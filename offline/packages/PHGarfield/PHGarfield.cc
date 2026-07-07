#include "PHGarfield.h"
#include <phool/phool.h>

#include <cdbobjects/CDBTTree.h>

#include <phfield/PHField3DCartesian.h>

#include <ffamodules/CDBInterface.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <TPolyLine3D.h>
#include <TVector3.h>
#include <TRotation.h>
#include <TFile.h>
#include <TH2.h>

#include <CLHEP/Units/SystemOfUnits.h>

#include <Garfield/ComponentUser.hh>
#include <Garfield/MediumMagboltz.hh>

#include <cmath>
#include <filesystem>
#include <functional>
#include <iostream>  // for basic_ostream, operat...
#include <vector>

PHGarfield::PHGarfield(const std::string& name,
                           const std::string& electricFieldMap,
                           double spaceChargeScale)
  : SubsysReco(name),
    //m_defaultGasfile("/sphenix/user/hemmick/gasfiles_20260624/Ar75_CF20_iso5.gas")
    m_defaultGasfile("/sphenix/user/hemmick/gasfiles_20260624"),
    m_electricFieldMap(electricFieldMap),
    m_spaceChargeScale(spaceChargeScale)
{
}


int PHGarfield::InitRun(PHCompositeNode* /*topNode*/)
{
  if (Verbosity() > 1)
  {
    std::cout << "PHGarfield::InitRun(PHCompositeNode *topNode) Initializing" << std::endl;
  }
  CDBInterface* m_cdb = CDBInterface::instance();

  //  Here we use the CDBInterface to set up the magnetic field map:
  std::string url = m_cdb->getUrl("FIELDMAP_TRACKING");
  m_field = new PHField3DCartesian(url, 1.0);

  //  Here we use the CDBInterface to set up the channel making of the TPC:
  std::string text = m_cdb->getUrl("TPC_FEE_CHANNEL_MAP");
  m_cdbTPCMAPttree = new CDBTTree(text);
  m_cdbTPCMAPttree->LoadCalibrations();


  // Load the optional axisymmetric space-charge field map.
  // Failure is non-fatal: Garfield then uses only the nominal 400 V/cm field.
  if (!m_electricFieldMap.empty())
  {
    if (!LoadElectricFieldCorrections(m_electricFieldMap))
    {
      std::cerr << PHWHERE << " Failed to load electric-field correction map: "
                << m_electricFieldMap << std::endl;
    }
  }


  //  Make the Garfield Component and register the methods that will interface to our fields...
  m_component = new Garfield::ComponentUser();
  m_component->SetMagneticField([this](double x, double y, double z, double& bx, double& by, double& bz)
                                { GetMagneticFieldTesla(x, y, z, bx, by, bz); });
  m_component->SetElectricField([this](double x, double y, double z, double& ex, double& ey, double& ez)
                                { GetElectricFieldVcm(x, y, z, ex, ey, ez); });
  
  // Here we fetch the gas from the CDB
  std::string gasfile = m_cdb->getUrl("PHGARFIELD_GAS");
  if (gasfile.empty() || !std::filesystem::exists(gasfile))
    {
      std::cerr << PHWHERE << " Missing CDB gasfile: " << gasfile << std::endl;
      std::cerr << PHWHERE << " Using default gasfile: " << m_defaultGasfile << std::endl;
      gasfile = m_defaultGasfile;
    }
  InitializeGas(gasfile);

  //InitializeGas("/direct/phenix+u/workarea/hemmick/code.sphenix/tkh/gas/gasfiles/");

  //  Diagnostic during code development...
  FillRadii();
  if (Verbosity() > 1)
  {
    PrintMaps();
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

void PHGarfield::FillRadii()
{
  //  Unload the pad map to get the radii in a handy location:
  for (unsigned int side = 0; side < 2; side++)
  {
    for (unsigned int sector = 0; sector < 12; sector++)
    {
      for (unsigned int fee = 0; fee < 26; fee++)
      {
        for (unsigned int channel = 0; channel < 256; channel++)
        {
          unsigned int key = (256 * (fee)) + channel;
          int layer = m_cdbTPCMAPttree->GetIntValue(key, "layer");
          double r = m_cdbTPCMAPttree->GetDoubleValue(key, "R") / CLHEP::cm;
          if (layer > 6)
          {
            radii[layer - 7] = r;
          }
        }
      }
    }
  }
}

void PHGarfield::PrintGarfield(double x, double y, double z) const
{
  double ex;
  double ey;
  double ez;
  double bx;
  double by;
  double bz;
  double vx;
  double vy;
  double vz;
  GetElectricFieldVcm(x, y, z, ex, ey, ez);
  GetMagneticFieldTesla(x, y, z, bx, by, bz);
  m_gas->ElectronVelocity(ex, ey, ez, bx, by, bz, vx, vy, vz);
  std::cout << " x:" << x
            << " y:" << y
            << " z:" << z
            << " ex:" << ex
            << " ey:" << ey
            << " ez:" << ez
            << " bx:" << bx
            << " by:" << by
            << " bz:" << bz
            << " vx:" << vx
            << " vy:" << vy
            << " vz:" << vz
            << std::endl;
}

void PHGarfield::PrintMaps() const
{
  //  Print out a few test points of the Garfield information
  PrintGarfield(0.0, 0.0, 0.1);
  PrintGarfield(0.0, 0.0, 100.0);
  PrintGarfield(0.0, 40.0, 100.1);
  PrintGarfield(0.0, 78.0, 010.1);

  //  Print out the pad coordinate map:
  int MAX = 10;
  int prints = 0;
  for (unsigned int side = 0; side < 2; side++)
  {
    for (unsigned int sector = 0; sector < 12; sector++)
    {
      for (unsigned int fee = 0; fee < 26; fee++)
      {
        for (unsigned int channel = 0; channel < 256; channel++)
        {
          unsigned int key = (256 * (fee)) + channel;
          int layer = m_cdbTPCMAPttree->GetIntValue(key, "layer");
          double phi = ((side == 1 ? 1 : -1) * (m_cdbTPCMAPttree->GetDoubleValue(key, "phi") - std::numbers::pi / 2.)) + ((sector % 12) * std::numbers::pi / 6);
          double r = m_cdbTPCMAPttree->GetDoubleValue(key, "R") / CLHEP::cm;

          phi = bounder(phi, PHI_MIN);

          if (layer > 6)
          {
            if (prints < MAX)
            {
              prints++;
              std::cout << " side: " << side;
              std::cout << " sector: " << sector;
              std::cout << " fee: " << fee;
              std::cout << " channel: " << channel;
              std::cout << " layer: " << layer;
              std::cout << " phi: " << phi;
              std::cout << " r: " << r;
              std::cout << std::endl;
            }
          }
        }
      }
    }
  }
}

void PHGarfield::MoveMagnet(double x, double y, double z){
magpos.SetXYZ(x,y,z);
return;
}
void PHGarfield::RotateMagnet(double theta_x, double theta_y, double theta_z){
magrot.RotateX(theta_x);
magrot.RotateY(theta_y);
magrot.RotateZ(theta_z);
return;
}
void PHGarfield::MoveTpc(double x, double y, double z){
printf("Translating TPC by (%f,%f,%f) cm\n",x,y,z);

tpcpos.SetXYZ(x,y,z);
return;
}
void PHGarfield::RotateTpc(double theta_x, double theta_y, double theta_z){
  printf("Rotating TPC by (%f,%f,%f) radians\n",theta_x,theta_y,theta_z);
tpcrot.RotateX(theta_x);
tpcrot.RotateY(theta_y);
tpcrot.RotateZ(theta_z);
return;
}

void PHGarfield::ConvertToLocal(double &x, double &y, double &z, TRotation rot, TVector3 trans) const{
  //printf("Convert to Local (%f,%f,%f)\n",x,y,z);

  //this assumes everything is in the same units!
  //convert coords from global coords in global axes
  //  to coords wrt tpc center (with global axes)
  TVector3 global;
  global.SetXYZ(x,y,z);
  TVector3 localRaw=global-trans;
  //rotate into the local axes:
  TRotation localRotInverse=rot.Inverse();
  TVector3 local=localRotInverse*localRaw;    
  x=local.X();
  y=local.Y();
  z=local.Z();
  return;
}
void PHGarfield::ConvertToGlobal(double &x, double &y, double &z, TRotation rot, TVector3 trans) const{
  printf("Convert to Global (%f,%f,%f)\n",x,y,z);
  //this assumes everything is in the same units!
  //inverse of the ConvertToLocal:
  TVector3 local;
  local.SetXYZ(x,y,z);
  //rotate back to global axes:
  TVector3 globalRaw=rot*local;
  TVector3 global=globalRaw+trans;
  x=global.X();
  y=global.Y();
  z=global.Z();
  return;
}

void PHGarfield::GetMagneticFieldTesla(double x_cm, double y_cm, double z_cm, double& bx_t, double& by_t, double& bz_t) const
{
cout <<PHWHERE << "enter";
  // NOTE:  Garfield uses  cm, V/cm, and Tesla.
  //        CLHEP    uses  mm, V/mm, and kiloTesla
  //        PHField3DCartesian follows the CLHEP conventions for magnetic fields.
/*
  //find the coordinates in global axes centered on magnet center:
  TVector3 raw;
  raw.SetXYZ(x_cm,y_cm,z_cm);
  TVector3 magraw=raw-magpos;
  //rotate into the magnet coordinates:
  TRotation magrotInverse=magrot.Inverse();
  TVector3 magcoord=magrotInverse*magraw;
  */
  double x=x_cm,y=y_cm, z=z_cm;
  ConvertToLocal(x,y,z,magrot,magpos);

  double point[4] =
      {
          x * CLHEP::cm,
          y * CLHEP::cm,
          z * CLHEP::cm,
          //(z_cm-20.0) * CLHEP::cm,
          0.0};

  double bfield[3] = {0.0, 0.0, 0.0};

  //  Get the magnetic field via the PHField3DCartesian object constructed usinf the CDB url reference.
  m_field->GetFieldValue(point, bfield);

  //bx_t = bfield[0] / CLHEP::tesla;
  //by_t = bfield[1] / CLHEP::tesla;
  //bz_t = bfield[2] / CLHEP::tesla;

TVector3 bfieldMag;
bfieldMag.SetXYZ(bfield[0],bfield[1],bfield[2]);
TVector3 bfieldGlobal=magrot*bfieldMag;
bx_t = bfieldGlobal.X() / CLHEP::tesla;
by_t = bfieldGlobal.Y() / CLHEP::tesla;
bz_t = bfieldGlobal.Z() / CLHEP::tesla;

cout <<PHWHERE << "exit";
return;
}


void PHGarfield::GetElectricFieldVcm(double x_cm, double y_cm, double z_cm, double& ex_vcm, double& ey_vcm, double& ez_vcm) const
{
cout <<PHWHERE << "enter";
  double x=x_cm,y=y_cm, z=z_cm;
  ConvertToLocal(x,y,z,tpcrot,tpcpos);


  double ex_loc,ey_loc,ez_loc;
  GetTpcFrameElectricFieldVcm(x,y,z,ex_loc,ey_loc,ez_loc);
  TVector3 fieldLocal;
  fieldLocal.SetXYZ(ex_loc,ey_loc,ez_loc); 
  //field is in TPC axes.  rotate into global axes:
  TVector3 fieldGlobal=tpcrot*fieldLocal;

  ex_vcm = fieldGlobal.X();
  ey_vcm = fieldGlobal.Y();
  ez_vcm = fieldGlobal.Z();
cout <<PHWHERE << "exit";

  return;
}

void PHGarfield::GetTpcFrameElectricFieldVcm(double x_cm, double y_cm, double z_cm, double& ex_vcm, double& ey_vcm, double& ez_vcm) const
{
cout <<PHWHERE << "enter";
  // NOTE:  Garfield uses  cm, V/cm, and Tesla.
  // The notebook maps use cm on their axes and V/m in their bins.
  // The map is produced for one TPC half using s = |z|, measured from
  // the central membrane toward the readout plane.

  const double r_cm = std::hypot(x_cm, y_cm);
  const double abs_z_cm = std::abs(z_cm);


  ex_vcm = 0.0;
  ey_vcm = 0.0;
  ez_vcm = z_cm > 0 ? -400.0 : 400.0;
  return;

  //Yuri's correction:
    if (!m_erCorrection || !m_ezCorrection || m_spaceChargeScale == 0.0)
  {
    return;
  }

  const double delta_er_vcm = m_spaceChargeScale *
      InterpolateCorrectionVcm(m_erCorrection, r_cm, abs_z_cm);
  const double delta_ez_local_vcm = m_spaceChargeScale *
      InterpolateCorrectionVcm(m_ezCorrection, r_cm, abs_z_cm);

  // Convert the cylindrical radial correction to Cartesian components.
  if (r_cm > 0.0)
  {
    ex_vcm += delta_er_vcm * x_cm / r_cm;
    ey_vcm += delta_er_vcm * y_cm / r_cm;
  }

  // hEzDefault is expressed along the local coordinate s = |z|.
  // Convert it to the global Cartesian z direction.
  ez_vcm += z_cm >= 0.0 ? delta_ez_local_vcm : -delta_ez_local_vcm;
cout <<PHWHERE << "exit";
  return;
}

bool PHGarfield::LoadElectricFieldCorrections(const std::string& filename)
{
  std::unique_ptr<TFile> input(TFile::Open(filename.c_str(), "READ"));
  if (!input || input->IsZombie())
  {
    std::cerr << PHWHERE << " Could not open electric-field map: "
              << filename << std::endl;
    return false;
  }

  auto* er = dynamic_cast<TH2*>(input->Get("QA/hErDefault"));
  auto* ez = dynamic_cast<TH2*>(input->Get("QA/hEzDefault"));

  // Also allow maps written at the ROOT-file top level.
  if (!er) er = dynamic_cast<TH2*>(input->Get("hErDefault"));
  if (!ez) ez = dynamic_cast<TH2*>(input->Get("hEzDefault"));

  if (!er || !ez)
  {
    std::cerr << PHWHERE
              << " Missing QA/hErDefault or QA/hEzDefault in "
              << filename << std::endl;
    return false;
  }

  delete m_erCorrection;
  delete m_ezCorrection;
  m_erCorrection = dynamic_cast<TH2*>(er->Clone("PHGarfield_ErCorrection"));
  m_ezCorrection = dynamic_cast<TH2*>(ez->Clone("PHGarfield_EzCorrection"));

  if (!m_erCorrection || !m_ezCorrection)
  {
    delete m_erCorrection;
    delete m_ezCorrection;
    m_erCorrection = nullptr;
    m_ezCorrection = nullptr;
    return false;
  }

  m_erCorrection->SetDirectory(nullptr);
  m_ezCorrection->SetDirectory(nullptr);

  std::cout << "Loaded axisymmetric electric-field corrections from "
            << filename << std::endl;
  std::cout << "  scale k_eff = " << m_spaceChargeScale << std::endl;
  std::cout << "  r range [cm] = ["
            << m_erCorrection->GetXaxis()->GetXmin() << ", "
            << m_erCorrection->GetXaxis()->GetXmax() << "]" << std::endl;
  std::cout << "  |z| range [cm] = ["
            << m_erCorrection->GetYaxis()->GetXmin() << ", "
            << m_erCorrection->GetYaxis()->GetXmax() << "]" << std::endl;

  return true;
}

double PHGarfield::InterpolateCorrectionVcm(const TH2* hist,
                                             double r_cm,
                                             double abs_z_cm) const
{
  if (!hist) return 0.0;

  const auto* xaxis = hist->GetXaxis();
  const auto* yaxis = hist->GetYaxis();
  if (r_cm < xaxis->GetXmin() || r_cm > xaxis->GetXmax() ||
      abs_z_cm < yaxis->GetXmin() || abs_z_cm > yaxis->GetXmax())
  {
    return 0.0;
  }

  // Notebook histograms store V/m; Garfield expects V/cm.
  return hist->Interpolate(r_cm, abs_z_cm) / 100.0;
}


void PHGarfield::InitializeGas(const std::string &dir)
{
  //  Create and fill the gas object so that we can trace particles through the gas...
  m_gas = new Garfield::MediumMagboltz();

  auto filename = [&](const int i)
  { return dir + "/PART_" + std::to_string(i) + ".gas"; };

  const std::string first = filename(0);
  if (!std::filesystem::exists(first))
  {
    std::cerr << "Missing first gas file: " << first << std::endl;
    return;
  }

  if (!m_gas->LoadGasFile(first))
  {
    std::cerr << "Failed to load " << first << std::endl;
    return;
  }

  for (int i = 1;; ++i)
  {
    const std::string file = filename(i);

    if (!std::filesystem::exists(file))
    {
      std::cout << "Stopping at first missing file: " << file << std::endl;
      break;
    }

    std::cout << "Merging " << file << std::endl;

    if (!m_gas->MergeGasFile(file, true))
    {
      std::cerr << "Failed to merge " << file << std::endl;
      return;
    }
  }
}

int PHGarfield::process_event(PHCompositeNode*)
{
  // Initial implementation doesn't do anything event-by-event.
  // Nonetheless, a future user might want do do something here...

  return Fun4AllReturnCodes::EVENT_OK;
}

double PHGarfield::bounder(double phi, double phi_min)
{
  double phi_max = phi_min + 2.0 * std::numbers::pi;
  while (phi < phi_min)
  {
    phi = phi + 2.0 * std::numbers::pi;
  }
  while (phi >= phi_max)
  {
    phi = phi - 2.0 * std::numbers::pi;
  }

  return phi;
}

TPolyLine3D* PHGarfield::ReverseDriftTpcCoords(double x_cm, double y_cm, double z_cm, double step_ns)
{
  printf("ReverseDrifting (TPC Coords) (%f,%f,%f, step=%f)\n",x_cm,y_cm,z_cm,step_ns);

  //x,y,z are denominated in tpc coordinate, so transform them to global
  double x=x_cm,y=y_cm, z=z_cm;
  ConvertToGlobal(x,y,z,tpcrot,tpcpos);

  TPolyLine3D* poly = ReverseDrift(x,y,z,step_ns);
  //polyline is in global coordinates, so transform it back, point by point.
  printf("ReverseDrifting (TPC Coords) Polyline has n=%d\n",poly->GetN());

  for (int i = 0; i < poly->GetN(); i++)
  {
    if (!(i%100)) printf("i=%d",i);
    double polyX=poly->GetP()[i*3];
    double polyY=poly->GetP()[i*3+1];
    double polyZ=poly->GetP()[i*3+2];
    ConvertToLocal(polyX,polyY,polyZ,tpcrot,tpcpos);
    
    poly->SetPoint(i, polyX,polyY,polyZ);
  }

  return poly;
}

TPolyLine3D* PHGarfield::ReverseDrift(double x, double y, double z, double step_ns)
{
  printf("ReverseDrifting (Global Coords) (%f,%f,%f, step=%f)\n",x,y,z,step_ns);

  std::vector<double> xlist;
  std::vector<double> ylist;
  std::vector<double> zlist;

  xlist.push_back(x);
  ylist.push_back(y);
  zlist.push_back(z);

  double ex;
  double ey;
  double ez;
  double bx;
  double by;
  double bz;
  double vx;
  double vy;
  double vz;

  double zPrevious = z;
  while (!StopHere(x, y, z, zPrevious))
  {
    zPrevious = z;
    GetMagneticFieldTesla(x, y, z, bx, by, bz);
    GetElectricFieldVcm(x, y, z, ex, ey, ez);
    m_gas->ElectronVelocity(ex, ey, ez, bx, by, bz, vx, vy, vz);

    x = x - vx * step_ns;
    y = y - vy * step_ns;
    z = z - vz * step_ns;

    xlist.push_back(x);
    ylist.push_back(y);
    zlist.push_back(z);
  }

  TPolyLine3D* poly = new TPolyLine3D(xlist.size() - 1);
  for (unsigned int i = 0; i < xlist.size() - 1; i++)
  {
    poly->SetPoint(i, xlist[i], ylist[i], zlist[i]);
  }

  return poly;
}

bool PHGarfield::StopHere(const double x, const double y, const double z,
                          const double zPrevious)
{
  const double r = std::hypot(x, y);

  if (r < 18.0)
  {
    return true;
  }
  if (r > 82.0)
  {
    return true;
  }
  if (z > 120.0)
  {
    return true;
  }
  if (z < -120.0)
  {
    return true;
  }

  // z crossed the central membrane.
  if (z * zPrevious < 0.0)
  {
    return true;
  }

  return false;
}
