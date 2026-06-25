//#include <QA.C>

#include <inttcalib/InttCalib.h>

#include <ffamodules/HeadReco.h>
#include <ffamodules/FlagHandler.h>
#include <ffamodules/SyncReco.h>
#include <ffamodules/CDBInterface.h>

#include <ffarawmodules/InttCheck.h>
#include <ffarawmodules/StreamingCheck.h>
#include <ffarawmodules/TpcCheck.h>

#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllUtils.h>

#include <fun4allraw/Fun4AllStreamingInputManager.h>
#include <fun4allraw/InputManagerType.h>
#include <fun4allraw/SingleGl1PoolInput.h>
#include <fun4allraw/SingleInttPoolInput.h>
#include <fun4allraw/SingleMicromegasPoolInput.h>
#include <fun4allraw/SingleMvtxPoolInput.h>
#include <fun4allraw/SingleTpcPoolInput.h>
#include <fun4allraw/SingleTpcTimeFrameInput.h>

#include <phool/recoConsts.h>
#include <GlobalVariables.C>

#include <phgarfield/PHGarfield.h>

#include <TPolyLine3D.h>
#include <TPolyLine.h>
#include <TGeoTube.h>
#include <TCanvas.h>
#include <TBox.h>
#include <TFile.h>
#include <cmath>
#include <vector>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libfun4allutils.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libfun4allraw.so)
R__LOAD_LIBRARY(libffarawmodules.so)
R__LOAD_LIBRARY(libcdbobjects.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libPHGarfield.so)

#define Nebdc 24
#define Nserver 2

// Global namespace to assist drawing...
TPolyLine3D *npoly3[48];
TPolyLine3D *spoly3[48];
TPolyLine   *npoly2[48];
TPolyLine   *spoly2[48];
TGeoTube    *tubby;
TCanvas     *canny;
TCanvas     *canny2;

TBox        *boxer1;
TBox        *boxer2;


static double WrapPhi(double phi)
{
  const double two_pi = 2.0 * std::acos(-1.0);
  while (phi < 0.0)
  {
    phi += two_pi;
  }
  while (phi >= two_pi)
  {
    phi -= two_pi;
  }
  return phi;
}

static TH3F* makeGuardedStandardTH3F(const std::string& name, const std::string& title,
                                     int nphi, int nr, int nz,
                                     double phi_min, double phi_max,
                                     const std::vector<double>& r_edges,
                                     double z_min, double z_max)
{
  const double dphi = (phi_max - phi_min) / nphi;
  const double dz = (z_max - z_min) / nz;
  std::vector<double> phi_bins(nphi + 3);
  std::vector<double> z_bins(nz + 3);
  for (int i = 0; i < nphi + 3; ++i)
  {
    phi_bins[i] = phi_min - dphi + i * dphi;
  }
  for (int i = 0; i < nz + 3; ++i)
  {
    z_bins[i] = z_min - dz + i * dz;
  }
  return new TH3F(name.c_str(), title.c_str(), nphi + 2, phi_bins.data(), nr + 2,
                  r_edges.data(), nz + 2, z_bins.data());
}

static bool InterpolatePathAtZ(TPolyLine3D* path, double z_target, double& r_out, double& phi_out)
{
  const int n = path->GetN();
  if (n < 2)
  {
    return false;
  }
  float* p = path->GetP();
  for (int i = 0; i < n - 1; i++)
  {
    const double x1 = p[3 * i + 0];
    const double y1 = p[3 * i + 1];
    const double z1 = p[3 * i + 2];
    const double x2 = p[3 * (i + 1) + 0];
    const double y2 = p[3 * (i + 1) + 1];
    const double z2 = p[3 * (i + 1) + 2];

    const double zmin = std::min(z1, z2);
    const double zmax = std::max(z1, z2);
    if (z_target < zmin || z_target > zmax)
    {
      continue;
    }

    const double frac = (z2 == z1) ? 0.0 : (z_target - z1) / (z2 - z1);
    const double x = x1 + frac * (x2 - x1);
    const double y = y1 + frac * (y2 - y1);
    r_out = std::hypot(x, y);
    phi_out = WrapPhi(std::atan2(y, x));
    return true;
  }
  return false;
}

static void fillGuardBins(TH3F* hist)
{
  //n bins including guards (but not under/overflow)
  const int nphi = hist->GetNbinsX() ;
  const int nr = hist->GetNbinsY() ;
  const int nz = hist->GetNbinsZ() ;

  // Step 1: Fill z guards (skip phi and r guard regions)
  //0=underflow, 1=guard, 2=first real
  // n+1=overflow, n=guard, n-1=last real
  for (int ip = 2; ip < nphi; ++ip)
  {
    for (int ir = 2; ir < nr; ++ir)
    {
      // zmin guard (iz=0) gets values from z=1
      hist->SetBinContent(ip, ir, 1, hist->GetBinContent(ip, ir, 2));
      // zmax guard (iz=nz+1) gets values from z=nz
      hist->SetBinContent(ip, ir, nz, hist->GetBinContent(ip, ir, nz-1));
    }
  }

  // Step 2: Fill r guards (skip phi guard, but include z guards)
  for (int ip = 2; ip < nphi; ++ip)
  {
    for (int iz = 1; iz <= nz; ++iz)
    {
      // rmin guard (ir=0) gets values from r=1
      hist->SetBinContent(ip, 1, iz, hist->GetBinContent(ip, 2, iz));
      // rmax guard (ir=nr+1) gets values from r=nr
      hist->SetBinContent(ip, nr, iz, hist->GetBinContent(ip, nr-1, iz));
    }
  }

  // Step 3: Fill phi guards (include all other guards )
  for (int ir = 1; ir <= nr + 1; ++ir)
  {
    for (int iz = 1; iz <= nz + 1; ++iz)
    {
      // phimin guard (ip=0) gets values from phimax (ip=nphi)
      hist->SetBinContent(1, ir, iz, hist->GetBinContent(nphi-1, ir, iz));
      // phimax guard (ip=nphi+1) gets values from phimin (ip=1)
      hist->SetBinContent(nphi, ir, iz, hist->GetBinContent(2, ir, iz));
    }
  }
}

void Fun4All_Garfield_Static()
{
  recoConsts* rc = recoConsts::instance();

  rc->set_StringFlag("CDB_GLOBALTAG","FieldMapTest");
  rc->set_uint64Flag("TIMESTAMP",1);

  auto cdb = CDBInterface::instance();
  std::string url = cdb->getUrl("FIELDMAP_TRACKING");
  std::cout << "Field map URL:\n" << url << std::endl;

  Fun4AllServer *se = Fun4AllServer::instance();

  //Enable::QA  = false;
  Enable::CDB = true;
  
  // Register a whole slew of input managers...
  // NOTE:  This depends upon the requested files being in frog.  Ribbit.
  char nextinput[500];
  char nextfile[500];
  Fun4AllInputManager* in[Nebdc]; 
  for (unsigned int ebdc=0; ebdc<24; ebdc++)
    {
      for (unsigned int server=0; server<2; server++)
	{
	  sprintf(nextinput,"ebdc%02d_%01d",ebdc,server);
	  //sprintf(nextfile,"DST_STREAMING_EVENT_ebdc%02d_%01d_run3line_laser_ana540_nocdbtag_v001-00064890-00000.root",ebdc,server);  // Line Laser
	  sprintf(nextfile,"DST_STREAMING_EVENT_ebdc%02d_%01d_run3auau_ana514_nocdbtag_v001-00075570-00000.root",ebdc,server);          // AuAu Zero Field
	  std::cout << nextfile << " " << nextinput << endl;
	  in[ebdc] = new Fun4AllDstInputManager(nextinput);
	  in[ebdc]->fileopen(nextfile);
	  se->registerInputManager(in[ebdc]);
	}
    }

  // Now register a flag handler because MAAABE it will make the CDB work correctly?
  //SubsysReco *fh = new FlagHandler();
  //se->registerSubsystem(fh);
  
  // Register Tom's Garfield analysis module.
  PHGarfield *phg = new PHGarfield();
  se->registerSubsystem(phg);

  se->run(4);

  const int nPhi = 48;
  const int nR = 48;
  const int nZ = 40;
  const double zRange = 102.0;
  const double pi = std::acos(-1.0);
  const double two_pi = 2.0 * pi;
  const double dphi = two_pi / static_cast<double>(nPhi);
  const double dz = zRange / static_cast<double>(nZ);

  std::vector<double> r_edges(nR + 3);
  const double r0 = phg->GetRadius(0);
  const double r1 = phg->GetRadius(1);
  const double r_last = phg->GetRadius(nR - 1);
  const double r_prev = phg->GetRadius(nR - 2);
  const double first_gap = (r1 - r0);
  const double last_gap = (r_last - r_prev);

  r_edges[1] = r0 - 0.5*(first_gap);
  for (int i = 1; i < nR; ++i)
  {
    r_edges[i + 1] = 0.5 * (phg->GetRadius(i - 1) + phg->GetRadius(i));
  }
  r_edges[nR + 1] = r_last + 0.5*(last_gap);
  r_edges[0] = r_edges[1] - first_gap;
  r_edges[nR + 2] = r_edges[nR + 1] + last_gap;

  const double phi_min = 0.0;
  const double phi_max = two_pi;
  double z_launch[]={zRange,-zRange};

  std::string side[2];
  side[0] = "posz";
  side[1] = "negz";
  std::string sepAxis[] = {"P", "R", "Z"};

  TH3F* hDistR[2];
  TH3F* hDistP[2];
  TH3F* hDistZ[2];

  // Create histograms for each side
  for (int s = 0; s < 2; ++s)
  {
    const std::string sideName = side[s];
    const double z_launch = (s == 0) ? zRange : -zRange;
    const double z_min=(s==0)? 0:-zRange;
    const double z_max=(s==0)? zRange:0;


    hDistP[s] = makeGuardedStandardTH3F(
        std::string("hIntDistortion" + sepAxis[0] + "_" + sideName),
        std::string("Integrated " + sepAxis[0] + " Deflection (rad) drifting from (phi,r,z) to z=endcap;phi (rad);r (cm);z(cm) (" + sideName + " side)"),
        nPhi, nR, nZ, phi_min, phi_max, r_edges, z_min, z_max);

    hDistR[s] = makeGuardedStandardTH3F(
        std::string("hIntDistortion" + sepAxis[1] + "_" + sideName),
        std::string("Integrated " + sepAxis[1] + " Deflection (cm) drifting from (phi,r,z) to z=endcap;phi (rad);r (cm);z(cm) (" + sideName + " side)"),
        nPhi, nR, nZ, phi_min, phi_max, r_edges, z_min, z_max);

    hDistZ[s] = makeGuardedStandardTH3F(
        std::string("hIntDistortion" + sepAxis[2] + "_" + sideName),
        std::string("Integrated " + sepAxis[2] + " Deflection (cm) drifting from (phi,r,z) to z=endcap;phi (rad);r (cm);z(cm) (" + sideName + " side)"),
        nPhi, nR, nZ, phi_min, phi_max, r_edges, z_min, z_max);

      for (int ir = 2; ir < nR+2; ++ir)
      {
        const double r_launch = hDistR[s]->GetYaxis()->GetBinCenter(ir);
        for (int ip = 2; ip < nPhi+2; ++ip)
        {
          const double phi_launch = hDistR[s]->GetXaxis()->GetBinCenter(ip);
          const double x_launch = r_launch * std::cos(phi_launch);
          const double y_launch = r_launch * std::sin(phi_launch);
          TPolyLine3D* path = phg->ReverseDrift(x_launch, y_launch, z_launch);
          if (!path)
          {
            // no path came back?
            continue;
          }

          for (int iz = 2; iz < nZ+2; ++iz)
          {
            const double z_center = hDistR[s]->GetZaxis()->GetBinCenter(iz);
            if ((z_launch > 0.0 && z_center <= 0.0) || (z_launch < 0.0 && z_center >= 0.0))
            {
              continue;
            }

            double r_actual = 0.0;
            double phi_actual = 0.0;
            if (!InterpolatePathAtZ(path, z_center, r_actual, phi_actual))
            {
              continue;
            }

            double dphi_val = phi_actual - phi_launch;
            if (dphi_val > pi)
            {
              dphi_val -= two_pi;
            }
            else if (dphi_val < -pi)
            {
              dphi_val += two_pi;
            }

            const double dr_val = r_actual - r_launch;
            //const double rdphi_val = r_launch * dphi_val;

            hDistR[s]->Fill(phi_launch, r_launch, z_center, dr_val);
            hDistP[s]->Fill(phi_launch, r_launch, z_center, dphi_val);
            hDistZ[s]->Fill(phi_launch, r_launch, z_center, 0.0);
          }

          delete path;
        }
      }
    }
  for (int s=0;s<2;s++){
    fillGuardBins(hDistR[s]);
    fillGuardBins(hDistP[s]);
    fillGuardBins(hDistZ[s]);
  }


  // Write all histograms to output file
  TFile output_file("GarfieldStaticCorrectionMaps.root", "RECREATE");
  for (int s = 0; s < 2; s++)
  {

    hDistR[s]->Write();
    hDistP[s]->Write();
    hDistZ[s]->Write();
    std::cout << "Wrote " << side[s] << " histograms with " << hDistR[s]->GetEntries() << " entries for dr." << std::endl;
  }
  output_file.Close();
  std::cout << "Wrote StaticCorrectionMaps.root" << std::endl;
}

