#include "TunableField.h"
#include "RosseggerReader.h"
#include "PadrowReader.h"

#include <format>
#include <TH3.h>
#include <TTree.h>
#include <TFile.h>

/*
A) Handling the primary space charge distribution:
1) get the primary space charge distribution from the fillSpaceChargeMaps output (I'll find a suitable pet file)
2) smooth by summing a bunch of them together and normalizing by the number of files we add.
3) smooth that averaged primary space charge by making it rotationally symmetric around z.
4) produce the electric field due to that smoothed charge by doing appropriate matrix multiplication with the Green's function (which is the rossegger data)

B) Handling the ion back flow (IBF) space charge distribution:
1) take our smoothed average primary space charge and 
smooth it in z.
2) make a decision about what to do with the charge that is not over an active area of the readout:
    i) divide those inactive regions into the parts that are closest to each padrow's active region, and assign the charge to the edge of that active region
    ii) simply remove the charge that is not over an active area of the readout
3) the resulting space charge may not be smooth any more.  For each padrow in one sector of the TPC, we compute the electric field due to the space charge we have assigned to that padrow.

C) Producing the total electric field:
1) We take in, for now, a small number of parameters:
    i) the primary ionization scale factor (which is the factor with which we scale the primary SC electric field
    ii) the ion back flow gain (this is the number of backflow ions per primary ion, and used to scale the IBF electric field for every padrow)
    iii) the external electric field (a z-directed electric field with no spatial variation except that it points toward z=0 (so that electrons travel away from z=0)
    iv) (and we have a function to knock-out any padrows we wish to set to zero IBF)

2) For each sector of the TPC, we create the sector IBF electric field by summing the IBF electric field for each padrow in a sector (unless knocked out), scaled by the ion back flow gain.
3) for each side of the TPC, we create the Side IBF electric field by summing the sector IBF electric fields for that side, rotating each by the appropriate amount in phi.
4) we create the total electric field by summing the primary SC electric field (scaled by the primary ionization scale factor) and the Side IBF electric field for each side.


*/


/*TODO:
- get the field-computer to do either single-phi (and then we enforce rotational symmetry manually) or the full phi range all-to-all (for the padrows).
- get RosseggerReader to do both positive and negative z by shifting the source and target points by the z extent if they're out of bounds.
- fix the handling of padrow boundaries.  We need to get the actual padrow boundaries from the PHG4TpcCylindericalGeomContainer
- implement a maximum distance in the field calculation, to skip source bins that are too far from the target point to contribute meaningfully to the field.  (make a little diagnostic test that shows the total field components as a function of the cutoff distance in the limiting case of perfectly uniform charge density?

-(auto) Implement automatic detection of the number of histograms present in the space charge file to ensure correct normalization in loadPrimarySpaceCharge and loadIBFSpaceCharge.
-(auto) Optimize computeField to pre-sum the Green's function azimuthal components when the source charge is rotationally symmetric, significantly reducing integration time.
-(auto) Generalize the frame bin skipping logic in makeIBFRotationallyPeriodic to be geometry-driven rather than relying on "magic" bin indices.
-(auto) Implement the PadrowReader logic to account for dead areas (e.g. support structures and module gaps) when calculating the active fraction of a bin.
*/

TunableField::TunableField()
{
  m_padrowReader = new PadrowReader();
  m_rosseggerReader = new RosseggerReader();

  m_primaryCharge = makeUnguardedStandardTH3F("h_primary_charge", "Primary Charge");
  m_ibfCharge = makeUnguardedStandardTH3F("h_ibf_charge", "IBF Charge");

  m_primaryFieldR = makeGuardedStandardTH3F("h_field_r", "E_{r}");
  m_primaryFieldP = makeGuardedStandardTH3F("h_field_p", "E_{\phi}");
  m_primaryFieldZ = makeGuardedStandardTH3F("h_field_z", "E_{z}");

  m_totalIbfFieldR = makeGuardedStandardTH3F("h_ibf_field_r", "E_{r} IBF");
  m_totalIbfFieldP = makeGuardedStandardTH3F("h_ibf_field_p", "E_{\phi} IBF");
  m_totalIbfFieldZ = makeGuardedStandardTH3F("h_ibf_field_z", "E_{z} IBF");
}

TunableField::~TunableField()
{
  delete m_padrowReader;
  delete m_rosseggerReader;
  delete m_primaryCharge;
  delete m_ibfCharge;
  delete m_primaryFieldR;
  delete m_primaryFieldP;
  delete m_primaryFieldZ;
  delete m_totalIbfFieldR;
  delete m_totalIbfFieldP;
  delete m_totalIbfFieldZ;
  clearPadrowFields();
  if (m_diagFile) m_diagFile->Close();
}

void TunableField::clearPadrowFields()
{
  for (TH3F* h : m_padrowCharges) delete h;
  for (TH3F* h : m_padrowFieldR) delete h;
  for (TH3F* h : m_padrowFieldP) delete h;
  for (TH3F* h : m_padrowFieldZ) delete h;
  m_padrowCharges.clear();
  m_padrowFieldR.clear();
  m_padrowFieldP.clear();
  m_padrowFieldZ.clear();
}

void TunableField::setDiagnosticFile(const std::string& filename)
{
  m_diagFile = TFile::Open(filename.c_str(), "RECREATE");
}

void TunableField::saveDiagnosticSpaceCharge(const std::string& name, TH3F* h)
{
  if (!m_diagFile) return;
  m_diagFile->cd();
  h->Write(name.c_str());
}

void TunableField::loadPrimarySpaceCharge()
{
  if (m_scFilename.empty()) return;

  TFile* f = TFile::Open(m_scFilename.c_str());
  if (!f) return;

  m_primaryCharge->Reset();
  float norm = 1.0 / 30.0;
  for (int i = 0; i < 30; ++i)
  {
    std::string hname = std::format("_h_SC_prim_{}", i);
    TH3* h = (TH3*) f->Get(hname.c_str());
    if (h)
    {
      for (int ir = 0; ir < m_nr; ++ir)
      {
        for (int ip = 0; ip < m_nphi; ++ip)
        {
          for (int iz = 0; iz < m_nz; ++iz)
          {
            double phi = m_primaryCharge->GetXaxis()->GetBinCenter(ip + 1);
            double r = m_primaryCharge->GetYaxis()->GetBinCenter(ir + 1);
            double z = m_primaryCharge->GetZaxis()->GetBinCenter(iz + 1);
            float q = h->GetBinContent(h->FindBin(phi, r, z));
            
            float currentQ = m_primaryCharge->GetBinContent(ip + 1, ir + 1, iz + 1);
            m_primaryCharge->SetBinContent(ip + 1, ir + 1, iz + 1, currentQ + q * norm);
          }
        }
      }
    }
  }
  f->Close();
  saveDiagnosticSpaceCharge("h_primary_charge_native", m_primaryCharge);
}

void TunableField::loadIBFSpaceCharge()
{
  if (m_scFilename.empty()) return;

  TFile* f = TFile::Open(m_scFilename.c_str());
  if (!f) return;

  m_ibfCharge->Reset();
  float norm = 1.0 / 30.0;
  for (int i = 0; i < 30; ++i)
  {
    std::string hname = std::format("_h_SC_ibf_{}", i);
    TH3* h = (TH3*) f->Get(hname.c_str());
    if (h)
    {
      for (int ir = 0; ir < m_nr; ++ir)
      {
        for (int ip = 0; ip < m_nphi; ++ip)
        {
          for (int iz = 0; iz < m_nz; ++iz)
          {
            double phi = m_ibfCharge->GetXaxis()->GetBinCenter(ip + 1);
            double r = m_ibfCharge->GetYaxis()->GetBinCenter(ir + 1);
            double z = m_ibfCharge->GetZaxis()->GetBinCenter(iz + 1);
            float q = h->GetBinContent(h->FindBin(phi, r, z));
            
            float currentQ = m_ibfCharge->GetBinContent(ip + 1, ir + 1, iz + 1);
            m_ibfCharge->SetBinContent(ip + 1, ir + 1, iz + 1, currentQ + q * norm);
          }
        }
      }
    }
  }
  f->Close();
  saveDiagnosticSpaceCharge("h_ibf_charge_native", m_ibfCharge);
}

void TunableField::loadRossegger()
{
  std::string filename = "/sphenix/user/rcorliss/rossegger/ross_phi1_sphenix_phislice_lookup_r26xp48xz40.root";
  if (m_rosseggerReader) m_rosseggerReader->Load(filename);
}

TVector3 TunableField::getInterpolatedGreen(float phi_src, float r_src, float z_src, float r_tgt, float z_tgt)
{
  if (m_rosseggerReader)
  {
    return m_rosseggerReader->GetInterpolatedField(phi_src, r_src, z_src, r_tgt, z_tgt);
  }
  return TVector3(0, 0, 0);
}

void TunableField::rotateAndAverage(TH3F* h)
{
  // make this charge distribution azimuthally symmetric by averaging over phi. 
  for (int ir = 1; ir <= m_nr; ++ir)
  {
    for (int iz = 1; iz <= m_nz; ++iz)
    {
      float phiSum = 0;
      for (int ip = 1; ip <= m_nphi; ++ip)
      {
        phiSum += h->GetBinContent(ip, ir, iz);
      }
      float phiAvg = phiSum / m_nphi;
      for (int ip = 1; ip <= m_nphi; ++ip)
      {
        h->SetBinContent(ip, ir, iz, phiAvg);
      }
    }
  }
}

void TunableField::makePrimaryInto2D()
{
  rotateAndAverage(m_primaryCharge);
  saveDiagnosticSpaceCharge("h_primary_charge_2D", m_primaryCharge);
}

void TunableField::makeIBFinto2D()
{
  // copy that native r and phi binning, and sum the ions across z. 
  // The ions in a given voxel are just the total charge in the column, 
  // divided by the volume of the column, time the volume of the particular voxel.
  for (int ip = 1; ip <= m_nphi; ++ip)
  {
    for (int ir = 1; ir <= m_nr; ++ir)
    {
      float zSum = 0;
      for (int iz = 1; iz <= m_nz; ++iz)
      {
        zSum += m_ibfCharge->GetBinContent(ip, ir, iz);
      }
      float zAvg = zSum / m_nz;
      for (int iz = 1; iz <= m_nz; ++iz)
      {
        m_ibfCharge->SetBinContent(ip, ir, iz, zAvg);
      }
    }
  }
  saveDiagnosticSpaceCharge("h_ibf_charge_2D", m_ibfCharge);
}

void TunableField::makeIBFRotationallyPeriodic()
{
  // make an azimuthally periodic copy by going over the input bins once 
  // and putting the charge in the phi=phi, phi+pi/6, phi+2pi/6, etc bins, 
  // (with wrapping) then divide contents by 12.

//magic rules for Evgeny's /sphenix/user/shulga/Work/IBF/DistortionMap/Files/Summary_hist_mdc2_UseFieldMaps_AA_event_0_bX10556072.root file:
int bin_period=17; //bins corresponding to 30 degrees in phi; phibins 0 and 17 are both frame bins.
int bin_wraparound_bonus=1; //the first and very last bin are the two halves of a single frame bin, so if we wrap around we have to add that bonus bin.


  
  TH3F* temp = (TH3F*)m_ibfCharge->Clone("h_ibf_periodic_temp");
  temp->Reset();//zero out the spacecharge to start.
  
  double dphi_sector = 2.0 * M_PI / 12.0;
  
  for (int ip = 1; ip <= m_nphi; ++ip)
  {
    double phi = m_ibfCharge->GetXaxis()->GetBinCenter(ip);
    double phispan = m_ibfCharge->GetXaxis()->GetBinWidth(ip);
    if (phispan<0.03) continue; //ignore charge in the frame bins per Tom.
    //this is hacky and matches Evgeny's phibinning.
    for (int ir = 1; ir <= m_nr; ++ir)
    {
      double rspan = m_ibfCharge->GetYaxis()->GetBinWidth(ir);
      if (rspan>7 && rspan<10) continue; ignore radial frame bins per Tom.
      //this is still hacky.  Need to get this from the geom properly, so we can get dead areas.
      for (int iz = 1; iz <= m_nz; ++iz)
      {
        float q = m_ibfCharge->GetBinContent(ip, ir, iz);
        if (q == 0) continue;
        
        for (int k = 0; k < 12; ++k)
        {
          double p_target = phi + k * dphi_sector;
          while (p_target >= 2.0 * M_PI) p_target -= 2.0 * M_PI;
          int ip_target = temp->GetXaxis()->FindBin(p_target);
          float q_prev = temp->GetBinContent(ip_target, ir, iz);
          temp->SetBinContent(ip_target, ir, iz, q_prev + q);
        }
      }
    }
  }
  
  // now divide by 12 and copy back
  m_ibfCharge->Reset();
  for (int ip = 1; ip <= m_nphi; ++ip)
    for (int ir = 1; ir <= m_nr; ++ir)
      for (int iz = 1; iz <= m_nz; ++iz)
        m_ibfCharge->SetBinContent(ip, ir, iz, temp->GetBinContent(ip, ir, iz) / 12.0);

  delete temp;
  saveDiagnosticSpaceCharge("h_ibf_charge_periodic", m_ibfCharge);
}

void TunableField::loadPadrowBoundaries()
{
  // get the geometry from the sPHENIX codebase and use it to determine the radial boundaries of each padrow, 
  // and store those boundaries in m_padrowBoundaries.
  //but for now, noop;
  return;
}

void TunableField::adjustIBFToActiveArea()
{
  // determine which (phi,r) locations are over active readout and which are not, 
  // and set the charge to zero where there is no active readout.

  //but for now, just noop;
  return;
  saveDiagnosticSpaceCharge("h_ibf_charge_active_trimmed", m_ibfCharge);
}

void TunableField::computeField(TH3F* sourceCharge, TH3F* fieldR, TH3F* fieldP, TH3F* fieldZ)
{
  const double eps0 = 8.854e-14; // C / (V * cm)
  const double epsinv = 1.0 / eps0;
  double field_phi = sourceCharge->GetXaxis()->GetBinCenter(1); 

  for (int ir = 0; ir < m_nr; ++ir)
  {
    double field_r = sourceCharge->GetYaxis()->GetBinCenter(ir + 1);
    for (int iz = 0; iz < m_nz; ++iz)
    {
      double field_z = sourceCharge->GetZaxis()->GetBinCenter(iz + 1);
      TVector3 field(0, 0, 0);
      for (int sr = 0; sr < m_nr; ++sr)
      {
        double source_r = sourceCharge->GetYaxis()->GetBinCenter(sr + 1);
        for (int sp = 0; sp < m_nphi; ++sp)
        {
          double source_phi = sourceCharge->GetXaxis()->GetBinCenter(sp + 1);
          for (int sz = 0; sz < m_nz; ++sz)
          {
            float q = sourceCharge->GetBinContent(sp + 1, sr + 1, sz + 1);
            if (q == 0) continue;
            double source_z = sourceCharge->GetZaxis()->GetBinCenter(sz + 1);

            TVector3 unitField;
            if (m_rosseggerReader)
            {
              unitField = m_rosseggerReader->GetInterpolatedField((float)source_phi, (float)source_r, (float)source_z, (float)field_r, (float)field_z);
            }
            else
            {
              printf("Warning: no RosseggerReader, returning zero field\n");
              unitField = TVector3(0, 0, 0);
            }
            field += unitField * q;
          }
        }
      }
      field *= epsinv;

      for (int ip = 0; ip < m_nphi; ++ip)
      {
        fieldR->SetBinContent(ip + 2, ir + 2, iz + 2, field.X());
        fieldP->SetBinContent(ip + 2, ir + 2, iz + 2, field.Y());
        fieldZ->SetBinContent(ip + 2, ir + 2, iz + 2, field.Z());
      }
    }
  }
  FillGuardBins(fieldR);
  FillGuardBins(fieldP);
  FillGuardBins(fieldZ);
}

void TunableField::computePadrowIBFFields()
{
  // for each padrow in one sector of the TPC, we compute the electric field 
  // due to the space charge we have assigned to that padrow.
  for (int i = 0; i < 48; ++i) {
    computeField(m_padrowCharges[i], m_padrowFieldR[i], m_padrowFieldP[i], m_padrowFieldZ[i]);
  }
}

void TunableField::computePrimaryField()
{
  /* Note(s):
  - we're only generating the primary field for one slice of phi, since by construction the primaries are rotationally symmetric.
  - I thought I might have more than one bullet. ;)
  - the Er (etc) functions from Rossegger return the (del V) gradient, so we need -1* to get the actual field (hence the subtraction, which I had to stare at :| )
  */
  computeField(m_primaryCharge, m_primaryFieldR, m_primaryFieldP, m_primaryFieldZ);
}

void TunableField::setPadrowKnockout(int padrow, bool knockout)
{
  if (padrow < 0) return;
  if (m_padrowKnockout.size() <= (size_t)padrow)
  {
    m_padrowKnockout.resize(padrow + 1, false);
  }
  m_padrowKnockout[padrow] = knockout;
}

void TunableField::FillGuardBins(TH3F* h)
{
  int nP = h->GetNbinsX();
  int nR = h->GetNbinsY();
  int nZ = h->GetNbinsZ();

  // Guard bins at z bounds (copy adjacent)
  for (int ip = 1; ip <= nP; ++ip)
  {
    for (int ir = 1; ir <= nR; ++ir)
    {
      h->SetBinContent(ip, ir, 1, h->GetBinContent(ip, ir, 2));
      h->SetBinContent(ip, ir, nZ, h->GetBinContent(ip, ir, nZ - 1));
    }
  }
  // Guard bins at r bounds (copy adjacent)
  //(since r/phi data at z bounds are valid now, this also fills the corners of the guard bins at r/z boundaries)
  for (int ip = 1; ip <= nP; ++ip)
  {
    for (int iz = 1; iz <= nZ; ++iz)
    {
      h->SetBinContent(ip, 1, iz, h->GetBinContent(ip, 2, iz));
      h->SetBinContent(ip, nR, iz, h->GetBinContent(ip, nR - 1, iz));
    }
  }
  // Guard bins at phi bounds (periodic wrap)
  for (int ir = 1; ir <= nR; ++ir)
  {
    for (int iz = 1; iz <= nZ; ++iz)
    {
      h->SetBinContent(1, ir, iz, h->GetBinContent(nP - 1, ir, iz));
      h->SetBinContent(nP, ir, iz, h->GetBinContent(2, ir, iz));
    }
  }
}

void TunableField::calculateFieldContributions()
{
//generate all the space charge distributions we need:

//Rossegger lookup (needed for field calculations)
loadRossegger();
setDiagnosticFile("tunableField_intermediates.hist.root");

//primary space charge (3D, one phi slice is all we use)
loadPrimarySpaceCharge();//load it in in its native binning
makePrimaryInto2D();//copy that native r and z binning, and sum the ions across phi.  The ions in a given voxel are just the total charge in the ring, divided by the volume of the ring, time the volume of the particular voxel. (this automatically makes us rotationally symmetric)

//IBF space charge (2D, no z dependent, 1/12 of azimuth is all we use)
loadIBFSpaceCharge();//load it in in its native binning
makeIBFinto2D();// copy that native r and phi binning, and sum the ions across z.  The ions in a given voxel are just the total charge in the column, divided by the volume of the column, time the volume of the particular voxel. (this automatically makes us smooth in z)
makeIBFRotationallyPeriodic();//make an azimuthally periodic copy by going over the input bins once and putting the charge in the phi=phi, phi+pi/6, phi+2pi/6, etc bins, (with wrapping) then divide contents by 12.

//get the geometry from the sPHENIX codebase, and use it to determine which (phi,r) locations are over active readout and which are not, and set the charge to zero where there is no active readout.
loadPadrowBoundaries();
adjustIBFToActiveArea();

//compute the electric field in the active half of the TPC for each distribution:
computePrimaryField();
computePadrowIBFFields();
m_componentFieldsAreReady = true;
}

void TunableField::generateTotalField()
{
  if(!m_componentFieldsAreReady){
    calculateFieldContributions();
  }

  // Sum the padrow fields into the total IBF field, respecting knockouts.
  m_totalIbfFieldR->Reset();
  m_totalIbfFieldP->Reset();
  m_totalIbfFieldZ->Reset();

  for (int i = 0; i < (int)m_padrowFieldR.size(); ++i)
  {
    bool knockout = (i < (int)m_padrowKnockout.size()) ? m_padrowKnockout[i] : false;
    if (!knockout)
    {
      m_totalIbfFieldR->Add(m_padrowFieldR[i]);
      m_totalIbfFieldP->Add(m_padrowFieldP[i]);
      m_totalIbfFieldZ->Add(m_padrowFieldZ[i]);
    }
  }

  // Diagnostic output for the assembled fields after knockouts but before overall scaling.
  if (m_diagFile)
  {
    m_diagFile->cd();
    m_primaryFieldR->Write("h_primary_field_R");
    m_primaryFieldP->Write("h_primary_field_P");
    m_primaryFieldZ->Write("h_primary_field_Z");

    m_totalIbfFieldR->Write("h_total_ibf_field_R");
    m_totalIbfFieldP->Write("h_total_ibf_field_P");
    m_totalIbfFieldZ->Write("h_total_ibf_field_Z");
  }

  m_fieldIsReady=true;
}

TVector3 TunableField::getFieldAt(const TVector3& pos)
{
if(!m_fieldIsReady){
  generateTotalField();
}

  TVector3 totalField(0, 0, (pos.Z() > 0 ? m_externalE : -m_externalE));

  double phi = pos.Phi();
  if (phi < 0) phi += 2.0 * M_PI;
  double r = pos.Perp();
  double z = std::abs(pos.Z());

  // Using ROOT's built-in trilinear interpolation for both components
  TVector3 scField(
      m_primaryFieldR->Interpolate(phi, r, z),
      m_primaryFieldP->Interpolate(phi, r, z),
      m_primaryFieldZ->Interpolate(phi, r, z)
  );
  
  TVector3 ibfField(
      m_totalIbfFieldR->Interpolate(phi, r, z),
      m_totalIbfFieldP->Interpolate(phi, r, z),
      m_totalIbfFieldZ->Interpolate(phi, r, z)
  );

  // Rotate the interpolated cylindrical components to the local phi position
  scField.RotateZ(phi);
  ibfField.RotateZ(phi);
  
  totalField += scField * m_primaryScale;
  totalField += ibfField * m_ibfGain;

  return totalField;
}

TH3F* TunableField::makeUnguardedStandardTH3F(const std::string& name, const std::string& title)
{
  double dphi = (2.0 * M_PI) / m_nphi;
  double dr = (m_rmax - m_rmin) / m_nr;
  double dz = m_zmax / m_nz;

  std::string fullTitle = title + ";#phi (rad);r (cm);z (cm)";
  return new TH3F(name.c_str(), fullTitle.c_str(), m_nphi, 0, 2 * M_PI, m_nr, m_rmin, m_rmax, m_nz, 0, m_zmax);
}

TH3F* TunableField::makeGuardedStandardTH3F(const std::string& name, const std::string& title)
{
  double dphi = (2.0 * M_PI) / m_nphi;
  double dr = (m_rmax - m_rmin) / m_nr;
  double dz = m_zmax / m_nz;
  // (adding 2 to each dimension for the guard bins...)

  std::string fullTitle = title + ";#phi (rad);r (cm);z (cm)";
  return new TH3F(name.c_str(), fullTitle.c_str(), m_nphi + 2, -dphi, 2 * M_PI + dphi, m_nr + 2, m_rmin - dr, m_rmax + dr, m_nz + 2, -dz, m_zmax + dz);
}
