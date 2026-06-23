#include "FieldMapBuilder.h"
#include "RosseggerReader.h"

#include <TH3.h>
#include <TVector3.h>

#include <cmath>
#include <cstdio>

FieldMapBuilder::FieldMapBuilder(RosseggerReader* rossegger, int nphi, int nr, int nz, float rmin, float rmax, float zmax)
  : m_rossegger(rossegger)
  , m_nphi(nphi)
  , m_nr(nr)
  , m_nz(nz)
  , m_rmin(rmin)
  , m_rmax(rmax)
  , m_zmax(zmax)
{
}

FieldMapBuilder::generateFieldHelpers()
{
  

}


float FieldMapBuilder::computeTotalChargeFromSamples(const std::vector<std::array<float,4>>& samples) const
{
  double total = 0.0;
  for (const auto &s : samples) total += s[3];
  return (float)total;
}

void FieldMapBuilder::computeFieldFromSamples(const std::vector<std::array<float,4>>& samples, TH3F* outR, TH3F* outP, TH3F* outZ)
{
  if (!outR || !outP || !outZ) return;

  const double eps0 = 8.854e-14; // C / (V * cm)
  const double epsinv = 1.0 / eps0;

  // iterate over target (r,z) grid and compute field once per (r,z), then copy into all phi bins
  for (int ir = 0; ir < m_nr; ++ir)
  {
    double field_r = outR->GetYaxis()->GetBinCenter(ir + 1);
    for (int iz = 0; iz < m_nz; ++iz)
    {
      double field_z = outR->GetZaxis()->GetBinCenter(iz + 1);
      TVector3 field(0,0,0);

      for (const auto &s : samples)
      {
        float source_phi = s[0];
        float source_r = s[1];
        float source_z = s[2];
        float q = s[3];
        if (q == 0) continue;

        TVector3 unitField(0,0,0);
        if (m_rossegger)
        {
          unitField = m_rossegger->GetInterpolatedField(source_phi, source_r, source_z, (float)field_r, (float)field_z);
        }
        field += unitField * q;
      }

      field *= epsinv;

      // store into phi bins (mirror TunableField behavior: interior bins start at +2)
      for (int ip = 0; ip < m_nphi; ++ip)
      {
        outR->SetBinContent(ip + 2, ir + 2, iz + 2, field.X());
        outP->SetBinContent(ip + 2, ir + 2, iz + 2, field.Y());
        outZ->SetBinContent(ip + 2, ir + 2, iz + 2, field.Z());
      }
    }
  }

  FillGuardBins(outR);
  FillGuardBins(outP);
  FillGuardBins(outZ);
}

void FieldMapBuilder::normalizeFieldByTotalCharge(TH3F* fieldR, TH3F* fieldP, TH3F* fieldZ, float totalCharge)
{
  if (!fieldR || !fieldP || !fieldZ) return;
  if (std::abs(totalCharge) < 1e-30) return;

  int nP = fieldR->GetNbinsX();
  int nR = fieldR->GetNbinsY();
  int nZ = fieldR->GetNbinsZ();

  for (int ip = 1; ip <= nP; ++ip)
    for (int ir = 1; ir <= nR; ++ir)
      for (int iz = 1; iz <= nZ; ++iz)
      {
        fieldR->SetBinContent(ip, ir, iz, fieldR->GetBinContent(ip, ir, iz) / totalCharge);
        fieldP->SetBinContent(ip, ir, iz, fieldP->GetBinContent(ip, ir, iz) / totalCharge);
        fieldZ->SetBinContent(ip, ir, iz, fieldZ->GetBinContent(ip, ir, iz) / totalCharge);
      }
}

void FieldMapBuilder::FillGuardBins(TH3F* h) const
{
  int nP = h->GetNbinsX();
  int nR = h->GetNbinsY();
  int nZ = h->GetNbinsZ();

  for (int ip = 1; ip <= nP; ++ip)
  {
    for (int ir = 1; ir <= nR; ++ir)
    {
      h->SetBinContent(ip, ir, 1, h->GetBinContent(ip, ir, 2));
      h->SetBinContent(ip, ir, nZ, h->GetBinContent(ip, ir, nZ - 1));
    }
  }
  for (int ip = 1; ip <= nP; ++ip)
  {
    for (int iz = 1; iz <= nZ; ++iz)
    {
      h->SetBinContent(ip, 1, iz, h->GetBinContent(ip, 2, iz));
      h->SetBinContent(ip, nR, iz, h->GetBinContent(ip, nR - 1, iz));
    }
  }
  for (int ir = 1; ir <= nR; ++ir)
  {
    for (int iz = 1; iz <= nZ; ++iz)
    {
      h->SetBinContent(1, ir, iz, h->GetBinContent(nP - 1, ir, iz));
      h->SetBinContent(nP, ir, iz, h->GetBinContent(2, ir, iz));
    }
  }
}
