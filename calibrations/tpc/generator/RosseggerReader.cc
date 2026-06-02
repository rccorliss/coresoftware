#include "RosseggerReader.h"

#include <TH3.h>
#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <format>
#include <iostream>

RosseggerReader::RosseggerReader()
  : m_nr(0)
  , m_rmin(0.0)
  , m_rmax(0.0)
  , m_nphi(0)
  , m_nz(0)
  , m_zmax(0.0)
  , m_dr(0.0)
  , m_dz(0.0)
{
}

RosseggerReader::~RosseggerReader()
{
  ClearLookup();
}

void RosseggerReader::ClearLookup()
{
  for (TH3F* h : m_greenFieldR) delete h;
  for (TH3F* h : m_greenFieldP) delete h;
  for (TH3F* h : m_greenFieldZ) delete h;
  m_greenFieldR.clear(); m_greenFieldP.clear(); m_greenFieldZ.clear();
}

TH3F* RosseggerReader::makeUnguardedStandardTH3F(const std::string& name, const std::string& title)
{
  std::string fullTitle = title + ";#phi (rad);r (cm);z (cm)";
  return new TH3F(name.c_str(), fullTitle.c_str(), m_nphi, 0, 2.0 * M_PI, m_nr, m_rmin, m_rmax, m_nz, 0, m_zmax);
}

void RosseggerReader::Load(const std::string& filename)
{
  TFile* input = TFile::Open(filename.c_str(), "READ");
  if (!input)
  {
    std::cout << "RosseggerReader::Load: could not open " << filename << std::endl;
    return;
  }

  TTree* tInfo = (TTree*) input->Get("info");
  if (!tInfo)
  {
    std::cout << "RosseggerReader::Load: 'info' tree not found" << std::endl;
    input->Close();
    return;
  }

  float f_rmin, f_rmax, f_zmin, f_zmax;
  int f_nr, f_np, f_nz;
  tInfo->SetBranchAddress("rmin", &f_rmin);
  tInfo->SetBranchAddress("rmax", &f_rmax);
  tInfo->SetBranchAddress("zmin", &f_zmin);
  tInfo->SetBranchAddress("zmax", &f_zmax);
  tInfo->SetBranchAddress("nr", &f_nr);
  tInfo->SetBranchAddress("nphi", &f_np);
  tInfo->SetBranchAddress("nz", &f_nz);
  tInfo->GetEntry(0);

  m_rmin = f_rmin;
  m_rmax = f_rmax;
  m_zmax = f_zmax;
  m_nr = f_nr;
  m_nphi = f_np;
  m_nz = f_nz;
  m_dr = (m_rmax - m_rmin) / (float) m_nr;
  m_dz = m_zmax / (float) m_nz;

  TTree* tLookup = (TTree*) input->Get("phislice");
  if (!tLookup) { input->Close(); return; }

  int ior, ifr, iophi, ioz, ifz;
  TVector3* unitf = nullptr;
  tLookup->SetBranchAddress("ir_source", &ior);
  tLookup->SetBranchAddress("ir_target", &ifr);
  tLookup->SetBranchAddress("ip_source", &iophi);
  tLookup->SetBranchAddress("iz_source", &ioz);
  tLookup->SetBranchAddress("iz_target", &ifz);
  tLookup->SetBranchAddress("Evec", &unitf);

  ClearLookup();
  m_greenFieldR.assign(m_nr * m_nz, nullptr);
  m_greenFieldP.assign(m_nr * m_nz, nullptr);
  m_greenFieldZ.assign(m_nr * m_nz, nullptr);

  long int nEntries = tLookup->GetEntries();
  for (long int i = 0; i < nEntries; i++)
  {
    tLookup->GetEntry(i);
    int target_idx = ifr + ifz * m_nr;
    
    if (m_greenFieldR[target_idx] == nullptr)
    {
      m_greenFieldR[target_idx] = makeUnguardedStandardTH3F(std::format("h_rR_{}_{}", ifr, ifz), "R");
      m_greenFieldP[target_idx] = makeUnguardedStandardTH3F(std::format("h_rP_{}_{}", ifr, ifz), "P");
      m_greenFieldZ[target_idx] = makeUnguardedStandardTH3F(std::format("h_rZ_{}_{}", ifr, ifz), "Z");
    }

    m_greenFieldR[target_idx]->SetBinContent(iophi + 1, ior + 1, ioz + 1, (float)unitf->X() * -1.0f);
    m_greenFieldP[target_idx]->SetBinContent(iophi + 1, ior + 1, ioz + 1, (float)unitf->Y() * -1.0f);
    m_greenFieldZ[target_idx]->SetBinContent(iophi + 1, ior + 1, ioz + 1, (float)unitf->Z() * -1.0f);
  }

  input->Close();
}

TVector3 RosseggerReader::GetInterpolatedField(float phi_src, float r_src, float z_src, float r_tgt, float z_tgt)
{
  // 1. Target interpolation weights (bilinear in r-z plane)
  float r_target_bin_f = (r_tgt - m_rmin) / m_dr - 0.5f;
  int ir0 = (int) std::floor(r_target_bin_f);
  int ir1 = ir0 + 1;
  float wr = r_target_bin_f - (float) ir0;

  float z_target_bin_f = (z_tgt - 0.0f) / m_dz - 0.5f;
  int iz0 = (int) std::floor(z_target_bin_f);
  int iz1 = iz0 + 1;
  float wz = z_target_bin_f - (float) iz0;

  if (ir0 < 0 || ir1 >= m_nr || iz0 < 0 || iz1 >= m_nz)
  {
    std::cout << "RosseggerReader::GetInterpolatedField: target out of bounds (r_tgt=" << r_tgt << ", z_tgt=" << z_tgt << ")" << " produced ir0=" << ir0 << ", ir1=" << ir1 << ", iz0=" << iz0 << ", iz1=" << iz1 << " (max=" << m_nr - 1 << ", " << m_nz - 1 << ")" << std::endl;
    return TVector3(0, 0, 0);
  }

  // Clamp target indices to grid range
  ir0 = std::max(0, std::min(m_nr - 1, ir0));
  ir1 = std::max(0, std::min(m_nr - 1, ir1));
  iz0 = std::max(0, std::min(m_nz - 1, iz0));
  iz1 = std::max(0, std::min(m_nz - 1, iz1));

  int idx[4] = {
      ir0 + iz0 * m_nr,
      ir1 + iz0 * m_nr,
      ir0 + iz1 * m_nr,
      ir1 + iz1 * m_nr
  };

  TVector3 fields[4];
  for (int i = 0; i < 4; i++)
  {
    if (m_greenFieldR[idx[i]])
    {
      fields[i].SetXYZ(
          m_greenFieldR[idx[i]]->Interpolate(phi_src, r_src, z_src),
          m_greenFieldP[idx[i]]->Interpolate(phi_src, r_src, z_src),
          m_greenFieldZ[idx[i]]->Interpolate(phi_src, r_src, z_src)
      );
    }
  }

  TVector3 z0field = fields[0] * (1.0f - wr) + fields[1] * wr;
  TVector3 z1field = fields[2] * (1.0f - wr) + fields[3] * wr;
  return z0field * (1.0f - wz) + z1field * wz;
}
