#ifndef FIELDMAPBUILDER_H
#define FIELDMAPBUILDER_H

#include <vector>
#include <array>

class RosseggerReader;
class TH3F;

class FieldMapBuilder
{
public:
  FieldMapBuilder(RosseggerReader* rossegger, int nphi, int nr, int nz, float rmin, float rmax, float zmax);

  // samples: each entry {phi, r, z, q}
  void computeFieldFromSamples(const std::vector<std::array<float,4>>& samples, TH3F* outR, TH3F* outP, TH3F* outZ);

  void normalizeFieldByTotalCharge(TH3F* fieldR, TH3F* fieldP, TH3F* fieldZ, float totalCharge);

  float computeTotalChargeFromSamples(const std::vector<std::array<float,4>>& samples) const;

private:
  RosseggerReader* m_rossegger;
  int m_nphi;
  int m_nr;
  int m_nz;
  float m_rmin;
  float m_rmax;
  float m_zmax;

  void FillGuardBins(TH3F* h) const;
};

#endif // FIELDMAPBUILDER_H
