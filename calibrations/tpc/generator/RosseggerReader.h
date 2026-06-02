#ifndef ROSSEGGERREADER_H
#define ROSSEGGERREADER_H

#include <TVector3.h>
#include <string>
#include <vector>

class TH3F;

/**
 * \class RosseggerReader
 * \brief Helper class to load and interpolate Rossegger Green's functions from a lookup table.
 */
class RosseggerReader
{
 public:
  RosseggerReader();
  ~RosseggerReader();

  void Load(const std::string& filename);
  
  /**
   * Performs a 5D interpolation of the Green's function.
   * Source interpolation is handled by TH3F::Interpolate.
   * Target interpolation is handled manually between target bins.
   */
  TVector3 GetInterpolatedField(float phi_src, float r_src, float z_src, float r_tgt, float z_tgt);

 private:
  void ClearLookup();
  TH3F* makeUnguardedStandardTH3F(const std::string& name, const std::string& title);

  int m_nr; float m_rmin; float m_rmax; int m_nphi; int m_nz; float m_zmax;
  float m_dr; float m_dz;

  std::vector<TH3F*> m_greenFieldR;
  std::vector<TH3F*> m_greenFieldP;
  std::vector<TH3F*> m_greenFieldZ;
};

#endif // ROSSEGGERREADER_H
