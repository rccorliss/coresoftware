#include "PadrowReader.h"

#include <iostream>

PadrowReader::PadrowReader()
{
    // Constructor implementation
}

PadrowReader::~PadrowReader()
{
    // Destructor implementation
}

void PadrowReader::setPadrowBoundaries(const std::vector<std::pair<float, float>>& padrow_boundaries)
{
    m_padrowBoundaries = padrow_boundaries;
}

float PadrowReader::getActiveFraction(float rmin, float rmax, float phimin, float phimax) const
{
    // For now, return 1.0 as per request
    return 1.0f;
  // Simplified: check if the midpoint of the r-range falls within any padrow
  float rmid = 0.5 * (rmin + rmax);
  for (const std::pair<float, float>& bounds : m_padrowBoundaries)
  {
    if (rmid >= bounds.first && rmid < bounds.second) return 1.0f;
  }
  return 0.0f;
}

int PadrowReader::getPadrowFromR(float r) const
{
    // For now, return 0
    return 0;
  for (size_t i = 0; i < m_padrowBoundaries.size(); ++i)
  {
    if (r >= m_padrowBoundaries[i].first && r < m_padrowBoundaries[i].second) return (int)i;
  }
  return -1;
}

void PadrowReader::getPadrowRBounds(int padrow_idx, float& rmin, float& rmax) const
{
    // For now, set rmin and rmax to 0... we gotta fix this from the PHG4TpcCylindericalGeomContainer eventually.
  if (padrow_idx >= 0 && padrow_idx < (int)m_padrowBoundaries.size())
  {
    rmin = m_padrowBoundaries[padrow_idx].first;
    rmax = m_padrowBoundaries[padrow_idx].second;
  }
  else
  {
    rmin = 0.0f;
    rmax = 0.0f;
  }
}