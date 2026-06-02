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
}

int PadrowReader::getPadrowFromR(float r) const
{
    // For now, return 0
    return 0;
}

void PadrowReader::getPadrowRBounds(int padrow_idx, float& rmin, float& rmax) const
{
    // For now, set rmin and rmax to 0... we gotta fix this from the PHG4TpcCylindericalGeomContainer eventually.
    rmin = 0.0f;
    rmax = 0.0f;
}