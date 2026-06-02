#ifndef PADROWREADER_H
#define PADROWREADER_H

#include <vector>
#include <utility> // For std::pair

/**
 * \class PadrowReader
 * \brief Helper class to provide information about TPC padrow geometry and active areas.
 */
class PadrowReader
{
public:
    PadrowReader();
    ~PadrowReader();

    void setPadrowBoundaries(const std::vector<std::pair<float, float>>& padrow_boundaries);

    // Returns the fraction of the given region that is in an active area.
    float getActiveFraction(float rmin, float rmax, float phimin, float phimax) const;

    // Returns the padrow index (1-48) for a given radial position r.
    int getPadrowFromR(float r) const;

    // For a given padrow index, it returns its rmin and rmax.
    void getPadrowRBounds(int padrow_idx, float& rmin, float& rmax) const;

private:
    std::vector<std::pair<float, float>> m_padrowBoundaries;
};

#endif // PADROWREADER_H