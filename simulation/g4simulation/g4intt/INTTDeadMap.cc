#include "INTTDeadMap.h"

#include <cassert>
#include <iostream>
#include <limits>

using namespace std;

int
    INTTDeadMap::s_wildCardID = -1;

const INTTDeadMap::Map&
INTTDeadMap::getDeadChannels(void) const
{
  static Map tmp_map;
  return tmp_map;
}

INTTDeadMap::Map&
INTTDeadMap::getDeadChannels(void)
{
  static Map tmp_map;
  return tmp_map;
}

void INTTDeadMap::addDeadChannel(const int layer, const int ieta, const int iphi)
{
  const PHG4CellDefs::keytype key = PHG4CellDefs::EtaPhiBinning::genkey(layer, ieta, iphi);
  addDeadChannel(key);
}

void INTTDeadMap::addDeadChannel(PHG4CellDefs::keytype key)
{
}

void INTTDeadMap::addDeadChannelINTT(const int layer,
                                     const int ladder_phi, const int ladder_z,
                                     const int strip_z, const int strip_phi)
{
  addDeadChannel(getINTTKey(layer,
                            ladder_phi, ladder_z,
                            strip_z, strip_phi));
}

bool INTTDeadMap::isDeadChannel(PHG4CellDefs::keytype key) const
{
  return false;
}

bool INTTDeadMap::isDeadChannel(const int layer, const int ieta, const int iphi) const
{
  const PHG4CellDefs::keytype key = PHG4CellDefs::EtaPhiBinning::genkey(layer, ieta, iphi);
  return isDeadChannel(key);
}

bool INTTDeadMap::isDeadChannelINTT(const int layer,
                                    const int ladder_phi, const int ladder_z,
                                    const int strip_z, const int strip_phi) const
{
  if (isDeadChannel(getINTTKey(layer,
                               ladder_phi, ladder_z,
                               strip_z, strip_phi)))
    return true;
  else if (isDeadChannel(getINTTKey(layer,
                                    ladder_phi, ladder_z,
                                    strip_z, s_wildCardID)))
    return true;
  else if (isDeadChannel(getINTTKey(layer,
                                    ladder_phi, ladder_z,
                                    s_wildCardID, s_wildCardID)))
    return true;
  else if (isDeadChannel(getINTTKey(layer,
                                    ladder_phi, s_wildCardID,
                                    s_wildCardID, s_wildCardID)))
    return true;
  else if (isDeadChannel(getINTTKey(layer,
                                    s_wildCardID, s_wildCardID,
                                    s_wildCardID, s_wildCardID)))
    return true;
  else
    return false;
}

int INTTDeadMap::isValid() const
{
  return size() > 0;
}

void INTTDeadMap::Reset()
{
}

void INTTDeadMap::identify(std::ostream& os) const
{
  os << "INTTDeadMap" << std::endl;
}

PHG4CellDefs::keytype INTTDeadMap::getINTTKey(int layer,
                                              int ladder_phi, int ladder_z,
                                              int strip_z, int strip_phi)
{
  static const int layer_bit = 8;
  static const int ladder_phi_bit = 16;
  static const int ladder_z_bit = 8;
  static const int strip_z_bit = 16;
  static const int strip_phi_bit = 16;

  bool wildcard = false;

  if (layer == s_wildCardID) wildcard = true;
  if (wildcard) layer = (1 << layer_bit) - 1;

  if (ladder_phi == s_wildCardID) wildcard = true;
  if (wildcard) ladder_phi = (1 << ladder_phi_bit) - 1;

  if (ladder_z == s_wildCardID) wildcard = true;
  if (wildcard) ladder_z = (1 << ladder_z_bit) - 1;

  if (strip_z == s_wildCardID) wildcard = true;
  if (wildcard) strip_z = (1 << strip_z_bit) - 1;

  if (strip_phi == s_wildCardID) wildcard = true;
  if (wildcard) strip_phi = (1 << strip_phi_bit) - 1;

  //  bit sum check
  assert(layer_bit + ladder_phi_bit + ladder_z_bit + strip_z_bit + strip_phi_bit == numeric_limits<PHG4CellDefs::keytype>::digits);

  //max range check
  assert(layer < (1 << layer_bit));
  assert(ladder_phi < (1 << ladder_phi_bit));
  assert(ladder_z < (1 << ladder_z_bit));
  assert(strip_z < (1 << strip_z_bit));
  assert(strip_phi < (1 << strip_phi_bit));

  //min range check
  assert(layer >= 0);
  assert(ladder_phi >= 0);
  assert(ladder_z >= 0);
  assert(strip_z >= 0);
  assert(strip_phi >= 0);

  PHG4CellDefs::keytype key = 0;

  key += layer;

  key <<= ladder_phi_bit;
  key += ladder_phi;

  key <<= ladder_z_bit;
  key += ladder_z;

  key <<= strip_z_bit;
  key += strip_z;

  key <<= strip_phi_bit;
  key += strip_phi;

  return key;
}
