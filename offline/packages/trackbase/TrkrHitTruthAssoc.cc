/**
 * @file trackbase/TrkrHitTruthAssoc.cc
 * @author D. McGlinchey
 * @date June 2018
 * @brief Implementation of TrkrHitTruthAssoc
 */

#include "TrkrHitTruthAssoc.h"
#include "TrkrDefs.h"

#include <g4main/PHG4HitDefs.h>

#include <algorithm>
#include <ostream>               // for operator<<, endl, basic_ostream, bas...

TrkrHitTruthAssoc::TrkrHitTruthAssoc()
  : m_map()
{
}

TrkrHitTruthAssoc::~TrkrHitTruthAssoc()
{
}

void
TrkrHitTruthAssoc::Reset()
{
m_map.clear();
}

void
TrkrHitTruthAssoc::identify(std::ostream &os) const
{
  os << "-----TrkrHitTruthAssoc-----" << std::endl;
  os << "Number of associations: " << m_map.size() << std::endl;

  for ( auto& entry : m_map )
  {
    int layer = TrkrDefs::getLayer(entry.first);
    os << "   hitset key: "  << entry.first << " layer " << layer
       << " hit key: " << entry.second.first
       << " g4hit key: " << entry.second.second
       << std::endl;
  }

  os << "------------------------------" << std::endl;

  return;
}

void
TrkrHitTruthAssoc::addAssoc(const TrkrDefs::hitsetkey hitsetkey, const TrkrDefs::hitkey hitkey, const PHG4HitDefs::keytype g4hitkey)
{
  // the association we want is between TrkrHit and PHG4Hit, but we need to know which TrkrHitSet the TrkrHit is in
  m_map.insert (std::make_pair(hitsetkey, std::make_pair(hitkey, g4hitkey)));
}

void
TrkrHitTruthAssoc::findOrAddAssoc(const TrkrDefs::hitsetkey hitsetkey, const TrkrDefs::hitkey hitkey, const PHG4HitDefs::keytype g4hitkey)
{
  // the association we want is between TrkrHit and PHG4Hit, but we need to know which TrkrHitSet the TrkrHit is in
  // check if this association already exists
  // We need all hitsets with this key

  const auto hitsetrange = m_map.equal_range(hitsetkey);
  if( std::any_of(
    hitsetrange.first,
    hitsetrange.second,
    [&hitkey, &g4hitkey]( const MMap::value_type& pair ) { return pair.second.first == hitkey && pair.second.second == g4hitkey; } ) )
  { return; }

  // Does not exist, create it
  const auto assoc = std::make_pair(hitkey, g4hitkey);
  m_map.insert (hitsetrange.second, std::make_pair(hitsetkey, std::move(assoc)));
}

void
TrkrHitTruthAssoc::removeAssoc(const TrkrDefs::hitsetkey hitsetkey, const TrkrDefs::hitkey hitkey)
{
  // remove all entries for this TrkrHit and its PHG4Hits, but we need to know which TrkrHitSet the TrkrHit is in
  // check if this association already exists
  // We need all hitsets with this key

  std::pair<MMap::iterator, MMap::iterator> hitsetrange = m_map.equal_range(hitsetkey);
  MMap::iterator mapiter = hitsetrange.first;
  for (mapiter = hitsetrange.first; mapiter != hitsetrange.second; ++mapiter)
    {
      if(mapiter->second.first == hitkey)
  {
    // exists, erase it
    //std::cout << "Association with hitsetkey " << hitsetkey << " hitkey " << hitkey
    //	    << " g4hitkey " << mapiter->second.second << "  exists, erasing it " << std::endl;
    m_map.erase (mapiter);
    return;
  }
    }

}

void
//TrkrHitTruthAssoc::ConstRange
TrkrHitTruthAssoc::getG4Hits(const TrkrDefs::hitsetkey hitsetkey, const unsigned int hidx, MMap &temp_map)
{
  //MMap temp_map;

  std::pair<MMap::iterator, MMap::iterator> hitsetrange = m_map.equal_range(hitsetkey);
  MMap::iterator mapiter = hitsetrange.first;

  for (mapiter = hitsetrange.first; mapiter != hitsetrange.second; ++mapiter)
    {
      if(mapiter->second.first == hidx)
  {
    //std::cout << "   Association with hitsetkey " << hitsetkey << " hitkey " << mapiter->second.first << " g4hitkey " << mapiter->second.second << " exists " << std::endl;
    // add this to the return object
    std::pair<TrkrDefs::hitkey, PHG4HitDefs::keytype> assoc = std::make_pair(mapiter->second.first, mapiter->second.second);
    temp_map.insert (std::make_pair(hitsetkey, assoc));
  }
    }

  //return std::make_pair(temp_map.begin(), temp_map.end());
}

