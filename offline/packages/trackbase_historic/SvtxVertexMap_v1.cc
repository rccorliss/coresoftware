#include "SvtxVertexMap_v1.h"

#include "SvtxVertex.h"

using namespace std;

ClassImp(SvtxVertexMap_v1)

    SvtxVertexMap_v1::SvtxVertexMap_v1()
  : _map()
{
}

SvtxVertexMap_v1::SvtxVertexMap_v1(const SvtxVertexMap_v1& vertexmap)
  : _map()
{
  for (ConstIter iter = vertexmap.begin();
       iter != vertexmap.end();
       ++iter)
  {
    const SvtxVertex* vertex = iter->second;
    _map.insert(make_pair(vertex->get_id(), vertex->Clone()));
  }
}

SvtxVertexMap_v1& SvtxVertexMap_v1::operator=(const SvtxVertexMap_v1& vertexmap)
{
  Reset();
  for (ConstIter iter = vertexmap.begin();
       iter != vertexmap.end();
       ++iter)
  {
    const SvtxVertex* vertex = iter->second;
    _map.insert(make_pair(vertex->get_id(), vertex->Clone()));
  }
  return *this;
}

SvtxVertexMap_v1::~SvtxVertexMap_v1()
{
  Reset();
}

void SvtxVertexMap_v1::Reset()
{
  for (Iter iter = _map.begin();
       iter != _map.end();
       ++iter)
  {
    SvtxVertex* vertex = iter->second;
    delete vertex;
  }
  _map.clear();
}

void SvtxVertexMap_v1::identify(ostream& os) const
{
  os << "SvtxVertexMap_v1: size = " << _map.size() << endl;
  return;
}

const SvtxVertex* SvtxVertexMap_v1::get(unsigned int id) const
{
  ConstIter iter = _map.find(id);
  if (iter == _map.end()) return NULL;
  return iter->second;
}

SvtxVertex* SvtxVertexMap_v1::get(unsigned int id)
{
  Iter iter = _map.find(id);
  if (iter == _map.end()) return NULL;
  return iter->second;
}

SvtxVertex* SvtxVertexMap_v1::insert(const SvtxVertex* vertex)
{
  unsigned int index = 0;
  if (!_map.empty()) index = _map.rbegin()->first + 1;
  _map.insert(make_pair(index, vertex->Clone()));
  _map[index]->set_id(index);
  return _map[index];
}
