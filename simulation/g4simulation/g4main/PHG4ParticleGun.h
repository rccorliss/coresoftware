// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4PARTICLEGUN_H
#define G4MAIN_PHG4PARTICLEGUN_H

#include "PHG4ParticleGeneratorBase.h"

#include <string>  // for string

class PHCompositeNode;

class PHG4ParticleGun : public PHG4ParticleGeneratorBase
{
 public:
  PHG4ParticleGun(const std::string &name = "PGUN");
  virtual ~PHG4ParticleGun() {}
  int process_event(PHCompositeNode *topNode);
};

#endif
