#ifndef PHG4TPCDIRECTLASER_H
#define PHG4TPCDIRECTLASER_H

#include <iostream>
#include <cmath>
#include <vector>
#include "TMath.h"
#include "TVector3.h"

//from phg4tpcsteppingaction.cc
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Hitv1.h>
//R__LOAD_LIBRARY(libphg4hit.so)


// all distances in mm, all angles in rad
// class that generates stripes and dummy hit coordinates
// stripes have width of one mm, length of one pad width, and are centered in middle of sector gaps

using namespace std;

class PHG4TpcDirectLaser {
public:
  PHG4TpcDirectLaser(); //default constructor

  double begin_CM, end_CM; // inner and outer radii of central membrane
  double ifc,ofc;
  
  vector<PHG4Hitv1*> PHG4Hits;
  
private:
  static const int nLasers = 4;
  const double mm = 1.0;
  const double cm = 10.0;

  TVector3 GetCmStrike(TVector3 start, TVector3 direction);
  TVector3 GetFieldcageStrike(TVector3 start, TVector3 direction);
  TVector3 GetCylinderStrike(TVector3 s, TVector3 v, float radius)
  
  int nElectrons;
 
  PHG4Hitv1* GenerateLaserHit(float theta, float phi, int laser);
};


#endif
