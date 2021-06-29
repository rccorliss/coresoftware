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
  int getSearchResult(double xcheck, double ycheck); // check if coords are in a stripe
  int getStripeID(double xcheck, double ycheck);
 
  int fullID;
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
  
  void CalculateVertices(int nStripes, int nPads, double R[], double spacing[], double x1a[][nRadii], double y1a[][nRadii], double x1b[][nRadii], double y1b[][nRadii], double x2a[][nRadii], double y2a[][nRadii], double x2b[][nRadii], double y2b[][nRadii], double x3a[][nRadii], double y3a[][nRadii], double x3b[][nRadii], double y3b[][nRadii], double padfrac, double str_width[][nRadii], double widthmod[], int nGoodStripes[], int keepUntil[], int nStripesIn[], int nStripesBefore[]);

  PHG4Hitv1* GetBotVerticesFromStripe(int moduleID, int radiusID, int stripeID);
  PHG4Hitv1* GetTopVerticesFromStripe(int moduleID, int radiusID, int stripeID);
  
  int SearchModule(int nStripes, double x1a[][nRadii], double x1b[][nRadii], double x2a[][nRadii], double x2b[][nRadii], double y1a[][nRadii], double y1b[][nRadii], double y2a[][nRadii], double y2b[][nRadii], double x3a[][nRadii], double y3a[][nRadii], double x3b[][nRadii], double y3b[][nRadii], double x, double y, int nGoodStripes[]);
  
  PHG4Hitv1* GetPHG4HitFromStripe(int petalID, int moduleID, int radiusID, int stripeID, int nElectrons);
};


#endif
