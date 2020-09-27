// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef PHMICROMEGASTPCTRACKMATCHING_H
#define PHMICROMEGASTPCTRACKMATCHING_H

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>

#include <string>
#include <set>
#include <vector>

class PHCompositeNode;
class SvtxTrackMap;
class SvtxTrack;
class SvtxVertexMap;
class TrkrClusterContainer;
class TrkrCluster;
class TrkrClusterHitAssoc;
class TrkrHitTruthAssoc;
class PHG4TruthInfoContainer;
class PHG4HitContainer;
class PHG4Particle;
class AssocInfoContainer;
class TF1;

#include <trackreco/PHTrackPropagating.h>

class PHMicromegasTpcTrackMatching : public PHTrackPropagating
{
 public:

  PHMicromegasTpcTrackMatching(const std::string &name = "PHMicromegasTpcTrackMatching");

  virtual ~PHMicromegasTpcTrackMatching();

  void set_rphi_search_window_lyr1(const double win){_rphi_search_win_1 = win;}
  void set_z_search_window_lyr1(const double win){_z_search_win_1 = win;}
  void set_rphi_search_window_lyr2(const double win){_rphi_search_win_2 = win;}
  void set_z_search_window_lyr2(const double win){_z_search_win_2 = win;}

  void set_field_dir(const double rescale)
  {
    _fieldDir = -1;
    if(rescale > 0)
      _fieldDir = 1;     
  }
  void set_field(const std::string field) { _field = field;}

 protected:
  int Setup(PHCompositeNode* topNode) override;

  int Process() override;

  int End() override;
  
 private:

  int GetNodes(PHCompositeNode* topNode);

  float line_circle_intersection(float x[], float y[], float z[], float radius);
  void CircleFitByTaubin (std::vector<TrkrCluster*> clusters, double &R, double &X0, double &Y0);
  void circle_circle_intersection(double r1, double r2, double x2, double y2, double &xplus, double &yplus, double &xminus, double &yminus);
  
  // default values, can be replaced from the macro, all in cm
  double _rphi_search_win_1 = 0.20; 
  double _z_search_win_1 = 13.0;
  double _rphi_search_win_2 = 13.0;  
  double _z_search_win_2 = 0.25;

  unsigned int _n_mm_layers = 2;
  double _mm_layer_radius[2] = { 82.2565, 82.6998};
  double _xplus[2] = {0, 0};
  double _yplus[2] = {0, 0};
  double _xminus[2] = {0, 0};
  double _yminus[2] = {0, 0};

  double _min_tpc_layer = 45;
  double _min_mm_layer = 55;
  
  SvtxTrack *_tracklet_tpc{nullptr};

  std::string _field;
  int _fieldDir = -1;

};

#endif // PHMICROMEGASTPCTRACKMATCHING_H
