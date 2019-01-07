// Tell emacs that this is a C++ source
// This file is really -*- C++ -*-.
#ifndef G4INTT_PHG4INTTSUBSYSTEM_H
#define G4INTT_PHG4INTTSUBSYSTEM_H

#include <g4detectors/PHG4DetectorGroupSubsystem.h>

#include <vector>

class PHG4INTTDetector;
class PHG4INTTSteppingAction;

class PHG4INTTSubsystem : public PHG4DetectorGroupSubsystem
{
 public:
  typedef std::vector<std::pair<int, int>> vpair;

  //! constructor
  PHG4INTTSubsystem(const std::string &name = "INTT", const vpair &layerconfig = vpair(0));

  //! destructor
  virtual ~PHG4INTTSubsystem(void)
  {
  }

  //! init
  /*!
  creates the detector_ object and place it on the node tree, under "DETECTORS" node (or whatever)
  reates the stepping action and place it on the node tree, under "ACTIONS" node
  creates relevant hit nodes that will be populated by the stepping action and stored in the output DST
  */
  int InitRunSubsystem(PHCompositeNode *);

  //! event processing
  /*!
  get all relevant nodes from top nodes (namely hit list)
  and pass that to the stepping action
  */
  int process_event(PHCompositeNode *);

  //! accessors (reimplemented)
  PHG4Detector *GetDetector(void) const;
  PHG4SteppingAction *GetSteppingAction(void) const { return m_SteppingAction; }
  void Print(const std::string &what = "ALL") const;

 private:
  void SetDefaultParameters();

  //! detector geometry
  /*! defives from PHG4Detector */
  PHG4INTTDetector *m_Detector;

  //! particle tracking "stepping" action
  /*! derives from PHG4SteppingActions */
  PHG4SteppingAction *m_SteppingAction;

  vpair m_LayerConfigVector;
  std::string m_DetectorType;
};

#endif
