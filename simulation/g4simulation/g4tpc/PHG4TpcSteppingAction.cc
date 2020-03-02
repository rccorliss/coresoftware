#include "PHG4TpcSteppingAction.h"
#include "PHG4TpcDetector.h"

#include <g4detectors/PHG4StepStatusDecode.h>

#include <phparameter/PHParameters.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hitv1.h>
#include <g4main/PHG4Shower.h>
#include <g4main/PHG4SteppingAction.h>         // for PHG4SteppingAction

#include <g4main/PHG4TrackUserInfoV1.h>

#include <phool/getClass.h>

#include <Geant4/G4ParticleDefinition.hh>
#include <Geant4/G4Step.hh>
#include <Geant4/G4StepPoint.hh>               // for G4StepPoint
#include <Geant4/G4StepStatus.hh>              // for fGeomBoundary, fPostSt...
#include <Geant4/G4String.hh>                  // for G4String
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>             // for G4ThreeVector
#include <Geant4/G4TouchableHandle.hh>         // for G4TouchableHandle
#include <Geant4/G4Track.hh>                   // for G4Track
#include <Geant4/G4TrackStatus.hh>             // for fStopAndKill
#include <Geant4/G4Types.hh>                   // for G4double
#include <Geant4/G4VPhysicalVolume.hh>         // for G4VPhysicalVolume
#include <Geant4/G4VTouchable.hh>              // for G4VTouchable
#include <Geant4/G4VUserTrackInformation.hh>   // for G4VUserTrackInformation

#include <cmath>                               // for isfinite
#include <cstdlib>                            // for exit
#include <iostream>
#include <string>                              // for operator<<, operator+

class PHCompositeNode;

using namespace std;
//____________________________________________________________________________..
PHG4TpcSteppingAction::PHG4TpcSteppingAction(PHG4TpcDetector* detector, const PHParameters* parameters)
  : PHG4SteppingAction(detector->GetName())
  , detector_(detector)
  , hits_(nullptr)
  , absorberhits_(nullptr)
  , hit(nullptr)
  , params(parameters)
  , savehitcontainer(nullptr)
  , saveshower(nullptr)
  , savevolpre(nullptr)
  , savevolpost(nullptr)
  , savetrackid(-1)
  , saveprestepstatus(-1)
  , savepoststepstatus(-1)
  , IsActive(params->get_int_param("active"))
  , IsBlackHole(params->get_int_param("blackhole"))
  , use_g4_steps(0)
{
  if (isfinite(params->get_double_param("steplimits")))
  {
    use_g4_steps = 1;
  }
  SetName(detector_->GetName());
}

PHG4TpcSteppingAction::~PHG4TpcSteppingAction()
{
  // if the last hit was a zero energie deposit hit, it is just reset
  // and the memory is still allocated, so we need to delete it here
  // if the last hit was saved, hit is a nullptr pointer which are
  // legal to delete (it results in a no operation)
  delete hit;
}
//____________________________________________________________________________..
bool PHG4TpcSteppingAction::UserSteppingAction(const G4Step* aStep, bool)
{
  G4TouchableHandle touch = aStep->GetPreStepPoint()->GetTouchableHandle();
  G4TouchableHandle touchpost = aStep->GetPostStepPoint()->GetTouchableHandle();
  // get volume of the current step
  G4VPhysicalVolume* volume = touch->GetVolume();

  // detector_->IsInTpc(volume)
  // returns
  //  0 is outside of Tpc
  //  1 is inside tpc gas
  // <0 is in tpc support structure (cage, endcaps,...)

  int whichactive = detector_->IsInTpc(volume);
  if (!whichactive)
  {
    return false;
  }
  unsigned int layer_id = 99;  // no layer number for the hit, use a non-existent one for now, replace it later
  // collect energy and track length step by step
  G4double edep = aStep->GetTotalEnergyDeposit() / GeV;
  G4double eion = (aStep->GetTotalEnergyDeposit() - aStep->GetNonIonizingEnergyDeposit()) / GeV;
  const G4Track* aTrack = aStep->GetTrack();

  // if this block stops everything, just put all kinetic energy into edep
  if (IsBlackHole)
  {
    edep = aTrack->GetKineticEnergy() / GeV;
    G4Track* killtrack = const_cast<G4Track*>(aTrack);
    killtrack->SetTrackStatus(fStopAndKill);
  }

  bool geantino = false;

  // the check for the pdg code speeds things up, I do not want to make
  // an expensive string compare for every track when we know
  // geantino or chargedgeantino has pid=0
  if (aTrack->GetParticleDefinition()->GetPDGEncoding() == 0 &&
      aTrack->GetParticleDefinition()->GetParticleName().find("geantino") != string::npos)
  {
    geantino = true;
  }
  G4StepPoint* prePoint = aStep->GetPreStepPoint();
  G4StepPoint* postPoint = aStep->GetPostStepPoint();
  int prepointstatus = prePoint->GetStepStatus();

  //       cout << "track id " << aTrack->GetTrackID() << endl;
  //       cout << "time prepoint: " << prePoint->GetGlobalTime() << endl;
  //       cout << "time postpoint: " << postPoint->GetGlobalTime() << endl;

  if ((use_g4_steps > 0 && whichactive > 0) ||
      prepointstatus == fGeomBoundary ||
      prepointstatus == fUndefined ||
      (prepointstatus == fPostStepDoItProc && savepoststepstatus == fGeomBoundary))
  {
    // this is for debugging weird occurances we have occasionally
    if (prepointstatus == fPostStepDoItProc && savepoststepstatus == fGeomBoundary)
    {
      cout << GetName() << ": New Hit for  " << endl;
      cout << "prestep status: " << PHG4StepStatusDecode::GetStepStatus(prePoint->GetStepStatus())
           << ", poststep status: " << PHG4StepStatusDecode::GetStepStatus(postPoint->GetStepStatus())
           << ", last pre step status: " << PHG4StepStatusDecode::GetStepStatus(saveprestepstatus)
           << ", last post step status: " << PHG4StepStatusDecode::GetStepStatus(savepoststepstatus) << endl;
      cout << "last track: " << savetrackid
           << ", current trackid: " << aTrack->GetTrackID() << endl;
      cout << "phys pre vol: " << volume->GetName()
           << " post vol : " << touchpost->GetVolume()->GetName() << endl;
      cout << " previous phys pre vol: " << savevolpre->GetName()
           << " previous phys post vol: " << savevolpost->GetName() << endl;
    }
    // if previous hit was saved, hit pointer was set to nullptr
    // and we have to make a new one
    if (!hit)
    {
      hit = new PHG4Hitv1();
    }
    hit->set_layer(layer_id);
    //here we set the entrance values in cm
    hit->set_x(0, prePoint->GetPosition().x() / cm);
    hit->set_y(0, prePoint->GetPosition().y() / cm);
    hit->set_z(0, prePoint->GetPosition().z() / cm);

    // momentum
    hit->set_px(0, prePoint->GetMomentum().x() / GeV);
    hit->set_py(0, prePoint->GetMomentum().y() / GeV);
    hit->set_pz(0, prePoint->GetMomentum().z() / GeV);

    // time in ns
    hit->set_t(0, prePoint->GetGlobalTime() / nanosecond);
    //set and save the track ID
    hit->set_trkid(aTrack->GetTrackID());
    savetrackid = aTrack->GetTrackID();
    //set the initial energy deposit
    hit->set_edep(0);
    if (whichactive > 0)  // return of IsInTpcDetector, > 0 hit in tpc gas volume, < 0 hit in support structures
    {
      hit->set_eion(0);
      // Now save the container we want to add this hit to
      savehitcontainer = hits_;
    }
    else
    {
      savehitcontainer = absorberhits_;
    }
    if (G4VUserTrackInformation* p = aTrack->GetUserInformation())
    {
      if (PHG4TrackUserInfoV1* pp = dynamic_cast<PHG4TrackUserInfoV1*>(p))
      {
        hit->set_trkid(pp->GetUserTrackId());
        hit->set_shower_id(pp->GetShower()->get_id());
        saveshower = pp->GetShower();
      }
    }
    // some sanity checks for inconsistencies
    // check if this hit was created, if not print out last post step status
    if (!hit || !isfinite(hit->get_x(0)))
    {
      cout << GetName() << ": hit was not created" << endl;
      cout << "prestep status: " << PHG4StepStatusDecode::GetStepStatus(prePoint->GetStepStatus())
           << ", poststep status: " << PHG4StepStatusDecode::GetStepStatus(postPoint->GetStepStatus())
           << ", last pre step status: " << PHG4StepStatusDecode::GetStepStatus(saveprestepstatus)
           << ", last post step status: " << PHG4StepStatusDecode::GetStepStatus(savepoststepstatus) << endl;
      cout << "last track: " << savetrackid
           << ", current trackid: " << aTrack->GetTrackID() << endl;
      cout << "phys pre vol: " << volume->GetName()
           << " post vol : " << touchpost->GetVolume()->GetName() << endl;
      cout << " previous phys pre vol: " << savevolpre->GetName()
           << " previous phys post vol: " << savevolpost->GetName() << endl;
      exit(1);
    }
    // check if track id matches the initial one when the hit was created
    if (aTrack->GetTrackID() != savetrackid)
    {
      cout << GetName() << ": hits do not belong to the same track" << endl;
      cout << "saved track: " << savetrackid
           << ", current trackid: " << aTrack->GetTrackID()
           << ", prestep status: " << prePoint->GetStepStatus()
           << ", previous post step status: " << savepoststepstatus
           << endl;

      exit(1);
    }
    saveprestepstatus = prePoint->GetStepStatus();
    savepoststepstatus = postPoint->GetStepStatus();
    savevolpre = volume;
    savevolpost = touchpost->GetVolume();
    // here we just update the exit values, it will be overwritten
    // for every step until we leave the volume or the particle
    // ceases to exist
    hit->set_x(1, postPoint->GetPosition().x() / cm);
    hit->set_y(1, postPoint->GetPosition().y() / cm);
    hit->set_z(1, postPoint->GetPosition().z() / cm);

    hit->set_px(1, postPoint->GetMomentum().x() / GeV);
    hit->set_py(1, postPoint->GetMomentum().y() / GeV);
    hit->set_pz(1, postPoint->GetMomentum().z() / GeV);

    hit->set_t(1, postPoint->GetGlobalTime() / nanosecond);

    //sum up the energy to get total deposited
    hit->set_edep(hit->get_edep() + edep);
    if (whichactive > 0)
    {
      hit->set_eion(hit->get_eion() + eion);
    }
    if (geantino)
    {
      hit->set_edep(-1);  // only energy=0 g4hits get dropped, this way geantinos survive the g4hit compression
      if (whichactive > 0)
      {
        hit->set_eion(-1);
      }
    }
    if (edep > 0)
    {
      if (G4VUserTrackInformation* p = aTrack->GetUserInformation())
      {
        if (PHG4TrackUserInfoV1* pp = dynamic_cast<PHG4TrackUserInfoV1*>(p))
        {
          pp->SetKeep(1);  // we want to keep the track
        }
      }
    }

    // if any of these conditions is true this is the last step in
    // this volume and we need to save the hit
    // postPoint->GetStepStatus() == fGeomBoundary: track leaves this volume
    // postPoint->GetStepStatus() == fWorldBoundary: track leaves this world
    // (happens when your detector goes outside world volume)
    // postPoint->GetStepStatus() == fAtRestDoItProc: track stops (typically
    // aTrack->GetTrackStatus() == fStopAndKill is also set)
    // aTrack->GetTrackStatus() == fStopAndKill: track ends
    if ((use_g4_steps > 0 && whichactive > 0) ||
        postPoint->GetStepStatus() == fGeomBoundary ||
        postPoint->GetStepStatus() == fWorldBoundary ||
        postPoint->GetStepStatus() == fAtRestDoItProc ||
        aTrack->GetTrackStatus() == fStopAndKill)
    {
      // save only hits with energy deposit (or -1 for geantino)
      if (hit->get_edep())
      {
        savehitcontainer->AddHit(layer_id, hit);
        if (saveshower)
        {
          saveshower->add_g4hit_id(savehitcontainer->GetID(), hit->get_hit_id());
        }

        double rin = sqrt(hit->get_x(0) * hit->get_x(0) + hit->get_y(0) * hit->get_y(0));
        double rout = sqrt(hit->get_x(1) * hit->get_x(1) + hit->get_y(1) * hit->get_y(1));
        if (Verbosity() > 10)
          if ((rin > 69.0 && rin < 70.125) || (rout > 69.0 && rout < 70.125))
          {
            cout << "Added Tpc g4hit with rin, rout = " << rin << "  " << rout
                 << " g4hitid " << hit->get_hit_id() << endl;
            cout << " xin " << hit->get_x(0)
                 << " yin " << hit->get_y(0)
                 << " zin " << hit->get_z(0)
                 << " rin " << rin
                 << endl;
            cout << " xout " << hit->get_x(1)
                 << " yout " << hit->get_y(1)
                 << " zout " << hit->get_z(1)
                 << " rout " << rout
                 << endl;
            cout << " xav " << (hit->get_x(1) + hit->get_x(0)) / 2.0
                 << " yav " << (hit->get_y(1) + hit->get_y(0)) / 2.0
                 << " zav " << (hit->get_z(1) + hit->get_z(0)) / 2.0
                 << " rav " << (rout + rin) / 2.0
                 << endl;
          }

        // ownership has been transferred to container, set to null
        // so we will create a new hit for the next track
        hit = nullptr;
      }
      else
      {
        // if this hit has no energy deposit, just reset it for reuse
        // this means we have to delete it in the dtor. If this was
        // the last hit we processed the memory is still allocated
        hit->Reset();
      }
    }
    // return true to indicate the hit was used
    return true;
  }
  else
  {
    return false;
  }
}

//____________________________________________________________________________..
void PHG4TpcSteppingAction::SetInterfacePointers(PHCompositeNode* topNode)
{
  string hitnodename;
  string absorbernodename;
  if (detector_->SuperDetector() != "NONE")
  {
    hitnodename = "G4HIT_" + detector_->SuperDetector();
    absorbernodename = "G4HIT_ABSORBER_" + detector_->SuperDetector();
  }
  else
  {
    hitnodename = "G4HIT_" + detector_->GetName();
    absorbernodename = "G4HIT_ABSORBER_" + detector_->GetName();
  }

  //now look for the map and grab a pointer to it.
  hits_ = findNode::getClass<PHG4HitContainer>(topNode, hitnodename.c_str());
  absorberhits_ = findNode::getClass<PHG4HitContainer>(topNode, absorbernodename.c_str());

  // if we do not find the node it's messed up.
  if (!hits_)
  {
    std::cout << "PHG4TpcSteppingAction::SetTopNode - unable to find " << hitnodename << std::endl;
  }
  if (!absorberhits_)
  {
    if (Verbosity() > 1)
    {
      cout << "PHG4HcalSteppingAction::SetTopNode - unable to find " << absorbernodename << endl;
    }
  }
}
