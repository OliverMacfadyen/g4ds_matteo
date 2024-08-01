#ifndef DSTRACKINGACTION_H
#define DSTRACKINGACTION_H

#include "G4ThreeVector.hh"
#include "G4UserTrackingAction.hh"

class DSTrackingAction : public G4UserTrackingAction {

 public:
  DSTrackingAction();
  virtual ~DSTrackingAction(){};
  virtual void PreUserTrackingAction(const G4Track*);
  virtual void PostUserTrackingAction(const G4Track*);

 private:
  int GetDS20kNChannel(G4ThreeVector);
};

#endif
