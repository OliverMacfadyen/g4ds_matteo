#ifndef DSStackingReD_h
#define DSStackingReD_h 1

#include "DSStackingAction.hh"
#include "DSVStackingAction.hh"
#include "G4ClassificationOfNewTrack.hh"
#include "G4UserStackingAction.hh"

class DSEventAction;
class G4StackManager;
class G4Track;

class DSStackingReD : public DSVStackingAction {
 public:
  DSStackingReD();  // EventAction*);
  virtual ~DSStackingReD();

 public:  // with description
  //---------------------------------------------------------------
  // vitual methods to be implemented by user
  //---------------------------------------------------------------

  virtual G4ClassificationOfNewTrack DSClassifyNewTrack(const G4Track* aTrack);
  //    The parent_ID of the track indicates the origin of it.
  //
  //    G4int parent_ID = aTrack->get_parentID();
  //
  //      parent_ID = 0 : primary particle
  //                > 0 : secondary particle
  //                < 0 : postponed from the previous event

  virtual void DSNewStage();
  //    This method is called by G4StackManager when the urgentStack
  //  becomes empty and contents in the waitingStack are transfered
  //  to the urgentStack.
  //    Note that this method is not called at the beginning of each
  //  event, but "PrepareNewEvent" is called.
  //
  //    In case re-classification of the stacked tracks is needed,
  //  use the following method to request to G4StackManager.
  //    stackManager->ReClassify();
  //  All of the stacked tracks in the waitingStack will be re-classified
  //  by "ClassifyNewTrack" method.
  //    To abort current event, use the following method.
  //    stackManager->clear();
  //
  //  Note that this way is valid and safe only for the case it is called
  //  from this user class. The more global way of event abortion is
  //
  //    G4UImanager * UImanager = G4UImanager::GetUIpointer();
  //    UImanager->ApplyCommand("/event/abort");

  virtual void DSPrepareNewEvent();
  //    This method is called by G4StackManager at the begining of
  //  each event.
  //    Be careful that the urgentStack and the waitingStack of
  //  G4StackManager are empty at this moment, because this method
  //  is called before accepting primary particles. Also, note that
  //  the postponeStack of G4StackManager may have some postponed
  //  tracks.

 private:
  //  DSEventAction* fEventAction;
};

#endif
