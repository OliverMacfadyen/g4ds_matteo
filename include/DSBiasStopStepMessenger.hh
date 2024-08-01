#ifndef DSBiasStopStepMessenger_h
#define DSBiasStopStepMessenger_h 1

#include "DSIO.hh"
#include "G4ParticleGun.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UImessenger.hh"

class DSBiasStopStep;

using namespace std;

class DSBiasStopStepMessenger : public G4UImessenger {

 public:
  DSBiasStopStepMessenger(DSBiasStopStep*);
  ~DSBiasStopStepMessenger();

  void SetNewValue(G4UIcommand*, G4String);

 private:
  DSBiasStopStep* stopper;
  G4UIdirectory* fDirectory;
  G4UIcmdWithAString* fBiasTypeCmd;
  G4UIcmdWithAString* fMaterialNameOutCmd;
  G4UIcmdWithAString* fMaterialNameInCmd;
  G4UIcmdWith3VectorAndUnit *fCenterCmd;
  G4UIcmdWithADoubleAndUnit *fRadiusCmd;
  G4UIcmdWithADoubleAndUnit *fSideCmd;
  G4UIcmdWithADoubleAndUnit *fHalfZCmd;


};

#endif


