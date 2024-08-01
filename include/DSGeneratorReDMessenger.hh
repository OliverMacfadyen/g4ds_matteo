#ifndef DSGeneratorReDMessenger_h
#define DSGeneratorReDMessenger_h 1

#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"
#include "G4UImessenger.hh"

class DSGeneratorReD;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithABool;
class G4UIcmdWithADoubleAndUnit;

class DSGeneratorReDMessenger : public G4UImessenger {

 public:
  DSGeneratorReDMessenger(DSGeneratorReD*);
  ~DSGeneratorReDMessenger();

  void SetNewValue(G4UIcommand*, G4String);

 private:
  DSGeneratorReD* fGenerator;
  G4UIdirectory* beamDirectory;
  G4UIdirectory* neutronDirectory;
  G4UIcmdWithADoubleAndUnit* beamEneCmd;
  G4UIcmdWithADoubleAndUnit* coneOpeningCmd;
  G4UIcmdWithAString* typeCmd;
  G4UIcmdWithADoubleAndUnit* TargetTPCDistanceCmd;
  G4UIcmdWithADoubleAndUnit* angleThetaXCmd;
  G4UIcmdWithADoubleAndUnit* anglePhiXCmd;
};

#endif
