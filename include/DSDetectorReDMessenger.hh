#ifndef DSDetectorReDMessenger_h
#define DSDetectorReDMessenger_h 1

#include "G4UImessenger.hh"

class DSDetectorReD;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithABool;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAnInteger;

class DSDetectorReDMessenger : public G4UImessenger {

 public:
  DSDetectorReDMessenger(DSDetectorReD*);
  ~DSDetectorReDMessenger();

  void SetNewValue(G4UIcommand*, G4String);

 private:
  DSDetectorReD* fDetector;
  G4UIdirectory* deteDirectory;
  G4UIdirectory* globalDirectory;
  G4UIcmdWithADoubleAndUnit* fbeamTPCDistanceCmd;
  G4UIcmdWithADoubleAndUnit* fangleThetaXCmd;
  G4UIcmdWithADoubleAndUnit* fanglePhiArCmd;
  // G4UIcmdWithAnInteger*        numberOfNeutronDetectorsCmd;
  // G4UIcmdWithADoubleAndUnit*   neutronDistanceCmd;
  // G4UIcmdWithADoubleAndUnit*   neutronEclipticCmd;
  // G4UIcmdWithADoubleAndUnit*   neutronPhiOffsetCmd;
};

#endif
