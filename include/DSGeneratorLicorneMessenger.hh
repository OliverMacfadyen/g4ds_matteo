//---------------------------------------------------------------------------//

#ifndef DSGeneratorLicorneMessenger_HH
#define DSGeneratorLicorneMessenger_HH

//---------------------------------------------------------------------------//

#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"
#include "G4UImessenger.hh"

#include "DSGeneratorLicorne.hh"
//---------------------------------------------------------------------------//

class DSGeneratorLicorne;

class DSGeneratorLicorneMessenger : public G4UImessenger {

 public:
  DSGeneratorLicorneMessenger(DSGeneratorLicorne*);
  ~DSGeneratorLicorneMessenger();

  void SetNewValue(G4UIcommand*, G4String);

 private:
  DSGeneratorLicorne* fGenerator;
  G4UIdirectory* fDirectory;
  // G4UIcmdWithAString* fPosition;
  G4UIcmdWith3VectorAndUnit* fPositionCmd;
  G4UIcmdWith3Vector* fDirectionCmd;
  G4UIcmdWithADoubleAndUnit* fRunTimeCmd;
  G4UIcmdWithADoubleAndUnit* fNeutronRateCmd;
  G4UIcmdWithADouble* fGammaToNeutronRatioCmd;
  G4UIcmdWithADoubleAndUnit* fPulsePeriodCmd;
  G4UIcmdWithADoubleAndUnit* fPulseWidthCmd;
  G4UIcmdWithABool* fPulseModeCmd;
};

#endif
