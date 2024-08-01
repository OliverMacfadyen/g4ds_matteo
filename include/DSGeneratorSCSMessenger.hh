//---------------------------------------------------------------------------//

#ifndef DSGeneratorSCSMessenger_HH
#define DSGeneratorSCSMessenger_HH

//---------------------------------------------------------------------------//

#include "DSGeneratorSCS.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"
#include "G4UImessenger.hh"

//---------------------------------------------------------------------------//

class DSGeneratorSCS;

class DSGeneratorSCSMessenger : public G4UImessenger {

 public:
  DSGeneratorSCSMessenger(DSGeneratorSCS*);
  ~DSGeneratorSCSMessenger();

  void SetNewValue(G4UIcommand*, G4String);

 private:
  DSGeneratorSCS* fGenerator;

  G4UIdirectory* fDirectory;
  G4UIcmdWithAString* fIsotopeCmd;

  G4UIcmdWithADoubleAndUnit* fSetEminCmd;
  G4UIcmdWithADoubleAndUnit* fSetEmaxCmd;
};

#endif
