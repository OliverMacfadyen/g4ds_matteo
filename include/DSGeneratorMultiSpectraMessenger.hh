//---------------------------------------------------------------------------//

#ifndef DSGENERATORMULTISPECRAMESSENGER_HH
#define DSGENERATORMULTISPECRAMESSENGER_HH

//---------------------------------------------------------------------------//

#include "G4UIcmdWithAnInteger.hh"
#include "G4UIdirectory.hh"
#include "G4UImessenger.hh"
//#include "G4UIcmdWithAString.hh"
//#include "G4UIcmdWithADoubleAndUnit.hh"
#include "DSGeneratorMultiSpectra.hh"

//---------------------------------------------------------------------------//

class DSGeneratorMultiSpectra;

class DSGeneratorMultiSpectraMessenger : public G4UImessenger {
 public:
  DSGeneratorMultiSpectraMessenger(DSGeneratorMultiSpectra*);
  ~DSGeneratorMultiSpectraMessenger();

  void SetNewValue(G4UIcommand*, G4String);

 private:
  DSGeneratorMultiSpectra* fGenerator;
  G4UIdirectory* fDirectory;

  G4UIcmdWithAnInteger* fNSpectraCmd;

  //    G4UIcmdWithAString*                  fIsotopeCmd;
  //    G4UIcmdWithADoubleAndUnit*           fSetEminCmd;
  //    G4UIcmdWithADoubleAndUnit*           fSetEmaxCmd;
};

#endif
