//---------------------------------------------------------------------------//
//    Code by pablo.garcia@ciemat.es to read files in HEPEVT format:
//---------------------------------------------------------------------------//
#ifndef DSGeneratorHEPevtMessenger_h

#define DSGeneratorHEPevtMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class DSGeneratorHEPevt;

class G4ParticleTable;
class G4UIcommand;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithADouble;
class G4UIcmdWithAString;

class DSGeneratorHEPevtMessenger : public G4UImessenger {

 public:
  DSGeneratorHEPevtMessenger(DSGeneratorHEPevt*);
  ~DSGeneratorHEPevtMessenger();

  void SetNewValue(G4UIcommand* command, G4String newValues);

 private:
  DSGeneratorHEPevt* generator;
  G4UIdirectory* fDirectory;
  G4UIcmdWithAString* fileNameCmd;
};

#endif
