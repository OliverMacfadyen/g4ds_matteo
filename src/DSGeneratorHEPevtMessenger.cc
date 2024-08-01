//---------------------------------------------------------------------------//
//    Code by pablo.garcia@ciemat.es to read files in HEPEVT format:
//---------------------------------------------------------------------------//
#include "DSGeneratorHEPevtMessenger.hh"
#include "DSGeneratorHEPevt.hh"
#include "DSLogger.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include <fstream>
#include <iomanip>
#include "G4ParticleTable.hh"
#include "G4ThreeVector.hh"
#include "G4Tokenizer.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"
#include "G4ios.hh"

DSGeneratorHEPevtMessenger::DSGeneratorHEPevtMessenger(DSGeneratorHEPevt* fPtclGun) {
  generator = fPtclGun;
  fDirectory = new G4UIdirectory("/ds/generator/hepevt/");
  fDirectory->SetGuidance("Control of HEPevt event generator");

  // name of the file containing the angular spectrum
  fileNameCmd = new G4UIcmdWithAString("/ds/generator/hepevt/filename", this);
  fileNameCmd->SetGuidance("Name of the file containing the events in HEPEVT format");
  fileNameCmd->SetGuidance("Format: text files only so far");
  fileNameCmd->SetGuidance("Default: events.hepevt");
}

DSGeneratorHEPevtMessenger::~DSGeneratorHEPevtMessenger() {
  delete fDirectory;
  delete fileNameCmd;
}

void DSGeneratorHEPevtMessenger::SetNewValue(G4UIcommand* command, G4String newValues) {
  if (command == fileNameCmd) {
    generator->SetFileName(newValues);
    DSLog(trace) << "File of events in HEPEVT format is " << generator->GetFileName() << endlog;
    DSLog(trace) << "Be sure that the format is correct (text file) !" << endlog;
  }
}
