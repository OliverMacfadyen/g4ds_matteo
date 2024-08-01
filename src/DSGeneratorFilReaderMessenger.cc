#include "DSGeneratorFilReaderMessenger.hh"
#include <stdio.h>
#include "DSLogger.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADouble.hh"

using namespace std;

DSGeneratorFilReaderMessenger::DSGeneratorFilReaderMessenger(DSGeneratorFilReader* gen) {
  generator = gen;
  fDirectory = new G4UIdirectory("/ds/generator/filreader/");
  fDirectory->SetGuidance("Control of DSFilReader generator");

  fFileCmd = new G4UIcmdWithAString("/ds/generator/filreader/file", this);
  fSkipEvents = new G4UIcmdWithADouble("/ds/generator/filreader/skip_events", this);
}

DSGeneratorFilReaderMessenger::~DSGeneratorFilReaderMessenger() {

  delete fDirectory;
  delete fFileCmd;
  delete fSkipEvents;
}

void DSGeneratorFilReaderMessenger::SetNewValue(G4UIcommand* cmd, G4String newValue) {
  if (cmd == fFileCmd) {
    DSIO::Get()->SetG4DSFile(newValue);
    DSIO::Get()->SetIsG4DS(true);
  } else if( cmd == fSkipEvents) {
    double n = fSkipEvents->ConvertToDouble(newValue);
    generator->SetNumberOfSkippedEvents(int (n));
    DSLog(routine) << "Number of skipped events " << n << endlog; 
  }
}

/*
 * Revision 1.1  2013/03/20 09:54:28  dfranco
 * New version of g4ds
 *
 *
 */
