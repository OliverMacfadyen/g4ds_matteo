#include "DSStackingOnlyOneParticleMessenger.hh"
#include "DSEventHandler.hh"
#include "DSLogger.hh"
#include "DSStackingOnlyOneParticle.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Tokenizer.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
// using namespace std;

DSStackingOnlyOneParticleMessenger::DSStackingOnlyOneParticleMessenger(DSStackingOnlyOneParticle* stack) {
  stacking = stack;

  fDirectory = new G4UIdirectory("/ds/stack/only_one_particle/");

  
  fParticleNameCmd = new G4UIcmdWithAString("/ds/stack/only_one_particle/name", this);
  

   fParticlePDGCmd = new G4UIcmdWithAnInteger("/ds/stack/only_one_particle/pdg", this);

}

DSStackingOnlyOneParticleMessenger::~DSStackingOnlyOneParticleMessenger() {

  delete fParticleNameCmd;
  delete fParticlePDGCmd;
}

void DSStackingOnlyOneParticleMessenger::SetNewValue(G4UIcommand* cmd, G4String newValue) {
  if (cmd == fParticlePDGCmd) {
    stacking->SetParticlePDG(fParticlePDGCmd->ConvertToInt(newValue));
    DSLog(routine) << "Store only particles with PDG code: " << newValue << endlog;
  } else if (cmd == fParticleNameCmd) {
     stacking->SetParticleName(newValue);
    DSLog(routine) << "Store only " << newValue << endlog;
 }
}
