
#include "DSStackingActionMessenger.hh"
#include "DSStackingAction.hh"
#include "DSStackingDefault.hh"
#include "DSStackingRDM.hh"
#include "DSStackingRDMChain.hh"
#include "DSStackingReD.hh"
#include "DSStackingNeutronCapture.hh"
#include "DSStackingOnlyOneParticle.hh"
#include "DSStorage.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"

using namespace std;

DSStackingActionMessenger::DSStackingActionMessenger(DSStackingAction* stack) {
  fStacking = stack;
  fDirectory = new G4UIdirectory("/ds/stack/");
  fDirectory->SetGuidance("Control commands for stack:");
  // /DS/generator/select command
  fSelectCmd = new G4UIcmdWithAString("/ds/stack/select", this);
  fSelectCmd->SetGuidance("Selects stack options");
  G4String candidates = "RDM RDMChain ReD Default NeutronCapture OnlyOneParticle";
  fSelectCmd->SetCandidates(candidates);
}

DSStackingActionMessenger::~DSStackingActionMessenger() {

  delete fDirectory;
  delete fSelectCmd;
  if (fStacking) delete fStacking;
}

void DSStackingActionMessenger::SetNewValue(G4UIcommand* command, G4String newValue) {
  if (command == fSelectCmd) {
    if (newValue == "RDM") {
      DSStorage::Get()->SetRDMDecay(true);
      fStacking->SetDSStackingAction(new DSStackingRDM);
    } else if (newValue == "RDMChain") {
      DSStorage::Get()->SetRDMDecay(true);
      DSStorage::Get()->SetIsDecayInVacuum(false);
      fStacking->SetDSStackingAction(new DSStackingRDMChain);
    } else if (newValue == "ReD") {
      DSStorage::Get()->SetReDStacking(true);
      fStacking->SetDSStackingAction(new DSStackingReD);
    } else if (newValue == "Default") {
      fStacking->SetDSStackingAction(new DSStackingDefault);
    } else if (newValue == "NeutronCapture") {
      fStacking->SetDSStackingAction(new DSStackingNeutronCapture);
   } else if (newValue == "OnlyOneParticle") {
      fStacking->SetDSStackingAction(new DSStackingOnlyOneParticle);
    }
  }
}

/*
 * $Log: DSStackingActionMessenger.cc,v $
 * Revision 1.1  2014/05/07 12:21:05  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.6  2013/07/07 09:42:16  dfranco
 * Added a default stacking for killing particles from long living nucleus
 * decays
 *
 *
 */
