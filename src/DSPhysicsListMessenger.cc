#include "DSPhysicsListMessenger.hh"
#include <G4UnitsTable.hh>
#include "DSIO.hh"
#include "DSLogger.hh"
#include "DSPhysicsList.hh"
#include "DSStorage.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

class DSPhysicsListMessenger;

using namespace std;

DSPhysicsListMessenger::DSPhysicsListMessenger(DSPhysicsList* phys) {
  fPhysicsList = phys;
  fDirectory = new G4UIdirectory("/ds/physics/");
  fDirectory->SetGuidance("Control of physical processes");

  fHadCmd = new G4UIcmdWithAString("/ds/physics/hadronic_list", this);
  G4String hadcandidates = "none HP QGSP_BERT_HP QGSP_BIC_HP FTF_BIC FTFP_BERT_HP Shielding FPFP_BERT";
  fHadCmd->SetCandidates(hadcandidates);

  fEMCmd = new G4UIcmdWithAString("/ds/physics/em_list", this);
  G4String emcandidates = "standard livermore";
  fEMCmd->SetCandidates(emcandidates);

  fHPRangeCutsCmd = new G4UIcmdWithABool("/ds/physics/HPRangeCuts", this);

  fDepositCutsCmd = new G4UIcmdWithABool("/ds/physics/DepositCuts", this);

  fKillS1S2Cmd = new G4UIcmdWithABool("/ds/physics/killS1S2", this);

  fKillS2Cmd = new G4UIcmdWithABool("/ds/physics/killS2", this);

  fKillS1Cmd = new G4UIcmdWithABool("/ds/physics/killS1", this);

  fScaleS2Cmd = new G4UIcmdWithADouble("/ds/physics/scaleS2", this);

  fTunedS1At200VCmd = new G4UIcmdWithABool("/ds/physics/tuned200V", this);

  fCherenkovCmd = new G4UIcmdWithABool("/ds/physics/cherenkov", this);
 
  fMatScintillationCmd = new G4UIcmdWithABool("/ds/physics/material_scintillation", this);

  fDriftFieldCmd = new G4UIcmdWithADoubleAndUnit("/ds/physics/DriftField", this);
  new G4UnitDefinition("V/cm", "V/cm", "Electric field", 100 * volt / m);
  new G4UnitDefinition("kV/cm", "kV/cm", "Electric field", 1e5 * volt / m);
  fDriftFieldCmd->SetUnitCategory("Electric field");
  fDriftFieldCmd->SetDefaultUnit("V/cm");
  fDriftFieldCmd->SetUnitCandidates("V/cm kV/cm");

  fExtractionFieldCmd = new G4UIcmdWithADoubleAndUnit("/ds/physics/ExtractionField", this);
  fExtractionFieldCmd->SetUnitCategory("Electric field");
  fExtractionFieldCmd->SetDefaultUnit("kV/cm");
  fExtractionFieldCmd->SetUnitCandidates("V/cm kV/cm");
}

DSPhysicsListMessenger::~DSPhysicsListMessenger() {
  delete fDirectory;
  delete fHadCmd;
  delete fEMCmd;
  delete fHPRangeCutsCmd;
  delete fKillS2Cmd;
  delete fKillS1S2Cmd;
  delete fScaleS2Cmd;
  delete fDriftFieldCmd;
  delete fExtractionFieldCmd;
  delete fTunedS1At200VCmd;
  delete fMatScintillationCmd;
}

void DSPhysicsListMessenger::SetNewValue(G4UIcommand* command, G4String newValue) {

  if (command == fHadCmd) {
    if (newValue == "none") fPhysicsList->SetHadronicList(-1);
    if (newValue == "HP") fPhysicsList->SetHadronicList(0);
    if (newValue == "QGSP_BERT_HP") fPhysicsList->SetHadronicList(1);
    if (newValue == "QGSP_BIC_HP") fPhysicsList->SetHadronicList(2);
    if (newValue == "FTF_BIC") fPhysicsList->SetHadronicList(3);
    if (newValue == "FTFP_BERT_HP") fPhysicsList->SetHadronicList(4);
    if (newValue == "Shielding") fPhysicsList->SetHadronicList(5);
    if (newValue == "FPFP_BERT") fPhysicsList->SetHadronicList(6);
    DSLog(routine) << " Hadronic List: " << newValue << endlog;
  } else if (command == fEMCmd) {
    if (newValue == "standard") fPhysicsList->SetEMList(1);
    if (newValue == "Livermore") fPhysicsList->SetEMList(2);
    DSLog(routine) << " EM List: " << newValue << endlog;
  } else if (command == fHPRangeCutsCmd) {
    fPhysicsList->SetHPRangeCuts(fHPRangeCutsCmd->ConvertToBool(newValue));
    DSLog(routine) << " High Precision Range Cuts in LAr: " << newValue << endlog;
  } else if (command == fDepositCutsCmd) {
    fPhysicsList->SetDepositCuts(fDepositCutsCmd->ConvertToBool(newValue));
    DSLog(routine) << " Deposit Physics Cuts:" << newValue << endlog;
  } else if (command == fKillS1S2Cmd) {
    DSStorage::Get()->SetKillS1S2(fKillS1S2Cmd->ConvertToBool(newValue));
    DSLog(routine) << " Kill S1 and S2 Light: " << newValue << endlog;
  } else if (command == fKillS2Cmd) {
    DSStorage::Get()->SetKillS2(fKillS2Cmd->ConvertToBool(newValue));
    DSLog(routine) << " Kill S2 Light: " << newValue << endlog;
  } else if (command == fKillS1Cmd) {
    DSStorage::Get()->SetKillS1(fKillS1Cmd->ConvertToBool(newValue));
    DSLog(routine) << " Kill S1 Light: " << newValue << endlog;
  } else if (command == fCherenkovCmd) {
    DSStorage::Get()->SetIsCherenkov(fCherenkovCmd->ConvertToBool(newValue));
    DSLog(routine) << " Cherenkov Process: " << newValue << endlog;
  } else if (command == fScaleS2Cmd) {
    DSStorage::Get()->SetScaleS2(fScaleS2Cmd->ConvertToDouble(newValue));
    DSLog(routine) << " Scale S2 Light by: " << newValue << endlog;
  } else if (command == fDriftFieldCmd) {
    DSStorage::Get()->SetDriftField(fDriftFieldCmd->ConvertToDimensionedDouble(newValue));
    DSStorage::Get()->SetTunedS1At200V(false);
    DSLog(routine) << " Drift Field: " << newValue << endlog;
  } else if (command == fExtractionFieldCmd) {
    DSStorage::Get()->SetExtractionField(fExtractionFieldCmd->ConvertToDimensionedDouble(newValue));
    DSLog(routine) << " Extraction Field: " << newValue << endlog;
  } else if (command == fTunedS1At200VCmd) {
    DSStorage::Get()->SetTunedS1At200V(fTunedS1At200VCmd->ConvertToBool(newValue));
    DSLog(routine) << " Tuned S1 and S2 at 200 V/cm: " << newValue << endlog;
  } else if (command == fMatScintillationCmd) {
    DSStorage::Get()->SetMatScintillation(command->ConvertToBool(newValue));
    DSLog(routine) << " ESR/TPB scintillation: " << newValue << endlog;
  }
}

