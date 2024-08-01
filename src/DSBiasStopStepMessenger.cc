#include "DSBiasStopStepMessenger.hh"
#include "DSEventHandler.hh"
#include "DSLogger.hh"
#include "DSBiasStopStep.hh"
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

DSBiasStopStepMessenger::DSBiasStopStepMessenger(DSBiasStopStep* _stopper) {
  stopper = _stopper;

  fDirectory = new G4UIdirectory("/ds/bias/");

  
  fBiasTypeCmd = new G4UIcmdWithAString("/ds/bias/type", this); 
  fBiasTypeCmd->SetCandidates("Material Sphere Cylinder Cube");
  
  fMaterialNameOutCmd = new G4UIcmdWithAString("/ds/bias/material_out", this);
  fMaterialNameInCmd  = new G4UIcmdWithAString("/ds/bias/material_in", this);
  fCenterCmd          = new G4UIcmdWith3VectorAndUnit("/ds/bias/center", this);
  fRadiusCmd          = new G4UIcmdWithADoubleAndUnit("/ds/bias/radius", this);
  fSideCmd            = new G4UIcmdWithADoubleAndUnit("/ds/bias/side", this);
  fHalfZCmd           = new G4UIcmdWithADoubleAndUnit("/ds/bias/half_z", this);
  


}

DSBiasStopStepMessenger::~DSBiasStopStepMessenger() {

  delete fBiasTypeCmd;
  delete fMaterialNameOutCmd;
  delete fMaterialNameInCmd;
  delete fCenterCmd;
  delete fRadiusCmd;
  delete fSideCmd;
  delete fHalfZCmd;
}

void DSBiasStopStepMessenger::SetNewValue(G4UIcommand* cmd, G4String newValue) {
  if (cmd == fBiasTypeCmd) {
    stopper->SetBiasGeometryType(newValue);
    DSLog(routine) << "Stop particles at: " << newValue << endlog;
  } else if (cmd == fMaterialNameOutCmd) {
    stopper->SetMaterialNameOut(newValue);
    DSLog(routine) << "Bias material-out: " << newValue << endlog; 
  } else if (cmd == fMaterialNameInCmd) {
    stopper->SetMaterialNameOut(newValue);
    DSLog(routine) << "Bias material-in: " << newValue << endlog; 
  } else if (cmd == fCenterCmd) {
    stopper->SetCenter(cmd->ConvertToDimensioned3Vector(newValue));
    DSLog(routine) << "Bias volume center: " << newValue << endlog; 
  } else if (cmd == fRadiusCmd) {
    stopper->SetRadius(cmd->ConvertToDimensionedDouble(newValue));
    DSLog(routine) << "Bias volume radius: " << newValue << endlog; 
  } else if (cmd == fSideCmd) {
    stopper->SetSide(cmd->ConvertToDimensionedDouble(newValue));
    DSLog(routine) << "Bias volume side: " << newValue << endlog; 
  } else if (cmd == fSideCmd) {
    stopper->SetHalfZ(cmd->ConvertToDimensionedDouble(newValue));
    DSLog(routine) << "Bias volume half Z: " << newValue << endlog; 
  }
}
