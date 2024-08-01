#include "DSGeneratorReDMessenger.hh"

#include "DSGeneratorReD.hh"

#include "DSLogger.hh"
#include "G4IonTable.hh"
#include "G4RunManager.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "globals.hh"

DSGeneratorReDMessenger::DSGeneratorReDMessenger(DSGeneratorReD* gen) {  // thePrim

  fGenerator = gen;
  //: fGeneratorPrimary(gen) //fPrimary(thePrim)

  beamDirectory = new G4UIdirectory("/ds/generator/red/PrimaryBeam/");
  beamDirectory->SetGuidance("ReD Primary Beam control");

  typeCmd = new G4UIcmdWithAString("/ds/generator/red/PrimaryBeam/source", this);
  typeCmd->SetGuidance("Select primary reaction: available: ");
  typeCmd->SetGuidance("Li7+p (Li7), D+D (DD), geantino");
  typeCmd->SetCandidates("Li7 DD geantino");

  beamEneCmd = new G4UIcmdWithADoubleAndUnit("/ds/generator/red/PrimaryBeam/energy", this);
  beamEneCmd->SetGuidance("Energy of the primary beam");
  /*
  TargetTPCDistanceCmd = new
  G4UIcmdWithADoubleAndUnit("/ds/generator/red/PrimaryBeam/position",this);
  TargetTPCDistanceCmd->SetGuidance("Distance between target and TPC's
  center."); TargetTPCDistanceCmd->SetUnitCandidates("cm");
  TargetTPCDistanceCmd->AvailableForStates(G4State_PreInit);
*/
  neutronDirectory = new G4UIdirectory("/ds/generator/red/NeutronBeam/");
  neutronDirectory->SetGuidance("ReD Neutron Beam control");

  coneOpeningCmd = new G4UIcmdWithADoubleAndUnit("/ds/generator/red/NeutronBeam/opening", this);
  coneOpeningCmd->SetGuidance("Opening angle of the neutrons's beam cone");
  coneOpeningCmd->SetGuidance("(default 2 deg) ");
  coneOpeningCmd->SetUnitCandidates("deg rad");
  coneOpeningCmd->AvailableForStates(G4State_Init, G4State_PreInit);

  angleThetaXCmd = new G4UIcmdWithADoubleAndUnit("/ds/generator/red/NeutronBeam/ThetaX", this);
  angleThetaXCmd->SetGuidance(
      "Angle between the HORIZONTAL plane (in the laboratory) containing the "
      "Li beam and the TPC's vector position "
      "(from target to center TPC).");
  // angleThetaXCmd->SetGuidance("Positive angles mean that the TPC is *below*
  // the Li beam");
  angleThetaXCmd->SetUnitCandidates("deg rad");
  angleThetaXCmd->AvailableForStates(G4State_Init, G4State_PreInit);

  anglePhiXCmd = new G4UIcmdWithADoubleAndUnit("/ds/generator/red/NeutronBeam/PhiX", this);
  anglePhiXCmd->SetGuidance(
      "Angle between the VERTICAL plane (in the laboratory) containing the Li "
      "beam and the TPC's vector position (from "
      "target to center TPC).");
  // fanglePhiXCmd->SetGuidance("Positive angles mean that the TPC is *on the
  // left* to the beam looking *forward* with respect to the target");
  anglePhiXCmd->SetUnitCandidates("deg rad");
  anglePhiXCmd->AvailableForStates(G4State_Init, G4State_PreInit);
}

DSGeneratorReDMessenger::~DSGeneratorReDMessenger() {
  delete beamDirectory;
  delete typeCmd;
  delete beamEneCmd;
  delete neutronDirectory;
  delete coneOpeningCmd;
  delete TargetTPCDistanceCmd;
  delete angleThetaXCmd;
  delete anglePhiXCmd;
}

void DSGeneratorReDMessenger::SetNewValue(G4UIcommand* command, G4String newValue) {
  if (command == typeCmd) fGenerator->SetReactionType(newValue);
  else if (command == beamEneCmd)
    fGenerator->SetBeamEnergy(beamEneCmd->GetNewDoubleValue(newValue));
  else if (command == coneOpeningCmd)
    fGenerator->SetConeOpening(coneOpeningCmd->GetNewDoubleValue(newValue));
  // else if(command == beamTPCDistanceCmd){
  // fGenerator->SetBeamToTPCDistance(beamTPCDistanceCmd->GetNewDoubleValue(newValue));
  // DSLog(routine) << "Target to TPC Distance: " << newValue << endlog ;
  // }
  else if (command == angleThetaXCmd) {
    fGenerator->SetNeutronBeamThetaX(angleThetaXCmd->GetNewDoubleValue(newValue));
    DSLog(routine) << "Neutron Beam to TPC Angle Theta: " << newValue << endlog;
  } else if (command == anglePhiXCmd) {
    fGenerator->SetNeutronBeamPhiX(anglePhiXCmd->GetNewDoubleValue(newValue));
    DSLog(routine) << "Neutron Beam to TPC Angle Phi: " << newValue << endlog;
  }
}
