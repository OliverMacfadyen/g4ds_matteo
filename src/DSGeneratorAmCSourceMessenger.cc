#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"

#include "DSEventHandler.hh"
#include "DSGeneratorAmCSource.hh"
#include "DSGeneratorAmCSourceMessenger.hh"
#include "DSLogger.hh"

using namespace std;

DSGeneratorAmCSourceMessenger::DSGeneratorAmCSourceMessenger(DSGeneratorAmCSource* gen) {

  generator = gen;
  fDirectory = new G4UIdirectory("/ds/generator/AmC/");
  fDirectory->SetGuidance("Control of DSAmCSource event generator");

  fAmCEnergyAngleCorrelationFileCmd = new G4UIcmdWithAString("/ds/generator/AmC/energyanglecorrelationfile", this);
  fAmCEnergyAngleCorrelationFileCmd->SetGuidance(
      "Select the correction between energy and angle distribution for AmC "
      "source");
  fSourceRotationCmd = new G4UIcmdWith3Vector("/ds/generator/AmC/rotation", this);
}

DSGeneratorAmCSourceMessenger::~DSGeneratorAmCSourceMessenger() {

  delete fDirectory;
  delete fSourceRotationCmd;
  delete fAmCEnergyAngleCorrelationFileCmd;
}

void DSGeneratorAmCSourceMessenger::SetNewValue(G4UIcommand* cmd, G4String newValue) {
  if (cmd == fAmCEnergyAngleCorrelationFileCmd) {
    DSLog(routine) << "AmC Energy Angle Correlation File Name  = " << newValue << endlog;
    generator->SetAmCEnergyAngleCorrelationFileName(newValue);
  } else if (cmd == fSourceRotationCmd) {
    generator->SetAmCSourceRotation(fSourceRotationCmd->ConvertTo3Vector(newValue));
  } else {
    DSLog(error) << "invalid value: " << newValue << ", should be one of AngleFile', 'EnergyFile'" << endlog;
  }
}
