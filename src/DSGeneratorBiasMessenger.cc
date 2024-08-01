#include "DSGeneratorBiasMessenger.hh"
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

DSGeneratorBiasMessenger::DSGeneratorBiasMessenger(DSGeneratorBias* gen) {
  generator = gen;
  fDirectory = new G4UIdirectory("/ds/generator/bias/");
  fDirectory->SetGuidance("Control of DSBias generator");

  fFileCmd                = new G4UIcmdWithAString("/ds/generator/bias/file", this);
  fSkipEventsCmd          = new G4UIcmdWithADouble("/ds/generator/bias/skip_events", this);
  fAmplificationFactorCmd = new G4UIcmdWithAnInteger("/ds/generator/bias/amplification", this); 
  fDirectionSmearingCmd   = new G4UIcmdWithADouble("/ds/generator/bias/smear_direction", this); 
  fEnergySmearingCmd      = new G4UIcmdWithADouble("/ds/generator/bias/smear_energy", this); 
  fMinimumEnergyCmd       = new G4UIcmdWithADoubleAndUnit("/ds/generator/bias/minimum_energy", this); 

}

DSGeneratorBiasMessenger::~DSGeneratorBiasMessenger() {

  delete fDirectory;
  delete fFileCmd;
  delete fSkipEventsCmd;
  delete fAmplificationFactorCmd;
  delete fDirectionSmearingCmd;
  delete fEnergySmearingCmd;
  delete fMinimumEnergyCmd;

}

void DSGeneratorBiasMessenger::SetNewValue(G4UIcommand* cmd, G4String newValue) {
  if (cmd == fFileCmd) {
    DSIO::Get()->SetG4DSFile(newValue);
    DSIO::Get()->SetIsG4DS(true);
  } else if( cmd == fSkipEventsCmd) {
    generator->SetNumberOfSkippedEvents(cmd->ConvertToInt(newValue));
    DSLog(routine) << "Number of skipped events " << newValue << endlog; 
  } else if( cmd == fAmplificationFactorCmd) {
    generator->SetAmplificationFactor(cmd->ConvertToInt(newValue));
    DSLog(routine) << "Bias amplification factor: " << cmd->ConvertToInt(newValue) << endlog; 
  } else if( cmd == fDirectionSmearingCmd) {
    generator->SetDirectionSmearing(cmd->ConvertToDouble(newValue));
    DSLog(routine) << "Direction smeared by " << newValue << " degree" << endlog; 
  } else if( cmd == fEnergySmearingCmd) {
    generator->SetEnergySmearing(cmd->ConvertToDouble(newValue)/100.);
    DSLog(routine) << "Energy smeared by " << newValue << " %" <<  endlog; 
  } else if( cmd == fMinimumEnergyCmd) {
    generator->SetMinimumEnergy(cmd->ConvertToDimensionedDouble(newValue));
    DSLog(routine) << "Do not track events with E > " << cmd->ConvertToDimensionedDouble(newValue)/keV << " keV" <<  endlog; 
  } 
}


