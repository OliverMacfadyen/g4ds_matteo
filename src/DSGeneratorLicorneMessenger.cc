#include "DSGeneratorLicorneMessenger.hh"
#include "DSLogger.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

DSGeneratorLicorneMessenger::DSGeneratorLicorneMessenger(DSGeneratorLicorne* gen) {

  fGenerator = gen;
  fDirectory = new G4UIdirectory("/ds/generator/licorne/");
  fDirectory->SetGuidance("Control of Licorne Neutron generator");

  fPositionCmd = new G4UIcmdWith3VectorAndUnit("/ds/generator/licorne/position", this);
  fPositionCmd->SetGuidance("Set the licorne beam position");
  fPositionCmd->SetUnitCategory("Length");
  fPositionCmd->SetDefaultUnit("cm");
  fPositionCmd->SetUnitCandidates("mm cm m");

  fDirectionCmd = new G4UIcmdWith3Vector("/ds/generator/licorne/direction", this);
  fDirectionCmd->SetGuidance("Set the gun direction");

  fPulseModeCmd = new G4UIcmdWithABool("/ds/generator/licorne/pulse_mode", this);

  fRunTimeCmd = new G4UIcmdWithADoubleAndUnit("/ds/generator/licorne/run_time", this);
  fRunTimeCmd->SetUnitCandidates("ns s microsecond ms");

  fNeutronRateCmd = new G4UIcmdWithADoubleAndUnit("/ds/generator/licorne/neutron_rate", this);
  fNeutronRateCmd->SetUnitCandidates("hertz kilohertz");

  fGammaToNeutronRatioCmd = new G4UIcmdWithADouble("/ds/generator/licorne/gamma_neutron_ratio", this);

  fPulsePeriodCmd = new G4UIcmdWithADoubleAndUnit("/ds/generator/licorne/pulse_period", this);
  fPulsePeriodCmd->SetUnitCandidates("ns microsecond ms s");

  fPulseWidthCmd = new G4UIcmdWithADoubleAndUnit("/ds/generator/licorne/pulse_width", this);
  fPulseWidthCmd->SetUnitCandidates("ns microsecond ms s");
}

DSGeneratorLicorneMessenger::~DSGeneratorLicorneMessenger() {
  delete fDirectory;
  delete fRunTimeCmd;
  delete fNeutronRateCmd;
  delete fGammaToNeutronRatioCmd;
  delete fPulseWidthCmd;
  delete fPulsePeriodCmd;
  delete fPulseModeCmd;
}

void DSGeneratorLicorneMessenger::SetNewValue(G4UIcommand* cmd, G4String newValue) {

  if (cmd == fPositionCmd) {
    fGenerator->SetVParticlePosition(fPositionCmd->ConvertToDimensioned3Vector(newValue));
  } else if (cmd == fDirectionCmd) {
    fGenerator->SetVParticleDirection(fDirectionCmd->ConvertTo3Vector(newValue));
  } else if (cmd == fPulseModeCmd) {
    DSLog(routine) << "Licorne Pulse Mode = " << newValue << endlog;
    fGenerator->SetPulseMode(fPulseModeCmd->ConvertToBool(newValue));
  } else if (cmd == fRunTimeCmd) {
    DSLog(routine) << "Acquisition Run Time = " << newValue << endlog;
    fGenerator->SetRunTime(fRunTimeCmd->ConvertToDimensionedDouble(newValue));
  } else if (cmd == fNeutronRateCmd) {
    DSLog(routine) << "Neutron Rate = " << newValue << endlog;
    fGenerator->SetNeutronRate(fNeutronRateCmd->ConvertToDimensionedDouble(newValue));
  } else if (cmd == fGammaToNeutronRatioCmd) {
    DSLog(routine) << "Gamma To Neutron Ratio = " << newValue << endlog;
    fGenerator->SetGammaToNeutronRatio(fGammaToNeutronRatioCmd->ConvertToDouble(newValue));
  } else if (cmd == fPulseWidthCmd) {
    DSLog(routine) << "Neutron Pulse Width = " << newValue << endlog;
    fGenerator->SetPulseWidth(fPulseWidthCmd->ConvertToDimensionedDouble(newValue));
  } else if (cmd == fPulsePeriodCmd) {
    DSLog(routine) << "Neutron Pulse Period = " << newValue << endlog;
    fGenerator->SetPulsePeriod(fPulsePeriodCmd->ConvertToDimensionedDouble(newValue));
  }
}
