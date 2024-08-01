#include "DSGeneratorSCSMessenger.hh"
#include "DSLogger.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

DSGeneratorSCSMessenger::DSGeneratorSCSMessenger(DSGeneratorSCS* gen) {

  fGenerator = gen;
  fDirectory = new G4UIdirectory("/ds/generator/scs/");
  fDirectory->SetGuidance("Control of Special Cross Sections generator");

  fIsotopeCmd = new G4UIcmdWithAString("/ds/generator/scs/isotope", this);
  G4String ioncandidates = "Ar39";
  fIsotopeCmd->SetCandidates(ioncandidates);

  fSetEminCmd = new G4UIcmdWithADoubleAndUnit("/ds/generator/scs/emin", this);
  fSetEminCmd->SetDefaultUnit("MeV");
  fSetEminCmd->SetUnitCandidates("eV keV MeV GeV");

  fSetEmaxCmd = new G4UIcmdWithADoubleAndUnit("/ds/generator/scs/emax", this);
  fSetEmaxCmd->SetDefaultUnit("MeV");
  fSetEmaxCmd->SetUnitCandidates("eV keV MeV GeV");
}

DSGeneratorSCSMessenger::~DSGeneratorSCSMessenger() {
  delete fDirectory;
  delete fIsotopeCmd;
  delete fSetEmaxCmd;
  delete fSetEminCmd;
}

void DSGeneratorSCSMessenger::SetNewValue(G4UIcommand* cmd, G4String newValue) {

  if (cmd == fIsotopeCmd) {
    fGenerator->SetIsotope(newValue);
  } else if (cmd == fSetEminCmd) {
    DSLog(routine) << "SCS Emin = " << newValue << endlog;
    fGenerator->SetEmin(fSetEminCmd->ConvertToDimensionedDouble(newValue));
  } else if (cmd == fSetEmaxCmd) {
    DSLog(routine) << "SCS Emax = " << newValue << endlog;
    fGenerator->SetEmax(fSetEmaxCmd->ConvertToDimensionedDouble(newValue));
  }
}
