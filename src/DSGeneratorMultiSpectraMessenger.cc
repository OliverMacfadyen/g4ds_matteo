#include "DSGeneratorMultiSpectraMessenger.hh"
#include "DSLogger.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

DSGeneratorMultiSpectraMessenger::DSGeneratorMultiSpectraMessenger(DSGeneratorMultiSpectra* gen) {
  fGenerator = gen;

  fDirectory = new G4UIdirectory("/ds/generator/multispectra/");
  fDirectory->SetGuidance("Control of Multi Spectra generator");

  fNSpectraCmd = new G4UIcmdWithAnInteger("/ds/generator/multispectra/nspectra", this);
  fNSpectraCmd->SetDefaultValue(0);

  //   fIsotopeCmd = new
  //   G4UIcmdWithAString("/ds/generator/multispectra/isotope",this); G4String
  //   ioncandidates = "Ar39"; fIsotopeCmd->SetCandidates( ioncandidates );
  //
  //   fSetEminCmd = new
  //   G4UIcmdWithADoubleAndUnit("/ds/generator/multispectra/emin",this);
  //   fSetEminCmd->SetDefaultUnit("MeV");
  //   fSetEminCmd->SetUnitCandidates("eV keV MeV GeV");
  //
  //   fSetEmaxCmd = new
  //   G4UIcmdWithADoubleAndUnit("/ds/generator/multispectra/emax",this);
  //   fSetEmaxCmd->SetDefaultUnit("MeV");
  //   fSetEmaxCmd->SetUnitCandidates("eV keV MeV GeV");
}

DSGeneratorMultiSpectraMessenger::~DSGeneratorMultiSpectraMessenger() {
  if (fDirectory) delete fDirectory;

  if (fNSpectraCmd) delete fNSpectraCmd;

  //   if (fIsotopeCmd)  delete fIsotopeCmd;
  //   if (fSetEmaxCmd)  delete fSetEmaxCmd;
  //   if (fSetEminCmd)  delete fSetEminCmd;
}

void DSGeneratorMultiSpectraMessenger::SetNewValue(G4UIcommand* cmd, G4String newValue) {
  if (cmd == fNSpectraCmd) {
    DSLog(routine) << "Number of Spectra = " << newValue << endlog;
    fGenerator->SetNSpectra(fNSpectraCmd->ConvertToInt(newValue));
  }
  //   else if(cmd == fIsotopeCmd)
  //   {
  //	  fGenerator->SetIsotope( newValue );
  //   }
  //   else if (cmd == fSetEminCmd)
  //   {
  //	  DSLog(routine) << "Spectra's Emin = " << newValue  << endlog ;
  //	  fGenerator->SetEmin(fSetEminCmd->ConvertToDimensionedDouble(newValue));
  //   }
  //   else if (cmd == fSetEmaxCmd)
  //   {
  //	  DSLog(routine) << "Spectra's Emax = " << newValue  << endlog ;
  //	  fGenerator->SetEmax(fSetEmaxCmd->ConvertToDimensionedDouble(newValue));
  //   }
}
