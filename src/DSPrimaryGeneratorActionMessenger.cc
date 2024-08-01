#include "DSPrimaryGeneratorActionMessenger.hh"
#include "DSGeneratorAmBeSource.hh"
#include "DSGeneratorAmCSource.hh"
#include "DSGeneratorCosmicRayMuons.hh"
#include "DSGeneratorEnergyDeposit.hh"
#include "DSGeneratorFLUKA.hh"
#include "DSGeneratorG4Gun.hh"
#include "DSGeneratorHEPevt.hh"
#include "DSGeneratorKr85m.hh"
#include "DSGeneratorLicorne.hh"
#include "DSGeneratorMultiEvent.hh"
#include "DSGeneratorBias.hh"
#include "DSGeneratorNeutronsAtGS.hh"
#include "DSGeneratorRDMDecayGun.hh"
#include "DSGeneratorReD.hh"
#include "DSGeneratorSCS.hh"
#include "DSLogger.hh"
#include "DSPrimaryGeneratorAction.hh"
#include "DSStorage.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "DSGeneratorFilReader.hh"

using namespace std;

DSPrimaryGeneratorActionMessenger::DSPrimaryGeneratorActionMessenger(DSPrimaryGeneratorAction* act) {
  fGeneratorPrimary = act;

  fDirectory = new G4UIdirectory("/ds/generator/");
  fDirectory->SetGuidance("Control commands for generators:");
  fDirectory->SetGuidance("/ds/generator/select: Select generator.");

  // /DS/generator/select command
  fSelectCmd = new G4UIcmdWithAString("/ds/generator/select", this);
  fSelectCmd->SetGuidance("Selects generator for events.");
  fSelectCmd->SetGuidance("Options are:");
  fSelectCmd->SetGuidance("G4gun: Standard G4 gun.");
  fSelectCmd->SetGuidance("CosmicRayMuons: cosmic-ray muons generator");
  fSelectCmd->SetGuidance("NeutronsAtGS: neutrons generator");
  fSelectCmd->SetGuidance("AmBeSource: AmBe generator");
  fSelectCmd->SetGuidance("AmCSource: AmC generator");
  fSelectCmd->SetGuidance("ReD: ReD generator");
  fSelectCmd->SetGuidance("HEPevt: read events from fine in HEPEVT format");
  fSelectCmd->SetGuidance("Bias: apply the biasing approach");

  G4String candidates =
      "G4Gun CosmicRayMuons NeutronsAtGS MultiEvent MultiSpectra RDM SCS "
      "EnergyDeposit AmBeSource Kr85m Licorne "
      "AmCSource FLUKA ReD HEPevt Bias FilReader";
  fSelectCmd->SetCandidates(candidates);

  fFLUKAfilenameCmd = new G4UIcmdWithAString("/ds/generator/fluka_filename", this);
}

DSPrimaryGeneratorActionMessenger::~DSPrimaryGeneratorActionMessenger() {

  delete fDirectory;
  delete fGeneratorPrimary;
  delete fFLUKAfilenameCmd;
}

void DSPrimaryGeneratorActionMessenger::SetNewValue(G4UIcommand* command, G4String newValue) {
  if (command == fSelectCmd) {

    DSLog(routine) << "Generator: " << newValue << endlog;

    if (newValue == "G4Gun") {
      fGeneratorPrimary->SetDSGenerator(new DSGeneratorG4Gun);
    } else if (newValue == "CosmicRayMuons") {
      fGeneratorPrimary->SetDSGenerator(new DSGeneratorCosmicRayMuons);
    } else if (newValue == "NeutronsAtGS") {
      fGeneratorPrimary->SetDSGenerator(new DSGeneratorNeutronsAtGS);
    } else if (newValue == "MultiEvent") {
      fGeneratorPrimary->SetDSGenerator(new DSGeneratorMultiEvent);
    } else if (newValue == "EnergyDeposit") {
      fGeneratorPrimary->SetDSGenerator(new DSGeneratorEnergyDeposit);
    } else if (newValue == "RDM") {
      fGeneratorPrimary->SetDSGenerator(new DSGeneratorRDMDecayGun);
    } else if (newValue == "SCS") {
      fGeneratorPrimary->SetDSGenerator(new DSGeneratorSCS);
    } else if (newValue == "AmBeSource") {
      fGeneratorPrimary->SetDSGenerator(new DSGeneratorAmBeSource);
    } else if (newValue == "AmCSource") {
      fGeneratorPrimary->SetDSGenerator(new DSGeneratorAmCSource);
    } else if (newValue == "Kr85m") {
      fGeneratorPrimary->SetDSGenerator(new DSGeneratorKr85m);
    } else if (newValue == "Licorne") {
      fGeneratorPrimary->SetDSGenerator(new DSGeneratorLicorne);
    } else if (newValue == "FLUKA") {
      fGeneratorPrimary->SetDSGenerator(new DSGeneratorFLUKA);
    } else if (newValue == "ReD") {
      fGeneratorPrimary->SetDSGenerator(new DSGeneratorReD);
    } else if (newValue == "HEPevt") {
      fGeneratorPrimary->SetDSGenerator(new DSGeneratorHEPevt);
    } else if (newValue == "Bias") {
      fGeneratorPrimary->SetDSGenerator(new DSGeneratorBias);
    } else if (newValue == "FilReader") {
      fGeneratorPrimary->SetDSGenerator(new DSGeneratorFilReader);
    }
  }

  // set the filename for the FLUKA generator (not worth to build a messenger
  // for just one command...?)
  if (command == fFLUKAfilenameCmd) { DSStorage::Get()->SetFLUKAfilename(newValue); }
}

/*
 * Revision 1.3  2013/03/22 13:50:09  dfranco
 * Added radioactive decay and radioactive decay chain generators. They require
 * the correspondent stacking actions. Two mac files are included as examples
 *
 */
