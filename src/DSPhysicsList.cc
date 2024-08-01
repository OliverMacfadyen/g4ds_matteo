#include "DSLogger.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "DSPhysicsList.hh"
#include "DSPhysicsListMessenger.hh"

#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "globals.hh"

#include "G4NeutronCaptureProcess.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleWithCuts.hh"

#include "G4ProductionCuts.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"

#include <iomanip>
#include "G4ios.hh"

#include "G4UserLimits.hh"
#include "G4Version.hh"

#include "DSOpticalPhysics.hh"
#include "G4DecayPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4HadronElasticPhysicsHP.hh"
#include "G4HadronPhysicsFTF_BIC.hh"
#include "G4HadronPhysicsQGSP_BERT_HP.hh"
#include "G4HadronPhysicsFTFP_BERT_HP.hh"
#include "G4HadronPhysicsFTFP_BERT.hh"
#include "G4HadronPhysicsQGSP_BIC_HP.hh"
#include "G4HadronPhysicsShielding.hh"
#include "G4IonPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4StepLimiterPhysics.hh"


#include "DSEventHandler.hh"
#include "DSStorage.hh"

// Constructor /////////////////////////////////////////////////////////////
DSPhysicsList::DSPhysicsList() : G4VModularPhysicsList() {
  VerboseLevel = 1;
  OpVerbLevel = 0;

  fHadronic = 0; //default: FPFT_BERT_HP
  fEM = 2;

  // DSStorage::Get()->SetIsCherenkov(false);

  SetVerboseLevel(VerboseLevel);
  isHPRangeCuts = true;
  isHadronicHP  = false;
  isDepositCuts = false;

  defaultCutValue = 1.0 * mm;  //

  // Standard EM
  fEMPhysics = nullptr;

  //Decay
  RegisterPhysics(new G4DecayPhysics());
  // Photo-nuclear and lepto-nuclear
  RegisterPhysics(new G4EmExtraPhysics());
  // Radioactive Decay Physics
  RegisterPhysics(new G4RadioactiveDecayPhysics());

  // OpticalPhysics
  fOptPhysics = new DSOpticalPhysics(OpVerbLevel);
  RegisterPhysics(fOptPhysics);
  fHadPhysics = nullptr;
  fMessenger = new DSPhysicsListMessenger(this);

  //Step Limiter (1/3)
  G4StepLimiterPhysics* stepLimitPhys = new G4StepLimiterPhysics();
  RegisterPhysics(stepLimitPhys);

}

// Destructor //////////////////////////////////////////////////////////////
DSPhysicsList::~DSPhysicsList() {
  if (fMessenger) delete fMessenger;
}

// Construct Processes //////////////////////////////////////////////////////
void DSPhysicsList::ConstructProcess() {

  G4VUserPhysicsList::AddTransportation();

  if (fEM == 1) {
    DSLog(routine) << "Standard EM Physics List Active" << endlog;
    fEMPhysics = new G4EmStandardPhysics(VerboseLevel);
  } else if (fEM == 2) {
    DSLog(routine) << "Livermore EM Physics List Active" << endlog;
    fEMPhysics = new G4EmLivermorePhysics(VerboseLevel);
  }
  if (fEMPhysics) fEMPhysics->ConstructProcess();

  fOptPhysics->SetVerboseLevel(OpVerbLevel);
  DSLog(routine) << "Optical Physics Active with verbosity = " << OpVerbLevel << endlog;
  fOptPhysics->ConstructProcess();

  // Warning fHadronic == 0 is not implemented!
  if (fHadronic == 0){ //Default physics list -> equal to fHadronic == 4
    // ConstructHad_Custom();
    DSLog(routine) << "Hadronic Physics Active: FTFP_BERT_HP model" << endlog;
    fHadPhysics = new G4HadronPhysicsFTFP_BERT_HP();
    isHadronicHP = true ;
  } else if (fHadronic == 1) {
    fHadPhysics = new G4HadronPhysicsQGSP_BERT_HP();
    isHadronicHP = true ;
    DSLog(routine) << "Hadronic Physics Active: QGSP_BERT_HP model" << endlog;
  } else if (fHadronic == 2) {
    fHadPhysics = new G4HadronPhysicsQGSP_BIC_HP();
    isHadronicHP = true ;
    DSLog(routine) << "Hadronic Physics Active: QSGP_BIC_HP model" << endlog;
  } else if (fHadronic == 3) {
    DSLog(routine) << "Hadronic Physics Active: FTF_BIC" << endlog;
    fHadPhysics = new G4HadronPhysicsFTF_BIC();
  } else if (fHadronic == 4) {
    DSLog(routine) << "Hadronic Physics Active: FTFP_BERT_HP model" << endlog;
    fHadPhysics = new G4HadronPhysicsFTFP_BERT_HP();
    isHadronicHP = true ;
  } else if (fHadronic == 5) {
    DSLog(routine) << "Hadronic Physics Active: Underground Shielding model" << endlog;
    fHadPhysics = new G4HadronPhysicsShielding("Shielding", true);
  } else if (fHadronic == 6)
    {
      DSLog(routine) << "Hadronic Physics Active: FPFP_BERT model (no HP!)" << endlog;
      fHadPhysics = new G4HadronPhysicsFTFP_BERT();
    }

  if (fHadPhysics)  // A valid choice was made for hadronic physics!
  {
    G4VPhysicsConstructor* hadElastic;
    if ( isHadronicHP == true )
      hadElastic = new G4HadronElasticPhysicsHP();
    else hadElastic = new G4HadronElasticPhysics();

    G4VPhysicsConstructor* ionPhysics = new G4IonPhysics();

    fHadPhysics->ConstructProcess();
    hadElastic->ConstructProcess();
    ionPhysics->ConstructProcess();

    DSLog(routine) << "Hadronic Elastic Processes ON" << endlog;
    DSLog(routine) << "Ion Physics List Active" << endlog;
  } else {
    DSLog(routine) << "Hadronic Physics DISABLED" << endlog;
  }

  // Decay
  G4VPhysicsConstructor* decayPhys = new G4DecayPhysics();
  decayPhys->ConstructProcess();

  // Photo-nuclear and lepto-nuclear
  G4VPhysicsConstructor* emextraPhys = new G4EmExtraPhysics();
  emextraPhys->ConstructProcess();

  // Radioactive Decay Physics
  G4VPhysicsConstructor* radDecayPhys = new G4RadioactiveDecayPhysics();
  radDecayPhys->ConstructProcess();

  // Step Limiter (2/3)
  G4VPhysicsConstructor* stepLimiterPhysics = new G4StepLimiterPhysics() ;
  stepLimiterPhysics->ConstructProcess();

  //  EnableAugerEffect();

  DumpProcessTable();

  OptionsForHadronicPhysics() ;
}

// Enable Auger Effect //////////////////////////////////////////////////////
#include "G4EmParameters.hh"
void DSPhysicsList::EnableAugerEffect() {
  DSLog(routine) << "Enabling Auger Effect" << endlog;
  G4EmParameters::Instance()->SetFluo(true);
  G4EmParameters::Instance()->SetAuger(true);
  G4EmParameters::Instance()->SetPixe(true);
}

// Options for the Hadronic Physics List
#include "G4ParticleHPManager.hh"
void DSPhysicsList::OptionsForHadronicPhysics() {
  DSLog(routine) << "Adding Options for Had Physics" << endlog;
  G4ParticleHPManager::GetInstance()->SetUseOnlyPhotoEvaporation( true );
  G4ParticleHPManager::GetInstance()->DumpSetting() ;
}

//
// DF:the following section is for liquid scintillators.
// I PROPOSE TO REMOVE IT

/*
#include "G4EmSaturation.hh"
#include "G4Scintillation.hh"
void DSPhysicsList::ConstructLSOp() {
  // G4 default scintillation process
  G4Scintillation* theScintProcessDef = new G4Scintillation("Scintillation");
  theScintProcessDef->SetTrackSecondariesFirst(true);
  theScintProcessDef->SetScintillationYieldFactor(1.0);      //
  theScintProcessDef->SetScintillationExcitationRatio(1.0);  //
  theScintProcessDef->SetVerboseLevel(OpVerbLevel);
  // Use Birks's correction in the scintillation process
  G4EmSaturation* emSaturation = G4LossTableManager::Instance()->EmSaturation();
  theScintProcessDef->AddSaturation(emSaturation);
  //theScintProcessDef->DumpPhysicsTable();

  auto theParticleIterator = GetParticleIterator();

  theParticleIterator->reset();
  while ((*theParticleIterator)()) {
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();

    if (theScintProcessDef->IsApplicable(*particle)) {
      pmanager->AddProcess(theScintProcessDef, ordDefault, ordInActive,
ordDefault); pmanager->SetProcessOrderingToLast(theScintProcessDef, idxAtRest);
      pmanager->SetProcessOrderingToLast(theScintProcessDef, idxPostStep);
    }
  }
}
*/

/*
//Custom model
// Elastic processes:
#include "G4ChipsElasticModel.hh"
#include "G4ElasticHadrNucleusHE.hh"
#include "G4HadronElasticProcess.hh"

// Inelastic processes:
#include "G4HadronInelasticProcess.hh"

#include "G4AlphaInelasticProcess.hh"
#include "G4AntiNeutronInelasticProcess.hh"
#include "G4AntiProtonInelasticProcess.hh"
#include "G4DeuteronInelasticProcess.hh"
#include "G4KaonMinusInelasticProcess.hh"
#include "G4KaonPlusInelasticProcess.hh"
#include "G4KaonZeroLInelasticProcess.hh"
#include "G4KaonZeroSInelasticProcess.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4PionMinusInelasticProcess.hh"
#include "G4PionPlusInelasticProcess.hh"
#include "G4ProtonInelasticProcess.hh"
#include "G4TritonInelasticProcess.hh"

//migration
// High energy FTFP model and Bertini cascade
#include "G4CascadeInterface.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4FTFModel.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4LundStringFragmentation.hh"
#include "G4PreCompoundModel.hh"
#include "G4TheoFSGenerator.hh"

// Cross sections
#include "G4AntiNuclElastic.hh"
#include "G4BGGNucleonInelasticXS.hh"
#include "G4BGGPionElasticXS.hh"
#include "G4ComponentAntiNuclNuclearXS.hh"
#include "G4CrossSectionDataSetRegistry.hh"
#include "G4CrossSectionElastic.hh"
#include "G4CrossSectionInelastic.hh"
#include "G4CrossSectionPairGG.hh"
#include "G4GGNuclNuclCrossSection.hh"
#include "G4HadronElastic.hh"
#include "G4NeutronCaptureProcess.hh"
#include "G4PiNuclearCrossSection.hh"
#include "G4VCrossSectionDataSet.hh"

// Neutron high-precision models: <20 MeV
#include "G4NeutronHPCapture.hh"
#include "G4NeutronHPCaptureData.hh"
#include "G4NeutronHPElastic.hh"
#include "G4NeutronHPElasticData.hh"
#include "G4NeutronHPInelastic.hh"
#include "G4NeutronHPInelasticData.hh"
//migration
//#include "G4LCapture.hh"

// Stopping processes
// migration
#include "G4AntiNeutronAnnihilationAtRest.hh"
#include "G4AntiProtonAbsorptionFritiof.hh"
#include "G4KaonMinusAbsorptionBertini.hh"
#include "G4PiMinusAbsorptionBertini.hh"


// ConstructHad()
// Makes discrete physics processes for the hadrons, at present limited
// to those particles with GHEISHA interactions (INTRC > 0).
// The processes are: Elastic scattering and Inelastic scattering.
// F.W.Jones  09-JUL-1998
void DSPhysicsList::ConstructHad_Custom()
{
//migration
//  G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
//  G4LElastic* theElasticModel = new G4LElastic;
//  theElasticProcess->RegisterMe(theElasticModel);

    //Elastic models
  const G4double elastic_elimitPi = 1.0*GeV;

  G4HadronElastic* elastic_lhep0 = new G4HadronElastic();
  G4HadronElastic* elastic_lhep1 = new G4HadronElastic();
  elastic_lhep1->SetMaxEnergy( elastic_elimitPi );
  G4ChipsElasticModel* elastic_chip = new G4ChipsElasticModel();
  G4ElasticHadrNucleusHE* elastic_he = new G4ElasticHadrNucleusHE();
  elastic_he->SetMinEnergy( elastic_elimitPi );


  // Inelastic scattering
  const G4double theFTFMin0 =    0.0*GeV;
  const G4double theFTFMin1 =    4.0*GeV;
  const G4double theFTFMax =   100.0*TeV;
  const G4double theBERTMin0 =   0.0*GeV;
  const G4double theBERTMin1 =  19.0*MeV;
  const G4double theBERTMax =    5.0*GeV;
  const G4double theHPMin =      0.0*GeV;
  const G4double theHPMax =     20.0*MeV;

  G4FTFModel * theStringModel = new G4FTFModel;
  G4ExcitedStringDecay * theStringDecay = new G4ExcitedStringDecay( new
G4LundStringFragmentation ); theStringModel->SetFragmentationModel(
theStringDecay ); G4PreCompoundModel * thePreEquilib = new G4PreCompoundModel(
new G4ExcitationHandler ); G4GeneratorPrecompoundInterface * theCascade = new
G4GeneratorPrecompoundInterface( thePreEquilib );

  G4TheoFSGenerator * theFTFModel0 = new G4TheoFSGenerator( "FTFP" );
  theFTFModel0->SetHighEnergyGenerator( theStringModel );
  theFTFModel0->SetTransport( theCascade );
  theFTFModel0->SetMinEnergy( theFTFMin0 );
  theFTFModel0->SetMaxEnergy( theFTFMax );

  G4TheoFSGenerator * theFTFModel1 = new G4TheoFSGenerator( "FTFP" );
  theFTFModel1->SetHighEnergyGenerator( theStringModel );
  theFTFModel1->SetTransport( theCascade );
  theFTFModel1->SetMinEnergy( theFTFMin1 );
  theFTFModel1->SetMaxEnergy( theFTFMax );

  G4CascadeInterface * theBERTModel0 = new G4CascadeInterface;
  theBERTModel0->SetMinEnergy( theBERTMin0 );
  theBERTModel0->SetMaxEnergy( theBERTMax );

  G4CascadeInterface * theBERTModel1 = new G4CascadeInterface;
  theBERTModel1->SetMinEnergy( theBERTMin1 );
  theBERTModel1->SetMaxEnergy( theBERTMax );

  G4VCrossSectionDataSet * thePiData = new G4CrossSectionPairGG( new
G4PiNuclearCrossSection, 91*GeV ); G4VCrossSectionDataSet * theAntiNucleonData =
new G4CrossSectionInelastic( new G4ComponentAntiNuclNuclearXS );
  G4VCrossSectionDataSet * theGGNuclNuclData =
G4CrossSectionDataSetRegistry::Instance()->
    GetCrossSectionDataSet(G4GGNuclNuclCrossSection::Default_Name());


  theParticleIterator->reset();
  while ((*theParticleIterator)())
    {
      G4ParticleDefinition* particle = theParticleIterator->value();
      G4ProcessManager* pmanager = particle->GetProcessManager();
      G4String particleName = particle->GetParticleName();

      if (particleName == "pi+")
        {
          // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new
G4HadronElasticProcess; theElasticProcess->AddDataSet( new G4BGGPionElasticXS(
particle ) ); theElasticProcess->RegisterMe( elastic_lhep1 );
          theElasticProcess->RegisterMe( elastic_he );
          pmanager->AddDiscreteProcess( theElasticProcess );
          //Inelastic scattering
          G4PionPlusInelasticProcess* theInelasticProcess =
            new G4PionPlusInelasticProcess("inelastic");
          theInelasticProcess->AddDataSet( thePiData );
          theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
          pmanager->AddDiscreteProcess( theInelasticProcess );
        }

      else if (particleName == "pi-")
        {
          // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new
G4HadronElasticProcess; theElasticProcess->AddDataSet( new G4BGGPionElasticXS(
particle ) ); theElasticProcess->RegisterMe( elastic_lhep1 );
          theElasticProcess->RegisterMe( elastic_he );
          pmanager->AddDiscreteProcess( theElasticProcess );
          //Inelastic scattering
          G4PionMinusInelasticProcess* theInelasticProcess =
            new G4PionMinusInelasticProcess("inelastic");
          theInelasticProcess->AddDataSet( thePiData );
          theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
          pmanager->AddDiscreteProcess( theInelasticProcess );
          //Absorption
          pmanager->AddRestProcess(new G4PiMinusAbsorptionBertini, ordDefault);
        }

      else if (particleName == "kaon+")
        {
          // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new
G4HadronElasticProcess; theElasticProcess->RegisterMe( elastic_lhep0 );
          pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering
          G4KaonPlusInelasticProcess* theInelasticProcess =
            new G4KaonPlusInelasticProcess("inelastic");
          theInelasticProcess->AddDataSet(
G4CrossSectionDataSetRegistry::Instance()->
                                           GetCrossSectionDataSet(G4ChipsKaonPlusInelasticXS::Default_Name()));
          theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
          pmanager->AddDiscreteProcess( theInelasticProcess );
        }

      else if (particleName == "kaon0S")
        {
          // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new
G4HadronElasticProcess; theElasticProcess->RegisterMe( elastic_lhep0 );
          pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering
          G4KaonZeroSInelasticProcess* theInelasticProcess =
            new G4KaonZeroSInelasticProcess("inelastic");
          theInelasticProcess->AddDataSet(
G4CrossSectionDataSetRegistry::Instance()->
                                           GetCrossSectionDataSet(G4ChipsKaonZeroInelasticXS::Default_Name()));
          theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
          pmanager->AddDiscreteProcess( theInelasticProcess );
        }

      else if (particleName == "kaon0L")
        {
          // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new
G4HadronElasticProcess; theElasticProcess->RegisterMe( elastic_lhep0 );
          pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering
          G4KaonZeroLInelasticProcess* theInelasticProcess =
            new G4KaonZeroLInelasticProcess("inelastic");
          theInelasticProcess->AddDataSet(
G4CrossSectionDataSetRegistry::Instance()->
                                           GetCrossSectionDataSet(G4ChipsKaonZeroInelasticXS::Default_Name()));
          theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
          pmanager->AddDiscreteProcess( theInelasticProcess );
        }

      else if (particleName == "kaon-")
        {
          // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new
G4HadronElasticProcess; theElasticProcess->RegisterMe( elastic_lhep0 );
          pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering
          G4KaonMinusInelasticProcess* theInelasticProcess =
            new G4KaonMinusInelasticProcess("inelastic");
          theInelasticProcess->AddDataSet(
G4CrossSectionDataSetRegistry::Instance()->
                                           GetCrossSectionDataSet(G4ChipsKaonMinusInelasticXS::Default_Name()));
          theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
          pmanager->AddDiscreteProcess( theInelasticProcess );
          pmanager->AddRestProcess(new G4KaonMinusAbsorptionBertini,
ordDefault);
        }

      else if (particleName == "proton")
        {
          // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new
G4HadronElasticProcess;
          theElasticProcess->AddDataSet(G4CrossSectionDataSetRegistry::Instance()->
                                        GetCrossSectionDataSet(G4ChipsProtonElasticXS::Default_Name()));
          theElasticProcess->RegisterMe( elastic_chip );
          pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering
          G4ProtonInelasticProcess* theInelasticProcess =
            new G4ProtonInelasticProcess("inelastic");
          theInelasticProcess->AddDataSet( new G4BGGNucleonInelasticXS(
G4Proton::Proton() ) ); theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
          pmanager->AddDiscreteProcess( theInelasticProcess );
        }
      else if (particleName == "anti_proton")
        {
          // Elastic scattering
          const G4double elastic_elimitAntiNuc = 100.0*CLHEP::MeV;
          G4AntiNuclElastic* elastic_anuc = new G4AntiNuclElastic();
          elastic_anuc->SetMinEnergy( elastic_elimitAntiNuc );
          G4CrossSectionElastic* elastic_anucxs = new G4CrossSectionElastic(
elastic_anuc->GetComponentCrossSection() ); G4HadronElastic* elastic_lhep2 = new
G4HadronElastic(); elastic_lhep2->SetMaxEnergy( elastic_elimitAntiNuc );
          G4HadronElasticProcess* theElasticProcess = new
G4HadronElasticProcess; theElasticProcess->AddDataSet( elastic_anucxs );
          theElasticProcess->RegisterMe( elastic_lhep2 );
          theElasticProcess->RegisterMe( elastic_anuc );
          pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering
          G4AntiProtonInelasticProcess* theInelasticProcess =
            new G4AntiProtonInelasticProcess("inelastic");
          theInelasticProcess->AddDataSet( theAntiNucleonData );
          theInelasticProcess->RegisterMe( theFTFModel0 );
          pmanager->AddDiscreteProcess( theInelasticProcess );
          // Absorption
          pmanager->AddRestProcess(new G4AntiProtonAbsorptionFritiof,
ordDefault);
        }

      else if (particleName == "neutron") {
        // elastic scattering
        G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
        theElasticProcess->AddDataSet(G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4ChipsNeutronElasticXS::Default_Name()));
        G4HadronElastic* elastic_neutronChipsModel = new G4ChipsElasticModel();
        elastic_neutronChipsModel->SetMinEnergy( 19.0*CLHEP::MeV );
        theElasticProcess->RegisterMe( elastic_neutronChipsModel );
        G4NeutronHPElastic * theElasticNeutronHP = new G4NeutronHPElastic;
        theElasticNeutronHP->SetMinEnergy( theHPMin );
        theElasticNeutronHP->SetMaxEnergy( theHPMax );
        theElasticProcess->RegisterMe( theElasticNeutronHP );
        theElasticProcess->AddDataSet( new G4NeutronHPElasticData );
        pmanager->AddDiscreteProcess( theElasticProcess );
        // inelastic scattering
        G4NeutronInelasticProcess* theInelasticProcess =
          new G4NeutronInelasticProcess("inelastic");
        theInelasticProcess->AddDataSet( new G4BGGNucleonInelasticXS(
G4Neutron::Neutron() ) ); theInelasticProcess->RegisterMe( theFTFModel1 );
        theInelasticProcess->RegisterMe( theBERTModel1 );
        G4NeutronHPInelastic * theNeutronInelasticHPModel = new
G4NeutronHPInelastic; theNeutronInelasticHPModel->SetMinEnergy( theHPMin );
        theNeutronInelasticHPModel->SetMaxEnergy( theHPMax );
        theInelasticProcess->RegisterMe( theNeutronInelasticHPModel );
        theInelasticProcess->AddDataSet( new G4NeutronHPInelasticData );
        pmanager->AddDiscreteProcess(theInelasticProcess);
        // capture
        G4NeutronCaptureProcess* theCaptureProcess =
          new G4NeutronCaptureProcess;
        G4NeutronHPCapture * theLENeutronCaptureModel = new G4NeutronHPCapture;
        theLENeutronCaptureModel->SetMinEnergy(theHPMin);
        theLENeutronCaptureModel->SetMaxEnergy(theHPMax);
        theCaptureProcess->RegisterMe(theLENeutronCaptureModel);
        theCaptureProcess->AddDataSet( new G4NeutronHPCaptureData);
        pmanager->AddDiscreteProcess(theCaptureProcess);

      }
      else if (particleName == "anti_neutron")
        {
          // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new
G4HadronElasticProcess; theElasticProcess->RegisterMe( elastic_lhep0 );
          pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering (include annihilation on-fly)
          G4AntiNeutronInelasticProcess* theInelasticProcess =
            new G4AntiNeutronInelasticProcess("inelastic");
          theInelasticProcess->AddDataSet( theAntiNucleonData );
          theInelasticProcess->RegisterMe( theFTFModel0 );
          pmanager->AddDiscreteProcess( theInelasticProcess );
        }

      else if (particleName == "deuteron")
        {
          // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new
G4HadronElasticProcess; theElasticProcess->RegisterMe( elastic_lhep0 );
          pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering
          G4DeuteronInelasticProcess* theInelasticProcess =
            new G4DeuteronInelasticProcess("inelastic");
          theInelasticProcess->AddDataSet( theGGNuclNuclData );
          theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
          pmanager->AddDiscreteProcess( theInelasticProcess );
        }

      else if (particleName == "triton")
        {
          // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new
G4HadronElasticProcess; theElasticProcess->RegisterMe( elastic_lhep0 );
          pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering
          G4TritonInelasticProcess* theInelasticProcess =
            new G4TritonInelasticProcess("inelastic");
          theInelasticProcess->AddDataSet( theGGNuclNuclData );
          theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
          pmanager->AddDiscreteProcess( theInelasticProcess );
        }
      else if (particleName == "alpha")
        {
          // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new
G4HadronElasticProcess; theElasticProcess->RegisterMe( elastic_lhep0 );
          pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering
          G4AlphaInelasticProcess* theInelasticProcess =
            new G4AlphaInelasticProcess("inelastic");
          theInelasticProcess->AddDataSet( theGGNuclNuclData );
          theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
          pmanager->AddDiscreteProcess( theInelasticProcess );
        }
    }
}
*/

// Cuts /////////////////////////////////////////////////////////////////////
void DSPhysicsList::SetCuts() {

  if (isDepositCuts) {
    DSLog(routine) << "Default Physics Cut for the DSGenneratorEnergyDeposit: 10 cm" << endlog;
    defaultCutValue = 10.0 * cm;  //
  }

  cutForGamma = defaultCutValue;
  cutForElectron = defaultCutValue;  // 1.0*nanometer;
  cutForPositron = defaultCutValue;

  DSLog(routine) << " Default Physics Cut Set to " << defaultCutValue / cm << " cm " << endlog;

  // special for low energy physics
  //   G4double lowlimit=250*eV;

  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma
  SetCutValue(cutForGamma, "gamma");
  SetCutValue(cutForElectron, "e-");
  SetCutValue(cutForPositron, "e+");

  SetCutValue(1 * cm, "opticalphoton");
  SetCutValue(defaultCutValue, "neutron");
  SetCutValue(1 * micrometer, "GenericIon");

  if (isDepositCuts) SetCutValue(10 * cm, "GenericIon");

  G4double lowlimit = 1 * keV;
  G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(lowlimit, 100. * GeV);

  // Set Range Cut to 1 mm for gammas, electrons and positrons
  SetCutsWithDefault();

  // If enabled, set the range cuts to 1 micron in the active LAr volume
  if (isHPRangeCuts) {
    double HPRangeCutEle = 1 * cm;
    double HPRangeCut = 1 * mm;
    double HPRangeCutNP = 1 * um;
    if (isDepositCuts) HPRangeCut = 10 * cm;

    G4Region* activeLArRegion = G4RegionStore::GetInstance()->GetRegion("LAr_Logic");
    G4ProductionCuts* cutsForLAr = new G4ProductionCuts;
    cutsForLAr->SetProductionCut(HPRangeCut, "gamma");
    cutsForLAr->SetProductionCut(HPRangeCutEle, "e-");
    cutsForLAr->SetProductionCut(HPRangeCutEle, "e+");
    cutsForLAr->SetProductionCut(HPRangeCutNP, "proton");
    cutsForLAr->SetProductionCut(HPRangeCutNP, "neutron");
    cutsForLAr->SetProductionCut(HPRangeCutNP, "GenericIon");
    activeLArRegion->SetProductionCuts(cutsForLAr);
  }

  if (verboseLevel > 0) DumpCutValuesTable();
}

void DSPhysicsList::DumpProcessTable() {

  int processmap[10000];
  for (int i = 0; i < 10000; i++) { processmap[i] = 0; }

  auto theParticleIterator = GetParticleIterator();

  theParticleIterator->reset();
  while ((*theParticleIterator)()) {
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();

    G4ProcessVector* temp = pmanager->GetProcessList();
    for (int i = 0; i < pmanager->GetProcessListLength(); i++) {
      int index = 1000 * (*temp)[i]->GetProcessType() + (*temp)[i]->GetProcessSubType();
      if (!processmap[index]) {
        processmap[index] = 1;
        DSLog(routine) << " Process Map Table - Index: " << index << " ; Process Name: " << (*temp)[i]->GetProcessName() << endlog;
      }
    }
  }
}
