#include "DSOpticalPhysics.hh"
#include "G4OpticalPhoton.hh"
#include "G4ProcessManager.hh"

#include "DSLight3.hh"
#include "DSOpWLS.hh"
#include "G4Cerenkov.hh"
#include "G4LossTableManager.hh"
#include "G4OpAbsorption.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4OpRayleigh.hh"
#include "G4OpWLS.hh"
#include "DSMatScintillation.hh"

#include "DSStorage.hh"

DSOpticalPhysics::DSOpticalPhysics(G4int verbose, const G4String& name) : G4VPhysicsConstructor(name) {
  verboseLevel = verbose;
}

DSOpticalPhysics::~DSOpticalPhysics() {
  ;
}

void DSOpticalPhysics::ConstructParticle() {
  // optical photon
  G4OpticalPhoton::OpticalPhotonDefinition();
}

void DSOpticalPhysics::ConstructProcess() {
  DSLight3* theScintProcessDefS1S2 = new DSLight3();

  // scintillation process for alpha:
  DSLight3* theScintProcessAlphaS1S2 = new DSLight3();

  // scintillation process for heavy nuclei
  DSLight3* theScintProcessNucS1S2 = new DSLight3();

  G4Cerenkov* theCerenkovProcessDef = new G4Cerenkov("Cherenkov");
  G4int MaxNumPhotons = 1000000;
  theCerenkovProcessDef->SetTrackSecondariesFirst(true);
  theCerenkovProcessDef->SetMaxNumPhotonsPerStep(MaxNumPhotons);


  DSMatScintillation* theScintillationProcessDef = new DSMatScintillation("MatScintillation");
  theScintillationProcessDef->SetTrackSecondariesFirst(true);
  theScintillationProcessDef->SetVerboseLevel(1);
  theScintillationProcessDef->DumpPhysicsTable();

  
  // default scintillation process
  //   theScintProcessDef->SetTrackSecondariesFirst(true);
  //   theScintProcessDef->SetScintillationYieldFactor(1.0); //
  //   theScintProcessDef->SetScintillationExcitationRatio(1.0); //
  // theScintProcessDef->SetVerboseLevel(OpVerbLevel);
  // Use Birks's correction in the scintillation process
  // G4EmSaturation* emSaturation =
  // G4LossTableManager::Instance()->EmSaturation();
  // theScintProcessDef->AddSaturation(emSaturation);
  // theScintProcessDef->DumpPhysicsTable();

  // optical processes
  G4OpAbsorption* theAbsorptionProcess = new G4OpAbsorption();
  G4OpRayleigh* theRayleighScatteringProcess = new G4OpRayleigh();
  G4OpBoundaryProcess* theBoundaryProcess = new G4OpBoundaryProcess();
  DSOpWLS* theWLSProcess = new DSOpWLS();
  //  theAbsorptionProcess->DumpPhysicsTable();
  //  theRayleighScatteringProcess->DumpPhysicsTable();
  theAbsorptionProcess->SetVerboseLevel(verboseLevel);
  // theRayleighScatteringProcess->SetVerboseLevel(OpVerbLevel);
  // theBoundaryProcess->SetVerboseLevel(OpVerbLevel);
  // G4OpticalSurfaceModel themodel = unified;
  // theBoundaryProcess->SetModel(themodel);
  theWLSProcess->UseTimeProfile("exponential");

  auto theParticleIterator = GetParticleIterator();

  theParticleIterator->reset();
  while ((*theParticleIterator)()) {
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();

    if (theScintProcessDefS1S2->IsApplicable(*particle)) {
      if (particle->GetParticleName() == "GenericIon") {
        pmanager->AddProcess(theScintProcessNucS1S2, ordDefault, ordInActive, ordDefault);  // AtRestDiscrete
        pmanager->SetProcessOrderingToLast(theScintProcessNucS1S2, idxAtRest);
        pmanager->SetProcessOrderingToLast(theScintProcessNucS1S2, idxPostStep);
        //pmanager->SetProcessOrdering(theScintProcessNucS1S2, idxPostStep);
     } else if (particle->GetParticleName() == "alpha") {
        pmanager->AddProcess(theScintProcessAlphaS1S2, ordDefault, ordInActive, ordDefault);
        pmanager->SetProcessOrderingToLast(theScintProcessAlphaS1S2, idxAtRest);
        pmanager->SetProcessOrderingToLast(theScintProcessAlphaS1S2, idxPostStep);
        //pmanager->SetProcessOrdering(theScintProcessAlphaS1S2, idxPostStep);
      } else {
        pmanager->AddProcess(theScintProcessDefS1S2, ordDefault, ordInActive, ordDefault);
        pmanager->SetProcessOrderingToLast(theScintProcessDefS1S2, idxAtRest);
        pmanager->SetProcessOrderingToLast(theScintProcessDefS1S2, idxPostStep);
        //pmanager->SetProcessOrdering(theScintProcessDefS1S2, idxPostStep);
      }
    }

    if (DSStorage::Get()->GetIsCherenkov() && theCerenkovProcessDef->IsApplicable(*particle)) {
      pmanager->AddProcess(theCerenkovProcessDef);
      pmanager->SetProcessOrdering(theCerenkovProcessDef, idxPostStep);
    }


    if (DSStorage::Get()->GetMatScintillation() && theScintillationProcessDef->IsApplicable(*particle)) {
      pmanager->AddProcess(theScintillationProcessDef);
      pmanager->SetProcessOrdering(theScintillationProcessDef, idxPostStep);
      //pmanager->SetProcessOrderingToLast(theScintillationProcessDef, idxPostStep);
      //pmanager->SetProcessOrdering(theScintillationProcessDef, idxPostStep);
    }

    if (particle == G4OpticalPhoton::Definition()) {
      pmanager->AddDiscreteProcess(theAbsorptionProcess);
      pmanager->AddDiscreteProcess(theRayleighScatteringProcess);
      pmanager->AddDiscreteProcess(theWLSProcess);
      pmanager->AddDiscreteProcess(theBoundaryProcess);
    }
  }
}
