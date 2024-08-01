#include "DSStackingOnlyOneParticle.hh"
#include "DSStackingOnlyOneParticleMessenger.hh"
#include "CLHEP/Random/RandExponential.h"
#include "DSEventHandler.hh"
#include "DSLogger.hh"
#include "DSStorage.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4StackManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4ios.hh"
#include "globals.hh"

using namespace std;



DSStackingOnlyOneParticle::DSStackingOnlyOneParticle() {
  fParticleTable = G4ParticleTable::GetParticleTable();
  fPDG = 2112 ;

  fMessenger = new DSStackingOnlyOneParticleMessenger(this);

  DSLog(routine) << "Kill particles if not " <<  fParticleTable->FindParticle(fPDG)->GetParticleName() << endlog;
  

}

DSStackingOnlyOneParticle::~DSStackingOnlyOneParticle() {
  ;
}


 
G4ClassificationOfNewTrack DSStackingOnlyOneParticle::DSClassifyNewTrack(const G4Track* aTrack) {
  if (aTrack->GetDefinition()->GetPDGEncoding() != fPDG) return fKill;

  return fUrgent;
}

void DSStackingOnlyOneParticle::DSNewStage() {}

void DSStackingOnlyOneParticle::DSPrepareNewEvent() {}

void DSStackingOnlyOneParticle::SetParticleName(string particle) {
  fPDG = fParticleTable->FindParticle(particle)->GetPDGEncoding();
}
