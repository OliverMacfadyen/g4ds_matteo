#include "DSStackingNeutronCapture.hh"
#include "CLHEP/Random/RandExponential.h"
#include "DSEventHandler.hh"
#include "DSLogger.hh"
#include "DSStorage.hh"
#include "G4AntiNeutrinoE.hh"
#include "G4AntiNeutrinoMu.hh"
#include "G4AntiNeutrinoTau.hh"
#include "G4GenericIon.hh"
#include "G4NeutrinoE.hh"
#include "G4NeutrinoMu.hh"
#include "G4NeutrinoTau.hh"
#include "G4OpticalPhoton.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4StackManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4ios.hh"
#include "globals.hh"

using namespace std;

DSStackingNeutronCapture::DSStackingNeutronCapture() {}

DSStackingNeutronCapture::~DSStackingNeutronCapture() {
  ;
}


 
G4ClassificationOfNewTrack DSStackingNeutronCapture::DSClassifyNewTrack(const G4Track* aTrack) {
  // neutron capture creator model  ID + Index = 25000 + 84
  int neutron_capture = 25084;
  if (aTrack->GetGlobalTime() > 100 * ms) return fKill;
  int process_id = aTrack->GetCreatorModelID () + aTrack->GetCreatorModelIndex (); 
  //int pdg        = aTrack->GetDefinition()->GetPDGEncoding();

  if (process_id == neutron_capture) {
    DSEventHandler::Get()->SetDepPID(aTrack->GetDefinition()->GetPDGEncoding());
    DSEventHandler::Get()->SetDepVolume(aTrack->GetVolume()->GetLogicalVolume()->GetMaterial()->GetIndex());
    DSEventHandler::Get()->SetDepEnergy(aTrack->GetKineticEnergy() / keV);
    DSEventHandler::Get()->SetDepStep(-1000);
    DSEventHandler::Get()->SetDepTime(aTrack->GetGlobalTime() / ns);
    DSEventHandler::Get()->SetDepPosition(aTrack->GetPosition() / cm);
    DSEventHandler::Get()->SetDeposits();
    return fKill;
  }

  return fUrgent;
}

void DSStackingNeutronCapture::DSNewStage() {}

void DSStackingNeutronCapture::DSPrepareNewEvent() {}


