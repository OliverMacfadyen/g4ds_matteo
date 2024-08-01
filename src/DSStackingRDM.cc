#include "DSStackingRDM.hh"
#include "CLHEP/Random/RandExponential.h"
#include "DSEventHandler.hh"
#include "DSLogger.hh"
#include "DSStackingRDMMessenger.hh"
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

DSStackingRDM::DSStackingRDM() {

  IsReclassify = false;
  fMessenger = new DSStackingRDMMessenger(this);
  fLevelZero = false;
  fPreTrack = 0;
  DSLog(routine) << "RDM Stacking Methode Active" << endlog;
}

DSStackingRDM::~DSStackingRDM() {
  ;
}

G4ClassificationOfNewTrack DSStackingRDM::DSClassifyNewTrack(const G4Track* aTrack) {
  // skip optical photons
  
  if (aTrack->GetDefinition()->GetPDGEncoding() == -22) return fUrgent;
  /*cout << "stack "
       << aTrack->GetDefinition()->GetParticleName() << " "
       << aTrack->GetGlobalTime()/ns << " "
       << aTrack->GetTrackID() << " "
       << aTrack->GetParentID() << " "
       << endl ;
  */
  // kill neutrinos
  if (aTrack->GetDefinition() == G4NeutrinoE::NeutrinoE()) return fKill;
  if (aTrack->GetDefinition() == G4NeutrinoMu::NeutrinoMu()) return fKill;
  if (aTrack->GetDefinition() == G4NeutrinoTau::NeutrinoTau()) return fKill;
  if (aTrack->GetDefinition() == G4AntiNeutrinoE::AntiNeutrinoE()) return fKill;
  if (aTrack->GetDefinition() == G4AntiNeutrinoMu::AntiNeutrinoMu()) return fKill;
  if (aTrack->GetDefinition() == G4AntiNeutrinoTau::AntiNeutrinoTau()) return fKill;




  // set lifetime of primary to 0 to oblige the decay 
  if (aTrack->GetDefinition()->GetAtomicNumber() > 2) {
    if (aTrack->GetParentID() == 0) { 
      aTrack->GetDefinition()->SetPDGLifeTime(0. * ns);
      fPIDs.push_back(aTrack->GetTrackID());
      return fUrgent ;
    } else {
      G4Ions* ion = (G4Ions*)aTrack->GetDefinition();
      fPIDs.push_back(aTrack->GetTrackID());
      if (ion->GetExcitationEnergy() == 0) return fKill;
    }
  } else {
    for (std::vector<int>::iterator it = fPIDs.begin(); it != fPIDs.end(); ++it) {
      if (aTrack->GetParentID() == *it) {
        

        // kill particles above a certain energy as defined by the user
        if((int)fPDGToBeKilled.size() >0) {
          for (G4int i = 0; i < G4int(fPDGToBeKilled.size()); i++)
            if (aTrack->GetDefinition()->GetPDGEncoding() == fPDGToBeKilled[i] && aTrack->GetKineticEnergy() < fLEnergyToBeKilled[i]) return fKill;
        }


        DSEventHandler::Get()->SetDId(int(DSEventHandler::Get()->GetVDaughters().size()));
        DSEventHandler::Get()->SetDPID(aTrack->GetParentID());
        DSEventHandler::Get()->SetDPDG(aTrack->GetDefinition()->GetPDGEncoding());
        DSEventHandler::Get()->SetDProcess(-1);
        DSEventHandler::Get()->SetDTime(aTrack->GetGlobalTime() / ns);
        DSEventHandler::Get()->SetDEnergy(aTrack->GetKineticEnergy() / keV);
        DSEventHandler::Get()->SetDPosition(aTrack->GetVertexPosition() / cm);
        DSEventHandler::Get()->SetDDirection(aTrack->GetVertexMomentumDirection());
        DSEventHandler::Get()->SetDaughters();
        break;
      } 
    }
  }
  

  return fUrgent;
}

void DSStackingRDM::DSNewStage() {
  DSStackClearAll();
  return;
}

void DSStackingRDM::DSPrepareNewEvent() {
  fCounter = 0;
  IsValid = false;
  stage = -1;
  fDaughters = 0;
  fAlpha = 0;
  fBeta = 0;
  fGamma = 0;

  fPIDs.clear();
  // cout << "------ NEW ------- " << endl ;
  // IsTheFirstDaughter = false ;
}

/*
 * Revision 1.1  2013/03/22 13:50:09  dfranco
 * Added radioactive decay and radioactive decay chain generators. They require
 * the correspondent stacking actions. Two mac files are included as examples
 *
 *
 */
