#include "DSStackingRDMChain.hh"
#include "CLHEP/Random/RandExponential.h"
#include "DSEventHandler.hh"
#include "DSLogger.hh"
#include "DSStackingRDMChainMessenger.hh"
#include "DSStorage.hh"
#include "G4AntiNeutrinoE.hh"
#include "G4AntiNeutrinoMu.hh"
#include "G4AntiNeutrinoTau.hh"
#include "G4GenericIon.hh"
#include "G4Ions.hh"
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

DSStackingRDMChain::DSStackingRDMChain() {
  changedParticle = new G4ParticleChange;
  fStopUntilIsotope = false ;
  fStopIsotope = 0 ;
  IsShort = false;
  fCounter = 0;
  DSStorage::Get()->SetRDMDecay(true);
  IsReclassify = false;
  fMaxLifeTime = 1.e30 * s;  // 10*24*3600*s ;
  isDaughter = false;
  fMessenger = new DSStackingRDMChainMessenger(this);
  isPostponed = false;
  isFirst = true;
  fGateTime = 100. * ns;

  //cout << "Stop at " << fStopIsotope << " " << fStopUntilIsotope << endl ;
  DSLog(routine) << "RDMChain Stacking Methode Active" << endlog;
}

DSStackingRDMChain::~DSStackingRDMChain() {
  ;
}

G4ClassificationOfNewTrack DSStackingRDMChain::DSClassifyNewTrack(const G4Track* aTrack) {

  if (fMaxLifeTime == 0) fMaxLifeTime = 1.e30 * day;

  // skip optical photons
  if (aTrack->GetDefinition()->GetPDGEncoding() == -22) return fUrgent;

  // kill neutrinos
  if (aTrack->GetDefinition() == G4NeutrinoE::NeutrinoE()) return fKill;
  if (aTrack->GetDefinition() == G4NeutrinoMu::NeutrinoMu()) return fKill;
  if (aTrack->GetDefinition() == G4NeutrinoTau::NeutrinoTau()) return fKill;
  if (aTrack->GetDefinition() == G4AntiNeutrinoE::AntiNeutrinoE()) return fKill;
  if (aTrack->GetDefinition() == G4AntiNeutrinoMu::AntiNeutrinoMu()) return fKill;
  if (aTrack->GetDefinition() == G4AntiNeutrinoTau::AntiNeutrinoTau()) return fKill;

  /*
  if (aTrack->GetDefinition()->GetPDGEncoding() > 1e6 &&
  aTrack->GetDefinition()->GetAtomicNumber () > 2) { G4Ions *ion =
  (G4Ions*)aTrack->GetDefinition(); cout <<
  aTrack->GetDefinition()->GetParticleName() << " "
    << aTrack->GetDefinition()->GetAtomicNumber() << " "
    <<  aTrack->GetDefinition()->GetAtomicMass() << " "
    <<  aTrack->GetParentID() << " "
    <<  "LT " << aTrack->GetDefinition()->GetPDGLifeTime()/s << " "
    << ion->GetExcitationEnergy ()
    <<  " ID: " <<aTrack->GetTrackID() << " "
    <<  " PID: " <<aTrack->GetParentID() << " "
    <<  " gtime: " << aTrack->GetGlobalTime()/s << " "
    <<  " life: " << aTrack->GetDefinition()->GetPDGLifeTime()/s<< " "
    <<  ion->GetExcitationEnergy () << " "
    << aTrack->GetDefinition()->GetPDGEncoding() << " "
    << mymap[aTrack->GetDefinition()->GetPDGEncoding()]
    <<  endl ;
  }*/

  //cout << aTrack->GetDefinition()->GetParticleName() << endl ;
  if (aTrack->GetDefinition()->GetPDGEncoding() > 1e6 && aTrack->GetDefinition()->GetAtomicNumber() > 2) {
    if (aTrack->GetDefinition()->GetPDGLifeTime() > 0)
      mymap[aTrack->GetDefinition()->GetPDGEncoding()] = aTrack->GetDefinition()->GetPDGLifeTime();
          // cout <<aTrack->GetDefinition()->GetPDGEncoding() << " "
          //      << aTrack->GetDefinition()->GetParticleName() << " "
          //      << aTrack->GetDefinition()->GetPDGLifeTime() << " "
          //      << endl ;

    int pid = aTrack->GetParentID();
    // int tid = aTrack->GetTrackID();
    G4Ions* ion = (G4Ions*)aTrack->GetDefinition();


    //cout << aTrack->GetDefinition()->GetParticleName() << " " << pid << " "  << endl ;

    // If stable and no excited level, kill the particle
    if (aTrack->GetDefinition()->GetPDGStable() && ion->GetExcitationEnergy() == 0) {
      // cout << "Uccido lo stabile: " <<
      // aTrack->GetDefinition()->GetParticleName() <<  endl ;
      isFirst = true;
      return fKill;
    }

    if (fStopUntilIsotope && aTrack->GetDefinition()->GetPDGEncoding()  == fStopIsotope) {
      //cout << aTrack->GetDefinition()->GetPDGEncoding() << endl ;
      isFirst = true;
      return fKill;
    }

    // If lifetime longer than threshold, kill it
    if (mymap[aTrack->GetDefinition()->GetPDGEncoding()] > fMaxLifeTime) {
      // cout << "Uccido il nucleo che vive decisamente troppo: " <<
      // aTrack->GetDefinition()->GetParticleName() <<  endl ;
      isFirst = true;
      return fKill;
    }

    if (pid == 0) {
      if (isFirst) {
        // cout << "lancio il padre " << pid << " " <<
        // aTrack->GetDefinition()->GetParticleName() <<endl ;
        aTrack->GetDefinition()->SetPDGLifeTime(0 * s);
        isFirst = false;
        DSEventHandler::Get()->SetPDG(aTrack->GetDefinition()->GetPDGEncoding());
        DSEventHandler::Get()->SetTime(CLHEP::RandExponential::shoot(mymap[aTrack->GetDefinition()->GetPDGEncoding()]));
        DSEventHandler::Get()->SetPosition(aTrack->GetPosition() / cm);
        DSEventHandler::Get()->SetDirection(aTrack->GetMomentumDirection());
        return fUrgent;
      } else {
        // cout << "Uccido il capostipite " <<endl ;
        if (fCounter > 2) {  // protezione contro strani isotopi eccitati che non decadono
          // cout << "Salvo il capostipite" << endl ;
          DSEventHandler::Get()->SetPDG(aTrack->GetDefinition()->GetPDGEncoding());
          DSEventHandler::Get()->SetTime(CLHEP::RandExponential::shoot(mymap[aTrack->GetDefinition()->GetPDGEncoding()]));
          DSEventHandler::Get()->SetPosition(aTrack->GetPosition() / cm);
          DSEventHandler::Get()->SetDirection(aTrack->GetMomentumDirection());

          return fUrgent;
        }
        fCounter++;
        return fKill;
      }
    }

    if (pid == -1 && !fDaughters) {
      fCounter = 0;
      // cout << "Lancio il figlio " << pid <<  " "
      // << aTrack->GetDefinition()->GetParticleName() << endl ;
      fDaughters = 1;
      DSEventHandler::Get()->SetPDG(aTrack->GetDefinition()->GetPDGEncoding());
      DSEventHandler::Get()->SetTime(CLHEP::RandExponential::shoot(mymap[aTrack->GetDefinition()->GetPDGEncoding()]));
      DSEventHandler::Get()->SetPosition(aTrack->GetPosition() / cm);
      DSEventHandler::Get()->SetDirection(aTrack->GetMomentumDirection());
      return fUrgent;
    }

    if (mymap[aTrack->GetDefinition()->GetPDGEncoding()] > fGateTime) {
      // cout << "Salvo il figlio " <<
      // aTrack->GetDefinition()->GetParticleName() << endl ;
      aTrack->GetDefinition()->SetPDGLifeTime(0 * s);
      isPostponed = true;
      return fPostpone;
    }

    if (isPostponed) {
      // cout << "Uccido il figlio del posticipato  " <<
      // aTrack->GetDefinition()->GetParticleName() << endl ;
      return fKill;
    }
  }

  // cout  << "Stacking Post "
  //   <<   aTrack->GetDefinition()->GetParticleName() << " "
  //   <<  " gtime: " << aTrack->GetGlobalTime() << " "
  //   <<  " ptime: " << aTrack->GetProperTime() << " "
  //   <<  " ltime: " << aTrack->GetLocalTime()/ns  << " "
  //   <<  " steps: " << aTrack->GetCurrentStepNumber() << " "
  //   <<  " E: " <<aTrack->GetKineticEnergy()/keV<< " "
  //   <<  " ID: " <<aTrack->GetTrackID() << " "
  //   <<  " PID: " <<aTrack->GetParentID() << " "
  //   << endl ;

  if (aTrack->GetParentID() <= 2 && aTrack->GetDefinition()->GetAtomicNumber() <= 4) {
    DSEventHandler::Get()->SetDId(int(DSEventHandler::Get()->GetVDaughters().size()));
    DSEventHandler::Get()->SetDPDG(aTrack->GetDefinition()->GetPDGEncoding());
    DSEventHandler::Get()->SetDProcess(0);
    DSEventHandler::Get()->SetDTime(aTrack->GetGlobalTime() / ns);
    DSEventHandler::Get()->SetDEnergy(aTrack->GetKineticEnergy() / keV);
    DSEventHandler::Get()->SetDPosition(aTrack->GetVertexPosition() / cm);
    DSEventHandler::Get()->SetDDirection(aTrack->GetVertexMomentumDirection());
    DSEventHandler::Get()->SetDaughters();
  }

  return fUrgent;
}

void DSStackingRDMChain::DSNewStage() {

  /* candidates
  DSStackAbort()
  DSStackClearAll()
  DSStackClearUrgentAndWaiting()
  DSStackClearWaiting()
  DSStackClearUrgent()
  DSStackClearPostponed()
  DSStackCheckStatus()
  DSStackReClassify()
  */
  // DSStackCheckStatus();
  // DSStackReClassify();
  // DSStackCheckStatus();

  return;
}

void DSStackingRDMChain::DSPrepareNewEvent() {
  // fCounter         = 0;
  IsValid = false;
  stage = 0;
  isSecondDaughter = false;
  fDaughters = 0;
  isPostponed = false;
  // isFirst          = false ;
  fAlpha = 0;
  fBeta = 0;
  fGamma = 0;
  fPIDs.clear();
}
