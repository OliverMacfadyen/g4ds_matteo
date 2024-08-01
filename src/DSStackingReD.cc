#include "DSStackingReD.hh"
#include "DSEventAction.hh"
#include "DSLogger.hh"

#include "G4Ions.hh"
#include "G4RunManager.hh"
#include "G4Track.hh"

//#include "DSStorage.hh"
//#include "G4SystemOfUnits.hh"
//#include "G4PhysicalConstants.hh"
//#include "DSEventHandler.hh"
//#include "G4StackManager.hh"
//#include "G4ios.hh"
//#include "G4ParticleTable.hh"
//#include "globals.hh"

//#include "G4VProcess.hh"
//#include "Primary.hh"

DSStackingReD::DSStackingReD() {  // :
                                  // fEventAction(G4RunManager::GetRunManager()->GetUserEventAction())
                                  // {
  DSLog(routine) << "ReD Stacking Methode Active" << endlog;
}

DSStackingReD::~DSStackingReD() {}

G4ClassificationOfNewTrack DSStackingReD::DSClassifyNewTrack(const G4Track* aTrack) {

  const G4ParticleDefinition* pDef = aTrack->GetParticleDefinition();
  const G4int pdg = pDef->GetPDGEncoding();
  //  if ( std::abs(pdg) < 23 ) // e-, e+, gamma
  //    return fUrgent;

  const G4int parentID = aTrack->GetParentID();
  const G4String name = pDef->GetParticleName();
  // const G4double T0 = aTrack->GetVertexKineticEnergy();
  const G4double T1 = aTrack->GetKineticEnergy();
  const G4ThreeVector pos0(aTrack->GetVertexPosition());
  const G4ThreeVector pos1(aTrack->GetPosition());
  const G4ThreeVector dir0(aTrack->GetVertexMomentumDirection());
  const G4ThreeVector dir1(aTrack->GetMomentumDirection());
  const G4VPhysicalVolume* volume(aTrack->GetVolume());
  // I could also cut on parentID != 0
  const G4String volName = volume ? volume->GetName() : "null";

  DSLog(routine) << "parentID " << parentID << "  " << name << "  "
                 << pdg
                 //                 << "    vertex: energy " << T0 << "  pos "
                 //                 << pos0 << "  dir " << dir0
                 << "    novertex: energy " << T1 << "  pos " << pos1 << "  dir " << dir1 << "    volume " << volName << endlog;

  // Attention: neutrons AFTER elastic scattering still keep the ParentID() = 0
  if (parentID == 0) {
    ; /*
    DSLog(routine) << "parentID " << parentID << "  " << name << "  " << pdg
                   << "vertex: energy " << T0 << "  pos " << pos0 << "  dir " <<
    dir0
                   << "novertex: energy " << T1 << "  pos " << pos1 << "  dir "
    << dir1
                   << endlog;
     */
  } else if (parentID > 0) {
    // first daughters
    G4String pname = aTrack->GetCreatorProcess()->GetProcessName();
    if (pname == "RadioactiveDecay") {
      ;
      //      pdg = pDef->GetPDGEncoding();
      //      energy = aTrack->GetDynamicParticle()->GetKineticEnergy();
    } else {
      // Se e' un nucleo secondario, nel rivelatore, NON generato da Radioactive
      // incrementa il numero di interazioni
      G4String tipo = pDef->GetParticleType();
      if (tipo == "nucleus") {
        G4int Z = ((G4Ions*)pDef)->GetAtomicNumber();
        if (Z > 3) {
          // G4bool checkVolume = aTrack->GetVolume()->GetName() == "ActiveLAr";
          //           const G4UserEventAction* myEventAction =
          //           G4RunManager::GetRunManager()->GetUserEventAction();
          //           DSEventAction* myEventAction =
          //           G4RunManager::GetRunManager()->GetUserEventAction();
          /*
          if (checkVolume) {
            //This is the theta wrt. the reference frame (= E field)
            G4ThreeVector mom = aTrack->GetMomentumDirection();
            //This is with respect to the primary neutron
            const Primary* thePrimary = static_cast<const
          Primary*>(G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
            G4ThreeVector neutronDir = thePrimary->FindNeutronDirection();
            //G4cout << "mom: " << mom << " neutronDir: " << neutronDir <<
          G4endl; G4double cosTh = mom.dot(neutronDir);
            //G4cout << cosTh << " " << std::acos(cosTh)/deg << " " <<
            //	(mom-neutronDir).theta()/deg << G4endl;
            fEventAction->AddRecoil(mom.theta(),mom.phi(),
          aTrack->GetKineticEnergy(), std::acos(cosTh));
          }
          */
        }
      }
    }
  }
  return fUrgent;
}

void DSStackingReD::DSNewStage() {}

void DSStackingReD::DSPrepareNewEvent() {}
