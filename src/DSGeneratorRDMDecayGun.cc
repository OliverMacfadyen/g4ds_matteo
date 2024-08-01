// --------------------------------------------------------------------------//
/**
 * AUTHOR: D. Franco
 * CONTACT: davide.franco@mi.infn.it
 */
// --------------------------------------------------------------------------//

#include "DSGeneratorRDMDecayGun.hh"
#include "DSEventHandler.hh"
#include "DSGeneratorRDMDecayGunMessenger.hh"
#include "DSLogger.hh"
#include "G4Event.hh"
#include "G4IonTable.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4RadioactiveDecay.hh"
#include "G4SPSRandomGenerator.hh"
#include "G4SystemOfUnits.hh"
#include "G4UImanager.hh"
//---------------------------------------------------------------------------//

DSGeneratorRDMDecayGun::DSGeneratorRDMDecayGun() : DSVGenerator("DSGeneratorRDMDecayGun") {
  fParticleGun = new G4ParticleGun;
  IsFirst = true;
  SetPositionFlag(1);
  fTheMessenger = new DSGeneratorRDMDecayGunMessenger(this);
  if (!fParticleGun) {
    DSLog(error) << "Could not allocate G4ParticleGun! Out of memory?" << endlog;
    DSLog(fatal) << endlog;
  }
  DSLog(routine) << "G4ParticleGun Constructed." << endlog;
}
//---------------------------------------------------------------------------//

// DSGeneratorRDMDecayGun::DSGeneratorRDMDecayGun(const DSGeneratorRDMDecayGun &
// other)
//{;}

//---------------------------------------------------------------------------//

void DSGeneratorRDMDecayGun::SetNucleus(DSGeneratorRDMNucleus theIon1) {

  fA = theIon1.GetA();
  fZ = theIon1.GetZ();
  fE = theIon1.GetE();
}

DSGeneratorRDMDecayGun::~DSGeneratorRDMDecayGun() {
  delete fTheMessenger;
  delete fParticleGun;
}

//---------------------------------------------------------------------------//

void DSGeneratorRDMDecayGun::DSGeneratePrimaries(G4Event* event) {

  if (IsFirst) {
    aIon = G4IonTable::GetIonTable()->GetIon(fZ, fA, fE);
    DSLog(routine) << "**************************" << endlog;
    DSLog(routine) << "  Isotope: " << aIon->GetParticleName() << endlog;
    DSLog(routine) << "**************************" << endlog;
    IsFirst = false;
  }

  fParticleGun->SetParticleDefinition(aIon);
  fParticleGun->SetParticlePosition(GetVParticlePosition());
  fParticleGun->SetParticleMomentumDirection(GetVParticleDirection());
  fParticleGun->SetParticleEnergy(GetVParticleEnergy(fParticleGun->GetParticleDefinition()));

  DSEventHandler::Get()->SetPDG(fParticleGun->GetParticleDefinition()->GetPDGEncoding());
  DSEventHandler::Get()->SetPosition(fParticleGun->GetParticlePosition() / cm);
  DSEventHandler::Get()->SetDirection(fParticleGun->GetParticleMomentumDirection());
  DSEventHandler::Get()->SetEnergy(fParticleGun->GetParticleEnergy() / keV);

  fParticleGun->GeneratePrimaryVertex(event);

  DSLog(development) << fParticleGun->GetParticleDefinition()->GetParticleName() << " "
                     << "Energy: " << fParticleGun->GetParticleEnergy() / keV << " keV; "
                     << "Position: " << fParticleGun->GetParticlePosition() / cm << " cm;  "
                     << "Direction: " << fParticleGun->GetParticleMomentumDirection() << " " << endlog;
}

/*
 * Revision 1.1  2013/03/22 13:50:09  dfranco
 * Added radioactive decay and radioactive decay chain generators. They require
 * the correspondent stacking actions. Two mac files are included as examples
 *
 */
