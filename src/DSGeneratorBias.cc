// --------------------------------------------------------------------------//
/**
 * AUTHOR: D. Franco (davide.franco@apc.in2p3.fr)
 * 
 * DATE: 16/11/2023
 */
// --------------------------------------------------------------------------//

#include "DSGeneratorBias.hh"
#include "DSEventHandler.hh"
#include "DSG4DSReader.hh"
#include "DSGeneratorBiasMessenger.hh"
#include "DSIO.hh"
#include "DSLogger.hh"
#include "DSStorage.hh"
#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SPSRandomGenerator.hh"
#include "G4SystemOfUnits.hh"
#include "G4IonTable.hh"
#include "CLHEP/Random/RandGauss.h"
#include "DSManager.hh"

//---------------------------------------------------------------------------//

DSGeneratorBias::DSGeneratorBias() : DSVGenerator("DSGeneratorBias") {

  fSkipEvents = 0;
  fAmplificationFactor = 1;
  fDirectionSmearing = 1; // degrees
  fEnergySmearing = 0.01; // fraction
  fMinimumEnergy = 0;
  fCounter = 0;

  fParticleTable = G4ParticleTable::GetParticleTable();

  DSStorage::Get()->SetSkipEventsWithNoDaughters(true);

  fMessenger = new DSGeneratorBiasMessenger(this);

  DSLog(routine) << "Bias Generatore Built" << endlog;
}
//---------------------------------------------------------------------------//

DSGeneratorBias::~DSGeneratorBias() {
  delete fMessenger;
}

//---------------------------------------------------------------------------//

void DSGeneratorBias::DSGeneratePrimaries(G4Event* event) {
  G4bool isRead     = false;
  G4bool isDaughter = false;

  if (fSkipEvents > 0) {
    DSG4DSReader::Get()->SkipEvents(fSkipEvents);
    fSkipEvents = 0;
  }

  
  if (fCounter == 0) { 
    while (!isDaughter) {
      DSG4DSReader::Get()->ClearAll();
      isRead = DSG4DSReader::Get()->ReadEvent();
      

      if (!isRead) {
        DSIO::Get()->CloseG4DSFile();
        DSStorage::Get()->SetAbortRun(true);
        fParticle = fParticleTable->FindParticle(0);
        G4PrimaryParticle* particle = new G4PrimaryParticle(fParticle);
        particle->SetKineticEnergy(0.1*eV);
        particle->SetMomentumDirection(G4ThreeVector(0,0,1));

        G4PrimaryVertex* vertex = new G4PrimaryVertex(G4ThreeVector(0,0,0), 0 * ns);

        vertex->SetPrimary(particle);
        event->AddPrimaryVertex(vertex);

        return ;
      }

      int ndau = int(DSG4DSReader::Get()->GetVDaughters().size());
      if (ndau > 0) {
        for(int i=0;i<ndau;++i) {
          if(DSG4DSReader::Get()->GetVDaughters()[i].Energy*keV > fMinimumEnergy) isDaughter = true;
        }
      }
    }
  }

  DSEventHandler::Get()->SetPDG(DSG4DSReader::Get()->GetEvent().PDG);
  DSEventHandler::Get()->SetPosition(G4ThreeVector(DSG4DSReader::Get()->GetEvent().Position[0],DSG4DSReader::Get()->GetEvent().Position[1],DSG4DSReader::Get()->GetEvent().Position[2]));
  DSEventHandler::Get()->SetDirection(G4ThreeVector(DSG4DSReader::Get()->GetEvent().Direction[0],DSG4DSReader::Get()->GetEvent().Direction[1],DSG4DSReader::Get()->GetEvent().Direction[2]));
  DSEventHandler::Get()->SetEnergy(DSG4DSReader::Get()->GetEvent().Energy);
  DSEventHandler::Get()->SetTime(DSG4DSReader::Get()->GetEvent().Time);



  for (int i = 0; i < G4int(DSG4DSReader::Get()->GetVDaughters().size()); ++i) {    
    int PDG = DSG4DSReader::Get()->GetVDaughters()[i].PDG ;
    if (PDG < 10000)
      fParticle = fParticleTable->FindParticle(PDG);
    else 
      fParticle = G4IonTable::GetIonTable()->GetIon(PDG);
    
    float energy = DSG4DSReader::Get()->GetVDaughters()[i].Energy*keV;

    if (fEnergySmearing > 0) energy =  CLHEP::RandGauss::shoot(energy, fEnergySmearing*energy);

    G4ThreeVector dir = G4ThreeVector(DSG4DSReader::Get()->GetVDaughters()[i].Direction[0], 
                                      DSG4DSReader::Get()->GetVDaughters()[i].Direction[1],
                                      DSG4DSReader::Get()->GetVDaughters()[i].Direction[2]);

    G4ThreeVector pos = G4ThreeVector(DSG4DSReader::Get()->GetVDaughters()[i].Position[0], 
                                      DSG4DSReader::Get()->GetVDaughters()[i].Position[1],
                                      DSG4DSReader::Get()->GetVDaughters()[i].Position[2]);

    
    if (fDirectionSmearing > 0) {
      double theta = CLHEP::RandGauss::shoot(pos.theta(),fDirectionSmearing*M_PI/180.); 
      double phi   = CLHEP::RandGauss::shoot(pos.phi(),fDirectionSmearing*M_PI/180.); 
      pos.setTheta(theta);
      pos.setPhi(phi);
    }

    G4PrimaryParticle* particle = new G4PrimaryParticle(fParticle);
    particle->SetKineticEnergy(DSG4DSReader::Get()->GetVDaughters()[i].Energy*keV);
    particle->SetMomentumDirection(dir/dir.mag());

    G4PrimaryVertex* vertex = new G4PrimaryVertex(pos*cm, DSG4DSReader::Get()->GetVDaughters()[i].Time * ns);

    vertex->SetPrimary(particle);
    event->AddPrimaryVertex(vertex);

    DSLog(development) << "bias " << i << " " 
         << DSG4DSReader::Get()->GetVDaughters()[i].PDG << " "
         << energy/keV << " "
         << pos/cm << endlog;


  }

  
  fCounter++;
  if (fCounter == fAmplificationFactor) {
    fCounter = 0;
    DSG4DSReader::Get()->ClearAll();
  }
}

