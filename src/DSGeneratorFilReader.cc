// --------------------------------------------------------------------------//
/**
 * AUTHOR: A. Meregaglia
 */
// --------------------------------------------------------------------------//

#include "DSGeneratorFilReader.hh"
#include "DSEventHandler.hh"
#include "DSG4DSReader.hh"
#include "DSGeneratorFilReaderMessenger.hh"
#include "DSIO.hh"
#include "DSLogger.hh"
#include "DSStorage.hh"
#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SPSRandomGenerator.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
//---------------------------------------------------------------------------//

DSGeneratorFilReader::DSGeneratorFilReader() : DSVGenerator("DSGeneratorFilReader") {

  fisRead = false;
  fisDau = false;
  fnDau = 0;
  fnumDau = 0;
  fRead = false;
  fSkipEvents = 0 ;
  fNSim = 0; 
  DSStorage::Get()->SetOverwriteCounter(false);
  G4SPSRandomGenerator* RndGen = new G4SPSRandomGenerator;
  fSPSPos = new G4SPSPosDistribution;
  fSPSAng = new G4SPSAngDistribution;
  fSPSAng->SetBiasRndm(RndGen);
  fSPSPos->SetBiasRndm(RndGen);
  fSPSAng->SetAngDistType("iso");
  fSPSAng->SetPosDistribution(fSPSPos);
  
  fSimulationIndex = 0 ; 
  fMultiplicator = DSStorage::Get()->GetMultiplicator() ;
  fParticleTable = G4ParticleTable::GetParticleTable();

  fMessenger = new DSGeneratorFilReaderMessenger(this);

  DSLog(routine) << "FilReader Generator Built" << endlog;
}
//---------------------------------------------------------------------------//

DSGeneratorFilReader::~DSGeneratorFilReader() {
  delete fMessenger;
}

//---------------------------------------------------------------------------//

void DSGeneratorFilReader::DSGeneratePrimaries(G4Event* event) {

 
/*
  if (fSkipEvents > 0) {
    DSG4DSReader::Get()->SkipEvents(fSkipEvents);
    fSkipEvents = 0;
  }
*/
    while (!fisDau) {
      fisRead = DSG4DSReader::Get()->ReadEvent();
      //cout << "fisRead " << fisRead << endl; // DO NOT REMOVE

      if (fisRead && fNSim<10000) {
        fnDau = int(DSG4DSReader::Get()->GetVDaughters().size());
        if (fnDau > 0) fisDau = true;
        fSimulationIndex = 0 ;
        fNSim++ ; 
      }
      else {
        DSLog(error) << " End of file reached. Next particle is ghost." << endlog;
        DSLog(routine) << " End of file reached. Next particle is ghost." << endlog;
        cout  << " End of file reached or 1000 ev without progress. Quitting. Next particle is ghost." << endl;
        
        DSIO::Get()->CloseG4DSFile();
        std::exit(0); 
        return;
        //break;
      }
    }

        
        G4int i = fnumDau;

        G4double ene = DSG4DSReader::Get()->GetVDaughters()[i].Energy;
        //cout << "Got daughter: " << i << " " << ene/keV << " " << fnDau  << " " << fMultiplicator  <<  " " << fSimulationIndex << endl;

        if ( fSimulationIndex < fMultiplicator  && ene > 300*keV ) {

          G4double xx = DSG4DSReader::Get()->GetVDaughters()[i].Position[0];
          G4double yy = DSG4DSReader::Get()->GetVDaughters()[i].Position[1];
          G4double zz = DSG4DSReader::Get()->GetVDaughters()[i].Position[2];

          G4double px = DSG4DSReader::Get()->GetVDaughters()[i].Direction[0];
          G4double py = DSG4DSReader::Get()->GetVDaughters()[i].Direction[1];
          G4double pz = DSG4DSReader::Get()->GetVDaughters()[i].Direction[2];
          G4double dir_len = sqrt(px*px + py*py + pz*pz);
          
          px = px/dir_len;
          py = py/dir_len;
          pz = pz/dir_len;

          G4double ene0 = DSG4DSReader::Get()->GetPrimaryEventEnergy();
          G4double weight = DSG4DSReader::Get()->GetVDaughters()[i].Time; //For biasing purpose, time is particle weight

          fPosition = G4ThreeVector(xx*cm, yy*cm, zz*cm);
          fDirection = G4ThreeVector(px, py, pz);
          fParticle = fParticleTable->FindParticle(DSG4DSReader::Get()->GetVDaughters()[i].PDG);

          //G4double phi  = G4RandUniform() * TWO_PI ; 
          G4double angle = G4RandGauss::shoot(0.0, 2.) *degree ; 
        

	  //G4double angle = 3.0 * degree;

	  G4ThreeVector randomAxis(G4UniformRand(), G4UniformRand(), G4UniformRand());
          randomAxis = randomAxis.unit();

	  G4RotationMatrix rotation;
	  rotation.rotate(angle, randomAxis);

	  fDirection = rotation * fDirection;


          G4PrimaryParticle* particle = new G4PrimaryParticle(fParticle);
          particle->SetKineticEnergy(ene);
          particle->SetMomentumDirection(fDirection);

          G4PrimaryVertex* vertex = new G4PrimaryVertex(fPosition, 0);

          vertex->SetWeight(weight);
          vertex->SetPrimary(particle);
          event->AddPrimaryVertex(vertex);

          DSEventHandler::Get()->SetPDG(DSG4DSReader::Get()->GetVDaughters()[i].PDG);
          DSEventHandler::Get()->SetPosition( DSG4DSReader::Get()->GetPrimaryPosition() );
          DSEventHandler::Get()->SetDirection(fDirection);
          DSEventHandler::Get()->SetEnergy(ene0);
          
          DSEventHandler::Get()->SetUserDouble(weight); //Weight is stored in UserDouble
          DSEventHandler::Get()->SetUserFloat1( ene / keV );
          
          DSEventHandler::Get()->SetUserFloat2( xx  ) ; 
          DSEventHandler::Get()->SetUserInt1( yy  ) ; 
          DSEventHandler::Get()->SetUserInt2( zz ) ; 
          
          DSEventHandler::Get()->SetUsers();

          DSEventHandler::Get()->SetRate(fMultiplicator);
          //cout << " made vertex " << fNSim <<  endl ;   

          fSimulationIndex++; 

        
          //cout << ene <<  " " << ene0 << " " << weight << endl ; 
        } else {
          
          fSimulationIndex = 0 ;
          fnumDau++;
        
        }


      if (fnumDau == fnDau) {
        fisDau = false; 
        fnumDau = 0;     
        fnDau = 0;
        DSG4DSReader::Get()->ClearAll();
        fNSim =0 ; 
      }
}

/*
 * Revision 1.1  2013/03/20 09:54:28  dfranco
 * New version of g4ds
 *
 *
 *      
 *
 */
