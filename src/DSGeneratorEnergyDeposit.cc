// --------------------------------------------------------------------------//
/**
 * AUTHOR: A. Meregaglia
 */
// --------------------------------------------------------------------------//

#include "DSGeneratorEnergyDeposit.hh"
#include "DSEventHandler.hh"
#include "DSG4DSReader.hh"
#include "DSGeneratorEnergyDepositMessenger.hh"
#include "DSIO.hh"
#include "DSLogger.hh"
#include "DSStorage.hh"
#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SPSRandomGenerator.hh"
#include "G4SystemOfUnits.hh"
#include "G4IonTable.hh"
//---------------------------------------------------------------------------//

DSGeneratorEnergyDeposit::DSGeneratorEnergyDeposit() : DSVGenerator("DSGeneratorEnergyDeposit") {

  fRead = false;
  fnumev = 0;
  fSkipEvents = 0 ;
  DSStorage::Get()->SetOverwriteCounter(false);
  G4SPSRandomGenerator* RndGen = new G4SPSRandomGenerator;
  fSPSPos = new G4SPSPosDistribution;
  fSPSAng = new G4SPSAngDistribution;
  fSPSAng->SetBiasRndm(RndGen);
  fSPSPos->SetBiasRndm(RndGen);
  fSPSAng->SetAngDistType("iso");
  fSPSAng->SetPosDistribution(fSPSPos);

  fParticleTable = G4ParticleTable::GetParticleTable();

  fMessenger = new DSGeneratorEnergyDepositMessenger(this);

  DSStorage::Get()->SetIsEnDepGenerator(true);

  DSLog(routine) << "EnergyDeposit Generatore Built" << endlog;
}
//---------------------------------------------------------------------------//

DSGeneratorEnergyDeposit::~DSGeneratorEnergyDeposit() {
  delete fMessenger;
}

//---------------------------------------------------------------------------//

void DSGeneratorEnergyDeposit::DSGeneratePrimaries(G4Event* event) {
  G4bool isRead = false;
  G4bool isNDepo = false;
//  int nSkip = DSStorage::Get()->GetSkipEvents();

  if (fSkipEvents > 0) {
    DSG4DSReader::Get()->SkipEvents(fSkipEvents);
    fSkipEvents = 0;
  }
  
  while (!isNDepo) {
    isRead = DSG4DSReader::Get()->ReadEvent();
    if (int(DSG4DSReader::Get()->GetVDeposits().size()) > 0) isNDepo = true;
  }
  if (!isRead) {
    DSIO::Get()->CloseG4DSFile();
    DSLog(routine) << "End of file reached" << endlog;
    return;
  }

  for (int i = 0; i < G4int(DSG4DSReader::Get()->GetVDeposits().size()); ++i) {

    G4double xx = DSG4DSReader::Get()->GetVDeposits()[i].Position[0] * cm;
    G4double yy = DSG4DSReader::Get()->GetVDeposits()[i].Position[1] * cm;
    G4double zz = DSG4DSReader::Get()->GetVDeposits()[i].Position[2] * cm;

    fPosition = G4ThreeVector(xx, yy, zz);
    
 
    if (DSG4DSReader::Get()->GetVDeposits()[i].PID < 10000)
      fParticle = fParticleTable->FindParticle(DSG4DSReader::Get()->GetVDeposits()[i].PID);
    else 
      fParticle = G4IonTable::GetIonTable()->GetIon(DSG4DSReader::Get()->GetVDeposits()[i].PID);
    
    G4ParticleMomentum theMomentum = fSPSAng->GenerateOne();

    G4PrimaryParticle* particle = new G4PrimaryParticle(fParticle);
    particle->SetKineticEnergy(DSG4DSReader::Get()->GetVDeposits()[i].Energy);
    particle->SetMomentumDirection(theMomentum);

    G4PrimaryVertex* vertex = new G4PrimaryVertex(fPosition, DSG4DSReader::Get()->GetVDeposits()[i].Time * ns);

    vertex->SetPrimary(particle);
    event->AddPrimaryVertex(vertex);
  }

  DSG4DSReader::Get()->ClearAll();
}

/*
 * Revision 1.1  2013/03/20 09:54:28  dfranco
 * New version of g4ds
 *
 *
 */
