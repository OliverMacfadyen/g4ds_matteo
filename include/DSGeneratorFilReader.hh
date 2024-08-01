// --------------------------------------------------------------------------//
/**
 * AUTHOR: A. Meregaglia
 */
// --------------------------------------------------------------------------//

#ifndef _DSGENERATORFilReader_HH
#define _DSGENERATORFilReader_HH

//---------------------------------------------------------------------------//

#include "DSVGenerator.hh"
#include "G4Event.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4SPSAngDistribution.hh"
#include "G4SPSPosDistribution.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"

using namespace std;

//---------------------------------------------------------------------------//
class DSGeneratorFilReaderMessenger;
class DSGeneratorFilReader : public DSVGenerator {
 public:
  // default constructor
  DSGeneratorFilReader();

  // destructor
  ~DSGeneratorFilReader();

  // public interface
  virtual void DSGeneratePrimaries(G4Event* event);

  void SetFile(string val) { fFile = val; }

  void SetNumberOfSkippedEvents(int val) { fSkipEvents = val; }


 private:
  DSGeneratorFilReaderMessenger* fMessenger;
  G4SPSAngDistribution* fSPSAng;
  G4SPSPosDistribution* fSPSPos;

  G4ParticleTable* fParticleTable;
  G4ParticleDefinition* fParticle;
  G4ThreeVector fPosition;
  G4ThreeVector fDirection;

  G4bool fRead;
  G4int fSkipEvents;
  G4int fNSim;
  G4String fFile;
 
  G4int fnumDau; 
  G4bool fisDau;
  G4bool fisRead;
  G4int fnDau;

  G4int fMultiplicator ; 
  G4int fSimulationIndex ; 
};
#endif
/*
 * Revision 1.1  2013/03/20 09:54:28  dfranco
 * New version of g4ds
 *
 *
 */
