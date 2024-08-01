// --------------------------------------------------------------------------//
/**
 * AUTHOR: A. Meregaglia
 */
// --------------------------------------------------------------------------//

#ifndef _DSGENERATORBIAS_HH
#define _DSGENERATORBIAS_HH

//---------------------------------------------------------------------------//

#include "DSVGenerator.hh"
#include "G4Event.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"

using namespace std;

//---------------------------------------------------------------------------//
class DSGeneratorBiasMessenger;
class DSGeneratorBias : public DSVGenerator {
 public:
  // default constructor
  DSGeneratorBias();

  // destructor
  ~DSGeneratorBias();

  // public interface
  virtual void DSGeneratePrimaries(G4Event* event);

  //void SetFile(string val) { fFile = val; }

  void SetNumberOfSkippedEvents(int val) { fSkipEvents = val; }
  void SetAmplificationFactor(int val)   { fAmplificationFactor = val; }
  void SetDirectionSmearing(float val)   { fDirectionSmearing = val; }
  void SetEnergySmearing(float val)      { fEnergySmearing = val; }
  void SetMinimumEnergy(float val)       { fMinimumEnergy = val; }


 private:
  DSGeneratorBiasMessenger* fMessenger;
  G4ParticleDefinition* fParticle; 
  G4ParticleTable* fParticleTable;
  
  G4int   fSkipEvents;
  G4int   fCounter;
  G4int   fAmplificationFactor;
  G4float fDirectionSmearing;
  G4float fEnergySmearing;
  G4float fMinimumEnergy;
};
#endif
