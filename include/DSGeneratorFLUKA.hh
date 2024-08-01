// --------------------------------------------------------------------------//
/**
 * AUTHOR: A. Meregaglia
 */
// --------------------------------------------------------------------------//

#ifndef _DSGENERATORFLUKA_HH
#define _DSGENERATORFLUKA_HH

//---------------------------------------------------------------------------//

#include <fstream>
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
// class DSGeneratorFLUKAMessenger;
class DSGeneratorFLUKA : public DSVGenerator {
 public:
  // default constructor
  DSGeneratorFLUKA();

  // destructor
  ~DSGeneratorFLUKA();

  // public interface
  virtual void DSGeneratePrimaries(G4Event* event);

  void SetFLUKAFileName(G4String val);

  // void SetFile (string val)   { fFile=val ; }

 private:
  int ConvertFlukaPDG(int);
  G4ThreeVector SmearDirection(G4ThreeVector , G4double );


  // DSGeneratorFLUKAMessenger*   fMessenger;
  G4SPSAngDistribution* fSPSAng;
  G4SPSPosDistribution* fSPSPos;

  G4ParticleTable* fParticleTable;
  G4ParticleDefinition* fParticle;
  G4ThreeVector fPosition;
  G4ThreeVector fDirection;
  G4bool fRead;
  G4int fnumev;
  ifstream fFile;
};
#endif
/*
 * Revision 1.1  2013/03/20 09:54:28  dfranco
 * New version of g4ds
 *
 *
 */
