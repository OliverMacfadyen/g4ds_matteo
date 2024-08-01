//---------------------------------------------------------------------------//

#ifndef _DSGENERATORMULTISPECTRA_HH
#define _DSGENERATORMULTISPECTRA_HH

//---------------------------------------------------------------------------//

#include <G4Event.hh>
#include <G4ParticleDefinition.hh>
#include <G4ParticleTable.hh>
#include <G4ThreeVector.hh>

#include "DSGeneratorMultiSpectraMessenger.hh"
#include "DSVGenerator.hh"

//---------------------------------------------------------------------------//

class DSGeneratorMultiSpectraMessenger;

class DSGeneratorMultiSpectra : public DSVGenerator {
 public:
  DSGeneratorMultiSpectra();
  ~DSGeneratorMultiSpectra();

  virtual void DSGeneratePrimaries(G4Event*);

  void SetNSpectra(G4int val) {
    fNSpectra = val;
    initArrays();
  }
  //   void                                       SetIsotope(int s, G4String
  //   val)      { if (s < fNSpectra && fIsotope)   fIsotope[s]   = val; } void
  //   SetMinEnergy(int s, double val)      { if (s < fNSpectra && fMinEnergy)
  //   fMinEnergy[s] = val; } void SetMaxEnergy(int s, double val)      { if (s
  //   < fNSpectra && fMaxEnergy) fMaxEnergy[s] = val; }

 private:
  void initArrays();
  G4int loadSpectrum(int);
  G4double pickAnEnergy(int);

  DSGeneratorMultiSpectraMessenger* fMessenger;

  G4int fNSpectra;
  G4String* fIsotope;
  G4String* fUnits;
  G4double* fMinEnergy;
  G4double* fMaxEnergy;
  // vector<G4double>                           *fSpectrumEne;
  vector<G4double>* fSpectrumCumul;

  G4ParticleDefinition* fParticle;
  G4ParticleTable* fParticleTable;

  G4ThreeVector fDirection;
  G4ThreeVector fPosition;
};

#endif
