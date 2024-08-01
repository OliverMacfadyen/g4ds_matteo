//---------------------------------------------------------------------------//

#ifndef _DSGENERATORKR85m_HH
#define _DSGENERATORKR85m_HH

//---------------------------------------------------------------------------//

#include <G4Event.hh>
#include <G4ParticleDefinition.hh>
#include <G4ParticleTable.hh>
#include <G4ThreeVector.hh>

#include "DSVGenerator.hh"

//---------------------------------------------------------------------------//

class DSGeneratorKr85mMessenger;

class DSGeneratorKr85m : public DSVGenerator {

 public:
  DSGeneratorKr85m();
  ~DSGeneratorKr85m();

  virtual void DSGeneratePrimaries(G4Event*);

 private:
  void LoadAr39CrossSection();
  G4double PickAnEnergy();

  G4bool fIsFirstEvent;

  G4ParticleDefinition* fParticle1;
  G4ParticleDefinition* fParticle2;
  G4ParticleTable* fParticleTable;

  G4double fSpectrumEne[1000];
  G4double fSpectrumCumul[1000];
  G4int fSamples;
};

#endif
