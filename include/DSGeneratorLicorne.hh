//---------------------------------------------------------------------------//

#ifndef _DSGENERATORLICORNE_HH
#define _DSGENERATORLICORNE_HH

//---------------------------------------------------------------------------//

#include <G4Event.hh>
#include <G4ParticleDefinition.hh>
#include <G4ParticleTable.hh>
#include <G4ThreeVector.hh>

#include "DSGeneratorLicorneMessenger.hh"
#include "DSVGenerator.hh"
#include "G4SPSAngDistribution.hh"
#include "G4SPSPosDistribution.hh"
#include "G4SPSRandomGenerator.hh"
//---------------------------------------------------------------------------//

class DSGeneratorLicorneMessenger;

class DSGeneratorLicorne : public DSVGenerator {

 public:
  DSGeneratorLicorne();
  ~DSGeneratorLicorne();

  virtual void DSGeneratePrimaries(G4Event*);

  void SetIsotope(string val) { fIsotope = val; }

  void SetGammaToNeutronRatio(double val) { fGammaToNeutronRatio = val; }
  void SetNeutronRate(double val) { fNeutronRate = val; }
  void SetRunTime(double val) { fRunTime = val; }
  void SetPulsePeriod(double val) { fPulsePeriod = val; }
  void SetPulseWidth(double val) { fPulseWidth = val; }
  void SetPulseMode(bool val) { fPulseMode = val; }

 private:
  void LoadNeutronEnergy();
  void LoadNeutronSpectra();
  void LoadNeutronAngle();
  void PickEneAngle(G4double&, G4double&);
  G4double PickAnEnergy();
  G4double PickAnAngle();

  DSGeneratorLicorneMessenger* fMessenger;

  string fIsotope;
  bool fIsFirstEvent;
  G4ParticleDefinition* fParticle;
  G4ParticleTable* fParticleTable;

  G4ThreeVector fDirection;
  G4ThreeVector fPosition;

  double fSpectrumEne[1000];
  double fSpectrumCumul[1000];
  int fSamples;

  double fAngle[1000];
  double fAngleCumul[1000];
  int fSamplesAngle;

  const static int fEnergy_tot_bin = 25;
  const static int fAngle_tot_bin = 700;
  double fLicorne_array[fEnergy_tot_bin][fAngle_tot_bin];
  double fMax_prob;
  double fEnergy_width;
  double fAngle_width;

  double fGammaToNeutronRatio;
  double fNeutronRate;
  double fRunTime;
  double fPulsePeriod;
  double fPulseWidth;
  bool fPulseMode;

  G4SPSAngDistribution* fAngDist;
  G4SPSPosDistribution* fPosDist;
  G4SPSRandomGenerator* fRandom;
};

#endif
