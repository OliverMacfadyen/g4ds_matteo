#ifndef DSLight3_h
#define DSLight3_h 1

#include <vector>
#include "G4DynamicParticle.hh"
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4OpticalPhoton.hh"
#include "G4ParticleMomentum.hh"
#include "G4PhysicsOrderedFreeVector.hh"
#include "G4Poisson.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4VRestDiscreteProcess.hh"
#include "Randomize.hh"
#include "globals.hh"
#include "templates.hh"
#include "G4SystemOfUnits.hh"

using namespace std;

class DSLight3 : public G4VRestDiscreteProcess {
 private:
 public:
  DSLight3(const G4String& processName = "DSLight3", G4ProcessType type = fElectromagnetic);
  ~DSLight3();

 public:
  G4bool IsApplicable(const G4ParticleDefinition& aParticleType);
  G4double GetMeanFreePath(const G4Track& aTrack, G4double, G4ForceCondition*);
  G4double GetMeanLifeTime(const G4Track& aTrack, G4ForceCondition*);
  G4VParticleChange* PostStepDoIt(const G4Track& aTrack, const G4Step& aStep);
  G4VParticleChange* AtRestDoIt(const G4Track& aTrack, const G4Step& aStep);
  float SaggedGasPocket(G4ThreeVector);
  float ExtractionField(float);
  float RCorrectionFactor(G4ThreeVector);
  void CreateS2(G4ThreeVector, G4double, G4ParticleChange*);
  void CreateClusters( G4double, G4double, G4double, G4ThreeVector, G4double, G4double); //, G4int);
  G4double GetLArDriftVelocity(G4double, G4double);
  G4float  GetRecombinationWrt200Vcm(G4float field);

  

  void SetTrackSecondariesFirst(G4bool val) { fTrackSecondariesFirst = val; }
  G4bool GetTrackSecondariesFirst() { return fTrackSecondariesFirst; }

  /*
  void SetScintillationYieldFactor(G4double val)            { fYieldFactorS1 =
  val; } G4double GetScintillationYieldFactor()                    { return
  fYieldFactorS1; }

  void SetScintillationYieldFactorS2(G4double val)          { fYieldFactorS2 =
  val; } G4double GetScintillationYieldFactorS2()                  { return
  fYieldFactorS2; }

  void SetScintillationExcitationRatio(G4double val)        { fExcitationRatio =
  val; } G4double GetScintillationExcitationRatio()                { return
  fExcitationRatio; }
*/

  double GetLArNuclearQuenching(G4double);
  static double GetRandomGamma(const double a, const double b);
  static unsigned int GetRandomNegativeBinomial(double mu, double sig);
  G4double GetPromptFraction(G4int, G4double);

 protected:
  int fTracker[10000];

 private:
  void ReadData();
  G4int BinomFluct(G4int, G4double);

  G4double Interpolator(double, vector<float>&, vector<float>&);

  G4double GetRecoProbAt200Vcm(double, bool);

  G4bool fTrackSecondariesFirst;

  // G4bool isFirstEnDepInLAr;
  G4bool IsARISQuenching;
  G4bool IsMeiQuenching;
  G4bool IsLindhardQuenching;

  G4double fRecoProb;
  G4double fQuenchingFactor;
  G4double fLightYield;
  G4double fS2Yield;
  G4double fS1ScaleFactor;
  G4double fMeanQuantaEnergy;
  G4double fLindhardFactor;

  G4double fExcitationRatioER;
  G4double fExcitationRatioNR;
  G4double fPhotonEneMean;
  G4double fPhotonEneWidth;

  G4double fLArGArBoundaryZ;
  //G4double fDriftFieldAboveGrid;
  G4double fGridSurfaceDistance;

  G4double fD_T;
  G4double fD_L;
  G4double fLArTauFast;
  G4double fLArTauSlow;
  G4double fTauReco;
  G4double fGArTauFast;
  G4double fGArTauSlow;

  G4double fTaueLAr;
  G4double fGArSingletToTriplet;
  G4double fSigma0_Sq;

  G4double fGASGAP;
  G4double fYieldScalingFactor;

  vector<float> fArgonEnergy;
  vector<float> fArgonELoss;

  vector<float> fElectronEnergy;
  vector<float> fElectronELoss;

  vector<float> fAlphaEnergy;
  vector<float> fAlphaELoss;

  vector<float> fRecombinationEnergy;
  vector<float> fRecombinationProbability;

  vector<float> fS2CorrectionPosition;
  vector<float> fS2CorrectionEfficiency;

  vector<float> recombination_ARIS_200_energies;
  vector<float> recombination_ARIS_200_values;

  vector<float> LEff_ARIS_energies;
  vector<float> LEff_ARIS_values;

  G4double fClustering_deltaR_max_cm = 20    ;
  G4double fClustering_deltaT_max_ns = 2000  ;  // ns
  G4double fClustering_dist_max_z_cm = 0.2   ;
  // char fLArPropertiesFile;
};

#endif

/*
 * $Log: DSLight3.hh,v $
 * Revision 1.9  2016/03/03 23:02:43  cz2
 * Added functions related to S2 radial correction.
 *
 * Revision 1.8  2016/02/26 10:48:43  pagnes
 * new calibration, used for MC paper plots (S1 ER from 37Ar,83mKr and 39Ar).
 * Mei NR quenching added
 *
 * Revision 1.7  2016/02/23 14:37:34  riffard
 * Modified S2 with a negative binomial
 *
 * Revision 1.6  2016/02/08 13:30:56  dfranco
 * added S2 radial correction to DS50
 *
 * Revision 1.5  2015/06/02 16:25:55  pagnes
 * added command to use a different LArScintillationProperties.txt file
 * (/ds/manager/lar_properties)
 *
 * Revision 1.4  2015/01/20 09:22:05  dfranco
 * do not produce light for neutral particles in DSLight3 and new model for f90
 *
 * Revision 1.3  2014/12/22 14:40:49  dfranco
 * added the option to activate the recombination probability at 200 V/cm
 * (/ds/physics/tunedS1); this option is by default true; selecting a specific
 * drift field automatically switch off the tunedS1 option
 *
 * Revision 1.2  2014/07/23 14:55:39  pagnes
 * S2 first tuning
 *
 * Revision 1.1  2014/06/27 14:54:32  perassos
 * Implementation of a new scintillation model
 *
 * Revision 1.1  2014/05/07 12:20:53  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.7  2014/01/07 14:14:19  perassos
 * Bug fixed in the storage of the energy deposits and major changes in the
 * application of the TI and DB models
 *
 * Revision 1.6  2013/08/16 15:33:51  perassos
 * Added S1 and S2 to DSLight3
 *
 * Revision 1.5  2013/08/06 13:58:18  dfranco
 * Added 3 variables to store the deposited energy in LAr, scintillator, and
 * water. The last two are not yet implemented. g4rooter has been updated with 3
 * new variables: tpcene, vetoene, and muene
 *
 * Revision 1.4  2013/08/02 15:46:00  dfranco
 * Further development on DSLight3
 *
 * Revision 1.3  2013/08/02 12:25:24  dfranco
 * Development of new DSLight class
 *
 * Revision 1.2  2013/08/01 14:28:09  dfranco
 * added energy loss data from SRIM/TRIM
 *
 * Revision 1.1  2013/07/25 09:55:59  dfranco
 * Added a second version of NEST, still in development, and actually not
 * working. The default version is DSLight
 *
 *
 */
