//---------------------------------------------------------------------------//
/**
 *
 * CLASS DECLARATION:  DSVGenerator.hh
 *
 *---------------------------------------------------------------------------//
 *
 * DESCRIPTION:
 *
 */
// Begin description of class here
/**
 *
 * Pure  base class for DS generators.
 *
 */
// End class description
//
/**
 * SPECIAL NOTES:
 *
 */
//
// --------------------------------------------------------------------------//
/**
 * AUTHOR: davide.franco@mi.infn.it
 */
// --------------------------------------------------------------------------//

#ifndef _DSVGENERATOR_HH
#define _DSVGENERATOR_HH

//---------------------------------------------------------------------------//

#include <iostream>
#include <string>
#include <vector>
#include "G4ParticleGun.hh"
#include "G4SPSAngDistribution.hh"
#include "G4SPSEneDistribution.hh"
#include "G4SPSPosDistribution.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

using namespace std;
//---------------------------------------------------------------------------//

class G4Event;
class G4UImessenger;
class G4Run;
class DSVGeneratorMessenger;
class G4SPSRandomGenerator;
//---------------------------------------------------------------------------//

class DSVGenerator {
 public:
  DSVGenerator(const G4String& myname);

  ~DSVGenerator();

  G4String GetGeneratorName() { return fGeneratorName; }
  void SetReportingFrequency(G4int freq) { fReportingFrequency = freq; }

  virtual void DSGeneratePrimaries(G4Event* event) = 0;

  void SetVerticalShift(G4double val) { fVerticalShift = val; }

 protected:
  G4UImessenger* fG4Messenger;  // G4Messenger for setting up generator.
  G4int fReportingFrequency;    // Generate report every fReportingFrequency events;

 public:
  G4ThreeVector GenerateInCryostats();
  G4ThreeVector GenerateInTeflon();
  G4ThreeVector GenerateInFusedSilica();
  G4ThreeVector GenerateInPMTPhotocathode();
  G4ThreeVector GenerateInPMTStem();
  G4ThreeVector GenerateInSiPM();
  G4ThreeVector GenerateInPDM();
  G4ThreeVector GenerateInSiPM_Veto();
  G4ThreeVector GenerateInGrid();
  G4ThreeVector GenerateInGridSupport();
  G4ThreeVector GenerateInCopperRings();
  G4ThreeVector GenerateInLiquidArgon();
  G4ThreeVector GenerateInLateralTPB();
  G4ThreeVector GenerateInLateralESR();
  G4ThreeVector GenerateInLateralAcrylicPanels();
  G4ThreeVector GenerateInGasPocket();
  G4ThreeVector GenerateInHolderSource();
  G4ThreeVector GenerateInLiquidScintillator();
  G4ThreeVector GenerateInVetoPMTs();
  G4ThreeVector GenerateInVetoSteelSphere();
  G4ThreeVector GenerateInDS20kSupportStructure();
  G4ThreeVector GenerateInPlasticScintillator();
  G4ThreeVector GenerateInArgonBufferInside();
  G4ThreeVector GenerateInArgonBufferOutside();
  G4ThreeVector GenerateIn7mSphere();
  
  G4ThreeVector GenerateInTPCBarrel();
  G4ThreeVector GenerateInOPsAcrylic();

  G4ThreeVector GenerateInArDMTestArgonCopper();
  G4ThreeVector GenerateInArDMTestArgonNS();
  G4ThreeVector GenerateInArDMTestArgonAcrylic();
  G4ThreeVector GenerateInArDMTestArgonUnderground();
  G4ThreeVector GenerateInArDMTestArgonReflector();
  G4ThreeVector GenerateInArDMTestArgonWLS();
  G4ThreeVector GenerateInArDMDartDet();
  G4ThreeVector GenerateInArDMDartExt();
  G4ThreeVector GenerateInArDMDartSiPM();
  G4ThreeVector GenerateProtoDUNEMembrane();
  G4ThreeVector GenerateProtoDUNEFoam();

  G4ThreeVector GenerateInGantryPipe();
  G4ThreeVector GenerateInTitaniumVessel();
  G4ThreeVector GenerateInAriserSourceHolder();

  G4ThreeVector GenerateInPhysVolume();
  G4ThreeVector GenerateInAnode();
  G4ThreeVector GenerateInCathode();

  void SphereToPrism(G4double, G4double, G4double, G4double, G4double, G4double&, G4double&, G4double&);

  G4bool CheckMaterial(G4ThreeVector, G4String);
  G4bool CheckMaterialIndex(G4ThreeVector, G4int);
  G4bool CheckVolumeName(G4ThreeVector, G4String);

  void SetVParticleDefinition(G4ParticleDefinition* part) { fParticle = part; }
  G4ParticleDefinition* GetVParticleDefinition() { return fParticle; }

  void SetVParticlePosition(G4ThreeVector vect) { fPosition = vect; }
  G4ThreeVector GetVParticlePosition();

  void SetVParticleDirection(G4ThreeVector dir) { fDirection = dir; }
  G4ThreeVector GetVParticleDirection();

  void SetVParticleEnergy(G4double ene) { fEnergy = ene; }
  G4double GetVParticleEnergy(G4ParticleDefinition* pDef);

  void SetIsVolumeDistribution(G4bool);
  G4bool GetIsVolumeDistribution() { return fVolumeFlag; }

  void SetIsEnergyDistribution(G4bool);
  G4bool GetIsEnergyDistribution() { return fEnergyFlag; }

  void SetVParticleCharge(G4double aCharge) { fCharge = aCharge; }
  G4double GetVParticleCharge() { return fCharge; }

  void SetVNumberOfParticles(G4int val) { fNumberOfParticles = val; }
  int GetVNumberOfParticles() { return fNumberOfParticles; }

  // Allows user to choose Point, Plane, Surface or Volume source
  // position distributions.
  void SetPosDisType(G4String string) { fSPSPos->SetPosDisType(string); }

  // Allows the user to choose the particular shape they wish for the
  // osition distribution. Choices are Square, Circle, Ellipse, Rectangle,
  // Sphere, Ellipsoid, Cylinder, Parallelepiped.
  void SetPosDisShape(G4String string) { fSPSPos->SetPosDisShape(string); }

  // Sets the co-ordinates of the centre of the position distribution.
  void SetCentreCoords(G4ThreeVector pos) { fSPSPos->SetCentreCoords(pos); }

  // Used to specify the co-ordinate system for the position distribution
  // along with SetPosRot2. SetPosRot1 sets the vector x' and need not be
  // a unit vector.
  void SetPosRot1(G4ThreeVector rot) { fSPSPos->SetPosRot1(rot); }

  // Used in connection with SetPosRot1. This sets a vector in the plane
  // x'y'. By a series of cross products x', y', z' are generated. Again
  // need not be a unit vector.
  void SetPosRot2(G4ThreeVector rot) { fSPSPos->SetPosRot2(rot); }

  // Sets the half length in x.
  void SetHalfX(G4double val) { fSPSPos->SetHalfX(val); }

  // Sets the half length in y.
  void SetHalfY(G4double val) { fSPSPos->SetHalfY(val); }

  // Sets the half length in z.
  void SetHalfZ(G4double val) { fSPSPos->SetHalfZ(val); }

  // Sets the radius where appropriate for source distribution shapes.
  void SetRadius(G4double radd) { fSPSPos->SetRadius(radd); }

  // Sets the inner radius where appropriate for source distribution shapes.
  void SetRadius0(G4double radd) { fSPSPos->SetRadius0(radd); }

  // Used to confine the start positions to a particular volume.
  void ConfineSourceToVolume(G4String vol) { fSPSPos->ConfineSourceToVolume(vol); }

  void SetEnergyDistribution(G4bool val) { fEnergyFlag = val; }

  void SetEnergyDisType(G4String fstring) { fSPSEne->SetEnergyDisType(fstring); }

  void SetEmin(G4double val) { fSPSEne->SetEmin(val); }

  void SetEmax(G4double val) { fSPSEne->SetEmax(val); }

  // Sets alpha for a power-law distribution.
  void SetAlpha(G4double val) { fSPSEne->SetAlpha(val); }

  // Sets Temperature for a Brem or BBody distributions.
  void SetTemp(G4double val) { fSPSEne->SetTemp(val); }

  // Sets Ezero for an exponential distribution.
  void SetEzero(G4double val) { fSPSEne->SetEzero(val); }

  // Sets gradient for a linear distribution.
  void SetGradient(G4double val) { fSPSEne->SetGradient(val); }

  // Sets intercept for a linear distribution.
  void SetInterCept(G4double val) { fSPSEne->SetInterCept(val); }

  // Set position flag
  void SetPositionFlag(G4int val) { fPositionFlag = val; }

  void SetUniformTPC();

  void SetUniformGasPocket();

  void SetTPCCenter();

  void SetEnergyFileName(string);

  void SetIsG2(G4bool val) { IsG2 = val; }

  void SetGenInPhysVolumeName(G4String volName) { fPhysVolumeName = volName; }
  G4String GetGenInPhysVolumeName() { return fPhysVolumeName; }
 
 protected:
  G4ThreeVector fPosition;
  G4ThreeVector fDirection;
  G4double fEnergy;
  G4ParticleDefinition* fParticle;

 private:
  G4Navigator* gNavigator;

  G4int fNumberOfHits;
  G4int fPositionFlag;
  G4String fGeneratorName;  // Name of Generator.
  DSVGeneratorMessenger* fMessenger;
  G4SPSPosDistribution* fSPSPos;
  G4SPSAngDistribution* fSPSAng;
  G4SPSEneDistribution* fSPSEne;
  G4bool fVolumeFlag;
  G4bool fEnergyFlag;
  G4SPSRandomGenerator* fRndGen;
  G4bool IsVolumetric;
  G4bool IsMyEnergyDistribution;
  G4bool IsG2;
  vector<double> fMyEne;
  vector<double> fMyProb;
  G4double fCharge;
  G4int fNumberOfParticles;
  G4double fVerticalShift;
  G4String fPhysVolumeName;
};
#endif
/*
 * $Log: DSVGenerator.hh,v $
 * Revision 1.8  2016/03/08 10:58:14  davini
 * command for uniform distibution in LSV (and DSDetectorNeutronVetoG3 cleanup)
 *
 * Revision 1.7  2015/08/31 09:15:58  pagnes
 * spatial generators for DS20K updated
 *
 * Revision 1.6  2015/01/14 16:58:42  dfranco
 * added source vial and correspondent commands by Laura and Bianca; manual
 * updated
 *
 * Revision 1.5  2014/10/22 14:04:20  dfranco
 * add method to generate events in Liqud Argon
 *
 * Revision 1.4  2014/07/18 13:54:52  dfranco
 * Added a new position generator in the Gas Pocket; calibrated the numnber of
 * photons generated per electron in S2; added a new command to generate
 * multiple particles in the same position
 *
 * Revision 1.3  2014/05/21 10:28:19  dfranco
 * added the possibility to shoot ions
 *
 * Revision 1.2  2014/05/08 11:00:42  pagnes
 * Added generator in PMTs stem
 *
 * Revision 1.1  2014/05/07 12:20:55  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.8  2014/04/22 09:52:34  perassos
 * Added the option to simulate bgd from the G2 cryostats
 *
 * Revision 1.7  2014/04/11 10:20:43  perassos
 * Added generation in materials
 *
 * Revision 1.6  2014/02/13 12:03:21  dfranco
 * added new commands for spatial distributions and fixed the manual
 *
 * Revision 1.5  2013/10/15 16:39:37  dfranco
 * added the possibility to define the energy spectrum reading an ascii file
 *
 * Revision 1.4  2013/05/06 10:25:25  dfranco
 * Fixed cylindrical spatial distribution
 *
 * Revision 1.3  2013/05/02 16:36:14  dfranco
 * Added several position distributions
 *
 * Revision 1.2  2013/03/22 14:09:40  dfranco
 * Added the cvs logger code at the end of each file
 *
 */
