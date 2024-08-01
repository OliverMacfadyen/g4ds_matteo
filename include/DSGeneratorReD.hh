#ifndef DSGeneratorReD_HH
#define DSGeneratorReD_HH

#include "DSGeneratorReDMessenger.hh"
#include "DSPrimaryGeneratorAction.hh"
#include "DSVGenerator.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleMomentum.hh"
#include "G4ThreeVector.hh"
#include "G4VPrimaryGenerator.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4Event;
class G4ParticleGun;
class DSGeneratorReDMessenger;

class DSGeneratorReD : public DSVGenerator  // G4VUserPrimaryGeneratorAction
{
 public:
  DSGeneratorReD();
  ~DSGeneratorReD();

 public:
  void DSGeneratePrimaries(G4Event* event);  // anEvent
  G4double GetNeutronAngle(G4ThreeVector momentum);
  G4double GetNeutronEnergy(G4double neutronAngle);
  G4double FindNeutronEnergy() const { return fNeutronEnergy; }
  G4ThreeVector FindNeutronDirection() const { return fNeutronMomentum; }
  G4ThreeVector FindMainDirection() const { return fMainDirection; }
  G4double GetChargedParticleEnergy(G4double neutronAngle);
  G4double GetChargedParticleAngle(G4double neutronAngle);

  G4double GetBeamEnergy() { return fBeamEnergy; }
  G4double FindBeamEnergy() const { return fBeamEnergy; }
  void SetBeamEnergy(G4double beam) { fBeamEnergy = beam; }

  void SetConeOpening(G4double angle) {
    fConeOpening = angle;
    fCosConeOpening = std::cos(fConeOpening);
  }
  G4double GetConeOpening() const { return fConeOpening; }

  void SetReactionType(G4String atype);
  G4String GetReactionType() { return fReactionType; }
  G4String FindReactionType() const { return fReactionType; }

  /*
    inline void         SetTargetToTPCDistance(G4double val)     {
    fTargetTPCDistance = val; } inline G4double     GetTargetToTPCDistance() {
    return fTargetTPCDistance; }

    inline void         SetBeamToTPCAngleThetaX(G4double val)  { fThetaArX =
    val; } inline G4double     GetBeamToTPCAngleThetaX()              { return
    fThetaArX; }

    inline void         SetBeamToTPCAnglePhiX(G4double val)    { fPhiArX = val;
    } inline G4double     GetBeamToTPCAnglePhiX()                { return
    fPhiArX; }
*/

  inline void SetNeutronBeamThetaX(G4double val) { fNeutronAngleThetaX = val; }
  inline G4double GetNeutronBeamThetaX() { return fNeutronAngleThetaX; }

  inline void SetNeutronBeamPhiX(G4double val) { fNeutronAnglePhiX = val; }
  inline G4double GetNeutronBeamPhiX() { return fNeutronAnglePhiX; }

 private:
  void CalculateKinematics(G4double neutronAngle);
  void CalculateKinematicsDD(G4double neutronAngle);

  G4GeneralParticleSource* particleGun;
  // G4ParticleGun* particleGun;

  G4double fBeamEnergy;
  // G4double fTPCAngleThetaX;
  // G4double fphiAr;
  G4ThreeVector fMainDirection;
  G4ThreeVector fNeutronMomentum;

  G4double fConeOpening;
  G4double fCosConeOpening;

  G4double fNeutronAngle;
  G4double fNeutronEnergy;
  G4double fChargedEnergy;
  G4double fChargedAngle;

  G4bool isCalculated;
  G4bool isCalculatedDD;

  G4String fReactionType;

  G4double fTargetTPCDistance;
  // G4double fangleThetaX;
  // G4double fanglePhiX;
  G4double fNeutronAngleThetaX;
  G4double fNeutronAnglePhiX;

  DSGeneratorReDMessenger* fMessenger;
};

#endif
