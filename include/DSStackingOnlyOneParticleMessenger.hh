#ifndef DSStackingOnlyOneParticleMessenger_h
#define DSStackingOnlyOneParticleMessenger_h 1

#include "DSIO.hh"
#include "G4ParticleGun.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UImessenger.hh"

class DSStackingOnlyOneParticle;

using namespace std;

class DSStackingOnlyOneParticleMessenger : public G4UImessenger {

 public:
  DSStackingOnlyOneParticleMessenger(DSStackingOnlyOneParticle*);
  ~DSStackingOnlyOneParticleMessenger();

  void SetNewValue(G4UIcommand*, G4String);

 private:
  DSStackingOnlyOneParticle* stacking;
  G4UIdirectory* fDirectory;
  G4UIcmdWithAString* fParticleNameCmd;
  G4UIcmdWithAnInteger* fParticlePDGCmd;
};

#endif
