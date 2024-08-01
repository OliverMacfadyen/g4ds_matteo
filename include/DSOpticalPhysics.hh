#ifndef DSOpticalPhysics_h
#define DSOpticalPhysics_h 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

class G4VProcess;

class DSOpticalPhysics : public G4VPhysicsConstructor {
 public:
  DSOpticalPhysics(G4int verbose = 0, const G4String& name = "DSOptical");
  virtual ~DSOpticalPhysics();

 protected:
  // construct particle and physics
  virtual void ConstructParticle();
  virtual void ConstructProcess();

 private:
  DSOpticalPhysics(const DSOpticalPhysics& right);
  DSOpticalPhysics& operator=(const DSOpticalPhysics& right);
};
#endif
