// --------------------------------------------------------------------------//
/**
 * AUTHOR: D. Franco
 * CONTACT: davide.franco@mi.infn.it
 */
// --------------------------------------------------------------------------//

#ifndef _BXGENERATORRDMDecayGun_HH
#define _BXGENERATORRDMDecayGun_HH

//---------------------------------------------------------------------------//

#include "DSGeneratorRDMNucleus.hh"
#include "DSVGenerator.hh"
#include "G4Event.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleDefinition.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"

//---------------------------------------------------------------------------//

class DSGeneratorRDMDecayGunMessenger;
class DSGeneratorRDMDecayGun : public DSVGenerator {
 public:
  // default constructor
  DSGeneratorRDMDecayGun();

  // copy constructor
  // DSGeneratorRDMDecayGun(const DSGeneratorRDMDecayGun &);

  // destructor
  virtual ~DSGeneratorRDMDecayGun();

  // Sets the isotope.
  void SetNucleus(DSGeneratorRDMNucleus theIon1);

  // Returns the specified isotope.
  DSGeneratorRDMNucleus GetNucleus() { return theIon; }

  // public interface
  virtual void DSGeneratePrimaries(G4Event* event);

  // protected members

  // private  members
 private:
  DSGeneratorRDMDecayGunMessenger* fTheMessenger;
  G4ParticleGun* fParticleGun;

  DSGeneratorRDMNucleus theIon;
  G4int fA;
  G4int fZ;
  G4double fE;
  G4bool IsFirst;
  G4ParticleDefinition* aIon;
};
#endif

/*
 * Revision 1.1  2013/03/22 13:50:09  dfranco
 * Added radioactive decay and radioactive decay chain generators. They require
 * the correspondent stacking actions. Two mac files are included as examples
 *
 */
