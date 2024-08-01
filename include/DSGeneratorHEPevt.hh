//---------------------------------------------------------------------------//
//    Code by pablo.garcia@ciemat.es to read files in HEPEVT format:
//---------------------------------------------------------------------------//
#ifndef DSGeneratorHEPevt_h
#define DSGeneratorHEPevt_h 1

#include "DSVGenerator.hh"
#include "G4DataVector.hh"
#include "G4Navigator.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleMomentum.hh"
#include "G4VPrimaryGenerator.hh"
#include "Randomize.hh"
#include "globals.hh"

#include "DSGeneratorHEPevtMessenger.hh"

#include <fstream>

class DSGeneratorHEPevt : public DSVGenerator {

 public:
  DSGeneratorHEPevt();
  virtual ~DSGeneratorHEPevt();
  virtual void DSGeneratePrimaries(G4Event* event);

  void SetFileName(G4String stri) { fileName = stri; };
  G4String GetFileName() { return fileName; };

  // Set methods
  // void SetParticleDefinition(G4ParticleDefinition * aParticleDefinition);

 private:
  DSGeneratorHEPevtMessenger* fTheMessenger;
  // G4Navigator* gNavigator;

  G4String fileName;
  G4bool firstEvent;
  ifstream HEPfile;

  // particle properties
  // G4int                  NumberOfParticlesToBeGenerated;
  // G4ParticleDefinition* particle_definition;
  G4ParticleMomentum particle_momentum;
  G4double particle_energy;
  G4ThreeVector particle_position;
  G4double particle_time;
  // G4double               particle_charge;
  // G4ThreeVector          particle_polarization;

  G4ParticleDefinition* fParticle;
  G4ParticleTable* fParticleTable;

  // nevhep - event number
  // nhep - number of entries in this event record
  // isthep(..) - status code
  //    0 - null
  //    1 - final state particle
  //    2 - intermediate state
  //    3 - documentation line
  // idhep(..) - particle ID, P.D.G. standard
  // jmohep(0,..) - position of mother particle in list
  // jmohep(1,..) - position of second mother particle in list
  // jdahep(0,..) - position of first daughter in list
  // jdahep(1,..) - position of last daughter in list
  // phep(0,..) - x momentum in GeV/c
  // phep(1,..) - y momentum in GeV/c
  // phep(2,..) - z momentum in GeV/c
  // phep(3,..) - energy in GeV
  // phep(4,..) - mass in GeV/c**2
  // vhep(0,..) - x vertex position in mm
  // vhep(1,..) - y vertex position in mm
  // vhep(2,..) - z vertex position in mm
  // vhep(3,..) - production time in mm/c

  G4int nevhep;
  G4int nhep, isthep, idhep, jmohep1, jmohep2, jdahep1, jdahep2;
  G4double phep[5], vhep[4];
};

#endif
