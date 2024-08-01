#ifndef _DSGENERATORAMCSOURCE_HH
#define _DSGENERATORAMCSOURCE_HH

//---------------------------------------------------------------------------//

#include <iostream>
#include <sstream>
#include <string>
#include "DSGeneratorAmCSourceMessenger.hh"
#include "DSVGenerator.hh"
#include "G4Event.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"
// #include "TFile.h"
// #include "TH2F.h"
//---------------------------------------------------------------------------//

class DSGeneratorAmCSource : public DSVGenerator {
 public:
  // default constructor
  DSGeneratorAmCSource();

  // destructor
  virtual ~DSGeneratorAmCSource();

  // public interface
  virtual void DSGeneratePrimaries(G4Event* event);

  void SetAmCSourceRotation(G4ThreeVector var) { fSourceRotation = var; }
  G4ThreeVector GetAmCSourceRotation() { return fSourceRotation; }
  void SetAmCEnergyAngleCorrelationFileName(string);

 private:
  DSGeneratorAmCSourceMessenger* fTheMessenger;

  G4ParticleTable* fParticleTable;
  G4ParticleDefinition* fParticle;
  G4ThreeVector fPosition;
  G4ThreeVector fDirection;
  G4ThreeVector fSourceRotation;

  G4bool isFirstTime;
  //TH2F* fenergyanglecorrelation2d;
  //linearized 2D histo. bincenter x, bincenter y, cumulative of the Z values
  double fenergyanglecorrelation2d_x[8000];
  double fenergyanglecorrelation2d_y[8000];
  double fenergyanglecorrelation2d_z[8000];

};
#endif
