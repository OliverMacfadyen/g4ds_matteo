#ifndef DSGeneratorBiasMessenger_HH
#define DSGeneratorBiasMessenger_HH 1

#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UImessenger.hh"

#include "DSGeneratorBias.hh"
#include "DSIO.hh"

class G4UIcmdWithADouble;

class DSGeneratorBias;

class DSGeneratorBiasMessenger : public G4UImessenger {

 public:
  DSGeneratorBiasMessenger(DSGeneratorBias*);
  ~DSGeneratorBiasMessenger();

  void SetNewValue(G4UIcommand*, G4String);

 private:
  DSGeneratorBias*           generator;
  G4UIdirectory*             fDirectory;
  G4UIcmdWithAString*        fFileCmd;
  G4UIcmdWithADouble*        fSkipEventsCmd;
  G4UIcmdWithAnInteger*      fAmplificationFactorCmd;
  G4UIcmdWithADouble*        fDirectionSmearingCmd;
  G4UIcmdWithADouble*        fEnergySmearingCmd;
  G4UIcmdWithADoubleAndUnit* fMinimumEnergyCmd;

};

#endif
/*
 * Revision 1.1  2013/03/20 09:54:28  dfranco
 * New version of g4ds
 *
 *
 */
