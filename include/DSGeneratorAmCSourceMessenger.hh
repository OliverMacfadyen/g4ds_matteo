#ifndef _DSGENERATORAMCSOURCEMESSENGER_H
#define _DSGENERATORAMCSOURCEMESSENGER_H

#include "DSIO.hh"

#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UImessenger.hh"

class DSPrimaryGeneratorAction;
class G4UIcmdWithAString;
class G4UIcmdWith3Vector;
class DSGeneratorAmCSource;

using namespace std;

class DSGeneratorAmCSourceMessenger : public G4UImessenger {

 public:
  // default constructor
  DSGeneratorAmCSourceMessenger(DSGeneratorAmCSource*);

  // destructor
  virtual ~DSGeneratorAmCSourceMessenger();

  void SetNewValue(G4UIcommand*, G4String);

 private:
  DSGeneratorAmCSource* generator;
  G4UIdirectory* fDirectory;
  G4UIcmdWithAString* fAmCEnergyAngleCorrelationFileCmd;
  G4UIcmdWith3Vector* fSourceRotationCmd;
};

#endif
