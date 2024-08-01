// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// *********--***********************************************************
//

#include "DSGeneratorG4GunMessenger.hh"
#include "DSGeneratorG4Gun.hh"
#include "DSLogger.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Tokenizer.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"

using namespace std;

DSGeneratorG4GunMessenger::DSGeneratorG4GunMessenger(DSGeneratorG4Gun* gen) {
  generator = gen;
  fShootIon = false;
  fAtomicNumber = 0;
  fAtomicMass = 0;
  fIonCharge = 0;
  fIonExciteEnergy = 0.0;
  fIonEnergyLevel = 0;

  fDirectory = new G4UIdirectory("/ds/generator/g4gun/");
  fDirectory->SetGuidance("Control of DSG4Gun event generator");

  fPositionCmd = new G4UIcmdWith3VectorAndUnit("/ds/generator/g4gun/position", this);
  fPositionCmd->SetGuidance("Set the gun position");
  fPositionCmd->SetUnitCategory("Length");
  fPositionCmd->SetDefaultUnit("cm");
  fPositionCmd->SetUnitCandidates("mm cm m");

  fDirectionCmd = new G4UIcmdWith3Vector("/ds/generator/g4gun/direction", this);
  fDirectionCmd->SetGuidance("Set the gun direction");

  fEnergyCmd = new G4UIcmdWithADoubleAndUnit("/ds/generator/g4gun/energy", this);
  fEnergyCmd->SetGuidance("Set the gun energy");
  fEnergyCmd->SetUnitCategory("Energy");
  fEnergyCmd->SetDefaultUnit("MeV");
  fEnergyCmd->SetUnitCandidates("eV keV MeV GeV");

  fIonCmd = new G4UIcommand("/ds/generator/g4gun/ion", this);
  fIonCmd->SetGuidance("Set properties of ion to be generated.");
  fIonCmd->SetGuidance("[usage] /ds/generator/g4gun/ion Z A Q");
  fIonCmd->SetGuidance("	 Z:(int) AtomicNumber");
  fIonCmd->SetGuidance("	 A:(int) AtomicMass");
  fIonCmd->SetGuidance("	 Q:(int) Charge of Ion (in unit of e)");
  fIonCmd->SetGuidance("	 E:(double) Excitation energy (in keV)");

  G4UIparameter* param;
  param = new G4UIparameter("Z", 'i', false);
  fIonCmd->SetParameter(param);
  param = new G4UIparameter("A", 'i', false);
  fIonCmd->SetParameter(param);
  param = new G4UIparameter("Q", 'i', true);
  param->SetDefaultValue(-1);
  fIonCmd->SetParameter(param);
  param = new G4UIparameter("E", 'd', true);
  param->SetDefaultValue(0.0);
  fIonCmd->SetParameter(param);

  fIonLvlCmd = new G4UIcommand("/ds/generator/g4gun/ionL", this);
  fIonLvlCmd->SetGuidance("Set properties of ion to be generated.");
  fIonLvlCmd->SetGuidance("[usage] /ds/generator/g4gun/ionL Z A Q I");
  fIonLvlCmd->SetGuidance("	    Z:(int) AtomicNumber");
  fIonLvlCmd->SetGuidance("	    A:(int) AtomicMass");
  fIonLvlCmd->SetGuidance("	    Q:(int) Charge of Ion (in unit of e)");
  fIonLvlCmd->SetGuidance("	    I:(int) Level number of metastable state (0 = ground)");

  G4UIparameter* paraml;
  paraml = new G4UIparameter("Z", 'i', false);
  fIonLvlCmd->SetParameter(paraml);
  paraml = new G4UIparameter("A", 'i', false);
  fIonLvlCmd->SetParameter(paraml);
  paraml = new G4UIparameter("Q", 'i', true);
  paraml->SetDefaultValue(-1);
  fIonLvlCmd->SetParameter(paraml);
  paraml = new G4UIparameter("I", 'i', true);
  paraml->SetDefaultValue("0");
  fIonLvlCmd->SetParameter(paraml);
}

DSGeneratorG4GunMessenger::~DSGeneratorG4GunMessenger() {

  delete fDirectory;
  delete fPositionCmd;
  delete fDirectionCmd;
  delete fEnergyCmd;
  delete fParticleCmd;
  delete fIonCmd;
  delete fIonLvlCmd;
}

void DSGeneratorG4GunMessenger::SetNewValue(G4UIcommand* cmd, G4String newValue) {
  newValue = "";
  if (cmd == fPositionCmd) {
    ;  // generator->SetParticlePosition(fPositionCmd->ConvertToDimensioned3Vector(newValue));
  } else if (cmd == fDirectionCmd) {
    ;  // generator->SetParticleMomentumDirection(fDirectionCmd->ConvertTo3Vector(newValue));
  } else if (cmd == fEnergyCmd) {
    ;  // generator->SetParticleEnergy(fEnergyCmd->ConvertToDimensionedDouble(newValue));
  } else if (cmd == fParticleCmd) {
    ;  //  generator->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle(newValue));
  }
}
