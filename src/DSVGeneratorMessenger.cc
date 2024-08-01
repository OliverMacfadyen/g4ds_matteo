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

#include <iostream>
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "DSLogger.hh"
#include "DSVGenerator.hh"
#include "DSVGeneratorMessenger.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
using namespace std;

DSVGeneratorMessenger::DSVGeneratorMessenger(DSVGenerator* gen) {
  generator = gen;
  fDirectory = new G4UIdirectory("/ds/generator/");
  fDirectory->SetGuidance("Control of DSG4Gun event generator");

  fPositionCmd = new G4UIcmdWith3VectorAndUnit("/ds/generator/position", this);
  fPositionCmd->SetGuidance("Set the gun position");
  fPositionCmd->SetUnitCategory("Length");
  fPositionCmd->SetDefaultUnit("cm");
  fPositionCmd->SetUnitCandidates("mm cm m");

  fDirectionCmd = new G4UIcmdWith3Vector("/ds/generator/direction", this);
  fDirectionCmd->SetGuidance("Set the gun direction");

  fEnergyCmd = new G4UIcmdWithADoubleAndUnit("/ds/generator/energy", this);
  fEnergyCmd->SetGuidance("Set the gun energy");
  fEnergyCmd->SetUnitCategory("Energy");
  fEnergyCmd->SetDefaultUnit("MeV");
  fEnergyCmd->SetUnitCandidates("eV keV MeV GeV");

  fSphereBulkCmd = new G4UIcmdWithADoubleAndUnit("/ds/generator/sphere_radius", this);
  fSphereBulkCmd->SetGuidance("Bulk radius");
  fSphereBulkCmd->SetUnitCategory("Length");
  fSphereBulkCmd->SetDefaultUnit("cm");
  fSphereBulkCmd->SetUnitCandidates("mm cm m");

  fSphereBulkRadMinCmd = new G4UIcmdWithADoubleAndUnit("/ds/generator/sphere_radius_min", this);
  fSphereBulkRadMinCmd->SetGuidance("Bulk minimum radius");
  fSphereBulkRadMinCmd->SetGuidance("Default: 0");
  fSphereBulkRadMinCmd->SetUnitCategory("Length");
  fSphereBulkRadMinCmd->SetDefaultUnit("cm");
  fSphereBulkRadMinCmd->SetUnitCandidates("mm cm m");

  fSphereSurfCmd = new G4UIcmdWithADoubleAndUnit("/ds/generator/surface_radius", this);
  fSphereSurfCmd->SetGuidance("Surface radius");
  fSphereSurfCmd->SetUnitCategory("Length");
  fSphereSurfCmd->SetDefaultUnit("cm");
  fSphereSurfCmd->SetUnitCandidates("mm cm m");

  fSphereCentreCmd = new G4UIcmdWith3VectorAndUnit("/ds/generator/set_center", this);
  fSphereCentreCmd->SetGuidance("Set the sphere position");
  fSphereCentreCmd->SetGuidance("Default: 0. 0. 0. cm");
  fSphereCentreCmd->SetUnitCategory("Length");
  fSphereCentreCmd->SetDefaultUnit("cm");
  fSphereCentreCmd->SetUnitCandidates("mm cm m");

  fEnergyDisTypeCmd = new G4UIcmdWithAString("/ds/generator/dist_energy", this);
  G4String candidates = "Mono Lin Pow Exp Gauss Brem BBody Cdg";
  fEnergyDisTypeCmd->SetCandidates(candidates);

  fSetEminCmd = new G4UIcmdWithADoubleAndUnit("/ds/generator/emin", this);
  fSetEminCmd->SetDefaultUnit("MeV");
  fSetEminCmd->SetUnitCandidates("eV keV MeV GeV");

  fSetEmaxCmd = new G4UIcmdWithADoubleAndUnit("/ds/generator/emax", this);
  fSetEmaxCmd->SetDefaultUnit("MeV");
  fSetEmaxCmd->SetUnitCandidates("eV keV MeV GeV");

  fSetAlphaCmd = new G4UIcmdWithADouble("/ds/generator/alpha", this);

  fSetTempCmd = new G4UIcmdWithADoubleAndUnit("/ds/generator/temp", this);
  fSetTempCmd->SetDefaultUnit("kelvin");
  fSetTempCmd->SetUnitCandidates("kelvin");

  fSetEzeroCmd = new G4UIcmdWithADoubleAndUnit("/ds/generator/ezero", this);
  fSetEzeroCmd->SetDefaultUnit("MeV");
  fSetEzeroCmd->SetUnitCandidates("eV keV MeV GeV");

  fSetGradientCmd = new G4UIcmdWithADouble("/ds/generator/gradient", this);

  fSetInterCeptCmd = new G4UIcmdWithADouble("/ds/generator/intercept", this);

  fParticleCmd = new G4UIcmdWithAString("/ds/generator/particle", this);
  fParticleCmd->SetGuidance("Set the gun type");

  // confine to volume
  fConfineCmd = new G4UIcmdWithAString("/ds/generator/confine", this);
  fConfineCmd->SetGuidance("Confine source to volume (NULL to unset).");
  fConfineCmd->SetGuidance("usage: confine VolName");
  fConfineCmd->SetParameterName("VolName", true, true);
  fConfineCmd->SetDefaultValue("NULL");

  fDisPosTypeCmd = new G4UIcmdWithAString("/ds/generator/postype", this);
  candidates = "Point Plane Surface Volume";
  fDisPosTypeCmd->SetCandidates(candidates);

  fDisPosShapeCmd = new G4UIcmdWithAString("/ds/generator/posshape", this);
  candidates =
      "Square Circle Ellipse Rectangle Sphere Ellipsoid Cylinder "
      "Parallelepiped";
  fDisPosShapeCmd->SetCandidates(candidates);

  fSetHalfXCmd = new G4UIcmdWithADoubleAndUnit("/ds/generator/sethalfX", this);
  fSetHalfXCmd->SetUnitCandidates("mm cm m");

  fSetHalfYCmd = new G4UIcmdWithADoubleAndUnit("/ds/generator/sethalfY", this);
  fSetHalfYCmd->SetUnitCandidates("mm cm m");

  fSetHalfZCmd = new G4UIcmdWithADoubleAndUnit("/ds/generator/sethalfZ", this);
  fSetHalfZCmd->SetUnitCandidates("mm cm m");

  fRadiusCmd = new G4UIcmdWithADoubleAndUnit("/ds/generator/setradius", this);
  fRadiusCmd->SetUnitCandidates("mm cm m");

  fRadius0Cmd = new G4UIcmdWithADoubleAndUnit("/ds/generator/setradius0", this);
  fRadius0Cmd->SetUnitCandidates("mm cm m");

  fEnergyFileCmd = new G4UIcmdWithAString("/ds/generator/energyfile", this);

  fUniformTPCCmd = new G4UIcmdWithABool("/ds/generator/tpcdistribution", this);

  fUniformGasPocketCmd = new G4UIcmdWithABool("/ds/generator/gaspocketdistribution", this);

  fTPCCenterCmd = new G4UIcmdWithABool("/ds/generator/tpccenter", this);

  fGenInCryostatsCmd = new G4UIcmdWithAnInteger("/ds/generator/bgd_cryostats", this);

  fGenInLiquidArgonCmd = new G4UIcmdWithAnInteger("/ds/generator/liquidargon", this);

  fGenInSiPMCmd = new G4UIcmdWithAnInteger("/ds/generator/bgd_sipm", this);

  // fGenInOptoLinkCmd = new
  // G4UIcmdWithAnInteger("/ds/generator/bgd_optolink",this);

  fGenInPDMCmd = new G4UIcmdWithAnInteger("/ds/generator/bgd_pdm", this);

  // fGenInCooperVesselCmd = new
  // G4UIcmdWithAnInteger("/ds/generator/bgd_copper_vessel",this);

  fGenInTitaniumVesselCmd = new G4UIcmdWithAnInteger("/ds/generator/bgd_titanium_vessel", this);

  fGenInSiPMVetoCmd = new G4UIcmdWithAnInteger("/ds/generator/bgd_veto_sipm", this);

  fGenInCopperRingsCmd = new G4UIcmdWithAnInteger("/ds/generator/bgd_rings", this);

  fGenInGridCmd = new G4UIcmdWithAnInteger("/ds/generator/bgd_grid", this);

  fGenInGridSupportCmd = new G4UIcmdWithAnInteger("/ds/generator/bgd_gridSupport", this);
  
  fGenInTeflonCmd = new G4UIcmdWithAnInteger("/ds/generator/bgd_reflector", this);

  fGenInLateralTPBCmd = new G4UIcmdWithAnInteger("/ds/generator/bgd_lateral_tpb", this);

  fGenInLateralESRCmd = new G4UIcmdWithAnInteger("/ds/generator/bgd_lateral_esr", this);

  fGenInLateralAcrylicCmd = new G4UIcmdWithAnInteger("/ds/generator/bgd_lateral_reflector_acrylic", this);

  fGenInFusedSilicaCmd = new G4UIcmdWithAnInteger("/ds/generator/bgd_windows", this);

  fGenInPMTPhotocathodeCmd = new G4UIcmdWithAnInteger("/ds/generator/bgd_pmt_photocathode", this);

  fGenInPMTStemCmd = new G4UIcmdWithAnInteger("/ds/generator/bgd_pmt_stem", this);

  fGenInG2Cryostats = new G4UIcmdWithABool("/ds/generator/is_G2_cryostat", this);

  fGenInHolderSourceCmd = new G4UIcmdWithAnInteger("/ds/generator/holderSource_on", this);

  fGenInLiquidScintillatorCmd = new G4UIcmdWithAnInteger("/ds/generator/liquidscintillator", this);

  fGenInPlasticScintillatorCmd = new G4UIcmdWithAnInteger("/ds/generator/plasticscintillator", this);

  fGenInTPCBarrelCmd = new G4UIcmdWithAnInteger("/ds/generator/TPCbarrel", this);

  fGenInOPsAcrylicCmd = new G4UIcmdWithAnInteger("/ds/generator/OPsAcrylic", this);

  fGenInArgonBufferInsideCmd = new G4UIcmdWithAnInteger("/ds/generator/argonbufferinside", this);

  fGenInArgonBufferOutsideCmd = new G4UIcmdWithAnInteger("/ds/generator/argonbufferoutside", this);

  fGenInVetoPMTsCmd = new G4UIcmdWithAnInteger("/ds/generator/veto_pmts", this);

  fGenInVetoSteelSphereCmd = new G4UIcmdWithAnInteger("/ds/generator/veto_sphere", this);

  fGenInDS20kSteelStructureCmd = new G4UIcmdWithAnInteger("/ds/generator/steel_structure", this);

  fGenIn7mSphereCmd = new G4UIcmdWithAnInteger("/ds/generator/sphere_7m", this);

  fGenInDuneCryostatFoamCmd = new G4UIcmdWithAnInteger("/ds/generator/dune_cryo_foam", this);

  fGenInDuneCryostatMembraneCmd = new G4UIcmdWithAnInteger("/ds/generator/dune_cryo_membrane", this);

  fNumberOfParticlesCmd = new G4UIcmdWithAnInteger("/ds/generator/numberofparticles", this);

  fGenInArDMTestArgonCopperCmd = new G4UIcmdWithAnInteger("/ds/generator/bgd_ardm_test_argon_copper", this);

  fGenInArDMTestArgonNSCmd = new G4UIcmdWithAnInteger("/ds/generator/bgd_ardm_test_argon_ns", this);

  fGenInArDMTestArgonAcrylicCmd = new G4UIcmdWithAnInteger("/ds/generator/bgd_ardm_test_argon_acrylic", this);

  fGenInArDMTestArgonUndergroundCmd = new G4UIcmdWithAnInteger("/ds/generator/bgd_ardm_test_argon_underground", this);

  fGenInArDMTestArgonReflectorCmd = new G4UIcmdWithAnInteger("/ds/generator/bgd_ardm_test_argon_reflector", this);

  fGenInArDMTestArgonWLSCmd = new G4UIcmdWithAnInteger("/ds/generator/bgd_ardm_test_argon_wls", this);

  fGenInArDMDartDetCmd = new G4UIcmdWithAnInteger("/ds/generator/bgd_ardm_dart_det", this);
  fGenInArDMDartExtCmd = new G4UIcmdWithAnInteger("/ds/generator/bgd_ardm_dart_external", this);
  fGenInArDMDartSiPMCmd = new G4UIcmdWithAnInteger("/ds/generator/bgd_ardm_dart_sipm", this);

  fGenInGantryPipeCmd = new G4UIcmdWithAnInteger("/ds/generator/gantrypipe", this);

  fAriserSourceHolderCmd = new G4UIcmdWithAnInteger("/ds/generator/ariser_source", this);

  fGenInPhysVolumeCmd = new G4UIcmdWithAString("/ds/generator/physvolume", this);

  fGenInAnodeCmd = new G4UIcmdWithAnInteger("/ds/generator/anode", this);

  fGenInCathodeCmd = new G4UIcmdWithAnInteger("/ds/generator/cathode", this);

}

DSVGeneratorMessenger::~DSVGeneratorMessenger() {

  delete fDirectory;
  delete fPositionCmd;
  delete fDirectionCmd;
  delete fEnergyCmd;
  delete fParticleCmd;
  delete fSphereBulkCmd;
  delete fSphereSurfCmd;
  delete fSphereCentreCmd;
  delete fConfineCmd;
  delete fEnergyDisTypeCmd;
  delete fSetEminCmd;
  delete fSetEmaxCmd;
  delete fSetAlphaCmd;
  delete fSetTempCmd;
  delete fSetEzeroCmd;
  delete fSetGradientCmd;
  delete fSetInterCeptCmd;
  delete fDisPosTypeCmd;
  delete fDisPosShapeCmd;
  delete fSetHalfXCmd;
  delete fSetHalfYCmd;
  delete fSetHalfZCmd;
  delete fRadiusCmd;
  delete fRadius0Cmd;
  delete fEnergyFileCmd;
  delete fGenInCryostatsCmd;
  delete fGenInLiquidArgonCmd;
  delete fGenInTeflonCmd;
  delete fGenInPMTStemCmd;
  delete fGenInGridCmd;
  delete fGenInGridSupportCmd;
  delete fGenInCopperRingsCmd;
  delete fGenInSiPMCmd;
  // delete fGenInOptoLinkCmd;
  delete fGenInPDMCmd;
  delete fGenInSiPMVetoCmd;
  delete fGenInFusedSilicaCmd;
  delete fGenInPMTPhotocathodeCmd;
  delete fGenInG2Cryostats;
  delete fGenInHolderSourceCmd;
  delete fGenInLiquidScintillatorCmd;
  delete fGenInPlasticScintillatorCmd;
  delete fGenInTPCBarrelCmd;
  delete fGenInOPsAcrylicCmd;
  delete fGenInArgonBufferInsideCmd;
  delete fGenInArgonBufferOutsideCmd;
  delete fNumberOfParticlesCmd;
  delete fUniformGasPocketCmd;
  delete fGenInVetoPMTsCmd;
  delete fGenInVetoSteelSphereCmd;
  delete fGenInArDMTestArgonAcrylicCmd;
  delete fGenInArDMTestArgonUndergroundCmd;
  delete fGenInArDMTestArgonNSCmd;
  delete fGenInArDMTestArgonReflectorCmd;
  delete fGenInArDMTestArgonWLSCmd;
  delete fGenInArDMTestArgonCopperCmd;
  delete fGenInArDMDartDetCmd;
  delete fGenInArDMDartExtCmd;
  delete fGenInArDMDartSiPMCmd;
  delete fGenInDS20kSteelStructureCmd;
  delete fGenInDuneCryostatFoamCmd;
  delete fGenInDuneCryostatMembraneCmd;
  delete fGenIn7mSphereCmd;
  delete fGenInGantryPipeCmd;
  // delete fGenInCooperVesselCmd;
  delete fGenInTitaniumVesselCmd;
  delete fAriserSourceHolderCmd;
  delete fGenInLateralTPBCmd;
  delete fGenInLateralESRCmd;
  delete fGenInLateralAcrylicCmd;
  delete fGenInPhysVolumeCmd;
  delete fGenInAnodeCmd;
  delete fGenInCathodeCmd;

}

void DSVGeneratorMessenger::SetNewValue(G4UIcommand* cmd, G4String newValue) {
  if (cmd == fPositionCmd) {
    generator->SetVParticlePosition(fPositionCmd->ConvertToDimensioned3Vector(newValue));
  } else if (cmd == fDirectionCmd) {
    generator->SetVParticleDirection(fDirectionCmd->ConvertTo3Vector(newValue));
  } else if (cmd == fEnergyCmd) {
    generator->SetVParticleEnergy(fEnergyCmd->ConvertToDimensionedDouble(newValue));
  } else if (cmd == fParticleCmd) {
    if (newValue == "Ar40") generator->SetVParticleDefinition(G4IonTable::GetIonTable()->GetIon(18, 40, 0, 0));
    else
      generator->SetVParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle(newValue));
  } else if (cmd == fSphereBulkCmd) {
    generator->SetPositionFlag(2);
    generator->SetIsVolumeDistribution(true);
    generator->SetPosDisType("Volume");
    generator->SetPosDisShape("Sphere");
    generator->SetRadius(fSphereBulkCmd->ConvertToDimensionedDouble(newValue));
  } else if (cmd == fSphereSurfCmd) {
    generator->SetPositionFlag(2);
    generator->SetIsVolumeDistribution(true);
    generator->SetPosDisType("Surface");
    generator->SetPosDisShape("Sphere");
    generator->SetRadius(fSphereSurfCmd->ConvertToDimensionedDouble(newValue));
  } else if (cmd == fSphereCentreCmd) {
    generator->SetCentreCoords(fSphereCentreCmd->ConvertToDimensioned3Vector(newValue));
  } else if (cmd == fSphereBulkRadMinCmd) {
    generator->SetRadius0(fSphereBulkRadMinCmd->ConvertToDimensionedDouble(newValue));
  } else if (cmd == fConfineCmd) {
    generator->ConfineSourceToVolume(newValue);
  } else if (cmd == fEnergyDisTypeCmd) {
    DSLog(routine) << "G4Gun Energy distribution " << newValue << endlog;
    generator->SetIsEnergyDistribution(true);
    generator->SetEnergyDisType(newValue);
  } else if (cmd == fSetEminCmd) {
    DSLog(routine) << "G4Gun Emin = " << newValue << endlog;
    generator->SetEmin(fSetEminCmd->ConvertToDimensionedDouble(newValue));
  } else if (cmd == fSetEmaxCmd) {
    DSLog(routine) << "G4Gun Emax = " << newValue << endlog;
    generator->SetEmax(fSetEmaxCmd->ConvertToDimensionedDouble(newValue));
  } else if (cmd == fSetAlphaCmd) {
    DSLog(routine) << "G4Gun Alpha = " << newValue << endlog;
    generator->SetAlpha(fSetAlphaCmd->ConvertToDouble(newValue));
  } else if (cmd == fSetTempCmd) {
    DSLog(routine) << "G4Gun Temperature = " << newValue << endlog;
    generator->SetTemp(fSetTempCmd->ConvertToDimensionedDouble(newValue));
  } else if (cmd == fSetEzeroCmd) {
    DSLog(routine) << "G4Gun Ezero = " << newValue << endlog;
    generator->SetEzero(fSetEzeroCmd->ConvertToDimensionedDouble(newValue));
  } else if (cmd == fSetGradientCmd) {
    DSLog(routine) << "G4Gun Gradient = " << newValue << endlog;
    generator->SetGradient(fSetGradientCmd->ConvertToDouble(newValue));
  } else if (cmd == fSetInterCeptCmd) {
    DSLog(routine) << "G4Gun Intercept = " << newValue << endlog;
    generator->SetInterCept(fSetInterCeptCmd->ConvertToDouble(newValue));
  } else if (cmd == fDisPosTypeCmd) {
    DSLog(routine) << "Position Distribution Type = " << newValue << endlog;
    generator->SetPositionFlag(2);
    generator->SetIsVolumeDistribution(true);
    generator->SetPosDisType(newValue);
  } else if (cmd == fDisPosShapeCmd) {
    DSLog(routine) << "Position Distribution Shape = " << newValue << endlog;
    generator->SetPositionFlag(2);
    generator->SetIsVolumeDistribution(true);
    generator->SetPosDisShape(newValue);
  } else if (cmd == fSetHalfXCmd) {
    DSLog(routine) << "HalfX = " << newValue << endlog;
    generator->SetHalfX(fSetHalfXCmd->ConvertToDimensionedDouble(newValue));
  } else if (cmd == fSetHalfYCmd) {
    DSLog(routine) << "HalfY = " << newValue << endlog;
    generator->SetHalfY(fSetHalfYCmd->ConvertToDimensionedDouble(newValue));
  } else if (cmd == fSetHalfZCmd) {
    DSLog(routine) << "Half Z = " << newValue << endlog;
    generator->SetHalfZ(fSetHalfZCmd->ConvertToDimensionedDouble(newValue));
  } else if (cmd == fRadiusCmd) {
    DSLog(routine) << "Radius = " << newValue << endlog;
    generator->SetRadius(fRadiusCmd->ConvertToDimensionedDouble(newValue));
  } else if (cmd == fRadius0Cmd) {
    DSLog(routine) << "Radius0 = " << newValue << endlog;
    generator->SetRadius0(fRadius0Cmd->ConvertToDimensionedDouble(newValue));
  } else if (cmd == fEnergyFileCmd) {
    DSLog(routine) << "Random Energy File Name  = " << newValue << endlog;
    generator->SetEnergyFileName(newValue);
  } else if (cmd == fUniformTPCCmd) {
    if (fUniformTPCCmd->ConvertToBool(newValue)) {
      generator->SetIsVolumeDistribution(true);
      DSLog(routine) << "Random Position Distribution in the TPC  = " << newValue << endlog;
      generator->SetUniformTPC();
    }
  // } else if (cmd == fUniformGasPocketCmd) {
  //   if (fUniformGasPocketCmd->ConvertToBool(newValue)) {
  //     generator->SetIsVolumeDistribution(true);
  //     DSLog(routine) << "Random Position Distribution in the Gas Pocket  = " << newValue << endlog;
  //     generator->SetUniformGasPocket();
  //   }
  } else if (cmd == fTPCCenterCmd) {
    if (fTPCCenterCmd->ConvertToBool(newValue)) {
      DSLog(routine) << "Set position or distribution in the TPC center  = " << newValue << endlog;
      generator->SetTPCCenter();
    }
  } else if (cmd == fGenInPhysVolumeCmd) {
    DSLog(routine) << "Generation confined to physical volume = " << newValue << endlog;
    generator->SetGenInPhysVolumeName(newValue);  
    generator->SetPositionFlag(100);
  }  else if (cmd == fGenInCryostatsCmd) {
    generator->SetPositionFlag(3);
  } else if (cmd == fGenInTeflonCmd) {
    generator->SetPositionFlag(4);
  } else if (cmd == fGenInFusedSilicaCmd) {
    generator->SetPositionFlag(5);
  } else if (cmd == fGenInPMTPhotocathodeCmd) {
    generator->SetPositionFlag(6);
  } else if (cmd == fGenInPMTStemCmd) {
    generator->SetPositionFlag(7);
  } else if (cmd == fGenInLiquidArgonCmd) {
    generator->SetIsVolumeDistribution(true);
    generator->SetPositionFlag(8);
  } else if (cmd == fGenInG2Cryostats) {
    generator->SetIsG2(fGenInG2Cryostats->ConvertToBool(newValue));
  } else if (cmd == fGenInHolderSourceCmd) {
    generator->SetPositionFlag(9);
  } else if (cmd == fGenInSiPMCmd) {
    generator->SetPositionFlag(10);
    // } else if (cmd == fGenInOptoLinkCmd){
    //       generator->SetPositionFlag(10);
    //       generator->SetVerticalShift(4*cm);
  } else if (cmd == fGenInGridCmd) {
    generator->SetPositionFlag(11);
  } else if (cmd == fGenInCopperRingsCmd) {
    generator->SetPositionFlag(12);
  } else if (cmd == fGenInLiquidScintillatorCmd) {
    generator->SetPositionFlag(13);
  } else if (cmd == fGenInVetoPMTsCmd) {
    generator->SetPositionFlag(14);
  } else if (cmd == fGenInVetoSteelSphereCmd) {
    generator->SetPositionFlag(15);
  } else if (cmd == fNumberOfParticlesCmd) {
    generator->SetVNumberOfParticles(fNumberOfParticlesCmd->ConvertToInt(newValue));
  } else if (cmd == fGenInArDMTestArgonCopperCmd) {
    generator->SetPositionFlag(16);
  } else if (cmd == fGenInArDMDartDetCmd) {
    generator->SetPositionFlag(17);
  } else if (cmd == fGenInArDMDartExtCmd) {
    generator->SetPositionFlag(18);
  } else if (cmd == fGenInArDMDartSiPMCmd) {
    generator->SetPositionFlag(19);
  } else if (cmd == fGenInArDMTestArgonNSCmd) {
    generator->SetPositionFlag(20);
  } else if (cmd == fGenInArDMTestArgonAcrylicCmd) {
    generator->SetPositionFlag(21);
  } else if (cmd == fGenInArDMTestArgonReflectorCmd) {
    generator->SetPositionFlag(22);
  } else if (cmd == fGenInArDMTestArgonWLSCmd) {
    generator->SetPositionFlag(23);
  } else if (cmd == fGenInArDMTestArgonUndergroundCmd) {
    generator->SetPositionFlag(24);
  } else if (cmd == fGenInDS20kSteelStructureCmd) {
    generator->SetPositionFlag(25);
  } else if (cmd == fGenIn7mSphereCmd) {
    generator->SetPositionFlag(26);
  } else if (cmd == fGenInPlasticScintillatorCmd) {
    generator->SetPositionFlag(27);
  } else if (cmd == fGenInArgonBufferInsideCmd) {
    generator->SetIsVolumeDistribution(true);
    generator->SetPositionFlag(28);
  } else if (cmd == fGenInArgonBufferOutsideCmd) {
    generator->SetIsVolumeDistribution(true);
    generator->SetPositionFlag(29);
  } else if (cmd == fGenInSiPMVetoCmd) {
    generator->SetPositionFlag(30);
  } else if (cmd == fGenInDuneCryostatFoamCmd) {
    generator->SetPositionFlag(31);
  } else if (cmd == fGenInDuneCryostatMembraneCmd) {
    generator->SetPositionFlag(32);
  } else if (cmd == fGenInGantryPipeCmd) {
    generator->SetIsVolumeDistribution(true);
    generator->SetPositionFlag(33);
  } else if (cmd == fGenInTitaniumVesselCmd) {
    generator->SetPositionFlag(34);
  } else if (cmd == fUniformGasPocketCmd) {
    generator->SetIsVolumeDistribution(true);
    DSLog(routine) << "Random Position Distribution in the Gas Pocket  = " << newValue << endlog;
    generator->SetPositionFlag(35);
    //generator->SetUniformGasPocket();
  } else if (cmd == fAriserSourceHolderCmd) {
    generator->SetIsVolumeDistribution(true);
    DSLog(routine) << "Random Position Distribution in the ARIS-ER source holder" << endlog;
    generator->SetPositionFlag(36);
  } else if (cmd == fGenInLateralTPBCmd) {
    generator->SetPositionFlag(37);
  } else if (cmd == fGenInLateralESRCmd) {
    generator->SetPositionFlag(38);
  } else if (cmd == fGenInLateralAcrylicCmd) {
    generator->SetPositionFlag(39);
  } else if (cmd == fGenInPDMCmd) {
    generator->SetPositionFlag(40);
  } else if (cmd == fGenInGridSupportCmd) {
    generator->SetPositionFlag(41);
  } else if (cmd == fGenInTPCBarrelCmd) {
    generator->SetPositionFlag(42);
  } else if (cmd == fGenInAnodeCmd) {
    generator->SetPositionFlag(43);
  } else if (cmd == fGenInCathodeCmd) {
    generator->SetPositionFlag(44);
  } else if (cmd == fGenInOPsAcrylicCmd) {
    generator->SetPositionFlag(45);
  }
  
}
