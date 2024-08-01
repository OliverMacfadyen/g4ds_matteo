//
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
// ********************************************************************

#include "DSDetectorConstructionMessenger.hh"
#include "DSDetectorConstruction.hh"
#include "DSIO.hh"
#include "DSLogger.hh"
#include "DSMaterial.hh"
#include "DSParameters.hh"
#include "DSStorage.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIdirectory.hh"
class DSDetectorConstructionMessenger;

using namespace std;

DSDetectorConstructionMessenger::DSDetectorConstructionMessenger(DSDetectorConstruction* det) {
  fDetectorConstruction = det;
  fDirectory = new G4UIdirectory("/ds/detector/");
  fDirectory->SetGuidance("Control detector gemetry and materials");

  fConfigurationCmd = new G4UIcmdWithAString("/ds/detector/configuration", this);
  G4String candidate_experiments = "ds50_full ds50_tpc_veto ds50_tpc ds20k ds20k_tpc dsproto dsproto0 red_tpc dart ardm licorne ariser tester";
  fConfigurationCmd->SetCandidates(candidate_experiments);

  fScintillatorCmd = new G4UIcmdWithAString("/ds/detector/scintillator", this);
  G4String candidates =
      "NatLi6Scintillator OrtoCarbScintillator "
      "Li6Scintillator BoronScintillator GdScintillator";
  fScintillatorCmd->SetCandidates(candidates);

  fCTFfillCmd = new G4UIcmdWithAString("/ds/detector/wt_material", this);
  G4String candidates2 = "BoronScintillator GdWater Water Air";
  fCTFfillCmd->SetCandidates(candidates2);

  fSourceCmd = new G4UIcmdWith3VectorAndUnit("/ds/detector/source_position", this);
  fSourceCmd->SetUnitCandidates("mm cm m");

  fExtLArScintillatingCmd = new G4UIcmdWithABool("/ds/detector/ExtLarScintillating", this);

  fCopperRingsCmd = new G4UIcmdWithABool("/ds/detector/real_field_rings", this);

  fVetoYieldFactorCmd = new G4UIcmdWithADouble("/ds/detector/vetoyieldfactor", this);

  fSourceHolderCenterCmd = new G4UIcmdWith3VectorAndUnit("/ds/detector/sourceholdercenter", this);
  fSourceHolderCenterCmd->SetGuidance("Set the center of source holder");
  fSourceHolderCenterCmd->SetDefaultUnit("cm");
  fSourceHolderCenterCmd->SetUnitCandidates("mm cm m");

  fSourceHolderThetaCmd = new G4UIcmdWithADoubleAndUnit("/ds/detector/sourceholdertheta", this);
  fSourceHolderThetaCmd->SetGuidance("Source Holder rotates around X axis");
  fSourceHolderThetaCmd->SetUnitCandidates("degree deg rad");
  fSourceHolderThetaCmd->SetDefaultUnit("deg");

  fSourceHolderPhiCmd = new G4UIcmdWithADoubleAndUnit("/ds/detector/sourceholderphi", this);
  fSourceHolderPhiCmd->SetGuidance("Source Holder rotates around Z axis");
  fSourceHolderPhiCmd->SetUnitCandidates("degree deg rad");
  fSourceHolderPhiCmd->SetDefaultUnit("deg");

  fSourceHolderLeadCmd = new G4UIcmdWithABool("/ds/detector/sourceholderlead", this);
  fSourceHolderLeadCmd->SetGuidance("Turn on/off the Lead Shield in Source Holder");

  fDS20kCryoWallThickCmd = new G4UIcmdWithADoubleAndUnit("/ds/detector/ds20cryo_thickness", this);
  fDS20kCryoWallThickCmd->SetUnitCandidates("mm cm");

  fArDMLeadHeightCmd = new G4UIcmdWithADoubleAndUnit("/ds/detector/ardm_lead", this);
  fArDMLeadHeightCmd->SetUnitCandidates("mm cm m");

  fDS20ReflectorCmd = new G4UIcmdWithAString("/ds/detector/reflector", this);
  G4String candidates3 = "baseline alternate";
  fDS20ReflectorCmd->SetCandidates(candidates3);

  fDS20CryoMaterialCmd = new G4UIcmdWithAnInteger("/ds/detector/ds20cryo_material", this);
  fDS20WindowMaterialCmd = new G4UIcmdWithAnInteger("/ds/detector/ds20window_material", this);

  fDS20kCryoCornerDistanceCmd = new G4UIcmdWithADoubleAndUnit("/ds/detector/ds20cryo_distance", this);
  fDS20kCryoCornerDistanceCmd->SetUnitCandidates("mm cm");

  fDS20kTPCheightCmd = new G4UIcmdWithADoubleAndUnit("/ds/detector/ds20cryo_tpcHeight", this);
  fDS20kTPCheightCmd->SetUnitCandidates("mm cm m");

  fDS20kTPCedgeCmd = new G4UIcmdWithADoubleAndUnit("/ds/detector/ds20cryo_tpcEdge", this);
  fDS20kTPCedgeCmd->SetUnitCandidates("mm cm m");

  fGdLayerThicknessCmd = new G4UIcmdWithADoubleAndUnit("/ds/detector/ds20k_gdFoil", this);
  fGdLayerThicknessCmd->SetUnitCandidates("mm cm");

  fDS20kfLArThicknessAboveTPCCmd = new G4UIcmdWithADoubleAndUnit("/ds/detector/ds20cryo_LArLevel", this);
  fDS20kfLArThicknessAboveTPCCmd->SetUnitCandidates("mm cm m");

  fIsDS20kLSVCmd = new G4UIcmdWithABool("/ds/detector/ds20lsv_detector", this);

  fIDiameterCmd = new G4UIcmdWithADoubleAndUnit("/ds/detector/ds20lsv_diameter", this);
  fIDiameterCmd->SetUnitCandidates("m cm");

  // Acrylic vessel for very first Veto design (LSV with PMTs in water)
  fAcrylicVesselDiameterCmd = new G4UIcmdWithADoubleAndUnit("/ds/detector/ds20acrylicVessel_diameter", this);
  fAcrylicVesselDiameterCmd->SetUnitCandidates("m cm");

  // set distance of PMTs layers from anode and cathode
  fDS20kSiPMOffsetCmd = new G4UIcmdWithADoubleAndUnit("/ds/detector/ds20sipm_offset", this);
  fDS20kSiPMOffsetCmd->SetUnitCandidates("mm cm");

  fDS20kTPCVesselThicknessCmd = new G4UIcmdWithADoubleAndUnit("/ds/detector/ds20_AcrylicWalls_Thick", this);
  fDS20kTPCVesselThicknessCmd->SetUnitCandidates("mm cm");

  fDS20kLArBufferThicknessCmd = new G4UIcmdWithADoubleAndUnit("/ds/detector/ds20_LArBuffers_Thick", this);
  fDS20kLArBufferThicknessCmd->SetUnitCandidates("mm cm");

  fDS20kGdPlasticThicknessCmd = new G4UIcmdWithADoubleAndUnit("/ds/detector/ds20_VetoShell_Thick", this);
  fDS20kGdPlasticThicknessCmd->SetUnitCandidates("mm cm");

  fDS20kHDPEShellThicknessCmd = new G4UIcmdWithADoubleAndUnit("/ds/detector/ds20_HDPEShell_Thick", this);
  fDS20kHDPEShellThicknessCmd->SetUnitCandidates("mm cm");

  fDS20kHDPEShellCapThicknessCmd = new G4UIcmdWithADoubleAndUnit("/ds/detector/ds20_HDPEShellCap_Thick", this);
  fDS20kHDPEShellCapThicknessCmd->SetUnitCandidates("mm cm");

  fDS20kTPBThicknessCmd = new G4UIcmdWithADoubleAndUnit("/ds/detector/ds20_TPB_Thick", this);
  fDS20kTPBThicknessCmd->SetUnitCandidates("mm cm");

  fDS20kTPBLayersCmd = new G4UIcmdWithAString("/ds/detector/ds20_TPB_Layers", this);

  fWTankRCmd = new G4UIcmdWithADoubleAndUnit("/ds/detector/wt_radius", this);
  fWTankRCmd->SetUnitCandidates("m cm");

  fWTankHeightCmd = new G4UIcmdWithADoubleAndUnit("/ds/detector/wt_height", this);
  fWTankHeightCmd->SetUnitCandidates("m cm");

  // G4UIdirectory *licorne_directory = fDirectory = new
  // G4UIdirectory("/ds/detector/licorne");

  fLicNDRadiusCmd = new G4UIcmdWithADoubleAndUnit("/ds/detector/licorne/distance", this);
  fLicNDRadiusCmd->SetUnitCandidates("m cm");

  fLicNDThetaCmd = new G4UIcmdWithADoubleAndUnit("/ds/detector/licorne/theta", this);
  fLicNDThetaCmd->SetUnitCandidates("degree deg rad");

  fLicNDPhi1Cmd = new G4UIcmdWithADoubleAndUnit("/ds/detector/licorne/phi1", this);
  fLicNDPhi1Cmd->SetUnitCandidates("degree deg rad");

  fLicNDPhi2Cmd = new G4UIcmdWithADoubleAndUnit("/ds/detector/licorne/phi2", this);
  fLicNDPhi2Cmd->SetUnitCandidates("degree deg rad");

  fLicNDNuclRecEnergyCmd = new G4UIcmdWithADoubleAndUnit("/ds/detector/licorne/nuclear_energy", this);
  fLicNDNuclRecEnergyCmd->SetUnitCandidates("keV MeV");

  fLicNeutronEnergyCmd = new G4UIcmdWithADoubleAndUnit("/ds/detector/licorne/neutron_energy", this);
  fLicNeutronEnergyCmd->SetUnitCategory("Energy");
  fLicNeutronEnergyCmd->SetDefaultUnit("MeV");
  fLicNeutronEnergyCmd->SetUnitCandidates("keV MeV");

  fLicActivateWallCmd = new G4UIcmdWithABool("/ds/detector/licorne/activate_wall", this);

  fLicWallThickCmd = new G4UIcmdWithADoubleAndUnit("/ds/detector/licorne/wall_thick", this);
  fLicWallThickCmd->SetUnitCandidates("m cm");

  fLicWallWindowLengthCmd = new G4UIcmdWithADoubleAndUnit("/ds/detector/licorne/wall_window_length", this);
  fLicWallWindowLengthCmd->SetUnitCandidates("m cm");

  fLicWallDistanceToTPCCmd = new G4UIcmdWithADoubleAndUnit("/ds/detector/licorne/wall_distance_to_TPC", this);
  fLicWallDistanceToTPCCmd->SetUnitCandidates("m cm");

  fLicExpHallCmd = new G4UIcmdWithABool("/ds/detector/licorne/exp_hall", this);

  fAriserSourceDistanceToCryoCmd = new G4UIcmdWithADoubleAndUnit("/ds/detector/ariser/source_distance_to_cryo", this);
  fAriserSourceDistanceToCryoCmd->SetUnitCategory("Length");

  fAriserBEGeDistanceToCryoCmd = new G4UIcmdWithADoubleAndUnit("/ds/detector/ariser/bege_distance_to_cryo", this);
  fAriserBEGeDistanceToCryoCmd->SetUnitCategory("Length");

  fAriserBEGeAngleCmd = new G4UIcmdWithADoubleAndUnit("/ds/detector/ariser/bege_angle", this);
  fAriserBEGeAngleCmd->SetUnitCategory("Angle");

  dart_directory = new G4UIdirectory("/ds/detector/dart/");
  fDSDartTeflonThickness = new G4UIcmdWithADoubleAndUnit("/ds/detector/dart/teflon_thickness", this);
  fDSDartTeflonThickness->SetUnitCandidates("m cm");
  fDSDartTeflonHeight = new G4UIcmdWithADoubleAndUnit("/ds/detector/dart/teflon_height", this);
  fDSDartTeflonHeight->SetUnitCandidates("m cm");
  fDSDartTeflonRadius = new G4UIcmdWithADoubleAndUnit("/ds/detector/dart/teflon_radius", this);
  fDSDartTeflonRadius->SetUnitCandidates("m cm");

  fDSDartDewarThickness = new G4UIcmdWithADoubleAndUnit("/ds/detector/dart/dewar_thickness", this);
  fDSDartDewarThickness->SetUnitCandidates("m cm");
  fDSDartDewarHeight = new G4UIcmdWithADoubleAndUnit("/ds/detector/dart/dewar_height", this);
  fDSDartDewarHeight->SetUnitCandidates("m cm");
  fDSDartDewarRadius = new G4UIcmdWithADoubleAndUnit("/ds/detector/dart/dewar_radius", this);
  fDSDartDewarRadius->SetUnitCandidates("m cm");

  red_directory = new G4UIdirectory("/ds/detector/red/");
  fReDConfigurationCmd = new G4UIcmdWithAnInteger("/ds/detector/red/configuration", this);
  fReDConfigurationCmd->SetGuidance("ReD geometry control");

  fReDNDSurveyModeCmd = new G4UIcmdWithAString("/ds/detector/red/ND_survey_mode", this);
  G4String SurveyCandidate = "xyz spherical";
  fReDNDSurveyModeCmd->SetCandidates(SurveyCandidate);

  fReDNDSurveyFileNameCmd = new G4UIcmdWithAString("/ds/detector/red/ND_survey_filename", this);

  fReDNDPositioningModeCmd = new G4UIcmdWithAString("/ds/detector/red/ND_positioning_mode", this);
  G4String PositioningCandidate = "std kinematic";
  fReDNDPositioningModeCmd->SetCandidates(PositioningCandidate);

  fTargetTPCDistanceCmd = new G4UIcmdWithADoubleAndUnit("/ds/detector/red/TPCPosition", this);
  fTargetTPCDistanceCmd->SetGuidance("Distance between target and TPC's center.");
  fTargetTPCDistanceCmd->SetUnitCandidates("cm");
  fTargetTPCDistanceCmd->AvailableForStates(G4State_PreInit);

  fangleThetaXCmd = new G4UIcmdWithADoubleAndUnit("/ds/detector/red/TPCThetaX", this);
  fangleThetaXCmd->SetGuidance(
      "Angle between the HORIZONTAL plane (in the laboratory) containing the "
      "beam and the TPC's vector position (from "
      "target to center TPC).");
  fangleThetaXCmd->SetGuidance("Positive angles mean that the TPC is *below* the beam");
  fangleThetaXCmd->SetUnitCandidates("deg rad");
  fangleThetaXCmd->AvailableForStates(G4State_PreInit);

  fanglePhiXCmd = new G4UIcmdWithADoubleAndUnit("/ds/detector/red/TPCPhiX", this);
  fanglePhiXCmd->SetGuidance(
      "Angle between the VERTICAL plane (in the laboratory) containing the "
      "beam and the TPC's vector position (from "
      "target to center TPC).");
  fanglePhiXCmd->SetGuidance(
      "Positive angles mean that the TPC is *on the left* to the beam looking "
      "*forward* with respect to the target");
  fanglePhiXCmd->SetUnitCandidates("deg rad");
  fanglePhiXCmd->AvailableForStates(G4State_PreInit);

  fProtoProtoBottomSiPMCmd = new G4UIcmdWithABool("/ds/detector/protoproto/bottomsipm", this);

  gantry_directory = new G4UIdirectory("/ds/detector/gantry/");
  fDSGantryConfiguration = new G4UIcmdWithABool("/ds/detector/gantry/configuration", this);
  fDSGantryPipeID = new G4UIcmdWithADoubleAndUnit("/ds/detector/gantry/pipeID", this);
  fDSGantryPipeID->SetUnitCandidates("m cm mm");
  fDSGantryPipeWall = new G4UIcmdWithADoubleAndUnit("/ds/detector/gantry/pipeWall", this);
  fDSGantryPipeWall->SetUnitCandidates("m cm mm");
  fDSGantryPipeBendingR = new G4UIcmdWithADoubleAndUnit("/ds/detector/gantry/pipeBendingR", this);
  fDSGantryPipeBendingR->SetUnitCandidates("m cm mm");
  fDSGantryDistanceTPCside = new G4UIcmdWithADoubleAndUnit("/ds/detector/gantry/distanceTPCside", this);
  fDSGantryDistanceTPCside->SetUnitCandidates("m cm mm");
  fDSGantryDistanceTPCbottom = new G4UIcmdWithADoubleAndUnit("/ds/detector/gantry/distanceTPCbottom", this);
  fDSGantryDistanceTPCbottom->SetUnitCandidates("m cm mm");
  fDSGantryDistanceBtwPipes = new G4UIcmdWithADoubleAndUnit("/ds/detector/gantry/distanceBtwPipes", this);
  fDSGantryDistanceBtwPipes->SetUnitCandidates("m cm mm");
  fDSGantrySurface = new G4UIcmdWithABool("/ds/detector/gantry/surface", this);
  fDSGantrySurfaceOption = new G4UIcmdWithAnInteger("/ds/detector/gantry/surfaceOption", this);
  fDSGantryShield = new G4UIcmdWithABool("/ds/detector/gantry/shield", this);
  fDSGantryShieldThickness = new G4UIcmdWithADoubleAndUnit("/ds/detector/gantry/shieldThickness", this);
  fDSGantryShieldThickness->SetUnitCandidates("cm mm");
  fDSGantryShieldHeight = new G4UIcmdWithADoubleAndUnit("/ds/detector/gantry/shieldHeight", this);
  fDSGantryShieldHeight->SetUnitCandidates("cm mm");
  fDSGantryShieldOffset = new G4UIcmdWithADoubleAndUnit("/ds/detector/gantry/shieldOffset", this);
  fDSGantryShieldOffset->SetUnitCandidates("cm mm");
  fDSGantryShieldMaterial = new G4UIcmdWithAnInteger("/ds/detector/gantry/shieldMaterial", this);

  fRefCmd = new G4UIcmdWithADouble("/ds/detector/ds20k_veto_esr_reflectivity", this);
  fRefCmd->SetGuidance("Reflectivity of ESR");
  fRefCmd->AvailableForStates(G4State_PreInit);

  fDS20kSiPMsCmd = new G4UIcmdWithABool("/ds/detector/ds20k_SiPMs", this);
  fDS20kSiPMsCmd->SetGuidance("Build and place SiPMs");
  fDS20kSiPMsCmd->AvailableForStates(G4State_PreInit);

  fDS20knSiPMsCmd = new G4UIcmdWithADouble("/ds/detector/ds20k_SiPMsNumber", this);
  fDS20knSiPMsCmd->SetGuidance("Required number of SiPMs");
  fDS20knSiPMsCmd->AvailableForStates(G4State_PreInit);

  fDS20kTPBOnSiPMCmd = new G4UIcmdWithABool("/ds/detector/ds20k_TPBOnSiPM", this);
  fDS20kTPBOnSiPMCmd->SetGuidance("TPB on SiPMs");
  fDS20kTPBOnSiPMCmd->AvailableForStates(G4State_PreInit);

  fDS20kSiPMUniformInsideSidesCmd = new G4UIcmdWithABool("/ds/detector/ds20k_SiPMUniformInsideSides", this);
  fDS20kSiPMUniformInsideSidesCmd->SetGuidance("SiPMs uniform configuration on IAB sides");
  fDS20kSiPMUniformInsideSidesCmd->AvailableForStates(G4State_PreInit);

  fDS20kSiPMUniformInsideCapsCmd = new G4UIcmdWithABool("/ds/detector/ds20k_SiPMUniformInsideCaps", this);
  fDS20kSiPMUniformInsideCapsCmd->SetGuidance("SiPMs uniform configuration on IAB caps");
  fDS20kSiPMUniformInsideCapsCmd->AvailableForStates(G4State_PreInit);

  fDS20kSiPMUniformOutsideSidesCmd = new G4UIcmdWithABool("/ds/detector/ds20k_SiPMUniformOutsideSides", this);
  fDS20kSiPMUniformOutsideSidesCmd->SetGuidance("SiPMs uniform configuration on OAB sides");
  fDS20kSiPMUniformOutsideSidesCmd->AvailableForStates(G4State_PreInit);

  fDS20kSiPMUniformOutsideCapsCmd = new G4UIcmdWithABool("/ds/detector/ds20k_SiPMUniformOutsideCaps", this);
  fDS20kSiPMUniformOutsideCapsCmd->SetGuidance("SiPMs uniform configuration on OAB caps");
  fDS20kSiPMUniformOutsideCapsCmd->AvailableForStates(G4State_PreInit);

  fDS20kbestLightCollectionCmd = new G4UIcmdWithABool("/ds/detector/ds20k_bestLightCollection", this);
  fDS20kbestLightCollectionCmd->SetGuidance("Best light collection configuration for SiPMs");
  fDS20kbestLightCollectionCmd->AvailableForStates(G4State_PreInit);

  fDS20kSiPMsAutoplacementCmd = new G4UIcmdWithABool("/ds/detector/ds20k_SiPMsAutoplacement", this);
  fDS20kSiPMsAutoplacementCmd->SetGuidance("SiPMs auto-placement from file");
  fDS20kSiPMsAutoplacementCmd->AvailableForStates(G4State_PreInit);

  fDS20kSiPMsAutoplacementFilenameCmd = new G4UIcmdWithAString("/ds/detector/ds20k_SiPMsAutoplacementFilename", this);
  fDS20kSiPMsAutoplacementFilenameCmd->SetGuidance("SiPMs auto-placement filename");
  fDS20kSiPMsAutoplacementFilenameCmd->AvailableForStates(G4State_PreInit);

  fDS20kCryoRCmd = new G4UIcmdWithADoubleAndUnit("/ds/detector/ds20cryo_cryostatR", this);
  fDS20kCryoHCmd = new G4UIcmdWithADoubleAndUnit("/ds/detector/ds20cryo_cryostatH", this);
  fDS20kCryoVacuumThCmd = new G4UIcmdWithADoubleAndUnit("/ds/detector/ds20cryo_vacuumTh", this);
  fDS20kCryoBottomCapCmd = new G4UIcmdWithADoubleAndUnit("/ds/detector/ds20cryo_bottomCapH", this);
  fDS20kCryoTopCapCmd = new G4UIcmdWithADoubleAndUnit("/ds/detector/ds20cryo_topCapH", this);
  fDS20kCryoTopOffsetCmd = new G4UIcmdWithADoubleAndUnit("/ds/detector/ds20cryo_topOffset", this);
  fDS20kCryoBottomOffsetCmd = new G4UIcmdWithADoubleAndUnit("/ds/detector/ds20cryo_bottomOffset", this);

  fDS20kGdConcentrationCmd = new G4UIcmdWithADouble("/ds/detector/ds20k_GdConcentration", this);

  fTPCConfigurationCmd = new G4UIcmdWithAString("/ds/detector/TPCconfiguration", this);
  G4String candidate_TPC = "GdPMMA PurePMMA Hybrid";
  fTPCConfigurationCmd->SetCandidates(candidate_TPC);

}

DSDetectorConstructionMessenger::~DSDetectorConstructionMessenger() {
  delete fDirectory;
  delete fConfigurationCmd;
  delete fDSGantryPipeID;
  delete fDSGantryPipeWall;
  delete fDSGantryPipeBendingR;
  delete fDSGantryDistanceTPCside;
  delete fDSGantryDistanceTPCbottom;
  delete fDSGantryDistanceBtwPipes;
  delete fDSGantrySurface;
  delete fDSGantrySurfaceOption;
  delete fDSGantryShield;
  delete fDSGantryShieldThickness;
  delete fDSGantryShieldHeight;
  delete fDSGantryShieldOffset;
  delete fDSGantryShieldMaterial;

  delete fScintillatorCmd;
  delete fSourceCmd;
  delete fArDMLeadHeightCmd;
  delete fExtLArScintillatingCmd;
  delete fCopperRingsCmd;
  delete fVetoYieldFactorCmd;
  delete fCTFfillCmd;
  delete fDS20CryoMaterialCmd;
  delete fDS20WindowMaterialCmd;
  delete fDS20kCryoWallThickCmd;
  delete fDS20kCryoCornerDistanceCmd;
  delete fDS20kTPCheightCmd;
  delete fDS20kTPCedgeCmd;
  delete fDS20kTPCVesselThicknessCmd;
  delete fDS20kTPBThicknessCmd;
  delete fDS20kHDPEShellThicknessCmd;
  delete fDS20kHDPEShellCapThicknessCmd;
  delete fDS20kTPBLayersCmd;

  delete fGdLayerThicknessCmd;
  delete fDS20kfLArThicknessAboveTPCCmd;
  delete fAcrylicVesselDiameterCmd;
  delete fIDiameterCmd;
  delete fWTankRCmd;
  delete fWTankHeightCmd;
  delete fIsDS20kLSVCmd;
  delete fLicNDRadiusCmd;
  delete fLicNDThetaCmd;
  delete fLicNDPhi1Cmd;
  delete fLicNDPhi2Cmd;
  delete fLicNDNuclRecEnergyCmd;
  delete fLicNeutronEnergyCmd;

  delete fLicActivateWallCmd;
  delete fLicWallThickCmd;
  delete fLicWallWindowLengthCmd;
  delete fLicWallDistanceToTPCCmd;

  delete fLicExpHallCmd;

  delete fAriserSourceDistanceToCryoCmd;
  delete fAriserBEGeDistanceToCryoCmd;
  delete fAriserBEGeAngleCmd;

  delete fSourceHolderCenterCmd;
  delete fSourceHolderThetaCmd;
  delete fSourceHolderPhiCmd;
  delete fSourceHolderLeadCmd;

  delete dart_directory;
  delete fDSDartTeflonThickness;
  delete fDSDartTeflonHeight;
  delete fDSDartTeflonRadius;

  delete red_directory;
  delete fReDConfigurationCmd;
  delete fReDNDSurveyModeCmd;
  delete fReDNDSurveyFileNameCmd;
  delete fReDNDPositioningModeCmd;
  delete fTargetTPCDistanceCmd;
  delete fangleThetaXCmd;
  delete fanglePhiXCmd;

  delete fProtoProtoBottomSiPMCmd;

  delete gantry_directory;
  delete fDSGantryConfiguration;

  delete fRefCmd;
  delete fDS20kSiPMsCmd;
  delete fDS20knSiPMsCmd;
  delete fDS20kTPBOnSiPMCmd;
  delete fDS20kbestLightCollectionCmd;
  delete fDS20kSiPMUniformInsideSidesCmd;
  delete fDS20kSiPMUniformInsideCapsCmd;
  delete fDS20kSiPMUniformOutsideSidesCmd;
  delete fDS20kSiPMUniformOutsideCapsCmd;
  delete fDS20ReflectorCmd;

  delete fDS20kbestLightCollectionCmd;

  delete fDS20kSiPMsAutoplacementCmd;
  delete fDS20kSiPMsAutoplacementFilenameCmd;

  // cryostat plan b
  delete fDS20kCryoRCmd;
  delete fDS20kCryoHCmd;
  delete fDS20kCryoVacuumThCmd;
  delete fDS20kCryoBottomCapCmd;
  delete fDS20kCryoTopCapCmd;
  delete fDS20kCryoTopOffsetCmd;
  delete fDS20kCryoBottomOffsetCmd;

  delete fDS20kGdConcentrationCmd;

  delete fTPCConfigurationCmd;
}

void DSDetectorConstructionMessenger::SetNewValue(G4UIcommand* command, G4String newValue) {

  if (command == fConfigurationCmd) {
    DSLog(routine) << "Detector Configuration: " << newValue << endlog;
    if (newValue == "ds50_full") {
      fDetectorConstruction->SetDetectorConfiguration(0);
    } else if (newValue == "ds50_tpc_veto") {
      fDetectorConstruction->SetDetectorConfiguration(1) ;
    } else if (newValue == "ds50_tpc") {
      fDetectorConstruction->SetDetectorConfiguration(2 );
    } else if (newValue == "tester") {
      fDetectorConstruction->SetDetectorConfiguration(3);
    } else if (newValue == "ds20k") {
      fDetectorConstruction->SetDetectorConfiguration(4 );
      DSStorage::Get()->Set20KGeometry(1);
    } else if (newValue == "ds20k_tpc") {
      fDetectorConstruction->SetDetectorConfiguration(5 );
      DSStorage::Get()->Set20KGeometry(1);
    } else if (newValue == "dsproto") {
      fDetectorConstruction->SetDetectorConfiguration(6 );
    } else if (newValue == "dsproto0") {
      fDetectorConstruction->SetDetectorConfiguration(7 );
    } else if (newValue == "red_tpc") {
      fDetectorConstruction->SetDetectorConfiguration(8 );
    } else if (newValue == "dart") {
      fDetectorConstruction->SetDetectorConfiguration(9 );
    } else if (newValue == "ardm") {
      fDetectorConstruction->SetDetectorConfiguration(10 );
    } else if (newValue == "licorne") {
      fDetectorConstruction->SetDetectorConfiguration(11);
    } else if (newValue == "ariser") {
      fDetectorConstruction->SetDetectorConfiguration(12);
    } else {
      DSLog(fatal) << "Fatal: detector geometry not recognized. Please check the spelling" << endlog;
    }

  } else if (command == fScintillatorCmd) {
    DSLog(routine) << "Choosen Scintillator: " << newValue << endlog;
    string volOrmass = " mass";
    if (newValue == "BoronScintillator") {
      DSStorage::Get()->SetScintillator(0);
      volOrmass = " volume";
    }
    if (newValue == "GdScintillator") DSStorage::Get()->SetScintillator(1);
    if (newValue == "Li6Scintillator") DSStorage::Get()->SetScintillator(2);
    if (newValue == "NatLi6Scintillator") DSStorage::Get()->SetScintillator(3);
    if (newValue == "OrtoCarbScintillator") DSStorage::Get()->SetScintillator(4);
    DSLog(routine) << newValue << " with " << DSStorage::Get()->GetTMBfraction() << volOrmass << " fraction" << endlog;
    DSLog(routine) << "...you can specify a new concentration with "
                      "/ds/manager/scint_mass_fraction command "
                   << endlog;
  } else if (command == fCTFfillCmd) {
    // G4String candidates = "BoronScintillator WaterGd Water";
    DSLog(routine) << "Choosen Material inside CTF: " << newValue << endlog;
    if (newValue == "Air") DSStorage::Get()->SetCTFfill(3);
    if (newValue == "BoronScintillator") DSStorage::Get()->SetCTFfill(2);
    if (newValue == "GdWater") DSStorage::Get()->SetCTFfill(1);
    if (newValue == "Water") DSStorage::Get()->SetCTFfill(0);
    /*
  }else  if (command == fSourceCmd) {
    fDetectorConstruction->SetIsSource(true);
    DSStorage::Get()->SetSourcePosition(fSourceCmd->ConvertToDimensioned3Vector(newValue));
    DSLog(routine) << "Source position: " << newValue << endlog ;
    */
  } else if (command == fDS20ReflectorCmd) {
    if (newValue == "baseline") {
      DSStorage::Get()->SetDS20kReflectorAlternate(0);
    } else if (newValue == "alternate") {
      DSStorage::Get()->SetDS20kReflectorAlternate(1);
    }
    DSLog(routine) << "DS20k TPC reflector configuration: " << newValue << endlog;
  } else if (command == fExtLArScintillatingCmd) {
    DSStorage::Get()->SetIsExternalLArScintillating(fExtLArScintillatingCmd->ConvertToBool(newValue));
    DSLog(routine) << "External Liquid Argon Scintillation: " << newValue << endlog;
  } else if (command == fCopperRingsCmd) {
    DSStorage::Get()->SetIsCopperRings(fCopperRingsCmd->ConvertToBool(newValue));
    DSLog(routine) << "True Copper Field Rings: " << newValue << endlog;
  } else if (command == fVetoYieldFactorCmd) {
    DSStorage::Get()->SetVetoYieldFactor(fVetoYieldFactorCmd->ConvertToDouble(newValue));
    DSLog(routine) << "Veto Yield Factor: " << newValue << endlog;
  } else if (command == fSourceHolderCenterCmd) {
    DSStorage::Get()->SetSourceHolderCenter(fSourceHolderCenterCmd->ConvertToDimensioned3Vector(newValue));
    DSStorage::Get()->SetSourceHolderFlag(true);
    DSLog(routine) << "Holder Radius: " << newValue << endlog;
  } else if (command == fSourceHolderThetaCmd) {
    DSStorage::Get()->SetSourceHolderTheta(fSourceHolderThetaCmd->ConvertToDimensionedDouble(newValue));
    DSStorage::Get()->SetSourceHolderFlag(true);
    DSLog(routine) << "Holder rotation Theta: " << newValue << endlog;
  } else if (command == fSourceHolderPhiCmd) {
    DSStorage::Get()->SetSourceHolderPhi(fSourceHolderPhiCmd->ConvertToDimensionedDouble(newValue));
    DSStorage::Get()->SetSourceHolderFlag(true);
    DSLog(routine) << "Holder rotation Phi: " << newValue << endlog;
  } else if (command == fSourceHolderLeadCmd) {
    DSStorage::Get()->SetSourceHolderLeadFlag(fSourceHolderLeadCmd->ConvertToBool(newValue));
    DSStorage::Get()->SetSourceHolderFlag(true);
  } else if (command == fArDMLeadHeightCmd) {
    DSStorage::Get()->SetArDMLeadHeight(fArDMLeadHeightCmd->ConvertToDimensionedDouble(newValue));
    DSLog(routine) << "ArDM lead height: " << newValue << endlog;
  } else if (command == fDS20kCryoWallThickCmd) {
    DSStorage::Get()->SetDS20kCryoWallThick(fDS20kCryoWallThickCmd->ConvertToDimensionedDouble(newValue));
    DSStorage::Get()->SetIsCustomCryostat(true);
    DSLog(routine) << "Custom cryo, DS20k CryoWallThickness: " << newValue << endlog;
  } else if (command == fDS20CryoMaterialCmd) {
    DSStorage::Get()->SetDS20kCryoMaterial(fDS20CryoMaterialCmd->ConvertToInt(newValue));
    // DSStorage::Get()->SetIsCustomCryostat(true);
    DSLog(routine) << "DS20k cryo material: " << newValue << endlog;
  } else if (command == fDS20WindowMaterialCmd) {
    DSStorage::Get()->SetDS20kWindowMaterial(fDS20WindowMaterialCmd->ConvertToInt(newValue));
    DSLog(routine) << "DS20k window material: " << newValue << endlog;
  } else if (command == fDS20kCryoCornerDistanceCmd) {
    DSStorage::Get()->SetDS20kCryoCornerDistance(fDS20kCryoCornerDistanceCmd->ConvertToDimensionedDouble(newValue));
    DSStorage::Get()->SetIsCustomCryostat(true);
    DSLog(routine) << "Custom cryo, DS20k CryoCornerDistance: " << newValue << endlog;
  } else if (command == fDS20kTPCheightCmd) {
    DSStorage::Get()->SetDS20kTPCheight(fDS20kTPCheightCmd->ConvertToDimensionedDouble(newValue));
    DSStorage::Get()->SetIsCustomCryostat(true);
    DSLog(routine) << "Custom cryo, DS20k TPC height: " << newValue << endlog;
  } else if (command == fDS20kTPCedgeCmd) {
    DSStorage::Get()->SetDS20kTPCedge(fDS20kTPCedgeCmd->ConvertToDimensionedDouble(newValue));
    DSStorage::Get()->SetIsCustomCryostat(true);
    DSLog(routine) << "Custom cryo, DS20k TPC edge : " << newValue << endlog;
  } else if (command == fGdLayerThicknessCmd) {
    DSStorage::Get()->SetGdLayerThickness(fGdLayerThicknessCmd->ConvertToDimensionedDouble(newValue));
    DSStorage::Get()->SetIsCustomCryostat(true);
    DSLog(routine) << "Custom cryo, DS20k GdLayer Thickness: " << newValue << endlog;
  } else if (command == fDS20kfLArThicknessAboveTPCCmd) {
    DSStorage::Get()->SetDS20kLArThicknessAboveTPC(fDS20kfLArThicknessAboveTPCCmd->ConvertToDimensionedDouble(newValue));
    DSLog(routine) << "DS20k, amount of NSLAr above the TPC: " << newValue << endlog;
  } else if (command == fDS20kSiPMOffsetCmd) {
    DSStorage::Get()->SetSiPMOffset(fDS20kSiPMOffsetCmd->ConvertToDimensionedDouble(newValue));
    DSLog(routine) << "DS20k SiPM Offset: " << newValue << endlog;
  } else if (command == fIDiameterCmd) {
    DSStorage::Get()->SetIDiameter(fIDiameterCmd->ConvertToDimensionedDouble(newValue));
    DSLog(routine) << "DS20k - LSV sphere diameter: " << newValue << endlog;
  } else if (command == fAcrylicVesselDiameterCmd) {
    DSStorage::Get()->SetAcrylicVesselDiameter(fAcrylicVesselDiameterCmd->ConvertToDimensionedDouble(newValue));
    DSLog(routine) << "DS20k - Acrylic Vessel diameter: " << newValue << endlog;
  } else if (command == fDS20kLArBufferThicknessCmd) {
    DSStorage::Get()->SetDS20kLArBufferThickness(fDS20kLArBufferThicknessCmd->ConvertToDimensionedDouble(newValue));
    DSLog(routine) << "DS20k, LAr Buffer Thickness: " << newValue << endlog;
  } else if (command == fDS20kTPCVesselThicknessCmd) {
    DSStorage::Get()->SetDS20kTPCVesselThickness(fDS20kTPCVesselThicknessCmd->ConvertToDimensionedDouble(newValue));
    DSLog(routine) << "DS20k, TPC Acrylic Vessel thickness: " << newValue << endlog;
  } else if (command == fDS20kHDPEShellThicknessCmd) {
    DSStorage::Get()->SetDS20kHDPEShellThickness(fDS20kHDPEShellThicknessCmd->ConvertToDimensionedDouble(newValue));
    DSLog(routine) << "DS20k, HDPE Shell thickness: " << newValue << endlog;
  } else if (command == fDS20kHDPEShellCapThicknessCmd) {
    DSStorage::Get()->SetDS20kHDPEShellCapThickness(fDS20kHDPEShellCapThicknessCmd->ConvertToDimensionedDouble(newValue));
    DSLog(routine) << "DS20k, HDPE Shell Cap thickness: " << newValue << endlog;
  } else if (command == fDS20kGdPlasticThicknessCmd) {
    DSStorage::Get()->SetDS20kGdPlasticThickness(fDS20kGdPlasticThicknessCmd->ConvertToDimensionedDouble(newValue));
    DSLog(routine) << "DS20k, Plastic Shell Thickness: " << newValue << endlog;
  } else if (command == fDS20kTPBThicknessCmd) {
    DSStorage::Get()->SetDS20kTPBThickness(fDS20kTPBThicknessCmd->ConvertToDimensionedDouble(newValue));
    DSLog(routine) << "DS20k, TPB Thickness: " << newValue << endlog;
  } else if (command == fDS20kTPBLayersCmd) {
    DSStorage::Get()->SetDS20kTPBLayers(newValue);
    DSLog(routine) << "DS20k, TPB Layers: " << newValue << endlog;
  } else if (command == fWTankRCmd) {
    DSStorage::Get()->SetWTankR(fWTankRCmd->ConvertToDimensionedDouble(newValue));
    DSLog(routine) << "Water Tank (CTF) Radius : " << newValue << endlog;
  } else if (command == fWTankHeightCmd) {
    DSStorage::Get()->SetWTankHeight(fWTankHeightCmd->ConvertToDimensionedDouble(newValue));
    DSLog(routine) << "Water Tank (CTF) Height : " << newValue << endlog;
  } else if (command == fIsDS20kLSVCmd) {
    DSStorage::Get()->SetIsDS20kLSV(fIsDS20kLSVCmd->ConvertToBool(newValue));
    DSLog(routine) << "DS20k - LSV sphere built : " << newValue << endlog;
  } else if (command == fLicNDRadiusCmd) {
    DSStorage::Get()->SetLicNDRadius(fLicNDRadiusCmd->ConvertToDimensionedDouble(newValue));
    DSLog(routine) << "Licorne neutron detector radius : " << newValue << endlog;
  } else if (command == fLicNDThetaCmd) {
    DSStorage::Get()->SetLicNDTheta(fLicNDThetaCmd->ConvertToDimensionedDouble(newValue));
    DSLog(routine) << "Licorne neutron detector theta : " << newValue << endlog;
  } else if (command == fLicNDPhi1Cmd) {
    DSStorage::Get()->SetLicNDPhi1(fLicNDPhi1Cmd->ConvertToDimensionedDouble(newValue));
    DSLog(routine) << "Licorne neutron detector no1 phi : " << newValue << endlog;
  } else if (command == fLicNDPhi2Cmd) {
    DSStorage::Get()->SetLicNDPhi2(fLicNDPhi2Cmd->ConvertToDimensionedDouble(newValue));
    DSLog(routine) << "Licorne neutron detector no2 phi : " << newValue << endlog;
  } else if (command == fLicNDNuclRecEnergyCmd) {
    DSStorage::Get()->SetLicNDNuclRecE(fLicNDNuclRecEnergyCmd->ConvertToDimensionedDouble(newValue));
    DSLog(routine) << "Licorne neutron detector theta to get Er = " << newValue << endlog;
  } else if (command == fLicNeutronEnergyCmd) {
    DSStorage::Get()->SetLicNeutronE(fLicNeutronEnergyCmd->ConvertToDimensionedDouble(newValue));
    DSLog(routine) << "Licorne neutron energy : " << newValue << endlog;
  } else if (command == fLicActivateWallCmd) {
    DSStorage::Get()->SetLicWallActivation(fLicActivateWallCmd->ConvertToBool(newValue));
    DSLog(routine) << "Licorne wall activation : " << newValue << endlog;
  } else if (command == fLicWallThickCmd) {
    DSStorage::Get()->SetLicWallThick(fLicWallThickCmd->ConvertToDimensionedDouble(newValue));
    DSLog(routine) << "Licorne wall thick : " << newValue << endlog;
  } else if (command == fLicWallWindowLengthCmd) {
    DSStorage::Get()->SetLicWallWindowLength(fLicWallWindowLengthCmd->ConvertToDimensionedDouble(newValue));
    DSLog(routine) << "Licorne wall window length : " << newValue << endlog;
  } else if (command == fLicWallDistanceToTPCCmd) {
    DSStorage::Get()->SetLicWallDistanceToTPC(fLicWallDistanceToTPCCmd->ConvertToDimensionedDouble(newValue));
    DSLog(routine) << "Licorne wall distance to TPC : " << newValue << endlog;
  } else if (command == fLicExpHallCmd) {
    DSStorage::Get()->SetLicExpHall(fLicExpHallCmd->ConvertToBool(newValue));
    DSLog(routine) << "Licorne experimental hall : " << newValue << endlog;
  } else if (command == fAriserSourceDistanceToCryoCmd) {
    DSStorage::Get()->SetAriserSourceDistanceToCryo(fAriserSourceDistanceToCryoCmd->ConvertToDimensionedDouble(newValue));
    DSLog(routine) << "ARIS-ER source distance to cryostat : " << newValue << endlog;
  } else if (command == fAriserBEGeDistanceToCryoCmd) {
    DSStorage::Get()->SetAriserBEGeDistanceToCryo(fAriserBEGeDistanceToCryoCmd->ConvertToDimensionedDouble(newValue));
    DSLog(routine) << "ARIS-ER BEGe detector distance to cryostat : " << newValue << endlog;
  } else if (command == fAriserBEGeAngleCmd) {
    DSStorage::Get()->SetAriserBEGeAngle(fAriserBEGeAngleCmd->ConvertToDimensionedDouble(newValue));
    DSLog(routine) << "ARIS-ER BEGe detector angle : " << newValue << endlog;
  } else if (command == fDSDartTeflonThickness) {
    DSStorage::Get()->SetDartTeflonThickness(fDSDartTeflonThickness->ConvertToDimensionedDouble(newValue));
    DSLog(routine) << "Dart Teflon  thickness : " << newValue << endlog;
  } else if (command == fDSDartTeflonHeight) {
    DSStorage::Get()->SetDartTeflonHeight(fDSDartTeflonHeight->ConvertToDimensionedDouble(newValue));
    DSLog(routine) << "Dart Teflon  height : " << newValue << endlog;
  } else if (command == fDSDartTeflonRadius) {
    DSStorage::Get()->SetDartTeflonRadius(fDSDartTeflonRadius->ConvertToDimensionedDouble(newValue));
    DSLog(routine) << "Dart Teflon  radius : " << newValue << endlog;
  } else if (command == fDSDartDewarThickness) {
    DSStorage::Get()->SetDartDewarThickness(fDSDartDewarThickness->ConvertToDimensionedDouble(newValue));
    DSLog(routine) << "Dart Dewar  thickness : " << newValue << endlog;
  } else if (command == fDSDartDewarHeight) {
    DSStorage::Get()->SetDartDewarHeight(fDSDartDewarHeight->ConvertToDimensionedDouble(newValue));
    DSLog(routine) << "Dart Dewar  height : " << newValue << endlog;
  } else if (command == fDSDartDewarRadius) {
    DSStorage::Get()->SetDartDewarRadius(fDSDartDewarRadius->ConvertToDimensionedDouble(newValue));
    DSLog(routine) << "Dart Dewar  radius : " << newValue << endlog;
  } else if (command == fReDConfigurationCmd) {
    DSStorage::Get()->SetReDConfiguration(fReDConfigurationCmd->ConvertToInt(newValue));
    DSLog(routine) << "ReD Detector Configuration : " << newValue << endlog;
  } else if (command == fReDNDSurveyModeCmd) {
    DSLog(routine) << "Choosen ReD Neutron Detector Survey mode: " << newValue << endlog;
    if (newValue == "xyz") DSStorage::Get()->SetReDNDSurveyMode(DSStorage::xyz);
    else if (newValue == "spherical")
      DSStorage::Get()->SetReDNDSurveyMode(DSStorage::spherical);
  } else if (command == fReDNDSurveyFileNameCmd) {
    DSStorage::Get()->SetReDNDSurveyFileName(newValue);
    DSLog(routine) << "ReD Neutron detector survey read from the file: " << newValue << endlog;
  } else if (command == fReDNDPositioningModeCmd) {
    DSLog(routine) << "Choosen ReD Neutron Detector Positioning Mode: " << newValue << endlog;
    if (newValue == "std") DSStorage::Get()->SetReDNDPositioningMode(DSStorage::std);
    else if (newValue == "kinematic")
      DSStorage::Get()->SetReDNDPositioningMode(DSStorage::kinematic);
  } else if (command == fTargetTPCDistanceCmd) {
    DSStorage::Get()->SetTargetToTPCDistance(fTargetTPCDistanceCmd->GetNewDoubleValue(newValue));
    DSLog(routine) << "TPC to Beam Distance: " << newValue << endlog;
  } else if (command == fangleThetaXCmd) {
    DSStorage::Get()->SetBeamToTPCAngleThetaX(fangleThetaXCmd->GetNewDoubleValue(newValue));
    DSLog(routine) << "Beam to TPC Angle Theta: " << newValue << endlog;
  } else if (command == fanglePhiXCmd) {
    DSStorage::Get()->SetBeamToTPCAnglePhiX(fanglePhiXCmd->GetNewDoubleValue(newValue));
    DSLog(routine) << "Beam to TPC Angle Phi: " << newValue << endlog;
  } else if (command == fProtoProtoBottomSiPMCmd) {
    DSStorage::Get()->SetIsProtoProtoBottomSiPM(fProtoProtoBottomSiPMCmd->ConvertToBool(newValue));
    DSLog(routine) << "DSProtoProto Bottom SiPMs: " << newValue << endlog;
  } else if (command == fDSGantryConfiguration) {
    DSStorage::Get()->SetGantryConfiguration(fDSGantryConfiguration->ConvertToBool(newValue));
    DSLog(routine) << "Gantry Configuration  " << newValue << endlog;
  } else if (command == fDSGantryPipeID) {
    DSStorage::Get()->SetGantryPipeID(fDSGantryPipeID->GetNewDoubleValue(newValue));
    DSLog(routine) << "Gantry Pipe ID  " << newValue << endlog;
  } else if (command == fDSGantryPipeWall) {
    DSStorage::Get()->SetGantryPipeWall(fDSGantryPipeWall->GetNewDoubleValue(newValue));
    DSLog(routine) << "Gantry Pipe Wall  " << newValue << endlog;
  } else if (command == fDSGantryPipeBendingR) {
    DSStorage::Get()->SetGantryPipeBendingR(fDSGantryPipeBendingR->GetNewDoubleValue(newValue));
    DSLog(routine) << "Gantry Pipe Bending R  " << newValue << endlog;
  } else if (command == fDSGantryDistanceTPCside) {
    DSStorage::Get()->SetGantryDistanceTPCside(fDSGantryDistanceTPCside->GetNewDoubleValue(newValue));
    DSLog(routine) << "Gantry Pipe Distance TPC side  " << newValue << endlog;
  } else if (command == fDSGantryDistanceTPCbottom) {
    DSStorage::Get()->SetGantryDistanceTPCbottom(fDSGantryDistanceTPCbottom->GetNewDoubleValue(newValue));
    DSLog(routine) << "Gantry Pipe Distance TPC bottom  " << newValue << endlog;
  } else if (command == fDSGantryDistanceBtwPipes) {
    DSStorage::Get()->SetGantryDistanceBtwPipes(fDSGantryDistanceBtwPipes->GetNewDoubleValue(newValue));
    DSLog(routine) << "Gantry Pipe Distance Btw Pipes  " << newValue << endlog;
  } else if (command == fDSGantrySurface) {
    DSStorage::Get()->SetGantrySurface(fDSGantrySurface->ConvertToBool(newValue));
    DSLog(routine) << "Gantry Surface Build " << newValue << endlog;
  } else if (command == fDSGantrySurfaceOption) {
    DSStorage::Get()->SetGantrySurfaceOption(fDSGantrySurfaceOption->ConvertToInt(newValue));
    DSLog(routine) << "Gantry Surface Option " << newValue << endlog;
  } else if (command == fDSGantryShield) {
    DSStorage::Get()->SetGantryShield(fDSGantryShield->ConvertToBool(newValue));
    DSLog(routine) << "Gantry Shield Build " << newValue << endlog;
  } else if (command == fDSGantryShieldThickness) {
    DSStorage::Get()->SetGantryShieldThickness(fDSGantryShieldThickness->GetNewDoubleValue(newValue));
    DSLog(routine) << "Gantry Shield Thickness " << newValue << endlog;
  } else if (command == fDSGantryShieldHeight) {
    DSStorage::Get()->SetGantryShieldHeight(fDSGantryShieldHeight->GetNewDoubleValue(newValue));
    DSLog(routine) << "Gantry Shield Height " << newValue << endlog;
  } else if (command == fDSGantryShieldOffset) {
    DSStorage::Get()->SetGantryShieldOffset(fDSGantryShieldOffset->GetNewDoubleValue(newValue));
    DSLog(routine) << "Gantry Shield Offset " << newValue << endlog;
  } else if (command == fDSGantryShieldMaterial) {
    DSStorage::Get()->SetGantryShieldMaterial(fDSGantryShieldMaterial->ConvertToInt(newValue));
    DSLog(routine) << "Gantry Shield Material " << newValue << endlog;
  } else if (command == fRefCmd) {
    DSParameters::Get()->SetTeflonTPBVisRef(fRefCmd->ConvertToDouble(newValue));
    DSLog(routine) << "ESR reflectivity of the ESR set to: " << newValue << endlog;
  } else if (command == fDS20kSiPMsCmd) {
    DSStorage::Get()->SetDS20kSiPMs(fDS20kSiPMsCmd->ConvertToBool(newValue));
    DSLog(routine) << "Build and place SiPMs: " << newValue << endlog;

  } else if (command == fDS20knSiPMsCmd) {
    DSStorage::Get()->SetDS20knSiPMs(fDS20knSiPMsCmd->ConvertToDouble(newValue));
    DSLog(routine) << "Required number of SiPMs: " << newValue << endlog;

  } else if (command == fDS20kTPBOnSiPMCmd) {
    DSStorage::Get()->SetDS20kTPBOnSiPM(fDS20kTPBOnSiPMCmd->ConvertToBool(newValue));
    DSLog(routine) << "TPB on SiPMs: " << newValue << endlog;

  } else if (command == fDS20kbestLightCollectionCmd) {
    DSStorage::Get()->SetDS20kbestLightCollection(fDS20kbestLightCollectionCmd->ConvertToBool(newValue));
    DSLog(routine) << "Best light collection configuration for SiPMs: " << newValue << endlog;

  } else if (command == fDS20kSiPMUniformInsideSidesCmd) {
    DSStorage::Get()->SetDS20kSiPMUniformInsideSides(fDS20kSiPMUniformInsideSidesCmd->ConvertToBool(newValue));
    DSLog(routine) << "SiPMs uniform configuration on IAB sides: " << newValue << endlog;
  } else if (command == fDS20kSiPMUniformInsideCapsCmd) {
    DSStorage::Get()->SetDS20kSiPMUniformInsideCaps(fDS20kSiPMUniformInsideCapsCmd->ConvertToBool(newValue));
    DSLog(routine) << "SiPMs uniform configuration on IAB caps: " << newValue << endlog;
  } else if (command == fDS20kSiPMUniformOutsideSidesCmd) {
    DSStorage::Get()->SetDS20kSiPMUniformOutsideSides(fDS20kSiPMUniformOutsideSidesCmd->ConvertToBool(newValue));
    DSLog(routine) << "SiPMs uniform configuration on OAB sides: " << newValue << endlog;
  } else if (command == fDS20kSiPMUniformOutsideCapsCmd) {
    DSStorage::Get()->SetDS20kSiPMUniformOutsideCaps(fDS20kSiPMUniformOutsideCapsCmd->ConvertToBool(newValue));
    DSLog(routine) << "SiPMs uniform configuration on OAB caps: " << newValue << endlog;
  }

  else if (command == fDS20kSiPMsAutoplacementCmd) {
    DSStorage::Get()->SetDS20kSiPMsAutoplacement(fDS20kSiPMsAutoplacementCmd->ConvertToBool(newValue));
    DSLog(routine) << "SiPMs auto-placement from file: " << newValue << endlog;
  }

  else if (command == fDS20kSiPMsAutoplacementFilenameCmd) {
    DSStorage::Get()->SetDS20kSiPMsAutoplacementFilename(fDS20kSiPMsAutoplacementFilenameCmd->ConvertToString(newValue));
    DSLog(routine) << "SiPMs auto-placement filename: " << newValue << endlog;
  } else if (command == fDS20kCryoRCmd) {
    DSStorage::Get()->SetDS20kCryoR(fDS20kCryoRCmd->ConvertToDouble(newValue));
    DSLog(routine) << "DS20kCryoR set to: " << newValue << endlog;
  } else if (command == fDS20kCryoHCmd) {
    DSStorage::Get()->SetDS20kCryoH(fDS20kCryoHCmd->ConvertToDouble(newValue));
    DSStorage::Get()->SetDS20KCryoBarrelH(fDS20kCryoHCmd->ConvertToDouble(newValue));
    DSLog(routine) << "DS20kCryoH set to: " << newValue << endlog;
  } else if (command == fDS20kCryoVacuumThCmd) {
    DSStorage::Get()->SetDS20kCryoVacuumTh(fDS20kCryoVacuumThCmd->ConvertToDouble(newValue));
    DSLog(routine) << "DS20kCryoVacuumTh set to: " << newValue << endlog;
  } else if (command == fDS20kCryoBottomCapCmd) {
    DSStorage::Get()->SetDS20kCryoBottomCap(fDS20kCryoBottomCapCmd->ConvertToDouble(newValue));
    DSLog(routine) << "DS20kCryoBottomCap set to: " << newValue << endlog;
  } else if (command == fDS20kCryoTopCapCmd) {
    DSStorage::Get()->SetDS20kCryoTopCap(fDS20kCryoTopCapCmd->ConvertToDouble(newValue));
    DSLog(routine) << "DS20kCryoTopCap set to: " << newValue << endlog;
  } else if (command == fDS20kCryoTopOffsetCmd) {
    DSStorage::Get()->SetDS20kCryoTopOffset(fDS20kCryoTopOffsetCmd->ConvertToDouble(newValue));
    DSLog(routine) << "DS20kCryoTopOffset set to: " << newValue << endlog;
  } else if (command == fDS20kCryoBottomOffsetCmd) {
    DSStorage::Get()->SetDS20kCryoBottomOffset(fDS20kCryoBottomOffsetCmd->ConvertToDouble(newValue));
    DSLog(routine) << "DS20kCryoTopOffset set to: " << newValue << endlog;
  } else if (command == fDS20kGdConcentrationCmd) {
    DSStorage::Get()->SetDS20kGdConcentration(fDS20kGdConcentrationCmd->ConvertToDouble(newValue));
    DSLog(routine) << "Gd concentration set to: " << newValue << endlog;
  } else if (command == fTPCConfigurationCmd) {
    DSLog(routine) << "TPC Configuration: " << newValue << endlog;
       if (newValue == "GdPMMA") { 
    }  else if (newValue == "PurePMMA") {
      DSStorage::Get()->SetPurePMMATPC(true);
    } else if (newValue == "Hybrid") {
      DSStorage::Get()->SetHybridTPC(true);
    } else {
      DSLog(fatal) << "Fatal: TPC configuration not recognized. Please check the spelling" << endlog;
    }
  
  }
}
