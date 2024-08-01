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
//

#ifndef DSDetectorConstructionMessenger_h
#define DSDetectorConstructionMessenger_h 1

#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UImessenger.hh"
#include "globals.hh"
using namespace std;

class G4UIcmdWith3VectorAndUnit;
class DSDetectorConstruction;
class G4UIcmdWithADoubleAndUnit;
class DSDetectorConstructionMessenger : public G4UImessenger {

 public:
  DSDetectorConstructionMessenger(DSDetectorConstruction*);
  ~DSDetectorConstructionMessenger();

  void SetNewValue(G4UIcommand*, G4String);

 private:
  DSDetectorConstruction* fDetectorConstruction;
  G4UIdirectory* fDirectory;
  G4UIcmdWithAString* fConfigurationCmd;
  G4UIcmdWithAString* fScintillatorCmd;
  G4UIcmdWithAString* fCTFfillCmd;
  G4UIcmdWith3VectorAndUnit* fSourceCmd;
  G4UIcmdWithABool* fExtLArScintillatingCmd;
  G4UIcmdWithAString* fDS20ReflectorCmd;
  G4UIcmdWithABool* fCopperRingsCmd;
  G4UIcmdWithADouble* fVetoYieldFactorCmd;

  G4UIcmdWithABool* fIsDS20kLSVCmd;
  G4UIcmdWithAnInteger* fDS20CryoMaterialCmd;
  G4UIcmdWithAnInteger* fDS20WindowMaterialCmd;
  G4UIcmdWithADoubleAndUnit* fDS20kfLArThicknessAboveTPCCmd;
  G4UIcmdWithADoubleAndUnit* fDS20kCryoWallThickCmd;
  G4UIcmdWithADoubleAndUnit* fDS20kCryoCornerDistanceCmd;
  G4UIcmdWithADoubleAndUnit* fDS20kTPCheightCmd;
  G4UIcmdWithADoubleAndUnit* fDS20kTPCedgeCmd;
  G4UIcmdWithADoubleAndUnit* fDS20kSiPMOffsetCmd;
  G4UIcmdWithADoubleAndUnit* fDS20kTPCVesselThicknessCmd;
  G4UIcmdWithADoubleAndUnit* fDS20kHDPEShellThicknessCmd;
  G4UIcmdWithADoubleAndUnit* fDS20kHDPEShellCapThicknessCmd;
  G4UIcmdWithADoubleAndUnit* fGdLayerThicknessCmd;
  G4UIcmdWithADoubleAndUnit* fDS20kGdPlasticThicknessCmd;
  G4UIcmdWithADoubleAndUnit* fDS20kLArBufferThicknessCmd;
  G4UIcmdWithADoubleAndUnit* fDS20kTPBThicknessCmd;
  G4UIcmdWithAString* fDS20kTPBLayersCmd;
  G4UIcmdWithADoubleAndUnit* fArDMLeadHeightCmd;

  G4UIcmdWithADoubleAndUnit* fIDiameterCmd;
  G4UIcmdWithADoubleAndUnit* fAcrylicVesselDiameterCmd;
  G4UIcmdWithADoubleAndUnit* fWTankRCmd;
  G4UIcmdWithADoubleAndUnit* fWTankHeightCmd;

  G4UIcmdWithADoubleAndUnit* fLicNDRadiusCmd;
  G4UIcmdWithADoubleAndUnit* fLicNDThetaCmd;
  G4UIcmdWithADoubleAndUnit* fLicNDPhi1Cmd;
  G4UIcmdWithADoubleAndUnit* fLicNDPhi2Cmd;
  G4UIcmdWithADoubleAndUnit* fLicNDNuclRecEnergyCmd;

  G4UIcmdWithADoubleAndUnit* fLicNeutronEnergyCmd;

  G4UIcmdWithABool* fLicActivateWallCmd;
  G4UIcmdWithADoubleAndUnit* fLicWallThickCmd;
  G4UIcmdWithADoubleAndUnit* fLicWallWindowLengthCmd;
  G4UIcmdWithADoubleAndUnit* fLicWallDistanceToTPCCmd;

  G4UIcmdWithABool* fLicExpHallCmd;

  G4UIcmdWithADoubleAndUnit* fAriserSourceDistanceToCryoCmd;
  G4UIcmdWithADoubleAndUnit* fAriserBEGeDistanceToCryoCmd;
  G4UIcmdWithADoubleAndUnit* fAriserBEGeAngleCmd;

  G4UIcmdWith3VectorAndUnit* fSourceHolderCenterCmd;
  G4UIcmdWithADoubleAndUnit* fSourceHolderThetaCmd;
  G4UIcmdWithADoubleAndUnit* fSourceHolderPhiCmd;
  G4UIcmdWithABool* fSourceHolderLeadCmd;

  G4UIdirectory* dart_directory;
  G4UIcmdWithADoubleAndUnit* fDSDartTeflonThickness;
  G4UIcmdWithADoubleAndUnit* fDSDartTeflonHeight;
  G4UIcmdWithADoubleAndUnit* fDSDartTeflonRadius;
  G4UIcmdWithADoubleAndUnit* fDSDartDewarThickness;
  G4UIcmdWithADoubleAndUnit* fDSDartDewarHeight;
  G4UIcmdWithADoubleAndUnit* fDSDartDewarRadius;

  G4UIdirectory* red_directory;
  G4UIcmdWithAnInteger* fReDConfigurationCmd;

  G4UIcmdWithAString* fReDNDSurveyModeCmd;
  G4UIcmdWithAString* fReDNDSurveyFileNameCmd;
  G4UIcmdWithAString* fReDNDPositioningModeCmd;

  G4UIcmdWithADoubleAndUnit* fTargetTPCDistanceCmd;
  G4UIcmdWithADoubleAndUnit* fangleThetaXCmd;
  G4UIcmdWithADoubleAndUnit* fanglePhiXCmd;

  G4UIcmdWithABool* fProtoProtoBottomSiPMCmd;

  G4UIdirectory* gantry_directory;
  G4UIcmdWithABool* fDSGantryConfiguration;
  G4UIcmdWithADoubleAndUnit* fDSGantryPipeID;
  G4UIcmdWithADoubleAndUnit* fDSGantryPipeWall;
  G4UIcmdWithADoubleAndUnit* fDSGantryPipeBendingR;
  G4UIcmdWithADoubleAndUnit* fDSGantryDistanceTPCside;
  G4UIcmdWithADoubleAndUnit* fDSGantryDistanceTPCbottom;
  G4UIcmdWithADoubleAndUnit* fDSGantryDistanceBtwPipes;
  G4UIcmdWithABool* fDSGantrySurface;
  G4UIcmdWithAnInteger* fDSGantrySurfaceOption;
  G4UIcmdWithABool* fDSGantryShield;
  G4UIcmdWithADoubleAndUnit* fDSGantryShieldThickness;
  G4UIcmdWithADoubleAndUnit* fDSGantryShieldHeight;
  G4UIcmdWithADoubleAndUnit* fDSGantryShieldOffset;
  G4UIcmdWithAnInteger* fDSGantryShieldMaterial;

  G4UIcmdWithADouble* fRefCmd;

  G4UIcmdWithABool* fDS20kSiPMsCmd;

  G4UIcmdWithADouble* fDS20knSiPMsCmd;

  G4UIcmdWithABool* fDS20kTPBOnSiPMCmd;

  G4UIcmdWithABool* fDS20kSiPMUniformInsideSidesCmd;
  G4UIcmdWithABool* fDS20kSiPMUniformInsideCapsCmd;
  G4UIcmdWithABool* fDS20kSiPMUniformOutsideSidesCmd;
  G4UIcmdWithABool* fDS20kSiPMUniformOutsideCapsCmd;

  G4UIcmdWithABool* fDS20kbestLightCollectionCmd;

  G4UIcmdWithABool* fDS20kSiPMsAutoplacementCmd;
  G4UIcmdWithAString* fDS20kSiPMsAutoplacementFilenameCmd;

  G4UIcmdWithADoubleAndUnit* fDS20kCryoRCmd;
  G4UIcmdWithADoubleAndUnit* fDS20kCryoHCmd;
  G4UIcmdWithADoubleAndUnit* fDS20kCryoVacuumThCmd;
  G4UIcmdWithADoubleAndUnit* fDS20kCryoBottomCapCmd;
  G4UIcmdWithADoubleAndUnit* fDS20kCryoTopCapCmd;
  G4UIcmdWithADoubleAndUnit* fDS20kCryoTopOffsetCmd;
  G4UIcmdWithADoubleAndUnit* fDS20kCryoBottomOffsetCmd;

  G4UIcmdWithADouble* fDS20kGdConcentrationCmd;

  G4UIcmdWithAString* fTPCConfigurationCmd;
};
#endif
/*
 * $Log: DSDetectorConstructionMessenger.hh,v $
 * Revision 1.16  2016/03/08 21:10:36  swesterd
 * Added source holder
 *
 * Revision 1.15  2016/03/03 15:07:02  davini
 * DS-20k PMTs for Neutron Veto (20inches and 8inches types)
 *
 * Revision 1.14  2015/12/08 15:04:03  cgiganti
 * add experimental hall to licorne
 *
 * Revision 1.13  2015/12/07 09:49:36  riffard
 * Licorne : wall shape modification
 *
 * Revision 1.12  2015/12/03 15:59:41  riffard
 * Licorne : wall addition
 *
 * Revision 1.11  2015/11/26 13:30:17  dfranco
 * licorne
 *
 * Revision 1.10  2015/10/31 11:26:58  pagnes
 * added commands to customize WT, DS20k-cryostats and DS20k-LSV geometries
 *
 * Revision 1.9  2015/10/28 16:07:36  pagnes
 * DS20k --- new custom cryostat geometry (construction and messengers)
 *
 * Revision 1.8  2015/10/12 16:10:52  pagnes
 * command to change CTF fill /ds/detector/CTFfill added
 *
 * Revision 1.7  2015/10/06 14:49:45  pagnes
 * realistec field rings (first try) added with /ds/detector/true_field_rings
 * (default is without)
 *
 * Revision 1.6  2015/01/14 16:58:42  dfranco
 * added source vial and correspondent commands by Laura and Bianca; manual
 * updated
 *
 * Revision 1.5  2014/11/21 10:19:06  dfranco
 * added a command to scale the veto scintillation yield factor and fixed the
 * visible energy variable in the veto
 *
 * Revision 1.4  2014/11/20 15:32:12  dfranco
 * added a command to remove scintillation process from liquid argon between TPC
 * and cryostat
 *
 * Revision 1.3  2014/11/06 17:39:51  dfranco
 * added source and fixed a bug: removed gaseous nitrogen from the neutron veto
 *
 * Revision 1.2  2014/05/07 14:27:26  dfranco
 * fixed some bugs and added GdScintillator
 *
 * Revision 1.1  2014/05/07 12:20:51  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.2  2013/03/22 14:09:39  dfranco
 * Added the cvs logger code at the end of each file
 *
 */
