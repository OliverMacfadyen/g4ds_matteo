#include "DSDetectorARISER.hh"

#include <iostream>
#include <stdio.h>

#include "DSIO.hh"
#include "DSLogger.hh"
#include "DSMaterial.hh"
#include "DSParameters.hh"
#include "DSStorage.hh"
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4Orb.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalConstants.hh"
#include "G4SubtractionSolid.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4VisAttributes.hh"

using namespace std;

////////////////////////////////////////////////
//////        Detector Description      ////////
////////////////////////////////////////////////

// Based on the description of the Licorne detector

// 24th April 2015
// schematics of the TPC for the LICORNE test beam detector, starting from scaling the 5K geometry (the cryostats in particular) and from some drawings from Hanguo
// OPTICS is directly derived from DS50
//
// Placements are realized with no shifts. The correct position is calculated when making the solids.
//
// GENERAL SCHEME:
// - nested single volumes, from the outside to the center: Vacuum Jacket (steel), Vacuum, Inner Cryo (ssteel).
// - the fill of the inner Volume is divided in two region at z = 1200, which correspond to the LAr-GAr interface (Inactive LAr and InactiveGar)
//   - the TPC lays inside the InactiveLAr region (Copper Rings, Reflector, TPBSide , Active LAr)
//   - the Gas Pocket lays in the upper part (reflector, Active GAr)
// - other volumes:
//   - the grid is placed 0.5 cm below the LAr surface
//   - Fused silica windows with SiPm on the outfacing surface and TPB on the innerfacing one
//   - Pseudo Argon Layer (condensation on the top window)

// addition note on TPBSide and TPBTop
// TPBSide is a volume sorrunding the LAr Active volume on the side and bottom. TPBTop is only between
// LArLAyer (Pseudo Ar) and the top window

DSDetectorARISER::DSDetectorARISER(G4VPhysicalVolume *myMotherVolume) {
  fMotherVolume = myMotherVolume;
  DSLog(routine) << "Constructing ARIS-ER Setup Geometry" << endlog;
  G4bool myCheckOverlap = DSStorage::Get()->GetCheckOverlap();

  G4ThreeVector myZeros(0., 0., 0.);

  // ---------------------------------------------------//
  // ------ Cryostats from dewar_double.pdf -------------//
  // ---------------------------------------------------//

  G4double cryoOuterR = 6.10 * 2.54 / 2 * cm;
  G4double cryoInnerR = 5.58 * 2.54 / 2 * cm;
  G4double cryoOuterThickness = 0.12 * 2.54 * cm;
  G4double cryoInnerThickness = 0.13 * 2.54 * cm;
  G4double cryoInnerH = 23.13 * 2.54 / 2 * cm;
  G4double cryoOuterH = 26.13 * 2.54 / 2 * cm;

  DSLog(routine) << "cryo outer " << cryoOuterR << " Inner " << cryoInnerR << " thickness " << cryoOuterThickness << "(mm)" << endlog;
  DSLog(routine) << "cryo H (external) " << cryoOuterH << endlog;
  DSLog(routine) << "cryo H (internal) " << cryoInnerH << endlog;
  // -Sleeve+MainBody+Inner Sleeve  from PTFE_sleeve.pdf
  //--> added Top Cap+sleeve+bottom cap for height (3.26+0.55+0.565)//

  G4double SleeveOuterR = GetCm(3.752 / 2) * cm;  // To be checked
  G4double SleeveOuterH = GetCm(4.375 / 2) * cm;  // GetCm((2.94+0.28)/2)*cm;
  G4double SleeveInnerR = GetCm(3.00 / 2) * cm;
  G4double SleeveInnerH = GetCm(4.375 / 2) * cm;  // GetCm(2.94/2)*cm;

  DSLog(routine) << "Teflon outer " << SleeveOuterR << " Inner " << SleeveInnerR << "(mm)" << endlog;
  DSLog(routine) << "Teflon outer H " << SleeveOuterH << endlog;

  G4double TPCLArOuterR = SleeveInnerR;
  G4double MainBodyH = GetCm(3.00 / 2) * cm;

  G4double TPCGArOuterH = GetCm(0.26 / 2) * cm;

  G4double FusedSilicaOuterR = GetCm(3.50 / 2) * cm;
  G4double FusedSilicaThickness = 0.1 * cm;  // for now assume 0.1 cm
  G4double GridThickness = 0.001 * cm;       /// GetCm(0.1/2)*cm; //TBC

  G4double TPCLArOuterH = MainBodyH - GridThickness - TPCGArOuterH;

  DSLog(routine) << "TPC outer R " << SleeveInnerR << " LAr H " << TPCLArOuterH << " GAr H " << TPCGArOuterH << "(mm)" << endlog;

  G4ThreeVector LArTPCShift(0, 0, -GetCm(0.26) * cm);

  G4double PMTBottomR = 7.75 / 2 * cm;
  G4double PMTBottomWindowH = 0.1 * cm;

  G4double PMTBottomBody1R = PMTBottomR;
  G4double PMTBottomBody1H = (SleeveOuterH - TPCLArOuterH - FusedSilicaThickness * 2 - PMTBottomWindowH * 2 + LArTPCShift.getZ()) / 2;
  G4double PMTBottomBody1R_vac = PMTBottomBody1R - 0.2 * cm;

  DSLog(routine) << PMTBottomBody1H << " " << SleeveOuterH << " " << TPCLArOuterH << " " << FusedSilicaThickness << " " << PMTBottomWindowH << endlog;

  G4double PMTBottomBody2R = 5.33 / 2 * cm;
  G4double PMTBottomBody2R_Vac = PMTBottomBody2R - 0.2 * cm;
  G4double PMTBottomBody2H = 11.4 / 2 * cm - PMTBottomBody1H;

  G4double PMTTopBody1H = (SleeveInnerH - MainBodyH - FusedSilicaThickness * 2 + 0.32 * cm) / 2;  // 0.32 should be better understood!
  G4double PMTTopBody2H = 2.825 * cm - PMTTopBody1H;
  G4double PMTTopBodyL = 2.05 / 2 * cm;
  G4double PMTTopBodyL_Vac = PMTTopBodyL - 0.2 * cm;
  G4double PMTTopPhotoCathodeH = 0.1 * cm;
  G4double PMTTopBody1H_Vac = PMTTopBody1H - PMTTopPhotoCathodeH;

  DSLog(routine) << "TPC " << TPCLArOuterR << " " << TPCLArOuterH << endlog;
  DSLog(routine) << "Fused silica " << FusedSilicaThickness << endlog;

  // shift materials for Z>0
  G4ThreeVector GridShift(0, 0, TPCLArOuterH + LArTPCShift.getZ() + GridThickness);
  G4ThreeVector GasPocketShift(0, 0, GridShift.getZ() + GridThickness + TPCGArOuterH);
  G4ThreeVector AnodeShift(0, 0, GasPocketShift.getZ() + TPCGArOuterH + FusedSilicaThickness);
  G4ThreeVector PMTTopBody1Shift(0, 0, AnodeShift.getZ() + FusedSilicaThickness + PMTTopBody1H);
  G4ThreeVector PMTTopBody2Shift(0, 0, PMTTopBody1Shift.getZ() + PMTTopBody1H + PMTTopBody2H);

  // shift materials with Z<0
  G4ThreeVector CathodeShift(0, 0, LArTPCShift.getZ() - TPCLArOuterH - FusedSilicaThickness);
  G4ThreeVector PMTBottomWindowShift(0, 0, CathodeShift.getZ() - FusedSilicaThickness - PMTBottomWindowH);
  G4ThreeVector PMTBottomBody1Shift(0, 0, PMTBottomWindowShift.getZ() - PMTBottomWindowH - PMTBottomBody1H);
  G4ThreeVector PMTBottomBody2Shift(0, 0, PMTBottomBody1Shift.getZ() - PMTBottomBody1H - PMTBottomBody2H);

  DSStorage::Get()->SetLArGArBoundaryPosZ(GridShift.getZ() + 100 * um);

  DSLog(routine) << "cathode shift " << CathodeShift.getZ() << " " << FusedSilicaThickness << endlog;
  DSLog(routine) << "lar tpc shift " << LArTPCShift.getZ() << " " << TPCLArOuterH << endlog;
  DSLog(routine) << "grid shift " << GridShift.getZ() << " " << GridThickness << endlog;
  DSLog(routine) << "gas pocket shift " << GasPocketShift.getZ() << " " << TPCGArOuterH << endlog;
  DSLog(routine) << "anode shift " << AnodeShift.getZ() << endlog;

  DSLog(routine) << LArTPCShift.getZ() << " " << TPCLArOuterH << " " << GridThickness * 2 << " " << TPCGArOuterH * 2 << endlog;

  //-------------------------//
  //      Outer VacCryo      //
  //-------------------------//

  G4Tubs *fSolidCryoOuter = new G4Tubs("CryoOuter_Solid",
                                       0.,
                                       cryoOuterR,
                                       cryoOuterH,
                                       0.,
                                       2 * M_PI);

  G4LogicalVolume *fLogicLicorne = new G4LogicalVolume(fSolidCryoOuter,
                                                       DSMaterial::Get()->GetStainlessSteel(),
                                                       "OuterCryostat_Logic");

  fPhysicLicorne = new G4PVPlacement(0,
                                     G4ThreeVector(0., 0., 0.),
                                     "OuterCryostat",
                                     fLogicLicorne,
                                     fMotherVolume,
                                     false,
                                     0,
                                     myCheckOverlap);

  //---------------------------//
  //         Vacuum            //
  //--------------------------//
  G4Tubs *fSolidCryoOuterVac = new G4Tubs("CryoOuterVac_Solid",
                                          0.,
                                          cryoOuterR - cryoOuterThickness,
                                          cryoOuterH - cryoOuterThickness,
                                          0.,
                                          2 * M_PI);

  G4LogicalVolume *fLogicLicorneVac = new G4LogicalVolume(fSolidCryoOuterVac,
                                                          DSMaterial::Get()->GetVacuum(),
                                                          "OuterCryostat_LogicVac");

  fPhysicLicorneVac = new G4PVPlacement(0,
                                        G4ThreeVector(0., 0., 0.),
                                        "OuterCryostatVac",
                                        fLogicLicorneVac,
                                        fPhysicLicorne,
                                        false,
                                        0,
                                        myCheckOverlap);

  //-------------------------//
  //      Inner Cryo         //
  //-------------------------//

  G4Tubs *fSolidInnerCryo = new G4Tubs("InnerCryo_Solid",
                                       0.,
                                       cryoInnerR,
                                       cryoInnerH,
                                       0.,
                                       2 * M_PI);

  G4LogicalVolume *fLogicInnerCryo = new G4LogicalVolume(fSolidInnerCryo,
                                                         DSMaterial::Get()->GetStainlessSteel(),
                                                         "SolidInnerCryo_Logic");

  fPhysicInnerCryo = new G4PVPlacement(0,
                                       myZeros,
                                       "InnerCryo",
                                       fLogicInnerCryo,
                                       fPhysicLicorneVac,
                                       false,
                                       0,
                                       myCheckOverlap);

  //-------------------------//
  //  Non scintillating LAr  //
  //-------------------------//

  G4Tubs *fSolidInactiveLar = new G4Tubs("InactiveLar_Solid",
                                         0,
                                         cryoInnerR - cryoInnerThickness,
                                         (cryoInnerH - cryoInnerThickness),
                                         0,
                                         2 * M_PI);

  G4LogicalVolume *fLogicInactiveLar = new G4LogicalVolume(fSolidInactiveLar,
                                                           DSMaterial::Get()->GetNSLiquidArgon(),
                                                           "SolidInactiveLar_Logic");

  // G4ThreeVector myInactiveLarShift (0,0,-(cryoInnerH-cryoInnerThickness)/4);

  fPhysicInactiveLar = new G4PVPlacement(0,
                                         myZeros,
                                         "InactiveLar",
                                         fLogicInactiveLar,
                                         fPhysicInnerCryo,
                                         false,
                                         0,
                                         myCheckOverlap);

  //-----------------------//
  //    Teflon Sleeve (1 piece for Sleeve, Main Body and Inner Sleeve            //
  //-----------------------//

  G4Tubs *fSolidExternalSleeve = new G4Tubs("ExternalSleeve_Solid",
                                            0,
                                            SleeveOuterR,
                                            SleeveOuterH,
                                            0,
                                            2 * M_PI);

  G4LogicalVolume *fLogicExternalSleeve = new G4LogicalVolume(fSolidExternalSleeve,
                                                              DSMaterial::Get()->GetTeflon(),
                                                              "SolidExternalSleeve_Logic");

  fPhysicExternalSleeve = new G4PVPlacement(0,
                                            myZeros,
                                            "ExternalSleeve",
                                            fLogicExternalSleeve,
                                            fPhysicInactiveLar,
                                            false,
                                            0,
                                            myCheckOverlap);

  //------------------------//
  //      TPB Layers        //
  //------------------------//

  // Layer surrounding GasPocket
  G4Tubs *fSolidTPBGAr = new G4Tubs("TPBLayerGAr_Solid",
                                    0,
                                    TPCLArOuterR,
                                    TPCGArOuterH,
                                    0,
                                    2 * M_PI);

  G4LogicalVolume *fLogicTPBGAr = new G4LogicalVolume(fSolidTPBGAr,
                                                      DSMaterial::Get()->GetTPB(),
                                                      "TPBGAr_Logic");

  fPhysicTPBGAr = new G4PVPlacement(0,
                                    GasPocketShift,
                                    "TPBGar",
                                    fLogicTPBGAr,
                                    fPhysicExternalSleeve,
                                    false,
                                    0,
                                    myCheckOverlap);

  // Layer surrounding ActiveLAr
  G4Tubs *fSolidTPBLAr = new G4Tubs("TPBLayerLAr_Solid",
                                    0,
                                    TPCLArOuterR,
                                    TPCLArOuterH,
                                    0,
                                    2 * M_PI);

  G4LogicalVolume *fLogicTPBLAr = new G4LogicalVolume(fSolidTPBLAr,
                                                      DSMaterial::Get()->GetTPB(),
                                                      "TPBLAr_Logic");

  fPhysicTPBLAr = new G4PVPlacement(0,
                                    LArTPCShift,
                                    "TPBLAr",
                                    fLogicTPBLAr,
                                    fPhysicExternalSleeve,
                                    false,
                                    0,
                                    myCheckOverlap);

  //-----------------------//
  //    Copper Rings        //
  //-----------------------//

  G4double RingOuterR = GetCm(3.44 / 2.) * cm;
  G4double RingInnerR = GetCm(3.42 / 2.) * cm;
  G4double RingH = GetCm(0.1 / 2.) * cm;
  // G4double SleeveH = GetCm(3.00 / 2.) * cm;
  G4double RingSpace = GetCm(0.14) * cm;

  G4Tubs *fSolidCopperRing = new G4Tubs("CopperRing_Solid",
                                        RingInnerR,
                                        RingOuterR,
                                        RingH,
                                        0,
                                        2 * M_PI);

  G4LogicalVolume *fLogicCopperRing = new G4LogicalVolume(fSolidCopperRing,
                                                          DSMaterial::Get()->GetMetalCopper(),
                                                          "CopperRing_Logic");

  // define 20 rings
  G4ThreeVector RingPlacement(0., 0., 0.);
  char ring_name[256];

  for (int j = 1; j < 20; j++) {
    snprintf(ring_name, 30, "CopperRing_%d", j);
    RingPlacement.setZ(TPCLArOuterH - RingSpace * j + LArTPCShift.getZ());
    // RingPlacement = RingPlacement - LArTPCShift;
    fPhysicCopperRing[j] = new G4PVPlacement(0,
                                       RingPlacement,
                                       ring_name,
                                       fLogicCopperRing,
                                       fPhysicExternalSleeve,
                                       false,
                                       0,
                                       myCheckOverlap);
  }

  //-----------------------//
  //    LAr TPC (Active_LAr        //
  //-----------------------//

  G4Tubs *fSolidActiveLArTPC = new G4Tubs("ActiveLArTPC_Solid",
                                          0,
                                          TPCLArOuterR - 0.1 * mm,
                                          TPCLArOuterH - 0.05 * mm,
                                          0,
                                          2 * M_PI);

  G4LogicalVolume *fLogicActiveLArTPC = new G4LogicalVolume(fSolidActiveLArTPC,
                                                            DSMaterial::Get()->GetLiquidArgon(),
                                                            "LAr_Logic");

  fActiveLArTPC = new G4PVPlacement(0,
                                    G4ThreeVector(0, 0, 0.05 * mm),
                                    "ActiveLAr",
                                    fLogicActiveLArTPC,
                                    fPhysicTPBLAr,
                                    false,
                                    0,
                                    myCheckOverlap);


  //-----------------------//
  //    Gas Pocket       //
  //-----------------------//
  G4Tubs *fSolidGasPocket = new G4Tubs("GasPocket_Solid",
                                       0,
                                       TPCLArOuterR - 0.1 * mm,
                                       TPCGArOuterH - 0.05 * mm,
                                       0,
                                       2 * M_PI);

  G4LogicalVolume *fLogicGasPocket = new G4LogicalVolume(fSolidGasPocket,
                                                         DSMaterial::Get()->GetGaseousArgon(),
                                                         "GasPocket_Logic");

  fPhysicGasPocket = new G4PVPlacement(0,
                                       G4ThreeVector(0, 0, -0.05 * mm),
                                       "GasPocket",
                                       fLogicGasPocket,
                                       fPhysicTPBGAr,
                                       false,
                                       0,
                                       myCheckOverlap);

  //-----------------------//
  //    Grid                //
  //-----------------------//

  G4Tubs *fSolidGrid = new G4Tubs("Grid_Solid",
                                  0,
                                  TPCLArOuterR - 0.1 * mm,
                                  GridThickness,
                                  0,
                                  2 * M_PI);

  G4LogicalVolume *fLogicGrid = new G4LogicalVolume(fSolidGrid,
                                                    DSMaterial::Get()->GetStainlessSteel(),
                                                    "Grid_Logic");

  fPhysicGrid = new G4PVPlacement(0,
                                  GridShift,
                                  "Grid",
                                  fLogicGrid,
                                  fPhysicExternalSleeve,
                                  false,
                                  0,
                                  myCheckOverlap);

  //-----------------------//
  //    Cathode Window       //
  //-----------------------//

  G4Tubs *fSolidCathodeWindow = new G4Tubs("CathodeWindow_Solid",
                                           0,
                                           FusedSilicaOuterR,
                                           FusedSilicaThickness,
                                           0,
                                           2 * M_PI);

  G4LogicalVolume *fLogicCathodeWindow = new G4LogicalVolume(fSolidCathodeWindow,
                                                             DSMaterial::Get()->GetFusedSilica(),
                                                             "CathodeWindow_Logic");

  fPhysicCathodeWindow = new G4PVPlacement(0,
                                           CathodeShift,
                                           "CathodeWindow",
                                           fLogicCathodeWindow,
                                           fPhysicExternalSleeve,
                                           false,
                                           0,
                                           myCheckOverlap);

  //-----------------------//
  //    Anode Window       //
  //-----------------------//
  G4Tubs *fSolidAnodeWindow = new G4Tubs("AnodeWindow_Solid",
                                         0,
                                         FusedSilicaOuterR,
                                         FusedSilicaThickness,
                                         0,
                                         2 * M_PI);

  G4LogicalVolume *fLogicAnodeWindow = new G4LogicalVolume(fSolidAnodeWindow,
                                                           DSMaterial::Get()->GetFusedSilica(),
                                                           "AnodeWindow_Logic");

  fPhysicAnodeWindow = new G4PVPlacement(0,
                                         AnodeShift,
                                         "AnodeWindow",
                                         fLogicAnodeWindow,
                                         fPhysicExternalSleeve,
                                         false,
                                         0,
                                         myCheckOverlap);

  //-----------------------//
  //    Upper PTFE ring     //
  //-----------------------//
  //

  //-----------------------//
  //    Top PMT arrays       //
  //-----------------------//

  //  PMT Window
  G4Box *fSolidTopPMTBody1 = new G4Box("TopPMTBody1_Solid",
                                       PMTTopBodyL,
                                       PMTTopBodyL,
                                       PMTTopBody1H);

  G4Box *fSolidTopPMTWindow = new G4Box("TopPMTBody1Window_Solid",
                                        PMTTopBodyL,
                                        PMTTopBodyL,
                                        PMTTopPhotoCathodeH);
  ;
  G4Box *fSolidTopPMTBody1_Vac = new G4Box("TopPMTBody1Vac_Solid",
                                           PMTTopBodyL_Vac,
                                           PMTTopBodyL_Vac,
                                           PMTTopBody1H_Vac);

  G4LogicalVolume *fLogicTopPMTBody1[7];
  G4LogicalVolume *fLogicTopPMTBody1_Window[7];
  G4LogicalVolume *fLogicTopPMTBody1_Vac[7];
  char logicTopPMT[256];
  char logicTopPMTWindow[256];
  char logicTopPMTVac[256];

  char physicTopPMT[256];
  char physicTopPMTWindow[256];
  char physicTopPMTVac[256];

  G4ThreeVector pmtwin_shift(0, 0, -PMTTopBody1H + PMTTopPhotoCathodeH);
  G4ThreeVector pmtvac_shift(0, 0, PMTTopPhotoCathodeH);  // 2*PMTTopPhotoCathodeH);

  for (int i = 0; i < 7; i++) {
    // shift 1.06 inches, 7 PMTs
    // pmt 0
    if (i == 0) {
      PMTTopBody1Shift.setX(0);
      PMTTopBody1Shift.setY(0);
    }
    // pmt 1
    else if (i == 1) {
      PMTTopBody1Shift.setX(GetCm(1.06) * cm);
      PMTTopBody1Shift.setY(0);
    }
    // pmt 2
    else if (i == 2) {
      PMTTopBody1Shift.setX(-GetCm(1.06) * cm);
      PMTTopBody1Shift.setY(0);
    }
    // pmt 3
    else if (i == 3) {
      PMTTopBody1Shift.setX(-GetCm(0.53) * cm);
      PMTTopBody1Shift.setY(GetCm(1.06) * cm);
    }
    // pmt 4
    else if (i == 4) {
      PMTTopBody1Shift.setX(GetCm(0.53) * cm);
      PMTTopBody1Shift.setY(GetCm(1.06) * cm);
    }
    // pmt 5
    else if (i == 5) {
      PMTTopBody1Shift.setX(GetCm(0.53) * cm);
      PMTTopBody1Shift.setY(-GetCm(1.06) * cm);
    }
    // pmt 6
    else if (i == 6) {
      PMTTopBody1Shift.setX(-GetCm(0.53) * cm);
      PMTTopBody1Shift.setY(-GetCm(1.06) * cm);
    }

    snprintf(logicTopPMT, 30, "TopPMTBody1_Logic_%d", i);
    snprintf(logicTopPMTWindow, 30, "Logic_TPMT_%d", i);
    snprintf(logicTopPMTVac, 30, "TopPMTBody1Vac_Logic_%d", i);

    fLogicTopPMTBody1[i] = new G4LogicalVolume(fSolidTopPMTBody1,
                                               DSMaterial::Get()->GetStainlessSteel(),
                                               logicTopPMT);

    fLogicTopPMTBody1_Window[i] = new G4LogicalVolume(fSolidTopPMTWindow,
                                                      DSMaterial::Get()->GetFakePhotocathode(),
                                                      logicTopPMT);

    fLogicTopPMTBody1_Vac[i] = new G4LogicalVolume(fSolidTopPMTBody1_Vac,
                                                   DSMaterial::Get()->GetVacuum(),
                                                   logicTopPMT);

    snprintf(physicTopPMT, 30, "TopPMTBody1_Physic_%d", i);
    // snprintf(physicTopPMTWindow,30, "TopPMTBody1Window_Physic_%d",i);
    snprintf(physicTopPMTWindow, 30, "TPMT_%d", i);
    snprintf(physicTopPMTVac, 30, "TopPMTBody1Vac_Physic_%d", i);

    fPhysicTopPMTBody1[i] = new G4PVPlacement(0,
                                              PMTTopBody1Shift,
                                              physicTopPMT,
                                              fLogicTopPMTBody1[i],
                                              fPhysicExternalSleeve,
                                              false,
                                              0,
                                              myCheckOverlap);

    fPhysicTopPMTBody1_Window[i] = new G4PVPlacement(0,
                                                     pmtwin_shift,
                                                     physicTopPMTWindow,
                                                     fLogicTopPMTBody1_Window[i],
                                                     fPhysicTopPMTBody1[i],
                                                     false,
                                                     0,
                                                     myCheckOverlap);

    fPhysicTopPMTBody1_Vac[i] = new G4PVPlacement(0,
                                                  pmtvac_shift,
                                                  physicTopPMTVac,
                                                  fLogicTopPMTBody1_Vac[i],
                                                  fPhysicTopPMTBody1[i],
                                                  false,
                                                  0,
                                                  myCheckOverlap);
  }

  //-----------------------//
  //    Bottom PMT         //
  //-----------------------//

  //  PMT Window
  G4Tubs *fSolidBottomPMTWindow = new G4Tubs("BottomPMTWindow_Solid",
                                             0,
                                             PMTBottomR,
                                             PMTBottomWindowH,
                                             0,
                                             2 * M_PI);

  G4LogicalVolume *fLogicBottomPMTWindow = new G4LogicalVolume(fSolidBottomPMTWindow,
                                                               DSMaterial::Get()->GetFakePhotocathode(),
                                                               "BottomPMTWindow_Logic");

  fPhysicBottomPMTWindow = new G4PVPlacement(0,
                                             PMTBottomWindowShift,
                                             "TPMT_7",
                                             fLogicBottomPMTWindow,
                                             fPhysicExternalSleeve,
                                             false,
                                             0,
                                             myCheckOverlap);

  /// defined in the region with Teflon
  G4Tubs *fSolidBottomPMTBody1 = new G4Tubs("BottomPMTBody1_Solid",
                                            0,
                                            PMTBottomBody1R,
                                            PMTBottomBody1H,
                                            0,
                                            2 * M_PI);

  G4LogicalVolume *fLogicBottomPMTBody1 = new G4LogicalVolume(fSolidBottomPMTBody1,
                                                              DSMaterial::Get()->GetStainlessSteel(),
                                                              "BottomPMTBody1_Logic");

  fPhysicBottomPMTBody1 = new G4PVPlacement(0,
                                            PMTBottomBody1Shift,
                                            "BottomPMTBody1",
                                            fLogicBottomPMTBody1,
                                            fPhysicExternalSleeve,
                                            false,
                                            0,
                                            myCheckOverlap);

  // vacuum
  G4Tubs *fSolidBottomPMTBody1Vac = new G4Tubs("BottomPMTBody1Vac_Solid",
                                               0,
                                               PMTBottomBody1R_vac,
                                               PMTBottomBody1H,
                                               0,
                                               2 * M_PI);

  G4LogicalVolume *fLogicBottomPMTBody1Vac = new G4LogicalVolume(fSolidBottomPMTBody1Vac,
                                                                 DSMaterial::Get()->GetVacuum(),
                                                                 "BottomPMTBody1Vac_Logic");

  fPhysicBottomPMTBody1Vac = new G4PVPlacement(0,
                                               myZeros,
                                               "BottomPMTBody1Vac",
                                               fLogicBottomPMTBody1Vac,
                                               fPhysicBottomPMTBody1,
                                               false,
                                               0,
                                               myCheckOverlap);

  ////defined in the region with NS LAr
  G4Tubs *fSolidBottomPMTBody2 = new G4Tubs("BottomPMTBody2_Solid",
                                            0,
                                            PMTBottomBody2R,
                                            PMTBottomBody2H,
                                            0,
                                            2 * M_PI);

  G4LogicalVolume *fLogicBottomPMTBody2 = new G4LogicalVolume(fSolidBottomPMTBody2,
                                                              DSMaterial::Get()->GetStainlessSteel(),
                                                              "BottomPMTBody2_Logic");

  fPhysicBottomPMTBody2 = new G4PVPlacement(0,
                                            PMTBottomBody2Shift,
                                            "BottomPMTBody2",
                                            fLogicBottomPMTBody2,
                                            fPhysicInactiveLar,
                                            false,
                                            0,
                                            myCheckOverlap);
  // vacuum
  G4Tubs *fSolidBottomPMTBody2Vac = new G4Tubs("BottomPMTBody2Vac_Solid",
                                               0,
                                               PMTBottomBody2R_Vac,
                                               PMTBottomBody2H,
                                               0,
                                               2 * M_PI);

  G4LogicalVolume *fLogicBottomPMTBody2Vac = new G4LogicalVolume(fSolidBottomPMTBody2Vac,
                                                                 DSMaterial::Get()->GetVacuum(),
                                                                 "BottomPMTBody2Vac_Logic");

  fPhysicBottomPMTBody2Vac = new G4PVPlacement(0,
                                               myZeros,
                                               "BottomPMTBody2Vac",
                                               fLogicBottomPMTBody2Vac,
                                               fPhysicBottomPMTBody2,
                                               false,
                                               0,
                                               myCheckOverlap);

  G4Region *fLArRegion = new G4Region("LAr_Logic");
  fLogicActiveLArTPC->SetRegion(fLArRegion);
  fLArRegion->AddRootLogicalVolume(fLogicActiveLArTPC);

  DefineSurfaces();
  DSStorage::Get()->SetPMTMaterialIndex(fLogicBottomPMTWindow->GetMaterial()->GetIndex());

  //--------------------------------------------------
  // Configure visualization
  //--------------------------------------------------
  G4Color colTransparent = G4Color(0, 0, 0, 0);
  G4Color colPMT = G4Color(0.75, 0.75, 0.75, 1);
  G4Color colCopper = G4Color(0.722, 0.451, 0.2, 0.3);
  G4Color colSteel = G4Color(0.533, 0.545, 0.553, 0.3);
  // G4Color colLAr = G4Color(0.952, 0.952, 0.956, 0.5);  // G4Color(0.831, 0.945, 0.976, 0.5)
  G4Color colLAr = colTransparent;
  G4Color colLArInactive = colLAr;
  G4Color colGAr = colLAr;
  G4Color colTeflon = G4Color(0.25, 0.25, 0.25, 0.3);
  G4Color colWindow = G4Color(0.58, 0.576, 0.69, 0.3);

  fMotherVolume->GetLogicalVolume()->SetVisAttributes(G4VisAttributes::GetInvisible());
  fLogicLicorne->SetVisAttributes(colSteel);
  fLogicLicorneVac->SetVisAttributes(colTransparent);
  fLogicInnerCryo->SetVisAttributes(colSteel);
  fLogicInactiveLar->SetVisAttributes(colLArInactive);
  fLogicExternalSleeve->SetVisAttributes(colTeflon);
  fLogicCopperRing->SetVisAttributes(colCopper);
  fLogicActiveLArTPC->SetVisAttributes(colLAr);
  fLogicGasPocket->SetVisAttributes(colGAr);
  fLogicGrid->SetVisAttributes(colSteel);
  fLogicCathodeWindow->SetVisAttributes(colWindow);
  fLogicAnodeWindow->SetVisAttributes(colWindow);
  for (int i = 0; i < 7; i++) {
    fLogicTopPMTBody1[i]->SetVisAttributes(colPMT);
    fLogicTopPMTBody1_Window[i]->SetVisAttributes(colPMT);
    fLogicTopPMTBody1_Vac[i]->SetVisAttributes(colPMT);
  }
  fLogicBottomPMTWindow->SetVisAttributes(colPMT);
  fLogicBottomPMTBody1->SetVisAttributes(colPMT);
  fLogicBottomPMTBody1Vac->SetVisAttributes(colPMT);
  fLogicBottomPMTBody2->SetVisAttributes(colPMT);
  fLogicBottomPMTBody2Vac->SetVisAttributes(colPMT);
}

DSDetectorARISER::~DSDetectorARISER() {
  ;  // delete fMessenger;
}

void DSDetectorARISER::DefineSurfaces() {
  ////////////////////////////////////////
  // Copied from DS DetectorDS50 - 24th april 2015
  ////////////////////////////////////////

  ////////////////////////////////////////
  // LAr -> Grid
  //  Note: in this model, all optical action takes place
  //  on entering the grid.
  ////////////////////////////////////////

  G4OpticalSurface *fOpGridLArSurface = new G4OpticalSurface("OpGridLArSurface");
  new G4LogicalBorderSurface("GridLArSurface", fActiveLArTPC, fPhysicGrid, fOpGridLArSurface);
  fOpGridLArSurface->SetType(dielectric_dielectric);
  fOpGridLArSurface->SetModel(glisur);
  // fOpGridLArSurface->SetModel( unified );
  fOpGridLArSurface->SetFinish(polished);
  // fOpGridLArSurface -> SetPolish(0.1);
  // fOpGridLArSurface->SetFinish( ground );
  // fOpGridLArSurface->SetSigmaAlpha(0.1);
  G4double GridLArENE[4] = {0.1 * eV, 8.0 * eV, 8.3 * eV, 20.0 * eV};
  G4double GRUV = DSParameters::Get()->GetLArGridUVRef();
  G4double GRVIS = DSParameters::Get()->GetLArGridVisRef();
  G4double GridLArREF[4] = {GRVIS, GRVIS, GRUV, GRUV};
  G4MaterialPropertiesTable *fGridLArSurfProp = new G4MaterialPropertiesTable();
  if (DSParameters::Get()->GetWithNewGridModel())
    // the grid model is described in DSStorage.cc, the surface is treated in G4OpBoundaryProcess.cc
    fGridLArSurfProp->AddConstProperty("DOGRID", 1, true);
  // Now use the following in old and new models.  By G4 convention, "reflectivity" is actually 1-absorption.
  fGridLArSurfProp->AddProperty("REFLECTIVITY", GridLArENE, GridLArREF, 4);
  //  else {
  //    fGridLArSurfProp->AddProperty("REFLECTIVITY", GridLArENE, GridLArREF, 4);
  //    }
  fOpGridLArSurface->SetMaterialPropertiesTable(fGridLArSurfProp);

  ////////////////////////////////////////
  // Grid -> LAr (keeping backward labeling convention)
  //  Note: in this model, all optical action takes place
  //  on entering the grid.  Exit action is to just continue
  //  in straight line.
  ////////////////////////////////////////

  if (DSParameters::Get()->GetWithNewGridModel()) {
    G4OpticalSurface *fOpLArGridSurface = new G4OpticalSurface("OpLArGridSurface");
    new G4LogicalBorderSurface("LArGridSurface", fPhysicGrid, fActiveLArTPC, fOpLArGridSurface);

    fOpLArGridSurface->SetType(dielectric_dielectric);
    //  fOpLArGridSurface->SetModel( glisur );
    //  fOpLArGridSurface->SetFinish( polished );
    G4MaterialPropertiesTable *fLArGridSurfProp = new G4MaterialPropertiesTable();
    // the grid model is described in DSStorage.cc, the surface is treated in G4OpBoundaryProcess.cc
    fLArGridSurfProp->AddConstProperty("DOGRIDEXIT", 1, true);
    fOpLArGridSurface->SetMaterialPropertiesTable(fLArGridSurfProp);
  }

  ////////////////////////////////////////
  // TPB <--> GAr and TPB <--> LAr
  // This surface will carry all the diffuse properties of the TPB
  // for both GAr and LAr.
  // Make this bi-directional
  ////////////////////////////////////////
  G4OpticalSurface *fOpTPBGArSurface = new G4OpticalSurface("OpTPBGArSurface");
  // Note: the following two work even when the "gas pocket" is LAr in no-pocket runs.
  new G4LogicalBorderSurface("GArTPBSurface", fPhysicGasPocket, fPhysicTPBGAr, fOpTPBGArSurface);
  new G4LogicalBorderSurface("TPBGArSurface", fPhysicTPBGAr, fPhysicGasPocket, fOpTPBGArSurface);
  new G4LogicalBorderSurface("LArTPBSurface", fActiveLArTPC, fPhysicTPBLAr, fOpTPBGArSurface);
  new G4LogicalBorderSurface("TPBLArSurface", fPhysicTPBLAr, fActiveLArTPC, fOpTPBGArSurface);
  // new G4LogicalBorderSurface("TPBLArLayerSurface", fPhysicTPB, fPhysicLArLayer, fOpTPBGArSurface );
  fOpTPBGArSurface->SetType(dielectric_dielectric);
  fOpTPBGArSurface->SetModel(unified);
  fOpTPBGArSurface->SetFinish(ground);
  fOpTPBGArSurface->SetSigmaAlpha(0.3);

  G4double VISTRAN = DSParameters::Get()->GetArTPBVisTran();
  const G4int NUM = 4;
  G4double pp[NUM] = {0.1 * eV, 8.0 * eV, 8.3 * eV, 20.0 * eV};  // vis, vis, UV, UV
  G4double specularlobe[NUM] = {0., 0., 0., 0.};                 //--
  G4double specularspike[NUM] = {0., 0., 0., 0.};                //----  gives all reflection to Lambertian lobe
  G4double backscatter[NUM] = {0., 0., 0., 0.};                  //--
  G4double reflectivity[NUM] = {1.0, 1.0, 1.0, 1.0};             //  To set 1-absorption
  G4double transmitivity[NUM] = {VISTRAN, VISTRAN, 1.0, 1.0};    //  To set reflection vs. transmission, overridding Fresnel
                                                                 //  For now, no angle dependence.

  G4MaterialPropertiesTable *fTPBGArSurfProp = new G4MaterialPropertiesTable();

  fTPBGArSurfProp->AddProperty("SPECULARLOBECONSTANT", pp, specularlobe, NUM);
  fTPBGArSurfProp->AddProperty("SPECULARSPIKECONSTANT", pp, specularspike, NUM);
  fTPBGArSurfProp->AddProperty("BACKSCATTERCONSTANT", pp, backscatter, NUM);
  fTPBGArSurfProp->AddProperty("REFLECTIVITY", pp, reflectivity, NUM);
  fTPBGArSurfProp->AddProperty("TRANSMITTANCE", pp, transmitivity, NUM);

  fTPBGArSurfProp->AddConstProperty("DOArTPB", 1, true);

  fOpTPBGArSurface->SetMaterialPropertiesTable(fTPBGArSurfProp);

  ////////////////////////////////////////
  // ITO /////
  ////////////////////////////////////////

  G4MaterialPropertiesTable *fITOSurfProp = new G4MaterialPropertiesTable();
  fITOSurfProp->AddConstProperty("DOITO", 1, true);

  ////////////////////////////////////////
  // AnodeWindow <--> TPB and CathodeWindow <--> TPB
  // In the current model, the diffuse nature of the TPB is handled
  // entirely at the TPB-GAr/LAr surface, not here.
  // Make this bi-directional.
  ////////////////////////////////////////
  G4OpticalSurface *fOpWindowTPBSurface = new G4OpticalSurface("OpWindowTPBSurface");
  new G4LogicalBorderSurface("BottomWindowTPBSurface", fPhysicCathodeWindow, fPhysicTPBLAr, fOpWindowTPBSurface);
  new G4LogicalBorderSurface("TPBBottomWindowSurface", fPhysicTPBLAr, fPhysicCathodeWindow, fOpWindowTPBSurface);
  new G4LogicalBorderSurface("TopWindowTPBSurface", fPhysicAnodeWindow, fPhysicTPBGAr, fOpWindowTPBSurface);
  new G4LogicalBorderSurface("TPBTopWindowSurface", fPhysicTPBGAr, fPhysicAnodeWindow, fOpWindowTPBSurface);
  fOpWindowTPBSurface->SetType(dielectric_dielectric);
  fOpWindowTPBSurface->SetModel(unified);
  //  fOpWindowTPBSurface->SetFinish( ground );
  fOpWindowTPBSurface->SetFinish(polished);
  // G4MaterialPropertiesTable *fWindowTPBSurfProp = new G4MaterialPropertiesTable();
  fOpWindowTPBSurface->SetMaterialPropertiesTable(fITOSurfProp);

  ///////////////
  // TPB --> Teflon (Reflector)
  //  Should be no Teflon --> TPB as surface is defined as dielectric_metal.
  ////////////////////////////////////////
  G4double TeflonTPBENE[4] = {0.1 * eV, 8.0 * eV, 8.3 * eV, 20.0 * eV};
  G4double TREFUV = DSParameters::Get()->GetTeflonTPBUVRef();
  G4double TREFVIS = DSParameters::Get()->GetTeflonTPBVisRef();
  G4double TeflonTPBREF[4] = {TREFVIS, TREFVIS, TREFUV, TREFUV};

  G4OpticalSurface *fOpTPBTeflonSurface = new G4OpticalSurface("OpTBPTeflonSurface");
  new G4LogicalBorderSurface("TPBLArTeflonSurface", fPhysicTPBLAr, fPhysicExternalSleeve, fOpTPBTeflonSurface);
  new G4LogicalBorderSurface("TPBGArTeflonSurface", fPhysicTPBGAr, fPhysicExternalSleeve, fOpTPBTeflonSurface);
  fOpTPBTeflonSurface->SetType(dielectric_metal);
  fOpTPBTeflonSurface->SetModel(unified);

  // PDM: though I can't see how, the following settings are giving the desired Lambertian reflection
  fOpTPBTeflonSurface->SetFinish(groundfrontpainted);
  // fOpTPBTeflonSurface->SetFinish(ground);
  fOpTPBTeflonSurface->SetSigmaAlpha(0.1);

  G4MaterialPropertiesTable *fTPBTeflonSurfProp = new G4MaterialPropertiesTable();
  fTPBTeflonSurfProp->AddProperty("REFLECTIVITY", TeflonTPBENE, TeflonTPBREF, 4);
  fOpTPBTeflonSurface->SetMaterialPropertiesTable(fTPBTeflonSurfProp);

  ////////////////////////////////////////
  // Teflon - LAr
  //  PDM: These refer to EXTERNAL surfaces of the TPC.  I don't think they will be executed
  //  and they haven't been debugged.  (Some are defined only for teflon --> LAr, which shouldn't happen.)
  ////////////////////////////////////////
  G4OpticalSurface *fOpLArTeflonSurface = new G4OpticalSurface("OpLArTeflonSurface");
  new G4LogicalBorderSurface("LArTeflonSurface", fPhysicExternalSleeve, fPhysicInactiveLar, fOpLArTeflonSurface);
  // new G4LogicalBorderSurface("LArPMTTopAssemblySurface", fPhysicPMTAssemblyTop, fPhysicInactiveLAr, fOpLArTeflonSurface);
  // new G4LogicalBorderSurface("LArPMTBottomAssemblySurface", fPhysicPMTAssemblyBottom, fPhysicInactiveLAr, fOpLArTeflonSurface);
  new G4LogicalBorderSurface("LArTeflonSupportSurface", fPhysicInactiveLar, fPhysicExternalSleeve, fOpLArTeflonSurface);
  fOpLArTeflonSurface->SetType(dielectric_metal);
  fOpLArTeflonSurface->SetModel(glisur);
  fOpLArTeflonSurface->SetFinish(ground);
  G4MaterialPropertiesTable *fLArTeflonSurfProp = new G4MaterialPropertiesTable();
  G4double TeflonLArENE[4] = {0.1 * eV, 8.0 * eV, 8.3 * eV, 20.0 * eV};
  G4double TAREFUV = DSParameters::Get()->GetTeflonLArUVRef();
  G4double TAREFVIS = DSParameters::Get()->GetTeflonLArVisRef();
  G4double TeflonLArREF[4] = {TAREFVIS, TAREFVIS, TAREFUV, TAREFUV};
  fLArTeflonSurfProp->AddProperty("REFLECTIVITY", TeflonLArENE, TeflonLArREF, 4);
  fOpLArTeflonSurface->SetMaterialPropertiesTable(fLArTeflonSurfProp);

  return;
}

double DSDetectorARISER::GetCm(double inches) { return inches * 2.54; }
