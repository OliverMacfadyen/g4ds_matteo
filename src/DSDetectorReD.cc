#include "G4Orb.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
#include "G4ThreeVector.hh"

#include "DSDetectorNeutronReD.hh"
#include "DSDetectorReD.hh"
#include "DSIO.hh"
#include "DSLogger.hh"
#include "DSMaterial.hh"
#include "DSParameters.hh"
#include "DSStorage.hh"
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4Cons.hh"
#include "G4Sphere.hh"
#include "G4VisAttributes.hh"

#include "G4NistManager.hh"

#include <iostream>
#include <stdio.h>

using namespace std;

DSDetectorReD::DSDetectorReD(G4VPhysicalVolume* myMotherVolume) {

  G4Colour myWhite(1.0, 1.0, 1.0);      // white
  G4Colour myGray(0.5, 0.5, 0.5);       // gray
  G4Colour myltGray(0.25, 0.25, 0.25);  // light gray
  G4Colour mydkGray(0.75, 0.75, 0.75);  // dark gray
  G4Colour myBlack(0.0, 0.0, 0.0);      // black
  G4Colour myRed(1.0, 0.0, 0.0);        // red
  G4Colour myGreen(0.0, 1.0, 0.0);      // green
  G4Colour myBlue(0.0, 0.0, 1.0);       // blue
  G4Colour myCyan(0.0, 1.0, 1.0);       // cyan
  G4Colour myMagenta(1.0, 0.0, 1.0);    // magenta
  G4Colour myYellow(1.0, 1.0, 0.0);     // yellow
  G4Colour myPink(0.94, 0.5, 0.5);      // pink

  fMotherVolume = myMotherVolume;
  DSLog(routine) << " Constructing ReD Setup Geometry" << endlog;

  // const double myTwoPi = 2 * M_PI * rad;
  //  Hey Michael Kuss, remember this for all times:
  //  THE OVERLAP CHECK OUTPUT GOES TO STDOUT, NOT TO THE LOG FILE!
  G4bool myCheckOverlap = DSStorage::Get()->GetCheckOverlap();
  G4ThreeVector myZeros(0., 0., 0.);

  // construct the array of LScis if a survey file is given
  if (DSStorage::Get()->GetReDNDSurveyFileName() != "") new DSDetectorNeutronReD(fMotherVolume);

  // ReD cryostat configuration
  //
  // I would have preferred to to do this as bitmap,
  // but unfortunately g4 doesn't accept easily hex or binary input to int.
  // And always thinking hex or binary and then writing decimal is error prone.
  //
  // So, I do it decimal from the beginning:
  // last two digits (0-99): code for the cryostat, right now:
  //                         0: old 65 cm Naples cryostat
  //                         1: new Criotec cryostat
  //                         99: ARIS cryostat
  // third digit:            0 always
  // digits 4 and 5:         vacuum flags:
  //                            0: all materials are what they are
  //                         1000: world volume material (air) -> vacuum
  //                         2000: stainless steel -> vacuum
  //                         4000: inactive LAr -> vacuum
  //                         8000: everything else (except for active LAr) ->
  //                         vacuum
  //                        16000: active LAr -> vacuum
  //
  // Example: setting ssteel and active LAr to vacuum for the ARIS-like cryostat
  // would be: 99 + 2000 + 16000 = 18099

  const G4int ReDConf = DSStorage::Get()->GetReDConfiguration();
  DSLog(routine) << "ReD cryostat configuration: " << ReDConf << endlog;
  G4int vacuumFlag = ReDConf / 1000;
  DSLog(routine) << "vacuum flag " << vacuumFlag << endlog;
  G4int cryostatFlag = ReDConf % 100;
  DSLog(routine) << "cryostat flag " << cryostatFlag << endlog;

  // ---------------------------------------------------//
  //                     the material                   //
  // ---------------------------------------------------//

  G4Material* Vacuum = DSMaterial::Get()->GetVacuum();  //  4

  G4Material* LiquidArgon = DSMaterial::Get()->GetLiquidArgon();      //  8
  G4Material* NSLiquidArgon = DSMaterial::Get()->GetNSLiquidArgon();  // 10
  G4Material* GaseousArgon = DSMaterial::Get()->GetGaseousArgon();    // 11

  G4Material* Steel = DSMaterial::Get()->GetSteel();                    // 14
  G4Material* StainlessSteel = DSMaterial::Get()->GetStainlessSteel();  // 15
  G4Material* MetalCopper = DSMaterial::Get()->GetMetalCopper();        // 24
  G4Material* Aluminum = DSMaterial::Get()->GetAluminum();              // 39

  G4Material* Acrylic = DSMaterial::Get()->GetAcrylic();            //  7
  G4Material* Teflon = DSMaterial::Get()->GetTeflon();              // 19
  G4Material* TPB = DSMaterial::Get()->GetTPB();                    // 28
  G4Material* ThreeMFoil = DSMaterial::Get()->GetThreeMFoil();      // 29
  G4Material* ITO = DSMaterial::Get()->GetITO();                    // 32
  G4Material* MetalSilicon = DSMaterial::Get()->GetMetalSilicon();  // 56

  // G4NistManager* manager = G4NistManager::Instance();  //added by Simone on
  // 25 Oct. 2017
  //   G4Material* Mylar = manager->FindOrBuildMaterial("G4_MYLAR");
  //   G4Material* Piombo = manager->FindOrBuildMaterial("G4_Pb");//MetalLead
  // G4Material* CH2 = manager->FindOrBuildMaterial("G4_PARAFFIN");

  // Visualization attributes
  //   TODO->SetVisAttributes(WhatVis);
  G4VisAttributes* MetalCopperVis = new G4VisAttributes(myPink);
  G4VisAttributes* StainlessSteelVis = new G4VisAttributes(myGray);
  G4VisAttributes* LiquidArgonVis = new G4VisAttributes(myRed);
  G4VisAttributes* NSLiquidArgonVis = new G4VisAttributes(myMagenta);
  G4VisAttributes* GaseousArgonVis = new G4VisAttributes(myYellow);

  // G4VisAttributes* allVis = new G4VisAttributes(mydkGray);
  //   allVis->SetForceSolid(true);
  // G4VisAttributes* targetVis = new G4VisAttributes(myBlue);
  //   targetVis->SetForceSolid(true);
  // G4VisAttributes* windowVis = new G4VisAttributes(myRed);
  //   windowVis->SetForceSolid(true);
  // G4VisAttributes* tubeVis = new G4VisAttributes(myYellow);
  //   tubeVis->SetForceSolid(true);

  if (vacuumFlag & 1) {
    DSLog(routine) << "Setting world volume material to vacuum" << endlog;
    G4Material* WorldMat = fMotherVolume->GetLogicalVolume()->GetMaterial();
    SetToVacuum(WorldMat);
    fMotherVolume->GetLogicalVolume()->SetMaterial(WorldMat);
  }

  if (vacuumFlag & 2) {
    SetToVacuum(StainlessSteel);
    SetToVacuum(Steel);
    SetToVacuum(Aluminum);
  }
  if (vacuumFlag & 4) { SetToVacuum(NSLiquidArgon); }
  if (vacuumFlag & 8) {
    SetToVacuum(MetalSilicon);
    SetToVacuum(MetalCopper);
    SetToVacuum(Teflon);
    SetToVacuum(TPB);
    SetToVacuum(ITO);
    SetToVacuum(GaseousArgon);
    SetToVacuum(Acrylic);
    SetToVacuum(ThreeMFoil);
  }
  if (vacuumFlag & 16) { SetToVacuum(LiquidArgon); }

  //  new CataniaBeamline(fMotherVolume,
  //  G4ThreeVector(-118.496*cm, 31.3385*cm, 48.3399*cm));

  // let's make this file even messier as it is ...
  // one day (soon) all these different cryostats (0, 1, 99) should have their
  // own class in their own files. right now, quick patching to get things done

  // Cryostat

  G4ThreeVector TPCshift(myZeros);
  G4PVPlacement* cryoOuterPhys = 0;
  if (cryostatFlag == 1) {

    // MK: vertical shift estimated
    TPCshift = G4ThreeVector(0, 0, -50 * mm);  // within inactiveLAr
    // MK: vertical shift estimated
    G4ThreeVector cryostatShift = TPCshift + G4ThreeVector(0, 0, 227.5 * mm);

    //------------------------//
    //      Outer Cryostat    //
    //------------------------//

    G4double cryoOuterH = 560.0 * mm;
    G4double cryoOuterR = 168.3 / 2.0 * mm;
    G4double cryoOuterW = 2.0 * mm;  // in the part list (4) it says x2, as well as for the end cap
                                     // (3).  But the measure below the end cap says x3!
    G4double cryoOuterWindowH = 280.0 * mm;
    G4double cryoOuterWindowW = 1.0 * mm;
    G4double cryoOuterTopRingH = 8.0 * mm;  // estimated, couldn't find measure
    G4double cryoOuterr = cryoOuterR - cryoOuterW;
    G4double cryoOuterWindowPosZ = (cryoOuterWindowH - cryoOuterH) / 2.0 + 10.0 * mm;
    G4double cryoOuterEndCapPosZ = -cryoOuterH / 2.0;  // a half sphere ("southern" hemisphere) has it's ref
                                                       // point anyway at 0,0,0
    G4double cryoOuterTopRingPosZ = cryoOuterH / 2.0 + cryoOuterTopRingH / 2.0;
    G4ThreeVector cryoOuterWindowPos(0, 0, cryoOuterWindowPosZ);
    G4ThreeVector cryoOuterEndCapPos(0, 0, cryoOuterEndCapPosZ);
    G4ThreeVector cryoOuterTopRingPos(0, 0, cryoOuterTopRingPosZ);

    // forward assignment
    G4double cryoInnerR = 131.0 / 2.0 * mm;
    /////////////////////

    G4Tubs* cryoOuterFullSolid = new G4Tubs("cryoOuterFullSolid", cryoOuterr, cryoOuterR, cryoOuterH / 2.0, 0.0, 2.0 * M_PI);
    // do the turning of the 1 mm entrance window (NOTA: without the extra 1 mm
    // for the outer radius it won't (may not) work
    G4Tubs* cryoOuterSolidTurnedOff = new G4Tubs("cryoOuterSolidTurnedOff", cryoOuterr + cryoOuterWindowW, cryoOuterR + 1 * mm, cryoOuterWindowH / 2.0, 0.0, 2.0 * M_PI);
    G4SubtractionSolid* cryoOuterSolid2 = new G4SubtractionSolid("cryoOuterSolid2", cryoOuterFullSolid, cryoOuterSolidTurnedOff, 0, cryoOuterWindowPos);
    // approximate the end cap by a sphere right now
    G4Sphere* cryoOuterEndCap = new G4Sphere("cryoOuterEndCap", cryoOuterr, cryoOuterR, 0, 2.0 * M_PI, 0.5 * M_PI, M_PI);
    G4UnionSolid* cryoOuterSolid1 = new G4UnionSolid("cryoOuterSolid2+cryoOuterEndCap", cryoOuterSolid2, cryoOuterEndCap, 0, cryoOuterEndCapPos);
    G4Tubs* cryoOuterTopRing = new G4Tubs("cryoOuterTopRing", cryoInnerR, cryoOuterR, cryoOuterTopRingH / 2.0, 0, 2.0 * M_PI);
    G4UnionSolid* cryoOuterSolid = new G4UnionSolid("cryoOuterSolid1+cryoOuterTopRing", cryoOuterSolid1, cryoOuterTopRing, 0, cryoOuterTopRingPos);
    G4LogicalVolume* cryoOuterLogic = new G4LogicalVolume(cryoOuterSolid, StainlessSteel, "cryoOuterLogic");
    cryoOuterLogic->SetVisAttributes(StainlessSteelVis);

    //------------------------//
    //      Vacuum            //
    //------------------------//

    G4double cryoVacH = cryoOuterH - cryoOuterTopRingH / 2.0;
    G4ThreeVector cryoVacEndCapPos(0, 0, cryoVacH / 2.0);

    G4Tubs* cryoVacSolid1 = new G4Tubs("cryoVacSolid1", cryoInnerR, cryoOuterr, cryoVacH / 2.0, 0, 2.0 * M_PI);
    // approximate the end cap by a sphere right now
    G4Sphere* cryoVacEndCapSolid1 = new G4Sphere("cryoVacEndCapSolid1", 0, cryoOuterr, 0, 2.0 * M_PI, 0.5 * M_PI, M_PI);

    // forward assignment
    G4double cryoInnerBotH = 5.0 * mm;
    G4Tubs* cryoInnerBot = new G4Tubs("cryoInnerBot", 0, cryoInnerR, cryoInnerBotH / 2.0, 0, 2.0 * M_PI);
    /////////////////////

    G4SubtractionSolid* cryoVacEndCapSolid = new G4SubtractionSolid("cryoVacEndCapSolid", cryoVacEndCapSolid1, cryoInnerBot);

    // If I do the following, running g4ds says:
    // ERROR: G4VSceneHandler::RequestPrimitives
    //   Polyhedron not available for cryoVacSolid1+cryoVacEndCap.
    //   This means it cannot be visualized on most systems.
    //   Contact the Visualization Coordinator.
    // So, I will avoid the union and place both parts of the vacuum separately
    //    G4UnionSolid* cryoVacSolid = new
    //    G4UnionSolid("cryoVacSolid1+cryoVacEndCap", cryoVacSolid1,
    //    cryoVacEndCap 0, cryoVacEndCapPos ); G4LogicalVolume* cryoVacLogic =
    //    new G4LogicalVolume(cryoVacSolid, Vacuum, "cryoVacLogic");

    G4LogicalVolume* cryoVacLogic = new G4LogicalVolume(cryoVacSolid1, Vacuum, "cryoVacLogic");
    G4LogicalVolume* cryoVacEndCapLogic = new G4LogicalVolume(cryoVacEndCapSolid, Vacuum, "cryoVacEndCapLogic");

    // push the cryostat in the world volume to keep the center of the TPC at
    // (0,0,0)
    cryoOuterPhys = new G4PVPlacement(0, cryostatShift, "cryoOuterPhys", cryoOuterLogic, fMotherVolume, false, 0, myCheckOverlap);
    G4PVPlacement* cryoVacPhys = new G4PVPlacement(0, cryostatShift, "cryoVacPhys", cryoVacLogic, fMotherVolume, false, 0, myCheckOverlap);
    G4PVPlacement* cryoVacEndCapPhys = new G4PVPlacement(0, cryostatShift + G4ThreeVector(0, 0, -cryoOuterH / 2.0), "cryoVacEndCapPhys", cryoVacEndCapLogic, fMotherVolume, false, 0, myCheckOverlap);
    Debug(cryoOuterPhys);
    Debug(cryoVacPhys);
    Debug(cryoVacEndCapPhys);

    //------------------------//
    //      Inner Cryostat    //
    //------------------------//

    //    G4double cryoInnerR = 131.0 / 2.0 * mm;  // initialized above
    G4double cryoInnerW = 0.7 * mm;
    G4double cryoInnerH = 597.0 * mm;
    G4double cryoInnerr = cryoInnerR - cryoInnerW;
    G4double cryoInnerPosZ = (cryoInnerH - cryoOuterH) / 2.0;  // no +cryoInnerBotH/2.0, the reference origin is cryoInnerSolid1
    G4ThreeVector cryoInnerPos(0, 0, cryoInnerPosZ);           // no +cryoInnerBotH/2.0, the reference origin is
                                                               // cryoInnerSolid1
    G4Tubs* cryoInnerSolid1 = new G4Tubs("cryoInnerSolid1", cryoInnerr, cryoInnerR, cryoInnerH / 2.0, 0, 2.0 * M_PI);
    G4UnionSolid* cryoInnerSolid = new G4UnionSolid("cryoInnerSolid1+cryoInnerBot", cryoInnerSolid1, cryoInnerBot, 0, G4ThreeVector(0, 0, -cryoInnerH / 2.0));
    G4LogicalVolume* cryoInnerLogic = new G4LogicalVolume(cryoInnerSolid, StainlessSteel, "cryoInnerLogic");
    cryoInnerLogic->SetVisAttributes(StainlessSteelVis);
    G4PVPlacement* cryoInnerPhys = new G4PVPlacement(0, cryostatShift + cryoInnerPos, "cryoInnerPhys", cryoInnerLogic, fMotherVolume, false, 0, myCheckOverlap);
    Debug(cryoInnerPhys);

    //-------------------------//
    //  Non scintillating LAr  //
    //-------------------------//

    G4double inactiveLArH = 300 * mm;  // guess by Michael, probable TODO
    G4Tubs* inactiveLArSolid = new G4Tubs("inactiveLArSolid", 0, cryoInnerr, inactiveLArH / 2.0, 0, 2 * M_PI);
    G4LogicalVolume* inactiveLArLogic = new G4LogicalVolume(inactiveLArSolid, NSLiquidArgon, "inactiveLArLogic");
    inactiveLArLogic->SetVisAttributes(NSLiquidArgonVis);
    //    G4ThreeVector inactiveLArPos(0, 0, (inactiveLArH-cryoInnerH)/2.0); //
    //    this was for placement in cryoInnerPhys
    G4ThreeVector inactiveLArPos = G4ThreeVector(0, 0, (inactiveLArH - cryoOuterH + cryoInnerBotH) / 2.0);
    fPhysicInactiveLAr = new G4PVPlacement(0, cryostatShift + inactiveLArPos, "inactiveLArPhys", inactiveLArLogic, fMotherVolume, false, 0, myCheckOverlap);
    Debug(fPhysicInactiveLAr);

    //-------------------------//
    //  Cryostat Skyline       //
    //-------------------------//

    G4double cryoConeInnerR = 254.5 / 2.0 * mm;
    G4double cryoConeW = 2.0 * mm;  // estimated, but I don't think it matters too much
    G4double cryoConeH = 100.0 * mm;
    G4double cryoConePosZ = cryoInnerPosZ + (cryoInnerH + cryoConeH) / 2.0;
    G4ThreeVector cryoConePos(0, 0, cryoConePosZ);
    G4Cons* cryoConeSolid = new G4Cons("cryoConeSolid", cryoInnerR - cryoConeW, cryoInnerR, cryoConeInnerR - cryoConeW, cryoConeInnerR, cryoConeH / 2.0, 0, 2.0 * M_PI);
    G4LogicalVolume* cryoConeLogic = new G4LogicalVolume(cryoConeSolid, StainlessSteel, "cryoConeLogic");
    cryoConeLogic->SetVisAttributes(StainlessSteelVis);
    G4PVPlacement* cryoConePhys = new G4PVPlacement(0, cryostatShift + cryoConePos, "cryoConePhys", cryoConeLogic, fMotherVolume, false, 0, myCheckOverlap);
    Debug(cryoConePhys);

    G4double cf250H = 26.0 * mm;
    G4double cf250OuterR = 304.8 / 2.0 * mm;
    G4double cf250InnerR = 0.0 * mm;  // I make them blind, I don't think there is a point in
                                      // detailed modelling thus far from the TPC
    G4double cf250aPosZ = cryoConePosZ + (cryoConeH + cf250H) / 2.0;
    G4double cf250bPosZ = cf250aPosZ + (cf250H + cf250H) / 2.0;  // weird, but I keep it this way
    G4ThreeVector cf250aPos(0, 0, cf250aPosZ);
    G4ThreeVector cf250bPos(0, 0, cf250bPosZ);
    G4Tubs* cf250Solid = new G4Tubs("cf250Solid", cf250InnerR, cf250OuterR, cf250H / 2.0, 0, 2.0 * M_PI);
    G4LogicalVolume* cf250Logic = new G4LogicalVolume(cf250Solid, StainlessSteel, "cf250Logic");
    cf250Logic->SetVisAttributes(StainlessSteelVis);
    G4PVPlacement* cf250aPhys = new G4PVPlacement(0, cryostatShift + cf250aPos, "cf250aPhys", cf250Logic, fMotherVolume, false, 0, myCheckOverlap);
    Debug(cf250aPhys);
    G4PVPlacement* cf250bPhys = new G4PVPlacement(0, cryostatShift + cf250bPos, "cf250bPhys", cf250Logic, fMotherVolume, false, 0, myCheckOverlap);
    Debug(cf250bPhys);

  } else {
    G4double CryoOuterR = 0.;
    G4double CryoOuterH = 0.;
    G4double CryoOuterThickness = 0.;
    G4double CryoInnerR = 0.;
    G4double CryoInnerH = 0.;
    G4double CryoInnerThickness = 0.;

    // Davide proposed classes to distiguish various configurations.
    // As quick-and-dirty patch I do it this way.
    if (cryostatFlag == 0) {         // Naples cryostat used for the first runs: ?????? - today
      CryoOuterR = 650. / 2 * mm;    // from Biagio
      CryoOuterH = 600. / 2 * mm;    // number invented by Michael
      CryoOuterThickness = 3. * mm;  // from Biagio
      CryoInnerR = 594. / 2 * mm;    // from Biagio
      CryoInnerH = 600. / 2 * mm;    // number invented by Michael
      CryoInnerThickness = 3. * mm;  // from Biagio
      G4double shiftToWall = 12.2 / 2 * cm;
      TPCshift = G4ThreeVector(-CryoInnerR + CryoInnerThickness + shiftToWall + 2. * mm, 0,
                               0);    // keep an arbitrary distance of 2 mm to the wall
    } else if (cryostatFlag == 99) {  // ARIS-like cryostat
      CryoOuterR = GetCm(6.10 / 2) * cm;
      CryoOuterH = GetCm(26.13 / 2) * cm;
      CryoOuterThickness = GetCm(0.048) * cm;
      CryoInnerR = GetCm(5.40 / 2) * cm;  // CryoOuterR - CryoOuterToInner; //should be 5.4/2 inches
      CryoInnerH = GetCm(23.13 / 2) * cm;
      CryoInnerThickness = GetCm(0.018) * cm;  // 0.13*2.54*cm;
    } else {
      DSLog(error) << "ReD cryostat configuration undefined: " << cryostatFlag << endlog;
      DSLog(fatal) << "Exiting ..." << endlog;
    }

    //------------------------//
    //      Outer Cryo        //
    //------------------------//

    G4Tubs* fSolidCryoOuter = new G4Tubs("CryoOuter_Solid", 0., CryoOuterR, CryoOuterH, 0., 2 * M_PI);

    G4LogicalVolume* fLogicCryoOuter = new G4LogicalVolume(fSolidCryoOuter, StainlessSteel, "OuterCryostat_Logic");
    fLogicCryoOuter->SetVisAttributes(StainlessSteelVis);

    // push the cryostat in the world volume to keep the center of the TPC at
    // (0,0,0)
    G4PVPlacement* fPhysicCryoOuter = new G4PVPlacement(0, -TPCshift, "OuterCryostat", fLogicCryoOuter, fMotherVolume, false, 0, myCheckOverlap);
    Debug(fPhysicCryoOuter);

    //--------------------------//
    //         Vacuum           //
    //--------------------------//

    G4Tubs* fSolidCryoVac = new G4Tubs("CryoVac_Solid", 0., CryoOuterR - CryoOuterThickness, CryoOuterH - CryoOuterThickness, 0., 2 * M_PI);

    G4LogicalVolume* fLogicCryoVac = new G4LogicalVolume(fSolidCryoVac, Vacuum, "CryostatVac_Logic");

    G4PVPlacement* fPhysicCryoVac = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), "CryostatVac", fLogicCryoVac, fPhysicCryoOuter, false, 0, myCheckOverlap);
    Debug(fPhysicCryoVac);

    //-------------------------//
    //      Inner Cryo         //
    //-------------------------//

    G4Tubs* fSolidCryoInner = new G4Tubs("CryoInner_Solid", 0., CryoInnerR, CryoInnerH, 0., 2 * M_PI);

    G4LogicalVolume* fLogicCryoInner = new G4LogicalVolume(fSolidCryoInner, StainlessSteel, "SolidCryoInner_Logic");
    fLogicCryoInner->SetVisAttributes(StainlessSteelVis);

    G4PVPlacement* fPhysicCryoInner = new G4PVPlacement(0, myZeros, "CryoInner", fLogicCryoInner, fPhysicCryoVac, false, 0, myCheckOverlap);
    Debug(fPhysicCryoInner);

    //-------------------------//
    //  Non scintillating LAr  //
    //-------------------------//

    G4Tubs* fSolidInactiveLAr = new G4Tubs("InactiveLAr_Solid", 0, CryoInnerR - CryoInnerThickness, (CryoInnerH - CryoInnerThickness), 0, 2 * M_PI);

    G4LogicalVolume* fLogicInactiveLAr = new G4LogicalVolume(fSolidInactiveLAr, NSLiquidArgon, "SolidInactiveLAr_Logic");
    fLogicInactiveLAr->SetVisAttributes(NSLiquidArgonVis);

    fPhysicInactiveLAr = new G4PVPlacement(0, myZeros, "InactiveLAr", fLogicInactiveLAr, fPhysicCryoInner, false, 0, myCheckOverlap);
    Debug(fPhysicInactiveLAr);
  }

  const int tpc = 1;
  // 0: old cylindrical TPC
  // 1: new UCLA quadratic TPC

  if (tpc == 0) {
    // ---------------------------------------------------//
    //           ReD cylindrical TPC dimensions           //
    // ---------------------------------------------------//

    // PTFE structure
    G4double PTFEOuterR = 7.4 / 2 * cm;
    G4double PTFEInnerR = 5.2 / 2 * cm;
    G4double PTFEWidth = 1.2 / 2 * cm;
    G4double PTFEOuterThickness = 0.3 * cm;
    G4double PTFETheta = extended_asin(PTFEWidth / PTFEOuterR);

    G4double PTFETop1OuterR = 12.2 / 2 * cm;
    G4double PTFETop1InnerR = 4.6 / 2 * cm;
    G4double PTFETop1H = 1. / 2 * cm;

    // define the helper volume "TPC"
    G4double TPCR = PTFETop1OuterR;  // fits precisely the large Teflon ring
    G4double TPCH = 200. / 2 * mm;   // approximate height

    G4double PTFETop2OuterR = 7.4 / 2 * cm;
    G4double PTFETop2InnerR = 4.6 / 2 * cm;
    G4double PTFETop2H = 0.8 / 2 * cm;

    G4double PTFETop3OuterR = 9.8 / 2 * cm;
    // G4double PTFETop3InnerR = 6.4 / 2 * cm;
    G4double PTFETop3H = 1.1 / 2 * cm;

    G4double TopSquareL = 6.4 / 2 * cm;
    G4double TopSquareH = PTFETop3H;

    G4double TopSIPML = 6.4 / 2 * cm;
    G4double TopSIPMH = 1. / 2 * mm;

    G4double PTFEBotOuterR = 9.8 / 2 * cm;
    // G4double PTFEBotInnerR = 6.4 / 2 * cm;
    G4double PTFEBotH = 1.45 / 2 * cm;

    G4double BotSquareL = 6.4 / 2 * cm;
    G4double BotSquareH = PTFEBotH;

    G4double BotSIPML = 6.4 / 2 * cm;
    G4double BotSIPMH = 1. / 2 * mm;

    // Copper rings
    G4double RingInnerR = PTFEOuterR;  //- PTFEOuterThickness;
    // Biagio writes the thickness is less than 0.2 mm
    G4double RingThickness = 0.2 * mm /* 0.1*cm */;
    G4double RingOuterR = RingInnerR + RingThickness;
    G4double RingLArOuterR = RingOuterR + PTFEOuterThickness;
    G4double RingLArInnerR = RingOuterR;

    // Biagio says there is no Teflon foil of 3 mm thickness, but a "thin"
    // reflector foil. I keep it Teflon, reduce it to 0.1 mm -> push it inwards
    // by 2.9 mm.
    const G4double reducedWidth = 2.9 * mm;
    G4double TeflonFoilOuterR = PTFEInnerR - reducedWidth;
    G4double TeflonThickness = 0.3 * cm - reducedWidth;
    G4double TeflonFoilInnerR = TeflonFoilOuterR - TeflonThickness;

    // TPBLayer
    G4double TPBLayerOuterR = TeflonFoilInnerR;
    G4double TPBLayerThickness = 0.1 * mm;
    G4double TPBLayerInnerR = TPBLayerOuterR - TPBLayerThickness;

    G4double TPCLArR = TPBLayerInnerR;

    G4double MainStructureH = 12.47 / 2 * cm;

    G4double PTFEOuterH = MainStructureH;
    // G4double PTFEInnerH = MainStructureH;
    G4double TeflonFoilH = MainStructureH;
    G4double TPBLayerH = MainStructureH;
    G4double TPCLArH = MainStructureH;
    G4double TPCGArH = 0.1 * cm;  // to be checked

    // G4double GridThickness = 0.001 * cm;  ///GetCm(0.1/2)*cm; //TBC

    G4double RingH = 0.78 / 2 * cm;
    G4double RingSpace = 0.89 * cm;
    // G4double RingSpaceFirst = RingH+RingSpace;

    DSStorage::Get()->SetLArGArBoundaryPosZ(TPCLArH - TPCGArH);
    DSLog(routine) << "TPC " << tpc << " gas pocket position: " << DSStorage::Get()->GetLArGArBoundaryPosZ() / mm << " TPCLArR " << TPCLArR / mm << " TPCLArH " << TPCLArH / mm << " TPCGArH " << TPCGArH / mm << " mm" << endlog;

    //--------------------------//
    //  Non-central TPC volume  //
    //--------------------------//
    // dummy volume to ease shift of the entire TPC

    G4Tubs* fSolidTPC = new G4Tubs("TPC_Solid", 0, TPCR, TPCH, 0, 2 * M_PI);

    G4LogicalVolume* fLogicTPC = new G4LogicalVolume(fSolidTPC, NSLiquidArgon, "TPC_Logic");
    fLogicTPC->SetVisAttributes(NSLiquidArgonVis);
    fLogicTPC->SetVisAttributes(false);

    // G4PVPlacement* fPhysicTPC =
    new G4PVPlacement(0, TPCshift, "TPC", fLogicTPC, fPhysicInactiveLAr, false, 0, myCheckOverlap);
    Debug(fPhysicTPC);

    //-----------------------//
    //    PTFE Structure     //
    //-----------------------//

    G4Tubs* fSolidExternalPTFE1 = new G4Tubs("ExternalPTFE1_Solid", PTFEInnerR, PTFEOuterR, PTFEOuterH, 45 * deg, PTFETheta);
    G4LogicalVolume* fLogicExternalPTFE1 = new G4LogicalVolume(fSolidExternalPTFE1, Teflon, "SolidExternalPTFE1_Logic");
    G4PVPlacement* fPTFEStructure1 = new G4PVPlacement(0, myZeros, "PTFEStructure_1", fLogicExternalPTFE1, fPhysicTPC, false, 0, myCheckOverlap);
    Debug(fPTFEStructure1);

    G4Tubs* fSolidExternalPTFE2 = new G4Tubs("ExternalPTFE2_Solid", PTFEInnerR, PTFEOuterR, PTFEOuterH, 165 * deg, PTFETheta);
    G4LogicalVolume* fLogicExternalPTFE2 = new G4LogicalVolume(fSolidExternalPTFE2, Teflon, "SolidExternalPTFE2_Logic");
    G4PVPlacement* fPTFEStructure2 = new G4PVPlacement(0, myZeros, "PTFEStructure_2", fLogicExternalPTFE2, fPhysicTPC, false, 0, myCheckOverlap);
    Debug(fPTFEStructure2);

    G4Tubs* fSolidExternalPTFE3 = new G4Tubs("ExternalPTFE3_Solid", PTFEInnerR, PTFEOuterR, PTFEOuterH, 285 * deg, PTFETheta);
    G4LogicalVolume* fLogicExternalPTFE3 = new G4LogicalVolume(fSolidExternalPTFE3, Teflon, "SolidExternalPTFE3_Logic");
    G4PVPlacement* fPTFEStructure3 = new G4PVPlacement(0, myZeros, "PTFEStructure_3", fLogicExternalPTFE3, fPhysicTPC, false, 0, myCheckOverlap);
    Debug(fPTFEStructure3);

    /////////////////////////////
    ///   Top and Bottom PTFE ///
    /////////////////////////////

    G4Tubs* fSolidTopPTFE1 = new G4Tubs("TopPTFE1_Solid", 0, PTFETop1OuterR, PTFETop1H, 0, 2 * M_PI);
    G4LogicalVolume* fLogicTopPTFE1 = new G4LogicalVolume(fSolidTopPTFE1, Teflon, "TopPTFE1_Logic");
    G4ThreeVector top1(0, 0, PTFEOuterH + PTFETop1H);
    G4PVPlacement* fPhysicTopPTFE1 = new G4PVPlacement(0, top1, "TopPTFE_1", fLogicTopPTFE1, fPhysicTPC, false, 0, myCheckOverlap);
    Debug(fPhysicTopPTFE1);

    ////////////////////////////////////////////////////////////////////////////
    ///gas pocket
    G4Tubs* fSolidTopInnerGAr = new G4Tubs("TopInnerGAr_Solid", 0, PTFETop1InnerR, PTFETop1H, 0, 2 * M_PI);
    G4LogicalVolume* fLogicTopGAr = new G4LogicalVolume(fSolidTopInnerGAr, GaseousArgon, "TopInnerGAr_Logic");
    fLogicTopGAr->SetVisAttributes(GaseousArgonVis);
    G4PVPlacement* fPhysicTopGAr1 = new G4PVPlacement(0, myZeros, "TopGAr_1", fLogicTopGAr, fPhysicTopPTFE1, false, 0, myCheckOverlap);
    Debug(fPhysicTopGAr1);
    ////////////////////////////////////////////////////////////////////////////
    ///gas pocket 1

    G4Tubs* fSolidTopPTFE2 = new G4Tubs("TopPTFE2_Solid", 0, PTFETop2OuterR, PTFETop2H, 0, 2 * M_PI);
    G4LogicalVolume* fLogicTopPTFE2 = new G4LogicalVolume(fSolidTopPTFE2, Teflon, "TopPTFE2_Logic");
    G4ThreeVector top2(0, 0, PTFEOuterH + PTFETop1H * 2 + PTFETop2H);
    G4PVPlacement* fPhysicTopPTFE2 = new G4PVPlacement(0, top2, "TopPTFE_2", fLogicTopPTFE2, fPhysicTPC, false, 0, myCheckOverlap);
    Debug(fPhysicTopPTFE2);

    ////////////////////////////////////////////////////////////////////////////
    ///gas pocket 2
    G4Tubs* fSolidTopGAr2 = new G4Tubs("TopGAr2_Solid", 0, PTFETop2InnerR, PTFETop2H, 0, 2 * M_PI);
    G4LogicalVolume* fLogicTopGAr2 = new G4LogicalVolume(fSolidTopGAr2, GaseousArgon, "TopGAr2_Logic");
    fLogicTopGAr2->SetVisAttributes(GaseousArgonVis);
    G4PVPlacement* fPhysicTopGAr2 = new G4PVPlacement(0, myZeros, "TopGAr_2", fLogicTopGAr2, fPhysicTopPTFE2, false, 0, myCheckOverlap);
    Debug(fPhysicTopGAr2);
    ////////////////////////////////////////////////////////////////////////////
    ///gas pocket 2

    G4Tubs* fSolidTopPTFE3 = new G4Tubs("TopPTFE3_Solid", 0, PTFETop3OuterR, PTFETop3H, 0, 2 * M_PI);
    G4LogicalVolume* fLogicTopPTFE3 = new G4LogicalVolume(fSolidTopPTFE3, Teflon, "TopPTFE_Logic");
    G4ThreeVector top3(0, 0, PTFEOuterH + PTFETop1H * 2 + PTFETop2H * 2 + PTFETop3H);
    G4PVPlacement* fPhysicTopPTFE = new G4PVPlacement(0, top3, "TopPTFE_3", fLogicTopPTFE3, fPhysicTPC, false, 0, myCheckOverlap);
    Debug(fPhysicTopPTFE);

    G4Box* fSolidTopSquarePTFE = new G4Box("TopPTFESquare_Solid", TopSquareL, TopSquareL, TopSquareH);
    G4LogicalVolume* fLogicTopSquarePTFE = new G4LogicalVolume(fSolidTopSquarePTFE, NSLiquidArgon, "TopSquarePTFE_Logic");
    fLogicTopSquarePTFE->SetVisAttributes(NSLiquidArgonVis);
    G4PVPlacement* fPhysicTopSquarePTFE = new G4PVPlacement(0, myZeros, "TopSquarePTFE", fLogicTopSquarePTFE, fPhysicTopPTFE, false, 0, myCheckOverlap);
    Debug(fPhysicTopSquarePTFE);

    G4Box* fSolidTopSIPM = new G4Box("TopSIPM_Solid", TopSIPML, TopSIPML, TopSIPMH);
    G4LogicalVolume* fLogicSiPMTop = new G4LogicalVolume(fSolidTopSIPM, MetalSilicon, "TopSIPM_Logic");
    G4ThreeVector top_sipm(0, 0, -PTFETop3H + TopSIPMH);
    G4PVPlacement* fTopSIPM = new G4PVPlacement(0, top_sipm, "TopSIPM", fLogicSiPMTop, fPhysicTopSquarePTFE, false, 0, myCheckOverlap);
    Debug(fTopSIPM);

    G4Tubs* fSolidBotPTFE = new G4Tubs("BotPTFE_Solid", 0, PTFEBotOuterR, PTFEBotH, 0, 2 * M_PI);
    G4LogicalVolume* fLogicBotPTFE = new G4LogicalVolume(fSolidBotPTFE, Teflon, "BotPTFE_Logic");
    G4ThreeVector bot1(0, 0, -PTFEOuterH - PTFEBotH);
    G4PVPlacement* fPhysicBotPTFE = new G4PVPlacement(0, bot1, "BotPTFE", fLogicBotPTFE, fPhysicTPC, false, 0, myCheckOverlap);
    Debug(fPhysicBotPTFE);

    G4Box* fSolidBotSquarePTFE = new G4Box("BotPTFESquare_Solid", BotSquareL, BotSquareL, BotSquareH);
    G4LogicalVolume* fLogicBotSquarePTFE = new G4LogicalVolume(fSolidBotSquarePTFE, NSLiquidArgon, "BotSquarePTFE_Logic");
    fLogicBotSquarePTFE->SetVisAttributes(NSLiquidArgonVis);
    G4PVPlacement* fPhysicBotSquarePTFE = new G4PVPlacement(0, myZeros, "BotSquarePTFE", fLogicBotSquarePTFE, fPhysicBotPTFE, false, 0, myCheckOverlap);
    Debug(fPhysicBotSquarePTFE);

    G4Box* fSolidBotSIPM = new G4Box("BotSIPM_Solid", BotSIPML, BotSIPML, BotSIPMH);
    G4LogicalVolume* fLogicSiPMBottom = new G4LogicalVolume(fSolidBotSIPM, MetalSilicon, "BotSIPM_Logic");
    G4ThreeVector bot_sipm(0, 0, PTFEBotH - BotSIPMH);
    G4PVPlacement* fBotSIPM = new G4PVPlacement(0, bot_sipm, "BotSIPM", fLogicSiPMBottom, fPhysicBotSquarePTFE, false, 0, myCheckOverlap);
    Debug(fBotSIPM);

    //-----------------------//
    //    Copper Rings       //
    //-----------------------//

    G4Tubs* fSolidCopperRing = new G4Tubs("CopperRing_Solid", RingInnerR, RingOuterR, RingH, 0, 2 * M_PI);
    G4LogicalVolume* fLogicCopperRing = new G4LogicalVolume(fSolidCopperRing, MetalCopper, "CopperRing_Logic");
    fLogicCopperRing->SetVisAttributes(MetalCopperVis);

    G4Tubs* fSolidLArRing = new G4Tubs("LArRing_Solid", RingLArInnerR, RingLArOuterR, RingH, 0, 2 * M_PI);
    G4LogicalVolume* fLogicLArRing = new G4LogicalVolume(fSolidLArRing, NSLiquidArgon, "LArRing_Logic");
    fLogicLArRing->SetVisAttributes(NSLiquidArgonVis);

    // define 8 rings
    G4ThreeVector RingPlacement(0., 0., 0.);
    char ring_name[256];
    char larring_name[256];
    G4PVPlacement* fRingCopper[8];
    G4PVPlacement* fRingLAr[8];

    for (int j = 0; j < 8; j++) {
      snprintf(ring_name,90, "CopperRing_%d", j);
      snprintf(larring_name,90, "LArRing_%d", j);
      RingPlacement.setZ(PTFEOuterH - RingH - (RingSpace + 2 * RingH) * j);

      fRingLAr[j] = new G4PVPlacement(0, RingPlacement, larring_name, fLogicLArRing, fPhysicTPC, false, 0, myCheckOverlap);
      Debug(fRingLAr[j]);

      fRingCopper[j] = new G4PVPlacement(0, RingPlacement, ring_name, fLogicCopperRing, fPhysicTPC, false, 0, myCheckOverlap);
      Debug(fRingCopper[j]);
    }

    //-----------------------//
    //    Teflon Foil        //
    //-----------------------//

    G4Tubs* fSolidTeflonFoil = new G4Tubs("TeflonFoil_Solid", 0, TeflonFoilOuterR, TeflonFoilH, 0, 2 * M_PI);
    G4LogicalVolume* fLogicTeflonFoil = new G4LogicalVolume(fSolidTeflonFoil, Teflon, "TeflonFoil_Logic");
    // G4PVPlacement* fPhysicTeflonFoil =
    new G4PVPlacement(0, myZeros, "TeflonFoil", fLogicTeflonFoil, fPhysicTPC, false, 0, myCheckOverlap);
    Debug(fPhysicTeflonFoil);

    //-----------------------//
    //    TPB Layer          //
    //-----------------------//

    G4Tubs* fSolidTPBLayer = new G4Tubs("TPBLayer_Solid", 0, TPBLayerOuterR, TPBLayerH, 0, 2 * M_PI);
    G4LogicalVolume* fLogicTPBLayer = new G4LogicalVolume(fSolidTPBLayer, TPB, "TPBLayer_Logic");
    G4PVPlacement* fPhysicTPBLayer = new G4PVPlacement(0, myZeros, "TPBLayer", fLogicTPBLayer, fPhysicTeflonFoil, false, 0, myCheckOverlap);
    Debug(fPhysicTPBLayer);

    //-----------------------------//
    //    LAr TPC (Active_LAr)     //
    //-----------------------------//

    G4Tubs* fSolidActiveLArTPC = new G4Tubs("ActiveLArTPC_Solid", 0, TPCLArR, TPCLArH, 0, 2 * M_PI);
    G4LogicalVolume* fLogicActiveLArTPC = new G4LogicalVolume(fSolidActiveLArTPC, LiquidArgon, "LAr_Logic");
    fLogicActiveLArTPC->SetVisAttributes(LiquidArgonVis);
    fActiveLArTPC = new G4PVPlacement(0, myZeros, "ActiveLAr", fLogicActiveLArTPC, fPhysicTPBLayer, false, 0, myCheckOverlap);
    Debug(fActiveLArTPC);

    //-----------------------//
    //      Gas Pocket       //
    //-----------------------//

    G4Tubs* fSolidGasPocket = new G4Tubs("GasPocket_Solid", 0, TPCLArR, TPCGArH, 0, 2 * M_PI);
    G4LogicalVolume* fLogicGasPocket = new G4LogicalVolume(fSolidGasPocket, GaseousArgon, "GasPocket_Logic");
    fLogicGasPocket->SetVisAttributes(GaseousArgonVis);
    G4ThreeVector gaspocket(0, 0, TPCLArH - TPCGArH);
    // G4PVPlacement* fPhysicGasPocket =
    new G4PVPlacement(0, gaspocket, "GasPocket", fLogicGasPocket, fActiveLArTPC, false, 0, myCheckOverlap);
    Debug(fPhysicGasPocket);

    //----------------------------------------------------------------------------

    G4Region* fLArRegion = new G4Region("LAr_Logic");
    fLogicActiveLArTPC->SetRegion(fLArRegion);
    fLArRegion->AddRootLogicalVolume(fLogicActiveLArTPC);

    DefineSurfaces();

    // Set the z coordinate of the LAr - GAr interface, necessary for S2
    // generation in DSLightX
    // DSStorage::Get()->SetLArGArBoundaryPosZ( LArGarIntefaceZ  + 1.0*um);

    // make SiPM as pe storing material
    DSStorage::Get()->SetPMTMaterialIndex(fLogicSiPMBottom->GetMaterial()->GetIndex());

    // DSStorage::Get()->SetPMTMaterialIndex(fLogicBottomPMTWindow->GetMaterial()->GetIndex());

  } else if (tpc == 1) {

    // ---------------------------------------------------//
    //           ReD quadratic TPC / UCLA                 //
    // ---------------------------------------------------//

    G4double coatingITOH = 25.0 * nm;
    G4double coatingTPBH = 1700.0 * nm;
    //    coatingITOH = 0.2 * mm;
    //    coatingTPBH = 0.15 * mm;
    G4double coatingITOH2 = coatingITOH / 2.0;
    G4double coatingTPBH2 = coatingTPBH / 2.0;

    //-----------------------//
    // TPC Container Volume  //
    //-----------------------//
    // NOW: the new cryostat has an outer radius of 65.5 mm and wall thickness
    // 0.7 * mm. Thus, a 92 mm x 92 mm x ~100 mm box won't fit. I HAVE TO FILLET
    // THE CORNERS!
    G4double TPCL2 = 46.1 * mm;  // real extent +- 46 mm
    G4double TPCH2 = 51.2 * mm;  // real extent -41.025 - +51.075 mm
    G4double filleterOuterR = TPCL2 * std::sqrt(2.0) + 0.1 * mm;
    G4double filleterInnerR = 65.5 * mm - 0.7 * mm - 0.1 * mm;
    G4Box* TPCboxSolid = new G4Box("TPCboxSolid", TPCL2, TPCL2, TPCH2);
    G4Tubs* TPCtubsSolid = new G4Tubs("TPCtubsSolid", filleterInnerR, filleterOuterR, TPCH2, 0, 2.0 * M_PI);
    G4SubtractionSolid* TPCSolid = new G4SubtractionSolid("TPCSolid", TPCboxSolid, TPCtubsSolid);
    G4LogicalVolume* TPCLogic = new G4LogicalVolume(TPCSolid, NSLiquidArgon, "TPCLogic");
    TPCLogic->SetVisAttributes(NSLiquidArgonVis);
    TPCLogic->SetVisAttributes(false);
    G4RotationMatrix TPCrot;
    //    TPCrot.rotateZ(22.5*deg);
    const G4Transform3D TPCtransform(TPCrot, TPCshift);
    // the object fPhysicInactiveLAr is created in the (future class) cryostat
    // G4PVPlacement* fPhysicTPC =
    new G4PVPlacement(TPCtransform, "TPC", TPCLogic, fPhysicInactiveLAr, false, 0, myCheckOverlap);
    Debug(fPhysicTPC);

    //----------------------------------------------------------------------------------//
    // Reflection_drift_side: // model the three volumes as staggered boxes,
    // with a final forth one of active LAr //
    //----------------------------------------------------------------------------------//
    //----------------------------------------//
    // Reflection_drift_side_back             //
    //----------------------------------------//
    G4double Reflection_drift_side_backOuterL2 = 28.05 * mm;
    G4double Reflection_drift_side_backH2 = 50.0 / 2.0 * mm;
    //    G4double Reflection_drift_side_backX       =  0.0       * mm;
    //    G4double Reflection_drift_side_backY       =  0.0       * mm;
    //    G4double Reflection_drift_side_backZ       =  0.0       * mm;
    //    G4ThreeVector
    //    Reflection_drift_side_backPos(Reflection_drift_side_backX,
    //    Reflection_drift_side_backY, Reflection_drift_side_backZ);
    G4Box* Reflection_drift_side_backSolid = new G4Box("Reflection_drift_side_backSolid", Reflection_drift_side_backOuterL2, Reflection_drift_side_backOuterL2, Reflection_drift_side_backH2);
    G4LogicalVolume* Reflection_drift_side_backLogic = new G4LogicalVolume(Reflection_drift_side_backSolid, Acrylic, "Reflection_drift_side_backLogic");
    G4PVPlacement* Reflection_drift_side_backPhys = new G4PVPlacement(0, myZeros, "Reflection_drift_side_backPhys", Reflection_drift_side_backLogic, fPhysicTPC, false, 0, myCheckOverlap);
    Debug(Reflection_drift_side_backPhys);
    //----------------------------------------//
    // Reflection_drift_side_ESR              //
    //----------------------------------------//
    G4double Reflection_drift_side_ESROuterL2 = 26.55 * mm;
    G4double Reflection_drift_side_ESRH2 = 50.0 / 2.0 * mm;
    G4Box* Reflection_drift_side_ESRSolid = new G4Box("Reflection_drift_side_ESRSolid", Reflection_drift_side_ESROuterL2, Reflection_drift_side_ESROuterL2, Reflection_drift_side_ESRH2);
    G4LogicalVolume* Reflection_drift_side_ESRLogic = new G4LogicalVolume(Reflection_drift_side_ESRSolid, ThreeMFoil, "Reflection_drift_side_ESRLogic");
    G4PVPlacement* Reflection_drift_side_ESRPhys = new G4PVPlacement(0, myZeros, "Reflection_drift_side_ESRPhys", Reflection_drift_side_ESRLogic, Reflection_drift_side_backPhys, false, 0, myCheckOverlap);
    Debug(Reflection_drift_side_ESRPhys);
    //----------------------------------------//
    // Reflection_drift_side_front            //
    //----------------------------------------//
    G4double Reflection_drift_side_frontOuterL2 = 26.5 * mm;
    G4double Reflection_drift_side_frontInnerL2 = 25.0 * mm;
    G4double Reflection_drift_side_frontH2 = 50.0 / 2.0 * mm;
    G4Box* Reflection_drift_side_frontSolid = new G4Box("Reflection_drift_side_frontSolid", Reflection_drift_side_frontOuterL2, Reflection_drift_side_frontOuterL2, Reflection_drift_side_frontH2);
    G4LogicalVolume* Reflection_drift_side_frontLogic = new G4LogicalVolume(Reflection_drift_side_frontSolid, Acrylic, "Reflection_drift_side_frontLogic");
    G4PVPlacement* Reflection_drift_side_frontPhys = new G4PVPlacement(0, myZeros, "Reflection_drift_side_frontPhys", Reflection_drift_side_frontLogic, Reflection_drift_side_ESRPhys, false, 0, myCheckOverlap);
    Debug(Reflection_drift_side_frontPhys);
    //----------------------------------------//
    // Reflection_drift_side_front TPB        //
    //----------------------------------------//
    G4double Reflection_drift_side_frontTPBOuterL2 = Reflection_drift_side_frontInnerL2;
    G4double Reflection_drift_side_frontTPBInnerL2 = Reflection_drift_side_frontTPBOuterL2 - coatingTPBH;  // needed for placement of activeLAr
    G4double Reflection_drift_side_frontTPBH2 = Reflection_drift_side_frontH2;
    G4Box* Reflection_drift_side_frontTPBSolid = new G4Box("Reflection_drift_side_frontTPBSolid", Reflection_drift_side_frontTPBOuterL2, Reflection_drift_side_frontTPBOuterL2, Reflection_drift_side_frontTPBH2);
    G4LogicalVolume* Reflection_drift_side_frontTPBLogic = new G4LogicalVolume(Reflection_drift_side_frontTPBSolid, TPB, "Reflection_drift_side_frontTPBLogic");
    G4PVPlacement* Reflection_drift_side_frontTPBPhys = new G4PVPlacement(0, myZeros, "Reflection_drift_side_frontTPBPhys", Reflection_drift_side_frontTPBLogic, Reflection_drift_side_frontPhys, false, 0, myCheckOverlap);
    Debug(Reflection_drift_side_frontTPBPhys);

    //----------------------------------------//
    // Field_cage_support_ex                  // simplified, a probable TODO
    //----------------------------------------//
    G4double Field_cage_support_exLX2 = 6.25 * mm;
    G4double Field_cage_support_exLY2 = 11.992 / 2.0 * mm;
    G4double Field_cage_support_exLZ2 = 55.008 / 2.0 * mm;
    G4double Field_cage_support_exX = 0.0 * mm;
    G4double Field_cage_support_exY = 39.004 * mm;
    G4double Field_cage_support_exZ = -4.028 * mm;
    G4Box* Field_cage_support_exFullSolid = new G4Box("Field_cage_support_exFullSolid", Field_cage_support_exLX2, Field_cage_support_exLY2, Field_cage_support_exLZ2);
    // the little box to be cut off: if I make it double larger than needed,
    // it's center will be at 0,LY2,LZ2!
    G4double Field_cage_support_exCutLX2 = 12.5 * mm;
    G4double Field_cage_support_exCutLY2 = 4.992 * mm;
    G4double Field_cage_support_exCutLZ2 = 8.564 * mm;
    G4Box* Field_cage_support_exCutSolid = new G4Box("Field_cage_support_exCutSolid", Field_cage_support_exCutLX2, Field_cage_support_exCutLY2, Field_cage_support_exCutLZ2);
    G4ThreeVector Field_cage_support_exCutPos(0, -Field_cage_support_exLY2, -Field_cage_support_exLZ2);
    G4SubtractionSolid* Field_cage_support_exSolid = new G4SubtractionSolid("Field_cage_support_exSolid", Field_cage_support_exFullSolid, Field_cage_support_exCutSolid, 0, Field_cage_support_exCutPos);
    G4LogicalVolume* Field_cage_support_exLogic = new G4LogicalVolume(Field_cage_support_exSolid, Teflon, "Field_cage_support_exLogic");
    G4PVPlacement* Field_cage_support_exPhys[4];
    for (int i = 0; i < 4; ++i) {  // I should learn to use repeated volumes
      char name[256];
      snprintf(name,90, "Field_cage_support_exPhys%d", i);
      G4double angle = i * 90.0 * deg;
      G4RotationMatrix rotZ;
      G4ThreeVector pos(Field_cage_support_exX, Field_cage_support_exY, Field_cage_support_exZ);
      rotZ.rotateZ(angle);
      pos.rotateZ(angle);
      const G4Transform3D transform(rotZ, pos);
      Field_cage_support_exPhys[i] = new G4PVPlacement(transform, name, Field_cage_support_exLogic, fPhysicTPC, false, 0, myCheckOverlap);
      Debug(Field_cage_support_exPhys[i]);
    }

    //----------------------------------------//
    // Field_cage_shaping_ring                //
    //----------------------------------------//
    G4double Field_cage_shaping_ringOuterL2 = 33.008 * mm;
    G4double Field_cage_shaping_ringD2 = 0.508 / 2.0 * mm;
    G4double Field_cage_shaping_ringOuterR = 3.008 * mm;
    // G4double Field_cage_shaping_ringInnerR = 2.5 * mm;
    G4double Field_cage_shaping_ringH2 = 4.0 / 2.0 * mm;
    G4double Field_cage_shaping_ringBarX = Field_cage_shaping_ringOuterL2 - Field_cage_shaping_ringD2;
    G4UnionSolid* Field_cage_shaping_ringSolid = G4QuadRingWithFillet(Field_cage_shaping_ringBarX, Field_cage_shaping_ringH2, Field_cage_shaping_ringD2, Field_cage_shaping_ringOuterR);
    G4LogicalVolume* Field_cage_shaping_ringLogic = new G4LogicalVolume(Field_cage_shaping_ringSolid, MetalCopper, "Field_cage_shaping_ringLogic");
    Field_cage_shaping_ringLogic->SetVisAttributes(MetalCopperVis);
    G4PVPlacement* Field_cage_shaping_ringPhys[9];
    for (int i = 0; i < 9; ++i) {  // I should learn to use repeated volumes
      char name[256];
      snprintf(name,90, "Field_cage_shaping_ringPhys%d", i);
      G4ThreeVector pos(Field_cage_shaping_ringBarX, 0, i * 5.0 * mm - 20.0 * mm);
      Field_cage_shaping_ringPhys[i] = new G4PVPlacement(0, pos, name, Field_cage_shaping_ringLogic, fPhysicTPC, false, 0, myCheckOverlap);
      Debug(Field_cage_shaping_ringPhys[i]);
    }

    //----------------------------------------//
    // Field_cage_support_in                  //HURZ
    //----------------------------------------//

    G4double Field_cage_support_inLX2 = 8.0 / 2.0 * mm;
    G4double Field_cage_support_inLY2 = 6.35 * mm;
    G4double Field_cage_support_inLZ2 = 25.0 * mm;
    G4double Field_cage_support_inX = 28.05 * mm;
    G4double Field_cage_support_inY = 28.05 * mm;
    G4double Field_cage_support_inZ = 0.0 * mm;
    G4double Field_cage_support_inCutBoxL2 = 6.0 * mm;  // 6 mm is sufficient, because all cutted edges have length
                                                        // sqrt(2)*4 mm
    // G4double Field_cage_support_inThickPlateH = 3.0 * mm;
    // G4double Field_cage_support_inThinPlateH = 1.0 * mm;
    //  ******** the "full" strut
    G4Box* Field_cage_support_inFullSolid = new G4Box("Field_cage_support_inFullSolid", Field_cage_support_inLX2, Field_cage_support_inLY2, Field_cage_support_inLZ2);
    // ******** the inner cut (full height)
    G4Box* Field_cage_support_inInnerCutSolid = new G4Box("Field_cage_support_inInnerCutSolid", Field_cage_support_inCutBoxL2, Field_cage_support_inCutBoxL2, 1.01 * Field_cage_support_inLZ2);
    // to be rotated by 45 deg, with a corner in the center of the full volume.
    G4RotationMatrix Field_cage_support_inRot;
    Field_cage_support_inRot.rotateZ(45.0 * deg);
    G4ThreeVector Field_cage_support_inPos(-Field_cage_support_inCutBoxL2 * std::sqrt(2.0), 0, 0);
    G4Transform3D Field_cage_support_inTrans(Field_cage_support_inRot, Field_cage_support_inPos);
    G4SubtractionSolid* Field_cage_support_inCutSolid = new G4SubtractionSolid("Field_cage_support_inCutSolid", Field_cage_support_inFullSolid, Field_cage_support_inInnerCutSolid, Field_cage_support_inTrans);
    // ******** the (two) outer cuts, not full height but between the two end
    // plates, cutting away the 8 intermediate plates too. Change: I noticed
    // that there is a collision between the Field_cage_support_in struts and
    // Field_cage_shaping_ring_thick_copper001 (the higher of the two). Thus,
    // right now I'll cut on the outer side completely.  Actually, I could use
    // the same solid for inner and outer cutting.
    //    G4Box* Field_cage_support_inOuterCutSolid = new
    //    G4Box("Field_cage_support_inOuterCutSolid",
    //    Field_cage_support_inCutBoxL2, Field_cage_support_inCutBoxL2,
    //    Field_cage_support_inLZ2-Field_cage_support_inThickPlateH);
    G4Box* Field_cage_support_inOuterCutSolid = Field_cage_support_inInnerCutSolid;
    // the extra shift of 4.45 mm accounts for the gap between
    // Reflection_drift_side_back and Field_cage_shaping_ring
    G4double Field_cage_support_inShift = (4.45 * mm + Field_cage_support_inCutBoxL2) / std::sqrt(2.0);
    Field_cage_support_inPos = G4ThreeVector(Field_cage_support_inShift, Field_cage_support_inShift, 0);
    Field_cage_support_inTrans = G4Transform3D(Field_cage_support_inRot, Field_cage_support_inPos);
    G4SubtractionSolid* Field_cage_support_inCutCutSolid = new G4SubtractionSolid("Field_cage_support_inCutCutSolid", Field_cage_support_inCutSolid, Field_cage_support_inOuterCutSolid, Field_cage_support_inTrans);
    Field_cage_support_inPos = G4ThreeVector(Field_cage_support_inShift, -Field_cage_support_inShift, 0);
    Field_cage_support_inTrans = G4Transform3D(Field_cage_support_inRot, Field_cage_support_inPos);
    G4SubtractionSolid* Field_cage_support_inSolid = new G4SubtractionSolid("Field_cage_support_inSolid", Field_cage_support_inCutCutSolid, Field_cage_support_inOuterCutSolid, Field_cage_support_inTrans);
    G4LogicalVolume* Field_cage_support_inLogic = new G4LogicalVolume(Field_cage_support_inSolid, Teflon, "Field_cage_support_inLogic");
    // There are still the little spacer plates missing.  TODO!
    G4PVPlacement* Field_cage_support_inPhys[4];
    for (int i = 0; i < 4; ++i) {  // I should learn to use repeated volumes
      char name[256];
      snprintf(name,90, "Field_cage_support_inPhys%d", i);
      G4RotationMatrix rot;
      rot.rotateZ((45.0 + i * 90.0) * deg);
      G4ThreeVector pos(Field_cage_support_inX, Field_cage_support_inY, Field_cage_support_inZ);
      pos.rotateZ(i * 90.0 * deg);
      G4Transform3D trans(rot, pos);
      Field_cage_support_inPhys[i] = new G4PVPlacement(trans, name, Field_cage_support_inLogic, fPhysicTPC, false, 0, myCheckOverlap);
      Debug(Field_cage_support_inPhys[i]);
    }

    //----------------------------------------//
    // Field_cage_shaping_ring_thin_ss        //
    //----------------------------------------//
    G4double Field_cage_shaping_ring_thin_ssOuterL2 = 46.0 * mm;
    G4double Field_cage_shaping_ring_thin_ssD2 = 10.0 / 2.0 * mm;
    G4double Field_cage_shaping_ring_thin_ssOuterR = 13.5 * mm;
    // G4double Field_cage_shaping_ring_thin_ssInnerR = 3.5 * mm;
    G4double Field_cage_shaping_ring_thin_ssH2 = 1.524 / 2.0 * mm;
    G4double Field_cage_shaping_ring_thin_ssZ[2] = {24.238 * mm, 25.812 * mm};
    G4double Field_cage_shaping_ring_thin_ssBarX = Field_cage_shaping_ring_thin_ssOuterL2 - Field_cage_shaping_ring_thin_ssD2;
    G4UnionSolid* Field_cage_shaping_ring_thin_ssSolid = G4QuadRingWithFillet(Field_cage_shaping_ring_thin_ssBarX, Field_cage_shaping_ring_thin_ssH2, Field_cage_shaping_ring_thin_ssD2, Field_cage_shaping_ring_thin_ssOuterR);
    G4LogicalVolume* Field_cage_shaping_ring_thin_ssLogic = new G4LogicalVolume(Field_cage_shaping_ring_thin_ssSolid, StainlessSteel, "Field_cage_shaping_ring_thin_ssLogic");
    Field_cage_shaping_ring_thin_ssLogic->SetVisAttributes(StainlessSteelVis);
    G4PVPlacement* Field_cage_shaping_ring_thin_ssPhys[2];
    for (int i = 0; i < 2; ++i) {  // I should learn to use repeated volumes
      char name[256];
      snprintf(name, 90,"Field_cage_shaping_ring_thin_ssPhys%d", i);
      G4ThreeVector pos(Field_cage_shaping_ring_thin_ssBarX, 0, Field_cage_shaping_ring_thin_ssZ[i]);
      Field_cage_shaping_ring_thin_ssPhys[i] = new G4PVPlacement(0, pos, name, Field_cage_shaping_ring_thin_ssLogic, fPhysicTPC, false, 0, myCheckOverlap);
      Debug(Field_cage_shaping_ring_thin_ssPhys[i]);
    }

    //----------------------------------------//
    // Extraction_grid                        //
    //----------------------------------------//
    G4double Extraction_gridOuterL2 = 75.5 / 2.0 * mm;  // estimated, it's surely smaller than 75.896 mm
    G4double Extraction_gridH2 = 0.05 / 2.0 * mm;
    G4double Extraction_gridX = 0.0 * mm;
    G4double Extraction_gridY = 0.0 * mm;
    G4double Extraction_gridZ = 25.025 * mm;
    G4Box* Extraction_gridSolid = new G4Box("Extraction_gridSolid", Extraction_gridOuterL2, Extraction_gridOuterL2, Extraction_gridH2);
    G4LogicalVolume* Extraction_gridLogic = new G4LogicalVolume(Extraction_gridSolid, NSLiquidArgon,
                                                                "Extraction_gridLogic");  // TODO: it's a ssteel mesh, have to create the
                                                                                          // material
    G4ThreeVector Extraction_gridPos(Extraction_gridX, Extraction_gridY, Extraction_gridZ);
    G4PVPlacement* Extraction_gridPhys = new G4PVPlacement(0, Extraction_gridPos, "Extraction_gridPhys", Extraction_gridLogic, fPhysicTPC, false, 0, myCheckOverlap);
    Debug(Extraction_gridPhys);

    //----------------------------------------//
    // Anode_plate                            //
    //----------------------------------------//
    G4double Anode_plateOuterL2 = 45.0 * mm;
    G4double Anode_plateOuterR = 12.7 * mm;
    G4double Anode_plateH = 4.5 * mm - 2.0 * coatingITOH - coatingTPBH;  // two layers of ITO and one (below) of TPB
    G4double Anode_plateH2 = Anode_plateH / 2.0;
    G4double Anode_plateX = 0.0 * mm;
    G4double Anode_plateY = 0.0 * mm;
    G4double Anode_plateZ = 37.3 * mm + coatingTPBH2;  // mean z is higher by half the TPB coating
    G4ThreeVector Anode_platePos(Anode_plateX, Anode_plateY, Anode_plateZ);
    G4SubtractionSolid* Anode_plateSolid = G4QuadPlateWithFillet(Anode_plateOuterL2, Anode_plateH2, Anode_plateOuterR);
    G4LogicalVolume* Anode_plateLogic = new G4LogicalVolume(Anode_plateSolid, Acrylic, "Anode_plateLogic");
    G4PVPlacement* Anode_platePhys = new G4PVPlacement(0, Anode_platePos, "Anode_platePhys", Anode_plateLogic, fPhysicTPC, false, 0, myCheckOverlap);
    Debug(Anode_platePhys);
    std::cout << "Anode H " << Anode_plateH / mm << " mm  IT0 " << coatingITOH / nm << " nm  TPB " << coatingTPBH / um << " um" << std::endl;
    //----------------------------------------//
    // Anode_plate ITO                        //
    //----------------------------------------//
    G4double Anode_plateITOX = Anode_plateX;
    G4double Anode_plateITOY = Anode_plateY;
    G4double Anode_plateITOZ[2] = {Anode_plateZ + Anode_plateH2 + coatingITOH2, Anode_plateZ - Anode_plateH2 - coatingITOH2};
    G4SubtractionSolid* Anode_plateITOSolid = G4QuadPlateWithFillet(Anode_plateOuterL2, coatingITOH2, Anode_plateOuterR);
    G4LogicalVolume* Anode_plateITOLogic = new G4LogicalVolume(Anode_plateITOSolid, ITO, "Anode_plateITOLogic");
    G4PVPlacement* Anode_plateITOPhys[2];
    for (int i = 0; i < 2; ++i) {  // I should learn to use repeated volumes
      char name[256];
      snprintf(name,90, "Anode_plateITOPhys%d", i);
      G4ThreeVector pos(Anode_plateITOX, Anode_plateITOY, Anode_plateITOZ[i]);
      Anode_plateITOPhys[i] = new G4PVPlacement(0, pos, name, Anode_plateITOLogic, fPhysicTPC, false, 0, myCheckOverlap);
      Debug(Anode_plateITOPhys[i]);
    }
    //----------------------------------------//
    // Anode_plate TPB                        //
    //----------------------------------------//
    G4double Anode_plateTPBX = Anode_plateX;
    G4double Anode_plateTPBY = Anode_plateY;
    G4double Anode_plateTPBZ = Anode_plateZ - Anode_plateH2 - coatingITOH - coatingTPBH2;  // below the plate
    G4ThreeVector Anode_plateTPBpos(Anode_plateTPBX, Anode_plateTPBY, Anode_plateTPBZ);
    G4SubtractionSolid* Anode_plateTPBSolid = G4QuadPlateWithFillet(Anode_plateOuterL2, coatingTPBH2, Anode_plateOuterR);
    G4LogicalVolume* Anode_plateTPBLogic = new G4LogicalVolume(Anode_plateTPBSolid, TPB, "Anode_plateTPBLogic");
    G4PVPlacement* Anode_plateTPBPhys = new G4PVPlacement(0, Anode_plateTPBpos, "Anode_plateTPBPhys", Anode_plateTPBLogic, fPhysicTPC, false, 0, myCheckOverlap);
    Debug(Anode_plateTPBPhys);

    //----------------------------------------//
    // Diving_bell_rim                        // simplified, a probable TODO
    //----------------------------------------//
    G4double Diving_bell_rimOuterL2 = 45.0 * mm;
    G4double Diving_bell_rimInnerL2 = 36.0 * mm;
    G4double Diving_bell_rimOuterR = 12.7 * mm;
    G4double Diving_bell_rimInnerR = 3.5 * mm;
    G4double Diving_bell_rimH2 = 6.952 / 2.0 * mm;
    G4double Diving_bell_rimX = 0.0 * mm;
    G4double Diving_bell_rimY = 0.0 * mm;
    G4double Diving_bell_rimZ = 30.05 * mm;
    G4ThreeVector Diving_bell_rimPos(Diving_bell_rimX, Diving_bell_rimY, Diving_bell_rimZ);
    G4SubtractionSolid* Diving_bell_rimSolid = G4QuadPlateWithFillet(Diving_bell_rimOuterL2, Diving_bell_rimH2, Diving_bell_rimOuterR, Diving_bell_rimInnerL2, Diving_bell_rimInnerR);
    G4LogicalVolume* Diving_bell_rimLogic = new G4LogicalVolume(Diving_bell_rimSolid, Teflon, "Diving_bell_rimLogic");
    G4PVPlacement* Diving_bell_rimPhys = new G4PVPlacement(0, Diving_bell_rimPos, "Diving_bell_rimPhys", Diving_bell_rimLogic, fPhysicTPC, false, 0, myCheckOverlap);
    Debug(Diving_bell_rimPhys);

    //----------------------------------------//
    // Anode_shaping_ring                     //  // I now notice that this ring
    // is identical to the two Field_cage_shaping_ring_thin_ss
    //----------------------------------------//
    G4double Anode_shaping_ringOuterL2 = 46.0 * mm;
    G4double Anode_shaping_ringD2 = 10.0 / 2.0 * mm;
    G4double Anode_shaping_ringOuterR = 13.5 * mm;
    // G4double Anode_shaping_ringInnerR = 3.5 * mm;
    G4double Anode_shaping_ringH2 = 1.524 / 2.0 * mm;
    G4double Anode_shaping_ringBarX = Anode_shaping_ringOuterL2 - Anode_shaping_ringD2;
    G4double Anode_shaping_ringX = Anode_shaping_ringBarX;
    G4double Anode_shaping_ringY = 0.0 * mm;
    G4double Anode_shaping_ringZ = 34.288 * mm;
    G4ThreeVector Anode_shaping_ringPos(Anode_shaping_ringX, Anode_shaping_ringY, Anode_shaping_ringZ);
    G4UnionSolid* Anode_shaping_ringSolid = G4QuadRingWithFillet(Anode_shaping_ringBarX, Anode_shaping_ringH2, Anode_shaping_ringD2, Anode_shaping_ringOuterR);
    G4LogicalVolume* Anode_shaping_ringLogic = new G4LogicalVolume(Anode_shaping_ringSolid, StainlessSteel, "Anode_shaping_ringLogic");
    Anode_shaping_ringLogic->SetVisAttributes(StainlessSteelVis);
    G4PVPlacement* Anode_shaping_ringPhys = new G4PVPlacement(0, Anode_shaping_ringPos, "Anode_shaping_ringPhys", Anode_shaping_ringLogic, fPhysicTPC, false, 0, myCheckOverlap);
    Debug(Anode_shaping_ringPhys);

    //----------------------------------------------------------------------------------//
    // Reflection_gas_side_back: // model the three volumes as staggered boxes,
    // the last one containing inactive LAr // and the gas pocket. //
    //----------------------------------------------------------------------------------//
    //----------------------------------------//
    // Reflection_gas_side_back               //
    //----------------------------------------//
    G4double Reflection_gas_side_backOuterL2 = 28.05 * mm;
    G4double Reflection_gas_side_backH2 = 10.0 / 2.0 * mm;
    G4double Reflection_gas_side_backX = 0.0 * mm;
    G4double Reflection_gas_side_backY = 0.0 * mm;
    G4double Reflection_gas_side_backZ = 30.05 * mm;
    G4ThreeVector Reflection_gas_side_backPos(Reflection_gas_side_backX, Reflection_gas_side_backY, Reflection_gas_side_backZ);
    G4Box* Reflection_gas_side_backSolid = new G4Box("Reflection_gas_side_backSolid", Reflection_gas_side_backOuterL2, Reflection_gas_side_backOuterL2, Reflection_gas_side_backH2);
    G4LogicalVolume* Reflection_gas_side_backLogic = new G4LogicalVolume(Reflection_gas_side_backSolid, Acrylic, "Reflection_gas_side_backLogic");
    G4PVPlacement* Reflection_gas_side_backPhys = new G4PVPlacement(0, Reflection_gas_side_backPos, "Reflection_gas_side_backPhys", Reflection_gas_side_backLogic, fPhysicTPC, false, 0, myCheckOverlap);
    Debug(Reflection_gas_side_backPhys);
    //----------------------------------------//
    // Reflection_gas_side_ESR                //
    //----------------------------------------//
    G4double Reflection_gas_side_ESROuterL2 = 26.55 * mm;
    G4double Reflection_gas_side_ESRH2 = 10.0 / 2.0 * mm;
    G4Box* Reflection_gas_side_ESRSolid = new G4Box("Reflection_gas_side_ESRSolid", Reflection_gas_side_ESROuterL2, Reflection_gas_side_ESROuterL2, Reflection_gas_side_ESRH2);
    G4LogicalVolume* Reflection_gas_side_ESRLogic = new G4LogicalVolume(Reflection_gas_side_ESRSolid, ThreeMFoil, "Reflection_gas_side_ESRLogic");
    G4PVPlacement* Reflection_gas_side_ESRPhys = new G4PVPlacement(0, myZeros, "Reflection_gas_side_ESRPhys", Reflection_gas_side_ESRLogic, Reflection_gas_side_backPhys, false, 0, myCheckOverlap);
    Debug(Reflection_gas_side_ESRPhys);
    //----------------------------------------//
    // Reflection_gas_side_front              //
    //----------------------------------------//
    G4double Reflection_gas_side_frontOuterL2 = 26.5 * mm;
    G4double Reflection_gas_side_frontInnerL2 = 25.0 * mm;
    G4double Reflection_gas_side_frontH2 = 10.0 / 2.0 * mm;
    G4Box* Reflection_gas_side_frontSolid = new G4Box("Reflection_gas_side_frontSolid", Reflection_gas_side_frontOuterL2, Reflection_gas_side_frontOuterL2, Reflection_gas_side_frontH2);
    G4LogicalVolume* Reflection_gas_side_frontLogic = new G4LogicalVolume(Reflection_gas_side_frontSolid, Acrylic, "Reflection_gas_side_frontLogic");
    G4PVPlacement* Reflection_gas_side_frontPhys = new G4PVPlacement(0, myZeros, "Reflection_gas_side_frontPhys", Reflection_gas_side_frontLogic, Reflection_gas_side_ESRPhys, false, 0, myCheckOverlap);
    Debug(Reflection_gas_side_frontPhys);
    //----------------------------------------//
    // Reflection_gas_side_front TPB        //
    //----------------------------------------//
    G4double Reflection_gas_side_frontTPBOuterL2 = Reflection_gas_side_frontInnerL2;
    G4double Reflection_gas_side_frontTPBInnerL2 = Reflection_gas_side_frontTPBOuterL2 - coatingTPBH;  // needed for placement of inactiveLAr and GAr
    G4double Reflection_gas_side_frontTPBH2 = Reflection_gas_side_frontH2;
    G4Box* Reflection_gas_side_frontTPBSolid = new G4Box("Reflection_gas_side_frontTPBSolid", Reflection_gas_side_frontTPBOuterL2, Reflection_gas_side_frontTPBOuterL2, Reflection_gas_side_frontTPBH2);
    G4LogicalVolume* Reflection_gas_side_frontTPBLogic = new G4LogicalVolume(Reflection_gas_side_frontTPBSolid, TPB, "Reflection_gas_side_frontTPBLogic");
    G4PVPlacement* Reflection_gas_side_frontTPBPhys = new G4PVPlacement(0, myZeros, "Reflection_gas_side_frontTPBPhys", Reflection_gas_side_frontTPBLogic, Reflection_gas_side_frontPhys, false, 0, myCheckOverlap);
    Debug(Reflection_gas_side_frontTPBPhys);

    //----------------------------------------//
    // Bottom_support                         // simplified, a probable TODO
    //----------------------------------------//
    G4double Bottom_supportOuterL2 = 45.0 * mm;
    G4double Bottom_supportInnerL2 = 31.0 * mm;
    G4double Bottom_supportOuterR = 12.7 * mm;
    G4double Bottom_supportInnerR = 0.0 * mm;
    G4double Bottom_supportH2 = 7.493 / 2.0 * mm;
    G4double Bottom_supportX = 0.0 * mm;
    G4double Bottom_supportY = 0.0 * mm;
    G4double Bottom_supportZ = -35.2785 * mm;
    G4ThreeVector Bottom_supportPos(Bottom_supportX, Bottom_supportY, Bottom_supportZ);
    G4SubtractionSolid* Bottom_supportSolid = G4QuadPlateWithFillet(Bottom_supportOuterL2, Bottom_supportH2, Bottom_supportOuterR, Bottom_supportInnerL2, Bottom_supportInnerR);
    G4LogicalVolume* Bottom_supportLogic = new G4LogicalVolume(Bottom_supportSolid, Teflon, "Bottom_supportLogic");
    G4PVPlacement* Bottom_supportPhys = new G4PVPlacement(0, Bottom_supportPos, "Bottom_supportPhys", Bottom_supportLogic, fPhysicTPC, false, 0, myCheckOverlap);
    Debug(Bottom_supportPhys);

    //----------------------------------------//
    // Field_cage_shaping_ring_thick_copper   //
    //----------------------------------------//
    G4double Field_cage_shaping_ring_thick_copperOuterL2 = 38.0 * mm;
    G4double Field_cage_shaping_ring_thick_copperD2 = 5.5 / 2.0 * mm;
    G4double Field_cage_shaping_ring_thick_copperOuterR = 8.0 * mm;
    // G4double Field_cage_shaping_ring_thick_copperInnerR = 2.5 * mm;
    G4double Field_cage_shaping_ring_thick_copperH2 = 2.032 / 2.0 * mm;
    G4double Field_cage_shaping_ring_thick_copperZ[2] = {-30.516 * mm, -23.984 * mm};
    G4double Field_cage_shaping_ring_thick_copperBarX = Field_cage_shaping_ring_thick_copperOuterL2 - Field_cage_shaping_ring_thick_copperD2;
    G4UnionSolid* Field_cage_shaping_ring_thick_copperSolid = G4QuadRingWithFillet(Field_cage_shaping_ring_thick_copperBarX, Field_cage_shaping_ring_thick_copperH2, Field_cage_shaping_ring_thick_copperD2, Field_cage_shaping_ring_thick_copperOuterR);
    G4LogicalVolume* Field_cage_shaping_ring_thick_copperLogic = new G4LogicalVolume(Field_cage_shaping_ring_thick_copperSolid, MetalCopper, "Field_cage_shaping_ring_thick_copperLogic");
    Field_cage_shaping_ring_thick_copperLogic->SetVisAttributes(MetalCopperVis);
    G4PVPlacement* Field_cage_shaping_ring_thick_copperPhys[2];
    for (int i = 0; i < 2; ++i) {  // I should learn to use repeated volumes
      char name[256];
      snprintf(name, 90,"Field_cage_shaping_ring_thick_copperPhys%d", i);
      G4ThreeVector pos(Field_cage_shaping_ring_thick_copperBarX, 0, Field_cage_shaping_ring_thick_copperZ[i]);
      Field_cage_shaping_ring_thick_copperPhys[i] = new G4PVPlacement(0, pos, name, Field_cage_shaping_ring_thick_copperLogic, fPhysicTPC, false, 0, myCheckOverlap);
      Debug(Field_cage_shaping_ring_thick_copperPhys[i]);
    }

    //----------------------------------------//
    // Cathode_plate                          //
    //----------------------------------------//
    G4double Cathode_plateOuterL2 = 38.0 * mm;
    G4double Cathode_plateOuterR = 8.0 * mm;
    G4double Cathode_plateH = 4.5 * mm - 2.0 * coatingITOH - coatingTPBH;  // two layers of ITO and one (above) of TPB
    G4double Cathode_plateH2 = Cathode_plateH / 2.0;
    G4double Cathode_plateX = 0.0 * mm;
    G4double Cathode_plateY = 0.0 * mm;
    G4double Cathode_plateZ = -27.25 * mm - coatingTPBH2;  // mean z is lower by half the TPB coating
    G4ThreeVector Cathode_platePos(Cathode_plateX, Cathode_plateY, Cathode_plateZ);
    G4SubtractionSolid* Cathode_plateSolid = G4QuadPlateWithFillet(Cathode_plateOuterL2, Cathode_plateH2, Cathode_plateOuterR);
    G4LogicalVolume* Cathode_plateLogic = new G4LogicalVolume(Cathode_plateSolid, Acrylic, "Cathode_plateLogic");
    G4PVPlacement* Cathode_platePhys = new G4PVPlacement(0, Cathode_platePos, "Cathode_platePhys", Cathode_plateLogic, fPhysicTPC, false, 0, myCheckOverlap);
    Debug(Cathode_platePhys);
    std::cout << "Cathode H " << Cathode_plateH / mm << " mm  IT0 " << coatingITOH / nm << " nm  TPB " << coatingTPBH / um << " um" << std::endl;
    //----------------------------------------//
    // Cathode_plate ITO                        //
    //----------------------------------------//
    G4double Cathode_plateITOX = Cathode_plateX;
    G4double Cathode_plateITOY = Cathode_plateY;
    G4double Cathode_plateITOZ[2] = {Cathode_plateZ + Cathode_plateH2 + coatingITOH2, Cathode_plateZ - Cathode_plateH2 - coatingITOH2};
    G4SubtractionSolid* Cathode_plateITOSolid = G4QuadPlateWithFillet(Cathode_plateOuterL2, coatingITOH2, Cathode_plateOuterR);
    G4LogicalVolume* Cathode_plateITOLogic = new G4LogicalVolume(Cathode_plateITOSolid, ITO, "Cathode_plateITOLogic");
    G4PVPlacement* Cathode_plateITOPhys[2];
    for (int i = 0; i < 2; ++i) {  // I should learn to use repeated volumes
      char name[256];
      snprintf(name,30, "Cathode_plateITOPhys%d", i);
      G4ThreeVector pos(Cathode_plateITOX, Cathode_plateITOY, Cathode_plateITOZ[i]);
      Cathode_plateITOPhys[i] = new G4PVPlacement(0, pos, name, Cathode_plateITOLogic, fPhysicTPC, false, 0, myCheckOverlap);
      Debug(Cathode_plateITOPhys[i]);
    }
    //----------------------------------------//
    // Cathode_plate TPB                        //
    //----------------------------------------//
    G4double Cathode_plateTPBX = Cathode_plateX;
    G4double Cathode_plateTPBY = Cathode_plateY;
    G4double Cathode_plateTPBZ = Cathode_plateZ + Cathode_plateH2 + coatingITOH + coatingTPBH2;  // above the plate
    G4ThreeVector Cathode_plateTPBpos(Cathode_plateTPBX, Cathode_plateTPBY, Cathode_plateTPBZ);
    G4SubtractionSolid* Cathode_plateTPBSolid = G4QuadPlateWithFillet(Cathode_plateOuterL2, coatingTPBH2, Cathode_plateOuterR);
    G4LogicalVolume* Cathode_plateTPBLogic = new G4LogicalVolume(Cathode_plateTPBSolid, TPB, "Cathode_plateTPBLogic");
    G4PVPlacement* Cathode_plateTPBPhys = new G4PVPlacement(0, Cathode_plateTPBpos, "Cathode_plateTPBPhys", Cathode_plateTPBLogic, fPhysicTPC, false, 0, myCheckOverlap);
    Debug(Cathode_plateTPBPhys);

    //----------------------------------------//
    // readout_bd                             //  180322: Bianca reported that
    // both readout boards are different! TODO
    //----------------------------------------//
    G4double readout_bdOuterL2 = 38.5 * mm;
    G4double readout_bdH2 = 2.0 / 2.0 * mm;
    G4double readout_bdX = 0.0 * mm;
    G4double readout_bdY = 0.0 * mm;
    G4double readout_bdZ[2] = {-40.025 * mm, 50.075 * mm};
    G4PVPlacement* readout_bdPhys[2];
    G4Box* readout_bdSolid = new G4Box("readout_bdSolid", readout_bdOuterL2, readout_bdOuterL2, readout_bdH2);
    G4LogicalVolume* readout_bdLogic = new G4LogicalVolume(readout_bdSolid, Teflon, "readout_bdLogic");  // TODO, wrong material
    for (int i = 0; i < 2; ++i) {                                                                        // I should learn to use repeated volumes
      char name[256];
      snprintf(name, 90,"readout_bdPhys%d", i);
      G4ThreeVector readout_bdPos(readout_bdX, readout_bdY, readout_bdZ[i]);
      readout_bdPhys[i] = new G4PVPlacement(0, readout_bdPos, name, readout_bdLogic, fPhysicTPC, false, 0, myCheckOverlap);
      Debug(readout_bdPhys[i]);
    }

    //----------------------------------------//
    // SiPM_Tile_Sub                          //
    //----------------------------------------//
    G4double SiPM_Tile_SubOuterL2 = 30.0 * mm;
    G4double SiPM_Tile_SubH2 = 1.64 / 2.0 * mm;
    G4double SiPM_Tile_SubX = 0.0 * mm;
    G4double SiPM_Tile_SubY = 0.0 * mm;
    G4double SiPM_Tile_SubZ[2] = {-31.845 * mm, 41.895 * mm};
    G4PVPlacement* SiPM_Tile_SubPhys[2];
    G4Box* SiPM_Tile_SubSolid = new G4Box("SiPM_Tile_SubSolid", SiPM_Tile_SubOuterL2, SiPM_Tile_SubOuterL2, SiPM_Tile_SubH2);
    G4LogicalVolume* SiPM_Tile_SubLogic = new G4LogicalVolume(SiPM_Tile_SubSolid, Teflon, "SiPM_Tile_SubLogic");
    for (int i = 0; i < 2; ++i) {  // I should learn to use repeated volumes
      char name[256];
      snprintf(name, 90,"SiPM_Tile_SubPhys%d", i);
      G4ThreeVector SiPM_Tile_SubPos(SiPM_Tile_SubX, SiPM_Tile_SubY, SiPM_Tile_SubZ[i]);
      SiPM_Tile_SubPhys[i] = new G4PVPlacement(0, SiPM_Tile_SubPos, name, SiPM_Tile_SubLogic, fPhysicTPC, false, 0, myCheckOverlap);
      Debug(SiPM_Tile_SubPhys[i]);
    }

    //----------------------------------------//
    // SiPM                                   //
    //----------------------------------------//
    // I could also make two monolithic SiPM plates, but two arrays of 32 each
    // look nicer
    G4double SiPMOuterLX2 = 5.95 / 2.0 * mm;
    G4double SiPMOuterLY2 = 11.95 / 2.0 * mm;
    G4double SiPMH2 = 0.5 / 2.0 * mm;
    G4PVPlacement* SiPMPhys[64];
    G4Box* SiPMSolid = new G4Box("SiPMSolid", SiPMOuterLX2, SiPMOuterLY2, SiPMH2);
    G4LogicalVolume* SiPMLogic = new G4LogicalVolume(SiPMSolid, MetalSilicon, "SiPMLogic");
    G4int SiPMcounter = 0;
    for (G4double SiPMZ = -30.675; SiPMZ < 80.0; SiPMZ += 30.675 + 40.725) {
      for (G4double SiPMY = -18.0; SiPMY < 19.0; SiPMY += 12.0) {
        for (G4double SiPMX = -21.0; SiPMX < 22.0; SiPMX += 6.0) {
          char name[256];
          snprintf(name,30, "SiPMPhys%d", SiPMcounter);
          G4ThreeVector SiPMPos(SiPMX, SiPMY, SiPMZ);
          SiPMPhys[SiPMcounter] = new G4PVPlacement(0, SiPMPos, name, SiPMLogic, fPhysicTPC, false, 0, myCheckOverlap);
          Debug(SiPMPhys[SiPMcounter]);
          ++SiPMcounter;
        }
      }
    }

    //----------------------------------------//
    // Top_SiPM_holder                        // simplified, a probable TODO
    //----------------------------------------//
    G4double Top_SiPM_holderOuterL2 = 45.0 * mm;
    G4double Top_SiPM_holderInnerL2 = 31.0 * mm;
    G4double Top_SiPM_holderOuterR = 12.7 * mm;
    G4double Top_SiPM_holderInnerR = 0.0 * mm;
    G4double Top_SiPM_holderH2 = 9.525 / 2.0 * mm;
    G4double Top_SiPM_holderX = 0.0 * mm;
    G4double Top_SiPM_holderY = 0.0 * mm;
    G4double Top_SiPM_holderZ = 44.3125 * mm;
    G4ThreeVector Top_SiPM_holderPos(Top_SiPM_holderX, Top_SiPM_holderY, Top_SiPM_holderZ);
    G4SubtractionSolid* Top_SiPM_holderSolid = G4QuadPlateWithFillet(Top_SiPM_holderOuterL2, Top_SiPM_holderH2, Top_SiPM_holderOuterR, Top_SiPM_holderInnerL2, Top_SiPM_holderInnerR);
    G4LogicalVolume* Top_SiPM_holderLogic = new G4LogicalVolume(Top_SiPM_holderSolid, Teflon, "Top_SiPM_holderLogic");
    G4PVPlacement* Top_SiPM_holderPhys = new G4PVPlacement(0, Top_SiPM_holderPos, "Top_SiPM_holderPhys", Top_SiPM_holderLogic, fPhysicTPC, false, 0, myCheckOverlap);
    Debug(Top_SiPM_holderPhys);

    //-----------------------------//
    // active LAr                  //
    //-----------------------------//
    // the active LAr is completely and entirely within
    // Reflection_drift_side_front
    G4double TPCactiveLArL2 = Reflection_drift_side_frontTPBInnerL2;
    G4double TPCactiveLArH2 = Reflection_drift_side_frontTPBH2;
    G4Box* TPCactiveLArSolid = new G4Box("TPCactiveLArSolid", TPCactiveLArL2, TPCactiveLArL2, TPCactiveLArH2);
    // The name of the logical volume has to be "LAr_Logic".
    // It will be used as a G4Region, and DSPhysicsList scans the G4RegionStore
    // for "LAr_Logic"!
    G4LogicalVolume* TPCactiveLArLogic = new G4LogicalVolume(TPCactiveLArSolid, LiquidArgon, "LAr_Logic");
    TPCactiveLArLogic->SetVisAttributes(LiquidArgonVis);
    // DSLightX requires that this volume is called "ActiveLAr", otherwise S2
    // won't be generated!
    fActiveLArTPC = new G4PVPlacement(0, myZeros, "ActiveLAr", TPCactiveLArLogic, Reflection_drift_side_frontTPBPhys, false, 0, myCheckOverlap);
    Debug(fActiveLArTPC);

    //
    // the space within Reflection_gas_side_front is completely filled with
    // inactive LAr or GAr (from above ~2 mm above the grid)
    //

    //-----------------------------//
    // inactive LAr                //
    //-----------------------------//
    G4double TPCinactiveLArL2 = Reflection_gas_side_frontTPBInnerL2;
    G4double TPCinactiveLArH2 = Reflection_gas_side_frontTPBH2;
    G4Box* TPCinactiveLArSolid = new G4Box("TPCinactiveLArSolid", TPCinactiveLArL2, TPCinactiveLArL2, TPCinactiveLArH2);
    G4LogicalVolume* TPCinactiveLArLogic = new G4LogicalVolume(TPCinactiveLArSolid, NSLiquidArgon, "TPCinactiveLArLogic");
    TPCinactiveLArLogic->SetVisAttributes(NSLiquidArgonVis);
    fPhysicInactiveLAr = new G4PVPlacement(0, myZeros, "TPCinactiveLArPhys", TPCinactiveLArLogic, Reflection_gas_side_frontTPBPhys, false, 0, myCheckOverlap);
    Debug(fPhysicInactiveLAr);

    //-----------------------//
    //      Gas Pocket       //
    //-----------------------//
    G4double TPCGArZshift = 2.0 * mm;  // the gas pocket starts about 2 mm above the extraction grid
    G4double TPCGArL2 = Reflection_gas_side_frontTPBInnerL2;
    G4double TPCGArH2 = TPCinactiveLArH2 - TPCGArZshift / 2.0;
    G4double TPCGArX = 0.0 * mm;
    G4double TPCGArY = 0.0 * mm;
    G4double TPCGArZ = TPCGArZshift / 2.0;
    G4ThreeVector TPCGArPos(TPCGArX, TPCGArY, TPCGArZ);
    G4Box* GasPocketSolid = new G4Box("GasPocketSolid", TPCGArL2, TPCGArL2, TPCGArH2);
    G4LogicalVolume* GasPocketLogic = new G4LogicalVolume(GasPocketSolid, GaseousArgon, "GasPocketLogic");
    GasPocketLogic->SetVisAttributes(GaseousArgonVis);
    G4PVPlacement* GasPocketPhys = new G4PVPlacement(0, TPCGArPos, "GasPocketPhys", GasPocketLogic, fPhysicInactiveLAr, false, 0, myCheckOverlap);
    Debug(GasPocketPhys);

    // Set the z coordinate of the LAr - GAr interface, necessary for S2
    // generation in DSLightX
    DSStorage::Get()->SetLArGArBoundaryPosZ(TPCactiveLArH2 + TPCGArZshift);
    DSLog(routine) << "TPC " << tpc << " gas pocket position: " << DSStorage::Get()->GetLArGArBoundaryPosZ() / mm << " TPCactiveLArH2 " << TPCactiveLArH2 / mm << " TPCGArZshift " << TPCGArZshift / mm << " mm" << endlog;

    //----------------------------------------------------------------------------
    // Sets the region (honestly, I don't understand this recursive stuff).
    // But that's how it is explained in the G4 v10.2 manual section 4.1.3.1.
    //----------------------------------------------------------------------------
    G4Region* LArRegion = new G4Region("LAr_Logic");
    TPCactiveLArLogic->SetRegion(LArRegion);
    LArRegion->AddRootLogicalVolume(TPCactiveLArLogic);

    DefineSurfaces();

    // make SiPM as pe storing material
    DSStorage::Get()->SetPMTMaterialIndex(SiPMLogic->GetMaterial()->GetIndex());
  }
}

DSDetectorReD::~DSDetectorReD() {
  ;  // delete fMessenger;
}

void DSDetectorReD::SetToVacuum(G4Material*& _m) {
  G4String _s = _m->GetName();
  _m = DSMaterial::Get()->GetVacuum();
  _m->SetName(_s);
  DSLog(routine) << _s << " set to " << _m->GetName() << endlog;
  DSLog(routine) << _m << endlog;
}

void DSDetectorReD::DefineSurfaces() {
  ////////////////////////////////////////
  // Copied from DS DetectorDS50 - 24th april 2015
  ////////////////////////////////////////

  ////////////////////////////////////////
  // LAr -> Grid
  //  Note: in this model, all optical action takes place
  //  on entering the grid.
  ////////////////////////////////////////

  G4OpticalSurface* fOpGridLArSurface = new G4OpticalSurface("OpGridLArSurface");
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
  G4MaterialPropertiesTable* fGridLArSurfProp = new G4MaterialPropertiesTable();
  if (DSParameters::Get()->GetWithNewGridModel())
    // the grid model is described in DSStorage.cc, the surface is treated in
    // G4OpBoundaryProcess.cc
    fGridLArSurfProp->AddConstProperty("DOGRID", 1);
  // Now use the following in old and new models.  By G4 convention,
  // "reflectivity" is actually 1-absorption.
  fGridLArSurfProp->AddProperty("REFLECTIVITY", GridLArENE, GridLArREF, 4);
  //  else {
  //    fGridLArSurfProp->AddProperty("REFLECTIVITY", GridLArENE, GridLArREF,
  //    4);
  //    }
  fOpGridLArSurface->SetMaterialPropertiesTable(fGridLArSurfProp);

  ////////////////////////////////////////
  // Grid -> LAr (keeping backward labeling convention)
  //  Note: in this model, all optical action takes place
  //  on entering the grid.  Exit action is to just continue
  //  in straight line.
  ////////////////////////////////////////

  if (DSParameters::Get()->GetWithNewGridModel()) {
    G4OpticalSurface* fOpLArGridSurface = new G4OpticalSurface("OpLArGridSurface");
    new G4LogicalBorderSurface("LArGridSurface", fPhysicGrid, fActiveLArTPC, fOpLArGridSurface);

    fOpLArGridSurface->SetType(dielectric_dielectric);
    //  fOpLArGridSurface->SetModel( glisur );
    //  fOpLArGridSurface->SetFinish( polished );
    G4MaterialPropertiesTable* fLArGridSurfProp = new G4MaterialPropertiesTable();
    // the grid model is described in DSStorage.cc, the surface is treated in
    // G4OpBoundaryProcess.cc
    fLArGridSurfProp->AddConstProperty("DOGRIDEXIT", 1);
    fOpLArGridSurface->SetMaterialPropertiesTable(fLArGridSurfProp);
  }

  ////////////////////////////////////////
  // TPB <--> GAr and TPB <--> LAr
  // This surface will carry all the diffuse properties of the TPB
  // for both GAr and LAr.
  // Make this bi-directional
  ////////////////////////////////////////
  G4OpticalSurface* fOpTPBGArSurface = new G4OpticalSurface("OpTPBGArSurface");
  // Note: the following two work even when the "gas pocket" is LAr in no-pocket
  // runs.
  new G4LogicalBorderSurface("GArTPBSurface", fPhysicGasPocket, fPhysicTPBGAr, fOpTPBGArSurface);
  new G4LogicalBorderSurface("TPBGArSurface", fPhysicTPBGAr, fPhysicGasPocket, fOpTPBGArSurface);
  new G4LogicalBorderSurface("LArTPBSurface", fActiveLArTPC, fPhysicTPBLAr, fOpTPBGArSurface);
  new G4LogicalBorderSurface("TPBLArSurface", fPhysicTPBLAr, fActiveLArTPC, fOpTPBGArSurface);
  // new G4LogicalBorderSurface("TPBLArLayerSurface", fPhysicTPB,
  // fPhysicLArLayer, fOpTPBGArSurface );
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

  G4MaterialPropertiesTable* fTPBGArSurfProp = new G4MaterialPropertiesTable();

  fTPBGArSurfProp->AddProperty("SPECULARLOBECONSTANT", pp, specularlobe, NUM);
  fTPBGArSurfProp->AddProperty("SPECULARSPIKECONSTANT", pp, specularspike, NUM);
  fTPBGArSurfProp->AddProperty("BACKSCATTERCONSTANT", pp, backscatter, NUM);
  fTPBGArSurfProp->AddProperty("REFLECTIVITY", pp, reflectivity, NUM);
  fTPBGArSurfProp->AddProperty("TRANSMITTANCE", pp, transmitivity, NUM);

  fTPBGArSurfProp->AddConstProperty("DOArTPB", 1);

  fOpTPBGArSurface->SetMaterialPropertiesTable(fTPBGArSurfProp);

  ////////////////////////////////////////
  // ITO /////
  ////////////////////////////////////////

  G4MaterialPropertiesTable* fITOSurfProp = new G4MaterialPropertiesTable();
  fITOSurfProp->AddConstProperty("DOITO", 1);

  ////////////////////////////////////////
  // AnodeWindow <--> TPB and CathodeWindow <--> TPB
  // In the current model, the diffuse nature of the TPB is handled
  // entirely at the TPB-GAr/LAr surface, not here.
  // Make this bi-directional.
  ////////////////////////////////////////
  G4OpticalSurface* fOpWindowTPBSurface = new G4OpticalSurface("OpWindowTPBSurface");
  new G4LogicalBorderSurface("BottomWindowTPBSurface", fPhysicCathodeWindow, fPhysicTPBLAr, fOpWindowTPBSurface);
  new G4LogicalBorderSurface("TPBBottomWindowSurface", fPhysicTPBLAr, fPhysicCathodeWindow, fOpWindowTPBSurface);
  new G4LogicalBorderSurface("TopWindowTPBSurface", fPhysicAnodeWindow, fPhysicTPBGAr, fOpWindowTPBSurface);
  new G4LogicalBorderSurface("TPBTopWindowSurface", fPhysicTPBGAr, fPhysicAnodeWindow, fOpWindowTPBSurface);
  fOpWindowTPBSurface->SetType(dielectric_dielectric);
  fOpWindowTPBSurface->SetModel(unified);
  //  fOpWindowTPBSurface->SetFinish( ground );
  fOpWindowTPBSurface->SetFinish(polished);
  // G4MaterialPropertiesTable *fWindowTPBSurfProp = new
  // G4MaterialPropertiesTable();
  fOpWindowTPBSurface->SetMaterialPropertiesTable(fITOSurfProp);

  ///////////////
  // TPB --> Teflon (Reflector)
  //  Should be no Teflon --> TPB as surface is defined as dielectric_metal.
  ////////////////////////////////////////
  G4double TeflonTPBENE[4] = {0.1 * eV, 8.0 * eV, 8.3 * eV, 20.0 * eV};
  G4double TREFUV = DSParameters::Get()->GetTeflonTPBUVRef();
  G4double TREFVIS = DSParameters::Get()->GetTeflonTPBVisRef();
  G4double TeflonTPBREF[4] = {TREFVIS, TREFVIS, TREFUV, TREFUV};

  G4OpticalSurface* fOpTPBTeflonSurface = new G4OpticalSurface("OpTBPTeflonSurface");
  new G4LogicalBorderSurface("TPBLArTeflonSurface", fPhysicTPBLAr, fPhysicTeflonFoil, fOpTPBTeflonSurface);
  new G4LogicalBorderSurface("TPBGArTeflonSurface", fPhysicTPBGAr, fPhysicTeflonFoil, fOpTPBTeflonSurface);
  fOpTPBTeflonSurface->SetType(dielectric_metal);
  fOpTPBTeflonSurface->SetModel(unified);

  // PDM: though I can't see how, the following settings are giving the desired
  // Lambertian reflection
  fOpTPBTeflonSurface->SetFinish(groundfrontpainted);
  // fOpTPBTeflonSurface->SetFinish(ground);
  fOpTPBTeflonSurface->SetSigmaAlpha(0.1);

  G4MaterialPropertiesTable* fTPBTeflonSurfProp = new G4MaterialPropertiesTable();
  fTPBTeflonSurfProp->AddProperty("REFLECTIVITY", TeflonTPBENE, TeflonTPBREF, 4);
  fOpTPBTeflonSurface->SetMaterialPropertiesTable(fTPBTeflonSurfProp);

  ////////////////////////////////////////
  // Teflon - LAr
  //  PDM: These refer to EXTERNAL surfaces of the TPC.  I don't think they will
  //  be executed and they haven't been debugged.  (Some are defined only for
  //  teflon --> LAr, which shouldn't happen.)
  ////////////////////////////////////////
  G4OpticalSurface* fOpLArTeflonSurface = new G4OpticalSurface("OpLArTeflonSurface");

  //  new G4LogicalBorderSurface("LArTeflonSurface", fPhysicTeflonFoil,
  //  fPhysicInactiveLAr, fOpLArTeflonSurface );
  new G4LogicalBorderSurface("LArTeflonSurface", fPhysicTeflonFoil, fPhysicTPC, fOpLArTeflonSurface);

  // new G4LogicalBorderSurface("LArPMTTopAssemblySurface",
  // fPhysicPMTAssemblyTop, fPhysicInactiveLAr, fOpLArTeflonSurface); new
  // G4LogicalBorderSurface("LArPMTBottomAssemblySurface",
  // fPhysicPMTAssemblyBottom, fPhysicInactiveLAr, fOpLArTeflonSurface);

  // MWK 170510: I don't understand why also this new object is saved in
  // fOpLArTeflonSurface. probably a mistake. Anyway, what does it mean that
  // vol1 and vol2 are reversed here?

  //  new G4LogicalBorderSurface("LArTeflonSupportSurface", fPhysicInactiveLAr,
  //  fPhysicTeflonFoil, fOpLArTeflonSurface);
  new G4LogicalBorderSurface("LArTeflonSupportSurface", fPhysicTPC, fPhysicTeflonFoil, fOpLArTeflonSurface);

  fOpLArTeflonSurface->SetType(dielectric_metal);
  fOpLArTeflonSurface->SetModel(glisur);
  fOpLArTeflonSurface->SetFinish(ground);
  G4MaterialPropertiesTable* fLArTeflonSurfProp = new G4MaterialPropertiesTable();
  G4double TeflonLArENE[4] = {0.1 * eV, 8.0 * eV, 8.3 * eV, 20.0 * eV};
  G4double TAREFUV = DSParameters::Get()->GetTeflonLArUVRef();
  G4double TAREFVIS = DSParameters::Get()->GetTeflonLArVisRef();
  G4double TeflonLArREF[4] = {TAREFVIS, TAREFVIS, TAREFUV, TAREFUV};
  fLArTeflonSurfProp->AddProperty("REFLECTIVITY", TeflonLArENE, TeflonLArREF, 4);
  fOpLArTeflonSurface->SetMaterialPropertiesTable(fLArTeflonSurfProp);

  return;
}

G4double DSDetectorReD::GetCm(G4double inches) {
  return inches * 2.54;
}

G4UnionSolid* DSDetectorReD::G4QuadRingWithFillet(G4double x, G4double h2, G4double d2, G4double R) {
  // The quadratic ring is composed of for straight segments, connected by arcs.
  // The parameters are:
  // x: distance of the center of each bar from the center of the ring
  // h2: half height of the ring
  // d2: half width of the ring
  // r: outer fillet radius of the ring

  G4double l2 = x + d2;
  G4double barL2 = l2 - R;

  G4Box* EdgeSolid = new G4Box("EdgeSolid", d2, barL2, h2);
  G4Tubs* ArcSolid = new G4Tubs("ArcSolid", R - 2.0 * d2, R, h2, 0, 0.5 * M_PI);
  G4ThreeVector shift(-R + d2, barL2, 0);
  G4UnionSolid* QuarterSolid = new G4UnionSolid("QuarterSolid", EdgeSolid, ArcSolid, 0, shift);
  G4RotationMatrix rotZ;
  rotZ.rotateZ(90.0 * deg);
  shift = G4ThreeVector(-x, x, 0);
  G4Transform3D transform(rotZ, shift);
  G4UnionSolid* HalfSolid = new G4UnionSolid("HalfSolid", QuarterSolid, QuarterSolid, transform);
  rotZ.rotateZ(90.0 * deg);
  shift = G4ThreeVector(-2.0 * x, 0, 0);
  transform = G4Transform3D(rotZ, shift);
  return new G4UnionSolid("Solid", HalfSolid, HalfSolid, transform);
}

G4SubtractionSolid* DSDetectorReD::G4QuadPlateWithFillet(G4double L2, G4double h2, G4double R, G4double l2, G4double _r) {
  // The quadratic plate with fillet and optional inner cut with optional inner
  // fillet. The parameters are: L2: outer half length h2: half height R:  outer
  // fillet radius l2: inner half length (if 0 means not cutted) r:  inner
  // fillet radius (if 0 means no radius)

  // FIXME: _r is an unused variable
  _r = 0;
  if (_r < 0) { cout << "Useless print, mandatory because _r is not used. This class should be fixed" << endl; };
  G4double x = L2 - R;

  G4Box* plate0Solid = new G4Box("plateSolid", L2, L2, h2);
  // the outer r should be sqrt(2)r, with 1.5r we are on the safe side
  // same for the height, make the arc a tad higher to avoid rounding problems
  G4Tubs* arcSolid = new G4Tubs("arcSolid", R, 1.5 * R, 1.01 * h2, 0, 0.5 * M_PI);
  G4RotationMatrix rotZ;
  G4ThreeVector shift(x, x, 0);
  G4Transform3D transform(rotZ, shift);
  G4SubtractionSolid* plate1Solid = new G4SubtractionSolid("plate1Solid", plate0Solid, arcSolid, transform);
  rotZ.rotateZ(90.0 * deg);
  shift = G4ThreeVector(-x, x, 0);
  transform = G4Transform3D(rotZ, shift);
  G4SubtractionSolid* plate2Solid = new G4SubtractionSolid("plate2Solid", plate1Solid, arcSolid, transform);
  rotZ.rotateZ(90.0 * deg);
  shift = G4ThreeVector(-x, -x, 0);
  transform = G4Transform3D(rotZ, shift);
  G4SubtractionSolid* plate3Solid = new G4SubtractionSolid("plate3Solid", plate2Solid, arcSolid, transform);
  rotZ.rotateZ(90.0 * deg);
  shift = G4ThreeVector(x, -x, 0);
  transform = G4Transform3D(rotZ, shift);
  G4SubtractionSolid* plate4Solid = new G4SubtractionSolid("plate4Solid", plate3Solid, arcSolid, transform);

  if (l2 == 0.0) return plate4Solid;

  // cut material

  // disabled, see below
  //  if ( r == 0.0 ) {  // just a plain plate
  G4Box* plateCutSolid = new G4Box("plateCutSolid", l2, l2, 1.01 * h2);
  return new G4SubtractionSolid("plate5Solid", plate4Solid, plateCutSolid);
  //  }

  // cut a plate with fillet (hopefully this won't result in a recursive chain
  // under certain conditions)
  // TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
  // TODO However, at run time this gives me: ERROR:
  // G4VSceneHandler::RequestPrimitives
  //   Polyhedron not available for plate5Solid.
  //   This means it cannot be visualized on most systems.
  //   Contact the Visualization Coordinator.
  // So, I don't fillet!
  /*
  G4SubtractionSolid* plateCutSolid = G4QuadPlateWithFillet(l2, h2, r);
  return new G4SubtractionSolid("plate5Solid", plate4Solid, plateCutSolid);
    */
}

//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
// Catania Beamline
//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------

CataniaBeamline::CataniaBeamline(G4VPhysicalVolume* motherVolume, const G4ThreeVector pos) {

  const double myTwoPi = 2 * M_PI * rad;
  G4bool myCheckOverlap = DSStorage::Get()->GetCheckOverlap();
  G4ThreeVector myZeros(0., 0., 0.);

  G4Material* Vacuum = DSMaterial::Get()->GetVacuum();      //  4
  G4Material* Steel = DSMaterial::Get()->GetSteel();        // 14
  G4Material* Aluminum = DSMaterial::Get()->GetAluminum();  // 39
  G4NistManager* manager = G4NistManager::Instance();       // added by Simone on 25 Oct. 2017
  //  G4Material* Mylar = manager->FindOrBuildMaterial("G4_MYLAR");
  //  G4Material* Piombo = manager->FindOrBuildMaterial("G4_Pb");//MetalLead
  G4Material* CH2 = manager->FindOrBuildMaterial("G4_PARAFFIN");

  G4Colour myGray(0.5, 0.5, 0.5);       // gray
  G4Colour mydkGray(0.75, 0.75, 0.75);  // dark gray
  G4Colour myBlue(0.0, 0.0, 1.0);       // blue

  G4VisAttributes* SteelVis = new G4VisAttributes(myGray);
  G4VisAttributes* allVis = new G4VisAttributes(mydkGray);
  G4VisAttributes* targetVis = new G4VisAttributes(myBlue);

  // ---------------------------------------------------//
  //                  ReD Scattering Chamber            //
  // ---------------------------------------------------//

  // Cylindrical Shape Scattering Chamber added by Simone on Feb. 2018
  G4double chamberInnerR = 29.4 * cm;
  G4double chamberOuterR = 30.0 * cm;
  G4double chamberH = 55.0 * cm;
  //  G4ThreeVector chamberPos(-118.496*cm, 31.3385*cm, 48.3399*cm);
  G4ThreeVector chamberPos = pos;
  G4double chamberTBH = 0.6 * cm;
  G4ThreeVector chamberTopPos(0, 0, (chamberH + chamberTBH) / 2.0);

  // the beam tube (in this model) has a length of 300 cm after and 40 cm before
  // the scattering chamber.
  G4double beamtubeLin = 40.0 * cm;
  G4double beamtubeLout = 300.0 * cm;
  G4double beamtubeL = beamtubeLin + beamtubeLout + 2.0 * chamberOuterR;
  G4double beamtubeInnerR = 4.4 * cm;
  G4double beamtubeOuterR = 5.0 * cm;
  // it should stick out 40 cm to the left
  G4ThreeVector beamtubePos(beamtubeL / 2.0 - chamberOuterR - beamtubeLin, 0.0,
                            0.0);  // it should be 20 cm below the top of the chamber

  // G4double Thickness = chamberOuterR - chamberInnerR;
  // G4double innRadius = 0. * cm;
  // G4double outRadius = 5. * cm;
  // G4double innRadius1 = 4.4 * cm;

  // CH2 target
  G4double targetX = 0.0025 * cm;  // spessore 2.5um
  G4double targetY = 4.0 * cm;     // base 4 cm
  G4double targetZ = 4.0 * cm;     // altezza 4 cm
  // MK 180306: right now I put the target at the entrance of the chamber, at a
  // distance of 1cm from the wall
  G4ThreeVector targetPos(-chamberInnerR + 1.0 * cm, 0,
                          0);  // relative to chamber center, to be placed in chamberVac

  // Aluminum beam stopper window
  G4double windowH = 0.4 * cm;
  G4ThreeVector windowPos = chamberPos + G4ThreeVector(chamberOuterR + beamtubeLout + windowH / 2.0, 0, 0);

  G4RotationMatrix* rot = new G4RotationMatrix();
  rot->rotateZ(-7. * deg);
  rot->rotateY(-109. * deg);
  rot->rotateX(-18.8 * deg);

  G4RotationMatrix* rot1 = new G4RotationMatrix();
  rot1->rotateY(-90. * deg);

  // G4double fThetaArX = DSStorage::Get()->GetBeamToTPCAngleThetaX();
  // // 18.8*deg wrt horizzontal G4double fphiArX =
  // DSStorage::Get()->GetBeamToTPCAnglePhiX();              //-12.75*deg wrt
  // vertical G4double fTargetTPCDistance =
  // DSStorage::Get()->GetTargetToTPCDistance();  // -138.496*cm

  // G4double argon_to_target_distance = - 150*cm;
  // G4double fThetaArX = 18.8*deg; //18.8 this is wrt the horizzontal (X axis)
  // G4double fphiArX = -12.75*deg; //-12.75 this is wrt the vertical (X axis)
  // G4double neutrondete_distance = 80*cm;

  // the beamline solid and the beamline vacuum

  G4Tubs* chamberSolid = new G4Tubs("chamber_Solid", 0.0, chamberOuterR, chamberH / 2., 0.0, myTwoPi);
  G4Tubs* chamberVacSolid = new G4Tubs("chamberVac_Solid", 0.0, chamberInnerR, chamberH / 2., 0.0, myTwoPi);

  G4Tubs* beamtubeSolid = new G4Tubs("beamtube_Solid", 0.0, beamtubeOuterR, beamtubeL / 2.0, 0.0, myTwoPi);
  G4Tubs* beamtubeVacSolid = new G4Tubs("beamtubeVac_Solid", 0.0, beamtubeInnerR, beamtubeL / 2.0, 0.0, myTwoPi);

  G4UnionSolid* beamlineSolid = new G4UnionSolid("chamber+beamtube_Solid", chamberSolid, beamtubeSolid, rot1, beamtubePos);
  G4UnionSolid* beamlineVacSolid = new G4UnionSolid("chamberVac+beamtubeVac_Solid", chamberVacSolid, beamtubeVacSolid, rot1, beamtubePos);

  G4LogicalVolume* beamlineLogic = new G4LogicalVolume(beamlineSolid, Steel, "beamline_Logic");
  G4LogicalVolume* beamlineVacLogic = new G4LogicalVolume(beamlineVacSolid, Vacuum, "beamlineVac_Logic");

  G4VPhysicalVolume* beamlinePhysics = new G4PVPlacement(0, chamberPos, "beamline", beamlineLogic, motherVolume, false, 0, myCheckOverlap);
  G4VPhysicalVolume* beamlineVacPhysics = new G4PVPlacement(0, myZeros, "beamlineVac", beamlineVacLogic, beamlinePhysics, false, 0, myCheckOverlap);

  Debug(beamlinePhysics);
  Debug(beamlineVacPhysics);

  // Top and bottom of the scattering chamber

  G4Tubs* chamberTopSolid = new G4Tubs("ChamberTop_Solid", 0., chamberOuterR, chamberTBH / 2.0, 0.0, myTwoPi);
  G4Tubs* chamberBotSolid = new G4Tubs("ChamberBot_Solid", 0., chamberOuterR, chamberTBH / 2.0, 0.0, myTwoPi);

  G4LogicalVolume* chamberTopLogic = new G4LogicalVolume(chamberTopSolid, Aluminum, "ChamberTop_Logic");
  G4LogicalVolume* chamberBotLogic = new G4LogicalVolume(chamberBotSolid, Aluminum, "ChamberBot_Logic");

  G4VPhysicalVolume* fPhysicChamberTop = new G4PVPlacement(0, chamberPos + chamberTopPos, "ChamberTop", chamberTopLogic, motherVolume, false, 0, myCheckOverlap);
  G4VPhysicalVolume* fPhysicChamberBot = new G4PVPlacement(0, chamberPos - chamberTopPos, "ChamberBot", chamberBotLogic, motherVolume, false, 0, myCheckOverlap);

  Debug(fPhysicChamberTop);
  Debug(fPhysicChamberBot);

  // CH2 target

  G4Box* targetSolid = new G4Box("target-Solid", targetX / 2.0, targetY / 2.0, targetZ / 2.0);
  G4LogicalVolume* targetLogic = new G4LogicalVolume(targetSolid, CH2, "target_Logic");
  G4VPhysicalVolume* fPhysicTarget = new G4PVPlacement(0, targetPos, "target", targetLogic, beamlineVacPhysics, false, 0, myCheckOverlap);
  Debug(fPhysicTarget);

  // Aluminum beam stopper window

  G4Tubs* windowSolid = new G4Tubs("window_Solid", 0.0, beamtubeOuterR, windowH / 2.0, 0.0, myTwoPi);
  G4LogicalVolume* windowLogic = new G4LogicalVolume(windowSolid, Aluminum, "window_Logic");
  G4VPhysicalVolume* fPhysicWindow = new G4PVPlacement(rot1, windowPos, "window", windowLogic, motherVolume, false, 0, myCheckOverlap);
  Debug(fPhysicWindow);

  // TODO should go out of this class
  // motherVolume->GetLogicalVolume()->SetVisAttributes(G4VisAttributes::Invisible);

  beamlineLogic->SetVisAttributes(SteelVis);
  chamberTopLogic->SetVisAttributes(allVis);
  chamberBotLogic->SetVisAttributes(allVis);
  targetLogic->SetVisAttributes(targetVis);
  windowLogic->SetVisAttributes(allVis);
}

G4double DSDetectorReD::extended_asin (G4double x) {
  if (x < -1) return - M_PI / 2. ;
  if (x >  1) return M_PI / 2. ;
  return asin (x);
}
