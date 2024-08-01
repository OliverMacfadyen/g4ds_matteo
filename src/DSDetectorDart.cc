/* -------------------------------------- */
/* DART geometry based on DS10 detector*/
/* Olivier Dadoun */
/* January 2017 */
/* -------------------------------------- */
#include <iostream>
#include <stdio.h>

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
#include "G4ThreeVector.hh"
#include "G4UIcommand.hh"

#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include "DSDetectorDart.hh"
#include "DSLogger.hh"
#include "DSMaterial.hh"
#include "DSStorage.hh"

#include "DSParameters.hh"

using namespace std;

DSDetectorDart::DSDetectorDart(G4VPhysicalVolume* myMotherVolume) {

  G4Colour myWhite(1.0, 1.0, 1.0);    // white
  G4Colour myGray(0.5, 0.5, 0.5);     // gray
  G4Colour myBlack(0.0, 0.0, 0.0);    // black
  G4Colour myRed(1.0, 0.0, 0.0);      // red
  G4Colour myGreen(0.0, 1.0, 0.0);    // green
  G4Colour myBlue(0.0, 0.0, 1.0);     // blue
  G4Colour myCyan(0.0, 1.0, 1.0);     // cyan
  G4Colour myMagenta(1.0, 0.0, 1.0);  // magenta
  G4Colour myYellow(1.0, 1.0, 0.0);   // yellow

  fMotherVolume = myMotherVolume;

  const double myTwoPi = 2 * M_PI * rad;
  G4bool myCheckOverlap = DSStorage::Get()->GetCheckOverlap();

  DSLog(routine) << " Constructing Dart Geometry " << endlog;

  // Copper Dewar + LAr bath
  G4double myDewarWall_t = 1 * cm;
  if (DSStorage::Get()->GetDartDewarThickness() != 0.) myDewarWall_t = DSStorage::Get()->GetDartDewarThickness();

  G4double myLArBath_h = 40. * cm;  // h e d sono stati invertiti per coerenza
  if (DSStorage::Get()->GetDartDewarHeight() != 0.) myLArBath_h = DSStorage::Get()->GetDartDewarHeight();

  G4double myLArBath_d = 40. * cm;
  if (DSStorage::Get()->GetDartDewarRadius() != 0.) myLArBath_d = 2. * DSStorage::Get()->GetDartDewarRadius();

  G4ThreeVector myZeros(0., 0.,
                        0.);  // TEMPORARY: At the center of the Mother Volume?

  fSolidDewar = new G4Tubs("Dewar_Solid", 0, myLArBath_d / 2. + myDewarWall_t, myLArBath_h / 2. + myDewarWall_t, 0, myTwoPi);
  fLogicDewar = new G4LogicalVolume(fSolidDewar, DSMaterial::Get()->GetStainlessSteel(), "Dewar_Logic");
  fPhysicDewar = new G4PVPlacement(0, myZeros, "Dewar", fLogicDewar, fMotherVolume, false, 0, myCheckOverlap);

  G4VisAttributes* LogVisAttCopper = new G4VisAttributes(myRed);
  fLogicDewar->SetVisAttributes(LogVisAttCopper);

  fSolidLArBath = new G4Tubs("LArBath_Solid", 0, myLArBath_d / 2., myLArBath_h / 2., 0, myTwoPi);
  fLogicLArBath = new G4LogicalVolume(fSolidLArBath, DSMaterial::Get()->GetNSLiquidArgon(), "LArBath_Logic");
  fPhysicLArBath = new G4PVPlacement(0, myZeros, "LArBath", fLogicLArBath, fPhysicDewar, false, 0, myCheckOverlap);

  G4VisAttributes* LogVisAttBath = new G4VisAttributes(myBlue);
  fLogicLArBath->SetVisAttributes(LogVisAttBath);

  G4double myTeflonWall_t = 1. * cm;
  if (DSStorage::Get()->GetDartTeflonThickness() != 0.) myTeflonWall_t = DSStorage::Get()->GetDartTeflonThickness();

  G4double myLArScint_h = 10. * cm;
  if (DSStorage::Get()->GetDartTeflonHeight() != 0.) myLArScint_h = DSStorage::Get()->GetDartTeflonHeight();

  G4double myLArScint_d = 10. * cm;
  if (DSStorage::Get()->GetDartTeflonRadius() != 0.) myLArScint_d = 2. * DSStorage::Get()->GetDartTeflonRadius();

  fSolidTeflon = new G4Tubs("Teflon_Solid", 0, myLArScint_d / 2. + myTeflonWall_t, myLArScint_h / 2. + myTeflonWall_t, 0, myTwoPi);
  fLogicTeflon = new G4LogicalVolume(fSolidTeflon, DSMaterial::Get()->GetTeflon(), "Teflon_Logic");
  fPhysicTeflon = new G4PVPlacement(0, myZeros, "Teflon", fLogicTeflon, fPhysicLArBath, false, 0, myCheckOverlap);

  G4VisAttributes* LogVisAttTeflon = new G4VisAttributes(myWhite);
  fLogicTeflon->SetVisAttributes(LogVisAttTeflon);
  //-----------------------//
  //    Copper Rings        //
  //-----------------------//
  G4double RingInnerR = myLArScint_d / 2. + myTeflonWall_t + 0.1 * cm;
  G4double RingOuterR = RingInnerR + 0.02 * cm;
  G4double RingH = 0.05 * cm;
  // G4double SleeveH = GetCm(3.00/2)*cm ;
  G4double RingSpace = 1. * cm;
  G4Tubs* fSolidCopperRing = new G4Tubs("CopperRing_Solid", RingInnerR, RingOuterR, RingH, 0, 2 * M_PI);
  G4LogicalVolume* fLogicCopperRing = new G4LogicalVolume(fSolidCopperRing, DSMaterial::Get()->GetMetalCopper(), "CopperRing_Logic");
  // define 20 rings
  G4ThreeVector RingPlacement(0., 0., 0.);
  char ring_name[256];
  //G4VPhysicalVolume *fRingCopper[20];

  for (int j = 0; j < 10; j++) {
    snprintf(ring_name, 30, "CopperRing_%d", j);
    RingPlacement.setZ(-myLArScint_h / 2 + RingSpace * j);
    // RingPlacement = RingPlacement - LArTPCShift;
    //fRingCopper[j] =
    new G4PVPlacement(0, RingPlacement, ring_name, fLogicCopperRing, fPhysicLArBath, false, 0, myCheckOverlap);
  }
  G4VisAttributes* LogVisAttScint = new G4VisAttributes(myRed);
  fLogicCopperRing->SetVisAttributes(LogVisAttScint);
  /* */
  fSolidLArScint = new G4Tubs("LArScint_Solid", 0, myLArScint_d / 2., myLArScint_h / 2., 0, myTwoPi);
  fLogicLArScint = new G4LogicalVolume(fSolidLArScint, DSMaterial::Get()->GetLiquidArgon(), "LArScint_Logic");
  // Named ActiveLAr mandatory for energy purpose
  fPhysicLArScint = new G4PVPlacement(0, myZeros, "ActiveLAr", fLogicLArScint, fPhysicTeflon, false, 0, myCheckOverlap);

  fLogicLArScint->SetVisAttributes(LogVisAttScint);

  // Gas Pocket and Inner LAr
  /*G4double myGAr_h = 2.0 * cm;
  G4ThreeVector myGasPocketPos( 0, 0,  myLArScint_h/2 - myGAr_h/2. );
  fSolidGasPocket  = new G4Tubs("GasPocket_Solid", 0, myLArScint_d/2.,
  myGAr_h/2., 0, myTwoPi ); fLogicGasPocket  = new
  G4LogicalVolume(fSolidGasPocket, DSMaterial::Get()->GetGaseousArgon(),
  "GasPocket_Logic" ); fPhysicGasPocket = new G4PVPlacement( 0, myGasPocketPos,
  "GasPocket", fLogicGasPocket, fPhysicLArScint, false, 0, myCheckOverlap );

  G4VisAttributes* LogVisAttGas= new G4VisAttributes(myMagenta);
  fLogicGasPocket->SetVisAttributes(LogVisAttGas);
   */

  // Set the z coordinate of the LAr - GAr interface, necessary for S2
  // generation in DSLightX
  // G4double myLArGArBoundaryPosZ = ( myGasPocketPos + myInnerVesselWallPos +
  // myLArBathPos ).z() - myGAr_h/2.; G4ThreeVector myGasPocketPos( 0, 0, 0.  );
  G4ThreeVector myInnerVesselWallPos(0, 0, 0);
  G4double myLArGArBoundaryPosZ = 0;
  DSStorage::Get()->SetLArGArBoundaryPosZ(myLArGArBoundaryPosZ);

  // The definition of the LAr region here is needed to set the range cuts for
  // this volume to a smaller value with respect to the rest of the detector (
  // see DSPhysicsList::SetCuts() )
  G4Region* fLArRegion = new G4Region("LAr_Logic");
  fLogicLArScint->SetRegion(fLArRegion);
  fLArRegion->AddRootLogicalVolume(fLogicLArScint);

  DefineSurfaces();
}

DSDetectorDart::~DSDetectorDart() {
  ;
}

void DSDetectorDart::DefineSurfaces() {

  ////////////////////////////////////////
  // GAr - LAr
  ////////////////////////////////////////
  G4OpticalSurface* fOpGArLArSurface = new G4OpticalSurface("OpGArLArSurface");
  new G4LogicalBorderSurface("GArLArSurface", fPhysicGasPocket, fPhysicInnerLAr, fOpGArLArSurface);
  fOpGArLArSurface->SetType(dielectric_dielectric);
  fOpGArLArSurface->SetModel(unified);
  fOpGArLArSurface->SetFinish(polished);
  G4MaterialPropertiesTable* fGArLArSurfProp = new G4MaterialPropertiesTable();
  fOpGArLArSurface->SetMaterialPropertiesTable(fGArLArSurfProp);

  ////////////////////////////////////////
  // LAr - StainlessSteel (supportRod + compression plates)
  ////////////////////////////////////////
  G4OpticalSurface* fOpLArSteelSurface = new G4OpticalSurface("OpLArSteelSurface");
  new G4LogicalBorderSurface("LArSteelSurfaceSupp0", fPhysicLArBath, fPhysicSupportRod0, fOpLArSteelSurface);
  new G4LogicalBorderSurface("LArSteelSurfaceSupp1", fPhysicLArBath, fPhysicSupportRod1, fOpLArSteelSurface);
  new G4LogicalBorderSurface("LArSteelSurfaceSupp2", fPhysicLArBath, fPhysicSupportRod2, fOpLArSteelSurface);
  new G4LogicalBorderSurface("LArSteelSurfaceSupp3", fPhysicLArBath, fPhysicSupportRod3, fOpLArSteelSurface);
  new G4LogicalBorderSurface("LArSteelSurfacePlateTop", fPhysicLArBath, fPhysicCompressionPlateTop, fOpLArSteelSurface);
  new G4LogicalBorderSurface("LArSteelSurfacePlateBot", fPhysicLArBath, fPhysicCompressionPlateBottom, fOpLArSteelSurface);
  fOpLArSteelSurface->SetType(dielectric_metal);
  fOpLArSteelSurface->SetModel(unified);
  fOpLArSteelSurface->SetFinish(polished);
  fOpLArSteelSurface->SetSigmaAlpha(0.1);
  G4MaterialPropertiesTable* fLArSteelSurfProp = new G4MaterialPropertiesTable();
  G4double LArSteelENE[2] = {0.1 * eV, 20.0 * eV};
  G4double LArSteelREF[2] = {0.99, 0.99};
  G4double LArSteelEFF[2] = {0.00, 0.00};
  fLArSteelSurfProp->AddProperty("REFLECTIVITY", LArSteelENE, LArSteelREF, 2);
  fLArSteelSurfProp->AddProperty("EFFICIENCY", LArSteelENE, LArSteelEFF, 2);
  fOpLArSteelSurface->SetMaterialPropertiesTable(fLArSteelSurfProp);

  ////////////////////////////////////////
  // LAr - GridSteel
  ////////////////////////////////////////
  G4OpticalSurface* fOpLArGridSurface = new G4OpticalSurface("OpLArGridSurface");
  new G4LogicalBorderSurface("LArGridSurfacePlateBot", fPhysicInnerLAr, fPhysicGrid, fOpLArGridSurface);
  fOpLArGridSurface->SetType(dielectric_dielectric);
  fOpLArGridSurface->SetModel(unified);
  fOpLArGridSurface->SetFinish(polished);
  // fOpLArGridSurface->SetSigmaAlpha(0.1);
  G4MaterialPropertiesTable* fLArGridSurfProp = new G4MaterialPropertiesTable();
  // G4double LArGridENE[2] = {0.1*eV, 20.0*eV};
  // G4double LArGridREF[2] = {0.11, 0.11};   // 89% absorption
  // G4double LArGridEFF[2] = {0.00, 0.00};
  // fLArGridSurfProp->AddProperty("REFLECTIVITY", LArGridENE, LArGridREF, 2);
  // fLArGridSurfProp->AddProperty("EFFICIENCY",   LArGridENE, LArGridEFF, 2);
  fOpLArGridSurface->SetMaterialPropertiesTable(fLArGridSurfProp);

  ////////////////////////////////////////
  // LArBath - Copper (dewar)
  ////////////////////////////////////////
  G4OpticalSurface* fOpLArCopperSurface = new G4OpticalSurface("OpLArCopperSurface");
  new G4LogicalBorderSurface("LArCopperSurface", fPhysicLArBath, fPhysicDewar, fOpLArCopperSurface);
  fOpLArCopperSurface->SetType(dielectric_metal);
  fOpLArCopperSurface->SetModel(unified);
  fOpLArCopperSurface->SetFinish(polished);
  fOpLArCopperSurface->SetSigmaAlpha(0.1);
  G4MaterialPropertiesTable* fLArCopperSurfProp = new G4MaterialPropertiesTable();
  G4double LArCopperENE[2] = {0.1 * eV, 20.0 * eV};
  G4double LArCopperREF[2] = {0.99, 0.99};
  G4double LArCopperEFF[2] = {0.00, 0.00};
  fLArCopperSurfProp->AddProperty("REFLECTIVITY", LArCopperENE, LArCopperREF, 2);
  fLArCopperSurfProp->AddProperty("EFFICIENCY", LArCopperENE, LArCopperEFF, 2);
  fOpLArCopperSurface->SetMaterialPropertiesTable(fLArCopperSurfProp);

  ////////////////////////////////////////
  // TPB - GAr
  ////////////////////////////////////////
  G4OpticalSurface* fOpTPBGArSurface = new G4OpticalSurface("OpTPBGArSurface");
  new G4LogicalBorderSurface("TPBGArSurface", fPhysicGasPocket, fPhysicTPBTopLateral, fOpTPBGArSurface);
  new G4LogicalBorderSurface("TPBGArSurface", fPhysicGasPocket, fPhysicTPBTopBase, fOpTPBGArSurface);
  fOpTPBGArSurface->SetType(dielectric_dielectric);
  fOpTPBGArSurface->SetModel(unified);
  fOpTPBGArSurface->SetFinish(ground);
  G4MaterialPropertiesTable* fTPBGArSurfProp = new G4MaterialPropertiesTable();
  fOpTPBGArSurface->SetMaterialPropertiesTable(fTPBGArSurfProp);

  ////////////////////////////////////////
  // LAr(bath) - Kapton
  ////////////////////////////////////////
  G4OpticalSurface* fOpLArKaptonSurface = new G4OpticalSurface("OpLArKaptonSurface");
  new G4LogicalBorderSurface("LArKaptonSurface", fPhysicLArBath, fPhysicKaptonBand, fOpLArKaptonSurface);
  fOpLArKaptonSurface->SetType(dielectric_metal);
  fOpLArKaptonSurface->SetModel(unified);
  fOpLArKaptonSurface->SetFinish(polished);
  fOpLArKaptonSurface->SetSigmaAlpha(0.1);
  G4MaterialPropertiesTable* fLArKaptonSurfProp = new G4MaterialPropertiesTable();
  G4double LArKaptonENE[2] = {0.1 * eV, 20.0 * eV};
  G4double LArKaptonREF[2] = {0.99, 0.99};
  G4double LArKaptonEFF[2] = {0.00, 0.00};
  fLArKaptonSurfProp->AddProperty("REFLECTIVITY", LArKaptonENE, LArKaptonREF, 2);
  fLArKaptonSurfProp->AddProperty("EFFICIENCY", LArKaptonENE, LArKaptonEFF, 2);
  fOpLArKaptonSurface->SetMaterialPropertiesTable(fLArKaptonSurfProp);

  ////////////////////////////////////////
  // TPB - Acrylic(inner vessel)
  ////////////////////////////////////////
  G4OpticalSurface* fOpTPBVesselSurface = new G4OpticalSurface("OpTPBVesselSurface");
  new G4LogicalBorderSurface("TPBVesselSurfaceTop", fPhysicTPBTopBase, fPhysicInnerVesselWindowTop, fOpTPBVesselSurface);
  new G4LogicalBorderSurface("TPBVesselSurfaceBot", fPhysicTPBBottomBase, fPhysicInnerVesselWindowBottom, fOpTPBVesselSurface);
  new G4LogicalBorderSurface("TPBVesselSurfaceTopLat", fPhysicTPBTopLateral, fPhysicInnerVesselWall, fOpTPBVesselSurface);
  new G4LogicalBorderSurface("TPBVesselSurfaceMidLat", fPhysicTPBMiddleLateral, fPhysicInnerVesselWall, fOpTPBVesselSurface);
  new G4LogicalBorderSurface("TPBVesselSurfaceLowLat", fPhysicTPBLowerLateral, fPhysicInnerVesselWall, fOpTPBVesselSurface);
  fOpTPBVesselSurface->SetType(dielectric_dielectric);
  fOpTPBVesselSurface->SetModel(unified);
  fOpTPBVesselSurface->SetFinish(ground);
  G4MaterialPropertiesTable* fTPBVesselSurfProp = new G4MaterialPropertiesTable();
  fOpTPBVesselSurface->SetMaterialPropertiesTable(fTPBVesselSurfProp);

  ////////////////////////////////////////
  // LAr - Acrylic(inner vessel)
  ////////////////////////////////////////
  G4OpticalSurface* fOpLArVesselSurface = new G4OpticalSurface("OpLArVesselSurface");
  new G4LogicalBorderSurface("LArVesselWall", fPhysicInnerLAr, fPhysicInnerVesselWall, fOpLArVesselSurface);
  new G4LogicalBorderSurface("LArVesselSurfaceTop", fPhysicLArBath, fPhysicInnerVesselWindowTop, fOpLArVesselSurface);
  new G4LogicalBorderSurface("LArVesselSurfaceWall", fPhysicLArBath, fPhysicInnerVesselWindowBottom, fOpLArVesselSurface);
  fOpLArVesselSurface->SetType(dielectric_dielectric);
  fOpLArVesselSurface->SetModel(unified);
  fOpLArVesselSurface->SetFinish(ground);
  G4MaterialPropertiesTable* fLArVesselSurfProp = new G4MaterialPropertiesTable();
  fOpLArVesselSurface->SetMaterialPropertiesTable(fLArVesselSurfProp);

  ////////////////////////////////////////
  // GAr - Acrylic(inner vessel)
  ////////////////////////////////////////
  G4OpticalSurface* fOpGArVesselSurface = new G4OpticalSurface("OpGArVesselSurface");
  new G4LogicalBorderSurface("GArVessel", fPhysicGasPocket, fPhysicInnerVesselWall, fOpGArVesselSurface);
  fOpGArVesselSurface->SetType(dielectric_dielectric);
  fOpGArVesselSurface->SetModel(unified);
  fOpGArVesselSurface->SetFinish(ground);
  G4MaterialPropertiesTable* fGArVesselSurfProp = new G4MaterialPropertiesTable();
  fOpGArVesselSurface->SetMaterialPropertiesTable(fGArVesselSurfProp);

  ////////////////////////////////////////
  // TPB - 3Mfoil (reflector)
  ////////////////////////////////////////
  G4OpticalSurface* fOpTPBThreeMSurface = new G4OpticalSurface("OpTPBThreeMSurface");
  new G4LogicalBorderSurface("TPBThreeMSurfaceTop", fPhysicTPBTopLateral, fPhysicThreeMTopFoil, fOpTPBThreeMSurface);
  new G4LogicalBorderSurface("TPBThreeMSurfaceMid", fPhysicTPBMiddleLateral, fPhysicThreeMMiddleFoil, fOpTPBThreeMSurface);
  new G4LogicalBorderSurface("TPBThreeMSurfaceBot", fPhysicTPBLowerLateral, fPhysicThreeMLowerFoil, fOpTPBThreeMSurface);
  fOpTPBThreeMSurface->SetType(dielectric_metal);
  fOpTPBThreeMSurface->SetModel(unified);
  fOpTPBThreeMSurface->SetFinish(ground);
  fOpTPBThreeMSurface->SetSigmaAlpha(0.5);
  G4double TPBThreeMENE[2] = {0.1 * eV, 20.0 * eV};
  G4double TPBThreeMREF[2] = {0.99, 0.99};
  G4double TPBThreeMEFF[2] = {0.00, 0.00};
  G4MaterialPropertiesTable* fTPBThreeMSurfProp = new G4MaterialPropertiesTable();
  fTPBThreeMSurfProp->AddProperty("REFLECTIVITY", TPBThreeMENE, TPBThreeMREF, 2);
  fTPBThreeMSurfProp->AddProperty("EFFICIENCY", TPBThreeMENE, TPBThreeMEFF, 2);
  fOpTPBThreeMSurface->SetMaterialPropertiesTable(fTPBThreeMSurfProp);

  ////////////////////////////////////////
  // LAr - 3Mfoil (reflector)
  ////////////////////////////////////////
  G4OpticalSurface* fOpLArThreeMSurface = new G4OpticalSurface("OpLArThreeMSurface");
  new G4LogicalBorderSurface("LArThreeMSurfaceTop", fPhysicInnerLAr, fPhysicThreeMTopFoil, fOpLArThreeMSurface);
  new G4LogicalBorderSurface("LArThreeMSurfaceMid", fPhysicInnerLAr, fPhysicThreeMMiddleFoil, fOpLArThreeMSurface);
  new G4LogicalBorderSurface("LArThreeMSurfaceBot", fPhysicInnerLAr, fPhysicThreeMLowerFoil, fOpLArThreeMSurface);
  fOpLArThreeMSurface->SetType(dielectric_metal);
  fOpLArThreeMSurface->SetModel(unified);
  fOpLArThreeMSurface->SetFinish(ground);
  fOpLArThreeMSurface->SetSigmaAlpha(0.5);
  G4double LArThreeMENE[2] = {0.1 * eV, 20.0 * eV};
  G4double LArThreeMREF[2] = {0.97, 0.97};
  G4double LArThreeMEFF[2] = {0.00, 0.00};
  G4MaterialPropertiesTable* fLArThreeMSurfProp = new G4MaterialPropertiesTable();
  fLArThreeMSurfProp->AddProperty("REFLECTIVITY", LArThreeMENE, LArThreeMREF, 2);
  fLArThreeMSurfProp->AddProperty("EFFICIENCY", LArThreeMENE, LArThreeMEFF, 2);
  fOpLArThreeMSurface->SetMaterialPropertiesTable(fLArThreeMSurfProp);

  ////////////////////////////////////////
  // Teflon - 3Mfoil   [NOT NEEDED ?]
  ////////////////////////////////////////
  G4OpticalSurface* fOpTeflonThreeMSurface = new G4OpticalSurface("OpTeflonThreeMSurface");
  new G4LogicalBorderSurface("TeflonThreeMSurfaceGas", fPhysicThreeMTopFoil, fPhysicTopGasSupportRing, fOpTeflonThreeMSurface);
  new G4LogicalBorderSurface("TeflonThreeMSurfaceLiq", fPhysicThreeMTopFoil, fPhysicTopLiqSupportRing, fOpTeflonThreeMSurface);
  new G4LogicalBorderSurface("TeflonThreeMSurfaceR0", fPhysicThreeMLowerFoil, fPhysicSupportRing0,
                             fOpTeflonThreeMSurface);  // da verificare !!!
  new G4LogicalBorderSurface("TeflonThreeMSurfaceR1", fPhysicThreeMMiddleFoil, fPhysicSupportRing1,
                             fOpTeflonThreeMSurface);  // da verificare !!!
  new G4LogicalBorderSurface("TeflonThreeMSurfaceR2", fPhysicThreeMTopFoil, fPhysicSupportRing2,
                             fOpTeflonThreeMSurface);  // da verificare !!!
  fOpTeflonThreeMSurface->SetType(dielectric_metal);
  fOpTeflonThreeMSurface->SetModel(unified);
  fOpTeflonThreeMSurface->SetFinish(ground);
  fOpTeflonThreeMSurface->SetSigmaAlpha(0.1);
  G4double TeflonThreeMENE[2] = {0.1 * eV, 20.0 * eV};
  G4double TeflonThreeMREF[2] = {0.99, 0.99};
  G4double TeflonThreeMEFF[2] = {0.00, 0.00};
  G4MaterialPropertiesTable* fTeflonThreeMSurfProp = new G4MaterialPropertiesTable();
  fTeflonThreeMSurfProp->AddProperty("REFLECTIVITY", TeflonThreeMENE, TeflonThreeMREF, 2);
  fTeflonThreeMSurfProp->AddProperty("EFFICIENCY", TeflonThreeMENE, TeflonThreeMEFF, 2);
  fOpTeflonThreeMSurface->SetMaterialPropertiesTable(fTeflonThreeMSurfProp);

  ////////////////////////////////////////
  // Acrylic - Copper
  ////////////////////////////////////////
  G4OpticalSurface* fOpAcrylicCopperSurface = new G4OpticalSurface("OpAcrylicCopperSurface");
  new G4LogicalBorderSurface("AcrylicCopperSurface", fPhysicInnerVesselWall, fPhysicFieldRings, fOpAcrylicCopperSurface);
  fOpAcrylicCopperSurface->SetType(dielectric_metal);
  fOpAcrylicCopperSurface->SetModel(unified);
  fOpAcrylicCopperSurface->SetFinish(ground);
  fOpAcrylicCopperSurface->SetSigmaAlpha(0.1);
  G4double AcrylicCopperENE[2] = {0.1 * eV, 20.0 * eV};
  G4double AcrylicCopperREF[2] = {0.99, 0.99};
  G4double AcrylicCopperEFF[2] = {0.00, 0.00};
  G4MaterialPropertiesTable* fAcrylicCopperSurfProp = new G4MaterialPropertiesTable();
  fAcrylicCopperSurfProp->AddProperty("REFLECTIVITY", AcrylicCopperENE, AcrylicCopperREF, 2);
  fAcrylicCopperSurfProp->AddProperty("EFFICIENCY", AcrylicCopperENE, AcrylicCopperEFF, 2);
  fOpAcrylicCopperSurface->SetMaterialPropertiesTable(fAcrylicCopperSurfProp);

  ////////////////////////////////////////
  // LAr - Teflon
  ////////////////////////////////////////
  G4OpticalSurface* fOpLArTeflonSurface = new G4OpticalSurface("OpLArTeflonSurface");
  new G4LogicalBorderSurface("LArTeflonSurfaceTop", fPhysicInnerLAr, fPhysicTopLiqSupportRing, fOpLArTeflonSurface);
  new G4LogicalBorderSurface("LArTeflonSurfaceR0", fPhysicInnerLAr, fPhysicSupportRing0, fOpLArTeflonSurface);
  new G4LogicalBorderSurface("LArTeflonSurfaceR1", fPhysicInnerLAr, fPhysicSupportRing1, fOpLArTeflonSurface);
  new G4LogicalBorderSurface("LArTeflonSurfaceR2", fPhysicInnerLAr, fPhysicSupportRing2, fOpLArTeflonSurface);
  new G4LogicalBorderSurface("LArTeflonSurfaceBub1", fPhysicLArBath, fPhysicBubblerTube1, fOpLArTeflonSurface);
  new G4LogicalBorderSurface("LArTeflonSurfaceBub2", fPhysicLArBath, fPhysicBubblerTube2, fOpLArTeflonSurface);
  fOpLArTeflonSurface->SetType(dielectric_metal);
  fOpLArTeflonSurface->SetModel(unified);
  fOpLArTeflonSurface->SetFinish(polished);
  fOpLArTeflonSurface->SetSigmaAlpha(0.1);
  G4MaterialPropertiesTable* fLArTeflonSurfProp = new G4MaterialPropertiesTable();
  G4double LArTeflonENE[2] = {0.1 * eV, 20.0 * eV};
  G4double LArTeflonREF[2] = {0.99, 0.99};
  G4double LArTeflonEFF[2] = {0.00, 0.00};
  fLArTeflonSurfProp->AddProperty("REFLECTIVITY", LArTeflonENE, LArTeflonREF, 2);
  fLArTeflonSurfProp->AddProperty("EFFICIENCY", LArTeflonENE, LArTeflonEFF, 2);
  fOpLArTeflonSurface->SetMaterialPropertiesTable(fLArTeflonSurfProp);

  ////////////////////////////////////////
  // GAr - Teflon (reflector)
  ////////////////////////////////////////
  G4OpticalSurface* fOpGArTeflonSurface = new G4OpticalSurface("OpGArTeflonSurface");
  new G4LogicalBorderSurface("GArTeflonSurfaceTop", fPhysicGasPocket, fPhysicTopGasSupportRing, fOpGArTeflonSurface);
  fOpGArTeflonSurface->SetType(dielectric_metal);
  fOpGArTeflonSurface->SetModel(unified);
  fOpGArTeflonSurface->SetFinish(polished);
  fOpGArTeflonSurface->SetSigmaAlpha(0.1);
  G4MaterialPropertiesTable* fGArTeflonSurfProp = new G4MaterialPropertiesTable();
  G4double GArTeflonENE[2] = {0.1 * eV, 20.0 * eV};
  G4double GArTeflonREF[2] = {0.99, 0.99};
  G4double GArTeflonEFF[2] = {0.00, 0.00};
  fGArTeflonSurfProp->AddProperty("REFLECTIVITY", GArTeflonENE, GArTeflonREF, 2);
  fGArTeflonSurfProp->AddProperty("EFFICIENCY", GArTeflonENE, GArTeflonEFF, 2);
  fOpGArTeflonSurface->SetMaterialPropertiesTable(fGArTeflonSurfProp);

  ////////////////////////////////////////
  // ITO /////
  ////////////////////////////////////////

  G4MaterialPropertiesTable* fITOSurfProp = new G4MaterialPropertiesTable();
  fITOSurfProp->AddConstProperty("DOITO", 1);

  ////////////////////////////////////////
  // LAr (inner) - TPB
  ////////////////////////////////////////
  G4OpticalSurface* fOpTPBLArSurface = new G4OpticalSurface("OpTBPLArSurface");
  new G4LogicalBorderSurface("TPBLArSurfaceInnerBotBas", fPhysicInnerLAr, fPhysicTPBBottomBase, fOpTPBLArSurface);
  new G4LogicalBorderSurface("TPBLArSurfaceInnerLowLat", fPhysicInnerLAr, fPhysicTPBLowerLateral, fOpTPBLArSurface);
  new G4LogicalBorderSurface("TPBLArSurfaceInnerMidLat", fPhysicInnerLAr, fPhysicTPBMiddleLateral, fOpTPBLArSurface);
  new G4LogicalBorderSurface("TPBLArSurfaceInnerTopLat", fPhysicInnerLAr, fPhysicTPBTopLateral,
                             fOpTPBLArSurface);  // not needed ?
  new G4LogicalBorderSurface("TPBLArSurfaceInnerTopBas", fPhysicInnerLAr, fPhysicTPBTopBase,
                             fOpTPBLArSurface);  // not needed ?
  fOpTPBLArSurface->SetType(dielectric_dielectric);
  fOpTPBLArSurface->SetModel(unified);
  fOpTPBLArSurface->SetFinish(ground);
  fOpTPBLArSurface->SetSigmaAlpha(0.1);
  G4MaterialPropertiesTable* fTPBLArSurfProp = new G4MaterialPropertiesTable();
  fOpTPBLArSurface->SetMaterialPropertiesTable(fTPBLArSurfProp);

  ////////////////////////////////////////
  // PMT surfaces
  ////////////////////////////////////////

  // Photocathode - LAr (bath)
  G4OpticalSurface* fOpPMTLArSurface = new G4OpticalSurface("OpPMTLArSurface");
  //G4LogicalBorderSurface *fPMTLArSurface[14];
  G4MaterialPropertiesTable* fPMTLArSurfProp;
  for (int i = 0; i < 14; i++)
    //fPMTLArSurface[i] =
    new G4LogicalBorderSurface("PMTLArSurface", fPhysicLArBath, fPhysicPMTWindow[i], fOpPMTLArSurface);
  fOpPMTLArSurface->SetType(dielectric_dielectric);
  fOpPMTLArSurface->SetModel(unified);
  fOpPMTLArSurface->SetFinish(polished);
  fPMTLArSurfProp = new G4MaterialPropertiesTable();
  fPMTLArSurfProp->AddConstProperty("REFLECTIVITY", 0.2);
  fPMTLArSurfProp->AddConstProperty("EFFICIENCY", 0.0);
  fOpPMTLArSurface->SetMaterialPropertiesTable(fPMTLArSurfProp);

  // LAr (bath) - Kovar (PMT body)
  G4OpticalSurface* fOpLArKovarSurface = new G4OpticalSurface("OpLArKovarSurface");
  for (int i = 0; i < 14; i++) new G4LogicalBorderSurface("LArKovarSurface", fPhysicLArBath, fPhysicPMTBody[i], fOpLArKovarSurface);
  fOpLArKovarSurface->SetType(dielectric_metal);
  fOpLArKovarSurface->SetModel(unified);
  fOpLArKovarSurface->SetFinish(polished);
  fOpLArKovarSurface->SetSigmaAlpha(0.1);
  G4MaterialPropertiesTable* fLArKovarSurfProp = new G4MaterialPropertiesTable();
  G4double LArKovarENE[2] = {0.1 * eV, 20.0 * eV};
  G4double LArKovarREF[2] = {1.00, 1.00};
  G4double LArKovarEFF[2] = {0.00, 0.00};
  fLArKovarSurfProp->AddProperty("REFLECTIVITY", LArKovarENE, LArKovarREF, 2);
  fLArKovarSurfProp->AddProperty("EFFICIENCY", LArKovarENE, LArKovarEFF, 2);
}
