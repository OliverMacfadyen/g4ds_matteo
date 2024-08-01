#include <iostream>
#include "G4Box.hh"
#include "G4Orb.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
#include "G4ThreeVector.hh"

#include "DSDetectorProto.hh"
#include "DSIO.hh"
#include "DSLogger.hh"
#include "DSMaterial.hh"
#include "DSParameters.hh"
#include "DSStorage.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

using namespace std;

////////////////////////////////////////////////
//////        Detector Description      ////////
////////////////////////////////////////////////

// 24th April 2015
// schematics of the 20k detector, starting from scaling the 5K geometry (the
// cryostats in particular) and from some drawings from Hanguo OPTICS is
// directly derived from DS50
//
// 9th October 2015
// change in the cryostats shape: we want to scale the tpc
//
// Placements are realized with no shifts. The correct position is calculated
// when making the solids.
//
// GENERAL SCHEME:
// - nested single volumes, from the outside to the center: Vacuum Jacket
// (steel), Vacuum, Inner Cryo (ssteel).
// - the fill of the inner Volume is divided in two region at z = 1200, which
// correspond to the LAr-GAr interface (Inactive LAr and InactiveGar)
//   - an additional layer of LAr is placed above 1200 cm (default is 5 cm)
//   - the TPC lays inside the InactiveLAr region (Copper Rings, Reflector,
//   TPBSide , Active LAr)
//   - the Gas Pocket lays in the upper part (reflector, Active GAr)
// - other volumes:
//   - the grid is placed 0.5 cm below the LAr surface
//   - Fused silica windows with SiPm on the outfacing surface and TPB on the
//   innerfacing one
//   - Acrylic windows (2 x 1.5) + sustaining plate (2.5 cm)
//   - Pseudo Argon Laryer (condensation on the top window)  - - - REMOVED
//   -

// addition note on TPBSide and TPBTop
// TPBSide is a volume sorrunding the LAr Active volume on the side and bottom.
// TPBTop is only between LArLAyer (Pseudo Ar) and the top window

// Oct'17
// - fixed dimensions of acrylic windows.
// - removed old cryostat (useless)

DSDetectorProto::DSDetectorProto(G4VPhysicalVolume* myMotherVolume) {

  bool IsAlternateDesign = false;
  // if (DSStorage::Get()->Get20KAlternateVeto()) IsAlternateDesign = true ;

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
  DSLog(routine) << " Constructing DS20k Geometry" << endlog;

  const double myTwoPi = 2 * M_PI * rad;

  G4RotationMatrix* myDefaultRotation2 = new G4RotationMatrix();
  G4double myDefaultAngle = myTwoPi / 16.;
  myDefaultRotation2->rotateZ(-myDefaultAngle);  // to correct for octagonal symmetry

  G4bool myCheckOverlap = DSStorage::Get()->GetCheckOverlap();
  G4ThreeVector myZeros(0., 0., 0.);

  // Parameters from drawings
  G4double myTPCHeight = DSStorage::Get()->GetDS20kTPCheight();  // default 120
                                                                 // cm
  G4double myTPCEdge = DSStorage::Get()->GetDS20kTPCedge();      // default 240 cm
  G4double myPMDSuppHeight = 7 * cm;                             // 0.03*myTPCEdge+3*cm ;
  G4double myPMDSuppEdge = myTPCEdge;                            // +5*cm;

  // Reflector
  double myReflectorThickness = 2.54 * cm;  // 1.0 inches

  // reduce thickness and change material if acrylic reflector
  bool IsAcrylicReflector = true;

  // build new (09/2018) baseline with no Cu vessel
  bool isNoCuVessel = true;

  if (IsAcrylicReflector) myReflectorThickness = 5 * mm;
  if (isNoCuVessel) {
    myReflectorThickness = DSStorage::Get()->GetDS20kTPCVesselThickness();
    if (DSStorage::Get()->GetDS20kCryoWallThick() > 0) myReflectorThickness = DSStorage::Get()->GetDS20kCryoWallThick();
  }

  // reduce thickness if DSProto
  if (myTPCEdge < 60 * cm) myReflectorThickness = 0.5 * cm;

  double myRingThickness = 0. * mm;  // 0.5*cm;
  double myTPCInnerRadius = GetOctagonInnerRadius(myTPCEdge);
  double myTeflonInnerRadius = myTPCInnerRadius + myReflectorThickness;
  double myRingOuterRadius = myTeflonInnerRadius + myRingThickness;

  // place the center of the TPC in the center of the ref system
  double myLArGArBoundaryPosZ = DSStorage::Get()->GetDS20kTPCheight() / 2.;
  // Set the z coordinate of the LAr - GAr interface, necessary for S2
  // generation in DSLightX
  DSStorage::Get()->SetLArGArBoundaryPosZ(myLArGArBoundaryPosZ + 1.0 * um);

  // Other parameters
  double myTPBThickness = 0.1 * mm;
  double myPArThickness = 0. * mm;
  double mySiPMOffset = DSStorage::Get()->GetSiPMOffset();  // default: 5 cm
  double mySiPmThickness = 1. * mm;
  double mySiPmBoardThickness = 5.0 * mm;
  double myAcrylicBoardThickness = 6.0 * mm;
  double myTopAcrylicThickness = 1.5 * cm;
  double myGasPocketThickness = 1 * cm;
  double myActiveGasLogicThickness = myGasPocketThickness;  // + myWindowsThickness + mySiPmThickness +
                                                            // mySiPmBoardThickness;

  // DS20k custom cryo specific
  // double myGdLayerThickness = DSStorage::Get()->GetGdLayerThickness();
  // //default 0 double myCryoToCornerDistance =
  // DSStorage::Get()->GetDS20kCryoCornerDistance();  //default 5 cm.

  // average steel thickness from SandroDeCecco
  double myCryoWallThickness = DSStorage::Get()->GetDS20kCryoWallThick();
  double myCryoWallThicknessI = DSStorage::Get()->GetDS20kCryoWallThick();
  if (myCryoWallThickness == 0 * cm) {
    myCryoWallThicknessI = 1.25 * cm;
    myCryoWallThickness = 1.75 * cm;
  }

  // cryostat material
  G4Material* myCryoMat = DSMaterial::Get()->GetStainlessSteel();
  if (DSStorage::Get()->GetDS20kCryoMaterial() == 1) myCryoMat = DSMaterial::Get()->GetMetalTitanium();
  else if (DSStorage::Get()->GetDS20kCryoMaterial() == 2)
    myCryoMat = DSMaterial::Get()->GetMetalCopperCryo();

  // windows material and thickness
  G4Material* myWindowMat = DSMaterial::Get()->GetAcrylic();
  double myWindowsThickness = 1.5 * cm;
  if (isNoCuVessel && myTPCEdge > 50 * cm) myWindowsThickness = DSStorage::Get()->GetDS20kTPCVesselThickness();

  if (DSStorage::Get()->GetDS20kWindowMaterial() == 1) {  // switch to fused silica
    myWindowMat = DSMaterial::Get()->GetFusedSilica();
    myWindowsThickness = 2.54 / 4 * cm;
    myTopAcrylicThickness = 0 * cm;
  }

  // set some z values for upper part of TPC (required here to set the total
  // height of the liquid)
  double myZPArLayerTop[2] = {myLArGArBoundaryPosZ + myGasPocketThickness, myLArGArBoundaryPosZ + myGasPocketThickness + myPArThickness};
  double myZTPBTop[2] = {myLArGArBoundaryPosZ + myGasPocketThickness + myPArThickness, myLArGArBoundaryPosZ + myGasPocketThickness + myPArThickness + myTPBThickness};
  double myZTopWindow[2] = {myLArGArBoundaryPosZ + myGasPocketThickness + myPArThickness + myTPBThickness, myLArGArBoundaryPosZ + myGasPocketThickness + myPArThickness + myTPBThickness + myWindowsThickness};
  double z1sipm = myZTopWindow[1] + mySiPMOffset;
  double z2sipm = myZTopWindow[1] + mySiPmThickness + mySiPMOffset;
  double myZSiPmTop[2] = {z1sipm, z2sipm};
  double myZSiPmTopBoard[2] = {z2sipm, z2sipm + mySiPmBoardThickness};
  double myZSiPmTopAcrylic[2] = {myZSiPmTopBoard[1], myZSiPmTopBoard[1] + myAcrylicBoardThickness};
  double myZTopAcrylic[2] = {myZSiPmTopAcrylic[1], myZSiPmTopAcrylic[1] + myTopAcrylicThickness};
  double myZBotWindow[2] = {myLArGArBoundaryPosZ - myTPCHeight - myWindowsThickness, myLArGArBoundaryPosZ - myTPCHeight};

  // set the amount of LAr above the TPC
  double myNSLarLevel = myZTopAcrylic[1] + DSStorage::Get()->GetDS20kLArThicknessAboveTPC() + myPMDSuppHeight;

  // --------------------------- //
  // ------   Cryostats   ------ //
  // --------------------------- //

  // Load cryostat profiles
  //  NOTE:  Each cryostat profile is evaluated in the ref frame with z = 0 at
  //  the CENTER of THE TPC
  // G4double myCryostatShiftZ = 0.;
  G4double LArGarIntefaceZ = myLArGArBoundaryPosZ;

  G4double myOuterCryostatZ[300];
  G4double myOuterCryostatRout[300];

  G4double myVacuumCryostatZ[300];
  G4double myVacuumCryostatRout[300];

  G4double myInnerCryostatZ[300];
  G4double myInnerCryostatRout[300];

  G4double myGasArgonZ[300];
  G4double myGasArgonRout[300];

  G4double myTopCryoFillerZ[300];
  G4double myTopCryoFillerRout[300];

  G4double myLiqArgonZ[300];
  G4double myLiqArgonRout[300];

  G4double myRmin[300];
  for (int ii = 0; ii < 300; ++ii) myRmin[ii] = 0.;

  G4int myNumPointsOuterCryo = 0;
  G4int myNumPointsVacCryo = 0;
  G4int myNumPointsInnerCryo = 0;
  // G4int myNumPointsLarGAs_tot = 0;
  G4int myNumPointsGasArgon = 0;
  G4int myNumPointsLiqArgon = 0;
  G4int myNumPointsTopCryoFiller = 0;

  // 2017-March all the hardcoded dimensions of the cryostat are changed with
  // parameters related to the size of the TPC, so that the cryostat scales when
  // DSProto is simulated. 2017-Oct removed reading of cryostat profile,
  // switched to function

  // double myOuterCryostatDistance  = .1958 * myTPCEdge
  // +myGdLayerThickness*1.04 + myCryoToCornerDistance;   // must be larger
  // than 23.5 cm
  double myDefaultR = 1.5667 * myTPCEdge / mm;  // 1880 ; //in mm without the unitmyRingOuterRadius +
                                                // myOuterCryostatDistance ;
  if (myTPCEdge / cm < 50)
    myDefaultR = 1.75 * myTPCEdge / mm;  // 1880 ; //in mm without the unitmyRingOuterRadius   +
                                         // myOuterCryostatDistance ;
  // double mycosalpha = sqrt(2) / 2.;
  // double mysinalpha = sqrt(2) / 2.;
  double myVacuumH = 0.08 * myTPCHeight / mm;      //  190 ; //in mm without the unit
  double myBottomCapH = 0.276 * myTPCHeight / mm;  // 550 + 112; //in mm without the unit
  double myTopCapH = 0.391 * myTPCHeight / mm;     // 700 + 240; //in mm without the unit
  double myDefaultRI = myDefaultR - myVacuumH;     // in mm without the unit //  myRingOuterRadius +
                                                   // myOuterCryostatDistance - 10*cm;

  double myTopOffset = 0.167 / 2. * myTPCHeight / mm + 120;  // 200 ;
  double myBotOffset = 0.125 / 2. * myTPCHeight / mm + 100;  // 150 ;

  PointColPtr pointCol;

  pointCol = DSDetectorProto::createGeometry(myDefaultR, myTPCHeight, 0., myBottomCapH, myTopCapH, 0., myTopOffset, myBotOffset);
  myNumPointsInnerCryo = pointCol->size();
  for (unsigned int i = 0; i < pointCol->size(); ++i) {
    myOuterCryostatZ[i] = pointCol->at(i).first * mm;
    myOuterCryostatRout[i] = pointCol->at(i).second * mm;
  }
  myNumPointsOuterCryo = pointCol->size();

  pointCol = DSDetectorProto::createGeometry(myDefaultR - myCryoWallThickness, myTPCHeight, 0., myBottomCapH - myCryoWallThickness, myTopCapH - myCryoWallThickness, 0., myTopOffset, myBotOffset);
  for (unsigned int i = 0; i < pointCol->size(); ++i) {
    myVacuumCryostatZ[i] = pointCol->at(i).first * mm;
    myVacuumCryostatRout[i] = pointCol->at(i).second * mm;
  }
  myNumPointsVacCryo = pointCol->size();

  pointCol = DSDetectorProto::createGeometry(myDefaultRI, myTPCHeight, 0., myBottomCapH - myVacuumH, myTopCapH - myVacuumH, 0., myTopOffset, myBotOffset);
  for (unsigned int i = 0; i < pointCol->size(); ++i) {
    myInnerCryostatZ[i] = pointCol->at(i).first * mm;
    myInnerCryostatRout[i] = pointCol->at(i).second * mm;
  }
  myNumPointsInnerCryo = pointCol->size();

  pointCol = DSDetectorProto::createGeometry(myDefaultRI - myCryoWallThicknessI, myTPCHeight, 0., myBottomCapH - myVacuumH - myCryoWallThicknessI, myTopCapH - myVacuumH - myCryoWallThicknessI, 0., myTopOffset, myBotOffset);
  int liqn = 0;
  int gasn = 0;
  for (unsigned int i = 0; i < pointCol->size(); ++i) {
    if (pointCol->at(i).first / mm < LArGarIntefaceZ) {
      myLiqArgonZ[i] = pointCol->at(i).first * mm;
      myLiqArgonRout[i] = pointCol->at(i).second * mm;
      ++liqn;
    } else {
      myGasArgonZ[i - liqn + 1] = pointCol->at(i).first * mm;
      myGasArgonRout[i - liqn + 1] = pointCol->at(i).second * mm;
      ++gasn;
    }
    myLiqArgonZ[liqn] = myLiqArgonZ[liqn - 1] - 1 * mm;
    myLiqArgonRout[liqn] = myLiqArgonRout[liqn - 1];
    myLiqArgonZ[liqn + 1] = myTPCHeight / 2.;
    myLiqArgonRout[liqn + 1] = myLiqArgonRout[liqn - 1];
    myLiqArgonZ[liqn + 2] = myTPCHeight / 2.;
    myLiqArgonRout[liqn + 2] = 0.;
    myGasArgonZ[0] = myTPCHeight / 2.;
    myGasArgonRout[0] = 0.;
    myGasArgonZ[1] = myTPCHeight / 2.;
    myGasArgonRout[1] = myLiqArgonRout[liqn];

    myNumPointsLiqArgon = liqn + 3;
    myNumPointsGasArgon = gasn + 1;

    // for (int i=0;i<liqn+3;++i) cout << myLiqArgonZ[i]/cm<<" " <<
    // myLiqArgonRout[i]/cm << endl ; for (int i=0;i<gasn+1;++i) cout <<
    // myGasArgonZ[i]/cm<<" " <<  myGasArgonRout[i]/cm << endl ;
  }
  for (int i = 0; i < gasn + 1; ++i) {
    if (myGasArgonZ[i] < myNSLarLevel) {
      myTopCryoFillerZ[i] = myGasArgonZ[i];
      myTopCryoFillerRout[i] = myGasArgonRout[i];
      ++myNumPointsTopCryoFiller;
    }
  }

  //-------------------------//
  //      Outer VacCryo      //
  //-------------------------//

  ///////////////////////////////////////////////////////////////////////////////
  if (!IsAlternateDesign) {
    ///////////////////////////////////////////////////////////////////////////////

    G4Polycone* fSolidDS20k = new G4Polycone("OuterCryostat_Solid", 0, myTwoPi, myNumPointsOuterCryo, myOuterCryostatZ, myRmin, myOuterCryostatRout);
    G4LogicalVolume* fLogicDS20k = new G4LogicalVolume(fSolidDS20k, myCryoMat, "OuterCryostat_Logic");
    G4double myVerticalShift = 0 * mm;
    // remember to shift a bit the cryostat placement compared to the acrylic
    // vessel if it is there
    if (DSStorage::Get()->GetAcrylicVesselDiameter()) myVerticalShift = DSStorage::Get()->GetDS20kTPCverticalShift();
    fPhysicDS20k = new G4PVPlacement(myDefaultRotation2, G4ThreeVector(0., 0., -myVerticalShift), "OuterCryostat", fLogicDS20k, fMotherVolume, false, 0, myCheckOverlap);

    fLogicDS20k->SetVisAttributes(myGreen);
    // fLogicDS20k->SetVisAttributes(G4VisAttributes::GetInvisible());

    G4Polycone* fSolidDS20kVac = new G4Polycone("OuterCryostat_SolidVac", 0, myTwoPi, myNumPointsVacCryo, myVacuumCryostatZ, myRmin, myVacuumCryostatRout);
    G4LogicalVolume* fLogicDS20kVac = new G4LogicalVolume(fSolidDS20kVac, DSMaterial::Get()->GetVacuum(), "OuterCryostat_LogicVac");
    fPhysicDS20kVac = new G4PVPlacement(0, myZeros, "OuterCryostatVac", fLogicDS20kVac, fPhysicDS20k, false, 0, myCheckOverlap);
    new G4VisAttributes(myWhite);
    fLogicDS20kVac->SetVisAttributes(myYellow);
    // fLogicDS20kVac->SetVisAttributes(G4VisAttributes::GetInvisible());

    //-------------------------//
    //      Inner Cryo         //
    //-------------------------//

    G4Polycone* fSolidInnerCryo = new G4Polycone("InnerCryo_Solid", 0, myTwoPi, myNumPointsInnerCryo, myInnerCryostatZ, myRmin, myInnerCryostatRout);
    G4LogicalVolume* fLogicInnerCryo = new G4LogicalVolume(fSolidInnerCryo, myCryoMat, "SolidInnerCryo_Logic");
    G4PVPlacement* fPhysicInnerCryo = new G4PVPlacement(0, myZeros, "InnerCryo", fLogicInnerCryo, fPhysicDS20kVac, false, 0, myCheckOverlap);
    fLogicInnerCryo->SetVisAttributes(myGreen);
    // fLogicInnerCryo->SetVisAttributes(G4VisAttributes::GetInvisible());

    // evrything must be divided in 2 pieces in Z for LAr and GAr regions
    G4Polycone* fSolidInactiveLar = new G4Polycone("SolidInactiveLar_Solid", 0, myTwoPi, myNumPointsLiqArgon, myLiqArgonZ, myRmin, myLiqArgonRout);
    G4LogicalVolume* fLogicInactiveLar = new G4LogicalVolume(fSolidInactiveLar, DSMaterial::Get()->GetNSLiquidArgon(), "SolidInactiveLar_Logic");
    fPhysicInactiveLar = new G4PVPlacement(0, myZeros, "SolidInactiveLar", fLogicInactiveLar, fPhysicInnerCryo, false, 0, myCheckOverlap);
    new G4VisAttributes(G4Colour(0.0, 0.0, 1.0, 0.9));
    // fLogicInactiveLar->SetVisAttributes(G4VisAttributes::GetInvisible());

    G4Polycone* fSolidInactiveGar = new G4Polycone("SolidInactiveGar_Solid", 0, myTwoPi, myNumPointsGasArgon, myGasArgonZ, myRmin, myGasArgonRout);
    // The GAr fills the upper part of the TPC.
    G4LogicalVolume* fLogicInactiveGar;
    fLogicInactiveGar = new G4LogicalVolume(fSolidInactiveGar, DSMaterial::Get()->GetGaseousArgon(), "SolidInactiveGar_Logic");
    fPhysicInactiveGar = new G4PVPlacement(0, myZeros, "SolidInactiveGar", fLogicInactiveGar, fPhysicInnerCryo, false, 0, myCheckOverlap);
    G4VisAttributes* myGarAttributes = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0, 0.9));
    fLogicInactiveGar->SetVisAttributes(myGarAttributes);

    // NS LiquidArgon Fills the GAr volume up to some level
    G4Polycone* fSolidTopCryoFiller = new G4Polycone("SolidTopCryoFiller_Solid", 0, myTwoPi, myNumPointsTopCryoFiller, myTopCryoFillerZ, myRmin, myTopCryoFillerRout);
    // fill the upper part of the cryostat with liquid argon  up to some
    // specific height (NSLArLevel)
    G4LogicalVolume* fLogicTopCryoFiller = new G4LogicalVolume(fSolidTopCryoFiller, DSMaterial::Get()->GetNSLiquidArgon(), "SolidTopCryoFiller_Logic");
    fPhysicTopCryoFiller = new G4PVPlacement(0, myZeros, "SolidTopCryoFiller", fLogicTopCryoFiller, fPhysicInactiveGar, false, 0, myCheckOverlap);

  } else if (IsAlternateDesign && !isNoCuVessel) {

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Cu Vessel
    G4double myInArBufinnerZ = (myTPCHeight / 2.) + 16 * cm + myPMDSuppHeight + mySiPMOffset;
    G4double myInArBufinnerR = myTPCInnerRadius + 12 * cm;

    double myCuVessZRingBottom[2] = {-myInArBufinnerZ, myInArBufinnerZ};
    double myCuVessInnerRadiusRing[2]{0 * cm, 0 * cm};
    double myCuVessOuterRadiusGdLayer[2] = {myInArBufinnerR, myInArBufinnerR};

    G4Polyhedra* fSolidCuVessel = new G4Polyhedra("fSolidCuVessel_Solid", 0, myTwoPi, 8, 2, myCuVessZRingBottom, myCuVessInnerRadiusRing, myCuVessOuterRadiusGdLayer);
    // Set this to 2 for DSVGenerator
    DSStorage::Get()->SetDS20kCryoMaterial(2);
    G4LogicalVolume* fLogicCuVessel = new G4LogicalVolume(fSolidCuVessel, DSMaterial::Get()->GetMetalCopperCryo(), "CuVess_Logic");
    fLogicCuVessel->SetVisAttributes(myRed);

    G4PVPlacement* fPhysicCuVessel = new G4PVPlacement(myDefaultRotation2, G4ThreeVector(0, 0, 0), "CuVess", fLogicCuVessel, myMotherVolume, false, 0, myCheckOverlap);

    G4double myCuVesselinnerZ = myInArBufinnerZ - 1 * cm;
    G4double myCuVesselinnerR = myInArBufinnerR - 1 * cm;

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // InActiveLAr filling the Cu vessel
    double myInactiveLArZRingBottom[2] = {-myCuVesselinnerZ, myCuVesselinnerZ};
    double myInactiveLArInnerRadiusRing[2]{0 * cm, 0 * cm};
    double myInactiveLArOuterRadiusGdLayer[2] = {myCuVesselinnerR, myCuVesselinnerR};

    G4Polyhedra* fSolidInactiveLArel = new G4Polyhedra("fSolidInactiveLArel_Solid", 0, myTwoPi, 8, 2, myInactiveLArZRingBottom, myInactiveLArInnerRadiusRing, myInactiveLArOuterRadiusGdLayer);

    G4LogicalVolume* fLogicInactiveLArel = new G4LogicalVolume(fSolidInactiveLArel, DSMaterial::Get()->GetNSLiquidArgon(), "InactiveLAr_Logic");
    //  fLogicInactiveLArel->SetVisAttributes(myYellow);

    fPhysicInactiveLar = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), "InactiveLAr", fLogicInactiveLArel,
                                           fPhysicCuVessel,  // fMotherVolume,
                                           false, 0, myCheckOverlap);

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // fPhysicTopCryoFiller is not really needed if the vessel is completely
    // filled with LAr
    fPhysicTopCryoFiller = fPhysicInactiveLar;

    ///////////////////////////////////////////////////////////////////////////////
  } else {  // just the acrylic vessel inside the veto volume, with

    // if (myTPCEdge/m<2.5)fMotherVolume->SetRotation (myDefaultRotation2)  ;
    fPhysicTopCryoFiller = fMotherVolume;
    fPhysicInactiveLar = fMotherVolume;
  }
  ///////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // PMD Support Frame, filled with non scintillating LAr
  double _myPMDSuppOuterRadius = GetOctagonInnerRadius(myPMDSuppEdge);
  double myPMDSuppZ[2] = {-myPMDSuppHeight / 2., myPMDSuppHeight / 2.};
  double myPMDSuppInnerRadius[2]{0 * cm, 0 * cm};
  double myPMDSuppOuterRadius[2] = {_myPMDSuppOuterRadius, _myPMDSuppOuterRadius};

  G4Polyhedra* fSolidPMDSupport = new G4Polyhedra("fSolidPMDSupport_Solid", 0, myTwoPi, 8, 2, myPMDSuppZ, myPMDSuppInnerRadius, myPMDSuppOuterRadius);

  G4LogicalVolume* fLogicPMDSupport = new G4LogicalVolume(fSolidPMDSupport, DSMaterial::Get()->GetSteel(), "PMDSupp_Logic");
  fLogicPMDSupport->SetVisAttributes(myBlue);

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // InLar inside PMD Support
  double PMDempty = 0.97976105;
  double _myInLarPMDSuppEdge = myPMDSuppEdge * PMDempty;
  double _myInLarPMDSuppHeight = myPMDSuppHeight * PMDempty;
  double _myInLarPMDSuppOuterRadius = GetOctagonInnerRadius(_myInLarPMDSuppEdge);

  double myInLarPMDSuppZ[2] = {-_myInLarPMDSuppHeight / 2, _myInLarPMDSuppHeight / 2};
  double myInLarPMDSuppInnerRadius[2]{0 * cm, 0 * cm};
  double myInLarPMDSuppOuterRadius[2] = {_myInLarPMDSuppOuterRadius, _myInLarPMDSuppOuterRadius};

  G4Polyhedra* fSolidInLarPMDSupport = new G4Polyhedra("fSolidInLarPMDSupport_Solid", 0, myTwoPi, 8, 2, myInLarPMDSuppZ, myInLarPMDSuppInnerRadius, myInLarPMDSuppOuterRadius);

  G4LogicalVolume* fLogicInLarPMDSupport = new G4LogicalVolume(fSolidInLarPMDSupport, DSMaterial::Get()->GetNSLiquidArgon(), "InLarPMDSupp_Logic");
  fLogicInLarPMDSupport->SetVisAttributes(myGreen);

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // 2 PMD support Steel Beams, in vertical, inside LAr
  double mySteelBeamPDM_x = 1.95 * GetOctagonInnerRadius(_myInLarPMDSuppEdge);
  double mySteelBeamPDM_y = 5 * mm;
  double mySteelBeamPDM_z = _myInLarPMDSuppHeight;

  G4Box* fSolidSteelBeamPDM = new G4Box("fSolidSteelBeamPDM", mySteelBeamPDM_x / 2, mySteelBeamPDM_y / 2, mySteelBeamPDM_z / 2);

  G4LogicalVolume* fLogicSteelBeamPDM = new G4LogicalVolume(fSolidSteelBeamPDM, DSMaterial::Get()->GetSteel(), "SteelBeamPDM_Logic");
  fLogicPMDSupport->SetVisAttributes(myBlack);

  // G4PVPlacement* fPhysicSteelBeamPDM;

  G4RotationMatrix* rm180 = new G4RotationMatrix();
  rm180->rotateZ(myTwoPi / 16.);
  G4ThreeVector v1(0, -_myInLarPMDSuppEdge / 3.5, 0);
  G4ThreeVector v2(0, +_myInLarPMDSuppEdge / 3.5, 0);
  v1.rotate(-myTwoPi / 16., G4ThreeVector(0, 0, 1));
  v2.rotate(-myTwoPi / 16., G4ThreeVector(0, 0, 1));

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // place Beam 1 inside the PDM support
  new G4PVPlacement(rm180, v1, fLogicSteelBeamPDM, "PMDSuppStBeam", fLogicInLarPMDSupport, false, 0, myCheckOverlap);

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // place Beam 2 inside the PDM support
  new G4PVPlacement(rm180, v2, fLogicSteelBeamPDM, "PMDSuppStBeam", fLogicInLarPMDSupport, false, 0, myCheckOverlap);

  // place InLAr and beams inside the PDM support
  new G4PVPlacement(0, G4ThreeVector(0, 0, 0), fLogicInLarPMDSupport, "PMDSuppInLAr", fLogicPMDSupport, false, 0, myCheckOverlap);

  /// PHYSICAL PLACEMENT of PDM support structure
  // Above the TPC
  new G4PVPlacement(0, G4ThreeVector(0, 0, myZTopAcrylic[1] + myPMDSuppHeight / 2), "PMDSuppTop", fLogicPMDSupport, fPhysicTopCryoFiller, false, 0, myCheckOverlap);

  // Below the TPC
  new G4PVPlacement(0, G4ThreeVector(0, 0, myZBotWindow[0] - 2 * cm - mySiPmBoardThickness - myPMDSuppHeight / 2.), "PMDSuppBot", fLogicPMDSupport, fPhysicInactiveLar, false, 0, myCheckOverlap);

  // 8 teflon pillars at each corner
  G4Tubs* pillar_solid = new G4Tubs("pillar_solid", 0 * cm, 1.5 * cm, abs(myZBotWindow[0]), 0, myTwoPi);
  G4LogicalVolume* pillar_logic = new G4LogicalVolume(pillar_solid, DSMaterial::Get()->GetTeflon(), "pillar_logic");
  G4ThreeVector v3(GetOctagonOuterRadius(myTPCEdge) + 4 * cm, 0, 0);

  // v3.rotate(myDefaultAngle, G4ThreeVector(0,0,1)) ;
  for (int i = 0; i < 8; ++i) {
    v3.rotate(i * myTwoPi / 8., G4ThreeVector(0, 0, 1));
    if (myTPCEdge > 1 * m && !isNoCuVessel) new G4PVPlacement(0, v3, "pillar", pillar_logic, fPhysicInactiveLar, false, 0, i == 0 ? myCheckOverlap : 0);
  }

  //-------------------------//
  //    Shaping Rings        //
  //-------------------------//
  // field shaping rings for now assume a constant wall of 5mm thick outside the
  // PTFE in average for MC first order estimation should be enough for now.
  //(6mm thick 2.2cm high ring spaced 2.5cm and 95 rings total).

  // LArGarIntefaceZ = 1200 *mm ;

  //////////////////////////////////////////////////////////////////////////////////////////////////
  double myZRingBottom[2] = {LArGarIntefaceZ - myTPCHeight + 0 * cm, LArGarIntefaceZ};
  double myInnerRadiusRing[2] = {0., 0.};
  double myOuterRadiusRing[2] = {myRingOuterRadius, myRingOuterRadius};
  // double myOuterRadiusGdLayer[2] = {myRingOuterRadius + myGdLayerThickness,
  // myRingOuterRadius + myGdLayerThickness};

  G4Polyhedra* fSolidRingBottom = new G4Polyhedra("RingBottom_Solid", 0, myTwoPi, 8, 2, myZRingBottom, myInnerRadiusRing, myOuterRadiusRing);

  G4LogicalVolume* fLogicRingBottom = new G4LogicalVolume(fSolidRingBottom, DSMaterial::Get()->GetMetalCopper(), "RingBottom_Logic");
  G4VPhysicalVolume* fPhysicRingBottom;
  fPhysicRingBottom = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), "FieldCageRingBottom", fLogicRingBottom, fPhysicInactiveLar, false, 0, myCheckOverlap);
  G4VisAttributes* myRingBottomAttributes = new G4VisAttributes(myRed);
  fLogicRingBottom->SetVisAttributes(myRingBottomAttributes);

  // G4VisAttributes *myRingTopAttributes = new G4VisAttributes(myBlue);
  // fLogicRingTop->SetVisAttributes(myRingTopAttributes);

  //-------------------------//
  //    Teflon Reflector       //
  //-------------------------//
  // outwards from the inner active wall, a 1-inch thick PTFE reflector
  // (octagonal shape).
  double myZTeflonBottom[2] = {LArGarIntefaceZ - myTPCHeight, LArGarIntefaceZ};
  // double myZTeflonTop[2]             = {LArGarIntefaceZ,LArGarIntefaceZ +
  // myGasPocketThickness};
  double myInnerRadiusTeflon[2] = {0, 0};
  double myOuterRadiusTeflon[2] = {myTeflonInnerRadius, myTeflonInnerRadius};

  // reflector bottom
  G4Polyhedra* fSolidTeflonBottom = new G4Polyhedra("TeflonBottom_Solid", 0, myTwoPi, 8, 2, myZTeflonBottom, myInnerRadiusTeflon, myOuterRadiusTeflon);

  G4LogicalVolume* fLogicTeflonBottom = new G4LogicalVolume(fSolidTeflonBottom, IsAcrylicReflector ? DSMaterial::Get()->GetAcrylic() : DSMaterial::Get()->GetTeflon(), "TeflonBottom_Logic");

  fPhysicTeflonBottom = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), "TeflonBottom", fLogicTeflonBottom, fPhysicRingBottom, false, 0, myCheckOverlap);
  //-------------------------//
  //          TPC            //
  //-------------------------//
  // TPC active volume section: Octagonal shape with edge 64.2cm when cold. (or
  // Inscribed circle radius = 77.5cm)
  double myZTPCBottom[2] = {LArGarIntefaceZ - myTPCHeight + myTPBThickness, LArGarIntefaceZ};
  double myZTPBSide[2] = {LArGarIntefaceZ - myTPCHeight, LArGarIntefaceZ};
  double myInnerRadiusTPC[2] = {0, 0};
  // double myOuterRadiusSiPMring[2]      = {myTPCInnerRadius,myTPCInnerRadius};
  // double myInnerRadiusSiPMring[2]      =
  // {myTPCInnerRadius-5*mm,myTPCInnerRadius-5*mm};
  double myOuterRadiusTPC[2] = {myTPCInnerRadius, myTPCInnerRadius};
  double myOuterRadiusTPBSide[2] = {myTPCInnerRadius + myTPBThickness, myTPCInnerRadius + myTPBThickness};

  // TPB layer (side and bottom)
  G4Polyhedra* fSolidTPBSide = new G4Polyhedra("TPBSide_Solid", 0, myTwoPi, 8, 2, myZTPBSide, myInnerRadiusTPC, myOuterRadiusTPBSide);

  G4LogicalVolume* fLogicTPBSide = new G4LogicalVolume(fSolidTPBSide, DSMaterial::Get()->GetTPB(), "TPBSide_Logic");
  fPhysicTPBSide = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), "TPBSide", fLogicTPBSide, fPhysicTeflonBottom, false, 0, myCheckOverlap);
  //////////////////////////////////////////////////////////////////////////////////////////////////
  // TPC Active Volume
  G4Polyhedra* fSolidActiveLAr = new G4Polyhedra("ActiveLAr_Solid", 0, myTwoPi, 8, 2, myZTPCBottom, myInnerRadiusTPC, myOuterRadiusTPC);

  G4LogicalVolume* fLogicActiveLAr = new G4LogicalVolume(fSolidActiveLAr, myTPCEdge / m < 2.5 ? DSMaterial::Get()->GetLiquidArgon() : DSMaterial::Get()->GetXeLiquidArgon(), "LAr_Logic");
  fPhysicActiveLAr = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), "ActiveLAr", fLogicActiveLAr,
                                       fPhysicTPBSide,  // fPhysicTPBSide,//fPhysicRingBottom,
                                       false, 0, myCheckOverlap);

  fLogicActiveLAr->SetVisAttributes(myBlue);
  G4Region* fLArRegion = new G4Region("LAr_Logic");
  fLogicActiveLAr->SetRegion(fLArRegion);
  fLArRegion->AddRootLogicalVolume(fLogicActiveLAr);

  DSStorage::Get()->SetLiquidArgonIndex((int)fLogicActiveLAr->GetMaterial()->GetIndex());
  //-------------------------------------------//
  //  Active GAr Region      //
  //-------------------------------------------//

  double myZTPCTop[2] = {LArGarIntefaceZ, LArGarIntefaceZ + myActiveGasLogicThickness};

  G4Polyhedra* fSolidGasPocket = new G4Polyhedra("GasPocket_Solid", 0, myTwoPi, 8, 2, myZTPCTop, myInnerRadiusTPC, myOuterRadiusTPC);

  G4LogicalVolume* fLogicGasPocket = new G4LogicalVolume(fSolidGasPocket, DSMaterial::Get()->GetGaseousArgon(), "GasPocket_Logic");

  fPhysicGasPocket = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), "GasPocket", fLogicGasPocket,
                                       fPhysicTopCryoFiller,  // fPhysicRingTop,
                                       false, 0, myCheckOverlap);
  //-------------------------------//
  //          Top and Bottom caps  //
  //-------------------------------//

  // Add offset of 5 mm
  z1sipm = myZBotWindow[0] - mySiPmThickness - 5 * mm;
  z2sipm = myZBotWindow[0] - 5 * mm;
  double myZSiPmBottom[2] = {z1sipm, z2sipm};
  double myZSiPmBottomBoard[2] = {z1sipm - mySiPmBoardThickness, z1sipm};
  double myZSiPmBottomAcrylic[2] = {myZSiPmBottomBoard[0] - myAcrylicBoardThickness, myZSiPmBottomBoard[0]};

  // top window
  G4Polyhedra* fSolidTopWindow = new G4Polyhedra("TopWindow_Solid", 0, myTwoPi, 8, 2, myZTopWindow, myInnerRadiusTPC, myOuterRadiusTPC);

  G4LogicalVolume* fLogicTopWindow = new G4LogicalVolume(fSolidTopWindow, myWindowMat, "TopWindow_Logic");
  fPhysicTopWindow = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), "TopWindow", fLogicTopWindow,
                                       fPhysicTopCryoFiller,  // fPhysicGasPocket,
                                       false, 0, myCheckOverlap);

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // TPB in between gar and window (top)
  G4Polyhedra* fSolidTPBTop = new G4Polyhedra("TPBTop_Solid", 0, myTwoPi, 8, 2, myZTPBTop, myInnerRadiusTPC, myOuterRadiusTPC);

  G4LogicalVolume* fLogicTPBTop = new G4LogicalVolume(fSolidTPBTop, DSMaterial::Get()->GetTPB(), "TPBTop_Logic");

  fPhysicTPBTop = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), "TPBTop", fLogicTPBTop,
                                    fPhysicTopCryoFiller,  // fPhysicGasPocket,
                                    false, 0, myCheckOverlap);
  if (myPArThickness > 0.00) {
    //////////////////////////////////////////////////////////////////////////////////////////////////
    // Pseudo Ar layer
    G4Polyhedra* fSolidLArLayer = new G4Polyhedra("LArLayer_Solid", 0, myTwoPi, 8, 2, myZPArLayerTop, myInnerRadiusTPC, myOuterRadiusTPC);

    G4LogicalVolume* fLogicLArLayer = new G4LogicalVolume(fSolidLArLayer, DSMaterial::Get()->GetPseudoArgon(), "LArLayer_Logic");
    // remove pseudo argon (WRONG NAME?)
    fPhysicTPBTop = new G4PVPlacement(0, G4ThreeVector(0, 0, 0.), "LArLayer", fLogicLArLayer, fPhysicTopCryoFiller, false, 0, myCheckOverlap);
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // bottom window
  G4Polyhedra* fSolidBotWindow = new G4Polyhedra("BotWindow_Solid", 0, myTwoPi, 8, 2, myZBotWindow, myInnerRadiusTPC, myOuterRadiusTeflon);

  G4LogicalVolume* fLogicBotWindow = new G4LogicalVolume(fSolidBotWindow, myWindowMat, "BotWindow_Logic");

  fPhysicBotWindow = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), "BotWindow", fLogicBotWindow,
                                       fPhysicInactiveLar,  // fPhysicRingTop,
                                       false, 0, myCheckOverlap);
  fLogicBotWindow->SetVisAttributes(G4Color(1.0, 1.0, 0.0, 0.6));
  fLogicTopWindow->SetVisAttributes(G4Color(1.0, 1.0, 0.0, 0.6));
  fLogicTeflonBottom->SetVisAttributes(G4Color(1.0, 1.0, 0.0, 0.6));

  /////////////////////////////////////////////////////////////////////
  // Top TPC belt/flange
  if (isNoCuVessel && myTPCEdge > 50 * cm && false) {
    ////////////////////////////////////////////////////////////////////
    double myZBelt[2] = {myZSiPmTop[0] - myReflectorThickness, myZSiPmTop[0]};
    double myOuterRadiusBelt[2] = {myOuterRadiusTPC[0] + myReflectorThickness + 5 * cm, myOuterRadiusTPC[1] + myReflectorThickness + 5 * cm};

    G4Polyhedra* fSolidBelt = new G4Polyhedra("Belt_Solid", 0, myTwoPi, 8, 2, myZBelt, myOuterRadiusTPC, myOuterRadiusBelt);

    G4LogicalVolume* fLogicBelt = new G4LogicalVolume(fSolidBelt, DSMaterial::Get()->GetAcrylic(), "Belt_Logic");

    new G4PVPlacement(0, G4ThreeVector(0, 0, 0), "Belt", fLogicBelt, fPhysicInactiveLar, false, 0, myCheckOverlap);

    fLogicBelt->SetVisAttributes(G4Color(1.0, 1.0, 0.0, 0.6));

    ////////////////////////////////////////////////////////////////////
  }

  //-----------------------------------//
  //      SiPM arrays and support      //
  //-----------------------------------//

  // Top Array
  G4Polyhedra* fSolidSiPMTop = new G4Polyhedra("SiPMTop_Solid", 0, myTwoPi, 8, 2, myZSiPmTop, myInnerRadiusTPC, myOuterRadiusTPC);

  G4LogicalVolume* fLogicSiPMTop = new G4LogicalVolume(fSolidSiPMTop, DSMaterial::Get()->GetMetalSilicon(), "SiPMTop_Logic");

  fPhysicSiPmTop = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), "SiPMTop", fLogicSiPMTop,
                                     fPhysicTopCryoFiller,  // fPhysicGasPocket,//fPhysicRingBottom,
                                     false, -56, myCheckOverlap);

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Top Array Board
  G4Polyhedra* fSolidSiPMTopBoard = new G4Polyhedra("SiPMTop_SolidBoard", 0, myTwoPi, 8, 2, myZSiPmTopBoard, myInnerRadiusTPC, myOuterRadiusTPC);

  G4LogicalVolume* fLogicSiPMTopBoard = new G4LogicalVolume(fSolidSiPMTopBoard, DSMaterial::Get()->GetMetalCopper(), "SiPMTop_LogicBoard");

  // G4PVPlacement* fPhysicSiPmTopBoard =
  new G4PVPlacement(0, G4ThreeVector(0, 0, 0), "SiPMTopBoard", fLogicSiPMTopBoard,
                    fPhysicTopCryoFiller,  // fPhysicGasPocket,//fPhysicRingBottom,
                    false, 0, myCheckOverlap);

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Top Array Acrylic layer to account for PDM acrylic
  G4Polyhedra* fSolidSiPMTopAcrylicBoard = new G4Polyhedra("SiPMTop_SolidBoard", 0, myTwoPi, 8, 2, myZSiPmTopAcrylic, myInnerRadiusTPC, myOuterRadiusTPC);

  G4LogicalVolume* fLogicSiPMTopAcrylicBoard = new G4LogicalVolume(fSolidSiPMTopAcrylicBoard, DSMaterial::Get()->GetAcrylic(), "SiPMTop_AcrylicBoard");

  // G4PVPlacement* fPhysicSiPmTopAcrylicBoard =
  new G4PVPlacement(0, G4ThreeVector(0, 0, 0), "SiPMTopAcrylicBoard", fLogicSiPMTopAcrylicBoard, fPhysicTopCryoFiller, false, 0, myCheckOverlap);

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Top Array Acrylic
  if (myZTopAcrylic[1] != myZTopAcrylic[0] && myWindowsThickness > 1 * cm) {
    G4Polyhedra* fSolidSiPMTopTopAcrylic = new G4Polyhedra("SiPMTop_SolidTopAcrylic", 0, myTwoPi, 8, 2, myZTopAcrylic, myInnerRadiusTPC, myOuterRadiusTPC);

    G4LogicalVolume* fLogicSiPMTopTopAcrylic = new G4LogicalVolume(fSolidSiPMTopTopAcrylic, DSMaterial::Get()->GetAcrylic(), "SiPMTop_LogicTopAcrylic");

    // G4PVPlacement* fPhysicSiPmTopTopAcrylic =
    new G4PVPlacement(0, G4ThreeVector(0, 0, 0), "SiPMTopTopAcrylic", fLogicSiPMTopTopAcrylic,
                      fPhysicTopCryoFiller,  // fPhysicGasPocket,//fPhysicRingBottom,
                      false, 0, myCheckOverlap);
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Bottom Array
  G4Polyhedra* fSolidSiPMBottom = new G4Polyhedra("SiPMBottom_Solid", 0, myTwoPi, 8, 2, myZSiPmBottom, myInnerRadiusTPC, myOuterRadiusTPC);

  G4LogicalVolume* fLogicSiPMBottom = new G4LogicalVolume(fSolidSiPMBottom, DSMaterial::Get()->GetMetalSilicon(), "SiPMBottom_Logic");

  fPhysicSiPmBottom = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), "SiPMBottom", fLogicSiPMBottom,
                                        fPhysicInactiveLar,  // fPhysicRingBottom,
                                        false, -56, myCheckOverlap);

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Bottom Array  Board
  G4Polyhedra* fSolidSiPMBottomBoard = new G4Polyhedra("SiPMBottomBoard_Solid", 0, myTwoPi, 8, 2, myZSiPmBottomBoard, myInnerRadiusTPC, myOuterRadiusTPC);

  G4LogicalVolume* fLogicSiPMBottomBoard = new G4LogicalVolume(fSolidSiPMBottomBoard, DSMaterial::Get()->GetMetalCopper(), "SiPMBottomBoard_Logic");

  // G4PVPlacement* fPhysicSiPmBottomBoard =
  new G4PVPlacement(0, G4ThreeVector(0, 0, 0), "SiPMBottomBoard", fLogicSiPMBottomBoard,
                    fPhysicInactiveLar,  // fPhysicRingBottomBoard,
                    false, 0, myCheckOverlap);

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Bottom Array Acrylic Board to account for PDM acrylic

  G4Polyhedra* fSolidSiPMBottomAcrylicBoard = new G4Polyhedra("SiPMBottomAcrylicBoard_Solid", 0, myTwoPi, 8, 2, myZSiPmBottomAcrylic, myInnerRadiusTPC, myOuterRadiusTPC);

  G4LogicalVolume* fLogicSiPMBottomAcrylicBoard = new G4LogicalVolume(fSolidSiPMBottomAcrylicBoard, DSMaterial::Get()->GetAcrylic(), "SiPMBottomAcrylicBoard_Logic");

  // G4PVPlacement* fPhysicSiPMBottomAcrylicBoard =
  new G4PVPlacement(0, G4ThreeVector(0, 0, 0), "SiPMBottomAcrylicBoard", fLogicSiPMBottomAcrylicBoard,
                    fPhysicInactiveLar,  // fPhysicRingBottomBoard,
                    false, 0, myCheckOverlap);

  ///////////////////////////////////////////////////////////////////////////////////////
  // Grid
  ///////////////////////////////////////////////////////////////////////////////////////

  G4double myGrid_h = 0.01 * cm;
  double myZGrid[2] = {LArGarIntefaceZ - .5 * cm, LArGarIntefaceZ - .5 * cm + myGrid_h};
  G4ThreeVector myGridPosition(0, 0, myLArGArBoundaryPosZ - 0.5 * cm);

  fSolidGrid = new G4Polyhedra("Grid_Solid", 0, myTwoPi, 8, 2, myZGrid, myInnerRadiusTPC, myOuterRadiusTPC);

  fLogicGrid = new G4LogicalVolume(fSolidGrid, DSMaterial::Get()->GetGridSteel(), "Grid_Logic");
  fPhysicGrid = new G4PVPlacement(0, myZeros, "Grid", fLogicGrid, fPhysicActiveLAr, false, 0, myCheckOverlap);

  DefineSurfaces();

  // make SiPM as pe storing material
  DSStorage::Get()->SetPMTMaterialIndex(fLogicSiPMBottom->GetMaterial()->GetIndex());
}

DSDetectorProto::~DSDetectorProto() {
  ;  // delete fMessenger;
}

void DSDetectorProto::DefineSurfaces() {
  ////////////////////////////////////////
  // Copied from DS DetectorDS50 - 24th april 2015
  ////////////////////////////////////////

  // TPC - BScint
  // warning: this only works if the current 20k design with LSV is constructed
  G4OpticalSurface* fOpElectropolishedStainlessSteelSurface = new G4OpticalSurface("OpElectropolishedStainlessSteelSurface");
  fOpElectropolishedStainlessSteelSurface->SetType(dielectric_metal);
  fOpElectropolishedStainlessSteelSurface->SetModel(unified);
  fOpElectropolishedStainlessSteelSurface->SetFinish(groundbackpainted);
  fOpElectropolishedStainlessSteelSurface->SetMaterialPropertiesTable(DSMaterial::Get()->GetElectropolishedStainlessSteelMPT());
  fDS20kOuterSurface = new G4LogicalBorderSurface("DS20kOuterSurface", fMotherVolume, fPhysicDS20k, fOpElectropolishedStainlessSteelSurface);

  // LAR-Grid
  G4OpticalSurface* fOpGridLArSurface = new G4OpticalSurface("OpGridLArSurface");
  new G4LogicalBorderSurface("GridLArSurface", fPhysicActiveLAr, fPhysicGrid, fOpGridLArSurface);
  fOpGridLArSurface->SetType(dielectric_dielectric);
  fOpGridLArSurface->SetModel(glisur);
  fOpGridLArSurface->SetFinish(polished);
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
  fOpGridLArSurface->SetMaterialPropertiesTable(fGridLArSurfProp);

  // Grid->LAR
  if (DSParameters::Get()->GetWithNewGridModel()) {
    G4OpticalSurface* fOpLArGridSurface = new G4OpticalSurface("OpLArGridSurface");
    new G4LogicalBorderSurface("LArGridSurface", fPhysicGrid, fPhysicActiveLAr, fOpLArGridSurface);
    cout << " With DS50 new grid model " << endl;
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
  G4OpticalSurface* fOpTPBGArSurface = new G4OpticalSurface("OpTPBGArSurface");
  new G4LogicalBorderSurface("GArTPBSurface", fPhysicGasPocket, fPhysicTPBTop, fOpTPBGArSurface);
  new G4LogicalBorderSurface("TPBGArSurface", fPhysicTPBTop, fPhysicGasPocket, fOpTPBGArSurface);
  new G4LogicalBorderSurface("LArTPBSurfaceSide", fPhysicActiveLAr, fPhysicTPBSide, fOpTPBGArSurface);
  new G4LogicalBorderSurface("TPBLArSurfaceSide", fPhysicTPBSide, fPhysicActiveLAr, fOpTPBGArSurface);
  //  new G4LogicalBorderSurface("LArTPBSurfaceBot", fPhysicActiveLAr,
  //  fPhysicTPBBottom, fOpTPBGArSurface ); new
  //  G4LogicalBorderSurface("TPBLArSurfaceBot", fPhysicTPBBottom,
  //  fPhysicActiveLAr, fOpTPBGArSurface );
  new G4LogicalBorderSurface("LArLayerTPBSurface", fPhysicLArLayer, fPhysicTPBTop, fOpTPBGArSurface);
  new G4LogicalBorderSurface("TPBLArLayerSurface", fPhysicTPBTop, fPhysicLArLayer, fOpTPBGArSurface);
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
  if (DSParameters::Get()->GetWithITO()) fITOSurfProp->AddConstProperty("DOITO", 1);

  ////////////////////////////////////////
  // BellTop (acrylic) <--> TPB and CathodeWindow <--> TPB (both with ITO)
  // In the current model, the diffuse nature of the TPB is handled
  // entirely at the TPB-GAr/LAr surface, not here.
  // Make this bi-directional.
  ////////////////////////////////////////
  G4OpticalSurface* fOpWindowTPBSurface = new G4OpticalSurface("OpWindowTPBSurface");
  new G4LogicalBorderSurface("TopWindowTPBSurface", fPhysicTopWindow, fPhysicTPBTop, fOpWindowTPBSurface);
  new G4LogicalBorderSurface("BottomWindowTPBSurface", fPhysicBotWindow, fPhysicTPBSide, fOpWindowTPBSurface);
  new G4LogicalBorderSurface("TPBTopWindowSurface", fPhysicTPBTop, fPhysicTopWindow, fOpWindowTPBSurface);
  new G4LogicalBorderSurface("TPBBottomWindowSurface", fPhysicTPBSide, fPhysicBotWindow, fOpWindowTPBSurface);
  fOpWindowTPBSurface->SetType(dielectric_dielectric);
  fOpWindowTPBSurface->SetModel(unified);
  //  fOpWindowTPBSurface->SetFinish( ground );
  fOpWindowTPBSurface->SetFinish(polished);
  // G4MaterialPropertiesTable *fWindowTPBSurfProp = new
  // G4MaterialPropertiesTable();
  fOpWindowTPBSurface->SetMaterialPropertiesTable(fITOSurfProp);

  // TODO 2016_12: one error to be checked with this surface
  /*
  ////////////////////////////////////////
  // BellTop (acrylic) <--> LAr (no ITO)         2017-02-16
  // Make this bi-directional.
  ////////////////////////////////////////
  G4OpticalSurface *fOpTopWindowArSurface     = new
  G4OpticalSurface("OpTopWindowArSurface"); new
  G4LogicalBorderSurface("TOpTopWindowLArSurfaceIn",   fPhysicTopCryoFiller,
  fPhysicTopWindow , fOpTopWindowArSurface ); new
  G4LogicalBorderSurface("TOpTopWindowLArSurfaceOut",   fPhysicTopWindow ,
  fPhysicTopCryoFiller, fOpTopWindowArSurface ); fOpTopWindowArSurface->SetType(
  dielectric_dielectric ); fOpTopWindowArSurface->SetModel( unified );
  //  fOpTopWindowArSurface->SetFinish( ground );
  fOpTopWindowArSurface->SetFinish( polished );
  //G4MaterialPropertiesTable *fWindowArSurfProp = new
  G4MaterialPropertiesTable();
  // fOpWindowArSurface->SetMaterialPropertiesTable( fITOSurfProp );     // NO
  ITO HERE
  */

  ////////////////////////////////////////
  // LAr and CathodeWindow <--> LAr (with ITO)    2017-02-16
  // Make this bi-directional.
  ////////////////////////////////////////
  G4OpticalSurface* fOpBotWindowArSurface = new G4OpticalSurface("OpBotWindowArSurface");
  new G4LogicalBorderSurface("BottomWindowLArSurfaceIn", fPhysicInactiveLar, fPhysicBotWindow, fOpBotWindowArSurface);
  new G4LogicalBorderSurface("BottomWindowLArSurfaceOut", fPhysicBotWindow, fPhysicInactiveLar, fOpBotWindowArSurface);
  fOpBotWindowArSurface->SetType(dielectric_dielectric);
  fOpBotWindowArSurface->SetModel(unified);
  //  fOpBotWindowArSurface->SetFinish( ground );
  fOpBotWindowArSurface->SetFinish(polished);
  // G4MaterialPropertiesTable *fWindowArSurfProp = new
  // G4MaterialPropertiesTable();
  fOpBotWindowArSurface->SetMaterialPropertiesTable(fITOSurfProp);  // With ITO HERE

  /*

  // not necessary
  ////////////////////////////////////////
  // BellTop (acrylic) <--> SiPM and CathodeWindow <--> SiPM (both with ITO)
  // Make this bi-directional.
  ////////////////////////////////////////
  G4OpticalSurface *fOpWindowSiPmSurface     = new
  G4OpticalSurface("OpWindowSiPmSurface"); new
  G4LogicalBorderSurface("TopWindowSiPmSurface",   fPhysicTopWindow ,
  fPhysicSiPmTop, fOpWindowSiPmSurface ); new
  G4LogicalBorderSurface("BottomWindowSiPmSurface", fPhysicBotWindow,
  fPhysicSiPmBottom, fOpWindowSiPmSurface );
  //new G4LogicalBorderSurface("SiPmTopWindowSurface",
  fPhysicSiPmTop,fPhysicTopWindow ,  fOpWindowSiPmSurface );
  //new G4LogicalBorderSurface("SiPmBottomWindowSurface", fPhysicSiPmBottom,
  fPhysicBotWindow, fOpWindowSiPmSurface ); fOpWindowSiPmSurface->SetType(
  dielectric_dielectric ); fOpWindowSiPmSurface->SetModel( unified );
  //  fOpWindowSiPmSurface->SetFinish( ground );
  fOpWindowSiPmSurface->SetFinish( polished );
  //G4MaterialPropertiesTable *fWindowSiPmSurfProp = new
  G4MaterialPropertiesTable(); fOpWindowSiPmSurface->SetMaterialPropertiesTable(
  fITOSurfProp );
  */

  ///////////////
  // TPB --> Teflon (Reflector)
  //  Should be no Teflon --> TPB as surface is defined as dielectric_metal.
  ////////////////////////////////////////
  G4double TeflonTPBENE[4] = {0.1 * eV, 8.0 * eV, 8.3 * eV, 20.0 * eV};
  G4double TREFUV = DSParameters::Get()->GetTeflonTPBUVRef();
  G4double TREFVIS = DSParameters::Get()->GetTeflonTPBVisRef();
  G4double TeflonTPBREF[4] = {TREFVIS, TREFVIS, TREFUV, TREFUV};
  G4OpticalSurface* fOpTPBTeflonSurface = new G4OpticalSurface("OpTBPTeflonSurface");
  new G4LogicalBorderSurface("TPBTeflonSurface", fPhysicTPBSide, fPhysicTeflonBottom, fOpTPBTeflonSurface);
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

  /*


  // TODO 2016_12: SiPM reflectivity (so far everything is absorbed, QE applied
  offline)

  ///////////////
  // LAr --> SiPM
  ////////////////////////////////////////
  G4double TeflonTPBENE2[4] = {0.1*eV, 8.0*eV, 8.3*eV, 20.0*eV};
  TREFUV  = 0; //DSParameters::Get()->GetTeflonTPBUVRef();
  TREFVIS = 0; //DSParameters::Get()->GetTeflonTPBVisRef();
  G4double TeflonTPBREF2[4] = {TREFVIS, TREFVIS ,TREFUV , TREFUV };
  G4OpticalSurface *fOpLArSiPMSurface = new
  G4OpticalSurface("OpLArSiPMSurface"); new
  G4LogicalBorderSurface("LArSiPMSurface", fPhysicInactiveGar, fPhysicSiPmTop,
  fOpLArSiPMSurface ); new G4LogicalBorderSurface("LArSiPMSurface",
  fPhysicInactiveLar, fPhysicSiPmBottom, fOpLArSiPMSurface );
  fOpLArSiPMSurface->SetType( dielectric_metal );
  fOpLArSiPMSurface->SetModel(unified);

  // PDM: though I can't see how, the following settings are giving the desired
  Lambertian reflection fOpLArSiPMSurface->SetFinish(groundfrontpainted);
  //fOpLArSiPMSurface->SetFinish(ground);
  fOpLArSiPMSurface->SetSigmaAlpha(0.1);

  G4MaterialPropertiesTable *fTPBTeflonSurfProp2 = new
  G4MaterialPropertiesTable(); fTPBTeflonSurfProp2->AddProperty("REFLECTIVITY",
  TeflonTPBENE2, TeflonTPBREF2, 4); G4double TeflonTPBREF3[4] = {1, 1 ,1 , 1 };
  fTPBTeflonSurfProp2->AddProperty("TRANSMITTANCE", TeflonTPBENE2,
  TeflonTPBREF3, 4); fOpLArSiPMSurface->SetMaterialPropertiesTable(
  fTPBTeflonSurfProp2 );
 */
  new G4LogicalBorderSurface("SSteelOuterSurface", fMotherVolume, fPhysicDS20kVac, fOpElectropolishedStainlessSteelSurface);
}

double DSDetectorProto::GetOctagonInnerRadius(double edge) {
  return edge / 2 * (1 + sqrt(2));
}
double DSDetectorProto::GetOctagonOuterRadius(double edge) {
  return edge / 2. * sqrt(4 + 2 * sqrt(2));
}

PointColPtr DSDetectorProto::createGeometry(double r0, double hc, double z0, double hb, double ht, double offset = 0, double mytopoff = 200, double mybotoff = 150) {
  // creates a cylindric shape with end caps in the form of rotation ellipsoids
  //
  // input:
  // r0: radius of cylinder
  // hc: height of cylinder
  // z0: symmetry plane of cylinder
  // hb: height of bottom cap
  // ht: height of top cap
  // offset: offset to the otherwise defined geometry

  // half height of cylinder
  const float hc2 = 0.5 * hc;

  const int nb = 35;
  const int nc = 36;
  const int nt = 35;

  PointColPtr points = new PointCol();

  // bottom cap

  // normal in -z direction, offset in direction -z

  // float mybotoff = 150 ;
  points->push_back(PairFF(z0 - hc2 - hb - offset - mybotoff, 0));
  for (int i = 1; i < nb; ++i) {
    const float angle = float(i) * 90.0 / float(nb);
    const float rads = M_PI / 180.0 * angle;
    float z = z0 - hc2 - hb * std::cos(rads) - mybotoff;
    float r = r0 * std::sin(rads);

    if (offset > 0) {
      // normal vector
      float n_z = -r0 * std::cos(rads);
      float n_r = hb * std::sin(rads);
      const float l = std::sqrt(n_z * n_z + n_r * n_r);
      // normal vector with length "offset"
      n_z *= offset / l;
      n_r *= offset / l;
      z += n_z;
      r += n_r;
    }
    points->push_back(PairFF(z, r));
  }

  // cylinder
  for (int i = 0; i < nc; ++i) {
    // normal in radial direction, offset in radius
    const float z = z0 + hc * (float(i) / float(nc - 1) - 0.5);
    const float r = r0 + offset;
    points->push_back(PairFF(z, r));
  }

  // top cap
  // float mytopoff = 200 ;
  for (int i = 0; i < mytopoff / 10; ++i) {
    float z = z0 + hc2 + 5 + i * 10;
    float r = r0 + offset;
    points->push_back(PairFF(z, r));
  }
  for (int i = 1; i < nt; ++i) {
    const float angle = float(i) * 90.0 / float(nb);
    const float rads = M_PI / 180.0 * angle;
    float z = z0 + hc2 + ht * std::sin(rads) + mytopoff;
    float r = r0 * std::cos(rads);
    if (offset > 0) {
      // normal vector
      float n_z = r0 * std::sin(rads);
      float n_r = ht * std::cos(rads);
      const float l = std::sqrt(n_z * n_z + n_r * n_r);
      // normal vector with length "offset"
      n_z *= offset / l;
      n_r *= offset / l;
      z += n_z;
      r += n_r;
    }
    points->push_back(PairFF(z, r));
  }
  // normal in +z direction, offset in direction +z
  points->push_back(PairFF(z0 + hc2 + ht + offset + mytopoff, 0));

  return points;
}
