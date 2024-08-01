#include <iostream>
#include "G4Box.hh"
#include "G4Orb.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalConstants.hh"
#include "G4Polycone.hh"
#include "G4SubtractionSolid.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4UserLimits.hh"


#include "DSDetectorDS20k.hh"
#include "DSIO.hh"
#include "DSLogger.hh"
#include "DSMaterial.hh"
#include "DSParameters.hh"
#include "DSStorage.hh"
#include "G4Colour.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4VisAttributes.hh"

using namespace std;

////////////////////////////////////////////////
//////        Detector Description      ////////
////////////////////////////////////////////////
// 2022 PA Merge to master branch and cleaning
// 2021 MR Implementaiton of planC geoemtry and veto SiPM placement
// May 2019
// schematics of the 20k detector
// OPTICS is directly derived from DS50
//
// Placements are realized with no shifts. The correct position is calculated
// when making the solids.

DSDetectorDS20k::DSDetectorDS20k(G4VPhysicalVolume* myMotherVolume) {

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

  vector<G4ThreeVector> SiPMPosVector;

  G4bool myCheckOverlap = DSStorage::Get()->GetCheckOverlap();

  G4ThreeVector myZeros(0., 0., 0.);

  // Parameters from drawings
  G4double myTPCHeight = DSStorage::Get()->GetDS20kTPCheight();
  G4double myTPCEdge = DSStorage::Get()->GetDS20kTPCedge();

  // Reflector and windows
  double myWindowsThickness = 5 * cm;                                            // default = 5 cm
  double myReflectorThickness = DSStorage::Get()->GetDS20kTPCVesselThickness();  // default 5 cm

  double myTPCInnerRadius = GetOctagonInnerRadius(myTPCEdge);
  double myReflectorInnerRadius = myTPCInnerRadius + myReflectorThickness;

  // place the center of the TPC at the center of the ref system
  double myLArGArBoundaryPosZ = DSStorage::Get()->GetDS20kTPCheight() / 2.;
  // Set the z coordinate of the LAr - GAr interface, necessary for S2
  // generation in DSLightX
  DSStorage::Get()->SetLArGArBoundaryPosZ(myLArGArBoundaryPosZ + 0.1 * um);

  // mass of material above and below the TPC
  //G4double mass_PCB_per_tile = 20 * g;    // Arlon
  //G4double mass_Invar_per_tile = 10 * g;  // Invar
  G4double mass_steel_plates = 810 * kg;  // Mass of each OP Steel structure

  bool remove_bonding_ring = DSStorage::Get()->GetDS20kRemoveBondingRing();

  //Step Limiter (3/3) - reduce step size for all particles in reflector volumes
  G4double maxStep = 0.1*mm;
  G4UserLimits * myStepLimit = new G4UserLimits(maxStep);

  // Other parameters
  double myTPBThickness = DSStorage::Get()->GetDS20kTPBThickness();
  double myPArThickness = 0. * mm;
  double mySiPMOffset_top = DSStorage::Get()->GetSiPMOffset();               // default: 5 cm
  double mySiPMOffset_bottom = DSStorage::Get()->GetSiPMOffset() + 1.5 * cm;  // default: 3 cm
  double mySiPMOffset_mismatch = mySiPMOffset_top - mySiPMOffset_bottom;
  double mySiPmThickness = 0.7 * mm;
  double myGasPocketThickness = DSStorage::Get()->GetDS20kGasPocketThickness();
  double myActiveGasLogicThickness = myGasPocketThickness;
  double myElectronicsSpace = DSStorage::Get()->GetDS20kLArThicknessAboveTPC();  // space above (below) top (bottom)
                                                                                 // window to accomodate electronics
                                                                                 // and mechanics
  double myElectronicsSpaceTop = myElectronicsSpace - 5.79 * cm;
  double myElectronicsSpaceBottom = myElectronicsSpace - 4.29 * cm;

  double myGdThicknessShield = 15 * cm;

  G4double LArGarIntefaceZ = myLArGArBoundaryPosZ;

  const G4double myPENThickness = DSStorage::Get()->GetDS20kWLSPENThickness();
  //const G4double myNSLArThickness = DSStorage::Get()->GetDS20kWLSLArThickness();
  const G4double myNSLArThickness = 1. * cm;
  
  // Define extension of the TPC vessel
  double top_vessel_4 = LArGarIntefaceZ + myGasPocketThickness + myWindowsThickness + 1.0 * cm; // extra cm added to anode thickness
  double top_vessel_3 = LArGarIntefaceZ;
  double top_vessel_2 = LArGarIntefaceZ - myWindowsThickness;
  double top_vessel_1 = LArGarIntefaceZ - myTPCHeight - myWindowsThickness;

  double rad_vessel_4 = myReflectorInnerRadius + myReflectorThickness;
  double rad_vessel_3 = myReflectorInnerRadius + myReflectorThickness;
  if (remove_bonding_ring) rad_vessel_3 -= myReflectorThickness;
  if (remove_bonding_ring) rad_vessel_4 -= myReflectorThickness;

  // double rad_vessel_2 = myReflectorInnerRadius;
  // double rad_vessel_1 = myReflectorInnerRadius;

  //-------------------------           //
  //   Skin of TPC Vessel               //
  //-------------------------           //

  double top_skin_6 = top_vessel_4 + myElectronicsSpaceTop + myGdThicknessShield + myPENThickness + 2. * mm;
  double top_skin_5 = top_vessel_4;
  double top_skin_4 = top_vessel_4;
  double top_skin_3 = top_vessel_3;
  double top_skin_2 = top_vessel_2;
  double top_skin_1 = top_vessel_1 - myElectronicsSpaceBottom - myGdThicknessShield - myPENThickness - 2. * mm;

  double myZSkin[6] = {top_skin_1, top_skin_2, top_skin_3, top_skin_4, top_skin_5, top_skin_6};
  double myInnerRadius[6] = {0, 0, 0, 0, 0, 0};

  double rad_skin_4 = myNSLArThickness + myTPBThickness + myPENThickness + rad_vessel_4 + 4.4 * cm;
  double rad_skin_3 = myNSLArThickness + myTPBThickness + myPENThickness + rad_vessel_3 + 4.4 * cm;
  double rad_skin_2 = myNSLArThickness + myTPBThickness + myPENThickness + myReflectorInnerRadius + 4.4 * cm;
  double rad_skin_1 = myNSLArThickness + myTPBThickness + myPENThickness + myReflectorInnerRadius + 4.4 * cm;
  double myOuterRadius[6] = {rad_skin_1, rad_skin_2, rad_skin_3, rad_skin_4, rad_skin_2, rad_skin_2};

  // double rad_TPB_4                     = myTPBThickness +  rad_vessel_4;
  // double rad_TPB_3                     = myTPBThickness +  rad_vessel_3;
  // double rad_TPB_2                     = myTPBThickness +
  // myReflectorInnerRadius  ;
  //double rad_TPB_1 = myTPBThickness + myReflectorInnerRadius;

  double top_GdAcrylic_6 = top_vessel_4 + myElectronicsSpaceTop + myGdThicknessShield;
  double top_GdAcrylic_5 = top_vessel_4;
  double top_GdAcrylic_4 = top_vessel_4;
  double top_GdAcrylic_3 = top_vessel_3;
  double top_GdAcrylic_2 = top_vessel_2;
  double top_GdAcrylic_1 = top_vessel_1 - myElectronicsSpaceBottom - myGdThicknessShield;
  double myZGdAcrylic[6] = {top_GdAcrylic_1, top_GdAcrylic_2, top_GdAcrylic_3, top_GdAcrylic_4, top_GdAcrylic_5, top_GdAcrylic_6};

  double rad_GdAcrylic_4 = myTPBThickness + rad_vessel_4 + 4.4 * cm;
  double rad_GdAcrylic_3 = myTPBThickness + rad_vessel_3 + 4.4 * cm;
  double rad_GdAcrylic_2 = myTPBThickness + myReflectorInnerRadius + 4.4 * cm;
  double rad_GdAcrylic_1 = myTPBThickness + myReflectorInnerRadius + 4.4 * cm;
  double myOuterGdAcrylicRadius[6] = {rad_GdAcrylic_1, rad_GdAcrylic_2, rad_GdAcrylic_3, rad_GdAcrylic_4, rad_GdAcrylic_2, rad_GdAcrylic_2};

  double top_NSLAr_6 = top_GdAcrylic_6 + 2. * mm;
  double top_NSLAr_5 = top_GdAcrylic_4;
  double top_NSLAr_4 = top_GdAcrylic_4;
  double top_NSLAr_3 = top_GdAcrylic_3;
  double top_NSLAr_2 = top_GdAcrylic_2;
  double top_NSLAr_1 = top_GdAcrylic_1 - 2. * mm;
  double myZNSLAr[6] = {top_NSLAr_1, top_NSLAr_2, top_NSLAr_3, top_NSLAr_4, top_NSLAr_5, top_NSLAr_6};

  double rad_NSLAr_4 = myNSLArThickness + rad_GdAcrylic_4;
  double rad_NSLAr_3 = myNSLArThickness + rad_GdAcrylic_3;
  double rad_NSLAr_2 = myNSLArThickness + rad_GdAcrylic_2;
  double rad_NSLAr_1 = myNSLArThickness + rad_GdAcrylic_1;
  double myOuterNSLArRadius[6] = {rad_NSLAr_1, rad_NSLAr_2, rad_NSLAr_3, rad_NSLAr_4, rad_NSLAr_2, rad_NSLAr_2};

  double myEpsilon = 0.0 * mm;

  ////////////////////////////////////////////

  ///////     Begin of TPC volumes     ///////

  ///////////////////////////////////////////

  // TPC PEN skin

  G4Polyhedra* fSolidPENSkin = new G4Polyhedra("PENSkin_Solid", 0, myTwoPi, 8, 6, myZSkin, myInnerRadius, myOuterRadius);

  G4LogicalVolume* fLogicPENSkin = new G4LogicalVolume(fSolidPENSkin, DSMaterial::Get()->GetPEN(), "PENSkin_Logic");

  myPENReflectorSkin = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), "PENSkin", fLogicPENSkin, fMotherVolume, false, 0, myCheckOverlap);

  // TPC NSLAr skin

  G4Polyhedra* fSolidNSLArSkin = new G4Polyhedra("NSLArSkin_Solid", 0, myTwoPi, 8, 6, myZNSLAr, myInnerRadius, myOuterNSLArRadius);

  G4LogicalVolume* fLogicNSLArSkin = new G4LogicalVolume(fSolidNSLArSkin, DSMaterial::Get()->GetNSLiquidArgon(), "NSLArSkin_Logic");

  myNSLArReflectorSkin = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), "NSLArSkin", fLogicNSLArSkin, myPENReflectorSkin, false, 0, myCheckOverlap);

  // TPC Gd-Acrylic walls

  G4Polyhedra* fSolidGdSkin = new G4Polyhedra("GdSkin_Solid", 0, myTwoPi, 8, 6, myZGdAcrylic, myInnerRadius, myOuterGdAcrylicRadius);

  G4LogicalVolume* fLogicGdSkin; 

  if (DSStorage::Get()->GetPurePMMATPC() || DSStorage::Get()->GetHybridTPC()) {
  	fLogicGdSkin = new G4LogicalVolume(fSolidGdSkin, DSMaterial::Get()->GetAcrylic(), "GdSkin_Logic");
  }
  else fLogicGdSkin = new G4LogicalVolume(fSolidGdSkin, DSMaterial::Get()->GetGdAcrylic(), "GdSkin_Logic");

  myGdReflectorSkin = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), "GdSkin", fLogicGdSkin, myNSLArReflectorSkin, false, 0, myCheckOverlap);

  // Thickness of the optical plane structure (inserted below)

  G4double myOPSteelRingThickness = 1.83 * cm; 

  // Placeholder for changing OPs material
  G4double myZPlaceholderTop[2] = {top_GdAcrylic_6 - myGdThicknessShield, top_GdAcrylic_6};
  G4double myOuterPlaceholderRadius[2] = {rad_GdAcrylic_1 - myOPSteelRingThickness, rad_GdAcrylic_1 - myOPSteelRingThickness};
  G4double myInnerPlaceholderRadius[2] = {0,0};
  G4Polyhedra* fSolidPlaceholderTop = new G4Polyhedra("PlaceholderTop_Solid", 0, myTwoPi, 8, 2, myZPlaceholderTop, myInnerPlaceholderRadius, myOuterPlaceholderRadius); 
  
  G4LogicalVolume* fLogicPlaceholderTop;
  if (DSStorage::Get()->GetPurePMMATPC()) {
  	fLogicPlaceholderTop = new G4LogicalVolume(fSolidPlaceholderTop, DSMaterial::Get()->GetAcrylic(), "TopPlaceholder_Logic");
  }
  else fLogicPlaceholderTop = new G4LogicalVolume(fSolidPlaceholderTop, DSMaterial::Get()->GetGdAcrylic(), "TopPlaceholder_Logic");

  G4VPhysicalVolume* fPhysPlaceholderTop = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), "PlaceholderTop", fLogicPlaceholderTop, myGdReflectorSkin, false, 0, myCheckOverlap);

  G4double myZPlaceholderBot[2] = {top_GdAcrylic_1, top_GdAcrylic_1 + myGdThicknessShield};
  G4Polyhedra* fSolidPlaceholderBot = new G4Polyhedra("PlaceholderBot_Solid", 0, myTwoPi, 8, 2, myZPlaceholderBot, myInnerPlaceholderRadius, myOuterPlaceholderRadius);
  
  G4LogicalVolume* fLogicPlaceholderBot;
  if (DSStorage::Get()->GetPurePMMATPC()) {
  fLogicPlaceholderBot = new G4LogicalVolume(fSolidPlaceholderBot, DSMaterial::Get()->GetAcrylic(), "BotPlaceholder_Logic");
  }
  else fLogicPlaceholderBot = new G4LogicalVolume(fSolidPlaceholderBot, DSMaterial::Get()->GetGdAcrylic(), "BotPlaceholder_Logic"); 

  G4VPhysicalVolume* fPhysPlaceholderBot = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), "PlaceholderBot", fLogicPlaceholderBot, myGdReflectorSkin, false, 0, myCheckOverlap);







































  //-------------------------           //
  //    TPC Vessel, including  windows //
  //-------------------------           //

  double myAddCathodeThickness = 1.5 * cm; // adds more thickness at the bottom-center of the cathode
  double myFlatBotCathodeExtent = 20. * cm; // Extent of the bottom flat (octagonal) surface of the cathode
  double top_vessel_0 = top_vessel_1 - myAddCathodeThickness;
  double myZAcrylicTPCVessel[5] = {top_vessel_0, top_vessel_1, top_vessel_2, top_vessel_3, top_vessel_4};
  double myInnerRadiusReflector[5] = {0, 0, 0, 0, 0};
  double myOuterRadiusTPBSide[2] = {myTPCInnerRadius + myTPBThickness, myTPCInnerRadius + myTPBThickness};
  double myOuterRadiusReflector[5] = {myFlatBotCathodeExtent, myOuterRadiusTPBSide[0] + 4. * mm + 4. * cm + myEpsilon, myOuterRadiusTPBSide[1] + 4. * mm + 4. * cm + myEpsilon, myOuterRadiusTPBSide[1] + 4. * mm + 4. * cm + myEpsilon, myOuterRadiusTPBSide[1] + 4. * mm + 4. * cm + myEpsilon};

  G4Polyhedra* fSolidAcrylicTPCVessel = new G4Polyhedra("AcrylicTPCVessel_Solid", 0, myTwoPi, 8, 5, myZAcrylicTPCVessel, myInnerRadiusReflector, myOuterRadiusReflector);

  G4LogicalVolume* fLogicAcrylicTPCVessel = new G4LogicalVolume(fSolidAcrylicTPCVessel, DSMaterial::Get()->GetAcrylic(), "AcrylicTPCVessel_Logic");

  fPhysicAcrylicTPCVessel = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), "AcrylicTPCVessel", fLogicAcrylicTPCVessel, myGdReflectorSkin, false, 0, myCheckOverlap);

  fLogicAcrylicTPCVessel->SetVisAttributes(myGreen);
  fLogicAcrylicTPCVessel->SetVisAttributes(myBlue);

  G4double myGastoGridDistance = 3.0 * mm;
  double myZTPBSide[2] = {LArGarIntefaceZ - myTPCHeight, LArGarIntefaceZ - myGastoGridDistance};
  double myInnerRadiusTPC[2] = {0, 0};
  
  // Pure acrylic top windows extension inside the Gd acrylic
  double ZTopWindowExtension[2] = {top_vessel_4, top_vessel_4 - 6 * cm};
  double ROutTopWindowExtension[2] = {rad_GdAcrylic_1, rad_GdAcrylic_1};
  double RinTopWindowExtension[2] = {myOuterRadiusReflector[1], myOuterRadiusReflector[1]};

  G4Polyhedra* fSolidAcrylicTopWindowExtension = new G4Polyhedra("AcrylicTopWindowExtension_Solid", 0, myTwoPi, 8, 2, ZTopWindowExtension, RinTopWindowExtension, ROutTopWindowExtension);

  G4LogicalVolume* fLogicAcrylicTopWindowExtension = new G4LogicalVolume(fSolidAcrylicTopWindowExtension, DSMaterial::Get()->GetAcrylic(), "AcrylicTopWindowExtension_Logic");

  new G4PVPlacement(0, G4ThreeVector(0, 0, 0), "AcrylicTopWindowExtension", fLogicAcrylicTopWindowExtension, myGdReflectorSkin, false, 0, myCheckOverlap);

  // Pure acrylic bottom window extension inside the Gd acrylic
  double ZBottomWindowExtension[2] = {top_vessel_1, top_vessel_1 + 63 * cm};
  double ROutBottomWindowExtension[2] = {myOuterRadiusReflector[1] + 5 * cm, myOuterRadiusReflector[1] + 5 * cm};
  double RinBottomWindowExtension[2] = {myOuterRadiusReflector[1], myOuterRadiusReflector[1]};

  G4Polyhedra* fSolidAcrylicBottomWindowExtension = new G4Polyhedra("AcrylicBottomWindowExtension_Solid", 0, myTwoPi, 8, 2, ZBottomWindowExtension, RinBottomWindowExtension, ROutBottomWindowExtension);

  G4LogicalVolume* fLogicAcrylicBottomWindowExtension = new G4LogicalVolume(fSolidAcrylicBottomWindowExtension, DSMaterial::Get()->GetAcrylic(), "AcrylicBottomWindowExtension_Logic");

  new G4PVPlacement(0, G4ThreeVector(0, 0, 0), "AcrylicBottomWindowExtension", fLogicAcrylicBottomWindowExtension, myGdReflectorSkin, false, 0, myCheckOverlap);

  //--------------------------------------------------//
  //     Acrylic sandwich reflector   - passive Ar    //
  //--------------------------------------------------//

  double myOuterRadiusTPBSideNSLAr[2] = {myOuterRadiusTPBSide[0] + 4 * mm + 4. * cm, myOuterRadiusTPBSide[1] + 4. * mm + 4. * cm};

  G4Polyhedra* fSolidLayerNSLAr = new G4Polyhedra("LayerNSLAr_Solid", 0, myTwoPi, 8, 2, myZTPBSide, myInnerRadiusTPC, myOuterRadiusTPBSideNSLAr);

  G4LogicalVolume* fLogicLayerNSLAr = new G4LogicalVolume(fSolidLayerNSLAr, DSMaterial::Get()->GetNSLiquidArgon(), "LayerNSLAr_Logic");

  fPhysicLayerNSLAr = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), "LayerNSLAr", fLogicLayerNSLAr, fPhysicAcrylicTPCVessel, false, 0, myCheckOverlap);


  if ( DSStorage::Get()->GetDS20kReflectorAlternate() == 1 )  {
    //--------------------------------------------------//
    //     Acrylic sandwich reflector   - ESR    //
    //--------------------------------------------------//

    double myOuterRadiusTPBSideESR[2]  = {myOuterRadiusTPBSide[0] + 4*mm + 1.06*mm , myOuterRadiusTPBSide[1] + 4.*mm + 1.06*mm };

    G4Polyhedra *fSolidLayerESR = new G4Polyhedra( "LayerESR_Solid",
                        0,
                        myTwoPi,
                        8,
                        2,
                        myZTPBSide,
                        myInnerRadiusTPC,
                        myOuterRadiusTPBSideESR );



    G4LogicalVolume *fLogicLayerESR= new G4LogicalVolume(fSolidLayerESR,
                                   DSMaterial::Get()->GetESR(), // ESR
                                   "LayerESR_Logic");
    fLogicLayerESR->SetUserLimits(myStepLimit);



    fPhysicLayerESR  = new G4PVPlacement(0,
                            G4ThreeVector(0,0,0),
                   "LayerESR",
                   fLogicLayerESR,
                   fPhysicLayerNSLAr,
                   false,
                   0,
                   myCheckOverlap);

    //--------------------------------------------------//
    //     Acrylic sandwich reflector   - Passive Ar 2  //
    //--------------------------------------------------//

    double myOuterRadiusTPBSidePassiveLAr2[2]  = {myOuterRadiusTPBSide[0] + 4*mm + 1.*mm , myOuterRadiusTPBSide[1] + 4.*mm + 1.0*mm };

    G4Polyhedra *fSolidLayerPassiveLAr2 = new G4Polyhedra( "LayerPassiveLAr2_Solid",
                        0,
                        myTwoPi,
                        8,
                        2,
                        myZTPBSide,
                        myInnerRadiusTPC,
                        myOuterRadiusTPBSidePassiveLAr2 );



    G4LogicalVolume *fLogicLayerPassiveLAr2= new G4LogicalVolume(fSolidLayerPassiveLAr2,
                                   DSMaterial::Get()->GetNSLiquidArgon(),
                                   "LayerPassiveLAr2_Logic");
    fLogicLayerPassiveLAr2->SetUserLimits(myStepLimit);


    fPhysicLayerPassiveLAr2  = new G4PVPlacement(0,
                            G4ThreeVector(0,0,0),
                   "LayerPassiveLAr2",
                   fLogicLayerPassiveLAr2,
                   fPhysicLayerESR,
                   false,
                   0,
                   myCheckOverlap);


  //------------------------------------------------------------//
  //     Acrylic sandwich reflector   -  TPB support acrylic         //
  //------------------------------------------------------------//

  double myOuterRadiusTPBSideAcrylicSand[2] = {myOuterRadiusTPBSide[0] + 4 * mm, myOuterRadiusTPBSide[1] + 4 * mm};
  G4Polyhedra* fSolidLayerAcrylicSand = new G4Polyhedra("LayerAcrylicSand_Solid", 0, myTwoPi, 8, 2, myZTPBSide, myInnerRadiusTPC, myOuterRadiusTPBSideAcrylicSand);

  G4LogicalVolume* fLogicLayerAcrylicSand = new G4LogicalVolume(fSolidLayerAcrylicSand, DSMaterial::Get()->GetDS20kPlasticScintillator(), "LayerAcrylicSand_Logic");
  fLogicLayerAcrylicSand->SetUserLimits(myStepLimit);
  fPhysicLayerAcrylicSand = new G4PVPlacement(0,  G4ThreeVector(0, 0, 0),
                                                  "LayerAcrylicSand",
                                                  fLogicLayerAcrylicSand,
                                                  fPhysicLayerPassiveLAr2,
                                                  false, 0, myCheckOverlap);

  } else {

    //baseline reflector

    //------------------------------------------------------------//
    //     Acrylic sandwich reflector   - passive acrylic         //
    //------------------------------------------------------------//

    double myOuterRadiusTPBSideAcrylicSand[2] = {myOuterRadiusTPBSide[0] + 5.06 * mm, myOuterRadiusTPBSide[1] + 5.06 * mm};
    G4Polyhedra* fSolidLayerAcrylicSand = new G4Polyhedra("LayerAcrylicSand_Solid", 0, myTwoPi, 8, 2, myZTPBSide, myInnerRadiusTPC, myOuterRadiusTPBSideAcrylicSand);

    G4LogicalVolume* fLogicLayerAcrylicSand = new G4LogicalVolume(fSolidLayerAcrylicSand, DSMaterial::Get()->GetDS20kPlasticScintillator(), "LayerAcrylicSand_Logic");
    fLogicLayerAcrylicSand->SetUserLimits(myStepLimit);

    fPhysicLayerAcrylicSand = new G4PVPlacement(0,  G4ThreeVector(0, 0, 0),
                                                    "LayerAcrylicSand",
                                                    fLogicLayerAcrylicSand,
                                                    fPhysicLayerNSLAr,
                                                    false, 0, myCheckOverlap);

    //--------------------------------------------------//
    //     Acrylic sandwich reflector   - Passive Ar 2  //
    //--------------------------------------------------//

    double myOuterRadiusTPBSidePassiveLAr2[2]  = {myOuterRadiusTPBSide[0] + 1.06*mm , myOuterRadiusTPBSide[1] + 1.06*mm };

    G4Polyhedra *fSolidLayerPassiveLAr2 = new G4Polyhedra( "LayerPassiveLAr2_Solid",
                        0,
                        myTwoPi,
                        8,
                        2,
                        myZTPBSide,
                        myInnerRadiusTPC,
                        myOuterRadiusTPBSidePassiveLAr2 );



    G4LogicalVolume *fLogicLayerPassiveLAr2= new G4LogicalVolume(fSolidLayerPassiveLAr2,
                                   DSMaterial::Get()->GetNSLiquidArgon(),
                                   "LayerPassiveLAr2_Logic");
    fLogicLayerPassiveLAr2->SetUserLimits(myStepLimit);

    fPhysicLayerPassiveLAr2  = new G4PVPlacement(0,
                            G4ThreeVector(0,0,0),
                   "LayerPassiveLAr2",
                   fLogicLayerPassiveLAr2,
                   fPhysicLayerAcrylicSand,
                   false,
                   0,
                   myCheckOverlap);

    //--------------------------------------------------//
    //     Acrylic sandwich reflector   - ESR    //
    //--------------------------------------------------//

    double myOuterRadiusTPBSideESR[2]  = {myOuterRadiusTPBSide[0] + .06*mm , myOuterRadiusTPBSide[1] + .06*mm };

    G4Polyhedra *fSolidLayerESR = new G4Polyhedra( "LayerESR_Solid",
                       0,
                       myTwoPi,
                       8,
                       2,
                       myZTPBSide,
                       myInnerRadiusTPC,
                       myOuterRadiusTPBSideESR );



    G4LogicalVolume *fLogicLayerESR= new G4LogicalVolume(fSolidLayerESR,
                                  DSMaterial::Get()->GetESR(),
                                  "LayerESR_Logic");
    fLogicLayerESR->SetUserLimits(myStepLimit);
    fPhysicLayerESR  = new G4PVPlacement(0,
                           G4ThreeVector(0,0,0),
                  "LayerESR",
                  fLogicLayerESR,
                  fPhysicLayerPassiveLAr2,
                  false,
                  0,
                  myCheckOverlap);


  }
  //-----------------------------------------------------//
  //         TPB coating on the lateral surfaces         //
  //-----------------------------------------------------//

  // TPB layer (side )
  G4Polyhedra* fSolidTPBSide = new G4Polyhedra("TPBSide_Solid", 0, myTwoPi, 8, 2, myZTPBSide, myInnerRadiusTPC, myOuterRadiusTPBSide);

  G4LogicalVolume* fLogicTPBSide = new G4LogicalVolume(fSolidTPBSide, DSMaterial::Get()->GetTPB(), "TPBSide_Logic");
  fLogicTPBSide->SetUserLimits(myStepLimit);
  fPhysicTPBSide = new G4PVPlacement(0, G4ThreeVector(0, 0, 0),
                                        "TPBSide",
                                        fLogicTPBSide,
                                        DSStorage::Get()->GetDS20kReflectorAlternate() == 0 ? fPhysicLayerESR : fPhysicLayerAcrylicSand,
                                        //fPhysicLayerAcrylicSand,
                                        false, 0, myCheckOverlap);

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // TPC Active Volume
  //////////////////////////////////////////////////////////////////////////////////////////////////

  double myZTPCBottom[2] = {LArGarIntefaceZ - myTPCHeight, LArGarIntefaceZ - myGastoGridDistance};
  double myOuterRadiusTPC[2] = {myTPCInnerRadius, myTPCInnerRadius};
  DSLog(debugging) << "UAr Inner Radius = " << myTPCInnerRadius << endlog;
  G4Polyhedra* fSolidActiveLAr = new G4Polyhedra("ActiveLAr_Solid", 0, myTwoPi, 8, 2, myZTPCBottom, myInnerRadiusTPC, myOuterRadiusTPC);

  G4LogicalVolume* fLogicActiveLAr = new G4LogicalVolume(fSolidActiveLAr, DSMaterial::Get()->GetLiquidArgon(), "LAr_Logic");
  //fLogicActiveLAr->SetUserLimits(myStepLimit);

  fPhysicActiveLAr = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), "ActiveLAr", fLogicActiveLAr, fPhysicTPBSide, false, 0, myCheckOverlap);

  fLogicActiveLAr->SetVisAttributes(myBlue);
  G4Region* fLArRegion = new G4Region("LAr_Logic");
  fLogicActiveLAr->SetRegion(fLArRegion);
  fLArRegion->AddRootLogicalVolume(fLogicActiveLAr);

  DSStorage::Get()->SetLiquidArgonIndex((int)fLogicActiveLAr->GetMaterial()->GetIndex());

  ///////////////////////////////////////////////////////////////////////////////////////
  // Grid
  ///////////////////////////////////////////////////////////////////////////////////////

  G4double myGrid_h = 0.02 * cm;
  double myZGrid[2] = {LArGarIntefaceZ - myGastoGridDistance, LArGarIntefaceZ - myGastoGridDistance + myGrid_h};
  
  fSolidGrid = new G4Polyhedra("Grid_Solid", 0, myTwoPi, 8, 2, myZGrid, myInnerRadiusTPC, myOuterRadiusTPBSideNSLAr);

  fLogicGrid = new G4LogicalVolume(fSolidGrid, DSMaterial::Get()->GetGridSteel(), "Grid_Logic");
  fPhysicGrid = new G4PVPlacement(0, myZeros, "Grid", fLogicGrid, fPhysicAcrylicTPCVessel, false, 0, myCheckOverlap);

  ///////////////////////////////////////////////////////////////////////////////////////
  // LAr above Grid
  ///////////////////////////////////////////////////////////////////////////////////////

  double myZLArAboveGrid[2] = {LArGarIntefaceZ - myGastoGridDistance + myGrid_h, LArGarIntefaceZ};

  G4Polyhedra*  fSolidLArAboveGrid = new G4Polyhedra("LArAboveGrid_Solid", 0, myTwoPi, 8, 2, myZLArAboveGrid, myInnerRadiusTPC, myOuterRadiusTPBSideNSLAr);

  G4LogicalVolume* fLogicGridLArAboveGrid = new G4LogicalVolume(fSolidLArAboveGrid, DSMaterial::Get()->GetLArAboveGrid(), "LArAboveGrid_Logic");
  fPhysicLArAboveGrid = new G4PVPlacement(0, myZeros, "LArAboveGrid", fLogicGridLArAboveGrid, fPhysicAcrylicTPCVessel, false, 0, myCheckOverlap);
  
  //-------------------------------------------//
  //  Active GAr Region                        //
  //-------------------------------------------//

  double myZTPCTop[2] = {LArGarIntefaceZ, LArGarIntefaceZ + myActiveGasLogicThickness};

  G4Polyhedra* fSolidGasPocket = new G4Polyhedra("GasPocket_Solid", 0, myTwoPi, 8, 2, myZTPCTop, myInnerRadiusTPC, myOuterRadiusTPBSideNSLAr);

  G4LogicalVolume* fLogicGasPocket = new G4LogicalVolume(fSolidGasPocket, DSMaterial::Get()->GetGaseousArgon(), "GasPocket_Logic");

  fPhysicGasPocket = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), "GasPocket", fLogicGasPocket, fPhysicAcrylicTPCVessel, false, 0, myCheckOverlap);

  double myZTPBTop[2] = {myLArGArBoundaryPosZ + myGasPocketThickness + myPArThickness, myLArGArBoundaryPosZ + myGasPocketThickness + myPArThickness + myTPBThickness};
  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Top window TPB coating
  //////////////////////////////////////////////////////////////////////////////////////////////////
  G4Polyhedra* fSolidTPBTop = new G4Polyhedra("TPBTop_Solid", 0, myTwoPi, 8, 2, myZTPBTop, myInnerRadiusTPC, myOuterRadiusTPBSideNSLAr);

  G4LogicalVolume* fLogicTPBTop = new G4LogicalVolume(fSolidTPBTop, DSMaterial::Get()->GetTPB(), "TPBTop_Logic");

  fPhysicTPBTop = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), "TPBTop", fLogicTPBTop, fPhysicAcrylicTPCVessel, false, 0, myCheckOverlap);

  double myZTPBBot[2] = {LArGarIntefaceZ - myTPCHeight - myTPBThickness, LArGarIntefaceZ - myTPCHeight};
  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Bottom window TPB coating
  //////////////////////////////////////////////////////////////////////////////////////////////////
  G4Polyhedra* fSolidTPBBot = new G4Polyhedra("TPBBot_Solid", 0, myTwoPi, 8, 2, myZTPBBot, myInnerRadiusTPC, myOuterRadiusTPC);

  G4LogicalVolume* fLogicTPBBot = new G4LogicalVolume(fSolidTPBBot, DSMaterial::Get()->GetTPB(), "TPBBot_Logic");

  fPhysicTPBBottom = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), "TPBBot", fLogicTPBBot, fPhysicAcrylicTPCVessel, false, 0, myCheckOverlap);

  ///////////////////////////////////
  ///////////////////////////////////
  //  END OF TPC
  ///////////////////////////////////
  ///////////////////////////////////

  //--------------------------------------------------//
  //     Passive argon above and below the TPC       //
  //--------------------------------------------------//
  // these volumes are containers of all the electronics and mechanics. Nothing
  // will be outside. The construction philosphy will also be slighly different:
  // the system is going to be symmetric: only one solid will be created, and
  // then placed above/below the TPC at the correct Z. The reference Z are
  // defined just below

  double myZpassiveAbove[2] = {top_vessel_4, top_vessel_4 + myElectronicsSpaceTop};
  double myZpassiveBelow[4] = {top_vessel_1 - myElectronicsSpaceBottom, top_vessel_1, top_vessel_1, top_vessel_0};

  double myInnerRadiusPassive[2] = {0, 0};
  double myInnerRadiusPassiveBelow[4] = {0, 0, 0, 0};
  double myOuterRadiusPassiveAbove[2] = {rad_GdAcrylic_1 - myOPSteelRingThickness, rad_GdAcrylic_1 - myOPSteelRingThickness};
  double myOuterRadiusPassiveBelow[4] = {rad_GdAcrylic_1 - myOPSteelRingThickness, rad_GdAcrylic_1 - myOPSteelRingThickness, myOuterRadiusTPBSide[0] + 4. * mm + 4. * cm + myEpsilon, myFlatBotCathodeExtent};

  G4Polyhedra* fSolidPassiveBot = new G4Polyhedra("PassiveBot_Solid", 0, myTwoPi, 8, 4, myZpassiveBelow, myInnerRadiusPassiveBelow, myOuterRadiusPassiveBelow);
  
  G4LogicalVolume* fLogicPassiveBot = new G4LogicalVolume(fSolidPassiveBot, DSMaterial::Get()->GetNSLiquidArgon(), "PassiveBot_Logic");

  fPhysicPassiveBottom = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), "PassiveBot", fLogicPassiveBot, myGdReflectorSkin, false, 0, myCheckOverlap);

  G4Polyhedra* fSolidPassiveTop = new G4Polyhedra("PassiveTop_Solid", 0, myTwoPi, 8, 2, myZpassiveAbove, myInnerRadiusPassive, myOuterRadiusPassiveAbove);

  G4LogicalVolume* fLogicPassiveTop = new G4LogicalVolume(fSolidPassiveTop, DSMaterial::Get()->GetNSLiquidArgon(), "PassiveTop_Logic");

  fPhysicPassiveTop = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), "PassiveTop", fLogicPassiveTop, myGdReflectorSkin, false, 0, myCheckOverlap);

  double myZSiPm[2] = {-0.5 * mySiPmThickness, 0.5 * mySiPmThickness};
  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Top Array of Silicon
  G4Polyhedra* fSolidSiPMTop = new G4Polyhedra("SiPMTop_Solid", 0, myTwoPi, 8, 2, myZSiPm, myInnerRadiusTPC, myOuterRadiusTPBSideNSLAr);

  G4LogicalVolume* fLogicSiPMTop = new G4LogicalVolume(fSolidSiPMTop, DSMaterial::Get()->GetMetalSilicon(), "SiPMTop_Logic");

  fPhysicSiPmTop = new G4PVPlacement(0, G4ThreeVector(0, 0, myZpassiveAbove[0] + mySiPMOffset_top + 0.5 * mySiPmThickness), "SiPMTop", fLogicSiPMTop, fPhysicPassiveTop, false, -56, myCheckOverlap);

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Bottom Array of Silicon
  fPhysicSiPmBottom = new G4PVPlacement(0, G4ThreeVector(0, 0, myZpassiveBelow[1] - mySiPMOffset_bottom - 0.5 * mySiPmThickness), "SiPMBottom", fLogicSiPMTop, fPhysicPassiveBottom, false, -56, myCheckOverlap);

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Arlon PCB
  // /////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////
  // PCB stackup comprise Arlon 55NT prepreg 1.1 mm total + 0.1 mm Copper (4 layers)
  // Implement a double stack of ArlonPrepreg + Copper to mimick the MB PCB - Tile PCB, in the middle 5 mm of NSLAr

  G4double myPCBCopperThickness = 0.1 * mm;
  G4double myPCBArlonPrepregThickness = 1.1 * mm;
  G4double myPCBLArThickness = 6.0 * mm;
  //  G4double PCB_thickness = mass_PCB_per_tile / DSMaterial::Get()->GetArlon()->GetDensity() / (25 * cm2);

  G4double myZPCBCopper[2] = {-0.5 * myPCBCopperThickness, 0.5 * myPCBCopperThickness};
  G4double Vertical_Offset = mySiPMOffset_top + mySiPmThickness + myPCBCopperThickness / 2.;  // compared to the reference Z: myZpassiveAbove[0] for
                                                                                              // top and myZpassiveBelow[1] for bot
  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Top Array of PCB Copper #1 (Tile)
  G4Polyhedra* fSolidCopperPCBs = new G4Polyhedra("CopperPCBs_Solid", 0, myTwoPi, 8, 2, myZPCBCopper, myInnerRadiusTPC, myOuterRadiusTPBSideNSLAr);

  G4LogicalVolume* fLogicCopperPCBs = new G4LogicalVolume(fSolidCopperPCBs, DSMaterial::Get()->GetMetalCopper(), "CopperPCBs_Logic");

  new G4PVPlacement(0, G4ThreeVector(0, 0, myZpassiveAbove[0] + Vertical_Offset), "CopperPCBTiles_top", fLogicCopperPCBs, fPhysicPassiveTop, false, 0, myCheckOverlap);

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Bottom Array of PCB Copper #1 (Tile)
  new G4PVPlacement(0, G4ThreeVector(0, 0, myZpassiveBelow[1] - Vertical_Offset + mySiPMOffset_mismatch), "CopperPCBTiles_bot", fLogicCopperPCBs, fPhysicPassiveBottom, false, 0, myCheckOverlap);

  //Arlon5NT prepreg layer
  G4double myZPCBArlon[2] = {-0.5 * myPCBArlonPrepregThickness, 0.5 * myPCBArlonPrepregThickness};

  Vertical_Offset += myPCBCopperThickness / 2. + myPCBArlonPrepregThickness / 2.; 
  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Top Array of PCB Arlon #1 (Tile)
  G4Polyhedra* fSolidArlonPCBs = new G4Polyhedra("ArlonPCBs_Solid", 0, myTwoPi, 8, 2, myZPCBArlon, myInnerRadiusTPC, myOuterRadiusTPBSideNSLAr);

  G4LogicalVolume* fLogicArlonPCBs = new G4LogicalVolume(fSolidArlonPCBs, DSMaterial::Get()->GetArlonPrepreg(), "ArlonPCBs_Logic");

  new G4PVPlacement(0, G4ThreeVector(0, 0, myZpassiveAbove[0] + Vertical_Offset), "ArlonPCBTiles_top", fLogicArlonPCBs, fPhysicPassiveTop, false, 0, myCheckOverlap);

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Bottom Array of PCB Arlon #1 (Tile)
  new G4PVPlacement(0, G4ThreeVector(0, 0, myZpassiveBelow[1] - Vertical_Offset + mySiPMOffset_mismatch), "ArlonPCBTiles_bot", fLogicArlonPCBs, fPhysicPassiveBottom, false, 0, myCheckOverlap);

  //Repeat for MB PCB adding NSLar in between

  Vertical_Offset += myPCBArlonPrepregThickness /2. + myPCBLArThickness + myPCBCopperThickness / 2.;  

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Top Array of PCB Copper #2 (MB)
  new G4PVPlacement(0, G4ThreeVector(0, 0, myZpassiveAbove[0] + Vertical_Offset), "CopperPCBMBs_top", fLogicCopperPCBs, fPhysicPassiveTop, false, 0, myCheckOverlap);

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Bottom Array of PCB Copper #2 (MB)
  new G4PVPlacement(0, G4ThreeVector(0, 0, myZpassiveBelow[1] - Vertical_Offset + mySiPMOffset_mismatch), "CopperPCBMB_bot", fLogicCopperPCBs, fPhysicPassiveBottom, false, 0, myCheckOverlap);

  //Arlon5NT prepreg layer
  Vertical_Offset += myPCBCopperThickness / 2. + myPCBArlonPrepregThickness / 2.; 
  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Top Array of PCB Arlon #2 (MB)
  new G4PVPlacement(0, G4ThreeVector(0, 0, myZpassiveAbove[0] + Vertical_Offset), "ArlonPCBMBs_top", fLogicArlonPCBs, fPhysicPassiveTop, false, 0, myCheckOverlap);

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Bottom Array of PCB Arlon #2 (MB)
  new G4PVPlacement(0, G4ThreeVector(0, 0, myZpassiveBelow[1] - Vertical_Offset + mySiPMOffset_mismatch), "ArlonPCBMBs_bot", fLogicArlonPCBs, fPhysicPassiveBottom, false, 0, myCheckOverlap);

  
    ////////////////////////////////////////////////////////////////////////////////////////////////
  // Stainless Steel optical plane structure   ////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////
  
  G4double myOPSteelRingTopHeight = myGdThicknessShield + myElectronicsSpaceTop;
  G4double myOPSteelRingBotHeight = myGdThicknessShield + myElectronicsSpaceBottom;
  G4double myZSteelTopRing[2] = {-myOPSteelRingTopHeight / 2., myOPSteelRingTopHeight / 2.};
  G4double myZSteelBotRing[2] = {-myOPSteelRingBotHeight / 2., myOPSteelRingBotHeight / 2.};
  G4double myOuterSteelRingRadius[2] = {rad_GdAcrylic_1, rad_GdAcrylic_1}; 
  G4double myInnerSteelRingRadius[2] = {rad_GdAcrylic_1 - myOPSteelRingThickness, rad_GdAcrylic_1 - myOPSteelRingThickness};

  ////////////////////////////////////////////////////////////////////////////////////
  // Steel Top ring
  G4Polyhedra* fSolidSteelTopRing = new G4Polyhedra("SteelTopRing_Solid", 0, myTwoPi, 8, 2, myZSteelTopRing, myInnerSteelRingRadius, myOuterSteelRingRadius);

  G4LogicalVolume* fLogicSteelTopRing = new G4LogicalVolume(fSolidSteelTopRing, DSMaterial::Get()->GetStainlessSteel(), "SteelTopRing_Logic");
		
  new G4PVPlacement(0, G4ThreeVector(0, 0, top_GdAcrylic_6 - myOPSteelRingTopHeight / 2.), "SteelTopRing", fLogicSteelTopRing, myGdReflectorSkin, false, 0, myCheckOverlap);

   ////////////////////////////////////////////////////////////////////////////////////
  // Steel Bottom ring
  G4Polyhedra* fSolidSteelBotRing = new G4Polyhedra("SteelBotRing_Solid", 0, myTwoPi, 8, 2, myZSteelBotRing, myInnerSteelRingRadius, myOuterSteelRingRadius);

  G4LogicalVolume* fLogicSteelBotRing = new G4LogicalVolume(fSolidSteelBotRing, DSMaterial::Get()->GetStainlessSteel(), "SteelBotRing_Logic");
		
  new G4PVPlacement(0, G4ThreeVector(0, 0, top_GdAcrylic_1 + myOPSteelRingBotHeight / 2.), "SteelBotRing", fLogicSteelBotRing, myGdReflectorSkin, false, 0, myCheckOverlap);

  G4double myMassSteelTopRing = DSMaterial::Get()->GetStainlessSteel()->GetDensity()  * 8. / (1 + sqrt(2.)) * myOPSteelRingThickness * myOPSteelRingTopHeight * (2. * myOuterSteelRingRadius[0] - myOPSteelRingThickness); // density times the volume of the octagonal ring

  //G4double myMassSteelBotRing = myMassSteelTopRing * myOPSteelRingBotHeight / myOPSteelRingTopHeight;

  G4double myMassSteelBeams = mass_steel_plates - myMassSteelTopRing;
     
  if (myMassSteelBeams < 0)DSLog(fatal) << "Mass of the steel ring is greater than total allowed mass =====>" << "Total mass(kg) = " << mass_steel_plates/kg << " Ring's mass(kg) = " <<  myMassSteelTopRing/kg << endlog;

  //cout << "Total mass(kg) = " << mass_steel_plates/kg << " Ring's mass(kg) = " <<  myMassSteelTopRing/kg << endl;
  
  G4double myLArSkinofSsBeam = 7.5 * mm; // twice the thickness of layers of NSLAr on each side of the beams
  G4int nOfBeams = 18;
  vector<double> mySsBeamLength;
  G4double myTotalSsBeamLength = 0.;
	
  G4ThreeVector p = {0, 0, 0};
  G4ThreeVector v = {0, 1, 0};
  v.rotateZ(22.5 * deg); // to match the coordinate system of the octagonal volume
  G4double myBeamSpacing = 2. * myInnerSteelRingRadius[0] / (nOfBeams + 1.);

  // this loop calculates the length of each beam and sum them
  for (int k = 0; k < nOfBeams; k++){
		  
    p ={myInnerSteelRingRadius[0] - (k + 1.)*myBeamSpacing, 0, 0};

    // rotate the vector to match the coordinate system of the octagonal volume
    p.rotateZ(22.5 * deg);
    mySsBeamLength.push_back(2. * fSolidSteelTopRing->DistanceToIn(p,v));
    myTotalSsBeamLength += mySsBeamLength[k];
  }

  // calculates the thickness of the beams given their total mass
  G4double mySsBeamThickness = (myMassSteelBeams) / (DSMaterial::Get()->GetStainlessSteel()->GetDensity() * myGdThicknessShield * myTotalSsBeamLength);

 
  G4RotationMatrix* rotSsBeam = new G4RotationMatrix();
  rotSsBeam->rotateZ(-22.5 * deg);
  
  // place the beams, each beam has a layer of LAr on both sides 
  for (int k = 0; k < nOfBeams; k++){
    
    p = {myInnerSteelRingRadius[0] - (k + 1.)*myBeamSpacing, 0, top_GdAcrylic_6 - myGdThicknessShield / 2.};
    
    if (p[0] > 0) {p[0] -=  (mySsBeamThickness + myLArSkinofSsBeam) / 2.;}
    else if (p[0] < 0){p[0] += (mySsBeamThickness + myLArSkinofSsBeam) / 2.;}
    
    G4Box* fSolidLArSkinBeam = new G4Box("LArSkinBeam_Solid", (mySsBeamThickness + myLArSkinofSsBeam) / 2., mySsBeamLength[k] / 2.,myGdThicknessShield / 2.);
    G4LogicalVolume* fLogicLArSkinBeam = new G4LogicalVolume(fSolidLArSkinBeam, DSMaterial::Get()->GetNSLiquidArgon(), "LArSkinBeam_Logic");
    G4Box* fSolidSsBeam = new G4Box("SsBeam_Solid", mySsBeamThickness / 2., mySsBeamLength[k] / 2.,myGdThicknessShield / 2.);
    G4LogicalVolume* fLogicSsBeam = new G4LogicalVolume(fSolidSsBeam, DSMaterial::Get()->GetStainlessSteel(), "SsBeam_Logic");

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Top steel beam
    new G4PVPlacement(0,G4ThreeVector(0, 0, 0), fLogicSsBeam, "SsBeam", fLogicLArSkinBeam, false, k, myCheckOverlap);
    
    new G4PVPlacement(rotSsBeam, p.rotateZ(22.5 * deg), "LArSkinBeamTop", fLogicLArSkinBeam, fPhysPlaceholderTop, false, k, myCheckOverlap);

     /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Bottom steel beam
    p.rotateZ(-22.5 * deg);
    p.setZ(top_GdAcrylic_1 + myGdThicknessShield / 2.);    
    new G4PVPlacement(rotSsBeam, p.rotateZ(22.5 * deg), "LArSkinBeamBot", fLogicLArSkinBeam, fPhysPlaceholderBot, false, k, myCheckOverlap);

  }

  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Stainless Steel grid support frame ////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////

  G4double myMassGridSupport = 126 * kg;
  G4double myOuterGridSupportRadius[2] = {rad_GdAcrylic_1 - 2.5 * cm, rad_GdAcrylic_1 - 2.5 * cm};
  G4double myInnerGridSupportRadius[2] = {myOuterRadiusReflector[1] + 2.5 * cm, myOuterRadiusReflector[1] + 2.5 * cm};

  G4double myGridSupportHeight = (myMassGridSupport / DSMaterial::Get()->GetStainlessSteel()->GetDensity()) * (1 + sqrt(2)) / (8 * (myOuterGridSupportRadius[0] * myOuterGridSupportRadius[0] - myInnerGridSupportRadius[0] * myInnerGridSupportRadius[0]));

  G4double myZGridSupport[2] = {LArGarIntefaceZ - myGastoGridDistance - myGridSupportHeight, LArGarIntefaceZ - myGastoGridDistance};
  
  G4Polyhedra* fSolidGridSupport = new G4Polyhedra("GridSupport_Solid", 0, myTwoPi, 8, 2, myZGridSupport, myInnerGridSupportRadius, myOuterGridSupportRadius);

    G4LogicalVolume* fLogicGridSupport = new G4LogicalVolume(fSolidGridSupport, DSMaterial::Get()->GetStainlessSteel(), "GridSupport_Logic");

    new G4PVPlacement(0, G4ThreeVector(0, 0, 0), "GridSupport", fLogicGridSupport, myGdReflectorSkin, false, 0, myCheckOverlap);

   ////////////////////////////////

  ///////  Inner VETO  SiPMs     //////
  // PA 10/23 - updated the geometry:
  // simplified. add single tile readout, changed channel index to [10000,11919]
  // shift down the first vPDU row by 10 cm
  ////////////////////////////////

  double z_sipm_above = myZSkin[5];   // upper edge of TPC
  double z_sipm_up = myZSkin[5] - 10*cm;   // upper edge of TPC
  double z_sipm_low = myZSkin[0];  // lower edge of TPC

  double r_insc_sipm = myOuterRadius[5];
  //double r_circ_sipm = 2. * myOuterRadius[5] / (sqrt(2. + sqrt(2.)));

  DSStorage::Get()->SetDS20kTPCHeight_top(myZSkin[5]);
  DSStorage::Get()->SetDS20kTPCHeight_bottom(myZSkin[0]);
  DSStorage::Get()->SetDS20kTPCInscribedR(myOuterRadius[5]);

  double h_sipm = z_sipm_up - z_sipm_low;
  double edge_sipm = r_insc_sipm / 1.207;

  //const G4double nsipm = DSStorage::Get()->GetDS20knSiPMs();
  const G4double SiPMSide = DSStorage::Get()->GetDS20kSiPMSide();
  const G4double SiPMHeight = DSStorage::Get()->GetDS20kSiPMHeight();
  const G4double vPDUSide   = 20.2 *cm ;

  // Teflon
  G4Material* teflon_mat = DSMaterial::Get()->GetTeflon();
  // Metal silicon
  G4Material* metalsilicon_mat = DSMaterial::Get()->GetMetalSilicon();

  const G4double epsilonSiPM = DSStorage::Get()->GetDS20kepsilonSiPM();

  /////////////////////////////////////////////////////////////////////////

  double fill_factor = 0.95 ;

  ///////////////////////////////////////////////////////////////////////
  // SiPMs on sides
  ///////////////////////////////////////////////////////////////////////
  G4Box* vPDU = new G4Box("vPDU", SiPMHeight / 2., vPDUSide / 2., vPDUSide / 2.);
  myLogicShapevPDU  = new G4LogicalVolume(vPDU, teflon_mat, "vPDU");

  G4Box* vQuadrant = new G4Box("vQuadrant", SiPMHeight/2. , vPDUSide / 4., vPDUSide / 4.);
  myLogicShapevQuadrant  = new G4LogicalVolume(vQuadrant, teflon_mat, "vQuadrant");

  // G4Box* SiPMSilica = new G4Box("SiPMSilica", 0.1 * cm, 0.5 * (SiPMSide / 2.)* fill_factor, 0.5 * (SiPMSide / 2.) * fill_factor);
  G4Box* SiPMSilica = new G4Box("SiPMSilica", 0.1 * cm, (SiPMSide / 2.)* fill_factor, (SiPMSide / 2.) * fill_factor);
  myLogicShapevSiPMSilica = new G4LogicalVolume(SiPMSilica, metalsilicon_mat, "SiPMSilica");

  //populate one vPDU quadrant with 4 SiPM tiles
  G4double shift = - 0.5 * vPDUSide / 4. ;
  for (int i = 0; i < 4; i++) {
    int index_1 = i ;
    new G4PVPlacement(0, G4ThreeVector(    - SiPMHeight / 2. + 1*mm  , shift  +  (vPDUSide /4.) * double(i%2), shift +  (vPDUSide /4.) *  int((i%4)>1)   ), myLogicShapevSiPMSilica, "SiPMSilica", myLogicShapevQuadrant,false, index_1, myCheckOverlap);
  }


  shift = - 0.5 * vPDUSide /2. ;
  for (int i = 0; i < 4; i++) {

    int index_2 = i % 4;       // 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15...
    // physSiPMSilica[i] = new
    new G4PVPlacement(0, G4ThreeVector(   0 , shift  +  (vPDUSide /2.) * double(i%2),  shift +  (vPDUSide /2.) *  int((i%4)>1) ),myLogicShapevQuadrant , "SiPMSilica", myLogicShapevPDU,false, index_2, myCheckOverlap);

  }
  /*
  G4double baseArea = 2. * (2. * sqrt(2.) * r_circ_sipm * r_circ_sipm);  // total area of bases
  G4double latArea = 8. * edge_sipm * h_sipm;                            // total area of the faces

  G4int nSiPMBase = floor(1.5 * (baseArea / (baseArea + latArea) * nsipm) * 0.5);  // number of sipm on each base (now hard coded for vPDU+)
  G4int nSiPMLat = floor((nsipm - nSiPMBase * 2)) / 8.;                            // number of sipm on each face
  */ // older version code which allowed a variable number of SiPMs depending on dimensions

  ///////////////////////////////////////////////////////////////////////
  //////////////// SiPMs on octagonal bases ///////////////////////////
  ///////////////////////////////////////////////////////////////////////

  G4RotationMatrix* rotSiPM_top = new G4RotationMatrix();
  rotSiPM_top->rotateY(-90. * deg);
  rotSiPM_top->rotateX(22.5 * deg);
  rotSiPM_top->rotateX(45 * deg);

  G4RotationMatrix* rotSiPM_bot = new G4RotationMatrix();
  rotSiPM_bot->rotateY(90. * deg);
  rotSiPM_bot->rotateX(-22.5 * deg);
  rotSiPM_bot->rotateX(45 * deg);

  int ChannelIDShift = 10000 ; //shift the chan id by 10000 units, so that there is no ambiguity with TPC channel IDs
  int CopyNumber_top = ChannelIDShift;
  int CopyNumber_bottom = ChannelIDShift + 80*4 ;
  int CopyNumber_side = ChannelIDShift + 160*4;

  vector<double> x_pos;
  vector<double> y_pos;

  int N = 5;

  x_pos.push_back(-164. * cm);
  y_pos.push_back(-41. * cm);

  x_pos.push_back(-82. * cm);
  y_pos.push_back(-41. * cm);

  x_pos.push_back(-41. * cm);
  y_pos.push_back(-82. * cm);

  x_pos.push_back(-41. * cm);
  y_pos.push_back(-164. * cm);

  x_pos.push_back(-123. * cm);
  y_pos.push_back(-123. * cm);

  int count_up = 0;
  int count_down = 0;

  for (int i = 0; i < N; i++) {

    G4ThreeVector p_up = G4ThreeVector(x_pos.at(i), y_pos.at(i), z_sipm_above + epsilonSiPM + SiPMHeight / 2.);
    G4ThreeVector p_down = G4ThreeVector(x_pos.at(i), y_pos.at(i), z_sipm_low - epsilonSiPM - SiPMHeight / 2.);

    for (int j = 0; j < 4; j++) {

      // Top cap
      p_up.rotateZ(22.5 * deg + j * 90. * deg);

      if (DSStorage::Get()->GetDS20kSiPMs() > 0) {
        physSiPMInternalUpperCap[count_up] = new G4PVPlacement(rotSiPM_top, p_up, "SiPMInternalUpperCap", myLogicShapevPDU, fMotherVolume, false, CopyNumber_top, myCheckOverlap);
      }

      SiPMPosVector.push_back(p_up.rotateZ(-22.5 * deg));
      count_up++;
      CopyNumber_top += 16;

      // Bottom cap
      p_down.rotateZ(22.5 * deg + j * 90. * deg);

      if (DSStorage::Get()->GetDS20kSiPMs() > 0) {
        physSiPMInternalLowerCap[count_down] = new G4PVPlacement(rotSiPM_bot, p_down, "SiPMInternalLowerCap", myLogicShapevPDU, fMotherVolume, false, CopyNumber_bottom, myCheckOverlap);
      }

      SiPMPosVector.push_back(p_down.rotateZ(-22.5 * deg));
      count_down++;
      CopyNumber_bottom += 16;
    }
  }

  //////////////// SiPMs on Faces ///////////////////////////////////////////
  /* double R = edge_sipm / h_sipm;

  int nrow = 1.;
  int ncol = 1.;

  while (nrow * ncol < nSiPMLat) {

    nrow += 1;
    ncol = floor(R * nrow);
  }
  */ // older version code which allowed a variable number of SiPMs depending on dimensions

  // 5 x 2 arrays of SiPMs on each lateral wall
  G4int nrow = 5;
  G4int ncol = 2;

  int count_side = 0;

  for (int k = 0; k < 8; k++) {

    G4RotationMatrix* rotSiPMSide = new G4RotationMatrix();
    rotSiPMSide->rotateZ(-22.5 * deg - k * 45. * deg);
    rotSiPMSide->rotateY(180. * deg);

    for (int i = 0; i < nrow; i++) {
      for (int j = 0; j < ncol; j++) {

        double step_edge = (edge_sipm - ncol * vPDUSide) / (ncol);
        double step_h = (h_sipm - nrow * vPDUSide) / (nrow);

        double x = r_insc_sipm + epsilonSiPM + SiPMHeight / 2.;
        double y = -edge_sipm / 2. + vPDUSide * (2 * j + 1) / 2. + step_edge * (2 * j + 1) / 2.;
        double z = z_sipm_low + vPDUSide * (2 * i + 1) / 2. + step_h * (2 * i + 1) / 2.;

        G4ThreeVector posside = G4ThreeVector(x, y, z);
        posside.rotateZ(k * 45 * deg);
        posside.rotateZ(22.5 * deg);

        if (DSStorage::Get()->GetDS20kSiPMs() > 0) {
          physSiPMInternalSide[count_side] = new G4PVPlacement(rotSiPMSide, posside, "SiPMInternalSide", myLogicShapevPDU, fMotherVolume, false, CopyNumber_side, myCheckOverlap);
        }

        CopyNumber_side += 16;
        count_side++;
        SiPMPosVector.push_back(posside.rotateZ(-22.5 * deg));
      }
    }
  }



  DSStorage::Get()->SetDS20kSiPMPosVector(SiPMPosVector);
  DefineSurfaces();
  // make SiPM as pe storing material
  DSStorage::Get()->SetPMTMaterialIndex(fLogicSiPMTop->GetMaterial()->GetIndex());

}

DSDetectorDS20k::~DSDetectorDS20k() {
  ;
}

void DSDetectorDS20k::DefineSurfaces() {

  ////////////////////////////////////////
  // Copied from DS DetectorDS50 - 24th april 2015
  ////////////////////////////////////////

  // LAR-Grid
  G4OpticalSurface* fOpGridLArSurface = new G4OpticalSurface("OpGridLArSurface");
  new G4LogicalBorderSurface("GridLArSurface", fPhysicActiveLAr, fPhysicGrid, fOpGridLArSurface);
  fOpGridLArSurface->SetType(dielectric_dielectric);
  fOpGridLArSurface->SetModel(glisur);
  fOpGridLArSurface->SetFinish(polished);

  // LAr above grid - Grid
  G4OpticalSurface* fOpGridLArAboveSurface = new G4OpticalSurface("OpGridLArAboveSurface");
  new G4LogicalBorderSurface("GridLArAboveSurface", fPhysicLArAboveGrid, fPhysicGrid, fOpGridLArAboveSurface);
  fOpGridLArAboveSurface->SetType(dielectric_dielectric);
  fOpGridLArAboveSurface->SetModel(glisur);
  fOpGridLArAboveSurface->SetFinish(polished);
  
  G4double GridLArENE[4] = {0.1 * eV, 8.0 * eV, 8.3 * eV, 20.0 * eV};
  G4double GRUV = DSParameters::Get()->GetLArGridUVRef();
  G4double GRVIS = DSParameters::Get()->GetLArGridVisRef();
  G4double GridLArREF[4] = {GRVIS, GRVIS, GRUV, GRUV};
  G4MaterialPropertiesTable* fGridLArSurfProp = new G4MaterialPropertiesTable();
  if (DSParameters::Get()->GetWithNewGridModel())
    // the grid model is described in DSStorage.cc, the surface is treated in
    // G4OpBoundaryProcess.cc
    fGridLArSurfProp->AddConstProperty("DOGRID", 1, true);
  // Now use the following in old and new models.  By G4 convention,
  // "reflectivity" is actually 1-absorption.
  fGridLArSurfProp->AddProperty("REFLECTIVITY", GridLArENE, GridLArREF, 4);
  fOpGridLArSurface->SetMaterialPropertiesTable(fGridLArSurfProp);
  fOpGridLArAboveSurface->SetMaterialPropertiesTable(fGridLArSurfProp);
  
  // Grid->LAR
  if (DSParameters::Get()->GetWithNewGridModel()) {
    G4OpticalSurface* fOpLArGridSurface = new G4OpticalSurface("OpLArGridSurface");
    new G4LogicalBorderSurface("LArGridSurface", fPhysicGrid, fPhysicActiveLAr, fOpLArGridSurface);

    //Grid - LAr above grid
    G4OpticalSurface* fOpLArAboveGridSurface = new G4OpticalSurface("OpLArAboveGridSurface");
    new G4LogicalBorderSurface("LArAboveGridSurface", fPhysicGrid, fPhysicLArAboveGrid, fOpLArAboveGridSurface);
    //cout << " With DS50 new grid model " << endl;
    fOpLArGridSurface->SetType(dielectric_dielectric);
    fOpLArAboveGridSurface->SetType(dielectric_dielectric);
    //  fOpLArGridSurface->SetModel( glisur );
    //  fOpLArGridSurface->SetFinish( polished );
    G4MaterialPropertiesTable* fLArGridSurfProp = new G4MaterialPropertiesTable();
    // the grid model is described in DSStorage.cc, the surface is treated in
    // G4OpBoundaryProcess.cc
    fLArGridSurfProp->AddConstProperty("DOGRIDEXIT", 1, true);
    fOpLArGridSurface->SetMaterialPropertiesTable(fLArGridSurfProp);
    fOpLArAboveGridSurface->SetMaterialPropertiesTable(fLArGridSurfProp);
  }

  // turning off the surface properties -MP 2019/06/2
  /*
  ////////////////////////////////////////
  // TPB <--> GAr and TPB <--> LAr
  G4OpticalSurface *fOpTPBGArSurface = new G4OpticalSurface("OpTPBGArSurface");
  new G4LogicalBorderSurface("GArTPBSurface", fPhysicGasPocket, fPhysicTPBTop,
  fOpTPBGArSurface ); new G4LogicalBorderSurface("TPBGArSurface", fPhysicTPBTop,
  fPhysicGasPocket, fOpTPBGArSurface ); new
  G4LogicalBorderSurface("LArTPBSurfaceSide", fPhysicActiveLAr, fPhysicTPBSide,
  fOpTPBGArSurface ); new G4LogicalBorderSurface("TPBLArSurfaceSide",
  fPhysicTPBSide, fPhysicActiveLAr, fOpTPBGArSurface ); new
  G4LogicalBorderSurface("LArTPBSurfaceBot", fPhysicActiveLAr, fPhysicTPBBottom,
  fOpTPBGArSurface ); new G4LogicalBorderSurface("TPBLArSurfaceBot",
  fPhysicTPBBottom, fPhysicActiveLAr, fOpTPBGArSurface );
  //new G4LogicalBorderSurface("LArLayerTPBSurface", fPhysicLArLayer,
  fPhysicTPBTop, fOpTPBGArSurface );
  //new G4LogicalBorderSurface("TPBLArLayerSurface", fPhysicTPBTop,
  fPhysicLArLayer, fOpTPBGArSurface ); fOpTPBGArSurface->SetType(
  dielectric_dielectric ); fOpTPBGArSurface->SetModel( unified );
  fOpTPBGArSurface->SetFinish( ground );
  fOpTPBGArSurface->SetSigmaAlpha(0.3);

  G4double VISTRAN = DSParameters::Get()->GetArTPBVisTran();

  const G4int NUM = 4;
  G4double pp[NUM] = {0.1*eV, 8.0*eV, 8.3*eV, 20.0*eV};  // vis, vis, UV, UV
  G4double specularlobe[NUM] = {0., 0., 0., 0.};         //--
  G4double specularspike[NUM] = {0., 0., 0., 0.};        //----  gives all
  reflection to Lambertian lobe G4double backscatter[NUM] = {0., 0., 0., 0.};
  //-- G4double reflectivity[NUM] = {1.0, 1.0, 1.0, 1.0};     //  To set
  1-absorption G4double transmitivity[NUM] = {VISTRAN, VISTRAN, 1.0, 1.0};    //
  To set reflection vs. transmission, overridding Fresnel
                                                         //  For now, no angle
  dependence. G4MaterialPropertiesTable *fTPBGArSurfProp = new
  G4MaterialPropertiesTable(); fTPBGArSurfProp ->
  AddProperty("SPECULARLOBECONSTANT",pp,specularlobe,NUM); fTPBGArSurfProp ->
  AddProperty("SPECULARSPIKECONSTANT",pp,specularspike,NUM); fTPBGArSurfProp ->
  AddProperty("BACKSCATTERCONSTANT",pp,backscatter,NUM); fTPBGArSurfProp ->
  AddProperty("REFLECTIVITY",pp,reflectivity,NUM); fTPBGArSurfProp ->
  AddProperty("TRANSMITTANCE",pp,transmitivity,NUM); fTPBGArSurfProp ->
  AddConstProperty("DOArTPB",1); fOpTPBGArSurface->SetMaterialPropertiesTable(
  fTPBGArSurfProp );
  */

  ////////////////////////////////////////
  // ITO (from DS50, kept as reference) //
  ////////////////////////////////////////

  // G4MaterialPropertiesTable *fITOSurfProp = new G4MaterialPropertiesTable();
  // if ( DSParameters::Get()->GetWithITO()  )
  // fITOSurfProp->AddConstProperty("DOITO",1);

  // do you want to simulate the actual refletor - LAr - acrylic - TPB sandwich
  // ?
  bool is_complex_reflector = DSStorage::Get()->GetDS20kReflectorAlternate();

  ////////////////////////////////////////////////////////////
  // TPB --> Reflector (Reflector)
  //  Should be no Reflector --> TPB as surface is defined as dielectric_metal.
  ///////////////////////////////////////////////////////////

  G4double ReflectorTPBENE[4] = {0.1 * eV, 8.0 * eV, 8.3 * eV, 20.0 * eV};
  G4double TREFUV = DSParameters::Get()->GetTeflonTPBUVRef();
  G4double TREFVIS = DSParameters::Get()->GetTeflonTPBVisRef();
  G4double ReflectorTPBREF[4] = {TREFVIS, TREFVIS, TREFUV, TREFUV};


  G4OpticalSurface* fOpTPBReflectorSurface = new G4OpticalSurface("OpTBPReflectorSurface");
  // TPC reflector
  // LL fPhysicLayerPassiveLAr2 , fPhysicLayerESR --> fPhysicLayerESR, fPhysicLayerNSLAr
  if ( is_complex_reflector ) new G4LogicalBorderSurface("TPBReflectorSurface", fPhysicLayerESR, fPhysicLayerNSLAr, fOpTPBReflectorSurface );
  else
    // LL fPhysicTPBSide, fPhysicLayerESR --> fPhysicLayerESR, fPhysicLayerPassiveLAr2
    new G4LogicalBorderSurface("TPBReflectorSurface", fPhysicLayerESR, fPhysicLayerPassiveLAr2, fOpTPBReflectorSurface);

  ///////////////
  // TPB Skin --> anything  (Reflector)
  //  Will do bi-directional to prevent light from the veto to enter the TPC
  // new G4LogicalBorderSurface("TPBReflectorSurface",  myTPBReflectorSkin ,
  // fMotherVolume , fOpTPBReflectorSurface );
  ////////////////////////////////////////
  fOpTPBReflectorSurface->SetType(dielectric_metal);
  fOpTPBReflectorSurface->SetModel(unified);

  // PDM: though I can't see how, the following settings are giving the desired
  // Lambertian reflection
  // fOpTPBReflectorSurface->SetFinish(groundfrontpainted);
  fOpTPBReflectorSurface->SetFinish(polished);
  // fOpTPBReflectorSurface->SetSigmaAlpha(0.1);

  G4MaterialPropertiesTable* fTPBReflectorSurfProp = new G4MaterialPropertiesTable();
  fTPBReflectorSurfProp->AddProperty("REFLECTIVITY", ReflectorTPBENE, ReflectorTPBREF, 4);
  fOpTPBReflectorSurface->SetMaterialPropertiesTable(fTPBReflectorSurfProp);

  G4String myTPBLayers = DSStorage::Get()->GetDS20kTPBLayers();
  if (myTPBLayers[3] == '0') {  // No TPB coating the inside of the veto --> kill
                                // VUV light and reflect visible
    G4double ForbibEnteringTPB[4] = {TREFVIS, TREFVIS, 0, 0};
    G4OpticalSurface* fOpTPBSkinSurface = new G4OpticalSurface("OpTBPSkinSurface");
    new G4LogicalBorderSurface("TPBSkinSurface", fMotherVolume, myTPBReflectorSkin, fOpTPBSkinSurface);
    fOpTPBSkinSurface->SetType(dielectric_metal);
    fOpTPBSkinSurface->SetModel(unified);
    fOpTPBSkinSurface->SetFinish(polished);
    G4MaterialPropertiesTable* fTPBSkinSurfProp = new G4MaterialPropertiesTable();
    fTPBSkinSurfProp->AddProperty("REFLECTIVITY", ReflectorTPBENE, ForbibEnteringTPB, 4);
    fOpTPBSkinSurface->SetMaterialPropertiesTable(fTPBSkinSurfProp);
  } else {  // TPB coating the inside of the veto --> let VUV light to enter and
            // be absorbed by TPB, reflect visible
    // G4double ForbibEnteringTPB[4] = {TREFVIS, TREFVIS, 0, 0};
    G4OpticalSurface* fOpTPBSkinSurface = new G4OpticalSurface("OpTBPSkinSurface");
    new G4LogicalBorderSurface("TPBSkinSurface", myTPBReflectorSkin, fPhysicAcrylicTPCVessel, fOpTPBSkinSurface);
    new G4LogicalBorderSurface("TPBSkinSurface", myTPBReflectorSkin, fPhysicPassiveTop, fOpTPBSkinSurface);
    new G4LogicalBorderSurface("TPBSkinSurface", myTPBReflectorSkin, fPhysicPassiveBottom, fOpTPBSkinSurface);
    fOpTPBSkinSurface->SetType(dielectric_metal);
    fOpTPBSkinSurface->SetModel(unified);
    // fOpTPBSkinSurface->SetFinish(polished);
    fOpTPBSkinSurface->SetFinish(ground);
    fOpTPBSkinSurface->SetSigmaAlpha(0.1);
    G4MaterialPropertiesTable* fTPBSkinSurfProp = new G4MaterialPropertiesTable();
    fTPBSkinSurfProp->AddProperty("REFLECTIVITY", ReflectorTPBENE, ReflectorTPBREF, 4);
    fOpTPBSkinSurface->SetMaterialPropertiesTable(fTPBSkinSurfProp);
  }
  /*
  ////////////////////////////////////////
  // Clevios /////
  ////////////////////////////////////////
  //surface either absorbs or transmits
  const G4int clevSize = 8;
  G4double clevEne[clevSize]   = {100., 375., 400., 420., 500., 600., 700.,
  800.}; //wl in nm G4double clevRefl[clevSize]  = {0., 0., 0.9457, 0.9553,
  0.9434, 0.9186, 0.8998, 0.8932}; //set (T+R) vs. absorption, i.e. ~15% chance
  of absorption G4double clevTrans[clevSize]; //set T vs. R for(int ij=0;
  ij<clevSize; ij++) { clevTrans[ij] = 1.; //all transmission and no reflection
    clevEne[ij] = 1240./clevEne[ij]*eV;
  }
  G4OpticalSurface *fOpCleviosSurface = new
  G4OpticalSurface("OpCleviosSurface"); new
  G4LogicalBorderSurface("TPBBottomAcrylicTPCVesselSurface", fPhysicTPBBottom,
  fPhysicAcrylicTPCVessel, fOpCleviosSurface ); new
  G4LogicalBorderSurface("AcrylicTPCVesselTPBBottomSurface",
  fPhysicAcrylicTPCVessel, fPhysicTPBBottom, fOpCleviosSurface ); new
  G4LogicalBorderSurface("TPBTopAcrylicTPCVesselSurface", fPhysicTPBTop,
  fPhysicAcrylicTPCVessel, fOpCleviosSurface ); new
  G4LogicalBorderSurface("AcrylicTPCVesselTPBTopSurface",
  fPhysicAcrylicTPCVessel, fPhysicTPBTop, fOpCleviosSurface ); new
  G4LogicalBorderSurface("PassiveBotAcrylicTPCVesselSurface",
  fPhysicAcrylicTPCVessel, fPhysicPassiveBottom, fOpCleviosSurface ); new
  G4LogicalBorderSurface("AcrylicTPCVesselPassiveBotSurface",
  fPhysicPassiveBottom, fPhysicAcrylicTPCVessel, fOpCleviosSurface );
  fOpCleviosSurface->SetType(dielectric_dielectric);
  fOpCleviosSurface->SetModel(unified);
  fOpCleviosSurface->SetFinish(polished);
  G4MaterialPropertiesTable *fCleviosSurfProp = new G4MaterialPropertiesTable();
  fCleviosSurfProp -> AddConstProperty("DOCLEVIOS",1);
  fCleviosSurfProp->AddProperty("REFLECTIVITY",clevEne,clevRefl,clevSize);
  fCleviosSurfProp->AddProperty("TRANSMITTANCE",clevEne,clevTrans,clevSize);
  fOpCleviosSurface->SetMaterialPropertiesTable( fCleviosSurfProp );
*/

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Neutron veto optics
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////
  // Teflon reflectivity //
  // PA teflon is a placeholder material, to prevent light hitting active
  // silicon from behind
  /////////////////////////

  G4double TefREF[4] = {0.2, 0.2, 0.01, 0.01};
  G4double Tefpp[4]  = {0.1 * eV, 8.0 * eV, 8.3 * eV, 20.0 * eV};

  G4OpticalSurface* fOpPDUTeflonSurface = new G4OpticalSurface("OpPDUTeflonSurface");
  fOpPDUTeflonSurface->SetType(dielectric_metal);
  fOpPDUTeflonSurface->SetModel(glisur);
  fOpPDUTeflonSurface->SetFinish(ground);

  G4MaterialPropertiesTable* fOpvPDUTeflonSurfProp = new G4MaterialPropertiesTable();
  fOpvPDUTeflonSurfProp->AddProperty("REFLECTIVITY", Tefpp, TefREF, 4);
  fOpPDUTeflonSurface->SetMaterialPropertiesTable(fOpvPDUTeflonSurfProp);

  // new G4LogicalSkinSurface("SiPMfusedsilicaoncap", myLogicShapeSiPMOnCap, fOpLArTeflonSurface);
  // new G4LogicalSkinSurface("SiPMfusedsilica", myLogicShapeSiPM, fOpLArTeflonSurface);
  new G4LogicalSkinSurface("OpvPDUTeflonSurface", myLogicShapevPDU, fOpPDUTeflonSurface);
  new G4LogicalSkinSurface("OpvPDUTeflonSurface", myLogicShapevQuadrant, fOpPDUTeflonSurface);

  ///////////////////////
  // SiPM reflectivity //
  ///////////////////////

  G4double myEnergy, myTra;
  G4int dim;

  G4double pp[916], SiPMTRN[916], SiPMREF[916];
  dim = 0;

  ifstream fSiPM_reflectivity("../data/detector/SiPM_transmittance_7deg_AstroCent.dat");

  while (!fSiPM_reflectivity.eof()) {
    fSiPM_reflectivity >> myEnergy >> myTra;
    if (fSiPM_reflectivity.eof()) break;
    pp[dim] = myEnergy;
    SiPMTRN[dim] = 1 - myTra;
    //        SiPMTRN[dim] = SiPMTRN[dim]*(0.82*0.98); // 0.82 is the lower
    //        limit of the correction factor for the reflectivity change for
    //        bare silicon between air and LAr
    SiPMTRN[dim] = SiPMTRN[dim] * (0.91 * 0.98);  // 0.91 is the central value of the correction factor for the
                                                  // reflectivity change for bare silicon between air and LAr
    //        SiPMTRN[dim] = SiPMTRN[dim]*(1.*0.98);   // 1.   is the upper
    //        limit of the correction factor for the reflectivity change for
    //        bare silicon between air and LAr
    SiPMTRN[dim] = 1 - SiPMTRN[dim];
    SiPMREF[dim] = 1.;
    dim++;
  }
  fSiPM_reflectivity.close();

  G4OpticalSurface* fOpLArSiPMSurface = new G4OpticalSurface("OpLArSiPMSurface");
  fOpLArSiPMSurface->SetType(dielectric_dielectric);
  fOpLArSiPMSurface->SetModel(unified);
  fOpLArSiPMSurface->SetFinish(polished);
  G4MaterialPropertiesTable* fLArSiPMSurfProp2 = new G4MaterialPropertiesTable();
  fLArSiPMSurfProp2->AddProperty("REFLECTIVITY", pp, SiPMREF, dim);
  fLArSiPMSurfProp2->AddProperty("TRANSMITTANCE", pp, SiPMTRN, dim);
  fOpLArSiPMSurface->SetMaterialPropertiesTable(fLArSiPMSurfProp2);

  new G4LogicalSkinSurface("SiPMfusedsilicaonside",  myLogicShapevSiPMSilica, fOpLArSiPMSurface);

  ///////////////////////////////////////////////////////////////////////////////
  // PEN --> Reflector (Reflector)
  // Should be no Reflector --> TPB as surface is defined as dielectric_metal.
  ///////////////////////////////////////////////////////////////////////////////

  // G4double TREFUV  = DSParameters::Get()->GetTeflonTPBUVRef();
  // G4double TREFVIS = DSParameters::Get()->GetTeflonTPBVisRef();

  // Choose 21 to use reflectivity of PEN air coupled to ESR or choose 22 to use
  // pure ESR reflectivity

  int ESR_reflectivity = 1;
  switch (ESR_reflectivity) {

    case 1:
      {
      //////////////////////////////////
      // PEN + air + ESR reflectivity //
      //////////////////////////////////
      G4double TeflonTPBENE[50] =
      {0.1, 2.067,2.175,2.214,2.255,2.340,2.385,2.431,2.436,2.531,2.583,2.638,2.696,2.725,2.756,2.787,2.818,2.884,2.918,2.952,2.988,3.024,3.039,3.047,3.054,3.062,3.069,3.077,3.085,3.092,3.100,3.108,3.116,3.123,3.131,3.139,3.147,3.155,3.163,3.171,3.179,3.188,3.196,3.204,3.212,3.221,3.263,
      8.0, 8.3, 20.0};
      for (int i=0; i<50; i++) TeflonTPBENE[i] *= eV;

      G4double TeflonTPBREF[50] = {100.130,
      100.130,99.995,99.856,99.681,99.659,99.569,99.351,99.306,99.018,98.652,98.415,98.283,98.018,97.856,97.606,97.457,97.134,96.928,96.827,96.247,95.737,95.359,95.197,95.048,94.876,94.684,94.463,94.055,93.650,93.147,92.562,91.812,90.904,89.807,88.506,86.957,85.242,83.156,80.678,77.811,74.615,71.004,67.089,62.924,58.670,20.000,20.000,
      TREFUV , TREFUV};
      for (int i=0; i<48; i++) TeflonTPBREF[i] *= TREFVIS/100.;

      G4OpticalSurface* fOpPENReflectorSurface = new G4OpticalSurface("OpTPBReflectorSurface");
      new G4LogicalBorderSurface("TPBReflectorSurface",  myPENReflectorSkin,  myNSLArReflectorSkin , fOpPENReflectorSurface );
      fOpPENReflectorSurface->SetType( dielectric_metal );
      fOpPENReflectorSurface->SetModel(unified);
      fOpPENReflectorSurface->SetFinish(polished);

      G4MaterialPropertiesTable *fOpTPBReflectorSurfProp = new G4MaterialPropertiesTable();
      fOpTPBReflectorSurfProp->AddProperty("REFLECTIVITY", TeflonTPBENE, TeflonTPBREF, 50); // Should be 50 for the arrays written above
      fOpPENReflectorSurface->SetMaterialPropertiesTable(fOpTPBReflectorSurfProp );
      }
    break;

    case 2:
      {
      ///////////////////////////
      // Pure ESR reflectivity //
      ///////////////////////////

      G4double myvalue, myene;
      G4double TeflonTPBENE[831], TeflonTPBREF[831];
      dim = 0;

      // ifstream
      // fESR_reflectivity("../data/detector/ESR_reflectivity_2_AstroCent.dat");
      //
      // if ( !fESR_reflectivity.is_open())
      //   DSLog(fatal) << "ERROR: Could not open ESR reflectivity file" <<endlog;
      while (DSIO::Get()->GetStreamDS20kESRreflectivity() >> myene >> myvalue) {
        // fESR_reflectivity >> myene >> myvalue ;
        // if(fESR_reflectivity.eof()) break;
        TeflonTPBENE[dim] = myene * eV;
        TeflonTPBREF[dim] = myvalue;
        // cout << "Energy: " << TeflonTPBENE[dim] << " Ref: " << TeflonTPBREF[dim]
        // << endl;

        // if(TeflonTPBENE[dim] > 3.27 || TeflonTPBENE[dim] < 2.05) {
        //       TeflonTPBREF[dim] = 0 ;
        // }
        // else{
        //       TeflonTPBREF[dim] = myvalue;
        // }
        dim++;
      }

      // fESR_reflectivity.close();

      for (int i = 0; i < dim; i++) TeflonTPBREF[i] *= TREFVIS / 100.;
      G4OpticalSurface* fOpPENReflectorSurface = new G4OpticalSurface("OpTPBReflectorSurface");
      new G4LogicalBorderSurface("TPBReflectorSurface", myNSLArReflectorSkin, myGdReflectorSkin, fOpPENReflectorSurface);

      fOpPENReflectorSurface->SetType(dielectric_metal);
      fOpPENReflectorSurface->SetModel(unified);
      fOpPENReflectorSurface->SetFinish(polished);
      G4MaterialPropertiesTable* fOpTPBReflectorSurfProp = new G4MaterialPropertiesTable();
      fOpTPBReflectorSurfProp->AddProperty("REFLECTIVITY", TeflonTPBENE, TeflonTPBREF,
                                           dim);  // Should be 50 for the arrays written above
      fOpPENReflectorSurface->SetMaterialPropertiesTable(fOpTPBReflectorSurfProp);
      }
    break;
  }//switch ends here --
}


// Use these functions to determine which octant a given (X,Y) coordinate is
double DSDetectorDS20k::Angle2D(double x1, double y1, double x2, double y2) {
  double dtheta, theta1, theta2;

  theta1 = atan2(y1, x1);
  theta2 = atan2(y2, x2);
  dtheta = theta2 - theta1;

  while (dtheta > M_PI) dtheta -= 2. * M_PI;
  while (dtheta < -M_PI) dtheta += 2. * M_PI;

  return (dtheta);
}

bool DSDetectorDS20k::InsideOctagon(float _radius, double pxx, double pyy)  //_radius is the circumscribed one
{

  double squareh[8], squarev[8];

  double sq = sqrt(2) / 2.;
  squareh[0] = 1;
  squarev[0] = 0;
  squareh[1] = sq;
  squarev[1] = sq;
  squareh[2] = 0;
  squarev[2] = 1;
  squareh[3] = -sq;
  squarev[3] = sq;
  squareh[4] = -1;
  squarev[4] = 0;
  squareh[5] = -sq;
  squarev[5] = -sq;
  squareh[6] = 0;
  squarev[6] = -1;
  squareh[7] = sq;
  squarev[7] = -sq;

  double angle = 0.;
  double p1h, p1v, p2h, p2v;

  double px = pxx * cos(M_PI / 8.) - pyy * sin(M_PI / 8);
  double py = pxx * sin(M_PI / 8.) + pyy * cos(M_PI / 8);

  for (int _i = 0; _i < 8; ++_i) {
    p1h = _radius * squareh[_i] - px;
    p1v = _radius * squarev[_i] - py;
    p2h = _radius * squareh[(_i + 1) % 8] - px;
    p2v = _radius * squarev[(_i + 1) % 8] - py;

    angle += DSDetectorDS20k::Angle2D(p1h, p1v, p2h, p2v);
  }
  if (abs(angle + M_PI / 19) < pi) return false;
  else
    return true;
}

/*
int DSDetectorDS20k::get_sector(double veto_x, double veto_y) {
  double ang;
  if (veto_x > 0) {
    if (veto_y >= 0) {
      ang = TMath::ATan((veto_y) / (veto_x));
    } else {
      ang = 2. * TMath::Pi() - abs(TMath::ATan((veto_y) / (veto_x)));
    }
  } else if (veto_x < 0) {
    ang = TMath::Pi() - TMath::ATan((veto_y) / (-veto_x));
  } else {
    if (veto_x == 0 && veto_y > 0) {
      ang = TMath::Pi() / 2.;
    } else if (veto_x == 0 && veto_y < 0) {
      ang = 3. * TMath::Pi() / 2.;
    } else {
      ang = 0;
    }
  }
  int half_sector = trunc(ang / (TMath::Pi() / 8));
  int sector = round(half_sector / 2.);
  if (half_sector == 15) { sector = 0; }
  return sector;
}
*/
double DSDetectorDS20k::GetOctagonInnerRadius(double edge) {
  return edge / 2 * (1 + sqrt(2));
}
double DSDetectorDS20k::GetOctagonOuterRadius(double edge) {
  return edge / 2. * sqrt(4 + 2 * sqrt(2));
}
