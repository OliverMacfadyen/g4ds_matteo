#include "DSDetectorPlasticVeto.hh"

#include "DSDetectorPMTNeutronVeto.hh"
#include "DSEventHandler.hh"
#include "DSIO.hh"
#include "DSLogger.hh"
#include "DSMaterial.hh"
#include "DSParameters.hh"
#include "DSStorage.hh"

#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4LogicalVolume.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4NistManager.hh"
#include "G4Orb.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4Polyhedra.hh"
#include "G4RunManager.hh"
#include "G4Sphere.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4UnionSolid.hh"
#include "G4VisAttributes.hh"
#include "G4VisExecutive.hh"

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

using namespace std;

DSDetectorPlasticVeto::DSDetectorPlasticVeto(G4VPhysicalVolume* myMotherVolume) {

  fMotherVolume = myMotherVolume;

  // Color definition for volumes

  const G4Colour myRed(1.0, 0.0, 0.0);
  const G4VisAttributes VisualRed = myRed;

  const G4Colour myBlue(0.0, 0.0, 1.0);
  const G4VisAttributes VisualBlue = myBlue;

  const G4Colour myGreen(0.0, 1.0, 0.0);
  const G4VisAttributes VisualGreen = myGreen;

  const G4Colour myYellow(1.0, 1.0, 0.0);
  const G4VisAttributes VisualYellow = myYellow;

  const G4Colour myPurple(1.0, 0.0, 1.0);
  const G4VisAttributes VisualPurple = myPurple;

  const G4Colour myCyan(0.0, 1.0, 1.0);
  const G4VisAttributes VisualCyan = myCyan;

  const G4Colour myBlack(0.0, 0.0, 0.0);
  const G4VisAttributes VisualBlack = myBlack;

  // Color definition for string in cout

  const std::string red("\033[0;31m");
  const std::string green("\033[1;32m");
  const std::string yellow("\033[1;33m");
  const std::string cyan("\033[0;36m");
  const std::string magenta("\033[0;35m");
  const std::string reset("\033[0m");

  // Mother volume setting
  if (myMotherVolume->GetLogicalVolume()->GetMaterial()->GetIndex() == 3) {  // check if the system is placed in side LAr or air

    G4Box* _fSolidWorld = new G4Box("World_Solid", 5 * m, 5 * m, 5 * m);

    G4LogicalVolume* _fLogicWorld = new G4LogicalVolume(_fSolidWorld, DSMaterial::Get()->GetOVLiquidArgon(), "World_Logic");
    G4Colour myWhite(1.0, 1.0, 1.0);  // white
    G4VisAttributes* LogVisAttWorld = new G4VisAttributes(myWhite);
    LogVisAttWorld->SetVisibility(false);
    _fLogicWorld->SetVisAttributes(LogVisAttWorld);

    G4VPhysicalVolume* LArVolume = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), "World", _fLogicWorld, myMotherVolume, false, 0);
    fMotherVolume = LArVolume;
  }

  // Overlaps checking
  G4bool myCheckOverlaps = DSStorage::Get()->GetCheckOverlap();

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////   MATERIALS SETTING
  ////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////

  // TPB
  G4Material* TPB_mat = DSMaterial::Get()->GetTPB();

  // Acrylic
  G4Material* acrylic_mat = DSMaterial::Get()->GetAcrylic();

  // Plastic
  G4Material* plastic_mat = DSMaterial::Get()->GetGdAcrylic();

  // LAr
  G4Material* PScintVetoLiquidArgon = DSMaterial::Get()->GetPScintVetoLiquidArgon();
  G4Material* VetoLiquidArgon = DSMaterial::Get()->GetOVLiquidArgon();
  G4Material* NSLiquidArgon = DSMaterial::Get()->GetNSLiquidArgon();

  // Copper
  G4Material* copper_mat = DSMaterial::Get()->GetMetalCopper();

  // Teflon
  G4Material* teflon_mat = DSMaterial::Get()->GetTeflon();

  // Metal silicon
  G4Material* metalsilicon_mat = DSMaterial::Get()->GetMetalSilicon();

  // Fused silica
  G4Material* fusedsilica_mat = DSMaterial::Get()->GetFusedSilica();

  /////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////   OFFSETS
  ///////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  // Offset of TPB panels on faces from inner plastic TPB
  const G4double internalOffset = DSStorage::Get()->GetDS20kinternalOffset();

  // Same offset but for TPB panels on edges
  const G4double intOff = internalOffset / (sin(67.5 * deg));  // From geom. considerations

  // Offset of TPB panels from outer plastic TPB
  const G4double externalOffset = DSStorage::Get()->GetDS20kexternalOffset();
  G4ThreeVector ExternalOffset = G4ThreeVector(externalOffset, 0, 0);

  // Same offset but for TPB panels on edges
  const G4double extOff = externalOffset / (sin(67.5 * deg));  // From geom. considerations
  G4ThreeVector ExtOff = G4ThreeVector(extOff, 0, 0);

  // Offset of vertical inner TPB panels from TPC
  const G4double verticalInnerPanelOffset = DSStorage::Get()->GetDS20kverticalInnerPanelOffset();
  ;
  G4ThreeVector VerticalInnerPanelOffset = G4ThreeVector(verticalInnerPanelOffset, 0, 0);

  // Same offset but for TPB panels on edges
  const G4double verticalInnerPanelOff = verticalInnerPanelOffset / (sin(67.5 * deg));
  G4ThreeVector VerticalInnerPanelOff = G4ThreeVector(verticalInnerPanelOff, 0, 0);

  // Offset of vertical outer TPB panels from copper TPB
  const G4double verticalOuterPanelOffset = DSStorage::Get()->GetDS20kverticalOuterPanelOffset();

  // Same offset but for TPB panels on edges
  const G4double verticalOuterPanelOff = verticalOuterPanelOffset / (sin(67.5 * deg));  // From geom. considerations

  // Offset for cilindrical inner volume on z axis
  // If too small, inner horizontal panels will overlap at the center!
  const G4double centerInnerOffset = DSStorage::Get()->GetDS20kcenterInnerOffset();
  G4ThreeVector CenterInnerOffset = G4ThreeVector(centerInnerOffset, 0, 0);

  // Offset for cilindrical outer volume on z axis
  // If too small, outer horizontal panels will overlap at the center!
  const G4double centerOuterOffset = DSStorage::Get()->GetDS20kcenterOuterOffset();
  G4ThreeVector CenterOuterOffset = G4ThreeVector(centerOuterOffset, 0, 0);

  // Upper horizontal inner TPB panel offset from electronics or, in general,
  // whatever volume is placed under them
  const G4double upperHorizontalInnerPanelOffset = DSStorage::Get()->GetDS20kupperHorizontalInnerPanelOffset();
  G4ThreeVector UpperHorizontalInnerPanelOffset = G4ThreeVector(0, 0, upperHorizontalInnerPanelOffset);

  // Lower horizontal inner TPB panel offset from electronics or, in general,
  // whatever volume is placed under them
  const G4double lowerHorizontalInnerPanelOffset = DSStorage::Get()->GetDS20klowerHorizontalInnerPanelOffset();
  G4ThreeVector LowerHorizontalInnerPanelOffset = G4ThreeVector(0, 0, lowerHorizontalInnerPanelOffset);

  // Upper horizontal outer TPB panel offset from copper TPB
  const G4double upperHorizontalOuterPanelOffset = DSStorage::Get()->GetDS20kupperHorizontalOuterPanelOffset();
  G4ThreeVector UpperHorizontalOuterPanelOffset = G4ThreeVector(0, 0, upperHorizontalOuterPanelOffset);

  // Lower horizontal outer TPB panel offset from copper TPB
  const G4double lowerHorizontalOuterPanelOffset = DSStorage::Get()->GetDS20klowerHorizontalOuterPanelOffset();
  G4ThreeVector LowerHorizontalOuterPanelOffset = G4ThreeVector(0, 0, lowerHorizontalOuterPanelOffset);

  ///////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////   PARAMETERS FROM DRAWINGS
  ////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////

  const G4double myTPCHeight = DSStorage::Get()->GetDS20kTPCheight();
  const G4double myTPCEdge = DSStorage::Get()->GetDS20kTPCedge();  // Internal edge without thickness
  const G4double myGasPocketThickness = DSStorage::Get()->GetDS20kGasPocketThickness();
  const G4double myElectronicsSpace = DSStorage::Get()->GetDS20kLArThicknessAboveTPC();
  const G4double myPlasticThickness = DSStorage::Get()->GetDS20kGdPlasticThickness();  // Active thickness wanted. In the
                                                                                       // middle there will be placed a layer
                                                                                       // of LAr

  const G4double myLArBufferThickness = DSStorage::Get()->GetDS20kLArBufferThickness();  // All the volume between the TPC and
                                                                                         // inner plastic TPB
                                                                                         // //All the volume between outer
                                                                                         // plastic TPB and copper TPB

  const G4double myTPBThickness = DSStorage::Get()->GetDS20kTPBThickness();
  const G4double myTPCVesselThickness = DSStorage::Get()->GetDS20kTPCVesselThickness();

  // TPB panels thickness
  const G4double myVerticalInnerTPBPanely = DSStorage::Get()->GetDS20kVerticalInnerTPBPanely();
  const G4double myHorizontalInnerTPBPanely = DSStorage::Get()->GetDS20kHorizontalInnerTPBPanely();
  const G4double myVerticalOuterTPBPanely = DSStorage::Get()->GetDS20kVerticalOuterTPBPanely();
  const G4double myHorizontalOuterTPBPanely = DSStorage::Get()->GetDS20kHorizontalOuterTPBPanely();

  //////////////  The following variables are referred to all the volume built
  ///in the TPC class  ///////////////////

  // Height of the entire volume built in TPC class. If one changes the
  // structure of the TPC, this variable must be manually changed.
  const G4double effectiveDS20kTPCheight = myTPCHeight + myGasPocketThickness + 2. * myTPCVesselThickness + 2. * myTPBThickness + 2. * myElectronicsSpace;

  const G4double inscribedTPCThick = myTPCEdge * ((1. + sqrt(2.)) / 2.) + myTPCVesselThickness + myTPBThickness;  // From geom. considerations

  const G4double apothemTPCThick = inscribedTPCThick;
  const G4ThreeVector aTPCThick = G4ThreeVector(apothemTPCThick, 0, 0);

  const G4double radiusTPCThick = apothemTPCThick / (((1. + sqrt(2.)) / 2.) * sqrt(2. - sqrt(2.)));  // From geom. considerations
  const G4ThreeVector rTPCThick = G4ThreeVector(radiusTPCThick, 0, 0);

  ///////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////   GLOBAL VARIABLES
  /////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////

  const G4double ang = 2 * (90. - 67.5) * deg;  // Default rotation angle

  cout << " " << endl;
  DSLog(routine) << red << "Constructing plastic veto geometry " << reset << endlog;
  cout << " " << endl;

  //######################################################################################
  //######################################################################################
  //######################################################################################

  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  // x axis is set to be orthogonal to a face of all octagonal solids
  //////////////////////////////////////////////////////////////////////
  // Only the most external mother volume is rotated, so that all the
  // daughters will be rotated too///////////////////////////////////////
  //////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////
  ///////////////////INNER TPB PANELS/////////////////////////////
  ////////////////////////////////////////////////////////////////

  /////////////////////////
  ///////Vertical TPB panels/////////////////////////////////////
  /////////////////////////

  // Vertical TPB panels on faces
  const G4double verticalInnerTPBPanelFacex = myLArBufferThickness - verticalInnerPanelOffset - internalOffset;  // Need to subtract panel gaps
  // Vertical TPB panels on edges
  const G4double verticalInnerTPBPanelEdgex = (verticalInnerTPBPanelFacex + verticalInnerPanelOffset + internalOffset) / (sin(67.5 * deg)) - myVerticalInnerTPBPanely * 0.5 * (1. / (tan(67.5 * deg))) - verticalInnerPanelOff - intOff;  // Panels on edges must be wider to fit the octagonal shape.
                                                                                                                                                                                                                                          // Calculus from geom. considerations

  const G4double verticalInnerTPBPanelz = effectiveDS20kTPCheight + 2. * myLArBufferThickness - 2. * internalOffset;  // 2.* for upper and lower spacing

  // Useful positions for vertical TPB panels on faces
  G4ThreeVector verticalInnerTPBPanelFaceX = G4ThreeVector(verticalInnerTPBPanelFacex / 2., 0, 0);  // Position of panel geom. center

  // Useful positions for vertical TPB panels on edges
  G4ThreeVector verticalInnerTPBPanelEdgeX = G4ThreeVector(verticalInnerTPBPanelEdgex / 2., 0, 0);  // Position of panel geom. center

  // Useful position for vertical TPB panels
  G4ThreeVector verticalInnerTPBPanelZ = G4ThreeVector(0, 0,
                                                       verticalInnerTPBPanelz / 2.);  // Vertical position of vertical inner panel geom. center

  // Solid for vertical TPB panels on faces
  G4Box* myVerticalInnerTPBPanelFace = new G4Box("verticalInnerTPBPanelFace", verticalInnerTPBPanelFacex / 2., myVerticalInnerTPBPanely / 2., verticalInnerTPBPanelz / 2.);
  // Solid for vertical TPB panels on edges
  G4Box* myVerticalInnerTPBPanelEdge = new G4Box("verticalInnerTPBPanelEdge", verticalInnerTPBPanelEdgex / 2., myVerticalInnerTPBPanely / 2., verticalInnerTPBPanelz / 2.);

  /////////////////////////
  ///////Horizontal TPB panels///////////////////////////////////////
  //////////////////////////

  // Horizontal TPB panels on faces
  const G4double horizontalInnerTPBPanelFacex = verticalInnerTPBPanelFacex + verticalInnerPanelOffset + apothemTPCThick - centerInnerOffset;
  // Horizontal panel on edges
  const G4double horizontalInnerTPBPanelEdgex = (horizontalInnerTPBPanelFacex + centerInnerOffset) / (((1. + sqrt(2.)) / 2.) * sqrt(2. - sqrt(2.))) - myHorizontalInnerTPBPanely * 0.5 * (1. / (tan(67.5 * deg))) - centerInnerOffset;  // Edges TPB panels must be wider to fit the octagonal
                                                                                                                                                                                                                                        // shape. Calculus from geom. considerations

  const G4double upperHorizontalInnerTPBPanelz = myLArBufferThickness - upperHorizontalInnerPanelOffset - internalOffset;  // Horizontal panels are vertically shorter
  const G4double lowerHorizontalInnerTPBPanelz = myLArBufferThickness - lowerHorizontalInnerPanelOffset - internalOffset;  // to leave space from electronics

  // Useful positions for upper horizontal TPB panels
  G4ThreeVector upperHorizontalInnerTPBPanelZ = G4ThreeVector(0, 0,
                                                              upperHorizontalInnerTPBPanelz / 2.);  // Position of panel geom. center

  // Useful positions for lower horizontal TPB panels
  G4ThreeVector lowerHorizontalInnerTPBPanelZ = G4ThreeVector(0, 0,
                                                              lowerHorizontalInnerTPBPanelz / 2.);  // Position of panel geom. center

  // Useful positions for TPB panels on faces
  G4ThreeVector horizontalInnerTPBPanelFaceX = G4ThreeVector(horizontalInnerTPBPanelFacex / 2., 0,
                                                             0);  // Position of panel geom. center

  // Useful positions for TPB panels on edges
  G4ThreeVector horizontalInnerTPBPanelEdgeX = G4ThreeVector(horizontalInnerTPBPanelEdgex / 2., 0,
                                                             0);  // Position of panel geom. center

  // Solid for upper horizontal TPB panels on faces
  G4Box* myUpperHorizontalInnerTPBPanelFace = new G4Box("upperHorizontalInnerTPBPanelFace", horizontalInnerTPBPanelFacex / 2., myHorizontalInnerTPBPanely / 2., upperHorizontalInnerTPBPanelz / 2.);

  // Solid for lower horizontal TPB panels on faces
  G4Box* myLowerHorizontalInnerTPBPanelFace = new G4Box("lowerHorizontalInnerTPBPanelFace", horizontalInnerTPBPanelFacex / 2., myHorizontalInnerTPBPanely / 2., lowerHorizontalInnerTPBPanelz / 2.);

  // Solid for upper horizontal TPB panels on edges
  G4Box* myUpperHorizontalInnerTPBPanelEdge = new G4Box("upperHorizontalInnerTPBPanelEdge", horizontalInnerTPBPanelEdgex / 2., myHorizontalInnerTPBPanely / 2., upperHorizontalInnerTPBPanelz / 2.);

  // Solid for lower horizontal TPB panels on edges
  G4Box* myLowerHorizontalInnerTPBPanelEdge = new G4Box("lowerHorizontalInnerTPBPanelEdge", horizontalInnerTPBPanelEdgex / 2., myHorizontalInnerTPBPanely / 2., lowerHorizontalInnerTPBPanelz / 2.);

  /////////////////////////////////////////////////////////////
  /////////Union solid and placement for inner TPB panels//////
  /////////////////////////////////////////////////////////////

  //
  /////////TPB panels on faces//////////////////////////////////////////
  //

  // Position of horizontal panel respect to vertical panel
  G4ThreeVector traslHorizUpperInnerTPBPanelFace = (verticalInnerTPBPanelZ - upperHorizontalInnerTPBPanelZ) - (horizontalInnerTPBPanelFaceX - verticalInnerTPBPanelFaceX);
  G4ThreeVector traslHorizLowerInnerTPBPanelFace = -(verticalInnerTPBPanelZ - lowerHorizontalInnerTPBPanelZ) - (horizontalInnerTPBPanelFaceX - verticalInnerTPBPanelFaceX);

  // Vertical TPB panel + upper horizontal TPB panel

  G4UnionSolid* tmpInnerTPBPanelFace = new G4UnionSolid("tmpInnerTPBPanelFace", myVerticalInnerTPBPanelFace, myUpperHorizontalInnerTPBPanelFace, 0, traslHorizUpperInnerTPBPanelFace);

  // Vertical TPB panel + upper horizontal TPB panel + lower horizontal TPB
  // panel
  G4UnionSolid* myInnerTPBPanelFace = new G4UnionSolid("innerTPBPanelFace", tmpInnerTPBPanelFace, myLowerHorizontalInnerTPBPanelFace, 0, traslHorizLowerInnerTPBPanelFace);

  G4LogicalVolume* myLogicShapeInnerTPBPanelFace[8];

  for (int i = 0; i < 8; i++) {

    myLogicShapeInnerTPBPanelFace[i] = new G4LogicalVolume(myInnerTPBPanelFace, TPB_mat, "innerTPBPanelFace");

    myLogicShapeInnerTPBPanelFace[i]->SetVisAttributes(VisualRed);
  }

  // Placement at the end of the file

  //
  /////////TPB panels on edges//////////////////////////////////////////
  //

  // Position of horizontal panel respect to vertical panel
  G4ThreeVector traslHorizUpperInnerTPBPanelEdge = (verticalInnerTPBPanelZ - upperHorizontalInnerTPBPanelZ) - (horizontalInnerTPBPanelEdgeX - verticalInnerTPBPanelEdgeX);
  G4ThreeVector traslHorizLowerInnerTPBPanelEdge = -(verticalInnerTPBPanelZ - lowerHorizontalInnerTPBPanelZ) - (horizontalInnerTPBPanelEdgeX - verticalInnerTPBPanelEdgeX);

  // Vertical TPB panel + upper horizontal TPB panel
  G4UnionSolid* tmpInnerTPBPanelEdge = new G4UnionSolid("tmpInnerTPBPanelEdge", myVerticalInnerTPBPanelEdge, myUpperHorizontalInnerTPBPanelEdge, 0, traslHorizUpperInnerTPBPanelEdge);

  // Vertical TPB panel + upper horizontal TPB panel + lower horizontal TPB
  // panel

  G4UnionSolid* myInnerTPBPanelEdge = new G4UnionSolid("InnerTPBPanelEdge", tmpInnerTPBPanelEdge, myLowerHorizontalInnerTPBPanelEdge, 0, traslHorizLowerInnerTPBPanelEdge);

  G4LogicalVolume* myLogicShapeInnerTPBPanelEdge[8];

  for (int i = 0; i < 8; i++) {

    myLogicShapeInnerTPBPanelEdge[i] = new G4LogicalVolume(myInnerTPBPanelEdge, TPB_mat, "innerTPBPanelEdge");

    myLogicShapeInnerTPBPanelEdge[i]->SetVisAttributes(VisualRed);
  }

  // Placement at the end of the file

  //////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////INNER ARYLIC PANELS////////////////////////////
  //////////////////////////////////////////////////////////////////////////////

  //
  //////////////////Vertical acrylic panels///////////////
  //

  // Dimensions defined from vertical TPB panels dimensions
  const G4double verticalAcrylicInnerPanely = myVerticalInnerTPBPanely - 2. * myTPBThickness;  // Panel thickness
  // Vertical acrylic panels on faces
  const G4double verticalAcrylicInnerPanelFacex = verticalInnerTPBPanelFacex - 2. * myTPBThickness;
  // Vertical acrylic panels on edges
  const G4double verticalAcrylicInnerPanelEdgex = verticalInnerTPBPanelEdgex - 2. * myTPBThickness;
  // Vertical acrylic panels height
  const G4double verticalAcrylicInnerPanelz = verticalInnerTPBPanelz - 2. * myTPBThickness;

  // Solid for vertical acrylic panels on faces
  G4Box* myVerticalAcrylicInnerPanelFace = new G4Box("verticalAcrylicInnerPanelFace", verticalAcrylicInnerPanelFacex / 2., verticalAcrylicInnerPanely / 2., verticalAcrylicInnerPanelz / 2.);

  // Solid for vertical acrylic panels on edges
  G4Box* myVerticalAcrylicInnerPanelEdge = new G4Box("verticalAcrylicInnerPanelEdge", verticalAcrylicInnerPanelEdgex / 2., verticalAcrylicInnerPanely / 2., verticalAcrylicInnerPanelz / 2.);

  // Useful positions for vertical acrylic panels on faces
  G4ThreeVector verticalAcrylicInnerPanelFaceX = G4ThreeVector(verticalAcrylicInnerPanelFacex / 2., 0,
                                                               0);  // Position of panel geom. center

  // Useful positions for vertical acrylic panels on edges
  G4ThreeVector verticalAcrylicInnerPanelEdgeX = G4ThreeVector(verticalAcrylicInnerPanelEdgex / 2., 0,
                                                               0);  // Position of panel geom. center

  // Useful position for vertical acrylic panels
  G4ThreeVector verticalAcrylicInnerPanelZ = G4ThreeVector(0, 0,
                                                           verticalAcrylicInnerPanelz / 2.);  // Vertical position of vertical inner panel geom. center

  //
  ////////////////Horizontal panels////////////////////
  //

  // Dimensions defined from horizontal TPB panels dimensions
  const G4double horizontalAcrylicInnerPanely = myHorizontalInnerTPBPanely - 2. * myTPBThickness;  // Panel thickness
  // Horizontal panels on faces
  const G4double horizontalAcrylicInnerPanelFacex = horizontalInnerTPBPanelFacex - 2. * myTPBThickness;
  // Horizontal panel on edges
  const G4double horizontalAcrylicInnerPanelEdgex = horizontalInnerTPBPanelEdgex - 2. * myTPBThickness;

  // Horizontal panels height
  const G4double upperHorizontalAcrylicInnerPanelz = upperHorizontalInnerTPBPanelz - 2. * myTPBThickness;
  const G4double lowerHorizontalAcrylicInnerPanelz = lowerHorizontalInnerTPBPanelz - 2. * myTPBThickness;

  // Solid for upper horizontal acrylic panels on faces
  G4Box* myUpperHorizontalAcrylicInnerPanelFace = new G4Box("upperHorizontalAcrylicInnerPanelFace", horizontalAcrylicInnerPanelFacex / 2., horizontalAcrylicInnerPanely / 2., upperHorizontalAcrylicInnerPanelz / 2.);

  // Solid for lower horizontal acrylic panels on faces
  G4Box* myLowerHorizontalAcrylicInnerPanelFace = new G4Box("lowerHorizontalAcrylicInnerPanelFace", horizontalAcrylicInnerPanelFacex / 2., horizontalAcrylicInnerPanely / 2., lowerHorizontalAcrylicInnerPanelz / 2.);

  // Solid for upper horizontal acrylic panels on edges
  G4Box* myUpperHorizontalAcrylicInnerPanelEdge = new G4Box("upperHorizontalAcrylicInnerPanelEdge", horizontalAcrylicInnerPanelEdgex / 2., horizontalAcrylicInnerPanely / 2., upperHorizontalAcrylicInnerPanelz / 2.);

  // Solid for lower horizontal acrylic panels on edges
  G4Box* myLowerHorizontalAcrylicInnerPanelEdge = new G4Box("lowerHorizontalAcrylicInnerPanelEdge", horizontalAcrylicInnerPanelEdgex / 2., horizontalAcrylicInnerPanely / 2., lowerHorizontalAcrylicInnerPanelz / 2.);

  // Useful positions for upper horizontal acrylic panels
  G4ThreeVector upperHorizontalAcrylicInnerPanelZ = G4ThreeVector(0, 0,
                                                                  upperHorizontalAcrylicInnerPanelz / 2.);  // Position of panel geom. center

  // Useful positions for lower horizontal acrylic panels
  G4ThreeVector lowerHorizontalAcrylicInnerPanelZ = G4ThreeVector(0, 0,
                                                                  lowerHorizontalAcrylicInnerPanelz / 2.);  // Position of panel geom. center

  // Useful positions for acrylic panels on faces
  G4ThreeVector horizontalAcrylicInnerPanelFaceX = G4ThreeVector(horizontalAcrylicInnerPanelFacex / 2., 0,
                                                                 0);  // Position of panel geom. center

  // Useful positions for acrylic panels on edges
  G4ThreeVector horizontalAcrylicInnerPanelEdgeX = G4ThreeVector(horizontalAcrylicInnerPanelEdgex / 2., 0,
                                                                 0);  // Position of panel geom. center

  /////////////////////////////////////////////////////////////
  /////Union solid and placement for inner acrylic panels//////
  /////////////////////////////////////////////////////////////

  //
  /////////Acrylic panels on faces//////////////////////////////////////////
  //

  // Position of horizontal panel respect to vertical panel
  G4ThreeVector traslHorizUpperAcrylicInnerPanelFace = (verticalAcrylicInnerPanelZ - upperHorizontalAcrylicInnerPanelZ) - (horizontalAcrylicInnerPanelFaceX - verticalAcrylicInnerPanelFaceX);
  G4ThreeVector traslHorizLowerAcrylicInnerPanelFace = -(verticalAcrylicInnerPanelZ - lowerHorizontalAcrylicInnerPanelZ) - (horizontalAcrylicInnerPanelFaceX - verticalAcrylicInnerPanelFaceX);

  // Vertical acrylic panel + upper horizontal acrylic panel
  G4UnionSolid* tmpAcrylicInnerPanelFace = new G4UnionSolid("tmpAcrylicInnerPanelFace", myVerticalAcrylicInnerPanelFace, myUpperHorizontalAcrylicInnerPanelFace, 0, traslHorizUpperAcrylicInnerPanelFace);

  // Vertical acrylic panel + upper horizontal acrylic panel + lower horizontal
  // acrylic panel
  G4UnionSolid* myAcrylicInnerPanelFace = new G4UnionSolid("innerAcrylicPanelFace", tmpAcrylicInnerPanelFace, myLowerHorizontalAcrylicInnerPanelFace, 0, traslHorizLowerAcrylicInnerPanelFace);

  G4LogicalVolume* myLogicShapeAcrylicInnerPanelFace[8];

  for (int i = 0; i < 8; i++) {

    myLogicShapeAcrylicInnerPanelFace[i] = new G4LogicalVolume(myAcrylicInnerPanelFace, acrylic_mat, "acrylicInnerPanelFace");

    myLogicShapeAcrylicInnerPanelFace[i]->SetVisAttributes(VisualBlue);
  }

  // Placement at the end of the file

  //
  /////////Acrylic panels on edges//////////////////////////////////////////
  //

  // Position of horizontal panel respect to vertical panel
  G4ThreeVector traslHorizUpperAcrylicInnerPanelEdge = (verticalAcrylicInnerPanelZ - upperHorizontalAcrylicInnerPanelZ) - (horizontalAcrylicInnerPanelEdgeX - verticalAcrylicInnerPanelEdgeX);
  G4ThreeVector traslHorizLowerAcrylicInnerPanelEdge = -(verticalAcrylicInnerPanelZ - lowerHorizontalAcrylicInnerPanelZ) - (horizontalAcrylicInnerPanelEdgeX - verticalAcrylicInnerPanelEdgeX);

  // Vertical acrylic panel + upper horizontal acrylic panel
  G4UnionSolid* tmpAcrylicInnerPanelEdge = new G4UnionSolid("tmpAcrylicInnerPanelEdge", myVerticalAcrylicInnerPanelEdge, myUpperHorizontalAcrylicInnerPanelEdge, 0, traslHorizUpperAcrylicInnerPanelEdge);

  // Vertical acrylic panel + upper horizontal acrylic panel + lower horizontal
  // acrylic panel
  G4UnionSolid* myAcrylicInnerPanelEdge = new G4UnionSolid("AcrylicInnerPanelEdge", tmpAcrylicInnerPanelEdge, myLowerHorizontalAcrylicInnerPanelEdge, 0, traslHorizLowerAcrylicInnerPanelEdge);

  G4LogicalVolume* myLogicShapeAcrylicInnerPanelEdge[8];

  for (int i = 0; i < 8; i++) {

    myLogicShapeAcrylicInnerPanelEdge[i] = new G4LogicalVolume(myAcrylicInnerPanelEdge, acrylic_mat, "innerAcrylicPanelEdge");

    myLogicShapeAcrylicInnerPanelEdge[i]->SetVisAttributes(VisualBlue);
  }

  // Placement at the end of the file

  //////////////////////////////////////////////////////////////////
  /////////CYLINDER CLOSING HORIZONTAL INNER PANELS CENTER/////////
  ////////////////////////////////////////////////////////////////////

  const G4double innerCylinderOffset = DSStorage::Get()->GetDS20kinnerCylinderOffset();  // Offset of inner TPB cylinders from
                                                                                         // horizontal panels
  // and from the extremes of the vertical space in which they are placed

  const G4double innerTPBCylinderThickness = DSStorage::Get()->GetDS20kinnerTPBCylinderThickness();  // Inner TPB cylinder
                                                                                                     // thickness. Must be lower
                                                                                                     // than cylinder outer radius

  //
  // INNER TPB CYLINDER
  //

  const G4double r2InnerTPBCylinder = centerInnerOffset - innerCylinderOffset;         // Inner cylinder external radius
  const G4double r1InnerTPBCylinder = r2InnerTPBCylinder - innerTPBCylinderThickness;  // Inner cylinder internal radius

  const G4double halfHeightUpperInnerTPBCylinder = myLArBufferThickness / 2. - innerCylinderOffset;  // Half upper cylinder height
  const G4double halfHeightLowerInnerTPBCylinder = myLArBufferThickness / 2. - innerCylinderOffset;  // Half lower cylinder height

  ////////////////////Upper inner cylinder///////////////////////////
  G4Tubs* myUpperInnerTPBCylinder = new G4Tubs("upperInnerTPBCylinder", r1InnerTPBCylinder, r2InnerTPBCylinder, halfHeightUpperInnerTPBCylinder, 0. * deg, 360. * deg);

  G4LogicalVolume* myLogicShapeUpperInnerTPBCylinder = new G4LogicalVolume(myUpperInnerTPBCylinder, TPB_mat, "UpperInnerTPBCylinder");
  myLogicShapeUpperInnerTPBCylinder->SetVisAttributes(VisualRed);

  // Upper inner cylinder position: at the center of the buffer with the two
  // offsets up and down
  G4ThreeVector posUpperInnerTPBCylinder = G4ThreeVector(0, 0, effectiveDS20kTPCheight / 2. + halfHeightUpperInnerTPBCylinder + innerCylinderOffset);

  // Placement at the and of the file

  //////////////////////////////Lower inner cylinder//////////////////
  G4Tubs* myLowerInnerTPBCylinder = new G4Tubs("lowerInnerTPBCylinder", r1InnerTPBCylinder, r2InnerTPBCylinder, halfHeightLowerInnerTPBCylinder, 0. * deg, 360. * deg);

  G4LogicalVolume* myLogicShapeLowerInnerTPBCylinder = new G4LogicalVolume(myLowerInnerTPBCylinder, TPB_mat, "LowerInnerTPBCylinder");
  myLogicShapeLowerInnerTPBCylinder->SetVisAttributes(VisualRed);

  // Lower inner cylinder position
  G4ThreeVector posLowerInnerTPBCylinder = -G4ThreeVector(0, 0, effectiveDS20kTPCheight / 2. + halfHeightLowerInnerTPBCylinder + innerCylinderOffset);

  // Placement at the end of the file

  //
  // INNER ACRYLIC CYLINDERS
  //

  const G4double r1InnerAcrylicCylinder = r1InnerTPBCylinder + myTPBThickness;  // Inner cylinder internal radius
  const G4double r2InnerAcrylicCylinder = r2InnerTPBCylinder - myTPBThickness;  // Inner cylinder external radius

  const G4double halfHeightUpperInnerAcrylicCylinder = halfHeightUpperInnerTPBCylinder;  // Half upper cylinder height
  const G4double halfHeightLowerInnerAcrylicCylinder = halfHeightLowerInnerTPBCylinder;  // Half lower cylinder height

  ////////////////////Upper inner cylinder///////////////////////////
  G4Tubs* myUpperInnerAcrylicCylinder = new G4Tubs("upperInnerAcrylicCylinder", r1InnerAcrylicCylinder, r2InnerAcrylicCylinder, halfHeightUpperInnerAcrylicCylinder, 0. * deg, 360. * deg);

  G4LogicalVolume* myLogicShapeUpperInnerAcrylicCylinder = new G4LogicalVolume(myUpperInnerAcrylicCylinder, acrylic_mat, "UpperInnerAcrylicCylinder");
  myLogicShapeUpperInnerAcrylicCylinder->SetVisAttributes(VisualBlue);

  // Placement at the and of the file

  //////////////////////////////Lower inner cylinder//////////////////
  G4Tubs* myLowerInnerAcrylicCylinder = new G4Tubs("lowerInnerAcrylicCylinder", r1InnerAcrylicCylinder, r2InnerAcrylicCylinder, halfHeightLowerInnerAcrylicCylinder, 0. * deg, 360. * deg);

  G4LogicalVolume* myLogicShapeLowerInnerAcrylicCylinder = new G4LogicalVolume(myLowerInnerAcrylicCylinder, acrylic_mat, "LowerInnerAcrylicCylinder");
  myLogicShapeLowerInnerAcrylicCylinder->SetVisAttributes(VisualBlue);

  // Placement at the end of the file

  ///////////////////////////////////////////////////////////////////
  ////////////////////////INNER LAr BUFFER///////////////////////////
  ///////////////////////////////////////////////////////////////////

  const G4double heightInnerLArBuffer = effectiveDS20kTPCheight + 2. * myLArBufferThickness;  // 2.* for upper and lower spacing

  const G4double inscribedInnerLArBuffer = inscribedTPCThick + myLArBufferThickness;  // Inscribed radius of inner LAr buffer

  G4double zInnerLArBuffer[2] = {-heightInnerLArBuffer / 2., heightInnerLArBuffer / 2.};
  G4double r1InnerLArBuffer[2] = {0, 0};
  G4double r2InnerLArBuffer[2] = {inscribedInnerLArBuffer, inscribedInnerLArBuffer};

  G4Polyhedra* InnerLArBuffer = new G4Polyhedra("InnerLArBuffer", 0. * deg, 360. * deg, 8, 2, zInnerLArBuffer, r1InnerLArBuffer, r2InnerLArBuffer);

  G4LogicalVolume* myLogicShapeInnerLArBuffer = new G4LogicalVolume(InnerLArBuffer, PScintVetoLiquidArgon, "InnerLArBuffer");
  myLogicShapeInnerLArBuffer->SetVisAttributes(VisualYellow);

  // Inner LAr buffer position
  G4ThreeVector posInnerLArBuffer = G4ThreeVector(0, 0, 0);

  // Placement at the end of the file

  ///////////////////////////////////////////////////////////////////
  /////////////////INNER PLASTIC TPB LAYER///////////////////////////
  ///////////////////////////////////////////////////////////////////

  const G4double heightInnerPlasticTPB = heightInnerLArBuffer + 2. * myTPBThickness;         // 2.* for upper and lower spacing
  const G4double inscribedInnerPlasticTPB = inscribedInnerLArBuffer;                         // Inscribed radius of TPB on plastic
  const G4double inscribedInnerPlasticTPBThick = inscribedInnerPlasticTPB + myTPBThickness;  // Inscribed radius of TPB on plastic + its thickness

  G4double zInnerPlasticTPB[2] = {-heightInnerPlasticTPB / 2., heightInnerPlasticTPB / 2.};
  G4double r1InnerPlasticTPB[2] = {0, 0};
  G4double r2InnerPlasticTPB[2] = {inscribedInnerPlasticTPBThick, inscribedInnerPlasticTPBThick};

  G4Polyhedra* innerPlasticTPB = new G4Polyhedra("innerPlasticTPB", 0. * deg, 360. * deg, 8, 2, zInnerPlasticTPB, r1InnerPlasticTPB, r2InnerPlasticTPB);

  G4LogicalVolume* myLogicShapeInnerPlasticTPB = new G4LogicalVolume(innerPlasticTPB, TPB_mat, "innerPlasticTPB");
  myLogicShapeInnerPlasticTPB->SetVisAttributes(VisualBlue);

  // Inner plastic TPB position
  G4ThreeVector posInnerPlasticTPB = G4ThreeVector(0, 0, 0);

  // Placement at the end of the file

  ///////////////////////////////////////////////////////////////////////////////
  //////////////////   PLASTIC AND LAr BUFFER IN PLASTIC //////////////////////
  ///////////////////////////////////////////////////////////////////////////////

  // Thickness of LAr buffer between two plastic solids
  const G4double plasticLArBufferThickness = DSStorage::Get()->GetDS20kplasticLArBufferThickness();  // Must be lower than plastic
                                                                                                     // thickness

  const G4double heightPlastic = heightInnerPlasticTPB + 2. * myPlasticThickness + 2. * plasticLArBufferThickness;  // 2.* for upper and lower spacing

  const G4double inscribedPlastic = inscribedInnerPlasticTPBThick;                                           // Inscribed radius of plastic
  const G4double inscribedPlasticThick = inscribedPlastic + myPlasticThickness + plasticLArBufferThickness;  // Inscribed radius of plastic + its thickness
                                                                                                             //  + LAr buffer thickness

  G4double zPlastic[2] = {-heightPlastic / 2., heightPlastic / 2.};
  G4double r1Plastic[2] = {0, 0};
  G4double r2Plastic[2] = {inscribedPlasticThick, inscribedPlasticThick};

  G4Polyhedra* plastic = new G4Polyhedra("plastic", 0. * deg, 360. * deg, 8, 2, zPlastic, r1Plastic, r2Plastic);

  G4LogicalVolume* myLogicShapePlastic = new G4LogicalVolume(plastic, plastic_mat, "plastic");
  myLogicShapePlastic->SetVisAttributes(VisualGreen);

  // Inner plastic position
  G4ThreeVector posPlastic = G4ThreeVector(0, 0, 0);

  // Placement at the end of the file

  ///////////////////////////////////////////////////////////////////////////////
  //////////////////////////////LAr BUFFER IN
  ///PLASTIC////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////

  // This volume is completely inside the plastic. We need to create a shell
  // (hollow) and place it in the plastic. The shell will be built from an outer
  // LAr buffer and an inner LAr buffer and their subtraction solid.

  const G4double heightOuterPlasticLArBuffer = heightInnerPlasticTPB + myPlasticThickness + 2. * plasticLArBufferThickness;  // Height outer plastic LAr buffer

  const G4double inscribedOuterPlasticLArBuffer = (inscribedPlastic + inscribedPlasticThick) / 2. - plasticLArBufferThickness / 2.;  // Inscribed radius of plastic LAr buffer
  const G4double inscribedOuterPlasticLArBufferThick = inscribedOuterPlasticLArBuffer + plasticLArBufferThickness;                   // Inscribed radius of LAr buffer + its
                                                                                                                                     // thickness

  G4double zOuterPlasticLArBuffer[2] = {-heightOuterPlasticLArBuffer / 2., heightOuterPlasticLArBuffer / 2.};
  G4double r1OuterPlasticLArBuffer[2] = {0, 0};
  G4double r2OuterPlasticLArBuffer[2] = {inscribedOuterPlasticLArBufferThick, inscribedOuterPlasticLArBufferThick};

  G4Polyhedra* outerPlasticLArBuffer = new G4Polyhedra("outerPlasticLArBuffer", 0. * deg, 360. * deg, 8, 2, zOuterPlasticLArBuffer, r1OuterPlasticLArBuffer, r2OuterPlasticLArBuffer);

  const G4double heightInnerPlasticLArBuffer = heightOuterPlasticLArBuffer - 2. * plasticLArBufferThickness;  // Height inner plastic LAr buffer

  const G4double inscribedInnerPlasticLArBuffer = inscribedOuterPlasticLArBuffer;  // Inscribed radius of inner LAr buffer

  G4double zInnerPlasticLArBuffer[2] = {-heightInnerPlasticLArBuffer / 2., heightInnerPlasticLArBuffer / 2.};
  G4double r1InnerPlasticLArBuffer[2] = {0, 0};
  G4double r2InnerPlasticLArBuffer[2] = {inscribedInnerPlasticLArBuffer, inscribedInnerPlasticLArBuffer};

  G4Polyhedra* innerPlasticLArBuffer = new G4Polyhedra("innerPlasticLArBuffer", 0. * deg, 360. * deg, 8, 2, zInnerPlasticLArBuffer, r1InnerPlasticLArBuffer, r2InnerPlasticLArBuffer);

  G4SubtractionSolid* plasticLArBuffer = new G4SubtractionSolid("plasticLArBuffer", outerPlasticLArBuffer, innerPlasticLArBuffer);

  G4LogicalVolume* myLogicShapePlasticLArBuffer = new G4LogicalVolume(plasticLArBuffer, NSLiquidArgon, "plasticLArBuffer");
  myLogicShapePlasticLArBuffer->SetVisAttributes(VisualPurple);

  // Plastic LAr buffer position
  G4ThreeVector posPlasticLArBuffer = G4ThreeVector(0, 0, 0);

  // Placement at the end of the file

  ///////////////////////////////////////////////////////////////////
  /////////////////OUTER PLASTIC TPB LAYER///////////////////////////
  ///////////////////////////////////////////////////////////////////

  const G4double heightOuterPlasticTPB = heightPlastic + 2. * myTPBThickness;  // 2.* for upper and lower spacing

  const G4double inscribedOuterPlasticTPB = inscribedPlasticThick;                           // Inscribed radius of TPB out of plastic
  const G4double inscribedOuterPlasticTPBThick = inscribedOuterPlasticTPB + myTPBThickness;  // Inscribed radius of TPB out of plastic + its thickness

  const G4double apothemOuterPlasticTPB = inscribedOuterPlasticTPBThick;
  const G4ThreeVector aOuterPlasticTPB = G4ThreeVector(apothemOuterPlasticTPB, 0, 0);

  const G4double radiusOuterPlasticTPB = apothemOuterPlasticTPB / (((1. + sqrt(2.)) / 2.) * sqrt(2. - sqrt(2.)));  // From geom. considerations
  const G4ThreeVector rOuterPlasticTPB = G4ThreeVector(radiusOuterPlasticTPB, 0, 0);

  G4double zOuterPlasticTPB[2] = {-heightOuterPlasticTPB / 2., heightOuterPlasticTPB / 2.};
  G4double r1OuterPlasticTPB[2] = {0, 0};
  G4double r2OuterPlasticTPB[2] = {inscribedOuterPlasticTPBThick, inscribedOuterPlasticTPBThick};

  G4Polyhedra* OuterPlasticTPB = new G4Polyhedra("OuterPlasticTPB", 0. * deg, 360. * deg, 8, 2, zOuterPlasticTPB, r1OuterPlasticTPB, r2OuterPlasticTPB);

  G4LogicalVolume* myLogicShapeOuterPlasticTPB = new G4LogicalVolume(OuterPlasticTPB, TPB_mat, "OuterPlasticTPB");
  myLogicShapeOuterPlasticTPB->SetVisAttributes(VisualBlue);

  // Outer plastic TPB position
  G4ThreeVector posOuterPlasticTPB = G4ThreeVector(0, 0, 0);

  // Placement at the end of the file

  ////////////////////////////////////////////////////////////////
  //////////////////////OUTER TPB PANELS//////////////////////////
  ////////////////////////////////////////////////////////////////

  ////////////////////////////////////////
  ///////Vertical outer TPB panels///////////////////////////////////////
  ////////////////////////////////////////

  // Vertical TPB panels on faces
  const G4double verticalOuterTPBPanelFacex = myLArBufferThickness - externalOffset - verticalOuterPanelOffset;  // Need to subtract panel gaps
  // Vertical TPB panels on edges
  const G4double verticalOuterTPBPanelEdgex = (verticalOuterTPBPanelFacex + externalOffset + verticalOuterPanelOffset) / (sin(67.5 * deg)) - myVerticalOuterTPBPanely * 0.5 * (1. / (tan(67.5 * deg))) - extOff - verticalOuterPanelOff;  // Panels on edges must be wider to fit the
                                                                                                                                                                                                                                          // octagonal shape. Calculus from geom.
                                                                                                                                                                                                                                          // considerations

  // Vertical TPB panels height
  const G4double verticalOuterTPBPanelz = heightOuterPlasticTPB + 2. * myLArBufferThickness - upperHorizontalOuterPanelOffset - lowerHorizontalOuterPanelOffset;  // 2.* for upper and lower spacing

  // Useful positions for vertical TPB panels on faces
  G4ThreeVector verticalOuterTPBPanelFaceX = G4ThreeVector(verticalOuterTPBPanelFacex / 2., 0, 0);  // Position of panel geom. center

  // Useful positions for vertical TPB panels on edges
  G4ThreeVector verticalOuterTPBPanelEdgeX = G4ThreeVector(verticalOuterTPBPanelEdgex / 2., 0, 0);  // Position of panel geom. center

  // Useful position for vertical TPB panels
  G4ThreeVector verticalOuterTPBPanelZ = G4ThreeVector(0, 0,
                                                       verticalOuterTPBPanelz / 2.);  // Vertical position of vertical outer panel geom. center

  // Solid for vertical TPB panels on faces
  G4Box* myVerticalOuterTPBPanelFace = new G4Box("verticalOuterTPBPanelFace", verticalOuterTPBPanelFacex / 2., myVerticalOuterTPBPanely / 2., verticalOuterTPBPanelz / 2.);
  // Solid for vertical TPB panels on edges
  G4Box* myVerticalOuterTPBPanelEdge = new G4Box("verticalOuterTPBPanelEdge", verticalOuterTPBPanelEdgex / 2., myVerticalOuterTPBPanely / 2., verticalOuterTPBPanelz / 2.);

  ////////////////////////////////////////
  ///////Horizontal outer TPB panels///////////////////////////////////////
  ////////////////////////////////////////

  // Horizontal TPB panels on faces
  const G4double horizontalOuterTPBPanelFacex = verticalOuterTPBPanelFacex + externalOffset + apothemOuterPlasticTPB - centerOuterOffset;
  // Horizontal TPB panels on edges
  const G4double horizontalOuterTPBPanelEdgex = (horizontalOuterTPBPanelFacex + centerOuterOffset) / ((((1. + sqrt(2.)) / 2.) * sqrt(2. - sqrt(2.)))) - myHorizontalOuterTPBPanely * 0.5 * (1. / (tan(67.5 * deg))) - centerOuterOffset;  // Edges panels must be wider to fit the octagonal
                                                                                                                                                                                                                                          // shape. Calculus from geom. considerations

  // Horizontal TPB panels height
  const G4double upperHorizontalOuterTPBPanelz = myLArBufferThickness - upperHorizontalOuterPanelOffset - externalOffset;
  const G4double lowerHorizontalOuterTPBPanelz = myLArBufferThickness - lowerHorizontalOuterPanelOffset - externalOffset;

  // Useful positions for upper horizontal TPB panels
  G4ThreeVector upperHorizontalOuterTPBPanelZ = G4ThreeVector(0, 0, upperHorizontalOuterTPBPanelz / 2.);

  // Useful positions for lower horizontal TPB panels
  G4ThreeVector lowerHorizontalOuterTPBPanelZ = G4ThreeVector(0, 0, lowerHorizontalOuterTPBPanelz / 2.);

  // Useful positions for TPB panels on faces
  G4ThreeVector horizontalOuterTPBPanelFaceX = G4ThreeVector(horizontalOuterTPBPanelFacex / 2., 0,
                                                             0);  // Position of panel geom. center respect to z axis

  // Solid for upper horizontal TPB panels on faces
  G4Box* myUpperHorizontalOuterTPBPanelFace = new G4Box("upperHorizontalOuterTPBPanelFace", horizontalOuterTPBPanelFacex / 2., myHorizontalOuterTPBPanely / 2., upperHorizontalOuterTPBPanelz / 2.);

  // Solid for lower horizontal TPB panels on faces
  G4Box* myLowerHorizontalOuterTPBPanelFace = new G4Box("lowerHorizontalOuterTPBPanelFace", horizontalOuterTPBPanelFacex / 2., myHorizontalOuterTPBPanely / 2., lowerHorizontalOuterTPBPanelz / 2.);

  // Useful positions for TPB panels on edges
  G4ThreeVector horizontalOuterTPBPanelEdgeX = G4ThreeVector(horizontalOuterTPBPanelEdgex / 2., 0,
                                                             0);  // Position of panel geom. center respect to z axis

  // Solid for upper horizontal TPB panels on edges
  G4Box* myUpperHorizontalOuterTPBPanelEdge = new G4Box("upperHorizontalOuterTPBPanelEdge", horizontalOuterTPBPanelEdgex / 2., myHorizontalOuterTPBPanely / 2., upperHorizontalOuterTPBPanelz / 2.);

  // Solid for lower horizontal TPB panels on edges
  G4Box* myLowerHorizontalOuterTPBPanelEdge = new G4Box("lowerHorizontalOuterTPBPanelEdge", horizontalOuterTPBPanelEdgex / 2., myHorizontalOuterTPBPanely / 2., lowerHorizontalOuterTPBPanelz / 2.);

  /////////////////////////////////////////////////////////////
  /////////////Union solid and placement for outer TPB panels//////
  /////////////////////////////////////////////////////////////

  //
  /////////TPB panels on faces//////////////////////////////////////////
  //

  // Position of horizontal TPB panel respect to vertical TPB panel
  G4ThreeVector traslHorizUpperOuterTPBPanelFace = (verticalOuterTPBPanelZ - upperHorizontalOuterTPBPanelZ) - (horizontalOuterTPBPanelFaceX - verticalOuterTPBPanelFaceX);
  G4ThreeVector traslHorizLowerOuterTPBPanelFace = -(verticalOuterTPBPanelZ - lowerHorizontalOuterTPBPanelZ) - (horizontalOuterTPBPanelFaceX - verticalOuterTPBPanelFaceX);

  // Vertical TPB panel + upper horizontal TPB panel
  G4UnionSolid* tmpOuterTPBPanelFace = new G4UnionSolid("tmpOuterTPBPanelFace", myVerticalOuterTPBPanelFace, myUpperHorizontalOuterTPBPanelFace, 0, traslHorizUpperOuterTPBPanelFace);

  // Vertical TPB panel + upper horizontal TPB panel + lower horizontal TPB
  // panel
  G4UnionSolid* myOuterTPBPanelFace = new G4UnionSolid("outerPanelFace", tmpOuterTPBPanelFace, myLowerHorizontalOuterTPBPanelFace, 0, traslHorizLowerOuterTPBPanelFace);

  G4LogicalVolume* myLogicShapeOuterTPBPanelFace[8];

  for (int i = 0; i < 8; i++) {

    myLogicShapeOuterTPBPanelFace[i] = new G4LogicalVolume(myOuterTPBPanelFace, TPB_mat, "OuterTPBPanelFace");

    myLogicShapeOuterTPBPanelFace[i]->SetVisAttributes(VisualRed);
  }

  // Placement at the end of the file

  //
  /////////TPB panels on edges//////////////////////////////////////////
  //

  // Position of horizontal TPB panel respect to vertical TPB panel
  G4ThreeVector traslHorizUpperOuterTPBPanelEdge = (verticalOuterTPBPanelZ - upperHorizontalOuterTPBPanelZ) - (horizontalOuterTPBPanelEdgeX - verticalOuterTPBPanelEdgeX);
  G4ThreeVector traslHorizLowerOuterTPBPanelEdge = -(verticalOuterTPBPanelZ - lowerHorizontalOuterTPBPanelZ) - (horizontalOuterTPBPanelEdgeX - verticalOuterTPBPanelEdgeX);

  // Vertical TPB panel + upper horizontal TPB panel
  G4UnionSolid* tmpOuterTPBPanelEdge = new G4UnionSolid("tmpOuterTPBPanelEdge", myVerticalOuterTPBPanelEdge, myUpperHorizontalOuterTPBPanelEdge, 0, traslHorizUpperOuterTPBPanelEdge);

  // Vertical TPB panel + upper horizontal TPB panel + lower horizontal TPB
  // panel
  G4UnionSolid* myOuterTPBPanelEdge = new G4UnionSolid("OuterTPBPanelEdge", tmpOuterTPBPanelEdge, myLowerHorizontalOuterTPBPanelEdge, 0, traslHorizLowerOuterTPBPanelEdge);

  G4LogicalVolume* myLogicShapeOuterTPBPanelEdge[8];

  for (int i = 0; i < 8; i++) {

    myLogicShapeOuterTPBPanelEdge[i] = new G4LogicalVolume(myOuterTPBPanelEdge, TPB_mat, "OuterTPBPanelEdge");

    myLogicShapeOuterTPBPanelEdge[i]->SetVisAttributes(VisualRed);
  }

  // Placement at the end of the file

  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////OUTER ARYLIC PANELS/////////////////////////////
  //////////////////////////////////////////////////////////////////////////////

  //
  //////////////////Vertical acrylic panels///////////////
  //

  // Dimensions defined from vertical TPB panels dimensions
  const G4double verticalAcrylicOuterPanely = myVerticalOuterTPBPanely - 2. * myTPBThickness;  // Panel thickness
  // Vertical acrylic panels on faces
  const G4double verticalAcrylicOuterPanelFacex = verticalOuterTPBPanelFacex - 2. * myTPBThickness;
  // Vertical acrylic panels on edges
  const G4double verticalAcrylicOuterPanelEdgex = verticalOuterTPBPanelEdgex - 2. * myTPBThickness;
  // Vertical acrylic panels height
  const G4double verticalAcrylicOuterPanelz = verticalOuterTPBPanelz - 2. * myTPBThickness;

  // Solid for vertical acrylic panels on faces
  G4Box* myVerticalAcrylicOuterPanelFace = new G4Box("verticalAcrylicOuterPanelFace", verticalAcrylicOuterPanelFacex / 2., verticalAcrylicOuterPanely / 2., verticalAcrylicOuterPanelz / 2.);

  // Solid for vertical acrylic panels on edges
  G4Box* myVerticalAcrylicOuterPanelEdge = new G4Box("verticalAcrylicOuterPanelEdge", verticalAcrylicOuterPanelEdgex / 2., verticalAcrylicOuterPanely / 2., verticalAcrylicOuterPanelz / 2.);

  // Useful positions for vertical acrylic panels on faces
  G4ThreeVector verticalAcrylicOuterPanelFaceX = G4ThreeVector(verticalAcrylicOuterPanelFacex / 2., 0,
                                                               0);  // Position of panel geom. center

  // Useful positions for vertical acrylic panels on edges
  G4ThreeVector verticalAcrylicOuterPanelEdgeX = G4ThreeVector(verticalAcrylicOuterPanelEdgex / 2., 0,
                                                               0);  // Position of panel geom. center

  // Useful position for vertical acrylic panels
  G4ThreeVector verticalAcrylicOuterPanelZ = G4ThreeVector(0, 0,
                                                           verticalAcrylicOuterPanelz / 2.);  // Vertical position of vertical inner panel geom. center

  //
  ////////////////Horizontal panels////////////////////
  //

  // Dimensions defined from horizontal TPB panels dimensions
  const G4double horizontalAcrylicOuterPanely = myHorizontalOuterTPBPanely - 2. * myTPBThickness;  // Panel thickness
  // Horizontal panels on faces
  const G4double horizontalAcrylicOuterPanelFacex = horizontalOuterTPBPanelFacex - 2. * myTPBThickness;
  // Horizontal panel on edges
  const G4double horizontalAcrylicOuterPanelEdgex = horizontalOuterTPBPanelEdgex - 2. * myTPBThickness;

  // Horizontal panels height
  const G4double upperHorizontalAcrylicOuterPanelz = upperHorizontalOuterTPBPanelz - 2. * myTPBThickness;
  const G4double lowerHorizontalAcrylicOuterPanelz = lowerHorizontalOuterTPBPanelz - 2. * myTPBThickness;

  // Solid for upper horizontal acrylic panels on faces
  G4Box* myUpperHorizontalAcrylicOuterPanelFace = new G4Box("upperHorizontalAcrylicOuterPanelFace", horizontalAcrylicOuterPanelFacex / 2., horizontalAcrylicOuterPanely / 2., upperHorizontalAcrylicOuterPanelz / 2.);

  // Solid for lower horizontal acrylic panels on faces
  G4Box* myLowerHorizontalAcrylicOuterPanelFace = new G4Box("lowerHorizontalAcrylicOuterPanelFace", horizontalAcrylicOuterPanelFacex / 2., horizontalAcrylicOuterPanely / 2., lowerHorizontalAcrylicOuterPanelz / 2.);

  // Solid for upper horizontal acrylic panels on edges
  G4Box* myUpperHorizontalAcrylicOuterPanelEdge = new G4Box("upperHorizontalAcrylicOuterPanelEdge", horizontalAcrylicOuterPanelEdgex / 2., horizontalAcrylicOuterPanely / 2., upperHorizontalAcrylicOuterPanelz / 2.);

  // Solid for lower horizontal acrylic panels on edges
  G4Box* myLowerHorizontalAcrylicOuterPanelEdge = new G4Box("lowerHorizontalAcrylicOuterPanelEdge", horizontalAcrylicOuterPanelEdgex / 2., horizontalAcrylicOuterPanely / 2., lowerHorizontalAcrylicOuterPanelz / 2.);

  // Useful positions for upper horizontal acrylic panels
  G4ThreeVector upperHorizontalAcrylicOuterPanelZ = G4ThreeVector(0, 0,
                                                                  upperHorizontalAcrylicOuterPanelz / 2.);  // Position of panel geom. center

  // Useful positions for lower horizontal acrylic panels
  G4ThreeVector lowerHorizontalAcrylicOuterPanelZ = G4ThreeVector(0, 0,
                                                                  lowerHorizontalAcrylicOuterPanelz / 2.);  // Position of panel geom. center

  // Useful positions for acrylic panels on faces
  G4ThreeVector horizontalAcrylicOuterPanelFaceX = G4ThreeVector(horizontalAcrylicOuterPanelFacex / 2., 0,
                                                                 0);  // Position of panel geom. center

  // Useful positions for acrylic panels on edges
  G4ThreeVector horizontalAcrylicOuterPanelEdgeX = G4ThreeVector(horizontalAcrylicOuterPanelEdgex / 2., 0,
                                                                 0);  // Position of panel geom. center

  /////////////////////////////////////////////////////////////
  /////Union solid and placement for outer acrylic panels//////
  /////////////////////////////////////////////////////////////

  //
  /////////Acrylic panels on faces//////////////////////////////////////////
  //

  // Position of horizontal panel respect to vertical panel
  G4ThreeVector traslHorizUpperAcrylicOuterPanelFace = (verticalAcrylicOuterPanelZ - upperHorizontalAcrylicOuterPanelZ) - (horizontalAcrylicOuterPanelFaceX - verticalAcrylicOuterPanelFaceX);
  G4ThreeVector traslHorizLowerAcrylicOuterPanelFace = -(verticalAcrylicOuterPanelZ - lowerHorizontalAcrylicOuterPanelZ) - (horizontalAcrylicOuterPanelFaceX - verticalAcrylicOuterPanelFaceX);

  // Vertical acrylic panel + upper horizontal acrylic panel
  G4UnionSolid* tmpAcrylicOuterPanelFace = new G4UnionSolid("tmpAcrylicOuterPanelFace", myVerticalAcrylicOuterPanelFace, myUpperHorizontalAcrylicOuterPanelFace, 0, traslHorizUpperAcrylicOuterPanelFace);

  // Vertical acrylic panel + upper horizontal acrylic panel + lower horizontal
  // acrylic panel
  G4UnionSolid* myAcrylicOuterPanelFace = new G4UnionSolid("innerAcrylicPanelFace", tmpAcrylicOuterPanelFace, myLowerHorizontalAcrylicOuterPanelFace, 0, traslHorizLowerAcrylicOuterPanelFace);

  G4LogicalVolume* myLogicShapeAcrylicOuterPanelFace[8];

  for (int i = 0; i < 8; i++) {

    myLogicShapeAcrylicOuterPanelFace[i] = new G4LogicalVolume(myAcrylicOuterPanelFace, acrylic_mat, "acrylicOuterPanelFace");

    myLogicShapeAcrylicOuterPanelFace[i]->SetVisAttributes(VisualBlue);
  }

  // Placement at the end of the file

  //
  /////////Acrylic panels on edges//////////////////////////////////////////
  //

  // Position of horizontal panel respect to vertical panel
  G4ThreeVector traslHorizUpperAcrylicOuterPanelEdge = (verticalAcrylicOuterPanelZ - upperHorizontalAcrylicOuterPanelZ) - (horizontalAcrylicOuterPanelEdgeX - verticalAcrylicOuterPanelEdgeX);
  G4ThreeVector traslHorizLowerAcrylicOuterPanelEdge = -(verticalAcrylicOuterPanelZ - lowerHorizontalAcrylicOuterPanelZ) - (horizontalAcrylicOuterPanelEdgeX - verticalAcrylicOuterPanelEdgeX);

  // Vertical acrylic panel + upper horizontal acrylic panel
  G4UnionSolid* tmpAcrylicOuterPanelEdge = new G4UnionSolid("tmpAcrylicOuterPanelEdge", myVerticalAcrylicOuterPanelEdge, myUpperHorizontalAcrylicOuterPanelEdge, 0, traslHorizUpperAcrylicOuterPanelEdge);

  // Vertical acrylic panel + upper horizontal acrylic panel + lower horizontal
  // acrylic panel
  G4UnionSolid* myAcrylicOuterPanelEdge = new G4UnionSolid("AcrylicOuterPanelEdge", tmpAcrylicOuterPanelEdge, myLowerHorizontalAcrylicOuterPanelEdge, 0, traslHorizLowerAcrylicOuterPanelEdge);

  G4LogicalVolume* myLogicShapeAcrylicOuterPanelEdge[8];

  for (int i = 0; i < 8; i++) {

    myLogicShapeAcrylicOuterPanelEdge[i] = new G4LogicalVolume(myAcrylicOuterPanelEdge, acrylic_mat, "innerAcrylicPanelEdge");

    myLogicShapeAcrylicOuterPanelEdge[i]->SetVisAttributes(VisualBlue);
  }

  // Placement at the end of the file

  //////////////////////////////////////////////////////////////////
  /////////CYLINDER CLOSING HORIZONTAL OUTER PANELS CENTER//////////
  //////////////////////////////////////////////////////////////////

  const G4double outerCylinderOffset = DSStorage::Get()->GetDS20kouterCylinderOffset();  // Offset of outer TPB cylinders from
                                                                                         // horizontal panels
  // and from the extremes of the vertical space in which they are placed

  const G4double outerTPBCylinderThickness = DSStorage::Get()->GetDS20kouterTPBCylinderThickness();  // TPB cylinder thickness. Must
                                                                                                     // be lower than cylinder outer
                                                                                                     // radius

  //
  //// OUTER TPB CYLINDERS
  //

  const G4double r2OuterTPBCylinder = centerOuterOffset - outerCylinderOffset;         // Outer cylinder external radius
  const G4double r1OuterTPBCylinder = r2OuterTPBCylinder - outerTPBCylinderThickness;  // Outer cylinder internal radius

  const G4double halfHeightUpperOuterTPBCylinder = myLArBufferThickness / 2. - outerCylinderOffset;  // Half upper cylinder height
  const G4double halfHeightLowerOuterTPBCylinder = myLArBufferThickness / 2. - outerCylinderOffset;  // Half lower cylinder height

  ////////////////////Upper outer cylinder///////////////////////////
  G4Tubs* myUpperOuterTPBCylinder = new G4Tubs("upperOuterTPBCylinder", r1OuterTPBCylinder, r2OuterTPBCylinder, halfHeightUpperOuterTPBCylinder, 0. * deg, 360. * deg);

  G4LogicalVolume* myLogicShapeUpperOuterTPBCylinder = new G4LogicalVolume(myUpperOuterTPBCylinder, TPB_mat, "UpperOuterTPBCylinder");
  myLogicShapeUpperOuterTPBCylinder->SetVisAttributes(VisualRed);

  // Upper inner cylinder position: at the center of the buffer with the two
  // offset up and down
  G4ThreeVector posUpperOuterTPBCylinder = G4ThreeVector(0, 0, heightOuterPlasticTPB / 2. + halfHeightUpperOuterTPBCylinder + outerCylinderOffset);

  // Placement at the end of the file

  //////////////////////////////Lower outer cylinder//////////////////
  G4Tubs* myLowerOuterTPBCylinder = new G4Tubs("lowerOuterTPBCylinder", r1OuterTPBCylinder, r2OuterTPBCylinder, halfHeightLowerOuterTPBCylinder, 0. * deg, 360. * deg);

  G4LogicalVolume* myLogicShapeLowerOuterTPBCylinder = new G4LogicalVolume(myLowerOuterTPBCylinder, TPB_mat, "LowerOuterTPBCylinder");
  myLogicShapeLowerOuterTPBCylinder->SetVisAttributes(VisualRed);

  // Lower inner cylinder position
  G4ThreeVector posLowerOuterTPBCylinder = -G4ThreeVector(0, 0, heightOuterPlasticTPB / 2. + halfHeightLowerOuterTPBCylinder + outerCylinderOffset);

  // Placement at the end of the file

  //
  // OUTER ACRYLIC CYLINDERS
  //

  const G4double r1OuterAcrylicCylinder = r1OuterTPBCylinder + myTPBThickness;  // Inner cylinder internal radius
  const G4double r2OuterAcrylicCylinder = r2OuterTPBCylinder - myTPBThickness;  // Inner cylinder external radius

  const G4double halfHeightUpperOuterAcrylicCylinder = halfHeightUpperOuterTPBCylinder;  // Half upper cylinder height
  const G4double halfHeightLowerOuterAcrylicCylinder = halfHeightLowerOuterTPBCylinder;  // Half lower cylinder height

  ////////////////////Upper inner cylinder///////////////////////////
  G4Tubs* myUpperOuterAcrylicCylinder = new G4Tubs("upperOuterAcrylicCylinder", r1OuterAcrylicCylinder, r2OuterAcrylicCylinder, halfHeightUpperOuterAcrylicCylinder, 0. * deg, 360. * deg);

  G4LogicalVolume* myLogicShapeUpperOuterAcrylicCylinder = new G4LogicalVolume(myUpperOuterAcrylicCylinder, acrylic_mat, "UpperOuterAcrylicCylinder");
  myLogicShapeUpperOuterAcrylicCylinder->SetVisAttributes(VisualBlue);

  // Placement at the and of the file

  //////////////////////////////Lower inner cylinder//////////////////
  G4Tubs* myLowerOuterAcrylicCylinder = new G4Tubs("lowerOuterAcrylicCylinder", r1OuterAcrylicCylinder, r2OuterAcrylicCylinder, halfHeightLowerOuterAcrylicCylinder, 0. * deg, 360. * deg);

  G4LogicalVolume* myLogicShapeLowerOuterAcrylicCylinder = new G4LogicalVolume(myLowerOuterAcrylicCylinder, acrylic_mat, "LowerOuterAcrylicCylinder");
  myLogicShapeLowerOuterAcrylicCylinder->SetVisAttributes(VisualBlue);

  // Placement at the end of the file

  ///////////////////////////////////////////////////////////////////
  ////////////////////////OUTER LAr BUFFER///////////////////////////
  ///////////////////////////////////////////////////////////////////

  const G4double heightOuterLArBuffer = heightOuterPlasticTPB + 2. * myLArBufferThickness;  // 2.* for upper and lower spacing

  const G4double inscribedOuterLArBuffer = inscribedOuterPlasticTPBThick + myLArBufferThickness;  // Inscribed radius of TPB on plastic

  G4double zOuterLArBuffer[2] = {-heightOuterLArBuffer / 2., heightOuterLArBuffer / 2.};
  G4double r1OuterLArBuffer[2] = {0, 0};
  G4double r2OuterLArBuffer[2] = {inscribedOuterLArBuffer, inscribedOuterLArBuffer};

  G4Polyhedra* OuterLArBuffer = new G4Polyhedra("OuterLArBuffer", 0. * deg, 360. * deg, 8, 2, zOuterLArBuffer, r1OuterLArBuffer, r2OuterLArBuffer);

  G4LogicalVolume* myLogicShapeOuterLArBuffer = new G4LogicalVolume(OuterLArBuffer, VetoLiquidArgon, "OuterLArBuffer");
  myLogicShapeOuterLArBuffer->SetVisAttributes(VisualYellow);

  // Outer LAr buffer position
  G4ThreeVector posOuterLArBuffer = G4ThreeVector(0, 0, 0);

  // Placement at the end of the file

  ///////////////////////////////////////////////////////////////////
  /////////////////EXTERNAL TPB LAYER///////////////////////////
  ///////////////////////////////////////////////////////////////////

  const G4double heightExternalTPB = heightOuterLArBuffer + 2. * myTPBThickness;  // 2.* for upper and lower spacing

  const G4double inscribedExternalTPB = inscribedOuterLArBuffer;                     // Inscribed radius of external TPB
  const G4double inscribedExternalTPBThick = inscribedExternalTPB + myTPBThickness;  // Inscribed radius of external TPB + its thickness

  G4double zExternalTPB[2] = {-heightExternalTPB / 2., heightExternalTPB / 2.};
  G4double r1ExternalTPB[2] = {0, 0};
  G4double r2ExternalTPB[2] = {inscribedExternalTPBThick, inscribedExternalTPBThick};

  G4Polyhedra* ExternalTPB = new G4Polyhedra("ExternalTPB", 0. * deg, 360. * deg, 8, 2, zExternalTPB, r1ExternalTPB, r2ExternalTPB);

  G4LogicalVolume* myLogicShapeExternalTPB = new G4LogicalVolume(ExternalTPB, TPB_mat, "ExternalTPB");
  myLogicShapeExternalTPB->SetVisAttributes(VisualBlue);

  // External TPB position
  G4ThreeVector posExternalTPB = G4ThreeVector(0, 0, 0);

  // Placement at the and of the file

  ///////////////////////////////////////////////////////////////////
  /////////////////EXTERNAL COPPER///////////////////////////
  ///////////////////////////////////////////////////////////////////

  // Thickness of copper shell
  const G4double myCopperThickness = DSStorage::Get()->GetDS20kCopperThickness();

  const G4double heightCopper = heightExternalTPB + 2. * myCopperThickness;  // 2.* for upper and lower spacing

  const G4double inscribedCopper = inscribedExternalTPBThick;                 // Inscribed radius of external copper
  const G4double inscribedCopperThick = inscribedCopper + myCopperThickness;  // Inscribed radius of copper + its thickness

  G4double zCopper[2] = {-heightCopper / 2., heightCopper / 2.};
  G4double r1Copper[2] = {0, 0};
  G4double r2Copper[2] = {inscribedCopperThick, inscribedCopperThick};

  G4Polyhedra* Copper = new G4Polyhedra("Copper", 0. * deg, 360. * deg, 8, 2, zCopper, r1Copper, r2Copper);

  G4LogicalVolume* myLogicShapeCopper = new G4LogicalVolume(Copper, copper_mat, "Copper");
  myLogicShapeCopper->SetVisAttributes(VisualPurple);

  // Copper position: traslation of 20 cm in height is consequence of the other
  // classes (TPC,cryostat).
  G4ThreeVector posCopper = G4ThreeVector(0, 0, -20. * cm);

  G4RotationMatrix* rotCopper = new G4RotationMatrix();  // Rotation i.o.t. have a side of all octagonal
                                                         // shapes othogonal to x axis
  rotCopper->rotateZ(ang / 2.);

  // Placement at the and of the file

  //////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////ALL VOLUMES
  ///PLACEMENT///////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////PHYSICAL VOLUMES DEFINED IN HEADER
  ///FILE//////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////

  // Copper placement, defined through physical mother volume

  physCopper = new G4PVPlacement(rotCopper, posCopper, "Copper", myLogicShapeCopper, fMotherVolume, false, 0, myCheckOverlaps);

  // External TPB placement

  physExternalTPB = new G4PVPlacement(0, posExternalTPB, myLogicShapeExternalTPB, "ExternalTPB", myLogicShapeCopper, false, 0, myCheckOverlaps);

  // Outer LAr buffer placement

  physOuterLArBuffer = new G4PVPlacement(0, posOuterLArBuffer, myLogicShapeOuterLArBuffer, "OuterLArBuffer", myLogicShapeExternalTPB, false, 0, myCheckOverlaps);

  // Outer panels on faces placement

  if (DSStorage::Get()->GetDS20kFacePanels()) {

    for (int i = 0; i < 8; i++) {
      G4RotationMatrix* rot = new G4RotationMatrix();
      rot->rotateZ(-i * ang - ang / 2.);  // Minus for rotation in the same direction as position
                                          // vector -ang/2. for positioning the first panel on face

      G4ThreeVector posOuterPanelFace = (aOuterPlasticTPB + ExternalOffset + verticalOuterTPBPanelFaceX + G4ThreeVector(0, 0,
                                                                                                                        -0.5 * upperHorizontalOuterPanelOffset + 0.5 * lowerHorizontalOuterPanelOffset)).rotateZ(i * ang + ang / 2.);  //+ang/2. for positioning the first panel on face

      string number;
      stringstream ss;
      ss << i;
      number = ss.str();

      physTPBOuterPanelFace[i] = new G4PVPlacement(rot, posOuterPanelFace, myLogicShapeOuterTPBPanelFace[i], "outerTPBPanelFace" + number, myLogicShapeOuterLArBuffer, false, 0, myCheckOverlaps);

      physAcrylicOuterPanelFace[i] = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), myLogicShapeAcrylicOuterPanelFace[i], "acrylicOuterPanelFace" + number, myLogicShapeOuterTPBPanelFace[i], false, 0, myCheckOverlaps);
    }
  }

  // Outer panels on edges placement

  if (DSStorage::Get()->GetDS20kEdgePanels()) {

    for (int i = 0; i < 8; i++) {
      G4RotationMatrix* rot = new G4RotationMatrix();
      rot->rotateZ(-i * ang);  // Minus for rotation in the same direction as position vector

      G4ThreeVector posOuterPanelEdge = (rOuterPlasticTPB + ExtOff + verticalOuterTPBPanelEdgeX + G4ThreeVector(0, 0, -0.5 * upperHorizontalOuterPanelOffset + 0.5 * lowerHorizontalOuterPanelOffset)).rotateZ(i * ang);

      string number;
      stringstream ss;
      ss << i;
      number = ss.str();

      physTPBOuterPanelEdge[i] = new G4PVPlacement(rot, posOuterPanelEdge, myLogicShapeOuterTPBPanelEdge[i], "outerTPBPanelEdge" + number, myLogicShapeOuterLArBuffer, false, 0, myCheckOverlaps);

      physAcrylicOuterPanelEdge[i] = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), myLogicShapeAcrylicOuterPanelEdge[i], "acrylicOuterPanelEdge" + number, myLogicShapeOuterTPBPanelEdge[i], false, 0, myCheckOverlaps);
    }
  }

  // Outer cylinders placement

  physUpperOuterTPBCylinder = new G4PVPlacement(0, posUpperOuterTPBCylinder, myLogicShapeUpperOuterTPBCylinder, "UpperOuterTPBCylinder", myLogicShapeOuterLArBuffer, false, 0, myCheckOverlaps);

  physUpperOuterAcrylicCylinder = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), myLogicShapeUpperOuterAcrylicCylinder, "UpperOuterAcrylicCylinder", myLogicShapeUpperOuterTPBCylinder, false, 0, myCheckOverlaps);

  physLowerOuterTPBCylinder = new G4PVPlacement(0, posLowerOuterTPBCylinder, myLogicShapeLowerOuterTPBCylinder, "LowerOuterTPBCylinder", myLogicShapeOuterLArBuffer, false, 0, myCheckOverlaps);

  physLowerOuterAcrylicCylinder = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), myLogicShapeLowerOuterAcrylicCylinder, "LowerOuterAcrylicCylinder", myLogicShapeLowerOuterTPBCylinder, false, 0, myCheckOverlaps);

  // Outer plastic TPB placement

  physOuterPlasticTPB = new G4PVPlacement(0, posOuterPlasticTPB, myLogicShapeOuterPlasticTPB, "OuterPlasticTPB", myLogicShapeOuterLArBuffer, false, 0, myCheckOverlaps);

  // Plastic placement

  physPlastic = new G4PVPlacement(0, posPlastic, myLogicShapePlastic, "plastic", myLogicShapeOuterPlasticTPB, false, 0, myCheckOverlaps);

  // Plastic LAr buffer placement

  physPlasticLArBuffer = new G4PVPlacement(0, posPlasticLArBuffer, myLogicShapePlasticLArBuffer, "plasticLArBuffer", myLogicShapePlastic, false, 0, myCheckOverlaps);

  // Inner plastic TPB placement

  physInnerPlasticTPB = new G4PVPlacement(0, posInnerPlasticTPB, myLogicShapeInnerPlasticTPB, "innerPlasticTPB", myLogicShapePlastic, false, 0, myCheckOverlaps);

  // Inner LAr buffer placement

  physInnerLArBuffer = new G4PVPlacement(0, posInnerLArBuffer, myLogicShapeInnerLArBuffer, "innerLArBuffer", myLogicShapeInnerPlasticTPB, false, 0, myCheckOverlaps);

  ///////////////////////////////////////////////////////////////////////////////////////
  fPhysicInnerLiquidArgon = physInnerLArBuffer;  // Needed for setting TPC mother volume
  ///////////////////////////////////////////////////////////////////////////////////////

  // Inner panels on faces placement

  if (DSStorage::Get()->GetDS20kFacePanels()) {

    for (int i = 0; i < 8; i++) {
      G4RotationMatrix* rot = new G4RotationMatrix();
      rot->rotateZ(-i * ang - ang / 2.);  // Minus for rotation in the same direction as position
                                          // vector -ang/2. for positioning the first panel on face

      G4ThreeVector posInnerPanelFace = (aTPCThick + VerticalInnerPanelOffset + verticalInnerTPBPanelFaceX).rotateZ(i * ang + ang / 2.);  //+ang/2. for positioning the first panel on face

      string number;
      stringstream ss;
      ss << i;
      number = ss.str();

      physTPBInnerPanelFace[i] = new G4PVPlacement(rot, posInnerPanelFace, myLogicShapeInnerTPBPanelFace[i], "innerTPBPanelFace" + number, myLogicShapeInnerLArBuffer, false, 0, myCheckOverlaps);

      physAcrylicInnerPanelFace[i] = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), myLogicShapeAcrylicInnerPanelFace[i], "acrylicInnerPanelFace" + number, myLogicShapeInnerTPBPanelFace[i], false, 0, myCheckOverlaps);
    }
  }

  // Inner panels on edges placement

  if (DSStorage::Get()->GetDS20kEdgePanels()) {

    for (int i = 0; i < 8; i++) {
      G4RotationMatrix* rot = new G4RotationMatrix();
      rot->rotateZ(-i * ang);  // Minus for rotation in the same direction as position vector

      G4ThreeVector posInnerPanelEdge = (rTPCThick + VerticalInnerPanelOff + verticalInnerTPBPanelEdgeX).rotateZ(i * ang);

      string number;
      stringstream ss;
      ss << i;
      number = ss.str();

      physTPBInnerPanelEdge[i] = new G4PVPlacement(rot, posInnerPanelEdge, myLogicShapeInnerTPBPanelEdge[i], "innerTPBPanelEdge" + number, myLogicShapeInnerLArBuffer, false, 0, myCheckOverlaps);

      physAcrylicInnerPanelEdge[i] = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), myLogicShapeAcrylicInnerPanelEdge[i], "acrylicInnerPanelEdge" + number, myLogicShapeInnerTPBPanelEdge[i], false, 0, myCheckOverlaps);
    }
  }

  // Inner cylinders placement

  physUpperInnerTPBCylinder = new G4PVPlacement(0, posUpperInnerTPBCylinder, myLogicShapeUpperInnerTPBCylinder, "UpperInnerTPBCylinder", myLogicShapeInnerLArBuffer, false, 0, myCheckOverlaps);

  physUpperInnerAcrylicCylinder = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), myLogicShapeUpperInnerAcrylicCylinder, "UpperInnerAcrylicCylinder", myLogicShapeUpperInnerTPBCylinder, false, 0, myCheckOverlaps);

  physLowerInnerTPBCylinder = new G4PVPlacement(0, posLowerInnerTPBCylinder, myLogicShapeLowerInnerTPBCylinder, "LowerInnerTPBCylinder", myLogicShapeInnerLArBuffer, false, 0, myCheckOverlaps);

  physLowerInnerAcrylicCylinder = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), myLogicShapeLowerInnerAcrylicCylinder, "LowerInnerAcrylicCylinder", myLogicShapeLowerInnerTPBCylinder, false, 0, myCheckOverlaps);

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  //////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////  SiPMs PLACEMENT
  /////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////

  const G4bool SiPMs = DSStorage::Get()->GetDS20kSiPMs();

  if (SiPMs) {
    DSLog(routine) << "SiPMs building and placement activated" << endlog;
  } else {
    DSLog(routine) << "SiPMs building and placement not activated" << endlog;
  }

  if (SiPMs) {

    // Creating SiPM objects

    const G4double nSiPMInside = DSStorage::Get()->GetDS20knSiPMInside();
    const G4double nSiPMOutside = DSStorage::Get()->GetDS20knSiPMOutside();

    /////   All the following dimensions and positions are derived from the
    ///previous code   /////

    const G4double SiPMSide = 5.2 * cm;
    const G4double SiPMHeight = 0.5 * cm;

    // Total number of SiPMs effectively placed
    G4int nSiPMsEffective = 0;

    // Copy number for SiPMs
    //
    // The syntax of the copy number will be the following: abcde , where:
    //  - a=1 for internal buffer, a=2 for external buffer
    //  - b=0,...7 for each sector
    //  - cde is the number of the SiPM in the sector (from 000 to xxx for SiPMs
    //  on sides, from 400 to 4xx for SiPMs on lower caps, from 500 to 5xx for
    //  SiPMs on upper caps)
    //
    G4int internalSidesSiPMCopyNumber;
    G4int externalSidesSiPMCopyNumber;
    G4int internalUpperCapSiPMCopyNumber;
    G4int internalLowerCapSiPMCopyNumber;
    G4int externalUpperCapSiPMCopyNumber;
    G4int externalLowerCapSiPMCopyNumber;

    // Output file for mapping SiPMs. Useful for root file
    G4String outputFilename = DSIO::Get()->GetBinaryFileName();
    G4String siPMOutputFilename = outputFilename.substr(0, outputFilename.length() - 4);
    G4String newSiPMOutputFilename = siPMOutputFilename.append("_SiPM_placement.txt");

    ofstream outputSiPMs;
    outputSiPMs.open(newSiPMOutputFilename.c_str());

    // Read SiPMs positions from file
    G4bool SiPMsAutoplacement = DSStorage::Get()->GetDS20kSiPMsAutoplacement();

    // Positions in input file must be in centimeters !!!
    ifstream inputSiPMs;
    G4String inputFilename = DSStorage::Get()->GetDS20kSiPMsAutoplacementFilename();

    // TPB on SiPMs or not
    G4bool TPBOnSiPM = DSStorage::Get()->GetDS20kTPBOnSiPM();

    // Uniform distribution for SiPMs or not
    G4bool SiPMUniformInsideSides = false;
    G4bool SiPMUniformInsideCaps = false;
    G4bool SiPMUniformOutsideSides = false;
    G4bool SiPMUniformOutsideCaps = false;

    // In DSStorage one can choose a different fraction of SiPMs on caps respect
    // to the total. This fraction is initialized to zero; if changed, the code
    // will use it.
    G4double SiPMFractionInside = 0;
    G4double SiPMFractionOutside = 0.4;

    // Setting parameters in case of best light collection is chosen.
    // The settings for UV absorption length and TPB visible reflectivity need
    // to be set by hand.

    const G4bool bestLightCollection = DSStorage::Get()->GetDS20kbestLightCollection();

    if (bestLightCollection && SiPMsAutoplacement) {

      SiPMUniformInsideSides = true;
      SiPMUniformInsideCaps = true;
      SiPMUniformOutsideSides = false;
      SiPMUniformOutsideCaps = false;

      SiPMFractionInside = 0;
      SiPMFractionOutside = 0.4;
      // Other properties of best light collection configuration are set later

    }

    else if (!bestLightCollection && SiPMsAutoplacement) {

      SiPMUniformInsideSides = DSStorage::Get()->GetDS20kSiPMUniformInsideSides();
      SiPMUniformInsideCaps = DSStorage::Get()->GetDS20kSiPMUniformInsideCaps();
      SiPMUniformOutsideSides = DSStorage::Get()->GetDS20kSiPMUniformOutsideSides();
      SiPMUniformOutsideCaps = DSStorage::Get()->GetDS20kSiPMUniformOutsideCaps();

      SiPMFractionInside = DSStorage::Get()->GetDS20kSiPMFractionInside();
      SiPMFractionOutside = DSStorage::Get()->GetDS20kSiPMFractionOutside();

    }

    else {

      inputSiPMs.open(inputFilename.c_str());
    }

    ////////////////////////////////////   On sides
    /////////////////////////////////////

    G4Box* SiPMOnSide = new G4Box("SiPMOnSide", SiPMHeight / 2., SiPMSide / 2., SiPMSide / 2.);
    G4LogicalVolume* myLogicShapeSiPMOnSide = new G4LogicalVolume(SiPMOnSide, teflon_mat, "SiPMOnSide");
    myLogicShapeSiPMOnSide->SetVisAttributes(myCyan);

    G4Box* SiPMActiveOnSide = new G4Box("SiPMActiveOnSide", 0.2 * cm, 2.5 * cm, 2.5 * cm);
    G4LogicalVolume* myLogicShapeSiPMActiveOnSide = new G4LogicalVolume(SiPMActiveOnSide, metalsilicon_mat, "SiPMActiveOnSide");
    myLogicShapeSiPMActiveOnSide->SetVisAttributes(myBlue);

    // Set TPB on SIPMs or not. If not the same volume will be made of fused
    // silica

    if (TPBOnSiPM) {

      G4Box* SiPMSilicaOnSide = new G4Box("SiPMSilicaOnSide", 250 * micrometer, 2.5 * cm, 2.5 * cm);
      G4LogicalVolume* myLogicShapeSiPMSilicaOnSide = new G4LogicalVolume(SiPMSilicaOnSide, fusedsilica_mat, "SiPMSilicaOnSide");
      myLogicShapeSiPMSilicaOnSide->SetVisAttributes(myYellow);

      G4Box* SiPMTPBOnSide = new G4Box("SiPMTPBOnSide", 50 * micrometer, 2.5 * cm, 2.5 * cm);
      G4LogicalVolume* myLogicShapeSiPMTPBOnSide = new G4LogicalVolume(SiPMTPBOnSide, TPB_mat, "SiPMTPBOnSide");
      myLogicShapeSiPMTPBOnSide->SetVisAttributes(myRed);

      physSiPMSilicaOnSide = new G4PVPlacement(0, G4ThreeVector(-0.2 * cm + 100 * micrometer + 250 * micrometer, 0, 0), myLogicShapeSiPMSilicaOnSide, "SiPMSilicaOnSide", myLogicShapeSiPMActiveOnSide, false, 0, myCheckOverlaps);

      physSiPMTPBOnSide = new G4PVPlacement(0, G4ThreeVector(-0.2 * cm + 50 * micrometer, 0, 0), myLogicShapeSiPMTPBOnSide, "SiPMTPBOnSide", myLogicShapeSiPMActiveOnSide, false, 0, myCheckOverlaps);

    }

    else {

      G4Box* SiPMSilicaOnSide = new G4Box("SiPMSilicaOnSide", 300 * micrometer, 2.5 * cm, 2.5 * cm);
      G4LogicalVolume* myLogicShapeSiPMSilicaOnSide = new G4LogicalVolume(SiPMSilicaOnSide, fusedsilica_mat, "SiPMSilicaOnSide");
      myLogicShapeSiPMSilicaOnSide->SetVisAttributes(myYellow);

      physSiPMSilicaOnSide = new G4PVPlacement(0, G4ThreeVector(-0.2 * cm + 300 * micrometer, 0, 0), myLogicShapeSiPMSilicaOnSide, "SiPMSilicaOnSide", myLogicShapeSiPMActiveOnSide, false, 0, myCheckOverlaps);
    }

    physSiPMActiveOnSide = new G4PVPlacement(0, G4ThreeVector(-SiPMHeight / 2. + 0.2 * cm, 0, 0), myLogicShapeSiPMActiveOnSide, "SiPMActiveOnSide", myLogicShapeSiPMOnSide, false, 0, myCheckOverlaps);

    //////////////////////////////////   On caps ///////////////////////////////

    G4Box* SiPMOnCap = new G4Box("SiPMOnCap", SiPMSide / 2., SiPMSide / 2., SiPMHeight / 2.);
    G4LogicalVolume* myLogicShapeSiPMOnCap = new G4LogicalVolume(SiPMOnCap, teflon_mat, "SiPMOnCap");
    myLogicShapeSiPMOnCap->SetVisAttributes(myCyan);

    G4Box* SiPMActiveOnCap = new G4Box("SiPMActiveOnCap", 2.5 * cm, 2.5 * cm, 0.2 * cm);
    G4LogicalVolume* myLogicShapeSiPMActiveOnCap = new G4LogicalVolume(SiPMActiveOnCap, metalsilicon_mat, "SiPMActiveOnCap");
    myLogicShapeSiPMActiveOnCap->SetVisAttributes(myBlue);

    // Set TPB on SIPMs or not. If not the same volume will be made of fused
    // silica

    if (TPBOnSiPM) {

      G4Box* SiPMSilicaOnCap = new G4Box("SiPMSilicaOnCap", 2.5 * cm, 2.5 * cm, 250. * micrometer);
      G4LogicalVolume* myLogicShapeSiPMSilicaOnCap = new G4LogicalVolume(SiPMSilicaOnCap, fusedsilica_mat, "SiPMSilicaOnCap");
      myLogicShapeSiPMSilicaOnCap->SetVisAttributes(myYellow);

      G4Box* SiPMTPBOnCap = new G4Box("SiPMTPBOnCap", 2.5 * cm, 2.5 * cm, 50. * micrometer);
      G4LogicalVolume* myLogicShapeSiPMTPBOnCap = new G4LogicalVolume(SiPMTPBOnCap, TPB_mat, "SiPMTPBOnCap");
      myLogicShapeSiPMTPBOnCap->SetVisAttributes(myRed);

      physSiPMSilicaOnCap = new G4PVPlacement(0, G4ThreeVector(0, 0, -0.2 * cm + 100 * micrometer + 250 * micrometer), myLogicShapeSiPMSilicaOnCap, "SiPMSilicaOnCap", myLogicShapeSiPMActiveOnCap, false, 0, myCheckOverlaps);

      physSiPMTPBOnCap = new G4PVPlacement(0, G4ThreeVector(0, 0, -0.2 * cm + 50 * micrometer), myLogicShapeSiPMTPBOnCap, "SiPMTPBOnCap", myLogicShapeSiPMActiveOnCap, false, 0, myCheckOverlaps);

    }

    else {

      G4Box* SiPMSilicaOnCap = new G4Box("SiPMSilicaOnCap", 2.5 * cm, 2.5 * cm, 300. * micrometer);
      G4LogicalVolume* myLogicShapeSiPMSilicaOnCap = new G4LogicalVolume(SiPMSilicaOnCap, fusedsilica_mat, "SiPMSilicaOnCap");
      myLogicShapeSiPMSilicaOnCap->SetVisAttributes(myYellow);

      physSiPMSilicaOnCap = new G4PVPlacement(0, G4ThreeVector(0, 0, -0.2 * cm + 300 * micrometer), myLogicShapeSiPMSilicaOnCap, "SiPMSilicaOnCap", myLogicShapeSiPMActiveOnCap, false, 0, myCheckOverlaps);
    }

    physSiPMActiveOnCap = new G4PVPlacement(0, G4ThreeVector(0, 0, -SiPMHeight / 2. + 0.2 * cm), myLogicShapeSiPMActiveOnCap, "SiPMActiveOnCap", myLogicShapeSiPMOnCap, false, 0, myCheckOverlaps);

    // If autoplacement is not activated, the code will skip to the manual
    // placement part
    if (SiPMsAutoplacement) {

      const G4double edgeInnerLArBuffer = inscribedInnerLArBuffer / ((1. + sqrt(2.)) / 2.);  // From geom. considerations

      const G4double edgeOuterPlasticTPBThick = inscribedOuterPlasticTPBThick / ((1. + sqrt(2.)) / 2.);  // From geom. considerations

      // SiPMs will not be placed on the face or caps of the inner LAr buffer
      // but a bit more internally SiPMs will not be placed on the face or caps
      // of the outer plastic TPB but a bit more externally
      const G4double epsilonSiPM = DSStorage::Get()->GetDS20kepsilonSiPM();

      // Ratio between lateral area and bases area of inner LAr buffer
      const G4double areaRatioInside = heightInnerLArBuffer / inscribedInnerLArBuffer;

      // Ratio between lateral area and bases area of outer plastic TPB
      const G4double areaRatioOutside = heightOuterPlasticTPB / inscribedOuterPlasticTPBThick;

      /////////////////////////////////////////////////////////////////////////////////////////////////
      //////////////////////////////      SiPMs on faces
      /////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////

      // Number of SiPMs per half face of inner LAr buffer, considering areas
      // ratio
      G4double nSiPMPerHalfFaceInside = nSiPMInside * (areaRatioInside / (areaRatioInside + 1)) * (1. / 8.) * (1 / 2.);

      // Number of SiPMs per half face of outer plastic TPB, considering areas
      // ratio
      G4double nSiPMPerHalfFaceOutside = nSiPMOutside * (areaRatioOutside / (areaRatioOutside + 1)) * (1. / 8.) * (1 / 2.);

      // Comparison of fraction of SiPMs on faces automatically calculated with
      // user defined via DSStorage.
      if (SiPMFractionInside != 0) { nSiPMPerHalfFaceInside = (nSiPMInside - (nSiPMInside * SiPMFractionInside)) * (1. / 8.) * (1. / 2.); }

      if (SiPMFractionOutside != 0) { nSiPMPerHalfFaceOutside = (nSiPMOutside - (nSiPMOutside * SiPMFractionOutside)) * (1. / 8.) * (1. / 2.); }

      // Ratio between effective half edge and height of a face of the inner LAr
      // buffer
      const G4double lengthRatioInside = (edgeInnerLArBuffer / 2. - myVerticalInnerTPBPanely / 2.) / heightInnerLArBuffer;

      // Ratio between effective half edge and height of a face of the outer
      // plastic TPB
      const G4double lengthRatioOutside = (edgeOuterPlasticTPBThick / 2. - myVerticalOuterTPBPanely / 2.) / heightOuterPlasticTPB;

      // The following numbers (number of SiPMs on faces, inside and outside)
      // are rounded up. This causes, depending on the theoretical number of
      // SiPMs in DSStorage, a possible excess on the number of SiPMs placed
      // respect to the theoretical one.

      // Number of SiPMs along the edge and the height of a face of the inner
      // LAr buffer
      G4double nSiPMHalfEdgeInside = 1;
      G4double nSiPMHeightInside = 1;

      // Number of SiPMs along the edge and the height of a face of the outer
      // plastic TPB
      G4double nSiPMHalfEdgeOutside = 1;
      G4double nSiPMHeightOutside = 1;

      // Considering lenghts ratio, one can calculate the optimal number of
      // SiPMs on faces
      while (nSiPMHalfEdgeInside * nSiPMHeightInside < nSiPMPerHalfFaceInside) {

        nSiPMHeightInside = nSiPMHeightInside + 1;
        nSiPMHalfEdgeInside = round(lengthRatioInside * nSiPMHeightInside);
      }

      while (nSiPMHalfEdgeOutside * nSiPMHeightOutside < nSiPMPerHalfFaceOutside) {

        nSiPMHeightOutside = nSiPMHeightOutside + 1;
        nSiPMHalfEdgeOutside = round(lengthRatioOutside * nSiPMHeightOutside);
      }

      ///////  Best light collection configuration property   ////////

      // In the IAB, the simple non uniform parabola distributions does not
      // work. The best result is obtained placing 5 less rows and using uniform
      // distribution both on caps and sides.

      if (bestLightCollection) { nSiPMHeightInside = nSiPMHeightInside - 5; }

      nSiPMsEffective = 16 * (nSiPMHalfEdgeInside * nSiPMHeightInside + nSiPMHalfEdgeOutside * nSiPMHeightOutside);

      // Since the placement of the SiPM begins a-step-distance from edge panels
      // along the face direction, unless the edge step is too small there will
      // not be overlap between SiPM and panels Even if SiPMs will not be placed
      // on the face of the LAr buffer but a bit more internally and on the face
      // of the outer plastic TPB but a bit more externally

      // We can find the optimal step:

      const G4double edgeStepInside = ((edgeInnerLArBuffer / 2. - myVerticalInnerTPBPanely / 2.) - nSiPMHalfEdgeInside * SiPMSide) / (nSiPMHalfEdgeInside + 1);
      const G4double heightStepInside = (heightInnerLArBuffer - nSiPMHeightInside * SiPMSide) / (nSiPMHeightInside + 1);

      const G4double edgeStepOutside = ((edgeOuterPlasticTPBThick / 2. - myVerticalOuterTPBPanely / 2.) - nSiPMHalfEdgeOutside * SiPMSide) / (nSiPMHalfEdgeOutside + 1);
      const G4double heightStepOutside = (heightOuterPlasticTPB - nSiPMHeightOutside * SiPMSide) / (nSiPMHeightOutside + 1);

      for (int k = 0; k < 8; k++) {

        G4RotationMatrix* rotSiPMInside = new G4RotationMatrix();

        rotSiPMInside->rotateZ(-k * ang - ang / 2.);  // Rotation i.o.t. have the side of
                                                      // SiPMs parallel to LAr buffer edge

        G4RotationMatrix* rotSiPMOutside = new G4RotationMatrix();

        rotSiPMOutside->rotateZ(-k * ang - ang / 2.);  // Rotation i.o.t. have the side of
                                                       // SiPMs parallel to LAr buffer edge

        rotSiPMOutside->rotateY(180. * deg);  // External SiPMs must be facing outside

        internalSidesSiPMCopyNumber = 10000 + k * 1000;
        externalSidesSiPMCopyNumber = 20000 + k * 1000;

        /////////////////////////////////
        ////Internal SiPMs
        /////////////////////////////////

        //////////////////////////////////////Uniform distribution for
        ///SiPMs/////////////////////////////

        if (SiPMUniformInsideSides) {

          // Columns
          for (int i = 1; i <= nSiPMHalfEdgeInside; i++) {

            // Rows
            for (int j = 1; j <= nSiPMHeightInside; j++) {

              // Constructing geometry as if X axis is in the direction of an
              // edge. Positions will be rotated later.
              G4ThreeVector posSiPMInside1 = G4ThreeVector(inscribedInnerLArBuffer - epsilonSiPM - SiPMHeight / 2., myVerticalInnerTPBPanely / 2. + i * edgeStepInside + (2 * i - 1) * SiPMSide / 2., -heightInnerLArBuffer / 2. + j * heightStepInside + (2 * j - 1) * (SiPMSide / 2.));

              G4ThreeVector posSiPMInside2 = G4ThreeVector(inscribedInnerLArBuffer - epsilonSiPM - SiPMHeight / 2., -(myVerticalInnerTPBPanely / 2. + i * edgeStepInside + (2 * i - 1) * SiPMSide / 2.), -heightInnerLArBuffer / 2. + j * heightStepInside + (2 * j - 1) * (SiPMSide / 2.));

              G4ThreeVector pos1 = posSiPMInside1.rotateZ(k * ang + ang / 2.);
              G4ThreeVector pos2 = posSiPMInside2.rotateZ(k * ang + ang / 2.);

              physSiPMInternalSide_1 = new G4PVPlacement(rotSiPMInside, pos1, myLogicShapeSiPMOnSide, "SiPMInternalSide", myLogicShapeInnerLArBuffer, false, internalSidesSiPMCopyNumber, myCheckOverlaps);

              physSiPMInside.push_back(physSiPMInternalSide_1);

              // Since one needs SiPMs coordinates respect to the world and not
              // respect to their mother volume (i.e. copper cage), and since
              //traslation along z axis is balaced, one needs only to de-rotate
              // since the copper cage is rotated respect to the world
              // (cryostat).
              outputSiPMs << internalSidesSiPMCopyNumber << " " << pos1.rotateZ(-ang / 2.) / cm << endl;

              internalSidesSiPMCopyNumber++;

              physSiPMInternalSide_2 = new G4PVPlacement(rotSiPMInside, pos2, myLogicShapeSiPMOnSide, "SiPMInternalSide", myLogicShapeInnerLArBuffer, false, internalSidesSiPMCopyNumber, myCheckOverlaps);

              physSiPMInside.push_back(physSiPMInternalSide_2);

              // Since one needs SiPMs coordinates respect to the world and not
              // respect to their mother volume (i.e. copper cage), and since
              //traslation along z axis is balaced, one needs only to de-rotate
              // since the copper cage is rotated respect to the world
              // (cryostat).
              outputSiPMs << internalSidesSiPMCopyNumber << " " << pos2.rotateZ(-ang / 2.) / cm << endl;

              internalSidesSiPMCopyNumber++;
            }
          }
        }

        ///////////////////////////////////Non uniform distribution for
        ///SiPMs////////////////////////////

        //*****   WARNING: in non uniform case one may visualise geometry
        //configuration to check if all volumes are
        //                 correctly placed.

        else {

          // Parameters from fitting of parabola shape from uniform distribution
          // (assuming zero the linear term) with TPB on SiPMs, free path 100 m
          // and refl= 0.98
          const G4double a = 0.0032 * (1. / (cm * cm));
          const G4double c = 4797;
          const G4double min = 4610;

          // Vector for z positions of SiPMs
          vector<G4double> v_zPosSiPM;

          G4double zPosSiPM;

          // Even number of raws of SiPMs
          if (int(nSiPMHeightInside) % 2 == 0) {

            const G4int nSiPMHalfHeightInside = (nSiPMHeightInside - 4) / 2.;  // The four raws subtracted will be placed later
            const G4int heightStep = (c - min) / (nSiPMHalfHeightInside + 1);

            for (int n = 1; n <= nSiPMHalfHeightInside; n++) {

              zPosSiPM = sqrt((n * heightStep) / a);

              v_zPosSiPM.push_back(zPosSiPM);
              v_zPosSiPM.push_back(-zPosSiPM);

              // Position of the four subtracted raws
              if (n == 1) {

                v_zPosSiPM.push_back(zPosSiPM * 0.2);
                v_zPosSiPM.push_back(-zPosSiPM * 0.2);
                v_zPosSiPM.push_back(zPosSiPM * 0.6);
                v_zPosSiPM.push_back(-zPosSiPM * 0.6);
              }
            }
          }

          // Odd number of raws of SiPMs
          if (int(nSiPMHeightInside) % 2 != 0) {

            const G4int nSiPMHalfHeightInside = (nSiPMHeightInside - 1.) / 2. - 1.;  // The removed raws will be placed later
            const G4double heightStep = (c - min) / (nSiPMHalfHeightInside + 1);

            // Since the number of raws is odd, one raw is surely at half height
            // of the face
            v_zPosSiPM.push_back(0);

            for (int n = 1; n <= nSiPMHalfHeightInside; n++) {

              zPosSiPM = sqrt((heightStep * n) / (a));

              v_zPosSiPM.push_back(zPosSiPM);
              v_zPosSiPM.push_back(-zPosSiPM);

              // Position of the removed raws
              if (n == 1) {

                v_zPosSiPM.push_back(zPosSiPM * 0.5);
                v_zPosSiPM.push_back(-zPosSiPM * 0.5);
              }
            }
          }

          // Raws
          for (int j = 1; j <= nSiPMHeightInside; j++) {

            // Columns
            for (int i = 1; i <= nSiPMHalfEdgeInside; i++) {

              // Constructing geometry as if X axis is in the direction of an
              // edge. Positions will be rotated later.
              G4ThreeVector posSiPMInside1 = G4ThreeVector(inscribedInnerLArBuffer - epsilonSiPM - SiPMHeight / 2., myVerticalInnerTPBPanely / 2. + i * edgeStepInside + (2 * i - 1) * SiPMSide / 2., v_zPosSiPM.at(j - 1));

              G4ThreeVector posSiPMInside2 = G4ThreeVector(inscribedInnerLArBuffer - epsilonSiPM - SiPMHeight / 2., -(myVerticalInnerTPBPanely / 2. + i * edgeStepInside + (2 * i - 1) * SiPMSide / 2.), v_zPosSiPM.at(j - 1));

              G4ThreeVector pos1 = posSiPMInside1.rotateZ(k * ang + ang / 2.);
              G4ThreeVector pos2 = posSiPMInside2.rotateZ(k * ang + ang / 2.);

              physSiPMInternalSide_1 = new G4PVPlacement(rotSiPMInside, pos1, myLogicShapeSiPMOnSide, "SiPMInternalSide", myLogicShapeInnerLArBuffer, false, internalSidesSiPMCopyNumber, myCheckOverlaps);

              physSiPMInside.push_back(physSiPMInternalSide_1);

              // Since one needs SiPMs coordinates respect to the world and not
              // respect to their mother volume (i.e. copper cage), and since
              //traslation along z axis is balaced, one needs only to de-rotate
              // since the copper cage is rotated respect to the world
              // (cryostat).
              outputSiPMs << internalSidesSiPMCopyNumber << " " << pos1.rotateZ(-ang / 2.) / cm << endl;

              internalSidesSiPMCopyNumber++;

              physSiPMInternalSide_2 = new G4PVPlacement(rotSiPMInside, pos2, myLogicShapeSiPMOnSide, "SiPMInternalSide", myLogicShapeInnerLArBuffer, false, internalSidesSiPMCopyNumber, myCheckOverlaps);

              physSiPMInside.push_back(physSiPMInternalSide_2);

              // Since one needs SiPMs coordinates respect to the world and not
              // respect to their mother volume (i.e. copper cage), and since
              //traslation along z axis is balaced, one needs only to de-rotate
              // since the copper cage is rotated respect to the world
              // (cryostat).
              outputSiPMs << internalSidesSiPMCopyNumber << " " << pos2.rotateZ(-ang / 2.) / cm << endl;

              internalSidesSiPMCopyNumber++;
            }
          }
        }

        ///////////////////////////////////
        ////External SiPMs
        ///////////////////////////////////

        //////////////////////////////////////Uniform distribution for
        ///SiPMs/////////////////////////////

        if (SiPMUniformOutsideSides) {

          // Columns
          for (int i = 1; i <= nSiPMHalfEdgeOutside; i++) {

            // Rows
            for (int j = 1; j <= nSiPMHeightOutside; j++) {

              // Constructing geometry as if X axis is in the direction of an
              // edge. Positions will be rotated later.
              G4ThreeVector posSiPMOutside1 = G4ThreeVector(inscribedOuterPlasticTPBThick + epsilonSiPM + SiPMHeight / 2., myVerticalOuterTPBPanely / 2. + i * edgeStepOutside + (2 * i - 1) * SiPMSide / 2., -heightOuterPlasticTPB / 2. + j * heightStepOutside + (2 * j - 1) * (SiPMSide / 2.));

              G4ThreeVector posSiPMOutside2 = G4ThreeVector(inscribedOuterPlasticTPBThick + epsilonSiPM + SiPMHeight / 2., -(myVerticalOuterTPBPanely / 2. + i * edgeStepOutside + (2 * i - 1) * SiPMSide / 2.), -heightOuterPlasticTPB / 2. + j * heightStepOutside + (2 * j - 1) * (SiPMSide / 2.));

              G4ThreeVector pos1 = posSiPMOutside1.rotateZ(k * ang + ang / 2.);
              G4ThreeVector pos2 = posSiPMOutside2.rotateZ(k * ang + ang / 2.);

              physSiPMExternalSide_1 = new G4PVPlacement(rotSiPMOutside, pos1, myLogicShapeSiPMOnSide, "SiPMExternalSide", myLogicShapeOuterLArBuffer, false, externalSidesSiPMCopyNumber, myCheckOverlaps);

              physSiPMOutside.push_back(physSiPMExternalSide_1);

              // Since one needs SiPMs coordinates respect to the world and not
              // respect to their mother volume (i.e. copper cage), and since
              //traslation along z axis is balaced, one needs only to de-rotate
              // since the copper cage is rotated respect to the world
              // (cryostat).
              outputSiPMs << externalSidesSiPMCopyNumber << " " << pos1.rotateZ(-ang / 2.) / cm << endl;

              externalSidesSiPMCopyNumber++;

              physSiPMExternalSide_2 = new G4PVPlacement(rotSiPMOutside, pos2, myLogicShapeSiPMOnSide, "SiPMExternalSide", myLogicShapeOuterLArBuffer, false, externalSidesSiPMCopyNumber, myCheckOverlaps);

              physSiPMOutside.push_back(physSiPMExternalSide_2);

              // Since one needs SiPMs coordinates respect to the world and not
              // respect to their mother volume (i.e. copper cage), and since
              //traslation along z axis is balaced, one needs only to de-rotate
              // since the copper cage is rotated respect to the world
              // (cryostat).
              outputSiPMs << externalSidesSiPMCopyNumber << " " << pos2.rotateZ(-ang / 2.) / cm << endl;

              externalSidesSiPMCopyNumber++;
            }
          }
        }

        ///////////////////////////////////Non uniform distribution for
        ///SiPMs////////////////////////////

        //*****   WARNING: in non uniform case one may visualise geometry
        //configuration to check if all volumes are
        //                 correctly placed.

        else {

          // Parameters from fitting of parabola shape from uniform distribution
          // (assuming zero the linear term) with TPB on SiPMs, free path 100 m
          // and refl= 0.98
          const G4double a = 0.00217 * (1. / (cm * cm));
          const G4double c = 881.856;
          const G4double min = 750.;

          // Vector for z positions of SiPMs
          vector<G4double> v_zPosSiPM;

          G4double zPosSiPM;

          // Even number of raws of SiPMs
          if (int(nSiPMHeightOutside) % 2 == 0) {

            const G4int nSiPMHalfHeightOutside = (nSiPMHeightOutside - 4) / 2.;  // The four raws subtracted will be placed later
            const G4int heightStep = (c - min) / (nSiPMHalfHeightOutside + 1);

            for (int n = 1; n <= nSiPMHalfHeightOutside; n++) {

              zPosSiPM = sqrt((n * heightStep) / a);

              v_zPosSiPM.push_back(zPosSiPM);
              v_zPosSiPM.push_back(-zPosSiPM);

              // Position of the four subtracted raws
              if (n == 1) {

                v_zPosSiPM.push_back(zPosSiPM * 0.2);
                v_zPosSiPM.push_back(-zPosSiPM * 0.2);
                v_zPosSiPM.push_back(zPosSiPM * 0.6);
                v_zPosSiPM.push_back(-zPosSiPM * 0.6);
              }
            }
          }

          // Odd number of raws of SiPMs
          if (int(nSiPMHeightOutside) % 2 != 0) {

            const G4int nSiPMHalfHeightOutside = (nSiPMHeightOutside - 1.) / 2. - 1.;  // The removed raws will be placed later
            const G4double heightStep = (c - min) / (nSiPMHalfHeightOutside + 1);

            // Since the number of raws is odd, one raw is surely at half height
            // of the face
            v_zPosSiPM.push_back(0);

            for (int n = 1; n <= nSiPMHalfHeightOutside; n++) {

              zPosSiPM = sqrt((heightStep * n) / (a));

              v_zPosSiPM.push_back(zPosSiPM);
              v_zPosSiPM.push_back(-zPosSiPM);

              // Position of the removed raws
              if (n == 1) {

                v_zPosSiPM.push_back(zPosSiPM * 0.5);
                v_zPosSiPM.push_back(-zPosSiPM * 0.5);
              }
            }
          }

          // Raws
          for (int j = 1; j <= nSiPMHeightOutside; j++) {

            // Columns
            for (int i = 1; i <= nSiPMHalfEdgeOutside; i++) {

              // Constructing geometry as if X axis is in the direction of an
              // edge. Positions will be rotated later.
              G4ThreeVector posSiPMOutside1 = G4ThreeVector(inscribedOuterPlasticTPBThick + epsilonSiPM + SiPMHeight / 2., myVerticalOuterTPBPanely / 2. + i * edgeStepOutside + (2 * i - 1) * SiPMSide / 2., v_zPosSiPM.at(j - 1));

              G4ThreeVector posSiPMOutside2 = G4ThreeVector(inscribedOuterPlasticTPBThick + epsilonSiPM + SiPMHeight / 2., -(myVerticalOuterTPBPanely / 2. + i * edgeStepOutside + (2 * i - 1) * SiPMSide / 2.), v_zPosSiPM.at(j - 1));

              G4ThreeVector pos1 = posSiPMOutside1.rotateZ(k * ang + ang / 2.);
              G4ThreeVector pos2 = posSiPMOutside2.rotateZ(k * ang + ang / 2.);

              physSiPMExternalSide_1 = new G4PVPlacement(rotSiPMOutside, pos1, myLogicShapeSiPMOnSide, "SiPMExternalSide", myLogicShapeOuterLArBuffer, false, externalSidesSiPMCopyNumber, myCheckOverlaps);

              physSiPMOutside.push_back(physSiPMExternalSide_1);

              // Since one needs SiPMs coordinates respect to the world and not
              // respect to their mother volume (i.e. copper cage), and since
              //traslation along z axis is balaced, one needs only to de-rotate
              // since the copper cage is rotated respect to the world
              // (cryostat).
              outputSiPMs << externalSidesSiPMCopyNumber << " " << pos1.rotateZ(-ang / 2.) / cm << endl;

              externalSidesSiPMCopyNumber++;

              physSiPMExternalSide_2 = new G4PVPlacement(rotSiPMOutside, pos2, myLogicShapeSiPMOnSide, "SiPMExternalSide", myLogicShapeOuterLArBuffer, false, externalSidesSiPMCopyNumber, myCheckOverlaps);

              physSiPMOutside.push_back(physSiPMExternalSide_2);

              // Since one needs SiPMs coordinates respect to the world and not
              // respect to their mother volume (i.e. copper cage), and since
              //traslation along z axis is balaced, one needs only to de-rotate
              // since the copper cage is rotated respect to the world
              // (cryostat).

              outputSiPMs << externalSidesSiPMCopyNumber << " " << pos2.rotateZ(-ang / 2.) / cm << endl;

              externalSidesSiPMCopyNumber++;
            }
          }
        }
      }

      ////////////////////////////////////////////////////////////////////////////////////////////////
      //////////////////////////////      SiPMs on caps
      /////////////////////////////////////////////
      ////////////////////////////////////////////////////////////////////////////////////////////////

      // The following numbers (number of SiPMs on caps, inside and outside) are
      // rounded up, as well as the number of SiPMs on faces. This causes,
      // depending on the theoretical number of SiPMs in DSStorage, a possible
      // excess on the number of SiPMs placed respect to the theoretical one.

      // Number of SiPMs per half sector of inner LAr buffer
      G4int nSiPMPerHalfSectorInside = round((nSiPMInside / (areaRatioInside + 1)) * (1. / 16.) * (1 / 2.));

      ///////  Best light collection configuration property   ////////

      // In the IAB there are too less SiPMs than the theoretical one, thus one
      // force the number of them on the caps One can choose +1.45% that leads
      // to 3024 total SiPMs or +1.35% that leads to 2864 total SiPMs
      if (bestLightCollection) {
        // nSiPMPerHalfSectorInside =
        // round((nSiPMInside/(areaRatioInside+1))*(1./16.)*(1/2.)*(1.35));
        nSiPMPerHalfSectorInside = round((nSiPMInside / (areaRatioInside + 1)) * (1. / 16.) * (1 / 2.) * (1.45));
      }

      // Number of SiPMs per cap of outer plastic TPB
      // G4int nSiPMPerHalfSectorOutside =
      // round(8.*(2.*nSiPMHalfEdgeOutside*nSiPMHeightOutside)/(32.*areaRatioOutside));
      G4int nSiPMPerHalfSectorOutside = round((nSiPMOutside / (areaRatioOutside + 1)) * (1. / 16.) * (1 / 2.));

      // Comparison of fraction of SiPMs on faces automatically calculated with
      // user defined via DSStorage. In DSStorage one can choose a different
      // fraction of SiPMs on caps respect to the total. This fraction is
      // initialized to zero; if changed, the code will use it.

      if (SiPMFractionInside != 0) { nSiPMPerHalfSectorInside = round(nSiPMInside * SiPMFractionInside * (1. / 2.) * (1. / 16.)); }

      if (SiPMFractionOutside != 0) { nSiPMPerHalfSectorOutside = round(nSiPMOutside * SiPMFractionOutside * (1. / 2.) * (1. / 16.)); }

      // Number of SiPms effectively placed on a sector of an internal cap
      G4int placedSiPMsOnInternalCapSector = 0;

      // Number of SiPms effectively placed on a sector of an external cap
      G4int placedSiPMsOnExternalCapSector = 0;

      G4double internalEpsilon = 0.05 * cm;  // Starting step between SiPMs
      const G4double internalAdd = 0.05 * cm;
      G4double internalIncr;

      G4double externalEpsilon = 0.05 * cm;  // Starting step between SiPMs
      const G4double externalAdd = 0.05 * cm;
      G4double externalIncr;

      G4bool recalculateInside = true;
      G4bool recalculateOutside = true;

      vector<G4ThreeVector> SiPMsUpperInternalCapSector1;
      vector<G4ThreeVector> SiPMsUpperInternalCapSector2;
      vector<G4ThreeVector> SiPMsLowerInternalCapSector1;
      vector<G4ThreeVector> SiPMsLowerInternalCapSector2;

      vector<G4ThreeVector> SiPMsUpperExternalCapSector1;
      vector<G4ThreeVector> SiPMsUpperExternalCapSector2;
      vector<G4ThreeVector> SiPMsLowerExternalCapSector1;
      vector<G4ThreeVector> SiPMsLowerExternalCapSector2;

      // SiPMs will be placed at some distance from panels, so that they can't
      // overlap.

      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ////////////////////////////////////////// Internal caps
      //////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      while (recalculateInside) {

        vector<G4int> nSiPMPerColumn;

        placedSiPMsOnInternalCapSector = 0;
        G4int i = 0;

        ///////  Best light collection configuration property   ////////

        // The number of the SiPMs on the internal caps have been increased but
        // there are too many SiPMs accumulated at the center. Thus one does not
        // place the central ones.
        if (bestLightCollection) { i = 2; }

        internalIncr = SiPMSide + internalEpsilon;  // This increment is
                                                    // arbitrary

        // Constructing geometry as if X axis is in the direction of an edge.
        // Positions will be rotated later.
        while (i * internalIncr + SiPMSide / 2. < inscribedInnerLArBuffer) {

          G4double y = (tan(22.5 * deg)) * i * internalIncr - myVerticalInnerTPBPanely / (2. * cos(22.5 * deg)) - myVerticalInnerTPBPanely / 2. - tan(22.5 * deg) * SiPMSide / 2.;  // Taking in account for panel thickness

          G4int n = floor((y + internalEpsilon) / (internalEpsilon + SiPMSide));

          if (n != 0) { nSiPMPerColumn.push_back(n); }

          for (int a = 1; a <= n; a++) {

            G4ThreeVector posSiPMUpperInternalCapSector1 = G4ThreeVector(i * internalIncr, myVerticalInnerTPBPanely / 2. + (y / (n + 1)) * a, heightInnerLArBuffer / 2. - epsilonSiPM - SiPMHeight / 2.);

            G4ThreeVector posSiPMUpperInternalCapSector2 = G4ThreeVector(i * internalIncr, -(myVerticalInnerTPBPanely / 2. + (y / (n + 1)) * a), heightInnerLArBuffer / 2. - epsilonSiPM - SiPMHeight / 2.);

            G4ThreeVector posSiPMLowerInternalCapSector1 = G4ThreeVector(i * internalIncr, myVerticalInnerTPBPanely / 2. + (y / (n + 1)) * a, -(heightInnerLArBuffer / 2. - epsilonSiPM - SiPMHeight / 2.));

            G4ThreeVector posSiPMLowerInternalCapSector2 = G4ThreeVector(i * internalIncr, -(myVerticalInnerTPBPanely / 2. + (y / (n + 1)) * a), -(heightInnerLArBuffer / 2. - epsilonSiPM - SiPMHeight / 2.));

            SiPMsUpperInternalCapSector1.push_back(posSiPMUpperInternalCapSector1);
            SiPMsUpperInternalCapSector2.push_back(posSiPMUpperInternalCapSector2);
            SiPMsLowerInternalCapSector1.push_back(posSiPMLowerInternalCapSector1);
            SiPMsLowerInternalCapSector2.push_back(posSiPMLowerInternalCapSector2);

            placedSiPMsOnInternalCapSector++;
            placedSiPMsOnInternalCapSector++;
          }

          i++;
        }

        if (placedSiPMsOnInternalCapSector * 8 > nSiPMPerHalfSectorInside * 16) {

          internalEpsilon = internalEpsilon + internalAdd;

          SiPMsUpperInternalCapSector1.clear();
          SiPMsUpperInternalCapSector2.clear();
          SiPMsLowerInternalCapSector1.clear();
          SiPMsLowerInternalCapSector2.clear();

          nSiPMPerColumn.clear();

        }

        else {

          recalculateInside = false;

          //*****   WARNING: in non uniform case one may visualise geometry
          //configuration to check if all volumes are
          //                 correctly placed.

          if (SiPMUniformInsideCaps == false) {

            int ncol = nSiPMPerColumn.size();

            // Parameters from fitting of parabola shape from non uniform
            // distribution on sides (assuming zero the linear term) with TPB on
            // SiPMs, free path 100 m and refl= 0.98.

            const G4double a = 0.0041 * (1. / (cm * cm));
            const G4double c = 4851;
            const G4double min = 4639;

            const G4double columnsStep = (c - min) / (ncol + 1);

            G4int doneSiPM = 0;

            for (int n = 1; n <= ncol; n++) {

              for (int nn = 0; nn < nSiPMPerColumn.at(n - 1); nn++) {

                G4double xPosSiPM = sqrt(((n)*columnsStep) / a) * (1 - 0.9 / (n * n));

                if (n == 1) { xPosSiPM = sqrt(((n)*columnsStep) / a) * (0.6); }

                G4ThreeVector posSiPM1 = SiPMsUpperInternalCapSector1.at(doneSiPM + nn);
                G4ThreeVector posSiPM2 = SiPMsUpperInternalCapSector2.at(doneSiPM + nn);
                G4ThreeVector posSiPM3 = SiPMsLowerInternalCapSector1.at(doneSiPM + nn);
                G4ThreeVector posSiPM4 = SiPMsLowerInternalCapSector2.at(doneSiPM + nn);

                posSiPM1.setX(xPosSiPM);
                posSiPM2.setX(xPosSiPM);
                posSiPM3.setX(xPosSiPM);
                posSiPM4.setX(xPosSiPM);

                SiPMsUpperInternalCapSector1[doneSiPM + nn] = posSiPM1;
                SiPMsUpperInternalCapSector2[doneSiPM + nn] = posSiPM2;
                SiPMsLowerInternalCapSector1[doneSiPM + nn] = posSiPM3;
                SiPMsLowerInternalCapSector2[doneSiPM + nn] = posSiPM4;
              }

              doneSiPM += nSiPMPerColumn.at(n - 1);
            }
          }

          for (int k = 0; k < 8; k++) {

            G4RotationMatrix* rotSiPMOnUpperCapInside = new G4RotationMatrix();

            rotSiPMOnUpperCapInside->rotateZ(-k * ang - ang / 2.);  // Rotation i.o.t. have the side of SiPMs
                                                                    // parallel to LAr buffer edge

            G4RotationMatrix* rotSiPMOnLowerCapInside = new G4RotationMatrix();

            rotSiPMOnLowerCapInside->rotateZ(-k * ang - ang / 2.);  // Rotation i.o.t. have the side of SiPMs
                                                                    // parallel to LAr buffer edge

            rotSiPMOnLowerCapInside->rotateY(180. * deg);  // Internal SiPMs on lower cap must be facing up

            internalUpperCapSiPMCopyNumber = 10000 + k * 1000 + 500;
            internalLowerCapSiPMCopyNumber = 10000 + k * 1000 + 400;

            for (int j = 0; j < placedSiPMsOnInternalCapSector / 2; j++) {

              G4ThreeVector pos1 = SiPMsUpperInternalCapSector1.at(j);
              pos1.rotateZ(k * ang + ang / 2.);

              G4ThreeVector pos2 = SiPMsUpperInternalCapSector2.at(j);
              pos2.rotateZ(k * ang + ang / 2.);

              G4ThreeVector pos3 = SiPMsLowerInternalCapSector1.at(j);
              pos3.rotateZ(k * ang + ang / 2.);

              G4ThreeVector pos4 = SiPMsLowerInternalCapSector2.at(j);
              pos4.rotateZ(k * ang + ang / 2.);

              physSiPMInternalUpperCap_1 = new G4PVPlacement(rotSiPMOnUpperCapInside, pos1, myLogicShapeSiPMOnCap, "SiPMInternalUpperCap", myLogicShapeInnerLArBuffer, false, internalUpperCapSiPMCopyNumber, myCheckOverlaps);

              physSiPMInside.push_back(physSiPMInternalUpperCap_1);

              // Since one needs SiPMs coordinates respect to the world and not
              // respect to their mother volume (i.e. copper cage), and since
              //traslation along z axis is balaced, one needs only to de-rotate
              // since the copper cage is rotated respect to the world
              // (cryostat).
              outputSiPMs << internalUpperCapSiPMCopyNumber << " " << pos1.rotateZ(-ang / 2.) / cm << endl;

              internalUpperCapSiPMCopyNumber++;

              physSiPMInternalUpperCap_2 = new G4PVPlacement(rotSiPMOnUpperCapInside, pos2, myLogicShapeSiPMOnCap, "SiPMInternalUpperCap", myLogicShapeInnerLArBuffer, false, internalUpperCapSiPMCopyNumber, myCheckOverlaps);

              physSiPMInside.push_back(physSiPMInternalUpperCap_2);

              // Since one needs SiPMs coordinates respect to the world and not
              // respect to their mother volume (i.e. copper cage), and since
              //traslation along z axis is balaced, one needs only to de-rotate
              // since the copper cage is rotated respect to the world
              // (cryostat).
              outputSiPMs << internalUpperCapSiPMCopyNumber << " " << pos2.rotateZ(-ang / 2.) / cm << endl;

              internalUpperCapSiPMCopyNumber++;

              physSiPMInternalLowerCap_1 = new G4PVPlacement(rotSiPMOnLowerCapInside, pos3, myLogicShapeSiPMOnCap, "SiPMInternalLowerCap", myLogicShapeInnerLArBuffer, false, internalLowerCapSiPMCopyNumber, myCheckOverlaps);

              physSiPMInside.push_back(physSiPMInternalLowerCap_1);

              // Since one needs SiPMs coordinates respect to the world and not
              // respect to their mother volume (i.e. copper cage), and since
              //traslation along z axis is balaced, one needs only to de-rotate
              // since the copper cage is rotated respect to the world
              // (cryostat).
              outputSiPMs << internalLowerCapSiPMCopyNumber << " " << pos3.rotateZ(-ang / 2.) / cm << endl;

              internalLowerCapSiPMCopyNumber++;

              physSiPMInternalLowerCap_2 = new G4PVPlacement(rotSiPMOnLowerCapInside, pos4, myLogicShapeSiPMOnCap, "SiPMInternalLowerCap", myLogicShapeInnerLArBuffer, false, internalLowerCapSiPMCopyNumber, myCheckOverlaps);

              physSiPMInside.push_back(physSiPMInternalLowerCap_2);

              // Since one needs SiPMs coordinates respect to the world and not
              // respect to their mother volume (i.e. copper cage), and since
              //traslation along z axis is balaced, one needs only to de-rotate
              // since the copper cage is rotated respect to the world
              // (cryostat).
              outputSiPMs << internalLowerCapSiPMCopyNumber << " " << pos4.rotateZ(-ang / 2.) / cm << endl;

              internalLowerCapSiPMCopyNumber++;
            }
          }

          nSiPMsEffective += placedSiPMsOnInternalCapSector * 8 * 2;
        }
      }

      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ////////////////////////////////////////// External caps
      //////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      while (recalculateOutside) {

        vector<G4int> nSiPMPerColumn;

        placedSiPMsOnExternalCapSector = 0;
        G4int i = 0;
        externalIncr = SiPMSide + externalEpsilon;  // This increment is
                                                    // arbitrary

        // Constructing geometry as if X axis is in the direction of an edge.
        // Positions will be rotated later.
        while (i * externalIncr + SiPMSide / 2. < inscribedOuterPlasticTPBThick) {

          G4double y = (tan(22.5 * deg)) * i * externalIncr - myVerticalOuterTPBPanely / (2. * cos(22.5 * deg)) - myVerticalOuterTPBPanely / 2. - tan(22.5 * deg) * SiPMSide / 2.;  // Taking in account for panel thickness

          G4int n = floor((y + externalEpsilon) / (externalEpsilon + SiPMSide));

          if (n != 0) { nSiPMPerColumn.push_back(n); }

          for (int a = 1; a <= n; a++) {

            G4ThreeVector posSiPMUpperExternalCapSector1 = G4ThreeVector(i * externalIncr, myVerticalOuterTPBPanely / 2. + (y / (n + 1)) * a, heightOuterPlasticTPB / 2. + epsilonSiPM + SiPMHeight / 2.);

            G4ThreeVector posSiPMUpperExternalCapSector2 = G4ThreeVector(i * externalIncr, -(myVerticalOuterTPBPanely / 2. + (y / (n + 1)) * a), heightOuterPlasticTPB / 2. + epsilonSiPM + SiPMHeight / 2.);

            G4ThreeVector posSiPMLowerExternalCapSector1 = G4ThreeVector(i * externalIncr, myVerticalOuterTPBPanely / 2. + (y / (n + 1)) * a, -(heightOuterPlasticTPB / 2. + epsilonSiPM + SiPMHeight / 2.));

            G4ThreeVector posSiPMLowerExternalCapSector2 = G4ThreeVector(i * externalIncr, -(myVerticalOuterTPBPanely / 2. + (y / (n + 1)) * a), -(heightOuterPlasticTPB / 2. + epsilonSiPM + SiPMHeight / 2.));

            SiPMsUpperExternalCapSector1.push_back(posSiPMUpperExternalCapSector1);
            SiPMsUpperExternalCapSector2.push_back(posSiPMUpperExternalCapSector2);
            SiPMsLowerExternalCapSector1.push_back(posSiPMLowerExternalCapSector1);
            SiPMsLowerExternalCapSector2.push_back(posSiPMLowerExternalCapSector2);

            placedSiPMsOnExternalCapSector++;
            placedSiPMsOnExternalCapSector++;
          }

          i++;
        }

        if (placedSiPMsOnExternalCapSector * 8 > nSiPMPerHalfSectorOutside * 16) {

          externalEpsilon = externalEpsilon + externalAdd;

          SiPMsUpperExternalCapSector1.clear();
          SiPMsUpperExternalCapSector2.clear();
          SiPMsLowerExternalCapSector1.clear();
          SiPMsLowerExternalCapSector2.clear();

          nSiPMPerColumn.clear();

        }

        else {

          recalculateOutside = false;

          //*****   WARNING: in non uniform case one may visualise geometry
          //configuration to check if all volumes are
          //                 correctly placed.

          if (SiPMUniformOutsideCaps == false) {

            int ncol = nSiPMPerColumn.size();

            // Parameters from fitting of parabola shape from non uniform
            // distribution on sides (assuming zero the linear term) with TPB on
            // SiPMs, free path 100 m and refl= 0.98.
            //
            //
            //
            //!!!!!Need to fit only until x value of plastic radius is
            //!reached!!!!!

            const G4double a = 0.004 * (1. / (cm * cm));
            const G4double c = 2418.56;
            const G4double min = 2180;

            const G4double columnsStep = (c - min) / (ncol + 1);

            G4int doneSiPM = 0;

            for (int n = 1; n <= ncol; n++) {

              for (int nn = 0; nn < nSiPMPerColumn.at(n - 1); nn++) {

                G4double xPosSiPM = sqrt(((n)*columnsStep) / a) * (1 - 0.9 / (n * n));

                if (n == 1) { xPosSiPM = sqrt(((n)*columnsStep) / a) * (0.6); }

                G4ThreeVector posSiPM1 = SiPMsUpperExternalCapSector1.at(doneSiPM + nn);
                G4ThreeVector posSiPM2 = SiPMsUpperExternalCapSector2.at(doneSiPM + nn);
                G4ThreeVector posSiPM3 = SiPMsLowerExternalCapSector1.at(doneSiPM + nn);
                G4ThreeVector posSiPM4 = SiPMsLowerExternalCapSector2.at(doneSiPM + nn);

                posSiPM1.setX(xPosSiPM);
                posSiPM2.setX(xPosSiPM);
                posSiPM3.setX(xPosSiPM);
                posSiPM4.setX(xPosSiPM);

                SiPMsUpperExternalCapSector1[doneSiPM + nn] = posSiPM1;
                SiPMsUpperExternalCapSector2[doneSiPM + nn] = posSiPM2;
                SiPMsLowerExternalCapSector1[doneSiPM + nn] = posSiPM3;
                SiPMsLowerExternalCapSector2[doneSiPM + nn] = posSiPM4;
              }

              doneSiPM += nSiPMPerColumn.at(n - 1);
            }
          }

          for (int k = 0; k < 8; k++) {

            G4RotationMatrix* rotSiPMOnUpperCapOutside = new G4RotationMatrix();

            rotSiPMOnUpperCapOutside->rotateZ(-k * ang - ang / 2.);  // Rotation i.o.t. have the side of SiPMs
                                                                     // parallel to outer plastic TPB edge

            rotSiPMOnUpperCapOutside->rotateY(180. * deg);  // External SiPMs on upper cap must be facing up

            G4RotationMatrix* rotSiPMOnLowerCapOutside = new G4RotationMatrix();

            rotSiPMOnLowerCapOutside->rotateZ(-k * ang - ang / 2.);  // Rotation i.o.t. have the side of SiPMs
                                                                     // parallel to outer plastic TPB edge

            externalUpperCapSiPMCopyNumber = 20000 + k * 1000 + 500;
            externalLowerCapSiPMCopyNumber = 20000 + k * 1000 + 400;

            for (int j = 0; j < placedSiPMsOnExternalCapSector / 2; j++) {

              G4ThreeVector pos1 = SiPMsUpperExternalCapSector1.at(j);
              pos1.rotateZ(k * ang + ang / 2.);

              G4ThreeVector pos2 = SiPMsUpperExternalCapSector2.at(j);
              pos2.rotateZ(k * ang + ang / 2.);

              G4ThreeVector pos3 = SiPMsLowerExternalCapSector1.at(j);
              pos3.rotateZ(k * ang + ang / 2.);

              G4ThreeVector pos4 = SiPMsLowerExternalCapSector2.at(j);
              pos4.rotateZ(k * ang + ang / 2.);

              physSiPMExternalUpperCap_1 = new G4PVPlacement(rotSiPMOnUpperCapOutside, pos1, myLogicShapeSiPMOnCap, "SiPMExternalUpperCap", myLogicShapeOuterLArBuffer, false, externalUpperCapSiPMCopyNumber, myCheckOverlaps);

              physSiPMInside.push_back(physSiPMExternalUpperCap_1);

              // Since one needs SiPMs coordinates respect to the world and not
              // respect to their mother volume (i.e. copper cage), and since
              //traslation along z axis is balaced, one needs only to de-rotate
              // since the copper cage is rotated respect to the world
              // (cryostat).
              outputSiPMs << externalUpperCapSiPMCopyNumber << " " << pos1.rotateZ(-ang / 2.) / cm << endl;

              externalUpperCapSiPMCopyNumber++;

              physSiPMExternalUpperCap_2 = new G4PVPlacement(rotSiPMOnUpperCapOutside, pos2, myLogicShapeSiPMOnCap, "SiPMExternalUpperCap", myLogicShapeOuterLArBuffer, false, externalUpperCapSiPMCopyNumber, myCheckOverlaps);

              physSiPMInside.push_back(physSiPMExternalUpperCap_2);

              // Since one needs SiPMs coordinates respect to the world and not
              // respect to their mother volume (i.e. copper cage), and since
              //traslation along z axis is balaced, one needs only to de-rotate
              // since the copper cage is rotated respect to the world
              // (cryostat).
              outputSiPMs << externalUpperCapSiPMCopyNumber << " " << pos2.rotateZ(-ang / 2.) / cm << endl;

              externalUpperCapSiPMCopyNumber++;

              physSiPMExternalLowerCap_1 = new G4PVPlacement(rotSiPMOnLowerCapOutside, pos3, myLogicShapeSiPMOnCap, "SiPMExternalLowerCap", myLogicShapeOuterLArBuffer, false, externalLowerCapSiPMCopyNumber, myCheckOverlaps);

              physSiPMInside.push_back(physSiPMExternalLowerCap_1);

              // Since one needs SiPMs coordinates respect to the world and not
              // respect to their mother volume (i.e. copper cage), and since
              //traslation along z axis is balaced, one needs only to de-rotate
              // since the copper cage is rotated respect to the world
              // (cryostat).
              outputSiPMs << externalLowerCapSiPMCopyNumber << " " << pos3.rotateZ(-ang / 2.) / cm << endl;

              externalLowerCapSiPMCopyNumber++;

              physSiPMExternalLowerCap_2 = new G4PVPlacement(rotSiPMOnLowerCapOutside, pos4, myLogicShapeSiPMOnCap, "SiPMExternalLowerCap", myLogicShapeOuterLArBuffer, false, externalLowerCapSiPMCopyNumber, myCheckOverlaps);

              physSiPMInside.push_back(physSiPMExternalLowerCap_2);

              // Since one needs SiPMs coordinates respect to the world and not
              // respect to their mother volume (i.e. copper cage), and since
              //traslation along z axis is balaced, one needs only to de-rotate
              // since the copper cage is rotated respect to the world
              // (cryostat).
              outputSiPMs << externalLowerCapSiPMCopyNumber << " " << pos4.rotateZ(-ang / 2.) / cm << endl;

              externalLowerCapSiPMCopyNumber++;
            }
          }

          nSiPMsEffective += placedSiPMsOnExternalCapSector * 8 * 2;
        }
      }

      cout << " " << endl;

      DSLog(routine) << "Activated SiPMs auto-placement" << endlog;

      cout << " " << endl;

      DSLog(routine) << "Required number of SiPMs in inner LAr buffer: " << nSiPMInside << endlog;
      DSLog(routine) << "Required number of SiPMs in outer LAr buffer: " << nSiPMOutside << endlog;

      cout << " " << endl;

      if (SiPMFractionInside != 0) {

        DSLog(routine) << "Required fraction of SiPMs on internal caps respect "
                          "to the total: "
                       << SiPMFractionInside << endlog;
      }

      if (SiPMFractionOutside != 0) {

        DSLog(routine) << "Required fraction of SiPMs on external caps respect "
                          "to the total: "
                       << SiPMFractionOutside << endlog;
      }

      cout << " " << endl;

      if (bestLightCollection) { DSLog(routine) << "Activated best light collection configuration for SiPMs " << endlog; }

      cout << " " << endl;

      DSLog(routine) << "TPB deposition on SiPMs: " << TPBOnSiPM << endl;

      cout << " " << endl;

      if (SiPMUniformInsideSides) {

        DSLog(routine) << "Uniform distribution on IAB sides" << endlog;

      }

      else {

        DSLog(routine) << "Non uniform distribution on IAB sides" << endlog;
      }

      if (SiPMUniformInsideCaps) {

        DSLog(routine) << "Uniform distribution on IAB caps" << endlog;

      }

      else {

        DSLog(routine) << "Non uniform distribution on IAB caps" << endlog;
      }

      if (SiPMUniformOutsideSides) {

        DSLog(routine) << "Uniform distribution on OAB sides" << endlog;

      }

      else {

        DSLog(routine) << "Non uniform distribution on OAB sides" << endlog;
      }

      if (SiPMUniformOutsideCaps) {

        DSLog(routine) << "Uniform distribution on OAB caps" << endlog;

      }

      else {

        DSLog(routine) << "Non uniform distribution on OAB caps" << endlog;
      }

      cout << " " << endl;

      DSLog(routine) << "Number of SiPMs placed on internal faces: " << (nSiPMHalfEdgeInside * nSiPMHeightInside) * 2 * 8 << endlog;

      DSLog(routine) << "Number of SiPMs placed on internal caps: " << (placedSiPMsOnInternalCapSector)*8 * 2 << endlog;

      DSLog(routine) << "Number of SiPMs placed on external faces: " << (nSiPMHalfEdgeOutside * nSiPMHeightOutside) * 2 * 8 << endlog;

      DSLog(routine) << "Number of SiPMs placed on external caps: " << (placedSiPMsOnExternalCapSector)*8 * 2 << endlog;

      DSLog(routine) << "Total number of SiPMs placed: " << nSiPMsEffective << endlog;

      cout << " " << endl;

      DSLog(routine) << "Fraction of SiPMs placed on internal caps respect to the total: " << (2. * placedSiPMsOnInternalCapSector * 8.) / (2. * placedSiPMsOnInternalCapSector * 8. + 2. * nSiPMHalfEdgeInside * nSiPMHeightInside * 8.) << endlog;

      DSLog(routine) << "Fraction of SiPMs placed on external caps respect to the total: " << (2. * placedSiPMsOnExternalCapSector * 8.) / (2. * placedSiPMsOnExternalCapSector * 8. + 2. * nSiPMHalfEdgeOutside * nSiPMHeightOutside * 8.) << endlog;

      cout << " " << endl;

    }  // End of autoplacement case

    /////////////////   SiPMs manual placement case  //////////////

    // Positions in input file must be in centimeters !!!

    else {

      G4int nSiPMInternalFaces = 0, nSiPMInternalCaps = 0, nSiPMExternalFaces = 0, nSiPMExternalCaps = 0;
      G4int SiPMsTotal = 0;

      G4int SiPMSector, SiPMPlace, SiPMBuffer;

      G4int SiPMCopyNumber;
      G4ThreeVector SiPMPosition, tmpSiPMPosition;

      while (inputSiPMs >> SiPMCopyNumber >> tmpSiPMPosition) {

        // FIXME
        // Error in output: istream ended before trying to input Hep3Vector
        //

        // Since one needs SiPMs coordinates respect to their mother volume and
        // not respect to the world and since traslation along z axis is
        // balaced, one needs only to rotate since the copper cage is rotated
        // respect to the world (cryostat).
        SiPMPosition = (tmpSiPMPosition * cm).rotateZ(ang / 2.);

        string SiPMManualName;

        SiPMsTotal++;

        SiPMBuffer = get_buffer(SiPMCopyNumber);  // 1 for inner buffer and 2 for outer buffer
        SiPMSector = get_sector(SiPMCopyNumber);  // From 0 to 7
        SiPMPlace = get_place(SiPMCopyNumber);    // 0 for faces, 4 for lower caps
                                                  // and 5 for upper caps

        G4RotationMatrix* rotSiPMManual = new G4RotationMatrix();
        rotSiPMManual->rotateZ(-SiPMSector * ang - ang / 2.);  // Rotation i.o.t. have the side of
                                                               // SiPMs parallel to LAr buffer edge

        G4LogicalVolume* SiPMManualMotherVolume;
        G4LogicalVolume* myLogicShapeSiPMManual;

        if (SiPMBuffer == 1) {
          SiPMManualName += "Internal";
          SiPMManualMotherVolume = myLogicShapeInnerLArBuffer;
          if (SiPMPlace == 4) {
            rotSiPMManual->rotateY(180. * deg);  // Internal SiPMs on lower cap must be facing up
            nSiPMInternalCaps++;
          } else if (SiPMPlace == 5) {
            nSiPMInternalCaps++;
          } else {
            nSiPMInternalFaces++;
          }

        } else {
          SiPMManualName += "External";
          SiPMManualMotherVolume = myLogicShapeOuterLArBuffer;
          if (SiPMPlace == 0) {
            rotSiPMManual->rotateY(180. * deg);  // SiPMs on outer faces must be facing outside
            nSiPMExternalFaces++;
          } else if (SiPMPlace == 5) {
            rotSiPMManual->rotateY(180. * deg);
            nSiPMExternalCaps++;
          } else {
            nSiPMExternalCaps++;
          }
        }

        if (SiPMPlace == 0) {
          SiPMManualName += "Side";
          myLogicShapeSiPMManual = myLogicShapeSiPMOnSide;
        } else {
          myLogicShapeSiPMManual = myLogicShapeSiPMOnCap;
          if (SiPMPlace == 4) {
            SiPMManualName += "LowerCap";
          } else {
            SiPMManualName += "UpperCap";
          }
        }

        physSiPMManual = new G4PVPlacement(rotSiPMManual, SiPMPosition, myLogicShapeSiPMManual, "SiPM" + SiPMManualName, SiPMManualMotherVolume, false, SiPMCopyNumber, myCheckOverlaps);

      }  // End of while loop

      DSLog(routine) << "Activated SiPMs manual placement from file: " << inputFilename << endlog;

      DSLog(routine) << "TPB deposition on SiPMs: " << TPBOnSiPM << endl;

      cout << " " << endl;

      DSLog(routine) << "Number of SiPMs placed on internal faces: " << nSiPMInternalFaces << endlog;

      DSLog(routine) << "Number of SiPMs placed on internal caps: " << nSiPMInternalCaps << endlog;

      DSLog(routine) << "Number of SiPMs placed on external faces: " << nSiPMExternalFaces << endlog;

      DSLog(routine) << "Number of SiPMs placed on external caps: " << nSiPMExternalCaps << endlog;

      DSLog(routine) << "Total number of SiPMs placed: " << SiPMsTotal << endlog;

      cout << " " << endl;

      inputSiPMs.close();

    }  // End of manual placement case

    outputSiPMs.close();

  }  // End of SiPMs placement, if SiPMs boolean variable is true

  DefineSurfaces();
}

DSDetectorPlasticVeto::~DSDetectorPlasticVeto() {
  ;  // delete fMessenger;
}

///////////////////////////////////////////////////////////////////////
////////////      REFLECTING SURFACES     /////////////////////////////
///////////////////////////////////////////////////////////////////////

void DSDetectorPlasticVeto::DefineSurfaces() {

  ////////////////////////////////////////
  // TPB --> teflon (Reflector)
  //  Should not be teflon --> TPB as surface is defined as dielectric_metal.
  ////////////////////////////////////////

  G4double TeflonTPBENE[4] = {0.1 * eV, 8.0 * eV, 8.3 * eV, 20.0 * eV};
  G4double TREFUV = DSParameters::Get()->GetTeflonTPBUVRef();
  G4double TREFVIS = DSParameters::Get()->GetTeflonTPBVisRef();
  G4double TeflonTPBREF[4] = {TREFVIS, TREFVIS, TREFUV, TREFUV};

  G4OpticalSurface* fOpTPBTeflonSurface = new G4OpticalSurface("OpTBPTeflonSurface");

  new G4LogicalBorderSurface("TPBTeflonSurface", physLowerInnerTPBCylinder, physLowerInnerAcrylicCylinder, fOpTPBTeflonSurface);

  new G4LogicalBorderSurface("TPBTeflonSurface", physUpperInnerTPBCylinder, physUpperInnerAcrylicCylinder, fOpTPBTeflonSurface);

  new G4LogicalBorderSurface("TPBTeflonSurface", physLowerOuterTPBCylinder, physLowerOuterAcrylicCylinder, fOpTPBTeflonSurface);

  new G4LogicalBorderSurface("TPBTeflonSurface", physUpperOuterTPBCylinder, physUpperOuterAcrylicCylinder, fOpTPBTeflonSurface);

  for (int i = 0; i < 8; i++) new G4LogicalBorderSurface("TPBTeflonSurface", physTPBInnerPanelFace[i], physAcrylicInnerPanelFace[i], fOpTPBTeflonSurface);

  for (int i = 0; i < 8; i++) new G4LogicalBorderSurface("TPBTeflonSurface", physTPBInnerPanelEdge[i], physAcrylicInnerPanelEdge[i], fOpTPBTeflonSurface);

  for (int i = 0; i < 8; i++) new G4LogicalBorderSurface("TPBTeflonSurface", physTPBOuterPanelFace[i], physAcrylicOuterPanelFace[i], fOpTPBTeflonSurface);

  for (int i = 0; i < 8; i++) new G4LogicalBorderSurface("TPBTeflonSurface", physTPBOuterPanelEdge[i], physAcrylicOuterPanelEdge[i], fOpTPBTeflonSurface);

  new G4LogicalBorderSurface("TPBTeflonSurface", physInnerPlasticTPB, physPlastic, fOpTPBTeflonSurface);

  new G4LogicalBorderSurface("TPBTeflonSurface", physOuterPlasticTPB, physPlastic, fOpTPBTeflonSurface);

  new G4LogicalBorderSurface("TPBTeflonSurface", physExternalTPB, physCopper, fOpTPBTeflonSurface);

  fOpTPBTeflonSurface->SetType(dielectric_metal);

  fOpTPBTeflonSurface->SetModel(unified);

  // fOpTPBTeflonSurface->SetFinish(groundfrontpainted);
  // fOpTPBTeflonSurface->SetFinish(ground);
  // fOpTPBTeflonSurface->SetSigmaAlpha(0.1);

  fOpTPBTeflonSurface->SetFinish(polished);

  G4MaterialPropertiesTable* fTPBTeflonSurfProp = new G4MaterialPropertiesTable();

  fTPBTeflonSurfProp->AddProperty("REFLECTIVITY", TeflonTPBENE, TeflonTPBREF, 4);

  fOpTPBTeflonSurface->SetMaterialPropertiesTable(fTPBTeflonSurfProp);
}

/////////  Functions for manual placement of SiPMs from file  /////////////

// Returns SiPM sector
G4int DSDetectorPlasticVeto::get_sector(G4int SiPMCopyNumber) {

  G4int sector;

  if (SiPMCopyNumber < 20000) {
    sector = floor((SiPMCopyNumber / 1000) - 10);
  } else {
    sector = floor((SiPMCopyNumber / 1000) - 20);
  }

  return sector;
}

// Returns 1 for inner buffer and 2 for outer buffer
G4int DSDetectorPlasticVeto::get_buffer(G4int SiPMCopyNumber) {

  G4int buffer = ceil(SiPMCopyNumber / 20000) + 1;

  return buffer;
}

// Returns 0 for faces, 4 for lower caps and 5 for upper caps
G4int DSDetectorPlasticVeto::get_place(G4int SiPMCopyNumber) {

  G4int place, tmpPlace;

  G4int buffer = ceil(SiPMCopyNumber / 20000) + 1;

  G4int sector;

  if (SiPMCopyNumber < 20000) {
    sector = floor(SiPMCopyNumber / 1000) - 10;
  } else {
    sector = floor(SiPMCopyNumber / 1000) - 20;
  }

  tmpPlace = (SiPMCopyNumber - buffer * 10000 - sector * 1000);

  if (tmpPlace < 400) {
    place = 0;
  }

  else if (tmpPlace >= 400 && tmpPlace < 500) {
    place = 4;
  }

  else {
    place = 5;
  }

  return place;
}
