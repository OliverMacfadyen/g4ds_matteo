#include <iostream>
#include "G4Box.hh"
#include "G4Orb.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
#include "G4ThreeVector.hh"

#include "DSDetectorProtoProto.hh"
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

// 25th January 2019
// based on drawing shown at https://agenda.infn.it/event/17949/

DSDetectorProtoProto::DSDetectorProtoProto(G4VPhysicalVolume* myMotherVolume) {

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
  DSLog(routine) << " Constructing DSProtoProto Geometry" << endlog;

  // const double myTwoPi = 2 * M_PI * rad;

  G4bool myCheckOverlap = DSStorage::Get()->GetCheckOverlap();
  G4ThreeVector myZeros(0., 0., 0.);

  //-------------------------//
  //      LAr bath           //
  //-------------------------//
  G4Box* solid_LArBath = new G4Box("solid_LArBath", 30 * cm, 30 * cm, 50 * cm);
  G4LogicalVolume* flogic_LArBath = new G4LogicalVolume(solid_LArBath, DSMaterial::Get()->GetNSLiquidArgon(), "LArBath_Logic");
  fPhysic_LArBath = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), "LArBath_Physic", flogic_LArBath, fMotherVolume, false, 0, myCheckOverlap);

  // Parameters from drawings
  G4double myTPCHeight = DSStorage::Get()->GetDS20kTPCheight();  // default 120
                                                                 // cm
  G4double myTPCEdge = DSStorage::Get()->GetDS20kTPCedge();      // default 240 cm
  // G4double myPMDSuppEdge = myTPCEdge + 5 * cm;

  // Reflector
  double myReflectorThickness = DSStorage::Get()->GetDS20kTPCVesselThickness();

  double myGasPocketThickness = 0.7 * cm;

  // place the center of the TPC in the center of the ref system
  double myLArGArBoundaryPosZ = DSStorage::Get()->GetDS20kTPCheight() / 2. - myGasPocketThickness;
  // Set the z coordinate of the LAr - GAr interface, necessary for S2
  // generation in DSLightX
  DSStorage::Get()->SetLArGArBoundaryPosZ(myLArGArBoundaryPosZ + 1.0 * um);

  // Other parameters
  // double myTPBThickness = 0.1 * mm;
  // double myPArThickness = 0. * mm;
  // double mySiPMOffset = DSStorage::Get()->GetSiPMOffset();  // default: 5 cm
  // double mySiPmThickness = 1. * mm;
  // double mySiPmBoardThickness = 5.0 * mm;
  // double myAcrylicBoardThickness = 6.0 * mm;
  // double myTopAcrylicThickness = 1.5 * cm;
  double myActiveGasLogicThickness = myGasPocketThickness;  // + myWindowsThickness + mySiPmThickness +
                                                            // mySiPmBoardThickness;

  if (DSStorage::Get()->GetGdLayerThickness() > 1 * mm) myActiveGasLogicThickness = DSStorage::Get()->GetGdLayerThickness();

  //-------------------------//
  //      Sealed Vessel      //
  //-------------------------//

  G4double vessel_external_x = myTPCEdge + 2 * myReflectorThickness;
  G4double vessel_external_z = myTPCHeight + 2 * myReflectorThickness;
  G4Box* solid_vessel = new G4Box("solid_vessel", vessel_external_x / 2., vessel_external_x / 2., vessel_external_z / 2.);
  G4LogicalVolume* flogic_vessel = new G4LogicalVolume(solid_vessel, DSMaterial::Get()->GetAcrylic(), "Vessel_Logic");
  fPhysic_vessel = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), "Vessel_Physic", flogic_vessel, fPhysic_LArBath, false, 0, myCheckOverlap);

  //-------------------------//
  //          TPC            //
  //-------------------------//

  //-------------------------//
  //         LAr             //
  G4double LAr_external_x = myTPCEdge;
  G4double LAr_external_z = myTPCHeight;
  G4Box* solid_LAr = new G4Box("solid_LAr", LAr_external_x / 2., LAr_external_x / 2., LAr_external_z / 2.);
  G4LogicalVolume* flogic_LAr = new G4LogicalVolume(solid_LAr, DSMaterial::Get()->GetLiquidArgon(), "LAr_Logic");
  fPhysic_LAr = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), "ActiveLAr", flogic_LAr, fPhysic_vessel, false, 0, myCheckOverlap);

  flogic_LAr->SetVisAttributes(myBlue);
  G4Region* fLArRegion = new G4Region("LAr_Logic");
  flogic_LAr->SetRegion(fLArRegion);
  fLArRegion->AddRootLogicalVolume(flogic_LAr);

  DSStorage::Get()->SetLiquidArgonIndex((int)flogic_LAr->GetMaterial()->GetIndex());

  //-------------------------//
  //         GAr             //
  G4double GAr_external_x = myTPCEdge;
  G4double GAr_external_z = myActiveGasLogicThickness;
  G4Box* solid_GAr = new G4Box("solid_GAr", GAr_external_x / 2., GAr_external_x / 2., GAr_external_z / 2.);
  G4LogicalVolume* flogic_GAr = new G4LogicalVolume(solid_GAr, DSMaterial::Get()->GetGaseousArgon(), "GAr_Logic");
  fPhysic_GAr = new G4PVPlacement(0, G4ThreeVector(0, 0, myTPCHeight / 2. - myActiveGasLogicThickness / 2.), "GAr_Physic", flogic_GAr, fPhysic_LAr, false, 0, myCheckOverlap);

  //-------------------------//
  //         TPB             //
  G4double TPB_external_x = myTPCEdge;
  G4double TPB_external_z = 100 * um;
  G4Box* solid_TPB = new G4Box("solid_TPB", TPB_external_x / 2., TPB_external_x / 2., TPB_external_z / 2.);
  G4LogicalVolume* flogic_TPB = new G4LogicalVolume(solid_TPB, DSMaterial::Get()->GetTPB(), "TPB_Logic");
  fPhysic_TPB = new G4PVPlacement(0, G4ThreeVector(0, 0, myActiveGasLogicThickness / 2. - TPB_external_z / 2.), "TPB_Physic", flogic_TPB, fPhysic_GAr, false, 0, myCheckOverlap);

  //---------------------//
  //      SiPM array     //
  //---------------------//

  // Top Array
  G4Box* fSolidSiPMTop = new G4Box("SolidSiPMTop", LAr_external_x / 2., LAr_external_x / 2., 0.5 * mm);
  G4LogicalVolume* fLogicSiPMTop = new G4LogicalVolume(fSolidSiPMTop, DSMaterial::Get()->GetMetalSilicon(), "SiPMTop_Logic");

  G4double vertical_coordinate = myTPCHeight / 2. + myReflectorThickness + 0.5 * mm + DSStorage::Get()->GetSiPMOffset() + 1 * mm;  // set minimum distance of 0.1 mm
  fPhysic_SiPmTop = new G4PVPlacement(0, G4ThreeVector(0, 0, vertical_coordinate), "SiPMTop", fLogicSiPMTop, fPhysic_LArBath, false, -56, myCheckOverlap);

  // Bottom Array
  if (DSStorage::Get()->GetIsProtoProtoBottomSiPM()) {
    G4Box* fSolidSiPMBottom = new G4Box("SolidSiPMBottom", LAr_external_x / 2., LAr_external_x / 2., 0.5 * mm);
    G4LogicalVolume* fLogicSiPMBottom = new G4LogicalVolume(fSolidSiPMBottom, DSMaterial::Get()->GetMetalSilicon(), "SiPMBottom_Logic");

    vertical_coordinate = -(myTPCHeight / 2. + myReflectorThickness + 0.5 * mm + DSStorage::Get()->GetSiPMOffset() + 1 * mm);  // set minimum distance of 0.1 mm
    fPhysic_SiPMBottom = new G4PVPlacement(0, G4ThreeVector(0, 0, vertical_coordinate), "SiPMBottom", fLogicSiPMBottom, fPhysic_LArBath, false, -56, myCheckOverlap);
  }

  // make SiPM as pe storing material
  DSStorage::Get()->SetPMTMaterialIndex(fLogicSiPMTop->GetMaterial()->GetIndex());

  DefineSurfaces();
}

DSDetectorProtoProto::~DSDetectorProtoProto() {
  ;  // delete fMessenger;
}

void DSDetectorProtoProto::DefineSurfaces() {

  ////////////////////////////////////////
  // TPB <--> GAr and TPB <--> LAr
  G4OpticalSurface* fOpTPBGArSurface = new G4OpticalSurface("OpTPBGArSurface");
  new G4LogicalBorderSurface("GArTPBSurface", fPhysic_GAr, fPhysic_TPB, fOpTPBGArSurface);
  new G4LogicalBorderSurface("TPBGArSurface", fPhysic_TPB, fPhysic_GAr, fOpTPBGArSurface);
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
  // Copied from DS DetectorDS50 - 24th april 2015
  ////////////////////////////////////////
  /*
  // TPC - BScint
  // warning: this only works if the current 20k design with LSV is constructed
  G4OpticalSurface *fOpElectropolishedStainlessSteelSurface = new
G4OpticalSurface("OpElectropolishedStainlessSteelSurface");
  fOpElectropolishedStainlessSteelSurface->SetType(dielectric_metal);
  fOpElectropolishedStainlessSteelSurface->SetModel(unified);
  fOpElectropolishedStainlessSteelSurface->SetFinish(groundbackpainted);
  fOpElectropolishedStainlessSteelSurface->SetMaterialPropertiesTable(DSMaterial::Get()->GetElectropolishedStainlessSteelMPT());
  fDS20kOuterSurface = new G4LogicalBorderSurface("DS20kOuterSurface",
fMotherVolume, fPhysicDS20k, fOpElectropolishedStainlessSteelSurface);

  // LAR-Grid
  G4OpticalSurface *fOpGridLArSurface = new
G4OpticalSurface("OpGridLArSurface"); new
G4LogicalBorderSurface("GridLArSurface", fPhysicActiveLAr, fPhysicGrid,
fOpGridLArSurface); fOpGridLArSurface->SetType( dielectric_dielectric );
  fOpGridLArSurface->SetModel( glisur );
  fOpGridLArSurface->SetFinish( polished );
  G4double GridLArENE[4] = {0.1*eV, 8.0*eV, 8.3*eV, 20.0*eV};
  G4double GRUV = DSParameters::Get()->GetLArGridUVRef();
  G4double GRVIS = DSParameters::Get()->GetLArGridVisRef();
  G4double GridLArREF[4] = {GRVIS, GRVIS, GRUV, GRUV};
  G4MaterialPropertiesTable *fGridLArSurfProp = new G4MaterialPropertiesTable();
  if (DSParameters::Get()->GetWithNewGridModel())
    // the grid model is described in DSStorage.cc, the surface is treated in
G4OpBoundaryProcess.cc fGridLArSurfProp->AddConstProperty("DOGRID",1);
  // Now use the following in old and new models.  By G4 convention,
"reflectivity" is actually 1-absorption.
  fGridLArSurfProp->AddProperty("REFLECTIVITY", GridLArENE, GridLArREF, 4);
  fOpGridLArSurface->SetMaterialPropertiesTable( fGridLArSurfProp );

  // Grid->LAR
  if (DSParameters::Get()->GetWithNewGridModel()) {
    G4OpticalSurface *fOpLArGridSurface = new
G4OpticalSurface("OpLArGridSurface"); new
G4LogicalBorderSurface("LArGridSurface", fPhysicGrid, fPhysicActiveLAr,
fOpLArGridSurface); cout << " With DS50 new grid model " << endl ;
    fOpLArGridSurface->SetType( dielectric_dielectric );
  //  fOpLArGridSurface->SetModel( glisur );
  //  fOpLArGridSurface->SetFinish( polished );
    G4MaterialPropertiesTable *fLArGridSurfProp = new
G4MaterialPropertiesTable();
    // the grid model is described in DSStorage.cc, the surface is treated in
G4OpBoundaryProcess.cc fLArGridSurfProp->AddConstProperty("DOGRIDEXIT",1);
    fOpLArGridSurface->SetMaterialPropertiesTable( fLArGridSurfProp );
  }

  ////////////////////////////////////////
  // TPB <--> GAr and TPB <--> LAr
  G4OpticalSurface *fOpTPBGArSurface = new G4OpticalSurface("OpTPBGArSurface");
  new G4LogicalBorderSurface("GArTPBSurface", fPhysicGasPocket, fPhysicTPBTop,
fOpTPBGArSurface ); new G4LogicalBorderSurface("TPBGArSurface", fPhysicTPBTop,
fPhysicGasPocket, fOpTPBGArSurface ); new
G4LogicalBorderSurface("LArTPBSurfaceSide", fPhysicActiveLAr, fPhysicTPBSide,
fOpTPBGArSurface ); new G4LogicalBorderSurface("TPBLArSurfaceSide",
fPhysicTPBSide, fPhysicActiveLAr, fOpTPBGArSurface );
//  new G4LogicalBorderSurface("LArTPBSurfaceBot", fPhysicActiveLAr,
fPhysicTPBBottom, fOpTPBGArSurface );
//  new G4LogicalBorderSurface("TPBLArSurfaceBot", fPhysicTPBBottom,
fPhysicActiveLAr, fOpTPBGArSurface ); new
G4LogicalBorderSurface("LArLayerTPBSurface", fPhysicLArLayer, fPhysicTPBTop,
fOpTPBGArSurface ); new G4LogicalBorderSurface("TPBLArLayerSurface",
fPhysicTPBTop, fPhysicLArLayer, fOpTPBGArSurface ); fOpTPBGArSurface->SetType(
dielectric_dielectric ); fOpTPBGArSurface->SetModel( unified );
  fOpTPBGArSurface->SetFinish( ground );
  fOpTPBGArSurface->SetSigmaAlpha(0.3);

  G4double VISTRAN = DSParameters::Get()->GetArTPBVisTran();

  const G4int NUM = 4;
  G4double pp[NUM] = {0.1*eV, 8.0*eV, 8.3*eV, 20.0*eV};  // vis, vis, UV, UV
  G4double specularlobe[NUM] = {0., 0., 0., 0.};         //--
  G4double specularspike[NUM] = {0., 0., 0., 0.};        //----  gives all
reflection to Lambertian lobe G4double backscatter[NUM] = {0., 0., 0., 0.}; //--
  G4double reflectivity[NUM] = {1.0, 1.0, 1.0, 1.0};     //  To set 1-absorption
  G4double transmitivity[NUM] = {VISTRAN, VISTRAN, 1.0, 1.0};    //  To set
reflection vs. transmission, overridding Fresnel
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

  ////////////////////////////////////////
  // ITO /////
  ////////////////////////////////////////

  G4MaterialPropertiesTable *fITOSurfProp = new G4MaterialPropertiesTable();
  if ( DSParameters::Get()->GetWithITO()  )
fITOSurfProp->AddConstProperty("DOITO",1);




  ////////////////////////////////////////
  // BellTop (acrylic) <--> TPB and CathodeWindow <--> TPB (both with ITO)
  // In the current model, the diffuse nature of the TPB is handled
  // entirely at the TPB-GAr/LAr surface, not here.
  // Make this bi-directional.
  ////////////////////////////////////////
  G4OpticalSurface *fOpWindowTPBSurface     = new
G4OpticalSurface("OpWindowTPBSurface"); new
G4LogicalBorderSurface("TopWindowTPBSurface",   fPhysicTopWindow ,
fPhysicTPBTop, fOpWindowTPBSurface ); new
G4LogicalBorderSurface("BottomWindowTPBSurface", fPhysicBotWindow,
fPhysicTPBSide, fOpWindowTPBSurface ); new
G4LogicalBorderSurface("TPBTopWindowSurface",    fPhysicTPBTop,fPhysicTopWindow
,  fOpWindowTPBSurface ); new G4LogicalBorderSurface("TPBBottomWindowSurface",
fPhysicTPBSide, fPhysicBotWindow, fOpWindowTPBSurface );
  fOpWindowTPBSurface->SetType( dielectric_dielectric );
  fOpWindowTPBSurface->SetModel( unified );
  //  fOpWindowTPBSurface->SetFinish( ground );
  fOpWindowTPBSurface->SetFinish( polished );
  //G4MaterialPropertiesTable *fWindowTPBSurfProp = new
G4MaterialPropertiesTable(); fOpWindowTPBSurface->SetMaterialPropertiesTable(
fITOSurfProp );

*/

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

  ////////////////////////////////////////
  // LAr and CathodeWindow <--> LAr (with ITO)    2017-02-16
  // Make this bi-directional.
  ////////////////////////////////////////
  G4OpticalSurface *fOpBotWindowArSurface     = new
  G4OpticalSurface("OpBotWindowArSurface"); new
  G4LogicalBorderSurface("BottomWindowLArSurfaceIn",fPhysicInactiveLar,
  fPhysicBotWindow,fOpBotWindowArSurface ); new
  G4LogicalBorderSurface("BottomWindowLArSurfaceOut", fPhysicBotWindow,
  fPhysicInactiveLar, fOpBotWindowArSurface ); fOpBotWindowArSurface->SetType(
  dielectric_dielectric ); fOpBotWindowArSurface->SetModel( unified );
  //  fOpBotWindowArSurface->SetFinish( ground );
  fOpBotWindowArSurface->SetFinish( polished );
  //G4MaterialPropertiesTable *fWindowArSurfProp = new
  G4MaterialPropertiesTable();
   fOpBotWindowArSurface->SetMaterialPropertiesTable( fITOSurfProp );     //
  With ITO HERE


  */

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

  ///////////////
  // TPB --> Teflon (Reflector)
  //  Should be no Teflon --> TPB as surface is defined as dielectric_metal.
  ////////////////////////////////////////
  G4double TeflonTPBENE[4] = {0.1*eV, 8.0*eV, 8.3*eV, 20.0*eV};
  G4double TREFUV  = DSParameters::Get()->GetTeflonTPBUVRef();
  G4double TREFVIS = DSParameters::Get()->GetTeflonTPBVisRef();
  G4double TeflonTPBREF[4] = {TREFVIS, TREFVIS ,TREFUV , TREFUV };
  G4OpticalSurface *fOpTPBTeflonSurface = new
  G4OpticalSurface("OpTBPTeflonSurface"); new
  G4LogicalBorderSurface("TPBTeflonSurface", fPhysicTPBSide,fPhysicTeflonBottom
  , fOpTPBTeflonSurface ); fOpTPBTeflonSurface->SetType( dielectric_metal );
  fOpTPBTeflonSurface->SetModel(unified);

  // PDM: though I can't see how, the following settings are giving the desired
  Lambertian reflection fOpTPBTeflonSurface->SetFinish(groundfrontpainted);
  //fOpTPBTeflonSurface->SetFinish(ground);
  fOpTPBTeflonSurface->SetSigmaAlpha(0.1);

  G4MaterialPropertiesTable *fTPBTeflonSurfProp = new
  G4MaterialPropertiesTable(); fTPBTeflonSurfProp->AddProperty("REFLECTIVITY",
  TeflonTPBENE, TeflonTPBREF, 4);
  fOpTPBTeflonSurface->SetMaterialPropertiesTable( fTPBTeflonSurfProp );

  */
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
  fTPBTeflonSurfProp2 ); new G4LogicalBorderSurface("SSteelOuterSurface",
                                  fMotherVolume,
                                fPhysicDS20kVac,
           fOpElectropolishedStainlessSteelSurface);

 */
}

/*
double DSDetectorProtoProto::GetOctagonInnerRadius(double edge) {  return edge
/2 *(1+sqrt(2)); } double DSDetectorProtoProto::GetOctagonOuterRadius(double
edge) {  return edge /2. *sqrt(4  + 2 *sqrt(2)) ; }

PointColPtr DSDetectorProtoProto::createGeometry(double r0, double hc, double
z0, double hb, double ht, double offset=0, double mytopoff=200, double
mybotoff=150) {
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

  //float mybotoff = 150 ;
  points->push_back(PairFF(z0-hc2-hb-offset-mybotoff,0));
  for ( int i=1; i<nb; ++i ) {
    const float angle = float(i) * 90.0 / float(nb);
    const float rads = M_PI / 180.0 * angle;
    float z = z0 - hc2 - hb * std::cos(rads)  - mybotoff ;
    float r =            r0 * std::sin(rads);

    if ( offset > 0 ) {
      // normal vector
      float n_z = -r0 * std::cos(rads);
      float n_r =  hb * std::sin(rads);
      const float l = std::sqrt(n_z*n_z+n_r*n_r);
      // normal vector with length "offset"
      n_z *= offset / l;
      n_r *= offset / l;
      z += n_z;
      r += n_r;
    }
    points->push_back(PairFF(z,r));
  }

  // cylinder
  for ( int i=0; i<nc; ++i ) {
    // normal in radial direction, offset in radius
    const float z = z0 + hc * ( float(i)/float(nc-1) - 0.5 );
    const float r = r0 + offset;
    points->push_back(PairFF(z,r));
  }

  // top cap
  //float mytopoff = 200 ;
  for (int i=0;i< mytopoff/10; ++i) {
    float z = z0 + hc2 +  5 + i*10 ;
    float r = r0 + offset;
    points->push_back(PairFF(z,r));
  }
  for ( int i=1; i<nt; ++i ) {
    const float angle = float(i) * 90.0 / float(nb);
    const float rads = M_PI / 180.0 * angle;
    float z = z0 + hc2 + ht * std::sin(rads) + mytopoff ;
    float r =            r0 * std::cos(rads);
    if ( offset > 0 ) {
      // normal vector
      float n_z = r0 * std::sin(rads);
      float n_r = ht * std::cos(rads);
      const float l = std::sqrt(n_z*n_z+n_r*n_r);
      // normal vector with length "offset"
      n_z *= offset / l;
      n_r *= offset / l;
      z += n_z;
      r += n_r;
    }
    points->push_back(PairFF(z,r));
  }
  // normal in +z direction, offset in direction +z
  points->push_back(PairFF(z0+hc2+ht+offset+mytopoff,0));

  return points;
}
*/
