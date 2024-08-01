#include "DSDetectorDS20kNeutronVeto.hh"
#include <fstream>
#include <iostream>
#include "DSEventHandler.hh"
#include "DSIO.hh"
#include "DSLogger.hh"
#include "DSMaterial.hh"
#include "DSParameters.hh"
#include "DSStorage.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalConstants.hh"
#include "G4Polycone.hh"
#include "G4Polyhedra.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4VisAttributes.hh"

using namespace std;

DSDetectorDS20kNeutronVeto::DSDetectorDS20kNeutronVeto(G4VPhysicalVolume* myMotherVolume) {
  // Apr 2022 PA. Merged to master, renamed and cleaned
  // 2021 AC: start from the G3 geometry to implement the planc geometry
  // Oct 2017. Finally commit the acrylic vessel design. default configuration
  // is still spherical veto

  fMotherVolume = myMotherVolume;
  // G4int volumeNumber = 10000;

  DSLog(routine) << " Constructing DS20k NeutronVeto Geometry" << endlog;

  const G4double myTwoPi = 2 * M_PI * rad;
  // const G4double myIDiameter         = DSStorage::Get()->GetIDiameter()  ;
  // //myAcrylicDiameter is the parameter that allows to distinguish between old
  // (no vessel) and new (with vessel) design
  // //the dimensions of the vessel are now (Jul'17) fixed at the current
  // proposal (3 m radius, 6.4 m height).
  // //The asymmetric cryostat is placed at the center of the vessel
  // const G4double myAcrylicDiameter   = 0 ; //
  // DSStorage::Get()->GetAcrylicVesselDiameter()  ;

  G4RotationMatrix* myDefaultRotation2 = new G4RotationMatrix();
  G4double myDefaultAngle = myTwoPi / 16.;
  myDefaultRotation2->rotateZ(-myDefaultAngle);  // to correct for octagonal symmetry

  G4bool myCheckOverlap = DSStorage::Get()->GetCheckOverlap();
  G4ThreeVector myZeros(0., 0., 0.);

  // ------------------------------------- //
  // ------   Vessel  dimensions ------ //
  // ------------------------------------- //

  // PA this part is taken from an old implementation of a double wall cryostat.
  // Only one membrane is used here

  // Vessel material
  G4Material* myCryoMat = DSMaterial::Get()->GetStainlessSteel();
  if (DSStorage::Get()->GetDS20kCryoMaterial() == 1) myCryoMat = DSMaterial::Get()->GetMetalTitanium();
  else if (DSStorage::Get()->GetDS20kCryoMaterial() == 2)
    myCryoMat = DSMaterial::Get()->GetMetalCopperCryo();

  //  const G4double myPlasticThickness   =
  //  DSStorage::Get()->GetDS20kGdPlasticThickness();
  const G4double myPlasticThickness = DSStorage::Get()->GetDS20kWLSPENThickness();
  // double myCryoToCornerDistance =
  // DSStorage::Get()->GetDS20kCryoCornerDistance();  //default 5 cm.

  // double myLArGArBoundaryPosZ = DSStorage::Get()->GetDS20kTPCheight() / 2.;
  //  Parameters from drawings
  // G4double myTPCHeight = DSStorage::Get()->GetDS20kTPCheight();  //default
  // 120 cm
  G4double myTPCEdge = DSStorage::Get()->GetDS20kTPCedge();  // default 240 cm

  double myCryoWallThickness = DSStorage::Get()->GetDS20kCryoWallThick();
  double myCryoWallThicknessI = DSStorage::Get()->GetDS20kCryoWallThick();
  double myCryoCapWallThickness = DSStorage::Get()->GetDS20kCryoCapWallThick();
  // G4cout<<" myCryoWallThickness "<<myCryoWallThickness<<"
  // myCryoWallThicknessI "<<myCryoWallThicknessI<<G4endl;
  if (myCryoWallThickness == 0 * cm) {
    myCryoWallThicknessI = 1.25 * cm;
    myCryoWallThickness = 1.75 * cm;
  }
  G4double myLArInterfaceThickness = DSStorage::Get()->GetDS20kWLSLArThickness();

  G4double myInnerCryostatZ[300];
  G4double myInnerCryostatRout[300];

  G4double myPlasticZ[300];
  G4double myPlasticRout[300];
  //G4double myPlasticUllageZ[300];
  //G4double myPlasticUllageRout[300];

  G4double myLiqArgonZ[300];
  G4double myLiqArgonRout[300];
  G4double myLiqArgonZI[300];
  G4double myLiqArgonRoutI[300];

  // Ullage
  //G4double myUllageZ[300];
  G4double myUllageZI[300];
  G4double myUllageRoutI[300];
  //G4double myUllageRout[300];

  G4double myRmin[300];
  for (int ii = 0; ii < 300; ++ii) myRmin[ii] = 0.;

  // G4int myNumPointsOuterCryo = 0;
  // G4int myNumPointsVacCryo = 0;
  G4int myNumPointsInnerCryo = 0;
  // G4int myNumPointsLarGAs_tot = 0;
  // G4int myNumPointsGasArgon = 0;
  G4int myNumPointsLiqArgonI = 0;
  G4int myNumPointsLiqArgon = 0;
  G4int myNumPointsPlastic = 0;
  // G4int myNumPointsTopCryoFiller = 0;
  G4int myNumPointsUllageI = 0;
  //G4int myNumPointsPlasticUllage = 0;
  //myNumPointsUllageG4int myNumPointsUllage = 0;
  // G4double maxr;//, maxz;

  // 2017-March all the hardcoded dimensions of the cryostat are changed with
  // parameters related to the size of the TPC, so that the cryostat scales when
  // DSProto is simulated. 2017-Oct removed reading of cryostat profile,
  // switched to function double myOuterCryostatDistance  = .1958 * myTPCEdge
  // +myGdLayerThickness*1.04 + myCryoToCornerDistance;   // must be larger
  // than 23.5 cm

  double myDefaultR = DSStorage::Get()->GetDS20kCryoR();  // 1880 ; //in mm without the unitmyRingOuterRadius
                                                          // + myOuterCryostatDistance ;
  if (myTPCEdge / cm < 50)
    myDefaultR = 1.75 * myTPCEdge / mm;  // 1880 ; //in mm without the unitmyRingOuterRadius   +

  double myBottomCapH = DSStorage::Get()->GetDS20kCryoBottomCap();
  double myTopCapH = DSStorage::Get()->GetDS20kCryoTopCap();
  double myDefaultRI = myDefaultR;

  double myTopOffset = DSStorage::Get()->GetDS20kCryoTopOffset();
  double myBotOffset = DSStorage::Get()->GetDS20kCryoBottomOffset();
  double myTCryostatBarrelH = DSStorage::Get()->GetDS20KCryoBarrelH();

  //set here the level of liquid UAr
  double myUllageLevel = 2350*mm;     // level calculated to have ~4% of gas with the current vessel geometry 


  PointColPtr pointCol;
  pointCol = DSDetectorDS20kNeutronVeto::createGeometry(myDefaultRI , myTCryostatBarrelH, 0.,myBottomCapH, myTopCapH, myTopOffset, myBotOffset,true,0);
  for ( unsigned int i=0; i<pointCol->size(); ++i ) {
    myInnerCryostatZ[i]    = pointCol->at(i).first  *mm ;
    myInnerCryostatRout[i] = pointCol->at(i).second *mm ;
  }
  myNumPointsInnerCryo = pointCol->size();

  // for (int i = 0; i<myNumPointsInnerCryo; i++)
  //         G4cout<<" Cryo : "<<myInnerCryostatZ[i]<<" "<<myInnerCryostatRout[i]<<G4endl;

  unsigned int temp_j =0;
  unsigned int temp_k =0;
  //LAr between PEN and ESR
  pointCol = DSDetectorDS20kNeutronVeto::createGeometry(myDefaultRI -myCryoWallThicknessI, myTCryostatBarrelH, 0.,myBottomCapH-myCryoCapWallThickness, myTopCapH-myCryoCapWallThickness, myTopOffset, myBotOffset,false,1);
  for ( unsigned int i=0; i<pointCol->size(); ++i ) {
    if(pointCol->at(i).first*mm<4000*mm){
      myLiqArgonZ[temp_k]    = pointCol->at(i).first  *mm ;
      myLiqArgonRout[temp_k] = pointCol->at(i).second *mm ;
      temp_k++;
    } /*else{
      myUllageZ[temp_j]= pointCol->at(i).first  *mm ;
      myUllageRout[temp_j] = pointCol->at(i).second *mm ;
      temp_j++;
    } */
  }

  //  myLiqArgonZ[temp_k-1]=myUllageLevel;
  //myLiqArgonRout[temp_k-1] = myUllageRout[0];
  // myUllageZ[0]=myUllageLevel;


  myNumPointsLiqArgon = temp_k;
  //myNumPointsUllage = temp_j;
  //for (int i = 0; i<myNumPointsLiqArgon; i++)
  //  G4cout<<" LAr : "<<myLiqArgonZ[i]<<" "<<myLiqArgonRout[i]<<G4endl;
  //for (int i = 0; i<myNumPointsUllage; i++)
  //  G4cout<<" Ullage : "<<myUllageZ[i]<<" "<<myUllageRout[i]<<G4endl;



   myNumPointsLiqArgon = pointCol->size();
  //PEN layer
  //pointCol = DSDetectorDS20kNeutronVeto::createGeometry(myDefaultRI - myCryoWallThicknessI-myLArInterfaceThickness, myTPCHeight, 0.,myBottomCapH-myVacuumH-myCryoWallThicknessI-myLArInterfaceThickness, myTopCapH-myVacuumH-myCryoWallThicknessI-myLArInterfaceThickness, 0., myTopOffset, myBotOffset);
   pointCol = DSDetectorDS20kNeutronVeto::createGeometry(myDefaultRI - myCryoWallThicknessI-myLArInterfaceThickness, myTCryostatBarrelH, 0.,myBottomCapH-myCryoCapWallThickness-myLArInterfaceThickness, myTopCapH-myCryoCapWallThickness-myLArInterfaceThickness, myTopOffset, myBotOffset,false,1);

  temp_j =0;
  temp_k =0;
  for ( unsigned int i=0; i<pointCol->size(); ++i ) {
    if(pointCol->at(i).first*mm<4000*mm){
      myPlasticZ[temp_k]    = pointCol->at(i).first  *mm ;
      myPlasticRout[temp_k] = pointCol->at(i).second *mm ;
      temp_k++;
      
    } /*else{
      myPlasticUllageZ[temp_j]    = pointCol->at(i).first  *mm ;
      myPlasticUllageRout[temp_j] = pointCol->at(i).second *mm ;
      temp_j++;
    } */
  }

  //  myPlasticZ[temp_k-1]=myUllageLevel;
  //  myPlasticRout[temp_k-1] = myPlasticUllageRout[0];
  // myPlasticUllageZ[0]=myUllageLevel;

  myNumPointsPlastic = pointCol->size();
  myNumPointsPlastic = temp_k;
  //myNumPointsPlasticUllage = temp_j;
  //for (int i = 0; i<myNumPointsPlastic; i++)
  //         G4cout<<" PEN in LAr : "<<myPlasticZ[i]<<" "<<myPlasticRout[i]<<G4endl;
  //for (int i = 0; i<myNumPointsUllage; i++)
  //         G4cout<<" PEN in GAr : "<<myPlasticUllageZ[i]<<" "<<myPlasticUllageRout[i]<<G4endl;


  //LAr inside PEN - below here, splitting between gas and liquid UAr

  //pointCol = DSDetectorDS20kNeutronVeto::createGeometry(myDefaultRI -myCryoWallThicknessI-myLArInterfaceThickness-myPlasticThickness, myTPCHeight, 0.,myBottomCapH-myVacuumH-myCryoWallThicknessI -myLArInterfaceThickness- myPlasticThickness, myTopCapH-myVacuumH-myCryoWallThicknessI -myLArInterfaceThickness- myPlasticThickness, 0., myTopOffset, myBotOffset);
  pointCol = DSDetectorDS20kNeutronVeto::createGeometry(myDefaultRI -myCryoWallThicknessI-myLArInterfaceThickness-myPlasticThickness, myTCryostatBarrelH, 0.,myBottomCapH-myCryoCapWallThickness -myLArInterfaceThickness- myPlasticThickness, myTopCapH-myCryoCapWallThickness -myLArInterfaceThickness- myPlasticThickness, myTopOffset, myBotOffset,false,1);

  temp_j =0;
  temp_k =0;
  for ( unsigned int i=0; i<pointCol->size(); ++i ) {
      myLiqArgonZI[temp_k]    = pointCol->at(i).first  *mm ;
      myLiqArgonRoutI[temp_k] = pointCol->at(i).second *mm ;
      temp_k++;
    if(pointCol->at(i).first*mm>myUllageLevel){
      myUllageZI[temp_j]= pointCol->at(i).first  *mm ;
      myUllageRoutI[temp_j] = pointCol->at(i).second *mm ;
      temp_j++;
    }
  }
  //myNumPointsLiqArgonI = pointCol->size();
  myNumPointsLiqArgonI = temp_k;
  myNumPointsUllageI = temp_j;

  //for (int i = 0; i<myNumPointsLiqArgonI; i++)
  // G4cout<<" LArI : "<<myLiqArgonZI[i]<<" "<<myLiqArgonRoutI[i]<<G4endl;
  // for (int i = 0; i<myNumPointsUllageI; i++)
  // G4cout<<" UllageI : "<<myUllageZI[i]<<" "<<myUllageRoutI[i]<<G4endl;

  myUllageZI[0]=myUllageLevel;

  //-------------------------//
  //      Outer Cryo         //
  //-------------------------//

  //remember to shift a bit the cryostat placement compared to the acrylic vessel if it is there
  //if ( DSStorage::Get()->GetAcrylicVesselDiameter()  )  myVerticalShift =  DSStorage::Get()->GetDS20kTPCverticalShift() ;
  //fPhysicDS20k  = new G4PVPlacement( myDefaultRotation2, G4ThreeVector(0.,0., - myVerticalShift), "OuterCryostat", fLogicDS20k,fMotherVolume, false, 0, myCheckOverlap );
  //G4VPhysicalVolume *fPhysicDS20k  = new G4PVPlacement( myDefaultRotation2, G4ThreeVector(0.,0., - myVerticalShift), "OuterCryostat", fLogicDS20k,fMotherVolume, false, 0, myCheckOverlap );
//  G4VPhysicalVolume *  fPhysicDS20k  = new G4PVPlacement( 0, G4ThreeVector(0.,0., 0), "OuterCryostat", fLogicDS20k,fMotherVolume, false, 0, myCheckOverlap );

  //fLogicDS20k->SetVisAttributes(myGreen);
  //fLogicDS20k->SetVisAttributes(G4VisAttributes::GetInvisible());

//  G4Polycone      *fSolidDS20kVac   = new G4Polycone( "OuterCryostat_SolidVac", 0, myTwoPi, myNumPointsVacCryo,myVacuumCryostatZ, myRmin, myVacuumCryostatRout  );
//  G4LogicalVolume *fLogicDS20kVac   = new G4LogicalVolume( fSolidDS20kVac, DSMaterial::Get()->GetVacuum(), "OuterCryostat_LogicVac" );
//  G4VPhysicalVolume *fPhysicDS20kVac  = new G4PVPlacement( 0, myZeros, "OuterCryostatVac", fLogicDS20kVac, fPhysicDS20k, false, 0, myCheckOverlap );
  //fLogicDS20kVac->SetVisAttributes(G4VisAttributes::GetInvisible());



  //-------------------------//
  //      Inner Vessel         //
  //-------------------------//

  G4Polycone      *fSolidInnerCryo   = new G4Polycone( "InnerCryo_Solid", 0, myTwoPi, myNumPointsInnerCryo,myInnerCryostatZ, myRmin,  myInnerCryostatRout );
  G4LogicalVolume *fLogicInnerCryo   = new G4LogicalVolume( fSolidInnerCryo,myCryoMat , "SolidInnerCryo_Logic" );
  fPhysicInnerCryo  = new G4PVPlacement( 0, G4ThreeVector(0,0,-20*cm), "InnerCryo", fLogicInnerCryo,fMotherVolume , false, 0, myCheckOverlap );
  //-------------------------//
  //LAr buffer inside PEN    //
  //-------------------------//

  // evrything must be divided in 2 pieces in Z for LAr and GAr regions
  G4Polycone      *fSolidLArBuffer   = new G4Polycone( "SoliLArBuffer_Solid", 0, myTwoPi, myNumPointsLiqArgon, myLiqArgonZ, myRmin, myLiqArgonRout  );
  G4LogicalVolume *fLogicLArBuffer   = new G4LogicalVolume( fSolidLArBuffer,  DSMaterial::Get()->GetNSLiquidArgon(), "SolidLArBuffer_Logic" );
  fPhysicLArBuffer  = new G4PVPlacement( 0, myZeros, "LArBuffer", fLogicLArBuffer, fPhysicInnerCryo, false, 0, myCheckOverlap );
/*
  G4Polycone      *fSolidGArBuffer   = new G4Polycone( "SoliGArBuffer_Solid", 0, myTwoPi, myNumPointsUllage, myUllageZ, myRmin, myUllageRout  );
  G4LogicalVolume *fLogicGArBuffer   = new G4LogicalVolume( fSolidGArBuffer,  DSMaterial::Get()->GetGaseousArgon(), "GArBuffer_Logic" );
  fPhysicGArBuffer  = new G4PVPlacement( 0, myZeros, "GArBuffer", fLogicGArBuffer, fPhysicInnerCryo, false, 0, myCheckOverlap );
*/
  //-------------------------//
  //      plastic layer      //
  //-------------------------//
  G4Polycone      *fSolidPlastic   = new G4Polycone( "SolidPlastic_Solid", 0, myTwoPi, myNumPointsPlastic, myPlasticZ, myRmin, myPlasticRout  );
  G4LogicalVolume *fLogicPlastic   = new G4LogicalVolume( fSolidPlastic,  DSMaterial::Get()->GetPEN(), "SolidPlastic_Logic" );
  fPhysicPlastic  = new G4PVPlacement( 0, myZeros, "SolidPlastic", fLogicPlastic, fPhysicLArBuffer, false, 0, myCheckOverlap );

/*
  G4Polycone      *fSolidPlasticUllage   = new G4Polycone( "SolidPlasticUllage_Solid", 0, myTwoPi, myNumPointsPlasticUllage, myPlasticUllageZ, myRmin, myPlasticUllageRout  );
  G4LogicalVolume *fLogicPlasticUllage   = new G4LogicalVolume( fSolidPlasticUllage,  DSMaterial::Get()->GetPEN(), "SolidPlasticUllage_Logic" );
  G4VPhysicalVolume  * fPhysicPlasticUllage  = new G4PVPlacement( 0, myZeros, "SolidPlasticUllage", fLogicPlasticUllage, fPhysicGArBuffer, false, 0, myCheckOverlap );
*/
  //fLogicPlastic->SetVisAttributes(G4VisAttributes::GetInvisible());


  // evrything must be divided in 2 pieces in Z for LAr and GAr regions
  G4RotationMatrix* rotUAr = new G4RotationMatrix();
  rotUAr->rotateZ(22.5*deg);

  // G4cout<<"UAr "<<G4endl;
  G4Polycone      *fSolidUAr   = new G4Polycone( "SolidUAr_Solid", 0, myTwoPi, myNumPointsLiqArgonI,myLiqArgonZI, myRmin, myLiqArgonRoutI );
  G4LogicalVolume *fLogicUAr   = new G4LogicalVolume( fSolidUAr, DSMaterial::Get()->GetIVLiquidArgon(), "SolidUAr_Logic" );
  fPhysicUAr  = new G4PVPlacement( rotUAr, myZeros, "SolidUAr", fLogicUAr, fPhysicPlastic, false, 0, myCheckOverlap );

  G4Polycone      *fSolidGUAr   = new G4Polycone( "SolidGUAr_Solid", 0, myTwoPi, myNumPointsUllageI,myUllageZI, myRmin, myUllageRoutI );
  G4LogicalVolume *fLogicGUAr   = new G4LogicalVolume( fSolidGUAr, DSMaterial::Get()->GetGaseousArgon(), "GUAr_Logic" );
  //fPhysicGUAr  = new G4PVPlacement(rotUAr, myZeros, "GUAr", fLogicGUAr, fPhysicPlasticUllage, false, 0, myCheckOverlap );
  fPhysicGUAr  = new G4PVPlacement(rotUAr, myZeros, "GUAr", fLogicGUAr, fPhysicUAr, false, 0, myCheckOverlap );
  // G4VisAttributes *myLarAttributes = new G4VisAttributes(G4Colour  (0.0, 0.0, 1.0, 0.9));

  //fLogicUAr->SetVisAttributes(G4VisAttributes::GetInvisible());

  //-------------------------//
  //      Plastic Shell      //
  //-------------------------//
  const G4int fPS_numZPlanes = 6 ;
  G4double fPS_zplanes[fPS_numZPlanes] = {0,0,0,0,0,0} ;
  G4double fPS_rInner[fPS_numZPlanes]  = {0,0,0,0,0,0} ;
  G4double fPS_rOuter[fPS_numZPlanes]  = {0,0,0,0,0,0} ;
  G4double fPSdistance  = 5*cm + 2*cm ; //rib thickness + distance
  G4double fPSthichness = DSStorage::Get()->GetDS20kHDPEShellThickness() ;
  G4double fPSbarrelHeight = 4740*mm ;
  G4double fPScapHeight = 840*mm ;
  G4double fPScapRadius = 120*cm / 2. ;
  G4double fPSinnerRadius = 5000*mm / 2. ;
  G4double fPScapAngle  = atan( (fPSinnerRadius-fPScapRadius) / fPScapHeight) ;
  //G4cout << fPScapAngle <<  " " << fPScapRadius - fPSthichness/sin(fPScapAngle) << G4endl ;

  fPS_zplanes[0] = - fPSbarrelHeight/2. - fPScapHeight;
  fPS_zplanes[1] = - fPSbarrelHeight/2. - 0.5*fPSthichness ;
  fPS_zplanes[2] = - fPSbarrelHeight/2. ;
  fPS_zplanes[3] = + fPSbarrelHeight/2. ;
  fPS_zplanes[4] = + fPSbarrelHeight/2. + 0.5*fPSthichness;
  fPS_zplanes[5] = + fPSbarrelHeight/2. + fPScapHeight;
  fPS_rInner[0] = fPScapRadius - fPSthichness/sin(fPScapAngle) - fPSthichness;
  fPS_rInner[1] = fPSinnerRadius  - fPSthichness/sin(fPScapAngle);
  fPS_rInner[2] = fPSinnerRadius + fPSdistance;
  fPS_rInner[3] = fPSinnerRadius + fPSdistance;
  fPS_rInner[4] = fPSinnerRadius - fPSthichness/sin(fPScapAngle);
  fPS_rInner[5] = fPScapRadius - fPSthichness/sin(fPScapAngle) - fPSthichness;
  fPS_rOuter[0] = fPScapRadius;
  fPS_rOuter[1] = fPSinnerRadius + fPSdistance + fPSthichness;
  fPS_rOuter[2] = fPSinnerRadius + fPSdistance + fPSthichness;
  fPS_rOuter[3] = fPSinnerRadius + fPSdistance+ fPSthichness;
  fPS_rOuter[4] = fPSinnerRadius + fPSdistance+ fPSthichness;
  fPS_rOuter[5] = fPScapRadius;


  if (fPSthichness > 0) {
    G4RotationMatrix* myDefaultRotation3 = new G4RotationMatrix();
    G4double myDefaultAngle3 = myTwoPi / 24.;
    myDefaultRotation3->rotateZ(-myDefaultAngle3);  // to correct for octagonal symmetry
    G4Polyhedra * fPlasticShellSolid      = new G4Polyhedra("PlasticShellSolid",0, myTwoPi, /*numSide*/12, /*numZPlanes*/ fPS_numZPlanes, fPS_zplanes, fPS_rInner, fPS_rOuter) ;
    G4LogicalVolume *fPlasticShellLogic   = new G4LogicalVolume( fPlasticShellSolid, DSMaterial::Get()->GetHDPE() , "PlasticShellLogic" );
    new G4PVPlacement( myDefaultRotation3  , G4ThreeVector(0,0,-20*cm), "PlasticShell", fPlasticShellLogic,fMotherVolume , false, 0, myCheckOverlap );

    G4double fPScap_positioning = fPScapHeight + fPSbarrelHeight/2. + fPSthichness/2. ;
    G4Tubs * fPScapSolid = new G4Tubs ("PScapSolid", 0, fPScapRadius, fPSthichness/2.,   0, myTwoPi);
    G4LogicalVolume * fPScapLogic   = new G4LogicalVolume( fPScapSolid, DSMaterial::Get()->GetHDPE() , "PScapLogic" );
    new G4PVPlacement( 0, G4ThreeVector(0,0,-20*cm-fPScap_positioning), "PSbottomCap", fPScapLogic,fMotherVolume , false, 0, myCheckOverlap );
    new G4PVPlacement( 0, G4ThreeVector(0,0,-20*cm+fPScap_positioning), "PStopCap", fPScapLogic,fMotherVolume , false, 0, myCheckOverlap );

    G4VisAttributes *myPSAttributes = new G4VisAttributes(G4Colour  (0.0, 0.0, 1.0, 1));
    fPScapLogic->SetVisAttributes(myPSAttributes) ;
    fPlasticShellLogic->SetVisAttributes(myPSAttributes) ;
  }


  //PA2021 remove surfaces
  //fPMTNeutronVeto20k->ConstructDefineSurfaces();
  DefineSurfaces();
}

DSDetectorDS20kNeutronVeto::~DSDetectorDS20kNeutronVeto() {
  ;  // delete fMessenger;
}

void DSDetectorDS20kNeutronVeto::DefineSurfaces() {

  ///////////////////////////////////////////////////////////////////////////////
  // PEN --> Reflector (Reflector)
  // Should be no Reflector --> TPB as surface is defined as dielectric_metal.
  ///////////////////////////////////////////////////////////////////////////////

  G4double TREFUV = DSParameters::Get()->GetTeflonTPBUVRef();
  G4double TREFVIS = DSParameters::Get()->GetTeflonTPBVisRef();

  // Choose 21 to use reflectivity of PEN air coupled to ESR or choose 22 to use
  // pure ESR reflectivity for PEN

  int PEN_reflectivity = 1;
   // PA original code uses a switch. Choose one default here
   switch (PEN_reflectivity) {

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

         G4OpticalSurface* fOpTPBReflectorSurface = new
         G4OpticalSurface("OpTPBReflectorSurface");

         new G4LogicalBorderSurface("TPBReflectorSurface", fPhysicPlastic,  fPhysicLArBuffer , fOpTPBReflectorSurface );

         fOpTPBReflectorSurface->SetType( dielectric_metal );

         fOpTPBReflectorSurface->SetModel(unified);

         fOpTPBReflectorSurface->SetFinish(polished);


         G4MaterialPropertiesTable *fOpTPBReflectorSurfProp = new
         G4MaterialPropertiesTable();

         fOpTPBReflectorSurfProp->AddProperty("REFLECTIVITY", TeflonTPBENE,
         TeflonTPBREF, 50); // Should be 50 for the arrays written above

         fOpTPBReflectorSurface->SetMaterialPropertiesTable(
         fOpTPBReflectorSurfProp );

       }
       break;
     case 2:
       {

  ///////////////////////////
  // Pure ESR reflectivity //
  ///////////////////////////

  G4double myvalue, myene;
  G4int dim;

  G4double TeflonTPBENE[831], TeflonTPBREF[831];
  dim = 0;
  // ifstream fESR_reflectivity = DSIO::GetStreamDS20kESRreflectivity();
  // while(!fESR_reflectivity.eof()) {
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

  for (int i = 0; i < dim; i++)
    TeflonTPBREF[i] *= TREFVIS / 100.;  // PA: rescaling by the value read in in
                                        // DSOpticsDS20k.dat? TODO: Confirm with Cenk

  G4OpticalSurface* fOpTPBReflectorSurface = new G4OpticalSurface("OpTPBReflectorSurface");
  new G4LogicalBorderSurface("TPBReflectorSurface", fPhysicLArBuffer, fPhysicInnerCryo, fOpTPBReflectorSurface);
  fOpTPBReflectorSurface->SetType(dielectric_metal);
  fOpTPBReflectorSurface->SetModel(unified);
  fOpTPBReflectorSurface->SetFinish(polished);

  G4MaterialPropertiesTable* fOpTPBReflectorSurfProp = new G4MaterialPropertiesTable();
  fOpTPBReflectorSurfProp->AddProperty("REFLECTIVITY", TeflonTPBENE, TeflonTPBREF,
                                       dim);  // Should be 50 for the arrays written above
  fOpTPBReflectorSurface->SetMaterialPropertiesTable(fOpTPBReflectorSurfProp);

  // switch terminates here
      }
      break;
  }
}

PointColPtr DSDetectorDS20kNeutronVeto::createGeometry(double r0, double hc, double z0, double hb, double ht, double mytopoffset, double mybotoffset,bool draw_ribs, int inner) {
  // G4cout<<r0<<" "<<hc<<" "<<z0<<" "<<hb<<" "<<ht<<" "<<offset<<" "<<mytopoff<<" "<<mybotoff<<G4endl;
  // creates a cylindric shape with end caps in the form of rotation ellipsoids
  //
  // input:
  // r0: radius of cylinder
  // hc: height of cylinder
  // z0: symmetry plane of cylinder
  // hb: height of bottom cap
  // ht: height of top cap
  // offset: offset to the otherwise defined geometry
  // mytopoffset: height of cylindrical shape between top cap and top rib
  // mybotoffset: height of cylindrical shape between bottom cap and bottom rib
  // draw_ribs:
  // inner: 
  // half height of cylinder
  const double hc2 = 0.5 * hc;
  
  // Flange and ribs thickness
  double flange_thick = 8.8*cm;
  double rib_thick = 3.5*cm;

  //const int nb = 35;
  //const int nc = 36;
 // const int nt = 35;
  const int nb = 90;
  //const int nc = 20;
  const int nt = 90;

  PointColPtr points = new PointCol();
// bottom cap

  for ( int i=1; i<nb; ++i ) {
    const double angle = double(i) * 90.0 / double(nb);
    const double rads = M_PI / 180.0 * angle;
    double z = z0 - hc2 - mybotoffset - hb * std::cos(rads);
    double r = (r0-0.4*cm*inner) * std::sin(rads);
    points->push_back(PairFF(z,r));
  }

  points->push_back(PairFF(z0 -hc2 - mybotoffset,r0-0.4*cm*inner)); //Last point of the cap's curve
  points->push_back(PairFF(z0 -hc2 - 4.4*cm,r0-0.4*cm*inner));          //Top flange-cap connection point - 44 = half flange thickness

  // Barrel
  int ncc = 2;

  for(int j=0; j<16;j++){

    if (j==0) {
      for ( int i=0; i<ncc; ++i ) {              // Bottom flange generation (2 points)
        const double z = z0 -4.4*cm- hc2  + (flange_thick)*(double(i)/double(ncc-1));
        if (inner==0){
          const double r = r0+6.7*cm;                // Flange positive radial offset
          points->push_back(PairFF(z,r));}
        else{
          const double r = r0-6.7*cm;
          points->push_back(PairFF(z,r));}       // Flange negative radial offset
      }  
    }
    
    for ( int i=0; i<ncc; ++i ) {                // Barrel wall points generation - 16 pieces (2 points each)
      const double z=z0 -hc2+4.4*cm + j*rib_thick + j*18.5*cm + (18.5*cm)*(i/(ncc-1));
      const double r=r0;
      points->push_back(PairFF(z,r));    
    }

    if(draw_ribs){                                // 15 Ribs Generation (2 points each)
      if(j<15){
        for ( int i=0; i<ncc; ++i ) {           
          const double z = z0 - hc2 + 4.4*cm +(j)*rib_thick + (j+1)*18.5*cm + (rib_thick)*(double(i)/double(ncc-1));
          const double r = r0+3.2*cm;               // 32 = Ribs radial offset
          points->push_back(PairFF(z,r));
        }
      }
    }

    if (j==15) {
      for ( int i=0; i<ncc; ++i ) {              // Top flange generation (2 points)
        const double z = z0 -4.4*cm- hc2 + 348.5*cm+8.8*cm+ (flange_thick)*(double(i)/double(ncc-1));
        if (inner==0){
          const double r = r0+6.7*cm;                // Flange positive radial offset
          points->push_back(PairFF(z,r));}
        else{
          const double r = r0-6.7*cm;
          points->push_back(PairFF(z,r));}       // Flange negative radial offset
      }  
    }
  }
  
  points->push_back(PairFF(z0 +4.4*cm +hc2,r0-0.4*cm*inner));   //Bottom flange-cap connection point - 44 = half flange thickness

 // top cap

  for ( int i=1; i<nt; ++i ) {
    const double angle = double(i) * 90.0 / double(nt);
    const double rads = M_PI / 180.0 * angle;
    double z = z0 + mytopoffset+ hc2  + ht * std::sin(rads);
    double r = (r0-0.4*cm*inner) * std::cos(rads);
    points->push_back(PairFF(z,r));
  }

  points->push_back(PairFF(z0+hc2+mytopoffset+ht,0));   // Last point of the top's curve

  return points;
}

/*
PointColPtr DSDetectorDS20kNeutronVeto::createGeometry(double r0, double hc, double z0, double hb, double ht, double , double mytopoffset, double mybotoffset, bool draw_ribs) {
//  G4cout<< "Pointcol param: " <<r0<<" "<<hc<<" "<<z0<<" "<<hb<<" "<<ht<<" "<<offset<<" "<<mytopoffset<<" "<<mybotoffset<<G4endl;
  // creates a cylindric shape with end caps in the form of rotation ellipsoids
  //
  // input:
  // r0: radius of cylinder
  // hc: height of cylinder
  // z0: symmetry plane of cylinder
  // hb: height of bottom cap
  // ht: height of top cap
  // offset: offset to the otherwise defined geometry
  // mytopoffset: height of cylindrical shape between top cap and top rib
  // mybotoffset: height of cylindrical shape between bottom cap and bottom rib

  // half height of cylinder
  const double hc2 = 0.5 * hc;
  double rib_thick = 10.*cm;
  double rib_small_thick = 3.*cm;

  //Correction to obtain correct caps height
  hb = hb - mybotoffset - rib_thick*0.5;
  ht = ht - mytopoffset - rib_thick*0.5;
  const int nb = 110;
  const int nc = 20;
  const int nt = 110;

  PointColPtr points = new PointCol();

  // bottom cap

  points->push_back(PairFF(z0-hc2-mybotoffset-hb,0));
  //cout << "offset in the cryo contruction : " << offset << endl;
  for ( int i=1; i<nb; ++i ) {
    const double angle = double(i) * 90.0 / double(nb);
    const double rads = M_PI / 180.0 * angle;
    double z = z0 - hc2 - mybotoffset - hb * std::cos(rads);
    double r =            r0 * std::sin(rads);

    if ( offset > 0 ) {
      // normal vector
      double n_z = -r0 * std::cos(rads);
      double n_r =  hb * std::sin(rads);
      const double l = std::sqrt(n_z*n_z+n_r*n_r);
      // normal vector with length "offset"
      n_z *= offset / l;
      n_r *= offset / l;
      z += n_z;
      r += n_r;
    }

  //cout << z << "  " << r << endl;
  points->push_back(PairFF(z,r));
  }

  // cylinder

double step = (hc - 5*rib_thick)/4;
double step_small = (step - rib_small_thick)*0.5;

if (!draw_ribs) {
  for ( int i=0; i<nc; ++i ) {
    // normal in radial direction, offset in radius
    const double z = z0 - hc2 - mybotoffset + (hc + mybotoffset + mytopoffset)*(double(i)/double(nc-1));
    const double r = r0;
    points->push_back(PairFF(z,r));
  }
}
else {

int ncc = 4;

    //cout << "Cylinder: " << endl;

  for ( int i=0; i<ncc; ++i ) {
    // normal in radial direction, offset in radius
    const double z = z0 - hc2 - mybotoffset + (mybotoffset)*(double(i)/double(ncc-1));
    const double r = r0;
    //cout << z << "  " << r << endl;
    points->push_back(PairFF(z,r));
    }

  for (int j=0;j<5;j++) {

     for ( int i=0; i<ncc; ++i ) {
     // normal in radial direction, offset in radius
     const double z = z0 - hc2 + j*rib_thick + j*step + (rib_thick)*(double(i)/double(ncc-1));
     const double r = r0 + 7.*cm; //7
     //cout << z << "  " << r << endl;
     points->push_back(PairFF(z,r));
     }

   if (j<4) {


     for ( int i=0; i<ncc; ++i ) {
     // normal in radial direction, offset in radius
     const double z = z0 - hc2 + (j+1)*rib_thick + j*step + (step_small)*(double(i)/double(ncc-1));
     const double r = r0;
     //cout << z << "  " << r << endl;

     points->push_back(PairFF(z,r));
     }

         for ( int i=0; i<ncc; ++i ) {
     // normal in radial direction, offset in radius
         const double z = z0 - hc2 + (j+1)*rib_thick + j*step + step_small + (rib_small_thick)*(double(i)/double(ncc-1));
         const double r = r0 + 5.*cm; //5
         // cout << z << "  " << r << endl;
         points->push_back(PairFF(z,r));
         }


         for ( int i=0; i<ncc; ++i ) {
     // normal in radial direction, offset in radius
         const double z = z0 - hc2 + (j+1)*rib_thick + j*step + step_small + rib_small_thick + (step_small)*(double(i)/double(ncc-1));
         const double r = r0;
         // cout << z << "  " << r << endl;
         points->push_back(PairFF(z,r));
         }
   }
  }

  for ( int i=0; i<ncc; ++i ) {
    // normal in radial direction, offset in radius
    const double z = z0 - hc2 + 5.*rib_thick + 4.*step + (mytopoffset)*(double(i)/double(ncc-1));
    const double r = r0;
    // cout << z << "  " << r << endl;
    points->push_back(PairFF(z,r));
    }
}

 // top cap

  //cout << "Top dome: " << endl;

  for ( int i=1; i<nt; ++i ) {
    const double angle = double(i) * 90.0 / double(nt);
    const double rads = M_PI / 180.0 * angle;
    double z = z0 + hc2 + mytopoffset + ht * std::sin(rads);
    double r =            r0 * std::cos(rads);
    if ( offset > 0 ) {
      // normal vector
      double n_z = r0 * std::sin(rads);
      double n_r = ht * std::cos(rads);
      const double l = std::sqrt(n_z*n_z+n_r*n_r);
      // normal vector with length "offset"
      n_z *= offset / l;
      n_r *= offset / l;
      z += n_z;
      r += n_r;
    }

    //cout << z << "  " << r << endl;
    points->push_back(PairFF(z,r));
  }
  // normal in +z direction, offset in direction +z
  points->push_back(PairFF(z0+hc2+mytopoffset+ht,0));


  return points;
}*/
