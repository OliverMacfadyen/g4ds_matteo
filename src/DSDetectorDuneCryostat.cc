#include "DSDetectorDuneCryostat.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "DSMaterial.hh"
#include "DSLogger.hh"
#include "DSStorage.hh"
#include <iostream>
#include <fstream>
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4UIcommand.hh"
#include "G4UnionSolid.hh"

using namespace std;

DSDetectorDuneCryostat::DSDetectorDuneCryostat(G4VPhysicalVolume *myMotherVolume) {

  fMotherVolume = myMotherVolume;

  G4Colour  myWhite   (1.0, 1.0, 1.0) ;  // white
  G4Colour  myGray    (0.5, 0.5, 0.5) ;  // gray
  G4Colour  myBlack   (0.0, 0.0, 0.0) ;  // black
  G4Colour  myRed     (1.0, 0.0, 0.0) ;  // red
  G4Colour  myGreen   (0.0, 1.0, 0.0) ;  // green
  G4Colour  myBlue    (0.0, 0.0, 1.0) ;  // blue
  G4Colour  myCyan    (0.0, 1.0, 1.0) ;  // cyan
  G4Colour  myMagenta (1.0, 0.0, 1.0) ;  // magenta
  G4Colour  myYellow  (1.0, 1.0, 0.0) ;  // yellow

  DSLog(routine) << " Constructing Dune Cryostat Geometry" << endlog ;

  myMotherVolume->GetLogicalVolume()  ->SetVisAttributes(G4VisAttributes::GetInvisible()) ;

  //Activate the placement of inter beam plastic/water shielding
  bool IsPlaceAdditionalShield = false ;

  //Plastic Shield is actually filled with Air
  G4double side   = 8548*mm +1.6*m + 1.6*m ;
  G4double height =7900*mm  +1.6*m + 1.6*m ;

  G4Material * mat = DSMaterial::Get()->GetAir() ;
  G4Box * fSolidPlasticShield  = new G4Box("PlasticShield_Solid",side/2., side/2., height/2.);
  G4LogicalVolume * fLogicPlasticShield  = new G4LogicalVolume(fSolidPlasticShield, mat, "PlasticShield_Logic");
  G4VPhysicalVolume *fPhysicPlasticShield = new G4PVPlacement(0,
         G4ThreeVector(0,0,20*cm),
         "PlasticShield",
         fLogicPlasticShield,
         fMotherVolume,
         false,
         0,
         DSStorage::Get()->GetCheckOverlap());

  fMotherVolume=fPhysicPlasticShield ;

  //Test net. Two bar dimensions.
  //
  //         22/13.5
  //      _________
  //     |         |    2.8 / 1.2
  //     ----  ----
  //        |  |
  //        |  |    1.7/0.8 thickness
  //        |  |
  //        |  |
  //        |  |
  //     ___|  |___
  //    |          |
  //     ----------         61./27.5 tall


  //increase by 30% the amount of material
  double increase = 1. ;

  G4Box * Vertical_I_cap   = new G4Box ("Vertical_I_cap", increase * 1.4 *cm , 11*cm , 10756/2.*mm ) ;
  G4Box * Vertical_I_body  = new G4Box ("Vertical_I_up", 61/2. *cm ,increase *1.4 *cm ,  10756/2.*mm ) ;
  G4UnionSolid* solid_half  = new G4UnionSolid("Vertical_half", Vertical_I_cap , Vertical_I_body, 0 , G4ThreeVector(-30.5*cm,0,0));
  G4UnionSolid* solid_big   = new G4UnionSolid("Vertical_big",Vertical_I_cap ,solid_half,  0 , G4ThreeVector(61*cm,0,0));


  G4Box * Vertical_I_cap_small   = new G4Box ("Vertical_I_cap_small", increase * 0.6*cm , 6.5*cm , 10756/2.*mm ) ;
  G4Box * Vertical_I_body_small  = new G4Box ("Vertical_I_up_small", 27.5/2. *cm ,increase *0.6 *cm ,  10756/2.*mm ) ;
  G4UnionSolid* solid_half_small  = new G4UnionSolid("Vertical_half_small", Vertical_I_cap_small , Vertical_I_body_small, 0 , G4ThreeVector(-27.5/2.*cm,0,0));
  G4UnionSolid* solid_small       = new G4UnionSolid("Vertical_small",Vertical_I_cap_small ,solid_half_small,  0 , G4ThreeVector(+27.5*cm,0,0));


  G4RotationMatrix * grot = new G4RotationMatrix ;
  grot->rotateX(90*deg) ;
  grot->rotateZ(90*deg) ;

  G4double small_half_dist   = 10.4/12.*m ;

  G4LogicalVolume * fLogicTest     = new G4LogicalVolume(solid_big, DSMaterial::Get()->GetStainlessSteel(), "test_steel_beam");
  G4LogicalVolume * fLogicTest_sma = new G4LogicalVolume(solid_small, DSMaterial::Get()->GetStainlessSteel(), "test_small_beam");

  //Volumes of materials among beams
  //G4Box * Body  = new G4Box ("Body", 61/2. *cm , small_half_dist-2.8*cm ,  10756/2.*mm ) ;
  //G4Box * Body  = new G4Box ("Body", 30/2. *cm , small_half_dist/2-1.4*cm ,  10756/2.*mm ) ;
  G4Box * BodyShield  = new G4Box ("BodyShield", 12.5 *cm , 14.5*cm ,  10756/2.*mm ) ;
  G4LogicalVolume * fLogicBodyShield     = new G4LogicalVolume(BodyShield, DSMaterial::Get()->GetHDPE(), "logic_body");


  for (int i=0;i<22;++i) {
  //-x

    if ( IsPlaceAdditionalShield )
      new G4PVPlacement(0,
            G4ThreeVector(-5.51*m,-4.5*m +small_half_dist/4 +small_half_dist/2*i  ,0),
            "phys_body",
            fLogicBodyShield,
            fMotherVolume,
            false,
            0,
            DSStorage::Get()->GetCheckOverlap());

    if (i<12) new G4PVPlacement(0,
            G4ThreeVector(-5.8*m,-4.5*m +small_half_dist*i  ,0),
            "Test_steel_beam",
            fLogicTest,
            fMotherVolume,
            false,
            0,
            DSStorage::Get()->GetCheckOverlap());


    if (i<11)   new G4PVPlacement(0,
            G4ThreeVector(-5.8*m + 15*cm,-4.5*m +small_half_dist*(i+0.5)  ,0),
            "Test_steel_beam_small",
            fLogicTest_sma,
            fMotherVolume,
            false,
            0,
            DSStorage::Get()->GetCheckOverlap());
  }


  for (int i=0;i<22;++i) {
  //+x

    if ( IsPlaceAdditionalShield )
      new G4PVPlacement(0,
                        G4ThreeVector(+5.49*m,-4.5*m +small_half_dist/4 +small_half_dist/2*i  ,0),
                        "Test_steel_beam",
                        fLogicBodyShield,
                        fMotherVolume,
                        false,
                        0,
                        DSStorage::Get()->GetCheckOverlap());

    if (i<12)   new G4PVPlacement(0,
            G4ThreeVector(+5.2*m,-4.5*m +small_half_dist*i  ,0),
            "Test_steel_beam",
            fLogicTest,
            fMotherVolume,
            false,
            0,
            DSStorage::Get()->GetCheckOverlap());

    if (i<11)   new G4PVPlacement(0,
            G4ThreeVector(+5.2*m +15*cm,-4.5*m +small_half_dist*(i+0.5)  ,0),
            "Test_steel_beam_small",
            fLogicTest_sma,
            fMotherVolume,
            false,
            0,
            DSStorage::Get()->GetCheckOverlap());
  }


  for (int i=0;i<22;++i) {
  //roof

    if ( IsPlaceAdditionalShield )
      new G4PVPlacement(grot,
                      G4ThreeVector(-4.7*m +small_half_dist/4 +small_half_dist/2*i ,0 ,5.09*m ),
                      "Test_steel_beam",
                      fLogicBodyShield,
                      fMotherVolume,
                      false,
                      0,
                      DSStorage::Get()->GetCheckOverlap());

    if (i<12)    new G4PVPlacement(grot,
            G4ThreeVector(-4.7*m +small_half_dist*i ,0 ,4.8*m ),
            "Test_steel_beam",
            fLogicTest,
            fMotherVolume,
            false,
            0,
            DSStorage::Get()->GetCheckOverlap());
    if (i<11)   new G4PVPlacement(grot,
            G4ThreeVector(-4.7*m+small_half_dist*(i+0.5),0,4.8*m+15*cm ),
            "Test_steel_beam_small",
            fLogicTest_sma,
            fMotherVolume,
            false,
            0,
            DSStorage::Get()->GetCheckOverlap());
  }


  for (int i=0;i<22;++i) {
  //floor

    if ( IsPlaceAdditionalShield )
      new G4PVPlacement(grot,
                        G4ThreeVector(-4.7*m + small_half_dist/4 +small_half_dist/2*i,0 ,-5.1*m),
                        "Test_steel_beam",
                        fLogicBodyShield,
                        fMotherVolume,
                        false,
                        0,
                        DSStorage::Get()->GetCheckOverlap());

    if (i<12)  new G4PVPlacement(grot,
              G4ThreeVector(-4.7*m +small_half_dist*i ,0 ,-5.4*m),
            "Test_steel_beam",
            fLogicTest,
            fMotherVolume,
            false,
            0,
            DSStorage::Get()->GetCheckOverlap());
    if (i<11)   new G4PVPlacement(grot,
            G4ThreeVector(-4.7*m+small_half_dist*(i+0.5),0, - 5.4*m + 16*cm),
            "Test_steel_beam_small",
            fLogicTest_sma,
            fMotherVolume,
            false,
            0,
            DSStorage::Get()->GetCheckOverlap());
  }



  G4RotationMatrix * grot2 = new G4RotationMatrix ;
  grot2->rotateZ(90*deg) ;
  grot2->rotateX(90*deg) ;
/*  */

  for (int i=0;i<22;++i) {
  //-y

    if ( IsPlaceAdditionalShield )
      new G4PVPlacement(grot2,
                        G4ThreeVector(0 , 5.51*m ,  -4.65*m +small_half_dist*0.97/4 +small_half_dist*0.97/2*i ),
                        "Test_steel_beam",
                        fLogicBodyShield,
                        fMotherVolume,
                        false,
                        0,
                        DSStorage::Get()->GetCheckOverlap());

    if (i<12)  new G4PVPlacement(grot2,
            G4ThreeVector(0 , 5.8*m ,  -4.65*m +small_half_dist*i*0.97 ),
            "Test_steel_beam",
            fLogicTest,
            fMotherVolume,
            false,
            0,
            DSStorage::Get()->GetCheckOverlap());

    if (i<11)   new G4PVPlacement(grot2,
            G4ThreeVector(0,5.8*m -15*cm, -4.65*m +small_half_dist*(i+0.5)*0.97 ),
            "Test_steel_beam_small",
            fLogicTest_sma,
            fMotherVolume,
            false,
            0,
            DSStorage::Get()->GetCheckOverlap());
  }

  for (int i=0;i<22;++i) {
  //+y

    if ( IsPlaceAdditionalShield )
      new G4PVPlacement(grot2,
                        G4ThreeVector(0 , -5.49*m ,  -4.65*m +small_half_dist*0.97/4 + small_half_dist*0.97/2*i ),
                        "Test_steel_beam",
                        fLogicBodyShield,
                        fMotherVolume,
                        false,
                        0,
                        DSStorage::Get()->GetCheckOverlap());

    if (i<12)   new G4PVPlacement(grot2,
            G4ThreeVector(0 , -5.2*m ,  -4.65*m +small_half_dist*i*0.97 ),
            "Test_steel_beam",
            fLogicTest,
            fMotherVolume,
            false,
            0,
            DSStorage::Get()->GetCheckOverlap());
    if (i<11)   new G4PVPlacement(grot2,
            G4ThreeVector(0,-5.2*m - 15*cm, -4.65*m +small_half_dist*0.97*(i+0.5) ),
            "Test_steel_beam_small",
            fLogicTest_sma,
            fMotherVolume,
            false,
            0,
            DSStorage::Get()->GetCheckOverlap());
  }


  /*
  bool isSlowDetailedGeometry = 0;
  //this works, but it is ultra slow!
  if (isSlowDetailedGeometry) {
    G4Box * Vertical_I_cap   = new G4Box ("Vertical_I_cap", 1.4 *cm , 11*cm , 10756/2.*mm ) ;
    G4Box * Vertical_I_body  = new G4Box ("Vertical_I_up", 61/2. *cm ,1.4 *cm ,  10756/2.*mm ) ;
    G4UnionSolid* solid_half  = new G4UnionSolid("Vertical_half", Vertical_I_cap , Vertical_I_body, 0 , G4ThreeVector(-30.5*cm,0,0));
    G4UnionSolid* solid_big   = new G4UnionSolid("Vertical_big",Vertical_I_cap ,solid_half,  0 , G4ThreeVector(61*cm,0,0));

    G4Box * Vertical_i_cap   = new G4Box ("Vertical_i_cap", 0.6 *cm , 13.5/2*cm , 10756/2.*mm ) ;
    G4Box * Vertical_i_body  = new G4Box ("Vertical_i_up", 27.5/2. *cm ,0.6*cm , 10756/2.*mm ) ;
    G4UnionSolid* small_half  = new G4UnionSolid("small_half", Vertical_i_cap , Vertical_i_body, 0 , G4ThreeVector(-275/2.*mm,0,0));
    G4UnionSolid* solid_small   = new G4UnionSolid("Vertical_bar",Vertical_i_cap ,small_half,  0 , G4ThreeVector(2*275/2.*mm,0,0));



    G4RotationMatrix * grot = new G4RotationMatrix ;
    grot->rotateX(90*deg) ;

    //make the first net.
    G4double small_half_dist   = 1.4/2.*m ;
    G4double small_half_dist_Z = 1.24/2.*m ;
    G4double big_dist = 1.3*m ;
    G4UnionSolid* intermediate_solid_vertical = new G4UnionSolid("int_vertical_small", solid_small , solid_big ,0,  G4ThreeVector(0,small_half_dist,0));
    //G4UnionSolid* intermediate_solid_vertical = new G4UnionSolid("int_vertical_small", solid_small , solid_small ,grot,  G4ThreeVector(0,+4.5*m,-4.5*m));

    G4double total_dist = 0 ;
    for (int i=0;i<13;++i) {
     cout << i <<" " << i%2 << " " <<  - 4.5 *m + (i+1)*small_half_dist_Z <<" " << (i+1)*small_half_dist << endl ;
     //if (i<11) continue ;

     if (i%2==1)  {
       intermediate_solid_vertical   = new G4UnionSolid("int_vertical_small", intermediate_solid_vertical , solid_small, 0 , G4ThreeVector(0, (i+1)*small_half_dist , 0));
      // intermediate_solid_vertical   = new G4UnionSolid("int_vertical_small", intermediate_solid_vertical , solid_small, grot , G4ThreeVector(0,  4.5 *m  ,  - 4.5 *m + (i+1)*small_half_dist_Z));

     } else  {
       intermediate_solid_vertical   = new G4UnionSolid("int_vertical_small", intermediate_solid_vertical , solid_big, 0 , G4ThreeVector(0, (i+1)*small_half_dist, 0));
     //  intermediate_solid_vertical   = new G4UnionSolid("int_vertical_small", intermediate_solid_vertical , solid_big, grot , G4ThreeVector(0, 4.5*m ,  - 4.5 *m + (i+1)*small_half_dist_Z));
       }
     }

    //G4UnionSolid* full_solid_vertical = new G4UnionSolid("int_vertical_small",  intermediate_solid_vertical, intermediate_solid_vertical ,grot,  G4ThreeVector(0, 4.5*m,-3.5*m));


    G4LogicalVolume * fLogicTest = new G4LogicalVolume(intermediate_solid_vertical, DSMaterial::Get()->GetStainlessSteel(), "Logical_steel_beam");

    fLogicTest ->SetVisAttributes(myRed) ;

    new G4PVPlacement(0,
           G4ThreeVector(  ( 8548 *mm + 2 * 80*cm  ) /2  + 2.4*cm , - 6.5 *small_half_dist ,0),
           "Phys_Steel_beam",
           fLogicTest,
           fMotherVolume,
           false,
           0,
           DSStorage::Get()->GetCheckOverlap());


  } else {
    ;
  }

    // so, use the following for the moment
    side   = 8548 *mm + 2 * 80*cm + 6.5*cm ;
    height = 7900 *mm + 2 * 80*cm + 6.5*cm;
    G4Box * fSolidSteelTank  = new G4Box("SteelTank_Solid",side/2., side/2., height/2.) ;
    G4LogicalVolume * fLogicSteelTank  = new G4LogicalVolume(fSolidSteelTank, DSMaterial::Get()->GetStainlessSteel(), "SteelTank_Logic");
    G4VPhysicalVolume * fPhysicSteelTank = new G4PVPlacement(0,
           G4ThreeVector(0,0,0),
           "SteelTank",
           fLogicSteelTank,
           fMotherVolume,
           false,
           0,
           DSStorage::Get()->GetCheckOverlap());
*/
  //}

  //start the sequence: from the outer wall:
  // 1 cm steel
  // 1 cm PlyWood
  // 50 cm insulating foam
  // 2 mm steel
  // 30 cm insulating foam
  // 1 cm PlyWood
  // 1.2 mm steel (membrane)


  side   = 8548 *mm + 2*80.1*cm + 2*11.2*mm;
  height = 7900 *mm + 2*80.1*cm + 2* 11.2*mm;
  fSolidSteelTank  = new G4Box("SteelTank_Solid",side/2., side/2., height/2.);
  fLogicSteelTank  = new G4LogicalVolume(fSolidSteelTank, DSMaterial::Get()->GetStainlessSteel(), "SteelTank_Logic");
  fPhysicSteelTank = new G4PVPlacement(0,
         G4ThreeVector(0,0,0),
         "SteelTank",
         fLogicSteelTank,
         fMotherVolume,
         false,
         0,
         DSStorage::Get()->GetCheckOverlap());

  G4VisAttributes * solidsteel= new G4VisAttributes(true);
  solidsteel->SetColour(1.0,1.0,1.0);  //white
  fLogicSteelTank->SetVisAttributes(solidsteel);

  side   = 8548 *mm + 2 + 2*79.1*cm + 2*11.2*mm ;
  height = 7900 *mm + 2 + 2*79.1*cm + 2*11.2*mm ;
  fSolidOuterPlyWood  = new G4Box("OuterPlyWood_Solid",side/2., side/2., height/2.);
  fLogicOuterPlyWood  = new G4LogicalVolume(fSolidOuterPlyWood, DSMaterial::Get()->GetPlyWood(), "OuterPlyWood_Logic");
  fPhysicOuterPlyWood = new G4PVPlacement(0,
         G4ThreeVector(0,0,0),
         "OuterPlyWood",
         fLogicOuterPlyWood,
         fPhysicSteelTank,
         false,
         0,
         DSStorage::Get()->GetCheckOverlap());

  side   = 8548 *mm + 2*78.1*cm + 2*11.2*mm   ;
  height = 7900 *mm + 2*78.1*cm + 2*11.2*mm ;
  fSolidOuterInsulatingFoam  = new G4Box("OuterInsulatingFoam_Solid",side/2., side/2., height/2.);
  fLogicOuterInsulatingFoam  = new G4LogicalVolume(fSolidOuterInsulatingFoam, DSMaterial::Get()->GetInsulatingFoam(), "OuterInsulatingFoam_Logic");
  fPhysicOuterInsulatingFoam = new G4PVPlacement(0,
         G4ThreeVector(0,0,0),
         "OuterInsulatingFoam",
         fLogicOuterInsulatingFoam,
         fPhysicOuterPlyWood,
         false,
         0,
         DSStorage::Get()->GetCheckOverlap());

  G4VisAttributes * InsulatingOuter= new G4VisAttributes(true);
  InsulatingOuter->SetColour(1.0,1.0,1.0);  //white
  fLogicOuterInsulatingFoam->SetVisAttributes(InsulatingOuter);

  //side   = 8548 *mm + 2*78.1*cm + 2*11.2*mm   ;
  side =  8548 *mm + 2*39.1*cm + 2*11.2*mm ;
  height = ( 781 - 391 )*mm  ;
  fSolidOuterInsulatingFoamRoof  = new G4Box("OuterInsulatingFoam_SolidRoof",side/2., side/2., height/2.);
  fLogicOuterInsulatingFoamRoof  = new G4LogicalVolume(fSolidOuterInsulatingFoamRoof, DSMaterial::Get()->GetInsulatingFoamRoof(), "OuterInsulatingFoam_LogicRoof");
  fPhysicOuterInsulatingFoamRoof = new G4PVPlacement(0,
         G4ThreeVector(0,0,((7900 *mm + 2*781*mm + 2*11.2*mm) - height)/2*mm),
         "OuterInsulatingFoamRoof",
         fLogicOuterInsulatingFoamRoof,
         fPhysicOuterInsulatingFoam,
         false,
         0,
         DSStorage::Get()->GetCheckOverlap());

 G4VisAttributes * InsulatingOuterRoof= new G4VisAttributes(true);
 InsulatingOuterRoof->SetColour(1.0,0.0,0.0);  //red
 fLogicOuterInsulatingFoamRoof->SetVisAttributes(InsulatingOuterRoof);

  side   = 8548 *mm + 2*39.1*cm + 2*11.2*mm   ;
  height = 7900 *mm + 2*39.1*cm + 2*11.2*mm ;
  fSolidInnerMembrane  = new G4Box("InnerMembrane_Solid",side/2., side/2., height/2.) ;
  fLogicInnerMembrane  = new G4LogicalVolume(fSolidInnerMembrane, DSMaterial::Get()->GetStainlessSteel(), "InnerMembrane_Logic");
  fPhysicInnerMembrane = new G4PVPlacement(0,
         G4ThreeVector(0,0,0),
         "InnerMembrane",
         fLogicInnerMembrane,
         fPhysicOuterInsulatingFoam,
         false,
         0,
         DSStorage::Get()->GetCheckOverlap());

 G4VisAttributes * InnerMembrane= new G4VisAttributes(true);
 InnerMembrane->SetColour(0.5,0.5,0.5);  //gray
 fLogicInnerMembrane->SetVisAttributes(InnerMembrane);

  side   = 8548 *mm + 2*39*cm + 2*11.2*mm ;
  height = 7900 *mm + 2*39*cm + 2*11.2*mm ;
  fSolidInsulatingFoam  = new G4Box("InsulatingFoam_Solid",side/2., side/2., height/2.);
  fLogicInsulatingFoam  = new G4LogicalVolume(fSolidInsulatingFoam, DSMaterial::Get()->GetInsulatingFoam(), "InsulatingFoam_Logic");
  fPhysicInsulatingFoam = new G4PVPlacement(0,
         G4ThreeVector(0,0,0),
         "InsulatingFoam",
         fLogicInsulatingFoam,
         fPhysicInnerMembrane,
         false,
         0,
         DSStorage::Get()->GetCheckOverlap());

  G4VisAttributes * Insulating= new G4VisAttributes(true);
  Insulating->SetColour(0.0,0.0,1.0);  //blue
  fLogicInsulatingFoam->SetVisAttributes(Insulating);

  side   = 8548 *mm + 2*39*cm + 2*11.2*mm ;
  height = 390*mm ;
  fSolidInsulatingFoamRoof  = new G4Box("InsulatingFoam_SolidRoof",side/2., side/2., height/2.);
  fLogicInsulatingFoamRoof  = new G4LogicalVolume(fSolidInsulatingFoamRoof, DSMaterial::Get()->GetInsulatingFoamRoof(), "InsulatingFoam_LogicRoof");
  fPhysicInsulatingFoamRoof = new G4PVPlacement(0,
         G4ThreeVector(0,0,((7900 *mm + 2*390*mm + 2*11.2*mm)-height)/2*mm),
         "InsulatingFoamRoof",
         fLogicInsulatingFoamRoof,
         fPhysicInsulatingFoam,
         false,
         0,
         DSStorage::Get()->GetCheckOverlap());

  G4VisAttributes * InsulatingRoof= new G4VisAttributes(true);
  InsulatingRoof->SetColour(0.0,1.0,0.0);  //green
  fLogicInsulatingFoamRoof->SetVisAttributes(InsulatingRoof);

  side   = 8548 *mm + 2 * 11.2*mm  ;
  height = 7900 *mm + 2 * 11.2*mm ;
  fSolidPlyWood  = new G4Box("PlyWood_Solid",side/2., side/2., height/2.);
  fLogicPlyWood  = new G4LogicalVolume(fSolidPlyWood, DSMaterial::Get()->GetPlyWood(), "PlyWood_Logic");
  fPhysicPlyWood = new G4PVPlacement(0,
         G4ThreeVector(0,0,0),
         "PlyWood",
         fLogicPlyWood,
         fPhysicInsulatingFoam,
         false,
         0,
         DSStorage::Get()->GetCheckOverlap());

  G4VisAttributes * PlyWood= new G4VisAttributes(true);
  PlyWood->SetColour(1.0,1.0,0.0);  //yellow
  fLogicPlyWood->SetVisAttributes(PlyWood);

  side   = 8548 *mm + 2 * 1.2*mm  ;
  height = 7900 *mm + 2 * 1.2*mm ;
  fSolidMembrane  = new G4Box("Membrane_Solid",side/2., side/2., height/2.);
  fLogicMembrane  = new G4LogicalVolume(fSolidMembrane, DSMaterial::Get()->GetStainlessSteel(), "Membrane_Logic");
  fPhysicMembrane = new G4PVPlacement(0,
         G4ThreeVector(0,0,0),
         "Membrane",
         fLogicMembrane,
         fPhysicPlyWood,
         false,
         0,
         DSStorage::Get()->GetCheckOverlap());



  G4double myPlasticThickness = 0 *cm ;
  if (myPlasticThickness/cm>0.0001) {  //add additional HDPE on the inner surface of cryo
    side   = 8548 *mm  ;
    height = 7900 *mm  ;
    fSolidAdditionalPlastic  = new G4Box("AdditionalPlastic_Solid", side/2., side/2., height/2.);
    fLogicAdditionalPlastic  = new G4LogicalVolume(fSolidAdditionalPlastic, DSMaterial::Get()->GetHDPE(), "AdditionalPlastic_Logic");
    fPhysicAdditionalPlastic = new G4PVPlacement(0,
           G4ThreeVector(0,0,0),
           "AdditionalPlastic",
           fLogicAdditionalPlastic,
           fPhysicMembrane,
           false,
           0,
           DSStorage::Get()->GetCheckOverlap());
  }



  side   = 8548 *mm  -2*myPlasticThickness;
  height = 7900 *mm  -2*myPlasticThickness;
  fSolidVetoArgon  = new G4Box("VetoArgon_Solid", side/2., side/2., height/2.);
  fLogicVetoArgon  = new G4LogicalVolume(fSolidVetoArgon, DSMaterial::Get()->GetOVLiquidArgon(), "VetoArgon_Logic");
  fPhysicVetoArgon = new G4PVPlacement(0,
         G4ThreeVector(0,0,0),
         "VetoArgon",
         fLogicVetoArgon,
         myPlasticThickness/mm>0.1?fPhysicAdditionalPlastic:fPhysicMembrane,
         false,
         0,
         DSStorage::Get()->GetCheckOverlap());



  //argon vapor sitting on top of the liquid
  side   = 8548 *mm -2*myPlasticThickness ;
  height = 400 *mm - myPlasticThickness ;
  fSolidArgonVapor  = new G4Box("ArgonVapor_Solid", side/2., side/2., height/2.);
  fLogicArgonVapor  = new G4LogicalVolume(fSolidArgonVapor, DSMaterial::Get()->GetGaseousArgon(), "ArgonVapor_Logic");
  fPhysicArgonVapor = new G4PVPlacement(0,
         G4ThreeVector(0,0,(7900-height)/2*mm-myPlasticThickness),
         "ArgonVapor",
         fLogicArgonVapor,
         fPhysicVetoArgon,
         false,
         0,
         DSStorage::Get()->GetCheckOverlap());



  fLogicArgonVapor      ->SetVisAttributes(G4VisAttributes::GetInvisible()) ;
  fLogicVetoArgon      ->SetVisAttributes(G4VisAttributes::GetInvisible()) ;
  fLogicMembrane      ->SetVisAttributes(G4VisAttributes::GetInvisible()) ;
  fLogicPlyWood        ->SetVisAttributes(G4VisAttributes::GetInvisible()) ;
  fLogicInsulatingFoam    ->SetVisAttributes(G4VisAttributes::GetInvisible()) ;
  fLogicInsulatingFoamRoof  ->SetVisAttributes(G4VisAttributes::GetInvisible()) ;
  fLogicInnerMembrane    ->SetVisAttributes(G4VisAttributes::GetInvisible()) ;
  fLogicOuterInsulatingFoam  ->SetVisAttributes(G4VisAttributes::GetInvisible()) ;
  fLogicOuterInsulatingFoamRoof  ->SetVisAttributes(G4VisAttributes::GetInvisible()) ;
  fLogicOuterPlyWood    ->SetVisAttributes(G4VisAttributes::GetInvisible()) ;
  //fLogicSteelTank      ->SetVisAttributes(G4VisAttributes::GetInvisible()) ;
  fLogicTest        ->SetVisAttributes(G4VisAttributes::GetInvisible()) ;
  fLogicTest_sma      ->SetVisAttributes(G4VisAttributes::GetInvisible()) ;

  fLogicPlasticShield->SetVisAttributes(G4VisAttributes::GetInvisible()) ;

  G4VisAttributes* beams_big = new G4VisAttributes(true);
  G4VisAttributes* beams_small = new G4VisAttributes(true);
  G4VisAttributes* body = new G4VisAttributes(true);

  beams_big->SetColour(1.0,0.0,0.0);
  beams_small->SetColour(1.0,1.0,0.0);
  body->SetColour(0.0,0.0,1.0);

  fLogicTest->SetVisAttributes(beams_big);
  fLogicTest_sma->SetVisAttributes(beams_small);
  fLogicBodyShield->SetVisAttributes(body);


}


DSDetectorDuneCryostat::~DSDetectorDuneCryostat(){
  ; //delete fMessenger;
}
