/* -------------------------------------- */
/* ArDM geometry based on DSDart detector*/
/* Olivier Dadoun */
/* March 2017 */
/* -------------------------------------- */
#include <iostream>
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4PVPlacement.hh"
#include "G4Sphere.hh"
#include "G4SubtractionSolid.hh"
#include "G4ThreeVector.hh"
#include "G4Torus.hh"
#include "G4Trd.hh"
#include "G4UIcommand.hh"
#include "G4UnionSolid.hh"

#include "G4Polyhedra.hh"

#include "G4Colour.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"

#include "DSDetectorArDM.hh"
#include "DSLogger.hh"
#include "DSMaterial.hh"
#include "DSStorage.hh"

#include "DSParameters.hh"
//#include "DSArDMpreparation.hh"

using namespace std;

DSDetectorArDM::DSDetectorArDM(G4VPhysicalVolume* myMotherVolume) {

  G4Colour myWhite(1.0, 1.0, 1.0);    // white
  G4Colour myGray(0.5, 0.5, 0.5);     // gray
  G4Colour myBlack(0.0, 0.0, 0.0);    // black
  G4Colour myRed(1.0, 0.0, 0.0);      // red
  G4Colour myGreen(0.0, 1.0, 0.0);    // green
  G4Colour myBlue(0.0, 0.0, 1.0);     // blue
  G4Colour myCyan(0.0, 1.0, 1.0);     // cyan
  G4Colour myMagenta(1.0, 0.0, 1.0);  // magenta
  G4Colour myYellow(1.0, 1.0, 0.0);   // yellow

  G4cout << " ************************************  " << G4endl;
  G4cout << "  Hello ArDM Detector Construction !!!" << G4endl;
  G4cout << " ************************************  " << G4endl;

  // Beginning import ArDM G4
  //  Init var
  ArDMvar.Init();

  // add PE shield
  addShield(myMotherVolume);
  addShieldLead(myMotherVolume);

  G4ThreeVector tankPos(0., 0., ArDMvar.TANK_CYLINDER_POS_Z);
  fTankMat = DSMaterial::Get()->GetStainlessSteel();
  fTankPhys = addTank(myMotherVolume, fTankMat, tankPos, 1);

  fTeflon = DSMaterial::Get()->GetTeflon();
  G4double z = ArDMvar.LARCOL_CYLINDER_POS_Z;
  G4ThreeVector pos(0., 0., z);
  G4Material* fLAr = DSMaterial::Get()->GetPseudoArgon();
  fPhysicLArBath = addLArColumn(myMotherVolume, fLAr, pos, 1);

  z = ArDMvar.GAR_COLUMN_Z;
  pos = G4ThreeVector(0., 0., z);
  G4Material* fGAr = DSMaterial::Get()->GetGaseousArgon();
  fGArColPhys = addGArColumn(myMotherVolume, fGAr, pos, 1);

  addWLS();

  addWLSSupport();
  addPMT("bottom");
  addPMT("top");

  addPMTSupport("bottom");
  addPMTSupport("top");

  addFieldShaperPillars();
  addFieldShaperRings();

  //#if ArDMvar.TURN_ON_SIDE_REFLECTOR
  addBtmSideReflector();
  addTopSideReflector();

  //#endif //ArDMvar.TURN_ON_SIDE_REFLECTOR
  // End import ArDM G4

  z = -ArDMvar.LARCOL_CYLINDER_POS_Z + DSStorage::Get()->GetArDMTestArgonHeight() / 2.;
  pos = G4ThreeVector(0., 0., z);
  G4Material* fArTPC = DSMaterial::Get()->GetLiquidArgon();
  fTestArgonPhys = addTestArgon(fPhysicLArBath, fArTPC, pos, 1);

  G4Region* fLArRegion = new G4Region("LAr_Logic");
  fLogicLArScint = fTestArgonPhys->GetLogicalVolume();
  fLogicLArScint->SetRegion(fLArRegion);
  fLArRegion->AddRootLogicalVolume(fLogicLArScint);

  // add DART SiPM plane
  //  addDartSiPM(fTestArgonPhys);
  DefineSurfaces();
}

void DSDetectorArDM::addShield(G4VPhysicalVolume* MotherVolume) {
  fMotherVolume = MotherVolume;
  if (!fMotherVolume) {
    cout << "in in DSDetectorArDM::addShield(), fMotherVolume = 0. exit." << endl;
    return;
  }

  G4ThreeVector fShieldPos(0., 0., 0.);
  fShieldMat = DSMaterial::Get()->GetHDPE();

  G4double innerR = ArDMvar.PE_SHIELD_INNER_RADIUS;
  G4double outerR = ArDMvar.PE_SHIELD_OUTER_RADIUS;
  G4double halfz = ArDMvar.PE_SHIELD_HALF_HEIGHT;

  // G4double topinnerR = 0.;
  // G4double topouterR = ArDMvar.PE_SHIELD_INNER_RADIUS;
  // G4double tophalfz = ArDMvar.PE_SHIELD_WIDTH;

  DSStorage::Get()->SetArDMShieldHeight(2. * halfz);
  DSStorage::Get()->SetArDMShieldRadius(outerR);

  // Vicente Pesudo (CIEMAT)
  //  G4double holeinnerR = 0.;
  //  G4double holeouterR = ArDMvar.PE_SHIELD_INNER_RADIUS;
  G4double holehalfz = ArDMvar.PE_SHIELD_HALF_HEIGHT - ArDMvar.PE_SHIELD_WIDTH;

  G4RotationMatrix* myDefaultRotation2 = new G4RotationMatrix();
  const double myTwoPi = 2 * M_PI * rad;
  G4double myDefaultAngle = myTwoPi / 16.;
  myDefaultRotation2->rotateZ(-myDefaultAngle);  // to correct for octagonal symmetry

  G4double zlimitSh8[2] = {-halfz, halfz};  // position of limiting surfaces of the prisma
  G4double rminSh8[2] = {0, 0};
  G4double rmaxSh8[2] = {outerR, outerR};  // outer radios at zlimit[i]

  G4Polyhedra* shield8FullSolid = new G4Polyhedra("Shield8Full", 0 * deg, 360. * deg, 8, 2, zlimitSh8, rminSh8, rmaxSh8);

  G4double zlimithole8[2] = {-holehalfz, holehalfz};  // position of limiting surfaces of the prisma
  G4double rminhole8[2] = {0, 0};
  G4double rmaxhole8[2] = {innerR, innerR};  // inner radios at zlimit[i]

  G4Polyhedra* shield8holeSolid = new G4Polyhedra("Shield8hole", 0 * deg, 360. * deg, 8, 2, zlimithole8, rminhole8, rmaxhole8);

  G4SubtractionSolid* shield8Solid = new G4SubtractionSolid("Shield8", shield8FullSolid, shield8holeSolid, 0, fShieldPos);

  G4LogicalVolume* fShield8Log = new G4LogicalVolume(shield8Solid, fShieldMat, "Shield8");
  // G4VPhysicalVolume* fShield8 =
  //     new G4PVPlacement(myDefaultRotation2, fShieldPos, "Shield8",
  //     fShield8Log, fMotherVolume, false, 0);

  G4VisAttributes* solidAtt = new G4VisAttributes(true);
  // G4VisAttributes* solidAtt = new G4VisAttributes(false);
  solidAtt->SetColour(1.0, 1.0, 0.0);
  solidAtt->SetForceAuxEdgeVisible(true);
  fShield8Log->SetVisAttributes(solidAtt);
}

// Edgar Sanchez (CIEMAT) DART Lead Test
void DSDetectorArDM::addShieldLead(G4VPhysicalVolume* MotherVolume) {
  fMotherVolume = MotherVolume;

  // By default, skip adding the ArDM lead rinh
  if (DSStorage::Get()->GetArDMLeadHeight() <= 0.) return;

  if (!fMotherVolume) {
    cout << "in in DSDetectorArDM::addShieldLead(), fMotherVolume = 0. exit." << endl;
    return;
  }

  G4ThreeVector fShieldLeadPos(0., 0., 0.);
  fShieldLeadMat = DSMaterial::Get()->GetMetalLead();

  G4double innerR = ArDMvar.PE_SHIELD_LEAD_INNER_RADIUS;
  G4double outerR = ArDMvar.PE_SHIELD_LEAD_OUTER_RADIUS;

  G4double halfz = (DSStorage::Get()->GetArDMLeadHeight()) / 2.;

  // G4double topinnerR = 0.;
  // G4double topouterR = ArDMvar.PE_SHIELD_LEAD_INNER_RADIUS;
  // G4double tophalfz = ArDMvar.PE_SHIELD_LEAD_WIDTH;

  DSStorage::Get()->SetArDMShieldLeadRadius(outerR);

  // G4double holehalfz = (DSStorage::Get()->GetArDMLeadHeight()) / 2. -
  // ArDMvar.PE_SHIELD_LEAD_WIDTH;

  G4double zlimitSh8[2] = {-halfz, halfz};  // position of limiting surfaces of the prisma
  G4double rminSh8[2] = {0, 0};
  G4double rmaxSh8[2] = {outerR, outerR};  // outer radios at zlimit[i]

  G4Polyhedra* shieldLead8FullSolid = new G4Polyhedra("ShieldLead8Full", 0 * deg, 360. * deg, 8, 2, zlimitSh8, rminSh8, rmaxSh8);

  G4double zlimithole8[2] = {-halfz, halfz};  // position of limiting surfaces of the prisma
  G4double rminhole8[2] = {0, 0};
  G4double rmaxhole8[2] = {innerR, innerR};  // inner radios at zlimit[i]

  G4Polyhedra* shieldLead8holeSolid = new G4Polyhedra("ShieldLead8hole", 0 * deg, 360. * deg, 8, 2, zlimithole8, rminhole8, rmaxhole8);

  G4SubtractionSolid* shieldLead8Solid = new G4SubtractionSolid("ShieldLead8", shieldLead8FullSolid, shieldLead8holeSolid, 0, fShieldLeadPos);

  G4LogicalVolume* fShieldLead8Log = new G4LogicalVolume(shieldLead8Solid, fShieldLeadMat, "ShieldLead8");
  // G4VPhysicalVolume* fShieldLead8 =
  //     new G4PVPlacement(0, fShieldLeadPos, "ShieldLead8", fShieldLead8Log,
  //     fMotherVolume, false, 0);

  G4VisAttributes* solidAtt = new G4VisAttributes(true);
  // G4VisAttributes* solidAtt = new G4VisAttributes(false);  ///No
  // visualization
  solidAtt->SetColour(1.0, 0.0, 0.0);
  solidAtt->SetForceAuxEdgeVisible(true);
  fShieldLead8Log->SetVisAttributes(solidAtt);
}

G4VPhysicalVolume* DSDetectorArDM::addTank(G4VPhysicalVolume* MotherVolume, G4Material* TankMat, G4ThreeVector fTankPos, int attribute) {
  fMotherVolume = MotherVolume;
  fTankMat = TankMat;
  if (!fMotherVolume) {
    cout << "in in DSDetectorArDM::addTank(), fMotherVolume = 0. exit." << endl;
    return NULL;
  }

  DSStorage::Get()->SetArDMTankHeight(2 * ArDMvar.TANK_CYLINDER_HALF_HEIGHT);
  DSStorage::Get()->SetArDMTotalTankHeight(2 * ArDMvar.TANK_HALF_HEIGHT);
  DSStorage::Get()->SetArDMTankRadius(ArDMvar.TANK_CYLINDER_OUTER_RADIUS);
  G4Tubs* cylinderSolid = new G4Tubs("tankCylinder", ArDMvar.TANK_CYLINDER_INNER_RADIUS, ArDMvar.TANK_CYLINDER_OUTER_RADIUS, ArDMvar.TANK_CYLINDER_HALF_HEIGHT, 0 * deg, 360 * deg);
  G4Sphere* r1010_arc_solid = new G4Sphere("tankR1010Arc", ArDMvar.TANK_BTM_PART_R1010_ARC_INNER_RADIUS, ArDMvar.TANK_BTM_PART_R1010_ARC_OUTER_RADIUS, 0 * deg, 360 * deg, 180 * deg - ArDMvar.TANK_BTM_PART_R1010_ARC_OPENING_ANGLE / 2, ArDMvar.TANK_BTM_PART_R1010_ARC_OPENING_ANGLE / 2);

  // the connecting part (r101-arc) is complicated !!
  // it's a part of a torus.
  //
  // how to build it logically ?
  // 1. you have an ("solid") arc of a ring with certain opening angle theta,
  // and innerRadius r1, outerRadius r2
  // 2. now you rotate this arc around an axis, which is perpendicular to the
  // normal vector of the arc, you'll get the connecting part. (yea, it's also a
  // ring. or rather a part of a torus (which you'd get if you rotate the whole
  // ring, and not just an arc of it, around a certain axis))
  //
  // in GEANT4, you can generate a whole torus, or a part of a torus with some
  // certain swept angle. but the cross section of the torus is always a ring,
  // and not an arc of a ring. GEANT4 doesn't have any standard class for that
  // shape.
  //
  //--> now, in order to get the connecting part, we have to cut the full torus
  //with some boxes.
  G4double sweptR = (ArDMvar.TANK_CYLINDER_INNER_RADIUS - ArDMvar.TANK_BTM_PART_R101_ARC_INNER_RADIUS);

  G4Torus* torus = new G4Torus("torus", ArDMvar.TANK_BTM_PART_R101_ARC_INNER_RADIUS, ArDMvar.TANK_BTM_PART_R101_ARC_OUTER_RADIUS, sweptR, 0 * deg, 360 * deg);
  // define boxes to cut the torus
  G4Box* upperBox = new G4Box("upperBox", ArDMvar.TANK_CYLINDER_OUTER_RADIUS + 1. * mm, ArDMvar.TANK_CYLINDER_OUTER_RADIUS + 1. * mm, ArDMvar.TANK_BTM_PART_R101_ARC_OUTER_RADIUS / 2 + 1. * mm);
  G4Cons* cone = new G4Cons("cone", ArDMvar.TANK_CYLINDER_INNER_RADIUS - ArDMvar.TANK_BTM_PART_R101_ARC_INNER_RADIUS - ArDMvar.TANK_BTM_PART_R101_ARC_OUTER_RADIUS,
                            ArDMvar.TANK_CYLINDER_INNER_RADIUS - ArDMvar.TANK_BTM_PART_R101_ARC_INNER_RADIUS + ArDMvar.TANK_BTM_PART_R101_ARC_OUTER_RADIUS * sin(ArDMvar.TANK_BTM_PART_R1010_ARC_OPENING_ANGLE / 2),
                            ArDMvar.TANK_CYLINDER_INNER_RADIUS - ArDMvar.TANK_BTM_PART_R101_ARC_INNER_RADIUS - ArDMvar.TANK_BTM_PART_R101_ARC_OUTER_RADIUS, ArDMvar.TANK_CYLINDER_INNER_RADIUS - ArDMvar.TANK_BTM_PART_R101_ARC_INNER_RADIUS,
                            ArDMvar.TANK_BTM_PART_R101_ARC_OUTER_RADIUS * cos(ArDMvar.TANK_BTM_PART_R1010_ARC_OPENING_ANGLE / 2) / 2, 0, 360. * deg);

  G4ThreeVector upperBoxPos(0, 0, upperBox->GetZHalfLength());
  G4ThreeVector conePos(0, 0, -cone->GetZHalfLength());

  G4SubtractionSolid* r101_arc_solid = new G4SubtractionSolid("tankR101Arc", torus, upperBox, 0, upperBoxPos);
  r101_arc_solid = new G4SubtractionSolid("tankR101Arc", r101_arc_solid, cone, 0, conePos);

  // all the parts:
  //  1. tank cylinder
  //  2. r101 connecting part
  //  3. r1010 bottom part
  //--> put these things together to get the shape for the dewar
  G4ThreeVector r101_pos(0, 0, ArDMvar.TANK_BTM_PART_R101_Z_RELATIVE_TO_TANK_CYLINDER);
  G4ThreeVector r1010_pos(0, 0, ArDMvar.TANK_BTM_PART_R1010_Z_RELATIVE_TO_TANK_CYLINDER);

  G4UnionSolid* tankSolid = new G4UnionSolid("tank", cylinderSolid, r101_arc_solid, 0, r101_pos);
  tankSolid = new G4UnionSolid("tank", tankSolid, r1010_arc_solid, 0, r1010_pos);

  G4LogicalVolume* fTankLog = new G4LogicalVolume(tankSolid, fTankMat, "tank");

  G4VisAttributes* tankAtt = new G4VisAttributes(attribute);
  // G4VisAttributes* tankAtt = new G4VisAttributes(false);
  tankAtt->SetColour(1., 1., 1.);  // white
  tankAtt->SetForceAuxEdgeVisible(true);
  fTankLog->SetVisAttributes(tankAtt);

  // put tank into world, positions of r101 and r1010 relative to the center of
  // the tank cylinder
  fTankPhys = new G4PVPlacement(0, fTankPos, "tank", fTankLog, fMotherVolume, false, 0);

  // addTopFlange();
  G4double innerR = ArDMvar.TOP_FLANGE_INNER_RADIUS;
  G4double outerR = ArDMvar.TOP_FLANGE_OUTER_RADIUS;
  G4double halfz = ArDMvar.TOP_FLANGE_HALF_HEIGHT_EFFECTIVE;

  G4Tubs* solid = new G4Tubs("topFlange", innerR, outerR, halfz, 0 * deg, 360 * deg);
  G4LogicalVolume* log = new G4LogicalVolume(solid, fTankMat, "topFlange");

  G4VisAttributes* lidAtt = new G4VisAttributes(true);
  //  G4VisAttributes* lidAtt = new G4VisAttributes(false);
  lidAtt->SetColour(0.0, 1.0, 1.0);  // cyan
  lidAtt->SetForceAuxEdgeVisible(true);
  log->SetVisAttributes(lidAtt);

  //
  G4ThreeVector pos(0, 0, ArDMvar.TOP_FLANGE_POS_Z);

  // G4VPhysicalVolume* fTopFlangePhys = new G4PVPlacement(0, pos, "topFlange",
  // log, fMotherVolume, false, 0); End addtopFlange

  return fTankPhys;
}

G4VPhysicalVolume* DSDetectorArDM::addLArColumn(G4VPhysicalVolume* MotherVolume, G4Material* fMat, G4ThreeVector fPos, int attribute) {
  fMotherVolume = MotherVolume;
  if (!fMotherVolume) {
    cout << "in DSDetectorArDM::addLArColumn: fMotherPhys = 0. exit." << endl;
    return NULL;
  }

  // LArCol = LAr cylinder + LAr btm part
  // LAr cylinder <--> tank_cylinder (or cylinderSolid) in constructTank
  // LAr btm part <--> LAr volume contained in the 'curved' btm part of the
  // tank, corresponding to r101 + r1010 arc together in constructTank the LAr
  // btm part is a union of :
  // 1. a part of a torus (<--> r101 arc)
  // 2. a part of a cone  (<--> r1010 arc + a part of r101 arc)
  G4double LArCylinder_innerR = 0;
  G4double LArCylinder_outerR = ArDMvar.LARCOL_CYLINDER_OUTER_RADIUS;
  G4double LArCylinder_halfz = ArDMvar.LARCOL_CYLINDER_HALFHEIGHT;

  G4Tubs* LArCylinder_solid = new G4Tubs("LArCylinder", LArCylinder_innerR, LArCylinder_outerR, LArCylinder_halfz, 0 * deg, 360 * deg);

  // a part of a torus (r101 arc)
  G4double LAr_r101_arc_innerR = 0;
  G4double LAr_r101_arc_outerR = ArDMvar.TANK_BTM_PART_R101_ARC_INNER_RADIUS;
  G4double sweptR = ArDMvar.TANK_CYLINDER_INNER_RADIUS - ArDMvar.TANK_BTM_PART_R101_ARC_INNER_RADIUS;
  G4Torus* torus = new G4Torus("LArTorus", LAr_r101_arc_innerR, LAr_r101_arc_outerR, sweptR, 0 * deg, 360 * deg);

  // define boxes / cones to cut the torus
  G4Box* upperBox = new G4Box("upperBox", ArDMvar.TANK_CYLINDER_INNER_RADIUS + 1. * mm, ArDMvar.TANK_CYLINDER_INNER_RADIUS + 1. * mm, ArDMvar.TANK_BTM_PART_R101_ARC_INNER_RADIUS / 2 + 1. * mm);

  // the cone used to cut the torus doesn't have to be this big,
  // but just to simplify things, we can also take a big cone.
  G4double cone_upper_innerR = 0;
  G4double cone_upper_outerR = ArDMvar.TANK_CYLINDER_INNER_RADIUS - ArDMvar.TANK_BTM_PART_R101_ARC_INNER_RADIUS;
  G4double cone_lower_innerR = 0;
  G4double cone_lower_outerR = ArDMvar.TANK_BTM_PART_R1010_ARC_OUTER_RADIUS * tan(ArDMvar.TANK_BTM_PART_R1010_ARC_OPENING_ANGLE / 2);
  G4double cone_halfz = (ArDMvar.DISTANCE_LOWER_EDGE_OF_TANK_CYLINDER_TO_TANK_LOWERMOST_POINT + ArDMvar.TANK_CYLINDER_THICKNESS) / 2;

  G4Cons* cone = new G4Cons("cone", cone_lower_innerR, cone_lower_outerR, cone_upper_innerR, cone_upper_outerR, cone_halfz, 0 * deg, 360 * deg);

  G4ThreeVector upperBoxPos(0, 0, upperBox->GetZHalfLength());
  G4ThreeVector conePos(0, 0, -cone->GetZHalfLength());

  G4SubtractionSolid* LAr_r101_solid = new G4SubtractionSolid("LArR101", torus, upperBox, 0, upperBoxPos);
  LAr_r101_solid = new G4SubtractionSolid("LArR101", LAr_r101_solid, cone, 0, conePos);

  // r101 + r1010 joint part
  // actually innerR = 0 mm, but we will later subtract a cone from this sphere,
  // it's better to approximate 0 mm by 0.001 mm
  G4double sphere_innerR = 0.001 * mm;
  G4double sphere_outerR = ArDMvar.TANK_BTM_PART_R1010_ARC_INNER_RADIUS;
  G4Sphere* r1010_r101 = new G4Sphere("sphere", sphere_innerR, sphere_outerR, 0, 360 * deg, 180 * deg - ArDMvar.TANK_BTM_PART_R1010_ARC_OPENING_ANGLE / 2, ArDMvar.TANK_BTM_PART_R1010_ARC_OPENING_ANGLE / 2);

  G4double r1010_r101_cone_lower_innerR = 0;
  G4double r1010_r101_cone_lower_outerR = (ArDMvar.TANK_BTM_PART_R1010_ARC_INNER_RADIUS - ArDMvar.TANK_BTM_PART_R101_ARC_INNER_RADIUS) * sin(ArDMvar.TANK_BTM_PART_R1010_ARC_OPENING_ANGLE / 2);

  G4double r1010_r101_cone_upper_innerR = 0;
  G4double r1010_r101_cone_upper_outerR = sphere_innerR;

  G4double r1010_r101_cone_halfz = (ArDMvar.TANK_BTM_PART_R1010_ARC_INNER_RADIUS - sphere_innerR - ArDMvar.TANK_BTM_PART_R101_ARC_INNER_RADIUS) * cos(ArDMvar.TANK_BTM_PART_R1010_ARC_OPENING_ANGLE / 2);

  G4Cons* r1010_r101_cone = new G4Cons("cone", r1010_r101_cone_lower_innerR, r1010_r101_cone_lower_outerR, r1010_r101_cone_upper_innerR, r1010_r101_cone_upper_outerR, r1010_r101_cone_halfz, 0 * deg, 360 * deg);

  G4ThreeVector r1010_r101_cone_pos(0, 0, -(ArDMvar.DISTANCE_R1010_ARC_CENTER_TO_LOWER_EDGE_OF_TANK_CYLINDER - r1010_r101_cone_halfz - sphere_innerR));
  G4SubtractionSolid* r1010_r101_solid = new G4SubtractionSolid("r1010_r101", r1010_r101, r1010_r101_cone, 0, r1010_r101_cone_pos);

  // put r101 and r1010_r101 joint part together
  G4ThreeVector r101_rel_r1010_r101_pos(0, 0, ArDMvar.DISTANCE_R1010_ARC_CENTER_TO_LOWER_EDGE_OF_TANK_CYLINDER);
  G4UnionSolid* btm_part_solid = new G4UnionSolid("btmPartSolid", LAr_r101_solid, r1010_r101_solid, 0, r101_rel_r1010_r101_pos);

  // put all pieces together to get the whole LAr volume
  // position vector relative to the center of LAr_cylinder_solid
  G4ThreeVector btmPart_pos(0, 0, -ArDMvar.LARCOL_CYLINDER_HALFHEIGHT);
  G4UnionSolid* LArCol_solid = new G4UnionSolid("LArCol_solid", LArCylinder_solid, btm_part_solid, 0, btmPart_pos);

  G4LogicalVolume* fLArColLog = new G4LogicalVolume(LArCol_solid, fMat, "LArCol");
  // G4VisAttributes* detAtt = new G4VisAttributes(false);
  G4VisAttributes* detAtt = new G4VisAttributes(attribute);
  detAtt->SetColour(1.0, 1.0, 0.0);  // yellow
  detAtt->SetForceAuxEdgeVisible(true);
  fLArColLog->SetVisAttributes(detAtt);
  // Named ActiveLAr mandatory for energy purpose
  // fLArColPhys = new
  // G4PVPlacement(0,fPos,"LArCol",fLArColLog,fMotherVolume,false,0);
  fLArColPhys = new G4PVPlacement(0, fPos, "ActiveLAr", fLArColLog, fMotherVolume, false, 0);

  return fLArColPhys;
}

G4VPhysicalVolume* DSDetectorArDM::addGArColumn(G4VPhysicalVolume* MotherVolume, G4Material* fMat, G4ThreeVector fPos, int attribute) {
  fMotherVolume = MotherVolume;
  if (!fMotherVolume) {
    cout << "in DSDetectorArDM::addGArColumn fMotherVolume = 0. exit." << endl;
    return NULL;
  }

  G4double z = ArDMvar.GAR_COLUMN_Z;
  G4ThreeVector pos(0., 0., z);

  // add detector (active volume) into tank
  G4double det_innerR = 0;
  G4double det_outerR = ArDMvar.TANK_CYLINDER_INNER_RADIUS;
  G4double det_half_height = ArDMvar.GAR_COLUMN_HALF_HEIGHT;

  G4Tubs* detSolid = new G4Tubs("GArCol", det_innerR, det_outerR, det_half_height, 0. * deg, 360. * deg);
  G4LogicalVolume* fGArColLog = new G4LogicalVolume(detSolid, fMat, "GArCol");

  G4VisAttributes* detAtt = new G4VisAttributes(attribute);
  // G4VisAttributes* detAtt     = new G4VisAttributes(false);
  detAtt->SetColour(1.0, 1.0, 1.0);  // white
  detAtt->SetForceAuxEdgeVisible(true);
  fGArColLog->SetVisAttributes(detAtt);

  fGArColPhys = new G4PVPlacement(0, fPos, "GArCol", fGArColLog, fMotherVolume, false, 0);

  return fGArColPhys;
}

//// Edgar DART new geometry two pipes up

G4VPhysicalVolume* DSDetectorArDM::addTestArgon(G4VPhysicalVolume* MotherVolume, G4Material* fMat, G4ThreeVector fPos, int attribute) {
  fMotherVolume = MotherVolume;
  if (!fMotherVolume) {
    cout << "in DSDetectorArDM::addTestArgon fMotherVolume = 0. exit." << endl;
    return NULL;
  }

  G4double det_innerR = 0.;
  G4double thickness_wall = 5. * mm;
  G4double thickness_flange = 8. * mm;

  G4double det_outerR = 65. * mm;
  G4double det_outerR_flange = 78. * mm;
  G4double det_half_height = 120. * mm;
  G4double det_inner_pipe = 6.35 * mm;
  G4double det_length_pipe = 150. * mm;

  G4double det_outerR_Acrylic = 41.6 * mm;
  G4double det_half_height_Acrylic = 100. * mm;

  G4double det_outerR_Reflector = 42.1 * mm;
  G4double det_half_height_Reflector = 100. * mm;

  G4double det_outerR_WLS = 35.6 * mm;

  G4double det_outerR_Ar = 35.5 * mm;
  G4double det_half_height_Ar = 100. * mm;

  G4double thickness_disk = 3. * mm;
  G4double thickness_TPB = 0.05 * mm;
  G4double thickness_reflector = 2. * mm;

  //// SiPM planes

  G4double x_axis = 25. * mm;
  G4double y_axis = 25. * mm;
  G4double thickness = 1. * mm;
  G4double half_height = 100. * mm;

  G4double x_axis_PCB = 35.5 * mm;
  G4double y_axis_PCB = 35.5 * mm;
  G4double thickness_PCB = 3. * mm;
  G4double half_height_PCB = 115. * mm;

  DSStorage::Get()->SetArDMTestArgonHeight(2 * (det_half_height + thickness_flange) + det_length_pipe);
  DSStorage::Get()->SetArDMTestArgonRadius(det_outerR_flange);

  G4Tubs* detCopperFlange = new G4Tubs("CopperFlangeTestArgon_solid", det_innerR, det_outerR_flange, thickness_flange, 0. * deg, 360. * deg);
  G4Tubs* detCopperFlangeDown = new G4Tubs("CopperFlangeDownTestArgon_solid", det_innerR, det_outerR + thickness_wall, thickness_flange / 2., 0. * deg, 360. * deg);
  G4Tubs* detCopperCylinder = new G4Tubs("CopperTestArgon_Cylinder", det_innerR, det_outerR + thickness_wall, det_half_height, 0. * deg, 360. * deg);

  G4Tubs* detPipeUp = new G4Tubs("PipeUp", det_inner_pipe, thickness_wall / 2. + det_inner_pipe, det_length_pipe, 0. * deg, 360. * deg);

  G4Tubs* detPipeDown = new G4Tubs("PipeUp", det_inner_pipe, thickness_wall / 2. + det_inner_pipe, det_length_pipe, 0. * deg, 360. * deg);

  G4ThreeVector det_pos(0, 0, det_half_height);
  G4ThreeVector det_pos_Down(0, 0, -(det_half_height + thickness_flange / 2.));
  G4ThreeVector det_pipe_up(-(det_outerR - 2.5 * thickness_wall), 0, det_half_height);
  G4ThreeVector det_pipe_down((det_outerR - 2.5 * thickness_wall), 0, det_half_height);

  G4UnionSolid* detCopperSolid_1 = new G4UnionSolid("CopperTestArgon_solid_1", detCopperCylinder, detCopperFlange, 0, det_pos);
  G4UnionSolid* detCopperSolid_2 = new G4UnionSolid("CopperTestArgon_solid_2", detCopperSolid_1, detCopperFlangeDown, 0, det_pos_Down);
  G4UnionSolid* detCopperSolid_3 = new G4UnionSolid("CopperTestArgon_solid_3", detCopperSolid_2, detPipeUp, 0, det_pipe_up);
  G4UnionSolid* detCopperSolid = new G4UnionSolid("CopperTestArgon_solid", detCopperSolid_3, detPipeDown, 0, det_pipe_down);

  G4LogicalVolume* fCopperTestArgonLog = new G4LogicalVolume(detCopperSolid, DSMaterial::Get()->GetMetalCopper(), "CopperTestArgon_log");
  fCopperTestArgonPhys = new G4PVPlacement(0, fPos, "CopperTestArgon_phys", fCopperTestArgonLog, fMotherVolume, false, 0);

  G4VisAttributes* detAttCopper = new G4VisAttributes(attribute);
  detAttCopper->SetColour(0.0, 1.0, 0.0);  // green
  detAttCopper->SetForceAuxEdgeVisible(true);
  fCopperTestArgonLog->SetVisAttributes(detAttCopper);

  G4Tubs* BatchArgonSolid = new G4Tubs("BatchArgon_solid", det_innerR, det_outerR, det_half_height, 0. * deg, 360. * deg);
  G4LogicalVolume* fBatchArgonLog = new G4LogicalVolume(BatchArgonSolid, DSMaterial::Get()->GetNSLiquidArgon(), "BatchTestArgon_log");

  G4VisAttributes* BatchArgonAtt = new G4VisAttributes(attribute);
  BatchArgonAtt->SetColour(0.5, 0.5, 0.5);  // grey
  BatchArgonAtt->SetForceAuxEdgeVisible(true);
  fBatchArgonLog->SetVisAttributes(BatchArgonAtt);

  G4ThreeVector pos(0., 0., 0.);
  fBatchArgonPhys = new G4PVPlacement(0, pos, "BatchTestArgon_phys", fBatchArgonLog, fCopperTestArgonPhys, false, 0);

  /// External acrylic inside DART

  G4Tubs* ExtAcrylicSolid = new G4Tubs("ExtAcrylic_solid", det_outerR_Reflector, det_outerR - 2. * mm, det_half_height, 0. * deg, 360. * deg);
  G4LogicalVolume* fExtAcrylicLog = new G4LogicalVolume(ExtAcrylicSolid, DSMaterial::Get()->GetAcrylic(), "ExtAcrylic_log");
  G4VisAttributes* ExtAcrylicAtt = new G4VisAttributes(attribute);
  ExtAcrylicAtt->SetColour(0.0, 3.0, 0.0);
  ExtAcrylicAtt->SetForceAuxEdgeVisible(true);
  fExtAcrylicLog->SetVisAttributes(ExtAcrylicAtt);

  // G4VPhysicalVolume* fExtAcrylicPhys =
  //     new G4PVPlacement(0, pos, "ExtAcrylic_phys", fExtAcrylicLog,
  //     fBatchArgonPhys, false, 0);

  /// PCB

  G4ThreeVector moveup_PCB(0, 0, half_height_PCB - thickness_PCB);

  G4Box* PCB = new G4Box("PCB", x_axis_PCB, y_axis_PCB, thickness_PCB);

  G4LogicalVolume* fPCBLog = new G4LogicalVolume(PCB, DSMaterial::Get()->GetAcrylic(), "PCB");

  G4VisAttributes* PCBAtt = new G4VisAttributes(true);
  PCBAtt->SetColour(1.0, 1.0, 1.0);
  PCBAtt->SetForceAuxEdgeVisible(true);
  fPCBLog->SetVisAttributes(PCBAtt);

  // G4VPhysicalVolume* fPCBup = new G4PVPlacement(0, moveup_PCB, "PCB",
  // fPCBLog, fBatchArgonPhys, false, 0); G4VPhysicalVolume* fPCBdw = new
  // G4PVPlacement(0, -moveup_PCB, "PCB", fPCBLog, fBatchArgonPhys, false, 1);

  /////////Pipe inside

  G4Tubs* PipeInsideSolid = new G4Tubs("PipeInsideTestArgon_solid", det_inner_pipe, thickness_wall / 2. + det_inner_pipe, det_half_height_Acrylic + thickness_flange + thickness_wall, 0. * deg, 360. * deg);
  G4LogicalVolume* fPipeInsideTestArgonLog = new G4LogicalVolume(PipeInsideSolid, DSMaterial::Get()->GetMetalCopper(), "PipeInsideTestArgon_log");

  G4VisAttributes* detAttPipeInside = new G4VisAttributes(attribute);

  detAttPipeInside->SetColour(0.0, 1.0, 0.0);  //
  detAttPipeInside->SetForceAuxEdgeVisible(true);
  fPipeInsideTestArgonLog->SetVisAttributes(detAttPipeInside);
  G4ThreeVector det_pipe_inside((det_outerR - 2.5 * thickness_wall), 0, thickness_flange + thickness_wall * 2.);

  // G4VPhysicalVolume* fPipeInsideTestArgonPhys = new G4PVPlacement(0,
  // det_pipe_inside, "PipeInsideTestArgon_phys",
  //                                                                 fPipeInsideTestArgonLog,
  //                                                                 fExtAcrylicPhys,
  //                                                                 false, 0);

  /// Reflector

  G4Tubs* ReflectorSolid = new G4Tubs("ReflectorTestArgon_solid", det_innerR, det_outerR_Reflector, det_half_height_Reflector + thickness_reflector, 0. * deg, 360. * deg);
  G4LogicalVolume* fReflectorTestArgonLog = new G4LogicalVolume(ReflectorSolid, DSMaterial::Get()->GetAluminum(), "ReflectorTestArgon_log");

  G4VisAttributes* detAttReflector = new G4VisAttributes(attribute);
  detAttReflector->SetColour(7.0, 0.0, 0.0);  //
  detAttReflector->SetForceAuxEdgeVisible(true);
  fReflectorTestArgonLog->SetVisAttributes(detAttReflector);

  fReflectorTestArgonPhys = new G4PVPlacement(0, pos, "ReflectorTestArgon_phys", fReflectorTestArgonLog, fBatchArgonPhys, false, 0);

  //////////Acrylic

  G4Tubs* AcrylicSolid = new G4Tubs("AcrylicTestArgon_solid", det_innerR, det_outerR_Acrylic, det_half_height_Acrylic, 0. * deg, 360. * deg);
  G4LogicalVolume* fAcrylicTestArgonLog = new G4LogicalVolume(AcrylicSolid, DSMaterial::Get()->GetAcrylicDART(), "AcrylicTestArgon_log");

  G4VisAttributes* detAttAcrylic = new G4VisAttributes(attribute);
  detAttAcrylic->SetColour(0.0, 3.0, 0.0);  //
  detAttAcrylic->SetForceAuxEdgeVisible(true);
  fAcrylicTestArgonLog->SetVisAttributes(detAttAcrylic);

  fAcrylicTestArgonPhys = new G4PVPlacement(0, pos, "AcrylicTestArgon_phys", fAcrylicTestArgonLog, fReflectorTestArgonPhys, false, 0);

  ///////////TPB

  G4Tubs* WLSSolid = new G4Tubs("WLSTestArgon_solid", det_outerR_Ar, det_outerR_WLS, half_height - 4 * thickness - 2 * thickness_disk, 0. * deg, 360. * deg);
  G4LogicalVolume* fWLSTestArgonLog = new G4LogicalVolume(WLSSolid, DSMaterial::Get()->GetTPB(), "WLSTestArgon_log");

  G4VisAttributes* detAttWLS = new G4VisAttributes(true);
  detAttWLS->SetColour(0.0, 5.0, 0.0);
  detAttWLS->SetForceAuxEdgeVisible(true);
  fWLSTestArgonLog->SetVisAttributes(detAttWLS);

  fWLSTestArgonPhys = new G4PVPlacement(0, pos, "WLSTestArgon_phys", fWLSTestArgonLog, fAcrylicTestArgonPhys, false, 0);

  /// Underground Argon

  G4Tubs* DeepArgonSolid = new G4Tubs("DeepArgonTestArgon_solid", det_innerR, det_outerR_Ar, det_half_height_Ar, 0. * deg, 360. * deg);
  G4LogicalVolume* fDeepArgonTestArgonLog = new G4LogicalVolume(DeepArgonSolid, fMat, "DeepArgonTestArgon_log");

  G4VisAttributes* detAttDeepArgon = new G4VisAttributes(true);
  detAttDeepArgon->SetColour(0.0, 6.0, 0.0);  //
  detAttDeepArgon->SetForceAuxEdgeVisible(true);
  fDeepArgonTestArgonLog->SetVisAttributes(detAttDeepArgon);

  fDeepArgonTestArgonPhys = new G4PVPlacement(0, pos, "ActiveLAr", fDeepArgonTestArgonLog, fAcrylicTestArgonPhys, false, 0);

  //// SiPM planes

  G4ThreeVector moveup(0, 0, half_height - thickness);
  G4Material* fSiPMMat = DSMaterial::Get()->GetMetalSilicon();
  G4Box* sipm = new G4Box("SiPM", x_axis, y_axis, thickness);

  G4LogicalVolume* fSiPMLog = new G4LogicalVolume(sipm, fSiPMMat, "SiPM");

  G4VisAttributes* solidAtt = new G4VisAttributes(true);
  solidAtt->SetColour(1.0, 0.0, 0.0);
  solidAtt->SetForceAuxEdgeVisible(true);
  fSiPMLog->SetVisAttributes(solidAtt);

  fSiPMup = new G4PVPlacement(0, moveup, "SiPM", fSiPMLog, fDeepArgonTestArgonPhys, false, 0);
  fSiPMdw = new G4PVPlacement(0, -moveup, "SiPM", fSiPMLog, fDeepArgonTestArgonPhys, false, 1);

  /// Acrylic Disks

  G4ThreeVector moveup_disk(0, 0, half_height - 4 * thickness - thickness_disk);
  G4Tubs* AcrylicDiskSolid = new G4Tubs("AcrylicDisk_solid", det_innerR, det_outerR_Ar, thickness_disk, 0. * deg, 360. * deg);

  G4LogicalVolume* fAcrylicDiskLog = new G4LogicalVolume(AcrylicDiskSolid, DSMaterial::Get()->GetAcrylicDART(), "AcrylicDisk");

  G4VisAttributes* AcrylicDiskAtt = new G4VisAttributes(true);
  AcrylicDiskAtt->SetColour(1.0, 0.0, 0.0);
  AcrylicDiskAtt->SetForceAuxEdgeVisible(true);
  fAcrylicDiskLog->SetVisAttributes(AcrylicDiskAtt);

  fDiskup = new G4PVPlacement(0, moveup_disk, "AcrylicDisk", fAcrylicDiskLog, fDeepArgonTestArgonPhys, false, 0);
  fDiskdw = new G4PVPlacement(0, -moveup_disk, "AcrylicDisk", fAcrylicDiskLog, fDeepArgonTestArgonPhys, false, 1);

  /// TPB Disks

  G4ThreeVector moveup_diskTPB(0, 0, half_height - 4 * thickness - 2 * thickness_disk - thickness_TPB);
  G4Tubs* TPBDiskSolid = new G4Tubs("TPBDisk_solid", det_innerR, det_outerR_Ar, thickness_TPB, 0. * deg, 360. * deg);

  G4LogicalVolume* fTPBDiskLog = new G4LogicalVolume(TPBDiskSolid, DSMaterial::Get()->GetTPB(), "TPBDisk");

  G4VisAttributes* TPBDiskAtt = new G4VisAttributes(true);
  TPBDiskAtt->SetColour(1.0, 0.0, 0.0);
  TPBDiskAtt->SetForceAuxEdgeVisible(true);
  fTPBDiskLog->SetVisAttributes(TPBDiskAtt);

  fDiskupTPB = new G4PVPlacement(0, moveup_diskTPB, "AcrylicDisk", fTPBDiskLog, fDeepArgonTestArgonPhys, false, 0);
  fDiskdwTPB = new G4PVPlacement(0, -moveup_diskTPB, "AcrylicDisk", fTPBDiskLog, fDeepArgonTestArgonPhys, false, 1);

  DSStorage::Get()->SetPMTMaterialIndex(fSiPMLog->GetMaterial()->GetIndex());

  ////////////////////////////
  ////DART Optical simulaton///
  ////////////////////////////

  // UAr <-> WLS
  G4OpticalSurface* fOpUArWLSSurface = new G4OpticalSurface("OpUArWLSSurface");
  fOpUArWLSSurface->SetType(dielectric_dielectric);
  fOpUArWLSSurface->SetModel(unified);
  fOpUArWLSSurface->SetFinish(ground);
  fOpUArWLSSurface->SetSigmaAlpha(0.3);
  new G4LogicalBorderSurface("UArWLSSurface", fDeepArgonTestArgonPhys, fWLSTestArgonPhys, fOpUArWLSSurface);
  new G4LogicalBorderSurface("WLSUArSurface", fWLSTestArgonPhys, fDeepArgonTestArgonPhys, fOpUArWLSSurface);
  /// Copy from DarkSide-20k
  G4double VISTRAN = DSParameters::Get()->GetArTPBVisTran();

  const G4int NUM = 4;
  G4double pp[NUM] = {0.1 * eV, 8.0 * eV, 8.3 * eV, 20.0 * eV};  // vis, vis, UV, UV
  G4double specularlobe[NUM] = {0., 0., 0., 0.};                 //--
  G4double specularspike[NUM] = {0., 0., 0., 0.};                //----  gives all reflection to Lambertian lobe
  G4double backscatter[NUM] = {0., 0., 0., 0.};                  //--
  G4double reflectivity[NUM] = {1.0, 1.0, 1.0, 1.0};             //  To set 1-absorption
  G4double transmitivity[NUM] = {VISTRAN, VISTRAN, 1.0, 1.0};    //  To set reflection vs. transmission, overridding Fresnel
                                                                 //  For now, no angle dependence.
  G4MaterialPropertiesTable* fUArWLSSurfProp = new G4MaterialPropertiesTable();
  fUArWLSSurfProp->AddProperty("SPECULARLOBECONSTANT", pp, specularlobe, NUM);
  fUArWLSSurfProp->AddProperty("SPECULARSPIKECONSTANT", pp, specularspike, NUM);
  fUArWLSSurfProp->AddProperty("BACKSCATTERCONSTANT", pp, backscatter, NUM);
  fUArWLSSurfProp->AddProperty("REFLECTIVITY", pp, reflectivity, NUM);
  fUArWLSSurfProp->AddProperty("TRANSMITTANCE", pp, transmitivity, NUM);
  fUArWLSSurfProp->AddConstProperty("DOArTPB", 1);
  fOpUArWLSSurface->SetMaterialPropertiesTable(fUArWLSSurfProp);

  /// UAR<-> Disk WLS

  G4OpticalSurface* fOpUArdiskWLSSurface = new G4OpticalSurface("OpUArdiskWLSSurface");
  fOpUArdiskWLSSurface->SetType(dielectric_dielectric);
  fOpUArdiskWLSSurface->SetModel(unified);
  fOpUArdiskWLSSurface->SetFinish(ground);
  fOpUArdiskWLSSurface->SetSigmaAlpha(0.3);
  new G4LogicalBorderSurface("UArdiskWLSSurface", fDeepArgonTestArgonPhys, fDiskupTPB, fOpUArWLSSurface);
  new G4LogicalBorderSurface("diskWLSUArSurface", fDiskupTPB, fDeepArgonTestArgonPhys, fOpUArWLSSurface);
  new G4LogicalBorderSurface("UArdiskWLSSurface", fDeepArgonTestArgonPhys, fDiskdwTPB, fOpUArWLSSurface);
  new G4LogicalBorderSurface("diskWLSUArSurface", fDiskdwTPB, fDeepArgonTestArgonPhys, fOpUArWLSSurface);

  ////WLS -> disk Acrylic

  G4OpticalSurface* fOpWLSdiskAcrylicSurface = new G4OpticalSurface("OpWLSdiskAcrylicSurface");
  fOpWLSdiskAcrylicSurface->SetType(dielectric_dielectric);
  fOpWLSdiskAcrylicSurface->SetModel(unified);
  fOpWLSdiskAcrylicSurface->SetFinish(polished);
  new G4LogicalBorderSurface("WLSdiskAcrylicSurface", fDiskupTPB, fDiskup, fOpWLSdiskAcrylicSurface);
  new G4LogicalBorderSurface("WLSdiskAcrylicSurface", fDiskdwTPB, fDiskdw, fOpWLSdiskAcrylicSurface);

  //  G4double reflectivity_acrylic[NUM] = {0.0, 0.0, 0.0, 0.0};     //  To set
  //  1-absorption G4double transmitivity_acrylic[NUM] = {1.0, 1.0, 1.0, 1.0};
  G4MaterialPropertiesTable* fTPBAcrylicSurfProp = new G4MaterialPropertiesTable();
  //  fTPBAcrylicSurfProp ->
  //  AddProperty("REFLECTIVITY",pp,reflectivity_acrylic,NUM);
  //  fTPBAcrylicSurfProp ->
  //  AddProperty("TRANSMITTANCE",pp,transmitivity_acrylic,NUM);
  fOpWLSdiskAcrylicSurface->SetMaterialPropertiesTable(fTPBAcrylicSurfProp);

  ////disk Acrylic -> UAr

  G4OpticalSurface* fOpdiskAcrylicUArSurface = new G4OpticalSurface("OpdiskAcrylicUArSurface");
  fOpdiskAcrylicUArSurface->SetType(dielectric_dielectric);
  fOpdiskAcrylicUArSurface->SetModel(unified);
  fOpdiskAcrylicUArSurface->SetFinish(polished);
  new G4LogicalBorderSurface("diskAcrylicUArSurface", fDiskup, fDeepArgonTestArgonPhys, fOpdiskAcrylicUArSurface);
  new G4LogicalBorderSurface("diskAcrylicUArSurface", fDiskdw, fDeepArgonTestArgonPhys, fOpdiskAcrylicUArSurface);

  //  G4double reflectivity_acrylic_UAr[NUM] = {0.0, 0.0, 0.0, 0.0};     //  To
  //  set 1-absorption G4double transmitivity_acrylic_UAr[NUM] =
  //  {1.0, 1.0, 1.0, 1.0};
  G4MaterialPropertiesTable* fAcrylicUArSurfProp = new G4MaterialPropertiesTable();
  //  fAcrylicUArSurfProp ->
  //  AddProperty("REFLECTIVITY",pp,reflectivity_acrylic_UAr,NUM);
  //  fAcrylicUArSurfProp ->
  //  AddProperty("TRANSMITTANCE",pp,transmitivity_acrylic_UAr,NUM);
  fOpdiskAcrylicUArSurface->SetMaterialPropertiesTable(fAcrylicUArSurfProp);

  //// UAr <->Acrylic

  G4OpticalSurface* fOpAcrylicUArSurface = new G4OpticalSurface("OAcrylicUArSurface");
  fOpAcrylicUArSurface->SetType(dielectric_dielectric);
  fOpAcrylicUArSurface->SetModel(unified);
  fOpAcrylicUArSurface->SetFinish(polished);
  new G4LogicalBorderSurface("AcrylicUArSurface", fAcrylicTestArgonPhys, fDeepArgonTestArgonPhys, fOpAcrylicUArSurface);
  new G4LogicalBorderSurface("UArAcrylicSurface", fDeepArgonTestArgonPhys, fAcrylicTestArgonPhys, fOpAcrylicUArSurface);

  //  G4double reflectivity_acrylic_UAr[NUM] = {0.0, 0.0, 0.0, 0.0};     //  To
  //  set 1-absorption G4double transmitivity_acrylic_UAr[NUM] =
  //  {1.0, 1.0, 1.0, 1.0};

  //  fAcrylicUArSurfProp ->
  //  AddProperty("REFLECTIVITY",pp,reflectivity_acrylic_UAr,NUM);
  //  fAcrylicUArSurfProp ->
  //  AddProperty("TRANSMITTANCE",pp,transmitivity_acrylic_UAr,NUM);
  fOpAcrylicUArSurface->SetMaterialPropertiesTable(fAcrylicUArSurfProp);

  ////WLS <-> Acrylic

  G4OpticalSurface* fOpWLSAcrylicSurface = new G4OpticalSurface("OpWLSAcrylicSurface");
  fOpWLSAcrylicSurface->SetType(dielectric_dielectric);
  fOpWLSAcrylicSurface->SetModel(unified);
  fOpWLSAcrylicSurface->SetFinish(polished);

  new G4LogicalBorderSurface("WLSAcrylicSurface", fWLSTestArgonPhys, fAcrylicTestArgonPhys, fOpWLSAcrylicSurface);

  G4OpticalSurface* fOpAcrylicWLSSurface = new G4OpticalSurface("OpAcrylicWLSSurface");
  fOpAcrylicWLSSurface->SetType(dielectric_dielectric);
  fOpAcrylicWLSSurface->SetModel(unified);
  fOpAcrylicWLSSurface->SetFinish(polished);

  new G4LogicalBorderSurface("AcrylicWLSSurface", fAcrylicTestArgonPhys, fWLSTestArgonPhys, fOpAcrylicWLSSurface);

  fOpWLSAcrylicSurface->SetMaterialPropertiesTable(fTPBAcrylicSurfProp);
  fOpAcrylicWLSSurface->SetMaterialPropertiesTable(fTPBAcrylicSurfProp);

  /// Acrylic -> Reflector

  G4double TeflonTPBENE[4] = {0.1 * eV, 8.0 * eV, 8.3 * eV, 20.0 * eV};
  G4double TeflonTPBREF[4] = {1., 1., 1., 1.};
  G4OpticalSurface* fOpAcrylicAlSurface = new G4OpticalSurface("OpAcrylicAlSurface");
  fOpAcrylicAlSurface->SetType(dielectric_metal);
  fOpAcrylicAlSurface->SetModel(unified);
  fOpAcrylicAlSurface->SetFinish(polished);
  new G4LogicalBorderSurface("AcrylicAlSurface", fAcrylicTestArgonPhys, fReflectorTestArgonPhys, fOpAcrylicAlSurface);

  G4MaterialPropertiesTable* fAcrylicAlSurfProp = new G4MaterialPropertiesTable();
  fAcrylicAlSurfProp->AddProperty("REFLECTIVITY", TeflonTPBENE, TeflonTPBREF, 4);
  fOpAcrylicAlSurface->SetMaterialPropertiesTable(fAcrylicAlSurfProp);

  /// UAr -> Reflector

  G4OpticalSurface* fOpArReflectorSurface = new G4OpticalSurface("OpArReflectorSurface");
  fOpArReflectorSurface->SetType(dielectric_metal);
  fOpArReflectorSurface->SetModel(unified);
  fOpArReflectorSurface->SetFinish(polished);
  new G4LogicalBorderSurface("ArReflectorSurface", fDeepArgonTestArgonPhys, fReflectorTestArgonPhys, fOpArReflectorSurface);

  G4MaterialPropertiesTable* fArReflectorSurfProp = new G4MaterialPropertiesTable();
  fArReflectorSurfProp->AddProperty("REFLECTIVITY", TeflonTPBENE, TeflonTPBREF, 4);
  fOpArReflectorSurface->SetMaterialPropertiesTable(fArReflectorSurfProp);

  /// UAr -> SiPM
  /*
  G4double SiPMREF[4] = {0., 0. ,0. , 0. };
  G4double SiPMTRN[4] = {1, 1 ,1 , 1 };

  G4OpticalSurface *fOpUArSiPMSurface = new
  G4OpticalSurface("OpUArSiPMSurface"); fOpUArSiPMSurface->SetType(
  dielectric_dielectric ); fOpUArSiPMSurface->SetModel(unified);
  fOpUArSiPMSurface->SetFinish(groundfrontpainted);
  fOpUArSiPMSurface->SetSigmaAlpha(0.1);
  new G4LogicalBorderSurface("UArSiPMSurface", fDeepArgonTestArgonPhys, fSiPMup,
  fOpUArSiPMSurface ); new G4LogicalBorderSurface("UArSiPMSurface",
  fDeepArgonTestArgonPhys, fSiPMdw, fOpUArSiPMSurface );

  G4MaterialPropertiesTable *fUArSiPMSurfProp2 = new
  G4MaterialPropertiesTable(); fUArSiPMSurfProp2->AddProperty("REFLECTIVITY",
  pp, SiPMREF, 4); fUArSiPMSurfProp2->AddProperty("TRANSMITTANCE", pp, SiPMTRN,
  4); fOpUArSiPMSurface->SetMaterialPropertiesTable( fUArSiPMSurfProp2 );
*/

  return fDeepArgonTestArgonPhys;  // fTestArgonPhys;
}

void DSDetectorArDM::addWLS() {
  // ring sector
  G4double Ring_innerR = ArDMvar.WLS_RINGSEC_INNER_RADIUS;
  G4double Ring_outerR = ArDMvar.WLS_RINGSEC_OUTER_RADIUS;
  G4double Ring_start_phi = ArDMvar.WLS_RINGSEC_START_PHI;
  G4double Ring_delta_phi = ArDMvar.WLS_RINGSEC_DELTA_PHI;

  G4double Ring_half_height_GAr = ArDMvar.WLS_RINGSEC_HALF_HEIGHT_GAR;
  G4double Ring_half_height_LAr = ArDMvar.WLS_RINGSEC_HALF_HEIGHT_LAR;

  G4Tubs* RingSolid_GAr = NULL;
  G4Tubs* RingSolid_LAr = NULL;

  if (Ring_half_height_GAr > 0) { RingSolid_GAr = new G4Tubs("wlsRing", Ring_innerR, Ring_outerR, Ring_half_height_GAr, Ring_start_phi, Ring_delta_phi); }

  if (Ring_half_height_LAr > 0) { RingSolid_LAr = new G4Tubs("wlsRing", Ring_innerR, Ring_outerR, Ring_half_height_LAr, Ring_start_phi, Ring_delta_phi); }

  // linear sector
  G4double Linsec_halfx = ArDMvar.WLS_LINSEC_HALF_X;
  G4double Linsec_halfy = ArDMvar.WLS_LINSEC_HALF_Y;
  G4double Linsec_halfz_GAr = ArDMvar.WLS_SUPPORT_LINSEC_HALF_HEIGHT_GAR;
  G4double Linsec_halfz_LAr = ArDMvar.WLS_SUPPORT_LINSEC_HALF_HEIGHT_LAR;

  G4Box* LinSecSolid_GAr = NULL;
  G4Box* LinSecSolid_LAr = NULL;

  if (Ring_half_height_GAr > 0) { LinSecSolid_GAr = new G4Box("wlsLinSec", Linsec_halfx, Linsec_halfy, Linsec_halfz_GAr); }
  if (Ring_half_height_LAr > 0) { LinSecSolid_LAr = new G4Box("wlsLinSec", Linsec_halfx, Linsec_halfy, Linsec_halfz_LAr); }

  // position of the linear sector relative to the ring sector
  G4RotationMatrix* Linsec_rotMat = new G4RotationMatrix;
  Linsec_rotMat->rotateX(0);
  Linsec_rotMat->rotateY(0);
  Linsec_rotMat->rotateZ(90. * deg);

  G4ThreeVector pos(ArDMvar.WLS_LINSEC_POS_X, ArDMvar.WLS_LINSEC_POS_Y,
                    0.);  // relative to the ring sector

  // build unionSolid of ring and linear sector
  G4UnionSolid* wls_GAr = NULL;
  G4UnionSolid* wls_LAr = NULL;

  if (RingSolid_GAr && LinSecSolid_GAr) {
    // wls_GAr = new
    // G4UnionSolid("WLS",RingSolid_GAr,LinSecSolid_LAr,Linsec_rotMat,pos); <---
    // should be a bug or not
    wls_GAr = new G4UnionSolid("WLS", RingSolid_GAr, LinSecSolid_GAr, Linsec_rotMat, pos);
  }
  if (RingSolid_LAr && LinSecSolid_LAr) { wls_LAr = new G4UnionSolid("WLS", RingSolid_LAr, LinSecSolid_LAr, Linsec_rotMat, pos); }

  G4LogicalVolume* wlsLog_GAr = NULL;
  G4LogicalVolume* wlsLog_LAr = NULL;

  fWLSMat = DSMaterial::Get()->GetTPB();

  if (wls_GAr) { wlsLog_GAr = new G4LogicalVolume(wls_GAr, fWLSMat, "WLS", 0, 0, 0); }
  if (wls_LAr) { wlsLog_LAr = new G4LogicalVolume(wls_LAr, fWLSMat, "WLS", 0, 0, 0); }

  G4VisAttributes* wlsAtt = new G4VisAttributes(true);
  // G4VisAttributes* wlsAtt = new G4VisAttributes(false);
  wlsAtt->SetColour(1.0, .0, 1.0);  // magenta
  wlsAtt->SetForceAuxEdgeVisible(true);

  if (wlsLog_GAr) { wlsLog_GAr->SetVisAttributes(wlsAtt); }
  if (wlsLog_LAr) { wlsLog_LAr->SetVisAttributes(wlsAtt); }

  G4ThreeVector wlsPos_GAr(0, 0, ArDMvar.WLS_RINGSEC_GAR_POS_Z_GARCOL);
  G4ThreeVector wlsPos_LAr(0, 0, ArDMvar.WLS_RINGSEC_LAR_POS_Z_LARCOL);

  fWLSGArPhys = NULL;
  fWLSLArPhys = NULL;

  if (wlsLog_GAr && fGArColPhys) { fWLSGArPhys = new G4PVPlacement(0, wlsPos_GAr, "WLS", wlsLog_GAr, fGArColPhys, false, 0); }

  if (wlsLog_LAr && fLArColPhys) { fWLSLArPhys = new G4PVPlacement(0, wlsPos_LAr, "WLS", wlsLog_LAr, fLArColPhys, false, 0); }

  return;
}

void DSDetectorArDM::addWLSSupport() {
  // ring sector
  G4double Ring_innerR = ArDMvar.WLS_SUPPORT_RINGSEC_INNER_RADIUS;
  G4double Ring_outerR = ArDMvar.WLS_SUPPORT_RINGSEC_OUTER_RADIUS;
  G4double Ring_start_phi = ArDMvar.WLS_SUPPORT_RINGSEC_START_PHI;
  G4double Ring_delta_phi = ArDMvar.WLS_SUPPORT_RINGSEC_DELTA_PHI;

  G4double Ring_half_height_GAr = ArDMvar.WLS_SUPPORT_RINGSEC_HALF_HEIGHT_GAR;
  G4double Ring_half_height_LAr = ArDMvar.WLS_SUPPORT_RINGSEC_HALF_HEIGHT_LAR;

  G4Tubs* RingSolid_GAr = NULL;
  G4Tubs* RingSolid_LAr = NULL;

  if (Ring_half_height_GAr > 0) { RingSolid_GAr = new G4Tubs("wlsSupportRing", Ring_innerR, Ring_outerR, Ring_half_height_GAr, Ring_start_phi, Ring_delta_phi); }
  if (Ring_half_height_LAr > 0) { RingSolid_LAr = new G4Tubs("wlsSupportRing", Ring_innerR, Ring_outerR, Ring_half_height_LAr, Ring_start_phi, Ring_delta_phi); }

  // linear sector
  G4double Linsec_halfx = ArDMvar.WLS_SUPPORT_LINSEC_HALF_X;
  G4double Linsec_halfy = ArDMvar.WLS_SUPPORT_LINSEC_HALF_Y;
  G4double Linsec_halfz_GAr = ArDMvar.WLS_SUPPORT_LINSEC_HALF_HEIGHT_GAR;
  G4double Linsec_halfz_LAr = ArDMvar.WLS_SUPPORT_LINSEC_HALF_HEIGHT_LAR;

  G4Box* LinSecSolid_GAr = NULL;
  G4Box* LinSecSolid_LAr = NULL;

  if (Linsec_halfz_GAr > 0) { LinSecSolid_GAr = new G4Box("wlsSupportLinSec", Linsec_halfx, Linsec_halfy, Linsec_halfz_GAr); }
  if (Linsec_halfz_LAr > 0) { LinSecSolid_LAr = new G4Box("wlsSupportLinSec", Linsec_halfx, Linsec_halfy, Linsec_halfz_LAr); }

  G4RotationMatrix* Linsec_rotMat = new G4RotationMatrix;
  Linsec_rotMat->rotateX(0);
  Linsec_rotMat->rotateY(0);
  Linsec_rotMat->rotateZ(90. * deg);

  G4ThreeVector pos(ArDMvar.WLS_SUPPORT_LINSEC_POS_X, ArDMvar.WLS_SUPPORT_LINSEC_POS_Y,
                    0.);  // relative to the ring sector

  // build unionSolid of ring and linear sector
  G4UnionSolid* wls_GAr = NULL;
  G4UnionSolid* wls_LAr = NULL;
  if (RingSolid_GAr && LinSecSolid_GAr) { wls_GAr = new G4UnionSolid("WLSSupport", RingSolid_GAr, LinSecSolid_GAr, Linsec_rotMat, pos); }
  if (RingSolid_LAr && LinSecSolid_LAr) { wls_LAr = new G4UnionSolid("WLSSupport", RingSolid_LAr, LinSecSolid_LAr, Linsec_rotMat, pos); }
  DSStorage::Get()->SetArDMTeflonHeight(Ring_half_height_LAr + Ring_half_height_GAr);
  DSStorage::Get()->SetArDMTeflonRadius(Ring_outerR);
  G4LogicalVolume* wlsLog_GAr = NULL;
  G4LogicalVolume* wlsLog_LAr = NULL;

  fWLSSupportMat = fTeflon;
  if (wls_GAr) { wlsLog_GAr = new G4LogicalVolume(wls_GAr, fWLSSupportMat, "WLSSupport", 0, 0, 0); }
  if (wls_LAr) { wlsLog_LAr = new G4LogicalVolume(wls_LAr, fWLSSupportMat, "WLSSupport", 0, 0, 0); }

  // place fWLSLog in fLArColLog
  G4ThreeVector wlsPos_GAr(0, 0, ArDMvar.WLS_SUPPORT_RINGSEC_GARCOL_POS_Z_GARCOL);
  G4ThreeVector wlsPos_LAr(0, 0, ArDMvar.WLS_SUPPORT_RINGSEC_LARCOL_POS_Z_LARCOL);

  if (wlsLog_GAr && fGArColPhys) { fWLSSupportGArPhys = new G4PVPlacement(0, wlsPos_GAr, "WLSSupport", wlsLog_GAr, fGArColPhys, false, 0); }
  if (wlsLog_LAr && fLArColPhys) { fWLSSupportLArPhys = new G4PVPlacement(0, wlsPos_LAr, "WLSSupport", wlsLog_LAr, fLArColPhys, false, 0); }

  G4VisAttributes* wlsAtt = new G4VisAttributes(true);
  // G4VisAttributes* wlsAtt = new G4VisAttributes(false);
  wlsAtt->SetColour(1.0, .0, 1.0);  // magenta
  wlsAtt->SetForceAuxEdgeVisible(true);
  if (wlsLog_GAr) { wlsLog_GAr->SetVisAttributes(wlsAtt); }
  if (wlsLog_LAr) { wlsLog_LAr->SetVisAttributes(wlsAtt); }

  return;
}

void DSDetectorArDM::addPMT(G4String bottom_or_top) {
  ostringstream name;
  name << "pmt";

  fPMTMat = DSMaterial::Get()->GetFakePhotocathode();
  G4LogicalVolume* pmtLog = constructPMT(fPMTMat, bottom_or_top, 1, name.str().c_str());
  G4VPhysicalVolume* motherPhys;
  G4double z;

  if (bottom_or_top == "top") {
    if (fGArColPhys) {
      motherPhys = fGArColPhys;
      z = ArDMvar.TOP_PMT_Z_GARCOL;
    } else {
      motherPhys = fMotherVolume;
      z = ArDMvar.TOP_PMT_Z;
    }
  } else {
    if (fLArColPhys) {
      motherPhys = fLArColPhys;
      z = ArDMvar.BTM_PMT_Z_LARCOL;
    } else {
      motherPhys = fMotherVolume;
      z = ArDMvar.BTM_PMT_Z;
    }
  }

  std::vector<G4ThreeVector> rPMT = ArDMvar.setPMTVector(bottom_or_top);
  for (unsigned int i = 0; i < rPMT.size(); i++) { rPMT[i].setZ(z); }

  std::vector<G4VPhysicalVolume*> pmtPhysArray;
  pmtPhysArray = placePMT(pmtLog, motherPhys, rPMT, bottom_or_top);

  if (bottom_or_top == "top") {
    fTopPMTArrayPhys = pmtPhysArray;
  } else {
    fBtmPMTArrayPhys = pmtPhysArray;
  }

  addPMTCoat(bottom_or_top);
  addPMTCathode(bottom_or_top);  //<-- this is the detector !

  // addPMTElectrode(bottom_or_top);
  // addPMTBase(bottom_or_top);

  return;
}

G4LogicalVolume* DSDetectorArDM::constructPMT(G4Material* PMTMat, const char* bottom_or_top, int attribute, const char* name) {
  fPMTMat = PMTMat;
  if (strcmp(bottom_or_top, "top") && strcmp(bottom_or_top, "bottom")) {
    cout << "in Geom_ROGeom.cc, constructPMT(): argument is run. exit." << endl;
    return NULL;
  }

  // middle cylinder
  G4double middle_cylinder_innerR = ArDMvar.PMT_MIDDLE_CYLINDER_INNER_RADIUS;
  G4double middle_cylinder_outerR = ArDMvar.PMT_MIDDLE_CYLINDER_OUTER_RADIUS;
  G4double middle_cylinder_halfz = ArDMvar.PMT_MIDDLE_CYLINDER_HALF_HEIGHT;
  G4Tubs* middle_cylinder_solid = new G4Tubs("pmtMiddleCylinder", middle_cylinder_innerR, middle_cylinder_outerR, middle_cylinder_halfz, 0 * deg, 360 * deg);

  // spherical part : top
  G4double spherical_part_innerR_top = ArDMvar.PMT_SPHERICAL_PART_INNER_RADIUS;
  G4double spherical_part_outerR_top = ArDMvar.PMT_SPHERICAL_PART_OUTER_RADIUS;
  G4double spherical_part_start_theta_top = 0;
  G4double spherical_part_delta_theta_top = ArDMvar.PMT_SPHERICAL_PART_OPENING_ANGLE / 2;

  G4Sphere* sph_part_solid_top = new G4Sphere("sphPartTop", spherical_part_innerR_top, spherical_part_outerR_top, 0 * deg, 360 * deg, spherical_part_start_theta_top, spherical_part_delta_theta_top);

  // spherical part : btm
  G4double spherical_part_innerR_btm = ArDMvar.PMT_SPHERICAL_PART_INNER_RADIUS;
  G4double spherical_part_outerR_btm = ArDMvar.PMT_SPHERICAL_PART_OUTER_RADIUS;
  G4double spherical_part_start_theta_btm = 180 * deg - ArDMvar.PMT_SPHERICAL_PART_OPENING_ANGLE / 2;
  G4double spherical_part_delta_theta_btm = (ArDMvar.PMT_SPHERICAL_PART_OPENING_ANGLE - ArDMvar.PMT_BTM_SPHERE_HOLE_OPENING_ANGLE) / 2;

  G4Sphere* sph_part_solid_btm = new G4Sphere("sphPartBtm", spherical_part_innerR_btm, spherical_part_outerR_btm, 0 * deg, 360 * deg, spherical_part_start_theta_btm, spherical_part_delta_theta_btm);

  // btm tube
  // connecting the pmt with pmt base
  G4double btmTube_innerR = ArDMvar.PMT_BTM_TUBE_INNER_RADIUS;
  G4double btmTube_outerR = ArDMvar.PMT_BTM_TUBE_OUTER_RADIUS;
  G4double btmTube_halfz = ArDMvar.PMT_BTM_TUBE_HALF_HEIGHT;
  G4Tubs* btmTube_solid = new G4Tubs("pmtBtmTube", btmTube_innerR, btmTube_outerR, btmTube_halfz, 0 * deg, 360 * deg);

  // relative positions of parts of pmt to pmt center
  // pmt center = center of middle cylinder
  G4ThreeVector topSphPos(0, 0, -ArDMvar.PMT_DISTANCE_PMT_CENTER_TO_TOP_SPHERE_CENTER);
  G4ThreeVector btmSphPos(0, 0, ArDMvar.PMT_DISTANCE_PMT_CENTER_TO_BTM_SPHERE_CENTER);
  G4ThreeVector btmTubePos(0, 0, -ArDMvar.PMT_DISTANCE_PMT_CENTER_TO_BTM_TUBE_CENTER);

  G4UnionSolid* pmtSolid = new G4UnionSolid("pmt", middle_cylinder_solid, sph_part_solid_top, 0, topSphPos);
  pmtSolid = new G4UnionSolid("pmt", pmtSolid, sph_part_solid_btm, 0, btmSphPos);
  pmtSolid = new G4UnionSolid("pmt", pmtSolid, btmTube_solid, 0, btmTubePos);

  G4LogicalVolume* pmtLog = new G4LogicalVolume(pmtSolid, fPMTMat, name);

  G4VisAttributes* pmtAtt = new G4VisAttributes(attribute);
  // G4VisAttributes* pmtAtt = new G4VisAttributes(false);
  pmtAtt->SetColour(0.0, 1.0, 0.0);  // green
  pmtAtt->SetForceAuxEdgeVisible(true);
  pmtLog->SetVisAttributes(pmtAtt);
  return pmtLog;
}

std::vector<G4VPhysicalVolume*> DSDetectorArDM::placePMT(G4LogicalVolume* fPMTLog, G4VPhysicalVolume* fMotherPhys, std::vector<G4ThreeVector> rPMT, const char* bottom_or_top) {
  G4RotationMatrix* rot = NULL;

  if (!strcmp(bottom_or_top, "top")) {
    rot = new G4RotationMatrix;
    rot->rotateX(180. * deg);
    rot->rotateY(0. * deg);
    rot->rotateZ(0. * deg);
  }

  std::vector<G4VPhysicalVolume*> physVolVec;
  ostringstream pmtname;

  for (unsigned int i = 0; i < rPMT.size(); i++) {
    pmtname.str("");
    if (!strcmp(bottom_or_top, "top")) {
      pmtname << fPMTLog->GetName() << i;
    } else if (!strcmp(bottom_or_top, "bottom")) {
      pmtname << fPMTLog->GetName() << i + 12;
    } else {
      pmtname << fPMTLog->GetName();
    }
    physVolVec.push_back((G4VPhysicalVolume*)(new G4PVPlacement(rot, rPMT[i], pmtname.str().c_str(), fPMTLog, fMotherPhys, false, 0)));
  }

  return physVolVec;
}
std::vector<G4VPhysicalVolume*> placePMT(std::vector<G4LogicalVolume*> fPMTLog, G4VPhysicalVolume* fMotherPhys, std::vector<G4ThreeVector> rPMT, const char* bottom_or_top) {
  G4RotationMatrix* rot = NULL;

  if (!strcmp(bottom_or_top, "top")) {
    rot = new G4RotationMatrix;
    rot->rotateX(180. * deg);
    rot->rotateY(0. * deg);
    rot->rotateZ(0. * deg);
  }

  std::vector<G4VPhysicalVolume*> physVolVec;
  ostringstream pmtname;

  for (unsigned int i = 0; i < rPMT.size(); i++) {
    pmtname.str("");

    if (!strcmp(bottom_or_top, "top")) {
      pmtname << fPMTLog[i]->GetName() << i;
    } else if (!strcmp(bottom_or_top, "bottom")) {
      pmtname << fPMTLog[i]->GetName() << i + 12;
    } else {
      pmtname << fPMTLog[i]->GetName();
    }

    physVolVec.push_back((G4VPhysicalVolume*)(new G4PVPlacement(rot, rPMT[i], pmtname.str().c_str(), fPMTLog[i], fMotherPhys, false, 0)));
  }
  return physVolVec;
}

G4LogicalVolume* DSDetectorArDM::constructPMTAux(G4Material* PMTMat, const char* bottom_or_top, int attribute, G4double pmt_inner_radius, G4double pmt_outer_radius, const char* name, G4double deltaTheta) {
  fPMTMat = PMTMat;
  if (strcmp(bottom_or_top, "top") && strcmp(bottom_or_top, "bottom")) {
    cout << "in DSDetectorArDM, constructPMTAux(): argument is wrong. exit." << endl;
    return NULL;
  }

  G4double startTheta = 0;
  G4double dTheta = deltaTheta;

  G4Sphere* pmtSolid = new G4Sphere("pmtAux", pmt_inner_radius, pmt_outer_radius, 0. * deg, 360. * deg, startTheta, dTheta);
  G4LogicalVolume* pmtLog = new G4LogicalVolume(pmtSolid, fPMTMat, name);

  G4VisAttributes* pmtAtt = new G4VisAttributes(attribute);
  // G4VisAttributes* pmtAtt = new G4VisAttributes(false);
  pmtAtt->SetColour(0.0, 1.0, 0.0);  // green
  pmtAtt->SetForceAuxEdgeVisible(true);
  pmtLog->SetVisAttributes(pmtAtt);
  return pmtLog;
}

void DSDetectorArDM::addPMTCoat(G4String bottom_or_top) {
  ostringstream name;
  name << "pmtCoat";

  G4double innerR;
  G4double outerR;
  G4double openAngle;

  if (bottom_or_top == "top") {
    innerR = ArDMvar.PMT_COATING_INNER_RADIUS;
    outerR = ArDMvar.PMT_COATING_OUTER_RADIUS;
    openAngle = ArDMvar.PMT_SPHERICAL_PART_OPENING_ANGLE;

    for (G4int i = 0; i < 12; i++) {
      // ((ArRecoPMT*)fAna->fRecoPMTs->At(i)) ->eConvEffPMT = PMTCONV;
    }
  } else {
    innerR = ArDMvar.PMT_COATING_INNER_RADIUS;
    outerR = ArDMvar.PMT_COATING_OUTER_RADIUS;
    openAngle = ArDMvar.PMT_SPHERICAL_PART_OPENING_ANGLE;
    for (G4int i = 12; i < 24; i++) {
      //((ArRecoPMT*)fAna->fRecoPMTs->At(i)) ->eConvEffPMT = PMTCONV;
    }
  }

  G4double deltaTheta = openAngle / 2;
  G4LogicalVolume* pmtLog = NULL;
  fWLSMatTop = DSMaterial::Get()->GetTPB();
  fWLSMat = DSMaterial::Get()->GetTPB();
  if (bottom_or_top == "top") {
    // to force no photon trapping in the top PMT WLS layer
    pmtLog = constructPMTAux(fWLSMatTop, bottom_or_top, 1, innerR, outerR, name.str().c_str(), deltaTheta);
  } else {
    pmtLog = constructPMTAux(fWLSMat, bottom_or_top, 1, innerR, outerR, name.str().c_str(), deltaTheta);
  }

  G4VPhysicalVolume* motherPhys;
  G4double z;
  if (bottom_or_top == "top") {
    if (fGArColPhys) {
      motherPhys = fGArColPhys;
      z = ArDMvar.TOP_PMT_COATING_Z_GARCOL;
    } else {
      motherPhys = fMotherVolume;
      z = ArDMvar.TOP_PMT_COATING_Z;
    }
  } else {
    if (fLArColPhys) {
      motherPhys = fLArColPhys;
      z = ArDMvar.BTM_PMT_COATING_Z_LARCOL;
    } else {
      motherPhys = fMotherVolume;
      z = ArDMvar.BTM_PMT_COATING_Z;
    }
  }

  vector<G4ThreeVector> rPMT = ArDMvar.setPMTVector(bottom_or_top);
  for (unsigned int i = 0; i < rPMT.size(); i++) { rPMT[i].setZ(z); }

  vector<G4VPhysicalVolume*> pmtPhysArray;
  pmtPhysArray = placePMT(pmtLog, motherPhys, rPMT, bottom_or_top);

  if (bottom_or_top == "top") {
    fTopPMTCoatArrayPhys = pmtPhysArray;
  } else {
    fBtmPMTCoatArrayPhys = pmtPhysArray;
  }

  return;
}

void DSDetectorArDM::addPMTCathode(G4String bottom_or_top) {
  ostringstream name;
  name << "pmtCathode";

  // define the shape of the PMT
  G4double innerR = ArDMvar.PMT_CATHODE_INNER_RADIUS;
  G4double outerR = ArDMvar.PMT_CATHODE_OUTER_RADIUS;
  G4double deltaTheta = ArDMvar.PMT_CATHODE_ACTIVE_RANGE;
  fPMTCathodeMat = fPMTMat;
  G4LogicalVolume* pmtCathodeLog = constructPMTAux(fPMTCathodeMat, bottom_or_top, 1, innerR, outerR, name.str().c_str(), deltaTheta);

  // define the behaviour of the sensitive region of the PMT
  //   G4VSensitiveDetector* sensPMT = getSD(bottom_or_top);
  //   G4SDManager*          SDman   = G4SDManager::GetSDMpointer();
  //   SDman->AddNewDetector(sensPMT);

  // mark it as sensitive, if not the whole PMT's surface is sensitive --> use
  // readout geometry to define the sensitive region !
  // pmtCathodeLog->SetSensitiveDetector(sensPMT);

  // place PMT in the detector
  G4VPhysicalVolume* motherPhys;
  G4double z;
  if (bottom_or_top == "top") {
    if (fGArColPhys) {
      motherPhys = fGArColPhys;
      z = ArDMvar.TOP_PMT_CATHODE_Z_GARCOL;
    } else {
      motherPhys = fMotherVolume;
      z = ArDMvar.TOP_PMT_CATHODE_Z;
    }
  } else {
    if (fLArColPhys) {
      motherPhys = fLArColPhys;
      z = ArDMvar.BTM_PMT_CATHODE_Z_LARCOL;
    } else {
      motherPhys = fMotherVolume;
      z = ArDMvar.BTM_PMT_CATHODE_Z;
    }
  }

  vector<G4ThreeVector> rPMT = ArDMvar.setPMTVector(bottom_or_top);
  for (unsigned int i = 0; i < rPMT.size(); i++) { rPMT[i].setZ(z); }

  vector<G4VPhysicalVolume*> pmtPhysArray;
  pmtPhysArray = placePMT(pmtCathodeLog, motherPhys, rPMT, bottom_or_top);

  if (bottom_or_top == "top") {
    fTopPMTCathodeArrayPhys = pmtPhysArray;
  } else {
    fBtmPMTCathodeArrayPhys = pmtPhysArray;
  }

  return;
}

void DSDetectorArDM::addPMTSupport(G4String bottom_or_top) {
  G4double innerR = 0;
  G4double outerR = ArDMvar.PMT_SUPPORT_RADIUS;
  G4double halfz = ArDMvar.PMT_SUPPORT_HALF_HEIGHT;

  G4Tubs* pmtSupportPlateSolid = new G4Tubs("pmtSupportPlate", innerR, outerR, halfz, 0 * deg, 360. * deg);

  // making holes for PMTs on the PMT support structure
  G4double pmtHole_innerR = 0;
  G4double pmtHole_outerR = ArDMvar.PMT_MIDDLE_CYLINDER_OUTER_RADIUS;
  G4double pmtHole_half_height = halfz + 0.1 * mm;  // 0.1*mm is here so that we actually see a hole in the
                                                    // plate when calling G4SubtractionSolid

  G4Tubs* pmtHoleSolid = new G4Tubs("pmtHole", pmtHole_innerR, pmtHole_outerR, pmtHole_half_height, 0 * deg, 360. * deg);

  vector<G4TwoVector> rpmt2D = ArDMvar.setPMTVector2D();

  G4BooleanSolid* pmtSupportSolid;
  pmtSupportSolid = new G4SubtractionSolid("pmtSupportSolid", pmtSupportPlateSolid, pmtHoleSolid, 0, G4ThreeVector(rpmt2D[0].x(), rpmt2D[0].y(), 0));

  for (unsigned int i = 1; i < rpmt2D.size(); i++) { pmtSupportSolid = new G4SubtractionSolid("pmtSupportSolid", pmtSupportSolid, pmtHoleSolid, 0, G4ThreeVector(rpmt2D[i].x(), rpmt2D[i].y(), 0)); }

  ostringstream pmtSupportName;
  pmtSupportName << bottom_or_top << "PMTSupport";
  fPMTSupportMat = fTeflon;
  G4LogicalVolume* fPMTSupportLog = new G4LogicalVolume(pmtSupportSolid, fPMTSupportMat, pmtSupportName.str().c_str());

  G4VisAttributes* pmtSupportAtt = new G4VisAttributes(true);
  // G4VisAttributes* pmtSupportAtt  = new G4VisAttributes(false);
  pmtSupportAtt->SetColour(G4Colour::Yellow());
  pmtSupportAtt->SetForceAuxEdgeVisible(true);
  fPMTSupportLog->SetVisAttributes(pmtSupportAtt);

  G4VPhysicalVolume* motherPhys;
  G4double z;
  if (bottom_or_top == "top") {
    if (fGArColPhys) {
      motherPhys = fGArColPhys;
      z = ArDMvar.TOP_PMT_SUPPORT_POS_Z_GARCOL;
    } else {
      motherPhys = fMotherVolume;
      z = ArDMvar.TOP_PMT_SUPPORT_POS_Z;
    }
  } else {
    if (fLArColPhys) {
      motherPhys = fLArColPhys;
      z = ArDMvar.BTM_PMT_SUPPORT_POS_Z_LARCOL;
    } else {
      motherPhys = fMotherVolume;
      z = ArDMvar.BTM_PMT_SUPPORT_POS_Z;
    }
  }

  G4ThreeVector pos(0, 0, z);
  G4VPhysicalVolume* fPMTSupportPhys = new G4PVPlacement(0, pos, pmtSupportName.str().c_str(), fPMTSupportLog, motherPhys, false, 0);

  if (bottom_or_top == "top") {
    fTopPMTSupportPhys = fPMTSupportPhys;
  } else {
    fBtmPMTSupportPhys = fPMTSupportPhys;
  }

  addPMTSupportCoat(bottom_or_top);

  return;
}

void DSDetectorArDM::addPMTSupportCoat(G4String bottom_or_top) {
  G4double innerR = 0;
  G4double ourterR = ArDMvar.PMT_SUPPORT_RADIUS;
  G4double halfz = ArDMvar.PMT_SUPPORT_COATING_HALF_HEIGHT;

  G4Tubs* pmtSupportCoatSolid = new G4Tubs("pmtSupportCoat", innerR, ourterR, halfz, 0 * deg, 360. * deg);

  // making holes for PMTs on the PMT support structure
  G4double pmtHole_innerR = 0;
  G4double pmtHole_outerR = ArDMvar.PMT_COATING_OUTER_RADIUS * sin(ArDMvar.PMT_ACTIVE_RANGE);
  G4double pmtHole_half_height = halfz + 0.01 * mm;
  // 0.01*mm is just there so that we actually see a hole in the plate when
  // calling G4SubtractionSolid

  G4Tubs* pmtHoleSolid = new G4Tubs("pmtHole", pmtHole_innerR, pmtHole_outerR, pmtHole_half_height, 0 * deg, 360. * deg);

  vector<G4TwoVector> rpmt2D = ArDMvar.setPMTVector2D();
  G4BooleanSolid* pmtSupportSolid;

  pmtSupportSolid = new G4SubtractionSolid("pmtSupportCoatSolid", pmtSupportCoatSolid, pmtHoleSolid, 0, G4ThreeVector(rpmt2D[0].x(), rpmt2D[0].y(), 0));

  for (unsigned int i = 1; i < rpmt2D.size(); i++) { pmtSupportSolid = new G4SubtractionSolid("pmtSupportCoatSolid", pmtSupportSolid, pmtHoleSolid, 0, G4ThreeVector(rpmt2D[i].x(), rpmt2D[i].y(), 0)); }

  ostringstream pmtSupportCoatName;
  pmtSupportCoatName << bottom_or_top << "PMTSupportCoat";
  G4LogicalVolume* fPMTSupportCoatLog = new G4LogicalVolume(pmtSupportSolid, fWLSMat, pmtSupportCoatName.str().c_str());

  G4VisAttributes* pmtSupportCoatAtt = new G4VisAttributes(true);
  // G4VisAttributes* pmtSupportCoatAtt  = new G4VisAttributes(false);
  pmtSupportCoatAtt->SetColour(G4Colour::Green());
  pmtSupportCoatAtt->SetForceAuxEdgeVisible(true);
  fPMTSupportCoatLog->SetVisAttributes(pmtSupportCoatAtt);

  G4VPhysicalVolume* motherPhys;
  G4double z;
  if (bottom_or_top == "top") {
    if (fGArColPhys) {
      motherPhys = fGArColPhys;
      z = ArDMvar.TOP_PMT_SUPPORT_COATING_POS_Z_GARCOL;
    } else {
      motherPhys = fMotherVolume;
      z = ArDMvar.TOP_PMT_SUPPORT_COATING_POS_Z;
    }
  } else {
    if (fLArColPhys) {
      motherPhys = fLArColPhys;
      z = ArDMvar.BTM_PMT_SUPPORT_COATING_POS_Z_LARCOL;
    } else {
      motherPhys = fMotherVolume;
      z = ArDMvar.BTM_PMT_SUPPORT_COATING_POS_Z;
    }
  }

  G4ThreeVector pos(0, 0, z);

  G4VPhysicalVolume* fPMTSupportCoatPhys = new G4PVPlacement(0, pos, pmtSupportCoatName.str().c_str(), fPMTSupportCoatLog, motherPhys, false, 0);

  if (bottom_or_top == "top") {
    fTopPMTSupportCoatPhys = fPMTSupportCoatPhys;
  } else {
    fBtmPMTSupportCoatPhys = fPMTSupportCoatPhys;
  }

  return;
}

void DSDetectorArDM::addFieldShaperPillars() {
  G4double innerR = ArDMvar.FIELD_SHAPER_PILLAR_INNER_RADIUS;
  G4double outerR = ArDMvar.FIELD_SHAPER_PILLAR_OUTER_RADIUS;
  G4double halfz = ArDMvar.FIELD_SHAPER_PILLAR_HALF_HEIGHT;

  G4Material* fPolyethylene = DSMaterial::Get()->GetHDPE();
  G4Tubs* solid = new G4Tubs("fieldShaperPillar", innerR, outerR, halfz, 0. * deg, 360. * deg);
  G4LogicalVolume* log = new G4LogicalVolume(solid, fPolyethylene, "fieldShaperPillar");

  G4VisAttributes* att = new G4VisAttributes(true);
  // G4VisAttributes* att   = new G4VisAttributes(false);
  att->SetColour(G4Colour::Magenta());
  att->SetForceAuxEdgeVisible(true);
  log->SetVisAttributes(att);

  G4double posz = fLArColPhys ? ArDMvar.FIELD_SHAPER_PILLAR_POS_Z_LAR : ArDMvar.FIELD_SHAPER_PILLAR_POS_Z;  // relative position to LAr column
  G4VPhysicalVolume* motherPhys = fLArColPhys ? fLArColPhys : fMotherVolume;

  // place the pillars outside of the field shaper rings, just touching rings.
  // there are 7 of them, placed 45 degrees from each other.
  // the linear sector of the reflector is perpendicular to the x-axis,
  // the 1st pillar is 45 degrees away from the x-axis.

  G4RotationMatrix* rot = NULL;
  G4double distance_to_the_det_center = ArDMvar.DISTANCE_FIELD_SHAPER_PILLAR_CENTER_TO_DETECTOR_CENTER;

  G4int count = 0;
  ostringstream name;

  for (G4double startphi = 45 * deg; startphi < 360 * deg; startphi += 45 * deg) {
    G4double x = distance_to_the_det_center * cos(startphi);
    G4double y = distance_to_the_det_center * sin(startphi);
    G4ThreeVector pos(x, y, posz);

    name.str("");
    name << log->GetName() << count;
    fFieldShaperPillarPhys.push_back(new G4PVPlacement(rot, pos, name.str().c_str(), log, motherPhys, false, 0));
    count++;
  }

  return;
}

void DSDetectorArDM::addFieldShaperRings() {
  // the field shaper ring has similar shape as the WLS.
  //--> 2 part : ring sector and linear sector.
  // ring sector = a torus
  // lin  sector = a tube cut at both ends
  // ring sector

  G4double ringsec_crossec_innerR = ArDMvar.FIELD_SHAPER_RING_RINGSEC_CROSSSECTION_INNER_RADIUS;
  G4double ringsec_crossec_outerR = ArDMvar.FIELD_SHAPER_RING_RINGSEC_CROSSSECTION_OUTER_RADIUS;
  G4double ringsec_torusR = ArDMvar.FIELD_SHAPER_RING_RINGSEC_TORUS_RADIUS;
  G4double ringsec_startphi = ArDMvar.FIELD_SHAPER_RING_RINGSEC_TORUS_START_PHI;
  G4double ringsec_delphi = ArDMvar.FIELD_SHAPER_RING_RINGSEC_TORUS_DELTA_PHI;

  G4Torus* ringsec_torus = new G4Torus("fieldShaperRingTorus", ringsec_crossec_innerR, ringsec_crossec_outerR, ringsec_torusR, ringsec_startphi, ringsec_delphi);

  // linear sector
  G4double linsec_innerR = ArDMvar.FIELD_SHAPER_LINSEC_INNER_RADIUS;
  G4double linsec_outerR = ArDMvar.FIELD_SHAPER_LINSEC_OUTER_RADIUS;
  G4double linsec_halfz = ArDMvar.FIELD_SHAPER_LINSEC_HALF_LENGTH;

  G4Tubs* tubeSolid = new G4Tubs("fieldShaperLinTube", linsec_innerR, linsec_outerR, linsec_halfz, 0. * deg, 360. * deg);

  // cut the tubes at both ends
  // define a trapezoid with the correct inclining angle for its sides,
  // then subtract the trapezoid from the tube, placing the trapezoid at the
  // correct positions. keep in mind that, at the moment, the tube is parallel
  // to the z-direction. it will be rotated later so that it'd be parallel to
  // y-axis and perpendicular to x-axis, and placed at z = 0.

  // the dimension of the trapezoid can be as big as you want,
  // it just needs to have the correct inclining angle on the side surface,
  // which will be used to cut the tube.

  // x,y,z have to be large enough to cut the whole xy cross section of the
  // tube,
  G4double trapezoid_z = 2 * ArDMvar.FIELD_SHAPER_LINSEC_OUTER_RADIUS + 10 * mm;  // 10*mm = dummy offset, just here to make
                                                                                  // sure that you cut enough off the tube
  G4double trapezoid_y_small = 2 * ArDMvar.FIELD_SHAPER_LINSEC_OUTER_RADIUS;
  G4double trapezoid_y_large = trapezoid_y_small + trapezoid_z * tan(90 * deg - ArDMvar.FIELD_SHAPER_RING_RINGSEC_TORUS_START_PHI);
  G4double trapezoid_x_small = 2 * ArDMvar.FIELD_SHAPER_LINSEC_OUTER_RADIUS;
  G4double trapezoid_x_large = trapezoid_x_small + trapezoid_z * tan(90 * deg - ArDMvar.FIELD_SHAPER_RING_RINGSEC_TORUS_START_PHI);

  G4Trd* trapezoidSolid = new G4Trd("trapezoid", trapezoid_x_small, trapezoid_x_large, trapezoid_y_small, trapezoid_y_large, trapezoid_z);

  /*************************************************************
    //subtracting the trapezoid from the tube
    //before subtracting, rotate the trapezoid by 90 degrees
    //clock-wise around the y-axis
    //
    // |\
    // |  \
    // |    \
    // |      \
    // |        \
    // |          \
    // |            \
    // |             |
    // |             |
    // |     A       |
    // |             |
    // |             |
    // |            /
    // |   __B___ /
    // |  |     /|
    // |  |  C/  |
    // |  | /    |
    // |  /      |
    // |/ |      |
    //    |      |
    //    |      |
    //    |  D   |
    //    |      |
    //    |      |
    //    |      |
    //    |      |
    //    |      |
    //    |______|
    //distance trapezoid's center <--> tube's center = AC + CD
    //AC = (small base + large base) / 2
    //CD = BD - CD, BD = half length of the tube
    //BC = tubeRadius * tan(ArDMvar.FIELD_SHAPER_RING_RINGSEC_TORUS_START_PHI)
    *************************************************************/
  G4RotationMatrix* trapezoidRotMat = new G4RotationMatrix;
  trapezoidRotMat->rotateY(ArDMvar.HANDEDNESS * 90 * deg);

  G4double AC = (trapezoid_x_small + trapezoid_x_large) / 2;
  G4double BD = linsec_halfz;
  G4double BC = linsec_outerR * tan(ArDMvar.FIELD_SHAPER_RING_RINGSEC_TORUS_START_PHI);
  G4double CD = BD - BC;
  G4double AB = AC - BC;
  G4double distance_trapezoidCenter_to_tubeCenter_upperEnd = AC + BD - CD;
  G4double distance_trapezoidCenter_to_tubeCenter_lowerEnd = distance_trapezoidCenter_to_tubeCenter_upperEnd - (2 * BD + AB);

  G4ThreeVector posUpperEnd(0, 0, distance_trapezoidCenter_to_tubeCenter_upperEnd);
  G4ThreeVector posLowerEnd(0, 0, distance_trapezoidCenter_to_tubeCenter_lowerEnd);

  G4SubtractionSolid* tubecut = new G4SubtractionSolid("fieldShaperRingLinTubeCut", tubeSolid, trapezoidSolid, trapezoidRotMat, posUpperEnd);
  tubecut = new G4SubtractionSolid("fieldShaperRingLinTubeCut", tubeSolid, trapezoidSolid, trapezoidRotMat, posLowerEnd);

  G4RotationMatrix* tubeRotMat = new G4RotationMatrix;
  tubeRotMat->rotateX(ArDMvar.HANDEDNESS * (-90) * deg);
  G4ThreeVector tubecutpos(ringsec_torusR * cos(ringsec_startphi), 0, 0);

  G4UnionSolid* solid = new G4UnionSolid("fieldShaperRing", ringsec_torus, tubecut, tubeRotMat, tubecutpos);
  G4LogicalVolume* log = new G4LogicalVolume(solid, fTankMat, "fieldShaperRing");

  G4VisAttributes* att = new G4VisAttributes(true);
  // G4VisAttributes* att   = new G4VisAttributes(false);
  att->SetColour(G4Colour::Blue());
  att->SetForceAuxEdgeVisible(true);
  log->SetVisAttributes(att);

  G4VPhysicalVolume* motherPhys = NULL;

  // place the pillars outside of the field shaper rings, just touching rings.
  // there are 7 of them, placed 45 degrees from each other.
  // the linear sector of the reflector is perpendicular to the x-axis,
  // the 1st pillar is 45 degrees away from the x-axis.
  G4int count = 0;
  ostringstream name;

  for (G4double ringz = ArDMvar.FIRST_FIELD_SHAPER_RING_POS_Z; ringz > ArDMvar.FIELD_SHAPER_PILLAR_POS_Z - ArDMvar.FIELD_SHAPER_PILLAR_HALF_HEIGHT; ringz -= ArDMvar.FIELD_SHAPER_RING_PITCH) {
    G4double arvol_z = 0;

    if (ringz > ArDMvar.LIQUID_SURFACE_POS_Z) {
      motherPhys = fGArColPhys ? fGArColPhys : fMotherVolume;
      arvol_z = fGArColPhys ? ArDMvar.GAR_COLUMN_Z : 0;
    } else {
      motherPhys = fLArColPhys ? fLArColPhys : fMotherVolume;
      arvol_z = fLArColPhys ? ArDMvar.LARCOL_CYLINDER_POS_Z : 0;
    }

    G4ThreeVector pos(0, 0, ringz - arvol_z);

    name.str("");
    name << log->GetName() << count;
    fFieldShaperRingsPhys.push_back(new G4PVPlacement(NULL, pos, name.str().c_str(), log, motherPhys, false, 0));

    count++;
  }

  return;
}

void DSDetectorArDM::addBtmSideReflector() {
  // bottom side reflector
  // a cone spanned by the lower edge of the reflector and the bottom pmtSupport
  // during gas test in mar/apr 2013
  G4double innerR1 = ArDMvar.BTM_SIDE_REFLECTOR_INNER_RADIUS_LOWER_EDGE;
  G4double outerR1 = ArDMvar.BTM_SIDE_REFLECTOR_OUTER_RADIUS_LOWER_EDGE;
  G4double innerR2 = ArDMvar.BTM_SIDE_REFLECTOR_INNER_RADIUS_UPPER_EDGE;
  G4double outerR2 = ArDMvar.BTM_SIDE_REFLECTOR_OUTER_RADIUS_UPPER_EDGE;
  G4double half_height = ArDMvar.BTM_SIDE_REFLECTOR_HALF_HEIGHT_GASTEST;
  G4double start_phi = ArDMvar.BTM_SIDE_REFLECTOR_START_PHI;
  G4double delta_phi = ArDMvar.BTM_SIDE_REFLECTOR_DELTA_PHI;

  G4Cons* coneSolid = new G4Cons("bottomSideReflector", innerR1, outerR1, innerR2, outerR2, half_height, start_phi, delta_phi);

  G4double cynPart_innerR = ArDMvar.BTM_SIDE_REFLECTOR_CYNLIDRICAL_PART_INNER_RADIUS;
  G4double cynPart_outerR = ArDMvar.BTM_SIDE_REFLECTOR_CYNLIDRICAL_PART_OUTER_RADIUS;
  G4double cynPart_half_height = ArDMvar.BTM_SIDE_REFLECTOR_CYNLIDRICAL_PART_HALF_HEIGHT;

  G4Tubs* tubeSolid = new G4Tubs("bottomSideReflector", cynPart_innerR, cynPart_outerR, cynPart_half_height, 0 * deg, 360. * deg);

  G4ThreeVector pos_tube_rel_cone(0., 0., ArDMvar.BTM_SIDE_REFLECTOR_CYNLIDRICAL_PART_POS_Z - ArDMvar.BTM_SIDE_REFLECTOR_POS_Z_GASTEST);

  G4UnionSolid* solid = new G4UnionSolid("bottomSideReflector", coneSolid, tubeSolid, 0, pos_tube_rel_cone);

  G4Material* fBtmSideReflectorMat;
  fBtmSideReflectorMat = fTeflon;
  fBtmSideReflectorLog = new G4LogicalVolume(solid, fBtmSideReflectorMat, "bottomSideReflector", 0, 0, 0);

  // place fBtmSideReflectorLog in fLArColLog
  G4ThreeVector pos(0, 0, ArDMvar.BTM_SIDE_REFLECTOR_POS_Z_GASTEST_LARCOL);

  // G4VPhysicalVolume* fBtmSideReflectorPhys =
  //     new G4PVPlacement(0, pos, "bottomSideReflector", fBtmSideReflectorLog,
  //     fLArColPhys, false, 0);

  G4VisAttributes* sideReflAtt = new G4VisAttributes(true);
  // G4VisAttributes* sideReflAtt = new G4VisAttributes(false);
  sideReflAtt->SetColour(G4Color::Green());
  sideReflAtt->SetForceAuxEdgeVisible(true);
  fBtmSideReflectorLog->SetVisAttributes(sideReflAtt);

  //#if ArDMvar.TURN_ON_SIDE_REFLECTOR_COATING
  addBtmSideReflectorCoat();
  //#endif //ArDMvar.TURN_ON_SIDE_REFLECTOR_COATING

  return;
}

void DSDetectorArDM::addTopSideReflector() {
  // top side refl. may span over GArCol-volume and LArCol-volume
  // hence, placing the top side refl. in only GArCol or LArCol will cause weird
  // things namely : assuming we place the top side refl. as daughter volume of
  // GArCol, then :
  //
  // 1. GEANT4 doesn't complain when the top side reflector sticks outside of
  // the GArCol-volume
  //    --> keyword "overlapping"
  //
  // 2. but when we build the optical surface between GArCol / LArCol and the
  // top side reflector, only the surface between GArCol and the top side
  // reflector works as we want it to. the surface between LArCol and the top
  // side reflector doesn't !
  //
  // 3. moreover, the optical photon can somehow just simply go through the part
  // of the top side reflector, which sticks into the LArCol volume
  //
  //--> solution :
  // devide the topside reflector into 2 parts
  // i.  the part in GArCol, and
  // ii. the part in LArCol

  G4double innerR1_GAr = ArDMvar.TOP_SIDE_REFLECTOR_INNER_RADIUS_UPPER_EDGE_IN_GARCOL;
  G4double outerR1_GAr = ArDMvar.TOP_SIDE_REFLECTOR_OUTER_RADIUS_UPPER_EDGE_IN_GARCOL;
  G4double innerR2_GAr = ArDMvar.TOP_SIDE_REFLECTOR_INNER_RADIUS_LOWER_EDGE_IN_GARCOL;
  G4double outerR2_GAr = ArDMvar.TOP_SIDE_REFLECTOR_OUTER_RADIUS_LOWER_EDGE_IN_GARCOL;
  G4double halfheight_GAr = ArDMvar.TOP_SIDE_REFLECTOR_HALF_HEIGHT_IN_GARCOL;

  G4double innerR1_LAr = ArDMvar.TOP_SIDE_REFLECTOR_INNER_RADIUS_UPPER_EDGE_IN_LARCOL;
  G4double outerR1_LAr = ArDMvar.TOP_SIDE_REFLECTOR_OUTER_RADIUS_UPPER_EDGE_IN_LARCOL;
  G4double innerR2_LAr = ArDMvar.TOP_SIDE_REFLECTOR_INNER_RADIUS_LOWER_EDGE_IN_LARCOL;
  G4double outerR2_LAr = ArDMvar.TOP_SIDE_REFLECTOR_OUTER_RADIUS_LOWER_EDGE_IN_LARCOL;
  G4double halfheight_LAr = ArDMvar.TOP_SIDE_REFLECTOR_HALF_HEIGHT_IN_LARCOL;

  G4double start_phi = ArDMvar.TOP_SIDE_REFLECTOR_START_PHI;
  G4double delta_phi = ArDMvar.TOP_SIDE_REFLECTOR_DELTA_PHI;

  G4Cons* ringSolid_GAr = NULL;
  G4Cons* ringSolid_LAr = NULL;

  if (halfheight_GAr > 0) { ringSolid_GAr = new G4Cons("topSideReflector", innerR2_GAr, outerR2_GAr, innerR1_GAr, outerR1_GAr, halfheight_GAr, start_phi, delta_phi); }

  if (halfheight_LAr > 0) { ringSolid_LAr = new G4Cons("topSideReflector", innerR2_LAr, outerR2_LAr, innerR1_LAr, outerR1_LAr, halfheight_LAr, start_phi, delta_phi); }

  G4LogicalVolume* fTopSideReflectorLog_GAr = NULL;
  G4LogicalVolume* fTopSideReflectorLog_LAr = NULL;

  G4Material* fTopSideReflectorMat;
  fTopSideReflectorMat = fTeflon;
  if (ringSolid_GAr) { fTopSideReflectorLog_GAr = new G4LogicalVolume(ringSolid_GAr, fTopSideReflectorMat, "topSideReflector", 0, 0, 0); }

  if (ringSolid_LAr) { fTopSideReflectorLog_LAr = new G4LogicalVolume(ringSolid_LAr, fTopSideReflectorMat, "topSideReflector", 0, 0, 0); }

  G4ThreeVector pos_GAr(0, 0, ArDMvar.TOP_SIDE_REFLECTOR_IN_GARCOL_POS_Z_GARCOL);
  G4ThreeVector pos_LAr(0, 0, ArDMvar.TOP_SIDE_REFLECTOR_IN_LARCOL_POS_Z_LARCOL);

  fTopSideReflectorGArPhys = NULL;
  fTopSideReflectorLArPhys = NULL;

  if (fTopSideReflectorLog_GAr && fGArColPhys) { fTopSideReflectorGArPhys = new G4PVPlacement(0, pos_GAr, "topSideReflector", fTopSideReflectorLog_GAr, fGArColPhys, false, 0); }

  if (fTopSideReflectorLog_LAr && fLArColPhys) { fTopSideReflectorLArPhys = new G4PVPlacement(0, pos_LAr, "topSideReflector", fTopSideReflectorLog_LAr, fLArColPhys, false, 0); }

  G4VisAttributes* sideReflAtt = new G4VisAttributes(true);
  // G4VisAttributes* sideReflAtt = new G4VisAttributes(false);
  sideReflAtt->SetColour(G4Color::Green());
  sideReflAtt->SetForceAuxEdgeVisible(true);
  if (fTopSideReflectorLog_GAr) fTopSideReflectorLog_GAr->SetVisAttributes(sideReflAtt);
  if (fTopSideReflectorLog_LAr) fTopSideReflectorLog_LAr->SetVisAttributes(sideReflAtt);

  //#if ArDMvar.TURN_ON_SIDE_REFLECTOR_COATING
  addTopSideReflectorCoat();
  //#endif //ArDMvar.TURN_ON_SIDE_REFLECTOR_COATING

  return;
}

void DSDetectorArDM::addBtmSideReflectorCoat() {
  G4double innerR1 = ArDMvar.BTM_SIDE_REFLECTOR_COAT_INNER_RADIUS_LOWER_EDGE;
  G4double outerR1 = ArDMvar.BTM_SIDE_REFLECTOR_COAT_OUTER_RADIUS_LOWER_EDGE;
  G4double innerR2 = ArDMvar.BTM_SIDE_REFLECTOR_COAT_INNER_RADIUS_UPPER_EDGE;
  G4double outerR2 = ArDMvar.BTM_SIDE_REFLECTOR_COAT_OUTER_RADIUS_UPPER_EDGE;
  G4double half_height = ArDMvar.BTM_SIDE_REFLECTOR_COAT_HALF_HEIGHT_GASTEST;
  G4double start_phi = ArDMvar.BTM_SIDE_REFLECTOR_START_PHI;
  G4double delta_phi = ArDMvar.BTM_SIDE_REFLECTOR_DELTA_PHI;

  G4Cons* coneSolid = new G4Cons("bottomSideReflector", innerR1, outerR1, innerR2, outerR2, half_height, start_phi, delta_phi);

  G4double cynPart_innerR = ArDMvar.BTM_SIDE_REFLECTOR_COAT_CYNLIDRICAL_PART_INNER_RADIUS;
  G4double cynPart_outerR = ArDMvar.BTM_SIDE_REFLECTOR_COAT_CYNLIDRICAL_PART_OUTER_RADIUS;
  G4double cynPart_half_height = ArDMvar.BTM_SIDE_REFLECTOR_COAT_CYNLIDRICAL_PART_HALF_HEIGHT;

  G4Tubs* tubeSolid = new G4Tubs("bottomSideReflector", cynPart_innerR, cynPart_outerR, cynPart_half_height, 0 * deg, 360. * deg);

  G4ThreeVector pos_tube_rel_cone(0, 0, ArDMvar.BTM_SIDE_REFLECTOR_COAT_CYNLIDRICAL_PART_POS_Z - ArDMvar.BTM_SIDE_REFLECTOR_COAT_POS_Z_GASTEST);

  G4UnionSolid* solid = new G4UnionSolid("bottomSideReflector", coneSolid, tubeSolid, 0, pos_tube_rel_cone);

  fBtmSideReflectorLog = new G4LogicalVolume(solid, fWLSMat, "bottomSideReflectorCoat", 0, 0, 0);

  // place fBtmSideReflectorLog in fLArColLog
  G4ThreeVector pos(0, 0, ArDMvar.BTM_SIDE_REFLECTOR_COAT_POS_Z_GASTEST_LARCOL);

  fBtmSideReflectorCoatPhys = new G4PVPlacement(0, pos, "bottomSideReflectorCoat", fBtmSideReflectorLog, fLArColPhys, false, 0);

  G4VisAttributes* sideReflAtt = new G4VisAttributes(true);
  // G4VisAttributes* sideReflAtt = new G4VisAttributes(false);
  sideReflAtt->SetColour(G4Color::Magenta());
  sideReflAtt->SetForceAuxEdgeVisible(true);
  fBtmSideReflectorLog->SetVisAttributes(sideReflAtt);

  return;
}
void DSDetectorArDM::addTopSideReflectorCoat() {
  G4double innerR1_GAr = ArDMvar.TOP_SIDE_REFLECTOR_COAT_INNER_RADIUS_UPPER_EDGE_IN_GARCOL;
  G4double outerR1_GAr = ArDMvar.TOP_SIDE_REFLECTOR_COAT_OUTER_RADIUS_UPPER_EDGE_IN_GARCOL;
  G4double innerR2_GAr = ArDMvar.TOP_SIDE_REFLECTOR_COAT_INNER_RADIUS_LOWER_EDGE_IN_GARCOL;
  G4double outerR2_GAr = ArDMvar.TOP_SIDE_REFLECTOR_COAT_OUTER_RADIUS_LOWER_EDGE_IN_GARCOL;
  G4double halfheight_GAr = ArDMvar.TOP_SIDE_REFLECTOR_COAT_HALF_HEIGHT_IN_GARCOL;

  G4double innerR1_LAr = ArDMvar.TOP_SIDE_REFLECTOR_COAT_INNER_RADIUS_UPPER_EDGE_IN_LARCOL;
  G4double outerR1_LAr = ArDMvar.TOP_SIDE_REFLECTOR_COAT_OUTER_RADIUS_UPPER_EDGE_IN_LARCOL;
  G4double innerR2_LAr = ArDMvar.TOP_SIDE_REFLECTOR_COAT_INNER_RADIUS_LOWER_EDGE_IN_LARCOL;
  G4double outerR2_LAr = ArDMvar.TOP_SIDE_REFLECTOR_COAT_OUTER_RADIUS_LOWER_EDGE_IN_LARCOL;
  G4double halfheight_LAr = ArDMvar.TOP_SIDE_REFLECTOR_COAT_HALF_HEIGHT_IN_LARCOL;

  G4double start_phi = ArDMvar.TOP_SIDE_REFLECTOR_COAT_START_PHI;
  G4double delta_phi = ArDMvar.TOP_SIDE_REFLECTOR_COAT_DELTA_PHI;

  G4Cons* ringSolid_GAr = NULL;
  G4Cons* ringSolid_LAr = NULL;

  if (halfheight_GAr) { ringSolid_GAr = new G4Cons("topSideReflectorCoat", innerR2_GAr, outerR2_GAr, innerR1_GAr, outerR1_GAr, halfheight_GAr, start_phi, delta_phi); }

  if (halfheight_LAr) { ringSolid_LAr = new G4Cons("topSideReflectorCoat", innerR2_LAr, outerR2_LAr, innerR1_LAr, outerR1_LAr, halfheight_LAr, start_phi, delta_phi); }

  G4LogicalVolume* fTopSideReflectorLog_GAr = NULL;
  G4LogicalVolume* fTopSideReflectorLog_LAr = NULL;

  if (ringSolid_GAr) { fTopSideReflectorLog_GAr = new G4LogicalVolume(ringSolid_GAr, fWLSMat, "topSideReflectorCoat", 0, 0, 0); }

  if (ringSolid_LAr) { fTopSideReflectorLog_LAr = new G4LogicalVolume(ringSolid_LAr, fWLSMat, "topSideReflectorCoat", 0, 0, 0); }

  G4ThreeVector pos_GAr(0, 0, ArDMvar.TOP_SIDE_REFLECTOR_COAT_IN_GARCOL_POS_Z_GARCOL);
  G4ThreeVector pos_LAr(0, 0, ArDMvar.TOP_SIDE_REFLECTOR_COAT_IN_LARCOL_POS_Z_LARCOL);

  if (fTopSideReflectorLog_GAr && fGArColPhys) { fTopSideReflectorCoatGArPhys = new G4PVPlacement(0, pos_GAr, "topSideReflectorCoat", fTopSideReflectorLog_GAr, fGArColPhys, false, 0); }

  if (fTopSideReflectorLog_LAr && fLArColPhys) { fTopSideReflectorCoatLArPhys = new G4PVPlacement(0, pos_LAr, "topSideReflectorCoat", fTopSideReflectorLog_LAr, fLArColPhys, false, 0); }

  G4VisAttributes* sideReflAtt = new G4VisAttributes(true);
  // G4VisAttributes* sideReflAtt = new G4VisAttributes(false);
  sideReflAtt->SetColour(G4Color::Magenta());
  sideReflAtt->SetForceAuxEdgeVisible(true);
  fTopSideReflectorLog_GAr->SetVisAttributes(sideReflAtt);
  fTopSideReflectorLog_LAr->SetVisAttributes(sideReflAtt);

  return;
}

DSDetectorArDM::~DSDetectorArDM() {
  ;
}
// took grom G4DSDetector
void DSDetectorArDM::DefineSurfaces() {

  /*
  ////////////////////////////////////////
  // GAr - LAr
  ////////////////////////////////////////
  G4OpticalSurface *fOpGArLArSurface = new G4OpticalSurface("OpGArLArSurface");
  new G4LogicalBorderSurface("GArLArSurface", fPhysicGasPocket, fPhysicInnerLAr,
fOpGArLArSurface); fOpGArLArSurface->SetType( dielectric_dielectric );
  fOpGArLArSurface->SetModel( unified );
  fOpGArLArSurface->SetFinish( polished );
  G4MaterialPropertiesTable *fGArLArSurfProp = new G4MaterialPropertiesTable();
  fOpGArLArSurface->SetMaterialPropertiesTable( fGArLArSurfProp );


  ////////////////////////////////////////
  // LAr - StainlessSteel (supportRod + compression plates)
  ////////////////////////////////////////
  G4OpticalSurface *fOpLArSteelSurface = new
G4OpticalSurface("OpLArSteelSurface"); new
G4LogicalBorderSurface("LArSteelSurfaceSupp0", fPhysicLArBath,
fPhysicSupportRod0, fOpLArSteelSurface); new
G4LogicalBorderSurface("LArSteelSurfaceSupp1", fPhysicLArBath,
fPhysicSupportRod1, fOpLArSteelSurface); new
G4LogicalBorderSurface("LArSteelSurfaceSupp2", fPhysicLArBath,
fPhysicSupportRod2, fOpLArSteelSurface); new
G4LogicalBorderSurface("LArSteelSurfaceSupp3", fPhysicLArBath,
fPhysicSupportRod3, fOpLArSteelSurface); new
G4LogicalBorderSurface("LArSteelSurfacePlateTop", fPhysicLArBath,
fPhysicCompressionPlateTop, fOpLArSteelSurface); new
G4LogicalBorderSurface("LArSteelSurfacePlateBot", fPhysicLArBath,
fPhysicCompressionPlateBottom, fOpLArSteelSurface); fOpLArSteelSurface->SetType(
dielectric_metal ); fOpLArSteelSurface->SetModel( unified );
  fOpLArSteelSurface->SetFinish( polished );
  fOpLArSteelSurface->SetSigmaAlpha(0.1);
  G4MaterialPropertiesTable *fLArSteelSurfProp = new
G4MaterialPropertiesTable(); G4double LArSteelENE[2] = {0.1*eV, 20.0*eV};
  G4double LArSteelREF[2] = {0.99, 0.99};
  G4double LArSteelEFF[2] = {0.00, 0.00};
  fLArSteelSurfProp->AddProperty("REFLECTIVITY", LArSteelENE, LArSteelREF, 2);
  fLArSteelSurfProp->AddProperty("EFFICIENCY",   LArSteelENE, LArSteelEFF, 2);
  fOpLArSteelSurface->SetMaterialPropertiesTable( fLArSteelSurfProp );


  ////////////////////////////////////////
  // LAr - GridSteel
  ////////////////////////////////////////
  G4OpticalSurface *fOpLArGridSurface = new
G4OpticalSurface("OpLArGridSurface"); new
G4LogicalBorderSurface("LArGridSurfacePlateBot", fPhysicInnerLAr,  fPhysicGrid,
fOpLArGridSurface); fOpLArGridSurface->SetType( dielectric_dielectric );
  fOpLArGridSurface->SetModel( unified );
  fOpLArGridSurface->SetFinish( polished );
  //fOpLArGridSurface->SetSigmaAlpha(0.1);
  G4MaterialPropertiesTable *fLArGridSurfProp = new G4MaterialPropertiesTable();
  //G4double LArGridENE[2] = {0.1*eV, 20.0*eV};
  //G4double LArGridREF[2] = {0.11, 0.11};   // 89% absorption
  //G4double LArGridEFF[2] = {0.00, 0.00};
  //fLArGridSurfProp->AddProperty("REFLECTIVITY", LArGridENE, LArGridREF, 2);
  //fLArGridSurfProp->AddProperty("EFFICIENCY",   LArGridENE, LArGridEFF, 2);
  fOpLArGridSurface->SetMaterialPropertiesTable( fLArGridSurfProp );


  ////////////////////////////////////////
  // LArBath - Copper (dewar)
  ////////////////////////////////////////
#if 0
  G4OpticalSurface *fOpLArCopperSurface = new
G4OpticalSurface("OpLArCopperSurface"); new
G4LogicalBorderSurface("LArCopperSurface", fPhysicLArBath, fPhysicDewar,
fOpLArCopperSurface); fOpLArCopperSurface->SetType( dielectric_metal );
  fOpLArCopperSurface->SetModel( unified );
  fOpLArCopperSurface->SetFinish( polished );
  fOpLArCopperSurface->SetSigmaAlpha(0.1);
  G4MaterialPropertiesTable *fLArCopperSurfProp = new
G4MaterialPropertiesTable(); G4double LArCopperENE[2] = {0.1*eV, 20.0*eV};
  G4double LArCopperREF[2] = {0.99, 0.99};
  G4double LArCopperEFF[2] = {0.00, 0.00};
  fLArCopperSurfProp->AddProperty("REFLECTIVITY", LArCopperENE, LArCopperREF,
2); fLArCopperSurfProp->AddProperty("EFFICIENCY",   LArCopperENE, LArCopperEFF,
2); fOpLArCopperSurface->SetMaterialPropertiesTable( fLArCopperSurfProp );
#endif

  ////////////////////////////////////////
  // TPB - GAr
  ////////////////////////////////////////
  G4OpticalSurface *fOpTPBGArSurface = new G4OpticalSurface("OpTPBGArSurface");
  new G4LogicalBorderSurface("TPBGArSurface", fPhysicGasPocket,
fPhysicTPBTopLateral, fOpTPBGArSurface ); new
G4LogicalBorderSurface("TPBGArSurface", fPhysicGasPocket, fPhysicTPBTopBase,
fOpTPBGArSurface ); fOpTPBGArSurface->SetType( dielectric_dielectric );
  fOpTPBGArSurface->SetModel( unified );
  fOpTPBGArSurface->SetFinish( ground );
  G4MaterialPropertiesTable *fTPBGArSurfProp = new G4MaterialPropertiesTable();
  fOpTPBGArSurface->SetMaterialPropertiesTable( fTPBGArSurfProp );


  ////////////////////////////////////////
  // LAr(bath) - Kapton
  ////////////////////////////////////////
  G4OpticalSurface *fOpLArKaptonSurface = new
G4OpticalSurface("OpLArKaptonSurface"); new
G4LogicalBorderSurface("LArKaptonSurface", fPhysicLArBath, fPhysicKaptonBand,
fOpLArKaptonSurface ); fOpLArKaptonSurface->SetType( dielectric_metal );
  fOpLArKaptonSurface->SetModel( unified );
  fOpLArKaptonSurface->SetFinish( polished );
  fOpLArKaptonSurface->SetSigmaAlpha(0.1);
  G4MaterialPropertiesTable *fLArKaptonSurfProp = new
G4MaterialPropertiesTable(); G4double LArKaptonENE[2] = {0.1*eV, 20.0*eV};
  G4double LArKaptonREF[2] = {0.99, 0.99};
  G4double LArKaptonEFF[2] = {0.00, 0.00};
  fLArKaptonSurfProp->AddProperty("REFLECTIVITY", LArKaptonENE, LArKaptonREF,
2); fLArKaptonSurfProp->AddProperty("EFFICIENCY",   LArKaptonENE, LArKaptonEFF,
2); fOpLArKaptonSurface->SetMaterialPropertiesTable( fLArKaptonSurfProp );


  ////////////////////////////////////////
  // TPB - Acrylic(inner vessel)
  ////////////////////////////////////////
  G4OpticalSurface *fOpTPBVesselSurface = new
G4OpticalSurface("OpTPBVesselSurface"); new
G4LogicalBorderSurface("TPBVesselSurfaceTop",    fPhysicTPBTopBase,
fPhysicInnerVesselWindowTop,    fOpTPBVesselSurface ); new
G4LogicalBorderSurface("TPBVesselSurfaceBot",    fPhysicTPBBottomBase,
fPhysicInnerVesselWindowBottom, fOpTPBVesselSurface ); new
G4LogicalBorderSurface("TPBVesselSurfaceTopLat", fPhysicTPBTopLateral,
fPhysicInnerVesselWall,         fOpTPBVesselSurface ); new
G4LogicalBorderSurface("TPBVesselSurfaceMidLat", fPhysicTPBMiddleLateral,
fPhysicInnerVesselWall,         fOpTPBVesselSurface ); new
G4LogicalBorderSurface("TPBVesselSurfaceLowLat", fPhysicTPBLowerLateral,
fPhysicInnerVesselWall,         fOpTPBVesselSurface );
  fOpTPBVesselSurface->SetType( dielectric_dielectric );
  fOpTPBVesselSurface->SetModel( unified );
  fOpTPBVesselSurface->SetFinish( ground );
  G4MaterialPropertiesTable *fTPBVesselSurfProp = new
G4MaterialPropertiesTable(); fOpTPBVesselSurface->SetMaterialPropertiesTable(
fTPBVesselSurfProp );


  ////////////////////////////////////////
  // LAr - Acrylic(inner vessel)
  ////////////////////////////////////////
  G4OpticalSurface *fOpLArVesselSurface = new
G4OpticalSurface("OpLArVesselSurface"); new
G4LogicalBorderSurface("LArVesselWall",        fPhysicInnerLAr,
fPhysicInnerVesselWall,         fOpLArVesselSurface ); new
G4LogicalBorderSurface("LArVesselSurfaceTop",  fPhysicLArBath,
fPhysicInnerVesselWindowTop,    fOpLArVesselSurface ); new
G4LogicalBorderSurface("LArVesselSurfaceWall", fPhysicLArBath,
fPhysicInnerVesselWindowBottom, fOpLArVesselSurface );
  fOpLArVesselSurface->SetType( dielectric_dielectric );
  fOpLArVesselSurface->SetModel( unified );
  fOpLArVesselSurface->SetFinish( ground );
  G4MaterialPropertiesTable *fLArVesselSurfProp = new
G4MaterialPropertiesTable(); fOpLArVesselSurface->SetMaterialPropertiesTable(
fLArVesselSurfProp );


  ////////////////////////////////////////
  // GAr - Acrylic(inner vessel)
  ////////////////////////////////////////
  G4OpticalSurface *fOpGArVesselSurface = new
G4OpticalSurface("OpGArVesselSurface"); new G4LogicalBorderSurface("GArVessel",
fPhysicGasPocket,   fPhysicInnerVesselWall,  fOpGArVesselSurface );
  fOpGArVesselSurface->SetType( dielectric_dielectric );
  fOpGArVesselSurface->SetModel( unified );
  fOpGArVesselSurface->SetFinish( ground );
  G4MaterialPropertiesTable *fGArVesselSurfProp = new
G4MaterialPropertiesTable(); fOpGArVesselSurface->SetMaterialPropertiesTable(
fGArVesselSurfProp );


  ////////////////////////////////////////
  // TPB - 3Mfoil (reflector)
  ////////////////////////////////////////
  G4OpticalSurface *fOpTPBThreeMSurface = new
G4OpticalSurface("OpTPBThreeMSurface"); new
G4LogicalBorderSurface("TPBThreeMSurfaceTop", fPhysicTPBTopLateral,
fPhysicThreeMTopFoil,    fOpTPBThreeMSurface ); new
G4LogicalBorderSurface("TPBThreeMSurfaceMid", fPhysicTPBMiddleLateral,
fPhysicThreeMMiddleFoil, fOpTPBThreeMSurface ); new
G4LogicalBorderSurface("TPBThreeMSurfaceBot", fPhysicTPBLowerLateral,
fPhysicThreeMLowerFoil,  fOpTPBThreeMSurface ); fOpTPBThreeMSurface->SetType(
dielectric_metal ); fOpTPBThreeMSurface->SetModel( unified );
  fOpTPBThreeMSurface->SetFinish( ground );
  fOpTPBThreeMSurface->SetSigmaAlpha(0.5);
  G4double TPBThreeMENE[2] = {0.1*eV, 20.0*eV};
  G4double TPBThreeMREF[2] = {0.99, 0.99};
  G4double TPBThreeMEFF[2] = {0.00, 0.00};
  G4MaterialPropertiesTable *fTPBThreeMSurfProp = new
G4MaterialPropertiesTable(); fTPBThreeMSurfProp->AddProperty("REFLECTIVITY",
TPBThreeMENE, TPBThreeMREF, 2); fTPBThreeMSurfProp->AddProperty("EFFICIENCY",
TPBThreeMENE, TPBThreeMEFF, 2); fOpTPBThreeMSurface->SetMaterialPropertiesTable(
fTPBThreeMSurfProp );


  ////////////////////////////////////////
  // LAr - 3Mfoil (reflector)
  ////////////////////////////////////////
  G4OpticalSurface *fOpLArThreeMSurface = new
G4OpticalSurface("OpLArThreeMSurface"); new
G4LogicalBorderSurface("LArThreeMSurfaceTop", fPhysicInnerLAr,
fPhysicThreeMTopFoil,    fOpLArThreeMSurface ); new
G4LogicalBorderSurface("LArThreeMSurfaceMid", fPhysicInnerLAr,
fPhysicThreeMMiddleFoil, fOpLArThreeMSurface ); new
G4LogicalBorderSurface("LArThreeMSurfaceBot", fPhysicInnerLAr,
fPhysicThreeMLowerFoil,  fOpLArThreeMSurface ); fOpLArThreeMSurface->SetType(
dielectric_metal ); fOpLArThreeMSurface->SetModel( unified );
  fOpLArThreeMSurface->SetFinish( ground );
  fOpLArThreeMSurface->SetSigmaAlpha(0.5);
  G4double LArThreeMENE[2] = {0.1*eV, 20.0*eV};
  G4double LArThreeMREF[2] = {0.97, 0.97};
  G4double LArThreeMEFF[2] = {0.00, 0.00};
  G4MaterialPropertiesTable *fLArThreeMSurfProp = new
G4MaterialPropertiesTable(); fLArThreeMSurfProp->AddProperty("REFLECTIVITY",
LArThreeMENE, LArThreeMREF, 2); fLArThreeMSurfProp->AddProperty("EFFICIENCY",
LArThreeMENE, LArThreeMEFF, 2); fOpLArThreeMSurface->SetMaterialPropertiesTable(
fLArThreeMSurfProp );


  ////////////////////////////////////////
  // Teflon - 3Mfoil   [NOT NEEDED ?]
  ////////////////////////////////////////
  G4OpticalSurface *fOpTeflonThreeMSurface = new
G4OpticalSurface("OpTeflonThreeMSurface"); new
G4LogicalBorderSurface("TeflonThreeMSurfaceGas", fPhysicThreeMTopFoil,
fPhysicTopGasSupportRing, fOpTeflonThreeMSurface ); new
G4LogicalBorderSurface("TeflonThreeMSurfaceLiq", fPhysicThreeMTopFoil,
fPhysicTopLiqSupportRing, fOpTeflonThreeMSurface ); new
G4LogicalBorderSurface("TeflonThreeMSurfaceR0",  fPhysicThreeMLowerFoil,
fPhysicSupportRing0,    fOpTeflonThreeMSurface );  // da verificare !!! new
G4LogicalBorderSurface("TeflonThreeMSurfaceR1",  fPhysicThreeMMiddleFoil,
fPhysicSupportRing1,   fOpTeflonThreeMSurface );// da verificare !!! new
G4LogicalBorderSurface("TeflonThreeMSurfaceR2",  fPhysicThreeMTopFoil,
fPhysicSupportRing2,      fOpTeflonThreeMSurface ); // da verificare !!!
  fOpTeflonThreeMSurface->SetType( dielectric_metal );
  fOpTeflonThreeMSurface->SetModel( unified );
  fOpTeflonThreeMSurface->SetFinish( ground );
  fOpTeflonThreeMSurface->SetSigmaAlpha(0.1);
  G4double TeflonThreeMENE[2] = {0.1*eV, 20.0*eV};
  G4double TeflonThreeMREF[2] = {0.99, 0.99};
  G4double TeflonThreeMEFF[2] = {0.00, 0.00};
  G4MaterialPropertiesTable *fTeflonThreeMSurfProp = new
G4MaterialPropertiesTable(); fTeflonThreeMSurfProp->AddProperty("REFLECTIVITY",
TeflonThreeMENE, TeflonThreeMREF, 2);
  fTeflonThreeMSurfProp->AddProperty("EFFICIENCY",   TeflonThreeMENE,
TeflonThreeMEFF, 2); fOpTeflonThreeMSurface->SetMaterialPropertiesTable(
fTeflonThreeMSurfProp );


  ////////////////////////////////////////
  // Acrylic - Copper
  ////////////////////////////////////////
  G4OpticalSurface *fOpAcrylicCopperSurface = new
G4OpticalSurface("OpAcrylicCopperSurface"); new
G4LogicalBorderSurface("AcrylicCopperSurface", fPhysicInnerVesselWall,
fPhysicFieldRings, fOpAcrylicCopperSurface ); fOpAcrylicCopperSurface->SetType(
dielectric_metal ); fOpAcrylicCopperSurface->SetModel( unified );
  fOpAcrylicCopperSurface->SetFinish( ground );
  fOpAcrylicCopperSurface->SetSigmaAlpha(0.1);
  G4double AcrylicCopperENE[2] = {0.1*eV, 20.0*eV};
  G4double AcrylicCopperREF[2] = {0.99, 0.99};
  G4double AcrylicCopperEFF[2] = {0.00, 0.00};
  G4MaterialPropertiesTable *fAcrylicCopperSurfProp = new
G4MaterialPropertiesTable(); fAcrylicCopperSurfProp->AddProperty("REFLECTIVITY",
AcrylicCopperENE, AcrylicCopperREF, 2);
  fAcrylicCopperSurfProp->AddProperty("EFFICIENCY",   AcrylicCopperENE,
AcrylicCopperEFF, 2); fOpAcrylicCopperSurface->SetMaterialPropertiesTable(
fAcrylicCopperSurfProp );


  ////////////////////////////////////////
  // LAr - Teflon
  ////////////////////////////////////////
  G4OpticalSurface *fOpLArTeflonSurface     = new
G4OpticalSurface("OpLArTeflonSurface"); new
G4LogicalBorderSurface("LArTeflonSurfaceTop",   fPhysicInnerLAr,
fPhysicTopLiqSupportRing, fOpLArTeflonSurface); new
G4LogicalBorderSurface("LArTeflonSurfaceR0",    fPhysicInnerLAr,
fPhysicSupportRing0,      fOpLArTeflonSurface); new
G4LogicalBorderSurface("LArTeflonSurfaceR1",    fPhysicInnerLAr,
fPhysicSupportRing1,      fOpLArTeflonSurface); new
G4LogicalBorderSurface("LArTeflonSurfaceR2",    fPhysicInnerLAr,
fPhysicSupportRing2,      fOpLArTeflonSurface); new
G4LogicalBorderSurface("LArTeflonSurfaceBub1",  fPhysicLArBath,
fPhysicBubblerTube1,      fOpLArTeflonSurface); new
G4LogicalBorderSurface("LArTeflonSurfaceBub2",  fPhysicLArBath,
fPhysicBubblerTube2,      fOpLArTeflonSurface); fOpLArTeflonSurface->SetType(
dielectric_metal ); fOpLArTeflonSurface->SetModel( unified );
  fOpLArTeflonSurface->SetFinish( polished );
  fOpLArTeflonSurface->SetSigmaAlpha(0.1);
  G4MaterialPropertiesTable *fLArTeflonSurfProp = new
G4MaterialPropertiesTable(); G4double LArTeflonENE[2] = {0.1*eV, 20.0*eV};
  G4double LArTeflonREF[2] = {0.99, 0.99};
  G4double LArTeflonEFF[2] = {0.00, 0.00};
  fLArTeflonSurfProp->AddProperty("REFLECTIVITY", LArTeflonENE, LArTeflonREF,
2); fLArTeflonSurfProp->AddProperty("EFFICIENCY",   LArTeflonENE, LArTeflonEFF,
2); fOpLArTeflonSurface->SetMaterialPropertiesTable( fLArTeflonSurfProp );


  ////////////////////////////////////////
  // GAr - Teflon (reflector)
  ////////////////////////////////////////
  G4OpticalSurface *fOpGArTeflonSurface     = new
G4OpticalSurface("OpGArTeflonSurface"); new
G4LogicalBorderSurface("GArTeflonSurfaceTop",   fPhysicGasPocket,
fPhysicTopGasSupportRing, fOpGArTeflonSurface); fOpGArTeflonSurface->SetType(
dielectric_metal ); fOpGArTeflonSurface->SetModel( unified );
  fOpGArTeflonSurface->SetFinish( polished );
  fOpGArTeflonSurface->SetSigmaAlpha(0.1);
  G4MaterialPropertiesTable *fGArTeflonSurfProp = new
G4MaterialPropertiesTable(); G4double GArTeflonENE[2] = {0.1*eV, 20.0*eV};
  G4double GArTeflonREF[2] = {0.99, 0.99};
  G4double GArTeflonEFF[2] = {0.00, 0.00};
  fGArTeflonSurfProp->AddProperty("REFLECTIVITY", GArTeflonENE, GArTeflonREF,
2); fGArTeflonSurfProp->AddProperty("EFFICIENCY",   GArTeflonENE, GArTeflonEFF,
2); fOpGArTeflonSurface->SetMaterialPropertiesTable( fGArTeflonSurfProp );


  ////////////////////////////////////////
  // ITO /////
  ////////////////////////////////////////

  G4MaterialPropertiesTable *fITOSurfProp = new G4MaterialPropertiesTable();
  fITOSurfProp->AddConstProperty("DOITO",1);


  ////////////////////////////////////////
  // LAr (inner) - TPB
  ////////////////////////////////////////
  G4OpticalSurface *fOpTPBLArSurface = new G4OpticalSurface("OpTBPLArSurface");
  new G4LogicalBorderSurface("TPBLArSurfaceInnerBotBas", fPhysicInnerLAr,
fPhysicTPBBottomBase,    fOpTPBLArSurface ); new
G4LogicalBorderSurface("TPBLArSurfaceInnerLowLat", fPhysicInnerLAr,
fPhysicTPBLowerLateral,  fOpTPBLArSurface ); new
G4LogicalBorderSurface("TPBLArSurfaceInnerMidLat", fPhysicInnerLAr,
fPhysicTPBMiddleLateral, fOpTPBLArSurface ); new
G4LogicalBorderSurface("TPBLArSurfaceInnerTopLat", fPhysicInnerLAr,
fPhysicTPBTopLateral,    fOpTPBLArSurface );   // not needed ? new
G4LogicalBorderSurface("TPBLArSurfaceInnerTopBas", fPhysicInnerLAr,
fPhysicTPBTopBase,       fOpTPBLArSurface );  // not needed ?
  fOpTPBLArSurface->SetType( dielectric_dielectric );
  fOpTPBLArSurface->SetModel( unified );
  fOpTPBLArSurface->SetFinish( ground );
  fOpTPBLArSurface->SetSigmaAlpha(0.1);
  G4MaterialPropertiesTable *fTPBLArSurfProp = new G4MaterialPropertiesTable();
  fOpTPBLArSurface->SetMaterialPropertiesTable( fTPBLArSurfProp );


  ////////////////////////////////////////
  // PMT surfaces
  ////////////////////////////////////////

  // Photocathode - LAr (bath)
  G4OpticalSurface *fOpPMTLArSurface = new G4OpticalSurface("OpPMTLArSurface");
  G4LogicalBorderSurface*     fPMTLArSurface[14];
  G4MaterialPropertiesTable*  fPMTLArSurfProp;
   for(int i = 0; i < 14; i++)  fPMTLArSurface[i] = new
G4LogicalBorderSurface("PMTLArSurface", fPhysicLArBath, fPhysicPMTWindow[i],
fOpPMTLArSurface ); fOpPMTLArSurface->SetType( dielectric_dielectric );
  fOpPMTLArSurface->SetModel( unified );
  fOpPMTLArSurface->SetFinish( polished );
  fPMTLArSurfProp = new G4MaterialPropertiesTable();
  fPMTLArSurfProp->AddConstProperty("REFLECTIVITY", 0.2);
  fPMTLArSurfProp->AddConstProperty("EFFICIENCY",   0.0);
  fOpPMTLArSurface->SetMaterialPropertiesTable( fPMTLArSurfProp );


  // LAr (bath) - Kovar (PMT body)
  G4OpticalSurface *fOpLArKovarSurface = new
G4OpticalSurface("OpLArKovarSurface"); G4LogicalBorderSurface*
fLArKovarSurface[14]; for(int i = 0; i < 14; i++)  fLArKovarSurface[i] = new
G4LogicalBorderSurface("LArKovarSurface", fPhysicLArBath, fPhysicPMTBody[i],
fOpLArKovarSurface ); fOpLArKovarSurface->SetType( dielectric_metal );
  fOpLArKovarSurface->SetModel( unified );
  fOpLArKovarSurface->SetFinish( polished );
  fOpLArKovarSurface->SetSigmaAlpha(0.1);
  G4MaterialPropertiesTable *fLArKovarSurfProp = new
G4MaterialPropertiesTable(); G4double LArKovarENE[2] = {0.1*eV, 20.0*eV};
  G4double LArKovarREF[2] = {1.00, 1.00};
  G4double LArKovarEFF[2] = {0.00, 0.00};
  fLArKovarSurfProp->AddProperty("REFLECTIVITY", LArKovarENE, LArKovarREF, 2);
  fLArKovarSurfProp->AddProperty("EFFICIENCY",   LArKovarENE, LArKovarEFF, 2);
*/
}

G4double ArDM::extended_asin (G4double x) {
  if (x < -1) return - M_PI / 2. ;
  if (x >  1) return M_PI / 2. ;
  return asin (x);
}


void ArDM::Init() {

  HANDEDNESS = -1;

  PI = M_PI ; // M_PI;
  BOLTZMAN_K = 1.380649e-23 ; // TMath::K();  // joule /kelvin              ; //boltzmann constant;

  NPMT = 12;
  NTOPPMT = NPMT;
  NBOTTOMPMT = NPMT;
  ELECTRIC_FIELD_STRENGTH = (0 * volt / cm);
  EXTRACTION_FIELD_STRENGTH = (4 * kilovolt / cm);
  DIELECTRIC_CONSTANT_GAR = (1.);
  DIELECTRIC_CONSTANT_LAR = (1.4);
  PHT_ENE_THRESHOLD = (6. * eV);  // if photon energy is larger than this
                                  // threshold --> it won't be detected !;
  BOILING_POINT_ARGON = (87. * kelvin);

  COMPONENT_RATIO_LAR_ELECTRON_LIKE = (0.23);
  COMPONENT_RATIO_LAR_NEUTRON_LIKE = (0.75);
  COMPONENT_RATIO_GAR_ELECTRON_LIKE = (0.4);
  COMPONENT_RATIO_GAR_NEUTRON_LIKE = (0.4);

  PMTCONV = (.727);
  REFLCONV = (.932);

  WORLD_HALF_SIZE = (4100. * mm);
  LID_RADIUS = (555. * mm);
  LID_HALF_HEIGHT = (20. * mm);
  TANK_HALF_HEIGHT = (2093. * mm / 2);
  TANK_CYLINDER_INNER_RADIUS = (500. * mm);
  TANK_CYLINDER_THICKNESS = (20. * mm);
  TANK_CYLINDER_OUTER_RADIUS = (TANK_CYLINDER_INNER_RADIUS + TANK_CYLINDER_THICKNESS);
  TANK_BTM_PART_R1010_ARC_INNER_RADIUS = (1010. * mm);
  TANK_BTM_PART_R1010_ARC_OUTER_RADIUS = (TANK_BTM_PART_R1010_ARC_INNER_RADIUS + TANK_CYLINDER_THICKNESS);
  TANK_BTM_PART_R101_ARC_INNER_RADIUS = (101. * mm);
  TANK_BTM_PART_R101_ARC_OUTER_RADIUS = (TANK_BTM_PART_R101_ARC_INNER_RADIUS + TANK_CYLINDER_THICKNESS);
  TANK_BTM_PART_R1010_ARC_CHORD = (2 * (TANK_CYLINDER_INNER_RADIUS - TANK_BTM_PART_R101_ARC_INNER_RADIUS));
  TANK_BTM_PART_R1010_ARC_OPENING_ANGLE = (2 * extended_asin(TANK_BTM_PART_R1010_ARC_CHORD / 2 / (TANK_BTM_PART_R1010_ARC_INNER_RADIUS - TANK_BTM_PART_R101_ARC_INNER_RADIUS)) / M_PI * 180 * deg);
  TANK_BTM_PART_R101_ARC_OPENING_ANGLE = (90 * deg - TANK_BTM_PART_R1010_ARC_OPENING_ANGLE / 2);
  DISTANCE_R1010_ARC_CENTER_TO_LOWER_EDGE_OF_TANK_CYLINDER = ((TANK_BTM_PART_R1010_ARC_INNER_RADIUS - TANK_BTM_PART_R101_ARC_INNER_RADIUS) * cos(TANK_BTM_PART_R1010_ARC_OPENING_ANGLE / 2));
  DISTANCE_LOWER_EDGE_OF_TANK_CYLINDER_TO_TANK_LOWERMOST_POINT = (TANK_BTM_PART_R1010_ARC_INNER_RADIUS - DISTANCE_R1010_ARC_CENTER_TO_LOWER_EDGE_OF_TANK_CYLINDER);
  TANK_CYLINDER_HALF_HEIGHT = (TANK_HALF_HEIGHT - DISTANCE_LOWER_EDGE_OF_TANK_CYLINDER_TO_TANK_LOWERMOST_POINT / 2);
  TANK_BTM_PART_R101_ARC_HEIGHT = (TANK_BTM_PART_R101_ARC_INNER_RADIUS * cos(TANK_BTM_PART_R1010_ARC_OPENING_ANGLE / 2));
  TANK_CYLINDER_POS_Z = (TANK_HALF_HEIGHT - TANK_CYLINDER_HALF_HEIGHT);
  TANK_BTM_PART_R1010_Z = (TANK_BTM_PART_R1010_ARC_INNER_RADIUS - TANK_HALF_HEIGHT);
  TANK_BTM_PART_R1010_Z_RELATIVE_TO_TANK_CYLINDER = (TANK_BTM_PART_R1010_Z - TANK_CYLINDER_POS_Z);
  DISTANCE_R1010_ARC_CENTER_TO_R101_ARC_CENTER = (DISTANCE_R1010_ARC_CENTER_TO_LOWER_EDGE_OF_TANK_CYLINDER);
  DISTANCE_TANK_LOWERMOST_POINT_TO_R101_ARC_CENTER = (TANK_BTM_PART_R1010_ARC_INNER_RADIUS - DISTANCE_R1010_ARC_CENTER_TO_R101_ARC_CENTER);
  DISTANCE_TANK_CYLINDER_CENTER_TO_R101_ARC_CENTER = (TANK_CYLINDER_HALF_HEIGHT);
  TANK_BTM_PART_R101_Z = (TANK_CYLINDER_POS_Z - DISTANCE_TANK_CYLINDER_CENTER_TO_R101_ARC_CENTER);
  TANK_BTM_PART_R101_Z_RELATIVE_TO_TANK_CYLINDER = (TANK_BTM_PART_R101_Z - TANK_CYLINDER_POS_Z);
  TOP_FLANGE_INNER_RADIUS = (0. * mm);
  TOP_FLANGE_OUTER_RADIUS = (555. * mm);
  TOP_FLANGE_HALF_HEIGHT = (40. * mm / 2);
  TOP_FLANGE_COMPENSATION_FACTOR_FOR_OMITTED_PILLARS = (0.9803485);
  TOP_FLANGE_HALF_HEIGHT_EFFECTIVE = (TOP_FLANGE_HALF_HEIGHT * TOP_FLANGE_COMPENSATION_FACTOR_FOR_OMITTED_PILLARS);
  TOP_FLANGE_POS_X = 0. * mm;
  TOP_FLANGE_POS_Y = 0. * mm;
  TOP_FLANGE_POS_Z = (TANK_CYLINDER_POS_Z + TANK_CYLINDER_HALF_HEIGHT + TOP_FLANGE_HALF_HEIGHT_EFFECTIVE);
  DISTANCE_TOP_FLANGE_TOP_PMT_CENTER = (249.2 * mm);
  DISTANCE_UPPER_EDGE_FIELD_SHAPER_LIQUID_SURFACE = (70. * mm);
  REFLECTOR_HALF_HEIGHT = (538.5 * mm);
  DISTANCE_TOP_FLANGE_UPPER_EDGE_FIELD_SHAPER = (473. * mm);
  DISTANCE_TOP_FLANGE_LIQUID_SURFACE = (DISTANCE_TOP_FLANGE_UPPER_EDGE_FIELD_SHAPER - DISTANCE_UPPER_EDGE_FIELD_SHAPER_LIQUID_SURFACE);
  LIQUID_SURFACE_POS_Z = (TOP_FLANGE_POS_Z - TOP_FLANGE_HALF_HEIGHT_EFFECTIVE - DISTANCE_TOP_FLANGE_LIQUID_SURFACE);
  DISTANCE_TOP_FLANGE_REFLECTOR_CENTER = (DISTANCE_TOP_FLANGE_UPPER_EDGE_FIELD_SHAPER + REFLECTOR_HALF_HEIGHT);
  DISTANCE_REFLECTOR_LOWER_EDGE_PROTECTION_GRID = (133. * mm);  //<-- check it ! ;
  DISTANCE_PROTECTION_GRID_TIP_OF_BTM_PMT = (21. * mm);         //<-- check it !;
  GAR_COLUMN_HALF_HEIGHT = (DISTANCE_TOP_FLANGE_LIQUID_SURFACE / 2);
  LAR_COLUMN_HALF_HEIGHT = (TANK_HALF_HEIGHT - GAR_COLUMN_HALF_HEIGHT);
  LAR_COLUMN_Z = (-(TANK_HALF_HEIGHT - LAR_COLUMN_HALF_HEIGHT));
  GAR_COLUMN_Z = (TANK_HALF_HEIGHT - GAR_COLUMN_HALF_HEIGHT);
  LARCOL_TOTAL_HEIGHT = (TANK_HALF_HEIGHT * 2 - DISTANCE_TOP_FLANGE_LIQUID_SURFACE);
  LARCOL_CYLINDER_HALFHEIGHT = ((LARCOL_TOTAL_HEIGHT - DISTANCE_LOWER_EDGE_OF_TANK_CYLINDER_TO_TANK_LOWERMOST_POINT) / 2);
  LARCOL_CYLINDER_OUTER_RADIUS = TANK_CYLINDER_INNER_RADIUS;
  LARCOL_CYLINDER_POS_Z = (-TANK_HALF_HEIGHT + DISTANCE_LOWER_EDGE_OF_TANK_CYLINDER_TO_TANK_LOWERMOST_POINT + LARCOL_CYLINDER_HALFHEIGHT);
  WLS_MEAN_ABSORPTION_LENGTH = (.01 * mm);
  WLS_THICKNESS_100_PERCENT_CONVERSION_EFFICIENCY = (0.2 * mm);
  WLS_RINGSEC_THICKNESS = getWLSThickness(REFLCONV);
  WLS_RINGSEC_INNER_RADIUS = (791. * mm / 2);
  WLS_RINGSEC_OUTER_RADIUS = (WLS_RINGSEC_INNER_RADIUS + WLS_RINGSEC_THICKNESS);
  WLS_RINGSEC_HALF_HEIGHT = REFLECTOR_HALF_HEIGHT;
  WLS_RINGSEC_START_PHI = (40. * deg);
  WLS_RINGSEC_DELTA_PHI = (280. * deg);
  WLS_RINGSEC_RADIUS = (WLS_RINGSEC_INNER_RADIUS + WLS_RINGSEC_OUTER_RADIUS) / 2;
  WLS_RINGSEC_POS_Z = (TANK_HALF_HEIGHT - DISTANCE_TOP_FLANGE_REFLECTOR_CENTER);
  WLS_RINGSEC_POS_Z_LARCOL = (WLS_RINGSEC_POS_Z - LARCOL_CYLINDER_POS_Z);
  WLS_LINSEC_HALF_X = (WLS_RINGSEC_RADIUS * sin(40. * deg));
  WLS_LINSEC_HALF_Y = (WLS_RINGSEC_OUTER_RADIUS - WLS_RINGSEC_INNER_RADIUS) * cos(40. * deg) / 2;
  WLS_LINSEC_HALF_Z = WLS_RINGSEC_HALF_HEIGHT;
  WLS_LINSEC_HALF_THICKNESS = WLS_LINSEC_HALF_Y;
  WLS_LINSEC_POS_X = WLS_RINGSEC_RADIUS * cos(40. * deg);
  WLS_LINSEC_POS_Y = (0. * mm);
  WLS_LINSEC_POS_Z = WLS_RINGSEC_POS_Z;
  WLS_SUPPORT_RINGSEC_INNER_RADIUS = WLS_RINGSEC_OUTER_RADIUS;
  WLS_SUPPORT_RINGSEC_THICKNESS = (.3 * mm);
  WLS_SUPPORT_RINGSEC_OUTER_RADIUS = (WLS_SUPPORT_RINGSEC_INNER_RADIUS + WLS_SUPPORT_RINGSEC_THICKNESS);
  WLS_SUPPORT_RINGSEC_HALF_HEIGHT = WLS_RINGSEC_HALF_HEIGHT;
  WLS_SUPPORT_RINGSEC_START_PHI = WLS_RINGSEC_START_PHI;
  WLS_SUPPORT_RINGSEC_DELTA_PHI = WLS_RINGSEC_DELTA_PHI;
  WLS_SUPPORT_RINGSEC_RADIUS = (WLS_SUPPORT_RINGSEC_INNER_RADIUS + WLS_SUPPORT_RINGSEC_OUTER_RADIUS) / 2;
  WLS_SUPPORT_RINGSEC_POS_Z = WLS_RINGSEC_POS_Z;
  WLS_SUPPORT_LINSEC_HALF_X = (WLS_SUPPORT_RINGSEC_RADIUS * sin(40. * deg));
  WLS_SUPPORT_LINSEC_HALF_Y = (WLS_SUPPORT_RINGSEC_OUTER_RADIUS - WLS_SUPPORT_RINGSEC_INNER_RADIUS) * cos(40. * deg) / 2;
  WLS_SUPPORT_LINSEC_HALF_Z = WLS_SUPPORT_RINGSEC_HALF_HEIGHT;
  WLS_SUPPORT_LINSEC_POS_X = (WLS_SUPPORT_RINGSEC_RADIUS * cos(40. * deg));
  WLS_SUPPORT_LINSEC_POS_Y = WLS_LINSEC_POS_Y;
  WLS_SUPPORT_LINSEC_POS_Z = WLS_LINSEC_POS_Z;
  WLS_UPPER_EDGE_POS_Z = (TOP_FLANGE_POS_Z - TOP_FLANGE_HALF_HEIGHT_EFFECTIVE - DISTANCE_TOP_FLANGE_UPPER_EDGE_FIELD_SHAPER);
  WLS_LOWER_EDGE_POS_Z = (WLS_UPPER_EDGE_POS_Z - 2 * WLS_RINGSEC_HALF_HEIGHT);
  WLS_RINGSEC_HALF_HEIGHT_GAR = ((WLS_UPPER_EDGE_POS_Z - LIQUID_SURFACE_POS_Z) / 2);
  WLS_RINGSEC_GAR_POS_Z = (WLS_UPPER_EDGE_POS_Z - WLS_RINGSEC_HALF_HEIGHT_GAR);
  WLS_RINGSEC_GAR_POS_Z_GARCOL = (WLS_RINGSEC_GAR_POS_Z - GAR_COLUMN_Z);
  WLS_LINSEC_HALF_HEIGHT_GAR = WLS_RINGSEC_HALF_HEIGHT_GAR;
  WLS_LINSEC_GAR_POS_X = WLS_LINSEC_POS_X;
  WLS_LINSEC_GAR_POS_Y = WLS_LINSEC_POS_Y;
  WLS_LINSEC_GAR_POS_Z = WLS_RINGSEC_GAR_POS_Z;
  WLS_LINSEC_GAR_POS_Z_GARCOL = WLS_RINGSEC_GAR_POS_Z_GARCOL;
  WLS_SUPPORT_RINGSEC_HALF_HEIGHT_GAR = WLS_RINGSEC_HALF_HEIGHT_GAR;
  WLS_SUPPORT_RINGSEC_GARCOL_POS_Z = WLS_RINGSEC_GAR_POS_Z;
  WLS_SUPPORT_RINGSEC_GARCOL_POS_Z_GARCOL = WLS_RINGSEC_GAR_POS_Z_GARCOL;
  WLS_SUPPORT_LINSEC_HALF_HEIGHT_GAR = WLS_SUPPORT_RINGSEC_HALF_HEIGHT_GAR;
  WLS_SUPPORT_LINSEC_GAR_POS_X = WLS_SUPPORT_LINSEC_HALF_X;
  WLS_SUPPORT_LINSEC_GAR_POS_Y = WLS_SUPPORT_LINSEC_HALF_Y;
  WLS_SUPPORT_LINSEC_GAR_POS_Z = WLS_SUPPORT_RINGSEC_GARCOL_POS_Z;
  WLS_SUPPORT_LINSEC_GAR_POS_Z_GARCOL = WLS_SUPPORT_RINGSEC_GARCOL_POS_Z_GARCOL;
  WLS_RINGSEC_HALF_HEIGHT_LAR = ((WLS_RINGSEC_HALF_HEIGHT_GAR > 0) ? (WLS_RINGSEC_HALF_HEIGHT - WLS_RINGSEC_HALF_HEIGHT_GAR) : WLS_RINGSEC_HALF_HEIGHT);
  WLS_RINGSEC_LAR_POS_Z = (WLS_LOWER_EDGE_POS_Z + WLS_RINGSEC_HALF_HEIGHT_LAR);
  WLS_RINGSEC_LAR_POS_Z_LARCOL = (WLS_RINGSEC_LAR_POS_Z - LARCOL_CYLINDER_POS_Z);
  WLS_LINSEC_HALF_HEIGHT_LAR = WLS_RINGSEC_HALF_HEIGHT_LAR;
  WLS_LINSEC_LAR_POS_X = WLS_LINSEC_POS_X;
  WLS_LINSEC_LAR_POS_Y = WLS_LINSEC_POS_Y;
  WLS_LINSEC_LAR_POS_Z = WLS_RINGSEC_LAR_POS_Z;
  WLS_LINSEC_LAR_POS_Z_LARCOL = WLS_RINGSEC_LAR_POS_Z_LARCOL;
  WLS_SUPPORT_RINGSEC_HALF_HEIGHT_LAR = WLS_RINGSEC_HALF_HEIGHT_LAR;
  WLS_SUPPORT_RINGSEC_LARCOL_POS_Z = WLS_RINGSEC_LAR_POS_Z;
  WLS_SUPPORT_RINGSEC_LARCOL_POS_Z_LARCOL = WLS_RINGSEC_LAR_POS_Z_LARCOL;
  WLS_SUPPORT_LINSEC_HALF_HEIGHT_LAR = WLS_SUPPORT_RINGSEC_HALF_HEIGHT_LAR;
  WLS_SUPPORT_LINSEC_LAR_POS_X = WLS_SUPPORT_LINSEC_HALF_X;
  WLS_SUPPORT_LINSEC_LAR_POS_Y = WLS_SUPPORT_LINSEC_HALF_Y;
  WLS_SUPPORT_LINSEC_LAR_POS_Z = WLS_SUPPORT_RINGSEC_LARCOL_POS_Z;
  WLS_SUPPORT_LINSEC_LAR_POS_Z_LARCOL = WLS_SUPPORT_RINGSEC_LARCOL_POS_Z_LARCOL;
  PMT_INNER_RADIUS = (131. * mm);
  PMT_THICK = (1. * mm);
  PMT_OUTER_RADIUS = (PMT_INNER_RADIUS + PMT_THICK);  // thickness of PMT = 1mm ??;
  PMT_ACTIVE_RANGE = (46.5 * deg);
  PMT_SPHERICAL_PART_INNER_RADIUS = (131. * mm);
  PMT_SPHERICAL_PART_OUTER_RADIUS = (PMT_SPHERICAL_PART_INNER_RADIUS + PMT_THICK);
  PMT_SPHERICAL_PART_OPENING_ANGLE = (2 * 50. * deg);
  PMT_MIDDLE_CYLINDER_INNER_RADIUS = (PMT_SPHERICAL_PART_INNER_RADIUS * sin(PMT_SPHERICAL_PART_OPENING_ANGLE / 2));
  PMT_MIDDLE_CYLINDER_OUTER_RADIUS = (PMT_MIDDLE_CYLINDER_INNER_RADIUS + PMT_THICK);
  PMT_MIDDLE_CYLINDER_HALF_HEIGHT = (27. * mm);
  PMT_SPHERICAL_PART_HEIGHT = (PMT_SPHERICAL_PART_OUTER_RADIUS - PMT_SPHERICAL_PART_INNER_RADIUS * cos(PMT_SPHERICAL_PART_OPENING_ANGLE / 2));
  PMT_BTM_TUBE_INNER_RADIUS = (42.25 * mm);
  PMT_BTM_TUBE_OUTER_RADIUS = (PMT_BTM_TUBE_INNER_RADIUS + PMT_THICK);
  PMT_BTM_TUBE_HALF_HEIGHT = (36. * mm);
  PMT_BTM_SPHERE_HOLE_OPENING_ANGLE = (2 * extended_asin(PMT_BTM_TUBE_INNER_RADIUS / PMT_SPHERICAL_PART_OUTER_RADIUS) / M_PI * 180 * deg);
  PMT_DISTANCE_PMT_CENTER_TO_SPHERE_CENTER = (fabs((PMT_SPHERICAL_PART_INNER_RADIUS * cos(PMT_SPHERICAL_PART_OPENING_ANGLE / 2) - PMT_MIDDLE_CYLINDER_HALF_HEIGHT)));
  PMT_DISTANCE_PMT_CENTER_TO_TOP_SPHERE_CENTER = PMT_DISTANCE_PMT_CENTER_TO_SPHERE_CENTER;
  PMT_DISTANCE_PMT_CENTER_TO_BTM_SPHERE_CENTER = PMT_DISTANCE_PMT_CENTER_TO_SPHERE_CENTER;
  PMT_DISTANCE_PMT_CENTER_TO_BTM_TUBE_CENTER = (PMT_SPHERICAL_PART_OUTER_RADIUS * cos(PMT_BTM_SPHERE_HOLE_OPENING_ANGLE / 2) - PMT_DISTANCE_PMT_CENTER_TO_TOP_SPHERE_CENTER + PMT_BTM_TUBE_HALF_HEIGHT);
  DISTANCE_TOP_FLANGE_UPPER_EDGE_OF_TOP_PMT_TUBE = (160. * mm);
  TOP_PMT_Z = (TOP_FLANGE_POS_Z - TOP_FLANGE_HALF_HEIGHT_EFFECTIVE - DISTANCE_TOP_FLANGE_UPPER_EDGE_OF_TOP_PMT_TUBE - PMT_BTM_TUBE_HALF_HEIGHT - PMT_DISTANCE_PMT_CENTER_TO_BTM_TUBE_CENTER);
  TOP_PMT_Z_GARCOL = (TOP_PMT_Z - GAR_COLUMN_Z);
  BTM_PMT_Z = (WLS_RINGSEC_POS_Z - WLS_RINGSEC_HALF_HEIGHT - DISTANCE_REFLECTOR_LOWER_EDGE_PROTECTION_GRID - DISTANCE_PROTECTION_GRID_TIP_OF_BTM_PMT - PMT_SPHERICAL_PART_HEIGHT - PMT_MIDDLE_CYLINDER_HALF_HEIGHT);
  BTM_PMT_Z_LARCOL = (BTM_PMT_Z - LARCOL_CYLINDER_POS_Z);
  TOP_PMT_CENTER_OF_TOP_SPHERE = (TOP_PMT_Z + PMT_DISTANCE_PMT_CENTER_TO_TOP_SPHERE_CENTER);
  TOP_PMT_CENTER_OF_TOP_SPHERE_GARCOL = (TOP_PMT_CENTER_OF_TOP_SPHERE - GAR_COLUMN_Z);
  BTM_PMT_CENTER_OF_TOP_SPHERE = (BTM_PMT_Z - PMT_DISTANCE_PMT_CENTER_TO_TOP_SPHERE_CENTER);
  BTM_PMT_CENTER_OF_TOP_SPHERE_LARCOL = (BTM_PMT_CENTER_OF_TOP_SPHERE - LARCOL_CYLINDER_POS_Z);
  PMT_BASE_RADIUS = (30. * mm);
  PMT_BASE_HALF_THICKNESS = (1.6 * mm);
  PMT_DISTANCE_LOWER_EDGE_OF_PMT_TUBE_TO_PMT_BASE_CENTER = (5. * mm);
  PMT_DISTANCE_PMT_CENTER_TO_PMT_BASE = (PMT_DISTANCE_PMT_CENTER_TO_BTM_TUBE_CENTER + PMT_BTM_TUBE_HALF_HEIGHT + PMT_DISTANCE_LOWER_EDGE_OF_PMT_TUBE_TO_PMT_BASE_CENTER);
  PMT_BASE_Z_TOP = (TOP_PMT_Z + PMT_DISTANCE_PMT_CENTER_TO_PMT_BASE);
  PMT_BASE_Z_TOP_GARCOL = (PMT_BASE_Z_TOP - GAR_COLUMN_Z);
  PMT_BASE_Z_BTM = (BTM_PMT_Z - PMT_DISTANCE_PMT_CENTER_TO_PMT_BASE);
  PMT_BASE_Z_BTM_LARCOL = (PMT_BASE_Z_BTM - LARCOL_CYLINDER_POS_Z);
  PMT_ELECTRODE_OUTER_RADIUS = PMT_BTM_TUBE_INNER_RADIUS;
  PMT_ELECTRODE_HALF_HEIGHT = PMT_BTM_TUBE_HALF_HEIGHT;
  PMT_ELECTRODE_VOLUME = ((197. / 8.00) * cm3);  // weight/density, material: stainless steel;
  PMT_ELECTRODE_INNER_RADIUS = sqrt((M_PI * 2 * PMT_ELECTRODE_HALF_HEIGHT * PMT_ELECTRODE_OUTER_RADIUS * PMT_ELECTRODE_OUTER_RADIUS - PMT_ELECTRODE_VOLUME) / M_PI / 2 / PMT_ELECTRODE_HALF_HEIGHT);
  PMT_ELECTRODE_Z_TOP = (TOP_PMT_Z + PMT_DISTANCE_PMT_CENTER_TO_BTM_TUBE_CENTER);
  PMT_ELECTRODE_Z_TOP_GARCOL = (PMT_ELECTRODE_Z_TOP - GAR_COLUMN_Z);
  PMT_ELECTRODE_Z_BTM = (BTM_PMT_Z - PMT_DISTANCE_PMT_CENTER_TO_BTM_TUBE_CENTER);
  PMT_ELECTRODE_Z_BTM_LARCOL = (PMT_ELECTRODE_Z_BTM - LARCOL_CYLINDER_POS_Z);
  HV_RESISTOR_BAR_INNER_RADIUS = (0. * mm);
  HV_RESISTOR_BAR_OUTER_RADIUS = (4.1 * mm);
  HV_RESISTOR_BAR_HALF_HEIGHT = (WLS_RINGSEC_HALF_HEIGHT);
  DISTANCE_BETWEEN_HV_RESISTOR_BARS = (100. * mm);
  HV_RESISTOR_BAR_1_POS_X = (WLS_SUPPORT_LINSEC_POS_X + WLS_LINSEC_HALF_THICKNESS + HV_RESISTOR_BAR_OUTER_RADIUS);
  HV_RESISTOR_BAR_1_POS_Y = (WLS_LINSEC_POS_Y + DISTANCE_BETWEEN_HV_RESISTOR_BARS / 2);
  HV_RESISTOR_BAR_1_POS_Z = WLS_LINSEC_POS_Z;
  HV_RESISTOR_BAR_1_POS_Z_LARCOL = (HV_RESISTOR_BAR_1_POS_Z - LARCOL_CYLINDER_POS_Z);
  HV_RESISTOR_BAR_2_POS_X = HV_RESISTOR_BAR_1_POS_X;
  HV_RESISTOR_BAR_2_POS_Y = (HV_RESISTOR_BAR_1_POS_Y - DISTANCE_BETWEEN_HV_RESISTOR_BARS);
  HV_RESISTOR_BAR_2_POS_Z = HV_RESISTOR_BAR_1_POS_Z;
  HV_RESISTOR_BAR_2_POS_Z_LARCOL = (HV_RESISTOR_BAR_2_POS_Z - LARCOL_CYLINDER_POS_Z);
  PMT_COATING_THICKNESS = getWLSThickness(PMTCONV);
  PMT_COATING_INNER_RADIUS = PMT_SPHERICAL_PART_OUTER_RADIUS;
  PMT_COATING_OUTER_RADIUS = (PMT_COATING_INNER_RADIUS + PMT_COATING_THICKNESS);
  TOP_PMT_COATING_Z = (TOP_PMT_Z + PMT_DISTANCE_PMT_CENTER_TO_SPHERE_CENTER);
  TOP_PMT_COATING_Z_GARCOL = (TOP_PMT_COATING_Z - GAR_COLUMN_Z);
  BTM_PMT_COATING_Z = (BTM_PMT_Z - PMT_DISTANCE_PMT_CENTER_TO_SPHERE_CENTER);
  BTM_PMT_COATING_Z_LARCOL = (BTM_PMT_COATING_Z - LARCOL_CYLINDER_POS_Z);
  PMT_CATHODE_OUTER_RADIUS = PMT_SPHERICAL_PART_INNER_RADIUS;
  PMT_CATHODE_THICKNESS = (0.001 * mm);
  PMT_CATHODE_INNER_RADIUS = (PMT_CATHODE_OUTER_RADIUS - PMT_CATHODE_THICKNESS);
  PMT_CATHODE_ACTIVE_RANGE = (46.5 * deg);
  PMT_ACTIVE_RANGE_RADIUS = (95.0240426 * mm);  //(PMT_OUTER_RADIUS*sin(PMT_CATHODE_ACTIVE_RANGE));
  TOP_PMT_CATHODE_Z = TOP_PMT_COATING_Z;
  TOP_PMT_CATHODE_Z_GARCOL = (TOP_PMT_CATHODE_Z - GAR_COLUMN_Z);
  BTM_PMT_CATHODE_Z = BTM_PMT_COATING_Z;
  BTM_PMT_CATHODE_Z_LARCOL = (BTM_PMT_CATHODE_Z - LARCOL_CYLINDER_POS_Z);
  PMT_SUPPORT_RADIUS = (450 * mm);
  PMT_SUPPORT_HALF_HEIGHT = (2.5 * mm);
  TOP_PMT_SUPPORT_POS_Z = (TOP_PMT_Z - PMT_MIDDLE_CYLINDER_HALF_HEIGHT + PMT_SUPPORT_HALF_HEIGHT);
  TOP_PMT_SUPPORT_POS_Z_GARCOL = (TOP_PMT_SUPPORT_POS_Z - GAR_COLUMN_Z);
  BTM_PMT_SUPPORT_POS_Z = (BTM_PMT_Z + PMT_MIDDLE_CYLINDER_HALF_HEIGHT - PMT_SUPPORT_HALF_HEIGHT);
  BTM_PMT_SUPPORT_POS_Z_LARCOL = (BTM_PMT_SUPPORT_POS_Z - LARCOL_CYLINDER_POS_Z);
  TOP_PMT_SUPPORT_POS_LOWER_EDGE = (TOP_PMT_SUPPORT_POS_Z - PMT_SUPPORT_HALF_HEIGHT);
  BTM_PMT_SUPPORT_POS_UPPER_EDGE = (BTM_PMT_SUPPORT_POS_Z + PMT_SUPPORT_HALF_HEIGHT);
  PMT_SUPPORT_CONVERSION_EFFICIENCY = (1.);
  PMT_SUPPORT_COATING_HALF_HEIGHT = (getWLSThickness(PMT_SUPPORT_CONVERSION_EFFICIENCY) / 2.);
  TOP_PMT_SUPPORT_COATING_POS_Z = (TOP_PMT_SUPPORT_POS_Z - PMT_SUPPORT_HALF_HEIGHT - PMT_SUPPORT_COATING_HALF_HEIGHT);
  TOP_PMT_SUPPORT_COATING_POS_Z_GARCOL = (TOP_PMT_SUPPORT_COATING_POS_Z - GAR_COLUMN_Z);
  BTM_PMT_SUPPORT_COATING_POS_Z = (BTM_PMT_SUPPORT_POS_Z + PMT_SUPPORT_HALF_HEIGHT + PMT_SUPPORT_COATING_HALF_HEIGHT);
  BTM_PMT_SUPPORT_COATING_POS_Z_LARCOL = (BTM_PMT_SUPPORT_COATING_POS_Z - LARCOL_CYLINDER_POS_Z);
  CATHODE_POS_Z = (WLS_RINGSEC_POS_Z - WLS_RINGSEC_HALF_HEIGHT - 1 * mm);
  CATHODE_POS_Z_LARCOL = (CATHODE_POS_Z - LARCOL_CYLINDER_POS_Z);
  CATHODE_PLATE_INNER_RADIUS = WLS_SUPPORT_RINGSEC_OUTER_RADIUS;
  CATHODE_PLATE_OUTER_RADIUS = (418. * mm);
  CATHODE_PLATE_HALF_THICKNESS = (13.4 * mm / 2);
  CATHODE_WIRE_INNER_RADIUS = (0. * mm);
  CATHODE_WIRE_OUTER_RADIUS = (.25 * mm);
  CATHODE_WIRE_PITCH = (20. * mm);
  PMMA_PLATE_OUTER_RADIUS = CATHODE_PLATE_OUTER_RADIUS;  // same as the cathode;
  PMMA_PLATE_HALF_THICKNESS = (10 * mm / 2);
  ITO_HALF_THICKNESS = (0.1 * mm / 2);
  BTM_PROTECTION_GRID_PLATE_INNER_RADIUS = WLS_SUPPORT_RINGSEC_OUTER_RADIUS;
  BTM_PROTECTION_GRID_PLATE_OUTER_RADIUS = PMT_SUPPORT_RADIUS;
  BTM_PROTECTION_GRID_PLATE_HALF_THICKNESS = (2.5 * mm);
  DISTANCE_LOWER_EDGE_CATHODE_PLATE_TO_UPPER_EDGE_BTM_PROTECTION_GRID = (120. * mm);
  BTM_PROTECTION_GRID_WIRE_INNER_RADIUS = CATHODE_WIRE_INNER_RADIUS;
  BTM_PROTECTION_GRID_WIRE_OUTER_RADIUS = CATHODE_WIRE_OUTER_RADIUS;
  BTM_PROTECTION_GRID_WIRE_PITCH = CATHODE_WIRE_PITCH;
  BTM_PROTECTION_GRID_POS_Z = (CATHODE_POS_Z - CATHODE_PLATE_HALF_THICKNESS - DISTANCE_LOWER_EDGE_CATHODE_PLATE_TO_UPPER_EDGE_BTM_PROTECTION_GRID - BTM_PROTECTION_GRID_PLATE_HALF_THICKNESS);
  BTM_PROTECTION_GRID_POS_Z_LARCOL = (BTM_PROTECTION_GRID_POS_Z - LARCOL_CYLINDER_POS_Z);
  SIDE_REFLECTOR_THICKNESS = (.254 * mm);
  BTM_SIDE_REFLECTOR_HALF_HEIGHT = (DISTANCE_REFLECTOR_LOWER_EDGE_PROTECTION_GRID / 2);
  BTM_SIDE_REFLECTOR_START_PHI = (0. * deg);
  BTM_SIDE_REFLECTOR_DELTA_PHI = (360. * deg);
  BTM_SIDE_REFLECTOR_POS_Z = (WLS_RINGSEC_POS_Z - WLS_RINGSEC_HALF_HEIGHT - BTM_SIDE_REFLECTOR_HALF_HEIGHT);  // relative to tank !!;
  BTM_SIDE_REFLECTOR_POS_Z_LARCOL = (BTM_SIDE_REFLECTOR_POS_Z - LARCOL_CYLINDER_POS_Z);                       // relative to LArColumn !;
  BTM_SIDE_REFLECTOR_INNER_RADIUS_UPPER_EDGE = CATHODE_PLATE_OUTER_RADIUS;
  BTM_SIDE_REFLECTOR_OUTER_RADIUS_UPPER_EDGE = (BTM_SIDE_REFLECTOR_INNER_RADIUS_UPPER_EDGE + SIDE_REFLECTOR_THICKNESS);
  BTM_SIDE_REFLECTOR_INNER_RADIUS_LOWER_EDGE = BTM_PROTECTION_GRID_PLATE_OUTER_RADIUS;
  BTM_SIDE_REFLECTOR_OUTER_RADIUS_LOWER_EDGE = (BTM_SIDE_REFLECTOR_INNER_RADIUS_LOWER_EDGE + SIDE_REFLECTOR_THICKNESS);
  BTM_SIDE_REFLECTOR_HALF_HEIGHT_GASTEST = (((CATHODE_POS_Z - CATHODE_PLATE_HALF_THICKNESS) - (BTM_PROTECTION_GRID_POS_Z + BTM_PROTECTION_GRID_PLATE_HALF_THICKNESS)) / 2);
  BTM_SIDE_REFLECTOR_POS_Z_GASTEST = (CATHODE_POS_Z - CATHODE_PLATE_HALF_THICKNESS - BTM_SIDE_REFLECTOR_HALF_HEIGHT_GASTEST);
  BTM_SIDE_REFLECTOR_POS_Z_GASTEST_LARCOL = (BTM_SIDE_REFLECTOR_POS_Z_GASTEST - LARCOL_CYLINDER_POS_Z);
  BTM_SIDE_REFLECTOR_CYNLIDRICAL_PART_INNER_RADIUS = BTM_PROTECTION_GRID_PLATE_OUTER_RADIUS;
  BTM_SIDE_REFLECTOR_CYNLIDRICAL_PART_OUTER_RADIUS = (BTM_SIDE_REFLECTOR_CYNLIDRICAL_PART_INNER_RADIUS + SIDE_REFLECTOR_THICKNESS);
  BTM_SIDE_REFLECTOR_CYNLIDRICAL_PART_HALF_HEIGHT = (((BTM_PROTECTION_GRID_POS_Z + BTM_PROTECTION_GRID_PLATE_HALF_THICKNESS) - (BTM_PMT_SUPPORT_POS_Z + PMT_SUPPORT_HALF_HEIGHT)) / 2);
  BTM_SIDE_REFLECTOR_CYNLIDRICAL_PART_POS_Z = ((BTM_PROTECTION_GRID_POS_Z - BTM_PROTECTION_GRID_PLATE_HALF_THICKNESS) - BTM_SIDE_REFLECTOR_CYNLIDRICAL_PART_HALF_HEIGHT);
  BTM_SIDE_REFLECTOR_CYNLIDRICAL_PART_POS_Z_LARCOL = (BTM_SIDE_REFLECTOR_CYNLIDRICAL_PART_POS_Z - LARCOL_CYLINDER_POS_Z);
  TOP_SIDE_REFLECTOR_INNER_RADIUS_LOWER_EDGE = WLS_SUPPORT_RINGSEC_OUTER_RADIUS;
  TOP_SIDE_REFLECTOR_OUTER_RADIUS_LOWER_EDGE = (TOP_SIDE_REFLECTOR_INNER_RADIUS_LOWER_EDGE + SIDE_REFLECTOR_THICKNESS);
  TOP_SIDE_REFLECTOR_INNER_RADIUS_UPPER_EDGE = PMT_SUPPORT_RADIUS;
  TOP_SIDE_REFLECTOR_OUTER_RADIUS_UPPER_EDGE = (TOP_SIDE_REFLECTOR_INNER_RADIUS_UPPER_EDGE + SIDE_REFLECTOR_THICKNESS);
  TOP_SIDE_REFLECTOR_HALF_HEIGHT = ((TOP_PMT_SUPPORT_POS_Z - PMT_SUPPORT_HALF_HEIGHT - (WLS_RINGSEC_POS_Z + WLS_RINGSEC_HALF_HEIGHT)) / 2);
  TOP_SIDE_REFLECTOR_START_PHI = (0. * deg);
  TOP_SIDE_REFLECTOR_DELTA_PHI = (360. * deg);
  TOP_SIDE_REFLECTOR_POS_Z = (WLS_RINGSEC_POS_Z + WLS_RINGSEC_HALF_HEIGHT + TOP_SIDE_REFLECTOR_HALF_HEIGHT);
  TOP_SIDE_REFLECTOR_POS_Z_GARCOL = (TOP_SIDE_REFLECTOR_POS_Z - GAR_COLUMN_Z);
  TOP_SIDE_REFLECTOR_HALF_HEIGHT_IN_GARCOL = ((TOP_PMT_SUPPORT_POS_Z - PMT_SUPPORT_HALF_HEIGHT - ((DISTANCE_UPPER_EDGE_FIELD_SHAPER_LIQUID_SURFACE > 0) ? LIQUID_SURFACE_POS_Z : WLS_UPPER_EDGE_POS_Z)) / 2);
  TOP_SIDE_REFLECTOR_HALF_HEIGHT_IN_LARCOL = (TOP_SIDE_REFLECTOR_HALF_HEIGHT - TOP_SIDE_REFLECTOR_HALF_HEIGHT_IN_GARCOL);
  TOP_SIDE_REFLECTOR_COTAN_TILTING_ANGLE = ((PMT_SUPPORT_RADIUS - WLS_SUPPORT_RINGSEC_OUTER_RADIUS) / (2 * TOP_SIDE_REFLECTOR_HALF_HEIGHT));
  TOP_SIDE_REFLECTOR_INNER_RADIUS_UPPER_EDGE_IN_GARCOL = PMT_SUPPORT_RADIUS;
  TOP_SIDE_REFLECTOR_OUTER_RADIUS_UPPER_EDGE_IN_GARCOL = (TOP_SIDE_REFLECTOR_INNER_RADIUS_UPPER_EDGE_IN_GARCOL + SIDE_REFLECTOR_THICKNESS);
  TOP_SIDE_REFLECTOR_INNER_RADIUS_LOWER_EDGE_IN_GARCOL = (TOP_SIDE_REFLECTOR_INNER_RADIUS_UPPER_EDGE_IN_GARCOL - 2 * TOP_SIDE_REFLECTOR_HALF_HEIGHT_IN_GARCOL * TOP_SIDE_REFLECTOR_COTAN_TILTING_ANGLE);
  TOP_SIDE_REFLECTOR_OUTER_RADIUS_LOWER_EDGE_IN_GARCOL = (TOP_SIDE_REFLECTOR_INNER_RADIUS_LOWER_EDGE_IN_GARCOL + SIDE_REFLECTOR_THICKNESS);
  TOP_SIDE_REFLECTOR_IN_GARCOL_POS_Z = (TOP_PMT_SUPPORT_POS_Z - PMT_SUPPORT_HALF_HEIGHT - TOP_SIDE_REFLECTOR_HALF_HEIGHT_IN_GARCOL);
  TOP_SIDE_REFLECTOR_IN_GARCOL_POS_Z_GARCOL = (TOP_SIDE_REFLECTOR_IN_GARCOL_POS_Z - GAR_COLUMN_Z);
  TOP_SIDE_REFLECTOR_INNER_RADIUS_UPPER_EDGE_IN_LARCOL = TOP_SIDE_REFLECTOR_INNER_RADIUS_LOWER_EDGE_IN_GARCOL;
  TOP_SIDE_REFLECTOR_OUTER_RADIUS_UPPER_EDGE_IN_LARCOL = (TOP_SIDE_REFLECTOR_INNER_RADIUS_UPPER_EDGE_IN_LARCOL + SIDE_REFLECTOR_THICKNESS);
  TOP_SIDE_REFLECTOR_INNER_RADIUS_LOWER_EDGE_IN_LARCOL = WLS_SUPPORT_RINGSEC_INNER_RADIUS;
  TOP_SIDE_REFLECTOR_OUTER_RADIUS_LOWER_EDGE_IN_LARCOL = (TOP_SIDE_REFLECTOR_INNER_RADIUS_LOWER_EDGE_IN_LARCOL + SIDE_REFLECTOR_THICKNESS);
  TOP_SIDE_REFLECTOR_IN_LARCOL_POS_Z = (TOP_SIDE_REFLECTOR_IN_GARCOL_POS_Z - TOP_SIDE_REFLECTOR_HALF_HEIGHT);
  TOP_SIDE_REFLECTOR_IN_LARCOL_POS_Z_LARCOL = (TOP_SIDE_REFLECTOR_IN_LARCOL_POS_Z - LARCOL_CYLINDER_POS_Z);
  SIDE_REFLECTOR_COAT_THICKNESS = (0.001 * mm);
  BTM_SIDE_REFLECTOR_COAT_OUTER_RADIUS_UPPER_EDGE = BTM_SIDE_REFLECTOR_INNER_RADIUS_UPPER_EDGE;
  BTM_SIDE_REFLECTOR_COAT_INNER_RADIUS_UPPER_EDGE = (BTM_SIDE_REFLECTOR_COAT_OUTER_RADIUS_UPPER_EDGE - SIDE_REFLECTOR_COAT_THICKNESS);
  BTM_SIDE_REFLECTOR_COAT_OUTER_RADIUS_LOWER_EDGE = BTM_SIDE_REFLECTOR_INNER_RADIUS_LOWER_EDGE;
  BTM_SIDE_REFLECTOR_COAT_INNER_RADIUS_LOWER_EDGE = (BTM_SIDE_REFLECTOR_COAT_OUTER_RADIUS_LOWER_EDGE - SIDE_REFLECTOR_COAT_THICKNESS);
  BTM_SIDE_REFLECTOR_COAT_HALF_HEIGHT_GASTEST = BTM_SIDE_REFLECTOR_HALF_HEIGHT_GASTEST;
  BTM_SIDE_REFLECTOR_COAT_POS_Z_GASTEST = BTM_SIDE_REFLECTOR_POS_Z_GASTEST;
  BTM_SIDE_REFLECTOR_COAT_POS_Z_GASTEST_LARCOL = (BTM_SIDE_REFLECTOR_COAT_POS_Z_GASTEST - LARCOL_CYLINDER_POS_Z);
  BTM_SIDE_REFLECTOR_COAT_START_PHI = BTM_SIDE_REFLECTOR_START_PHI;
  BTM_SIDE_REFLECTOR_COAT_DELTA_PHI = BTM_SIDE_REFLECTOR_DELTA_PHI;
  BTM_SIDE_REFLECTOR_COAT_CYNLIDRICAL_PART_OUTER_RADIUS = BTM_SIDE_REFLECTOR_CYNLIDRICAL_PART_INNER_RADIUS;
  BTM_SIDE_REFLECTOR_COAT_CYNLIDRICAL_PART_INNER_RADIUS = (BTM_SIDE_REFLECTOR_COAT_CYNLIDRICAL_PART_OUTER_RADIUS - SIDE_REFLECTOR_COAT_THICKNESS);
  BTM_SIDE_REFLECTOR_COAT_CYNLIDRICAL_PART_HALF_HEIGHT = BTM_SIDE_REFLECTOR_CYNLIDRICAL_PART_HALF_HEIGHT;
  BTM_SIDE_REFLECTOR_COAT_CYNLIDRICAL_PART_POS_Z = BTM_SIDE_REFLECTOR_CYNLIDRICAL_PART_POS_Z;
  BTM_SIDE_REFLECTOR_COAT_CYNLIDRICAL_PART_POS_Z_LARCOL = (BTM_SIDE_REFLECTOR_COAT_CYNLIDRICAL_PART_POS_Z - LARCOL_CYLINDER_POS_Z);
  TOP_SIDE_REFLECTOR_COAT_OUTER_RADIUS_UPPER_EDGE_IN_GARCOL = TOP_SIDE_REFLECTOR_INNER_RADIUS_UPPER_EDGE_IN_GARCOL;
  TOP_SIDE_REFLECTOR_COAT_INNER_RADIUS_UPPER_EDGE_IN_GARCOL = (TOP_SIDE_REFLECTOR_COAT_OUTER_RADIUS_UPPER_EDGE_IN_GARCOL - SIDE_REFLECTOR_COAT_THICKNESS);
  TOP_SIDE_REFLECTOR_COAT_OUTER_RADIUS_LOWER_EDGE_IN_GARCOL = TOP_SIDE_REFLECTOR_INNER_RADIUS_LOWER_EDGE_IN_GARCOL;
  TOP_SIDE_REFLECTOR_COAT_INNER_RADIUS_LOWER_EDGE_IN_GARCOL = (TOP_SIDE_REFLECTOR_COAT_OUTER_RADIUS_LOWER_EDGE_IN_GARCOL - SIDE_REFLECTOR_COAT_THICKNESS);
  TOP_SIDE_REFLECTOR_COAT_HALF_HEIGHT_IN_GARCOL = TOP_SIDE_REFLECTOR_HALF_HEIGHT_IN_GARCOL;
  TOP_SIDE_REFLECTOR_COAT_IN_GARCOL_POS_Z = TOP_SIDE_REFLECTOR_IN_GARCOL_POS_Z;
  TOP_SIDE_REFLECTOR_COAT_IN_GARCOL_POS_Z_GARCOL = (TOP_SIDE_REFLECTOR_COAT_IN_GARCOL_POS_Z - GAR_COLUMN_Z);
  TOP_SIDE_REFLECTOR_COAT_OUTER_RADIUS_UPPER_EDGE_IN_LARCOL = TOP_SIDE_REFLECTOR_INNER_RADIUS_UPPER_EDGE_IN_LARCOL;
  TOP_SIDE_REFLECTOR_COAT_INNER_RADIUS_UPPER_EDGE_IN_LARCOL = (TOP_SIDE_REFLECTOR_COAT_OUTER_RADIUS_UPPER_EDGE_IN_LARCOL - SIDE_REFLECTOR_COAT_THICKNESS);
  TOP_SIDE_REFLECTOR_COAT_OUTER_RADIUS_LOWER_EDGE_IN_LARCOL = TOP_SIDE_REFLECTOR_INNER_RADIUS_LOWER_EDGE_IN_LARCOL;
  TOP_SIDE_REFLECTOR_COAT_INNER_RADIUS_LOWER_EDGE_IN_LARCOL = (TOP_SIDE_REFLECTOR_COAT_OUTER_RADIUS_LOWER_EDGE_IN_LARCOL - SIDE_REFLECTOR_COAT_THICKNESS);
  TOP_SIDE_REFLECTOR_COAT_HALF_HEIGHT_IN_LARCOL = TOP_SIDE_REFLECTOR_HALF_HEIGHT_IN_LARCOL;
  TOP_SIDE_REFLECTOR_COAT_IN_LARCOL_POS_Z = TOP_SIDE_REFLECTOR_IN_LARCOL_POS_Z;
  TOP_SIDE_REFLECTOR_COAT_IN_LARCOL_POS_Z_LARCOL = (TOP_SIDE_REFLECTOR_COAT_IN_LARCOL_POS_Z - LARCOL_CYLINDER_POS_Z);
  TOP_SIDE_REFLECTOR_COAT_START_PHI = TOP_SIDE_REFLECTOR_START_PHI;
  TOP_SIDE_REFLECTOR_COAT_DELTA_PHI = TOP_SIDE_REFLECTOR_DELTA_PHI;
  FIELD_SHAPER_RING_RINGSEC_CROSSSECTION_OUTER_RADIUS = (3. * mm);
  FIELD_SHAPER_RING_RINGSEC_CROSSSECTION_THICKNESS = (1. * mm);
  FIELD_SHAPER_RING_RINGSEC_CROSSSECTION_INNER_RADIUS = (FIELD_SHAPER_RING_RINGSEC_CROSSSECTION_OUTER_RADIUS - FIELD_SHAPER_RING_RINGSEC_CROSSSECTION_THICKNESS);
  FIELD_SHAPER_RING_RINGSEC_TORUS_RADIUS = (400. * mm);
  FIELD_SHAPER_RING_RINGSEC_TORUS_START_PHI = WLS_RINGSEC_START_PHI;
  FIELD_SHAPER_RING_RINGSEC_TORUS_DELTA_PHI = WLS_RINGSEC_DELTA_PHI;
  FIELD_SHAPER_LINSEC_INNER_RADIUS = FIELD_SHAPER_RING_RINGSEC_CROSSSECTION_INNER_RADIUS;
  FIELD_SHAPER_LINSEC_OUTER_RADIUS = FIELD_SHAPER_RING_RINGSEC_CROSSSECTION_OUTER_RADIUS;
  FIELD_SHAPER_LINSEC_HALF_LENGTH = ((FIELD_SHAPER_RING_RINGSEC_TORUS_RADIUS + FIELD_SHAPER_RING_RINGSEC_CROSSSECTION_OUTER_RADIUS) * sin(FIELD_SHAPER_RING_RINGSEC_TORUS_START_PHI));
  FIELD_SHAPER_PILLAR_INNER_RADIUS = (0. * mm);
  FIELD_SHAPER_PILLAR_OUTER_RADIUS = (20. * mm);
  FIELD_SHAPER_PILLAR_HALF_HEIGHT = (1177. / 2 * mm);
  FIELD_SHAPER_PILLAR_POS_Z = (CATHODE_POS_Z + CATHODE_PLATE_HALF_THICKNESS + FIELD_SHAPER_PILLAR_HALF_HEIGHT);
  FIELD_SHAPER_PILLAR_POS_Z_LAR = (FIELD_SHAPER_PILLAR_POS_Z - LARCOL_CYLINDER_POS_Z);
  DISTANCE_FIELD_SHAPER_PILLAR_CENTER_TO_DETECTOR_CENTER = (FIELD_SHAPER_RING_RINGSEC_TORUS_RADIUS + FIELD_SHAPER_RING_RINGSEC_CROSSSECTION_OUTER_RADIUS + FIELD_SHAPER_PILLAR_OUTER_RADIUS);
  FIELD_SHAPER_RING_PITCH = (40. * mm);
  DISTANCE_PILLAR_UPPER_EDGE_TO_CENTER_OF_FIRST_FIELD_SHAPER_RING = (103. * mm);
  FIRST_FIELD_SHAPER_RING_POS_Z = (FIELD_SHAPER_PILLAR_POS_Z + FIELD_SHAPER_PILLAR_HALF_HEIGHT - DISTANCE_PILLAR_UPPER_EDGE_TO_CENTER_OF_FIRST_FIELD_SHAPER_RING);
  PERLITE_COLUMN_INNER_RADIUS = (0. * mm);
  PERLITE_COLUMN_OUTER_RADIUS = (718. * mm);
  PERLITE_COLUMN_HALF_HEIGHT = (376. * mm);
  DISTANCE_PERLITE_COLUMN_UPPER_EDGE_TO_TOP_FLANGE_LOWER_EDGE = (279.5 * mm);
  PERLITE_COLUMN_UPPER_EDGE_POS_Z = (TOP_FLANGE_POS_Z - TOP_FLANGE_HALF_HEIGHT + DISTANCE_PERLITE_COLUMN_UPPER_EDGE_TO_TOP_FLANGE_LOWER_EDGE);
  PERLITE_COLUMN_POS_Z = (PERLITE_COLUMN_UPPER_EDGE_POS_Z - PERLITE_COLUMN_HALF_HEIGHT);
  NEUTRON_SHIELD_THICKNESS = (90. * mm);
  NEUTRON_SHIELD_DISTANCE_TO_REFLECTOR = (10. * mm);
  NEUTRON_SHIELD_INNER_RADIUS = (WLS_SUPPORT_RINGSEC_OUTER_RADIUS + NEUTRON_SHIELD_DISTANCE_TO_REFLECTOR);
  NEUTRON_SHIELD_OUTER_RADIUS = (NEUTRON_SHIELD_INNER_RADIUS + NEUTRON_SHIELD_THICKNESS);
  NEUTRON_SHIELD_HALF_HEIGHT = (WLS_SUPPORT_RINGSEC_HALF_HEIGHT + BTM_SIDE_REFLECTOR_HALF_HEIGHT);
  NEUTRON_SHIELD_POS_Z = (WLS_SUPPORT_RINGSEC_POS_Z - BTM_SIDE_REFLECTOR_HALF_HEIGHT);
  NEUTRON_SHIELD_POS_Z_LAR = (NEUTRON_SHIELD_POS_Z - LARCOL_CYLINDER_POS_Z);

  PE_SHIELD_WIDTH = (600. * mm);
  PE_SHIELD_INNER_RADIUS = TOP_FLANGE_OUTER_RADIUS + (45. * mm);
  PE_SHIELD_OUTER_RADIUS = PE_SHIELD_INNER_RADIUS + PE_SHIELD_WIDTH;
  PE_SHIELD_HALF_HEIGHT = (4470. * mm / 2.);
  PE_SHIELD_POS_Z = 0.;

  PE_SHIELD_LEAD_WIDTH = (50. * mm);
  PE_SHIELD_LEAD_INNER_RADIUS = TANK_CYLINDER_OUTER_RADIUS;
  PE_SHIELD_LEAD_OUTER_RADIUS = PE_SHIELD_LEAD_INNER_RADIUS + PE_SHIELD_LEAD_WIDTH;
  PE_SHIELD_LEAD_POS_Z = 2093. * mm / 2;
}
