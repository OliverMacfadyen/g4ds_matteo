#include "DSDetectorARISERGamma.hh"

#include <assert.h>

#include <iostream>

#include "DSDetectorARISER.hh"
#include "DSIO.hh"
#include "DSLogger.hh"
#include "DSMaterial.hh"
#include "DSParameters.hh"
#include "DSStorage.hh"
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalConstants.hh"
#include "G4SubtractionSolid.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4VisAttributes.hh"

using namespace std;

/**
 * Description.
 *
 * The Na-22 source is located along the +x axis, aligned with the z-center of the TPC.
 * The BaF2 detector which detects the back-to-back 511 keV gamma is on the +x side of the source holder.
 * The BaF2 detector which detects the isotropic 1270 keV gamma is on the +z side of the source holder.
 */

DSDetectorARISERGamma::DSDetectorARISERGamma(G4VPhysicalVolume *myMotherVolume) {
  fMotherVolume = myMotherVolume;
  DSLog(routine) << "Constructing ARIS-ER Gamma Detector Geometries" << endlog;
  G4bool myCheckOverlap = DSStorage::Get()->GetCheckOverlap();

  //--------------------------------------------------------------
  // Positioning
  //--------------------------------------------------------------
  // Na22 source and holder
  const G4double sourceDistanceToCryo = DSStorage::Get()->GetAriserSourceDistanceToCryo();

  // BaF2 detectors
  const G4double bafDistanceToSource511 = 3 * cm;
  const G4double bafDistanceToSource1270 = bafDistanceToSource511;

  // Germanium detector
  const G4double begeDistanceToCryo = DSStorage::Get()->GetAriserBEGeDistanceToCryo();
  const G4double begeAngle = DSStorage::Get()->GetAriserBEGeAngle();

  DSLog(routine) << "Na22 source distance to cryostat: " << sourceDistanceToCryo / mm << "(mm)" << endlog;
  DSLog(routine) << "BaF2 detector (511 keV gamma) distance to source: " << bafDistanceToSource511 / mm << "(mm)" << endlog;
  DSLog(routine) << "BaF2 detector (1270 keV gamma) distance to source: " << bafDistanceToSource1270 / mm << "(mm)" << endlog;
  DSLog(routine) << "BEGe detector distance to cryostat: " << begeDistanceToCryo / mm << "(mm)" << endlog;
  DSLog(routine) << "BEGe detector angle: " << begeAngle / deg << "(deg)" << endlog;

  //--------------------------------------------------------------
  // Fixed properties
  //--------------------------------------------------------------
  // TPC parameters
  const G4double cryoOuterR = 6.10 * 2.54 / 2. * cm;
  const G4ThreeVector LArTPCShift(0, 0, -GetCm(0.26) * cm);

  // Na22 source and holder
  const G4double sourceActiveLength = 1.5 * mm;
  const G4double sourceActiveDiameter = 2.1 * mm;
  const G4double sourceHolderLength = 3.0 * mm;
  const G4double sourceHolderDiameter = 3.0 * mm;
  const G4double sourceHolderWindowThickness = 0.225 * mm;
  const G4ThreeVector sourceHCeramicPosition = G4ThreeVector(0., 0., 0.5 * sourceHolderLength - 0.5 * sourceActiveLength - sourceHolderWindowThickness);

  // BaF2 detectors
  const G4double bafDetectorHeight = 30. / 2. * mm;
  const G4double bafDetectorRadiusMin = 25. / 2. * mm;
  const G4double bafDetectorRadiusMax = 38. / 2. * mm;
  const G4double bafShieldHeight = 32. / 2. * mm;
  const G4double bafShieldRadiusMin = 27. / 2. * mm;
  const G4double bafShieldRadiusMax = 40.87 / 2. * mm;

  // Germanium detector
  const G4double begeCrystalRadius = 51.3 * mm / 2.;
  const G4double begeCrystalThickness = 21.7 * mm;
  const G4double begeVacuumThickness = 4.7 * mm;
  const G4double begeShieldThickness = 0.1 * mm;
  const G4double begeWindowThickness = 0.6 * mm;

  //--------------------------------------------------------------
  // Na22 source and holder
  //--------------------------------------------------------------
  DSLog(routine) << "Building Na-22 source and holder." << endlog;
  G4ThreeVector sourcePosition(sourceDistanceToCryo + cryoOuterR, 0., 0.);
  sourcePosition += LArTPCShift;
  DSStorage::Get()->SetAriserSourcePosition(sourcePosition);
  DSLog(routine) << "Na-22 source position"
                 << ": X " << sourcePosition.getX() << ", Y " << sourcePosition.getY() << ", Z " << sourcePosition.getZ() << endlog;
  G4RotationMatrix *sourceRotation = new G4RotationMatrix();
  sourceRotation->rotateY(0.5 * M_PI);

  G4Tubs *fSolidSourceHolder = new G4Tubs("SourceHolder_Solid", 0, sourceHolderDiameter / 2., sourceHolderLength / 2., 0, 2 * M_PI);
  G4LogicalVolume *fLogicSourceHolder = new G4LogicalVolume(fSolidSourceHolder, DSMaterial::Get()->GetStainlessSteel(), "SourceHolder_Logic");
  fPhysicSourceHolder = new G4PVPlacement(sourceRotation, sourcePosition, "SourceHolder", fLogicSourceHolder, fMotherVolume, false, 0, myCheckOverlap);

  G4Tubs *fSolidSourceActive = new G4Tubs("SourceActive_Solid", 0, sourceActiveDiameter / 2., sourceActiveLength / 2., 0, 2 * M_PI);
  G4LogicalVolume *fLogicSourceActive = new G4LogicalVolume(fSolidSourceActive, DSMaterial::Get()->GetCeramic(), "SourceActive_Logic");
  fPhysicSourceActive = new G4PVPlacement(0, sourceHCeramicPosition, "SourceActive", fLogicSourceActive, fPhysicSourceHolder, false, 0, myCheckOverlap);

  //--------------------------------------------------------------
  // BaF2 gamma detectors
  //--------------------------------------------------------------
  DSLog(routine) << "Building BaF2 gamma detectors." << endlog;
  G4RotationMatrix *bafBaseRotation = new G4RotationMatrix();
  bafBaseRotation->rotateY(M_PI);  // Account for direction of conical shape

  // Detector for 511 keV gamma
  G4ThreeVector bafPosition511 = G4ThreeVector(bafDistanceToSource511 + 0.5 * bafShieldHeight, 0., 0.);
  bafPosition511 += sourcePosition;
  G4RotationMatrix *bafRotation511 = new G4RotationMatrix();  // Make sure BaF2 is facing the source
  bafRotation511->rotateZ(-(bafPosition511 - sourcePosition).phi());
  bafRotation511->rotateY(-(bafPosition511 - sourcePosition).theta());
  *bafRotation511 = (*bafBaseRotation) * (*bafRotation511);
  DSLog(routine) << "BaF Detector for 511 keV gamma"
                 << ": X " << bafPosition511.getX() << ", Y " << bafPosition511.getY() << ", Z " << bafPosition511.getZ() << endlog;

  // Detector for 1270 keV gamma
  G4ThreeVector bafPosition1270 = G4ThreeVector(0., 0., bafDistanceToSource1270 + 0.5 * bafShieldHeight);
  bafPosition1270 += sourcePosition;
  G4RotationMatrix *bafRotation1270 = new G4RotationMatrix();  // Make sure BaF2 is facing the source
  bafRotation1270->rotateZ(-(bafPosition1270 - sourcePosition).phi());
  bafRotation1270->rotateY(-(bafPosition1270 - sourcePosition).theta());
  *bafRotation1270 = (*bafBaseRotation) * (*bafRotation1270);
  DSLog(routine) << "BaF Detector for 1270 keV gamma"
                 << ": X " << bafPosition1270.getX() << ", Y " << bafPosition1270.getY() << ", Z " << bafPosition1270.getZ() << endlog;

  G4Cons *fSolidBafShield = new G4Cons("BafShield_solid", 0, bafShieldRadiusMax, 0, bafShieldRadiusMin, bafShieldHeight, 0, 2 * pi);
  G4LogicalVolume *fLogicBafShield511 = new G4LogicalVolume(fSolidBafShield, DSMaterial::Get()->GetAluminum(), "BafShield511_logic");
  G4LogicalVolume *fLogicBafShield1270 = new G4LogicalVolume(fSolidBafShield, DSMaterial::Get()->GetAluminum(), "BafShield1270_logic");
  fPhysicBafShield511 = new G4PVPlacement(bafRotation511, bafPosition511, "BafShield511", fLogicBafShield511, fMotherVolume, false, 0, myCheckOverlap);
  fPhysicBafShield1270 = new G4PVPlacement(bafRotation1270, bafPosition1270, "BafShield1270", fLogicBafShield1270, fMotherVolume, false, 0, myCheckOverlap);

  G4Cons *fSolidBafDetector = new G4Cons("BafDetector_solid", 0, bafDetectorRadiusMax, 0, bafDetectorRadiusMin, bafDetectorHeight, 0, 2 * pi);
  G4LogicalVolume *fLogicBafDetector511 = new G4LogicalVolume(fSolidBafDetector, DSMaterial::Get()->GetBariumFluoride1(), "BafDetector511_logic");
  G4LogicalVolume *fLogicBafDetector1270 = new G4LogicalVolume(fSolidBafDetector, DSMaterial::Get()->GetBariumFluoride2(), "BafDetector1270_logic");
  fPhysicBafDetector511 = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), fLogicBafDetector511, "BafDetector511", fLogicBafShield511, false, 0, myCheckOverlap);
  fPhysicBafDetector1270 = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), fLogicBafDetector1270, "BafDetector1270", fLogicBafShield1270, false, 0, myCheckOverlap);

  //--------------------------------------------------
  // Germanium detector (for Compton scattered gamma)
  //--------------------------------------------------
  DSLog(routine) << "Building germanium gamma detector." << endlog;
  G4double begeDetectorHeight = begeWindowThickness + begeVacuumThickness + begeCrystalThickness;
  G4double begeDetectorRadius = begeShieldThickness + begeVacuumThickness + begeCrystalRadius;

  G4ThreeVector begePosition = G4ThreeVector(-cryoOuterR - begeDistanceToCryo  - 0.5 * begeDetectorHeight, 0., 0.);
  begePosition.rotateZ(-begeAngle);
  G4RotationMatrix *begeRotation = new G4RotationMatrix();  // Make sure BEGe detector is facing the TPC
  // begeRotation->rotateZ(begePosition.phi());
  // begeRotation->rotateY(begePosition.theta());
  begeRotation->rotateY(-M_PI / 2.);  // Face towards the TPC
  begeRotation->rotateX(-begeAngle);
  begePosition += LArTPCShift;  // Now add the vertical offset after using the position to calculate the rotation
  DSLog(routine) << "BEGe detector" << ": X " << begePosition.getX() << ", Y " << begePosition.getY() << ", Z " << begePosition.getZ() << endlog;

  G4Tubs *fSolidBegeShield = new G4Tubs("BegeShield_solid", 0., begeDetectorRadius, begeDetectorHeight / 2., 0., 2 * M_PI);
  G4LogicalVolume *fLogicBegeShield = new G4LogicalVolume(fSolidBegeShield, DSMaterial::Get()->GetStainlessSteel(), "BegeShield_logic");
  fPhysicBegeShield = new G4PVPlacement(begeRotation, begePosition, "BegeShield", fLogicBegeShield, fMotherVolume, false, 0, myCheckOverlap);

  G4Tubs *fSolidBegeWindow = new G4Tubs("BegeWindow_solid", 0., begeDetectorRadius - begeShieldThickness, begeDetectorHeight / 2., 0., 2 * M_PI);
  G4LogicalVolume *fLogicBegeWindow = new G4LogicalVolume(fSolidBegeWindow, DSMaterial::Get()->GetCarbon(), "BegeWindow_logic");
  fPhysicBegeWindow = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), "BegeWindow", fLogicBegeWindow, fPhysicBegeShield, false, 0, myCheckOverlap);

  G4Tubs *fSolidBegeVacuum = new G4Tubs("BegeVacuum_solid", 0., begeDetectorRadius - begeShieldThickness, (begeDetectorHeight - begeWindowThickness) / 2., 0., 2 * M_PI);
  G4LogicalVolume *fLogicBegeVacuum = new G4LogicalVolume(fSolidBegeVacuum, DSMaterial::Get()->GetVacuum(), "BegeVacuum_logic");
  fPhysicBegeVacuum = new G4PVPlacement(0, G4ThreeVector(0, 0, -begeWindowThickness / 2.), "BegeVacuum", fLogicBegeVacuum, fPhysicBegeWindow, false, 0., myCheckOverlap);

  G4Tubs *fSolidBegeCrystal = new G4Tubs("BegeCrystal_solid", 0., begeCrystalRadius, begeCrystalThickness / 2., 0., 2 * M_PI);
  G4LogicalVolume *fLogicBegeCrystal = new G4LogicalVolume(fSolidBegeCrystal, DSMaterial::Get()->GetGermanium(), "BegeCrystal_Logic");
  fPhysicBegeCrystal = new G4PVPlacement(0, G4ThreeVector(0., 0., -begeVacuumThickness / 2.), "BegeCrystal", fLogicBegeCrystal, fPhysicBegeVacuum, false, 0., myCheckOverlap);

  DefineSurfaces();

  //--------------------------------------------------
  // Configure visualization
  //--------------------------------------------------
  G4Color colGray = G4Color(0.5, 0.5, 0.5, 0.5);
  G4Color colTransparent = G4Color(0, 0, 0, 0);
  fLogicSourceHolder->SetVisAttributes(colGray);
  fLogicSourceActive->SetVisAttributes(G4Color(1, 0, 0, 0.5));
  fLogicBafShield511->SetVisAttributes(colGray);
  fLogicBafShield1270->SetVisAttributes(colGray);
  fLogicBafDetector511->SetVisAttributes(G4Color(0, 0, 1, 0.5));
  fLogicBafDetector1270->SetVisAttributes(G4Color(0, 1, 1, 0.5));
  fLogicBegeShield->SetVisAttributes(colGray);
  fLogicBegeWindow->SetVisAttributes(colGray);
  fLogicBegeVacuum->SetVisAttributes(colTransparent);
  fLogicBegeCrystal->SetVisAttributes(G4Color(1, 0, 1, 0.5));
}

DSDetectorARISERGamma::~DSDetectorARISERGamma() {
  ;  // delete fMessenger;
}

void DSDetectorARISERGamma::DefineSurfaces() { return; }
