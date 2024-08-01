//////////////////////////////////////////
//                     _                //
//    __ _  __ _ _ __ | |_ _ __ _   _   //
//   / _` |/ _` | '_ \| __| '__| | | |  //
//  | (_| | (_| | | | | |_| |  | |_| |  //
//   \__, |\__,_|_| |_|\__|_|   \__, |  //
//   |___/                      |___/   //
//                                      //
//////////////////////////////////////////
// Gantry geometry devise for calibration
// dadoun@in2p3.fr, last update October 2019
// Re-written by AK 200701 (akish@cern.ch) to have automatic build with only few
// input dimensions

#include "DSDetectorGantryCalibration.hh"
#include "DSMaterial.hh"
#include "DSParameters.hh"
#include "DSStorage.hh"
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4Sphere.hh"
#include "G4ThreeVector.hh"
#include "G4Torus.hh"
#include "G4Transform3D.hh"
#include "G4Tubs.hh"
#include "G4UnionSolid.hh"
#include "G4VisAttributes.hh"

using namespace std;

G4Colour myRed(1.0, 0.0, 0.0);
G4Colour myGreen(0.0, 1.0, 0.0);
G4Colour myBlue(0.0, .0, 1.0);
G4Colour myYellow(1.0, 1.0, 0.0);
G4Colour myMagenta(0.0, 1.0, 1.0);
G4Colour myWhite(1.0, 1.0, 1.0);
G4Colour myGrey(0.5, 0.5, 0.5);

DSDetectorGantryCalibration::DSDetectorGantryCalibration(G4VPhysicalVolume* myMotherVolume) {
  fMotherVolume = myMotherVolume;

  GantryConfiguration = DSStorage::Get()->GetGantryConfiguration();
  if (GantryConfiguration == 1) {
    G4cout << "GantryConfiguration = " << GantryConfiguration << G4endl;
    GantrySurface = DSStorage::Get()->GetGantrySurface();
    GantryShield = DSStorage::Get()->GetGantryShield();

    // build gantry pipes
    ConstructGantry();

    // build optical surfaces for the pipes
    G4cout << "GantrySurface = " << GantrySurface << G4endl;
    if (GantrySurface == 1) { DefineGantrySurfaces(); }

    // build lead shield for YBe source
    G4cout << "GantryShield = " << GantryShield << G4endl;
    if (GantryShield == 1) { ConstructGantryShield(); }
  }
}

void DSDetectorGantryCalibration::ConstructGantry() {

  matSLS = DSMaterial::Get()->GetStainlessSteel();
  matNitrogen = DSMaterial::Get()->GetGaseousNitrogen();
  TwoPi = 2 * M_PI * rad;
  DefaultAngle = TwoPi / 16.;  // rotation due to octogone stuff, see DSDetectorPlasticVeto
  fCheckOverlaps = DSStorage::Get()->GetCheckOverlap();
  half_tpc_edge = DSStorage::Get()->GetDS20kTPCedge() / 2;
  VetoBufferThickness = DSStorage::Get()->GetDS20kLArBufferThickness();
  TPCWallThickness = DSStorage::Get()->GetDS20kTPCVesselThickness();
  TPBThickness = DSStorage::Get()->GetDS20kTPBThickness();
  ElectronicsThickness = DSStorage::Get()->GetDS20kLArThicknessAboveTPC();  // this is not below the TPC, but
                                                                            // has the same dimension

  Rinner = half_tpc_edge / tan(DefaultAngle);
  Router = half_tpc_edge / sin(DefaultAngle);

  PipeID = DSStorage::Get()->GetGantryPipeID();
  PipeWall = DSStorage::Get()->GetGantryPipeWall();
  PipeOD = PipeID + 2 * PipeWall;
  PipeIR = PipeID / 2;
  PipeOR = PipeID / 2 + PipeWall;
  PipeBendingR = DSStorage::Get()->GetGantryPipeBendingR();
  DistanceTPCside = DSStorage::Get()->GetGantryDistanceTPCside();
  DistanceTPCbottom = DSStorage::Get()->GetGantryDistanceTPCbottom();
  DistanceBtwPipes = DSStorage::Get()->GetGantryDistanceBtwPipes();

  // if the shield to be built, check dimensions
  if (GantryShield == 1) {
    GantryShieldThickness = DSStorage::Get()->GetGantryShieldThickness();
    GantryShieldHeight = DSStorage::Get()->GetGantryShieldHeight();
    GantryShieldOffset = DSStorage::Get()->GetGantryShieldOffset();
    GantryShieldMaterial = DSStorage::Get()->GetGantryShieldMaterial();
  }

  Rinner_ext = Rinner + TPCWallThickness + TPBThickness;
  Router_ext = Rinner_ext / cos(DefaultAngle);
  half_tpc_edge_ext = Router_ext * tan(DefaultAngle);

  TPCheight = DSStorage::Get()->GetDS20kTPCheight();
  Pipe1HeightTotal = TPCheight + 2 * TPCWallThickness + 2 * ElectronicsThickness + DistanceTPCbottom + VetoBufferThickness * 0.8;  // the only adhoc scaling factor to not hit the top veto panel
  Pipe2HeightTotal = Pipe1HeightTotal + DistanceBtwPipes + PipeOR;
  if (GantryShield == 1) {  // 2nd pipe must be longer as it goes lower through the shield
    Pipe2HeightTotal += GantryShieldHeight + GantryShieldOffset;
  }
  Pipe1VertHeight = Pipe1HeightTotal - PipeBendingR;
  Pipe2VertHeight = Pipe2HeightTotal - PipeBendingR;

  PipeVertR = Rinner_ext + DistanceTPCside + PipeOR;
  PipeHorLength = (PipeVertR - PipeBendingR) * 2;

  G4double PipeBendAngle = 90.;                                     // degree
  G4double PipeBendAngle_rad = (M_PI / 180) * PipeBendAngle * rad;  // in rad

  // Pipe 1
  G4Tubs* Pipe1Vert = new G4Tubs("GantryPipe1Vert", 0, PipeOR, Pipe1VertHeight / 2, 0., 2 * M_PI * rad);
  G4Tubs* Pipe2Vert = new G4Tubs("GantryPipe2Vert", 0, PipeOR, Pipe2VertHeight / 2, 0., 2 * M_PI * rad);
  G4Tubs* PipeHor = new G4Tubs("GantryPipeHor", 0, PipeOR, PipeHorLength / 2, 0., 2 * M_PI * rad);
  G4Torus* PipeBend = new G4Torus("GantryPipeBend", 0, PipeOR, PipeBendingR, 0, PipeBendAngle_rad);

  G4Tubs* Pipe1NitroVert = new G4Tubs("GantryPipe1NitroVert", 0, PipeIR, Pipe1VertHeight / 2, 0., 2 * M_PI * rad);
  G4Tubs* Pipe2NitroVert = new G4Tubs("GantryPipe2NitroVert", 0, PipeIR, Pipe2VertHeight / 2, 0., 2 * M_PI * rad);
  G4Tubs* PipeNitroHor = new G4Tubs("GantryPipeNitroHor", 0, PipeIR, PipeHorLength / 2, 0., 2 * M_PI * rad);
  G4Torus* PipeNitroBend = new G4Torus("GantryPipeNitroBend", 0, PipeIR, PipeBendingR, 0, PipeBendAngle_rad);

  G4ThreeVector PosPipe1Bend1(0, 0, 0);
  PosPipe1Bend1 = G4ThreeVector(0, PipeBendingR, -Pipe1VertHeight / 2);
  G4ThreeVector PosPipe2Bend1(0, 0, 0);
  PosPipe2Bend1 = G4ThreeVector(0, PipeBendingR, -Pipe2VertHeight / 2);
  G4RotationMatrix* RotationBend1 = new G4RotationMatrix();
  RotationBend1->rotateX(M_PI);
  RotationBend1->rotateY(M_PI / 2);
  G4Transform3D* transformPipe1Bend1 = new G4Transform3D(*RotationBend1, PosPipe1Bend1);
  G4Transform3D* transformPipe2Bend1 = new G4Transform3D(*RotationBend1, PosPipe2Bend1);

  G4ThreeVector PosPipe1Hor(0, 0, 0);
  PosPipe1Hor = G4ThreeVector(0, PipeHorLength / 2 + PipeBendingR, -Pipe1VertHeight / 2 - PipeBendingR);
  G4ThreeVector PosPipe2Hor(0, 0, 0);
  PosPipe2Hor = G4ThreeVector(0, PipeHorLength / 2 + PipeBendingR, -Pipe2VertHeight / 2 - PipeBendingR);
  G4RotationMatrix* RotationHor = new G4RotationMatrix();
  RotationHor->rotateX(M_PI / 2);
  G4Transform3D* transformPipe1Hor = new G4Transform3D(*RotationHor, PosPipe1Hor);
  G4Transform3D* transformPipe2Hor = new G4Transform3D(*RotationHor, PosPipe2Hor);

  G4ThreeVector PosPipe1Bend2(0, 0, 0);
  PosPipe1Bend2 = G4ThreeVector(0, PipeHorLength + PipeBendingR, -Pipe1VertHeight / 2);
  G4ThreeVector PosPipe2Bend2(0, 0, 0);
  PosPipe2Bend2 = G4ThreeVector(0, PipeHorLength + PipeBendingR, -Pipe2VertHeight / 2);
  G4RotationMatrix* RotationBend2 = new G4RotationMatrix();
  RotationBend2->rotateX(M_PI);
  RotationBend2->rotateY(M_PI / 2);
  RotationBend2->rotateZ(M_PI);
  G4Transform3D* transformPipe1Bend2 = new G4Transform3D(*RotationBend2, PosPipe1Bend2);
  G4Transform3D* transformPipe2Bend2 = new G4Transform3D(*RotationBend2, PosPipe2Bend2);

  G4ThreeVector PosPipeVert(0, 0, 0);
  PosPipeVert = G4ThreeVector(0, PipeHorLength + 2 * PipeBendingR, 0);
  G4RotationMatrix* RotationVert = new G4RotationMatrix();  // no rotation
  G4Transform3D* transformPipeVert = new G4Transform3D(*RotationVert, PosPipeVert);

  G4UnionSolid* Pipe1Full = new G4UnionSolid("Pipe1Full", Pipe1Vert, PipeBend, *transformPipe1Bend1);
  Pipe1Full = new G4UnionSolid("Pipe1Full", Pipe1Full, PipeHor, *transformPipe1Hor);
  Pipe1Full = new G4UnionSolid("Pipe1Full", Pipe1Full, PipeBend, *transformPipe1Bend2);
  Pipe1Full = new G4UnionSolid("Pipe1Full", Pipe1Full, Pipe1Vert, *transformPipeVert);

  G4UnionSolid* Pipe1NitroFull = new G4UnionSolid("Pipe1NitroFull", Pipe1NitroVert, PipeNitroBend, *transformPipe1Bend1);
  Pipe1NitroFull = new G4UnionSolid("Pipe1NitroFull", Pipe1NitroFull, PipeNitroHor, *transformPipe1Hor);
  Pipe1NitroFull = new G4UnionSolid("Pipe1NitroFull", Pipe1NitroFull, PipeNitroBend, *transformPipe1Bend2);
  Pipe1NitroFull = new G4UnionSolid("Pipe1NitroFull", Pipe1NitroFull, Pipe1NitroVert, *transformPipeVert);

  G4UnionSolid* Pipe2Full = new G4UnionSolid("Pipe2Full", Pipe2Vert, PipeBend, *transformPipe2Bend1);
  Pipe2Full = new G4UnionSolid("Pipe2Full", Pipe2Full, PipeHor, *transformPipe2Hor);
  Pipe2Full = new G4UnionSolid("Pipe2Full", Pipe2Full, PipeBend, *transformPipe2Bend2);
  Pipe2Full = new G4UnionSolid("Pipe2Full", Pipe2Full, Pipe2Vert, *transformPipeVert);

  G4UnionSolid* Pipe2NitroFull = new G4UnionSolid("Pipe2NitroFull", Pipe2NitroVert, PipeNitroBend, *transformPipe2Bend1);
  Pipe2NitroFull = new G4UnionSolid("Pipe2NitroFull", Pipe2NitroFull, PipeNitroHor, *transformPipe2Hor);
  Pipe2NitroFull = new G4UnionSolid("Pipe2NitroFull", Pipe2NitroFull, PipeNitroBend, *transformPipe2Bend2);
  Pipe2NitroFull = new G4UnionSolid("Pipe2NitroFull", Pipe2NitroFull, Pipe2NitroVert, *transformPipeVert);

  G4double Pipe1Vert_z = -TPCheight / 2 - TPCWallThickness - ElectronicsThickness - TPBThickness - DistanceTPCbottom - PipeOR + PipeBendingR / 2 + Pipe1HeightTotal / 2;                                  // works for PipeID=5cm, BendingR=20cm
  G4double Pipe2Vert_z = -TPCheight / 2 - TPCWallThickness - ElectronicsThickness - TPBThickness - DistanceTPCbottom - PipeOR + PipeBendingR / 2 + Pipe2HeightTotal / 2 - DistanceBtwPipes - 2 * PipeOR;  // works for PipeID=5cm, BendingR=20cm
  if (GantryShield == 1) {                                                                                                                                                                                // move the 2nd pipe lower to fit the shielding cylinder
    Pipe2Vert_z += -GantryShieldHeight / 2 - GantryShieldOffset;
  }

  G4ThreeVector PosPipe1(0, 0, 0);
  PosPipe1 = G4ThreeVector(-PipeVertR * sin(DefaultAngle), -PipeVertR * cos(DefaultAngle), Pipe1Vert_z);

  G4RotationMatrix* RotationPipe1 = new G4RotationMatrix();
  RotationPipe1->rotateZ(-DefaultAngle);
  G4Transform3D* transformPipe1 = new G4Transform3D(*RotationPipe1, PosPipe1);

  G4ThreeVector PosPipe2(0, 0, 0);
  PosPipe2 = G4ThreeVector(PipeVertR * cos(DefaultAngle), -PipeVertR * sin(DefaultAngle), Pipe2Vert_z);

  G4RotationMatrix* RotationPipe2 = new G4RotationMatrix();
  RotationPipe2->rotateZ(M_PI / 2 - DefaultAngle);
  G4Transform3D* transformPipe2 = new G4Transform3D(*RotationPipe2, PosPipe2);

  Pipe1log = new G4LogicalVolume(Pipe1Full, matSLS, "Pipe1log");
  Pipe2log = new G4LogicalVolume(Pipe2Full, matSLS, "Pipe2log");
  Pipe1phys = new G4PVPlacement(*transformPipe1, "Pipe1phys", Pipe1log, fMotherVolume, false, 0, fCheckOverlaps);
  Pipe2phys = new G4PVPlacement(*transformPipe2, "Pipe2phys", Pipe2log, fMotherVolume, false, 0, fCheckOverlaps);

  PipeNitro1log = new G4LogicalVolume(Pipe1NitroFull, matNitrogen, "PipeNitro1log");
  PipeNitro2log = new G4LogicalVolume(Pipe2NitroFull, matNitrogen, "PipeNitro2log");
  PipeNitro1phys = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), "PipeNitro1phys", PipeNitro1log, Pipe1phys, false, 0, fCheckOverlaps);
  PipeNitro2phys = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), "PipeNitro2phys", PipeNitro2log, Pipe2phys, false, 0, fCheckOverlaps);

  G4VisAttributes* pipe1_vis = new G4VisAttributes(myBlue);
  pipe1_vis->SetVisibility(true);
  pipe1_vis->SetForceSolid(true);
  Pipe1log->SetVisAttributes(pipe1_vis);

  G4VisAttributes* pipe2_vis = new G4VisAttributes(myRed);
  pipe2_vis->SetVisibility(true);
  pipe2_vis->SetForceSolid(true);
  Pipe2log->SetVisAttributes(pipe2_vis);

  G4VisAttributes* nitro_vis = new G4VisAttributes(myWhite);
  nitro_vis->SetVisibility(true);
  nitro_vis->SetForceSolid(true);
  PipeNitro1log->SetVisAttributes(nitro_vis);
  PipeNitro2log->SetVisAttributes(nitro_vis);

  // get mass of the pipes
  G4double mass_pipe1 = 0.;
  mass_pipe1 += Pipe1log->GetMass(false, false) / CLHEP::kg;
  G4double mass_pipe2 = 0.;
  mass_pipe2 += Pipe2log->GetMass(false, false) / CLHEP::kg;

  // verbosity
  G4cout << "Pipe ID " << PipeID << "mm,  wall " << PipeWall << "mm,  from TPC by " << DistanceTPCside << " mm" << endl;
  G4cout << "Radial source position at  " << PipeVertR << " mm" << endl;
  G4cout << "Mass of the 1st pipe: " << mass_pipe1 << " kg" << G4endl;
  G4cout << "Mass of the 2nd pipe: " << mass_pipe2 << " kg" << G4endl;
  G4cout << "Total mass of piping steel: " << mass_pipe1 + mass_pipe2 << " kg" << G4endl;
}

DSDetectorGantryCalibration::~DSDetectorGantryCalibration() {}

void DSDetectorGantryCalibration::DefineGantrySurfaces() {

  G4cout << "Gantry pipes optics ON" << G4endl;

  GantrySurfaceOption = DSStorage::Get()->GetGantrySurfaceOption();

  if (GantrySurfaceOption == 1) {
    G4cout << "Gantry pipes surface: Untreated Stainless Steel" << G4endl;
    G4OpticalSurface* fOpUntreatedStainlessSteelSurface = new G4OpticalSurface("OpUntreatedStainlessSteelSurface");
    fOpUntreatedStainlessSteelSurface->SetType(dielectric_metal);
    fOpUntreatedStainlessSteelSurface->SetModel(unified);
    fOpUntreatedStainlessSteelSurface->SetFinish(groundbackpainted);
    fOpUntreatedStainlessSteelSurface->SetMaterialPropertiesTable(DSMaterial::Get()->GetUntreatedStainlessSteelMPT());

    new G4LogicalBorderSurface("Pipe1Surface", fMotherVolume, Pipe1phys, fOpUntreatedStainlessSteelSurface);
    new G4LogicalBorderSurface("Pipe2Surface", fMotherVolume, Pipe2phys, fOpUntreatedStainlessSteelSurface);
  } else if (GantrySurfaceOption == 2) {
    G4cout << "Gantry pipes surface: Electropolished Stainless Steel" << G4endl;
    G4OpticalSurface* fOpElectropolishedStainlessSteelSurface = new G4OpticalSurface("OpElectropolishedStainlessSteelSurface");
    fOpElectropolishedStainlessSteelSurface->SetType(dielectric_metal);
    fOpElectropolishedStainlessSteelSurface->SetModel(unified);
    fOpElectropolishedStainlessSteelSurface->SetFinish(groundbackpainted);
    fOpElectropolishedStainlessSteelSurface->SetMaterialPropertiesTable(DSMaterial::Get()->GetElectropolishedStainlessSteelMPT());

    new G4LogicalBorderSurface("Pipe1Surface", fMotherVolume, Pipe1phys, fOpElectropolishedStainlessSteelSurface);
    new G4LogicalBorderSurface("Pipe2Surface", fMotherVolume, Pipe2phys, fOpElectropolishedStainlessSteelSurface);
  } else if (GantrySurfaceOption == 3) {
    G4cout << "Gantry pipes surface: TPB-coated Stainless Steel" << G4endl;
    G4double SSTPBENE[4] = {0.1 * eV, 8.0 * eV, 8.3 * eV, 20.0 * eV};
    G4double SSTREFUV = DSParameters::Get()->GetTeflonTPBUVRef();
    G4double SSTREFVIS = DSParameters::Get()->GetTeflonTPBVisRef();
    G4double SSTPBREF[4] = {SSTREFVIS, SSTREFVIS, SSTREFUV, SSTREFUV};

    G4OpticalSurface* fOpSSTPBSurface = new G4OpticalSurface("OpSSTBPSurface");
    fOpSSTPBSurface->SetType(dielectric_metal);
    fOpSSTPBSurface->SetModel(unified);
    fOpSSTPBSurface->SetFinish(polished);

    G4MaterialPropertiesTable* fSSTPBSurfProp = new G4MaterialPropertiesTable();
    fSSTPBSurfProp->AddProperty("REFLECTIVITY", SSTPBENE, SSTPBREF, 4);
    fOpSSTPBSurface->SetMaterialPropertiesTable(fSSTPBSurfProp);

    new G4LogicalBorderSurface("Pipe1Surface", fMotherVolume, Pipe1phys, fOpSSTPBSurface);
    new G4LogicalBorderSurface("Pipe2Surface", fMotherVolume, Pipe2phys, fOpSSTPBSurface);
  } else if (GantrySurfaceOption) {
    G4cout << "ERROR: Pipes Optical Surface not defined." << G4endl;
  }
}

void DSDetectorGantryCalibration::ConstructGantryShield() {

  G4cout << "Building calibration gantry shield" << G4endl;

  G4cout << " GantryShieldOffset = " << GantryShieldOffset << G4endl;

  if (GantryShieldMaterial == 1) {
    matShield = DSMaterial::Get()->GetMetalLead();
  } else if (GantryShieldMaterial == 2) {
    matShield = DSMaterial::Get()->GetMetalTungsten();
  } else if (GantryShieldMaterial == 3) {
    matShield = DSMaterial::Get()->GetMetalTantalum();
  } else if (GantryShieldMaterial) {
    G4cout << "ERROR: Gantry Shield material not defined." << G4endl;
  }

  G4double GantryShieldPipeGap = 1.0 * mm;  // radial gap between the pipe and the shield

  G4Tubs* GantryShieldCyl = new G4Tubs("GantryPipe1Vert", PipeOR + GantryShieldPipeGap, PipeOR + GantryShieldPipeGap + GantryShieldThickness, GantryShieldHeight / 2, 0., 2 * M_PI * rad);

  G4double GantryShield_z = -TPCheight / 2 - TPCWallThickness - ElectronicsThickness - TPBThickness - GantryShieldOffset - GantryShieldHeight / 2;
  G4ThreeVector PosGantryShield(0, 0, 0);
  PosGantryShield = G4ThreeVector(PipeVertR * cos(DefaultAngle), -PipeVertR * sin(DefaultAngle), GantryShield_z);

  G4RotationMatrix* RotationGantryShield = new G4RotationMatrix();  // no rotation
  G4Transform3D* transformGantryShield = new G4Transform3D(*RotationGantryShield, PosGantryShield);

  GantryShieldLog = new G4LogicalVolume(GantryShieldCyl, matShield, "GantryShieldLog");
  GantryShieldPhys = new G4PVPlacement(*transformGantryShield, "GantryShieldPhys", GantryShieldLog, fMotherVolume, false, 0, fCheckOverlaps);

  G4VisAttributes* GantryShield_vis = new G4VisAttributes(myGrey);
  GantryShield_vis->SetVisibility(true);
  GantryShield_vis->SetForceSolid(true);
  GantryShieldLog->SetVisAttributes(GantryShield_vis);

  // get mass of the pipes
  G4double mass_GantryShield = 0.;
  mass_GantryShield += GantryShieldLog->GetMass(false, false) / CLHEP::kg;
  G4cout << "Mass of the Gantry Shield: " << mass_GantryShield << " kg" << G4endl;
}
