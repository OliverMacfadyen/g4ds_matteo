#ifndef DSDetectorGantryCalibration_H
#define DSDetectorGantryCalibration_H

#include "DSDetectorConstruction.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4Material.hh"
#include "G4OpticalSurface.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalConstants.hh"
#include "G4RotationMatrix.hh"
#include "G4Sphere.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Torus.hh"
#include "G4Transform3D.hh"
#include "G4Tubs.hh"
#include "G4UnionSolid.hh"
#include "G4VSolid.hh"
#include "G4VisAttributes.hh"
#include "globals.hh"

class DSDetectorGantryCalibration {

 public:
  DSDetectorGantryCalibration(G4VPhysicalVolume*);
  ~DSDetectorGantryCalibration();

  void ConstructGantry();
  void DefineGantrySurfaces();
  void ConstructGantryShield();

 private:
  G4VPhysicalVolume* fMotherVolume;
  G4bool fCheckOverlaps;
  G4bool GantryConfiguration;
  G4bool GantrySurface;
  G4int GantrySurfaceOption;
  G4bool GantryShield;
  G4int GantryShieldMaterial;

  G4Material* matSLS;
  G4Material* matNitrogen;
  G4Material* matShield;
  G4double TwoPi;
  G4double DefaultAngle;
  G4double half_tpc_edge;
  G4double TPCWallThickness;
  G4double TPBThickness;
  G4double ElectronicsThickness;
  G4double VetoBufferThickness;
  G4double TPCheight;
  G4double Pipe1HeightTotal;
  G4double Pipe2HeightTotal;
  G4double PipeVertR;
  G4double Pipe1VertHeight;
  G4double Pipe2VertHeight;
  G4double PipeHorLength;

  G4double Rinner;
  G4double Router;
  G4double Rinner_ext;
  G4double Router_ext;
  G4double half_tpc_edge_ext;
  G4double PipeID;
  G4double PipeWall;
  G4double PipeOD;
  G4double PipeBendingR;
  G4double PipeIR;
  G4double PipeOR;
  G4double DistanceTPCside;
  G4double DistanceTPCbottom;
  G4double DistanceBtwPipes;

  G4double GantryShieldThickness;
  G4double GantryShieldHeight;
  G4double GantryShieldOffset;

  G4LogicalVolume* Pipe1log;
  G4LogicalVolume* Pipe2log;
  G4PVPlacement* Pipe1phys;
  G4PVPlacement* Pipe2phys;

  G4LogicalVolume* PipeNitro1log;
  G4LogicalVolume* PipeNitro2log;
  G4PVPlacement* PipeNitro1phys;
  G4PVPlacement* PipeNitro2phys;

  G4LogicalVolume* GantryShieldLog;
  G4PVPlacement* GantryShieldPhys;

  // G4RotationMatrix* DefaultRotation;
};

#endif
