#ifndef DSDetectorReD_H
#define DSDetectorReD_H

#include "G4LogicalVolume.hh"
#include "G4Polycone.hh"
#include "G4SubtractionSolid.hh"
#include "G4Tubs.hh"
#include "G4UnionSolid.hh"
#include "G4VPhysicalVolume.hh"

#include "G4LogicalBorderSurface.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4OpticalSurface.hh"

#include "DSLogger.hh"

class DSMaterial;

class DSDetectorReD {

 public:
  DSDetectorReD(G4VPhysicalVolume*);
  ~DSDetectorReD();

  // G4VPhysicalVolume*          GetDetectorComponent()  { return
  // fPhysicOuterLiqArgon; }
 G4double extended_asin(G4double ang);
 private:
  void DefineSurfaces();
  void SetToVacuum(G4Material*&);
  void Debug(G4VPhysicalVolume* vol) { DSLog(routine) << vol->GetName() << " pos " << vol->GetObjectTranslation() << " in " << vol->GetMotherLogical()->GetName() << endlog; }

  G4double GetCm(G4double);
  G4UnionSolid* G4QuadRingWithFillet(G4double x, G4double h2, G4double d2, G4double R);
  G4SubtractionSolid* G4QuadPlateWithFillet(G4double L2, G4double h2, G4double R, G4double l2 = 0, G4double r = 0);

  G4VPhysicalVolume* fMotherVolume;
  // G4VPhysicalVolume* fPhysicReDCryo;
  // G4VPhysicalVolume* fPhysicPMTTopCap;
  // G4VPhysicalVolume* fPhysicPMTTop;
  // G4VPhysicalVolume* fPhysicPMTBottomCap;
  // G4VPhysicalVolume* fPhysicPMTBottom;
  G4VPhysicalVolume* fPhysicGasPocket;
  G4VPhysicalVolume* fPhysicCathodeWindow;
  // G4VPhysicalVolume* fPhysicDivingBell;
  G4VPhysicalVolume* fPhysicAnodeWindow;
  G4VPhysicalVolume* fPhysicGrid;
  // G4VPhysicalVolume* fPhysicGridRing;
  // G4VPhysicalVolume* fPhysicFieldCage;
  // G4VPhysicalVolume* fPhysicFieldCageUpper;
  // G4VPhysicalVolume* fPhysicExternalSleeve;
  // G4VPhysicalVolume* fPhysicInnerSleeve;
  // G4VPhysicalVolume* fPhysicInnerSleeveGas;
  // G4VPhysicalVolume* fPhysicCopperRing;
  G4VPhysicalVolume* fPhysicTeflonFoil;

  // TPB Layer
  G4VPhysicalVolume* fPhysicTPBGAr;
  G4VPhysicalVolume* fPhysicTPBLAr;

  // Active LAr
  G4VPhysicalVolume* fActiveLArTPC;
  G4VPhysicalVolume* fPhysicInactiveLAr;
  G4VPhysicalVolume* fPhysicTPC;
};

//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
// Catania Beamline
//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------

class CataniaBeamline {
 public:
  CataniaBeamline(G4VPhysicalVolume*, const G4ThreeVector);
  ~CataniaBeamline();

 private:
  void Debug(G4VPhysicalVolume* vol) { DSLog(routine) << vol->GetName() << " pos " << vol->GetObjectTranslation() << " in " << vol->GetMotherLogical()->GetName() << endlog; }
};

#endif
