#ifndef DSDetectorLicorne_H
#define DSDetectorLicorne_H

#include "G4LogicalBorderSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4OpticalSurface.hh"
#include "G4Polycone.hh"
#include "G4SubtractionSolid.hh"
#include "G4Tubs.hh"
#include "G4UnionSolid.hh"
#include "G4VPhysicalVolume.hh"

class DSMaterial;

class DSDetectorLicorne {
 public:
  DSDetectorLicorne(G4VPhysicalVolume*);
  ~DSDetectorLicorne();

  // G4VPhysicalVolume*          GetDetectorComponent()  { return  fPhysicOuterLiqArgon; }

 private:
  void DefineSurfaces();

  double GetCm(double);

  G4VPhysicalVolume* fMotherVolume;

  // Experimental room
  G4VPhysicalVolume* fPhysicLicorneFloor;
  G4VPhysicalVolume* fPhysicLeadWall[6];

  // Cryostat and vacuum
  G4VPhysicalVolume* fPhysicLicorne;
  G4VPhysicalVolume* fPhysicLicorneVac;
  G4VPhysicalVolume* fPhysicInnerCryo;

  // TPC
  G4VPhysicalVolume* fPhysicInactiveLar;
  G4VPhysicalVolume* fPhysicExternalSleeve;
  G4VPhysicalVolume* fPhysicCopperRing[20];
  G4VPhysicalVolume* fActiveLArTPC;
  G4VPhysicalVolume* fPhysicGasPocket;
  G4VPhysicalVolume* fPhysicGrid;
  G4VPhysicalVolume* fPhysicCathodeWindow;
  G4VPhysicalVolume* fPhysicAnodeWindow;

  // TPB Layer
  G4VPhysicalVolume* fPhysicTPBGAr;
  G4VPhysicalVolume* fPhysicTPBLAr;

  // PMTs
  G4VPhysicalVolume* fPhysicTopPMTBody1[7];
  G4VPhysicalVolume* fPhysicTopPMTBody1_Window[7];
  G4VPhysicalVolume* fPhysicTopPMTBody1_Vac[7];
  G4VPhysicalVolume* fPhysicBottomPMTWindow;
  G4VPhysicalVolume* fPhysicBottomPMTBody1;
  G4VPhysicalVolume* fPhysicBottomPMTBody1Vac;
  G4VPhysicalVolume* fPhysicBottomPMTBody2;
  G4VPhysicalVolume* fPhysicBottomPMTBody2Vac;
};

#endif
