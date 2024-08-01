#ifndef DSDetectorARISERGamma_H
#define DSDetectorARISERGamma_H

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

class DSDetectorARISERGamma {
 public:
  DSDetectorARISERGamma(G4VPhysicalVolume*);
  ~DSDetectorARISERGamma();
  inline double GetCm(double inches) { return inches * 2.54; }

 private:
  G4VPhysicalVolume* fMotherVolume;

  G4VPhysicalVolume* fPhysicSourceHolder;
  G4VPhysicalVolume* fPhysicSourceActive;

  G4VPhysicalVolume* fPhysicBafShield511;
  G4VPhysicalVolume* fPhysicBafShield1270;
  G4VPhysicalVolume* fPhysicBafDetector511;
  G4VPhysicalVolume* fPhysicBafDetector1270;

  G4VPhysicalVolume* fPhysicBegeShield;
  G4VPhysicalVolume* fPhysicBegeWindow;
  G4VPhysicalVolume* fPhysicBegeVacuum;
  G4VPhysicalVolume* fPhysicBegeCrystal;

 private:
  void DefineSurfaces();
};

#endif
