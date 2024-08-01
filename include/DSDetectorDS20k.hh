#ifndef DSDetectorDS20k_H
#define DSDetectorDS20k_H

#include "G4Ellipsoid.hh"
#include "G4IntersectionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4Polycone.hh"
#include "G4Polyhedra.hh"
#include "G4SubtractionSolid.hh"
#include "G4Tubs.hh"
#include "G4UnionSolid.hh"
#include "G4VPhysicalVolume.hh"

#include "G4LogicalBorderSurface.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4OpticalSurface.hh"

class DSMaterial;

class DSDetectorDS20k {

 public:
  DSDetectorDS20k(G4VPhysicalVolume*);
  ~DSDetectorDS20k();
  int get_sector(double, double);
  double Angle2D(double, double, double, double);
  bool InsideOctagon(float, double, double);

 private:
  void DefineSurfaces();

  double GetOctagonInnerRadius(double);
  double GetOctagonOuterRadius(double);

  G4Polyhedra* fSolidGrid;
  G4LogicalVolume* fLogicGrid;

  G4VPhysicalVolume* fMotherVolume;
  G4VPhysicalVolume* myTPBReflectorSkin;
  G4VPhysicalVolume* myNSLArReflectorSkin;
  G4VPhysicalVolume* myPENReflectorSkin;
  G4VPhysicalVolume* myGdReflectorSkin;
  G4VPhysicalVolume* fPhysicActiveLAr;
  G4VPhysicalVolume* fPhysicSiPmTop;
  G4VPhysicalVolume* fPhysicSiPmBottom;
  G4VPhysicalVolume* fPhysicGrid;
  G4VPhysicalVolume* fPhysicGasPocket;
  G4VPhysicalVolume* fPhysicTPBSide;
  G4VPhysicalVolume* fPhysicTPBTop;
  G4VPhysicalVolume* fPhysicTPBBottom;
  G4VPhysicalVolume* fPhysicAcrylicTPCVessel;
  G4VPhysicalVolume* fPhysicLayerAcrylicSand;
  G4VPhysicalVolume* fPhysicLayerNSLAr;
  G4VPhysicalVolume* fPhysicPassiveTop;
  G4VPhysicalVolume* fPhysicPassiveBottom;
  G4VPhysicalVolume* fPhysicLayerESR;
  G4VPhysicalVolume* fPhysicLayerPassiveLAr2;
  G4VPhysicalVolume* fPhysicLArAboveGrid;
  
  G4LogicalVolume * myLogicShapevQuadrant;
  G4LogicalVolume * myLogicShapevPDU;
  G4LogicalVolume* myLogicShapevSiPMSilica;

  G4VPhysicalVolume* physSiPMSilicaOnSide[1280];
  // G4VPhysicalVolume* physSiPMActiveOnSide[80];
  G4VPhysicalVolume* physSiPMSilicaOnCap[640];
  // G4VPhysicalVolume* physSiPMActiveOnCap[40];
  G4VPhysicalVolume* physSiPMInternalUpperCap[20];
  G4VPhysicalVolume* physSiPMInternalLowerCap[20];
  G4VPhysicalVolume* physSiPMInternalSide[80];
};

#endif
