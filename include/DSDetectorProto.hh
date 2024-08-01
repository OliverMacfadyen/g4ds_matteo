#ifndef DSDetectorProto_H
#define DSDetectorProto_H

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

typedef std::pair<double, double> PairFF;
typedef std::vector<std::pair<double, double>> PointCol;
typedef std::vector<std::pair<double, double>>::iterator PointColItr;
typedef std::vector<std::pair<double, double>>* PointColPtr;

class DSMaterial;

class DSDetectorProto {

 public:
  DSDetectorProto(G4VPhysicalVolume*);
  ~DSDetectorProto();

  //  G4VPhysicalVolume*          GetDetectorComponent()  { return
  //  fPhysicOuterLiqArgon; }
  // G4VPhysicalVolume*          GetDetectorComponent()  { return
  // fPhysicCryoOuter; }

 private:
  void DefineSurfaces();
  PointColPtr createGeometry(double, double, double, double, double, double, double, double);

  double GetOctagonInnerRadius(double);
  double GetOctagonOuterRadius(double);

  G4Polyhedra* fSolidGrid;
  G4LogicalVolume* fLogicGrid;

  G4VPhysicalVolume* fMotherVolume;
  G4VPhysicalVolume* fPhysicDS20kVac;
  G4VPhysicalVolume* fPhysicActiveLAr;
  G4VPhysicalVolume* fPhysicSiPmTop;
  G4VPhysicalVolume* fPhysicSiPmBottom;
  G4VPhysicalVolume* fPhysicGrid;
  G4VPhysicalVolume* fPhysicGasPocket;
  G4VPhysicalVolume* fPhysicTPBSide;
  G4VPhysicalVolume* fPhysicTPBTop;
  // G4VPhysicalVolume* fPhysicTPBBottom;
  G4VPhysicalVolume* fPhysicLArLayer;
  G4VPhysicalVolume* fPhysicTopWindow;
  // G4VPhysicalVolume* fPhysicAcrylic;
  G4VPhysicalVolume* fPhysicBotWindow;
  G4VPhysicalVolume* fPhysicTeflonBottom;

  G4VPhysicalVolume* fPhysicDS20k;

  // optical boundary: TPC cryostat - scintillator
  G4LogicalBorderSurface* fDS20kOuterSurface;

  G4VPhysicalVolume* fPhysicTopCryoFiller;

  G4VPhysicalVolume* fPhysicInactiveLar;
  G4VPhysicalVolume* fPhysicInactiveGar;
};

#endif
