#ifndef DSDetectorProtoProto_H
#define DSDetectorProtoProto_H

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

class DSDetectorProtoProto {

 public:
  DSDetectorProtoProto(G4VPhysicalVolume*);
  ~DSDetectorProtoProto();

 private:
  void DefineSurfaces();

  // G4Polyhedra* fSolidGrid;
  // G4LogicalVolume* fLogicGrid;

  G4VPhysicalVolume* fMotherVolume;
  G4VPhysicalVolume* fPhysic_LAr;
  G4VPhysicalVolume* fPhysic_GAr;
  G4VPhysicalVolume* fPhysic_vessel;
  G4VPhysicalVolume* fPhysic_LArBath;
  G4VPhysicalVolume* fPhysic_SiPmTop;
  G4VPhysicalVolume* fPhysic_SiPMBottom;
  G4VPhysicalVolume* fPhysic_TPB;

  // G4VPhysicalVolume* fPhysicDS20kVac;
  // G4VPhysicalVolume* fPhysicActiveLAr;
  // G4VPhysicalVolume* fPhysicSiPmBottom;
  // G4VPhysicalVolume* fPhysicGrid;
  // G4VPhysicalVolume* fPhysicGasPocket;
  // G4VPhysicalVolume* fPhysicTPBSide;
  // G4VPhysicalVolume* fPhysicTPBBottom;
  // G4VPhysicalVolume* fPhysicLArLayer;
  // G4VPhysicalVolume* fPhysicTopWindow;
  // G4VPhysicalVolume* fPhysicAcrylic;
  // G4VPhysicalVolume* fPhysicBotWindow;
  // G4VPhysicalVolume* fPhysicTeflonBottom;
  // G4VPhysicalVolume* fPhysicDS20k;
};

#endif
