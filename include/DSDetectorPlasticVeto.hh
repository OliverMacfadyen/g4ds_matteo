#ifndef DSDetectorPlasticVeto_H
#define DSDetectorPlasticVeto_H

#include <vector>
#include "G4Box.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4OpticalSurface.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Tubs.hh"
#include "G4UnionSolid.hh"
#include "G4VPhysicalVolume.hh"

class DSMaterial;

class DSDetectorPlasticVeto {

 public:
  DSDetectorPlasticVeto(G4VPhysicalVolume*);
  ~DSDetectorPlasticVeto();

  G4VPhysicalVolume* GetDetectorComponent() { return fPhysicInnerLiquidArgon; }

 private:
  void DefineSurfaces();

  // Functions useful for manual SiPMs placement from file
  G4int get_sector(G4int);
  G4int get_buffer(G4int);
  G4int get_place(G4int);

  G4int CalculateCoverage(G4bool isTop, G4double percent, G4double R, G4double Z, G4double SiPM_Area);
  void PlaceSiPMsInCaps(G4double, G4double, G4double, G4int, G4LogicalVolume*, G4VPhysicalVolume*, int);
  void PlaceSiPMsInSides(G4double, G4double, G4double, G4int, G4LogicalVolume*, G4VPhysicalVolume*, int);
  int placedSiPM;

  DSMaterial* dsmaterial;
  G4LogicalVolume* smallSiPMLogic;
  G4LogicalVolume* CellSiPMLogic;
  G4LogicalVolume* ActiveSiPMLogic;
  G4LogicalVolume* smallTPBLogic;
  G4LogicalVolume* smallSilicaLogic;
  G4VPhysicalVolume* CoatedActiveSiPM;
  G4VPhysicalVolume* FullActiveSiPM;
  G4VPhysicalVolume* SilicaActiveSiPM;
  G4VPhysicalVolume* fMotherVolume;
  G4VPhysicalVolume* fPhysicInnerLiquidArgon;

  G4VPhysicalVolume* physSiPMSilicaOnSide;
  G4VPhysicalVolume* physSiPMTPBOnSide;
  G4VPhysicalVolume* physSiPMActiveOnSide;
  G4VPhysicalVolume* physSiPMSilicaOnCap;
  G4VPhysicalVolume* physSiPMTPBOnCap;
  G4VPhysicalVolume* physSiPMActiveOnCap;
  G4VPhysicalVolume* physSiPMInternalSide_1;
  G4VPhysicalVolume* physSiPMInternalSide_2;
  G4VPhysicalVolume* physSiPMExternalSide_1;
  G4VPhysicalVolume* physSiPMExternalSide_2;
  G4VPhysicalVolume* physSiPMInternalUpperCap_1;
  G4VPhysicalVolume* physSiPMInternalUpperCap_2;
  G4VPhysicalVolume* physSiPMInternalLowerCap_1;
  G4VPhysicalVolume* physSiPMInternalLowerCap_2;
  G4VPhysicalVolume* physSiPMExternalUpperCap_1;
  G4VPhysicalVolume* physSiPMExternalUpperCap_2;
  G4VPhysicalVolume* physSiPMExternalLowerCap_1;
  G4VPhysicalVolume* physSiPMExternalLowerCap_2;

  G4VPhysicalVolume* physSiPMManual;

  G4VPhysicalVolume* physCopper;
  G4VPhysicalVolume* physExternalTPB;
  G4VPhysicalVolume* physOuterLArBuffer;
  G4VPhysicalVolume* physTPBOuterPanelFace[8];
  G4VPhysicalVolume* physAcrylicOuterPanelFace[8];
  G4VPhysicalVolume* physTPBOuterPanelEdge[8];
  G4VPhysicalVolume* physAcrylicOuterPanelEdge[8];
  G4VPhysicalVolume* physUpperOuterTPBCylinder;
  G4VPhysicalVolume* physUpperOuterAcrylicCylinder;
  G4VPhysicalVolume* physLowerOuterTPBCylinder;
  G4VPhysicalVolume* physLowerOuterAcrylicCylinder;
  G4VPhysicalVolume* physOuterPlasticTPB;
  G4VPhysicalVolume* physPlasticLArBuffer;
  G4VPhysicalVolume* physPlastic;
  G4VPhysicalVolume* physInnerPlasticTPB;
  G4VPhysicalVolume* physInnerLArBuffer;
  G4VPhysicalVolume* physTPBInnerPanelFace[8];
  G4VPhysicalVolume* physAcrylicInnerPanelFace[8];
  G4VPhysicalVolume* physTPBInnerPanelEdge[8];
  G4VPhysicalVolume* physAcrylicInnerPanelEdge[8];
  G4VPhysicalVolume* physUpperInnerTPBCylinder;
  G4VPhysicalVolume* physUpperInnerAcrylicCylinder;
  G4VPhysicalVolume* physLowerInnerTPBCylinder;
  G4VPhysicalVolume* physLowerInnerAcrylicCylinder;

  std::vector<G4VPhysicalVolume*> physSiPMInside;
  std::vector<G4VPhysicalVolume*> physSiPMOutside;
};
#endif
