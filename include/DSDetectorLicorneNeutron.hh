#ifndef DSDetectorLicorneNeutron_H
#define DSDetectorLicorneNeutron_H

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

class DSDetectorLicorneNeutron {
 public:
  DSDetectorLicorneNeutron(G4VPhysicalVolume*);
  ~DSDetectorLicorneNeutron();

  // G4VPhysicalVolume*          GetDetectorComponent()  { return  fPhysicOuterLiqArgon; }

 private:
  static G4double AngleFromNuclearRecoilEnergy(G4double G4E_ni, G4double G4E_af);

 private:
  void DefineSurfaces();

  //  double 			GetCm(double );

  G4VPhysicalVolume* fMotherVolume;
};

#endif
