#ifndef DSDetectorDuneCryostat_H
#define DSDetectorDuneCryostat_H

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4Para.hh"
#include "G4Trap.hh"
#include "G4Tubs.hh"
#include "G4VPhysicalVolume.hh"

class DSMaterial;

class DSDetectorDuneCryostat {

 public:
  DSDetectorDuneCryostat(G4VPhysicalVolume*);
  ~DSDetectorDuneCryostat();

  G4VPhysicalVolume* GetDetectorComponent() { return fPhysicVetoArgon; }

 private:
  // DSMaterial* dsmaterial;
  // G4Material* fWorldMaterial;

  G4VPhysicalVolume* fMotherVolume;

  G4Box* fSolidSteelTank;
  G4LogicalVolume* fLogicSteelTank;
  G4VPhysicalVolume* fPhysicSteelTank;

  G4Box* fSolidOuterPlyWood;
  G4LogicalVolume* fLogicOuterPlyWood;
  G4VPhysicalVolume* fPhysicOuterPlyWood;

  G4Box* fSolidOuterInsulatingFoam;
  G4LogicalVolume* fLogicOuterInsulatingFoam;
  G4VPhysicalVolume* fPhysicOuterInsulatingFoam;

  G4Box* fSolidOuterInsulatingFoamRoof;
  G4LogicalVolume* fLogicOuterInsulatingFoamRoof;
  G4VPhysicalVolume* fPhysicOuterInsulatingFoamRoof;

  G4Box* fSolidInsulatingFoam;
  G4LogicalVolume* fLogicInsulatingFoam;
  G4VPhysicalVolume* fPhysicInsulatingFoam;

  G4Box* fSolidInsulatingFoamRoof;
  G4LogicalVolume* fLogicInsulatingFoamRoof;
  G4VPhysicalVolume* fPhysicInsulatingFoamRoof;

  G4Box* fSolidMembrane;
  G4LogicalVolume* fLogicMembrane;
  G4VPhysicalVolume* fPhysicMembrane;

  G4Box* fSolidInnerMembrane;
  G4LogicalVolume* fLogicInnerMembrane;
  G4VPhysicalVolume* fPhysicInnerMembrane;

  G4Box* fSolidPlyWood;
  G4LogicalVolume* fLogicPlyWood;
  G4VPhysicalVolume* fPhysicPlyWood;

  G4Box* fSolidVetoArgon;
  G4LogicalVolume* fLogicVetoArgon;
  G4VPhysicalVolume* fPhysicVetoArgon;

  G4Box* fSolidArgonVapor;
  G4LogicalVolume* fLogicArgonVapor;
  G4VPhysicalVolume* fPhysicArgonVapor;

  G4Box* fSolidAdditionalPlastic;
  G4LogicalVolume* fLogicAdditionalPlastic;
  G4VPhysicalVolume* fPhysicAdditionalPlastic;
};

#endif
