#ifndef DSDetectorPET_H
#define DSDetectorPET_H

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4Para.hh"
#include "G4Trap.hh"
#include "G4Tubs.hh"
#include "G4VPhysicalVolume.hh"

class DSMaterial;

class DSDetectorPET {

 public:
  DSDetectorPET(G4VPhysicalVolume*);
  ~DSDetectorPET();

  G4VPhysicalVolume* GetDetectorComponent() { return fPhysicSteelTank[0]; }

 private:
  void DefineSurfaces();

  G4VPhysicalVolume* fMotherVolume;

  G4Tubs* fSolidSteelTank;
  G4LogicalVolume* fLogicSteelTank;
  G4VPhysicalVolume* fPhysicSteelTank[2];

  G4Tubs* fSolidTeflon;
  G4LogicalVolume* fLogicTeflon;
  G4VPhysicalVolume* fPhysicTeflon;

  G4Tubs* fSolidTPB;
  G4LogicalVolume* fLogicTPB;
  G4VPhysicalVolume* fPhysicTPB;

  G4Tubs* fSolidLAr;
  G4LogicalVolume* fLogicLAr;
  G4VPhysicalVolume* fPhysicLAr;

  G4Tubs* fSolidSiPMTop;
  G4LogicalVolume* fLogicSiPMTop;
  G4VPhysicalVolume* fPhysicSiPMTop;

  G4Tubs* fSolidSiPMBot;
  G4LogicalVolume* fLogicSiPMBot;
  G4VPhysicalVolume* fPhysicSiPMBot;

  G4Tubs* fSolidTPBTop;
  G4LogicalVolume* fLogicTPBTop;
  G4VPhysicalVolume* fPhysicTPBTop;

  G4Tubs* fSolidTPBBot;
  G4LogicalVolume* fLogicTPBBot;
  G4VPhysicalVolume* fPhysicTPBBot;
};

#endif
/*
 * $Log: DSDetectorPET.hh,v $
 * Revision 1.1  2016/02/05 16:43:47  pagnes
 * basic PET geometry added
 *
 *
 */
