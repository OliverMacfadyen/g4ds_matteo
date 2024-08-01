#include "DSDetectorPET.hh"
#include <fstream>
#include <iostream>
#include "DSLogger.hh"
#include "DSMaterial.hh"
#include "DSStorage.hh"
#include "G4Colour.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4UIcommand.hh"
#include "G4VisAttributes.hh"

using namespace std;

DSDetectorPET::DSDetectorPET(G4VPhysicalVolume* myMotherVolume) {

  fMotherVolume = myMotherVolume;

  DSLog(routine) << " Constructing PET Geometry" << endlog;
  G4double myLArRadius = 5.0 * cm;
  G4double myLArHeight = 20.0 * cm;

  G4double mySteelThickness = 1.0 * mm;
  G4double myTeflonThickness = 5.0 * mm;
  G4double myTPBThickness = .1 * mm;
  G4double mySiPMThickness = 2.0 / 2. * mm;

  G4double myTPBRadius = myLArRadius + myTPBThickness;
  G4double myTPBHeight = myLArHeight + myTPBThickness;

  G4double myTeflonRadius = myTPBRadius + myTeflonThickness;
  G4double myTeflonHeight = myLArHeight + myTeflonThickness;

  G4double mySteelRadius = myTeflonRadius + mySteelThickness;
  G4double mySteelHeight = myTeflonHeight + mySteelThickness;

  G4double mySiPMRadius = myLArRadius;
  G4double mySiPMHeight = mySiPMThickness;

  G4double myTPBTopRadius = myLArRadius;
  G4double myTPBTopHeight = myTPBThickness;
  G4double myTPBTopShift = myLArHeight - 2 * mySiPMThickness - myTPBTopHeight;

  G4double myTPBBotRadius = myLArRadius;
  G4double myTPBBotHeight = myTPBThickness;
  G4double myTPBBotShift = -myLArHeight + 2 * mySiPMThickness + myTPBTopHeight;

  G4double myPETShift = 50 * cm;

  // Steel Cryostat
  fSolidSteelTank = new G4Tubs("Steel_Solid", 0, mySteelRadius, mySteelHeight, 0, twopi * rad);
  fLogicSteelTank = new G4LogicalVolume(fSolidSteelTank, DSMaterial::Get()->GetSteel(), "Steel_Logic");
  fPhysicSteelTank[0] = new G4PVPlacement(0, G4ThreeVector(0, 0, myPETShift), "Steel1", fLogicSteelTank, fMotherVolume, false, 1, DSStorage::Get()->GetCheckOverlap());

  fPhysicSteelTank[1] = new G4PVPlacement(0, G4ThreeVector(0, 0, -myPETShift), "Steel1", fLogicSteelTank, fMotherVolume, false, 1, DSStorage::Get()->GetCheckOverlap());

  // Teflon volume
  fSolidTeflon = new G4Tubs("Teflon_Solid", 0, myTeflonRadius, myTeflonHeight, 0, twopi * rad);
  fLogicTeflon = new G4LogicalVolume(fSolidTeflon, DSMaterial::Get()->GetTeflon(), "Teflon_Logic");
  fPhysicTeflon = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), "Teflon", fLogicTeflon, fPhysicSteelTank[0], false, 1, DSStorage::Get()->GetCheckOverlap());

  // TPB volume
  fSolidTPB = new G4Tubs("TPB_Solid", myLArRadius, myTPBRadius, myTPBHeight, 0, twopi * rad);
  fLogicTPB = new G4LogicalVolume(fSolidTPB, DSMaterial::Get()->GetTPB(), "TPB_Logic");
  fPhysicTPB = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), "TPB", fLogicTPB, fPhysicTeflon, false, 1, DSStorage::Get()->GetCheckOverlap());

  // LAr volume
  fSolidLAr = new G4Tubs("LAr_Solid", 0, myLArRadius, myLArHeight, 0, twopi * rad);
  fLogicLAr = new G4LogicalVolume(fSolidLAr, DSMaterial::Get()->GetLiquidArgon(), "ActiveLAr_Logic");
  fPhysicLAr = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), "LAr", fLogicLAr, fPhysicTeflon, false, 1, DSStorage::Get()->GetCheckOverlap());

  // SiPM Top
  fSolidSiPMTop = new G4Tubs("SiPMTop_Solid", 0, mySiPMRadius, mySiPMHeight, 0, twopi * rad);
  fLogicSiPMTop = new G4LogicalVolume(fSolidSiPMTop, DSMaterial::Get()->GetMetalSilicon(), "SiPMTop_Logic");
  fPhysicSiPMTop = new G4PVPlacement(0, G4ThreeVector(0, 0, +myLArHeight - mySiPMHeight), "SiPMTop", fLogicSiPMTop, fPhysicLAr, false, 1, DSStorage::Get()->GetCheckOverlap());

  // SiPM Bot
  fSolidSiPMBot = new G4Tubs("SiPMBot_Solid", 0, mySiPMRadius, mySiPMHeight, 0, twopi * rad);
  fLogicSiPMBot = new G4LogicalVolume(fSolidSiPMBot, DSMaterial::Get()->GetMetalSilicon(), "SiPMBot_Logic");
  fPhysicSiPMBot = new G4PVPlacement(0, G4ThreeVector(0, 0, -myLArHeight + mySiPMHeight), "SiPMBot", fLogicSiPMBot, fPhysicLAr, false, 1, DSStorage::Get()->GetCheckOverlap());

  // TPB Top
  fSolidTPBTop = new G4Tubs("TPBTop_Solid", 0, myTPBTopRadius, myTPBTopHeight, 0, twopi * rad);
  fLogicTPBTop = new G4LogicalVolume(fSolidTPBTop, DSMaterial::Get()->GetTPB(), "TPBTop_Logic");
  fPhysicTPBTop = new G4PVPlacement(0, G4ThreeVector(0, 0, myTPBTopShift), "TPBTop", fLogicTPBTop, fPhysicLAr, false, 1, DSStorage::Get()->GetCheckOverlap());

  // TPB Bot
  fSolidTPBBot = new G4Tubs("TPBBot_Solid", 0, myTPBBotRadius, myTPBBotHeight, 0, twopi * rad);
  fLogicTPBBot = new G4LogicalVolume(fSolidTPBBot, DSMaterial::Get()->GetTPB(), "TPBBot_Logic");
  fPhysicTPBBot = new G4PVPlacement(0, G4ThreeVector(0, 0, myTPBBotShift), "TPBBot", fLogicTPBBot, fPhysicLAr, false, 1, DSStorage::Get()->GetCheckOverlap());

  G4Region* fLArRegion = new G4Region("LAr_Logic");
  fLogicLAr->SetRegion(fLArRegion);
  fLArRegion->AddRootLogicalVolume(fLogicLAr);

  // make SiPM as pe storing material
  DSStorage::Get()->SetPMTMaterialIndex(fLogicSiPMBot->GetMaterial()->GetIndex());

  G4VisAttributes* mySiPMAttributes = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0, 0.9));
  G4VisAttributes* myLArAttributes = new G4VisAttributes(G4Colour(.0, 4.0, 2.0, 0.1));
  fLogicSiPMBot->SetVisAttributes(mySiPMAttributes);
  fLogicSiPMTop->SetVisAttributes(mySiPMAttributes);
  fLogicLAr->SetVisAttributes(myLArAttributes);

  DefineSurfaces();
}
DSDetectorPET::~DSDetectorPET() {
  ;
}
void DSDetectorPET::DefineSurfaces() {

  // Steel - Teflon   OK?
  G4OpticalSurface* fOpSteelTeflonSurface = new G4OpticalSurface("OpSteelTeflonSurface");
  new G4LogicalBorderSurface("OpSteelTeflonBorderSurface0", fPhysicTeflon, fPhysicSteelTank[0], fOpSteelTeflonSurface);
  new G4LogicalBorderSurface("OpSteelTeflonBorderSurface1", fPhysicTeflon, fPhysicSteelTank[1], fOpSteelTeflonSurface);
  fOpSteelTeflonSurface->SetType(dielectric_metal);
  fOpSteelTeflonSurface->SetModel(glisur);
  fOpSteelTeflonSurface->SetFinish(ground);
  G4MaterialPropertiesTable* fMTBSteelTeflon = new G4MaterialPropertiesTable();
  fMTBSteelTeflon->AddConstProperty("REFLECTIVITY", 1.0);
  fMTBSteelTeflon->AddConstProperty("EFFICIENCY", 0.0);
  fOpSteelTeflonSurface->SetMaterialPropertiesTable(fMTBSteelTeflon);

  // Teflon - TPB  OK?
  G4OpticalSurface* fOpTeflonTPBSurface = new G4OpticalSurface("OpTeflonTPBSurface");
  new G4LogicalBorderSurface("OpTeflonTPBBorderSurface", fPhysicTPB, fPhysicTeflon, fOpTeflonTPBSurface);
  fOpTeflonTPBSurface->SetType(dielectric_metal);
  fOpTeflonTPBSurface->SetModel(unified);
  fOpTeflonTPBSurface->SetFinish(groundfrontpainted);
  fOpTeflonTPBSurface->SetSigmaAlpha(0.1);
  G4MaterialPropertiesTable* fMTBTeflonTPB = new G4MaterialPropertiesTable();
  G4double TeflonTPBENE[2] = {0.1 * eV, 20.0 * eV};
  G4double TeflonTPBREF[2] = {0.90, 0.90};
  G4double TeflonTPBEFF[2] = {0.00, 0.00};
  fMTBTeflonTPB->AddProperty("REFLECTIVITY", TeflonTPBENE, TeflonTPBREF, 2);
  fMTBTeflonTPB->AddProperty("EFFICIENCY", TeflonTPBENE, TeflonTPBEFF, 2);

  // SiPM - Teflon
  G4OpticalSurface* fOpSiPMTeflonSurface = new G4OpticalSurface("OpSiPMTeflonSurface");
  new G4LogicalBorderSurface("OpSiPMTeflonBorderSurface", fPhysicSiPMBot, fPhysicTeflon, fOpSiPMTeflonSurface);
  fOpSiPMTeflonSurface->SetType(dielectric_dielectric);
  fOpSiPMTeflonSurface->SetModel(unified);
  fOpSiPMTeflonSurface->SetFinish(polished);
  fOpSiPMTeflonSurface->SetSigmaAlpha(0.0);
  G4MaterialPropertiesTable* fMTBSiPMTeflon = new G4MaterialPropertiesTable();
  // fMTBLArTeflon->AddConstProperty("EFFICIENCY", 0.0);
  // fMTBLArTeflon->AddConstProperty("REFLECTIVITY", 0.1);
  // fMTBLArTeflon->AddConstProperty("SPECULARLOBECONSTANT",  1.0);
  // fMTBLArTeflon->AddConstProperty("BACKSCATTERCONSTANT",   1.0);
  // fMTBLArTeflon->AddConstProperty("SPECULARSPIKECONSTANT", 1.0);
  fOpSiPMTeflonSurface->SetMaterialPropertiesTable(fMTBSiPMTeflon);

  // SiPM - TPB
  G4OpticalSurface* fOpSiPMTPBSurface = new G4OpticalSurface("OpSiPMTPBSurface");
  new G4LogicalBorderSurface("OpSiPMTPBTopBorderSurface", fPhysicSiPMTop, fPhysicTPBTop, fOpSiPMTPBSurface);
  new G4LogicalBorderSurface("OpSiPMTPBBotBorderSurface", fPhysicSiPMBot, fPhysicTPBBot, fOpSiPMTPBSurface);
  fOpSiPMTPBSurface->SetType(dielectric_dielectric);
  fOpSiPMTPBSurface->SetModel(unified);
  fOpSiPMTPBSurface->SetFinish(polished);
  fOpSiPMTPBSurface->SetSigmaAlpha(0.0);
  G4MaterialPropertiesTable* fMTBSiPMTPB = new G4MaterialPropertiesTable();
  // fMTBLArTPB->AddConstProperty("EFFICIENCY", 0.0);
  // fMTBLArTPB->AddConstProperty("REFLECTIVITY", 0.1);
  // fMTBLArTPB->AddConstProperty("SPECULARLOBECONSTANT",  1.0);
  // fMTBLArTPB->AddConstProperty("BACKSCATTERCONSTANT",   1.0);
  // fMTBLArTPB->AddConstProperty("SPECULARSPIKECONSTANT", 1.0);
  fOpSiPMTPBSurface->SetMaterialPropertiesTable(fMTBSiPMTPB);

  // TPB - LAr
  G4OpticalSurface* fOpLArTPBSurface = new G4OpticalSurface("OpLArTPBSurface");
  new G4LogicalBorderSurface("OpLArTPBBorderSurface", fPhysicTPB, fPhysicLAr, fOpLArTPBSurface);
  new G4LogicalBorderSurface("OpLArTPBTopBorderSurface", fPhysicTPBTop, fPhysicLAr, fOpLArTPBSurface);
  new G4LogicalBorderSurface("OpLArTPBBotBorderSurface", fPhysicTPBBot, fPhysicLAr, fOpLArTPBSurface);
  fOpLArTPBSurface->SetType(dielectric_dielectric);
  fOpLArTPBSurface->SetModel(unified);
  fOpLArTPBSurface->SetFinish(ground);
  fOpLArTPBSurface->SetSigmaAlpha(0.3);
  G4MaterialPropertiesTable* fMTBLArTPB = new G4MaterialPropertiesTable();
  // fMTBLArTPB->AddConstProperty("EFFICIENCY", 0.0);
  // fMTBLArTPB->AddConstProperty("REFLECTIVITY", 0.1);
  // fMTBLArTPB->AddConstProperty("SPECULARLOBECONSTANT",  1.0);
  // fMTBLArTPB->AddConstProperty("BACKSCATTERCONSTANT",   1.0);
  // fMTBLArTPB->AddConstProperty("SPECULARSPIKECONSTANT", 1.0);
  fOpLArTPBSurface->SetMaterialPropertiesTable(fMTBLArTPB);

  // Steel - LAr
  // G4OpticalSurface *fOpSteelLArSurface     = new
  // G4OpticalSurface("OpSteelLArSurface"); new
  // G4LogicalBorderSurface("OpSteelLArBorderSurface", fPhysicLAr,
  // fPhysicSteelTank, fOpSteelLArSurface ); fOpSteelLArSurface->SetType(
  // dielectric_dielectric ); fOpSteelLArSurface->SetModel( glisur );
  // fOpSteelLArSurface->SetFinish( polished );
  // G4MaterialPropertiesTable *fMTBSteelLAr = new G4MaterialPropertiesTable();
  // fMTBSteelLAr->AddConstProperty("REFLECTIVITY", 0.0);
  // fMTBSteelLAr->AddConstProperty("EFFICIENCY", 0.0);
  // fOpSteelLArSurface->SetMaterialPropertiesTable( fMTBSteelLAr );
}

/*
 * $Log: DSDetectorPET.cc,v $
 * Revision 1.1  2016/02/05 16:43:44  pagnes
 * basic PET geometry added
 *
 *
 */
