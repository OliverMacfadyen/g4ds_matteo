#include "DSMaterial.hh"
#include <fstream>
#include "DSEventHandler.hh"
#include "DSLogger.hh"
#include "DSParameters.hh"
#include "DSStorage.hh"
#include "G4NistManager.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

using namespace std;

DSMaterial* DSMaterial::me = 0;

DSMaterial::DSMaterial() {
  fTPCMaterialIndex = -1;
  fIVMaterialIndex = -1;
  fOVMaterialIndex = -1;
  DefineMaterials();
  DefineBoundaries();
  DefineProperties();
  SetMaterialIndexes();
}

DSMaterial::~DSMaterial() {}

DSMaterial* DSMaterial::Get() {
  if (!me) me = new DSMaterial();

  return me;
}

void DSMaterial::DefineMaterials() {

  G4String name, symbol;
  G4int z, ncomponents, natoms;
  G4double a, density, fractionmass, temperature, pressure;

  G4Element* H = new G4Element(name = "Hydrogen", symbol = "H", z = 1, a = 1.00794 * g / mole);
  G4Element* Li = new G4Element(name = "Lithium", symbol = "Li", z = 3, a = 6.941 * g / mole);
  G4Element* C = new G4Element(name = "Carbon", symbol = "C", z = 6, a = 12.0107 * g / mole);
  G4Element* Be = new G4Element(name = "Beryllium", symbol = "Be", z = 4, a = 9.01 * g / mole);
  // G4Element* B =  new G4Element(name="Boron",     symbol="B",  z=5,
  // a=10.811*g/mole  );
  G4Element* N = new G4Element(name = "Nitrogen", symbol = "N", z = 7, a = 14.00067 * g / mole);
  G4Element* O = new G4Element(name = "Oxygen", symbol = "O", z = 8, a = 15.9994 * g / mole);
  G4Element* Fl = new G4Element(name = "Fluorine", symbol = "F", z = 9, a = 18.9984 * g / mole);
  G4Element* Na = new G4Element(name = "Sodium", symbol = "Na", z = 11, a = 22.9898 * g / mole);
  G4Element* Mg = new G4Element(name = "Magnesium", symbol = "Mg", z = 12, a = 24.305 * g / mole);
  G4Element* Al = new G4Element(name = "Aluminum", symbol = "Al", z = 13, a = 26.9815 * g / mole);
  G4Element* Si = new G4Element(name = "Silicium", symbol = "Si", z = 14, a = 28.0855 * g / mole);
  G4Element* S = new G4Element(name = "Sulfur", symbol = "S", z = 16, a = 32.065 * g / mole);
  G4Element* Ar = new G4Element(name = "Argon", symbol = "Ar", z = 18, a = 39.948 * g / mole);
  G4Element* K = new G4Element(name = "Potassium", symbol = "K", z = 19, a = 39.0989 * g / mole);
  G4Element* Ca = new G4Element(name = "Calcium", symbol = "Ca", z = 20, a = 40.078 * g / mole);
  // G4Element* Ti = new G4Element(name="Titanium",  symbol="Ti", z=22,
  // a=47.867*g/mole  );
  G4Element* Cr = new G4Element(name = "Chromium", symbol = "Cr", z = 24, a = 51.9961 * g / mole);
  G4Element* Mn = new G4Element(name = "Manganese", symbol = "Mn", z = 25, a = 54.93805 * g / mole);
  G4Element* Fe = new G4Element(name = "Iron", symbol = "Fe", z = 26, a = 55.845 * g / mole);
  G4Element* Co = new G4Element(name = "Cobalt", symbol = "Co", z = 27, a = 58.9331 * g / mole);
  G4Element* Ni = new G4Element(name = "Nickel", symbol = "Ni", z = 28, a = 58.70 * g / mole);
  // G4Element* Cu = new G4Element(name="Copper",    symbol="Cu", z=29,
  // a=63.546*g/mole  );
  G4Element* Ge = new G4Element(name="Germanium", symbol="Ge", z=32, a=72.630*g/mole);
  G4Element* P = new G4Element(name = "Phosphorus", symbol = "P", z = 15, a = 30.97398 * g / mole);
  G4Element* Rb = new G4Element(name = "Rubidium", symbol = "Rb", z = 37, a = 85.4678 * g / mole);
  G4Element* In = new G4Element(name = "Indium", symbol = "In", z = 49, a = 114.818 * g / mole);
  G4Element* Sn = new G4Element(name = "Selenium", symbol = "Sn", z = 50,
                                a = 118.710 * g / mole);  // Tin????
  G4Element* Sb = new G4Element(name = "Antimony", symbol = "Sb", z = 51, a = 121.760 * g / mole);
  G4Element* Pb = new G4Element(name = "Lead", symbol = "Pb", z = 82, a = 207.2 * g / mole);
  G4Element* I = new G4Element(name = "Iodine", symbol = "I", z = 53, a = 126.9045 * g / mole);
  G4Element* Xe = new G4Element(name = "Xenon", symbol = "Xe", z = 54, a = 131.293 * g / mole);
  G4Element* Cs = new G4Element(name = "Cesium", symbol = "Cs", z = 55, a = 132.91 * g / mole);
  G4Element* Ba = new G4Element(name = "Barium", symbol = "Ba", z = 56, a = 137.327 * g / mole);
  // G4Element* Cl = new G4Element(name="Clorine",   symbol="Cl", z=17,
  // a=35.453*g/mole  ); G4Element* Ta  = new G4Element(name="Tallium",
  // symbol="Ta", z=73, a=180.94*g/mole  ); G4Element* W   = new
  // G4Element(name="Tungsten",  	symbol="W",  z=74, a=182.30*g/mole  );
  G4Element* Ta = G4NistManager::Instance()->FindOrBuildElement("Ta");
  G4Element* W = G4NistManager::Instance()->FindOrBuildElement("W");

  G4Isotope* B10 = new G4Isotope(name = "B10", 5, 10, a = 10.013 * g / mole);
  G4Isotope* B11 = new G4Isotope(name = "B11", 5, 11, a = 11.009 * g / mole);
  G4Element* NewB = new G4Element(name = "NewB", symbol = "NewB", 2);
  NewB->AddIsotope(B10, 19.9 * perCent);
  NewB->AddIsotope(B11, 80.1 * perCent);

  G4Isotope* Li6 = new G4Isotope(name = "Li6", 3, 6, a = 6 * g / mole);
  G4Isotope* Li7 = new G4Isotope(name = "Li7", 3, 7, a = 7 * g / mole);
  G4Element* EnrLi = new G4Element(name = "EnrLi", symbol = "EnrLi", 2);
  EnrLi->AddIsotope(Li6, 95 * perCent);
  EnrLi->AddIsotope(Li7, 5 * perCent);

  G4Element* natLi = new G4Element(name = "natLi", symbol = "natLi", 2);
  natLi->AddIsotope(Li6, 7.5 * perCent);
  natLi->AddIsotope(Li7, 92.5 * perCent);

  G4Isotope* Cl35 = new G4Isotope(name = "Cl35", 17, 35, a = 34.453 * g / mole);
  G4Isotope* Cl37 = new G4Isotope(name = "Cl37", 17, 37, a = 37.453 * g / mole);
  G4Element* natCl = new G4Element(name = "natCl", symbol = "natCl", 2);
  natCl->AddIsotope(Cl35, 75.8 * perCent);
  natCl->AddIsotope(Cl37, 24.2 * perCent);

  G4Isotope* Cu63 = new G4Isotope(name = "Cu63", 29, 63, a = 62.930 * g / mole);
  G4Isotope* Cu65 = new G4Isotope(name = "Cu65", 29, 65, a = 64.928 * g / mole);
  G4Element* Cu = new G4Element(name = "Cu", symbol = "Cu", 2);
  Cu->AddIsotope(Cu63, 69.15 * perCent);
  Cu->AddIsotope(Cu65, 30.85 * perCent);

  G4Isotope* Ti46 = new G4Isotope(name = "Ti46", 22, 46, a = 45.952 * g / mole);
  G4Isotope* Ti47 = new G4Isotope(name = "Ti47", 22, 47, a = 46.951 * g / mole);
  G4Isotope* Ti48 = new G4Isotope(name = "Ti48", 22, 48, a = 47.948 * g / mole);
  G4Isotope* Ti49 = new G4Isotope(name = "Ti49", 22, 49, a = 49.947 * g / mole);
  G4Isotope* Ti50 = new G4Isotope(name = "Ti50", 22, 50, a = 49.994 * g / mole);
  G4Element* Ti = new G4Element(name = "Ti", symbol = "Ti", 5);
  Ti->AddIsotope(Ti46, 8.25 * perCent);
  Ti->AddIsotope(Ti47, 7.44 * perCent);
  Ti->AddIsotope(Ti48, 73.72 * perCent);
  Ti->AddIsotope(Ti49, 5.41 * perCent);
  Ti->AddIsotope(Ti50, 5.18 * perCent);

  G4Isotope* Gd152 = new G4Isotope("Gd152", 64, 152, 152.0 * g / mole);
  G4Isotope* Gd154 = new G4Isotope("Gd154", 64, 154, 154.0 * g / mole);
  G4Isotope* Gd155 = new G4Isotope("Gd155", 64, 155, 155.0 * g / mole);
  G4Isotope* Gd156 = new G4Isotope("Gd156", 64, 156, 156.0 * g / mole);
  G4Isotope* Gd157 = new G4Isotope("Gd157", 64, 157, 157.0 * g / mole);
  G4Isotope* Gd158 = new G4Isotope("Gd158", 64, 158, 158.0 * g / mole);
  G4Isotope* Gd160 = new G4Isotope("Gd160", 64, 160, 160.0 * g / mole);

  G4Element* Gd = new G4Element("Gadolinium", "Gd", 7);
  Gd->AddIsotope(Gd152, 0.2 * perCent);
  Gd->AddIsotope(Gd154, 2.2 * perCent);
  Gd->AddIsotope(Gd155, 14.9 * perCent);  // beware: it is abundance,
  Gd->AddIsotope(Gd156, 20.6 * perCent);  //        not fractionMass
  Gd->AddIsotope(Gd157, 15.7 * perCent);
  Gd->AddIsotope(Gd158, 24.7 * perCent);
  Gd->AddIsotope(Gd160, 21.7 * perCent);

  //---------------------//
  // Material Definition //
  //---------------------//

  // PC
  const G4double densityPC = 0.877 * g / cm3;  // from Borexino
  fPC = new G4Material(name = "PC", densityPC, ncomponents = 2);
  fPC->AddElement(C, natoms = 9);
  fPC->AddElement(H, natoms = 12);

  // TMB
  const G4double densityTMB = 0.932 * g / cm3;  // from Wikipedia
  fTMB = new G4Material(name = "TMB", densityTMB, ncomponents = 4);
  fTMB->AddElement(H, natoms = 9);
  fTMB->AddElement(O, natoms = 3);
  fTMB->AddElement(C, natoms = 3);
  fTMB->AddElement(NewB, natoms = 1);

  const G4double fraction = DSStorage::Get()->GetTMBfraction();
  // fraction = 0.2 for 5us capture time in DS20k
  if ((fraction < 0.) || (fraction > 1.)) DSLog(fatal) << " The scintillator volume fraction must be between 0 and 1" << endlog;
  density = (1. - fraction) * densityPC + fraction * densityTMB;
  fBoronScintillator = new G4Material(name = "BoronScintillator", density, ncomponents = 2);
  fBoronScintillator->AddMaterial(fTMB, fractionmass = fraction);
  fBoronScintillator->AddMaterial(fPC, fractionmass = (1. - fraction));

  G4cout << fBoronScintillator << endl;

  // Air
  density = 1.290 * mg / cm3;
  fAir = new G4Material(name = "Air", density, ncomponents = 2);
  fAir->AddElement(N, fractionmass = 0.7);
  fAir->AddElement(O, fractionmass = 0.3);

  // Vacuum
  density = 1.29e-20 * g / cm3;
  fVacuum = new G4Material(name = "Vacuum", density, ncomponents = 2);
  fVacuum->AddElement(N, fractionmass = 0.7);
  fVacuum->AddElement(O, fractionmass = 0.3);

  // PPO
  density = 1.000 * g / cm3;  // should be corrected
  fPPO = new G4Material(name = "PPO", density, ncomponents = 4);
  fPPO->AddElement(C, natoms = 15);
  fPPO->AddElement(H, natoms = 11);
  fPPO->AddElement(N, natoms = 1);
  fPPO->AddElement(O, natoms = 1);

  // DMP
  density = 1.000 * g / cm3;  // should be corrected
  fDMP = new G4Material(name = "DMP", density, ncomponents = 4);
  fDMP->AddElement(C, natoms = 10);
  fDMP->AddElement(H, natoms = 10);
  fDMP->AddElement(N, natoms = 1);
  fDMP->AddElement(O, natoms = 1);

  // Acrylic
  density = 1.19 * g / cm3;
  fAcrylic = new G4Material(name = "Acrylic", density, ncomponents = 3, kStateSolid);
  fAcrylic->AddElement(C, 5);
  fAcrylic->AddElement(O, 2);
  fAcrylic->AddElement(H, 8);

  // Liquid Argon
  density = 1.40 * g / cm3;
  temperature = 87 * kelvin;
  fLiquidArgon = new G4Material(name = "LiquidArgon", density, ncomponents = 1, kStateLiquid, temperature);
  fLiquidArgon->AddElement(Ar, 1.0);

  // Pseudo Argon
  density = 1.40 * g / cm3;
  temperature = 87 * kelvin;
  fPseudoArgon = new G4Material(name = "PseudoArgon", density, ncomponents = 1, kStateLiquid, temperature);
  fPseudoArgon->AddElement(Ar, 1.0);

  // Non Scintillating Liquid Argon
  density = 1.40 * g / cm3;
  temperature = 87 * kelvin;
  fNSLiquidArgon = new G4Material(name = "NSLiquidArgon", density, ncomponents = 1, kStateLiquid, temperature);
  fNSLiquidArgon->AddElement(Ar, 1.0);

  // Gaseous Argon
  density = 5.4 * mg / cm3;
  temperature = 87.00 * kelvin;
  pressure = 0.93 * atmosphere;
  fGaseousArgon = new G4Material(name = "GaseousArgon", density, ncomponents = 1, kStateGas, temperature, pressure);
  fGaseousArgon->AddElement(Ar, 1.0);

  // Liquid Xenon
  density = 3.057 * g / cm3;
  temperature = 165 * kelvin;
  fLiquidXenon = new G4Material(name = "LiquidXenon", density, ncomponents = 1, kStateLiquid, temperature);
  fLiquidXenon->AddElement(Xe, 1.0);

  // Gaseous Xenon
  density = 5.89 * mg / cm3;
  temperature = 165 * kelvin;
  pressure = 0.93 * atmosphere;
  fGaseousXenon = new G4Material(name = "GaseousXenon", density, ncomponents = 1, kStateGas, temperature, pressure);
  fGaseousXenon->AddElement(Xe, 1.0);

  // Stainless steel (Medical Physics, Vol 25, No 10, Oct 1998)
  density = 8.290 * g / cm3;
  fSteel = new G4Material(name = "Steel", density, ncomponents = 5);
  fSteel->AddElement(Mn, 0.02);
  fSteel->AddElement(Si, 0.01);
  fSteel->AddElement(Cr, 0.19);
  fSteel->AddElement(Ni, 0.10);
  fSteel->AddElement(Fe, 0.68);

  // Stainless Steel
  density = 7.7 * g / cm3;
  fStainlessSteel = new G4Material(name = "StainlessSteel", density, ncomponents = 3);
  fStainlessSteel->AddElement(Cr, fractionmass = 0.20);
  fStainlessSteel->AddElement(Fe, fractionmass = 0.68);
  fStainlessSteel->AddElement(Ni, fractionmass = 0.12);

  // Grid Steel
  density = 7.7 * g / cm3;
  fGridSteel = new G4Material(name = "GridSteel", density, ncomponents = 3);
  fGridSteel->AddElement(Cr, fractionmass = 0.20);
  fGridSteel->AddElement(Fe, fractionmass = 0.68);
  fGridSteel->AddElement(Ni, fractionmass = 0.12);

  // Water
  density = 1.00 * g / cm3;
  fWater = new G4Material(name = "Water", density, ncomponents = 2);
  fWater->AddElement(H, natoms = 2);
  fWater->AddElement(O, natoms = 1);

  // Metal Lead
  density = 11.34 * g / cm3;
  fMetalLead = new G4Material(name = "MetalLead", density, ncomponents = 1);
  fMetalLead->AddElement(Pb, 1.0);

  //  Teflon
  density = 2.165 * g / cm3;
  fTeflon = new G4Material(name = "Teflon", density = 2.165 * g / cm3, ncomponents = 2);
  fTeflon->AddElement(Fl, fractionmass = 0.76);
  fTeflon->AddElement(C, fractionmass = 0.24);

  //  Rock: hep-ex/0312050v2 Wulandari et al.
  density = 2.71 * g / cm3;  //
  fRock = new G4Material(name = "Rock", density, ncomponents = 7);
  fRock->AddElement(C, fractionmass = 0.1188);
  fRock->AddElement(O, fractionmass = 0.4892);  // <-- here I added 1.01% to obtain a fractionmass = 1
  fRock->AddElement(Mg, fractionmass = 0.0558);
  fRock->AddElement(Al, fractionmass = 0.0103);
  fRock->AddElement(Si, fractionmass = 0.0127);
  fRock->AddElement(K, fractionmass = 0.0103);
  fRock->AddElement(Ca, fractionmass = 0.3029);

  //  Norite (SNOlab rock) Density from Thomas Cement AB (Thomas Concrete Group)
  density = 3.15 * g / cm3;
  fNorite = new G4Material(name = "Norite", density, ncomponents = 11, kStateSolid);
  fNorite->AddElement(H, fractionmass = 0.15 * perCent);
  fNorite->AddElement(C, fractionmass = 0.04 * perCent);
  fNorite->AddElement(O, fractionmass = 46.0 * perCent);
  fNorite->AddElement(Na, fractionmass = 2.2 * perCent);
  fNorite->AddElement(Mn, fractionmass = 3.4 * perCent);
  fNorite->AddElement(Al, fractionmass = 9.0 * perCent);
  fNorite->AddElement(Si, fractionmass = 26.2 * perCent);
  fNorite->AddElement(K, fractionmass = 1.2 * perCent);
  fNorite->AddElement(Ca, fractionmass = 5.2 * perCent);
  //fNorite->AddElement(Mn, fractionmass = 0.1 * perCent);
  fNorite->AddElement(Fe, fractionmass = 6.2 * perCent);
  fNorite->AddElement(Ti, fractionmass = 0.5 * perCent);

  // High Density PolyEthylene
  density = 0.942 * g / cm3;
  fHDPE = new G4Material(name = "HDPE", density, ncomponents = 2);
  fHDPE->AddElement(H, 0.135);
  fHDPE->AddElement(C, 0.865);

  // Fake Styrofoam
  density = 0.0996996029 * g / cm3;
  fStyrofoam = new G4Material(name = "Styrofoam", density, ncomponents = 2);
  fStyrofoam->AddElement(H, 0.135);
  fStyrofoam->AddElement(C, 0.865);

  // Metal Copper
  density = 8.96 * g / cm3;
  fMetalCopper = new G4Material(name = "MetalCopper", density, ncomponents = 1);
  fMetalCopper->AddElement(Cu, 1.0);

  // Metal Copper
  density = 8.96 * g / cm3;
  fMetalCopperCryo = new G4Material(name = "MetalCopperCryo", density, ncomponents = 1);
  fMetalCopperCryo->AddElement(Cu, 1.0);

  // Metal Titanium
  density = 4.506 * g / cm3;
  fMetalTitanium = new G4Material(name = "MetalTitanium", density, ncomponents = 1);
  fMetalTitanium->AddElement(Ti, 1.0);
  // Ship Steel
  density = 7.85 * g / cm3;
  fShipSteel = new G4Material(name = "ShipSteel", density, ncomponents = 2);
  fShipSteel->AddElement(C, fractionmass = 0.002);
  fShipSteel->AddElement(Fe, fractionmass = 0.998);

  // TPB (Tetra-Phenyl Butadiene)
  // Density copied from PseudoCumene
  density = 0.9 * g / cm3;
  fTPB = new G4Material(name = "TPB", density, ncomponents = 2, kStateSolid);
  fTPB->AddElement(C, 28);
  fTPB->AddElement(H, 22);

  // ThreeMFoil (Polyethylene terephthalate)
  // http://multimedia.3m.com/mws/mediawebserver?mwsId=66666UuZjcFSLXTtM8T6oXTVEVuQEcuZgVs6EVs6E666666--
  // Density and composition from wiki
  density = 1.4 * g / cm3;
  fThreeMFoil = new G4Material(name = "ThreeMFoil", density, ncomponents = 3);
  fThreeMFoil->AddElement(C, 10);
  fThreeMFoil->AddElement(H, 8);
  fThreeMFoil->AddElement(O, 4);

  // Fused Silica (same composition as quartz)
  // Density from Wikipedia
  density = 2.203 * g / cm3;
  fFusedSilica = new G4Material(name = "FusedSilica", density, ncomponents = 2, kStateSolid);
  fFusedSilica->AddElement(Si, 1);
  fFusedSilica->AddElement(O, 2);

  // Sodium Iodide
  density = 3.67 * g / cm3;
  fSodiumIodide = new G4Material(name = "SodiumIodide", density, ncomponents = 2);
  fSodiumIodide->AddElement(Na, natoms = 1);
  fSodiumIodide->AddElement(I, natoms = 1);

  // ITO
  // Indium corporation product data sheet
  density = 7.2 * g / cm3;
  fITO = new G4Material(name = "ITO", density, ncomponents = 3);
  fITO->AddElement(In, fractionmass = 0.7444);
  fITO->AddElement(Sn, fractionmass = 0.0788);
  fITO->AddElement(O, fractionmass = 0.1768);

  // Liquid Nitrogen
  density = 0.807 * g / cm3;
  fLiquidNitrogen = new G4Material(name = "LiquidNitrogen", density, ncomponents = 1, kStateLiquid, temperature = 77 * kelvin);
  fLiquidNitrogen->AddElement(N, 1.0);

  // Gaseous Nitrogen
  density = 0.000125 * g / cm3;
  fGaseousNitrogen = new G4Material(name = "GaseousNitrogen", density, ncomponents = 1, kStateGas, temperature = 298 * kelvin);
  fGaseousNitrogen->AddElement(N, 1.0);

  // Kovar
  density = 8.0 * g / cm3;
  fKovar = new G4Material(name = "Kovar", density, ncomponents = 3);
  fKovar->AddElement(Fe, fractionmass = 0.54);
  fKovar->AddElement(Ni, fractionmass = 0.29);
  fKovar->AddElement(Co, fractionmass = 0.17);

  // Bialkali (Wiki: Photocathode)
  density = 3.39 * g / cm3;
  fBialkali = new G4Material(name = "Bialkali", density, ncomponents = 3);
  fBialkali->AddElement(Sb, fractionmass = 0.34);
  fBialkali->AddElement(Rb, fractionmass = 0.33);
  fBialkali->AddElement(Cs, fractionmass = 0.33);

  // Kapton (Wiki: Kapton)
  density = 1.0 * g / cm3;
  fKapton = new G4Material(name = "Kapton", density, ncomponents = 4);
  fKapton->AddElement(C, 22);
  fKapton->AddElement(N, 2);
  fKapton->AddElement(O, 5);
  fKapton->AddElement(H, 10);

  // Concrete: hep-ex/0312050v2 Wulandari et al.
  density = 2.4 * g / cm3;  //
  fConcrete = new G4Material(name = "Concrete", density, ncomponents = 13);
  fConcrete->AddElement(H, fractionmass = 0.0089);
  fConcrete->AddElement(C, fractionmass = 0.0799);
  fConcrete->AddElement(O, fractionmass = 0.4971);  // <-- here I add 1.28% to obtain a fractionmass = 1
  fConcrete->AddElement(Na, fractionmass = 0.0006);
  fConcrete->AddElement(Mg, fractionmass = 0.0085);
  fConcrete->AddElement(Al, fractionmass = 0.0009);
  fConcrete->AddElement(Si, fractionmass = 0.0386);
  fConcrete->AddElement(P, fractionmass = 0.0004);
  fConcrete->AddElement(S, fractionmass = 0.0016);
  fConcrete->AddElement(K, fractionmass = 0.0054);
  fConcrete->AddElement(Ca, fractionmass = 0.3534);  // <-- here I add 1.28% to obtain a fractionmass = 1
  fConcrete->AddElement(Ti, fractionmass = 0.0004);
  fConcrete->AddElement(Fe, fractionmass = 0.0043);

  // Borex materials

  // Aluminum
  density = 2.700 * g / cm3;
  fAluminum = new G4Material(name = "Aluminum", density, ncomponents = 1);
  fAluminum->AddElement(Al, fractionmass = 1.);

  // Glass
  density = 2.2 * g / cm3;  // should be corrected
  fGlass = new G4Material(name = "Glass", density, ncomponents = 2);
  fGlass->AddElement(Si, natoms = 1);
  fGlass->AddElement(O, natoms = 2);

  // Nylon
  density = 1.1 * g / cm3;  // should be corrected
  fNylon = new G4Material(name = "Nylon", density, ncomponents = 3);
  fNylon->AddElement(H, fractionmass = 0.08);
  fNylon->AddElement(C, fractionmass = 0.60);
  fNylon->AddElement(O, fractionmass = 0.32);

  // NylonL
  fNylonL = new G4Material(name = "NylonL", density, ncomponents = 3);
  fNylonL->AddElement(H, fractionmass = 0.08);
  fNylonL->AddElement(C, fractionmass = 0.60);
  fNylonL->AddElement(O, fractionmass = 0.32);

  // Bialkali
  density = 4.000 * g / cm3;  // should be corrected, but not important
  fOldBialkali = new G4Material(name = "Bialkali", density, ncomponents = 2);
  fOldBialkali->AddElement(K, fractionmass = 0.5);
  fOldBialkali->AddElement(Cs, fractionmass = 0.5);

  // Paint
  density = 4.000 * g / cm3;  // should be corrected, but not important
  fPaint = new G4Material(name = "Paint", density, ncomponents = 2);
  fPaint->AddElement(K, fractionmass = 0.5);
  fPaint->AddElement(Cs, fractionmass = 0.5);

  // fBorexScintillator
  density = 0.8774 * g / cm3;
  fBorexScintillator = new G4Material(name = "BorexScintillator", density, ncomponents = 2);
  fBorexScintillator->AddMaterial(fPC, fractionmass = 99.83 * perCent);
  fBorexScintillator->AddMaterial(fPPO, fractionmass = 0.17 * perCent);

  // fBorexScintillator1
  density = 0.8774 * g / cm3;
  fBorexScintillator1 = new G4Material(name = "BorexScintillator1", density, ncomponents = 2);
  fBorexScintillator1->AddMaterial(fPC, fractionmass = 99.83 * perCent);
  fBorexScintillator1->AddMaterial(fPPO, fractionmass = 0.17 * perCent);

  // fBorexScintillator2
  density = 0.8774 * g / cm3;
  fBorexScintillator2 = new G4Material(name = "BorexScintillator2", density, ncomponents = 2);
  fBorexScintillator2->AddMaterial(fPC, fractionmass = 99.83 * perCent);
  fBorexScintillator2->AddMaterial(fPPO, fractionmass = 0.17 * perCent);

  // DMPBuffer
  // Density of the buffer reported by Frank Calaprice is by 0.2 kg/m^3 higher
  // than the scintillator 2.5 g/l DMP - October 2009
  density = 0.8776 * g / cm3;
  fDMPbuffer = new G4Material(name = "DMPbuffer", density, ncomponents = 2);
  fDMPbuffer->AddMaterial(fPC, fractionmass = 99.66 * perCent);
  fDMPbuffer->AddMaterial(fDMP, fractionmass = 0.34 * perCent);

  // Quartz
  density = 2.200 * g / cm3;
  fQuartz = new G4Material(name = "Quartz", density, ncomponents = 2);
  fQuartz->AddElement(Si, 1);
  fQuartz->AddElement(O, 2);

  //  Beryllium
  density = 1.848 * g / cm3;
  fBe = new G4Material(name = "Beryllium", density, ncomponents = 1);
  fBe->AddElement(Be, fractionmass = 1.);

  // DeLRin
  density = 1.47 * g / cm3;
  fDerlin = new G4Material(name = "Derlin", density, ncomponents = 2);
  fDerlin->AddElement(H, fractionmass = 0.067);
  fDerlin->AddElement(C, fractionmass = 0.933);

  // BlackHole
  density = 100 * g / cm3;
  fBlackHole = new G4Material(name = "BlackHole", density, ncomponents = 1);
  fBlackHole->AddElement(C, 2);

  // Gadolinium Oxide
  density = 7.41 * g / cm3;
  fGdOxide = new G4Material(name = "fGdOxide", density, ncomponents = 2);
  fGdOxide->AddElement(Gd, natoms = 2);
  fGdOxide->AddElement(O, natoms = 3);
  // fGdOxide->AddElement(S, natoms=3);

  // GdScintillator
  density = 0.8776 * g / cm3;
  fGdScintillator = new G4Material(name = "GdScintillator", density, ncomponents = 2);
  // fGdScintillator->AddMaterial(fBorexScintillator,
  // fractionmass=(877.6-1.)/877.6); fGdScintillator->AddMaterial(fGdOxide,
  // fractionmass=1./877.6 );
  //  fraction = 0.0014 (1 g/l) for 35us capture time in DS20k
  fGdScintillator->AddMaterial(fGdOxide, fractionmass = fraction);
  fGdScintillator->AddMaterial(fBorexScintillator, fractionmass = 1 - fraction);

  // GdWater
  density = 1.0 * g / cm3;
  fGdWater = new G4Material(name = "GdWater", density, ncomponents = 2);
  fGdWater->AddMaterial(fWater, fractionmass = (1000. - 2.) / 1000);
  fGdWater->AddMaterial(fGdOxide, fractionmass = 2. / 1000.);

  // Metal Silicon
  density = 2.33 * g / cm3;
  fMetalSilicon = new G4Material(name = "MetalSilicon", density, ncomponents = 1);
  fMetalSilicon->AddElement(Si, 1.0);

  // BoroSilicate glass
  density = 2.23 * g / cm3;
  fBoroSiliGlass = new G4Material(name = "BoroSilicateGlass", density, ncomponents = 7);
  fBoroSiliGlass->AddElement(Li, 0.0047);
  fBoroSiliGlass->AddElement(NewB, .0529);
  fBoroSiliGlass->AddElement(Na, .046);
  fBoroSiliGlass->AddElement(Si, .3243);
  fBoroSiliGlass->AddElement(Al, .0318);
  fBoroSiliGlass->AddElement(O, .5215);
  fBoroSiliGlass->AddElement(Ba, .0179);

  // Fake Photocatode - - - As Fused Silica
  density = 2.203 * g / cm3;
  fFakePhotocathode = new G4Material(name = "FakePhotocathode", density, ncomponents = 2, kStateSolid);
  fFakePhotocathode->AddElement(Si, 1);
  fFakePhotocathode->AddElement(O, 2);

  density = 2.07 * g / cm3;
  G4Material* ClLi = new G4Material(name = "ClLi", density, ncomponents = 2);
  ClLi->AddElement(EnrLi, natoms = 1);
  ClLi->AddElement(natCl, natoms = 1);

  density = 2.07 * g / cm3;
  G4Material* NatClLi = new G4Material(name = "NatClLi", density, ncomponents = 2);
  NatClLi->AddElement(natLi, natoms = 1);
  NatClLi->AddElement(natCl, natoms = 1);

  density = 0.95 * g / cm3;
  G4Material* OrtoCarb = new G4Material(name = "OrtoCarb", density, ncomponents = 3);
  OrtoCarb->AddElement(C, natoms = 2);
  OrtoCarb->AddElement(H, natoms = 12);
  OrtoCarb->AddElement(NewB, natoms = 10);

  density = 0.877 * g / cm3;  // same as PC
  fLi6Scintillator = new G4Material(name = "Li6Scintillator", density, ncomponents = 2);
  double ClLimassFract = fraction;  // 7.33e-2 ; // 0.010 * 42.39 / 6.94; //not 0.0015 as Prospect
                                    // because of capture time Changed by B.Bottino 17March
  // fraction =  7.33e-3 for 40us capture time in DS20k
  fLi6Scintillator->AddMaterial(ClLi, fractionmass = ClLimassFract);
  fLi6Scintillator->AddMaterial(fPC, fractionmass = 1 - ClLimassFract);

  density = 0.877 * g / cm3;  // same as PC
  fNatLi6Scintillator = new G4Material(name = "NatLi6Scintillator", density, ncomponents = 2);
  double NatClLimassFract = fraction;  // 0.011 * 42.39 / 6.94; // Changed by B.Bottino 16March
  fNatLi6Scintillator->AddMaterial(NatClLi, fractionmass = NatClLimassFract);
  fNatLi6Scintillator->AddMaterial(fPC, fractionmass = 1 - ClLimassFract);

  density = 0.877 * g / cm3;  // same as PC
  fOrtoCarbScintillator = new G4Material(name = "OrtoCarbScintillator", density, ncomponents = 2);
  double OrtoCarbmassFract = fraction;  // 0.003; // Changed by B.Bottino from
  // fraction =  0.085 for 3us capture time in DS20k
  fOrtoCarbScintillator->AddMaterial(OrtoCarb, fractionmass = OrtoCarbmassFract);
  fOrtoCarbScintillator->AddMaterial(fPC, fractionmass = 1 - OrtoCarbmassFract);

  // GdAcrylic
  G4double GdConcentration = DSStorage::Get()->GetDS20kGdConcentration();
  G4double GdAcrylicConcentration = GdConcentration * (1. / 0.87);
  density = 1.2 * g / cm3;
  fGdAcrylic = new G4Material(name = "GdAcrylic", density, ncomponents = 2, kStateSolid);
  fGdAcrylic->AddMaterial(fAcrylic, fractionmass = 1. - GdAcrylicConcentration);
  fGdAcrylic->AddMaterial(fGdOxide, fractionmass = GdAcrylicConcentration);

  // Veto - Non Scintillating Liquid Argon
  density = 1.40 * g / cm3;
  temperature = 87 * kelvin;
  fIVLiquidArgon = new G4Material(name = "IVLiquidArgon", density, ncomponents = 1, kStateLiquid, temperature);
  fIVLiquidArgon->AddElement(Ar, 1.0);

  // Naphthalene
  density = 1.0 * g / cm3;
  G4Material* Naphthalene = new G4Material(name = "Naphthalene", density, ncomponents = 2);
  Naphthalene->AddElement(C, 10);
  Naphthalene->AddElement(H, 8);

  // Altustipe - the naphthalene can be removed
  density = 1.2 * g / cm3;
  temperature = 290 * kelvin;
  fAltustipe = new G4Material(name = "Altustipe", density, ncomponents = 2, kStateSolid, temperature);
  fAltustipe->AddMaterial(fAcrylic, fractionmass = (1 - fraction));
  fAltustipe->AddMaterial(Naphthalene, fractionmass = fraction);
  // fAltustipe->AddMaterial(fTMB,fractionmass =0.12 );
  // fAltustipe->AddMaterial(fGdOxide,fractionmass =0.12 );

  // BC454
  density = 1.026 * g / cm3;
  temperature = 83 * kelvin;
  fBC454 = new G4Material(name = "BC454", density, ncomponents = 3, kStateSolid, temperature);
  fBC454->AddElement(C, 32);
  fBC454->AddElement(H, 37);
  fBC454->AddElement(NewB, 2);

  // PlyWood
  density = 0.7 * g / cm3;
  fPlyWood = new G4Material(name = "PlyWood", density, ncomponents = 4, kStateSolid);
  fPlyWood->AddElement(C, fractionmass = 0.496);
  fPlyWood->AddElement(O, fractionmass = 0.44);
  fPlyWood->AddElement(H, fractionmass = 0.063);
  fPlyWood->AddElement(N, fractionmass = 0.001);

  // Insulating Polyurethane Foam
  density = 0.09 * g / cm3;
  fInsulatingFoam = new G4Material(name = "InsulatingFoam", density, ncomponents = 4, kStateSolid);
  fInsulatingFoam->AddElement(C, fractionmass = 0.4286);
  fInsulatingFoam->AddElement(O, fractionmass = 0.3265);
  fInsulatingFoam->AddElement(H, fractionmass = 0.1428);
  fInsulatingFoam->AddElement(N, fractionmass = 0.1021);

  // Outer Veto - Non Scintillating Liquid Argon
  density = 1.40 * g / cm3;
  temperature = 87 * kelvin;
  fPScintVetoLiquidArgon = new G4Material(name = "LiquidArgonPlasticVeto", density, ncomponents = 1, kStateLiquid, temperature);
  fPScintVetoLiquidArgon->AddElement(Ar, 1.0);

  // Acrylic--DART
  density = 1.19 * g / cm3;
  fAcrylicDART = new G4Material(name = "AcrylicDART", density, ncomponents = 3, kStateSolid);
  fAcrylicDART->AddElement(C, 5);
  fAcrylicDART->AddElement(O, 2);
  fAcrylicDART->AddElement(H, 8);

  // Acrylic -- Final candidate for DS20k veto plastic scintillator.
  // use std acrylic for the moment
  density = 1.19 * g / cm3;
  fDS20kPlasticScintillator = new G4Material(name = "ReflectorAcrylic", density, ncomponents = 3, kStateSolid);
  fDS20kPlasticScintillator->AddElement(C, 5);
  fDS20kPlasticScintillator->AddElement(O, 2);
  fDS20kPlasticScintillator->AddElement(H, 8);

  // Methane
  density = 1.4 * g / cm3;
  fMethane = new G4Material(name = "Methane", density, ncomponents = 2, kStateLiquid);
  fMethane->AddElement(C, 1);
  fMethane->AddElement(H, 4);

  // Veto Mixture
  density = 1.4 * g / cm3;
  double CH4_Xe_fraction = 0;
  fOVLiquidArgon = new G4Material(name = "OVLiquidArgon", density, ncomponents = CH4_Xe_fraction > 1E-3 ? 3 : 1, kStateLiquid);
  if (CH4_Xe_fraction > 0.001) {
    fOVLiquidArgon->AddMaterial(fMethane, fractionmass = fraction * (1. - 0.00116));
    fOVLiquidArgon->AddMaterial(fLiquidXenon, fractionmass = 0.00116);  // mass
                                                                      // frac
  }
  fOVLiquidArgon->AddMaterial(fLiquidArgon, fractionmass = CH4_Xe_fraction > 1E-3 ? (1 - fraction) * (1. - 0.00116) : 1);

  // Xe doped Liquid Argon
  density = 1.40 * g / cm3;
  temperature = 87 * kelvin;
  fXeLiquidArgon = new G4Material(name = "XeLiquidArgon", density, ncomponents = 2, kStateLiquid, temperature);
  fXeLiquidArgon->AddMaterial(fLiquidArgon, fractionmass = 1.0 - fraction);
  fXeLiquidArgon->AddMaterial(fLiquidXenon, fractionmass = fraction);
  G4cout << fXeLiquidArgon << endl;

  // invar
  density = 8.1 * g / cm3;
  temperature = 87 * kelvin;
  fInvar = new G4Material(name = "Invar", density, ncomponents = 2, kStateSolid, temperature);
  fInvar->AddElement(Fe, fractionmass = 0.64);
  fInvar->AddElement(Ni, fractionmass = 0.36);

  // arlon
  density = 2.3 * g / cm3;
  temperature = 87 * kelvin;
  fArlon = new G4Material(name = "Arlon", density, ncomponents = 1, kStateSolid, temperature);
  fArlon->AddMaterial(fAcrylic, fractionmass = 1);

  // tungsten
  density = 19.25 * g / cm3;
  temperature = 87 * kelvin;
  fMetalTungsten = new G4Material(name = "MetalTungsten", density, ncomponents = 1, kStateSolid, temperature);
  fMetalTungsten->AddElement(W, 1.0);

  // tantalum
  density = 16.69 * g / cm3;
  temperature = 87 * kelvin;
  fMetalTantalum = new G4Material(name = "MetalTantalum", density, ncomponents = 1, kStateSolid, temperature);
  fMetalTantalum->AddElement(Ta, 1.0);

  // PEN
  density = 1.36 * g / cm3;
  fPEN = new G4Material(name = "PEN", density, ncomponents = 3, kStateSolid);
  fPEN->AddElement(C, 14);
  fPEN->AddElement(H, 10);
  fPEN->AddElement(O, 4);

  //Insulating Polyurethane Foam Low Density for the roof
  density      = 0.04*g/cm3;
  fInsulatingFoamRoof = new G4Material(name="InsulatingFoamRoof",density,ncomponents=4,kStateSolid);
  fInsulatingFoamRoof->AddElement(C, fractionmass=0.4286);
  fInsulatingFoamRoof->AddElement(O, fractionmass=0.3265);
  fInsulatingFoamRoof->AddElement(H, fractionmass=0.1428);
  fInsulatingFoamRoof->AddElement(N, fractionmass=0.1021);

  // Barium fluoride scintillator (no optical properties yet)
  //   [0]=generic material, [1]=for ARIS-ER 511 keV detector, [1]=for ARIS-ER 1270 keV detector
  density = 4.88*g/cm3;
  for (int i=0; i<3; i++) {
    string matname = "BariumFluoride";
    if (i != 0)  matname += "_" + std::to_string(i);
    fBariumFluorideList[i] = new G4Material(matname, density, ncomponents=2, kStateSolid);
    fBariumFluorideList[i]->AddElement(Ba, 1);
    fBariumFluorideList[i]->AddElement(Fl, 2);
  }

  // Ceramic
  density = 4*g/cm3;
  fCeramic = new G4Material(name="Ceramic", density, ncomponents=2, kStateSolid);
  fCeramic->AddElement(Al, 2);
  fCeramic->AddElement(O, 3);

  // Carbon
  density = 1.55 * g / cm3;
  fCarbon = new G4Material(name="Carbon", density, ncomponents=1);
  fCarbon->AddElement(C, 1.0);

  // Germanium
  density = 5.5 * g / cm3;
  fGermanium = new G4Material(name="Germanium", density, ncomponents=1);
  fGermanium->AddElement(Ge, 1.0);

  // ESR
  density = 1.4 * g / cm3; // as if PET
  fESR = new G4Material(name = "ESR", density, ncomponents = 3);
  fESR->AddElement(H, 8);
  fESR->AddElement(C, 10);
  fESR->AddElement(O, 4);

  // arlon 55NT prepreg for PCB
  density = 1.3 * g / cm3;
  temperature = 87 * kelvin;
  fArlonPrepreg = new G4Material(name = "ArlonPrepreg", density, ncomponents = 4, kStateSolid, temperature);
  fArlonPrepreg->AddElement(H,fractionmass = 0.065);
  fArlonPrepreg->AddElement(C,fractionmass = 0.668);
  fArlonPrepreg->AddElement(O,fractionmass = 0.209);
  fArlonPrepreg->AddElement(N,fractionmass = 0.058);

  // Liquid Argon above Grid
  density = 1.40 * g / cm3;
  temperature = 87 * kelvin;
  fLArAboveGrid = new G4Material(name = "LArAboveGrid", density, ncomponents = 1, kStateLiquid, temperature);
  fLArAboveGrid->AddElement(Ar, 1.0);
  
}

G4Material* DSMaterial::GetNDReDScintillator() {
  // as of Borexino scintillator

  G4int ncomponents;
  G4double density, fractionmass;

  density = 0.8774 * g / cm3;
  int idBorex = fNDReDScintillatorList.size();

  string name = "BorexScintillator_";
  ostringstream oss_idBorex;
  oss_idBorex << idBorex;

  name += oss_idBorex.str();

  G4Material* BorexScintillator = new G4Material(name, density, ncomponents = 2);
  BorexScintillator->AddMaterial(fPC, fractionmass = 99.83 * perCent);
  BorexScintillator->AddMaterial(fPPO, fractionmass = 0.17 * perCent);

  fNDReDScintillatorList.push_back(BorexScintillator);

  DSLog(development) << "----------------------------------------" << endlog;
  DSLog(development) << "           New Material indexes         " << endlog;
  DSLog(development) << "----------------------------------------" << endlog;
  DSLog(routine) << BorexScintillator->GetName() << ": " << (int)BorexScintillator->GetIndex() << endlog;
  DSLog(development) << "----------------------------------------" << endlog;

  int borex_info = ((int)fNDReDScintillatorList.size()) << 8;
  borex_info = borex_info | (int)BorexScintillator->GetIndex();

  DSStorage::Get()->SetBoronScintillatorIndex(borex_info);
  DSEventHandler::Get()->SetScintillatorIndex(borex_info);

  return BorexScintillator;
}

void DSMaterial::DefineProperties() {
  //---------------------------------------------------------------------------------
  // Scintillator
  //---------------------------------------------------------------------------------

  const G4double hc = h_Planck * c_light / eV / nm;  // in eV*nm, ~ 1239.84172

  // local variables
  G4double myvalue, myene, myene2, mywl;
  G4int dim, dim1, dim2;
  const G4int BScint_NUMENTRIES = 401;
  const G4int PPOScint_NUMENTRIES = 401;
  G4double BScint_Energy[BScint_NUMENTRIES];
  G4double BScint_SCINT[BScint_NUMENTRIES];
  G4double ppo_Energy[PPOScint_NUMENTRIES];
  G4double ppo_SCINT[BScint_NUMENTRIES];
  // G4double BScint_const_rind = 1.4; // Choosing 1.4, since it is between 1.3
  // (TMB) and 1.5 (PC)

  // Fill in values for SS_RIND now, since it needs the same wavelengths as the
  // scint
  //G4double StainlessSteel_RIND[BScint_NUMENTRIES];
  //G4double StainlessSteel_const_rind = 1.4; // 1.56;

  // Load BScint properties, as measured by Aldo
  const char bscintFileName[] = "../data/detector/bscintSpectrum.dat";
  std::ifstream bscint_data(bscintFileName);
  G4int il = 0;
  if (!bscint_data.is_open()) DSLog(fatal) << "ERROR: COULD NOT OPEN BSCINT SCINT FILE " << bscintFileName << endlog;
  while (bscint_data >> mywl >> myvalue) {
    const int j = BScint_NUMENTRIES - 1 - il;
    // Convert wavelength from nm to energy
    BScint_Energy[j] = (hc / mywl) * eV;
    BScint_SCINT[j] = myvalue;
    // DSLog(debugging) << "BSCINT ENERGY = " << BScint_Energy[j]/eV << "
    // eV\tSCINT = " << BScint_SCINT[j] << endlog;
    il++;
  }
  bscint_data.close();

  // Load PPO Properties
  const char ppoFileName[] = "../data/detector/ppoEmissionSpectrum.dat";
  std::ifstream ppo_data(ppoFileName);
  il = 0;
  if (!ppo_data.is_open()) DSLog(fatal) << "ERROR: COULD NOT OPEN PPO DATA FILE : " << ppoFileName << endlog;
  while (ppo_data >> mywl >> myvalue) {
    const int j = BScint_NUMENTRIES - 1 - il;
    // Convert wavelength from nm to energy
    ppo_Energy[j] = (hc / mywl) * eV;
    ppo_SCINT[j] = myvalue;
    //StainlessSteel_RIND[j] = StainlessSteel_const_rind;
    il++;
  }
  ppo_data.close();

  // LS Absorption length
  // Combining PC and TMB absorption lenght
  // See e.g Knoll
  // 1 / lambda_tot = Sum_i ( f_i /lambda_i )
  // where f_i is the volume fraction

  // Hereby I am assuming that all attenuation length is absorption length
  // This is a very conservative assumption on the light yield
  // e.g. Rayleigh scattering becomes relevant for lambda > 400 nm
  // see NIM A 400 (2000) 360-371

  vector<G4double> myLSAbsEnergy;
  vector<G4double> myLSAbsLength;
  myLSAbsEnergy.reserve(1000);
  myLSAbsLength.reserve(1000);
  const vector<G4double> myPCAttEnergy = DSParameters::Get()->GetPCAttEnergy();
  const G4double TMBfrac = DSStorage::Get()->GetTMBfraction();
  const G4double PCfrac = 1. - TMBfrac;
  for (unsigned i = 0; i < myPCAttEnergy.size(); i++) {
    const G4double lpc = DSParameters::Get()->GetPCAttLength()[i];
    const G4double ltmb = DSParameters::Get()->GetTMBAttLength()[i];
    const G4double lambda = 1. / (PCfrac / lpc + TMBfrac / ltmb);
    myLSAbsEnergy.push_back(myPCAttEnergy[i]);
    myLSAbsLength.push_back(lambda);
  }

  // Enter Material Properties into table
  G4MaterialPropertiesTable* BScint_MPT = new G4MaterialPropertiesTable();
  BScint_MPT->AddProperty("FASTCOMPONENT", BScint_Energy, BScint_SCINT, BScint_NUMENTRIES, true);
  BScint_MPT->AddProperty("SLOWCOMPONENT", BScint_Energy, BScint_SCINT, BScint_NUMENTRIES, true);
  BScint_MPT->AddProperty("RINDEX", &(DSParameters::Get()->GetPCRefrEnergy())[0], &(DSParameters::Get()->GetPCRefrIndex())[0], DSParameters::Get()->GetPCRefrEnergy().size(), true);
  BScint_MPT->AddProperty("ABSLENGTH", &myLSAbsEnergy[0], &myLSAbsLength[0], myLSAbsEnergy.size(), true);
  const G4double scintYield = (DSParameters::Get()->GetLSYield()) * (DSParameters::Get()->GetVPmtMaxQe_adjusted());
  BScint_MPT->AddConstProperty("SCINTILLATIONYIELD", scintYield);
  BScint_MPT->AddConstProperty("YIELDRATIO", DSParameters::Get()->GetLSYieldRatio(), true);
  BScint_MPT->AddConstProperty("RESOLUTIONSCALE", DSParameters::Get()->GetLSResolutionScale(), true);
  BScint_MPT->AddConstProperty("FASTTIMECONSTANT", DSParameters::Get()->GetLSFastTime(), true);
  BScint_MPT->AddConstProperty("SLOWTIMECONSTANT", DSParameters::Get()->GetLSSlowTime(), true);
  BScint_MPT->AddProperty("WLSCOMPONENT", ppo_Energy, ppo_SCINT, BScint_NUMENTRIES);
  BScint_MPT->AddProperty("WLSABSLENGTH", &(DSParameters::Get()->GetPPOAttEnergy())[0], &(DSParameters::Get()->GetPPOAttLength())[0], DSParameters::Get()->GetPPOAttEnergy().size());
  BScint_MPT->AddConstProperty("WLSTIMECONSTANT", DSParameters::Get()->GetLSWLSTime());
  BScint_MPT->AddConstProperty("WLSEFFICIENCY", DSParameters::Get()->GetLSWLSEfficiency(), true);
  fBoronScintillator->SetMaterialPropertiesTable(BScint_MPT);
  fGdScintillator->SetMaterialPropertiesTable(BScint_MPT);
  // TODO for the moment, use boronscintillator properties
  fLi6Scintillator->SetMaterialPropertiesTable(BScint_MPT);

  const G4double pcBirk = DSParameters::Get()->GetLSBirksBeta();
  fBoronScintillator->GetIonisation()->SetBirksConstant(pcBirk);
  fGdScintillator->GetIonisation()->SetBirksConstant(pcBirk);
  // TODO for the moment, use boronscintillator properties
  fLi6Scintillator->GetIonisation()->SetBirksConstant(pcBirk);

  //---------------------------------------------------------------------------------
  // Water and GdWater
  //---------------------------------------------------------------------------------

  G4double myWEnergy[2] = {0.1 * eV, 15. * eV};
  G4double myWRI[2] = {1.33, 1.33};
  G4MaterialPropertiesTable* water_MPT = new G4MaterialPropertiesTable();
  water_MPT->AddProperty("RINDEX", myWEnergy, myWRI, 2);
  fWater->SetMaterialPropertiesTable(water_MPT);
  fGdWater->SetMaterialPropertiesTable(water_MPT);
  // Temporary attribute water properties to air
  fAir->SetMaterialPropertiesTable(water_MPT);

  //---------------------------------------------------------------------------------
  // Stainless Steel
  //---------------------------------------------------------------------------------
  G4double mySSEnergy[2] = {0.1 * eV, 10 * eV};
  G4double mySSABSL[2] = {0.1 * nm, 0.1 * nm};
  // G4double mySSRI[2]     = {1.47, 1.47};
  G4MaterialPropertiesTable* ss_MPT = new G4MaterialPropertiesTable();
  // ss_MPT->AddProperty("RINDEX",mySSEnergy,mySSRI,2)->SetSpline(true);
  ss_MPT->AddProperty("ABSLENGTH", mySSEnergy, mySSABSL, 2);
  fStainlessSteel->SetMaterialPropertiesTable(ss_MPT);

  //---------------------------------------------------------------------------------
  // Gaseous Nitrogen
  //---------------------------------------------------------------------------------
  G4double myGN2Energy[2] = {0.1 * eV, 10. * eV};
  G4double myGN2ABSL[2] = {100. * m, 100. * nm};
  G4double myGN2RI[2] = {1., 1.};
  G4MaterialPropertiesTable* gN2_MPT = new G4MaterialPropertiesTable();
  gN2_MPT->AddProperty("RINDEX", myGN2Energy, myGN2RI, 2);
  gN2_MPT->AddProperty("ABSLENGTH", myGN2Energy, myGN2ABSL, 2);
  fGaseousNitrogen->SetMaterialPropertiesTable(gN2_MPT);

  //---------------------------------------------------------------------------------
  // Fused Silica
  //---------------------------------------------------------------------------------
  G4MaterialPropertiesTable* myFusedSilica = new G4MaterialPropertiesTable();
  G4double myFSEnergy[4], myFSRI[4], myFSAbs[4], myRLength[4];
  myFSEnergy[0] = 1.0 * eV;
  myFSEnergy[1] = 8.0 * eV;
  myFSRI[0] = DSParameters::Get()->GetFusedSilicaVisRind();
  myFSAbs[0] = DSParameters::Get()->GetFusedSilicaVisAbs() * m;
  myRLength[0] = DSParameters::Get()->GetFSilicaRaylVisLength() * m;
  myFSRI[1] = myFSRI[0];
  myFSAbs[1] = myFSAbs[0];
  myRLength[1] = myRLength[0];
  myFSEnergy[2] = 8.3 * eV;
  myFSEnergy[3] = 20.0 * eV;
  myFSRI[2] = DSParameters::Get()->GetFusedSilicaUVRind();
  myFSAbs[2] = DSParameters::Get()->GetFusedSilicaUVAbs() * m;
  ;
  myRLength[2] = DSParameters::Get()->GetFSilicaRaylUVLength() * m;
  myFSRI[3] = myFSRI[2];
  myFSAbs[3] = myFSAbs[2];
  myRLength[3] = myRLength[2];
  myFusedSilica->AddProperty("RINDEX", myFSEnergy, myFSRI, 4);
  myFusedSilica->AddProperty("ABSLENGTH", myFSEnergy, myFSAbs, 4);
  myFusedSilica->AddProperty("RAYLEIGH", myFSEnergy, myRLength, 4);
  fFusedSilica->SetMaterialPropertiesTable(myFusedSilica);

  //---------------------------------------------------------------------------------
  // Teflon
  //---------------------------------------------------------------------------------

  G4MaterialPropertiesTable* myTeflon = new G4MaterialPropertiesTable;
  G4double myTeflonEnergy[3], myTeflonAbs[3];
  myTeflonEnergy[0] = 1.0 * eV;
  //myTeflonRI[0] = 1.4;
  myTeflonAbs[0] = .1 * mm;
  myTeflonEnergy[1] = 5.0 * eV;
  //myTeflonRI[1] = 1.4;
  myTeflonAbs[1] = .1 * mm;
  myTeflonEnergy[2] = 10.0 * eV;
  //myTeflonRI[2] = 1.4;
  myTeflonAbs[2] = .1 * mm;
  // myTeflon->AddProperty("RINDEX", myTeflonEnergy,    myTeflonRI, 3 );
  myTeflon->AddProperty("ABSLENGTH", myTeflonEnergy, myTeflonAbs, 3);
  fTeflon->SetMaterialPropertiesTable(myTeflon);

  //---------------------------------------------------------------------------------
  // Aluminum
  //---------------------------------------------------------------------------------

  G4MaterialPropertiesTable* myAlumi = new G4MaterialPropertiesTable;
  G4double myAlumiEnergy[3], myAlumiAbs[3];
  myAlumiEnergy[0] = 1.0 * eV;
  //myAlumiRI[0] = 1.;
  myAlumiAbs[0] = .1 * mm;
  myAlumiEnergy[1] = 5.0 * eV;
  //myAlumiRI[1] = 1.;
  myAlumiAbs[1] = .1 * mm;
  myAlumiEnergy[2] = 10.0 * eV;
  //myAlumiRI[2] = 1.;
  myAlumiAbs[2] = .1 * mm;
  //  myAlumi->AddProperty("RINDEX", myAlumiEnergy,    myAlumiRI, 3 );
  myAlumi->AddProperty("ABSLENGTH", myAlumiEnergy, myAlumiAbs, 3);
  fAluminum->SetMaterialPropertiesTable(myAlumi);

  //---------------------------------------------------------------------------------
  // Bialkali
  //---------------------------------------------------------------------------------
  // this is used instead of Bialkali for the cathode as fake material to stop
  // photons FusedSilica can not be used, because of the window disk

  G4MaterialPropertiesTable* myBialkali = new G4MaterialPropertiesTable();
  G4double myBialkaliEnergy[4], myBialkaliRI[4], myBialkaliAbs[4];
  myBialkaliEnergy[0] = 1.0 * eV;
  myBialkaliRI[0] = DSParameters::Get()->GetPhotocathodeVisRind();
  myBialkaliAbs[0] = 0.1 * nm;
  myBialkaliEnergy[1] = 8.0 * eV;
  myBialkaliRI[1] = DSParameters::Get()->GetPhotocathodeVisRind();
  ;
  myBialkaliAbs[1] = 0.1 * nm;
  myBialkaliEnergy[2] = 8.3 * eV;
  myBialkaliRI[2] = DSParameters::Get()->GetPhotocathodeUVRind();
  myBialkaliAbs[2] = 0.1 * nm;
  myBialkaliEnergy[3] = 20.0 * eV;
  myBialkaliRI[3] = DSParameters::Get()->GetPhotocathodeUVRind();
  myBialkaliAbs[3] = 0.1 * nm;
  myBialkali->AddProperty("RINDEX", myBialkaliEnergy, myBialkaliRI, 4);
  myBialkali->AddProperty("ABSLENGTH", myBialkaliEnergy, myBialkaliAbs, 4);
  fBialkali->SetMaterialPropertiesTable(myBialkali);

  // FakePhotocathode - - - same properties as bialkali
  fFakePhotocathode->SetMaterialPropertiesTable(myBialkali);

  //---------------------------------------------------------------------------------
  // Metal Silicon
  //---------------------------------------------------------------------------------
  const int mySiliconRindDim = 46;  // from https://refractiveindex.info/?shelf=main&book=Si&page=Aspnes
  double mySiliconRindEne[mySiliconRindDim] = {0.2066, 0.2101, 0.2138, 0.2175, 0.2214, 0.2254, 0.2296, 0.2339, 0.2384, 0.2431, 0.2480, 0.2530, 0.2583, 0.2638, 0.2695, 0.2755, 0.2818, 0.2883, 0.2952, 0.3024, 0.3100, 0.3179, 0.3263,
                                               0.3351, 0.3444, 0.3542, 0.3647, 0.3757, 0.3875, 0.3999, 0.4133, 0.4275, 0.4428, 0.4592, 0.4769, 0.4959, 0.5166, 0.5391, 0.5636, 0.5904, 0.6199, 0.6525, 0.6888, 0.7293, 0.7749, 0.8266};  // um
  double mySiliconRind[mySiliconRindDim] = {1.010, 1.083, 1.133, 1.186, 1.247, 1.340, 1.471, 1.579, 1.589, 1.571, 1.570, 1.597, 1.658, 1.764, 1.988, 2.452, 3.120, 4.087, 4.888, 5.020, 5.010, 5.016, 5.065,
                                            5.156, 5.296, 5.610, 6.522, 6.709, 6.062, 5.570, 5.222, 4.961, 4.753, 4.583, 4.442, 4.320, 4.215, 4.123, 4.042, 3.969, 3.906, 3.847, 3.796, 3.752, 3.714, 3.673};  // rindex
  for (int ij = 0; ij < mySiliconRindDim; ij++) { mySiliconRindEne[ij] = h_Planck * c_light / (1000. * mySiliconRindEne[ij] * nm); }
  // Now reverse the two arrays because energy must be increasing
  G4double mySiliconRindEneReversed[mySiliconRindDim];
  G4double mySiliconRindReversed[mySiliconRindDim];
  for (G4int i = 0; i < mySiliconRindDim; i++) {
    mySiliconRindEneReversed[i] = mySiliconRindEne[mySiliconRindDim - 1 - i];
    mySiliconRindReversed[i] = mySiliconRind[mySiliconRindDim - 1 - i];
  }

  // same goal as bialkali
  G4MaterialPropertiesTable* myMetalSilicon = new G4MaterialPropertiesTable();
  // http://refractiveindex.info/?shelf=main&book=Si&page=Pierce
  // myMetalSilicon->AddProperty("IMAGINARYRINDEX",  mySiEnergy  ,mySiImRind ,
  // 6); myMetalSilicon->AddProperty("RINDEX",    myBialkaliEnergy,
  // myBialkaliRI,  4);
  myMetalSilicon->AddProperty("RINDEX", mySiliconRindEneReversed, mySiliconRindReversed, mySiliconRindDim);
  myMetalSilicon->AddProperty("ABSLENGTH", myBialkaliEnergy, myBialkaliAbs, 4);
  fMetalSilicon->SetMaterialPropertiesTable(myMetalSilicon);
  fMetalCopper->SetMaterialPropertiesTable(myMetalSilicon);

  //---------------------------------------------------------------------------------
  // Liquid Argon
  //---------------------------------------------------------------------------------
  G4MaterialPropertiesTable* myLiquidArgon           = new G4MaterialPropertiesTable();
  G4MaterialPropertiesTable* myLArAboveGrid          = new G4MaterialPropertiesTable();
  G4MaterialPropertiesTable* myVetoLiquidArgon       = new G4MaterialPropertiesTable();
  G4MaterialPropertiesTable* myPseudoArgon           = new G4MaterialPropertiesTable();


  dim = 0;
  G4double LArRefIndex[1000], LArRefIndexEne[1000], LArRayLength[1000], LArAbsLength[1000];
  G4double PArRefIndex[1000];
   for (int ij = 0; ij < 55; ij++) {
    LArRefIndexEne[ij] = (1.0 + ij * 0.2) * eV;
    const G4double lambda = h_Planck * c_light / LArRefIndexEne[ij];
    LArRefIndex[ij]  = GetLArRefIndex(lambda / nm);
    LArRayLength[ij] = GetLArRayLength(lambda / nm) * DSParameters::Get()->GetLArRayleighScale();
    PArRefIndex[ij]  = DSParameters::Get()->GetPArRind();
    if (lambda < 150 * nm) LArAbsLength[ij] = DSParameters::Get()->GetLiquidArgonUVAbs() * m;
    if (lambda > 150 * nm) LArAbsLength[ij] = DSParameters::Get()->GetLiquidArgonVisAbs() * m;
    dim = ij;
  }
  int LArRefIndexDim = dim + 1;

  myLiquidArgon->AddProperty("RINDEX",        LArRefIndexEne, LArRefIndex,  dim + 1);
  myLiquidArgon->AddProperty("RAYLEIGH",      LArRefIndexEne, LArRayLength, dim + 1);
  myLiquidArgon->AddProperty("ABSLENGTH",     LArRefIndexEne, LArAbsLength, dim + 1);

  myLArAboveGrid->AddProperty("RINDEX",       LArRefIndexEne, LArRefIndex,  dim + 1);
  myLArAboveGrid->AddProperty("RAYLEIGH",     LArRefIndexEne, LArRayLength, dim + 1);
  myLArAboveGrid->AddProperty("ABSLENGTH",    LArRefIndexEne, LArAbsLength, dim + 1);
  
  myVetoLiquidArgon->AddProperty("RINDEX",    LArRefIndexEne, LArRefIndex,  dim + 1);
  myVetoLiquidArgon->AddProperty("RAYLEIGH",  LArRefIndexEne, LArRayLength, dim + 1);
  myVetoLiquidArgon->AddProperty("ABSLENGTH", LArRefIndexEne, LArAbsLength, dim + 1);

  myPseudoArgon->AddProperty("RINDEX",        LArRefIndexEne, PArRefIndex,  dim + 1);
  myPseudoArgon->AddProperty("RAYLEIGH",      LArRefIndexEne, LArRayLength, dim + 1);
  myPseudoArgon->AddProperty("ABSLENGTH",     LArRefIndexEne, LArAbsLength, dim + 1);

  myLiquidArgon->AddConstProperty( "ELECTRICFIELD", DSStorage::Get()->GetDriftField(), true);   
  myLArAboveGrid->AddConstProperty("ELECTRICFIELD", DSStorage::Get()->GetExtractionField(), true);  
  myVetoLiquidArgon->AddConstProperty("ELECTRICFIELD", 0 * volt / cm, true); 
  myPseudoArgon->AddConstProperty("ELECTRICFIELD", 0 * volt / cm, true); 

  
  //myLiquidArgon->AddConstProperty("TOTALNUM_INT_SITES", -1, true);

  fLiquidArgon->SetMaterialPropertiesTable(myLiquidArgon);
  fXeLiquidArgon->SetMaterialPropertiesTable(myLiquidArgon);

  fNSLiquidArgon->SetMaterialPropertiesTable(myVetoLiquidArgon);
  fIVLiquidArgon->SetMaterialPropertiesTable(myVetoLiquidArgon); // IV 
  fPScintVetoLiquidArgon->SetMaterialPropertiesTable(myVetoLiquidArgon);
  fOVLiquidArgon->SetMaterialPropertiesTable(myVetoLiquidArgon); // OV

  fPseudoArgon->SetMaterialPropertiesTable(myPseudoArgon);
  // Liquid Argon above grid
  fLArAboveGrid->SetMaterialPropertiesTable(myLArAboveGrid);

  //---------------------------------------------------------------------------------
  // Liquid Argon
  //---------------------------------------------------------------------------------
  // G4MaterialPropertiesTable* myLiquidArgon = new G4MaterialPropertiesTable();
  // dim = 0;
  // G4double LArRefIndex[1000], LArRefIndexEne[1000], LArRayLength[1000], LArAbsLength[1000];
  // for (int ij = 0; ij < 55; ij++) {
  //   LArRefIndexEne[ij] = (1.0 + ij * 0.2) * eV;
  //   const G4double lambda = h_Planck * c_light / LArRefIndexEne[ij];
  //   LArRefIndex[ij] = GetLArRefIndex(lambda / nm);
  //   LArRayLength[ij] = GetLArRayLength(lambda / nm) * DSParameters::Get()->GetLArRayleighScale();
  //   // L700 why no to modify LAr absorption length for VUV light?
  //   if (lambda < 150 * nm) LArAbsLength[ij] = DSParameters::Get()->GetLiquidArgonUVAbs() * m;
  //   if (lambda > 150 * nm) LArAbsLength[ij] = DSParameters::Get()->GetLiquidArgonVisAbs() * m;
  //   dim = ij;
  // }
  // myLiquidArgon->AddProperty("RINDEX", LArRefIndexEne, LArRefIndex, dim + 1);
  // myLiquidArgon->AddProperty("RAYLEIGH", LArRefIndexEne, LArRayLength, dim + 1);
  // myLiquidArgon->AddProperty("ABSLENGTH", LArRefIndexEne, LArAbsLength, dim + 1);

  // // field for NEST
  // myLiquidArgon->AddConstProperty("ELECTRICFIELD", 1000 * volt / cm,
  //                                 true);  // for missed nooks and crannies
  // myLiquidArgon->AddConstProperty("TOTALNUM_INT_SITES", -1, true);

  // fLiquidArgon->SetMaterialPropertiesTable(myLiquidArgon);



  // fXeLiquidArgon->SetMaterialPropertiesTable(myLiquidArgon);





  // //---------------------------------------------------------------------------------
  // // Veto Liquid Argon
  // //---------------------------------------------------------------------------------
  // G4MaterialPropertiesTable* myVetoLiquidArgon = new G4MaterialPropertiesTable();
  // dim = 0;

  // for (int ij = 0; ij < 55; ij++) {
  //   LArRefIndexEne[ij] = (1.0 + ij * 0.2) * eV;
  //   const G4double lambda = h_Planck * c_light / LArRefIndexEne[ij];
  //   LArRefIndex[ij] = GetLArRefIndex(lambda / nm);
  //   LArRayLength[ij] = GetLArRayLength(lambda / nm) * DSParameters::Get()->GetLArRayleighScale();
  //   // L700 why no to modify LAr absorption length for VUV light?
  //   if (lambda < 150 * nm) LArAbsLength[ij] = DSParameters::Get()->GetLiquidArgonUVAbs() * m;
  //   if (lambda > 150 * nm) LArAbsLength[ij] = DSParameters::Get()->GetLiquidArgonVisAbs() * m;
  //   dim = ij;
  // }
  // myVetoLiquidArgon->AddProperty("RINDEX", LArRefIndexEne, LArRefIndex, dim + 1);
  // myVetoLiquidArgon->AddProperty("RAYLEIGH", LArRefIndexEne, LArRayLength, dim + 1);
  // myVetoLiquidArgon->AddProperty("ABSLENGTH", LArRefIndexEne, LArAbsLength, dim + 1);

  // // field for DSLight3
  // myVetoLiquidArgon->AddConstProperty("ELECTRICFIELD", 0 * volt / cm,
  //                                     true);  // for missed nooks and crannies

  // fNSLiquidArgon->SetMaterialPropertiesTable(myVetoLiquidArgon);
  // fIVLiquidArgon->SetMaterialPropertiesTable(myVetoLiquidArgon);
  // fPScintVetoLiquidArgon->SetMaterialPropertiesTable(myVetoLiquidArgon);
  // fOVLiquidArgon->SetMaterialPropertiesTable(myVetoLiquidArgon);

  // //---------------------------------------------------------------------------------
  // // Pseudo Argon
  // //---------------------------------------------------------------------------------
  // G4MaterialPropertiesTable* myPseudoArgon = new G4MaterialPropertiesTable();
  // dim = 0;
  // G4double PArRefIndex[1000], PArRefIndexEne[1000], PArRayLength[1000], PArAbsLength[1000];
  // for (int ij = 0; ij < 55; ij++) {
  //   PArRefIndexEne[ij] = (1.0 + ij * 0.2) * eV;
  //   const G4double lambda = h_Planck * c_light / PArRefIndexEne[ij];
  //   PArRefIndex[ij] = DSParameters::Get()->GetPArRind();
  //   PArRayLength[ij] = GetLArRayLength(lambda / nm) * DSParameters::Get()->GetLArRayleighScale();
  //   // L700 why no to modify LAr absorption length for VUV light?
  //   if (lambda < 150 * nm) PArAbsLength[ij] = DSParameters::Get()->GetLiquidArgonUVAbs() * m;
  //   if (lambda > 150 * nm) PArAbsLength[ij] = DSParameters::Get()->GetLiquidArgonVisAbs() * m;
  //   dim = ij;
  // }
  // myPseudoArgon->AddProperty("RINDEX", PArRefIndexEne, PArRefIndex, dim + 1);
  // myPseudoArgon->AddProperty("RAYLEIGH", PArRefIndexEne, PArRayLength, dim + 1);
  // myPseudoArgon->AddProperty("ABSLENGTH", PArRefIndexEne, PArAbsLength, dim + 1);
  // int LArRefIndexDim = dim + 1;  // store array dimension for TPB

  // // field for NEST
  // myPseudoArgon->AddConstProperty("ELECTRICFIELD", 3000 * volt / cm,
  //                                 true);  // for missed nooks and crannies
  // myPseudoArgon->AddConstProperty("TOTALNUM_INT_SITES", -1, true);

  // fPseudoArgon->SetMaterialPropertiesTable(myPseudoArgon);
  // // Liquid Argon above grid
  // fLArAboveGrid->SetMaterialPropertiesTable(myPseudoArgon);
  
  //---------------------------------------------------------------------------------
  // Gaseous Argon
  //---------------------------------------------------------------------------------

  G4MaterialPropertiesTable* myGaseousArgon = new G4MaterialPropertiesTable();
  G4double GArRefIndex[1000], GArRefIndexEne[1000], GArRayLength[1000], GArAbsLength[1000];
  dim = 0;
  for (int ij = 0; ij < 55; ij++) {
    GArRefIndexEne[ij] = (1.4 + ij * 0.2) * eV;
    const G4double lambda = h_Planck * c_light / GArRefIndexEne[ij];
    //    GArRefIndex[ij]    =
    //    GetGArRefIndex(lambda/nm)*DSParameters::Get()->GetGArRindexScale();
    GArRefIndex[ij] = GetGArRefIndex(lambda / nm);
    GArRayLength[ij] = GetGArRayLength(lambda / nm);
    if (lambda < 150 * nm) GArAbsLength[ij] = DSParameters::Get()->GetGaseousArgonUVAbs() * m;
    if (lambda > 150 * nm) GArAbsLength[ij] = DSParameters::Get()->GetGaseousArgonVisAbs() * m;
    dim = ij;
  }

  myGaseousArgon->AddProperty("RINDEX", GArRefIndexEne, GArRefIndex, dim + 1);
  myGaseousArgon->AddProperty("RAYLEIGH", GArRefIndexEne, GArRayLength, dim + 1);
  myGaseousArgon->AddProperty("ABSLENGTH", GArRefIndexEne, GArAbsLength, dim + 1);

  // field for NEST

  myGaseousArgon->AddConstProperty("ELECTRICFIELD", 3000 * volt / cm,
                                   true);  // for missed nooks and crannies
  myGaseousArgon->AddConstProperty("TOTALNUM_INT_SITES", -1, true);

  fGaseousArgon->SetMaterialPropertiesTable(myGaseousArgon);

  //---------------------------------------------------------------------------------
  // TPB
  //---------------------------------------------------------------------------------

  G4MaterialPropertiesTable* myTPB = new G4MaterialPropertiesTable();

  G4double TPB_ENE[1000], TPB_EMISSION_VAL[1000], TPB_ABSORPTION_VAL[1000];
  dim = 0;
  G4double TPBNorma = 0;
  ifstream ftpb_emission("../data/detector/tpb_emission_spectrum.dat");
  // ifstream
  // ftpb_emission("../data/detector/tpb_emission_spectrum_tempdep.dat");
  if (!ftpb_emission.is_open()) DSLog(fatal) << "ERROR: Could not open TPB emission file" << endlog;
  while (!ftpb_emission.eof()) {
    ftpb_emission >> myene >> myvalue;
    if (ftpb_emission.eof()) break;
    TPB_ENE[dim] = myene * eV;
    if (hc / (TPB_ENE[dim] / nm) > 600) {
      TPB_EMISSION_VAL[dim] = 0;
    } else {
      TPB_EMISSION_VAL[dim] = myvalue;
      TPBNorma += myvalue;
    }
    dim++;
  }
  ftpb_emission.close();

  // PA normalise the emission spectrum here
  for (int ii = 0; ii < dim; ++ii) { TPB_EMISSION_VAL[ii] /= TPBNorma; }

  G4double TPB_ENE2[1000], TPB_WLS_ABSORPTION_VAL[1000];
  dim = 0;
  ifstream ftpb_absorption("../data/detector/tpb_absorption_length.dat");
  if (!ftpb_absorption.is_open()) DSLog(fatal) << "ERROR: Could not open TPB absorption file" << endlog;
  while (!ftpb_absorption.eof()) {
    ftpb_absorption >> myene >> myvalue;
    TPB_ENE2[dim] = myene * eV;
    if (ftpb_absorption.eof()) break;
    TPB_WLS_ABSORPTION_VAL[dim] = myvalue * m * DSParameters::Get()->GetWLSAbsorptionFactor();
    if (myene > 8.3) {
      TPB_ABSORPTION_VAL[dim] = DSParameters::Get()->GetTPBUVAbs() * m;
      //TPB_RINDEX[dim] = DSParameters::Get()->GetTPBUVRind();
    } else {
      TPB_ABSORPTION_VAL[dim] = DSParameters::Get()->GetTPBVisAbs() * m;
      //TPB_RINDEX[dim] = DSParameters::Get()->GetTPBVisRind();
    }
    dim++;
  }
  ftpb_absorption.close();

  myTPB->AddProperty("WLSCOMPONENT", TPB_ENE, TPB_EMISSION_VAL, dim);
  myTPB->AddProperty("WLSABSLENGTH", TPB_ENE2, TPB_WLS_ABSORPTION_VAL, dim);
  myTPB->AddProperty("ABSLENGTH", TPB_ENE2, TPB_ABSORPTION_VAL, dim);
  // myTPB->AddProperty("RINDEX",      TPB_ENE, TPB_RINDEX,             dim);
  myTPB->AddProperty("RINDEX", LArRefIndexEne, LArRefIndex, LArRefIndexDim);  // set TPB RINDEX to same as LAr
  myTPB->AddConstProperty("WLSMEANNUMBERPHOTONS", DSParameters::Get()->GetWLSMeanNumberPhotons());
  myTPB->AddConstProperty("WLSTIMECONSTANT", DSParameters::Get()->GetWLSTimeConstant_ns() * ns);
  myTPB->AddConstProperty("WLSEFFICIENCY", DSParameters::Get()->GetWLSEfficiency(), true);

  G4double myscintYield_TPB[3], myscintTimeconstant_TPB[3];

  myscintYield_TPB[0]= 0.19;
  myscintTimeconstant_TPB[0]= 5.2e-3 * us;

  //myscintYield_TPB[1]= 0.41 * (9.6/5.87) / 0.195 * 0;
  //myscintTimeconstant_TPB[1]= 0.9 * us ;

  myscintYield_TPB[1]=  0.81;
  myscintTimeconstant_TPB[1]= 2.3e3 * us;

  myscintYield_TPB[2]= 0.41 * (9.6/5.87) / 0.195; // in nph / keV
  myTPB->AddConstProperty("MATSCINTILLATIONYIELD1", myscintYield_TPB[0], true);
  myTPB->AddConstProperty("MATSCINTILLATIONTIMECONSTANT1", myscintTimeconstant_TPB[0], true);
  myTPB->AddProperty("MATSCINTILLATIONCOMPONENT1", TPB_ENE, TPB_EMISSION_VAL, dim, true); 
  
  myTPB->AddConstProperty("MATSCINTILLATIONYIELD2", myscintYield_TPB[1], true);
  myTPB->AddConstProperty("MATSCINTILLATIONTIMECONSTANT2", myscintTimeconstant_TPB[1], true);
  myTPB->AddProperty("MATSCINTILLATIONCOMPONENT2", TPB_ENE, TPB_EMISSION_VAL, dim, true); 


  myTPB->AddConstProperty("MATSCINTILLATIONYIELD", myscintYield_TPB[2]/ keV, true); // in nph / keV
  myTPB->AddConstProperty("MATRESOLUTIONSCALE", 1, true); // in nph / keV

  fTPB->SetMaterialPropertiesTable(myTPB);



  //---------------------------------------------------------------------------------
  // ESR
  //---------------------------------------------------------------------------------



  G4MaterialPropertiesTable* myESR = new G4MaterialPropertiesTable();

  G4double ESR_ENE[1000], ESR_EMISSION_VAL[1000], ESR_ABSORPTION_VAL[1000];
  dim1 = 0;
  G4double ESRNorma = 0;
  ifstream fESR_emission("../data/detector/pen_emission_spectrum_tempdep.dat");
  if (!fESR_emission.is_open()) DSLog(fatal) << "ERROR: Could not open PEN emission file" << endlog;
  while (!fESR_emission.eof()) {
    fESR_emission >> myene >> myvalue;
    if (fESR_emission.eof()) break;
    ESR_ENE[dim1] = myene * eV;
    if (hc / (ESR_ENE[dim1] / nm) > 600) {
      ESR_EMISSION_VAL[dim1] = 0;
    } else {
      ESR_EMISSION_VAL[dim1] = myvalue;
      ESRNorma += myvalue;
    }
    dim1++;
  }
  fESR_emission.close();

  for (int ii = 0; ii < dim1; ++ii) ESR_EMISSION_VAL[ii] /= ESRNorma;

  // LL commented ESR_RINDEX
  G4double ESR_ENE2[1000], ESR_WLS_ABSORPTION_VAL[1000], /*ESR_RINDEX[1000],*/ ESR_RAYL[1000];
  dim2 = 0;
  ifstream fESR_absorption("../data/detector/pen_absorption_length.dat");
  if (!fESR_absorption.is_open()) DSLog(fatal) << "ERROR: Could not open PEN absorption file" << endlog;
  while (!fESR_absorption.eof()) {
    fESR_absorption >> myene2 >> myvalue;
    if (fESR_absorption.eof()) break;
    ESR_ENE2[dim2] = myene2 * eV;
    ESR_WLS_ABSORPTION_VAL[dim2] = myvalue * DSParameters::Get()->GetWLSAbsorptionFactor() * m / 4.;
    // cout << "PEN WLS Absorption Factor: " << ESR_WLS_ABSORPTION_VAL[dim2] <<
    // endl;
    if (myene2 > 8.3) {
      ESR_ABSORPTION_VAL[dim2] = DSParameters::Get()->GetPENUVAbs() * m;
    // LL commented the following line
      //ESR_RINDEX[dim2] = DSParameters::Get()->GetPENUVRind();
      ESR_RAYL[dim2] = DSParameters::Get()->GetPENUVRaylLength() * m;
    } else {
      ESR_ABSORPTION_VAL[dim2] = DSParameters::Get()->GetPENVisAbs() * m;
    // LL commented the following line
      //ESR_RINDEX[dim2] = DSParameters::Get()->GetPENVisRind();
      ESR_RAYL[dim2] = DSParameters::Get()->GetPENVisRaylLength() * m;
    }
    dim2++;
  }
  fESR_absorption.close();

  myESR->AddProperty("WLSCOMPONENT", ESR_ENE, ESR_EMISSION_VAL, dim1);
  myESR->AddProperty("WLSABSLENGTH", ESR_ENE2, ESR_WLS_ABSORPTION_VAL, dim2);
  myESR->AddProperty("ABSLENGTH", ESR_ENE2, ESR_ABSORPTION_VAL, dim2);
  // LL commented the following line
  //myESR->AddProperty("RINDEX", ESR_ENE2, ESR_RINDEX, dim2);
  myESR->AddProperty("RAYLEIGH", ESR_ENE2, ESR_RAYL, dim2);
  //  myESR->AddProperty("RAYLEIGH",    ESR_ENE2,
  //  DSParameters::Get()->GetPENRaylLength()*um);
  // LL uncommented the following line
  myESR->AddProperty("RINDEX", LArRefIndexEne, LArRefIndex, LArRefIndexDim); //set PEN RINDEX to same as LAr
  myESR->AddConstProperty("WLSMEANNUMBERPHOTONS", DSParameters::Get()->GetWLSMeanNumberPhotons());
  myESR->AddConstProperty("WLSTIMECONSTANT", DSParameters::Get()->GetPENTimeConstant_ns() * ns);
  // LL "DSParameters::Get()->GetPENEfficiency() --> DSParameters::Get()->GetPENEfficiency()*0.22"
  myESR->AddConstProperty("WLSEFFICIENCY", DSParameters::Get()->GetPENEfficiency()*0.22, true);
 

  //G4double myscintYield[4]= { 9.26e-1, 1.96e-3, 6.40e-2, 0.200*1.77};  
  //G4double myscintTimeconstant[3]= {16.5 * ns, 2110 * ns,  111 * ns};  
  // parameters provided by Marco 
  //double probs[4] = {0.52912342, 0.23811609, 0.13964312, 0.093117377};
  //double taus[4] = {16.79,      110.3,       2113.,      47980.}; 
  // skip slowest compoment
  G4double myscintYield[4]= { 0.52912342, 0.23811609, 0.13964312, 0.200*1.77};  
  G4double myscintTimeconstant[3]= {16.79 * ns, 110.3 * ns,  2113 * ns};  

  myESR->AddConstProperty("MATSCINTILLATIONYIELD1", myscintYield[0],true);
  myESR->AddConstProperty("MATSCINTILLATIONTIMECONSTANT1", myscintTimeconstant[0],true);
  myESR->AddProperty("MATSCINTILLATIONCOMPONENT1", ESR_ENE, ESR_EMISSION_VAL, dim1, true); 
  
  myESR->AddConstProperty("MATSCINTILLATIONYIELD2", myscintYield[1], true);
  myESR->AddConstProperty("MATSCINTILLATIONTIMECONSTANT2", myscintTimeconstant[1], true);
  myESR->AddProperty("MATSCINTILLATIONCOMPONENT2", ESR_ENE, ESR_EMISSION_VAL, dim1, true); 

  myESR->AddConstProperty("MATSCINTILLATIONYIELD3", myscintYield[2], true);
  myESR->AddConstProperty("MATSCINTILLATIONTIMECONSTANT3", myscintTimeconstant[2], true);
  myESR->AddProperty("MATSCINTILLATIONCOMPONENT3", ESR_ENE, ESR_EMISSION_VAL, dim1, true); 

  myESR->AddConstProperty("MATSCINTILLATIONYIELD", myscintYield[3]/ keV, true); // in nph / keV
  myESR->AddConstProperty("MATRESOLUTIONSCALE", 1, true); // in nph / keV
  

  fESR->SetMaterialPropertiesTable(myESR);


  //---------------------------------------------------------------------------------
  // PEN
  //---------------------------------------------------------------------------------

  G4MaterialPropertiesTable* myPEN = new G4MaterialPropertiesTable();

  G4double PEN_ENE[1000], PEN_EMISSION_VAL[1000], PEN_ABSORPTION_VAL[1000];
  dim1 = 0;
  G4double PENNorma = 0;
  ifstream fpen_emission("../data/detector/pen_emission_spectrum_tempdep.dat");
  if (!fpen_emission.is_open()) DSLog(fatal) << "ERROR: Could not open PEN emission file" << endlog;
  while (!fpen_emission.eof()) {
    fpen_emission >> myene >> myvalue;
    if (fpen_emission.eof()) break;
    PEN_ENE[dim1] = myene * eV;
    if (hc / (PEN_ENE[dim1] / nm) > 600) {
      PEN_EMISSION_VAL[dim1] = 0;
    } else {
      PEN_EMISSION_VAL[dim1] = myvalue;
      PENNorma += myvalue;
    }
    dim1++;
  }
  fpen_emission.close();

  for (int ii = 0; ii < dim1; ++ii) PEN_EMISSION_VAL[ii] /= PENNorma;

  G4double PEN_ENE2[1000], PEN_WLS_ABSORPTION_VAL[1000], PEN_RINDEX[1000], PEN_RAYL[1000];
  dim2 = 0;
  ifstream fpen_absorption("../data/detector/pen_absorption_length.dat");
  if (!fpen_absorption.is_open()) DSLog(fatal) << "ERROR: Could not open PEN absorption file" << endlog;
  while (!fpen_absorption.eof()) {
    fpen_absorption >> myene2 >> myvalue;
    if (fpen_absorption.eof()) break;
    PEN_ENE2[dim2] = myene2 * eV;
    PEN_WLS_ABSORPTION_VAL[dim2] = myvalue * DSParameters::Get()->GetWLSAbsorptionFactor() * m / 4.;
    // cout << "PEN WLS Absorption Factor: " << PEN_WLS_ABSORPTION_VAL[dim2] <<
    // endl;
    if (myene2 > 8.3) {
      PEN_ABSORPTION_VAL[dim2] = DSParameters::Get()->GetPENUVAbs() * m;
      PEN_RINDEX[dim2] = DSParameters::Get()->GetPENUVRind();
      PEN_RAYL[dim2] = DSParameters::Get()->GetPENUVRaylLength() * m;
    } else {
      PEN_ABSORPTION_VAL[dim2] = DSParameters::Get()->GetPENVisAbs() * m;
      PEN_RINDEX[dim2] = DSParameters::Get()->GetPENVisRind();
      PEN_RAYL[dim2] = DSParameters::Get()->GetPENVisRaylLength() * m;
    }
    dim2++;
  }
  fpen_absorption.close();

  myPEN->AddProperty("WLSCOMPONENT", PEN_ENE, PEN_EMISSION_VAL, dim1);
  myPEN->AddProperty("WLSABSLENGTH", PEN_ENE2, PEN_WLS_ABSORPTION_VAL, dim2);
  myPEN->AddProperty("ABSLENGTH", PEN_ENE2, PEN_ABSORPTION_VAL, dim2);
  myPEN->AddProperty("RINDEX", PEN_ENE2, PEN_RINDEX, dim2);
  myPEN->AddProperty("RAYLEIGH", PEN_ENE2, PEN_RAYL, dim2);
  //  myPEN->AddProperty("RAYLEIGH",    PEN_ENE2,
  //  DSParameters::Get()->GetPENRaylLength()*um); myPEN->AddProperty("RINDEX",
  //  LArRefIndexEne, LArRefIndex, LArRefIndexDim); //set PEN RINDEX to same as
  //  LAr
  myPEN->AddConstProperty("WLSMEANNUMBERPHOTONS", DSParameters::Get()->GetWLSMeanNumberPhotons());
  myPEN->AddConstProperty("WLSTIMECONSTANT", DSParameters::Get()->GetPENTimeConstant_ns() * ns);
  myPEN->AddConstProperty("WLSEFFICIENCY", DSParameters::Get()->GetPENEfficiency(), true);
  fPEN->SetMaterialPropertiesTable(myPEN);

  //---------------------------------------------------------------------------------
  // GridSteel
  //---------------------------------------------------------------------------------
  G4MaterialPropertiesTable* myGridSteel = new G4MaterialPropertiesTable();
  // G4double myGridSteelEnergy[4],  myGridSteelAbs[4], myGridRefIndex[4];
  dim = 0;
  G4double GridRefIndex[1000], GridRefIndexEne[1000], GridAbsLength[1000];
  for (int ij = 0; ij < 55; ij++) {
    GridRefIndexEne[ij] = (1.0 + ij * 0.2) * eV;
    const G4double lambda = h_Planck * c_light / GridRefIndexEne[ij];
    GridRefIndex[ij] = GetLArRefIndex(lambda / nm) * DSParameters::Get()->GetGridSteelRindScale();
    if (lambda < 150 * nm) GridAbsLength[ij] = DSParameters::Get()->GetGridSteelUVAbs() * m;
    else
      GridAbsLength[ij] = DSParameters::Get()->GetGridSteelVisAbs() * m;
    dim = ij;
  }

  myGridSteel->AddProperty("RINDEX", GridRefIndexEne, GridRefIndex, dim + 1);
  myGridSteel->AddProperty("ABSLENGTH", GridRefIndexEne, GridAbsLength, dim + 1);
  fGridSteel->SetMaterialPropertiesTable(myGridSteel);

  //---------------------------------------------------------------------------------
  // Kovar
  //---------------------------------------------------------------------------------
  G4MaterialPropertiesTable* myKovar = new G4MaterialPropertiesTable();
  G4double myKovarEnergy[3], myKovarAbs[3];
  myKovarEnergy[0] = 0.1 * eV;
  myKovarAbs[0] = 0.01 * nm;
  myKovarEnergy[1] = 5.0 * eV;
  myKovarAbs[1] = 0.01 * nm;
  myKovarEnergy[2] = 20.0 * eV;
  myKovarAbs[2] = 0.01 * nm;
  myKovar->AddProperty("ABSLENGTH", myKovarEnergy, myKovarAbs, 3);
  fKovar->SetMaterialPropertiesTable(myKovar);

  //---------------------------------------------------------------------------------
  // Black Hole
  //---------------------------------------------------------------------------------
  G4MaterialPropertiesTable* myBlackHole = new G4MaterialPropertiesTable();
  G4double myBHEnergy[3], myBHRI[3], myBHAbs[3];
  myBHEnergy[0] = 0.1 * eV;
  myBHRI[0] = 1.0;
  myBHAbs[0] = 0.001 * mm;
  myBHEnergy[1] = 5.0 * eV;
  myBHRI[1] = 1.0;
  myBHAbs[1] = 0.001 * mm;
  myBHEnergy[2] = 20.0 * eV;
  myBHRI[2] = 1.0;
  myBHAbs[2] = 0.001 * mm;
  myBlackHole->AddProperty("RINDEX", myBHEnergy, myBHRI, 3);
  myBlackHole->AddProperty("ABSLENGTH", myBHEnergy, myBHAbs, 3);
  fBlackHole->SetMaterialPropertiesTable(myBlackHole);
  fArlon->SetMaterialPropertiesTable(myBlackHole);

  //---------------------------------------------------------------------------------
  // Acrylic
  //---------------------------------------------------------------------------------
  G4MaterialPropertiesTable* myAcrylic = new G4MaterialPropertiesTable();
  // G4double myAcrEnergy[18], myAcrRI[18];

  /*
  myAcrEnergy[0] = 0.1*eV  ; myAcrRI[0] = 1.5;  myAcrAbs[0] = 10*m;
  myAcrEnergy[1] = 5.0*eV ;  myAcrRI[1] = 1.5;  myAcrAbs[1] = 10*m;
  myAcrEnergy[2] = 5.1*eV ;  myAcrRI[2] = 1.9;  myAcrAbs[2] = 1*m;
  myAcrEnergy[3] = 20.0*eV ; myAcrRI[3] = 1.9;  myAcrAbs[3] = 1*nm;
  */

  // Jul'17 -- added acrylic abs length in Genoa - based on
  //  G4double myAcrLambda[18] = {  200.0, 220.0, 240.0, 260.0, 280.0, 300.0,
  //  310.0, 320.0, 330.0, 340.0, 350.0, 360.0, 380.0, 400.0, 420.0, 440.0,
  //  450.0, 800.0 } ; G4double myAcrAbs[18]    = { 0.249165, 0.249165,
  //  0.249165, 0.249165, 2.39349, 22.9885, 69.9301, 100.0, 140.056, 200.0,
  //  270.27, 380.228, 751.88, 1129.94, 1227.6, 1343.72, 1410.44, 1410.44} ;

  const int myAcrSize = 504;
  G4double myAcrAbsEnergy[myAcrSize] = {60,  200, 280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297, 298, 299, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333,
                                        334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 375, 376, 377, 378, 379, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389,
                                        390, 391, 392, 393, 394, 395, 396, 397, 398, 399, 400, 401, 402, 403, 404, 405, 406, 407, 408, 409, 410, 411, 412, 413, 414, 415, 416, 417, 418, 419, 420, 421, 422, 423, 424, 425, 426, 427, 428, 429, 430, 431, 432, 433, 434, 435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445,
                                        446, 447, 448, 449, 450, 451, 452, 453, 454, 455, 456, 457, 458, 459, 460, 461, 462, 463, 464, 465, 466, 467, 468, 469, 470, 471, 472, 473, 474, 475, 476, 477, 478, 479, 480, 481, 482, 483, 484, 485, 486, 487, 488, 489, 490, 491, 492, 493, 494, 495, 496, 497, 498, 499, 500, 501,
                                        502, 503, 504, 505, 506, 507, 508, 509, 510, 511, 512, 513, 514, 515, 516, 517, 518, 519, 520, 521, 522, 523, 524, 525, 526, 527, 528, 529, 530, 531, 532, 533, 534, 535, 536, 537, 538, 539, 540, 541, 542, 543, 544, 545, 546, 547, 548, 549, 550, 551, 552, 553, 554, 555, 556, 557,
                                        558, 559, 560, 561, 562, 563, 564, 565, 566, 567, 568, 569, 570, 571, 572, 573, 574, 575, 576, 577, 578, 579, 580, 581, 582, 583, 584, 585, 586, 587, 588, 589, 590, 591, 592, 593, 594, 595, 596, 597, 598, 599, 600, 601, 602, 603, 604, 605, 606, 607, 608, 609, 610, 611, 612, 613,
                                        614, 615, 616, 617, 618, 619, 620, 621, 622, 623, 624, 625, 626, 627, 628, 629, 630, 631, 632, 633, 634, 635, 636, 637, 638, 639, 640, 641, 642, 643, 644, 645, 646, 647, 648, 649, 650, 651, 652, 653, 654, 655, 656, 657, 658, 659, 660, 661, 662, 663, 664, 665, 666, 667, 668, 669,
                                        670, 671, 672, 673, 674, 675, 676, 677, 678, 679, 680, 681, 682, 683, 684, 685, 686, 687, 688, 689, 690, 691, 692, 693, 694, 695, 696, 697, 698, 699, 700, 701, 702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 718, 719, 720, 721, 722, 723, 724, 725,
                                        726, 727, 728, 729, 730, 731, 732, 733, 734, 735, 736, 737, 738, 739, 740, 741, 742, 743, 744, 745, 746, 747, 748, 749, 750, 751, 752, 753, 754, 755, 756, 757, 758, 759, 760, 761, 762, 763, 764, 765, 766, 767, 768, 769, 770, 771, 772, 773, 774, 775, 776, 777, 778, 779, 780, 800};
  G4double myAcrAbs[myAcrSize] = {
      0,       0,       1.32284, 1.29312, 1.30753, 1.30943, 1.26653, 1.32372, 1.29873, 1.3092,  1.30371, 1.30577, 1.32406, 1.31243, 1.3171,  1.28569, 1.30494, 1.33508, 1.30388, 1.29962, 1.29142, 1.30584, 1.28672, 1.31488, 1.30954, 1.31865, 1.29685, 1.28557, 1.3144,  1.32283, 1.31425, 1.27107, 1.302,   1.29515,
      1.29604, 1.31832, 1.27713, 1.33705, 1.31362, 1.34722, 1.34076, 1.33069, 1.35506, 1.35238, 1.33769, 1.33566, 1.35435, 1.32866, 1.35414, 1.34693, 1.40201, 1.3828,  1.39343, 1.31766, 1.30547, 1.34563, 1.31455, 1.34873, 1.32728, 1.281,   1.35277, 1.31982, 1.30124, 1.28945, 1.29276, 1.2827,  1.28362, 1.29432,
      1.25258, 1.21353, 1.27814, 1.26223, 1.22418, 1.23389, 1.19746, 1.20099, 1.23181, 1.2066,  1.21319, 1.21229, 1.22824, 1.2102,  1.20421, 1.17812, 1.19581, 1.23024, 1.15775, 1.21839, 1.20776, 1.21069, 1.1799,  1.20016, 1.1395,  1.19189, 1.16832, 1.13224, 1.24984, 1.50017, 1.77389, 1.78046, 2.09596, 1.93195,
      2.85164, 3.44997, 2.62686, 4.97381, 5.97033, 14.0861, 18.7694, 22.3729, 23.5664, 24.2125, 24.918,  27.6302, 32.0202, 38.1495, 45.9201, 55.8979, 68.4572, 83.8411, 102.648, 125.635, 153.884, 187.777, 229.237, 278.717, 339.222, 410.576, 495.219, 597.237, 707.002, 839.051, 980.058, 1143.41, 1321,    1510.09,
      1705.82, 1932.72, 2165.05, 2359.05, 2580.08, 2771.07, 3033.05, 3167.62, 3409.31, 3479.23, 3707.19, 3851.96, 3897.3,  4225.34, 4205.58, 4315.72, 4419.21, 4657.92, 4658.72, 4759.31, 4806.68, 4802.77, 5000.05, 5121.6,  5309.11, 5333.83, 5338.58, 5643.52, 5467.12, 5694.31, 5779.82, 5594.7,  5879.86, 5981.84,
      5787.35, 6262.51, 6254.24, 6416.64, 6962.87, 6881.92, 7197.05, 7142.11, 7293.65, 7876.71, 8372.81, 7681.55, 7740.36, 7990.87, 8104.59, 8363.88, 8592.66, 8136.1,  8629.33, 8565.14, 8687.34, 9020.28, 9114.16, 8924.67, 9189.4,  9277.85, 9158.99, 9574.41, 10854.5, 9460.01, 9664.68, 10045.7, 10371.2, 9946.08,
      10623.3, 11275.5, 10485.8, 11033.6, 11213.6, 11160.2, 11517.3, 11466.1, 11184,   11771.9, 11520.6, 12506,   11235.2, 12208.5, 12694.5, 11598.4, 13119.1, 11726.9, 12782.3, 12282.6, 13379.2, 14005.3, 12249.8, 12855.6, 12727.9, 12390.6, 13738.5, 14171.4, 13698.8, 14075.3, 14030.4, 15346.9, 13350.7, 14807,
      14633.6, 17460.4, 15319,   15551.3, 15958,   14883.8, 14163.4, 15687.7, 14723.8, 15681.3, 16460.4, 15074.8, 15287.5, 15808.1, 16863.4, 15432.6, 16162.6, 16028.7, 17023.9, 16573.2, 14573.7, 15343.5, 16263,   14356.3, 14277.8, 15545.9, 16271.4, 14849.5, 14983.5, 15728.8, 15986.8, 14954.6, 16380.3, 14752.6,
      14499.7, 16345.6, 18455.6, 13926.8, 17307.8, 15290.8, 15605.9, 16007.6, 17090.7, 18384.3, 16610.8, 16896,   17407.5, 19442.5, 19325.7, 20769.7, 18066.3, 17648.5, 18863.4, 22099.9, 20030.3, 23828.9, 19256,   20513.2, 26932.9, 20529.7, 18444.9, 20507.2, 20803.8, 23749.3, 20872.5, 21161.2, 20857.3, 21194.2,
      22768.2, 20472.9, 18433.5, 21932.4, 19893.6, 21459.4, 19480.1, 19474.4, 22644.7, 25062.3, 20185.8, 21164,   20515.7, 21679.3, 23751.5, 21387,   24228.8, 22827.1, 22586.3, 21646.6, 19780.9, 18207.4, 17148.8, 17317.4, 14297.8, 15752.4, 15363.6, 14105.6, 12582.5, 12050,   11068.4, 10582.3, 10052.6, 9658.73,
      9156.67, 9611.99, 8740.63, 9141.28, 8545.84, 8675.79, 8096.21, 8927.55, 8966.51, 9515.41, 8602.26, 9440.08, 10283.8, 10265.4, 10888.3, 10528.8, 12258.1, 12131.1, 14388.4, 15769.1, 16420.3, 17624.4, 15974.7, 19390.5, 21176.3, 19823.3, 24277.3, 20107.8, 21017.3, 26572.9, 24786.3, 23206.9, 25619.9, 29196.8,
      22952.4, 27438.1, 22454,   30586.4, 20477.8, 22649.7, 22041.4, 20035.7, 17730.3, 17692.6, 20030.4, 19537,   16425.3, 16144.3, 19343.1, 15152.6, 14117.4, 14447.6, 15711.6, 13035.8, 14160.7, 15208.4, 14528.4, 13168.9, 13102.1, 20340.7, 14727,   14351,   13941.1, 13676.7, 14920.8, 14700,   13243,   12718.4,
      15343.7, 13408.6, 13693.9, 12562.9, 13435.8, 13291.9, 12706.8, 13482.3, 12490.9, 12314.1, 11564,   10328.3, 10482,   9638.82, 9286.56, 8971.19, 8252.47, 7627.67, 7387.04, 6738.8,  6368.22, 6173.61, 5617.03, 5641.01, 5206.5,  4730.11, 4560.32, 4190.11, 3948.63, 3640.44, 3263.39, 3022.28, 2781.11, 2542.86,
      2304.26, 2125.74, 1978.46, 1820.81, 1705.02, 1617.41, 1552.42, 1469.53, 1420.46, 1389.11, 1373.85, 1360.76, 1356.43, 1357.77, 1376.89, 1393.14, 1407.4,  1442.59, 1487.05, 1538.95, 1622.01, 1698.69, 1792.43, 1935.75, 2089.58, 2274.06, 2455.86, 2639.89, 2934.61, 3027.92, 3422.37, 3714.02, 4134.33, 4542.13,
      5160.33, 5363.87, 5977.23, 6442.42, 7013.28, 7784.87, 8675.49, 8942.09, 9518.33, 9760.12, 10260.8, 9922.78, 10979.6, 10000.5, 10391.3, 9134.34, 9521.78, 9595.54, 9179.65, 8279.66, 8223.8,  8020.19, 6951.69, 6875.43, 6592.61, 6159.45, 6023.03, 6023.03};
  for (int i = 0; i < myAcrSize; i++) {
    // cout << myAcrAbsEnergy[i] << " " << myAcrAbs[i] << endl;
    myAcrAbsEnergy[i] = h_Planck * c_light / (nm * myAcrAbsEnergy[i]);
    myAcrAbs[i] *= mm;
  }

  G4double RINDEX_value1[62] = {60,  200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 310, 320, 330, 340, 350, 360, 370, 380, 390, 400, 410, 420, 430, 440, 450, 460, 470, 480, 490,
                                500, 510, 520, 530, 540, 550, 560, 570, 580, 590, 600, 610, 620, 630, 640, 650, 660, 670, 680, 690, 700, 710, 720, 730, 740, 750, 760, 770, 780, 790, 800};
  G4double RINDEX_value2[62] = {1.597, 1.597, 1.584, 1.573, 1.564, 1.556, 1.550, 1.544, 1.539, 1.534, 1.531, 1.527, 1.524, 1.521, 1.519, 1.516, 1.514, 1.512, 1.510, 1.509, 1.507, 1.506, 1.505, 1.503, 1.502, 1.501, 1.500, 1.499, 1.499, 1.498, 1.497,
                                1.496, 1.496, 1.495, 1.494, 1.494, 1.493, 1.493, 1.492, 1.492, 1.491, 1.491, 1.490, 1.490, 1.490, 1.489, 1.489, 1.488, 1.488, 1.488, 1.488, 1.487, 1.487, 1.487, 1.486, 1.486, 1.486, 1.486, 1.485, 1.485, 1.485, 1.485};
  /*
  for (int i=17;i>=0;--i) {
    myAcrLambda[i] *= nm ;
    myAcrEnergy[i]   = h_Planck*c_light / myAcrLambda[i]  ;
    myAcrRI[i]     = 1.5 ;
    myAcrAbs[i]    *= mm ;
  }
  */

  // Now reverse the arrays because energy must be increasing
  G4double myacr_RINDEX_ene[62];
  G4double myacr_RINDEX[62];
  for (G4int i = 0; i < 62; i++) {
    myacr_RINDEX_ene[i] = h_Planck * c_light / (RINDEX_value1[61 - i] * nm);
    myacr_RINDEX[i] = RINDEX_value2[61 - i];
  }
  G4double myAcrAbsEnergyReversed[myAcrSize];
  G4double myAcrAbsReversed[myAcrSize];
  for (G4int i = 0; i < myAcrSize; i++) {
    myAcrAbsEnergyReversed[i] = myAcrAbsEnergy[myAcrSize - 1 - i];
    myAcrAbsReversed[i] = myAcrAbs[myAcrSize - 1 - i];
  }

  myAcrylic->AddProperty("RINDEX", myacr_RINDEX_ene, myacr_RINDEX, 62);
  // myAcrylic->AddProperty("ABSLENGTH", myAcrEnergy, myAcrAbs , 18);
  myAcrylic->AddProperty("ABSLENGTH", myAcrAbsEnergyReversed, myAcrAbsReversed, myAcrSize);
  fAcrylic->SetMaterialPropertiesTable(myAcrylic);

  //---------------------------------------------------------------------------------
  // Set the properties of final scintillator for DS20k veto.
  // (use the acrylic properties for the moment)
  //---------------------------------------------------------------------------------
  fDS20kPlasticScintillator->SetMaterialPropertiesTable(myAcrylic);

  //---------------------------------------------------------------------------------
  // Acrylic--DART
  //---------------------------------------------------------------------------------

  G4MaterialPropertiesTable* myAcrylicDART = new G4MaterialPropertiesTable();
  G4double myAcrDARTEnergy[4], myAcrDARTRI[4], myAcrDARTAbs[4];

  myAcrDARTEnergy[0] = 0.1 * eV;
  myAcrDARTRI[0] = 1.5;
  myAcrDARTAbs[0] = 1. * m;
  myAcrDARTEnergy[1] = 5.0 * eV;
  myAcrDARTRI[1] = 1.5;
  myAcrDARTAbs[1] = 1. * m;
  myAcrDARTEnergy[2] = 5.1 * eV;
  myAcrDARTRI[2] = 1.9;
  myAcrDARTAbs[2] = 1 * m;
  myAcrDARTEnergy[3] = 20.0 * eV;
  myAcrDARTRI[3] = 1.9;
  myAcrDARTAbs[3] = 1 * nm;

  myAcrylicDART->AddProperty("RINDEX", myAcrDARTEnergy, myAcrDARTRI, 4);
  myAcrylicDART->AddProperty("ABSLENGTH", myAcrDARTEnergy, myAcrDARTAbs, 4);
  fAcrylicDART->SetMaterialPropertiesTable(myAcrylicDART);

  //---------------------------------------------------------------------------------
  // Altustipe and BC454
  //---------------------------------------------------------------------------------
  G4MaterialPropertiesTable* myAltustipe = new G4MaterialPropertiesTable();
  G4double myAltustipeEnergy[3], myAltustipeRI[3], myAltustipeAbs[3];
  myAltustipeEnergy[0] = 0.1 * eV;
  myAltustipeRI[0] = 1.6;
  myAltustipeAbs[0] = 120 * cm;
  myAltustipeEnergy[1] = 5.0 * eV;
  myAltustipeRI[1] = 1.6;
  myAltustipeAbs[1] = 1200 * mm;
  myAltustipeEnergy[2] = 20.0 * eV;
  myAltustipeRI[2] = 1.6;
  myAltustipeAbs[2] = 1200 * mm;
  myAltustipe->AddProperty("RINDEX", myAltustipeEnergy, myAltustipeRI, 3);
  myAltustipe->AddProperty("ABSLENGTH", myAltustipeEnergy, myAltustipeAbs, 3);
  fAltustipe->SetMaterialPropertiesTable(myAltustipe);
  fBC454->SetMaterialPropertiesTable(myAltustipe);
}
//---------------------------------------------------------------------------------------
void DSMaterial::DefineBoundaries() {

  //------------------------------------------
  // Lumirror
  //------------------------------------------

  std::vector<G4double> myLumirrorEnergy(DSParameters::Get()->GetLumirrorEnergy());
  std::vector<G4double> myLumirrorRefl(DSParameters::Get()->GetLumirrorReflec());
  const size_t mysize(myLumirrorEnergy.size());
  std::vector<G4double> myLumirrorEff(mysize, DSParameters::Get()->GetLumirrorEff());
  std::vector<G4double> myLumirrorSpecLobe(mysize, DSParameters::Get()->GetLumirrorSpecLobe());
  std::vector<G4double> myLumirrorSpecSpike(mysize, DSParameters::Get()->GetLumirrorSpecSpike());
  std::vector<G4double> myLumirrorBackscat(mysize, DSParameters::Get()->GetLumirrorBackscat());
  std::vector<G4double> myLumirrorRindex(mysize, DSParameters::Get()->GetLumirrorRindex());

  // Load properties into material properties table
  fLumirrorMPT = new G4MaterialPropertiesTable();
  fLumirrorMPT->AddProperty("SPECULARLOBECONSTANT", &myLumirrorEnergy[0], &myLumirrorSpecLobe[0], mysize);
  fLumirrorMPT->AddProperty("SPECULARSPIKECONSTANT", &myLumirrorEnergy[0], &myLumirrorSpecSpike[0], mysize);
  fLumirrorMPT->AddProperty("BACKSCATTERCONSTANT", &myLumirrorEnergy[0], &myLumirrorBackscat[0], mysize);
  fLumirrorMPT->AddProperty("REFLECTIVITY", &myLumirrorEnergy[0], &myLumirrorRefl[0], mysize);
  fLumirrorMPT->AddProperty("EFFICIENCY", &myLumirrorEnergy[0], &myLumirrorEff[0], mysize);
  fLumirrorMPT->AddProperty("RINDEX", &myLumirrorEnergy[0], &myLumirrorRindex[0], mysize);

  //--------------------------------------------------------
  // ElectroPolished StainlessSteel (TPC cryostat from LS)
  //--------------------------------------------------------

  // Using the same energy vector of Lumirror to fill constant properties
  std::vector<G4double> myEpssEnergy(myLumirrorEnergy);
  std::vector<G4double> myEpssRefl(mysize, DSParameters::Get()->GetEpssRefl());
  std::vector<G4double> myEpssEff(mysize, DSParameters::Get()->GetEpssEff());
  std::vector<G4double> myEpssSpecLobe(mysize, DSParameters::Get()->GetEpssSpecLobe());
  std::vector<G4double> myEpssSpecSpike(mysize, DSParameters::Get()->GetEpssSpecSpike());
  std::vector<G4double> myEpssBackscat(mysize, DSParameters::Get()->GetEpssBackscat());
  std::vector<G4double> myEpssRindex(mysize, DSParameters::Get()->GetEpssRindex());

  std::vector<G4double> myPMTBackEnergy(myLumirrorEnergy);
  std::vector<G4double> myPMTBackRefl(mysize, DSParameters::Get()->GetPMTBackRefl());
  std::vector<G4double> myPMTBackEff(mysize, DSParameters::Get()->GetPMTBackEff());
  std::vector<G4double> myPMTBackSpecLobe(mysize, DSParameters::Get()->GetPMTBackSpecLobe());
  std::vector<G4double> myPMTBackSpecSpike(mysize, DSParameters::Get()->GetPMTBackSpecSpike());
  std::vector<G4double> myPMTBackBackscat(mysize, DSParameters::Get()->GetPMTBackBackscat());
  std::vector<G4double> myPMTBackRindex(mysize, DSParameters::Get()->GetPMTBackRindex());

  std::vector<G4double> myUnssEnergy(myLumirrorEnergy);
  std::vector<G4double> myUnssRefl(mysize, DSParameters::Get()->GetUnssRefl());
  std::vector<G4double> myUnssEff(mysize, DSParameters::Get()->GetUnssEff());
  std::vector<G4double> myUnssSpecLobe(mysize, DSParameters::Get()->GetUnssSpecLobe());
  std::vector<G4double> myUnssSpecSpike(mysize, DSParameters::Get()->GetUnssSpecSpike());
  std::vector<G4double> myUnssBackscat(mysize, DSParameters::Get()->GetUnssBackscat());
  std::vector<G4double> myUnssRindex(mysize, DSParameters::Get()->GetUnssRindex());

  fElectropolishedStainlessSteelMPT = new G4MaterialPropertiesTable();
  fElectropolishedStainlessSteelMPT->AddProperty("SPECULARLOBECONSTANT", &myEpssEnergy[0], &myEpssSpecLobe[0], mysize);
  fElectropolishedStainlessSteelMPT->AddProperty("SPECULARSPIKECONSTANT", &myEpssEnergy[0], &myEpssSpecSpike[0], mysize);
  fElectropolishedStainlessSteelMPT->AddProperty("BACKSCATTERCONSTANT", &myEpssEnergy[0], &myEpssBackscat[0], mysize);
  fElectropolishedStainlessSteelMPT->AddProperty("REFLECTIVITY", &myEpssEnergy[0], &myEpssRefl[0], mysize);
  fElectropolishedStainlessSteelMPT->AddProperty("EFFICIENCY", &myEpssEnergy[0], &myEpssEff[0], mysize);
  fElectropolishedStainlessSteelMPT->AddProperty("RINDEX", &myEpssEnergy[0], &myEpssRindex[0], mysize);

  fPMTBackMPT = new G4MaterialPropertiesTable();
  fPMTBackMPT->AddProperty("SPECULARLOBECONSTANT", &myPMTBackEnergy[0], &myPMTBackSpecLobe[0], mysize);
  fPMTBackMPT->AddProperty("SPECULARSPIKECONSTANT", &myPMTBackEnergy[0], &myPMTBackSpecSpike[0], mysize);
  fPMTBackMPT->AddProperty("BACKSCATTERCONSTANT", &myPMTBackEnergy[0], &myPMTBackBackscat[0], mysize);
  fPMTBackMPT->AddProperty("REFLECTIVITY", &myPMTBackEnergy[0], &myPMTBackRefl[0], mysize);
  fPMTBackMPT->AddProperty("EFFICIENCY", &myPMTBackEnergy[0], &myPMTBackEff[0], mysize);
  fPMTBackMPT->AddProperty("RINDEX", &myPMTBackEnergy[0], &myPMTBackRindex[0], mysize);

  fUntreatedStainlessSteelMPT = new G4MaterialPropertiesTable();
  fUntreatedStainlessSteelMPT->AddProperty("SPECULARLOBECONSTANT", &myUnssEnergy[0], &myUnssSpecLobe[0], mysize);
  fUntreatedStainlessSteelMPT->AddProperty("SPECULARSPIKECONSTANT", &myUnssEnergy[0], &myUnssSpecSpike[0], mysize);
  fUntreatedStainlessSteelMPT->AddProperty("BACKSCATTERCONSTANT", &myUnssEnergy[0], &myUnssBackscat[0], mysize);
  fUntreatedStainlessSteelMPT->AddProperty("REFLECTIVITY", &myUnssEnergy[0], &myUnssRefl[0], mysize);
  fUntreatedStainlessSteelMPT->AddProperty("EFFICIENCY", &myUnssEnergy[0], &myUnssEff[0], mysize);
  fUntreatedStainlessSteelMPT->AddProperty("RINDEX", &myUnssEnergy[0], &myUnssRindex[0], mysize);

  //-------------------
  // Veto PMT Photo Cathode
  // ------------------

  std::vector<G4double> myVCathodeEnergy(DSParameters::Get()->GetVCathodeEnergy());
  std::vector<G4double> myVCathodeReflec(DSParameters::Get()->GetVCathodeReflec());
  const size_t vcatsize(myVCathodeEnergy.size());
  std::vector<G4double> myVCathodeEff(vcatsize, DSParameters::Get()->GetVCathodeEff());
  std::vector<G4double> myVCathodeSpecLobe(vcatsize, DSParameters::Get()->GetVCathodeSpecLobe());
  std::vector<G4double> myVCathodeSpecSpike(vcatsize, DSParameters::Get()->GetVCathodeSpecSpike());
  std::vector<G4double> myVCathodeBackscat(vcatsize, DSParameters::Get()->GetVCathodeBackscat());
  std::vector<G4double> myVCathodeRindex(vcatsize, DSParameters::Get()->GetVCathodeRindex());

  fVPhotocathodeMPT = new G4MaterialPropertiesTable();
  fVPhotocathodeMPT->AddConstProperty("DOTRANSMISSION", 1, true);
  fVPhotocathodeMPT->AddProperty("REFLECTIVITY", &myVCathodeEnergy[0], &myVCathodeReflec[0], vcatsize);
  fVPhotocathodeMPT->AddProperty("EFFICIENCY", &myVCathodeEnergy[0], &myVCathodeEff[0], vcatsize);
  fVPhotocathodeMPT->AddProperty("SPECULARLOBECONSTANT", &myVCathodeEnergy[0], &myVCathodeSpecLobe[0], vcatsize);
  fVPhotocathodeMPT->AddProperty("SPECULARSPIKECONSTANT", &myVCathodeEnergy[0], &myVCathodeSpecSpike[0], vcatsize);
  fVPhotocathodeMPT->AddProperty("BACKSCATTERCONSTANT", &myVCathodeEnergy[0], &myVCathodeBackscat[0], vcatsize);
  fVPhotocathodeMPT->AddProperty("RINDEX", &myVCathodeEnergy[0], &myVCathodeRindex[0], vcatsize);

  //------------------------------
  // TPC PMT Photo Cathode
  //------------------------------

  fPhotocathodeMPT = new G4MaterialPropertiesTable();
  fPhotocathodeMPT->AddProperty("SPECULARLOBECONSTANT", &DSParameters::Get()->GetPhotocathodeEnergy()[0], &DSParameters::Get()->GetPhotocathodeSpecularLobe()[0], DSParameters::Get()->GetQENumEntries());
  fPhotocathodeMPT->AddProperty("SPECULARSPIKECONSTANT", &DSParameters::Get()->GetPhotocathodeEnergy()[0], &DSParameters::Get()->GetPhotocathodeSpecularSpike()[0], DSParameters::Get()->GetQENumEntries());
  fPhotocathodeMPT->AddProperty("BACKSCATTERCONSTANT", &DSParameters::Get()->GetPhotocathodeEnergy()[0], &DSParameters::Get()->GetPhotocathodeBackscatter()[0], DSParameters::Get()->GetQENumEntries());
  fPhotocathodeMPT->AddProperty("REFLECTIVITY", &DSParameters::Get()->GetPhotocathodeEnergy()[0], &DSParameters::Get()->GetPhotocathodeReflectance()[0], DSParameters::Get()->GetQENumEntries());
  fPhotocathodeMPT->AddProperty("EFFICIENCY", &DSParameters::Get()->GetPhotocathodeEnergy()[0], &DSParameters::Get()->GetPhotocathodeEfficiency()[0], DSParameters::Get()->GetQENumEntries());
  fPhotocathodeMPT->AddProperty("RINDEX", &DSParameters::Get()->GetPhotocathodeEnergy()[0], &DSParameters::Get()->GetPhotocathodeRindex()[0], DSParameters::Get()->GetQENumEntries());
}

//---------------------------------------------------------------------------------------

void DSMaterial::SetMaterialIndexes() {
  DSLog(development) << "----------------------------------------" << endlog;
  DSLog(development) << "              Material indexes          " << endlog;
  DSLog(development) << "----------------------------------------" << endlog;
  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
  for (int i = 0; i < (int)G4Material::GetNumberOfMaterials(); ++i) {
    G4Material* aMaterial = (*theMaterialTable)[i];
    DSLog(routine) << aMaterial->GetName() << ": " << (int)aMaterial->GetIndex() << endlog;

    if (aMaterial->GetName() == "BoronScintillator") DSStorage::Get()->SetBoronScintillatorIndex((int)aMaterial->GetIndex());
    else if (aMaterial->GetName() == "LiquidArgon") {
      DSStorage::Get()->SetLiquidArgonIndex((int)aMaterial->GetIndex());
      fTPCMaterialIndex = (int)aMaterial->GetIndex();
    }
    else if (aMaterial->GetName() == "IVLiquidArgon")   fIVMaterialIndex = (int)aMaterial->GetIndex();
    else if (aMaterial->GetName() == "OVLiquidArgon")   fOVMaterialIndex = (int)aMaterial->GetIndex();
    else if (aMaterial->GetName() == "XeLiquidArgon" && DSStorage::Get()->GetDS20kTPCedge() > 2.5 * m) DSStorage::Get()->SetLiquidArgonIndex((int)aMaterial->GetIndex());
    else if (aMaterial->GetName() == "GaseousArgon")    DSStorage::Get()->SetGaseousArgonIndex((int)aMaterial->GetIndex());
    else if (aMaterial->GetName() == "LArAboveGrid")    DSStorage::Get()->SetLArAboveGridIndex((int)aMaterial->GetIndex());

  }
  DSLog(development) << "----------------------------------------" << endlog;


  DSEventHandler::Get()->SetLArIndex(DSStorage::Get()->GetLiquidArgonIndex());
  if (DSStorage::Get()->GetScintillator() == 0) DSEventHandler::Get()->SetScintillatorIndex(GetBoronScintillator()->GetIndex());
  else if (DSStorage::Get()->GetScintillator() == 1)
    DSEventHandler::Get()->SetScintillatorIndex(GetGdScintillator()->GetIndex());
  else if (DSStorage::Get()->GetScintillator() == 2)
    DSEventHandler::Get()->SetScintillatorIndex(GetLi6Scintillator()->GetIndex());
  else if (DSStorage::Get()->GetScintillator() == 3)
    DSEventHandler::Get()->SetScintillatorIndex(GetNatLi6Scintillator()->GetIndex());
  else if (DSStorage::Get()->GetScintillator() == 4)
    DSEventHandler::Get()->SetScintillatorIndex(GetOrtoCarbScintillator()->GetIndex());
}

//---------------------------------------------------------------------------------------

G4double DSMaterial::GetLArRefIndex(G4double lambda) {

  G4double epsilon;
  const G4double LArRho = 1.390;
  const G4double ArRho = 0.001784;
  if (lambda <= 107.05) return 1.0e4;   // lambda MUST be > 107.05 nm
  epsilon = lambda / 1000.0;            // switch to micrometers
  epsilon = 1.0 / (epsilon * epsilon);  // 1 / (lambda)^2
  epsilon = 1.2055e-2 * (0.2075 / (91.012 - epsilon) + 0.0415 / (87.892 - epsilon) + 4.3330 / (214.02 - epsilon));
  epsilon *= (8. / 12.);        // Bideau-Sellmeier -> Clausius-Mossotti
  epsilon *= (LArRho / ArRho);  // density correction (Ar gas -> LAr liquid)
  if ((epsilon < 0.0) || (epsilon > 0.999999)) return 4.0e6;
  epsilon = (1.0 + 2.0 * epsilon) / (1.0 - epsilon);  // solve Clausius-Mossotti

  return sqrt(epsilon);
}
//---------------------------------------------------------------------------------------

G4double DSMaterial::GetGArRefIndex(G4double lambda) {

  G4double epsilon;

  if (lambda <= 107.05) return 1.0e4;   // lambda MUST be > 107.05 nm
  epsilon = lambda / 1000.0;            // switch to micrometers
  epsilon = 1.0 / (epsilon * epsilon);  // 1 / (lambda)^2
  epsilon = 1.2055e-2 * (0.2075 / (91.012 - epsilon) + 0.0415 / (87.892 - epsilon) + 4.3330 / (214.02 - epsilon));
  epsilon *= (8. / 12.);  // Bideau-Sellmeier -> Clausius-Mossotti
  if ((epsilon < 0.0) || (epsilon > 0.999999)) return 4.0e6;
  epsilon = (1.0 + 2.0 * epsilon) / (1.0 - epsilon);  // solve Clausius-Mossotti

  return sqrt(epsilon);
}

//---------------------------------------------------------------------------------------

G4double DSMaterial::GetLArRayLength(G4double lambda) {
  const G4double LArT = 87.0;                                    // the actual temperature of LAr in detector
  const G4double LArKT = 2.18e-10 * cm * cm / (1.e-5 * newton);  // LAr isothermal compressibility
  const G4double k = 1.380658e-23;                               // the Boltzmann constant

  G4double h = GetLArEpsilon(lambda);
  if (h < 1.00000001) h = 1.00000001;                // just a precaution
  h = (h - 1.0) * (h + 2.0);                         // the "dielectric constant" dependance
  h *= h;                                            // take the square
  h *= LArKT * LArT * k;                             // compressibility * temperature * Boltzmann constant
  h /= lambda * lambda * lambda * lambda * 1.0e-36;  // (lambda)^4
  h *= 9.18704494231105429;                          // (2 * Pi / 3)^3
  //   if ( h < (1.0 / (10.0 * km)) ) h = 1.0 / (10.0 * km); // just a
  //   precaution if ( h > (1.0 / (0.1 * nanometer)) ) h = 1.0 / (0.1 *
  //   nanometer); // just a precaution
  //  return ( 100.0 / h );
  return (100. / h) * um;
}

//---------------------------------------------------------------------------------------

G4double DSMaterial::GetGArRayLength(G4double lambda) {
  const G4double ArT = 87.0;                                   // the actual temperature of LAr in detector
  const G4double ArKT = 2.18e-5 * cm * cm / (1.e-5 * newton);  // LAr isothermal compressibility
  const G4double k = 1.380658e-23;                             // the Boltzmann constant

  G4double h = GetGArEpsilon(lambda);
  if (h < 1.00000001) h = 1.00000001;                // just a precaution
  h = (h - 1.0) * (h + 2.0);                         // the "dielectric constant" dependance
  h *= h;                                            // take the square
  h *= ArKT * ArT * k;                               // compressibility * temperature * Boltzmann constant
  h /= lambda * lambda * lambda * lambda * 1.0e-36;  // (lambda)^4
  h *= 9.18704494231105429;                          // (2 * Pi / 3)^3
  //   if ( h < (1.0 / (10.0 * km)) ) h = 1.0 / (10.0 * km); // just a
  //   precaution if ( h > (1.0 / (0.1 * nanometer)) ) h = 1.0 / (0.1 *
  //   nanometer); // just a precaution
  return (100.0 / h) * um;
}

//---------------------------------------------------------------------------------------
G4double DSMaterial::GetLArEpsilon(G4double lambda) {
  const G4double LArRho = 1.390;
  const G4double ArRho = 0.001784;
  G4double epsilon;

  if (lambda <= 107.05) return 1.0e4;   // lambda MUST be > 107.05 nm
  epsilon = lambda / 1000.0;            // switch to micrometers
  epsilon = 1.0 / (epsilon * epsilon);  // 1 / (lambda)^2
  epsilon = 1.2055e-2 * (0.2075 / (91.012 - epsilon) + 0.0415 / (87.892 - epsilon) + 4.3330 / (214.02 - epsilon));
  epsilon *= (8. / 12.);        // Bideau-Sellmeier -> Clausius-Mossotti
  epsilon *= (LArRho / ArRho);  // density correction (Ar gas -> LAr liquid)
  if ((epsilon < 0.0) || (epsilon > 0.999999)) return 4.0e6;
  epsilon = (1.0 + 2.0 * epsilon) / (1.0 - epsilon);  // solve Clausius-Mossotti

  return epsilon;
}
//---------------------------------------------------------------------------------------
G4double DSMaterial::GetGArEpsilon(G4double lambda) {

  G4double epsilon;

  if (lambda <= 107.05) return 1.0e4;   // lambda MUST be > 107.05 nm
  epsilon = lambda / 1000.0;            // switch to micrometers
  epsilon = 1.0 / (epsilon * epsilon);  // 1 / (lambda)^2
  epsilon = 1.2055e-2 * (0.2075 / (91.012 - epsilon) + 0.0415 / (87.892 - epsilon) + 4.3330 / (214.02 - epsilon));
  epsilon *= (8. / 12.);  // Bideau-Sellmeier -> Clausius-Mossotti
  if ((epsilon < 0.0) || (epsilon > 0.999999)) return 4.0e6;
  epsilon = (1.0 + 2.0 * epsilon) / (1.0 - epsilon);  // solve Clausius-Mossotti

  return epsilon;
}
