#ifndef DSDetectorDS20kNeutronVeto_H
#define DSDetectorDS20kNeutronVeto_H

//#include "DSDetectorPMTNeutronVeto20k.hh"
#include "G4Box.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4OpticalSurface.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4SubtractionSolid.hh"
#include "G4Tubs.hh"
#include "G4UnionSolid.hh"
#include "G4VPhysicalVolume.hh"

typedef std::pair<double, double> PairFF;
typedef std::vector<std::pair<double, double>> PointCol;
typedef std::vector<std::pair<double, double>>::iterator PointColItr;
typedef std::vector<std::pair<double, double>>* PointColPtr;

class DSMaterial;

class DSDetectorDS20kNeutronVeto {

 public:
  DSDetectorDS20kNeutronVeto(G4VPhysicalVolume*);
  ~DSDetectorDS20kNeutronVeto();

  // G4VPhysicalVolume* GetDetectorComponent()  { return  fPhysicBScintillator;}
  G4VPhysicalVolume* GetDetectorComponent() { return fPhysicUAr; }

 private:
  void DefineSurfaces();

  // DSMaterial* dsmaterial;
  PointColPtr createGeometry(double, double, double, double, double, double, double, bool,int);

  G4VPhysicalVolume* fMotherVolume;

  G4VPhysicalVolume* fPhysicUAr;
  G4VPhysicalVolume* fPhysicGUAr;
  G4VPhysicalVolume* fPhysicLArBuffer;
  // G4VPhysicalVolume* fPhysicGArBuffer;
  G4VPhysicalVolume* fPhysicInnerCryo;
  G4VPhysicalVolume* fPhysicPlastic;

};

#endif
