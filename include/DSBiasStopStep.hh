#ifndef DSBIASSTOPSTEP_H
#define DSBIASSTOPSTEP_H

#include "G4ThreeVector.hh"
#include "G4UserTrackingAction.hh"

class G4Step;

class DSBiasStopStepMessenger ;

enum class BiasGeometryEnum { None=0, Material=1, Cylinder=2, Sphere=3, Cube=4 }; 

class DSBiasStopStep  {

 public:
  DSBiasStopStep();
  ~DSBiasStopStep(){};


  G4bool CheckIfCrosses(const G4Step*);

  G4bool CheckIfMaterial(const G4Step*);
  G4bool CheckIfSphere(const G4Step*);
  G4bool CheckIfCube(const G4Step*);
  G4bool CheckIfCylinder(const G4Step*);

  G4ThreeVector GetInteceptedPosition(const G4Step*);


  G4ThreeVector InterceptSphere(const G4Step*);
  G4ThreeVector InterceptCube(const G4Step*);
  G4ThreeVector InterceptCylinder(const G4Step*);


  void SetMaterialNameOut(G4String val) {fMaterialNameOut = val;}
  G4String GetMaterialNameOut() {return fMaterialNameOut;}

  void SetMaterialNameIn(G4String val) {fMaterialNameIn = val;}
  G4String GetMaterialNameIn() {return fMaterialNameIn;}

  void SetRadius(G4float val) {fRadius = val;}
  G4float GetRadius() {return fRadius;}

  void SetCenter(G4ThreeVector val) {fCenter = val;}
  G4ThreeVector& GetCenter() {return fCenter;}

  void SetSide(G4float val) {fSide = val;}
  G4float GetSide() {return fSide;}

  void SetHalfZ(G4float val) {fHalfZ = val;}
  G4float GetHalfZ() {return fHalfZ;}


  void SetBiasGeometryType(G4String);
  
  
 private:
  DSBiasStopStepMessenger* fMessenger ;
  BiasGeometryEnum fBiasType ;

  G4String  fMaterialNameOut;
  G4String  fMaterialNameIn;
  G4int     fMaterialIdxOut;
  G4int     fMaterialIdxIn;
  G4float   fRadius;
  G4float   fSide;
  G4float   fHalfZ;
  G4ThreeVector fCenter;
  G4double  fEPSILON;

};

#endif
