#include "DSBiasStopStep.hh"
#include "DSBiasStopStepMessenger.hh"
#include "DSLogger.hh"
#include "G4SystemOfUnits.hh"
#include "G4Step.hh"
#include <string.h>
#include "DSMaterial.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"

using namespace std;

////////////////////////////////////////////////////////////////////////
//
// Author: Davide Franco (davide.franco@apc.in2p3.fr)
//
// Date: 16/11/23
////////////////////////////////////////////////////////////////////////





DSBiasStopStep::DSBiasStopStep() {
  fBiasType         = BiasGeometryEnum::None;
  fMaterialNameOut  = "TMB";
  fMaterialNameIn   = "ITO";
  fMaterialIdxOut   = -999 ;
  fMaterialIdxIn    = -999 ;
  fRadius           = 2*m;
  fSide             = 2*m;
  fHalfZ            = 2*m;
  fCenter           = G4ThreeVector(0,0,0);

  fEPSILON          = 1e-5;


  fMessenger = new DSBiasStopStepMessenger(this);

  // need to call a material (whatever) to construct the material table
  //DSMaterial::Get()->GetPPO();

  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();

  for (int i = 0; i < (int) G4Material::GetNumberOfMaterials(); ++i) {
    G4Material* aMaterial = (*theMaterialTable)[i];
    if (aMaterial->GetName() == fMaterialNameIn)  fMaterialIdxIn  = (int) aMaterial->GetIndex();
    if (aMaterial->GetName() == fMaterialNameOut) fMaterialIdxOut = (int) aMaterial->GetIndex();
  }

}

G4bool DSBiasStopStep::CheckIfMaterial(const G4Step *step) {
  bool isGood = false ;
  if ( (int) step->GetPostStepPoint()->GetMaterial()->GetIndex() == fMaterialIdxOut) isGood = true ;
  if (  fMaterialIdxIn != -999 &&  (int) step->GetPreStepPoint()->GetMaterial()->GetIndex() != fMaterialIdxIn) isGood = false;
  return isGood ;
}

G4bool DSBiasStopStep::CheckIfSphere(const G4Step *step) {
  float _r1 = (step->GetPreStepPoint()->GetPosition()-fCenter).mag();
  float _r2 = (step->GetPostStepPoint()->GetPosition()-fCenter).mag();

  return ((_r1 > fRadius && _r2 < fRadius) || (_r1 < fRadius && _r2 > fRadius));
}

G4bool DSBiasStopStep::CheckIfCube(const G4Step *step) {
  G4ThreeVector _v1 = step->GetPreStepPoint()->GetPosition() - fCenter ;
  G4ThreeVector _v2 = step->GetPostStepPoint()->GetPosition() - fCenter ;
  bool _v1_in  = fabs(_v1.x()) < fSide/2 && fabs(_v1.y()) < fSide/2 &&  fabs(_v1.z()) < fSide/2;
  bool _v1_out = fabs(_v1.x()) > fSide/2 or fabs(_v1.y()) > fSide/2 or  fabs(_v1.z()) > fSide/2;
  bool _v2_in  = fabs(_v2.x()) < fSide/2 && fabs(_v2.y()) < fSide/2 &&  fabs(_v2.z()) < fSide/2;
  bool _v2_out = fabs(_v2.x()) > fSide/2 or fabs(_v2.y()) > fSide/2 or  fabs(_v2.z()) > fSide/2;
  return (_v1_in && _v2_out) || (_v2_in && _v1_out);
}

G4bool DSBiasStopStep::CheckIfCylinder(const G4Step *step) {
  G4ThreeVector _v1 = step->GetPreStepPoint()->GetPosition() - fCenter ;
  G4ThreeVector _v2 = step->GetPostStepPoint()->GetPosition() - fCenter ;
  bool _v1_in  = pow(_v1.x(),2) + pow(_v1.y(),2) < fRadius*fRadius && fabs(_v1.z()) < fHalfZ;
  bool _v2_in  = pow(_v2.x(),2) + pow(_v2.y(),2) < fRadius*fRadius && fabs(_v2.z()) < fHalfZ;
  bool _v1_out = pow(_v1.x(),2) + pow(_v1.y(),2) > fRadius*fRadius || fabs(_v1.z()) > fHalfZ;
  bool _v2_out = pow(_v2.x(),2) + pow(_v2.y(),2) > fRadius*fRadius || fabs(_v2.z()) > fHalfZ;
  return (_v1_in && _v2_out) || (_v2_in && _v1_out);
}


G4bool DSBiasStopStep::CheckIfCrosses(const G4Step *step) {
  bool intercepted = false;
  switch (fBiasType) {
    case BiasGeometryEnum::None:
      break ;
    case BiasGeometryEnum::Material:
      intercepted = CheckIfMaterial(step) ;
      break ;
    case BiasGeometryEnum::Sphere:
      intercepted = CheckIfSphere(step) ;
      break ;
    case BiasGeometryEnum::Cube:
      intercepted = CheckIfCube(step) ;
      break ;
    case BiasGeometryEnum::Cylinder:
      intercepted = CheckIfCylinder(step) ;
      break ;
  }
  return intercepted ;

}


G4ThreeVector DSBiasStopStep::GetInteceptedPosition(const G4Step *step) {
  G4ThreeVector pos = step->GetPostStepPoint()->GetPosition();
  switch (fBiasType) {
    case BiasGeometryEnum::None:
      break;
    case BiasGeometryEnum::Material:
      break;
    case BiasGeometryEnum::Sphere:
      pos = InterceptSphere(step) ;
      break ;
    case BiasGeometryEnum::Cube:
      pos = InterceptCube(step) ;
      break ;
    case BiasGeometryEnum::Cylinder:
      pos = InterceptCylinder(step) ;
      break ;
  }
  return pos ;

}


G4ThreeVector DSBiasStopStep::InterceptSphere(const G4Step* step) {

  G4ThreeVector v1 =  step->GetPreStepPoint()->GetPosition() - fCenter;
  G4ThreeVector v2 =  step->GetPostStepPoint()->GetPosition() - fCenter;
  G4ThreeVector dist = v2-v1;

  G4ThreeVector out = v1 ;

  float a = dist.mag2();
  float b = dist.dot(v1);
  float delta = pow(2*b,2) - 4*a*(v1.mag2() - fRadius*fRadius);
  if (delta >= 0) {
    float u1 = (-2 * b - sqrt(delta))/(2*a);
    float u2 = (-2 * b + sqrt(delta))/(2*a);
    G4ThreeVector sol1 = v1 + u1*dist;
    G4ThreeVector sol2 = v1 + u2*dist;
    float dist1 = (sol1-v1).mag2() + (sol1-v2).mag2();
    float dist2 = (sol2-v1).mag2() + (sol2-v2).mag2();

    if (dist1 < dist2) out = sol1 ;
    else out = sol2 ;
  }

  return out + fCenter;
}




G4ThreeVector DSBiasStopStep::InterceptCube(const G4Step* step) {

  G4ThreeVector v1 =  step->GetPreStepPoint()->GetPosition()  - fCenter;
  G4ThreeVector v2 =  step->GetPostStepPoint()->GetPosition() - fCenter;
  G4ThreeVector out = G4ThreeVector(0,0,0);

  float _m  = 0;
  for(int i=0;i<3;++i) {
    if (((v1[i] < fSide/2.) && (v2[i] > fSide/2.)) || ((v1[i] > fSide/2.) && (v2[i] < fSide/2.))) {

      _m = (fSide/2. - v1[i])/(v2[i] - v1[i]);
      out = _m*(v2-v1) + v1;

      int intercepted = 1;
      for(int j=0;j<3;++j) {
        if(i==j) continue ;
        if (fabs(out[j]) > fSide/2.) intercepted = 0;
      }
      if (intercepted ==1)  break;


    } else if  (((v1[i] < -fSide/2.) && (v2[i] > -fSide/2.)) || ((v1[i] > -fSide/2.) && (v2[i] < -fSide/2.))) {

      _m = (-fSide/2. - v1[i])/(v2[i] - v1[i]);
      out = _m*(v2-v1) + v1;
      int intercepted = 1;
      for(int j=0;j<3;++j) {
        if(i==j) continue ;
        if (fabs(out[j]) > fSide/2.) intercepted = 0;
      }
      if (intercepted ==1)  break;
    }
  }

  return out + fCenter;
}



G4ThreeVector DSBiasStopStep::InterceptCylinder(const G4Step* step) {

  G4ThreeVector v1   = step->GetPreStepPoint()->GetPosition()  - fCenter ;
  G4ThreeVector v2   = step->GetPostStepPoint()->GetPosition() - fCenter ;
  G4ThreeVector dist = v2-v1;
  G4ThreeVector out  = G4ThreeVector(0,0,0);

  float _m  = 0;

  if ((v1.z() > fHalfZ && v2.z() < fHalfZ ) || (v1.z() < fHalfZ && v2.z() > fHalfZ )) {
    _m = (fHalfZ -  v1.z())/(v2.z() - v1.z());
  } else if ((v1.z() > -fHalfZ && v2.z() < -fHalfZ ) || (v1.z() < -fHalfZ && v2.z() > -fHalfZ )) {
    _m = (-fHalfZ -  v1.z())/(v2.z() - v1.z());
  }
  out = G4ThreeVector(_m*(v2.x()-v1.x()) + v1.x(), _m*(v2.y()-v1.y()) + v1.y(), _m*(v2.z()-v1.z())+ v1.z());
  if(_m != 0 && pow(out.x(),2)+pow(out.y(),2) < fRadius*fRadius) return out + fCenter;

  G4ThreeVector v1_0 = G4ThreeVector(v1.x(),v1.y(),0);
  G4ThreeVector v2_0 = G4ThreeVector(v2.x(),v2.y(),0);
  float t = v1_0.dot(v1_0-v2_0)/ (v1_0-v2_0).dot(v1_0-v2_0);
  float d2 = v1_0.mag2() - pow(t,2)*(v1_0-v2_0).dot(v1_0-v2_0);
  float delta2  = (fRadius * fRadius -d2)/(v1_0-v2_0).dot(v1_0-v2_0);

  G4ThreeVector sol1 = v1 + (t+sqrt(delta2))*(v2-v1);
  G4ThreeVector sol2 = v1 + (t-sqrt(delta2))*(v2-v1);

  float dist1 = (sol1-v1).mag2() + (sol1-v2).mag2();
  float dist2 = (sol2-v1).mag2() + (sol2-v2).mag2();

  if (dist1 < dist2) out = sol1 ;
  else out = sol2 ;

  return out + fCenter;

}


void DSBiasStopStep::SetBiasGeometryType(G4String val) {

  if      (strcmp(val, "Sphere")   == 0) fBiasType = BiasGeometryEnum::Sphere;
  else if (strcmp(val, "Cylinder") == 0) fBiasType = BiasGeometryEnum::Cylinder;
  else if (strcmp(val, "Cube")     == 0) fBiasType = BiasGeometryEnum::Cube;
  else if (strcmp(val, "Material") == 0) fBiasType = BiasGeometryEnum::Material ;
  else fBiasType = BiasGeometryEnum::None ;

  DSLog(routine) << "Bias: stop particles at " << val <<  " (" << int(fBiasType) << ")" << endlog ;

}
