#ifndef DSDetectorNeutronReD_H
#define DSDetectorNeutronReD_H

#include "G4LogicalVolume.hh"
#include "G4Polycone.hh"
#include "G4SubtractionSolid.hh"
#include "G4Tubs.hh"
#include "G4UnionSolid.hh"
#include "G4VPhysicalVolume.hh"

#include "G4LogicalBorderSurface.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4OpticalSurface.hh"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
// class LSci: describes position and pointing, as taken from survey file     //
////////////////////////////////////////////////////////////////////////////////

class LSci {
 private:
  G4ThreeVector pos;
  G4String type;
  G4RotationMatrix* rot;

 public:
  LSci() : pos(G4ThreeVector()), type(""), rot(0){};
  void SetPosition(G4ThreeVector p) { pos = p; }
  void SetType(G4String t) { type = t; }
  void SetRotation(G4RotationMatrix* r) { rot = r; }
  const G4ThreeVector GetPosition() const { return pos; }
  const G4String GetType() const { return type; }
  /*const*/ G4RotationMatrix* GetRotation() const { return rot; }
};
inline std::ostream& operator<<(std::ostream& out, const LSci& lsci) {
  out << "pos " << lsci.GetPosition() << " type ";
  if (lsci.GetType() != "") out << lsci.GetType();
  else
    out << "-";
  out << " rot ";
  const G4RotationMatrix* dummy = lsci.GetRotation();
  if (dummy)  // pointer exists
    out << *dummy;
  else  // pointer 0
    out << dummy;
  return out;
}

////////////////////////////////////////////////////////////////////////////////
// class DSDetectorNeutronReD: reads the survey file and prepared the setup   //
////////////////////////////////////////////////////////////////////////////////

class DSDetectorNeutronReD {

 public:
  DSDetectorNeutronReD(G4VPhysicalVolume*);
  ~DSDetectorNeutronReD() { /*delete fMessenger;*/
  }

 private:
  void parseSurveyFile(G4String surveyFileName);
  //  void getPositionFromSphericalSurvey(G4String surveyFileName,
  //  vector<G4ThreeVector>& ND_position);
  void DefineSurfaces() {
    ////////////////////////////////////////
    // Copied from DS DetectorDS50 - 24th april 2015
    ////////////////////////////////////////
    return;
  }

  G4VPhysicalVolume* fMotherVolume;
  std::vector<LSci> lsciCol;
};

////////////////////////////////////////////////////////////////////////////////

#endif
