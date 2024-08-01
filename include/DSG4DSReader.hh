//
//	DSG4DSReader.hh
//

#ifndef DSG4DSReader_h
#define DSG4DSReader_h 1
//#include "DSInputVertex.hh"
//#include "DSOutputStructure.hh"
#include <stdio.h>
#include <iostream>
#include "DSEventStructure.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

using namespace std;

class DSG4DSReader {
 private:
  DSG4DSReader();

 public:
  ~DSG4DSReader();
  static DSG4DSReader* Get();
  void ReadHeader();
  G4bool ReadEvent();
  void DumpHeader();
  void SkipEvents(G4int);
  void DumpEvent();
  void SetG4DSFile();

  G4double GetPrimaryEventEnergy () ; 
  G4ThreeVector GetPrimaryPosition () ; 
  G4int GetPrimEvID () ; 

  void ClearAll();

  // vector<DepositStructure>&           GetVDeposits() { return theDeposits ; }
  vector<DepositStructure> GetVDeposits() { return theDeposits; }
  vector<DaughterStructure> GetVDaughters() { return theDaughters; }

  DepositStructure& GetDeposits() { return theDepositStructure; }
  DaughterStructure& GetDaughters() { return theDaughterStructure; }

  void SetDeposit(DepositStructure val) { theDepositStructure = val; }
  void SetDeposits() { theDeposits.push_back(theDepositStructure); }

  void SetDaughter(DaughterStructure val) { theDaughterStructure = val; }
  void SetDaughters() { theDaughters.push_back(theDaughterStructure); }
  
  EventStructureDiskFormat GetEvent() { return fEvent; }

 private:
  static DSG4DSReader* me;

  void SkipEvents();

  EventStructureDiskFormat fEvent;
  DaughterStructure fDaughter;
  DepositStructure fDeposit;
  UserStructure fUser;
  PhotonStructure fPhoton;
  PhotoElectronStructure fPhotoElectron;
  PhotoElectronStructure fMuPhotoElectron;
  PhotoElectronStructure fVetoPhotoElectron;

  DepositStructure  theDepositStructure;
  DaughterStructure theDaughterStructure;
  vector<DepositStructure> theDeposits;
  vector<DaughterStructure> theDaughters;
};

#endif
/*
 * $Log: DSG4DSReader.hh,v $
 * Revision 1.1  2014/05/07 12:20:52  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.1  2014/03/11 09:56:25  meregaglia
 * Added generator starting from energy deposits
 *
 * Revision 1.1  2009-10-20 12:37:38  dfranco
 * Added a new generator in order to read a g4ds output file and begin a
 * simulation from the energy deposits. Useful for external background.
 *
 *
 */
