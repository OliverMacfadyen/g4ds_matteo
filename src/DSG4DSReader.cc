//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
//
//
//	DSG4DSReader.cc
//
#include "DSG4DSReader.hh"
#include "DSEventHandler.hh"
#include "DSEventStructure.hh"
#include "DSIO.hh"
#include "DSLogger.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
//#include "DSEventHandler.hh"
using namespace std;

DSG4DSReader* DSG4DSReader::me = 0;

// singleton
DSG4DSReader::DSG4DSReader() {}

DSG4DSReader* DSG4DSReader::Get() {
  if (!me) me = new DSG4DSReader();
  return me;
}

G4bool DSG4DSReader::ReadEvent() {

  int BuffDimension;
  DSIO::Get()->GetG4DSFile().read(reinterpret_cast<char*>(&BuffDimension), sizeof(int));

  DSIO::Get()->GetG4DSFile().read(reinterpret_cast<char*>(&fEvent), sizeof(EventStructureDiskFormat));
  //DSEventHandler::Get()->SetEventID(fEvent.EventID);
  //DSEventHandler::Get()->SetEnergy(fEvent.Energy);
  //DSEventHandler::Get()->SetTime(fEvent.Time);
  //DSEventHandler::Get()->SetPDG(fEvent.PDG);
 
  for (int i = 0; i < fEvent.NDaughters; i++) {
    DSIO::Get()->GetG4DSFile().read(reinterpret_cast<char*>(&fDaughter), sizeof(DaughterStructure));
    fDaughter.Energy *= keV;
    DSEventHandler::Get()->SetDaughter(fDaughter);
    DSEventHandler::Get()->SetDaughters();
    SetDaughter(fDaughter);
    SetDaughters();
  }

  for (int i = 0; i < fEvent.NDeposits; i++) {
    DSIO::Get()->GetG4DSFile().read(reinterpret_cast<char*>(&fDeposit), sizeof(DepositStructure));
    fDeposit.Energy *= keV;
    SetDeposit(fDeposit);
    SetDeposits();
    // cout<<fDeposit.Energy*MeV<<" "<<fDeposit.Time<<endl;
    // cout<<" size "<<theDeposits.size()<<endl;
  }

  for (int i = 0; i < fEvent.NUsers; i++) {
    DSIO::Get()->GetG4DSFile().read(reinterpret_cast<char*>(&fUser), sizeof(UserStructure));
    DSEventHandler::Get()->SetUser(fUser);
    DSEventHandler::Get()->SetUsers();
  }

  for (int i = 0; i < fEvent.NPH; i++) {
    DSIO::Get()->GetG4DSFile().read(reinterpret_cast<char*>(&fPhoton), sizeof(PhotonStructure));
    DSEventHandler::Get()->SetPhoton(fPhoton);
    DSEventHandler::Get()->SetPhotons();
  }

  for (int i = 0; i < fEvent.NPE; i++) {
    DSIO::Get()->GetG4DSFile().read(reinterpret_cast<char*>(&fPhotoElectron), sizeof(PhotoElectronStructure));
    DSEventHandler::Get()->SetPhotoElectron(fPhotoElectron);
    DSEventHandler::Get()->SetPhotoElectrons();
  }

  for (int i = 0; i < fEvent.VetoNPE; i++) {
    DSIO::Get()->GetG4DSFile().read(reinterpret_cast<char*>(&fMuPhotoElectron), sizeof(PhotoElectronStructure));
    DSEventHandler::Get()->SetMuPhotoElectron(fMuPhotoElectron);
    DSEventHandler::Get()->SetMuPhotoElectrons();
  }

  for (int i = 0; i < fEvent.MuNPE; i++) {
    DSIO::Get()->GetG4DSFile().read(reinterpret_cast<char*>(&fVetoPhotoElectron), sizeof(PhotoElectronStructure));
    DSEventHandler::Get()->SetMuPhotoElectron(fVetoPhotoElectron);
    DSEventHandler::Get()->SetMuPhotoElectrons();
  }

  int BuffDimension2;
  DSIO::Get()->GetG4DSFile().read(reinterpret_cast<char*>(&BuffDimension2), sizeof(int));
  
  //if (DSLogger::GetSeverity() == DSLogger::trace) DSLog(trace) << BuffDimension <<  " !! " << BuffDimension2 << endlog ;
  
  if (BuffDimension > 0  && BuffDimension == BuffDimension2) return true;
  return false;
}
void DSG4DSReader::ReadHeader() {

  int BuffDimension;
  DSIO::Get()->GetG4DSFile().read(reinterpret_cast<char*>(&BuffDimension), sizeof(int));
  DSIO::Get()->GetG4DSFile().ignore(BuffDimension);
  int BuffDimension2;
  DSIO::Get()->GetG4DSFile().read(reinterpret_cast<char*>(&BuffDimension2), sizeof(int));
  if (BuffDimension != BuffDimension2) DSLog(error) << "Error in reading the G4DS file header!" << endl;
}

void DSG4DSReader::SkipEvents(G4int nevents) {
  int size1, size2;
  for(int i=0;i<nevents;++i) {
    DSIO::Get()->GetG4DSFile().read(reinterpret_cast<char*>(&size1), sizeof(int));
    DSIO::Get()->GetG4DSFile().ignore(size1);
    DSIO::Get()->GetG4DSFile().read(reinterpret_cast<char*>(&size2), sizeof(int));
  }

}

void DSG4DSReader::DumpHeader() {}

void DSG4DSReader::DumpEvent() {}

void DSG4DSReader::ClearAll() {
  theDeposits.clear();
  theDaughters.clear();
}

G4double DSG4DSReader::GetPrimaryEventEnergy () {
  //if (fEvent.NUsers > 0 ) return fUser.UserFloat1; 
  return fEvent.Energy  ; 
}
 
G4ThreeVector DSG4DSReader::GetPrimaryPosition () {
  G4ThreeVector vec  (fEvent.Position [0], fEvent.Position [1], fEvent.Position [2] ) ; 
  return vec  ;
}

G4int GetPrimEvID () {
  return 0; // fEvent.fEventID ; 
}


