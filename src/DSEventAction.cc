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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DSEventAction.hh"
#include <fstream>
#include <iostream>
#include <string>
#include "DSEventActionMessenger.hh"
#include "DSEventHandler.hh"
#include "DSIO.hh"
#include "DSLogger.hh"
#include "DSMaterial.hh"
#include "DSStorage.hh"
#include "DSManager.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4PhysicalConstants.hh"
#include "G4SDManager.hh"
#include "G4String.hh"
#include "G4SystemOfUnits.hh"
#include "G4Timer.hh"
#include "G4Trajectory.hh"
#include "G4TrajectoryContainer.hh"
#include "G4UImanager.hh"
#include "G4UnitsTable.hh"
#include "G4VHitsCollection.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"
#include "Randomize.hh"

#include <algorithm>
#include <vector>

using namespace std;

DSEventAction::DSEventAction() {
  fTotNPE = 0;
  fTimer = new G4Timer;
  fTimer->Start();
  fMessenger = new DSEventActionMessenger(this);
}

DSEventAction::~DSEventAction() {
  delete fMessenger;
}

void DSEventAction::BeginOfEventAction(const G4Event*) {
  // FIX ME: evt is an unused variable
  ;
}

bool DSEventAction::ArDMWriteFilter() {

  G4Material* matAir = DSMaterial::Get()->GetAir();
  G4Material* matVac = DSMaterial::Get()->GetVacuum();
  G4Material* matLAr = DSMaterial::Get()->GetLiquidArgon();
  G4Material* matPAr = DSMaterial::Get()->GetPseudoArgon();

  bool inArDM = false;

  for (int i = 0; i < DSEventHandler::Get()->GetNDeposits(); i++) {
    DepositStructure& dep = DSEventHandler::Get()->GetVDeposits()[i];

    // if deposit in ArDM or DART
    if (dep.Volume == int(matPAr->GetIndex()) || dep.Volume == int(matLAr->GetIndex())) {
      if (!inArDM) DSStorage::Get()->IncrArDMeventsShieldPE();
      if (dep.Energy > 0) return true;  // event accepted, will be written on file
      inArDM = true;
    }

    // events crossing the PE shield of ArDM
    if (!inArDM && dep.Volume != int(matAir->GetIndex()) && dep.Volume != int(matVac->GetIndex())) {
      inArDM = true;
      DSStorage::Get()->IncrArDMeventsShieldPE();
    }
  }

  return false;  // event rejected, will not to be written on file
}

bool DSEventAction::AriserWriteFilter() {
  G4Material* matBaf511 = DSMaterial::Get()->GetBariumFluoride1();
  G4Material* matBaf1270 = DSMaterial::Get()->GetBariumFluoride2();
  G4Material* matGe = DSMaterial::Get()->GetGermanium();
  G4Material* matLAr = DSMaterial::Get()->GetLiquidArgon();

  bool is_keep = false;
  for (int i = 0; i < DSEventHandler::Get()->GetNDeposits(); i++) {
    DepositStructure& dep = DSEventHandler::Get()->GetVDeposits()[i];
    if (dep.Volume == int(matBaf511->GetIndex()) || dep.Volume == int(matBaf1270->GetIndex()) || dep.Volume == int(matGe->GetIndex()) || dep.Volume == int(matLAr->GetIndex())) {
      is_keep = true;
      break;
    }
  }
  return is_keep;  // event rejected, will not to be written on file
}

bool DSEventAction::TPCEnergyWriteFilter() {

  // true: the event passes the selection and should be written
  // GetTPCDepEnergy is a float (in default units of keV);
  // GetWriteFilterMinimumEnergy are G4double with units
  if (DSEventHandler::Get()->GetTPCDepEnergy() > DSStorage::Get()->GetWriteFilterMinimumEnergy() / keV && DSEventHandler::Get()->GetTPCDepEnergy() < DSStorage::Get()->GetWriteFilterMaximumEnergy() / keV) return true;
  else
    return false;
}

bool DSEventAction::NoDaughtersWriteFilter() {

  // true: the event passes the selection and should be written
  // GetTPCDepEnergy is a float (in default units of keV);
  // GetWriteFilterMinimumEnergy are G4double with units
  if (DSEventHandler::Get()->GetNDaughters() == 0) return true;
  else
    return false;
}


void DSEventAction::EndOfEventAction(const G4Event* evt) {

  // Save event infos
  if (not DSStorage::Get()->GetOverwriteCounter()) DSEventHandler::Get()->SetEventID(evt->GetEventID());
  DSEventHandler::Get()->SetNPE(int(DSEventHandler::Get()->GetVPhotoElectrons().size()));
  DSEventHandler::Get()->SetMuNPE(int(DSEventHandler::Get()->GetVMuPhotoElectrons().size()));
  DSEventHandler::Get()->SetVetoNPE(int(DSEventHandler::Get()->GetVVetoPhotoElectrons().size()));
  DSEventHandler::Get()->SetNPH(int(DSEventHandler::Get()->GetVPhotons().size()));
  DSEventHandler::Get()->SetNDeposits(int(DSEventHandler::Get()->GetVDeposits().size()));
  DSEventHandler::Get()->SetNDaughters(int(DSEventHandler::Get()->GetVDaughters().size()));
  DSEventHandler::Get()->SetNUsers(int(DSEventHandler::Get()->GetVUsers().size()));
  DSEventHandler::Get()->SetNClusters(int(DSEventHandler::Get()->GetVClusters().size()));  //edited here to try and add a clusters size element 

  // Write info on standard output and log file
  int EventCounter = DSStorage::Get()->GetEventCounter();
  fTotNPE += DSEventHandler::Get()->GetNPE();
  fTotNPE += DSEventHandler::Get()->GetVetoNPE();
  fTotNPE += DSEventHandler::Get()->GetMuNPE();
  int OverAllNPE = DSEventHandler::Get()->GetNPE() + DSEventHandler::Get()->GetVetoNPE() + DSEventHandler::Get()->GetMuNPE();
  int ClusterSize = DSEventHandler::Get()->GetNClusters(); 

  //if (!evt->GetEventID()) {
  //  fTimer->Stop();
  //  fTimer->Start();
  //}
  if ((evt->GetEventID() % EventCounter) == 0) {
    fTimer->Stop();
    if (evt->GetEventID() == 0) {
      DSLog(routine) << ">>> Event " << evt->GetEventID() << ";  NPE = " << OverAllNPE << ";  NPE/event = " << G4float(fTotNPE) << "; Cluster Size = " << ClusterSize << ";" << endlog;
    } else {
      DSLog(routine) << ">>> Event " << evt->GetEventID() << ";  NPE = " << OverAllNPE << ";  NPE/event = " << G4float(fTotNPE) / G4float(EventCounter) << "; Cluster Size = " << ClusterSize << ";  CPUTime/event = " << fTimer->GetRealElapsed() / G4float(EventCounter) << " s" << endlog;
    }
    DSLog(trace) << "    Starting Position: " << DSEventHandler::Get()->GetPosition() << " cm" << endlog;
    DSLog(trace) << "    Energy           : " << DSEventHandler::Get()->GetEnergy() << " keV" << endlog;
    DSLog(trace) << "    PDG              : " << DSEventHandler::Get()->GetPDG() << endlog;
    if (DSStorage::Get()->GetVerbosity() > 0 && DSStorage::Get()->GetNumberOfHits()) DSLog(trace) << "    Fraction of true positions:  " << float(evt->GetEventID()) / DSStorage::Get()->GetNumberOfHits() << endlog;
    DSLog(trace) << endlog;
    // if (DSStorage::Get()->GetWriteFilter()) {
    // DSLog(trace)   << "    Events filtered: " <<
    // DSStorage::Get()->GetArDMeventsShieldPE() <<  endlog;
    DSLog(trace) << endlog;
    fTimer->Start();
    fTotNPE = 0;
  }


  // Extra info on standard output
  if (DSLogger::GetSeverity() <= DSLogger::trace) {
    if (DSStorage::Get()->GetVerbosity() > 0) DSEventHandler::Get()->DumpHeader();
    if (DSStorage::Get()->GetVerbosity() > 0) DSEventHandler::Get()->DumpEvent();
    if (DSStorage::Get()->GetVerbosity() > 1) DSEventHandler::Get()->DumpDaughter();
    if (DSStorage::Get()->GetVerbosity() > 2) DSEventHandler::Get()->DumpDeposit();
    if (DSStorage::Get()->GetVerbosity() > 2) DSEventHandler::Get()->DumpCluster();
    if (DSStorage::Get()->GetVerbosity() > 2) DSEventHandler::Get()->DumpUser();
    if (DSStorage::Get()->GetVerbosity() > 3) DSEventHandler::Get()->DumpPhotoElectron();
    if (DSStorage::Get()->GetVerbosity() > 3) DSEventHandler::Get()->DumpMuPhotoElectron();
    if (DSStorage::Get()->GetVerbosity() > 3) DSEventHandler::Get()->DumpVetoPhotoElectron();
    if (DSStorage::Get()->GetVerbosity() > 4) DSEventHandler::Get()->DumpPhoton();
  }

  // (PABLOG) Skip events if filter is active and condition is not met
  // (/ds/manager/writefilter 1) (PA22) extend this feature to all the
  // geometries
  // if (DSStorage::Get()->GetWriteFilter()) {
  //   if ((DSStorage::Get()->GetArDMGeometry() && !ArDMWriteFilter())  // ArDM keeps his own filter
  //       || (DSStorage::Get()->GetAriserGeometry() && !AriserWriteFilter())    // ARIS-ER keeps its own filter
  //       || (!TPCEnergyWriteFilter())) {
  //     // DSLog(trace)   << "    WriteFilter failed: skip event" << endlog;
  //     DSEventHandler::Get()->ClearAll();
  //     return;
  //   }
  // }

  // if (DSStorage::Get()->SkipEventsWithNoDaughters() && DSEventHandler::Get()->GetNDaughters() == 0) {
  //   DSEventHandler::Get()->ClearAll();
  //   return ;
  // }

  if (DSStorage::Get()->IsAbortRun()) DSManager::Get()->AbortRun(true);

  if(DSManager::Get()->EventSurvivesFiltering() == false) {
    DSEventHandler::Get()->ClearAll();
    return ;
  }

  // Write output on binary file
  int SIZE = sizeof(EventStructureDiskFormat) + DSEventHandler::Get()->GetNDaughters() * sizeof(DaughterStructure) + DSEventHandler::Get()->GetNDeposits() * sizeof(DepositStructure) + DSEventHandler::Get()->GetNUsers() * sizeof(UserStructure) + DSEventHandler::Get()->GetNPH() * sizeof(PhotonStructure) +
             DSEventHandler::Get()->GetNPE() * sizeof(PhotoElectronStructure) + DSEventHandler::Get()->GetVetoNPE() * sizeof(PhotoElectronStructure) + DSEventHandler::Get()->GetMuNPE() * sizeof(PhotoElectronStructure);

  DSIO::Get()->GetBinaryFile().write(reinterpret_cast<char*>(&SIZE), sizeof(int));
  DSIO::Get()->GetBinaryFile().write(reinterpret_cast<char*>(&DSEventHandler::Get()->GetEvent()), sizeof(EventStructureDiskFormat));
  for (int i = 0; i < DSEventHandler::Get()->GetNDaughters(); i++) DSIO::Get()->GetBinaryFile().write(reinterpret_cast<char*>(&DSEventHandler::Get()->GetVDaughters()[i]), sizeof(DaughterStructure));
  for (int i = 0; i < DSEventHandler::Get()->GetNDeposits(); i++) DSIO::Get()->GetBinaryFile().write(reinterpret_cast<char*>(&DSEventHandler::Get()->GetVDeposits()[i]), sizeof(DepositStructure));
  for (int i = 0; i < DSEventHandler::Get()->GetNUsers(); i++) DSIO::Get()->GetBinaryFile().write(reinterpret_cast<char*>(&DSEventHandler::Get()->GetVUsers()[i]), sizeof(UserStructure));
  for (int i = 0; i < DSEventHandler::Get()->GetNPH(); i++) DSIO::Get()->GetBinaryFile().write(reinterpret_cast<char*>(&DSEventHandler::Get()->GetVPhotons()[i]), sizeof(PhotonStructure));
  for (int i = 0; i < DSEventHandler::Get()->GetNPE(); i++) DSIO::Get()->GetBinaryFile().write(reinterpret_cast<char*>(&DSEventHandler::Get()->GetVPhotoElectrons()[i]), sizeof(PhotoElectronStructure));
  for (int i = 0; i < DSEventHandler::Get()->GetMuNPE(); i++) DSIO::Get()->GetBinaryFile().write(reinterpret_cast<char*>(&DSEventHandler::Get()->GetVMuPhotoElectrons()[i]), sizeof(PhotoElectronStructure));
  for (int i = 0; i < DSEventHandler::Get()->GetVetoNPE(); i++) DSIO::Get()->GetBinaryFile().write(reinterpret_cast<char*>(&DSEventHandler::Get()->GetVVetoPhotoElectrons()[i]), sizeof(PhotoElectronStructure));
  for (int i = 0; i < DSEventHandler::Get()->GetNClusters(); i++) DSIO::Get()->GetBinaryFile().write(reinterpret_cast<char*>(&DSEventHandler::Get()->GetVClusters()[i]), sizeof(ClusterStructure));
  DSIO::Get()->GetBinaryFile().write(reinterpret_cast<char*>(&SIZE), sizeof(int));

  DSEventHandler::Get()->ClearAll();



}
