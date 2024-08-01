#include "DSSteppingAction.hh"
#include "DSEventHandler.hh"
#include "DSLogger.hh"
#include "DSMaterial.hh"
#include "DSStorage.hh"
#include "G4PhysicalConstants.hh"
#include "G4Step.hh"
#include "G4SteppingManager.hh"
#include "DSBiasStopStep.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4UnitsTable.hh"
#include "G4ios.hh"
#include "G4Ions.hh"

using namespace std;

DSSteppingAction::DSSteppingAction() {
  fStopStep =  new DSBiasStopStep;
  //DSLog(routine) <<  fStopStep->GetStopMaterial() << endlog ;
}

void DSSteppingAction::UserSteppingAction(const G4Step* theStep) {

  bool IsTherePostStepPoint;
  if (theStep->GetTrack()->GetVolume()->GetLogicalVolume()->GetMaterial()->GetIndex() == 3) IsTherePostStepPoint = false;
  else
    IsTherePostStepPoint = true;

  //  if(DSLogger::GetSeverity() == DSLogger::development &&
  //  theStep->GetTrack()->GetVolume()->GetLogicalVolume()->GetMaterial()->Get
  //     && theStep->GetPostStepPoint()) {
  //    DSLog(development)
  if (DSLogger::GetSeverity() == DSLogger::development && theStep->GetPostStepPoint() && IsTherePostStepPoint)
    cout << " " << theStep->GetTrack()->GetDefinition()->GetParticleName() << " "
         << " " << theStep->GetTrack()->GetDefinition()->GetPDGEncoding() << " "
         << " E: " << theStep->GetTrack()->GetKineticEnergy() / keV << " keV; "
         << " Edep: " << theStep->GetTotalEnergyDeposit() / keV
         << " "
         //<<  " "
         //<<theStep->GetTrack()->GetVolume()->GetLogicalVolume()->GetMaterial()->GetName()<<
         //" "
         << " " << theStep->GetPreStepPoint()->GetPhysicalVolume()->GetLogicalVolume()->GetMaterial()->GetName() << " "
         << " " << theStep->GetPostStepPoint()->GetPhysicalVolume()->GetLogicalVolume()->GetMaterial()->GetName()
         << " "
         //<<  " " <<theStep->GetTrack()->GetVolume()->GetName()<< " "
         << " " << theStep->GetPostStepPoint()->GetPosition() / cm << " "
         << " step " << theStep->GetStepLength() / um << " "
         << " ID: " << theStep->GetTrack()->GetTrackID() << " "
         << " Parent ID: " << theStep->GetTrack()->GetParentID() << " "
         << " gtime: " << theStep->GetTrack()->GetGlobalTime()
         << " "
         // <<  " ltime: " <<theStep->GetTrack()->GetLocalTime()/ns  << " "
         // <<  " steps: " <<theStep->GetTrack()->GetCurrentStepNumber() << " "
         << theStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()
         //    << endlog ;
         << endl;
  //  }

  if (fStopStep->CheckIfCrosses(theStep) && theStep->GetPostStepPoint()->GetKineticEnergy () > 0) {

    G4ThreeVector versor = theStep->GetPostStepPoint()->GetPosition() - theStep->GetPreStepPoint()->GetPosition();
    versor /= versor.mag();

    DSEventHandler::Get()->SetDId(int(DSEventHandler::Get()->GetVDaughters().size()));
    DSEventHandler::Get()->SetDPID(theStep->GetTrack()->GetParentID());
    DSEventHandler::Get()->SetDPDG(theStep->GetTrack()->GetDefinition()->GetPDGEncoding());
    DSEventHandler::Get()->SetDProcess(0);
    DSEventHandler::Get()->SetDTime(0); // aTrack->GetGlobalTime() / ns
    DSEventHandler::Get()->SetDEnergy(theStep->GetPostStepPoint()->GetKineticEnergy ()/keV);
    DSEventHandler::Get()->SetDPosition((G4ThreeVector) fStopStep->GetInteceptedPosition(theStep) / cm);
    DSEventHandler::Get()->SetDDirection(versor);
    DSEventHandler::Get()->SetDaughters();

    theStep->GetTrack()->SetTrackStatus(fStopAndKill);
  }

  //
  if (theStep->GetTrack()->GetGlobalTime() > DSStorage::Get()->GetTimeCut()) theStep->GetTrack()->SetTrackStatus(fKillTrackAndSecondaries);

  if (theStep->GetTrack()->GetDefinition()->GetPDGEncoding() == 5 || theStep->GetTrack()->GetDefinition()->GetPDGEncoding() == 6 || theStep->GetTrack()->GetDefinition()->GetPDGEncoding() == 27 || theStep->GetTrack()->GetDefinition()->GetPDGEncoding() == 28 ||
      theStep->GetTrack()->GetDefinition()->GetPDGEncoding() == 43 || theStep->GetTrack()->GetDefinition()->GetPDGEncoding() == 44)
    theStep->GetTrack()->SetTrackStatus(fStopAndKill);

  // Write total deposited energy

  // Write deposits

  if (theStep->GetTrack()->GetDefinition()->GetPDGEncoding() != -22 && theStep->GetTotalEnergyDeposit() > 0) {

    int mat_idx = (int) theStep->GetTrack()->GetVolume()->GetLogicalVolume()->GetMaterial()->GetIndex();

    if (mat_idx == DSMaterial::Get()->GetTPCMaterialIndex() || mat_idx == DSStorage::Get()->GetLArAboveGridIndex())
          DSEventHandler::Get()->SetTPCDepEnergy(DSEventHandler::Get()->GetTPCDepEnergy() + theStep->GetTotalEnergyDeposit() / keV);
    else if (mat_idx == DSMaterial::Get()->GetIVMaterialIndex())
          DSEventHandler::Get()->SetVetoDepEnergy(DSEventHandler::Get()->GetVetoDepEnergy() + theStep->GetTotalEnergyDeposit() / keV);
    else if (mat_idx == DSMaterial::Get()->GetOVMaterialIndex())
          DSEventHandler::Get()->SetMuDepEnergy(DSEventHandler::Get()->GetMuDepEnergy() + theStep->GetTotalEnergyDeposit() / keV);

    if (DSStorage::Get()->GetWriteDeposits()) {
      // protection against writing excited atoms 
      bool isExcitedAtom = theStep->GetTrack()->GetDefinition()->GetPDGEncoding() > 1e9 && ((G4Ions*)theStep->GetTrack()->GetDefinition())->GetExcitationEnergy() > 0;
 
      if (not isExcitedAtom) { 
        // check copy number
        DSEventHandler::Get()->SetDepPID(theStep->GetTrack()->GetDefinition()->GetPDGEncoding());
        DSEventHandler::Get()->SetDepVolume(theStep->GetTrack()->GetVolume()->GetLogicalVolume()->GetMaterial()->GetIndex());
        DSEventHandler::Get()->SetDepEnergy(theStep->GetTotalEnergyDeposit() / keV);

        DSEventHandler::Get()->SetDepStep(theStep->GetStepLength() / um);
        DSEventHandler::Get()->SetDepTime(theStep->GetTrack()->GetGlobalTime() / ns);
        DSEventHandler::Get()->SetDepPosition(theStep->GetPostStepPoint()->GetPosition() / cm);

        DSEventHandler::Get()->SetDeposits();
      }
    }
  }
}


/*
 * $Log: DSSteppingAction.cc,v $
 * Revision 1.8  2015/10/31 12:18:57  pagnes
 * bug fixed
 *
 * Revision 1.7  2015/10/31 11:25:10  pagnes
 * fixed error with step info cout. Now working with ds/manager/log development
 *
 * Revision 1.6  2015/10/15 09:24:37  dfranco
 * all units updated: position in cm and energy in keV
 *
 * Revision 1.5  2015/04/28 11:46:43  pagnes
 * tpcene (sum of the deposits in ActiveLAr volume) un-commented
 *
 * Revision 1.4  2014/11/21 10:19:00  dfranco
 * added a command to scale the veto scintillation yield factor and fixed the
 * visible energy variable in the veto
 *
 * Revision 1.3  2014/11/13 16:47:05  dfranco
 * removed variables which were creating conflicts with the previous version of
 * g4ds10
 *
 * Revision 1.2  2014/10/13 18:43:57  swesterd
 * fixed veto PMT positions
 *
 * Revision 1.1  2014/05/07 12:21:05  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.13  2013/08/27 07:13:14  swesterd
 * add visible energy for the neutron veto
 *
 * Revision 1.12  2013/08/06 13:58:20  dfranco
 * Added 3 variables to store the deposited energy in LAr, scintillator, and
 * water. The last two are not yet implemented. g4rooter has been updated with 3
 * new variables: tpcene, vetoene, and muene
 *
 * Revision 1.11  2013/07/24 09:49:02  dfranco
 * Added S1 and S2 equivalent energies, the material index for each deposit, the
 * command killS1S2 to kill photons and electrons generated by DSLight (after
 * storing the equivalent energies)
 *
 * Revision 1.10  2013/07/07 09:52:39  dfranco
 * removed a cout
 *
 * Revision 1.9  2013/05/31 13:43:51  dfranco
 * Change the DSLog(development) to add pre and post step volumes
 *
 * Revision 1.8  2013/05/29 10:55:58  dfranco
 * Change DSLog level from debugging to development
 *
 * Revision 1.7  2013/05/27 11:28:20  dfranco
 * Segmentation fault bug fixed, for the TPC case only. Useless files removed.
 * Fixed the RINDEX, which was not defined in the whole range, for the
 * scintillator.
 *
 *
 */
