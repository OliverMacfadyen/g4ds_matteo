#include "DSTrackingAction.hh"
#include <vector>
#include "DSEventHandler.hh"
#include "DSLogger.hh"
#include "DSParameters.hh"
#include "DSStorage.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"
#include "G4TrackStatus.hh"
#include "G4TrackingManager.hh"

using namespace std;

DSTrackingAction::DSTrackingAction() {
  ;
}

void DSTrackingAction::PreUserTrackingAction(const G4Track* aTrack) {
  /*
  if(aTrack->GetDefinition()->GetPDGEncoding() != -22)
  cout << "Track " << aTrack->GetDefinition()->GetParticleName() << " "
       << aTrack->GetVertexKineticEnergy ()/keV << " "
   << aTrack->GetTrackID() << " "
   << aTrack->GetParentID() << " "
       <<endl ;


  if(   aTrack->GetDefinition()->GetPDGEncoding() == 5
     || aTrack->GetDefinition()->GetPDGEncoding() == 6
     || aTrack->GetDefinition()->GetPDGEncoding() == 27
     || aTrack->GetDefinition()->GetPDGEncoding() == 28
     || aTrack->GetDefinition()->GetPDGEncoding() == 43
     || aTrack->GetDefinition()->GetPDGEncoding() == 44
  ) aTrack->SetTrackStatus(fStopAndKill );



  //cout << fKillTrackAndSecondaries << endl ;
  //if(aTrack->GetGlobalTime () > 1.*ms) aTrack->SetTrackStatus(  3);

 */
  if (DSStorage::Get()->GetWriteDaughters() &&
      DSStorage::Get()->GetRDMDecay() == 0
      //      && aTrack->GetParentID()  <= DSStorage::Get()->GetNDaughters() //
      //      default = 1
      && aTrack->GetTrackID() > 1 && DSStorage::Get()->GetIsEnDepGenerator() == false && aTrack->GetDefinition()->GetPDGEncoding() != -22) {

    int index = aTrack->GetCreatorProcess()->GetProcessType() * 1000 + aTrack->GetCreatorProcess()->GetProcessSubType();
    DSEventHandler::Get()->SetDId(aTrack->GetTrackID());    // int(DSEventHandler::Get()->GetVDaughters().size()));
    DSEventHandler::Get()->SetDPID(aTrack->GetParentID());  // aTrack->GetTrackID());
    DSEventHandler::Get()->SetDPDG(aTrack->GetDefinition()->GetPDGEncoding());
    DSEventHandler::Get()->SetDProcess(index);
    DSEventHandler::Get()->SetDTime(aTrack->GetGlobalTime() / ns);
    DSEventHandler::Get()->SetDEnergy(aTrack->GetVertexKineticEnergy() / keV);
    DSEventHandler::Get()->SetDPosition(aTrack->GetVertexPosition() / m);
    DSEventHandler::Get()->SetDDirection(aTrack->GetVertexMomentumDirection());
    DSEventHandler::Get()->SetDaughters();

    if (aTrack->GetDefinition()->GetPDGEncoding() != -22) DSLog(development) << " Save daughters: " << aTrack->GetDefinition()->GetParticleName() << " " << aTrack->GetTrackID() << " " << aTrack->GetParentID() << " " << aTrack->GetVertexKineticEnergy() / keV << " " << aTrack->GetGlobalTime() / ns << endlog;
  }
}

/////////////////////////////////////////////////////////////////////////

void DSTrackingAction::PostUserTrackingAction(const G4Track* aTrack) {

  int nchannel = -10;

  if (aTrack->GetTrackStatus() == fStopAndKill && aTrack->GetDefinition()->GetPDGEncoding() == -22 && DSLogger::GetSeverity() == DSLogger::development)
    DSLog(development) << " EoT " << aTrack->GetDefinition()->GetParticleName() << " " << aTrack->GetStep()->GetPostStepPoint()->GetPosition() / cm << " " << aTrack->GetVolume()->GetName() << endlog;

  // save tpc photoelectrons
  if (DSStorage::Get()->GetPMTMaterialIndex() == (int)aTrack->GetVolume()->GetLogicalVolume()->GetMaterial()->GetIndex() && aTrack->GetDefinition()->GetPDGEncoding() == -22 && aTrack->GetTrackStatus() == fStopAndKill) {

    if (DSStorage::Get()->Get20KGeometry()) {
      // Fill PE structure
      nchannel = -10;
      // select ph reaching the SiPM planes in the TPC (the copy No is manually
      // set to -56, the same as the material index (x -1)
      if (aTrack->GetVolume()->GetCopyNo() == -56) nchannel = GetDS20kNChannel((G4ThreeVector)aTrack->GetStep()->GetPostStepPoint()->GetPosition());
      else
        nchannel = -1;  // Veto SiPMs

      // compute binomial PDE selection only once
      G4double detect = G4UniformRand();

      // TPC PE's
      if (detect < DSParameters::Get()->GetSiPMPDE(nm * h_Planck * c_light * 1E12 / (aTrack->GetKineticEnergy() / eV)) /*PDE*/
          && nchannel > -1) {
        DSEventHandler::Get()->SetPhotoElectronPMT(nchannel);
        DSEventHandler::Get()->SetPhotoElectronTime(aTrack->GetGlobalTime() / ns);
        DSEventHandler::Get()->SetPhotoElectrons();
      }

      // Veto PE's
      if (nchannel == -1) {
        if (detect < DSParameters::Get()->GetSiPMPDE(nm * h_Planck * c_light * 1E12 / (aTrack->GetKineticEnergy() / eV))) {

          int myCopyNumber = 0;
          //PA 10/23 - tile index + 4*quadrant_index + PDU index
          myCopyNumber = aTrack->GetTouchableHandle()->GetCopyNumber(2) + 4 * aTrack->GetTouchableHandle()->GetCopyNumber(1) + aTrack->GetTouchableHandle()->GetCopyNumber(0);

          DSEventHandler::Get()->SetVetoPhotoElectronPMT(myCopyNumber);
          DSEventHandler::Get()->SetVetoPhotoElectronTime(aTrack->GetGlobalTime() / ns);
          DSEventHandler::Get()->SetVetoPhotoElectrons();
        }
      }
      /*  */
      /*
  //for debugging purposees, store everything at this stage
  DSEventHandler::Get()->SetPhotonVolumeID(my_sum_cpy_n);
  //DSEventHandler::Get()->SetPhotonVolumeID(aTrack->GetVolume()->GetCopyNo ()
  ); DSEventHandler::Get()->SetPhotonTime(aTrack->GetGlobalTime()/ns);
  DSEventHandler::Get()->SetPhotonPID(detect>0.4?0:1);
  DSEventHandler::Get()->SetPhotonWavelength(h_Planck*c_light*1E12/(aTrack->GetKineticEnergy()/eV));
  DSEventHandler::Get()->SetPhotonPosition(aTrack->GetStep()->GetPostStepPoint()->GetPosition()/cm);
  DSEventHandler::Get()->SetPhotons();
  */

    } else if (DSStorage::Get()->GetArDMGeometry()) {
      DSEventHandler::Get()->SetPhotoElectronTime(aTrack->GetGlobalTime() / ns);
      DSEventHandler::Get()->SetPhotoElectrons();
    } else if (DSStorage::Get()->GetReDGeometry()) {
      DSEventHandler::Get()->SetPhotonVolumeID(aTrack->GetVolume()->GetLogicalVolume()->GetMaterial()->GetIndex());
      // DSEventHandler::Get()->SetPhotonVolumeID(aTrack->GetVolume()->GetCopyNo
      // ()  );
      DSEventHandler::Get()->SetPhotonTime(aTrack->GetGlobalTime() / ns);
      DSEventHandler::Get()->SetPhotonPID(aTrack->GetParentID());
      DSEventHandler::Get()->SetPhotonWavelength(h_Planck * c_light * 1E12 / (aTrack->GetKineticEnergy() / eV));
      DSEventHandler::Get()->SetPhotonPosition(aTrack->GetStep()->GetPostStepPoint()->GetPosition() / cm);
      DSEventHandler::Get()->SetPhotons();
    } else if (DSStorage::Get()->GetPETGeometry()) {

      if (aTrack->GetKineticEnergy() < 5 * eV && G4UniformRand() < 0.4) {
        DSEventHandler::Get()->SetPhotonTime(aTrack->GetGlobalTime() / ns);
        DSEventHandler::Get()->SetPhotonVolumeID(aTrack->GetVolume()->GetLogicalVolume()->GetMaterial()->GetIndex());
        DSEventHandler::Get()->SetPhotonPID(aTrack->GetParentID());
        DSEventHandler::Get()->SetPhotonWavelength(h_Planck * c_light * 1E12 / (aTrack->GetKineticEnergy() / eV));
        DSEventHandler::Get()->SetPhotonPosition(aTrack->GetStep()->GetPostStepPoint()->GetPosition() / cm);
        DSEventHandler::Get()->SetPhotons();
      }

    } else {
      G4String pmtname = aTrack->GetVolume()->GetName();
      if (!pmtname.find("TPMT_")) {
        // apply the QE cut
        if (G4UniformRand() < DSParameters::Get()->GetTPCQE(h_Planck * c_light * 1E12 / (aTrack->GetKineticEnergy() / eV))) {
          pmtname.erase(0, 5);
          DSEventHandler::Get()->SetPhotoElectronPMT(atoi(pmtname.c_str()));
          DSEventHandler::Get()->SetPhotoElectronTime(aTrack->GetGlobalTime() / ns);
          DSEventHandler::Get()->SetPhotoElectrons();
          if (DSStorage::Get()->GetVerbosity() > 3) {
            DSLog(development) << "TPC PhotoElectron: "
                               << "PMT " << DSEventHandler::Get()->GetPhotoElectronPMT() << "T   " << DSEventHandler::Get()->GetPhotoElectronTime() << " ns " << aTrack->GetVolume()->GetLogicalVolume()->GetMaterial()->GetName() << " " << aTrack->GetVolume()->GetName() << " " << aTrack->GetVolume()->GetCopyNo() << " "
                               << aTrack->GetVolume()->GetLogicalVolume()->GetName() << " " << aTrack->GetStep()->GetPostStepPoint()->GetPosition() / m << " " << aTrack->GetStep()->GetPostStepPoint()->GetPosition().mag() / m << " " << endlog;
          }
        }
      }
    }
  }

  // save veto photoelectrons
  if (DSStorage::Get()->GetVetoPMTMaterialIndex() == (int)aTrack->GetVolume()->GetLogicalVolume()->GetMaterial()->GetIndex() && aTrack->GetDefinition()->GetPDGEncoding() == -22 && aTrack->GetTrackStatus() == fStopAndKill) {
    DSLog(development) << "in veto PMT" << endlog;
    G4String pmtname = aTrack->GetVolume()->GetName();
    if (!pmtname.find("VPMT_")) {

      // Compute the photon wavelength
      const G4double myPhotonWL = nm * h_Planck * c_light * 1E12 / (aTrack->GetKineticEnergy() / eV);
      // apply the QE cut
      if (G4UniformRand() < DSParameters::Get()->GetNVQE(myPhotonWL)) {
        pmtname.erase(0, 5);
        DSEventHandler::Get()->SetVetoPhotoElectronPMT(atoi(pmtname.c_str()));
        DSEventHandler::Get()->SetVetoPhotoElectronTime(aTrack->GetGlobalTime() / ns);
        DSEventHandler::Get()->SetVetoPhotoElectrons();

        if (DSStorage::Get()->GetVerbosity() > 3) {
          DSLog(development) << "Veto PhotoElectron: "
                             << "PMT " << DSEventHandler::Get()->GetVetoPhotoElectronPMT() << "T   " << DSEventHandler::Get()->GetVetoPhotoElectronTime() << " ns " << aTrack->GetVolume()->GetLogicalVolume()->GetMaterial()->GetName() << " " << aTrack->GetVolume()->GetName() << " " << aTrack->GetVolume()->GetCopyNo()
                             << " " << aTrack->GetVolume()->GetLogicalVolume()->GetName() << " " << aTrack->GetStep()->GetPostStepPoint()->GetPosition() / m << " " << aTrack->GetStep()->GetPostStepPoint()->GetPosition().mag() / m << " " << endlog;
        }
      }
    }
  }

  // save water tank photoelectrons
  if (DSStorage::Get()->GetMuPMTMaterialIndex() == (int)aTrack->GetVolume()->GetLogicalVolume()->GetMaterial()->GetIndex() && aTrack->GetDefinition()->GetPDGEncoding() == -22 && aTrack->GetTrackStatus() == fStopAndKill) {
    G4String pmtname = aTrack->GetVolume()->GetName();
    if (!pmtname.find("WPMT_")) {
      pmtname.erase(0, 5);
      DSEventHandler::Get()->SetMuPhotoElectronPMT(atoi(pmtname.c_str()));
      DSEventHandler::Get()->SetMuPhotoElectronTime(aTrack->GetGlobalTime() / ns);
      DSEventHandler::Get()->SetMuPhotoElectrons();

      if (DSStorage::Get()->GetVerbosity() > 3) {
        DSLog(development) << "WT PhotoElectron: "
                           << "PMT " << DSEventHandler::Get()->GetMuPhotoElectronPMT() << "T   " << DSEventHandler::Get()->GetMuPhotoElectronTime() << " ns " << aTrack->GetVolume()->GetLogicalVolume()->GetMaterial()->GetName() << " " << aTrack->GetVolume()->GetName() << " " << aTrack->GetVolume()->GetCopyNo() << " "
                           << aTrack->GetVolume()->GetLogicalVolume()->GetName() << " " << aTrack->GetStep()->GetPostStepPoint()->GetPosition() / m << " " << aTrack->GetStep()->GetPostStepPoint()->GetPosition().mag() / m << " " << endlog;
      }
    }
  }

  // save photons
  if (DSStorage::Get()->GetWritePhotons() && aTrack->GetTrackStatus() == fStopAndKill && aTrack->GetDefinition()->GetPDGEncoding() == -22
      /*   && !DSStorage::Get()->Get20KGeometry() */) {
    if (DSStorage::Get()->GetVerbosity() > 4) {
      DSLog(development) << "Photon: " << aTrack->GetVolume()->GetLogicalVolume()->GetMaterial()->GetName() << " " << aTrack->GetVolume()->GetLogicalVolume()->GetName() << " " << aTrack->GetStep()->GetPostStepPoint()->GetPosition() / m << " " << aTrack->GetStep()->GetPostStepPoint()->GetPosition().mag() / m << " "
                         << endlog;
    }
    DSEventHandler::Get()->SetPhotonVolumeID(aTrack->GetVolume()->GetLogicalVolume()->GetMaterial()->GetIndex());
    // DSEventHandler::Get()->SetPhotonVolumeID(aTrack->GetVolume()->GetCopyNo
    // ()  );
    DSEventHandler::Get()->SetPhotonTime(aTrack->GetGlobalTime() / ns);
    DSEventHandler::Get()->SetPhotonPID(aTrack->GetParentID());
    // DSEventHandler::Get()->SetPhotonPID(nchannel);
    DSEventHandler::Get()->SetPhotonWavelength(h_Planck * c_light * 1E12 / (aTrack->GetKineticEnergy() / eV));
    DSEventHandler::Get()->SetPhotonPosition(aTrack->GetStep()->GetPostStepPoint()->GetPosition() / cm);

    DSEventHandler::Get()->SetPhotons();
  }
}

int DSTrackingAction::GetDS20kNChannel(G4ThreeVector position) {

  // everything in cm

  G4double z = position.z() / cm;
  G4double x = position.x() / cm;
  G4double y = position.y() / cm;

  G4double MB_heigh = 5 * 5 /*5 PMDs*/ + 2 * 0.25 /*two times half the spacing */;
  G4double Maximal_coordinate = 7 * MB_heigh;

  if (DSStorage::Get()->GetDS20kTPCedge() < 0.6 * m) {
    x += MB_heigh / 2.;
    y += MB_heigh / 2.;
  }

  if (abs(x) > Maximal_coordinate) return -2;
  if (abs(y) > Maximal_coordinate) return -2;

  x += Maximal_coordinate;
  y += Maximal_coordinate;

  int nchannel = 0;
  int noffset = 0;

  int MB_x = 0;
  int MB_y = 0;

  while (x > MB_x * MB_heigh) ++MB_x;
  while (y > MB_y * MB_heigh) ++MB_y;

  G4double MB_bot_left_x = (MB_x - 1) * MB_heigh;
  G4double MB_bot_left_y = (MB_y - 1) * MB_heigh;

  nchannel = (MB_x - 1) * 25 + (MB_y - 1) * 25 * 14;

  G4double xtile = x - 1. * MB_bot_left_x;
  G4double ytile = y - 1. * MB_bot_left_y;

  if (xtile < 0.25) return -3;  // cut 0.25 cm from edges
  if (ytile < 0.25) return -3;
  if (xtile > 5 * 5 + 0.25) return -4;  // cut 0.25 cm from edges
  if (ytile > 5 * 5 + 0.25) return -4;

  int tile_x = 0;
  int tile_y = 0;

  while (xtile > tile_x * 5 + 0.25) ++tile_x;
  while (ytile > tile_y * 5 + 0.25) ++tile_y;

  nchannel += ((tile_x - 1) * 5 + tile_y - 1);

  if (z > 0.) noffset = 196 * 25;
  nchannel += noffset;

  return nchannel;
}
