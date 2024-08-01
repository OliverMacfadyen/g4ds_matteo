//---------------------------------------------------------------------------//
//    Code by pablo.garcia@ciemat.es to read files in HEPEVT format:
//---------------------------------------------------------------------------//

#include "DSGeneratorHEPevt.hh"
#include "DSLogger.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include <math.h>
#include "DSEventHandler.hh"
#include "G4Event.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4PrimaryParticle.hh"
#include "G4RunManager.hh"
#include "G4TransportationManager.hh"
#include "G4VPhysicalVolume.hh"

DSGeneratorHEPevt::DSGeneratorHEPevt() : DSVGenerator("DSGeneratorHEPevt") {

  // initialization
  fileName = "events.hepevt";
  firstEvent = true;

  fParticleTable = G4ParticleTable::GetParticleTable();

  fTheMessenger = new DSGeneratorHEPevtMessenger(this);
  // gNavigator =
  // G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();

  DSLog(debugging) << "HEPevt Generator Built" << endlog;
}

DSGeneratorHEPevt::~DSGeneratorHEPevt() {
  HEPfile.close();
  delete fTheMessenger;
}

/***
void DSGeneratorHEPevt::SetParticleDefinition (G4ParticleDefinition*
aParticleDefinition) { particle_definition = aParticleDefinition;
  // particle_charge = particle_definition->GetPDGCharge();
}
***/

void DSGeneratorHEPevt::DSGeneratePrimaries(G4Event* event) {

  if (event) {}  //////////////// temporary hack to avoid warnings ////////////////////

  DSLog(debugging) << "Inside HEPevt generator !!" << endlog;
  DSLog(debugging) << "HEPevt filename: " << fileName << endlog;

  int dummy;

  if (firstEvent) {
    firstEvent = false;
    nevhep = -1;

    HEPfile.open(fileName);

    if (!HEPfile.is_open()) {
      DSLog(error) << "Unable to open input file!" << endlog;
      G4RunManager::GetRunManager()->AbortRun(true);
      return;
    }
  }

  // read one event each time
  HEPfile >> dummy;

  if (!HEPfile) {
    G4RunManager::GetRunManager()->AbortRun(true);
    return;
  }

  nevhep++;

  HEPfile >> nhep;

  DSLog(debugging) << "Event number " << nevhep << " - particles " << nhep << endlog;

  G4int particle_counter = 0;

  G4ThreeVector position = GetVParticlePosition();

  for (int i = 0; i < nhep; i++) {
    HEPfile >> isthep >> idhep >> jmohep1 >> jmohep2 >> jdahep1 >> jdahep2;
    for (int j = 0; j < 5; j++) HEPfile >> phep[j];
    for (int j = 0; j < 4; j++) HEPfile >> vhep[j];

    DSLog(debugging) << "Particle # " << i + 1 << " - status " << isthep << " - PDG " << idhep << endlog;

    // do the following for stable particles (istat==1) other than neutrinos
    // (12, 14, 16)
    if (isthep == 1 && abs(idhep) != 12 && abs(idhep) != 14 && abs(idhep) != 16) {
      DSLog(debugging) << "px     " << phep[0] << " GeV" << endlog;
      DSLog(debugging) << "py     " << phep[1] << " GeV" << endlog;
      DSLog(debugging) << "pz     " << phep[2] << " GeV" << endlog;
      DSLog(debugging) << "energy " << phep[3] << " GeV" << endlog;
      DSLog(debugging) << "mass   " << phep[4] << " GeV" << endlog;

      // position from event generator
      //    particle_position.setX(vhep[0]*mm);
      //    particle_position.setY(vhep[1]*mm);
      //    particle_position.setZ(vhep[2]*mm);

      // position from mac file definition
      particle_position = position;
      DSLog(debugging) << "position " << position << " mm" << endlog;

      particle_time = vhep[3] / c_light;  // units: ns
      DSLog(debugging) << "time " << particle_time << " ns" << endlog;

      particle_momentum.setX(phep[0] * GeV);
      particle_momentum.setY(phep[1] * GeV);
      particle_momentum.setZ(phep[2] * GeV);

      particle_energy = phep[3] * GeV;

      fParticle = fParticleTable->FindParticle(idhep);

      G4PrimaryVertex* vertex = new G4PrimaryVertex(particle_position, particle_time);
      G4PrimaryParticle* particle = new G4PrimaryParticle(fParticle, particle_momentum.getX(), particle_momentum.getY(), particle_momentum.getZ());

      vertex->SetPrimary(particle);
      event->AddPrimaryVertex(vertex);

      DSEventHandler::Get()->SetDId(particle_counter);
      DSEventHandler::Get()->SetDPDG(idhep);
      DSEventHandler::Get()->SetDPosition(particle_position);
      DSEventHandler::Get()->SetDDirection(particle_momentum / keV);
      DSEventHandler::Get()->SetDEnergy(particle_energy / keV);
      DSEventHandler::Get()->SetDTime(particle_time / ns);
      DSEventHandler::Get()->SetDaughters();

      particle_counter++;
    }
  }
}
