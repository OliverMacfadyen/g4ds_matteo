#include "DSGeneratorAmCSource.hh"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include "DSEventHandler.hh"
#include "DSGeneratorAmCSourceMessenger.hh"
#include "DSLogger.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SPSRandomGenerator.hh"
#include "G4SystemOfUnits.hh"
// #include "TFile.h"
// #include "TString.h"
using namespace std;

DSGeneratorAmCSource::DSGeneratorAmCSource() : DSVGenerator("DSGeneratorAmCSource") {

  isFirstTime = true;
  fPosition = G4ThreeVector(0., 0., 0.);
  fDirection = G4ThreeVector(0., 0., 0.);
  fSourceRotation = G4ThreeVector(0., 0., 0.);
  fParticleTable = G4ParticleTable::GetParticleTable();
  fTheMessenger = new DSGeneratorAmCSourceMessenger(this);
  DSLog(routine) << "AmC Source Generator Built" << endlog;
}

DSGeneratorAmCSource::~DSGeneratorAmCSource() {
  delete fTheMessenger;
}

void DSGeneratorAmCSource::DSGeneratePrimaries(G4Event* event) {
  if (isFirstTime) { isFirstTime = false; }

  fPosition = GetVParticlePosition();
  fParticle = fParticleTable->FindParticle(2112);
  G4int pdg = fParticle->GetPDGEncoding();
  G4double mass = fParticle->GetPDGMass();
  G4double KinE = 0;
  G4double theta = 0;

  // fenergyanglecorrelation2d->GetRandom2(theta, KinE);
  G4double sample_z_cumulative = G4UniformRand() ;
  G4int idx = 0 ;
  while (sample_z_cumulative > fenergyanglecorrelation2d_z[idx]) ++idx ;
  KinE  = fenergyanglecorrelation2d_y[idx];
  theta = fenergyanglecorrelation2d_x[idx];

  KinE *= MeV;
  theta *= radian;

  G4double energy = KinE + mass;
  G4double phi = twopi * G4UniformRand();  // phi uniform in [0, 2*pi]
  G4ThreeVector ur(std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta));
  G4ThreeVector dir = GetAmCSourceRotation();
  ur.rotateY(dir.theta());
  ur.rotateZ(dir.phi());
  fDirection = ur;

  G4double pmom = std::sqrt(energy * energy - mass * mass);
  G4double px = pmom * fDirection.x();
  G4double py = pmom * fDirection.y();
  G4double pz = pmom * fDirection.z();

  DSEventHandler::Get()->SetEnergy(KinE / keV);
  DSEventHandler::Get()->SetPosition(fPosition / cm);
  DSEventHandler::Get()->SetDirection(fDirection);
  DSEventHandler::Get()->SetPDG(pdg);
  DSEventHandler::Get()->SetTime(0);

  G4PrimaryVertex* vertex = new G4PrimaryVertex(fPosition, 0);
  G4PrimaryParticle* particle = new G4PrimaryParticle(fParticle, px, py, pz);

  vertex->SetPrimary(particle);
  event->AddPrimaryVertex(vertex);

  DSEventHandler::Get()->SetDId(1);
  DSEventHandler::Get()->SetDPDG(pdg);
  DSEventHandler::Get()->SetDTime(0.);
  DSEventHandler::Get()->SetDPosition(fPosition / cm);
  DSEventHandler::Get()->SetDEnergy(KinE / keV);
  DSEventHandler::Get()->SetDDirection(fDirection);
  DSEventHandler::Get()->SetDaughters();
}

/*
void DSGeneratorAmCSource::SetAmCEnergyAngleCorrelationFileName(string filename) {
  //  string filename="../data/physics/amcNoFoil_out.root";
  TFile* fFile = new TFile(filename.c_str());
  DSLog(routine) << "data file for AmC source:  " << filename << endlog;
  if (fFile->IsZombie()) {
    DSLog(error) << "data file for AmC source could not be opened." << endlog;
    DSLog(fatal) << endlog;
  }
  fenergyanglecorrelation2d = (TH2F*)fFile->Get("thetaNVsEnH");
}
*/

void DSGeneratorAmCSource::SetAmCEnergyAngleCorrelationFileName(string filename) {
  ifstream fin (filename.c_str() ) ;
  int idx = 0 ;
  while (1)  {
    if (idx >= 8000) break  ;
    fin >> fenergyanglecorrelation2d_x[idx] >> fenergyanglecorrelation2d_y[idx] >> fenergyanglecorrelation2d_z[idx] ;
    ++idx ;
  }

}
