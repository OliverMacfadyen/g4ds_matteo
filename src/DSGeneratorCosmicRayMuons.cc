//---------------------------------------------------------------------------//
//    Adapted by davide.franco@mi.infn.it from the MaGe code:
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// bb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nu//
//                                                                           //
//                                                                           //
//                         MaGe Simulation                                   //
//                                                                           //
//      This code implementation is the intellectual property of the         //
//      MAJORANA and Gerda Collaborations. It is based on Geant4, an         //
//      intellectual property of the RD44 GEANT4 collaboration.              //
//                                                                           //
//                        *********************                              //
//                                                                           //
//    Neither the authors of this software system, nor their employing       //
//    institutes, nor the agencies providing financial support for this      //
//    work  make  any representation or  warranty, express or implied,       //
//    regarding this software system or assume any liability for its use.    //
//    By copying, distributing or modifying the Program (or any work based   //
//    on on the Program) you indicate your acceptance of this statement,     //
//    and all its terms.                                                     //
//                                                                           //
// bb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nu//
//---------------------------------------------------------------------------//

#include "DSGeneratorCosmicRayMuons.hh"
#include "DSLogger.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
//#include "DSManager.hh"
//#include "DSDataCollection.hh"

#include <math.h>
#include "DSEventHandler.hh"
#include "G4Event.hh"
#include "G4MuonMinus.hh"
#include "G4MuonPlus.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4PrimaryParticle.hh"
#include "G4TransportationManager.hh"
#include "G4VPhysicalVolume.hh"
#include "Randomize.hh"

DSGeneratorCosmicRayMuons::DSGeneratorCosmicRayMuons() : DSVGenerator("DSGeneratorCosmicRayMuons") {

  NumberOfParticlesToBeGenerated = 1;
  particle_definition = G4MuonMinus::MuonMinusDefinition();
  G4ThreeVector zero(0., 0., 0.);
  particle_momentum_direction = G4ParticleMomentum(0., 0., -1.);
  particle_energy = 500 * GeV;
  particle_position = zero;
  particle_time = 0.0;
  particle_polarization = zero;
  particle_charge = particle_definition->GetPDGCharge();

  // geometry initialization
  halfz = 8.1 * m;
  Radius = 8.0 * m;
  spectralIndex = 3.7;
  rockDepth = 3700.0 * m;  // m.w.e.

  fileName = "zenith_azimuth.dat";
  samplerInitialized = false;
  fileFound = false;

  // limits for the energy sampling
  energysup = 10 * TeV;
  energyinf = 1 * GeV;

  fTheMessenger = new DSGeneratorCosmicRayMuonsMessenger(this);
  gNavigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
  // DSOutputVertex::Get()->SetRagMin(halfz);
  // DSOutputVertex::Get()->SetRagMax(Radius);

  angularpdf = new G4DataVector();
  costheta_pdf = new G4DataVector();
  phi_pdf = new G4DataVector;
  nx = 0;
  ny = 0;
}

DSGeneratorCosmicRayMuons::~DSGeneratorCosmicRayMuons() {
  delete fTheMessenger;
  delete angularpdf;
  delete costheta_pdf;
  delete phi_pdf;
}

void DSGeneratorCosmicRayMuons::SampleInitialPosition() {
  G4ThreeVector RandPos;
  G4double z = 0.;

  // Positions are generated uniformly from a circle
  G4double r = sqrt(Radius * Radius * G4UniformRand());
  G4double phi = twopi * G4UniformRand();
  RandPos.setX(r * cos(phi));
  RandPos.setY(r * sin(phi));
  RandPos.setZ(z);
  G4ThreeVector CentreCoords(0., 0., halfz);
  particle_position = CentreCoords + RandPos;
}

void DSGeneratorCosmicRayMuons::GenerateAngularSpectrum() {
  G4double costheta = 0.;
  G4double phi = 0;
  G4double sintheta = 0;
  size_t j = 0;

  if (!fileFound) {
    // According to Phys. Rev. D 44 (1991) 3543, the angular spectrum
    // underground is proportional to 1/cos(theta) --> with respect to the
    // horizontal
    G4double costhetamin = 0.00001;
    costheta = pow(costhetamin, (1.0 - G4UniformRand()));
    phi = G4UniformRand() * 360 * degree;
    sintheta = sqrt(1.0 - costheta * costheta);
    particle_momentum_direction.setX(cos(phi) * costheta);
    particle_momentum_direction.setY(sin(phi) * costheta);
    particle_momentum_direction.setZ(-1.0 * sintheta);
  } else {  // sampling from file
    // warning: the spectrum MUST be evenly-spaced in costheta and in phi
    G4double step1 = 1. / ((G4double)nx);
    G4double step2 = 360. / ((G4double)ny);

    // explicit sampling of costheta and phi
    G4double xrand = G4UniformRand();
    // this is slow. Something more efficienct would be good
    for (j = 0; j < costheta_pdf->size(); j++) {
      if (xrand > (*angularpdf)[j] && xrand <= (*angularpdf)[j + 1]) break;
    }
    costheta = step1 * (G4UniformRand() - 0.5) + (*costheta_pdf)[j];
    phi = (step2 * (G4UniformRand() - 0.5) + (*phi_pdf)[j]) * degree;
    if (costheta < 0) costheta = 0;
    if (costheta > 1) costheta = 1;
    sintheta = sqrt(1.0 - costheta * costheta);

    // DSLog(debugging) << xrand << " " << j << " " << costheta << " " <<
    // phi/degree << endlog;

    particle_momentum_direction.setX(cos(phi) * sintheta);
    particle_momentum_direction.setY(sin(phi) * sintheta);
    particle_momentum_direction.setZ(-1.0 * costheta);
  }
}

void DSGeneratorCosmicRayMuons::ShootEnergy() {
  G4double step = (log10(energysup) - log10(energyinf)) / ((G4double)nbin - 1);
  G4int j = 0;

  // initialize the sampler: logarithmic sampling between 1 GeV and 10 TeV
  // angular distribution is read from file (MACRO) Astrop. Journ. 412: 310-311
  // 1993
  if (!samplerInitialized) {
    G4double log_energy = 0;
    G4double true_energy = 0;
    userpdf[0] = 0.;
    for (j = 1; j < nbin; j++) {
      log_energy = j * step;
      true_energy = energyinf * std::pow(10, log_energy);
      G4double base = energyinf * std::pow(10, log_energy) * (1.0 - std::pow(10, -1.0 * step));
      userpdf[j] = userpdf[j - 1] + MuonSpectrum(true_energy) * base;
    }
    G4double xnormalization = userpdf[nbin - 1];
    for (j = 0; j < nbin; j++) { userpdf[j] = userpdf[j] / xnormalization; }
    DSLog(trace) << "Energy sampler initialized." << endlog;
    DSLog(trace) << "Spectrum from Lipari and Stanev, Phys. Rev. D 44 (1991) 3543" << endlog;

    // Read from file the MACRO angular spectrum
    size_t jj = 0;
    G4String pathString;
    char* path = getenv("DSDATA");
    if (path) pathString = path;
    else  // if (!path)
    {
      DSLog(warning) << "DADATA environment variable not set!" << endlog;
      pathString = "../data/physics";
    }
    G4String pathFile = pathString + "/" + fileName;
    std::ifstream file(pathFile);
    if (file.is_open()) {
      G4double bb = 0;
      G4double aa = 0;
      G4double cc = 0;
      angularpdf->push_back(bb);

      do {
        aa = 0;
        bb = 0;
        cc = 0;
        file >> aa >> cc >> bb;  // format: costheta, phi, distribution
        //	DSLog(debugging) << aa << " " << bb << endlog;
        if (aa != 0 || bb != 0 || cc != 0)  // prevents the last value to be inserted twice
        {
          costheta_pdf->push_back(aa);
          phi_pdf->push_back(cc);
          bb += (*angularpdf)[angularpdf->size() - 1];
          angularpdf->push_back(bb);
        }
      } while (!file.eof());

      xnormalization = (*angularpdf)[angularpdf->size() - 1];
      for (jj = 0; jj < angularpdf->size(); jj++) {
        (*angularpdf)[jj] = (*angularpdf)[jj] / xnormalization;
        // DSLog(debugging) << "angularpdf[" << jj << "] = " <<
        // (*angularpdf)[jj] << endlog;
      }
      // Finds the number of bins in x (costheta) and y (phi)
      nx = 0;
      ny = 0;
      for (jj = 0; jj < costheta_pdf->size(); jj++) {
        if ((*costheta_pdf)[0] == (*costheta_pdf)[jj]) ny++;
        if ((*phi_pdf)[0] == (*phi_pdf)[jj]) nx++;
      }
      DSLog(trace) << "Angular sampler initialized. " << endlog;
      DSLog(trace) << "Data from MACRO, Astrop. Journ. 412 (1993) 310-311" << endlog;
      DSLog(debugging) << "Total binning: " << nx << " in costheta" << endlog;
      DSLog(debugging) << "Total binning: " << ny << " in phi" << endlog;
      file.close();
      fileFound = true;
    } else {
      DSLog(warning) << "Data file " + pathFile + " not found!" << endlog;
      DSLog(warning) << "Using 1/(costheta) angolar distribution " << endlog;
    }

    samplerInitialized = true;
  }
  // The spectrum generated in this way is coincident with the one reported in
  // Fig.7 of Astrop. Phys. 6 (1997) 129 for the Gran Sasso depth

  G4double xrand = G4UniformRand();
  for (j = 0; j < (nbin - 1); j++) {
    if (xrand > userpdf[j] && xrand <= userpdf[j + 1]) break;
  }
  G4double sampled_bin = (G4double)j;
  if (userpdf[j + 1] != userpdf[j]) sampled_bin += ((userpdf[j + 1] - xrand) / (userpdf[j + 1] - userpdf[j]));
  particle_energy = energyinf * pow(10, sampled_bin * step);
  return;
}

void DSGeneratorCosmicRayMuons::SetParticleDefinition(G4ParticleDefinition* aParticleDefinition) {
  particle_definition = aParticleDefinition;
  particle_charge = particle_definition->GetPDGCharge();
}

void DSGeneratorCosmicRayMuons::DSGeneratePrimaries(G4Event* event) {

  // DSLog(debugging) << "Start generation of primary vertex" << endlog;

  if (particle_definition == NULL) {
    DSLog(error) << "No particle has been defined!" << endlog;
    return;
  }

  SampleInitialPosition();

  // Energy stuff
  ShootEnergy();

  // Angular stuff
  GenerateAngularSpectrum();

  // create a new vertex
  G4PrimaryVertex* vertex = new G4PrimaryVertex(particle_position, particle_time);

  // create new primaries and set them to the vertex
  G4double mass = particle_definition->GetPDGMass();
  G4double energy = particle_energy + mass;
  G4double pmom = sqrt(energy * energy - mass * mass);
  G4double px = pmom * particle_momentum_direction.x();
  G4double py = pmom * particle_momentum_direction.y();
  G4double pz = pmom * particle_momentum_direction.z();

  if (G4UniformRand() > 0.5) particle_definition = G4MuonMinus::MuonMinusDefinition();
  else
    particle_definition = G4MuonPlus::MuonPlusDefinition();

  DSLog(debugging) << "Particle name: " << particle_definition->GetParticleName() << endlog;
  DSLog(debugging) << "       Energy: " << particle_energy / GeV << " GeV" << endlog;
  DSLog(debugging) << "     Position: " << particle_position / m << " m" << endlog;
  DSLog(debugging) << "    Direction: " << particle_momentum_direction << endlog;
  //   DSLog(debugging) << " NumberOfParticlesToBeGenerated: "
  // 		   << NumberOfParticlesToBeGenerated << endlog;

  DSEventHandler::Get()->SetPosition(particle_position / cm);
  DSEventHandler::Get()->SetEnergy(particle_energy / keV);
  DSEventHandler::Get()->SetDirection(particle_momentum_direction);
  DSEventHandler::Get()->SetPDG(particle_definition->GetPDGEncoding());

  for (G4int i = 0; i < NumberOfParticlesToBeGenerated; i++) {
    G4PrimaryParticle* particle = new G4PrimaryParticle(particle_definition, px, py, pz);
    particle->SetMass(mass);
    particle->SetCharge(particle_charge);
    particle->SetPolarization(particle_polarization.x(), particle_polarization.y(), particle_polarization.z());
    vertex->SetPrimary(particle);
  }

  event->AddPrimaryVertex(vertex);

  G4String materialName = gNavigator->LocateGlobalPointAndSetup(particle_position)->GetName();
  DSLog(debugging) << "Created primary vertex in " << materialName << endlog;

  // theDataCollection->SetEnergy(particle_energy);
  // theDataCollection->SetPosition(particle_position);
  // theDataCollection->SetDirection(particle_momentum_direction);
  // theDataCollection->SetShootDist(sqrt(particle_momentum_direction.x()*particle_momentum_direction.x()+
  //                                    particle_momentum_direction.y()*particle_momentum_direction.y()));
}

G4double DSGeneratorCosmicRayMuons::MuonSpectrum(G4double ene) {
  // The parametrization is taken from P. Lipari and T. Stanev
  // Phys. Rev. D 44 (1991) 3543

  // The energy must be expressed in TeV
  // Depth in km w.e.

  G4double beta, epsilon;
  if (spectralIndex == 2.0) {
    beta = 0.465;
    epsilon = 0.569;
  } else if (spectralIndex == 2.7) {
    beta = 0.418;
    epsilon = 0.557;
  } else  // this is the default case
  {
    beta = 0.383;
    epsilon = 0.618;
  }

  G4double part1 = epsilon * (1 - std::exp(-1.0 * beta * rockDepth / km));
  G4double part2 = std::pow((ene / TeV) + part1, -1 * spectralIndex);
  G4double part3 = std::exp(beta * (rockDepth / km) * (1 - spectralIndex));
  return part2 * part3;
}

/*
 * Revision 1.2  2013/03/20 09:54:28  dfranco
 * New version of g4ds
 *

 */
