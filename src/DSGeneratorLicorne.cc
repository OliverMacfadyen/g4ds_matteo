#include <fstream>
#include "G4PhysicalConstants.hh"
#include "G4Poisson.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include "G4ParticleGun.hh"

#include "DSEventHandler.hh"
#include "DSGeneratorLicorne.hh"

DSGeneratorLicorne::DSGeneratorLicorne() : DSVGenerator("DSGeneratorLicorne") {

  fPulseMode = false;
  fIsFirstEvent = true;
  fParticleTable = G4ParticleTable::GetParticleTable();

  fDirection = G4ThreeVector(0., 0., 1.);
  fPosition = G4ThreeVector(0., 0., 0.);

  fSamples = 0;
  fGammaToNeutronRatio = 0.;
  fNeutronRate = 100 * kilohertz;
  fRunTime = 100 * microsecond;
  fPulsePeriod = 800 * ns;
  fPulseWidth = 1.5 * ns;

  fAngDist = new G4SPSAngDistribution;
  fPosDist = new G4SPSPosDistribution;
  fPosDist->SetPosDisType("Point");
  fAngDist->SetPosDistribution(fPosDist);
  fRandom = new G4SPSRandomGenerator;
  fAngDist->SetAngDistType("iso");
  fAngDist->SetBiasRndm(fRandom);

  fMessenger = new DSGeneratorLicorneMessenger(this);

  DSLog(routine) << "Licorne Generator Constructed." << endlog;
}

DSGeneratorLicorne::~DSGeneratorLicorne() {
  ;
}

void DSGeneratorLicorne::DSGeneratePrimaries(G4Event* event) {

  if (fIsFirstEvent) {
    LoadNeutronSpectra();
    fIsFirstEvent = false;
    DSLog(routine) << " Mean number of gammas per neutron: " << fGammaToNeutronRatio << endlog;
    DSLog(routine) << " Neutron Rate: " << fNeutronRate / hertz << " Hz" << endlog;
    DSLog(routine) << " Run Time: " << fRunTime / microsecond << " us" << endlog;

    DSEventHandler::Get()->SetRate(fNeutronRate / kilohertz);

    if (fPulseMode) {
      DSLog(routine) << " Mean number of neutrons per Run: " << fNeutronRate * fRunTime << endlog;
      DSLog(routine) << " Mean number of neutrons per pulse: " << fNeutronRate * fPulseWidth << endlog;
    }
  }

  if (!fPulseMode) {  // Single Neutron Mode

    fParticle = fParticleTable->FindParticle(2112);
    fPosition = GetVParticlePosition();

    double myKinEne = 0;
    double myAngle = 0;
    PickEneAngle(myKinEne, myAngle);

    double phi = twopi * G4UniformRand();
    double sinPhi = std::sin(phi);
    double cosPhi = std::cos(phi);

    fDirection.setZ(std::sin(myAngle) * cosPhi);
    fDirection.setY(std::sin(myAngle) * sinPhi);
    fDirection.setX(std::cos(myAngle));

    double myMass = fParticle->GetPDGMass();
    double myEnergy = myKinEne + myMass;
    double myTime = 0.;

    double myPmon = std::sqrt(myEnergy * myEnergy - myMass * myMass);
    double myPx = myPmon * fDirection.x();
    double myPy = myPmon * fDirection.y();
    double myPz = myPmon * fDirection.z();

    G4PrimaryVertex* myVertex = new G4PrimaryVertex(fPosition, myTime);  // Position, time
    G4PrimaryParticle* myParticle = new G4PrimaryParticle(fParticle, myPx, myPy, myPz);

    myVertex->SetPrimary(myParticle);
    event->AddPrimaryVertex(myVertex);

    DSEventHandler::Get()->SetEnergy(myKinEne);
    DSEventHandler::Get()->SetPosition(fPosition / cm);
    DSEventHandler::Get()->SetDirection(fDirection);
    DSEventHandler::Get()->SetPDG(fParticle->GetPDGEncoding());
    DSEventHandler::Get()->SetTime(myTime / ns);

    int gamma_counter = 0;

    int nGammas = G4Poisson(fGammaToNeutronRatio);
    for (int i = 0; i < nGammas; ++i) {

      fParticle = fParticleTable->FindParticle(22);
      G4double gEnergy = 478 * keV;
      G4double gTime = G4UniformRand() * fRunTime;
      G4ThreeVector gDirection = fAngDist->GenerateOne();
      G4double gPx = gEnergy * gDirection.x();
      G4double gPy = gEnergy * gDirection.y();
      G4double gPz = gEnergy * gDirection.z();

      G4PrimaryVertex* gVertex = new G4PrimaryVertex(fPosition, gTime);  // Position, time
      G4PrimaryParticle* gParticle = new G4PrimaryParticle(fParticle, gPx, gPy, gPz);
      gVertex->SetPrimary(gParticle);
      event->AddPrimaryVertex(gVertex);

      DSEventHandler::Get()->SetDId(gamma_counter);
      DSEventHandler::Get()->SetDPosition(fPosition / cm);
      DSEventHandler::Get()->SetDDirection(gDirection);
      DSEventHandler::Get()->SetDPDG(fParticle->GetPDGEncoding());
      DSEventHandler::Get()->SetDTime(gTime);
      DSEventHandler::Get()->SetDEnergy(gEnergy / keV);
      DSEventHandler::Get()->SetDaughters();
      gamma_counter++;
    }

  } else {  // Pulse Mode
    fPosition = GetVParticlePosition();

    int nPulses = int(fRunTime / (fPulsePeriod) + 0.5);
    double mean_neutron_per_pulse = fNeutronRate * fPulsePeriod;

    fParticle = fParticleTable->FindParticle(2112);

    int counter = 0;

    for (int i = 0; i < nPulses; ++i) {
      int nNeutrons = G4Poisson(mean_neutron_per_pulse);
      double myTime = i * fPulsePeriod;
      for (int k = 0; k < nNeutrons; ++k) {

        double myKinEne = 0;
        double myAngle = 0;
        PickEneAngle(myKinEne, myAngle);

        double phi = twopi * G4UniformRand();
        double sinPhi = std::sin(phi);
        double cosPhi = std::cos(phi);

        fDirection.setZ(std::sin(myAngle) * cosPhi);
        fDirection.setY(std::sin(myAngle) * sinPhi);
        fDirection.setX(std::cos(myAngle));

        double myMass = fParticle->GetPDGMass();
        double myEnergy = myKinEne + myMass;

        double myPmon = std::sqrt(myEnergy * myEnergy - myMass * myMass);
        double myPx = myPmon * fDirection.x();
        double myPy = myPmon * fDirection.y();
        double myPz = myPmon * fDirection.z();

        double myLocalTime = myTime + 1.5 * ns * G4UniformRand();

        G4PrimaryVertex* myVertex = new G4PrimaryVertex(fPosition, myLocalTime);
        G4PrimaryParticle* myParticle = new G4PrimaryParticle(fParticle, myPx, myPy, myPz);

        myVertex->SetPrimary(myParticle);
        event->AddPrimaryVertex(myVertex);

        DSEventHandler::Get()->SetDId(counter);
        DSEventHandler::Get()->SetDPosition(fPosition / cm);
        DSEventHandler::Get()->SetDDirection(fDirection);
        DSEventHandler::Get()->SetDPDG(fParticle->GetPDGEncoding());
        DSEventHandler::Get()->SetDTime(myLocalTime / ns);
        DSEventHandler::Get()->SetDEnergy(myKinEne / keV);
        DSEventHandler::Get()->SetDaughters();
        counter++;
      }
    }

    fParticle = fParticleTable->FindParticle(22);
    double mean_n_gamma = fNeutronRate * fRunTime * fGammaToNeutronRatio;
    int nGammas = G4Poisson(mean_n_gamma);
    for (int i = 0; i < nGammas; ++i) {

      G4double gEnergy = 478 * keV;
      G4double gTime = G4UniformRand() * fRunTime;
      G4ThreeVector gDirection = fAngDist->GenerateOne();
      G4double gPx = gEnergy * gDirection.x();
      G4double gPy = gEnergy * gDirection.y();
      G4double gPz = gEnergy * gDirection.z();

      G4PrimaryVertex* gVertex = new G4PrimaryVertex(fPosition, gTime);  // Position, time
      G4PrimaryParticle* gParticle = new G4PrimaryParticle(fParticle, gPx, gPy, gPz);
      gVertex->SetPrimary(gParticle);
      event->AddPrimaryVertex(gVertex);

      DSEventHandler::Get()->SetDId(counter);
      DSEventHandler::Get()->SetDPosition(fPosition / cm);
      DSEventHandler::Get()->SetDDirection(gDirection);
      DSEventHandler::Get()->SetDPDG(22);
      DSEventHandler::Get()->SetDTime(gTime);
      DSEventHandler::Get()->SetDEnergy(gEnergy / keV);
      DSEventHandler::Get()->SetDaughters();
      counter++;
    }
  }

  DSEventHandler::Get()->SetUserFloat1(fPulseWidth / ns);
  DSEventHandler::Get()->SetUserFloat2(fPulsePeriod / ns);
  DSEventHandler::Get()->SetUserDouble(fRunTime / ns);
  DSEventHandler::Get()->SetUsers();

  if (event->GetNumberOfPrimaryVertex() == 0) event->SetEventAborted();
}

void DSGeneratorLicorne::LoadNeutronSpectra() {

  // from Licorne simulation
  //  angle from 0 to 7 degrees (700 bins)
  //  energy from 1 to 2 MeV   (25 bins)
  //  probability already normalized in text file

  ifstream mySpectrum("../data/physics/licorneSpectra.dat");
  if (!mySpectrum.is_open()) {
    DSLog(error) << "Couldn't load the Licorne spectrum. File not found." << endl;
    DSLog(fatal) << "" << endlog;
  }

  double ene;
  double angle;
  double weight;
  int ene_bin;
  int angle_bin;
  fMax_prob = 0;
  fEnergy_width = 0.04;
  fAngle_width = 0.01;

  while (mySpectrum >> angle >> ene >> weight) {
    ene_bin = int(ene * fEnergy_tot_bin - 25. + 0.00001);
    angle_bin = int(angle * 100 + 0.00001);  // fAngle_tot_bin - 7;
    if (weight > fMax_prob) fMax_prob = weight;
    fLicorne_array[ene_bin][angle_bin] = weight;
  }

  mySpectrum.close();
}

void DSGeneratorLicorne::PickEneAngle(G4double& ene1, G4double& angle1) {

  G4double angle_pick = G4UniformRand();
  G4int angle_bin = angle_pick * fAngle_tot_bin;
  G4double ene_pick = G4UniformRand();
  G4int ene_bin = ene_pick * fEnergy_tot_bin;
  G4double prob = G4UniformRand();

  while (prob > fLicorne_array[ene_bin][angle_bin]) {
    angle_pick = G4UniformRand();
    ene_pick = G4UniformRand();
    angle_bin = angle_pick * fAngle_tot_bin;
    ene_bin = ene_pick * fEnergy_tot_bin;
    prob = G4UniformRand();
  }

  ene1 = (ene_pick + 1.) * MeV;
  angle1 = angle_pick * 7 * degree;

  return;
}
