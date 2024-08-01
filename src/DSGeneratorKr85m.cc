#include <Randomize.hh>
#include <fstream>
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "DSEventHandler.hh"
#include "DSGeneratorKr85m.hh"

DSGeneratorKr85m::DSGeneratorKr85m() : DSVGenerator("DSGeneratorKr85m") {

  fIsFirstEvent = true;
  fSamples = 0;

  DSLog(routine) << "Kr85m Generator Constructed." << endlog;
}

DSGeneratorKr85m::~DSGeneratorKr85m() {}

void DSGeneratorKr85m::DSGeneratePrimaries(G4Event* event) {

  if (fIsFirstEvent) {

    fParticle1 = fParticleTable->FindParticle("e-");
    fParticle2 = fParticleTable->FindParticle("gamma");
    LoadAr39CrossSection();
    fIsFirstEvent = false;
  }

  G4ThreeVector fIsotopePosition = GetVParticlePosition();

  DSEventHandler::Get()->SetEnergy(0);
  DSEventHandler::Get()->SetPosition(fIsotopePosition / cm);
  DSEventHandler::Get()->SetDirection(G4ThreeVector(0, 0, 0));
  DSEventHandler::Get()->SetPDG(1000360850);
  DSEventHandler::Get()->SetTime(0);

  // Beta Kr85m
  fDirection = GetVParticleDirection();

  G4int myPDG = fParticle1->GetPDGEncoding();

  G4double myKinEne = PickAnEnergy() * keV;
  G4double myMass = fParticle1->GetPDGMass();
  G4double myEnergy = myKinEne + myMass;
  G4double myTime = 0.;

  G4double myPmon = std::sqrt(myEnergy * myEnergy - myMass * myMass);
  G4double myPx = myPmon * fDirection.x();
  G4double myPy = myPmon * fDirection.y();
  G4double myPz = myPmon * fDirection.z();

  G4PrimaryVertex* myVertex = new G4PrimaryVertex(fIsotopePosition, myTime);  // Position, time
  G4PrimaryParticle* myParticle = new G4PrimaryParticle(fParticle1, myPx, myPy, myPz);

  myVertex->SetPrimary(myParticle);
  event->AddPrimaryVertex(myVertex);

  DSEventHandler::Get()->SetDEnergy(myKinEne / keV);
  DSEventHandler::Get()->SetDPosition(fIsotopePosition / cm);
  DSEventHandler::Get()->SetDDirection(fDirection);
  DSEventHandler::Get()->SetDPDG(myPDG);
  DSEventHandler::Get()->SetDTime(0);
  DSEventHandler::Get()->SetDaughters();

  // Gamma Kr85m

  fDirection = GetVParticleDirection();
  myPDG = fParticle2->GetPDGEncoding();

  myKinEne = 514.0067 * keV;
  myEnergy = myKinEne + myMass;
  myTime = CLHEP::RandExponential::shoot(1.015 / log(2) * microsecond);

  myPx = myKinEne * fDirection.x();
  myPy = myKinEne * fDirection.y();
  myPz = myKinEne * fDirection.z();

  G4PrimaryVertex* myVertex2 = new G4PrimaryVertex(fIsotopePosition, myTime);  // Position, time
  G4PrimaryParticle* myParticle2 = new G4PrimaryParticle(fParticle2, myPx, myPy, myPz);

  myVertex2->SetPrimary(myParticle2);
  event->AddPrimaryVertex(myVertex2);

  DSEventHandler::Get()->SetDEnergy(myKinEne / keV);
  DSEventHandler::Get()->SetDPosition(fIsotopePosition / cm);
  DSEventHandler::Get()->SetDDirection(fDirection);
  DSEventHandler::Get()->SetDPDG(myPDG);
  DSEventHandler::Get()->SetDTime(myTime);
  DSEventHandler::Get()->SetDaughters();
}

void DSGeneratorKr85m::LoadAr39CrossSection() {

  ifstream fSpectrum("../data/physics/kr85m_spectrum.dat");

  if (!fSpectrum.is_open()) {
    DSLog(error) << "Couldn't load the Kr85m spectrum. File not found." << endl;
    DSLog(fatal) << "" << endlog;
  }

  while (fSpectrum >> fSpectrumEne[fSamples] >> fSpectrumCumul[fSamples]) {
    if (fSamples == 0) fSpectrumCumul[fSamples] *= fSpectrumEne[fSamples];
    else
      fSpectrumCumul[fSamples] = fSpectrumCumul[fSamples - 1] + fSpectrumCumul[fSamples] * (fSpectrumEne[fSamples] - fSpectrumEne[fSamples - 1]);
    fSamples++;
  }
  fSpectrum.close();

  for (int i = 0; i < fSamples; i++) fSpectrumCumul[i] /= fSpectrumCumul[fSamples - 1];
}

G4double DSGeneratorKr85m::PickAnEnergy() {

  G4double p = G4UniformRand();

  for (int i = 0; i < fSamples; i++) {

    if (fSpectrumCumul[i] > p) {

      G4double E0 = 0.;
      G4double dE = fSpectrumEne[0];
      G4double dp = fSpectrumCumul[0];

      if (i > 0) {
        E0 = fSpectrumEne[i - 1];
        dE = fSpectrumEne[i] - fSpectrumEne[i - 1];
        dp = fSpectrumCumul[i] - fSpectrumCumul[i - 1];
      }

      G4double E = E0 + dE * (p - fSpectrumCumul[i - 1]) / dp;
      return E;
    }
  }

  return 0.;
}
