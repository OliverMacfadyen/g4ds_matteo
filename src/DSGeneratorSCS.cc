#include <Randomize.hh>
#include <fstream>
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "DSEventHandler.hh"
#include "DSGeneratorSCS.hh"

DSGeneratorSCS::DSGeneratorSCS() : DSVGenerator("DSGeneratorSCS") {

  fIsotope = "Ar39";
  fIsFirstEvent = true;
  fParticleTable = G4ParticleTable::GetParticleTable();

  fDirection = G4ThreeVector(0., 0., 1.);
  fPosition = G4ThreeVector(0., 0., 0.);

  fMinEnergy = 0 * keV;
  fMaxEnergy = 1000 * keV;

  fSamples = 0;

  fMessenger = new DSGeneratorSCSMessenger(this);

  DSLog(routine) << "SCS Generator Constructed." << endlog;
}

DSGeneratorSCS::~DSGeneratorSCS() {}

void DSGeneratorSCS::DSGeneratePrimaries(G4Event* event) {

  if (fIsFirstEvent) {

    if (fIsotope == "Ar39") {
      fParticle = fParticleTable->FindParticle("e-");
      LoadAr39CrossSection();
    }

    fIsFirstEvent = false;
  }

  fDirection = GetVParticleDirection();
  fPosition = GetVParticlePosition();

  G4double myPDG = fParticle->GetPDGEncoding();

  G4double myKinEne = PickAnEnergy() * MeV;
  while (myKinEne < fMinEnergy || myKinEne > fMaxEnergy) myKinEne = PickAnEnergy() * MeV;

  G4double myMass = fParticle->GetPDGMass();
  G4double myEnergy = myKinEne + myMass;
  G4double myTime = 0.;

  G4double myPmon = std::sqrt(myEnergy * myEnergy - myMass * myMass);
  G4double myPx = myPmon * fDirection.x();
  G4double myPy = myPmon * fDirection.y();
  G4double myPz = myPmon * fDirection.z();

  G4PrimaryVertex* myVertex = new G4PrimaryVertex(fPosition, myTime);  // Position, time
  G4PrimaryParticle* myParticle = new G4PrimaryParticle(fParticle, myPx, myPy, myPz);

  myVertex->SetPrimary(myParticle);
  event->AddPrimaryVertex(myVertex);

  DSEventHandler::Get()->SetEnergy(myKinEne / keV);
  DSEventHandler::Get()->SetPosition(fPosition / cm);
  DSEventHandler::Get()->SetDirection(fDirection);
  DSEventHandler::Get()->SetPDG(myPDG);
  DSEventHandler::Get()->SetTime(myTime / ns);
}

void DSGeneratorSCS::LoadAr39CrossSection() {

  ifstream fSpectrum("../data/physics/Ar39spc-highstatistics.dat");

  if (!fSpectrum.is_open()) {
    DSLog(error) << "Couldn't load the Ar39 spectrum. File not found." << endl;
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

G4double DSGeneratorSCS::PickAnEnergy() {

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
      if (E >= 0) return E;
      else
        return -E;
    }
  }

  return 0.;
}
