#include <Randomize.hh>
#include <fstream>
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "DSEventHandler.hh"
#include "DSGeneratorMultiSpectra.hh"

DSGeneratorMultiSpectra::DSGeneratorMultiSpectra() : DSVGenerator("DSGeneratorMultiSpectra") {
  fNSpectra = 0;

  fIsotope = fUnits = 0;
  fMinEnergy = fMaxEnergy = 0;
  fSpectrumCumul = 0;

  fParticle = 0;
  fParticleTable = G4ParticleTable::GetParticleTable();

  fDirection = G4ThreeVector(0., 0., 1.);
  fPosition = G4ThreeVector(0., 0., 0.);

  fMessenger = new DSGeneratorMultiSpectraMessenger(this);

  DSLog(routine) << "Multiple Spectra Generator constructed. Set desired "
                    "number of spectra by "
                    "/ds/generator/multispectra/nspectra N"
                 << endlog;
}

DSGeneratorMultiSpectra::~DSGeneratorMultiSpectra() {
  if (fIsotope) delete[] fIsotope;
  if (fUnits) delete[] fUnits;
  if (fMinEnergy) delete[] fMinEnergy;
  if (fMaxEnergy) delete[] fMaxEnergy;
  // if (fSpectrumEne)   delete[] fSpectrumEne;
  if (fSpectrumCumul) delete[] fSpectrumCumul;
}

void DSGeneratorMultiSpectra::DSGeneratePrimaries(G4Event* event) {
  double ku = 1.;

  DSLog(routine) << "Generated particles in each vertex: " << endlog;

  for (int n = 0; n < fNSpectra; ++n) {
    fParticle = fParticleTable->FindParticle(fIsotope[n]);

    fDirection = GetVParticleDirection();
    fPosition = GetVParticlePosition();

    G4String fu = fUnits[n];
    if (fu == "eV") ku = eV;
    else if (fu == "keV")
      ku = keV;
    else if (fu == "MeV")
      ku = MeV;
    else if (fu == "GeV")
      ku = GeV;
    else if (fu == "TeV")
      ku = TeV;
    else if (fu == "PeV")
      ku = PeV;
    else {
      DSLog(error) << "DSGeneratorMultiSpectra::DSGeneratePrimaries: units not "
                      "recognized for spectrum "
                   << n << "; skipping vertex generation." << endl;
      continue;
    }

    G4double myPDG = fParticle->GetPDGEncoding();
    G4double myKinEne = pickAnEnergy(n) * ku;

    G4double myMass = fParticle->GetPDGMass();
    G4double myEnergy = myKinEne + myMass;
    G4double myTime = 0.;

    G4double myPmon = std::sqrt(myEnergy * myEnergy - myMass * myMass);
    G4double myPx = myPmon * fDirection.x();
    G4double myPy = myPmon * fDirection.y();
    G4double myPz = myPmon * fDirection.z();

    G4PrimaryVertex* myVertex = new G4PrimaryVertex(fPosition, myTime);  // Position, time
    G4PrimaryParticle* myParticle = new G4PrimaryParticle(fParticle, myPx, myPy, myPz);

    DSLog(routine) << n << " " << fParticle->GetParticleName() << " " << myKinEne << " MeV" << endlog;

    myVertex->SetPrimary(myParticle);
    event->AddPrimaryVertex(myVertex);

    DSEventHandler::Get()->SetEnergy(myKinEne / keV);
    DSEventHandler::Get()->SetPosition(fPosition / cm);
    DSEventHandler::Get()->SetDirection(fDirection);
    DSEventHandler::Get()->SetPDG(myPDG);
    DSEventHandler::Get()->SetTime(myTime / ns);
  }
}

G4int DSGeneratorMultiSpectra::loadSpectrum(int n) {
  char dir[2000];
  snprintf(dir, 50, "../data/physics/MultiSpectra/spectrum_%.3d.dat", n);

  ifstream inpFile(dir);

  if (!inpFile.is_open()) {
    DSLog(fatal) << "Couldn't load spectrum " << dir << "; terminating." << endl;
    return -1;
  }

  G4double val, sum = 0.;

  // header
  inpFile >> fIsotope[n] >> fMinEnergy[n] >> fMaxEnergy[n] >> fUnits[n];

  // content
  while (inpFile >> val) {
    sum += val;

    // fSpectrumEne[n].push_back(val);
    fSpectrumCumul[n].push_back(sum);
  }

  inpFile.close();

  // normalize
  for (vector<G4double>::iterator it = fSpectrumCumul[n].begin(); it != fSpectrumCumul[n].end(); ++it) *it /= sum;

  return n;
}

G4double DSGeneratorMultiSpectra::pickAnEnergy(int n) {
  assert(fSpectrumCumul[n].size());
  vector<G4double>::iterator it_cum = fSpectrumCumul[n].begin();

  G4double E0 = fMinEnergy[n];
  G4double dE = (fMaxEnergy[n] - fMinEnergy[n]) / fSpectrumCumul[n].size();

  G4double p0 = 0., dp;
  G4double p = p0;

  do p = G4UniformRand();
  while (!p);

  while (it_cum != fSpectrumCumul[n].end()) {
    dp = *it_cum - p0;

    if (*it_cum > p) return E0 + dE * (p - p0) / dp;

    E0 += dE;
    p0 = *it_cum;
    ++it_cum;
  }

  return 0.;
}

void DSGeneratorMultiSpectra::initArrays() {
  if (fNSpectra > 0) {
    if (fIsotope) delete[] fIsotope;
    fIsotope = new G4String[fNSpectra];

    if (fUnits) delete[] fUnits;
    fUnits = new G4String[fNSpectra];

    if (fMinEnergy) delete[] fMinEnergy;
    fMinEnergy = new G4double[fNSpectra];

    if (fMaxEnergy) delete[] fMaxEnergy;
    fMaxEnergy = new G4double[fNSpectra];

    //	  if (fSpectrumEne)
    //	  {
    //		 delete[] fSpectrumEne;
    //	  }
    //	  fSpectrumEne = new vector<G4double>[fNSpectra];

    if (fSpectrumCumul) { delete[] fSpectrumCumul; }
    fSpectrumCumul = new vector<G4double>[fNSpectra];

    for (int i = 0; i < fNSpectra; ++i) loadSpectrum(i);
  }
}
