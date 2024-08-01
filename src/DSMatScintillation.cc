// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// Scintillation Light Class Implementation
//
// File:        DSMatScintillation.cc
// Description: RestDiscrete Process - Generation of Scintillation Photons
// Version:     1.0
// Created:     1998-11-07
// Author:      Peter Gumplinger
// Updated:     2010-10-20 Allow the scintillation yield to be a function
//              of energy deposited by particle type
//              Thanks to Zach Hartwig (Department of Nuclear
//              Science and Engineeering - MIT)
//              2010-09-22 by Peter Gumplinger
//              > scintillation rise time included, thanks to
//              > Martin Goettlich/DESY
//              2005-08-17 by Peter Gumplinger
//              > change variable name MeanNumPhotons -> MeanNumberOfPhotons
//              2005-07-28 by Peter Gumplinger
//              > add G4ProcessType to constructor
//              2004-08-05 by Peter Gumplinger
//              > changed StronglyForced back to Forced in GetMeanLifeTime
//              2002-11-21 by Peter Gumplinger
//              > change to use G4Poisson for small MeanNumberOfPhotons
//              2002-11-07 by Peter Gumplinger
//              > now allow for fast and slow scintillation component
//              2002-11-05 by Peter Gumplinger
//              > now use scintillation constants from G4Material
//              2002-05-09 by Peter Gumplinger
//              > use only the PostStepPoint location for the origin of
//                scintillation photons when energy is lost to the medium
//                by a neutral particle
//              2000-09-18 by Peter Gumplinger
//              > change: aSecondaryPosition=x0+rand*aStep.GetDeltaPosition();
//                        aSecondaryTrack->SetTouchable(0);
//              2001-09-17, migration of Materials to pure STL (mma)
//              2003-06-03, V.Ivanchenko fix compilation warnings
//
 
#include "DSMatScintillation.hh"
 
#include "globals.hh"
#include "G4DynamicParticle.hh"
#include "G4EmProcessSubType.hh"
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4MaterialPropertyVector.hh"
#include "G4OpticalParameters.hh"
#include "G4ParticleMomentum.hh"
#include "G4ParticleTypes.hh"
#include "G4PhysicalConstants.hh"
#include "G4PhysicsFreeVector.hh"
#include "G4PhysicsTable.hh"
#include "G4Poisson.hh"
#include "G4ScintillationTrackInformation.hh"
#include "G4StepPoint.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"
#include "G4PhysicsModelCatalog.hh"
 

 using namespace std ;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DSMatScintillation::DSMatScintillation(const G4String& processName,
                                 G4ProcessType type)
  : G4VRestDiscreteProcess(processName, type)
  , fIntegralTable1(nullptr)
  , fIntegralTable2(nullptr)
  , fIntegralTable3(nullptr)
  , fIntegralTable4(nullptr)
  , fEmSaturation(nullptr)
  , fNumPhotons(0)
{
  secID = G4PhysicsModelCatalog::GetModelID("model_Scintillation");
  SetProcessSubType(fScintillation);
 
#ifdef G4DEBUG_SCINTILLATION
  ScintTrackEDep  = 0.;
  ScintTrackYield = 0.;
#endif
 
  if(verboseLevel > 0)
  {
    G4cout << GetProcessName() << " is created " << G4endl;
  }

  Initialise();
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DSMatScintillation::~DSMatScintillation()
{
  if(fIntegralTable1 != nullptr)
  {
    fIntegralTable1->clearAndDestroy();
    delete fIntegralTable1;
  }
  if(fIntegralTable2 != nullptr)
  {
    fIntegralTable2->clearAndDestroy();
    delete fIntegralTable2;
  }
  if(fIntegralTable3 != nullptr)
  {
    fIntegralTable3->clearAndDestroy();
    delete fIntegralTable3;
  }
  if(fIntegralTable4 != nullptr)
  {
    fIntegralTable4->clearAndDestroy();
    delete fIntegralTable4;
  }
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DSMatScintillation::ProcessDescription(std::ostream& out) const
{
  out << "Scintillation simulates production of optical photons produced\n"
         "by a high energy particle traversing matter.\n"
         "Various material properties need to be defined.\n";
  G4VRestDiscreteProcess::DumpInfo();
 
  G4OpticalParameters* params = G4OpticalParameters::Instance();
  out << "Track secondaries first: " << params->GetScintTrackSecondariesFirst();
  out << "Finite rise time: " << params->GetScintFiniteRiseTime();
  out << "Save track information: " << params->GetScintTrackInfo();
  out << "Stack photons: " << params->GetScintStackPhotons();
  out << "Verbose level: " << params->GetScintVerboseLevel();
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool DSMatScintillation::IsApplicable(const G4ParticleDefinition& aParticleType)
{
  if(aParticleType.GetParticleName() == "opticalphoton")
    return false;
  if(aParticleType.IsShortLived())
    return false;
  return true;
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DSMatScintillation::PreparePhysicsTable(const G4ParticleDefinition&)
{
  Initialise();
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DSMatScintillation::Initialise()
{
  G4OpticalParameters* params = G4OpticalParameters::Instance();
  SetTrackSecondariesFirst(params->GetScintTrackSecondariesFirst());
  SetFiniteRiseTime(params->GetScintFiniteRiseTime());
  SetScintillationTrackInfo(params->GetScintTrackInfo());
  SetStackPhotons(params->GetScintStackPhotons());
  SetVerboseLevel(params->GetScintVerboseLevel());
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DSMatScintillation::BuildPhysicsTable(const G4ParticleDefinition&)
{
  if(fIntegralTable1 != nullptr)
  {
    fIntegralTable1->clearAndDestroy();
    delete fIntegralTable1;
    fIntegralTable1 = nullptr;
  }
  if(fIntegralTable2 != nullptr)
  {
    fIntegralTable2->clearAndDestroy();
    delete fIntegralTable2;
    fIntegralTable2 = nullptr;
  }
  if(fIntegralTable3 != nullptr)
  {
    fIntegralTable3->clearAndDestroy();
    delete fIntegralTable3;
    fIntegralTable3 = nullptr;
  }
  if(fIntegralTable4 != nullptr)
  {
    fIntegralTable4->clearAndDestroy();
    delete fIntegralTable4;
    fIntegralTable4 = nullptr;
  } 
  const G4MaterialTable* materialTable = G4Material::GetMaterialTable();
  size_t numOfMaterials                = G4Material::GetNumberOfMaterials();
 
  // create new physics table
  if(!fIntegralTable1)
    fIntegralTable1 = new G4PhysicsTable(numOfMaterials);
  if(!fIntegralTable2)
    fIntegralTable2 = new G4PhysicsTable(numOfMaterials);
  if(!fIntegralTable3)
    fIntegralTable3 = new G4PhysicsTable(numOfMaterials);
  if(!fIntegralTable4)
    fIntegralTable4 = new G4PhysicsTable(numOfMaterials);
 
  for(size_t i = 0; i < numOfMaterials; ++i)
  {
    G4PhysicsFreeVector* vector1 = new G4PhysicsFreeVector();
    G4PhysicsFreeVector* vector2 = new G4PhysicsFreeVector();
    G4PhysicsFreeVector* vector3 = new G4PhysicsFreeVector();
    G4PhysicsFreeVector* vector4 = new G4PhysicsFreeVector();
 
    // Retrieve vector of scintillation wavelength intensity for
    // the material from the material's optical properties table.
    G4MaterialPropertiesTable* MPT =
      ((*materialTable)[i])->GetMaterialPropertiesTable();
 
    if(MPT)
    {
      G4MaterialPropertyVector* MPV =
        MPT->GetProperty("MATSCINTILLATIONCOMPONENT1");
      if(MPV)
      {
        // Retrieve the first intensity point in vector
        // of (photon energy, intensity) pairs
        G4double currentIN = (*MPV)[0];
        if(currentIN >= 0.0)
        {
          // Create first (photon energy, Scintillation Integral pair
          G4double currentPM  = MPV->Energy(0);
          G4double currentCII = 0.0;
          vector1->InsertValues(currentPM, currentCII);
 
          // Set previous values to current ones prior to loop
          G4double prevPM  = currentPM;
          G4double prevCII = currentCII;
          G4double prevIN  = currentIN;
 
          // loop over all (photon energy, intensity)
          // pairs stored for this material
          for(size_t ii = 1; ii < MPV->GetVectorLength(); ++ii)
          {
            currentPM = MPV->Energy(ii);
            currentIN = (*MPV)[ii];
            currentCII =
              prevCII + 0.5 * (currentPM - prevPM) * (prevIN + currentIN);
 
            vector1->InsertValues(currentPM, currentCII);
 
            prevPM  = currentPM;
            prevCII = currentCII;
            prevIN  = currentIN;
          }
        }
      }
 
      MPV = MPT->GetProperty("MATSCINTILLATIONCOMPONENT2");
      if(MPV)
      {
        // Retrieve the first intensity point in vector
        // of (photon energy, intensity) pairs
        G4double currentIN = (*MPV)[0];
        if(currentIN >= 0.0)
        {
          // Create first (photon energy, Scintillation Integral pair
          G4double currentPM  = MPV->Energy(0);
          G4double currentCII = 0.0;
          vector2->InsertValues(currentPM, currentCII);
 
          // Set previous values to current ones prior to loop
          G4double prevPM  = currentPM;
          G4double prevCII = currentCII;
          G4double prevIN  = currentIN;
 
          // loop over all (photon energy, intensity)
          // pairs stored for this material
          for(size_t ii = 1; ii < MPV->GetVectorLength(); ++ii)
          {
            currentPM = MPV->Energy(ii);
            currentIN = (*MPV)[ii];
            currentCII =
              prevCII + 0.5 * (currentPM - prevPM) * (prevIN + currentIN);
 
            vector2->InsertValues(currentPM, currentCII);
 
            prevPM  = currentPM;
            prevCII = currentCII;
            prevIN  = currentIN;
          }
        }
      }
      MPV = MPT->GetProperty("MATSCINTILLATIONCOMPONENT3");
      if(MPV)
      {
        // Retrieve the first intensity point in vector
        // of (photon energy, intensity) pairs
        G4double currentIN = (*MPV)[0];
        if(currentIN >= 0.0)
        {
          // Create first (photon energy, Scintillation Integral pair
          G4double currentPM  = MPV->Energy(0);
          G4double currentCII = 0.0;
          vector3->InsertValues(currentPM, currentCII);
 
          // Set previous values to current ones prior to loop
          G4double prevPM  = currentPM;
          G4double prevCII = currentCII;
          G4double prevIN  = currentIN;
 
          // loop over all (photon energy, intensity)
          // pairs stored for this material
          for(size_t ii = 1; ii < MPV->GetVectorLength(); ++ii)
          {
            currentPM = MPV->Energy(ii);
            currentIN = (*MPV)[ii];
            currentCII =
              prevCII + 0.5 * (currentPM - prevPM) * (prevIN + currentIN);
 
            vector3->InsertValues(currentPM, currentCII);
 
            prevPM  = currentPM;
            prevCII = currentCII;
            prevIN  = currentIN;
          }
        }
      }
      MPV = MPT->GetProperty("MATSCINTILLATIONCOMPONENT4");
      if(MPV)
      {
        // Retrieve the first intensity point in vector
        // of (photon energy, intensity) pairs
        G4double currentIN = (*MPV)[0];
        if(currentIN >= 0.0)
        {
          // Create first (photon energy, Scintillation Integral pair
          G4double currentPM  = MPV->Energy(0);
          G4double currentCII = 0.0;
          vector4->InsertValues(currentPM, currentCII);
 
          // Set previous values to current ones prior to loop
          G4double prevPM  = currentPM;
          G4double prevCII = currentCII;
          G4double prevIN  = currentIN;
 
          // loop over all (photon energy, intensity)
          // pairs stored for this material
          for(size_t ii = 1; ii < MPV->GetVectorLength(); ++ii)
          {
            currentPM = MPV->Energy(ii);
            currentIN = (*MPV)[ii];
            currentCII =
              prevCII + 0.5 * (currentPM - prevPM) * (prevIN + currentIN);
 
            vector3->InsertValues(currentPM, currentCII);
 
            prevPM  = currentPM;
            prevCII = currentCII;
            prevIN  = currentIN;
          }
        }
      }
    }

    fIntegralTable1->insertAt(i, vector1);
    fIntegralTable2->insertAt(i, vector2);
    fIntegralTable3->insertAt(i, vector3);
    fIntegralTable4->insertAt(i, vector4);
  }
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VParticleChange* DSMatScintillation::AtRestDoIt(const G4Track& aTrack,
                                               const G4Step& aStep)
// This routine simply calls the equivalent PostStepDoIt since all the
// necessary information resides in aStep.GetTotalEnergyDeposit()
{
  return DSMatScintillation::PostStepDoIt(aTrack, aStep);
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VParticleChange* DSMatScintillation::PostStepDoIt(const G4Track& aTrack,
                                                 const G4Step& aStep)
// This routine is called for each tracking step of a charged particle
// in a scintillator. A Poisson/Gauss-distributed number of photons is
// generated according to the scintillation yield formula, distributed
// evenly along the track segment and uniformly into 4pi.
{
  aParticleChange.Initialize(aTrack);
  fNumPhotons = 0;
 
  const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
  const G4Material* aMaterial        = aTrack.GetMaterial();
 
  G4StepPoint* pPreStepPoint  = aStep.GetPreStepPoint();
  G4StepPoint* pPostStepPoint = aStep.GetPostStepPoint();
 
  G4ThreeVector x0 = pPreStepPoint->GetPosition();
  G4ThreeVector p0 = aStep.GetDeltaPosition().unit();
  G4double t0      = pPreStepPoint->GetGlobalTime();
 
  G4double TotalEnergyDeposit = aStep.GetTotalEnergyDeposit();
 
  G4MaterialPropertiesTable* MPT = aMaterial->GetMaterialPropertiesTable();
  if(!MPT)
    return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);
 
  G4int N_timeconstants = 1;
 
  if(MPT->GetProperty("MATSCINTILLATIONCOMPONENT4"))
    N_timeconstants = 4;
  else if(MPT->GetProperty("MATSCINTILLATIONCOMPONENT3"))
    N_timeconstants = 3;
  else if(MPT->GetProperty("MATSCINTILLATIONCOMPONENT2"))
    N_timeconstants = 2;
  else if(!(MPT->GetProperty("MATSCINTILLATIONCOMPONENT1")))
  {
    // no components were specified
    return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);
  }
 
  G4double ResolutionScale = MPT->GetConstProperty("MATRESOLUTIONSCALE");
  G4double MeanNumberOfPhotons;
 
  G4double yield1     = 0.;
  G4double yield2     = 0.;
  G4double yield3     = 0.;
  G4double yield4     = 0.;
  G4double sum_yields = 0.;
 
  yield1 = MPT->ConstPropertyExists("MATSCINTILLATIONYIELD1")
             ? MPT->GetConstProperty("MATSCINTILLATIONYIELD1")
             : 1.;
  yield2 = MPT->ConstPropertyExists("MATSCINTILLATIONYIELD2")
             ? MPT->GetConstProperty("MATSCINTILLATIONYIELD2")
             : 0.;
  yield3 = MPT->ConstPropertyExists("MATSCINTILLATIONYIELD3")
             ? MPT->GetConstProperty("MATSCINTILLATIONYIELD3")
             : 0.;
  yield4 = MPT->ConstPropertyExists("MATSCINTILLATIONYIELD4")
             ? MPT->GetConstProperty("MATSCINTILLATIONYIELD4")
             : 0.;
  // The default linear scintillation process
  // Units: [# scintillation photons / MeV]
  MeanNumberOfPhotons = MPT->GetConstProperty("MATSCINTILLATIONYIELD");
  // Birk's correction via fEmSaturation and specifying scintillation by
  // by particle type are physically mutually exclusive
  if(fEmSaturation)
    MeanNumberOfPhotons *=
      (fEmSaturation->VisibleEnergyDepositionAtAStep(&aStep));
  else
    MeanNumberOfPhotons *= TotalEnergyDeposit;
  
  sum_yields = yield1 + yield2 + yield3 + yield4;


  if(MeanNumberOfPhotons > 10.)
  {
    G4double sigma = ResolutionScale * std::sqrt(MeanNumberOfPhotons);
    fNumPhotons = G4int(G4RandGauss::shoot(MeanNumberOfPhotons, sigma) + 0.5);
  }
  else
  {
    fNumPhotons = G4int(G4Poisson(MeanNumberOfPhotons));
  }
 
  if(fNumPhotons <= 0 || !fStackingFlag)
  {
    // return unchanged particle and no secondaries
    aParticleChange.SetNumberOfSecondaries(0);
    return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);
  }
  
  // DF: 
  //    +10 photons is introduced to increase the buffer
  //    this is needed because binomial shoots can sometime 
  //    fluctuate NumPhot > fNumPhotons 
  aParticleChange.SetNumberOfSecondaries(fNumPhotons+10);
 
  if(fTrackSecondariesFirst)
  {
    if(aTrack.GetTrackStatus() == fAlive)
      aParticleChange.ProposeTrackStatus(fSuspend);
  }
 
  G4int materialIndex = aMaterial->GetIndex();
 
  // Retrieve the Scintillation Integral for this material
  // new G4PhysicsFreeVector allocated to hold CII's
  size_t numPhot                     = fNumPhotons;
  G4double scintTime                 = 0.;
  G4double riseTime                  = 0.;
  G4PhysicsFreeVector* scintIntegral = nullptr;
  G4ScintillationType scintType      = Slow;
  fGeneratedPhotons                  = 0 ;
  
  for(G4int scnt = 0; scnt < N_timeconstants; ++scnt)
  {
    // if there is 1 time constant it is #1, etc.
    if(scnt == 0)
    {
      if(N_timeconstants == 1)
      {
        numPhot = fNumPhotons;
      }
      else
      {
        numPhot = CLHEP::RandBinomial::shoot(fNumPhotons, yield1 / sum_yields);
      }
      scintTime = MPT->GetConstProperty("MATSCINTILLATIONTIMECONSTANT1");
      if(fFiniteRiseTime)
      {
        riseTime = MPT->GetConstProperty("MATSCINTILLATIONRISETIME1");
      }
      scintType = Fast;
      scintIntegral =
        (G4PhysicsFreeVector*) ((*fIntegralTable1)(materialIndex));
      fGeneratedPhotons += numPhot;
    }
    else if(scnt == 1)
    {
      // to be consistent with old version (due to double->int conversion)
      if(N_timeconstants == 2)
      {
        if (fNumPhotons - fGeneratedPhotons > 0)
          numPhot = fNumPhotons - fGeneratedPhotons;
        else numPhot = 0;
      }
      else
      {
        //numPhot = yield2 / sum_yields * fNumPhotons;
        numPhot = CLHEP::RandBinomial::shoot(fNumPhotons, yield2 / sum_yields);
      }
      scintTime = MPT->GetConstProperty("MATSCINTILLATIONTIMECONSTANT2");
      if(fFiniteRiseTime)
      {
        riseTime = MPT->GetConstProperty("MATSCINTILLATIONRISETIME2");
      }
      scintType = Medium;
      scintIntegral =
        (G4PhysicsFreeVector*) ((*fIntegralTable2)(materialIndex));
     fGeneratedPhotons += numPhot;
    }
    else if(scnt == 2)
    {
      if(N_timeconstants == 3)
      {
        if (fNumPhotons - fGeneratedPhotons > 0)
          numPhot = fNumPhotons - fGeneratedPhotons;
        else numPhot = 0;
      }
      else
      {
        //numPhot = yield2 / sum_yields * fNumPhotons;
        numPhot = CLHEP::RandBinomial::shoot(fNumPhotons, yield3 / sum_yields);
      }
      scintTime = MPT->GetConstProperty("MATSCINTILLATIONTIMECONSTANT3");
      if(fFiniteRiseTime)
      {
        riseTime = MPT->GetConstProperty("MATSCINTILLATIONRISETIME3");
      }
      scintType = Slow;
      scintIntegral =
        (G4PhysicsFreeVector*) ((*fIntegralTable3)(materialIndex));
     fGeneratedPhotons += numPhot;
    }
    else if(scnt == 3)
    {
      if (fNumPhotons - fGeneratedPhotons > 0)
        numPhot = fNumPhotons - fGeneratedPhotons;
      else numPhot = 0;

      scintTime = MPT->GetConstProperty("MATSCINTILLATIONTIMECONSTANT4");
      if(fFiniteRiseTime)
      {
        riseTime = MPT->GetConstProperty("MATSCINTILLATIONRISETIME4");
      }
      scintType = Slow;
      scintIntegral =
        (G4PhysicsFreeVector*) ((*fIntegralTable4)(materialIndex));

     fGeneratedPhotons += numPhot;
    }
 
    if(!scintIntegral)
      continue;
    

    G4double CIImax = scintIntegral->GetMaxValue();
    for(size_t i = 0; i < numPhot; ++i)
    {
      // Determine photon energy
      G4double CIIvalue      = G4UniformRand() * CIImax;
      G4double sampledEnergy = scintIntegral->GetEnergy(CIIvalue);
 
      if(verboseLevel > 1)
      {
        G4cout << "sampledEnergy = " << sampledEnergy << G4endl;
        G4cout << "CIIvalue =        " << CIIvalue << G4endl;
      }
 
      // Generate random photon direction
      G4double cost = 1. - 2. * G4UniformRand();
      G4double sint = std::sqrt((1. - cost) * (1. + cost));
      G4double phi  = twopi * G4UniformRand();
      G4double sinp = std::sin(phi);
      G4double cosp = std::cos(phi);
      G4ParticleMomentum photonMomentum(sint * cosp, sint * sinp, cost);
 
      // Determine polarization of new photon
      G4ThreeVector photonPolarization(cost * cosp, cost * sinp, -sint);
      G4ThreeVector perp = photonMomentum.cross(photonPolarization);
      phi                = twopi * G4UniformRand();
      sinp               = std::sin(phi);
      cosp               = std::cos(phi);
      photonPolarization = (cosp * photonPolarization + sinp * perp).unit();
 
      // Generate a new photon:
      G4DynamicParticle* scintPhoton =
        new G4DynamicParticle(opticalphoton, photonMomentum);
      scintPhoton->SetPolarization(photonPolarization);
      scintPhoton->SetKineticEnergy(sampledEnergy);
 
      // Generate new G4Track object:
      G4double rand = G4UniformRand();
      if(aParticle->GetDefinition()->GetPDGCharge() == 0)
      {
        rand = 1.0;
      }
 
      // emission time distribution
      G4double delta = rand * aStep.GetStepLength();
      G4double deltaTime =
        delta /
        (pPreStepPoint->GetVelocity() +
         rand * (pPostStepPoint->GetVelocity() - pPreStepPoint->GetVelocity()) /
           2.);
      if(riseTime == 0.0)
      {
        deltaTime -= scintTime * std::log(G4UniformRand());
      }
      else
      {
        deltaTime += sample_time(riseTime, scintTime);
      }
 
      G4double secTime          = t0 + deltaTime;
      G4ThreeVector secPosition = x0 + rand * aStep.GetDeltaPosition();
 
      G4Track* secTrack = new G4Track(scintPhoton, secTime, secPosition);
      secTrack->SetTouchableHandle(
        aStep.GetPreStepPoint()->GetTouchableHandle());
      secTrack->SetParentID(aTrack.GetTrackID());
      secTrack->SetCreatorModelID(secID);
      if(fScintillationTrackInfo)
        secTrack->SetUserInformation(
          new G4ScintillationTrackInformation(scintType));
      aParticleChange.AddSecondary(secTrack);
    }
  }
 
  if(verboseLevel > 1)
  {
    G4cout << "\n Exiting from DSMatScintillation::DoIt -- NumberOfSecondaries = "
           << aParticleChange.GetNumberOfSecondaries() << G4endl;
  }
 
  return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double DSMatScintillation::GetMeanFreePath(const G4Track&, G4double,
                                          G4ForceCondition* condition)
{
  *condition = StronglyForced;
  return DBL_MAX;
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double DSMatScintillation::GetMeanLifeTime(const G4Track&,
                                          G4ForceCondition* condition)
{
  *condition = Forced;
  return DBL_MAX;
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double DSMatScintillation::sample_time(G4double tau1, G4double tau2)
{
  // tau1: rise time and tau2: decay time
  // Loop checking, 07-Aug-2015, Vladimir Ivanchenko
  while(true)
  {
    G4double ran1 = G4UniformRand();
    G4double ran2 = G4UniformRand();
 
    // exponential distribution as envelope function: very efficient
    G4double d = (tau1 + tau2) / tau2;
    // make sure the envelope function is
    // always larger than the bi-exponential
    G4double t  = -1.0 * tau2 * std::log(1. - ran1);
    G4double gg = d * single_exp(t, tau2);
    if(ran2 <= bi_exp(t, tau1, tau2) / gg)
      return t;
  }
  return -1.0;
}
 

 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DSMatScintillation::DumpPhysicsTable() const
{
  if(fIntegralTable1)
  {
    for(size_t i = 0; i < fIntegralTable1->entries(); ++i)
    {
      ((G4PhysicsFreeVector*) (*fIntegralTable1)[i])->DumpValues();
    }
  }
  if(fIntegralTable2)
  {
    for(size_t i = 0; i < fIntegralTable2->entries(); ++i)
    {
      ((G4PhysicsFreeVector*) (*fIntegralTable2)[i])->DumpValues();
    }
  }
  if(fIntegralTable3)
  {
    for(size_t i = 0; i < fIntegralTable3->entries(); ++i)
    {
      ((G4PhysicsFreeVector*) (*fIntegralTable3)[i])->DumpValues();
    }
  }
  if(fIntegralTable4)
  {
    for(size_t i = 0; i < fIntegralTable4->entries(); ++i)
    {
      ((G4PhysicsFreeVector*) (*fIntegralTable4)[i])->DumpValues();
    }
  }
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DSMatScintillation::SetTrackSecondariesFirst(const G4bool state)
{
  fTrackSecondariesFirst = state;
  G4OpticalParameters::Instance()->SetScintTrackSecondariesFirst(
    fTrackSecondariesFirst);
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DSMatScintillation::SetFiniteRiseTime(const G4bool state)
{
  fFiniteRiseTime = state;
  G4OpticalParameters::Instance()->SetScintFiniteRiseTime(fFiniteRiseTime);
}
 

 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DSMatScintillation::SetScintillationTrackInfo(const G4bool trackType)
{
  fScintillationTrackInfo = trackType;
  G4OpticalParameters::Instance()->SetScintTrackInfo(fScintillationTrackInfo);
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DSMatScintillation::SetStackPhotons(const G4bool stackingFlag)
{
  fStackingFlag = stackingFlag;
  G4OpticalParameters::Instance()->SetScintStackPhotons(fStackingFlag);
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DSMatScintillation::SetVerboseLevel(G4int verbose)
{
  verboseLevel = verbose;
  G4OpticalParameters::Instance()->SetScintVerboseLevel(verboseLevel);
}
