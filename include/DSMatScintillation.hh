//
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
// Scintillation Light Class Definition
//
// File:        DSMatScintillation.hh
// Description: Discrete Process - Generation of Scintillation Photons
// Version:     1.0
// Created:     1998-11-07
// Author:      Peter Gumplinger
// Updated:     2010-10-20 Allow the scintillation yield to be a function
//                         of energy deposited by particle type
//                         Thanks to Zach Hartwig (Department of Nuclear
//                         Science and Engineeering - MIT)
//              2005-07-28 add G4ProcessType to constructor
//              2002-11-21 change to user G4Poisson for small MeanNumPotons
//              2002-11-07 allow for fast and slow scintillation
//              2002-11-05 make use of constant material properties
//              2002-05-16 changed to inherit from VRestDiscreteProcess
//              2002-05-09 changed IsApplicable method
//              1999-10-29 add method and class descriptors
//
//
 
#ifndef DSMatScintillation_h
#define DSMatScintillation_h 1
 
#include "globals.hh"
#include "G4EmSaturation.hh"
#include "G4OpticalPhoton.hh"
#include "G4VRestDiscreteProcess.hh"
 
class G4PhysicsTable;
class G4Step;
class G4Track;
 
// Class Description:
// RestDiscrete Process - Generation of Scintillation Photons.
// Class inherits publicly from G4VRestDiscreteProcess.
// Class Description - End:
 
class DSMatScintillation : public G4VRestDiscreteProcess
{
 public:
  explicit DSMatScintillation(const G4String& processName = "Scintillation",
                           G4ProcessType type          = fElectromagnetic);
  ~DSMatScintillation();
 
  DSMatScintillation(const DSMatScintillation& right) = delete;
  DSMatScintillation& operator=(const DSMatScintillation& right) = delete;
 
  // DSMatScintillation Process has both PostStepDoIt (for energy
  // deposition of particles in flight) and AtRestDoIt (for energy
  // given to the medium by particles at rest)
 
  G4bool IsApplicable(const G4ParticleDefinition& aParticleType) override;
  // Returns true -> 'is applicable', for any particle type except
  // for an 'opticalphoton' and for short-lived particles
 
  void ProcessDescription(std::ostream&) const override;
  void DumpInfo() const override {ProcessDescription(G4cout);};
 
  void BuildPhysicsTable(const G4ParticleDefinition& aParticleType) override;
  // Build table at the right time
 
  void PreparePhysicsTable(const G4ParticleDefinition& part) override;
  void Initialise();
 
  G4double GetMeanFreePath(const G4Track& aTrack, G4double,
                           G4ForceCondition*) override;
  // Returns infinity; i. e. the process does not limit the step,
  // but sets the 'StronglyForced' condition for the DoIt to be
  // invoked at every step.
 
  G4double GetMeanLifeTime(const G4Track& aTrack, G4ForceCondition*) override;
  // Returns infinity; i. e. the process does not limit the time,
  // but sets the 'StronglyForced' condition for the DoIt to be
  // invoked at every step.
 
  G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
                                  const G4Step& aStep) override;
  G4VParticleChange* AtRestDoIt(const G4Track& aTrack,
                                const G4Step& aStep) override;
 
  
  void SetTrackSecondariesFirst(const G4bool state);
  // If set, the primary particle tracking is interrupted and any
  // produced scintillation photons are tracked next. When all
  // have been tracked, the tracking of the primary resumes.
 
  G4bool GetTrackSecondariesFirst() const;
  // Returns the boolean flag for tracking secondaries first.
 
  void SetFiniteRiseTime(const G4bool state);
  // If set, the DSMatScintillation process expects the user to have
  // set the constant material property SCINTILLATIONRISETIME{1,2,3}.
 
  G4bool GetFiniteRiseTime() const;
  // Returns the boolean flag for a finite scintillation rise time.
 
  G4PhysicsTable* GetIntegralTable1() const;
  // Returns the address of scintillation integral table #1.
 
  G4PhysicsTable* GetIntegralTable2() const;
  // Returns the address of scintillation integral table #2.
 
  G4PhysicsTable* GetIntegralTable3() const;
  // Returns the address of scintillation integral table #3.
 
  G4PhysicsTable* GetIntegralTable4() const;
  // Returns the address of scintillation integral table #3.
 
  void AddSaturation(G4EmSaturation* sat);
  // Adds Birks Saturation to the process.
 
  void RemoveSaturation();
  // Removes the Birks Saturation from the process.
 
  G4EmSaturation* GetSaturation() const;
  // Returns the Birks Saturation.
 
 
 
  void SetScintillationTrackInfo(const G4bool trackType);
  // Call by the user to set the DSMatScintillationTrackInformation
  // to scintillation photon track
 
  G4bool GetScintillationTrackInfo() const;
  // Return the boolean for whether or not the
  // DSMatScintillationTrackInformation is set to the scint. photon track
 
  void SetStackPhotons(const G4bool);
  // Call by the user to set the flag for stacking the scint. photons
 
  G4bool GetStackPhotons() const;
  // Return the boolean for whether or not the scint. photons are stacked
 
  G4int GetNumPhotons() const;
  // Returns the current number of scint. photons (after PostStepDoIt)
 
  void DumpPhysicsTable() const;
  // Prints the fast and slow scintillation integral tables.
 
  void SetVerboseLevel(G4int);
  // sets verbosity
 
 private:
 
  G4PhysicsTable* fIntegralTable1;
  G4PhysicsTable* fIntegralTable2;
  G4PhysicsTable* fIntegralTable3;
  G4PhysicsTable* fIntegralTable4
  ;
 
  G4EmSaturation* fEmSaturation;
  const G4ParticleDefinition* opticalphoton =
    G4OpticalPhoton::OpticalPhotonDefinition();
 
  G4int fNumPhotons;

  G4int fGeneratedPhotons;
  
  G4bool fScintillationTrackInfo;
  G4bool fStackingFlag;
  G4bool fTrackSecondariesFirst;
  G4bool fFiniteRiseTime;
 
#ifdef G4DEBUG_SCINTILLATION
  G4double ScintTrackEDep, ScintTrackYield;
#endif
 
  G4double single_exp(G4double t, G4double tau2);
  G4double bi_exp(G4double t, G4double tau1, G4double tau2);
 
  // emission time distribution when there is a finite rise time
  G4double sample_time(G4double tau1, G4double tau2);
 
  G4int secID = -1;  // creator modelID
 
};
 
// Inline methods
 
inline G4bool DSMatScintillation::GetTrackSecondariesFirst() const
{
  return fTrackSecondariesFirst;
}
 
inline G4bool DSMatScintillation::GetFiniteRiseTime() const
{
  return fFiniteRiseTime;
}
 
inline G4PhysicsTable* DSMatScintillation::GetIntegralTable1() const
{
  return fIntegralTable1;
}
 
inline G4PhysicsTable* DSMatScintillation::GetIntegralTable2() const
{
  return fIntegralTable2;
}
 
inline G4PhysicsTable* DSMatScintillation::GetIntegralTable3() const
{
  return fIntegralTable3;
}
 
inline void DSMatScintillation::AddSaturation(G4EmSaturation* sat)
{
  fEmSaturation = sat;
}
 
inline void DSMatScintillation::RemoveSaturation() { fEmSaturation = nullptr; }
 
inline G4EmSaturation* DSMatScintillation::GetSaturation() const
{
  return fEmSaturation;
}
 
 
inline G4bool DSMatScintillation::GetScintillationTrackInfo() const
{
  return fScintillationTrackInfo;
}
 
inline G4bool DSMatScintillation::GetStackPhotons() const { return fStackingFlag; }
 
inline G4int DSMatScintillation::GetNumPhotons() const { return fNumPhotons; }
 
inline G4double DSMatScintillation::single_exp(G4double t, G4double tau2)
{
  return std::exp(-1.0 * t / tau2) / tau2;
}
 
inline G4double DSMatScintillation::bi_exp(G4double t, G4double tau1,
                                        G4double tau2)
{
  return std::exp(-1.0 * t / tau2) * (1 - std::exp(-1.0 * t / tau1)) / tau2 /
         tau2 * (tau1 + tau2);
}
 
#endif /* DSMatScintillation_h */
