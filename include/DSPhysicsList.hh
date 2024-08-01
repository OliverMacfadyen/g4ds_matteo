#ifndef DSPhysicsList_h
#define DSPhysicsList_h 1

#include "DSParameters.hh"
#include "G4VModularPhysicsList.hh"
#include "globals.hh"

class DSPhysicsListMessenger;

class DSPhysicsList : public G4VModularPhysicsList {
 public:
  DSPhysicsList();
  ~DSPhysicsList();
  virtual void SetCuts() override;

  void SetHPRangeCuts(G4bool value) { isHPRangeCuts = value; }
  G4bool GetHPRangeCuts() { return isHPRangeCuts; }

  void SetDepositCuts(G4bool value) { isDepositCuts = value; }
  G4bool GetDepositCuts() { return isDepositCuts; }

  void SetHadronicList(G4int value) { fHadronic = value; }
  G4int GetHadronicList() { return fHadronic; }

  void SetEMList(G4int value) { fEM = value; }
  G4int GetEMList() { return fEM; }

  //! \brief Construct physics process
  virtual void ConstructProcess() override;

 private:
  void EnableAugerEffect();
  void DumpProcessTable();
  void ConstructHad_Custom();
  void OptionsForHadronicPhysics();
  // void ConstructLSOp();
  // void ConstructLSOp2();

  DSPhysicsListMessenger* fMessenger;
  G4VPhysicsConstructor* fOptPhysics;
  G4VPhysicsConstructor* fHadPhysics;
  G4VPhysicsConstructor* fEMPhysics;

  G4int VerboseLevel;
  G4int OpVerbLevel;

  G4double cutForGamma;
  G4double cutForElectron;
  G4double cutForPositron;
  G4double cutForProton;
  G4double cutForAlpha;
  G4double cutForGenericIon;

  G4bool isHPRangeCuts;
  G4bool isDepositCuts;
  G4bool isHadronicHP;

  G4int fHadronic;
  G4int fEM;
};

#endif
