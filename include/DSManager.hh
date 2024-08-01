#ifndef _BXMANAGER_HH
#define _BXMANAGER_HH

#include "G4RunManager.hh"
//---------------------------------------------------------------------------//

class DSManagerMessenger;

//---------------------------------------------------------------------------//

enum class FilterEnum { None=0, TPCEnergy=1, IVEnergy=2, Daughters=3, Deposits=4, ArDM=5, ARIS=6 };


class DSManager : public G4RunManager {
 private:
  // default constructor
  DSManager();

 public:
  static DSManager* Get();

  // copy constructor
  DSManager(const DSManager&);

  G4bool EventSurvivesFiltering(); 
  void SetFilter(G4String);
  void SetMinEnergy(G4double val) { fMinEnergy = val;}
  void SetMaxEnergy(G4double val) { fMaxEnergy = val;}

  // public interface

  // protected members
 protected:
  // private  members
 private:
  // destructor
  virtual ~DSManager();

  static DSManager* Manager;
  // Pointers to main objects.
  // G4RunManager                    *fG4RunManager;
  DSManagerMessenger* fDSMessenger;

  FilterEnum fFilterType ;

  G4bool FilterTPCEnergy();
  G4bool FilterIVEnergy();
  G4bool FilterDaughters();
  G4bool FilterDeposits();
  G4bool FilterArDM();
  G4bool FilterARIS();

  G4double fMinEnergy;
  G4double fMaxEnergy;



};
#endif
/*
 * $Log: DSManager.hh,v $
 * Revision 1.1  2014/05/07 12:20:53  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.2  2013/03/22 14:09:40  dfranco
 * Added the cvs logger code at the end of each file
 *
 */
