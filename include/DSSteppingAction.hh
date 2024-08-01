#ifndef DSSteppingAction_h
#define DSSteppingAction_h 1

#include "DSParameters.hh"
#include "G4UserSteppingAction.hh"

class  DSBiasStopStep;


class DSSteppingAction : public G4UserSteppingAction {
 public:
  DSSteppingAction();
  ~DSSteppingAction(){};

  void UserSteppingAction(const G4Step*);

 private:
  DSBiasStopStep *fStopStep ;

};

#endif
/*
 * $Log: DSSteppingAction.hh,v $
 * Revision 1.1  2014/05/07 12:20:55  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.3  2013/08/27 07:13:12  swesterd
 * add visible energy for the neutron veto
 *
 * Revision 1.2  2013/03/22 14:09:40  dfranco
 * Added the cvs logger code at the end of each file
 *
 */
