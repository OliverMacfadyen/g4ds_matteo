#ifndef DSManagerMessenger_h
#define DSManagerMessenger_h 1

#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UImessenger.hh"
#include "globals.hh"

using namespace std;

class DSManager;
class G4UIcmdWithAString;

class DSManagerMessenger : public G4UImessenger {

 public:
  DSManagerMessenger(DSManager*);
  ~DSManagerMessenger();

  void SetNewValue(G4UIcommand*, G4String);

 private:
  DSManager *fManager; 

  G4UIdirectory* fDirectory;
  G4UIcmdWithAString* fDSLogCmd;
  G4UIcmdWithAString* fDSOpticsCmd;
  G4UIcmdWithAString* fLArPropertiesCmd;
  G4UIcmdWithABool* fOverLapCmd;
  G4UIcmdWithABool* fGDMLCmd;
  G4UIcmdWithABool* fWritePhotonsCmd;
  G4UIcmdWithAnInteger* fEventCounterCmd;
  G4UIcmdWithAnInteger* fVerbosityCmd;
  G4UIcmdWithABool* fWriteDaughtersCmd;
  G4UIcmdWithABool* fWriteDepositsCmd;
  G4UIcmdWithABool* fWriteThermalElectronsCmd;
  G4UIcmdWithAnInteger* fDaughterDepthCmd;
  G4UIcmdWithADouble* fTMBfractionCmd;
  G4UIcmdWithADoubleAndUnit* fTimeCutCmd;
  G4UIcmdWithADoubleAndUnit* fWriteFilterMinECmd;
  G4UIcmdWithADoubleAndUnit* fWriteFilterMaxECmd;
  G4UIcmdWithAnInteger* fFastSimulationCmd;
  G4UIcmdWithABool* fStoreEmittedCherenkovCmd;
  G4UIcmdWithABool* fSkipEventsWithNoDaughtersCmd;
  G4UIcmdWithABool* fWriteFilterCmd;
  G4UIcmdWithAString* fChooseFilterCmd; 
 
  G4UIcmdWithADouble* fMultiplicatorCmd;
};

#endif
/*
 * $Log: DSManagerMessenger.hh,v $
 * Revision 1.5  2015/10/22 10:18:45  dfranco
 * added Cherenov process (by default off)
 *
 * Revision 1.4  2015/06/02 16:25:54  pagnes
 * added command to use a different LArScintillationProperties.txt file
 * (/ds/manager/lar_properties)
 *
 * Revision 1.3  2014/07/23 14:52:38  pagnes
 * write thermal e- and kill S1 commands added
 *
 * Revision 1.2  2014/07/16 08:23:13  pagnes
 * QE scaling to 1.0 added (/ds/manager/fast_simulation xxx)
 *
 * Revision 1.1  2014/05/07 12:20:53  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.5  2014/04/11 12:33:35  pagnes
 * command to set TMB/PC ratio inside veto added
 *
 * Revision 1.4  2014/03/19 16:37:36  dfranco
 * update external configuration reader for optics tuning
 *
 * Revision 1.3  2013/03/22 16:24:40  dfranco
 * added a command to set the "daughter depth" level to store in the output file
 *
 * Revision 1.2  2013/03/22 14:09:40  dfranco
 * Added the cvs logger code at the end of each file
 *
 */
