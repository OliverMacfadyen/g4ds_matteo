//---------------------------------------------------------------------------//
//                                                                           //
//                                                                           //
//                         G4DS Simulation                                   //
//                                                                           //
// --------------------------------------------------------------------------//

#include "DSDetectorConstruction.hh"
#include "DSEventAction.hh"
#include "DSEventHandler.hh"
#include "DSLogger.hh"
#include "DSManager.hh"
#include "DSPhysicsList.hh"
#include "DSPrimaryGeneratorAction.hh"
#include "DSRunAction.hh"
#include "DSStackingAction.hh"
#include "DSSteppingAction.hh"
#include "DSStorage.hh"
#include "DSTrackingAction.hh"

#include "DSIO.hh"

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIsession.hh"
#include "G4UItcsh.hh"
#include "G4UIterminal.hh"

#include "DSVisManager.hh"
#include "G4TrajectoryDrawByParticleID.hh"
#include "G4VisExecutive.hh"

#include "G4UIExecutive.hh"

#ifdef G4UI_USE_ROOT
#include "G4UIRoot.hh"
#endif

/*#ifdef G4UI_USE_XM
#include "G4UIXm.hh"
#endif
*/

#ifdef G4LIB_USE_GDML
#include "G4GDMLParser.hh"
#endif

#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <string>
#include "G4ios.hh"

//Global variable to count the amount of photons
DSManager* runManager;
DSDetectorConstruction* theDetector;
DSPhysicsList* thePhysicsList;
// Functions called from main().
void PrintHeader();
void PrintUsage(void);

using namespace std;

int main(int argc, char** argv) {

  PrintHeader();
  G4bool IsStack = false;

  //---------------------------------------------------------------------------//
  // Definition of promary parameters
  DSLogger::SetSeverity(DSLogger::routine);
  DSIO::Get()->CheckFileName("output");
  string filename;
  if (argc > 1) {
    ifstream ifs(argv[1]);
    string _s;
    while (getline(ifs, _s)) {
      if (!_s.find("/DSlog")) {
        DSLogger::SetSeverity(DSLogger::toEnum(_s.substr(_s.find(" ") + 1)));
      } else if (!_s.find("/run/filename")) {
        filename = _s.substr(_s.find(" ") + 1);
      } else if (!_s.find("/ds/stack")) {
        IsStack = true;
      } else if (!_s.find("/run/beamOn")) {
        int nEventsToGenerate = std::atoi(_s.substr(_s.find(" ") + 1).c_str());
        DSEventHandler::Get()->SetEvents(nEventsToGenerate);
      }
    }
    ifs.close();
  }
  //if(argc==3)  filename = argv[2];
  DSIO::Get()->CheckFileName(filename);

  DSLog(routine) << "Output file name: " << DSIO::Get()->GetFileName() << endlog;
  // Copy fo the input file in the log file
  DSIO::Get()->OpenLogFiles();
  DSIO::Get()->GetStreamLogFile() << endl;
  DSIO::Get()->GetStreamLogFile() << "------------------------------------------" << endl;
  DSIO::Get()->GetStreamLogFile() << "                   G4DS                   " << endl;
  DSIO::Get()->GetStreamLogFile() << endl;
  DSIO::Get()->GetStreamLogFile() << "     The Dark Side Geant4 Simulator        " << endl;
  DSIO::Get()->GetStreamLogFile() << "------------------------------------------" << endl;
  DSIO::Get()->GetStreamLogFile() << endl;
  DSIO::Get()->GetStreamLogFile() << endl;
  DSIO::Get()->GetStreamLogFile() << "------------------------------------------" << endl;
  DSIO::Get()->GetStreamLogFile() << "                   Input                  " << endl;
  DSIO::Get()->GetStreamLogFile() << "------------------------------------------" << endl;
  DSIO::Get()->GetStreamLogFile() << endl;
  if (argc > 1) {
    ifstream ifs(argv[1]);
    string _s;
    while (getline(ifs, _s)) DSIO::Get()->GetStreamLogFile() << _s << endl;
    DSIO::Get()->GetStreamLogFile() << endl;
    DSIO::Get()->GetStreamLogFile() << endl;
    DSIO::Get()->GetStreamLogFile() << "------------------------------------------" << endl;
    DSIO::Get()->GetStreamLogFile() << "                   Output                 " << endl;
    DSIO::Get()->GetStreamLogFile() << "------------------------------------------" << endl;
    DSIO::Get()->GetStreamLogFile() << endl;

    ifs.close();
  }

  //---------------------------------------------------------------------------//

  DSLog(trace) << "Creating G4 Run Manager" << endlog;
  runManager = DSManager::Get();

  // Register detector geometry and materials.
  DSLog(trace) << "Creating and registering G4 geometry" << endlog;
  DSDetectorConstruction* myDetector = new DSDetectorConstruction;
  runManager->SetUserInitialization(myDetector);

  // Register Geant4 physics processes
  DSLog(trace) << "Creating and registering G4  physics processes" << endlog;
  DSPhysicsList* myPhysicsList = new DSPhysicsList;
  runManager->SetUserInitialization(myPhysicsList);

  // Register event generator.
  DSLog(trace) << "Creating and registering event generator" << endlog;
  DSPrimaryGeneratorAction* myGenerator = new DSPrimaryGeneratorAction;
  runManager->SetUserAction(myGenerator);

  // Register run action. What to do at beginning and end of each run.
  DSLog(trace) << "Registering G4 run action." << endlog;
  DSRunAction* myRunAction = new DSRunAction;
  runManager->SetUserAction(myRunAction);

  // Register event action, ie. what to save/compute for each event.
  DSLog(trace) << "Registering G4 event action." << endlog;
  DSEventAction* myEventAction = new DSEventAction;
  runManager->SetUserAction(myEventAction);

  //Register tracking action
  DSLog(trace) << "Registering G4 tracking action." << endlog;
  DSTrackingAction* myTrackAction = new DSTrackingAction;
  runManager->SetUserAction(myTrackAction);

  // Register stepping action, ie. what to save/compute for each step.
  DSLog(trace) << "Registering G4 stepping action." << endlog;
  DSSteppingAction* myStepAction = new DSSteppingAction;
  runManager->SetUserAction(myStepAction);

  if (IsStack) {
    // Register stacking action, ie. what to save/compute for each step.
    DSLog(trace) << "Registering G4 stacking action." << endlog;
    DSStackingAction* myStackAction = new DSStackingAction;
    runManager->SetUserAction(myStackAction);
  }

  // visualization manager
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
  G4TrajectoryDrawByParticleID* model = new G4TrajectoryDrawByParticleID;

  model->SetDefault("cyan");
  model->Set("gamma", "red");
  model->Set("neutron", "blue");
  model->Set("opticalphoton", "green");
  model->Set("e+", "magenta");
  model->Set("e-", G4Colour(0.3, 0.3, 0.3));

  visManager->RegisterModel(model);

  visManager->SelectTrajectoryModel(model->Name());

  if (argc > 1 && !strcmp(argv[1], "-h")) {
    PrintUsage();
  } else if (argc > 2 && !strcmp(argv[2], "Qt")) {
    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    G4UIExecutive* ui = new G4UIExecutive(argc, argv);
    G4cout << "Entering Qt visu mode...\n";
    G4cout << argv[1] << " mac file have been loaded " << '\n';
    G4String fileName = argv[1];
    G4String command = "/control/execute ";
    UImanager->ApplyCommand(command + fileName);
    ui->SessionStart();
    delete ui;
  } else if (argc > 2 && !strcmp(argv[2], "session")) {
    G4cout << "Entering batch mode...\n";
    G4cout << "Executing script file from command line: " << argv[1] << '\n';
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    G4UImanager::GetUIpointer()->ApplyCommand(command + fileName);
    G4UIsession* session = new G4UIterminal;
    session->SessionStart();
    delete session;
  } else if (argc > 2 && !strcmp(argv[2], "geo")) {

    G4cout << "Entering geo mode...\n";
    G4cout << "Executing script file from command line: " << argv[1] << '\n';
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    G4UImanager::GetUIpointer()->ApplyCommand(command + fileName);
    string geofilename = DSIO::Get()->GetFileName();
    ;
    geofilename += ".gdml";
    G4cout << "Write geometry to " << geofilename << '\n';

#ifdef G4LIB_USE_GDML
    G4GDMLParser parser;
    parser.Write(geofilename.c_str());
#endif

  } else {
    G4cout << "Entering batch mode...\n";
    G4cout << "Executing script file from command line: " << argv[1] << '\n';
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    G4UImanager::GetUIpointer()->ApplyCommand(command + fileName);
  }
}

//---------------------------------------------------------------------------//

void PrintHeader(void) {
  G4cout << "Dark Side Monte Carlo Simulation" << G4endl;
  G4cout << "------------------------------------------" << G4endl;
  G4cout << "Version: 2.00" << G4endl;
  G4cout << "based on Geant4.10.00.p01" << G4endl;
}

//---------------------------------------------------------------------------//

void PrintUsage(void) {
  G4cout << "Usage:" << G4endl;
  G4cout << "g4DS -h : Displays this message" << G4endl;
  G4cout << "g4DS <filename> : Executes script <filename>" << G4endl;
  G4cout << "g4DS <filename> session: execute G4DS interactively after having executed the script <filename>" << G4endl
         << G4endl;
  G4cout << "g4DS <filename> geo: execute G4DS and produce a gdml geometry file" << G4endl << G4endl;
}

//---------------------------------------------------------------------------//
