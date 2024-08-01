#include <iostream>
#include <stdio.h>

#include "G4Orb.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"


#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
#include "G4ThreeVector.hh"

#include "DSDetectorNeutronReD.hh"
#include "DSIO.hh"
#include "DSLogger.hh"
#include "DSMaterial.hh"
#include "DSParameters.hh"
#include "DSStorage.hh"
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4Cons.hh"
#include "G4VisAttributes.hh"

using namespace std;

////////////////////////////////////////////////
//////        Detector Description      ////////
////////////////////////////////////////////////
DSDetectorNeutronReD::DSDetectorNeutronReD(G4VPhysicalVolume* myMotherVolume) {

  DSLog(routine) << " Constructing ReD Setup Geometry : " << endlog;

  G4Colour myWhite(1.0, 1.0, 1.0);      // white
  G4Colour myGray(0.5, 0.5, 0.5);       // gray
  G4Colour myltGray(0.25, 0.25, 0.25);  // light gray
  G4Colour mydkGray(0.75, 0.75, 0.75);  // dark gray
  G4Colour myBlack(0.0, 0.0, 0.0);      // black
  G4Colour myRed(1.0, 0.0, 0.0);        // red
  G4Colour mydkRed(0.5, 0.0, 0.0);      // dark red
  G4Colour myGreen(0.0, 1.0, 0.0);      // green
  G4Colour myBlue(0.0, 0.0, 1.0);       // blue
  G4Colour myCyan(0.0, 1.0, 1.0);       // cyan
  G4Colour myMagenta(1.0, 0.0, 1.0);    // magenta
  G4Colour myYellow(1.0, 1.0, 0.0);     // yellow

  fMotherVolume = myMotherVolume;

  G4bool myCheckOverlap = DSStorage::Get()->GetCheckOverlap();

  //  G4ThreeVector ND_pointing_toward (0, 0, 0);

  //----------------------------------------------------------------------------
  // Get the neutron detector positions
  // The neutron detector position is now refered to the window position
  //----------------------------------------------------------------------------
  //  vector<G4ThreeVector> ND_position;
  vector<G4RotationMatrix*> ND_rotation;

  if (DSStorage::Get()->GetReDNDPositioningMode() == DSStorage::std || DSStorage::Get()->GetReDNDPositioningMode() == DSStorage::kinematic) {
    /*
    ND_pointing_toward.setX( 0 );
    ND_pointing_toward.setY( 0 );
    ND_pointing_toward.setZ( 0 );
    */
    switch (DSStorage::Get()->GetReDNDSurveyMode()) {
      case DSStorage::xyz:
        parseSurveyFile(DSStorage::Get()->GetReDNDSurveyFileName() /*, ND_position*/);
        //      getPositionFromXYZSurvey(DSStorage::Get()->GetReDNDSurveyFileName(),
        //      ND_position);
        DSLog(routine) << "Define neutron detector position from xyz coordinates" << endlog;
        break;

      case DSStorage::spherical:
        //      getPositionFromSphericalSurvey(DSStorage::Get()->GetReDNDSurveyFileName(),
        //      ND_position);
        DSLog(routine) << "Define neutron detector position from spherical coordinates" << endlog;
        DSLog(error) << "Code is broken.  If you really need it, ask Michael "
                        "Kuss to fix it!"
                     << endlog;
        std::exit(1);
        break;

      default:
        DSLog(error) << "Unrecognized survey mode." << endlog;
        DSLog(error) << "Exiting !!" << endlog;
        exit(1);
        break;
    }

    for (std::vector<LSci>::iterator it = lsciCol.begin(); it != lsciCol.end(); ++it) {

      if (it->GetType() == "")  // type is unset, use the default 3in
        it->SetType("3in");

      if (!it->GetRotation()) {  // rotation matrix pointer is not set
        G4RotationMatrix* myNeutronRotation = new G4RotationMatrix;
        const G4ThreeVector center(0 * cm, 0 * cm, 0 * cm);
        const G4ThreeVector dir = center - it->GetPosition();
        const G4double phi = 180.0 * degree - dir.phi();
        const G4double theta = dir.theta();
        switch (DSStorage::Get()->GetReDNDPositioningMode()) {
          case DSStorage::std:
            // Davide pointed implicitely to (0,0,0)
            //          myNeutronRotation->rotateZ(it->GetPosition().phi());
            //          myNeutronRotation->rotateY(180*degree-it->GetPosition().theta());
            myNeutronRotation->rotateZ(phi);
            myNeutronRotation->rotateY(theta);
            break;
          case DSStorage::kinematic:
            myNeutronRotation->rotateY(90 * degree);
            break;
          default:
            DSLog(error) << "Unrecognized positioning mode." << endlog;
            DSLog(error) << "Exiting !!" << endlog;
            exit(1);
            break;
        }
        it->SetRotation(myNeutronRotation);
        //      ND_rotation.push_back(myNeutronRotation);
      }
    }
  }

  DSLog(routine) << "There are " << lsciCol.size() << " neutron detectors to construct." << endlog;

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  // Define the neutron detector construction
  //----------------------------------------------------------------------------

  DSLog(routine) << "Building the Scionix EJ309 detectors." << endlog;

  // Loop through number of neutron detectors
  for (std::vector<LSci>::const_iterator it = lsciCol.begin(); it != lsciCol.end(); ++it) {

    // The dummy type just places dummy volumes to have some means of debugging
    // when placing LScis. It doesn't have any other function.
    if (it->GetType() == "dummy") {
      G4Orb* dummySolid = new G4Orb("dummySolid", 1 * cm);
      G4LogicalVolume* dummyLogic = new G4LogicalVolume(dummySolid, DSMaterial::Get()->GetAir(), "dummyLogic");
      dummyLogic->SetVisAttributes(myBlack);
      // I notice that in the entire NeutronReD class I use G4PVPlacement
      // without assigning to a variable
      new G4PVPlacement(NULL, it->GetPosition(), "dummyPhys", dummyLogic, fMotherVolume, false, 0, myCheckOverlap);
      continue;
      // dummy volume placed?  Let's go!
    }

    // dimensions from Rescigno's drawing VS-0499-91

    // detector volume envelope
    G4double LSciR = 101.6 * mm / 2.0;
    G4double LSciH = (325. + 15.) * mm / 2.0;  // 15 mm includes also connectors, for beauty
    G4double LSciZ = LSciH;
    // ATTENTION: +LSciZ is close, -LSciZ far from the origin.

    // Liquid scintillator cell
    // cell thin part
    G4double cellR = 79.25 * mm / 2.0;
    G4double cellH = 70.9 * mm / 2.0;
    G4double cellD = 1.52 * mm;
    // cell thick part
    G4double cellRR = 86. * mm / 2.0;
    G4double cellHH = 21.3 * mm / 2.0;  // ( 98.6 - 6.4 - 70.9 ) / 2
    //  G4double cellDD = cellRR - cellR + cellD;
    // cell entrance window
    G4double cellWR = cellR;
    G4double cellWH = cellD / 2.0;
    G4double cellWZ = LSciZ - cellWH;
    // cell cylinder
    G4double cellCr = cellR - cellD;
    G4double cellCR = cellR;
    G4double cellCH = cellH - cellWH;
    G4double cellCZ = cellWZ - cellWH - cellCH;
    // cell flange
    G4double cellFr = cellCr;
    G4double cellFR = cellRR;
    G4double cellFH = cellHH;
    G4double cellFZ = cellCZ - cellCH - cellFH;

    // PMT
    G4double pmtD = 1.0 * mm / 2.0;  // measured/estimated from drawing
    G4double pmtR1 = 101.6 * mm / 2.0;
    G4double pmtR2 = 83.2 * mm / 2.0;
    G4double pmtR3 = 58.8 * mm / 2.0;
    // flange (actually I am merging cell and PMT flange here!)
    G4double pmtFR = pmtR1;
    G4double pmtFr = pmtR2 - pmtD;
    G4double pmtFH = 15.4 * mm / 2.0;  // ( 107.6 - 98.6 ) + 6.4
    G4double pmtFZ = cellFZ - cellFH - pmtFH;
    // cylinder
    G4double pmtCR = pmtR2;
    G4double pmtCr = pmtCR - pmtD;
    G4double pmtCH = 50. * mm / 2.0;
    G4double pmtCZ = pmtFZ - pmtFH - pmtCH;
    // seam
    G4double pmtSr1 = pmtCr;
    G4double pmtSR1 = pmtCR;
    G4double pmtSR2 = pmtR3;
    G4double pmtSr2 = pmtR3 - pmtD;
    G4double pmtSH = 10.4 * mm / 2.0;  // measured/estimated from drawing
    G4double pmtSZ = pmtCZ - pmtCH - pmtSH;
    // housing
    G4double pmtHr = pmtSr2;
    G4double pmtHR = pmtSR2;
    G4double pmtHH = (157. - pmtD) * mm / 2.0;  // measured/estimated from the drawing, pmtD
                                                // subtracts thickness of an end plate
    G4double pmtHZ = pmtSZ - pmtSH - pmtHH;
    // magnetic shield
    G4double pmtMr = pmtHR;
    G4double pmtMR = pmtMr + 0.64 * mm / 2.0;
    G4double pmtMH = 102. * mm / 2.0;        // measured/estimated from the drawing
    G4double pmtMZ = pmtSZ - pmtSH - pmtMH;  // reference is the seam
    G4Material* pmtMmat = DSMaterial::Get()->GetSteel();
    // closing plate
    G4double pmtPR = pmtR3;
    G4double pmtPH = pmtD / 2.0;
    G4double pmtPZ = pmtHZ - pmtHH - pmtPH;
    // connectors (just for beauty)
    G4double connR = 5. * mm;
    G4double connH = 15. * mm / 2.0;
    G4double connD = 14. * mm;
    G4double connZ = pmtPZ - pmtPH - connH;
    // active volume
    G4double actR = cellCr;
    G4double actH = cellCH + cellFH;         // spans over cylinder and flange ring
    G4double actZ = cellWZ - cellWH - actH;  // reference to cell window
    // glass window (in active volume)
    G4double winR = 32. * mm;                                    // measured/estimated from the drawing
    G4double winH = 16. * mm / 2.0;                              // measured/estimated from the drawing
    G4double winZ = -actH + winH;                                // reference to active volume, but place it above lower
                                                                 // end which means inside (standard is below)
    G4Material* winMat = DSMaterial::Get()->GetBoroSiliGlass();  // should be BK7, I guess
                                                                 // BoroSiliGlass is fine

    if (it->GetType() == "5in") {
      // dimensions from Rescigno's drawing VS-0813-150

      // detector volume envelope
      LSciR = 176. * mm / 2.0;
      LSciH = (415.9 + 66.5 + 15.) * mm / 2.0;  // 15 mm includes also connectors, for beauty
      LSciZ = LSciH;
      // ATTENTION: +LSciZ is close, -LSciZ far from the origin.

      // Liquid scintillator cell
      // cell thin part
      cellR = 130. * mm / 2.0;
      cellH = 145.8 * mm / 2.0;
      cellD = 1.52 * mm;
      // cell thick part
      cellRR = 141. * mm / 2.0;
      cellHH = (165.4 - 2 * cellH - 7.2) * mm / 2.0;
      //  G4double cellDD = cellRR - cellR + cellD;
      // cell entrance window
      cellWR = cellR;
      cellWH = cellD / 2.0;
      cellWZ = LSciZ - cellWH;
      // cell cylinder
      cellCr = cellR - cellD;
      cellCR = cellR;
      cellCH = cellH - cellWH;
      cellCZ = cellWZ - cellWH - cellCH;
      // cell flange
      cellFr = cellCr;
      cellFR = cellRR;
      cellFH = cellHH;
      cellFZ = cellCZ - cellCH - cellFH;

      // PMT
      pmtD = 1. * mm / 2.0;  // measured/estimated from drawing
      pmtR1 = 160. * mm / 2.0;
      pmtR2 = 139.3 * mm / 2.0;
      pmtR3 = 58.3 * mm / 2.0;
      // flange (actually I am merging cell and PMT flange here!)
      // WHICH IS NOT OK FOR THE 5"
      // LSCI!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! OPTICAL
      // COUPLING IS MISSING IN
      // BOTH!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      pmtFR = pmtR1;
      pmtFr = pmtR2 - pmtD;
      pmtFH = 7.2 * mm / 2.0;  // I AM NEGLECTING THE RED RING!
      pmtFZ = cellFZ - cellFH - pmtFH;
      // cylinder
      pmtCR = pmtR2;
      pmtCr = pmtCR - pmtD;
      pmtCH = 167. * mm / 2.0;
      pmtCZ = pmtFZ - pmtFH - pmtCH;
      // seam
      pmtSr1 = pmtCr;
      pmtSR1 = pmtCR;
      pmtSR2 = pmtR3;
      pmtSr2 = pmtR3 - pmtD;
      pmtSH = 40. * mm / 2.0;  // measured/estimated from drawing
      pmtSZ = pmtCZ - pmtCH - pmtSH;
      // housing
      pmtHr = pmtSr2;
      pmtHR = pmtSR2;
      pmtHH = (43. + 66.5 - pmtD) * mm / 2.0;  // measured/estimated from the drawing, pmtD subtracts
                                               // thickness of an end plate
      pmtHZ = pmtSZ - pmtSH - pmtHH;
      // magnetic shield
      // difference:
      // for 3" it's on the housing of the tube,
      // for 5" on the cylinder part
      pmtMr = pmtCR;
      pmtMR = pmtMr + 0.64 * mm / 2.0;  // I put it around the cylinder!
      pmtMH = pmtCH;
      pmtMZ = pmtCZ;
      pmtMmat = DSMaterial::Get()->GetSteel();
      // closing plate
      pmtPR = pmtR3;
      pmtPH = pmtD / 2.0;
      pmtPZ = pmtHZ - pmtHH - pmtPH;
      // connectors (just for beauty)
      connR = 5. * mm;
      connH = 15. * mm / 2.0;
      connD = 14. * mm;
      connZ = pmtPZ - pmtPH - connH;
      // active volume
      actR = cellCr;
      actH = cellCH + cellFH;         // spans over cylinder and flange ring
      actZ = cellWZ - cellWH - actH;  // reference to cell window
      // glass window (in active volume)
      winR = 110. * mm / 2.0;
      winH = 32. * mm / 2.0;
      winZ = -actH + winH;                             // reference to active volume, but place it above
                                                       // lower end which means inside (standard is below)
      winMat = DSMaterial::Get()->GetBoroSiliGlass();  // should be BK7, I guess
                                                       // BoroSili glass is fine
    }

    // reference point (shift between reference point and the center of the
    // envelope solid)
    //    G4double refZ = actZ;         // active volume
    G4double refZ = actZ + winH;  // "real" active volume.  The BK7 window is in the active
                                  // volume, thus replaces some liquid.
    //    G4double refZ = LSciH;        // entrance window

    DSLog(routine) << "LSci detector extension "
                   << "height " << 2. * LSciH << " z " << LSciZ << " ( " << cellWZ + cellWH << " ) to " << connZ - connH << " reference z " << refZ << endlog;

    //--------------------------------//
    //     the "envelope" volume      //
    //--------------------------------//
    //    G4Box* refSolid = new G4Box("refSolid", 1.1*LSciR, 1.1*LSciR, 1.0*mm);
    G4Box* refSolid = new G4Box("refSolid", 1.0 * mm, 1.0 * mm, 1.0 * mm);
    G4Tubs* frontSolid = new G4Tubs("frontSolid", 0.0, cellR, cellH, 0.0, 2 * M_PI);
    G4Tubs* backSolid = new G4Tubs("backSolid", 0.0, LSciR, LSciH - cellH, 0.0, 2 * M_PI);
    // I have found no other way to redefine the reference point for a volume.
    // Seems geant4 decides that (for symmetric volumes it's the center), and
    // doesn't allow to change it. refSolid defines the center, and solidLSci
    // it's dimensions.
    G4UnionSolid* extraSolid = new G4UnionSolid("refSolid+frontSolid", refSolid, frontSolid, 0, G4ThreeVector(0, 0, LSciH - cellH - refZ));
    G4UnionSolid* LSciSolid = new G4UnionSolid("LSciSolid", extraSolid, backSolid, 0, G4ThreeVector(0, 0, -cellH - refZ));
    G4LogicalVolume* logicLSci = new G4LogicalVolume(LSciSolid, DSMaterial::Get()->GetAir(), "LSciLogic");
    logicLSci->SetVisAttributes(G4Color(1, 1, 1, 0.1));
    // logicLSci->SetVisAttributes(G4VisAttributes::Invisible);
    G4ThreeVector pos = it->GetPosition();
    char name[256];
    G4int counter = it - lsciCol.begin();
    snprintf(name,30, "LSciPhys%d", counter);

    DSLog(routine) << "LSci neutron detector " << counter << " placement at " << pos << " r " << pos.r() << endlog;
    new G4PVPlacement(it->GetRotation(), pos, name, logicLSci, fMotherVolume, false, 0, myCheckOverlap);

    //----------------------------------------------------------------------------
    // Solid Geometries
    //----------------------------------------------------------------------------

    // make all this a polycone one day ... maybe
    G4Tubs* solidLSciCellWin = new G4Tubs("LSciCellWin_Solid", 0.0, cellWR, cellWH, 0.0, 2 * M_PI);

    G4Tubs* solidLSciCellCyl = new G4Tubs("LSciCellCyl_Solid", cellCr, cellCR, cellCH, 0.0, 2 * M_PI);

    G4Tubs* solidLSciCellFlange = new G4Tubs("LSciCellFlange_Solid", cellFr, cellFR, cellFH, 0.0, 2 * M_PI);

    G4Tubs* solidLSciPMTFlange = new G4Tubs("LSciPMTFlange_Solid", pmtFr, pmtFR, pmtFH, 0.0, 2 * M_PI);

    G4Tubs* solidLSciPMTCell = new G4Tubs("LSciPMTCell_Solid", pmtCr, pmtCR, pmtCH, 0.0, 2 * M_PI);

    G4Cons* solidLSciPMTSeam = new G4Cons("LSciPMTSeam_Solid", pmtSr2, pmtSR2, pmtSr1, pmtSR1, pmtSH, 0.0, 2 * M_PI);

    G4Tubs* solidLSciPMTHouse = new G4Tubs("LSciPMTHousing_Solid", pmtHr, pmtHR, pmtHH, 0.0, 2 * M_PI);

    G4Tubs* solidLSciPMTMag = new G4Tubs("LSciPMTMagShield_Solid", pmtMr, pmtMR, pmtMH, 0.0, 2 * M_PI);

    G4Tubs* solidLSciPMTend = new G4Tubs("LSciPMTend_Solid", 0.0, pmtPR, pmtPH, 0.0, 2 * M_PI);

    G4Tubs* solidLSciConn = new G4Tubs("LSciConn_Solid", 0.0, connR, connH, 0.0, 2 * M_PI);

    G4Tubs* solidLSciActive = new G4Tubs("LSciActive_Solid", 0., actR, actH, 0., 2 * M_PI);

    G4Tubs* solidLSciWin = new G4Tubs("LSciWin_Solid", 0., winR, winH, 0., 2 * M_PI);

    G4LogicalVolume* logicLSciCellWin = new G4LogicalVolume(solidLSciCellWin, DSMaterial::Get()->GetAluminum(), "LSciCellWin_Logic");
    G4VisAttributes* hsciVis = new G4VisAttributes(myRed);
    hsciVis->SetVisibility(true);
    hsciVis->SetForceSolid(true);
    logicLSciCellWin->SetVisAttributes(hsciVis);

    new G4PVPlacement(NULL, G4ThreeVector(0, 0, cellWZ - refZ), logicLSciCellWin, "cellWindow", logicLSci, false, 0, myCheckOverlap);

    G4LogicalVolume* logicLSciCellCyl = new G4LogicalVolume(solidLSciCellCyl, DSMaterial::Get()->GetAluminum(), "LSciCellCyl_Logic");
    logicLSciCellCyl->SetVisAttributes(hsciVis);

    new G4PVPlacement(NULL, G4ThreeVector(0, 0, cellCZ - refZ), logicLSciCellCyl, "cellCylinder", logicLSci, false, 0, myCheckOverlap);

    G4LogicalVolume* logicLSciCellFlange = new G4LogicalVolume(solidLSciCellFlange, DSMaterial::Get()->GetAluminum(), "LSciCellFlange_Logic");
    G4VisAttributes* flangeVis = new G4VisAttributes(myBlue);

    flangeVis->SetVisibility(true);
    flangeVis->SetForceSolid(true);
    logicLSciCellFlange->SetVisAttributes(flangeVis);
    new G4PVPlacement(NULL, G4ThreeVector(0, 0, cellFZ - refZ), logicLSciCellFlange, "cellCylinderFlange", logicLSci, false, 0, myCheckOverlap);

    G4LogicalVolume* logicLSciPMTFlange = new G4LogicalVolume(solidLSciPMTFlange, DSMaterial::Get()->GetAluminum(), "LSciPMTFlange_Logic");
    logicLSciPMTFlange->SetVisAttributes(flangeVis);

    new G4PVPlacement(NULL, G4ThreeVector(0, 0, pmtFZ - refZ), logicLSciPMTFlange, "pmtFlange", logicLSci, false, 0, myCheckOverlap);

    G4LogicalVolume* logicLSciPMTCell = new G4LogicalVolume(solidLSciPMTCell, DSMaterial::Get()->GetAluminum(), "LSciPMTCell_Logic");
    G4VisAttributes* allVis = new G4VisAttributes(mydkGray);
    allVis->SetVisibility(true);
    allVis->SetForceSolid(true);
    logicLSciPMTCell->SetVisAttributes(allVis);

    new G4PVPlacement(NULL, G4ThreeVector(0, 0, pmtCZ - refZ), logicLSciPMTCell, "pmtCell", logicLSci, false, 0, myCheckOverlap);

    G4LogicalVolume* logicLSciPMTSeam = new G4LogicalVolume(solidLSciPMTSeam, DSMaterial::Get()->GetAluminum(), "LSciPMTSeam_Logic");
    logicLSciPMTSeam->SetVisAttributes(allVis);

    new G4PVPlacement(NULL, G4ThreeVector(0, 0, pmtSZ - refZ), logicLSciPMTSeam, "pmtSeam", logicLSci, false, 0, myCheckOverlap);

    G4LogicalVolume* logicLSciPMTHouse = new G4LogicalVolume(solidLSciPMTHouse, DSMaterial::Get()->GetAluminum(), "LSciPMTHouse_Logic");
    logicLSciPMTHouse->SetVisAttributes(allVis);

    new G4PVPlacement(NULL, G4ThreeVector(0, 0, pmtHZ - refZ), logicLSciPMTHouse, "pmtHouse", logicLSci, false, 0, myCheckOverlap);

    G4LogicalVolume* logicLSciPMTMag = new G4LogicalVolume(solidLSciPMTMag, pmtMmat, "LSciPMTMagShield_Logic");
    logicLSciPMTMag->SetVisAttributes(allVis);

    new G4PVPlacement(NULL, G4ThreeVector(0, 0, pmtMZ - refZ), logicLSciPMTMag, "pmtMagShield", logicLSci, false, 0, myCheckOverlap);

    G4LogicalVolume* logicLSciPMTend = new G4LogicalVolume(solidLSciPMTend, DSMaterial::Get()->GetAluminum(), "LSciPMTendPlate_Logic");
    logicLSciPMTend->SetVisAttributes(allVis);

    new G4PVPlacement(NULL, G4ThreeVector(0, 0, pmtPZ - refZ), logicLSciPMTend, "pmtEndPlate", logicLSci, false, 0, myCheckOverlap);

    for (int j = 0; j < 3; ++j) {  // could these be parametrized volumes, or replicated?)
      G4LogicalVolume* logicLSciConn = new G4LogicalVolume(solidLSciConn, DSMaterial::Get()->GetStainlessSteel(), "LSciConn_Logic");
      G4VisAttributes* connVis = new G4VisAttributes(myGray);
      connVis->SetVisibility(true);
      connVis->SetForceSolid(true);
      logicLSciConn->SetVisAttributes(connVis);

      new G4PVPlacement(NULL, G4ThreeVector(std::cos(j * M_PI / 2.0) * connD, std::sin(j * M_PI / 2.0) * connD, connZ - refZ), logicLSciConn, "conn", logicLSci, false, 0, myCheckOverlap);
    }

    G4LogicalVolume* logicLSciActive = new G4LogicalVolume(solidLSciActive,
                                                           DSMaterial::Get()->GetNDReDScintillator(),  // should be EJ309
                                                           "LSciActive_Logic");
    // logicLSciActive->SetVisAttributes(G4Color(1, 0, 0, 1));
    new G4PVPlacement(NULL, G4ThreeVector(0, 0, actZ - refZ), logicLSciActive, "ActiveReDLSci", logicLSci, false, 0, myCheckOverlap);

    G4LogicalVolume* logicLSciWin = new G4LogicalVolume(solidLSciWin, winMat, "LSciWin_Logic");
    new G4PVPlacement(NULL, G4ThreeVector(0, 0, winZ), logicLSciWin, "BK7win",
                      logicLSciActive,  // in active volume!
                      false, 0, myCheckOverlap);
  }  // end loop over number of detectors
  //  DefineSurfaces(); commented because the code is empty
}

#include <G4UIcommand.hh>

void DSDetectorNeutronReD::parseSurveyFile(G4String surveyFileName) {

  ifstream fileSurvey(surveyFileName);
  if (!fileSurvey) {
    fileSurvey.close();
    DSLog(error) << "The survey file " << surveyFileName << " was not found " << endlog;
    exit(0);
  }

  std::string _s;
  while (std::getline(fileSurvey, _s)) {
    if (_s[0] == '#')  // commented line
      continue;

    // extract the elements of the input line
    std::istringstream iss(_s);
    G4String v;
    std::vector<G4String> vec;
    while (iss >> v)     // yields true if extraction succeeded
      vec.push_back(v);  // Attention: vec can have variable length

    if (vec.size() == 0)  // empty line
      continue;

    // fill the element of class "LSci"
    LSci lsci;
    G4String x, y, z;

    // the first six elements are always the position (with units)
    x = vec[0] + vec[1];
    y = vec[2] + vec[3];
    z = vec[4] + vec[5];
    const G4ThreeVector pos(G4UIcommand::ConvertToDimensionedDouble(x), G4UIcommand::ConvertToDimensionedDouble(y), G4UIcommand::ConvertToDimensionedDouble(z));
    lsci.SetPosition(pos);

    // the (optional) 7th element is the LSci type
    if (vec.size() > 6) lsci.SetType(vec[6]);

    if (vec.size() == 10 || vec.size() == 13) {  // we can construct a rotation matrix
      G4ThreeVector dir;
      if (vec.size() == 10) {  // elements 8 - 10 are a direction
        x = vec[7];
        y = vec[8];
        z = vec[9];
        dir = G4ThreeVector(G4UIcommand::ConvertToDouble(x), G4UIcommand::ConvertToDouble(y), G4UIcommand::ConvertToDouble(z));
      } else if (vec.size() == 13) {  // elements 8 - 13 are a position (with
                                      // units) to point to
        x = vec[7] + vec[8];
        y = vec[9] + vec[10];
        z = vec[11] + vec[12];
        G4ThreeVector posTo(G4UIcommand::ConvertToDimensionedDouble(x), G4UIcommand::ConvertToDimensionedDouble(y), G4UIcommand::ConvertToDimensionedDouble(z));
        dir = posTo - pos;
      }
      G4RotationMatrix* rot = new G4RotationMatrix;
      const G4double phi = 180.0 * degree - dir.phi();
      const G4double theta = dir.theta();
      rot->rotateZ(phi);
      rot->rotateY(theta);
      lsci.SetRotation(rot);
    }
    lsciCol.push_back(lsci);
  }
}

/*
void DSDetectorNeutronReD::getPositionFromSphericalSurvey(G4String
surveyFileName,vector<G4ThreeVector> &ND_position){

  ND_position.clear();

  ifstream fileSurvey(surveyFileName);
  if(!fileSurvey){
    fileSurvey.close();
    DSLog(error) << "The survey file" << surveyFileName<<" was not found " <<
endlog; exit(0);
  }
  G4String R,R_unit,theta,theta_unit,phi,phi_unit;
  while(!fileSurvey.eof()){

    fileSurvey>>R>>R_unit>>theta>>theta_unit>>phi>>phi_unit;

    if( fileSurvey.eof() ) break;

    R     += " " + R_unit;
    theta += " " + theta_unit;
    phi   += " " + phi_unit;


    G4ThreeVector v(0,0,0);
    v.setRThetaPhi(G4UIcommand::ConvertToDimensionedDouble(R),
                   G4UIcommand::ConvertToDimensionedDouble(theta),
                   G4UIcommand::ConvertToDimensionedDouble(phi));

    ND_position.push_back(v);
  }
}
*/
