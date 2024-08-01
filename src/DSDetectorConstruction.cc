#include <fstream>
#include <iostream>
#include "G4PhysicalConstants.hh"
#include "G4String.hh"
#include "G4SystemOfUnits.hh"

#include "G4GeometryManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4SolidStore.hh"

#include "G4CSGSolid.hh"
#include "G4PVPlacement.hh"
#include "G4VPhysicalVolume.hh"

#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Sphere.hh"
#include "G4Tubs.hh"

// migration : canceled out BREPPolycone and added GenerciPolycone
//#include  "G4GenericPolycone.hh"
#include "G4Polycone.hh"

#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4OpticalSurface.hh"

#include "G4Point3D.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"

#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4ios.hh"

#include "G4IntersectionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"

#include "DSDetectorConstruction.hh"
#include "DSDetectorConstructionMessenger.hh"
#include "DSLogger.hh"
#include "DSMaterial.hh"

#include "DSDetectorDS10.hh"
#include "DSDetectorDS20k.hh"
#include "DSDetectorDS20kNeutronVeto.hh"
#include "DSDetectorDS50.hh"
#include "DSDetectorDuneCryostat.hh"
#include "DSDetectorNeutronVeto.hh"
#include "DSDetectorPET.hh"
#include "DSDetectorPlasticVeto.hh"
#include "DSDetectorProto.hh"
#include "DSDetectorProtoProto.hh"
#include "DSDetectorReD.hh"
#include "DSDetectorWaterTank.hh"
#include "DSDetectorArDM.hh"
#include "DSDetectorLicorne.hh"
#include "DSDetectorLicorneNeutron.hh"
#include "DSDetectorARISER.hh"
#include "DSDetectorARISERGamma.hh"
#include "DSDetectorCalibrationDevice.hh"
#include "DSDetectorDart.hh"
#include "DSDetectorGantryCalibration.hh"
#include "DSDetectorSource.hh"
#include "DSDetectorTester.hh"
#include "DSEventHandler.hh"
#include "DSStorage.hh"

#include <math.h>
#include <time.h>
#include <sstream>

using namespace std;

extern G4int thePMTconfig;

G4int NumberOfVetoPMT;

DSDetectorConstruction::DSDetectorConstruction() {

  time_t rawtime;
  time(&rawtime);
  DSLog(routine) << ctime(&rawtime) << endlog;

  fDetectorConfiguration = 0;

  //  fIsSource  = false ;

  fMessenger = new DSDetectorConstructionMessenger(this);
}

DSDetectorConstruction::~DSDetectorConstruction() {
  delete fMessenger;
  // delete gantry;
}

G4VPhysicalVolume* DSDetectorConstruction::Construct() {

  return ConstructDetector();
}

void DSDetectorConstruction::DefineSizes() {
  ;
}

G4VPhysicalVolume* DSDetectorConstruction::ConstructDetector() {
  // World
  fSolidWorld = new G4Box("World_Solid", 20 * m, 20 * m, 20 * m);
  fLogicWorld = new G4LogicalVolume(fSolidWorld, DSMaterial::Get()->GetAir(), "World_Logic");
  G4Colour myWhite(1.0, 1.0, 1.0);  // white
  // G4VisAttributes* LogVisAttWorld = new G4VisAttributes(myWhite);
  // LogVisAttWorld->SetVisibility(false);
  // fLogicWorld->SetVisAttributes(LogVisAttWorld);

  fPhysicWorld = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), "World", fLogicWorld, NULL, false, 0);

  DSEventHandler::Get()->SetDetectorFlag(fDetectorConfiguration);

  if (fDetectorConfiguration == 0) {  // entire DS50 detector
    DSLog(routine) << " Detector Configuration - TPC+NV+WT: " << fDetectorConfiguration << endlog;
    // Water Tank
    DSDetectorWaterTank* WTVolume = new DSDetectorWaterTank(fPhysicWorld);
    // Neutron Veto
    DSDetectorNeutronVeto* NVVolume = new DSDetectorNeutronVeto(WTVolume->GetDetectorComponent());
    // TPC
    new DSDetectorDS50(NVVolume->GetDetectorComponent());

    // if(fIsSource)  new DSDetectorSource( NVVolume->GetDetectorComponent() ) ;
    if (DSStorage::Get()->GetSourceHolderFlag()) new DSDetectorSource(NVVolume->GetDetectorComponent());

  } else if (fDetectorConfiguration == 1) {  // DS50 TPC + Veto
    DSLog(routine) << " Detector Configuration - TPC+NV: " << fDetectorConfiguration << endlog;
    // Neutron Veto
    DSDetectorNeutronVeto* NVVolume = new DSDetectorNeutronVeto(fPhysicWorld);
    // TPC
    new DSDetectorDS50(NVVolume->GetDetectorComponent());

  } else if (fDetectorConfiguration == 2) {  // DS50 - only TPC
    DSLog(routine) << " Detector Configuration - TPC: " << fDetectorConfiguration << endlog;
    // TPC
    new DSDetectorDS50(fPhysicWorld);

  } else if (fDetectorConfiguration == 3) {  // WT + NV + DS50 TPC + Collimator - 808 to leave a gap with
                                             // respect to "real" detector configurations. This is for
                                             // testing purposes
    DSLog(routine) << " Detector Configuration - TPC+NV+WT+Collimator: " << fDetectorConfiguration << endlog;
    // Water Tank
    DSDetectorWaterTank* WTVolume = new DSDetectorWaterTank(fPhysicWorld);
    // Neutron Veto
    DSDetectorNeutronVeto* NVVolume = new DSDetectorNeutronVeto(WTVolume->GetDetectorComponent());
    // TPC
    new DSDetectorDS50(NVVolume->GetDetectorComponent());
    // Collimator
    new DSDetectorCalibrationDevice(NVVolume->GetDetectorComponent());

  } else if (fDetectorConfiguration == 4) {  // 20k - with GdWater Veto inside CTF

    DSLog(routine) << " Detector Configuration - DS20k: " << fDetectorConfiguration << endlog;

    // recast the rate header value to store the TPC size.
    DSEventHandler::Get()->SetRate(DSStorage::Get()->GetDS20kTPCheight() / cm * 1000 + DSStorage::Get()->GetDS20kTPCedge() / cm);

    DSDetectorDuneCryostat* ProtoDuneCryo = new DSDetectorDuneCryostat(fPhysicWorld);
    DSDetectorDS20kNeutronVeto* VetoVessel = new DSDetectorDS20kNeutronVeto(ProtoDuneCryo->GetDetectorComponent());
    /* DSDetectorDS20k *ds20k =  */
    new DSDetectorDS20k(VetoVessel->GetDetectorComponent());

    // if(DSStorage::Get()->GetGantryConfiguration()){
    //   gantry = new
    //   DSDetectorGantryCalibration(PlasticScintillator->GetDetectorComponent());
    // }

  } else if (fDetectorConfiguration == 5) {  // DS20k tpc only

    DSLog(routine) << " Detector Configuration - DS20k (TPC only): " << fDetectorConfiguration << endlog;
    // recast the rate header value to store the TPC size.
    DSEventHandler::Get()->SetRate(DSStorage::Get()->GetDS20kTPCheight() / cm * 1000 + DSStorage::Get()->GetDS20kTPCedge() / cm);
    /* DSDetectorDS20k *ds20k = */ new DSDetectorDS20k(fPhysicWorld);

  } else if (fDetectorConfiguration == 6) {  // ds Proto
    DSLog(routine) << "Detector configuration - ds Proto  " << fDetectorConfiguration << endlog;
    new DSDetectorProto(fPhysicWorld);

  } else if (fDetectorConfiguration == 7) {  // ds Proto0
    DSLog(routine) << "Detector configuration - ds Proto0 " << fDetectorConfiguration << endlog;
    new DSDetectorProtoProto(fPhysicWorld);

  } else if (fDetectorConfiguration == 8) {  // TPC for ReD setup"
    DSLog(routine) << "Detector configuration - TPC for ReD - 2 " << fDetectorConfiguration << endlog;
    DSStorage::Get()->SetReDGeometry(true);
    new DSDetectorReD(fPhysicWorld);
    // DSDetectorNeutronReD *NeutronDetector = new
    // DSDetectorNeutronReD(fPhysicWorld);
    //  180129 MK:
    //  I removed the neutron detector construction from here to DSDetectorReD.
    //  I wanted to have it optional.
    //  In DSDetectorReD I can check for the specification of
    //  ND_survey_filename. If set, I'll load the LScis, otherwise not. Same for
    //  CT beamline and various scattering chambers and cryostats. We should
    //  control our geometry by flag in the messenger, and then construct only
    //  what is needed.
    //
    //  end 180129 MK
  }

  else if (fDetectorConfiguration == 9) {  // DART geometry
    DSLog(routine) << "Detector configuration (configuration 314)  - DART -  " << fDetectorConfiguration << endlog;
    new DSDetectorDart(fPhysicWorld);
    DSStorage::Get()->SetDartGeometry(true);
  }

  else if (fDetectorConfiguration == 10) {  // ArDM geometry
    DSLog(routine) << "Detector configuration (configuration 421)  - ArDM -  " << fDetectorConfiguration << endlog;
    new DSDetectorArDM(fPhysicWorld);
    DSStorage::Get()->SetArDMGeometry(true);
  }

  else if (fDetectorConfiguration == 11) {  // Licorne (ARIS) geometry
    new DSDetectorLicorne(fPhysicWorld);
    new DSDetectorLicorneNeutron(fPhysicWorld);
  }

  else if (fDetectorConfiguration == 12) {  // ARIS-ER (UC Davis) geometry
    new DSDetectorARISER(fPhysicWorld);
    new DSDetectorARISERGamma(fPhysicWorld);
    DSStorage::Get()->SetAriserGeometry(true);
  }

  DSLog(development) << " LAr - GAr boundary z coordinate set to: " << DSStorage::Get()->GetLArGArBoundaryPosZ() / cm << " cm" << endlog;

  return fPhysicWorld;
}

void DSDetectorConstruction::PrintDetectorParameters() {}

void DSDetectorConstruction::UpdateGeometry() {
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructDetector());
}
