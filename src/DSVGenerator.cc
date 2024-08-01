//---------------------------------------------------------------------------//
/**
 *
 * CLASS DECLARATION:  DSVGenerator.cc
 *
 *---------------------------------------------------------------------------//
 *
 * DESCRIPTION:
 *
 */
// Begin description of class here
/**
 *
 * Pure virtual base class for DS generators.
 *
 */
// End class description
//
/**
 * SPECIAL NOTES:
 *
 */
//
// --------------------------------------------------------------------------//
/**
 * AUTHOR: davide.franco@mi.infn.it
 */
// --------------------------------------------------------------------------//

#include "DSVGenerator.hh"
#include <cstdlib>
#include <fstream>
#include "DSEventHandler.hh"
#include "DSLogger.hh"
#include "DSStorage.hh"
#include "DSVGeneratorMessenger.hh"
#include "G4Material.hh"
#include "G4PhysicalConstants.hh"
#include "G4SPSRandomGenerator.hh"
#include "G4SystemOfUnits.hh"
#include "G4TransportationManager.hh"
#include "Randomize.hh"
//---------------------------------------------------------------------------//

DSVGenerator::DSVGenerator(const G4String& myname) : fG4Messenger(0), fReportingFrequency(1000) {
  fGeneratorName = myname;
  fVolumeFlag = false;
  fEnergyFlag = false;
  fNumberOfHits = 0;
  G4ThreeVector zero(0., 0., 0.);
  fRndGen = new G4SPSRandomGenerator;
  IsMyEnergyDistribution = false;
  fCharge = 0.0;
  fNumberOfParticles = 1;
  gNavigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
  fPositionFlag = 1;
  fVerticalShift = 0 * cm;
  IsG2 = false;
  // SetRadius0(0.);
  // SetCentreCoords(zero);

  IsVolumetric = false;
  fMessenger = new DSVGeneratorMessenger(this);

  fPhysVolumeName = "";
}

//---------------------------------------------------------------------------//

// DSVGenerator::DSVGenerator(const DSVGenerator & other)
//{;}

//---------------------------------------------------------------------------//

DSVGenerator::~DSVGenerator() {
  delete fSPSPos;
  delete fSPSAng;
  delete fSPSEne;
}

G4ThreeVector DSVGenerator::GetVParticleDirection() {

  if (fPositionFlag != 1) {
    //  if(fVolumeFlag) {
    fSPSAng->SetAngDistType("iso");
    fSPSAng->SetPosDistribution(fSPSPos);
    return fSPSAng->GenerateOne();
  }

  return fDirection;
}

G4ThreeVector DSVGenerator::GetVParticlePosition() {

  G4ThreeVector thePos(0., 0., 0.);

  if (fPositionFlag == 1) return fPosition;
  else if (fPositionFlag == 2)
    return fSPSPos->GenerateOne();
  else if (fPositionFlag == 3)
    thePos = GenerateInCryostats();
  else if (fPositionFlag == 4)
    thePos = GenerateInTeflon();
  else if (fPositionFlag == 5)
    thePos = GenerateInFusedSilica();
  else if (fPositionFlag == 6)
    thePos = GenerateInPMTPhotocathode();
  else if (fPositionFlag == 7)
    thePos = GenerateInPMTStem();
  else if (fPositionFlag == 8)
    thePos = GenerateInLiquidArgon();
  else if (fPositionFlag == 9)
    thePos = GenerateInHolderSource();
  else if (fPositionFlag == 10)
    thePos = GenerateInSiPM();
  else if (fPositionFlag == 11)
    thePos = GenerateInGrid();
  else if (fPositionFlag == 12)
    thePos = GenerateInCopperRings();
  else if (fPositionFlag == 13)
    thePos = GenerateInLiquidScintillator();
  else if (fPositionFlag == 14)
    thePos = GenerateInVetoPMTs();
  else if (fPositionFlag == 15)
    thePos = GenerateInVetoSteelSphere();
  else if (fPositionFlag == 16)
    thePos = GenerateInArDMTestArgonCopper();
  else if (fPositionFlag == 17)
    thePos = GenerateInArDMDartDet();
  else if (fPositionFlag == 18)
    thePos = GenerateInArDMDartExt();
  else if (fPositionFlag == 19)
    thePos = GenerateInArDMDartSiPM();
  else if (fPositionFlag == 20)
    thePos = GenerateInArDMTestArgonNS();
  else if (fPositionFlag == 21)
    thePos = GenerateInArDMTestArgonAcrylic();
  else if (fPositionFlag == 22)
    thePos = GenerateInArDMTestArgonReflector();
  else if (fPositionFlag == 23)
    thePos = GenerateInArDMTestArgonWLS();
  else if (fPositionFlag == 24)
    thePos = GenerateInArDMTestArgonUnderground();
  else if (fPositionFlag == 25)
    thePos = GenerateInDS20kSupportStructure();
  else if (fPositionFlag == 26)
    thePos = GenerateIn7mSphere();
  else if (fPositionFlag == 27)
    thePos = GenerateInPlasticScintillator();
  else if (fPositionFlag == 28)
    thePos = GenerateInArgonBufferInside();
  else if (fPositionFlag == 29)
    thePos = GenerateInArgonBufferOutside();
  else if (fPositionFlag == 30)
    thePos = GenerateInSiPM_Veto();
  else if (fPositionFlag == 31)
    thePos = GenerateProtoDUNEFoam();
  else if (fPositionFlag == 32)
    thePos = GenerateProtoDUNEMembrane();
  else if (fPositionFlag == 33)
    thePos = GenerateInGantryPipe();
  else if (fPositionFlag == 34)
    thePos = GenerateInTitaniumVessel();
  else if (fPositionFlag == 35)
    thePos = GenerateInGasPocket();
  else if (fPositionFlag == 36)
    thePos = GenerateInAriserSourceHolder();
  else if (fPositionFlag == 37)
    thePos = GenerateInLateralTPB();
  else if (fPositionFlag == 38)
    thePos = GenerateInLateralESR();
  else if (fPositionFlag == 39)
    thePos = GenerateInLateralAcrylicPanels();
  else if (fPositionFlag == 40)
    thePos = GenerateInPDM();
  else if (fPositionFlag == 41)
    thePos = GenerateInGridSupport();
  else if (fPositionFlag == 42)  
    thePos = GenerateInTPCBarrel();
  else if (fPositionFlag == 43)
    thePos = GenerateInAnode(); 
  else if (fPositionFlag == 44)
    thePos = GenerateInCathode();  
  else if (fPositionFlag == 45)
    thePos = GenerateInOPsAcrylic(); 

  else if (fPositionFlag == 100)
    thePos = GenerateInPhysVolume();
  return thePos;
}

G4double DSVGenerator::GetVParticleEnergy(G4ParticleDefinition* pDef) {
  if (fEnergyFlag) return fSPSEne->GenerateOne(pDef);
  if (IsMyEnergyDistribution) {
    double val = G4UniformRand();

    // add a protection to negative energies
    while (val < fMyProb[0]) val = G4UniformRand();

    for (int i = 1; i < int(fMyProb.size()); ++i) {
      if (fMyProb[i] > val) {
        double x2 = fMyProb[i];
        double x1 = fMyProb[i - 1];
        double y2 = fMyEne[i];
        double y1 = fMyEne[i - 1];
        if (x2 - x1 == 0) continue;
        double _m = (y2 - y1) / (x2 - x1);
        double q = y2 - x2 * _m;
        fEnergy = _m * val + q;
        if (fEnergy < 0)
          cout << " EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE "
               << " " << i << " " << fMyProb[i] << " " << fMyProb[i - 1] << " " << val << endl;
        if (fEnergy < fMyEne[0]) cout << " ZZZZZZZZZZZZZZZZZZZZZZZZZZZz " << endl;
        return _m * val + q;
      }
    }
  }

  return fEnergy;
}

void DSVGenerator::SetIsVolumeDistribution(G4bool val) {
  fVolumeFlag = val;
  if (fVolumeFlag && !IsVolumetric) {
    DSLog(routine) << " Random spatial distribution activated " << endlog;
    fSPSAng = new G4SPSAngDistribution();
    fSPSPos = new G4SPSPosDistribution();
    fSPSPos->SetBiasRndm(fRndGen);
    fSPSAng->SetBiasRndm(fRndGen);
    IsVolumetric = true;
  }
}

void DSVGenerator::SetIsEnergyDistribution(G4bool val) {
  fEnergyFlag = val;
  if (fEnergyFlag) {
    fSPSEne = new G4SPSEneDistribution;
    fSPSEne->SetBiasRndm(fRndGen);
  }
}

void DSVGenerator::SetEnergyFileName(string filename) {
  IsMyEnergyDistribution = true;
  int count = 0;
  double x, y;
  double prob = 0;
  ifstream ff(filename.c_str());
  while (!ff.eof()) {
    ff >> x >> y;
    if (ff.eof()) break;
    prob += y;
    fMyEne.push_back(x * keV);
    fMyProb.push_back(prob);
    count++;
  }
  ff.close();
  for (int i = 0; i < count; ++i) fMyProb[i] /= prob;
}

void DSVGenerator::SetUniformTPC() {
  fPositionFlag = 2;
  SetIsVolumeDistribution(true);
  SetPosDisType("Volume");
  SetPosDisShape("Cylinder");
  SetCentreCoords(G4ThreeVector(0, 0, -36.5 * mm));
  SetRadius(200. * mm);
  SetHalfZ(200. * mm);

  if (DSStorage::Get()->Get5KGeometry()) {
    SetCentreCoords(G4ThreeVector(0, 0, 0 * mm));
    SetRadius(840. * mm);
    SetHalfZ(908. * mm);
  }
  if (DSStorage::Get()->Get20KGeometry()) {
    SetCentreCoords(G4ThreeVector(0, 0, 0 * mm));
    SetRadius(DSStorage::Get()->GetDS20kTPCedge() / 2. * sqrt(4 + 2 * sqrt(2)) + 10 * cm);
    SetHalfZ(DSStorage::Get()->GetDS20kTPCheight() / 2. + 3 * cm);
  }
  if (DSStorage::Get()->GetDartGeometry()) {
    SetCentreCoords(G4ThreeVector(0, 0, 0));
    SetHalfZ(DSStorage::Get()->GetDartTeflonHeight() / 2. + DSStorage::Get()->GetDartTeflonThickness());
    SetRadius(DSStorage::Get()->GetDartTeflonRadius() + DSStorage::Get()->GetDartTeflonThickness());
    /*G4cout << "MSG From DSVGenerator " << " "
           << " Teflon specs Height : " <<
       DSStorage::Get()->GetDartTeflonHeight()/cm << " cm "
           << " Radius : " << DSStorage::Get()->GetDartTeflonRadius()/cm << " cm
       "
           << " Thickness : " << DSStorage::Get()->GetDartTeflonThickness()/cm
       << " cm " << G4endl;  */
  }
  ConfineSourceToVolume("ActiveLAr");
}

void DSVGenerator::SetUniformGasPocket() {
  fPositionFlag = 2;
  SetIsVolumeDistribution(true);
  SetPosDisType("Volume");
  SetPosDisShape("Cylinder");
  SetCentreCoords(G4ThreeVector(0, 0, +150.0 * mm));
  SetRadius(200. * mm);
  SetHalfZ(6. * mm);
  if (DSStorage::Get()->Get5KGeometry()) {
    SetCentreCoords(G4ThreeVector(0, 0, 905 * mm));
    SetRadius(840. * mm);
    SetHalfZ(10. * mm);
  }
  ConfineSourceToVolume("GasPocket");
}

void DSVGenerator::SetTPCCenter() {
  SetVParticlePosition(G4ThreeVector(0, 0, -36.54 * mm));
  if (fVolumeFlag) SetCentreCoords(G4ThreeVector(0, 0, -36.54 * mm));
  if (DSStorage::Get()->Get5KGeometry()) SetCentreCoords(G4ThreeVector(0, 0, 0));
}

G4ThreeVector DSVGenerator::GenerateInCryostats() {

  G4double center_z = IsG2 ? 0. * mm : 161.5 * mm;
  G4double height = IsG2 ? 2231 * mm : 1400. * mm;
  G4double radius = IsG2 ? 910. * mm : 356. * mm;

  if (DSStorage::Get()->Get5KGeometry()) {

    center_z = 130. * mm;
    height = 2774 * mm;
    radius = 996 * mm;

  } else if (DSStorage::Get()->Get20KGeometry()) {

    // center_z = 25.*mm;
    // height   = 8400*mm + ( DSStorage::Get()->GetDS20kTPCheight() - 240*cm);
    // radius   = ( DSStorage::Get()-> GetDS20kTPCedge()/2. *sqrt(4  + 2
    // *sqrt(2))  +  DSStorage::Get()->GetDS20kCryoCornerDistance() + 50*cm);

    center_z = 0. * cm;
    height = 7. * m;
    radius = 2.7 * m;

  } else if (DSStorage::Get()->GetDartGeometry()) {

    center_z = 0. * cm;
    radius = DSStorage::Get()->GetDartDewarRadius() + DSStorage::Get()->GetDartDewarThickness();
    height = DSStorage::Get()->GetDartDewarHeight() + DSStorage::Get()->GetDartDewarThickness();

  } else if (DSStorage::Get()->GetArDMGeometry()) {

    center_z = 0 * mm;
    height = DSStorage::Get()->GetArDMTotalTankHeight();
    radius = DSStorage::Get()->GetArDMTankRadius();
  }

  SetIsVolumeDistribution(true);
  SetPosDisType("Volume");
  SetPosDisShape("Cylinder");
  SetCentreCoords(G4ThreeVector(0, 0, center_z));
  SetHalfZ(height / 2.);
  SetRadius(radius);

  string mymat = "StainlessSteel";
  if (DSStorage::Get()->Get20KGeometry()) {
    if (DSStorage::Get()->GetDS20kCryoMaterial() == 0) mymat = "StainlessSteel";
    else if (DSStorage::Get()->GetDS20kCryoMaterial() == 1)
      mymat = "MetalTitanium";
    else if (DSStorage::Get()->GetDS20kCryoMaterial() == 2)
      mymat = "MetalCopperCryo";
  }

  G4ThreeVector myPos(0., 0., 0.);
  myPos = fSPSPos->GenerateOne();
  while (!CheckMaterial(myPos, mymat.c_str())) myPos = fSPSPos->GenerateOne();

  return myPos;
}

G4ThreeVector DSVGenerator::GenerateInLiquidArgon() {

  string mymat = "LiquidArgon";

  G4double center_z = IsG2 ? 0. * mm : 161.5 * mm;
  G4double height = IsG2 ? 2231 * mm : 1400. * mm;
  G4double radius = IsG2 ? 910. * mm : 356. * mm;

  if (DSStorage::Get()->Get5KGeometry()) {
    center_z = 0 * mm;
    height = 1820 * mm;
    radius = 990 * mm;

  } else if (DSStorage::Get()->Get20KGeometry()) {
    center_z = 0 * mm;
    radius = DSStorage::Get()->GetDS20kTPCedge() / 2. * sqrt(4 + 2 * sqrt(2)) + 10 * cm;
    height = DSStorage::Get()->GetDS20kTPCheight() + 10 * cm;

  } else if (DSStorage::Get()->GetArDMGeometry()) {
    // vpesudo
    center_z = 0 * mm;
    height = DSStorage::Get()->GetArDMTotalTankHeight();
    radius = DSStorage::Get()->GetArDMTankRadius();
    mymat = "PseudoArgon";
  }

  SetIsVolumeDistribution(true);
  SetPosDisType("Volume");
  SetPosDisShape("Cylinder");
  SetCentreCoords(G4ThreeVector(0, 0, center_z));
  SetHalfZ(height / 2.);
  SetRadius(radius);

  G4ThreeVector myPos(0., 0., 0.);
  myPos = fSPSPos->GenerateOne();

  if (DSStorage::Get()->GetArDMGeometry()) {
    while (!(CheckMaterial(myPos, mymat.c_str()))) myPos = fSPSPos->GenerateOne();
  } else {
    while (!(CheckMaterial(myPos, mymat.c_str()) || CheckMaterial(myPos, "NSLiquidArgon"))) myPos = fSPSPos->GenerateOne();
  }

  // std::cout<<myPos.x()<<" "<<myPos.y()<<" "<<myPos.z()<<std::endl;

  return myPos;
}
G4ThreeVector DSVGenerator::GenerateInGasPocket() {

  string mymat = "GaseousArgon";


  // DS20k
  float radius = DSStorage::Get()->GetDS20kTPCedge() / 2. * sqrt(4 + 2 * sqrt(2)) + 10 * cm;
  float height = DSStorage::Get()->GetDS20kTPCheight() /2. ;

  SetIsVolumeDistribution(true);
  SetPosDisType("Volume");
  SetPosDisShape("Cylinder");
  SetCentreCoords(G4ThreeVector(0, 0, height));
  SetHalfZ(2.5*cm/2.);
  SetRadius(radius);

  G4ThreeVector myPos(0., 0., 0.);
  myPos = fSPSPos->GenerateOne();

  while (!(CheckMaterial(myPos, mymat.c_str()))) myPos = fSPSPos->GenerateOne();

  return myPos;
}


G4ThreeVector DSVGenerator::GenerateInLateralTPB() {

  SetIsVolumeDistribution(true);
  SetPosDisType("Volume");
  SetPosDisShape("Para");

  SetCentreCoords(G4ThreeVector(176* cm, 0, 0. ));
  SetHalfZ(175.25 * cm);
  SetHalfX(1 * cm);
  SetHalfY(75 * cm);

  G4ThreeVector myPos(0., 0., 0.);
  myPos = fSPSPos->GenerateOne();
  while (!CheckMaterial(myPos, "TPB") ) myPos = fSPSPos->GenerateOne();

  int val = G4UniformRand() * 8  ;
  double angle = val * M_PI/4. ;
  double oldX = myPos.getX() ;
  double oldY = myPos.getY() ;
  myPos.setX ( oldX* cos ( angle ) -  oldY * sin ( angle ) )   ;
  myPos.setY ( oldX* sin ( angle ) + oldY * cos (angle ) )   ;

  return myPos;

}

G4ThreeVector DSVGenerator::GenerateInLateralESR() {

//sample from a small volume and then randomly rotate on other faces
  SetIsVolumeDistribution(true);
  SetPosDisType("Volume");
  SetPosDisShape("Para");

  SetCentreCoords(G4ThreeVector(176* cm, 0, 0. ));
  SetHalfZ(175.25 * cm);
  SetHalfX(2 * cm);
  SetHalfY(75 * cm);

  G4ThreeVector myPos(0., 0., 0.);
  myPos = fSPSPos->GenerateOne();
  while (!CheckMaterial(myPos, "ESR") ) myPos = fSPSPos->GenerateOne();

  int val = G4UniformRand() * 8  ;
  double angle = val * M_PI/4. ;
  double oldX = myPos.getX() ;
  double oldY = myPos.getY() ;
  myPos.setX ( oldX* cos ( angle ) -  oldY * sin ( angle ) )   ;
  myPos.setY ( oldX* sin ( angle ) + oldY * cos (angle ) )   ;

  return myPos;

}

G4ThreeVector DSVGenerator::GenerateInLateralAcrylicPanels() {

  SetIsVolumeDistribution(true);
  SetPosDisType("Volume");
  SetPosDisShape("Para");

  SetCentreCoords(G4ThreeVector(176* cm, 0, 0. ));
  SetHalfZ(175.25 * cm);
  SetHalfX(2 * cm);
  SetHalfY(75 * cm);

  G4ThreeVector myPos(0., 0., 0.);
  myPos = fSPSPos->GenerateOne();
  while (!CheckMaterial(myPos, "ReflectorAcrylic") ) myPos = fSPSPos->GenerateOne();

  int val = G4UniformRand() * 8  ;
  double angle = val * M_PI/4. ;
  double oldX = myPos.getX() ;
  double oldY = myPos.getY() ;
  myPos.setX ( oldX* cos ( angle ) -  oldY * sin ( angle ) )   ;
  myPos.setY ( oldX* sin ( angle ) + oldY * cos (angle ) )   ;

  return myPos;

}
G4ThreeVector DSVGenerator::GenerateInTeflon() {

  SetIsVolumeDistribution(true);
  SetPosDisType("Volume");
  SetPosDisShape("Cylinder");
  SetCentreCoords(G4ThreeVector(0, 0, 160. * mm));
  SetHalfZ(860. / 2. * mm);
  SetRadius(235. * mm);
  if (DSStorage::Get()->Get5KGeometry()) {
    SetCentreCoords(G4ThreeVector(0, 0, 0));
    SetHalfZ(905 * mm);
    SetRadius(867. * mm);
  }
  if (DSStorage::Get()->Get20KGeometry()) {
    SetCentreCoords(G4ThreeVector(0, 0, 0));
    SetHalfZ(DSStorage::Get()->GetDS20kTPCheight() / 2 + DSStorage::Get()->GetDS20kTPCVesselThickness());
    SetRadius(DSStorage::Get()->GetDS20kTPCedge() / 2. * sqrt(4 + 2 * sqrt(2)) + DSStorage::Get()->GetDS20kTPCVesselThickness());
  }
  if (DSStorage::Get()->GetDartGeometry()) {
    SetCentreCoords(G4ThreeVector(0, 0, 0));
    SetHalfZ(DSStorage::Get()->GetDartTeflonHeight() / 2. + DSStorage::Get()->GetDartTeflonThickness());
    SetRadius(DSStorage::Get()->GetDartTeflonRadius() + DSStorage::Get()->GetDartTeflonThickness());
    /*G4cout << "MSG From DSVGenerator " << " "
           << " Teflon specs Height : " <<
       DSStorage::Get()->GetDartTeflonHeight()/cm << " cm "
           << " Radius : " << DSStorage::Get()->GetDartTeflonRadius()/cm << " cm
       "
           << " Thickness : " << DSStorage::Get()->GetDartTeflonThickness()/cm
       << " cm " << G4endl;  */
  }
  if (DSStorage::Get()->GetArDMGeometry()) {
    SetCentreCoords(G4ThreeVector(0, 0, 0));
    SetHalfZ(DSStorage::Get()->GetArDMTeflonHeight() / 2);
    SetRadius(DSStorage::Get()->GetArDMTeflonRadius());
  }

  G4ThreeVector myPos(0., 0., 0.);
  myPos = fSPSPos->GenerateOne();
  while (!(CheckMaterial(myPos, "Teflon") || CheckMaterial(myPos, "Acrylic"))) myPos = fSPSPos->GenerateOne();

  return myPos;
}

G4ThreeVector DSVGenerator::GenerateInArDMTestArgonCopper() {
  G4double center_z = 0.;
  G4double height = 0.;
  G4double radius = 0.;

  if (DSStorage::Get()->GetArDMGeometry()) {
    center_z = 0. * cm;
    height = DSStorage::Get()->GetArDMTestArgonHeight();
    radius = DSStorage::Get()->GetArDMTestArgonRadius();
  }

  SetIsVolumeDistribution(true);
  SetPosDisType("Volume");
  SetPosDisShape("Cylinder");
  SetCentreCoords(G4ThreeVector(0, 0, center_z));
  SetHalfZ(height / 2.);
  SetRadius(radius);

  string mymat = "MetalCopper";

  G4ThreeVector myPos(0., 0., 0.);
  myPos = fSPSPos->GenerateOne();
  while (!CheckMaterial(myPos, mymat.c_str())) myPos = fSPSPos->GenerateOne();

  return myPos;
}

G4ThreeVector DSVGenerator::GenerateInArDMTestArgonAcrylic() {
  G4double center_z = 0.;
  G4double height = 0.;
  G4double radius = 0.;

  if (DSStorage::Get()->GetArDMGeometry()) {
    center_z = 0. * cm;
    height = DSStorage::Get()->GetArDMTestArgonHeight();
    radius = DSStorage::Get()->GetArDMTestArgonRadius();
  }

  SetIsVolumeDistribution(true);
  SetPosDisType("Volume");
  SetPosDisShape("Cylinder");
  SetCentreCoords(G4ThreeVector(0, 0, center_z));
  SetHalfZ(height / 2.);
  SetRadius(radius);

  string mymat = "Acrylic";

  G4ThreeVector myPos(0., 0., 0.);
  myPos = fSPSPos->GenerateOne();
  while (!CheckMaterial(myPos, mymat.c_str())) myPos = fSPSPos->GenerateOne();

  return myPos;
}

G4ThreeVector DSVGenerator::GenerateInArDMTestArgonUnderground() {
  G4double center_z = 0.;
  G4double height = 0.;
  G4double radius = 0.;

  if (DSStorage::Get()->GetArDMGeometry()) {
    center_z = 0. * cm;
    height = DSStorage::Get()->GetArDMTestArgonHeight();
    radius = DSStorage::Get()->GetArDMTestArgonRadius();
  }

  SetIsVolumeDistribution(true);
  SetPosDisType("Volume");
  SetPosDisShape("Cylinder");
  SetCentreCoords(G4ThreeVector(0, 0, center_z));
  SetHalfZ(height / 2.);
  SetRadius(radius);

  string mymat = "LiquidArgon";

  G4ThreeVector myPos(0., 0., 0.);
  myPos = fSPSPos->GenerateOne();
  while (!CheckMaterial(myPos, mymat.c_str())) myPos = fSPSPos->GenerateOne();

  return myPos;
}

G4ThreeVector DSVGenerator::GenerateInArDMTestArgonNS() {
  G4double center_z = 0.;
  G4double height = 0.;
  G4double radius = 0.;

  if (DSStorage::Get()->GetArDMGeometry()) {
    center_z = 0. * cm;
    height = DSStorage::Get()->GetArDMTestArgonHeight();
    radius = DSStorage::Get()->GetArDMTestArgonRadius();
  }

  SetIsVolumeDistribution(true);
  SetPosDisType("Volume");
  SetPosDisShape("Cylinder");
  SetCentreCoords(G4ThreeVector(0, 0, center_z));
  SetHalfZ(height / 2.);
  SetRadius(radius);

  string mymat = "NSLiquidArgon";

  G4ThreeVector myPos(0., 0., 0.);
  myPos = fSPSPos->GenerateOne();
  while (!CheckMaterial(myPos, mymat.c_str())) myPos = fSPSPos->GenerateOne();

  return myPos;
}

G4ThreeVector DSVGenerator::GenerateInArDMTestArgonReflector() {
  G4double center_z = 0.;
  G4double height = 0.;
  G4double radius = 0.;

  if (DSStorage::Get()->GetArDMGeometry()) {
    center_z = 0. * cm;
    height = DSStorage::Get()->GetArDMTestArgonHeight();
    radius = DSStorage::Get()->GetArDMTestArgonRadius();
  }

  SetIsVolumeDistribution(true);
  SetPosDisType("Volume");
  SetPosDisShape("Cylinder");
  SetCentreCoords(G4ThreeVector(0, 0, center_z));
  SetHalfZ(height / 2.);
  SetRadius(radius);

  string mymat = "Aluminum";

  G4ThreeVector myPos(0., 0., 0.);
  myPos = fSPSPos->GenerateOne();
  while (!CheckMaterial(myPos, mymat.c_str())) myPos = fSPSPos->GenerateOne();

  return myPos;
}

G4ThreeVector DSVGenerator::GenerateInArDMTestArgonWLS() {
  G4double center_z = 0.;
  G4double height = 0.;
  G4double radius = 0.;

  if (DSStorage::Get()->GetArDMGeometry()) {
    center_z = 0. * cm;
    height = DSStorage::Get()->GetArDMTestArgonHeight();
    radius = DSStorage::Get()->GetArDMTestArgonRadius();
  }

  SetIsVolumeDistribution(true);
  SetPosDisType("Volume");
  SetPosDisShape("Cylinder");
  SetCentreCoords(G4ThreeVector(0, 0, center_z));
  SetHalfZ(height / 2.);
  SetRadius(radius);

  string mymat = "TPB";

  G4ThreeVector myPos(0., 0., 0.);
  myPos = fSPSPos->GenerateOne();
  while (!CheckMaterial(myPos, mymat.c_str())) myPos = fSPSPos->GenerateOne();

  return myPos;
}

G4ThreeVector DSVGenerator::GenerateInArDMDartDet() {
  G4double center_z = 0.;
  G4double height = 0.;
  G4double radius = 0.;

  if (DSStorage::Get()->GetArDMGeometry()) {
    center_z = 0. * cm;
    height = DSStorage::Get()->GetArDMTestArgonHeight();
    radius = DSStorage::Get()->GetArDMTestArgonRadius();
  }

  SetIsVolumeDistribution(true);
  SetPosDisType("Volume");
  SetPosDisShape("Cylinder");
  SetCentreCoords(G4ThreeVector(0, 0, center_z));
  SetHalfZ(height / 2.);
  SetRadius(radius);

  //  vpesudo
  string mymat = "LiquidArgon";

  G4ThreeVector myPos(0., 0., 0.);
  myPos = fSPSPos->GenerateOne();
  while (!CheckMaterial(myPos, mymat.c_str())) myPos = fSPSPos->GenerateOne();

  return myPos;
}

G4ThreeVector DSVGenerator::GenerateInArDMDartSiPM() {
  G4double center_z = 0.;
  G4double height = 0.;
  G4double radius = 0.;

  if (DSStorage::Get()->GetArDMGeometry()) {
    center_z = 0. * cm;
    height = DSStorage::Get()->GetArDMTestArgonHeight();
    radius = DSStorage::Get()->GetArDMTestArgonRadius();
  }

  SetIsVolumeDistribution(true);
  SetPosDisType("Volume");
  SetPosDisShape("Cylinder");
  SetCentreCoords(G4ThreeVector(0, 0, center_z));
  SetHalfZ(height / 2.);
  SetRadius(radius);

  // double val = G4UniformRand();
  // if ( val<0.5 ) SetCentreCoords( G4ThreeVector(0,0,-center_z) );
  // else           SetCentreCoords( G4ThreeVector(0,0, center_z) );

  string mymat = "MetalSilicon";

  G4ThreeVector myPos(0., 0., 0.);
  myPos = fSPSPos->GenerateOne();
  while (!CheckMaterial(myPos, mymat.c_str())) myPos = fSPSPos->GenerateOne();

  return myPos;
}

void DSVGenerator::SphereToPrism(G4double faces, G4double hp, G4double rp,  // input: prism parameters
                                 G4double theta, G4double phi,              // input: spherical coordinates
                                 G4double& xp, G4double& yp,
                                 G4double& zp) {  // output: coordinates on the prism surface
                                                  //
                                                  // SphereToPrism.C      (pablo.garcia@ciemat.es, 11 may 2018)
                                                  //
  // This function projects a sphere onto the surface of a regular prism,
  // returning the spatial coordinates (x, y, z) of the projection for the two
  // spherical angles (theta, phi). Conventional spherical coordinates are used:
  //     - faces: number of faces of the prism
  //     - hp:    height of the prism
  //     - rp:    radius of the prism (center of the lateral face)
  //     - theta: polar angle from 0 to pi
  //     - phi:   azimuthal angle from 0 to 2*pi
  //
  // The prism position and orientation convention is as follows:
  //     - (0,0,0): geometrical centre of the prism ("centre of mass")
  //     - z axis:  from (0,0,0) to the centre of the upper polygonal base
  //     - x axis:  from (0,0,0) to the center of the parallelogram (side face)
  //     - y axis:  orthogonal to the XZ plane
  //

  if (faces < 3) {  // not a prism
    xp = 0.0;
    yp = 0.0;
    zp = 0.0;
    return;
  }

  // one sector of the prism (radians)
  G4double alpha = 360. / faces *  ( M_PI / 180.);

  // convert phi range from [-pi,pi] to [0,2pi]
  phi += M_PI ;

  // coordinates in the prism sides
  G4double xp1, yp1, l1;
  G4double xbar, ybar, b;

  G4int sector = (G4int)((phi + alpha / 2.) / alpha);
  G4double phibar = phi - alpha * sector;

  xbar = rp;
  ybar = rp * tan(phibar);
  b = phi - phibar;

  xp1 = xbar * cos(b) - ybar * sin(b);
  yp1 = xbar * sin(b) + ybar * cos(b);

  l1 = sqrt(xp1 * xp1 + yp1 * yp1);

  // coordinates in the prism base
  G4double rz, lz;

  G4double ang  = M_PI/2. - theta;
  G4double sign = (ang > 0) - (ang < 0) ;
  zp = hp / 2. * sign ;
  rz = zp / cos(theta);
  xp = rz * sin(theta) * cos(phi);
  yp = rz * sin(theta) * sin(phi);

  lz = sqrt(xp * xp + yp * yp);

  // select base or side
  if (lz > l1) {
    xp = xp1;
    yp = yp1;
    zp = l1 / tan(theta);
  }

  DSLog(debugging) << "SphereToPrism: " << xp << " " << yp << " " << zp << endlog;

  return;
}

G4ThreeVector DSVGenerator::GenerateInArDMDartExt() {

  // This generator projects a point in a sphere (theta, phi) onto
  // an octogonal prism, returning the x, y, z, coordinates on the
  // surface of the prism. The code SphereToPrism can be reused
  // for any n-sided prism just modifying the 'faces' parameter.

  G4ThreeVector myPos(0., 0., 0.);

  if (!DSStorage::Get()->GetArDMGeometry()) return myPos;

  G4double xp, yp, zp;      // coordinates on the surface of the prism
  G4double height, radius;  // number of faces, height and radius of the prism

  height = DSStorage::Get()->GetArDMShieldHeight();
  radius = DSStorage::Get()->GetArDMShieldRadius();

  SetIsVolumeDistribution(true);
  SetPosDisType("Volume");
  SetPosDisShape("Sphere");
  SetCentreCoords(G4ThreeVector(0, 0, 0));
  SetRadius(1 * cm);

  // string mymat = "HDPE" ;

  myPos = fSPSPos->GenerateOne();  // point in a spheric volume

  SphereToPrism(8, height, radius, myPos.getTheta(), myPos.getPhi(), xp, yp, zp);

  myPos.set(xp, yp, zp);

  return myPos;
}

G4ThreeVector DSVGenerator::GenerateInCopperRings() {

  SetIsVolumeDistribution(true);
  SetPosDisType("Volume");
  SetPosDisShape("Cylinder");
  /*  SetCentreCoords(G4ThreeVector(0,0,160.*mm));
  SetHalfZ(860./2.*mm);
  SetRadius(235.*mm);
  if (DSStorage::Get()->Get5KGeometry()) {
    SetCentreCoords(G4ThreeVector(0,0,0));
    SetHalfZ(905*mm);
    SetRadius(867.*mm);
  }
*/
  if (DSStorage::Get()->Get20KGeometry()) {
    SetCentreCoords(G4ThreeVector(0, 0, 0));
    SetRadius(DSStorage::Get()->GetDS20kTPCedge() / 2. * sqrt(4 + 2 * sqrt(2)) + 10 * cm);
    SetHalfZ(DSStorage::Get()->GetDS20kTPCheight() / 2. + 3 * cm);
  }

  G4ThreeVector myPos(0., 0., 0.);
  myPos = fSPSPos->GenerateOne();
  while (!CheckMaterial(myPos, "MetalCopper")) myPos = fSPSPos->GenerateOne();

  return myPos;
}

G4ThreeVector DSVGenerator::GenerateInGrid() {

  SetIsVolumeDistribution(true);
  SetPosDisType("Volume");
  SetPosDisShape("Cylinder");
  /*  SetCentreCoords(G4ThreeVector(0,0,160.*mm));
  SetHalfZ(860./2.*mm);
  SetRadius(235.*mm);
  if (DSStorage::Get()->Get5KGeometry()) {
    SetCentreCoords(G4ThreeVector(0,0,0));
    SetHalfZ(905*mm);
    SetRadius(867.*mm);
  }
*/
  if (DSStorage::Get()->Get20KGeometry()) {
    SetCentreCoords(G4ThreeVector(0, 0, DSStorage::Get()->GetDS20kTPCheight() / 2. - 0.3 * cm));
    SetRadius(DSStorage::Get()->GetDS20kTPCedge() / 2. * sqrt(4 + 2 * sqrt(2)) + 10 * cm);
    SetHalfZ(10. * mm);
  }

  G4ThreeVector myPos(0., 0., 0.);
  myPos = fSPSPos->GenerateOne();
  while (!CheckMaterial(myPos, "GridSteel")) myPos = fSPSPos->GenerateOne();

  return myPos;
}

G4ThreeVector DSVGenerator::GenerateInGridSupport() {

  SetIsVolumeDistribution(true);
  SetPosDisType("Volume");
  SetPosDisShape("Cylinder");

  G4double radius = DSStorage::Get()->GetDS20kTPCedge() / 2. * (1 + sqrt(2)) + DSStorage::Get()->GetDS20kTPCVesselThickness() + DSStorage::Get()->GetDS20kTPBThickness() + 4.4 * cm;
  radius *= 1.1; // to get outer radius of Gd-Acrylic
  
  SetCentreCoords(G4ThreeVector(0, 0, DSStorage::Get()->GetDS20kTPCheight() / 2. - 1.5 * cm));
  SetRadius(radius);
  SetHalfZ(6. * cm);
  

  G4ThreeVector myPos(0., 0., 0.);
  myPos = fSPSPos->GenerateOne();
  while (!CheckMaterial(myPos, "StainlessSteel")) myPos = fSPSPos->GenerateOne();
  
  return myPos;
}

G4ThreeVector DSVGenerator::GenerateInPDM() {

  G4double pos_z;
  G4double half_z;
  G4double radius;

  pos_z = -DSStorage::Get()->GetDS20kTPCheight() / 2. - 5 * cm - DSStorage::Get()->GetSiPMOffset() - 3. * cm;
  half_z = 50. * mm;
  radius = (DSStorage::Get()->GetDS20kTPCedge()) / (sqrt(2 - sqrt(2))) + 5. * cm;

  if (G4UniformRand() > 0.5) {
    pos_z = DSStorage::Get()->GetDS20kTPCheight() / 2. + 5. * cm + 0.7 * cm + DSStorage::Get()->GetSiPMOffset();
  }
  SetIsVolumeDistribution(true);
  SetPosDisType("Volume");
  SetPosDisShape("Cylinder");

  SetCentreCoords(G4ThreeVector(0, 0, pos_z));
  SetHalfZ(half_z);
  SetRadius(radius);

  G4ThreeVector myPos(0., 0., 0.);
  myPos = fSPSPos->GenerateOne();
  while (!(CheckMaterial(myPos, "ArlonPrepreg") || CheckMaterial(myPos, "MetalCopper"))) myPos = fSPSPos->GenerateOne();

  return myPos;
}

G4ThreeVector DSVGenerator::GenerateInSiPM() {
  //  G4double pos_z  = 160.5 * mm ;
  //  G4double half_z = 10./2. * mm ;
  //  G4double radius = 213.*mm ;

  G4double pos_z;
  G4double half_z;
  G4double radius;

  /*
    if (DSStorage::Get()->Get5KGeometry()) {
    pos_z  = 913*mm ;
    half_z = 13. * mm ;
    radius = 840.*mm ;
    }
  */
  // int sign = -1;
  //  if (DSStorage::Get()->Get20KGeometry()) {
  pos_z = -DSStorage::Get()->GetDS20kTPCheight() / 2. - 5 * cm - DSStorage::Get()->GetSiPMOffset() - 3. * cm;
  half_z = 50. * mm;
  // radius = DSStorage::Get()->GetDS20kTPCedge()/2. *sqrt(4  + 2 *sqrt(2))  +
  // 10*cm ;
  radius = (DSStorage::Get()->GetDS20kTPCedge()) / (sqrt(2 - sqrt(2))) + 1. * cm;
  // }

  if (G4UniformRand() > 0.5) {
    /*
          if (DSStorage::Get()->Get5KGeometry()) {
          pos_z  = -903*mm ;
          } else
    */
    //sign = 1;
    //  if (DSStorage::Get()->Get20KGeometry()) {
    pos_z = DSStorage::Get()->GetDS20kTPCheight() / 2. + 5. * cm + 0.7 * cm + DSStorage::Get()->GetSiPMOffset();
    // } else {
    //  pos_z = -221. * mm;
    //  half_z = 10./2. * mm;
    // }
  }
  SetIsVolumeDistribution(true);
  SetPosDisType("Volume");
  SetPosDisShape("Cylinder");

  SetCentreCoords(G4ThreeVector(0, 0, pos_z));
  SetHalfZ(half_z);
  SetRadius(radius);

  G4ThreeVector myPos(0., 0., 0.);
  myPos = fSPSPos->GenerateOne();
  while (!CheckMaterial(myPos, "MetalSilicon")) myPos = fSPSPos->GenerateOne();

  //myPos.setZ(myPos.getZ() + sign * offset);

  return myPos;
}

G4ThreeVector DSVGenerator::GenerateInSiPM_Veto() {

  vector<G4ThreeVector> posVector = DSStorage::Get()->GetDS20kSiPMPosVector();

  SetIsVolumeDistribution(true);
  SetPosDisType("Volume");
  SetPosDisShape("Sphere");

  G4double random = round((G4UniformRand()) * (posVector.size() - 1));

  SetCentreCoords(posVector.at(random));
  SetRadius((DSStorage::Get()->GetDS20kSiPMSide()) * 0.5 * sqrt(2) + 0.2 * cm);

  G4ThreeVector myPos(0., 0., 0.);

  myPos = fSPSPos->GenerateOne();
  while (!CheckMaterial(myPos, "Teflon")) myPos = fSPSPos->GenerateOne();

  return myPos;
}

G4ThreeVector DSVGenerator::GenerateInFusedSilica() {

  G4double pos_z = -35.5 * mm;
  G4double half_z = 22 * cm;
  G4double radius = 223. * mm;
  if (DSStorage::Get()->Get5KGeometry()) {
    pos_z = 913 * mm;
    half_z = 13. * mm;
    radius = 840. * mm;
  }
  if (DSStorage::Get()->Get20KGeometry()) {
    pos_z = -DSStorage::Get()->GetDS20kTPCheight() / 2. - 1 * cm;
    half_z = 60. * mm;
    radius = DSStorage::Get()->GetDS20kTPCedge() / 2. * sqrt(4 + 2 * sqrt(2)) + 10 * cm;
  }

  if (G4UniformRand() > 1.5 / 4.) {
    // in he current design,  there are three acrylic layer (1.5 at the bottom,
    // and 1.5 + 1 cm at the top)
    if (DSStorage::Get()->Get5KGeometry()) {
      pos_z = -903 * mm;
    } else if (DSStorage::Get()->Get20KGeometry()) {
      pos_z = DSStorage::Get()->GetDS20kTPCheight() / 2. + 5 * cm;
    }
  }

  SetIsVolumeDistribution(true);
  SetPosDisType("Volume");
  SetPosDisShape("Cylinder");

  SetCentreCoords(G4ThreeVector(0, 0, pos_z));
  SetHalfZ(half_z);
  SetRadius(radius);

  G4ThreeVector myPos(0., 0., 0.);
  myPos = fSPSPos->GenerateOne();

  string thematerial = "FusedSilica";
  if (DSStorage::Get()->Get20KGeometry()) thematerial = "Acrylic";

  while (!CheckMaterial(myPos, thematerial.c_str())) myPos = fSPSPos->GenerateOne();

  return myPos;
}

G4ThreeVector DSVGenerator::GenerateInPMTPhotocathode() {

  SetIsVolumeDistribution(true);
  SetPosDisType("Volume");

  string thematerial = "FakePhotocathode";

  if (DSStorage::Get()->Get20KGeometry()) {

    SetPosDisShape("Sphere");
    SetCentreCoords(G4ThreeVector(0, 0, 0));
    SetRadius(DSStorage::Get()->GetIDiameter() / 2.);
    SetRadius0(DSStorage::Get()->GetIDiameter() / 2. - 40 * cm);

    thematerial = "Bialkali";

  } else if (DSStorage::Get()->GetArDMGeometry()) {

    SetPosDisShape("Cylinder");
    SetRadius(45. * cm);
    SetHalfZ(15. * cm);

    double val = G4UniformRand();
    if (val < 0.5) SetCentreCoords(G4ThreeVector(0, 0, -75. * cm));
    else
      SetCentreCoords(G4ThreeVector(0, 0, 75. * cm));

  } else {

    // default: DS50
    G4double pos_z = 170.5 * mm;
    G4double half_z = 12. / 2. * mm;
    if (G4UniformRand() > 0.5) {
      pos_z = -220.5 * mm;
      half_z = 12. / 2. * mm;
    }

    SetPosDisShape("Cylinder");
    SetCentreCoords(G4ThreeVector(0, 0, pos_z));
    SetHalfZ(half_z);
    SetRadius(235. * mm);
  }

  G4ThreeVector myPos(0., 0., 0.);
  myPos = fSPSPos->GenerateOne();
  while (!CheckMaterial(myPos, thematerial.c_str())) myPos = fSPSPos->GenerateOne();

  return myPos;
}

G4ThreeVector DSVGenerator::GenerateInPMTStem() {

  G4double pos_z = 287.6 * mm;
  G4double half_z = 1. * mm;
  if (G4UniformRand() > 0.5) {
    pos_z = -344.86 * mm;
    half_z = 1. * mm;
  }

  SetIsVolumeDistribution(true);
  SetPosDisType("Volume");
  SetPosDisShape("Cylinder");
  SetCentreCoords(G4ThreeVector(0, 0, pos_z));
  SetHalfZ(half_z);
  SetRadius(235. * mm);

  G4ThreeVector myPos(0., 0., 0.);
  myPos = fSPSPos->GenerateOne();
  while (!CheckMaterial(myPos, "BoroSilicateGlass")) myPos = fSPSPos->GenerateOne();

  return myPos;
}

G4ThreeVector DSVGenerator::GenerateInHolderSource() {

  G4ThreeVector center = DSStorage::Get()->GetSourceHolderCenter();
  SetIsVolumeDistribution(true);
  SetPosDisType("Volume");
  SetPosDisShape("Sphere");
  SetCentreCoords(center);
  SetRadius(40. * mm);

  G4ThreeVector myPos(0., 0., 0.);
  myPos = fSPSPos->GenerateOne();
  //  while( !CheckMaterial( myPos, "Teflon" ) ) myPos = fSPSPos->GenerateOne();
  while (!CheckMaterial(myPos, "Air")) myPos = fSPSPos->GenerateOne();

  return myPos;
}

G4ThreeVector DSVGenerator::GenerateInLiquidScintillator() {

  static const G4double radius = DSStorage::Get()->GetIDiameter() / 2 + 10 * cm;

  SetIsVolumeDistribution(true);
  SetPosDisType("Volume");
  SetPosDisShape("Sphere");
  SetCentreCoords(G4ThreeVector(0., 0., 0.));
  SetRadius(radius);

  G4ThreeVector myPos(0., 0., 0.);
  myPos = fSPSPos->GenerateOne();
  while (!CheckMaterialIndex(myPos, DSEventHandler::Get()->GetScintillatorIndex())) myPos = fSPSPos->GenerateOne();

  return myPos;
}

G4ThreeVector DSVGenerator::GenerateInPlasticScintillator() {

  static const G4double radius = 4 * m;

  SetIsVolumeDistribution(true);
  SetPosDisType("Volume");
  SetPosDisShape("Sphere");
  SetCentreCoords(G4ThreeVector(0., 0., 0.));
  SetRadius(radius);

  G4ThreeVector myPos(0., 0., 0.);
  myPos = fSPSPos->GenerateOne();
  while (!(CheckMaterial(myPos, "GdAcrylic") || CheckMaterial(myPos, "VetoPScint"))) myPos = fSPSPos->GenerateOne();

  return myPos;
}

G4ThreeVector DSVGenerator::GenerateInTPCBarrel() {

  static const G4double radius = 4 * m;

  SetIsVolumeDistribution(true);
  SetPosDisType("Volume");
  SetPosDisShape("Sphere");
  SetCentreCoords(G4ThreeVector(0., 0., 0.));
  SetRadius(radius);

  string mymat = "GdAcrylic";
  if (DSStorage::Get()->GetPurePMMATPC() || DSStorage::Get()->GetHybridTPC()) mymat = "Acrylic";

  G4ThreeVector myPos(0., 0., 0.);
  myPos = fSPSPos->GenerateOne();
  while (!(CheckVolumeName(myPos, "GdSkin") && CheckMaterial(myPos, mymat))) myPos = fSPSPos->GenerateOne();

  return myPos;
}

G4ThreeVector DSVGenerator::GenerateInOPsAcrylic() {

  static const G4double radius = 4 * m;

  SetIsVolumeDistribution(true);
  SetPosDisType("Volume");
  SetPosDisShape("Sphere");
  SetCentreCoords(G4ThreeVector(0., 0., 0.));
  SetRadius(radius);

  string mymat = "GdAcrylic";
  if (DSStorage::Get()->GetPurePMMATPC()) mymat = "Acrylic";

  G4ThreeVector myPos(0., 0., 0.);
  myPos = fSPSPos->GenerateOne();
  while (!( (CheckVolumeName(myPos, "PlaceholderTop") || CheckVolumeName(myPos, "PlaceholderBot")) && CheckMaterial(myPos, mymat))) myPos = fSPSPos->GenerateOne();

  return myPos;
}

G4ThreeVector DSVGenerator::GenerateInArgonBufferInside() {  // inside

  static const G4double radius = 5 * m;

  SetIsVolumeDistribution(true);
  SetPosDisType("Volume");
  SetPosDisShape("Sphere");
  SetCentreCoords(G4ThreeVector(0., 0., 0.));
  SetRadius(radius);

  G4ThreeVector myPos(0., 0., 0.);
  myPos = fSPSPos->GenerateOne();
  while (!CheckMaterial(myPos, "LiquidArgonPlasticVeto")) myPos = fSPSPos->GenerateOne();

  return myPos;
}

G4ThreeVector DSVGenerator::GenerateInArgonBufferOutside() {  // outside

  static const G4double radius = 5 * m;

  SetIsVolumeDistribution(true);
  SetPosDisType("Volume");
  SetPosDisShape("Sphere");
  SetCentreCoords(G4ThreeVector(0., 0., 0.));
  SetRadius(radius);

  G4ThreeVector myPos(0., 0., 0.);
  myPos = fSPSPos->GenerateOne();
  while (!CheckMaterial(myPos, "IVLiquidArgon")) myPos = fSPSPos->GenerateOne();

  return myPos;
}

G4ThreeVector DSVGenerator::GenerateInDS20kSupportStructure() {
 
  G4double pos_z;
  G4double half_z;
  G4double radius;

  pos_z = -DSStorage::Get()->GetDS20kTPCheight() / 2. - DSStorage::Get()->GetDS20kGasPocketThickness() - DSStorage::Get()->GetSiPMOffset() - 6.5 * cm - 7 * cm;
  half_z = 10. * cm;
  radius = DSStorage::Get()->GetDS20kTPCedge() / 2. * (1 + sqrt(2)) + DSStorage::Get()->GetDS20kTPCVesselThickness() + DSStorage::Get()->GetDS20kWLSPENThickness() + DSStorage::Get()->GetDS20kTPBThickness() + 5.4 * cm;
  radius *= 1.1; // to get outer radius

  if (G4UniformRand() > 0.5) {
    pos_z = DSStorage::Get()->GetDS20kTPCheight() / 2. + DSStorage::Get()->GetDS20kGasPocketThickness() + DSStorage::Get()->GetSiPMOffset() + 8. * cm + 7. * cm;
  }
  SetIsVolumeDistribution(true);
  SetPosDisType("Volume");
  SetPosDisShape("Cylinder");

  SetCentreCoords(G4ThreeVector(0, 0, pos_z));
  SetHalfZ(half_z);
  SetRadius(radius);

  G4ThreeVector myPos(0., 0., 0.);
  myPos = fSPSPos->GenerateOne();
  while (!CheckMaterial(myPos, "StainlessSteel")) myPos = fSPSPos->GenerateOne();
  
  return myPos;
}

G4ThreeVector DSVGenerator::GenerateIn7mSphere() {

  static const G4double radius = 3.5 * m;

  SetIsVolumeDistribution(true);
  SetPosDisType("Volume");
  SetPosDisShape("Sphere");
  SetCentreCoords(G4ThreeVector(0., 0., 0.));
  SetRadius(radius);

  G4ThreeVector myPos(0., 0., 0.);
  myPos = fSPSPos->GenerateOne();

  G4ThreeVector null(0., 0., 0.);
  G4ThreeVector* ptr;
  ptr = &null;

  G4VPhysicalVolume* theVolume;
  theVolume = gNavigator->LocateGlobalPointAndSetup(myPos, ptr, true);
  G4Material* amat = theVolume->GetLogicalVolume()->GetMaterial();
  G4int theMatIndex = amat->GetIndex();

  DSEventHandler::Get()->SetUserInt1(theMatIndex);
  DSEventHandler::Get()->SetUsers();

  return myPos;
}

G4ThreeVector DSVGenerator::GenerateInVetoPMTs() {

  static const G4double radius = DSStorage::Get()->GetIDiameter() / 2 + 10 * cm;

  SetIsVolumeDistribution(true);
  SetPosDisType("Volume");
  SetPosDisShape("Sphere");
  SetCentreCoords(G4ThreeVector(0., 0., 0.));
  SetRadius(radius);

  G4ThreeVector myPos(0., 0., 0.);
  myPos = fSPSPos->GenerateOne();
  while (!CheckMaterial(myPos, "FakePhotocathode")) myPos = fSPSPos->GenerateOne();

  return myPos;
}

G4ThreeVector DSVGenerator::GenerateInVetoSteelSphere() {

  static const G4double radius = DSStorage::Get()->GetIDiameter() / 2 + 10 * cm;

  SetIsVolumeDistribution(true);
  SetPosDisType("Volume");
  SetPosDisShape("Sphere");
  SetCentreCoords(G4ThreeVector(0., 0., 0.));
  SetRadius(radius);

  G4ThreeVector myPos(0., 0., 0.);
  myPos = fSPSPos->GenerateOne();

  while (!((CheckMaterial(myPos, "StainlessSteel") || CheckMaterial(myPos, "Acrylic")) && sqrt(myPos.x() * myPos.x() + myPos.y() * myPos.y() + myPos.z() * myPos.z()) > radius - 50 * cm)) { myPos = fSPSPos->GenerateOne(); }
  return myPos;
}

G4ThreeVector DSVGenerator::GenerateProtoDUNEFoam() {

  static const G4double radius = 9 * m;

  SetIsVolumeDistribution(true);
  SetPosDisType("Volume");
  SetPosDisShape("Cylinder");
  SetCentreCoords(G4ThreeVector(0., 0., 0.));
  SetRadius(radius);

  SetHalfZ(7 * m);
  G4ThreeVector myPos(0., 0., 0.);
  myPos = fSPSPos->GenerateOne();
  while (!CheckMaterial(myPos, "InsulatingFoam")) myPos = fSPSPos->GenerateOne();

  return myPos;
}

G4ThreeVector DSVGenerator::GenerateProtoDUNEMembrane() {

  // static const G4double radius = 9 * m;

  SetIsVolumeDistribution(true);
  SetPosDisType("Surface");
  SetPosDisShape("Parallelepiped");
  SetCentreCoords(G4ThreeVector(0., 0., 0.));

  double maxZ = 7.9 / 2. * m + 1. * mm;
  double maxXY = 8.548 / 2. * m + 1 * mm;

  double myZ, myX, myY;

  if (G4UniformRand() < 1. / 6) {  // ceiling and floor
    if (G4UniformRand() < 0.5)     // top
      myZ = maxZ;
    else
      myZ = -maxZ;
    myX = -maxXY + 2 * maxXY * G4UniformRand();
    myY = -maxXY + 2 * maxXY * G4UniformRand();
  } else {
    double mywall = G4UniformRand();
    if (mywall < 0.25) {  // one wall
      myX = maxXY;
      myZ = -maxZ + 2 * maxZ * G4UniformRand();
      myY = -maxXY + 2 * maxXY * G4UniformRand();
    } else if (mywall >= .25 && mywall < .5) {
      myX = -maxXY;
      myZ = -maxZ + 2 * maxZ * G4UniformRand();
      myY = -maxXY + 2 * maxXY * G4UniformRand();
    } else if (mywall >= .5 && mywall < .75) {
      myY = maxXY;
      myZ = -maxZ + 2 * maxZ * G4UniformRand();
      myX = -maxXY + 2 * maxXY * G4UniformRand();
    } else {
      myY = -maxXY;
      myZ = -maxZ + 2 * maxZ * G4UniformRand();
      myX = -maxXY + 2 * maxXY * G4UniformRand();
    }
  }
  // SetHalfZ( 3*m  );
  // SetHalfX (3*m ) ;
  // SetHalfY (3*m );

  G4ThreeVector myPos(myX, myY, myZ);
  // myPos = fSPSPos->GenerateOne();
  // while( !CheckMaterial( myPos, "InsulatingFoam" ) ) myPos =
  // fSPSPos->GenerateOne();

  return myPos;
}

G4ThreeVector DSVGenerator::GenerateInGantryPipe() {
  static const G4double radius = 2.2 * m;
  static const G4double height = 2.2 * m;

  SetIsVolumeDistribution(true);
  SetPosDisType("Volume");
  SetPosDisShape("Cylinder");

  SetCentreCoords(G4ThreeVector(0., 0., 0.));
  SetRadius(radius);
  SetHalfZ(height);

  G4ThreeVector myPos(0., 0., 0.);
  myPos = fSPSPos->GenerateOne();
  // later might need to create a separate material for gantry steel in the
  // DSMaterial
  while (!CheckMaterial(myPos, "StainlessSteel")) myPos = fSPSPos->GenerateOne();

  return myPos;
}

G4bool DSVGenerator::CheckMaterial(G4ThreeVector pos, G4String MatName) {

  G4bool found = false;
  fNumberOfHits++;
  DSStorage::Get()->SetNumberOfHits(fNumberOfHits);

  G4ThreeVector null(0., 0., 0.);
  G4ThreeVector* ptr;
  ptr = &null;

  G4VPhysicalVolume* theVolume;
  theVolume = gNavigator->LocateGlobalPointAndSetup(pos, ptr, true);
  G4Material* amat = theVolume->GetLogicalVolume()->GetMaterial();
  G4String theMatName = amat->GetName();
  if (theMatName == MatName) found = true;

  return found;
}

G4bool DSVGenerator::CheckMaterialIndex(G4ThreeVector pos, G4int MatIndex) {

  G4bool found = false;
  fNumberOfHits++;
  DSStorage::Get()->SetNumberOfHits(fNumberOfHits);

  G4ThreeVector null(0., 0., 0.);
  G4ThreeVector* ptr;
  ptr = &null;

  G4VPhysicalVolume* theVolume;
  theVolume = gNavigator->LocateGlobalPointAndSetup(pos, ptr, true);
  G4Material* amat = theVolume->GetLogicalVolume()->GetMaterial();
  G4int theMatIndex = amat->GetIndex();
  if (theMatIndex == MatIndex) found = true;

  return found;
}

G4bool DSVGenerator::CheckVolumeName(G4ThreeVector pos, G4String MatName) {

  G4bool found = false;
  fNumberOfHits++;
  DSStorage::Get()->SetNumberOfHits(fNumberOfHits);

  G4ThreeVector null(0., 0., 0.);
  G4ThreeVector* ptr;
  ptr = &null;

  G4VPhysicalVolume* theVolume;
  theVolume = gNavigator->LocateGlobalPointAndSetup(pos, ptr, true);
  G4String theMatName = theVolume->GetName();
  if (theMatName == MatName) found = true;
  return found;
}


/*
G4ThreeVector DSVGenerator::GenerateInCopperVessel(){

  G4double gap_copper = 2.*cm;

  G4double radius = (222.*cm)*2/(sqrt(2+sqrt(2))) + gap_copper;
  G4double height = 230.*cm + 2.*gap_copper;

  SetIsVolumeDistribution(true);
  SetPosDisType("Volume");
  SetPosDisShape("Cylinder");

  SetCentreCoords(G4ThreeVector(0.,0.,0.));
  SetRadius(radius);
  SetHalfZ(height);

  G4ThreeVector myPos(0.,0.,0.);
  myPos = fSPSPos->GenerateOne();
   while( !CheckMaterial( myPos, "MetalCopper" ) ) myPos =
fSPSPos->GenerateOne(); return myPos;
}
*/

G4ThreeVector DSVGenerator::GenerateInTitaniumVessel() {

  G4double gap_titanium = 20. * cm;

  G4double radius = (232.5 * cm) * 2 / (sqrt(2 + sqrt(2))) + gap_titanium;
  G4double height_up = 269.285 * cm + gap_titanium;
  G4double height_down = 269.285 * cm + gap_titanium;

  SetIsVolumeDistribution(true);
  SetPosDisType("Volume");
  SetPosDisShape("Cylinder");

  SetCentreCoords(G4ThreeVector(0., 0., (height_up - height_down) * 0.5));
  SetRadius(radius);
  SetHalfZ((height_up + height_down) * 0.5);

  G4double z_min = 200. * cm;
  G4double r_min = 220. * cm;

  string mymat = "StainlessSteel";
  if (DSStorage::Get()->Get20KGeometry()) {
    if (DSStorage::Get()->GetDS20kCryoMaterial() == 0) mymat = "StainlessSteel";
    else if (DSStorage::Get()->GetDS20kCryoMaterial() == 1)
      mymat = "MetalTitanium";
    else if (DSStorage::Get()->GetDS20kCryoMaterial() == 2)
      mymat = "MetalCopperCryo";
  }

  G4ThreeVector myPos(0., 0., 0.);
  myPos = fSPSPos->GenerateOne();
  while (!CheckMaterial(myPos, mymat) || (abs(myPos.getZ()) < z_min && sqrt(myPos.getX() * myPos.getX() + myPos.getY() * myPos.getY()) < r_min)) myPos = fSPSPos->GenerateOne();
  return myPos;
}

G4ThreeVector DSVGenerator::GenerateInAriserSourceHolder(){
  G4ThreeVector center = DSStorage::Get()->GetAriserSourcePosition();
  G4double sourceActiveDiameter = 2.1 * mm;

  SetIsVolumeDistribution(true);
  SetPosDisType("Volume");
  SetPosDisShape("Sphere");
  SetCentreCoords(center);
  SetRadius(1.05 * sourceActiveDiameter / 2.);

  G4ThreeVector myPos(0.,0.,0.);
  myPos = fSPSPos->GenerateOne();
  while( !CheckMaterial(myPos, "Ceramic") ) myPos = fSPSPos->GenerateOne();
  // DSLog(development) << "Na22 Decay Position (cm): " << myPos.getX()/cm << "  " << myPos.getY()/cm << "  " << myPos.getZ()/cm << endlog;
  return myPos;
}

G4ThreeVector DSVGenerator::GenerateInPhysVolume() {

  G4String VolName = GetGenInPhysVolumeName();
 
  static const G4double radius = 4 * m;

  SetIsVolumeDistribution(true);
  SetPosDisType("Volume");
  SetPosDisShape("Sphere");
  SetCentreCoords(G4ThreeVector(0., 0., 0.));
  SetRadius(radius);

  G4ThreeVector myPos(0., 0., 0.);
  myPos = fSPSPos->GenerateOne();
  while (!(CheckVolumeName(myPos, VolName))) myPos = fSPSPos->GenerateOne();

  return myPos;
 
}

G4ThreeVector DSVGenerator::GenerateInAnode() {

  static const G4double radius = 4 * m;

  SetIsVolumeDistribution(true);
  SetPosDisType("Volume");
  SetPosDisShape("Sphere");
  SetCentreCoords(G4ThreeVector(0., 0., 0.));
  SetRadius(radius);

  G4ThreeVector myPos(0., 0., 0.);
  myPos = fSPSPos->GenerateOne();
  while (!( (CheckVolumeName(myPos, "AcrylicTPCVessel") && myPos.getZ() > 0) || CheckVolumeName(myPos, "AcrylicTopWindowExtension"))) myPos = fSPSPos->GenerateOne();

  return myPos;

}

G4ThreeVector DSVGenerator::GenerateInCathode() {

  static const G4double radius = 4 * m;

  SetIsVolumeDistribution(true);
  SetPosDisType("Volume");
  SetPosDisShape("Sphere");
  SetCentreCoords(G4ThreeVector(0., 0., 0.));
  SetRadius(radius);

  G4ThreeVector myPos(0., 0., 0.);
  myPos = fSPSPos->GenerateOne();
  while (!( (CheckVolumeName(myPos, "AcrylicTPCVessel") && myPos.getZ() < 0) || CheckVolumeName(myPos, "AcrylicBottomWindowExtension"))) myPos = fSPSPos->GenerateOne();

  return myPos;

}
