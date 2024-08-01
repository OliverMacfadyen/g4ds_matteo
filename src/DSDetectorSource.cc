#include "DSDetectorSource.hh"
#include <iostream>
#include "DSIO.hh"
#include "DSLogger.hh"
#include "DSMaterial.hh"
#include "DSParameters.hh"
#include "DSStorage.hh"
#include "G4DisplacedSolid.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalConstants.hh"
#include "G4SubtractionSolid.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VisAttributes.hh"

using namespace std;

////////////////////////////////////////////////
//////        Detector Description      ////////
////////////////////////////////////////////////
/*

From the center (0,0,0) to the top

LiquidArgon     ActiveLAr          27.537 cm
GaseousArgon    GasPocket          28.474 cm
TPB             TPB                28.484 cm
FusedSilica     BellTop            28.802 cm
LiquidArgon     InnerLiqArgon      29.304 cm
LiquidArgon     PMTLArDisk_0       29.404 cm
Bialkali        TPMT_0             29.7627 cm

Back of the PMT:

Vacuum          PMTHeadVac_0       32.853 cm
Vacuum          PMTJoinVac_0       32.908 cm
Vacuum          PMTBodyVac_0       36.365 cm
Vacuum          PMTTopVac_0        41.592 cm
Kovar           PMTTop_0           41.647 cm
LiquidArgon     PMTAssemblyTop     43.528 cm
LiquidArgon     OuterLiquidArgon   44.4 cm
GaseousArgon    GaseousArgon       62.141 cm
StainlessSteel  InnerCryostat      62.6 cm
Vacuum          VacuumCryostat     63.4 cm
GaseousArgon    TrunkAr            123.4 cm
Air             World              126.301 cm
StainlessSteel  Trunk6             126.584 cm
GaseousArgon    TrunkAr            131.406 cm
StainlessSteel  Trunk6             131.689 cm
Air             World              2000 cm


From the center (0,0,0) to the left

LiquidArgon     ActiveLAr          17.77 cm
TPB             TPB                17.78 cm
Teflon          Reflector          20.32 cm
MetalCopper     FieldRings         20.7367 cm
LiquidArgon     InnerLiqArgon      21.59 cm
Teflon          TeflonSupport      23.495 cm
LiquidArgon     OuterLiquidArgon   25.197 cm
StainlessSteel  InnerCryostat      25.65 cm
Vacuum          VacuumCryostat     31.6983 cm
StainlessSteel  OuterCryostat      32.1 cm


*/
////////////////////////////////////////////////

DSDetectorSource::DSDetectorSource(G4VPhysicalVolume* myMotherVolume) {

  fMotherVolume = myMotherVolume;
  G4int volumeNumber = 20160209;
  /*
  const double myTwoPi = 2*M_PI*rad;
  bool   myCheckOverlap   = DSStorage::Get()->GetCheckOverlap();

  fPosition = DSStorage::Get()->GetSourcePosition();


  DSLog(routine) << " Constructing Source Geometry" << endlog ;

  G4RotationMatrix* rotZ90 = new G4RotationMatrix;
  rotZ90->rotateZ( M_PI/2.*rad );

  G4Tubs *mySolidVial            = new G4Tubs ("Vial_Solid", 0, 14.5*mm, 5*mm,
  0, myTwoPi); G4LogicalVolume *myLogicVial   = new G4LogicalVolume(mySolidVial,
  DSMaterial::Get()->GetStainlessSteel(), "Vial_Logic"); G4VPhysicalVolume
  *myPhysicVial = new G4PVPlacement( rotZ90, fPosition, "Vial", myLogicVial,
  fMotherVolume, true, 0, myCheckOverlap );


  G4Tubs *mySolidSource             = new G4Tubs ("Source_Solid", 0, 12.5*mm,
  3*mm, 0, myTwoPi); G4LogicalVolume *myLogicSource    = new
  G4LogicalVolume(mySolidSource, DSMaterial::Get()->GetTeflon(),
  "Source_Logic"); new G4PVPlacement( 0, G4ThreeVector(0,0,0), "Source",
  myLogicSource, myPhysicVial, true, 1, myCheckOverlap );
  */

  // Source holder
  if (DSStorage::Get()->GetSourceHolderFlag()) {
    DSLog(routine) << "Constructing Source Holder inside DS50 Neutron Veto: " << volumeNumber + 1 << endlog;

    // geometry of Source Holder Stainless steel
    G4double sourceholderouterradius = 12.5 * mm;  // total outer radius 25mm
    G4double sourceholderlength = 75 * mm;
    // geometry of Vacuum inside Source Holder
    G4double sourceholdervacuumradius = 10.5 * mm;  // gap between statinless steel and vacuum is 2mm
    G4double sourceholdervacuumlength = 64 * mm;
    G4double sourceholderuppergap = 10 * mm;  // z gap between source holder and vacuum

    G4ThreeVector mysourceholdercenter = DSStorage::Get()->GetSourceHolderCenter();
    G4double theta = DSStorage::Get()->GetSourceHolderTheta();
    G4double phi = DSStorage::Get()->GetSourceHolderPhi();
    G4RotationMatrix* myHolderRotation = new G4RotationMatrix;
    myHolderRotation->rotateZ(phi);
    myHolderRotation->rotateX(theta);

    mySolidSteelHolder = new G4Tubs("SteelHolder_Solid", 0, sourceholderouterradius, sourceholderlength / 2.0, 0, twopi * rad);
    myLogicSteelHolder = new G4LogicalVolume(mySolidSteelHolder, DSMaterial::Get()->GetStainlessSteel(), "SteelHolder_Logic");
    myPhysicalSteelHolder = new G4PVPlacement(myHolderRotation, mysourceholdercenter, "SteelHolder", myLogicSteelHolder, fMotherVolume, false, volumeNumber++, DSStorage::Get()->GetCheckOverlap());
    G4VisAttributes* rot = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0));
    myLogicSteelHolder->SetVisAttributes(rot);
    DSLog(routine) << "Source holder Rotation Theta: " << myPhysicalSteelHolder->GetObjectRotation()->getTheta() << " Phi: " << myPhysicalSteelHolder->GetObjectRotation()->getPhi() << endlog;
    DSLog(routine) << "Source holder Center: " << myPhysicalSteelHolder->GetObjectTranslation() << endlog;
    G4double sourceholdervacuumshifter = sourceholderlength / 2. - sourceholderuppergap - sourceholdervacuumlength / 2.;  // in terms of source holder
    mySolidVacuumHolder = new G4Tubs("VacuumHolder_Solid", 0, sourceholdervacuumradius, sourceholdervacuumlength / 2., 0, twopi * rad);
    myLogicVacuumHolder = new G4LogicalVolume(mySolidVacuumHolder, DSMaterial::Get()->GetVacuum(), "VacuumHolder_Logic");
    myPhysicalVacuumHolder = new G4PVPlacement(0, G4ThreeVector(0, 0, sourceholdervacuumshifter), "VacuumHolder", myLogicVacuumHolder, myPhysicalSteelHolder, false, volumeNumber++, DSStorage::Get()->GetCheckOverlap());
    G4VisAttributes* rosa = new G4VisAttributes(G4Colour(1.0, 0.0, 1.0));
    myLogicVacuumHolder->SetVisAttributes(rosa);
    DSLog(routine) << "Source holder Vacuum Center: " << myPhysicalVacuumHolder->GetObjectTranslation() << endlog;
    if (DSStorage::Get()->GetSourceHolderLeadFlag()) {  // turn on the Lead
                                                        // shield
      // geometry of lead shield(if it's on)
      G4double sourceholderleadradius = 8.5 * mm;  // 17mm/2.
      G4double sourceholderleadlength = 6.35 * mm;
      G4double sourceholderleadthickness = 2 * mm;
      G4double sourceholdercontainerradius = sourceholderleadradius - sourceholderleadthickness;
      G4double sourceholdercontainerlength = sourceholderleadlength - sourceholderleadthickness * 2;

      G4double sourceholderleadshieldshifter = sourceholderleadlength / 2. - sourceholdervacuumlength / 2.;  // in terms of vacuum

      mySolidLeadShield = new G4Tubs("LeadShield_Solid", 0, sourceholderleadradius, sourceholderleadlength / 2., 0, twopi * rad);
      myLogicLeadShield = new G4LogicalVolume(mySolidLeadShield, DSMaterial::Get()->GetMetalLead(), "LeadShield_Logic");
      myPhysicalLeadShield = new G4PVPlacement(0, G4ThreeVector(0, 0, sourceholderleadshieldshifter), "LeadShield", myLogicLeadShield, myPhysicalVacuumHolder, false, volumeNumber++, DSStorage::Get()->GetCheckOverlap());
      myLogicVacuumHolder->SetVisAttributes(G4Colour(1.0, 1.0, 0));
      DSLog(routine) << "Source holder Lead Shield Center: " << myPhysicalLeadShield->GetObjectTranslation() << endlog;
      mySolidSourceContainer = new G4Tubs("SourceContainer_Solid", 0, sourceholdercontainerradius, sourceholdercontainerlength / 2., 0, twopi * rad);
      // myLogicSourceContainer      = new
      // G4LogicalVolume(mySolidSourceContainer,DSMaterial::Get()->GetTeflon(),"SourceContainer_Logic");
      myLogicSourceContainer = new G4LogicalVolume(mySolidSourceContainer, DSMaterial::Get()->GetAir(), "SourceContainer_Logic");
      myPhysicalSourceContainer = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), "SourceContainer", myLogicSourceContainer, myPhysicalLeadShield, false, volumeNumber++, DSStorage::Get()->GetCheckOverlap());
      myLogicSourceContainer->SetVisAttributes(G4Colour(1.0, 1.0, 0));
      DSLog(routine) << "Source holder Source Container Center: " << myPhysicalSourceContainer->GetObjectTranslation() << endlog;

    } else {
      // geometry of teflon(lead shield is off)
      G4double sourceholderteflonradius = sourceholdervacuumradius;
      G4double sourceholderteflonlength = 3 * mm;
      G4double sourceholderteflonshifter = sourceholderteflonlength / 2. - sourceholdervacuumlength / 2.;  // in terms of vaccum

      mySolidSourceDisk = new G4Tubs("SourceDisk_Solid", 0, sourceholderteflonradius, sourceholderteflonlength / 2., 0, twopi * rad);
      // myLogicSourceDisk      = new G4LogicalVolume(mySolidSourceDisk,
      // DSMaterial::Get()->GetTeflon(), "SourceDisk_Logic");
      myLogicSourceDisk = new G4LogicalVolume(mySolidSourceDisk, DSMaterial::Get()->GetAir(), "SourceDisk_Logic");
      myPhysicalSourceDisk = new G4PVPlacement(0, G4ThreeVector(0, 0, sourceholderteflonshifter), "SourceDisk", myLogicSourceDisk, myPhysicalVacuumHolder, false, volumeNumber++, DSStorage::Get()->GetCheckOverlap());
      G4VisAttributes* green = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0));
      myLogicSourceDisk->SetVisAttributes(green);
      DSLog(routine) << "Source holder Teflon Center: " << myPhysicalSourceDisk->GetObjectTranslation() << endlog;
    }
  }

  DefineSurfaces();
}

DSDetectorSource::~DSDetectorSource() {
  ;  // delete fMessenger;
}

void DSDetectorSource::DefineSurfaces() {
  ;
}

/*
 * $Log: DSDetectorSource.cc,v $
 * Revision 1.2  2016/03/08 21:09:36  swesterd
 * Added source holder
 *
 * Revision 1.1  2014/11/06 17:39:45  dfranco
 * added source and fixed a bug: removed gaseous nitrogen from the neutron veto
 *
 *
 */
