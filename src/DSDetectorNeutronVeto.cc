#include "DSDetectorNeutronVeto.hh"
#include <fstream>
#include <iostream>
#include "DSDetectorPMTNeutronVeto.hh"
#include "DSLogger.hh"
#include "DSMaterial.hh"
#include "DSParameters.hh"
#include "DSStorage.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4VisAttributes.hh"
using namespace std;

DSDetectorNeutronVeto::DSDetectorNeutronVeto(G4VPhysicalVolume* myMotherVolume) {

  fMotherVolume = myMotherVolume;

  G4int volumeNumber = 10000;

  DSLog(routine) << " Constructing NeutronVeto Geometry" << endlog;

  // Stainless Steel Vessel

  fSolidSteelVessel = new G4Orb("SteelVessel_Solid", 2008. * mm);
  fLogicSteelVessel = new G4LogicalVolume(fSolidSteelVessel, DSMaterial::Get()->GetStainlessSteel(), "SteelVessel_Logic");
  fPhysicSteelVessel = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), "SteelVessel", fLogicSteelVessel, fMotherVolume, false, volumeNumber++, DSStorage::Get()->GetCheckOverlap());

  DSLog(routine) << fPhysicSteelVessel->GetName() << " = " << fPhysicSteelVessel->GetCopyNo() << endlog;

  G4Material* myFillMaterial = 0;
  if (DSStorage::Get()->GetScintillator() == 0) {
    myFillMaterial = DSMaterial::Get()->GetBoronScintillator();
    DSLog(routine) << "Borate Scintillator in the DS20k neutron veto" << endlog;
  } else if (DSStorage::Get()->GetScintillator() == 1) {
    myFillMaterial = DSMaterial::Get()->GetGdScintillator();
    DSLog(routine) << "Gd Scintillator in the DS20k neutron veto" << endlog;
  } else if (DSStorage::Get()->GetScintillator() == 2) {
    myFillMaterial = DSMaterial::Get()->GetLi6Scintillator();
    DSLog(routine) << "Li6-enriched Scintillator in the DS20k neutron veto" << endlog;
  } else if (DSStorage::Get()->GetScintillator() == 3) {
    myFillMaterial = DSMaterial::Get()->GetNatLi6Scintillator();
    DSLog(routine) << "Li6-natural Scintillator in the DS20k neutron veto" << endlog;
    G4cout << myFillMaterial << endl;
  } else if (DSStorage::Get()->GetScintillator() == 4) {
    myFillMaterial = DSMaterial::Get()->GetOrtoCarbScintillator();
    DSLog(routine) << "Orto-Carbo Scintillator in the DS20k neutron veto" << endlog;
  } else {
    DSLog(fatal) << "Scintillator index " << DSStorage::Get()->GetScintillator() << " not defined" << endlog;
  }

  fSolidBScintillator = new G4Orb("BoronScintillator_Solid", 2000. * mm);
  fLogicBScintillator = new G4LogicalVolume(fSolidBScintillator, myFillMaterial, "BoronScintillator_Logic");
  fPhysicBScintillator = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), "BoronScintillatorVolume", fLogicBScintillator, fPhysicSteelVessel, false, volumeNumber++, DSStorage::Get()->GetCheckOverlap());

  DSLog(routine) << fPhysicBScintillator->GetName() << " = " << fPhysicBScintillator->GetCopyNo() << endlog;

  // PMT
  new DSDetectorPMTNeutronVeto(fPhysicBScintillator);
  /*
  //Lumirror Sheath outside flat part of cryostat
  G4double sheathGap = 1.*nm;
  G4double cryoSheathMaxZ = DSParameters::Get()->GetCryoSheathZ()[0];
  G4double cryoSheathMinZ = DSParameters::Get()->GetCryoSheathZ()[1];
  G4double cryoSheathR = DSParameters::Get()->GetCryoSheathR()[0]+sheathGap;
  G4double sheathThickness = 1.*nm;
  fSolidCryoSheath = new G4Tubs("CryoSheath_Solid",
                                cryoSheathR,
                                cryoSheathR+sheathThickness,
                                (cryoSheathMaxZ - cryoSheathMinZ)/2.,
                                0, 2*M_PI);
  fLogicCryoSheath = new G4LogicalVolume(fSolidCryoSheath,
  DSMaterial::Get()->GetBoronScintillator(), "CryoSheath_Logic");
  fPhysicCryoSheath = new G4PVPlacement(0,
                                        G4ThreeVector(0,0,(cryoSheathMaxZ+cryoSheathMinZ)/2.),
                                        "CryoSheath",
                                        fLogicCryoSheath,
                                        fPhysicBScintillator,
                                        false,
                                        volumeNumber++,
                                        DSStorage::Get()->GetCheckOverlap());
  DSLog(routine) << fPhysicCryoSheath->GetName() << " = " <<
  fPhysicCryoSheath->GetCopyNo() << endlog;

  fLogicCryoSheath->SetVisAttributes(new G4VisAttributes(G4Colour(0,.5,.5)));
  */

  // Cover the top flange with untreated SS
  G4double topFlangeRadius = 45. * cm;
  G4double topFlangeAngle = std::asin(topFlangeRadius / (200. * cm));
  G4double topFlangeThickness = 8. * mm;
  fSolidTopFlange = new G4Sphere("TopFlange_Solid", 200. * cm - topFlangeThickness, 200. * cm, 0, 2 * M_PI, 0, topFlangeAngle);
  fLogicTopFlange = new G4LogicalVolume(fSolidTopFlange, DSMaterial::Get()->GetStainlessSteel(), "TopFlange_Logic");
  fPhysicTopFlange = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), "TopFlangeVolume", fLogicTopFlange, fPhysicBScintillator, false, volumeNumber++, DSStorage::Get()->GetCheckOverlap());

  DSLog(routine) << fPhysicTopFlange->GetName() << " = " << fPhysicTopFlange->GetCopyNo() << endlog;

  /*
//Source holder
  if(DSStorage::Get()->GetSourceHolderFlag()){
    DSLog(routine) << "Constructing Source Holder inside DS50 Neutron Veto: " <<
volumeNumber+1 << endlog;

    //geometry of Source Holder Stainless steel
    G4double sourceholderouterradius = 12.5*mm; //total outer radius 25mm
    G4double sourceholderlength = 75*mm;
    //geometry of Vacuum inside Source Holder
    G4double sourceholdervacuumradius = 10.5*mm; //gap between statinless steel
and vacuum is 2mm G4double sourceholdervacuumlength = 64*mm; G4double
sourceholderuppergap = 10*mm; //z gap between source holder and vacuum

    G4ThreeVector mysourceholdercenter =
DSStorage::Get()->GetSourceHolderCenter(); G4double theta =
DSStorage::Get()->GetSourceHolderTheta(); G4double phi   =
DSStorage::Get()->GetSourceHolderPhi(); G4RotationMatrix *myHolderRotation = new
G4RotationMatrix; myHolderRotation->rotateZ(phi);
    myHolderRotation->rotateX(theta);

    mySolidSteelHolder  = new
G4Tubs("SteelHolder_Solid",0,sourceholderouterradius,sourceholderlength/2.0,0,twopi*rad);
    myLogicSteelHolder  = new G4LogicalVolume(mySolidSteelHolder,
DSMaterial::Get()->GetStainlessSteel(), "SteelHolder_Logic");
    myPhysicalSteelHolder  = new G4PVPlacement(myHolderRotation,
                                               mysourceholdercenter,
                                               "SteelHolder",
                                               myLogicSteelHolder,
                                               fPhysicBScintillator,
                                               false,
                                               volumeNumber++,
                                               DSStorage::Get()->GetCheckOverlap());
    G4VisAttributes* rot       = new G4VisAttributes(G4Colour(1.0,0.0,0.0));
    myLogicSteelHolder->SetVisAttributes(rot);
    DSLog(routine)<<"Source holder Rotation Theta:
"<<myPhysicalSteelHolder->GetObjectRotation()->getTheta()<<" Phi:
"<<myPhysicalSteelHolder->GetObjectRotation()->getPhi()<<endlog;
    DSLog(routine)<<"Source holder Center:
"<<myPhysicalSteelHolder->GetObjectTranslation()<<endlog; G4double
sourceholdervacuumshifter =
sourceholderlength/2.-sourceholderuppergap-sourceholdervacuumlength/2.;//in
terms of source holder mySolidVacuumHolder    = new G4Tubs("VacuumHolder_Solid",
0,sourceholdervacuumradius,sourceholdervacuumlength/2.,0,twopi*rad);
    myLogicVacuumHolder    = new G4LogicalVolume(mySolidVacuumHolder,
DSMaterial::Get()->GetVacuum(), "VacuumHolder_Logic"); myPhysicalVacuumHolder =
new G4PVPlacement(0, G4ThreeVector(0,0,sourceholdervacuumshifter),
                                               "VacuumHolder",
                                               myLogicVacuumHolder,
                                               myPhysicalSteelHolder,
                                               false,
                                               volumeNumber++,
                                               DSStorage::Get()->GetCheckOverlap());
    G4VisAttributes* rosa      = new G4VisAttributes(G4Colour(1.0,0.0,1.0));
    myLogicVacuumHolder->SetVisAttributes(rosa);
    DSLog(routine)<<"Source holder Vacuum Center:
"<<myPhysicalVacuumHolder->GetObjectTranslation()<<endlog;

    if(DSStorage::Get()->GetSourceHolderLeadFlag()){ //turn on the Lead shield
      //geometry of lead shield(if it's on)
      G4double sourceholderleadradius = 8.5*mm; //17mm/2.
      G4double sourceholderleadlength = 6.35*mm;
      G4double sourceholderleadthickness = 2*mm;
      G4double sourceholdercontainerradius =
sourceholderleadradius-sourceholderleadthickness; G4double
sourceholdercontainerlength =
sourceholderleadlength-sourceholderleadthickness*2;

      G4double sourceholderleadshieldshifter =
sourceholderleadlength/2.-sourceholdervacuumlength/2.;//in terms of vacuum

      mySolidLeadShield      = new
G4Tubs("LeadShield_Solid",0,sourceholderleadradius,sourceholderleadlength/2.,0,twopi*rad);
      myLogicLeadShield      = new
G4LogicalVolume(mySolidLeadShield,DSMaterial::Get()->GetMetalLead(),"LeadShield_Logic");
      myPhysicalLeadShield   = new G4PVPlacement(0,
                                                 G4ThreeVector(0,0,sourceholderleadshieldshifter),
                                                 "LeadShield",
                                                 myLogicLeadShield,
                                                 myPhysicalVacuumHolder,
                                                 false,
                                                 volumeNumber++,
                                                 DSStorage::Get()->GetCheckOverlap());
      myLogicVacuumHolder->SetVisAttributes(G4Colour(1.0,1.0,0));
      DSLog(routine)<<"Source holder Lead Shield Center:
"<<myPhysicalLeadShield->GetObjectTranslation()<<endlog; mySolidSourceContainer
= new
G4Tubs("SourceContainer_Solid",0,sourceholdercontainerradius,sourceholdercontainerlength/2.,0,twopi*rad);
      // myLogicSourceContainer      = new
G4LogicalVolume(mySolidSourceContainer,DSMaterial::Get()->GetTeflon(),"SourceContainer_Logic");
      myLogicSourceContainer      = new
G4LogicalVolume(mySolidSourceContainer,DSMaterial::Get()->GetAir(),"SourceContainer_Logic");
      myPhysicalSourceContainer   = new G4PVPlacement(0,
                                                      G4ThreeVector(0,0,0),
                                                      "SourceContainer",
                                                      myLogicSourceContainer,
                                                      myPhysicalLeadShield,
                                                      false,
                                                      volumeNumber++,
                                                      DSStorage::Get()->GetCheckOverlap());
      myLogicSourceContainer->SetVisAttributes(G4Colour(1.0,1.0,0));
      DSLog(routine)<<"Source holder Source Container Center:
"<<myPhysicalSourceContainer->GetObjectTranslation()<<endlog;

    }else{
      //geometry of teflon(lead shield is off)
      G4double sourceholderteflonradius = sourceholdervacuumradius;
      G4double sourceholderteflonlength = 3*mm;
      G4double sourceholderteflonshifter =
sourceholderteflonlength/2.-sourceholdervacuumlength/2.; //in terms of vaccum

      mySolidSourceDisk      = new G4Tubs("SourceDisk_Solid",
0,sourceholderteflonradius,sourceholderteflonlength/2.,0,twopi*rad);
      myLogicSourceDisk      = new G4LogicalVolume(mySolidSourceDisk,
DSMaterial::Get()->GetTeflon(), "SourceDisk_Logic"); myPhysicalSourceDisk   =
new G4PVPlacement(0, G4ThreeVector(0,0,sourceholderteflonshifter), "SourceDisk",
                                                 myLogicSourceDisk,
                                                 myPhysicalVacuumHolder,
                                                 false,
                                                 volumeNumber++,
                                                 DSStorage::Get()->GetCheckOverlap());
      G4VisAttributes* green     = new G4VisAttributes(G4Colour(0.0,1.0,0.0));
      myLogicSourceDisk->SetVisAttributes(green);
      DSLog(routine)<<"Source holder Teflon Center:
"<<myPhysicalSourceDisk->GetObjectTranslation()<<endlog;

    }
  }

*/
  DefineSurfaces();
}

DSDetectorNeutronVeto::~DSDetectorNeutronVeto() {
  ;  // delete fMessenger;
}
void DSDetectorNeutronVeto::DefineSurfaces() {
  fOpUntreatedStainlessSteelSurface = new G4OpticalSurface("OpUntreatedStainlessSteelSurface");
  fOpUntreatedStainlessSteelSurface->SetType(dielectric_metal);
  fOpUntreatedStainlessSteelSurface->SetModel(unified);
  fOpUntreatedStainlessSteelSurface->SetFinish(groundbackpainted);
  fOpUntreatedStainlessSteelSurface->SetMaterialPropertiesTable(DSMaterial::Get()->GetUntreatedStainlessSteelMPT());

  fOpElectropolishedStainlessSteelSurface = new G4OpticalSurface("OpElectropolishedStainlessSteelSurface");
  fOpElectropolishedStainlessSteelSurface->SetType(dielectric_metal);
  fOpElectropolishedStainlessSteelSurface->SetModel(unified);
  fOpElectropolishedStainlessSteelSurface->SetFinish(groundbackpainted);
  fOpElectropolishedStainlessSteelSurface->SetMaterialPropertiesTable(DSMaterial::Get()->GetElectropolishedStainlessSteelMPT());

  fOpLumirrorSurface = new G4OpticalSurface("OpLumirrorSurface");
  fOpLumirrorSurface->SetType(dielectric_metal);
  fOpLumirrorSurface->SetModel(unified);
  fOpLumirrorSurface->SetFinish(groundfrontpainted);
  fOpLumirrorSurface->SetMaterialPropertiesTable(DSMaterial::Get()->GetLumirrorMPT());

  fSteelInnerSurface = new G4LogicalBorderSurface("SSteelInnerSurface", fPhysicBScintillator, fPhysicSteelVessel, fOpLumirrorSurface);

  fSteelInnerSurfaceFlip = new G4LogicalBorderSurface("SSteelInnerSurface", fPhysicSteelVessel, fPhysicBScintillator, fOpLumirrorSurface);

  fSteelOuterSurface = new G4LogicalBorderSurface("SSteelOuterSurface", fPhysicSteelVessel, fMotherVolume, fOpUntreatedStainlessSteelSurface);
  /*
  fCryoSheathSurface = new G4LogicalBorderSurface("CryoSheathSurface",
                                                  fPhysicBScintillator,
                                                  fPhysicCryoSheath,
                                                  fOpLumirrorSurface);
  */

  fTopFlangeSurface = new G4LogicalBorderSurface("TopFlangeSurface", fPhysicBScintillator, fPhysicTopFlange, fOpLumirrorSurface);
}

/*
 * $Log: DSDetectorNeutronVeto.cc,v $
 * Revision 1.6  2016/03/08 21:09:36  swesterd
 * Added source holder
 *
 * Revision 1.5  2015/01/14 16:58:36  dfranco
 * added source vial and correspondent commands by Laura and Bianca; manual
 * updated
 *
 * Revision 1.4  2014/11/06 17:39:45  dfranco
 * added source and fixed a bug: removed gaseous nitrogen from the neutron veto
 *
 * Revision 1.3  2014/10/13 18:43:48  swesterd
 * fixed veto PMT positions
 *
 * Revision 1.2  2014/05/07 14:27:31  dfranco
 * fixed some bugs and added GdScintillator
 *
 * Revision 1.1  2014/05/07 12:21:01  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.13  2014/04/18 16:21:42  swesterd
 * fixed an overlap in the G2 detector
 *
 * Revision 1.12  2013/08/27 04:07:01  swesterd
 * some fine tuning of the boron scintillator kB and scint yield, and some
 * modifications to the DSG2 geometry
 *
 * Revision 1.11  2013/08/05 03:13:52  swesterd
 * some fine tuning of bscint and veto parameters
 *
 * Revision 1.10  2013/06/06 23:00:19  swesterd
 * moved veto PMT numbers from copy number to name and made a separate logical
 * volume for each veto PMT. Assigned each physical volume in the veto a unique
 * copy number 1wxyz, where w=1 for PMT bulbs, w=2 for PMT bases, w=3 for
 * photocathodes, and w=0 for everything else
 *
 * Revision 1.9  2013/06/05 23:03:32  swesterd
 * moved optical boundary MPTs to DSMaterial and gave the trunks optical
 * boundary properties consistent with untreated stainless steel
 *
 * Revision 1.8  2013/06/04 01:02:29  swesterd
 * other than the optical boundary of the trunks, the veto optics appear to be
 * complete and up and running...modulo whatever I may have missed...
 *
 * Revision 1.7  2013/05/27 23:59:02  swesterd
 * added a (currently commented out) Lumirror sheath to the cryostat and
 * introduced DSOpBoundaryProcess to try to figure out why the boundaries are
 * being screwy, with some edits so that it can handle constant and vector
 * properties with freaking out
 *
 * Revision 1.6  2013/05/14 04:20:22  swesterd
 * Fixed the veto PMT geometry
 *
 * Revision 1.5  2013/05/07 23:06:29  swesterd
 * added optical boundaries and Lumirror in the veto
 *
 * Revision 1.4  2013/05/07 16:15:45  swesterd
 * Optical processes now seem to be working in the boron-loaded scintillator
 *
 * Revision 1.3  2013/05/07 09:44:05  dfranco
 * Changed logical names of BoronScintillator volumes
 *
 * Revision 1.2  2013/03/22 14:09:40  dfranco
 * Added the cvs logger code at the end of each file
 *
 */
