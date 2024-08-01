#include "DSDetectorLicorneNeutron.hh"

#include <iostream>

#include "DSIO.hh"
#include "DSLogger.hh"
#include "DSMaterial.hh"
#include "DSParameters.hh"
#include "DSStorage.hh"
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4Orb.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalConstants.hh"
#include "G4SubtractionSolid.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4VisAttributes.hh"

using namespace std;

////////////////////////////////////////////////
//////        Detector Description      ////////
////////////////////////////////////////////////

// 24th April 2015
// schematics of the TPC for the LICORNE test beam detector, starting from scaling the 5K geometry (the cryostats in particular) and from some drawings from Hanguo
// OPTICS is directly derived from DS50
//
// Placements are realized with no shifts. The correct position is calculated when making the solids.
//
// GENERAL SCHEME:
// - nested single volumes, from the outside to the center: Vacuum Jacket (steel), Vacuum, Inner Cryo (ssteel).
// - the fill of the inner Volume is divided in two region at z = 1200, which correspond to the LAr-GAr interface (Inactive LAr and InactiveGar)
//   - the TPC lays inside the InactiveLAr region (Copper Rings, Reflector, TPBSide , Active LAr)
//   - the Gas Pocket lays in the upper part (reflector, Active GAr)
// - other volumes:
//   - the grid is placed 0.5 cm below the LAr surface
//   - Fused silica windows with SiPm on the outfacing surface and TPB on the innerfacing one
//   - Pseudo Argon Laryer (condensation on the top window)

// addition note on TPBSide and TPBTop
// TPBSide is a volume sorrunding the LAr Active volume on the side and bottom. TPBTop is only between
// LArLAyer (Pseudo Ar) and the top window

DSDetectorLicorneNeutron::DSDetectorLicorneNeutron(G4VPhysicalVolume *myMotherVolume) {
  fMotherVolume = myMotherVolume;
  DSLog(routine) << " Constructing Licorne Neutron Detector Geometry" << endlog;

  G4bool myCheckOverlap = DSStorage::Get()->GetCheckOverlap();
  G4ThreeVector myZeros(0., 0., 0.);

  //----------------------------------------------------------------------------
  // Get Info for neutron detector construction
  //----------------------------------------------------------------------------
  const int neutronDetectorNumber = 2;

  G4double detectorRadius = DSStorage::Get()->GetLicNDRadius();

  G4double detectorPhi1 = DSStorage::Get()->GetLicNDPhi1();
  G4double detectorPhi2 = DSStorage::Get()->GetLicNDPhi2();
  G4double detectorNuclRecEnergy = DSStorage::Get()->GetLicNDNuclRecE();

  G4double neutronEnergy = DSStorage::Get()->GetLicNeutronE();

  G4double detectorTheta;
  if (DSStorage::Get()->GetLicNdAngleByAng() == true)
    detectorTheta = DSStorage::Get()->GetLicNDTheta();
  else
    detectorTheta = AngleFromNuclearRecoilEnergy(neutronEnergy, detectorNuclRecEnergy);

  DSLog(routine) << "    with : " << endlog;
  DSLog(routine) << "       - detectorRadius= " << detectorRadius / cm << " cm" << endlog;
  DSLog(routine) << "       - detectorPhi1= " << detectorPhi1 / deg << " degree" << endlog;
  DSLog(routine) << "       - detectorPhi2= " << detectorPhi2 / deg << " degree" << endlog;
  DSLog(routine) << "       - detectorNuclRecEnergy= " << detectorNuclRecEnergy / keV << " keV" << endlog;
  DSLog(routine) << "       - neutronEnergy= " << neutronEnergy / keV << " keV" << endlog;
  DSLog(routine) << "       - detectorTheta= " << detectorTheta / deg << " degree" << endlog;

  //----------------------------------------------------------------------------
  // Set the neutron detector position in R/theta/phi
  //----------------------------------------------------------------------------
  G4ThreeVector detectorPosition[neutronDetectorNumber];

  G4double x1 = detectorRadius * cos(detectorTheta);
  G4double y1 = detectorRadius * sin(detectorTheta) * sin(detectorPhi1);
  G4double z1 = detectorRadius * sin(detectorTheta) * cos(detectorPhi1);

  detectorPosition[0].set(x1, y1, z1);

  G4double x2 = detectorRadius * cos(detectorTheta);
  G4double y2 = detectorRadius * sin(detectorTheta) * sin(detectorPhi2);
  G4double z2 = detectorRadius * sin(detectorTheta) * cos(detectorPhi2);

  detectorPosition[1].set(x2, y2, z2);

  /*  G4double detectorPhi[ neutronDetectorNumber ];
    detectorPhi[0] = detectorPhi1;
    detectorPhi[1] = detectorPhi2;
  */

  G4double neutronVesselHigh = 10 * cm;
  G4double neutronVesselThick = 0.5 * cm;

  G4double neutronActiveRadius = 10 * cm;
  G4double neutronActiveHigh = 2.5 * cm;

  // G4ThreeVector neutronSourcePos (-150. , 0., 0.);

  // G4double PMTBottomCylinderH = 11.4/2 * cm;
  //-------------------------//
  //       Outer VacCryo      //
  //-------------------------//
  G4Tubs *solidNDVesse = new G4Tubs("NeutronVessel_Solid",
                                    0.,
                                    neutronActiveRadius + neutronVesselThick,
                                    0.5 * neutronVesselHigh,
                                    0.,
                                    2 * M_PI);

  G4Tubs *fSolidNeutronActive = new G4Tubs("NeutronActive_Solid",
                                           0.,
                                           neutronActiveRadius,
                                           0.5 * neutronActiveHigh,
                                           0.,
                                           2 * M_PI);

  vector<G4Material *> BorecMat;
  BorecMat.push_back(DSMaterial::Get()->GetBorexScintillator1());
  BorecMat.push_back(DSMaterial::Get()->GetBorexScintillator2());

  // Loop over the detector number
  for (int i = 0; i < neutronDetectorNumber; i++) {
    G4LogicalVolume *fLogicNeutronLicorneVessel = new G4LogicalVolume(solidNDVesse,
                                                                      DSMaterial::Get()->GetStainlessSteel(),
                                                                      "NeutronLicorne_Logic");

    G4LogicalVolume *fLogicNeutronActiveLicorne = new G4LogicalVolume(fSolidNeutronActive,
                                                                      BorecMat[i],
                                                                      "NeutronActiveLicorne_Logic");

    new G4PVPlacement(NULL,
                      G4ThreeVector(0, 0, -0.5 * neutronVesselHigh + 0.5 * neutronActiveHigh + neutronVesselThick),
                      fLogicNeutronActiveLicorne,
                      "ActiveVolume",
                      fLogicNeutronLicorneVessel,
                      false,
                      0,
                      myCheckOverlap);

    G4RotationMatrix *myNeutronRotation = new G4RotationMatrix;

    myNeutronRotation->rotateZ(-detectorPosition[i].phi());
    myNeutronRotation->rotateY(-detectorPosition[i].theta());

    new G4PVPlacement(myNeutronRotation,
                      detectorPosition[i],
                      "NeutronPosition",
                      fLogicNeutronLicorneVessel,
                      fMotherVolume,
                      false,
                      0,
                      myCheckOverlap);
  }

  // G4PVPlacement *fPhysicNeutronLicorne[5];
  // G4PVPlacement *fPhysicNeutronActiveLicorne[5];

  // G4double neutronSourceR  = 5 * cm;
  // G4double neutronSourceH  = 5 * cm;

  /*G4double myTheta_min = -60*deg;
  G4double myTheta_step = 30*deg;

  G4double myRadiusNeutron = 100*cm ;

  G4ThreeVector myNeutronActivePos (0,0,neutronBodyH-neutronActiveH);

 for(int i=0;i<5;i++){

      G4RotationMatrix* myNeutronRotation = new G4RotationMatrix;

      G4double myTheta = myTheta_min + i * myTheta_step;
      G4double myPhi = 0*deg;

      myNeutronRotation->rotateZ( -myPhi );
      myNeutronRotation->rotateY( -(M_PI + myTheta) );

      G4ThreeVector myUnitaryVect( sin( myTheta) * cos( myPhi ), sin( myTheta ) * sin( myPhi ), cos( myTheta ) );
      G4ThreeVector myNeutronPos    = myRadiusNeutron * myUnitaryVect;

      fPhysicNeutronLicorne[i] = new G4PVPlacement(myNeutronRotation,
                               myNeutronPos,
                         "NeutronPosition",
                         fLogicNeutronLicorne,
                         fMotherVolume,
                         false,
                         0,
                         myCheckOverlap);
      fPhysicNeutronActiveLicorne[i] = new G4PVPlacement(0,
                           myNeutronActivePos,
                                 "NeutronActivePosition",
                           fLogicNeutronActiveLicorne,
                           fPhysicNeutronLicorne[i],
                           false,
                               0,
                     myCheckOverlap);


  }
*/

  /*
  G4Tubs *fSolidNeutronSource = new G4Tubs("NeutronSource_Solid",
                          0.,
                      neutronSourceR,
                      neutronSourceH,
                      0.,
                      2*M_PI);

                      */

  /*
   G4LogicalVolume   *fLogicNeutronSource = new G4LogicalVolume(fSolidNeutronSource,
                       DSMaterial::Get()->GetStainlessSteel(),
                 "NeutronSource_Logic");


   G4ThreeVector myNeutronSourcePos (-150., 0, 0);


   G4PVPlacement *fPhysicNeutronSource = new G4PVPlacement(0,
                       myNeutronSourcePos,
                 "NeutronSourcePosition",
                 fLogicNeutronSource,
                                         fMotherVolume,
                 false,
                             0,
                 myCheckOverlap);

   */

  DefineSurfaces();

  // Set the z coordinate of the LAr - GAr interface, necessary for S2 generation in DSLightX
  // DSStorage::Get()->SetLArGArBoundaryPosZ( LArGarIntefaceZ  + 1.0*um);

  // make SiPM as pe storing material
  // DSStorage::Get()->SetPMTMaterialIndex(fLogicSiPMBottom->GetMaterial()->GetIndex());
}

DSDetectorLicorneNeutron::~DSDetectorLicorneNeutron() {
  ;  // delete fMessenger;
}

void DSDetectorLicorneNeutron::DefineSurfaces() {
  ////////////////////////////////////////
  // Copied from DS DetectorDS50 - 24th april 2015
  ////////////////////////////////////////

  return;
}

G4double DSDetectorLicorneNeutron::AngleFromNuclearRecoilEnergy(G4double G4E_ni, G4double G4E_af) {
  double m_a = 38.962383 * 931.5 * 1.e6;  // Mass of 40Ar nucleus in eV
  double m_n = 1.008 * 931.5 * 1.e6;      // Mass of neutron in eV (1.008 amu * conversion factor to MeV * convert to eV)

  double E_ni = G4E_ni / eV;  // Initial neutron energy in eV
  double E_af = G4E_af / eV;  // Final nuclear recoil energy in eV
  double E_nf = E_ni - E_af;

  double v_ni = sqrt(2 * E_ni / m_n);  // initial velocity magnitude of neutron
  double v_af = sqrt(2 * E_af / m_a);  // final velocity magnitude of nucleus
  double v_nf = sqrt(2 * E_nf / m_n);  // final velocity magnitude of neutron

  double num = m_a * m_a * v_af * v_af + m_n * m_n * v_ni * v_ni - m_n * m_n * v_nf * v_nf;
  double den = 2 * m_n * m_a * v_ni * v_af;
  double phi = acos(num / den);  // angle from the neutron beam axis of the nuclear recoil

  // Transverse component
  // double theta = asin((m_a * v_af * sin(phi)) / (m_n * v_nf)) ; // angle from the neutron beam axis of the scattered neutron

  // Longitudianl component
  double theta = acos((m_n * v_ni - m_a * v_af * cos(phi)) / (m_n * v_nf));  // angle from the neutron beam axis of the scattered neutron

  G4double theta_deg = theta * (180. / 3.14159265359) * degree;

  return theta_deg;
}
