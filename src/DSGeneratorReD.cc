#include "DSGeneratorReD.hh"
#include <fstream>
#include "DSDetectorReD.hh"
#include "DSEventHandler.hh"
#include "DSGeneratorReDMessenger.hh"
#include "DSStorage.hh"
#include "G4Electron.hh"
#include "G4Event.hh"
#include "G4Gamma.hh"
#include "G4Geantino.hh"
#include "G4GeneralParticleSource.hh"
#include "G4IonTable.hh"
#include "G4MuonMinus.hh"
#include "G4MuonPlus.hh"
#include "G4Neutron.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"
#include "globals.hh"

DSGeneratorReD::DSGeneratorReD() : DSVGenerator("DSGeneratorReD") {

  fNeutronEnergy = -1.;
  fChargedEnergy = -1.;
  fChargedAngle = -1.;
  isCalculated = false;
  isCalculatedDD = false;

  fReactionType = "Li7";  //"DD" or "Li7"; default

  particleGun = new G4GeneralParticleSource();

  if (fReactionType == "Li7") {
    fBeamEnergy = 20 * MeV;
  } else if (fReactionType == "DD") {
    fBeamEnergy = 50 * keV;
  }

  fNeutronAngleThetaX = 0. * deg;
  fNeutronAnglePhiX = 0. * deg;

  fNeutronAngleThetaX = DSStorage::Get()->GetNeutronBeamThetaX();  // 18.8*deg wrt horizzontal
  fNeutronAnglePhiX = DSStorage::Get()->GetNeutronBeamPhiX();      //-12.75*deg wrt vertical

  fMainDirection = G4ThreeVector(0, 0, 1);
  fMainDirection.setTheta(pi / 2 + fNeutronAngleThetaX);
  fMainDirection.setPhi(fNeutronAnglePhiX);

  // fConeOpening = 2.*deg;
  fCosConeOpening = std::cos(fConeOpening);
  // G4cout << "CosConeOpening " << fCosConeOpening << G4endl;

  // Set the defaults
  particleGun->SetParticleDefinition(G4Neutron::Definition());
  // energy
  particleGun->GetCurrentSource()->GetEneDist()->SetMonoEnergy(1 * MeV);
  G4ThreeVector momentum(0, 0, 1);
  momentum.setTheta(pi / 2 + fNeutronAngleThetaX);
  momentum.setPhi(fNeutronAnglePhiX);
  // direction
  particleGun->GetCurrentSource()->GetAngDist()->SetParticleMomentumDirection(momentum);
  // position
  particleGun->GetCurrentSource()->GetPosDist()->SetCentreCoords(G4ThreeVector(0., 0., 0.));

  fMessenger = new DSGeneratorReDMessenger(this);
  /*
  if(!particleGun) {
    DSLog(error) << "Could not allocate G4ParticleGun! Out of memory?"<<endlog;
    DSLog(fatal) << endlog;
  }
*/

  DSLog(routine) << "ReD Generator Constructed." << endlog;
}

DSGeneratorReD::~DSGeneratorReD() {
  delete particleGun;
  delete fMessenger;
}

G4double DSGeneratorReD::GetNeutronAngle(G4ThreeVector momentum) {
  G4ThreeVector v_asseX(1, 0, 0);
  G4ThreeVector v_asseY(0, 1, 0);
  if (fReactionType == "Li7") {
    if (!isCalculated) fNeutronAngle = momentum.angle(v_asseX);
  } else if (fReactionType == "DD") {
    if (!isCalculatedDD) fNeutronAngle = momentum.angle(v_asseY);
  }
  return fNeutronAngle;
}

void DSGeneratorReD::DSGeneratePrimaries(G4Event* event)  // anEvent
{
  // Sample neutron angle
  isCalculated = false;
  isCalculatedDD = false;

  static bool firstTime = false;

  if (!firstTime) {
    // const DSDetectorReD* theDete = static_cast<const DSDetectorReD*>
    //(G4RunManager::GetRunManager()->GetUserDetectorConstruction());

    // fNeutronAngleThetaX = theDete->GetTPCAngleThetaX();
    // fNeutronAnglePhiX = theDete->GetTPCAnglePhiX();
    fMainDirection = G4ThreeVector(0, 0, 1);
    fMainDirection.setTheta(pi / 2 + fNeutronAngleThetaX);
    fMainDirection.setPhi(fNeutronAnglePhiX);
  }

  // Direction: cone of 2 deg opening
  G4double cosTheta = fCosConeOpening + (1. - fCosConeOpening) * G4UniformRand();

  G4double sinTheta = std::sqrt(1 - cosTheta * cosTheta);
  G4double phi = twopi * G4UniformRand();
  G4ThreeVector momentum(sin(phi) * sinTheta, cos(phi) * sinTheta, cosTheta);
  momentum.rotateUz(fMainDirection);

  // theta() is with respect to the z-axis, subtract 90 deg to make it with
  // respect to the beam direction == x-axis

  G4double theta = GetNeutronAngle(momentum);

  G4double fThetaArX = DSStorage::Get()->GetBeamToTPCAngleThetaX();  // 18.8*deg wrt horizzontal
  G4double fphiArX = DSStorage::Get()->GetBeamToTPCAnglePhiX();      //-12.75*deg wrt vertical
  fTargetTPCDistance = DSStorage::Get()->GetTargetToTPCDistance();   // -138.496*cm

  G4double Angolo_DD = std::acos(std::cos(fThetaArX) * std::sin(fphiArX)) / deg;
  G4double Angolo_7Li = std::acos(std::cos(fThetaArX) * std::cos(fphiArX)) / deg;

  if (!firstTime) {
    G4cout << "Reaction type: " << fReactionType << G4endl;
    G4cout << "Particle gun particle definition: " << particleGun->GetParticleDefinition()->GetParticleName() << G4endl;
    G4cout << "Generating primaries over a cone of opening: " << fConeOpening / deg << " deg" << G4endl;
    G4cout << "CosConeOpening: " << fCosConeOpening << G4endl;
    G4cout << "Solid angle: " << (1 - fCosConeOpening) * twopi << " srad" << G4endl;
    G4cout << "Beam energy: " << fBeamEnergy / MeV << " MeV" << G4endl;
    G4cout << "Target to TPC distance: " << fTargetTPCDistance / cm << " cm" << G4endl;
    G4cout << "Angle of TPC's center's vector wrt X axis (axis of 7Li's "
              "primary beam): "
           << Angolo_7Li << G4endl;
    G4cout << "Angle of TPC's center's vector wrt Y axis (axis of DD's primary "
              "beam): "
           << Angolo_DD << G4endl;
    G4cout << "Neutron Main Direction: " << fMainDirection << G4endl;
    G4cout << "Neutron Main Direction (theta): " << fMainDirection.getTheta() / deg << G4endl;
    G4cout << "Neutron Main Direction (phi): " << fMainDirection.getPhi() / deg << G4endl;
    G4cout << "Neutron Scattering Angle: " << fNeutronAngle / deg << G4endl;
    firstTime = true;
  }

  fNeutronMomentum = momentum;

  // direction
  particleGun->GetCurrentSource()->GetAngDist()->SetParticleMomentumDirection(momentum);

  // energy
  if (GetNeutronEnergy(theta) > 0) {
    particleGun->GetCurrentSource()->GetEneDist()->SetMonoEnergy(fNeutronEnergy);
    // G4cout << " Shooting energy: " << fNeutronEnergy/MeV << " MeV at " <<
    // theta/deg << " degrees" << G4endl;
  }

  particleGun->GeneratePrimaryVertex(event);  // anEvent

  /*
  DSLog(routine) << particleGun->GetParticleDefinition()->GetPDGEncoding() << '
  '
                 << particleGun->GetParticlePosition() << ' '
                 << particleGun->GetParticleMomentumDirection() << ' '
                 << particleGun->GetParticleEnergy() << ' '
                 << particleGun->GetParticleTime()
                 << endlog;
  */
  DSEventHandler::Get()->SetPDG(particleGun->GetParticleDefinition()->GetPDGEncoding());
  DSEventHandler::Get()->SetPosition(particleGun->GetParticlePosition());
  DSEventHandler::Get()->SetDirection(particleGun->GetParticleMomentumDirection());
  DSEventHandler::Get()->SetEnergy(particleGun->GetParticleEnergy());
  // MK 180303: not sure that the time is used somewhere
  //  DSEventHandler::Get()->SetTime(particleGun->GetParticleTime());
}

G4double DSGeneratorReD::GetNeutronEnergy(G4double neutronAngle) {
  if (fReactionType == "Li7") {
    if (!isCalculated) CalculateKinematics(neutronAngle);
  } else if (fReactionType == "DD") {
    if (!isCalculatedDD) CalculateKinematicsDD(neutronAngle);
  }
  return fNeutronEnergy;
}

G4double DSGeneratorReD::GetChargedParticleEnergy(G4double neutronAngle) {
  if (fReactionType == "Li7") {
    if (!isCalculated) CalculateKinematics(neutronAngle);
  } else if (fReactionType == "DD") {
    if (!isCalculatedDD) CalculateKinematicsDD(neutronAngle);
  }
  return fChargedEnergy;
}

G4double DSGeneratorReD::GetChargedParticleAngle(G4double neutronAngle) {
  if (fReactionType == "Li7") {
    if (!isCalculated) CalculateKinematics(neutronAngle);
  } else if (fReactionType == "DD") {
    if (!isCalculatedDD) CalculateKinematicsDD(neutronAngle);
  }
  return fChargedAngle;
}

void DSGeneratorReD::CalculateKinematics(G4double neutronAngle) {
  // G4cout << "In CalculateKinematics - neutronAngle/deg = " <<
  // neutronAngle/deg <<  G4endl; Inputs
  G4double beamEnergy = GetBeamEnergy();

  // Calculations for 7Li(p,n)
  G4double u = 931.49432 * MeV;      // atomic mass unit (unit of MeV/c^2)
  G4double xmbe = 7.016928 * u;      // 7Be mass
  G4double xmneu = 1.008664916 * u;  // neutron mass
  G4double xmprot = 1.007825 * u;    // proton mass
  G4double xmli = 7.01600 * u;       // 7Li mass

  /*G4double Q,TH3,B,E1,P1,P12,Etot,cos3,sin3,topolino,pippo,pluto,P3;
  G4double En,Ebe,P4,prova,param,FI4;
  */
  G4double Q = xmli + xmprot - xmneu - xmbe;  // reaction q-value
  // G4double B=xmprot/(xmprot+xmli);
  G4double Eth = -Q * ((xmli + xmprot + xmbe + xmneu) / (2. * xmprot));  // energy threshold
  if (beamEnergy < Eth) {
    G4cout << "Beam energy below threshold:" << Eth / MeV << " MeV" << G4endl;
    isCalculated = true;
    fNeutronEnergy = -1. * MeV;
    return;
  }

  G4double P1 = sqrt(2. * xmli * beamEnergy);  // 7Li impulse
  G4double P12 = P1 * P1;
  G4double Etot = beamEnergy + Q;  // total energy of outgoing particles
  G4double cos3 = std::cos(neutronAngle);
  G4double sin3 = sin(neutronAngle);

  G4double topolino = cos3 * cos3;
  G4double pippo = (4. * P12 * topolino) - (4. * (1. + (xmbe / xmneu))) * (P12 - (2. * xmbe * Etot));
  G4double pluto = 2. * (1. + (xmbe / xmneu));

  if (pippo < 0) {
    G4cout << "Kinematic condition not allowed" << G4endl;
    isCalculated = true;
    fNeutronEnergy = -1. * MeV;
    return;
  }

  G4double P3 = ((2. * P1 * cos3) + sqrt(pippo)) / pluto;  // neutron impulse
  G4double En = (P3 * P3) / (2. * xmneu);                  // neutron energy

  G4double Ebe = Etot - En;             // beryllium energy
  G4double P4 = sqrt(2. * xmbe * Ebe);  // beryllium impulse

  // G4double FI4=asin((P3/P4)*sin3);                  //beryllium recoil angle

  // Outputs
  fNeutronEnergy = (P3 * P3) / (2. * xmneu);
  // G4cout << "END: " << fNeutronEnergy/MeV << " " << P3 << " " << xmneu << "
  // MeV" << G4endl;
  fChargedEnergy = Etot - En;
  fChargedAngle = std::asin((P3 / P4) * sin3);

  isCalculated = true;
  return;
}

void DSGeneratorReD::CalculateKinematicsDD(G4double neutronAngle) {

  // Inputs
  G4double beamEnergy = GetBeamEnergy();

  // Calculations for D(d,n)3He
  G4double uma = 931494.0954 * keV;  // keV === da pdg == 931.494 0954(57) MeV/c2
  G4double mD = 2.014102 * uma;      // deuton mass
  G4double m3He = 3.016029 * uma;    // 3He mass
  G4double mN = 1.008664916 * uma;   // neutron mass

  G4double Q = (2 * mD) - (mN + m3He);
  G4double P1 = sqrt(2. * mD * beamEnergy);  // modulo impulso D incidente
  G4ThreeVector p1(P1, 0, 0);                // vettore impulso D incidente

  G4double cos3 = std::cos(neutronAngle);
  G4double sin3 = std::sin(neutronAngle);

  G4double sqrtNEnergy = (sqrt(mD * mN * beamEnergy) * cos3 + sqrt(mD * mN * beamEnergy * cos3 * cos3 + (m3He + mN) * (m3He * Q + (m3He - mD) * beamEnergy))) / (m3He + mN);
  // Eccomi!
  fNeutronEnergy = sqrtNEnergy * sqrtNEnergy;

  G4double P3 = sqrt(2. * mN * fNeutronEnergy);  // modulo impulso n
  G4ThreeVector p3(P3 * cos3, P3 * sin3, 0);     // vettore impulso n

  G4ThreeVector p4 = p1 - p3;  // vettore impulso 3He
  G4double P4 = p4.mag();      // modulo impulso 3He

  fChargedEnergy = P4 * P4 / (2. * m3He);
  fChargedAngle = p4.angle(p1);

  isCalculatedDD = true;
  return;
}

void DSGeneratorReD::SetReactionType(G4String atype) {

  fReactionType = atype;

  if (fReactionType == "Li7") G4cout << "Reaction type set to Li7+p" << G4endl;
  else if (fReactionType == "DD") {
    fBeamEnergy = 50 * keV;
    G4cout << "Reaction type set to DD" << G4endl;
    G4cout << "Beam energy reset at 50 keV" << G4endl;
  } else if (fReactionType == "geantino")
    particleGun->SetParticleDefinition(G4Geantino::Definition());
  G4cout << "Particle type set to " << particleGun->GetParticleDefinition()->GetParticleName() << G4endl;

  return;
}
