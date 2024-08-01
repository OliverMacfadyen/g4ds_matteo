// --------------------------------------------------------------------------//
/**
 * AUTHOR: A. Meregaglia
 */
// --------------------------------------------------------------------------//

#include "DSGeneratorFLUKA.hh"
#include "DSEventHandler.hh"
#include "DSIO.hh"
#include "DSLogger.hh"
#include "DSStorage.hh"
#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SPSRandomGenerator.hh"
#include "G4SystemOfUnits.hh"
//#include "DSGeneratorFLUKAMessenger.hh"
#include "DSG4DSReader.hh"
#include "G4StateManager.hh"
#include "G4ApplicationState.hh"
#include "G4ExceptionSeverity.hh"

//---------------------------------------------------------------------------//

DSGeneratorFLUKA::DSGeneratorFLUKA() : DSVGenerator("DSGeneratorFLUKA") {
  fParticleTable = G4ParticleTable::GetParticleTable();
  fFile.open(DSStorage::Get()->GetFLUKAfilename().c_str());
  DSLog(routine) << "FLUKA Generator Built"
                 << " " << DSStorage::Get()->GetFLUKAfilename() << endlog;
}
//---------------------------------------------------------------------------//

DSGeneratorFLUKA::~DSGeneratorFLUKA() {
  ;  // delete fMessenger;
}

//---------------------------------------------------------------------------//

void DSGeneratorFLUKA::DSGeneratePrimaries(G4Event* event) {

  int nparticles;

  bool IsSmearDirecitons = true ;

  if (!fFile) {
    DSLog(routine) << "Please set fluka input file name with /ds/generator/fluka_filename " << endlog;
    DSLog(fatal) << "Please set fluka input file name with /ds/generator/fluka_filename " << endlog;
  }

  fFile >> nparticles;
  if (fFile.eof()) {
     G4Exception("End of FLUKA file","DSGeneratorFLUKA", RunMustBeAborted, "End of FLUKA file");
    return;
  }


  int pdg, g4ds_pdg;
  double ene, t0, xx, yy, zz, px, py, pz, outer, inner, tpc;
  for (int i = 0; i < nparticles; ++i) {
    fFile >> pdg >> ene >> t0 >> xx >> yy >> zz >> px >> py >> pz >> outer >> inner >> tpc;
    if (fFile.eof()) {
      G4Exception("End of FLUKA file","DSGeneratorFLUKA", RunMustBeAborted, "End of FLUKA file");
      return;
    }


    double vunit = sqrt(pow(px, 2) + pow(py, 2) + pow(pz, 2));

    px /= vunit;
    py /= vunit;
    pz /= vunit;
    ene   *= GeV;
    outer *= GeV;
    inner *= GeV;
    tpc   *= GeV;
    t0    *= s ;



    fPosition = G4ThreeVector(xx * cm, yy * cm, zz * cm);
    g4ds_pdg = ConvertFlukaPDG(pdg);

    fParticle = fParticleTable->FindParticle(g4ds_pdg);

    G4ThreeVector dir = G4ThreeVector(px, py, pz);
    if (IsSmearDirecitons) dir = SmearDirection(dir, 2); // 2 degrees


    G4PrimaryParticle* particle = new G4PrimaryParticle(fParticle);
    particle->SetKineticEnergy(ene);
    //particle->SetMomentumDirection(theMomentum);
    particle->SetMomentumDirection(dir);

    G4PrimaryVertex* vertex = new G4PrimaryVertex(fPosition, t0);
    vertex->SetPrimary(particle);
    event->AddPrimaryVertex(vertex);

    DSEventHandler::Get()->SetUserFloat1(tpc / MeV);
    DSEventHandler::Get()->SetUserFloat2(inner / MeV);
    DSEventHandler::Get()->SetUserDouble(outer / MeV);
    DSEventHandler::Get()->SetUsers();

    DSEventHandler::Get()->SetDId(i);
    DSEventHandler::Get()->SetDPosition(fPosition);
    DSEventHandler::Get()->SetDDirection(dir);
    DSEventHandler::Get()->SetDPDG(g4ds_pdg);
    DSEventHandler::Get()->SetDTime(t0 / ns);
    DSEventHandler::Get()->SetDEnergy(ene / MeV);
    DSEventHandler::Get()->SetDaughters();

    DSLog(development) << pdg << " " << g4ds_pdg << " " << ene << " " << t0 << " " << xx << " " << yy << " " << zz << " " << px << " " << py << " " << pz << " Done" << endlog;
  }
}

G4ThreeVector DSGeneratorFLUKA::SmearDirection(G4ThreeVector dir, G4double sigma_alpha) {


    G4double alpha = G4RandGauss::shoot(0.0,sigma_alpha*twopi/360.);

    G4double phi = G4UniformRand()*twopi;

    G4double SinAlpha = std::sin(alpha);
    G4double CosAlpha = std::cos(alpha);
    G4double SinPhi = std::sin(phi);
    G4double CosPhi = std::cos(phi);

    G4double unit_x = SinAlpha * CosPhi;
    G4double unit_y = SinAlpha * SinPhi;
    G4double unit_z = CosAlpha;

    G4ThreeVector FacetNormal;

    FacetNormal.setX(unit_x);
    FacetNormal.setY(unit_y);
    FacetNormal.setZ(unit_z);

    G4ThreeVector tmpNormal = dir;

    FacetNormal.rotateUz(tmpNormal);

    return FacetNormal;
}

int DSGeneratorFLUKA::ConvertFlukaPDG(int pdg) {

  if (pdg == -6) return 1000020400;  // 4-HELIUM (1)
  else if (pdg == -5)
    return 1000020300;  //  3-HELIUM (1)
  else if (pdg == -4)
    return 1000010300;  //  TRITON   (1)
  else if (pdg == -3)
    return 1000010200;  //  DEUTERON (1)
  else if (pdg == -1)
    return 50;  //  OPTIPHOT
  else if (pdg == 1)
    return 2212;  //	 PROTON
  else if (pdg == 2)
    return -2212;  //	 APROTON
  else if (pdg == 3)
    return 11;  //	 ELECTRON
  else if (pdg == 4)
    return -11;  //	 POSITRON
  else if (pdg == 5)
    return 12;  //	 NEUTRIE
  else if (pdg == 6)
    return -12;  //	 ANEUTRIE
  else if (pdg == 7)
    return 22;  //	 PHOTON
  else if (pdg == 8)
    return 2112;  //	 NEUTRON
  else if (pdg == 9)
    return -2112;  //	 ANEUTRON
  else if (pdg == 10)
    return -13;  //	  MUON+
  else if (pdg == 11)
    return 13;  //	  MUON-
  else if (pdg == 12)
    return 130;  //	  KAONLONG
  else if (pdg == 13)
    return 211;  //	   PION+
  else if (pdg == 14)
    return -211;  //	   PION-
  else if (pdg == 15)
    return 321;  //	   KAON+
  else if (pdg == 16)
    return -321;  //	  KAON-
  else if (pdg == 17)
    return 3122;  //	LAMBDA
  else if (pdg == 18)
    return -3122;  //	ALAMBDA
  else if (pdg == 19)
    return 310;  //	KAONSHRT
  else if (pdg == 20)
    return 3112;  //	     SIGMA-
  else if (pdg == 21)
    return 3222;  //	  SIGMA+
  else if (pdg == 22)
    return 3212;  //	  SIGMAZER
  else if (pdg == 23)
    return 111;  //	     PIZERO
  else if (pdg == 24)
    return 311;  //	     KAONZERO
  else if (pdg == 25)
    return -311;  //	     AKAONZER
  else if (pdg == 27)
    return 14;  //	     NEUTRIM
  else if (pdg == 28)
    return -14;  //	  ANEUTRIM
  else if (pdg == 31)
    return -3222;  //	     ASIGMA-
  else if (pdg == 32)
    return -3212;  //	     ASIGMAZE
  else if (pdg == 33)
    return -3112;  //	     ASIGMA+
  else if (pdg == 34)
    return 3322;  //	     XSIZERO
  else if (pdg == 35)
    return -3322;  //	     AXSIZERO
  else if (pdg == 36)
    return 3312;  //	     XSI-
  else if (pdg == 37)
    return -3312;  //	     AXSI+
  else if (pdg == 38)
    return 3334;  //	     OMEGA-
  else if (pdg == 39)
    return -3334;  //	     AOMEGA+
  else if (pdg == 41)
    return -15;  //	  TAU+
  else if (pdg == 42)
    return 15;  //	     TAU-
  else if (pdg == 43)
    return 16;  //	     NEUTRIT
  else if (pdg == 44)
    return -16;  //	  ANEUTRIT
  else if (pdg == 45)
    return 411;  //	  D+
  else if (pdg == 46)
    return -411;  //	  D-
  else if (pdg == 47)
    return 421;  //	  D0
  else if (pdg == 48)
    return -421;  //	  D0BAR
  else if (pdg == 49)
    return 431;  //	  DS+
  else if (pdg == 50)
    return -431;  //	  DS-
  else if (pdg == 51)
    return 4122;  //	  LAMBDAC+
  else if (pdg == 52)
    return 4232;  //	  XSIC+
  else if (pdg == 53)
    return 4112;  //	  XSIC0
  else if (pdg == 54)
    return 4322;  //	  XSIPC+
  else if (pdg == 55)
    return 4312;  //	  XSIPC0
  else if (pdg == 56)
    return 4332;  //	  OMEGAC0
  else if (pdg == 57)
    return -4122;  //	  ALAMBDC-
  else if (pdg == 58)
    return -4232;  //	  AXSIC-
  else if (pdg == 59)
    return -4132;  //	  AXSIC0
  else if (pdg == 60)
    return -4322;  //	  AXSIPC-
  else if (pdg == 61)
    return -4312;  //	  AXSIPC0
  else if (pdg == 62)
    return -4332;  //	  AOMEGAC0
  else
    return 0;
}
