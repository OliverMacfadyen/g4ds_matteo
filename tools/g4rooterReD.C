// --------------------------------------------------------------------------//
/**
 * AUTHOR: D. Franco
 * CONTACT: dfranco@in2p3.fr
 *
 * Generate a root file reading the binary file from the g4ds output
 *
 * Modified by M. Kuss (Michael.Kuss@pi.infn.it) for the ReD TPC Naples setup.
 *
 * Functions S1quench and QScintillator, written for Licorne, have to be
 * revisited.  Are the S1 quenching (probably) and light yield (probably not)
 * the same at Naples?
 *
 * Also check f90TPC, recoEneNeu..., etc.
 *
 */
// --------------------------------------------------------------------------//
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TObject.h"
#include "TRandom.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TH1.h"
#include "TRandom.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TMath.h"
#include "string.h"
#include "TMinuit.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "algorithm"
#include "../include/DSEventStructure.hh"

using namespace std;

const int MAXDAUGHTER = 70000;
const int MAXDEPOSIT  = 10000;
const int MAXNPH      = 100000;
const int MAXNPE      = 2000000;
const int MAXUSER     = 1000;
// 8 x 3in + 1 x 5in, ordered by id assigned during geometry construction
// check if # of LScis could be written to file header
const int MAXLSC      = 9;
// id's of LScis 65 - 72
//const int LSC0        = 65;
// after Simone's changes the first id is 73
// after AcrylicDART: 74
// probably also this could be written to the file header
const int LSC0        = 74;

// constants
Double_t keV = 1.0;
Double_t MeV = 1E3 * keV;
Double_t GeV = 1E3 * MeV;
Double_t cm = 1.0;
Double_t mm = 0.1 * cm;
Double_t n_m    =   939.5653451     * MeV;
Double_t Ar40_m =    37.21552324198 * GeV;

vector<Float_t> vpdf, vpdf_time, vtheta, vlength, vjitter_prob, vjitter_time;

HeaderStructure                theHeader;
EventStructureDiskFormat       theEvent;
vector<DepositStructure>       theDeposits;
vector<DaughterStructure>      theDaughters;
vector<UserStructure>          theUsers ;
vector<PhotonStructure>        thePhotons;
vector<PhotoElectronStructure> thePhotoElectrons ;
vector<PhotoElectronStructure> theVetoPhotoElectrons ;
vector<PhotoElectronStructure> theMuPhotoElectrons ;

struct cmp_photoelectron{
  bool operator() (const PhotoElectronStructure& a, const PhotoElectronStructure& b) {
    return a.Time<b.Time ;
  }
};
struct cmp_photon{
  bool operator() (const PhotonStructure& a, const PhotonStructure& b) {
    return a.Time<b.Time ;
  }
};
struct cmp_deposit{
  bool operator() (const DepositStructure& a, const DepositStructure& b) {
    return a.Time<b.Time ;
  }
};

bool _readHeader (ifstream *file) {
  int event_size;
  int event_size2;
  file->read( (char *)(&event_size),  sizeof(int) );
  file->read( (char *)(&theHeader),   sizeof(HeaderStructure) );
  file->read( (char *)(&event_size2), sizeof(int) );
  if ( !(*file) ) return false;
  return true;
}

DepositStructure _readDeposit (ifstream *file) {
  DepositStructure theDeposit ;
  file->read ((char *)(&theDeposit), sizeof( DepositStructure) );
  return theDeposit;
}

DaughterStructure _readDaughter (ifstream *file) {
  DaughterStructure theDaughter ;
  file->read ((char *)(&theDaughter), sizeof( DaughterStructure) );
  return theDaughter;
}

UserStructure _readUser (ifstream *file) {
  UserStructure theUser ;
  file->read ((char *)(&theUser), sizeof( UserStructure) );
  return theUser;
}

PhotonStructure _readPhoton (ifstream *file) {
  PhotonStructure thePhoton ;
  file->read ((char *)(&thePhoton), sizeof( PhotonStructure) );
  return thePhoton;
}

PhotoElectronStructure _readPhotoElectron (ifstream *file) {
  PhotoElectronStructure thePhotoElectron ;
  file->read ((char *)(&thePhotoElectron), sizeof( PhotoElectronStructure) );
  return thePhotoElectron;
}

PhotoElectronStructure _readVetoPhotoElectron (ifstream *file) {
  PhotoElectronStructure theVetoPhotoElectron ;
  file->read ((char *)(&theVetoPhotoElectron), sizeof( PhotoElectronStructure) );
  return theVetoPhotoElectron;
}

PhotoElectronStructure _readMuPhotoElectron (ifstream *file) {
  PhotoElectronStructure theMuPhotoElectron ;
  file->read ((char *)(&theMuPhotoElectron), sizeof( PhotoElectronStructure) );
  return theMuPhotoElectron;
}



bool _readEvent (ifstream *file) {
  theDeposits.clear();
  theDaughters.clear();
  theUsers.clear() ;
  thePhotons.clear();
  thePhotoElectrons.clear() ;
  theVetoPhotoElectrons.clear() ;
  theMuPhotoElectrons.clear() ;

  if (!(*file)) return false;

  int event_size;
  int event_size2;
  file->read ((char *)(&event_size), sizeof (int));
  file->read ((char *)(&theEvent), sizeof( EventStructureDiskFormat)  );
  if(file->eof()) return false ;

  if(theEvent.NDaughters > MAXDAUGHTER) { cout << "Fatal: NDaughters > MAXDAUGHTER : " << theEvent.NDaughters << " > " << MAXDAUGHTER << endl ; exit(0) ;}
  if(theEvent.NDeposits  > MAXDEPOSIT ) { cout << "Fatal: NDeposits = " << theEvent.NDeposits << " > MAXDEPOSIT" << endl ; exit(0) ;}
  if(theEvent.NUsers     > MAXUSER)     { cout << "Fatal: NUsers > MAXUSER" << endl ; exit(0) ;}
  if(theEvent.NPH        > MAXNPH)      { cout << "Fatal: NPH > MAXNPH" << endl ; exit(0) ;}
  if(theEvent.NPE        > MAXNPE)      { cout << "Fatal: NPE > MAXNPE" << endl ; exit(0) ;}
  if(theEvent.VetoNPE    > MAXNPE)      { cout << "Fatal: VetoNPE > MAXNPE" << endl ; exit(0) ;}
  if(theEvent.MuNPE      > MAXNPE)      { cout << "Fatal: MuNPE > MAXNPE" << endl ; exit(0) ;}


  for(int i=0; i<theEvent.NDaughters; i++) theDaughters.push_back(_readDaughter(file));
  for(int i=0; i<theEvent.NDeposits; i++)  theDeposits.push_back(_readDeposit(file));
  for(int i=0; i<theEvent.NUsers; i++)     theUsers.push_back(_readUser(file));
  for(int i=0; i<theEvent.NPH; i++)        thePhotons.push_back(_readPhoton(file));
  for(int i=0; i<theEvent.NPE; i++)        thePhotoElectrons.push_back(_readPhotoElectron(file));
  for(int i=0; i<theEvent.VetoNPE; i++)    theVetoPhotoElectrons .push_back(_readMuPhotoElectron(file));
  for(int i=0; i<theEvent.MuNPE; i++)      theMuPhotoElectrons.push_back(_readMuPhotoElectron(file));
  file->read ((char *)(&event_size2), sizeof (int));
  if(file->eof()) return false ;
  if(event_size != event_size2) return false ;
  //  cout << "Problem at event " << theEvent.EventID << endl;
  //}


  return true;

}

bool _skipEvent (ifstream *file) {
  int event_size;
  int event_size2;
  file->read ((char *)(&event_size), sizeof (int));
  char * buffer = new char [event_size];
  file->read (buffer,event_size);
  file->read ((char *)(&event_size2), sizeof (int));
  if(file->eof()) return false ;
  return true;
}



///quenching for TPC
double S1quench(double ene, double alpha) {

  double W = 19.5e-3; // keV
  double epsilon = 1;

  double p0=  0.292032;
  double p1=  1.37809;
  double p2= -0.155008;
  double p3= -0.00146758;
  double p4=  0.916802;
  double p5=  0.621216;
  double myRecoProb = p0 * ( 1 - p1 * exp(p2*ene) ) * TMath::Exp( p3 * pow(ene,p4) ) + p5;

  double myNumQuanta   = ene / W;
  double myNumIons     = myNumQuanta / ( 1 + alpha );
  double myNumExcitons = myNumQuanta - myNumIons;

  double myphotons = myNumExcitons + myNumIons * myRecoProb;
  double out = myphotons * W / ene;

  return out;
}

double QScintillator(double ene, int pdg) {
  if(pdg < 100) return ene ;
  else if(pdg <4000) return ene/3. ;
  else return ene/10 ;

  return 0;
}

Double_t sqr(Double_t x) {
  return x * x;
}

Double_t Erec(Double_t ene0, TVector3 gun_pos, TVector3 tpc_pos, TVector3 LSci_pos) {
  // taken from http://skisickness.com/2010/04/25/
  // Reaction X1(X2,X3)X4 , i.e. 40Ar(n,n')40Ar'
  Double_t x2_T = ene0;
  Double_t x2_M = n_m;

  // the projectile
  Double_t x2_E = x2_M + x2_T;
  Double_t x2_P = TMath::Sqrt( sqr(x2_E) - sqr(x2_M) );
  TVector3 x2_Vect_Unit = (tpc_pos-gun_pos).Unit();
  TLorentzVector X2(x2_Vect_Unit*x2_P, x2_E);
  //  std::cout << "---------------------------------------------------------------- X2 (neutron)" << std::endl;
  //  X2.Print();
  //  std::cout << "M " << X2.M()/MeV << " MeV  M2 " << X2.M2()/MeV/MeV << " MeV^2  P " << X2.P()/MeV << " MeV  Beta " << X2.Beta() << "  E " << X2.E()/MeV << " MeV  E-M " << (X2.E()-X2.M())/MeV << " MeV" << std::endl;

  // the target
  TLorentzVector X1(TVector3(), Ar40_m);
  //  std::cout << "---------------------------------------------------------------- X1 (40Ar)" << std::endl;
  //  X1.Print();
  //  std::cout << "M " << X1.M()/GeV << " GeV  M2 " << X1.M2()/GeV/GeV << " GeV^2  P " << X1.P()/keV << " keV  Beta " << X1.Beta() << "  E " << X1.E()/GeV << " GeV  E-M " << (X1.E()-X1.M())/keV << " keV" << std::endl;

  // the cm system
  TLorentzVector cm = X2 + X1;
  Double_t s = cm.M2();
  Double_t p_cm = TMath::Sqrt( ( sqr(s-X1.M2()-X2.M2()) - 4.0*X1.M2()*X2.M2() ) / 4.0 / s );
  Double_t p_cm2 = sqr(p_cm);
  Double_t rap_cm = TMath::Log( ( p_cm + TMath::Sqrt(X1.M2()+p_cm2) ) / X1.M() );
  /*
  std::cout << "---------------------------------------------------------------- cm" << std::endl;
  cm.Print();
  std::cout << "M " << cm.M()/GeV << " GeV  M2 " << cm.M2()/GeV/GeV << " GeV^2  P " << cm.P()/MeV << " MeV  Beta " << cm.Beta() << "  E " << cm.E()/GeV << " GeV  E-M " << (cm.E()-cm.M())/MeV << " MeV" << std::endl;
  std::cout << " particle momentum " << 1E-3*p_cm << " rapidity " << rap_cm << std::endl;
  */

  // the ejectile
  //  std::cout << "---------------------------------------------------------------- X3 (neutron)" << std::endl;
  TVector3 x3_Vect_Unit = (LSci_pos-tpc_pos).Unit();
  Double_t th3 = x2_Vect_Unit.Angle(x3_Vect_Unit);
  //  std::cout << "scattering Angle " << TMath::RadToDeg() * th3 << " deg" << std::endl;
  // of the ejectile we know only it's direction, and it's mass.  We have to reconstruct the 4mom vector.
  Double_t x3_M2 = X2.M2();
  Double_t sinHrap = TMath::SinH(rap_cm);
  Double_t sinHrap2 = sqr(sinHrap);
  Double_t x3_P = ( TMath::Sqrt(x3_M2+p_cm2) * TMath::Cos(th3) * sinHrap + TMath::CosH(rap_cm) * TMath::Sqrt(p_cm2-x3_M2*sqr(th3)*sinHrap2) ) / ( 1.0 + sqr(th3) + sinHrap2 );
  //  std::cout << "ejectile momentum " << x3_P/MeV << " MeV" << std::endl;
  TVector3 x3_Vect = x3_Vect_Unit * x3_P;
  Double_t x3_E = TMath::Sqrt(x3_M2+x3_Vect.Mag2());
  TLorentzVector X3(x3_Vect, x3_E);
  //  X3.Print();
  // note that for M() we have rounding errors.  M is Sqrt(E*E-P*P), with E^2 and P^2 differing by about 7 orders of magnitude!
  //  std::cout << "M " << X3.M()/MeV << " MeV  M2 " << X3.M2()/MeV/MeV << " MeV^2  P " << X3.P()/MeV << " MeV  Beta " << X3.Beta() << "  E " << X3.E()/MeV << " MeV  E-M " << (X3.E()-X3.M())/MeV << " MeV" << std::endl;

  // the recoil
  TLorentzVector X4 = cm - X3;
  /*
  std::cout << "---------------------------------------------------------------- X4 (40Ar)" << std::endl;
  X4.Print();
  std::cout << "M " << X4.M()/GeV << " GeV  M2 " << X4.M2()/GeV/GeV << " GeV^2  P " << X4.P()/MeV << " MeV  Beta " << X4.Beta() << "  E " << X4.E()/GeV << " GeV  E-M " << (X4.E()-X4.M())/keV << " keV" << std::endl;
  */

  Double_t X4_T = X4.E() - X4.M();
  return X4_T;
}

//--------------------------------------------
//                    main
//--------------------------------------------
int main (int argc, char *argv[]) {

  if ( argc == 1 || ( argc > 1 && !string(argv[1]).find("help") ) ) {
  //if(file.find("help") < 10) {
    cout << "Usage: " << argv[0] << " [FILE] [OPTIONS] [OUTPUT]" <<endl ;
    cout <<endl ;
    cout << " Options: " << endl ;
    cout << " nevents=N:  max number (N) of events to process (default: 1G)" <<endl ;
    cout << " TPC:        save only events with a signal in the TPC" <<endl ;
    cout << " skipN=N:    skip the first N events (default: 0)" <<endl ;
    cout << " kB=xxx:     set the Birks parameter (default: 0.012 cm/MeV)" <<endl ;
    cout << " LY=xxx:     set the light yield (default: 500 p.e./MeV)" <<endl ;
    cout << endl;
    cout << " Output: " << endl ;
    cout << " filename.root (default: FILE.root)" << endl ;
    cout << endl;
    cout << "Version v783.r4" << endl ;
    cout << "..bye ;) DF  (dfranco@in2p3.fr)" << endl ;
    return 0 ;
  }

  string file = argv[1] ;

  if(file.find(".fil") == string::npos) {
    cout << "file " <<  file << " not found.... Bye!" << endl;
    cout << "ps: the input file needs the .fil extension"  << endl;
    return 0;
  }
  string rootfile = file;
  rootfile.replace(file.find(".fil"),4,".root");

  int nevents     = 1000000000;
  int skipNevents = 0;
  Float_t kB        = 0.012;
  Float_t LY        = 500.;
  bool  TPCstore  = 0;
  bool verbose = 0;
  bool debug = 0;
  int loop=2;
  while(argv[loop]) {
    string argument = argv[loop];
    if(!argument.find("nevents=")) {
      argument.erase(0,8);
      nevents = (int)  atoi(argument.c_str()) ;
      cout << "   max number of events to process: " << nevents << endl;
    }
    if(!argument.find("skipN=")) {
      argument.erase(0,6);
      skipNevents = atoi(argument.c_str()) ;
      cout << "   events to skip: " << skipNevents << endl;
    }
    if(!argument.find("kB=")) {
      argument.erase(0,3);
      kB = atof(argument.c_str()) ;
      cout << "   kB: " << kB << " cm/MeV" << endl;
    }
    if(!argument.find("LY=")) {
      argument.erase(0,3);
      kB = atof(argument.c_str()) ;
      cout << "   LY: " << LY << " p.e./MeV" << endl;
    }
    if(!argument.find("TPC")) {
      TPCstore = true ;
      cout << "   Save only events with a signal in the TPC" << endl;
    }
    if(!argument.find("verbose")) {
      verbose = true;
      cout << "   Print LSci deposit information" << endl;
    }
    if(!argument.find("debug")) {
      debug = true;
      cout << "   Print Event deposit information" << endl;
    }
    size_t found;
    found=argument.rfind(".root");
    if(found!=string::npos) {
      rootfile = argument;
      cout << "   output root file: " << argument << endl ;
    }

    loop++;
  }

  TFile *ff = new TFile(rootfile.c_str(),"recreate");

  // event extra variables
  Float_t radius = 0;

  // daughter variables
  int    Did[MAXDAUGHTER],
         Dpdg[MAXDAUGHTER],
         Dpid[MAXDAUGHTER],
         Dprocess[MAXDAUGHTER],
         Dmat[MAXDAUGHTER];
  Float_t  Dene[MAXDAUGHTER], Dx[MAXDAUGHTER], Dy[MAXDAUGHTER], Dz[MAXDAUGHTER],
         Dr[MAXDAUGHTER], Dpx[MAXDAUGHTER], Dpy[MAXDAUGHTER], Dpz[MAXDAUGHTER];
  double Dtime[MAXDAUGHTER];

  // deposit variables
  int    *dep_pdg   = new int[MAXDEPOSIT];
  int    *dep_mat   = new int[MAXDEPOSIT];
  Float_t  *dep_ene   = new Float_t[MAXDEPOSIT];
  Float_t  *dep_x     = new Float_t[MAXDEPOSIT];
  Float_t  *dep_y     = new Float_t[MAXDEPOSIT];
  Float_t  *dep_z     = new Float_t[MAXDEPOSIT];
  Float_t  *dep_r     = new Float_t[MAXDEPOSIT];
  Float_t  *dep_step  = new Float_t[MAXDEPOSIT];
  double *dep_time  = new double[MAXDEPOSIT];

  // user variables
  int    INT1[MAXUSER], INT2[MAXUSER];
  Float_t  FLOAT1[MAXUSER], FLOAT2[MAXUSER];
  double DOUBLE[MAXUSER];

  // pe variables
  int    *pe_pmt = new int[MAXNPE];
  double *pe_time= new double[MAXNPE];
  int NPES1, NPES2 ;

  // veto pe variables
  int    *veto_pe_pmt  = new int[MAXNPE];
  double *veto_pe_time = new double[MAXNPE];

  // mu pe variables
  int    *mu_pe_pmt  = new int[MAXNPE];
  double *mu_pe_time = new double[MAXNPE];

  // ph variables
  int    ph_volume[MAXNPH], ph_pid[MAXNPH];
  Float_t  ph_wl[MAXNPH], ph_x[MAXNPH], ph_y[MAXNPH], ph_z[MAXNPH];
  double ph_time[MAXNPH];

  Float_t s1ene,s2ene, veto_visene, mu_visene, qene, qnpe, tpcene, vetoene, muene ;

  int LArIndex;

  // ReD variables (adapted from Licorne)
  // Strange: I can't define these vars in front of the Branch statements. It triggers seg faults at run time!
  Double_t  timeTPC;
  Float_t   eneTPC;
  Float_t   f90TPC;
  Int_t     ndepTPC;
  Int_t     nAr40;
  Float_t   Ar40eneCalc;
  Float_t   Ar40ene;
  Float_t   Ar40x;
  Float_t   Ar40y;
  Float_t   Ar40z;
  Float_t   Ar40px;
  Float_t   Ar40py;
  Float_t   Ar40pz;
  Int_t*    ndepLSc   = new Int_t[MAXLSC];
  Double_t* timeLSc   = new Double_t[MAXLSC];
  Float_t*  x0LSc     = new Float_t[MAXLSC];
  Float_t*  y0LSc     = new Float_t[MAXLSC];
  Float_t*  z0LSc     = new Float_t[MAXLSC];
  Float_t*  xLSc      = new Float_t[MAXLSC];
  Float_t*  yLSc      = new Float_t[MAXLSC];
  Float_t*  zLSc      = new Float_t[MAXLSC];
  Float_t*  xeneLSc   = new Float_t[MAXLSC];
  Float_t*  yeneLSc   = new Float_t[MAXLSC];
  Float_t*  zeneLSc   = new Float_t[MAXLSC];
  Float_t*  eneLScRaw = new Float_t[MAXLSC];
  Float_t*  eneLSc    = new Float_t[MAXLSC];
  Float_t*  eneLScTOF = new Float_t[MAXLSC];
  Int_t     nTrgLSc;
  TRandom3* rand = new TRandom3(0);
  const Float_t sx = 11.95 * mm / TMath::Sqrt(12.0);
  const Float_t sy =  5.95 * mm / TMath::Sqrt(12.0);
  const Float_t sz =  1.0  * mm;

  TTree *dstree = new TTree("dstree","The G4DS Root Tree");
  dstree->SetMaxVirtualSize(100000);
  dstree->SetMaxVirtualSize(10000000000);

  dstree->Branch("ev",             &theEvent.EventID,         "ev/I");
  dstree->Branch("pdg",            &theEvent.PDG,             "pdg/I");
  dstree->Branch("ene0",           &theEvent.Energy,          "ene0/F");
  dstree->Branch("s1ene",          &theEvent.S1Energy,        "s1ene/F");
  dstree->Branch("s2ene",          &theEvent.S2Energy,        "s2ene/F");
  dstree->Branch("veto_visene",    &theEvent.VetoVisEnergy,   "veto_visene/F");
  dstree->Branch("mu_visene",      &theEvent.MuVisEnergy,     "mu_visene/F");

  dstree->Branch("tpcene",         &theEvent.TPCDepEnergy,   "tpcene/F");
  dstree->Branch("vetoene",        &theEvent.VetoDepEnergy, "vetoene/F");
  dstree->Branch("muene",          &theEvent.MuDepEnergy,     "muene/F");
  dstree->Branch("x",              &theEvent.Position[0],     "x/F");
  dstree->Branch("y",              &theEvent.Position[1],     "y/F");
  dstree->Branch("z",              &theEvent.Position[2],     "z/F");
  dstree->Branch("r",              &radius,		      "radius/F");
  dstree->Branch("px",             &theEvent.Direction[0],    "px/F");
  dstree->Branch("py",             &theEvent.Direction[1],    "py/F");
  dstree->Branch("pz",             &theEvent.Direction[2],    "pz/F");
  dstree->Branch("bx",             &theEvent.CenterOfMass[0], "bx/F");
  dstree->Branch("by",             &theEvent.CenterOfMass[1], "by/F");
  dstree->Branch("bz",             &theEvent.CenterOfMass[2], "bz/F");

  dstree->Branch("npeS1",          &NPES1 ,	              "npeS1/I");
  dstree->Branch("npeS2",          &NPES2 ,        	      "npeS2/I");
  dstree->Branch("npe",            &theEvent.NPE ,	      "npe/I");
  dstree->Branch("munpe" ,         &theEvent.MuNPE ,	      "munpe/I");
  dstree->Branch("vnpe",           &theEvent.VetoNPE ,        "vnpe/I");
  dstree->Branch("nph",            &theEvent.NPH,	      "nph/I");
  dstree->Branch("ndaughters",     &theEvent.NDaughters,      "ndaughters/I");
  dstree->Branch("ndeposits" ,     &theEvent.NDeposits,       "ndeposits/I");
  dstree->Branch("nusers",         &theEvent.NUsers,	      "nusers/I");
  dstree->Branch("dau_id",         Did,                  "Did[ndaughters]/I");
  dstree->Branch("dau_pdg",        Dpdg,                 "Dpdg[ndaughters]/I");
  dstree->Branch("dau_pid",        Dpid,                 "Dpid[ndaughters]/I");
  dstree->Branch("dau_process",    Dprocess,             "Dprocess[ndaughters]/I");
  dstree->Branch("dau_time",       Dtime,                "Dtime[ndaughters]/D");
  dstree->Branch("dau_ene",        Dene,                 "Dene[ndaughters]/F");
  dstree->Branch("dau_x",          Dx,                   "Dx[ndaughters]/F");
  dstree->Branch("dau_y",          Dy,                   "Dy[ndaughters]/F");
  dstree->Branch("dau_z",          Dz,                   "Dz[ndaughters]/F") ;
  dstree->Branch("dau_r",          Dr,                   "Dr[ndaughters]/F");
  dstree->Branch("dau_px",         Dpx,                  "Dpx[ndaughters]/F");
  dstree->Branch("dau_py",         Dpy,                  "Dpy[ndaughters]/F") ;
  dstree->Branch("dau_pz",         Dpz,                  "Dpz[ndaughters]/F");
  dstree->Branch("dep_pdg",     dep_pdg,                  "dep_pdg[ndeposits]/I");
  dstree->Branch("dep_mat",     dep_mat,                  "dep_mat[ndeposits]/I");
  dstree->Branch("dep_time",    dep_time,                 "dep_time[ndeposits]/D");
  dstree->Branch("dep_ene",     dep_ene,                  "dep_ene[ndeposits]/F")  ;
  dstree->Branch("dep_step",    dep_step,                 "dep_step[ndeposits]/F")  ;
  dstree->Branch("dep_x",       dep_x,                    "dep_x[ndeposits]/F");
  dstree->Branch("dep_y",       dep_y,                    "dep_y[ndeposits]/F");
  dstree->Branch("dep_z",       dep_z,                    "dep_z[ndeposits]/F") ;
  dstree->Branch("dep_r",       dep_r,                    "dep_r[ndeposits]/F") ;

  dstree->Branch("userint1",    INT1,                     "int1[nusers]/I");
  dstree->Branch("userint2",    INT2,                     "int2[nusers]/I");
  dstree->Branch("userfloat1",  FLOAT1,                   "float1[nusers]/F")  ;
  dstree->Branch("userfloat2",  FLOAT2,                   "float2[nusers]/F");
  dstree->Branch("userdouble0", DOUBLE,                   "double0[nusers]/D");

  dstree->Branch("pe_time",     pe_time,                  "pe_time[npe]/D");
  dstree->Branch("pe_pmt",      pe_pmt,                   "pe_pmt[npe]/I");

  dstree->Branch("vpe_time",  veto_pe_time,               "veto_pe_time[vnpe]/D");
  dstree->Branch("vpe_pmt",   veto_pe_pmt,                "veto_pe_pmt[vnpe]/I");

  dstree->Branch("mupe_time",   mu_pe_time,               "mu_pe_time[munpe]/D");
  dstree->Branch("mupe_pmt",    mu_pe_pmt,                "mu_pe_pmt[munpe]/I");

  dstree->Branch("ph_volume",    ph_volume,               "ph_volume[nph]/I");
  dstree->Branch("ph_pid",       ph_pid,                  "ph_pid[nph]/I");
  dstree->Branch("ph_wl",        ph_wl,                   "ph_wl[nph]/F");
  dstree->Branch("ph_x",         ph_x,                    "ph_x[nph]/F");
  dstree->Branch("ph_y",         ph_y,                    "ph_y[nph]/F");
  dstree->Branch("ph_z",         ph_z,                    "ph_z[nph]/F");
  dstree->Branch("ph_time",      ph_time,                 "ph_time[nph]/D");

  // ReD (adapted from Licorne)
  dstree->Branch("ndepTPC",    &ndepTPC,   "ndepTPC/I");
  dstree->Branch("timeTPC",    &timeTPC,   "timeTPC/D");       // TPC time
  dstree->Branch("eneTPC",     &eneTPC,    "eneTPC/F");        // ene TPC
  dstree->Branch("f90TPC",     &f90TPC,    "f90TPC/F");        // f90 TPC, calculated from dep_elec and dep_nucl

  // Ar40 recoil
  dstree->Branch("nAr40",       &nAr40,       "nAr40/I");       // number of Ar40 recoils
  dstree->Branch("Ar40eneCalc", &Ar40eneCalc, "Ar40eneCalc/F"); // Ar40 recoil calculated from n kinematics
  dstree->Branch("Ar40ene",     &Ar40ene,     "Ar40ene/F");     // Ar40 recoil energy (only last recoil ... if more than one) is saved)
  dstree->Branch("Ar40x",       &Ar40x,       "Ar40x/F");       // Ar40 recoil position
  dstree->Branch("Ar40y",       &Ar40y,       "Ar40y/F");
  dstree->Branch("Ar40z",       &Ar40z,       "Ar40z/F");
  // the momentum dir has to be filled in g4ds
  dstree->Branch("Ar40px",      &Ar40px,      "Ar40px/F");      // Ar40 recoil momentum
  dstree->Branch("Ar40py",      &Ar40py,      "Ar40py/F");
  dstree->Branch("Ar40pz",      &Ar40pz,      "Ar40pz/F");
 
  // MK: again, read MAXLSc, don't hardcode here
  dstree->Branch("ndepLSc",     ndepLSc,   "ndepLSc[9]/I");    // number of hits in neutron detector
  dstree->Branch("timeLSc",     timeLSc,   "timeLSc[9]/D");    // time of flight from TPC
  dstree->Branch("x0LSc",       x0LSc,     "x0LSc[9]/F");      // position of first hit
  dstree->Branch("y0LSc",       y0LSc,     "y0LSc[9]/F");
  dstree->Branch("z0LSc",       z0LSc,     "z0LSc[9]/F");
  dstree->Branch("xLSc",        xLSc,      "xLSc[9]/F");       // position of energy centroid
  dstree->Branch("yLSc",        yLSc,      "yLSc[9]/F");
  dstree->Branch("zLSc",        zLSc,      "zLSc[9]/F");
  dstree->Branch("eneLScRaw",   eneLScRaw, "eneLScRaw[9]/F");  // energy
  dstree->Branch("eneLSc",      eneLSc,    "eneLSc[9]/F");     // energy, parsed through QScintillator
  dstree->Branch("eneLScTOF",   eneLScTOF, "eneLScTOF[9]/F");  // neutron energy calculated from ToF
  dstree->Branch("nTrgLSc",     &nTrgLSc,  "nTrgLSc/I");       // number of neutron detectors with energy deposits

  // Open the binary file
  ifstream *_bin_fstream;
  _bin_fstream = new ifstream (file.c_str(), std::ios::binary);
  cout << endl ;
  cout << "Binary File: " << file << endl;
  cout << endl ;
  if (!(*_bin_fstream)) {
    std::cerr << "Cannot open file. Exiting...\n";
    exit(1);
  }

  Int_t counter = 0;

  // Read Header
  _readHeader(_bin_fstream);
  // dump the entire header
  std::cout << "theHeader.Events: "            << theHeader.Events << std::endl;       
  std::cout << "theHeader.Run: "               << theHeader.Run << std::endl;       
  std::cout << "theHeader.PDG: "               << theHeader.PDG << std::endl;        
  std::cout << "theHeader.LArIndex: "          << theHeader.LArIndex << std::endl;        
  std::cout << "theHeader.ScintillatorIndex: " << theHeader.ScintillatorIndex << std::endl;    
  std::cout << "theHeader.Rate: "              << theHeader.Rate << std::endl;       
  std::cout << "theHeader.DetectorFlag: "      << theHeader.DetectorFlag << std::endl; 

  LArIndex = theHeader.LArIndex;

  // try the method from g4licorne_pulsed (doesn't work for ReD)
  size_t found_cm = file.find("cm");
  std::cout << "found_cm: " << found_cm << std::endl;
  string dist = file.substr(0,found_cm);
  std::cout << "dist: " << dist << std::endl;
  dist.erase(0,dist.find("_")+1);
  std::cout << "dist: " << dist << std::endl;

  // Loop over the events
  for(Int_t _i=0; _i<nevents; _i++) {

    // Skip Events
    if(_i < skipNevents) {
      _skipEvent(_bin_fstream);
      continue;
    }
    // Read Event
    if(!_readEvent(_bin_fstream)) break ;

    // Close the binary file if the end is reached
    if(_bin_fstream->eof()) break ;

    // Print counter
    if (!(counter % 10000) && counter>0) std::cout << counter <<" processed events " << " (event id = " <<  theEvent.EventID << ")"<< std::endl;
    counter++ ;

    // sort photon and photoelectron vectors by time
    std::sort(theDeposits.begin(),theDeposits.end(), cmp_deposit());
    std::sort(thePhotons.begin(),thePhotons.end(), cmp_photon());
    std::sort(thePhotoElectrons.begin(),thePhotoElectrons.end(), cmp_photoelectron());
    std::sort(theVetoPhotoElectrons.begin(),theVetoPhotoElectrons.end(), cmp_photoelectron());
    std::sort(theMuPhotoElectrons.begin(),theMuPhotoElectrons.end(), cmp_photoelectron());

    // Fill daughter variables
    for(Int_t i=0;i<theEvent.NDaughters;++i) {
      Did[i]      = theDaughters[i].Id ;
      Dpdg[i]     = theDaughters[i].PDG ;
      Dpid[i]     = theDaughters[i].PID ;
      Dprocess[i] = theDaughters[i].Process ;
      Dtime[i]    = theDaughters[i].Time ;
      Dene[i]     = theDaughters[i].Energy ;
      Dx[i]       = theDaughters[i].Position[0] ;
      Dy[i]       = theDaughters[i].Position[1] ;
      Dz[i]       = theDaughters[i].Position[2] ;
      Dpx[i]      = theDaughters[i].Direction[0] ;
      Dpy[i]      = theDaughters[i].Direction[1] ;
      Dpz[i]      = theDaughters[i].Direction[2] ;
      Dr[i]       = sqrt(pow(Dx[i],2)+pow(Dy[i],2)+pow(Dz[i],2));
    }

    // ReD (adapted from Licorne)

    const Float_t neutronDetectorRadius = 80.0;
    bool firstTPC = true;
    bool firstLSc[MAXLSC];

    ndepTPC = 0;
    eneTPC  = 0;
    nAr40   = 0;
    nTrgLSc = 0;
    //    Double_t timeLScMin[MAXLSC];

    for ( Int_t i=0; i<MAXLSC; ++i ) {
      firstLSc[i] = true;
      ndepLSc[i] = 0;
      eneLScRaw[i] = 0;
      eneLSc[i] = 0;
      eneLScTOF[i] = -1;
      timeLSc[i] = -1;
      x0LSc[i] = -1;
      y0LSc[i] = -1;
      z0LSc[i] = -1;
      xLSc[i] = -1;
      yLSc[i] = -1;
      zLSc[i] = -1;
      xeneLSc[i] = 0;
      yeneLSc[i] = 0;
      zeneLSc[i] = 0;
    }
    f90TPC  = -1;
    timeTPC = -1;

    Float_t qarAr   = 0.25;             // quenching of Ar in LAr
    Float_t depElec = 0;
    Float_t depNucl = 0;

    // violent debugging
    if ( debug ) {
      TString name = TString::Format("events/%07d.log", theEvent.EventID);
      system("mkdir -p events");
      FILE* out = fopen(name, "w");
      fprintf(out, "%d\n", theEvent.EventID);
      for ( Int_t i=0; i<theEvent.NDeposits; ++i ) {
        fprintf(out, "%10d  %2d  %8.3f keV  ( %8.2f , %8.2f, %8.2f ) cm  %7.3f ns  %f cm\n",  theDeposits[i].PID, theDeposits[i].Volume, theDeposits[i].Energy, theDeposits[i].Position[0], theDeposits[i].Position[1], theDeposits[i].Position[2], theDeposits[i].Time, theDeposits[i].Step);
      }
      fclose(out);
    }

    // Fill deposits variables
    for ( Int_t i=0; i<theEvent.NDeposits; ++i ) {
      dep_pdg[i]  = theDeposits[i].PID;
      dep_mat[i]  = theDeposits[i].Volume;
      //      std::cout << "event " << theEvent.EventID << " dep_mat[" << i << "] " << dep_mat[i] << std::endl;
      dep_ene[i]  = theDeposits[i].Energy;
      dep_step[i] = theDeposits[i].Step;
      dep_x[i]    = theDeposits[i].Position[0];
      dep_y[i]    = theDeposits[i].Position[1];
      dep_z[i]    = theDeposits[i].Position[2];
      dep_time[i] = theDeposits[i].Time;
      dep_r[i]    = sqrt(pow(dep_x[i],2)+pow(dep_y[i],2)+pow(dep_z[i],2));

      //--------------------------------------------------------------------------------
      // TPC
      //--------------------------------------------------------------------------------
      if ( dep_mat[i] == LArIndex ) {
        if ( firstTPC ) {  // this should not be necessary, because the deposits are time ordered
          firstTPC = false;
          timeTPC  = dep_time[i];
        }
        ++ndepTPC;
	eneTPC += dep_ene[i];
        if ( dep_pdg[i] < 100 )
          depElec += S1quench(dep_ene[i], 0.21);
        else
          depNucl += S1quench(dep_ene[i], 1) * 0.25;

        // Ar40
        Ar40ene     = 0.0;
        Ar40eneCalc = 0.0;
        Ar40x = 0.0;
        Ar40y = 0.0;
        Ar40z = 0.0;
        Ar40px = 0.0;
        Ar40py = 0.0;
        Ar40pz = 0.0;

        if ( dep_pdg[i] == 1000180400 ) {
          ++nAr40;
          Ar40ene = dep_ene[i];
          Ar40x   = dep_x[i];
          Ar40y   = dep_y[i];
          Ar40z   = dep_z[i];
          // momenta to be assigned
        }
      }

      //--------------------------------------------------------------------------------
      // Neutron Detectors
      //--------------------------------------------------------------------------------
      const Int_t num = dep_mat[i] - LSC0;
      //      std::cout << dep_mat[i] << std::endl;
      if ( 0 <= num && num < MAXLSC ) {
        if ( firstLSc[num] ) {
          firstLSc[num] = false;
          timeLSc[num]  = dep_time[i];
          // position of first hit
          x0LSc[num] = dep_x[i];
          y0LSc[num] = dep_y[i];
          z0LSc[num] = dep_z[i];
        }
        ++ndepLSc[num];
        eneLScRaw[num] += dep_ene[i];
        xeneLSc[num] += dep_ene[i] * dep_x[i];
        yeneLSc[num] += dep_ene[i] * dep_y[i];
        zeneLSc[num] += dep_ene[i] * dep_z[i];
        eneLSc[num] += QScintillator(dep_ene[i],dep_pdg[i]);
        if ( verbose )
          std::cout << theEvent.EventID << ' ' << num << ' ' << i << ' ' << dep_time[i] << ' ' << dep_pdg[i] << ' ' << dep_ene[i] << std::endl;
      }
    }

    //--------------------------------------------------------------------------------
    // LSci: do some analysis AFTER all deposits have been parsed
    //--------------------------------------------------------------------------------
    for ( Int_t num=0; num<MAXLSC; ++num ) {
      if ( timeLSc[num] > 0 ) {
        ++nTrgLSc;  // this LSci has energy deposit
        timeLSc[num] -= timeTPC;  // ToF from TPC
        if ( eneLScRaw[num] ) {
          xLSc[num] = xeneLSc[num] / eneLScRaw[num];
          yLSc[num] = yeneLSc[num] / eneLScRaw[num];
          zLSc[num] = zeneLSc[num] / eneLScRaw[num];
        }
        else
          std::cout << "To have timeLSc but 0 eneLSc should not happen" << std::endl;
        Double_t v = 1e9 * neutronDetectorRadius / timeLSc[num];
        v /= 3e8;
        eneLScTOF[num] = 0.5 * 1.00866554916 * 931.5 * pow(v,2);
        /*
        std::cout << "num " << num
                  << " timeLSc " << timeLSc[num]
                  << " timeTPC " << timeTPC
                  << " eneLScTOF " << eneLScTOF[num]
                  << std::endl;
        */
      }
    }

    if(depNucl + depElec > 0) f90TPC = depNucl / ( depNucl + depElec );

    // Fill user variables
    for(int i=0;i<theEvent.NUsers;++i) {
      INT1[i]    = theUsers[i].UserInt1;
      INT2[i]    = theUsers[i].UserInt2;
      FLOAT1[i]  = theUsers[i].UserFloat1;
      FLOAT2[i]  = theUsers[i].UserFloat2;
      DOUBLE[i]  = theUsers[i].UserDouble;
    }
    NPES1 = 0 ;
    NPES2 = 0 ;
    // Fill pe variables
    for(int i=0;i<theEvent.NPE;++i) {
      pe_pmt[i]  = thePhotoElectrons[i].PMT;
      pe_time[i] = thePhotoElectrons[i].Time;
      if (pe_time[i] < 10000) NPES1++ ;
      else NPES2++ ;
    }
    // Fill veto pe variables
    for(int i=0;i<theEvent.VetoNPE;++i) {
      veto_pe_pmt[i]  = theVetoPhotoElectrons[i].PMT;
      veto_pe_time[i] = theVetoPhotoElectrons[i].Time;
    }
     // Fill mu pe variables
    for(int i=0;i<theEvent.MuNPE;++i) {
      mu_pe_pmt[i]  = theMuPhotoElectrons[i].PMT;
      mu_pe_time[i] = theMuPhotoElectrons[i].Time;
    }
    // Fill photon variables
    for(int i=0;i<theEvent.NPH;++i) {
      ph_volume[i] = thePhotons[i].VolumeID;
      ph_pid[i]    = thePhotons[i].PID;
      ph_wl[i]     = thePhotons[i].Wavelength;
      ph_x[i]      = thePhotons[i].Position[0];
      ph_y[i]      = thePhotons[i].Position[1];
      ph_z[i]      = thePhotons[i].Position[2];
      ph_time[i]   = thePhotons[i].Time;
    }

    // fill the calculated recoil:
    // if there is a hit in the ONE LSci (UNFORTUNATELY I HAVE TO HARDCODE THE POSITION HERE)
    // if there is a 40Ar recoil in the TPC (in MC I can choose exactly one recoil, won't work for data)
    //    std::cout << "hurz" << std::endl;
    if ( ndepLSc[0] && nAr40 == 1 ) {
      Float_t ene0 = theEvent.Energy * MeV;  // in the fil it's stored in MeV!!!!

      // Gun position: currently without position uncertaincy
      TVector3 gun_pos(theEvent.Position[0]*mm, theEvent.Position[1]*mm, theEvent.Position[2]*mm);  // in the fil it's stored in mm!!!!  

      // Hit position in the TPC:
      //      TVector3 tpc_pos(Ar40x, Ar40y, Ar40z);  // cm, precise MC position
      TVector3 tpc_pos(rand->Gaus(Ar40x,sx), rand->Gaus(Ar40y,sy), rand->Gaus(Ar40z,sz));  // cm, using realistic position resolution

      // Hit position in the LSci:
      //      TVector3 LSci_pos(x0LSc[0], y0LSc[0], z0LSc[0]);  // cm, MC position, [0] works only if we have one LSci only, TODO
      TVector3 LSci_pos(198.39*cm, 25.33*cm, 0.00*cm);  // LSci center position (we don't know better)

      Ar40eneCalc = Erec(ene0, gun_pos, tpc_pos, LSci_pos) / keV;
      /*
      std::cout << "ene0 " << ene0/keV << " keV" << std::endl;
      std::cout << "gun_pos     "; gun_pos.Print();
      std::cout << "tpc_pos     "; tpc_pos.Print();
      std::cout << "LSci_pos    "; LSci_pos.Print();
      std::cout << "Ar40 energy " << Ar40eneCalc/keV << " keV (" << Ar40ene << ")" << std::endl;
      */
    }

    // Fill the event
    if ( TPCstore ) {
      if ( theEvent.S1Energy > 0 )
        dstree->Fill();
    }
    else {
      dstree->Fill();
    }
  }

  // Write the tree
  dstree->Write(); 

  // Close the root file
  ff->Close();

  // Close the binary file
  _bin_fstream->close();
  cout << "Rootfile " << rootfile.c_str() << " created! " << endl ;
  cout << "Bye...!" << endl ;
  return 0 ;
}
