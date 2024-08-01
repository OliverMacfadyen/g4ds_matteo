// --------------------------------------------------------------------------//
/**
 * AUTHOR: D. Franco
 * CONTACT: dfranco@in2p3.fr
 *
 * Generate a root file reading the binary file from the g4ds output
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
#include "string.h"
#include "TMinuit.h"
#include "TVector3.h"
#include "algorithm"
#include "../include/DSEventStructure.hh"


using namespace std;

int MAXDAUGHTER = 70000;
int MAXDEPOSIT  = 2000000;
int MAXNPH      = 0;
int MAXNPE      = 2000000;
int MAXUSER     = 1000;
int MAXSIGNAL   = 1000;

struct ParisSignal{

  int    ndep ;
  int    quality;                    // 0: good, 1: fake, 2: contaminated
  float  eneTPC;                // True energy of the event
  float  s1ene;               // equiv to S1 ene
  float  s2ene;               // equiv to S1 ene
  double timeTPC;                // Time of the event in the TPC

  float  eneTOFDet1;                // Reconstructed energy of the nuclear recoil form TOF
  float  eneTOFDet2;                // Reconstructed energy of the nuclear recoil form TOF

  float  eneDet1;
  float  eneDet2;

  float  psdDet1;
  float  psdDet2;

  double timeDet1;
  double timeDet2;

  float f90;
};

ParisSignal mysignal ;

vector<float> vpdf, vpdf_time, vtheta, vlenght, vjitter_prob, vjitter_time;

HeaderStructure                theHeader;
EventStructureDiskFormat       theEvent;
vector<DepositStructure>       theDeposits;
vector<DaughterStructure>      theDaughters;
vector<UserStructure>          theUsers ;
vector<PhotonStructure>        thePhotons;
vector<PhotoElectronStructure> thePhotoElectrons ;
vector<PhotoElectronStructure> theVetoPhotoElectrons ;
vector<PhotoElectronStructure> theMuPhotoElectrons ;
vector<ParisSignal>            theSignals;

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
  file->read ((char *)(&event_size), sizeof (int));
  file->read ((char *)(&theHeader), sizeof( HeaderStructure)  );
  file->read ((char *)(&event_size2), sizeof (int));
  if (!(*file)) return false;
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
  theSignals.clear();

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


void reset_signal() {

  mysignal.ndep       = 0 ;
  mysignal.quality    = -1 ;
  mysignal.eneTPC     = 0;
  mysignal.s1ene      = 0;
  mysignal.s2ene      = 0;
  mysignal.timeTPC    = 0;
  mysignal.eneTOFDet1 = 0;
  mysignal.eneTOFDet2 = 0;
  mysignal.eneDet1    = 0;
  mysignal.eneDet2    = 0;
  mysignal.psdDet1    = -1;
  mysignal.psdDet2    = -1;
  mysignal.timeDet1   = 0;
  mysignal.timeDet2   = 0;
  mysignal.f90        = -1;
}



///quenching for TPC
double S1quench(double ene, double alpha) {

  double W = 19.5e-3; // keV
  double epsilon =1;


  double p0=    2.92032e-01 ;
  double p1=    1.37809e+00 ;
  double p2=   -1.55008e-01 ;
  double p3=   -1.46758e-02 ;
  double p4=    9.16802e-01 ;
  double p5=    6.21216e-01  ;
  double myRecoProb = p0 * (1 - p1*exp( p2 * ene )) *
            TMath::Exp( (p3) * pow( ene, (p4) ) ) +
            p5 ;

  double myNumQuanta   = ene /W ;
  double myNumIons     = myNumQuanta / ( 1 + alpha ) ;
  double myNumExcitons = myNumQuanta - myNumIons;

  double myphotons   = myNumExcitons + myNumIons*myRecoProb ;
  double out = myphotons*W ;

  return out;

}
double S2quench(double ene, double alpha) {

  double W = 19.5e-3; // keV
  double epsilon =1;


  double p0=    2.92032e-01 ;
  double p1=    1.37809e+00 ;
  double p2=   -1.55008e-01 ;
  double p3=   -1.46758e-02 ;
  double p4=    9.16802e-01 ;
  double p5=    6.21216e-01  ;
  double myRecoProb = p0 * (1 - p1*exp( p2 * ene )) *
            TMath::Exp( (p3) * pow( ene, (p4) ) ) +
            p5 ;

  double myNumQuanta   = ene /W ;
  double myNumIons     = myNumQuanta / ( 1 + alpha )*(1 - myRecoProb) ;
  double out           = myNumIons*W ;

  return out;

}


double QScintillator(double ene, int pdg) {
  if(pdg < 100) return ene ;
  else if(pdg <4000) return ene/3. ;
  else return ene/10 ;

  return 0;
}

//--------------------------------------------
//                    main
//--------------------------------------------
int main (int argc, char *argv[]) {

  if(argc == 1 || (argc > 1 && !string(argv[1]).find("help")) )  {
  //if(file.find("help") < 10) {
    cout << "Usage: g4rooter [FILE] [OPTIONS] [OUTPUT]" <<endl ;
    cout <<endl ;
    cout << " Options: " << endl ;
    cout << " nevents=N:  max number (N) of events to process (default: 10000)" <<endl ;
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
  int nevents     = 100000000;
  int skipNevents = 0;
  float kB        = 0.012;
  float LY        = 500.;

  float GateTime     = 1e4;
  string rootfile    = file;

  int loop=2;
  while(argv[loop]) {
    string argument = argv[loop];
    if(!argument.find("nevents=")) {
      argument.erase(0,8);
      nevents = (int)  atoi(argument.c_str()) ;
      cout << "   max number of events to process: " << nevents << endl;
    }
    if(!argument.find("gate=")) {
      argument.erase(0,5);
      GateTime = (float)  atof(argument.c_str()) ;
      cout << "   Acquisition Gate Time: " << GateTime << " ns" << endl;
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
    size_t found;
    found=argument.rfind(".root");
    if(found!=string::npos) {
      rootfile = argument;
      cout << "   output root file: " << argument << endl ;
    }

    loop++;
  }

  char extname[20];
  sprintf(extname,"_%dus.root",int(GateTime/1000.));


  //rootfile.replace(file.find(".fil"),4,".root");

  rootfile.replace(file.find(".fil"),4,extname);

  cout <<"Output Root File: "<< rootfile << endl ;


  TFile *ff = new TFile(rootfile.c_str(),"recreate");

  // event extra variables
  float radius = 0;

  // daughter variables
  int    Did[MAXDAUGHTER], Dpdg[MAXDAUGHTER], Dpid[MAXDAUGHTER], Dprocess[MAXDAUGHTER],
         Dmat[MAXDAUGHTER];
  float  Dene[MAXDAUGHTER], Dx[MAXDAUGHTER], Dy[MAXDAUGHTER], Dz[MAXDAUGHTER],
         Dr[MAXDAUGHTER], Dpx[MAXDAUGHTER], Dpy[MAXDAUGHTER], Dpz[MAXDAUGHTER];
  double Dtime[MAXDAUGHTER];

  int    *dep_pdg   = new int[MAXDEPOSIT];
  int    *dep_mat   = new int[MAXDEPOSIT];
  float  *dep_ene   = new float[MAXDEPOSIT];
  float  *dep_x     = new float[MAXDEPOSIT];
  float  *dep_y     = new float[MAXDEPOSIT];
  float  *dep_z     = new float[MAXDEPOSIT];
  float  *dep_r     = new float[MAXDEPOSIT];
  float  *dep_step  = new float[MAXDEPOSIT];
  double *dep_time  = new double[MAXDEPOSIT];
  int    *dep_track = new int[MAXDEPOSIT];

  // Paris Signals

  int    *s_id        = new int[MAXSIGNAL];
  int    *s_ndep      = new int[MAXSIGNAL];
  int    *s_quality   = new int[MAXSIGNAL];                    // 0: good, 1: fake, 2: contaminated
  float  *s_eneTPC    = new float[MAXSIGNAL];                // True energy of the event
  float  *s_s1ene     = new float[MAXSIGNAL];               // equiv to S1 ene
  float  *s_s2ene     = new float[MAXSIGNAL];               // equiv to S1 ene
  double *s_timeTPC   = new double[MAXSIGNAL];                // Time of the event in the TPC
  float  *s_eneTOFDet1= new float[MAXSIGNAL];                // Reconstructed energy of the nuclear recoil form TOF
  float  *s_eneTOFDet2= new float[MAXSIGNAL];                // Reconstructed energy of the nuclear recoil form TOF
  float  *s_eneDet1   = new float[MAXSIGNAL];
  float  *s_eneDet2   = new float[MAXSIGNAL];
  float  *s_psdDet1   = new float[MAXSIGNAL];
  float  *s_psdDet2   = new float[MAXSIGNAL];
  double *s_timeDet1  = new double[MAXSIGNAL];
  double *s_timeDet2  = new double[MAXSIGNAL];
  float  *s_f90       = new float[MAXSIGNAL];


  // user variables
  int    INT1[MAXUSER], INT2[MAXUSER];
  float  FLOAT1[MAXUSER], FLOAT2[MAXUSER];
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
  float  ph_wl[MAXNPH], ph_x[MAXNPH], ph_y[MAXNPH], ph_z[MAXNPH];
  double ph_time[MAXNPH];

  // licorne variables
  int ndepTPC;
  float f90TPC, eneDet1, eneDet2, recoEne1, recoEne2, psdDet1, psdDet2;
  double timeDet1,timeDet2, timeTPC;
  float recoEneNeuTOFDet1,recoEneNeuTOFDet2;

  int LArIndex;
  int nsignals ;

  float  thePulsePeriod;
  float  thePulseWidth;
  double theRunTime ;

  float s1ene,s2ene, veto_visene, mu_visene, ene, qene, qnpe, tpcene, vetoene, muene ;
  TTree *dstree = new TTree("dstree","The G4DS Root Tree");
  dstree->SetMaxVirtualSize(10000000);
  //dstree->SetMaxVirtualSize(10000000000);

  dstree->Branch("ev",             &theEvent.EventID,         "ev/I");
  dstree->Branch("pdg",            &theEvent.PDG,             "pdg/I");
  dstree->Branch("ene0",           &theEvent.Energy,          "ene0/F");
  dstree->Branch("s1ene",          &theEvent.S1Energy,        "s1ene/F");
  dstree->Branch("s2ene",          &theEvent.S2Energy,        "s2ene/F");
  dstree->Branch("veto_visene",    &theEvent.VetoVisEnergy,   "veto_visene/F");
  dstree->Branch("mu_visene",      &theEvent.MuVisEnergy,     "mu_visene/F");
  dstree->Branch("t0",             &theEvent.Time,            "t0/D");

  dstree->Branch("tpcene",          &theEvent.TPCDepEnergy,   "tpcene/F");
  dstree->Branch("vetoene",          &theEvent.VetoDepEnergy, "vetoene/F");
  dstree->Branch("muene",          &theEvent.MuDepEnergy,     "muene/F");
  dstree->Branch("ene",            &ene,		      "ene/F");
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

  dstree->Branch("vpe_time",    veto_pe_time,               "veto_pe_time[vnpe]/D");
  dstree->Branch("vpe_pmt",     veto_pe_pmt,                "veto_pe_pmt[vnpe]/I");

  dstree->Branch("mupe_time",   mu_pe_time,               "mu_pe_time[munpe]/D");
  dstree->Branch("mupe_pmt",    mu_pe_pmt,                "mu_pe_pmt[munpe]/I");

  dstree->Branch("ph_volume",    ph_volume,               "ph_volume[nph]/I");
  dstree->Branch("ph_pid",       ph_pid,                  "ph_pid[nph]/I");
  dstree->Branch("ph_wl",        ph_wl,                   "ph_wl[nph]/F");
  dstree->Branch("ph_x",         ph_x,                    "ph_x[nph]/F");
  dstree->Branch("ph_y",         ph_y,                    "ph_y[nph]/F");
  dstree->Branch("ph_z",         ph_z,                    "ph_z[nph]/F");
  dstree->Branch("ph_time",      ph_time,                 "ph_time[nph]/D");

  // PARIS
  dstree->Branch("nsignals",        &nsignals,        "nsignals/I") ;
  dstree->Branch("s_id",            s_id,             "s_id[nsignals]/I") ;
  dstree->Branch("s_ndep",          s_ndep,           "s_ndep[nsignals]/I") ;
  dstree->Branch("s_quality",       s_quality,        "s_quality[nsignals]/I") ;
  dstree->Branch("s_eneTPC",        s_eneTPC,         "s_eneTPC[nsignals]/F") ;
  dstree->Branch("s_s1ene",         s_s1ene,          "s_s1ene[nsignals]/F") ;
  dstree->Branch("s_s2ene",         s_s2ene,          "s_s2ene[nsignals]/F") ;
  dstree->Branch("s_timeTPC",       s_timeTPC,        "s_timeTPC[nsignals]/D") ;
  dstree->Branch("s_eneTOFDet1",    s_eneTOFDet1,     "s_eneTOFDet1[nsignals]/F") ;
  dstree->Branch("s_eneTOFDet2",    s_eneTOFDet2,     "s_eneTOFDet2[nsignals]/F") ;
  dstree->Branch("s_eneDet1",       s_eneDet1,        "s_eneDet1[nsignals]/F") ;
  dstree->Branch("s_eneDet2",       s_eneDet2,        "s_eneDet2[nsignals]/F") ;
  dstree->Branch("s_psdDet1",       s_psdDet1,        "s_psdDet1[nsignals]/F") ;
  dstree->Branch("s_psdDet2",       s_psdDet2,        "s_psdDet2[nsignals]/F") ;
  dstree->Branch("s_timeDet1",      s_timeDet1,       "s_timeDet1[nsignals]/D") ;
  dstree->Branch("s_timeDet2",      s_timeDet2,       "s_timeDet2[nsignals]/D") ;
  dstree->Branch("s_f90",           s_f90,            "s_f90[nsignals]/F") ;

  dstree->Branch("RunTime",         &theRunTime,      "RunTime/D") ;
  dstree->Branch("PulsePeriod",     &thePulsePeriod,  "PulsePeriod/F") ;
  dstree->Branch("PulseWidth",      &thePulseWidth,   "PulseWidth/F") ;
  dstree->Branch("GateTime",        &GateTime,        "GateTime/F") ;


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

  int counter = 0;

  // Read Header
  _readHeader(_bin_fstream);
  LArIndex = theHeader.LArIndex ;
  float theNeutronRate = theHeader.Rate ;

  float neutronDetectorRadius;

  if(theHeader.Rate < 0.0001){

    size_t found_cm = file.find("cm");

    string dist = file.substr(0,found_cm);
    dist.erase(0,dist.find("_")+1);

    if(found_cm == string::npos ){
      cout<<"warning : neutronDetectorRadius could not be determine"<<endl;
      neutronDetectorRadius = 0.1* 0.8 ;
    }else{
      neutronDetectorRadius = atof(dist.c_str()) * 0.1;
    }

  }else {
    neutronDetectorRadius = theHeader.Rate;
  }



  // Loop over the events
  for(int _i=0; _i<nevents; _i++) {

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

    thePulseWidth  =  theUsers[0].UserFloat1; // ns
    thePulsePeriod =  theUsers[0].UserFloat2; // ns
    theRunTime     =  theUsers[0].UserDouble; // ns


    // initialize variables
    ene         =  0;

    // licorne

    ndepTPC     =  0;
    eneDet1     =  0;
    eneDet2     =  0;
    recoEne1    =  0;
    recoEne2    =  0;
    f90TPC      = -1;
    psdDet1     = -1;
    psdDet2     = -1;
    timeDet1    = -1;
    timeDet2    = -1;
    timeTPC     = -1;
    recoEneNeuTOFDet1 = -1;
    recoEneNeuTOFDet2 = -1;

    float qarAr   = 0.25 ;             // quenching of Ar in LAr
    float depElec = 0 ;
    float depNucl = 0 ;


    // Fill daughter variables
    for(int i=0;i<theEvent.NDaughters;++i) {
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




    bool has_first_tpcHit  = false ;
    bool firstDet1 = false ;
    bool firstDet2 = false ;
    depNucl        = 0;
    depElec        = 0;

    reset_signal();

    double tpcTime0 = 0;


    // Fill deposits variables
    for(int i=0;i<theEvent.NDeposits;++i) {
      dep_pdg[i]  = theDeposits[i].PID;
      dep_mat[i]  = theDeposits[i].Volume;
      dep_ene[i]  = theDeposits[i].Energy;
      dep_step[i] = theDeposits[i].Step;
      dep_x[i]    = theDeposits[i].Position[0];
      dep_y[i]    = theDeposits[i].Position[1];
      dep_z[i]    = theDeposits[i].Position[2];
      dep_time[i] = theDeposits[i].Time;
      dep_r[i]    = sqrt(pow(dep_x[i],2)+pow(dep_y[i],2)+pow(dep_z[i],2));



      if(has_first_tpcHit) {
        if(dep_time[i] - tpcTime0 < GateTime) {
          if(dep_mat[i] == LArIndex) {
            mysignal.ndep++;
            if(dep_pdg[i] < 100) {
              mysignal.s1ene  += S1quench(dep_ene[i],0.21);
              mysignal.s2ene  += S2quench(dep_ene[i],0.21);
              depElec         += S1quench(dep_ene[i],0.21);
            } else {
              mysignal.s1ene  += S1quench(dep_ene[i],1)*0.25;
              mysignal.s2ene  += S2quench(dep_ene[i],1)*0.25;
              depNucl         += S1quench(dep_ene[i],1)*0.25;
            }
            mysignal.eneTPC += dep_ene[i];
            mysignal.quality = 1 ;
          } else if(dep_mat[i] == 46) { // neutron detector 1
            mysignal.eneDet1 += QScintillator(dep_ene[i],dep_pdg[i]);
            if(!firstDet1) {
              mysignal.timeDet1 = dep_time[i] - tpcTime0;
              firstDet1 = true ;
            }
          } else if(dep_mat[i] == 47) { // neutron detector 2
            mysignal.eneDet2 += QScintillator(dep_ene[i],dep_pdg[i]);
            if(!firstDet2) {
              mysignal.timeDet2 = dep_time[i] - tpcTime0;
              firstDet2 = true ;
            }

          }
        } else {
          has_first_tpcHit = false ;
          firstDet1 = false ;
          firstDet2 = false ;
          if(depElec+depNucl > 0)   mysignal.f90 = depNucl/(depElec+depNucl);
          if(mysignal.timeDet1 > 0 && mysignal.timeTPC > 0) {
            float v = 1e9 * neutronDetectorRadius/(mysignal.timeDet1 - mysignal.timeTPC );
            v /= 2.998e8;
            mysignal.eneTOFDet1 = 0.5 * 1.0086649160 * 931.494061 * pow(v,2);
          }
          if(mysignal.timeDet2 > 0 && mysignal.timeTPC > 0) {
            float v = 1e9 * neutronDetectorRadius/(mysignal.timeDet2 - mysignal.timeTPC );
            v /= 2.998e8;
            mysignal.eneTOFDet2 = 0.5 * 1.0086649160 * 931.494061 * pow(v,2);
          }

          theSignals.push_back(mysignal);
          depNucl  = 0;
          depElec  = 0;
          tpcTime0 = 0 ;
          reset_signal();
        }
      }

      // first hit in LAr
      if(!has_first_tpcHit && dep_mat[i] == LArIndex) {
        tpcTime0  = int(dep_time[i]/thePulsePeriod)*thePulsePeriod ;
        has_first_tpcHit  = true ;
        mysignal.timeTPC = dep_time[i] - tpcTime0;
        mysignal.ndep++;
        if(dep_pdg[i] < 100) {
          mysignal.s1ene  += S1quench(dep_ene[i],0.21);
          mysignal.s2ene  += S2quench(dep_ene[i],0.21);
          depElec         += S1quench(dep_ene[i],0.21);
        } else {
          mysignal.s1ene  += S1quench(dep_ene[i],1)*0.25;
          mysignal.s2ene  += S2quench(dep_ene[i],1)*0.25;
          depNucl         += S1quench(dep_ene[i],1)*0.25;
        }
        mysignal.eneTPC += dep_ene[i];
      }
    }

    nsignals = int(theSignals.size());

    for(int i=0;i<nsignals;++i) {
      s_id[i]       = i;
      s_ndep[i]     = theSignals[i].ndep;
      s_eneTPC[i]   = theSignals[i].eneTPC;
      s_s1ene[i]    = theSignals[i].s1ene;
      s_s2ene[i]    = theSignals[i].s2ene;
      s_timeTPC[i]  = theSignals[i].timeTPC;
      s_eneTOFDet1[i]   = theSignals[i].eneTOFDet1;
      s_eneTOFDet1[i]   = theSignals[i].eneTOFDet1;
      s_eneDet1[i]  = theSignals[i].eneDet1;
      s_eneDet2[i]  = theSignals[i].eneDet2;
      s_psdDet1[i]  = theSignals[i].psdDet1;
      s_psdDet2[i]  = theSignals[i].psdDet2;
      s_timeDet1[i] = theSignals[i].timeDet1;
      s_timeDet2[i] = theSignals[i].timeDet2;
      s_f90[i]      = theSignals[i].f90;


      s_quality[i] = 0 ;
      if(s_ndep[i] == 1 && (s_eneDet1[i] > 0 || s_eneDet2[i] > 0)) s_quality[i] = 1 ;
    }

/*

      // TPC
      if(dep_mat[i] == LArIndex)  {
        if(firstTPC) {
          firstTPC = false ;
          mysignal.timeTPC  = dep_time[i] ;
        }
        ndepTPC++ ;
        if(dep_pdg[i] < 100)  depElec += S1quench(dep_ene[i],0.21);
        else depNucl += S1quench(dep_ene[i],1)*0.25 ;
      }

      // Neutron Detector 1
      if(dep_mat[i] == 46)  {
        eneDet1 += QScintillator(dep_ene[i],dep_pdg[i]);
        if(firstDet1) {
          firstDet1 = false ;
          timeDet1  = dep_time[i] ;
        }
      }

      // Neutron Detector 2
      if(dep_mat[i] == 47)  {
        eneDet2 += QScintillator(dep_ene[i],dep_pdg[i]);
        if(firstDet2) {
          firstDet2 = false ;
          timeDet2  = dep_time[i] ;
        }
      }
    }

    if(timeDet1 > 0){
      double v = 1e9 * neutronDetectorRadius / (timeDet1 - timeTPC);
      v /= 3e8;
      recoEneNeuTOFDet1 = 0.5 * 1.00866554916 * 931.5 * pow(v,2);
    }

    if(timeDet2 > 0){
      double v = 1e9 * neutronDetectorRadius / (timeDet2 - timeTPC);
      v /= 3e8;
      recoEneNeuTOFDet2 = 0.5 * 1.00866554916 * 931.5 * pow(v,2);
    }



    if(depNucl + depElec > 0) f90TPC =  depNucl/(depNucl + depElec);

*/

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

    // Fill the event
    dstree->Fill();

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


/*
 * $Log: g4licorne_pulsed.C,v $
 * Revision 1.3  2015/12/02 15:02:35  riffard
 * Addition of the TOF energy reconstruction
 *
 * Revision 1.2  2015/11/29 13:23:11  dfranco
 * changed output name in licorne pulsed reco
 *
 * Revision 1.1  2015/11/28 12:55:51  dfranco
 * licorne generantor and analysis updated
 *
 * Revision 1.2  2015/11/27 09:04:49  riffard
 * addition of the TOF neutron energy reconstruction
 *
 * Revision 1.1  2015/11/26 13:30:31  dfranco
 * licorne
 *
 *
 */
