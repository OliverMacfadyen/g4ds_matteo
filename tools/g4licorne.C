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
int MAXDEPOSIT  = 10000;
int MAXNPH      = 100000;
int MAXNPE      = 2000000;
int MAXUSER     = 1000;


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
  double out = myphotons*W/ene ;

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

  int nevents     = 100000000;
  int skipNevents = 0;
  float kB        = 0.012;
  float LY        = 500.;
  bool  TPCstore  = 0;
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
  float radius = 0;

  // daughter variables
  int    Did[MAXDAUGHTER], Dpdg[MAXDAUGHTER], Dpid[MAXDAUGHTER], Dprocess[MAXDAUGHTER],
         Dmat[MAXDAUGHTER];
  float  Dene[MAXDAUGHTER], Dx[MAXDAUGHTER], Dy[MAXDAUGHTER], Dz[MAXDAUGHTER],
         Dr[MAXDAUGHTER], Dpx[MAXDAUGHTER], Dpy[MAXDAUGHTER], Dpz[MAXDAUGHTER];
  double Dtime[MAXDAUGHTER];

  // deposit variables
/*  int    dep_pdg[MAXDEPOSIT], dep_mat[MAXDEPOSIT], dep_track[MAXDEPOSIT];
  float  dep_ene[MAXDEPOSIT], dep_x[MAXDEPOSIT], dep_y[MAXDEPOSIT],
         dep_z[MAXDEPOSIT], dep_r[MAXDEPOSIT], dep_step[MAXDEPOSIT];
  double dep_time[MAXDEPOSIT];
*/
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
  float eneTPC, f90TPC, eneDet1, eneDet2, recoEne1, recoEne2, psdDet1, psdDet2;
  double timeDet1,timeDet2, timeTPC;
  float recoEneNeuTOFDet1,recoEneNeuTOFDet2;

  int LArIndex;

  float s1ene,s2ene, veto_visene, mu_visene, ene, qene, qnpe, tpcene, vetoene, muene ;
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

  // Licorne
  dstree->Branch("ndepTPC",         &ndepTPC,         "ndepTPC/I") ;
  dstree->Branch("eneTPC",          &eneTPC,          "eneTPC/F") ;       // f90 TPC
  dstree->Branch("f90TPC",          &f90TPC,          "f90TPC/F") ;       // f90 TPC
  dstree->Branch("eneDet1",         &eneDet1,         "eneDet1/F") ;      // energy of neutron detector 1
  dstree->Branch("eneDet2",         &eneDet2,         "eneDet2/F") ;      // energy of neutron detector 2
  dstree->Branch("timeDet1",        &timeDet1,        "timeDet1/D") ;     // time of neutron detector 1
  dstree->Branch("timeDet2",        &timeDet2,        "timeDet2/D") ;     // time of neutron detector 2
  dstree->Branch("timeTPC",         &timeTPC,         "timeTPC/D") ;       // TPC time
  dstree->Branch("recoEne1",        &recoEne1,        "recoEne1/F") ;     // reco energy of neutron detector 1
  dstree->Branch("recoEne2",        &recoEne2,        "recoEne2/F") ;     // reco energy of neutron detector 2
  dstree->Branch("psdDet1",         &psdDet1,         "psdDet1/F") ;     // PSD of neutron detector 1
  dstree->Branch("psdDet2",         &psdDet2,         "psdDet2/F") ;     // PSD of neutron detector 2
  dstree->Branch("recoEneNeuTOFDet1", &recoEneNeuTOFDet1, "recoEneNeuTOFDet1/F") ; // Neutron energy from TOF
  dstree->Branch("recoEneNeuTOFDet2", &recoEneNeuTOFDet2, "recoEneNeuTOFDet2/F") ; // Neutron energy from TOF



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

  const float neutronDetectorRadius = 0.1*theHeader.Rate;


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

    // initialize variables
    ene         =  0;

    // licorne

    ndepTPC     =  0;
    eneTPC      =  0;
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




    bool firstTPC  = true ;
    bool firstDet1 = true ;
    bool firstDet2 = true ;

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

      // TPC
      if(dep_mat[i] == LArIndex)  {
        if(firstTPC) {
          firstTPC = false ;
          timeTPC  = dep_time[i] ;
        }
        ndepTPC++ ;
	eneTPC += dep_ene[i] ;
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
    if(TPCstore) {
      if(theEvent.S1Energy > 0) dstree->Fill();
    } else { dstree->Fill();}

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
 * $Log: g4licorne.C,v $
 * Revision 1.4  2015/12/04 09:57:24  dfranco
 * added option to write only good events
 *
 * Revision 1.3  2015/12/03 17:56:31  dfranco
 * g4licorne update
 *
 * Revision 1.2  2015/11/27 09:04:49  riffard
 * addition of the TOF neutron energy reconstruction
 *
 * Revision 1.1  2015/11/26 13:30:31  dfranco
 * licorne
 *
 *
 */
