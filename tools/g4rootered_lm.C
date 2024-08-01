// --------------------------------------------------------------------------//
/*
 * g4rooter specialised for DS20k S2only searches
 * paolo.agnes@cern.ch
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
#include "TRandom3.h"
#include "TMath.h"
#include "TCanvas.h"
#include "string.h"
#include "TMinuit.h"
#include "TVector3.h"
#include "algorithm"
#include "../include/DSEventStructure.hh"


using namespace std;

TRandom3 * ran = new TRandom3();

bool isLightTree = true ;
const int MAXDAUGHTER = 10000;
const int MAXDEPOSIT  = 100000;
const int MAXNPH      = 1;
const int MAXNPE      = 1;
const int MAXUSER     = 100;

//clustering
bool isToBeAdded;
// cluser params
float dist_max_xy = 500. ;    //  cm
float dist_max_z  = 0.2 ;  //  cm
float deltaT_max  = dist_max_z*10000.  ;  // ns

//SiPM paramters
double SiPM_time_jitter_stddev = 10 ;
double tile_time_jitter_stddev = 10 ;
double SiPM_QE                 = 0.4 ;
double SiPM_dark_cts_in_7us    = 1.;
double fXX_spe_resolution      = 0.42 ;

//Veto paramters
float VetoLY                 = 0.52 ; //pe/keV
float VetoOpticsFano         = 1 ;
double VetoAcquisitionWindow = 200e3 ; // ns

int ScintillatorIndex ;
int LArIndex;

struct Cluster {

  double       genPartPDG ;   // PDG of the particle generating the energy deposit
  double       genPartZ ;     // Z of the particle generating the energy deposit
  double       Length ;
  double       Radius ;
  double       Energy ;
  double       kinEne ;       // kinetic energy of the particle generating the energy deposit
  double       dEdx ;
  double       nucl;
  double       elec;
  double       T0 ;
  float        X0 ;
  float        Y0 ;
  float        Z0 ;
  float        X1 ;
  float        Y1 ;
  float        Z1 ;
  int          npe ;
  int         nDeposits;
  int         nElectrons ;
  int         nPhotons ;
  int         nExcitons ;

};


HeaderStructure                theHeader;
EventStructureDiskFormat       theEvent;
vector<DepositStructure>       theDeposits;
vector<DaughterStructure>      theDaughters;
vector<UserStructure>          theUsers ;
vector<PhotonStructure>        thePhotons;
vector<PhotoElectronStructure> thePhotoElectrons ;
vector<PhotoElectronStructure> theVetoPhotoElectrons ;
vector<PhotoElectronStructure> theMuPhotoElectrons ;
vector<Cluster>                Clusters ;


// tuned to reproduce data
float ZToTDrift (float Z)  {
  return -(Z-175.)/350.5 * 3.75e3;  //assuming 3.76 ms maximum drift time in DS20k
} // cm to us

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
  if(theEvent.NPE        > MAXNPE)      { cout << "Fatal: NPE > MAXNPE = "<<theEvent.NPE << endl ; exit(0) ;}
  if(theEvent.VetoNPE    > MAXNPE)      { cout << "Fatal: VetoNPE > MAXNPE" << endl ; exit(0) ;}
  if(theEvent.MuNPE      > MAXNPE)      { cout << "Fatal: MuNPE > MAXNPE" << endl ; exit(0) ;}


  for(int i=0; i<theEvent.NDaughters; i++) theDaughters.push_back(_readDaughter(file));
  for(int i=0; i<theEvent.NDeposits; i++)  theDeposits.push_back(_readDeposit(file));
  for(int i=0; i<theEvent.NUsers; i++)     theUsers.push_back(_readUser(file));
  for(int i=0; i<theEvent.NPH; i++)        thePhotons.push_back(_readPhoton(file));
  for(int i=0; i<theEvent.NPE; i++)        thePhotoElectrons.push_back(_readPhotoElectron(file));
  for(int i=0; i<theEvent.VetoNPE; i++)    theVetoPhotoElectrons.push_back(_readMuPhotoElectron(file));
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


double get_ds20k_s1_corr(double x,double y) {
/*   double  p0   =       -1.05133e-03   ; //3.58187e-03   1.01122e-02   1.75984e-01
   double  p1 =         -3.2242-05  ;// 2.42217e-03   8.08249e-03  -1.08269e-01
   double  p2 =         -6.07911e-8  ;// 1.10745e-04   1.16607e-04  -2.05519e+00
   double  p3 =          2.0442e-8  ;// 7.23814e-07   1.25878e-06  -2.58920e+01
   double  p4 =          9.71755e-9 ;//  7.14966e-07   1.00530e-06  -1.22457e+03
   double  p5  =        -2.1363e-10 ;//  6.13616e-09   1.08507e-08  -4.31370e+04
   double  p6   =        9.7782e-01 ;//  1.97079e-01   7.00541e-01   6.77370e-0
*/
  double  p0   =       -9.92391e-04   ; //3.58187e-03   1.01122e-02   1.75984e-01
  double  p1 =         -7.86036e-05  ;// 2.42217e-03   8.08249e-03  -1.08269e-01
  double  p2 =         -1.52499e-06  ;// 1.10745e-04   1.16607e-04  -2.05519e+00
  double  p3 =    1.23053e-08  ;// 7.23814e-07   1.25878e-06  -2.58920e+01
  double  p4 =    1.74786e-08 ;//  7.14966e-07   1.00530e-06  -1.22457e+03
  double  p5  =        -1.32184e-10 ;//  6.13616e-09   1.08507e-08  -4.31370e+04
  double  p6   =  9.90096e-01 ;//  1.97079e-01   7.00541e-01   6.77370e-0

  return p6+p0*x+p1*y+p2*x*y+p3*x*x*y+p4*y*y*x+p5*x*x*y*y ;
}


double interpolator(double ene, int dim, double *x, double *y) {
  for(int i=0;i<dim;++i) {
    if(x[i] > ene) {
      double m = (y[i] - y[i-1])/(x[i] - x[i-1]);
      double q = y[i] - m*x[i];
      return m*ene + q ;
    }
  }

  return 0;
}


//__________________________________________________________________________________
////////////////////////////////////////////////////////////////////////////////////
//                    main
////////////////////////////////////////////////////////////////////////////////////

int main (int argc, char *argv[]) {

  if(argc == 1 || (argc > 1 && !string(argv[1]).find("help")) )  {
  //if(file.find("help") < 10) {
    cout << "Usage: g4rooter_lm [FILE] [OPTIONS] [OUTPUT]" <<endl ;
    cout <<endl ;
    cout << " Options: " << endl ;
    cout << " nevents=N:  max number (N) of events to process (default: 10000)" <<endl ;
    cout << " skipN=N:    skip the first N events (default: 0)" <<endl ;
    cout << endl;
    cout << " Output: " << endl ;
    cout << " filename.root (default: FILE.root)" << endl ;
    cout << endl;
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

  int loop=2;
  while(argv[loop]) {
    string argument = argv[loop];
    if(!argument.find("nevents=")) {
      argument.erase(0,8);
      nevents = (int)  atoi(argument.c_str()) ;
      cout << "   max number of events to process: " << nevents << endl;
    }
    if(!argument.find("clz=")) {
      argument.erase(0,4);
      dist_max_z = (float)  atof(argument.c_str()) ;
      cout << "   clustering in z set to: " << dist_max_z << endl;
    }
    if(!argument.find("skipN=")) {
      argument.erase(0,6);
      skipNevents = atoi(argument.c_str()) ;
      cout << "   events to skip: " << skipNevents << endl;
    }
    if(!argument.find("ds20k_noise=")) {
      argument.erase(0,12);
      SiPM_dark_cts_in_7us = atof(argument.c_str()) ;
      cout << "   ds20k_noise: " << SiPM_dark_cts_in_7us << " cts/7us" << endl;
    }
    if(!argument.find("ds20k_SiPM_res=")) {
      argument.erase(0,15);
      SiPM_time_jitter_stddev = atof(argument.c_str()) ;
      cout << "   ds20k_SiPM_time_res: " << SiPM_time_jitter_stddev << " ns" << endl;
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
  int counter_ene = 0;

  // Read Header
  _readHeader(_bin_fstream);
  LArIndex = theHeader.LArIndex ;
  ScintillatorIndex = theHeader.ScintillatorIndex ;


  cout << "nBeamON: " << int ( theHeader.Events )  << endl ;
  cout << "DetectorFlag: " << theHeader.DetectorFlag << endl ;
  cout << "Rate: " << theHeader.Rate << endl ;

  if (theHeader.DetectorFlag >= 10 ) fXX_spe_resolution = 0.1 ;

  // event extra variables
  float radius = 0;

  // daughter variables
  int    Did[MAXDAUGHTER], Dpdg[MAXDAUGHTER], Dpid[MAXDAUGHTER], Dprocess[MAXDAUGHTER];
  float  Dene[MAXDAUGHTER], Dx[MAXDAUGHTER], Dy[MAXDAUGHTER], Dz[MAXDAUGHTER],
         Dr[MAXDAUGHTER], Dpx[MAXDAUGHTER], Dpy[MAXDAUGHTER], Dpz[MAXDAUGHTER];
  double Dtime[MAXDAUGHTER];

  // deposit variables
  int    * dep_pdg = new int[MAXDEPOSIT];
  int    * dep_mat = new int[MAXDEPOSIT];
  int    * dep_id = new int[MAXDEPOSIT];
  float  * dep_ene = new float[MAXDEPOSIT] ;
  float  * dep_qene = new float[MAXDEPOSIT] ;
  float  * dep_x = new float[MAXDEPOSIT] ;
  float  * dep_y = new float[MAXDEPOSIT] ;
  float  * dep_z = new float[MAXDEPOSIT] ;
  float  * dep_r= new float[MAXDEPOSIT] ;
  //float * dep_Prompt = new float[MAXDEPOSIT]  ;
  //float * dep_Del = new float[MAXDEPOSIT]  ;
  float  * dep_step = new float[MAXDEPOSIT];
  double * dep_time = new double[MAXDEPOSIT];

  // user variables
  int    INT1[MAXUSER], INT2[MAXUSER];
  float  FLOAT1[MAXUSER], FLOAT2[MAXUSER];
  double DOUBLE[MAXUSER];


  // pe variables
  int    *pe_pmt = new int[MAXNPE];
  int    pe_ch[38];
  double *pe_time = new double[MAXNPE];
  float  s1_max_frac ;
  // veto pe variables
  int    *veto_pe_pmt  = new int[MAXNPH];
  double *veto_pe_time = new double [MAXNPH];

  // mu pe variables
  int    *mu_pe_pmt  = new int[MAXNPH];
  double *mu_pe_time = new double[MAXNPH];


  // ph variables
  int    *ph_volume = new int[MAXNPH];
  int    *ph_pid    = new int[MAXNPH];
  float  *ph_wl     = new float[MAXNPH];
  float  *ph_x      = new float[MAXNPH] ;
  float  *ph_y      = new float[MAXNPH] ;
  float  *ph_z      = new float[MAXNPH];
  double *ph_time   = new double[MAXNPH];

  int ndepoTPC;
  int ds20npe ;
  float total_s1, total_s2;
  int ndeptpc, nclus, cl_ndep[MAXDEPOSIT], ab_mat, cl_npe[MAXDEPOSIT];
  float s1ene,s2ene, veto_visene, mu_visene, ene, qene, qnpe, tpcene, vetoene,
    muene , depTPCTot , depVeto, npeTPC, npeTPC400, depTPC400 ,
    depVeto70, eneTPC, eneVeto, eneslideVeto, timeslideVeto, npeslideVeto, depTot;
  double ab_time, timeVeto;
  float f90like , f90npe, f120like, f150like , f260like, f200like, cl_ene[MAXDEPOSIT], cl_true_ene[MAXDEPOSIT],cl_x[MAXDEPOSIT], cl_y[MAXDEPOSIT], cl_z[MAXDEPOSIT], cl_t[MAXDEPOSIT];
  int fall_tile ;
  float f90like_tile , f120like_tile, f150like_tile ,f260like_tile,   f200like_tile ;
  float f90like_noise , f120like_noise, f150like_noise ,f260like_noise,   f200like_noise ;
  float cl_nucl[MAXDEPOSIT], cl_elec[MAXDEPOSIT];

  double tpromptVeto, tdelayedVeto,tlateVeto;
  //int    epromptVeto, edelayedVeto,elateVeto;
  int     prompt_npeVeto    ;
  int     prompt_npeNoise    ;
  double  prompt_timeVeto   ;
  double  prompt_zVeto   ;
  double  prompt_rVeto   ;

  int     late_npeVeto       ;
  double  late_timeVeto    ;
  int     del_npeVeto       ;
  double  del_timeVeto       ;

  int s1npe, s2npe;

  float qdepVeto ;

  float tdrift , s1_corr ;

  float  mat_energy_fraction[100];

  int nphc;
  TTree *dstree  = new TTree("dstree","The G4DS Root Tree");
  TTree *header  = new TTree("header","The G4DS Root Tree Header");
  dstree->SetMaxVirtualSize(100000);

  header->Branch ("nBeamOn",       &theHeader.Events,          "nBeamOn/I" ) ;
  header->Branch ("nSaved",        &counter,                 "nSaved/I") ;
  header->Branch ("detector_flag", &theHeader.DetectorFlag,  "detector_flag/I") ;

  dstree->Branch("ev",             &theEvent.EventID,         "ev/I");
  dstree->Branch("pdg",            &theEvent.PDG,             "pdg/I");
  dstree->Branch("ene0",           &theEvent.Energy,          "ene0/F");
  if (!isLightTree )  {
    dstree->Branch("s1ene",          &theEvent.S1Energy,        "s1ene/F");
    dstree->Branch("s2ene",       &theEvent.S2Energy,  "s2ene/F");
    dstree->Branch("veto_visene",    &theEvent.VetoVisEnergy,  "veto_visene/F");
    dstree->Branch("mu_visene",      &theEvent.MuVisEnergy,  "mu_visene/F");
    dstree->Branch("vetoene",        &theEvent.VetoDepEnergy, "vetoene/F");
    dstree->Branch("muene",      &theEvent.MuDepEnergy,     "muene/F");
  }
  dstree->Branch("tpcene",          &theEvent.TPCDepEnergy,   "tpcene/F");
  dstree->Branch("x",       &theEvent.Position[0],     "x/F");
  dstree->Branch("y",       &theEvent.Position[1],     "y/F");
  dstree->Branch("z",       &theEvent.Position[2],     "z/F");
  if (!isLightTree ) {
    dstree->Branch("ene",            &ene,          "ene/F");
    dstree->Branch("r",              &radius,          "radius/F");
    dstree->Branch("px",             &theEvent.Direction[0],    "px/F");
    dstree->Branch("py",             &theEvent.Direction[1],    "py/F");
    dstree->Branch("pz",             &theEvent.Direction[2],    "pz/F");
    //dstree->Branch("bx",             &theEvent.CenterOfMass[0], "bx/F");
    //dstree->Branch("by",             &theEvent.CenterOfMass[1], "by/F");
    //dstree->Branch("bz",             &theEvent.CenterOfMass[2], "bz/F");
  }
  dstree->Branch("tdrift",         &tdrift ,                "tdrift/F");

  dstree->Branch("npe",            &theEvent.NPE ,        "npe/I");
  dstree->Branch("munpe" ,         &theEvent.MuNPE ,        "munpe/I");
  dstree->Branch("vnpe",           &theEvent.VetoNPE ,        "vnpe/I");
  dstree->Branch("nph",            &theEvent.NPH,        "nph/I");
  dstree->Branch("nusers",         &theEvent.NUsers,        "nusers/I");

  // dstree->Branch("s1npe",            &s1npe ,        "s1npe/I");
  // dstree->Branch("s2npe",            &s2npe ,        "s2npe/I");
  // dstree->Branch("s1",            &total_s1 ,        "s1/F");
  // dstree->Branch("s2",            &total_s2 ,        "s2/F");
  // dstree->Branch("s1_corr",        &s1_corr ,                "s1_corr/F");
  // if (theHeader.DetectorFlag == 10) {
  //   dstree->Branch("s1_tile",        &fall_tile ,                "s1_tile/I");
  //   //dstree->Branch("nphqe",            &ds20npe ,        "nphqe/I");
  // }

  dstree->Branch("ndeposits" ,     &theEvent.NDeposits,       "ndeposits/I");
  dstree->Branch("dep_id",      dep_id,                   "dep_id[ndeposits]/I");
  dstree->Branch("dep_pdg",     dep_pdg,                  "dep_pdg[ndeposits]/I");
  dstree->Branch("dep_mat",     dep_mat,                  "dep_mat[ndeposits]/I");
  dstree->Branch("dep_time",    dep_time,                 "dep_time[ndeposits]/D");
  dstree->Branch("dep_ene",     dep_ene,                  "dep_ene[ndeposits]/F")  ;
  dstree->Branch("dep_qene",     dep_qene,                  "dep_qene[ndeposits]/F")  ;
  dstree->Branch("dep_step",    dep_step,                 "dep_step[ndeposits]/F")  ;
  dstree->Branch("dep_x",       dep_x,                    "dep_x[ndeposits]/F");
  dstree->Branch("dep_y",       dep_y,                    "dep_y[ndeposits]/F");
  dstree->Branch("dep_z",       dep_z,                    "dep_z[ndeposits]/F") ;
  dstree->Branch("dep_r",       dep_r,                    "dep_r[ndeposits]/F") ;

  if (!isLightTree ) {
    dstree->Branch("ndaughters",     &theEvent.NDaughters,      "ndaughters/I");
    dstree->Branch("ndepositsTPC" ,  &ndepoTPC,                 "ndepositsTPC/I");
    dstree->Branch("nphc",         &nphc,        "nphc/I");

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


    dstree->Branch("mat_fraction",      &mat_energy_fraction[0],            "mat_fraction[100]/F") ;
  }

  // dstree->Branch("depVeto",   &depVeto,       "depVeto/F") ;
  // //should be equal to vetoene
  // dstree->Branch("qdepVeto",  &qdepVeto,       "qdepVeto/F") ;
  // dstree->Branch("prompt_npeVeto",  &prompt_npeVeto,      "prompt_npeVeto/I");
  // dstree->Branch("prompt_npeNoise",   &prompt_npeNoise,        "prompt_npeNoise/I");
  // dstree->Branch("prompt_timeVeto",  &prompt_timeVeto,      "prompt_timeVeto/D");
  // dstree->Branch("prompt_zVeto",       &prompt_zVeto,        "prompt_zVeto/D");
  // dstree->Branch("prompt_rVeto",       &prompt_rVeto,        "prompt_rVeto/D");
  // dstree->Branch("late_npeVeto",        &late_npeVeto,          "late_npeVeto/I");
  // dstree->Branch("late_timeVeto",       &late_timeVeto,          "late_timeVeto/D");
  dstree->Branch("userint1",    INT1,      "int1[nusers]/I");
  dstree->Branch("userint2",    INT2,      "int2[nusers]/I");
  dstree->Branch("userfloat1",  FLOAT1,      "float1[nusers]/F")  ;
  dstree->Branch("userfloat2",  FLOAT2,      "float2[nusers]/F");
  dstree->Branch("userdouble0", DOUBLE,      "double0[nusers]/D");
  // dstree->Branch("pe_time",     pe_time,      "pe_time[npe]/D");
  // dstree->Branch("pe_pmt",      pe_pmt,      "pe_pmt[npe]/I");

  dstree->Branch("nclus",         &nclus,            "nclus/I") ;
  dstree->Branch("cl_ene",        cl_ene,            "cl_ene[nclus]/F") ;
  dstree->Branch("cl_true_ene",   cl_true_ene,       "cl_true_ene[nclus]/F") ;
  dstree->Branch("cl_ndep",       cl_ndep,            "cl_ndep[nclus]/I") ;
  dstree->Branch("cl_x",          cl_x,                     "cl_x[nclus]/F") ;
  dstree->Branch("cl_y",          cl_y,                     "cl_y[nclus]/F") ;
  dstree->Branch("cl_z",          cl_z,                     "cl_z[nclus]/F") ;
  dstree->Branch("cl_t",          cl_t,                    "cl_t[nclus]/F") ;
  dstree->Branch("cl_npe",        cl_npe,                    "cl_npe[nclus]/I") ;
  dstree->Branch("cl_nucl",        cl_nucl,                    "cl_nucl[nclus]/F") ;
  dstree->Branch("cl_elec",        cl_elec,                    "cl_elec[nclus]/F") ;


  if (!isLightTree ) {
    //dstree->Branch("pe_time",     pe_time,                  "pe_time[npe]/D");
    //dstree->Branch("pe_pmt",      pe_pmt,                   "pe_pmt[npe]/I");

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

  }

  ab_time=0;
  ab_mat = -1 ;

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
    ene         = 0;
    for(int i=0;i<100;++i) mat_energy_fraction[i] = 0;

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

    // Fill deposits variables
    float qarAr   = 0.25 ;             // Ar in LAr
    float quench ;

    depTPCTot = 0;
    depVeto=0;
    depTot=0;
    npeTPC=0;
    npeTPC400=0;
    depVeto70=0;
    ndeptpc = 0 ;
    f90like = -1;
    f90npe = -1;
    f120like = -1;
    f150like = -1;
    f260like = -1;
    f200like = -1;
    eneTPC = 0;
    eneVeto= 0;
    timeVeto = 1e12;
    eneslideVeto = 0;
    timeslideVeto = 1e12;
    npeslideVeto = 0;
    tdrift = -1. ;
    s1_corr = 0. ;
    ds20npe   = 0 ;
    total_s1   = 0 ;
    total_s2   = 0 ;
    s1npe   = 0 ;
    s2npe   = 0 ;
    nphc = 0 ;

    // For the Clustering
    Cluster LocalClus;
    Clusters.clear();

    float  nuceneT = 0 , lepeneT = 0;
    float nucene, lepene, TotEne;
    float dep_dist_xy = 0,   dep_dist_z  = 0,  deltaT      = 0;
    isToBeAdded = false ;

    int neutron_id = 0 ;
    double tmp_maxpe = 0, tmp_time_maxpe = 0, t0 = 0;

    bool first = true;

    ndepoTPC = 0;
    float tot_dep_ene = 0;
    qdepVeto = 0 ;
    for(int i=0;i<theEvent.NDeposits;++i) {

      dep_id[i]   = i+1 ;
      dep_pdg[i]  = theDeposits[i].PID;
      dep_mat[i]  = theDeposits[i].Volume;
      dep_ene[i]  = theDeposits[i].Energy;
      dep_step[i] = theDeposits[i].Step;
      dep_x[i]    = theDeposits[i].Position[0];
      dep_y[i]    = theDeposits[i].Position[1];
      dep_z[i]    = theDeposits[i].Position[2];
      dep_time[i] = theDeposits[i].Time;
      dep_r[i]    = sqrt(pow(dep_x[i],2)+pow(dep_y[i],2));
      dep_qene[i] = 0;

      mat_energy_fraction[dep_mat[i]]+= dep_ene[i];
      tot_dep_ene += dep_ene[i];

      if(dep_mat[i] == 8) depTPCTot += dep_ene[i]; // TPC
      if(dep_mat[i] == 66) { // UAr veto
        depVeto += dep_ene[i];
        if(dep_pdg[i]< 30) qdepVeto += dep_ene[i];
        else if(dep_pdg[i] >= 30 && dep_pdg[i] < 3000)   qdepVeto += dep_ene[i] * 0.25;
      }
      if(dep_mat[i] == 8 ) ndepoTPC++;

      // CLUSTERING ALGORITHM
      if (dep_mat[i] != LArIndex ) continue ;
      // continue if DS50 geometry and deposits in LAr are outside the active volume
      if (theHeader.DetectorFlag <= 2  && (dep_x[i]*dep_x[i]+dep_y[i]*dep_y[i] > 17.77*17.77 || dep_z[i] > 14.7 ||  dep_z[i]<-21.439))  continue ;

      double alpha = 0;

      double qdepo       = dep_ene[i];
      dep_qene[i] = qdepo;

      isToBeAdded = true;
      for( int j = 0; j < int( Clusters.size() ); j++) {

        dep_dist_z =   abs( Clusters[j].Z0 - dep_z[i] ) ;
        deltaT     =   abs( Clusters[j].T0 - dep_time[i] );
        // clustering conditions
        if( dep_dist_z < dist_max_z &&  deltaT < deltaT_max ) {

        quench = 0 ;
        if (dep_pdg[i] < 30  )
          Clusters[j].elec += dep_qene[i] ;
        else
          Clusters[j].nucl += dep_qene[i] ;

         //weighted average
         Clusters[j].Z0 = ( Clusters[j].Z0*Clusters[j].Energy + dep_z[i]*dep_qene[i] ) / (dep_qene[i]+Clusters[j].Energy)  ;
         Clusters[j].X0 = ( Clusters[j].X0*Clusters[j].Energy + dep_x[i]*dep_qene[i] ) / (dep_qene[i]+Clusters[j].Energy)  ;
         Clusters[j].Y0 = ( Clusters[j].Y0*Clusters[j].Energy + dep_y[i]*dep_qene[i] ) / (dep_qene[i]+Clusters[j].Energy)  ;

         Clusters[j].Energy   += dep_qene[i] ;
         Clusters[j].nDeposits++;

         isToBeAdded = false;
         break;
       }
      }
      // Otherwise, new cluster
      if( isToBeAdded ){
        quench = 0 ;
        if (dep_pdg[i] < 30  ) {
          LocalClus.elec = dep_qene[i] ;
          LocalClus.nucl = 0;
        } else  {
          LocalClus.nucl = dep_qene[i] ;
          LocalClus.elec = 0;
        }

        LocalClus.genPartPDG =  dep_pdg[i] ;
        LocalClus.Energy     =  dep_qene[i];
        LocalClus.T0         =  dep_time[i];
        LocalClus.X0         =  dep_x[i] ;
        LocalClus.Y0         =  dep_y[i] ;
        LocalClus.Z0         =  dep_z[i] ;
        LocalClus.nDeposits  =  1;
        Clusters.push_back( LocalClus ) ;
      }
    }

    for(int i=0;i<100;++i) if(tot_dep_ene > 0) mat_energy_fraction[i] /= tot_dep_ene ;




    nclus = int(Clusters.size());
    //if(nclus != 1) continue ;
    //if(nclus == 0) continue ;

    double T0;
    if(nclus!=0)
      T0 = Clusters[0].T0 ;
    else
      T0 = 0;

    //write the clusters
    for (int i=0; i<Clusters.size() ; ++i) {
      cl_true_ene[i]   = Clusters[i].Energy;
      cl_x[i]          = Clusters[i].X0;
      cl_y[i]          = Clusters[i].Y0;
      cl_z[i]          = Clusters[i].Z0;
      cl_t[i]          = Clusters[i].T0;
      // QUENCHING is NOT accounted for
      if (Clusters[i].nucl >0.5 ) cl_nucl[i]     = Clusters[i].nucl;
      else cl_nucl[i]  = 0 ;
      //////////////////////////////////////////////
      if (Clusters[i].elec >0.5 ) cl_elec[i]   = Clusters[i].elec ;
      else cl_elec[i]  = 0 ;
      cl_ene[i]        = cl_elec[i]  +  cl_nucl[i] ;
      cl_ndep[i]       = Clusters[i].nDeposits;
    }

    tdrift  = ZToTDrift(cl_z[0]) ;

    // Fill user variables
    for(int i=0;i<theEvent.NUsers;++i) {
      INT1[i]    = theUsers[i].UserInt1;
      INT2[i]    = theUsers[i].UserInt2;
      FLOAT1[i]  = theUsers[i].UserFloat1;
      FLOAT2[i]  = theUsers[i].UserFloat2;
      DOUBLE[i]  = theUsers[i].UserDouble;
    }

    for(int i=0;i<theEvent.NPE;++i) {
      pe_pmt[i]  = thePhotoElectrons[i].PMT;
      pe_time[i] = thePhotoElectrons[i].Time;
    }

    if ( theHeader.DetectorFlag >= 10 )
      s1_corr = 1./get_ds20k_s1_corr(theEvent.Position[2],sqrt(pow(theEvent.Position[0], 2)+pow(theEvent.Position[1],2))) ;
    else  {
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

    double fX_t0 = 0. ;
    /////////////////////////////////////////////////////////////////////////////////
    // Fill photon variables
    // USED as for photoelectrons in DS20k
    /////////////////////////////////////////////////////////////////////////////////
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

  //TFile *ff = new TFile(rootfile.c_str(),"recreate");
  // Write the trees
  header->Fill() ;
  header->Write();

  dstree->Write();



  // Close the root file
  ff->Close();

  // Close the binary file
  _bin_fstream->close();
  cout << "Rootfile " << rootfile.c_str() << " created! " << endl ;
  cout << "Bye...!" << endl ;
  return 0 ;
}
