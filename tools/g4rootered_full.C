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

bool isLightTree = false ; 
const int MAXDAUGHTER = 100;
const int MAXDEPOSIT  = 20000;
const int MAXNPH      = 10000000;
const int MAXNPE      = 1000000;
const int MAXUSER     = 1000;


vector<float> vpdf, vpdf_time, vtheta, vlenght, vjitter_prob, vjitter_time;

static const int ntiles = 9800 ; // 14x14 motherboards for DS20k, 25 tiles per MB, top and bottom arrays
static const double MB_heigh = 5*5 + 2*0.25; //5 SiPMs + two times half the edge spacing
static const double Maximal_coordinate = 7*MB_heigh; //TPC has radius of 7 motherboards

static const int h_ntiles = 72 ; 
static const float tiles_lowedge = -180, tiles_maxedge = 180 ; 
TH2D * h_tile_top  = new TH2D("h_tile_top", "h_tile_top", h_ntiles,tiles_lowedge,tiles_maxedge, h_ntiles,tiles_lowedge,tiles_maxedge);
TH2D * h_tile_bot  = new TH2D("h_tile_bot", "h_tile_bot", h_ntiles,tiles_lowedge,tiles_maxedge, h_ntiles,tiles_lowedge,tiles_maxedge);

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
  int	       nDeposits;
  int	       nElectrons ;
  int	       nPhotons ;
  int	       nExcitons ;

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
  if ( theHeader.DetectorFlag < 4 )  return -(Z-14.13677)/35.570*376.   ;
  else return -(Z-120)/240 * 2.5e3;  //assuming 2.5 ms maximum drift time in DS20k  
} // cm to us

// official analysis
double s1_correction( double t_drift, double t_drift_max=376.) {
  Double_t z = t_drift/(0.5*t_drift_max); // note normalization is to 0.5*t_drift_max
  // looked at Kr peaks in 15us t_drift windows (Run5330+5340), and fit these to [0]*z^5 + [1]*z^4 + [2]*z^3+[3]*z^2+[4]*z+[5].
  Double_t fit_par0 = 0.0407;
  Double_t fit_par1 = -0.206;
  Double_t fit_par2 = 0.407;
  Double_t fit_par3 = -0.389;
  Double_t fit_par4 = 0.247;
  Double_t fit_par5 = 0.898;
  // normalizing all points on fitted curve to expected Kr peak at t_drift_max/2
  Double_t exp_Kr_peak_at_half_t_drift_max = fit_par0 + fit_par1 + fit_par2 + fit_par3 + fit_par4 + fit_par5;
  Double_t exp_Kr_peak_at_t_drift = fit_par0*z*z*z*z*z + fit_par1*z*z*z*z + fit_par2*z*z*z + fit_par3*z*z + fit_par4*z + fit_par5;
  return exp_Kr_peak_at_half_t_drift_max/exp_Kr_peak_at_t_drift; // s1 correction factor
}





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
  if(event_size != event_size2) return false ;
  return true;
} 

////////////////////////////////////////////////////////////////////////////////////
// Read noise from real veto baseline
////////////////////////////////////////////////////////////////////////////////////
bool readNoise(TH1F *hh){
  double npe[12],cum[12];
  double norma  = 0 ;
  double norma2 = 0 ;
  double tmp = 0 ;
 
  ifstream fin("noise_2.dat");
  for(int i=0;i<11;++i) {
    fin >> npe[i] >> tmp ;

    if(i ==0 )    cum[i] = tmp;
    else          cum[i] = cum[i-1] + tmp ;  
    
    npe[i] += 5 ;
    norma  += tmp ;
  }
  fin.close();

  for(int i=0;i<11;++i) {
    cum[i] /= norma ;
    //cout << npe[i] << " " << cum[i] << endl;
  }

  for(int i=0;i<1e6;++i) {
    double myran = ran->Uniform();
    for(int k=0;k<11;++k) 
      
      if(cum[k] > myran ) {
        double mm ;
        double qq ;
        if(k == 0) {
          mm = (npe[k] - 0)/(cum[k] - 0);
	  qq = npe[k] - mm*cum[k]; 
	} else {
          mm = (npe[k] - npe[k-1])/(cum[k] - cum[k-1]);
	  qq = npe[k] - mm*cum[k];
	}
        hh->Fill(myran*mm + qq);
        break ;
      }
  
  }

  return true ; 

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
  if(event_size != event_size2) return false ;
  return true;
}

////////////////////////////////////////////////////////////////////////////////////
// Quenching for TPC - we need to take into account the increase of Ly at low energy (due to higher recombination).
// With this function we calculate the amount of energy (in keV) that goes in S1 
// The recombination is supposed to be the same for both ER and NR, but the excitation/ionization ratio is different
////////////////////////////////////////////////////////////////////////////////////
double S1quench(double ene, double alpha) {

  double W = 19.5e-3; // keV
  double epsilon =1;

  // new model obtained with 37Ar and 83mKr constraint - oct '15 6th - as from DSLight3
  double p0=   2.96766e-01 ;
  double p1=   3.95496e+00 ;
  double p2=  -5.17812e-01 ;
  double p3=  -1.38485e-02 ;
  double p4=   9.12436e-01 ;
  double p5=   6.61046e-01 ;
  double myRecoProb = p0 * (1 - p1*exp( p2 * ene )) * 
            TMath::Exp( (p3) * pow( ene, (p4) ) ) + 
            p5 ;  

  double myNumQuanta   = ene /W ;
  double myNumIons     = myNumQuanta / ( 1 + alpha ) ;
  double myNumExcitons = myNumQuanta - myNumIons;
  
  double myphotons   = myNumExcitons + myNumIons*myRecoProb ;
  double out = myphotons*W/ene ;
  return out;//(epsilon*ene/W/(1+alpha) *(alpha + reco))*W/ene ;

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
  double  p3 =  	1.23053e-08  ;// 7.23814e-07   1.25878e-06  -2.58920e+01
  double  p4 =  	1.74786e-08 ;//  7.14966e-07   1.00530e-06  -1.22457e+03
  double  p5  =        -1.32184e-10 ;//  6.13616e-09   1.08507e-08  -4.31370e+04
  double  p6   =	9.90096e-01 ;//  1.97079e-01   7.00541e-01   6.77370e-0

  return p6+p0*x+p1*y+p2*x*y+p3*x*x*y+p4*y*y*x+p5*x*x*y*y ; 
}

////////////////////////////////////////////////////////////////////////////////////
// Quenching factors for VETO
////////////////////////////////////////////////////////////////////////////////////
double MeV = 1000.0;

// 18
double proton_energy[] = {0*MeV,0.029*MeV,0.094*MeV,0.2*MeV,0.34*MeV,
                          0.52*MeV, 0.72*MeV, 0.94*MeV,2*MeV,3*MeV,
                          4*MeV,6*MeV,10*MeV,20*MeV,30*MeV,
                          40*MeV,60*MeV,100*MeV};
// 18  
double proton_quenching[] = {0*MeV,0.003*MeV,0.005*MeV,0.009*MeV,0.020*MeV,
                             0.041*MeV,0.071*MeV,0.127*MeV,0.6*MeV,1*MeV,
			     1.6*MeV,3*MeV,6*MeV,13*MeV,20*MeV,30*MeV,45*MeV,70*MeV};
    
// 13
double alpha_energy[] = {0*MeV,.7*MeV,.85*MeV,1*MeV,1.2*MeV,
				       2*MeV,3*MeV,4*MeV,6*MeV,
				       8*MeV,10*MeV,20*MeV,30*MeV};
// 13
double alpha_quenching[] = {0*MeV,.017*MeV,.02*MeV,.025*MeV,.036*MeV,
				.06*MeV,.15*MeV,.22*MeV,.5*MeV,
				1*MeV,1.8*MeV,7*MeV,12*MeV};
    
//quenching for carbon nuclear recoils
//AW from Astropart. Phys. 16 (2002) 333-338 and NIM 33 (1965) 131-135
// 11
double carbon_energy[] = {0*MeV,0.046*MeV, 0.111*MeV, 0.229*MeV,
				       0.368*MeV, 0.500*MeV, 1.2*MeV, 2.0*MeV,
				       3.0*MeV, 4.0*MeV, 5.0*MeV};
// 11
double carbon_quenching[] = {0*MeV,0.0022*MeV, 0.0026*MeV, 0.0032*MeV,
				0.0044*MeV, 0.005*MeV, 0.007*MeV, 0.011*MeV,
				0.018*MeV, 0.025*MeV, 0.035*MeV};


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

double qalpha(double ene) {
   return interpolator(ene,13,alpha_energy, alpha_quenching)/ene * 0.548 ;
}
double qproton(double ene) {
   return interpolator(ene,18,proton_energy, proton_quenching)/ene  ;
} 
double qcarbon(double ene) {
   return interpolator(ene,11,carbon_energy, carbon_quenching)/ene  ;
}
double qelectron (double ene)/*keV*/ {
  //RNS
  //kB = 0.012 cm/MeV from Borexino
  //Birks Quenching parameterized as in
  //"The ionization quench factor in liquid-scintillation counting standardizations
  //Malonda, Carles
  //Applied Radiation and Isotopes 51 (1999) 183-188

  double A1 = 0.32903;
  double A2 = 0.14404;
  double A3 = 0.08059;
  double A4 = -0.04536E-3;
  double A5 = 0.15623;
  double A6 = 0.06611;
  return ( (A1 + A2*TMath::Log(ene) + A3*TMath::Log(ene)*TMath::Log(ene) + A4*TMath::Log(ene)*TMath::Log(ene)*TMath::Log(ene))/
	   (1 + A5*TMath::Log(ene) + A6*TMath::Log(ene)*TMath::Log(ene) + A4*TMath::Log(ene)*TMath::Log(ene)*TMath::Log(ene))  );
}


double getqfactor(double ene,int pdg) {
  //return 1;
  
  if(pdg - 1e9 == 20040) return qalpha(ene); 
  //if(pdg - 1e9 == 20040) return 1./50.; ///alpha from B has 1.47 MeV, quenching between 30 and 60 keV
  else if(pdg == 11 || pdg == -11 || pdg == 22 ) return qelectron(ene);
  else if(pdg == 2212 || pdg == 2112) return qproton(ene);
  else if(pdg - 1e9 > 30000 ) return qcarbon(ene);
  //else if(pdg - 1e9 > 30000 ) return 0;
  
  
  // else 
  //cout<<"quenching not computed, pdg = "<<" "<<pdg<<endl;

  return 1 ;

}



//__________________________________________________________________________________
////////////////////////////////////////////////////////////////////////////////////
//                    main 
////////////////////////////////////////////////////////////////////////////////////

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
  string rootfile = file;
  rootfile.replace(file.find(".fil"),4,".root");
  
  int nevents     = 100000000;
  int skipNevents = 0;
  float kB        = 0.012;
  float LY        = 500.;

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
    if(!argument.find("ds20k_tile_res=")) {  
      argument.erase(0,15);  
      tile_time_jitter_stddev = atof(argument.c_str()) ;
      cout << "   ds20k_SiPM_tile_res: " << tile_time_jitter_stddev << " ns" << endl;
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
  
   cout << "DetectorFlag: " << theHeader.DetectorFlag << endl ; 
   if (theHeader.DetectorFlag == 10 ) fXX_spe_resolution = 0.1 ; 
  
  

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
  
  //SiPM tile variables
  int tile_ch[ntiles];
  float tile_x[ntiles], tile_y[ntiles];
  int tile_s1npe[ntiles], tile_s2npe[ntiles], tile_s1max;
  float tile_bary_x, tile_bary_y, tile_s1maxfrac;
  int numTiles = ntiles;
  
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
  TTree *dstree = new TTree("dstree","The G4DS Root Tree");
  dstree->SetMaxVirtualSize(100000);
    
  dstree->Branch("ev",             &theEvent.EventID,         "ev/I");
  dstree->Branch("pdg",            &theEvent.PDG,             "pdg/I");
  dstree->Branch("ene0",           &theEvent.Energy,          "ene0/F");  
  if (!isLightTree )  { 
    dstree->Branch("s1ene",          &theEvent.S1Energy,        "s1ene/F");     
    dstree->Branch("s2ene",	     &theEvent.S2Energy,	"s2ene/F");	
    dstree->Branch("veto_visene",    &theEvent.VetoVisEnergy,	"veto_visene/F");    
    dstree->Branch("mu_visene",      &theEvent.MuVisEnergy,	"mu_visene/F");      
    dstree->Branch("vetoene",	      &theEvent.VetoDepEnergy, "vetoene/F");	 
    dstree->Branch("muene",	    &theEvent.MuDepEnergy,     "muene/F");     
  }
  dstree->Branch("tpcene",          &theEvent.TPCDepEnergy,   "tpcene/F");     
  dstree->Branch("x",		   &theEvent.Position[0],     "x/F");		     
  dstree->Branch("y",		   &theEvent.Position[1],     "y/F");		     
  dstree->Branch("z",		   &theEvent.Position[2],     "z/F");		     
  if (!isLightTree ) { 
    dstree->Branch("ene",            &ene,		      "ene/F");        
    dstree->Branch("r",              &radius,		      "radius/F");     
    dstree->Branch("px",             &theEvent.Direction[0],    "px/F");	       
    dstree->Branch("py",             &theEvent.Direction[1],    "py/F");	       
    dstree->Branch("pz",             &theEvent.Direction[2],    "pz/F");	       
    //dstree->Branch("bx",             &theEvent.CenterOfMass[0], "bx/F");	       
    //dstree->Branch("by",             &theEvent.CenterOfMass[1], "by/F");	       
    //dstree->Branch("bz",             &theEvent.CenterOfMass[2], "bz/F");	       
  }
  dstree->Branch("tdrift",         &tdrift ,	              "tdrift/F");	   
     
  dstree->Branch("npe",            &theEvent.NPE ,	      "npe/I");        
  dstree->Branch("munpe" ,         &theEvent.MuNPE ,	      "munpe/I");      
  dstree->Branch("vnpe",           &theEvent.VetoNPE ,        "vnpe/I");       
  dstree->Branch("nph",            &theEvent.NPH,	      "nph/I");        

  dstree->Branch("s1npe",            &s1npe ,	      "s1npe/I");        
  dstree->Branch("s2npe",            &s2npe ,	      "s2npe/I");        
  dstree->Branch("s1",            &total_s1 ,	      "s1/F");        
  dstree->Branch("s2",            &total_s2 ,	      "s2/F");        
  dstree->Branch("s1_corr",        &s1_corr ,	              "s1_corr/F");	
  if (theHeader.DetectorFlag == 10) {
    dstree->Branch("s1_tile",        &fall_tile ,	              "s1_tile/I");	  
    //dstree->Branch("nphqe",            &ds20npe ,	      "nphqe/I");        
  }
    
  if (!isLightTree ) { 
    dstree->Branch("ndaughters",     &theEvent.NDaughters,      "ndaughters/I"); 
    dstree->Branch("ndeposits" ,     &theEvent.NDeposits,       "ndeposits/I");  
    dstree->Branch("ndepositsTPC" ,  &ndepoTPC,                 "ndepositsTPC/I");
    dstree->Branch("nusers",         &theEvent.NUsers,	      "nusers/I");    
    dstree->Branch("nphc",         &nphc,	      "nphc/I");    

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

    dstree->Branch("mat_fraction",      &mat_energy_fraction[0],	          "mat_fraction[100]/F") ;
  }
  
  dstree->Branch("depVeto",	 &depVeto,	     "depVeto/F") ;
  //should be equal to vetoene
  dstree->Branch("qdepVeto",	&qdepVeto,	     "qdepVeto/F") ;
  dstree->Branch("prompt_npeVeto",	&prompt_npeVeto,		  "prompt_npeVeto/I"); 
  dstree->Branch("prompt_npeNoise",	 &prompt_npeNoise,		    "prompt_npeNoise/I"); 
  dstree->Branch("prompt_timeVeto",	&prompt_timeVeto,		  "prompt_timeVeto/D");      
  dstree->Branch("prompt_zVeto",       &prompt_zVeto,		    "prompt_zVeto/D");      
  dstree->Branch("prompt_rVeto",       &prompt_rVeto,		    "prompt_rVeto/D");      
  dstree->Branch("late_npeVeto",        &late_npeVeto,		      "late_npeVeto/I"); 
  dstree->Branch("late_timeVeto",       &late_timeVeto,		      "late_timeVeto/D");      
  dstree->Branch("userint1",    INT1,			"int1[nusers]/I");  
  dstree->Branch("userint2",    INT2,			"int2[nusers]/I");
  dstree->Branch("userfloat1",  FLOAT1,			"float1[nusers]/F")  ; 
  dstree->Branch("userfloat2",  FLOAT2,			"float2[nusers]/F");      	   
  dstree->Branch("userdouble0", DOUBLE,			"double0[nusers]/D");  
  dstree->Branch("pe_time",     pe_time,  		"pe_time[npe]/D");   	   
  dstree->Branch("pe_pmt",      pe_pmt,			"pe_pmt[npe]/I");   

  if (theHeader.DetectorFlag == 10) {
    
    dstree->Branch("f180like",      &f90like,	          "f180/F") ;
    dstree->Branch("f200like",      &f200like,	          "f200/F") ;
    dstree->Branch("f220like",      &f120like,	          "f220/F") ;
    //dstree->Branch("f240like",      &f150like,	          "f240/F") ;
    //dstree->Branch("f260like",      &f260like,	          "f260/F") ;
    dstree->Branch("f180like_noise",      &f90like_noise,	                  "f180_noise/F") ;
    dstree->Branch("f200like_noise",      &f200like_noise,	          "f200_noise/F") ;
    dstree->Branch("f220like_noise",      &f120like_noise,	          "f220_noise/F") ;
    //dstree->Branch("f240like_noise",      &f150like_noise,	          "f240_noise/F") ;
    //dstree->Branch("f260like_noise",      &f260like_noise,	          "f260_noise/F") ;
    
    dstree->Branch("f180like_tile",      &f90like_tile,	                  "f180_tile/F") ;
    dstree->Branch("f200like_tile",      &f200like_tile,	          "f200_tile/F") ;
    dstree->Branch("f220like_tile",      &f120like_tile,	          "f220_tile/F") ;
    //dstree->Branch("f240like_tile",      &f150like_tile,	          "f240_tile/F") ;
    //dstree->Branch("f260like_tile",      &f260like_tile,	          "f260_tile/F") ;
    
  } else {
    dstree->Branch("f90",         &f90like,	          "f90/F") ;
    if (!isLightTree )  dstree->Branch("f90npe",      &f90npe,	          "f90npe/F") ;
    
  }
  //dstree->Branch("ab_time",    &ab_time,                   "ab_time/D");   
  //dstree->Branch("ab_mat",    &ab_mat,                     "ab_mat/I");   
  dstree->Branch("nclus",         &nclus,	          "nclus/I") ;
  if (!isLightTree ) {
    dstree->Branch("cl_ene",        cl_ene,	          "cl_ene[nclus]/F") ;
    dstree->Branch("cl_true_ene",   cl_true_ene,   	  "cl_true_ene[nclus]/F") ;
    dstree->Branch("cl_ndep",       cl_ndep,	          "cl_ndep[nclus]/I") ;
    dstree->Branch("cl_x",          cl_x,	                   "cl_x[nclus]/F") ;
    dstree->Branch("cl_y",          cl_y,	                   "cl_y[nclus]/F") ;
    dstree->Branch("cl_z",          cl_z,	                   "cl_z[nclus]/F") ;
    dstree->Branch("cl_t",          cl_t,	                  "cl_t[nclus]/F") ;
    dstree->Branch("cl_npe",        cl_npe,	                  "cl_npe[nclus]/I") ;
    dstree->Branch("cl_nucl",        cl_nucl,	                  "cl_nucl[nclus]/F") ;
    dstree->Branch("cl_elec",        cl_elec,	                  "cl_elec[nclus]/F") ;
  
 
    //dstree->Branch("pe_time",     pe_time,                  "pe_time[npe]/D");
    //dstree->Branch("pe_pmt",      pe_pmt,                   "pe_pmt[npe]/I");   
    if (theHeader.DetectorFlag != 10)  dstree->Branch("pe_ch",       pe_ch,                    "pe_ch[38]/I");   
    if (theHeader.DetectorFlag != 10)  dstree->Branch("s1_max_frac", &s1_max_frac,             "s1_max_frac/F");   
    
    dstree->Branch("vpe_time",  veto_pe_time,               "veto_pe_time[vnpe]/D");   
    dstree->Branch("vpe_pmt",   veto_pe_pmt,                "veto_pe_pmt[vnpe]/I"); 
    
    dstree->Branch("mupe_time",   mu_pe_time,               "mu_pe_time[munpe]/D");   
    dstree->Branch("mupe_pmt",    mu_pe_pmt,                "mu_pe_pmt[munpe]/I");   
    
    //dstree->Branch("ph_volume",    ph_volume,               "ph_volume[nph]/I");
    //dstree->Branch("ph_pid",       ph_pid,                  "ph_pid[nph]/I");
    dstree->Branch("ph_wl",        ph_wl,                   "ph_wl[nph]/F");
    dstree->Branch("ph_x",         ph_x,                    "ph_x[nph]/F");
    dstree->Branch("ph_y",         ph_y,                    "ph_y[nph]/F");
    dstree->Branch("ph_z",         ph_z,                    "ph_z[nph]/F");
    dstree->Branch("ph_time",      ph_time,                 "ph_time[nph]/D");
    
    ///channel number for a given tile (matches algorithm in DSTrackingAction
    ///for a given PE, its tile number can be obtained with the pe_pmt variable
    ///index of tile array is same as channel number
    ///  bottom array channels: [0,4899]      top array channels: [4900,9799]
    //dstree->Branch("ntiles",         &numTiles,       "ntiles/I");
    //dstree->Branch("tile_ch",        tile_ch,         "tile_ch[ntiles]/I");     //tile channel number, matches index of tile_arrays
    //dstree->Branch("tile_x",         tile_x,          "tile_x[ntiles]/F");      // center position of tile in cm
    //dstree->Branch("tile_y",         tile_y,          "tile_y[ntiles]/F");      
    //dstree->Branch("tile_s1npe",     tile_s1npe,      "tile_s1npe[ntiles]/I");  //number of pe in S1 for a given tile
    dstree->Branch("tile_s1max",     &tile_s1max,     "tile_s1max/I");            //index of tile with maximum number of S1 PEs
    dstree->Branch("tile_s1maxfrac", &tile_s1maxfrac, "tile_s1maxfrac/F");        //fraction of PE in tile with maximum # of S1 PEs
    //dstree->Branch("tile_s2npe",     tile_s2npe,      "tile_s2npe[ntiles]/I");  //number of pe in S2 for a given tile
    dstree->Branch("tile_bary_x",    &tile_bary_x,    "tile_bary_x/F");           //charge barycenter calculated using pe in a given channel in cm
    dstree->Branch("tile_bary_y",    &tile_bary_y,    "tile_bary_y/F");
  }
  
  for(int my=0; my<14; my++) {  //loop over motherboard indicies for 14x14 MB grid
    for(int mx=0; mx<14; mx++) { 
      for(int tx=0; tx<5; tx++) {  //loop over tile indicies for 5x5 tiles per MB
        for(int ty=0; ty<5; ty++) {
          int chBot = mx*25 + my*25*14 + tx*5 + ty; 
          int chTop = mx*25 + my*25*14 + tx*5 + ty + 196*25; //adjust channel number for top SiPM array
          tile_ch[chBot] = chBot;
          tile_ch[chTop] = chTop;
          
          //  pos = bottom left of MB + dist to center of tile + shift to center of tile + left edge of MB - shift to centered on (0,0)
          tile_x[chBot] = mx*MB_heigh + tx*5 + 2.5 + 0.25 - Maximal_coordinate;
          tile_x[chTop] = mx*MB_heigh + tx*5 + 2.5 + 0.25 - Maximal_coordinate;
          tile_y[chBot] = my*MB_heigh + ty*5 + 2.5 + 0.25 - Maximal_coordinate;
          tile_y[chTop] = my*MB_heigh + ty*5 + 2.5 + 0.25 - Maximal_coordinate;
          
          if (theHeader.DetectorFlag==12)  { //if ProtoProto, adjust coord so tile is centered on (0,0)
            tile_x[chBot] -= MB_heigh/2.0 ; 
            tile_x[chTop] -= MB_heigh/2.0 ; 
            tile_y[chBot] -= MB_heigh/2.0 ; 
            tile_y[chTop] -= MB_heigh/2.0 ; 
          }
        }
      }
    }
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
    
    //if(s1ene < 0.00000001) continue ;
    
    
    // sort photon and photoelectron vectors by time
    std::sort(theDeposits.begin(),theDeposits.end(), cmp_deposit());
    std::sort(thePhotons.begin(),thePhotons.end(), cmp_photon());
    std::sort(thePhotoElectrons.begin(),thePhotoElectrons.end(), cmp_photoelectron());
    std::sort(theVetoPhotoElectrons.begin(),theVetoPhotoElectrons.end(), cmp_photoelectron());
    std::sort(theMuPhotoElectrons.begin(),theMuPhotoElectrons.end(), cmp_photoelectron());
    
    // initialize variables 
    ene         = 0;
    //TPC
    double tp0 = 2.65965e-01;
    double tp1 = 4.92038e-01;
    double tp2 = 1.58933e-04;   
    double tp3 = -3.46480e-02;
    double tp4 = 7.67703e-01;
    double tp5 = 6.31385e-01;
    double tp6 = -5.64532e-02;
    double tp7 = -7.86833e-06;
    double fDriftField = 200;
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
    float TPCLY   = 7.100;   //photons per keV
    float qvalpha = 0.027;             // alpha in veto (40keVee for 1.47MeV alphas)
    float qvli    = 0.0 ;               // Li in veto  -  TODO         
    float qvp     = 0.3333;            // p in veto 
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

      if(dep_mat[i] == 8) depTPCTot += dep_ene[i];
      if(dep_mat[i] == 17 || dep_mat[i] == 2 || dep_mat[i] == 52) {
        depVeto += dep_ene[i];
	if(dep_pdg[i]< 30) qdepVeto += dep_ene[i];
	else if(dep_pdg[i] >= 30 && dep_pdg[i] < 3000)   qdepVeto += dep_ene[i]/3.5;
	else qdepVeto += dep_ene[i]/25.;
      }
      if(dep_mat[i] == 8 ) ndepoTPC++;
      
      //Some clusters are produced (especially by n captures) after the end of the aquisition window. We should account for this time difference
      //if(dep_time[i] > 400e3) break  ; 
      
      
      // CLUSTERING ALGORITHM
      if (dep_mat[i] != LArIndex ) continue ; 
      // continue if DS50 geometry and deposits in LAr are outside the active volume
      if (theHeader.DetectorFlag <= 2  && (dep_x[i]*dep_x[i]+dep_y[i]*dep_y[i] > 17.77*17.77 || dep_z[i] > 14.7 ||  dep_z[i]<-21.439))  continue ;
      
      

      double alpha = 0;
/*   // This should be done later, once the total cluster energy is defined
      // TPC quenching
      if(dep_pdg[i] - 1e9 > 100000)
	alpha = 1;
      else
	alpha = 0.21;

      double s1_dep = S1quench(dep_ene[i],alpha);
      double qdepo  = s1_dep * dep_ene[i];
  
      if(dep_pdg[i] -1e9 >  100000) 
	qdepo *= qarAr;     // Lindhard factor
*/	
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
    
    
    
    
    //select a specific number of clusters
    nclus = int(Clusters.size());
    //if(nclus != 1) continue ;
    //if(nclus == 0) continue ;
  
    double T0;
    if(nclus!=0)
      T0 = Clusters[0].T0 ;
    else
      T0 = 0;
      

    bool   prompt_firstVeto  = true ;

    double tmp_maxpe_late    = 0; 
    double t0_late           = -100;

    double tmp_maxpe_delayed = 0 ;
    double t0_delayed        = -100;
    
    double elateVeto        = 0;
    double edelayedVeto      = 0;
    double enelateVeto       = 0 ;
    double enedelayedVeto    = 0 ;
    double delnpeVeto        = 0 ;
    double latenpeVeto        = 0 ;
  
    prompt_npeVeto   = 0;
    prompt_npeNoise  = 0;
    prompt_timeVeto  = -10000;
    prompt_zVeto     = -10000;
    prompt_rVeto     = -10000;
    late_npeVeto     = 0;
    late_timeVeto    = -10000;
    del_npeVeto      = 0;
    del_timeVeto     = -10000;



    for(int k=0;k<theEvent.NDeposits;++k) {
    
      if(dep_mat[k] != ScintillatorIndex) continue ;
      //same as before! some events are otside the acquisition window. How to treat them? 
      if(dep_time[k] > VetoAcquisitionWindow) break ;
      
      double deltaT = dep_time[k] - T0 ;
      double quench = getqfactor(dep_ene[k],dep_pdg[k]);
      // prompt cut 
      if(deltaT > -10 && deltaT  < 200) {
        if(prompt_firstVeto) {
          prompt_timeVeto = dep_time[k] - T0;
          prompt_firstVeto = false ;
	  prompt_zVeto = dep_z[k];
	  prompt_rVeto = dep_r[k];
        }
        int vnpe = 0 ; 
	double vene = quench*VetoLY*dep_ene[k];
        if ( vene < 20 ) vnpe = ran->Poisson(vene)  ;
        else vnpe =  int(ran->Gaus( vene, VetoOpticsFano * sqrt(vene)) + 0.5) ;
	prompt_npeVeto += vnpe ;
	
      }
      //if(deltaT > 8800) {
        //tmp_maxpe_late	   = quench*dep_ene[k];
        //t0_late  	   = dep_time[k];
      //}
      // old delayed
      //if(deltaT < 200000 && deltaT > 200) {
      //  tmp_maxpe_delayed  = quench*dep_ene[k];
       // t0_delayed         = dep_time[k];
      //the old delayed becomes the late
      if(deltaT < VetoAcquisitionWindow && deltaT > 200) {
        tmp_maxpe_late  = quench*dep_ene[k];
        t0_late         = dep_time[k];
      }

      for(int kk=k+1;kk<theEvent.NDeposits;++kk){
	
	if(dep_mat[kk] != ScintillatorIndex) continue ;
	
	// late
	//if(dep_time[kk] - T0 > 8800)  {
	//if( dep_time[kk] - t0_late < 300.){
	//tmp_maxpe_late += getqfactor(dep_ene[kk],dep_pdg[kk])*dep_ene[kk];
	//} else break;
	//}
	
	// delayed
	//if(dep_time[kk] - T0 < 200000 && dep_time[kk] - T0 > 200)  {
	//  if( dep_time[kk] - t0_delayed < 300.){
	//    tmp_maxpe_delayed += getqfactor(dep_ene[kk],dep_pdg[kk])*dep_ene[kk];
	//  } else break;
	
	//the old delayed becomes the late 
	if(dep_time[kk] - T0 <VetoAcquisitionWindow&& dep_time[kk] - T0 > 200)  {
	  if( dep_time[kk] - t0_late < 300.){
	    tmp_maxpe_late += getqfactor(dep_ene[kk],dep_pdg[kk])*dep_ene[kk];
	  } else break;
	
	}
	
      }
      
      //sliding window simulation
      if(tmp_maxpe_late > elateVeto){
	elateVeto  = tmp_maxpe_late;
	late_timeVeto  = t0_late - T0;
      }  

      //if(tmp_maxpe_delayed > edelayedVeto){
      //edelayedVeto  = tmp_maxpe_delayed;
      //del_timeVeto  = t0_delayed - T0;
      //}
    }
    
    //if(edelayedVeto == 0)              del_npeVeto = 0 ;
    //else if(edelayedVeto*VetoLY < 20 ) del_npeVeto = ran->Poisson(edelayedVeto*VetoLY);
    //else                               del_npeVeto = int(ran->Gaus(edelayedVeto*VetoLY,sqrt(edelayedVeto*VetoLY))+0.5 );

    if(elateVeto == 0)                 late_npeVeto = 0;
    else if(elateVeto*VetoLY < 20 )    late_npeVeto = ran->Poisson(elateVeto*VetoLY);
    else                               late_npeVeto = int(ran->Gaus(elateVeto*VetoLY,VetoOpticsFano*sqrt(elateVeto*VetoLY))+0.5 );


    //cout<<"npe veto prompt  "<<prompt_npeVeto<<" late "<<del_npeVeto<<endl;
    
    ////// Cleaning of the clusters
    //int count_cl_erased = 0 ;
    //int ClToBeErased[200];
    //for(int i=0; i< int(Clusters.size()); ++i) {
    //  if(Clusters[i].Energy*TPCLY < 20) Clusters[i].npe =  ran->Poisson(Clusters[i].Energy*TPCLY)  ;
    //  else                              Clusters[i].npe =  int(ran->Gaus(Clusters[i].Energy*TPCLY , sqrt(Clusters[i].Energy*TPCLY)) + 0.5) ;
    //  
    //  if(Clusters[i].npe < 5) {
    //    ClToBeErased[count_cl_erased] = i ;
    //    count_cl_erased++;
    //  }
    //}
    //////remove very low energy clusters
    //for(int i=count_cl_erased-1;i>=0;--i) Clusters.erase(Clusters.begin() + ClToBeErased[i]);
    
    //write the clusters
    for (int i=0; i<Clusters.size() ; ++i) {
      cl_true_ene[i]   = Clusters[i].Energy; 
      cl_x[i]          = Clusters[i].X0; 
      cl_y[i]          = Clusters[i].Y0; 
      cl_z[i]          = Clusters[i].Z0; 
      cl_t[i]          = Clusters[i].T0; 
      // QUENCHING is accounted for with qAr = 0.25  ////////////
      // and the recombination too (S1quench) 
      if (Clusters[i].nucl >0.5 ) cl_nucl[i]     = Clusters[i].nucl*S1quench(Clusters[i].nucl,1.) * qarAr; 
      else cl_nucl[i]  = 0 ; 
      //////////////////////////////////////////////
      if (Clusters[i].elec >0.5 ) cl_elec[i]	 = Clusters[i].elec*S1quench(Clusters[i].elec,.21)   ; 
      else cl_elec[i]  = 0 ; 
      cl_ene[i]        = cl_elec[i]  +  cl_nucl[i] ; 
      cl_ndep[i]       = Clusters[i].nDeposits;  
      if(Clusters[i].Energy*TPCLY < 20) cl_npe[i]        = ran->Poisson(cl_ene[i]*TPCLY)  ;
      else  cl_npe[i] =  int(ran->Gaus(cl_ene[i]*TPCLY , sqrt(cl_ene[i]*TPCLY)) + 0.5) ;

    } 


    /*
    ////////////////////////////////////// 
    //    ! corrected after 
    if (nclus>0) { 
      tdrift  = ZToTDrift(cl_z[0]) ; 
        if ( theHeader.DetectorFlag <= 2 )  s1_corr = theEvent.NPE * s1_correction(tdrift) ;
        else if ( theHeader.DetectorFlag == 10 ) s1_corr = 1./get_ds20k_s1_corr(cl_z[0], sqrt(cl_x[0]*cl_x[0] + cl_y[0]*cl_y[0]));
    }
    */
    
    tdrift  = ZToTDrift(cl_z[0]) ; 

    // Fill user variables
    for(int i=0;i<theEvent.NUsers;++i) {
      INT1[i]    = theUsers[i].UserInt1;
      INT2[i]    = theUsers[i].UserInt2;
      FLOAT1[i]  = theUsers[i].UserFloat1;
      FLOAT2[i]  = theUsers[i].UserFloat2; 
      DOUBLE[i]  = theUsers[i].UserDouble; 
    }
    // Fill pe variables
    
    double f90_tile = 0.;
    double f120_tile = 0.;
    double f150_tile = 0.;
    double f200_tile = 0.;
    double f260_tile = 0.;
    fall_tile = 0.;
    double f90_light = 0.;
    double f120_light = 0.;
    double f150_light = 0.;
    double f200_light = 0.;
    double f260_light = 0.;
    double fall_light = 0.;
    
    for (int i=0;i<38;++i) pe_ch[i] =0 ; 
    s1_max_frac =  -1. ; 
    
    //SiPM channel assignment reproduced from DSTrackingAction
    for(int i=0; i<ntiles; i++) { 
      tile_s1npe[i]=0;
      tile_s2npe[i]=0;
    }
    tile_s1max=-1;
    tile_s1maxfrac=0;
    tile_bary_x=0;
    tile_bary_y=0;
    int sumBaryPE=0;
    
    for(int i=0;i<theEvent.NPE;++i) {
      pe_pmt[i]  = thePhotoElectrons[i].PMT;    
      pe_time[i] = thePhotoElectrons[i].Time; 
      if (theHeader.DetectorFlag < 3) ++pe_ch[pe_pmt[i]]  ; 
      if(pe_time[i]/1000. > 8) {
	      s2npe++;
        tile_s2npe[pe_pmt[i]]++;
        //compute charge barycenter using s2 of top tiles
        if(pe_pmt[i]>=4900) {
          sumBaryPE++;
          tile_bary_x += tile_x[pe_pmt[i]];
          tile_bary_y += tile_y[pe_pmt[i]];
        }
      }else{
	      s1npe++;
        tile_s1npe[pe_pmt[i]]++;
        if(tile_s1max<0)  tile_s1max = pe_pmt[i]; 
        else if(tile_s1npe[tile_s1max]<tile_s1npe[pe_pmt[i]])  tile_s1max = pe_pmt[i];
      }
      //compute f90
      if (theHeader.DetectorFlag < 3) { 
	if(pe_time[i]-pe_time[0]<90)
          f90_light+=1;
	if(pe_time[i]-pe_time[0]<7000)
          fall_light+=1;
	if (pe_time[i]>8000) total_s2+=1. ; 
      }  
    }	 
    tile_bary_x = tile_bary_x/(double)sumBaryPE;
    tile_bary_y = tile_bary_y/(double)sumBaryPE;
    tile_s1maxfrac = (double)tile_s1npe[tile_s1max]/(double)s1npe;
    
    double Prompt= 0 , Late=0  ; 
    if(theEvent.NPE!=0 && fall_light>0) { 
      f90npe = f90_light/fall_light;
      Prompt =  ran->Gaus(f90_light, fXX_spe_resolution * sqrt(f90_light)) ; 
      Late  =  ran->Gaus( fall_light - f90_light , fXX_spe_resolution * sqrt(fall_light - f90_light) ) ; 
      f90like =  Prompt/(Prompt+Late) ; 
    }
    
    if ( theHeader.DetectorFlag == 10 )
      s1_corr = 1./get_ds20k_s1_corr(theEvent.Position[2],sqrt(pow(theEvent.Position[0], 2)+pow(theEvent.Position[1],2))) ; 
    else  { 
      total_s1 = (Prompt+Late) ; 
      s1_corr = (Prompt+Late) * s1_correction(tdrift) ; 
    }
    
    float pe_frac[38] ; 
    for (int i=0;i<38;++i) { 
      pe_frac[i] = pe_ch[i] / float (theEvent.NPE) ; 
      if (s1_max_frac <pe_frac[i] )  s1_max_frac = pe_frac[i]  ; 
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
    
    bool isFirstPEFound = false ; 
    int FirstPEFound_index  = 0 ; 
    
    // PSD parameter for tiles 
    h_tile_top->Reset() ; 
    h_tile_bot->Reset() ;

    if (SiPM_time_jitter_stddev > 0.00000001 ) {
      for(int i=0;i<theEvent.NPH;++i) {
        thePhotons[i].Time += ran->Gaus(0,SiPM_time_jitter_stddev) ;
      }
      //re-do the ph time sorting
      std::sort(thePhotons.begin(),thePhotons.end(), cmp_photon());
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
      //number of cerenkov photons produced in water or gdwater
      if(ph_volume[i] == 17 || ph_volume[i] == 52) nphc+= thePhotons[i].PID;
      //compute f90
      if (theHeader.DetectorFlag ==10 && ph_volume[i] == 56) { 
        if (ran->Uniform () > SiPM_QE ) continue ; 
	
	
	int binx = h_tile_bot->GetXaxis()->FindBin(ph_x[i]) ;
	int biny = h_tile_bot->GetYaxis()->FindBin(ph_y[i]) ;
	if (ph_z[i] < 0) { 
          if (h_tile_bot->GetBinContent(binx, biny) == 0 ) h_tile_bot->SetBinContent( binx, biny , ran->Gaus(ph_time[i] ,tile_time_jitter_stddev) ) ; 	} else {
          if (h_tile_top->GetBinContent(binx, biny) == 0 ) h_tile_top->SetBinContent( binx, biny , ran->Gaus(ph_time[i] ,tile_time_jitter_stddev)) ; 
	}
	
	if (!isFirstPEFound)  FirstPEFound_index = i ; 
	isFirstPEFound = true ; 

        ++ds20npe; 
        total_s1+=1. ; 
	if(ph_time[i]-ph_time[FirstPEFound_index]<180)
          f90_light+=1;
	if(ph_time[i]-ph_time[FirstPEFound_index]<240)
            f150_light += 1; 
	if(ph_time[i]-ph_time[FirstPEFound_index]<200)
            f200_light+=1 ;
	if(ph_time[i]-ph_time[FirstPEFound_index]<220)
            f120_light+=1;

	if(ph_time[i]-ph_time[FirstPEFound_index]<260)
            f260_light+=1 ; 

	if(ph_time[i]-ph_time[FirstPEFound_index]<7000)
          fall_light+=1;
	else total_s2+=1 ;   
      }
     }
     
     
     // DC noise 
     int noise = 0 ; 
     if (theHeader.DetectorFlag >2 && theEvent.NPH>0) { 
       noise = ran->Poisson(SiPM_dark_cts_in_7us) ; 

       // Adding (if the tile didn't already fired) the noise to a random tile
       for (int i=0;i<noise;++i){
	 int xx = int ( ran->Uniform()*72 +0.5 ); 
	 int yy = int ( ran->Uniform()*72 +0.5 ); 
	 while (xx*xx + yy*yy > 36*36) {
           xx = int ( ran->Uniform()*72 +0.5 ); 
           yy = int ( ran->Uniform()*72 +0.5 ); 
	 }
	 if (ran->Uniform()>=0.5) {
           if (h_tile_bot->GetBinContent(xx, yy) == 0 ) h_tile_bot->SetBinContent( xx, yy , ran->Uniform()*7e-6  ) ; 
	 } else {
           if (h_tile_top->GetBinContent(xx, yy) == 0 ) h_tile_top->SetBinContent( xx, yy , ran->Uniform()*7e-6 ) ; 
	 }
       }

       //loop over the tiles to find the first pe
       if(theEvent.NPH!=0 && theHeader.DetectorFlag >2 ) {
	 for (int i=1;i<=h_ntiles;++i){
           for (int j=1;j<=h_ntiles;++j){
	     double tt = 0 ; 
	     for (int m=0;m<2;++m) {
	       if (m==0) tt=h_tile_top->GetBinContent(i,j); 
	       else      tt=h_tile_bot->GetBinContent(i,j); 
	       if (tt==0.) continue ; 
	       if (tt < fX_t0 ) fX_t0 = tt ; 
             }  
           }
	 }
       }

       //loop over the tiles to compute the variables
       if(theEvent.NPH!=0 && theHeader.DetectorFlag >2 ) {
	 for (int i=1;i<=h_ntiles;++i){
           for (int j=1;j<=h_ntiles;++j){
	     double tt = 0 ; 
	     for (int m=0;m<2;++m) {

	       if (m==0) tt=h_tile_top->GetBinContent(i,j); 
	       else      tt=h_tile_bot->GetBinContent(i,j); 

	       if (tt==0 || tt-fX_t0 > 8000) continue ; 
	       else ++fall_tile ; 
	       if (tt-fX_t0< 220 )  f120_tile+=1. ; 
	       if (tt-fX_t0< 240 )  f150_tile+=1. ; 
	       if (tt-fX_t0< 180 )   f90_tile+=1. ; 
	       if (tt-fX_t0< 200 )  f200_tile+=1. ; 
	       if (tt-fX_t0< 260 )  f260_tile+=1. ; 
	     }
           }
	 }
       }
       f90like_tile = f90_tile/fall_tile ; 
       f120like_tile = f120_tile/fall_tile ; 
       f150like_tile = f150_tile/fall_tile ; 
       f200like_tile = f200_tile/fall_tile ; 
       f260like_tile = f260_tile/fall_tile ; 

       double fall_light_sm ;       
       double fall_light_sm_noise;  

       if(theEvent.NPH!=0 && theHeader.DetectorFlag >2 ) {
	 fall_light_sm = ran->Gaus(fall_light,   fXX_spe_resolution * sqrt(fall_light) ) ; 
	 fall_light_sm_noise = fall_light_sm + ran->Gaus(noise, fXX_spe_resolution * sqrt(noise) ) ; 
	 double f90_light_sm  = ran->Gaus(f90_light ,   fXX_spe_resolution * sqrt(f90_light ) ) ; 
	 double f120_light_sm = ran->Gaus(f120_light,   fXX_spe_resolution*   sqrt(f120_light) ) ; 
	 double f150_light_sm = ran->Gaus(f150_light,   fXX_spe_resolution *  sqrt(f150_light) ) ; 
	 double f200_light_sm = ran->Gaus(f200_light,   fXX_spe_resolution *  sqrt(f200_light) ) ; 
	 double f260_light_sm = ran->Gaus(f260_light,   fXX_spe_resolution *  sqrt(f260_light) ) ; 

	 f90like = f90_light_sm/fall_light_sm;
	 f120like = f120_light_sm/fall_light_sm;
	 f150like = f150_light_sm/fall_light_sm;
	 f200like = f200_light_sm/fall_light_sm;
	 f260like = f260_light_sm/fall_light_sm;
	 f90like_noise = f90_light_sm/(fall_light_sm_noise);
	 f120like_noise = f120_light_sm/(fall_light_sm_noise);
	 f150like_noise = f150_light_sm/(fall_light_sm_noise);
	 f200like_noise = f200_light_sm/(fall_light_sm_noise);
	 f260like_noise = f260_light_sm/(fall_light_sm_noise);
	 total_s1 = fall_light_sm_noise ; 
	 total_s2 =  ran->Gaus(total_s2,   fXX_spe_resolution *  sqrt(total_s2) ); 
       }
     double tmp = (total_s1)*s1_corr ; // TODO get_ds20k_s1_corr(cl_z[0], sqrt(cl_x[0]*cl_x[0] + cl_y[0]*cl_y[0]));
     s1_corr     = tmp ; 
     ds20npe += noise ; 
   }
    
    // Fill the event
    dstree->Fill();
    
  }
  
  //TFile *ff = new TFile(rootfile.c_str(),"recreate");
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
 * $Log: g4rootered_full.C,v $
 * Revision 1.18  2015/10/16 17:05:43  dfranco
 * bug on neutron fixed
 *
 * Revision 1.17  2015/10/15 10:39:53  pagnes
 * useless variables removed from the output tree
 *
 * Revision 1.16  2015/10/15 09:24:39  dfranco
 * all units updated: position in cm and energy in keV
 *
 * Revision 1.15  2015/10/14 10:40:51  dfranco
 * added total energy deposit variables
 *
 * Revision 1.14  2015/10/14 09:10:09  dfranco
 * fixed two bugs in the clustering
 *
 * Revision 1.13  2015/10/13 10:28:03  dfranco
 * fix nuclear and electron recoils fractions in the clusters
 *
 * Revision 1.12  2015/10/07 14:07:49  pagnes
 * forgotten cout removed
 *
 * Revision 1.11  2015/10/07 13:50:19  pagnes
 * bug fixed in tdrift calculation
 *
 * Revision 1.10  2015/10/07 09:55:35  pagnes
 * data-like variables added: s1_max_frac, s1_corr and tdrift (for evts with 1 clust or more)
 *
 * Revision 1.9  2015/09/03 10:09:40  pagnes
 * change z clustering max dist from command line (clz=###)
 *
 * Revision 1.8  2015/01/22 09:56:46  dfranco
 * add f90
 *
 * Revision 1.7  2015/01/20 13:01:29  dfranco
 * add f90
 *
 * Revision 1.6  2015/01/20 10:57:12  dfranco
 * compute f90
 *
 * Revision 1.5  2014/11/07 15:22:32  dfranco
 * move TPC to right position for clustering
 *
 * Revision 1.4  2014/11/03 15:37:10  dfranco
 * update the clustering algorithm
 *
 * Revision 1.3  2014/10/23 11:40:20  dfranco
 * fix bug in clustering
 *
 * Revision 1.2  2014/07/09 13:06:12  pagnes
 * Generators in materials fixed
 *
 * Revision 1.1  2014/05/08 10:59:17  pagnes
 * Scintillator Index added in binary header
 *
 * Revision 1.14  2013/10/20 16:30:20  swesterd
 * updated the veto scintillator optical properties based on Aldos measurements
 *
 * Revision 1.13  2013/10/01 06:26:25  swesterd
 * added waveform averager
 *
 * Revision 1.12  2013/08/27 04:07:02  swesterd
 * some fine tuning of the boron scintillator kB and scint yield, and some modifications to the DSG2 geometry
 *
 * Revision 1.11  2013/08/20 03:25:53  swesterd
 * added G2 TPC geoemtry (not complete) and added monoenergetic energy distribution to generator
 *
 * Revision 1.10  2013/08/06 13:58:22  dfranco
 * Added 3 variables to store the deposited energy in LAr, scintillator, and water. The last two are not yet implemented. g4rooter has been updated with 3 new variables: tpcene, vetoene, and muene
 *
 * Revision 1.9  2013/07/24 09:49:03  dfranco
 * Added S1 and S2 equivalent energies, the material index for each deposit, the command killS1S2 to kill photons and electrons generated by DSLight (after storing the equivalent energies)
 *
 * Revision 1.8  2013/06/13 10:04:58  dfranco
 * update g4rooter variable names
 *
 * Revision 1.7  2013/06/04 01:02:31  swesterd
 * other than the optical boundary of the trunks, the veto optics appear to be complete and up and running...modulo whatever I may have missed...
 *
 * Revision 1.6  2013/05/07 23:06:32  swesterd
 * added optical boundaries and Lumirror in the veto
 *
 * Revision 1.5  2013/04/04 09:22:54  dfranco
 * added variables to dstree
 *
 * Revision 1.4  2013/04/04 09:07:47  dfranco
 * added deposit step length, visibile (quenched) energy, original energy to dstree
 *
 * Revision 1.3  2013/03/22 13:23:21  dfranco
 * deposit times sorted
 *
 * Revision 1.2  2013/03/20 09:54:28  dfranco
 * New version of g4ds
 *
 */
