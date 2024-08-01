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

double squareh[8], squarev[8] , xnsp[1000], ynsp[1000]; 

bool isLightTree = 1 ; 
//constants 
const double PI    = TMath::Pi(); 
const double TWOPI = 2*TMath::Pi();  ; 
const int MAXDAUGHTER = 10000;
const int MAXDEPOSIT  = 100000;
const int MAXNPH      = 100000;
const int MAXNPE      = 100000;
const int MAXUSER     = 1000;


vector<float> vpdf, vpdf_time, vtheta, vlenght, vjitter_prob, vjitter_time;

//clustering 
bool isToBeAdded;  
int isHcap, isCucap, isArcap, isCcap, isSicap, isOcap, isGdcap, pdgCap;
double cap_time      = 0. ; 
// cluser params
float dist_max_xy = 500. ;    //  cm 
float dist_max_z  = 0.2 ;  //  cm 
float deltaT_max  = dist_max_z*10000.  ;  // ns 
double VetoAcquisitionWindow = 2000e3 ; // ns

int ScintillatorIndex ; 
int LArIndex; 
bool reverse (int i,int j) { return (i>j); }
float TPCHeight ; 
float TPCEdge   ; 

////////////////////////////////////////////////////////
//compute theta angle
////////////////////////////////////////////////////////
double Angle2D(double x1, double y1, double x2, double y2)
{
   double dtheta,theta1,theta2;

   theta1 = atan2(y1,x1);
   theta2 = atan2(y2,x2);
   dtheta = theta2 - theta1;
   
   while (dtheta > PI)
      dtheta -= TWOPI;
   while (dtheta < -PI)
      dtheta += TWOPI;

   return(dtheta);
}

////////////////////////////////////////////////////////
//check if the (x,y) point is inside the octagon
////////////////////////////////////////////////////////
bool InsideOctagon(float _radius , double pxx, double pyy)
{
   int i;
   double angle=0.;
   double p1h, p1v,p2h, p2v;

   double px = pxx *cos (PI/8.) - pyy * sin(PI/8) ; 
   double py = pxx *sin (PI/8.) + pyy * cos(PI/8) ; 

   for (int _i=0;_i<8;++_i) {
      p1h = _radius*squareh[_i] - px;
      p1v = _radius*squarev[_i] - py;
      p2h = _radius*squareh[(_i+1)%8] - px;
      p2v = _radius*squarev[(_i+1)%8] - py;
      
      angle += Angle2D(p1h,p1v,p2h,p2v);
   }
   if (fabs(angle+PI/19) < PI)
      return kFALSE;
   else
      return kTRUE;
}


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
  else return   -(Z - TPCHeight/2. )/240  * 2.5e3  ;  //assuming 2.5 ms maximum drift time in DS20k  
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

////////////////////////////////////////////////////////////////////////////////////
// Quenching for TPC - we need to take into account the increase of Ly at low energy (due to higher recombination).
// With this function we calculate the amount of energy (in keV) that goes in S1 
// The recombination is supposed to be the same for both ER and NR, but the excitation/ionization ratio is different
// update: fraction of S1 that goes in scintillation, assumng ER recombination probability with a user defined alpha = Nex/Ni
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


// the quenching in liquid scintillator uses the above fuctions
double getqfactor(double ene,int pdg) {
  //return 1;
  
  if(pdg - 1e9 == 20040) return qalpha(ene); 
  else if(pdg == 11 || pdg == -11 || pdg == 22 )  qelectron(ene);
  else if(pdg == 2212 || pdg == 2112) return   qproton(ene);
  else if(pdg - 1e9 > 30000 ) return  qcarbon(ene);
  else if(pdg - 1e9 > 30000 ) return 0;
  
  return 1 ;

}


// the quenching in plastic is approximated as follows
double getqfactor_PS(double ene,int pdg) {
  //return 1;
  
  if(pdg - 1e9 == 20040) return 1./25; 
  else if(pdg == 11 || pdg == -11 || pdg == 22 ) return  1; 
  else if(pdg == 2212 || pdg == 2112) return  1/3.5 ;  
  else if(pdg - 1e9 > 30000 ) return 0 ; 

  return 1 ;

}

// the quenching in LAr (outside the TPC) is approximated as follows
double getqfactor_LAr(double ene,int pdg) {
  
  if(pdg - 1e9 == 20040) return 1.;  
  else if(pdg == 11 || pdg == -11 || pdg == 22 ) return  1; 
  else if(pdg == 2212 || pdg == 2112) return  1/4. ;  
  else if(pdg - 1e9 > 30000 ) return 1./4. ; 
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


  float myR, myZ;

  bool FVcut = false;
  bool NRcut = false;

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
    if (!argument.find("myR=")) {
    argument.erase(0,4);
    myR = (float) atof(argument.c_str());
    cout << "   myR set to: " << myR << endl;
    }   
    if (!argument.find("myZ=")) {
    argument.erase(0,4);
    myZ = (float) atof(argument.c_str());
    cout << "   myZ set to: " << myZ << endl;
    }
    if (!argument.find("FVcut=")) {
    argument.erase(0,6);
    FVcut = (bool)atoi(argument.c_str());
    cout << "   FVcut (nclusNR && isFV30) set to: " << FVcut << endl;
    }
    if (!argument.find("NRcut=")) {
    argument.erase(0,6);
    NRcut = (bool)atoi(argument.c_str());
    cout << "   NRcut (nclusNR) set to: " << NRcut << endl;
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
  TPCHeight = int (theHeader.Rate * 1e-3) ; 
  TPCEdge   = int (theHeader.Rate - TPCHeight*1e3 ) ; 
  
  TPCHeight = 350.5 ; 
  TPCEdge   = 145.1 ; 

   cout << "DetectorFlag: " << theHeader.DetectorFlag << endl ; 
   if (theHeader.DetectorFlag == 10 || theHeader.DetectorFlag == 11) { 
     cout << "TPC dimensions: " << TPCHeight << " cm tall, " << TPCEdge <<" cm edge" << endl ;      
     cout << "Scintillator material: " << ScintillatorIndex << endl ; 
   } 
  // event extra variables
  float radius = 0;

  // daughter variables
  int  * Did = new int [MAXDAUGHTER]; 
  int * Dpdg = new int [MAXDAUGHTER]; 
  int * Dpid = new int [MAXDAUGHTER]; 
  int * Dprocess = new int [MAXDAUGHTER];
  float *Dene = new float [MAXDAUGHTER]; 
  float *  Dx = new float [MAXDAUGHTER]; 
  float *  Dy = new float [MAXDAUGHTER]; 
  float *  Dz = new float [MAXDAUGHTER]; 
  float *  Dr = new float [MAXDAUGHTER]; 
  float *  Dpx= new float [MAXDAUGHTER]; 
  float *  Dpy= new float [MAXDAUGHTER]; 
  float *  Dpz= new float [MAXDAUGHTER];
  double * Dtime= new double [MAXDAUGHTER];
  
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


//bar variables                                                                                                                                                                         
  int nbars = 0;
  int nbars_prompt = 0;
  int nbars_prompt10 = 0;
  int nbars_prompt100 = 0;
  int nbars_late = 0;
  int nbars_late10 = 0;
  int nbars_late100 = 0;
  float bar_tot_energy =0;
  float bar_late_energy =0;
  double tdep_min = 0;
  double tdep_max=0;
  float * bar_prompt_ene = new float[MAXDEPOSIT];                                                                                 
  float * bar_late_ene   = new float[MAXDEPOSIT];                                                                                 
  float * bar_step = new float[MAXDEPOSIT];

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
  int *isAcrylicVeto = new int[MAXDEPOSIT] ; 
  int *isOuterVeto = new int[MAXDEPOSIT] ; 
  int ds20npe ; 
  float total_s1, total_s2;
  int ndeptpc, nclus, nclus_thr, cl_ndep[MAXDEPOSIT], ab_mat, cl_npe[MAXDEPOSIT];
  float s1ene,s2ene, veto_visene, mu_visene, ene, qene, qnpe, tpcene, vetoene,
    muene , depTPCTot , depVeto, npeTPC, npeTPC400, depTPC400 ,
    depVeto70, eneTPC, eneVeto, eneslideVeto, timeslideVeto, npeslideVeto, depTot; 
  double ab_time, timeVeto;
  float cl_ene[MAXDEPOSIT], cl_true_ene[MAXDEPOSIT],cl_x[MAXDEPOSIT], cl_y[MAXDEPOSIT], cl_z[MAXDEPOSIT], cl_t[MAXDEPOSIT]; 
  int fall_tile ; 
  float cl_nucl[MAXDEPOSIT], cl_elec[MAXDEPOSIT];
  float cap_y, cap_x, cap_z , energyinNS;

  int isFV30[MAXDEPOSIT], isFV5[MAXDEPOSIT]; 
  int isFV10 , isFV15 , isFV20, isFV35; 
  
  double sq  = sqrt(2)/2. ; 
  squareh[0] = 1 ; 
  squarev[0] = 0 ; 
  squareh[1] = sq ; 
  squarev[1] = sq ; 
  squareh[2] = 0 ; 
  squarev[2] = 1 ; 
  squareh[3] = -sq ; 
  squarev[3] = sq ; 
  squareh[4] = -1 ; 
  squarev[4] = 0 ; 
  squareh[5] = -sq ; 
  squarev[5] = -sq; 
  squareh[6] = 0 ; 
  squarev[6] = -1 ; 
  squareh[7] = sq ; 
  squarev[7] = -sq ; 

  double tpromptVeto, tdelayedVeto,tlateVeto;
  //int    epromptVeto, edelayedVeto,elateVeto;
  float tdrift , s1_corr ; 
  float   prompt_qdepMat[85] , late_qdepMat[85]	 , prompt_depMat[85]   , late_depMat[85]	 ,late_timeMat[85]	 ; 
  float   prompt_eneVeto_PS ,prompt_eneVeto_LS , prompt_eneVeto_Ar , prompt_eneVeto  ;   
  float late_eneVeto_PS ; 
  float late_eneVeto_LS ; 
  float late_eneVeto_Ar ; 
  float late_eneVeto	;  

  float ene_shell_0, ene_shell_1, ene_shell_2, ene_shell_3, ene_shell_4, ene_shell_5, ene_shell_6, ene_shell_7;

 
  float  mat_energy_fraction[100];
  float cap_gamma_ene ; 
  int cap_gamma_mult ; 
  int nphc;
  TTree *nentree = new TTree("nentree","The G4DS Tree Counting Events");
  nentree->Branch("nentries",             &counter,         "nentries/I");

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
//  if (!isLightTree ) { 
    dstree->Branch("ene",            &ene,		      "ene/F");        
    dstree->Branch("r",              &radius,		      "radius/F");     
    dstree->Branch("px",             &theEvent.Direction[0],    "px/F");	       
    dstree->Branch("py",             &theEvent.Direction[1],    "py/F");	       
    dstree->Branch("pz",             &theEvent.Direction[2],    "pz/F");	       
    //dstree->Branch("bx",             &theEvent.CenterOfMass[0], "bx/F");	       
    //dstree->Branch("by",             &theEvent.CenterOfMass[1], "by/F");	       
    //dstree->Branch("bz",             &theEvent.CenterOfMass[2], "bz/F");	       
    dstree->Branch("munpe" ,         &theEvent.MuNPE ,	      "munpe/I");      
    dstree->Branch("vnpe",           &theEvent.VetoNPE ,        "vnpe/I");       
    dstree->Branch("nph",            &theEvent.NPH,	      "nph/I");        
    dstree->Branch("s1",            &total_s1 ,	      "s1/F");        
    dstree->Branch("s2",            &total_s2 ,	      "s2/F");        
    dstree->Branch("s1_corr",        &s1_corr ,	              "s1_corr/F");	
    dstree->Branch("npe",            &theEvent.NPE ,	      "npe/I");        
  //}
  dstree->Branch("tdrift",         &tdrift ,	              "tdrift/F");	   
    
  dstree->Branch("ndaughters",     &theEvent.NDaughters,      "ndaughters/I"); 
  dstree->Branch("ndeposits" ,     &theEvent.NDeposits,       "ndeposits/I");  
  dstree->Branch("ndepositsTPC" ,  &ndepoTPC,                 "ndepositsTPC/I");
  dstree->Branch("nusers",         &theEvent.NUsers,	      "nusers/I");
  /*
   * dstree->Branch("nbars" ,         &nbars,                     "nbars/I");    
  dstree->Branch("nbars_prompt" ,         &nbars_prompt,                     "nbars_prompt/I");    
  dstree->Branch("nbars_prompt10" ,         &nbars_prompt10,                     "nbars_prompt10/I");    
  dstree->Branch("nbars_prompt100" ,         &nbars_prompt100,                     "nbars_prompt100/I");    
  dstree->Branch("nbars_late" ,         &nbars_late,                     "nbars_late/I");    
  dstree->Branch("nbars_late10" ,         &nbars_late10,                     "nbars_late10/I");    
  dstree->Branch("nbars_late100" ,         &nbars_late100,                     "nbars_late100/I");    
  //dstree->Branch("nphc",         &nphc,	      "nphc/I");    
  //dstree->Branch("isAcrylicVeto",            isAcrylicVeto,	      "isAcrylicVeto[ndeposits]/I");        
  //dstree->Branch("isOuterVeto",            isOuterVeto,	      "isOuterVeto[ndeposits]/I");        
  */
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
  if (!isLightTree )  {
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
  

  } 

  dstree->Branch("mat_fraction",  &mat_energy_fraction[0],	          "mat_fraction[100]/F") ;
  
  dstree->Branch("prompt_qdepMat",  prompt_qdepMat  ,	          "prompt_qdepMat[85]/F") ;
  dstree->Branch("late_qdepMat",   late_qdepMat      ,	          "late_qdepMat[85]/F") ;
  dstree->Branch("prompt_depMat",  prompt_depMat    ,	          "prompt_depMat[85]/F") ;
  dstree->Branch("late_depMat",    late_depMat        ,	          "late_depMat[85]/F") ;
  dstree->Branch("late_timeMat",   late_timeMat      ,	          "late_timeMat[85]/F") ;

  
  //dstree->Branch("prompt_eneVeto_PS",       &prompt_eneVeto_PS  ,	          "prompt_eneVeto_PS/F") ;
  //dstree->Branch("prompt_eneVeto_LS",       &prompt_eneVeto_LS  ,	          "prompt_eneVeto_LS/F") ;
  dstree->Branch("prompt_eneVeto_Ar",       &prompt_eneVeto_Ar  ,	          "prompt_eneVeto_Ar/F") ;
  //dstree->Branch("prompt_eneVeto",       &prompt_eneVeto  ,	          "prompt_eneVeto/F") ;
  //dstree->Branch("late_eneVeto_PS",       &late_eneVeto_PS  ,	          "late_eneVeto_PS/F") ;
  //dstree->Branch("late_eneVeto_LS",       &late_eneVeto_LS  ,	          "late_eneVeto_LS/F") ;
  dstree->Branch("late_eneVeto_Ar",       &late_eneVeto_Ar  ,	          "late_eneVeto_Ar/F") ;
  //dstree->Branch("late_eneVeto",       &late_eneVeto  ,	          "late_eneVeto/F") ;
    
  dstree->Branch("pdgCap",	 &pdgCap,	     "pdgCap/I") ;
  dstree->Branch("cap_x",       &cap_x,                    "cap_x/F");      
  dstree->Branch("cap_y",       &cap_y,                    "cap_y/F");      
  dstree->Branch("cap_z",       &cap_z,                    "cap_z/F") ;   
  dstree->Branch("cap_gamma_ene",	   &cap_gamma_ene,			 "cap_gamma_ene/F") ;   
  dstree->Branch("cap_gamma_mult",	   &cap_gamma_mult,			 "cap_gamma_mult/I") ;   
  dstree->Branch("cap_time",       &cap_time,                    "cap_time/D") ;     

  dstree->Branch("depVeto",	 &depVeto,	     "depVeto/F") ;

  int nclusNR, nclusER, goodNRidx;   // nclusNR ==1 && isFV30[goodNRidx] 
  float energyER = 0. ; 
  float energyERinel = 0;
  dstree->Branch("goodNRidx",         &goodNRidx,	          "goodNRidx/I") ;
  dstree->Branch("nclusNR",         &nclusNR,	          "nclusNR/I") ;
  dstree->Branch("nclusER",         &nclusER,	          "nclusER/I") ;
  dstree->Branch("energyER",         &energyER,	          "energyER/F") ;
  dstree->Branch("energyERinel",         &energyERinel,           "energyERinel/F") ;

  dstree->Branch("nclus",         &nclus,                 "nclus/I") ;

  dstree->Branch("isFV5",         &isFV5,	          "isFV5[nclus]/I") ;
  dstree->Branch("isFV10",        &isFV10,	          "isFV10/I") ;
  dstree->Branch("isFV15",         &isFV15,	          "isFV15/I") ;
  dstree->Branch("isFV20",         &isFV20,	          "isFV20/I") ;
  dstree->Branch("isFV30",         isFV30,	          "isFV30[nclus]/I") ;
  dstree->Branch("isFV35",         &isFV35,	          "isFV35/I") ;
   
  dstree->Branch("userint1",    INT1,			"int1[nusers]/I");  
  dstree->Branch("userint2",    INT2,			"int2[nusers]/I");
  dstree->Branch("userfloat1",  FLOAT1,			"float1[nusers]/F")  ; 
  dstree->Branch("userfloat2",  FLOAT2,			"float2[nusers]/F");      	   
  dstree->Branch("userdouble0", DOUBLE,			"double0[nusers]/D");  

  //dstree->Branch("nclus",         &nclus,	          "nclus/I") ;
  dstree->Branch("nclus_thr",         &nclus_thr,	          "nclus_thr/I") ;
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
  
    
  if (!isLightTree ) {
    dstree->Branch("pe_time",     pe_time,  		"pe_time[npe]/D");   	   
    dstree->Branch("pe_pmt",      pe_pmt,			"pe_pmt[npe]/I");   
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

    pdgCap  = 0 ; 
    cap_z   = -1e5 ; 
    cap_x   = -1e5 ; 
    cap_y   = -1e5 ; 
    cap_time   = 0. ; 
    cap_gamma_mult = 0 ; 
    cap_gamma_ene   = 0. ; 
    
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
      
      if (Dprocess[i] == 4131 ) {  

	if (Dpdg[i] > 1e9 ) pdgCap = int( Dpdg[i] - 1e9 )/ 1e4; 
	cap_x = 1.*Dx[i] ; 
	cap_y = 1.*Dy[i] ; 
	cap_z = 1.*Dz[i] ;
	cap_time = theDaughters[i].Time ; 
	
	if (Dpdg[i] == 22 ) { 
	  cap_gamma_mult+= 1 ; 
	  cap_gamma_ene+=Dene[i] ; 
	} 
      }
    }

    // Fill deposits variables
    float TPCLY   = 7.100;   //photons per keV
    float qvalpha = 0.027;             // alpha in veto (40keVee for 1.47MeV alphas)
    float qvli    = 0.0 ;               // Li in veto  -  TODO         
    float qvp     = 0.3333;            // p in veto 
    float qarAr   = 0.25 ;             // Ar in LAr
    float quench ; 
  
    isFV10 = 0 ; 
    isFV15 = 0 ; 
    isFV20 = 0 ; 
    isFV35 = 0 ; 
       
 
    depTPCTot = 0;
    depVeto=0;
    depTot=0;
    npeTPC=0;
    npeTPC400=0;
    depVeto70=0;
    ndeptpc = 0 ;
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
    energyinNS = 0 ; 
    
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
      tot_dep_ene  +=  dep_ene[i];
      if(dep_mat[i] == 8) depTPCTot += dep_ene[i];
      if(dep_mat[i] == 8  ) ndepoTPC++;
      
      isAcrylicVeto[i] = ( (dep_mat[i]==7 || dep_mat[i]== 65 ||  dep_mat[i]== 74) && ((sqrt(dep_x[i]*dep_x[i]+dep_y[i]*dep_y[i])>150) || (abs(dep_z[i]) > 130))) ; 
      isOuterVeto[i]   = (dep_mat[i]==66 && ((sqrt(dep_x[i]*dep_x[i]+dep_y[i]*dep_y[i])<330) && (abs(dep_z[i]) < 305))) ; 
      
      //Passive argon, non-scintillating, around the TPC
      if(dep_mat[i] == 10 ) {
	if(dep_pdg[i]< 30) { 
	  energyinNS += dep_ene[i];
	}
	else { 
	  energyinNS += dep_ene[i] * 0.25;
	}
      }
      
      //Some clusters are produced (especially by n captures) after the end of the aquisition window. We should account for this time difference
      //if(dep_time[i] > 400e3) break  ; 
      
      
      // CLUSTERING ALGORITHM
      if (dep_mat[i] != LArIndex ) continue ;     
      
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
    
    /////////////////////////////////////////////////////////////////////////////////////
    // add a cut on the cluster energy
    /////////////////////////////////////////////////////////////////////////////////////
    int ClToBeErased[200];
    int n_to_be_erased = 0 ; 
  
    nclusNR      = 0 ; 
    nclusER      = 0 ; 
    goodNRidx    = 0 ; 
    energyER     = 0 ;
    energyERinel = 0 ;
    for(int i=0; i< int(Clusters.size()); ++i) {
      if (Clusters[i].nucl>Clusters[i].elec ) { // NR
        nclusNR += 1  ; 
    	if (Clusters[i].Energy < 0.4  ) { 
    	    ClToBeErased[i] = 1 ; 
	    ++n_to_be_erased ; 
    	} else { 
          ClToBeErased[i] = 0 ; 
          goodNRidx       = i ; 
        }
      } else {  //ER 
        nclusER += 1 ; 
	if ( Clusters[i].T0 > cap_time - 2.  ) energyER += Clusters[i].Energy;
	else energyERinel += Clusters[i].Energy;
      }
    }

    /////////////////////////////////////////////////////////////////////////////////////
    // to be discussed before application 
    /////////////////////////////////////////////////////////////////////////////////////
    //for(int i=0; i< int(Clusters.size()); ++i) {
    //    if ( ClToBeErased[i] )  Clusters.erase(Clusters.begin() + i);
    //}     
    
    nclus   = int(Clusters.size());
    nclusNR = nclusNR - n_to_be_erased ;
  
    double T0;
    if(nclus!=0)
      T0 = Clusters[0].T0 ;
    else
      T0 = 0;
    
    
    //empty material list  
    for (int i=0;i<85;++i) { 
      prompt_qdepMat[i] = 0 ; 
      late_qdepMat[i] = 0 ; 
      prompt_depMat[i] = 0 ; 
      late_depMat[i] = 0 ; 
      late_timeMat[i] = -1e5 ;       
    }


      ene_shell_0 = 0;
      ene_shell_1 = 0;
      ene_shell_2 = 0;
      ene_shell_3 = 0;
      ene_shell_4 = 0;
      ene_shell_5 = 0;
      ene_shell_6 = 0;
      ene_shell_7 = 0;


      //compute energy deposits in all the materials, summing energy in prompt and late windows
      for(int k=0;k<theEvent.NDeposits;++k) {

      double deltaT = dep_time[k] - T0 ;
      //same as before! some events are otside the acquisition window. How to treat them? 
      if( deltaT > VetoAcquisitionWindow ) break ;
      
      //set quenching factors for energy deposits in active materials
      double quench = 1. ; 
      if (dep_mat[k] == 8 || dep_mat[k]==10 || dep_mat[k]==72 ||  dep_mat[k]==66) 
        quench  = getqfactor_LAr(dep_ene[k], dep_pdg[k]) ; 
      else if (dep_mat[k] == 2 || dep_mat[k] == 62 || dep_mat[k] == 63 || dep_mat[k] == 64 ) 
        quench  = getqfactor(dep_ene[k], dep_pdg[k]) ; 
      else if (dep_mat[k] == 7 || dep_mat[k] == 65 || dep_mat[k] == 68 || dep_mat[k] == 69|| dep_mat[k] == 74)  
        quench  = getqfactor_PS(dep_ene[k], dep_pdg[k]) ; 
      else quench = 0 ;  //if passive material, set to 0
            
      //prompt window is arbitrarly chosen to be [-200,200] ns around the TPC signal
      if (deltaT < 200 && deltaT > -200) { 
        prompt_depMat[dep_mat[k]]  += dep_ene[k] ;; 
        prompt_qdepMat[dep_mat[k]] += quench * dep_ene[k] ; 	
      }
      
      //late window is [200,MatAcquisitionWindow ] ns after the TPC signal
      //late time is the time of the last deposit (after a capture, should be very close to the capture time) 
      if (deltaT > 200 ) { 
        late_depMat[dep_mat[k]]  += dep_ene[k] ;; 
        late_qdepMat[dep_mat[k]] += quench * dep_ene[k] ; 	
	late_timeMat[dep_mat[k]] = dep_time [k] ; 
      }

    
     if (dep_mat[k] == 66) {

       ene_shell_0 += dep_ene[k];

       if (sqrt(dep_x[k]*dep_x[k] + dep_y[k]*dep_y[k]) < myR + 10 &&  abs(dep_z[k]) < myZ +10) {       
       ene_shell_1 += dep_ene[k];
       }
       if (sqrt(dep_x[k]*dep_x[k] + dep_y[k]*dep_y[k]) >= myR + 10 && sqrt(dep_x[k]*dep_x[k] + dep_y[k]*dep_y[k]) < myR + 15 && abs(dep_z[k]) >= myZ + 10 && abs(dep_z[k]) < myZ +15) {
       ene_shell_2 += dep_ene[k];
       }
       if (sqrt(dep_x[k]*dep_x[k] + dep_y[k]*dep_y[k]) >= myR + 15 && sqrt(dep_x[k]*dep_x[k] + dep_y[k]*dep_y[k]) < myR + 20 && abs(dep_z[k]) >= myZ + 15 && abs(dep_z[k]) < myZ +20) {
       ene_shell_3 += dep_ene[k];
       }
       if (sqrt(dep_x[k]*dep_x[k] + dep_y[k]*dep_y[k]) >= myR + 20 && sqrt(dep_x[k]*dep_x[k] + dep_y[k]*dep_y[k]) < myR + 25 && abs(dep_z[k]) >= myZ + 20 && abs(dep_z[k]) < myZ +25) {
       ene_shell_4 += dep_ene[k];
       }
       if (sqrt(dep_x[k]*dep_x[k] + dep_y[k]*dep_y[k]) >= myR + 25 && sqrt(dep_x[k]*dep_x[k] + dep_y[k]*dep_y[k]) < myR + 30 && abs(dep_z[k]) >= myZ + 25 && abs(dep_z[k]) < myZ +30) {
       ene_shell_5 += dep_ene[k];
       }
       if (sqrt(dep_x[k]*dep_x[k] + dep_y[k]*dep_y[k]) >= myR + 30 && sqrt(dep_x[k]*dep_x[k] + dep_y[k]*dep_y[k]) < myR + 35 && abs(dep_z[k]) >= myZ + 30 && abs(dep_z[k]) < myZ +35) {
       ene_shell_6 += dep_ene[k];
       }
       if (sqrt(dep_x[k]*dep_x[k] + dep_y[k]*dep_y[k]) >= myR + 35 && sqrt(dep_x[k]*dep_x[k] + dep_y[k]*dep_y[k]) < myR + 40 && abs(dep_z[k]) >= myZ + 35 && abs(dep_z[k]) < myZ +40) {
       ene_shell_7 += dep_ene[k];
       }

     }
     
    } //end of deposits loop
    
    prompt_eneVeto_PS = 0 ; 
    prompt_eneVeto_LS = 0 ;  
    prompt_eneVeto_Ar = 0 ;  
    prompt_eneVeto  = 0 ;  
    late_eneVeto_PS =  0 ; 
    late_eneVeto_LS =  0 ; 
    late_eneVeto_Ar =  0 ; 
    late_eneVeto    =   0 ; 
  
    //set some aliases
    prompt_eneVeto_PS = prompt_qdepMat[65] + prompt_qdepMat[68] + prompt_qdepMat[69] +prompt_qdepMat[74]   ; 
    prompt_eneVeto_LS = prompt_qdepMat[2] + prompt_qdepMat[62] + prompt_qdepMat[63] + prompt_qdepMat[64] ; 
    prompt_eneVeto_Ar = prompt_qdepMat[72] + prompt_qdepMat[66] ; 
    prompt_eneVeto    = prompt_eneVeto_PS + prompt_eneVeto_LS + prompt_eneVeto_Ar ; 
    late_eneVeto_PS = late_qdepMat[65] + late_qdepMat[68] + late_qdepMat[69] +late_qdepMat[74]; 
    late_eneVeto_LS = late_qdepMat[2] + late_qdepMat[62] + late_qdepMat[63]+ late_qdepMat[64] ; 
    late_eneVeto_Ar = late_qdepMat[72] + late_qdepMat[66] ; 
    late_eneVeto    = late_eneVeto_PS + late_eneVeto_LS + late_eneVeto_Ar ; 
    
    
    
    //write the clusters
    for (int i=0; i<Clusters.size() ; ++i) {
      cl_true_ene[i]   = Clusters[i].Energy; 
      cl_x[i]          = Clusters[i].X0; 
      cl_y[i]          = Clusters[i].Y0; 
      cl_z[i]          = Clusters[i].Z0; 
      cl_t[i]          = Clusters[i].T0; 
      // QUENCHING is accounted for with qAr = 0.25  ////////////
      // and the recombination too (S1quench) 
      cl_ene[i]        = Clusters[i].nucl * qarAr + Clusters[i].elec;  //restore energy and not S1 energy 
      if (Clusters[i].nucl >0.5 ) cl_nucl[i]     = Clusters[i].nucl*S1quench(Clusters[i].nucl,1.) * qarAr; 
      else cl_nucl[i]  = 0 ; 
      //////////////////////////////////////////////
      if (Clusters[i].elec >0.5 ) cl_elec[i]	 = Clusters[i].elec*S1quench(Clusters[i].elec,.21)   ; 
      else cl_elec[i]  = 0 ; 
      cl_ndep[i]       = Clusters[i].nDeposits;  
      float cl_ene_reco = cl_nucl[i] + cl_elec[i] ; 
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
    
    ///////////////////////////////////////////////////////////////
    // bar multiplicity analysis 
    ///////////////////////////////////////////////////////////////
    
    vector<float> dep_step_cpy;
    double tmin =10000;
    double tmax =0;
                                                                                                      
    for(int i=0; i<theEvent.NDeposits;++i){
      //if(dep_mat[i]==ScintillatorIndex){
      if(dep_mat[i]==66 || dep_mat[i]==72 ){
        if (dep_step[i]==0 ) continue ; 
	dep_step_cpy.push_back(dep_step[i]);
	if(dep_time[i]<tmin){
	  tmin = dep_time[i];
	}
	if(dep_time[i]>tmax){
	  tmax = dep_time[i];
	}
     }
    }
    std::sort(dep_step_cpy.begin(), dep_step_cpy.end());
    //get the number of bars with at least one deposit
    std::vector<float>::iterator it;
    it = std::unique(dep_step_cpy.begin(),dep_step_cpy.end() );
    dep_step_cpy.resize(std::distance(dep_step_cpy.begin(),it) );
    nbars = dep_step_cpy.size();
       
    for(int i=0;i<nbars;i++){
      bar_step[i] = dep_step_cpy[i];
      tdep_min = tmin;
      tdep_max = tmax;
    }
    
    //intialize all the variables
    bar_tot_energy = 0 ; 
    bar_late_energy = 0 ; 
    for (int i=0;i<20;++i) { 
      bar_prompt_ene[i] = 0 ; 
      bar_late_ene[i] = 0 ;       
    } 
    
    
    //loop to fill bars 
    for(int i=0; i<theEvent.NDeposits;++i){
      //if(dep_mat[i]==ScintillatorIndex && dep_step[i] !=0 ){
      if( 1  ){
	quench  = 1 ; //  getqfactor_PS(dep_ene[i], dep_pdg[i]);
        bool isMatched=false;
	bar_tot_energy += ( dep_ene[i] * quench ) ;
	int j=0;
        while(!isMatched && j<nbars){
          if(bar_step[j] == dep_step[i]){
	    if (dep_time[i] - T0 <= 200)  bar_prompt_ene[j] += dep_ene[i] * quench;
	    if (dep_time[i] - T0 >  200)  bar_late_ene[j]   += dep_ene[i] * quench;
            isMatched=true;
          }
          j++;
        }
      }
    }
    
    //copy bar info in ad hoc vectors
    vector<float> bse_prompt, bse_late ; 
    vector<float> bse_prompt10, bse_late10 ; 
    vector<float> bse_prompt100, bse_late100 ; 
    for (int i=0;i<nbars;++i) {
      if (bar_prompt_ene[i] > 1 ) bse_prompt.push_back(bar_prompt_ene[i] ) ; 
      if (bar_late_ene[i]   > 1 ) bse_late.push_back(bar_late_ene[i] ) ; 
      if (bar_prompt_ene[i] > 10 ) bse_prompt10.push_back(bar_prompt_ene[i] ) ; 
      if (bar_late_ene[i]   > 10 ) bse_late10.push_back(bar_late_ene[i] ) ; 
      if (bar_prompt_ene[i] > 100 ) bse_prompt100.push_back(bar_prompt_ene[i] ) ; 
      if (bar_late_ene[i]   > 100 ) bse_late100.push_back(bar_late_ene[i] ) ; 
    }
    //cout << nbars << endl ; 
    std::sort(bse_late.begin(), bse_late.end())   ; 
    std::sort(bse_prompt.begin(), bse_prompt.end())   ;        

    //re-initialized arrays
    for (int i=0;i<20;++i) { 
      bar_prompt_ene[i] = 0 ; 
      bar_late_ene[i] = 0 ;       
    } 
    
    for (int i=0;i<bse_prompt.size();++i) { 
      //cout <<"P " << i <<"/" << bse_prompt.size()  <<" " << bse_prompt[i] << endl ; 
      bar_prompt_ene[i] = bse_prompt[bse_prompt.size()-1-i] ; 
    } 
    for (int i=0;i<bse_late.size();++i) { 
      bar_late_ene[i] = bse_late[bse_late.size()-1-i] ; 
      //cout <<"L " << i <<"/" << bse_late.size()  <<" " << bse_late[i] << endl ; 
      bar_late_energy += bar_late_ene[i] ; 
    } 
    
    nbars_prompt = bse_prompt.size()   ; 
    nbars_late = bse_late.size()  ; 
    nbars_prompt10 = bse_prompt10.size()   ; 
    nbars_late10   = bse_late10.size()  ; 
    nbars_prompt100= bse_prompt100.size()   ; 
    nbars_late100  = bse_late100.size()  ; 

    // Fill user variables
    for(int i=0;i<theEvent.NUsers;++i) {
      INT1[i]    = theUsers[i].UserInt1;
      INT2[i]    = theUsers[i].UserInt2;
      FLOAT1[i]  = theUsers[i].UserFloat1;
      FLOAT2[i]  = theUsers[i].UserFloat2; 
      DOUBLE[i]  = theUsers[i].UserDouble; 
    }
    // Fill pe variables    
    for (int i=0;i<38;++i) pe_ch[i] =0 ; 
    s1_max_frac =  -1. ; 
    
    for(int i=0;i<theEvent.NPE;++i) {
      pe_pmt[i]  = thePhotoElectrons[i].PMT;    
      pe_time[i] = thePhotoElectrons[i].Time; 
      if (theHeader.DetectorFlag < 3 ) ++pe_ch[pe_pmt[i]]  ; 
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
    
   //if (nclus==1) { 
   for (int i=0; i<Clusters.size() ; ++i) {
     double apothem = TPCEdge * 1.207 ; 
     double cosPi8  = 1./cos(PI/8.) ; 
     double HalfZ   = TPCHeight/2. ;

     isFV30[i]=0;
     isFV5[i]=0;
 
     if ( InsideOctagon((apothem - 5)*cosPi8   , cl_x[i], cl_y[i] ) && cl_z[i] > -HalfZ + 5  && cl_z[i]  < HalfZ - 5 ) isFV5[i]  =1 ; 
     if ( InsideOctagon((apothem - 10)*cosPi8  , cl_x[0], cl_y[0] ) && cl_z[0] > -HalfZ + 10 && cl_z[0]  < HalfZ - 10) isFV10 =1 ; 
     if ( InsideOctagon((apothem - 15)*cosPi8  , cl_x[0], cl_y[0] ) && cl_z[0] > -HalfZ + 15 && cl_z[0]  < HalfZ - 15) isFV15 =1 ; 
     if ( InsideOctagon((apothem - 20)*cosPi8  , cl_x[0], cl_y[0] ) && cl_z[0] > -HalfZ + 20 && cl_z[0]  < HalfZ - 20) isFV20 =1 ; 
     if ( InsideOctagon((apothem - 30)*cosPi8  , cl_x[i], cl_y[i] ) && cl_z[i] > -HalfZ + 30 && cl_z[i]  < HalfZ - 30) isFV30[i] =1 ; 
     if ( InsideOctagon((apothem - 35)*cosPi8  , cl_x[0], cl_y[0] ) && cl_z[0] > -HalfZ + 35 && cl_z[0]  < HalfZ - 35) isFV35 =1 ; 
   }

   // Fill the event
   
  if (FVcut == true) {
    if (nclusNR==1 && isFV30[goodNRidx]==1) dstree->Fill();
    }
  else if (NRcut == true) {
  if (nclusNR==1) dstree->Fill();
  }
  else dstree->Fill();

  }
  
  nentree->Fill() ; 

  // Write the tree
  dstree->Write();
  nentree->Write() ; 
  

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

