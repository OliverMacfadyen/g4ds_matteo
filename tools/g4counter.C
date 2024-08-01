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


int _NEVENTS ;
int _COUNTER         = 0 ;
int _GLOBAL_COUNTER  = 0 ;
bool _NPE            = 0 ;
bool _S1ENE          = 0 ;
ofstream lout;

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
  _NEVENTS =   theHeader.Events ;
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
  
  
  int size_structs =   theEvent.NDaughters*sizeof(DaughterStructure) 
                     + theEvent.NDeposits *sizeof(DepositStructure)
                     + theEvent.NUsers *sizeof(UserStructure)
                     + theEvent.NPH *sizeof(PhotonStructure)
                     + theEvent.NPE *sizeof(PhotoElectronStructure)
                     + theEvent.VetoNPE *sizeof(PhotoElectronStructure)
                     + theEvent.MuNPE *sizeof(PhotoElectronStructure);
  
  file->ignore(size_structs);
  //file->seekg(size_structs+file->tellg());
  //for(int i=0; i<theEvent.NDaughters; i++)  _readDaughter(file);
  //for(int i=0; i<theEvent.NDeposits; i++)  _readDeposit(file);
  //for(int i=0; i<theEvent.NUsers; i++)     _readUser(file);
  //for(int i=0; i<theEvent.NPH; i++)        _readPhoton(file);
  //for(int i=0; i<theEvent.NPE; i++)        _readPhotoElectron(file);
  //for(int i=0; i<theEvent.VetoNPE; i++)    _readPhotoElectron(file);
  //for(int i=0; i<theEvent.MuNPE; i++)      _readPhotoElectron(file);
  
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

//--------------------------------------------
//                    exe 
//--------------------------------------------
void exe (string filename, int first) {
  ifstream *_bin_fstream = new ifstream (filename.c_str(), std::ios::binary);
  _readHeader(_bin_fstream); 


  for(int i=0; i<_NEVENTS; i++) { 
    if(!_readEvent(_bin_fstream))          break ; 
    if(!(i % 100000)) cout << i << " " << float(i)/_NEVENTS << endl ;
    if(_bin_fstream->eof())                break ;
    ++_GLOBAL_COUNTER ;  
    _COUNTER++;
  } 

  cout << _GLOBAL_COUNTER<<"=" <<float(_GLOBAL_COUNTER)/_NEVENTS << endl ;
}


//--------------------------------------------
//                    main 
//--------------------------------------------
int main (int argc, char *argv[]) {

  if(argc == 1 || (argc > 1 && !string(argv[1]).find("help")) )  { 
  //if(file.find("help") < 10) {
    cout << "Usage: g4counter <filename>" <<endl ;
    cout <<endl ;
    return 0 ;
  }

  std::string delimiter = "/";
  string filename = argv[1];
  
  string s = filename;

  size_t pos = 0;
  std::string token;
  while ((pos = s.find(delimiter)) != std::string::npos) {
      token = s.substr(0, pos);
      s.erase(0, pos + delimiter.length());
  }
  string log = s; 
  for(int i=0;i<4;++i) log.erase( log.end()-1 );
  log += ".counter";
  lout.open(log,ofstream::out);
  if(filename.find(".fil") != string::npos) {
    cout << "filename: " << filename << endl ;
    exe(filename, 1);
    lout << filename << endl ;
  }
  lout <<  _GLOBAL_COUNTER << " " << _NEVENTS << endl ;
  lout <<  float(_GLOBAL_COUNTER)/_NEVENTS << endl ;
  lout.close();
  
  
  cout << "number of written events: " << _COUNTER << endl ;
  return 0;
  
}

