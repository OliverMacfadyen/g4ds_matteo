
#include "G4Event.hh"
#include "G4PhysicalConstants.hh"
#include "G4RunManager.hh"
#include "G4String.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include "DSIO.hh"
#include "DSLogger.hh"
#include "DSStorage.hh"
#include "G4String.hh"

using namespace std;

DSIO* DSIO::me = 0;

// singleton
DSIO::DSIO() {

  fIsBinary = 0;
  fBinaryFileName = "outtest";
  fStreamDSGeometryFileName = "../data/detector/DSGeometry.dat";
  fStreamDSCryoProfFileName = "../data/detector/DSCryostatProfiles.dat";
  fStreamDSG3CryoProfFileName = "../data/detector/DSG3CryostatProfilesFromTechDrawings.dat";
  fStreamDS50VPMTGeometryFileName = "../data/detector/VPMTGeometry.dat";
  fStreamDS20kVPMTGeometryFileName = "../data/detector/DS20kPMTNeutronVeto.dat";
  fStreamDSOpticsFileName = "../data/detector/DSOptics.dat";
  fStreamLSVOpticsFileName = "../data/detector/LSVOptics.dat";
  fStreamDS50VetoPmtQEFileName = "../data/detector/DS50VetoPmtQE.dat";

  fStreamLArPropertiesFileName = "../data/physics/LArScintillatioProperties.txt";

  fStreamPCRindexFileName = "../data/physics/PCRefractionIndex.dat";
  fStreamPCAttLengthFileName = "../data/physics/PCAttenuationLength.dat";
  fStreamPPOAttLengthFileName = "../data/physics/PPOAttenuationLength.dat";
  fStreamTMBAttLengthFileName = "../data/physics/TMBAttenuationLength.dat";

  fStreamLumirrorReflectivityFileName = "../data/detector/LumirrorReflectivity.dat";

  // DS20k
  fStreamDS20kESRreflectivityFileName = "../data/detector/ESR_reflectivity_2_AstroCent.dat";
}

DSIO* DSIO::Get() {
  if (!me) me = new DSIO();

  return me;
}

const G4String& DSIO::CheckFileName(const G4String& mys) {
  std::ifstream file;
  std::ostringstream newfilename;
  newfilename << mys << ".log";
  file.open(newfilename.str().c_str(), std::ios::in);
  if (!file) {
    ChangeName(mys);
    newfilename.str("");
    newfilename << mys;
  } else {
    file.close();
    for (int i = 1; i < 1E4; i++) {
      newfilename.str("");
      newfilename << mys << "_v" << i << ".log";
      file.open(newfilename.str().c_str(), std::ios::in);
      if (!file) {
        newfilename.str("");
        newfilename << mys << "_v" << i;
        ChangeName(newfilename.str());
        DSLog(trace) << "Output files already exist! Name changed: " << newfilename.str() << endl;
        break;
      }
      file.close();
    }
  }
  fFileName = newfilename.str();
  return fFileName;
}

void DSIO::ChangeName(const G4String& mys) {
  std::ostringstream newfilename;

  newfilename << mys << ".log";
  SetLogFileName(newfilename.str());

  newfilename.str("");
  newfilename << mys << ".fil";

  SetBinaryFileName(newfilename.str());
}

//--------------------------------------------------------------
// Binary File
//--------------------------------------------------------------

void DSIO::OpenBinaryFile() {
  fBinaryFile.open(fBinaryFileName.c_str(), ofstream::out | ofstream::binary);
}

void DSIO::CloseBinaryFile() {
  // fBinaryFile.flush();
  fBinaryFile.close();
  DSLog(routine) << "Binary File Closed" << endlog;
}

//--------------------------------------------------------------
// Log Files
//--------------------------------------------------------------

void DSIO::OpenLogFiles() {
  if (!fStreamLogFile.is_open()) {
    DSLog(routine) << "Log file created " << endl;

    fStreamLogFile.open(fLogFileName, ofstream::out);
    fStreamLogFile << "#############  G4DS Log File  #############" << endl;
  }
}

void DSIO::CloseLogFiles() {
  fStreamLogFile.close();
  DSLog(routine) << "Log Files Flushed And Closed" << endlog;
}

//--------------------------------------------------------------
// G4DS Files
//--------------------------------------------------------------

void DSIO::OpenG4DSFile() {
  fG4DSFile.open(fG4DSFileName.c_str(), ifstream::in | ifstream::binary);
}

void DSIO::CloseG4DSFile() {
  fG4DSFile.close();
  DSLog(routine) << "G4DS File Closed" << endlog;
}

// DSGeometry
ifstream& DSIO::GetStreamDSGeometry() {
  if (!fStreamDSGeometry.is_open()) { fStreamDSGeometry.open(fStreamDSGeometryFileName, ifstream::in); }
  if (fStreamDSGeometry.eof()) {
    fStreamDSGeometry.close();
    fStreamDSGeometry.clear();
    fStreamDSGeometry.open(fStreamDSGeometryFileName, ifstream::in);
  }
  return fStreamDSGeometry;
}

void DSIO::CloseStreamDSGeometry() {

  if (fStreamDSGeometry.is_open()) fStreamDSGeometry.close();
  DSLog(routine) << "Geometry Input File Closed" << endlog;
}

// DSOptics
ifstream& DSIO::GetStreamDSOptics() {
  if (DSStorage::Get()->Get20KGeometry()) fStreamDSOpticsFileName = "../data/detector/DSOpticsDS20k.dat";

  if (!fStreamDSOptics.is_open()) { fStreamDSOptics.open(fStreamDSOpticsFileName, ifstream::in); }
  if (fStreamDSOptics.eof()) {
    fStreamDSOptics.close();
    fStreamDSOptics.clear();
    fStreamDSOptics.open(fStreamDSOpticsFileName, ifstream::in);
  }
  return fStreamDSOptics;
}

void DSIO::CloseStreamDSOptics() {
  if (fStreamDSOptics.is_open()) fStreamDSOptics.close();
  DSLog(routine) << "Optics Input File Closed" << endlog;
}

// LSVOptics
ifstream& DSIO::GetStreamLSVOptics() {
  if (!fStreamLSVOptics.is_open()) {
    fStreamLSVOptics.open(fStreamLSVOpticsFileName, ifstream::in);
    if (!fStreamLSVOptics.is_open()) DSLog(fatal) << "Error. File " << fStreamPCRindexFileName << " not found." << endlog;
  }
  if (fStreamLSVOptics.eof()) {
    fStreamLSVOptics.close();
    fStreamLSVOptics.clear();
    fStreamLSVOptics.open(fStreamLSVOpticsFileName, ifstream::in);
  }
  return fStreamLSVOptics;
}

void DSIO::CloseStreamLSVOptics() {
  if (fStreamLSVOptics.is_open()) fStreamLSVOptics.close();
  DSLog(routine) << "LSV Optics Input File Closed" << endlog;
}

// LArProperties
ifstream& DSIO::GetStreamLArProperties() {
  if (!fStreamLArProperties.is_open()) { fStreamLArProperties.open(fStreamLArPropertiesFileName, ifstream::in); }
  if (fStreamLArProperties.eof()) {
    fStreamLArProperties.close();
    fStreamLArProperties.clear();
    fStreamLArProperties.open(fStreamLArPropertiesFileName, ifstream::in);
  }
  return fStreamLArProperties;
}

void DSIO::CloseStreamLArProperties() {

  if (fStreamLArProperties.is_open()) fStreamLArProperties.close();
  DSLog(routine) << "Optics Input File Closed" << endlog;
}

ifstream& DSIO::GetStreamDSCryostatProfile() {
  if (!fStreamDSCryoProfile.is_open()) { fStreamDSCryoProfile.open(fStreamDSCryoProfFileName, ifstream::in); }
  if (fStreamDSCryoProfile.eof()) {
    fStreamDSCryoProfile.close();
    fStreamDSCryoProfile.clear();
    fStreamDSCryoProfile.open(fStreamDSCryoProfFileName, ifstream::in);
  }
  return fStreamDSCryoProfile;
}

// Veto PMT Geometry
ifstream& DSIO::GetStreamDSVPMTGeometry() {

  if (!fStreamDSVPMTGeometry.is_open()) {
    const G4bool is20k = DSStorage::Get()->GetIsDS20kLSV();
    // choose Geometry File
    // ** TEMPORARY -> should be independent of DSStorage
    if (is20k) {
      fStreamDSVPMTGeometryFileName = fStreamDS20kVPMTGeometryFileName;
    } else
      fStreamDSVPMTGeometryFileName = fStreamDS50VPMTGeometryFileName;
    fStreamDSVPMTGeometry.open(fStreamDSVPMTGeometryFileName, ifstream::in);
  }
  if (fStreamDSVPMTGeometry.eof()) {
    fStreamDSVPMTGeometry.close();
    fStreamDSVPMTGeometry.clear();
    fStreamDSVPMTGeometry.open(fStreamDSVPMTGeometryFileName, ifstream::in);
  }
  return fStreamDSVPMTGeometry;
}

void DSIO::CloseStreamDSVPMTGeometry() {

  if (fStreamDSVPMTGeometry.is_open()) fStreamDSVPMTGeometry.close();
  DSLog(routine) << "Veto PMT Geometry Input File Closed" << endlog;
}

// DS G3 Cryostat
ifstream& DSIO::GetStreamDSG3CryostatProfile() {

  if (!fStreamDSG3CryoProfile.is_open()) { fStreamDSG3CryoProfile.open(fStreamDSG3CryoProfFileName, ifstream::in); }
  if (fStreamDSG3CryoProfile.eof()) {
    fStreamDSG3CryoProfile.close();
    fStreamDSG3CryoProfile.clear();
    fStreamDSG3CryoProfile.open(fStreamDSG3CryoProfFileName, ifstream::in);
  }
  return fStreamDSG3CryoProfile;
}

void DSIO::CloseStreamDSCryostatProfile() {

  if (fStreamDSCryoProfile.is_open()) fStreamDSCryoProfile.close();
  DSLog(routine) << "Cryostat Profile Input File Closed" << endlog;
}

void DSIO::CloseStreamDSG3CryostatProfile() {

  if (fStreamDSG3CryoProfile.is_open()) fStreamDSG3CryoProfile.close();
  DSLog(routine) << "G3 Cryostat Profile Input File Closed" << endlog;
}

// LS: PC Refraction Index
ifstream& DSIO::GetStreamPCRefractionIndex() {
  if (!fStreamPCRindex.is_open()) {
    fStreamPCRindex.open(fStreamPCRindexFileName, ifstream::in);
    if (!fStreamPCRindex) DSLog(fatal) << "Error. File " << fStreamPCRindexFileName << " not found." << endlog;
  }
  if (fStreamPCRindex.eof()) {
    fStreamPCRindex.close();
    fStreamPCRindex.clear();
    fStreamPCRindex.open(fStreamPCRindexFileName, ifstream::in);
  }
  return fStreamPCRindex;
}

void DSIO::CloseStreamPCRefractionIndex() {
  if (fStreamPCRindex.is_open()) fStreamPCRindex.close();
  DSLog(routine) << "PC Refraction Index Input File Closed" << endlog;
}

// LS: PC Attenuation length
ifstream& DSIO::GetStreamPCAttLength() {
  if (!fStreamPCAttLength.is_open()) {
    fStreamPCAttLength.open(fStreamPCAttLengthFileName, ifstream::in);
    if (!fStreamPCAttLength) DSLog(fatal) << "Error. File " << fStreamPCAttLengthFileName << " not found." << endlog;
  }
  if (fStreamPCAttLength.eof()) {
    fStreamPCAttLength.close();
    fStreamPCAttLength.clear();
    fStreamPCAttLength.open(fStreamPCAttLengthFileName, ifstream::in);
  }
  return fStreamPCAttLength;
}

void DSIO::CloseStreamPCAttLength() {
  if (fStreamPCAttLength.is_open()) fStreamPCAttLength.close();
  DSLog(routine) << "PC Attenuation Length Input File Closed" << endlog;
}

// LS: PPO Attenuation length
ifstream& DSIO::GetStreamPPOAttLength() {
  if (!fStreamPPOAttLength.is_open()) {
    fStreamPPOAttLength.open(fStreamPPOAttLengthFileName, ifstream::in);
    if (!fStreamPPOAttLength) DSLog(fatal) << "Error. File " << fStreamPPOAttLengthFileName << " not found." << endlog;
  }
  if (fStreamPPOAttLength.eof()) {
    fStreamPPOAttLength.close();
    fStreamPPOAttLength.clear();
    fStreamPPOAttLength.open(fStreamPPOAttLengthFileName, ifstream::in);
  }
  return fStreamPPOAttLength;
}

void DSIO::CloseStreamPPOAttLength() {
  if (fStreamPPOAttLength.is_open()) fStreamPPOAttLength.close();
  DSLog(routine) << "PPO Attenuation Length Input File Closed" << endlog;
}

// LS: TMB Attenuation length
ifstream& DSIO::GetStreamTMBAttLength() {
  if (!fStreamTMBAttLength.is_open()) {
    fStreamTMBAttLength.open(fStreamTMBAttLengthFileName, ifstream::in);
    if (!fStreamTMBAttLength) DSLog(fatal) << "Error. File " << fStreamTMBAttLengthFileName << " not found." << endlog;
  }
  if (fStreamTMBAttLength.eof()) {
    fStreamTMBAttLength.close();
    fStreamTMBAttLength.clear();
    fStreamTMBAttLength.open(fStreamTMBAttLengthFileName, ifstream::in);
  }
  return fStreamTMBAttLength;
}

void DSIO::CloseStreamTMBAttLength() {
  if (fStreamTMBAttLength.is_open()) fStreamTMBAttLength.close();
  DSLog(routine) << "TMB Attenuation Length Input File Closed" << endlog;
}

// LSV: Lumirror Relfectivity
ifstream& DSIO::GetStreamLumirrorReflectivity() {
  if (!fStreamLumirrorReflectivity.is_open()) {
    fStreamLumirrorReflectivity.open(fStreamLumirrorReflectivityFileName, ifstream::in);
    if (!fStreamLumirrorReflectivity) DSLog(fatal) << "Error. File " << fStreamLumirrorReflectivityFileName << " not found." << endlog;
  }
  if (fStreamLumirrorReflectivity.eof()) {
    fStreamLumirrorReflectivity.close();
    fStreamLumirrorReflectivity.clear();
    fStreamLumirrorReflectivity.open(fStreamLumirrorReflectivityFileName, ifstream::in);
  }
  return fStreamLumirrorReflectivity;
}

void DSIO::CloseStreamLumirrorReflectivity() {
  if (fStreamLumirrorReflectivity.is_open()) fStreamLumirrorReflectivity.close();
  DSLog(routine) << "Lumirror Reflectivity Input File Closed" << endlog;
}

// LSV: PMT Quantum Effinciency
ifstream& DSIO::GetStreamDS50VetoPmtQE() {
  if (!fStreamDS50VetoPmtQE.is_open()) {
    fStreamDS50VetoPmtQE.open(fStreamDS50VetoPmtQEFileName, ifstream::in);
    if (!fStreamDS50VetoPmtQE) DSLog(fatal) << "Error. File " << fStreamDS50VetoPmtQEFileName << " not found." << endlog;
  }
  if (fStreamDS50VetoPmtQE.eof()) {
    fStreamDS50VetoPmtQE.close();
    fStreamDS50VetoPmtQE.clear();
    fStreamDS50VetoPmtQE.open(fStreamDS50VetoPmtQEFileName, ifstream::in);
  }
  return fStreamDS50VetoPmtQE;
}

void DSIO::CloseStreamDS50VetoPmtQE() {
  if (fStreamDS50VetoPmtQE.is_open()) fStreamDS50VetoPmtQE.close();
  DSLog(routine) << "DS50 Veto PMT Quantum Efficiency Input File Closed" << endlog;
}

// DS20k - ESR reflectivity for the LAr->ESR surface in the DS20k veto. Measured
// at AstroCent
ifstream& DSIO::GetStreamDS20kESRreflectivity() {
  if (!fStreamDS20kESRreflectivity.is_open()) {
    fStreamDS20kESRreflectivity.open(fStreamDS20kESRreflectivityFileName, ifstream::in);
    if (!fStreamDS20kESRreflectivity) DSLog(fatal) << "Error. File " << fStreamDS20kESRreflectivityFileName << " not found." << endlog;
  }
  if (fStreamDS20kESRreflectivity.eof()) {
    fStreamDS20kESRreflectivity.close();
    fStreamDS20kESRreflectivity.clear();
    fStreamDS20kESRreflectivity.open(fStreamDS20kESRreflectivityFileName, ifstream::in);
  }
  return fStreamDS20kESRreflectivity;
}

/*
 * $Log: DSIO.cc,v $
 * Revision 1.4  2015/06/02 16:25:58  pagnes
 * added command to use a different LArScintillationProperties.txt file
 * (/ds/manager/lar_properties)
 *
 * Revision 1.3  2014/10/13 18:43:57  swesterd
 * fixed veto PMT positions
 *
 * Revision 1.2  2014/07/25 14:07:20  perassos
 * Implementation of the DSG3 (TPC + NV). Configuration #8
 *
 * Revision 1.1  2014/05/07 12:21:03  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.5  2014/03/11 16:49:56  dfranco
 * read optical parameters from external file (DSOptics.dat)
 *
 * Revision 1.4  2014/03/11 09:54:38  meregaglia
 * Added generator starting from energy deposits
 *
 * Revision 1.3  2013/03/22 14:09:40  dfranco
 * Added the cvs logger code at the end of each file
 *
 */
