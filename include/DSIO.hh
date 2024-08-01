#ifndef _BXIO_HH
#define _BXIO_HH 1

#include <fstream>
#include <iostream>
#include "G4String.hh"

//using namespace std;

class DSIO {
 private:
  DSIO();

 public:
  static DSIO* Get();

  virtual ~DSIO() {}

  const G4String& CheckFileName(const G4String&);
  const G4String& GetFileName() const { return fFileName; }
  void SetFileName(const G4String& a) { fFileName = a; }

  // Binary File
  std::ofstream& GetBinaryFile() { return fBinaryFile; }
  const G4String& GetBinaryFileName() const { return fBinaryFileName; }
  void OpenBinaryFile();
  void CloseBinaryFile();
  G4bool GetIsBinary() const { return fIsBinary; }
  void SetIsBinary(G4bool a) { fIsBinary = a; }

  // Log Files
  const G4String& GetLogFileName() const { return fLogFileName; }
  std::ofstream& GetStreamLogFile() { return fStreamLogFile; }
  void OpenLogFiles();
  void CloseLogFiles();

  // G4DS
  std::ifstream& GetG4DSFile() { return fG4DSFile; }
  const G4String& GetG4DSFileName() const { return fG4DSFileName; }
  void OpenG4DSFile();
  void CloseG4DSFile();
  void SetIsG4DS(G4bool a) { fIsG4DS = a; }
  G4bool IsG4DS() const { return fIsG4DS; }

  void SetDSGeometry(const G4String& val) { fStreamDSGeometryFileName = val; }
  const G4String& GetDSGeometry() const { return fStreamDSGeometryFileName; }
  std::ifstream& GetStreamDSGeometry();
  std::ifstream& GetStreamDSCryostatProfile();
  std::ifstream& GetStreamDSG3CryostatProfile();

  void CloseStreamDSGeometry();
  void CloseStreamDSCryostatProfile();
  void CloseStreamDSG3CryostatProfile();

  void SetDSVPMTGeometry(const G4String& val) { fStreamDSVPMTGeometryFileName = val; }
  const G4String& GetDSVPMTGeometry() const { return fStreamDSVPMTGeometryFileName; }
  std::ifstream& GetStreamDSVPMTGeometry();
  void CloseStreamDSVPMTGeometry();

  std::ifstream& GetStreamDSOptics();
  void SetDSOpticsFileName(const G4String& val) { fStreamDSOpticsFileName = val; }
  void CloseStreamDSOptics();

  std::ifstream& GetStreamLSVOptics();
  void SetLSVOpticsFileName(const G4String& val) { fStreamLSVOpticsFileName = val; }
  void CloseStreamLSVOptics();

  std::ifstream& GetStreamLArProperties();
  void SetLArPropertiesFileName(const G4String& val) { fStreamLArPropertiesFileName = val; }
  void CloseStreamLArProperties();

  // PC refraction index
  std::ifstream& GetStreamPCRefractionIndex();
  void SetPCRefractionIndexFileName(const G4String& val) { fStreamPCRindexFileName = val; }
  void CloseStreamPCRefractionIndex();

  // PC attenuation length
  std::ifstream& GetStreamPCAttLength();
  void SetPCAttLengthFileName(const G4String& val) { fStreamPCAttLengthFileName = val; }
  void CloseStreamPCAttLength();

  // PPO attenuation length
  std::ifstream& GetStreamPPOAttLength();
  void SetPPOAttLengthFileName(const G4String& val) { fStreamPPOAttLengthFileName = val; }
  void CloseStreamPPOAttLength();

  // TMB attenuation length
  std::ifstream& GetStreamTMBAttLength();
  void SetTMBAttLengthFileName(const G4String& val) { fStreamTMBAttLengthFileName = val; }
  void CloseStreamTMBAttLength();

  // Lumirror Reflectivity
  std::ifstream& GetStreamLumirrorReflectivity();
  void SetLumirrorReflectivityFileName(const G4String& val) { fStreamLumirrorReflectivityFileName = val; }
  void CloseStreamLumirrorReflectivity();

  // DS50 Veto PMT QE
  std::ifstream& GetStreamDS50VetoPmtQE();
  void SetDS50VetoPmtQEFileName(const G4String& val) { fStreamDS50VetoPmtQEFileName = val; }
  void CloseStreamDS50VetoPmtQE();

  // DS20k ESR reflectivity measured at Astrocent for neutron veto
  std::ifstream& GetStreamDS20kESRreflectivity();
  // void       SetDS20kESRreflectivityFileName(const string& val) {
  // fStreamDS20kESRreflectivityFileName = val ;} void
  // CloseStreamDS20kESRreflectivity();

  void SetG4DSFile(const G4String& name) { fG4DSFileName = name; }

 private:
  static DSIO* me;

  std::ifstream fStreamDSGeometry;
  G4String fStreamDSGeometryFileName;

  // optics tuning
  std::ifstream fStreamDSOptics;
  G4String fStreamDSOpticsFileName;

  // LSV optics tuning
  std::ifstream fStreamLSVOptics;
  G4String fStreamLSVOpticsFileName;

  std::ifstream fStreamLArProperties;
  G4String fStreamLArPropertiesFileName;

  std::ifstream fStreamDSCryoProfile;
  G4String fStreamDSCryoProfFileName;

  std::ifstream fStreamDSG3CryoProfile;
  G4String fStreamDSG3CryoProfFileName;

  // Veto PMTs
  std::ifstream fStreamDSVPMTGeometry;
  G4String fStreamDSVPMTGeometryFileName;
  G4String fStreamDS50VPMTGeometryFileName;
  G4String fStreamDS20kVPMTGeometryFileName;

  // PseudoCumene (PC) refraction index
  std::ifstream fStreamPCRindex;
  G4String fStreamPCRindexFileName;

  // PseudoCumene (PC) attenuation length
  std::ifstream fStreamPCAttLength;
  G4String fStreamPCAttLengthFileName;

  // PPO attenuation length
  std::ifstream fStreamPPOAttLength;
  G4String fStreamPPOAttLengthFileName;

  // TMB attenuation length
  std::ifstream fStreamTMBAttLength;
  G4String fStreamTMBAttLengthFileName;

  // Lumirror reflectivity
  std::ifstream fStreamLumirrorReflectivity;
  G4String fStreamLumirrorReflectivityFileName;

  // DS50 Veto PMT Quantum Efficiency
  std::ifstream fStreamDS50VetoPmtQE;
  G4String fStreamDS50VetoPmtQEFileName;

  std::ifstream fStreamDS20kESRreflectivity;
  G4String fStreamDS20kESRreflectivityFileName;

  std::ofstream fBinaryFile;
  std::ofstream fStreamLogFile;

  G4String fFileName;
  G4String fLogFileName;
  G4String fBinaryFileName;
  G4String fG4DSFileName;
  std::ifstream fG4DSFile;
  G4bool fIsBinary;
  G4bool fIsG4DS;

  void SetLogFileName(const G4String& name) { fLogFileName = name; }
  void SetBinaryFileName(const G4String& name) { fBinaryFileName = name; }

  void ChangeName(const G4String&);
};

#endif
/*
 * $Log: DSIO.hh,v $
 * Revision 1.4  2015/06/02 16:25:56  pagnes
 * added command to use a different LArScintillationProperties.txt file
 * (/ds/manager/lar_properties)
 *
 * Revision 1.3  2014/10/13 18:43:46  swesterd
 * fixed veto PMT positions
 *
 * Revision 1.2  2014/07/25 14:07:25  perassos
 * Implementation of the DSG3 (TPC + NV). Configuration #8
 *
 * Revision 1.1  2014/05/07 12:20:53  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.5  2014/03/19 16:44:43  dfranco
 * update external configuration reader for optics tuning
 *
 * Revision 1.4  2014/03/11 16:50:00  dfranco
 * read optical parameters from external file (DSOptics.dat)
 *
 * Revision 1.3  2014/03/11 09:56:25  meregaglia
 * Added generator starting from energy deposits
 *
 * Revision 1.2  2013/03/22 14:09:40  dfranco
 * Added the cvs logger code at the end of each file
 *
 */
