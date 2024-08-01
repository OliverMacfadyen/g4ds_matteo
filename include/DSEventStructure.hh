#ifndef DSEventStructure_h
#define DSEventStructure_h 1
#include <vector>

//using namespace std;

#pragma pack(1)
struct HeaderStructure {
  int Events;
  int Run;
  int PDG;
  int LArIndex;
  int ScintillatorIndex;
  float Rate;
  int DetectorFlag;
};

struct PhotoElectronStructure {
  int PMT;
  double Time;
};

struct PhotonStructure {
  int VolumeID;
  int PID;
  float Wavelength;
  float Position[3];
  double Time;
};

struct DepositStructure {
  int PID;
  int Volume;
  float Energy;
  float Step;
  float Position[3];
  double Time;
};

struct ClusterStructure {
  float  RecoilID;
  //int    Material;
  float  Energy;
  float  S1Energy;
  float  S2Energy;
  float  Position[3];
  double Time;
};

struct DaughterStructure {
  int Id;
  int PDG;
  int PID;
  int Process;
  double Time;
  float Energy;
  float Position[3];
  float Direction[3];
};

struct UserStructure {
  int UserInt1;
  int UserInt2;
  float UserFloat1;
  float UserFloat2;
  double UserDouble;
};

struct EventStructureDiskFormat {
  int EventID;            // Event ID
  int PDG;                // Particle Data Group Code
  double Time;            // Time
  float Energy;           // Kinetic energy
  float S1Energy;         // S1 visible kinetic energy
  float S2Energy;         // S2 visible kinetic energy
  float VetoVisEnergy;    // Neutron Veto visible kinetic energy
  float MuVisEnergy;      // Muon Veto visible kinetic energy
  float TPCDepEnergy;     // Energy deposited in the TPC (LAr only)
  float VetoDepEnergy;    // Energy deposited in scintillator
  float MuDepEnergy;      // Energy deposited in water
  float Position[3];      // Position
  float Direction[3];     // Direction
  float CenterOfMass[3];  // Center of mass
  int NPE;                // Number of photoelectrons in the TPC  (size of the vector
                          // thePhotonElectrons)
  int VetoNPE;            // Number of photoelectrons in the Veto (size of the vector
                          // theVetoPhotonElectrons)
  int MuNPE;              // Number of photoelectrons in the Tank (size of the vector
                          // theMuPhotonElectrons)
  int NPH;                // Number of photoelectrons in the TPC  (size of the vector
                          // thePhotons)
  int NDaughters;         // Number of daughters (size of the vector theDaughters)
  int NDeposits;          // Number of deposits  (size of the vector theDeposits)
  int NUsers;             // Size of the vector theUsers
  int NClusters;       // Size of the vector theClusters  
};

struct EventStructure : public EventStructureDiskFormat {
  std::vector<DaughterStructure> theDaughters;
  std::vector<DepositStructure> theDeposits;
  std::vector<ClusterStructure> theClusters;
  std::vector<PhotoElectronStructure> thePhotoElectrons;
  std::vector<PhotoElectronStructure> theVetoPhotoElectrons;
  std::vector<PhotoElectronStructure> theMuPhotoElectrons;
  std::vector<PhotonStructure> thePhotons;
  std::vector<UserStructure> theUsers;
};

#pragma pack()

#endif

/*
 * $Log: DSEventStructure.hh,v $
 * Revision 1.4  2014/11/13 16:47:09  dfranco
 * removed variables which were creating conflicts with the previous version of
 * g4ds10
 *
 * Revision 1.3  2014/10/13 18:43:46  swesterd
 * fixed veto PMT positions
 *
 * Revision 1.2  2014/05/08 10:59:21  pagnes
 * Scintillator Index added in binary header
 *
 * Revision 1.1  2014/05/07 12:20:52  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.6  2013/08/06 13:58:18  dfranco
 * Added 3 variables to store the deposited energy in LAr, scintillator, and
 * water. The last two are not yet implemented. g4rooter has been updated with 3
 * new variables: tpcene, vetoene, and muene
 *
 * Revision 1.5  2013/07/24 09:48:57  dfranco
 * Added S1 and S2 equivalent energies, the material index for each deposit, the
 * command killS1S2 to kill photons and electrons generated by DSLight (after
 * storing the equivalent energies)
 *
 * Revision 1.4  2013/04/04 09:04:17  dfranco
 * added step length info to the deposit structure
 *
 * Revision 1.3  2013/03/22 14:09:40  dfranco
 * Added the cvs logger code at the end of each file
 *
 */