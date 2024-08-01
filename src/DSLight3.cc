////////////////////////////////////////////////////////////////////////
//
// Author: D. Franco (dfranco@in2p3.fr)
//
// S1 S2  Scintillation Light Class Implementation
////////////////////////////////////////////////////////////////////////

#include "DSLight3.hh"
#include "DSEventHandler.hh"
#include "DSIO.hh"
#include "DSLogger.hh"
#include "DSMaterial.hh"
#include "DSParameters.hh"
#include "DSStorage.hh"
#include "G4EmProcessSubType.hh"
#include "G4EventManager.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleTypes.hh"
#include "G4PhysicalConstants.hh"
#include "G4SteppingManager.hh"
#include "G4TrackVector.hh"

using namespace std;

/*
* some important parameters:

* LAr scintillation properties are read from a specific file (see ReadData() )
* The drift time calculation is based on the fLArGArBoundaryZ value (set in
DSDetectorConstruction.cc )
*/



///////////////////////////////////////////////////////////////////
DSLight3::DSLight3(const G4String& processName, G4ProcessType type) : G4VRestDiscreteProcess(processName, type) {

  SetProcessType(fUserDefined);

  fTrackSecondariesFirst = false;

  // isFirstEnDepInLAr = true;

  fRecoProb = 0;
  fLArGArBoundaryZ     = DSStorage::Get()->GetLArGArBoundaryPosZ();
  fGridSurfaceDistance = DSStorage::Get()->GetGridSurfaceDistance();

  //if (DSStorage::Get()->GetTunedS1At200V()) DriftField = 200.0 * volt / cm;

  if (DSStorage::Get()->GetKillS1S2()) {
    DSStorage::Get()->SetKillS1(true);
    DSStorage::Get()->SetKillS2(true);
  }

  
  fGASGAP = DSParameters::Get()->GetGasPocketThickness();
  if (DSStorage::Get()->Get20KGeometry()) fGASGAP = DSStorage::Get()->GetDS20kGasPocketThickness();

  fYieldScalingFactor = 1 ;
  if (DSStorage::Get()->GetFastSimulation()) {
    fYieldScalingFactor = DSParameters::Get()->GetRealPmtMaxQe();
    if (DSStorage::Get()->Get20KGeometry())  fYieldScalingFactor = DSParameters::Get()->GetSiPMMaxPDE();
  }

  fArgonEnergy.clear();
  fArgonELoss.clear();

  fElectronEnergy.clear();
  fElectronELoss.clear();

  fAlphaEnergy.clear();
  fAlphaELoss.clear();

  fS2CorrectionPosition.clear();
  fS2CorrectionEfficiency.clear();

  recombination_ARIS_200_energies.clear();
  recombination_ARIS_200_values.clear();

  LEff_ARIS_energies.clear();
  LEff_ARIS_values.clear();

  ReadData();
  IsMeiQuenching = false;
  IsLindhardQuenching = false;
  IsARISQuenching = true;
  DSLog(debugging) << GetProcessName() << " is created " << endlog;

}

DSLight3::~DSLight3() {}

G4VParticleChange* DSLight3::AtRestDoIt(const G4Track& aTrack, const G4Step& aStep) {
  return DSLight3::PostStepDoIt(aTrack, aStep);
}

///////////////////////////////////////////////////////////////////
G4VParticleChange* DSLight3::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep) {

  aParticleChange.Initialize(aTrack);
  // cout <<
  // G4EventManager::GetEventManager()->GetStackManager()->GetNTotalTrack () <<
  // endl ;

  // At the event begin
  // if(aTrack.GetTrackID() == 1 && aTrack.GetCurrentStepNumber () == 1)  {
  //  isFirstEnDepInLAr = true;
  //  DSStorage::Get()->SetDSLightTrackSecondaries(1);
  //}

  const G4ParticleDefinition* myParticle = aTrack.GetParticleDefinition();
  const G4Material* myMaterial = aStep.GetPreStepPoint()->GetMaterial();
  
  // continue if optical photon
  if (myParticle->GetPDGEncoding() == -22) return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);
  // continue if not an argon
  if ((G4int)myMaterial->GetElement(0)->GetZ() != 18) return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);
  
  // continue if zero deposit
  if (aStep.GetTotalEnergyDeposit() <= 0) return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);

  // continue if not in LAr TPC or above grid
  int mat_idx = (int) aTrack.GetVolume()->GetLogicalVolume()->GetMaterial()->GetIndex();
  if (mat_idx != DSMaterial::Get()->GetTPCMaterialIndex() && mat_idx != DSStorage::Get()->GetLArAboveGridIndex() && mat_idx != DSMaterial::Get()->GetIVMaterialIndex()) 
    return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);
  


  // DF: remove following commented code (old) if above is validated 


  // // if optical photon...
  // if (myParticle->GetPDGEncoding() == -22 ||
  //     // ... or if not in LAr...
  //     // DSStorage::Get()->GetLiquidArgonIndex() != (G4int)
  //     // myMaterial->GetIndex()  ||
  //     18 != (G4int)myMaterial->GetElement(0)->GetZ() ||  // generalize to other LAr volumes, but negelect
  //                                                        // NSLiquidArgon (by def) and OVLiquidArgon (outside)
  //                                                        // for the moment
  //     myMaterial->GetName() == "NSLiquidArgon" || myMaterial->GetName() == "OVLiquidArgon" || myMaterial->GetName() == "GaseousArgon" ||
  //     // ... or if deposited energy is null...
  //     aStep.GetTotalEnergyDeposit() <= 0 )
  //     // ... or if the charge is zero
  //     //(myParticle->GetPDGCharge() == 0 && myParticle->GetPDGEncoding() != 22 && myParticle->GetPDGEncoding() != 2112)
  //   // ... return unchanged track.
  //   return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);

  // G4double      myPDG    = myParticle->GetPDGEncoding() ;
  G4double myZ       = myParticle->GetAtomicNumber();
  G4double myT0      = aStep.GetPreStepPoint()->GetGlobalTime();
  G4ThreeVector myX0 = aStep.GetPreStepPoint()->GetPosition();
  G4double myT1      = aStep.GetPostStepPoint()->GetGlobalTime();
  G4ThreeVector myX1 = aStep.GetPostStepPoint()->GetPosition();
  G4double myDepEne  = aStep.GetTotalEnergyDeposit();
  // G4double      myKinEne = aTrack.GetKineticEnergy () +
  // aStep.GetTotalEnergyDeposit (); G4double myKinEne =
  // aStep.GetPreStepPoint()->GetKineticEnergy();

  // G4double myFanoFactor = 0.2;
  G4double myExcitationRatio = myZ == 18 ?  fExcitationRatioNR : fExcitationRatioER;
  if (myParticle->GetPDGEncoding() == 22 )  myExcitationRatio = fExcitationRatioER;
  if (myParticle->GetPDGEncoding() == 2112) myExcitationRatio = fExcitationRatioNR;

  //////////////////////////
  //                      //
  //   Quenching Factor   //
  //                      //
  //////////////////////////

  fQuenchingFactor = (myZ == 18 || myParticle->GetPDGEncoding() == 2112) ? GetLArNuclearQuenching(aTrack.GetVertexKineticEnergy()) : 1;

  ///////////////////////////////////
  //                               //
  //   Recombination Probability   //
  //                               //
  ///////////////////////////////////

  G4double InitialKinEne = fQuenchingFactor * aTrack.GetVertexKineticEnergy();

  G4MaterialPropertiesTable* aMaterialPropertiesTable = myMaterial->GetMaterialPropertiesTable();
  G4double DriftField = aMaterialPropertiesTable->GetConstProperty("ELECTRICFIELD");

  if (DriftField == 0) fRecoProb = 1;
  else
    fRecoProb = GetRecoProbAt200Vcm(InitialKinEne, bool(myZ == 18))*GetRecombinationWrt200Vcm(DriftField);
  
  G4double myNumQuanta = fQuenchingFactor * myDepEne / fMeanQuantaEnergy;
  G4double myNumIons = myNumQuanta / (1 + myExcitationRatio);
  G4double myNumExcitons = myNumQuanta - myNumIons;
  // cout << myDepEne / keV <<" " << InitialKinEne/keV << endl ;

  G4double myphotons = myNumExcitons + myNumIons * fRecoProb;
  G4double myelectrons = myNumIons * (1 - fRecoProb);

  //compute the true prompt fraction for PSD
  G4double mySingletTripletRatio = GetPromptFraction(myZ, InitialKinEne);


  bool isInActiveLAr = (mat_idx == DSMaterial::Get()->GetTPCMaterialIndex()) ||  (mat_idx == DSStorage::Get()->GetLArAboveGridIndex());
  //bool isInActiveLAr = (aStep.GetPreStepPoint()->GetPhysicalVolume()->GetName() == "ActiveLAr") || (aStep.GetPreStepPoint()->GetPhysicalVolume()->GetName() == "LArAboveGrid") 

  // ? 1 : 0;

  // Add clustering in Active LAr -- active by default
  if (DSStorage::Get()->GetClusterWritingActivated() &&  isInActiveLAr==1 ) {
    CreateClusters(myphotons, myelectrons, myT0,  myX1 , myDepEne , mySingletTripletRatio); //,  int (aStep.GetPreStepPoint()->GetMaterial()->GetIndex()) ); // myX0 + (myX1-myX0)/2. ) ;
    //cout << LocalClus.X0 << " " << isInActiveLAr <<" " <<  myphotons<< " "<< myT0<< " "<< myX0<< " "<< myX1<< " "<<  myX0 + (myX1-myX0)/2. <<endl ;
  }


  // Kill S1 and S2 if required   (see the constructor - we want to save
  // electrons!)
  // if( DSStorage::Get()->GetKillS1S2() )
  //  return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);


  G4int myNumPhotons = 0;
  G4int myNumElectrons = 0;

  if (myphotons < 20) myNumPhotons = int(G4Poisson(myphotons) + 0.5);
  else
    myNumPhotons = int(G4RandGauss::shoot(myphotons, sqrt(myphotons)) + 0.5);
  if (myNumPhotons < 0) myNumPhotons = 0;

  if (myelectrons < 20) myNumElectrons = int(G4Poisson(myelectrons) + 0.5);
  else
    myNumElectrons = int(G4RandGauss::shoot(myelectrons, sqrt(myelectrons)) + 0.5);
  if (myNumElectrons < 0) myNumElectrons = 0;

  // G4int myNumPhotons   = int ( myNumExcitons + myNumIons*fRecoProb + 0.5);
  // G4int myNumElectrons = int ( myNumQuanta + 0.5 ) - myNumPhotons;        //
  // + 0.5 for the rounding
  // G4cout << "Debug from Ligth3 " <<
  // aStep.GetPreStepPoint()->GetPhysicalVolume()->GetName() << G4endl;
  if (isInActiveLAr) {
    DSEventHandler::Get()->SetS1Energy(DSEventHandler::Get()->GetS1Energy() + (float(myNumPhotons) * fMeanQuantaEnergy) / keV);
    DSEventHandler::Get()->SetS2Energy(DSEventHandler::Get()->GetS2Energy() + (float(myNumElectrons) * fMeanQuantaEnergy) / keV);
  }

  // scale by a calibration factor set in LArScintillationProperties ;
  myNumPhotons *= fS1ScaleFactor;

  //scale by a factor equal to MaxQE or MaxPDE for fast simulations
  myNumPhotons *= fYieldScalingFactor;
  

  aParticleChange.SetNumberOfSecondaries(myNumPhotons + 1000 * myNumElectrons);
  if (fTrackSecondariesFirst) {
    if (aTrack.GetTrackStatus() == fAlive) aParticleChange.ProposeTrackStatus(fSuspend);
  }

  // if KillS1 or KillS1S2, No S1 photons
  if (isInActiveLAr) {
    if (DSStorage::Get()->GetKillS1()) myNumPhotons = 0;
  } else
    myNumPhotons *= DSStorage::Get()->GetVetoYieldFactor();
  // Generate S1

  for (int j = 0; j < (int)myNumPhotons; j++) {

    // Momentum
    G4double cost = 1 - 2 * G4UniformRand();
    G4double sint = sqrt(1 - cost * cost);

    G4double phi = 2 * M_PI * G4UniformRand();
    G4double cosp = cos(phi);
    G4double sinp = sin(phi);

    G4double px = sint * cosp;
    G4double py = sint * sinp;
    G4double pz = cost;

    G4ThreeVector myPhotonMomentum(px, py, pz);

    // Polarization
    G4double sx = cost * cosp;
    G4double sy = cost * sinp;
    G4double sz = -sint;

    G4ThreeVector myPhotonPolarization(sx, sy, sz);
    G4ThreeVector myPerpendicular = myPhotonMomentum.cross(myPhotonPolarization);

    phi = 2 * M_PI * G4UniformRand();
    cosp = cos(phi);
    sinp = sin(phi);

    myPhotonPolarization = cosp * myPhotonPolarization + sinp * myPhotonMomentum;
    myPhotonPolarization = myPhotonPolarization.unit();

    // Energy
    G4double mySampledEnergy = 0;
    while (mySampledEnergy <= 0) mySampledEnergy = G4RandGauss::shoot(fPhotonEneMean, fPhotonEneWidth);

    // Position
    G4ThreeVector myPhotonPosition = myX0 + G4UniformRand() * (myX1 - myX0);

    // Time
    //   photonTime = Global Time ( + Recombination Time ) + Singlet/Triplet
    //   Time
    G4double myPhotonTime = myT0 + G4UniformRand() * (myT1 - myT0);

    // G4double mySingletTripletRatioEx   = 0.;        // Singlet to triplet
    // ratio for excitons G4double mySingletTripletRatioReco = 0.;        //
    // Singlet to triplet ratio for recombining electrons
    G4UniformRand() < mySingletTripletRatio ? myPhotonTime -= fLArTauFast * log(G4UniformRand()) : myPhotonTime -= fLArTauSlow * log(G4UniformRand());

    // Create a new Photon
    G4DynamicParticle* myScintillationPhoton = new G4DynamicParticle(G4OpticalPhoton::OpticalPhoton(), myPhotonMomentum);
    myScintillationPhoton->SetPolarization(myPhotonPolarization.x(), myPhotonPolarization.y(), myPhotonPolarization.z());
    myScintillationPhoton->SetKineticEnergy(mySampledEnergy);

    // Create a new track
    G4Track* aSecondaryTrack = new G4Track(myScintillationPhoton, myPhotonTime, myPhotonPosition);
    aParticleChange.AddSecondary(aSecondaryTrack);
  }

  // Generate S2
  if (isInActiveLAr && DriftField > 0 && fLArGArBoundaryZ != -100 * m) {

    bool isAboveGrid = (mat_idx == DSStorage::Get()->GetLArAboveGridIndex());

    for (int j = 0; j < myNumElectrons; j++) {
      // Position
      G4ThreeVector aSecondaryPosition = myX0 + G4UniformRand() * (myX1 - myX0);

      G4double myEnDepZ = aSecondaryPosition.z();

      G4double myVDrift = GetLArDriftVelocity(myMaterial->GetTemperature(), DriftField);

      G4double mySigmaDT = sqrt(fSigma0_Sq + 2 * fD_T * fabs(fLArGArBoundaryZ - myEnDepZ) / myVDrift);  // cm
      G4double dr = fabs(G4RandGauss::shoot(0., mySigmaDT));
      G4double myPhi = 2 * M_PI * G4UniformRand();

      aSecondaryPosition[0] += cos(myPhi) * dr;
      aSecondaryPosition[1] += sin(myPhi) * dr;
      G4ThreeVector aSecondaryPosition_old = aSecondaryPosition;
      aSecondaryPosition[2] = fLArGArBoundaryZ + 1 * um;
      // Time
      if (DriftField < 100 * volt / cm && myMaterial->GetState() == kStateLiquid) fD_L = 8 * cm2 / s;

      G4double mySigmaDL = sqrt(fSigma0_Sq + 2 * fD_L * fabs(fLArGArBoundaryZ - myEnDepZ) / myVDrift);  // cm
      G4double dt = G4RandGauss::shoot(0., (mySigmaDL / myVDrift));
      G4double myTDrift = 0;
      if (isAboveGrid) {
        myTDrift=  (fLArGArBoundaryZ + fGridSurfaceDistance - myEnDepZ) / myVDrift;   
      } else {
        myTDrift = (fLArGArBoundaryZ - myEnDepZ) / myVDrift + fGridSurfaceDistance / GetLArDriftVelocity(myMaterial->GetTemperature(), DSStorage::Get()->GetExtractionField()); 
      }

      G4double aSecondaryTime = myT0 + G4UniformRand() * (myT1 - myT0) + myTDrift + dt;

      // the aSecondaryPosition  spread may cause some electrons
      // to be generated inside the gas pocket (for boundary evts).
      // we will these electrons
      if (aSecondaryTime < 0.0) continue;

      // Kill Thermal Electrons according to purity...
      G4double prob = exp(-(myTDrift / ns) / (fTaueLAr / ns));
      G4double myTestVal = G4UniformRand();
      if (myTestVal > prob) continue;

      else {
        if (DSStorage::Get()->GetKillS2() == false) { CreateS2(aSecondaryPosition, aSecondaryTime, &aParticleChange); }
        if (DSStorage::Get()->GetWriteThermalElectrons()) {
          DSEventHandler::Get()->SetDepPID(-1);
          DSEventHandler::Get()->SetDepVolume(0);
          DSEventHandler::Get()->SetDepEnergy(0);
          DSEventHandler::Get()->SetDepStep(0);
          DSEventHandler::Get()->SetDepTime(aSecondaryTime / ns);
          DSEventHandler::Get()->SetDepPosition(aSecondaryPosition_old / cm);
          DSEventHandler::Get()->SetDeposits();
        }
      }
    }
  }
  // ** --------------------

  return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);
}

//////////////////////////////////////////////
//  GetMeanFreePath
//////////////////////////////////////////////
G4double DSLight3::GetMeanFreePath(const G4Track&, G4double, G4ForceCondition* condition) {
  *condition = StronglyForced;

  return DBL_MAX;
}

G4double DSLight3::GetMeanLifeTime(const G4Track&, G4ForceCondition* condition) {
  *condition = Forced;
  return DBL_MAX;
}

/////////////////////////////////////////////
//        BinomFluct
/////////////////////////////////////////////
G4int DSLight3::BinomFluct(G4int N0, G4double prob) {
  G4int N1 = 0;
  if (prob == 0.00) return N1;
  if (prob == 1.00) return N0;

  if (N0 < 10) {
    for (G4int i = 0; i < N0; i++) {
      if (G4UniformRand() < prob) N1++;
    }
  } else {
    G4double mean = N0 * prob;
    G4double sigma = sqrt(N0 * prob * (1 - prob));
    N1 = G4int(G4RandGauss::shoot(mean, sigma) + 0.5);
    if (N1 > N0) N1 = N0;
    if (N1 < 0) N1 = 0;
  }
  return N1;
}

//////////////////////////////////////////////
//  Read Energy Loss Data (form SRIM/TRIM)
//////////////////////////////////////////////
void DSLight3::ReadData() {
  float ene, pos, val;
  ifstream ff;

  ff.open("../data/physics/electron_eloss.dat");
  while (!ff.eof()) {
    ff >> ene >> val;
    if (ff.eof()) break;
    fElectronEnergy.push_back(ene * keV);
    fElectronELoss.push_back(val * 10);
  }
  ff.close();

  ff.open("../data/physics/argon_eloss.dat");
  while (!ff.eof()) {
    ff >> ene >> val;
    if (ff.eof()) break;
    fArgonEnergy.push_back(ene * keV);
    fArgonELoss.push_back(val * 10);
  }
  ff.close();

  ff.open("../data/physics/alpha_eloss.dat");
  while (!ff.eof()) {
    ff >> ene >> val;
    if (ff.eof()) break;
    fAlphaEnergy.push_back(ene * keV);
    fAlphaELoss.push_back(val * keV / micrometer);
  }
  ff.close();

  ff.open("../data/detector/s2_radial_correction.dat");
  while (ff >> pos >> val) {
    fS2CorrectionPosition.push_back(pos * cm);
    fS2CorrectionEfficiency.push_back(val);
  }
  ff.close();

  ff.open("../data/physics/NR_LEff_ARIS.dat");
  while (ff >> pos >> val) {
    LEff_ARIS_energies.push_back(pos);
    LEff_ARIS_values.push_back(val);
  }
  ff.close();

  ff.open("../data/physics/NR_recombination_200Vpcm_ARIS.dat");
  while (ff >> pos >> val) {
    recombination_ARIS_200_energies.push_back(pos);
    recombination_ARIS_200_values.push_back(val);
  }
  ff.close();

  // READ LAr Scintillation Properties

  // Default values, Real ones are Read in ReadData() -->
  // ../data/physics/LArScintillationProperties.txt
  fLindhardFactor = 0.25;  // not used anymore after tuning with AmBe
  fMeanQuantaEnergy = (19.5 * eV);
  fExcitationRatioER = 0.21;
  fExcitationRatioNR = 1.;
  fPhotonEneMean = 9.69 * eV;   // 9.81*eV;    lambda = 126.8 nm
  fPhotonEneWidth = 0.22 * eV;  // 0.60*eV;    FWHM = 7.8 nm
  fS1ScaleFactor = 3.5;
  fD_T = 4.8 * cm2 / s;  // Transverse Diffusion from ICARUS NIM A527 (2004) 329
  fD_L = 18 * cm2 / s;   // Longitudinal Diffusion from atrazhev  Timoshkim theory at 1kV/cm
  fD_L = 5 * cm2 / s;    // PA updated, based on DS50 electroluminescence paper
  fLArTauFast = 6 * ns;
  fLArTauSlow = 1600 * ns;
  fGArTauFast = 11 * ns;
  fGArTauSlow = 3200 * ns;
  fGArSingletToTriplet = 0.09;
  fSigma0_Sq = 0.00;
  fTauReco = 0.8 * ns;   // Kubota Phys Rev B 20 (1979)
  fTaueLAr = 15.8 * ms;  // electron lifetime
  fLightYield = 7000;    // to be tuned!!!! ph/MeV
  // end of default values

  string mys;
  while (getline(DSIO::Get()->GetStreamLArProperties(), mys)) {
    if (!mys.find("MeanQuantaEnergy")) {
      fMeanQuantaEnergy = atof(mys.substr(mys.find("=") + 1).c_str()) * eV;
      DSLog(routine) << " MeanQuantaEnergy: " << fMeanQuantaEnergy << endlog;
    }
    if (!mys.find("S1ScaleFactor")) {
      fS1ScaleFactor = atof(mys.substr(mys.find("=") + 1).c_str());
      DSLog(routine) << " S1ScaleFactor: " << fS1ScaleFactor << endlog;
    }
    if (!mys.find("S2Yield")) {
      fS2Yield = atof(mys.substr(mys.find("=") + 1).c_str());
      DSLog(routine) << " S2Yield: " << fS2Yield << endlog;
    }
    if (!mys.find("ExcitationRatioER")) {
      fExcitationRatioER = atof(mys.substr(mys.find("=") + 1).c_str());
      DSLog(routine) << " ExcitationRatioER: " << fExcitationRatioER << endlog;
    }
    if (!mys.find("ExcitationRatioNR")) {
      fExcitationRatioNR = atof(mys.substr(mys.find("=") + 1).c_str());
      DSLog(routine) << " ExcitationRatioNR: " << fExcitationRatioNR << endlog;
    }
    if (!mys.find("LindhardFactor")) {
      fLindhardFactor = atof(mys.substr(mys.find("=") + 1).c_str());
      DSLog(routine) << " LindhardFactor: " << fLindhardFactor << endlog;
    }
    if (!mys.find("PhotonEneMean")) {
      fPhotonEneMean = atof(mys.substr(mys.find("=") + 1).c_str()) * eV;
      DSLog(routine) << " PhotonEneMean: " << fPhotonEneMean << endlog;
    }
    if (!mys.find("PhotonEneWidth")) {
      fPhotonEneWidth = atof(mys.substr(mys.find("=") + 1).c_str()) * eV;
      DSLog(routine) << " PhotonEneWidth: " << fPhotonEneWidth << endlog;
    }
    if (!mys.find("D_T")) {
      fD_T = atof(mys.substr(mys.find("=") + 1).c_str()) * cm2 / s;
      DSLog(routine) << " D_T: " << fD_T << endlog;
    }
    if (!mys.find("D_L")) {
      fD_L = atof(mys.substr(mys.find("=") + 1).c_str()) * cm2 / s;
      DSLog(routine) << " D_L: " << fD_L << endlog;
    }
    if (!mys.find("Sigma0_Sq")) {
      fSigma0_Sq = atof(mys.substr(mys.find("=") + 1).c_str()) * mm2;
      DSLog(routine) << " Sigma0_Sq: " << fD_L << endlog;
    }
    if (!mys.find("LArTauFast")) {
      fLArTauFast = atof(mys.substr(mys.find("=") + 1).c_str()) * ns;
      DSLog(routine) << " LArTauFast: " << fLArTauFast << endlog;
    }
    if (!mys.find("LArTauSlow")) {
      fLArTauSlow = atof(mys.substr(mys.find("=") + 1).c_str()) * ns;
      DSLog(routine) << " LArTauSlow: " << fLArTauSlow << endlog;
    }
    if (!mys.find("GArTauFast")) {
      fGArTauFast = atof(mys.substr(mys.find("=") + 1).c_str()) * ns;
      DSLog(routine) << " GArTauFast: " << fGArTauFast << endlog;
    }
    if (!mys.find("GArTauSlow")) {
      fGArTauSlow = atof(mys.substr(mys.find("=") + 1).c_str()) * ns;
      DSLog(routine) << " GArTauSlow: " << fGArTauSlow << endlog;
    }
    if (!mys.find("TauReco")) {
      fTauReco = atof(mys.substr(mys.find("=") + 1).c_str()) * ns;
      DSLog(routine) << " TauReco: " << fTauReco << endlog;
    }
    if (!mys.find("TaueLAr")) {
      fTaueLAr = atof(mys.substr(mys.find("=") + 1).c_str()) * ms;
      DSLog(routine) << " TaueLAr: " << fTaueLAr << endlog;
    }
    if (!mys.find("GArSingletToTriplet")) {
      fTaueLAr = atof(mys.substr(mys.find("=") + 1).c_str());
      DSLog(routine) << " GArSingletToTriplet: " << fGArSingletToTriplet << endlog;
    }
  }
}

//////////////////////////////////////////////
//  Energy Loss Interpolator
//////////////////////////////////////////////
G4double DSLight3::Interpolator(double ene, vector<float>& v1, vector<float>& v2) {

  if (ene < v1.front()) return v2.front();
  if (ene > v1.back()) return v2.back();

  int ndim = int(v1.size());

  for (int i = 0; i < ndim; i++) {
    if (v1[i] > ene) {
      double _mm = (v2[i - 1] - v2[i]) / (v1[i - 1] - v1[i]);
      double q = v2[i] - _mm * v1[i];
      return (_mm * ene + q);
    }
  }
  return 0.;
}

G4bool DSLight3::IsApplicable(const G4ParticleDefinition& aParticleType) {
  aParticleType.IsShortLived();
  // DF:
  // return always true
  // this is crucial to understand the origin and the end of the event
  // DO NOT MODIFY IT
  return true;
}

//////////////////////////////////////////////////////////
//     Return Gaspocket Height Based on Sagging Model
//////////////////////////////////////////////////////////
float DSLight3::SaggedGasPocket(G4ThreeVector position) {
  G4ThreeVector x1 = position;
  float tpcRadi = 179.45;
  float r = sqrt(x1[0] * x1[0] + x1[1] * x1[1]);
  float par[] = {9.83263e-01, -2.16438e-10, 1.19817, -5.97500e-5, -1.15250e-4};
  // float GASGAP = DSParameters::Get()->GetGasPocketThickness();
  float GASGAP = 10 * mm;
  if (r > tpcRadi) { r = par[2] * tpcRadi; }
  float gasgap_sagging = GASGAP * (par[3] * x1[0] + par[4] * x1[1] + par[0] + par[1] * pow(tpcRadi * par[2] * tpcRadi * par[2] - r * r, 2));
  return gasgap_sagging;
}

//////////////////////////////////////////////////////////
//     Return Extraction Field Based on Sagging Model
//////////////////////////////////////////////////////////
float DSLight3::ExtractionField(float gasgap_sagging) {
  float eps_l = 1.52;
  float dl = 0.5;
  float ext = 2.8;
  // float GASGAP = DSParameters::Get()->GetGasPocketThickness();
  float GASGAP = 10 * mm;

  return (GASGAP * 0.1 * eps_l + dl) * ext / (eps_l * gasgap_sagging * 0.1 + dl);
}

//////////////////////////////////////////////////////////////////
//     S2 Correction Factor Based on Sagging Model
//////////////////////////////////////////////////////////////////
float DSLight3::RCorrectionFactor(G4ThreeVector position) {
  float gasgap_sagging = SaggedGasPocket(position);
  float extField = ExtractionField(gasgap_sagging);
  float par1[] = {2.985363, -2.974136, 0.863177, -5.456808e-2};
  // float GASGAP = DSParameters::Get()->GetGasPocketThickness();
  float GASGAP = 10 * mm;
  float factor = (par1[0] + par1[1] * extField + par1[2] * pow(extField, 2) + par1[3] * pow(extField, 3)) * (gasgap_sagging / GASGAP);
  return factor;
}

//////////////////////////////
//       Create S2
//////////////////////////////
void DSLight3::CreateS2(G4ThreeVector position, G4double time, G4ParticleChange* apc) {

  G4ThreeVector x1 = position;
  G4double myT1 = time;

  // set to gas
  /*

  const G4Material* aMaterial = DSMaterial::Get()->GetGaseousArgon();
  G4double Pressure = aMaterial->GetPressure();
  G4MaterialPropertiesTable* aMaterialPropertiesTable =
  aMaterial->GetMaterialPropertiesTable(); G4double ExtractionField =
  DSStorage::Get()->GetExtractionField();

  // in Ar at 3kV/cm all extracted //http://arxiv.org/pdf/astro-ph/0701286.pdf
  //otherwise a probability has to be set

  G4int numberofsecondary;//from http://arxiv.org/pdf/1207.2292.pdf
  G4double phgasA = 0.0813;
  G4double phgasB = 139;
  G4double phgasG = 30.6;
  G4double GASGAP = 0.25*cm;
  G4double phpercm  = phgasA * ExtractionField/(volt/cm) - phgasB*Pressure/bar -
  phgasG; numberofsecondary = G4int(
  G4RandGauss::shoot(phpercm*GASGAP/cm,sqrt(phpercm*GASGAP/cm)) + 0.5 );  //
  floor(YieldFactor) removed because always = 1 // to check numberofsecondary =
  G4int(numberofsecondary/DSStorage::Get()->GetScaleS2() + 0.5) ;
  */

  float meannumberofphotons;
  if (DSParameters::Get()->GetWithSaggingModel()) {
    meannumberofphotons = fS2Yield * RCorrectionFactor(x1);
  } else {
    meannumberofphotons = fS2Yield;
    // If DS50, apply the radial correction to S2
    if (DSEventHandler::Get()->GetDetectorFlag() < 3) {
      double _xyradius = sqrt(pow(x1.x(), 2) + pow(x1.y(), 2));
      meannumberofphotons *= Interpolator(_xyradius, fS2CorrectionPosition, fS2CorrectionEfficiency);
    }
  }

  // PA: scale deafult S2 yield by a factor (g2 tuning)
  meannumberofphotons *= DSStorage::Get()->GetScaleS2();

  //scale by a factor equal to MaxQE or MaxPDE for fast simulations
  meannumberofphotons *= fYieldScalingFactor ;

  double numberofsecondary_mean = meannumberofphotons;
  double numberofsecondary_sig = 1.101 * sqrt(numberofsecondary_mean);
  G4int numberofsecondary = GetRandomNegativeBinomial(numberofsecondary_mean, numberofsecondary_sig);

  if (numberofsecondary < 0) numberofsecondary = 0;

  DSLog(debugging) << "\n Exiting from DS2Light::DoIt -- "
                   << "NumberOfSecondaries = "
                   << " " << numberofsecondary << endlog;

  numberofsecondary = int(numberofsecondary);

  for (int k = 0; k < numberofsecondary; k++) {
    // start particle creation
    G4double sampledEnergy;
    G4DynamicParticle* aQuantum;

    // Generate random direction
    G4double cost = 1. - 2. * G4UniformRand();
    G4double sint = std::sqrt((1. - cost) * (1. + cost));
    G4double phi = twopi * G4UniformRand();
    G4double sinp = std::sin(phi);
    G4double cosp = std::cos(phi);
    G4double px = sint * cosp;
    G4double py = sint * sinp;
    G4double pz = cost;

    // Create momentum direction vector
    G4ParticleMomentum photonMomentum(px, py, pz);

    // Determine polarization of new photon
    G4double sx = cost * cosp;
    G4double sy = cost * sinp;
    G4double sz = -sint;
    G4ThreeVector photonPolarization(sx, sy, sz);
    G4ThreeVector perp = photonMomentum.cross(photonPolarization);
    phi = twopi * G4UniformRand();
    sinp = std::sin(phi);
    cosp = std::cos(phi);
    photonPolarization = cosp * photonPolarization + sinp * perp;
    photonPolarization = photonPolarization.unit();

    // Generate a new photon:
    G4double E_eVLAr = 9.7;
    sampledEnergy = G4RandGauss::shoot(E_eVLAr * eV, 0.2 * eV);
    aQuantum = new G4DynamicParticle(G4OpticalPhoton::OpticalPhoton(), photonMomentum);
    aQuantum->SetPolarization(photonPolarization.x(), photonPolarization.y(), photonPolarization.z());

    // assign energy to make particle real
    aQuantum->SetKineticEnergy(sampledEnergy);
    // G4cout << "sampledEnergy = " << sampledEnergy << G4endl;

    // electroluminesence emission time distribution
    G4double aSecondaryTime = myT1, SingTripRatio = fGArSingletToTriplet;

    if (G4UniformRand() < SingTripRatio) aSecondaryTime -= fGArTauFast * log(G4UniformRand());  // singlet
    else
      aSecondaryTime -= fGArTauSlow * log(G4UniformRand());  // triplet

    G4ThreeVector aSecondaryPosition = x1;
    aSecondaryPosition[2] += (G4UniformRand() * fGASGAP / mm) * mm;
    // if (DSStorage::Get()->Get20KGeometry())  aSecondaryPosition[2] -=
    // GASGAP/2. ;  //PA this is no more required after the planc merge to
    // master

    // add some time, based on electron drift velocity across gap
    G4double eMobilitygAr = 475.0 * cm * cm / s / volt;  // http://journals.aps.org/pr/pdf/10.1103/PhysRev.166.871
    G4double gArDriftField = 2800.0 * volt / cm;         // needs to be updated
    G4double egArDriftVel = eMobilitygAr * gArDriftField;
    aSecondaryTime += (aSecondaryPosition - x1).mag() / egArDriftVel;

    // GEANT4 business: stuff you need to make a new track
    G4Track* aSecondaryTrack = new G4Track(aQuantum, aSecondaryTime, aSecondaryPosition);
    apc->AddSecondary(aSecondaryTrack);
  }
}

//////////////////////////////
//  GetLArDriftVelocity
//////////////////////////////
//G4double DSLight3::GetElectronDriftVelocity(G4float field, G4float T=87.7){ 
G4double DSLight3::GetLArDriftVelocity(G4double tempinput, G4double efieldinput) {
  //
  // Get the electron drift velocity (from https://lar.bnl.gov/properties/trans.html)
  //Args:
  //  field: electric field 
  //  T: temperature in Kelvin
  // Returns:
  //  float: drift velocity in mm/us.
  
  //
  G4float pars[6] = {551.6, 7158.3, 4440.43, 4.29, 43.63, 0.2053};
  G4float T0 =  89*kelvin;
  tempinput  /= kelvin;
  
  G4float field = efieldinput/(kilovolt / cm) ;
  G4double mu = (pars[0] + pars[1]*field + pars[2]*pow(field,3/2.) + pars[3]*pow(field,5./2)) / (1 + pars[1]/pars[0]*field + pars[4]*pow(field,2) + pars[5]*pow(field,3));
  mu *= pow(tempinput/T0,-3./2);
  return mu*field/100*mm/us;
}
// G4double DSLight3::GetLArDriftVelocity(G4double tempinput, G4double efieldinput) {
//   // from ATLAS internal note LARG-NO-058

//   if ( DSStorage::Get()->GetTunedS1At200V() ) return 93000 * cm / s;  // from DS50

//   double p1 = -0.016863;
//   double p2 = -0.0083412;
//   double p3 = 0.18088;
//   double p4 = 8.9751;
//   double p5 = 1.4614;
//   double p6 = 0.32891;
//   double T = tempinput;
//   double T0 = 92.91;
//   double E = efieldinput / (kilovolt / cm);  // in kV/cm

//   double vdrift = (p1 * (T - T0) + 1) * (p3 * E * log(1 + p4 / E) + p5 * pow(E, p6)) + p2 * (T - T0);
//   vdrift *= 1e5;  // cm/s
//   return vdrift*cm/s;


// }



//////////////////////////////
//  Prompt fraction
//////////////////////////////

G4double DSLight3::GetPromptFraction (G4int myZ, G4double InitialKinEne)  {

  G4double mySingletTripletRatio = 0.;


  // -- Electron Recoil
  if (myZ != 18 && myZ != 2) {
    //////////////////////////////
    // f90 from Hinkley fit, as used to compare g4ds to data in the MC paper
    //////////////////////////////

    G4double ene = InitialKinEne / keV  ;
    double f90p0 = 0.249488;
    double f90p1 = 0.146597;
    double f90p2 = -2.89121;
    mySingletTripletRatio =  f90p0 * (1 + exp((-f90p1 * (ene - f90p2))));

  }
  // -- Alphas
  //else if (myZ == 2) {
  //  mySingletTripletRatio = (-0.065492 + 1.9996 * exp(-myDepEne / MeV)) / (1 + 0.082154 / pow(myDepEne / MeV, 2.)) + 2.1811;  // check
  //}
  // -- Nuclear Recoils
  else {

    // here the result fo the fit to true energy (assuming Mei Qunching!)
    // it is ok until npe < ~250
    G4double myTrueRecEne =  InitialKinEne / fQuenchingFactor / keV;  // true recoil energy. InitialKinEne is scaled by fQuenchingFactor above
    if (myTrueRecEne > 180) myTrueRecEne = 180;


    G4double  ratio_p0 = 0.513575;      ///-   0.00935121
    G4double  ratio_p1 = 0.00664834;    ///-   0.000618988
    G4double  ratio_p2 = -7.26861e-05;  ///-   1.22681e-05
    G4double  ratio_p3 = 3.69379e-07;   ///-   9.10873e-08
    G4double  ratio_p4 = -7.04932e-10;  ///-   2.24589e-10

    mySingletTripletRatio = ratio_p0 - 0.045 + ratio_p1 * myTrueRecEne + ratio_p2 * pow(myTrueRecEne, 2) + ratio_p3 * pow(myTrueRecEne, 3) + ratio_p4 * pow(myTrueRecEne, 4);

  }
  return mySingletTripletRatio ;
}


//////////////////////////////
//  Recombination Probability
//////////////////////////////
G4double DSLight3::GetRecoProbAt200Vcm(double InitialKinEne, bool isNR) {
  //  - Extracted from Ar39 spectrum in DS10 data
  //  - 9 parameters (more detailed modelization of the field dependence)
  //  - Improved agreement with DS10 data in the region of the recombination
  //  peak
  double myRecoProb = 0;

  if (!isNR) {
    double p0 = 2.96766e-01;   //  7.50908e-03   0.00000e+00  -2.61950e-01
    double p1 = 3.95496e+00;   // 2.77124e+00  -0.00000e+00  -1.52624e+01
    double p2 = -5.17812e-01;  // 2.69001e-01   0.00000e+00  -1.88046e+02
    double p3 = -1.38485e-02;  // 1.13375e-03  -0.00000e+00   4.47780e+03
    double p4 = 9.12436e-01;   // 1.34073e-02  -0.00000e+00  -9.04345e+02
    double p5 = 6.61046e-01;   // 1.59214e-04  -0.00000e+00   2.03157e+01

    myRecoProb = p0 * (1 - p1 * exp(p2 * InitialKinEne / keV)) * exp((p3)*pow(InitialKinEne / keV, (p4))) + p5;
  } else {
    // recombination for NR from ARIS, at 200 V/cm, as function of visible
    // energy (keVee). measured between 7.1 keVNR (1.7 keVee) and 117 keVNR (41
    // keVee). Extrapolated otherwise from analytical expression in Phys. Rev. D
    // 97, 112005 (2018)
    myRecoProb = Interpolator(InitialKinEne / keV, recombination_ARIS_200_energies, recombination_ARIS_200_values);
  }

  return myRecoProb;
}

//////////////////////////////
//  Nuclear Quenching Factor
//////////////////////////////

G4double DSLight3::GetLArNuclearQuenching(double myene) {
  myene /= keV;

  if (IsMeiQuenching) {
    // model suggested by Mei arXiv:0712.2470
    // tuned on AmBe data (4.4 MeV coincidence in Veto)
    // limit at high energy: 0.28

    if (myene > 250) myene = 250;
    double p0 = 0.172706;
    double p1 = 0.0254609;
    double p2 = -0.000178219;
    double p3 = 0.179962;
    return (p1 * std::log(myene + p0) + p2 * myene + p3);
  } else if (IsLindhardQuenching) {
    // pure Lindhard

    double plin[10];
    plin[0] = 0.191995;      ///-   0.00156303
    plin[1] = 0.00853716;    ///-   0.000255132
    plin[2] = -0.000205518;  ///-   1.23408e-05
    plin[3] = 3.16629e-06;   ///-   2.59109e-07
    plin[4] = -2.76908e-08;  ///-   2.76198e-09
    plin[5] = 1.35493e-10;   ///-   1.56351e-11
    plin[6] = -3.4568e-13;   ///-   4.47693e-14
    plin[7] = 3.57839e-16;   ///-   5.09874e-17
    plin[8] = 0.407796;      ///-   0.00469415
    plin[9] = 0.00056365;    ///-   1.77868e-05

    double qf = 0.;
    if (myene < 250)
      for (int i = 0; i < 8; ++i) qf += plin[i] * pow(myene, i);
    else if (myene < 400)
      qf = plin[8] + myene * plin[9];
    else
      qf = plin[8] + 400 * plin[9];
    return qf;
  } else if (IsARISQuenching) {
    // LEff from ARIS. Measured between 7.1 and 120 keVNR. Extrapolated
    // otherwise according to the expression in Phys. Rev. D 97, 112005 (2018)
    return Interpolator(myene, LEff_ARIS_energies, LEff_ARIS_values);

  } else
    return 0.25;
}

//////////////////////////////
// Get random number from the Gamma distribution (gsl)
//////////////////////////////
double DSLight3::GetRandomGamma(const double a, const double b) {
  /* assume a > 0 */
  if (a < 1) {
    double u = G4RandFlat::shoot(0., 1.);
    return GetRandomGamma(1.0 + a, b) * pow(u, 1.0 / a);
  }
  {
    double x, v, u;
    double d = a - 1.0 / 3.0;
    double c = (1.0 / 3.0) / sqrt(d);

    while (1) {
      do {
        x = G4RandGauss::shoot(0, 1.0);
        v = 1.0 + c * x;
      } while (v <= 0);

      v = v * v * v;
      u = G4RandFlat::shoot(0., 1.);

      if (u < 1 - 0.0331 * x * x * x * x) break;
      if (log(u) < 0.5 * x * x + d * (1 - v + log(v))) break;
    }
    return b * d * v;
  }
}


void DSLight3::CreateClusters (G4double myphotons, G4double myelectrons, G4double mytime, G4ThreeVector myPos, G4double myTrueDepEne, G4double pfraction /*, G4int matIdx*/) {

  //read all the values in units of cm, ns, keV

  G4double dep_x    = myPos[0] / cm;
  G4double dep_y    = myPos[1] / cm;
  G4double dep_z    = myPos[2] / cm;

  G4double dep_S1ene    = (myphotons) * fMeanQuantaEnergy / keV ;
  G4double dep_S2ene    = (myelectrons) * fMeanQuantaEnergy / keV ;
  //G4double dep_ene      = (myphotons + myelectrons) * fMeanQuantaEnergy / keV ;
  G4double dep_true_ene = myTrueDepEne / keV ;

  G4double dep_time = mytime / ns  ;
  //cout << " __ mytime: "<< mytime << "; with dep_time being: " << dep_time << "__.";
  G4bool fMakeNewCluster = true ;
  cout << " size of the VClusters: " << int( DSEventHandler::Get()->GetVClusters().size()  ) << ". " << endl ; 

  for( int j = 0; j < int( DSEventHandler::Get()->GetVClusters().size()  ); j++){
    cout << "Loop entered for j = " << j << endl ;
    G4double deltaZ     =   fabs( DSEventHandler::Get()->GetVClusters()[j].Position[2] - dep_z ) ;
    G4double deltaT     =   fabs( DSEventHandler::Get()->GetVClusters()[j].Time - dep_time );
    G4double deltaR     =   sqrt ( pow (  DSEventHandler::Get()->GetVClusters()[j].Position[0]-dep_x, 2) +  pow (  DSEventHandler::Get()->GetVClusters()[j].Position[1]-dep_y, 2)   ) ;
    cout << "deposit " << j << ": deltaZ = " << deltaZ << "; deltaT = " << deltaT << "; deltaR = " << deltaR << "." << endl ;
    // clustering conditions
    if( deltaZ < fClustering_dist_max_z_cm &&  deltaT < fClustering_deltaT_max_ns  && deltaR < fClustering_deltaR_max_cm ) {
      //weighted average

      DSEventHandler::Get()->GetVClusters()[j].Position[2] = ( DSEventHandler::Get()->GetVClusters()[j].Position[2]*DSEventHandler::Get()->GetVClusters()[j].S2Energy + dep_z*dep_S2ene ) / (dep_S2ene+DSEventHandler::Get()->GetVClusters()[j].S2Energy)  ;
      DSEventHandler::Get()->GetVClusters()[j].Position[0] = ( DSEventHandler::Get()->GetVClusters()[j].Position[0]*DSEventHandler::Get()->GetVClusters()[j].S2Energy + dep_x*dep_S2ene ) / (dep_S2ene+DSEventHandler::Get()->GetVClusters()[j].S2Energy)  ;
      DSEventHandler::Get()->GetVClusters()[j].Position[1] = ( DSEventHandler::Get()->GetVClusters()[j].Position[1]*DSEventHandler::Get()->GetVClusters()[j].S2Energy + dep_y*dep_S2ene ) / (dep_S2ene+DSEventHandler::Get()->GetVClusters()[j].S2Energy)  ;

      G4double cl_prompt = DSEventHandler::Get()->GetVClusters()[j].RecoilID * DSEventHandler::Get()->GetVClusters()[j].S1Energy ;
      G4double cl_late   = DSEventHandler::Get()->GetVClusters()[j].S1Energy  - cl_prompt ;
      cl_prompt          += dep_S1ene * pfraction ;
      cl_late            += (1-pfraction) * dep_S1ene ;

      DSEventHandler::Get()->GetVClusters()[j].RecoilID  =  cl_prompt / (cl_prompt+cl_late)  ;
      DSEventHandler::Get()->GetVClusters()[j].S1Energy  +=  dep_S1ene ;
      DSEventHandler::Get()->GetVClusters()[j].S2Energy  +=  dep_S2ene ;
      DSEventHandler::Get()->GetVClusters()[j].Energy    +=  dep_true_ene ;
      cout << "REUSED cluster for deposit " << j << " and loop exited " << endl ;
      fMakeNewCluster = false;
      break;
    }
  }
  // Otherwise, new cluster
  if( fMakeNewCluster ){

    DSEventHandler::Get()->SetClusterS1Energy  (   myphotons* fMeanQuantaEnergy / keV ) ;
    DSEventHandler::Get()->SetClusterS2Energy  (   myelectrons* fMeanQuantaEnergy / keV  ) ;
    DSEventHandler::Get()->SetClusterEnergy  (   myTrueDepEne   / keV ) ;
    DSEventHandler::Get()->SetClusterTime    (  dep_time ) ;
    DSEventHandler::Get()->SetClusterPosition(myPos / cm )  ;
    //DSEventHandler::Get()->SetClusterMaterial(matIdx)  ;
    DSEventHandler::Get()->SetClusterRecoilID(pfraction)  ;
    DSEventHandler::Get()->SetClusters() ;
    cout << "NEW cluster  " << endl ;
  }

}




//////////////////////////////
// Get random number from the negative binomial distribution (gsl)
//   mu = mean of the distribution
//   sig = rms of the distribution
//////////////////////////////
unsigned int DSLight3::GetRandomNegativeBinomial(double mu, double sig) {

  double p = mu / sig / sig;
  double n = mu * mu / (sig * sig - mu);

  double X = GetRandomGamma(n, 1.0);
  unsigned int k = G4Poisson(X * (1 - p) / p);
  return k;
}

//////////////////////////////
// Return recombination probability scale factor wrt to 200 V/cm 
//    (from the fit of fig. 3 in https://doi.org/10.1016/j.nima.2003.11.423)
// Args:
//   field : electric field in kV/cm
// Returns:
//   float : scale factor
// 
//////////////////////////////


G4float DSLight3:: GetRecombinationWrt200Vcm(G4float field) { 
  if (field == 0.2*kilovolt / cm) return 1;
  float pars[3] = {1.890, 1.2455,  205.1e-3};
  field /= (kilovolt / cm);
  float field0 = 0.2 ;
  float reco_at_200Vcm = pars[0]*field0/(1+pars[0]*field0)/pars[1] + pars[2];
  float reco = pars[0]*field/(1+pars[0]*field)/pars[1] + pars[2];
  return reco_at_200Vcm/reco;
}



