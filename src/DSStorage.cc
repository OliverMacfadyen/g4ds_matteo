#include "DSStorage.hh"
#include <complex>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include "DSLogger.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

using namespace std;

DSStorage* DSStorage::me = 0;

// singleton
DSStorage::DSStorage() {

  fMultiplicator = 1 ;
  fIsDecayInVacuum  = false ;
  fIsEnDepGenerator = false;
  fOverlap = 1;
  fEventCounter = 1000;
  fOverwriteCounter = false ;
  fExportGDML = 0;
  fWritePhotons = 0;
  fVerbosity = -1;
  fWriteDeposits = 0;
  fWriteThermalElectrons = 0;
  fWriteDaughters = 0;
  fWriteFilter = 0;
  fArDMeventsShieldPE = 0;
  fPMT = 1000000;
  fVetoPMT = 1000000;
  fMuPMT = 1000000;
  fRDMDecay = false;
  fReDStacking = false;
  fRealPDGMeanLife = 0;
  fPreAbsTime = 0;
  fNDaughters = 1;
  fRDMChain = false;
  fKillS1S2 = false;
  fKillS2 = false;
  fKillS1 = false;
  fScaleS2 = 1.0;
  fDriftField = 200*volt/cm;       // V/cm
  fExtractionField = 3500*volt/cm;  // V/cm
  fTMBfraction = 0.1;
  fTimeCut = 10 * ms;
  fScintillator = 0;  // 0: BorateScintillator  ; 1 Gd Scintillator ; 2 Li6Scintillator ; 3
                      // NatLi6Scintillator ; 4 OrtoCarbonScintillator
  fCTFfill = 0;       // 0 Water  ; 1 Gd loaded water ; 2 BorateScintillator
  fFastSimulation = 1;
  fIsExternalLArScintillating = false;
  fIsCopperRings = false;
  fVetoYieldFactor = 1.0;
  fTunedS1At200V = true;

  fSourceHolderCenter = G4ThreeVector(-33.85 * cm, 3.75 * cm, -3.65 * cm);
  fSourceHolderTheta = 90 * deg;
  fSourceHolderPhi = 0. * deg;
  fSourceHolderFlag = false;
  fSourceHolderLeadFlag = false;

  fClusterWritingActivated = true;

  fIs5KGeometry = false;
  fIs20KGeometry = false;
  fIsArDMGeometry = false;
  fIsDartGeometry = false;
  fIsReDGeometry = false;
  fIsPETGeometry = false;
  fIsDS20kLSV = false;
  fDS20kTPCverticalShift = 15 * cm;

  fIsCherenkov = false;
  fIsEmittedCherenkov = false;

  fIsCustomCryostat = false;
  fDS20kReflectorAlternate = 0;
  fDS20kCryoCornerDistance = 5 * cm;  // deprecated
  fDS20kTPCheight = 350.5 * cm;
  fDS20kTPCedge = 145.1 * cm;
  fArDMLeadHeight = -1 * cm;              // Lead Height: negative by default (NO LEAD RING)
  fDS20kWindowMaterial = 0;               // 0: acrylic, 1: quartz
  fDS20kLArThicknessAboveTPC = 8.5 * cm;  // amount of NonScint LAr above the TPC
  fDS20kGasPocketThickness = 0.7 * cm;
  fSiPMOffset = 1 * cm;  // updated to planc HA version 1.*cm;
  fDS20kTPCVesselThickness = 15 * cm;
  fDS20kHDPEShellThickness = 5*cm ;
  fDS20kHDPEShellCapThickness = 5.*cm ; 
  fDS20kTPBThickness = 3 * um;
  fDS20kTPBLayers = "1111";  // Where the layers of TPB go, left most digit must
                             // be a 1, the next digit is closet to vessel, right
                             // most is outermost layer. 1= true, 0 = false
  fGridSurfaceDistance = 3*mm;
  fGdLayerThickness = 0 * cm;
  fDS20kRemoveBondingRing = 1;

  fDS20kWLSLArThickness = 1000 * um;
  fDS20kWLSPENThickness = 25 * um;

  fWriteFilterMinimumEnergy = 0;    // keV
  fWriteFilterMaximumEnergy = 1e9;  // keV

  // Build and place SiPMs or not
  fDS20kSiPMs = true;

  // Gd concentration - N.B: Gd concentration, not Gd2O3 one.
  fDS20kGdConcentration = 0.01;

  // Number of vPDUs to be placed. NOT to be changed.
  // Using an input of 100, the code automatically places 120 vPDUs
  fDS20knSiPMs = 100;

  // SiPMs will not be placed on the outer surfaces of the TPC but at some
  // distance
  fDS20kepsilonSiPM = 0.05 * cm;

  // SiPM side and height
  fDS20kSiPMSide = 5.0 * cm;
  fDS20kSiPMHeight = 0.5 * cm;

  fG3WTankR = 5 * m;
  fG3IDiameter = 800 * cm;
  fG3WTankHeight = 9 * m;
  fAcrylicVesselDiameter = 0 * cm;

  fIsLicNdAngleByAng = true;
  fIsLicNdAngleByEnergy = false;

  fLicNDNuclRecE = 1 * MeV;
  fLicNDTheta = 60 * degree;
  fLicNDPhi1 = 30 * degree;
  fLicNDPhi2 = 210 * degree;
  fLicNDRadius = 1 * m;

  fLicNeutronE = 1.5 * MeV;

  fLicWallActivated = false;
  fLicWallThick = 5 * cm;
  fLicWallWindowLength = 10 * cm;
  fLicWallDistanceToTPC = 5 * cm;

  fDartTeflonThickness = 0.;
  fDartTeflonHeight = 0.;
  fDartTeflonRadius = 0.;

  fReDConfiguration = 0;  // 0: Napoli 65 cm cryostat, 99: ARIS like cryostat

  fIsProtoProtoBottomSiPM = 0;

  fLicExpHall = false;

  // ARIS-ER parameters
  fIsAriserGeometry = false;  // Set  in DSDetectorARISERGamma.cc
  fAriserSourceDistanceToCryo = 5 * cm;
  fAriserBEGeDistanceToCryo = 3 * cm;
  fAriserBEGeAngle = 20 * deg;
  fAriserSourcePosition = G4ThreeVector(0, 0, 0);  // Set in DSDetectorARISERGamma.cc

  fFLUKAfilename = "../data/physics/fluka_output_example.dat";

  fWL[0] = 300.0;
  fRI[0] = 2.375;
  fEC[0] = 0.247;
  fWL[1] = 350.0;
  fRI[1] = 2.292;
  fEC[1] = 0.082;
  fWL[2] = 400.0;
  fRI[2] = 2.182;
  fEC[2] = 0.045;
  fWL[3] = 450.0;
  fRI[3] = 2.1;
  fEC[3] = 0.021;
  fWL[4] = 500.0;
  fRI[4] = 2.06;
  fEC[4] = 0.016;
  fWL[5] = 550.0;
  fRI[5] = 2.05;
  fEC[5] = 0.014;
  fWL[6] = 600.0;
  fRI[6] = 2.04;
  fEC[6] = 0.012;
  fWL[7] = 650.0;
  fRI[7] = 2.03;
  fEC[7] = 0.011;
  fWL[8] = 700.0;
  fRI[8] = 2.02;
  fEC[8] = 0.0105;
  fWL[9] = 750.0;
  fRI[9] = 2.01;
  fEC[9] = 0.0105;
  fWL[10] = 800.0;
  fRI[10] = 1.914;
  fEC[10] = 0.01;

  fSourcePosition = G4ThreeVector(60.0 * cm, 0, 0);

  ifstream fin("../data/detector/GridTransm.dat");

  double xx, tt;
  for (int i = 0; i < 100; ++i) {
    fin >> xx >> tt;
    fAng.push_back(xx / 360 * 2. * M_PI) ; //  TMath::Pi());
    // PDM 5/21/2014: just store the value directly.
    //    fTrans.push_back(1-(1-tt)/2.);
    fTrans.push_back(tt);
  }
  fin.close();

  ifstream feff("../data/detector/tpb_efficiency.dat");

  for (int iy = 0; iy < 100; iy++) {
    for (int ix = 0; ix < 100; ix++) { feff >> tpb_eff[ix][iy]; }
  }
  feff.close();

  ////////////////////////////////////////////
  // Calibration pipes
  fGantryConfiguration = 0;
  fGantryGlobalPosition = G4ThreeVector();
  fGantryPipeID = 3.0 * cm;
  fGantryPipeWall = 1.5 * mm;
  fGantryPipeBendingR = 40.0 * cm;
  fGantryDistanceTPCside = 3.0 * cm;
  fGantryDistanceTPCbottom = 8.25 * cm;
  fGantryDistanceBtwPipes = 3.0 * mm;
  fGantrySurface = 0;
  fGantrySurfaceOption = 2;
  fGantryShield = 0;
  fGantryShieldThickness = 10.0 * cm;
  fGantryShieldHeight = 10.0 * cm;
  fGantryShieldOffset = 1.0 * cm;
  fGantryShieldMaterial = 1;

  // plan c vessel - names refers to a cryostat
  fDS20kCryoWallThick = 0.8 * cm;
  fDS20kCryoCapWallThick = 1.2 * cm;
  fDS20kCryoMaterial = 0;  // 0: stainless steel, 1: titanium
  fDS20kCryoR = 2333 * mm;
  fDS20kCryoH = 3573 * mm;
  fDS20KCryoBarrelH = 3573 * mm;
  fDS20kCryoBottomCap = 906 * mm;
  fDS20kCryoTopCap = 906 * mm;
  fDS20kCryoTopOffset = 138 * mm;
  fDS20kCryoBottomOffset = 108 * mm;
  // old from double wall cryostat
  fDS20kCryoBottomCapOut = 806 * mm;
  fDS20kCryoTopCapOut = 806 * mm;
  fDS20kCryoTopOffsetOut = 275.5 * mm;
  fDS20kCryoBottomOffsetOut = 275.5 * mm;
  fDS20kCryoVacuumTh = 80 * mm;

  //////////////////////////////////////////////

  ///////     OLD - Plan A variables     ///////

  // Number of SiPMs to be placed.
  fDS20knSiPMInside = 2000;
  fDS20knSiPMOutside = 1000;

  fDS20kSiPMsAutoplacement = true;
  fDS20kSiPMsAutoplacementFilename = "manual_SiPM_Placement.txt";

  // Fraction of SiPMs on caps respect to the total. If different from zero, the
  // code will use them instead of the ones calculated.
  fDS20kSiPMFractionInside = 0;
  fDS20kSiPMFractionOutside = 0.4;

  // TPB on SiPMs or not
  fDS20kTPBOnSiPM = false;

  // Uniform distribution or not for SiPMs
  fDS20kSiPMUniformInsideSides = true;
  fDS20kSiPMUniformInsideCaps = true;
  fDS20kSiPMUniformOutsideSides = false;
  fDS20kSiPMUniformOutsideCaps = false;

  // Best configuration for uniform light collection
  // The settings for UV absorption length and TPB visible reflectivity need to
  // be set by hand.
  fDS20kbestLightCollection = true;

  // TPB panels thichness
  fDS20kVerticalInnerTPBPanely = 1. * cm;
  fDS20kHorizontalInnerTPBPanely = 1. * cm;
  fDS20kVerticalOuterTPBPanely = 1. * cm;
  fDS20kHorizontalOuterTPBPanely = 1. * cm;

  // Plastic veto offsets
  fDS20kinternalOffset = 0.5 * mm;  // Offset of TPB panels on faces from inner plastic TPB
  fDS20kexternalOffset = 0.5 * mm;  // Offset of TPB panels on faces from outer plastic TPB

  fDS20kcenterInnerOffset = 3. * cm;  // Offset for cilindrical inner volume on z axis
                                      // If too small, inner horizontal panels will overlap at the center!

  fDS20kcenterOuterOffset = 3. * cm;  // Offset for cilindrical outer volume on z axis
                                      // If too small, outer horizontal panels will overlap at the center!

  fDS20kupperHorizontalInnerPanelOffset = 0.5 * mm;  // Upper horizontal inner TPB panel offset from electronics or,
                                                     // in general,  whatever volume is placed under them
  fDS20klowerHorizontalInnerPanelOffset = 0.5 * mm;  // Lower horizontal inner TPB panel offset from electronics or,
                                                     // in general,  whatever volume is placed under them
  fDS20kupperHorizontalOuterPanelOffset = 0.5 * mm;  // Upper horizontal outer TPB panel offset from copper TPB
  fDS20klowerHorizontalOuterPanelOffset = 0.5 * mm;  // Lower horizontal outer TPB panel offset from copper TPB

  fDS20kverticalInnerPanelOffset = 0.5 * mm;  // Offset of vertical inner TPB panels from TPC
  fDS20kverticalOuterPanelOffset = 0.5 * mm;  // Offset of vertical outer TPB panels from copper TPB

  fDS20kplasticLArBufferThickness = 1. * cm;  // Thickness of LAr buffer between two plastic solids. Must be
                                              // lower than plastic thickness

  fDS20kinnerCylinderOffset = 0.5 * mm;       // Offset of inner TPB cylinders from horizontal panels
  fDS20kouterCylinderOffset = 0.5 * mm;       // Offset of outer TPB cylinders from horizontal panels
  fDS20kinnerTPBCylinderThickness = 1. * cm;  // Inner TPB cylinder thickness. Must be lower than cylinder
                                              // outer radius
  fDS20kouterTPBCylinderThickness = 1. * cm;  // Outer TPB cylinder thickness. Must be lower than cylinder
                                              // outer radius

  fDS20kCopperThickness = 1. * mm;  // Thickness of copper shell

  fDS20kFacePanels = false;  // To set panels on faces
  fDS20kEdgePanels = true;   // To set panels on edges

  fDS20kBufferPercentageCoverageInside = 1;  // percentage coverage of LAr buffer covered by SiPMs
  fDS20kBufferPercentageCoverageOutside = 1;

  fDS20kLArBufferThickness = 40 * cm;
  fDS20kGdPlasticThickness = 15 * cm;

  fSkipEventsWithNoDaughters = false ;
  fAbortRun = false ;
  fMatScintillation = false ;

  fPurePMMATPC = false;
  fHybridTPC = false;
}

DSStorage* DSStorage::Get() {
  if (!me) me = new DSStorage();

  return me;
}

G4double DSStorage::GetTPBEfficiency(G4double x, G4double y) {

  int bx, by;
  bx = (int)(-y * 0.25 + 50);
  by = (int)(x * 0.25 + 50);
  return tpb_eff[bx][by];
}

G4double DSStorage::GetGridParameters(G4double th0) {

  // PDM 5/21/2014: No redefinition of angle of incidence needed.
  //  th0 -= TMath::Pi()/2.;
  double inc = 0;
  inc = th0;

  // PDM: having trouble with normal incidence: returns xa=xb=0.
  if (inc == 0) return 1.0;

  // Interpolate from table
  std::vector<double>::iterator it = fAng.begin();
  std::vector<double>::iterator tr = fTrans.begin();
  while (*it < inc) {
    ++it;
    ++tr;
  }
  float xa, xb, ya, yb;
  xb = *it;
  --it;
  xa = *it;
  yb = *tr;
  --tr;
  ya = *tr;

  double slope = (yb - ya) / (xb - xa);
  fGridParameter = slope * inc + ya - xa * slope;
  // pdm  std::cout<<"***DSStorage: xa,xb,ya,yb= "<<xa<<"  "<<xb<<"  "<<ya<<"
  // "<<yb<<std::endl;
  return fGridParameter;
}

G4double DSStorage::GetPENEfficiency(G4double x, G4double y) {
  int bx, by;
  bx = (int)(-y * 0.25 + 50);
  by = (int)(x * 0.25 + 50);
  return pen_eff[bx][by];
}

/*
G4double* DSStorage::GetITOParameters(G4double n1, G4double n3, G4double d, G4double lambda, G4double th0) {

  // PDM 5/21/2014: No redefinition of angle of incidence needed.
  //  th0 -= TMath::Pi()/2.;
  // n1 materiale 1
  // n3 materiale 2
  // spessore in nm
  // lambda wl luce incidente
  // th0 theta angolo incidenza in radianti rispetto alla normale

  //  Double_t i_n1, i_n3, i_th0, i_d, i_lambda;

  G4double n = 0;
  G4double k = 0;

  TComplex n2;
  TComplex th1, th2, th3;
  TComplex ct1, st1, ct2, st2, ct3, st3;
  TComplex delta;

  // i_n1 = 1.233; i_n3 = 1.49;
  // i_n=2.182; i_k=0.045;
  // i_d=100.0; i_lambda = 420;
  // i_th0 = 0.0;

  //  std::cout << n << "  " << k << std::endl;

  int np = 11;
  if (lambda < fWL[0] || lambda > fWL[np - 1]) {
    fITOParameters[0] = -1;
    fITOParameters[1] = -1;
    fITOParameters[2] = -1;
    return fITOParameters;
  }

  for (int i = 0; i < np; i++) {
    if ((lambda >= fWL[i]) && (lambda < fWL[i + 1])) {
      n = fRI[i + 1] + (fRI[i] - fRI[i + 1]) * (lambda - fWL[i + 1]) / (fWL[i] - fWL[i + 1]);
      k = fEC[i + 1] + (fEC[i] - fEC[i + 1]) * (lambda - fWL[i + 1]) / (fWL[i] - fWL[i + 1]);
      break;
    }
  }

  n2(n, k);

  th1(th0, 0.0);
  // PDM 7/25/14: Debugging to make this code give same answers as stand-alone
  // Matlab code it was based on.
  //  Main problem was in definition of ct3, which is imaginary when there is
  //  total internal reflection. Secondary problem was in deriving ct2 and ct3
  //  from st2 and st3, where there is apparently a difference in the sign
  //  convention between the root TComplex class and Matlab.  I got around this
  //  by an empirical fix (note the Conjugate in ct3 but NOT in ct2), that I
  //  compared to the Matlab result for many cases.

  //  th2 = TComplex::ASin((n1*TMath::Sin(th0))/n2);
  //  th3(TMath::ASin((n1*TMath::Sin(th0))/n3), 0.0);

  ct1 = TComplex::Cos(th1);
  st1 = TComplex::Sin(th1);

  //  ct2 = TComplex::Cos(th2);
  //  st2 = TComplex::Sin(th2);
  st2 = st1 * (n1 / n2);
  ct2 = TComplex::Sqrt(1. - st2 * st2);

  //  ct3 = TComplex::Cos(th3);
  //  st3 = TComplex::Sin(th3);
  st3 = st1 * (n1 / n3);
  ct3 = TComplex::Sqrt(1. - st3 * st3);
  ct3 = TComplex::Conjugate(ct3);

  delta = (2.0 * TMath::Pi() * d / lambda) * (n2 * ct2);

  TComplex t12p, t23p, r12p, r23p;
  TComplex t12s, t23s, r12s, r23s;

  t12s = (2.0 * n1 * ct1) / (n1 * ct1 + n2 * ct2);
  t23s = (2.0 * n2 * ct2) / (n2 * ct2 + n3 * ct3);

  t12p = (2.0 * n1 * ct1) / (n2 * ct1 + n1 * ct2);
  t23p = (2.0 * n2 * ct2) / (n3 * ct2 + n2 * ct3);

  r12s = (n1 * ct1 - n2 * ct2) / (n1 * ct1 + n2 * ct2);
  r23s = (n2 * ct2 - n3 * ct3) / (n2 * ct2 + n3 * ct3);

  r12p = (n2 * ct1 - n1 * ct2) / (n2 * ct1 + n1 * ct2);
  r23p = (n3 * ct2 - n2 * ct3) / (n3 * ct2 + n2 * ct3);

  TComplex M11p, M21p;
  TComplex Rptotal, Tptotal;

  M11p = TComplex::Exp(-TComplex::I() * delta) + r12p * r23p * TComplex::Exp(TComplex::I() * delta);
  M11p /= t12p * t23p;

  M21p = r12p * TComplex::Exp(-TComplex::I() * delta) + r23p * TComplex::Exp(TComplex::I() * delta);
  M21p /= t12p * t23p;

  Rptotal = M21p / M11p;
  Tptotal = 1.0 / M11p;

  TComplex M11s, M21s;
  TComplex Rstotal, Tstotal;

  M11s = TComplex::Exp(-TComplex::I() * delta) + r12s * r23s * TComplex::Exp(TComplex::I() * delta);
  M11s /= t12s * t23s;

  M21s = r12s * TComplex::Exp(-TComplex::I() * delta) + r23s * TComplex::Exp(TComplex::I() * delta);
  M21s /= t12s * t23s;

  Rstotal = M21s / M11s;
  Tstotal = 1.0 / M11s;

  Double_t Rs, Rp, Ts, Tp;

  Rs = Rstotal.Rho2();
  Rp = Rptotal.Rho2();
  Ts = Tstotal.Rho2() * (n3 * ct3.Re()) / (n1 * ct1.Re());
  Tp = Tptotal.Rho2() * (n3 * TComplex::Conjugate(ct3).Re()) / (n1 * ct1.Re());

  G4double RR = (Rs + Rp) / 2.0;  // riflettivita'   -> riflessione speculare
  G4double TT = (Ts + Tp) / 2.0;  // trasmissivita'  -> stessa direzione
  G4double AA = 1.0 - (TT + RR);  // aa assorbimento -> muore

  fITOParameters[0] = RR;
  fITOParameters[1] = TT;
  fITOParameters[2] = AA;
  return fITOParameters;
}

*/

G4double *DSStorage::GetITOParameters(G4double n1, G4double n3, G4double d,
G4double lambda, G4double th0) {

  //PA see the function above for description.
  //Restoring the version which doesn't use TComplex. It was found in the past
  //a possible reduction of the LY (DS50). Running the two functions for a couple
  //of values returns however the same numerical outputs.

  //  Double_t i_n1, i_n3, i_th0, i_d, i_lambda;

  G4double n = 0;
  G4double k = 0;

  complex<G4double> n2;
  complex<G4double> th1; //th2,th3;
  complex<G4double> ct1,st1,ct2,st2,ct3,st3;
  complex<G4double> delta;


  // i_n1 = 1.233; i_n3 = 1.49;
  // i_n=2.182; i_k=0.045;
  // i_d=100.0; i_lambda = 420;
  // i_th0 = 0.0;

  //  std::cout << n << "  " << k << std::endl;

  int np = 11;
  if(lambda<fWL[0] || lambda>fWL[np-1]) {
    fITOParameters[0] = -1;
    fITOParameters[1] = -1;
    fITOParameters[2] = -1;
    return fITOParameters;
  }

  for(int i=0; i<np; i++) {
    if( (lambda>=fWL[i]) && (lambda<fWL[i+1]) ) {
      n = fRI[i+1]+(fRI[i]-fRI[i+1])*(lambda-fWL[i+1])/(fWL[i]-fWL[i+1]);
      k = fEC[i+1]+(fEC[i]-fEC[i+1])*(lambda-fWL[i+1])/(fWL[i]-fWL[i+1]);
      break;
    }
  }


  n2  = complex<double>(n, k);

  th1 = complex<double>( th0 , 0.0 );

  ct1 = cos(th1);
  st1 = sin(th1);

  //  ct2 = cos(th2);
  //  st2 = sin(th2);
  st2 = st1*(n1/n2);
  ct2 = sqrt(complex<G4double>(1.,0) - st2*st2);

  //  ct3 = cos(th3);
  //  st3 = sin(th3);
  st3 = st1*(n1/n3);
  ct3 = sqrt( complex<G4double>( 1.,0) -st3*st3);
  ct3 = conj(ct3);

  delta = (2.0*M_PI*d/lambda) * (n2 * ct2);

  complex<G4double> t12p,t23p,r12p,r23p;
  complex<G4double> t12s,t23s,r12s,r23s;

  t12s = (2.0 * n1 * ct1) / (n1*ct1 + n2*ct2);
  t23s = (2.0 * n2 * ct2) / (n2*ct2 + n3*ct3);

  t12p = (2.0 * n1 * ct1) / (n2*ct1 + n1*ct2);
  t23p = (2.0 * n2 * ct2) / (n3*ct2 + n2*ct3);

  r12s = (n1*ct1 - n2*ct2) / (n1*ct1 + n2*ct2);
  r23s = (n2*ct2 - n3*ct3) / (n2*ct2 + n3*ct3);

  r12p = (n2*ct1 - n1*ct2) / (n2*ct1 + n1*ct2);
  r23p = (n3*ct2 - n2*ct3) / (n3*ct2 + n2*ct3);

  complex<G4double> M11p,M21p;
  complex<G4double> Rptotal, Tptotal;

  M11p = exp(-complex<G4double> (0,1)*delta) + r12p*r23p*exp(complex<G4double>
(0,1)*delta); M11p /= t12p*t23p;

  M21p = r12p*exp(-complex<G4double> (0,1) *delta) + r23p*exp(complex<G4double>
(0,1) *delta); M21p /= t12p*t23p;

  Rptotal = M21p/M11p;
  Tptotal = 1.0/M11p;

  complex<G4double> M11s,M21s;
  complex<G4double> Rstotal, Tstotal;

  M11s = exp(-complex<G4double> (0,1) *delta) + r12s*r23s*exp(complex<G4double>
(0,1) *delta); M11s /= t12s*t23s;

  M21s = r12s*exp(-complex<G4double> (0,1) *delta) + r23s*exp(complex<G4double>
(0,1) *delta); M21s /= t12s*t23s;

  Rstotal = M21s/M11s;
  Tstotal = 1.0/M11s;

  G4double Rs,Rp,Ts,Tp;

  Rs = (Rstotal * conj(Rstotal)).real() ;
  Rp = (Rptotal * conj(Rptotal)).real() ;
  Ts = (Tstotal * conj(Tstotal)).real() * (n3*ct3.real())/(n1*ct1.real());
  Tp = (Tptotal * conj(Tptotal)).real() * (n3*conj(ct3).real())/(n1*ct1.real());

  G4double RR = (Rs+Rp)/2.0;    // riflettivita'   -> riflessione speculare
  G4double TT = (Ts+Tp)/2.0;    // trasmissivita'  -> stessa direzione
  G4double AA = 1.0 - (TT+RR);  // aa assorbimento -> muore


  fITOParameters[0] = RR;
  fITOParameters[1] = TT;
  fITOParameters[2] = AA;
  return fITOParameters;

}



/*
 * $Log: DSStorage.cc,v $
 * Revision 1.27  2016/03/08 21:09:36  swesterd
 * Added source holder
 *
 * Revision 1.26  2016/03/03 23:19:42  cz2
 * Added a function to read the tpb_efficiency_map for defect study.
 *
 * Revision 1.25  2016/02/05 16:48:22  pagnes
 * basic PET geometry added
 *
 * Revision 1.24  2015/12/08 15:03:33  cgiganti
 * add experimental hall to licorne
 *
 * Revision 1.23  2015/12/07 09:49:41  riffard
 * Licorne : wall shape modification
 *
 * Revision 1.22  2015/12/03 15:59:29  riffard
 * Licorne : wall addition
 *
 * Revision 1.21  2015/11/26 13:30:29  dfranco
 * licorne
 *
 * Revision 1.20  2015/10/31 11:27:00  pagnes
 * added commands to customize WT, DS20k-cryostats and DS20k-LSV geometries
 *
 * Revision 1.19  2015/10/28 15:46:03  pagnes
 * added variables to build a custom cryostat
 *
 * Revision 1.18  2015/10/22 10:18:41  dfranco
 * added Cherenov process (by default off)
 *
 * Revision 1.17  2015/10/15 09:24:37  dfranco
 * all units updated: position in cm and energy in keV
 *
 * Revision 1.16  2015/10/12 16:11:31  pagnes
 * default time cut changed from 1 ms to 10 ms
 *
 * Revision 1.15  2015/10/06 14:49:41  pagnes
 * realistec field rings (first try) added with /ds/detector/true_field_rings
 * (default is without)
 *
 * Revision 1.14  2015/07/07 12:11:44  dfranco
 * re-committed the old ITO function; the exchange between complex and TComplex
 * reduces the LY
 *
 * Revision 1.13  2015/07/07 07:12:09  dfranco
 * TComplex substituted with complex class; gridparameter function cleaned from
 * unused parameters
 *
 * Revision 1.12  2015/04/23 14:04:02  pagnes
 * DS20K geometry added (config 10)
 *
 * Revision 1.11  2015/01/14 16:58:35  dfranco
 * added source vial and correspondent commands by Laura and Bianca; manual
 * updated
 *
 * Revision 1.10  2014/12/22 14:40:43  dfranco
 * added the option to activate the recombination probability at 200 V/cm
 * (/ds/physics/tunedS1); this option is by default true; selecting a specific
 * drift field automatically switch off the tunedS1 option
 *
 * Revision 1.9  2014/11/21 10:18:59  dfranco
 * added a command to scale the veto scintillation yield factor and fixed the
 * visible energy variable in the veto
 *
 * Revision 1.8  2014/11/20 15:32:05  dfranco
 * added a command to remove scintillation process from liquid argon between TPC
 * and cryostat
 *
 * Revision 1.7  2014/11/06 17:39:45  dfranco
 * added source and fixed a bug: removed gaseous nitrogen from the neutron veto
 *
 * Revision 1.6  2014/07/25 17:15:39  meyers
 * Fixed a bug in ITO code.  Now gives same results as the standalone code it
 * was based on
 *
 * Revision 1.5  2014/07/23 14:52:41  pagnes
 * write thermal e- and kill S1 commands added
 *
 * Revision 1.4  2014/07/16 08:23:04  pagnes
 * QE scaling to 1.0 added (/ds/manager/fast_simulation xxx)
 *
 * Revision 1.3  2014/06/03 13:31:35  meyers
 * Migrate TPC grid and ITO optics updates to g4ds10
 *
 * Revision 1.17  2014/05/28 14:33:40  meyers
 * Remove 90-degree rotation in ITO angle of incidence
 *
 * Revision 1.16  2014/05/23 17:55:49  meyers
 * Update grid optical model
 *
 * Revision 1.15  2014/04/11 12:33:29  pagnes
 * command to set TMB/PC ratio inside veto added
 *
 * Revision 1.14  2014/03/19 16:37:27  dfranco
 * update external configuration reader for optics tuning
 *
 * Revision 1.13  2014/03/11 09:54:38  meregaglia
 * Added generator starting from energy deposits
 *
 * Revision 1.12  2014/01/29 13:13:41  perassos
 * Update of the electric field handling and of the Nuclear Recoils generator
 *
 * Revision 1.11  2014/01/07 14:10:36  perassos
 * Added the commands to set in the macfile the electric field and the
 * Thomas-Imel parameters
 *
 * Revision 1.10  2013/11/19 10:33:20  perassos
 * Added methods to handle the electric field and the liquid/gas interface z
 * coordinate
 *
 * Revision 1.9  2013/07/24 09:49:02  dfranco
 * Added S1 and S2 equivalent energies, the material index for each deposit, the
 * command killS1S2 to kill photons and electrons generated by DSLight (after
 * storing the equivalent energies)
 *
 * Revision 1.8  2013/06/15 07:16:33  dfranco
 * Added spacename TMath for using TComplex for the ITO optical surface
 *
 * Revision 1.7  2013/06/11 22:48:33  dfranco
 * Added ITO optical boundary. The Sernelius function is defined in DSStorage,
 * and called by G4OpBoundaryProcess. New ITO bool variable added to
 * DSDetectorDS50.cc as surface property (G4MaterialPropertyTable)
 *
 * Revision 1.6  2013/06/10 14:15:39  dfranco
 * Added two commands: /ds/physics/killS2 and /ds/physics/scaleS2 to kill or
 * scale the S2 light
 *
 * Revision 1.5  2013/04/03 10:14:25  dfranco
 * Fixed bugs with RDM and RDMChain staking actions. The logic of Geant4 is
 * changed. Different excited states of a nucleus correspond to new particles
 * (trackID). Code adapted.
 *
 * Revision 1.4  2013/03/22 14:09:40  dfranco
 * Added the cvs logger code at the end of each file
 *
 */
