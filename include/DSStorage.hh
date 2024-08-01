#ifndef _DSStorage_HH
#define _DSStorage_HH 1

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "DSLogger.hh"
#include "G4String.hh"
#include "G4ThreeVector.hh"

using namespace std;

class DSStorage {
 private:
  DSStorage();

 public:
  static DSStorage* Get();

  virtual ~DSStorage() {}

 private:
  static DSStorage* me;

 public:
  enum NDPosMode { std, kinematic, efficiency };

  enum NDSurveyMode { xyz, spherical };

  inline void SetIsEnDepGenerator(G4bool val) { fIsEnDepGenerator = val; }
  inline G4int GetIsEnDepGenerator() { return fIsEnDepGenerator; }

  inline void SetCheckOverlap(G4int val) { fOverlap = val; }
  inline G4int GetCheckOverlap() { return fOverlap; }

  inline void SetEventCounter(G4int val) { fEventCounter = val; }
  inline G4int GetEventCounter() { return fEventCounter; }

  inline void SetExportGDML(G4bool val) { fExportGDML = val; }
  inline G4bool GetExportGDML() { return fExportGDML; }

  inline void SetWritePhotons(G4bool val) { fWritePhotons = val; }
  inline G4bool GetWritePhotons() { return fWritePhotons; }

  inline void SetVerbosity(G4int val) { fVerbosity = val; }
  inline G4int GetVerbosity() { return fVerbosity; }

  inline void SetBoronScintillatorIndex(G4int val) { fBoronScintillatorIndex = val; }
  inline G4int GetBoronScintillatorIndex() { return fBoronScintillatorIndex; }

  inline void SetLiquidArgonIndex(G4int val) { fLiquidArgonIndex = val; }
  inline G4int GetLiquidArgonIndex() { return fLiquidArgonIndex; }

  inline void SetDS20kReflectorAlternate(G4int val) { fDS20kReflectorAlternate = val; }
  inline G4int GetDS20kReflectorAlternate() { return fDS20kReflectorAlternate; }

  inline void SetGaseousArgonIndex(G4int val) { fGaseousArgonIndex = val; }
  inline G4int GetGaseousArgonIndex() { return fGaseousArgonIndex; }

  inline void SetLArAboveGridIndex(G4int val) { fLArAboveGridIndex = val; }
  inline G4int GetLArAboveGridIndex() { return fLArAboveGridIndex; }

  inline void SetPMTMaterialIndex(G4int val) { fPMT = val; }
  inline G4int GetPMTMaterialIndex() { return fPMT; }

  inline void SetVetoPMTMaterialIndex(G4int val) { fVetoPMT = val; }
  inline G4int GetVetoPMTMaterialIndex() { return fVetoPMT; }

  inline void SetMuPMTMaterialIndex(G4int val) { fMuPMT = val; }
  inline G4int GetMuPMTMaterialIndex() { return fMuPMT; }

  inline void SetWriteDeposits(G4int val) { fWriteDeposits = val; }
  inline G4int GetWriteDeposits() { return fWriteDeposits; }

  inline void SetWriteThermalElectrons(G4int val) { fWriteThermalElectrons = val; }
  inline G4int GetWriteThermalElectrons() { return fWriteThermalElectrons; }

  inline void SetWriteDaughters(G4int val) { fWriteDaughters = val; }
  inline G4int GetWriteDaughters() { return fWriteDaughters; }

  inline void SetWriteFilter(G4int val) { fWriteFilter = val; }
  inline G4int GetWriteFilter() { return fWriteFilter; }

  inline void SetWriteFilterMinimumEnergy(G4double val) { fWriteFilterMinimumEnergy = val; }
  inline G4double GetWriteFilterMinimumEnergy() { return fWriteFilterMinimumEnergy; }

  inline void SetWriteFilterMaximumEnergy(G4double val) { fWriteFilterMaximumEnergy = val; }
  inline G4double GetWriteFilterMaximumEnergy() { return fWriteFilterMaximumEnergy; }

  inline void IncrArDMeventsShieldPE() { fArDMeventsShieldPE++; }
  inline G4int GetArDMeventsShieldPE() { return fArDMeventsShieldPE; }

  inline void SetRDMDecay(G4bool val) { fRDMDecay = val; }
  inline G4bool GetRDMDecay() { return fRDMDecay; }

  inline void SetReDStacking(G4bool val) { fReDStacking = val; }
  inline G4bool GetReDStacking() { return fReDStacking; }

  inline void SetRealPDGMeanLife(G4double val) { fRealPDGMeanLife = val; }
  inline G4double GetRealPDGMeanLife() { return fRealPDGMeanLife; }

  inline void SetPreAbsTime(G4double val) { fPreAbsTime = val; }
  inline G4double GetPreAbsTime() { return fPreAbsTime; }

  inline void SetNDaughters(G4int val) { fNDaughters = val; }
  inline G4int GetNDaughters() { return fNDaughters; }

  inline void SetRDMChain(G4bool val) { fRDMChain = val; }
  inline G4bool GetRDMChain() { return fRDMChain; }

  inline void SetKillS1S2(G4bool val) { fKillS1S2 = val; }
  inline G4bool GetKillS1S2() { return fKillS1S2; }

  inline void SetKillS2(G4bool val) { fKillS2 = val; }
  inline G4bool GetKillS2() { return fKillS2; }

  inline void SetKillS1(G4bool val) { fKillS1 = val; }
  inline G4bool GetKillS1() { return fKillS1; }

  inline void SetScaleS2(G4double val) { fScaleS2 = val; }
  inline G4double GetScaleS2() { return fScaleS2; }

  inline void SetTMBfraction(G4double val) { fTMBfraction = val; }
  inline G4double GetTMBfraction() { return fTMBfraction; }

  inline void SetTimeCut(G4double val) { fTimeCut = val; }
  inline G4double GetTimeCut() { return fTimeCut; }

  inline void SetTrackCounter(G4int val) { fTrackCounter = val; }
  inline G4int GetTrackCounter() { return fTrackCounter; }

  inline void SetEndOfEvent(G4bool val) { fEndOfEvent = val; }
  inline G4bool GetEndOfEvent() { return fEndOfEvent; }

  inline void SetDSLightTrackSecondaries(G4bool val) { fDSLightTrackSecondaries = val; }
  inline G4bool GetDSLightTrackSecondaries() { return fDSLightTrackSecondaries; }

  inline void SetLArGArBoundaryPosZ(G4double val) { fLArGArBoundaryPositionZ = val; }
  inline G4double GetLArGArBoundaryPosZ() { return fLArGArBoundaryPositionZ; }

  inline void SetGridSurfaceDistance(G4double val) { fGridSurfaceDistance = val; }
  inline G4double GetGridSurfaceDistance() { return fGridSurfaceDistance; }

  inline void SetDriftField(G4double val) { fDriftField = val; }
  inline G4double GetDriftField() { return fDriftField; }

  inline void SetExtractionField(G4double val) { fExtractionField = val; }
  inline G4double GetExtractionField() { return fExtractionField; }

  inline void SetThomasImelNullField(G4double val) { fThomasImelNullField = val; }
  inline G4double GetThomasImelNullField() { return fThomasImelNullField; }

  inline void SetThomasImelEp0(G4double val) { fThomasImelEp0 = val; }
  inline G4double GetThomasImelEp0() { return fThomasImelEp0; }

  inline void SetThomasImelEp1(G4double val) { fThomasImelEp1 = val; }
  inline G4double GetThomasImelEp1() { return fThomasImelEp1; }

  inline void SetDokeBirksNFp1(G4double val) { fDokeBirksNFp1 = val; }
  inline G4double GetDokeBirksNFp1() { return fDokeBirksNFp1; }

  inline void SetDokeBirksNFp3(G4double val) { fDokeBirksNFp3 = val; }
  inline G4double GetDokeBirksNFp3() { return fDokeBirksNFp3; }

  inline void SetDokeBirksEp1(G4double val) { fDokeBirksEp1 = val; }
  inline G4double GetDokeBirksEp1() { return fDokeBirksEp1; }

  inline void SetDokeBirksEp2(G4double val) { fDokeBirksEp2 = val; }
  inline G4double GetDokeBirksEp2() { return fDokeBirksEp2; }

  inline void SetDokeBirksEp3(G4double val) { fDokeBirksEp3 = val; }
  inline G4double GetDokeBirksEp3() { return fDokeBirksEp3; }

  G4double* GetITOParameters(G4double, G4double, G4double, G4double, G4double);
  G4double GetGridParameters(G4double);
  G4double GetTPBEfficiency(G4double, G4double);
  G4double GetPENEfficiency(G4double, G4double);

  inline void SetNumberOfHits(G4int val) { fNumberOfHits = val; }
  inline G4int GetNumberOfHits() { return fNumberOfHits; }

  inline void SetScintillator(G4int val) { fScintillator = val; }
  inline G4int GetScintillator() { return fScintillator; }

  inline void SetCTFfill(G4int val) { fCTFfill = val; }
  inline G4int GetCTFfill() { return fCTFfill; }

  inline void SetFastSimulation(G4int val) { fFastSimulation = val; }
  inline G4int GetFastSimulation() { return fFastSimulation; }

  inline void SetSourcePosition(G4ThreeVector val) { fSourcePosition = val; }
  inline G4ThreeVector GetSourcePosition() { return fSourcePosition; }

  inline void SetIsExternalLArScintillating(G4bool val) { fIsExternalLArScintillating = val; }
  inline G4bool GetIsExternalLArScintillating() { return fIsExternalLArScintillating; }

  inline void SetIsCopperRings(G4bool val) { fIsCopperRings = val; }
  inline G4bool GetIsCopperRings() { return fIsCopperRings; }

  inline void SetVetoYieldFactor(G4double val) { fVetoYieldFactor = val; }
  inline G4double GetVetoYieldFactor() { return fVetoYieldFactor; }

  inline void SetTunedS1At200V(G4bool val) { fTunedS1At200V = val; }
  inline G4bool GetTunedS1At200V() { return fTunedS1At200V; }

  inline void SetSourceHolderCenter(G4ThreeVector val) { fSourceHolderCenter = val; }
  inline G4ThreeVector GetSourceHolderCenter() { return fSourceHolderCenter; }

  inline void SetSourceHolderTheta(G4double val) { fSourceHolderTheta = val; }
  inline G4double GetSourceHolderTheta() { return fSourceHolderTheta; }

  inline void SetSourceHolderPhi(G4double val) { fSourceHolderPhi = val; }
  inline G4double GetSourceHolderPhi() { return fSourceHolderPhi; }

  inline void SetSourceHolderFlag(G4bool val) { fSourceHolderFlag = val; }
  inline G4bool GetSourceHolderFlag() { return fSourceHolderFlag; }

  inline void SetSourceHolderLeadFlag(G4bool val) { fSourceHolderLeadFlag = val; }
  inline G4bool GetSourceHolderLeadFlag() { return fSourceHolderLeadFlag; }

  inline void SetClusterWritingActivated(G4bool val) { fClusterWritingActivated = val; }
  inline G4bool GetClusterWritingActivated() { return fClusterWritingActivated; }

  inline void Set5KGeometry(G4bool val) { fIs5KGeometry = val; }
  inline G4bool Get5KGeometry() { return fIs5KGeometry; }

  inline void Set20KGeometry(G4bool val) { fIs20KGeometry = val; }
  inline G4bool Get20KGeometry() { return fIs20KGeometry; }

  inline void Set20KAlternateVeto(G4bool val) { fIs20KAlternateVeto = val; }
  inline G4bool Get20KAlternateVeto() { return fIs20KAlternateVeto; }

  inline void SetDartGeometry(G4bool val) { fIsDartGeometry = val; }
  inline G4bool GetDartGeometry() { return fIsDartGeometry; }

  inline void SetDartTeflonThickness(G4double val) { fDartTeflonThickness = val; }
  inline G4double GetDartTeflonThickness() { return fDartTeflonThickness; }
  inline void SetDartTeflonHeight(G4double val) { fDartTeflonHeight = val; }
  inline G4double GetDartTeflonHeight() { return fDartTeflonHeight; }
  inline void SetDartTeflonRadius(G4double val) { fDartTeflonRadius = val; }
  inline G4double GetDartTeflonRadius() { return fDartTeflonRadius; }

  inline void SetDartDewarThickness(G4double val) { fDartDewarThickness = val; }
  inline G4double GetDartDewarThickness() { return fDartDewarThickness; }
  inline void SetDartDewarHeight(G4double val) { fDartDewarHeight = val; }
  inline G4double GetDartDewarHeight() { return fDartDewarHeight; }
  inline void SetDartDewarRadius(G4double val) { fDartDewarRadius = val; }
  inline G4double GetDartDewarRadius() { return fDartDewarRadius; }

  inline void SetArDMTeflonHeight(G4double val) { fArDMTeflonHeight = val; }
  inline G4double GetArDMTeflonHeight() { return fArDMTeflonHeight; }
  inline void SetArDMTeflonRadius(G4double val) { fArDMTeflonRadius = val; }
  inline G4double GetArDMTeflonRadius() { return fArDMTeflonRadius; }

  inline void SetArDMTankHeight(G4double val) { fArDMTankHeight = val; }
  inline G4double GetArDMTankHeight() { return fArDMTankHeight; }
  inline void SetArDMTankRadius(G4double val) { fArDMTankRadius = val; }
  inline G4double GetArDMTankRadius() { return fArDMTankRadius; }
  // vpesudo
  inline void SetArDMTotalTankHeight(G4double val) { fArDMTotalTankHeight = val; }
  inline G4double GetArDMTotalTankHeight() { return fArDMTotalTankHeight; }

  inline void SetArDMShieldHeight(G4double val) { fArDMShieldHeight = val; }
  inline G4double GetArDMShieldHeight() { return fArDMShieldHeight; }
  inline void SetArDMShieldRadius(G4double val) { fArDMShieldRadius = val; }
  inline G4double GetArDMShieldRadius() { return fArDMShieldRadius; }

  inline void SetArDMLeadHeight(G4double val) { fArDMLeadHeight = val; }
  inline G4double GetArDMLeadHeight() { return fArDMLeadHeight; }

  inline void SetArDMShieldLeadRadius(G4double val) { fArDMShieldLeadRadius = val; }
  inline G4double GetArDMShieldLeadRadius() { return fArDMShieldLeadRadius; }

  inline void SetArDMTestArgonHeight(G4double val) { fArDMTestArgonHeight = val; }
  inline G4double GetArDMTestArgonHeight() { return fArDMTestArgonHeight; }

  inline void SetArDMTestArgonRadius(G4double val) { fArDMTestArgonRadius = val; }
  inline G4double GetArDMTestArgonRadius() { return fArDMTestArgonRadius; }

  inline void SetArDMGeometry(G4bool val) { fIsArDMGeometry = val; }
  inline G4bool GetArDMGeometry() { return fIsArDMGeometry; }

  inline void SetReDGeometry(G4bool val) { fIsReDGeometry = val; }
  inline G4bool GetReDGeometry() { return fIsReDGeometry; }

  inline void SetPETGeometry(G4bool val) { fIsPETGeometry = val; }
  inline G4bool GetPETGeometry() { return fIsPETGeometry; }

  inline void SetIsDS20kLSV(G4bool val) { fIsDS20kLSV = val; }
  inline G4bool GetIsDS20kLSV() const { return fIsDS20kLSV; }

  inline void SetIsCherenkov(G4bool val) { fIsCherenkov = val; }
  inline G4bool GetIsCherenkov() { return fIsCherenkov; }

  inline void SetIsEmittedCherenkov(G4bool val) { fIsEmittedCherenkov = val; }
  inline G4bool GetIsEmittedCherenkov() { return fIsEmittedCherenkov; }

  inline void SetIsCustomCryostat(G4bool val) { fIsCustomCryostat = val; }
  inline G4bool GetCustomCryostat() { return fIsCustomCryostat; }

  inline void SetDS20kCryoWallThick(G4double val) { fDS20kCryoWallThick = val; }
  inline G4double GetDS20kCryoWallThick() { return fDS20kCryoWallThick; }

  inline void SetDS20kCryoCapWallThick(G4double val) { fDS20kCryoCapWallThick = val; }
  inline G4double GetDS20kCryoCapWallThick() { return fDS20kCryoCapWallThick; }
  
  inline void SetDS20kCryoCornerDistance(G4double val) { fDS20kCryoCornerDistance = val; }
  inline G4double GetDS20kCryoCornerDistance() { return fDS20kCryoCornerDistance; }

  inline void SetDS20kTPCheight(G4double val) { fDS20kTPCheight = val; }
  inline G4double GetDS20kTPCheight() { return fDS20kTPCheight; }

  inline void SetDS20kTPCedge(G4double val) { fDS20kTPCedge = val; }
  inline G4double GetDS20kTPCedge() { return fDS20kTPCedge; }

  inline void SetDS20kHDPEShellThickness(G4double val) { fDS20kHDPEShellThickness = val; };
  inline G4double GetDS20kHDPEShellThickness() { return fDS20kHDPEShellThickness; };

  inline void SetDS20kHDPEShellCapThickness(G4double val) { fDS20kHDPEShellCapThickness = val; };
  inline G4double GetDS20kHDPEShellCapThickness() { return fDS20kHDPEShellCapThickness; };

  inline void SetDS20kCryoMaterial(G4int val) { fDS20kCryoMaterial = val; }
  inline G4int GetDS20kCryoMaterial() { return fDS20kCryoMaterial; }

  inline void SetDS20kWindowMaterial(G4int val) { fDS20kWindowMaterial = val; }
  inline G4int GetDS20kWindowMaterial() { return fDS20kWindowMaterial; }

  inline void SetDS20kRemoveBondingRing(G4int val) { fDS20kRemoveBondingRing = val; }
  inline G4int GetDS20kRemoveBondingRing() { return fDS20kRemoveBondingRing; }

  inline void SetDS20kLArThicknessAboveTPC(G4int val) { fDS20kLArThicknessAboveTPC = val; }
  inline G4int GetDS20kLArThicknessAboveTPC() { return fDS20kLArThicknessAboveTPC; }

  inline void SetSiPMOffset(G4double val) { fSiPMOffset = val; };
  inline G4double GetSiPMOffset() { return fSiPMOffset; };

  inline void SetDS20kTPCVesselThickness(G4double val) { fDS20kTPCVesselThickness = val; };
  inline G4double GetDS20kTPCVesselThickness() { return fDS20kTPCVesselThickness; };

  inline void SetDS20kLArBufferThickness(G4double val) { fDS20kLArBufferThickness = val; };
  inline G4double GetDS20kLArBufferThickness() { return fDS20kLArBufferThickness; };

  inline void SetDS20kGdPlasticThickness(G4double val) { fDS20kGdPlasticThickness = val; };
  inline G4double GetDS20kGdPlasticThickness() { return fDS20kGdPlasticThickness; };

  inline void SetDS20kGasPocketThickness(G4double val) { fDS20kGasPocketThickness = val; };

  inline G4double GetDS20kGasPocketThickness() { return fDS20kGasPocketThickness; };

  inline void SetDS20kFacePanels(G4bool val) { fDS20kFacePanels = val; };
  inline G4bool GetDS20kFacePanels() { return fDS20kFacePanels; }

  inline void SetDS20kEdgePanels(G4bool val) { fDS20kEdgePanels = val; };
  inline G4bool GetDS20kEdgePanels() { return fDS20kEdgePanels; }

  inline void SetDS20kSiPMs(G4bool val) { fDS20kSiPMs = val; };
  inline G4bool GetDS20kSiPMs() { return fDS20kSiPMs; }

  inline void SetDS20kSiPMSide(G4double val) { fDS20kSiPMSide = val; };
  inline G4double GetDS20kSiPMSide() { return fDS20kSiPMSide; };

  inline void SetDS20kSiPMHeight(G4double val) { fDS20kSiPMHeight = val; };
  inline G4double GetDS20kSiPMHeight() { return fDS20kSiPMHeight; };

  inline void SetDS20knSiPMs(G4double val) { fDS20knSiPMs = val; };
  inline G4double GetDS20knSiPMs() { return fDS20knSiPMs; };

  inline void SetDS20knSiPMInside(G4double val) { fDS20knSiPMInside = val; };
  inline G4double GetDS20knSiPMInside() { return fDS20knSiPMInside; };

  inline void SetDS20knSiPMOutside(G4double val) { fDS20knSiPMOutside = val; };
  inline G4double GetDS20knSiPMOutside() { return fDS20knSiPMOutside; };

  inline void SetDS20kSiPMsAutoplacement(G4bool val) { fDS20kSiPMsAutoplacement = val; };
  inline G4bool GetDS20kSiPMsAutoplacement() { return fDS20kSiPMsAutoplacement; }

  inline void SetDS20kSiPMsAutoplacementFilename(G4String val) { fDS20kSiPMsAutoplacementFilename = val; };
  inline G4String GetDS20kSiPMsAutoplacementFilename() { return fDS20kSiPMsAutoplacementFilename; };

  inline void SetDS20kepsilonSiPM(G4double val) { fDS20kepsilonSiPM = val; };
  inline G4double GetDS20kepsilonSiPM() { return fDS20kepsilonSiPM; };

  inline void SetDS20kSiPMFractionInside(G4double val) { fDS20kSiPMFractionInside = val; };
  inline G4double GetDS20kSiPMFractionInside() { return fDS20kSiPMFractionInside; }

  inline void SetDS20kSiPMFractionOutside(G4double val) { fDS20kSiPMFractionOutside = val; };
  inline G4double GetDS20kSiPMFractionOutside() { return fDS20kSiPMFractionOutside; }

  inline void SetDS20kTPBOnSiPM(G4bool val) { fDS20kTPBOnSiPM = val; };
  inline G4bool GetDS20kTPBOnSiPM() { return fDS20kTPBOnSiPM; }

  inline void SetDS20kSiPMUniformInsideSides(G4bool val) { fDS20kSiPMUniformInsideSides = val; };
  inline G4bool GetDS20kSiPMUniformInsideSides() { return fDS20kSiPMUniformInsideSides; }

  inline void SetDS20kSiPMUniformInsideCaps(G4bool val) { fDS20kSiPMUniformInsideCaps = val; };
  inline G4bool GetDS20kSiPMUniformInsideCaps() { return fDS20kSiPMUniformInsideCaps; }

  inline void SetDS20kSiPMUniformOutsideSides(G4bool val) { fDS20kSiPMUniformOutsideSides = val; };
  inline G4bool GetDS20kSiPMUniformOutsideSides() { return fDS20kSiPMUniformOutsideSides; }

  inline void SetDS20kSiPMUniformOutsideCaps(G4bool val) { fDS20kSiPMUniformOutsideCaps = val; };
  inline G4bool GetDS20kSiPMUniformOutsideCaps() { return fDS20kSiPMUniformOutsideCaps; }

  inline void SetDS20kbestLightCollection(G4bool val) { fDS20kbestLightCollection = val; };
  inline G4bool GetDS20kbestLightCollection() { return fDS20kbestLightCollection; }

  inline void SetDS20kVerticalInnerTPBPanely(G4double val) { fDS20kVerticalInnerTPBPanely = val; };
  inline G4double GetDS20kVerticalInnerTPBPanely() { return fDS20kVerticalInnerTPBPanely; }

  inline void SetDS20kHorizontalInnerTPBPanely(G4double val) { fDS20kHorizontalInnerTPBPanely = val; };
  inline G4double GetDS20kHorizontalInnerTPBPanely() { return fDS20kHorizontalInnerTPBPanely; }

  inline void SetDS20kVerticalOuterTPBPanely(G4double val) { fDS20kVerticalOuterTPBPanely = val; };
  inline G4double GetDS20kVerticalOuterTPBPanely() { return fDS20kVerticalOuterTPBPanely; }

  inline void SetDS20kHorizontalOuterTPBPanely(G4double val) { fDS20kHorizontalOuterTPBPanely = val; };
  inline G4double GetDS20kHorizontalOuterTPBPanely() { return fDS20kHorizontalOuterTPBPanely; }

  inline void SetDS20kinternalOffset(G4double val) { fDS20kinternalOffset = val; };
  inline G4double GetDS20kinternalOffset() { return fDS20kinternalOffset; }

  inline void SetDS20kexternalOffset(G4double val) { fDS20kexternalOffset = val; };
  inline G4double GetDS20kexternalOffset() { return fDS20kexternalOffset; }

  inline void SetDS20kcenterInnerOffset(G4double val) { fDS20kcenterInnerOffset = val; };
  inline G4double GetDS20kcenterInnerOffset() { return fDS20kcenterInnerOffset; }

  inline void SetDS20kcenterOuterOffset(G4double val) { fDS20kcenterOuterOffset = val; };
  inline G4double GetDS20kcenterOuterOffset() { return fDS20kcenterOuterOffset; }

  inline void SetDS20kupperHorizontalInnerPanelOffset(G4double val) { fDS20kupperHorizontalInnerPanelOffset = val; };
  inline G4double GetDS20kupperHorizontalInnerPanelOffset() { return fDS20kupperHorizontalInnerPanelOffset; }

  inline void SetDS20klowerHorizontalInnerPanelOffset(G4double val) { fDS20klowerHorizontalInnerPanelOffset = val; };
  inline G4double GetDS20klowerHorizontalInnerPanelOffset() { return fDS20klowerHorizontalInnerPanelOffset; }

  inline void SetDS20kupperHorizontalOuterPanelOffset(G4double val) { fDS20kupperHorizontalOuterPanelOffset = val; };
  inline G4double GetDS20kupperHorizontalOuterPanelOffset() { return fDS20kupperHorizontalOuterPanelOffset; }

  inline void SetDS20klowerHorizontalOuterPanelOffset(G4double val) { fDS20klowerHorizontalOuterPanelOffset = val; };
  inline G4double GetDS20klowerHorizontalOuterPanelOffset() { return fDS20klowerHorizontalOuterPanelOffset; }

  inline void SetDS20kverticalInnerPanelOffset(G4double val) { fDS20kverticalInnerPanelOffset = val; };
  inline G4double GetDS20kverticalInnerPanelOffset() { return fDS20kverticalInnerPanelOffset; }

  inline void SetDS20kverticalOuterPanelOffset(G4double val) { fDS20kverticalOuterPanelOffset = val; };
  inline G4double GetDS20kverticalOuterPanelOffset() { return fDS20kverticalOuterPanelOffset; }

  inline void SetDS20kplasticLArBufferThickness(G4double val) { fDS20kplasticLArBufferThickness = val; };
  inline G4double GetDS20kplasticLArBufferThickness() { return fDS20kplasticLArBufferThickness; }

  inline void SetDS20kCopperThickness(G4double val) { fDS20kCopperThickness = val; };
  inline G4double GetDS20kCopperThickness() { return fDS20kCopperThickness; }

  inline void SetDS20kinnerCylinderOffset(G4double val) { fDS20kinnerCylinderOffset = val; };
  inline G4double GetDS20kinnerCylinderOffset() { return fDS20kinnerCylinderOffset; }

  inline void SetDS20kouterCylinderOffset(G4double val) { fDS20kouterCylinderOffset = val; };
  inline G4double GetDS20kouterCylinderOffset() { return fDS20kouterCylinderOffset; }

  inline void SetDS20kinnerTPBCylinderThickness(G4double val) { fDS20kinnerTPBCylinderThickness = val; };
  inline G4double GetDS20kinnerTPBCylinderThickness() { return fDS20kinnerTPBCylinderThickness; }

  inline void SetDS20kouterTPBCylinderThickness(G4double val) { fDS20kouterTPBCylinderThickness = val; };
  inline G4double GetDS20kouterTPBCylinderThickness() { return fDS20kouterTPBCylinderThickness; }

  inline void SetDS20kTPBThickness(G4double val) { fDS20kTPBThickness = val; };
  inline G4double GetDS20kTPBThickness() { return fDS20kTPBThickness; }

  inline void SetDS20kTPBLayers(G4String val) { fDS20kTPBLayers = val; };
  inline G4String GetDS20kTPBLayers() { return fDS20kTPBLayers; };

  inline G4double GetDS20kTPCverticalShift() { return fDS20kTPCverticalShift; }

  inline void SetGdLayerThickness(G4double val) { fGdLayerThickness = val; }
  inline G4double GetGdLayerThickness() { return fGdLayerThickness; }

  inline void SetWTankR(G4int val) { fG3WTankR = val; }
  inline G4double GetWTankR() { return fG3WTankR; }

  inline void SetAcrylicVesselDiameter(G4int val) { fAcrylicVesselDiameter = val; }
  inline G4double GetAcrylicVesselDiameter() { return fAcrylicVesselDiameter; }

  inline void SetIDiameter(G4int val) { fG3IDiameter = val; }
  inline G4double GetIDiameter() { return fG3IDiameter; }

  inline void SetWTankHeight(G4int val) { fG3WTankHeight = val; }
  inline G4double GetWTankHeight() { return fG3WTankHeight; }

  inline void SetLicNDNuclRecE(G4double val) {
    fLicNDNuclRecE = val;
    SetLicNdAngleByEnergy();
  };
  inline G4double GetLicNDNuclRecE() { return fLicNDNuclRecE; }

  inline void SetLicNDTheta(G4double val) {
    fLicNDTheta = val;
    SetLicNdAngleByAng();
  };
  inline G4double GetLicNDTheta() { return fLicNDTheta; }

  inline void SetLicNDPhi1(G4double val) { fLicNDPhi1 = val; };
  inline G4double GetLicNDPhi1() { return fLicNDPhi1; }

  inline void SetLicNDPhi2(G4double val) { fLicNDPhi2 = val; };
  inline G4double GetLicNDPhi2() { return fLicNDPhi2; }

  inline void SetLicNDRadius(G4double val) { fLicNDRadius = val; };
  inline G4double GetLicNDRadius() { return fLicNDRadius; }

  inline G4bool GetLicNdAngleByEnergy() { return fIsLicNdAngleByEnergy; };
  inline G4bool GetLicNdAngleByAng() { return fIsLicNdAngleByAng; };

  inline void SetLicNeutronE(G4double val) { fLicNeutronE = val; }
  inline G4double GetLicNeutronE() { return fLicNeutronE; };

  inline void SetLicWallActivation(G4bool val) { fLicWallActivated = val; };
  inline G4bool GetLicWallActivation() { return fLicWallActivated; };

  inline void SetLicWallThick(G4double val) { fLicWallThick = val; };
  inline G4double GetLicWallThick() { return fLicWallThick; };

  inline void SetLicWallWindowLength(G4double val) { fLicWallWindowLength = val; };
  inline G4double GetLicWallWindowLength() { return fLicWallWindowLength; };

  inline void SetLicWallDistanceToTPC(G4double val) { fLicWallDistanceToTPC = val; };
  inline G4double GetLicWallDistanceToTPC() { return fLicWallDistanceToTPC; };

  inline void SetLicExpHall(G4double val) { fLicExpHall = val; };
  inline G4double GetLicExpHall() { return fLicExpHall; };

  inline void SetAriserGeometry(G4bool val) { fIsAriserGeometry = val; }
  inline G4bool GetAriserGeometry() { return fIsAriserGeometry; }

  inline void SetAriserSourceDistanceToCryo(G4double val) { fAriserSourceDistanceToCryo = val; };
  inline G4double GetAriserSourceDistanceToCryo() { return fAriserSourceDistanceToCryo; };

  inline void SetAriserBEGeDistanceToCryo(G4double val) { fAriserBEGeDistanceToCryo = val; };
  inline G4double GetAriserBEGeDistanceToCryo() { return fAriserBEGeDistanceToCryo; };

  inline void SetAriserBEGeAngle(G4double val) { fAriserBEGeAngle = val; };
  inline G4double GetAriserBEGeAngle() { return fAriserBEGeAngle; };

  inline void SetAriserSourcePosition(G4ThreeVector val) { fAriserSourcePosition = val; };
  inline G4ThreeVector GetAriserSourcePosition() { return fAriserSourcePosition; };

  inline void SetFLUKAfilename(G4String val) { fFLUKAfilename = val; };
  inline G4String GetFLUKAfilename() { return fFLUKAfilename; };

  inline void SetReDConfiguration(G4int val) { fReDConfiguration = val; }
  inline G4int GetReDConfiguration() const { return fReDConfiguration; }

  inline void SetReDNDSurveyMode(NDSurveyMode val) { fReDNDSurveyMode = val; };
  inline NDSurveyMode GetReDNDSurveyMode() { return fReDNDSurveyMode; };

  inline void SetReDNDSurveyFileName(G4String val) { fReDNDSurveyFilename = val; };
  inline G4String GetReDNDSurveyFileName() { return fReDNDSurveyFilename; };

  inline void SetReDNDPositioningMode(NDPosMode val) { fReDNDPositioningMode = val; };
  inline NDPosMode GetReDNDPositioningMode() { return fReDNDPositioningMode; };

  inline void SetTargetToTPCDistance(G4double val) { fTargetTPCDistance = val; };
  inline G4double GetTargetToTPCDistance() { return fTargetTPCDistance; };

  inline void SetBeamToTPCAngleThetaX(G4double val) { fThetaArX = val; };
  inline G4double GetBeamToTPCAngleThetaX() { return fThetaArX; };

  inline void SetBeamToTPCAnglePhiX(G4double val) { fPhiArX = val; };
  inline G4double GetBeamToTPCAnglePhiX() { return fPhiArX; };

  inline void SetNeutronBeamThetaX(G4double val) { fNeutronAngleThetaX = val; };
  inline G4double GetNeutronBeamThetaX() { return fNeutronAngleThetaX; };

  inline void SetNeutronBeamPhiX(G4double val) { fNeutronAnglePhiX = val; };
  inline G4double GetNeutronBeamPhiX() { return fNeutronAnglePhiX; };

  inline void SetIsProtoProtoBottomSiPM(G4bool val) { fIsProtoProtoBottomSiPM = val; };
  inline G4bool GetIsProtoProtoBottomSiPM() { return fIsProtoProtoBottomSiPM; };

  /*  inline void         SetGantryConfiguration(G4int val)      {
    fGantryConfiguration = val; }; inline G4int        GetGantryConfiguration()
    { return fGantryConfiguration; }; inline void         SetGantrySource(G4int
    val)             { fGantryConfiguration = val; }; inline G4int
    GetGantrySource()                      { return fGantryConfiguration; };
    inline void SetGantryRelativeSourcePosition(G4ThreeVector val)
    {fGantryRelativeSourcePos = val;}; inline G4ThreeVector
    GetGantryRelativeSourcePosition()            {return
    fGantryRelativeSourcePos;}; inline void
    SetGantryGlobalPosition(G4ThreeVector val) {fGantryGlobalPosition=val;};
    inline G4ThreeVector GetGantryGlobalPosition(){return
    fGantryGlobalPosition;}; inline void SetGantrySourceMaterial(G4String val)
    {fGantrySourceMat=val;}; inline G4String GetGantrySourceMaterial(){return
    fGantrySourceMat;}; inline G4double GetGantrySourceRadius(){return
    fGantrySourceRadius;}; inline void SetGantrySourceRadius(G4double
    val){fGantrySourceRadius=val;};
*/
  // for calibration gantry
  inline void SetGantryConfiguration(G4bool val) { fGantryConfiguration = val; };
  inline G4bool GetGantryConfiguration() { return fGantryConfiguration; };
  inline void SetGantryGlobalPosition(G4ThreeVector val) { fGantryGlobalPosition = val; };
  inline G4ThreeVector GetGantryGlobalPosition() { return fGantryGlobalPosition; };
  inline void SetGantryPipeID(G4double val) { fGantryPipeID = val; };
  inline G4double GetGantryPipeID() { return fGantryPipeID; };
  inline void SetGantryPipeWall(G4double val) { fGantryPipeWall = val; };
  inline G4double GetGantryPipeWall() { return fGantryPipeWall; };
  inline void SetGantryPipeBendingR(G4double val) { fGantryPipeBendingR = val; };
  inline G4double GetGantryPipeBendingR() { return fGantryPipeBendingR; };
  inline void SetGantryDistanceTPCside(G4double val) { fGantryDistanceTPCside = val; };
  inline G4double GetGantryDistanceTPCside() { return fGantryDistanceTPCside; };
  inline void SetGantryDistanceTPCbottom(G4double val) { fGantryDistanceTPCbottom = val; };
  inline G4double GetGantryDistanceTPCbottom() { return fGantryDistanceTPCbottom; };
  inline void SetGantryDistanceBtwPipes(G4double val) { fGantryDistanceBtwPipes = val; };
  inline G4double GetGantryDistanceBtwPipes() { return fGantryDistanceBtwPipes; };
  inline void SetGantrySurface(G4bool val) { fGantrySurface = val; };
  inline G4bool GetGantrySurface() { return fGantrySurface; };
  inline void SetGantrySurfaceOption(G4int val) { fGantrySurfaceOption = val; };
  inline G4int GetGantrySurfaceOption() { return fGantrySurfaceOption; };
  inline void SetGantryShield(G4bool val) { fGantryShield = val; };
  inline G4bool GetGantryShield() { return fGantryShield; };
  inline void SetGantryShieldThickness(G4double val) { fGantryShieldThickness = val; };
  inline G4double GetGantryShieldThickness() { return fGantryShieldThickness; };
  inline void SetGantryShieldHeight(G4double val) { fGantryShieldHeight = val; };
  inline G4double GetGantryShieldHeight() { return fGantryShieldHeight; };
  inline void SetGantryShieldOffset(G4double val) { fGantryShieldOffset = val; };
  inline G4double GetGantryShieldOffset() { return fGantryShieldOffset; };
  inline void SetGantryShieldMaterial(G4int val) { fGantryShieldMaterial = val; };
  inline G4int GetGantryShieldMaterial() { return fGantryShieldMaterial; };

  inline void SetDS20kWLSLArThickness(G4double val) { fDS20kWLSLArThickness = val; };
  inline G4double GetDS20kWLSLArThickness() { return fDS20kWLSLArThickness; }
  inline void SetDS20kWLSPENThickness(G4double val) { fDS20kWLSPENThickness = val; };
  inline G4double GetDS20kWLSPENThickness() { return fDS20kWLSPENThickness; }

  // plan b cryostat
  inline void SetDS20kCryoR(G4double val) { fDS20kCryoR = val; };
  inline G4double GetDS20kCryoR() { return fDS20kCryoR; }
  inline void SetDS20kCryoH(G4double val) { fDS20kCryoH = val; };
  inline G4double GetDS20kCryoH() { return fDS20kCryoH; }
  inline void SetDS20kCryoVacuumTh(G4double val) { fDS20kCryoVacuumTh = val; };
  inline G4double GetDS20kCryoVacuumTh() { return fDS20kCryoVacuumTh; }
  inline void SetDS20kCryoBottomCap(G4double val) { fDS20kCryoBottomCap = val; };
  inline G4double GetDS20kCryoBottomCap() { return fDS20kCryoBottomCap; }
  inline void SetDS20kCryoTopCap(G4double val) { fDS20kCryoTopCap = val; };
  inline G4double GetDS20kCryoTopCap() { return fDS20kCryoTopCap; }
  inline void SetDS20kCryoTopOffset(G4double val) { fDS20kCryoTopOffset = val; };
  inline G4double GetDS20kCryoTopOffset() { return fDS20kCryoTopOffset; }
  inline void SetDS20kCryoBottomOffset(G4double val) { fDS20kCryoBottomOffset = val; };
  inline G4double GetDS20kCryoBottomOffset() { return fDS20kCryoBottomOffset; }
  inline void SetDS20kCryoBottomCapOut(G4double val) { fDS20kCryoBottomCapOut = val; };
  inline G4double GetDS20kCryoBottomCapOut() { return fDS20kCryoBottomCapOut; }
  inline void SetDS20kCryoTopCapOut(G4double val) { fDS20kCryoTopCapOut = val; };
  inline G4double GetDS20kCryoTopCapOut() { return fDS20kCryoTopCapOut; }
  inline void SetDS20kCryoTopOffsetOut(G4double val) { fDS20kCryoTopOffsetOut = val; };
  inline G4double GetDS20kCryoTopOffsetOut() { return fDS20kCryoTopOffsetOut; }
  inline void SetDS20kCryoBottomOffsetOut(G4double val) { fDS20kCryoBottomOffsetOut = val; };
  inline G4double GetDS20kCryoBottomOffsetOut() { return fDS20kCryoBottomOffsetOut; }
  inline void SetDS20KCryoBarrelH(G4double val) { fDS20KCryoBarrelH = val; };
  inline G4double GetDS20KCryoBarrelH() { return fDS20KCryoBarrelH; }

  inline void SetDS20kTPCHeight_top(G4double val) { fDS20kTPCHeight_top = val; };
  inline G4double GetDS20kTPCHeight_top() { return fDS20kTPCHeight_top; }

  inline void SetDS20kTPCHeight_bottom(G4double val) { fDS20kTPCHeight_bottom = val; };
  inline G4double GetDS20kTPCHeight_bottom() { return fDS20kTPCHeight_bottom; }

  inline void SetDS20kTPCInscribedR(G4double val) { fDS20kTPCInscribedR = val; };
  inline G4double GetDS20kTPCInscribedR() { return fDS20kTPCInscribedR; }

  inline void SetDS20kSiPMPosVector(vector<G4ThreeVector> val) { fDS20kSiPMPosVector = val; };
  inline vector<G4ThreeVector> GetDS20kSiPMPosVector() { return fDS20kSiPMPosVector; }

  inline void SetDS20kGdConcentration(G4double val) { fDS20kGdConcentration = val; };
  inline G4double GetDS20kGdConcentration() { return fDS20kGdConcentration; }

  inline void SetOverwriteCounter(G4bool val) { fOverwriteCounter = val; };
  inline G4double GetOverwriteCounter() { return fOverwriteCounter; }

  inline void SetIsDecayInVacuum(G4bool val) { fOverwriteCounter = val; };
  inline G4bool IsDecayInVacuum() { return fOverwriteCounter; }

  inline void SetSkipEventsWithNoDaughters(G4bool val) {fSkipEventsWithNoDaughters = val;}
  inline G4bool SkipEventsWithNoDaughters() { return fSkipEventsWithNoDaughters; }

  inline void SetAbortRun(G4bool val) {fAbortRun = val;}
  inline G4bool IsAbortRun() { return fAbortRun; }

  inline void SetMatScintillation(G4bool val) {fMatScintillation = val;}
  inline G4bool GetMatScintillation() { return fMatScintillation; }

  inline void SetPurePMMATPC(G4bool val) {fPurePMMATPC = val;}
  inline G4bool GetPurePMMATPC() { return fPurePMMATPC; }

  inline void SetHybridTPC(G4bool val) {fHybridTPC = val;}
  inline G4bool GetHybridTPC() { return fHybridTPC; }

  inline void SetMultiplicator(G4double val) { fMultiplicator = val; };
  inline G4double GetMultiplicator() { return fMultiplicator; }


 private:
  inline void SetLicNdAngleByEnergy() {
    fIsLicNdAngleByAng = false;
    fIsLicNdAngleByEnergy = true;
  };

  inline void SetLicNdAngleByAng() {
    fIsLicNdAngleByAng = true;
    fIsLicNdAngleByEnergy = false;
  };

 private:
  G4bool fIsEnDepGenerator;
  G4int fOverlap;
  G4int fEventCounter;
  G4bool fExportGDML;
  G4bool fWritePhotons;
  G4int fVerbosity;
  G4int fBoronScintillatorIndex;
  G4int fLiquidArgonIndex;
  G4int fGaseousArgonIndex;
  G4int fLArAboveGridIndex;
  G4int fWriteThermalElectrons;
  G4int fWriteDeposits;
  G4int fWriteDaughters;
  G4int fWriteFilter;
  G4int fArDMeventsShieldPE;
  G4int fPMT;
  G4int fVetoPMT;
  G4int fMuPMT;
  G4bool fRDMDecay;
  G4bool fReDStacking;
  G4double fRealPDGMeanLife;
  G4double fPreAbsTime;
  G4int fNDaughters;
  G4bool fRDMChain;
  G4bool fKillS1S2;
  G4bool fKillS2;
  G4bool fKillS1;
  G4double fScaleS2;
  G4double fITOParameters[3];
  G4double fGridParameter;
  vector<double> fAng;
  vector<double> fTrans;
  G4double tpb_eff[100][100];
  G4double fTMBfraction;
  G4double fTimeCut;
  G4double pen_eff[100][100];

  G4double fWL[11];
  G4double fRI[11];
  G4double fEC[11];
  G4int fTrackCounter;
  G4bool fEndOfEvent;
  G4bool fDSLightTrackSecondaries;
  G4double fLArGArBoundaryPositionZ;  // z coordinate of the LAr-GAr interface -
  G4double fGridSurfaceDistance;  // z coordinate of the LAr-GAr interface -
                                      // detector configuration dependent
  G4double fDriftField;               // Electric field (V/cm)
  G4double fExtractionField;          // Electric field (V/cm)

  G4double fThomasImelNullField;  // Thomas Imel parameter at null field
  G4double fThomasImelEp0;        // Thomas Imel parameter at non-null field = p0 * E**p1
  G4double fThomasImelEp1;

  G4double fDokeBirksNFp1;
  G4double fDokeBirksNFp3;
  G4double fDokeBirksEp1;
  G4double fDokeBirksEp2;
  G4double fDokeBirksEp3;

  G4int fNumberOfHits;
  G4int fScintillator;
  G4int fFastSimulation;

  G4ThreeVector fSourcePosition;
  G4bool fIsExternalLArScintillating;
  G4bool fIsCopperRings;
  G4double fVetoYieldFactor;

  G4bool fTunedS1At200V;

  G4ThreeVector fSourceHolderCenter;
  G4double fSourceHolderTheta;
  G4double fSourceHolderPhi;
  G4bool fSourceHolderFlag;
  G4bool fSourceHolderLeadFlag;

  G4bool fClusterWritingActivated;
  G4int fCTFfill;
  G4bool fIs5KGeometry;
  G4bool fIs20KGeometry;
  G4bool fIs20KAlternateVeto;
  G4bool fIsDartGeometry;
  G4bool fIsArDMGeometry;
  G4bool fIsReDGeometry;
  G4bool fIsPETGeometry;
  G4bool fIsDS20kLSV;
  G4double fDS20kTPCverticalShift;
  G4double fDS20kTPCVesselThickness;
  G4bool fIsCherenkov;
  G4bool fIsEmittedCherenkov;

  G4double fDartTeflonThickness;
  G4double fDartTeflonHeight;
  G4double fDartTeflonRadius;

  G4double fDartDewarThickness;
  G4double fDartDewarHeight;
  G4double fDartDewarRadius;

  G4double fArDMTeflonHeight;
  G4double fArDMTeflonRadius;
  G4double fArDMTankHeight;
  G4double fArDMTankRadius;
  G4double fArDMTestArgonHeight;
  G4double fArDMTestArgonRadius;
  G4double fArDMShieldHeight;
  G4double fArDMShieldRadius;

  G4double fArDMLeadHeight;
  G4double fArDMShieldLeadRadius;

  G4double fArDMTotalTankHeight;

  G4bool fIsCustomCryostat;
  G4double fDS20kCryoWallThick;
  G4double fDS20kCryoCapWallThick;
  G4double fDS20kCryoCornerDistance;
  G4double fDS20kTPCheight;
  G4double fDS20kTPCedge;
  G4int fDS20kCryoMaterial;
  G4int fDS20kWindowMaterial;
  G4int fDS20kRemoveBondingRing;
  G4double fSiPMOffset;
  G4int fDS20kReflectorAlternate;

  G4double fDS20kLArBufferThickness;
  G4double fDS20kGdPlasticThickness;
  G4double fDS20kBufferPercentageCoverageInside;
  G4double fDS20kBufferPercentageCoverageOutside;
  G4double fDS20kTPBThickness;
  G4String fDS20kTPBLayers;
  G4double fDS20kGasPocketThickness;

  G4bool fDS20kFacePanels;
  G4bool fDS20kEdgePanels;

  G4bool fDS20kSiPMs;

  G4double fDS20kSiPMSide;
  G4double fDS20kSiPMHeight;

  G4double fDS20knSiPMInside;
  G4double fDS20knSiPMOutside;

  G4double fDS20knSiPMs;

  G4bool fDS20kSiPMsAutoplacement;
  G4String fDS20kSiPMsAutoplacementFilename;

  G4double fDS20kepsilonSiPM;

  G4double fDS20kSiPMFractionInside;
  G4double fDS20kSiPMFractionOutside;

  G4bool fDS20kTPBOnSiPM;

  G4bool fDS20kSiPMUniformInsideSides;
  G4bool fDS20kSiPMUniformInsideCaps;
  G4bool fDS20kSiPMUniformOutsideSides;
  G4bool fDS20kSiPMUniformOutsideCaps;

  G4bool fDS20kbestLightCollection;

  G4double fDS20kVerticalInnerTPBPanely;
  G4double fDS20kHorizontalInnerTPBPanely;
  G4double fDS20kVerticalOuterTPBPanely;
  G4double fDS20kHorizontalOuterTPBPanely;

  G4double fDS20kinternalOffset;
  G4double fDS20kexternalOffset;
  G4double fDS20kcenterInnerOffset;
  G4double fDS20kcenterOuterOffset;
  G4double fDS20kupperHorizontalInnerPanelOffset;
  G4double fDS20klowerHorizontalInnerPanelOffset;
  G4double fDS20kupperHorizontalOuterPanelOffset;
  G4double fDS20klowerHorizontalOuterPanelOffset;

  G4double fDS20kverticalInnerPanelOffset;
  G4double fDS20kverticalOuterPanelOffset;

  G4double fDS20kplasticLArBufferThickness;

  G4double fDS20kCopperThickness;

  G4double fDS20kinnerCylinderOffset;
  G4double fDS20kouterCylinderOffset;
  G4double fDS20kinnerTPBCylinderThickness;
  // G4double fDS20kinnerAcrylicCylinderThickness;
  G4double fDS20kouterTPBCylinderThickness;
  // G4double fDS20kouterAcrylicCylinderThickness;

  G4double fDS20kHDPEShellThickness;
  G4double fDS20kHDPEShellCapThickness; 
  G4double fGdLayerThickness;
  G4double fDS20kLArThicknessAboveTPC;

  G4double fG3WTankR;
  G4double fG3IDiameter;
  G4double fAcrylicVesselDiameter;
  G4double fG3WTankHeight;

  G4bool fIsLicNdAngleByAng;
  G4bool fIsLicNdAngleByEnergy;
  G4double fLicNDNuclRecE;
  G4double fLicNDTheta;
  G4double fLicNDPhi1;
  G4double fLicNDPhi2;
  G4double fLicNDRadius;

  G4double fLicNeutronE;

  G4double fWriteFilterMinimumEnergy;
  G4double fWriteFilterMaximumEnergy;

  G4bool fLicWallActivated;
  G4double fLicWallThick;
  G4double fLicWallWindowLength;
  G4double fLicWallDistanceToTPC;

  G4bool fLicExpHall;

  // ARIS-ER parameters
  G4bool fIsAriserGeometry;
  G4double fAriserSourceDistanceToCryo;
  G4double fAriserBEGeDistanceToCryo;
  G4double fAriserBEGeAngle;
  G4ThreeVector fAriserSourcePosition;

  G4String fFLUKAfilename;

  G4int fReDConfiguration;

  NDSurveyMode fReDNDSurveyMode;
  G4String fReDNDSurveyFilename;
  NDPosMode fReDNDPositioningMode;

  G4double fTargetTPCDistance;
  G4double fThetaArX;
  G4double fPhiArX;
  G4double fNeutronAngleThetaX;
  G4double fNeutronAnglePhiX;

  G4bool fIsProtoProtoBottomSiPM;

  G4bool fGantryConfiguration;
  G4ThreeVector fGantryGlobalPosition;
  G4double fGantryPipeID;
  G4double fGantryPipeWall;
  G4double fGantryPipeBendingR;
  G4double fGantryDistanceTPCside;
  G4double fGantryDistanceTPCbottom;
  G4double fGantryDistanceBtwPipes;
  G4bool fGantrySurface;
  G4int fGantrySurfaceOption;
  G4bool fGantryShield;
  G4double fGantryShieldThickness;
  G4double fGantryShieldHeight;
  G4double fGantryShieldOffset;
  G4int fGantryShieldMaterial;
  G4double fDS20kWLSLArThickness;
  G4double fDS20kWLSPENThickness;

  G4double fDS20kCryoR;
  G4double fDS20kCryoVacuumTh;
  G4double fDS20kCryoBottomCap;
  G4double fDS20kCryoTopCap;
  G4double fDS20kCryoBottomCapOut;
  G4double fDS20kCryoTopCapOut;
  G4double fDS20kCryoTopOffset;
  G4double fDS20kCryoBottomOffset;
  G4double fDS20kCryoTopOffsetOut;
  G4double fDS20kCryoBottomOffsetOut;
  G4double fDS20kCryoH;
  G4double fDS20KCryoBarrelH;

  G4double fDS20kTPCHeight_top;
  G4double fDS20kTPCHeight_bottom;
  G4double fDS20kTPCInscribedR;

  vector<G4ThreeVector> fDS20kSiPMPosVector;

  G4double fDS20kGdConcentration;

  G4bool fOverwriteCounter ;

  G4bool fIsDecayInVacuum;

  G4bool fSkipEventsWithNoDaughters;
  G4bool fAbortRun;
  G4bool fMatScintillation;

  G4bool fPurePMMATPC;
  G4bool fHybridTPC;

  G4double fMultiplicator ;

};

#endif
