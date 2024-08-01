// --------------------------------------------------------------------------//
/**
 * AUTHOR: Davide Franco
 *
 */
// --------------------------------------------------------------------------//

#ifndef DSParametersH_H
#define DSParametersH_H
#include <iostream>
#include <vector>
#include "G4EmSaturation.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4ThreeVector.hh"
#include "math.h"

//---------------------------------------------------------------------------//

using namespace std;
class DSParameters {
 private:
  DSParameters();

 public:
  static DSParameters* Get();

  virtual ~DSParameters() {}

  void ReadDetectorGeometry();
  void ReadDetectorProperties();
  void ReadVetoPMTGeometry();
  void ReadOpticsLSV();
  void ReadOpticsTuning();
  void AdjustSiPMPDE() ; 

  G4double GetSiPMPDE(G4double);
  G4double GetTPCQE(G4double) const;
  G4double GetNVQE(G4double);
  G4double GetWorldSizeX() const { return fWorldSizeX; }
  G4double GetWorldSizeY() const { return fWorldSizeY; }
  G4double GetWorldSizeZ() const { return fWorldSizeZ; }

  // Tank
  G4double GetTankRmax() const { return fTankRmax; }
  G4double GetTankHeight() const { return fTankHeight; }

  G4double GetSSSRadius() const { return fSSSRadius; }

  // Cryostat
  G4double GetCryostatShiftZ() const { return fCryostatShiftZ; }

  G4double GetTrunkShiftZ() const { return fTrunkShiftZ; }
  G4double GetTrunkDiameter() const { return fTrunkDiameter; }
  G4double GetTrunkThickness() const { return fTrunkThickness; }
  G4double GetTrunkTopBottomOffsetZ() const { return fTrunkTopBottomOffsetZ; }
  G4double GetTrunkTopBottomOffsetX() const { return fTrunkTopBottomOffsetX; }
  G4double GetTrunkBottomHeight() const { return fTrunkBottomHeight; }
  G4double GetTrunkMiddleHeight() const { return fTrunkMiddleHeight; }
  G4double GetTrunkDistToCryostatAxis() const { return fTrunkDistToCryostatAxis; }

  // TPC
  G4double GetTPCShiftZ() const { return fTPCShiftZ; }

  G4double GetTeflonSupportDiameter() const { return fTeflonSupportDiameter; }
  G4double GetTeflonSupportThickness() const { return fTeflonSupportThickness; }
  G4double GetTeflonSupportHeight() const { return fTeflonSupportHeight; }

  G4double GetPMTAssemblyHeight() const { return fPMTAssemblyHeight; }
  G4double GetTeflonCapDiameter() const { return fTeflonCapDiameter; }
  G4double GetTeflonCapHeight() const { return fTeflonCapHeight; }

  G4double GetLArBottomLayerThickness() const { return fLArBottomLayerThickness; }
  G4double GetCathodeWindowDiameter() const { return fCathodeWindowDiameter; }
  G4double GetCathodeWindowHeight() const { return fCathodeWindowHeight; }
  G4double GetITOThickness() const { return fITOThickness; }

  G4double GetReflectorInnerDiameter() const { return fReflectorInnerDiameter; }
  G4double GetReflectorOuterDiameter() const { return fReflectorOuterDiameter; }
  G4double GetReflectorHeight() const { return fReflectorHeight; }
  G4double GetAboveGridInnerDiameter() const { return fAboveGridInnerDiameter; }
  G4double GetAboveGridOuterDiameter() const { return fAboveGridOuterDiameter; }
  G4double GetAboveGridHeight() const { return fAboveGridHeight; }
  G4double GetGasPocketThickness() const { return fGasPocketThickness; }
  G4double GetTPBThickness() const { return fTPBThickness; }
  G4double GetPENThickness() const { return fPENThickness; }

  G4double GetDivingBellHeight() const { return fDivingBellHeight; }
  G4double GetDivingBellTopHeight() const { return fDivingBellTopHeight; }
  G4double GetDivingBellOuterDiameter() const { return fDivingBellOuterDiameter; }
  G4double GetDivingBellInnerDiameter() const { return fDivingBellInnerDiameter; }

  G4double GetFieldRingsHeight() const { return fFieldRingsHeight; }
  G4double GetFieldRingsThickness() const { return fFieldRingsThickness; }

  // TPC PMTs
  G4double GetPMTBodyDiameter() const { return fPMTBodyDiameter; }
  G4double GetPMTBodyHeight() const { return fPMTBodyHeight; }
  G4double GetPMTHeadDiameter() const { return fPMTHeadDiameter; }
  G4double GetPMTHeadHeight() const { return fPMTHeadHeight; }
  G4double GetPMTWallThickness() const { return fPMTWallThickness; }
  G4double GetPMTWindowThickness() const { return fPMTWindowThickness; }
  G4double GetPMTSpacing() const { return fPMTSpacing; }
  G4double GetPMTOffset() const { return fPMTOffset; }

  // Veto PMT Positions
  G4int GetNVPMTs() const { return fVPMTNum.size(); }
  const vector<G4double>& GetVPMTNumber() const { return fVPMTNum; }
  const vector<G4double>& GetVPMTTheta() const { return fVPMTTheta; }
  const vector<G4double>& GetVPMTPhi() const { return fVPMTPhi; }

  // Properties
  G4double GetPmtMaxQe() const { return fPmtMaxQe; }
  G4double GetSiPMMaxPDE() const { return fSiPMMaxPDE; }
  G4double GetRealPmtMaxQe() const { return fRealPmtMaxQe; }
  G4double GetPmtMaxQe_adjusted() const { return fPmtMaxQe_adjusted; }
  const vector<G4double>& GetPmtQeWl() const { return fPmtQeWl; }
  const vector<G4double>& GetPmtQe() const { return fPmtQe; }
  G4double GetVPmtMaxQe() const { return fVPmtMaxQe; }
  G4double GetVPmtMaxQe_adjusted() const { return fVPmtMaxQe_adjusted; }
  const vector<G4double>& GetVPmtQeWl() const { return fVPmtQeWl; }
  const vector<G4double>& GetVPmtQe() const { return fVPmtQe; }

  inline G4double nmToEnergy(G4double wl);
  std::vector<G4double> ConvertTo4DoubleVector(const char* st);
  int binSearch(const G4double* list, const G4double& x, int minIndex, int maxIndex);
  int search(const G4double* list, const G4double& x, int listSize);
  G4double interpolate(const G4double* ys, const G4double* xs, const G4double& x, int size);

  // TPC Optics Tuning
  G4double GetGridSteelRindScale() { return fGridSteelRindScale; }
  G4double GetPhotocathodeUVRind() { return fPhotocathodeUVRind; }
  G4double GetPhotocathodeVisRind() { return fPhotocathodeVisRind; }
  G4double GetFusedSilicaUVRind() { return fFusedSilicaUVRind; }
  G4double GetFusedSilicaVisRind() { return fFusedSilicaVisRind; }
  G4double GetTPBUVRind() { return fTPBUVRind; }
  G4double GetTPBVisRind() { return fTPBVisRind; }
  G4double GetPENUVRind() { return fPENUVRind; }
  G4double GetPENVisRind() { return fPENVisRind; }
  G4double GetGaseousArgonUVAbs() { return fGaseousArgonUVAbs; }
  G4double GetGaseousArgonVisAbs() { return fGaseousArgonVisAbs; }
  G4double GetLiquidArgonUVAbs() { return fLiquidArgonUVAbs; }
  G4double GetLiquidArgonVisAbs() { return fLiquidArgonVisAbs; }
  G4double GetGridSteelUVAbs() { return fGridSteelUVAbs; }
  G4double GetGridSteelVisAbs() { return fGridSteelVisAbs; }
  G4double GetTPBUVAbs() { return fTPBUVAbs; }
  G4double GetTPBVisAbs() { return fTPBVisAbs; }
  G4double GetPENUVAbs() { return fPENUVAbs; }
  G4double GetPENVisAbs() { return fPENVisAbs; }
  G4double GetFusedSilicaUVAbs() { return fFusedSilicaUVAbs; }
  G4double GetFusedSilicaVisAbs() { return fFusedSilicaVisAbs; }
  G4double GetWLSAbsorptionFactor() { return fWLSAbsorptionFactor; }
  G4double GetWLSMeanNumberPhotons() { return fWLSMeanNumberPhotons; }
  G4double GetWLSTimeConstant_ns() { return fWLSTimeConstant_ns; }
  G4double GetWLSEfficiency() { return fWLSEfficiency; }
  G4int GetWithITO() { return fWithITO; }
  G4int GetWithGasPocket() { return fWithGasPocket; }
  G4int GetWithNewGridModel() { return fWithNewGridModel; }
  G4double GetGridNormalTransparency() { return fGridNormalTransparency; }
  G4double GetGArRindexScale() { return fGArRindexScale; }
  G4double GetLArRayleighScale() { return fLArRayleighScale; }
  G4double GetFSilicaRaylVisLength() { return fFSilicaRaylVisLength; }
  G4double GetFSilicaRaylUVLength() { return fFSilicaRaylUVLength; }
  G4double GetLArGridUVRef() { return fLArGridUVRef; }
  G4double GetLArGridVisRef() { return fLArGridVisRef; }
  G4double GetTeflonTPBUVRef() { return fTeflonTPBUVRef; }
  G4double GetTeflonTPBVisRef() { return fTeflonTPBVisRef; }
  G4double GetPENTimeConstant_ns() { return fPENTimeConstant_ns; }
  G4double GetPENEfficiency() { return fPENEfficiency; }
  G4double GetPMTLArUVRef() { return fPMTLArUVRef; }
  G4double GetPMTLArVisRef() { return fPMTLArVisRef; }
  G4double GetTeflonLArUVRef() { return fTeflonLArUVRef; }
  G4double GetTeflonLArVisRef() { return fTeflonLArVisRef; }
  G4double GetArTPBVisTran() { return fArTPBVisTran; }
  G4double GetPArRind() { return fPArRind; }
  G4double GetWithDefect() { return fWithDefect; }
  G4double GetWithSaggingModel() { return fWithSaggingModel; }
  G4double GetPENUVRaylLength() { return fPENUVRaylLength; }
  G4double GetPENVisRaylLength() { return fPENVisRaylLength; }
  void SetTeflonTPBVisRef(G4double val) { fTeflonTPBVisRef = val; }

  // LS Optics
  G4double GetLSYield() const { return fLSYield; }
  G4double GetLSBirksBeta() const { return fLSBirksBeta; }
  G4double GetLSBirksProton() const { return fLSBirksProton; }
  G4double GetLSBirksAlpha() const { return fLSBirksAlpha; }
  G4double GetLSBirks2Beta() const { return fLSBirks2Beta; }
  G4double GetLSBirks2Proton() const { return fLSBirks2Proton; }
  G4double GetLSBirks2Alpha() const { return fLSBirks2Alpha; }
  G4double GetLSResolutionScale() const { return fLSResolutionScale; }
  const std::vector<G4double>& GetLSBetaScintTime() const { return fLSBetaScintTime; }
  const std::vector<G4double>& GetLSBetaScintWeight() const { return fLSBetaScintWeight; }
  const std::vector<G4double>& GetLSAlphaScintTime() const { return fLSAlphaScintTime; }
  const std::vector<G4double>& GetLSAlphaScintWeight() const { return fLSAlphaScintWeight; }
  G4double GetLSFastTime() const { return fLSFastTime; }
  G4double GetLSSlowTime() const { return fLSSlowTime; }
  G4double GetLSYieldRatio() const { return fLSYieldRatio; }
  G4double GetLSWLSTime() const { return fLSWLSTime; }
  G4double GetLSWLSEfficiency() const { return fLSWLSEfficiency; }

  // LSV Optical Properties: constants
  G4double GetLumirrorSpecLobe() const { return fLumirrorSpecLobe; }
  G4double GetLumirrorSpecSpike() const { return fLumirrorSpecSpike; }
  G4double GetLumirrorBackscat() const { return fLumirrorBackscat; }
  // G4double GetLumirrorRefl()		const { return fLumirrorRefl; }
  G4double GetLumirrorEff() const { return fLumirrorEff; }
  G4double GetLumirrorRindex() const { return fLumirrorRindex; }
  G4double GetEpssSpecLobe() const { return fEpssSpecLobe; }
  G4double GetEpssSpecSpike() const { return fEpssSpecSpike; }
  G4double GetEpssBackscat() const { return fEpssBackscat; }
  G4double GetEpssRefl() const { return fEpssRefl; }
  G4double GetEpssEff() const { return fEpssEff; }
  G4double GetEpssRindex() const { return fEpssRindex; }
  G4double GetVCathodeEff() const { return fVCathodeEff; }
  G4double GetVCathodeSpecLobe() const { return fVCathodeSpecLobe; }
  G4double GetVCathodeSpecSpike() const { return fVCathodeSpecSpike; }
  G4double GetVCathodeBackscat() const { return fVCathodeBackscat; }
  G4double GetVCathodeRindex() const { return fVCathodeRindex; }
  G4double GetPMTBackSpecLobe() const { return fPMTBackSpecLobe; }
  G4double GetPMTBackSpecSpike() const { return fPMTBackSpecSpike; }
  G4double GetPMTBackBackscat() const { return fPMTBackBackscat; }
  G4double GetPMTBackRefl() const { return fPMTBackRefl; }
  G4double GetPMTBackEff() const { return fPMTBackEff; }
  G4double GetPMTBackRindex() const { return fPMTBackRindex; }
  G4double GetUnssSpecLobe() const { return fUnssSpecLobe; }
  G4double GetUnssSpecSpike() const { return fUnssSpecSpike; }
  G4double GetUnssBackscat() const { return fUnssBackscat; }
  G4double GetUnssRefl() const { return fUnssRefl; }
  G4double GetUnssEff() const { return fUnssEff; }
  G4double GetUnssRindex() const { return fUnssRindex; }

  // LS Optical Properties: vectors
  const vector<G4double>& GetPCAttEnergy() const { return fPCAttEnergy; }
  const vector<G4double>& GetPCAttLength() const { return fPCAttLength; }
  vector<G4double> GetPPOAttEnergy() const { return fPPOAttEnergy; }
  vector<G4double> GetPPOAttLength() const { return fPPOAttLength; }
  const vector<G4double>& GetTMBAttEnergy() const { return fTMBAttEnergy; }
  const vector<G4double>& GetTMBAttLength() const { return fTMBAttLength; }
  vector<G4double> GetPCRefrEnergy() const { return fPCRefrEnergy; }
  vector<G4double> GetPCRefrIndex() const { return fPCRefrIndex; }

  const vector<G4double>& GetLumirrorEnergy() const { return fLumirrorEnergy; }
  const vector<G4double>& GetLumirrorReflec() const { return fLumirrorReflec; }

  // Veto PMT Photo Cathode
  const vector<G4double>& GetVCathodeEnergy() const { return fVCathodeEnergy; }
  const vector<G4double>& GetVCathodeReflec() const { return fVCathodeReflec; }

  // TPC PhotoCathode
  G4int GetQENumEntries() const { return fPhotocathode_Energy.size(); }
  vector<G4double> GetPhotocathodeEnergy() { return fPhotocathode_Energy; }
  vector<G4double> GetPhotocathodeSpecularLobe() { return fPhotocathode_SpecularLobe; }
  vector<G4double> GetPhotocathodeSpecularSpike() { return fPhotocathode_SpecularSpike; }
  vector<G4double> GetPhotocathodeBackscatter() { return fPhotocathode_Backscatter; }
  vector<G4double> GetPhotocathodeReflectance() { return fPhotocathode_Reflectance; }
  vector<G4double> GetPhotocathodeEfficiency() { return fPhotocathode_Efficiency; }
  vector<G4double> GetPhotocathodeRindex() { return fPhotocathode_Rindex; }

  // Cryostat Sheath Properties
  G4double* GetCryoSheathR() { return fCryoSheathR; }
  G4double* GetCryoSheathZ() { return fCryoSheathZ; }

  // EM Saturation
  void AddSaturation(G4EmSaturation* sat) { emSaturation = sat; }
  void RemoveSaturation() { emSaturation = NULL; }
  G4EmSaturation* GetSaturation() { return emSaturation; }

 private:
  static DSParameters* me;

  G4String fDSGeometryFileName;

  G4float fWorldSizeX;
  G4float fWorldSizeY;
  G4float fWorldSizeZ;
  G4float fTankRmax;
  G4float fTankHeight;

  G4double fSSSRadius;

  // Cryostat
  G4double fCryostatShiftZ;

  G4double fTrunkShiftZ;
  G4double fTrunkDiameter;
  G4double fTrunkThickness;
  G4double fTrunkTopBottomOffsetZ;
  G4double fTrunkTopBottomOffsetX;
  G4double fTrunkBottomHeight;
  G4double fTrunkMiddleHeight;
  G4double fTrunkDistToCryostatAxis;

  // TPC
  G4double fTPCShiftZ;

  G4double fTeflonSupportDiameter;
  G4double fTeflonSupportThickness;
  G4double fTeflonSupportHeight;

  G4double fPMTAssemblyHeight;
  G4double fTeflonCapDiameter;
  G4double fTeflonCapHeight;

  G4double fLArBottomLayerThickness;
  G4double fCathodeWindowDiameter;
  G4double fCathodeWindowHeight;
  G4double fITOThickness;

  G4double fReflectorInnerDiameter;
  G4double fReflectorOuterDiameter;
  G4double fReflectorHeight;
  G4double fAboveGridInnerDiameter;
  G4double fAboveGridOuterDiameter;
  G4double fAboveGridHeight;
  G4double fGasPocketThickness;
  G4double fTPBThickness;
  G4double fPENThickness;

  G4double fDivingBellHeight;
  G4double fDivingBellTopHeight;
  G4double fDivingBellOuterDiameter;
  G4double fDivingBellInnerDiameter;

  G4double fFieldRingsHeight;
  G4double fFieldRingsThickness;

  // TPC PMTs
  G4double fPMTBodyDiameter;
  G4double fPMTBodyHeight;
  G4double fPMTHeadDiameter;
  G4double fPMTHeadHeight;
  G4double fPMTWallThickness;
  G4double fPMTWindowThickness;
  G4double fPMTSpacing;
  G4double fPMTOffset;

  // Veto PMT Positions
  vector<G4double> fVPMTNum;
  vector<G4double> fVPMTTheta;
  vector<G4double> fVPMTPhi;

  // TPC Optics Tuning
  G4double fGridSteelRindScale;
  G4double fPhotocathodeUVRind;
  G4double fPhotocathodeVisRind;
  G4double fFusedSilicaUVRind;
  G4double fFusedSilicaVisRind;
  G4double fTPBUVRind;
  G4double fTPBVisRind;
  G4double fPENUVRind;
  G4double fPENVisRind;
  G4double fGaseousArgonUVAbs;
  G4double fGaseousArgonVisAbs;
  G4double fLiquidArgonUVAbs;
  G4double fLiquidArgonVisAbs;
  G4double fGridSteelUVAbs;
  G4double fGridSteelVisAbs;
  G4double fFusedSilicaUVAbs;
  G4double fFusedSilicaVisAbs;
  G4double fTPBUVAbs;
  G4double fTPBVisAbs;
  G4double fPENUVAbs;
  G4double fPENVisAbs;
  G4double fWLSAbsorptionFactor;
  G4double fWLSMeanNumberPhotons;
  G4double fWLSTimeConstant_ns;
  G4double fWLSEfficiency;
  G4double fPENTimeConstant_ns;
  G4double fPENEfficiency;
  G4int fWithITO;
  G4int fWithGasPocket;
  G4int fWithNewGridModel;
  G4double fGridNormalTransparency;
  G4double fGArRindexScale;
  G4double fLArRayleighScale;
  G4double fFSilicaRaylVisLength;
  G4double fFSilicaRaylUVLength;
  G4double fLArGridUVRef;
  G4double fLArGridVisRef;
  G4double fTeflonTPBUVRef;
  G4double fTeflonTPBVisRef;
  G4double fPENUVRaylLength;
  G4double fPENVisRaylLength;
  G4double fPMTLArUVRef;
  G4double fPMTLArVisRef;
  G4double fTeflonLArUVRef;
  G4double fTeflonLArVisRef;
  G4double fArTPBVisTran;
  G4double fPArRind;
  G4int fWithDefect;
  G4int fWithSaggingModel;

  // LSV Optics
  G4double fLSYield;
  G4double fLSBirksBeta;
  G4double fLSBirksProton;
  G4double fLSBirksAlpha;
  G4double fLSBirks2Beta;
  G4double fLSBirks2Proton;
  G4double fLSBirks2Alpha;
  G4double fLSResolutionScale;
  std::vector<G4double> fLSBetaScintTime;
  std::vector<G4double> fLSBetaScintWeight;
  std::vector<G4double> fLSAlphaScintTime;
  std::vector<G4double> fLSAlphaScintWeight;
  G4double fLSFastTime;
  G4double fLSSlowTime;
  G4double fLSYieldRatio;
  G4double fLSWLSTime;
  G4double fLSWLSEfficiency;
  // Lumirror boundary properties
  G4double fLumirrorEff;
  G4double fLumirrorSpecLobe;
  G4double fLumirrorSpecSpike;
  G4double fLumirrorBackscat;
  G4double fLumirrorRindex;
  // EPSS (ElectroPolished Stainless Steel) boundary properties: boundary
  // between LS and TPC cryostat
  G4double fEpssRefl;
  G4double fEpssEff;
  G4double fEpssSpecLobe;
  G4double fEpssSpecSpike;
  G4double fEpssBackscat;
  G4double fEpssRindex;
  // LSV PMT PhotoCathode
  G4double fVCathodeRefl;
  G4double fVCathodeEff;
  G4double fVCathodeSpecLobe;
  G4double fVCathodeSpecSpike;
  G4double fVCathodeBackscat;
  G4double fVCathodeRindex;
  // LSV PMT Back
  G4double fPMTBackRefl;
  G4double fPMTBackEff;
  G4double fPMTBackSpecLobe;
  G4double fPMTBackSpecSpike;
  G4double fPMTBackBackscat;
  G4double fPMTBackRindex;
  // USS (Untreated StainlesSteel) boundary properties: boundary between LS and
  // flanges
  G4double fUnssRefl;
  G4double fUnssEff;
  G4double fUnssSpecLobe;
  G4double fUnssSpecSpike;
  G4double fUnssBackscat;
  G4double fUnssRindex;

  vector<G4double> fPCAttEnergy;
  vector<G4double> fPCAttLength;
  vector<G4double> fPPOAttEnergy;
  vector<G4double> fPPOAttLength;
  vector<G4double> fTMBAttEnergy;
  vector<G4double> fTMBAttLength;
  vector<G4double> fPCRefrEnergy;
  vector<G4double> fPCRefrIndex;

  // LSV Optical Boundary Properties
  vector<G4double> fLumirrorEnergy;
  vector<G4double> fLumirrorReflec;
  vector<G4double> fVCathodeEnergy;
  vector<G4double> fVCathodeReflec;

  // TPC PMT Cathode surface (are they used??)
  vector<G4double> fPhotocathode_Energy;
  vector<G4double> fPhotocathode_Reflectance;
  vector<G4double> fPhotocathode_Efficiency;
  vector<G4double> fPhotocathode_SpecularLobe;
  vector<G4double> fPhotocathode_SpecularSpike;
  vector<G4double> fPhotocathode_Backscatter;
  vector<G4double> fPhotocathode_Rindex;

  // Properties
  G4double fPmtMaxQe;
  G4double fSiPMMaxPDE;
  G4double fRealPmtMaxQe;
  G4double fPmtMaxQe_adjusted;
  G4double fVPmtMaxQe;
  G4double fVPmtMaxQe_adjusted;

  std::vector<G4double> fPmtQeWl;
  std::vector<G4double> fPmtQe;
  std::vector<G4double> fVPmtQeWl;
  std::vector<G4double> fVPmtQe;

  //DS20k SiPM PDE from reference 
  bool IsSipmPDEAdjusted = false; 
  G4double myPDEWL[16] = {285.87, 301.14, 320.99, 340.83, 365.26, 380.53, 391.22, 400.38, 411.06, 420.22, 435.49, 449.23, 464.50, 499.61, 588.16, 698.09};  // in nm
  G4double myPDE[16] = {18.39, 37.24, 41.87, 42.43, 45.03, 48.12, 52.72, 52.75, 52.77, 52.29, 48.77, 46.26, 40.71, 38.76, 28.82, 18.93};                    // PDE %

  // Dimensions for the lumirror sheath around the outside of the cryostat
  G4double fCryoSheathR[2];
  G4double fCryoSheathZ[2];

  G4EmSaturation* emSaturation;
};

#endif
/*
 * $Log: DSParameters.hh,v $
 * Revision 1.10  2016/03/11 13:19:57  davini
 * returning const vector<G4double>& instead of vector<G4double> where possible
 *
 * Revision 1.9  2016/03/11 11:41:32  davini
 * GetVPmtQe GetVPMTQeWl returns const vector<G4double>& instead of
 * vector<G4double>. Enforced use of const double* (instead of double*) for
 * passing arrays that should not be modyfied by functions
 *
 * Revision 1.8  2016/03/03 23:06:15  cz2
 * Added parameters related to TPB defects and S2 radial correction.
 *
 * Revision 1.7  2015/01/17 11:31:45  pagnes
 * PAr model added form optical tuning
 *
 * Revision 1.6  2015/01/07 16:45:59  pagnes
 * changed veto optical properties format from arrays to vectors
 *
 * Revision 1.5  2014/11/05 15:47:09  pagnes
 * temporary optics tuning
 *
 * Revision 1.4  2014/10/13 18:43:46  swesterd
 * fixed veto PMT positions
 *
 * Revision 1.3  2014/07/16 08:23:12  pagnes
 * QE scaling to 1.0 added (/ds/manager/fast_simulation xxx)
 *
 * Revision 1.2  2014/06/03 13:31:33  meyers
 * Migrate TPC grid and ITO optics updates to g4ds10
 *
 * Revision 1.1  2014/05/07 12:20:54  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.17  2014/03/19 16:37:36  dfranco
 * update external configuration reader for optics tuning
 *
 * Revision 1.16  2014/03/11 16:50:00  dfranco
 * read optical parameters from external file (DSOptics.dat)
 *
 * Revision 1.15  2013/08/27 07:13:12  swesterd
 * add visible energy for the neutron veto
 *
 * Revision 1.14  2013/06/19 18:35:25  swesterd
 * added DSScintCelll and made tpc PMTs' QE and reflections work like veto PMTs
 *
 * Revision 1.13  2013/06/05 23:03:28  swesterd
 * moved optical boundary MPTs to DSMaterial and gave the trunks optical
 * boundary properties consistent with untreated stainless steel
 *
 * Revision 1.12  2013/06/04 23:38:46  swesterd
 * Changed the length of the Lumirror reflectance array to be an adjustable
 * size, set it to 60 elements
 *
 * Revision 1.11  2013/06/04 14:11:38  dfranco
 * Added a function returning the TPC QE in DSParameters, applied in the
 * tracking action
 *
 * Revision 1.10  2013/05/27 23:59:00  swesterd
 * added a (currently commented out) Lumirror sheath to the cryostat and
 * introduced DSOpBoundaryProcess to try to figure out why the boundaries are
 * being screwy, with some edits so that it can handle constant and vector
 * properties with freaking out
 *
 * Revision 1.9  2013/05/25 07:58:22  swesterd
 * Got the veto PMT optical boundaries all working along with photocathode
 * optical properties, added PMT quantum efficiency to DSTrackingAction, and
 * added a function to DSTrackingAction that locates and quadratically
 * interpolates points in data, for getting useful QEs
 *
 * Revision 1.8  2013/05/07 23:06:26  swesterd
 * added optical boundaries and Lumirror in the veto
 *
 * Revision 1.7  2013/05/06 14:59:52  perassos
 * Updates on the TPC surface properties and geometry
 *
 * Revision 1.6  2013/05/01 08:20:23  swesterd
 * added boron-loaded scintillator optical properties
 *
 * Revision 1.5  2013/04/26 16:17:19  perassos
 * TPC PMT QE added
 *
 * Revision 1.4  2013/04/19 16:10:09  perassos
 * Added ITO and TPB layers
 *
 * Revision 1.3  2013/03/22 14:09:40  dfranco
 * Added the cvs logger code at the end of each file
 *
 */

