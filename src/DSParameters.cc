// --------------------------------------------------------------------------//
/**
 * AUTHOR: Davide Franco
 *
 */
// --------------------------------------------------------------------------//

//---------------------------------------------------------------------------//

#include "DSParameters.hh"  //Present DS Class Headers
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//---------------------------------------------------------------------------//

#include <algorithm>
#include <fstream>
#include <sstream>
#include "DSIO.hh"
#include "DSLogger.hh"
#include "DSStorage.hh"

using namespace std;

DSParameters* DSParameters::me = 0;

// singleton
DSParameters::DSParameters() {
  ReadDetectorGeometry();
  ReadDetectorProperties();
  ReadVetoPMTGeometry();
  ReadOpticsLSV();
  ReadOpticsTuning();
}

DSParameters* DSParameters::Get() {
  if (!me) me = new DSParameters();
  return me;
}

void DSParameters::ReadDetectorGeometry() {
  G4String mys;

  //DSIO::Get()->GetStreamLogFile() << endl;
  //DSIO::Get()->GetStreamLogFile() << "########## Detector Geometry ###########" << endl;

  while (getline(DSIO::Get()->GetStreamDSGeometry(), mys)) {

    //DSIO::Get()->GetStreamLogFile() << mys << endl;

    if (!mys.find("WorldSizeX")) {
      fWorldSizeX = atof(mys.substr(mys.find("=") + 1).c_str()) * mm;
    } else if (!mys.find("WorldSizeY")) {
      fWorldSizeY = atof(mys.substr(mys.find("=") + 1).c_str()) * mm;
    } else if (!mys.find("WorldSizeZ")) {
      fWorldSizeZ = atof(mys.substr(mys.find("=") + 1).c_str()) * mm;
    } else if (!mys.find("DSTankHeight")) {
      fTankHeight = atof(mys.substr(mys.find("=") + 1).c_str()) * mm;
    } else if (!mys.find("DSTankRMax")) {
      fTankRmax = atof(mys.substr(mys.find("=") + 1).c_str()) * mm;
    } else if (!mys.find("SSSRadius")) {
      fSSSRadius = atof(mys.substr(mys.find("=") + 1).c_str()) * mm;
    } else if (!mys.find("CryostatShiftZ")) {
      fCryostatShiftZ = atof(mys.substr(mys.find("=") + 1).c_str()) * mm;
    } else if (!mys.find("TrunkShiftZ")) {
      fTrunkShiftZ = atof(mys.substr(mys.find("=") + 1).c_str()) * mm;
    } else if (!mys.find("TrunkDiameter")) {
      fTrunkDiameter = atof(mys.substr(mys.find("=") + 1).c_str()) * mm;
    } else if (!mys.find("TrunkThickness")) {
      fTrunkThickness = atof(mys.substr(mys.find("=") + 1).c_str()) * mm;
    } else if (!mys.find("TrunkTopBottomOffsetZ")) {
      fTrunkTopBottomOffsetZ = atof(mys.substr(mys.find("=") + 1).c_str()) * mm;
    } else if (!mys.find("TrunkTopBottomOffsetX")) {
      fTrunkTopBottomOffsetX = atof(mys.substr(mys.find("=") + 1).c_str()) * mm;
    } else if (!mys.find("TrunkBottomHeight")) {
      fTrunkBottomHeight = atof(mys.substr(mys.find("=") + 1).c_str()) * mm;
    } else if (!mys.find("TrunkMiddleHeight")) {
      fTrunkMiddleHeight = atof(mys.substr(mys.find("=") + 1).c_str()) * mm;
    } else if (!mys.find("TrunkDistToCryostatAxis")) {
      fTrunkDistToCryostatAxis = atof(mys.substr(mys.find("=") + 1).c_str()) * mm;
    } else if (!mys.find("TPCShiftZ")) {
      fTPCShiftZ = atof(mys.substr(mys.find("=") + 1).c_str()) * mm;
    } else if (!mys.find("TeflonSupportDiameter")) {
      fTeflonSupportDiameter = atof(mys.substr(mys.find("=") + 1).c_str()) * mm;
    } else if (!mys.find("TeflonSupportThickness")) {
      fTeflonSupportThickness = atof(mys.substr(mys.find("=") + 1).c_str()) * mm;
    } else if (!mys.find("TeflonSupportHeight")) {
      fTeflonSupportHeight = atof(mys.substr(mys.find("=") + 1).c_str()) * mm;
    } else if (!mys.find("TeflonSupportHeight")) {
      fTeflonSupportHeight = atof(mys.substr(mys.find("=") + 1).c_str()) * mm;
    } else if (!mys.find("PMTAssemblyHeight")) {
      fPMTAssemblyHeight = atof(mys.substr(mys.find("=") + 1).c_str()) * mm;
    } else if (!mys.find("TeflonCapDiameter")) {
      fTeflonCapDiameter = atof(mys.substr(mys.find("=") + 1).c_str()) * mm;
    } else if (!mys.find("TeflonCapHeight")) {
      fTeflonCapHeight = atof(mys.substr(mys.find("=") + 1).c_str()) * mm;
    } else if (!mys.find("LArBottomLayerThickness")) {
      fLArBottomLayerThickness = atof(mys.substr(mys.find("=") + 1).c_str()) * mm;
    } else if (!mys.find("CathodeWindowDiameter")) {
      fCathodeWindowDiameter = atof(mys.substr(mys.find("=") + 1).c_str()) * mm;
    } else if (!mys.find("CathodeWindowHeight")) {
      fCathodeWindowHeight = atof(mys.substr(mys.find("=") + 1).c_str()) * mm;
    } else if (!mys.find("ITOThickness")) {
      fITOThickness = atof(mys.substr(mys.find("=") + 1).c_str()) * mm;
    } else if (!mys.find("ReflectorInnerDiameter")) {
      fReflectorInnerDiameter = atof(mys.substr(mys.find("=") + 1).c_str()) * mm;
    } else if (!mys.find("ReflectorOuterDiameter")) {
      fReflectorOuterDiameter = atof(mys.substr(mys.find("=") + 1).c_str()) * mm;
    } else if (!mys.find("ReflectorHeight")) {
      fReflectorHeight = atof(mys.substr(mys.find("=") + 1).c_str()) * mm;
    } else if (!mys.find("AboveGridInnerDiameter")) {
      fAboveGridInnerDiameter = atof(mys.substr(mys.find("=") + 1).c_str()) * mm;
    } else if (!mys.find("AboveGridOuterDiameter")) {
      fAboveGridOuterDiameter = atof(mys.substr(mys.find("=") + 1).c_str()) * mm;
    } else if (!mys.find("AboveGridHeight")) {
      fAboveGridHeight = atof(mys.substr(mys.find("=") + 1).c_str()) * mm;
    } else if (!mys.find("GasPocketThickness")) {
      fGasPocketThickness = atof(mys.substr(mys.find("=") + 1).c_str()) * mm;
    } else if (!mys.find("TPBThickness")) {
      fTPBThickness = atof(mys.substr(mys.find("=") + 1).c_str()) * mm;
    } else if (!mys.find("PENThickness")) {
      fPENThickness = atof(mys.substr(mys.find("=") + 1).c_str()) * mm;
    } else if (!mys.find("DivingBellHeight")) {
      fDivingBellHeight = atof(mys.substr(mys.find("=") + 1).c_str()) * mm;
    } else if (!mys.find("DivingBellTopHeight")) {
      fDivingBellTopHeight = atof(mys.substr(mys.find("=") + 1).c_str()) * mm;
    } else if (!mys.find("DivingBellOuterDiameter")) {
      fDivingBellOuterDiameter = atof(mys.substr(mys.find("=") + 1).c_str()) * mm;
    } else if (!mys.find("DivingBellInnerDiameter")) {
      fDivingBellInnerDiameter = atof(mys.substr(mys.find("=") + 1).c_str()) * mm;
    } else if (!mys.find("FieldRingsHeight")) {
      fFieldRingsHeight = atof(mys.substr(mys.find("=") + 1).c_str()) * mm;
    } else if (!mys.find("FieldRingsThickness")) {
      fFieldRingsThickness = atof(mys.substr(mys.find("=") + 1).c_str()) * mm;
    } else if (!mys.find("PMTBodyDiameter")) {
      fPMTBodyDiameter = atof(mys.substr(mys.find("=") + 1).c_str()) * mm;
    } else if (!mys.find("PMTBodyHeight")) {
      fPMTBodyHeight = atof(mys.substr(mys.find("=") + 1).c_str()) * mm;
    } else if (!mys.find("PMTHeadDiameter")) {
      fPMTHeadDiameter = atof(mys.substr(mys.find("=") + 1).c_str()) * mm;
    } else if (!mys.find("PMTHeadHeight")) {
      fPMTHeadHeight = atof(mys.substr(mys.find("=") + 1).c_str()) * mm;
    } else if (!mys.find("PMTWallThickness")) {
      fPMTWallThickness = atof(mys.substr(mys.find("=") + 1).c_str()) * mm;
    } else if (!mys.find("PMTWindowThickness")) {
      fPMTWindowThickness = atof(mys.substr(mys.find("=") + 1).c_str()) * mm;
    } else if (!mys.find("PMTSpacing")) {
      fPMTSpacing = atof(mys.substr(mys.find("=") + 1).c_str()) * mm;
    } else if (!mys.find("PMTOffset")) {
      fPMTOffset = atof(mys.substr(mys.find("=") + 1).c_str()) * mm;
    }
  }

  DSIO::Get()->CloseStreamDSGeometry();
  DSLog(trace) << "Detector Geometry loaded" << endlog;
}

void DSParameters::ReadVetoPMTGeometry() {
  DSLog(routine) << "Reading LSV PMT Geometry from file: " << DSIO::Get()->GetDSVPMTGeometry() << endlog;
  G4double pmtnum, theta, phi;  // theta = polar angle, phi = azimuthal angle
  //DSIO::Get()->GetStreamLogFile() << endl;
  //DSIO::Get()->GetStreamLogFile() << "########## Veto PMT Positions ###########" << endl;
  fVPMTNum.clear();
  fVPMTTheta.clear();
  fVPMTPhi.clear();

  int index = 0;
  if (!(DSStorage::Get()->GetIsDS20kLSV())) {  // DS-50
    DSLog(routine) << "Reading LSV PMT Geometry for DS-50 " << endlog;
    while (DSIO::Get()->GetStreamDSVPMTGeometry() >> pmtnum >> theta >> phi) {
      if (index >= 110) break;
      //DSIO::Get()->GetStreamLogFile() << index << ", pmt " << pmtnum << "\ttheta = " << theta << "\tphi = " << phi << endl;
      fVPMTNum.push_back(pmtnum);
      fVPMTTheta.push_back(theta);
      fVPMTPhi.push_back(phi);
      index++;
    }
  } else {  // DS-20k
    DSLog(routine) << "Reading LSV PMT Geometry for DS-20k " << endlog;
    while (DSIO::Get()->GetStreamDSVPMTGeometry() >> theta >> phi) {
      pmtnum = index;
      //DSIO::Get()->GetStreamLogFile() << index << ", pmt " << pmtnum << "\ttheta = " << theta << "\tphi = " << phi << endl;
      fVPMTNum.push_back(pmtnum);
      fVPMTTheta.push_back(theta);
      fVPMTPhi.push_back(phi);
      index++;
    }
  }
  DSIO::Get()->CloseStreamDSVPMTGeometry();
}

void DSParameters::ReadOpticsTuning() {
  G4String mys;
  //DSIO::Get()->GetStreamLogFile() << endl;
  //DSIO::Get()->GetStreamLogFile() << "########## Detector Optics Tuning ###########" << endl;

  while (getline(DSIO::Get()->GetStreamDSOptics(), mys)) {
    //DSIO::Get()->GetStreamLogFile() << mys << endl;

    if (!mys.find("GridSteelRindScale")) {
      fGridSteelRindScale = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << "Scale factor Grid Steel RINDEX: " <<
      // fGridSteelRindScale << endlog;
    } else if (!mys.find("PhotocathodeUVRind")) {
      fPhotocathodeUVRind = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << "Photocathode RINDEX UV light: " << fPhotocathodeUVRind
      // << endlog;
    } else if (!mys.find("PhotocathodeVisRind")) {
      fPhotocathodeVisRind = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << "Photocathode RINDEX Visible light: " <<
      // fPhotocathodeVisRind << endlog;
    } else if (!mys.find("FusedSilicaUVRind")) {
      fFusedSilicaUVRind = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << "Fused Silica RINDEX UV light: " << fFusedSilicaUVRind
      // << endlog;
    } else if (!mys.find("FusedSilicaVisRind")) {
      fFusedSilicaVisRind = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << "Fused Silica RINDEX visible light: " <<
      // fFusedSilicaVisRind << endlog;
    } else if (!mys.find("TPBUVRind")) {
      fTPBUVRind = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << "TPB RINDEX UV light: " << fTPBUVRind << endlog;
    } else if (!mys.find("TPBVisRind")) {
      fTPBVisRind = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << "TPB RINDEX Visible light: " << fTPBVisRind << endlog;
    } else if (!mys.find("PENUVRind")) {
      fPENUVRind = atof(mys.substr(mys.find("=") + 1).c_str());
    } else if (!mys.find("PENVisRind")) {
      fPENVisRind = atof(mys.substr(mys.find("=") + 1).c_str());
    } else if (!mys.find("GaseousArgonUVAbs")) {
      fGaseousArgonUVAbs = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << "Gaseous Ar Absorption Length UV light: " <<
      // fGaseousArgonUVAbs << endlog;
    } else if (!mys.find("GaseousArgonVisAbs")) {
      fGaseousArgonVisAbs = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << "Gaseous Ar Absorption Length Vis light: " <<
      // fGaseousArgonVisAbs << endlog;
    } else if (!mys.find("LiquidArgonUVAbs")) {
      fLiquidArgonUVAbs = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << "Liquid Ar Absorption Length UV light: " <<
      // fLiquidArgonUVAbs << endlog;
    } else if (!mys.find("LiquidArgonVisAbs")) {
      fLiquidArgonVisAbs = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << "Liquid Ar Absorption Length Visible light: " <<
      // fLiquidArgonVisAbs << endlog;
    } else if (!mys.find("GridSteelUVAbs")) {
      fGridSteelUVAbs = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << "GridSteel Absorption Length UV light: " <<
      // fGridSteelUVAbs << endlog;
    } else if (!mys.find("GridSteelVisAbs")) {
      fGridSteelVisAbs = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << "GridSteel Absorption Length Visible light: " <<
      // fGridSteelVisAbs << endlog;
    } else if (!mys.find("FusedSilicaUVAbs")) {
      fFusedSilicaUVAbs = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << "FusedSilica Absorption Length UV light: " <<
      // fFusedSilicaUVAbs << endlog;
    } else if (!mys.find("FusedSilicaVisAbs")) {
      fFusedSilicaVisAbs = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << "FusedSilica Absorption Length Visible light: " <<
      // fFusedSilicaVisAbs << endlog;
    } else if (!mys.find("TPBUVAbs")) {
      fTPBUVAbs = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << "TPB Absorption Length UV light: " << fTPBUVAbs <<
      // endlog;
    } else if (!mys.find("TPBVisAbs")) {
      fTPBVisAbs = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << "TPB Absorption Length Visible light: " << fTPBVisAbs
      // << endlog;
    } else if (!mys.find("WLSAbsorptionFactor")) {
      fWLSAbsorptionFactor = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << "TPB W.L.Shifting-Length scaling factor: " <<
      // fWLSAbsorptionFactor << endlog;
    } else if (!mys.find("WLSMeanNumberPhotons")) {
      fWLSMeanNumberPhotons = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << "TPB mean number of photons: " << fWLSMeanNumberPhotons
      // << endlog;
    } else if (!mys.find("WLSTimeConstant_ns")) {
      fWLSTimeConstant_ns = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << "TPB time constant (ns): " << fWLSTimeConstant_ns <<
      // endlog;
    } else if (!mys.find("WLSEfficiency")) {
      fWLSEfficiency = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << "TPB Efficiency: " << fWLSEfficiency << endlog;
    } else if (!mys.find("PENTimeConstant_ns")) {
      fPENTimeConstant_ns = atof(mys.substr(mys.find("=") + 1).c_str());
    } else if (!mys.find("PENEfficiency")) {
      fPENEfficiency = atof(mys.substr(mys.find("=") + 1).c_str());
    } else if (!mys.find("PENUVAbs")) {
      fPENUVAbs = atof(mys.substr(mys.find("=") + 1).c_str());
    } else if (!mys.find("PENVisAbs")) {
      fPENVisAbs = atof(mys.substr(mys.find("=") + 1).c_str());
    } else if (!mys.find("PENUVRaylLength")) {
      fPENUVRaylLength = atof(mys.substr(mys.find("=") + 1).c_str());
    } else if (!mys.find("PENVisRaylLength")) {
      fPENVisRaylLength = atof(mys.substr(mys.find("=") + 1).c_str());
    } else if (!mys.find("WLSAbsorptionFactor")) {
      fWLSAbsorptionFactor = atof(mys.substr(mys.find("=") + 1).c_str());
    } else if (!mys.find("GArRindexScale")) {
      fGArRindexScale = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << "GAr Rayleigh scattering scaling: " <<
      // fGArRayleighScale << endlog;
    } else if (!mys.find("LArRayleighScale")) {
      fLArRayleighScale = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << "LAr Rayleigh scattering scaling: " <<
      // fLArRayleighScale << endlog;
    } else if (!mys.find("FSilicaRaylVisLength")) {
      fFSilicaRaylVisLength = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << "LAr Rayleigh scattering scaling: " <<
      // fLArRayleighScale << endlog;
    } else if (!mys.find("FSilicaRaylUVLength")) {
      fFSilicaRaylUVLength = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << "LAr Rayleigh scattering scaling: " <<
      // fLArRayleighScale << endlog;
    } else if (!mys.find("WithITO")) {
      fWithITO = atoi(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << "ITO placement (bool): " << fWithITO << endlog;
    } else if (!mys.find("WithGasPocket")) {
      fWithGasPocket = atoi(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << "GasPocket placement (bool): " << fWithGasPocket <<
      // endlog;
    } else if (!mys.find("WithNewGridModel")) {
      fWithNewGridModel = atoi(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << "GasPocket placement (bool): " << fWithNewGridModel <<
      // endlog;
    } else if (!mys.find("LArGridUVRef")) {
      fLArGridUVRef = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << "LAr - Grid Surf Reflectivity (1 - abs), UV light: " <<
      // fLArGridUVRef << endlog;
    } else if (!mys.find("LArGridVisRef")) {
      fLArGridVisRef = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << "LAr - Grid Surf Reflectivity (1 - abs), vis light: "
      // << fLArGridVisRef << endlog;
    }
    // P. Meyers: eliminated GridReflection, added GridNormalTransparency
    else if (!mys.find("GridNormalTransparency")) {
      fGridNormalTransparency = atof(mys.substr(mys.find("=") + 1).c_str());
    } else if (!mys.find("TeflonTPBUVRef")) {
      fTeflonTPBUVRef = atof(mys.substr(mys.find("=") + 1).c_str());
    } else if (!mys.find("TeflonTPBVisRef")) {
      fTeflonTPBVisRef = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << "Teflon-TPB Surf Reflectivity (1 - abs), vis light: "
      // << fTeflonTPBVisRef << endlog;
    } else if (!mys.find("TeflonLArUVRef")) {
      fTeflonLArUVRef = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << "Teflon-LAr Surf Reflectivity (1 - abs), UV light: " <<
      // fTeflonLArUVRef << endlog;
    } else if (!mys.find("TeflonLArVisRef")) {
      fTeflonLArVisRef = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << "Teflon-LAr Surf Reflectivity (1 - abs), vis light: "
      // << fTeflonLArVisRef << endlog;
    } else if (!mys.find("PMTLArUVRef")) {
      fPMTLArUVRef = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << "PMT-LAr Surf Reflectivity (1 - abs), UV light: " <<
      // fPMTLArUVRef << endlog;
    } else if (!mys.find("PMTLArVisRef")) {
      fPMTLArVisRef = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << "PMT-LAr Surf Reflectivity (1 - abs), vis light: " <<
      // fPMTLArVisRef << endlog;
    } else if (!mys.find("ArTPBVisTran")) {
      fArTPBVisTran = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << "PMT-LAr Surf Reflectivity (1 - abs), vis light: " <<
      // fPMTLArVisRef << endlog;
    } else if (!mys.find("PArRind")) {
      fPArRind = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << "Index of refraction of the pseudo argon: " << PArRind
      // << endlog;
    } else if (!mys.find("WithDefect")) {
      fWithDefect = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << "Flag determines whether to apply TPB defect: " <<
      // WithDefect << endlog;
    } else if (!mys.find("WithSaggingModel")) {
      fWithSaggingModel = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << "Flag determines whether to use sagging model: " <<
      // WithSaggingModel << endlog;
    }
  }
  DSIO::Get()->CloseStreamDSOptics();

  DSLog(trace) << "#######    Detector Optics loaded      #######" << endlog;
}

void DSParameters::ReadDetectorProperties() {

  // Load Detector Properties
  ifstream myPropertiesFile("../data/detector/DSProperties.dat");
  string myString;

  if (!myPropertiesFile.is_open()) DSLog(fatal) << "Couldn't load Detector Properties" << endlog;
  while (getline(myPropertiesFile, myString)) {

    if (!myString.find("PmtMaxQe")) {
      fPmtMaxQe = atof(myString.substr(myString.find("=") + 1).c_str());
      fRealPmtMaxQe = fPmtMaxQe;
    }
    if (!myString.find("SiPMMaxPDE")) {
      fSiPMMaxPDE = atof(myString.substr(myString.find("=") + 1).c_str());
    }
  }
  myPropertiesFile.close();
  DSLog(trace) << "Detector Properties loaded" << endlog;

  // Load TPC PMT Q.E.
  ifstream myQeFile("../data/detector/DSPmtQeDS50.dat");
  if (!myQeFile.is_open()) DSLog(fatal) << "Couldn't load the TPC PMT Q.E." << endlog;
  while (getline(myQeFile, myString)) {

    fPmtQeWl.push_back(atof(myString.substr(0, myString.find(",")).c_str()));
    fPmtQe.push_back(atof(myString.substr(myString.find(",") + 1).c_str()));
  }
  myQeFile.close();
  DSLog(trace) << "TPC PMT Q.E. loaded." << endlog;

  // TPC PMT QEs
  // Is this section deprecated? SD
  const G4double hc = h_Planck * c_light / eV / nm;  // eV*nm
  // Define constant properties
  G4double tpc_const_refl = 0.2;
  G4double tpc_const_eff = 0.;
  G4double tpc_const_speclobe = 1.;
  G4double tpc_const_specspike = 0.;
  G4double tpc_const_backscat = 0.;
  G4double tpc_const_rindex = 1.;
  G4double wavelength, qe, qe_adjusted, refl;
  fPmtMaxQe_adjusted = fPmtMaxQe / (1 - tpc_const_refl * (1 - fPmtMaxQe));
  // Load properties into arrays
  const char tqeFileName[] = "../data/detector/PMTTPC.dat";  // "../data/detector/DSPmtQeDS50.dat";
  myPropertiesFile.open(tqeFileName);
  if (!myPropertiesFile.is_open()) DSLog(fatal) << "WARNING : COULD NOT OPEN TPC PMT DATA FILE " << tqeFileName << endlog;
  while (myPropertiesFile >> wavelength >> qe) {
    qe_adjusted = fPmtMaxQe * qe / (1 - tpc_const_refl * (1 - fPmtMaxQe * qe));
    qe_adjusted /= fPmtMaxQe_adjusted;
    refl = (1 - qe) * tpc_const_refl;
    refl *= fPmtMaxQe_adjusted;

    // fPmtQe.push_back(qe_adjusted);
    // fPmtQeWl.push_back(wavelength*nm);
    fPhotocathode_Energy.push_back((hc / wavelength) * eV);
    fPhotocathode_Reflectance.push_back(refl);
    fPhotocathode_Efficiency.push_back(tpc_const_eff);
    fPhotocathode_SpecularLobe.push_back(tpc_const_speclobe);
    fPhotocathode_SpecularSpike.push_back(tpc_const_specspike);
    fPhotocathode_Backscatter.push_back(tpc_const_backscat);
    fPhotocathode_Rindex.push_back(tpc_const_rindex);

    DSLog(debugging) << "TPC PMT QE : WL = " << wavelength << " nm\tQE = " << qe << "\tQE' = " << qe_adjusted << "\tr = " << refl << endlog;
  }
  std::reverse(fPhotocathode_Energy.begin(), fPhotocathode_Energy.end());
  std::reverse(fPhotocathode_Reflectance.begin(), fPhotocathode_Reflectance.end());
  std::reverse(fPhotocathode_Efficiency.begin(), fPhotocathode_Efficiency.end());
  std::reverse(fPhotocathode_SpecularLobe.begin(), fPhotocathode_SpecularLobe.end());
  std::reverse(fPhotocathode_SpecularSpike.begin(), fPhotocathode_SpecularSpike.end());
  std::reverse(fPhotocathode_Backscatter.begin(), fPhotocathode_Backscatter.end());
  std::reverse(fPhotocathode_Rindex.begin(), fPhotocathode_Rindex.end());

  myPropertiesFile.close();

  // Get the dimensions of the Lumirror sheath around the cryostat
  G4int indexMaxR = 0;
  fCryoSheathR[0] = 0;
  fCryoSheathR[1] = 0;
  fCryoSheathZ[0] = 0;
  fCryoSheathZ[1] = 0;

  G4int myNumPointsOuterCryo;
  G4double myOuterCryostatZ[30], myOuterCryostatR[30];
  DSIO::Get()->GetStreamDSCryostatProfile() >> myNumPointsOuterCryo;
  for (G4int i = 0; i < myNumPointsOuterCryo; i++) {
    DSIO::Get()->GetStreamDSCryostatProfile() >> myOuterCryostatZ[i] >> myOuterCryostatR[i];
    myOuterCryostatZ[i] += fCryostatShiftZ;
    if (i > 0 && myOuterCryostatZ[i] < myOuterCryostatZ[i - 1] && myOuterCryostatR[i] < myOuterCryostatR[i - 1] && myOuterCryostatR[i] > fCryoSheathR[0]) {
      indexMaxR = i;
      fCryoSheathR[0] = myOuterCryostatR[i];
    }
  }
  DSIO::Get()->CloseStreamDSCryostatProfile();
  fCryoSheathR[0] = myOuterCryostatR[indexMaxR];
  fCryoSheathR[1] = myOuterCryostatR[indexMaxR + 1];
  fCryoSheathZ[0] = myOuterCryostatZ[indexMaxR];
  fCryoSheathZ[1] = myOuterCryostatZ[indexMaxR + 1];
  DSLog(debugging) << "Cryostat Sheath Dimensions : " << endlog;
  DSLog(debugging) << "\tfCryoSheathR[0] = " << fCryoSheathR[0] << endlog;
  DSLog(debugging) << "\tfCryoSheathR[1] = " << fCryoSheathR[1] << endlog;
  DSLog(debugging) << "\tfCryoSheathZ[0] = " << fCryoSheathZ[0] << endlog;
  DSLog(debugging) << "\tfCryoSheathZ[1] = " << fCryoSheathZ[1] << endlog;
}

// LSV Optics
void DSParameters::ReadOpticsLSV() {

  G4String mys;
  //DSIO::Get()->GetStreamLogFile() << "########## LSV Optics  ###########" << endl;

  while (getline(DSIO::Get()->GetStreamLSVOptics(), mys)) {
    //DSIO::Get()->GetStreamLogFile() << mys << endl;

    if (!mys.find("ScintillationYield")) {
      fLSYield = atof(mys.substr(mys.find("=") + 1).c_str()) / MeV;
      // DSLog(trace) << " LS Scintillation Yield : " << fLSYield * MeV << "
      // /MeV" << endlog;
    } else if (!mys.find("BirksBeta")) {
      fLSBirksBeta = atof(mys.substr(mys.find("=") + 1).c_str()) * cm / MeV;
      DSLog(trace) << " LS Birks Constant for Betas: " << fLSBirksBeta * MeV / cm << " cm/MeV" << endlog;
    } else if (!mys.find("BirksAlpha")) {
      fLSBirksAlpha = atof(mys.substr(mys.find("=") + 1).c_str()) * cm / MeV;
      DSLog(trace) << " LS Birks Constant for Alphas: " << fLSBirksAlpha * MeV / cm << " cm/MeV" << endlog;
    } else if (!mys.find("BirksProton")) {
      fLSBirksProton = atof(mys.substr(mys.find("=") + 1).c_str()) * cm / MeV;
      DSLog(trace) << " LS Birks Constant for Protons: " << fLSBirksProton * MeV / cm << " cm/MeV" << endlog;
    } else if (!mys.find("Birks2Beta")) {
      fLSBirks2Beta = atof(mys.substr(mys.find("=") + 1).c_str()) * cm * cm / MeV / MeV;
      DSLog(trace) << " LS Second order Birks Constant for Betas: " << fLSBirks2Beta * MeV * MeV / cm / cm << " (cm/MeV)^2" << endlog;
    } else if (!mys.find("Birks2Alpha")) {
      fLSBirks2Alpha = atof(mys.substr(mys.find("=") + 1).c_str()) * cm * cm / MeV / MeV;
      DSLog(trace) << " LS Second order Birks Constant for Alphas: " << fLSBirks2Alpha * MeV * MeV / cm / cm << " (cm/MeV)^2" << endlog;
    } else if (!mys.find("Birks2Proton")) {
      fLSBirks2Proton = atof(mys.substr(mys.find("=") + 1).c_str()) * cm * cm / MeV / MeV;
      DSLog(trace) << " LS Second order Birks Constant for Protons: " << fLSBirks2Proton * MeV * MeV / cm / cm << " (cm/MeV)^2" << endlog;
    } else if (!mys.find("ResolutionScale")) {
      fLSResolutionScale = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << " LS Resolution Scale : " << fLSResolutionScale << ""
      // << endlog;
    } else if (!mys.find("BetaScintTime")) {
      fLSBetaScintTime.push_back(atof(mys.substr(mys.find("=") + 1).c_str()) * ns);
      // DSLog(trace) << " LS Beta Scintillation Time " <<
      // fLSBetaScintTime.size() << " : " << fLSBetaScintTime.back() / ns << "
      // ns" << endlog;
    } else if (!mys.find("BetaScintWeight")) {
      fLSBetaScintWeight.push_back(atof(mys.substr(mys.find("=") + 1).c_str()));
      // DSLog(trace) << " LS Beta Scintillation Weight " <<
      // fLSBetaScintWeight.size() << " : " << fLSBetaScintWeight.back() << ""
      // << endlog;
    } else if (!mys.find("AlphaScintTime")) {
      fLSAlphaScintTime.push_back(atof(mys.substr(mys.find("=") + 1).c_str()) * ns);
      // DSLog(trace) << " LS Alpha Scintillation Time " <<
      // fLSAlphaScintTime.size() << " : " << fLSAlphaScintTime.back() / ns << "
      // ns" << endlog;
    } else if (!mys.find("AlphaScintWeight")) {
      fLSAlphaScintWeight.push_back(atof(mys.substr(mys.find("=") + 1).c_str()));
      // DSLog(trace) << " LS Alpha Scintillation Weight " <<
      // fLSAlphaScintWeight.size() << " : " << fLSAlphaScintWeight.back() << ""
      // << endlog;
    } else if (!mys.find("FastTime")) {
      fLSFastTime = atof(mys.substr(mys.find("=") + 1).c_str()) * ns;
      // DSLog(trace) << " LS Fast Time : " << fLSFastTime / ns << " ns" <<
      // endlog;
    } else if (!mys.find("SlowTime")) {
      fLSSlowTime = atof(mys.substr(mys.find("=") + 1).c_str()) * ns;
      // DSLog(trace) << " LS Slow Time : " << fLSSlowTime / ns << " ns" <<
      // endlog;
    } else if (!mys.find("YieldRatio")) {
      fLSYieldRatio = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << " LS Yield Ratio : " << fLSYieldRatio << "" << endlog;
    } else if (!mys.find("WLSTime")) {
      fLSWLSTime = atof(mys.substr(mys.find("=") + 1).c_str()) * ns;
      // DSLog(trace) << " LS WLS Time : " << fLSWLSTime / ns << " ns" <<
      // endlog;
    } else if (!mys.find("WLSEfficiency")) {
      fLSWLSEfficiency = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << "LS WLS Efficiency : " << fLSWLSEfficiency << "" <<
      // endlog;
    }
    // else if (!mys.find("LumirrorRefl")) {
    // fLumirrorRefl = atof(mys.substr(mys.find("=")+1).c_str());
    // DSLog(trace) << "Lumirror Reflectivity : " << fLumirrorRefl << "" <<
    // endlog;
    // }
    else if (!mys.find("LumirrorEff")) {
      fLumirrorEff = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << "Lumirror Efficiency : " << fLumirrorEff << "" <<
      // endlog;
    } else if (!mys.find("LumirrorSpecLobe")) {
      fLumirrorSpecLobe = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << "Lumirror Specular Lobe : " << fLumirrorSpecLobe << ""
      // << endlog;
    } else if (!mys.find("LumirrorSpecSpike")) {
      fLumirrorSpecSpike = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << "Lumirror Specular Spike : " << fLumirrorSpecSpike <<
      // "" << endlog;
    } else if (!mys.find("LumirrorBackscat")) {
      fLumirrorBackscat = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << "Lumirror Backscattering : " << fLumirrorBackscat << ""
      // << endlog;
    } else if (!mys.find("LumirrorRindex")) {
      fLumirrorRindex = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << "Lumirror Refractive Index : " << fLumirrorRindex << ""
      // << endlog;
    } else if (!mys.find("EpssRefl")) {
      fEpssRefl = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << " EPSS (TPC cryostat) Reflectivity : " << fEpssRefl <<
      // "" << endlog;
    } else if (!mys.find("EpssEff")) {
      fEpssEff = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << "EPSS (TPC cryostat) Efficiency : " << fEpssEff << ""
      // << endlog;
    } else if (!mys.find("EpssSpecLobe")) {
      fEpssSpecLobe = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << " EPSS (TPC cryostat) Specular Lobe : " <<
      // fEpssSpecLobe << "" << endlog;
    } else if (!mys.find("EpssSpecSpike")) {
      fEpssSpecSpike = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << "EPSS (TPC cryostat) Specular Spike : " <<
      // fEpssSpecSpike << "" << endlog;
    } else if (!mys.find("EpssBackscat")) {
      fEpssBackscat = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << "EPSS (TPC cryostat) Backscattering : " <<
      // fEpssBackscat << "" << endlog;
    } else if (!mys.find("EpssRindex")) {
      fEpssRindex = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << "EPSS (TPC cryostat) Refractive Index : " <<
      // fEpssRindex << "" << endlog;
    } else if (!mys.find("VCathodeRefl")) {
      fVCathodeRefl = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << " PMT Back Reflectivity : " << fVCathodeRefl << "" <<
      // endlog;
    } else if (!mys.find("VCathodeEff")) {
      fVCathodeEff = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << " PMT Back Efficiency : " << fVCathodeEff << "" <<
      // endlog;
    } else if (!mys.find("VCathodeSpecLobe")) {
      fVCathodeSpecLobe = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << " PMT Back Specular Lobe : " << fVCathodeSpecLobe << ""
      // << endlog;
    } else if (!mys.find("VCathodeSpecSpike")) {
      fVCathodeSpecSpike = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << " PMT Back Specular Spike : " << fVCathodeSpecSpike <<
      // "" << endlog;
    } else if (!mys.find("VCathodeBackscat")) {
      fVCathodeBackscat = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << " PMT Back Backscattering : " << fVCathodeBackscat <<
      // "" << endlog;
    } else if (!mys.find("VCathodeRindex")) {
      fVCathodeRindex = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << " PMT Back Refractive Index : " << fVCathodeRindex <<
      // "" << endlog;
    } else if (!mys.find("PMTBackRefl")) {
      fPMTBackRefl = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << " PMT Back Reflectivity : " << fPMTBackRefl << "" <<
      // endlog;
    } else if (!mys.find("PMTBackEff")) {
      fPMTBackEff = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << "PMT Back Efficiency : " << fPMTBackEff << "" <<
      // endlog;
    } else if (!mys.find("PMTBackSpecLobe")) {
      fPMTBackSpecLobe = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << "PMT Back Specular Lobe : " << fPMTBackSpecLobe << ""
      // << endlog;
    } else if (!mys.find("PMTBackSpecSpike")) {
      fPMTBackSpecSpike = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << "PMT Back Specular Spike : " << fPMTBackSpecSpike << ""
      // << endlog;
    } else if (!mys.find("PMTBackBackscat")) {
      fPMTBackBackscat = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << "PMT Back Backscattering : " << fPMTBackBackscat << ""
      // << endlog;
    } else if (!mys.find("PMTBackRindex")) {
      fPMTBackRindex = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << "PMT Back Refractive Index : " << fPMTBackRindex << ""
      // << endlog;
    } else if (!mys.find("UnssRefl")) {
      fUnssRefl = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << "Untreated Stainless Steel Reflectivity : " <<
      // fUnssRefl << "" << endlog;
    } else if (!mys.find("UnssEff")) {
      fUnssEff = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << "Untreated Stainless Steel Efficiency : " << fUnssEff
      // << "" << endlog;
    } else if (!mys.find("UnssSpecLobe")) {
      fUnssSpecLobe = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << "Untreated Stainless Steel Specular Lobe : " <<
      // fUnssSpecLobe << "" << endlog;
    } else if (!mys.find("UnssSpecSpike")) {
      fUnssSpecSpike = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << "Untreated Stainless Steel Specular Spike : " <<
      // fUnssSpecSpike << "" << endlog;
    } else if (!mys.find("UnssBackscat")) {
      fUnssBackscat = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << "Untreated Stainless Steel Backscattering : " <<
      // fUnssBackscat << "" << endlog;
    } else if (!mys.find("UnssRindex")) {
      fUnssRindex = atof(mys.substr(mys.find("=") + 1).c_str());
      // DSLog(trace) << "Untreated Stainless Steel Refractive Index : " <<
      // fUnssRindex << "" << endlog;
    }
  }
  DSIO::Get()->CloseStreamLSVOptics();

  // Check that the sizes of fLS(..)ScintTime and fLS(..)ScintWeight vectors are
  // the same
  if (fLSBetaScintTime.size() != fLSBetaScintWeight.size())
    DSLog(error) << " Sizes of Beta Scintillation Time and Weight vectors are "
                    "not the same: "
                 << fLSBetaScintTime.size() << " != " << fLSBetaScintWeight.size() << endlog;
  if (fLSAlphaScintTime.size() != fLSAlphaScintWeight.size())
    DSLog(error) << " Sizes of Alpha Scintillation Time and Weight vectors are "
                    "not the same: "
                 << fLSAlphaScintTime.size() << " != " << fLSAlphaScintWeight.size() << endlog;

  // Get PC, PPO, and TMB attenuation lengths
  // WARNING: the wavelenght vectors should be equal in order to combine the
  // attenuation lengths
  G4double mywl, myvalue;
  // PC
  fPCAttEnergy.reserve(1000);
  fPCAttLength.reserve(1000);
  while (DSIO::Get()->GetStreamPCAttLength() >> mywl >> myvalue) {
    fPCAttEnergy.push_back(nmToEnergy(mywl));
    fPCAttLength.push_back(myvalue * m);
  }
  DSIO::Get()->CloseStreamPCAttLength();
  // reversing the vector in order to have low energy (high wavelength) at the
  // begining
  std::reverse(fPCAttEnergy.begin(), fPCAttEnergy.end());
  std::reverse(fPCAttLength.begin(), fPCAttLength.end());

  // PPO
  fPPOAttEnergy.reserve(1000);
  fPPOAttLength.reserve(1000);
  while (DSIO::Get()->GetStreamPPOAttLength() >> mywl >> myvalue) {
    fPPOAttEnergy.push_back(nmToEnergy(mywl));
    fPPOAttLength.push_back(myvalue * m);
  }
  DSIO::Get()->CloseStreamPPOAttLength();
  // reversing the vector in order to have low energy (high wavelength) at the
  // begining
  std::reverse(fPPOAttEnergy.begin(), fPPOAttEnergy.end());
  std::reverse(fPPOAttLength.begin(), fPPOAttLength.end());
  if (fPCAttEnergy.size() != fPPOAttEnergy.size())
    DSLog(warning) << "WARNING: PC and PPO attenuation length vector differs "
                      "in size. The vectors should be equal in "
                      "size to combine them."
                   << endlog;

  // TMB
  fTMBAttEnergy.reserve(1000);
  fTMBAttLength.reserve(1000);
  while (DSIO::Get()->GetStreamTMBAttLength() >> mywl >> myvalue) {
    fTMBAttEnergy.push_back(nmToEnergy(mywl));
    fTMBAttLength.push_back(myvalue * m);
  }
  DSIO::Get()->CloseStreamTMBAttLength();
  // reversing the vector in order to have low energy (high wavelength) at the
  // begining
  std::reverse(fTMBAttEnergy.begin(), fTMBAttEnergy.end());
  std::reverse(fTMBAttLength.begin(), fTMBAttLength.end());
  if (fPCAttEnergy.size() != fTMBAttEnergy.size())
    DSLog(warning) << "WARNING: PC and TMB attenuation length vector differs "
                      "in size. The vectors should be equal in "
                      "size to combine them."
                   << endlog;

  // PC Refractive Index
  fPCRefrIndex.reserve(1000);
  fPCRefrEnergy.reserve(1000);
  while (DSIO::Get()->GetStreamPCRefractionIndex() >> mywl >> myvalue) {
    fPCRefrEnergy.push_back(nmToEnergy(mywl));
    fPCRefrIndex.push_back(myvalue);
  }
  DSIO::Get()->CloseStreamPCRefractionIndex();
  std::reverse(fPCRefrEnergy.begin(), fPCRefrEnergy.end());
  std::reverse(fPCRefrIndex.begin(), fPCRefrIndex.end());

  // Lumirror Reflectivity
  fLumirrorReflec.reserve(1000);
  fLumirrorEnergy.reserve(1000);
  while (DSIO::Get()->GetStreamLumirrorReflectivity() >> mywl >> myvalue) {
    fLumirrorEnergy.push_back(nmToEnergy(mywl));
    fLumirrorReflec.push_back(myvalue / 100.);
  }
  DSIO::Get()->CloseStreamLumirrorReflectivity();
  // NO reverse: Lumirror wavelenghts are sorted in decreasing order in the file
  // std::reverse(fLumirrorReflec.begin(), fLumirrorReflec.end());
  // std::reverse(fLumirrorEnergy.begin(), fLumirrorEnergy.end());

  // Veto PMT QEs and reflectivity
  vector<G4double> myqe;
  myqe.reserve(110);
  fVPmtMaxQe = 0.;
  while (DSIO::Get()->GetStreamDS50VetoPmtQE() >> mywl >> myvalue) {
    // *temporary patch* for 20k PMT 20i type: rescale the 20i QE by the ratio
    // of the max QE (20i/8i) .SD.
    if (DSStorage::Get()->GetIsDS20kLSV()) myvalue *= (26. / 37.);
    // end of patch
    myqe.push_back(myvalue);
    fVPmtMaxQe = std::max(myvalue, fVPmtMaxQe);
    fVPmtQeWl.push_back(mywl * nm);
    fVCathodeEnergy.push_back(nmToEnergy(mywl));
    // DSLog(debugging) << "VETO PMT QE : WL = " << wavelength << " nm\tQE = "
    // << qe << "\tQE' = " << qe_adjusted << "\tr = " << refl << endlog;
  }
  DSIO::Get()->CloseStreamDS50VetoPmtQE();

  fVPmtMaxQe_adjusted = fVPmtMaxQe / (1 - fVCathodeRefl * (1 - fVPmtMaxQe));

  for (unsigned i = 0; i < myqe.size(); i++) {
    G4double qe = myqe[i];
    G4double qe_adjusted = qe / (1 - fVCathodeRefl * (1 - qe));
    qe /= fVPmtMaxQe;
    qe_adjusted /= fVPmtMaxQe_adjusted;
    G4double refl = (1 - qe) * fVCathodeRefl;
    refl *= fVPmtMaxQe_adjusted;

    fVPmtQe.push_back(qe_adjusted);
    fVCathodeReflec.push_back(refl);
  }

  // Revert all vectors: energy should be in increasing order
  std::reverse(fVCathodeEnergy.begin(), fVCathodeEnergy.end());
  std::reverse(fVCathodeReflec.begin(), fVCathodeReflec.end());

  DSLog(routine) << "LSV Optics Properties loaded" << endlog;
}

// Method to convert a wavelength in nm to G4 energy unit
inline G4double DSParameters::nmToEnergy(G4double wavelength) {
  const G4double hc = h_Planck * c_light / nm;
  return (hc / wavelength);
}

/////////////////////////////////////////////////////////
//// Method returning the SiPM PDE , given the WL in nm
/////////////////////////////////////////////////////////

void DSParameters::AdjustSiPMPDE () { 
  const int myPDEsize = 16;
  
  for (int i = 0; i < myPDEsize; i++) {
    myPDE[i] *= 0.01;     // convert to fraction
    myPDE[i] *= fSiPMMaxPDE / 0.5277;  // scale so maximum PDE is 0.40

    if (DSStorage::Get()->GetFastSimulation() ) { 
      myPDE[i] /= fSiPMMaxPDE ; 
    }

    myPDEWL[i] *= nm;
  }
  IsSipmPDEAdjusted = true ; 
}


G4double DSParameters::GetSiPMPDE(G4double myPhotonWL) {
  
  if (!IsSipmPDEAdjusted) DSParameters::AdjustSiPMPDE () ; 
  const int myPDEsize = 16;

  if (myPhotonWL < myPDEWL[0]) return 0;
  if (myPhotonWL > myPDEWL[myPDEsize - 1]) return 0;

  G4int myIndex = 0;
  G4double _WL0 = myPDEWL[0];
  G4double _WL1 = myPDEWL[1];
  for (int it = 0; it < myPDEsize; ++it) {
    _WL1 = myPDEWL[it];
    if (myPDEWL[it] > myPhotonWL) break;
    _WL0 = myPDEWL[it];
    myIndex++;
  }

  const G4double _QE0 = myPDE[myIndex - 1];
  const G4double _QE1 = myPDE[myIndex];
  // interpolation
  const G4double _m = (_QE0 - _QE1) / (_WL0 - _WL1);
  const G4double _q = _QE1 - _m * _WL1;
  // cout << endl << myPhotonWL/nm << " " << _QE0 << " " << _QE1 << " " << _m*myPhotonWL + _q << endl << endl;
  return _m * myPhotonWL + _q;
}

/////////////////////////////////////////////////////////
//// Method returning the TPC QE , given the WL in nm
/////////////////////////////////////////////////////////
G4double DSParameters::GetTPCQE(G4double myPhotonWL) const {

  //  if( myPhotonWL < fPmtQeWl[1] )       return fPmtQe[1];
  if (myPhotonWL < fPmtQeWl[1]) return 0;  // return 0 for the QE if WL<200 nm
  if (myPhotonWL > fPmtQeWl.back()) return fPmtQe.back();

  G4int myIndex = 0;
  G4double _WL0 = 0;
  G4double _WL1 = 0;

  for (std::vector<G4double>::const_iterator it = fPmtQeWl.begin(); it != fPmtQeWl.end(); ++it) {
    _WL1 = *it;
    if (*it > myPhotonWL) break;
    _WL0 = *it;
    myIndex++;
  }

  const G4double _QE0 = fPmtQe[myIndex - 1] * fPmtMaxQe;
  const G4double _QE1 = fPmtQe[myIndex] * fPmtMaxQe;
  // interpolation
  const G4double _m = (_QE0 - _QE1) / (_WL0 - _WL1);
  const G4double _q = _QE1 - _m * _WL1;

  return _m * myPhotonWL + _q;
}

/////////////////////////////////////////////////////////
//// Method returning the Neutron Veto QE , given the WL in nm
/////////////////////////////////////////////////////////
G4double DSParameters::GetNVQE(G4double myPhotonWL) {

  G4double myQE = interpolate(&fVPmtQe[0], &fVPmtQeWl[0], myPhotonWL, (int)fVPmtQe.size());
  // *temporary patch* for 20k PMT 20i type: rescale the 8i QE by the ratio of
  // the max QE (20i/8i) .SD.
  // const static G4int myPMTtype = DSStorage::Get()->GetDS20kLSVPMT();
  // if (myPMTtype==2) myQE*=(26./37.);
  // end of patch

  // DSLog(development) << "WL = " << myPhotonWL << "\tQE = " << myQE << endlog;

  return myQE;
}

int DSParameters::binSearch(const G4double* list, const G4double& x, int minIndex, int maxIndex) {
  if (x < list[minIndex + 1] && x >= list[maxIndex - 1]) return minIndex;
  if (x < list[(minIndex + maxIndex) / 2]) return binSearch(list, x, minIndex, (minIndex + maxIndex) / 2);
  return binSearch(list, x, (minIndex + maxIndex) / 2, maxIndex);
}

int DSParameters::search(const G4double* list, const G4double& x, int listSize) {
  if (x <= list[0]) return 0;
  if (x >= list[listSize - 1]) return listSize - 1;
  return binSearch(list, x, 0, listSize - 1);
}

// quadratically interpolate a point in an array
G4double DSParameters::interpolate(const G4double* ys, const G4double* xs, const G4double& x, int size) {
  // Find the greatest element <= x
  int myIndex = search(xs, x, size);
  // Check if it is below or above the range. If so, return endpoint
  if (myIndex == 0) return ys[0];
  if (myIndex == size - 1) return ys[size - 1];
  // Check if x falls on a known point and no interpolation is needed
  if (xs[myIndex] == x) return ys[myIndex];

  const G4double denom1 = (xs[myIndex - 1] - xs[myIndex]) * (xs[myIndex - 1] - xs[myIndex + 1]);
  const G4double denom2 = (xs[myIndex - 1] - xs[myIndex]) * (xs[myIndex] - xs[myIndex + 1]);
  const G4double denom3 = (xs[myIndex - 1] - xs[myIndex + 1]) * (xs[myIndex] - xs[myIndex + 1]);

  const G4double b = ys[myIndex - 1] / denom1 - ys[myIndex] / denom2 + ys[myIndex + 1] / denom3;
  const G4double c = -(xs[myIndex] + xs[myIndex + 1]) * ys[myIndex - 1] / denom1 + (xs[myIndex - 1] + xs[myIndex + 1]) * ys[myIndex] / denom2 - (xs[myIndex - 1] + xs[myIndex]) * ys[myIndex + 1] / denom3;
  const G4double d = xs[myIndex] * xs[myIndex + 1] * ys[myIndex - 1] / denom1 - xs[myIndex - 1] * xs[myIndex + 1] * ys[myIndex] / denom2 + xs[myIndex - 1] * xs[myIndex] * ys[myIndex + 1] / denom3;

  return b * x * x + c * x + d;
}

std::vector<G4double> DSParameters::ConvertTo4DoubleVector(const char* st) {
  G4double v1;
  G4double v2;
  G4double v3;
  G4double v4;
  std::istringstream is(st);
  is >> v1 >> v2 >> v3 >> v4;
  std::vector<G4double> tmp;
  tmp.push_back(v1);
  tmp.push_back(v2);
  tmp.push_back(v3);
  tmp.push_back(v4);
  return tmp;
}

/*
 * $Log: DSParameters.cc,v $
 * Revision 1.11  2016/03/03 23:17:35  cz2
 * Added two parameters related to TPB defects and model for S2 radial
 * correction factor.
 *
 * Revision 1.10  2015/01/17 11:31:53  pagnes
 * PAr model added form optical tuning
 *
 * Revision 1.9  2015/01/07 16:45:53  pagnes
 * changed veto optical properties format from arrays to vectors
 *
 * Revision 1.8  2014/12/22 22:42:04  jbrodsky
 * Fixed bug in filling veto PMT locations from data file. The loop was going
 * one over the limit of the relevant arrays (trying to fill entry 110 in an
 * array of size 110) which caused bad codebehavior and segfaults.
 *
 * Revision 1.7  2014/11/20 17:08:55  dfranco
 * QE of PMTs = 0 if wavelength<200 nm
 *
 * Revision 1.6  2014/11/13 16:47:04  dfranco
 * removed variables which were creating conflicts with the previous version of
 * g4ds10
 *
 * Revision 1.5  2014/11/05 15:47:13  pagnes
 * temporary optics tuning
 *
 * Revision 1.4  2014/10/13 18:43:57  swesterd
 * fixed veto PMT positions
 *
 * Revision 1.3  2014/07/16 08:23:04  pagnes
 * QE scaling to 1.0 added (/ds/manager/fast_simulation xxx)
 *
 * Revision 1.2  2014/06/03 13:31:35  meyers
 * Migrate TPC grid and ITO optics updates to g4ds10
 *
 * Revision 1.1  2014/05/07 12:21:04  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.21  2014/03/19 16:37:26  dfranco
 * update external configuration reader for optics tuning
 *
 * Revision 1.20  2014/03/11 16:49:56  dfranco
 * read optical parameters from external file (DSOptics.dat)
 *
 * Revision 1.19  2013/08/05 12:25:13  swesterd
 * Added G4OpBoundaryProcess.hh so that light can pass through some optical
 * boundaries
 *
 * Revision 1.18  2013/08/05 03:13:52  swesterd
 * some fine tuning of bscint and veto parameters
 *
 * Revision 1.17  2013/06/24 13:43:29  dfranco
 * cout removed
 *
 * Revision 1.16  2013/06/24 13:05:55  dfranco
 * TPC QE values were filled twice: once in the standard way, the second
 * deconvoluting the reflections. The second has been commented
 *
 * Revision 1.15  2013/06/22 07:21:21  dfranco
 * Fixed a bug in the photoelectron absorption in DS50. Back to the previous QE
 * method
 *
 * Revision 1.14  2013/06/19 18:35:28  swesterd
 * added DSScintCelll and made tpc PMTs' QE and reflections work like veto PMTs
 *
 * Revision 1.13  2013/06/05 23:03:32  swesterd
 * moved optical boundary MPTs to DSMaterial and gave the trunks optical
 * boundary properties consistent with untreated stainless steel
 *
 * Revision 1.12  2013/06/04 23:38:49  swesterd
 * Changed the length of the Lumirror reflectance array to be an adjustable
 * size, set it to 60 elements
 *
 * Revision 1.11  2013/06/04 14:11:36  dfranco
 * Added a function returning the TPC QE in DSParameters, applied in the
 * tracking action
 *
 * Revision 1.10  2013/05/27 23:59:02  swesterd
 * added a (currently commented out) Lumirror sheath to the cryostat and
 * introduced DSOpBoundaryProcess to try to figure out why the boundaries are
 * being screwy, with some edits so that it can handle constant and vector
 * properties with freaking out
 *
 * Revision 1.9  2013/05/25 07:58:23  swesterd
 * Got the veto PMT optical boundaries all working along with photocathode
 * optical properties, added PMT quantum efficiency to DSTrackingAction, and
 * added a function to DSTrackingAction that locates and quadratically
 * interpolates points in data, for getting useful QEs
 *
 * Revision 1.8  2013/05/07 23:06:30  swesterd
 * added optical boundaries and Lumirror in the veto
 *
 * Revision 1.7  2013/05/06 14:59:52  perassos
 * Updates on the TPC surface properties and geometry
 *
 * Revision 1.6  2013/05/01 08:20:24  swesterd
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
