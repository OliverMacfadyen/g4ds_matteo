/* -------------------------------------- */
/* ArDM geometry based on DSDart detector*/
/* Olivier Dadoun */
/* March 2017 */
/* -------------------------------------- */

#ifndef DSDetectorArDM_H
#define DSDetectorArDM_H

#include "globals.hh"
#include "vector"

#include "G4LogicalVolume.hh"
#include "G4Tubs.hh"
#include "G4VPhysicalVolume.hh"

#include "G4LogicalBorderSurface.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4OpticalSurface.hh"

//#include "ArDM_Analysis.hh"

struct ArDM {
  inline std::vector<G4TwoVector> setPMTVector2D();
  inline std::vector<G4ThreeVector> setPMTVector(const char* top_or_bottom = "bottom");
  inline G4double getWLSThickness(G4double convEff);

  G4double extended_asin(G4double ang) ; // redefine TMath::Asin
  void Init();


  double HANDEDNESS;

  double PI;
  double BOLTZMAN_K;  // joule /kelvin              ; //boltzmann constant;

  int NPMT;
  int NTOPPMT;
  int NBOTTOMPMT;
  double ELECTRIC_FIELD_STRENGTH;
  double EXTRACTION_FIELD_STRENGTH;
  double DIELECTRIC_CONSTANT_GAR;
  double DIELECTRIC_CONSTANT_LAR;
  double PHT_ENE_THRESHOLD;  // if photon energy is larger than this threshold
                             // --> it won't be detected !;
  double BOILING_POINT_ARGON;

  double COMPONENT_RATIO_LAR_ELECTRON_LIKE;
  double COMPONENT_RATIO_LAR_NEUTRON_LIKE;
  double COMPONENT_RATIO_GAR_ELECTRON_LIKE;
  double COMPONENT_RATIO_GAR_NEUTRON_LIKE;

  double PMTCONV;
  double REFLCONV;

  double WORLD_HALF_SIZE;
  double LID_RADIUS;
  double LID_HALF_HEIGHT;
  double TANK_HALF_HEIGHT;
  double TANK_CYLINDER_INNER_RADIUS;
  double TANK_CYLINDER_THICKNESS;
  double TANK_CYLINDER_OUTER_RADIUS;
  double TANK_BTM_PART_R1010_ARC_INNER_RADIUS;
  double TANK_BTM_PART_R1010_ARC_OUTER_RADIUS;
  double TANK_BTM_PART_R101_ARC_INNER_RADIUS;
  double TANK_BTM_PART_R101_ARC_OUTER_RADIUS;
  double TANK_BTM_PART_R1010_ARC_CHORD;
  double TANK_BTM_PART_R1010_ARC_OPENING_ANGLE;
  double TANK_BTM_PART_R101_ARC_OPENING_ANGLE;
  double DISTANCE_R1010_ARC_CENTER_TO_LOWER_EDGE_OF_TANK_CYLINDER;
  double DISTANCE_LOWER_EDGE_OF_TANK_CYLINDER_TO_TANK_LOWERMOST_POINT;
  double TANK_CYLINDER_HALF_HEIGHT;
  double TANK_BTM_PART_R101_ARC_HEIGHT;
  double TANK_CYLINDER_POS_Z;
  double TANK_BTM_PART_R1010_Z;
  double TANK_BTM_PART_R1010_Z_RELATIVE_TO_TANK_CYLINDER;
  double DISTANCE_R1010_ARC_CENTER_TO_R101_ARC_CENTER;
  double DISTANCE_TANK_LOWERMOST_POINT_TO_R101_ARC_CENTER;
  double DISTANCE_TANK_CYLINDER_CENTER_TO_R101_ARC_CENTER;
  double TANK_BTM_PART_R101_Z;
  double TANK_BTM_PART_R101_Z_RELATIVE_TO_TANK_CYLINDER;
  double TOP_FLANGE_INNER_RADIUS;
  double TOP_FLANGE_OUTER_RADIUS;
  double TOP_FLANGE_HALF_HEIGHT;
  double TOP_FLANGE_COMPENSATION_FACTOR_FOR_OMITTED_PILLARS;
  double TOP_FLANGE_HALF_HEIGHT_EFFECTIVE;
  double TOP_FLANGE_POS_X;
  double TOP_FLANGE_POS_Y;
  double TOP_FLANGE_POS_Z;
  double DISTANCE_TOP_FLANGE_TOP_PMT_CENTER;
  double DISTANCE_UPPER_EDGE_FIELD_SHAPER_LIQUID_SURFACE;
  double REFLECTOR_HALF_HEIGHT;
  double DISTANCE_TOP_FLANGE_UPPER_EDGE_FIELD_SHAPER;
  double DISTANCE_TOP_FLANGE_LIQUID_SURFACE;
  double LIQUID_SURFACE_POS_Z;
  double DISTANCE_TOP_FLANGE_REFLECTOR_CENTER;
  double DISTANCE_REFLECTOR_LOWER_EDGE_PROTECTION_GRID;  //<-- check it ! ;
  double DISTANCE_PROTECTION_GRID_TIP_OF_BTM_PMT;        //<-- check it !;
  double GAR_COLUMN_HALF_HEIGHT;
  double LAR_COLUMN_HALF_HEIGHT;
  double LAR_COLUMN_Z;
  double GAR_COLUMN_Z;
  double LARCOL_TOTAL_HEIGHT;
  double LARCOL_CYLINDER_HALFHEIGHT;
  double LARCOL_CYLINDER_OUTER_RADIUS;
  double LARCOL_CYLINDER_POS_Z;
  double WLS_MEAN_ABSORPTION_LENGTH;
  double WLS_THICKNESS_100_PERCENT_CONVERSION_EFFICIENCY;
  double WLS_RINGSEC_THICKNESS;
  double WLS_RINGSEC_INNER_RADIUS;
  double WLS_RINGSEC_OUTER_RADIUS;
  double WLS_RINGSEC_HALF_HEIGHT;
  double WLS_RINGSEC_START_PHI;
  double WLS_RINGSEC_DELTA_PHI;
  double WLS_RINGSEC_RADIUS;
  double WLS_RINGSEC_POS_Z;
  double WLS_RINGSEC_POS_Z_LARCOL;
  double WLS_LINSEC_HALF_X;
  double WLS_LINSEC_HALF_Y;
  double WLS_LINSEC_HALF_Z;
  double WLS_LINSEC_HALF_THICKNESS;
  double WLS_LINSEC_POS_X;
  double WLS_LINSEC_POS_Y;
  double WLS_LINSEC_POS_Z;
  double WLS_SUPPORT_RINGSEC_INNER_RADIUS;
  double WLS_SUPPORT_RINGSEC_THICKNESS;
  double WLS_SUPPORT_RINGSEC_OUTER_RADIUS;
  double WLS_SUPPORT_RINGSEC_HALF_HEIGHT;
  double WLS_SUPPORT_RINGSEC_START_PHI;
  double WLS_SUPPORT_RINGSEC_DELTA_PHI;
  double WLS_SUPPORT_RINGSEC_RADIUS;
  double WLS_SUPPORT_RINGSEC_POS_Z;
  double WLS_SUPPORT_LINSEC_HALF_X;
  double WLS_SUPPORT_LINSEC_HALF_Y;
  double WLS_SUPPORT_LINSEC_HALF_Z;
  double WLS_SUPPORT_LINSEC_POS_X;
  double WLS_SUPPORT_LINSEC_POS_Y;
  double WLS_SUPPORT_LINSEC_POS_Z;
  double WLS_UPPER_EDGE_POS_Z;
  double WLS_LOWER_EDGE_POS_Z;
  double WLS_RINGSEC_HALF_HEIGHT_GAR;
  double WLS_RINGSEC_GAR_POS_Z;
  double WLS_RINGSEC_GAR_POS_Z_GARCOL;
  double WLS_LINSEC_HALF_HEIGHT_GAR;
  double WLS_LINSEC_GAR_POS_X;
  double WLS_LINSEC_GAR_POS_Y;
  double WLS_LINSEC_GAR_POS_Z;
  double WLS_LINSEC_GAR_POS_Z_GARCOL;
  double WLS_SUPPORT_RINGSEC_HALF_HEIGHT_GAR;
  double WLS_SUPPORT_RINGSEC_GARCOL_POS_Z;
  double WLS_SUPPORT_RINGSEC_GARCOL_POS_Z_GARCOL;
  double WLS_SUPPORT_LINSEC_HALF_HEIGHT_GAR;
  double WLS_SUPPORT_LINSEC_GAR_POS_X;
  double WLS_SUPPORT_LINSEC_GAR_POS_Y;
  double WLS_SUPPORT_LINSEC_GAR_POS_Z;
  double WLS_SUPPORT_LINSEC_GAR_POS_Z_GARCOL;
  double WLS_RINGSEC_HALF_HEIGHT_LAR;
  double WLS_RINGSEC_LAR_POS_Z;
  double WLS_RINGSEC_LAR_POS_Z_LARCOL;
  double WLS_LINSEC_HALF_HEIGHT_LAR;
  double WLS_LINSEC_LAR_POS_X;
  double WLS_LINSEC_LAR_POS_Y;
  double WLS_LINSEC_LAR_POS_Z;
  double WLS_LINSEC_LAR_POS_Z_LARCOL;
  double WLS_SUPPORT_RINGSEC_HALF_HEIGHT_LAR;
  double WLS_SUPPORT_RINGSEC_LARCOL_POS_Z;
  double WLS_SUPPORT_RINGSEC_LARCOL_POS_Z_LARCOL;
  double WLS_SUPPORT_LINSEC_HALF_HEIGHT_LAR;
  double WLS_SUPPORT_LINSEC_LAR_POS_X;
  double WLS_SUPPORT_LINSEC_LAR_POS_Y;
  double WLS_SUPPORT_LINSEC_LAR_POS_Z;
  double WLS_SUPPORT_LINSEC_LAR_POS_Z_LARCOL;
  double PMT_INNER_RADIUS;
  double PMT_THICK;
  double PMT_OUTER_RADIUS;  // thickness of PMT = 1mm ??;
  double PMT_ACTIVE_RANGE;
  double PMT_SPHERICAL_PART_INNER_RADIUS;
  double PMT_SPHERICAL_PART_OUTER_RADIUS;
  double PMT_SPHERICAL_PART_OPENING_ANGLE;
  double PMT_MIDDLE_CYLINDER_INNER_RADIUS;
  double PMT_MIDDLE_CYLINDER_OUTER_RADIUS;
  double PMT_MIDDLE_CYLINDER_HALF_HEIGHT;
  double PMT_SPHERICAL_PART_HEIGHT;
  double PMT_BTM_TUBE_INNER_RADIUS;
  double PMT_BTM_TUBE_OUTER_RADIUS;
  double PMT_BTM_TUBE_HALF_HEIGHT;
  double PMT_BTM_SPHERE_HOLE_OPENING_ANGLE;
  double PMT_DISTANCE_PMT_CENTER_TO_SPHERE_CENTER;
  double PMT_DISTANCE_PMT_CENTER_TO_TOP_SPHERE_CENTER;
  double PMT_DISTANCE_PMT_CENTER_TO_BTM_SPHERE_CENTER;
  double PMT_DISTANCE_PMT_CENTER_TO_BTM_TUBE_CENTER;
  double DISTANCE_TOP_FLANGE_UPPER_EDGE_OF_TOP_PMT_TUBE;
  double TOP_PMT_Z;
  double TOP_PMT_Z_GARCOL;
  double BTM_PMT_Z;
  double BTM_PMT_Z_LARCOL;
  double TOP_PMT_CENTER_OF_TOP_SPHERE;
  double TOP_PMT_CENTER_OF_TOP_SPHERE_GARCOL;
  double BTM_PMT_CENTER_OF_TOP_SPHERE;
  double BTM_PMT_CENTER_OF_TOP_SPHERE_LARCOL;
  double PMT_BASE_RADIUS;
  double PMT_BASE_HALF_THICKNESS;
  double PMT_DISTANCE_LOWER_EDGE_OF_PMT_TUBE_TO_PMT_BASE_CENTER;
  double PMT_DISTANCE_PMT_CENTER_TO_PMT_BASE;
  double PMT_BASE_Z_TOP;
  double PMT_BASE_Z_TOP_GARCOL;
  double PMT_BASE_Z_BTM;
  double PMT_BASE_Z_BTM_LARCOL;
  double PMT_ELECTRODE_OUTER_RADIUS;
  double PMT_ELECTRODE_HALF_HEIGHT;
  double PMT_ELECTRODE_VOLUME;  // weight/density, material: stainless steel;
  double PMT_ELECTRODE_INNER_RADIUS;
  double PMT_ELECTRODE_Z_TOP;
  double PMT_ELECTRODE_Z_TOP_GARCOL;
  double PMT_ELECTRODE_Z_BTM;
  double PMT_ELECTRODE_Z_BTM_LARCOL;
  double HV_RESISTOR_BAR_INNER_RADIUS;
  double HV_RESISTOR_BAR_OUTER_RADIUS;
  double HV_RESISTOR_BAR_HALF_HEIGHT;
  double DISTANCE_BETWEEN_HV_RESISTOR_BARS;
  double HV_RESISTOR_BAR_1_POS_X;
  double HV_RESISTOR_BAR_1_POS_Y;
  double HV_RESISTOR_BAR_1_POS_Z;
  double HV_RESISTOR_BAR_1_POS_Z_LARCOL;
  double HV_RESISTOR_BAR_2_POS_X;
  double HV_RESISTOR_BAR_2_POS_Y;
  double HV_RESISTOR_BAR_2_POS_Z;
  double HV_RESISTOR_BAR_2_POS_Z_LARCOL;
  double PMT_COATING_THICKNESS;
  double PMT_COATING_INNER_RADIUS;
  double PMT_COATING_OUTER_RADIUS;
  double TOP_PMT_COATING_Z;
  double TOP_PMT_COATING_Z_GARCOL;
  double BTM_PMT_COATING_Z;
  double BTM_PMT_COATING_Z_LARCOL;
  double PMT_CATHODE_OUTER_RADIUS;
  double PMT_CATHODE_THICKNESS;
  double PMT_CATHODE_INNER_RADIUS;
  double PMT_CATHODE_ACTIVE_RANGE;
  double PMT_ACTIVE_RANGE_RADIUS;  //);
  double TOP_PMT_CATHODE_Z;
  double TOP_PMT_CATHODE_Z_GARCOL;
  double BTM_PMT_CATHODE_Z;
  double BTM_PMT_CATHODE_Z_LARCOL;
  double PMT_SUPPORT_RADIUS;
  double PMT_SUPPORT_HALF_HEIGHT;
  double TOP_PMT_SUPPORT_POS_Z;
  double TOP_PMT_SUPPORT_POS_Z_GARCOL;
  double BTM_PMT_SUPPORT_POS_Z;
  double BTM_PMT_SUPPORT_POS_Z_LARCOL;
  double TOP_PMT_SUPPORT_POS_LOWER_EDGE;
  double BTM_PMT_SUPPORT_POS_UPPER_EDGE;
  double PMT_SUPPORT_CONVERSION_EFFICIENCY;
  double PMT_SUPPORT_COATING_HALF_HEIGHT;
  double TOP_PMT_SUPPORT_COATING_POS_Z;
  double TOP_PMT_SUPPORT_COATING_POS_Z_GARCOL;
  double BTM_PMT_SUPPORT_COATING_POS_Z;
  double BTM_PMT_SUPPORT_COATING_POS_Z_LARCOL;
  double CATHODE_POS_Z;
  double CATHODE_POS_Z_LARCOL;
  double CATHODE_PLATE_INNER_RADIUS;
  double CATHODE_PLATE_OUTER_RADIUS;
  double CATHODE_PLATE_HALF_THICKNESS;
  double CATHODE_WIRE_INNER_RADIUS;
  double CATHODE_WIRE_OUTER_RADIUS;
  double CATHODE_WIRE_PITCH;
  double PMMA_PLATE_OUTER_RADIUS;
  double PMMA_PLATE_HALF_THICKNESS;
  double ITO_HALF_THICKNESS;
  double BTM_PROTECTION_GRID_PLATE_INNER_RADIUS;
  double BTM_PROTECTION_GRID_PLATE_OUTER_RADIUS;
  double BTM_PROTECTION_GRID_PLATE_HALF_THICKNESS;
  double DISTANCE_LOWER_EDGE_CATHODE_PLATE_TO_UPPER_EDGE_BTM_PROTECTION_GRID;
  double BTM_PROTECTION_GRID_WIRE_INNER_RADIUS;
  double BTM_PROTECTION_GRID_WIRE_OUTER_RADIUS;
  double BTM_PROTECTION_GRID_WIRE_PITCH;
  double BTM_PROTECTION_GRID_POS_Z;
  double BTM_PROTECTION_GRID_POS_Z_LARCOL;
  double SIDE_REFLECTOR_THICKNESS;
  double BTM_SIDE_REFLECTOR_HALF_HEIGHT;
  double BTM_SIDE_REFLECTOR_START_PHI;
  double BTM_SIDE_REFLECTOR_DELTA_PHI;
  double BTM_SIDE_REFLECTOR_POS_Z;         // relative to tank !!;
  double BTM_SIDE_REFLECTOR_POS_Z_LARCOL;  // relative to LArColumn !;
  double BTM_SIDE_REFLECTOR_INNER_RADIUS_UPPER_EDGE;
  double BTM_SIDE_REFLECTOR_OUTER_RADIUS_UPPER_EDGE;
  double BTM_SIDE_REFLECTOR_INNER_RADIUS_LOWER_EDGE;
  double BTM_SIDE_REFLECTOR_OUTER_RADIUS_LOWER_EDGE;
  double BTM_SIDE_REFLECTOR_HALF_HEIGHT_GASTEST;
  double BTM_SIDE_REFLECTOR_POS_Z_GASTEST;
  double BTM_SIDE_REFLECTOR_POS_Z_GASTEST_LARCOL;
  double BTM_SIDE_REFLECTOR_CYNLIDRICAL_PART_INNER_RADIUS;
  double BTM_SIDE_REFLECTOR_CYNLIDRICAL_PART_OUTER_RADIUS;
  double BTM_SIDE_REFLECTOR_CYNLIDRICAL_PART_HALF_HEIGHT;
  double BTM_SIDE_REFLECTOR_CYNLIDRICAL_PART_POS_Z;
  double BTM_SIDE_REFLECTOR_CYNLIDRICAL_PART_POS_Z_LARCOL;
  double TOP_SIDE_REFLECTOR_INNER_RADIUS_LOWER_EDGE;
  double TOP_SIDE_REFLECTOR_OUTER_RADIUS_LOWER_EDGE;
  double TOP_SIDE_REFLECTOR_INNER_RADIUS_UPPER_EDGE;
  double TOP_SIDE_REFLECTOR_OUTER_RADIUS_UPPER_EDGE;
  double TOP_SIDE_REFLECTOR_HALF_HEIGHT;
  double TOP_SIDE_REFLECTOR_START_PHI;
  double TOP_SIDE_REFLECTOR_DELTA_PHI;
  double TOP_SIDE_REFLECTOR_POS_Z;
  double TOP_SIDE_REFLECTOR_POS_Z_GARCOL;
  double TOP_SIDE_REFLECTOR_HALF_HEIGHT_IN_GARCOL;
  double TOP_SIDE_REFLECTOR_HALF_HEIGHT_IN_LARCOL;
  double TOP_SIDE_REFLECTOR_COTAN_TILTING_ANGLE;
  double TOP_SIDE_REFLECTOR_INNER_RADIUS_UPPER_EDGE_IN_GARCOL;
  double TOP_SIDE_REFLECTOR_OUTER_RADIUS_UPPER_EDGE_IN_GARCOL;
  double TOP_SIDE_REFLECTOR_INNER_RADIUS_LOWER_EDGE_IN_GARCOL;
  double TOP_SIDE_REFLECTOR_OUTER_RADIUS_LOWER_EDGE_IN_GARCOL;
  double TOP_SIDE_REFLECTOR_IN_GARCOL_POS_Z;
  double TOP_SIDE_REFLECTOR_IN_GARCOL_POS_Z_GARCOL;
  double TOP_SIDE_REFLECTOR_INNER_RADIUS_UPPER_EDGE_IN_LARCOL;
  double TOP_SIDE_REFLECTOR_OUTER_RADIUS_UPPER_EDGE_IN_LARCOL;
  double TOP_SIDE_REFLECTOR_INNER_RADIUS_LOWER_EDGE_IN_LARCOL;
  double TOP_SIDE_REFLECTOR_OUTER_RADIUS_LOWER_EDGE_IN_LARCOL;
  double TOP_SIDE_REFLECTOR_IN_LARCOL_POS_Z;
  double TOP_SIDE_REFLECTOR_IN_LARCOL_POS_Z_LARCOL;
  double SIDE_REFLECTOR_COAT_THICKNESS;
  double BTM_SIDE_REFLECTOR_COAT_OUTER_RADIUS_UPPER_EDGE;
  double BTM_SIDE_REFLECTOR_COAT_INNER_RADIUS_UPPER_EDGE;
  double BTM_SIDE_REFLECTOR_COAT_OUTER_RADIUS_LOWER_EDGE;
  double BTM_SIDE_REFLECTOR_COAT_INNER_RADIUS_LOWER_EDGE;
  double BTM_SIDE_REFLECTOR_COAT_HALF_HEIGHT_GASTEST;
  double BTM_SIDE_REFLECTOR_COAT_POS_Z_GASTEST;
  double BTM_SIDE_REFLECTOR_COAT_POS_Z_GASTEST_LARCOL;
  double BTM_SIDE_REFLECTOR_COAT_START_PHI;
  double BTM_SIDE_REFLECTOR_COAT_DELTA_PHI;
  double BTM_SIDE_REFLECTOR_COAT_CYNLIDRICAL_PART_OUTER_RADIUS;
  double BTM_SIDE_REFLECTOR_COAT_CYNLIDRICAL_PART_INNER_RADIUS;
  double BTM_SIDE_REFLECTOR_COAT_CYNLIDRICAL_PART_HALF_HEIGHT;
  double BTM_SIDE_REFLECTOR_COAT_CYNLIDRICAL_PART_POS_Z;
  double BTM_SIDE_REFLECTOR_COAT_CYNLIDRICAL_PART_POS_Z_LARCOL;
  double TOP_SIDE_REFLECTOR_COAT_OUTER_RADIUS_UPPER_EDGE_IN_GARCOL;
  double TOP_SIDE_REFLECTOR_COAT_INNER_RADIUS_UPPER_EDGE_IN_GARCOL;
  double TOP_SIDE_REFLECTOR_COAT_OUTER_RADIUS_LOWER_EDGE_IN_GARCOL;
  double TOP_SIDE_REFLECTOR_COAT_INNER_RADIUS_LOWER_EDGE_IN_GARCOL;
  double TOP_SIDE_REFLECTOR_COAT_HALF_HEIGHT_IN_GARCOL;
  double TOP_SIDE_REFLECTOR_COAT_IN_GARCOL_POS_Z;
  double TOP_SIDE_REFLECTOR_COAT_IN_GARCOL_POS_Z_GARCOL;
  double TOP_SIDE_REFLECTOR_COAT_OUTER_RADIUS_UPPER_EDGE_IN_LARCOL;
  double TOP_SIDE_REFLECTOR_COAT_INNER_RADIUS_UPPER_EDGE_IN_LARCOL;
  double TOP_SIDE_REFLECTOR_COAT_OUTER_RADIUS_LOWER_EDGE_IN_LARCOL;
  double TOP_SIDE_REFLECTOR_COAT_INNER_RADIUS_LOWER_EDGE_IN_LARCOL;
  double TOP_SIDE_REFLECTOR_COAT_HALF_HEIGHT_IN_LARCOL;
  double TOP_SIDE_REFLECTOR_COAT_IN_LARCOL_POS_Z;
  double TOP_SIDE_REFLECTOR_COAT_IN_LARCOL_POS_Z_LARCOL;
  double TOP_SIDE_REFLECTOR_COAT_START_PHI;
  double TOP_SIDE_REFLECTOR_COAT_DELTA_PHI;
  double FIELD_SHAPER_RING_RINGSEC_CROSSSECTION_OUTER_RADIUS;
  double FIELD_SHAPER_RING_RINGSEC_CROSSSECTION_THICKNESS;
  double FIELD_SHAPER_RING_RINGSEC_CROSSSECTION_INNER_RADIUS;
  double FIELD_SHAPER_RING_RINGSEC_TORUS_RADIUS;
  double FIELD_SHAPER_RING_RINGSEC_TORUS_START_PHI;
  double FIELD_SHAPER_RING_RINGSEC_TORUS_DELTA_PHI;
  double FIELD_SHAPER_LINSEC_INNER_RADIUS;
  double FIELD_SHAPER_LINSEC_OUTER_RADIUS;
  double FIELD_SHAPER_LINSEC_HALF_LENGTH;
  double FIELD_SHAPER_PILLAR_INNER_RADIUS;
  double FIELD_SHAPER_PILLAR_OUTER_RADIUS;
  double FIELD_SHAPER_PILLAR_HALF_HEIGHT;
  double FIELD_SHAPER_PILLAR_POS_Z;
  double FIELD_SHAPER_PILLAR_POS_Z_LAR;
  double DISTANCE_FIELD_SHAPER_PILLAR_CENTER_TO_DETECTOR_CENTER;
  double FIELD_SHAPER_RING_PITCH;
  double DISTANCE_PILLAR_UPPER_EDGE_TO_CENTER_OF_FIRST_FIELD_SHAPER_RING;
  double FIRST_FIELD_SHAPER_RING_POS_Z;
  double PERLITE_COLUMN_INNER_RADIUS;
  double PERLITE_COLUMN_OUTER_RADIUS;
  double PERLITE_COLUMN_HALF_HEIGHT;
  double DISTANCE_PERLITE_COLUMN_UPPER_EDGE_TO_TOP_FLANGE_LOWER_EDGE;
  double PERLITE_COLUMN_UPPER_EDGE_POS_Z;
  double PERLITE_COLUMN_POS_Z;
  double NEUTRON_SHIELD_THICKNESS;
  double NEUTRON_SHIELD_DISTANCE_TO_REFLECTOR;
  double NEUTRON_SHIELD_INNER_RADIUS;
  double NEUTRON_SHIELD_OUTER_RADIUS;
  double NEUTRON_SHIELD_HALF_HEIGHT;
  double NEUTRON_SHIELD_POS_Z;
  double NEUTRON_SHIELD_POS_Z_LAR;
  double PE_SHIELD_INNER_RADIUS;
  double PE_SHIELD_OUTER_RADIUS;
  double PE_SHIELD_HALF_HEIGHT;
  double PE_SHIELD_POS_Z;
  double PE_SHIELD_WIDTH;
  double PE_SHIELD_LEAD_INNER_RADIUS;  // Test for DART
  double PE_SHIELD_LEAD_OUTER_RADIUS;
  double PE_SHIELD_LEAD_HALF_HEIGHT;
  double PE_SHIELD_LEAD_POS_Z;
  double PE_SHIELD_LEAD_WIDTH;
};

inline G4double ArDM::getWLSThickness(G4double convEff) {
  return ((convEff) > 0 && (convEff) < 1) ? (-std::log(1 - (convEff)) * (WLS_MEAN_ABSORPTION_LENGTH)) : (WLS_THICKNESS_100_PERCENT_CONVERSION_EFFICIENCY);
}

inline std::vector<G4TwoVector> ArDM::setPMTVector2D() {
  G4double R1 = 122.4 * mm;
  G4double R2 = 244.8 * mm;
  G4double R3 = 323.84 * mm;

  G4TwoVector r[12];  //[NPMT];

  // first calculate the coordinates of the PMTs
  // assuming that the x-axis is the symmetry axis of the linear sector of the
  // reflector the positive x-direction is from the center of the cross section
  // of the reflector towards the linear sector the 3 first PMTs are located on
  // the innermost circle with radius R1 120 degrees shifted from each other

  // 1st PMT
  r[0] = G4TwoVector(-R1, 0);

  G4TwoVector r0copy = r[0];
  // 2nd and 3rd PMT

  r[1] = r[0];
  r[1].rotate (120 * deg);

  r[2] = r[0];
  r[2].rotate(240 * deg);

  // PMT 3,4,5 are located on the 2nd circle with radius R2
  // relative configuration is the same as for PMT 0,1,2
  // but rotated around the z-axis by 180 degrees
  r[3] = r[0];
  r[3].rotate(180 * deg);
  r[3] *=  R2 / R1;

  r[4] = r[3];
  r[4].rotate(120 * deg);

  r[5] = r[3];
  r[5].rotate(240 * deg);

  // PMT 6,7,8,9,10,11 are on the outermost circle with radius R3
  // 6, 7, 8  are 120   degrees shifted from each other
  // 9,10,11  are 120   degrees shifted from each other
  // 6 is -79.1 degrees away from 3
  // 8 and 9 are 38.21 degrees away from each other
  r[6] = r[3];
  r[6].rotate(-79.1 * deg);
  r[6] *=  R3 / R2;

  r[7] = r[6];
  r[7].rotate(120 * deg);

  r[8] = r[6];
  r[8].rotate(240 * deg);

  r[9] = r[8];
  r[9].rotate(38.21 * deg);
  r[10] = r[9];
  r[10].rotate(120 * deg);
  r[11] = r[9];
  r[11].rotate(240 * deg);

  std::vector<G4TwoVector> rPMT;
  for (int i = 0; i < NPMT; i++) { rPMT.push_back(r[i]); }

  return rPMT;
}
inline std::vector<G4ThreeVector> ArDM::setPMTVector(const char* top_or_bottom) {
  std::vector<G4ThreeVector> rPMT;
  G4double z;

  if (!strcmp(top_or_bottom, "bottom")) {
    z = BTM_PMT_Z;
  } else {
    z = TOP_PMT_Z;
  }

  std::vector<G4TwoVector> rPMT2D = setPMTVector2D();

  for (unsigned int i = 0; i < rPMT2D.size(); i++) { rPMT.push_back(G4ThreeVector(rPMT2D[i].x(), rPMT2D[i].y(), z)); }

  return rPMT;
}

class DSDetectorArDM {

 public:
  DSDetectorArDM(G4VPhysicalVolume*);
  ~DSDetectorArDM();

 private:
  void DefineSurfaces();

  G4VPhysicalVolume* fShield;
  G4VPhysicalVolume* fDartSiPM;
  G4VPhysicalVolume* fTankPhys;
  G4VPhysicalVolume* fMotherVolume;
  G4VPhysicalVolume* fPhysicLArBath;
  G4LogicalVolume* fLogicLArScint;
  G4VPhysicalVolume* fPhysicTank;
  G4VPhysicalVolume* fGArColPhys;

  G4VPhysicalVolume* fTestArgonPhys;

  ArDM ArDMvar;

 protected:
  // add PE shield as a plain cylinder. When it works set the right geometry.
  void addShield(G4VPhysicalVolume* fMotherVolume);
  // esanchez
  // add PE shield as a plain cylinder. When it works set the right geometry.
  void addShieldLead(G4VPhysicalVolume* fMotherVolume);
  // add DART detector plane (fake SiPMs).
  void addDartSiPM(G4VPhysicalVolume* fMotherVolume);

  G4VPhysicalVolume* addTank(G4VPhysicalVolume* fMotherVolume, G4Material* fTankMat, G4ThreeVector fTankPos,
                             int attribute);  //** new ** tank = cylinder + curved bottom part
  G4VPhysicalVolume* addLArColumn(G4VPhysicalVolume* fMotherVolume, G4Material* fMat, G4ThreeVector fPos, int attribute);
  G4VPhysicalVolume* addGArColumn(G4VPhysicalVolume* fMotherVolume, G4Material* fMat, G4ThreeVector fPos, int attribute);

  // Test the Argon batch for DARKSIDE
  G4VPhysicalVolume* addTestArgon(G4VPhysicalVolume* fMotherVolume, G4Material* fMat, G4ThreeVector fPos, int attribute);

  void addWLS();
  G4Material* fWLSMat;
  void addWLSSupport();
  void addPMT(G4String bottom_or_top);
  void addPMTCoat(G4String bottom_or_top);
  void addPMTCathode(G4String bottom_or_top);
  void addPMTSupport(G4String bottom_or_top);
  void addPMTSupportCoat(G4String bottom_or_top);
  void addFieldShaperPillars();
  void addFieldShaperRings();
  // G4VSensitiveDetector* getSD(G4String bottom_or_top);
  void addBtmSideReflector();
  void addTopSideReflector();
  void addBtmSideReflectorCoat();
  void addTopSideReflectorCoat();

  // ArDM_Analysis* fAna;
  G4Material* fShieldMat;
  G4Material* fShieldLeadMat;
  G4Material* fTankMat;
  G4Material* fTeflon;
  G4Material* fWLSSupportMat;
  G4Material* fWLSMatTop;
  G4Material* fPMTMat;
  G4Material* fPMTCathodeMat;
  G4Material* fPMTSupportMat;

  G4LogicalVolume* constructPMT(G4Material* fPMTMat, const char* bottom_or_top, int attribute, const char* name);
  G4LogicalVolume* constructPMTAux(G4Material* fPMTMat, const char* bottom_or_top, int attribute, G4double pmt_inner_radius, G4double pmt_outer_radius, const char* name, G4double deltaTheta);
  std::vector<G4VPhysicalVolume*> fTopPMTCoatArrayPhys;
  std::vector<G4VPhysicalVolume*> fBtmPMTCoatArrayPhys;
  std::vector<G4VPhysicalVolume*> fBtmPMTCathodeArrayPhys;
  std::vector<G4VPhysicalVolume*> fTopPMTCathodeArrayPhys;
  // physical volume of the polyethylene pillars, there are 7 of them around +
  // outside the main reflector
  std::vector<G4VPhysicalVolume*> fFieldShaperPillarPhys;

  G4VPhysicalVolume* fWLSGArPhys;
  G4VPhysicalVolume* fWLSLArPhys;
  G4VPhysicalVolume* fLArColPhys;

  G4VPhysicalVolume* fWLSSupportGArPhys;
  G4VPhysicalVolume* fWLSSupportLArPhys;

  G4VPhysicalVolume* fTopPMTSupportPhys;
  G4VPhysicalVolume* fBtmPMTSupportPhys;

  G4VPhysicalVolume* fTopPMTSupportCoatPhys;  // coating of support tructure for PMT
  G4VPhysicalVolume* fBtmPMTSupportCoatPhys;  // coating of support tructure for PMT

  G4LogicalVolume* fBtmSideReflectorLog;

  G4VPhysicalVolume* fTopSideReflectorGArPhys;
  G4VPhysicalVolume* fTopSideReflectorLArPhys;

  G4VPhysicalVolume* fBtmSideReflectorCoatPhys;
  G4VPhysicalVolume* fTopSideReflectorCoatGArPhys;
  G4VPhysicalVolume* fTopSideReflectorCoatLArPhys;

  std::vector<G4VPhysicalVolume*> placePMT(G4LogicalVolume* fPMTLog, G4VPhysicalVolume* fMotherPhys, std::vector<G4ThreeVector> rPMT, const char* bottom_or_top);
  std::vector<G4VPhysicalVolume*> placePMT(std::vector<G4LogicalVolume*> fPMTLog, G4VPhysicalVolume* fMotherLog, std::vector<G4ThreeVector> rPMT, const char* name);
  std::vector<G4VPhysicalVolume*> fBtmPMTArrayPhys;
  std::vector<G4VPhysicalVolume*> fTopPMTArrayPhys;

  // physical volume of the field shaper rings, there are 27 of them
  std::vector<G4VPhysicalVolume*> fFieldShaperRingsPhys;
  // DEfine surface stuff ........

  // Rods
  G4Tubs* fSolidSupportRod;
  G4LogicalVolume* fLogicSupportRod;
  G4VPhysicalVolume* fPhysicSupportRod0;
  G4VPhysicalVolume* fPhysicSupportRod1;
  G4VPhysicalVolume* fPhysicSupportRod2;
  G4VPhysicalVolume* fPhysicSupportRod3;

  // Vessel
  G4Tubs* fSolidInnerVesselWall;
  G4LogicalVolume* fLogicInnerVesselWall;
  G4VPhysicalVolume* fPhysicInnerVesselWall;

  G4Tubs* fSolidInnerVesselWindow;
  G4LogicalVolume* fLogicInnerVesselWindow;
  G4VPhysicalVolume* fPhysicInnerVesselWindowTop;
  G4VPhysicalVolume* fPhysicInnerVesselWindowBottom;

  // Gas Pocket and Active LAr
  G4Tubs* fSolidGasPocket;
  G4LogicalVolume* fLogicGasPocket;
  G4VPhysicalVolume* fPhysicGasPocket;

  /// Edgar ArDM
  G4VPhysicalVolume* fDeepArgonTestArgonPhys;
  G4VPhysicalVolume* fWLSTestArgonPhys;
  G4VPhysicalVolume* fAcrylicTestArgonPhys;
  G4VPhysicalVolume* fReflectorTestArgonPhys;
  G4VPhysicalVolume* fSiPMup;
  G4VPhysicalVolume* fSiPMdw;
  G4VPhysicalVolume* fBatchArgonPhys;
  G4VPhysicalVolume* fCopperTestArgonPhys;
  G4VPhysicalVolume* fDiskupTPB;
  G4VPhysicalVolume* fDiskdwTPB;
  G4VPhysicalVolume* fDiskup;
  G4VPhysicalVolume* fDiskdw;

  G4Tubs* fSolidInnerLAr;
  G4LogicalVolume* fLogicInnerLAr;
  G4VPhysicalVolume* fPhysicInnerLAr;

  // Grid
  G4Tubs* fSolidGrid;
  G4LogicalVolume* fLogicGrid;
  G4VPhysicalVolume* fPhysicGrid;

  // Active Volume Boundaries
  G4Tubs* fSolidThreeMLowerFoil;
  G4Tubs* fSolidThreeMMiddleFoil;
  G4Tubs* fSolidThreeMTopFoil;
  G4LogicalVolume* fLogicThreeMLowerFoil;
  G4LogicalVolume* fLogicThreeMMiddleFoil;
  G4LogicalVolume* fLogicThreeMTopFoil;
  G4VPhysicalVolume* fPhysicThreeMLowerFoil;
  G4VPhysicalVolume* fPhysicThreeMMiddleFoil;
  G4VPhysicalVolume* fPhysicThreeMTopFoil;

  // 3M foil support rings
  G4Tubs* fSolidSupportRing;
  G4Tubs* fSolidTopLiqSupportRing;
  G4Tubs* fSolidTopGasSupportRing;
  G4LogicalVolume* fLogicSupportRing;
  G4LogicalVolume* fLogicTopLiqSupportRing;
  G4LogicalVolume* fLogicTopGasSupportRing;
  G4VPhysicalVolume* fPhysicSupportRing0;
  G4VPhysicalVolume* fPhysicSupportRing1;
  G4VPhysicalVolume* fPhysicSupportRing2;
  G4VPhysicalVolume* fPhysicTopLiqSupportRing;
  G4VPhysicalVolume* fPhysicTopGasSupportRing;

  // TPB
  G4Tubs* fSolidTPBLowerLateral;
  G4Tubs* fSolidTPBMiddleLateral;
  G4Tubs* fSolidTPBTopLateral;
  G4Tubs* fSolidTPBBases;
  G4LogicalVolume* fLogicTPBLowerLateral;
  G4LogicalVolume* fLogicTPBMiddleLateral;
  G4LogicalVolume* fLogicTPBTopLateral;
  G4LogicalVolume* fLogicTPBBottomBase;
  G4LogicalVolume* fLogicTPBTopBase;

  G4VPhysicalVolume* fPhysicTPBLowerLateral;
  G4VPhysicalVolume* fPhysicTPBMiddleLateral;
  G4VPhysicalVolume* fPhysicTPBTopLateral;
  G4VPhysicalVolume* fPhysicTPBBottomBase;
  G4VPhysicalVolume* fPhysicTPBTopBase;

  // Field Rings and Kapton cover
  G4Tubs* fSolidFieldRing;
  G4LogicalVolume* fLogicFieldRing;
  G4VPhysicalVolume* fPhysicFieldRings;

  G4Tubs* fSolidKaptonBand;
  G4LogicalVolume* fLogicKaptonBand;
  G4VPhysicalVolume* fPhysicKaptonBand;

  // Compression plate
  G4Tubs* fSolidCompressionPlateTmp;
  G4Tubs* fSolidPMTMold;
  G4SubtractionSolid* fSolidCompressionPlate[7];
  G4LogicalVolume* fLogicCompressionPlate;
  G4VPhysicalVolume* fPhysicCompressionPlateTop;
  G4VPhysicalVolume* fPhysicCompressionPlateBottom;

  // Rods between the Compression plates
  G4Tubs* fSolidSteelRod;
  G4LogicalVolume* fLogicSteelRod;
  G4VPhysicalVolume* fPhysicThisSteelRod[12];

  // PMTs
  G4Tubs* fSolidPMTWindow;
  G4Tubs* fSolidPMTBody;
  G4Tubs* fSolidPMTVacuum;
  G4LogicalVolume* fLogicPMTBody;
  G4LogicalVolume* fLogicPMTVacuum;
  G4LogicalVolume* fLogicPMTWindow[14];
  G4VPhysicalVolume* fPhysicPMTWindow[14];
  G4VPhysicalVolume* fPhysicPMTBody[14];
  G4VPhysicalVolume* fPhysicPMTVacuum;

  // Bubbler Tubes
  G4Tubs* fSolidBubblerTube;
  G4LogicalVolume* fLogicBubblerTube;
  G4VPhysicalVolume* fPhysicBubblerTube1;
  G4VPhysicalVolume* fPhysicBubblerTube2;
};

#endif
