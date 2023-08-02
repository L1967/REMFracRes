//----------------------------fracResFileKeyWord.h---------------------------
//
//
//Description:	fracResFileKeyWordclassusestoredataparameterFracResfile
//					keyword
//
//Documents:		fracResFileKeyWord.h
//
//Creationdate:20/11/2022
//
//Author:			PascalSiegel-RabahNamar
//
//Copyright(C)2022byin2earth
//
//----------------------------fracResFileKeyWord.h---------------------------
#pragma once

#ifndef fracResFileKeyWord_h
#define fracResFileKeyWord_h

//localincludes

//standardincludes
#include<fstream>
#include<iostream>
#include<iomanip>
#include<memory>
#include<vector>

#include<utility>



namespace remfracres{

using namespace std;
/**
*fracResFileKeyWordcontainsthekeywordparametersforFracResdataParameter
*
*@authorPascalSiegel-RabahNamar-OlivierJaquet
*@date2022
*/

class FracResFileKeyWord
{


//variables

public:


	  const double ndvalue =-9999.0;

	 const pair< string ,string > distribution_min={"a","Min"};

	 const pair< string ,string > distribution_max ={"b","Max"};
	 const pair< string ,string > distribution_min_triangular ={"min","Min"};
	 const pair< string ,string > distribution_max_triangular ={"max","Max"};
	 const pair< string ,string > distribution_mode ={"mode","Mode"};
	 const pair< string ,string > distribution_alpha ={"alpha","Alpha"};
	 const pair< string ,string > distribution_theta ={"theta","Theta"};
	 const pair< string ,string > distribution_mean ={"m","Mean"};
	 const pair< string ,string > distribution_mean_log_normal ={"mean","Mean"};
	 const pair< string ,string > distribution_std_dev ={"sigma","StdDev"};
	 const pair< string ,string > distribution_std_dev_log_normal ={"s","StdDev"};
	 const string  distribution_key_val_separator = ":";

	 const string new_line = "\n";
	 const string separator = ";";
	 const string value_separator = ",";
	 const string vale_name = "TEST";
	 const string fracres_model_name = "MODEL_NAME";
	 const string struct_model_entity_id = "STRUCT_MODEL_ENTITY_ID";

	 const string struct_model_entity_name = "STRUCT_MODEL_ENTITY_NAME";
	 const string struct_model_path = "STRUCTURAL_MODEL_PATH";
	 const string data_support_entity_id = "DATA_SUPPORT_ENTITY_ID";
	 const string data_support_name = "DATA_SUPPORT_NAME";
	 const string manage_bed_interface = "MANAGE_BED_INTERFACE";
	 const string manage_bed_parallel = "MANAGE_BED_PARALLEL";

	//Units
	 const string unit_begin = "UNIT_BEGIN";
	 const string unit_age = "UNIT_AGE";
	 const string unit_duration = "UNIT_DURATION";
	 const string unit_distance = "UNIT_DISTANCE";
	 const string unit_altitude = "UNIT_ALTITUDE";
	 const string unit_production_duration = "UNIT_PRODUCTION_DURATION";
	 const string unit_erosion_duration = "UNIT_EROSION_DURATION";
	 const string unit_subsidence_duration = "UNIT_SUBSIDENCE_DURATION";
	 const string unit_river_sediment_volume = "UNIT_RIVER_SEDIMENT_VOLUME";
	 const string unit_river_sediment_duration = "UNIT_RIVER_SEDIMENT_DURATION";
	 const string unit_temperature = "UNIT_TEMPERATURE";
	 const string unit_slope = "UNIT_SLOPE";
	 const string unit_end = "UNIT_END";

	//stage
	 const string stage_begin = "STAGE_BEGIN";
	 const string stage_name = "STAGE_NAME";
	 const string stage_index = "STAGE_INDEX";
	 const string stage_begin_stage = "STAGE_BEDDING_STAGE";
	 const string stage_end = "STAGE_END";

	//Fracture
	 const string fracture_begin = "FRACTURE_BEGIN";
	 const string fracture_id = "FRACTURE_ID";//Onlyforsave/load,toknowwitchfracturetolink
	 const string fracture_name = "FRACTURE_NAME";
	 const string fracture_category = "FRACTURE_CATEGORY";
	 const string fracture_type = "FRACTURE_TYPE";
	 const string fracture_end = "FRACTURE_END";

	//Hypothesis
	 const string hyp_begin = "HYP_BEGIN";
	 const string hyp_name = "HYP_NAME";
	 const string hyp_index = "HYP_INDEX";
	 const string hyp_bed_is_region = "HYP_BED_IS_REGION";
	 const string hyp_bed_region = "HYP_BED_REGION";
	 const string hyp_bed_region_name = "HYP_BED_REGION_NAME";
	 const string hyp_bed_prop_discrete_prop = "HYP_BED_PROP_DISCRETE_PROP";
	 const string hyp_bed_prop_discrete_name = "HYP_BED_PROP_DISCRETE_PROP_NAME";
	 const string hyp_end = "HYP_END";

	//beddingaperture
	 const string bed_aperture_begin = "BED_APERTURE_BEGIN";
	 const string bed_aperture_name = "BED_APERTURE_NAME";
	 const string bed_aperture_parameter = "BED_APERTURE_PARAMETER";
	 const string bed_aperture_unit = "BED_APERTURE_UNIT";
	 const string bed_aperture_data_type = "BED_APERTURE_DATA_TYPE";
	 const string bed_aperture_value = "BED_APERTURE_VALUE";
	 const string bed_aperture_property = "BED_APERTURE_PROPERTY";
	 const string bed_aperture_property_name = "BED_APERTURE_PROPERTY_NAME";
	 const string bed_aperture_distribution = "BED_APERTURE_DISTRIBUTION";
	 const string bed_aperture_distribution_str = "BED_APERTURE_DISTRIBUTION_STR";
	 const string bed_aperture_discrete_prop_values = "BED_APERTURE_DISCRETE_PROP_VALUES";
	 const string bed_aperture_end = "BED_APERTURE_END";

	//Fractureproperties
	 const string frac_prop_begin = "FRAC_PROP_BEGIN";
	 const string frac_prop_selected = "FRAC_PROP_SELECTED";
	 const string frac_prop_fracture_id = "FRAC_PROP_FRACTURE_ID";//Onlyuseforsave/load,toretrievethefracture
	 const string frac_prop_index = "FRAC_PROP_INDEX";
	 const string frac_prop_clusters_type = "FRAC_PROP_CLUSTERS_TYPE";
	 const string frac_prop_use_previous = "FRAC_PROP_USE_PREVIOUS";
	 const string frac_prop_height_correl_type = "FRAC_PROP_HEIGHT_CORREL_TYPE";
	 const string frac_prop_height_correl_value = "FRAC_PROP_HEIGHT_CORREL_VALUE";
	 const string frac_prop_height_correl_variable = "FRAC_PROP_HEIGHT_CORREL_VARIABLE";
	 const string frac_prop_height_correl_property = "FRAC_PROP_HEIGHT_CORREL_PROPERTY";
	 const string frac_prop_height_correl_property_name = "FRAC_PROP_HEIGHT_CORREL_PROPERTY_NAME";

	 const string frac_prop_aperture_correl_type = "FRAC_PROP_APERTURE_CORREL_TYPE";
	 const string frac_prop_aperture_correl_value = "FRAC_PROP_APERTURE_CORREL_VALUE";
	 const string frac_prop_aperture_correl_variable = "FRAC_PROP_APERTURE_CORREL_VARIABLE";
	 const string frac_prop_aperture_correl_property = "FRAC_PROP_APERTURE_CORREL_PROPERTY";
	 const string frac_prop_aperture_correl_property_name = "FRAC_PROP_APERTURE_CORREL_PROPERTY_NAME";

	 const string frac_prop_fault_core_type = "FRAC_PROP_FAULT_CORE_TYPE";
	 const string frac_prop_fault_core_value = "FRAC_PROP_FAULT_CORE_VALUE";
	 const string frac_prop_fault_core_property = "FRAC_PROP_FAULT_CORE_PROPERTY";
	 const string frac_prop_fault_core_property_name = "FRAC_PROP_FAULT_CORE_PROPERTY_NAME";

	 const string frac_prop_fault_core_distribution = "FRAC_PROP_FAULT_CORE_DISTRIBUTION";
	 const string frac_prop_fault_core_distribution_str = "FRAC_PROP_FAULT_CORE_DISTRIBUTION_STR";

	 const string frac_prop_corridor_fractures = "FRAC_PROP_CORRIDOR_FRACTURES";
	 const string frac_prop_rate_reactivate = "FRAC_PROP_RATE_REACTIVATE";
	 const string frac_prop_end = "FRAC_PROP_END";

	//BeDParallelslip
	 const string bed_parrallel_begin = "BED_PARALLEL_BEGIN";
	 const string bed_parrallel_is_region = "BED_PARALLEL_IS_REGION";
	 const string bed_parrallel_region = "BED_PARALLEL_REGION";
	 const string bed_parrallel_region_name = "BED_PARALLEL_REGION_NAME";
	 const string bed_parrallel_prop_discrete_prop = "BED_PARALLEL_PROP_DISCRETE_PROP";
	 const string bed_parrallel_prop_discrete_prop_name = "BED_PARALLEL_PROP_DISCRETE_PROP_NAME";
	 const string bed_parrallel_slope_parameter = "BED_PARALLEL_SLOPE_PARAMETER";
	 const string bed_parrallel_slope_unit = "BED_PARALLEL_SLOPE_UNIT";
	 const string bed_parrallel_slope_data_type = "BED_PARALLEL_SLOPE_DATA_TYPE";
	 const string bed_parrallel_slope_value = "BED_PARALLEL_SLOPE_VALUE";
	 const string bed_parrallel_slope_prop = "BED_PARALLEL_SLOPE_PROP";
	 const string bed_parrallel_slope_prop_name = "BED_PARALLEL_SLOPE_PROP_NAME";
	 const string bed_parrallel_slope_distribution = "BED_PARALLEL_SLOPE_DISTRIBUTION";
	 const string bed_parrallel_slope_distribution_str = "BED_PARALLEL_SLOPE_DISTRIBUTION_STR";
	 const string bed_parrallel_slope_ratio = "BED_APERTURE_SLOPE_RATIO";

	 const string bed_parrallel_end = "BED_PARALLEL_END";

	//Intersection
	 const string intersection_begin = "INTERSECTION_BEGIN";
	 const string intersection_fracture_type = "INTERSECTION_FRACTURE_TYPE";
	 const string intersection_type = "INTERSECTION_TYPE";
	 const string intersection_type_modifiable = "INTERSECTION_TYPE_MODIFIABLE";
	 const string intersection_end = "INTERSECTION_END";

	//Distribution
	 const string distribution_begin = "DISTRIBUTION_BEGIN";
	 const string distribution_parameter_type = "DISTRIBUTION_PARAMETER_TYPE";
	 const string distribution_array_type = "DISTRIBUTION_ARRAY_TYPE";
	 const string distribution_data_type = "DISTRIBUTION_DATA_TYPE";
	 const string distribution_unit = "DISTRIBUTION_UNIT";
	 const string distribution_value = "DISTRIBUTION_VALUE";
	 const string distribution_property = "DISTRIBUTION_PROPERTY";
	 const string distribution_property_name = "DISTRIBUTION_PROPERTY_NAME";
	 const string distribution_distribution = "DISTRIBUTION_DISTRIBUTION";
	 const string distribution_distribution_str = "DISTRIBUTION_DISTRIBUTION_STR";
	 const string distribution_end = "DISTRIBUTION_END";

	//Geometry
	 const string geometry_begin = "GEOMETRY_BEGIN";
	 const string geometry_parameter_type = "GEOMETRY_PARAMETER_TYPE";
	 const string geometry_array_type = "GEOMETRY_ARRAY_TYPE";

	 const string geometry_data_type = "GEOMETRY_DATA_TYPE";
	 const string geometry_unit = "GEOMETRY_UNIT";
	 const string geometry_value = "GEOMETRY_VALUE";
	 const string geometry_property = "GEOMETRY_PROPERTY";
	 const string geometry_property_name = "GEOMETRY_PROPERTY_NAME";
	 const string geometry_distribution = "GEOMETRY_DISTRIBUTION";
	 const string geometry_distribution_str = "GEOMETRY_DISTRIBUTION_STR";
	 const string geometry_normal_distribution_min = "GEOMETRY_NORMAL_DISTRIBUTION_MIN";
	 const string geometry_normal_distribution_max = "GEOMETRY_NORMAL_DISTRIBUTION_MAX";
	 const string geometry_end = "GEOMETRY_END";
	//Representation&Upscaling
	 const string rep_hypothesis_begin = "REP_HYPOTHESIS_BEGIN";
	 const string rep_hypothesis_name = "REP_HYPOTHESIS_NAME";
	 const string rep_hypothesis_index = "REP_HYPOTHESIS_INDEX";
	 const string rep_hypothesis_grid_id = "REP_HYPOTHESIS_GRID_ID";
	 const string rep_hypothesis_grid_name = "REP_HYPOTHESIS_GRID_NAME";
	 const string rep_hypothesis_upscaling_mode = "REP_HYPOTHESIS_UPSCALING_MODE";
	 const string rep_hypothesis_upscaling_method = "REP_HYPOTHESIS_UPSCALING_METHOD";
	 const string rep_hypothesis_boundary = "REP_HYPOTHESIS_BOUNDARY";
	 const string rep_hypothesis_resolution = "REP_HYPOTHESIS_RESOLUTION";
	 const string rep_hypothesis_end = "REP_HYPOTHESIS_END";

	//Reptable
	 const string representation_begin = "REPRESENTATION_BEGIN";
	 const string representation_fracture_type = "REPRESENTATION_FRACTURE_TYPE";
	 const string representation_rep = "REPRESENTATION_REP";
	 const string representation_medium = "REPRESENTATION_MEDIUM";
	 const string representation_cut_off_param = "REPRESENTATION_CUT_OFF_PARAM";
	 const string representation_cut_off_val = "REPRESENTATION_CUT_OFF_VAL";
	 const string representation_end = "REPRESENTATION_END";

	//Scenario
	 const string scenario_begin = "SCENARIO_BEGIN";
	 const string scenario_is_modeling = "SCENARIO_IS_MODELING";
	 const string scenario_is_qc = "SCENARIO_IS_QC";
	 const string scenario_nb_real = "SCENARIO_NB_REAL";
	 const string scenario_rep_hyp_index = "SCENARIO_REP_HYP_INDEX";
	 const string scenario_nb_estimate_fractures = "SCENARIO_NB_ESTIMATE_FRACTURES";

	 const string scenario_name = "SCENARIO_NAME";
	 const string scenario_index = "SCENARIO_INDEX";
	 const string scenario_stage_hyp_link = "-";
	 const string scenario_stage_hyp_separator = ",";
	 const string scenario_stages_hyp_map = "SCENARIO_STAGES_HYPS_MAP";
	 const string scenario_end = "SCENARIO_END";
		//RUN param
	 const string run_param_begin = "RUN_PARAM_BEGIN";
	 const string run_param_seed = "RUN_PARAM_SEED";
	 const string run_param_save_upscaled = "RUN_PARAM_SAVE_UPSCALED";
	 const string run_param_save_dfn = "RUN_PARAM_SAVE_DFN";
	 const string run_param_target_mesh_res = "RUN_PARAM_TARGET_MESH_RES";
	 const string run_param_save_density = "RUN_PARAM_SAVE_DENSITY";
	 const string run_param_denisty_enum = "RUN_PARAM_DENSITY_ENUM";
	 const string run_param_end = "RUN_PARAM_END";

		//MATRIX prop
	 const string matrix_prop_begin = "MATRIX_PROP_BEGIN";
	 const string matrix_prop_init_poro_name = "MATRIX_PROP_INIT_PORO_NAME";
	 const string matrix_prop_init_poro_selected = "MATRIX_PROP_INIT_PORO_SELECTED";
	 const string matrix_prop_init_xperm_name = "MATRIX_PROP_INIT_X_PERM_NAME";
	 const string matrix_prop_init_xperm_selected = "MATRIX_PROP_INIT_X_PERM_SELECTED";
	 const string matrix_prop_init_yperm_name = "MATRIX_PROP_INIT_Y_PERM_NAME";
	 const string matrix_prop_init_yperm_selected = "MATRIX_PROP_INIT_Y_PERM_SELECTED";
	 const string matrix_prop_init_zperm_name = "MATRIX_PROP_INIT_Z_PERM_NAME";
	 const string matrix_prop_init_zperm_selected = "MATRIX_PROP_INIT_Z_PERM_SELECTED";
	 const string matrix_prop_end = "MATRIX_PROP_END";
	 // OUTPUTS
	 const string outputs_begin = "OUTPUTS_BEGIN";
	 const string outputs_path = "OUTPUTS_PATH";
	 const string outputs_density= "OUTPUTS_DENSITY";
	 const string outputs_density_selected = "OUTPUTS_DENSITY_SELECTED";
	 const string outputs_first_second_exchange = "OUTPUTS_FIRS_SECOND_EXCHANGE";
	 const string outputs_first_second_exchange_selected = "OUTPUTS_FIRS_SECOND_EXCHANGE_SELECTED";
	 const string outputs_second_third_exchange = "OUTPUTS_SECOND_THIRD_EXCHANGE";
	 const string outputs_second_third_exchange_selected = "OUTPUTS_SECOND_THIRD_EXCHANGE_SELECTED";
	 const string outputs_first_third_exchange = "OUTPUTS_FIRST_THIRD_EXCHANGE";
	 const string outputs_first_third_exchange_selected = "OUTPUTS_FIRST_THIRD_EXCHANGE_SELECTED";
	 const string outputs_fractures = "OUTPUTS_FRACTURES";
	 const string outputs_fractures_selected = "OUTPUTS_FRACTURES_SELECTED";
	 const string outputs_aperture = "OUTPUTS_APERTURE";
	 const string outputs_aperture_selected = "OUTPUTS_APERTURE_SELECTED";
	 const string outputs_length = "OUTPUTS_LENGTH";
	 const string outputs_length_selected = "OUTPUTS_LENGTH_SELECTED";
	 const string outputs_height = "OUTPUTS_HEIGHT";
	 const string outputs_height_selected = "OUTPUTS_HEIGHT_SELECTED";
	 const string outputs_dip_azimut = "OUTPUTS_DIP_AZIMUTH";
	 const string outputs_dip_azimut_selected = "OUTPUTS_DIP_AZIMUTH_SELECTED";
	 const string outputs_dip_angle = "OUTPUTS_DIP_ANGLE";
	 const string outputs_dip_angle_selected = "OUTPUTS_DIP_ANGLE_SELECTED";
	 const string outputs_mult_t = "OUTPUTS_MULT_T";
	 const string outputs_mult_t_selected = "OUTPUTS_MULT_T_SELECTED";
	 const string outputs_clusters = "OUTPUTS_CLUSTERS";
	 const string outputs_clusters_selected = "OUTPUTS_CLUSTERS_SELECTED";
	 const string outputs_clusters_length = "OUTPUTS_CLUSTERS_LENGTH";
	 const string outputs_clusters_length_selected = "OUTPUTS_CLUSTERS_LENGTH_SELECTED";
	 const string outputs_clusters_width = "OUTPUTS_CLUSTERS_WIDTH";
	 const string outputs_clusters_width_selected = "OUTPUTS_CLUSTERS_WIDTH_SELECTED";
	 const string outputs_clusters_height = "OUTPUTS_CLUSTERS_HEIGHT";
	 const string outputs_clusters_height_selected = "OUTPUTS_CLUSTERS_HEIGHT_SELECTED";
	 const string outputs_end = "OUTPUTS_END";
public :
	 /**
	 	 *      Default constructor
	 	 */

	 FracResFileKeyWord(){

	 };


	 	/**
	 	 *		Destructor
	 	 */

	 ~FracResFileKeyWord(){};


};
//inlinemethods


}
#endif//fracResFileKeyWord_h
