/*[
 * Godiag Software.
 *
 * Copyright, TotalFinaElf and Colenco Power Engineering, 
 * 2002-2003, All rights reserved.
 *
 * The software and information contained herein are proprietary to, and
 * comprise valuable trade secrets of, TotalFinaElf and Colenco Power Engineering
 * which intend to preserve as trade secrets such software and information.
 * This software is furnished pursuant to a written license agreement and
 * may be used, copied, transmitted, and stored only in accordance with
 * the terms of such license and with the inclusion of the above copyright
 * notice.  This software and information or any other copies thereof may
 * not be provided or otherwise made available to any other person.
 *
]*/


//----------------------------FractureSet.cpp  ---------------------------
// 
//  FractureSet contains the parameters for the simulation 
//  of a Boolean model and the generated fractures
//  
//
//  Author:  Pascal Siegel - Rabah Namar - Olivier Jaquet 
//
//	Date: 2002-2003
//
//----------------------------  FractureSet.cpp  ---------------------------


// local includes


#include <remfracres/fractureset.h>

#include <memory>
#include <string>
#include <array>

// lib includes 

#include <geode/mesh/core/tetrahedral_solid.h>
#include <geode/geometry/basic_objects/tetrahedron.h>
#include <geode/geometry/mensuration.h>
#include <geode/mesh/helpers/convert_surface_mesh.h>

#include <geode/mesh/core/geode/geode_triangulated_surface.h>
#include <geode/mesh/core/geode/geode_polygonal_surface.h>
#include <geode/mesh/core/polygonal_surface.h>
#include <geode/mesh/core/triangulated_surface.h>
#include <geode/mesh/io/triangulated_surface_input.h>
#include <geode/mesh/io/triangulated_surface_output.h>
#include <geode/mesh/io/polygonal_surface_input.h>
#include <geode/mesh/io/polygonal_surface_output.h>
#include <geode/mesh/builder/geode/geode_triangulated_surface_builder.h>
#include <geode/mesh/builder/geode/geode_polygonal_surface_builder.h>
#include <geode/mesh/builder/triangulated_surface_builder.h>
#include <geode/mesh/builder/polygonal_surface_builder.h>

//#include <boost/random.hpp>
#include <boost/random.hpp>
#include <boost/math/distributions/uniform.hpp>
#include <boost/math/distributions.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/triangle_distribution.hpp>
#include <boost/random/lognormal_distribution.hpp>
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/math/distributions/lognormal.hpp>
#include <boost/math/distributions/complement.hpp>
#include <boost/math/distributions/pareto.hpp>
// system includes                       

namespace remfracres {

using namespace geode;
using namespace std;
using namespace fractures_intersect;



//-----------------------------------------------
// default constructor of the class FractureSet
//-----------------------------------------------

FractureSet::FractureSet() 
{
	name_= "set0";

	length_1_ = "";

	length_2_  = "";
	
	dip_azimut_  = "";

	dip_angle_  = "";
	
	aperture_name_  = "";

	density_ = 1e-6;

	complementary_fracture_set_active_ = false;

	density_function_=NULL;

	coef_a_=0.0; 
		
	seed_=1337;

	fractures_number_ = 0;
	fractures_number_foot_wall_ = 0;

	fractures_number_hanging_wall_=0;

	fractures_activable_ = 0;

	fractures_inserted_ = 0;

	evaluate_fractures_number_ = 0;
	evaluate_fractures_number_foot_wall_ = 0;
	evaluate_fractures_number_hanging_wall_ = 0;

	fractures_activable_hanging_wall_ = 0;
	fractures_activable_foot_wall_ = 0;

	//fractures_ = NULL;

	density_property_name_ = "none";
	constraint_property_name_ = "none"; 

	shadow_zone_distribution_name_ = "";
	shadow_zone_property_name_ = "none";

	is_edit_ = false;

    is_imported_ = false;
	is_created_ = false;
	distribution_length1_type_index_ = 0;
	distribution_length2_type_index_ = 0;
	distribution_azimut_type_index_ = 0;
	distribution_angle_type_index_ = 0;
	distribution_aperture_type_index_ = 0;
	distribution_shadow_zone_type_index_ =0;
	density_type_hanging_wall_index_ = 0;
	density_type_foot_wall_index_ = 0;

	distr_length1_val_ = 100.0;
	distr_length2_val_ = 50;
	distr_azimut_val_ = 180.0;;
	distr_angle_val_= 45.0;
	distr_aperture_val_= 1.0;
	distr_shadow_zone_val_=0.5;
	distr_hanging_wall_val_ = 1.0;
	distr_foot_wall_val_=1.0;
	density_type_index_ = 0;
	density_law_type_ =0;

	to_generate_ = false;
	density_property_name_ = "density";
	length1_type_ = distribution_type::uniform;
	length2_type_= distribution_type::uniform;
	azimuth_type_= distribution_type::uniform;
	dip_angle_type_= distribution_type::uniform;
	aperture_type_= distribution_type::uniform;
	density_type_ = distribution_type::uniform;
	shadow_zone_type_ = distribution_type::uniform;
	density_hanging_wall_type_ = distribution_type::uniform;
	density_foot_wall_type_ = distribution_type::uniform;


	minmaxmode_density_values_.resize(3);
	minmaxmode_density_values_[0]=0;
	minmaxmode_density_values_[1]=1;
	minmaxmode_density_values_[2]=0;

	minmaxmode_shadow_zone_values_.resize(3);
	minmaxmode_shadow_zone_values_[0]=0;
	minmaxmode_shadow_zone_values_[1]=1;
	minmaxmode_shadow_zone_values_[2]=0;

	minmaxmode_values_.resize(15);

	for(int i= 0 ; i < 5; i++){
	 minmaxmode_values_[3*i]=0;
	 minmaxmode_values_[3*i + 1]=1;
	 minmaxmode_values_[3*i + 2]=0;
	}
	model_fractures_ = fractures_intersect::model();
	fracture_set_index_=1;
	default_stress_zone_max_ = 0.0;
	directory_path_ = remfracres::data_path;
	is_shadow_zone_active_ =true;
	success_ = false;

	cylinder_step_ = 12;
	tab_type_intersection_.resize(8, fractures_intersect::FracResIntersectionTypeEnum::RANDOM);
	type_category_ = fractures_intersect::FracResIntersectionCategoryEnum::FAULTS;
	type_intersection_= fractures_intersect::FracResIntersectionTypeEnum::RANDOM;
	array_index_ =0;
	type_array_ = array_type::unkown;
	echelon_space_ = 0.1;

	max_distance_connectivity_ = 0.2;
	nb_fractures_connection_relay_ =0.0;

	width_stress_zone_law_type_=0;

	width_stress_zone_distribution_name_="None";

	width_stress_zone_min_max_mode_.resize(3, 0.0);
	width_stress_zone_min_max_mode_[1]=1.0;

	width_stress_zone_property_name_="None";

	density_distribution_name_ ="";
	density_min_max_.resize(3, 0);
	density_min_max_[1]=1.0;


	minmaxmode_density_hanging_wall_values_.resize(3);
	minmaxmode_density_hanging_wall_values_[0]=0;
	minmaxmode_density_hanging_wall_values_[1]=1;
	minmaxmode_density_hanging_wall_values_[2]=0;

	minmaxmode_density_foot_wall_values_.resize(3);
	minmaxmode_density_foot_wall_values_[0]=0;
	minmaxmode_density_foot_wall_values_[1]=1;
	minmaxmode_density_foot_wall_values_[2]=0;


	distribution_echelon_type_index_ =0;
	echelon_space_type_= distribution_type::uniform;
	minmaxmode_echelon_space_values_.resize(3);
	minmaxmode_echelon_space_values_[0]=0;
	minmaxmode_echelon_space_values_[1]=1;
	minmaxmode_echelon_space_values_[2]=0;
	echelon_space_distribution_name_ = "";

    is_active_ = false;
    model_geo_cont_index_ =0;
    output_prefixe_name_ = "";
	total_simulation_time_ = 0.0;
	total_writting_output_time_ = 0.0;
	length_max_ = 0;

}


//---------------------------------------------
//  copy constructor of the class FractureSet
//---------------------------------------------

FractureSet::FractureSet(const FractureSet& from)
{		
	name_= from.name_;

	length_1_ = from.length_1_;

	length_2_  = from.length_2_;
	
	dip_azimut_  = from.dip_azimut_;

	dip_angle_  = from.dip_angle_;
	
	aperture_name_  = from.aperture_name_;

	density_ = from.density_;

	complementary_fracture_set_active_ = from.complementary_fracture_set_active_;

	density_function_= from.density_function_;

	coef_a_= from.coef_a_; 
		
	seed_= from.seed_;

	fractures_number_ = from.fractures_number_;

	evaluate_fractures_number_ = from.evaluate_fractures_number_;

	fractures_activable_ = from.fractures_activable_;

	fractures_inserted_ = from.fractures_inserted_;

	//fractures_ = from.fractures_;

	generation_box_ = from.generation_box_;

	oriented_generation_box_ = from.oriented_generation_box_;

	oriented_generation_box_damage_zone_ = from.oriented_generation_box_damage_zone_;

	density_property_name_ = from.density_property_name_;

	constraint_property_name_ = from.constraint_property_name_;

	is_edit_ = from.is_edit_;

    is_imported_ = from.is_imported_;
	is_created_ = from.is_created_;


	distribution_length1_type_index_ = from.distribution_length1_type_index_;
	distribution_length2_type_index_ = from.distribution_length2_type_index_;
	distribution_azimut_type_index_ = from.distribution_azimut_type_index_;
	distribution_angle_type_index_ = from.distribution_angle_type_index_;
	distribution_aperture_type_index_ = from.distribution_aperture_type_index_;
	
	distr_length1_val_ = from.distr_length1_val_;
	distr_length2_val_ = from.distr_length2_val_;
	distr_azimut_val_ = from.distr_azimut_val_;;
	distr_angle_val_ = from.distr_angle_val_;
	distr_aperture_val_ = from.distr_aperture_val_;
	

	density_type_index_ = from.density_type_index_;
	density_law_type_ = from.density_law_type_;

	to_generate_ = from.to_generate_;
	density_property_name_ = from.density_property_name_;
	length1_type_ = from.length1_type_;
	length2_type_= from.length2_type_;
	azimuth_type_= from.azimuth_type_;
	dip_angle_type_= from.dip_angle_type_;
	aperture_type_= from.aperture_type_;
	density_type_ = from.density_type_;
	minmaxmode_density_values_ = from.minmaxmode_density_values_;
	minmaxmode_shadow_zone_values_ = from.minmaxmode_shadow_zone_values_;

	minmaxmode_values_ = from.minmaxmode_values_;

	model_fractures_ = from.model_fractures_;

	fracture_set_index_=from.fracture_set_index_;

	default_stress_zone_max_ = from.default_stress_zone_max_;

	directory_path_ = from.directory_path_;

	is_shadow_zone_active_ = from.is_shadow_zone_active_;

	success_ = from.success_;

	cylinder_step_ = from.cylinder_step_;

	tab_type_intersection_ = from.tab_type_intersection_;
	type_category_ = from.type_category_;
	type_intersection_ = from.type_intersection_;
	array_index_ = from.array_index_;
	type_array_ = from.type_array_;

	echelon_space_ = from.echelon_space_;
	max_distance_connectivity_ = from.max_distance_connectivity_;
	nb_fractures_connection_relay_ = from.nb_fractures_connection_relay_;

	width_stress_zone_law_type_= from.width_stress_zone_law_type_;

	width_stress_zone_distribution_name_= from.width_stress_zone_distribution_name_;

	width_stress_zone_min_max_mode_ = from.width_stress_zone_min_max_mode_;

	width_stress_zone_property_name_ = from.width_stress_zone_property_name_;

	density_distribution_name_ = from.density_distribution_name_;
	density_min_max_ = from.density_min_max_;

	shadow_zone_distribution_name_ = from.shadow_zone_distribution_name_;
	shadow_zone_property_name_ = from.shadow_zone_property_name_;
	distribution_shadow_zone_type_index_ = from.distribution_shadow_zone_type_index_;
	distr_shadow_zone_val_ = from.distr_shadow_zone_val_;
	shadow_zone_val_ = from.shadow_zone_val_;
	shadow_zone_type_ = from.shadow_zone_type_;
	randShadowZone_ = from.randShadowZone_;

	randHangingWallDensity_ = from.randHangingWallDensity_;
	randFootWallDensity_ = from.randFootWallDensity_;

	fractures_number_hanging_wall_ = from.fractures_number_hanging_wall_;
	evaluate_fractures_number_hanging_wall_ = from.evaluate_fractures_number_hanging_wall_;
	fractures_activable_hanging_wall_ = from.fractures_activable_hanging_wall_;
	fractures_number_foot_wall_ = from.fractures_number_foot_wall_;
	evaluate_fractures_number_foot_wall_ = from.evaluate_fractures_number_foot_wall_;
	fractures_activable_foot_wall_ = from.fractures_activable_foot_wall_;

	density_hanging_wall_type_ = from.density_hanging_wall_type_;
	density_foot_wall_type_ = from.density_foot_wall_type_;

	density_hanging_wall_law_type_ = from.density_hanging_wall_law_type_;
	density_foot_wall_law_type_ = from.density_foot_wall_law_type_;

	minmaxmode_density_hanging_wall_values_ = from.minmaxmode_density_hanging_wall_values_;
	minmaxmode_density_foot_wall_values_ = from.minmaxmode_density_foot_wall_values_;
	distr_hanging_wall_val_ = from.distr_hanging_wall_val_;
	distr_foot_wall_val_ = from.distr_foot_wall_val_;

	density_hanging_wall_val_ = from.density_hanging_wall_val_;
	density_foot_wall_val_ = from.density_foot_wall_val_;

	density_hanging_wall_property_name_ = from.density_hanging_wall_property_name_;
	density_foot_wall_property_name_ = from.density_foot_wall_property_name_;
	density_hanging_wall_distribution_name_ = from.density_hanging_wall_distribution_name_;
	density_foot_wall_distribution_name_ = from.density_foot_wall_distribution_name_;

	distribution_echelon_type_index_ = from.distribution_echelon_type_index_;
	echelon_space_type_= from.echelon_space_type_;
	minmaxmode_echelon_space_values_ = from.minmaxmode_echelon_space_values_;
	randEchelon_space_ = from.randEchelon_space_;
	echelon_space_val_ = from.echelon_space_val_;
	echelon_space_distribution_name_ =  from.echelon_space_distribution_name_;

	is_active_ = from.is_active_;
	model_geo_cont_index_ = from.model_geo_cont_index_;
	output_prefixe_name_ = from.output_prefixe_name_;
	total_simulation_time_ = from.total_simulation_time_;
	total_writting_output_time_ = from.total_writting_output_time_;
	length_max_ = from.length_max_;
	randDistPowerLength1_ = from.randDistPowerLength1_;
	randDistPowerLength2_ = from.randDistPowerLength2_;
	randDistPowerAzimuth_ = from.randDistPowerAzimuth_;
	randDistPowerDipAngle_ = from.randDistPowerDipAngle_;
	randDistPowerAperture_ = from.randDistPowerAperture_;
	randDistPowerDensity_ = from.randDistPowerDensity_;
	randDistPowerMinimumSpace_ = from.randDistPowerMinimumSpace_;
	randDistPowerHangingWall_ = from.randDistPowerHangingWall_;
	randDistPowerFootWall_ = from.randDistPowerFootWall_;
 	randDistPowerEchelonSpace_ = from.randDistPowerEchelonSpace_;
 	randDistPowerShadowZone_ = from.randDistPowerShadowZone_;
 	randDistPowerAzimuthAnastomosing_ = from.randDistPowerAzimuthAnastomosing_;
}


//-------------------------------------------
// destructor of the class FractureSet
//-------------------------------------------

FractureSet::~FractureSet() 
{
}


//-----------------------------------------------
// assignment operator of the class FractureSet
//----------------------------------------------

FractureSet& FractureSet::operator=(const FractureSet& from)
{ 
	name_= from.name_;

	length_1_ = from.length_1_;

	length_2_  = from.length_2_;
	
	dip_azimut_  = from.dip_azimut_;

	dip_angle_  = from.dip_angle_;
	
	aperture_name_  = from.aperture_name_;

	density_ = from.density_;

	complementary_fracture_set_active_ = from.complementary_fracture_set_active_;

	density_function_= from.density_function_;

	coef_a_= from.coef_a_; 
		
	seed_= from.seed_;

	fractures_number_ = from.fractures_number_;

	evaluate_fractures_number_ = from.evaluate_fractures_number_;

	fractures_activable_ = from.fractures_activable_;

	fractures_inserted_ = from.fractures_inserted_;

//	fractures_ = from.fractures_;

	generation_box_ = from.generation_box_;

	oriented_generation_box_ = from.oriented_generation_box_;

	oriented_generation_box_damage_zone_ = from.oriented_generation_box_damage_zone_;

	density_property_name_ = from.density_property_name_;

	constraint_property_name_ = from.constraint_property_name_;

	is_edit_ = from.is_edit_;

    is_imported_ = from.is_imported_;

	is_created_ = from.is_created_;


	distribution_length1_type_index_ = from.distribution_length1_type_index_;
	distribution_length2_type_index_ = from.distribution_length2_type_index_;
	distribution_azimut_type_index_ = from.distribution_azimut_type_index_;
	distribution_angle_type_index_ = from.distribution_angle_type_index_;
	distribution_aperture_type_index_ = from.distribution_aperture_type_index_;
	density_type_index_ = from.density_type_index_;
	density_law_type_ = from.density_law_type_;
	distr_length1_val_ = from.distr_length1_val_;
	distr_length2_val_ = from.distr_length2_val_;
	distr_azimut_val_ = from.distr_azimut_val_;;
	distr_angle_val_ = from.distr_angle_val_;
	distr_aperture_val_ = from.distr_aperture_val_;
	to_generate_ = from.to_generate_;

	density_property_name_ = from.density_property_name_;
	length1_type_ = from.length1_type_;
	length2_type_= from.length2_type_;
	azimuth_type_= from.azimuth_type_;
	dip_angle_type_= from.dip_angle_type_;
	aperture_type_= from.aperture_type_;
	density_type_ = from.density_type_;
	minmaxmode_density_values_ = from.minmaxmode_density_values_;
	minmaxmode_shadow_zone_values_ = from.minmaxmode_shadow_zone_values_;

	minmaxmode_values_ = from.minmaxmode_values_;

	model_fractures_ = from.model_fractures_;

	fracture_set_index_=from.fracture_set_index_;

	default_stress_zone_max_ = from.default_stress_zone_max_;

	directory_path_ = from.directory_path_;

	is_shadow_zone_active_ = from.is_shadow_zone_active_;

	success_ = from.success_;

	cylinder_step_ = from.cylinder_step_;

	tab_type_intersection_ = from.tab_type_intersection_;

	type_category_ = from.type_category_;

	type_intersection_ = from.type_intersection_;
	array_index_ = from.array_index_;
	type_array_ = from.type_array_;
	echelon_space_ = from.echelon_space_;
	max_distance_connectivity_ = from.max_distance_connectivity_;
	nb_fractures_connection_relay_ = from.nb_fractures_connection_relay_;
	width_stress_zone_law_type_= from.width_stress_zone_law_type_;

	width_stress_zone_distribution_name_= from.width_stress_zone_distribution_name_;

	width_stress_zone_min_max_mode_ = from.width_stress_zone_min_max_mode_;

	width_stress_zone_property_name_ = from.width_stress_zone_property_name_;

	density_distribution_name_ = from.density_distribution_name_;

	density_min_max_ = from.density_min_max_;
	shadow_zone_distribution_name_ = from.shadow_zone_distribution_name_;
	shadow_zone_property_name_ = from.shadow_zone_property_name_;
	distribution_shadow_zone_type_index_ = from.distribution_shadow_zone_type_index_;
	distr_shadow_zone_val_ = from.distr_shadow_zone_val_;
	shadow_zone_val_ = from.shadow_zone_val_;
	shadow_zone_type_ = from.shadow_zone_type_;
	randShadowZone_ = from.randShadowZone_;

	fractures_number_hanging_wall_ = from.fractures_number_hanging_wall_;
	evaluate_fractures_number_hanging_wall_ = from.evaluate_fractures_number_hanging_wall_;
	fractures_activable_hanging_wall_ = from.fractures_activable_hanging_wall_;
	fractures_number_foot_wall_ = from.fractures_number_foot_wall_;
	evaluate_fractures_number_foot_wall_ = from.evaluate_fractures_number_foot_wall_;
	fractures_activable_foot_wall_ = from.fractures_activable_foot_wall_;

	density_hanging_wall_type_ = from.density_hanging_wall_type_;
	density_foot_wall_type_ = from.density_foot_wall_type_;

	density_hanging_wall_law_type_ = from.density_hanging_wall_law_type_;
	density_foot_wall_law_type_ = from.density_foot_wall_law_type_;

	minmaxmode_density_hanging_wall_values_ = from.minmaxmode_density_hanging_wall_values_;
	minmaxmode_density_foot_wall_values_ = from.minmaxmode_density_foot_wall_values_;
	distr_hanging_wall_val_ = from.distr_hanging_wall_val_;
	distr_foot_wall_val_ = from.distr_foot_wall_val_;

	density_hanging_wall_val_ = from.density_hanging_wall_val_;
	density_foot_wall_val_ = from.density_foot_wall_val_;

	density_hanging_wall_property_name_ = from.density_hanging_wall_property_name_;
	density_foot_wall_property_name_ = from.density_foot_wall_property_name_;
	density_hanging_wall_distribution_name_ = from.density_hanging_wall_distribution_name_;
	density_foot_wall_distribution_name_ = from.density_foot_wall_distribution_name_;
	randHangingWallDensity_ = from.randHangingWallDensity_;
	randFootWallDensity_ = from.randFootWallDensity_;

	distribution_echelon_type_index_ = from.distribution_echelon_type_index_;
	echelon_space_type_= from.echelon_space_type_;
	minmaxmode_echelon_space_values_ = from.minmaxmode_echelon_space_values_;
	randEchelon_space_ = from.randEchelon_space_;
	echelon_space_val_ = from.echelon_space_val_;
	echelon_space_distribution_name_ =  from.echelon_space_distribution_name_;
	is_active_ = from.is_active_;
	model_geo_cont_index_ = from.model_geo_cont_index_;
	output_prefixe_name_ = from.output_prefixe_name_;
	total_simulation_time_ = from.total_simulation_time_;
	total_writting_output_time_ = from.total_writting_output_time_;
	length_max_ = from.length_max_;
	randDistPowerLength1_ = from.randDistPowerLength1_;
	randDistPowerLength2_ = from.randDistPowerLength2_;
	randDistPowerAzimuth_ = from.randDistPowerAzimuth_;
	randDistPowerDipAngle_ = from.randDistPowerDipAngle_;
	randDistPowerAperture_ = from.randDistPowerAperture_;
	randDistPowerDensity_ = from.randDistPowerDensity_;
	randDistPowerMinimumSpace_ = from.randDistPowerMinimumSpace_;
	randDistPowerHangingWall_ = from.randDistPowerHangingWall_;
	randDistPowerFootWall_ = from.randDistPowerFootWall_;
 	randDistPowerEchelonSpace_ = from.randDistPowerEchelonSpace_;
 	randDistPowerShadowZone_ = from.randDistPowerShadowZone_;
 	randDistPowerAzimuthAnastomosing_ = from.randDistPowerAzimuthAnastomosing_;

	return *this; 
}


//--------------------------------------------------------
//      copy function 
//--------------------------------------------------------
  
void FractureSet::copy(const FractureSet& from)
{

	length_1_ = from.length_1_;

	length_2_  = from.length_2_;
	
	dip_azimut_  = from.dip_azimut_;

	dip_angle_  = from.dip_angle_;
	
	aperture_name_  = from.aperture_name_;

	density_ = from.density_;

	complementary_fracture_set_active_ = from.complementary_fracture_set_active_;

	density_function_= from.density_function_;

	coef_a_= from.coef_a_; 
		
	seed_= from.seed_;

	fractures_number_ = from.fractures_number_;

	fractures_activable_ = from.fractures_number_;

	fractures_inserted_ = 0;


	// tsurf copy 
	

	generation_box_ = from.generation_box_;

	oriented_generation_box_ = from.oriented_generation_box_;

	oriented_generation_box_damage_zone_ = from.oriented_generation_box_damage_zone_;

	density_property_name_ = from.density_property_name_;

	constraint_property_name_ = from.constraint_property_name_;

    is_imported_ = from.is_imported_;

	is_created_ = from.is_created_;


	distribution_length1_type_index_ = from.distribution_length1_type_index_;
	distribution_length2_type_index_ = from.distribution_length2_type_index_;
	distribution_azimut_type_index_ = from.distribution_azimut_type_index_;
	distribution_angle_type_index_ = from.distribution_angle_type_index_;
	distribution_aperture_type_index_ = from.distribution_aperture_type_index_;

	density_type_index_ = from.density_type_index_;
	density_law_type_ = from.density_law_type_;
	distr_length1_val_ = from.distr_length1_val_;
	distr_length2_val_ = from.distr_length2_val_;
	distr_azimut_val_ = from.distr_azimut_val_;
	distr_angle_val_ = from.distr_angle_val_;
	distr_aperture_val_ = from.distr_aperture_val_;
	to_generate_ = from.to_generate_;

	density_property_name_ = from.density_property_name_;
	length1_type_ = from.length1_type_;
	length2_type_= from.length2_type_;
	azimuth_type_= from.azimuth_type_;
	dip_angle_type_= from.dip_angle_type_;
	aperture_type_= from.aperture_type_;
	density_type_ = from.density_type_;
	minmaxmode_density_values_ = from.minmaxmode_density_values_;
	minmaxmode_shadow_zone_values_ = from.minmaxmode_shadow_zone_values_;

	minmaxmode_values_ = from.minmaxmode_values_;

	model_fractures_ = from.model_fractures_;

	fracture_set_index_=from.fracture_set_index_;
	default_stress_zone_max_ = from.default_stress_zone_max_;
	directory_path_ = from.directory_path_;
	is_shadow_zone_active_ = from.is_shadow_zone_active_;
	success_ = from.success_;
	cylinder_step_ = from.cylinder_step_;
	tab_type_intersection_ = from.tab_type_intersection_;
	type_category_ = from.type_category_;

	type_intersection_ = from.type_intersection_;

	array_index_ = from.array_index_;
	type_array_ = from.type_array_;
	echelon_space_ = from.echelon_space_;
	max_distance_connectivity_ = from.max_distance_connectivity_;
	nb_fractures_connection_relay_ = from.nb_fractures_connection_relay_;
	width_stress_zone_law_type_= from.width_stress_zone_law_type_;

	width_stress_zone_distribution_name_= from.width_stress_zone_distribution_name_;

	width_stress_zone_min_max_mode_ = from.width_stress_zone_min_max_mode_;

	width_stress_zone_property_name_ = from.width_stress_zone_property_name_;

	density_distribution_name_ = from.density_distribution_name_;
	density_min_max_ = from.density_min_max_;

	shadow_zone_distribution_name_ = from.shadow_zone_distribution_name_;
	shadow_zone_property_name_ = from.shadow_zone_property_name_;
	distribution_shadow_zone_type_index_ = from.distribution_shadow_zone_type_index_;
	distr_shadow_zone_val_ = from.distr_shadow_zone_val_;
	shadow_zone_val_ = from.shadow_zone_val_;
	shadow_zone_type_ = from.shadow_zone_type_;
	randShadowZone_ = from.randShadowZone_;

	fractures_number_hanging_wall_ = from.fractures_number_hanging_wall_;
	evaluate_fractures_number_hanging_wall_ = from.evaluate_fractures_number_hanging_wall_;
	fractures_activable_hanging_wall_ = from.fractures_activable_hanging_wall_;
	fractures_number_foot_wall_ = from.fractures_number_foot_wall_;
	evaluate_fractures_number_foot_wall_ = from.evaluate_fractures_number_foot_wall_;
	fractures_activable_foot_wall_ = from.fractures_activable_foot_wall_;

	density_hanging_wall_type_ = from.density_hanging_wall_type_;
	density_foot_wall_type_ = from.density_foot_wall_type_;

	density_hanging_wall_law_type_ = from.density_hanging_wall_law_type_;
	density_foot_wall_law_type_ = from.density_foot_wall_law_type_;

	minmaxmode_density_hanging_wall_values_ = from.minmaxmode_density_hanging_wall_values_;
	minmaxmode_density_foot_wall_values_ = from.minmaxmode_density_foot_wall_values_;
	distr_hanging_wall_val_ = from.distr_hanging_wall_val_;
	distr_foot_wall_val_ = from.distr_foot_wall_val_;

	density_hanging_wall_val_ = from.density_hanging_wall_val_;
	density_foot_wall_val_ = from.density_foot_wall_val_;

	density_hanging_wall_property_name_ = from.density_hanging_wall_property_name_;
	density_foot_wall_property_name_ = from.density_foot_wall_property_name_;
	density_hanging_wall_distribution_name_ = from.density_hanging_wall_distribution_name_;
	density_foot_wall_distribution_name_ = from.density_foot_wall_distribution_name_;
	randHangingWallDensity_ = from.randHangingWallDensity_;
	randFootWallDensity_ = from.randFootWallDensity_;

	distribution_echelon_type_index_ = from.distribution_echelon_type_index_;
	echelon_space_type_= from.echelon_space_type_;
	minmaxmode_echelon_space_values_ = from.minmaxmode_echelon_space_values_;
	randEchelon_space_ = from.randEchelon_space_;
	echelon_space_val_ = from.echelon_space_val_;
	echelon_space_distribution_name_ =  from.echelon_space_distribution_name_;
	is_active_ = from.is_active_;
	model_geo_cont_index_ = from.model_geo_cont_index_;
	output_prefixe_name_ = from.output_prefixe_name_;
	total_simulation_time_ = from.total_simulation_time_;
	total_writting_output_time_ = from.total_writting_output_time_;
	length_max_ = from.length_max_;
	randDistPowerLength1_ = from.randDistPowerLength1_;
	randDistPowerLength2_ = from.randDistPowerLength2_;
	randDistPowerAzimuth_ = from.randDistPowerAzimuth_;
	randDistPowerDipAngle_ = from.randDistPowerDipAngle_;
	randDistPowerAperture_ = from.randDistPowerAperture_;
	randDistPowerDensity_ = from.randDistPowerDensity_;
	randDistPowerMinimumSpace_ = from.randDistPowerMinimumSpace_;
	randDistPowerHangingWall_ = from.randDistPowerHangingWall_;
	randDistPowerFootWall_ = from.randDistPowerFootWall_;
 	randDistPowerEchelonSpace_ = from.randDistPowerEchelonSpace_;
 	randDistPowerShadowZone_ = from.randDistPowerShadowZone_;
 	randDistPowerAzimuthAnastomosing_ = from.randDistPowerAzimuthAnastomosing_;

}

//--------------------------------------------------------
//      distribution from index
//--------------------------------------------------------

void FractureSet::set_distribution_type_from_index(int index, std::string name)
{
	//Switch
	    switch(index){
	       case 0:
	    	   length1_type_ = get_distribution_type(name);
	            break;
	       case 1:
	    	   length2_type_ =  get_distribution_type(name);
	            break;
	       case 2:
	    	   azimuth_type_ =  get_distribution_type(name);
	            break;
	       case 3:
	    	   dip_angle_type_ =  get_distribution_type(name);
	            break;
	       case 4:
	    	   aperture_type_ =  get_distribution_type(name);
	            break;
	       case 5:
	    	   density_type_ =  get_distribution_type(name);
	            break;
	       case 6:
	    	   shadow_zone_type_ =  get_distribution_type(name);
	            break;
	       case 7:
	    	   density_hanging_wall_type_ =  get_distribution_type(name);
	            break;
	       case 8:
	    	   density_foot_wall_type_ =  get_distribution_type(name);
	            break;
	       case 9:
	    	   echelon_space_type_ =  get_distribution_type(name);
	            break;
	        default:
	            break;
	    }

}

FractureSet::distribution_type FractureSet::get_distribution_type(std::string name)
{
	distribution_type  type=distribution_type::uniform;

	if(name == "triangular"){

		type =  distribution_type::triangulated;
	}else if( name == "normal"){
		type =  distribution_type::normal;
	}else if( name == "lognormal"){
		type =  distribution_type::lognormal;
	}else if( name == "power"){
		type =  distribution_type::power;
	}


	return type;

}
//------------------------------------------------------
// set distribution
//------------------------------------------------------
/*
void FractureSet::set_distribution(int index,boost::mt19937&  center_gen)
{
    //Switch
    switch(index){
       case 0:
    	   randLength1_ = get_distribution_from_bound(length1_type_,minmaxmode_values_[0],minmaxmode_values_[1],minmaxmode_values_[2],center_gen);
            break;
       case 1:
    	   randLength2_ = get_distribution_from_bound(length2_type_,minmaxmode_values_[3],minmaxmode_values_[4],minmaxmode_values_[5],center_gen);
            break;
       case 2:
    	   randAzimuth_ = get_distribution_from_bound(azimuth_type_,minmaxmode_values_[6],minmaxmode_values_[7],minmaxmode_values_[8],center_gen);
            break;
       case 3:
    	   randDipAngle_ = get_distribution_from_bound(dip_angle_type_,minmaxmode_values_[9],minmaxmode_values_[10],minmaxmode_values_[11],center_gen);
            break;
       case 4:
    	   randAperture_ = get_distribution_from_bound(aperture_type_,minmaxmode_values_[12],minmaxmode_values_[13],minmaxmode_values_[14],center_gen);
            break;
        default:
            break;
    }
}
*/
//------------------------------------------------------
// register in the geobase the TSurf of the FractureSet
//------------------------------------------------------

/*  boost::variate_generator<boost::mt19937&, boost::uniform_real< double > > FractureSet :: get_uniform_dist(double min, double max, boost::mt19937& center_gen)
{
	 boost::uniform_real< double > test_l1(min, max);
	 static boost::variate_generator<boost::mt19937&, boost::uniform_real< double > > rand(center_gen, test_l1);
    return rand;
}
*/

//------------------------------------------------------
// register in the geobase the TSurf of the FractureSet
//------------------------------------------------------

void FractureSet::show_graphics()
{	
}


//--------------------------
// set the domain limits
//--------------------------

void FractureSet::set_generation_box(geode::BoundingBox3D box)
{
	generation_box_ = box;
}
//--------------------------
// set the domain limits
//--------------------------

void FractureSet::set_generation_oriented_box(fractures_intersect::oriented_bounding_box box)
{
	oriented_generation_box_ = box;
}



void FractureSet::update_bounding_box_for_fault_zone( double damage_zone_width)
{

	oriented_generation_box_damage_zone_ = new fractures_intersect::oriented_bounding_box(oriented_generation_box_.center_,oriented_generation_box_.R_,oriented_generation_box_.S_,oriented_generation_box_.T_,oriented_generation_box_.extension_);

	oriented_generation_box_damage_zone_->extension_[2]= damage_zone_width;

	volume_ = oriented_generation_box_damage_zone_->volume();
	area_ = oriented_generation_box_damage_zone_->area();

}
/**
 * @fn void initialize_shadow_zone_distribution(boost::mt19937&)
 * @brief
 *
 * @param center_gen
 */
void FractureSet::initialize_shadow_zone_distribution(boost::mt19937 &center_gen) {
	if (shadow_zone_type_ == distribution_type::uniform) {
	randShadowZone_ = boost::bind(boost::random::uniform_real_distribution<>(minmaxmode_shadow_zone_values_[0], minmaxmode_shadow_zone_values_[1]),center_gen);
	}
	else if (shadow_zone_type_ == distribution_type::triangulated) {
		randShadowZone_ = boost::bind(
				boost::random::triangle_distribution<>(minmaxmode_shadow_zone_values_[0],
						minmaxmode_shadow_zone_values_[2], minmaxmode_shadow_zone_values_[1]),
				center_gen);
	} else if (shadow_zone_type_ == distribution_type::normal) {
		randShadowZone_ = boost::bind(
				boost::random::normal_distribution<>(minmaxmode_shadow_zone_values_[0],
						minmaxmode_shadow_zone_values_[1]), center_gen);
	} else if (shadow_zone_type_ == distribution_type::lognormal) {
		randShadowZone_ = boost::bind(
				boost::random::lognormal_distribution<>(minmaxmode_shadow_zone_values_[0],
						minmaxmode_shadow_zone_values_[1]), center_gen);
	}else if (shadow_zone_type_ == distribution_type::power){
		randDistPowerShadowZone_ = DistributionPower(minmaxmode_shadow_zone_values_[0],minmaxmode_shadow_zone_values_[1]);
	}
}
/**
 * @fn void initialize_density_distribution(boost::mt19937&)
 * @brief
 *
 * @param center_gen
 */
void FractureSet::initialize_density_distribution(boost::mt19937 &center_gen) {
	if (density_type_ == distribution_type::uniform) {
	randDensity_ = boost::bind(boost::random::uniform_real_distribution<>(minmaxmode_density_values_[0], minmaxmode_density_values_[1]),center_gen);
	}
	else if (density_type_ == distribution_type::triangulated) {
		randDensity_ = boost::bind(
				boost::random::triangle_distribution<>(minmaxmode_density_values_[0],
						minmaxmode_density_values_[2], minmaxmode_density_values_[1]),
				center_gen);
	} else if (density_type_ == distribution_type::normal) {
		randDensity_ = boost::bind(
				boost::random::normal_distribution<>(minmaxmode_density_values_[0],
						minmaxmode_density_values_[1]), center_gen);
	} else if (density_type_ == distribution_type::lognormal) {
		randDensity_ = boost::bind(
				boost::random::lognormal_distribution<>(minmaxmode_density_values_[0],
						minmaxmode_density_values_[1]), center_gen);
	}else if (density_type_ == distribution_type::power){
		randDistPowerDensity_ = DistributionPower(minmaxmode_density_values_[0],minmaxmode_density_values_[1]);
	}
}

/**
 * @fn void initialize_echelon_density_distribution(boost::mt19937&)
 * @brief
 *
 * @param center_gen
 */
void FractureSet::initialize_echelon_space_density_distribution(boost::mt19937 &center_gen) {
	if (echelon_space_type_ == distribution_type::uniform) {
	randEchelon_space_ = boost::bind(boost::random::uniform_real_distribution<>(minmaxmode_echelon_space_values_[0], minmaxmode_echelon_space_values_[1]),center_gen);
	} else if (echelon_space_type_ == distribution_type::triangulated) {
		randEchelon_space_ = boost::bind(
				boost::random::triangle_distribution<>(minmaxmode_echelon_space_values_[0],
						minmaxmode_echelon_space_values_[2], minmaxmode_echelon_space_values_[1]),
				center_gen);
	} else if (echelon_space_type_ == distribution_type::normal) {
		randEchelon_space_ = boost::bind(
				boost::random::normal_distribution<>(minmaxmode_echelon_space_values_[0],
						minmaxmode_echelon_space_values_[1]), center_gen);
	} else if (echelon_space_type_ == distribution_type::lognormal) {
		randEchelon_space_ = boost::bind(
				boost::random::lognormal_distribution<>(minmaxmode_echelon_space_values_[0],
						minmaxmode_echelon_space_values_[1]), center_gen);
	}else if (echelon_space_type_ == distribution_type::power){
		randDistPowerEchelonSpace_ = DistributionPower(minmaxmode_echelon_space_values_[0],minmaxmode_echelon_space_values_[1]);
	}
}

/**
 * @fn void initialize_density_distribution(boost::mt19937&)
 * @brief
 *
 * @param center_gen
 */
void FractureSet::initialize_hanging_wall_density_distribution(boost::mt19937 &center_gen) {
	if (density_hanging_wall_type_ == distribution_type::uniform) {
	randHangingWallDensity_ = boost::bind(boost::random::uniform_real_distribution<>(minmaxmode_density_hanging_wall_values_[0], minmaxmode_density_hanging_wall_values_[1]),center_gen);
	} else if (density_hanging_wall_type_ == distribution_type::triangulated) {
		randHangingWallDensity_ = boost::bind(
				boost::random::triangle_distribution<>(minmaxmode_density_hanging_wall_values_[0],
						minmaxmode_density_hanging_wall_values_[2], minmaxmode_density_hanging_wall_values_[1]),
				center_gen);
	} else if (density_hanging_wall_type_ == distribution_type::normal) {
		randHangingWallDensity_ = boost::bind(
				boost::random::normal_distribution<>(minmaxmode_density_hanging_wall_values_[0],
						minmaxmode_density_hanging_wall_values_[1]), center_gen);
	} else if (density_hanging_wall_type_ == distribution_type::lognormal) {
		randHangingWallDensity_ = boost::bind(
				boost::random::lognormal_distribution<>(minmaxmode_density_hanging_wall_values_[0],
						minmaxmode_density_hanging_wall_values_[1]), center_gen);
	}else if (density_hanging_wall_type_ == distribution_type::power){
		randDistPowerHangingWall_ = DistributionPower(minmaxmode_density_hanging_wall_values_[0],minmaxmode_density_hanging_wall_values_[1]);
	}
}



/**
 * @fn void initialize_density_distribution(boost::mt19937&)
 * @brief
 *
 * @param center_gen
 */
void FractureSet::initialize_foot_wall_density_distribution(boost::mt19937 &center_gen) {
	if (density_foot_wall_type_ == distribution_type::uniform) {
	randFootWallDensity_ = boost::bind(boost::random::uniform_real_distribution<>(minmaxmode_density_foot_wall_values_[0], minmaxmode_density_foot_wall_values_[1]),center_gen);
	} else if (density_foot_wall_type_ == distribution_type::triangulated) {
		randFootWallDensity_ = boost::bind(
				boost::random::triangle_distribution<>(minmaxmode_density_foot_wall_values_[0],
						minmaxmode_density_foot_wall_values_[2], minmaxmode_density_foot_wall_values_[1]),
				center_gen);
	} else if (density_foot_wall_type_ == distribution_type::normal) {
		randFootWallDensity_ = boost::bind(
				boost::random::normal_distribution<>(minmaxmode_density_foot_wall_values_[0],
						minmaxmode_density_foot_wall_values_[1]), center_gen);
	} else if (density_foot_wall_type_ == distribution_type::lognormal) {
		randFootWallDensity_ = boost::bind(
				boost::random::lognormal_distribution<>(minmaxmode_density_foot_wall_values_[0],
						minmaxmode_density_foot_wall_values_[1]), center_gen);
	}else if (density_foot_wall_type_ == distribution_type::power){
		randDistPowerFootWall_ = DistributionPower(minmaxmode_density_foot_wall_values_[0],minmaxmode_density_foot_wall_values_[1]);
	}
}



/**
 * @fn void initialize_length1_distribution(boost::mt19937&)
 * @brief
 *
 * @param center_gen
 */
void FractureSet::initialize_length1_distribution(boost::mt19937 &center_gen) {
	if (length1_type_ == distribution_type::uniform){
		randLength1_ = boost::bind(boost::random::uniform_real_distribution<>(minmaxmode_values_[0], minmaxmode_values_[1]),center_gen);
	}
	else if (length1_type_ == distribution_type::triangulated) {
		randLength1_ = boost::bind(
				boost::random::triangle_distribution<>(minmaxmode_values_[0],
						minmaxmode_values_[2], minmaxmode_values_[1]),
				center_gen);
	} else if (length1_type_ == distribution_type::normal) {
		randLength1_ = boost::bind(
				boost::random::normal_distribution<>(minmaxmode_values_[0],
						minmaxmode_values_[1]), center_gen);
	} else if (length1_type_ == distribution_type::lognormal) {
		randLength1_ = boost::bind(
				boost::random::lognormal_distribution<>(minmaxmode_values_[0],
						minmaxmode_values_[1]), center_gen);
	}else if (length1_type_ == distribution_type::power){
		randDistPowerLength1_ = DistributionPower(minmaxmode_values_[0],minmaxmode_values_[1]);
	}
}
/**
 * @fn void initialize_lenght2_distribution(boost::mt19937&)
 * @brief
 *
 * @param center_gen
 */
void FractureSet::initialize_lenght2_distribution( boost::mt19937 &center_gen) {
	if (length2_type_ == distribution_type::uniform){
		randLength2_ = boost::bind(
			boost::random::uniform_real_distribution<>(minmaxmode_values_[3], minmaxmode_values_[4]),
			center_gen);
   }else if (length2_type_ == distribution_type::triangulated) {
		randLength2_ = boost::bind(
				boost::random::triangle_distribution<>(minmaxmode_values_[3],
						minmaxmode_values_[5], minmaxmode_values_[4]),
				center_gen);
	} else if (length2_type_ == distribution_type::normal) {
		randLength2_ = boost::bind(
				boost::random::normal_distribution<>(minmaxmode_values_[3],
						minmaxmode_values_[4]), center_gen);
	} else if (length2_type_ == distribution_type::lognormal) {
		randLength2_ = boost::bind(
				boost::random::lognormal_distribution<>(minmaxmode_values_[3],
						minmaxmode_values_[4]), center_gen);
	}else if (length2_type_ == distribution_type::power){
		randDistPowerLength2_ = DistributionPower(minmaxmode_values_[3],minmaxmode_values_[4]);
	}
}
/**
 * @fn void initialize_azimut_distribution(boost::mt19937&)
 * @brief
 *
 * @param center_gen
 */
void FractureSet::initialize_azimut_distribution( boost::mt19937 &center_gen) {


	if(type_array_ == array_type::anastomosing){

		if(minmaxmode_values_[6] < 10) minmaxmode_values_[6]=10;
		if(minmaxmode_values_[7] <= (minmaxmode_values_[6] + 10) ) minmaxmode_values_[7] = 20;
		if (azimuth_type_ == distribution_type::uniform) {
		randAzimut_ = boost::bind(boost::random::uniform_real_distribution<>(minmaxmode_values_[6], minmaxmode_values_[7]),center_gen);
		randAzimutAnastomosing_= boost::bind(boost::random::uniform_real_distribution<>(360-minmaxmode_values_[7], 350),center_gen);

		} else if (azimuth_type_ == distribution_type::triangulated) {
				randAzimut_ = boost::bind(boost::random::triangle_distribution<>(minmaxmode_values_[6],minmaxmode_values_[8], minmaxmode_values_[7]),center_gen);
				randAzimutAnastomosing_= boost::bind(boost::random::triangle_distribution<>(360-minmaxmode_values_[7],minmaxmode_values_[8], 350),center_gen);
		} else if (azimuth_type_ == distribution_type::normal) {
			randAzimut_ = boost::bind(boost::random::normal_distribution<>(minmaxmode_values_[6],minmaxmode_values_[7]), center_gen);
			randAzimutAnastomosing_=boost::bind(boost::random::normal_distribution<>(360-minmaxmode_values_[7],350), center_gen);
		} else if (azimuth_type_ == distribution_type::lognormal) {
			randAzimut_ = boost::bind(boost::random::lognormal_distribution<>(minmaxmode_values_[6],minmaxmode_values_[7]), center_gen);
			randAzimutAnastomosing_= boost::bind(boost::random::lognormal_distribution<>(360-minmaxmode_values_[7],350), center_gen);
		}else if (azimuth_type_ == distribution_type::power){
			randDistPowerAzimuth_ = DistributionPower(minmaxmode_values_[6],minmaxmode_values_[7]);
			randDistPowerAzimuthAnastomosing_ = DistributionPower(360-minmaxmode_values_[7],350);
		}

	}if(type_array_ == array_type::corridor){


		if (azimuth_type_ == distribution_type::uniform) {
		randAzimut_ = boost::bind(boost::random::uniform_real_distribution<>(minmaxmode_values_[6], minmaxmode_values_[7]),center_gen);
		}


	}else{
		if (azimuth_type_ == distribution_type::uniform){
			randAzimut_ = boost::bind(boost::random::uniform_real_distribution<>(minmaxmode_values_[6], minmaxmode_values_[7]),center_gen);
		}
		else if (azimuth_type_ == distribution_type::triangulated) {
			randAzimut_ = boost::bind(boost::random::triangle_distribution<>(minmaxmode_values_[6],minmaxmode_values_[8], minmaxmode_values_[7]),center_gen);
		} else if (azimuth_type_ == distribution_type::normal) {
			randAzimut_ = boost::bind(boost::random::normal_distribution<>(minmaxmode_values_[6],minmaxmode_values_[7]), center_gen);
		} else if (azimuth_type_ == distribution_type::lognormal) {
			randAzimut_ = boost::bind(boost::random::lognormal_distribution<>(minmaxmode_values_[6],minmaxmode_values_[7]), center_gen);
		}else if (azimuth_type_ == distribution_type::power){
			randDistPowerAzimuth_ = DistributionPower(minmaxmode_values_[6],minmaxmode_values_[7]);
		}
	}



}
/**
 * @fn void initialize_diap_angle_distribution(boost::mt19937&)
 * @brief
 *
 * @param center_gen
 */
void FractureSet::initialize_dip_angle_distribution( boost::mt19937 &center_gen) {
	if (dip_angle_type_ == distribution_type::uniform) {
	randDip_ = boost::bind(
			boost::random::uniform_real_distribution<>(minmaxmode_values_[9],
					minmaxmode_values_[10]), center_gen);
	}
	else if (dip_angle_type_ == distribution_type::triangulated) {
		randDip_ = boost::bind(
				boost::random::triangle_distribution<>(minmaxmode_values_[9],
						minmaxmode_values_[11], minmaxmode_values_[10]),
				center_gen);
	} else if (dip_angle_type_ == distribution_type::normal) {
		randDip_ = boost::bind(
				boost::random::normal_distribution<>(minmaxmode_values_[9],
						minmaxmode_values_[10]), center_gen);
	} else if (dip_angle_type_ == distribution_type::lognormal) {
		randDip_ = boost::bind(
				boost::random::lognormal_distribution<>(minmaxmode_values_[9],
						minmaxmode_values_[10]), center_gen);
	}else if (dip_angle_type_ == distribution_type::power){
		randDistPowerDipAngle_ = DistributionPower(minmaxmode_values_[9],minmaxmode_values_[10]);
	}
}
/**
 * @fn void initialize_aperture_distribution(boost::mt19937&)
 * @brief
 *
 * @param center_gen
 */
void FractureSet::initialize_aperture_distribution( boost::mt19937 &center_gen) {
	if (aperture_type_ == distribution_type::uniform) {
	randAperture_ = boost::bind(
			boost::random::uniform_real_distribution<>(minmaxmode_values_[12],
					minmaxmode_values_[13]), center_gen);
	}
	else if (aperture_type_ == distribution_type::triangulated) {
		randAperture_ = boost::bind(
				boost::random::triangle_distribution<>(minmaxmode_values_[12],
						minmaxmode_values_[14], minmaxmode_values_[13]),
				center_gen);
	} else if (aperture_type_ == distribution_type::normal) {
		randAperture_ = boost::bind(
				boost::random::normal_distribution<>(minmaxmode_values_[12],
						minmaxmode_values_[13]), center_gen);
	} else if (aperture_type_ == distribution_type::lognormal) {
		randAperture_ = boost::bind(
				boost::random::lognormal_distribution<>(minmaxmode_values_[12],
						minmaxmode_values_[13]), center_gen);
	}else if (aperture_type_ == distribution_type::power){
		randDistPowerAperture_ = DistributionPower(minmaxmode_values_[12],minmaxmode_values_[13]);
	}
}

void FractureSet::initialize_fracture_set_distribution(boost::mt19937& center_gen) {
	initialize_length1_distribution(center_gen);
	initialize_lenght2_distribution(center_gen);
	initialize_azimut_distribution(center_gen);
	initialize_dip_angle_distribution(center_gen);
	initialize_aperture_distribution(center_gen);
	initialize_density_distribution(center_gen);
	initialize_echelon_space_density_distribution(center_gen);
	initialize_hanging_wall_density_distribution(center_gen);
	initialize_foot_wall_density_distribution(center_gen);
	initialize_shadow_zone_distribution(center_gen);
}

void FractureSet::clean_vertices_and_triangles_array() {
	vertices_.clear();
	triangles_.clear();
	triangles_l1_.clear();
	triangles_l2_.clear();
	triangles_orientation_.clear();
	triangles_dip_.clear();
	triangles_aperture_.clear();
}

float FractureSet::evaluate_fracture_number_from_oriented_box(const geode::StructuralModel &model,std::vector<CellId>& cells,bool density_variable_in_grid){

	float depthmax=0.0;
	float depth;
	float ntot;
	float ntot_P31;
	float theta_max=0.0;
	float theta_max_P31=0.0;
	float theta_depth;
	float theta_depth_P31;
	theta_max = 1.0;
	theta_max_P31 = 1.0;

	depthmax = oriented_generation_box_.extension_.z();

	// determine the mean density, the volume and the area of the model


	int nb_frac;

	nb_frac = evaluate_fractures_number_;



	double l1_mean = (minmaxmode_values_[1] + minmaxmode_values_[0])/2;
	double l2_mean = (minmaxmode_values_[3] + minmaxmode_values_[4])/2;

	double surface_mean = l1_mean*l2_mean;

	// case with density variable in a grid
	if (density_variable_in_grid == true)
	{
		generate_fractures_density_in_array(density_variable_in_grid,density_,model,cells,surface_mean,density_law_type_);
		if(density_law_type_ ==1){
			ntot = density_;
		}else if(density_law_type_ ==2){

		    ntot_P31 = density_;
		}
	}
	//  other cases
	else
	{

		// determine the mean thickness of the sgrid
		float mean_thickness= volume_/area_;
		float mean_thickness_P31= volume_/area_;

		if( coef_a_ < 0.0)
		{
			theta_depth= mean_thickness *coef_a_ + density_;
			theta_depth_P31 =  mean_thickness *coef_a_ + density_/surface_mean;

			if( theta_depth < 0.0) mean_thickness=-density_/coef_a_;
			if( theta_depth_P31 < 0.0) mean_thickness_P31=-density_/(surface_mean*coef_a_);
		}

		ntot=((coef_a_/2)*((mean_thickness)*(mean_thickness)) + density_*(mean_thickness))*area_;
		ntot_P31 = ((coef_a_/2)*((mean_thickness_P31)*(mean_thickness_P31)) + (density_/surface_mean)*(mean_thickness_P31))*area_;


		if (coef_a_ > 0.0)
		{
			theta_max = coef_a_*depthmax + density_;
			theta_max_P31 = coef_a_*depthmax + density_/surface_mean;
		}
		else
		{
			theta_max=density_;
			theta_max_P31 = density_/surface_mean;
		}


		float min_v = std::numeric_limits<int>::min();
		float max_v = std::numeric_limits<int>::max();

		if(density_law_type_==1){
			if (ntot < min_v || ntot > max_v  )
			{
				if (ntot > max_v) evaluate_fractures_number_ = std::numeric_limits<int>::max();
				if (ntot < min_v) evaluate_fractures_number_ = std::numeric_limits<int>::min();
				return theta_max;
			}
		}else if(density_law_type_==2) {
			if (ntot_P31 < min_v || ntot_P31 > max_v  )
			{
				if (ntot_P31 > max_v) evaluate_fractures_number_ = std::numeric_limits<int>::max();
				if (ntot_P31 < min_v) evaluate_fractures_number_ = std::numeric_limits<int>::min();
				return theta_max_P31;
			}
		}
	}


	if(density_law_type_==1){
	  evaluate_fractures_number_ = int(ntot);

	} else if(density_law_type_==2) {
		 evaluate_fractures_number_ = int(ntot_P31);
		 theta_max = theta_max_P31;
	}else{
		theta_max=1.0;
	}


	return theta_max;
}

float FractureSet::evaluate_fracture_number_from_oriented_box_foot_wall_zone(const geode::StructuralModel &model,std::vector<CellId>& cells,bool density_variable_in_grid, double ratio_footwall){

	float depthmax;
	float depth;
	float ntot;
	float ntot_P31;

	float theta_max;
	float theta_max_P31;

	float theta_depth;
	float theta_depth_P31;


	theta_max = 1.0;
	theta_max_P31= 1.0;

	depthmax = oriented_generation_box_.extension_.z()*ratio_footwall;

	// determine the mean density, the volume and the area of the model


	int nb_frac;

	nb_frac = evaluate_fractures_number_foot_wall_;



	double l1_mean = (minmaxmode_density_foot_wall_values_[1] + minmaxmode_density_foot_wall_values_[0])/2;
	double l2_mean = (minmaxmode_density_foot_wall_values_[3] + minmaxmode_density_foot_wall_values_[4])/2;


	double surface_mean = l1_mean*l2_mean;

	// case with density variable in a grid
	if (density_variable_in_grid == true)
	{
		generate_fractures_density_in_array(density_variable_in_grid,density_,model,cells,surface_mean,density_law_type_);
		if(density_law_type_ ==1){
			ntot = int(volume_*ratio_footwall * density_);
		}else if(density_law_type_ ==2){

			ntot_P31 = int(volume_ *ratio_footwall* density_);
		}



	}
	//  other cases
	else
	{

		// determine the mean thickness of the sgrid
		float mean_thickness= volume_*ratio_footwall/area_;
		float mean_thickness_P31= volume_*ratio_footwall/area_;

		if( coef_a_ < 0.0)
		{
			theta_depth= mean_thickness *coef_a_ + density_foot_wall_;
			theta_depth_P31 =  mean_thickness_P31 *coef_a_ + density_foot_wall_/surface_mean;

			if( theta_depth < 0.0) mean_thickness=-density_foot_wall_/coef_a_;
			if( theta_depth_P31 < 0.0) mean_thickness_P31=-density_foot_wall_/(surface_mean*coef_a_);


		}

		ntot=((coef_a_/2)*((mean_thickness)*(mean_thickness)) + density_foot_wall_*(mean_thickness))*area_;
		ntot_P31 = ((coef_a_/2)*((mean_thickness_P31)*(mean_thickness_P31)) + (density_foot_wall_/surface_mean)*(mean_thickness_P31))*area_;


		if (coef_a_ > 0.0)
		{
			theta_max = coef_a_*depthmax + density_foot_wall_;
			theta_max_P31 = coef_a_*depthmax + density_foot_wall_/surface_mean;


		}
		else
		{
			theta_max=density_foot_wall_;
			theta_max_P31 = density_foot_wall_/surface_mean;

		}

		float min_v = std::numeric_limits<int>::min();
		float max_v = std::numeric_limits<int>::max();

		if(density_foot_wall_law_type_==1){
			if (ntot < min_v || ntot > max_v  )
			{
				if (ntot > max_v) evaluate_fractures_number_foot_wall_ = std::numeric_limits<int>::max();
				if (ntot < min_v) evaluate_fractures_number_foot_wall_ = std::numeric_limits<int>::min();
				return theta_max;
			}



		}else if(density_foot_wall_law_type_==2) {
			if (ntot_P31 < min_v || ntot_P31 > max_v  )
			{
				if (ntot_P31 > max_v) evaluate_fractures_number_foot_wall_ = std::numeric_limits<int>::max();
				if (ntot_P31 < min_v) evaluate_fractures_number_foot_wall_ = std::numeric_limits<int>::min();
				return theta_max_P31;
			}
		}




	}


	if(density_foot_wall_law_type_==1){
		evaluate_fractures_number_foot_wall_ = int(ntot);

	} else if(density_foot_wall_law_type_==2) {
		evaluate_fractures_number_foot_wall_ = int(ntot_P31);
		 theta_max = theta_max_P31;
	}else{
		theta_max=1.0;
	}


	return theta_max;
}

float FractureSet::evaluate_fracture_number_from_oriented_box_hanging_wall_zone(const geode::StructuralModel &model,std::vector<CellId>& cells,bool density_variable_in_grid, double ratio_hanging){

	float depthmax;
	float depth;
	float ntot;
	float ntot_P31;

	float theta_max;
	float theta_max_P31;

	float theta_depth;
	float theta_depth_P31;


	theta_max = 1.0;
	theta_max_P31= 1.0;

	depthmax = oriented_generation_box_.extension_.z()*ratio_hanging;

	// determine the mean density, the volume and the area of the model


	int nb_frac;

	nb_frac = evaluate_fractures_number_hanging_wall_;



	double l1_mean = (minmaxmode_values_[1] + minmaxmode_values_[0])/2;
	double l2_mean = (minmaxmode_values_[3] + minmaxmode_values_[4])/2;


	double surface_mean = l1_mean*l2_mean;

	// case with density variable in a grid
	if (density_variable_in_grid == true)
	{
		generate_fractures_density_in_array(density_variable_in_grid,density_,model,cells,surface_mean,density_law_type_);
		ntot = int(volume_*ratio_hanging * density_hanging_wall_);
		ntot_P31 = int(volume_ *ratio_hanging* density_hanging_wall_/surface_mean);

	}
	//  other cases
	else
	{

		// determine the mean thickness of the sgrid
		float mean_thickness= volume_*ratio_hanging/area_;
		float mean_thickness_P31= volume_*ratio_hanging/area_;

		if( coef_a_ < 0.0)
		{
			theta_depth= mean_thickness *coef_a_ + density_hanging_wall_;
			theta_depth_P31 =  mean_thickness_P31 *coef_a_ + density_hanging_wall_/surface_mean;

			if( theta_depth < 0.0) mean_thickness=-density_hanging_wall_/coef_a_;
			if( theta_depth_P31 < 0.0) mean_thickness_P31=-density_hanging_wall_/(surface_mean*coef_a_);


		}

		ntot=((coef_a_/2)*((mean_thickness)*(mean_thickness)) + density_hanging_wall_*(mean_thickness))*area_;
		ntot_P31 = ((coef_a_/2)*((mean_thickness_P31)*(mean_thickness_P31)) + (density_hanging_wall_/surface_mean)*(mean_thickness_P31))*area_;


		if (coef_a_ > 0.0)
		{
			theta_max = coef_a_*depthmax + density_hanging_wall_;
			theta_max_P31 = coef_a_*depthmax + density_hanging_wall_/surface_mean;


		}
		else
		{
			theta_max=density_hanging_wall_;
			theta_max_P31 = density_hanging_wall_/surface_mean;

		}

		float min_v = std::numeric_limits<int>::min();
		float max_v = std::numeric_limits<int>::max();

		if(density_hanging_wall_law_type_==1){
			if (ntot < min_v || ntot > max_v  )
			{
				if (ntot > max_v) evaluate_fractures_number_hanging_wall_ = std::numeric_limits<int>::max();
				if (ntot < min_v) evaluate_fractures_number_hanging_wall_ = std::numeric_limits<int>::min();
				return theta_max;
			}



		}else if(density_hanging_wall_law_type_==2) {
			if (ntot_P31 < min_v || ntot_P31 > max_v  )
			{
				if (ntot_P31 > max_v) evaluate_fractures_number_hanging_wall_ = std::numeric_limits<int>::max();
				if (ntot_P31 < min_v) evaluate_fractures_number_hanging_wall_ = std::numeric_limits<int>::min();
				return theta_max_P31;
			}
		}




	}


	if(density_hanging_wall_law_type_==1){
		evaluate_fractures_number_hanging_wall_ = int(ntot);

	} else if(density_hanging_wall_law_type_==2) {
		evaluate_fractures_number_hanging_wall_ = int(ntot_P31);
		 theta_max = theta_max_P31;
	}else{
		theta_max=1.0;
	}


	return theta_max;
}

float FractureSet::evaluate_fracture_number(const geode::StructuralModel &model,std::vector<CellId>& cells,bool density_variable_in_grid){
	
	float depthmax=0.0;
	float depth;
	float ntot;
	float ntot_P31;
	float theta_max=0.0;
	float theta_max_P31=0.0;
	float theta_depth;
	float theta_depth_P31;
	theta_max = 1.0;
	theta_max_P31 = 1.0;
	geode::Point3D  min_box = generation_box_.min();
	geode::Point3D  max_box = generation_box_.max();

	geode::Vector3D box_length(min_box, max_box);
	depthmax = box_length.value(2);

	// determine the mean density, the volume and the area of the model


	int nb_frac;

	nb_frac = evaluate_fractures_number_;



	double l1_mean = (minmaxmode_values_[1] + minmaxmode_values_[0])/2;
	double l2_mean = (minmaxmode_values_[3] + minmaxmode_values_[4])/2;



	if (distribution_length1_type_index_ == 1) {
		l1_mean = distr_length1_val_;
	}
	if (distribution_length2_type_index_ == 1) {
		l2_mean = distr_length2_val_;
	}


	double surface_mean = l1_mean*l2_mean;

	// case with density variable in a grid 	
	if (density_variable_in_grid == true) 
	{
		generate_fractures_density(density_variable_in_grid,density_,model,cells,surface_mean,density_law_type_);

		ntot = int(volume_ * density_);
		ntot_P31 = int(volume_ * density_/surface_mean);
	}
	//  other cases
	else
	{

		// determine the mean thickness of the sgrid 
		float mean_thickness= volume_/area_;
		float mean_thickness_P31= volume_/area_;

		if( coef_a_ < 0.0)
		{
			theta_depth= mean_thickness *coef_a_ + density_;
			theta_depth_P31 =  mean_thickness *coef_a_ + density_/surface_mean;

			if( theta_depth < 0.0) mean_thickness=-density_/coef_a_;
			if( theta_depth_P31 < 0.0) mean_thickness_P31=-density_/(surface_mean*coef_a_);
		}

		ntot=((coef_a_/2)*((mean_thickness)*(mean_thickness)) + density_*(mean_thickness))*area_;
		ntot_P31 = ((coef_a_/2)*((mean_thickness_P31)*(mean_thickness_P31)) + (density_/surface_mean)*(mean_thickness_P31))*area_;


		if (coef_a_ > 0.0)
		{
			theta_max = coef_a_*depthmax + density_;
			theta_max_P31 = coef_a_*depthmax + density_/surface_mean;
		}
		else
		{
			theta_max=density_;
			theta_max_P31 = density_/surface_mean;
		}


		float min_v = std::numeric_limits<int>::min();
		float max_v = std::numeric_limits<int>::max();

		if(density_law_type_==1){
			if (ntot < min_v || ntot > max_v  )
			{
				if (ntot > max_v) evaluate_fractures_number_ = std::numeric_limits<int>::max();
				if (ntot < min_v) evaluate_fractures_number_ = std::numeric_limits<int>::min();
				return theta_max;
			}
		}else if(density_law_type_==2) {
			if (ntot_P31 < min_v || ntot_P31 > max_v  )
			{
				if (ntot_P31 > max_v) evaluate_fractures_number_ = std::numeric_limits<int>::max();
				if (ntot_P31 < min_v) evaluate_fractures_number_ = std::numeric_limits<int>::min();
				return theta_max_P31;
			}
		}
	}


	if(density_law_type_==1){
	  evaluate_fractures_number_ = int(ntot);

	} else if(density_law_type_==2) {
		 evaluate_fractures_number_ = int(ntot_P31);
		 theta_max = theta_max_P31;
	}else{
		theta_max=1.0;
	}


	return theta_max;
}

void FractureSet::evaluate_final_fracture_number(boost::uniform_01<boost::mt19937>& center_uni){
	// generation of number_frac fractures according to a Poisson 
	// distribution law where ntot
	// is the Poisson parameter equal to mean number of fractures
	int number_mean_frac= int(-evaluate_fractures_number_);
	float sum=0;
	int count=0;
	int number_frac=0;
	//  Initialisation of random parameter from uniform law





	while ( sum > number_mean_frac)
	{
		sum+=log(center_uni());
		count++;

	}	
	
	number_frac=std::max(count-1,0);

	set_fractures_number(number_frac);
}


void FractureSet::evaluate_final_fracture_number_hanging_wall(boost::uniform_01<boost::mt19937>& center_uni){
	// generation of number_frac fractures according to a Poisson
	// distribution law where ntot
	// is the Poisson parameter equal to mean number of fractures
	int number_mean_frac= int(-evaluate_fractures_number_hanging_wall_);
	float sum=0;
	int count=0;
	int number_frac=0;
	//  Initialisation of random parameter from uniform law





	while ( sum > number_mean_frac)
	{
		sum+=log(center_uni());
		count++;

	}

	number_frac=std::max(count-1,0);

	set_fractures_number_hanging_wall(number_frac);
}

void FractureSet::evaluate_final_fracture_number_foot_wall(boost::uniform_01<boost::mt19937>& center_uni){
	// generation of number_frac fractures according to a Poisson
	// distribution law where ntot
	// is the Poisson parameter equal to mean number of fractures
	int number_mean_frac= int(-evaluate_fractures_number_foot_wall_);
	float sum=0;
	int count=0;
	int number_frac=0;
	//  Initialisation of random parameter from uniform law





	while ( sum > number_mean_frac)
	{
		sum+=log(center_uni());
		count++;

	}

	number_frac=std::max(count-1,0);

	set_fractures_number_foot_wall(number_frac);
}



void FractureSet::initialize_geometrical_fractures_data() {
	// test if asociated tsurf object exist or not
	// if( fractures_ == NULL ) return;
	length1_.resize(fractures_number_, 0);
	length2_.resize(fractures_number_, 0);
	Orientation_.resize(fractures_number_, 0);
	Dip_.resize(fractures_number_, 0);
	aperture_.resize(fractures_number_, 0);
	density_val_.resize(fractures_number_, 0);
	shadow_zone_val_.resize(fractures_number_, 0);
	density_hanging_wall_val_.resize(fractures_number_, 0);
	density_foot_wall_val_.resize(fractures_number_, 0);
	//echelon_space_val_.resize(fractures_number_, 0);
}


void FractureSet::update_values_of_fracture_index_from_property(CellId id_cell,const geode::StructuralModel &model, int ifrac){
	//cout << randLength1_() << endl;
	const auto &block_mesh = model.block(id_cell.blockId).mesh<geode::TetrahedralSolid3D>();

	if (distribution_length1_type_index_ > 1) {

		bool testPoperty =block_mesh.polyhedron_attribute_manager().attribute_exists(length1_property_name_);

			if (!testPoperty) {
				cout << " The Property " << length1_property_name_
						<< " doesn't exist for block " << model.block(id_cell.blockId).name() << endl;

			}else{
			// get petrophysique properties from polyhedron
			const auto &attributPropertyLenght1_ =block_mesh.polyhedron_attribute_manager().find_or_create_attribute<geode::VariableAttribute, double>(length1_property_name_,ndvalue_);


		     length1_[ifrac] = attributPropertyLenght1_->value(id_cell.polyhedronIndex); //randLength1_();
			}
	}
	if (distribution_length2_type_index_ >1) {
		bool testPoperty =block_mesh.polyhedron_attribute_manager().attribute_exists(length2_property_name_);

			if (!testPoperty) {
				cout << " The Property " << length2_property_name_
						<< " doesn't exist for block " << model.block(id_cell.blockId).name() << endl;

			}else{
			// get petrophysique properties from polyhedron
			const auto &attributPropertyLength2_ =block_mesh.polyhedron_attribute_manager().find_or_create_attribute<geode::VariableAttribute, double>(length2_property_name_,ndvalue_);

		      length2_[ifrac] = attributPropertyLength2_->value(id_cell.polyhedronIndex);
			}
	}

	if (type_array_ == array_type::anastomosing) {
		if (distribution_azimut_type_index_ > 1) {
		    bool testPoperty =block_mesh.polyhedron_attribute_manager().attribute_exists(azimut_property_name_);

			if (!testPoperty) {
				cout << " The Property " << azimut_property_name_
						<< " doesn't exist for block " << model.block(id_cell.blockId).name() << endl;

			}else{
			// get petrophysique properties from polyhedron
			const auto &attributPropertyAzimut_ =block_mesh.polyhedron_attribute_manager().find_or_create_attribute<geode::VariableAttribute, double>(azimut_property_name_,ndvalue_);
			if (ifrac & 1) {

					Orientation_[ifrac] = 270 - attributPropertyAzimut_->value(id_cell.polyhedronIndex);


			} else {
				Orientation_[ifrac] = attributPropertyAzimut_->value(id_cell.polyhedronIndex);
			}
			}
		}

	} else if (type_array_ == array_type::corridor) {
		if (distribution_azimut_type_index_ > 1) {
		    bool testPoperty =block_mesh.polyhedron_attribute_manager().attribute_exists(azimut_property_name_);

			if (!testPoperty) {
				cout << " The Property " << azimut_property_name_
						<< " doesn't exist for block " << model.block(id_cell.blockId).name() << endl;

			}else {
			// get petrophysique properties from polyhedron
			const auto &attributPropertyAzimut_ =block_mesh.polyhedron_attribute_manager().find_or_create_attribute<geode::VariableAttribute, double>(azimut_property_name_,ndvalue_);
			Orientation_[ifrac] = attributPropertyAzimut_->value(id_cell.polyhedronIndex);
			}
		}
	} else {
		if (distribution_azimut_type_index_ > 1) {
		    bool testPoperty =block_mesh.polyhedron_attribute_manager().attribute_exists(azimut_property_name_);

			if (!testPoperty) {
				cout << " The Property " << azimut_property_name_
						<< " doesn't exist for block " << model.block(id_cell.blockId).name() << endl;

			}else{
			// get petrophysique properties from polyhedron
				const auto &attributPropertyAzimut_ =block_mesh.polyhedron_attribute_manager().find_or_create_attribute<geode::VariableAttribute, double>(azimut_property_name_,ndvalue_);
				Orientation_[ifrac] = attributPropertyAzimut_->value(id_cell.polyhedronIndex);
			}
		}
	}
	if (distribution_angle_type_index_ > 1) {
	    bool testPoperty =block_mesh.polyhedron_attribute_manager().attribute_exists(dip_angle_property_name_);

		if (!testPoperty) {
			cout << " The Property " << dip_angle_property_name_
					<< " doesn't exist for block " << model.block(id_cell.blockId).name() << endl;

		}else{
			// get petrophysique properties from polyhedron
			const auto &attributPropertyDipAngle =block_mesh.polyhedron_attribute_manager().find_or_create_attribute<geode::VariableAttribute, double>(dip_angle_property_name_,ndvalue_);
			Dip_[ifrac] = attributPropertyDipAngle->value(id_cell.polyhedronIndex);
		}
	}
	if (distribution_aperture_type_index_ > 1) {
		bool testPoperty =block_mesh.polyhedron_attribute_manager().attribute_exists(aperture_property_name_);

		if (!testPoperty) {
			cout << " The Property " << aperture_property_name_
					<< " doesn't exist for block " << model.block(id_cell.blockId).name() << endl;
		}else{
			// get petrophysique properties from polyhedron
			const auto &attributPropertyAperture =block_mesh.polyhedron_attribute_manager().find_or_create_attribute<geode::VariableAttribute, double>(aperture_property_name_,ndvalue_);

			aperture_[ifrac] = attributPropertyAperture->value(id_cell.polyhedronIndex);
		}
	}
	if (density_law_type_ > 0) {
		if (density_type_index_ > 1) {
			bool testPoperty =block_mesh.polyhedron_attribute_manager().attribute_exists(density_property_name_);

			if (!testPoperty) {
				cout << " The Property " << density_property_name_
						<< " doesn't exist for block " << model.block(id_cell.blockId).name() << endl;
			}else{
			// get petrophysique properties from polyhedron
			const auto &attributPropertyDensity =block_mesh.polyhedron_attribute_manager().find_or_create_attribute<geode::VariableAttribute, double>(density_property_name_,ndvalue_);

			density_val_[ifrac] = attributPropertyDensity->value(id_cell.polyhedronIndex);
			}
		}
	}

	if (distribution_shadow_zone_type_index_ > 1) {

			bool testPoperty =block_mesh.polyhedron_attribute_manager().attribute_exists(width_stress_zone_property_name_);

			if (!testPoperty) {
				cout << " The Property " << width_stress_zone_property_name_
						<< " doesn't exist for block " << model.block(id_cell.blockId).name() << endl;

			}else{
				// get petrophysique properties from polyhedron
				const auto &attributPropertyShadow_zone =block_mesh.polyhedron_attribute_manager().find_or_create_attribute<geode::VariableAttribute, double>(width_stress_zone_property_name_,ndvalue_);
				shadow_zone_val_[ifrac] = attributPropertyShadow_zone->value(id_cell.polyhedronIndex);
			}
	}


}
bool FractureSet::update_values_from_fracture_index(bool complementary,
		int &ifrac, double &length_max_l1, double &length_max_l2,boost::mt19937& center_gen) {
	//cout << randLength1_() << endl;
	if (distribution_length1_type_index_ == 0) {
		if(length1_type_ == distribution_type::power){
			length1_[ifrac] = randDistPowerLength1_.GcInvPower(center_gen);
		}else{
			 length1_[ifrac] = randLength1_();
		}
	} else {
		length1_[ifrac] = distr_length1_val_;
	}
	if (distribution_length2_type_index_ == 0) {
		if(length2_type_ == distribution_type::power){
			length2_[ifrac] = randDistPowerLength2_.GcInvPower(center_gen);
		}else{
			length2_[ifrac] = randLength2_();
		}
	} else {
		length2_[ifrac] = distr_length2_val_;
	}
	length_max_l1 = std::max((double) ((length_max_l1)), length1_[ifrac]);
	length_max_l2 = std::max((double) ((length_max_l2)), length2_[ifrac]);
	if (type_array_ == array_type::anastomosing) {
		if (ifrac & 1) {
			if (distribution_azimut_type_index_ == 0) {
				//Orientation_[ifrac] = randAzimutAnastomosing_();
				if(azimuth_type_ == distribution_type::power){
				  Orientation_[ifrac] = randDistPowerAzimuthAnastomosing_.GcInvPower(center_gen);
				}else{
				 Orientation_[ifrac] = randAzimutAnastomosing_();
				}
			} else {
				Orientation_[ifrac] = 270 - distr_azimut_val_;
			}
		} else {
			if (distribution_azimut_type_index_ == 0) {
				//Orientation_[ifrac] = randAzimut_();
				if(azimuth_type_ == distribution_type::power){
				  Orientation_[ifrac] = randDistPowerAzimuth_.GcInvPower(center_gen);
				}else{
				 Orientation_[ifrac] = randAzimut_();
				}
			} else {
				Orientation_[ifrac] = distr_azimut_val_;
			}
		}
	} else if (type_array_ == array_type::corridor) {
		Orientation_[ifrac] = randAzimut_();
	} else {
		if (distribution_azimut_type_index_ == 0) {
			//Orientation_[ifrac] = randAzimut_();
			if(azimuth_type_ == distribution_type::power){
			  Orientation_[ifrac] = randDistPowerAzimuth_.GcInvPower(center_gen);
			}else{
			 Orientation_[ifrac] = randAzimut_();
			}
		} else {
			Orientation_[ifrac] = distr_azimut_val_;
		}
	}

	// generate the complementary set
	if (complementary_fracture_set_active_ == true) {
		if (complementary == true) {
			Orientation_[ifrac] += 180.0;
			if (Orientation_[ifrac] > 360)
				Orientation_[ifrac] = Orientation_[ifrac] - 360;

			complementary = false;
		} else
			complementary = false;
	}
	if (distribution_angle_type_index_ == 0) {
		//Dip_[ifrac] = randDip_();
		if(dip_angle_type_ == distribution_type::power){
			Dip_[ifrac] = randDistPowerDipAngle_.GcInvPower(center_gen);
		}else{
		Dip_[ifrac] = randDip_();
		}
	} else {
		Dip_[ifrac] = distr_angle_val_;
	}
	if (distribution_aperture_type_index_ == 0) {
		//aperture_[ifrac] = randAperture_();
		if(aperture_type_ == distribution_type::power){
		 aperture_[ifrac] = randDistPowerAperture_.GcInvPower(center_gen);
		}else{
		 aperture_[ifrac] = randAperture_();
		}
	} else {
		aperture_[ifrac] = distr_aperture_val_;
	}
	if (density_law_type_ > 0) {
		if (density_type_index_ == 0) {
			//density_val_[ifrac] = randDensity_();
			if(density_type_ == distribution_type::power){
			 density_val_[ifrac] = randDistPowerDensity_.GcInvPower(center_gen);
			}else{
			 density_val_[ifrac] = randDensity_();
			}
		} else {
			density_val_[ifrac] = density_;
		}
	}
	if (density_hanging_wall_type_ > 0) {
		if (density_type_hanging_wall_index_ == 0) {
			//density_hanging_wall_val_[ifrac] = randHangingWallDensity_();
			if(density_hanging_wall_type_ == distribution_type::power){
				density_hanging_wall_val_[ifrac] = randDistPowerHangingWall_.GcInvPower(center_gen);
			}else{
				density_hanging_wall_val_[ifrac] = randHangingWallDensity_();
			}
		} else {
			density_hanging_wall_val_[ifrac] = density_hanging_wall_;
		}
	}
	if (density_foot_wall_type_ > 0) {
		if (density_type_foot_wall_index_ == 0) {
			//density_foot_wall_val_[ifrac] = randFootWallDensity_();
			if(density_foot_wall_type_ == distribution_type::power){
				density_foot_wall_val_[ifrac] = randDistPowerFootWall_.GcInvPower(center_gen);
			}else{
				density_foot_wall_val_[ifrac] = randFootWallDensity_();
			}
		} else {
			density_foot_wall_val_[ifrac] = density_foot_wall_;
		}
	}
	if (distribution_shadow_zone_type_index_ == 1) {
		//shadow_zone_val_[ifrac] = randShadowZone_();
		if(shadow_zone_type_ == distribution_type::power){
			shadow_zone_val_[ifrac] = randDistPowerShadowZone_.GcInvPower(center_gen);
		}else{
			shadow_zone_val_[ifrac] = randShadowZone_();
		}
	} else {
		shadow_zone_val_[ifrac] = default_stress_zone_max_;
	}
	return complementary;
}

double FractureSet::update_data_from_distribution(boost::mt19937& center_gen) {
	initialize_fracture_set_distribution(center_gen);
	// test if asociated tsurf object exist or not
	// if( fractures_ == NULL ) return;
	initialize_geometrical_fractures_data();

	double length_max_l1 = 0.0;
	double length_max_l2 = 0.0;

	//	SmartPointer<PercentageJobTracer> bjob = new PercentageJobTracer(Format::tr("Fractures Generation..."), 2*fractures_number_, true);
	int job_evolution = 0;
	bool complementary = false;
	//randLength1_ = boost::bind(boost::random::uniform_real_distribution<>(minmaxmode_values_[0], minmaxmode_values_[1]),center_gen);
	// determined the fractures parameters
	for (int ifrac = 0; ifrac < fractures_number_; ifrac++) {
		//cout << randLength1_() << endl;
		complementary = update_values_from_fracture_index(complementary, ifrac,
				length_max_l1, length_max_l2,center_gen);
		job_evolution++;
	}
  return std::max(length_max_l1,length_max_l2);
}





void FractureSet::fractures_evaluate_generation(const geode::StructuralModel &model, std::shared_ptr< geode::AABBTree3D >& tree,std::shared_ptr< BoxAABBEvalDistance3D >& disteval, std::vector<CellId> cells,std::shared_ptr<fractures_intersect::fault_array>& test_array,  bool estim_volume)
{
    // init the tsurf
    init_tsurf(name_);

	bool density_variable_in_grid = false;

	if (density_property_name_ != "none" && density_type_index_ !=0)
	{
		density_variable_in_grid = true;

	}

	boost::mt19937  center_gen;
	boost::uint32_t seed;
	seed=seed_*( array_index_ + fracture_set_index_+1);
    center_gen.seed(seed);

    float theta_max = evaluate_fracture_number(model,cells,density_variable_in_grid);
}

//-------------------------------------------
// generation of the fractures
//-------------------------------------------

void FractureSet::fractures_generation(const geode::StructuralModel &model, std::shared_ptr< geode::AABBTree3D >& tree,std::shared_ptr< BoxAABBEvalDistance3D >& disteval, std::vector<CellId> cells,std::shared_ptr<fractures_intersect::fault_array>& test_array,  bool estim_volume)
{
    // init the tsurf
    init_tsurf(name_);

	bool density_variable_in_grid = false;

	if (density_property_name_ != "none" && density_type_index_ !=0)
	{
		density_variable_in_grid = true;

	}

	boost::mt19937  center_gen;
	boost::uint32_t seed;
	seed=seed_*( array_index_ + fracture_set_index_+1);
    center_gen.seed(seed);

    float theta_max = evaluate_fracture_number(model,cells,density_variable_in_grid);

	// return if just a estimation of the fractures number
	if (estim_volume == true) return;

	// generation of number_frac fractures according to a Poisson
	// distribution law where ntot
	// is the Poisson parameter equal to mean number of fractures
	boost::uniform_01<boost::mt19937> center_uni(center_gen);

	if(density_law_type_>0) evaluate_final_fracture_number(center_uni);
	
	geode::Point3D center_p, p0_uv;
		
	double length_max=update_data_from_distribution(center_gen);

	float l_w= std::abs(generation_box_.max().value(2)-generation_box_.min().value(2));
	geode::Vector3D box_length(generation_box_.min(), generation_box_.max());
	int index = 0 ;

	float* densityTab = NULL;
	// case with density variable in a grid
	if (density_variable_in_grid == true) 
	{

		index = calculation_fractures_from_property(model,tree,disteval,cells,test_array,estim_volume,center_uni);



	}
	// other cases
	else
	{


		float z_init= generation_box_.min().value(2)-(length_max/2); // utiliser la bouding box de la grille :

		float z_bloc= (l_w + length_max);

		p0_uv.set_value(0, 0.5);
		p0_uv.set_value(1, 0.5);
		p0_uv.set_value(2,0);


		int security = 0; 


		// loop until all the fractures are generated 
		while ( index < fractures_number_ && security < fractures_number_*1000)
		{
			// random allocation of fracture center 
			// in real x,y,z 
			center_p.set_value(0, generation_box_.min().value(0)+center_uni()*box_length.value(0));
			center_p.set_value(1,generation_box_.min().value(1)+center_uni()*box_length.value(1));
			center_p.set_value(2,z_init + center_uni()* z_bloc);
			fractures_intersect::Point3D germe_pt = fractures_intersect::Point3D(center_p.value(0),center_p.value(1),center_p.value(2));

			if(!test_array->segment_point_array_center_intersect_fracture(germe_pt)) continue;
			// get the topo at center position
			float topo = generation_box_.max().value(2);
			float bottom =  generation_box_.min().value(2);

			double depth;

			// verify that the topo exist and the center is 
			// in the z generation domain
			if ( center_p.value(2) > topo )  depth = 0.0;
			else depth= topo - center_p.value(2);



			float theta_uni=theta_max*center_uni();

			float theta_depth= coef_a_* depth + density_;
			if (theta_depth <0.0) theta_depth=0.0;

			
			if ((center_p.value(2) > topo + 0.5) || (center_p.value(2) < bottom - 0.5) ) theta_depth=0.0;

			if  ( theta_uni < theta_depth)
			{
				geode::index_t box_id;
				geode::Point3D nearest_point;
				double distance;

				std::tie( box_id, nearest_point, distance ) = tree->closest_element_box( center_p, *disteval );

				CellId id_cell = cells[box_id];

				update_values_of_fracture_index_from_property(id_cell,model,index);


				if( add_fracture(center_p,p0_uv, length1_[index] , length2_[index], Orientation_[index],
    						 Dip_[index], aperture_[index], false, NULL,NULL) == true)
				{
					fractures_intersect::Point3D germe_pt = fractures_intersect::Point3D(center_p.value(0),center_p.value(1),center_p.value(2));
					std::shared_ptr<fractures_intersect::Cfaille> faille = std::shared_ptr<Cfaille>( new fractures_intersect::Cfaille(fracture_set_index_,germe_pt,Orientation_[index],Dip_[index],length1_[index] , length2_[index],aperture_[index],default_stress_zone_max_));
					//failles_vector_.push_back(faille);
					index++;


				}
				


			}

			security++; 
		
		} // end loop over the fractures


		//if (security == fractures_number_*1000) activated_security = true;


	} // end other cases 


	if(index == fractures_number_) {
		success_ = true;
		cout << " Number of created fractures for set " << name_ << " is : " << index << endl;

	}

	// update the phase and the activable fracture

	fractures_activable_ = fractures_number_;

	fractures_inserted_ = 0;




}

std::vector< fractures_intersect::Point3D  > FractureSet::get_echelon_position_from_initial_germ(fractures_intersect::Point3D pt, int ifrac){

	std::vector< fractures_intersect::Point3D  > vect_point_echelon;
	std::vector< fractures_intersect::Point3D  > vect_point_echelon_moins;
	std::vector< fractures_intersect::Point3D  > vect_point_echelon_plus;
	vect_point_echelon.clear();
	vect_point_echelon_plus.clear();
	vect_point_echelon_moins.clear();
	double space = echelon_space_val_[ifrac]/oriented_generation_box_.extension_.x();
	double x = pt.x();
	vect_point_echelon_plus.push_back(pt);
	while( x <= 1.0 + space){

		fractures_intersect::Point3D new_pt = fractures_intersect::Point3D(x+space,pt.y(),pt.z());
		vect_point_echelon_plus.push_back(new_pt);
		x +=space;
	}

	x = pt.x();
	while( x >= -space){

		fractures_intersect::Point3D new_pt = fractures_intersect::Point3D(x-space,pt.y(),pt.z());
		vect_point_echelon_moins.push_back(new_pt);
		x -=space;
	}


	int size=(vect_point_echelon_moins.size()-1);

	for(int i=size ; i >= 0; i--){
		vect_point_echelon.push_back(vect_point_echelon_moins[i]);
	}
	for(int j=0 ; j < vect_point_echelon_plus.size(); j++){
		vect_point_echelon.push_back(vect_point_echelon_plus[j]);
	}
	vect_point_echelon_plus.clear();
	vect_point_echelon_moins.clear();

	return vect_point_echelon;

}

void FractureSet::fractures_generation_from_shadow_zone_in_array(const geode::StructuralModel &model, std::shared_ptr< geode::AABBTree3D >& tree,
	    std::shared_ptr< BoxAABBEvalDistance3D >& disteval,  std::vector<CellId> cells,std::shared_ptr<fractures_intersect::fault_array>& test_array, bool estim_volume,fractures_intersect::Point3D pt_ref, int indice_array_loc)
{
    // init the tsurf
    init_tsurf(name_);

	bool density_variable_in_grid = false;

	if (density_property_name_ != "none" && density_type_index_ !=0)
	{
		density_variable_in_grid = true;

	}

	boost::mt19937  center_gen;
	boost::uint32_t seed;
	seed=seed_*( array_index_ + fracture_set_index_+1);
    center_gen.seed(seed);

    float theta_max = evaluate_fracture_number_from_oriented_box(model,cells,density_variable_in_grid);

	// return if just a estimation of the fractures number
	if (estim_volume == true) return;

	// generation of number_frac fractures according to a Poisson
	// distribution law where ntot
	// is the Poisson parameter equal to mean number of fractures
	boost::uniform_01<boost::mt19937> center_uni(center_gen);

	if(density_law_type_>0) evaluate_final_fracture_number(center_uni);

	geode::Point3D center_p, p0_uv;

	double length_max=update_data_from_distribution(center_gen);

	float l_w= std::abs(oriented_generation_box_.extension_.z());

	int index = 0 ;

	float* densityTab = NULL;
	std::vector<std::shared_ptr<fractures_intersect::Cfaille> > tab_fractures;
     if(is_shadow_zone_active_){
    	 type_intersection_= fractures_intersect::FracResIntersectionTypeEnum::CROSSING;
     }

     stress_zone_set_ = std::shared_ptr<fractures_intersect::fracture_set_geo >( new fractures_intersect::fracture_set_geo(fracture_set_index_, type_category_, seed_,tab_type_intersection_));

	// case with density variable in a grid
	if (density_variable_in_grid == true)
	{

		index = calculation_array_fractures_from_property(model,tree,disteval,cells,test_array,estim_volume,center_uni);



	}
	// other cases
	else
	{


		float z_init=  pt_ref.z() -(length_max/2); // utiliser la bouding box de la grille :

		float z_bloc= (l_w + length_max);

		p0_uv.set_value(0, 0.5);
		p0_uv.set_value(1, 0.5);
		p0_uv.set_value(2,0);


		int security = 0;

		center_p.set_value(0, center_uni());
		center_p.set_value(1, center_uni());
		center_p.set_value(2, center_uni());
		fractures_intersect::Point3D germe_pt ;

		std::vector< fractures_intersect::Point3D  >  vect_pt_echelon;


		if(type_array_== array_type::echelon){
			fractures_intersect::Point3D  local_germ(center_p.value(0),center_p.value(1),center_p.value(2));

			 vect_pt_echelon = get_echelon_position_from_initial_germ(local_germ,indice_array_loc);

			fractures_number_ = vect_pt_echelon.size();
			double length_max=update_data_from_distribution(center_gen);
		}



		// loop until all the fractures are generated
		while ( index < fractures_number_ && security < fractures_number_*100)
		{
			// random allocation of fracture center
			// in real x,y,z




			// get the topo at center position
			//float topo = pt_ref.z() + oriented_generation_box_.extension_.z();
			//float bottom =  pt_ref.z();

			double depth;

			// verify that the topo exist and the center is
			// in the z generation domain
			//if ( center_p.value(2) > topo )  depth = 0.0;
			//else depth= topo - center_p.value(2);



			float theta_uni=theta_max*center_uni();

			//float theta_depth= coef_a_* depth + density_;
			//if (theta_depth <0.0) theta_depth=0.0;


			//if ((center_p.value(2) > topo + length_max*0.5) || (center_p.value(2) < bottom - length_max*0.5) ) theta_depth=0.0;

			//if  ( theta_uni < theta_depth)
			//{
			    if(type_array_ != array_type::echelon){
			    	center_p.set_value(0, center_uni());
					center_p.set_value(1, center_uni());
					center_p.set_value(2, center_uni());
				    germe_pt = oriented_generation_box_.generate_point(center_p.value(0),center_p.value(1),center_p.value(2));
			    }else{
			    	if(index > vect_pt_echelon.size() ) continue;
			    	germe_pt = oriented_generation_box_.generate_point(vect_pt_echelon[index].coord_(0),vect_pt_echelon[index].coord_(1),vect_pt_echelon[index].coord_(2));
			    }
			    if(!test_array->segment_point_array_center_intersect_fracture(germe_pt)) continue;
				geode::index_t box_id;
				geode::Point3D nearest_point;
				double distance;
				geode::Point3D query({germe_pt.x(),germe_pt.y(),germe_pt.z()});
				std::tie( box_id, nearest_point, distance ) = tree->closest_element_box( query, *disteval );
				CellId id_cell = cells[box_id];

				update_values_of_fracture_index_from_property(id_cell,model,index);


				std::shared_ptr<fractures_intersect::Cfaille> faille = std::shared_ptr<fractures_intersect::Cfaille>(new fractures_intersect::Cfaille(fracture_set_index_,germe_pt,Orientation_[index],Dip_[index],length1_[index] , length2_[index],aperture_[index],shadow_zone_val_[index] > 0));


				if(type_array_==array_type::anastomosing){
					if (index & 1){

					 faille->familly_index_=2;
					}else{
					  faille->familly_index_=1;
					}
					stress_zone_set_->add_fracture_and_abbuting(faille);

					index++;
				}else{


					if(shadow_zone_val_[index] > 0){
										is_shadow_zone_active_ = true;
							if ( !stress_zone_set_->in_all_stress_zone__abutting_3d(faille,cylinder_step_))
							{
								//stress_zone_set_->add_fracture_and_create_cylinder(faille);

								index++;

							}


						}else{
							is_shadow_zone_active_ = false;
							stress_zone_set_->add_fracture(faille);

							index++;
						}

				}

			//}

			security++;

		} // end loop over the fractures



		//if (security == fractures_number_*1000) activated_security = true;


	} // end other cases


	int index_connection_fault =0;
	if(type_array_==array_type::relay){
		std::vector<int> global_number;
		std::vector<int> local_number(index);
		std::vector<std::pair<int,int>> correspondance;

		global_number.push_back(0);
		local_number[0]=1;
		int indice=-1;
		int jfaille=0;
		if(index > 1){
			do{

				double minDist=1E9;

				for(int ifaille=0; ifaille < index; ifaille++){
					if(ifaille == jfaille || local_number[ifaille] >0) continue;
					double  dist = stress_zone_set_->tab_fractures_[jfaille]->distance_to_plan(stress_zone_set_->tab_fractures_[ifaille]->Get_point_ref());
					if(dist < minDist  ){
						minDist=dist;
						indice=ifaille;
					}
				}
				global_number.push_back(indice);
				if(indice >= 0){
				std::pair<int,int> pairSet ={jfaille,indice};
				correspondance.push_back(pairSet);
				}
				jfaille=indice;
				local_number[indice]=1;


			}while(global_number.size() != index);

		}
		for(int jfaille=0; jfaille < index; jfaille++){

			for(int ifaille = jfaille+1; ifaille < index; ifaille++){

			 std::shared_ptr<fractures_intersect::Cfaille>  newFaille = getRelayConnectionFaille(fracture_set_index_,stress_zone_set_->tab_fractures_[global_number[jfaille]],stress_zone_set_->tab_fractures_[global_number[ifaille]] );

			 if(newFaille != NULL){
				 newFaille->connection_relay_=1;
				 stress_zone_set_->add_fracture_and_abbuting(newFaille);
				 index_connection_fault++;
			 }
			}

		}

	}
	nb_fractures_connection_relay_ = index_connection_fault;



	if(index == fractures_number_){

		success_ = true;
		cout << " Number of created fractures for set " << name_ << " is : " << index << endl;
	}
	// update the phase and the activable fracture

	fractures_activable_ = fractures_number_;

	fractures_inserted_ = 0;

	if(is_shadow_zone_active_){

		if(type_intersection_== fractures_intersect::FracResIntersectionTypeEnum::ABUTTING){
			cout<<" Fracture set "<< name_ << " have shadow zone with abbuting" <<endl;
		}else {
			cout<<" Fracture set "<< name_ << " have shadow zone with no abbuting" <<endl;
		}
	}




}

void FractureSet::fractures_generation_from_shadow_zone_in_corridor(const geode::StructuralModel &model, std::shared_ptr< geode::AABBTree3D >& tree,
	    std::shared_ptr< BoxAABBEvalDistance3D >& disteval, std::vector<CellId> cells, std::shared_ptr<fractures_intersect::fault_array>& test_array, bool estim_volume,fractures_intersect::Point3D pt_ref)
{
    // init the tsurf
    init_tsurf(name_);

	bool density_variable_in_grid = false;

	if (density_property_name_ != "none" && density_type_index_ !=0)
	{
		density_variable_in_grid = true;

	}

	boost::mt19937  center_gen;
	boost::uint32_t seed;
	seed=seed_*( array_index_ + fracture_set_index_+1);
    center_gen.seed(seed);

    float theta_max = evaluate_fracture_number_from_oriented_box(model,cells,density_variable_in_grid);

	// return if just a estimation of the fractures number
	if (estim_volume == true) return;

	// generation of number_frac fractures according to a Poisson
	// distribution law where ntot
	// is the Poisson parameter equal to mean number of fractures
	boost::uniform_01<boost::mt19937> center_uni(center_gen);

	if(density_law_type_>0) evaluate_final_fracture_number(center_uni);

	geode::Point3D center_p, p0_uv;

	double length_max=update_data_from_distribution(center_gen);

	float l_w= std::abs(oriented_generation_box_.extension_.z());

	int index = 0 ;

	float* densityTab = NULL;
	std::vector<std::shared_ptr<fractures_intersect::Cfaille> > tab_fractures;
     if(is_shadow_zone_active_){
    	 type_intersection_= fractures_intersect::FracResIntersectionTypeEnum::CROSSING;
     }

     stress_zone_set_ = std::shared_ptr<fractures_intersect::fracture_set_geo >( new fractures_intersect::fracture_set_geo(fracture_set_index_, type_category_, seed_,tab_type_intersection_));

	// case with density variable in a grid
	if (density_variable_in_grid == true)
	{



		index = calculation_array_fractures_from_property(model,tree,disteval,cells,test_array,estim_volume,center_uni);

	}
	// other cases
	else
	{


		float z_init=  pt_ref.z() -(length_max/2); // utiliser la bouding box de la grille :

		float z_bloc= (l_w + length_max);

		p0_uv.set_value(0, 0.5);
		p0_uv.set_value(1, 0.5);
		p0_uv.set_value(2,0);


		int security = 0;

		center_p.set_value(0, center_uni());
		center_p.set_value(1, center_uni());
		center_p.set_value(2, center_uni());
		fractures_intersect::Point3D germe_pt ;

		// loop until all the fractures are generated
		while ( index < fractures_number_ && security < fractures_number_*100)
		{
			// random allocation of fracture center
			// in real x,y,z




			// get the topo at center position
			//float topo = pt_ref.z() + oriented_generation_box_.extension_.z();
			//float bottom =  pt_ref.z();

			double depth;

			// verify that the topo exist and the center is
			// in the z generation domain
			//if ( center_p.value(2) > topo )  depth = 0.0;
			//else depth= topo - center_p.value(2);



			float theta_uni=theta_max*center_uni();

			//float theta_depth= coef_a_* depth + density_;
			//if (theta_depth <0.0) theta_depth=0.0;


			//if ((center_p.value(2) > topo + length_max*0.5) || (center_p.value(2) < bottom - length_max*0.5) ) theta_depth=0.0;

			//if  ( theta_uni < theta_depth)
			//{

				center_p.set_value(0, center_uni());
				center_p.set_value(1, center_uni());
				center_p.set_value(2, center_uni());


				geode::index_t box_id;
				geode::Point3D nearest_point;
				double distance;
			    std::tie( box_id, nearest_point, distance ) = tree->closest_element_box( center_p, *disteval );
			    germe_pt = oriented_generation_box_.generate_point(center_p.value(0),center_p.value(1),center_p.value(2));

			    if(!test_array->segment_point_array_center_intersect_fracture(germe_pt)) continue;
				CellId id_cell = cells[box_id];

				update_values_of_fracture_index_from_property(id_cell,model,index);

				std::shared_ptr<fractures_intersect::Cfaille> faille = std::shared_ptr<fractures_intersect::Cfaille>(new fractures_intersect::Cfaille(fracture_set_index_,germe_pt,Orientation_[index],Dip_[index],length1_[index] , length2_[index],aperture_[index],shadow_zone_val_[index]));



				if(shadow_zone_val_[index] > 0){
					is_shadow_zone_active_=true;
					if ( !stress_zone_set_->in_all_stress_zone__abutting_3d(faille,cylinder_step_))
					{
						//stress_zone_set_->add_fracture_and_create_cylinder(faille);

						index++;

					}


				}else{
					is_shadow_zone_active_=false;
					stress_zone_set_->add_fracture(faille);

					index++;
				}


			security++;

		} // end loop over the fractures



		//if (security == fractures_number_*1000) activated_security = true;


	} // end other cases


	int index_connection_fault =0;
	if(type_array_==array_type::relay){
		std::vector<int> global_number;
		std::vector<int> local_number(index);
		std::vector<std::pair<int,int>> correspondance;

		global_number.push_back(0);
		local_number[0]=1;
		int indice=-1;
		int jfaille=0;
		do{

			double minDist=1E9;

			for(int ifaille=0; ifaille < index; ifaille++){
				if(ifaille == jfaille || local_number[ifaille] >0) continue;
				double  dist = stress_zone_set_->tab_fractures_[jfaille]->distance_to_plan(stress_zone_set_->tab_fractures_[ifaille]->Get_point_ref());
				if(dist < minDist  ){
					minDist=dist;
					indice=ifaille;
				}
			}
			global_number.push_back(indice);
			if(indice >= 0){
			std::pair<int,int> pairSet ={jfaille,indice};
			correspondance.push_back(pairSet);
			}
			jfaille=indice;
			local_number[indice]=1;


		}while(global_number.size() != index);


		for(int jfaille=0; jfaille < index; jfaille++){

			for(int ifaille = jfaille+1; ifaille < index; ifaille++){

			 std::shared_ptr<fractures_intersect::Cfaille>  newFaille = getRelayConnectionFaille(fracture_set_index_,stress_zone_set_->tab_fractures_[global_number[jfaille]],stress_zone_set_->tab_fractures_[global_number[ifaille]] );

			 if(newFaille != NULL){
				 newFaille->connection_relay_=1;
				 stress_zone_set_->add_fracture_and_abbuting(newFaille);
				 index_connection_fault++;
			 }
			}

		}

	}
	nb_fractures_connection_relay_ = index_connection_fault;



	if(index == fractures_number_){

		success_ = true;
		cout << " Number of created fractures for set " << name_ << " is : " << index << endl;
	}
	// update the phase and the activable fracture

	fractures_activable_ = fractures_number_;

	fractures_inserted_ = 0;

	if(is_shadow_zone_active_){

		if(type_intersection_== fractures_intersect::FracResIntersectionTypeEnum::ABUTTING){
			cout<<" Fracture set "<< name_ << " have shadow zone with abbuting" <<endl;
		}else {
			cout<<" Fracture set "<< name_ << " have shadow zone with no abbuting" <<endl;
		}
	}




}



void FractureSet::fractures_generation_from_shadow_zone_in_damage_fault_zone(const geode::StructuralModel &model, std::shared_ptr< geode::AABBTree3D >& tree,
	    std::shared_ptr< BoxAABBEvalDistance3D >& disteval, std::vector<CellId> cells, std::shared_ptr<fractures_intersect::fault_array>& test_array,  bool estim_volume,double damage_hanging_wall_zone, double damage_foot_wall_zone)
{
    // init the tsurf
    init_tsurf(name_);

	bool density_variable_in_grid = false;

	if (density_property_name_ != "none" && density_type_index_ !=0)
	{
		density_variable_in_grid = true;

	}

	boost::mt19937  center_gen;
	boost::uint32_t seed;
	seed=seed_*( array_index_ + fracture_set_index_+1);
    center_gen.seed(seed);

    float theta_max_hangingwall = evaluate_fracture_number_from_oriented_box_hanging_wall_zone(model,cells,density_variable_in_grid, damage_hanging_wall_zone);
    float theta_max_footwall = evaluate_fracture_number_from_oriented_box_foot_wall_zone(model,cells,density_variable_in_grid, damage_foot_wall_zone);

	// return if just a estimation of the fractures number
	if (estim_volume == true) return;

	// generation of number_frac fractures according to a Poisson
	// distribution law where ntot
	// is the Poisson parameter equal to mean number of fractures
	boost::uniform_01<boost::mt19937> center_uni(center_gen);

	if(density_hanging_wall_law_type_>0) evaluate_final_fracture_number_hanging_wall(center_uni);
	if(density_foot_wall_law_type_>0) evaluate_final_fracture_number_foot_wall(center_uni);

	fractures_number_ = fractures_number_hanging_wall_ + fractures_number_foot_wall_;
	geode::Point3D center_p, p0_uv;

	double length_max=update_data_from_distribution(center_gen);

	float l_w_hanging= std::abs(oriented_generation_box_.extension_.z()*damage_hanging_wall_zone);
	float l_w_footwall= std::abs(oriented_generation_box_.extension_.z()*damage_foot_wall_zone);


	int index = 0 ;
	int index_foot_wall = 0 ;
	int index_hanging_wall = 0 ;

	float* densityTab = NULL;
	std::vector<std::shared_ptr<fractures_intersect::Cfaille> > tab_fractures;
     if(is_shadow_zone_active_){
    	 type_intersection_= fractures_intersect::FracResIntersectionTypeEnum::CROSSING;
     }

     stress_zone_set_ = std::shared_ptr<fractures_intersect::fracture_set_geo >( new fractures_intersect::fracture_set_geo(fracture_set_index_, type_category_, seed_,tab_type_intersection_));

	// case with density variable in a grid
	if (density_variable_in_grid == true)
	{

		index_foot_wall = calculation_foot_wall_fractures_from_property(model,tree,disteval,cells,test_array,estim_volume,center_uni,damage_foot_wall_zone);
		int index_hanging = calculation_hanging_wall_fractures_from_property(model,tree,disteval,cells,test_array,estim_volume,center_uni,damage_foot_wall_zone,damage_hanging_wall_zone);

		index_hanging_wall=index_foot_wall + index_hanging;

	}
	// other cases
	else
	{


		//float z_init= pt_ref.z() -(length_max/2); // utiliser la bouding box de la grille :

		float z_bloc_hanging= (l_w_hanging + length_max);

		float z_bloc_footwall= (l_w_footwall + length_max);

		p0_uv.set_value(0, 0.5);
		p0_uv.set_value(1, 0.5);
		p0_uv.set_value(2,0);


		int security = 0;

		// loop until all the fractures are generated
		while ( index_foot_wall < fractures_number_foot_wall_  && security < fractures_number_foot_wall_*100)
		{
			center_p.set_value(0, center_uni());
			center_p.set_value(1, center_uni());
			center_p.set_value(2, center_uni()*damage_foot_wall_zone);
			fractures_intersect::Point3D germe_pt = oriented_generation_box_.generate_point(center_p.value(0),center_p.value(1),center_p.value(2));
			if(!test_array->segment_point_array_center_intersect_fracture(germe_pt)) continue;
			geode::index_t box_id;
			geode::Point3D nearest_point;
			double distance;

			std::tie( box_id, nearest_point, distance ) = tree->closest_element_box( center_p, *disteval );


				if(center_p.value(2)/damage_foot_wall_zone >= center_uni() ){


					CellId id_cell = cells[box_id];

					update_values_of_fracture_index_from_property(id_cell,model,index);

					std::shared_ptr<fractures_intersect::Cfaille> faille = std::shared_ptr<fractures_intersect::Cfaille>(new fractures_intersect::Cfaille(fracture_set_index_,germe_pt,Orientation_[index_foot_wall],Dip_[index_foot_wall],length1_[index_foot_wall] , length2_[index_foot_wall],aperture_[index_foot_wall],shadow_zone_val_[index]));

					if(shadow_zone_val_[index] > 0){
						is_shadow_zone_active_=true;
						if ( !stress_zone_set_->in_all_stress_zone__abutting_3d(faille,cylinder_step_))
						{
							index_foot_wall++;
						}

					}else{
						is_shadow_zone_active_=false;
						stress_zone_set_->add_fracture(faille);
						index_foot_wall++;
					}
				}

			security++;

		}

		security = 0;

		index_hanging_wall +=index_foot_wall;
		// loop until all the fractures are generated
		while ( (index_hanging_wall-index_foot_wall) < fractures_number_hanging_wall_  && security < fractures_number_hanging_wall_*100)
		{
			center_p.set_value(0, center_uni());
			center_p.set_value(1, center_uni());
			center_p.set_value(2, damage_foot_wall_zone + center_uni()*damage_hanging_wall_zone);
			fractures_intersect::Point3D germe_pt = oriented_generation_box_.generate_point(center_p.value(0),center_p.value(1),center_p.value(2));
			if(!test_array->segment_point_array_center_intersect_fracture(germe_pt)) continue;
			geode::index_t box_id;
			geode::Point3D nearest_point;
			double distance;

			std::tie( box_id, nearest_point, distance ) = tree->closest_element_box( center_p, *disteval );


				if( center_uni() <= ( 1 - (center_p.value(2) - damage_foot_wall_zone)/damage_hanging_wall_zone)){



					CellId id_cell = cells[box_id];

					update_values_of_fracture_index_from_property(id_cell,model,index);

					std::shared_ptr<fractures_intersect::Cfaille> faille = std::shared_ptr<fractures_intersect::Cfaille>(new fractures_intersect::Cfaille(fracture_set_index_,germe_pt,Orientation_[index_hanging_wall],Dip_[index_hanging_wall],length1_[index_hanging_wall] , length2_[index_hanging_wall],aperture_[index_hanging_wall],shadow_zone_val_[index]));

				if(shadow_zone_val_[index] > 0){
					is_shadow_zone_active_=true;
					if ( !stress_zone_set_->in_all_stress_zone__abutting_3d(faille,cylinder_step_))
					{
						index_hanging_wall++;
					}

				}else{
					is_shadow_zone_active_=false;
					stress_zone_set_->add_fracture(faille);
					index_hanging_wall++;
				}



				}
			security++;

		 } // end loop over the fractures


	} // end other cases



	if(index_foot_wall == fractures_number_foot_wall_){

		success_ = true;
		cout << " Number of created foot wall fractures for set " << name_ << " is : " << index_foot_wall << endl;
	}
	// update the phase and the activable fracture

	fractures_activable_foot_wall_ = fractures_number_foot_wall_;


	if((index_hanging_wall-index_foot_wall) == fractures_number_hanging_wall_){

		success_ = true;
		cout << " Number of created hanging wall fractures for set " << name_ << " is : " << index_hanging_wall-index_foot_wall << endl;
	}
	// update the phase and the activable fracture

	fractures_activable_hanging_wall_ = fractures_number_hanging_wall_;

	fractures_number_ = fractures_number_foot_wall_ + fractures_number_hanging_wall_;
	fractures_activable_ = fractures_activable_foot_wall_ + fractures_activable_hanging_wall_;

	if(is_shadow_zone_active_){

		if(type_intersection_== fractures_intersect::FracResIntersectionTypeEnum::ABUTTING){
			cout<<" Fracture set "<< name_ << " have shadow zone with abbuting" <<endl;
		}else {
			cout<<" Fracture set "<< name_ << " have shadow zone with no abbuting" <<endl;
		}
	}




}



/**
 * @fn std::shared_ptr<fractures_intersect::Cfaille> getRelayConnectionFaille(int, std::shared_ptr<fractures_intersect::Cfaille>, std::shared_ptr<fractures_intersect::Cfaille>)
 * @brief
 *
 * @param index_fractureSet
 * @param faille1
 * @param faille2
 * @return
 */

std::shared_ptr<fractures_intersect::Cfaille> FractureSet::getRelayConnectionFaille( int index_fractureSet,std::shared_ptr<fractures_intersect::Cfaille> faille1, std::shared_ptr<fractures_intersect::Cfaille> faille2){



	fractures_intersect::Point3D germe_pt;

	germe_pt[0]= (faille1->Get_point_ref().x() + faille2->Get_point_ref().x())/2;
	germe_pt[1]= (faille1->Get_point_ref().y() + faille2->Get_point_ref().y())/2;
	germe_pt[2]= (faille1->Get_point_ref().z() + faille2->Get_point_ref().z())/2;

	double l2 = std::max(faille1->ext2L_selon_lpgp, faille2->ext2L_selon_lpgp);

	double l1 = distance(faille1->Get_point_ref(), faille2->Get_point_ref());

	double fl1= faille1->ext2L_selon_lpgp/2;
	double fl2= faille2->ext2L_selon_lpgp/2;
	double dist_min = std::min(faille1->ext2L_selon_lpgp/2, faille2->ext2L_selon_lpgp/2)*std::sqrt(2);


	bool inFaille1 = faille1->in_Plan(germe_pt);

	bool inFaille2 = faille2->in_Plan(germe_pt);



	if(!inFaille1 || !inFaille2 ) return NULL;


	double azimut = std::max(faille1->Get_azimut(), faille2->Get_azimut()) + 90.0;

	double dip = std::max(faille1->pendage, faille2->pendage);
	double aperture = (faille1->aperture_+ faille2->aperture_ )/2;

	std::shared_ptr<fractures_intersect::Cfaille> failleConnection = std::shared_ptr< fractures_intersect::Cfaille> (new fractures_intersect::Cfaille(index_fractureSet,germe_pt,azimut,dip,l1,l2,aperture,0));


	return failleConnection;


}
//-------------------------------------------
// generation of the fractures
//-------------------------------------------

void FractureSet::fractures_generation_from_shadow_zone(const geode::StructuralModel &model, std::shared_ptr< geode::AABBTree3D >& tree,
	    std::shared_ptr< BoxAABBEvalDistance3D >& disteval,  std::vector<CellId> cells, std::shared_ptr<fractures_intersect::fault_array>& test_array, bool estim_volume)
{
    // init the tsurf
    init_tsurf(name_);

	bool density_variable_in_grid = false;

	if (density_property_name_ != "none" && density_type_index_ !=0)
	{
		density_variable_in_grid = true;

	}

	boost::mt19937  center_gen;
	boost::uint32_t seed;
	seed=seed_*( array_index_ + fracture_set_index_+1);
    center_gen.seed(seed);

    float theta_max = evaluate_fracture_number(model,cells,density_variable_in_grid);

	// return if just a estimation of the fractures number
	if (estim_volume == true) return;

	// generation of number_frac fractures according to a Poisson
	// distribution law where ntot
	// is the Poisson parameter equal to mean number of fractures
	boost::uniform_01<boost::mt19937> center_uni(center_gen);

	if(density_law_type_>0) evaluate_final_fracture_number(center_uni);

	geode::Point3D center_p, p0_uv;

	double length_max=update_data_from_distribution(center_gen);

	float l_w= std::abs(generation_box_.max().value(2)-generation_box_.min().value(2));
	geode::Vector3D box_length(generation_box_.min(), generation_box_.max());
	int index = 0 ;

	float* densityTab = NULL;

     if(is_shadow_zone_active_){
    	 type_intersection_= fractures_intersect::FracResIntersectionTypeEnum::CROSSING;
     }

     type_intersection_ = tab_type_intersection_[type_category_];

     stress_zone_set_ = std::shared_ptr<fractures_intersect::fracture_set_geo >( new fractures_intersect::fracture_set_geo(fracture_set_index_, type_category_, seed_,tab_type_intersection_));


 	std::string fracname =  absl::StrCat(name_, "_", fracture_set_index_);
 	std::string surface_final = absl::StrCat(this->output_prefixe_name_,  fracname);
 	std::string filename = absl::StrCat(directory_path_,surface_final, ".dat");
    std::fstream fout(filename,fout.out);
	// case with density variable in a grid
	if (density_variable_in_grid == true)
	{

		index = calculation_fractures_from_property(model,tree,disteval,cells,test_array, estim_volume,center_uni);


	}
	// other cases
	else
	{

        double pas_dist = 0.5;
        if(length_max*pas_dist > l_w){
        	pas_dist = 0.03125;
        }
		float z_init= generation_box_.min().value(2)-(length_max*pas_dist); // utiliser la bouding box de la grille :

		float z_bloc= (l_w + length_max*2*pas_dist);

		p0_uv.set_value(0, 0.5);
		p0_uv.set_value(1, 0.5);
		p0_uv.set_value(2,0);


		int security = 0;

		geode::index_t box_id;
		geode::Point3D nearest_point;
		double distance;

		 double ndv = -9999.00;

		 if(fout.is_open()){
			 fout<< fractures_number_ << endl;
		 }
		// loop until all the fractures are generated
		while ( index < fractures_number_ && security < fractures_number_*100)
		{
			// random allocation of fracture center
			// in real x,y,z
			center_p.set_value(0, generation_box_.min().value(0)+center_uni()*box_length.value(0));
			center_p.set_value(1,generation_box_.min().value(1)+center_uni()*box_length.value(1));
			center_p.set_value(2,z_init + center_uni()* z_bloc);


			if(index ==21){
				int ducon=0;
				ducon++;
			}

			// get the topo at center position
			float topo = generation_box_.max().value(2);
			float bottom =  generation_box_.min().value(2);

			double depth;

			// verify that the topo exist and the center is
			// in the z generation domain
			if ( center_p.value(2) > topo )  depth = 0.0;
			else depth= topo - center_p.value(2);



			float theta_uni=theta_max*center_uni();

			float theta_depth= coef_a_* depth + density_;
			if (theta_depth <0.0) theta_depth=0.0;


			if ((center_p.value(2) > topo + length_max*pas_dist) || (center_p.value(2) < bottom - length_max*pas_dist) ) theta_depth=0.0;

			if((theta_uni < theta_depth &&  density_law_type_ > 0 ) || density_law_type_==0)
			{
				fractures_intersect::Point3D germe_pt = fractures_intersect::Point3D(center_p.value(0),center_p.value(1),center_p.value(2));
				if(test_array != NULL && !test_array->segment_point_array_center_intersect_fracture(germe_pt)) continue;

				geode::Point3D query({germe_pt.x(), germe_pt.y(),germe_pt.z()});
				std::tie( box_id, nearest_point, distance ) = tree->closest_element_box( query, *disteval );
				CellId id_cell = cells[box_id];

				update_values_of_fracture_index_from_property(id_cell,model,index);


				std::shared_ptr<fractures_intersect::Cfaille> faille = std::shared_ptr<fractures_intersect::Cfaille>(new fractures_intersect::Cfaille(fracture_set_index_,germe_pt,Orientation_[index],Dip_[index],length1_[index] , length2_[index],aperture_[index],shadow_zone_val_[index]));
				if(fout.is_open()){faille->save(fout);}
				if(shadow_zone_val_[index] > 0){
					is_shadow_zone_active_ = true;
					if ( !stress_zone_set_->in_all_stress_zone__abutting_3d(faille,cylinder_step_))
					{
						//stress_zone_set_->add_fracture_and_create_cylinder(faille);


						index++;


					}

				}else{
					is_shadow_zone_active_ = false;
					//if(index==1){
					stress_zone_set_->add_fracture(faille);
					//}

					index++;
				}

			}

			security++;

		} // end loop over the fractures


		//if (security == fractures_number_*1000) activated_security = true;


	} // end other cases

	fout.close();

	if(index == fractures_number_){
		success_ = true;
		cout << " Number of created fractures for set " << name_ << " is : " << index << endl;
	}
	// update the phase and the activable fracture

	fractures_activable_ = fractures_number_;

	fractures_inserted_ = 0;

	if(is_shadow_zone_active_){

		if(type_intersection_== fractures_intersect::FracResIntersectionTypeEnum::ABUTTING){
			cout<<" Fracture set "<< name_ << " have shadow zone with abbuting" <<endl;
		}else {
			cout<<" Fracture set "<< name_ << " have shadow zone with no abbuting" <<endl;
		}
	}




}



//----------------------------------------------
// create a TFace and add it in the Tsurf
// if no intersection with conditionnal info
//----------------------------------------------

bool FractureSet::add_fracture(geode::Point3D p0_xyz, geode::Point3D p0_uv , double l1 , double l2, double orientation,
    						   double dip, double aperture, bool check_intersect ,float* constraint_data, float* density_data)
{
	geode::Vector3D u,v;



	u.set_value(0, l1);
	u.set_value(1, 0.0);
	u.set_value(2, 0.0);
	v.set_value(0, 0.0);
	v.set_value(1, l2);
	v.set_value(2, 0.0);

		
	double theta, alpha;
	double rad=2*M_PI/360;

	alpha = M_PI- orientation*rad;
	theta = dip*rad;

	geode::Vector3D R,S;

	if (dip==0.0)
	{
		R.set_value(0,1.0);
		R.set_value(1,0.0);
		R.set_value(2,0.0);
		S.set_value(0,0.0);
		S.set_value(1,1.0);
		S.set_value(2,0.0);
	}
    else
	{
		R.set_value(0,cos(alpha));
		R.set_value(1,sin(alpha));
		R.set_value(2, 0.0);

		S.set_value(0,-cos(theta)*sin(alpha));
		S.set_value(1,cos(theta)*cos(alpha));
		S.set_value(2,sin(theta));
	}

	geode::Vector3D u1,v1;
	geode::Point3D p0,p1;

    p0.set_value(0,-l1 * p0_uv.value(0));
	p0.set_value(1,-l2 * p0_uv.value(1));
	p0.set_value(2,0.0);

	p1.set_value(0,R.value(0)*p0.value(0)+ S.value(0)*p0.value(1)+p0_xyz.value(0));
	p1.set_value(1,R.value(1)*p0.value(0)+ S.value(1)*p0.value(1)+p0_xyz.value(1));
	p1.set_value(2,R.value(2)*p0.value(0)+ S.value(2)*p0.value(1)+p0_xyz.value(2));
        
	u1.set_value(0,R.value(0)*u.value(0)+ S.value(0)*u.value(1));
	u1.set_value(1,R.value(1)*u.value(0)+ S.value(1)*u.value(1));
	u1.set_value(2,R.value(2)*u.value(0)+ S.value(2)*u.value(1));

	v1.set_value(0,R.value(0)*v.value(0)+ S.value(0)*v.value(1));
	v1.set_value(1,R.value(1)*v.value(0)+ S.value(1)*v.value(1));
	v1.set_value(2,R.value(2)*v.value(0)+ S.value(2)*v.value(1));





	vertices_.push_back(p1);

	vertices_.push_back(p1 + v1);

	vertices_.push_back(p1 + v1+ u1);

	vertices_.push_back(p1 + u1);


	std::array< geode::index_t,3 > tr1;
	tr1[0] = vertices_.size()-4;
	tr1[1] = vertices_.size()-2;
	tr1[2] = vertices_.size()-3;
	triangles_.push_back(tr1);
	triangles_l1_.push_back(l1);
	triangles_l2_.push_back(l2);
	triangles_orientation_.push_back(orientation);
	triangles_dip_.push_back(dip);
	triangles_aperture_.push_back(aperture);

	std::array< geode::index_t,3 > tr2;
	tr2[0] = vertices_.size()-4;
	tr2[1] = vertices_.size()-1;
	tr2[2] = vertices_.size()-2;
	triangles_.push_back(tr2);

	triangles_l1_.push_back(l1);
	triangles_l2_.push_back(l2);
	triangles_orientation_.push_back(orientation);
	triangles_dip_.push_back(dip);
	triangles_aperture_.push_back(aperture);


	return true; 
}


//-------------------------------------------------
// get the TSurf Object representing the fractures
//-------------------------------------------------


//----------------------------------------------
// get the number of fractures
//----------------------------------------------

int FractureSet::get_fractures_number()
{
	return fractures_number_;
}

//----------------------------------------------
// set the number of fractures
//----------------------------------------------

void FractureSet::set_fractures_number(int n_frac)
{
	
	fractures_activable_ += n_frac - fractures_number_;
	if (fractures_activable_ < 0) fractures_activable_ = 0;

	fractures_number_ = n_frac;
	evaluate_fractures_number_ = n_frac;
}


//----------------------------------------------
// set the number of fractures
//----------------------------------------------

void FractureSet::set_fractures_number_hanging_wall(int n_frac)
{

	fractures_activable_hanging_wall_ += n_frac - fractures_number_hanging_wall_;
	if (fractures_activable_hanging_wall_ < 0) fractures_activable_hanging_wall_ = 0;

	fractures_number_hanging_wall_ = n_frac;
	evaluate_fractures_number_hanging_wall_ = n_frac;
}

//----------------------------------------------
// set the number of fractures
//----------------------------------------------

void FractureSet::set_fractures_number_foot_wall(int n_frac)
{

	fractures_activable_foot_wall_ += n_frac - fractures_number_foot_wall_;
	if (fractures_activable_foot_wall_ < 0) fractures_activable_foot_wall_ = 0;

	fractures_number_foot_wall_ = n_frac;
	evaluate_fractures_number_foot_wall_ = n_frac;
}
//----------------------------------------------
// get the evaluate number of fractures
//----------------------------------------------

int FractureSet::get_evaluate_fractures_number()
{
	return evaluate_fractures_number_;
}

//----------------------------------------------
// get the number of activable fractures
//----------------------------------------------

int FractureSet::get_fractures_activables()
{
	return fractures_activable_;
}


//----------------------------------------------
// desactivate a given number of fractures
//----------------------------------------------

int FractureSet::desactivate_fractures(int n_fracture)
{
	int n = fractures_activable_;

	fractures_activable_ = fractures_activable_ + n_fracture;

	if (fractures_activable_ > fractures_number_ ) 
	{
		fractures_activable_= fractures_number_;
		return fractures_number_ - n;
	}
	else 
	{
		return n_fracture;
	}

}

//----------------------------------------------
// activate a given number of fractures
//----------------------------------------------

int FractureSet::activate_fractures(int n_fracture)
{
	int n = fractures_activable_;

	fractures_activable_ = fractures_activable_ - n_fracture;

	if (fractures_activable_ > 0) 
	{
		return n_fracture;
	}
	else 
	{
		fractures_activable_= 0;
		return n;
	}

}
//---------------------------------------------------
// estimate the fractures density for the karstsim sgrid
//---------------------------------------------------

void FractureSet::initialize_volume_from_global_grid_of_structural_model( const geode::StructuralModel &model)
{
	// get the total grid from geode structural model
	const auto& solidTot = geode::convert_brep_into_solid( model);
    const auto att_value =solidTot->polyhedron_attribute_manager().find_attribute<geode::uuid_from_conversion_attribute_type>(geode::uuid_from_conversion_attribute_name);
    const auto& box = solidTot->bounding_box();
    double lengthX = box.max().value(0) - box.min().value(0);
    double lengthY = box.max().value(1) - box.min().value(1);
    double lengthZ = box.max().value(2) - box.min().value(2);

    area_ = lengthX*lengthY;
    volume_=area_*lengthZ;

}


int FractureSet::calculation_fractures_from_property(const geode::StructuralModel &model, std::shared_ptr< geode::AABBTree3D >& tree,
		std::shared_ptr< BoxAABBEvalDistance3D >& disteval, std::vector<CellId> cells,std::shared_ptr<fractures_intersect::fault_array>& test_array, bool estim_volume, boost::uniform_01<boost::mt19937>& center_uni){


	geode::Point3D p0_uv, center_p;
	p0_uv.set_value(0, 0.5);
	p0_uv.set_value(1, 0.5);
	p0_uv.set_value(2,0);

	int security = 0;
	int index =0;
	geode::Vector3D box_length(generation_box_.min(), generation_box_.max());

	// loop until all the fractures are generated
	while ( index < fractures_number_ && security < fractures_number_*100)
	{

		center_p.set_value(0, generation_box_.min().value(0)+center_uni()*box_length.value(0));
		center_p.set_value(1, generation_box_.min().value(1)+center_uni()*box_length.value(1));
		center_p.set_value(2, generation_box_.min().value(2)+center_uni()*box_length.value(2));


			geode::index_t box_id;
			geode::Point3D nearest_point;
			double distance;
			fractures_intersect::Point3D germe_pt = oriented_generation_box_.generate_point(center_p.value(0),center_p.value(1),center_p.value(2));

			if(test_array != NULL && !test_array->segment_point_array_center_intersect_fracture(germe_pt)) continue;
			geode::Point3D query({center_p.value(0), center_p.value(1),center_p.value(2)});
			std::tie( box_id, nearest_point, distance ) = tree->closest_element_box( query, *disteval );


			// karstsim grid index
			// find nearest cell in grid to get cell density

			float center_density=cells_density_[box_id];


			if ((max_density_*center_uni() < center_density) )
			{



				CellId id_cell = cells[box_id];

				update_values_of_fracture_index_from_property(id_cell,model,index);

				fractures_intersect::Point3D germe_pt = fractures_intersect::Point3D(center_p.value(0),center_p.value(1),center_p.value(2));
				std::shared_ptr<fractures_intersect::Cfaille> faille = std::shared_ptr<fractures_intersect::Cfaille>(new fractures_intersect::Cfaille(fracture_set_index_,germe_pt,Orientation_[index],Dip_[index],length1_[index] , length2_[index],aperture_[index],shadow_zone_val_[index]));

				if(shadow_zone_val_[index] > 0){
					is_shadow_zone_active_ = true;
					if ( !stress_zone_set_->in_all_stress_zone__abutting_3d(faille,cylinder_step_))
					{
						stress_zone_set_->add_fracture_and_create_cylinder(faille);

						index++;

					}

				}else{
					is_shadow_zone_active_ = false;
					stress_zone_set_->add_fracture(faille);

					index++;
				}




			}

			security++;



	}

return index;


}


int FractureSet::calculation_array_fractures_from_property(const geode::StructuralModel &model, std::shared_ptr< geode::AABBTree3D >& tree,
		std::shared_ptr< BoxAABBEvalDistance3D >& disteval, std::vector<CellId> cells,std::shared_ptr<fractures_intersect::fault_array>& test_array, bool estim_volume, boost::uniform_01<boost::mt19937>& center_uni){


	geode::Point3D p0_uv, center_p;
	p0_uv.set_value(0, 0.5);
	p0_uv.set_value(1, 0.5);
	p0_uv.set_value(2,0);

	int security = 0;
	int index =0;
	geode::Vector3D box_length(generation_box_.min(), generation_box_.max());

	// loop until all the fractures are generated
	while ( index < fractures_number_ && security < fractures_number_*100)
	{

		center_p.set_value(0, center_uni());
		center_p.set_value(1, center_uni());
		center_p.set_value(2, center_uni());


			geode::index_t box_id;
			geode::Point3D nearest_point;
			double distance;
			fractures_intersect::Point3D germe_pt = oriented_generation_box_.generate_point(center_p.value(0),center_p.value(1),center_p.value(2));
			if(!test_array->segment_point_array_center_intersect_fracture(germe_pt)) continue;
			geode::Point3D query({germe_pt.x(), germe_pt.y(),germe_pt.z()});
			std::tie( box_id, nearest_point, distance ) = tree->closest_element_box( query, *disteval );

			// karstsim grid index
			// find nearest cell in grid to get cell density

			float center_density=cells_density_[box_id];


			if ((max_density_*center_uni() < center_density) )
			{

				CellId id_cell = cells[box_id];

				update_values_of_fracture_index_from_property(id_cell,model,index);

				std::shared_ptr<fractures_intersect::Cfaille> faille = std::shared_ptr<fractures_intersect::Cfaille>(new fractures_intersect::Cfaille(fracture_set_index_,germe_pt,Orientation_[index],Dip_[index],length1_[index] , length2_[index],aperture_[index],shadow_zone_val_[index]));

				if(shadow_zone_val_[index] > 0){
					is_shadow_zone_active_ = true;
					if ( !stress_zone_set_->in_all_stress_zone__abutting_3d(faille,cylinder_step_))
					{
						stress_zone_set_->add_fracture_and_create_cylinder(faille);

						index++;

					}

				}else{
					is_shadow_zone_active_ = false;
					stress_zone_set_->add_fracture(faille);

					index++;
				}




			}

			security++;



	}

return index;


}
int FractureSet::calculation_hanging_wall_fractures_from_property(const geode::StructuralModel &model, std::shared_ptr< geode::AABBTree3D >& tree,
		std::shared_ptr< BoxAABBEvalDistance3D >& disteval, std::vector<CellId> cells,std::shared_ptr<fractures_intersect::fault_array>& test_array, bool estim_volume, boost::uniform_01<boost::mt19937>& center_uni, double damage_foot_wall_zone,double damage_hanging_wall_zone){


	geode::Point3D p0_uv, center_p;
	p0_uv.set_value(0, 0.5);
	p0_uv.set_value(1, 0.5);
	p0_uv.set_value(2,0);

	int security = 0;
	int index =0;
	geode::Vector3D box_length(generation_box_.min(), generation_box_.max());

	// loop until all the fractures are generated
	while ( index < fractures_number_hanging_wall_  && security < fractures_number_hanging_wall_*100)
	{

		center_p.set_value(0, center_uni());
		center_p.set_value(1, center_uni());
		center_p.set_value(2, damage_foot_wall_zone + center_uni()*damage_hanging_wall_zone);


		if( center_uni() <= ( 1 - (center_p.value(2) - damage_foot_wall_zone)/damage_hanging_wall_zone)){

			fractures_intersect::Point3D germe_pt = oriented_generation_box_.generate_point(center_p.value(0),center_p.value(1),center_p.value(2));
			if(!test_array->segment_point_array_center_intersect_fracture(germe_pt)) continue;
			geode::index_t box_id;
			geode::Point3D nearest_point;
			double distance;
			geode::Point3D query({germe_pt.x(),germe_pt.y(),germe_pt.z()});
			std::tie( box_id, nearest_point, distance ) = tree->closest_element_box( query, *disteval );



			// karstsim grid index
			// find nearest cell in grid to get cell density

			float center_density=cells_density_[box_id];


			if ((max_density_*center_uni() < center_density) )
			{




				CellId id_cell = cells[box_id];

				update_values_of_fracture_index_from_property(id_cell,model,index);


				std::shared_ptr<fractures_intersect::Cfaille> faille = std::shared_ptr<fractures_intersect::Cfaille>(new fractures_intersect::Cfaille(fracture_set_index_,germe_pt,Orientation_[index],Dip_[index],length1_[index] , length2_[index],aperture_[index],shadow_zone_val_[index]));

				if(shadow_zone_val_[index] > 0){
					is_shadow_zone_active_ = true;
					if ( !stress_zone_set_->in_all_stress_zone__abutting_3d(faille,cylinder_step_))
					{
						stress_zone_set_->add_fracture_and_create_cylinder(faille);

						index++;

					}

				}else{
					is_shadow_zone_active_ = false;
					stress_zone_set_->add_fracture(faille);

					index++;
				}




			}

			security++;
		}


	}

	return index;

}

//
int FractureSet::calculation_foot_wall_fractures_from_property(const geode::StructuralModel &model, std::shared_ptr< geode::AABBTree3D >& tree,
		std::shared_ptr< BoxAABBEvalDistance3D >& disteval, std::vector<CellId> cells,std::shared_ptr<fractures_intersect::fault_array>& test_array, bool estim_volume, boost::uniform_01<boost::mt19937>& center_uni, double damage_foot_wall_zone){


	geode::Point3D p0_uv, center_p;
	p0_uv.set_value(0, 0.5);
	p0_uv.set_value(1, 0.5);
	p0_uv.set_value(2,0);

	int security = 0;
	int index =0;
	geode::Vector3D box_length(generation_box_.min(), generation_box_.max());

	// loop until all the fractures are generated
	while ( index < fractures_number_foot_wall_  && security < fractures_number_foot_wall_*100)
	{

		center_p.set_value(0, center_uni());
		center_p.set_value(1, center_uni());
		center_p.set_value(2, center_uni()*damage_foot_wall_zone);


		if(center_p.value(2)/damage_foot_wall_zone >= center_uni() ){

			geode::index_t box_id;
			geode::Point3D nearest_point;
			double distance;
			fractures_intersect::Point3D germe_pt = oriented_generation_box_.generate_point(center_p.value(0),center_p.value(1),center_p.value(2));
			if(!test_array->segment_point_array_center_intersect_fracture(germe_pt)) continue;
			geode::Point3D query({germe_pt.x(),germe_pt.y(),germe_pt.z()});
			std::tie( box_id, nearest_point, distance ) = tree->closest_element_box( query, *disteval );

			// karstsim grid index
			// find nearest cell in grid to get cell density

			float center_density=cells_density_[box_id];


			if ((max_density_*center_uni() < center_density) )
			{




				CellId id_cell = cells[box_id];

				update_values_of_fracture_index_from_property(id_cell,model,index);


				std::shared_ptr<fractures_intersect::Cfaille> faille = std::shared_ptr<fractures_intersect::Cfaille>(new fractures_intersect::Cfaille(fracture_set_index_,germe_pt,Orientation_[index],Dip_[index],length1_[index] , length2_[index],aperture_[index],shadow_zone_val_[index]));

				if(shadow_zone_val_[index] > 0){
					is_shadow_zone_active_ = true;
					if ( !stress_zone_set_->in_all_stress_zone__abutting_3d(faille,cylinder_step_))
					{
						stress_zone_set_->add_fracture_and_create_cylinder(faille);

						index++;

					}

				}else{
					is_shadow_zone_active_ = false;
					stress_zone_set_->add_fracture(faille);

					index++;
				}




			}

			security++;
		}


	}

  return index;
}
//---------------------------------------------------
// estimate the fractures density for the karstsim sgrid
//---------------------------------------------------

void FractureSet::initialize_density_P31_from_structural_model( const geode::StructuralModel &model,std::vector<CellId>& cells)
{
    geode_BlocName_to_BlockId_.clear();
    geode_BlocName_to_Density_exist_.clear();
    geode_BlocName_to_volume_.clear();
    geode_BlocName_to_area_.clear();
    geode_BlocName_to_density_.clear();
    double density_loc =0;
    double volume_loc=0;
    double density_tot=0;
    bool testDensity =true;

    for (const auto &block : model.blocks()){
    	const auto name = block.name();
    	const auto& mesh = block.mesh<geode::TetrahedralSolid3D>();
    	bool testDensity = mesh.polyhedron_attribute_manager().attribute_exists(density_property_name_);

		if(!testDensity){
				cout << " The Property " << density_property_name_<< " doesn't exist for block " << block.name() << endl;
		 }

    }
    cells_density_.resize(cells.size(), 0.0);
    for( const auto& cell_id : cells)
    {
		const auto name = model.block(cell_id.blockId).name();
		const auto& id = cell_id.blockId;
		geode_BlocName_to_BlockId_.insert({(std::string)name,cell_id.blockId});

		const auto& mesh = model.block(cell_id.blockId).mesh<geode::TetrahedralSolid3D>();
		bool testDensity = mesh.polyhedron_attribute_manager().attribute_exists(density_property_name_);
		volume_loc = geode::tetrahedron_volume(mesh.tetrahedron(cell_id.polyhedronIndex));
		if(!testDensity || volume_loc == 0){

			cells_density_[cell_id.polyhedronIndex] = 0.0;
		}else{

			//double lenghtX = mesh.bounding_box().max().value(0) - mesh.bounding_box().min().value(0);
			//double lenghtY = mesh.bounding_box().max().value(1) - mesh.bounding_box().min().value(1);
			//double area = lenghtX*lenghtY;
			//double lenghtZ = mesh.bounding_box().max().value(2) - mesh.bounding_box().min().value(2);
			//geode_BlocName_to_area_.insert({(std::string)name,area});
			//geode_BlocName_to_volume_.insert({(std::string)name,area*lenghtZ});

			//geode_BlocName_to_Density_exist_.insert({(std::string)name,testDensity});

			const auto& attributDensity_ = mesh.polyhedron_attribute_manager().find_or_create_attribute< geode::VariableAttribute, double >(density_property_name_, ndvalue_ );

			density_loc = volume_loc*attributDensity_->value(cell_id.polyhedronIndex);
;
			cells_density_[cell_id.polyhedronIndex] = density_loc ;
			double valeur = attributDensity_->value(cell_id.polyhedronIndex);
			//if(valeur > 0) std::cout << " density = " << valeur << " volume= " << volume_loc << " density= " << density_loc << endl;
			max_density_ = std::max(max_density_,density_loc);

			//geode_BlocName_to_density_.insert({(std::string)name,density_loc});
			density_+=density_loc;

		}


		volume_ += volume_loc;

    }
   if(volume_ > 0)  density_ /= volume_;


}
//---------------------------------------------------
// estimate the fractures density for the karstsim sgrid
//---------------------------------------------------

void FractureSet::initialize_density_P32_from_structural_model( const geode::StructuralModel &model, std::vector<CellId>& cells, double lfrac_mean)
{
    geode_BlocName_to_BlockId_.clear();
    geode_BlocName_to_Density_exist_.clear();
    geode_BlocName_to_volume_.clear();
    geode_BlocName_to_area_.clear();
    geode_BlocName_to_density_.clear();
    double density_loc =0;
    double volume_loc=0;
    double density_tot=0;

    for (const auto &block : model.blocks()){
    	const auto name = block.name();
    	const auto& mesh = block.mesh<geode::TetrahedralSolid3D>();
    	bool testDensity = mesh.polyhedron_attribute_manager().attribute_exists(density_property_name_);

		if(!testDensity){
				cout << " The Property " << density_property_name_<< " doesn't exist for block " << block.name() << endl;
		 }

    }
    cells_density_.resize(cells.size(), 0.0);
    for( const auto& cell_id : cells)
    {
		const auto name = model.block(cell_id.blockId).name();
		const auto& id = cell_id.blockId;
		geode_BlocName_to_BlockId_.insert({(std::string)name,cell_id.blockId});

		const auto& mesh = model.block(cell_id.blockId).mesh<geode::TetrahedralSolid3D>();
		bool testDensity = mesh.polyhedron_attribute_manager().attribute_exists(density_property_name_);

		volume_loc = geode::tetrahedron_volume(mesh.tetrahedron(cell_id.polyhedronIndex));
		if(!testDensity || volume_loc == 0){
			cells_density_[cell_id.polyhedronIndex] = 0.0;
		}else{

			//double lenghtX = mesh.bounding_box().max().value(0) - mesh.bounding_box().min().value(0);
			//double lenghtY = mesh.bounding_box().max().value(1) - mesh.bounding_box().min().value(1);
			//double area = lenghtX*lenghtY;
			//double lenghtZ = mesh.bounding_box().max().value(2) - mesh.bounding_box().min().value(2);
			//geode_BlocName_to_area_.insert({(std::string)name,area});
			//geode_BlocName_to_volume_.insert({(std::string)name,area*lenghtZ});

			//geode_BlocName_to_Density_exist_.insert({(std::string)name,testDensity});

			const auto& attributDensity_ = mesh.polyhedron_attribute_manager().find_or_create_attribute< geode::VariableAttribute, double >(density_property_name_, ndvalue_ );

			density_loc = volume_loc*attributDensity_->value(cell_id.polyhedronIndex)/lfrac_mean;
			cells_density_[cell_id.polyhedronIndex] =density_loc;
			double valeur = attributDensity_->value(cell_id.polyhedronIndex);
			//if(valeur > 0) std::cout << " density = " << valeur << " volume= " << volume_loc << " density= " << density_loc << endl;
			max_density_ = std::max(max_density_,density_loc);

			//geode_BlocName_to_density_.insert({(std::string)name,density_loc});
			density_+=density_loc;

		}


		volume_ += volume_loc;

    }
    if(volume_ > 0) {
    	density_ /= volume_;
    	density_ *=lfrac_mean;
    }
}
//---------------------------------------------------
// estimate the fractures density for the karstsim sgrid
//---------------------------------------------------

void FractureSet::initialize_density_P31_from_structural_model_in_array( const geode::StructuralModel &model,std::vector<CellId>& cells)
{
    geode_BlocName_to_BlockId_.clear();
    geode_BlocName_to_Density_exist_.clear();
    geode_BlocName_to_volume_.clear();
    geode_BlocName_to_area_.clear();
    geode_BlocName_to_density_.clear();
    double density_loc =0;
    double volume_loc=0;
    double density_tot=0;
    for (const auto &block : model.blocks()){
    	const auto name = block.name();
    	const auto& mesh = block.mesh<geode::TetrahedralSolid3D>();
    	bool testDensity = mesh.polyhedron_attribute_manager().attribute_exists(density_property_name_);

		if(!testDensity){
				cout << " The Property " << density_property_name_<< " doesn't exist for block " << block.name() << endl;
		 }

    }
    cells_density_.resize(cells.size(), 0.0);
    for( const auto& cell_id : cells)
    {
		const auto name = model.block(cell_id.blockId).name();
		const auto& id = cell_id.blockId;
		geode_BlocName_to_BlockId_.insert({(std::string)name,cell_id.blockId});

		const auto& mesh = model.block(cell_id.blockId).mesh<geode::TetrahedralSolid3D>();

		geode::Point3D barycenter = mesh.polyhedron_barycenter(cell_id.polyhedronIndex);
		fractures_intersect::Point3D pt(barycenter.value(0),barycenter.value(1),barycenter.value(2));
		if(!oriented_generation_box_.contains(pt)) {
			cells_density_[cell_id.polyhedronIndex] = 0.0;
			continue;
		}
		bool testDensity = mesh.polyhedron_attribute_manager().attribute_exists(density_property_name_);
		volume_loc = geode::tetrahedron_volume(mesh.tetrahedron(cell_id.polyhedronIndex));
		if(!testDensity || volume_loc==0){
			cells_density_[cell_id.polyhedronIndex] = 0.0;
		}else{

			//double lenghtX = mesh.bounding_box().max().value(0) - mesh.bounding_box().min().value(0);
			//double lenghtY = mesh.bounding_box().max().value(1) - mesh.bounding_box().min().value(1);
			//double area = lenghtX*lenghtY;
			//double lenghtZ = mesh.bounding_box().max().value(2) - mesh.bounding_box().min().value(2);
			//geode_BlocName_to_area_.insert({(std::string)name,area});
			//geode_BlocName_to_volume_.insert({(std::string)name,area*lenghtZ});

			//geode_BlocName_to_Density_exist_.insert({(std::string)name,testDensity});

			const auto& attributDensity_ = mesh.polyhedron_attribute_manager().find_or_create_attribute< geode::VariableAttribute, double >(density_property_name_, ndvalue_ );

			density_loc = volume_loc*attributDensity_->value(cell_id.polyhedronIndex);

			//cells_density_.push_back(density_loc);
			cells_density_[cell_id.polyhedronIndex] = density_loc;
			double valeur = attributDensity_->value(cell_id.polyhedronIndex);
			//if(valeur > 0) std::cout << " density = " << valeur << " volume= " << volume_loc << " density= " << density_loc << endl;
			max_density_ = std::max(max_density_,density_loc);

			//geode_BlocName_to_density_.insert({(std::string)name,density_loc});
			density_+=density_loc;

		}


		volume_ += volume_loc;

    }
    if(volume_ > 0)  density_ /= volume_;
}
//---------------------------------------------------
// estimate the fractures density for the karstsim sgrid
//---------------------------------------------------

void FractureSet::initialize_density_P32_from_structural_model_in_array( const geode::StructuralModel &model, std::vector<CellId>& cells, double lfrac_mean)
{
    geode_BlocName_to_BlockId_.clear();
    geode_BlocName_to_Density_exist_.clear();
    geode_BlocName_to_volume_.clear();
    geode_BlocName_to_area_.clear();
    geode_BlocName_to_density_.clear();
    double density_loc =0;
    double volume_loc=0;
    double density_tot=0;
    cells_density_.clear();
    for (const auto &block : model.blocks()){
    	const auto name = block.name();
    	const auto& mesh = block.mesh<geode::TetrahedralSolid3D>();
    	bool testDensity = mesh.polyhedron_attribute_manager().attribute_exists(density_property_name_);


		if(!testDensity){
				cout << " The Property " << density_property_name_<< " doesn't exist for block " << block.name() << endl;
		 }

    }
    cells_density_.resize(cells.size(), 0.0);

    for( const auto& cell_id : cells)
    {
		const auto name = model.block(cell_id.blockId).name();
		const auto& id = cell_id.blockId;
		geode_BlocName_to_BlockId_.insert({(std::string)name,cell_id.blockId});

		const auto& mesh = model.block(cell_id.blockId).mesh<geode::TetrahedralSolid3D>();
		geode::Point3D barycenter = mesh.polyhedron_barycenter(cell_id.polyhedronIndex);
		fractures_intersect::Point3D pt(barycenter.value(0),barycenter.value(1),barycenter.value(2));
		if(!oriented_generation_box_.contains(pt)) {
			cells_density_[cell_id.polyhedronIndex]=0.0;
			continue;
		}

		bool testDensity = mesh.polyhedron_attribute_manager().attribute_exists(density_property_name_);
		volume_loc = geode::tetrahedron_volume(mesh.tetrahedron(cell_id.polyhedronIndex));
		if(!testDensity || volume_loc == 0){
			cells_density_[cell_id.polyhedronIndex]=0.0;

		}else{

			//double lenghtX = mesh.bounding_box().max().value(0) - mesh.bounding_box().min().value(0);
			//double lenghtY = mesh.bounding_box().max().value(1) - mesh.bounding_box().min().value(1);
			//double area = lenghtX*lenghtY;
			//double lenghtZ = mesh.bounding_box().max().value(2) - mesh.bounding_box().min().value(2);
			//geode_BlocName_to_area_.insert({(std::string)name,area});
			//geode_BlocName_to_volume_.insert({(std::string)name,area*lenghtZ});

			//geode_BlocName_to_Density_exist_.insert({(std::string)name,testDensity});

			const auto& attributDensity_ = mesh.polyhedron_attribute_manager().find_or_create_attribute< geode::VariableAttribute, double >(density_property_name_, ndvalue_ );

			density_loc = volume_loc*attributDensity_->value(cell_id.polyhedronIndex)/lfrac_mean;
			//cells_density_.push_back(density_loc);
			cells_density_[cell_id.polyhedronIndex] = density_loc;
			double valeur = attributDensity_->value(cell_id.polyhedronIndex);
			//if(valeur > 0) std::cout << " density = " << valeur << " volume= " << volume_loc << " density= " << density_loc << endl;
			max_density_ = std::max(max_density_,density_loc);

			//geode_BlocName_to_density_.insert({(std::string)name,density_loc});
			density_+=density_loc;

		}


		volume_ += volume_loc;

    }
    if(volume_ > 0) {
    	density_ /= volume_;
    	density_ *=lfrac_mean;
    }
}
//---------------------------------------------------
// estimate the fractures density for the karstsim sgrid 
//---------------------------------------------------

void FractureSet::generate_fractures_density( bool density_variable_in_grid, double density_val,const geode::StructuralModel &model,std::vector< remfracres::CellId >& cells,double lfrac_mean, int density_type)
{
		
    if (density_variable_in_grid == true) density_ = 0.0 ;

	density_ = density_val;
	max_density_=density_val;

	if(density_type==1){

	   initialize_density_P31_from_structural_model(model,cells);

	}else{

	   initialize_density_P32_from_structural_model(model,cells,lfrac_mean);

	}




}


//---------------------------------------------------
// estimate the fractures density for the karstsim sgrid
//---------------------------------------------------

void FractureSet::generate_fractures_density_in_array( bool density_variable_in_grid, double density_val,const geode::StructuralModel &model,std::vector< remfracres::CellId >& cells,double lfrac_mean, int density_type)
{

    if (density_variable_in_grid == true) density_ = 0.0 ;

	density_ = density_val;
	max_density_=density_val;

	if(density_type==1){

	   initialize_density_P31_from_structural_model_in_array(model,cells);

	}else{

	   initialize_density_P32_from_structural_model_in_array(model,cells,lfrac_mean);

	}




}


    
//-------------------------------------------------------------
//  init the tsurf property    
//-------------------------------------------------------------


void FractureSet::init_tsurf(std::string name)
{

    name_ = name;


}   


    
//-------------------------------------------------------------
//  fill the fracture set with data of gocad fracture network    
//-------------------------------------------------------------


bool FractureSet::import_gocad_fracture_network(std::string filename)
{

    is_imported_ =true;



    std::fstream is(filename, std::ios_base::in);



    const char* code ="";

    int n_property = 0;
    int aperture_property = -1;
    float no_data_value_aperture = -99999;

    /*
    std::string class_names;

    while( code=catalog->code(in)) 
	{
        // read property info 
		if( Str::equal(code,"GOCAD"))
        {
            char class_name[80];
            in.width(80) ;
            in >> class_name;
            
            class_names = CString(class_name);

            if (class_names != "FRACTURENETWORK") return false;
            
        }

		// read property info 
		if( Str::equal(code,"HEADER") ||  Str::equal(code,"HEADER{" ))
        {
            bool find_name = false;

            catalog->eat_line(in);

            char c;
            while( catalog->get(in,c) ) 
            {
                while(isspace(c) && catalog->get(in,c) );
                
                char ident[80];
                register char * p=0;

                if( c == 'n' ) 
                {
                    ident[0]= c;
                    for( p = ident; catalog->get(in,c) && c != ':' && c != '\n'; *++p = c );
                } 
                
                *++p = '\0';

                if(Str::equal(ident,"name")) 
                {
                    char fract_name[80];
                    in >> fract_name;
            
                    init_tsurf(fract_name,karst_sim); 

                    hcs = fractures_->get_coordinate_system();

                    find_name = true;

                    break;                    

                }

                if( c == '}' ) break;

            }
        
            if (find_name == false) return false;
          
        }

        if( Str::equal(code,"PROPERTIES"))
        {
            CopyCString property_name;

            char c;
            while( catalog->get(in,c) ) 
            {

                while( c != '\n' && isspace(c) && catalog->get(in,c) );
                if( c == '\n' ) break;

                char ident[80];
                register char * p;

                if( c == '"' ) 
                {
                    for( p=ident; catalog->get(in,c) && c != '"' && c != '\n'; *p++ = c );
                } 
                else 
                {
                    for( p=ident, *p=c; !isspace(c) && catalog->get(in,c); *++p = c );
                }
                
                *p = '\0';

                if( p != ident ) 
                {
                    n_property ++;
                    if (CString(ident) == "APERTURE" || CString(ident) == "aperture" || CString(ident) == "Aperture")  aperture_property = n_property-1;
          
                }

                if( c == '\n' ) break;

            }
            

        }

        if( Str::equal(code,"NO_DATA_VALUES"))
        {
           
            float no_value; 

            int i;

            for (i = 0;i<n_property; i++)
            {
			    in >> no_value; 
                if (i == aperture_property) no_data_value_aperture = no_value;
            }
            
        }


        if( Str::equal(code,"SUBFRACTURENETWORK"))
        {
            // not taking into account for the moment 
        }



        if( Str::equal(code,"APFRACT"))
        {
            int id_fract;
            double x,y,z;

            in >> id_fract >> x >> y >> z;  // read fracture id and coordinates
            Point3d center(x,y,z);

            double x1, y1, z1, x2, y2, z2 ;
            in >> x1 >> y1 >> z1 >> x2 >> y2 >> z2 ;
           
            Vector3d v1(x1,y1,z1); 
            Vector3d v2(x2,y2,z2);

			double aperture = 1e-20;

            for (int i = 0;i<n_property; i++) 
            {
                in  >> aperture ;
                if (i == aperture_property) break; 
            }

            center=  hcs->transform_into(hcs->reference(), center);
            v1=  hcs->transform_into(hcs->reference(), v1);
            v2=  hcs->transform_into(hcs->reference(), v2);
            
            add_fracture_with_vector(center,v1,v2,aperture);
   
            fractures_number_++;
        }

        if( Str::equal(code,"PFRACT"))
        {
            int id_fract;
            double x,y,z;

            in >> id_fract >> x >> y >> z;  // read fracture id and coordinates
            Point3d center(x,y,z);

            double l1,l2,orientation,dip;
            in >> orientation >> dip >>  l1 >> l2 ;

            double aperture=1e-20;
            
            for (int i = 0;i<n_property; i++) 
            {
                in  >> aperture ;
                if (i == aperture_property) break; 
            }

			Point3d p0_uv(0.5, 0.5,0.0);

            center=  hcs->transform_into(hcs->reference(), center);
            // for the moment l1 et l2 are not changed 
           
            add_fracture(center, p0_uv, l1 , l2, orientation, dip, aperture, false,NULL,NULL,-1,NULL,NULL);

            fractures_number_++;

        }
           
            
    }


    if (fractures_ != NULL) 
    {


        evaluate_fractures_number_ = fractures_number_;
        fractures_activable_ = fractures_number_;
		return true;
    }
	else return false; 
*/

    return true;
}



//----------------------------------------------
// create a TFace and add it in the Tsurf
//----------------------------------------------

void FractureSet::add_fracture_with_vector(geode::Point3D center, geode::Vector3D v_major , geode::Vector3D v_minor,float aperture)
{

	// create the Tface for a given fracture  
   // TFace * tface = fractures_->create_element();

	// create the 4 atom of the Tface 
		
   /* Atom* a0 = tface-> create_atom(center-v_major-v_minor);
	tface-> add_atom(a0);
	Atom* a1 = tface-> create_atom(center-v_major+v_minor);
	tface-> add_atom(a1);
	Atom* a2 = tface-> create_atom(center+v_major+v_minor);
	tface-> add_atom(a2);
	Atom* a3 = tface-> create_atom(center+v_major-v_minor);
	tface-> add_atom(a3);

	// create the 2 triangles of the Tface
    Trgl * trgl;
    trgl = tface-> create_simplex(a0,a2,a1); 
	tface-> add_simplex(trgl);
	trgl = tface-> create_simplex(a0,a3,a2); 
	tface-> add_simplex(trgl);
		  
	  
	// add the fracture properties to each Tface atom
    for( PtrListItr<Atom> j(tface->atoms()); j.more(); j.next() ) 
	{
		Atom * a = j.cur() ;
		a->value(index_aper_) = aperture ; // initial value
		  
	}
	
	// add the Tface  to the fracture set TSurf
	fractures_->add(tface);*/
}

//----------------------------------------------
// Initialise the fractureset property to no datat 
// value 
//----------------------------------------------

void FractureSet::initialise_property_to_nd_value(std::string property_name)
{

	return;
  //  if (fractures_ == NULL) return ;

	/*
    int property_index = fractures_->datapack().find_property_record_index(property_name);

    if (property_index == -1) return;

    Property* property =  fractures_->datapack().find_property(property_name);

   NoDataValue nd_value = property->no_data_value();
   float valNd = nd_value.default_value();

    // iteration over the TSurf triangles 
	for( TSurfTrglsItr it(TSurfAPI::get_trgls(fractures_)); it.more(); it.next()) 
	{
	
		// get the current triangle
		Trgl * t = it.cur();

		// get the triangle atoms 
		Atom *a[3] = {NULL,NULL,NULL};
		t->atoms(a);
  
        a[0]->value(property_index) = valNd;
        a[1]->value(property_index) = valNd;
        a[2]->value(property_index) = valNd;
    }
	*/
}

//----------------------------------------------
// change de phase property for the next fractures  
//----------------------------------------------

void FractureSet::change_phase_property(int frac_number,int phase_number)
{

	return;
   // if (fractures_ == NULL) return ;

    /*
    int property_index = fractures_->datapack().find_property_record_index("Phase_number");

    if (property_index == -1) return;

	int begin = 2*fractures_inserted_;
	int end = begin + 2*frac_number;


    int size= -1;

    // iteration over the TSurf triangles 
	for( TSurfTrglsItr it(TSurfAPI::get_trgls(fractures_)); it.more(); it.next()) 
	{
	    size ++; 

        if (size < begin || size > end-1 ) continue;

		// get the current triangle
		Trgl * t = it.cur();

		// get the triangle atoms 
		Atom *a[3] = {NULL,NULL,NULL};
		t->atoms(a);
  
        a[0]->value(property_index) = phase_number;
        a[1]->value(property_index) = phase_number;
        a[2]->value(property_index) = phase_number;
    }
		*/
    fractures_inserted_ += frac_number;
}

//---------------------------------------------------
//   save  the fracture set data for the karstsim 
//---------------------------------------------------
	
 bool FractureSet::save(std::string filename)
 {

	 std::fstream fout(filename,fout.out);

	 if(fout.is_open()){

	fout << "FRACTURE_SET_NAME" << ";" << name_<< endl;
	fout << "FRACTURES_NUMBER" << ";" << get_fractures_number()<<endl;

	fout << "IS_IMPORTED" << ";" << is_imported_ <<endl;

	fout << "WELL_DATA_PROPERTY" <<  ";" << constraint_property_name_<<'"' << endl;

	fout << "DENSITY_PROPERTY"   <<  ";" << density_property_name_<<'"' << endl;
	
	fout.setf(std::ios::scientific);
    fout.precision(12);

	fout << "DENSITY_VALUE" <<  ";"  <<	density_<< endl;
	fout << "DENSITY_FUNCTION_TYPE"<<  ";"   << 0 << endl;
	fout << "DENSITY_COEF_A"  <<  ";"  << coef_a_ << endl;
	fout << "SEED "  <<  ";" <<	seed_ << endl;

	fout << "LENGTH1_INDEX_TYPE " <<  ";" << distribution_length1_type_index_<<  ";" << distr_length1_val_ <<  endl;

	fout << "DISTRIBUTION_NAME_LENGTH1 " <<  ";" << length_1_<<  ";" << minmaxmode_values_[0] <<  ";" << minmaxmode_values_[1]<<  ";" << minmaxmode_values_[2]  << endl;

	fout << "LENGTH2_INDEX_TYPE " <<  ";" << distribution_length2_type_index_<<  ";" << distr_length2_val_ <<  endl;

	fout << "DISTRIBUTION_NAME_LENGTH2 "<<  ";" << length_2_<<  ";" << minmaxmode_values_[3] <<  ";" << minmaxmode_values_[4]<<  ";" << minmaxmode_values_[5] << endl;

	fout << "DIP_AZIMUTH_INDEX_TYPE " <<  ";" << distribution_azimut_type_index_<<  ";" << distr_azimut_val_ <<  endl;
	
	fout << "DISTRIBUTION_NAME_DIP_AZIMUT " <<  ";" <<	dip_azimut_<<  ";" << minmaxmode_values_[6] <<  ";" << minmaxmode_values_[7]<<  ";" << minmaxmode_values_[8] << endl;

	fout << "DIP_INDEX_TYPE " <<  ";" << distribution_angle_type_index_<<  ";" << distr_angle_val_ <<  endl;
	
	fout << "DISTRIBUTION_NAME_DIP " <<  ";" << dip_angle_<<  ";" << minmaxmode_values_[9] <<  ";" << minmaxmode_values_[10]<<  ";" << minmaxmode_values_[11] << endl;

	fout << "APERTURE_INDEX_TYPE " <<  ";" << distribution_aperture_type_index_<<  ";" << distr_aperture_val_ <<  endl;

	fout << "DISTRIBUTION_NAME_APERTURE "  <<  ";" <<	aperture_name_<<  ";" << minmaxmode_values_[12] <<  ";" << minmaxmode_values_[13]<<  ";" << minmaxmode_values_[14] << endl;

	fout << "COMPLEMENTARY_SET_ACTIVE " << " " <<complementary_fracture_set_active_ <<endl;

	fout<<endl;

	fout.close();
      return true;
	 }
    return false;
}




 //---------------------------------------------------
 //   open a fracture set data for the karstsim
 //---------------------------------------------------
 bool FractureSet::open_loc(std::fstream& fin)
 {

      if(fin.is_open()){

 	//const char* code ;
 	std::string comment;
 	std::string separator;
 	bool read_ok = false;;

 	int index;

 	std::string name;

 	fin >> comment >> separator >> name;

 	name_ = name;

 	//fin >> comment >> comment >> index;

 	//is_shadow_zone_active_ =(bool) index;

 	fin >> comment >> comment >> index;

 	type_category_= static_cast <fractures_intersect::FracResIntersectionCategoryEnum>(index);

 	//fin >> comment >> comment >> index;

 	//type_intersection_= static_cast <fractures_intersect::FracResIntersectionTypeEnum>(index);


    for(int i = 0 ; i < 8 ; i++){

    	fin >> name >> comment >>  index ;
    	tab_type_intersection_[i]= static_cast <fractures_intersect::FracResIntersectionTypeEnum>(index);
    }

 	//directory_path_ = name;

 	fin >> comment >> separator >> default_stress_zone_max_;

 	is_shadow_zone_active_= false;
 	if(default_stress_zone_max_ > 0) is_shadow_zone_active_=true;


 	fin >> comment >> separator >> density_law_type_;

 	int n_frac;

 	fin >> comment >> separator >> n_frac;

 	set_fractures_number(n_frac);

     fin >> comment >> separator >> is_imported_;

 	fin >> comment >> separator >> name;
 	//catalog->read_quoted_strfing(fin,name);
 	constraint_property_name_ = name;

 	fin >> comment >> separator >> name;
 	//catalog->read_quoted_string(in,name);
 	density_property_name_ = name;

 	fin >> comment >> separator >> density_;
 	fin >> comment >> separator >> density_type_index_ ;
 	fin >> comment >> separator >> coef_a_;
 	fin >> comment >> separator >> seed_ ;

 	double min,max,mode, space_echelon;

 	fin >> comment >> separator >>distribution_length1_type_index_ >> separator >> distr_length1_val_;

 	fin >> comment >> separator >> name>> separator >> min >> separator >> max >> separator >> mode;
 	//catalog->read_quoted_string(in,name);
 	length_1_ = name;
 	minmaxmode_values_[0]=min;
 	minmaxmode_values_[1]=max;
 	minmaxmode_values_[2]=mode;

 	fin >> comment >> separator >>distribution_length2_type_index_ >> separator >> distr_length2_val_;

 	fin >> comment >> separator >> name >> separator >> min >> separator >> max >> separator >> mode;
 	//catalog->read_quoted_string(in,name);
 	length_2_ = name;
 	minmaxmode_values_[3]=min;
 	minmaxmode_values_[4]=max;
 	minmaxmode_values_[5]=mode;

 	fin >> comment >> comment >>distribution_azimut_type_index_ >> comment >> distr_azimut_val_;

 	fin >> comment >> comment >> name >> comment >> min >> comment >> max >> comment >> mode;
 	//catalog->read_quoted_string(in,name);
 	dip_azimut_ = name;
 	minmaxmode_values_[6]=min;
 	minmaxmode_values_[7]=max;
 	minmaxmode_values_[8]=mode;

 	fin >> comment >> comment >>distribution_angle_type_index_ >> comment >> distr_angle_val_;

 	fin >> comment >> comment >> name >> comment >> min >> comment >> max >> comment >> mode;
 	//catalog->read_quoted_string(in,name);
 	dip_angle_ = name;
 	minmaxmode_values_[9]=min;
 	minmaxmode_values_[10]=max;
 	minmaxmode_values_[11]=mode;

 	fin >> comment >> comment >>distribution_aperture_type_index_ >> comment >> distr_aperture_val_;

 	fin >> comment >> comment >> name >> comment >> min >> comment >> max >> comment >> mode;
 	//catalog->read_quoted_string(in,name);
 	aperture_name_ = name;
 	minmaxmode_values_[12]=min;
 	minmaxmode_values_[13]=max;
 	minmaxmode_values_[14]=mode;

 	fin >> comment >> comment >> index;
 	complementary_fracture_set_active_ = bool(index);

 	fin >> comment >> comment >> index;
    type_array_= array_type(index);
    fin >> comment >> comment >> space_echelon;
    echelon_space_ = space_echelon;



 	read_ok = true;

     return read_ok;
      }
      return false;
 }


//---------------------------------------------------
//   open a fracture set data for the karstsim
//---------------------------------------------------
bool FractureSet::open(std::string filename)
{

	 std::fstream fin(filename,fin.in);

     if(fin.is_open()){

    	 open_loc(fin);
     }
     return false;
}

void FractureSet::initialize_dip(std::vector<double>& new_dip){
	Dip_=new_dip;
}
void FractureSet::initialize_orientation(std::vector<double>& new_orientation){
	Orientation_= new_orientation;
}
void FractureSet::initialize_length1(std::vector<double>& new_length1){
	length1_ = new_length1;
}
void FractureSet::initialize_length2(std::vector<double>& new_length2){
	length2_ = new_length2;
}
void FractureSet::initialize_aperture(std::vector<double>& new_aperture){
	aperture_ = new_aperture;
}
void FractureSet::init_surface_geometry( int indice,const geode::TriangulatedSurface3D& surface , geode::TriangulatedSurfaceBuilder3D& builder )
{



	for(const auto pt : vertices_){
		builder.create_point(pt);
	}

	for(const auto tri : triangles_){
	    builder.create_triangle( tri );
	}

	std::shared_ptr< geode::VariableAttribute< double > > attribut =
			surface.polygon_attribute_manager()
			.find_or_create_attribute< geode::VariableAttribute, double >(
					"Identifier", ndvalue_ );
	std::shared_ptr< geode::VariableAttribute< double > > attributL1 =
			surface.polygon_attribute_manager()
			.find_or_create_attribute< geode::VariableAttribute, double >(
					"Length1", ndvalue_ );
	std::shared_ptr< geode::VariableAttribute< double > > attributL2 =
			surface.polygon_attribute_manager()
			.find_or_create_attribute< geode::VariableAttribute, double >(
					"Length2", ndvalue_ );
	std::shared_ptr< geode::VariableAttribute< double > > attributDip =
			surface.polygon_attribute_manager()
			.find_or_create_attribute< geode::VariableAttribute, double >(
					"Dip", ndvalue_ );
	std::shared_ptr< geode::VariableAttribute< double > > attributOrientation =
			surface.polygon_attribute_manager()
			.find_or_create_attribute< geode::VariableAttribute, double >(
					"Orientation", ndvalue_ );
	std::shared_ptr< geode::VariableAttribute< double > > attributAperture =
			surface.polygon_attribute_manager()
			.find_or_create_attribute< geode::VariableAttribute, double >(
					"Aperture", ndvalue_ );

	for(int i=0; i < triangles_.size() ; i++){
		attribut->set_value(i , indice);
		if(triangles_l1_.size() >0)  attributL1->set_value(i , triangles_l1_[i]);
		if(triangles_l2_.size() >0)attributL2->set_value(i , triangles_l2_[i]);
		if(triangles_orientation_.size() >0)attributOrientation->set_value(i , triangles_orientation_[i]);
		if(triangles_dip_.size() >0) attributDip->set_value(i , triangles_dip_[i]);
		if(triangles_aperture_.size() >0) attributAperture->set_value(i , triangles_aperture_[i]);
	}

}
void FractureSet::save_trinagulate_surface_to_ts(std::string filename, std::string name_surface){

	 std::ofstream fout(filename,fout.out);

     int z_sign{ 1 };
     double ndv = -9999.00;

	 if(fout.is_open()){
         fout << "GOCAD TSurf 1" << endl;
         fout << "HEADER {" << endl;
         fout << "name: " << name_surface << endl;
         fout << "mesh: " << "on" << endl;
         fout << "}" << endl;
         fout << "GOCAD_ORIGINAL_COORDINATE_SYSTEM" << endl;
         fout << "NAME " << " \" Pdgm_Epic_Local \" " << endl;
         fout << "PROJECTION " << "Unknown"<< endl;
         fout << "DATUM " << "Unknown"<< endl;
         fout << "AXIS_NAME " <<  "X" << " " << "Y" << " " <<"Z" << endl;
         fout << "AXIS_UNIT " << "m" << " " << "m" << " " << "m" << endl;
         fout << "ZPOSITIVE " << ( z_sign == 1 ? "Elevation" : "Depth" )<< endl;
         fout << "END_ORIGINAL_COORDINATE_SYSTEM" << endl;
         fout << "GEOLOGICAL_FEATURE "<< name_surface << endl;
         fout << "PROPERTIES"<< " "<<  "SnS_data_mismatch" << endl;
         fout << "PROP_LEGAL_RANGES" << " " << "**none**"<< " " <<"**none**"  << endl;
         fout << "NO_DATA_VALUES"<< " "<< ndv<< endl;
         fout << "PROPERTY_CLASSES"<< " " <<"sns_data_mismatch"<< endl;
         fout << "PROPERTY_KINDS"<< " " << "Length" << endl;
         fout << "PROPERTY_SUBCLASSES"<< " " <<"LINEARFUNCTION" << " "<< "Float" << " "<< 1 << " "<< 0 <<endl;
         fout<< "ESIZES"<<" "<<1<<endl;
         fout << "UNITS"<<" "<<"m"<<endl;
         fout << "PROPERTY_CLASS_HEADER" << " " << "X" << " "<< "{" << endl;
         fout << "kind:" << "X" << endl;
         fout<< "unit:" << "m" << endl;
         fout<< "pclip:" << 99 << endl;
         fout << "}" << endl;
         fout << "PROPERTY_CLASS_HEADER" << " " << "Y" << " "<< "{" << endl;
         fout << "kind:" << "Y" << endl;
         fout<< "unit:" << "m" << endl;
         fout<< "pclip:" << 99 << endl;
         fout << "}" << endl;
         fout << "PROPERTY_CLASS_HEADER" << " " << "Z" << " "<< "{" << endl;
         fout << "kind:" << "Z" << endl;
         fout<< "unit:" << "m" << endl;
         fout<< "pclip:" << 99 << endl;
         fout << "is_Z: on" << endl;
         fout << "}" << endl;
         fout << "PROPERTY_CLASS_HEADER" << " " << "sns_data_mismatch" << " "<< "{" << endl;
         fout << "kind:" << "Length" << endl;
         fout<< "unit:" << "m" << endl;
         fout<< "pclip:" << 99 << endl;
         fout << "}" << endl;
         fout<<"TFACE"<<endl;
         int node_count=0;
         for(const auto pt : vertices_){
        	 fout<<std::string("PVRTX ")<< node_count <<" " << pt.value(0)<< " " << pt.value(1) << " " << pt.value(2) << " " << 1<<std::endl;
        	 node_count++;
         }
         int tri_count=0;
         for(const auto tri : triangles_){
        	 fout<<std::string("TRGL ")<<  tri[0]<< " " << tri[1] << " " << tri[2] << std::endl;
        	 tri_count++;
         }
         fout<<"END"<<endl;
	 }
	 fout.close();
}

/**
 * @fn bool save_fracture_set(std::string)
 * @brief
 *
 * @param surface_filename
 * @return
 */
bool FractureSet::save_fracture_set(std::string surface_filename){
	auto surface = geode::TriangulatedSurface3D::create(
	geode::OpenGeodeTriangulatedSurface3D::impl_name_static() );
	auto builder = geode::TriangulatedSurfaceBuilder3D::create( *surface );
	init_surface_geometry(fracture_set_index_,*surface, *builder);
	std::string surface_final = absl::StrCat( surface_filename,name_);

	const auto output_file_native = absl::StrCat( directory_path_,surface_final, ".vtp");
	const auto output_file_surface_ts = absl::StrCat( directory_path_,surface_final, ".ts");

	save_trinagulate_surface_to_ts(output_file_surface_ts, name_);
	//geode::save_triangulated_surface( *surface, output_file_native );
	return true;
}

bool FractureSet::test_All_fractures_boolean() {

	geode::OpenGeodeModelLibrary::initialize();
	geode::OpenGeodeMeshLibrary::initialize();
	geode::GeosciencesExplicitLibrary::initialize();
	geode::GeosciencesIOModelLibrary::initialize();
	geode::GeosciencesIOMeshLibrary::initialize();
    geode::IOModelLibrary::initialize();

	std::string filename = "LS2.";
	std::string prefixe_fracture = "FractureSet_";
	std::string surface_filename =  absl::StrCat(prefixe_fracture,filename);
	auto surface = geode::TriangulatedSurface3D::create(
	geode::OpenGeodeTriangulatedSurface3D::impl_name_static() );
	auto builder = geode::TriangulatedSurfaceBuilder3D::create( *surface );

	const geode::StructuralModel& model = geode::load_structural_model(absl::StrCat(directory_path_, filename, "lso"));

	fractures_intersect::model frac_intersec_model;
	set_fractures_number(2);
   	fracture_set_index_=1;
	density_=0.01;
	for(int i=0; i < 10;i++){
		length1_.push_back(0.5*(i+1));
		length2_.push_back(0.1*(i+1));
		Dip_.push_back(10*(i+1));
		Orientation_.push_back(180*(i+1));
		aperture_.push_back(0.5*(i+1));
	}
	clean_vertices_and_triangles_array();
	set_generation_box(model.bounding_box());
	//fractures_generation(model,false);
	init_surface_geometry(fracture_set_index_,*surface, *builder);

	const auto output_file_native = absl::StrCat( directory_path_,surface_filename, "vtp");

	const auto output_file_surface_native = absl::StrCat( directory_path_,surface_filename, "dxf");

	// auto poly_surface = geode::convert_surface_mesh_into_polygonal_surface( *surface );


	//geode::save_polygonal_surface(*poly_surface, output_file_surface_native);

	geode::save_triangulated_surface( *surface, output_file_native );


	return true;
}

void FractureSet::set_volume(float val){
	volume_=val;
}

void FractureSet::set_area(float val){
	area_= val;
}




}

/*int main(int argc, char **argv) {
	using namespace geode;
	using namespace std;
	using namespace remfracres;
	using namespace fractures_intersect;
	std::cout << "Hello world";
	if (argc < 2) {

		try {
			detail::initialize_geosciences_io();
			detail::initialize_model_io();
			detail::initialize_mesh_io();
			std::string filename = "LS2.";
			std::string prefixe_fracture = "FractureSet_";
			std::string surface_filename =  absl::StrCat(prefixe_fracture,filename);
			std::string filenameParameter = absl::StrCat(absl::StrCat(remfracres::data_path, "LS2_parameter_fractures.txt"));

			auto surface = geode::TriangulatedSurface3D::create(
			geode::OpenGeodeTriangulatedSurface3D::impl_name_static() );
			auto builder = geode::TriangulatedSurfaceBuilder3D::create( *surface );

			const geode::StructuralModel& model = geode::load_structural_model(absl::StrCat(remfracres::data_path, filename, "lso"));
			remfracres::FractureSet obj;
			fractures_intersect::model frac_intersec_model;

			obj.open(filenameParameter);
			obj.clean_vertices_and_triangles_array();
			obj.set_generation_box(model.bounding_box());
			obj.fractures_generation(model,false);

			obj.init_surface_geometry(obj.fracture_set_index_,*surface, *builder);
			obj.directory_path_ = remfracres::data_path;

			const auto output_file_native = absl::StrCat( obj.directory_path_,surface_filename, "vtp");

			const auto output_file_surface_ts = absl::StrCat( obj.directory_path_,surface_filename, "ts");

			obj.save_trinagulate_surface_to_ts(output_file_surface_ts, obj.name_);


			geode::save_triangulated_surface( *surface, output_file_native );

			Logger::info("TEST SUCCESS");

			return 0;
		} catch (...) {
			return geode::geode_lippincott();
		}
	}
	return 0;
}*/

