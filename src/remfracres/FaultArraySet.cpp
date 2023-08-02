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


//----------------------------FaultArraySet.cpp  ---------------------------
// 
//  FaultArraySet contains the parameters for the simulation
//  of a Boolean model and the generated fractures
//  
//
//  Author:  Pascal Siegel - Rabah Namar - Olivier Jaquet 
//
//	Date: 2002-2003
//
//----------------------------  FaultArraySet.cpp  ---------------------------


// local includes


#include <remfracres/FaultArraySet.h>

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
#include <boost/math/distributions/lognormal.hpp>
#include <boost/math/distributions/complement.hpp>
// system includes                       

namespace remfracres {

using namespace geode;
using namespace std;
using namespace fractures_intersect;



//-----------------------------------------------
// default constructor of the class FractureSet
//-----------------------------------------------

FaultArraySet::FaultArraySet()
{
	std::string name_= "FaultArraySet0";

	length_1_ = "";

	length_2_  = "";
	
	dip_azimut_  = "";

	dip_angle_  = "";
	
	aperture_name_  = "";


	density_ = 1e-6;

	complementary_fracture_set_active_ = true;

	density_function_=NULL;

	coef_a_=0.0; 
		
	seed_=1337;

	faultArray_number_ = 0;

	faultArray_activable_ = 0;

	faultArray_inserted_ = 0;

	evaluate_faultArray_number_ = 0;

	//fractures_ = NULL;

	density_property_name_ = "none";
	constraint_property_name_ = "none"; 
	damage_zone_hanging_property_name_="none";
	damage_zone_footwall_property_name_="none";

	is_edit_ = false;

    is_imported_ = false;
	is_created_ = false;
	distribution_length1_type_index_ = 0;
	distribution_length2_type_index_ = 0;
	distribution_azimut_type_index_ = 0;
	distribution_angle_type_index_ = 0;
	distribution_aperture_type_index_ = 0;
	distribution_damage_zone_hanging_index_ = 0;
	distribution_damage_zone_footwall_index_ = 0;

	distr_length1_val_ = 100.0;
	distr_length2_val_ = 50;
	distr_azimut_val_ = 180.0;
	distr_angle_val_= 0.0;
	distr_aperture_val_= 1.0;
	distr_damage_zone_hanging_val_= 1.0;
	distr_damage_zone_footwall_val_= 1.0;

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
	minimum_space_type_ = distribution_type::uniform;
	damage_zone_hanging_type_ = distribution_type::uniform;
	damage_zone_footwall_type_ = distribution_type::uniform;

	minmaxmode_density_values_.resize(3);
	minmaxmode_density_values_[0]=0;
	minmaxmode_density_values_[1]=1;
	minmaxmode_density_values_[2]=0;


	minmaxmode_minimum_space_values_.resize(3, 0);
	minmaxmode_minimum_space_values_[1]=1;

	minmaxmode_damage_zone_hanging_values_.resize(3);
	minmaxmode_damage_zone_hanging_values_[0]=0;
	minmaxmode_damage_zone_hanging_values_[1]=1;
	minmaxmode_damage_zone_hanging_values_[2]=0;

	minmaxmode_damage_zone_footwall_values_.resize(3);
	minmaxmode_damage_zone_footwall_values_[0]=0;
	minmaxmode_damage_zone_footwall_values_[1]=1;
	minmaxmode_damage_zone_footwall_values_[2]=0;

	minmaxmode_values_.resize(15);

	for(int i= 0 ; i < 5; i++){
	 minmaxmode_values_[3*i]=0;
	 minmaxmode_values_[3*i + 1]=1;
	 minmaxmode_values_[3*i + 2]=0;
	}
	model_fractures_ = fractures_intersect::model();
	fault_array_set_index_=1;
	default_stress_zone_max_ = 0.0;
	directory_path_ = remfracres::data_path;
	is_shadow_zone_active_ =false;
	success_ = false;

	cylinder_step_ = 12;
	tab_type_intersection_.resize(9, fractures_intersect::FracResIntersectionTypeEnum::RANDOM);
	type_category_ = fractures_intersect::FracResIntersectionCategoryEnum::FAULTS;
	type_intersection_= fractures_intersect::FracResIntersectionTypeEnum::RANDOM;

	nb_failles_=0;
	nb_fracture_sets_=0;
	nb_fault_array_sets_ = 0;
	fracResClusterType_ = FracResClustersFaultTypeEnum::Null;
	usePreviousFracture_ = false;
	typeFamilyArray_ = array_type::unkown;
 	heightCorrelationType_ = FracResGeometryCorrelationTypeEnum::INDEPENDANT;
 	heightCorrelationValue_ = 0.0;
 	heightCorrelationVariable_ = false;
 	heightCorrelationProperty_name_ = "";

 	apertureCorrelationType_ = FracResGeometryCorrelationTypeEnum::INDEPENDANT;
    apertureCorrelationValue_ = 0.0;
 	apertureCorrelationVariable_ = false;
 	apertureCorrelationProperty_name_ = "";

 	faultCoreType_ = FracResDataTypeEnum::CONSTANT;
    faultCoreValue_ = 0.0;
 	faultCoreDistribution_name_ = "";
 	faultCoreProperty_name_ = "";

 	rate_reactive_ = 0.0;



	density_min_max_.resize(3, 0);


	density_distribution_name_ = "";
	density_min_max_[1]=1.0;

	minimum_space_distribution_name_ = "";
	minimum_space_property_name_ = "None";
	distribution_minimum_space_type_index_ = 0;
	distr_minimum_space_val_ = 0.5;

	typeCluster_ = cluster_type::no_cluster;
	damage_zone_hanging_distribution_name_ = "";
	damage_zone_footwall_distribution_name_ = "";
	is_active_ = false;
	model_geo_cont_index_ =0;
	output_prefixe_name_ ="";
	total_simulation_time_ = 0.0;
	total_writting_output_time_ = 0.0;
	length_max_ = 0.0;

}


//---------------------------------------------
//  copy constructor of the class FractureSet
//---------------------------------------------

FaultArraySet::FaultArraySet(const FaultArraySet& from)
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

	faultArray_number_ = from.faultArray_number_;

	evaluate_faultArray_number_ = from.evaluate_faultArray_number_;

	faultArray_activable_ = from.faultArray_activable_;

	faultArray_inserted_ = from.faultArray_inserted_;

	//fractures_ = from.fractures_;

	generation_box_ = from.generation_box_;

	density_property_name_ = from.density_property_name_;


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
	minmaxmode_values_ = from.minmaxmode_values_;

	model_fractures_ = from.model_fractures_;

	fault_array_set_index_=from.fault_array_set_index_;

	default_stress_zone_max_ = from.default_stress_zone_max_;

	directory_path_ = from.directory_path_;

	is_shadow_zone_active_ = from.is_shadow_zone_active_;

	success_ = from.success_;

	cylinder_step_ = from.cylinder_step_;

	tab_type_intersection_ = from.tab_type_intersection_;
	type_category_ = from.type_category_;
	type_intersection_ = from.type_intersection_;

	nb_failles_= from.nb_failles_;

	nb_fracture_sets_= from.nb_fracture_sets_;
	nb_fault_array_sets_ = from.nb_fault_array_sets_;
	fracturesSet_list_ = from.fracturesSet_list_;
	fracResClusterType_ = from.fracResClusterType_;
	usePreviousFracture_ = from.usePreviousFracture_;
	typeFamilyArray_ = from.typeFamilyArray_;
	heightCorrelationType_ = from.heightCorrelationType_;
 	heightCorrelationValue_ = from.heightCorrelationValue_;
 	heightCorrelationVariable_ = from.heightCorrelationVariable_;
 	heightCorrelationProperty_name_ = from.heightCorrelationProperty_name_;

 	apertureCorrelationType_ = from.apertureCorrelationType_;
    apertureCorrelationValue_ = from.apertureCorrelationValue_;
 	apertureCorrelationVariable_ = from.apertureCorrelationVariable_;
 	apertureCorrelationProperty_name_ = from.apertureCorrelationProperty_name_;

 	faultCoreType_ = from.faultCoreType_;
    faultCoreValue_ = from.faultCoreValue_;
 	faultCoreDistribution_name_ = from.faultCoreDistribution_name_;
 	faultCoreProperty_name_ = from.faultCoreProperty_name_;

 	rate_reactive_ = from.rate_reactive_;

 	minimum_space_between_array_law_type_ = from.minimum_space_between_array_law_type_;

	minimum_space_between_array_law_value_ = from.minimum_space_between_array_law_value_;

	minimum_space_between_array_distribution_name_ = from.minimum_space_between_array_distribution_name_;

    minimum_space_between_array_min_max_mode_ = from.minimum_space_between_array_min_max_mode_;

	minimum_space_between_array_property_name_ = from.minimum_space_between_array_property_name_;

	density_distribution_name_ = from.density_distribution_name_;
	density_min_max_=from.density_min_max_;


	minimum_space_distribution_name_ = from.minimum_space_distribution_name_;
	minimum_space_property_name_ = from.minimum_space_property_name_;
	distribution_minimum_space_type_index_ = from.distribution_minimum_space_type_index_;
	distr_minimum_space_val_ = from.distr_minimum_space_val_;
	minimum_space_type_= from.minimum_space_type_;
	minmaxmode_minimum_space_values_ = from.minmaxmode_minimum_space_values_;

	damage_zone_hanging_distribution_name_ = from.damage_zone_hanging_distribution_name_;
	damage_zone_hanging_val_ = from.damage_zone_hanging_val_;
	damage_zone_hanging_property_name_ = from.damage_zone_hanging_property_name_;
	distribution_damage_zone_hanging_index_ = from.distribution_damage_zone_hanging_index_;
	distr_damage_zone_hanging_val_ = from.distr_damage_zone_hanging_val_;
	minmaxmode_damage_zone_hanging_values_ = from.minmaxmode_damage_zone_hanging_values_;

	damage_zone_footwall_distribution_name_ = from.damage_zone_footwall_distribution_name_;
	damage_zone_footwall_val_ = from.damage_zone_footwall_val_;
	damage_zone_footwall_property_name_ = from.damage_zone_footwall_property_name_;
	distribution_damage_zone_footwall_index_ = from.distribution_damage_zone_footwall_index_;
	distr_damage_zone_footwall_val_ = from.distr_damage_zone_footwall_val_;
	minmaxmode_damage_zone_footwall_values_ = from.minmaxmode_damage_zone_footwall_values_;
	typeCluster_ = from.typeCluster_;
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

}


//-------------------------------------------
// destructor of the class FractureSet
//-------------------------------------------

FaultArraySet::~FaultArraySet()
{
}


//-----------------------------------------------
// assignment operator of the class FractureSet
//----------------------------------------------

FaultArraySet& FaultArraySet::operator=(const FaultArraySet& from)
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

	faultArray_number_ = from.faultArray_number_;

	evaluate_faultArray_number_ = from.evaluate_faultArray_number_;

	faultArray_activable_ = from.faultArray_activable_;

	faultArray_inserted_ = from.faultArray_inserted_;

//	fractures_ = from.fractures_;

	generation_box_ = from.generation_box_;

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

	minmaxmode_values_ = from.minmaxmode_values_;

	model_fractures_ = from.model_fractures_;

	fault_array_set_index_=from.fault_array_set_index_;

	default_stress_zone_max_ = from.default_stress_zone_max_;

	directory_path_ = from.directory_path_;

	is_shadow_zone_active_ = from.is_shadow_zone_active_;

	success_ = from.success_;

	cylinder_step_ = from.cylinder_step_;

	tab_type_intersection_ = from.tab_type_intersection_;

	type_category_ = from.type_category_;

	type_intersection_ = from.type_intersection_;

	nb_failles_= from.nb_failles_;

	nb_fracture_sets_= from.nb_fracture_sets_;
	nb_fault_array_sets_ = from.nb_fault_array_sets_;
	fracturesSet_list_ = from.fracturesSet_list_;
	fracResClusterType_ = from.fracResClusterType_;
	usePreviousFracture_ = from.usePreviousFracture_;
	typeFamilyArray_ = from.typeFamilyArray_;
	heightCorrelationType_ = from.heightCorrelationType_;
 	heightCorrelationValue_ = from.heightCorrelationValue_;
 	heightCorrelationVariable_ = from.heightCorrelationVariable_;
 	heightCorrelationProperty_name_ = from.heightCorrelationProperty_name_;

 	apertureCorrelationType_ = from.apertureCorrelationType_;
    apertureCorrelationValue_ = from.apertureCorrelationValue_;
 	apertureCorrelationVariable_ = from.apertureCorrelationVariable_;
 	apertureCorrelationProperty_name_ = from.apertureCorrelationProperty_name_;

 	faultCoreType_ = from.faultCoreType_;
    faultCoreValue_ = from.faultCoreValue_;
 	faultCoreDistribution_name_ = from.faultCoreDistribution_name_;
 	faultCoreProperty_name_ = from.faultCoreProperty_name_;

 	rate_reactive_ = from.rate_reactive_;
 	minimum_space_between_array_law_type_ = from.minimum_space_between_array_law_type_;

	minimum_space_between_array_law_value_ = from.minimum_space_between_array_law_value_;

	minimum_space_between_array_distribution_name_ = from.minimum_space_between_array_distribution_name_;

    minimum_space_between_array_min_max_mode_ = from.minimum_space_between_array_min_max_mode_;

	minimum_space_between_array_property_name_ = from.minimum_space_between_array_property_name_;

	density_distribution_name_ = from.density_distribution_name_;

	density_min_max_=from.density_min_max_;

	minimum_space_distribution_name_ = from.minimum_space_distribution_name_;
	minimum_space_property_name_ = from.minimum_space_property_name_;
	distribution_minimum_space_type_index_ = from.distribution_minimum_space_type_index_;
	distr_minimum_space_val_ = from.distr_minimum_space_val_;
	minimum_space_type_= from.minimum_space_type_;
	minmaxmode_minimum_space_values_ = from.minmaxmode_minimum_space_values_;

	constraint_property_name_ = from.constraint_property_name_;

	damage_zone_hanging_distribution_name_ = from.damage_zone_hanging_distribution_name_;
	damage_zone_hanging_val_ = from.damage_zone_hanging_val_;
	damage_zone_hanging_property_name_ = from.damage_zone_hanging_property_name_;
	distribution_damage_zone_hanging_index_ = from.distribution_damage_zone_hanging_index_;
	distr_damage_zone_hanging_val_ = from.distr_damage_zone_hanging_val_;
	minmaxmode_damage_zone_hanging_values_ = from.minmaxmode_damage_zone_hanging_values_;

	damage_zone_footwall_distribution_name_ = from.damage_zone_footwall_distribution_name_;
	damage_zone_footwall_val_ = from.damage_zone_footwall_val_;
	damage_zone_footwall_property_name_ = from.damage_zone_footwall_property_name_;
	distribution_damage_zone_footwall_index_ = from.distribution_damage_zone_footwall_index_;
	distr_damage_zone_footwall_val_ = from.distr_damage_zone_footwall_val_;
	minmaxmode_damage_zone_footwall_values_ = from.minmaxmode_damage_zone_footwall_values_;
	typeCluster_ = from.typeCluster_;
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
	return *this; 
}


//--------------------------------------------------------
//      copy function 
//--------------------------------------------------------
  
void FaultArraySet::copy(const FaultArraySet& from)
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

	faultArray_number_ = from.faultArray_number_;

	faultArray_activable_ = from.faultArray_number_;

	faultArray_inserted_ = 0;


	// tsurf copy 
	

	generation_box_ = from.generation_box_;

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

	minmaxmode_values_ = from.minmaxmode_values_;

	model_fractures_ = from.model_fractures_;

	fault_array_set_index_=from.fault_array_set_index_;
	default_stress_zone_max_ = from.default_stress_zone_max_;
	directory_path_ = from.directory_path_;
	is_shadow_zone_active_ = from.is_shadow_zone_active_;
	success_ = from.success_;
	cylinder_step_ = from.cylinder_step_;
	tab_type_intersection_ = from.tab_type_intersection_;
	type_category_ = from.type_category_;

	type_intersection_ = from.type_intersection_;
	nb_failles_= from.nb_failles_;

	nb_fracture_sets_= from.nb_fracture_sets_;
	nb_fault_array_sets_ = from.nb_fault_array_sets_;
	fracturesSet_list_ = from.fracturesSet_list_;
	fracResClusterType_ = from.fracResClusterType_;
	usePreviousFracture_ = from.usePreviousFracture_;
	typeFamilyArray_ = from.typeFamilyArray_;
	heightCorrelationType_ = from.heightCorrelationType_;
 	heightCorrelationValue_ = from.heightCorrelationValue_;
 	heightCorrelationVariable_ = from.heightCorrelationVariable_;
 	heightCorrelationProperty_name_ = from.heightCorrelationProperty_name_;

 	apertureCorrelationType_ = from.apertureCorrelationType_;
    apertureCorrelationValue_ = from.apertureCorrelationValue_;
 	apertureCorrelationVariable_ = from.apertureCorrelationVariable_;
 	apertureCorrelationProperty_name_ = from.apertureCorrelationProperty_name_;

 	faultCoreType_ = from.faultCoreType_;
    faultCoreValue_ = from.faultCoreValue_;
 	faultCoreDistribution_name_ = from.faultCoreDistribution_name_;
 	faultCoreProperty_name_ = from.faultCoreProperty_name_;

 	rate_reactive_ = from.rate_reactive_;

 	minimum_space_between_array_law_type_ = from.minimum_space_between_array_law_type_;

	minimum_space_between_array_law_value_ = from.minimum_space_between_array_law_value_;

	minimum_space_between_array_distribution_name_ = from.minimum_space_between_array_distribution_name_;

    minimum_space_between_array_min_max_mode_ = from.minimum_space_between_array_min_max_mode_;

	minimum_space_between_array_property_name_ = from.minimum_space_between_array_property_name_;

	density_distribution_name_ = from.density_distribution_name_;

	density_min_max_=from.density_min_max_;

	minimum_space_distribution_name_ = from.minimum_space_distribution_name_;
	minimum_space_property_name_ = from.minimum_space_property_name_;
	distribution_minimum_space_type_index_ = from.distribution_minimum_space_type_index_;
	distr_minimum_space_val_ = from.distr_minimum_space_val_;
	minimum_space_type_= from.minimum_space_type_;
	minmaxmode_minimum_space_values_ = from.minmaxmode_minimum_space_values_;

	constraint_property_name_ = from.constraint_property_name_;

	damage_zone_hanging_distribution_name_ = from.damage_zone_hanging_distribution_name_;
	damage_zone_hanging_val_ = from.damage_zone_hanging_val_;
	damage_zone_hanging_property_name_ = from.damage_zone_hanging_property_name_;
	distribution_damage_zone_hanging_index_ = from.distribution_damage_zone_hanging_index_;
	distr_damage_zone_hanging_val_ = from.distr_damage_zone_hanging_val_;
	minmaxmode_damage_zone_hanging_values_ = from.minmaxmode_damage_zone_hanging_values_;

	damage_zone_footwall_distribution_name_ = from.damage_zone_footwall_distribution_name_;
	damage_zone_footwall_val_ = from.damage_zone_footwall_val_;
	damage_zone_footwall_property_name_ = from.damage_zone_footwall_property_name_;
	distribution_damage_zone_footwall_index_ = from.distribution_damage_zone_footwall_index_;
	distr_damage_zone_footwall_val_ = from.distr_damage_zone_footwall_val_;
	minmaxmode_damage_zone_footwall_values_ = from.minmaxmode_damage_zone_footwall_values_;
	typeCluster_ = from.typeCluster_;
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
}

//--------------------------------------------------------
//      distribution from index
//--------------------------------------------------------

void FaultArraySet::set_distribution_type_from_index(int index, std::string name)
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
	    	   minimum_space_type_ =  get_distribution_type(name);
	            break;
	        default:
	            break;
	    }

}

FaultArraySet::distribution_type FaultArraySet::get_distribution_type(std::string name)
{
	distribution_type  type=distribution_type::uniform;

	if(name == "triangular"){

		type =  distribution_type::triangulated;
	}else if( name == "normal"){
		type =  distribution_type::normal;
	}else if( name == "lognormal"){
		type =  distribution_type::lognormal;
	}else if(name == "power"){
		type= distribution_type::power;
	}


	return type;

}


//------------------------------------------------------
// register in the geobase the TSurf of the FractureSet
//------------------------------------------------------

void FaultArraySet::show_graphics()
{	
}


//--------------------------
// set the domain limits
//--------------------------

void FaultArraySet::set_generation_box(geode::BoundingBox3D box)
{
	generation_box_ = box;
}
/**
 * @fn void initialize_damage_zone_distribution(boost::mt19937&)
 * @brief
 *
 * @param center_gen
 */
void FaultArraySet::initialize_damage_zone_hanging_distribution(boost::mt19937 &center_gen) {
	randDamageHangingZone_ = boost::bind(boost::random::uniform_real_distribution<>(minmaxmode_damage_zone_hanging_values_[0], minmaxmode_damage_zone_hanging_values_[1]),center_gen);
	if (damage_zone_hanging_type_ == distribution_type::triangulated) {
		randDamageHangingZone_ = boost::bind(
				boost::random::triangle_distribution<>(minmaxmode_damage_zone_hanging_values_[0],
						minmaxmode_damage_zone_hanging_values_[2], minmaxmode_damage_zone_hanging_values_[1]),
				center_gen);
	} else if (damage_zone_hanging_type_ == distribution_type::normal) {
		randDamageHangingZone_ = boost::bind(
				boost::random::normal_distribution<>(minmaxmode_damage_zone_hanging_values_[0],
						minmaxmode_damage_zone_hanging_values_[1]), center_gen);
	} else if (damage_zone_hanging_type_ == distribution_type::lognormal) {
		randDamageHangingZone_ = boost::bind(
				boost::random::lognormal_distribution<>(minmaxmode_damage_zone_hanging_values_[0],
						minmaxmode_damage_zone_hanging_values_[1]), center_gen);
	}else if (damage_zone_hanging_type_ == distribution_type::power){
		randDistPowerHangingWall_ = DistributionPower(minmaxmode_damage_zone_hanging_values_[0],minmaxmode_damage_zone_hanging_values_[1]);
	}
}

/**
 * @fn void initialize_damage_zone_distribution(boost::mt19937&)
 * @brief
 *
 * @param center_gen
 */
void FaultArraySet::initialize_damage_zone_footwall_distribution(boost::mt19937 &center_gen) {
	randDamageFootwallZone_ = boost::bind(boost::random::uniform_real_distribution<>(minmaxmode_damage_zone_footwall_values_[0], minmaxmode_damage_zone_footwall_values_[1]),center_gen);
	if (damage_zone_footwall_type_ == distribution_type::triangulated) {
		randDamageFootwallZone_ = boost::bind(
				boost::random::triangle_distribution<>(minmaxmode_damage_zone_footwall_values_[0],
						minmaxmode_damage_zone_footwall_values_[2], minmaxmode_damage_zone_footwall_values_[1]),
				center_gen);
	} else if (damage_zone_footwall_type_ == distribution_type::normal) {
		randDamageFootwallZone_ = boost::bind(
				boost::random::normal_distribution<>(minmaxmode_damage_zone_footwall_values_[0],
						minmaxmode_damage_zone_footwall_values_[1]), center_gen);
	} else if (damage_zone_footwall_type_ == distribution_type::lognormal) {
		randDamageFootwallZone_ = boost::bind(
				boost::random::lognormal_distribution<>(minmaxmode_damage_zone_footwall_values_[0],
						minmaxmode_damage_zone_footwall_values_[1]), center_gen);
	}else if (damage_zone_footwall_type_ == distribution_type::power){
		randDistPowerFootWall_ = DistributionPower(minmaxmode_damage_zone_footwall_values_[0],minmaxmode_damage_zone_footwall_values_[1]);
	}
}
/**
 * @fn void initialize_length1_distribution(boost::mt19937&)
 * @brief
 *
 * @param center_gen
 */
void FaultArraySet::initialize_minimum_space_distribution(boost::mt19937 &center_gen) {
	randMinimumSpace_ = boost::bind(boost::random::uniform_real_distribution<>(minmaxmode_minimum_space_values_[0], minmaxmode_minimum_space_values_[1]),center_gen);
	if (minimum_space_type_ == distribution_type::triangulated) {
		randMinimumSpace_ = boost::bind(
				boost::random::triangle_distribution<>(minmaxmode_minimum_space_values_[0],
						minmaxmode_minimum_space_values_[2], minmaxmode_minimum_space_values_[1]),
				center_gen);
	} else if (minimum_space_type_ == distribution_type::normal) {
		randMinimumSpace_ = boost::bind(
				boost::random::normal_distribution<>(minmaxmode_minimum_space_values_[0],
						minmaxmode_minimum_space_values_[1]), center_gen);
	} else if (minimum_space_type_ == distribution_type::lognormal) {
		randMinimumSpace_ = boost::bind(
				boost::random::lognormal_distribution<>(minmaxmode_minimum_space_values_[0],
						minmaxmode_minimum_space_values_[1]), center_gen);
	}else if (minimum_space_type_ == distribution_type::power){
		randDistPowerMinimumSpace_ = DistributionPower(minmaxmode_minimum_space_values_[0],minmaxmode_minimum_space_values_[1]);
	}
}



/**
 * @fn void initialize_length1_distribution(boost::mt19937&)
 * @brief
 *
 * @param center_gen
 */
void FaultArraySet::initialize_length1_distribution(boost::mt19937 &center_gen) {
	randLength1_ = boost::bind(boost::random::uniform_real_distribution<>(minmaxmode_values_[0], minmaxmode_values_[1]),center_gen);
	if (length1_type_ == distribution_type::triangulated) {
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
void FaultArraySet::initialize_lenght2_distribution( boost::mt19937 &center_gen) {
	randLength2_ = boost::bind(
			boost::random::uniform_real_distribution<>(minmaxmode_values_[3], minmaxmode_values_[4]),
			center_gen);
	if (length2_type_ == distribution_type::triangulated) {
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
void FaultArraySet::initialize_azimut_distribution( boost::mt19937 &center_gen) {
	randAzimut_ = boost::bind(boost::random::uniform_real_distribution<>(minmaxmode_values_[6], minmaxmode_values_[7]),center_gen);
	if (azimuth_type_ == distribution_type::triangulated) {
		randAzimut_ = boost::bind(boost::random::triangle_distribution<>(minmaxmode_values_[6],minmaxmode_values_[8], minmaxmode_values_[7]),center_gen);
	} else if (azimuth_type_ == distribution_type::normal) {
		randAzimut_ = boost::bind(boost::random::normal_distribution<>(minmaxmode_values_[6],minmaxmode_values_[7]), center_gen);
	} else if (azimuth_type_ == distribution_type::lognormal) {
		randAzimut_ = boost::bind(boost::random::lognormal_distribution<>(minmaxmode_values_[6],minmaxmode_values_[7]), center_gen);
	}else if (azimuth_type_ == distribution_type::power){
		randDistPowerAzimuth_ = DistributionPower(minmaxmode_values_[6],minmaxmode_values_[7]);
	}
}
/**
 * @fn void initialize_diap_angle_distribution(boost::mt19937&)
 * @brief
 *
 * @param center_gen
 */
void FaultArraySet::initialize_dip_angle_distribution( boost::mt19937 &center_gen) {
	randDip_ = boost::bind(
			boost::random::uniform_real_distribution<>(minmaxmode_values_[9],
					minmaxmode_values_[10]), center_gen);
	if (dip_angle_type_ == distribution_type::triangulated) {
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
void FaultArraySet::initialize_aperture_distribution( boost::mt19937 &center_gen) {
	randAperture_ = boost::bind(
			boost::random::uniform_real_distribution<>(minmaxmode_values_[12],
					minmaxmode_values_[13]), center_gen);
	if (aperture_type_ == distribution_type::triangulated) {
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
/**
 * @fn void initialize_density_distribution(boost::mt19937&)
 * @brief
 *
 * @param center_gen
 */
void FaultArraySet::initialize_density_distribution( boost::mt19937 &center_gen) {
	randDensity_ = boost::bind(boost::random::uniform_real_distribution<>(minmaxmode_density_values_[0], minmaxmode_density_values_[1]),center_gen);
		if (density_type_ == distribution_type::triangulated) {
			randDensity_ = boost::bind(
					boost::random::triangle_distribution<>(minmaxmode_density_values_[0],
							minmaxmode_density_values_[2], minmaxmode_density_values_[1]),
					center_gen);
		} else if (length1_type_ == distribution_type::normal) {
			randDensity_ = boost::bind(
					boost::random::normal_distribution<>(minmaxmode_density_values_[0],
							minmaxmode_density_values_[1]), center_gen);
		} else if (length1_type_ == distribution_type::lognormal) {
			randDensity_ = boost::bind(
					boost::random::lognormal_distribution<>(minmaxmode_density_values_[0],
							minmaxmode_density_values_[1]), center_gen);
		}else if (density_type_ == distribution_type::power){
			randDistPowerDensity_ = DistributionPower(minmaxmode_density_values_[0],minmaxmode_density_values_[1]);
		}
}
void FaultArraySet::initialize_fault_array_set_distribution(boost::mt19937& center_gen) {
	initialize_length1_distribution(center_gen);
	initialize_lenght2_distribution(center_gen);
	initialize_azimut_distribution(center_gen);
	initialize_dip_angle_distribution(center_gen);
	initialize_aperture_distribution(center_gen);
	initialize_density_distribution(center_gen);
	initialize_minimum_space_distribution(center_gen);
	initialize_damage_zone_hanging_distribution(center_gen);
	initialize_damage_zone_footwall_distribution(center_gen);
	fracturesSet_list_[0].initialize_echelon_space_density_distribution(center_gen);
}

float FaultArraySet::evaluate_fault_array_number(const geode::StructuralModel &model,std::vector<CellId> cells,bool density_variable_in_grid){
	
	float depthmax=0.0;

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

	nb_frac = evaluate_faultArray_number_;



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
		generate_fault_array_density(density_variable_in_grid,density_,model, cells,surface_mean,density_law_type_);
		ntot = int(volume_ * density_);
		ntot_P31 = int(volume_ * density_);
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
				if (ntot > max_v) evaluate_faultArray_number_ = std::numeric_limits<int>::max();
				if (ntot < min_v) evaluate_faultArray_number_ = std::numeric_limits<int>::min();
				return theta_max;
			}
		}else if(density_law_type_==2) {
			if (ntot_P31 < min_v || ntot_P31 > max_v  )
			{
				if (ntot_P31 > max_v) evaluate_faultArray_number_ = std::numeric_limits<int>::max();
				if (ntot_P31 < min_v) evaluate_faultArray_number_ = std::numeric_limits<int>::min();
				return theta_max_P31;
			}
		}
	}


	if(density_law_type_==1){
	  evaluate_faultArray_number_ = int(ntot);

	} else if(density_law_type_==2) {
		 evaluate_faultArray_number_ = int(ntot_P31);
		 theta_max = theta_max_P31;
	}else{
		theta_max=1.0;
	}


	return theta_max;
}

void FaultArraySet::evaluate_final_fault_array_number(boost::uniform_01<boost::mt19937>& center_uni){
	// generation of number_frac fractures according to a Poisson 
	// distribution law where ntot
	// is the Poisson parameter equal to mean number of fractures
	int number_mean_frac= int(-evaluate_faultArray_number_);
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

	set_fault_array_number(number_frac);
}

void FaultArraySet::initialize_geometrical_fractures_data() {
	// test if asociated tsurf object exist or not
	// if( fractures_ == NULL ) return;
	length1_.resize(faultArray_number_, 0);
	length2_.resize(faultArray_number_, 0);
	Orientation_.resize(faultArray_number_, 0);
	Dip_.resize(faultArray_number_, 0);
	aperture_.resize(faultArray_number_, 0);
	density_val_.resize(faultArray_number_, 0);
	minimum_space_val_.resize(faultArray_number_, 0);
	damage_zone_hanging_val_.resize(faultArray_number_, 0);
	damage_zone_footwall_val_.resize(faultArray_number_, 0);
	minimum_space_val_.resize(faultArray_number_, 0);
	fracturesSet_list_[0].echelon_space_val_.resize(faultArray_number_, 0);
}

double FaultArraySet::update_data_from_distribution(boost::mt19937& center_gen) {
	initialize_fault_array_set_distribution(center_gen);
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
	for (int ifrac = 0; ifrac < faultArray_number_; ifrac++) {
		//cout << randLength1_() << endl;
		if (distribution_length1_type_index_ == 1) {
			length1_[ifrac] = distr_length1_val_;
			 //randLength1_();
		} else {
			if(length1_type_ == distribution_type::power){
				length1_[ifrac] = randDistPowerLength1_.GcInvPower(center_gen);
			}else{
			     length1_[ifrac] = randLength1_();
			}
		}
		if (distribution_length2_type_index_ == 1) {
			length2_[ifrac] = distr_length2_val_;
		} else {
			if(length2_type_ == distribution_type::power){
							length2_[ifrac] = randDistPowerLength2_.GcInvPower(center_gen);
			}else{
			length2_[ifrac] = randLength2_();
			}
		}
		length_max_l1 = std::max((double) (length_max_l1), length1_[ifrac]);
		length_max_l2 = std::max((double) (length_max_l2), length2_[ifrac]);
		if (distribution_azimut_type_index_ == 1) {
			Orientation_[ifrac] = 90+distr_azimut_val_;
		} else {
			if(azimuth_type_ == distribution_type::power){
			  Orientation_[ifrac] = 90+ randDistPowerAzimuth_.GcInvPower(center_gen);
			}else{
			 Orientation_[ifrac] = 90+randAzimut_();
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
			Dip_[ifrac] = distr_angle_val_;
		} else {
			if(dip_angle_type_ == distribution_type::power){
				Dip_[ifrac] = randDistPowerDipAngle_.GcInvPower(center_gen);
			}else{
			Dip_[ifrac] = randDip_();
			}
		}
		if (distribution_aperture_type_index_ == 1) {
			aperture_[ifrac] = distr_aperture_val_;
		} else {
			if(aperture_type_ == distribution_type::power){
			 aperture_[ifrac] = randDistPowerAperture_.GcInvPower(center_gen);
			}else{
			 aperture_[ifrac] = randAperture_();
			}

		}
		if(density_law_type_ > 0){

			if (density_type_index_ == 0) {
				density_val_[ifrac] = density_;

			} else {
				if(density_type_ == distribution_type::power){
				 density_val_[ifrac] = randDistPowerDensity_.GcInvPower(center_gen);
				}else{
				 density_val_[ifrac] = randDensity_();
				}
			}
		}

		if (distribution_minimum_space_type_index_ == 0) {
			minimum_space_val_[ifrac] = distr_minimum_space_val_;

		} else {
			if(minimum_space_type_ == distribution_type::power){
			 minimum_space_val_[ifrac] = randDistPowerMinimumSpace_.GcInvPower(center_gen);
			}else{
			 minimum_space_val_[ifrac] = randMinimumSpace_();
			}
		}

		if (distribution_damage_zone_hanging_index_ == 0) {
			if(damage_zone_hanging_type_ == distribution_type::power){
				damage_zone_hanging_val_[ifrac] = randDistPowerHangingWall_.GcInvPower(center_gen);
			}else{
			  damage_zone_hanging_val_[ifrac] = randDamageHangingZone_();
			}
		} else {
			damage_zone_hanging_val_[ifrac] = distr_damage_zone_hanging_val_;
		}

		if (distribution_damage_zone_footwall_index_ == 0) {
			if(damage_zone_footwall_type_ == distribution_type::power){
			 damage_zone_footwall_val_[ifrac] = randDistPowerFootWall_.GcInvPower(center_gen);
			}else{
			 damage_zone_footwall_val_[ifrac] = randDamageFootwallZone_();
			}
		} else {
			damage_zone_footwall_val_[ifrac] = distr_damage_zone_footwall_val_;
		}

		if(typeFamilyArray_ == FaultArraySet::echelon){
			if (fracturesSet_list_[0].distribution_echelon_type_index_ == 0) {
				fracturesSet_list_[0].echelon_space_val_[ifrac] = fracturesSet_list_[0].echelon_space_;

			} else {
				if(fracturesSet_list_[0].echelon_space_type_ == distribution_type::power){
					fracturesSet_list_[0].echelon_space_val_[ifrac] = fracturesSet_list_[0].randDistPowerEchelonSpace_.GcInvPower(center_gen);
				}else{
					fracturesSet_list_[0].echelon_space_val_[ifrac] = fracturesSet_list_[0].randEchelon_space_();
				}
			}
		}

		job_evolution++;
	}
  return std::max(length_max_l1,length_max_l2);
}



int FaultArraySet::fault_array_fractures_evaluation(const geode::StructuralModel &model, std::shared_ptr< geode::AABBTree3D >& tree, std::shared_ptr< BoxAABBEvalDistance3D >& disteval,  std::vector<CellId> cells,bool estim_volume){
	bool density_variable_in_grid = false;

	if (density_property_name_ != "none" && density_type_index_ !=0)
	{
		density_variable_in_grid = true;

	}

	boost::mt19937  center_gen;
	boost::uint32_t seed;
	seed=seed_* (fault_array_set_index_+1);
    center_gen.seed(seed);

    float theta_max = evaluate_fault_array_number(model,cells,density_variable_in_grid);

    FractureSet& setFrac = fracturesSet_list_[0];
    setFrac.seed_ = seed_;
	setFrac.fracture_set_index_ += 0;
	setFrac.set_volume(volume_);
	setFrac.set_area(area_);
	std::shared_ptr<fault_array> test_array = NULL;
	setFrac.fractures_evaluate_generation(model, tree, disteval,cells, test_array, estim_volume);
	return faultArray_number_*setFrac.fractures_number_;
}

//----------------------------------------------------------
// generation array faults
//----------------------------------------------------------
/**
 * @fn void fault_array_generation(const geode::StructuralModel&, bool=false)
 * @brief
 *
 * @param model
 * @param estim_volume
 */
void FaultArraySet::fault_array_generation(const geode::StructuralModel &model, std::shared_ptr< geode::AABBTree3D >& tree,

	    std::shared_ptr< BoxAABBEvalDistance3D >& disteval, std::vector<CellId> cells,fractures_intersect::Ccontrainte_Geo& geo_contraint, bool estim_volume ){
	bool density_variable_in_grid = false;

	if (density_property_name_ != "none" && density_type_index_ !=0)
	{
		density_variable_in_grid = true;

	}

	boost::mt19937  center_gen;
	boost::uint32_t seed;
	seed=seed_* (fault_array_set_index_+1);
    center_gen.seed(seed);

    float theta_max = evaluate_fault_array_number(model,cells,density_variable_in_grid);

	// return if just a estimation of the fractures number
	if (estim_volume == true) return;

	// generation of number_frac fractures according to a Poisson
	// distribution law where ntot
	// is the Poisson parameter equal to mean number of fractures
	boost::uniform_01<boost::mt19937> center_uni(center_gen);

	if(density_law_type_>0) evaluate_final_fault_array_number(center_uni);

	geode::Point3D center_p, p0_uv;

	double length_max=update_data_from_distribution(center_gen);

	float l_w= std::abs(generation_box_.max().value(2)-generation_box_.min().value(2));
	geode::Vector3D box_length(generation_box_.min(), generation_box_.max());
	int index = 0 ;

	float* densityTab = NULL;
	// case with density variable in a grid
	if (density_variable_in_grid == true)
	{

		index= calculation_array_from_property(model,tree,disteval,cells,geo_contraint,estim_volume,center_uni);


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
		while ( index < faultArray_number_ && security < faultArray_number_*1000)
		{
			// random allocation of fracture center
			// in real x,y,z
			center_p.set_value(0, generation_box_.min().value(0)+center_uni()*box_length.value(0));
			center_p.set_value(1,generation_box_.min().value(1)+center_uni()*box_length.value(1));
			center_p.set_value(2,z_init + center_uni()* z_bloc);



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


			if ((center_p.value(2) > topo + length_max*0.5) || (center_p.value(2) < bottom - length_max*0.5) ) theta_depth=0.0;

			if  ((theta_uni < theta_depth &&  density_law_type_ >0 ) || density_law_type_==0)
			{
				//if( add_fault_array(center_p,p0_uv, length1_[index] , length2_[index], Orientation_[index],
    			//			 Dip_[index], aperture_[index], false, NULL,NULL) == true)
				//{
					fractures_intersect::Point3D germe_pt = fractures_intersect::Point3D(center_p.value(0),center_p.value(1),center_p.value(2));

					if(!fault_array_geo_set_->check_array_center(germe_pt)){
						security++;
						continue;
					}
					geode::index_t box_id;
					geode::Point3D nearest_point;
					double distance;
					geode::Point3D query({germe_pt.x(),germe_pt.y(),germe_pt.z()});
					std::tie( box_id, nearest_point, distance ) = tree->closest_element_box( query, *disteval );

					CellId id_cell = cells[box_id];

					update_values_of_array_index_from_property(id_cell,model,index);


					std::shared_ptr<fault_array> test_array = std::shared_ptr<fault_array>(new fault_array(fault_array_set_index_,germe_pt,Orientation_[index],Orientation_[index], length1_[index] , length2_[index],aperture_[index],minimum_space_val_[index]));


					for( std::shared_ptr<fracture_set_geo>  fractureset : geo_contraint.tab_fractures_set_){

						for( std::shared_ptr<Cfaille> faille : fractureset->tab_fractures_){

							test_array->check_and_intersected_fracture(faille);
						}


					}



					cout<<" Fracture set  on  Array N° "<< index + 1 << " creation " <<endl;
					int nb_success=0;
					bool relayChoice=false;
					//test_array->check_and_intersected_fracture()
					for(FractureSet& setFrac : fracturesSet_list_){

						setFrac.array_index_ = index + fault_array_set_index_;
						setFrac.seed_ = seed_;
						setFrac.fracture_set_index_ += 0;
						setFrac.set_volume(test_array->oriented_bb_->volume());
						setFrac.set_area(test_array->oriented_bb_->area());
						setFrac.set_generation_oriented_box(*test_array->oriented_bb_);
						setFrac.complementary_fracture_set_active_ = false;

						setFrac.fractures_generation_from_shadow_zone_in_array(model, tree,disteval,cells,test_array,estim_volume, fractures_intersect::Point3D(0,0,0),index);

						if(setFrac.is_shadow_zone_active_){
							setFrac.type_intersection_=fractures_intersect::FracResIntersectionTypeEnum::CROSSING;
							cout<<" Fracture set "<< setFrac.name_ << " is crossing with shadow zone" <<endl;

						}else{

							if(setFrac.type_array_== FractureSet::array_type::unkown){
								if(setFrac.type_intersection_ == fractures_intersect::FracResIntersectionTypeEnum::CROSSING){
									cout<<" Fracture set "<< setFrac.name_ << " is crossing with no shadow zone" <<endl;
								}else if(setFrac.type_intersection_ == fractures_intersect::FracResIntersectionTypeEnum::ABUTTING){
									cout<<" Fracture set "<< setFrac.name_ << " is abbuting with no shadow zone" <<endl;
								}else{
									cout<<" Fracture set "<< setFrac.name_ << " is randomly crossing and abbuting with no shadow zone" <<endl;
								}
							}else if(setFrac.type_array_== FractureSet::array_type::anastomosing){
								cout<<" Fracture set "<< setFrac.name_ << " is anastomosing" <<endl;
							}else if(setFrac.type_array_== FractureSet::array_type::echelon){
								cout<<" Fracture set "<< setFrac.name_ << " is echelon" <<endl;
							}else if(setFrac.type_array_== FractureSet::array_type::relay){
								cout<<" Fracture set "<< setFrac.name_ << " is relay" <<endl;
							}
						}

						relayChoice=setFrac.type_array_== FractureSet::array_type::relay;

						setFrac.fractures_activable_ = setFrac.fractures_number_;

						setFrac.fractures_inserted_ = 0;

						if(setFrac.success_ ){
							nb_success++;

						}
						test_array->tab_fractures_set_.push_back(setFrac.stress_zone_set_);

					}

					if(relayChoice) test_array->intersect_all_fracture_set(10);

					fault_array_geo_set_->add_array(test_array);

					index++;


			//	}



			}

			security++;

		} // end loop over the fractures


		if (security == faultArray_number_*1000) {
			if(faultArray_number_ > 0){
			std::cout <<" germe point of Array N° " << index + 1 <<" is contain inside available arrays and reach security number so calculation stop." << endl;
			}else{
				std::cout <<" The number of  Array is equal to zero. change density value" << endl;
			}

		}


	} // end other cases


	if(index == faultArray_number_) {
		success_ = true;
		cout << " Number of created array for set " << name_ << " is : " << index << endl;

	}

	// update the phase and the activable fracture

	faultArray_activable_ = faultArray_number_;

	faultArray_inserted_ = 0;
}



//---------------------------------------------------------------------
// generation corridor faults
//---------------------------------------------------------------------

/**
 * @fn void fault_array_generation(const geode::StructuralModel&, bool=false)
 * @brief
 *
 * @param model
 * @param estim_volume
 */
void FaultArraySet::fault_corridor_generation(const geode::StructuralModel &model, std::shared_ptr< geode::AABBTree3D >& tree,

	    std::shared_ptr< BoxAABBEvalDistance3D >& disteval,  std::vector<CellId> cells,fractures_intersect::Ccontrainte_Geo& geo_contraint,bool estim_volume ){
	bool density_variable_in_grid = false;

	if (density_property_name_ != "none" && density_type_index_ !=0)
	{
		density_variable_in_grid = true;

	}

	boost::mt19937  center_gen;
	boost::uint32_t seed;
	seed=seed_* (fault_array_set_index_+1);
    center_gen.seed(seed);

    float theta_max = evaluate_fault_array_number(model,cells,density_variable_in_grid);

	// return if just a estimation of the fractures number
	if (estim_volume == true) return;

	// generation of number_frac fractures according to a Poisson
	// distribution law where ntot
	// is the Poisson parameter equal to mean number of fractures
	boost::uniform_01<boost::mt19937> center_uni(center_gen);

	if(density_law_type_>0)evaluate_final_fault_array_number(center_uni);

	geode::Point3D center_p, p0_uv;

	double length_max=update_data_from_distribution(center_gen);

	float l_w= std::abs(generation_box_.max().value(2)-generation_box_.min().value(2));
	geode::Vector3D box_length(generation_box_.min(), generation_box_.max());
	int index = 0 ;

	float* densityTab = NULL;
	// case with density variable in a grid
	if (density_variable_in_grid == true)
	{

		index= calculation_array_corridor_from_property(model,tree,disteval,cells,geo_contraint,estim_volume,center_uni);


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
		while ( index < faultArray_number_ && security < faultArray_number_*1000)
		{
			// random allocation of fracture center
			// in real x,y,z
			center_p.set_value(0, generation_box_.min().value(0)+center_uni()*box_length.value(0));
			center_p.set_value(1,generation_box_.min().value(1)+center_uni()*box_length.value(1));
			center_p.set_value(2,z_init + center_uni()* z_bloc);



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


			if ((center_p.value(2) > topo + length_max*0.5) || (center_p.value(2) < bottom - length_max*0.5) ) theta_depth=0.0;

			if  ((theta_uni < theta_depth &&  density_law_type_ >0 ) || density_law_type_==0)
			{
				//if( add_fault_array(center_p,p0_uv, length1_[index] , length2_[index], Orientation_[index],
    			//			 Dip_[index], aperture_[index], false, NULL,NULL) == true)
				//{
					fractures_intersect::Point3D germe_pt = fractures_intersect::Point3D(center_p.value(0),center_p.value(1),center_p.value(2));

					if(!fault_array_geo_set_->check_array_center(germe_pt)){
						security++;
						continue;
					}

					geode::index_t box_id;
					geode::Point3D nearest_point;
					double distance;
					geode::Point3D query({germe_pt.x(),germe_pt.y(),germe_pt.z()});
					std::tie( box_id, nearest_point, distance ) = tree->closest_element_box( query, *disteval );

					CellId id_cell = cells[box_id];

					update_values_of_array_index_from_property(id_cell,model,index);

					std::shared_ptr<fault_array> test_array = std::shared_ptr<fault_array>(new fault_array(fault_array_set_index_,germe_pt,Orientation_[index],90.0
							, length1_[index] , length2_[index],aperture_[index],minimum_space_val_[index]));


					cout<<" Fracture set  on  cluster corridor N° "<< index + 1 << " creation " <<endl;
					int nb_success=0;
					bool relayChoice=false;
					for(FractureSet& setFrac : fracturesSet_list_){

						setFrac.array_index_ = index + fault_array_set_index_;
						setFrac.seed_ = seed_;

						setFrac.fracture_set_index_ += 0;
						setFrac.set_volume(test_array->oriented_bb_->volume());
						setFrac.set_area(area_);
						setFrac.set_generation_oriented_box(*test_array->oriented_bb_);

						setFrac.distribution_azimut_type_index_=0;
						setFrac.azimuth_type_ = FractureSet::distribution_type::uniform;
						setFrac.minmaxmode_values_[6] = test_array->azimut_ -10;
						setFrac.minmaxmode_values_[7] = test_array->azimut_ +10;

						setFrac.distribution_angle_type_index_=0;
						setFrac.dip_angle_type_ = FractureSet::distribution_type::uniform;
						setFrac.minmaxmode_values_[9] = 80;
						setFrac.minmaxmode_values_[10] = 100;
						setFrac.complementary_fracture_set_active_ = false;
						setFrac.fractures_generation_from_shadow_zone_in_corridor(model, tree,disteval,cells,test_array,estim_volume, fractures_intersect::Point3D(0,0,0));



						if(setFrac.is_shadow_zone_active_){
							setFrac.type_intersection_=fractures_intersect::FracResIntersectionTypeEnum::CROSSING;
							cout<<" Fracture set "<< setFrac.name_ << " is crossing with shadow zone" <<endl;

						}else{

							if(setFrac.type_array_== FractureSet::array_type::unkown){
								if(setFrac.type_intersection_ == fractures_intersect::FracResIntersectionTypeEnum::CROSSING){
									cout<<" Fracture set "<< setFrac.name_ << " is crossing with no shadow zone" <<endl;
								}else if(setFrac.type_intersection_ == fractures_intersect::FracResIntersectionTypeEnum::ABUTTING){
									cout<<" Fracture set "<< setFrac.name_ << " is abbuting with no shadow zone" <<endl;
								}else{
									cout<<" Fracture set "<< setFrac.name_ << " is randomly crossing and abbuting with no shadow zone" <<endl;
								}
							}
						}



						setFrac.fractures_activable_ = setFrac.fractures_number_;

						setFrac.fractures_inserted_ = 0;

						if(setFrac.success_ ){
							nb_success++;

						}
						test_array->tab_fractures_set_.push_back(setFrac.stress_zone_set_);

					}


					fault_array_geo_set_->add_array(test_array);

					index++;


			//	}



			}

			security++;

		} // end loop over the fractures


		if (security == faultArray_number_*1000) {
			if(faultArray_number_ > 0){
			std::cout <<" germe point of Array N° " << index + 1 <<" is contain inside available arrays and reach security number so calculation stop." << endl;
			}else{
				std::cout <<" The number of  Array is equal to zero. change density value" << endl;
			}

		}


	} // end other cases


	if(index == faultArray_number_) {
		success_ = true;
		cout << " Number of created array for set " << name_ << " is : " << index << endl;

	}

	// update the phase and the activable fracture

	faultArray_activable_ = faultArray_number_;

	faultArray_inserted_ = 0;
}


//---------------------------------------------------------------------
// generation corridor faults
//---------------------------------------------------------------------

/**
 * @fn void fault_zone_generation(const geode::StructuralModel&, bool=false)
 * @brief
 *
 * @param model
 * @param estim_volume
 */
void FaultArraySet::fault_zone_generation(const geode::StructuralModel &model, std::shared_ptr< geode::AABBTree3D >& tree,

	    std::shared_ptr< BoxAABBEvalDistance3D >& disteval, std::vector<CellId> cells,fractures_intersect::Ccontrainte_Geo& geo_contraint, bool estim_volume ){
	bool density_variable_in_grid = false;

	if (density_property_name_ != "none" && density_type_index_ !=0)
	{
		density_variable_in_grid = true;

	}

	boost::mt19937  center_gen;
	boost::uint32_t seed;
	seed=seed_* (fault_array_set_index_+1);
    center_gen.seed(seed);

    float theta_max = evaluate_fault_array_number(model,cells,density_variable_in_grid);

	// return if just a estimation of the fractures number
	if (estim_volume == true) return;

	// generation of number_frac fractures according to a Poisson
	// distribution law where ntot
	// is the Poisson parameter equal to mean number of fractures
	boost::uniform_01<boost::mt19937> center_uni(center_gen);

	if(density_law_type_>0)evaluate_final_fault_array_number(center_uni);

	geode::Point3D center_p, p0_uv;

	double length_max=update_data_from_distribution(center_gen);

	float l_w= std::abs(generation_box_.max().value(2)-generation_box_.min().value(2));
	geode::Vector3D box_length(generation_box_.min(), generation_box_.max());
	int index = 0 ;

	float* densityTab = NULL;
	// case with density variable in a grid
	if (density_variable_in_grid == true)
	{

		index= calculation_array_fault_zone_from_property(model,tree,disteval,cells,geo_contraint,estim_volume,center_uni);



	}
	// other cases
	else
	{


		float z_init= generation_box_.min().value(2)-(length_max/2); // utiliser la bouding box de la grille :

		float z_bloc= (l_w + length_max/2);

		p0_uv.set_value(0, 0.5);
		p0_uv.set_value(1, 0.5);
		p0_uv.set_value(2,0);


		int security = 0;


		// loop until all the fractures are generated
		while ( index < faultArray_number_ && security < faultArray_number_*1000)
		{
			// random allocation of fracture center
			// in real x,y,z
			center_p.set_value(0, generation_box_.min().value(0)+center_uni()*box_length.value(0));
			center_p.set_value(1,generation_box_.min().value(1)+center_uni()*box_length.value(1));
			center_p.set_value(2,z_init + center_uni()* z_bloc);



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


			if ((center_p.value(2) > topo + length_max*0.5) || (center_p.value(2) < bottom - length_max*0.5) ) theta_depth=0.0;

			if  ((theta_uni < theta_depth &&  density_law_type_ >0 ) || density_law_type_==0)
			{
				//if( add_fault_array(center_p,p0_uv, length1_[index] , length2_[index], Orientation_[index],
    			//			 Dip_[index], aperture_[index], false, NULL,NULL) == true)
				//{
					fractures_intersect::Point3D germe_pt = fractures_intersect::Point3D(center_p.value(0),center_p.value(1),center_p.value(2));

					if(!fault_array_geo_set_->check_array_center(germe_pt)){
						security++;
						continue;
					}
					geode::index_t box_id;
					geode::Point3D nearest_point;
					double distance;
					geode::Point3D query({germe_pt.x(),germe_pt.y(),germe_pt.z()});
					std::tie( box_id, nearest_point, distance ) = tree->closest_element_box( query, *disteval );

					CellId id_cell = cells[box_id];

					update_values_of_array_index_from_property(id_cell,model,index);


					std::shared_ptr<fault_array> test_array = std::shared_ptr<fault_array>(new fault_array(fault_array_set_index_,germe_pt,Orientation_[index],Dip_[index]
							, length1_[index] , length2_[index],aperture_[index],minimum_space_val_[index]));

					test_array->ratio_hangingwall_ = damage_zone_hanging_val_[index];
					test_array->ratio_footwall_ = damage_zone_footwall_val_[index];

					cout<<" Fracture set  on  cluster fault zone N° "<< index + 1 << " creation " <<endl;
					int nb_success=0;

					for(FractureSet& setFrac : fracturesSet_list_){

						setFrac.array_index_ = index + fault_array_set_index_;
						setFrac.seed_ = seed_;

						setFrac.fracture_set_index_ += 0;
						setFrac.set_generation_oriented_box(*test_array->oriented_bb_);

						//setFrac.update_bounding_box_for_fault_zone(test_array->thickness_);

						setFrac.distribution_azimut_type_index_=0;
						setFrac.azimuth_type_ = FractureSet::distribution_type::uniform;
						setFrac.minmaxmode_values_[6] = test_array->azimut_ - 5.0;//65.0;
						setFrac.minmaxmode_values_[7] = test_array->azimut_ + 5.0;//85.0; //85
						setFrac.distr_azimut_val_ = test_array->azimut_;
						setFrac.distribution_angle_type_index_=0;
						setFrac.dip_angle_type_ = FractureSet::distribution_type::uniform;
						setFrac.minmaxmode_values_[9] = test_array->pendage_ - 5.0;//20.0;
						setFrac.minmaxmode_values_[10] = test_array->pendage_ + 5.0;//45.0; //45
						setFrac.distr_angle_val_ = test_array->pendage_;
                        setFrac.complementary_fracture_set_active_ = false;
						setFrac.fractures_generation_from_shadow_zone_in_damage_fault_zone(model, tree,disteval,cells,test_array,estim_volume, test_array->ratio_hangingwall_, test_array->ratio_footwall_);



						if(setFrac.is_shadow_zone_active_){
							setFrac.type_intersection_=fractures_intersect::FracResIntersectionTypeEnum::CROSSING;
							cout<<" Fracture set "<< setFrac.name_ << " is crossing with shadow zone" <<endl;

						}else{

							 if(setFrac.type_array_== FractureSet::array_type::fault_zone){
								cout<<" Fracture set "<< setFrac.name_ << " is fault zone" <<endl;
							}
						}


						setFrac.fractures_activable_ = setFrac.fractures_number_;

						setFrac.fractures_inserted_ = 0;

						if(setFrac.success_ ){
							nb_success++;

						}
						test_array->tab_fractures_set_.push_back(setFrac.stress_zone_set_);

					}

					test_array->intersect_all_fracture_set(10);
					test_array->thickness_=0.0;

					test_array->update_box_from_damage_zone();

					fault_array_geo_set_->add_array(test_array);

					index++;


			//	}



			}

			security++;

		} // end loop over the fractures


		if (security == faultArray_number_*1000) {
			if(faultArray_number_ > 0){
			std::cout <<" germe point of Array N° " << index + 1 <<" is contain inside available arrays and reach security number so calculation stop." << endl;
			}else{
				std::cout <<" The number of  Array is equal to zero. change density value" << endl;
			}

		}


	} // end other cases


	if(index == faultArray_number_) {
		success_ = true;
		cout << " Number of created array for set " << name_ << " is : " << index << endl;

	}

	// update the phase and the activable fracture

	faultArray_activable_ = faultArray_number_;

	faultArray_inserted_ = 0;
}










//----------------------------------------------
// create a TFace and add it in the Tsurf
// if no intersection with conditionnal info
//----------------------------------------------

bool FaultArraySet::add_fault_array(geode::Point3D p0_xyz, geode::Point3D p0_uv , double l1 , double l2, double orientation,
    						   double dip, double aperture, bool check_intersect ,float* constraint_data)
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





	return true; 
}


//-------------------------------------------------
// get the TSurf Object representing the fractures
//-------------------------------------------------


//----------------------------------------------
// get the number of fractures
//----------------------------------------------

int FaultArraySet::get_fault_array_number()
{
	return faultArray_number_;
}

//----------------------------------------------
// set the number of fractures
//----------------------------------------------

void FaultArraySet::set_fault_array_number(int n_frac)
{
	
	faultArray_activable_ += n_frac - faultArray_number_;
	if (faultArray_activable_ < 0) faultArray_activable_ = 0;

	faultArray_number_ = n_frac;
	evaluate_faultArray_number_ = n_frac;
}


//----------------------------------------------
// get the evaluate number of fractures
//----------------------------------------------

int FaultArraySet::get_evaluate_fault_array_number()
{
	return evaluate_faultArray_number_;
}

//----------------------------------------------
// get the number of activable fractures
//----------------------------------------------

int FaultArraySet::get_fault_array_activables()
{
	return faultArray_activable_;
}


//----------------------------------------------
// desactivate a given number of fractures
//----------------------------------------------

int FaultArraySet::desactivate_fault_array(int n_fracture)
{
	int n = faultArray_activable_;

	faultArray_activable_ = faultArray_activable_ + n_fracture;

	if (faultArray_activable_ > faultArray_number_ )
	{
		faultArray_activable_= faultArray_number_;
		return faultArray_number_ - n;
	}
	else 
	{
		return n_fracture;
	}

}

//----------------------------------------------
// activate a given number of fractures
//----------------------------------------------

int FaultArraySet::activate_fault_array(int n_fracture)
{
	int n = faultArray_activable_;

	faultArray_activable_ = faultArray_activable_ - n_fracture;

	if (faultArray_activable_ > 0)
	{
		return n_fracture;
	}
	else 
	{
		faultArray_activable_= 0;
		return n;
	}
}


int FaultArraySet::calculation_array_from_property(const geode::StructuralModel &model, std::shared_ptr< geode::AABBTree3D >& tree,

			    std::shared_ptr< BoxAABBEvalDistance3D >& disteval, std::vector<CellId> cells,fractures_intersect::Ccontrainte_Geo& geo_contraint,bool estim_volume, boost::uniform_01<boost::mt19937>& center_uni){


	int security = 0;

	int index=0;
	geode::Point3D center_p;
	geode::Vector3D box_length(generation_box_.min(), generation_box_.max());
	// loop until all the fractures are generated
	while ( index < faultArray_number_  && security < faultArray_number_*1000)
	{

		// random allocation of fracture center
		// in real x,y,z
		center_p.set_value(0, generation_box_.min().value(0)+center_uni()*box_length.value(0));
		center_p.set_value(1, generation_box_.min().value(1)+center_uni()*box_length.value(1));
		center_p.set_value(2, generation_box_.min().value(2)+center_uni()*box_length.value(2));
		fractures_intersect::Point3D germe_pt = fractures_intersect::Point3D(center_p.value(0),center_p.value(1),center_p.value(2));

		if(!fault_array_geo_set_->check_array_center(germe_pt)){
			security++;
			continue;
		}
		geode::index_t box_id;
		geode::Point3D nearest_point;
		double distance;

		std::tie( box_id, nearest_point, distance ) = tree->closest_element_box( center_p, *disteval );
		// karstsim grid index
		// find nearest cell in grid to get cell density

		float center_density= cells_density_[box_id];



		if (max_density_*center_uni() < center_density)
		{

			CellId id_cell = cells[box_id];

			update_values_of_array_index_from_property(id_cell,model,index);


			std::shared_ptr<fault_array> test_array = std::shared_ptr<fault_array>(new fault_array(fault_array_set_index_,germe_pt,Orientation_[index],Orientation_[index], length1_[index] , length2_[index],aperture_[index],minimum_space_val_[index]));


				cout<<" Fracture set  on  Array N° "<< index + 1 << " creation " <<endl;
				int nb_success=0;
				bool relayChoice=false;
				//test_array->check_and_intersected_fracture()
				for(FractureSet& setFrac : fracturesSet_list_){

					setFrac.array_index_ = index + fault_array_set_index_;
					setFrac.seed_ = seed_;
					setFrac.fracture_set_index_ += 0;
					setFrac.set_volume(test_array->oriented_bb_->volume());
					setFrac.set_area(test_array->oriented_bb_->area());
					setFrac.set_generation_oriented_box(*test_array->oriented_bb_);
					setFrac.complementary_fracture_set_active_ = false;

					setFrac.fractures_generation_from_shadow_zone_in_array(model, tree,disteval,cells,test_array,estim_volume, fractures_intersect::Point3D(0,0,0),index);

					if(setFrac.is_shadow_zone_active_){
						setFrac.type_intersection_=fractures_intersect::FracResIntersectionTypeEnum::CROSSING;
						cout<<" Fracture set "<< setFrac.name_ << " is crossing with shadow zone" <<endl;

					}else{

						if(setFrac.type_array_== FractureSet::array_type::unkown){
							if(setFrac.type_intersection_ == fractures_intersect::FracResIntersectionTypeEnum::CROSSING){
								cout<<" Fracture set "<< setFrac.name_ << " is crossing with no shadow zone" <<endl;
							}else if(setFrac.type_intersection_ == fractures_intersect::FracResIntersectionTypeEnum::ABUTTING){
								cout<<" Fracture set "<< setFrac.name_ << " is abbuting with no shadow zone" <<endl;
							}else{
								cout<<" Fracture set "<< setFrac.name_ << " is randomly crossing and abbuting with no shadow zone" <<endl;
							}
						}else if(setFrac.type_array_== FractureSet::array_type::anastomosing){
							cout<<" Fracture set "<< setFrac.name_ << " is anastomosing" <<endl;
						}else if(setFrac.type_array_== FractureSet::array_type::echelon){
							cout<<" Fracture set "<< setFrac.name_ << " is echelon" <<endl;
						}else if(setFrac.type_array_== FractureSet::array_type::relay){
							cout<<" Fracture set "<< setFrac.name_ << " is relay" <<endl;
						}
					}

					relayChoice=setFrac.type_array_== FractureSet::array_type::relay;

					setFrac.fractures_activable_ = setFrac.fractures_number_;

					setFrac.fractures_inserted_ = 0;

					if(setFrac.success_ ){
						nb_success++;

					}
					test_array->tab_fractures_set_.push_back(setFrac.stress_zone_set_);

				}

				if(relayChoice) test_array->intersect_all_fracture_set(10);

				fault_array_geo_set_->add_array(test_array);

				index++;



		}

		security++;


	}
	return index;
}

int FaultArraySet::calculation_array_fault_zone_from_property(const geode::StructuralModel &model, std::shared_ptr< geode::AABBTree3D >& tree,

			    std::shared_ptr< BoxAABBEvalDistance3D >& disteval, std::vector<CellId> cells,fractures_intersect::Ccontrainte_Geo& geo_contraint,bool estim_volume, boost::uniform_01<boost::mt19937>& center_uni){


	int security = 0;

	int index=0;
	geode::Point3D center_p;
	geode::Vector3D box_length(generation_box_.min(), generation_box_.max());
	// loop until all the fractures are generated
	while ( index < faultArray_number_  && security < faultArray_number_*1000)
	{

		// random allocation of fracture center
		// in real x,y,z
		center_p.set_value(0, generation_box_.min().value(0)+center_uni()*box_length.value(0));
		center_p.set_value(1, generation_box_.min().value(1)+center_uni()*box_length.value(1));
		center_p.set_value(2, generation_box_.min().value(2)+center_uni()*box_length.value(2));
		fractures_intersect::Point3D germe_pt = fractures_intersect::Point3D(center_p.value(0),center_p.value(1),center_p.value(2));

		if(!fault_array_geo_set_->check_array_center(germe_pt)){
			security++;
			continue;
		}
		geode::index_t box_id;
		geode::Point3D nearest_point;
		double distance;

		std::tie( box_id, nearest_point, distance ) = tree->closest_element_box( center_p, *disteval );
		// karstsim grid index
		// find nearest cell in grid to get cell density

		float center_density= cells_density_[box_id];



		if (max_density_*center_uni() < center_density)
		{


			CellId id_cell = cells[box_id];

			update_values_of_array_index_from_property(id_cell,model,index);


			std::shared_ptr<fault_array> test_array = std::shared_ptr<fault_array>(new fault_array(fault_array_set_index_,germe_pt,Orientation_[index],Orientation_[index], length1_[index] , length2_[index],aperture_[index],minimum_space_val_[index]));


			test_array->ratio_hangingwall_ = damage_zone_hanging_val_[index];
			test_array->ratio_footwall_ = damage_zone_footwall_val_[index];

			cout<<" Fracture set  on  cluster fault zone N° "<< index + 1 << " creation " <<endl;
			int nb_success=0;

			for(FractureSet& setFrac : fracturesSet_list_){

				setFrac.array_index_ = index + fault_array_set_index_;
				setFrac.seed_ = seed_;

				setFrac.fracture_set_index_ += 0;
				setFrac.set_generation_oriented_box(*test_array->oriented_bb_);

				//setFrac.update_bounding_box_for_fault_zone(test_array->thickness_);

				setFrac.distribution_azimut_type_index_=0;
				setFrac.azimuth_type_ = FractureSet::distribution_type::uniform;
				setFrac.minmaxmode_values_[6] = test_array->azimut_ - 5.0;//65.0;
				setFrac.minmaxmode_values_[7] = test_array->azimut_ + 5.0;//85.0; //85
				setFrac.distr_azimut_val_ = test_array->azimut_;
				setFrac.distribution_angle_type_index_=0;
				setFrac.dip_angle_type_ = FractureSet::distribution_type::uniform;
				setFrac.minmaxmode_values_[9] = test_array->pendage_ - 5.0;//20.0;
				setFrac.minmaxmode_values_[10] = test_array->pendage_ + 5.0;//45.0; //45
				setFrac.distr_angle_val_ = test_array->pendage_;
                setFrac.complementary_fracture_set_active_ = false;
				setFrac.fractures_generation_from_shadow_zone_in_damage_fault_zone(model, tree,disteval,cells,test_array,estim_volume, test_array->ratio_hangingwall_, test_array->ratio_footwall_);



				if(setFrac.is_shadow_zone_active_){
					setFrac.type_intersection_=fractures_intersect::FracResIntersectionTypeEnum::CROSSING;
					cout<<" Fracture set "<< setFrac.name_ << " is crossing with shadow zone" <<endl;

				}else{

					 if(setFrac.type_array_== FractureSet::array_type::fault_zone){
						cout<<" Fracture set "<< setFrac.name_ << " is fault zone" <<endl;
					}
				}


				setFrac.fractures_activable_ = setFrac.fractures_number_;

				setFrac.fractures_inserted_ = 0;

				if(setFrac.success_ ){
					nb_success++;

				}
				test_array->tab_fractures_set_.push_back(setFrac.stress_zone_set_);

			}

			test_array->intersect_all_fracture_set(10);
			test_array->thickness_=0.0;

			test_array->update_box_from_damage_zone();

			fault_array_geo_set_->add_array(test_array);


				index++;



		}

		security++;


	}
	return index;
}


int FaultArraySet::calculation_array_corridor_from_property(const geode::StructuralModel &model, std::shared_ptr< geode::AABBTree3D >& tree,

			    std::shared_ptr< BoxAABBEvalDistance3D >& disteval, std::vector<CellId> cells,fractures_intersect::Ccontrainte_Geo& geo_contraint,bool estim_volume, boost::uniform_01<boost::mt19937>& center_uni){


	int security = 0;

	int index=0;
	geode::Point3D center_p;
	geode::Vector3D box_length(generation_box_.min(), generation_box_.max());
	// loop until all the fractures are generated
	while ( index < faultArray_number_  && security < faultArray_number_*1000)
	{

		// random allocation of fracture center
		// in real x,y,z
		center_p.set_value(0, generation_box_.min().value(0)+center_uni()*box_length.value(0));
		center_p.set_value(1, generation_box_.min().value(1)+center_uni()*box_length.value(1));
		center_p.set_value(2, generation_box_.min().value(2)+center_uni()*box_length.value(2));
		fractures_intersect::Point3D germe_pt = fractures_intersect::Point3D(center_p.value(0),center_p.value(1),center_p.value(2));

		if(!fault_array_geo_set_->check_array_center(germe_pt)){
			security++;
			continue;
		}
		geode::index_t box_id;
		geode::Point3D nearest_point;
		double distance;

		std::tie( box_id, nearest_point, distance ) = tree->closest_element_box( center_p, *disteval );
		// karstsim grid index
		// find nearest cell in grid to get cell density

		float center_density= cells_density_[box_id];



		if (max_density_*center_uni() < center_density)
		{


			CellId id_cell = cells[box_id];

			update_values_of_array_index_from_property(id_cell,model,index);


			std::shared_ptr<fault_array> test_array = std::shared_ptr<fault_array>(new fault_array(fault_array_set_index_,germe_pt,Orientation_[index],Orientation_[index], length1_[index] , length2_[index],aperture_[index],minimum_space_val_[index]));



			cout<<" Fracture set  on  cluster corridor N° "<< index + 1 << " creation " <<endl;
			int nb_success=0;
			bool relayChoice=false;
			for(FractureSet& setFrac : fracturesSet_list_){

				setFrac.array_index_ = index + fault_array_set_index_;
				setFrac.seed_ = seed_;

				setFrac.fracture_set_index_ += 0;
				setFrac.set_volume(test_array->oriented_bb_->volume());
				setFrac.set_area(area_);
				setFrac.set_generation_oriented_box(*test_array->oriented_bb_);

				setFrac.distribution_azimut_type_index_=0;
				setFrac.azimuth_type_ = FractureSet::distribution_type::uniform;
				setFrac.minmaxmode_values_[6] = test_array->azimut_ -10;
				setFrac.minmaxmode_values_[7] = test_array->azimut_ +10;

				setFrac.distribution_angle_type_index_=0;
				setFrac.dip_angle_type_ = FractureSet::distribution_type::uniform;
				setFrac.minmaxmode_values_[9] = 80;
				setFrac.minmaxmode_values_[10] = 100;
				setFrac.complementary_fracture_set_active_ = false;
				setFrac.fractures_generation_from_shadow_zone_in_corridor(model, tree,disteval,cells,test_array,estim_volume, fractures_intersect::Point3D(0,0,0));



				if(setFrac.is_shadow_zone_active_){
					setFrac.type_intersection_=fractures_intersect::FracResIntersectionTypeEnum::CROSSING;
					cout<<" Fracture set "<< setFrac.name_ << " is crossing with shadow zone" <<endl;

				}else{

					if(setFrac.type_array_== FractureSet::array_type::unkown){
						if(setFrac.type_intersection_ == fractures_intersect::FracResIntersectionTypeEnum::CROSSING){
							cout<<" Fracture set "<< setFrac.name_ << " is crossing with no shadow zone" <<endl;
						}else if(setFrac.type_intersection_ == fractures_intersect::FracResIntersectionTypeEnum::ABUTTING){
							cout<<" Fracture set "<< setFrac.name_ << " is abbuting with no shadow zone" <<endl;
						}else{
							cout<<" Fracture set "<< setFrac.name_ << " is randomly crossing and abbuting with no shadow zone" <<endl;
						}
					}
				}


				setFrac.fractures_activable_ = setFrac.fractures_number_;

				setFrac.fractures_inserted_ = 0;

				if(setFrac.success_ ){
					nb_success++;

				}
				test_array->tab_fractures_set_.push_back(setFrac.stress_zone_set_);

			}



			fault_array_geo_set_->add_array(test_array);

			index++;




		}

		security++;


	}
	return index;
}











//---------------------------------------------------
// estimate the fractures density for the karstsim sgrid
//---------------------------------------------------

void FaultArraySet::initialize_density_P31_from_structural_model( const geode::StructuralModel &model,std::vector<CellId>& cells)
{
    geode_BlocName_to_BlockId_.clear();
    geode_BlocName_to_Density_exist_.clear();
    geode_BlocName_to_volume_.clear();
    geode_BlocName_to_area_.clear();
    geode_BlocName_to_density_.clear();
    double density_loc =0;
    double volume_loc=0;
    double density_tot=0;

    for( const auto& cell_id : cells)
    {
		const auto name = model.block(cell_id.blockId).name();
		const auto& id = cell_id.blockId;
		geode_BlocName_to_BlockId_.insert({(std::string)name,cell_id.blockId});

		const auto& mesh = model.block(cell_id.blockId).mesh<geode::TetrahedralSolid3D>();
		bool testDensity = mesh.polyhedron_attribute_manager().attribute_exists(density_property_name_);

		if(!testDensity){
			cells_density_.push_back(0.0);
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
			cells_density_.push_back(density_loc);

			max_density_ = std::max(max_density_,(float)attributDensity_->value(cell_id.polyhedronIndex));

			//geode_BlocName_to_density_.insert({(std::string)name,density_loc});
			density_+=density_loc;

		}

		volume_loc = geode::tetrahedron_volume(mesh.tetrahedron(cell_id.polyhedronIndex));
		volume_ += volume_loc;

    }

}
//---------------------------------------------------
// estimate the fractures density for the karstsim sgrid
//---------------------------------------------------

void FaultArraySet::initialize_density_P32_from_structural_model( const geode::StructuralModel &model, std::vector<CellId>& cells, double lfrac_mean)
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

    for( const auto& cell_id : cells)
    {
		const auto name = model.block(cell_id.blockId).name();
		const auto& id = cell_id.blockId;
		geode_BlocName_to_BlockId_.insert({(std::string)name,cell_id.blockId});

		const auto& mesh = model.block(cell_id.blockId).mesh<geode::TetrahedralSolid3D>();
		bool testDensity = mesh.polyhedron_attribute_manager().attribute_exists(density_property_name_);

		if(!testDensity){
			cells_density_.push_back(0.0);
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
			cells_density_.push_back(density_loc);

			max_density_ = std::max(max_density_,(float)attributDensity_->value(cell_id.polyhedronIndex));

			//geode_BlocName_to_density_.insert({(std::string)name,density_loc});
			density_+=density_loc;

		}

		volume_loc = geode::tetrahedron_volume(mesh.tetrahedron(cell_id.polyhedronIndex));
		volume_ += volume_loc;

    }

}
//---------------------------------------------------
// estimate the fractures density for the karstsim sgrid 
//---------------------------------------------------

void FaultArraySet::generate_fault_array_density( bool density_variable_in_grid, double density_val,const geode::StructuralModel &model,std::vector< remfracres::CellId >& cells,double lfrac_mean, int density_type)
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


//----------------------------------------------
// Initialise the fractureset property to no datat 
// value 
//----------------------------------------------

void FaultArraySet::initialise_property_to_nd_value(std::string property_name)
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

void FaultArraySet::change_phase_property(int frac_number,int phase_number)
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
    faultArray_inserted_ += frac_number;
}

//---------------------------------------------------
//   save  the fracture set data for the karstsim 
//---------------------------------------------------
	
 bool FaultArraySet::save(std::string filename)
 {

	 std::fstream fout(filename,fout.out);

	 if(fout.is_open()){


	fout << "FAULT_ARRAY_SET_NAME" << "; " << name_<< endl;

	fout << "CATEGORY_TYPE " << "; " << type_category_<< endl;
	fout << "MECHANIAL_UNITS "<< "; " << tab_type_intersection_[0]<< endl;
	fout << "BEDDINGS "<< "; " <<  tab_type_intersection_[0]<< endl;
	fout << "JOINTS  "<< "; " <<  tab_type_intersection_[0]<< endl;
	fout << "FAULTS  "<< "; " <<  tab_type_intersection_[0]<< endl;
	fout << "FAULT_ARRAYS  "<< "; " <<  tab_type_intersection_[0]<< endl;
	fout << "FRACTURE_CLUSTERS_JOINT_SWARM  "<< "; " << tab_type_intersection_[0]<< endl;
	fout << "FRACTURE_CLUSTERS_CORRIDORS_ZONES "<< "; " <<  tab_type_intersection_[0]<< endl;
	fout << "DEFAULT_ARRAY_ZONE_MAX " << "; " << is_shadow_zone_active_<< endl;

	fout << "P_LAW_TYPE" << ";" << density_type_index_<<endl;
	fout << "FAULT_ARRAY_NUMBER" << ";" << get_fault_array_number()<<endl;

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
 bool FaultArraySet::open_loc(std::fstream& fin)
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


    for(int i = 0 ; i < 7 ; i++){

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

 	set_fault_array_number(n_frac);

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

 	double min,max,mode;

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

	fin >> comment >> separator >> nb_fracture_sets_;
	//fin >> comment >> separator >> index;
	//intersect_inside_set_ = bool(index);

	for(int i =0; i < nb_fracture_sets_;i++ ){

		remfracres::FractureSet obj;
		obj.open_loc(fin);
		obj.fracture_set_index_=i;
		obj.cylinder_step_=cylinder_step_;
		fracturesSet_list_.push_back(obj);
	}
 	read_ok = true;

     return read_ok;
      }
      return false;
 }


//---------------------------------------------------
//   open a fracture set data for the karstsim
//---------------------------------------------------
bool FaultArraySet::open(std::string filename)
{

	 std::fstream fin(filename,fin.in);

     if(fin.is_open()){

    	 open_loc(fin);
     }
     return false;
}

void FaultArraySet::initialize_dip(std::vector<double>& new_dip){
	Dip_=new_dip;
}
void FaultArraySet::initialize_orientation(std::vector<double>& new_orientation){
	Orientation_= new_orientation;
}
void FaultArraySet::initialize_length1(std::vector<double>& new_length1){
	length1_ = new_length1;
}
void FaultArraySet::initialize_length2(std::vector<double>& new_length2){
	length2_ = new_length2;
}
void FaultArraySet::initialize_aperture(std::vector<double>& new_aperture){
	aperture_ = new_aperture;
}
void FaultArraySet::init_surface_geometry( int indice,const geode::TriangulatedSurface3D& surface , geode::TriangulatedSurfaceBuilder3D& builder )
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
		attributL1->set_value(i , triangles_l1_[i]);
		attributL2->set_value(i , triangles_l2_[i]);
		attributOrientation->set_value(i , triangles_orientation_[i]);
		attributDip->set_value(i , triangles_dip_[i]);
		attributAperture->set_value(i , triangles_aperture_[i]);
	}

}


void FaultArraySet::set_volume(float val){
	volume_=val;
}

void FaultArraySet::set_area(float val){
	area_= val;
}



void FaultArraySet::update_values_of_array_index_from_property(CellId id_cell,const geode::StructuralModel &model, int ifrac){
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

	if (distribution_damage_zone_hanging_index_ > 1) {

		bool testPoperty =block_mesh.polyhedron_attribute_manager().attribute_exists(damage_hanging_wall_property_name_);

		if (!testPoperty) {
			cout << " The Property " << damage_hanging_wall_property_name_
					<< " doesn't exist for block " << model.block(id_cell.blockId).name() << endl;
		}else{
		// get petrophysique properties from polyhedron
		const auto &attributPropertyHanging =block_mesh.polyhedron_attribute_manager().find_or_create_attribute<geode::VariableAttribute, double>(damage_hanging_wall_property_name_,ndvalue_);
		damage_zone_hanging_val_[ifrac] = attributPropertyHanging->value(id_cell.polyhedronIndex);
		}
	}
	if (distribution_damage_zone_footwall_index_ > 1) {

		bool testPoperty =block_mesh.polyhedron_attribute_manager().attribute_exists(damage_foot_wall_property_name_);
		if (!testPoperty) {
			cout << " The Property " << damage_foot_wall_property_name_
					<< " doesn't exist for block " << model.block(id_cell.blockId).name() << endl;
		}else{
			// get petrophysique properties from polyhedron
			const auto &attributPropertyFootwall =block_mesh.polyhedron_attribute_manager().find_or_create_attribute<geode::VariableAttribute, double>(damage_foot_wall_property_name_,ndvalue_);
			damage_zone_footwall_val_[ifrac] = attributPropertyFootwall->value(id_cell.polyhedronIndex);
		}
	}

	if (distribution_minimum_space_type_index_ > 1) {

		bool testPoperty =block_mesh.polyhedron_attribute_manager().attribute_exists(minimum_space_property_name_);
		if (!testPoperty) {
			cout << " The Property " << minimum_space_property_name_
					<< " doesn't exist for block " << model.block(id_cell.blockId).name() << endl;
		}else{
			// get petrophysique properties from polyhedron
			const auto &attributPropertyFootwall =block_mesh.polyhedron_attribute_manager().find_or_create_attribute<geode::VariableAttribute, double>(minimum_space_property_name_,ndvalue_);
			minimum_space_val_[ifrac] = attributPropertyFootwall->value(id_cell.polyhedronIndex);
		}
	}

}

}


