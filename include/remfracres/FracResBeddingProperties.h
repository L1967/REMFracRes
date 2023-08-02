//----------------------------  FracResBeddingProperties.h  ---------------------------
// 
//
//  Description:	FracResBeddingProperties class use for build hypothesis fracturing
//					information
// 
//  Documents:		FracResBeddingProperties.h
//
//  Creation date  : 01/11/2022
//
//  Author:			Pascal Siegel - Rabah Namar
//
//  Copyright (C) 2022 by in2earth
//
//---------------------------- FracResBeddingProperties.h  ---------------------------
#pragma once

#ifndef FracResBeddingProperties_h
#define FracResBeddingProperties_h

// local includes
#include <remfracres/FracResBeddingAperture.h>
#include <remfracres/FracResDataTypeEnum.h>
#include <numeric/model.h>
// standard includes
#include <fstream>
#include <iostream>
#include <iomanip>
#include <memory>
#include <vector>
#include <string>

// boost includes 
#include <boost/random.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <Eigen/Core>
#include <Eigen/Geometry>


namespace remfracres{

/**
 *    FracResBeddingProperties.h contains the parameters for the build
 *    hypothesis fracturing
 * 
 *    @author Pascal Siegel - Rabah Namar - Olivier Jaquet 
 *    @date 2023
 */
 
 class FracResBeddingProperties
{

// variables

public:

	 std::string stage_name_;

	 std::string fracres_hypothesis_name_;

	 std::string name_;

	 std::string slope_parameter_name_;

	 FracResDataTypeEnum slope_data_type_;

	 std::string slope_unit_name_;

	 double slope_data_value_;

	 double slope_ratio_data_value_;

	 std::string bed_parrallel_slope_distribution_name_;

	 std::string bed_parrallel_slope_distribution_type_name_;

	 std::vector< double > min_mode_max_;

	 std::string bed_parrallel_discrete_property_name_;

	 std::string bed_parrallel_slope_prop_name_;

	 int propertIndex_;

	 bool is_region_;

	 bool bedParallelPropertyActive_;

	 std::string  region_name_;

    std::vector< FracResBeddingAperture* > bedding_list_;
	
    int index_fracture_set_geo_begin_;

    int nb_fracture_set_geo_;

    std::string output_directory_path_;

    std::string output_prefixe_name_;
	
// functions

public: 

	/**
	 *      Default constructor 
	 */

     FracResBeddingProperties();

	/**
	 *		Copy constructor
	 * 
	 *      @param from The value to copy to this object
	 */

     FracResBeddingProperties(const FracResBeddingProperties& from);

	/**
	 *		Destructor 
	 */
        


	/** 
     *      Assignment operator 
	 *
     *      @param from the value to assign to this object 
     * 
     *      @return A reference to this object 
     */
  
	FracResBeddingProperties& operator=(const FracResBeddingProperties& from);

	/** 
     *      copy function 
	 *
     *      @param from the value to assign to this object 
     */
  
	void copy(const FracResBeddingProperties& from);
	
	
};
// inline methods


 }
#endif  // FracResBeddingProperties_h
