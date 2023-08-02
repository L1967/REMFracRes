//----------------------------  FracResBeddingAperture.h  ---------------------------
// 
//
//  Description:	FracResBeddingAperture class use for build hypothesis fracturing
//					information
// 
//  Documents:		FracResBeddingAperture.h
//
//  Creation date  : 01/11/2022
//
//  Author:			Pascal Siegel - Rabah Namar
//
//  Copyright (C) 2022 by in2earth
//
//---------------------------- FracResBeddingAperture.h  ---------------------------
#pragma once

#ifndef FracResBeddingAperture_h
#define FracResBeddingAperture_h

// local includes
#include <remfracres/FracResDataTypeEnum.h>

// standard includes
#include <fstream>
#include <iostream>
#include <iomanip>
#include <memory>
#include <vector>
#include <string>



namespace remfracres{

/**
 *    FracResBeddingAperture.h contains the parameters for the build
 *    hypothesis fracturing
 * 
 *    @author Pascal Siegel - Rabah Namar - Olivier Jaquet 
 *    @date 2023
 */
 
 class FracResBeddingAperture
{

// variables

public:

	 std::string aperture_name_;

	 std::string parameter_name_;

	 std::string name_;

	 std::string unit_name_;

	 std::string discrete_property_name_;

	 FracResDataTypeEnum dataType_;

	 double value_;

	 std::string distribution_name_;

	 std::vector< std::pair< int,int > > discrete_property_selected_values_;

     int index_;

     std::vector< double > min_mode_max_;
	bool is_active_;


    std::string output_directory_path_;

    std::string output_prefixe_name_;

    int index_begin_;

    int index_end_;
	
// functions

public: 

	/**
	 *      Default constructor 
	 */

     FracResBeddingAperture();

	/**
	 *		Copy constructor
	 * 
	 *      @param from The value to copy to this object
	 */

     FracResBeddingAperture(const FracResBeddingAperture& from);

	/**
	 *		Destructor 
	 */
        
	~FracResBeddingAperture(){};

	/** 
     *      Assignment operator 
	 *
     *      @param from the value to assign to this object 
     * 
     *      @return A reference to this object 
     */
  
	FracResBeddingAperture& operator=(const FracResBeddingAperture& from);

	/** 
     *      copy function 
	 *
     *      @param from the value to assign to this object 
     */
  
	void copy(const FracResBeddingAperture& from);
	
	
};
// inline methods


 }
#endif  // FracResBeddingAperture_h
