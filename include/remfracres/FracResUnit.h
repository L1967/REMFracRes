//----------------------------  FracResUnit.h  ---------------------------
// 
//
//  Description:	FracResUnit class use for build DFN fracturing
//					information
// 
//  Documents:		FracResUnit.h
//
//  Creation date  : 01/11/2022
//
//  Author:			Pascal Siegel - Rabah Namar
//
//  Copyright (C) 2022 by in2earth
//
//---------------------------- FracResUnit.h  ---------------------------
#pragma once

#ifndef FracResUnit_h
#define FracResUnit_h

// local includes
#include <remfracres/FracResTypeEnum.h>

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
 *    FracResUnit.h contains the parameters for the build
 *    hypothesis fracturing
 * 
 *    @author Pascal Siegel - Rabah Namar - Olivier Jaquet 
 *    @date 2023
 */
 
 class FracResUnit
{

// variables

public:

	 std::string name_;
	 std::string unit_age_;
	 std::string unit_duration_;
	 std::string unit_distance_;
	 std::string unit_altitude_;
	 std::string unit_production_duration_;
	 std::string unit_erosion_duration_;
	 std::string unit_subsidence_duration_;
	 std::string unit_river_sediment_duration_;
	 std::string unit_river_sediment_volume_;
	 std::string unit_temperature_;
	 std::string unit_slope_;
	
// functions

public: 

	/**
	 *      Default constructor 
	 */

     FracResUnit();

	/**
	 *		Copy constructor
	 * 
	 *      @param from The value to copy to this object
	 */

     FracResUnit(const FracResUnit& from);

	/**
	 *		Destructor 
	 */
        
	~FracResUnit(){};

	/** 
     *      Assignment operator 
	 *
     *      @param from the value to assign to this object 
     * 
     *      @return A reference to this object 
     */
  
	FracResUnit& operator=(const FracResUnit& from);

	/** 
     *      copy function 
	 *
     *      @param from the value to assign to this object 
     */
  
	void copy(const FracResUnit& from);
	
	
};
// inline methods


 }
#endif  // FracResUnit_h
