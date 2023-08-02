//----------------------------  FracResFracture.h  ---------------------------
// 
//
//  Description:	FracResFracture class use for build hypothesis fracturing
//					information
// 
//  Documents:		FracResFracture.h
//
//  Creation date  : 01/11/2022
//
//  Author:			Pascal Siegel - Rabah Namar
//
//  Copyright (C) 2022 by in2earth
//
//---------------------------- FracResFracture.h  ---------------------------
#pragma once

#ifndef FracResFracture_h
#define FracResFracture_h

// local includes
#include <remfracres/FracResTypeEnum.h>
#include <remfracres/FracResCategoryEnum.h>

// standard includes
#include <fstream>
#include <iostream>
#include <iomanip>
#include <memory>
#include <vector>
#include <string>



namespace remfracres{

/**
 *    FracResFracture.h contains the parameters for the build
 *    hypothesis fracturing
 * 
 *    @author Pascal Siegel - Rabah Namar - Olivier Jaquet 
 *    @date 2023
 */
 
 class FracResFracture
{

// variables

public:

	 std::string stage_name_;
	 std::string name_;
     int index_Id_;
     FracResTypeEnum dataType_;
     fractures_intersect::FracResIntersectionCategoryEnum categoryType_;
	
	
// functions

public: 

	/**
	 *      Default constructor 
	 */

     FracResFracture();

	/**
	 *		Copy constructor
	 * 
	 *      @param from The value to copy to this object
	 */

     FracResFracture(const FracResFracture& from);

	/**
	 *		Destructor 
	 */
        
	~FracResFracture(){};

	/** 
     *      Assignment operator 
	 *
     *      @param from the value to assign to this object 
     * 
     *      @return A reference to this object 
     */
  
	FracResFracture& operator=(const FracResFracture& from);

	/** 
     *      copy function 
	 *
     *      @param from the value to assign to this object 
     */
  
	void copy(const FracResFracture& from);
	
	
};
// inline methods


 }
#endif  // FracResFracture_h
