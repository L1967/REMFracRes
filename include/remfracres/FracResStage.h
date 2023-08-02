//----------------------------  FracResStage.h  ---------------------------
// 
//
//  Description:	FracResStage class use for build array fault
//					information
// 
//  Documents:		FracResStage.h
//
//  Creation date  : 01/11/2022
//
//  Author:			Pascal Siegel - Rabah Namar
//
//  Copyright (C) 2022 by in2earth
//
//---------------------------- FracResStage.h  ---------------------------
#pragma once

#ifndef FracResStage_h
#define FracResStage_h

// local includes
#include <remfracres/FracResHypothesis.h>
#include <remfracres/FracResFracture.h>

// standard includes
#include <fstream>
#include <iostream>
#include <iomanip>
#include <memory>
#include <vector>
#include <string>




namespace remfracres{

/**
 *    FracResStage contains the parameters for the simulation
 *    stage fracturing
 * 
 *    @author Pascal Siegel - Rabah Namar - Olivier Jaquet 
 *    @date 2022
 */
 
 class FracResStage
{

// variables

public:

	 std::string name_;

	 int stageIndex_;

	 bool beddingStage_;

     std::vector<FracResHypothesis> hypotheisisList_;
	 std::map< int, FracResFracture > fractureList_;
	 std::vector<int>	evaluate_fractures_number_total_per_hypothesis_ ;
	
// functions

public: 

	/**
	 *      Default constructor 
	 */

     FracResStage();

	/**
	 *		Copy constructor
	 * 
	 *      @param from The value to copy to this object
	 */

     FracResStage(const FracResStage& from);

	/**
	 *		Destructor 
	 */
        
	~FracResStage(){};

	/** 
     *      Assignment operator 
	 *
     *      @param from the value to assign to this object 
     * 
     *      @return A reference to this object 
     */
  
	FracResStage& operator=(const FracResStage& from);

	/** 
     *      copy function 
	 *
     *      @param from the value to assign to this object 
     */
  
	void copy(const FracResStage& from);
	
	
};
// inline methods


 }
#endif  // FracResStage_h
