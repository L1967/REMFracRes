//----------------------------  FracResHypothesis.h  ---------------------------
// 
//
//  Description:	FracResHypothesis class use for build hypothesis fracturing
//					information
// 
//  Documents:		FracResHypothesis.h
//
//  Creation date  : 01/11/2022
//
//  Author:			Pascal Siegel - Rabah Namar
//
//  Copyright (C) 2022 by in2earth
//
//---------------------------- FracResHypothesis.h  ---------------------------
#pragma once

#ifndef FracResHypothesis_h
#define FracResHypothesis_h

// local includes
#include <remfracres/FracResBeddingProperties.h>
#include <remfracres/FaultArraySet.h>
#include <remfracres/fractureset.h>



// standard includes
#include <fstream>
#include <iostream>
#include <iomanip>
#include <memory>
#include <vector>
#include <string>



namespace remfracres{

/**
 *    FracResHypothesis.h contains the parameters for the build
 *    hypothesis fracturing
 * 
 *    @author Pascal Siegel - Rabah Namar - Olivier Jaquet 
 *    @date 2023
 */
 
 class FracResHypothesis
{

// variables

public:

	 std::string stage_name_;

	 std::string name_;

	 int hypothesisIndex_;

	 bool beddingIsRegion_;

	 std::string beddingRegionName_;

	 std::string beddingPropertyDiscreteName_;

	 std::vector< FaultArraySet* > faultArrayList_;

	 std::vector< FaultArraySet* > faultCorridorList_;

	 std::vector< FaultArraySet* > faultZoneList_;

	 std::vector< FractureSet* > fractureList_;

	 std::vector< FracResBeddingAperture* > beddingList_;

	 std::vector< FracResBeddingProperties* > bed_parallel_List_;
	
	 std::vector<int> type_dfn_;

	 std::vector<int> index_dfn_;

	 int evaluate_fractures_number_total_;

	 int nb_bed_parallel_active_;
	 int nb_bed_interface_active_;
	 int nb_single_fault_active_;
	 int nb_fault_array_active_;
	 int nb_cluster_corridors_active_;
	 int nb_cluster_fault_zone_active_;


// functions

public: 

	/**
	 *      Default constructor 
	 */

     FracResHypothesis();


 	/**
 	 *		Destructor
 	 */

 	~FracResHypothesis(){};

	/**
	 *		Copy constructor
	 * 
	 *      @param from The value to copy to this object
	 */

     FracResHypothesis(const FracResHypothesis& from);


	/** 
     *      Assignment operator 
	 *
     *      @param from the value to assign to this object 
     * 
     *      @return A reference to this object 
     */
  
	FracResHypothesis& operator =(const FracResHypothesis& from);

	/** 
     *      copy function 
	 *
     *      @param from the value to assign to this object 
     */
  
	void copy(const FracResHypothesis& from);
	


	void update_active_dfn_number();


};
// inline methods


 }
#endif  // FracResHypothesis_h
