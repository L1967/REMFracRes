//----------------------------  FracResBeddingLocalisation.h  ---------------------------
// 
//
//  Description:	FracResBeddingLocalisation class use for bedding localisation
//					information
// 
//  Documents:		FracResBeddingLocalisation.h
//
//  Creation date  : 01/06/2023
//
//  Author:			Pascal Siegel - Rabah Namar
//
//  Copyright (C) 2023 by in2earth
//
//---------------------------- FracResBeddingLocalisation.h  ---------------------------
#pragma once

#ifndef FracResBeddingLocalisation_h
#define FracResBeddingLocalisation_h



// standard includes
#include <fstream>
#include <iostream>
#include <iomanip>
#include <memory>
#include <vector>
#include <string>
// geode includes

#include <geode/basic/uuid.h>


namespace fractures_intersect{

/**
 *    FracResBeddingLocalisation.h contains the parameters for the bedding localisation
 *
 * 
 *    @author Pascal Siegel - Rabah Namar - Olivier Jaquet 
 *    @date 2023
 */
 
 class FracResBeddingLocalisation
{

// variables

public:


	 int index_;

	 int index_polyedre_1_;

	 int index_polyedre_2_;

	 geode::uuid  block_id_;
	
// functions

public: 

	/**
	 *      Default constructor 
	 */

	 FracResBeddingLocalisation();


	 FracResBeddingLocalisation( int index, int polyhedron_index1, int polyhedron_index2, geode::uuid id_block );

	/**
	 *		Copy constructor
	 * 
	 *      @param from The value to copy to this object
	 */

	 FracResBeddingLocalisation(const FracResBeddingLocalisation& from);

	/**
	 *		Destructor 
	 */
        
	~FracResBeddingLocalisation(){};

	/** 
     *      Assignment operator 
	 *
     *      @param from the value to assign to this object 
     * 
     *      @return A reference to this object 
     */
  
	FracResBeddingLocalisation& operator=(const FracResBeddingLocalisation& from);

	/** 
     *      copy function 
	 *
     *      @param from the value to assign to this object 
     */
  
	void copy(const FracResBeddingLocalisation& from);
	
	
};
// inline methods


 }
#endif  // FracResBeddingLocalisation_h
