//----------------------------  FracResScenario.h  ---------------------------
// 
//
//  Description:	FracResScenario class use for build Scenario fracturing
//					information
// 
//  Documents:		FracResScenario.h
//
//  Creation date  : 01/11/2022
//
//  Author:			Pascal Siegel - Rabah Namar
//
//  Copyright (C) 2022 by in2earth
//
//---------------------------- FracResScenario.h  ---------------------------
#pragma once

#ifndef FracResScenario_h
#define FracResScenario_h

// local includes
#include <remfracres/FracResStage.h>
#include <remfracres/FracResHypothesis.h>



// standard includes
#include <fstream>
#include <iostream>
#include <iomanip>
#include <memory>
#include <vector>
#include <string>

#include <numeric/model.h>

namespace remfracres{

/**
 *    FracResScenario.h contains the parameters for the build
 *    Scenario fracturing
 * 
 *    @author Pascal Siegel - Rabah Namar - Olivier Jaquet 
 *    @date 2023
 */
 
 class FracResScenario
{

// variables

public:

	std::string scenario_file_name_;
	std::string name_;
	int scenarioIndex_;
	int estimate_Fracture_Number_ ;
	bool modeling_ ;
	int numberOfRealisation_ ;
	std::vector< std::pair<int, int> > scenario_stage_hypothesis_index_vector_;
	std::map<int, std::vector<  bool > > success_list_;
	std::string path_result_directory_;
	std::vector< fractures_intersect::model > model_fractures_list_;
	int nb_failles_;
	std::vector< double > total_simulation_time_list_;
	std::vector< double > total_writting_output_time_list_;
	std::vector< int > evaluate_fractures_number_list_;


// functions

public: 

	/**
	 *      Default constructor 
	 */

     FracResScenario();


 	/**
 	 *		Destructor
 	 */

 	~FracResScenario(){};

	/**
	 *		Copy constructor
	 * 
	 *      @param from The value to copy to this object
	 */

     FracResScenario(const FracResScenario& from);


	/** 
     *      Assignment operator 
	 *
     *      @param from the value to assign to this object 
     * 
     *      @return A reference to this object 
     */
  
	FracResScenario& operator =(const FracResScenario& from);

	/** 
     *      copy function 
	 *
     *      @param from the value to assign to this object 
     */
  
	void copy(const FracResScenario& from);
	

};
// inline methods


 }
#endif  // FracResScenario_h
