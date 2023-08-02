//----------------------------  FracResCategoryEnum.h  ---------------------------
// 
//
//  Description:	FracResCategoryEnum class use for enulmerate Fractures categories
//					information
// 
//  Documents:		FracResCategoryEnum.h
//
//  Creation date  : 01/11/2022
//
//  Author:			Pascal Siegel - Rabah Namar
//
//  Copyright (C) 2022 by in2earth
//
//---------------------------- FracResCategoryEnum.h  ---------------------------
#pragma once

#ifndef FracResCategoryEnum_h
#define FracResCategoryEnum_h

// local includes


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
 *    FracResCategoryEnum.h contains the enum parameters for fracture dfn category
 *
 * 
 *    @author Pascal Siegel - Rabah Namar - Olivier Jaquet 
 *    @date 2023
 */
 
 enum  FracResCategoryEnum
{
   BEDDINGS = 0,
   JOINTS = 1,
   FAULTS = 2,
   FAULTS_ARRAY = 3,
   FRACTURE_CLUSTER = 4

	
	
};
// inline methods


 }
#endif  // FracResCategoryEnum
