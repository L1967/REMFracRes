//----------------------------  FracResTypeEnum.h  ---------------------------
// 
//
//  Description:	FracResTypeEnum class use for build hypothesis fracturing
//					information
// 
//  Documents:		FracResTypeEnum.h
//
//  Creation date  : 01/11/2022
//
//  Author:			Pascal Siegel - Rabah Namar
//
//  Copyright (C) 2022 by in2earth
//
//---------------------------- FracResTypeEnum.h  ---------------------------
#pragma once

#ifndef FracResTypeEnum_h
#define FracResTypeEnum_h

// local includes
#include <remfracres/FracResCategoryEnum.h>

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
 *    FracResTypeEnum.h contains the enum parameters for fracture dfn type
 *
 * 
 *    @author Pascal Siegel - Rabah Namar - Olivier Jaquet 
 *    @date 2023
 */
 
 enum  FracResTypeEnum
{
   BED_PARALLEL_SLIP = 0,
   SYSTEMATIC_JOINTS = 1,
   CONJUGATE_JOINTS = 2,
   SINGLE_FAULTS = 3,
   CONJUGATE_FAULTS = 4,
   ANASTOMOSING_FAULTS = 5,
   RELAY_FAULTS = 6,
   EN_ECHELON_FAULTS = 7,
   JOINTS_SWARM = 8,
   FRACTURE_CORRIDORS = 9,
   FAULT_ZONES = 10

	
	
};
// inline methods


 }
#endif  // FracResTypeEnum
