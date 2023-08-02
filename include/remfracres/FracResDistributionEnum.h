//----------------------------  FracResDistributionEnum.h  ---------------------------
// 
//
//  Description:	FracResDistributionEnum class use for enulmerate Fractures categories
//					information
// 
//  Documents:		FracResDistributionEnum.h
//
//  Creation date  : 01/11/2022
//
//  Author:			Pascal Siegel - Rabah Namar
//
//  Copyright (C) 2022 by in2earth
//
//---------------------------- FracResDistributionEnum.h  ---------------------------
#pragma once

#ifndef FracResDistributionEnum_h
#define FracResDistributionEnum_h

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
 *    FracResDistributionEnum.h contains the enum parameters for fracture dfn category
 *
 * 
 *    @author Pascal Siegel - Rabah Namar - Olivier Jaquet 
 *    @date 2023
 */
 
 enum  FracResDistributionEnum
{
   FRACTURE_CORRIDORS_DENSITY = 0,
   FAULT_ZONE_DENSITY = 1,
   FRACTURE_DENSITY = 2,
   FAULT_DENSITY = 3,
   MINIMUM_SPACE_BETWEEN_ARRAYS = 4,
   ARRAY_DENSITY = 5,
   FAULT_SEGMENT_DENSITY = 6,
   DISTANCE_BETWEEN_SEG = 7,
   WIDTH_STRESS_SHADOW_FAMILY1 = 8,
   WIDTH_STRESS_SHADOW_FAMILY2 = 9,
   FRACTURE_DENSITY_CORRIDORS = 10,
   WIDTH_STRESS_SHADOW_ZONE = 11,
   MINIMUM_SPACING_BETWEEN_FAULT_ZONES = 12,
   MINIMUM_SPACING_BETWEEN_FRACTURE_CORRIDORS = 13


	
	
};
// inline methods


 }
#endif  // FracResDistributionEnum
