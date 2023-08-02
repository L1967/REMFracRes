//----------------------------  FracResGeometryEnum.h  ---------------------------
// 
//
//  Description:	FracResGeometryEnum class use for enulmerate Fractures categories
//					information
// 
//  Documents:		FracResGeometryEnum.h
//
//  Creation date  : 01/11/2022
//
//  Author:			Pascal Siegel - Rabah Namar
//
//  Copyright (C) 2022 by in2earth
//
//---------------------------- FracResGeometryEnum.h  ---------------------------
#pragma once

#ifndef FracResGeometryEnum_h
#define FracResGeometryEnum_h

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
 *    FracResGeometryEnum.h contains the enum parameters for fracture dfn category
 *
 * 
 *    @author Pascal Siegel - Rabah Namar - Olivier Jaquet 
 *    @date 2023
 */
 
 enum  FracResGeometryEnum
{
   ANGLE_WITH_FAUL_ARRAY = 0,
   AZIMUTH = 1,
   DIP_AZIMUTH = 2,
   DIP_ANGLE = 3,
   LENGTH = 4,
   HEIGHT = 5,
   APERTURE = 6,
   WIDTH = 7,
   CORE_ZONE_DIP_AZIMUTH = 8,
   CORE_ZONE_DIP_ANGLE = 9,
   CORE_ZONE_LENGTH = 10,
   CORE_ZONE_HEIGHT = 11,
   CORE_ZONE_APERTURE = 12,
   DAMAGE_ZONE_WIDTH = 13,
   HANGING_WALL_WIDTH = 14,
   FOOT_WALL_WIDTH = 15,
   ANGLE_ARRAY_DIRECTION = 16,
   MAX_DISTANCE_CONNECTIVITY = 17


	
	
};
// inline methods


 }
#endif  // FracResGeometryEnum
