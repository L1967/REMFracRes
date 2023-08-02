//----------------------------  FracResDistributionArrayTypeEnum.h  ---------------------------
// 
//
//  Description:	FracResDistributionArrayTypeEnum class use for enulmerate Fractures categories
//					information
// 
//  Documents:		FracResDistributionArrayTypeEnum.h
//
//  Creation date  : 01/11/2022
//
//  Author:			Pascal Siegel - Rabah Namar
//
//  Copyright (C) 2022 by in2earth
//
//---------------------------- FracResDistributionArrayTypeEnum.h  ---------------------------
#pragma once

#ifndef FracResDistributionArrayTypeEnum_h
#define FracResDistributionArrayTypeEnum_h

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
 *    FracResDistributionArrayTypeEnum.h contains the enum parameters for fracture dfn category
 *
 * 
 *    @author Pascal Siegel - Rabah Namar - Olivier Jaquet 
 *    @date 2023
 */
 
 enum  FracResDistributionArrayTypeEnum
{
   SINGLE_FAULT = 0,
   FAULT_ARRAY_ARRAY = 1,
   FAULT_ARRAY_SEGMENT = 2,
   FRACTURE_CLUSTERS_ZONE_CLUSTERS = 3,
   FRACTURE_CLUSTERS_ZONE_FRACTURES_HANGING = 4,
   FRACTURE_CLUSTERS_ZONE_FRACTURES_FOOTWALL = 5,
   FRACTURE_CLUSTERS_ZONE_FRACTURES = 6,
   FRACTURE_CLUSTERS_CORRIDOR_FRACTURES_CORRIDORS = 7,
   FRACTURE_CLUSTERS_CORRIDOR_FRACTURES = 8
	
};
// inline methods


 }
#endif  // FracResDistributionArrayTypeEnum
