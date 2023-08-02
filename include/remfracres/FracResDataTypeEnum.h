//----------------------------  FracResDataTypeEnum.h  ---------------------------
// 
//
//  Description:	FracResDataTypeEnum class use for enulmerate Fractures categories
//					information
// 
//  Documents:		FracResDataTypeEnum.h
//
//  Creation date  : 01/11/2022
//
//  Author:			Pascal Siegel - Rabah Namar
//
//  Copyright (C) 2022 by in2earth
//
//---------------------------- FracResDataTypeEnum.h  ---------------------------
#pragma once

#ifndef FracResDataTypeEnum_h
#define FracResDataTypeEnum_h

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
 *    FracResDataTypeEnum.h.h contains the enum parameters for fracture dfn category
 *
 * 
 *    @author Pascal Siegel - Rabah Namar - Olivier Jaquet 
 *    @date 2023
 */
 
 enum  FracResDataTypeEnum
{
   CONSTANT = 0,
   DISTRIBUTION = 1,
   PROPERTY = 2,
   NORMAL_DISTRIBUTION = 3

	
	
};
// inline methods


 }
#endif  // FracResDataTypeEnum.h
