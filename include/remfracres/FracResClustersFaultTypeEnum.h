//----------------------------  FracResClustersFaultTypeEnum.h  ---------------------------
// 
//
//  Description:	FracResClustersFaultTypeEnum class use for enulmerate Fractures categories
//					information
// 
//  Documents:		FracResClustersFaultTypeEnum.h
//
//  Creation date  : 01/11/2022
//
//  Author:			Pascal Siegel - Rabah Namar
//
//  Copyright (C) 2022 by in2earth
//
//---------------------------- FracResClustersFaultTypeEnum.h  ---------------------------
#pragma once

#ifndef FracResClustersFaultTypeEnum_h
#define FracResClustersFaultTypeEnum_h

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
 *    FracResClustersFaultTypeEnum.h contains the enum parameters for fracture dfn category
 *
 * 
 *    @author Pascal Siegel - Rabah Namar - Olivier Jaquet 
 *    @date 2023
 */
 
 enum FracResClustersFaultTypeEnum
{
   Null = 0,
   DIP_SLIP_NORMAL = 1,
   DIP_SLIP_REVERSE =2,
   STRIKE_SLIP = 3

};
// inline methods


 }
#endif  // FracResClustersFaultTypeEnum
