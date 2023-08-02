//----------------------------  FracResGeometryCorrelationTypeEnum.h  ---------------------------
// 
//
//  Description:	FracResGeometryCorrelationTypeEnum class use for enumerate Fractures categories
//					information
// 
//  Documents:		FracResGeometryCorrelationTypeEnum.h
//
//  Creation date  : 01/11/2022
//
//  Author:			Pascal Siegel - Rabah Namar
//
//  Copyright (C) 2022 by in2earth
//
//---------------------------- FracResGeometryCorrelationTypeEnum.h  ---------------------------
#pragma once

#ifndef FracResGeometryCorrelationTypeEnum_h
#define FracResGeometryCorrelationTypeEnum_h

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
 *    FracResGeometryCorrelationTypeEnum.h.h contains the enum parameters for fracture dfn category
 *
 * 
 *    @author Pascal Siegel - Rabah Namar - Olivier Jaquet 
 *    @date 2023
 */
 
 enum  FracResGeometryCorrelationTypeEnum
{
   INDEPENDANT= 0,
   CORRELATED=1

	
	
};
// inline methods


 }
#endif  // FracResGeometryCorrelationTypeEnum.h
