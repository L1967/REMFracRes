//----------------------------  model.h  ---------------------------
// 
//
//  Description:	Base class for all the model objects. 
//                  mesh, variables, boundary , material data , numeric solver 
//
//  Documents:		model.h
//
//  Creation date  : 01/2006
//      
//
//  Author:			 Pascal Siegel
//
//  Copyright (C) 2006 by Colenco Power Engineering Ltd.
//
//---------------------------- model.h  ---------------------------

#ifndef model_h
#define model_h


// system includes
#include <vector>
#include <map>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <primitive_mesh/contrainte_geo.h>

namespace fractures_intersect {

// local includes 



/**
 *    model description : Base class for all the model objects.
 *						is a template.
 *						 mesh, variables, boundary , material data , numeric solver 
 *   
 * 
 *    @author Siegel pascal And Rabah Namar
 *    @date 06/2010
 */


class model
{

// variables


public:


	/**
	 *	geometric constraint 	 
	 */
	Ccontrainte_Geo geo_cont_; 



private:


// functions


public : 


	/**
	 *      Default constructor 
	 */
    model();



	/**
	 *	Destructor 
	 */       
	~model()
	{
	}; 


	/**
	 *	read the fractures data   
	 */

	bool  read_fractures(std::string filename);

	/**
	 *	read the fractures data   
	 */

	bool  save_meshed_plan(std::string filename);





};

	bool save_meshed_vtk_plan(std::string filename, std::vector<Point3D>& tab_point, std::vector<std::vector<int>>& tab_triangle,
			std::vector<property<float>>& tab_triangle_face_property, bool is_binary = false);

	template <typename T>
	void SwapEnd(T& var){
	    char* varArray = reinterpret_cast<char*>(&var);
	    for(long i = 0; i < static_cast<long>(sizeof(var)/2); i++)
	        std::swap(varArray[sizeof(var) - 1 - i],varArray[i]);
	}

}







#endif  // model_h
