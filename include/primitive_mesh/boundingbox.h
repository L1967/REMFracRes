//----------------------------  BoundingBox.h  ---------------------------
// 
//
//  Description :		This class is for the bounding box of the Primitives 
//						and other objects (parts of a Mesh,...)
//						A BoundingBox is constituted by the RealDomain on which it extends.
//  
//
//  Documents:			boundingbox.h
//						boundingbox.cpp
//
//
//  Creation date  :	10/01/20013
//
//  Author:				Pascal Siegel
//
//  Copyright (C) 2013 by in2earth
//
//---------------------------- BoundingBox.h  ---------------------------

#ifndef BoundingBox_h
#define BoundingBox_h


// local includes 
#include <primitive_mesh/point.h>


// project includes 


// system includes


// forward references 


namespace fractures_intersect {


/**
 *    BoundingBox description : This class is for the bounding box of the Primitives
 *								and other objects (parts of a Mesh,...)
 *								A BoundingBox is constituted by the RealDomain on which it extends.
 *   
 * 
 *    @author Nicolas Hubschwerlen
 *    @date 05/2003
 */

class BoundingBox {


	// variables

public:

	/**
	 *		point min
	 */
	Point3D pt_min_;

	/**
	 *		point max 
	 */
	Point3D pt_max_;


// functions


public : 


	/**
	 *      Default constructor 
	 */
	BoundingBox() {};


	/**
	 *      Default constructor 
	 */
    BoundingBox(Point3D& pt1, Point3D& pt2)
	{
		pt_min_.x() = std::min(pt1.x(),pt2.x());
		pt_min_.y() = std::min(pt1.y(),pt2.y()); 
		pt_min_.z() = std::min(pt1.z(),pt2.z()); 

		pt_max_.x() = std::max(pt1.x(),pt2.x());
		pt_max_.y() = std::max(pt1.y(),pt2.y()); 
		pt_max_.z() = std::max(pt1.z(),pt2.z()); 
	}


	/**
	 *	Destructor 
	 */        
	~BoundingBox() {};

	/** 
     *      add point
	 *
	 *		@param new Point
	 *
	 *		@return a bool true if the tested point is inside the BoundingBox
     */
	void add_point(Point3D& pt)
	{
		pt_min_.x() = std::min(pt_min_.x(),pt.x());
		pt_min_.y() = std::min(pt_min_.y(),pt.y()); 
		pt_min_.z() = std::min(pt_min_.z(),pt.z()); 

		pt_max_.x() = std::max(pt_max_.x(),pt.x());
		pt_max_.y() = std::max(pt_max_.y(),pt.y()); 
		pt_max_.z() = std::max(pt_max_.z(),pt.z()); 
	}

	/** 
     *	Extends the bounding box : adds "margins" of a width epsilon. 
     *  This is useful to avoid some truncature / computing errors.
     */

	void extendBoundingBox(double epsilon)
	{
		pt_min_.x() -= epsilon ;
		pt_min_.y() -= epsilon ;
		pt_min_.z() -= epsilon ;

		pt_max_.x() += epsilon ;
		pt_max_.y() += epsilon ;
		pt_max_.z() += epsilon ;
		
	}

	
	bool intersection(const BoundingBox &box)
	{
		Point3D min(
			std::max(pt_min_.x(), box.pt_min_.x()),
			std::max(pt_min_.y(), box.pt_min_.y()),
			std::max(pt_min_.z(), box.pt_min_.z())
		 );

		Point3D max(
			std::min(pt_max_.x(), box.pt_max_.x()),
			std::min(pt_max_.y(), box.pt_max_.y()),
			std::min(pt_max_.z(), box.pt_max_.z())
		 );

		 if( max.x() < min.x() || max.y() < min.y() || max.z() < min.z() )
		 {
			return false;
		 } 
		 else return true;
		 
	}


	
	
};

}




#endif  // BoundingBox_h
