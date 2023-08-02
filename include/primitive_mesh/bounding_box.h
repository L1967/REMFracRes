//----------------------------  bounding_box.h  ---------------------------
// 
//
//  Description :		This class is for the bounding box of the Primitives 
//						and other objects (parts of a Mesh,...)
//						A bounding_box is constituted by the RealDomain on which it extends.
//  
//
//  Documents:			bounding_box.h
//						bounding_box.cpp
//
//
//  Creation date  :	10/01/20013
//
//  Author:				Pascal Siegel
//
//  Copyright (C) 2013 by in2earth
//
//---------------------------- bounding_box.h  ---------------------------

#ifndef bounding_box_h
#define bounding_box_h


// local includes 
#include <limits>
#include <primitive_mesh/point.h>


// project includes 


// system includes


// forward references 


namespace fractures_intersect {


/**
 *    bounding_box description : This class is for the bounding box of the Primitives
 *								and other objects (parts of a Mesh,...)
 *								A bounding_box is constituted by the RealDomain on which it extends.
 *   
 */

class bounding_box {


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

	bool  bounding_set_; // 1 si on dejas calcul� la bounding box


// functions


public : 


	/**
	 *      Default constructor 
	 */
	bounding_box() {bounding_set_ = false;


	pt_min_ = {std::numeric_limits<double>::max(), std::numeric_limits< double >::max(),std::numeric_limits< double >::max()};
	pt_max_ = {std::numeric_limits<double>::lowest(), std::numeric_limits< double >::lowest(),std::numeric_limits< double >::lowest()};

	};


	/**
	 *      Default constructor 
	 */
    bounding_box(Point3D& pt_min, Point3D& pt_max)
	{
		pt_min_ = pt_min;
		pt_max_ = pt_max;

		bounding_set_ = true;
	}


	/**
	 *	Destructor 
	 */        
	~bounding_box() {};

	/** 
     *      add point
	 *
	 *		@param new Point
	 *
	 *		@return a bool true if the tested point is inside the bounding_box
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

	void set_point(Point3D& pt)
	{
		pt_min_.coord_ = pt.coord_;
		pt_max_= pt.coord_;
	}

	double length() {
		return (pt_max_.coord_- -pt_min_.coord_).norm();
	}

	void add_bounding_box(const bounding_box &box)
	{
		pt_min_.x() = std::min(pt_min_.x(),box.pt_min_.x());
		pt_min_.y() = std::min(pt_min_.y(),box.pt_min_.y());
		pt_min_.z() = std::min(pt_min_.z(),box.pt_min_.z());

		pt_max_.x() = std::max(pt_max_.x(),box.pt_max_.x());
		pt_max_.y() = std::max(pt_max_.y(),box.pt_max_.y());
		pt_max_.z() = std::max(pt_max_.z(),box.pt_max_.z());
	}

	/** 
     *	Extends the bounding box : adds "margins" of a width epsilon. 
     *  This is useful to avoid some truncature / computing errors.
     */

	void extend(double epsilon)
	{
		pt_min_.x() -= epsilon ;
		pt_min_.y() -= epsilon ;
		pt_min_.z() -= epsilon ;

		pt_max_.x() += epsilon ;
		pt_max_.y() += epsilon ;
		pt_max_.z() += epsilon ;
		
	}

	void extend(Eigen::Matrix<double, 3, 1>& vec)
	{
		pt_min_.coord_ = pt_min_.coord_ - vec ;
		pt_max_.coord_ = pt_max_.coord_ + vec ;
	}

	bool contains(Point3D& pt) {
		if (pt.x()<pt_min_.x() || pt.x()>pt_max_.x()) return false;
		if (pt.y()<pt_min_.y() || pt.y()>pt_max_.y()) return false;
		if (pt.z()<pt_min_.z() || pt.z()>pt_max_.z()) return false;

		return true;
	}


	bool contains(const bounding_box &box) {
		return (box.pt_max_.x() < pt_max_.x() && box.pt_min_.x() > pt_min_.x()   &&
				box.pt_max_.y() < pt_max_.y() && box.pt_min_.y() > pt_min_.y()   &&
				box.pt_max_.z() < pt_max_.z() && box.pt_min_.z() > pt_min_.z());
	}

	bool iscontained(const bounding_box &box) {
		return (pt_max_.x() < box.pt_max_.x() && pt_min_.x() > box.pt_min_.x()   &&
				pt_max_.y() < box.pt_max_.y() && pt_min_.y() > box.pt_min_.y()   &&
				pt_max_.z() < box.pt_max_.z() && pt_min_.z() > box.pt_min_.z());
	}

	
	bool intersection(const bounding_box &box)
	{

		if (pt_max_.x() < box.pt_min_.x() || pt_min_.x() > box.pt_max_.x()) return false;
		if (pt_max_.y() < box.pt_min_.y() || pt_min_.y() > box.pt_max_.y()) return false;
		if (pt_max_.z() < box.pt_min_.z() || pt_min_.z() > box.pt_max_.z()) return false;
		// Overlapping on all axes means AABBs are intersecting
		return true;

	}


	bool intersectionAndContains(const bounding_box &box)
	{
		if (intersection(box)) return true;

		if (contains(box) ||  iscontained(box)) return true ;

		return false;

	}

	Vector3D diagonal()
	{
	    return (pt_max_.coord_ - pt_min_.coord_);
	}

	Point3D center() const
	{
	    return Point3D((pt_max_.coord_ + pt_min_.coord_) / 2.);
	}

	double signed_distance(Point3D& pt) {

		bool inside{ true };
		Vector3D result;

		for(int i = 0; i < 3; i++ )
		{
			double value = pt[i];

			if( value < pt_min_[i] )
			{
				inside = false;
				result[i]  =  value - pt_min_[i];
			}
			else if( value > pt_max_[i] )
			{
				inside = false;
				result[i]  =  value - pt_max_[i];
			}

		}

		if( !inside )
		{
			return result.norm();
		}
		else {
			result = pt.coord_ -center().coord_;
			return result.norm();
		}
	}
	
	
};

}




#endif  // bounding_box_h
