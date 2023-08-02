//----------------------------  oriented_bounding_box.h  ---------------------------
// 
//
//  Description :		This class is for the bounding box of the Primitives 
//						and other objects (parts of a Mesh,...)
//						A oriented_bounding_box is constituted by the RealDomain on which it extends.
//  
//
//  Documents:			oriented_bounding_box.h
//						oriented_bounding_box.cpp
//
//
//  Creation date  :	10/01/20013
//
//  Author:				Pascal Siegel
//
//  Copyright (C) 2013 by in2earth
//
//---------------------------- oriented_bounding_box.h  ---------------------------

#ifndef oriented_bounding_box_h
#define oriented_bounding_box_h


// local includes 
#include <primitive_mesh/point.h>
#include <primitive_mesh/faille.h>


// project includes 


// system includes


// forward references 


namespace fractures_intersect {


/**
 *    oriented_bounding_box description : This class is for the oriented bounding box of the Primitives
 *								and other objects (parts of a Mesh,...)
 *								A oriented_bounding_box is constituted by the RealDomain on which it extends.
 *   
 */

class oriented_bounding_box {


	// variables

public:

	Vector3D	R_,S_,T_;	// axis unit vector : T  normal vector
	Point3D	center_;	// axis center


	Eigen::Matrix<double, 3, 3> rotation_matrix_;
	Eigen::Matrix<double, 3, 3> rotation_matrix_tr_;

	/**
	 *		extension of box in each axis
	 */
	Vector3D extension_;

	bool  bounding_set_; // ture if the box is generated

	double  volume_; // volume of bounding box


	double  area_; // area
	// functions


public : 


	/**
	 *      Default constructor 
	 */
	oriented_bounding_box() {

		R_ = {1.0, 0.0,0.0};
		S_ = {0.0, 1.0,0.0};
		T_ = {0.0, 0.0,1.0};

		extension_ = {1.0, 1.0,1.0};

		bounding_set_ = false;

	};


	/**
	 *      Default constructor 
	 */
	oriented_bounding_box(Point3D& center, Vector3D&	R, Vector3D& S,Vector3D& T, Vector3D& extension)
	{
		R_ = R;
		S_ = S;
		T_ = T;

		center_ = center;

		extension_ = extension;

		bounding_set_ = true;
	};

	/**
	 *      Default constructor
	 */
	oriented_bounding_box(Cfaille& fault)
	{
		R_ = fault.Get_R();
		S_ = fault.Get_S();
		T_ = fault.Get_vec_normal();

		center_ = fault.Get_point_ref();

		extension_ = {fault.ext2L_selon_lpgp, fault.ext2l_orthog_lpgp, fault.aperture_};

		bounding_set_ = true;

		update_rotation_matrix();
	};

	/**
	 *      Default constructor
	 */
	oriented_bounding_box(bounding_box& box)
	{
		R_ = {1.0, 0.0,0.0};
		S_ = {0.0, 1.0,0.0};
		T_ = {0.0, 0.0,1.0};

		extension_ = box.pt_max_.coord_ - box.pt_min_.coord_;

		bounding_set_ = true;

		update_rotation_matrix();
	};


	/**
	 *	Destructor 
	 */        
	~oriented_bounding_box() {};


	void update_rotation_matrix()
	{
		rotation_matrix_tr_.col(0) = R_;
		rotation_matrix_tr_.col(1) = S_;
		rotation_matrix_tr_.col(2) = T_;

		rotation_matrix_ = rotation_matrix_tr_.transpose();
	};


	Point3D local_to_global(Point3D& local) {
		return Point3D(rotation_matrix_tr_*local.coord_ + center_.coord_);
	};

	Point3D global_to_local(Point3D& global) {
		return Point3D(rotation_matrix_*(global.coord_ - center_.coord_));
	};

	Point3D generate_point(double random_x, double random_y, double random_z) {
		random_x = 2*random_x -1.0;
		random_y = 2*random_y -1.0;
		random_z = 2*random_z -1.0;

		Point3D localCoord (0.5*random_x*extension_.x(),  0.5*random_y*extension_.y(), 0.5*random_z*extension_.z());
		return local_to_global(localCoord);
	};
	Point3D generate_point(double random_x, double random_y, double random_z, double zinit, double ratio) {
		random_x = 2*random_x -1.0;
		random_y = 2*random_y -1.0;
		random_z = (random_z/ratio) -1.0;

		Point3D localCoord (0.5*random_x*extension_.x(),  0.5*random_y*extension_.y(), zinit+ ratio*random_z*extension_.z());
		return local_to_global(localCoord);
	};

	/** 
     *	Extends the bounding box : adds "margins" of a width epsilon. 
     *  This is useful to avoid some truncature / computing errors.
     */

	void extend(double epsilon)
	{
		 extension_[0] += 2.0 * epsilon ;
		 extension_[1] += 2.0 * epsilon ;
		 extension_[2] += 2.0 * epsilon ;
	};

	void extend(Eigen::Matrix<double, 3, 1>& vec)
	{
		extension_ += 2.0*vec;
	};

	bool contains(Point3D& pt) {

		Point3D local = global_to_local(pt);

		if (local.x()< extension_.x()/(-2.0) || local.x()> extension_.x()/2.0) return false;
		if (local.y()< extension_.y()/(-2.0) || local.y()> extension_.y()/2.0) return false;
		if (local.z()< extension_.z()/(-2.0) || local.z()> extension_.z()/2.0) return false;

		return true;
	}

	double surface()
	{
		return  extension_.x() *  extension_.y();
	}

	
	bool intersection(const oriented_bounding_box &box)
	{
		return true;
	};


	bool intersectionAndContains(const oriented_bounding_box &box)
	{
		return false;
	};

	double volume(){

		volume_ = extension_.x()*extension_.y()*extension_.z();
		return volume_;
	};

	double area(){

		area_ = extension_.x()*extension_.y();
		return area_;
	};
	double diagonal()
	{
		return  extension_.norm();
	};
};

}




#endif  // oriented_bounding_box_h
