//----------------------------  point.h  ---------------------------
// 
//
//  Description: Definition of the class point
//               
//
//  Documents: VPAC modelling repport 
//
//  Creation date  : 01/2013
//       modified  : 
//
//  Author: Siegel Pascal
//
//  Copyright (C) 2013 by in2earth Ltd.
//
//
//---------------------------- point.h  ---------------------------

#ifndef point_h
#define point_h


// system includes

// project includes 

// project includes 
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cmath>

// local includes 

namespace fractures_intersect {

/**
 *    coordinates point  
 *   
 *    @author Pascal Siegel
 *    @date 05/2013
 */

extern double precis,TOL,pi,pt_tol;

static double MIN_THICKNESS = precis;

typedef Eigen::Matrix<double, 3, 1> Vector3D;


class Point3D
{

// variables

public: 

    /**
	 *	vector of coordinates 
     */

    Eigen::Matrix<double, 3, 1>  coord_;


private:


// functions


public : 

	/**
	 *	vector of coordinates 
     */

    bool calculated_;


	/**
	 *      Default constructor 
	 */
    
    Point3D() {coord_[0] = 0.0; coord_[1] = 0.0; coord_[2] = 0.0; calculated_=false;};

    /**
	 *      Default constructor 
	 */
    
    Point3D(double x, double y , double z = 0 ) {
		coord_[0] = x; coord_[1] = y; coord_[2] = z; calculated_ =false;};


    Point3D(Eigen::Matrix<double, 3, 1>  vect) {
        coord_ = vect; calculated_ =false;};

    Point3D(Point3D  pt, bool calc) {
        coord_ = pt.coord_; calculated_ = calc;};

    Point3D(Eigen::Matrix<double, 3, 1>  vect, bool calc) {
        coord_ = vect; calculated_ =calc;};


	/**
	 *	Destructor 
	 */
        
    ~Point3D() {}; 

	/** 
     *	Assignment operator 
	 *
     *  @param pt The value to assign to this object
     * 
     *  @return A reference to this object 
     */
  
    Point3D& operator=(const Point3D& pt)  {coord_ = pt.coord_; calculated_=pt.calculated_; return *this;};
		

    /** 
     *		Get the x coordinate of the point      
	 *
     *      @return an object of template type T
     */

	double  x() const { return coord_[0]; };

	
	/** 
     *		Get the y coordinate of the point     
	 *
     *      @return an object of template type T
     */

	double y() const { return coord_[1]; }; 


    /** 
     *		Get the z coordinate of the point     
	 *
     *      @return an object of template type T
     */

	double z() const { return coord_[2]; };

	/** 
     *		Get the x coordinate of the point      
	 *
     *      @return an object of template type T
     */

	double& x() { return coord_[0]; };

	
	/** 
     *		Get the y coordinate of the point     
	 *
     *      @return an object of template type T
     */

	double& y()  { return coord_[1]; }; 


    /** 
     *		Get the z coordinate of the point     
	 *
     *      @return an object of template type T
     */

	double& z() { return coord_[2]; };


    /** 
     *		Get the vector coordinates of the point     
	 *
     *      @return an vector of template type T
     */

    Eigen::Matrix<double, 3, 1>& vector() { return coord_ ;};

	double& operator[](int i) 
	{
		 return coord_[i]; 
	}

	void arrondi_selon_precis()
	{
		if (std::abs(coord_[0]) < TOL) coord_[0] = 0.;
		if (std::abs(coord_[1]) < TOL) coord_[1] = 0.;
		if (std::abs(coord_[2]) < TOL) coord_[2] = 0.;
	}

	
	/** 
	 *		Equality operator
	 */
	bool operator==(const Point3D& p) const 
	{
		if ( (std::abs(p.x()-coord_(0)) < MIN_THICKNESS) && (std::abs(p.y() - coord_(1)) < MIN_THICKNESS) && (std::abs(p.z()-coord_(2)) < MIN_THICKNESS) )return true;
		else return false;
	};

	/** 
     *      Not Equal operator
	 *
     *      @return a bool
     */
	bool operator!=(const Point3D& p) const 
	{
		return !operator==(p);
	};

	/** 
     *	operateur <
	 */

	bool operator<(Point3D p)
	{
		// Test de l'egalite avec tolerance ("precis")
		if (operator==(p) == true ) return false;

		// on est sur que les points ne sont pas "proches"
		if (x() < p.x() && y() < p.y() && z() < p.z()) return true;
		return false;
	}

	/** 
     *	operateur <=
	 */
	bool operator <=(Point3D p)
	{
		// on est sur que les points ne sont pas "proches"
		if (x() <= p.x() && y() <= p.y() && z() <= p.z()) return true;
		return false;
	}

	


};

std::ostream& operator<<(std::ostream& os, Point3D& pt);
std::istream& operator>>(std::istream& is, Point3D& pt);


//-----------------------------------------------------
// Distance au carre entre deux points 
//-----------------------------------------------------
double distance_carre(Point3D& a, Point3D& b);

//-----------------------------------------------------
// Distance entre deux points 
//-----------------------------------------------------
double distance(Point3D& a, Point3D& b);

//-----------------------------------------------------
// max et min entre deux point 
//-----------------------------------------------------

Point3D pmax(Point3D& a,Point3D& b);

Point3D pmin(Point3D& a,Point3D& b);

// inline methods

double tri_surface(Point3D& p1, Point3D& p2, Point3D& p3);

double tri_surface_grad(Point3D& p1, Point3D& p2, Point3D& p3, Point3D& elem_center, Eigen::Matrix<double, 3, 1>& gradient);

/** 
*		Get the cross product of 2 vectors    
*
*      @return a vector of type double
*/
Eigen::Matrix<double, 3, 1> prod_vec(Eigen::Matrix<double, 3, 1>& v, Eigen::Matrix<double, 3, 1>& w);

double prod_mixte(Eigen::Matrix<double, 3, 1> vec1, Eigen::Matrix<double, 3, 1> vec2,Eigen::Matrix<double, 3, 1> vec3);

bool iscolinear(Eigen::Matrix<double, 3, 1>& v, Eigen::Matrix<double, 3, 1>& w);

Point3D RandomRay(double r);
//----------------------------------------
void normalize(Eigen::Matrix<double, 3, 1>& v);

}


#endif 
