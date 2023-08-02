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


// system includes

// project includes 

// project includes 
#include<primitive_mesh/point.h>
#include <time.h>
#include <iomanip>
#include <fstream>

// local includes 

namespace fractures_intersect {

/**
 *    coordinates point  
 *   
 *    @author Pascal Siegel
 *    @date 05/2013
 */


//-----------------------------------------------------
// operateur d'ecriture pour la structure Tpoint3D
//-----------------------------------------------------
std::ostream& operator<<(std::ostream& os, Point3D& pt)
{
	os<<std::setiosflags(std::ios::fixed)<<std::setprecision(6)<<std::setiosflags(std::ios::showpoint);
    os<<std::setw(20)<<pt.x()<<" "<<std::setw(20)<<pt.y()<<" " <<std::setw(20)<<pt.z();
	return os;
} 

//-----------------------------------------------------
// operateur de lecture pour la structure Tpoint3D
//-----------------------------------------------------
std::istream& operator>>(std::istream& is, Point3D& pt)
{
	int cal;
	is>>pt.x()>>pt.y()>>pt.z()>>cal;
	pt.calculated_ = cal; 
	return is;
} 



//-----------------------------------------------------
// Distance au carre entre deux points 
//-----------------------------------------------------
double distance_carre(Point3D& a, Point3D& b)
{
	double dist = (a.vector()-b.vector()).norm();
	return dist*dist;
	
}

//-----------------------------------------------------
// Distance entre deux points 
//-----------------------------------------------------
double distance(Point3D& a, Point3D& b)
{
	double dist = (a.vector()-b.vector()).norm();
	return dist;
}


//-----------------------------------------------------
// max et min entre deux point 
//-----------------------------------------------------

Point3D pmax(Point3D& a,Point3D& b)
{
	Point3D pt;

	pt.x()=std::max(a.x(),b.x());
	pt.y()=std::max(a.y(),b.y());
	pt.z()=std::max(a.z(),b.z());

	return pt;
}


Point3D pmin(Point3D& a,Point3D& b)
{
	Point3D pt;

	pt.x()=std::min(a.x(),b.x());
	pt.y()=std::min(a.y(),b.y());
	pt.z()=std::min(a.z(),b.z());

	return pt;
}



// inline methods

double tri_surface(Point3D& p1, Point3D& p2, Point3D& p3)
{
	Eigen::Matrix<double, 3, 1> v = p2.vector()-p1.vector();
	Eigen::Matrix<double, 3, 1> w = p3.vector()-p1.vector();
	Eigen::Matrix<double, 3, 1> p = v.cross(w);

	return 0.5* p.norm();
};

double tri_surface_grad(Point3D& p1, Point3D& p2, Point3D& p3, Point3D& elem_center, Eigen::Matrix<double, 3, 1>& gradient)
{
	Eigen::Matrix<double, 3, 1> v = p2.vector()-p1.vector();
	Eigen::Matrix<double, 3, 1> w = p3.vector()-p1.vector();
	Eigen::Matrix<double, 3, 1> p = v.cross(w);

	Eigen::Matrix<double, 3, 1>  face_center = (p1.vector()+p2.vector()+p3.vector())/3.0;

	Eigen::Matrix<double, 3, 1> externe_direction = face_center - elem_center.vector();

	if (externe_direction.dot(p) < 0 ) p = -p;

	return 0.5*gradient.dot(p);
};


/** 
*		Get the cross product of 2 vectors    
*
*      @return a vector of type double
*/
Eigen::Matrix<double, 3, 1> prod_vec(Eigen::Matrix<double, 3, 1>& v, Eigen::Matrix<double, 3, 1>& w)
{
	return v.cross(w);
	/*Eigen::Matrix<double, 3, 1> p;
	p(0)=  v(1)*w(2)-v(2)*w(1);
	p(1)=  v(2)*w(0)-v(0)*w(2);
	p(2)=  v(0)*w(1)-v(1)*w(0);
	return p;*/
};

//----------------------------------------
// Produit mixte (vec1,vec2,vec3) 
// (= determinant 3x3 de la matrice
// formee par les trois vecteurs mis en colonne)
//----------------------------------------
double prod_mixte(Eigen::Matrix<double, 3, 1> vec1, Eigen::Matrix<double, 3, 1> vec2,Eigen::Matrix<double, 3, 1> vec3)
{
	normalize(vec1);
    normalize(vec2);
	normalize(vec3);

	double det = vec1.cross(vec2).dot(vec3);
	if (std::abs(det) < TOL) det = 0.;
	return det;
}



//----------------------------------------
// Verification si deux vecteurs sont colineaires
// true = colineaires
//----------------------------------------
bool iscolinear(Eigen::Matrix<double, 3, 1>& v, Eigen::Matrix<double, 3, 1>& w)
{
	return (v.cross(w).norm()<TOL);
}


//----------------------------------------
// Vecteur normalise
//----------------------------------------
void normalize(Eigen::Matrix<double, 3, 1>& v)
{
	double x = v.norm();
	if (x != 0.)
	{
		v = v/x;
		if (std::abs(v[0]) < TOL) v[0] = 0.;
		if (std::abs(v[1]) < TOL) v[1] = 0.;
		if (std::abs(v[2]) < TOL) v[2] = 0.;
	}
}

//-----------------------------------------------------
// Ray tracing : on genere un point de maniere aleatoire
// dans l'espace en le placant a une distance r donnee
// (adapte de la fonction de Joseph O'Rourke, avril 1998)
//-----------------------------------------------------
Point3D RandomRay(double r)
{
	srand((unsigned)time(NULL));

	// Generate a random point on a sphere of radius 1.
	// the sphere is sliced at z, and a random point at angle t
	// generated on the circle of intersection.
	double	z = 2.0 * (double)rand()/RAND_MAX - 1.0,
			t = 2.0 * pi * (double)rand() / RAND_MAX,
			w = std::sqrt(1. - z*z),
			x = r * w * std::cos(t),
			y = r * w * std::sin(t);

	z *= r;
	
	return Point3D(x,y,z);
}
}



