// Version 1.0 - Pascal SIEGEL - 2000

#ifndef __VARGLOB_H__
#define __VARGLOB_H__

#include <cmath>

namespace fractures_intersect {

double			precis = 1.e-4, // faire dialog pour optimiser
				pt_tol = 1.e-4,
				TOL = 1.e-6,
				pi = std::acos(double(-1)),
				rad = pi/180.0,
				irad = 180.0/pi;

}

#endif //__VARGLOB_H__
