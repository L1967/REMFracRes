#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include<primitive_mesh/point.h>

namespace fractures_intersect {

 void exactinit();

 double orient2d(double *pa, double *pb, double *pc);   
 
 double orient3d(double *pa, double *pb, double *pc, double *pd);  
 
 double incircle(double *pa, double *pb, double *pc, double *pd);  
 
 double insphere(double *pa, double *pb, double *pc, double *pd, double *pe);

 void tricircumcenter(double *a, double *b, double *c, double *circumcenter);

 void tricircumcenter3d(double *a,double *b, double *c, double *circumcenter);

 void tetcircumcenter(double *a, double *b, double *c, double *d, double *circumcenter);

 double Orient3D(Point3D tpa,Point3D tpb,Point3D tpc,Point3D tpd);
 
 double Orient2D(Point3D tpa,Point3D tpb,Point3D tpc);

 double uniformdoublerand();

}

