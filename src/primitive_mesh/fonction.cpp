#include <stdlib.h>
#include <float.h>
#include <primitive_mesh/fonction.h>

namespace fractures_intersect {
//--------------------------------------------
// retourne un nombre aleatoire entre 0 et n
//---------------------------------------------

int random(int n)
{
	int i;
	double div;

	div=(double) n/ (double) (RAND_MAX+1);
	i=(int) (rand()*div);

	return i;
}

//--------------------------------------------
// qui arrondi un double a la precision tol
//--------------------------------------------

double dround(double x, double tol)
{
      bool flag;
      static volatile double X;

	  x/=tol;

      if (true == (flag = (x < 0.0)))
            X = -x;
      else  X = x;


			
      X +=  1.0/ (DBL_EPSILON);
      X -=  1.0/ (DBL_EPSILON);

	  X*=tol;


      return ((flag) ? -X : X);
}

}
