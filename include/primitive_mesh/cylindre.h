//-------------------------------------------
// Version 1.0 - Pascal SIEGEL - 2000
//-------------------------------------------

#if !defined(__CYLINDRE_H__)
	#define __CYLINDRE_H__

#include <primitive_mesh/polyedre.h>
#include <primitive_mesh/faille.h>


namespace fractures_intersect {

//--------------------------------------------------------------------
// class Ccylindre : public Cpolyedre
// 
// cylindres circulaires droits uniquement
// definition : directrice du cylindre delimite par un segment 
// (2 points 3D : c0 = origine, c1 = extremite)
// et un rayon (constant)
// pour la discretisation ulterieure ou l'intersection non analytique,
// cette classe herite du polyedre (facetisation sous forme de plans)
//--------------------------------------------------------------------

class Ccylindre : public Cpolyedre
{
	// Variables
public:
	std::shared_ptr<Cfaille> faille_support;// points du segment de directrice

	Point3D	    ptExtrem[2];// points du segment de directrice

	double		a_rayon;		// rayon du cylindre a de l'elipse
	double		b_rayon;		// rayon du cylindre b de l'elipse

	int			nb_facettes;
				
public:
	// Construction/destruction
    Ccylindre(std::shared_ptr<Cfaille> faille_support, int n = 12);
    virtual ~Ccylindre();
	
	// Attributs
	

	// Operations
	bool		Actualise_facettes(int n);
};

}
#endif //__CYLINDRE_H__
