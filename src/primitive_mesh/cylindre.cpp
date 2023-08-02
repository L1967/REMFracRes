//-------------------------------------------
// Version 1.0 - Pascal SIEGEL - 2000
//-------------------------------------------


#include <cmath>
#include <primitive_mesh/cylindre.h>

namespace fractures_intersect {

const long MAXCAR = 120;
//////////////////////////////////////////////////////////////////
//-------------------------------------------
// constructeur/destructeur de la classe Ccylindre
//-------------------------------------------

Ccylindre::Ccylindre(std::shared_ptr<Cfaille> faille, int n):
Cpolyedre("cyl",Cylindre)
{	
	faille_support = faille;// points du segment de directrice


	ptExtrem[0].coord_ = (faille->tab_primitives[0]->tab_point[0].coord_ + faille->tab_primitives[0]->tab_point[1].coord_)/2.0;

	ptExtrem[1].coord_ = (faille->tab_primitives[0]->tab_point[2].coord_ + faille->tab_primitives[0]->tab_point[3].coord_)/2.0;


	double length = faille->ext2L_selon_lpgp;
	a_rayon = faille->ext2L_selon_lpgp/2.0;
	b_rayon = faille->stress_zone_max_;

	/*if(a_rayon > b_rayon ){

		if(a_rayon/b_rayon > 5) a_rayon = b_rayon;

	}else{

		if(a_rayon/b_rayon < 5) a_rayon = b_rayon ;

	}*/

	Actualise_facettes(n);
}

Ccylindre::~Ccylindre() {}





//-----------------------------------------------------
// Determination des n facettes (plans) du cylindre
// stockees dans tab_plan. Plan 0 = face origine,
// Plan 1 = face extremite
//-----------------------------------------------------
bool Ccylindre::Actualise_facettes(int n)
{
	nb_facettes = n;
	ptExtrem[0].coord_ = (faille_support->tab_primitives[0]->tab_point[0].coord_ + faille_support->tab_primitives[0]->tab_point[1].coord_)/2.0;

	ptExtrem[1].coord_ = (faille_support->tab_primitives[0]->tab_point[2].coord_ + faille_support->tab_primitives[0]->tab_point[3].coord_)/2.0;


	double length = faille_support->ext2L_selon_lpgp;
	a_rayon = faille_support->ext2L_selon_lpgp/2.0;
	b_rayon = faille_support->stress_zone_max_;



	if (a_rayon == 0.0 || b_rayon == 0.0) return false;

	n = n+2;

	// Nombre de facettes insuffisant pour decrire le cylindre
	if (n < 5) return false;
	
	Vector3D T = ptExtrem[1].coord_ - ptExtrem[0].coord_;
	T.normalize();

	if (T.norm() == 0.) return false;
	
	int nplan = tab_plan.size();
	
	if (nplan)
	{
		tab_plan.clear();
	}

	tab_plan.resize(n);
	Point3D	pt,ptres,pt0;
	Vector3D	vn = -T;
	double		pi2 = 2.*pi,dteta = pi2/(n-2),teta = 0.,eps = 1.;
	std::string sz = nom + "_plan0";

	tab_plan[0] = std::shared_ptr<Cplan>(new Cplan(sz,vn,ptExtrem[0]));
	tab_plan[0]->update_RS(faille_support->Get_R());
	tab_plan[0]->tab_primitives.resize(1);
	tab_plan[0]->tab_primitives[0] = std::shared_ptr<Cpolygone>(new Cpolygone("plan0_contour"));
	
	sz = nom + "_plan1";
	tab_plan[1] = std::shared_ptr<Cplan>(new Cplan(sz,T,ptExtrem[1]));
	tab_plan[1]->update_RS(faille_support->Get_R());
	tab_plan[1]->Set_taille_element(taille_element);
	tab_plan[1]->Set_ind_materiaux(ind_materiaux);
	tab_plan[1]->tab_primitives.resize(1);
	tab_plan[1]->tab_primitives[0] = std::shared_ptr<Cpolygone>(new Cpolygone("plan1_contour"));

	std::shared_ptr<Cprimitive>	prim0 = tab_plan[0]->tab_primitives[0],
				prim1 = tab_plan[1]->tab_primitives[0];
	
	prim0->Set_taille_element(taille_element);
	prim0->Set_ind_materiaux(ind_materiaux);
	prim1->Set_taille_element(taille_element);
	prim1->Set_ind_materiaux(ind_materiaux);
	


	for(int i = 0 ; i < n-2 ; i++)
	{
		eps *= -0.01;
		teta = i*dteta;
		pt.coord_[0] = (a_rayon + precis*eps)*cos(teta);
		pt.coord_[1] = (b_rayon + precis*eps)*sin(teta);
		pt.coord_[2] = 0.;
		ptres = tab_plan[0]->Change_coord_3D(pt);
		prim0->tab_point.push_back(ptres);


	}
	prim0->tab_point.push_back(prim0->tab_point[0]);
	
	Vector3D vec = ptExtrem[1].coord_ - ptExtrem[0].coord_;
	prim1->tab_point.resize(n-1);
	for(int i = 0 ; i < n-1 ; i++)
	{
		ptres.coord_ = prim0->tab_point[i].coord_ + vec;
		prim1->tab_point[n-i-2] = ptres;
	}
	

	Cpolyedre::Actualise_plans();	
	Cpolyedre::Set_Bounding_box();
	return true;
}

}

