
#include <primitive_mesh/boite.h>
#include <boost/random.hpp>
#include <boost/format.hpp>
#include <map>

namespace fractures_intersect {
 
extern double dround(double,double);
extern double pt_tol,TOL;
const long MAXCAR = 120;
//-------------------------------------------------------
//
// fonctions de la classe classe Cpolyedre :
// primitive volumique 
// limite par un ensemble d'objet plan
// un forme un volume ferme
//
// Version 1.0 Pascal SIEGEL - 2000
//-------------------------------------------------------

//-------------------------------------------
// constructeur/destructeur de la classe Cpolyedre
//-------------------------------------------
Cpolyedre::Cpolyedre(std::string sz,Tprimitive nv_type):
Cbase_primitive(sz,nv_type),type_poly(Zone_affinement)
{}

Cpolyedre::~Cpolyedre() 
{
/*	for(int i=0;i<tab_plan.size();i++)
	{
		delete tab_plan[i];
	}
	tab_plan.RemoveAll();*/
	
	Del_plans();
}

//---------------------------------------------
// Attributs
//---------------------------------------------
Tpolygone Cpolyedre::Get_type_poly() const { return type_poly; }
void Cpolyedre::Set_type_poly(Tpolygone nv_type) { type_poly = nv_type; }
void Cpolyedre::Set_type_poly(bool typ) 
{
	if (typ==true) type_poly=Zone_affinement;
	else type_poly=Trou;
}

void Cpolyedre::Set_nom(std::string nv_nom)
{

	std::string		szn,sz;

	if (!nv_nom.empty()) nom = nv_nom;

	for(size_t i=0;i<tab_plan.size();i++)
	{
		sz  =  boost::str(boost::format("plan%1%") % (i) );
		szn = nom + "_" + sz;
		tab_plan[i]->Set_nom(szn);
	}

}

//-------------------------------------------
// Determination de la "bounding box" du polyedre
//-------------------------------------------	
void Cpolyedre::Set_Bounding_box()
{
	if (boundingBox_.bounding_set_== false)
	{
	
		int n = tab_plan.size();
		if (!n) return;
		
		for (int i = 0 ; i < n ; i++)
		{
			tab_plan[i]->Set_Bounding_box();
			boundingBox_.add_bounding_box(tab_plan[i]->get_boundingBox());
		}

		boundingBox_.bounding_set_ = true;
	}
}

//---------------------------------------------
// fonction qui retourne oui si le point pt
// est dans le polyedre
// si strict = true (par defaut) : la limite du
// polyedre est consideree comme externe (point
// strictement a l'interieur)
//---------------------------------------------
bool Cpolyedre::Point_ds_polyedre(Point3D pt,bool strict)
{
	int n = tab_plan.size();

	// Nombre de plans insuffisant : le polyedre n'est pas ferme
	if (n < 4) return false;

	// On actualise la boite englobante du polyedre
	Set_Bounding_box();
	
	// On verifie que le point est situe dedans
	if (!In_Bounding_box(pt)) return false;
	
	// Si oui, on verifie que le point n'est pas sur une face
	// du polyedre
	for (int i = 0 ; i < n ; i++)
	{
		if (!tab_plan[i]->tab_primitives.size()) return false;
		if (tab_plan[i]->Point_ds_plan(pt,pt_tol) == true)
		{
			return !strict;
		}
	}


	// Finalement, on fait du "Ray tracing"
	Point3D	pti,ptj;
	Vector3D	nj,v;
	Cpolyligne	seg;
	double		ray = distance(boundingBox_.pt_min_,boundingBox_.pt_max_),
				num = 0.,denom = 0.,d = 0.,u = 0.;
	int			intersect = 0,degenere = 0;
	
	int i = 0;

    Ray_degenere:
	while (i++ < n)
	{
		intersect = 0;
		pti.vector() = pt.vector() + RandomRay(ray).vector();
		seg = Cpolyligne(pt,pti);
		v = pti.vector() - pt.vector();
		degenere = 0;
		for (int j = 0 ; j < n ; j++)
		{
			std::shared_ptr<Cpolygone> poly = std::dynamic_pointer_cast<Cpolygone>(tab_plan[j]->tab_primitives[0]);
			if (poly->Intersect_Bounding_box(&seg) == true)
			{
				degenere++;
//				ptj = poly[j].tab_point[0]->pt;
				ptj = tab_plan[j]->Get_point_ref();
				nj = tab_plan[j]->Get_vec_normal();
				d = nj.dot(ptj.vector());
				num = d - nj.dot(pt.vector());
			    denom = v.dot(nj);

				// Segment parallele au plan
			    if (denom == 0.0)
				{  
					// pti est sur le plan => mauvais rayon, on recommence
					if (num == 0.0)
					{
						degenere = 0;
						goto Ray_degenere;
					}
				}
				else u = num / denom;
				u = dround(u,TOL);

				if (0.0 < u && u < 1.0) 
				{
					// Point d'intersection entre le plan et
					// le segment
					ptj.vector() = pt.vector() + v*u;

					// Si ptj est dans le polygone contour,
					// on incremente intersect
					if (poly->Point_ds_polygone(tab_plan[j].get(),ptj,pt_tol,strict) == true)
						intersect++;
				}
				else if (u == 1.0 || u == 0.0)
				{
					degenere = 0;
					goto Ray_degenere;
				}
			}
		}
		if (degenere > 0) break;
	}

	if(degenere > 0 && (intersect % 2) == 1) return true;
	return false;
}

//---------------------------------------------
// Operations
//---------------------------------------------
void Cpolyedre::Actualise_plans()
{
	int			n = tab_plan.size();

	if (n < 5) return;
	
	std::shared_ptr<Cprimitive> prim0 = tab_plan[0]->tab_primitives[0];
	std::shared_ptr<Cprimitive> prim1 = tab_plan[1]->tab_primitives[0];
	std::string		szn,sz;
	
	for(int i = 2, i1 = 0 , i2 = 0 ; i < n ; i++)
	{
		sz = boost::str(boost::format("plan%1%") % (i) );
		szn = nom + "_" + sz;
		tab_plan[i] = std::shared_ptr<Cplan>(new Cplan());
		tab_plan[i]->Set_nom(szn);
		tab_plan[i]->tab_primitives.resize(1);
		sz += "_contour";
		tab_plan[i]->tab_primitives[0] = std::shared_ptr<Cpolygone>(new Cpolygone(sz));
		std::shared_ptr<Cprimitive> prim = tab_plan[i]->tab_primitives[0];
		prim->Set_taille_element(taille_element);
		prim->Set_ind_materiaux(ind_materiaux);
		i1 = i-1;
		i2 = n-i;
		if (i == 2) i2 = 0;
		if (i == n-1) i1 = 0;
		prim->tab_point.push_back(prim0->tab_point[i1]);
		prim->tab_point.push_back(prim0->tab_point[i-2]);
		prim->tab_point.push_back(prim1->tab_point[i2]);
		prim->tab_point.push_back(prim1->tab_point[n-i-1]);
		prim->tab_point.push_back(prim0->tab_point[i1]);
	
		tab_plan[i]->Set_point_ref(prim->tab_point[0]);
		tab_plan[i]->Set_taille_element(taille_element);
		tab_plan[i]->Set_ind_materiaux(ind_materiaux);
		if (tab_plan[i]->Actualise_repere(true) == true)
			tab_plan[i]->Set_equation_plan();


	}

	for(size_t j=0;j<tab_plan.size();j++) {
		tab_plan[j]->Set_ind_materiaux(ind_materiaux);
		std::dynamic_pointer_cast<Cpolygone>(tab_plan[j]->tab_primitives[0])->update_point2D(tab_plan[j].get());
	}
}

//---------------------------------------------
// On regarde dans les tableaux de primitives
// des plans si des polylignes sont discontinues
// ou non
//---------------------------------------------namespace fractures_intersect {
void Cpolyedre::scinde_lignes()
{
	for (size_t i = 0 ; i < tab_plan.size() ; i++)
	{
			tab_plan[i]->scinde_lignes();
	}
}

//---------------------------------------------
// Fusion des primitives identiques
// issues des intersections
//---------------------------------------------
void Cpolyedre::Fusionne_primitives()
{
	for (size_t i = 0 ; i < tab_plan.size() ; i++)
	{
			tab_plan[i]->Fusionne_primitives();
	}
}

//---------------------------------------------
// Intersection polyedre/primitive de base
//---------------------------------------------
void Cpolyedre::intersect(Cbase_primitive *prim)
{
	for (size_t i = 0 ; i < tab_plan.size() ; i++)
	{
			tab_plan[i]->intersect(prim);
	}
}


//---------------------------------------------
// Intersection polyedre/primitive de base
//---------------------------------------------
bool Cpolyedre::check_possible_intersect_polyedre(Cbase_primitive *prim)
{

	Cpolyedre* poly = dynamic_cast<Cpolyedre*>(prim);
	if (poly == 0) return false;

	if (Intersect_Bounding_box(prim) ==  false) return false;

	for (size_t i = 0 ; i < tab_plan.size() ; i++)
	{
		if (tab_plan[i]->Intersect_Bounding_box(poly) ==  false) continue;

		for (size_t j = 0 ; j < poly->tab_plan.size() ; j++)
		{
			if( tab_plan[i]->check_intersect_plan(poly->tab_plan[j].get(),pt_tol) == true) return true;
		}

	}

	return false;
}

//---------------------------------------------
// Intersection abutting fracture polyedre
//---------------------------------------------
bool Cpolyedre::abutting_fracture(Cbase_primitive *prim)
{

	Cfaille* faille = dynamic_cast<Cfaille*>(prim);
	if (faille == 0) return false;

	if (Intersect_Bounding_box(prim) ==  false) return false;


	std::vector<std::pair<int, double> > vec;

	bool intersect = false;

	for (size_t i = 0 ; i < tab_plan.size() ; i++)
	{
		if (tab_plan[i]->Intersect_Bounding_box(faille) ==  false) continue;

		double dist = 0.0;

		if (faille->intersect_plan_limitated_dist(tab_plan[i].get(),pt_tol,dist)) {
			std::pair<int, double> plan_dist(i, dist);
			vec.push_back(plan_dist);
		}
	}

	if (vec.empty()) return false;

	std::sort(vec.begin(), vec.end(), [](const std::pair<int,double> &left, const std::pair<int,int> &right) {
	    return left.second > right.second;
	});

	for (size_t i = 0 ; i < vec.size() ; i++)
	{
		int index_plan = vec[i].first;

		if (faille->intersect_plan_limitated(tab_plan[index_plan].get(),pt_tol)) {
			//std::cout << "simplify to box poly from index plan " << index_plan << std::endl;
			faille->simplify_to_box_poly();
			intersect = true;
		}


	}


	return intersect;
}


//-------------------------------------------
// Elimination des primitives issues des
// intersections
//-------------------------------------------	
void Cpolyedre::Del_intersect()
{
	int n = tab_plan.size();
	if (!n) return;
	
	for (int i = n-1 ; i >= 0 ; i--)
	{
		if (!tab_plan[i]->tab_primitives.size())
			Del_plans(i);
		
		else if (tab_plan[i] != NULL)
			tab_plan[i]->Del_intersect();
	}
}

//-------------------------------------------
// Elimination d'un plan du tableau
// des plans
//-------------------------------------------
void Cpolyedre::Del_plans(int i)
{
	Del_plan(tab_plan,i);
}



//-----------------------------------------------------------
// Indice materiau du polyedre
//----------------------------------------------------------- 

void Cpolyedre::Set_ind_materiaux(int nv_ind)  
{ 
	ind_materiaux=nv_ind;
	for(size_t j=0;j<tab_plan.size();j++) tab_plan[j]->Set_ind_materiaux(nv_ind);
}




//-------------------------------------------
// Elimination d'un polyedre du tableau
// des polyedres
//-------------------------------------------	
void Del_polyedre(std::vector<std::shared_ptr<Cpolyedre> >& tab,int i)
{
	int n = tab.size();
	if (!n || i > n) return;
	
	if (i >= 0)
	{
		tab.erase(tab.begin()+i);
	}
	else if (i == -1)
	{
		tab.clear();
	}
}

//-----------------------------------------------------
// write operator for plan mesh triangle
//-----------------------------------------------------
std::ostream& operator<<(std::ostream& os, Cpolyedre* poly)
{

	int n = poly->tab_plan.size();

	for (int i = 0 ; i < n ; i++)
	{
		os<<poly->tab_plan[i];
	}


	return os;

}

}

