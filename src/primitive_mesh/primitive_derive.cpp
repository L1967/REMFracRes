
#include <primitive_mesh/plan.h>
#include <primitive_mesh/primitive_derive.h>

namespace fractures_intersect {

extern double pt_tol,TOL;
//-------------------------------------------------------
// Version 1.0 Pascal SIEGEL - 2000
//-------------------------------------------------------

//-------------------------------------------------------
// Classe Cpolyligne : public Cprimitive
// tableau de plusieurs Point3D 
//-------------------------------------------------------

//-------------------------------------------
// constructeurs/destructeur de la classe Cpolyligne
//-------------------------------------------

Cpolyligne::Cpolyligne(std::string sz,Tprimitive nv_type):
Cprimitive(sz,nv_type) {}

Cpolyligne::Cpolyligne(Point3D a,Point3D b):
Cprimitive("segment",Polyligne)
{
	tab_point.push_back(a);
	tab_point.push_back(b);
	tab_connec.push_back(true);
}

Cpolyligne::Cpolyligne(const Cpolyligne& poly):
Cprimitive((const Cprimitive&)poly) 
{
	Del_connec();
	
	int n = poly.tab_connec.size();
	if (n > 0)
	{
		tab_connec.resize(n);
		for (int i = 0 ; i < n ; i++)
			tab_connec[i] = poly.tab_connec[i];
	}
}

Cpolyligne::~Cpolyligne() 
{ 
	Del_connec();
}

//-------------------------------------------
// operateur = entre deux polylignes
//-------------------------------------------	
Cpolyligne& Cpolyligne::operator=(const Cpolyligne& poly)
{
	*((Cprimitive*)this) = (const Cprimitive&)poly;
	Del_connec();
	
	int n = poly.tab_connec.size();
	if (n > 0)
	{
		tab_connec.resize(n);
		for (int i = 0 ; i < n ; i++)
			tab_connec[i] = poly.tab_connec[i];
	}
	return *this;
}

//-------------------------------------------
// Insertion d'un point de maniere ordonnee
// sur un segment de polyligne (utile si d'autres
// points existent dans le segment)
// p0 = point origine du segment
// p1 = point extremite
// pt = point a inserer
// connec : booleen indique si on met a jour
// le tableau des connexions
//-------------------------------------------	
bool Cpolyligne::Insert(Point3D p0,Point3D p1,Point3D p,bool connec)
{
	if (p == p0 || p == p1) return false;

	Vector3D	vec= p1.vector() - p0.vector();
	double denom = vec.norm();
	if (denom == 0.0) return false;
	
	int i = Find_point(p0,0,pt_tol);
	if (i == -1) return false;
	
	int j = Find_point(p1,i,pt_tol);
	if (j == -1) return false;

	Vector3D vp= p.vector() - p0.vector();
	Point3D	point;
	double		ui = 0.,u = vp.norm()/denom;

	if (u > 1.) return false;

	for (int ij = i+1 ; ij <= j ; ij++)
	{
		point = tab_point[ij];
		if (p == point) return false;
		vp = point.vector() - p0.vector();
		ui = vp.norm()/denom;
		if (u < ui)
		{
			if (connec == true &&  tab_connec.size() > 0) tab_connec.insert(tab_connec.begin()+ij-1,true);
			tab_point.insert(tab_point.begin()+ij,p);
			return true;
		}
	}
	return false;
}

//-------------------------------------------
// Elimination d'une connexion dans le tableau
// des connexions
//-------------------------------------------	
void Cpolyligne::Del_connec(int i)
{
	int n = tab_connec.size();
	if (!n || i > n) return;
	
	if (i >= 0)
		tab_connec.erase(tab_connec.begin()+i);

	else if (i == -1)
	{
		// On detruit tout (valeur par defaut)
		tab_connec.clear();
	}
}

//-------------------------------------------
// Fonction testant si une polyligne
// est incluse dans un polygone
//-------------------------------------------	
bool Cpolyligne::Primitive_incluse(Cprimitive *prim,bool strict,bool cal_inter)
{
	if (prim->Get_type() != Polygone) return false;
//	if (!prim->In_Bounding_box(this)) return false;
	
	Cpolygone	*poly = (Cpolygone*)prim;
	int		n = tab_point.size();
	
	if (type == Polygone) n--;
	
	for (int i = 0  ; i < n ; i++)
	{
		if (poly->Point_ds_polygone(tab_point[i],pt_tol,strict,cal_inter) == false)
		return false;
	}
	
	return true;
}


//-------------------------------------------
// Fonction testant si le segment j de la polyligne p
// est inclus dans le segment i de la polyligne this
//-------------------------------------------	
bool Cpolyligne::Segment_inclus(Cpolyligne *p,int i,int j)
{
	if (i >= tab_point.size()-1) return false;
	if (j >= p->tab_point.size()-1) return false;

	Cpolyligne poly(tab_point[i],tab_point[i+1]);

	
	bool inclus = (poly.Point_ds_polyligne(p->tab_point[j],pt_tol) == true &&
				poly.Point_ds_polyligne(p->tab_point[j+1],pt_tol) == true);
	return inclus;
}


//-------------------------------------------
// Fonction testant tous les segments de la polyligne
// pour voir si un segment est inclus dans un autre
//-------------------------------------------	
void Cpolyligne::Test_segment_inclus()
{
	int	n = tab_point.size();
	if (n < 3) return;

	int	*flag = new int[n-1],j = 0;
	
	for (int i = 0 ; i < n-1 ; i++) flag[i] = -1;
	
	for (int i = 0 ; i < n-1 ; i++)
	{
		if (flag[i] != -1 || !tab_connec[i]) continue;

		bool	inclus = false;
		
		j = i+1;
		while (j < n-1 && !inclus)
		{
			if (flag[j] != -1 || !tab_connec[j]) break;
			
			if (Segment_inclus(this,i,j) == true)
			{
				flag[j] = i;
			}
			else if ((inclus = Segment_inclus(this,j,i)) == true)
			{
				flag[i] = j;
			}
			j++;
		}
	}
	
	for (int i = 0 ; i < n-1 ; i++)
	{
		j = flag[i];
		if (j != -1)
		{
			// le segment i est inclus dans le segment j
			Point3D	pj0 = tab_point[j],
					pj1 = tab_point[j+1];

			Insert(pj0,pj1,tab_point[i]);
			Insert(pj0,pj1,tab_point[i+1]);
			tab_connec[i] = false;
		}
	}

}

//-------------------------------------------
// Fonction testant si un point appartient
// a une polyligne
//-------------------------------------------	
bool Cpolyligne::Point_ds_polyligne(Point3D pt,double tol,bool cal_inter)
{
	if (!In_Bounding_box(pt)) return false;
	
	Point3D	pt1 = tab_point[0],pt2,pti;
	Vector3D	vl,v;
	double		den = 0.,u = 0.;

	for (size_t i = 1 ; i < tab_point.size() ; i++)
	{
		pt2 = tab_point[i];
		vl = pt2.vector() - pt1.vector();
		v = pt.vector() - pt1.vector();
		den = vl.norm()*vl.norm();
		u = v.dot(vl);
	
		// Calcul de la distance entre le point et le segment
		if (den > tol)
		{
			u /= den;
			
			double u_tol=pt_tol/vl.norm();

			if (u<0 && u>(-u_tol)) u=0.0;
		    if (u>1.0 && (u-1.0)<u_tol) u=1.0;

			if (u >= 0. && u <= 1.)
			{				
				pti.vector() = pt1.vector() + vl*u;
				if (distance(pti,pt) < tol)
				{
					// le point pt appartient au segment
					if (cal_inter == true)
					{
						// on est dans le cas de calcul d'intersection
						// on enrichit la polyligne avec le point
//						return Insert(pt1,pt2,pt);
						Insert(pt1,pt2,Point3D(pt,true));
					}
					return true;
				}
			}
		}
		// pt1 et pt2 sont coincidents a tol pres
		else
		{
			if (distance(pt1,pt) < tol || distance(pt2,pt) < tol)
			{
				// le point pt est proche de pt1 et/ou pt2 (a tol pres)
				if (cal_inter == true)
				{
					// on est dans le cas de calcul d'intersection
					// on enrichit la polyligne avec le point
					return Insert(pt1,pt2,Point3D(pt,true));
				}
				return true;
			}
		}
		pt1 = pt2;
	}
	return false;
}

//-------------------------------------------
// Fonction qui cree une nouvelle polyligne
// des la premiere discontinuite rencontree
//-------------------------------------------	
bool Cpolyligne::scinde_ligne(Cpolyligne *p)
{
	int n = tab_point.size();
	if (!n || !tab_connec.size()) return false;

	int j = -1; 

	for (int i = 0 ; i < n-1 ; i++)
	{
		if (tab_connec[i] == false) 
		{
			j = i;
			break;
		}
	}
	
	if (j == -1) return false;
	
	// la ligne est discontinue et ne possede qu'un segment
	// on elimine la ligne
	if (n == 2)
	{
		Del_point();
		Del_connec();
		return false;
	}

	// la discontinuite concerne le 1er segment
	// on elimine donc le 1er point de tab_point
	if (!j)
	{
		Del_point(0);
		Del_connec(0);
		return true;
	}
	
	// la discontinuite concerne le dernier segmentnamespace fractures_intersect {
	// on elimine donc le dernier point de tab_point
	if (j == n-2)
	{
		Del_point(n-1);
		Del_connec(j);
		return true;
	}
	
	// la nouvelle polyligne contient les points situes
	// apres la discontinuite (on les eliminera donc
	// de tab_point)
	for (int i = j+1 ; i < n ; i++)
	{
		p->tab_point.push_back(tab_point[i]);
		if (i < n-1) p->tab_connec.push_back(tab_connec[i]);
	}
	
	// on ne garde dans tab_point que les points
	// situes avant la discontinuite 
	for (int i = n-1 ; i > j ; i--)
	{
		Del_point(i);
		Del_connec(i-1);
	}
	return true;
}

//-------------------------------------------
// Fonction qui connecte une ligne a une autre
// si les deux lignes ont un point commun
//-------------------------------------------	
bool Cpolyligne::connecte_ligne(Cpolyligne *p)
{
	size_t n = tab_point.size(),
		np = p->tab_point.size();
	
	if (n < 2 || np < 2) return false;
	if (!calcule || !p->calcule) return false;
	
	int connec = 0;

	int j = -1; 

	for (int i = 0  ; i < np ; i += np-1)
	{
		if (tab_point[0] == p->tab_point[i])
		{
			j = i;
			break;
		}
	}
	
	if (j == -1)
	{
		for (size_t i = 0  ; i < np ; i += np-1)
		{
			if (tab_point[1] == p->tab_point[i])
			{
				j = i;
				connec++;
				break;
			}
		}
	}
	
	// les lignes ne sont pas connectees
	if (j == -1) return false;
		
	// Si les tableau des connexions n'existe pas,namespace fractures_intersect {
	// on le cree pour pouvoir l'utiliser dans les boucles
	// generales
	bool	connec1 = (tab_connec.size() > 0),
			connec2 = (p->tab_connec.size() > 0);

	if (!connec1)
	{
		tab_connec.resize(n-1);
		for (int i = 0 ; i < n-1 ; i++)
		{
			tab_connec[i] = true;
		}
	}

	if (!connec2)
	{
		p->tab_connec.resize(np-1);
		for (int i = 0 ; i < np-1 ; i++)
		{
			p->tab_connec[i] = true;
		}
	}

	if (connec == 1)
	{
		if (!j)
		{
			for (int i = 1 ; i < np ; i++)
			{
				tab_point.push_back(p->tab_point[i]);
				tab_connec.push_back(p->tab_connec[i-1]);
			}
		}
		else
		{
			for (int i = np-2 ; i >= 0 ; i--)
			{
				tab_point.push_back(p->tab_point[i]);
				tab_connec.push_back(p->tab_connec[i]);
			}
		}
	}
	else
	{
		Cpolyligne pl;	
		if (j == np-1)
		{
			for (int i = 0 ; i < np ; i++)
			{
				pl.tab_point.push_back(p->tab_point[i]);
				if (i < np-1) pl.tab_connec.push_back(p->tab_connec[i]);
			}
			for (int i = 1 ; i < n ; i++)
			{
				pl.tab_point.push_back(tab_point[i]);
				pl.tab_connec.push_back(tab_connec[i-1]);
			}
		}
		else
		{
			for (int i = np-1 ; i > 0 ; i--)
			{
				pl.tab_point.push_back(p->tab_point[i]);
				pl.tab_connec.push_back(p->tab_connec[i-1]);
			}
			for (int i = 0 ; i < n ; i++)
			{
				pl.tab_point.push_back(tab_point[i]);
				if (i < n-1) pl.tab_connec.push_back(tab_connec[i]);
			}
		}
		Del_point();
		n += np-1;
		for (int i = 0 ; i < n ; i++)
		{
			tab_point.push_back(pl.tab_point[i]);
			if (i < n-1) tab_connec.push_back(tab_connec[i]);
		}
	}
	return true;
}

//-------------------------------------------
// Intersection
//-------------------------------------------	
bool Cpolyligne::intersect(Cprimitive *prim,bool)
{
	if (!Intersect_Bounding_box(prim)) return false;
	
	bool	inter = false;
	switch (prim->Get_type())
	{
		case Polyligne:
//				if (!prim->In_Bounding_box(this)) break;
				inter = intersect_lignes(this,(Cpolyligne*)prim,true);
				if (inter == true)
				{
					// on fusionne les points coincidents eventuels
					Fusionne_points(pt_tol);
					prim->Fusionne_points(pt_tol);
				}
			break;
		
		case Point:
				inter = Point_ds_polyligne(prim->tab_point[0],
											pt_tol,true);
				if (inter == true)
				{
					// on fusionne les points coincidents eventuels
					Fusionne_points(pt_tol);
				}
			break;

		default:break;
	}
	return inter;
}

//-------------------------------------------------------
// Classe Cpolygone : public Cprimitive 
// tableau de Point3D + taille des elements
// de fonction de localisation ds un plan
//-------------------------------------------------------

//-------------------------------------------
// constructeur/destructeur de la classe Cpolygone
//-------------------------------------------
Cpolygone::Cpolygone(std::string sz) : Cpolyligne(sz,Polygone),
type_poly(Zone_affinement) {}

Cpolygone::Cpolygone(const Cpolygone& poly):
Cpolyligne((const Cpolyligne&)poly)
{
	type_poly = poly.type_poly;
}

Cpolygone::Cpolygone(const Cpolyligne& poly):
Cpolyligne(poly)
{
	type = Polygone;
	type_poly = Zone_affinement;
	nom = "polygone";
}

Cpolygone::~Cpolygone() {}

//-----------------------------------------------------------
//  Attributs
//----------------------------------------------------------- 
Tpolygone  Cpolygone::Get_type_poly() const { return type_poly; }
void Cpolygone::Set_type_poly(Tpolygone nv_type) { type_poly=nv_type; }

//-----------------------------------------------------------
//  Operations
//----------------------------------------------------------- 
Vector3D Cpolygone::Get_vec_normal()
{
	int n = tab_point.size();
	if (n < 3) return  Point3D().coord_;
	
	Point3D	p1 = tab_point[1];
	Vector3D	v,v3;
	Vector3D v1 =  p1.vector()-tab_point[0].vector();
	Vector3D vn;
	int			i = 2;
	
	normalize(v1);
	while (i < n )
	{
		Vector3D v2= tab_point[i].vector() - p1.vector();
		normalize(v2);
		v = v1.cross(v2);
		if (v.norm() > pt_tol) vn = v;
		i++;
	}
	return vn;
}


//-----------------------------------------------------------
// On recupere le plan contenant le polygone
//----------------------------------------------------------- 

Cplan* Cpolygone::get_plane()
{

	Vector3D vec = Get_vec_normal();
	if (vec.norm() == 0 || !tab_point.size()) return NULL;
	Cplan* plane = new Cplan("plan",vec,tab_point[0]) ;
	// primitive 0 du plan cree = polygone contour
	plane->tab_primitives.push_back(std::shared_ptr<Cpolygone>(new Cpolygone(*this)));
	std::dynamic_pointer_cast<Cpolygone>(plane->tab_primitives[0])->update_point2D(plane);

	return plane;
}

//-------------------------------------------
// Fonction testant si un point 3D
// est inclus dans un polygone (strictement si strict = true)
//-------------------------------------------	
bool Cpolygone::Point_ds_polygone(Point3D pt,double tol,bool strict,bool cal_inter)
{
	if (!In_Bounding_box(pt)) return false;

	std::shared_ptr<Cplan>	pl =  std::shared_ptr<Cplan>(get_plane());
	bool	in = pl->Point_ds_polygone(pl->tab_primitives[0].get(),pt);
	if (!in && !strict)
		in = Point_ds_polyligne(pt,tol,cal_inter);
	return in;
}

//-------------------------------------------
// Fonction testant si un point 3D
// est inclus dans un polygone (strictement si strict = true)
//-------------------------------------------
bool Cpolygone::Point_ds_polygone(Cplan* plane, Point3D pt,double tol,bool strict,bool cal_inter)
{
	if (!In_Bounding_box(pt)) return false;

	bool	in = plane->Point_ds_polygone(this,pt);
	if (!in && !strict)
		in = Point_ds_polyligne(pt,tol,cal_inter);
	return in;
}


void Cpolygone::update_point2D(Cplan* plane) {

	tab_point_2d_.clear();
	for (size_t i = 0 ; i <tab_point.size() ; i++)
	{
		tab_point_2d_.push_back(plane->Change_coord_2D(tab_point[i]));
	}
}


//-------------------------------------------
// Intersection
//-------------------------------------------	
bool Cpolygone::intersect(Cprimitive *prim,bool benrichi)
{
	if (!Intersect_Bounding_box(prim)) return false;
	
	bool	inter = false;
	switch (prim->Get_type())
	{
		case Polygone:
				inter = intersect_polygones(this,(Cpolygone*)prim);
			break;

		case Polyligne:
				inter = intersect_polygone_ligne(this,(Cpolyligne*)prim,benrichi);
			break;

		case Point:
				inter = Point_ds_polyligne(prim->tab_point[0],
											pt_tol,true);
				if (inter == true)
				{
					// on fusionne les points coincidents eventuels
					Fusionne_points(pt_tol);
				}
			break;

		default:break;
	}
	return inter;
}


//----------------------------------------------------
//
// Fonctions de gestion des intersections elementaires
//
//----------------------------------------------------

//----------------------------------------------------
// Intersection polyligne/polyligne 
//----------------------------------------------------
bool intersect_lignes(Cpolyligne *lin1,Cpolyligne *lin2,
					  bool blin)
{
	int	nlin1 = lin1->tab_point.size(),
		nlin2 = lin2->tab_point.size();

	if (nlin1 < 2 || nlin2 < 2) return false;
	
	Point3D		p0,p1,pl0 = lin2->tab_point[0],pl1,pt;
	Vector3D		vl1,vl2,vec;
	double			d = 0.,u1 = 0.,u2 = 0.;
	int				inter = 0;
	bool			connec1 = (lin1->tab_connec.size() > 0),
					connec2 = (lin2->tab_connec.size() > 0);
	
	// Si les tableau des connexions n'existe pas,
	// on le cree pour pouvoir l'utiliser dans les boucles
	// generales
	if (!connec1)
	{
		lin1->tab_connec.resize(nlin1-1);
		for (int i = 0 ; i < nlin1-1 ; i++)
		{
			lin1->tab_connec[i] = true;
		}
	}

	if (!connec2)
	{
		lin2->tab_connec.resize(nlin2-1);
		for (int i = 0 ; i < nlin2-1 ; i++)
		{
			lin2->tab_connec[i] = true;
		}
	}
	
	for (int i = 1, ij = 1 ; i < nlin2 ; i++)
	{
			pl1 = lin2->tab_point[ij];

		// test d'intersection uniquement si les points
		// du segment sont connectes
		if (lin2->tab_connec[ij-1] == true)
		{
			vl2 = pl1.vector() - pl0.vector();
			p0 = lin1->tab_point[0];
			nlin1 = lin1->tab_point.size();
			for (int j = 1 , jj = 1 ; j < nlin1 ; j++)
			{
					p1 = lin1->tab_point[jj];
				if (lin1->tab_connec[jj-1] == true)
				{
					vl1 = p1.vector() - p0.vector();
 					vec = pl0.vector() - p0.vector();
					d = prod_mixte(vec,vl2,vl1);
					// si d != 0, pas d'intersection possible
					// (segments non coplanaires)
					if (d == 0.0)
					{
						Point3D pp;
						pp.calculated_ = true;
						if (!iscolinear(vl2,vl1)) // lignes non parallelles
						{
							// On travaille dans le plan des segments
							vec = vl2.cross(vl1);
							Cplan *planlin = new Cplan("plan",vec,p0);
							Vector3D vl1p(vl1),vl2p(vl2);
							vec = pl0.vector() - p0.vector();
							vec = planlin->Projection_sur(vec);
							vl1p = planlin->Projection_sur(vl1);
							vl2p = planlin->Projection_sur(vl2);
							d = vl2p[1]*vl1p[0] - vl2p[0]*vl1p[1];
							u2 = (vl1p[1]*vec[0] - vl1p[0]*vec[1])/d;
							u1 = (vl2p[1]*vec[0] - vl2p[0]*vec[1])/d;
							

							double u1_tol=pt_tol/vl1.norm();
							double u2_tol=pt_tol/vl2.norm();

							if (u1<0 && u1>(-u1_tol)) u1=0.0;
					        if (u1>1.0 && (u1-1.0)<u1_tol) u1=1.0;

							if (u2<0 && u2>(-u2_tol)) u2=0.0;
					        if (u2>1.0 && (u2-1.0)<u2_tol) u2=1.0;

							if (u1 >= 0. && u1 <= 1. && u2 >= 0. && u2 <= 1.)
							{
								// le point intersection des 2 segments existe
								pt.vector() = p0.vector() + vl1*u1;
								pp = pt;
								if (blin == true)
								{
									if (lin1->Insert(p0,p1,pp) == true) jj++;
								}
								if (lin2->Insert(pl0,pl1,pp) == true) ij++;
								inter++;
							}

							delete planlin;
						}
						
							// lignes paralleles ou confondues
							Cpolyligne l1(p0,p1);
							
							if (l1.Point_ds_polyligne(pl0,pt_tol) == true)
							{
								if (blin == true)
								{
									pp = pl0;
									if (lin1->Insert(p0,p1,pp) == true) jj++;
								}
								inter++;
							}
							if (l1.Point_ds_polyligne(pl1,pt_tol) == true)
							{
								if (blin == true)
								{
									pp = pl1;
									if (lin1->Insert(p0,p1,pp) == true) jj++;
								}
								inter++;
							}

							Cpolyligne l2(pl0,pl1);
							if (l2.Point_ds_polyligne(p0,pt_tol) == true)
							{
								pp = p0;
								if (lin2->Insert(pl0,pl1,pp) == true) ij++;
								inter++;
							}
							if (l2.Point_ds_polyligne(p1,pt_tol) == true)
							{
								pp = p1;
								if (lin2->Insert(pl0,pl1,pp) == true) ij++;
								inter++;
							}
					
					}
				}
				p0 = p1;
				jj++;
			}
		}
		pl0 = pl1;
		ij++;
	}
	
	// Si les tableau des connexions n'existait pas,
	// on le nettoie
	if (!connec1)
		lin1->Del_connec();

	if (!connec2)
		lin2->Del_connec();

	return (inter != 0);
}

//-------------------------------------------
// Intersection polygone/ligne : revient a une
// intersection polyligne/polyligne 
//-------------------------------------------	
bool intersect_polygone_ligne(Cpolygone *poly,Cpolyligne *lin,bool bpoly)
{
	// On nettoie le tableau des connexions du polygone
	// (pas de discontinuites possibles)
	poly->Del_connec();
	
	// On calcule les points intersections entre les
	// segments du polygone et de la polyligne
	if (intersect_lignes((Cpolyligne*)poly,lin,bpoly) == true)
	{
		// On verifie qu'il n'y a pas de points coincidents
		// dans le polygone et la polyligne enrichis
		// par les intersections
		// la tolerance est pt_tol
		poly->Fusionne_points(pt_tol);
		lin->Fusionne_points(pt_tol);
		return true;
	}
	return false;
}

//-------------------------------------------
// Intersection polygone/polygone
// retour = true si l'intersection existe
//-------------------------------------------	
bool intersect_polygones(Cpolygone *poly1,Cpolygone *poly2)
{
	// On nettoie les tableaux des connexions des polygones
	// (pas de discontinuites possibles)
	poly1->Del_connec();
	poly2->Del_connec();
	
	// On calcule les points intersections entre les
	// segments des polygones (passage aux polylignes)
	if (intersect_lignes((Cpolyligne*)poly1,(Cpolyligne*)poly2,true) == true)
	{
		// On verifie qu'il n'y a pas de points coincidents
		// dans les polygones enrichis par les intersections
		// la tolerance est pt_tol
		poly1->Fusionne_points(pt_tol);
		poly2->Fusionne_points(pt_tol);
		return true;
	}
	return false;
}

}
