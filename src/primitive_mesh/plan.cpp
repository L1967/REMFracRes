
#include <time.h>
#include <primitive_mesh/polyedre.h>
#include <primitive_mesh/faille.h>
#include <primitive_mesh/plan.h>
#include <boost/random.hpp>
#include <boost/format.hpp>

namespace fractures_intersect {

extern double pt_tol,TOL;

//-------------------------------------------
//
// fonctions de la classe Cplan
// 
// equation du plan + limites externes
// et internes : contour + tableau de 
// primitives du plan 
//
// Version 1.0 Pascal SIEGEL - 2000
//-------------------------------------------

//-------------------------------------------
// constructeur/destructeur de la classe Cplan
//-------------------------------------------

Cplan::Cplan() : Cbase_primitive("plan",Plan) {
}

Cplan::Cplan(std::string sz,Vector3D& vec,Point3D& pt):
Cbase_primitive(sz,Plan)
{	
	
	bequation = false;
	pt0_ = pt;
	if (vec.norm() > 0.0)
	{
		T = vec;
		normalize(T);
		Set_equation_plan();
		Actualise_repere();
	}
}

Cplan::Cplan(const Cplan& plan):
Cbase_primitive((const Cbase_primitive&)plan)
{
	
	bequation = false;
	pt0_ = plan.pt0_;
	T = plan.T;
	R = plan.R;
	S = plan.S;
	Set_equation_plan();
}

Cplan::~Cplan() 
{
}

//-------------------------------------------
// Elimination d'une primitive du tableau
// des primitives
//-------------------------------------------
void Cplan::Del_primitives(int i)
{
	Del_primitive(tab_primitives,i);
}

//-------------------------------------------
// Fusionne les primitives identiques
// issues d'intersection
//-------------------------------------------
void Cplan::Fusionne_primitives()
{
	Fusionne(tab_primitives);
}

//---------------------------------------------
// On scinde les polylignes discontinues
// dans le tableau des primitives
//---------------------------------------------
void Cplan::scinde_lignes()
{
	scinde_ligne(tab_primitives);
}

//-------------------------------------------
// Equation du plan connaissant le vecteur normal
// et un point du plan
//-------------------------------------------
void Cplan::Set_equation_plan()
{ 
	for (int i = 0 ; i < 3 ; i++)
		a[i] = T[i];
	a[3] = -Vector3D(pt0_.vector()).dot(T);

	bequation = true;
}

//-------------------------------------------
// Repere o.n.d. du plan (R,S,T) connaissant 
// le vecteur normal T
//-------------------------------------------
bool Cplan::Actualise_repere(bool cal)
{
	bool ok = true; 
	if (cal == true) ok = Set_vec_normal();
	if (ok == true)
	{
		// On initialise R et S (vecteurs nuls)
		R = Vector3D(Point3D().vector()); S = Vector3D(Point3D().vector());
		Cbase_primitive::Actualise_repere(R,S,T);
		Point3D R_pt(R);
		Point3D S_pt(S);

		if (R_pt == Point3D() || S_pt == Point3D()) ok = Set_RS();
	}
	return ok;
}

//-------------------------------------------
// Equation du plan
//-------------------------------------------
bool Cplan::Actualise_equation()
{	
	Point3D T_pt(T);

	if (T_pt == Point3D())
	{
		if (Actualise_repere(true) == false) 
			return false;
	}
	
	if (bequation == false) Set_equation_plan();
	return true;
}

//-------------------------------------------
// fonctions d'acces au variables
//-------------------------------------------
bool Cplan::Get_bequation() const { return bequation; }

Point3D& Cplan::Get_point_ref() { return pt0_; }
void Cplan::Set_point_ref(Point3D pt) { pt0_=pt; }

double Cplan::Get_coef(int i) const { return a[i]; }
void Cplan::Set_coef(double x,int i) { a[i] = x; }

Vector3D& Cplan::Get_R() { return R; }
Vector3D& Cplan::Get_S() { return S; }
Vector3D& Cplan::Get_vec_normal() { return T; }

bool Cplan::Set_vec_normal()
{ 
	// On calcule le vecteur normal a partir de 3 points du plan
	// (points donnes dans le polygone contour = tab_primitives[0])
	if (!tab_primitives.size() || 
		tab_primitives[0]->tab_point.size() < 3) return false;
	
	std::shared_ptr<Cpolygone> prim = std::dynamic_pointer_cast<Cpolygone>(tab_primitives[0]);
	// Points numerotes pour que la normale soit 
	// exterieure lorsqu'on tourne dans le sens direct
	Vector3D vn = prim->Get_vec_normal();

	pt0_=prim->tab_point[0];

	if (vn.norm() > 0.0)
	{
		T = vn;
		normalize(T);
		return true;
	}
	return false;
}

bool Cplan::Set_RS()
{ 
	if (!tab_primitives.size() || 
		tab_primitives[0]->tab_point.size() < 3)
	{
		// On calcule le vecteur R de maniere aleatoire
		Vector3D vecd = {0.,0.,0.};
		int			n = 0;
	
		if (T[2] != 0.) n = 2;
		else if (T[1] != 0.) n = 1;
		double	u = T[n];

		std::srand(std::time(NULL));
		for (int i = 0 ; i < 3 ; i++) 
		{
			if (i != n) vecd[i] = static_cast<double>(std::rand()) / RAND_MAX;
		}
		vecd[n] = -T.dot(Vector3D(vecd))/u;
		R = vecd;
	}
	else
	{
		// On calcule le vecteur R a partir du premier segment
		// du polygone contour = tab_primitives[0]
		std::shared_ptr<Cprimitive>	prim = tab_primitives[0];
		// Points numerotes pour que la normale soit 
		// exterieure lorsqu'on tourne dans le sens direct
		R = prim->tab_point[1].vector() - prim->tab_point[0].vector();
	}
	if (R.norm() == 0.) return false;
	normalize(R);
	S = T.cross(R);
	if (S.norm() == 0.) return false;
	return true;
}

bool Cplan::update_RS(Vector3D R_fixed)
{
	R = R_fixed;
	S = T.cross(R);
	if (S.norm() == 0.) return false;
	return true;
}




bool operator==(Cplan p1,Cplan p2)
{
	if (p1.R == p2.R && p1.S == p2.S && p1.T == p2.T)
		return true;
	return false;
}

bool operator!=(Cplan p1,Cplan p2)
{
	return !(p1 == p2);
}

//---------------------------------------------
// Transformation des coordonnees 3D d'un point 
// du plan en coordonnes R,S ds le plan
//---------------------------------------------
Point3D Cplan::Change_coord_2D(Point3D& pt)
{
	Vector3D vec = pt.vector()- pt0_.vector();
	return Point3D(vec.dot(R), vec.dot(S),0.0);
}

//---------------------------------------------
// Transformation des coordonnees 2D R,S du plan
// d'un point du plan en coordonnes x,y,z 3D
//---------------------------------------------

Point3D Cplan::Change_coord_3D(Point3D& pt)
{
	Point3D res_pt;

	res_pt.x()=R[0]*pt.x()+S[0]*pt.y();
	res_pt.y()=R[1]*pt.x()+S[1]*pt.y();
	res_pt.z()=R[2]*pt.x()+S[2]*pt.y();

	res_pt.vector()=res_pt.vector()+pt0_.vector();

	return res_pt;

}

//---------------------------------------------
// Projection d'un vecteur 3D
// dans le plan (R,S)
//---------------------------------------------
Vector3D Cplan::Projection_sur(Vector3D& vec)
{
	Vector3D res_vec;
	res_vec[0] = vec.dot(R);
	res_vec[1] = vec.dot(S);
	res_vec[2] = 0.0;

	return res_vec;
}

//---------------------------------------------
// fonction qui projecte orthogonalement
// le point pt sur le plan
// cf. bourke
//---------------------------------------------

Point3D Cplan::Projection_sur(Point3D& pt)
{

	double u=T.dot(pt0_.vector()- pt.vector());

	return Point3D(pt.vector() + T*u);

	
}

//---------------------------------------------
// Coordonnees 3D xmax et ymax et zmax
//---------------------------------------------
Point3D Cplan::Get_max_coord2D()
{
	Point3D	pt,pt_max(-1e20,-1e20,0.0);
	std::shared_ptr<Cprimitive> prim = tab_primitives[0]; // countour du plan

	for(size_t i=0;i<prim->tab_point.size();i++)
	{
		pt=Change_coord_2D(prim->tab_point[i]); // changement de coordonnees 3D->2D
		pt_max.x()=std::max(pt.x(),pt_max.x());
		pt_max.y()=std::max(pt.y(),pt_max.y());
	}
	return pt_max;
}
		
//---------------------------------------------
// Coordonnees 2D xmin et ymin
//---------------------------------------------
Point3D Cplan::Get_min_coord2D()
{
	Point3D	pt,pt_min(1e20,1e20,0.0);
	std::shared_ptr<Cprimitive> prim = tab_primitives[0]; // countour du plan

	for(size_t i=0;i<prim->tab_point.size();i++)
	{
		pt=Change_coord_2D(prim->tab_point[i]); // changement de coordonnees 3D->2D
		pt_min.x()=std::min(pt.x(),pt_min.x());
		pt_min.y()=std::min(pt.y(),pt_min.y());
	}
	return pt_min;
}

//-------------------------------------------
// Determination de la "bounding box" du plan
// (comprant toutes les primitives)
//-------------------------------------------	
void Cplan::Set_Bounding_box()
{
	
	if (boundingBox_.bounding_set_== false)
	{
		Cbase_primitive::Set_Bounding_box();

		int n = tab_primitives.size();
		if (!n) return;
		
		for (int i = 0 ; i < n ; i++)
		{
			tab_primitives[i]->Set_Bounding_box();
			if (i == 0)
			{
				boundingBox_ = tab_primitives[i]->get_boundingBox();
			}
			else
			{
				boundingBox_.add_bounding_box(tab_primitives[i]->get_boundingBox());
			}
		}

		boundingBox_.bounding_set_ = true;
	}
}

//-------------------------------------------
// Determination de la "bounding box" d'une
// primitive iprim du plan
//-------------------------------------------	
void Cplan::Set_Bounding_box(int iprim)
{
	int n = tab_primitives.size();
	if (!n || iprim > n) return;
	tab_primitives[iprim]->Set_Bounding_box();
}

//---------------------------------------------
// Test si un point est strictement dans le 
// polygone poly (=true)
// la limite du polygone est consideree
// comme externe
//---------------------------------------------
bool Cplan::Point_ds_polygone(Cprimitive* poly,Point3D& pt_3d)
{
	if (poly->Get_type() != Polygone) return false;
	
	// On verifie que le point appartient a la bounding box
	// du polygone
	poly->Set_Bounding_box();
	if (!poly->In_Bounding_box(pt_3d)) return false;
	
	Point3D pt = Change_coord_2D(pt_3d); // coordonnees planaires
	
	Point3D	pt1(poly->tab_point_2d_[0].vector() - pt.vector());
    Point3D	pt2;
	int			intersect = 0; // nombre d'intersection avec le rayon

	for (size_t i = 1 ; i < poly->tab_point_2d_.size() ; i++)
	{
		pt2.vector() = poly->tab_point_2d_[i].vector()- pt.vector();
	
		// pt appartient a la limite du polygone
		if (pt1.x() == 0.0 && pt1.y() == 0.0)
			return false;
		
		if((pt1.y() > 0.0) != (pt2.y() > 0.0))
		{
			double x = (pt1.x()*pt2.y() - pt2.x()*pt1.y())/(pt2.y() - pt1.y());
		    if (x > 0.0) intersect++;
		}
		pt1 = pt2;
	}
 	if ((intersect % 2) == 1) return true;
	return false;
}

//---------------------------------------------
// Test si un point est proche d'un plan
// (distance < tolerance donnee)
// cal = true : calcul d'intersections (par defaut)
//---------------------------------------------
bool Cplan::Point_ds_plan(Point3D& pt,double tol,bool cal)
{	
	if (Actualise_equation() == false) 
		return false;

	double dist = a[0]*pt.x() + a[1]*pt.y() + a[2]*pt.z()+ a[3];
	dist = std::abs(dist)/T.norm();

	if (dist < tol)
	{
		// le point est dans le plan
		// si le plan n'est pas limite (pas de polygone contour)
		// ou pas de calcul d'intersections, on renvoie true
		if (!tab_primitives.size() || !cal) return true;
		// sinon on verifie qu'il appartient au polygone contour
//		return Point_ds_polygone(tab_primitives[0],pt);
		return std::dynamic_pointer_cast<Cpolygone>(tab_primitives[0])->Point_ds_polygone(pt,tol,false);
	}
	return false;
}


//---------------------------------------------
// distance
//---------------------------------------------
double Cplan::distance_to_plan(Point3D& pt)
{
	double dist = 1e+10;

	if (Actualise_equation() == false) return dist;

	dist = a[0]*pt.x() + a[1]*pt.y() + a[2]*pt.z()+ a[3];
	dist = std::abs(dist)/T.norm();

	return dist ;

}

//---------------------------------------------
// distance
//---------------------------------------------
double Cplan::oriented_distance_to_plan(Point3D& pt)
{
	double dist = 1e+10;
	if (Actualise_equation() == false) return dist;

	return T.dot(pt.coord_ -pt0_.coord_);
}



//---------------------------------------------
// Test si deux plans sont confondus
//---------------------------------------------
bool Cplan::Plan_confondu_avec(Cplan *plan,double tol)
{
	Vector3D vec = T.cross(plan->Get_vec_normal());
	Point3D vec_pt(vec);
	if (vec_pt != Point3D()) return false;
	
	// Les plans sont paralleles
	if (Point_ds_plan(plan->Get_point_ref(),tol,false) == true)
	{
		// Les plans sont confondus (distance < tol)
		// On utilise la meme base (R,S,T) pour les 2 plans
		R = plan->Get_R();
		S = plan->Get_S();
		T = plan->Get_vec_normal();
		return true;
	}
	return false;
}

//---------------------------------------------
// Intersections plan/plan et plan/primitive 
// (polygone,polyligne,point)
//---------------------------------------------
void Cplan::intersect(Cbase_primitive *prim,int i0)
{
	if (!Intersect_Bounding_box(prim)) return;
	
	Tprimitive	ptype = prim->Get_type();

	if (ptype == Polyedre || ptype == Cylindre 
		|| ptype == Boite_Rect)
	{
		Cpolyedre* poly = (Cpolyedre*)(prim);
		for (size_t i = 0 ; i < poly->tab_plan.size() ; i++)
		{	
			intersect_plan(poly->tab_plan[i].get(),pt_tol);
		}
	}
	else if (ptype == Plan || ptype == Faille)
	{
		intersect_plan((Cplan*)prim,pt_tol);
	}
	else
	{
		Cprimitive*	pri = (Cprimitive*) prim;
		bool		inter = false;
		Tprimitive	ptypi;
		int			nprim = tab_primitives.size();

		for (int i = 0 ; i < nprim ; i++)
		{
			std::shared_ptr<Cprimitive> primi = tab_primitives[i];
			ptypi = primi->Get_type();
			// intersection d'un polygone avec une polyligne
			// et non le contraire
			if ((ptypi == Polyligne && ptype == Polygone)
				|| (ptypi == Point && ptype != Point))

				inter = pri->intersect(primi.get());
			
			else
				inter = primi->intersect(pri);
		}

		if (inter == false)
		{
			// Pas d'intersection detectee
			// on verifie que ce n'est pas un cas 
			// particulier (primitive incluse strictement
			// dans un polygone)
			if (pri->Primitive_incluse(tab_primitives[0].get()) == true)
			{
				int j = tab_primitives.size();
				switch (ptype)
				{
					case Polygone: {

						std::shared_ptr<Cpolygone> new_polygone = std::shared_ptr<Cpolygone>(new Cpolygone(*(Cpolygone*)prim));
						tab_primitives.push_back(new_polygone); }
					break;

					case Polyligne: 
						tab_primitives.push_back(std::shared_ptr<Cpolyligne>(new Cpolyligne(*(Cpolyligne*)prim)));
					break;

					case Point: 
						tab_primitives.push_back(std::shared_ptr<Cprimitive>(new Cprimitive(*pri)));
						break;
				}
				if (j != -1) tab_primitives[j]->Set_calcule(true);
			}
		}
		else if (!i0 && ptype == Polyligne)
		{
			std::shared_ptr<Cpolygone>	poly = std::dynamic_pointer_cast<Cpolygone>(tab_primitives[0]);
			Cpolyligne	*lin = (Cpolyligne*)prim;
			Point3D	p0 = lin->tab_point[0],p1,pt;
			
			for (size_t j = 1 ; j < lin->tab_point.size() ; j++)
			{
				p1 = lin->tab_point[j];
				Vector3D v01 = p1.vector() - p0.vector(); 

				if (std::abs(T.dot(v01)) < TOL)
				{
					// Le segment est coplanaire
					Cpolyligne pl(p0,p1);
					if (pl.Primitive_incluse(poly.get(),false,false) == true)
					{
							// on determine son point milieu
						pt.vector() = (p0.vector() + p1.vector())*0.5;
						if (poly->Point_ds_polygone(this,pt,pt_tol) == true)
						{
							pl.Set_calcule(true);
							tab_primitives.push_back(std::shared_ptr<Cpolyligne>(new Cpolyligne(pl)));
						}
					}
				}
				p0 = p1;
			}
		}
	
	}
}


//---------------------------------------------
// Intersection plan/plan uniquement
//---------------------------------------------
void Cplan::intersect_plan(Cplan *plan,double tol)
{
	double dist;

	if (!Intersect_Bounding_box(plan)) return;

	Point3D	p1,p0,pt;
	bool pas_intersection;

	if (Plan_confondu_avec(plan,tol) == true)
	{
		for (size_t i = 0 ; i < plan->tab_primitives.size() ; i++)
		{
			intersect(plan->tab_primitives[i].get());
		}

	}
	else
	{
		if (Actualise_equation() == false ||
			plan->Actualise_equation() == false) return;

		std::shared_ptr<Cpolygone> poly1 = std::shared_ptr<Cpolygone>(new Cpolygone(*std::dynamic_pointer_cast<Cpolygone>(tab_primitives[0]).get()));
		std::shared_ptr<Cpolygone> poly2 = std::shared_ptr<Cpolygone>(new Cpolygone(*std::dynamic_pointer_cast<Cpolygone>(plan->tab_primitives[0]).get()));

		if (!poly1->Intersect_Bounding_box(poly2.get())) return;

		pas_intersection=true;;

		dist=std::max(poly1->Distance_max(),poly2->Distance_max());

		double	a0 = plan->Get_coef(0),
				a1 = plan->Get_coef(1),
				a2 = plan->Get_coef(2),
				a3 = plan->Get_coef(3);
				
		// Vecteur directeur de la droite 3D
		Vector3D	vl = T.cross(plan->Get_vec_normal());

		normalize(vl);
		if (vl.norm() == 0.0) return;
		int n = poly1->tab_point.size();
		
		if (n > 2)
		{		
			// Polyligne issue de l'intersection entre la droite 3D
			// et le polygone contour de ce plan
			
			Point3D	p1,p0 = poly1->tab_point[0],pt;
			Vector3D	v;
			double		u = 0., den = 0.;
			
			p0 = poly1->tab_point[0];
			for (int i = 1 ; i < n ; i++)
			{
				p1 = poly1->tab_point[i];
				v = p1.vector() - p0.vector();
				if (!iscolinear(v,vl)) // lignes non parallelles
				{
					den = a0*v[0] + a1*v[1] + a2*v[2];
					u = -(a0*p0.x() + a1*p0.y() + a2*p0.z() + a3)/den;

					double u_tol=pt_tol/v.norm();
				

					if (u<0 && u>(-u_tol)) u=0.0;
					if (u>1.0 && (u-1.0)<u_tol) u=1.0;

					if (u >= 0 && u <= 1.0)
					{
						// Le point intersection entre le segment du
						// polygone poly1 et la droite 3D existe
						pt.vector() = p0.vector() + v*u;

						pas_intersection=false;

						break;

						// on enrichit la ligne
						//lin->tab_point.push_back(new Cpoint3D(pt,true));
					}
				}
				p0 = p1;
			}


			if (pas_intersection==true) return;


			std::shared_ptr<Cpolyligne> lin = std::shared_ptr<Cpolyligne>(new Cpolyligne);
			lin->Set_calcule(true);


			lin->tab_point.push_back(Point3D(pt.vector()+vl*(2.0*dist),true));
			lin->tab_point.push_back(Point3D(pt.vector()-vl*(2.0*dist),true));

			poly1->intersect(lin.get());
			poly2->intersect(lin.get());

			lin->Fusionne_points(pt_tol);

			lin->tab_point.erase(lin->tab_point.begin()+0);
			lin->tab_point.erase(lin->tab_point.begin() + lin->tab_point.size()-1);


			if (lin->tab_point.size() > 0)
			{

				n = lin->tab_point.size();
				
				bool	bpoly1 = false,bpoly2 = false;

				// On va enrichir les polygones contour poly1 et poly2 si necessaire
				if (n > 1)
				{
					lin->tab_connec.resize(n-1);
					p0 = lin->tab_point[0];
					for (int i = 1 ; i < n ; i++)
					{
						p1 = lin->tab_point[i];
						v = p1.vector() - p0.vector();
						// point milieu du segment
						pt.vector() = p0.vector() + v*0.5;
						bpoly1 = poly1->Point_ds_polygone(this,pt,pt_tol,false);
						bpoly2 = poly2->Point_ds_polygone(plan, pt,pt_tol,false);
						if (bpoly1 == true && bpoly2 == true)
						{
							lin->tab_connec[i-1] = true;
						}
						else
						{
							lin->tab_connec[i-1] = false;
							if (bpoly1 == true)
							{
								if (poly2->Point_ds_polyligne(p0,pt_tol))
								{

									Cprimitive point1(Point3D(p0,true));
									intersect(&point1,1);
								}
								if (poly2->Point_ds_polyligne(p1,pt_tol))
								{

									Cprimitive point1(Point3D(p1,true));
									intersect(&point1,1);
								}

							}
							if (bpoly2 == true)
							{
							if (poly1->Point_ds_polyligne(p0,pt_tol))
								{

									Cprimitive point1(Point3D(p0,true));
									plan->intersect(&point1,1);
								}
								if (poly1->Point_ds_polyligne(p1,pt_tol))
								{

									Cprimitive point1(Point3D(p1,true));
									plan->intersect(&point1,1);
								}
							}

						}
						p0 = p1;
					}
					Cpolyligne	pij;
					std::string		sz,szp = lin->Get_nom();
					int			k = 1,nij = 0;
					bool		cont = true;
		
					do
					{
						pij = Cpolyligne(sz);
						cont = lin->scinde_ligne(&pij);
						nij = pij.tab_point.size();
						if (nij > 0)
						{
                            pij.Set_Bounding_box();
							pij.Set_calcule(true);
							intersect(&pij,1);
							plan->intersect(&pij,1);
							if (nij > 1) 
							{
								tab_primitives.push_back(std::shared_ptr<Cpolyligne>(new Cpolyligne(pij)));
								plan->tab_primitives.push_back(std::shared_ptr<Cpolyligne>(new Cpolyligne(pij)));
							}
							
							if (cont == true) sz  = szp + boost::str(boost::format("%1%") % (++k) );
								
						}
					}
					while (cont == true);
					n = lin->tab_point.size();
				}
				
				if (n == 1)
				{
					p0 = lin->tab_point[0];
					bpoly1 = poly1->Point_ds_polygone(this,p0,pt_tol,false);
					bpoly2 = poly2->Point_ds_polygone(plan,p0,pt_tol,false);
					
					if (bpoly1 == true && bpoly2 == true)
					{
						// on enrichit eventuellement les polygones
						poly1->Point_ds_polyligne(p0,pt_tol,true);
						poly2->Point_ds_polyligne(p0,pt_tol,true);

						Cprimitive	point1(lin->tab_point[0]),
									point2(lin->tab_point[0]);
					
						intersect(&point1,1);
						plan->intersect(&point2,1);
					}
				}
				else if (n > 1)
				{
					std::shared_ptr<Cpolyligne> lin1 = std::shared_ptr<Cpolyligne>(new Cpolyligne(*lin));
								
					lin1->Set_Bounding_box();

					intersect(lin1.get());
					plan->intersect(lin1.get());
					
					tab_primitives.push_back(std::shared_ptr<Cpolyligne>(new Cpolyligne(*lin1.get())));
					plan->tab_primitives.push_back(std::shared_ptr<Cpolyligne>(new Cpolyligne(*lin1.get())));


				}
			}
		}
	}
}

//-------------------------------------------
// Elimination des primitives issues des
// intersections
//-------------------------------------------	
void Cplan::Del_intersect()
{
	int n = tab_primitives.size();
	if (!n) return;
	
	for (int i = n-1 ; i >= 0 ; i--)
	{
		if (!tab_primitives[i]->tab_point.size())
			Del_primitives(i);
		
		else if (tab_primitives[i] != NULL)
		{
			if (tab_primitives[i]->Get_calcule() == true)
			{
				Del_primitives(i);
			}
			else
			{
				tab_primitives[i]->Del_intersect();
				if (!tab_primitives[i]->tab_point.size())
					Del_primitives(i);
			}
		}
	}
}




//---------------------------------------------
// On regarde dans le tableau des primitives
// si des polylignes sont discontinues
// si oui, on les scinde pour ne garder que 
// des polylignes continues
//---------------------------------------------
void scinde_ligne(std::vector<std::shared_ptr<Cprimitive> >& tab)
{
	if (!tab.size()) return;

	for (int i = tab.size()-1 ; i > 0 ; i--)
	{
		if (tab[i]->Get_type() == Polyligne)
		{
			std::shared_ptr<Cpolyligne>	pi = std::dynamic_pointer_cast<Cpolyligne>(tab[i]);
			std::string		sz,szp = pi->Get_nom();
			int			j = 1;
			bool		cont = true;
			
			do
			{
				Cpolyligne pij = Cpolyligne(sz);
				cont = pi->scinde_ligne(&pij);
				if (pij.tab_point.size() > 0)
				{
					pij.Set_calcule(true);
					tab.push_back(std::shared_ptr<Cprimitive>(new Cpolyligne(pij)));
					if (cont == true) sz  = szp + boost::str(boost::format("%1%") % (++j) );
				}
			}
			while (cont == true);
		}
		if (tab[i]->tab_point.empty())
		{
			Del_primitive(tab,i);
		}
	}
}

//-------------------------------------------
// Elimination d'une primitive du tableau
// des primitives
//-------------------------------------------	
void Del_primitive(std::vector<std::shared_ptr<Cprimitive>>& tab,int i)
{
	int n = tab.size();
	if (!n || i > n) return;
	
	if (i >= 0)
	{
		tab.erase(tab.begin() + i);
	}
	else if (i == -1)
	{
		tab.clear();
	}
}

//-------------------------------------------
// Elimination d'un plan du tableau
// des plans
//-------------------------------------------	
void Del_plan(std::vector<std::shared_ptr<Cplan>>& tab,int i)
{
	int n = tab.size();
	if (!n || i > n) return;
	
	if (i >= 0)
	{

		tab.erase(tab.begin() + i);
	}
	else if (i == -1)
	{
		tab.clear();
	}
}

//-------------------------------------------
// Fusionne les primitives identiques ou 
// Connecte ou supprime les polylignes
// redondantes issues d'intersection
//-------------------------------------------
void Fusionne(std::vector<std::shared_ptr<Cprimitive>>& tab)
{

	int	n,np,np1;
	Point3D p1,p2;
	bool kill;

	// on supprime pour chaque primitive les points doubles cr�es 
    // lors des intersections

	for (size_t i = 1 ; i < tab.size() ; i++)
	{
		np=tab[i]->tab_point.size();

		if (np>1)  
		{
			tab[i]->Fusionne_points(pt_tol);

        // si un polygone ou une polyligne c'est transforme en un point 

			if (tab[i]->tab_point.size()==1)
			{
					
				tab.push_back(std::shared_ptr<Cprimitive>(new Cprimitive(*tab[i])));
				tab[i]->Set_type(Point);
				tab.erase(tab.begin() + i);
			}
		}
		else if(np==1)
		{
			if(tab[i]->Get_calcule()==true)
			{
				p1= tab[i]->tab_point[0];
				kill=false;
				for(size_t j =0 ; j < tab.size(); j++)
				{
					if (i!=j)
					{
						if(tab[j]->Find_point(p1,0,pt_tol))
						{
							kill=true;
							break;
						}
					}
				}
				if(kill==true) 
				{
					tab.erase(tab.begin() + i);
					i--;
				}
			}
		}

	}


	for( size_t i = 0 ; i < tab.size(); i++)
	{
		np=tab[i]->tab_point.size();
		
		if(np>1) 
		{
			for(size_t j = i+1 ; j < tab.size(); j++)
			{
				kill=false;
				
				if(tab[j]->Get_calcule()==true)
				{
					np1=tab[j]->tab_point.size();
					
					if (np1==2) 
					{
						p1=tab[j]->tab_point[0];
						p2=tab[j]->tab_point[1];
						if(tab[i]->Segment_inclus(p1,p2))
						{		
							kill=true;	
						}
					}
						
				}
				if(kill==true) 
				{
					tab.erase(tab.begin() + j);
					j--;
				}
			}

		}

	}

		
	n=tab.size();

	if (n <= 1) return;
	
	bool	*flag = new bool[n];
	for (int i = 0 ; i < n ; i++) flag[i] = false;

	for (int i = 1 ; i < n ; i++)
	{
		if (!flag[i])
		{
			if (tab[i]->Get_type() == Polyligne)
			{
				int nn = tab[i]->tab_point.size();
				if (!nn) flag[i] = true;
				else if (tab[i]->tab_point[0] == tab[i]->tab_point[nn-1])
				{
					flag[i] = true;
					int j = tab.size();
					tab.push_back(std::shared_ptr<Cpolygone>(new Cpolygone(*(std::dynamic_pointer_cast<Cpolyligne>(tab[i])))));
					tab[j]->Set_calcule(true);
					std::string sz  = boost::str(boost::format("polygone_fusion %1% ") % (i) );
					tab[j]->Set_nom(sz);
				}
			}
		}
	}
	

	for (int i = n-1 ; i >= 0 ; i--)
	{
		if (flag[i] == true)
		{
			Del_primitive(tab,i);
		}
	}
	delete [] flag;


}

//---------------------------------------------
// Intersection plan/plan with surface limitation
//---------------------------------------------
bool Cplan::intersect_plan_limitated(Cplan *plan,double tol) {


	if (!Intersect_Bounding_box(plan)) return false;

	bool pas_intersection;

	if (Plan_confondu_avec(plan,tol) == true) return false;

	if (Actualise_equation() == false ||
			plan->Actualise_equation() == false) return false; ;

	std::shared_ptr<Cpolygone> poly1 = std::shared_ptr<Cpolygone>(new Cpolygone(*std::dynamic_pointer_cast<Cpolygone>(tab_primitives[0]).get()));
	std::shared_ptr<Cpolygone> poly2 = std::shared_ptr<Cpolygone>(std::dynamic_pointer_cast<Cpolygone>(plan->tab_primitives[0]));

	pas_intersection=true;;

	double dist= std::max(poly1->Distance_max(), poly2->Distance_max());

	double	a0 = plan->Get_coef(0),
			a1 = plan->Get_coef(1),
			a2 = plan->Get_coef(2),
			a3 = plan->Get_coef(3);

	Vector3D	vl = T.cross(plan->Get_vec_normal());

	normalize(vl);
	if (vl.norm() == 0.0) return false;
	int n = poly1->tab_point.size();

	if (n < 3) return false;

	// Polyligne issue de l'intersection entre la droite 3D
	// et le polygone contour de ce plan

	Vector3D	v;
	double		u = 0., den = 0.;

	Point3D pt;

	Point3D p0 = poly1->tab_point[0];
	for (int i = 1 ; i < n ; i++)
	{
		Point3D p1 = poly1->tab_point[i];
		v = p1.vector() - p0.vector();
		if (!iscolinear(v,vl)) // lignes non parallelles
		{
			den = a0*v[0] + a1*v[1] + a2*v[2];
			u = -(a0*p0.x() + a1*p0.y() + a2*p0.z() + a3)/den;

			double u_tol=pt_tol/v.norm();

			if (u<0 && u>(-u_tol)) u=0.0;
			if (u>1.0 && (u-1.0)<u_tol) u=1.0;

			if (u >= 0 && u <= 1.0)
			{
				// Le point intersection entre le segment du
				// polygone poly1 et la droite 3D existe
				pt.vector() = p0.vector() + v*u;

				pas_intersection=false;
				break;
			}
		}

		p0 = p1;
	}


	if(pas_intersection==true) return false;



	std::shared_ptr<Cpolyligne>	lin2 = std::shared_ptr<Cpolyligne>(new Cpolyligne);
	lin2->Set_calcule(true);

	lin2->tab_point.push_back(Point3D(pt.vector()+vl*(2.0*dist),true));
	lin2->tab_point.push_back(Point3D(pt.vector()-vl*(2.0*dist),true));

	poly2->intersect(lin2.get(), false);

	lin2->Fusionne_points(pt_tol);

	lin2->tab_point.erase(lin2->tab_point.begin());
	lin2->tab_point.erase(lin2->tab_point.begin() + lin2->tab_point.size()-1);

	std::shared_ptr<Cpolyligne> lin1 = std::shared_ptr<Cpolyligne>(new Cpolyligne);
	lin1->Set_calcule(true);

	lin1->tab_point.push_back(Point3D(pt.vector()+vl*(2.0*dist),true));
	lin1->tab_point.push_back(Point3D(pt.vector()-vl*(2.0*dist),true));

	poly1->intersect(lin1.get());

	lin1->Fusionne_points(pt_tol);

	lin1->tab_point.erase(lin1->tab_point.begin());
	lin1->tab_point.erase(lin1->tab_point.begin() + lin1->tab_point.size()-1);


	if (lin2->tab_point.empty()) return false;


	poly1->update_point2D(this);
	int count = 0;
	for (size_t i = 0; i < lin2->tab_point.size(); i++)
	{

		if (poly1->Point_ds_polygone(this,lin2->tab_point[i],pt_tol,false)) {
			count++;
		}

	}
	for (size_t i = 0; i < lin1->tab_point.size(); i++)
	{

		if (poly2->Point_ds_polygone(plan,lin1->tab_point[i],pt_tol,false)) {
			count++;
		}

	}

	if (count < 2) return false;



	if (lin1->tab_point.size() == 2)
	{
		std::shared_ptr<Cpolygone> new_poly1 =std::shared_ptr<Cpolygone>(new Cpolygone("poly1"));
		std::shared_ptr<Cpolygone> new_poly2 =std::shared_ptr<Cpolygone>(new Cpolygone("poly2"));

		Point3D p0 = lin1->tab_point[0];
		Point3D p1 = lin1->tab_point[1];

		int pos_pt_0 = poly1->Find_point(p0,1,pt_tol);
		int pos_pt_1 = poly1->Find_point(p1,1,pt_tol);

		int min_pos = std::min(pos_pt_0,pos_pt_1);
		int max_pos = std::max(pos_pt_0,pos_pt_1);

        int count_1 = 0;
		for (int pli = 0; pli < min_pos+1; pli++)
		{
			new_poly1->tab_point.push_back(poly1->tab_point[pli]);
			count_1++;
		}
		new_poly1->tab_point[min_pos].calculated_ = true;

		int count_2 = 0;
		for (int pli = min_pos; pli < max_pos+1; pli++)
		{
			new_poly2->tab_point.push_back(poly1->tab_point[pli]);
			count_2++;
		}
		new_poly2->tab_point[0].calculated_ = true;
		new_poly2->tab_point[count_2-1].calculated_ = true;

		new_poly2->tab_point.push_back(poly1->tab_point[min_pos]);
		new_poly2->tab_point[count_2].calculated_ = true;

		for (size_t pli = max_pos; pli < poly1->tab_point.size(); pli++)
		{
			new_poly1->tab_point.push_back(poly1->tab_point[pli]);
		}

		new_poly1->tab_point[count_1].calculated_ = true;

		new_poly1->update_point2D(this);

		if (new_poly1->Point_ds_polygone(this, this->pt0_,pt_tol,false)) {
			tab_primitives[0] = new_poly1;
		}
		else {
			tab_primitives[0] = new_poly2;
			new_poly2->update_point2D(this);
		}



		return true;
	}

	return false;
}

//---------------------------------------------
// Intersection plan/plan with surface limitation
//---------------------------------------------
bool Cplan::intersect_plan_limitated_dist(Cplan *plan,double tol, double& dist_min)
{

	if (!Intersect_Bounding_box(plan)) return false;

	bool pas_intersection;

	if (Plan_confondu_avec(plan,tol) == true) return false;

	if (Actualise_equation() == false ||
			plan->Actualise_equation() == false) return false; ;

	std::shared_ptr<Cpolygone> poly1 = std::shared_ptr<Cpolygone>(std::dynamic_pointer_cast<Cpolygone>(tab_primitives[0]));
	std::shared_ptr<Cpolygone> poly2 = std::shared_ptr<Cpolygone>(std::dynamic_pointer_cast<Cpolygone>(plan->tab_primitives[0]));

	pas_intersection=true;;

	double dist= std::max(poly1->Distance_max(), poly2->Distance_max());

	double	a0 = plan->Get_coef(0),
			a1 = plan->Get_coef(1),
			a2 = plan->Get_coef(2),
			a3 = plan->Get_coef(3);

	Vector3D	vl = T.cross(plan->Get_vec_normal());

	normalize(vl);
	if (vl.norm() == 0.0) return false;
	int n = poly1->tab_point.size();

	if (n < 3) return false;

	// Polyligne issue de l'intersection entre la droite 3D
	// et le polygone contour de ce plan

	Vector3D	v;
	double		u = 0., den = 0.;

	Point3D pt;

	Point3D p0 = poly1->tab_point[0];
	for (int i = 1 ; i < n ; i++)
	{
		Point3D p1 = poly1->tab_point[i];
		v = p1.vector() - p0.vector();
		if (!iscolinear(v,vl)) // lignes non parallelles
		{
			den = a0*v[0] + a1*v[1] + a2*v[2];
			u = -(a0*p0.x() + a1*p0.y() + a2*p0.z() + a3)/den;

			double u_tol=pt_tol/v.norm();

			if (u<0 && u>(-u_tol)) u=0.0;
			if (u>1.0 && (u-1.0)<u_tol) u=1.0;

			if (u >= 0 && u <= 1.0)
			{
				// Le point intersection entre le segment du
				// polygone poly1 et la droite 3D existe
				pt.vector() = p0.vector() + v*u;

				pas_intersection=false;
				break;
			}
		}

		p0 = p1;
	}


	if(pas_intersection==true) return false;


	std::shared_ptr<Cpolyligne>	lin2 = std::shared_ptr<Cpolyligne>(new Cpolyligne);
	lin2->Set_calcule(true);

	lin2->tab_point.push_back(Point3D(pt.vector()+vl*(2.0*dist),true));
	lin2->tab_point.push_back(Point3D(pt.vector()-vl*(2.0*dist),true));

	poly2->intersect(lin2.get(), false);

	lin2->Fusionne_points(pt_tol);

	lin2->tab_point.erase(lin2->tab_point.begin()+0);
	lin2->tab_point.erase(lin2->tab_point.begin() + lin2->tab_point.size()-1);

	std::shared_ptr<Cpolyligne> lin1 = std::shared_ptr<Cpolyligne>(new Cpolyligne);
	lin1->Set_calcule(true);

	lin1->tab_point.push_back(Point3D(pt.vector()+vl*(2.0*dist),true));
	lin1->tab_point.push_back(Point3D(pt.vector()-vl*(2.0*dist),true));

	poly1->intersect(lin1.get(),false);

	lin1->Fusionne_points(pt_tol);

	lin1->tab_point.erase(lin1->tab_point.begin()+0);
	lin1->tab_point.erase(lin1->tab_point.begin() + lin1->tab_point.size()-1);


	if (lin2->tab_point.empty()) return false;

	int count = 0;
	for (size_t i = 0; i < lin2->tab_point.size(); i++)
	{

		if (poly1->Point_ds_polygone(this,lin2->tab_point[i],pt_tol,false)) {
			count++;
		}

	}
	for (size_t i = 0; i < lin1->tab_point.size(); i++)
	{

		if (poly2->Point_ds_polygone(plan,lin1->tab_point[i],pt_tol,false)) {
			count++;
		}

	}

	if (count< 2) return false;


	if (lin1->tab_point.size() == 2)
	{
		Point3D pt_center;
		pt_center.coord_ = (lin1->tab_point[0].coord_+lin1->tab_point[1].coord_)/2.0;

		dist_min = distance(pt_center, pt0_);

		return true;
	}

	return false;
}


//---------------------------------------------
// Intersection plan/plan uniquement
//---------------------------------------------
bool Cplan::check_intersect_plan(Cplan *plan,double tol)
{
	double dist;

	if (Intersect_Bounding_box(plan) == false) return false;

	Point3D	p1,p0,pt;
	bool pas_intersection;

	if (Plan_confondu_avec(plan,tol) == true) return false;

	if (Actualise_equation() == false ||
		plan->Actualise_equation() == false) return false;

	std::shared_ptr<Cpolygone> poly1 = std::shared_ptr<Cpolygone>(std::dynamic_pointer_cast<Cpolygone>(tab_primitives[0]));
	std::shared_ptr<Cpolygone> poly2 = std::shared_ptr<Cpolygone>(std::dynamic_pointer_cast<Cpolygone>(plan->tab_primitives[0]));


	if (!poly1->Intersect_Bounding_box(poly2.get())) return false;

	pas_intersection=true;;

	dist=std::max(poly1->Distance_max(),poly2->Distance_max());

	double	a0 = plan->Get_coef(0),
			a1 = plan->Get_coef(1),
			a2 = plan->Get_coef(2),
			a3 = plan->Get_coef(3);

	// Vecteur directeur de la droite 3D
	Vector3D	vl = T.cross(plan->Get_vec_normal());

	normalize(vl);
	if (vl.norm() == 0.0) return false;
	int n = poly1->tab_point.size();

	if (n < 3) return false;

	Vector3D	v;
	double		u = 0., den = 0.;

	p0 = poly1->tab_point[0];
	for (int i = 1 ; i < n ; i++)
	{
		p1 = poly1->tab_point[i];
		v = p1.vector() - p0.vector();
		if (!iscolinear(v,vl)) // lignes non parallelles
		{
			den = a0*v[0] + a1*v[1] + a2*v[2];
			u = -(a0*p0.x() + a1*p0.y() + a2*p0.z() + a3)/den;

			double u_tol=pt_tol/v.norm();

			if (u<0 && u>(-u_tol)) u=0.0;
			if (u>1.0 && (u-1.0)<u_tol) u=1.0;

			if (u >= 0 && u <= 1.0)
			{
				// Le point intersection entre le segment du
				// polygone poly1 et la droite 3D existe
				pt.vector() = p0.vector() + v*u;

				pas_intersection=false;

				break;

				// on enrichit la ligne
				//lin->tab_point.push_back(new Cpoint3D(pt,true));
			}
		}
			p0 = p1;
	}


	if(pas_intersection==true) return false;


	std::shared_ptr<Cpolyligne>	lin = std::shared_ptr<Cpolyligne>(new Cpolyligne);
	lin->Set_calcule(true);


	lin->tab_point.push_back(Point3D(pt.vector()+vl*(2.0*dist),true));
	lin->tab_point.push_back(Point3D(pt.vector()-vl*(2.0*dist),true));

	poly1->intersect(lin.get(),false);

	n = lin->tab_point.size();

	poly2->intersect(lin.get(),false);

	n = lin->tab_point.size();

	lin->Fusionne_points(pt_tol);

	lin->tab_point.erase(lin->tab_point.begin()+0);
	lin->tab_point.erase(lin->tab_point.begin() + lin->tab_point.size()-1);

	int count = 0;

	n = lin->tab_point.size();

	p0 = lin->tab_point[0];
	for (int i = 1 ; i < n ; i++)
	{
		p1 = lin->tab_point[i];
		v = p1.vector() - p0.vector();
						// point milieu du segment
		pt.vector() = p0.vector() + v*0.5;
		bool bpoly1 = poly1->Point_ds_polygone(this,pt,pt_tol,false);
		bool bpoly2 = poly2->Point_ds_polygone(plan,pt,pt_tol,false);
		if (bpoly1 == true && bpoly2 == true) count++;
		p0 = p1;
	}


    return  (count > 0);

}


//-----------------------------------------------------
// write operator for plan mesh triangle
//-----------------------------------------------------
void Cplan::add_mesh_data(std::vector<Point3D>& tab_point, std::vector<std::vector<int>>& tab_triangle, std::vector<std::vector<float>>& tab_triangle_fracture_index,
		std::vector<float>& property)
{

	int index = tab_point.size();

	int node_size = tab_primitives[0]->tab_point.size()-1;

	Point3D center;

	for (int i = 0; i < node_size; i++)
	{
		tab_point.push_back(tab_primitives[0]->tab_point[i]);
		center.coord_ += tab_primitives[0]->tab_point[i].coord_;
	}

	if (node_size > 0) center.coord_ = center.coord_ /node_size;

	if (node_size == 3) {

		std::vector<int> triangle1 = {0+index, 1 + index , 2 + index};
		tab_triangle.push_back(triangle1);
		tab_triangle_fracture_index.push_back(property);
	}
	else if (node_size == 4) {
		std::vector<int> triangle1 = {0+index, 1 + index , 2 + index};
		tab_triangle.push_back(triangle1);
		tab_triangle_fracture_index.push_back(property);
		std::vector<int> triangle2 = {0+index, 2 + index , 3 + index};
		tab_triangle.push_back(triangle2);
		tab_triangle_fracture_index.push_back(property);
	}
	else {
		tab_point.push_back(center.coord_);

		for (int i = 0; i < node_size; i++)
		{
			int ni = i+1;
			if (ni == node_size) ni = 0;
			std::vector<int> triangle1 = {i+index, ni + index , node_size + index};
			tab_triangle_fracture_index.push_back(property);
			tab_triangle.push_back(triangle1);

		}

	}

}


//-----------------------------------------------------
// write operator for plan mesh triangle
//-----------------------------------------------------
std::ostream& operator<<(std::ostream& os, Cplan* plan)
{
	int node_size = plan->tab_primitives[0]->tab_point.size()-1;

	os<<std::string("TFACE")<<std::endl;

	Point3D center;

	for (int i = 0; i < node_size; i++)
	{
		os<<std::string("PVRTX ")<< i <<" " << plan->tab_primitives[0]->tab_point[i]<<std::endl;
		center.coord_ += plan->tab_primitives[0]->tab_point[i].coord_;
	}

	if (node_size > 0) center.coord_ = center.coord_ /node_size;

	if (node_size == 3) {

		os<<std::string("TRGL");
		os<<" " << " 0 1 2" <<std::endl;

	}
	else if (node_size == 4) {
		os<<std::string("TRGL");
		os<<" " << " 0 1 2" <<std::endl;
		os<<std::string("TRGL");
		os<<" " << " 0 2 3" <<std::endl;
	}
	else {
		os<<std::string("PVRTX ")<< node_size <<" " << center <<std::endl;

		for (int i = 0; i < node_size; i++)
		{
			os<<std::string("TRGL");
			int ni = i+1;
			if (ni == node_size) ni = 0;
			os<<" " << i <<" "<< ni <<" " << node_size <<std::endl;
		}

	}


	return os;

}



}





