
#include <primitive_mesh/primitive.h>
#include <primitive_mesh/boite.h>
#include <string>
#include <cmath>

namespace fractures_intersect {


extern double precis,pt_tol,TOL;
//-------------------------------------------
// Version 1.0 Pascal SIEGEL - 2000
//-------------------------------------------

//-------------------------------------------
// classe Cbase_primitive
//-------------------------------------------

//------------------------------------------------------
// construction/destruction de la classe Cbase_primitive
//------------------------------------------------------
Cbase_primitive::Cbase_primitive(std::string sz, Tprimitive nv_type,bool cal):
nom(sz),type(nv_type),calcule(cal)
{
	Selection = false;
	taille_element=1000;
	ind_materiaux=0;
	active = true;
}

Cbase_primitive::Cbase_primitive(const Cbase_primitive& prim)
{
	Selection = prim.Selection;
	boundingBox_ = prim.boundingBox_;
	taille_element = prim.taille_element;
	ind_materiaux = prim.ind_materiaux; 
	type = prim.type;
	if (prim.nom.empty() == false) nom = prim.nom;
	calcule = prim.calcule;
	active = prim.active;


}

Cbase_primitive::~Cbase_primitive() 
{
}


void Cbase_primitive::operator =(Cbase_primitive prim)
{
	Selection = prim.Selection;
	boundingBox_ = prim.boundingBox_;
	taille_element = prim.taille_element;
	ind_materiaux = prim.ind_materiaux; 
	type = prim.type;
	if (prim.nom.empty() == false) nom = prim.nom;
	calcule = prim.calcule;
	active = prim.active;

}

//--------------------------------------------------------
//  Attributs de la classe Cbase_primitive
//--------------------------------------------------------
Tprimitive Cbase_primitive::Get_type() const
{ return type; }

void Cbase_primitive::Set_type(Tprimitive nv_type)
{ 
	type = nv_type; 
}


std::string Cbase_primitive::Get_nom() const { return nom; }

void Cbase_primitive::Set_nom(std::string  nv_nom)
{
	if (!nv_nom.empty()) nom = nv_nom;
}


bounding_box& Cbase_primitive::get_boundingBox() { return boundingBox_; }

//-----------------------------------------------------------
// Indice materiau du polyedre
//----------------------------------------------------------- 
int Cbase_primitive::Get_ind_materiaux() const 
{ return ind_materiaux; }

void Cbase_primitive::Set_ind_materiaux(int nv_ind)  
{ ind_materiaux=nv_ind; }

//--------------------------------------------------------
//  calcule = true si une primitive de base est issue des
//  calculs d'intersection et sinon false
//--------------------------------------------------------

bool Cbase_primitive::Get_calcule() const { return calcule; }
void Cbase_primitive::Set_calcule(bool cal) { calcule = cal; }

//--------------------------------------------------------
//  Etat de selection de la primitive  
//--------------------------------------------------------
void Cbase_primitive::Selectionne_Moi(bool sel) { Selection = sel; }
bool Cbase_primitive::Selectionne_Moi() const {	return Selection; }

//-------------------------------------------
// Operations WARNING Pas implement�
//-------------------------------------------	
void Cbase_primitive::Set_taille_element(double taille) 
{
	taille_element = taille;
}

	double Cbase_primitive::Get_taille_element() const
{ return taille_element; }



//-------------------------------------------
// Determination d'un repere o.n.d (e0,e1,e2)
// connaissant e2
//-------------------------------------------	
void Cbase_primitive::Actualise_repere(Vector3D& e0,Vector3D& e1,
									Vector3D& e2)
{
	// Cette fonction est utile si la base o.n.d cherchee est 
	// suivant les axes du repere global (x,y,z)
	Vector3D x(1,0,0);
	Vector3D y(0,1,0);
	Vector3D z(0,0,1);

	double sign = 0.;

	if (e2.cross(x).norm() == 0.)
	{
		sign = (e2.dot(x) > 0. ? 1 : -1);
		// e2 = x*sign
		e0 = y*sign;
		e1 = z;
	}
	
	else if (e2.cross(y).norm() == 0.)
	{
		sign = (e2.dot(y) > 0. ? 1 : -1);
		// e2 = y*sign
		e0 = z*sign;
		e1 = x;
	}

	else if (e2.cross(z).norm() == 0.)
	{
		sign = (e2.dot(z) > 0. ? 1 : -1);
		// e2 = z*sign
		e0 = x*sign;
		e1 = y;
	}
}

//-------------------------------------------
// Determination de la boite englobante de la primitive
// fonction virtuelle implementee ulterieurement
//-------------------------------------------	
void Cbase_primitive::Set_Bounding_box() 
{
}

//-------------------------------------------
// Fonction verifiant si un point est dans la
// boite englobante de la primitive
//-------------------------------------------	
bool Cbase_primitive::In_Bounding_box(Point3D& pt)
{
	Set_Bounding_box();
	return boundingBox_.contains(pt);
}

//-------------------------------------------
// Fonction verifiant si une autre primitive de base
// est dans la boite englobante de celle-ci
//-------------------------------------------	
bool Cbase_primitive::In_Bounding_box(Cbase_primitive *p)
{

	Set_Bounding_box();
	p->Set_Bounding_box();

	return boundingBox_.contains(p->get_boundingBox());
}

//-------------------------------------------
// Fonction verifiant si les boites englobantes
// des deux primitives s'intersectent
//-------------------------------------------	
bool Cbase_primitive::Intersect_Bounding_box(Cbase_primitive *p)
{


	Set_Bounding_box();
	p->Set_Bounding_box();
	
	return boundingBox_.intersectionAndContains(p->get_boundingBox());
}

//-------------------------------------------
// Fonction calculant la plus grande distance 
// separant deux points de la primitive
//-------------------------------------------
double Cbase_primitive::Distance_max()
{
	Set_Bounding_box();
	return distance(boundingBox_.pt_min_,boundingBox_.pt_max_);
}


//-------------------------------------------
// classe Cprimitive : public Cbase_primitive
// coordonn�es du noeud + variables associ�es 
//-------------------------------------------

//-------------------------------------------------
// construction/destruction de la classe Cprimitive
//-------------------------------------------------

Cprimitive::Cprimitive(std::string sz,Tprimitive nv_type,bool cal):
Cbase_primitive(sz,nv_type,cal)
{}

Cprimitive::Cprimitive(const Cprimitive& prim):
Cbase_primitive((const Cbase_primitive&)prim)
{
	for(size_t i = 0 ; i < prim.tab_point.size() ; i++)
	{
		tab_point.push_back(prim.tab_point[i]);
	}
}

Cprimitive::Cprimitive(Point3D pt): Cbase_primitive()
{
	tab_point.push_back(pt);
	calcule = pt.calculated_;
}

Cprimitive::~Cprimitive()
{
}

Cprimitive& Cprimitive::operator=(const Cprimitive& prim)
{
	Selection = prim.Selection;
	taille_element = prim.taille_element;
	ind_materiaux = prim.ind_materiaux; 
	type = prim.type;
	if (prim.nom.empty() == false) nom = prim.nom;
	calcule = prim.calcule;

	Del_point();
	for(size_t i = 0 ; i < prim.tab_point.size() ; i++)
	{
		tab_point.push_back(prim.tab_point[i]);
	}

	tab_point_2d_.clear();
	for(size_t i = 0 ; i < prim.tab_point_2d_.size() ; i++)
	{
		tab_point_2d_.push_back(prim.tab_point_2d_[i]);
	}


	return *this;
}

//-------------------------------------------
// Elimination d'un point dans le tableau des points
//-------------------------------------------	
void Cprimitive::Del_point(int i)
{
	int n = tab_point.size();
	if (!n || i > n) return;
	
	if (i >= 0)
	{
		tab_point.erase(tab_point.begin() + i);
	}
	else if (i == -1)
	{
		tab_point.clear();
	}
}

//-------------------------------------------
// Recherche d'un point dans le tableau des points
// (a tol pres et a partir de l'indice i0)
//-------------------------------------------	
int Cprimitive::Find_point(Point3D pt,int i0,double tol)
{
	int	n = tab_point.size();
	if (!n) return -1;

	Point3D pti;

	for (int i = i0 ; i < n ; i++)
	{
		pti = tab_point[i];
		if (pti == pt || distance(pti,pt) < tol) 
			return i;
	}
	return -1;
}

//-------------------------------------------
// La distance de fusionnement des points est definie
// par tol (par defaut, on peut prendre tol = pt_tol)
//-------------------------------------------	
void Cprimitive::Fusionne_points(double tol)
{
	
	if (tab_point.size()<1) return;
	
	int		n = 0;
	Point3D pt0,pt1;
	Cpolyligne	*p = NULL;

	if (type == Polyligne)
	{
		p = (Cpolyligne*)this;
	}

	// Pour un polygone, le premier et le dernier
	// points sont identiques (polyligne fermee)
	if (type == Polygone) n=1;
	

	pt0=tab_point[0];

	for(size_t i=1;i<tab_point.size()-n;i++)
	{
			pt1=tab_point[i];
			
			if (tab_point[i].calculated_==true && pt1==pt0)
			{
				Del_point(i);
				p->Del_connec(i-1);
			}
			else pt0=pt1;
	}

}

//-------------------------------------------
// Determination de la "bounding box" de la primitive
//-------------------------------------------	
void Cprimitive::Set_Bounding_box()
{
	
	if (boundingBox_.bounding_set_== false)
	{
		Cbase_primitive::Set_Bounding_box();

		int n = tab_point.size();
		
		if (!n) return;
		if (type == Polygone) n--;

		boundingBox_.set_point(tab_point[0]);
		
		Point3D pti,pt=Point3D(pt_tol,pt_tol,pt_tol);

		for (int i = 1 ; i < n ; i++)
		{
			boundingBox_.add_point(tab_point[i]);
		}

		Eigen::Matrix<double, 3, 1> vec = boundingBox_.pt_max_.coord_-boundingBox_.pt_min_.coord_;
		double max_size =vec.maxCoeff();

		pt.coord_ *= max_size;

		boundingBox_.extend(pt.coord_);
		

		boundingBox_.bounding_set_ = true;
	}

}

//-------------------------------------------
// Fonction testant si une primitive (point,
// polyligne ou polygone) est incluse
// dans un polygone
//-------------------------------------------	
bool Cprimitive::Primitive_incluse(Cprimitive *prim,bool strict,bool)
{
	if (type != Point || prim->Get_type() != Polygone)
		return false;
	
	return ((Cpolygone*)prim)->Point_ds_polygone(tab_point[0],pt_tol,strict);
}

//-------------------------------------------
// Fonction testant si le segment pt1 pt2
// est inclus dans la polyligne this
//-------------------------------------------	
bool Cprimitive::Segment_inclus(Point3D pt1,Point3D pt2)
{
	Point3D p1,p2;

	p1=tab_point[0];

	for (size_t i=1;i<tab_point.size();i++)
	{
		p2=tab_point[i];
		if((p1==pt1 && p2==pt2) || (p1==pt2 && p2==pt1)) return true;
		p1=p2;
	}

	return false;
}


//----------------------------------------
// operateur == : test sur les tableaux de points
// uniquement (utile pour les calculs intersections)
//----------------------------------------
bool operator==(Cprimitive prim1,Cprimitive prim2)
{
	int	n1 = prim1.tab_point.size(),
		n2 = prim2.tab_point.size();

	if (n1 != n2) return false;

	for (int i = 0 ; i < n1 ; i++)
	{
//		if (prim2.Find_point(prim1.tab_point[i]->pt,i,precis) != i)
		if (prim2.Find_point(prim1.tab_point[i],i,precis) == -1)
		{
			return false;
		}
	}
	
	return true;

}

//-------------------------------------------
// Elimination de la primitive issue d'une
// intersection
//-------------------------------------------	
void Cprimitive::Del_intersect()
{
	int n = tab_point.size();
	if (!n) return;
	
	for (int i = n-1 ; i >= 0  ; i--)
	{
		if (tab_point[i].calculated_ == true)
		{
			Del_point(i);
		}
	}
}

//-------------------------------------------
// Intersection
//-------------------------------------------	
bool Cprimitive::intersect(Cprimitive*,bool)
{
	return false;
}

}



