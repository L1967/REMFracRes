//---------------------------
// Version 1.0 - Pascal SIEGEL - 2000
//-------------------------------------------

#include <primitive_mesh/boite.h>

namespace fractures_intersect {

extern double Orient3D(Point3D,Point3D,Point3D,Point3D);
extern double TOL;

const long MAXCAR = 120;
//-----------------------------------------------------
// Constructeur de la structure TBoite
//-----------------------------------------------------
TBoite::TBoite(Point3D pt,double a,double b,double c):
pt0(pt),dx(a),dy(b),dz(c) {}

//----------------------------------------
// operateur == avec un double
//----------------------------------------
bool TBoite::operator==(Point3D a)
{
	Point3D p(dx,dy,dz);
	return (p == a);
}

//-----------------------------------------------------
// operateur d'ecriture pour la structure TBoite
//-----------------------------------------------------
std::ostream& operator<<(std::ostream& os, TBoite& bo)
{   
	os << std::string("Point reference/Origine : ") << bo.pt0 << "\n";
	os << std::string("dx, dy, dz : ") << bo.dx << std::string(" ");
	os << bo.dy << std::string(" ") << bo.dz << "\n";
	return os;
} 

//-----------------------------------------------------
// operateur de lecture pour la structure TBoite
//-----------------------------------------------------
std::istream& operator>>(std::istream& is, TBoite& bo)
{
	char d[MAXCAR];
	
	is.getline(d,MAXCAR,':');
	is >> bo.pt0;
	is.getline(d,MAXCAR,'\n');
	is.getline(d,MAXCAR,':');
	is >> bo.dx >> bo.dy >> bo.dz;
	is.getline(d,MAXCAR,'\n');
	return is;
} 



//-------------------------------------------
// classe Cboite_rect : public Cpolyedre
//
// constructeur/destructeur de la classe Cboite_rect
//-------------------------------------------

Cboite_rect::Cboite_rect(TBoite bo):
Cpolyedre("boite_rect",Boite_Rect),box(bo) {}

Cboite_rect::~Cboite_rect() {}

//-----------------------------------------------------
// Attributs
//-----------------------------------------------------
	
Point3D Cboite_rect::Get_point_ref() const{ return box.pt0; }
void Cboite_rect::Set_point_ref(Point3D pt) { box.pt0 = pt; }

TBoite Cboite_rect::Get_box() const{ return box; }
void Cboite_rect::Set_box(TBoite bo) { box = bo; }

double Cboite_rect::Get_dx() const{ return box.dx; }
void Cboite_rect::Set_dx(double val) { box.dx = val; }

double Cboite_rect::Get_dy() const{ return box.dy; }
void Cboite_rect::Set_dy(double val) { box.dy = val; }

double Cboite_rect::Get_dz() const{ return box.dz; }
void Cboite_rect::Set_dz(double val) { box.dz = val; }

Point3D Cboite_rect::Get_point_coin(int i) const
{ 
	Point3D pti = box.pt0;
	switch (i)
	{
		default: break;

		case 7 : pti.y() += box.dy;
		case 4 : pti.z() += box.dz;
		case 1 : pti.x() += box.dx; break;

		case 6 : pti.x() += box.dx;
		case 2 : pti.y() += box.dy; break;
		
		case 5 : pti.y() += box.dy;
		case 3 : pti.z() += box.dz; break;
	}
	return pti;
}

//-----------------------------------------------------
// Operations
//-----------------------------------------------------
void Cboite_rect::Actualise_plans()
{
	int nplan = tab_plan.size();
	if (nplan) Del_plans();
	
	tab_plan.resize(6);
	Point3D	pt1(0.,box.dy,0.0),pt;
	Vector3D	vn(0,-1,0); // Normale exterieure au 1er plan

	std::string		sz = nom + "_plan0";

	tab_plan[0] = std::shared_ptr<Cplan>(new Cplan(sz,vn,box.pt0));
	tab_plan[0]->Set_taille_element(taille_element);
	tab_plan[0]->Set_ind_materiaux(ind_materiaux);
	tab_plan[0]->tab_primitives.resize(1);
	tab_plan[0]->tab_primitives[0] = std::shared_ptr<Cpolygone> (new Cpolygone("plan0_contour"));
	
	pt.vector() = box.pt0.vector() + pt1.vector();
	vn = -vn;
	sz = nom + "_plan1";
	tab_plan[1] = std::shared_ptr<Cplan>(new Cplan(sz,vn,pt));
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
	pt = box.pt0;
	prim0->tab_point.push_back(pt);
	pt.x() += box.dx;
	prim0->tab_point.push_back(pt);
	pt.z() += box.dz;
	prim0->tab_point.push_back(pt);
	pt.x() = box.pt0.x();
	prim0->tab_point.push_back(pt);
	pt = box.pt0;
	prim0->tab_point.push_back(pt);
	pt.vector() = pt.vector() + pt1.vector();
	prim1->tab_point.push_back(pt);
	pt.z() += box.dz;
	prim1->tab_point.push_back(pt);
	pt.x() += box.dx;
	prim1->tab_point.push_back(pt);
	pt.z() = box.pt0.z();
	prim1->tab_point.push_back(pt);
	pt.x() = box.pt0.x();
	prim1->tab_point.push_back(pt);

	Cpolyedre::Actualise_plans();
}

//--------------------------------------------------------------------
// classe Cboite : public Cpolyedre
//
// constructeur/destructeur de la classe Cboite
//-------------------------------------------

Cboite::Cboite(): Cpolyedre("couche",Boite)
{
	epaisseur = 0.;
	bfaille = true;
}

Cboite::Cboite(Point3D pt_ref,double azim,double pend,
double ext_selon_lpgp,double ext_ortho_lpgp,double ep):
Cpolyedre("couche",Boite)
{
	plan_ref = Cfaille(0,pt_ref,azim,pend,ext_selon_lpgp,ext_ortho_lpgp, ep, ep*10);
	epaisseur = ep;
	bfaille = true;
	Actualise_boite();
}

Cboite::Cboite(Point3D pt[4][2],TBoite *box):
Cpolyedre("boite",Boite)
{
	// point_def[][j] : j = 0 => points du plan top
	//					j = 1 => points du plan bottom
	// on verifiera ulterieurement avec la fonction
	// Verifie_points_def();

	for (int i = 0 ; i < 4 ; i++)
	{
		for (int j = 0 ; j < 2 ; j++)
		{
			point_def[i][j] = pt[i][j];
		}
	}
	bfaille = false;
	Actualise_boite(box);
}

Cboite::~Cboite() {}

//-----------------------------------------------------
// Attributs
//-----------------------------------------------------
	
Point3D Cboite::Get_point_def(int i,int j) const
{ 
	if (i < 4 && j < 2) return point_def[i][j];
	return Point3D();
}

void Cboite::Set_point_def(int i,int j,Point3D pt)
{
	if (i < 4 && j < 2) point_def[i][j] = pt; 
}

double Cboite::Get_epaisseur() const{ return epaisseur; }
void Cboite::Set_epaisseur(double val) { epaisseur = val; }

Cfaille *Cboite::Get_plan_ref() { return &plan_ref; }

void Cboite::Set_plan_ref(const Cfaille& f) 
{ plan_ref = f; }

//-----------------------------------------------------
// Operations
//-----------------------------------------------------
bool Cboite::Verifie_points_def()
{
	if (bfaille == true) return false;

	bool stop = false;
	
	for (int i = 0 ; i < 4 && !stop ; i++)
	{
		if (point_def[i][0] < point_def[i][1])
			stop = true;
	}
	if (stop == true) return false;
	
	for (int j = 0 ; j < 2 ; j++)
	{
		if (Orient3D(point_def[0][j],point_def[1][j],
					point_def[2][j],point_def[3][j]) != 0)
		{
			// les 4 points donnes pour top (j = 0) ou 
			// bottom (j = 1) ne sont pas coplanaires
			// => on recalcule le z du 4e point
			// pour qu'ils le soient
			Cpolygone poly;
			for (int i = 0 ; i < 3 ; i++)
				poly.tab_point.push_back(point_def[i][j]);
			poly.tab_point.push_back(point_def[0][j]);
			
			std::shared_ptr<Cplan>	pl =  std::shared_ptr<Cplan>(poly.get_plane());
			double	c = pl->Get_coef(2);
			
			if (pl == NULL || fabs(c) < TOL) return false;
			
			double	a = pl->Get_coef(0),
					b = pl->Get_coef(1),
					d = pl->Get_coef(3);

			point_def[3][j].z() = -(a*point_def[3][j].x() + b*point_def[3][j].y() + d)/c;
		}
			
		// Les 4 points sont maintenant coplanaires
		// On verifie qu'ils forment un polygone valide
		// (pas de croisements) sinon, on change l'ordre
		// de definition dans le tableau (on ne le fait que pour
		// j = 0, car les points top/bottom marchent par paires)
		Cpolyligne	p1(point_def[0][0],point_def[3][0]),
					p2(point_def[1][0],point_def[2][0]),
					p3(point_def[0][0],point_def[1][0]),
					p4(point_def[2][0],point_def[3][0]);
		Point3D	pt;

		if (intersect_lignes(&p1,&p2,false) == true)
		{
			// Croisement entre segments [0,3] et [1,2]
			// On inverse 2 et 3
			for (j = 0 ; j < 2 ; j++)
			{
				pt = point_def[3][j];
				point_def[3][j] = point_def[2][j];
				point_def[2][j] = pt;
			}
		}
		else if (intersect_lignes(&p3,&p4,false) == true)
		{
			// Croisement entre segments [0,1] et [2,3]
			// On inverse 1 et 2
			for (j = 0 ; j < 2 ; j++)
			{
				pt = point_def[2][j];
				point_def[2][j] = point_def[1][j];
				point_def[1][j] = pt;
			}
		}
	}
	return true;
}

bool Cboite::Actualise_boite(TBoite *box)
{
	// box = boite englobante (domaine d'etude)
	if (tab_plan.size()) Del_plans();
	tab_plan.resize(6);
	
	std::string		sz;

	if (bfaille == false)
	{
		// Creation de la boite par 8 points
		// On verifie que les points sont valides a ce niveau
		if (!Verifie_points_def()) return false;

		Cpolygone	poly[2],pl0,pl1;
		Vector3D	vz(0.,0.,1.);
		
		sz = nom + "_top";		
		for (int i = 0 ; i < 4 ; i++)
		{
			// Plan top
			pl0.tab_point.push_back(point_def[i][0]);
		}
		Vector3D vn = pl0.Get_vec_normal();
		
		if ( (vn).dot(vz) <= 0.)
		{
			for (int i = 3 ; i >= 0 ; i--)
			{
				poly[0].tab_point.push_back(pl0.tab_point[i]);
			}
		}
		else
		{
			for (int i = 0 ; i < 4 ; i++)
			{
				poly[0].tab_point.push_back(pl0.tab_point[i]);
			}
			vn = poly[0].Get_vec_normal();
		}
		tab_plan[0] = std::shared_ptr<Cplan>(new Cplan(sz,vn,point_def[0][0]));
		
		for (int i = 0 ; i < 4 ; i++)
		{
			// Plan bottom
			pl1.tab_point.push_back(point_def[i][1]);
		}
		vn = pl1.Get_vec_normal();
		
		if ((vn).dot(vz) >= 0.)
		{
			for (int i = 3 ; i >= 0 ; i--)
			{
				poly[1].tab_point.push_back(pl1.tab_point[i]);
			}
			vn = poly[1].Get_vec_normal();
		}
		else
		{
			for (int i = 0 ; i < 4 ; i++)
			{
				poly[1].tab_point.push_back(pl1.tab_point[i]);
			}
		}
		sz = nom + "_bottom";
		tab_plan[1] = std::shared_ptr<Cplan>(new Cplan(sz,vn,point_def[0][1]));
		
		for (int j = 0 ; j < 2 ; j++)
		{
			poly[j].tab_point.push_back(point_def[0][j]);
			poly[j].Set_taille_element(taille_element);
			poly[j].Set_ind_materiaux(ind_materiaux);
			tab_plan[j]->Set_taille_element(taille_element);
			tab_plan[j]->Set_ind_materiaux(ind_materiaux);
			tab_plan[j]->tab_primitives.resize(1);
		}
		tab_plan[0]->tab_primitives[0] = std::shared_ptr<Cpolygone>(new Cpolygone(poly[0]));
		tab_plan[0]->tab_primitives[0]->Set_nom("top_contour");

		tab_plan[1]->tab_primitives[0] = std::shared_ptr<Cpolygone>(new Cpolygone(poly[1]));
		tab_plan[1]->tab_primitives[0]->Set_nom("bottom_contour");

		
		Cpolyedre::Actualise_plans();
	}
	else
	{


		Point3D	    ptExtrem[2];// points du segment de directrice
		ptExtrem[0].coord_ = plan_ref.Get_point_ref().coord_ - 0.5*epaisseur * plan_ref.Get_vec_normal();
		ptExtrem[1].coord_ = plan_ref.Get_point_ref().coord_ + 0.5*epaisseur * plan_ref.Get_vec_normal();

		int nplan = tab_plan.size();

		if (nplan) tab_plan.clear();
		tab_plan.resize(nplan);


		Vector3D	vn = -plan_ref.Get_vec_normal();
		std::string sz = nom + "_plan0";

		tab_plan[0] = std::shared_ptr<Cplan>(new Cplan(sz,vn,ptExtrem[0]));
		tab_plan[0]->update_RS(plan_ref.Get_R());
		tab_plan[0]->tab_primitives.resize(1);
		tab_plan[0]->tab_primitives[0] = std::shared_ptr<Cpolygone>(new Cpolygone("plan0_contour"));

		sz = nom + "_plan1";
		tab_plan[1] = std::shared_ptr<Cplan>(new Cplan(sz,plan_ref.Get_vec_normal(),ptExtrem[1]));
		tab_plan[1]->update_RS(plan_ref.Get_R());
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

		for(int i = 0 ; i < plan_ref.tab_primitives[0]->tab_point.size(); i++)
		{
			Point3D	pt;
			pt.coord_= plan_ref.tab_primitives[0]->tab_point[i].coord_ - 0.5*epaisseur * plan_ref.Get_vec_normal();
			prim0->tab_point.push_back(pt);

		}

		for(int i =  plan_ref.tab_primitives[0]->tab_point.size()-1 ; i > -1; i--)
		{
			Point3D	pt_top;
			pt_top.coord_ = plan_ref.tab_primitives[0]->tab_point[i].coord_ + 0.5*epaisseur * plan_ref.Get_vec_normal();

			prim1->tab_point.push_back(pt_top);
		}

		Cpolyedre::Actualise_plans();
		Cpolyedre::Set_Bounding_box();

	}


	return true;
}


void Cboite::Actualise_boite_damage_zone(float damage_foot_wall, float damage_hanging_wall)
{
	// box = boite englobante (domaine d'etude)
	if (tab_plan.size()) Del_plans();
	tab_plan.resize(6);



		Point3D	    ptExtrem[2];// points du segment de directrice
		ptExtrem[0].coord_ = plan_ref.Get_point_ref().coord_ - damage_hanging_wall*epaisseur * plan_ref.Get_vec_normal();
		ptExtrem[1].coord_ = plan_ref.Get_point_ref().coord_ + damage_foot_wall*epaisseur * plan_ref.Get_vec_normal();

		int nplan = tab_plan.size();

		if (nplan) tab_plan.clear();
		tab_plan.resize(nplan);


		Vector3D	vn = -plan_ref.Get_vec_normal();
		std::string sz = nom + "_plan0";

		tab_plan[0] = std::shared_ptr<Cplan>(new Cplan(sz,vn,ptExtrem[0]));
		tab_plan[0]->update_RS(plan_ref.Get_R());
		tab_plan[0]->tab_primitives.resize(1);
		tab_plan[0]->tab_primitives[0] = std::shared_ptr<Cpolygone>(new Cpolygone("plan0_contour"));

		sz = nom + "_plan1";
		tab_plan[1] = std::shared_ptr<Cplan>(new Cplan(sz,plan_ref.Get_vec_normal(),ptExtrem[1]));
		tab_plan[1]->update_RS(plan_ref.Get_R());
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

		for(int i = 0 ; i < plan_ref.tab_primitives[0]->tab_point.size(); i++)
		{
			Point3D	pt;
			pt.coord_= plan_ref.tab_primitives[0]->tab_point[i].coord_ - damage_hanging_wall*epaisseur * plan_ref.Get_vec_normal();
			prim0->tab_point.push_back(pt);

		}

		for(int i =  plan_ref.tab_primitives[0]->tab_point.size()-1 ; i > -1; i--)
		{
			Point3D	pt_top;
			pt_top.coord_ = plan_ref.tab_primitives[0]->tab_point[i].coord_ + damage_foot_wall*epaisseur * plan_ref.Get_vec_normal();

			prim1->tab_point.push_back(pt_top);
		}

		Cpolyedre::Actualise_plans();
		Cpolyedre::Set_Bounding_box();


}


}
