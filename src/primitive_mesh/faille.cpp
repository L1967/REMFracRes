
#include <primitive_mesh/boite.h>
#include <primitive_mesh/faille.h>
#include <boost/random.hpp>
#include <boost/format.hpp>
#include <cmath>
#include <fstream>
#include <iostream>

namespace fractures_intersect {

extern double dround(double,double);
extern double rad,pi,pt_tol,TOL;
////////////////////////////////////////////////////////////////
//-------------------------------------------
//
// fonctions de la classe Cfaille : public Cplan
// 
// equation du plan + limites externes
// et internes : contour + tableau de 
// primitives du plan 
// 
// + pour la faille : azimut et pendage
//
// Version 1.0 Pascal SIEGEL - 2000
//-------------------------------------------

//-------------------------------------------
// constructeur/destructeur de la classe Cfaille
//-------------------------------------------

Cfaille::Cfaille(int fracture_set_index , Point3D pt_ref,double azim,double pend,
double ext_selon_lpgp,double ext_ortho_lpgp, double aperture , double stress_zone_max):
Cplan()
{
	pt0_ = pt_ref;
	fracture_set_index_ = fracture_set_index;
	azimut = azim;
	pendage = pend;
	ext2L_selon_lpgp = ext_selon_lpgp;
	ext2l_orthog_lpgp = ext_ortho_lpgp;
	limite = false;
	type = Faille;
	stress_zone_max_= stress_zone_max;
	aperture_ = aperture;
	connection_relay_ =0;
	ratio_perm_ = 1.0;
	slope_azimuth_ = 0.0;
	familly_index_ =1;

	Creation_faille_finie();


}

Cfaille::Cfaille(const Cfaille& f): Cplan((const Cplan&)f)
{		
	pt0_= f.pt0_;
	fracture_set_index_ = f.fracture_set_index_;
	azimut = f.azimut;
	pendage = f.pendage;
	ext2L_selon_lpgp = f.ext2L_selon_lpgp;
	ext2l_orthog_lpgp = f.ext2l_orthog_lpgp;
	limite = false;
	type = Faille;
	stress_zone_max_= 0.0;
	aperture_ = f.aperture_;
	connection_relay_ =f.connection_relay_;
	ratio_perm_ = f.ratio_perm_;
	slope_azimuth_ = f.slope_azimuth_;
	familly_index_ = f.familly_index_;
	Creation_faille_finie();
}

// constructeur, destructeur
Cfaille::Cfaille(Cpolygone* polygone)
{
	tab_primitives.resize(1);
	tab_primitives[0] = std::shared_ptr<Cpolygone>(polygone);

	fracture_set_index_ = 0;
	aperture_ = 0.0;

	stress_zone_max_= 0.0;
	azimut = 0.0;
	pendage = 0.0;
	ext2L_selon_lpgp = 0.0;
	ext2l_orthog_lpgp = 0.0;

	limite = true;
	connection_relay_ =0;
	ratio_perm_ = 1.0;
	slope_azimuth_ = 0.0;
	familly_index_ =1;
	Set_point_ref(polygone->tab_point[0]);

	Set_vec_normal();
	Set_RS();

	Set_equation_plan();
	std::dynamic_pointer_cast<Cpolygone>(tab_primitives[0])->update_point2D(this);
}


// constructeur, destructeur
Cfaille::Cfaille(int fracture_set_index, double aperture, std::vector<Point3D>& poly, std::string name)
{
	tab_primitives.resize(1);
	tab_primitives[0] = std::shared_ptr<Cpolygone>(new Cpolygone(name));

	if (poly.empty()) return;

	for (size_t i = 0; i < poly.size(); i++) {
		tab_primitives[0]->tab_point.push_back(poly[i]);
	}

	tab_primitives[0]->tab_point.push_back(poly[0]);

	Set_point_ref(tab_primitives[0]->tab_point[0]);

	Set_vec_normal();
	Set_RS();

	Set_equation_plan();

	azimut = 0.0;
	pendage = 0.0;
	ext2L_selon_lpgp = 0.0;
	ext2l_orthog_lpgp = 0.0;

	fracture_set_index_ = fracture_set_index;
	aperture_ = aperture;

	stress_zone_max_= 0.0;
	connection_relay_ =0;

	limite = true;
	ratio_perm_ = 1.0;
	slope_azimuth_ = 0.0;
	familly_index_ = 1;

	std::dynamic_pointer_cast<Cpolygone>(tab_primitives[0])->update_point2D(this);
}



Cfaille::~Cfaille() {}

//-------------------------------------------
// fonctions d'acces au variables
//-------------------------------------------

double Cfaille::Get_azimut() const { return azimut; }
void Cfaille::Set_azimut(double azim) { azimut = azim; }
double Cfaille::Get_SlopeAzimut() const { return slope_azimuth_; }
void Cfaille::Set_SlopeAzimut(double azim) { slope_azimuth_ = azim; }
double Cfaille::Get_RatioPermeability() const { return ratio_perm_; }
void Cfaille::Set_RatioPermeability(double ratio) { ratio_perm_ = ratio; }
double Cfaille::Get_aperture() const { return aperture_; }
void Cfaille::Set_aperture(double ap) { aperture_ = ap; }
double Cfaille::Get_pendage() const { return pendage; }
void Cfaille::Set_pendage(double pend) { pendage = pend; }

bool Cfaille::Get_limite() const { return limite; }
void   Cfaille::Set_limite(bool limit) { limite=limit; }

double Cfaille::Get_ext2L_selon_lpgp() const
{ return ext2L_selon_lpgp; }

void Cfaille::Set_ext2L_selon_lpgp(double l)
{ ext2L_selon_lpgp = l; }

double Cfaille::Get_ext2l_orthog_lpgp() const
{ return ext2l_orthog_lpgp; }

void Cfaille::Set_ext2l_orthog_lpgp(double l)
{ ext2l_orthog_lpgp = l; }

//-------------------------------------------
// fonctions qui remet a jour les vecteur 
// R,S et T qui definissent la faille
//-------------------------------------------
bool Cfaille::Actualise_repere(bool /*cal*/)
{	
	double tetha,alpha;

	tetha = pendage*rad;
	alpha = - azimut*rad;

	if (pendage == 0.0 && azimut == 0.0) // failles horizontales
	{
		R[0] = 1.0;
		R[1] = 0.0;
		R[2] = 0.0;

		S[0] = 0.0;
		S[1] = 1.0;
		S[2] = 0.0;
	}
	else	// autres failles
    	{
	
		R[0] = std::cos(alpha);
		R[1] = std::sin(alpha);
		R[2] = 0.0;
		Point3D ptr(R); 
		ptr.arrondi_selon_precis();
		R = ptr.vector(); 

		S[0] = +std::cos(tetha)*std::sin(alpha);
		S[1] = -std::cos(tetha)*std::cos(alpha);
		S[2] = std::sin(tetha);
		Point3D pts(S); 
		pts.arrondi_selon_precis();
		S = pts.vector(); 
		
	}
	T = R.cross(S);
	Point3D ptt(T);
	ptt.arrondi_selon_precis();
	T = ptt.vector(); 
	return true;
}

//-------------------------------------------
// Creation/actualisation de la faille
//-------------------------------------------
void Cfaille::Creation_faille_finie()
{	
	Actualise_repere();

	if (ext2l_orthog_lpgp == 0. && ext2L_selon_lpgp == 0.) return;
	limite = true;
	
	// on ne modifie que le contour de la faille (primitive 0)
	// si d'autres primitives existent, on n'y touche pas
	Del_primitives();			
	tab_primitives.resize(1);

	tab_primitives[0] = std::shared_ptr<Cpolygone>(new Cpolygone("faille_contour"));
	Point3D pt = Change_coord_2D(pt0_),ptres,pt1;
	


	pt.x() -= ext2L_selon_lpgp/2.;
	pt.y() -= ext2l_orthog_lpgp/2.;

	ptres = Change_coord_3D(pt);
	pt1 = ptres;


	tab_primitives[0]->tab_point.push_back(ptres);
	pt.x() += ext2L_selon_lpgp;;
	ptres = Change_coord_3D(pt);

	tab_primitives[0]->tab_point.push_back(ptres);
	pt.y() += ext2l_orthog_lpgp;
	ptres = Change_coord_3D(pt);

	tab_primitives[0]->tab_point.push_back(ptres);
	pt.x() -= ext2L_selon_lpgp;
	ptres = Change_coord_3D(pt);

	tab_primitives[0]->tab_point.push_back(ptres);
	tab_primitives[0]->tab_point.push_back(pt1);

	Set_equation_plan();

	std::dynamic_pointer_cast<Cpolygone>(tab_primitives[0])->update_point2D(this);

}

bool Cfaille::Set_faille_contour(TBoite *box)
{
	if (limite == true || box == NULL) return false;
	Actualise_repere();
	Set_equation_plan();

	// On cree le polygone contour de la faille (primitive 0)
	// en determinant l'intersection du plan avec le domaine box
	Del_primitives();			

	std::shared_ptr<Cpolyligne>	lin = std::shared_ptr<Cpolyligne>(new Cpolyligne);
	Cboite_rect	boite(*box);
	boite.Actualise_plans();

	for (int i = 0 ; i < 6 ; i++)
	{
		std::shared_ptr<Cplan> plan = boite.tab_plan[i];

		// Vecteur directeur de la droite 3D
		Vector3D vec = T.cross(plan->Get_vec_normal());

		// Les plans sont paralleles
		if (vec.norm() == 0.0)
		{
			// Les plans sont confondus (a pt_tol pres)
			if (plan->Point_ds_plan(pt0_,pt_tol,false) == true)
			{
				// On prend le polygone contour du plan
				// de la boite comme contour de faille
				std::shared_ptr<Cpolygone> p = std::dynamic_pointer_cast<Cpolygone>(plan->tab_primitives[0]);
				std::shared_ptr<Cpolygone> new_polygone = std::shared_ptr<Cpolygone>(new Cpolygone(*p.get()));
				tab_primitives.push_back(new_polygone);
				tab_primitives[0]->Set_nom("faille_contour");

				return true;
			}
		}
		else
		{
			std::shared_ptr<Cpolygone> contour = std::dynamic_pointer_cast<Cpolygone>(plan->tab_primitives[0]);
			Point3D	p1,p0 = contour->tab_point[0],pt;
			Vector3D	v;
			double		u = 0., den = 0.;
			
			for (size_t j = 1 ; j < contour->tab_point.size() ;j++)
			{
				p1 = contour->tab_point[j];
				v = p1.vector() - p0.vector();
				if (!iscolinear(v,vec)) // lignes non parallelles
				{
					den = T.dot(v);
					u = -(T.dot(p0.vector()) + a[3]);
					den = dround(den,TOL);
					
					if (den != 0.)
					{
						u /= den;
						if (u >= 0. && u <= 1.)
						{
							// Le point intersection entre le segment du
							// polygone contour et la droite 3D existe
							pt.vector() = p0.vector() + v*u;
							// on enrichit la ligne
							lin->tab_point.push_back(pt);
						}
					}
				}
				p0 = p1;
			}
		}
	}
		
	if (!lin->tab_point.size())
	{
		return false;
	}

	lin->Fusionne_points(pt_tol);
	lin->tab_point.push_back(lin->tab_point[0]);
	Cpolygone poly(*lin);
	if ((poly.Get_vec_normal()).dot(T) < 0.)
	{
		Point3D pi;
		int n = lin->tab_point.size();
		for (int i = 1 ; i < n/2 ; i++)
		{
			pi = lin->tab_point[n-i-1];
			lin->tab_point[n-i-1] = lin->tab_point[i];
			lin->tab_point[i] = pi;
		}
	}
	tab_primitives.push_back(std::shared_ptr<Cpolygone>(new Cpolygone(*lin)));
	tab_primitives[0]->Set_nom("faille_contour");
	std::dynamic_pointer_cast<Cpolygone>(tab_primitives[0])->update_point2D(this);
	return true;
}

//---------------------------------------------
// Transformation des coordonnees 3D d'un point
// du plan en coordonnes R,S ds le plan
//---------------------------------------------
bool Cfaille::in_Plan(Point3D& pt)
{



	Point3D projPlanPt = Projection_sur(pt);

	if ( (std::dynamic_pointer_cast<Cpolygone>(tab_primitives[0]))->Point_ds_polygone(this,projPlanPt,pt_tol,false) == false) return false;

	return true;


}

//---------------------------------------------
// Transformation des coordonnees 3D d'un point
// du plan en coordonnes R,S ds le plan
//---------------------------------------------
bool Cfaille::in_stress_zone(Point3D& pt)
{

	double dist = distance_to_plan(pt);
	if (dist > 	stress_zone_max_) return false;

	Point3D projPlanPt = Projection_sur(pt);

	if ( (std::dynamic_pointer_cast<Cpolygone>(tab_primitives[0]))->Point_ds_polygone(this,projPlanPt,pt_tol,false) == false) return false;

	Point3D pt2D = Change_coord_2D(pt);
	double x= pt2D.x();

	double a_elipse = ext2L_selon_lpgp/2.0;

	double stress_zone = (stress_zone_max_/a_elipse)*std::sqrt(a_elipse*a_elipse -  x*x);

	if (dist <= stress_zone) return true;
	else return false;


}


void Cfaille::simplify_to_box_poly() {

	std::shared_ptr<Cpolygone>	poly = std::dynamic_pointer_cast<Cpolygone>(tab_primitives[0]);

	double neg_x = 1e20;
	double pos_x = 1e20;
	double neg_y = 1e20;
	double pos_y = 1e20;

	Point3D pt_ref = Change_coord_2D(pt0_);

	int tab_size = tab_primitives[0]->tab_point.size();

	for (int i = 0 ; i < tab_size -1; i++)
	{

		Point3D pt = Change_coord_2D(tab_primitives[0]->tab_point[i]);
		Point3D pt1 = Change_coord_2D(tab_primitives[0]->tab_point[i+1]);


		if (tab_primitives[0]->tab_point[i].calculated_  && tab_primitives[0]->tab_point[i+1].calculated_ && tab_size > 5) {

			bool suppress_one = false;
			bool take_x = false;
			bool first_point = true;

			if (pt.x() * pt1.x() > 0 && std::abs(pt.x()-pt1.x()) < ext2L_selon_lpgp/2) {
				suppress_one = true;
				take_x = true;

				if (std::abs(pt.x()) < std::abs(pt1.x())) {
					first_point = true;
				} else {
					first_point = false;
				}
			}
			if (pt.y() * pt1.y() > 0 && std::abs(pt.y()-pt1.y()) < ext2l_orthog_lpgp/2) {

				if ( (take_x == true  && std::abs(pt.y()-pt1.y()) < std::abs(pt.x()-pt1.x())) || take_x == false)  {

					suppress_one = true;

					if (std::abs(pt.y()) < std::abs(pt1.y())) {
						first_point = true;
					} else {
						first_point = false;
					}
				}
			}

			if (suppress_one == true) {
				if (first_point == false) pt.coord_ = pt1.coord_;
			}
			else {
				pt.coord_ = (pt.coord_ + pt1.coord_)/2.0;
			}

			i++;
		}


		double x = pt.x();
		double y = pt.y();

		if (x < 0) neg_x = std::min(neg_x,-x);
		else if (x > 0) pos_x = std::min(pos_x,x);

		if (y < 0) neg_y = std::min(neg_y,-y);
		else if (y > 0)  pos_y = std::min(pos_y,y);
	}

	if (neg_x > 1e19) neg_x = 0.0;
	if (pos_x > 1e19) pos_x = 0.0;

	if (neg_y > 1e19) neg_y = 0.0;
	if (pos_y > 1e19) pos_y = 0.0;

	Point3D pt_new(pos_x - (pos_x + neg_x)/2, pos_y-(pos_y + neg_y)/2, 0.0);
	pt0_ = Change_coord_3D(pt_new);


	ext2L_selon_lpgp = pos_x + neg_x;
	ext2l_orthog_lpgp = pos_y + neg_y;

	Creation_faille_finie();
}








//-------------------------------------------
// operateur=
//-------------------------------------------
Cfaille& Cfaille::operator=(const Cfaille& f)
{

	pt0_ = f.pt0_;
	azimut = f.azimut;
	pendage = f.pendage;
	ext2L_selon_lpgp = f.ext2L_selon_lpgp;
	ext2l_orthog_lpgp = f.ext2l_orthog_lpgp;
	limite = false;
	type = Faille;
	stress_zone_max_= f.stress_zone_max_;
	aperture_ = f.aperture_;
	ratio_perm_ = f.ratio_perm_;
	slope_azimuth_ = f.slope_azimuth_;



	Creation_faille_finie();
	return *this;
}


void Cfaille::save(std::fstream & fout){
	if(fout.is_open()){
	         fout << pt0_.x() << " " << pt0_.y() << " " << pt0_.z() << " "<< ext2L_selon_lpgp << " "<< ext2l_orthog_lpgp << " "<< aperture_ << " "<< azimut << " "<< pendage << " "<< std::endl;
	}
}


}
