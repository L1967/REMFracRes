//-------------------------------------------
// Version 1.0 - Pascal SIEGEL - 2000
//-------------------------------------------

#if !defined(__BOITE_H__)
	#define __BOITE_H__

#include<primitive_mesh/point.h>
#include <primitive_mesh/polyedre.h>
#include <primitive_mesh/faille.h>
#include <primitive_mesh/oriented_bounding_box.h>


namespace fractures_intersect {

//-----------------------------------------------------
// structure TBoite : coordonn�es du point d'origine
// + dx,dy,dz
//------------------------------------------------------
struct TBoite
{
	Point3D pt0;
	double dx,dy,dz;
	TBoite(Point3D pt = Point3D(),double a = 0.,
					double b = 0.,double c = 0.);
	
	bool	operator==(Point3D a);
};


//--------------------------------------------------------------------
// classe Cboite_rect : public Cpolyedre
// 
// boite rectangulaire heritant du polyedre 
// (ensemble de plans, utile pour discretiser ou 
// pour les calculs d'intersection)
//--------------------------------------------------------------------

class Cboite_rect : public Cpolyedre
{
	// Variables
protected:
	TBoite		box;
				
public:
	// Construction/destruction
    Cboite_rect(TBoite bo = TBoite());
	~Cboite_rect();
	
	// Attributs
	TBoite		Get_box() const;
	Point3D	Get_point_ref() const;
	Point3D	Get_point_coin(int i) const;
	double		Get_dx() const;
	double		Get_dy() const;
	double		Get_dz() const;

	void		Set_box(TBoite bo),
				Set_point_ref(Point3D pt),
				Set_dx(double val),
				Set_dy(double val),
				Set_dz(double val);
	
	// Operations
	virtual void	Actualise_plans();
};

std::ostream& operator<<(std::ostream& os, TBoite& bo);
std::istream& operator>>(std::istream& is, TBoite& bo);

//--------------------------------------------------------------------
// classe Cboite : public Cpolyedre
// 
// boite quelconque heritant du polyedre 
// volume defini par 8 points ou par azimut,pendage et extensions
// utile pour creer les facies ou couches geologiques
//--------------------------------------------------------------------

class Cboite : public Cpolyedre
{
	// Variables
protected:
	Cfaille		plan_ref;// cas ou la couche est definie avec azimut, 
						 // pendage, extensions et epaisseur
	double		epaisseur;// epaisseur de la couche
	Point3D	point_def[4][2];// points limites definissant le facies

private:
	bool		bfaille;	// mode de creation de la boite
				
public:
	// Construction/destruction
    Cboite();
    Cboite(Point3D pt_ref,double azim,double pend,
		double ext_selon_lpgp,double ext_ortho_lpgp,
		double ep = 0.);
    Cboite(Point3D pt[4][2],TBoite *box = NULL);
	~Cboite();
	
	// Attributs
	Point3D	Get_point_def(int i,int j = 0) const;
	double		Get_epaisseur() const;
	Cfaille		*Get_plan_ref();

	void		Set_plan_ref(const Cfaille& f),
				Set_point_def(int i,int j,Point3D pt),
				Set_epaisseur(double val);
	void        Actualise_boite_damage_zone(float damage_foot_wall, float damage_hanging_wall);
	// Operations
	bool	Actualise_boite(TBoite *box = NULL),
			Verifie_points_def();
};

}

#endif //__BOITE_H__
