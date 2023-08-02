#include <primitive_mesh/primitive_derive.h>
#include <vector>
#include <memory>

#if !defined(LePlan)
#define LePlan

namespace fractures_intersect {

//-------------------------------------------
// Version 1.0 Pascal SIEGEL - 2000
//-------------------------------------------

class Cmaillage_Plan_trian;
class Cmaillage_Plan;

//-------------------------------------------
//
// classe Cplan : public Cbase_primitive
// 
// equation du plan + limites externes
// et internes : contour + tableau de 
// primitives du plan 
//-------------------------------------------

class Cplan : public Cbase_primitive
{
	// variables	
public : 

	std::vector<std::shared_ptr<Cprimitive>> tab_primitives;	// tableau des primitives
	// primitive 0 = contour
private:
	bool		bequation; // booleen indiquant si l'equation du plan a ete calculee
protected:

	Vector3D	R,S,T;	// vecteurs du repere : T vecteur normal
	Point3D	pt0_;	// un point du plan
	double		a[4];	// coefficients de l'equation du plan 
	// a0x+a1y+a2z+a3=0


public:

	// constructeur, destructeur
	Cplan(std::string sz , Vector3D& vec ,Point3D& pt);
	Cplan();
	Cplan(const Cplan& plan);
	virtual ~Cplan();

	// Attributs
	bool		Get_bequation() const;
	Point3D&	    Get_point_ref();
	double		Get_coef(int i) const;
	Vector3D&	Get_R();
	Vector3D&	Get_S();
	Vector3D&	Get_vec_normal();
	void		Set_coef(double x,int i);
	void		Set_point_ref(Point3D pt);

	// Equation du plan (a,b,c,d) calculee a partir de R,S,T et pt0
	void		Set_equation_plan();

	// Vecteur normal determine en fonction de 3 points du plan
	bool		Set_vec_normal();
	bool		Set_RS();
	bool        update_RS(Vector3D R_fixed);

	// Calcul du repere R,S,T
	virtual bool Actualise_repere(bool cal = false);
	virtual bool Actualise_equation();

	// fonctions liees au maillage
	Cmaillage_Plan_trian* Get_maillage_trian();     // retourne le maillage_plan triangule 
	void Set_maillage_trian(Cmaillage_Plan_trian* m_plan);

	Cmaillage_Plan* Get_maillage();     // retourne le maillage_plan
	void Set_maillage(Cmaillage_Plan* m_plan);

	// Operations
	Point3D	Change_coord_2D(Point3D& pt); // changement de coordonnes 3D->2D plan
	Point3D	Change_coord_3D(Point3D& pt); // changement de coordonnes 2D plan->3D
	Point3D	Get_max_coord2D();  // en coordonn�es planaires
	Point3D	Get_min_coord2D();

	Point3D	Projection_sur(Point3D& pt); // projection du pt sur le plan
	Vector3D Projection_sur(Vector3D& vec); // Projection d'un vecteur sur le plan

	double distance_to_plan(Point3D& pt);
	double oriented_distance_to_plan(Point3D& pt);

	bool Point_ds_polygone(Cprimitive* poly,Point3D& pt);
	bool Point_ds_plan(Point3D& pt,double tol,bool cal = true);
	bool Plan_confondu_avec(Cplan *plan,double tol);

	void intersect(Cbase_primitive *prim,int i0 = 0);
	void intersect_plan(Cplan *plan,double tol);
	bool intersect_plan_limitated(Cplan *plan,double tol);
	bool intersect_plan_limitated_dist(Cplan *plan,double tol, double& dist_min);
	bool check_intersect_plan(Cplan *plan,double tol);

	void Fusionne_primitives();
	void Del_primitives(int i = -1);
	void scinde_lignes();
	void Set_Bounding_box(int iprim);
	void Del_intersect();

	virtual void Set_Bounding_box();

	void Select_noeud_face();

	friend bool	operator==(Cplan p1,Cplan p2);
	friend bool	operator!=(Cplan p1,Cplan p2);

	void add_mesh_data(std::vector<Point3D>& tab_point, std::vector<std::vector<int>>& tab_triangle, std::vector<std::vector<float>>& tab_triangle_fracture_index,
			std::vector<float>& property);


};
void scinde_ligne(std::vector<std::shared_ptr<Cprimitive> >& tab);
void Del_primitive(std::vector<std::shared_ptr<Cprimitive> >& tab,int i);
void Del_plan(std::vector<std::shared_ptr<Cplan> >& tab,int i);
void Fusionne(std::vector<std::shared_ptr<Cprimitive> >& tab);

std::ostream& operator<<(std::ostream& os, Cplan* plan);

}

#endif
