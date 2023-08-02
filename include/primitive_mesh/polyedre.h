#include <primitive_mesh/plan.h>


#if !defined(LesPolyedres )
	#define LesPolyedres

namespace fractures_intersect {

//-------------------------------------------
// Version 1.0 Pascal SIEGEL - 2000
//-------------------------------------------

//-------------------------------------------
//
// classe Cpolyedre : public Cbase_primitive
// primitive volumique limitee par un ensemble 
// d'objets plan en forme de volume ferme
//-------------------------------------------

class Cpolyedre : public Cbase_primitive
{
	// variables
protected:
	Tpolygone	type_poly;

public :

	std::vector< std::shared_ptr<Cplan> > tab_plan; // tableau de primitive plan

	// constructeur, destructeur
	Cpolyedre(std::string sz = "poly",Tprimitive nv_type = Polyedre);
	virtual ~Cpolyedre();

	// Attributs
	Tpolygone		Get_type_poly() const;
	
	void			Set_type_poly(Tpolygone nv_type);
	void            Set_type_poly(bool typ);
	
	// Operations
	bool			Point_ds_polyedre(Point3D pt,bool strict = true);
	void			Fusionne_primitives(),
					Del_plans(int i = -1),
					scinde_lignes();
	
	virtual void	intersect(Cbase_primitive *prim);
	virtual void	Del_intersect();
	virtual void	Actualise_plans();
	virtual void	Set_Bounding_box();
	virtual void	Set_nom(std::string nv_nom);
	virtual void    Set_ind_materiaux(int nv_ind);


	void Select_surface();
	bool check_possible_intersect_polyedre(Cbase_primitive *prim);
	bool abutting_fracture(Cbase_primitive *prim);
	
};
void Del_polyedre(std::vector<std::shared_ptr<Cpolyedre> >& tab,int i);
std::ostream& operator<<(std::ostream& os, Cpolyedre* poly);

}

#endif
