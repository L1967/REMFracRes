#include <primitive_mesh/boite.h>
#include <primitive_mesh/cylindre.h>
#include <primitive_mesh/faille.h>
#include <primitive_mesh/fracture_set_geo.h>

#include <numeric/property.h>




#if !defined( LesContrainte_Geo)
	#define LesContrainte_Geo

namespace fractures_intersect {

//-------------------------------------------------------
// Version 1.0 modifiee - Pascal SIEGEL - 2000
//-------------------------------------------------------

//-------------------------------------------------------------
// classe Ccontrainte_Geo
// regroupe les contraintes geometriques et geologiques
// discretis�es ou non 
//
// compos� d'un tableau de polyedres
// + tableau de  plan 
// + tableau de primitives
//
// Version 1.0 Pascal SIEGEL - 2000
//-------------------------------------------------------------

class Ccontrainte_Geo
{
public:

	// variables

	std::vector<std::shared_ptr<Cpolyedre> > tab_polyedres;	//  tableau des polyedres
	std::vector<std::shared_ptr<Cplan> > tab_plans;				//  tableau des plans
	std::vector<std::shared_ptr<Cprimitive> > tab_primitives; // tableau des primitives 3D non inclue ds les polyedres ou ds les plans
	std::vector<std::shared_ptr<fracture_set_geo> > tab_fractures_set_; // tableau des primitives 3D non inclue


protected:

	TBoite domaine;							// domaine d'etude;


public:
	
	// constructeur, destructeur
    Ccontrainte_Geo();
    ~Ccontrainte_Geo(); 

	// Attributs
	void	Set_domaine(TBoite dom);
	TBoite	Get_domaine() const;

		// Operations
	std::shared_ptr<Cbase_primitive> Recherche_par_nom(std::string nom);
	std::shared_ptr<Cpolyedre>		 Get_countour();
	void				intersect_all(),
						intersect_all_fracture_set(int number_of_thread),
						Fusionne_primitives(),
						Del_primitives(int i = -1),
						Del_plans(int i = -1),
						Del_polyedres(int i = -1),
						scinde_lignes(),
						Del_intersect(),
						Remove_by_name(std::string nom),
						Add_primitive(std::shared_ptr<Cbase_primitive> prim),
						Change_contour(std::shared_ptr<Cpolyedre> poly),
						Set_countour(std::shared_ptr<Cpolyedre> poly),
						Set_primitive(Cbase_primitive* prim);

	void Remove_all();
	void add_fracture_set(std::shared_ptr<fracture_set_geo> fracture_set);


	void get_fracture_set_surface_mesh(int fracture_set_index, std::vector<Point3D>& tab_point, std::vector<std::vector<int>>& tab_triangle, std::vector<property<float>>&  tab_triangle_fracture_property) ;

	// fonctions de lecture ecriture 
	friend std::ostream& operator<<(std::ostream& os, Ccontrainte_Geo* geo_imp);
	friend std::istream& operator>>(std::istream& is, Ccontrainte_Geo* geo_imp); 
	
	/**
	 *	read the data discretisation + boundary  
	 */

	bool read_data(std::string filename);

	/**
	 *	read the fractures data   
	 */

	bool  read_fractures(std::string filename);




};

}


#endif
