#include <primitive_mesh/primitive.h>
#include <vector>

#if !defined(LesPrimitives_derive)
	#define LesPrimitives_derive

namespace fractures_intersect {

//-------------------------------------------
// Version 1.0 Pascal SIEGEL - 2000
//-------------------------------------------

//-------------------------------------------
// classe Cpolyligne : public Cprimitive 
// tableau de Point3D 
//-------------------------------------------

class Cpolyligne : public Cprimitive
{
	// variables
public:
	// indicateurs servant pour les lignes discontinues
	std::vector<bool> tab_connec;
						

	// constructeur, destructeur
	Cpolyligne(std::string sz = "ligne",
		Tprimitive nv_type = Polyligne);
	Cpolyligne(Point3D a, Point3D b);
	Cpolyligne(const Cpolyligne& poly);
	virtual ~Cpolyligne();
	
	// Operations
	Cpolyligne&	operator=(const Cpolyligne& poly);

	bool			Point_ds_polyligne(Point3D pt,double tol,
										bool cal_inter = false),
					scinde_ligne(Cpolyligne *p),
					connecte_ligne(Cpolyligne *p),
					Segment_inclus(Cpolyligne *p,int i,int j),
					Insert(Point3D p0,Point3D p1, Point3D p,bool connec = true);
	void			Del_connec(int i = -1);
	void			Test_segment_inclus();

	virtual bool	intersect(Cprimitive *prim,bool benrichi = true);
    virtual bool	Primitive_incluse(Cprimitive *prim,bool strict = true,bool cal_inter = true);
};

//-------------------------------------------
// classe Cpolygone : public Cpolyligne
// polyligne fermee + taille des elements
// fonction de localisation ds un plan
//-------------------------------------------
class Cplan;
class Cpolygone : public Cpolyligne
{

	// variables
protected:
	Tpolygone	type_poly;		// type de polygone



public:

	// constructeur, destructeur
	Cpolygone(std::string sz = "poly");
	Cpolygone(const Cpolygone& poly);
	Cpolygone(const Cpolyligne& poly);
	~Cpolygone();

	// Attributs 
	Tpolygone	Get_type_poly() const;
	void		Set_type_poly(Tpolygone nv_type);

	// Operations
	Vector3D		Get_vec_normal();
	Cplan* get_plane();

	bool			Point_ds_polygone(Point3D pt,double tol,bool strict = true,
									bool cal_inter = false);
	bool			Point_ds_polygone(Cplan* plane, Point3D pt,double tol,bool strict = true,
									bool cal_inter = false);

	void			update_point2D(Cplan* plane);


	virtual bool	intersect(Cprimitive *prim,bool benrichi = true);


};
// Operations d'intersection
bool intersect_lignes(Cpolyligne *lin1,Cpolyligne *lin2,bool blin);
bool intersect_polygone_ligne(Cpolygone *poly,Cpolyligne *lin,bool bpoly);
bool intersect_polygones(Cpolygone *poly1,Cpolygone *poly2);

}

#endif
