#include <primitive_mesh/bounding_box.h>
#include<primitive_mesh/point.h>
#include<Eigen/Core>
#include<string>
#include<vector>
#include<memory>

#if !defined(LesPrimitives )
	#define LesPrimitives

namespace fractures_intersect {




//-------------------------------------------
// Version 1.0 Pascal SIEGEL - 2000
//-------------------------------------------

// type d'objets primitives 2D/3D
enum Tprimitive {Point,Polyligne,Polygone,Plan,Faille,
				Polyedre,Cylindre,Boite_Rect,Boite};

// type d'objets polygones/polyedres
enum Tpolygone {Trou,Zone_affinement};


//-------------------------------------------
// classe Cbase_primitive
// classe g�n�rique de base utilis�e pour definir les
// primitives, les plans, polyedres, etc
//
//-------------------------------------------
class Cbase_primitive
{
	// Variables membres
protected:
    
	bounding_box		boundingBox_;// points min et max de la bounding box
	Tprimitive		type;		// type de primitive de base
	std::string		nom;		// nom de la primitive de base
	double			taille_element;// taille des elements
	int				ind_materiaux; // indice du materiau
	bool			calcule,	// si true primitive issue d'intersections
					Selection,  // Selection o/n de la primitive
					active;	// indique si la primitive est active ou non
							// si oui, on la prend en compte dans le calcul
							// des intersections ou le maillage


public:
	
	// Constructeur/destructeur
	Cbase_primitive(std::string sz = "point",Tprimitive nv_type = Point,
					bool cal = false);
	Cbase_primitive(const Cbase_primitive& prim);
	virtual ~Cbase_primitive();

	void operator = (Cbase_primitive prim);
	
	// Attributs
	Tprimitive		Get_type() const;
	std::string		Get_nom() const;
	double			Get_taille_element() const;
	int				Get_ind_materiaux() const;
	bool			Get_calcule() const;
	bool            Bounding_box_existe();
	bool			Selectionne_Moi() const;
    void			Set_type(Tprimitive nv_type),
					Set_taille_element(double taille),
					Set_calcule(bool cal),
					Selectionne_Moi(bool sel);
	virtual void    Set_ind_materiaux(int nv_ind);
	virtual void    Set_nom(std::string nv_nom);

	bounding_box&    get_boundingBox();
	             
	// Operations
	virtual void	Set_Bounding_box();
	
	bool			In_Bounding_box(Point3D& pt),
					In_Bounding_box(Cbase_primitive *p),
					Intersect_Bounding_box(Cbase_primitive *p);



	double          Distance_max();
					

	void			Actualise_repere(Vector3D& e0,Vector3D& e1,
									Vector3D& e2);

	
};

//-------------------------------------------
// classe Cprimitive : public Cbase_primitive
// forme de base utilis� pour definir des contours
// polygones, polylignes, points....
// par defaut, la classe Cprimitive represente
// un point
//-------------------------------------------
class Cprimitive : public Cbase_primitive
{
	// variables

public : 

	std::vector<Point3D> tab_point; // tableau de points 3D
	std::vector<Point3D> tab_point_2d_; // tableau de points 3D
	
	// constructeur, destructeur

	Cprimitive(std::string sz = "point",Tprimitive nv_type = Point,
		bool cal = false);
	Cprimitive(Point3D pt);
	Cprimitive(const Cprimitive& prim);
	~Cprimitive();

	// Operations
	virtual void	Set_Bounding_box();
	virtual void	Fusionne_points(double tol);
	virtual void	Del_intersect();
	virtual bool	intersect(Cprimitive *prim,bool benrichi = true);
	virtual bool	Primitive_incluse(Cprimitive *prim,bool strict = true,bool cal_inter = true);
	virtual int		Find_point(Point3D pt,int i0,double tol);
	bool            Segment_inclus(Point3D pt1,Point3D pt2);
	
	Cprimitive& operator=(const Cprimitive& prim);
	
	void			Del_point(int i = -1);

	friend bool operator==(Cprimitive prim1,Cprimitive prim2);

	
};
}

#endif
