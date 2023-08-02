#include <primitive_mesh/boite.h>
#include <primitive_mesh/cylindre.h>
#include <primitive_mesh/faille.h>
#include <primitive_mesh/FracResBeddingLocalisation.h>
#include <numeric/property.h>


#if !defined(Fracture_set_geo)
#define Fracture_set_geo

namespace fractures_intersect {

enum FracResIntersectionCategoryEnum {
  BED_INTERFACE = 0,
  MECHANIAL_UNITS = 1,
  BED_PARALLEL = 2,
  JOINTS = 3,
  FAULTS = 4,
  FAULT_ARRAYS =5,
  FRACTURE_CLUSTERS_JOINT_SWARM = 6,
  FRACTURE_CLUSTERS = 7,
  FRACTURE_CLUSTERS_CORRIDORS_ZONES = 8,
  FRACTURE_CLUSTERS_FAULT_ZONES = 9



};

enum FracResIntersectionTypeEnum {
  ABUTTING = 0,
  CROSSING = 1,
  RANDOM = 2,
  ABBUTING_CROSSING=3
};

//-------------------------------------------------------------
// class fracture_set_geo
//-------------------------------------------------------------

class fracture_set_geo
{
public:

	// variables

	std::vector<std::shared_ptr<Cfaille> > tab_fractures_;
	std::vector<std::shared_ptr<Ccylindre> > tab_stress_zone_cylinders_;

	int index_;
	int seed_;
	FracResIntersectionCategoryEnum set_type_;

	std::vector< FracResIntersectionTypeEnum > tab_type_intersection_;
	std::vector< std::shared_ptr<FracResBeddingLocalisation> > tab_bedding_localisation_;



public:
	
	// constructeur, destructeur
    fracture_set_geo(int index, FracResIntersectionCategoryEnum set_type, int seed, std::vector<FracResIntersectionTypeEnum>& tab_type_intersection ,
    		 std::vector<std::shared_ptr<Cfaille> >& tab_fractures);
    fracture_set_geo(int index, FracResIntersectionCategoryEnum set_type, int seed, std::vector<FracResIntersectionTypeEnum>& tab_type_intersection);
    fracture_set_geo();
    ~fracture_set_geo();

	void add_fracture(std::shared_ptr<Cfaille> fracture);
	void add_fracture_and_abbuting(std::shared_ptr<Cfaille> fracture_add);
	void add_fracture_and_create_cylinder(std::shared_ptr<Cfaille> fracture, int n_plan = 12);
	void add_bedding_localisation(std::shared_ptr<FracResBeddingLocalisation>  localisation);
	void clear_stress_zone_cylinder();

	bool in_all_stress_zone(Point3D& pt);
	bool in_all_stress_zone_3d(std::shared_ptr<Cfaille> faille);
	bool in_all_stress_zone__abutting_3d(std::shared_ptr<Cfaille> faille,int n_plan);
	void get_fracture_set_surface_mesh(std::vector<Point3D>& tab_point, std::vector<std::vector<int>>& tab_triangle, std::vector<property<float>>&  tab_triangle_fracture_property) ;
	void get_fracture_set_protection_zone_surface_mesh(std::vector<Point3D>& tab_point, std::vector<std::vector<int>>& tab_triangle ,  std::vector<property<float>>& tab_triangle_poly_index);
	void abbuting_all_fractures();

};

}


#endif
