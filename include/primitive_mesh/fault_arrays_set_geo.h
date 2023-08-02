#include <primitive_mesh/boite.h>
#include <primitive_mesh/cylindre.h>
#include <primitive_mesh/fault_array.h>
#include <numeric/property.h>

#if !defined(Fault_arrays_set_geo)
#define Fault_array_set_geo

namespace fractures_intersect {

//-------------------------------------------------------------
// class fault_array_set_geo
//-------------------------------------------------------------

class fault_arrays_set_geo
{
public:

	// variables

	std::vector<std::shared_ptr<fault_array> > tab_fault_arrays_;

	int index_;
	int seed_;



public:
	
	// constructeur, destructeur
	fault_arrays_set_geo(int index, int seed,std::vector<std::shared_ptr<fault_array> >& tab_fault_arrays);
	fault_arrays_set_geo(int index, int seed);
	fault_arrays_set_geo();
    ~fault_arrays_set_geo();

	void add_array(std::shared_ptr<fault_array> fault_array);

	bool check_array_center(Point3D& pt);

	void get_fracture_set_surface_mesh(std::vector<Point3D>& tab_point, std::vector<std::vector<int>>& tab_triangle, std::vector<property<float>>&  tab_triangle_fracture_property) ;
	void get_array_surface_mesh(std::vector<Point3D>& tab_point, std::vector<std::vector<int>>& tab_triangle ,  std::vector<property<float>>& tab_triangle_poly_index);
    void get_fracture_set_protection_zone_surface_mesh(std::vector<Point3D>& tab_point, std::vector<std::vector<int>>& tab_triangle,std::vector<property<float>>& tab_triangle_fault_array_property);
	
};

}

#endif



