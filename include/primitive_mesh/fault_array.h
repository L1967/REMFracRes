#include <primitive_mesh/boite.h>
#include <primitive_mesh/faille.h>
#include <primitive_mesh/fracture_set_geo.h>
#include <numeric/property.h>




#if !defined(Fault_array)
#define Fault_array

namespace fractures_intersect {

//-------------------------------------------------------------
// class fault_array
//-------------------------------------------------------------

class fault_array
{
public:

	// variables

	std::vector< std::shared_ptr<fracture_set_geo> > tab_fractures_set_; // tableau des primitives 3D non inclue
	std::shared_ptr<Cboite > box_3d_;
	std::shared_ptr<oriented_bounding_box> oriented_bb_;
	std::shared_ptr<oriented_bounding_box> oriented_bb_protection_;
	std::vector<std::shared_ptr<Cfaille> > tab_intersected_fractures_;
	double azimut_;
	double pendage_;
	double ext2L_selon_lpgp_;
	double ext2l_orthog_lpgp_ ;
	double thickness_;
	Point3D pt_ref_ ;

	int array_set_index_;

	double protection_distance_;

	int index_;
	int seed_;
	double min_space_distance_;
	double damage_zone_width_;
	double ratio_footwall_;
	double ratio_hangingwall_;


public:
	
	fault_array(int array_set_index = 0, Point3D pt_ref = Point3D(), double azim = 0.,
			double pend = 0.,double ext_selon_lpgp = 0.,
			double ext_ortho_lpgp = 0., double thickness = 0.0, double protection_distance_ = 0.0);

    fault_array() {};
    ~fault_array() {};

    Point3D generate_point(double random_x, double random_y, double random_z);
    bool contains_in_protection_zone(Point3D& pt);

	void get_fracture_set_surface_mesh(std::vector<Point3D>& tab_point, std::vector<std::vector<int>>& tab_triangle, std::vector<property<float>>&  tab_triangle_fracture_property) ;
	void get_array_box_surface_mesh(std::vector<Point3D>& tab_point, std::vector<std::vector<int>>& tab_triangle , std::vector<property<float>>& tab_triangle_poly_index);
	void intersect_all_fracture_set(int number_of_thread);
	void set_minimum_space_between_array( double l);
	void set_damage_width( double l);
	void set_ratio_footwall( double l);
	void set_ratio_hangingwall( double l);
	void check_and_intersected_fracture(std::shared_ptr<Cfaille> fracture);
	bool segment_point_array_center_intersect_fracture(Point3D& pt);
	void update_box_from_damage_zone();
	
};

}


#endif
