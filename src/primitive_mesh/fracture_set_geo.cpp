
#include <primitive_mesh/fracture_set_geo.h>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <memory>

namespace fractures_intersect {

using namespace std;
const long MAXCAR = 120;
//-------------------------------------------------------------
// class fracture_set_geo
//
// fracture set data
//-------------------------------------------------------


fracture_set_geo::fracture_set_geo()
{
	index_ = -1;
	set_type_ = FracResIntersectionCategoryEnum::FAULTS;
	seed_ = 12548;
}

fracture_set_geo::fracture_set_geo(int index, FracResIntersectionCategoryEnum set_type, int seed, std::vector<FracResIntersectionTypeEnum>& tab_type_intersection)
{
	index_ = index;
	seed_ = seed;
	set_type_ = set_type;
	tab_type_intersection_ = tab_type_intersection;
}

fracture_set_geo::fracture_set_geo(int index, FracResIntersectionCategoryEnum set_type, int seed, std::vector<FracResIntersectionTypeEnum>& tab_type_intersection ,
    		std::vector<std::shared_ptr<Cfaille> >& tab_fractures) {
	index_ = index;
	seed_ = seed;
	set_type_ = set_type;
	tab_type_intersection_ = tab_type_intersection;
	tab_fractures_ = tab_fractures;
}

fracture_set_geo::~fracture_set_geo()
{
}
void fracture_set_geo::add_bedding_localisation(std::shared_ptr<FracResBeddingLocalisation>  localisation){

	tab_bedding_localisation_.push_back(localisation);
}
void fracture_set_geo::add_fracture(std::shared_ptr<Cfaille> fracture) {
	tab_fractures_.push_back(fracture);
}

void fracture_set_geo::add_fracture_and_create_cylinder(std::shared_ptr<Cfaille> fracture,int n_plan) {
	tab_fractures_.push_back(fracture);
	tab_stress_zone_cylinders_.push_back(std::shared_ptr<Ccylindre>(new Ccylindre(fracture, n_plan)));
	
}

void fracture_set_geo::add_fracture_and_abbuting(std::shared_ptr<Cfaille> fracture_add) {



	for (size_t fj = 0 ; fj < tab_fractures_.size() ; fj++)
	{

		if(fj==14 && tab_fractures_.size()==17){
			int ducon=0;
			ducon++;
		}
		Cfaille* fracture= tab_fractures_[fj].get();
		bool intersection = fracture_add->intersect_plan_limitated(fracture,pt_tol);


		string result = intersection == 0 ? "false" : "true";

		if(intersection > 0) std::cout<< "fracture " << fj <<" intersection with " << tab_fractures_.size() << " is " << result << std::endl;
	}



	tab_fractures_.push_back(fracture_add);
}



void fracture_set_geo::clear_stress_zone_cylinder() {
	tab_stress_zone_cylinders_.clear();
}

//-------------------------------------------
// check all plan stress  zone
//-------------------------------------------
bool fracture_set_geo::in_all_stress_zone(Point3D& pt)
{

	for (size_t i = 0 ; i < tab_fractures_.size(); i++)
	{
		if (tab_fractures_[i]->in_stress_zone(pt) == true) return true;
	}

	return false;
}


//-------------------------------------------
// check all plan stress  zone
//-------------------------------------------
bool fracture_set_geo::in_all_stress_zone_3d(std::shared_ptr<Cfaille> faille)
{

	std::shared_ptr<Ccylindre> cyl_faille = std::shared_ptr<Ccylindre>(new Ccylindre(faille,12));

	for (size_t i = 0 ; i < tab_stress_zone_cylinders_.size() ; i++)
	{
		if (tab_stress_zone_cylinders_[i]->faille_support->in_stress_zone(cyl_faille->faille_support->Get_point_ref()) == true) return true;
		if (tab_stress_zone_cylinders_[i]->check_possible_intersect_polyedre(cyl_faille.get()) == true) return true;
	}

	return false;
}

//-------------------------------------------
// check all plan stress  zone
//-------------------------------------------
bool fracture_set_geo::in_all_stress_zone__abutting_3d(std::shared_ptr<Cfaille> faille, int n_plan)
{


	for (size_t i = 0 ; i < tab_stress_zone_cylinders_.size() ; i++)
	{
		if (tab_stress_zone_cylinders_[i]->faille_support->in_stress_zone(faille->Get_point_ref()) == true) return true;
		if (faille->intersect_plan_limitated(tab_stress_zone_cylinders_[i]->faille_support.get(),pt_tol) ) {
			faille->simplify_to_box_poly();
		}

		tab_stress_zone_cylinders_[i]->abutting_fracture(faille.get());
	}


	std::shared_ptr<Ccylindre> new_cylinder = std::shared_ptr<Ccylindre>(new Ccylindre(faille, n_plan));
	//new_cylinder->Actualise_facettes(new_cylinder->nb_facettes);

	for (size_t i = 0 ; i < tab_stress_zone_cylinders_.size() ; i++)
	{
		std::shared_ptr<Cfaille> faille_cyl =  tab_stress_zone_cylinders_[i]->faille_support;
		if (new_cylinder->abutting_fracture(faille_cyl.get()) == true) {
			tab_stress_zone_cylinders_[i]->Actualise_facettes(tab_stress_zone_cylinders_[i]->nb_facettes);

		}
	}




	tab_fractures_.push_back(faille);
	tab_stress_zone_cylinders_.push_back(new_cylinder);


	return false;
}



void fracture_set_geo::get_fracture_set_surface_mesh(std::vector<Point3D>& tab_point, std::vector<std::vector<int>>& tab_triangle,
		std::vector<property<float>>& tab_triangle_fracture_property)
{

	int fracture_count = 0;

	std::vector<std::vector<float>> tab_triangle_properties;

	for (size_t i = 0 ; i < tab_fractures_.size(); i++)
	{
		std::vector<float> property;
		property.push_back(index_);
		property.push_back(fracture_count);
		property.push_back(tab_fractures_[i]->azimut);
		property.push_back(tab_fractures_[i]->pendage);
		property.push_back(tab_fractures_[i]->ext2L_selon_lpgp);
		property.push_back(tab_fractures_[i]->ext2l_orthog_lpgp);
		property.push_back(tab_fractures_[i]->aperture_);
		property.push_back(tab_fractures_[i]->stress_zone_max_);
		property.push_back(tab_fractures_[i]->familly_index_);

		tab_fractures_[i]->add_mesh_data(tab_point, tab_triangle, tab_triangle_properties, property);
		fracture_count++;
	}

	size_t tr_size = tab_triangle.size();


	tab_triangle_fracture_property.push_back(property<float>("fracture_set_index",tr_size));
	tab_triangle_fracture_property.push_back(property<float>("fracture_index",tr_size));
	tab_triangle_fracture_property.push_back(property<float>("azimut",tr_size));
	tab_triangle_fracture_property.push_back(property<float>("dip",tr_size));
	tab_triangle_fracture_property.push_back(property<float>("length_1",tr_size));
	tab_triangle_fracture_property.push_back(property<float>("length_2",tr_size));
	tab_triangle_fracture_property.push_back(property<float>("aperture",tr_size));
	tab_triangle_fracture_property.push_back(property<float>("stress_zone",tr_size));
	tab_triangle_fracture_property.push_back(property<float>("familly_index",tr_size));


	for (size_t i = 0 ; i < tr_size ; i++)
	{
		for (size_t j = 0 ; j < tab_triangle_fracture_property.size() ; j++)
		{
			tab_triangle_fracture_property[j].set_value(i,tab_triangle_properties[i][j]);
		}

	}





}

void fracture_set_geo::get_fracture_set_protection_zone_surface_mesh(std::vector<Point3D>& tab_point,
		std::vector<std::vector<int>>& tab_triangle, std::vector<property<float>>& tab_triangle_poly_index)
{

	std::vector<std::vector<float>> tab_triangle_properties;

	int poly_count = 0;
	for (size_t i = 0 ; i < tab_stress_zone_cylinders_.size() ; i++)
	{
		std::vector<float> property;
		property.push_back(index_);
		property.push_back(poly_count);
		property.push_back(tab_stress_zone_cylinders_[i]->faille_support->azimut);
		property.push_back(tab_stress_zone_cylinders_[i]->faille_support->pendage);
		property.push_back(tab_stress_zone_cylinders_[i]->faille_support->ext2L_selon_lpgp);
		property.push_back(tab_stress_zone_cylinders_[i]->faille_support->ext2l_orthog_lpgp);
		property.push_back(tab_stress_zone_cylinders_[i]->faille_support->aperture_);
		property.push_back(tab_stress_zone_cylinders_[i]->faille_support->stress_zone_max_);
		property.push_back(tab_stress_zone_cylinders_[i]->faille_support->familly_index_);


		for (size_t j = 0; j < tab_stress_zone_cylinders_[i]->tab_plan.size() ; j++)
		{
			tab_stress_zone_cylinders_[i]->tab_plan[j]->add_mesh_data(tab_point, tab_triangle, tab_triangle_properties, property);
		}
		poly_count++;
	}

	size_t tr_size = tab_triangle.size();

	tab_triangle_poly_index.push_back(property<float>("fracture_set_index",tr_size));
	tab_triangle_poly_index.push_back(property<float>("fracture_index",tr_size));
	tab_triangle_poly_index.push_back(property<float>("azimut",tr_size));
	tab_triangle_poly_index.push_back(property<float>("dip",tr_size));
	tab_triangle_poly_index.push_back(property<float>("length_1",tr_size));
	tab_triangle_poly_index.push_back(property<float>("length_2",tr_size));
	tab_triangle_poly_index.push_back(property<float>("aperture",tr_size));
	tab_triangle_poly_index.push_back(property<float>("stress_zone",tr_size));
	tab_triangle_poly_index.push_back(property<float>("familly_index",tr_size));



	for (size_t i = 0 ; i < tr_size ; i++)
	{
		for (size_t j = 0 ; j < tab_triangle_poly_index.size() ; j++)
		{
			tab_triangle_poly_index[j].set_value(i,tab_triangle_properties[i][j]);
		}

	}
}



}


