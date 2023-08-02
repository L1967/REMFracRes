

#include <primitive_mesh/fault_arrays_set_geo.h>
#include <iomanip>
#include <fstream>
#include <memory>

namespace fractures_intersect {
//-------------------------------------------------------------
// class fault_array_set_geo
//
// fracture set data
//-------------------------------------------------------


fault_arrays_set_geo::fault_arrays_set_geo()
{
	index_ = -1;
	seed_ = 12548;
}

fault_arrays_set_geo::fault_arrays_set_geo(int index, int seed)
{
	index_ = index;
	seed_ = seed;
}

fault_arrays_set_geo::fault_arrays_set_geo(int index, int seed, std::vector<std::shared_ptr<fault_array> >& tab_fault_arrays) {
	index_ = index;
	seed_ = seed;
	tab_fault_arrays_ = tab_fault_arrays;
}

fault_arrays_set_geo::~fault_arrays_set_geo()
{
}

void fault_arrays_set_geo::add_array(std::shared_ptr<fault_array> fault_array) {
	tab_fault_arrays_.push_back(fault_array);
}

bool fault_arrays_set_geo::check_array_center(Point3D& pt){

	for (size_t fa = 0 ; fa < tab_fault_arrays_.size(); fa++)
	{
		if (tab_fault_arrays_[fa]->oriented_bb_protection_->contains(pt)) return false;
	}
	return true;
}

void fault_arrays_set_geo::get_fracture_set_surface_mesh(std::vector<Point3D>& tab_point, std::vector<std::vector<int>>& tab_triangle,
		std::vector<property<float>>& tab_triangle_fault_array_property)
{

	std::vector<std::vector<float>> tab_triangle_properties;

	for (size_t fa = 0 ; fa < tab_fault_arrays_.size(); fa++)
	{

		for (size_t fi = 0 ; fi < tab_fault_arrays_[fa]->tab_fractures_set_.size(); fi++)
		{
			int fracture_count = 0;
			for (size_t i = 0 ; i < tab_fault_arrays_[fa]->tab_fractures_set_[fi]->tab_fractures_.size(); i++)
			{
				std::vector<float> property;
				property.push_back(fa);
				property.push_back(fi);
				property.push_back(fracture_count);
				property.push_back(tab_fault_arrays_[fa]->tab_fractures_set_[fi]->tab_fractures_[i]->azimut-90);
				property.push_back(tab_fault_arrays_[fa]->tab_fractures_set_[fi]->tab_fractures_[i]->pendage);
				property.push_back(tab_fault_arrays_[fa]->tab_fractures_set_[fi]->tab_fractures_[i]->ext2L_selon_lpgp);
				property.push_back(tab_fault_arrays_[fa]->tab_fractures_set_[fi]->tab_fractures_[i]->ext2l_orthog_lpgp);
				property.push_back(tab_fault_arrays_[fa]->tab_fractures_set_[fi]->tab_fractures_[i]->aperture_);
				property.push_back(tab_fault_arrays_[fa]->tab_fractures_set_[fi]->tab_fractures_[i]->stress_zone_max_);
				property.push_back(tab_fault_arrays_[fa]->tab_fractures_set_[fi]->tab_fractures_[i]->connection_relay_);
				property.push_back(tab_fault_arrays_[fa]->tab_fractures_set_[fi]->tab_fractures_[i]->familly_index_);

				tab_fault_arrays_[fa]->tab_fractures_set_[fi]->tab_fractures_[i]->add_mesh_data(tab_point, tab_triangle, tab_triangle_properties, property);
				fracture_count++;
			}
		}
	}

	size_t tr_size = tab_triangle.size();


	tab_triangle_fault_array_property.push_back(property<float>("fault_array_set_index",tr_size));
	tab_triangle_fault_array_property.push_back(property<float>("fracture_set_index",tr_size));
	tab_triangle_fault_array_property.push_back(property<float>("fault_index",tr_size));
	tab_triangle_fault_array_property.push_back(property<float>("azimut",tr_size));
	tab_triangle_fault_array_property.push_back(property<float>("dip",tr_size));
	tab_triangle_fault_array_property.push_back(property<float>("length_1",tr_size));
	tab_triangle_fault_array_property.push_back(property<float>("length_2",tr_size));
	tab_triangle_fault_array_property.push_back(property<float>("aperture",tr_size));
	tab_triangle_fault_array_property.push_back(property<float>("stress_zone",tr_size));
	tab_triangle_fault_array_property.push_back(property<float>("connection_relay",tr_size));
	tab_triangle_fault_array_property.push_back(property<float>("familly_index",tr_size));



	for (size_t i = 0 ; i < tr_size ; i++)
	{
		for (size_t j = 0 ; j < tab_triangle_fault_array_property.size() ; j++)
		{
			tab_triangle_fault_array_property[j].set_value(i,tab_triangle_properties[i][j]);
		}

	}





}




void fault_arrays_set_geo::get_fracture_set_protection_zone_surface_mesh(std::vector<Point3D>& tab_point, std::vector<std::vector<int>>& tab_triangle,
		std::vector<property<float>>& tab_triangle_fault_array_property)
{

	std::vector<std::vector<float>> tab_triangle_properties;

	for (size_t fa = 0 ; fa < tab_fault_arrays_.size(); fa++)
	{

		for (size_t fi = 0 ; fi < tab_fault_arrays_[fa]->tab_fractures_set_.size(); fi++)
		{
			int fracture_count = 0;
			for (size_t i = 0 ; i < tab_fault_arrays_[fa]->tab_fractures_set_[fi]->tab_stress_zone_cylinders_.size(); i++)
			{
				std::vector<float> property;
				property.push_back(fa);
				property.push_back(fi);
				property.push_back(fracture_count);
				property.push_back(tab_fault_arrays_[fa]->tab_fractures_set_[fi]->tab_stress_zone_cylinders_[i]->faille_support->azimut - 90);
				property.push_back(tab_fault_arrays_[fa]->tab_fractures_set_[fi]->tab_stress_zone_cylinders_[i]->faille_support->pendage);
				property.push_back(tab_fault_arrays_[fa]->tab_fractures_set_[fi]->tab_stress_zone_cylinders_[i]->faille_support->ext2L_selon_lpgp);
				property.push_back(tab_fault_arrays_[fa]->tab_fractures_set_[fi]->tab_stress_zone_cylinders_[i]->faille_support->ext2l_orthog_lpgp);
				property.push_back(tab_fault_arrays_[fa]->tab_fractures_set_[fi]->tab_stress_zone_cylinders_[i]->faille_support->aperture_);
				property.push_back(tab_fault_arrays_[fa]->tab_fractures_set_[fi]->tab_stress_zone_cylinders_[i]->faille_support->stress_zone_max_);
				property.push_back(tab_fault_arrays_[fa]->tab_fractures_set_[fi]->tab_stress_zone_cylinders_[i]->faille_support->familly_index_);


				for (size_t j = 0; j < tab_fault_arrays_[fa]->tab_fractures_set_[fi]->tab_stress_zone_cylinders_[i]->tab_plan.size() ; j++)
				{
					tab_fault_arrays_[fa]->tab_fractures_set_[fi]->tab_stress_zone_cylinders_[i]->tab_plan[j]->add_mesh_data(tab_point, tab_triangle, tab_triangle_properties, property);
				}


				fracture_count++;
			}
		}
	}

	size_t tr_size = tab_triangle.size();


	tab_triangle_fault_array_property.push_back(property<float>("fault_array_set_index",tr_size));
	tab_triangle_fault_array_property.push_back(property<float>("fracture_set_index",tr_size));
	tab_triangle_fault_array_property.push_back(property<float>("fault_index",tr_size));
	tab_triangle_fault_array_property.push_back(property<float>("azimut",tr_size));
	tab_triangle_fault_array_property.push_back(property<float>("dip",tr_size));
	tab_triangle_fault_array_property.push_back(property<float>("length_1",tr_size));
	tab_triangle_fault_array_property.push_back(property<float>("length_2",tr_size));
	tab_triangle_fault_array_property.push_back(property<float>("aperture",tr_size));
	tab_triangle_fault_array_property.push_back(property<float>("stress_zone",tr_size));
	tab_triangle_fault_array_property.push_back(property<float>("familly_index",tr_size));



	for (size_t i = 0 ; i < tr_size ; i++)
	{
		for (size_t j = 0 ; j < tab_triangle_fault_array_property.size() ; j++)
		{
			tab_triangle_fault_array_property[j].set_value(i,tab_triangle_properties[i][j]);
		}

	}





}












void fault_arrays_set_geo::get_array_surface_mesh(std::vector<Point3D>& tab_point,
		std::vector<std::vector<int>>& tab_triangle, std::vector<property<float>>& tab_triangle_poly_index)
{

	std::vector<std::vector<float>> tab_triangle_properties;


	for (size_t fa = 0 ; fa < tab_fault_arrays_.size(); fa++)
	{

		std::vector<float> property;
		property.push_back(index_);
		property.push_back(fa);
		property.push_back(tab_fault_arrays_[fa]->azimut_-90);
		property.push_back(tab_fault_arrays_[fa]->pendage_);
		property.push_back(tab_fault_arrays_[fa]->ext2L_selon_lpgp_);
		property.push_back(tab_fault_arrays_[fa]->ext2l_orthog_lpgp_);
		property.push_back(tab_fault_arrays_[fa]->thickness_);
		property.push_back(tab_fault_arrays_[fa]->protection_distance_);


		for (size_t j = 0; j < tab_fault_arrays_[fa]->box_3d_->tab_plan.size() ; j++)
		{
			tab_fault_arrays_[fa]->box_3d_->tab_plan[j]->add_mesh_data(tab_point, tab_triangle, tab_triangle_properties, property);
		}

		tab_fault_arrays_[fa]->box_3d_->Get_plan_ref()->add_mesh_data(tab_point, tab_triangle, tab_triangle_properties, property);
	}

	size_t tr_size = tab_triangle.size();

	tab_triangle_poly_index.push_back(property<float>("fault_array_set_index",tr_size));
	tab_triangle_poly_index.push_back(property<float>("fault_array_index",tr_size));
	tab_triangle_poly_index.push_back(property<float>("azimut",tr_size));
	tab_triangle_poly_index.push_back(property<float>("dip",tr_size));
	tab_triangle_poly_index.push_back(property<float>("length_1",tr_size));
	tab_triangle_poly_index.push_back(property<float>("length_2",tr_size));
	tab_triangle_poly_index.push_back(property<float>("thickness",tr_size));
	tab_triangle_poly_index.push_back(property<float>("protection_distance",tr_size));



	for (size_t i = 0 ; i < tr_size ; i++)
	{
		for (size_t j = 0 ; j < tab_triangle_poly_index.size() ; j++)
		{
			tab_triangle_poly_index[j].set_value(i,tab_triangle_properties[i][j]);
		}

	}
}



}


