
#include <primitive_mesh/fault_array.h>
#include <iomanip>
#include <fstream>
#include <memory>

namespace fractures_intersect {


//-------------------------------------------------------------
// class fault_array
//-------------------------------------------------------


fault_array::fault_array(int array_set_index, Point3D pt_ref,double azim, double pend,double ext_selon_lpgp, double ext_ortho_lpgp,
		double thickness, double protection_distance)
{
	array_set_index_ = array_set_index;
	azimut_ = azim;
	pendage_ = pend;
	ext2L_selon_lpgp_ = ext_selon_lpgp;
	ext2l_orthog_lpgp_ = ext_ortho_lpgp;
	thickness_ =thickness;
	protection_distance_ = protection_distance;

	pt_ref_ = pt_ref;

	box_3d_ = std::shared_ptr<Cboite >(new Cboite(pt_ref_,azim,pend,
		ext_selon_lpgp,ext_ortho_lpgp, thickness));


	Cfaille faille(0, pt_ref_,azim,pend,ext_selon_lpgp,ext_ortho_lpgp, thickness);

	oriented_bb_ = std::shared_ptr<oriented_bounding_box>( new oriented_bounding_box(faille));
	//	std::shared_ptr<oriented_bounding_box> oriented_bb_protection_;

	oriented_bb_protection_ = std::shared_ptr<oriented_bounding_box>(new oriented_bounding_box(faille));


	oriented_bb_protection_->extend(std::max(ext_selon_lpgp, ext_ortho_lpgp)/2.0 + protection_distance);

	min_space_distance_ = 0.1;

	damage_zone_width_ = 0.5;
	ratio_footwall_ = 0.5;
	ratio_hangingwall_=0.5;
}

void fault_array::update_box_from_damage_zone(){
	box_3d_->Actualise_boite_damage_zone(ratio_footwall_,ratio_hangingwall_);
}
void fault_array::set_damage_width( double l){
	damage_zone_width_ =l;
}
void fault_array::set_ratio_footwall( double l){
	ratio_footwall_ = l;
}
void fault_array::set_ratio_hangingwall( double l){
	ratio_hangingwall_ = l;
}
void fault_array::set_minimum_space_between_array(double l){
	min_space_distance_ = l;
}


Point3D fault_array::generate_point(double random_x, double random_y, double random_z) {
	return oriented_bb_->generate_point(random_x, random_y, random_z);
}

bool fault_array::contains_in_protection_zone(Point3D& pt) {
	return oriented_bb_protection_->contains(pt);
}
void fault_array::check_and_intersected_fracture(std::shared_ptr<Cfaille> fracture)
{

	if (fracture->distance_to_plan(pt_ref_) > oriented_bb_->diagonal()/2.0) return;

	if (box_3d_->Intersect_Bounding_box(fracture.get()) == false) return;

	for( Point3D pt : fracture->tab_primitives[0]->tab_point) {
		if (oriented_bb_->contains(pt)) {
			tab_intersected_fractures_.push_back(fracture);
			return;
		}
	}

	bool poss = false;
	bool neg = false;
	bool intersect = false;

	for (int i = 0; i < 2 ; i++) {
		for(Point3D pt: box_3d_->tab_plan[i]->tab_primitives[0]->tab_point) {
			double dist = fracture->oriented_distance_to_plan(pt);
			if (dist > 0) poss = true;
			if (dist < 0) neg = true;
			if (poss && neg) {
				intersect = true;
				break;
			}
		}
		if (intersect) break;

	}

	if (intersect == false) return;


	for(std::shared_ptr<Cplan> plane: box_3d_->tab_plan) {

		if( plane->check_intersect_plan(fracture.get(),pt_tol) == true) {
			tab_intersected_fractures_.push_back(fracture);
			return;
		}
	}

}


bool fault_array::segment_point_array_center_intersect_fracture(Point3D& pt) {
	return true;

	for(std::shared_ptr<Cfaille> fracture: tab_intersected_fractures_) {

		double d1 = fracture->oriented_distance_to_plan(pt_ref_);
		double d2 = fracture->oriented_distance_to_plan(pt);

		if ((d1 > 0 && d2 > 0) || (d1 < 0 && d2 < 0)) continue;
		double t = d1 / (d1 - d2);

		Point3D intersection_pt(pt_ref_.coord_ + t * (pt.coord_ - pt_ref_.coord_));

		std::shared_ptr<Cpolygone> polygon = std::shared_ptr<Cpolygone>(std::dynamic_pointer_cast<Cpolygone>(fracture->tab_primitives[0]));

		if (polygon->Point_ds_polygone(fracture.get(),intersection_pt,pt_tol,false)) return true;
	}

	return false;
}
//-------------------------------------------
// Operations relatives aux intersections
//-------------------------------------------
void fault_array::intersect_all_fracture_set(int number_of_thread)
{

	std::uniform_real_distribution<double> distribution(0.0,1.0);

	for (size_t si = 0 ; si < tab_fractures_set_.size() ; si++)
	{

		std::mt19937 generator(tab_fractures_set_[si]->seed_);

//		#pragma  omp  parallel for schedule(dynamic,5)  num_threads(number_of_thread) if (tab_fractures_set_[si]->tab_fractures_.size() > 10)
		for (size_t fi  = 0 ; fi < tab_fractures_set_[si]->tab_fractures_.size() ; fi++)
		{
			Cfaille* fracture_ref  = tab_fractures_set_[si]->tab_fractures_[fi].get();

			for (size_t sj = 0 ; sj < si ; sj++)
			{
				bool random = false;

				if (tab_fractures_set_[si]->tab_type_intersection_[tab_fractures_set_[sj]->set_type_] == CROSSING) continue;
				if (tab_fractures_set_[si]->tab_type_intersection_[tab_fractures_set_[sj]->set_type_] == RANDOM) {
					random = true;
				}

				for (size_t fj = 0 ; fj < tab_fractures_set_[sj]->tab_fractures_.size() ; fj++)
				{
					Cfaille* fracture= tab_fractures_set_[sj]->tab_fractures_[fj].get();
					double prob = 1.0;
					if (random) prob = distribution(generator);

					if (prob > 0.5) {
						fracture_ref->intersect_plan_limitated(fracture,pt_tol);
					}
				}
			}

		}

	}
}

void fault_array::get_fracture_set_surface_mesh(std::vector<Point3D>& tab_point, std::vector<std::vector<int>>& tab_triangle,
		std::vector<property<float>>& tab_triangle_fracture_property)
{
	std::vector<std::vector<float>> tab_triangle_properties;

	for (size_t fi = 0 ; fi < tab_fractures_set_.size(); fi++)
	{
		int fracture_count = 0;
		for (size_t i = 0 ; i < tab_fractures_set_[fi]->tab_fractures_.size(); i++)
		{
			std::vector<float> property;
			property.push_back(fi);
			property.push_back(fracture_count);
			property.push_back(tab_fractures_set_[fi]->tab_fractures_[i]->azimut);
			property.push_back(tab_fractures_set_[fi]->tab_fractures_[i]->pendage);
			property.push_back(tab_fractures_set_[fi]->tab_fractures_[i]->ext2L_selon_lpgp);
			property.push_back(tab_fractures_set_[fi]->tab_fractures_[i]->ext2l_orthog_lpgp);
			property.push_back(tab_fractures_set_[fi]->tab_fractures_[i]->aperture_);
			property.push_back(tab_fractures_set_[fi]->tab_fractures_[i]->stress_zone_max_);

			tab_fractures_set_[fi]->tab_fractures_[i]->add_mesh_data(tab_point, tab_triangle, tab_triangle_properties, property);
			fracture_count++;
		}
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



	for (size_t i = 0 ; i < tr_size ; i++)
	{
		for (size_t j = 0 ; j < tab_triangle_fracture_property.size() ; j++)
		{
			tab_triangle_fracture_property[j].set_value(i,tab_triangle_properties[i][j]);
		}

	}

}

}


