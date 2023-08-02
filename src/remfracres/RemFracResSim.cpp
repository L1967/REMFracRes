/*[
 * RemFracRes Software.
 *
 * Copyright, TotalEnergies, 
 * 2002-2022, All rights reserved.
 *
 * The software and information contained herein are proprietary to, and
 * comprise valuable trade secrets of, TotalEnergies 
 * which intend to preserve as trade secrets such software and information.
 * This software is furnished pursuant to a written license agreement and
 * may be used, copied, transmitted, and stored only in accordance with
 * the terms of such license and with the inclusion of the above copyright
 * notice.  This software and information or any other copies thereof may
 * not be provided or otherwise made available to any other person.
 *
 ]*/

//----------------------------RemFracResSim.cpp  ---------------------------
// 
//  RemFracResSim contains the parameters for the simulation 
//  of a Boolean Abutting and Anstomose model and the generated fractures
//  
//
//  Author:  Pascal Siegel - Rabah Namar - Olivier Jaquet 
//
//	Date: 2002-2022
//
//----------------------------  RemFracResSim.cpp  ---------------------------

// local includes
#include <remfracres/RemFracResSim.h>
#include <primitive_mesh/varglob.h>
#include <primitive_mesh/contrainte_geo.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <memory>
#include <string>
#include <array>
#include <chrono>
#include <ctime>
#include <cmath>
#include <sys/stat.h>


// lib includes 

#include <geode/mesh/core/solid_mesh.h>
#include <geode/mesh/core/tetrahedral_solid.h>
#include <geode/geometry/basic_objects/tetrahedron.h>
#include <geode/geometry/basic_objects/triangle.h>
#include <geode/model/mixin/core/block.h>
#include <geode/geometry/mensuration.h>
#include <geode/mesh/helpers/convert_surface_mesh.h>
#include <geode/basic/attribute_manager.h>
#include <geode/mesh/builder/solid_edges_builder.h>
#include <geode/mesh/builder/solid_facets_builder.h>
#include <geode/mesh/core/geode/geode_tetrahedral_solid.h>
#include <geode/mesh/core/solid_edges.h>
#include <geode/mesh/core/solid_facets.h>
#include <geode/mesh/core/geode/geode_triangulated_surface.h>
#include <geode/geometry/vector.h>
#include <geode/geometry/mensuration.h>
#include <geode/geometry/aabb.h>
#include <geode/mesh/helpers/aabb_surface_helpers.h>
#include <geode/mesh/helpers/ray_tracing.h>
#include <geode/mesh/builder/surface_mesh_builder.h>
#include <geode/mesh/helpers/convert_surface_mesh.h>
#include <geode/basic/logger.h>
#include <geode/basic/range.h>
#include <geode/geometry/bounding_box.h>
#include <geode/model/helpers/component_mesh_polygons.h>
#include <geode/model/helpers/convert_to_mesh.h>
#include <geode/mesh/helpers/convert_solid_mesh.h>
#include <geode/mesh/io/tetrahedral_solid_output.h>

#include <geode/mesh/core/geode/geode_polygonal_surface.h>
#include <geode/mesh/core/polygonal_surface.h>
#include <geode/mesh/core/triangulated_surface.h>
#include <geode/mesh/io/triangulated_surface_input.h>
#include <geode/mesh/io/triangulated_surface_output.h>
#include <geode/mesh/io/polygonal_surface_input.h>
#include <geode/mesh/io/polygonal_surface_output.h>
#include <geode/mesh/builder/geode/geode_triangulated_surface_builder.h>
#include <geode/mesh/builder/geode/geode_polygonal_surface_builder.h>
#include <geode/mesh/builder/triangulated_surface_builder.h>
#include <geode/mesh/builder/polygonal_surface_builder.h>
//absl
#include <absl/container/flat_hash_map.h>
#include <absl/container/flat_hash_set.h>
#include <absl/strings/str_split.h>
#include <absl/strings/string_view.h>
//#include <boost/random.hpp>
#include <boost/random.hpp>
#include <boost/math/distributions/uniform.hpp>
#include <boost/math/distributions.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/triangle_distribution.hpp>
#include <boost/random/lognormal_distribution.hpp>
#include <boost/math/distributions/lognormal.hpp>
#include <boost/math/distributions/complement.hpp>
// system includes                       

namespace remfracres {

using namespace geode;
using namespace std;
using namespace fractures_intersect;
using namespace absl;

//-----------------------------------------------
// default constructor of the class RemFracResSim
//-----------------------------------------------
/**
 * @fn  RemFracResSim()
 * @brief
 *
 */
RemFracResSim::RemFracResSim() {
	name_ = "simultion_0";
	bed_interface_property_name_ = "None";
	outpout_file_name_ ="output";
	method_type_ = 0;
	volume_ = 0;
	area_ = 0;
	model_fractures_ = fractures_intersect::model();
	fracturesSet_list_.clear();
	use_global_outpout_directory_path_ = false;
	outpout_directory_path_ = remfracres::data_path;
	intersect_inside_set_ = false;
	success_ = false;
	cylinder_step_ = 12;
	nb_failles_ = 0;
	seed_ = 12345;
	nb_fault_arrays_set_ = 0;
	manage_bed_interface_ = false;
	nb_total_bed_ = 0;
	nb_realisation_=1;
	fracresFileKeyWord_ = new FracResFileKeyWord();

	fileManagement_ = new RemFracResFileManagement();
	fracResUnit_ = new FracResUnit();
	fracResStageList_.clear();
	fracResScenarioList_.clear();
	stress_geo_list_.clear();
	stress_geo_parallel_list_.clear();
	bed_name_vect_.clear();
	bed_name_parallel_vect_.clear();
	by_scenario_ = false;
	total_simulation_time_ = 0.0;
	total_writting_output_time_ = 0.0;
	active_discrete_pair_array_.clear();


}

//---------------------------------------------
//  copy constructor of the class RemFracResSim
//---------------------------------------------
/**
 * @fn  RemFracResSim(const RemFracResSim&)
 * @brief
 *
 * @param from
 */
RemFracResSim::RemFracResSim(const RemFracResSim &from) {
	name_ = from.name_;
	method_type_ = from.method_type_;
	volume_ = from.volume_;
	area_ = from.area_;
	model_fractures_ = from.model_fractures_;
	fracturesSet_list_ = from.fracturesSet_list_;
	use_global_outpout_directory_path_ =
			from.use_global_outpout_directory_path_;
	outpout_directory_path_ = from.outpout_directory_path_;
	intersect_inside_set_ = from.intersect_inside_set_;
	success_ = from.success_;
	cylinder_step_ = from.cylinder_step_;
	nb_failles_ = from.nb_failles_;
	seed_ = from.seed_;
	nb_fault_arrays_set_ = from.nb_fault_arrays_set_;
	faultArraySet_list_ = from.faultArraySet_list_;
	fault_array_zone_set_list_ = from.fault_array_zone_set_list_;
	fault_array_type_intersection_ = from.fault_array_type_intersection_;

	manage_bed_interface_ = from.manage_bed_interface_;
	fracresFileKeyWord_ = from.fracresFileKeyWord_;
	fracResUnit_ = from.fracResUnit_;
	fracResStageList_.clear();
	fracResScenarioList_.clear();
	stress_geo_list_ = from.stress_geo_list_;
	stress_geo_parallel_list_ = from.stress_geo_parallel_list_;
	bed_interface_property_name_ = from.bed_interface_property_name_;
	fileManagement_ = from.fileManagement_;
	nb_total_bed_ = from.nb_total_bed_;
	bed_name_vect_ = from.bed_name_vect_;
	bed_name_parallel_vect_ = from.bed_name_parallel_vect_;
	nb_realisation_= from.nb_realisation_;
	fracResScenarioList_ = from.fracResScenarioList_;
	by_scenario_ = from.by_scenario_;
	total_simulation_time_ = from.total_simulation_time_;
	total_writting_output_time_ = from.total_writting_output_time_;
	outpout_file_name_ = from.outpout_file_name_;
	active_discrete_pair_array_ = from.active_discrete_pair_array_;

}

//-------------------------------------------
// destructor of the class RemFracResSim
//-------------------------------------------
/**
 * @fn  ~RemFracResSim()
 * @brief
 *
 */
RemFracResSim::~RemFracResSim() {
	fracResStageList_.clear();
	fracResScenarioList_.clear();

}

//-----------------------------------------------
// assignment operator of the class RemFracResSim
//----------------------------------------------
/**
 * @fn RemFracResSim operator =&(const RemFracResSim&)
 * @brief
 *
 * @param from
 * @return
 */
RemFracResSim& RemFracResSim::operator=(const RemFracResSim &from) {
	name_ = from.name_;

	method_type_ = from.method_type_;
	volume_ = from.volume_;
	area_ = from.area_;
	model_fractures_ = fractures_intersect::model();
	fracturesSet_list_ = from.fracturesSet_list_;
	use_global_outpout_directory_path_ =
			from.use_global_outpout_directory_path_;
	outpout_directory_path_ = from.outpout_directory_path_;
	intersect_inside_set_ = from.intersect_inside_set_;
	success_ = from.success_;
	cylinder_step_ = from.cylinder_step_;
	nb_failles_ = from.nb_failles_;
	seed_ = from.seed_;
	nb_fault_arrays_set_ = from.nb_fault_arrays_set_;
	faultArraySet_list_ = from.faultArraySet_list_;
	fault_array_zone_set_list_ = from.fault_array_zone_set_list_;
	fault_array_type_intersection_ = from.fault_array_type_intersection_;

	manage_bed_interface_ = from.manage_bed_interface_;
	fracresFileKeyWord_ = from.fracresFileKeyWord_;
	fracResUnit_ = from.fracResUnit_;
	stress_geo_list_ = from.stress_geo_list_;
	stress_geo_parallel_list_ = from.stress_geo_parallel_list_;
	bed_interface_property_name_ = from.bed_interface_property_name_;
	fileManagement_ = from.fileManagement_;
	nb_total_bed_ = from.nb_total_bed_;
	bed_name_vect_ = from.bed_name_vect_;
	bed_name_parallel_vect_ = from.bed_name_parallel_vect_;
	nb_realisation_= from.nb_realisation_;
	fracResScenarioList_ = from.fracResScenarioList_;
	by_scenario_ = from.by_scenario_;
	total_simulation_time_ = from.total_simulation_time_;
	total_writting_output_time_ = from.total_writting_output_time_;
	outpout_file_name_ = from.outpout_file_name_;
	active_discrete_pair_array_ = from.active_discrete_pair_array_;
	return *this;
}

//--------------------------------------------------------
//      copy function 
//--------------------------------------------------------
/**
 * @fn void copy(const RemFracResSim&)
 * @brief
 *
 * @param from
 */
void RemFracResSim::copy(const RemFracResSim &from) {

	method_type_ = from.method_type_;
	volume_ = from.volume_;
	area_ = from.area_;
	model_fractures_ = fractures_intersect::model();
	fracturesSet_list_ = from.fracturesSet_list_;
	use_global_outpout_directory_path_ =
			from.use_global_outpout_directory_path_;
	outpout_directory_path_ = from.outpout_directory_path_;
	intersect_inside_set_ = from.intersect_inside_set_;
	success_ = from.success_;
	cylinder_step_ = from.cylinder_step_;
	nb_failles_ = from.nb_failles_;
	seed_ = from.seed_;
	nb_fault_arrays_set_ = from.nb_fault_arrays_set_;
	faultArraySet_list_ = from.faultArraySet_list_;
	fault_array_zone_set_list_ = from.fault_array_zone_set_list_;
	fault_array_type_intersection_ = from.fault_array_type_intersection_;
	manage_bed_interface_ = from.manage_bed_interface_;
	fracresFileKeyWord_ = from.fracresFileKeyWord_;
	fracResUnit_ = from.fracResUnit_;
	stress_geo_list_ = from.stress_geo_list_;
	stress_geo_parallel_list_ = from.stress_geo_parallel_list_;
	bed_interface_property_name_ = from.bed_interface_property_name_;
	fileManagement_ = from.fileManagement_;
	nb_total_bed_ = from.nb_total_bed_;
	bed_name_vect_ = from.bed_name_vect_;
	bed_name_parallel_vect_ = from.bed_name_parallel_vect_;
	nb_realisation_= from.nb_realisation_;
	fracResScenarioList_ = from.fracResScenarioList_;
	by_scenario_ = from.by_scenario_;
	total_simulation_time_ = from.total_simulation_time_;
	total_writting_output_time_ = from.total_writting_output_time_;
	outpout_file_name_ = from.outpout_file_name_;
	active_discrete_pair_array_ = from.active_discrete_pair_array_;
}

//--------------------------
// set the domain limits
//--------------------------
/**
 * @fn void set_generation_box(geode::BoundingBox3D)
 * @brief
 *
 * @param box
 */
void RemFracResSim::set_generation_box(geode::BoundingBox3D box) {
	generation_box_ = box;
}

//--------------------------
// set the model domain limits
//--------------------------

/**
 * @fn void set_model_fractures_bounding_box(const geode::StructuralModel&)
 * @brief
 *
 * @param model
 */
void RemFracResSim::set_model_fractures_bounding_box(
		const geode::StructuralModel &model) {

	set_generation_box(model.bounding_box());
	geode::Point3D min_box = generation_box_.min();
	geode::Point3D max_box = generation_box_.max();

	geode::Vector3D box_length(min_box, max_box);
	area_ = box_length.value(0) * box_length.value(1);
	volume_ = area_ * box_length.value(2);

	bbox_ = std::shared_ptr < fractures_intersect::Cboite_rect
			> (new fractures_intersect::Cboite_rect());

	fractures_intersect::Point3D pt = fractures_intersect::Point3D(
			min_box.value(0), min_box.value(1), min_box.value(2));
	bbox_->Set_dx(box_length.value(0));
	bbox_->Set_dy(box_length.value(1));
	bbox_->Set_dz(box_length.value(2));
	bbox_->Set_ind_materiaux(-1);
	bbox_->Set_type_poly(Zone_affinement);
	bbox_->Set_taille_element(10);
	bbox_->Set_nom("VOI");
	bbox_->Actualise_plans();
	model_fractures_.geo_cont_.Add_primitive(bbox_);
}
geode::index_t RemFracResSim::closest_element_index(geode::Point3D query, const geode::AABBTree3D& tree,const BoxAABBEvalDistance3D& disteval){



    geode::index_t box_id;
     geode::Point3D nearest_point;
     double distance;
   //std::tie( box_id, nearest_point, distance ) = tree.closest_element_box( query, disteval );
	return box_id;

}
geode::BoundingBox3D RemFracResSim::computePolyhedronBoundingBox(const geode::SolidMesh3D &mesh, geode::index_t polyhedronIndex) {
    geode::BoundingBox3D result;
    for (geode::index_t vertexIndex : mesh.polyhedron_vertices(polyhedronIndex)) {
        result.add_point(mesh.point(vertexIndex));
    }
    return result;
}
std::vector<geode::BoundingBox3D> RemFracResSim::create_bounding_box(const geode::StructuralModel &model){
	  geode::Logger::info("Building aabb tree from structural model ", model.id().string());

	    // Compute the cell center of each cells, and the (flat) bounding box around each cell center.
	    std::vector<geode::BoundingBox3D> boundingBoxes;
	    for (const geode::Block3D &block : model.blocks()) {
	        const geode::SolidMesh3D &mesh = block.mesh();
	        for (geode::index_t polyhedronIndex : geode::Range(mesh.nb_polyhedra())) {
	            cells_.emplace_back(CellId{ block.id(), polyhedronIndex });
	            boundingBoxes.emplace_back(computePolyhedronBoundingBox(mesh, polyhedronIndex));
	        }
	    }

	    // Build the tree with all the bounding boxes.
	    // Warning: the order of elements in `boundingBoxes` must correspond to the order of `cells_`.

		return boundingBoxes;

}


 void RemFracResSim::get_property_attribute_names_from_block_model(const geode::StructuralModel &model){
	 for (const geode::Block3D &block : model.blocks()) {
		 const geode::SolidMesh3D &mesh = block.mesh();
		 absl::FixedArray< absl::string_view > array_names = mesh.polyhedron_attribute_manager().attribute_names();
		 std::cout << "Block : " << block.name() << endl;
		 if(array_names.size() > 0){
			 for(auto& attr_name : array_names){
				 std::cout << attr_name << endl;
			 }
		 }else{
			 std::cout << " no properties " << endl;
		 }

	 }

}

void RemFracResSim::create_fault_zone_dfn(const FracResStage &stage,
		const FracResHypothesis &fracResHypothesis, int ireal,
		const geode::StructuralModel &model, bool estim_volume,
		const std::string &realisation_directory, FaultArraySet *setArrayZone,
		fractures_intersect::model &model_fractures_loc) {
	std::cout
			<< " --------------  FAULT ZONE CREATION IN PROGRESS ------------------"
			<< endl;
	setArrayZone->fault_array_set_index_ += 1000 * stage.stageIndex_
			+ 100 * fracResHypothesis.hypothesisIndex_;
	//model_fractures_.geo_cont_.tab_fractures_set_.size();
	setArrayZone->fault_array_set_index_ *= (ireal + 1);
	setArrayZone->seed_ = seed_;
	setArrayZone->volume_ = volume_;
	setArrayZone->fault_array_geo_set_ = std::shared_ptr
			< fractures_intersect::fault_arrays_set_geo
			> (new fractures_intersect::fault_arrays_set_geo(
					setArrayZone->fault_array_set_index_, setArrayZone->seed_));
	setArrayZone->area_ = area_;
	setArrayZone->set_generation_box(generation_box_);
	setArrayZone->fault_zone_generation(model, tree_, disteval_, cells_,
			model_fractures_loc.geo_cont_, estim_volume);
	setArrayZone->directory_path_ = realisation_directory;
}

void RemFracResSim::create_corridors_dfn(const FracResStage &stage,
		const FracResHypothesis &fracResHypothesis, int ireal,
		const std::string &realisation_directory,
		const geode::StructuralModel &model, bool estim_volume,
		FaultArraySet *setArrayCorridor,
		fractures_intersect::model &model_fractures_loc) {
	std::cout
			<< " --------------  FAULT CORRIDOR CREATION IN PROGRESS ------------------"
			<< endl;
	setArrayCorridor->fault_array_set_index_ += 1000 * stage.stageIndex_
			+ 100 * fracResHypothesis.hypothesisIndex_;
	setArrayCorridor->fault_array_set_index_ *= (ireal + 1);
	//model_fractures_.geo_cont_.tab_fractures_set_.size();
	setArrayCorridor->seed_ = seed_;
	setArrayCorridor->volume_ = volume_;
	setArrayCorridor->fault_array_geo_set_ = std::shared_ptr
			< fractures_intersect::fault_arrays_set_geo
			> (new fractures_intersect::fault_arrays_set_geo(
					setArrayCorridor->fault_array_set_index_,
					setArrayCorridor->seed_));
	setArrayCorridor->area_ = area_;
	setArrayCorridor->set_generation_box(generation_box_);
	setArrayCorridor->directory_path_ = realisation_directory;
	setArrayCorridor->fault_corridor_generation(model, tree_, disteval_, cells_,
			model_fractures_loc.geo_cont_, estim_volume);
}

int RemFracResSim::create_bed_parallel_slip_dfn( const geode::StructuralModel &model,fractures_intersect::model& model_fractures_loc,FracResBeddingProperties *bed_parallel, std::string path,std::string output_prefixe_file_name) {
	int nb_failles =0;
	bool BedInterfacePropertyIsSameAsBedParallelProperty = false;
	if (bed_interface_property_name_
			== bed_parallel->bed_parrallel_discrete_property_name_) {
		BedInterfacePropertyIsSameAsBedParallelProperty = true;
	}
	if (BedInterfacePropertyIsSameAsBedParallelProperty) {
		if (nb_total_bed_ > 0) {
			update_bed_interface(model, model_fractures_loc,seed_, bed_parallel);
		}
	} else {
		create_bed_parallel(model, seed_, bed_parallel);
		bed_parallel->index_fracture_set_geo_begin_ = model_fractures_loc.geo_cont_.tab_fractures_set_.size() - 1;
		for ( auto stress_geo_parallel : stress_geo_parallel_list_) {

			std::shared_ptr<fractures_intersect::fracture_set_geo> geo_parallel = stress_geo_parallel.second;
			model_fractures_loc.geo_cont_.add_fracture_set(geo_parallel);
			nb_failles++;
		}
	}
	bed_parallel->output_directory_path_ = path;
	bed_parallel->output_prefixe_name_ = output_prefixe_file_name;

	bed_parallel->nb_fracture_set_geo_ = stress_geo_parallel_list_.size();

/*
	int begin = bed_parallel->index_fracture_set_geo_begin_;
	int end = bed_parallel->index_fracture_set_geo_begin_+ bed_parallel->nb_fracture_set_geo_;

	for (int ifaille = begin; ifaille < end; ifaille++) {

		std::string methode_suffixe = absl::StrCat("_BEDDINGS_","PARALLEL_SLIP_");
		std::string surface_filename = absl::StrCat(output_prefixe_file_name, methode_suffixe);
		std::string surface_final = absl::StrCat(surface_filename, bed_parallel->name_);
		const auto output_file_native = absl::StrCat(path,surface_final, ".vtk");
		save_fracture_surface_mesh_from_model(ifaille,model_fractures_loc,output_file_native);
	}

        */


	return nb_failles;
}

void RemFracResSim::create_fault_array_dfn(const FracResStage &stage,
		const FracResHypothesis &fracResHypothesis, int ireal,
		const std::string &realisation_directory,
		const geode::StructuralModel &model, bool estim_volume,
		FaultArraySet *setArrayZone,
		fractures_intersect::model &model_fractures_loc) {
	std::cout
			<< " --------------  FAULT ARRAY CREATION IN PROGRESS ------------------"
			<< endl;
	setArrayZone->fault_array_set_index_ += 1000 * stage.stageIndex_
			+ 100 * fracResHypothesis.hypothesisIndex_;
	setArrayZone->fault_array_set_index_ *= (ireal + 1);
	//model_fractures_.geo_cont_.tab_fractures_set_.size();
	setArrayZone->seed_ = seed_;
	setArrayZone->volume_ = volume_;
	setArrayZone->fault_array_geo_set_ = std::shared_ptr
			< fractures_intersect::fault_arrays_set_geo
			> (new fractures_intersect::fault_arrays_set_geo(
					setArrayZone->fault_array_set_index_, setArrayZone->seed_));
	setArrayZone->area_ = area_;
	setArrayZone->set_generation_box(generation_box_);
	setArrayZone->directory_path_ = realisation_directory;
	setArrayZone->fault_array_generation(model, tree_, disteval_, cells_,
			model_fractures_loc.geo_cont_, estim_volume);
}

void RemFracResSim::create_single_fault_dfn(const FracResStage &stage,
		const FracResHypothesis &fracResHypothesis, int ireal,
		fractures_intersect::model &model_fractures_loc,
		const std::string &realisation_directory,
		const geode::StructuralModel &model, bool estim_volume,
		FractureSet *setFrac) {
	std::cout
			<< " --------------  SINGLE FAULT CREATION IN PROGRESS ------------------"
			<< endl;
	setFrac->fracture_set_index_ += 1000 * stage.stageIndex_
			+ 100 * fracResHypothesis.hypothesisIndex_;
	setFrac->fracture_set_index_ *= (ireal + 1);
	setFrac->seed_ = seed_;
	setFrac->set_volume(volume_);
	setFrac->model_fractures_.geo_cont_ = model_fractures_loc.geo_cont_;
	setFrac->set_area(area_);
	setFrac->set_generation_box(generation_box_);
	setFrac->directory_path_ = realisation_directory;
	std::shared_ptr<fault_array> test_array = NULL;
	setFrac->fractures_generation_from_shadow_zone(model, tree_, disteval_,
			cells_, test_array, estim_volume);
	if (setFrac->stress_zone_set_->tab_fractures_.size() > 0) {
		model_fractures_loc.geo_cont_.add_fracture_set(
				setFrac->stress_zone_set_);
		if (setFrac->is_shadow_zone_active_) {
			setFrac->type_intersection_ =
					fractures_intersect::FracResIntersectionTypeEnum::CROSSING;
			cout << " Fracture set " << setFrac->name_
					<< " is crossing with shadow zone" << endl;
		} else {
			if (setFrac->type_intersection_
					== fractures_intersect::FracResIntersectionTypeEnum::CROSSING) {
				cout << " Fracture set " << setFrac->name_
						<< " is crossing with no shadow zone" << endl;
			} else if (setFrac->type_intersection_
					== fractures_intersect::FracResIntersectionTypeEnum::ABUTTING) {
				cout << " Fracture set " << setFrac->name_
						<< " is abbuting with no shadow zone" << endl;
			} else {
				cout << " Fracture set " << setFrac->name_
						<< " is randomly crossing and abbuting with no shadow zone"
						<< endl;
			}
		}
		setFrac->fractures_activable_ = setFrac->fractures_number_;
		setFrac->fractures_inserted_ = 0;
		setFrac->model_geo_cont_index_= model_fractures_loc.geo_cont_.tab_fractures_set_.size()-1;
	}
}

void RemFracResSim::time_creation_message(const std::string &nameDfn,std::chrono::time_point < std::chrono::system_clock>& total_StartTime) {
	std::chrono::time_point < std::chrono::system_clock> total_EndTime = std::chrono::system_clock::now();
	double write_miliseconde = std::chrono::duration_cast< std::chrono::milliseconds> (total_EndTime - total_StartTime).count();
	double write_second = write_miliseconde / 1000;
	cout << " Temps total de simulation du dfn " << nameDfn << endl;
	cout << " nb miliseconds " << write_miliseconde << " nb second "<< write_second << '\n';
	cout << endl;
}

void RemFracResSim::create_bed_interface_slip_dfn(const geode::StructuralModel &model,const FracResHypothesis &fracResHypothesis,fractures_intersect::model &model_fractures_loc,
		std::string path,std::string output_prefixe_file_name, int index) {
	bed_interface_property_name_ = fracResHypothesis.beddingPropertyDiscreteName_;
	FracResBeddingAperture *bed_interface = fracResHypothesis.beddingList_[index];
	if(index > 0){
		bed_interface->index_begin_ = fracResHypothesis.beddingList_[index -1]->index_end_;

	}

	create_bed_surface(model, fracResHypothesis.beddingList_,fracResHypothesis.beddingPropertyDiscreteName_, index);
	nb_total_bed_ = 0;


	for (int j = 0;j < bed_interface->discrete_property_selected_values_.size();j++) {
		std::shared_ptr<fractures_intersect::fracture_set_geo> stress_geo =  stress_geo_list_[bed_interface->discrete_property_selected_values_[j]];
		cout << "nb fractures for ( " << bed_interface->discrete_property_selected_values_[j].first<< " , "<< bed_interface->discrete_property_selected_values_[j].second << " )= "<< stress_geo->tab_fractures_.size() << endl;
		nb_total_bed_ += stress_geo->tab_fractures_.size();
		model_fractures_loc.geo_cont_.add_fracture_set(stress_geo);
		bed_interface->index_end_++;

	}
}

void RemFracResSim::message_bedding_interface(const FracResStage &stage) {
	if (stage.beddingStage_) {
		std::cout
				<< " -------------- bedding : Bed_interface Hypothesis -------------------"
				<< endl;
	} else {
		std::cout
				<< " ...............Bedding : No bed interface ---------------------------"
				<< endl;
	}
}





int RemFracResSim::message_hypothesis_dfn_status(FracResHypothesis& fracResHypothesis){
	int nb_total_dfn = 0;

	int nb_bedding_interface = fracResHypothesis.beddingList_.size();
	int nb_bedding_parallel =fracResHypothesis.bed_parallel_List_.size();
	fracResHypothesis.update_active_dfn_number();


	int nb_single_fault = fracResHypothesis.nb_single_fault_active_;
	int nb_array_fault = fracResHypothesis.nb_fault_array_active_;
	int nb_corridor_fault = fracResHypothesis.nb_cluster_corridors_active_;
	int nb_fault_zone = fracResHypothesis.nb_cluster_fault_zone_active_;
	int nb_total = 0;
	if (nb_bedding_parallel > 0) {
		std::cout << "number of bed parallel : " << nb_bedding_parallel<< endl;
		nb_total += nb_bedding_parallel;
	}
	if (nb_bedding_interface > 0) {
		std::cout << "number of bed interface : "<< nb_bedding_interface << endl;
		nb_total += nb_bedding_interface;
	}
	if (nb_single_fault > 0) {
		std::cout << "number of single fault : " << nb_single_fault<< endl;
		nb_total += nb_single_fault;
	}
	if (nb_array_fault > 0) {
		std::cout << "number of array fault : " << nb_array_fault<< endl;
		nb_total += nb_array_fault;
	}

	if (nb_corridor_fault > 0) {
		std::cout << "number of clusters corridors  : "<< nb_corridor_fault << endl;
		nb_total += nb_corridor_fault;
	}
	if (nb_fault_zone > 0) {
		std::cout << "number of clusters fault zone  : "<< nb_fault_zone << endl;
		nb_total += nb_fault_zone;
	}
	std::cout << endl;
	std::cout << "number of dfn  : " << nb_total << " " << endl;
	nb_total_dfn += nb_total;
	return nb_total_dfn;
}

void RemFracResSim::create_bed_interface_from_hypothesis(const FracResHypothesis &fracResHypothesis, int indice_bed_interface,const geode::StructuralModel &model,
		fractures_intersect::model &model_fractures_loc,std::string nameDfn, std::string path,std::string output_prefixe_file_name) {
	FracResBeddingAperture *bed_interface =fracResHypothesis.beddingList_[indice_bed_interface];


	bed_interface->output_prefixe_name_=output_prefixe_file_name;
	bed_interface->output_directory_path_ = path;




	std::chrono::time_point < std::chrono::system_clock > total_StartTime =std::chrono::system_clock::now();
	create_bed_interface_slip_dfn(model, fracResHypothesis,model_fractures_loc,path,output_prefixe_file_name,indice_bed_interface);
	if(indice_bed_interface ==0){
		fracResBeddingzero_->index_begin_=fracResHypothesis.beddingList_[indice_bed_interface]->discrete_property_selected_values_.size();
		fracResBeddingzero_->index_end_=fracResHypothesis.beddingList_[indice_bed_interface]->discrete_property_selected_values_.size();
		fracResBeddingzero_->output_prefixe_name_=output_prefixe_file_name;
		fracResBeddingzero_->output_directory_path_ = path;
		for (int j = 0;j < fracResBeddingzero_->discrete_property_selected_values_.size();j++) {
				std::shared_ptr<fractures_intersect::fracture_set_geo> stress_geo =  stress_geo_list_[fracResBeddingzero_->discrete_property_selected_values_[j]];
				cout << "nb fractures for ( " << fracResBeddingzero_->discrete_property_selected_values_[j].first<< " , "<< fracResBeddingzero_->discrete_property_selected_values_[j].second << " )= "<< stress_geo->tab_fractures_.size() << endl;
				nb_total_bed_ += stress_geo->tab_fractures_.size();
				model_fractures_loc.geo_cont_.add_fracture_set(stress_geo);
				fracResBeddingzero_->index_end_++;

			}

		bed_interface->index_end_ += fracResBeddingzero_->discrete_property_selected_values_.size();
	}
	nameDfn = "";
	time_creation_message(nameDfn, total_StartTime);


}

void RemFracResSim::save_bed_interface_from_hypothesis(const FracResHypothesis &fracResHypothesis, int indice_bed_interface,fractures_intersect::model &model_fractures_loc) {
	FracResBeddingAperture *bed_interface =fracResHypothesis.beddingList_[indice_bed_interface];


	if (bed_interface->is_active_) {
		int ifaille=0;


		for (int j = bed_interface->index_begin_;j < bed_interface->index_begin_ + bed_interface->discrete_property_selected_values_.size();j++) {



			std::vector<fractures_intersect::Point3D> tab_point;
			std::vector<std::vector<int>> tab_triangle;
			std::vector<fractures_intersect::property<float>> tab_triangle_face_property;

			model_fractures_loc.geo_cont_.get_fracture_set_surface_mesh(j, tab_point, tab_triangle,tab_triangle_face_property);

			std::string methode_suffixe = absl::StrCat("_BEDDINGS_","INTERFACE_SLIP_");

			std::string surface_filename = absl::StrCat(bed_interface->output_prefixe_name_, methode_suffixe);
			std::pair<int, int> val = bed_interface->discrete_property_selected_values_[j-bed_interface->index_begin_];
			std::string surface_final = absl::StrCat(surface_filename, bed_name_vect_[val]);

			std::string discrete_pair = std::to_string(val.first)+ "_" + std::to_string(val.second);
			surface_final = absl::StrCat(surface_filename,discrete_pair);
			const auto output_file_native = absl::StrCat(bed_interface->output_directory_path_,surface_final, ".vtk");
			save_meshed_vtk_plan(output_file_native, tab_point,tab_triangle, tab_triangle_face_property);
			if (tab_point.size() <= 0) {ifaille++;}
		}

		if(indice_bed_interface == 0){
			for (int j = fracResBeddingzero_->index_begin_;j < fracResBeddingzero_->index_end_;j++) {



						std::vector<fractures_intersect::Point3D> tab_point;
						std::vector<std::vector<int>> tab_triangle;
						std::vector<fractures_intersect::property<float>> tab_triangle_face_property;

						model_fractures_loc.geo_cont_.get_fracture_set_surface_mesh(j, tab_point, tab_triangle,tab_triangle_face_property);

						std::string methode_suffixe = absl::StrCat("_BEDDINGS_","INTERFACE_SLIP_");

						std::string surface_filename = absl::StrCat(fracResBeddingzero_->output_prefixe_name_, methode_suffixe);
						std::pair<int, int> val = fracResBeddingzero_->discrete_property_selected_values_[j-fracResBeddingzero_->index_begin_];
						std::string surface_final = absl::StrCat(surface_filename, bed_name_vect_[val]);

						std::string discrete_pair = std::to_string(val.first)+ "_" + std::to_string(val.second);
						surface_final = absl::StrCat(surface_filename,discrete_pair);
						const auto output_file_native = absl::StrCat(fracResBeddingzero_->output_directory_path_,surface_final, ".vtk");
						save_meshed_vtk_plan(output_file_native, tab_point,tab_triangle, tab_triangle_face_property);
						if (tab_point.size() <= 0) {ifaille++;}
					}
		}

	}
}

int RemFracResSim::get_evaluate_fracture_single_faults_from_hypothesis(const FracResHypothesis &fracResHypothesis, int indice_single_fault,
		const FracResStage &stage, const geode::StructuralModel &model) {
	 int nb_fracture =0;
	 FractureSet* setFrac= fracResHypothesis.fractureList_ [indice_single_fault];
	// for( FractureSet* setFrac : fracResHypothesis.fractureList_ ){
		if (setFrac->is_active_) {
			//cout << " the single fault " << setFrac->name_ << " is active " << endl;
			evaluate_single_fault_dfn(stage, fracResHypothesis, model, false,setFrac);

			nb_fracture = setFrac->evaluate_fractures_number_;

		}
	// }

	return nb_fracture;
}

void RemFracResSim::evaluate_single_fault_dfn(const FracResStage &stage,
		const FracResHypothesis &fracResHypothesis,const geode::StructuralModel &model, bool estim_volume,FractureSet *setFrac) {
	//std::cout
	//		<< " --------------  SINGLE FAULT EVALUATION FRACTURES NUMBER IN PROGRESS ------------------"
	//		<< endl;
	setFrac->seed_ = seed_;
	setFrac->set_volume(volume_);
	setFrac->set_area(area_);
	setFrac->set_generation_box(generation_box_);
	std::shared_ptr<fault_array> test_array = NULL;
	setFrac->fractures_evaluate_generation(model, tree_, disteval_,cells_, test_array, estim_volume);

}

int RemFracResSim::create_single_faults_from_hypothesis(const FracResHypothesis &fracResHypothesis, int indice_single_fault,
		const FracResStage &stage, int ireal,fractures_intersect::model &model_fractures_loc,
		const std::string &realisation_directory,const geode::StructuralModel &model, bool estim_volume,std::string nameDfn,std::string output_prefixe_name) {
	 int nb_success =0;
	FractureSet *setFrac = fracResHypothesis.fractureList_[indice_single_fault];
	setFrac->output_prefixe_name_=output_prefixe_name;
	if (setFrac->is_active_) {
		cout << " the single fault " << setFrac->name_ << " is active " << endl;
		std::chrono::time_point < std::chrono::system_clock > total_StartTime =std::chrono::system_clock::now();
		create_single_fault_dfn(stage, fracResHypothesis, ireal,model_fractures_loc, realisation_directory, model, estim_volume,setFrac);
		nameDfn = setFrac->name_;
		//model_fractures_loc.geo_cont_.intersect_all_fracture_set(10);
		time_creation_message(setFrac->name_, total_StartTime);
		if (setFrac->success_) {
			nb_success++;
		}

		//saveFractureSet(model_fractures_loc,output_prefixe_name, setFrac);
	}
	return nb_success;	return nb_success;

}

bool RemFracResSim::save_output_single_faults_from_hypothesis(fractures_intersect::model &model_fractures_loc,const FracResHypothesis &fracResHypothesis, int indice_single_fault){
	FractureSet *setFrac = fracResHypothesis.fractureList_[indice_single_fault];
	if (setFrac->is_active_) {

		std::string methode_suffixe = absl::StrCat("_", fracResHypothesis.stage_name_,"_");
		std::string surface_filename = absl::StrCat(setFrac->output_prefixe_name_,
				methode_suffixe);
	saveFractureSet(model_fractures_loc,surface_filename, setFrac);
	return true;
	}
	return false;
}

bool RemFracResSim::save_output_array_faults_from_hypothesis(std::string methode_suffixe,const FracResHypothesis &fracResHypothesis, int indice_fault_array, int type_array){

	FaultArraySet *setArrayZone = NULL;
	if(type_array == 0){
		setArrayZone =fracResHypothesis.faultArrayList_[indice_fault_array];
	}else if(type_array == 1){
		setArrayZone =fracResHypothesis.faultCorridorList_[indice_fault_array];
	}else if(type_array == 2){
		setArrayZone =fracResHypothesis.faultZoneList_[indice_fault_array];
	}
	if (setArrayZone != NULL && setArrayZone->is_active_) {

		std::string fileNameNew = absl::StrCat(setArrayZone->output_prefixe_name_,"_",fracResHypothesis.stage_name_);

		methode_suffixe = absl::StrCat("_", enumCategoryToString(setArrayZone->type_category_),
							"_");


	saveArrayFault(methode_suffixe,fileNameNew, setArrayZone);
	return true;
	}
	return false;
}


int RemFracResSim::create_fault_array_from_hypothesis(const FracResHypothesis &fracResHypothesis, int indice_fault_array,
		const FracResStage &stage, int ireal,const std::string &realisation_directory,
		const geode::StructuralModel &model, bool estim_volume,std::string nameDfn,std::string output_prefixe_name,fractures_intersect::model &model_fractures_loc) {
	 int nb_success=0;
	FaultArraySet *setArrayZone =fracResHypothesis.faultArrayList_[indice_fault_array];
	setArrayZone->output_prefixe_name_=output_prefixe_name;
	if (setArrayZone->is_active_) {
		cout << " the fault array " << setArrayZone->name_ << " is active "<< endl;
		std::chrono::time_point < std::chrono::system_clock > total_StartTime =std::chrono::system_clock::now();
		create_fault_array_dfn(stage, fracResHypothesis, ireal,realisation_directory, model, estim_volume, setArrayZone,model_fractures_loc);
		nameDfn = setArrayZone->name_;
		//model_fractures_loc.geo_cont_.intersect_all_fracture_set(10);
		time_creation_message(setArrayZone->name_, total_StartTime);
		if (setArrayZone->success_) {
			nb_success++;
		}

		std::string methode_suffixe = "_FAULTS_ARRAY_";
		//saveArrayFault(methode_suffixe,output_prefixe_name, setArrayZone);
	}
	return nb_success;
}


int RemFracResSim::get_evaluate_fracture_fault_array_from_hypothesis(const FracResHypothesis &fracResHypothesis, int indice_fault_array,
		const FracResStage &stage, const geode::StructuralModel &model) {
	 int nb_fracture =0;
	FaultArraySet *setArrayZone =fracResHypothesis.faultArrayList_[indice_fault_array];
	// for( FaultArraySet* setArrayZone : fracResHypothesis.faultArrayList_ ){
		if (setArrayZone->is_active_) {
			//cout << " the fault array  " << setArrayZone->name_ << " is active " << endl;
			nb_fracture = evaluate_fault_array_dfn(stage, fracResHypothesis, model, false,setArrayZone);
	       // break;
		}
	// }
	return nb_fracture;
}

int RemFracResSim::evaluate_fault_array_dfn(const FracResStage &stage,
		const FracResHypothesis &fracResHypothesis,const geode::StructuralModel &model, bool estim_volume,FaultArraySet *setArray) {
	//std::cout
	//		<< " --------------  FAULT ARRAY EVALUATION FRACTURES NUMBER IN PROGRESS ------------------"
	//		<< endl;
	setArray->seed_ = seed_;
	setArray->set_volume(volume_);
	setArray->set_area(area_);
	setArray->set_generation_box(generation_box_);
	int valeur_fracture = setArray->fault_array_fractures_evaluation(model, tree_, disteval_,cells_,estim_volume);
    return valeur_fracture;
}
int RemFracResSim::get_evaluate_fracture_cluster_corridor_from_hypothesis(const FracResHypothesis &fracResHypothesis, int indice_fault_array,
		const FracResStage &stage, const geode::StructuralModel &model) {
	 int nb_fracture =0;
	 FaultArraySet *setArrayCorridor =fracResHypothesis.faultCorridorList_[indice_fault_array];
    //for( FaultArraySet* setArrayCorridor : fracResHypothesis.faultCorridorList_){
		if (setArrayCorridor->is_active_) {
			//cout << " the cluster fault corridor  " << setArrayCorridor->name_ << " is active " << endl;
			nb_fracture = evaluate_fault_array_dfn(stage, fracResHypothesis, model, false,setArrayCorridor);
			//break;
		}
    //}
	return nb_fracture;
}

int RemFracResSim::get_evaluate_fracture_cluster_fault_zone_from_hypothesis(const FracResHypothesis &fracResHypothesis, int indice_fault_array,
		const FracResStage &stage, const geode::StructuralModel &model) {
	 int nb_fracture =0;
	 FaultArraySet *setArrayFaultZone =fracResHypothesis.faultZoneList_[indice_fault_array];
	 //for( FaultArraySet* setArrayFaultZone : fracResHypothesis.faultZoneList_){
		if (setArrayFaultZone->is_active_) {
			//cout << " the cluster fault zone  " << setArrayFaultZone->name_ << " is active " << endl;
			nb_fracture = evaluate_fault_array_dfn(stage, fracResHypothesis, model, false,setArrayFaultZone);
			//break;
		}
	// }
	return nb_fracture;
}





int RemFracResSim::create_bed_parallel_slip_from_hypothesis(const FracResHypothesis &fracResHypothesis, int indice_bed_parallel,const geode::StructuralModel &model,
		fractures_intersect::model &model_fractures_loc,std::string nameDfn,std::string path,std::string output_prefixe_file_name) {
	FracResBeddingProperties *bed_parallel =fracResHypothesis.bed_parallel_List_[indice_bed_parallel];
	bed_parallel->output_prefixe_name_=output_prefixe_file_name;
	bed_parallel->output_directory_path_ = path;

	std::chrono::time_point < std::chrono::system_clock > total_StartTime =std::chrono::system_clock::now();
	int nb_faille_loc = create_bed_parallel_slip_dfn(model, model_fractures_loc,bed_parallel,path,output_prefixe_file_name);
	nameDfn = bed_parallel->name_;
	time_creation_message(nameDfn, total_StartTime);
	return nb_faille_loc;
}

void RemFracResSim::save_bed_parallel_slip_from_hypothesis(const FracResHypothesis &fracResHypothesis, int indice_bed_parallel,
		fractures_intersect::model &model_fractures_loc) {
	FracResBeddingProperties *bed_parallel =fracResHypothesis.bed_parallel_List_[indice_bed_parallel];
	bed_parallel->nb_fracture_set_geo_ = stress_geo_parallel_list_.size();
	bed_parallel->index_fracture_set_geo_begin_ = model_fractures_loc.geo_cont_.tab_fractures_set_.size() - 1;

	int begin = bed_parallel->index_fracture_set_geo_begin_;
	int end = bed_parallel->index_fracture_set_geo_begin_+ bed_parallel->nb_fracture_set_geo_;


	for (int ifaille = begin; ifaille < end; ifaille++) {

		std::string methode_suffixe = absl::StrCat("_BEDDINGS_","PARALLEL_SLIP_");
		std::string surface_filename = absl::StrCat(bed_parallel->output_prefixe_name_, methode_suffixe);
		std::string surface_final = absl::StrCat(surface_filename, bed_parallel->name_);
		const auto output_file_native = absl::StrCat(bed_parallel->output_directory_path_,surface_final, ".vtk");
		save_fracture_surface_mesh_from_model(ifaille,model_fractures_loc,output_file_native);
	}

}



int RemFracResSim::create_clusters_corridors_from_hypothesis(const FracResHypothesis &fracResHypothesis, int indice_fault_corridor,const FracResStage &stage, int ireal,
		const std::string &realisation_directory,const geode::StructuralModel &model, bool estim_volume,std::string nameDfn,std::string output_prefixe_name,fractures_intersect::model &model_fractures_loc) {
	int nb_success=0;
	FaultArraySet *setArrayCorridor =fracResHypothesis.faultCorridorList_[indice_fault_corridor];
	setArrayCorridor->output_prefixe_name_=output_prefixe_name;
	if (setArrayCorridor->is_active_) {
		cout << " the cluster fault corridor " << setArrayCorridor->name_<< " is active " << endl;
		std::chrono::time_point < std::chrono::system_clock > total_StartTime =std::chrono::system_clock::now();
		create_corridors_dfn(stage, fracResHypothesis, ireal,realisation_directory, model, estim_volume, setArrayCorridor,model_fractures_loc);
		nameDfn = setArrayCorridor->name_;
		//model_fractures_loc.geo_cont_.intersect_all_fracture_set(10);
		time_creation_message(setArrayCorridor->name_, total_StartTime);
		if (setArrayCorridor->success_) {
			nb_success++;
		}

		//std::string methode_suffixe ="_CLUSTERS_FAULTS_CORRIDORS_";
		//saveArrayFault(methode_suffixe,output_prefixe_name, setArrayCorridor);

	}
	return nb_success;
}

int RemFracResSim::create_clusters_fault_zone_from_hypothesis(const FracResHypothesis &fracResHypothesis, int indice_fault_zone,const FracResStage &stage, int ireal,
		const geode::StructuralModel &model, bool estim_volume,const std::string &realisation_directory, std::string nameDfn,std::string output_prefixe_name,fractures_intersect::model &model_fractures_loc) {
	int nb_success=0;
	FaultArraySet *setArrayZone =fracResHypothesis.faultZoneList_[indice_fault_zone];
	setArrayZone->output_prefixe_name_=output_prefixe_name;
	if (setArrayZone->is_active_) {
		cout << " the cluster fault zone " << setArrayZone->name_<< " is active " << endl;
		std::chrono::time_point < std::chrono::system_clock > total_StartTime =std::chrono::system_clock::now();
		create_fault_zone_dfn(stage, fracResHypothesis, ireal, model,estim_volume, realisation_directory, setArrayZone,model_fractures_loc);
		nameDfn = setArrayZone->name_;
		//model_fractures_loc.geo_cont_.intersect_all_fracture_set(10);
		time_creation_message(nameDfn, total_StartTime);
		if (setArrayZone->success_) {
			nb_success++;
		}

		//std::string methode_suffixe = "_CLUSTERS_FAULTS_ZONES_";
		//saveArrayFault(methode_suffixe,output_prefixe_name, setArrayZone);

	}
	return nb_success;
}

void RemFracResSim::run_simulation_from_stage_and_hypothesesis(const FracResHypothesis &fracResHypothesis,const geode::StructuralModel &model, const std::string &nameDfn,
		const FracResStage &stage, int ireal,const std::string &realisation_directory, bool estim_volume,fractures_intersect::model &model_fractures_loc,FracResScenario &scenario) {

	int indice_bed_interface = 0;
	int indice_bed_parallel = 0;
	int indice_single_fault = 0;
	int indice_fault_array = 0;
	int indice_fault_corridor = 0;
	int indice_fault_zone = 0;
	std::string output_file_name_prefixe = scenario.scenario_file_name_ + "_Real0" + std::to_string(ireal+1);

	stress_geo_list_.clear();
	bed_name_vect_.clear();
	active_discrete_pair_array_.clear();


	for (FracResBeddingAperture *fracTest :fracResHypothesis.beddingList_) {

		for (int j = 0; j < fracTest->discrete_property_selected_values_.size();j++) {
			active_discrete_pair_array_.insert({fracTest->discrete_property_selected_values_[j], true});
		}
	}
	fracResBeddingzero_ = std::shared_ptr< remfracres::FracResBeddingAperture > (new FracResBeddingAperture());


	for (int type_dfn : fracResHypothesis.type_dfn_) {
		int nb_success = 0;
		switch (type_dfn) {
		case 0: {
			create_bed_interface_from_hypothesis(fracResHypothesis,indice_bed_interface, model, model_fractures_loc, nameDfn,realisation_directory,output_file_name_prefixe);
			scenario.nb_failles_ =model_fractures_loc.geo_cont_.tab_fractures_set_.size();

			if (scenario.nb_failles_ > 0) indice_bed_interface++;
			break;
		}
		case 1: {
			nb_success = create_single_faults_from_hypothesis(fracResHypothesis,indice_single_fault, stage, ireal, model_fractures_loc,
					realisation_directory, model, estim_volume, nameDfn,output_file_name_prefixe);
			if (nb_success >= 0) indice_single_fault++;
			scenario.nb_failles_ += nb_success;
			break;
		}
		case 2: {
			nb_success = create_fault_array_from_hypothesis(fracResHypothesis,indice_fault_array, stage, ireal, realisation_directory,
					model, estim_volume, nameDfn, output_file_name_prefixe,model_fractures_loc);
			if (nb_success >= 0) indice_fault_array++;
			scenario.nb_failles_ += nb_success;
			break;
		}
		case 3: {
			int nb_faille_loc = create_bed_parallel_slip_from_hypothesis(fracResHypothesis, indice_bed_parallel, model,model_fractures_loc, nameDfn,realisation_directory,output_file_name_prefixe);
			scenario.nb_failles_ += nb_faille_loc;
			if (nb_faille_loc >= 0) indice_bed_parallel++;
			break;
		}
		case 4: {
			nb_success = create_clusters_corridors_from_hypothesis(fracResHypothesis, indice_fault_corridor, stage, ireal,realisation_directory,
					model, estim_volume, nameDfn,output_file_name_prefixe,model_fractures_loc);
			if (nb_success >= 0) indice_fault_corridor++;
			scenario.nb_failles_ += nb_success;
			break;
		}
		case 5: {
			nb_success = create_clusters_fault_zone_from_hypothesis(fracResHypothesis, indice_fault_zone, stage, ireal, model,
					estim_volume, realisation_directory, nameDfn,output_file_name_prefixe,model_fractures_loc);
			if (nb_success >= 0) indice_fault_zone++;
			scenario.nb_failles_ += nb_success;
			break;
		}
		default: {
		}
		} // switch
	}

	// loop type_dfn}
}




void RemFracResSim::save_output_from_stage_and_hypothesesis(const FracResHypothesis &fracResHypothesis,const FracResStage &stage, fractures_intersect::model &model_fractures_loc,FracResScenario &scenario) {

	int indice_bed_interface = 0;
	int indice_bed_parallel = 0;
	int indice_single_fault = 0;
	int indice_fault_array = 0;
	int indice_fault_corridor = 0;
	int indice_fault_zone = 0;


	for (int type_dfn : fracResHypothesis.type_dfn_) {
		int nb_success = 0;
		switch (type_dfn) {
		case 0: {

			save_bed_interface_from_hypothesis(fracResHypothesis,indice_bed_interface,model_fractures_loc);
			indice_bed_interface++;
			break;
		}
		case 1: {
			save_output_single_faults_from_hypothesis(model_fractures_loc,fracResHypothesis,indice_single_fault);
			indice_single_fault++;
			break;
		}
		case 2: {



			std::string methode_suffixe = "_FAULTS_ARRAY_";
			save_output_array_faults_from_hypothesis(methode_suffixe,fracResHypothesis,indice_fault_array,0);
			indice_fault_array++;

			break;
		}
		case 3: {
			save_bed_parallel_slip_from_hypothesis(fracResHypothesis, indice_bed_parallel,model_fractures_loc);
			indice_bed_parallel++;
			break;
		}
		case 4: {
			std::string methode_suffixe = "_CLUSTERS_FAULTS_CORRIDORS_";
			save_output_array_faults_from_hypothesis(methode_suffixe,fracResHypothesis,indice_fault_corridor,1);
			indice_fault_corridor++;
			break;
		}
		case 5: {
			std::string methode_suffixe = "_CLUSTERS_FAULTS_ZONES_";
			save_output_array_faults_from_hypothesis(methode_suffixe,fracResHypothesis,indice_fault_zone,2);
			indice_fault_zone++;
			break;
		}
		default: {
		}
		} // switch
	}
	// loop type_dfn}
}




std::map<std::pair<int,int>,std::string> RemFracResSim::fractures_number_evaluation_for_scenario(const geode::StructuralModel &model){
	int nb_scenario_ = fracResScenarioList_.size();
	set_model_fractures_bounding_box(model);
	nb_failles_ = 0;
	std::map<std::pair<int,int>,std::string> mapElementNumber;

	for (FracResScenario& scenario : fracResScenarioList_) {
		std::pair<int,int> evaluate_number_pair;
		evaluate_number_pair.first = scenario.scenarioIndex_;
		int nb_totalfracture_scenario = 0;
		for( std::pair<int,int> stagePair : scenario.scenario_stage_hypothesis_index_vector_){
			FracResStage& stage = 	fracResStageList_[stagePair.first];
			FracResHypothesis& fracResHypothesis = stage.hypotheisisList_[stagePair.second-1];
			int nb_fractures = getEvaluateFracturesNumberPerStageHypothesis(model,stage,fracResHypothesis);
			fracResHypothesis.evaluate_fractures_number_total_ = nb_fractures;
			stage.evaluate_fractures_number_total_per_hypothesis_.push_back(nb_fractures);
			nb_totalfracture_scenario += nb_fractures;
		}
		evaluate_number_pair.second = nb_totalfracture_scenario;
		std::string message = " Evaluate fracture number of scenario " + scenario.name_ + " ";
		mapElementNumber.insert({evaluate_number_pair, message});
	}
	return mapElementNumber;
}


void RemFracResSim::affiche_fracture_number_evaluation(std::map<std::pair<int,int>,std::string>&  mapNumber){
	for (const auto& elem : mapNumber) {

		Logger::info(elem.second +  std::to_string(elem.first.second));

	}

}
//-------------------------------------------
// generation of the fractures
//-------------------------------------------
/**
 * @fn void fractures_set_generation_for_scenario(const geode::StructuralModel&, bool)
 * @brief
 *
 * @param model
 * @param estim_volume
 */
void RemFracResSim::fractures_set_generation_for_scenario(const geode::StructuralModel &model,std::string nameFile, bool estim_volume) {

	int nb_scenario_ = fracResScenarioList_.size();
	set_model_fractures_bounding_box(model);
	nb_failles_ = 0;
    std::string nameDfn;

	for (FracResScenario& scenario : fracResScenarioList_) {
		int nb_total_dfn = 0;
		std::cout << " -------------- Scernario N°" << scenario.scenarioIndex_ << " named " << scenario.name_ << "------------------" << endl;
		int nb_stage = scenario.scenario_stage_hypothesis_index_vector_.size();
		std::vector<bool>success_list;
		success_list.clear();
		if(!scenario.modeling_) continue;
		std::string scenario_directory = outpout_directory_path_ + scenario.name_ + "/";
		create_local_directory_output(scenario_directory);
		scenario.scenario_file_name_ = nameFile + "_"+ scenario.name_;
		for( int ireal =0; ireal < scenario.numberOfRealisation_; ireal ++){

			std::string realisation_directory = outpout_directory_path_ + scenario.name_ + "/Real0" + std::to_string(ireal+1) + "/";
			create_local_directory_output(realisation_directory);
			fractures_intersect::model model_fractures_loc;
			model_fractures_loc.geo_cont_.Add_primitive(bbox_);


			std::cout << "Realisation : " << ireal+1 << endl;
			double simul_second_tot=0;
			for( std::pair<int,int> stagePair : scenario.scenario_stage_hypothesis_index_vector_){

				FracResStage& stage = 	fracResStageList_[stagePair.first];
				FracResHypothesis& fracResHypothesis = stage.hypotheisisList_[stagePair.second-1];
				std::cout<< " stage : " << stage.name_ << " and hypothesis : " << fracResHypothesis.name_ << std::endl;
				message_bedding_interface(stage);
				nb_total_dfn += message_hypothesis_dfn_status(fracResHypothesis);

				std::chrono::time_point < std::chrono::system_clock > total_StartTime;
				std::chrono::time_point < std::chrono::system_clock > total_EndTime;
			    total_StartTime = std::chrono::system_clock::now();


				run_simulation_from_stage_and_hypothesesis(fracResHypothesis, model, nameDfn, stage, ireal, realisation_directory, estim_volume, model_fractures_loc, scenario);


				total_EndTime = std::chrono::system_clock::now();
				double write_miliseconde = std::chrono::duration_cast< std::chrono::milliseconds> (total_EndTime - total_StartTime).count();
				double write_second = write_miliseconde / 1000;
				cout << " Temps total de simulation et d'ecriture du résultat de l'hypothese" << endl;
				cout << " nb miliseconds " << write_miliseconde << " nb second "<< write_second << '\n';
				cout << endl;
				simul_second_tot +=write_miliseconde;
				if (scenario.nb_failles_ == (nb_total_dfn)){
					success_list.push_back(true);
				}else{
					success_list.push_back(false);
				}
			} //loop pair a effectuer

			model_fractures_loc.geo_cont_.intersect_all_fracture_set(10);

			for( std::pair<int,int> stagePair : scenario.scenario_stage_hypothesis_index_vector_){

				FracResStage& stage = 	fracResStageList_[stagePair.first];
				FracResHypothesis& fracResHypothesis = stage.hypotheisisList_[stagePair.second-1];
				save_output_from_stage_and_hypothesesis(fracResHypothesis,stage, model_fractures_loc, scenario);

			}
			scenario.total_simulation_time_list_.push_back(simul_second_tot);
			scenario.success_list_.insert({ireal,success_list});
			scenario.model_fractures_list_.push_back(model_fractures_loc);
		} // loop realisation

	} // loop scenario



}

//-------------------------------------------
// generation of the fractures
//-------------------------------------------
/**
 * @fn void stage_hypothesis_fractures_set_generation(const geode::StructuralModel&, bool)
 * @brief
 *
 * @param model
 * @param estim_volume
 */
void RemFracResSim::stage_hypothesis_fractures_set_generation(
		const geode::StructuralModel &model, bool estim_volume) {

	int nb_stage_ = fracResStageList_.size();
	set_model_fractures_bounding_box(model);
	nb_failles_ = 0;
	int nb_success = 0;
	int nb_total_dfn = 0;
	for (FracResStage stage : fracResStageList_) {

		std::cout << " -------------- Stage N°" << stage.stageIndex_
				<< " named " << stage.name_ << "------------------" << endl;
		int nb_hypothesis = stage.hypotheisisList_.size();
		if (stage.beddingStage_) {

			std::cout
					<< " -------------- bedding : Bed_interface Hypothesis -------------------"
					<< endl;
			std::cout << "number of bedding hypothesis : " << nb_hypothesis
					<< endl;

		} else {
			std::cout
					<< " ...............Bedding : No bed interface ---------------------------"
					<< endl;
			std::cout << "number of  hypothesis : " << nb_hypothesis << endl;
		}

		for (FracResHypothesis& fracResHypothesis : stage.hypotheisisList_) {

			int nb_bedding_interface = fracResHypothesis.beddingList_.size();
			int nb_bedding_parallel =fracResHypothesis.bed_parallel_List_.size();
			fracResHypothesis.update_active_dfn_number();
			int nb_single_fault = fracResHypothesis.nb_single_fault_active_;
			int nb_array_fault = fracResHypothesis.nb_fault_array_active_;
			int nb_corridor_fault = fracResHypothesis.nb_cluster_corridors_active_;
			int nb_fault_zone = fracResHypothesis.nb_cluster_fault_zone_active_;
			int nb_total = 0;
			if (nb_bedding_parallel > 0) {
				std::cout << "number of bed parallel : " << nb_bedding_parallel
						<< endl;
				nb_total += nb_bedding_parallel;
			}
			if (nb_bedding_interface > 0) {
				std::cout << "number of bed interface : "
						<< nb_bedding_interface << endl;
				nb_total += nb_bedding_interface;
			}
			if (nb_single_fault > 0) {
				std::cout << "number of single fault : " << nb_single_fault
						<< endl;
				nb_total += nb_single_fault;
			}
			if (nb_array_fault > 0) {
				std::cout << "number of array fault : " << nb_array_fault
						<< endl;
				nb_total += nb_array_fault;
			}

			if (nb_corridor_fault > 0) {
				std::cout << "number of clusters corridors  : "
						<< nb_corridor_fault << endl;
				nb_total += nb_corridor_fault;
			}
			if (nb_fault_zone > 0) {
				std::cout << "number of clusters fault zone  : "
						<< nb_fault_zone << endl;
				nb_total += nb_fault_zone;
			}
			std::cout << endl;

			std::cout << "number of dfn  : " << nb_total << " " << endl;
			nb_total_dfn += nb_total;

			int indice_bed_interface = 0;
			int indice_bed_parallel = 0;
			int indice_single_fault = 0;
			int indice_fault_array = 0;
			int indice_fault_corridor = 0;
			int indice_fault_zone = 0;
			for (int type_dfn : fracResHypothesis.type_dfn_) {
				std::chrono::time_point < std::chrono::system_clock
						> total_StartTime;
				std::chrono::time_point < std::chrono::system_clock
						> total_EndTime;

				switch (type_dfn) {
				case 0: {

					FracResBeddingAperture *bed_interface =
							fracResHypothesis.beddingList_[indice_bed_interface];
					if (!bed_interface->is_active_) {
						cout << " the bed interface " << bed_interface->name_
								<< "is inactive " << endl;
						break;
					}
					total_StartTime = std::chrono::system_clock::now();
					bed_interface_property_name_ =
							fracResHypothesis.beddingPropertyDiscreteName_;
					create_bed_surface(model, fracResHypothesis.beddingList_,
							fracResHypothesis.beddingPropertyDiscreteName_, indice_bed_interface);

					nb_total_bed_ = 0;
					for (auto const &stress_geo : stress_geo_list_) {

						cout << "nb fractures for ( " << stress_geo.first.first
								<< " , " << stress_geo.first.second << " )= "
								<< stress_geo.second->tab_fractures_.size()
								<< endl;

						nb_total_bed_ +=
								stress_geo.second->tab_fractures_.size();

						model_fractures_.geo_cont_.add_fracture_set(
								stress_geo.second);
					}
					nb_failles_ =model_fractures_.geo_cont_.tab_fractures_set_.size();
					/*total_EndTime = std::chrono::system_clock::now();
					double write_miliseconde = std::chrono::duration_cast< std::chrono::milliseconds> (total_EndTime - total_StartTime).count();
					double write_second = write_miliseconde / 1000;
					cout << " Temps total de simulation des bed interface slip "<< endl; //<<bed_interface->name_<< endl;
					cout << " nb miliseconds " << write_miliseconde<< " nb second " << write_second << '\n';
					cout << endl;*/
					time_creation_message("",total_StartTime);
					indice_bed_interface++;
					break;
				}
				case 1: {
					FractureSet *setFrac =fracResHypothesis.fractureList_[indice_single_fault];
					if (!setFrac->is_active_) {
						cout << " the single fault " << setFrac->name_<< "is inactive " << endl;
						break;
					}
					total_StartTime = std::chrono::system_clock::now();
					std::cout
							<< " --------------  SINGLE FAULT CREATION IN PROGRESS ------------------"
							<< endl;
					setFrac->fracture_set_index_ += 1000 * stage.stageIndex_ + 100*fracResHypothesis.hypothesisIndex_;
                    setFrac->model_geo_cont_index_ = model_fractures_.geo_cont_.tab_fractures_set_.size();
					setFrac->seed_ = seed_;
					setFrac->set_volume(volume_);
					setFrac->model_fractures_.geo_cont_ =
							model_fractures_.geo_cont_;
					setFrac->set_area(area_);
					setFrac->set_generation_box(generation_box_);
					setFrac->directory_path_ = outpout_directory_path_;

					std::shared_ptr<fault_array> test_array = NULL;
					setFrac->fractures_generation_from_shadow_zone(model,tree_,disteval_,cells_,test_array, estim_volume);

					if (setFrac->stress_zone_set_->tab_fractures_.size() > 0) {
						model_fractures_.geo_cont_.add_fracture_set(
								setFrac->stress_zone_set_);
						if (setFrac->is_shadow_zone_active_) {
							setFrac->type_intersection_ =
									fractures_intersect::FracResIntersectionTypeEnum::CROSSING;
							cout << " Fracture set " << setFrac->name_
									<< " is crossing with shadow zone" << endl;

						} else {
							if (setFrac->type_intersection_
									== fractures_intersect::FracResIntersectionTypeEnum::CROSSING) {
								cout << " Fracture set " << setFrac->name_
										<< " is crossing with no shadow zone"
										<< endl;
							} else if (setFrac->type_intersection_
									== fractures_intersect::FracResIntersectionTypeEnum::ABUTTING) {
								cout << " Fracture set " << setFrac->name_
										<< " is abbuting with no shadow zone"
										<< endl;
							} else {
								cout << " Fracture set " << setFrac->name_
										<< " is randomly crossing and abbuting with no shadow zone"
										<< endl;
							}
						}
						setFrac->fractures_activable_ =
								setFrac->fractures_number_;

						setFrac->fractures_inserted_ = 0;

						if (setFrac->success_) {
							nb_success++;

						}

					}

					/*total_EndTime = std::chrono::system_clock::now();
					double write_miliseconde = std::chrono::duration_cast< std::chrono::milliseconds> (total_EndTime - total_StartTime).count();
					double write_second = write_miliseconde / 1000;
					cout << " Temps total de simulation du dfn "<< setFrac->name_ << endl;
					cout << " nb miliseconds " << write_miliseconde<< " nb second " << write_second << '\n';
					cout << endl;*/
					time_creation_message(setFrac->name_,total_StartTime);
					indice_single_fault++;
					nb_failles_++;
					break;
				}
				case 2: {

					FaultArraySet *setArrayZone =
							fracResHypothesis.faultArrayList_[indice_fault_array];
					if (!setArrayZone->is_active_) {
						cout << " the fault array " << setArrayZone->name_
								<< "is inactive " << endl;
						break;
					}
					total_StartTime = std::chrono::system_clock::now();
					std::cout
							<< " --------------  FAULT ARRAY CREATION IN PROGRESS ------------------"
							<< endl;

					setArrayZone->fault_array_set_index_ += 1000 * stage.stageIndex_ + 100*fracResHypothesis.hypothesisIndex_;
							//model_fractures_.geo_cont_.tab_fractures_set_.size();
					setArrayZone->model_geo_cont_index_ = model_fractures_.geo_cont_.tab_fractures_set_.size();
					setArrayZone->seed_ = seed_;
					setArrayZone->volume_ = volume_;
					setArrayZone->fault_array_geo_set_ = std::shared_ptr
							< fractures_intersect::fault_arrays_set_geo
							> (new fractures_intersect::fault_arrays_set_geo(
									setArrayZone->fault_array_set_index_,
									setArrayZone->seed_));
					setArrayZone->area_ = area_;
					setArrayZone->set_generation_box(generation_box_);
					setArrayZone->fault_array_generation(model,tree_,disteval_, cells_,model_fractures_.geo_cont_,estim_volume);
					setArrayZone->directory_path_ = outpout_directory_path_;
					if (setArrayZone->success_) {
						nb_success++;

					}
					/*total_EndTime = std::chrono::system_clock::now();
					double write_miliseconde = std::chrono::duration_cast
							< std::chrono::milliseconds
							> (total_EndTime - total_StartTime).count();
					double write_second = write_miliseconde / 1000;
					cout << " Temps total de simulation du dfn "
							<< setArrayZone->name_ << endl;
					cout << " nb miliseconds " << write_miliseconde
							<< " nb second " << write_second << '\n';*/
					time_creation_message(setArrayZone->name_,total_StartTime);
					indice_fault_array++;
					nb_failles_++;
					break;
				}
				case 3: {

					FracResBeddingProperties *bed_parallel =
							fracResHypothesis.bed_parallel_List_[indice_bed_parallel];

					total_StartTime = std::chrono::system_clock::now();
					bool BedInterfacePropertyIsSameAsBedParallelProperty = false;
					if (bed_interface_property_name_
							== bed_parallel->bed_parrallel_discrete_property_name_) {
						BedInterfacePropertyIsSameAsBedParallelProperty = true;
					}

					if (BedInterfacePropertyIsSameAsBedParallelProperty) {
						if (nb_total_bed_ > 0) {
							update_bed_interface(model,model_fractures_, seed_, bed_parallel);
						}
					} else {
						create_bed_parallel(model, seed_, bed_parallel);
						for (auto const &stress_geo_parallel : stress_geo_parallel_list_) {
							model_fractures_.geo_cont_.add_fracture_set(
									stress_geo_parallel.second);
							nb_failles_++;

						}
					}
					bed_parallel->nb_fracture_set_geo_ =
							stress_geo_parallel_list_.size();
					bed_parallel->index_fracture_set_geo_begin_ =
							model_fractures_.geo_cont_.tab_fractures_set_.size()
									- 1;
					/*total_EndTime = std::chrono::system_clock::now();
					double write_miliseconde = std::chrono::duration_cast
							< std::chrono::milliseconds
							> (total_EndTime - total_StartTime).count();
					double write_second = write_miliseconde / 1000;
					cout << " Temps total de simulation du dfn "
							<< bed_parallel->name_ << endl;
					cout << " nb miliseconds " << write_miliseconde
							<< " nb second " << write_second << '\n';*/
					time_creation_message(bed_parallel->name_ ,total_StartTime);
					indice_bed_parallel++;
					break;
				}
				case 4: {
					FaultArraySet *setArrayCorridor =
							fracResHypothesis.faultCorridorList_[indice_fault_corridor];
					if (!setArrayCorridor->is_active_) {
						cout << " the cluster fault corridor "
								<< setArrayCorridor->name_ << "is inactive "
								<< endl;
						break;
					}
					total_StartTime = std::chrono::system_clock::now();
					std::cout
							<< " --------------  FAULT CORRIDOR CREATION IN PROGRESS ------------------"
							<< endl;
					setArrayCorridor->fault_array_set_index_ +=1000 * stage.stageIndex_ + 100*fracResHypothesis.hypothesisIndex_;
							//model_fractures_.geo_cont_.tab_fractures_set_.size();
					setArrayCorridor->model_geo_cont_index_ = model_fractures_.geo_cont_.tab_fractures_set_.size();
					setArrayCorridor->seed_ = seed_;
					setArrayCorridor->volume_ = volume_;
					setArrayCorridor->fault_array_geo_set_ = std::shared_ptr< fractures_intersect::fault_arrays_set_geo> (new fractures_intersect::fault_arrays_set_geo(
									setArrayCorridor->fault_array_set_index_,setArrayCorridor->seed_));
					setArrayCorridor->area_ = area_;
					setArrayCorridor->set_generation_box(generation_box_);
					setArrayCorridor->fault_corridor_generation(model,tree_,disteval_,cells_,model_fractures_.geo_cont_,
							estim_volume);
					setArrayCorridor->directory_path_ = outpout_directory_path_;
					if (setArrayCorridor->success_) {
						nb_success++;

					}
					/*total_EndTime = std::chrono::system_clock::now();
					double write_miliseconde = std::chrono::duration_cast
							< std::chrono::milliseconds
							> (total_EndTime - total_StartTime).count();
					double write_second = write_miliseconde / 1000;
					cout << " Temps total de simulation du dfn "
							<< setArrayCorridor->name_ << endl;
					cout << " nb miliseconds " << write_miliseconde
							<< " nb second " << write_second << '\n';*/
					time_creation_message(setArrayCorridor->name_,total_StartTime);
					indice_fault_corridor++;
					nb_failles_++;
					break;
				}
				case 5: {
					FaultArraySet *setArrayZone =
							fracResHypothesis.faultZoneList_[indice_fault_zone];
					if (!setArrayZone->is_active_) {
						cout << " the cluster fault zone "
								<< setArrayZone->name_ << "is inactive "
								<< endl;
						break;
					}
					total_StartTime = std::chrono::system_clock::now();
					std::cout
							<< " --------------  FAULT ZONE CREATION IN PROGRESS ------------------"
							<< endl;

					setArrayZone->fault_array_set_index_ +=1000 * stage.stageIndex_ + 100*fracResHypothesis.hypothesisIndex_;
							//model_fractures_.geo_cont_.tab_fractures_set_.size();
					setArrayZone->model_geo_cont_index_ = model_fractures_.geo_cont_.tab_fractures_set_.size();
					setArrayZone->seed_ = seed_;
					setArrayZone->volume_ = volume_;
					setArrayZone->fault_array_geo_set_ = std::shared_ptr
							< fractures_intersect::fault_arrays_set_geo
							> (new fractures_intersect::fault_arrays_set_geo(
									setArrayZone->fault_array_set_index_,
									setArrayZone->seed_));
					setArrayZone->area_ = area_;
					setArrayZone->set_generation_box(generation_box_);
					setArrayZone->fault_zone_generation(model, tree_,disteval_,cells_,model_fractures_.geo_cont_,estim_volume);
					setArrayZone->directory_path_ = outpout_directory_path_;
					if (setArrayZone->success_) {
						nb_success++;

					}
					/*total_EndTime = std::chrono::system_clock::now();
					double write_miliseconde = std::chrono::duration_cast
							< std::chrono::milliseconds
							> (total_EndTime - total_StartTime).count();
					double write_second = write_miliseconde / 1000;
					cout << " Temps total de simulation du dfn "
							<< setArrayZone->name_ << endl;
					cout << " nb miliseconds " << write_miliseconde
							<< " nb second " << write_second << '\n';*/
					time_creation_message(setArrayZone->name_,total_StartTime);
					indice_fault_zone++;
					//nb_failles_++;
					break;
				}
				default: {

				}

				} // switch

			} // loop type_dfn

		} //loop Hypothesis
	} // loop stage

	if (nb_success == (nb_total_dfn))
		success_ = true;

}







/**
 * @fn std::vector<std::shared_ptr<fractures_intersect::fracture_set_geo>> create_bed_surface(const geode::StructuralModel&, FracResBeddingAperture*)
 * @brief
 *
 * @param model
 * @param bed_interface
 * @return
 */

//std::map< std::pair<int,int>, std::shared_ptr<fractures_intersect::fracture_set_geo> >
void RemFracResSim::create_bed_surface(const geode::StructuralModel &model,
		std::vector<FracResBeddingAperture*> bed_interface_list,
		std::string beddingPropertyName, int indice_bed_interface) {

	std::vector<std::shared_ptr<fractures_intersect::fracture_set_geo> > stress_sets;

	//std::map< std::pair<int,int> ,std::shared_ptr<fractures_intersect::fracture_set_geo> > stress_geo_list_;

	std::map<std::pair<int, int>, double> aperture_val_map;
	std::map<std::pair<int, int>, string> bed_name_map;

	std::map<std::pair<int, int>, int> bed_interface_map;

	//std::vector< fractures_intersect::FracResIntersectionTypeEnum> tab_type_intersection;
	std::vector<fractures_intersect::FracResIntersectionTypeEnum> tab_type_intersection(
			7, fractures_intersect::FracResIntersectionTypeEnum::CROSSING);
	//stress_zone_sets.resize(bed_interface->discrete_property_selected_values_.size());

	int nb_interface_pair_active = 0;
	int nb_interface_active = 0;
	aperture_val_map.clear();
	//for (FracResBeddingAperture *bed_interface : bed_interface_list) {


	FracResBeddingAperture *bed_interface = bed_interface_list[indice_bed_interface];
	if (bed_interface->is_active_) {
			cout << " the bed interface " << bed_interface->name_ << " is active "
					<< endl;
	}

		for (int j = 0;j < bed_interface->discrete_property_selected_values_.size();j++) {

			if (beddings_type_intersection_[bed_interface->name_].size() <= 0)beddings_type_intersection_[bed_interface->name_] =tab_type_intersection;
			std::shared_ptr<fractures_intersect::fracture_set_geo> stress_zone_set =
					std::shared_ptr < fractures_intersect::fracture_set_geo
							> (new fractures_intersect::fracture_set_geo(j,
									fractures_intersect::FracResIntersectionCategoryEnum::BED_INTERFACE,
									seed_,
									beddings_type_intersection_[bed_interface->name_]));
			//stress_sets.push_back(stress_zone_set);
			map<std::pair<int, int>,std::shared_ptr<fractures_intersect::fracture_set_geo>>::iterator it =stress_geo_list_.find(bed_interface->discrete_property_selected_values_[j]);
			if (it == stress_geo_list_.end()) {
				stress_geo_list_.insert({ bed_interface->discrete_property_selected_values_[j],stress_zone_set });
				bed_name_vect_.insert({ bed_interface->discrete_property_selected_values_[j],bed_interface->name_});
			}else{
				stress_geo_list_[bed_interface->discrete_property_selected_values_[j]] = stress_zone_set;
				bed_name_vect_[bed_interface->discrete_property_selected_values_[j]] = bed_interface->name_;
			}
			aperture_val_map.insert({ bed_interface->discrete_property_selected_values_[j],bed_interface->value_ });
			bed_name_map.insert({ bed_interface->discrete_property_selected_values_[j],bed_interface->name_ });
			bed_interface_map.insert({ bed_interface->discrete_property_selected_values_[j],nb_interface_active });
			nb_interface_pair_active++;
		}
		//nb_interface_active=indice_bed_interface +1;
	//}
	boost::mt19937 center_gen;
	center_gen.seed(seed_);
	boost::function<double()> randSlope;

	int block_index = 0;
	int nb_bedinterface = nb_interface_pair_active;
	for (const auto &block : model.blocks()) {

		absl::flat_hash_map<geode::PolyhedronFacet, bool> map_polyhedronfacet_active;

		// store name and uuid of block object

		// get block mesh object type as TetrahedralSolid3D not Polyedral
		const auto &block_mesh = block.mesh<geode::TetrahedralSolid3D>();

		// enable edges for edge computation and create property for checking edge

		//const auto attr_edge_from = block_mesh.edges().edge_attribute_manager().find_or_create_attribute< geode::VariableAttribute,double >( "edge_id", ndvalue_ );

		// check petrophysique properties exist
		bool testPoperty =block_mesh.polyhedron_attribute_manager().attribute_exists(beddingPropertyName);

		if (!testPoperty) {
			cout << " The discrete Property " << beddingPropertyName
					<< " doesn't exist for block " << block.name() << endl;
			continue;
		}
		// get petrophysique properties from polyhedron
		const auto &attributProperty_ =block_mesh.polyhedron_attribute_manager().find_or_create_attribute<geode::VariableAttribute, double>(beddingPropertyName,ndvalue_);

		if (!block_mesh.are_facets_enabled())block_mesh.enable_facets();

		for (const auto facet_id : geode::Range {block_mesh.facets().nb_facets() }) {

			const PolyhedraAroundFacet polyAround =block_mesh.polyhedra_from_facet_vertices(block_mesh.facets().facet_vertices(facet_id));

			if (polyAround.size() < 2) continue;
			int fault_val = attributProperty_->value(polyAround[0].polyhedron_id);
			int fault_adj_val = attributProperty_->value(polyAround[1].polyhedron_id);


			if ((fault_val >= 0 && fault_adj_val >= 0) && (fault_val != fault_adj_val)) {



				std::pair<int, int> pairVal(std::min(fault_val, fault_adj_val),std::max(fault_val, fault_adj_val));


				map<std::pair<int, int>,bool>::iterator i =active_discrete_pair_array_.find(pairVal);
				map<std::pair<int, int>,std::string>::iterator ibed =bed_name_vect_.find(pairVal);
				if (i == active_discrete_pair_array_.end() && ibed == bed_name_vect_.end()) {

					fracResBeddingzero_->value_=0.0;

					fracResBeddingzero_->discrete_property_selected_values_.push_back(pairVal);
					std::shared_ptr<fractures_intersect::fracture_set_geo> stress_zone_set =std::shared_ptr< fractures_intersect::fracture_set_geo> (new fractures_intersect::fracture_set_geo(nb_bedinterface,
											fractures_intersect::FracResIntersectionCategoryEnum::BED_INTERFACE,
											seed_, tab_type_intersection));
					aperture_val_map.insert( { pairVal, 0 });
					bed_name_map.insert( { pairVal, "bed_interface" });
					stress_geo_list_.insert( { pairVal, stress_zone_set });
					bed_name_vect_.insert({ pairVal,"bed_interface"});
					//bed_interface->index_end_++;
					nb_bedinterface++;

				}


				if(bed_name_map.find(pairVal) == bed_name_map.end()){
					continue;
				}


				const geode::Triangle3D tri = block_mesh.triangle(polyAround[0]);
				geode::Point3D p0 = tri.vertices()[0].get();
				geode::Point3D p1 = tri.vertices()[1].get();
				geode::Point3D p2 = tri.vertices()[2].get();

				fractures_intersect::Point3D germe_pt0 =fractures_intersect::Point3D(p0.value(0), p0.value(1),p0.value(2));
				fractures_intersect::Point3D germe_pt1 =fractures_intersect::Point3D(p1.value(0), p1.value(1),p1.value(2));
				fractures_intersect::Point3D germe_pt2 =fractures_intersect::Point3D(p2.value(0), p2.value(1),p2.value(2));

				std::vector<fractures_intersect::Point3D> nodes;
				nodes.push_back(germe_pt0);
				nodes.push_back(germe_pt1);
				nodes.push_back(germe_pt2);


				if (bed_name_vect_.find(pairVal)!= bed_name_vect_.end()) {

					FracResBeddingAperture* bed_interface =bed_interface_list[indice_bed_interface];
					aperture_val_map[pairVal] = bed_interface->value_;

					if (bed_interface->dataType_== FracResDataTypeEnum::DISTRIBUTION) {
						initialize_bed_aperture_distribution(center_gen,randSlope, bed_interface->min_mode_max_,bed_interface->distribution_name_);

						aperture_val_map[pairVal] = randSlope();
					} else if (bed_interface->dataType_== FracResDataTypeEnum::PROPERTY) {
						std::shared_ptr<geode::VariableAttribute<double> > attributApertureProperty;

						// check petrophysique properties exist
						bool testPopertyAperture =block_mesh.polyhedron_attribute_manager().attribute_exists(bed_interface->discrete_property_name_);
						if (!testPopertyAperture) continue;

						// get petrophysique properties from polyhedron
						attributApertureProperty =block_mesh.polyhedron_attribute_manager().find_or_create_attribute<geode::VariableAttribute, double>(bed_interface->discrete_property_name_,ndvalue_);

						aperture_val_map[pairVal] =(attributApertureProperty->value(polyAround[0].polyhedron_id)+ attributApertureProperty->value(polyAround[1].polyhedron_id))/ 2;

					}

				}
				double value = aperture_val_map[pairVal];
				if(value > 0.5){
					int ducon= 0;
					ducon++;
				}

				std::shared_ptr<fractures_intersect::Cfaille> faille =std::shared_ptr < fractures_intersect::Cfaille> (new fractures_intersect::Cfaille(stress_geo_list_[pairVal]->tab_fractures_.size(),value, nodes, bed_name_map[pairVal])); //stress_geo_list_[pairVal]->index_

				std::shared_ptr<fractures_intersect::FracResBeddingLocalisation> localisation =std::shared_ptr< fractures_intersect::FracResBeddingLocalisation> (new fractures_intersect::FracResBeddingLocalisation(stress_geo_list_[pairVal]->tab_fractures_.size()
												- 1,
										polyAround[0].polyhedron_id,
										polyAround[1].polyhedron_id, block.id()));

				stress_geo_list_[pairVal]->add_bedding_localisation(localisation);
				stress_geo_list_[pairVal]->add_fracture(faille);

			}

		}

		/*
		 for( const auto p : geode::Range{ block_mesh.nb_polyhedra() } )
		 {

		 int fault_val = attributProperty_->value(p);
		 // initialize attribut Facies value if neede
		 if(fault_val < 0) continue;


		 int indice=0;
		 for( std::pair<int,int> pairSet : bed_interface->discrete_property_selected_values_){

		 for(auto faceId : geode::LRange(block_mesh.nb_polyhedron_facets(p))){

		 geode::PolyhedronFacet facet{p,faceId};

		 if(!map_polyhedronfacet_active[facet]) map_polyhedronfacet_active.insert({facet, false});

		 const auto opt_polyhedron_adj = block_mesh.polyhedron_adjacent( facet );
		 if( !opt_polyhedron_adj )
		 {
		 int ducon=0;
		 continue;

		 }
		 if( const auto adj = block_mesh.polyhedron_adjacent(facet) ){


		 geode::PolyhedronFacet adjfact = block_mesh.polyhedron_adjacent_facet( facet ).value();

		 if(!map_polyhedronfacet_active[adjfact]) map_polyhedronfacet_active.insert({adjfact, false});

		 int fault_val_adj = attributProperty_->value(adjfact.polyhedron_id);

		 if(fault_val==4 && fault_val_adj == 16){
		 int ducon=0;
		 ducon++;
		 }

		 if(fault_val == pairSet.first && fault_val_adj == pairSet.second){

		 map_polyhedronfacet_active[facet]= true;
		 map_polyhedronfacet_active[adjfact]= true;

		 const geode::Triangle3D tri = block_mesh.triangle(adjfact);

		 geode::Point3D p0= tri.vertices()[0].get() ;
		 geode::Point3D p1= tri.vertices()[1].get() ;
		 geode::Point3D p2= tri.vertices()[2].get() ;

		 fractures_intersect::Point3D germe_pt0 = fractures_intersect::Point3D(p0.value(0),p0.value(1),p0.value(2));
		 fractures_intersect::Point3D germe_pt1 = fractures_intersect::Point3D(p1.value(0),p1.value(1),p1.value(2));
		 fractures_intersect::Point3D germe_pt2 = fractures_intersect::Point3D(p2.value(0),p2.value(1),p2.value(2));

		 std::vector< fractures_intersect::Point3D > nodes ;
		 nodes.push_back(germe_pt0);
		 nodes.push_back(germe_pt1);
		 nodes.push_back(germe_pt2);
		 double value = bed_interface->value_;

		 if (bed_interface->dataType_ == FracResDataTypeEnum::DISTRIBUTION){
		 initialize_bed_aperture_distribution( center_gen,randSlope, bed_interface->min_mode_max_ ,bed_interface->distribution_name_);

		 value = randSlope();
		 }else if (bed_interface->dataType_ == FracResDataTypeEnum::PROPERTY){
		 value = (attributApertureProperty->value(adjfact.polyhedron_id)+  attributApertureProperty->value(p))/2;

		 }





		 std::shared_ptr<fractures_intersect::Cfaille> faille = std::shared_ptr<fractures_intersect::Cfaille>(new fractures_intersect::Cfaille(nb_failles_+  stress_sets[indice]->tab_fractures_.size(),value,nodes,bed_interface->name_));

		 stress_sets[indice]->add_fracture(faille);	;
		 }
		 }

		 }

		 indice++;
		 }
		 }*/

		block_index++;

	}

	//return stress_sets;*/

	//return stress_geo_list_;

}
/**
 * @fn void initialize_slope_distribution(boost::mt19937&, boost::function<double ()>&, std::vector<double>, std::string)
 * @brief
 *
 * @param center_gen
 * @param rand
 * @param minmaxmode_values
 * @param type_name
 */

void RemFracResSim::initialize_aperture_slope_distribution(
		boost::mt19937 &center_gen, boost::function<double()> &rand,
		std::vector<double> minmaxmode_values, std::string type_name) {

	FractureSet::distribution_type type_slope =
			FractureSet::distribution_type::uniform;

	if (type_name == "TRIANGULAR") {

		type_slope = FractureSet::distribution_type::triangulated;
	} else if (type_name == "NORMAL") {
		type_slope = FractureSet::distribution_type::normal;
	} else if (type_name == "LOGNORMAL") {
		type_slope = FractureSet::distribution_type::lognormal;
	}

	rand = boost::bind(
			boost::random::uniform_real_distribution<>(minmaxmode_values[0],
					minmaxmode_values[1]), center_gen);
	if (type_slope == FractureSet::distribution_type::triangulated) {
		rand = boost::bind(
				boost::random::triangle_distribution<>(minmaxmode_values[0],
						minmaxmode_values[2], minmaxmode_values[1]),
				center_gen);
	} else if (type_slope == FractureSet::distribution_type::normal) {
		rand = boost::bind(
				boost::random::normal_distribution<>(minmaxmode_values[0],
						minmaxmode_values[1]), center_gen);
	} else if (type_slope == FractureSet::distribution_type::lognormal) {
		rand = boost::bind(
				boost::random::lognormal_distribution<>(minmaxmode_values[0],
						minmaxmode_values[1]), center_gen);
	}
}

void RemFracResSim::initialize_bed_aperture_distribution(
		boost::mt19937 &center_gen, boost::function<double()> &rand,
		std::vector<double> minmaxmode_values, std::string type_name) {

	FractureSet::distribution_type type_bed_aperture =
			FractureSet::distribution_type::uniform;

	if (type_name == "TRIANGULAR") {

		type_bed_aperture = FractureSet::distribution_type::triangulated;
	} else if (type_name == "NORMAL") {
		type_bed_aperture = FractureSet::distribution_type::normal;
	} else if (type_name == "LOGNORMAL") {
		type_bed_aperture = FractureSet::distribution_type::lognormal;
	}

	rand = boost::bind(
			boost::random::uniform_real_distribution<>(minmaxmode_values[0],
					minmaxmode_values[1]), center_gen);
	if (type_bed_aperture == FractureSet::distribution_type::triangulated) {
		rand = boost::bind(
				boost::random::triangle_distribution<>(minmaxmode_values[0],
						minmaxmode_values[2], minmaxmode_values[1]),
				center_gen);
	} else if (type_bed_aperture == FractureSet::distribution_type::normal) {
		rand = boost::bind(
				boost::random::normal_distribution<>(minmaxmode_values[0],
						minmaxmode_values[1]), center_gen);
	} else if (type_bed_aperture == FractureSet::distribution_type::lognormal) {
		rand = boost::bind(
				boost::random::lognormal_distribution<>(minmaxmode_values[0],
						minmaxmode_values[1]), center_gen);
	}
}

/**
 * @fn std::vector<std::vector<std::shared_ptr<fractures_intersect::fracture_set_geo>>> create_bed_parallel(const geode::StructuralModel&, boost::uint32_t, FracResBeddingProperties*)
 * @brief
 *
 * @param model
 * @param seed
 * @param bed_parallel
 * @return
 */

void RemFracResSim::create_bed_parallel(const geode::StructuralModel &model,
		boost::uint32_t seed, FracResBeddingProperties *bed_parallel) {
	std::map<std::pair<int, int>, double> aperture_val_map;
	std::map<std::pair<int, int>, string> bed_name_map;

	std::map<std::pair<int, int>, int> bed_interface_map;
	bed_name_parallel_vect_.clear();

	//std::vector< fractures_intersect::FracResIntersectionTypeEnum> tab_type_intersection;
	std::vector<fractures_intersect::FracResIntersectionTypeEnum> tab_type_intersection(
			7, fractures_intersect::FracResIntersectionTypeEnum::CROSSING);
	//stress_zone_sets.resize(bed_interface->discrete_property_selected_values_.size());

	int nb_interface_pair_active = 0;
	int nb_interface_active = 0;

	double slope_azimut = bed_parallel->slope_data_value_;
	double ratio_perm = 1.0;
	boost::mt19937 center_gen;
	center_gen.seed(seed);

	ratio_perm = bed_parallel->slope_ratio_data_value_;
	boost::function<double()> randSlope;

	std::string slope_property_name = "";

	for (FracResBeddingAperture *fracTest : bed_parallel->bedding_list_) {

		for (int j = 0; j < fracTest->discrete_property_selected_values_.size();
				j++) {

			if (beddings_type_intersection_[fracTest->name_].size() <= 0)
				beddings_type_intersection_[fracTest->name_] =
						tab_type_intersection;

			std::shared_ptr<fractures_intersect::fracture_set_geo> stress_zone_set =
					std::shared_ptr < fractures_intersect::fracture_set_geo
							> (new fractures_intersect::fracture_set_geo(j,
									fractures_intersect::FracResIntersectionCategoryEnum::BED_PARALLEL,
									seed_,
									beddings_type_intersection_[fracTest->name_]));
			stress_geo_parallel_list_.insert(
					{ fracTest->discrete_property_selected_values_[j],
							stress_zone_set });
			aperture_val_map.insert(
					{ fracTest->discrete_property_selected_values_[j],
							fracTest->value_ });
			bed_name_map.insert(
					{ fracTest->discrete_property_selected_values_[j],
							fracTest->name_ });
			bed_interface_map.insert(
					{ fracTest->discrete_property_selected_values_[j],
							nb_interface_active });
			nb_interface_pair_active++;
			bed_name_parallel_vect_.insert(
					{ fracTest->discrete_property_selected_values_[j],
							fracTest->name_ });
		}
		nb_interface_active++;

	}

	int block_index = 0;
	int nb_bedinterface = nb_interface_pair_active;

	for (const auto &block : model.blocks()) {

		absl::flat_hash_map<geode::PolyhedronFacet, bool> map_polyhedronfacet_active;

		// store name and uuid of block object

		// get block mesh object type as TetrahedralSolid3D not Polyedral
		const auto &block_mesh = block.mesh<geode::TetrahedralSolid3D>();

		if (bed_parallel->is_region_
				&& block.name() != bed_parallel->region_name_)
			continue;

		// check petrophysique properties exist

		bool testPoperty =
				block_mesh.polyhedron_attribute_manager().attribute_exists(
						bed_parallel->bed_parrallel_discrete_property_name_);
		if (!testPoperty)
			continue;
		// get petrophysique properties from polyhedron
		const auto &attributProperty =
				block_mesh.polyhedron_attribute_manager().find_or_create_attribute<
						geode::VariableAttribute, double>(
						bed_parallel->bed_parrallel_discrete_property_name_,
						ndvalue_);

		// get petrophysique properties from polyhedron
		std::shared_ptr<geode::VariableAttribute<double> > attributSlopeProperty;
		if (bed_parallel->slope_data_type_ == FracResDataTypeEnum::PROPERTY) {
			bool testSlopePoperty =
					block_mesh.polyhedron_attribute_manager().attribute_exists(
							bed_parallel->bed_parrallel_slope_prop_name_);
			const auto &attributSlopeProperty =
					block_mesh.polyhedron_attribute_manager().find_or_create_attribute<
							geode::VariableAttribute, double>(
							bed_parallel->bed_parrallel_slope_prop_name_,
							ndvalue_);
			if (!testSlopePoperty)
				continue;
		}

		if (!block_mesh.are_facets_enabled())
			block_mesh.enable_facets();

		for (const auto facet_id : geode::Range {
				block_mesh.facets().nb_facets() }) {
			const PolyhedraAroundFacet polyAround =
					block_mesh.polyhedra_from_facet_vertices(
							block_mesh.facets().facet_vertices(facet_id));

			if (polyAround.size() < 2)
				continue;
			int fault_val = attributProperty->value(
					polyAround[0].polyhedron_id);
			int fault_adj_val = attributProperty->value(
					polyAround[1].polyhedron_id);

			if ((fault_val >= 0 && fault_adj_val > 0)
					&& (fault_val != fault_adj_val)) {

				std::pair<int, int> pairVal(std::min(fault_val, fault_adj_val),
						std::max(fault_val, fault_adj_val));

				FracResBeddingAperture *fracAperture =
						bed_parallel->bedding_list_[bed_interface_map[pairVal]];
				std::shared_ptr<geode::VariableAttribute<double> > attributApertureProperty;
				if (fracAperture->dataType_ == FracResDataTypeEnum::PROPERTY) {
					// check petrophysique properties exist
					bool testPopertyAperture =
							block_mesh.polyhedron_attribute_manager().attribute_exists(
									fracAperture->discrete_property_name_);
					if (!testPopertyAperture)
						continue;

					// get petrophysique properties from polyhedron
					attributApertureProperty =
							block_mesh.polyhedron_attribute_manager().find_or_create_attribute<
									geode::VariableAttribute, double>(
									fracAperture->discrete_property_name_,
									ndvalue_);
				}

				const geode::Triangle3D tri = block_mesh.triangle(
						polyAround[0]);
				geode::Point3D p0 = tri.vertices()[0].get();
				geode::Point3D p1 = tri.vertices()[1].get();
				geode::Point3D p2 = tri.vertices()[2].get();

				fractures_intersect::Point3D germe_pt0 =
						fractures_intersect::Point3D(p0.value(0), p0.value(1),
								p0.value(2));
				fractures_intersect::Point3D germe_pt1 =
						fractures_intersect::Point3D(p1.value(0), p1.value(1),
								p1.value(2));
				fractures_intersect::Point3D germe_pt2 =
						fractures_intersect::Point3D(p2.value(0), p2.value(1),
								p2.value(2));

				std::vector<fractures_intersect::Point3D> nodes;
				nodes.push_back(germe_pt0);
				nodes.push_back(germe_pt1);
				nodes.push_back(germe_pt2);
				double value = aperture_val_map[pairVal];

				if (bed_interface_map.find(pairVal)
						!= bed_interface_map.end()) {

					FracResBeddingAperture *bed_interface =
							bed_parallel->bedding_list_[bed_interface_map[pairVal]];

					if (bed_interface->dataType_
							== FracResDataTypeEnum::DISTRIBUTION) {
						initialize_bed_aperture_distribution(center_gen,
								randSlope, fracAperture->min_mode_max_,
								fracAperture->distribution_name_);

						aperture_val_map[pairVal] = randSlope();
					} else if (bed_interface->dataType_
							== FracResDataTypeEnum::PROPERTY) {
						std::shared_ptr<geode::VariableAttribute<double> > attributApertureProperty;

						// check petrophysique properties exist
						bool testPopertyAperture =
								block_mesh.polyhedron_attribute_manager().attribute_exists(
										fracAperture->discrete_property_name_);
						if (!testPopertyAperture)
							continue;

						// get petrophysique properties from polyhedron
						attributApertureProperty =
								block_mesh.polyhedron_attribute_manager().find_or_create_attribute<
										geode::VariableAttribute, double>(
										fracAperture->discrete_property_name_,
										ndvalue_);

					}

				}

				value = fracAperture->value_;

				if (fracAperture->dataType_
						== FracResDataTypeEnum::DISTRIBUTION) {
					initialize_bed_aperture_distribution(center_gen, randSlope,
							fracAperture->min_mode_max_,
							fracAperture->distribution_name_);

					value = randSlope();
				} else if (fracAperture->dataType_
						== FracResDataTypeEnum::PROPERTY) {
					value = (attributApertureProperty->value(
							polyAround[0].polyhedron_id)
							+ attributApertureProperty->value(
									polyAround[1].polyhedron_id)) / 2;

				}
				if (bed_parallel->slope_data_type_
						== FracResDataTypeEnum::DISTRIBUTION) {
					initialize_aperture_slope_distribution(center_gen,
							randSlope, bed_parallel->min_mode_max_,
							bed_parallel->bed_parrallel_slope_distribution_type_name_);
					slope_azimut = randSlope();
				} else if (bed_parallel->slope_data_type_
						== FracResDataTypeEnum::PROPERTY) {
					slope_azimut = (attributSlopeProperty->value(
							polyAround[0].polyhedron_id)
							+ attributSlopeProperty->value(
									polyAround[1].polyhedron_id)) / 2;

				}
				std::shared_ptr<fractures_intersect::Cfaille> faille =
						std::shared_ptr < fractures_intersect::Cfaille
								> (new fractures_intersect::Cfaille(
										stress_geo_parallel_list_[pairVal]->index_,
										value, nodes, fracAperture->name_));

				faille->Set_SlopeAzimut(slope_azimut);
				faille->Set_RatioPermeability(ratio_perm);
				stress_geo_parallel_list_[pairVal]->add_fracture(faille);

			}

		}

		block_index++;

	}

}
/**
 * @fn void update_bed_interface(const geode::StructuralModel&, boost::uint32_t, FracResBeddingProperties*)
 * @brief
 *
 * @param model
 * @param seed
 * @param bed_parallel
 * @return
 */

void RemFracResSim::update_bed_interface(const geode::StructuralModel &model,fractures_intersect::model& model_fractures_loc,
		boost::uint32_t seed, FracResBeddingProperties *bed_parallel) {

	double slope_azimut = 0.0;
	double ratio_perm = 1.0;
	boost::mt19937 center_gen;
	center_gen.seed(seed);

	ratio_perm = bed_parallel->slope_ratio_data_value_;
	boost::function<double()> randSlope;

	std::string slope_property_name = "";
	std::shared_ptr<geode::VariableAttribute<double> > attributApertureProperty;
	std::shared_ptr<geode::VariableAttribute<double> > attributSlopeProperty;


	map<std::pair<int, int>,std::shared_ptr<fractures_intersect::fracture_set_geo>>::iterator i_stress_geo ;

	for (FracResBeddingAperture *fracTest : bed_parallel->bedding_list_) {

		slope_azimut = bed_parallel->slope_data_value_;

		for (int j = fracTest->index_begin_;j < fracTest->index_begin_ + fracTest->discrete_property_selected_values_.size();j++) {

			std::pair<int, int> pairVal =fracTest->discrete_property_selected_values_[j];


			i_stress_geo =stress_geo_list_.find(pairVal);
			if (i_stress_geo == stress_geo_list_.end()) {
				continue;
			}

			for (std::shared_ptr<fractures_intersect::FracResBeddingLocalisation> localisation : model_fractures_loc.geo_cont_.tab_fractures_set_[j] ->tab_bedding_localisation_) {

				int indice_poly1 = localisation->index_polyedre_1_;
				int indice_poly2 = localisation->index_polyedre_2_;
				const auto &block_mesh =
						model.block(localisation->block_id_).mesh<
								geode::TetrahedralSolid3D>();

				double aperture = fracTest->value_;
				if (fracTest->dataType_ == FracResDataTypeEnum::PROPERTY) {
					// check petrophysique properties exist
					bool testPopertyAperture =
							block_mesh.polyhedron_attribute_manager().attribute_exists(
									fracTest->discrete_property_name_);
					if (!testPopertyAperture)
						continue;

					// get petrophysique properties from polyhedron
					attributApertureProperty =
							block_mesh.polyhedron_attribute_manager().find_or_create_attribute<
									geode::VariableAttribute, double>(
									fracTest->discrete_property_name_,
									ndvalue_);
					aperture = (attributApertureProperty->value(indice_poly1)
							+ attributApertureProperty->value(indice_poly1))
							/ 2;

				} else if (fracTest->dataType_
						== FracResDataTypeEnum::DISTRIBUTION) {
					initialize_bed_aperture_distribution(center_gen, randSlope,
							fracTest->min_mode_max_,
							fracTest->distribution_name_);

					aperture = randSlope();
				}

				if (bed_parallel->slope_data_type_
						== FracResDataTypeEnum::PROPERTY) {
					bool testSlopePoperty =
							block_mesh.polyhedron_attribute_manager().attribute_exists(
									bed_parallel->bed_parrallel_slope_prop_name_);
					const auto &attributSlopeProperty =
							block_mesh.polyhedron_attribute_manager().find_or_create_attribute<
									geode::VariableAttribute, double>(
									bed_parallel->bed_parrallel_slope_prop_name_,
									ndvalue_);
					if (!testSlopePoperty)
						continue;
					slope_azimut = (attributSlopeProperty->value(indice_poly1)
							+ attributSlopeProperty->value(indice_poly2)) / 2;

				} else if (bed_parallel->slope_data_type_
						== FracResDataTypeEnum::DISTRIBUTION) {
					initialize_aperture_slope_distribution(center_gen,
							randSlope, bed_parallel->min_mode_max_,
							bed_parallel->bed_parrallel_slope_distribution_type_name_);
					slope_azimut = randSlope();
				}

				int indice_faille = localisation->index_;

				std::shared_ptr<fractures_intersect::Cfaille> faille =
						model_fractures_loc.geo_cont_.tab_fractures_set_[j]->tab_fractures_[indice_faille];

				model_fractures_loc.geo_cont_.tab_fractures_set_[j]->tab_fractures_[indice_faille]->Set_aperture(
						std::max(aperture, faille->Get_aperture()));
				model_fractures_loc.geo_cont_.tab_fractures_set_[j]->tab_fractures_[indice_faille]->Set_SlopeAzimut(
						slope_azimut);
				model_fractures_loc.geo_cont_.tab_fractures_set_[j]->tab_fractures_[indice_faille]->Set_RatioPermeability(
						ratio_perm);

			}

		}

	}

}
//-------------------------------------------
// generation of the fractures
//-------------------------------------------
/**
 * @fn void fractures_set_generation(const geode::StructuralModel&, bool)
 * @brief
 *
 * @param model
 * @param estim_volume
 */
void RemFracResSim::fractures_set_generation(const geode::StructuralModel &model,bool estim_volume) {

	nb_fracture_sets_ = fracturesSet_list_.size();

	set_model_fractures_bounding_box(model);
	nb_failles_ = 0;

	initialize_beddings_discontinuity(model);
	initalize_mechanical_discontinuity(model);

	std::cout << " --------------  SINGLE FAULT CREATION ------------------"
			<< endl;
	int nb_success = 0;
	for (FractureSet &setFrac : fracturesSet_list_) {

		setFrac.fracture_set_index_ += nb_failles_;
		setFrac.set_volume(volume_);
		setFrac.model_fractures_.geo_cont_ = model_fractures_.geo_cont_;
		setFrac.set_area(area_);
		setFrac.set_generation_box(generation_box_);
		std::shared_ptr<fault_array> test_array = NULL;
		setFrac.fractures_generation_from_shadow_zone(model,  tree_,disteval_,cells_,test_array,estim_volume);

		model_fractures_.geo_cont_.add_fracture_set(setFrac.stress_zone_set_);

		if (setFrac.is_shadow_zone_active_) {
			setFrac.type_intersection_ =
					fractures_intersect::FracResIntersectionTypeEnum::CROSSING;
			cout << " Fracture set " << setFrac.name_
					<< " is crossing with shadow zone" << endl;

		} else {
			if (setFrac.type_intersection_
					== fractures_intersect::FracResIntersectionTypeEnum::CROSSING) {
				cout << " Fracture set " << setFrac.name_
						<< " is crossing with no shadow zone" << endl;
			} else if (setFrac.type_intersection_
					== fractures_intersect::FracResIntersectionTypeEnum::ABUTTING) {
				cout << " Fracture set " << setFrac.name_
						<< " is abbuting with no shadow zone" << endl;
			} else {
				cout << " Fracture set " << setFrac.name_
						<< " is randomly crossing and abbuting with no shadow zone"
						<< endl;
			}
		}

		setFrac.fractures_activable_ = setFrac.fractures_number_;

		setFrac.fractures_inserted_ = 0;

		if (setFrac.success_) {
			nb_success++;

		}

	}

	std::cout
			<< " --------------  FAULT ARRAY CREATION IN PROGRESS ------------------"
			<< endl;
	for (FaultArraySet &setArrayZone : faultArraySet_list_) {

		setArrayZone.fault_array_set_index_ += nb_failles_
				+ fracturesSet_list_.size();
		setArrayZone.volume_ = volume_;
		setArrayZone.fault_array_geo_set_ =
				std::shared_ptr < fractures_intersect::fault_arrays_set_geo
						> (new fractures_intersect::fault_arrays_set_geo(
								setArrayZone.fault_array_set_index_,
								setArrayZone.seed_));
		setArrayZone.area_ = area_;
		setArrayZone.set_generation_box(generation_box_);
		setArrayZone.fault_array_generation(model,  tree_,disteval_,cells_,model_fractures_.geo_cont_,estim_volume);

		if (setArrayZone.success_) {
			nb_success++;

		}

	}

	if (nb_success == (nb_fracture_sets_ + faultArraySet_list_.size()))
		success_ = true;

}

/**
 * @fn absl::flat_hash_map<geode::uuid,std::vector<geode::uuid>> map_surface_from_horizons(const geode::StructuralModel&)
 * @brief Identify model surface and horizons object
 *
 * @pre
 * @post
 * @param model
 * @return
 */

absl::flat_hash_map<geode::uuid, std::vector<geode::uuid> > RemFracResSim::map_beddings_from_horizons(
		const geode::StructuralModel &model) {
	absl::flat_hash_map<geode::uuid, std::vector<geode::uuid> > map_surface_to_horizons;

	index_t nbhorizons = model.nb_horizons();
	if (nbhorizons > 0) {
		std::vector<geode::uuid> horizons_uuids;
		horizons_uuids.reserve(model.nb_horizons());
		for (const auto &horizon : model.horizons()) {
			horizons_uuids.push_back(horizon.id());
			std::vector<geode::uuid> item_uuids;
			std::vector<geode::uuid> item_bounds;
			std::vector<geode::uuid> item_horizons;
			std::string name_horizon = (std::string) horizon.name();

			if (beddings_name_.find(name_horizon) != beddings_name_.end()) {
				for (const auto &item : model.horizon_items(horizon)) {
					item_uuids.push_back(item.id());
					if (model.nb_collections(item.id()) == 1) {
						item_horizons.push_back(item.id());
						discontinuity_map_to_surface_id_is_active_.insert( {
								item.id(),
								beddings_surface_is_active_[name_horizon] });

					} else if (model.nb_collections(item.id()) > 1) {

						item_bounds.push_back(item.id());
					}
				}
				map_surface_to_horizons.insert( { horizon.id(), item_uuids });

				if (item_horizons.size() > 0) {
					discontinuity_map_to_surface_id_.insert( { name_horizon,
							item_horizons });

				}

				if (item_bounds.size() > 0)
					border_map_to_surface_id_.insert( { name_horizon,
							item_bounds });
			}
		}
	}
	return map_surface_to_horizons;
}

/**
 * @fn absl::flat_hash_map<geode::uuid,std::vector<geode::uuid>> map_surface_from_horizons(const geode::StructuralModel&)
 * @brief Identify model surface and horizons object
 *
 * @pre
 * @post
 * @param model
 * @return
 */

absl::flat_hash_map<geode::uuid, std::vector<geode::uuid> > RemFracResSim::map_mechanical_from_horizons(
		const geode::StructuralModel &model) {
	absl::flat_hash_map<geode::uuid, std::vector<geode::uuid> > map_mechanical_to_horizons;

	index_t nbhorizons = model.nb_horizons();
	if (nbhorizons > 0) {
		std::vector<geode::uuid> horizons_uuids;
		horizons_uuids.reserve(model.nb_horizons());
		for (const auto &horizon : model.horizons()) {
			horizons_uuids.push_back(horizon.id());
			std::vector<geode::uuid> item_uuids;
			std::vector<geode::uuid> item_bounds;
			std::vector<geode::uuid> item_horizons;
			std::string name_horizon = (std::string) horizon.name();

			if (mechanical_name_.find(name_horizon) != mechanical_name_.end()) {
				for (const auto &item : model.horizon_items(horizon)) {
					item_uuids.push_back(item.id());
					if (model.nb_collections(item.id()) == 1) {
						item_horizons.push_back(item.id());
						mechanical_discontinuity_map_to_surface_id_is_active_.insert(
								{ item.id(),
										mechanical_surface_is_active_[name_horizon] });

					} else if (model.nb_collections(item.id()) > 1) {

						item_bounds.push_back(item.id());
					}
				}
				map_mechanical_to_horizons.insert(
						{ horizon.id(), item_uuids });

				if (item_horizons.size() > 0) {
					mechanical_discontinuity_map_to_surface_id_.insert( {
							name_horizon, item_horizons });

				}

				if (item_bounds.size() > 0)
					mechanical_border_map_to_surface_id_.insert( { name_horizon,
							item_bounds });
			}
		}
	}
	return map_mechanical_to_horizons;
}

/**
 * @fn void simulate_intersection(bool)
 * @brief
 *
 * @param active
 */

void RemFracResSim::simulate_intersection(bool active) {

	//fractures_intersect::Point3D pt(generation_box_.center().value(0),generation_box_.center().value(1),generation_box_.center().value(2));

	//Cfaille* faille = (Cfaille*) model_fractures_.geo_cont_.tab_plans[0];

	//bool ducon =  model_fractures_.geo_cont_.in_all_stress_zone(pt);

	// ducon =  model_fractures_.geo_cont_.in_all_stress_zone_3d(faille);

	model_fractures_.geo_cont_.intersect_all_fracture_set(10);

}

//---------------------------------------------------
// estimate the fractures density for the karstsim sgrid
//---------------------------------------------------
/**
 * @fn void initialize_volume_and_area_from_global_grid_of_structural_model(const geode::StructuralModel&)
 * @brief
 *
 * @param model
 */
void RemFracResSim::initialize_volume_and_area_from_global_grid_of_structural_model(
		const geode::StructuralModel &model) {
	// get the total grid from geode structural model
	const auto &solidTot = geode::convert_brep_into_solid(model);
	const auto att_value =
			solidTot->polyhedron_attribute_manager().find_attribute
					< geode::uuid_from_conversion_attribute_type
					> (geode::uuid_from_conversion_attribute_name);
	const auto &box = solidTot->bounding_box();
	double lengthX = box.max().value(0) - box.min().value(0);
	double lengthY = box.max().value(1) - box.min().value(1);
	double lengthZ = box.max().value(2) - box.min().value(2);

	area_ = lengthX * lengthY;
	volume_ = area_ * lengthZ;

}

//---------------------------------------------------
//   update outpout directory path for the FracResSim
//---------------------------------------------------
/**
 * @fn void update_fracture_set_outpout_directory_path(std::string)
 * @brief
 *
 * @param path
 */
void RemFracResSim::update_fracture_set_outpout_directory_path(std::string path) {
	for (FractureSet &setFrac : fracturesSet_list_) {
		setFrac.directory_path_ = path;
	}
	for (FaultArraySet &setArray : faultArraySet_list_) {
		setArray.directory_path_ = path;
	}
}

/**
 * @fn std::string enumCategoryToString(fractures_intersect::FracResIntersectionCategoryEnum)
 * @brief
 *
 * @param cat
 * @return
 */

std::string RemFracResSim::enumCategoryToString(
		fractures_intersect::FracResIntersectionCategoryEnum cat) {

	std::string buffer = "BEDDINGS";
	switch (cat) {
	case fractures_intersect::FracResIntersectionCategoryEnum::MECHANIAL_UNITS:
		buffer = "MECHANICAL_UNITS";
		break;
	case fractures_intersect::FracResIntersectionCategoryEnum::BED_PARALLEL:
		buffer = "BED_PARALLEL";
		break;
	case fractures_intersect::FracResIntersectionCategoryEnum::BED_INTERFACE:
		buffer = "BEDDINGS";
		break;
	case fractures_intersect::FracResIntersectionCategoryEnum::FAULTS:
		buffer = "SingleFault";
		break;
	case fractures_intersect::FracResIntersectionCategoryEnum::FAULT_ARRAYS:
		buffer = "ArrayFault";
		break;
	case fractures_intersect::FracResIntersectionCategoryEnum::JOINTS:
		buffer = "JOINTS";
		break;
	case fractures_intersect::FracResIntersectionCategoryEnum::FRACTURE_CLUSTERS_CORRIDORS_ZONES:
		buffer = "CorridorFault";
		break;
	case fractures_intersect::FracResIntersectionCategoryEnum::FRACTURE_CLUSTERS_FAULT_ZONES:
		buffer = "FaultZone";
		break;
	case fractures_intersect::FracResIntersectionCategoryEnum::FRACTURE_CLUSTERS_JOINT_SWARM:
		buffer = "FRACTURE_CLUSTERS_JOINT_SWARM";
		break;
	default:
		break;
	}
	return buffer;
}
/**
 * @fn std::string enumIntersectionToString(fractures_intersect::FracResIntersectionTypeEnum)
 * @brief
 *
 * @param type
 * @return
 */
std::string RemFracResSim::enumIntersectionToString(
		fractures_intersect::FracResIntersectionTypeEnum type) {

	std::string buffer = "CROSSING";
	switch (type) {
	case fractures_intersect::FracResIntersectionTypeEnum::CROSSING:
		buffer = "CROSSING";
		break;
	case fractures_intersect::FracResIntersectionTypeEnum::ABUTTING:
		buffer = "ABUTTING";
		break;
	case fractures_intersect::FracResIntersectionTypeEnum::RANDOM:
		buffer = "RANDOM";
		break;
	default:
		break;
	}
	return buffer;
}
/**
 * @fn void saveFractureSet(std::string, const std::string&, FractureSet*, int)
 * @brief
 *
 * @param methode_suffixe
 * @param surface_filename2
 * @param setFrac
 * @param index
 */
void RemFracResSim::saveFractureSet(fractures_intersect::model& model,const std::string &surface_filename2, FractureSet *setFrac) {
	std::vector<fractures_intersect::Point3D> tab_point;
	std::vector<std::vector<int> > tab_triangle;
	std::vector<fractures_intersect::property<float> > tab_triangle_face_property;
	model.geo_cont_.get_fracture_set_surface_mesh(setFrac->model_geo_cont_index_, tab_point,
			tab_triangle, tab_triangle_face_property);
	setFrac->clean_vertices_and_triangles_array();
	for (int j = 0; j < tab_point.size(); j++) {
		geode::Point3D pt;
		pt.set_value(0, tab_point[j].x());
		pt.set_value(1, tab_point[j].y());
		pt.set_value(2, tab_point[j].z());
		setFrac->vertices_.push_back(pt);
	}
	for (int k = 0; k < tab_triangle.size(); k++) {
		std::array<geode::index_t, 3> tri;
		tri[0] = tab_triangle[k][0];
		tri[1] = tab_triangle[k][1];
		tri[2] = tab_triangle[k][2];
		setFrac->triangles_.push_back(tri);
		int ij = k / 2;
	}
	for (int l = 0; l < tab_triangle_face_property.size(); l++) {
		if (tab_triangle_face_property[l].size() > 6) {
			setFrac->triangles_l1_.push_back(
					tab_triangle_face_property[l].value(4));
			setFrac->triangles_l2_.push_back(
					tab_triangle_face_property[l].value(5));
			setFrac->triangles_orientation_.push_back(
					tab_triangle_face_property[l].value(2));
			setFrac->triangles_dip_.push_back(
					tab_triangle_face_property[l].value(3));
			setFrac->triangles_aperture_.push_back(
					tab_triangle_face_property[l].value(6));
		}
	}
	/*std::string methode_suffixe = absl::StrCat(
			absl::StrCat("_", enumCategoryToString(setFrac->type_category_),
					"_"), enumIntersectionToString(setFrac->type_intersection_),
			"_");
	std::string surface_filename = absl::StrCat(surface_filename2,
			methode_suffixe);*/
	setFrac->save_fracture_set(surface_filename2);
	std::string methode_suffixe = absl::StrCat("_", enumCategoryToString(setFrac->type_category_),
						"_");

	std::string fracname =  absl::StrCat(setFrac->name_, methode_suffixe, setFrac->fracture_set_index_);

	std::string surface_final = absl::StrCat(surface_filename2,  fracname);
	const auto output_file_native = absl::StrCat(setFrac->directory_path_,
			surface_final, ".vtk");
	fractures_intersect::save_meshed_vtk_plan(output_file_native, tab_point,
			tab_triangle, tab_triangle_face_property);
	auto surface = geode::TriangulatedSurface3D::create(
			geode::OpenGeodeTriangulatedSurface3D::impl_name_static());
	auto builder = geode::TriangulatedSurfaceBuilder3D::create(*surface);
	setFrac->init_surface_geometry(setFrac->fracture_set_index_, *surface, *builder);
	const auto output_file_native_ts = absl::StrCat(setFrac->directory_path_,
			surface_final, ".vtp");
	geode::save_triangulated_surface(*surface, output_file_native_ts);
	if (setFrac->is_shadow_zone_active_) {
		tab_point.clear();
		tab_triangle.clear();
		tab_triangle_face_property.clear();
		std::string surface_final = absl::StrCat(
				absl::StrCat(surface_filename2, "_cyl_"), fracname);
		const auto output_file_cycl_native = absl::StrCat(
				setFrac->directory_path_, surface_final, ".vtk");
		setFrac->stress_zone_set_->get_fracture_set_protection_zone_surface_mesh(
				tab_point, tab_triangle, tab_triangle_face_property);
		fractures_intersect::save_meshed_vtk_plan(output_file_cycl_native,
				tab_point, tab_triangle, tab_triangle_face_property, true);
	}
}
/**
 * @fn void saveArrayFault(std::string, const std::string&, FaultArraySet*, int)
 * @brief
 *
 * @param methode_suffixe
 * @param surface_filename2
 * @param setArrayZone
 * @param index
 */
void RemFracResSim::saveArrayFault(std::string methode_suffixe,
		const std::string &surface_filename2, FaultArraySet *setArrayZone) {
	std::vector<fractures_intersect::Point3D> tab_point;
	std::vector<std::vector<int> > tab_triangle;
	std::vector<fractures_intersect::property<float> > tab_triangle_face_property;


	std::string surface_filename = absl::StrCat(surface_filename2,"_",setArrayZone->name_);

	std::string surface_filename_suffixe = absl::StrCat(surface_filename,"_",methode_suffixe);

	std::string surface_final = absl::StrCat(surface_filename_suffixe,"_",setArrayZone->fault_array_set_index_);

	const auto output_file_native = absl::StrCat(setArrayZone->directory_path_,surface_final, ".vtk");
	setArrayZone->fault_array_geo_set_->get_array_surface_mesh(tab_point,tab_triangle, tab_triangle_face_property);
	fractures_intersect::save_meshed_vtk_plan(output_file_native, tab_point,tab_triangle, tab_triangle_face_property, true);
	cout << " Fichier Array output " << output_file_native << endl;
	tab_point.clear();
	tab_triangle.clear();
	tab_triangle_face_property.clear();
	methode_suffixe += "FRACTURE_SET_";
	surface_filename = absl::StrCat(surface_filename2, methode_suffixe);
	surface_final = absl::StrCat(surface_filename,setArrayZone->fracturesSet_list_[0].name_);
	std::string output_file_native_frac = absl::StrCat(setArrayZone->directory_path_, surface_final, ".vtk");
	setArrayZone->fault_array_geo_set_->get_fracture_set_surface_mesh(tab_point,tab_triangle, tab_triangle_face_property);
	fractures_intersect::save_meshed_vtk_plan(output_file_native_frac,tab_point, tab_triangle, tab_triangle_face_property, true);
	cout << " Fichier fracture set output " << output_file_native_frac << endl;
	if (setArrayZone->is_shadow_zone_active_) {
		tab_point.clear();
		tab_triangle.clear();
		tab_triangle_face_property.clear();
		surface_filename = absl::StrCat(surface_filename2,"_",setArrayZone->name_);

		std::string surface_final_suffixe = absl::StrCat(surface_filename,"_",methode_suffixe);

		std::string surface_final = absl::StrCat(surface_final_suffixe,"_",setArrayZone->fault_array_set_index_);

		surface_final = absl::StrCat(surface_final, "_cyl");
		const auto output_file_cycl_native = absl::StrCat(setArrayZone->directory_path_, surface_final, ".vtk");
		setArrayZone->fault_array_geo_set_->get_fracture_set_protection_zone_surface_mesh(tab_point, tab_triangle, tab_triangle_face_property);
		fractures_intersect::save_meshed_vtk_plan(output_file_cycl_native,tab_point, tab_triangle, tab_triangle_face_property, true);
		cout << " Fichier fracture set output with shadow zone "<< output_file_cycl_native << endl;
	}
}





void RemFracResSim::save_fracture_surface_mesh_from_model(int ifaille,fractures_intersect::model& model,std::string output_file_native) {
	std::vector<fractures_intersect::Point3D> tab_point;
	std::vector<std::vector<int> > tab_triangle;
	std::vector<fractures_intersect::property<float> > tab_triangle_face_property;
	model.geo_cont_.get_fracture_set_surface_mesh(ifaille, tab_point,
			tab_triangle, tab_triangle_face_property);
	save_meshed_vtk_plan(output_file_native, tab_point, tab_triangle,
			tab_triangle_face_property);
}

//---------------------------------------------------
//   save  the fracture set data for the FracResSim
//---------------------------------------------------
/**
 * @fn bool saveFromStage(std::string, std::string)
 * @brief
 *
 * @param surface_filename2
 * @param path
 * @return
 */
bool RemFracResSim::saveFromStage(std::string surface_filename2,
		std::string path) {
	std::unordered_map<std::string, std::vector<geode::uuid> >::iterator it =
			discontinuity_map_to_surface_id_.begin();

	int ifaille = 0;
	std::string methode_suffixe;
	for (it = discontinuity_map_to_surface_id_.begin();
			it != discontinuity_map_to_surface_id_.end(); it++) {
		std::vector<fractures_intersect::Point3D> tab_point;
		std::vector<std::vector<int>> tab_triangle;
		std::vector<fractures_intersect::property<float>> tab_triangle_face_property;

		if (beddings_surface_is_active_[it->first]) {
			model_fractures_.geo_cont_.get_fracture_set_surface_mesh(ifaille,
					tab_point, tab_triangle, tab_triangle_face_property);

			std::string bed_type_name = "CROSSING_";

			methode_suffixe = absl::StrCat("_BEDDINGS_", bed_type_name);

			std::string surface_filename = absl::StrCat(surface_filename2,
					methode_suffixe);
			std::string surface_final = absl::StrCat(surface_filename,
					it->first);
			const auto output_file_native = absl::StrCat(path, surface_final,
					".vtk");
			save_meshed_vtk_plan(output_file_native, tab_point, tab_triangle,
					tab_triangle_face_property);
			ifaille++;
		}

	}

	std::unordered_map<std::string, std::vector<geode::uuid> >::iterator itm =
			mechanical_discontinuity_map_to_surface_id_.begin();

	for (itm = mechanical_discontinuity_map_to_surface_id_.begin();
			itm != mechanical_discontinuity_map_to_surface_id_.end(); itm++) {
		std::vector<fractures_intersect::Point3D> tab_point;
		std::vector<std::vector<int>> tab_triangle;
		std::vector<fractures_intersect::property<float>> tab_triangle_face_property;

		if (mechanical_surface_is_active_[itm->first]) {
			model_fractures_.geo_cont_.get_fracture_set_surface_mesh(ifaille,
					tab_point, tab_triangle, tab_triangle_face_property);
			std::string bed_type_name = "CROSSING_";

			methode_suffixe = absl::StrCat("_MECHANICAL_UNITS_", bed_type_name);

			std::string surface_filename3 = absl::StrCat(surface_filename2,
					methode_suffixe);
			std::string surface_final3 = absl::StrCat(surface_filename3,
					itm->first);
			const auto output_file_native = absl::StrCat(path, surface_final3,
					".vtk");
			save_meshed_vtk_plan(output_file_native, tab_point, tab_triangle,
					tab_triangle_face_property);
			ifaille++;
		}

	}
	int nb_stage_ = fracResStageList_.size();

	nb_failles_ = 0;
	int nb_success = 0;
	int nb_total_dfn = 0;
	for (FracResStage stage : fracResStageList_) {

		std::cout << " -------------- Stage N°" << stage.stageIndex_
				<< " named " << stage.name_ << "------------------" << endl;
		int nb_hypothesis = stage.hypotheisisList_.size();
		if (stage.beddingStage_) {

			std::cout
					<< " -------------- bedding : Bed_interface Hypothesis -------------------"
					<< endl;
			std::cout << "number of bedding hypothesis : " << nb_hypothesis
					<< endl;

		} else {
			std::cout
					<< " ...............Bedding : No bed interface ---------------------------"
					<< endl;
			std::cout << "number of  hypothesis : " << nb_hypothesis << endl;
		}
		int indice_bed_interface = 0;
		int indice_bed_parallel = 0;
		int indice_single_fault = 0;
		int indice_fault_array = 0;
		int indice_fault_corridor = 0;
		int indice_fault_zone = 0;

		for (FracResHypothesis fracResHypothesis : stage.hypotheisisList_) {

			for (int type_dfn : fracResHypothesis.type_dfn_) {

				switch (type_dfn) {
				case 0: {
					FracResBeddingAperture *bed_interface =fracResHypothesis.beddingList_[indice_bed_interface];
					if (!bed_interface->is_active_) {
						cout << " bed interface set " << bed_interface->name_<< " is ianctive" << endl;
						break;
					}

					int ifaille = 0;
					for (auto const &stress_geo : stress_geo_list_) {
						std::vector<fractures_intersect::Point3D> tab_point;
						std::vector<std::vector<int>> tab_triangle;
						std::vector<fractures_intersect::property<float>> tab_triangle_face_property;

						model_fractures_.geo_cont_.get_fracture_set_surface_mesh(ifaille, tab_point, tab_triangle,tab_triangle_face_property);

						if (tab_point.size() <= 0) {
							ifaille++;
							continue;
						}
						std::string bed_type_name = "INTERFACE_SLIP_";

						methode_suffixe = absl::StrCat("_BEDDINGS_",bed_type_name);

						std::string surface_filename = absl::StrCat(surface_filename2, methode_suffixe);
						std::pair<int, int> val = stress_geo.first;
						std::string surface_final = absl::StrCat(surface_filename, bed_name_vect_[val]);

						std::string discrete_pair = std::to_string(val.first)+ "_" + std::to_string(val.second);
						surface_final = absl::StrCat(surface_filename,discrete_pair);
						const auto output_file_native = absl::StrCat(path,surface_final, ".vtk");
						save_meshed_vtk_plan(output_file_native, tab_point,tab_triangle, tab_triangle_face_property);
						ifaille++;
					}

					indice_bed_interface++;
					break;
				}
				case 1: {
					FractureSet *setFrac =fracResHypothesis.fractureList_[indice_single_fault];
					if (!setFrac->is_active_) {
						cout << " single fault set " << setFrac->name_
								<< " is ianctive" << endl;
						break;
					}
					std::cout
							<< " --------------  SINGLE FAULT save ------------------"
							<< endl;
					saveFractureSet(model_fractures_,surface_filename2, setFrac);
					indice_single_fault++;
					break;
				}
				case 2: {
					FaultArraySet *setArrayZone =
							fracResHypothesis.faultArrayList_[indice_fault_array];
					if (!setArrayZone->is_active_) {
						cout << " fault array set " << setArrayZone->name_
								<< " is ianctive" << endl;
						break;
					}
					std::cout
							<< " --------------  FAULT ARRAY save ------------------"
							<< endl;

					methode_suffixe = "_FAULTS_ARRAY_";
					saveArrayFault(methode_suffixe, surface_filename2,
							setArrayZone);
					indice_fault_array++;
					break;
				}
				case 3: {
					FracResBeddingProperties *bed_parallel =fracResHypothesis.bed_parallel_List_[indice_bed_parallel];

					int begin = bed_parallel->index_fracture_set_geo_begin_;
					int end = bed_parallel->index_fracture_set_geo_begin_+ bed_parallel->nb_fracture_set_geo_;

					for (int ifaille = begin; ifaille < end; ifaille++) {



						std::string bed_type_name = "PARALLEL_SLIP_";
						methode_suffixe = absl::StrCat("_BEDDINGS_",bed_type_name);
						std::string surface_filename = absl::StrCat(surface_filename2, methode_suffixe);
						std::string surface_final = absl::StrCat(surface_filename, bed_parallel->name_);
						const auto output_file_native = absl::StrCat(path,surface_final, ".vtk");
						save_fracture_surface_mesh_from_model(ifaille,model_fractures_,output_file_native);

						ifaille++;

					}

					indice_bed_parallel++;
					break;
				}
				case 4: {
					FaultArraySet *setArrayCorridor =
							fracResHypothesis.faultCorridorList_[indice_fault_corridor];
					if (!setArrayCorridor->is_active_) {
						cout << " cluster fault corridor set "
								<< setArrayCorridor->name_ << " is ianctive"
								<< endl;
						break;
					}
					std::cout
							<< " --------------  FAULT Corridor save ------------------"
							<< endl;

					methode_suffixe = "_CLUSTERS_FAULTS_CORRIDORS_";
					saveArrayFault(methode_suffixe, surface_filename2,
							setArrayCorridor);
					indice_fault_corridor++;
					break;
				}
				case 5: {
					FaultArraySet *setArrayZone =
							fracResHypothesis.faultZoneList_[indice_fault_zone];
					if (!setArrayZone->is_active_) {
						cout << " cluster fault zone set "
								<< setArrayZone->name_ << " is ianctive"
								<< endl;
						break;
					}
					std::cout
							<< " --------------  FAULT Zone save ------------------"
							<< endl;

					methode_suffixe = "_CLUSTERS_FAULTS_ZONES_";
					saveArrayFault(methode_suffixe, surface_filename2,
							setArrayZone);
					indice_fault_zone++;
					break;
				}
				default: {
					return false;
				}

				} // switch

			} // loop type_dfn

		} //loop Hypothesis
	} // loop stage

	return true;

}
bool RemFracResSim::saveMechanicalUnit(std::string surface_filename2,std::string path) {

	std::unordered_map<std::string, std::vector<geode::uuid> >::iterator it =discontinuity_map_to_surface_id_.begin();

	int ifaille = 0;
	std::string methode_suffixe;
	for (it = discontinuity_map_to_surface_id_.begin();it != discontinuity_map_to_surface_id_.end(); it++) {
		std::vector<fractures_intersect::Point3D> tab_point;
		std::vector<std::vector<int>> tab_triangle;
		std::vector<fractures_intersect::property<float>> tab_triangle_face_property;

		if (beddings_surface_is_active_[it->first]) {
			model_fractures_.geo_cont_.get_fracture_set_surface_mesh(ifaille,tab_point, tab_triangle, tab_triangle_face_property);
			std::string bed_type_name = "CROSSING_";
			methode_suffixe = absl::StrCat("_BEDDINGS_", bed_type_name);
			std::string surface_filename = absl::StrCat(surface_filename2,methode_suffixe);
			std::string surface_final = absl::StrCat(surface_filename,it->first);
			const auto output_file_native = absl::StrCat(path, surface_final,".vtk");
			save_meshed_vtk_plan(output_file_native, tab_point, tab_triangle,tab_triangle_face_property);
			ifaille++;
		}

	}

	std::unordered_map<std::string, std::vector<geode::uuid> >::iterator itm =mechanical_discontinuity_map_to_surface_id_.begin();

	for (itm = mechanical_discontinuity_map_to_surface_id_.begin();itm != mechanical_discontinuity_map_to_surface_id_.end(); itm++) {
		std::vector<fractures_intersect::Point3D> tab_point;
		std::vector<std::vector<int>> tab_triangle;
		std::vector<fractures_intersect::property<float>> tab_triangle_face_property;

		if (mechanical_surface_is_active_[itm->first]) {
			model_fractures_.geo_cont_.get_fracture_set_surface_mesh(ifaille,tab_point, tab_triangle, tab_triangle_face_property);
			std::string bed_type_name = "CROSSING_";
			methode_suffixe = absl::StrCat("_MECHANICAL_UNITS_", bed_type_name);
			std::string surface_filename3 = absl::StrCat(surface_filename2,methode_suffixe);
			std::string surface_final3 = absl::StrCat(surface_filename3,itm->first);
			const auto output_file_native = absl::StrCat(path, surface_final3,".vtk");
			save_meshed_vtk_plan(output_file_native, tab_point, tab_triangle,tab_triangle_face_property);
			ifaille++;
		}

	}

	return true;

}
//---------------------------------------------------
//   save  the fracture set data for the FracResSim
//---------------------------------------------------
/**
 * @fn bool save(std::string, std::string)
 * @brief
 *
 * @param surface_filename2
 * @param path
 * @return
 */
bool RemFracResSim::save(std::string surface_filename2, std::string path) {

	std::unordered_map<std::string, std::vector<geode::uuid> >::iterator it =
			discontinuity_map_to_surface_id_.begin();

	int ifaille = 0;
	std::string methode_suffixe;
	for (it = discontinuity_map_to_surface_id_.begin();
			it != discontinuity_map_to_surface_id_.end(); it++) {
		std::vector<fractures_intersect::Point3D> tab_point;
		std::vector<std::vector<int>> tab_triangle;
		std::vector<fractures_intersect::property<float>> tab_triangle_face_property;

		if (beddings_surface_is_active_[it->first]) {
			model_fractures_.geo_cont_.get_fracture_set_surface_mesh(ifaille,
					tab_point, tab_triangle, tab_triangle_face_property);

			std::string bed_type_name = "CROSSING_";

			methode_suffixe = absl::StrCat("_BEDDINGS_", bed_type_name);

			std::string surface_filename = absl::StrCat(surface_filename2,
					methode_suffixe);
			std::string surface_final = absl::StrCat(surface_filename,
					it->first);
			const auto output_file_native = absl::StrCat(path, surface_final,
					".vtk");
			save_meshed_vtk_plan(output_file_native, tab_point, tab_triangle,
					tab_triangle_face_property);
			ifaille++;
		}

	}

	std::unordered_map<std::string, std::vector<geode::uuid> >::iterator itm =
			mechanical_discontinuity_map_to_surface_id_.begin();

	for (itm = mechanical_discontinuity_map_to_surface_id_.begin();
			itm != mechanical_discontinuity_map_to_surface_id_.end(); itm++) {
		std::vector<fractures_intersect::Point3D> tab_point;
		std::vector<std::vector<int>> tab_triangle;
		std::vector<fractures_intersect::property<float>> tab_triangle_face_property;

		if (mechanical_surface_is_active_[itm->first]) {
			model_fractures_.geo_cont_.get_fracture_set_surface_mesh(ifaille,
					tab_point, tab_triangle, tab_triangle_face_property);
			std::string bed_type_name = "CROSSING_";

			methode_suffixe = absl::StrCat("_MECHANICAL_UNITS_", bed_type_name);

			std::string surface_filename3 = absl::StrCat(surface_filename2,
					methode_suffixe);
			std::string surface_final3 = absl::StrCat(surface_filename3,
					itm->first);
			const auto output_file_native = absl::StrCat(path, surface_final3,
					".vtk");
			save_meshed_vtk_plan(output_file_native, tab_point, tab_triangle,
					tab_triangle_face_property);
			ifaille++;
		}

	}

	for (FractureSet &setFrac : fracturesSet_list_) {

		std::vector<fractures_intersect::Point3D> tab_point;
		std::vector<std::vector<int>> tab_triangle;
		std::vector<fractures_intersect::property<float>> tab_triangle_face_property;

		model_fractures_.geo_cont_.get_fracture_set_surface_mesh(
				setFrac.fracture_set_index_, tab_point, tab_triangle,
				tab_triangle_face_property);

		setFrac.clean_vertices_and_triangles_array();

		for (int j = 0; j < tab_point.size(); j++) {
			geode::Point3D pt;
			pt.set_value(0, tab_point[j].x());
			pt.set_value(1, tab_point[j].y());
			pt.set_value(2, tab_point[j].z());
			setFrac.vertices_.push_back(pt);
		}

		for (int k = 0; k < tab_triangle.size(); k++) {
			std::array<geode::index_t, 3> tri;
			tri[0] = tab_triangle[k][0];
			tri[1] = tab_triangle[k][1];
			tri[2] = tab_triangle[k][2];

			setFrac.triangles_.push_back(tri);
			int ij = k / 2;

		}
		for (int l = 0; l < tab_triangle_face_property.size(); l++) {
			if (tab_triangle_face_property[l].size() > 6) {
				setFrac.triangles_l1_.push_back(
						tab_triangle_face_property[l].value(4));
				setFrac.triangles_l2_.push_back(
						tab_triangle_face_property[l].value(5));
				setFrac.triangles_orientation_.push_back(
						tab_triangle_face_property[l].value(2));
				setFrac.triangles_dip_.push_back(
						tab_triangle_face_property[l].value(3));
				setFrac.triangles_aperture_.push_back(
						tab_triangle_face_property[l].value(6));
			}

		}

		methode_suffixe = absl::StrCat(
				absl::StrCat("_", enumCategoryToString(setFrac.type_category_),
						"_"),
				enumIntersectionToString(setFrac.type_intersection_), "_");

		std::string surface_filename = absl::StrCat(surface_filename2,
				methode_suffixe);
		setFrac.save_fracture_set(surface_filename);
		std::string surface_final = absl::StrCat(surface_filename,
				setFrac.name_);

		const auto output_file_native = absl::StrCat(setFrac.directory_path_,
				surface_final, ".vtk");

		fractures_intersect::save_meshed_vtk_plan(output_file_native, tab_point,
				tab_triangle, tab_triangle_face_property);

		auto surface = geode::TriangulatedSurface3D::create(
				geode::OpenGeodeTriangulatedSurface3D::impl_name_static());
		auto builder = geode::TriangulatedSurfaceBuilder3D::create(*surface);
		setFrac.init_surface_geometry(setFrac.fracture_set_index_, *surface,
				*builder);
		const auto output_file_native_ts = absl::StrCat(setFrac.directory_path_,
				surface_final, ".vtp");
		geode::save_triangulated_surface(*surface, output_file_native_ts);

		if (setFrac.is_shadow_zone_active_) {
			tab_point.clear();
			tab_triangle.clear();
			tab_triangle_face_property.clear();

			std::string surface_final = absl::StrCat(
					absl::StrCat(surface_filename, "_cyl_"), setFrac.name_);
			const auto output_file_cycl_native = absl::StrCat(
					setFrac.directory_path_, surface_final, ".vtk");

			setFrac.stress_zone_set_->get_fracture_set_protection_zone_surface_mesh(
					tab_point, tab_triangle, tab_triangle_face_property);

			fractures_intersect::save_meshed_vtk_plan(output_file_cycl_native,
					tab_point, tab_triangle, tab_triangle_face_property, true);
		}
	}

	for (FaultArraySet &setFrac : faultArraySet_list_) {

		std::vector<fractures_intersect::Point3D> tab_point;
		std::vector<std::vector<int>> tab_triangle;
		std::vector<fractures_intersect::property<float>> tab_triangle_face_property;

		methode_suffixe = "_FAULTS_ARRAY_";
		std::string surface_filename = absl::StrCat(surface_filename2,
				methode_suffixe);

		std::string surface_final = absl::StrCat(surface_filename,
				setFrac.name_);
		const auto output_file_native = absl::StrCat(setFrac.directory_path_,
				surface_final, ".vtk");
		setFrac.fault_array_geo_set_->get_array_surface_mesh(tab_point,
				tab_triangle, tab_triangle_face_property);
		fractures_intersect::save_meshed_vtk_plan(output_file_native, tab_point,
				tab_triangle, tab_triangle_face_property, true);
		cout << " Fichier fault Array output " << output_file_native << endl;

		tab_point.clear();
		tab_triangle.clear();
		tab_triangle_face_property.clear();

		methode_suffixe = "_FAULTS_ARRAY_FRACTURE_SET_";
		surface_filename = absl::StrCat(surface_filename2, methode_suffixe);
		std::string fracName = absl::StrCat(setFrac.name_,"_", setFrac.fault_array_set_index_);
		surface_final = absl::StrCat(surface_filename, fracName);
		std::string output_file_native_frac = absl::StrCat(
				setFrac.directory_path_, surface_final, ".vtk");
		setFrac.fault_array_geo_set_->get_fracture_set_surface_mesh(tab_point,
				tab_triangle, tab_triangle_face_property);
		fractures_intersect::save_meshed_vtk_plan(output_file_native_frac,
				tab_point, tab_triangle, tab_triangle_face_property, true);
		cout << " Fichier Array fracture set output " << output_file_native
				<< endl;
		if (setFrac.is_shadow_zone_active_) {
			tab_point.clear();
			tab_triangle.clear();
			tab_triangle_face_property.clear();
			surface_filename = absl::StrCat(surface_filename2, methode_suffixe);
			surface_final = absl::StrCat(absl::StrCat(surface_filename, "cyl_"),
					setFrac.name_);
			const auto output_file_cycl_native = absl::StrCat(
					setFrac.directory_path_, surface_final, ".vtk");

			setFrac.fault_array_geo_set_->get_fracture_set_protection_zone_surface_mesh(
					tab_point, tab_triangle, tab_triangle_face_property);

			fractures_intersect::save_meshed_vtk_plan(output_file_cycl_native,
					tab_point, tab_triangle, tab_triangle_face_property, true);
			cout << " Fichier Array fracture set output "
					<< output_file_cycl_native << endl;
		}

	}
}

void RemFracResSim::initialize_beddings_discontinuity(
		const geode::StructuralModel &model) {
	std::unordered_map<std::string, std::vector<geode::uuid> >::iterator it =
			discontinuity_map_to_surface_id_.begin();

	while (it != discontinuity_map_to_surface_id_.end()) {

		std::cout << it->first << " :: ";
		std::cout << std::endl;
		std::shared_ptr<fractures_intersect::fracture_set_geo> stress_zone_set =
				std::shared_ptr < fractures_intersect::fracture_set_geo
						> (new fractures_intersect::fracture_set_geo(
								nb_failles_,
								fractures_intersect::FracResIntersectionCategoryEnum::BED_INTERFACE,
								seed_, beddings_type_intersection_[it->first]));

		for (auto &uuid : it->second) {

			std::cout << uuid.string() << " :: ";
			std::cout << std::endl;

			if (discontinuity_map_to_surface_id_is_active_[uuid]) {
				const auto &mesh_surf = model.surface(uuid).mesh<
						geode::TriangulatedSurface3D>();
				std::shared_ptr<geode::VariableAttribute<double> > attribut =
						mesh_surf.polygon_attribute_manager().find_or_create_attribute<
								geode::VariableAttribute, double>(
								beddings_name_[it->first], ndvalue_);

				for (const auto &polyId : geode::Range { mesh_surf.nb_polygons() }) {

					const geode::Triangle3D tri = mesh_surf.triangle(polyId);

					geode::Point3D p0 = tri.vertices()[0].get();
					geode::Point3D p1 = tri.vertices()[1].get();
					geode::Point3D p2 = tri.vertices()[2].get();

					fractures_intersect::Point3D germe_pt0 =
							fractures_intersect::Point3D(p0.value(0),
									p0.value(1), p0.value(2));
					fractures_intersect::Point3D germe_pt1 =
							fractures_intersect::Point3D(p1.value(0),
									p1.value(1), p1.value(2));
					fractures_intersect::Point3D germe_pt2 =
							fractures_intersect::Point3D(p2.value(0),
									p2.value(1), p2.value(2));

					std::vector<fractures_intersect::Point3D> nodes;
					nodes.push_back(germe_pt0);
					nodes.push_back(germe_pt1);
					nodes.push_back(germe_pt2);
					double value = beddings_aperture_default_value_[it->first];
					std::shared_ptr<fractures_intersect::Cfaille> faille =
							std::shared_ptr < fractures_intersect::Cfaille
									> (new fractures_intersect::Cfaille(
											nb_failles_, value, nodes,
											it->first));
					stress_zone_set->add_fracture(faille);

				}
			}
		}
		if (beddings_surface_is_active_[it->first]) {
			model_fractures_.geo_cont_.add_fracture_set(stress_zone_set);
			nb_failles_++;
		}
		it++;
	}

}

/**
 * @fn void initalize_mechanical_discontinuity(const geode::StructuralModel&)
 * @brief
 *
 * @param model
 */

void RemFracResSim::initalize_mechanical_discontinuity(
		const geode::StructuralModel &model) {
	std::unordered_map<std::string, std::vector<geode::uuid> >::iterator it =
			mechanical_discontinuity_map_to_surface_id_.begin();

	while (it != mechanical_discontinuity_map_to_surface_id_.end()) {

		std::cout << it->first << " :: ";
		std::cout << std::endl;
		std::shared_ptr<fractures_intersect::fracture_set_geo> stress_zone_set =
				std::shared_ptr < fractures_intersect::fracture_set_geo
						> (new fractures_intersect::fracture_set_geo(
								nb_failles_,
								fractures_intersect::FracResIntersectionCategoryEnum::MECHANIAL_UNITS,
								seed_, mechanical_type_intersection_[it->first]));

		for (auto &uuid : it->second) {

			std::cout << uuid.string() << " :: ";
			std::cout << std::endl;

			if (mechanical_discontinuity_map_to_surface_id_is_active_[uuid]) {
				const auto &mesh_surf = model.surface(uuid).mesh<
						geode::TriangulatedSurface3D>();
				std::shared_ptr<geode::VariableAttribute<double> > attribut =
						mesh_surf.polygon_attribute_manager().find_or_create_attribute<
								geode::VariableAttribute, double>(
								mechanical_name_[it->first], ndvalue_);

				for (const auto &polyId : geode::Range { mesh_surf.nb_polygons() }) {

					const geode::Triangle3D tri = mesh_surf.triangle(polyId);

					geode::Point3D p0 = tri.vertices()[0].get();
					geode::Point3D p1 = tri.vertices()[1].get();
					geode::Point3D p2 = tri.vertices()[2].get();

					fractures_intersect::Point3D germe_pt0 =
							fractures_intersect::Point3D(p0.value(0),
									p0.value(1), p0.value(2));
					fractures_intersect::Point3D germe_pt1 =
							fractures_intersect::Point3D(p1.value(0),
									p1.value(1), p1.value(2));
					fractures_intersect::Point3D germe_pt2 =
							fractures_intersect::Point3D(p2.value(0),
									p2.value(1), p2.value(2));

					std::vector<fractures_intersect::Point3D> nodes;
					nodes.push_back(germe_pt0);
					nodes.push_back(germe_pt1);
					nodes.push_back(germe_pt2);
					double value = mechanical_aperture_default_value_[it->first];
					std::shared_ptr<fractures_intersect::Cfaille> faille =
							std::shared_ptr < fractures_intersect::Cfaille
									> (new fractures_intersect::Cfaille(
											nb_failles_, value, nodes,
											it->first));
					stress_zone_set->add_fracture(faille);

				}
			}
		}
		if (mechanical_surface_is_active_[it->first]) {
			model_fractures_.geo_cont_.add_fracture_set(stress_zone_set);
			nb_failles_++;
		}
		it++;
	}

}

//---------------------------------------------------
//   open a fracture set data for the karstsim 
//---------------------------------------------------
/**
 * @fn bool open(std::string)
 * @brief
 *
 * @param filename
 * @return
 */
bool RemFracResSim::open(std::string filename) {

	std::fstream fin(filename, fin.in);

	if (fin.is_open()) {
		//const char* code ;
		std::string comment;
		std::string separator;
		bool read_ok = false;

		int index;
		int nb_bound;

		std::string name;
		std::string nameloc;

		fracturesSet_list_.clear();
		fin >> comment >> separator >> name_;

		fin >> comment >> separator >> seed_;

		fin >> comment >> separator >> index;

		use_global_outpout_directory_path_ = bool(index);
		fin >> comment >> separator >> outpout_directory_path_;
		fin >> comment >> separator >> cylinder_step_;

		fin >> name >> comment >> nb_bound;

		std::cout << name << " number: " << nb_bound;
		std::cout << std::endl;

		bool is_acif = false;
		double value = 0;
		int type;
		for (int k = 0; k < nb_bound; k++) {

			fin >> name >> comment >> is_acif >> comment >> nameloc >> comment
					>> value;
			beddings_name_.emplace(name, nameloc);
			beddings_surface_is_active_.emplace(name, is_acif);
			beddings_aperture_default_value_.emplace(name, value);
			std::vector<fractures_intersect::FracResIntersectionTypeEnum> beddings_intersection;
			beddings_intersection.resize(7,
					fractures_intersect::FracResIntersectionTypeEnum::RANDOM);
			for (int i = 0; i < 7; i++) {
				fin >> nameloc >> comment >> type;
				beddings_intersection[i] =
						static_cast<fractures_intersect::FracResIntersectionTypeEnum>(type);
			}
			beddings_type_intersection_.emplace(name, beddings_intersection);

		}

		fin >> name >> comment >> nb_bound;
		std::cout << name << " number: " << nb_bound;
		std::cout << std::endl;

		is_acif = false;
		value = 0;

		for (int k = 0; k < nb_bound; k++) {

			fin >> name >> comment >> is_acif >> comment >> nameloc >> comment
					>> value;
			mechanical_name_.emplace(name, nameloc);
			mechanical_surface_is_active_.emplace(name, is_acif);
			mechanical_aperture_default_value_.emplace(name, value);
			std::vector<fractures_intersect::FracResIntersectionTypeEnum> mechanical_intersection;
			mechanical_intersection.resize(7,
					fractures_intersect::FracResIntersectionTypeEnum::RANDOM);
			for (int i = 0; i < 7; i++) {
				fin >> nameloc >> comment >> type;
				mechanical_intersection[i] =
						static_cast<fractures_intersect::FracResIntersectionTypeEnum>(type);
			}
			mechanical_type_intersection_.emplace(name,
					mechanical_intersection);
		}

		fin >> comment >> separator >> nb_fracture_sets_;
		//fin >> comment >> separator >> index;
		//intersect_inside_set_ = bool(index);

		for (int i = 0; i < nb_fracture_sets_; i++) {

			remfracres::FractureSet obj;
			obj.open_loc(fin);
			obj.fracture_set_index_ = i;
			obj.cylinder_step_ = cylinder_step_;
			fracturesSet_list_.push_back(obj);
		}

		nb_fault_arrays_set_ = 0;

		if (fin.eof())
			return read_ok;

		fin >> comment >> separator >> nb_fault_arrays_set_;
		//fin >> comment >> separator >> index;
		//intersect_inside_set_ = bool(index);

		for (int i = 0; i < nb_fault_arrays_set_; i++) {

			remfracres::FaultArraySet obj;
			obj.open_loc(fin);
			obj.fault_array_set_index_ = i;
			faultArraySet_list_.push_back(obj);
		}

		fin.close();
		return read_ok;
	}
	return false;
}

//---------------------------------------------------
//   open a fracture set data for the karstsim
//---------------------------------------------------

/**
 * @fn bool openFracres(std::string)
 * @brief
 *
 * @param filename
 * @return
 */
bool RemFracResSim::openFracres(std::string filename) {

	std::fstream fin(filename, fin.in);

	if (fin.is_open()) {
		//const char* code ;
		bool read_ok = false;

		read_ok = manageLine(fin);

		fin.close();
		return read_ok;
	}
	return false;
}
/**
 * @fn std::vector<std::string> splitLine(std::istringstream&, char)
 * @brief
 *
 * @param my_stream
 * @param delim
 * @return
 */
std::vector<std::string> RemFracResSim::splitLine(std::istringstream &my_stream,
		char delim) {
	std::vector<std::string> tokenVector;
	while (my_stream) {
		std::string substring;
		if (!std::getline(my_stream, substring, delim))
			break;

		tokenVector.push_back(substring);
	}
	return tokenVector;
}
/**
 * @fn bool manageLine(std::fstream&)
 * @brief
 *
 * @param fin
 * @return
 */
bool RemFracResSim::manageLine(std::fstream &fin) {
	// Traverse till stream is valid
	// To store the stream string
	std::string contenu_ligne;
	string token;
	string value;
	string separator;
	std::istringstream my_stream;
	while (std::getline(fin, contenu_ligne)) // tant que l'on peut mettre la ligne dans "contenu"
	{

		my_stream.str(contenu_ligne);
		cout << contenu_ligne << endl;  // on l'affiche

		std::vector<std::string> tokenVector = splitLine(my_stream, ';');

		if (tokenVector[0] == fracresFileKeyWord_->fracres_model_name) {
			name_ = tokenVector[1];
		} else if (tokenVector[0] == fracresFileKeyWord_->data_support_name) {
			//name_ =tokenVector[1];
		} else if (tokenVector[0]== fracresFileKeyWord_->manage_bed_interface) {
			bool b = (tokenVector[1] == "true");
			manage_bed_interface_ = b;

		} else if (tokenVector[0] == fracresFileKeyWord_->unit_begin) {
			manageUnit(fin);
		} else if (tokenVector[0] == fracresFileKeyWord_->stage_begin) {
			FracResStage fracResStage = manageStage(fin);
			fracResStageList_.push_back(fracResStage);
		} else if (tokenVector[0] == fracresFileKeyWord_->scenario_begin) {
			FracResScenario fracResScenario = manageScenario(fin);
			fracResScenarioList_.push_back(fracResScenario);
		} else if (tokenVector[0] == fracresFileKeyWord_->run_param_begin) {
			manageRunParameter(fin);
		} else if (tokenVector[0] == fracresFileKeyWord_->matrix_prop_begin) {
			manageMatrixProp(fin);
		} else if (tokenVector[0] == fracresFileKeyWord_->outputs_begin) {
			manageOutputs(fin);
		}
		my_stream.clear();
	}
	return true;
}
/**
 * @fn void manageUnit(std::fstream&)
 * @brief
 *
 * @param fin
 */
void RemFracResSim::manageUnit(std::fstream &fin) {

	std::string contenu_ligne;
	std::string comment;
	std::string separator;
	std::istringstream my_stream;

	while (getline(fin, contenu_ligne)) // tant que l'on peut mettre la ligne dans "contenu"
	{

		cout << contenu_ligne << endl;  // on l'affiche
		my_stream.str(contenu_ligne);

		std::vector<std::string> tokenVector = splitLine(my_stream, ';');

		if (tokenVector[0] == fracresFileKeyWord_->unit_age) {
			fracResUnit_->unit_age_ = tokenVector[1];
		} else if (tokenVector[0] == fracresFileKeyWord_->unit_duration) {
			fracResUnit_->unit_duration_ = tokenVector[1];
		} else if (tokenVector[0] == fracresFileKeyWord_->unit_distance) {
			fracResUnit_->unit_distance_ = tokenVector[1];
		} else if (tokenVector[0] == fracresFileKeyWord_->unit_altitude) {
			fracResUnit_->unit_altitude_ = tokenVector[1];
		} else if (tokenVector[0]
				== fracresFileKeyWord_->unit_production_duration) {
			fracResUnit_->unit_production_duration_ = tokenVector[1];
		} else if (tokenVector[0]
				== fracresFileKeyWord_->unit_erosion_duration) {
			fracResUnit_->unit_erosion_duration_ = tokenVector[1];
		} else if (tokenVector[0]
				== fracresFileKeyWord_->unit_subsidence_duration) {
			fracResUnit_->unit_subsidence_duration_ = tokenVector[1];

		} else if (tokenVector[0]
				== fracresFileKeyWord_->unit_river_sediment_duration) {
			fracResUnit_->unit_river_sediment_duration_ = tokenVector[1];
		} else if (tokenVector[0]
				== fracresFileKeyWord_->unit_river_sediment_volume) {
			fracResUnit_->unit_river_sediment_volume_ = tokenVector[1];
		} else if (tokenVector[0] == fracresFileKeyWord_->unit_temperature) {
			fracResUnit_->unit_temperature_ = tokenVector[1];
		} else if (tokenVector[0] == fracresFileKeyWord_->unit_slope) {
			fracResUnit_->unit_slope_ = tokenVector[1];
		} else if (tokenVector[0] == fracresFileKeyWord_->unit_end) {

			break;
		}

		my_stream.clear();
	}

}
/**
 * @fn FracResStage manageStage(std::fstream&)
 * @brief
 *
 * @param fin
 * @return
 */
FracResStage RemFracResSim::manageStage(std::fstream &fin) {
	std::string contenu_ligne;
	std::string comment;
	std::string separator;

	FracResStage fracResStage = FracResStage();
	std::istringstream my_stream;
	while (getline(fin, contenu_ligne)) // tant que l'on peut mettre la ligne dans "contenu"
	{

		cout << contenu_ligne << endl;  // on l'affiche
		my_stream.str(contenu_ligne);

		std::vector<std::string> tokenVector = splitLine(my_stream, ';');

		if (tokenVector[0] == fracresFileKeyWord_->stage_name) {
			fracResStage.name_ = tokenVector[1];
		} else if (tokenVector[0] == fracresFileKeyWord_->stage_index) {
			fracResStage.stageIndex_ = std::atoi(tokenVector[1].c_str());
		} else if (tokenVector[0] == fracresFileKeyWord_->stage_begin_stage) {
			fracResStage.beddingStage_ = (tokenVector[1] == "true");
		} else if (tokenVector[0] == fracresFileKeyWord_->fracture_begin) {
			FracResFracture fracresFracture = FracResFracture();
			fracresFracture = manageFracture(fin);
			fracresFracture.stage_name_ = fracResStage.name_;
			fracResStage.fractureList_.insert( { fracresFracture.index_Id_,
					fracresFracture });

		} else if (tokenVector[0] == fracresFileKeyWord_->hyp_begin) {
			FracResHypothesis fracResHypothesis = FracResHypothesis();

			if (fracResStage.beddingStage_) {
				manageBeddingStageHypothesis(fin, fracResHypothesis);
			} else {
				manageClassicStageHypothesis(fin, fracResStage,
						fracResHypothesis);
			}
			fracResHypothesis.stage_name_ = fracResStage.name_;
			fracResStage.hypotheisisList_.push_back(fracResHypothesis);
		} else if (tokenVector[0] == fracresFileKeyWord_->stage_end) {

			break;
		}

		my_stream.clear();
	}
	return fracResStage;
}
/**
 * @fn FracResHypothesis manageRepHypothesis(std::fstream&)
 * @brief
 *
 * @param fin
 * @return
 */
FracResHypothesis RemFracResSim::manageRepHypothesis(std::fstream &fin) {
	std::string contenu_ligne;
	std::string comment;
	std::string separator;
	FracResHypothesis fracResHypothesis = FracResHypothesis();
	std::istringstream my_stream;
	while (getline(fin, contenu_ligne)) // tant que l'on peut mettre la ligne dans "contenu"
	{

		cout << contenu_ligne << endl;  // on l'affiche
		my_stream.str(contenu_ligne);

		std::vector<std::string> tokenVector = splitLine(my_stream, ';');

		if (tokenVector[0] == fracresFileKeyWord_->hyp_name) {
			fracResHypothesis.name_ = tokenVector[1];
		} else if (tokenVector[0] == fracresFileKeyWord_->hyp_index) {
			fracResHypothesis.hypothesisIndex_ = stoi(tokenVector[1]);
		} else if (tokenVector[0] == fracresFileKeyWord_->hyp_bed_is_region) {
			fracResHypothesis.beddingIsRegion_ = (tokenVector[1] == "true");
		} else if (tokenVector[0] == fracresFileKeyWord_->hyp_bed_region_name) {
			fracResHypothesis.beddingRegionName_ = tokenVector[1];
		} else if (tokenVector[0]
				== fracresFileKeyWord_->hyp_bed_prop_discrete_name) {
			fracResHypothesis.beddingPropertyDiscreteName_ = tokenVector[1];
		} else if (tokenVector[0] == fracresFileKeyWord_->hyp_end) {

			break;
		}
		my_stream.clear();
	}
	return fracResHypothesis;
}

/**
 * @fn void manageMatrixProp(std::fstream&)
 * @brief
 *
 * @param fin
 */
void RemFracResSim::manageMatrixProp(std::fstream &fin) {
	std::string contenu_ligne;
	std::string comment;
	std::string separator;

	std::istringstream my_stream;
	while (getline(fin, contenu_ligne)) // tant que l'on peut mettre la ligne dans "contenu"
	{

		cout << contenu_ligne << endl;  // on l'affiche
		my_stream.str(contenu_ligne);

		std::vector<std::string> tokenVector = splitLine(my_stream, ';');

		if (tokenVector[0] == fracresFileKeyWord_->matrix_prop_init_poro_name) {

		} else if (tokenVector[0]
				== fracresFileKeyWord_->matrix_prop_init_poro_selected) {
			bool test = (tokenVector[1] == "true");
		} else if (tokenVector[0]
				== fracresFileKeyWord_->matrix_prop_init_xperm_name) {

		} else if (tokenVector[0]
				== fracresFileKeyWord_->matrix_prop_init_xperm_selected) {
			bool testx = (tokenVector[1] == "true");
		} else if (tokenVector[0]
				== fracresFileKeyWord_->matrix_prop_init_yperm_name) {

		} else if (tokenVector[0]
				== fracresFileKeyWord_->matrix_prop_init_yperm_selected) {
			bool testy = (tokenVector[1] == "true");
		} else if (tokenVector[0]
				== fracresFileKeyWord_->matrix_prop_init_zperm_name) {

		} else if (tokenVector[0]
				== fracresFileKeyWord_->matrix_prop_init_zperm_selected) {
			bool testz = (tokenVector[1] == "true");
		} else if (tokenVector[0] == fracresFileKeyWord_->matrix_prop_end) {

			break;
		}
		my_stream.clear();
	}

}

/**
 * @fn void manageOutputs(std::fstream&)
 * @brief
 *
 * @param fin
 */
void RemFracResSim::manageOutputs(std::fstream &fin) {
	std::string contenu_ligne;
	std::string comment;
	std::string separator;

	std::istringstream my_stream;
	while (getline(fin, contenu_ligne)) // tant que l'on peut mettre la ligne dans "contenu"
	{

		cout << contenu_ligne << endl;  // on l'affiche
		my_stream.str(contenu_ligne);

		std::vector<std::string> tokenVector = splitLine(my_stream, ';');
		//if(tokenVector.size() == 1 && tokenVector[0] != fracresFileKeyWord_->outputs_end) continue;

		if (tokenVector[0] == fracresFileKeyWord_->outputs_density) {

		} else if (tokenVector[0]
				== fracresFileKeyWord_->outputs_density_selected) {
			bool test = (tokenVector[1] == "true");
		} else if (tokenVector[0]
				== fracresFileKeyWord_->outputs_first_second_exchange) {

		} else if (tokenVector[0]
				== fracresFileKeyWord_->outputs_first_second_exchange_selected) {
			bool test2 = (tokenVector[1] == "true");
		} else if (tokenVector[0]
				== fracresFileKeyWord_->outputs_second_third_exchange) {

		} else if (tokenVector[0]
				== fracresFileKeyWord_->outputs_second_third_exchange_selected) {
			bool test3 = (tokenVector[1] == "true");
		} else if (tokenVector[0]
				== fracresFileKeyWord_->outputs_first_third_exchange) {

		} else if (tokenVector[0]
				== fracresFileKeyWord_->outputs_first_third_exchange_selected) {
			bool test4 = (tokenVector[1] == "true");
		} else if (tokenVector[0] == fracresFileKeyWord_->outputs_fractures) {

		} else if (tokenVector[0]
				== fracresFileKeyWord_->outputs_fractures_selected) {
			bool test5 = (tokenVector[1] == "true");
		} else if (tokenVector[0] == fracresFileKeyWord_->outputs_aperture) {

		} else if (tokenVector[0]
				== fracresFileKeyWord_->outputs_aperture_selected) {
			bool test5 = (tokenVector[1] == "true");
		} else if (tokenVector[0] == fracresFileKeyWord_->outputs_length) {

		} else if (tokenVector[0]
				== fracresFileKeyWord_->outputs_length_selected) {
			bool test6 = (tokenVector[1] == "true");
		} else if (tokenVector[0] == fracresFileKeyWord_->outputs_height) {

		} else if (tokenVector[0]
				== fracresFileKeyWord_->outputs_height_selected) {
			bool test7 = (tokenVector[1] == "true");
		} else if (tokenVector[0] == fracresFileKeyWord_->outputs_dip_azimut) {

		} else if (tokenVector[0]
				== fracresFileKeyWord_->outputs_dip_azimut_selected) {
			bool test8 = (tokenVector[1] == "true");
		} else if (tokenVector[0] == fracresFileKeyWord_->outputs_dip_angle) {

		} else if (tokenVector[0]
				== fracresFileKeyWord_->outputs_dip_angle_selected) {
			bool test9 = (tokenVector[1] == "true");
		} else if (tokenVector[0] == fracresFileKeyWord_->outputs_mult_t) {

		} else if (tokenVector[0]
				== fracresFileKeyWord_->outputs_mult_t_selected) {
			bool test10 = (tokenVector[1] == "true");
		} else if (tokenVector[0] == fracresFileKeyWord_->outputs_clusters) {

		} else if (tokenVector[0]
				== fracresFileKeyWord_->outputs_clusters_selected) {
			bool test11 = (tokenVector[1] == "true");
		} else if (tokenVector[0]
				== fracresFileKeyWord_->outputs_clusters_length) {

		} else if (tokenVector[0]
				== fracresFileKeyWord_->outputs_clusters_length_selected) {
			bool test12 = (tokenVector[1] == "true");
		} else if (tokenVector[0]
				== fracresFileKeyWord_->outputs_clusters_width) {

		} else if (tokenVector[0]
				== fracresFileKeyWord_->outputs_clusters_width_selected) {
			bool test13 = (tokenVector[1] == "true");
		} else if (tokenVector[0]
				== fracresFileKeyWord_->outputs_clusters_height) {

		} else if (tokenVector[0]
				== fracresFileKeyWord_->outputs_clusters_height_selected) {
			bool test14 = (tokenVector[1] == "true");
		} else if (tokenVector[0] == fracresFileKeyWord_->outputs_path) {
			outpout_directory_path_ = tokenVector[1];
		} else if (tokenVector[0] == fracresFileKeyWord_->outputs_end) {

			break;
		}
		my_stream.clear();
	}

}

/**
 * @fn void manageScenario(std::fstream&)
 * @brief
 *
 * @param fin
 */
FracResScenario RemFracResSim::manageScenario(std::fstream &fin) {
	std::string contenu_ligne;
	std::string comment;
	std::string separator;
	FracResScenario fracResScenario = FracResScenario();
	std::istringstream my_stream;
	while (getline(fin, contenu_ligne)) // tant que l'on peut mettre la ligne dans "contenu"
	{


		cout << contenu_ligne << endl;  // on l'affiche
		my_stream.str(contenu_ligne);

		std::vector<std::string> tokenVector = splitLine(my_stream, ';');
		if (tokenVector[0] == fracresFileKeyWord_->outputs_density) {

		} else if (tokenVector[0]== fracresFileKeyWord_->scenario_is_modeling) {
			fracResScenario.modeling_ = (tokenVector[1] == "true");
		} else if (tokenVector[0] == fracresFileKeyWord_->scenario_is_qc) {

			bool testQc = (tokenVector[1] == "true");

		} else if (tokenVector[0] == fracresFileKeyWord_->scenario_nb_real) {
			fracResScenario.numberOfRealisation_ = stoi(tokenVector[1]);
		} else if (tokenVector[0]== fracresFileKeyWord_->scenario_rep_hyp_index) {
			int test0 = stoi(tokenVector[1]);
		} else if (tokenVector[0]== fracresFileKeyWord_->scenario_nb_estimate_fractures) {
			fracResScenario.estimate_Fracture_Number_ = stoi(tokenVector[1]);
		} else if (tokenVector[0] == fracresFileKeyWord_->scenario_name) {
			fracResScenario.name_= tokenVector[1];
		} else if (tokenVector[0] == fracresFileKeyWord_->scenario_index) {
			fracResScenario.scenarioIndex_ = stoi(tokenVector[1]);
		} else if (tokenVector[0]== fracresFileKeyWord_->scenario_stage_hyp_link) {

		} else if (tokenVector[0]== fracresFileKeyWord_->scenario_stages_hyp_map) {
			std::vector<std::pair<int, int>> elements;
			split(tokenVector[1], elements);
			fracResScenario.scenario_stage_hypothesis_index_vector_=elements;
		} else if (tokenVector[0] == fracresFileKeyWord_->scenario_end) {

			break;
		}
		my_stream.clear();
	}

	return fracResScenario;
}
/**
 * @fn void manageStageHypScenarioKeyVal(std::fstream&)
 * @brief
 *
 * @param fin
 */
void RemFracResSim::manageStageHypScenarioKeyVal(std::fstream &fin) {
	std::string contenu_ligne;
	std::string comment;
	std::string separator;

	while (getline(fin, contenu_ligne)) // tant que l'on peut mettre la ligne dans "contenu"
	{
		std::istringstream my_stream;

		cout << contenu_ligne << endl;  // on l'affiche
		my_stream.str(contenu_ligne);

		// To store the stream string
		string token;
		string value;

		// Traverse till stream is valid
		my_stream >> token >> separator >> value;
	}

}
/**
 * @fn void manageUpscalingRep(std::fstream&)
 * @brief
 *
 * @param fin
 */
void RemFracResSim::manageUpscalingRep(std::fstream &fin) {
	std::string contenu_ligne;
	std::string comment;
	std::string separator;

	while (getline(fin, contenu_ligne)) // tant que l'on peut mettre la ligne dans "contenu"
	{
		std::istringstream my_stream;

		cout << contenu_ligne << endl;  // on l'affiche
		my_stream.str(contenu_ligne);

		// To store the stream string
		string token;
		string value;

		// Traverse till stream is valid
		my_stream >> token >> separator >> value;
	}

}
/**
 * @fn void manageRepresentationObject(std::fstream&)
 * @brief
 *
 * @param fin
 */
void RemFracResSim::manageRepresentationObject(std::fstream &fin) {
	std::string contenu_ligne;
	std::string comment;
	std::string separator;

	while (getline(fin, contenu_ligne)) // tant que l'on peut mettre la ligne dans "contenu"
	{
		std::istringstream my_stream;

		cout << contenu_ligne << endl;  // on l'affiche
		my_stream.str(contenu_ligne);

		// To store the stream string
		string token;
		string value;

		// Traverse till stream is valid
		my_stream >> token >> separator >> value;
	}

}

/**
 * @fn FracResFracture manageFracture(std::fstream&)
 * @brief
 *
 * @param fin
 * @return
 */
void RemFracResSim::manageRunParameter(std::fstream &fin) {
	std::string contenu_ligne;
	std::string comment;
	std::string separator;
	std::istringstream my_stream;

	while (getline(fin, contenu_ligne)) // tant que l'on peut mettre la ligne dans "contenu"
	{

		cout << contenu_ligne << endl;  // on l'affiche
		my_stream.str(contenu_ligne);
		std::vector<std::string> tokenVector = splitLine(my_stream, ';');
		if (tokenVector.size() == 1
				&& tokenVector[0] != fracresFileKeyWord_->run_param_end)
			continue;

		if (tokenVector[0] == fracresFileKeyWord_->run_param_seed) {
			seed_ = std::strtod(tokenVector[1].c_str(), NULL);
		} else if (tokenVector[0]
				== fracresFileKeyWord_->run_param_save_upscaled) {

		} else if (tokenVector[0] == fracresFileKeyWord_->run_param_save_dfn) {

		} else if (tokenVector[0]
				== fracresFileKeyWord_->run_param_target_mesh_res) {

		} else if (tokenVector[0]
				== fracresFileKeyWord_->run_param_save_density) {

		} else if (tokenVector[0]
				== fracresFileKeyWord_->run_param_denisty_enum) {

		} else if (tokenVector[0] == fracresFileKeyWord_->run_param_end) {

			break;
		}
		my_stream.clear();
	}

}

/**
 * @fn FracResFracture manageFracture(std::fstream&)
 * @brief
 *
 * @param fin
 * @return
 */
FracResFracture RemFracResSim::manageFracture(std::fstream &fin) {
	std::string contenu_ligne;
	std::string comment;
	std::string separator;
	std::istringstream my_stream;
	FracResFracture fracresFracture = FracResFracture();
	while (getline(fin, contenu_ligne)) // tant que l'on peut mettre la ligne dans "contenu"
	{

		cout << contenu_ligne << endl;  // on l'affiche
		my_stream.str(contenu_ligne);
		std::vector<std::string> tokenVector = splitLine(my_stream, ';');
		if (tokenVector.size() == 1
				&& tokenVector[0] != fracresFileKeyWord_->fracture_end)
			continue;
		if (tokenVector[0] == fracresFileKeyWord_->fracture_name) {
			fracresFracture.name_ = tokenVector[1];
		} else if (tokenVector[0] == fracresFileKeyWord_->fracture_id) {
			fracresFracture.index_Id_ = std::atoi(tokenVector[1].c_str());
		} else if (tokenVector[0] == fracresFileKeyWord_->fracture_type) {
			fracresFracture.dataType_ = (FracResTypeEnum) getFractureType(
					tokenVector[1]);
			if (fracresFracture.dataType_ == FracResTypeEnum::FAULT_ZONES) {
				fracresFracture.categoryType_ =
						fractures_intersect::FracResIntersectionCategoryEnum::FRACTURE_CLUSTERS_FAULT_ZONES;
			} else if (fracresFracture.dataType_
					== FracResTypeEnum::FRACTURE_CORRIDORS) {
				fracresFracture.categoryType_ =
						fractures_intersect::FracResIntersectionCategoryEnum::FRACTURE_CLUSTERS_CORRIDORS_ZONES;
			}
		} else if (tokenVector[0] == fracresFileKeyWord_->fracture_category) {
			fracresFracture.categoryType_ =
					(fractures_intersect::FracResIntersectionCategoryEnum) getCategoryType(
							tokenVector[1]);
		} else if (tokenVector[0] == fracresFileKeyWord_->fracture_end) {

			break;
		}
		my_stream.clear();
	}
	return fracresFracture;

}
/**
 * @fn FracResHypothesis manageHypothesis(std::fstream&)
 * @brief
 *
 * @param fin
 * @return
 */
FracResHypothesis RemFracResSim::manageHypothesis(std::fstream &fin) {
	std::string contenu_ligne;
	std::string comment;
	std::string separator;
	FracResHypothesis fracResHypothesis = FracResHypothesis();
	std::istringstream my_stream;
	while (getline(fin, contenu_ligne)) // tant que l'on peut mettre la ligne dans "contenu"
	{

		cout << contenu_ligne << endl;  // on l'affiche
		my_stream.str(contenu_ligne);

		std::vector<std::string> tokenVector = splitLine(my_stream, ';');
		if (tokenVector.size() == 1
				&& tokenVector[0] != fracresFileKeyWord_->hyp_end)
			continue;
		if (tokenVector[0] == fracresFileKeyWord_->hyp_end) {

			break;
		}
		my_stream.clear();
	}
	return fracResHypothesis;
}
/**
 * @fn void manageBeddingStageHypothesis(std::fstream&, FracResHypothesis&)
 * @brief
 *
 * @param fin
 * @param fracResHypothesis
 */
void RemFracResSim::manageBeddingStageHypothesis(std::fstream &fin,
		FracResHypothesis &fracResHypothesis) {
	std::string contenu_ligne;
	std::string comment;
	std::string separator;

	std::istringstream my_stream;
	while (getline(fin, contenu_ligne)) // tant que l'on peut mettre la ligne dans "contenu"
	{

		cout << contenu_ligne << endl;  // on l'affiche
		my_stream.str(contenu_ligne);

		std::vector<std::string> tokenVector = splitLine(my_stream, ';');

		if (tokenVector[0] == fracresFileKeyWord_->hyp_name) {
			fracResHypothesis.name_ = tokenVector[1];
		} else if (tokenVector[0] == fracresFileKeyWord_->hyp_index) {
			fracResHypothesis.hypothesisIndex_ = std::atoi(
					tokenVector[1].c_str());
		} else if (tokenVector[0] == fracresFileKeyWord_->hyp_bed_is_region) {
			fracResHypothesis.beddingIsRegion_ = (tokenVector[1] == "true");
		} else if (tokenVector[0] == fracresFileKeyWord_->hyp_bed_region_name) {
			std::string property_name =tokenVector[1];
			 std::transform(property_name.begin(), property_name.end(), property_name.begin(), my_toupper());

			fracResHypothesis.beddingRegionName_ = property_name;
		} else if (tokenVector[0]
				== fracresFileKeyWord_->hyp_bed_prop_discrete_name) {
			std::string property_name =tokenVector[1];
			 std::transform(property_name.begin(), property_name.end(), property_name.begin(), my_toupper());
			fracResHypothesis.beddingPropertyDiscreteName_ = property_name;
		} else if (tokenVector[0] == fracresFileKeyWord_->bed_aperture_begin) {
			FracResBeddingAperture *fracResBeddingAperture_ =
					manageBeddingAperture(fin);
			fracResHypothesis.type_dfn_.push_back(0);
			fracResHypothesis.beddingList_.push_back(fracResBeddingAperture_);
		} else if (tokenVector[0] == fracresFileKeyWord_->hyp_end) {

			break;
		}
		my_stream.clear();
	}

}
/**
 * @fn FracResBeddingAperture manageBeddingAperture*(std::fstream&)
 * @brief
 *
 * @param fin
 * @return
 */
FracResBeddingAperture* RemFracResSim::manageBeddingAperture(
		std::fstream &fin) {
	std::string contenu_ligne;
	std::string comment;
	std::string separator;
	std::istringstream my_stream;

	FracResBeddingAperture *fracResBeddingAperture =
			new FracResBeddingAperture();
	while (getline(fin, contenu_ligne)) // tant que l'on peut mettre la ligne dans "contenu"
	{

		cout << contenu_ligne << endl;  // on l'affiche
		my_stream.str(contenu_ligne);

		std::vector<std::string> tokenVector = splitLine(my_stream, ';');

		if (tokenVector[0] == fracresFileKeyWord_->bed_aperture_name) {
			fracResBeddingAperture->name_ = tokenVector[1];
		} else if (tokenVector[0]
				== fracresFileKeyWord_->bed_aperture_parameter) {
			fracResBeddingAperture->aperture_name_ = tokenVector[1];
		} else if (tokenVector[0] == fracresFileKeyWord_->bed_aperture_unit) {
			fracResBeddingAperture->unit_name_ = tokenVector[1];
		} else if (tokenVector[0]
				== fracresFileKeyWord_->bed_aperture_data_type) {
			fracResBeddingAperture->dataType_ =
					(FracResDataTypeEnum) getDataType(tokenVector[1]);
		} else if (tokenVector[0] == fracresFileKeyWord_->bed_aperture_value) {
			fracResBeddingAperture->value_ = std::strtod(tokenVector[1].c_str(),
					NULL);
		} else if (tokenVector[0]
				== fracresFileKeyWord_->bed_aperture_distribution) {
			fracResBeddingAperture->distribution_name_ = tokenVector[1];
		} else if (tokenVector[0]
				== fracresFileKeyWord_->bed_aperture_property_name) {
			std::string property_name =tokenVector[1];
			 std::transform(property_name.begin(), property_name.end(), property_name.begin(), my_toupper());

			fracResBeddingAperture->discrete_property_name_ = property_name;
		} else if (tokenVector[0]
				== fracresFileKeyWord_->bed_aperture_property_name) {
			std::string property_name =tokenVector[1];
			 std::transform(property_name.begin(), property_name.end(), property_name.begin(), my_toupper());

			fracResBeddingAperture->discrete_property_name_ = property_name;
		} else if (tokenVector[0]
				== fracresFileKeyWord_->bed_aperture_discrete_prop_values) {
			std::vector<std::pair<int, int>> elements;
			split(tokenVector[1], elements);
			fracResBeddingAperture->discrete_property_selected_values_ =
					elements;
		} else if (tokenVector[0] == fracresFileKeyWord_->bed_aperture_end) {

			break;
		}
		my_stream.clear();
	}
	return fracResBeddingAperture;

}
/**
 * @fn void split(const string&, std::vector<std::pair<int,int>>&)
 * @brief
 *
 * @param chaine
 * @param elements
 */

void RemFracResSim::split(const string &chaine,
		std::vector<std::pair<int, int> > &elements) {

	std::string sousChaine;

	std::vector<std::string> tokenVector;
	std::istringstream my_stream;
	my_stream.str(chaine);
	while (my_stream) {
		std::string substring;
		if (!std::getline(my_stream, substring, ','))
			break;

		tokenVector.push_back(substring);
	}

	for (int j = 0; j < tokenVector.size(); j++) {
		sousChaine = tokenVector[j];
		size_t pos, begin = 0;
		vector<int> list;
		std::pair<int, int> listPair;
		do {
			pos = sousChaine.find('-', begin);
			int index = std::atoi(
					sousChaine.substr(begin, pos - begin).c_str());
			list.push_back(index);
			begin = pos + 1;
		} while (pos != string::npos);
		listPair = { list[0], list[1] };
		elements.push_back(listPair);
	}

}
/**
 * @fn void manageClassicStageHypothesis(std::fstream&, FracResStage&, FracResHypothesis&)
 * @brief
 *
 * @param fin
 * @param fracresStage
 * @param fracResHypothesis
 */
void RemFracResSim::manageClassicStageHypothesis(std::fstream &fin,
		FracResStage &fracresStage, FracResHypothesis &fracResHypothesis) {
	std::string contenu_ligne;
	std::string comment;
	std::string separator;

	std::istringstream my_stream;
	while (getline(fin, contenu_ligne)) // tant que l'on peut mettre la ligne dans "contenu"
	{

		cout << contenu_ligne << endl;  // on l'affiche
		my_stream.str(contenu_ligne);
		std::vector<std::string> tokenVector = splitLine(my_stream, ';');

		if (tokenVector[0] == fracresFileKeyWord_->hyp_name) {
			fracResHypothesis.name_ = tokenVector[1];
		} else if (tokenVector[0] == fracresFileKeyWord_->hyp_index) {
			fracResHypothesis.hypothesisIndex_ = std::atoi(
					tokenVector[1].c_str());
		} else if (tokenVector[0] == fracresFileKeyWord_->frac_prop_begin) {
			manageFractureProperty(fin, fracresStage, fracResHypothesis);
		} else if (tokenVector[0] == fracresFileKeyWord_->hyp_end) {

			break;
		}
		my_stream.clear();
	}

}
/**
 * @fn void manageFractureProperty(std::fstream&, FracResStage&, FracResHypothesis&)
 * @brief
 *
 * @param fin
 * @param fracresStage
 * @param fracResHypothesis
 */
void RemFracResSim::manageFractureProperty(std::fstream &fin,
		FracResStage &fracresStage, FracResHypothesis &fracResHypothesis) {
	std::string contenu_ligne;
	std::string comment;
	std::string separator;

	bool active_dfn = false;
	int id_dfn = -1;
	std::string dfnName = "";
	fractures_intersect::FracResIntersectionCategoryEnum category =
			fractures_intersect::FracResIntersectionCategoryEnum::FAULTS;
	FracResTypeEnum type = FracResTypeEnum::SINGLE_FAULTS;
	FracResBeddingAperture *fracResBedding = NULL;
	FaultArraySet *faultArray = NULL;
	FractureSet *faultSingle = NULL;
	while (getline(fin, contenu_ligne)) // tant que l'on peut mettre la ligne dans "contenu"
	{
		std::istringstream my_stream;
		cout << contenu_ligne << endl;  // on l'affiche
		my_stream.str(contenu_ligne);

		std::vector<std::string> tokenVector = splitLine(my_stream, ';');

		if (tokenVector[0] == fracresFileKeyWord_->frac_prop_selected) {
			active_dfn = (tokenVector[1] == "true");

		} else if (tokenVector[0]
				== fracresFileKeyWord_->frac_prop_fracture_id) {
			id_dfn = std::atoi(tokenVector[1].c_str());
			category = fracresStage.fractureList_[id_dfn].categoryType_;
			type = fracresStage.fractureList_[id_dfn].dataType_;
			dfnName = fracresStage.fractureList_[id_dfn].name_;
			if (category
					== fractures_intersect::FracResIntersectionCategoryEnum::FAULTS) {
				faultSingle = new FractureSet();
				faultSingle->type_category_ = category;
				faultSingle->type_array_ = faultSingle->array_type::unkown;
				faultSingle->name_ = dfnName;
			} else if (category
					== fractures_intersect::FracResIntersectionCategoryEnum::FAULT_ARRAYS) {
				faultArray = new FaultArraySet();
				faultArray->name_ = dfnName;
				faultArray->type_category_ = category;
				if (type == FracResTypeEnum::ANASTOMOSING_FAULTS) {

					faultArray->typeFamilyArray_ = FaultArraySet::anastomosing;
					FractureSet frac = FractureSet();
					frac.name_ = faultArray->name_ + "_" + "Dfn_Anastomosing";
					frac.minmaxmode_values_[6] = 10;
					frac.minmaxmode_values_[7] = 20;
					frac.minmaxmode_values_[8] = 1;

					frac.type_array_ = frac.array_type::anastomosing;
					frac.type_category_ =
							fractures_intersect::FracResIntersectionCategoryEnum::FAULTS;
					faultArray->fracturesSet_list_.push_back(frac);
				} else if (type == FracResTypeEnum::EN_ECHELON_FAULTS) {
					faultArray->typeFamilyArray_ = FaultArraySet::echelon;
					FractureSet frac = FractureSet();
					frac.name_ = faultArray->name_ + "_" + "Dfn_Echelon";
					frac.type_array_ = frac.array_type::echelon;
					frac.type_category_ =
							fractures_intersect::FracResIntersectionCategoryEnum::FAULTS;
					faultArray->fracturesSet_list_.push_back(frac);

				} else if (type == FracResTypeEnum::RELAY_FAULTS) {
					faultArray->typeFamilyArray_ = FaultArraySet::relay;
					FractureSet frac = FractureSet();
					frac.name_ = faultArray->name_ + "_" + "Dfn_Relay";
					frac.type_category_ =
							fractures_intersect::FracResIntersectionCategoryEnum::FAULTS;
					frac.type_array_ = frac.array_type::relay;
					faultArray->fracturesSet_list_.push_back(frac);
				} else if (type == FracResTypeEnum::FRACTURE_CORRIDORS) {
					faultArray->typeFamilyArray_ = FaultArraySet::corridor;
					FractureSet frac = FractureSet();
					frac.name_ = faultArray->name_ + "_" + "Dfn_Corridor";
					frac.type_array_ = frac.array_type::corridor;
					frac.type_category_ =
							fractures_intersect::FracResIntersectionCategoryEnum::FAULTS;
					faultArray->fracturesSet_list_.push_back(frac);
				} else if (type == FracResTypeEnum::FAULT_ZONES) {
					faultArray->typeFamilyArray_ = FaultArraySet::fault_zone;
					FractureSet frac = FractureSet();
					frac.name_ = faultArray->name_ + "_" + "Dfn_Fault_Zone";
					frac.type_array_ = frac.array_type::fault_zone;
					frac.type_category_ =
							fractures_intersect::FracResIntersectionCategoryEnum::FAULTS;
					faultArray->fracturesSet_list_.push_back(frac);
				}

			} else if (category
					== fractures_intersect::FracResIntersectionCategoryEnum::FRACTURE_CLUSTERS_CORRIDORS_ZONES) {
				faultArray = new FaultArraySet();
				faultArray->name_ = dfnName;
				faultArray->typeFamilyArray_ = FaultArraySet::corridor;
				FractureSet frac = FractureSet();
				frac.name_ = faultArray->name_ + "_" + "Dfn_Corridor";
				frac.type_array_ = frac.array_type::corridor;
				frac.type_category_ =
						fractures_intersect::FracResIntersectionCategoryEnum::FRACTURE_CLUSTERS_CORRIDORS_ZONES;
				faultArray->fracturesSet_list_.push_back(frac);

			} else if (category
					== fractures_intersect::FracResIntersectionCategoryEnum::FRACTURE_CLUSTERS_FAULT_ZONES) {
				faultArray = new FaultArraySet();
				faultArray->name_ = dfnName;
				faultArray->typeFamilyArray_ = FaultArraySet::fault_zone;
				FractureSet frac = FractureSet();
				frac.name_ = faultArray->name_ + "_" + "Dfn_Fault_Zone";
				frac.type_array_ = frac.array_type::fault_zone;
				frac.type_category_ =
						fractures_intersect::FracResIntersectionCategoryEnum::FRACTURE_CLUSTERS_FAULT_ZONES;
				faultArray->fracturesSet_list_.push_back(frac);

			}
			if (fracResBedding != NULL) {
				fracResBedding->is_active_ = active_dfn;
			} else if (faultSingle != NULL) {
				faultSingle->is_active_ = active_dfn;
			} else if (faultArray != NULL) {
				faultArray->is_active_ = active_dfn;
			}
		} else if (tokenVector[0] == fracresFileKeyWord_->frac_prop_index) {

			int index = std::atoi(tokenVector[1].c_str());
			if (fracResBedding != NULL) {
				fracResBedding->index_ = index;
			} else if (faultSingle != NULL) {
				faultSingle->fracture_set_index_ = index;
			} else if (faultArray != NULL) {
				faultArray->fault_array_set_index_ = index;
			}
		} else if (tokenVector[0]
				== fracresFileKeyWord_->frac_prop_clusters_type) {
			if (faultArray != NULL) {
				faultArray->fracResClusterType_ =
						(FracResClustersFaultTypeEnum) getClustersFaultType(
								tokenVector[1]);
			}
		} else if (tokenVector[0]
				== fracresFileKeyWord_->frac_prop_use_previous) {
			if (faultArray != NULL) {
				faultArray->usePreviousFracture_ = (tokenVector[1] == "true");
			}
		} else if (tokenVector[0]
				== fracresFileKeyWord_->frac_prop_corridor_fractures) {
			if (faultArray != NULL && faultArray->usePreviousFracture_) {

				faultArray->typeFamilyArray_ = FaultArraySet::corridor;
			}
		} else if (tokenVector[0]
				== fracresFileKeyWord_->frac_prop_height_correl_type) {
			if (faultArray != NULL) {

				faultArray->heightCorrelationType_ =
						(FracResGeometryCorrelationTypeEnum) getGeometryCorrelationType(
								tokenVector[1]);
			}
		} else if (tokenVector[0]
				== fracresFileKeyWord_->frac_prop_height_correl_value) {
			if (faultArray != NULL) {

				faultArray->heightCorrelationValue_ = std::strtod(
						tokenVector[1].c_str(), NULL);
			}
		} else if (tokenVector[0]
				== fracresFileKeyWord_->frac_prop_height_correl_variable) {
			if (faultArray != NULL) {

				faultArray->heightCorrelationVariable_ = (tokenVector[1]
						== "true");
			}
		} else if (tokenVector[0]
				== fracresFileKeyWord_->frac_prop_height_correl_property_name) {
			if (faultArray != NULL) {
				std::string property_name =tokenVector[1];
				 std::transform(property_name.begin(), property_name.end(), property_name.begin(), my_toupper());

				faultArray->heightCorrelationProperty_name_ = property_name;
			}
		} else if (tokenVector[0]
				== fracresFileKeyWord_->frac_prop_aperture_correl_type) {
			if (faultArray != NULL) {

				faultArray->apertureCorrelationType_ =
						(FracResGeometryCorrelationTypeEnum) getGeometryCorrelationType(
								tokenVector[1]);
			}
		} else if (tokenVector[0]
				== fracresFileKeyWord_->frac_prop_aperture_correl_value) {
			if (faultArray != NULL) {

				faultArray->apertureCorrelationValue_ = std::strtod(
						tokenVector[1].c_str(), NULL);
			}
		} else if (tokenVector[0]
				== fracresFileKeyWord_->frac_prop_aperture_correl_variable) {
			if (faultArray != NULL) {

				faultArray->apertureCorrelationVariable_ = (tokenVector[1]
						== "true");
			}
		} else if (tokenVector[0]
				== fracresFileKeyWord_->frac_prop_aperture_correl_property_name) {
			if (faultArray != NULL) {
				std::string property_name =tokenVector[1];
				 std::transform(property_name.begin(), property_name.end(), property_name.begin(), my_toupper());

				faultArray->apertureCorrelationProperty_name_ = property_name;
			}
		} else if (tokenVector[0]
				== fracresFileKeyWord_->frac_prop_fault_core_type) {
			if (faultArray != NULL) {

				faultArray->faultCoreType_ = (FracResDataTypeEnum) getDataType(
						tokenVector[1]);
			}
		} else if (tokenVector[0]
				== fracresFileKeyWord_->frac_prop_fault_core_value) {
			if (faultArray != NULL) {

				faultArray->faultCoreValue_ = std::strtod(
						tokenVector[1].c_str(), NULL);
			}
		} else if (tokenVector[0]
				== fracresFileKeyWord_->frac_prop_fault_core_property_name) {
			if (faultArray != NULL) {
				std::string property_name =tokenVector[1];
				 std::transform(property_name.begin(), property_name.end(), property_name.begin(), my_toupper());

				faultArray->faultCoreProperty_name_ = property_name;
			}
		} else if (tokenVector[0]
				== fracresFileKeyWord_->frac_prop_fault_core_distribution) {
			if (faultArray != NULL) {

				faultArray->faultCoreDistribution_name_ = tokenVector[1];
			}
		} else if (tokenVector[0]
				== fracresFileKeyWord_->frac_prop_rate_reactivate) {
			if (faultArray != NULL) {

				faultArray->rate_reactive_ = std::strtod(tokenVector[1].c_str(),
						NULL);
			}
		} else if (tokenVector[0] == fracresFileKeyWord_->intersection_begin) {

			int indice = -1;
			fractures_intersect::FracResIntersectionTypeEnum tab_type_intersection =
					manageIntersection(fin, indice);
			if (indice < 0)
				continue;
			if (faultSingle != NULL) {
				if (indice >= 0)
					faultSingle->tab_type_intersection_[indice] =
							tab_type_intersection;
			} else if (faultArray != NULL) {
				if (indice >= 0)
					faultArray->tab_type_intersection_[indice] =
							tab_type_intersection;
			}
		} else if (tokenVector[0] == fracresFileKeyWord_->distribution_begin) {

			if (faultSingle != NULL || faultArray != NULL) {
				manageTotalDistribution(fin, faultArray, faultSingle);
			}
		} else if (tokenVector[0] == fracresFileKeyWord_->geometry_begin) {
			if (faultSingle != NULL || faultArray != NULL) {
				manageTotalGeometry(fin, faultArray, faultSingle);
			}
		} else if (tokenVector[0] == fracresFileKeyWord_->bed_parrallel_begin) {
			manageBedParallel(fin, fracresStage, fracResHypothesis);
		} else if (tokenVector[0] == fracresFileKeyWord_->frac_prop_end) {

			int index = fracResHypothesis.bed_parallel_List_.size();
			if (index > 0 && !dfnName.empty()) {
				fracResHypothesis.bed_parallel_List_[index - 1]->name_ =
						dfnName;
				fracResHypothesis.bed_parallel_List_[index - 1]->fracres_hypothesis_name_ =
						fracResHypothesis.name_;
				fracResHypothesis.bed_parallel_List_[index - 1]->stage_name_ =
						fracresStage.name_;
			}

			if (faultSingle != NULL) {
				fracResHypothesis.fractureList_.push_back(faultSingle);
				fracResHypothesis.type_dfn_.push_back(1);
			}
			if (faultArray != NULL) {

				if (faultArray->typeFamilyArray_ == FaultArraySet::corridor) {
					fracResHypothesis.type_dfn_.push_back(4);
					fracResHypothesis.faultCorridorList_.push_back(faultArray);
				} else if (faultArray->typeFamilyArray_
						== FaultArraySet::fault_zone) {
					fracResHypothesis.type_dfn_.push_back(5);
					fracResHypothesis.faultZoneList_.push_back(faultArray);
				} else {
					fracResHypothesis.type_dfn_.push_back(2);
					fracResHypothesis.faultArrayList_.push_back(faultArray);
				}
			}
			break;
		}
		my_stream.clear();
	}

}
// SINGLE FAULT or FractureSet
/**
 * @fn void setGeometryDipAzimuth(int, double, std::vector<std::string>&, const std::vector<double>&, FractureSet*)
 * @brief
 *
 * @param geometry_data_type
 * @param geometry_value
 * @param name
 * @param min_max_mode
 * @param faultSingle
 */
void RemFracResSim::setGeometryDipAzimuth(int geometry_data_type,
		double geometry_value, std::vector<std::string> &name,
		const std::vector<double> &min_max_mode, FractureSet *faultSingle) {

	if (std::isnan(geometry_value))
		geometry_value = 1e-6;
	if (geometry_data_type == 0) {
		faultSingle->distribution_azimut_type_index_ = 1;
		faultSingle->distr_azimut_val_ = geometry_value;
	} else if (geometry_data_type == 1) {
		faultSingle->distribution_azimut_type_index_ = 0;
		if (name[1] == "NORMAL") {
			faultSingle->azimuth_type_= FractureSet::distribution_type::normal;
		} else if (name[1] == " TRIANGULAR") {
			faultSingle->azimuth_type_= FractureSet::distribution_type::triangulated;
			faultSingle->minmaxmode_values_[8] = min_max_mode[1];
		} else if (name[1] == "UNIFORM") {
			faultSingle->azimuth_type_= FractureSet::distribution_type::uniform;
		} else if (name[1] == "LOGNORMAL") {
			faultSingle->azimuth_type_= FractureSet::distribution_type::lognormal;
		}

		faultSingle->minmaxmode_values_[6] = min_max_mode[0];
		faultSingle->minmaxmode_values_[7] = min_max_mode[2];
	} else if (geometry_data_type == 2) {
		std::string property_name =name[0];
		 std::transform(property_name.begin(), property_name.end(), property_name.begin(), my_toupper());

		faultSingle->distribution_azimut_type_index_ = 2;
		faultSingle->azimut_property_name_=property_name;
	}else if (geometry_data_type == 3) {
		faultSingle->distribution_azimut_type_index_ = 0;
		faultSingle->azimuth_type_= FractureSet::distribution_type::uniform;
		faultSingle->minmaxmode_values_[6] = min_max_mode[0];
		faultSingle->minmaxmode_values_[7] = min_max_mode[2];
	}
	else {
		faultSingle->distribution_azimut_type_index_ = geometry_data_type;
	}
}
/**
 * @fn void setGeometryDipAngle(int, double, std::vector<std::string>&, const std::vector<double>&, FractureSet*)
 * @brief
 *
 * @param geometry_data_type
 * @param geometry_value
 * @param name
 * @param min_max_mode
 * @param faultSingle
 */
void RemFracResSim::setGeometryDipAngle(int geometry_data_type,
		double geometry_value, std::vector<std::string> &name,
		const std::vector<double> &min_max_mode, FractureSet *faultSingle) {
	if (std::isnan(geometry_value))
		geometry_value = 90;
	if (geometry_data_type == 0) {
		faultSingle->distribution_angle_type_index_ = 1;
		faultSingle->distr_angle_val_ = geometry_value;
	} else if (geometry_data_type == 1) {
		faultSingle->distribution_angle_type_index_ = 0;

		if (name[1] == "NORMAL") {
			faultSingle->dip_angle_type_ =FractureSet::distribution_type::normal;
		} else if (name[1] == "TRIANGULAR") {

			faultSingle->dip_angle_type_ =FractureSet::distribution_type::triangulated;
			faultSingle->minmaxmode_values_[11] = min_max_mode[1];

		} else if (name[1] == "UNIFORM") {
			faultSingle->dip_angle_type_ =FractureSet::distribution_type::uniform;
		} else if (name[1] == "LOGNORMAL") {
			faultSingle->dip_angle_type_ =FractureSet::distribution_type::lognormal;
		}

		faultSingle->minmaxmode_values_[9] = min_max_mode[0];
		faultSingle->minmaxmode_values_[10] = min_max_mode[2];

	}else if (geometry_data_type == 2) {
		faultSingle->distribution_angle_type_index_ = 2;
		std::string property_name =name[0];
		 std::transform(property_name.begin(), property_name.end(), property_name.begin(), my_toupper());

		faultSingle->dip_angle_property_name_=property_name;
	} else {
		faultSingle->distribution_angle_type_index_ = geometry_data_type;
	}

}
/**
 * @fn void setGeometryLength(int, double, std::vector<std::string>&, const std::vector<double>&, FractureSet*)
 * @brief
 *
 * @param geometry_data_type
 * @param geometry_value
 * @param name
 * @param min_max_mode
 * @param faultSingle
 */
void RemFracResSim::setGeometryLength(int geometry_data_type,
		double geometry_value, std::vector<std::string> &name,
		const std::vector<double> &min_max_mode, FractureSet *faultSingle) {
	if (std::isnan(geometry_value))
		geometry_value = 1e-6;
	if (geometry_data_type == 0) {
		faultSingle->distribution_length1_type_index_ = 1;
		faultSingle->distr_length1_val_ = geometry_value;
	} else if (geometry_data_type == 1) {
		faultSingle->distribution_length1_type_index_ = 0;

		if (name[1] == "NORMAL") {
			faultSingle->length1_type_ = FractureSet::distribution_type::normal;
		} else if (name[1] == "TRIANGULAR") {

			faultSingle->length1_type_ =
					FractureSet::distribution_type::triangulated;
			faultSingle->minmaxmode_values_[2] = min_max_mode[1];

		} else if (name[1] == "UNIFORM") {
			faultSingle->length1_type_ =FractureSet::distribution_type::uniform;
		} else if (name[1] == "LOGNORMAL") {
			faultSingle->length1_type_ =FractureSet::distribution_type::lognormal;
		}

		faultSingle->minmaxmode_values_[0] = min_max_mode[0];
		faultSingle->minmaxmode_values_[1] = min_max_mode[2];

	}else if (geometry_data_type == 2) {
		faultSingle->distribution_length1_type_index_ = 2;
		std::string property_name =name[0];
		 std::transform(property_name.begin(), property_name.end(), property_name.begin(), my_toupper());

		faultSingle->length1_property_name_ = property_name;
	} else {
		faultSingle->distribution_length1_type_index_ = geometry_data_type;
	}
}
/**
 * @fn void setGeometryHeight(int, double, std::vector<std::string>&, const std::vector<double>&, FractureSet*)
 * @brief
 *
 * @param geometry_data_type
 * @param geometry_value
 * @param name
 * @param min_max_mode
 * @param faultSingle
 */
void RemFracResSim::setGeometryHeight(int geometry_data_type,
		double geometry_value, std::vector<std::string> &name,
		const std::vector<double> &min_max_mode, FractureSet *faultSingle) {
	if (std::isnan(geometry_value))
		geometry_value = 1e-6;
	if (geometry_data_type == 0) {
		faultSingle->distribution_length2_type_index_ = 1;
		faultSingle->distr_length2_val_ = geometry_value;
	} else if (geometry_data_type == 1) {
		faultSingle->distribution_length2_type_index_ = 0;

		if (name[1] == "NORMAL") {
			faultSingle->length2_type_ = FractureSet::distribution_type::normal;
		} else if (name[1] == "TRIANGULAR") {

			faultSingle->length2_type_ =FractureSet::distribution_type::triangulated;
			faultSingle->minmaxmode_values_[5] = min_max_mode[1];

		} else if (name[1] == "UNIFORM") {
			faultSingle->length2_type_ =FractureSet::distribution_type::uniform;
		} else if (name[1] == "LOGNORMAL") {
			faultSingle->length2_type_ =FractureSet::distribution_type::lognormal;
		}

		faultSingle->minmaxmode_values_[3] = min_max_mode[0];
		faultSingle->minmaxmode_values_[4] = min_max_mode[2];

	}else if (geometry_data_type == 2) {
		faultSingle->distribution_length2_type_index_ = 2;
		std::string property_name =name[0];
		 std::transform(property_name.begin(), property_name.end(), property_name.begin(), my_toupper());

		faultSingle->length2_property_name_ = property_name;
	} else {
		faultSingle->distribution_length2_type_index_ = geometry_data_type;
	}
}
/**
 * @fn void setGeometryAperture(int, double, std::vector<std::string>&, const std::vector<double>&, FractureSet*)
 * @brief
 *
 * @param geometry_data_type
 * @param geometry_value
 * @param name
 * @param min_max_mode
 * @param faultSingle
 */
void RemFracResSim::setGeometryAperture(int geometry_data_type,
		double geometry_value, std::vector<std::string> &name,
		const std::vector<double> &min_max_mode, FractureSet *faultSingle) {
	if (std::isnan(geometry_value))
		geometry_value = 1e-6;
	if (geometry_data_type == 0) {
		faultSingle->distribution_aperture_type_index_ = 1;
		faultSingle->distr_aperture_val_ = geometry_value;
	} else if (geometry_data_type == 1) {
		faultSingle->distribution_aperture_type_index_ = 0;

		if (name[1] == "NORMAL") {
			faultSingle->aperture_type_ =FractureSet::distribution_type::normal;
		} else if (name[1] == "TRIANGULAR") {

			faultSingle->aperture_type_ =FractureSet::distribution_type::triangulated;
			faultSingle->minmaxmode_values_[14] = min_max_mode[1];

		} else if (name[1] == "UNIFORM") {
			faultSingle->aperture_type_ =FractureSet::distribution_type::uniform;
		} else if (name[1] == "LOGNORMAL") {
			faultSingle->aperture_type_ =FractureSet::distribution_type::lognormal;
		}

		faultSingle->minmaxmode_values_[12] = min_max_mode[0];
		faultSingle->minmaxmode_values_[13] = min_max_mode[2];

	}else if (geometry_data_type == 2) {
		faultSingle->distribution_aperture_type_index_ = 2;
		std::string property_name =name[0];
		 std::transform(property_name.begin(), property_name.end(), property_name.begin(), my_toupper());

		faultSingle->aperture_property_name_ = property_name;
	} else {
		faultSingle->distribution_aperture_type_index_ = geometry_data_type;
	}
}
/**
 * @fn void setGeometrySinglefault(int, int, double, std::vector<std::string>&, const std::vector<double>&, FractureSet*)
 * @brief
 *
 * @param geometry_parameter_type
 * @param geometry_data_type
 * @param geometry_value
 * @param name
 * @param min_max_mode
 * @param faultSingle
 */
void RemFracResSim::setGeometrySinglefault(int geometry_parameter_type,
		int geometry_data_type, double geometry_value,
		std::vector<std::string> &name, const std::vector<double> &min_max_mode,
		FractureSet *faultSingle) {

	if (geometry_parameter_type == 2 || geometry_parameter_type == 16) {
		if (std::isnan(geometry_value))
			geometry_value = 180;
		setGeometryDipAzimuth(geometry_data_type, geometry_value, name,
				min_max_mode, faultSingle);
	} else if (geometry_parameter_type == 3) {
		if (std::isnan(geometry_value))
			geometry_value = 45;
		setGeometryDipAngle(geometry_data_type, geometry_value, name,
				min_max_mode, faultSingle);
	} else if (geometry_parameter_type == 4) {
		if (std::isnan(geometry_value))
			geometry_value = 1;
		setGeometryLength(geometry_data_type, geometry_value, name,
				min_max_mode, faultSingle);
	} else if (geometry_parameter_type == 5) {
		if (std::isnan(geometry_value))
			geometry_value = 1;
		setGeometryHeight(geometry_data_type, geometry_value, name,
				min_max_mode, faultSingle);
	} else if (geometry_parameter_type == 6) {
		if (std::isnan(geometry_value))
			geometry_value = 0.1;
		setGeometryAperture(geometry_data_type, geometry_value, name,
				min_max_mode, faultSingle);
	}
}
/**
 * @fn void setGeometrySingleFaultCorridor(int, int, double, std::vector<std::string>&, const std::vector<double>&, FractureSet*)
 * @brief
 *
 * @param geometry_parameter_type
 * @param geometry_data_type
 * @param geometry_value
 * @param name
 * @param min_max_mode
 * @param faultSingle
 */
void RemFracResSim::setGeometrySingleFaultCorridor(int geometry_parameter_type,
		int geometry_data_type, double geometry_value,
		std::vector<std::string> &name, const std::vector<double> &min_max_mode,
		FractureSet *faultSingle) {

	if (geometry_parameter_type == 4) {
		if (std::isnan(geometry_value))
			geometry_value = 1;
		setGeometryLength(geometry_data_type, geometry_value, name,
				min_max_mode, faultSingle);
	} else if (geometry_parameter_type == 5) {
		if (std::isnan(geometry_value))
			geometry_value = 1;
		setGeometryHeight(geometry_data_type, geometry_value, name,
				min_max_mode, faultSingle);
	} else if (geometry_parameter_type == 6) {
		if (std::isnan(geometry_value))
			geometry_value = 0.1;
		setGeometryAperture(geometry_data_type, geometry_value, name,
				min_max_mode, faultSingle);
	}
}

/**
 * @fn void setGeometrySingleFaultZone(int, int, double, std::vector<std::string>&, const std::vector<double>&, FractureSet*)
 * @brief
 *
 * @param geometry_parameter_type
 * @param geometry_data_type
 * @param geometry_value
 * @param name
 * @param min_max_mode
 * @param faultSingle
 */
void RemFracResSim::setGeometrySingleFaultZone(int geometry_parameter_type,
		int geometry_data_type, double geometry_value,
		std::vector<std::string> &name, const std::vector<double> &min_max_mode,
		FaultArraySet *faultArray) {

	if (geometry_parameter_type == 4) {
		if (std::isnan(geometry_value))
			geometry_value = 1;
		setGeometryLength(geometry_data_type, geometry_value, name,
				min_max_mode, &faultArray->fracturesSet_list_[0]);
	} else if (geometry_parameter_type == 5) {
		if (std::isnan(geometry_value))
			geometry_value = 1;
		setGeometryHeight(geometry_data_type, geometry_value, name,
				min_max_mode, &faultArray->fracturesSet_list_[0]);
	} else if (geometry_parameter_type == 6) {
		if (std::isnan(geometry_value))
			geometry_value = 1;
		setGeometryAperture(geometry_data_type, geometry_value, name,
				min_max_mode, &faultArray->fracturesSet_list_[0]);
	}
}

// ARRAY
/**
 * @fn void setGeometryAzimuthArray(int, double, std::vector<std::string>&, const std::vector<double>&, FaultArraySet*)
 * @brief
 *
 * @param geometry_data_type
 * @param geometry_value
 * @param name
 * @param min_max_mode
 * @param faultArray
 */
void RemFracResSim::setGeometryAzimuthArray(int geometry_data_type,
		double geometry_value, std::vector<std::string> &name,
		const std::vector<double> &min_max_mode, FaultArraySet *faultArray) {

	if (geometry_data_type == 0) {
		if (std::isnan(geometry_value))
			geometry_value = 180;
		faultArray->distribution_azimut_type_index_ = 1;
		faultArray->distr_azimut_val_ = geometry_value;
	} else if (geometry_data_type == 1) {
		faultArray->distribution_azimut_type_index_ = 0;
		if (name[1] == "NORMAL") {
			faultArray->azimuth_type_ =
					FaultArraySet::distribution_type::normal;
		} else if (name[1] == "TRIANGULAR") {
			faultArray->azimuth_type_ =
					FaultArraySet::distribution_type::triangulated;
			faultArray->minmaxmode_values_[8] = min_max_mode[1];
		} else if (name[1] == "UNIFORM") {
			faultArray->azimuth_type_ =
					FaultArraySet::distribution_type::uniform;
		} else if (name[1] == "LOGNORMAL") {
			faultArray->azimuth_type_ =
					FaultArraySet::distribution_type::lognormal;
		}

		faultArray->minmaxmode_values_[6] = min_max_mode[0];
		faultArray->minmaxmode_values_[7] = min_max_mode[2];
	}else if (geometry_data_type == 2) {
		faultArray->distribution_azimut_type_index_ = 2;
		std::string property_name =name[0];
		 std::transform(property_name.begin(), property_name.end(), property_name.begin(), my_toupper());

		faultArray->azimut_property_name_ = property_name;
	} else {
		faultArray->distribution_azimut_type_index_ = geometry_data_type;
	}
}
/**
 * @fn void setGeometryAngleArray(int, double, std::vector<std::string>&, const std::vector<double>&, FaultArraySet*)
 * @brief
 *
 * @param geometry_data_type
 * @param geometry_value
 * @param name
 * @param min_max_mode
 * @param faultArray
 */
void RemFracResSim::setGeometryAngleArray(int geometry_data_type,
		double geometry_value, std::vector<std::string> &name,
		const std::vector<double> &min_max_mode, FaultArraySet *faultArray) {
	if (std::isnan(geometry_value))
		geometry_value = 45;
	if (geometry_data_type == 0) {
		faultArray->distribution_angle_type_index_ = 1;
		faultArray->distr_angle_val_ = geometry_value;
	} else if (geometry_data_type == 1 || geometry_data_type == 3) {
		faultArray->distribution_angle_type_index_ = 0;

		if (name[1] == "NORMAL") {
			faultArray->dip_angle_type_ =
					FaultArraySet::distribution_type::normal;
		} else if (name[1] == "TRIANGULAR") {

			faultArray->dip_angle_type_ =
					FaultArraySet::distribution_type::triangulated;
			faultArray->minmaxmode_values_[11] = min_max_mode[1];

		} else if (name[1] == "UNIFORM") {
			faultArray->dip_angle_type_ =
					FaultArraySet::distribution_type::uniform;
		} else if (name[1] == "LOGNORMAL") {
			faultArray->dip_angle_type_ =
					FaultArraySet::distribution_type::lognormal;
		}

		faultArray->minmaxmode_values_[9] = min_max_mode[0];
		faultArray->minmaxmode_values_[10] = min_max_mode[2];

	} else if (geometry_data_type == 2) {
		faultArray->distribution_angle_type_index_ = 2;
		std::string property_name =name[0];
		 std::transform(property_name.begin(), property_name.end(), property_name.begin(), my_toupper());

		faultArray->dip_angle_property_name_ = property_name;
	}else {
		faultArray->distribution_angle_type_index_ = geometry_data_type;
	}

}
/**
 * @fn void setGeometryLengthArray(int, double, std::vector<std::string>&, const std::vector<double>&, FaultArraySet*)
 * @brief
 *
 * @param geometry_data_type
 * @param geometry_value
 * @param name
 * @param min_max_mode
 * @param faultArray
 */
void RemFracResSim::setGeometryLengthArray(int geometry_data_type,
		double geometry_value, std::vector<std::string> &name,
		const std::vector<double> &min_max_mode, FaultArraySet *faultArray) {
	if (std::isnan(geometry_value))
		geometry_value = 1;
	if (geometry_data_type == 0) {
		faultArray->distribution_length1_type_index_ = 1;
		faultArray->distr_length1_val_ = geometry_value;
	} else if (geometry_data_type == 1) {
		faultArray->distribution_length1_type_index_ = 0;

		if (name[1] == "NORMAL") {
			faultArray->length1_type_ =
					FaultArraySet::distribution_type::normal;
		} else if (name[1] == "TRIANGULAR") {

			faultArray->length1_type_ =
					FaultArraySet::distribution_type::triangulated;
			faultArray->minmaxmode_values_[2] = min_max_mode[1];

		} else if (name[1] == "UNIFORM") {
			faultArray->length1_type_ =
					FaultArraySet::distribution_type::uniform;
		} else if (name[1] == "LOGNORMAL") {
			faultArray->length1_type_ =
					FaultArraySet::distribution_type::lognormal;
		}

		faultArray->minmaxmode_values_[0] = min_max_mode[0];
		faultArray->minmaxmode_values_[1] = min_max_mode[2];

	}else if (geometry_data_type == 2) {
		faultArray->distribution_length1_type_index_ = 2;
		std::string property_name =name[0];
		 std::transform(property_name.begin(), property_name.end(), property_name.begin(), my_toupper());

		faultArray->length1_property_name_ = property_name;
	} else {
		faultArray->distribution_length1_type_index_ = geometry_data_type;
	}
}
/**
 * @fn void setGeometryHeightArray(int, double, std::vector<std::string>&, const std::vector<double>&, FaultArraySet*)
 * @brief
 *
 * @param geometry_data_type
 * @param geometry_value
 * @param name
 * @param min_max_mode
 * @param faultArray
 */
void RemFracResSim::setGeometryHeightArray(int geometry_data_type,
		double geometry_value, std::vector<std::string> &name,
		const std::vector<double> &min_max_mode, FaultArraySet *faultArray) {
	if (std::isnan(geometry_value))
		geometry_value = 1;
	if (geometry_data_type == 0) {
		faultArray->distribution_length2_type_index_ = 1;
		faultArray->distr_length2_val_ = geometry_value;
	} else if (geometry_data_type == 1) {
		faultArray->distribution_length2_type_index_ = 0;

		if (name[1] == "NORMAL") {
			faultArray->length2_type_ =
					FaultArraySet::distribution_type::normal;

		} else if (name[1] == "TRIANGULAR") {

			faultArray->length2_type_ =
					FaultArraySet::distribution_type::triangulated;
			faultArray->minmaxmode_values_[5] = min_max_mode[1];

		} else if (name[1] == "UNIFORM") {
			faultArray->length2_type_ =
					FaultArraySet::distribution_type::uniform;
		} else if (name[1] == "LOGNORMAL") {
			faultArray->length2_type_ =
					FaultArraySet::distribution_type::lognormal;
		}

		faultArray->minmaxmode_values_[3] = min_max_mode[0];
		faultArray->minmaxmode_values_[4] = min_max_mode[2];

	}else if (geometry_data_type == 2) {
		faultArray->distribution_length2_type_index_ = 2;
		std::string property_name =name[0];
		 std::transform(property_name.begin(), property_name.end(), property_name.begin(), my_toupper());

		faultArray->length2_property_name_ = property_name;
	} else {
		faultArray->distribution_length2_type_index_ = geometry_data_type;
	}
}
/**
 * @fn void setGeometryDamageHangingArray(int, double, std::vector<std::string>&, const std::vector<double>&, FaultArraySet*)
 * @brief
 *
 * @param geometry_data_type
 * @param geometry_value
 * @param name
 * @param min_max_mode
 * @param faultArray
 */

void RemFracResSim::setGeometryDamageHangingArray(int geometry_data_type,
		double geometry_value, std::vector<std::string> &name,
		const std::vector<double> &min_max_mode, FaultArraySet *faultArray) {
	if (std::isnan(geometry_value))
		geometry_value = 1e-6;
	if (geometry_data_type == 0) {
		faultArray->distribution_damage_zone_hanging_index_ = 1;
		faultArray->distr_damage_zone_hanging_val_ = geometry_value;
	} else if (geometry_data_type == 1) {
		faultArray->distribution_damage_zone_hanging_index_ = 0;

		if (name[1] == "NORMAL") {
			faultArray->damage_zone_hanging_type_ =
					FaultArraySet::distribution_type::normal;
		} else if (name[1] == "TRIANGULAR") {

			faultArray->damage_zone_hanging_type_ =
					FaultArraySet::distribution_type::triangulated;
			faultArray->minmaxmode_damage_zone_hanging_values_[1] =
					min_max_mode[1];

		} else if (name[1] == "UNIFORM") {
			faultArray->damage_zone_hanging_type_ =
					FaultArraySet::distribution_type::uniform;
		} else if (name[1] == "LOGNORMAL") {
			faultArray->damage_zone_hanging_type_ =
					FaultArraySet::distribution_type::lognormal;
		}

		faultArray->minmaxmode_damage_zone_hanging_values_[0] = min_max_mode[0];
		faultArray->minmaxmode_damage_zone_hanging_values_[2] = min_max_mode[2];

	}else if (geometry_data_type == 2) {
		faultArray->distribution_damage_zone_hanging_index_ = 2;
		std::string property_name =name[0];
		 std::transform(property_name.begin(), property_name.end(), property_name.begin(), my_toupper());

		faultArray->damage_hanging_wall_property_name_ = property_name;
	} else {
		faultArray->distribution_damage_zone_hanging_index_ =
				geometry_data_type;
	}
}
/**
 * @fn void setGeometryDamageFootwallArray(int, double, std::vector<std::string>&, const std::vector<double>&, FaultArraySet*)
 * @brief
 *
 * @param geometry_data_type
 * @param geometry_value
 * @param name
 * @param min_max_mode
 * @param faultArray
 */
void RemFracResSim::setGeometryDamageFootwallArray(int geometry_data_type,
		double geometry_value, std::vector<std::string> &name,
		const std::vector<double> &min_max_mode, FaultArraySet *faultArray) {
	if (std::isnan(geometry_value))
		geometry_value = 1e-6;
	if (geometry_data_type == 0) {
		faultArray->distribution_damage_zone_footwall_index_ = 1;
		faultArray->distr_damage_zone_footwall_val_ = geometry_value;
	} else if (geometry_data_type == 1) {
		faultArray->distribution_damage_zone_footwall_index_ = 0;

		if (name[1] == "NORMAL") {
			faultArray->damage_zone_footwall_type_ =
					FaultArraySet::distribution_type::normal;
		} else if (name[1] == "TRIANGULAR") {

			faultArray->damage_zone_footwall_type_ =
					FaultArraySet::distribution_type::triangulated;
			faultArray->minmaxmode_damage_zone_footwall_values_[1] =
					min_max_mode[1];

		} else if (name[1] == "UNIFORM") {
			faultArray->damage_zone_footwall_type_ =
					FaultArraySet::distribution_type::uniform;
		} else if (name[1] == "LOGNORMAL") {
			faultArray->damage_zone_footwall_type_ =
					FaultArraySet::distribution_type::lognormal;
		}

		faultArray->minmaxmode_damage_zone_footwall_values_[0] =
				min_max_mode[0];
		faultArray->minmaxmode_damage_zone_footwall_values_[2] =
				min_max_mode[2];

	}else if (geometry_data_type == 2) {
		faultArray->distribution_damage_zone_footwall_index_ = 2;
		std::string property_name =name[0];
		 std::transform(property_name.begin(), property_name.end(), property_name.begin(), my_toupper());

		faultArray->damage_foot_wall_property_name_ = property_name ;
	} else {
		faultArray->distribution_damage_zone_footwall_index_ =
				geometry_data_type;
	}
}
/**
 * @fn void setGeometryApertureArray(int, double, std::vector<std::string>&, const std::vector<double>&, FaultArraySet*)
 * @brief
 *
 * @param geometry_data_type
 * @param geometry_value
 * @param name
 * @param min_max_mode
 * @param faultArray
 */
void RemFracResSim::setGeometryApertureArray(int geometry_data_type,
		double geometry_value, std::vector<std::string> &name,
		const std::vector<double> &min_max_mode, FaultArraySet *faultArray) {
	if (std::isnan(geometry_value))
		geometry_value = 1e-6;
	if (geometry_data_type == 0) {
		faultArray->distribution_aperture_type_index_ = 1;
		faultArray->distr_aperture_val_ = geometry_value;
	} else if (geometry_data_type == 1) {
		faultArray->distribution_aperture_type_index_ = 0;

		if (name[1] == "NORMAL") {
			faultArray->aperture_type_ =
					FaultArraySet::distribution_type::normal;
		} else if (name[1] == "TRIANGULAR") {

			faultArray->aperture_type_ =
					FaultArraySet::distribution_type::triangulated;
			faultArray->minmaxmode_values_[14] = min_max_mode[1];

		} else if (name[1] == "UNIFORM") {
			faultArray->aperture_type_ =
					FaultArraySet::distribution_type::uniform;
		} else if (name[1] == "LOGNORMAL") {
			faultArray->aperture_type_ =
					FaultArraySet::distribution_type::lognormal;
		}

		faultArray->minmaxmode_values_[12] = min_max_mode[0];
		faultArray->minmaxmode_values_[13] = min_max_mode[2];

	} else if (geometry_data_type == 2) {
		faultArray->distribution_aperture_type_index_=2;
		std::string property_name =name[0];
		 std::transform(property_name.begin(), property_name.end(), property_name.begin(), my_toupper());

		faultArray->aperture_property_name_ = property_name;
	}else {
		faultArray->distribution_aperture_type_index_ = geometry_data_type;
	}
}

/**
 * @fn void setGeometryArray(int, int, double, std::vector<std::string>&, const std::vector<double>&, FaultArraySet*)
 * @brief
 *
 * @param geometry_parameter_type
 * @param geometry_data_type
 * @param geometry_value
 * @param name
 * @param min_max_mode
 * @param faultArray
 */
void RemFracResSim::setGeometryArray(int geometry_parameter_type,
		int geometry_data_type, double geometry_value,
		std::vector<std::string> &name, const std::vector<double> &min_max_mode,
		FaultArraySet *faultArray) {

	if (geometry_parameter_type == 1) {
		setGeometryAzimuthArray(geometry_data_type, geometry_value, name,
				min_max_mode, faultArray);
	} else if (geometry_parameter_type == 3) {
		setGeometryAngleArray(geometry_data_type, geometry_value, name,
				min_max_mode, faultArray);
	} else if (geometry_parameter_type == 4) {
		setGeometryLengthArray(geometry_data_type, geometry_value, name,
				min_max_mode, faultArray);
	} else if (geometry_parameter_type == 7) {
		setGeometryHeightArray(geometry_data_type, geometry_value, name,
				min_max_mode, faultArray);
	} else if (geometry_parameter_type == 5 || geometry_parameter_type == 6) {
		setGeometryApertureArray(geometry_data_type, geometry_value, name,
				min_max_mode, faultArray);
	}
}
/**
 * @fn void setGeometryCorridor(int, int, double, std::vector<std::string>&, const std::vector<double>&, FaultArraySet*)
 * @brief
 *
 * @param geometry_parameter_type
 * @param geometry_data_type
 * @param geometry_value
 * @param name
 * @param min_max_mode
 * @param faultArray
 */
void RemFracResSim::setGeometryCorridor(int geometry_parameter_type,
		int geometry_data_type, double geometry_value,
		std::vector<std::string> &name, const std::vector<double> &min_max_mode,
		FaultArraySet *faultArray) {

	if (geometry_parameter_type == 1) {
		setGeometryAzimuthArray(geometry_data_type, geometry_value, name,
				min_max_mode, faultArray);
	} else if (geometry_parameter_type == 3) {
		faultArray->distribution_angle_type_index_ = 1;
		faultArray->distr_angle_val_ = 90;
	} else if (geometry_parameter_type == 4) {
		setGeometryLengthArray(geometry_data_type, geometry_value, name,
				min_max_mode, faultArray);
	} else if (geometry_parameter_type == 5) {
		setGeometryHeightArray(geometry_data_type, geometry_value, name,
						min_max_mode, faultArray);
		//faultArray->distribution_length2_type_index_ = 1;
		//faultArray->distr_length2_val_ = geometry_value;

	} else if (geometry_parameter_type == 6 || geometry_parameter_type == 7) {
		setGeometryApertureArray(geometry_data_type, geometry_value, name,
				min_max_mode, faultArray);
	}
}
/**
 * @fn void setGeometryFaultZone(int, int, double, std::vector<std::string>&, const std::vector<double>&, FaultArraySet*)
 * @brief
 *
 * @param geometry_parameter_type
 * @param geometry_data_type
 * @param geometry_value
 * @param name
 * @param min_max_mode
 * @param faultArray
 */

void RemFracResSim::setGeometryFaultZone(int geometry_parameter_type,
		int geometry_data_type, double geometry_value,
		std::vector<std::string> &name, const std::vector<double> &min_max_mode,
		FaultArraySet *faultArray) {

	if (geometry_parameter_type == 8) {
		setGeometryAzimuthArray(geometry_data_type, geometry_value, name,
				min_max_mode, faultArray);
	} else if (geometry_parameter_type == 9) {
		setGeometryAngleArray(geometry_data_type, geometry_value, name,
				min_max_mode, faultArray);
	} else if (geometry_parameter_type == 10) {
		setGeometryLengthArray(geometry_data_type, geometry_value, name,
				min_max_mode, faultArray);
	} else if (geometry_parameter_type == 11) {

		setGeometryHeightArray(geometry_data_type, geometry_value, name,
				min_max_mode, faultArray);

	} else if (geometry_parameter_type == 12 || geometry_parameter_type == 13) {
		setGeometryApertureArray(geometry_data_type, geometry_value, name,
				min_max_mode, faultArray);
	} else if (geometry_parameter_type == 14) {
		if (std::isnan(geometry_value))
			geometry_value = 0.5;
		setGeometryDamageHangingArray(geometry_data_type, geometry_value, name,
				min_max_mode, faultArray);
	} else if (geometry_parameter_type == 15) {
		if (std::isnan(geometry_value))
			geometry_value = 0.5;
		setGeometryDamageFootwallArray(geometry_data_type, geometry_value, name,
				min_max_mode, faultArray);
	}
}

/**
 * @fn void manageTotalGeometry(std::fstream&, FaultArraySet*, FractureSet*)
 * @brief
 *
 * @param fin
 * @param faultArray
 * @param faultSingle
 */

void RemFracResSim::manageTotalGeometry(std::fstream &fin,
		FaultArraySet *faultArray, FractureSet *faultSingle) {
	int geometry_parameter_type = -1;
	int geometry_array_type = -1;
	int geometry_data_type = -1;
	std::string unit = "M";
	double geometry_value = -1.0;
	std::vector<double> min_max_mode;
	std::vector<std::string> name = manageGeometry(fin, geometry_parameter_type,
			geometry_array_type, geometry_data_type, geometry_value,
			min_max_mode);
	if (faultSingle != NULL && geometry_array_type == 0) {
		setGeometrySinglefault(geometry_parameter_type, geometry_data_type,
				geometry_value, name, min_max_mode, faultSingle);
	}
	if (faultArray != NULL) {
		if (geometry_array_type == 1) {  // array

			setGeometryArray(geometry_parameter_type, geometry_data_type,
					geometry_value, name, min_max_mode, faultArray);

		} else if (geometry_array_type == 2) {   // array_common

			setGeometrySinglefault(geometry_parameter_type, geometry_data_type,
					geometry_value, name, min_max_mode,
					&faultArray->fracturesSet_list_[0]);

		} else if (geometry_array_type == 3) {    // array_familly1

			setGeometrySinglefault(geometry_parameter_type, geometry_data_type,
					geometry_value, name, min_max_mode,
					&faultArray->fracturesSet_list_[0]);
		} else if (geometry_array_type == 4) {    // array_familly2
			setGeometrySinglefault(geometry_parameter_type, geometry_data_type,
					geometry_value, name, min_max_mode,
					&faultArray->fracturesSet_list_[0]);

		} else if (geometry_array_type == 5) { // FRACTURE_CLUSTERS_ZONE_CLUSTERS
			setGeometryFaultZone(geometry_parameter_type, geometry_data_type,
					geometry_value, name, min_max_mode, faultArray);

		} else if (geometry_array_type == 6) { // FRACTURE_CLUSTERS_ZONE_FRACTURES
			setGeometrySingleFaultZone(geometry_parameter_type,
					geometry_data_type, geometry_value, name, min_max_mode,
					faultArray);

		} else if (geometry_array_type == 7) { // FRACTURE_CLUSTERS_CORRIDOR_CLUSTERS
			setGeometryCorridor(geometry_parameter_type, geometry_data_type,
					geometry_value, name, min_max_mode, faultArray);

		} else if (geometry_array_type == 8) { // FRACTURE_CLUSTERS_CORRIDOR_FRACTURE
			setGeometrySingleFaultCorridor(geometry_parameter_type,
					geometry_data_type, geometry_value, name, min_max_mode,
					&faultArray->fracturesSet_list_[0]);

		}
	}

}
/**
 * @fn void setDistributionSingleFault(int, int, double, std::vector<std::string>&, std::vector<double>&, FractureSet*)
 * @brief
 *
 * @param distribution_parameter_type
 * @param distribution_data_type
 * @param distribution_value
 * @param name
 * @param min_max_mode
 * @param faultSingle
 */
void RemFracResSim::setDistributionSingleFault(int distribution_parameter_type,
		int distribution_data_type, double distribution_value,
		std::vector<std::string> &name, std::vector<double> &min_max_mode,
		FractureSet *faultSingle) {

	if (distribution_parameter_type == 3) {
		if (std::isnan(distribution_value))
			distribution_value = 1e-6;

		if (distribution_data_type == 0) {
			faultSingle->density_ = distribution_value;
		} else if (distribution_data_type == 1) {
			faultSingle->density_distribution_name_ = name[0];

			if (name[1] == "NORMAL") {
				faultSingle->density_type_
						== FractureSet::distribution_type::normal;
			} else if (name[1] == "TRIANGULAR") {

				faultSingle->density_type_ =
						FractureSet::distribution_type::triangulated;
				faultSingle->minmaxmode_density_values_[2] = min_max_mode[1];

			} else if (name[1] == "UNIFORM") {
				faultSingle->density_type_ =
						FractureSet::distribution_type::uniform;
			} else if (name[1] == "LOGNORMAL") {
				faultSingle->density_type_ =
						FractureSet::distribution_type::lognormal;
			}

			faultSingle->minmaxmode_density_values_[0] = min_max_mode[0];
			faultSingle->minmaxmode_density_values_[1] = min_max_mode[2];

		} else {
			faultSingle->density_type_index_ = 2;
			std::string property_name =name[0];
			 std::transform(property_name.begin(), property_name.end(), property_name.begin(), my_toupper());

			faultSingle->density_property_name_ = property_name;
		}
	} else if (distribution_parameter_type == 11) {
		if (std::isnan(distribution_value))
			distribution_value = 0.1;

		if (distribution_data_type == 0) {
			faultSingle->distribution_shadow_zone_type_index_ = 0;
			faultSingle->default_stress_zone_max_ = distribution_value;
		} else if (distribution_data_type == 1) {
			faultSingle->width_stress_zone_distribution_name_ = name[0];

			faultSingle->distribution_shadow_zone_type_index_ = 1;
			if (name[1] == "NORMAL") {
				faultSingle->shadow_zone_type_
						== FractureSet::distribution_type::normal;
			} else if (name[1] == "TRIANGULAR") {

				faultSingle->shadow_zone_type_ =
						FractureSet::distribution_type::triangulated;
				faultSingle->minmaxmode_shadow_zone_values_[2] =
						min_max_mode[1];

			} else if (name[1] == "UNIFORM") {
				faultSingle->shadow_zone_type_ =
						FractureSet::distribution_type::uniform;
			} else if (name[1] == "LOGNORMAL") {
				faultSingle->shadow_zone_type_ =
						FractureSet::distribution_type::lognormal;
			}

			faultSingle->minmaxmode_shadow_zone_values_[0] = min_max_mode[0];
			faultSingle->minmaxmode_shadow_zone_values_[1] = min_max_mode[2];

		} else {
			faultSingle->distribution_shadow_zone_type_index_ = 2;
			std::string property_name =name[0];
			 std::transform(property_name.begin(), property_name.end(), property_name.begin(), my_toupper());

			faultSingle->width_stress_zone_property_name_ = property_name;
		}
	}
}
/**
 * @fn void setDistributionHangingWall(int, int, double, std::vector<std::string>&, std::vector<double>&, FaultArraySet*)
 * @brief
 *
 * @param distribution_parameter_type
 * @param distribution_data_type
 * @param distribution_value
 * @param name
 * @param min_max_mode
 * @param faultArray
 */
void RemFracResSim::setDistributionHangingWall(int distribution_parameter_type,
		int distribution_data_type, double distribution_value,
		std::vector<std::string> &name, std::vector<double> &min_max_mode,
		FaultArraySet *faultArray) {
	if (distribution_parameter_type == 2) {
		if (std::isnan(distribution_value))
			distribution_value = 1e-6;
		faultArray->distribution_damage_zone_hanging_index_ =
				distribution_data_type;
		if (distribution_data_type == 0) {
			faultArray->distr_damage_zone_hanging_val_ = distribution_value;

		} else if (distribution_data_type == 1) {
			faultArray->damage_zone_hanging_distribution_name_ = name[0];
			faultArray->distribution_damage_zone_hanging_index_ = 0;
			if (name[1] == "NORMAL") {
				faultArray->damage_zone_hanging_type_
						== FaultArraySet::distribution_type::normal;
			} else if (name[1] == "TRIANGULAR") {

				faultArray->damage_zone_hanging_type_ =
						FaultArraySet::distribution_type::triangulated;
				faultArray->minmaxmode_damage_zone_hanging_values_[2] =
						min_max_mode[1];

			} else if (name[1] == "UNIFORM") {
				faultArray->damage_zone_hanging_type_ =
						FaultArraySet::distribution_type::uniform;
			} else if (name[1] == "LOGNORMAL") {
				faultArray->damage_zone_hanging_type_ =
						FaultArraySet::distribution_type::lognormal;
			}

			faultArray->minmaxmode_damage_zone_hanging_values_[0] =
					min_max_mode[0];
			faultArray->minmaxmode_damage_zone_hanging_values_[1] =
					min_max_mode[2];

		} else {
			faultArray->distribution_damage_zone_hanging_index_ = 1;
			std::string property_name =name[0];
			 std::transform(property_name.begin(), property_name.end(), property_name.begin(), my_toupper());

			faultArray->damage_zone_hanging_property_name_ = property_name;
		}
	}
}
/**
 * @fn void setDistributionFootWall(int, int, double, std::vector<std::string>&, std::vector<double>&, FaultArraySet*)
 * @brief
 *
 * @param distribution_parameter_type
 * @param distribution_data_type
 * @param distribution_value
 * @param name
 * @param min_max_mode
 * @param faultArray
 */
void RemFracResSim::setDistributionFootWall(int distribution_parameter_type,
		int distribution_data_type, double distribution_value,
		std::vector<std::string> &name, std::vector<double> &min_max_mode,
		FaultArraySet *faultArray) {
	if (distribution_parameter_type == 2) {
		if (std::isnan(distribution_value))
			distribution_value = 1e-6;
		faultArray->distribution_damage_zone_footwall_index_ =
				distribution_data_type;
		if (distribution_data_type == 0) {
			faultArray->distr_damage_zone_footwall_val_ = distribution_value;
		} else if (distribution_data_type == 1) {
			faultArray->damage_zone_footwall_distribution_name_ = name[0];
			faultArray->distribution_damage_zone_footwall_index_ = 0;
			if (name[1] == "NORMAL") {
				faultArray->damage_zone_footwall_type_
						== FaultArraySet::distribution_type::normal;
			} else if (name[1] == "TRIANGULAR") {

				faultArray->damage_zone_footwall_type_ =
						FaultArraySet::distribution_type::triangulated;
				faultArray->minmaxmode_damage_zone_footwall_values_[2] =
						min_max_mode[1];

			} else if (name[1] == "UNIFORM") {
				faultArray->damage_zone_footwall_type_ =
						FaultArraySet::distribution_type::uniform;
			} else if (name[1] == "LOGNORMAL") {
				faultArray->damage_zone_footwall_type_ =
						FaultArraySet::distribution_type::lognormal;
			}

			faultArray->minmaxmode_damage_zone_footwall_values_[0] =
					min_max_mode[0];
			faultArray->minmaxmode_damage_zone_footwall_values_[1] =
					min_max_mode[2];

		} else {
			faultArray->distribution_damage_zone_footwall_index_ = 1;
			std::string property_name =name[0];
			 std::transform(property_name.begin(), property_name.end(), property_name.begin(), my_toupper());

			faultArray->damage_zone_footwall_property_name_ = property_name;
		}
	}
}
/**
 * @fn void setDistributionArray(int, int, double, std::vector<std::string>&, std::vector<double>&, FaultArraySet*)
 * @brief
 *
 * @param distribution_parameter_type
 * @param distribution_data_type
 * @param distribution_value
 * @param name
 * @param min_max_mode
 * @param faultArray
 */
void RemFracResSim::setDistributionArray(int distribution_parameter_type,
		int distribution_data_type, double distribution_value,
		std::vector<std::string> &name, std::vector<double> &min_max_mode,
		FaultArraySet *faultArray) {
	if (distribution_parameter_type == 5 || distribution_parameter_type == 1
			|| distribution_parameter_type == 0) {
		if (std::isnan(distribution_value))
			distribution_value = 1e-6;

		if (distribution_data_type == 0) {
			faultArray->density_ = distribution_value;
		} else if (distribution_data_type == 1) {
			faultArray->density_distribution_name_ = name[0];

			if (name[1] == "NORMAL") {
				faultArray->density_type_
						== FaultArraySet::distribution_type::normal;
			} else if (name[1] == "TRIANGULAR") {

				faultArray->density_type_ =
						FaultArraySet::distribution_type::triangulated;
				faultArray->minmaxmode_density_values_[2] = min_max_mode[1];

			} else if (name[1] == "UNIFORM") {
				faultArray->density_type_ =
						FaultArraySet::distribution_type::uniform;
			} else if (name[1] == "LOGNORMAL") {
				faultArray->density_type_ =
						FaultArraySet::distribution_type::lognormal;
			}

			faultArray->minmaxmode_density_values_[0] = min_max_mode[0];
			faultArray->minmaxmode_density_values_[1] = min_max_mode[2];

		} else {
			std::string property_name =name[0];
			 std::transform(property_name.begin(), property_name.end(), property_name.begin(), my_toupper());

			faultArray->density_property_name_ = property_name;
		}
	} else if (distribution_parameter_type == 4
			|| distribution_parameter_type == 12
			|| distribution_parameter_type == 13) {
		if (std::isnan(distribution_value))
			distribution_value = 1e-6;
		faultArray->minimum_space_between_array_law_type_ =
				distribution_data_type;
		if (distribution_data_type == 0) {
			faultArray->minimum_space_between_array_law_value_ =
					distribution_value;
		} else if (distribution_data_type == 1) {
			faultArray->minimum_space_between_array_distribution_name_ =
					name[0];
			faultArray->distribution_minimum_space_type_index_ = 0;
			if (name[1] == "NORMAL") {
				faultArray->minimum_space_type_
						== FaultArraySet::distribution_type::normal;
			} else if (name[1] == "TRIANGULAR") {

				faultArray->minimum_space_type_ =
						FaultArraySet::distribution_type::triangulated;
				faultArray->minmaxmode_minimum_space_values_[2] =
						min_max_mode[1];

			} else if (name[1] == "UNIFORM") {
				faultArray->minimum_space_type_ =
						FaultArraySet::distribution_type::uniform;
			} else if (name[1] == "LOGNORMAL") {
				faultArray->minimum_space_type_ =
						FaultArraySet::distribution_type::lognormal;
			}

			faultArray->minmaxmode_minimum_space_values_[0] = min_max_mode[0];
			faultArray->minmaxmode_minimum_space_values_[1] = min_max_mode[2];

		} else {
			faultArray->distribution_minimum_space_type_index_ = 1;
			std::string property_name =name[0];
			 std::transform(property_name.begin(), property_name.end(), property_name.begin(), my_toupper());
			faultArray->minimum_space_between_array_property_name_ = property_name;
		}
	}
}

/**
 * @fn void setDistributionArrayFault(int, int, double, std::vector<std::string>&, std::vector<double>&, FaultArraySet*)
 * @brief
 *
 * @param distribution_density_law_type
 * @param distribution_data_type
 * @param distribution_value
 * @param name
 * @param min_max_mode
 * @param faultArray
 */
void RemFracResSim::setDistributionArrayFault(int distribution_density_law_type,
		int distribution_data_type, double distribution_value,
		std::vector<std::string> &name, std::vector<double> &min_max_mode,
		FaultArraySet *faultArray) {

	faultArray->fracturesSet_list_[0].density_law_type_ =
			distribution_density_law_type;
	if (distribution_density_law_type == 0) {
		if (std::isnan(distribution_value))
			distribution_value = 5.0;
		faultArray->fracturesSet_list_[0].set_fractures_number(
				distribution_value);
	} else {
		if (distribution_data_type == 0) {
			if (std::isnan(distribution_value))
				distribution_value = 1e-6;
			faultArray->fracturesSet_list_[0].density_ = distribution_value;

		} else if (distribution_data_type == 1) {
			faultArray->fracturesSet_list_[0].density_distribution_name_ =
					name[0];
			if (name[1] == "NORMAL") {
				faultArray->fracturesSet_list_[0].density_type_
						== FractureSet::distribution_type::normal;
			} else if (name[1] == "TRIANGULAR") {

				faultArray->fracturesSet_list_[0].density_type_ =
						FractureSet::distribution_type::triangulated;
				faultArray->fracturesSet_list_[0].minmaxmode_density_values_[2] =
						min_max_mode[1];

			} else if (name[1] == "UNIFORM") {
				faultArray->fracturesSet_list_[0].density_type_ =
						FractureSet::distribution_type::uniform;
			} else if (name[1] == "LOGNORMAL") {
				faultArray->fracturesSet_list_[0].density_type_ =
						FractureSet::distribution_type::lognormal;
			}

			faultArray->fracturesSet_list_[0].minmaxmode_density_values_[0] =
					min_max_mode[0];
			faultArray->fracturesSet_list_[0].minmaxmode_density_values_[1] =
					min_max_mode[2];

		} else {
			std::string property_name =name[0];
			 std::transform(property_name.begin(), property_name.end(), property_name.begin(), my_toupper());

			faultArray->fracturesSet_list_[0].density_property_name_ = property_name;
		}
	}
}

/**
 * @fn void setDistributionArrayStressZoneFault(int, double, std::vector<std::string>&, std::vector<double>&, FaultArraySet*)
 * @brief
 *
 * @param distribution_data_type
 * @param distribution_value
 * @param name
 * @param min_max_mode
 * @param faultArray
 */
void RemFracResSim::setDistributionArrayStressZoneFault(
		int distribution_data_type, double distribution_value,
		std::vector<std::string> &name, std::vector<double> &min_max_mode,
		FaultArraySet *faultArray) {
	if (std::isnan(distribution_value))
		distribution_value = 0.5;
	faultArray->fracturesSet_list_[0].width_stress_zone_law_type_ =
			distribution_data_type;
	if (distribution_data_type == 0) {
		faultArray->fracturesSet_list_[0].default_stress_zone_max_ =
				distribution_value;
	} else if (distribution_data_type == 1) {
		faultArray->fracturesSet_list_[0].width_stress_zone_distribution_name_ =
				name[0];

		if (name[1] == "NORMAL") {
			faultArray->fracturesSet_list_[0].shadow_zone_type_
					== FractureSet::distribution_type::normal;
		} else if (name[1] == "TRIANGULAR") {

			faultArray->fracturesSet_list_[0].shadow_zone_type_ =
					FractureSet::distribution_type::triangulated;
			faultArray->fracturesSet_list_[0].minmaxmode_shadow_zone_values_[2] =
					min_max_mode[1];

		} else if (name[1] == "UNIFORM") {
			faultArray->fracturesSet_list_[0].shadow_zone_type_ =
					FractureSet::distribution_type::uniform;
		} else if (name[1] == "LOGNORMAL") {
			faultArray->fracturesSet_list_[0].shadow_zone_type_ =
					FractureSet::distribution_type::lognormal;
		}

		faultArray->fracturesSet_list_[0].minmaxmode_shadow_zone_values_[0] =
				min_max_mode[0];
		faultArray->fracturesSet_list_[0].minmaxmode_shadow_zone_values_[1] =
				min_max_mode[2];
	} else {
		std::string property_name =name[0];
		 std::transform(property_name.begin(), property_name.end(), property_name.begin(), my_toupper());

		faultArray->fracturesSet_list_[0].width_stress_zone_property_name_ =
				property_name;
	}
}
/**
 * @fn void setDistributionArrayStressZoneFault(int, double, std::vector<std::string>&, std::vector<double>&, FaultArraySet*)
 * @brief
 *
 * @param distribution_data_type
 * @param distribution_value
 * @param name
 * @param min_max_mode
 * @param faultArray
 */
void RemFracResSim::setDistributionArrayEchelonSpaceFault(
		int distribution_data_type, double distribution_value,
		std::vector<std::string> &name, std::vector<double> &min_max_mode,
		FaultArraySet *faultArray) {
	if (std::isnan(distribution_value))
		distribution_value = 0.5;
	faultArray->fracturesSet_list_[0].distribution_echelon_type_index_ =
			distribution_data_type;
	if (distribution_data_type == 0) {
		faultArray->fracturesSet_list_[0].echelon_space_ = distribution_value;
	} else if (distribution_data_type == 1) {
		faultArray->fracturesSet_list_[0].echelon_space_distribution_name_ =
				name[0];

		if (name[1] == "NORMAL") {
			faultArray->fracturesSet_list_[0].echelon_space_type_
					== FractureSet::distribution_type::normal;
		} else if (name[1] == "TRIANGULAR") {

			faultArray->fracturesSet_list_[0].echelon_space_type_ =
					FractureSet::distribution_type::triangulated;
			faultArray->fracturesSet_list_[0].minmaxmode_echelon_space_values_[2] =
					min_max_mode[1];

		} else if (name[1] == "UNIFORM") {
			faultArray->fracturesSet_list_[0].echelon_space_type_ =
					FractureSet::distribution_type::uniform;
		} else if (name[1] == "LOGNORMAL") {
			faultArray->fracturesSet_list_[0].echelon_space_type_ =
					FractureSet::distribution_type::lognormal;
		}

		faultArray->fracturesSet_list_[0].minmaxmode_echelon_space_values_[0] =
				min_max_mode[0];
		faultArray->fracturesSet_list_[0].minmaxmode_echelon_space_values_[1] =
				min_max_mode[2];
	} else {
		faultArray->fracturesSet_list_[0].echelon_space_distribution_name_ =
				name[0];
	}
}
/**
 * @fn void setDistributionArrayFaultFamili1(int, int, double, std::vector<std::string>&, std::vector<double>&, FaultArraySet*)
 * @brief
 *
 * @param distribution_density_law_type
 * @param distribution_data_type
 * @param distribution_value
 * @param name
 * @param min_max_mode
 * @param faultArray
 */
void RemFracResSim::setDistributionArrayFaultFamili1(
		int distribution_density_law_type, int distribution_data_type,
		double distribution_value, std::vector<std::string> &name,
		std::vector<double> &min_max_mode, FaultArraySet *faultArray) {

	faultArray->fracturesSet_list_[0].density_law_type_ =
			distribution_density_law_type;
	if (distribution_density_law_type == 0) {
		if (std::isnan(distribution_value))
			distribution_value = 5.0;

		faultArray->fracturesSet_list_[0].set_fractures_number(
				distribution_value);
	} else {
		if (std::isnan(distribution_value))
			distribution_value = 1e-6;
		faultArray->fracturesSet_list_[0].density_type_index_ =
				distribution_data_type;
		if (distribution_data_type == 0) {
			faultArray->fracturesSet_list_[0].density_ = distribution_value;
		} else if (distribution_data_type == 1) {

			faultArray->fracturesSet_list_[0].density_distribution_name_ =
					name[0];
			if (name[1] == "NORMAL") {
				faultArray->fracturesSet_list_[0].density_type_
						== FractureSet::distribution_type::normal;
			} else if (name[1] == "TRIANGULAR") {

				faultArray->fracturesSet_list_[0].density_type_ =
						FractureSet::distribution_type::triangulated;
				faultArray->fracturesSet_list_[0].minmaxmode_density_values_[2] =
						min_max_mode[1];

			} else if (name[1] == "UNIFORM") {
				faultArray->fracturesSet_list_[0].density_type_ =
						FractureSet::distribution_type::uniform;
			} else if (name[1] == "LOGNORMAL") {
				faultArray->fracturesSet_list_[0].density_type_ =
						FractureSet::distribution_type::lognormal;
			}

			faultArray->fracturesSet_list_[0].minmaxmode_density_values_[0] =
					min_max_mode[0];
			faultArray->fracturesSet_list_[0].minmaxmode_density_values_[1] =
					min_max_mode[2];

		} else {
			std::string property_name =name[0];
			 std::transform(property_name.begin(), property_name.end(), property_name.begin(), my_toupper());

			faultArray->fracturesSet_list_[0].density_property_name_ = property_name;
		}
	}
}
/**
 * @fn void setDistributionArray(int, int, int, int, double, std::vector<std::string>&, std::vector<double>&, FaultArraySet*)
 * @brief
 *
 * @param distribution_density_law_type
 * @param distribution_array_type
 * @param distribution_parameter_type
 * @param distribution_data_type
 * @param distribution_value
 * @param name
 * @param min_mode_max
 * @param faultArray
 */
void RemFracResSim::setDistributionArray(int distribution_density_law_type,
		int distribution_array_type, int distribution_parameter_type,
		int distribution_data_type, double distribution_value,
		std::vector<std::string> &name, std::vector<double> &min_mode_max,
		FaultArraySet *faultArray) {
	if (distribution_array_type == 1) {   // fault array
		setDistributionArray(distribution_parameter_type,
				distribution_data_type, distribution_value, name, min_mode_max,
				faultArray);
	} else if (distribution_array_type == 2) { // fault array segment

		if (distribution_parameter_type == 6) {   // fault segment density
			setDistributionArrayFault(distribution_density_law_type,
					distribution_data_type, distribution_value, name,
					min_mode_max, faultArray);
		} else if (distribution_parameter_type == 11) { // stress shadow zone
			setDistributionArrayStressZoneFault(distribution_data_type,
					distribution_value, name, min_mode_max, faultArray);
		} else if (distribution_parameter_type == 11
				&& faultArray->typeFamilyArray_
						== FaultArraySet::anastomosing) {
			setDistributionArrayFaultFamili1(distribution_density_law_type,
					distribution_data_type, distribution_value, name,
					min_mode_max, faultArray);
		} else if (distribution_parameter_type == 11
				&& faultArray->typeFamilyArray_
						== FaultArraySet::anastomosing) {
			setDistributionArrayFaultFamili1(distribution_density_law_type,
					distribution_data_type, distribution_value, name,
					min_mode_max, faultArray);
		}
	} else if (distribution_array_type == 3) {  //fault zone clusters
		setDistributionArray(distribution_parameter_type,
				distribution_data_type, distribution_value, name, min_mode_max,
				faultArray);
	} else if (distribution_array_type == 4) {  // fault zone structural
	} else if (distribution_array_type == 5) {   // fault zone hanging
		setDistributionHangingWall(distribution_density_law_type,
				distribution_data_type, distribution_value, name, min_mode_max,
				faultArray);
	} else if (distribution_array_type == 6) {   // fault zone footwall
		setDistributionFootWall(distribution_density_law_type,
				distribution_data_type, distribution_value, name, min_mode_max,
				faultArray);
	} else if (distribution_array_type == 7) {   // fault zone fractures
	} else if (distribution_array_type == 8) {   // corridor clusters
		setDistributionArrayFaultFamili1(distribution_density_law_type,
				distribution_data_type, distribution_value, name, min_mode_max,
				faultArray);
	} else if (distribution_array_type == 9) {   // corridor fractures
		setDistributionArray(distribution_parameter_type,
				distribution_data_type, distribution_value, name, min_mode_max,
				faultArray);
	}
}

/**
 * @fn void manageTotalDistribution(std::fstream&, FaultArraySet*, FractureSet*)
 * @brief
 *
 * @param fin
 * @param faultArray
 * @param faultSingle
 */
void RemFracResSim::manageTotalDistribution(std::fstream &fin,
		FaultArraySet *faultArray, FractureSet *faultSingle) {
	int distribution_parameter_type = -1;
	int distribution_array_type = -1;
	int distribution_data_type = -1;
	int distribution_density_law_type = -1;
	double distribution_value = -1.0;
	std::vector<double> min_max_mode;
	std::vector<std::string> name = manageDistribution(fin,
			distribution_parameter_type, distribution_array_type,
			distribution_data_type, distribution_density_law_type,
			distribution_value, min_max_mode);
	//if(distribution_density_law_type < 0) distribution_density_law_type =0;

	if (faultSingle != NULL && distribution_array_type == 0) {
		if (distribution_density_law_type == 0) {
			if (std::isnan(distribution_value))
				distribution_value = 5.0;

			faultSingle->set_fractures_number(distribution_value);
			if (distribution_parameter_type == 3)
				faultSingle->density_law_type_ = distribution_density_law_type;
		} else {

			setDistributionSingleFault(distribution_parameter_type,
					distribution_data_type, distribution_value, name,
					min_max_mode, faultSingle);
			if (distribution_parameter_type == 3)
				faultSingle->density_law_type_ = distribution_density_law_type;
		}
	} else if (faultArray != NULL) {
		if (distribution_density_law_type == 0) {
			if (std::isnan(distribution_value))
				distribution_value = 5.0;

			if ((distribution_array_type == 1
					&& distribution_parameter_type == 5)
					|| (distribution_array_type == 3
							&& (distribution_parameter_type == 16
									|| distribution_parameter_type == 1))
					|| (distribution_array_type == 7
							&& distribution_parameter_type == 16)) {
				faultArray->density_law_type_ = distribution_density_law_type;
				faultArray->set_fault_array_number(distribution_value);
			} else if (distribution_parameter_type == 10
					&& (distribution_array_type == 8)) {
				faultArray->fracturesSet_list_[0].density_law_type_ =
						distribution_density_law_type;
				faultArray->fracturesSet_list_[0].set_fractures_number(
						distribution_value);
			} else if (distribution_parameter_type == 2
					&& (distribution_array_type == 5)) {
				faultArray->fracturesSet_list_[0].density_hanging_wall_law_type_ =
						distribution_density_law_type;
				faultArray->fracturesSet_list_[0].set_fractures_number_hanging_wall(
						distribution_value);
			} else if (distribution_parameter_type == 2
					&& (distribution_array_type == 6)) {
				faultArray->fracturesSet_list_[0].density_foot_wall_law_type_ =
						distribution_density_law_type;
				faultArray->fracturesSet_list_[0].set_fractures_number_foot_wall(
						distribution_value);
			} else if (distribution_array_type == 9
					&& distribution_parameter_type == 0) {
				faultArray->density_law_type_ = distribution_density_law_type;
				faultArray->set_fault_array_number(distribution_value);
			} else if (distribution_array_type == 2
					&& distribution_parameter_type == 6) {
				faultArray->fracturesSet_list_[0].density_law_type_ =
						distribution_density_law_type;
				faultArray->fracturesSet_list_[0].set_fractures_number(
						distribution_value);
			}
		} else {
			if ((distribution_array_type == 1)
					|| (distribution_array_type == 3
							&& (distribution_parameter_type == 0
									|| distribution_parameter_type == 1
									|| distribution_parameter_type == 14))
					|| (distribution_array_type == 7
							&& (distribution_parameter_type == 0
									|| distribution_parameter_type == 1
									|| distribution_parameter_type == 16))) {
				//faultArray->density_law_type_ = distribution_density_law_type;
				setDistributionArray(distribution_density_law_type,
						distribution_array_type, distribution_parameter_type,
						distribution_data_type, distribution_value, name,
						min_max_mode, faultArray);
			} else if (((distribution_parameter_type == 1
					|| distribution_parameter_type == 12)
					&& distribution_array_type == 3)) {
				setDistributionArray(distribution_density_law_type,
						distribution_array_type, distribution_parameter_type,
						distribution_data_type, distribution_value, name,
						min_max_mode, faultArray);

			} else if (((distribution_parameter_type == 2)
					&& (distribution_array_type == 5
							|| distribution_array_type == 6))) {
				setDistributionArray(distribution_density_law_type,
						distribution_array_type, distribution_parameter_type,
						distribution_data_type, distribution_value, name,
						min_max_mode, faultArray);

			} else if (((distribution_parameter_type == 13)
					&& (distribution_array_type == 9))) {
				setDistributionArray(distribution_density_law_type,
						distribution_array_type, distribution_parameter_type,
						distribution_data_type, distribution_value, name,
						min_max_mode, faultArray);

			} else if ((distribution_parameter_type == 6)
					&& (distribution_array_type == 2)) {
				setDistributionArrayFault(distribution_density_law_type,
						distribution_data_type, distribution_value, name,
						min_max_mode, faultArray);

			} else if ((distribution_parameter_type == 11)
					&& (distribution_array_type == 2)) {
				setDistributionArrayStressZoneFault(distribution_data_type,
						distribution_value, name, min_max_mode, faultArray);
			} else if ((distribution_parameter_type == 7)
					&& (distribution_array_type == 2)) {
				//distance between segement echelon case
				setDistributionArrayEchelonSpaceFault(distribution_data_type,
						distribution_value, name, min_max_mode, faultArray);
			}

		}
	}
}

/**
 * @fn fractures_intersect::FracResIntersectionTypeEnum manageIntersection(std::fstream&, int&)
 * @brief
 *
 * @param fin
 * @param index
 * @return
 */

fractures_intersect::FracResIntersectionTypeEnum RemFracResSim::manageIntersection(
		std::fstream &fin, int &index) {
	std::string contenu_ligne;
	std::string comment;
	std::string separator;
	std::istringstream my_stream;
	fractures_intersect::FracResIntersectionTypeEnum tab_type_intersection;

	while (getline(fin, contenu_ligne)) // tant que l'on peut mettre la ligne dans "contenu"
	{

		cout << contenu_ligne << endl;  // on l'affiche
		my_stream.str(contenu_ligne);
		std::vector<std::string> tokenVector = splitLine(my_stream, ';');

		if (tokenVector[0] == fracresFileKeyWord_->intersection_type) {
			tab_type_intersection = getIntersectionType(tokenVector[1]);
		} else if (tokenVector[0]
				== fracresFileKeyWord_->intersection_fracture_type) {

			index = getIntersectionCategoryType(tokenVector[1]);
		} else if (tokenVector[0] == fracresFileKeyWord_->intersection_end) {

			break;
		}
		my_stream.clear();
	}

	return tab_type_intersection;
}

/**
 * @fn std::vector<std::string> manageDistribution(std::fstream&, int&, int&, int&, int&, double&, std::vector<double>&)
 * @brief
 *
 * @param fin
 * @param distribution_parameter_type
 * @param distribution_array_type
 * @param distribution_data_type
 * @param distribution_density_law_type
 * @param distribution_value
 * @param minmax
 * @return
 */
std::vector<std::string> RemFracResSim::manageDistribution(std::fstream &fin,
		int &distribution_parameter_type, int &distribution_array_type,
		int &distribution_data_type, int &distribution_density_law_type,
		double &distribution_value, std::vector<double> &minmax) {
	std::string contenu_ligne;
	std::string comment;
	std::string separator;

	std::istringstream my_stream;

	std::vector<std::string> nameList_;
	minmax.resize(3, 0);
	minmax[2] = 1.0;

	while (getline(fin, contenu_ligne)) // tant que l'on peut mettre la ligne dans "contenu"
	{

		cout << contenu_ligne << endl;  // on l'affiche
		my_stream.str(contenu_ligne);
		// To store the stream string
		std::vector<std::string> tokenVector = splitLine(my_stream, ';');

		if (tokenVector[0]
				== fracresFileKeyWord_->distribution_parameter_type) {
			distribution_parameter_type = getDistributionType(tokenVector[1]);
		} else if (tokenVector[0]
				== fracresFileKeyWord_->distribution_array_type) {
			distribution_array_type = getDistributionArrayType(tokenVector[1]);
		} else if (tokenVector[0]
				== fracresFileKeyWord_->distribution_data_type) {
			distribution_data_type = getDataType(tokenVector[1]);
		} else if (tokenVector[0] == fracresFileKeyWord_->distribution_unit) {
			distribution_density_law_type = getDensityLawType(tokenVector[1]);
		} else if (tokenVector[0] == fracresFileKeyWord_->distribution_value) {
			distribution_value = std::strtod(tokenVector[1].c_str(), NULL);
		} else if (tokenVector[0]
				== fracresFileKeyWord_->distribution_property_name) {
			nameList_.push_back(tokenVector[1]);
			nameList_.push_back("None");
		} else if (tokenVector[0]
				== fracresFileKeyWord_->distribution_distribution_str) {

			std::istringstream my_stream2;
			my_stream2.str(tokenVector[1]);
			std::vector<std::string> tokenVectorVirgule = splitLine(my_stream2,
					',');
			nameList_.push_back(tokenVectorVirgule[0]);
			nameList_.push_back(tokenVectorVirgule[1]);

			std::istringstream my_stream3;
			my_stream3.str(tokenVectorVirgule[2]);
			std::vector<std::string> tokenVectorValues = splitLine(my_stream3,
					':');

			minmax[0] = std::strtod(tokenVectorValues[1].c_str(), NULL);
			minmax[1] = std::strtod(tokenVectorValues[3].c_str(), NULL);
			minmax[2] = std::strtod(tokenVectorValues[5].c_str(), NULL);

		} else if (tokenVector[0] == fracresFileKeyWord_->distribution_end) {

			break;
		}
		my_stream.clear();
	}
	return nameList_;
}

/**
 * @fn std::vector<std::string> manageGeometry(std::fstream&, int&, int&, int&, double&, std::vector<double>&)
 * @brief
 *
 * @param fin
 * @param geometry_parameter_type
 * @param geometry_array_type
 * @param geometry_data_type
 * @param distribution_value
 * @param min_max_mode
 * @return
 */
std::vector<std::string> RemFracResSim::manageGeometry(std::fstream &fin,
		int &geometry_parameter_type, int &geometry_array_type,
		int &geometry_data_type, double &distribution_value,
		std::vector<double> &min_max_mode) {
	std::string contenu_ligne;
	std::string comment;
	std::string separator;
	std::istringstream my_stream;

	std::vector<std::string> name;
	name.resize(2, "None");
	min_max_mode.resize(3, 0.0);
	min_max_mode[2] = 1.0;

	while (getline(fin, contenu_ligne)) // tant que l'on peut mettre la ligne dans "contenu"
	{

		cout << contenu_ligne << endl;  // on l'affiche
		my_stream.str(contenu_ligne);

		std::vector<std::string> tokenVector = splitLine(my_stream, ';');

		if (tokenVector[0] == fracresFileKeyWord_->geometry_parameter_type) {
			geometry_parameter_type = getGeometryType(tokenVector[1]);
		} else if (tokenVector[0] == fracresFileKeyWord_->geometry_array_type) {
			geometry_array_type = getGeometryArrayType(tokenVector[1]);
		} else if (tokenVector[0] == fracresFileKeyWord_->geometry_data_type) {
			geometry_data_type = getDataType(tokenVector[1]);

			if (geometry_data_type == 3)
				name[1] = "NORMAL";
		} else if (tokenVector[0] == fracresFileKeyWord_->geometry_unit) {
			std::string unit = tokenVector[1];
		} else if (tokenVector[0] == fracresFileKeyWord_->geometry_value) {
			distribution_value = std::strtod(tokenVector[1].c_str(), NULL);
		} else if (tokenVector[0]
				== fracresFileKeyWord_->geometry_property_name) {
			name[0] = tokenVector[1];

		} else if (tokenVector[0]
				== fracresFileKeyWord_->geometry_distribution_str) {
			std::istringstream my_stream2;
			my_stream2.str(tokenVector[1]);
			std::vector<std::string> tokenVectorVirgule = splitLine(my_stream2,
					',');
			name[0] = (tokenVectorVirgule[0]);
			name[1] = (tokenVectorVirgule[1]);

			if (name[1] == "POWER") {
				std::istringstream my_stream3;
				my_stream3.str(tokenVectorVirgule[2]);
				std::vector<std::string> tokenVectorValues = splitLine(
						my_stream3, ':');
				min_max_mode[0] = std::strtod(tokenVectorValues[1].c_str(),
						NULL);
				my_stream3.clear();
				my_stream3.str(tokenVectorVirgule[3]);
				std::vector<std::string> tokenVectorValues2 = splitLine(
						my_stream3, ':');
				min_max_mode[2] = std::strtod(tokenVectorValues2[1].c_str(),
						NULL);
			} else {
				std::istringstream my_stream3;
				my_stream3.str(tokenVectorVirgule[2]);
				std::vector<std::string> tokenVectorValuesMin = splitLine(
						my_stream3, ':');
				min_max_mode[0] = std::strtod(tokenVectorValuesMin[1].c_str(),
						NULL);

				if(tokenVectorVirgule.size() > 4){
				my_stream3.clear();
				my_stream3.str(tokenVectorVirgule[3]);
				std::vector<std::string> tokenVectorValuesMode = splitLine(
										my_stream3, ':');
				min_max_mode[1] = std::strtod(tokenVectorValuesMode[1].c_str(),
						NULL);
				my_stream3.clear();
				my_stream3.str(tokenVectorVirgule[4]);
				std::vector<std::string> tokenVectorValuesMax = splitLine(
										my_stream3, ':');
				min_max_mode[2] = std::strtod(tokenVectorValuesMax[1].c_str(),
						NULL);
				}else{
					my_stream3.clear();
					my_stream3.str(tokenVectorVirgule[3]);
					std::vector<std::string> tokenVectorValuesMax = splitLine(
											my_stream3, ':');
					min_max_mode[2] = std::strtod(tokenVectorValuesMax[1].c_str(),
							NULL);
				}
			}

		} else if (tokenVector[0]
				== fracresFileKeyWord_->geometry_normal_distribution_min) {
			min_max_mode[0] = std::strtod(tokenVector[1].c_str(), NULL);
		} else if (tokenVector[0]
				== fracresFileKeyWord_->geometry_normal_distribution_max) {
			min_max_mode[2] = std::strtod(tokenVector[1].c_str(), NULL);
		} else if (tokenVector[0] == fracresFileKeyWord_->geometry_end) {

			break;
		}
		my_stream.clear();
	}
	return name;
}
/**
 * @fn void manageBedParallel(std::fstream&, FracResStage&, FracResHypothesis&)
 * @brief
 *
 * @param fin
 * @param fracresStage
 * @param fracResHypothesis
 */
void RemFracResSim::manageBedParallel(std::fstream &fin,
		FracResStage &fracresStage, FracResHypothesis &fracResHypothesis) {
	std::string contenu_ligne;
	std::string comment;
	std::string separator;
	std::istringstream my_stream;
	FracResBeddingProperties *fracResBeddingProperties =
			new FracResBeddingProperties();
	if (fracResBeddingProperties->min_mode_max_.size() <= 0) {
		fracResBeddingProperties->min_mode_max_.resize(3, 0.0);
		fracResBeddingProperties->min_mode_max_[1] = 1.0;
	}
	while (getline(fin, contenu_ligne)) // tant que l'on peut mettre la ligne dans "contenu"
	{
		cout << contenu_ligne << endl;  // on l'affiche
		my_stream.str(contenu_ligne);

		std::vector<std::string> tokenVector = splitLine(my_stream, ';');

		if (tokenVector[0]
				== fracresFileKeyWord_->bed_parrallel_slope_parameter) {
			fracResBeddingProperties->slope_parameter_name_ = tokenVector[1];

		} else if (tokenVector[0]
				== fracresFileKeyWord_->bed_parrallel_slope_data_type) {
			fracResBeddingProperties->slope_data_type_ =
					(FracResDataTypeEnum) getDataType(tokenVector[1]);

		} else if (tokenVector[0]
				== fracresFileKeyWord_->bed_parrallel_slope_unit) {
			fracResBeddingProperties->slope_unit_name_ = tokenVector[1];
		} else if (tokenVector[0]
				== fracresFileKeyWord_->bed_parrallel_slope_value) {
			fracResBeddingProperties->slope_data_value_ = std::strtod(
					tokenVector[1].c_str(), NULL);

		} else if (tokenVector[0]
				== fracresFileKeyWord_->bed_parrallel_prop_discrete_prop_name) {
			std::string property_name =tokenVector[1];
			 std::transform(property_name.begin(), property_name.end(), property_name.begin(), my_toupper());

			fracResBeddingProperties->bed_parrallel_discrete_property_name_ =
					property_name;
		} else if (tokenVector[0]
				== fracresFileKeyWord_->bed_parrallel_is_region) {
			fracResBeddingProperties->is_region_ = (tokenVector[1] == "true");
		} else if (tokenVector[0]
				== fracresFileKeyWord_->bed_parrallel_region_name) {
			if (fracResBeddingProperties->is_region_)
				fracResBeddingProperties->region_name_ = tokenVector[1];

		} else if (tokenVector[0]
				== fracresFileKeyWord_->bed_parrallel_slope_ratio) {
			fracResBeddingProperties->slope_ratio_data_value_ = std::strtod(
					tokenVector[1].c_str(), NULL);

		} else if (tokenVector[0]
				== fracresFileKeyWord_->bed_parrallel_slope_distribution_str) {
			fracResBeddingProperties->bed_parrallel_slope_distribution_name_ =
					tokenVector[1];
			std::istringstream my_stream2;
			my_stream2.str(tokenVector[1]);
			std::vector<std::string> tokenVectorVirgule = splitLine(my_stream2,
					',');
			fracResBeddingProperties->bed_parrallel_slope_distribution_name_ =
					(tokenVectorVirgule[0]);
			fracResBeddingProperties->bed_parrallel_slope_distribution_type_name_ =
					(tokenVectorVirgule[1]);

			std::istringstream my_stream3;
			my_stream3.str(tokenVectorVirgule[2]);
			std::vector<std::string> tokenVectorValues = splitLine(my_stream3,
					':');

			fracResBeddingProperties->min_mode_max_[0] = std::strtod(
					tokenVectorValues[1].c_str(), NULL);
			fracResBeddingProperties->min_mode_max_[2] = std::strtod(
					tokenVectorValues[3].c_str(), NULL);
			fracResBeddingProperties->min_mode_max_[1] = std::strtod(
					tokenVectorValues[5].c_str(), NULL);

		} else if (tokenVector[0]
				== fracresFileKeyWord_->bed_parrallel_slope_prop_name) {
			fracResBeddingProperties->bed_parrallel_slope_prop_name_ =
					tokenVector[1];
		} else if (tokenVector[0] == fracresFileKeyWord_->bed_aperture_begin) {
			FracResBeddingAperture *fracResBeddingAperture_ =
					manageBeddingAperture(fin);
			fracResBeddingProperties->bedding_list_.push_back(
					fracResBeddingAperture_);
		} else if (tokenVector[0] == fracresFileKeyWord_->bed_parrallel_end) {

			fracResHypothesis.bed_parallel_List_.push_back(
					fracResBeddingProperties);
			fracResHypothesis.type_dfn_.push_back(3);
			break;
		}

		my_stream.clear();

	}

}

/**
 * @fn std::vector<common_solver::Point3D> write_vertices(const geode::StructuralModel&)
 * @brief write vertices on a common-solver::point data structure for global mesh
 *
 * @pre
 * @post
 * @param model
 * @return
 */
void RemFracResSim::write_vertices_geode(const geode::StructuralModel &model) {

	static constexpr geode::index_t OFFSET_START { 1 };

	auto nb_vertices = model.nb_unique_vertices();
	grid_vertices_.reserve(nb_vertices);
	std::vector<geode::index_t> atoms;
	std::vector<absl::flat_hash_map<geode::ComponentMeshVertex, geode::index_t> > vertices_;
	vertices_.resize(nb_vertices);
	std::vector<geode::index_t> atoms_index;
	std::cout << "nb unique vertex = " << nb_vertices << "";
	std::cout << std::endl;
	for (const auto v : geode::Range { nb_vertices }) {
		const auto block_vertices = model.component_mesh_vertices(v,
				geode::Block3D::component_type_static());
		auto &mcvs = vertices_[v];
		mcvs.reserve(block_vertices.size());
		// std::cout <<" block_vertice size = " << block_vertices.size() << "" ;
		//	std::cout <<  std::endl;
		geode::index_t id { 0 };
		auto first_uuid = block_vertices.front().component_id.id();

		//std::cout <<" first_uuid = " << first_uuid.string() << "" ;
		//	std::cout <<  std::endl;

		for (const auto i : geode::Range { 1, block_vertices.size() }) {
			const auto &uuid = block_vertices[i].component_id.id();
			if (uuid < first_uuid) {
				first_uuid = uuid;
				id = i;
			}
		}

		const auto first = id;

		if (std::find(atoms_index.begin(), atoms_index.end(), v + OFFSET_START)
				== atoms_index.end())
			atoms_index.push_back(v + OFFSET_START);

		const auto &first_vertex = block_vertices[first];
		mcvs.emplace(first_vertex, v + OFFSET_START);

		const auto &block = model.block(first_vertex.component_id.id());
		const auto &mesh = block.mesh<geode::TetrahedralSolid3D>();
		// std::cout << "VRTX " << count++ << " "<< mesh.point( first_vertex.vertex ).string() ;
		//		std::cout <<  std::endl;

		geode::Point3D geode_pt = mesh.point(first_vertex.vertex);

		grid_vertices_.push_back(geode_pt);
		for (const auto i : geode::Indices { block_vertices }) {
			if (i == first) {
				continue;
			}
			atoms.push_back(v + OFFSET_START);
			mcvs.emplace(block_vertices[i], nb_vertices + OFFSET_START);

			nb_vertices++;
		}
	}

	std::cout << "vertices size " << grid_vertices_.size() << " ";
	std::cout << std::endl;

}
/**
 * @fn void save_vtk_tetra(std::string, std::string)
 * @brief
 *
 * @param filename
 * @param path
 */

void RemFracResSim::save_vtk_tetra(std::string filename, std::string path) {

	const auto output_file_native_vtk = absl::StrCat(path, filename, ".vtk");
	bool success = save_meshed_vtk_grid_tetra(output_file_native_vtk);

}
//-----------------------------------------------------
//save the meshed fractures data
//-----------------------------------------------------
/**
 * @fn bool save_meshed_vtk_grid_tetra(std::string)
 * @brief
 *
 * @param filename
 * @return
 */
bool RemFracResSim::save_meshed_vtk_grid_tetra(std::string filename) {
	std::fstream os(filename.data(), std::ios_base::out);

	if (os.is_open() == false) {
		std::cout << "Could not save vtk the mesh plan file" << filename
				<< std::endl;
		return false;
	}

	os << "# vtk DataFile Version 2.0" << std::endl;
	os << "Cube faille" << std::endl;
	os << "ASCII" << std::endl;
	os << "DATASET UNSTRUCTURED_GRID" << std::endl;

	os << "POINTS " << grid_vertices_.size() << " double" << std::endl;
	for (int i = 0; i < grid_vertices_.size(); i++) {
		os << grid_vertices_[i].value(0) << " " << grid_vertices_[i].value(1)
				<< " " << grid_vertices_[i].value(2) << std::endl;
	}

	os << "CELLS " << tetra_component_.size() << " "
			<< tetra_component_.size() * 5 << std::endl;

	for (int i = 0; i < tetra_component_.size(); i++) {
		os << "4 ";
		for (int j = 0; j < 4; j++)
			os << tetra_component_[i]->t_nodes_[j] << " ";
		os << std::endl;
	}

	os << "CELL_TYPES " << tetra_component_.size() << std::endl;
	for (int i = 0; i < tetra_component_.size(); i++) {
		os << "10 " << std::endl;
	}

	os << "CELL_DATA " << tetra_component_.size() << std::endl;
	os << "SCALARS cell_index float 1" << std::endl;
	os << "LOOKUP_TABLE default" << std::endl;

	for (int k = 0; k < tetra_component_.size(); k++) {
		float val = tetra_component_[k]->active_index_ * 1.0 + 1.0;
		os << val << " ";
	}
	os << std::endl;

	os << "SCALARS block_index float 1" << std::endl;
	os << "LOOKUP_TABLE default" << std::endl;

	for (int l = 0; l < tetra_component_.size(); l++) {
		float val = tetra_component_[l]->geode_block_index_ * 1.0 + 1.0;
		os << val << " ";
	}
	os << std::endl;

	os.close();
	return true;

}

/**
 * @fn void initialize_block_tetra(const geode::StructuralModel&)
 * @brief
 *
 * @param model
 */
void RemFracResSim::initialize_block_tetra(
		const geode::StructuralModel &model) {

	int indice_block = 0;
	int cell_index = 0;
	//block_global_index_[indice_block]=cell_index;
	for (const auto &block : model.blocks()) {
		const auto name = block.name();
		const auto &id = block.component_id();
		const auto &mesh = block.mesh<geode::TetrahedralSolid3D>();
		mesh.enable_edges();
		const auto attr_edge_from =
				mesh.edges().edge_attribute_manager().find_or_create_attribute<
						geode::VariableAttribute, double>("edge_id", ndvalue_);

		for (const auto p : geode::Range { mesh.nb_polyhedra() }) {

			const bool test_border = mesh.is_polyhedron_on_border(p);
			int n1 = model.unique_vertex(
					{ id, mesh.polyhedron_vertex( { p, 0 }) });
			int n2 = model.unique_vertex(
					{ id, mesh.polyhedron_vertex( { p, 1 }) });
			int n3 = model.unique_vertex(
					{ id, mesh.polyhedron_vertex( { p, 2 }) });
			int n4 = model.unique_vertex(
					{ id, mesh.polyhedron_vertex( { p, 3 }) });

			//std::cout <<" "<< n1 <<" "<< n2 << " "<< n3 <<" "<<n4;
			//std::cout <<  std::endl;

			fractures_intersect::tetra *tet = new fractures_intersect::tetra(n1,
					n2, n3, n4, p, cell_index);

			tet->border_ = test_border;
			tet->t_nodes_[0] = n1;
			tet->geode_block_index_ = indice_block;
			tet->geode_polyhedron_id_ = p;
			tet->active_index_ = cell_index;

			tetra_component_.push_back(tet);

			cell_index++;

		}
		indice_block++;

	}

	std::cout << " tetra component size " << tetra_component_.size();
	std::cout << std::endl;

}


int  RemFracResSim::runsimulation(std::string fileParameter, std::string fileStrm){


	try {

		geode::OpenGeodeModelLibrary::initialize();
		geode::OpenGeodeMeshLibrary::initialize();
		geode::GeosciencesExplicitLibrary::initialize();
		geode::GeosciencesIOModelLibrary::initialize();
		geode::GeosciencesIOMeshLibrary::initialize();
		geode::IOModelLibrary::initialize();

		const geode::StructuralModel &model = geode::load_structural_model(
				fileStrm);

		remfracres::RemFracResSim obj;

		obj.boundingBoxes_ = obj.create_bounding_box(model);

		obj.tree_ = std::shared_ptr < geode::AABBTree3D > (new geode::AABBTree3D(obj.boundingBoxes_));

		obj.disteval_= std::shared_ptr < BoxAABBEvalDistance3D > (new BoxAABBEvalDistance3D(obj.boundingBoxes_));


		std::string path = remfracres::data_path;

		if (model.nb_blocks() > 1) {
			obj.write_vertices_geode(model);
			obj.initialize_block_tetra(model);
		}


		size_t posFirst = fileParameter.find_last_of(".");

		std::string	pathParam = fileParameter.substr(0, posFirst + 1);
		obj.openFracres(fileParameter);

		absl::flat_hash_map<geode::uuid, std::vector<geode::uuid> > map_beddings_to_horizons =
				obj.map_beddings_from_horizons(model);
		absl::flat_hash_map<geode::uuid, std::vector<geode::uuid> > map_mechanical_to_horizons =
				obj.map_mechanical_from_horizons(model);

		//if (!obj.use_global_outpout_directory_path_)path = obj.outpout_directory_path_;


		 posFirst = path.find_last_of("/");

		path = path.substr(0, posFirst + 1) + "Result_" + pathParam;

		int result = mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);

		if (result == 0) {
			std::cout << "Répertoire " << path <<" have been created" << std::endl;
		} else {
			rmdir(path.c_str());
			mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
			std::cout << "Répertoire " << path <<" is still exist." << std::endl;
		}

		path +="/";
		obj.outpout_directory_path_ = path;

		obj.update_fracture_set_outpout_directory_path(path);

		float seed = obj.seed_;
		obj.stage_hypothesis_fractures_set_generation(model, false);

		//obj.fractures_set_generation(model,false);

		obj.simulate_intersection(true);


		std::string namefile_save1 = absl::StrCat("LS2", "_");
		std::string namefile_save = absl::StrCat(namefile_save1, obj.name_);

		std::string tetra_filename = absl::StrCat("tetra_", namefile_save);
		//obj.save(namefile_save, path);
		obj.saveFromStage(namefile_save, path);
		obj.save_vtk_tetra(tetra_filename, path);
		if (obj.success_) {

			Logger::info("TEST SUCCESS");
		} else {
			Logger::info("TEST FAILED - Check Fractures generation");
		}

		return 0;
	} catch (...) {
		return geode::geode_lippincott();
	}



}


void RemFracResSim::create_main_directory_output(std::string &path) {
	int result = mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
	if (result == 0) {
		std::cout << "Répertoire " << path << " have been created" << std::endl;
	} else {
		std::string pathnew = "rm -rf " + path ;

		system(pathnew.c_str());

		mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
		std::cout << "Répertoire " << path << " have been removed and recreated." << std::endl;

	}
}

void RemFracResSim::create_local_directory_output(std::string &path) {
	int result = mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
	if (result == 0) {
		std::cout << "Répertoire " << path << " have been created" << std::endl;
	} else {

		std::cout << "Répertoire " << path << " is created yet." << std::endl;

	}
}

void RemFracResSim::create_result_directory(const std::string &fileParameter, size_t posFirst,
		const std::string &namefile, std::string &path) {
	path = fileParameter.substr(0, posFirst + 1) + "RESULT_" + namefile;
	create_main_directory_output(path);
	path += "/";
}

void RemFracResSim::simulationAllScenario(const geode::StructuralModel &model, std::string nameFile) {
	std::cout<<""<<endl;
	std::cout << "" << endl;

	fractures_set_generation_for_scenario(model,nameFile, false);
	std::cout<<""<<endl;
	std::cout << "" << endl;

}


void RemFracResSim::simulationAllStage(const geode::StructuralModel &model){
	std::cout<<""<<endl;
	std::cout << "" << endl;

	std::chrono::time_point < std::chrono::system_clock > total_StartTime;
	std::chrono::time_point < std::chrono::system_clock > total_EndTime;

	total_StartTime = std::chrono::system_clock::now();
	stage_hypothesis_fractures_set_generation(model, false);
	simulate_intersection(true);
	total_EndTime = std::chrono::system_clock::now();
	double write_miliseconde = std::chrono::duration_cast< std::chrono::milliseconds> (total_EndTime - total_StartTime).count();
	total_simulation_time_ = write_miliseconde ;
	cout << " Temps total de simulation" << endl;
	cout << " nb miliseconds " << write_miliseconde << " nb second "<< total_simulation_time_ << '\n';
	cout << endl;

	total_StartTime = std::chrono::system_clock::now();

	std::string namefile_save = absl::StrCat(name_, "_");
	std::string tetra_filename = absl::StrCat("tetra_", namefile_save);
	//obj.save(namefile_save, path);
	saveFromStage(namefile_save, outpout_directory_path_);
	save_vtk_tetra(tetra_filename, outpout_directory_path_);
	total_EndTime = std::chrono::system_clock::now();
	write_miliseconde = std::chrono::duration_cast < std::chrono::milliseconds> (total_EndTime - total_StartTime).count();
	total_writting_output_time_ = write_miliseconde;
	cout << " Temps total d ecriture" << endl;
	cout << " nb miliseconds " << write_miliseconde << " nb second "<< total_writting_output_time_ << '\n';
	cout << endl;
	std::cout<<""<<endl;
	std::cout << "" << endl;



}
std::string RemFracResSim::getPathAndFileName(const std::string &filename,std::string &namefile) {
	size_t posFirst0 = filename.find_last_of("/");
	std::string path = filename.substr(0, posFirst0 + 1);
	size_t posLast = filename.find_last_of(".");
	namefile = filename.substr(posFirst0 + 1, posLast - posFirst0);
	std::cout << "filename1test= " << filename;
	std::cout << std::endl;
	std::cout << " namefile= " << namefile << " :: " << path << " :: "
			<< posFirst0 << " :: " << posLast << " :: " << filename.length();
	std::cout << std::endl;
	return path;
}

size_t RemFracResSim::getParameterSuffixeName(const std::string &fileParameter,std::string &namefile,std::string& path ) {
	size_t posFirst = fileParameter.find_last_of("/");
	path = fileParameter.substr(0, posFirst + 1);
	size_t posLast = fileParameter.find_last_of(".");
	namefile.clear();
	namefile = fileParameter.substr(posFirst + 1, posLast-posFirst-1);


	return posFirst;
}

void RemFracResSim::run_fracturing_simulation(const geode::StructuralModel &model,const std::string &fileParameter) {

	boundingBoxes_ = create_bounding_box(model);
	tree_ = std::shared_ptr < geode::AABBTree3D> (new geode::AABBTree3D(boundingBoxes_));
	disteval_ = std::shared_ptr < BoxAABBEvalDistance3D> (new BoxAABBEvalDistance3D(boundingBoxes_));
	if (model.nb_blocks() > 1) {
		write_vertices_geode(model);
		initialize_block_tetra(model);
	}else{
		 std::string tetra_filename = absl::StrCat(outpout_directory_path_,"Grille");
		 absl::string_view filename_without_ext{ tetra_filename };
		 const auto output_filename_vtu =absl::StrCat( filename_without_ext, "_output.vtu" );

		for (const auto &block : model.blocks()){
			const auto name = block.name();
			const auto& mesh = block.mesh<geode::TetrahedralSolid3D>();

			geode::save_tetrahedral_solid(mesh, output_filename_vtu);
		}


	}
	absl::flat_hash_map<geode::uuid, std::vector<geode::uuid> > map_beddings_to_horizons = map_beddings_from_horizons(model);
	absl::flat_hash_map<geode::uuid, std::vector<geode::uuid> > map_mechanical_to_horizons = map_mechanical_from_horizons(model);

	openFracres(fileParameter);

	by_scenario_ = fracResScenarioList_.size() > 0;

	if(by_scenario_){
		std::map<std::pair<int,int>,std::string > mapEvaluateFracture = fractures_number_evaluation_for_scenario(model);
		affiche_fracture_number_evaluation(mapEvaluateFracture);
		simulationAllScenario(model,outpout_file_name_);
	}else{
		getEvaluateFracturesNumberPerAllStage(model);
			getEvaluateFracturesNumberMessage();
		simulationAllStage(model);
	}

	getSuccessMessage();

}

int RemFracResSim::getEvaluateFracturesNumberPerStageHypothesis(const geode::StructuralModel &model,const FracResStage &stage,const FracResHypothesis &fracResHypothesis) {

	int evaluate_fracture_number =0;
	int indice_bed_interface = 0;
	int indice_bed_parallel = 0;
	int indice_single_fault = 0;
	int indice_fault_array = 0;
	int indice_fault_corridor = 0;
	int indice_fault_zone = 0;

	for (int type_dfn : fracResHypothesis.type_dfn_) {

		switch (type_dfn) {
		case 0: {
			evaluate_fracture_number += fracResHypothesis.beddingList_.size();
			break;
		}
		case 1: {
			evaluate_fracture_number += get_evaluate_fracture_single_faults_from_hypothesis(fracResHypothesis,indice_single_fault, stage, model);
			indice_single_fault++;
			break;
		}
		case 2: {
			evaluate_fracture_number += get_evaluate_fracture_fault_array_from_hypothesis(fracResHypothesis,indice_fault_array, stage, model);
			indice_fault_array++;
			break;
		}
		case 3: {
			evaluate_fracture_number += fracResHypothesis.bed_parallel_List_.size();
			break;
		}
		case 4: {
			evaluate_fracture_number += get_evaluate_fracture_cluster_corridor_from_hypothesis(fracResHypothesis, indice_fault_corridor, stage,model);
			indice_fault_corridor++;
			break;
		}
		case 5: {
			evaluate_fracture_number += get_evaluate_fracture_cluster_fault_zone_from_hypothesis(fracResHypothesis, indice_fault_zone, stage,  model);
			indice_fault_zone++;
			break;
		}
		default: {
			break;
		}
		} // switch
	}
	return evaluate_fracture_number;
	// loop type_dfn}
}
 void RemFracResSim::getEvaluateFracturesNumberPerAllStage(const geode::StructuralModel &model){

	set_model_fractures_bounding_box(model);
	for (FracResStage& stage : fracResStageList_) {

		std::cout << " -------------- Stage N°" << stage.stageIndex_<< " named " << stage.name_ << "------------------" << endl;

		for (FracResHypothesis& fracResHypothesis : stage.hypotheisisList_) {
			std::cout << " -------------- Hypothesis N°" << fracResHypothesis.hypothesisIndex_<< " named " << fracResHypothesis.name_ << "------------------" << endl;
			int nb_fractures = getEvaluateFracturesNumberPerStageHypothesis(model,stage,fracResHypothesis);
			fracResHypothesis.evaluate_fractures_number_total_ = nb_fractures;
			stage.evaluate_fractures_number_total_per_hypothesis_.push_back(nb_fractures);
		}

	}
}

 void RemFracResSim::getEvaluateFracturesNumberMessage(){


 	for (FracResStage stage : fracResStageList_) {

 		if(stage.evaluate_fractures_number_total_per_hypothesis_.size() <=0) continue;
 		int nb_hyp=0;
 		for (FracResHypothesis fracResHypothesis : stage.hypotheisisList_) {

 			int number_of_fracture = stage.evaluate_fractures_number_total_per_hypothesis_[nb_hyp];

 			Logger::info(" Evaluated number of fractures for stage " + stage.name_  + " hypothesis " + fracResHypothesis.name_ + " : " + std::to_string(number_of_fracture));

 			nb_hyp++;
 		}

 	}
 }


void RemFracResSim::getTimeMessage() {
	if (by_scenario_) {
		for (FracResScenario scenario : fracResScenarioList_) {
			Logger::info("Scenario : " + scenario.name_);
			if(!scenario.modeling_) continue;
			for (int ireal = 0; ireal < scenario.numberOfRealisation_;
					ireal++) {
				Logger::info(
						"Total simulation and writting time [s] for realization "
								+ std::to_string(ireal + 1) + ": "
								+ getFormatTime(
										scenario.total_simulation_time_list_[ireal]));

				for(int i=0; i < scenario.success_list_[ireal].size(); i++){
					success_ = scenario.success_list_[ireal][i];
				}
			}


		}
	} else {
		Logger::info(
				"Total simulation time [s] : "
						+ getFormatTime(total_simulation_time_));
		Logger::info(
				"Total writting time [s] : "
						+ getFormatTime(total_writting_output_time_));
	}
}

void RemFracResSim::getSuccessMessage() {
	getTimeMessage();
	if (success_) {
		Logger::info("TEST SUCCESS");
	} else {
		Logger::info("TEST FAILED - Check Fractures generation");
	}


}

std::string  RemFracResSim::getFormatTime( int time_second){


	int  time_days = time_second / 86400;
	time_second %=24;

	int time_hours= time_second / 3600;
	time_second %=3600;

	int time_minutes= time_second / 60;
	time_second %=60;

	int time_seconds = time_second/ 60;

	time_second %=60;

	int time_miliseconds = time_second;



	std::string result = std::to_string(time_days) + " j: " + std::to_string(time_hours) + " h: " + std::to_string(time_minutes) + " m: " + std::to_string(time_seconds)  + " s: " + std::to_string(time_miliseconds)  + " ms";

	return result;
}

}



/**
 * @fn int main(int, char**)
 * @brief
 *
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char **argv) {
	using namespace geode;
	using namespace std;
	using namespace remfracres;
	using namespace fractures_intersect;
	//using std::filesystem::current_path;

	if (argc < 2) {

		try {

			std::string filename = absl::StrCat(absl::StrCat(remfracres::data_path,"LS2.lso"));
			filename = absl::StrCat(absl::StrCat(remfracres::data_path,"Test_FRACRES_output.og_strm"));

			//filename = "fast_grid_properties_new.";
			std::string prefixe_fracture = "FractureSet_";
			std::string fractureName = "b1_f1_fracture_arrays_relay.dat";
			std::string nameFracShadow =
					"LS2_parameter_fractures_FRACRES_v5_WO_ShadowZone_TEST.txt";
			std::string nameFracShadow2 =
					"LS2_parameter_fractures_FRACRES_v5_With_ShadowZone_TEST.txt";

			std::string surface_filename = absl::StrCat(prefixe_fracture,filename);
			std::string filenameParameter = absl::StrCat(absl::StrCat(remfracres::data_path, nameFracShadow2));
			//std::string model_filename = "reference_grid_test_output.og_strm";
			std::string nameFracParam = "FracResSettingsTestV1_Antoine.txt";
			nameFracParam = "FRACRES_parameter_set2.txt";
			//nameFracParam = "FRACRES_v1102_SingleFaults_TestEL01_CstValues_Bis.txt";
			nameFracParam = "FRACRES_v1003_SingleFaults_TestEL01_DipAzimuth_DistribUniform.txt";
			nameFracParam = "FRACRES_v1003_SingleFaults_TestEL01_GeometryDistrib_WO_Correlation.txt";
			nameFracParam ="FRACRES_v1103_SingleFaults_TestEL04_Height_DistribLogNormal.txt";
			nameFracParam ="FRACRES_v1106_2Sets_SingleFaults_TestEL01b_Abutting_With_ShadowZone_DistribGeometricParameters.txt";
			nameFracParam ="FRACRES_v1106_ArraysClusters_MultiScenarios_TestRN01.txt";

			nameFracParam = "FRACRES_v1108_REM2520_TestEL01.txt";
			nameFracParam = "FRACRES_v1109_REM2532_TestEL01.txt";
			nameFracParam = "FRACRES_v1109_REM2532_TestEL02.txt";

			nameFracParam ="FRACRES_v1106_REM2395_TestEL02.txt";
			nameFracParam = "Test_ra_1234_V2.txt";

			nameFracParam = "FRACRES_parameter_set_V1110_AL_01.txt";
			nameFracParam ="FRACRES_48523_REM2058_2102_TestEL01_run.txt";
			nameFracParam = "FRACRES_V1112_TestProperty_TestEL01.txt";
			nameFracParam = "FRACRES_V1112_DensityVariable_TestEL01.txt";
			nameFracParam = "FRACRES_V1112_TestProperty1_TestEL01.txt";

			std::string fileParameter = absl::StrCat(absl::StrCat(remfracres::data_path, nameFracParam));
			remfracres::RemFracResSim obj;
			geode::OpenGeodeModelLibrary::initialize();
			geode::OpenGeodeMeshLibrary::initialize();
			geode::GeosciencesExplicitLibrary::initialize();
			geode::GeosciencesIOModelLibrary::initialize();
			geode::GeosciencesIOMeshLibrary::initialize();
			geode::IOModelLibrary::initialize();
			const geode::StructuralModel &model = geode::load_structural_model(filename);
			obj.get_property_attribute_names_from_block_model(model);
			std::string namefile;
			std::string path = obj.getPathAndFileName(filename, namefile);
			size_t posFirst = obj.getParameterSuffixeName(fileParameter, namefile,path);
			obj.create_result_directory(fileParameter, posFirst, namefile, path);
			obj.outpout_directory_path_ = path;
			obj.outpout_file_name_ = namefile;

			obj.run_fracturing_simulation(model, fileParameter);



			return 0;
		} catch (...) {
			return geode::geode_lippincott();
		}
	} else {


		std::string filename = argv[1];

		std::string fileParameter = argv[2];


		remfracres::RemFracResSim obj;
		geode::OpenGeodeModelLibrary::initialize();
		geode::OpenGeodeMeshLibrary::initialize();
		geode::GeosciencesExplicitLibrary::initialize();
		geode::GeosciencesIOModelLibrary::initialize();
		geode::GeosciencesIOMeshLibrary::initialize();
		geode::IOModelLibrary::initialize();
		const geode::StructuralModel &model = geode::load_structural_model(filename);
		obj.get_property_attribute_names_from_block_model(model);
		std::string namefile;
		std::string path = obj.getPathAndFileName(filename, namefile);
		size_t posFirst = obj.getParameterSuffixeName(fileParameter, namefile,path);
		obj.create_result_directory(fileParameter, posFirst, namefile, path);
		obj.outpout_directory_path_ = path;
		obj.outpout_file_name_ = namefile;
		obj.run_fracturing_simulation(model, fileParameter);


	}
	return 0;
}

