//----------------------------  remfracres_sim.h  ---------------------------
// 
//
//  Description:	remfracres_sim class use for build fractures sets
//					information
// 
//  Documents:		RemFracResSim.h
//
//  Creation date  : 20/08/2022
//
//  Author:			Rabah Namar
//
//  Copyright (C) 2022 by in2earth
//
//---------------------------- remfracres_sim.h  ---------------------------
#pragma once

#ifndef RemFracResSim_h
#define RemFracResSim_h

// local includes
#include <remfracres/FracResStage.h>
#include <remfracres/FracResScenario.h>
#include <remfracres/FracResUnit.h>
#include <remfracres/RemFracResFileManagement.h>
#include <remfracres/CellId.h>
#include <remfracres/FacetId.h>
#include <remfracres/BoxAABBEvalDistance.h>

//#include <remfracres/FaultArraySet.h>
#include <remfracres/fracResFileKeyWord.h>
#include <primitive_mesh/tetra.h>
#include <numeric/model.h>
// standard includes
#include <iomanip>
#include <vector>
#include <array>
#include <list>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <utility>
#include <chrono>
#include <ctime>
#include <cctype>
#include <algorithm>
// boost includes 
#include <boost/random.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <Eigen/Core>
#include <Eigen/Geometry>

// geode include
#include <geode/geometry/basic_objects/plane.h>
#include <geode/geometry/basic_objects/segment.h>
#include <geode/geometry/basic_objects/triangle.h>
#include <geode/geometry/point.h>
#include <geode/geometry/bounding_box.h>
#include <geode/io/mesh/common.h>
#include <geode/io/model/common.h>
#include <geode/model/mixin/core/block.h>
#include <geode/model/mixin/core/model_boundary.h>
#include <geode/model/mixin/core/surface.h>
#include <geode/geosciences/explicit/representation/core/structural_model.h>
#include <geode/geosciences/explicit/representation/io/structural_model_output.h>
#include <geode/geosciences/explicit/representation/io/structural_model_input.h>
#include <geode/geosciences/explicit/representation/builder/structural_model_builder.h>
#include <geode/geosciences/explicit/representation/core/structural_model.h>
#include <geode/geosciences/explicit/common.h>
#include <geode/geosciences_io/mesh/common.h>
#include <geode/geosciences_io/model/common.h>


#include <geode/basic/types.h>
#include <geode/basic/uuid.h>
#include <geode/geometry/aabb.h>
#include <geode/geometry/distance.h>
#include <geode/mesh/core/solid_mesh.h>
#include <geode/mesh/core/tetrahedral_solid.h>
#include <geode/geometry/basic_objects/triangle.h>
#include <geode/geometry/basic_objects/tetrahedron.h>
#include <geode/model/mixin/core/block.h>
#include <geode/model/mixin/core/model_boundary.h>
#include <geode/model/mixin/core/surface.h>
#include <geode/basic/assert.h>
#include <geode/basic/attribute_manager.h>
#include <geode/basic/filename.h>
#include <geode/basic/logger.h>
#include <geode/mesh/builder/solid_edges_builder.h>
#include <geode/mesh/builder/solid_facets_builder.h>
#include <geode/mesh/core/geode/geode_tetrahedral_solid.h>
#include <geode/mesh/core/solid_edges.h>
#include <geode/mesh/core/solid_facets.h>
#include <geode/mesh/core/surface_mesh.h>
#include <geode/model/helpers/convert_model_meshes.h>
#include <geode/model/representation/core/brep.h>
#include <geode/model/helpers/convert_to_mesh.h>
#include <geode/model/helpers/component_mesh_polygons.h>
#include <geode/model/helpers/component_mesh_edges.h>
#include <geode/model/helpers/component_mesh_vertices.h>

#include <geode/mesh/core/geode/geode_triangulated_surface.h>
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
#include <geode/mesh/builder/surface_edges_builder.h>

#include <geode/mesh/core/surface_edges.h>
#include <geode/basic/assert.h>
#include <geode/basic/logger.h>
#include <absl/container/flat_hash_map.h>
#include <absl/container/flat_hash_set.h>
#include <absl/strings/str_split.h>
#include <absl/strings/string_view.h>



namespace remfracres{

/**
 *    FractureSet contains the parameters for the simulation 
 *    of a Boolean model and the generated fractures
 * 
 *    @author Pascal Siegel - Rabah Namar - Olivier Jaquet 
 *    @date 2022
 */
 
 class RemFracResSim
{


	 struct my_tolower
	 {
	     char operator()(char c) const
	     {
	         return std::tolower(static_cast<unsigned char>(c));
	     }
	 };

	 struct my_toupper
	 {
	     char operator()(char c) const
	     {
	         return std::toupper(static_cast<unsigned char>(c));
	     }
	 };


// variables

public:


	const double ndvalue_ = -9999.0;


	int nb_failles_;

	/**
	 *      name of RemFracResSim
	 */

	std::string name_;

	/**
	 *      fracture set number
	 */

	int nb_fracture_sets_;


	/**
	 *     stage number
	 */

	int nb_stage_;

	/**
	 *
	 */
	bool manage_bed_interface_;

	std::vector< FracResStage > fracResStageList_;

	std::vector< FracResScenario > fracResScenarioList_;


	FracResUnit* fracResUnit_;

	absl::flat_hash_map< std::string,geode::uuid >  geode_BlocName_to_BlockId_;

	absl::flat_hash_map< std::string,bool >  geode_BlocName_to_Density_exist_;

	absl::flat_hash_map< std::string,double >  geode_BlocName_to_volume_;

	absl::flat_hash_map< std::string,double >  geode_BlocName_to_area_;

	absl::flat_hash_map< std::string,double >  geode_BlocName_to_density_;
	
	/**
	 *      type method for fractures sets generation
	 *      0 : cross (boolean) 1: Abbutting 2 : Anastomose
	 */
	
	int method_type_;
	
	/**
	 *      generation domaine limits 
	 *      (grid limits + fracture max length/2)
	 */

	geode::BoundingBox3D generation_box_;

	/**
	 *      volume of the grid where the generation take place 
	 */

	float volume_;

	/**
	 *		area of the grid where the generation take place 
	 */

	float area_;

	/**
	 *		model fractures set
	 */
     fractures_intersect::model model_fractures_;

 	/**
 	 *		bounding box model
 	 */

     std::shared_ptr< fractures_intersect::Cboite_rect> bbox_;



 	/**
 	 *		list of  fractures set element
 	 */
     std::vector< FractureSet > fracturesSet_list_;

 	/**
  	 *		boolean for outpout dircetory same as input file
  	 */
      bool use_global_outpout_directory_path_;

      /**
      */
      std::string outpout_directory_path_;

      /**
       */
       std::string outpout_file_name_;


      /**
       *
      */
      bool intersect_inside_set_;

      /**
       *
      */

      bool success_;



      /**
       *
      */
      int cylinder_step_;


  	/**
  	 *      general seed value for intersection and fracture set
  	 */

  	int seed_;


      /**
  	 *    vector of surface name for discontinuity aperture
  	 *
  	 */

      std::map< std::string, std::string > beddings_name_;

		/**
		* default aperture values
		*/
      std::map< std::string,double> beddings_aperture_default_value_;

      /**
        *
       */
       std::unordered_map<  std::string, bool > beddings_surface_is_active_;


       std::unordered_map<  std::string, std::vector< fractures_intersect::FracResIntersectionTypeEnum > >beddings_type_intersection_;

       /**
    	 *    vector of surface name for discontinuity aperture
    	 *
    	 */

        std::map< std::string,std::string > mechanical_name_;

  		/**
  		* default aperture values
  		*/
        std::map< std::string,double> mechanical_aperture_default_value_;

        /**
          *
         */
         std::unordered_map<  std::string, bool > mechanical_surface_is_active_;


         std::unordered_map<  std::string, std::vector< fractures_intersect::FracResIntersectionTypeEnum > > mechanical_type_intersection_;


         /**
      	  * map for horizon ou fault  name (string) -> vecteur of geode::uuid
      	  */
           std::unordered_map< std::string,std::vector< geode::uuid > > mechanical_discontinuity_map_to_surface_id_;

           /**
            *
           */
           std::unordered_map<  geode::uuid , bool > mechanical_discontinuity_map_to_surface_id_is_active_;


           /**
        	  * map for boundary name (string) -> vecteur of geode::uuid
        	  */
          std::unordered_map< std::string,std::vector< geode::uuid > > mechanical_border_map_to_surface_id_;








      /**
   	  * map for horizon ou fault  name (string) -> vecteur of geode::uuid
   	  */
        std::unordered_map< std::string,std::vector< geode::uuid > > discontinuity_map_to_surface_id_;

        /**
         *
        */
        std::unordered_map<  geode::uuid , bool > discontinuity_map_to_surface_id_is_active_;


        /**
     	  * map for boundary name (string) -> vecteur of geode::uuid
     	  */
       std::unordered_map< std::string,std::vector< geode::uuid > > border_map_to_surface_id_;


       /**
     	 *    vector of surface name for boundary condition inflow
     	 *
     	 */

         std::vector< std::pair<std::string,std::string>  > surface_boundaries_name_in_;


         /**
          * vector of tetra on global mesh
          */
         std::vector< fractures_intersect::tetra *> tetra_component_;


         std::vector< geode::Point3D > grid_vertices_;



         std::vector< std::shared_ptr<fractures_intersect::fault_arrays_set_geo > >fault_array_zone_set_list_;


         std::shared_ptr< remfracres::FracResBeddingAperture > fracResBeddingzero_;
      	/**
      	 *		list of  fractures set element
      	 */
          std::vector< FaultArraySet > faultArraySet_list_;

         int nb_fault_arrays_set_;

         std::unordered_map<  std::string, std::vector< fractures_intersect::FracResIntersectionTypeEnum > > fault_array_type_intersection_;



         std::map< std::pair<int,int>, std::shared_ptr<fractures_intersect::fracture_set_geo> > stress_geo_list_;

         std::map< std::pair<int,int>,string > bed_name_vect_;

         std::map< std::pair<int,int>,std::shared_ptr<fractures_intersect::fracture_set_geo> > stress_geo_parallel_list_;
         std::map< std::pair<int,int>,string > bed_name_parallel_vect_;
         std::map< std::pair<int,int>,bool > active_discrete_pair_array_;

		FracResFileKeyWord*  fracresFileKeyWord_;


		std::string bed_interface_property_name_;


		RemFracResFileManagement* fileManagement_;

		int nb_total_bed_;


		std::vector<CellId> cells_;

		std::vector<geode::BoundingBox3D> boundingBoxes_;


		std::shared_ptr< geode::AABBTree3D > tree_;

		std::shared_ptr< BoxAABBEvalDistance3D > disteval_;

		int nb_realisation_;

		bool by_scenario_;


		double total_simulation_time_;
		double total_writting_output_time_;
// functions

public: 

	/**
	 *      Default constructor 
	 */

     RemFracResSim();

	/**
	 *		Copy constructor
	 * 
	 *      @param from The value to copy to this object
	 */

     RemFracResSim(const RemFracResSim& from);

	/**
	 *		Destructor 
	 */
        
	~RemFracResSim();

	/** 
     *      Assignment operator 
	 *
     *      @param from the value to assign to this object 
     * 
     *      @return A reference to this object 
     */
  
	RemFracResSim& operator=(const RemFracResSim& from);

	/** 
     *      copy function 
	 *
     *      @param from the value to assign to this object 
     */
  
	void copy(const RemFracResSim& from);
	


	/**
	 *      set_generation_box: set the domain limits
	 *                       
	 *		@param box: 3D box, delimiting the domain 
	 */
	
	void set_generation_box(geode::BoundingBox3D box);

	/**
	 * @fn void set_model_fractures_bounding_box(const geode::StructuralModel&)
	 * @brief
	 *
	 * @param model
	 */

	void set_model_fractures_bounding_box(const geode::StructuralModel& model);

	/**
	 *      fractures_generation: generation of the fractures 
	 */
	
	void fractures_set_generation(const geode::StructuralModel& model,bool estim_volume = false);


	/**
	 *      save the fracture set data
	 *      for the karstsim 
	 *
	 *		@param out: ostreamb 
     * 
	 */

     bool save(std::string filename, std::string path);

	 /**
	 *      open a fracture set data   
	 *      for the karstsim 
	 *
	 *		@param in: istreamb 
	 *
	 *		@param catalog: pointer to an AsciiCatalog
	 */

      bool open(std::string filename);

	/**
	 * @fn absl::flat_hash_map<geode::uuid,std::vector<geode::uuid>> map_surface_from_horizons(const geode::StructuralModel&)
	 * @brief
	 *
	 * @param model
	 * @return
	 */

      absl::flat_hash_map< geode::uuid,std::vector<geode::uuid> > map_beddings_from_horizons(const geode::StructuralModel &model);
	/**
	 * @fn void initialize_volume_and_area_from_global_grid_of_structural_model(const geode::StructuralModel&)
	 * @brief
	 *
	 * @param model
	 */
      void initialize_volume_and_area_from_global_grid_of_structural_model(const geode::StructuralModel& model);

		/**
		 * @fn void simulate_intersection(bool)
		 * @brief
		 *
		 * @param active
		 */
      void simulate_intersection(bool active);

      /**
       * @fn void update_fracture_set_outpout_directory_path(std::string)
       * @brief
       *
       * @param path
       */

      void update_fracture_set_outpout_directory_path(std::string path);
      /**
       * @fn int initialize_beddings_discontinuity(const geode::StructuralModel&)
       * @brief
       *
       * @param model
       * @return
       */
      void initialize_beddings_discontinuity(const geode::StructuralModel &model);


      /**
       * @fn void initalize_mechanical_discontinuity(const geode::StructuralModel&)
       * @brief
       *
       * @param model
       */
      void initalize_mechanical_discontinuity(const geode::StructuralModel &model);


  	/**
  	 * @fn absl::flat_hash_map<geode::uuid,std::vector<geode::uuid>> map_surface_from_horizons(const geode::StructuralModel&)
  	 * @brief
  	 *
  	 * @param model
  	 * @return
  	 */
      absl::flat_hash_map< geode::uuid,std::vector<geode::uuid> > map_mechanical_from_horizons(const geode::StructuralModel &model);


      std::string enumCategoryToString(fractures_intersect::FracResIntersectionCategoryEnum cat);

      std::string enumIntersectionToString(fractures_intersect::FracResIntersectionTypeEnum cat);


      /**
       * @fn bool save_vtk(std::string)
       * @brief
       *
       * @param filename
       * @return
       */

      	void save_vtk_tetra(std::string filename, std::string path);
      /**
       *
       * @fn std::vector<geode::Point3D> write_vertices_geode(const geode::StructuralModel&)
       * @brief
       *
       * @param model
       * @return
       */


      	void  write_vertices_geode(const geode::StructuralModel &model);

      	/**
      	 * @fn bool save_meshed_vtk_grid_tetra(std::string, std::vector<geode::Point3D>&, std::vector<std::vector<int>>&, std::vector<std::vector<float>>&)
      	 * @brief
      	 *
      	 * @param filename
      	 * @param tab_point
      	 * @param tab_element
      	 * @param tab_element_property
      	 * @return
      	 */
      		bool save_meshed_vtk_grid_tetra(std::string filename);

   /**
    * @fn void initialize_block_tetra(const geode::StructuralModel&)
    * @brief
    *
    * @param model
    */

      	    void initialize_block_tetra(const geode::StructuralModel &model);




	/**
	*      open a fracture set data
	*      for the karstsim
	*
	*		@param in: istreamb
	*
	*		@param catalog: pointer to an AsciiCatalog
	*/

	bool openFracres(std::string filename);
/**
 * @fn bool manageLine(std::fstream&)
 * @brief
 *
 * @param fin
 * @return
 */
	bool manageLine(std::fstream& fin);
/**
 * @fn void manageUnit(std::fstream&)
 * @brief
 *
 * @param fin
 */
	void manageUnit(std::fstream& fin);
/**
 * @fn FracResStage manageStage(std::fstream&)
 * @brief
 *
 * @param fin
 * @return
 */
	FracResStage manageStage(std::fstream& fin);
/**
 * @fn FracResHypothesis manageRepHypothesis(std::fstream&)
 * @brief
 *
 * @param fin
 * @return
 */
	FracResHypothesis  manageRepHypothesis(std::fstream& fin);

	FracResScenario manageScenario(std::fstream& fin);
	void manageRunParameter(std::fstream& fin);

	void manageMatrixProp(std::fstream& fin);
	void manageOutputs(std::fstream& fin);

    void manageStageHypScenarioKeyVal(std::fstream& fin);
	void manageUpscalingRep(std::fstream& fin);
	void manageRepresentationObject(std::fstream& fin);
	FracResFracture manageFracture(std::fstream& fin);

	FracResHypothesis manageHypothesis(std::fstream& fin);

	void manageBeddingStageHypothesis(std::fstream& fin, FracResHypothesis& fracResHypothesis );
	FracResBeddingAperture* manageBeddingAperture(std::fstream& fin);
	void manageClassicStageHypothesis(std::fstream& fin, FracResStage& fracresStage, FracResHypothesis& fracResHypothesis );
	void manageFractureProperty(std::fstream& fin, FracResStage& fracresStage, FracResHypothesis& fracResHypothesis);
	fractures_intersect::FracResIntersectionTypeEnum  manageIntersection(std::fstream& fin, int& index);
	std::vector< std::string > manageDistribution(std::fstream& fin,int& distribution_parameter_type,int& distribution_array_type, int& distribution_data_type, int& distribution_density_law_type,double& distribution_value,std::vector<double> &min_max_mode);
	std::vector< std::string >  manageGeometry(std::fstream& fin, int& distribution_parameter_type,int& distribution_array_type, int& distribution_data_type,double& distribution_value,std::vector< double>& min_max_mode);
	void manageBedParallel(std::fstream& fin, FracResStage& fracresStage, FracResHypothesis& fracResHypothesis);
	void split(const string& chaine,std::vector< std::pair<int,int> >& elements);
	void manageTotalDistribution(std::fstream& fin, FaultArraySet* faultArray,FractureSet* faultSingle);
	void manageTotalGeometry(std::fstream& fin, FaultArraySet* faultArray,FractureSet* faultSingle);

	void setGeometryDipAzimuth(int geometry_data_type,double geometry_value, std::vector< std::string >& name,const std::vector<double> &min_max_mode, FractureSet *faultSingle);
	void setGeometryDipAngle(int geometry_data_type,double geometry_value, std::vector< std::string >& name,const std::vector<double> &min_max_mode, FractureSet *faultSingle);
	void setGeometryHeight(int geometry_data_type,double geometry_value, std::vector< std::string >& name,const std::vector<double> &min_max_mode, FractureSet *faultSingle);
	void setGeometryLength(int geometry_data_type,double geometry_value, std::vector< std::string >& name,const std::vector<double> &min_max_mode, FractureSet *faultSingle);
	void setGeometryAperture(int geometry_data_type,double geometry_value, std::vector< std::string >& name,const std::vector<double> &min_max_mode, FractureSet *faultSingle);

	void setGeometryAzimuthArray(int geometry_data_type,double geometry_value, std::vector< std::string >& name,const std::vector<double> &min_max_mode, FaultArraySet *faultSingle);
	void setGeometryAngleArray(int geometry_data_type,double geometry_value, std::vector< std::string >& name,const std::vector<double> &min_max_mode, FaultArraySet *faultSingle);
	void setGeometryHeightArray(int geometry_data_type,double geometry_value, std::vector< std::string >& name,const std::vector<double> &min_max_mode, FaultArraySet *faultSingle);
	void setGeometryLengthArray(int geometry_data_type,double geometry_value, std::vector< std::string >& name,const std::vector<double> &min_max_mode, FaultArraySet *faultSingle);
	void setGeometryApertureArray(int geometry_data_type,double geometry_value, std::vector< std::string >& name,const std::vector<double> &min_max_mode, FaultArraySet *faultSingle);
	void setGeometryDamageHangingArray(int geometry_data_type,double geometry_value, std::vector< std::string >& name,const std::vector<double> &min_max_mode, FaultArraySet *faultSingle);
	void setGeometryDamageFootwallArray(int geometry_data_type,double geometry_value, std::vector< std::string >& name,const std::vector<double> &min_max_mode, FaultArraySet *faultSingle);

	void setGeometrySingleFaultCorridor(int geometry_parameter_type,int geometry_data_type, double geometry_value, std::vector<std::string>& name,const std::vector<double> &min_max_mode, FractureSet *faultSingle);
	void setGeometrySingleFaultZone(int geometry_parameter_type,int geometry_data_type, double geometry_value, std::vector<std::string>& name,const std::vector<double> &min_max_mode, FaultArraySet *faultSingle);


	void stage_hypothesis_fractures_set_generation(const geode::StructuralModel &model,bool estim_volume);
	void setGeometrySinglefault(int geometry_parameter_type,
			int geometry_data_type, double geometry_value,
			std::vector< std::string > &name, const std::vector<double> &min_max_mode,
			FractureSet *faultSingle);

	void setGeometryArray(int geometry_parameter_type,
			int geometry_data_type, double geometry_value,
			std::vector< std::string > &name, const std::vector<double> &min_max_mode,
			FaultArraySet *faultArray);
	void setGeometryCorridor(int geometry_parameter_type,
		int geometry_data_type, double geometry_value, std::vector<std::string> &name,
		const std::vector<double> &min_max_mode, FaultArraySet *faultArray);
	void setGeometryFaultZone(int geometry_parameter_type,
		int geometry_data_type, double geometry_value, std::vector<std::string> &name,
		const std::vector<double> &min_max_mode, FaultArraySet *faultArray);
	void setDistributionSingleFault(int distribution_parameter_type,
			int distribution_density_law_type, double distribution_value,
			std::vector< std::string > &name, std::vector< double >& min_mode_max, FractureSet *faultSingle);
	void setDistributionArray(int distribution_parameter_type,
			int distribution_density_law_type, double distribution_value,
			std::vector< std::string > &name, std::vector< double >& min_mode_max, FaultArraySet *faultArray);
	void setDistributionArrayFault(int distribution_density_law_type,int distribution_data_type,
			double distribution_value, std::vector< std::string > &name, std::vector< double >& min_mode_max,
			FaultArraySet *faultArray);
	void setDistributionArrayStressZoneFault(int distribution_density_law_type,
			double distribution_value, std::vector< std::string > &name, std::vector< double >& min_mode_max,
			FaultArraySet *faultArray);
	void setDistributionArrayFaultFamili1(int distribution_density_law_type,int distribution_data_type,
			double distribution_value, std::vector< std::string > &name, std::vector< double >& min_mode_max,
			FaultArraySet *faultArray);

	void setDistributionHangingWall(int distribution_density_law_type,int distribution_data_type,
			double distribution_value, std::vector< std::string > &name, std::vector< double >& min_mode_max,
			FaultArraySet *faultArray);

	void setDistributionFootWall(int distribution_density_law_type,int distribution_data_type,
			double distribution_value, std::vector< std::string > &name, std::vector< double >& min_mode_max,
			FaultArraySet *faultArray);

	void setDistributionArray(int distribution_density_law_type,int distribution_array_type,
			int distribution_parameter_type, int distribution_data_type,
			double distribution_value, std::vector< std::string > &name, std::vector< double >& min_mode_max,
			FaultArraySet *faultArray);
	std::vector<std::string> splitLine(std::istringstream& my_stream, char delim);
	bool saveFromStage(std::string surface_filename2, std::string path);
	void saveFractureSet(fractures_intersect::model& model,const std::string &surface_filename2, FractureSet *setFrac);
	void saveArrayFault(std::string methode_suffixe,const std::string &surface_filename2, FaultArraySet* setArrayZone);
	//std::map< std::pair<int,int>, std::shared_ptr<fractures_intersect::fracture_set_geo> >
	void create_bed_surface(const geode::StructuralModel &model,std::vector< FracResBeddingAperture* > bed_interface_list, std::string beddingPropertyName, int index);
	void create_bed_parallel(const geode::StructuralModel &model,boost::uint32_t seed,FracResBeddingProperties* bed_interface);
	void update_bed_interface(const geode::StructuralModel &model,fractures_intersect::model& model_fractures_loc,boost::uint32_t seed,FracResBeddingProperties* bed_interface);

	void initialize_aperture_slope_distribution( boost::mt19937 &center_gen, boost::function< double()>& rand , std::vector< double > minmaxmode_values, std::string type_name);
	void setDistributionArrayEchelonSpaceFault(int distribution_data_type, double distribution_value,std::vector< std::string >& name,std::vector< double >& min_max_mode,FaultArraySet *faultArray) ;
	void initialize_bed_aperture_distribution( boost::mt19937 &center_gen, boost::function< double()>& rand , std::vector< double > minmaxmode_values, std::string type_name);

	geode::BoundingBox3D computePolyhedronBoundingBox(const geode::SolidMesh3D &mesh, geode::index_t polyhedronIndex);
	geode::index_t closest_element_index(geode::Point3D pt,const geode::AABBTree3D& tree,const BoxAABBEvalDistance3D& disteval);
	std::vector<geode::BoundingBox3D> create_bounding_box(const geode::StructuralModel &model);

	int runsimulation(std::string fileParameter, std::string fileStrm);
	void fractures_set_generation_for_scenario(const geode::StructuralModel &model,std::string nameFile, bool estim_volume) ;
	void create_result_directory(const std::string &fileParameter, size_t posFirst,const std::string &namefile, std::string &path);
	void run_fracturing_simulation(const geode::StructuralModel &model,const std::string &fileParameter);
	void create_main_directory_output(std::string &path);
	void create_local_directory_output(std::string &path);
	void simulationAllScenario(const geode::StructuralModel &model,std::string nameFile);
	std::string getPathAndFileName(const std::string &filename,std::string &namefile);
	size_t getParameterSuffixeName(const std::string &fileParameter,std::string &namefile,std::string& path );
	bool saveMechanicalUnit(std::string surface_filename2,std::string path);

	void create_fault_zone_dfn(const FracResStage &stage,
			const FracResHypothesis &fracResHypothesis, int ireal,
			const geode::StructuralModel &model, bool estim_volume,
			const std::string &realisation_directory,
			FaultArraySet *setArrayZone,
			fractures_intersect::model &model_fractures_loc);
	void create_corridors_dfn(const FracResStage &stage,
			const FracResHypothesis &fracResHypothesis, int ireal,
			const std::string &realisation_directory,
			const geode::StructuralModel &model, bool estim_volume,
			FaultArraySet *setArrayCorridor,
			fractures_intersect::model &model_fractures_loc);
	int create_bed_parallel_slip_dfn( const geode::StructuralModel &model,fractures_intersect::model& model_fractures_loc,
			FracResBeddingProperties *bed_parallel,std::string path,std::string output_prefixe_file_name);
	void create_fault_array_dfn(const FracResStage &stage,
			const FracResHypothesis &fracResHypothesis, int ireal,
			const std::string &realisation_directory,
			const geode::StructuralModel &model, bool estim_volume,
			FaultArraySet *setArrayZone,
			fractures_intersect::model &model_fractures_loc);
	void create_single_fault_dfn(const FracResStage &stage,
			const FracResHypothesis &fracResHypothesis, int ireal,
			fractures_intersect::model &model_fractures_loc,
			const std::string &realisation_directory,
			const geode::StructuralModel &model, bool estim_volume,
			FractureSet *setFrac);
	void time_creation_message(const std::string &nameDfn,std::chrono::time_point < std::chrono::system_clock>& start_point);
	void create_bed_interface_slip_dfn(const geode::StructuralModel &model,
			const FracResHypothesis &fracResHypothesis,
			fractures_intersect::model &model_fractures_loc,std::string path,std::string output_prefixe_file_name, int index);
	void message_bedding_interface(const FracResStage &stage);
	int message_hypothesis_dfn_status(FracResHypothesis& hypothesis);
	void create_bed_interface_from_hypothesis(
			const FracResHypothesis &fracResHypothesis,
			int indice_bed_interface, const geode::StructuralModel &model,
			fractures_intersect::model &model_fractures_loc,
			std::string nameDfn,std::string path,std::string output_prefixe_file_name);
	int create_single_faults_from_hypothesis(
			const FracResHypothesis &fracResHypothesis, int indice_single_fault,
			const FracResStage &stage, int ireal,
			fractures_intersect::model &model_fractures_loc,
			const std::string &realisation_directory,
			const geode::StructuralModel &model, bool estim_volume,
			std::string nameDfn, std::string output_prefixe_name);
	int create_fault_array_from_hypothesis(
			const FracResHypothesis &fracResHypothesis, int indice_fault_array,
			const FracResStage &stage, int ireal,
			const std::string &realisation_directory,
			const geode::StructuralModel &model, bool estim_volume,
			std::string nameDfn, std::string output_prefixe_name,
			fractures_intersect::model &model_fractures_loc);
	int create_bed_parallel_slip_from_hypothesis(
			const FracResHypothesis &fracResHypothesis, int indice_bed_parallel,
			const geode::StructuralModel &model,
			fractures_intersect::model &model_fractures_loc,
			std::string nameDfn,std::string path,std::string output_prefixe_file_name);

	int create_clusters_corridors_from_hypothesis(
			const FracResHypothesis &fracResHypothesis,
			int indice_fault_corridor, const FracResStage &stage, int ireal,
			const std::string &realisation_directory,
			const geode::StructuralModel &model, bool estim_volume,
			std::string nameDfn,std::string output_prefixe_name,fractures_intersect::model &model_fractures_loc);

	int create_clusters_fault_zone_from_hypothesis(
			const FracResHypothesis &fracResHypothesis, int indice_fault_zone,
			const FracResStage &stage, int ireal,
			const geode::StructuralModel &model, bool estim_volume,
			const std::string &realisation_directory, std::string nameDfn,
			std::string output_prefixe_name,fractures_intersect::model &model_fractures_loc);

	void run_simulation_from_stage_and_hypothesesis(
			const FracResHypothesis &fracResHypothesis,
			const geode::StructuralModel &model,
			const std::string &nameDfn,
			const FracResStage &stage, int ireal,
			const std::string &realisation_directory, bool estim_volume,
			fractures_intersect::model &model_fractures_loc,
			FracResScenario &scenario);


	void save_fracture_surface_mesh_from_model(int ifaille,fractures_intersect::model& model,std::string output_file_native);

	void simulationAllStage(const geode::StructuralModel &model);

	void getSuccessMessage();
	std::string getFormatTime(int time_second);
	void getEvaluateFracturesNumberMessage();

	int getEvaluateFracturesNumberPerStageHypothesis(const geode::StructuralModel &model,const FracResStage &stage,
			const FracResHypothesis &fracResHypothesis);
	void getEvaluateFracturesNumberPerAllStage( const geode::StructuralModel &model);


	void evaluate_single_fault_dfn(const FracResStage &stage,
			const FracResHypothesis &fracResHypothesis,
			const geode::StructuralModel &model, bool estim_volume,
			FractureSet *setFrac);
	int get_evaluate_fracture_single_faults_from_hypothesis(const FracResHypothesis &fracResHypothesis, int indice_single_fault,
			const FracResStage &stage, const geode::StructuralModel &model) ;


	int evaluate_fault_array_dfn(const FracResStage &stage,
			const FracResHypothesis &fracResHypothesis,
			const geode::StructuralModel &model, bool estim_volume,
			FaultArraySet *setFrac);
	int get_evaluate_fracture_fault_array_from_hypothesis(const FracResHypothesis &fracResHypothesis, int indice_fault_array,
			const FracResStage &stage, const geode::StructuralModel &model) ;

	int get_evaluate_fracture_cluster_fault_zone_from_hypothesis(const FracResHypothesis &fracResHypothesis, int indice_fault_array,
			const FracResStage &stage, const geode::StructuralModel &model) ;
	int get_evaluate_fracture_cluster_corridor_from_hypothesis(const FracResHypothesis &fracResHypothesis, int indice_fault_array,
			const FracResStage &stage, const geode::StructuralModel &model) ;

	bool save_output_single_faults_from_hypothesis(fractures_intersect::model &model_fractures_loc,const FracResHypothesis &fracResHypothesis, int indice_single_fault);
	bool save_output_array_faults_from_hypothesis(std::string methode_suffixe,const FracResHypothesis &fracResHypothesis, int indice_fault_array, int type);
	void save_output_from_stage_and_hypothesesis(const FracResHypothesis &fracResHypothesis,const FracResStage &stage, fractures_intersect::model &model_fractures_loc,FracResScenario &scenario);

	void save_bed_parallel_slip_from_hypothesis(const FracResHypothesis &fracResHypothesis, int indice_bed_parallel,fractures_intersect::model &model_fractures_loc);
	void save_bed_interface_from_hypothesis(const FracResHypothesis &fracResHypothesis, int indice_bed_interface,fractures_intersect::model &model_fractures_loc) ;

	std::map<std::pair<int,int>,std::string> fractures_number_evaluation_for_scenario(const geode::StructuralModel &model);

	void affiche_fracture_number_evaluation(std::map<std::pair<int,int>,std::string>&  mapNumber);
	void get_property_attribute_names_from_block_model(const geode::StructuralModel &model);
/**
 *
 * @fn int getCategoryType(std::string)
 * @brief
 *
 * @param value
 * @return
 */

	int getCategoryType(std::string value){

		int cat_index = -1;
		if(value== "BEDDING" || value== "BED_INTERFACE"){
			cat_index = 0;
		}else if(value== "JOINTS"){
			cat_index = 3;
		}else if(value== "FAULTS"){
			cat_index = 4;
		}else if(value== "FAULTS_ARRAYS"){
			cat_index = 5;
		}else if(value== "FRACTURE_CLUSTERS"){
			cat_index = 7;
		}else if(value== "BED_PARALLEL"){
			cat_index = 2;
		}

		return cat_index;

	};
/**
 * @fn int getFractureType(std::string)
 * @brief
 *
 * @param value
 * @return
 */
	int getFractureType( std::string value){

		int cat_index = -1;

		if(value == "BED_PARALLEL_SLIP"){
			cat_index = 0;
		}else if(value == "SYSTEMATIC_JOINTS"){
			cat_index = 1;
		}else if(value == "CONJUGATE_JOINTS"){
			cat_index = 2;
		}else if(value == "SINGLE_FAULTS"){
			cat_index = 3;
		}else if(value == "CONJUGATE_FAULTS"){
			cat_index = 4;
		}else if(value == "ANASTOMOSING_FAULTS"){
			cat_index = 5;
		}else if(value == "RELAY_FAULTS"){
			cat_index = 6;
		}else if(value == "EN_ECHELON_FAULTS"){
			cat_index = 7;
		}else if(value == "JOINTS_SWARM"){
			cat_index = 8;
		}else if(value == "FRACTURE_CORRIDORS"){
			cat_index = 9;
		}else if(value == "FAULT_ZONES"){
			cat_index = 10;
		}

		return cat_index;
	};
/**
 * @fn int getIntersectionCategoryType(std::string)
 * @brief
 *
 * @param value
 * @return
 */
	int getIntersectionCategoryType( std::string value){
	  int cat_index=-1;
	if(value=="BEDDINGS" || value=="BED_INTERFACE"){
		cat_index = 0;
	}else if(value=="MECHANICAL_UNITS"){
		cat_index = 1;
	}else if(value=="BED_PARALLEL"){
	   cat_index = 2;
	}else if(value=="JOINTS"){
		cat_index = 3;
	}else if(value=="FAULTS"){
		cat_index = 4;
	}else if(value=="FAULTS_ARRAYS"){
		cat_index = 5;
	}else if(value=="FRACTURE_CLUSTERS_JOINT_SWARM"){
		cat_index = 6;
	}else if(value=="FRACTURE_CLUSTERS_CORRIDORS_ZONES"){
		cat_index = 7;
	}

	return cat_index;
	};
/**
 * @fn fractures_intersect::FracResIntersectionTypeEnum getIntersectionType(std::string)
 * @brief
 *
 * @param value
 * @return
 */
	fractures_intersect::FracResIntersectionTypeEnum getIntersectionType(std::string value){
		fractures_intersect::FracResIntersectionTypeEnum  tab_type_intersection;
	   if(value=="ABUTTING_CROSSING"){
		   tab_type_intersection =   fractures_intersect::FracResIntersectionTypeEnum::ABBUTING_CROSSING;
	   }else if(value=="ABUTTING"){
		   tab_type_intersection =   fractures_intersect::FracResIntersectionTypeEnum::ABUTTING;
	   }else if(value=="CROSSING"){
		   tab_type_intersection =   fractures_intersect::FracResIntersectionTypeEnum::CROSSING;
	   }else if(value=="RANDOM"){
		   tab_type_intersection =   fractures_intersect::FracResIntersectionTypeEnum::RANDOM;
	   }
	   return tab_type_intersection;
	};
/**
 * @fn int getDistributionType(std::string)
 * @brief
 *
 * @param value
 * @return
 */
	int getDistributionType( std::string value){
	  int cat_index=-1;
	if(value=="FRACTURE_CORRIDORS_DENSITY"){
		cat_index = 0;
	}else if(value=="FAULT_ZONE_DENSITY"){
		cat_index = 1;
	}else if(value=="FAULT_SEGMENTS_DENSITY"){
	   cat_index = 2;
	}else if(value=="FAULT_DENSITY"){
		cat_index = 3;
	}else if(value=="MINIMUM_SPACE_BETWEEN_ARRAYS"){
		cat_index = 4;
	}else if(value=="ARRAY_DENSITY"){
		cat_index = 5;
	}else if(value=="FAULT_SEGMENT_DENSITY"){
		cat_index = 6;
	}else if(value=="DISTANCE_BETWEEN_SEG"){
		cat_index = 7;
	}else if(value=="WIDTH_STRESS_SHADOW_FAMILY1"){
		cat_index = 8;
	}else if(value=="WIDTH_STRESS_SHADOW_FAMILY2"){
		cat_index = 9;
	}else if(value=="FAULT_SEGMENTS_DENSITY_CORRIDORS"){
		cat_index = 10;
	}else if(value=="WIDTH_STRESS_SHADOW_ZONE"){
		cat_index = 11;
	}else if(value=="MINIMUM_SPACING_BETWEEN_FAULT_ZONES"){
		cat_index = 12;
	}else if(value=="MINIMUM_SPACING_BETWEEN_FRACTURE_CORRIDORS"){
		cat_index = 13;
	}else if(value=="CLUSTER_DENSITY"){
		cat_index = 14;
	}else if(value=="CLUSTER_NUMBER"){
		cat_index = 15;
	}else if (value == "FAULT_SEGMENTS_DENSITY_CORRIDORS"){
		cat_index = 16;
	}

	return cat_index;
	};
/**
 * @fn int getDistributionArrayType(std::string)
 * @brief
 *
 * @param value
 * @return
 */
	int getDistributionArrayType( std::string value){
	  int cat_index=-1;
	if(value=="SINGLE_FAULT"){
		cat_index = 0;
	}else if(value=="FAULT_ARRAY_ARRAY"){
		cat_index = 1;
	}else if(value=="FAULT_ARRAY_SEGMENT"){
	   cat_index = 2;
	}else if(value=="FRACTURE_CLUSTERS_ZONE_CLUSTERS"){
		cat_index = 3;
	}else if(value=="FRACTURE_CLUSTERS_ZONE_CLUSTERS_STRUCTURAL"){
		cat_index = 4;
	}else if(value=="FRACTURE_CLUSTERS_ZONE_FAULT_SEGMENTS_HANGING"){
		cat_index = 5;
	}else if(value=="FRACTURE_CLUSTERS_ZONE_FAULT_SEGMENTS_FOOTWALL"){
		cat_index = 6;
	}else if(value=="FRACTURE_CLUSTERS_ZONE_FAULT_SEGMENTS"){
		cat_index = 7;
	}else if(value=="FRACTURE_CLUSTERS_CORRIDOR_FAULT_SEGMENTS"){
		cat_index = 8;
	}else if(value=="FRACTURE_CLUSTERS_CORRIDOR_FRACTURE_CORRIDORS"){
		cat_index = 9;
	}

	return cat_index;
	};

	int getGeometryType( std::string value){
	  int cat_index=-1;
	if(value=="ANGLE_WITH_FAUL_ARRAY"){
		cat_index = 0;
	}else if(value=="AZIMUTH"){
		cat_index = 1;
	}else if(value=="DIP_AZIMUTH"){
	   cat_index = 2;
	}else if(value=="DIP_ANGLE"){
		cat_index = 3;
	}else if(value=="LENGTH"){
		cat_index = 4;
	}else if(value=="HEIGHT"){
		cat_index = 5;
	}else if(value=="APERTURE"){
		cat_index = 6;
	}else if(value=="WIDTH"){
		cat_index = 7;
	}else if(value=="CORE_ZONE_DIP_AZIMUTH"){
		cat_index = 8;
	}else if(value=="CORE_ZONE_DIP_ANGLE"){
		cat_index = 9;
	}else if(value=="CORE_ZONE_LENGTH"){
		cat_index = 10;
	}else if(value=="CORE_ZONE_HEIGHT"){
		cat_index = 11;
	}else if(value=="CORE_ZONE_APERTURE"){
		cat_index = 12;
	}else if(value=="DAMAGE_ZONE_WIDTH"){
		cat_index = 13;
	}else if(value=="HANGING_WALL_WIDTH"){
		cat_index = 14;
	}else if(value=="FOOTWALL_WIDTH"){
		cat_index = 15;
	}else if(value=="ANGLE_ARRAY_DIRECTION"){
		cat_index = 16;
	}else if(value=="MAX_DISTANCE_CONNECTIVITY"){
		cat_index = 17;
	}

	return cat_index;
	};




	int getGeometryFaultZoneType( std::string value){
	  int cat_index=-1;
	if(value=="DIP_SLIP_NORMAL"){
		cat_index = 0;
	}else if(value=="DIP_SLIP_REVERSE"){
		cat_index = 1;
	}else if(value=="STRIKE_SLIP"){
	   cat_index = 2;
	}

	return cat_index;
	};




	int getGeometryArrayType( std::string value){
	  int cat_index=-1;
	if(value=="SINGLE_FAULT"){
		cat_index = 0;
	}else if(value=="FAULT_ARRAY_ARRAYS"){
		cat_index = 1;
	}else if(value=="FAULT_SEGMENTS_COMMON"){
	   cat_index = 2;
	}else if(value=="FAULT_SEGMENTS_FAMILY_1"){
		cat_index = 3;
	}else if(value=="FAULt_SEGMENTS_FAMILY_2"){
		cat_index = 4;
	}else if(value=="FRACTURE_CLUSTERS_ZONE_CLUSTERS"){
		cat_index = 5;
	}else if(value=="FRACTURE_CLUSTERS_ZONE_FAULT_SEGMENTS"){
		cat_index = 6;
	}else if(value=="FRACTURE_CLUSTERS_CORRIDOR_CLUSTERS"){
		cat_index = 7;
	}else if(value=="FRACTURE_CLUSTERS_CORRIDOR_FAULT_SEGMENTS"){
		cat_index = 8;
	}

	return cat_index;
	};

	int getDataType( std::string value){
	  int cat_index=-1;
	if(value=="CONSTANT"){
		cat_index = 0;
	}else if(value=="DISTRIBUTION"){
		cat_index = 1;
	}else if(value=="PROPERTY"){
	   cat_index = 2;
	}else if(value=="NORMAL_DISTRIBUTION"){
		cat_index = 3;
	}
	return cat_index;
	};

	int getClustersFaultType( std::string value){
	  int cat_index=-1;
	if(value=="Null"){
		cat_index = 0;
	}else if(value=="DIP_SLIP_NORMAL"){
		cat_index = 1;
	}else if(value=="DIP_SLIP_REVERSE"){
	   cat_index = 2;
	}else if(value=="STRIKE_SLIP"){
		cat_index = 3;
	}
	return cat_index;
	};

	int getGeometryCorrelationType( std::string value){
		  int cat_index=-1;
		if(value=="INDEPENDANT"){
			cat_index = 0;
		}else if(value=="CORRELATED"){
			cat_index = 1;
		}
		return cat_index;
	};

	int getDensityLawType( std::string value){
	  int cat_index=-1;
	  int pos1=0;
	  int pos2 = value.find("P00", pos1);
	  if ( pos2 != value.npos )
	  {
		  cat_index = 0;

	  }else{

	  pos2 = value.find("P30", pos1);
	  if ( pos2 != value.npos )
	  {
		  cat_index = 1;

	  }else{
		   pos2 = value.find("P32", pos1);
		  if ( pos2 != value.npos )
		  {
			  cat_index = 2;

		  }
		  pos2 = value.find("Number of", pos1);
		  if ( pos2 != value.npos )
		  {
			  cat_index = 0;

		  }
	  }
	  }
	   return cat_index;
	};
    /*!
       * Number of fine cells.
       */
      geode::index_t cellCount() const {
          return cells_.size();
	}

	void getTimeMessage();


};

// inline methods


 }
#endif  // RemFracResSim_h
