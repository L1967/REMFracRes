//----------------------------  graph_sim.h  ---------------------------
// 
//
//  Description:	graph_sim class use for build graph from geode StructuralModel
//					information
// 
//  Documents:		fractureset.h
//
//  Creation date  : 10/07/2022
//
//  Author:			Rabah Namar
//
//  Copyright (C) 2022 by in2earth
//
//---------------------------- fractureset.h  ---------------------------
#pragma once

#ifndef fractureset_h
#define fractureset_h

// local includes
#include <remfracres/REMFracRes.h>
#include <remfracres/CellId.h>
#include <remfracres/BoxAABBEvalDistance.h>
#include <remfracres/DistributionPower.h>
#include <numeric/model.h>
#include <primitive_mesh/fault_array.h>
// standard includes
#include <fstream>
#include <iostream>
#include <iomanip>
#include <memory>
#include <vector>
#include <string>

// boost includes 
#include <boost/random.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <Eigen/Core>
#include <Eigen/Geometry>

// geode include
#include <geode/basic/uuid.h>
#include <geode/geometry/aabb.h>
#include <geode/geometry/distance.h>
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
 
 class FractureSet
{

// variables

public:


	const double ndvalue_ = -9999.0;



	/**
	 *  array type of fracture set if single unknown
	*/

	enum array_type { unkown, bedding, anastomosing, echelon, relay, corridor, fault_zone, clusters};

	/**
	 *  Distribution type
	*/

	enum distribution_type { uniform, triangulated, normal, lognormal, power};
	/**
	 *      name of FractureSet  
	 */

	std::string name_;
	/**
	 *      index of FractureSet
	 */
	int fracture_set_index_;

	/**
	 *      Distribution name used for  
	 *		the determination of length_1
	 */

	std::string length_1_;

	/**
	 *      Distribution name used for  
	 *		the determination of length_2
	 */

	std::string length_2_;

	/**
	 *      Distribution name used for  
	 *		the determination of dip azimut
	 */
	
	std::string dip_azimut_;

	/**
	 *      Distribution name used for  
	 *		the determination of dip angle
	 */

	std::string dip_angle_;

	/**
	 *      Distribution name used for  
	 *		the determination of the fracture aperture
	 */
	
	std::string aperture_name_;

	/**
	 *      name of property used for  
	 *		the determination of the fracture density
	 */
	
	std::string density_property_name_;


	std::string density_distribution_name_;

	/**
	 *      name of property used for
	 *		the determination of the hanging wall fracture density
	 */

	std::string density_hanging_wall_property_name_;


	std::string density_hanging_wall_distribution_name_;


	/**
	 *      name of property used for
	 *		the determination of the foot wall fracture density
	 */

	std::string density_foot_wall_property_name_;


	std::string density_foot_wall_distribution_name_;


	std::string shadow_zone_distribution_name_;

	std::string shadow_zone_property_name_;




	/**
	 *
	 * 		density min et max
	 *
	*/

	std::vector< double > density_min_max_;

	/**
	 *      name of property used for  
	 *		constraining of the fracture 
	 */
	
	std::string constraint_property_name_;

	/**
	 *      fracture density , mean density for property grid
	 */
	
	double density_;

	double density_hanging_wall_;

	double density_foot_wall_;

	/**
	 *      true if the complementary fracture set is active 
	 */
	
	bool complementary_fracture_set_active_;


	/**
	 *      pointer function  for the deterministic function       
	 *      (density variation with depth ) used in the generation of fracture  
	 *      centers (using a regionalised Poisson process)
	 */

	float (*density_function_) (float);

	/**
	 *      coefficient a of determinitic function used in generation of fracture centers
	 */

	double coef_a_; 

	/**
	 *      seed value for random generator  
	 */
		
	long int seed_;

	/**
	 *		number of generated fractures  
	 */

	int fractures_number_;

	/**
	 *		number of possible generated fractures  
	 */

	int evaluate_fractures_number_;


	/**
	 *		number of generated fractures in hanging  wall damage zone
	 */

	int fractures_number_hanging_wall_;

	/**
	 *		number of possible generated fractures in hanging wall damage zone
	 */

	int evaluate_fractures_number_hanging_wall_;



	/**
	 *		number of generated fractures  in foot wall damage zone
	 */

	int fractures_number_foot_wall_;

	/**
	 *		number of possible generated fractures   in hanging damage zone
	 */

	int evaluate_fractures_number_foot_wall_;


	/**
	 *		number of inserted fractures
	 */

	int fractures_inserted_;

	/**
	 *		true if the fracture is editing  
	 */

	bool is_edit_;

    /**
	 *		is imported, true if is imported from fraca
	 */

	bool is_imported_;
	/**
	 *		true if the fracture is created
	 */

	bool is_created_;

	/**
		 *		true if the fracture is created
		 */

	int distribution_length1_type_index_;

	int distribution_length2_type_index_;

	int distribution_azimut_type_index_;

	int distribution_angle_type_index_;

	int distribution_aperture_type_index_;

	int distribution_shadow_zone_type_index_;




	/**
		 *		true if the fracture is created
		 */

	int density_type_index_;
	int density_type_hanging_wall_index_;
	int density_type_foot_wall_index_;

	float distr_length1_val_;
	float distr_length2_val_;
	float distr_azimut_val_;
	float distr_angle_val_;
	float distr_aperture_val_;
	float distr_shadow_zone_val_;

	float distr_hanging_wall_val_;
	float distr_foot_wall_val_;


	bool to_generate_;



	absl::flat_hash_map< std::string,geode::uuid >  geode_BlocName_to_BlockId_;

	absl::flat_hash_map< std::string,bool >  geode_BlocName_to_Density_exist_;

	absl::flat_hash_map< std::string,double >  geode_BlocName_to_volume_;

	absl::flat_hash_map< std::string,double >  geode_BlocName_to_area_;

	absl::flat_hash_map< std::string,double >  geode_BlocName_to_density_;

	std::vector< double > Dip_;
	std::vector< double > Orientation_;
	std::vector< double > length1_;
	std::vector< double > length2_;
	std::vector< double > aperture_;
	std::vector< double > density_val_;
	std::vector< double > density_hanging_wall_val_;
	std::vector< double > density_foot_wall_val_;
	std::vector< double > shadow_zone_val_;

	std::vector< geode::Point3D> vertices_;
	std::vector< std::array< geode::index_t,3 > > triangles_;
	std::vector< double > triangles_l1_;
	std::vector< double > triangles_l2_;
	std::vector< double > triangles_orientation_;
	std::vector< double > triangles_dip_;
	std::vector< double > triangles_aperture_;


	distribution_type length1_type_;
	distribution_type length2_type_;
	distribution_type azimuth_type_;
	distribution_type dip_angle_type_;
	distribution_type aperture_type_;
	distribution_type density_type_;
	distribution_type density_hanging_wall_type_;
	distribution_type density_foot_wall_type_;
	distribution_type shadow_zone_type_;

	std::string length1_property_name_;
	std::string length2_property_name_;
	std::string azimut_property_name_;
	std::string dip_angle_property_name_;
	std::string aperture_property_name_;



	std::vector< double > minmaxmode_values_;
	std::vector< double > minmaxmode_density_values_;
	std::vector< double > minmaxmode_density_hanging_wall_values_;
	std::vector< double > minmaxmode_density_foot_wall_values_;

	std::vector< double > minmaxmode_shadow_zone_values_;

	boost::function< double()> randLength1_;
	boost::function< double()> randLength2_;
	boost::function< double()> randAzimut_;
	boost::function< double()> randAzimutAnastomosing_;
	boost::function< double()> randDip_;
	boost::function< double()> randAperture_;
	boost::function< double()> randDensity_;
	boost::function< double()> randHangingWallDensity_;
	boost::function< double()> randFootWallDensity_;
	boost::function< double()> randShadowZone_;

	int distribution_echelon_type_index_;
	distribution_type echelon_space_type_;
	std::vector< double > minmaxmode_echelon_space_values_;
	boost::function< double()> randEchelon_space_;
	std::vector< double > echelon_space_val_;
	std::string echelon_space_distribution_name_;

	std::shared_ptr< fractures_intersect::fracture_set_geo  > stress_zone_set_;

	std::vector< fractures_intersect::FracResIntersectionTypeEnum > tab_type_intersection_;

	fractures_intersect::FracResIntersectionCategoryEnum type_category_;

	fractures_intersect::FracResIntersectionTypeEnum type_intersection_;

	double default_stress_zone_max_;

	std::string directory_path_;


    /**
     *  boolean value for shadows zone active;
     *
     */

	bool is_shadow_zone_active_;


	bool success_;


	int cylinder_step_;

	int array_index_;

	array_type type_array_;

	double echelon_space_;

	double max_distance_connectivity_;

	int nb_fractures_connection_relay_;

	int width_stress_zone_law_type_;

	std::string width_stress_zone_distribution_name_;

	std::vector< double >  width_stress_zone_min_max_mode_;

	std::string width_stress_zone_property_name_;

	float mult_z;

	bool is_active_;

	int model_geo_cont_index_;

	std::string output_prefixe_name_;
	double total_simulation_time_;
	double total_writting_output_time_;

	double length_max_;



public:
	
	
	/**
	 *      number of activable fracture (fracture not yet used in 
	 *      the karstification phases)   
	 */
	
	int fractures_activable_;
	
	int fractures_activable_hanging_wall_;

	int fractures_activable_foot_wall_;

	/**
	 *      Pointer of TSurf object (graphics representation 
	 *      of the fractures) a Tsurf contains all the fractures 
	 *      repesentation, a fractures = 2 TFaces
	 */

	//TSurf* fractures_;
	
	/**
	 *      generation domaine limits 
	 *      (grid limits + fracture max length/2)
	 */

	geode::BoundingBox3D generation_box_;


	fractures_intersect::oriented_bounding_box oriented_generation_box_;

	fractures_intersect::oriented_bounding_box* oriented_generation_box_damage_zone_;






	/**
	 *      volume of the grid where the generation take place 
	 */

	float volume_;

	/**
	 *		area of the grid where the generation take place 
	 */

	float area_;

	/**
	 *      maximun fracture density in property grid
	 */
	
	double max_density_;

	/**
	 *		index number of each property for include data property values
	 */ 

     int index_leng1_, index_leng2_	, index_orient_ ,index_dip_ ,index_aper_ , index_density_, index_shadow_zone_;


     fractures_intersect::model model_fractures_;

     fractures_intersect::Cboite_rect* bbox_;


     /**
      * Density law type
      * 0 : POO   nb_fractures
      * 1 : P32   density (m3) * volume
      * 2 : P31   density (m2/m3) / surface_moyenne * volume;
      */

     int density_law_type_;

     int density_hanging_wall_law_type_;


     int density_foot_wall_law_type_;

     fractures_intersect::fracture_set_geo fracture_set_geo_;



     std::vector< double > cells_density_;

 	DistributionPower randDistPowerLength1_;
 	DistributionPower randDistPowerLength2_;
 	DistributionPower randDistPowerAzimuth_;
 	DistributionPower randDistPowerAzimuthAnastomosing_;
 	DistributionPower randDistPowerDipAngle_;
 	DistributionPower randDistPowerAperture_;
 	DistributionPower randDistPowerDensity_;
 	DistributionPower randDistPowerMinimumSpace_;
 	DistributionPower randDistPowerHangingWall_;
 	DistributionPower randDistPowerFootWall_;
 	DistributionPower randDistPowerEchelonSpace_;
 	DistributionPower randDistPowerShadowZone_;

	
// functions

public: 

	/**
	 *      Default constructor 
	 */

    FractureSet();

	/**
	 *		Copy constructor
	 * 
	 *      @param from The value to copy to this object
	 */

    FractureSet(const FractureSet& from);

	/**
	 *		Destructor 
	 */
        
	~FractureSet(); 

	/** 
     *      Assignment operator 
	 *
     *      @param from the value to assign to this object 
     * 
     *      @return A reference to this object 
     */
  
	FractureSet& operator=(const FractureSet& from); 

	/** 
     *      copy function 
	 *
     *      @param from the value to assign to this object 
     */
  
	void copy(const FractureSet& from); 
	
	/** 
     *      show_graphics: register in the geobase the TSurf of the FractureSet  
	 *		to show the distribution in the KartSim grid  
     */
	
	void show_graphics();

	/**
	 *      set_generation_box: set the domain limits
	 *                       
	 *		@param box: 3D box, delimiting the domain 
	 */
	
	void set_generation_box(geode::BoundingBox3D box);



	/**
	 *      set_generation_oriented_box: set the domain limits
	 *
	 *		@param box: 3D box, delimiting the domain
	 */

	void set_generation_oriented_box(fractures_intersect::oriented_bounding_box box);

	/**
	 *      fractures_generation: generation of the fractures 
	 */
	
	void fractures_generation(const geode::StructuralModel &model, std::shared_ptr< geode::AABBTree3D >& tree,

		    std::shared_ptr< BoxAABBEvalDistance3D >& disteval, std::vector<CellId> cells,std::shared_ptr<fractures_intersect::fault_array>& test_array,bool estim_volume = false);


	/**
	 *      add_fracture: create a TFace and add it in the Tsurf
	 *      if no intersection with conditional info return true
	 *
	 *		@param center: center of the fracture 
	 *		@param l1: lenght1 in m
	 *		@param l2: lenght1 in m 
	 *		@param azimut: dip azimut 
	 *		@param angle: dip angle
	 *		@param aperture: aperture in m
	 */

	bool add_fracture(geode::Point3D p0_xyz, geode::Point3D p0_uv , double l1 , double l2, double orientation,
    						   double dip, double aperture, bool check_intersect ,float* constraint_data, float* density_data);

	/**
	 *      get_fractures: get the TSurf Object representing 
	 *      the fractures 
	 *
	 *		@return A pointer to a TSurf
	 */

	//TSurf* get_fractures();

	/**
	 *      set_fractures: set the TSurf Object representing 
	 *      the fractures 
	 *
	 *		@param tsurf A pointer to a TSurf
	 */

	// void set_fractures(TSurf* tsurf);

	/**
	 *      get_fractures_number: get the number of fractures 
	 *
	 *		@return  The number of fractures
	 */
	
	int get_fractures_number();
	

	/**
	 *		set the number of fractures
	 * 
	 *		@param n_frac the number of fractures
	 */

	void set_fractures_number(int n_frac);


	/**
	 *      get_fractures_number: get the evaluate number of fractures 
	 *
	 *		@return  The number of fractures
	 */
	
	int get_evaluate_fractures_number();
	


	/**
	 *      get_fractures_activables: get the number of activable fractures
	 *
	 *		@return The number of activable fractures
	 */
	
	int get_fractures_activables();
	
	/**
	 *      desactivate_fractures: desactivate a given number of fractures
	 *                       
	 *		@param n_fracture: number of fracture to desactivate 
	 *
     *      @return the number of deactivated fractures 
	 */

	int  desactivate_fractures(int n_fracture);
	 
	/**
	 *      dactivate_fractures: activate a given number of fractures
	 *                       
	 *		@param n_fracture: number of fracture to activate 
	 *
     *      @return the number of activated fractures 
	 */
	
	int activate_fractures(int n_fracture);

	/**
	 *      generate_fractures_density: estimate the fractures density   
	 *      for the grid
	 */

	void generate_fractures_density (bool density_variable_in_grid,double density_val,const geode::StructuralModel &model,std::vector< remfracres::CellId >& cells,double lfrac_mean, int density_type);

	/**
	 *      generate_fractures_density: estimate the fractures density
	 *      for the grid
	 */

	void generate_fractures_density_in_array (bool density_variable_in_grid,double density_val,const geode::StructuralModel &model,std::vector< remfracres::CellId >& cells,double lfrac_mean, int density_type);

	/**
	 *      fill the fracture set with data of gocad fracture network    
	 *      for the karstsim sgrid
	 *
	 *		@param filename: the file name 
     * 
	 */

     bool import_gocad_fracture_network(std::string filename);

	/**
	 *      init tsurf associated to the fractureset     
	 *
	 *		@param name: fracture set name  
     * 
	 */

     void init_tsurf(std::string name);

    /**
	 *      add_fracture: create a TFace and add it in the Tsurf with vector
	 *
	 *		@param center: center of the fracture 
	 *		@param v_major: major vector
	 *		@param v_minor: minor vector 
	 *		@param aperture: aperture in m
	 */

     void add_fracture_with_vector(geode::Point3D center, geode::Vector3D v_major ,geode:: Vector3D v_minor,float aperture);

    /**
	 *      Initialise the fractureset property to no datat 
     *      value 
     *
     *		@param property_name: name of the property
	 */

    void initialise_property_to_nd_value(std::string property_name);

    /**
	 *      change de phase property for the next fractures  
     *
     *		@param frac_number
     *      @param phase_number
	 */

    void change_phase_property(int frac_number,int phase_number);

	/**
	 *      save the fracture set data
	 *      for the karstsim 
	 *
	 *		@param out: ostreamb 
     * 
	 */

     bool save(std::string filename);

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
 	 *      open a fracture set data
 	 *      from
 	 *
 	 *		@param fin: fstream
 	 *
 	 *		@param catalog: pointer to an AsciiCatalog
 	 */

       bool open_loc(std::fstream& fin);

       /**
               * @fn void initialize_density_from_structural_model(const geode::StructuralModel&)
        * @brief
        *
        * @param model
        */
      void initialize_density_P32_from_structural_model(const geode::StructuralModel& model,std::vector<CellId>& cells,double lfrac_mean);


      /**
              * @fn void initialize_density_from_structural_model(const geode::StructuralModel&)
       * @brief
       *
       * @param model
       */
     void initialize_density_P31_from_structural_model_in_array(const geode::StructuralModel& model,std::vector<CellId>& cells);

     /**
             * @fn void initialize_density_from_structural_model(const geode::StructuralModel&)
      * @brief
      *
      * @param model
      */
    void initialize_density_P32_from_structural_model_in_array(const geode::StructuralModel& model,std::vector<CellId>& cells,double lfrac_mean);


    /**
            * @fn void initialize_density_from_structural_model(const geode::StructuralModel&)
     * @brief
     *
     * @param model
     */
   void initialize_density_P31_from_structural_model(const geode::StructuralModel& model,std::vector<CellId>& cells);
      /**
       * @fn void initialize_volume_from_global_grid_of_structural_model(const geode::StructuralModel&)
       * @brief
       *
       * @param model
       */
      void initialize_volume_from_global_grid_of_structural_model(const geode::StructuralModel& model);
      /**
       * @fn float evaluate_fracture_number(const geode::StructuralModel&, bool)
       * @brief
       *
       * @param model
       * @param density_variable_in_grid
       * @return
       */

      float evaluate_fracture_number(const geode::StructuralModel& model, std::vector<CellId>& cells,bool density_variable_in_grid);
      /**
       * @fn void evaluate_final_fracture_number(boost::uniform_01<boost::mt19937>&)
       * @brief
       *
       * @param center_uni
       */

       void evaluate_final_fracture_number( boost::uniform_01<boost::mt19937>& center_uni);
       /**
        * @fn void initialize_dip(std::vector<double>&)
        * @brief
        *
        * @param new_dip
        */

      void initialize_dip(std::vector<double>& new_dip);
      /**
       * @fn void initialize_orientation(std::vector<double>&)
       * @brief
       *
       * @param new_orientation
       */
      void initialize_orientation(std::vector<double>& new_orientation);
      /**
       * @fn void initialize_length1(std::vector<double>&)
       * @brief
       *
       * @param new_length1
       */
      void initialize_length1(std::vector<double>& new_length1);
      /**
       * @fn void initialize_length2(std::vector<double>&)
       * @brief
       *
       * @param new_length2
       */
      void initialize_length2(std::vector<double>& new_length2);
      /**
       * @fn void initialize_aperture(std::vector<double>&)
       * @brief
       *
       * @param new_aperture
       */
      void initialize_aperture(std::vector<double>& new_aperture);
      /**
       * @fn void init_surface_geometry(int, const geode::TriangulatedSurface3D&, geode::TriangulatedSurfaceBuilder3D&)
       * @brief
       *
       * @param indice
       * @param edged_curve
       * @param builder
       */
      void init_surface_geometry( int indice,const geode::TriangulatedSurface3D& edged_curve,geode::TriangulatedSurfaceBuilder3D& builder );
		/**
		 * @fn void set_distribution_type_from_index(int, std::string)
		 * @brief
		 *
		 * @param index
		 * @param name
		 */
      void set_distribution_type_from_index(int index,std::string name);
      /**
       * @fn distribution_type get_distribution_type(std::string)
       * @brief
       *
       * @param name
       * @return
       */
      distribution_type get_distribution_type(std::string name);
      /**
       * @fn void initialize_length1_distribution(boost::mt19937&)
       * @brief
       *
       * @param center_gen
       */
	  void initialize_length1_distribution(boost::mt19937 &center_gen);
	  /**
	   * @fn void initialize_lenght2_distribution(boost::mt19937&)
	   * @brief
	   *
	   * @param center_gen
	   */
	  void initialize_lenght2_distribution(boost::mt19937 &center_gen);
	  /**
	   * @fn void initialize_azimut_distribution(boost::mt19937&)
	   * @brief
	   *
	   * @param center_gen
	   */
	  void initialize_azimut_distribution(boost::mt19937 &center_gen);
	  /**
	   * @fn void initialize_dip_angle_distribution(boost::mt19937&)
	   * @brief
	   *
	   * @param center_gen
	   */
	  void initialize_dip_angle_distribution(boost::mt19937 &center_gen);
	  /**
	   * @fn void initialize_aperture_distribution(boost::mt19937&)
	   * @brief
	   *
	   * @param center_gen
	   */
	  void initialize_aperture_distribution(boost::mt19937 &center_gen);
	  /**
	   * @fn void initialize_fracture_set_distribution(boost::mt19937&)
	   * @brief
	   *
	   * @param center_gen
	   */
	  void initialize_fracture_set_distribution(boost::mt19937& center_gen);

	  /**
	   * @fn void save_trinagulate_surface_to_ts(std::string, std::string)
	   * @brief
	   *
	   * @param filename
	   * @param name_surface
	   */
	  void save_trinagulate_surface_to_ts(std::string filename, std::string name_surface);

	  /**
	   * @fn bool test_All_fractures_boolean()
	   * @brief
	   *
	   * @return
	   */
	  bool test_All_fractures_boolean();

	  /**
	   * @fn void clean_vertices_and_triangles_array()
	   * @brief
	   *
	   */
	  void clean_vertices_and_triangles_array();
	  /**
	   * @fn void initialize_geometrical_fractures_data()
	   * @brief
	   *
	   */

	  void initialize_geometrical_fractures_data();
	  /**
	   * @fn double update_data_from_distribution(boost::mt19937&)
	   * @brief
	   *
	   * @param center_gen
	   * @return
	   */

	  double update_data_from_distribution(boost::mt19937 &center_gen);

	  /**
	  * @fn bool save_fracture_set(std::string)
	   * @brief
	   *
	   * @param surface_filename
	   * @return
	   */
	  bool save_fracture_set(std::string surface_filename);
	/**
	 * @fn void set_volume(float)
	 * @brief
	 *
	 * @param val
	 */
	  void set_volume(float val);
	/**
	 * @fn void set_area(float)
	 * @brief
	 *
	 * @param val
	 */
	  void set_area(float val);


	  /**
	   * @fn void update_bounding_box_for_fault_zone(double)
	   * @brief
	   *
	   * @param damage_zone_width
	   */


	  void update_bounding_box_for_fault_zone(double damage_zone_width);
	  /**
	  	   * @fn void fractures_generation_from_shadow_zone(const geode::StructuralModel&, bool)
	   * @brief
	   *
	   * @param model
	   * @param estim_volume
	   */
	  void fractures_generation_from_shadow_zone(const geode::StructuralModel &model, std::shared_ptr< geode::AABBTree3D >& tree,

			    std::shared_ptr< BoxAABBEvalDistance3D >& disteval, std::vector<CellId> cells, std::shared_ptr<fractures_intersect::fault_array>& test_array,  bool estim_volume);


	  std::shared_ptr< fractures_intersect::Cfaille> getRelayConnectionFaille( int indexFractureSet, std::shared_ptr<fractures_intersect::Cfaille> faille1, std::shared_ptr<fractures_intersect::Cfaille> facille2);
	  /**
	  	   * @fn void fractures_generation_from_shadow_zone(const geode::StructuralModel&, bool)
	   * @brief
	   *
	   * @param model
	   * @param estim_volume
	   */
	  void fractures_generation_from_shadow_zone_in_array(const geode::StructuralModel &model, std::shared_ptr< geode::AABBTree3D >& tree,

			    std::shared_ptr< BoxAABBEvalDistance3D >& disteval, std::vector<CellId> cells,std::shared_ptr<fractures_intersect::fault_array>& test_array,  bool estim_volume, fractures_intersect::Point3D pt_ref, int indice_array_loc);

	  void fractures_generation_from_shadow_zone_in_corridor(const geode::StructuralModel &model, std::shared_ptr< geode::AABBTree3D >& tree,

			    std::shared_ptr< BoxAABBEvalDistance3D >& disteval, std::vector<CellId> cells,std::shared_ptr<fractures_intersect::fault_array>& test_array,  bool estim_volume, fractures_intersect::Point3D pt_ref);

	  float evaluate_fracture_number_from_oriented_box(const geode::StructuralModel &model, std::vector<CellId>& cells,bool density_variable_in_grid);

	  void fractures_generation_from_shadow_zone_in_damage_fault_zone(const geode::StructuralModel &model, std::shared_ptr< geode::AABBTree3D >& tree,

			    std::shared_ptr< BoxAABBEvalDistance3D >& disteval, std::vector<CellId> cells, std::shared_ptr<fractures_intersect::fault_array>& test_array, bool estim_volume, double damage_hanging_wall_zone, double damage_foot_wall_zone);


	  std::vector< fractures_intersect::Point3D  > get_echelon_position_from_initial_germ(fractures_intersect::Point3D, int ifrac);


	  void initialize_density_distribution(boost::mt19937 &center_gen) ;

	  void initialize_shadow_zone_distribution(boost::mt19937 &center_gen) ;


	  float evaluate_fracture_number_from_oriented_box_hanging_wall_zone(const geode::StructuralModel &model,std::vector<CellId>& cells,bool density_variable_in_grid,double ratio_hanging);

	  float evaluate_fracture_number_from_oriented_box_foot_wall_zone(const geode::StructuralModel &model,std::vector<CellId>& cells,bool density_variable_in_grid,double ratio_footwall);

	  void set_fractures_number_hanging_wall(int n_frac);
	  void set_fractures_number_foot_wall(int n_frac);

	  void evaluate_final_fracture_number_hanging_wall(boost::uniform_01<boost::mt19937>& center_uni);
	  void evaluate_final_fracture_number_foot_wall(boost::uniform_01<boost::mt19937>& center_uni);
	  void initialize_foot_wall_density_distribution(boost::mt19937 &center_gen) ;
	  void initialize_hanging_wall_density_distribution(boost::mt19937 &center_gen) ;

	  void initialize_echelon_space_density_distribution(boost::mt19937 &center_gen) ;
	  bool update_values_from_fracture_index(bool complementary, int &ifrac,
			double &length_max_l1, double &length_max_l2,boost::mt19937& center_gen);

	  void update_values_of_fracture_index_from_property( remfracres::CellId id_cell, const geode::StructuralModel &model,int ifrac);

	  int calculation_array_fractures_from_property(const geode::StructuralModel &model, std::shared_ptr< geode::AABBTree3D >& tree,

	 			    std::shared_ptr< BoxAABBEvalDistance3D >& disteval, std::vector<CellId> cells,std::shared_ptr<fractures_intersect::fault_array>& test_array,  bool estim_volume,boost::uniform_01<boost::mt19937>& center_uni);


	  int calculation_fractures_from_property(const geode::StructuralModel &model, std::shared_ptr< geode::AABBTree3D >& tree,

			    std::shared_ptr< BoxAABBEvalDistance3D >& disteval, std::vector<CellId> cells, std::shared_ptr<fractures_intersect::fault_array>& test_array, bool estim_volume,boost::uniform_01<boost::mt19937>& center_uni);

	  int calculation_hanging_wall_fractures_from_property(const geode::StructuralModel &model, std::shared_ptr< geode::AABBTree3D >& tree,

			    std::shared_ptr< BoxAABBEvalDistance3D >& disteval, std::vector<CellId> cells, std::shared_ptr<fractures_intersect::fault_array>& test_array, bool estim_volume,boost::uniform_01<boost::mt19937>& center_uni, double damage_foot_wall_zone,double damage_hanging_wall_zone);

	  int calculation_foot_wall_fractures_from_property(const geode::StructuralModel &model, std::shared_ptr< geode::AABBTree3D >& tree,

			    std::shared_ptr< BoxAABBEvalDistance3D >& disteval, std::vector<CellId> cells,std::shared_ptr<fractures_intersect::fault_array>& test_array,  bool estim_volume,boost::uniform_01<boost::mt19937>& center_uni,double damage_foot_wall_zone);

	  void fractures_evaluate_generation(const geode::StructuralModel &model, std::shared_ptr< geode::AABBTree3D >& tree,std::shared_ptr< BoxAABBEvalDistance3D >& disteval, std::vector<CellId> cells,std::shared_ptr<fractures_intersect::fault_array>& test_array,  bool estim_volume);
};
// inline methods


 }
#endif  // fractureset_h
