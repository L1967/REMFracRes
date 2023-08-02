//----------------------------  FaultArraySet.h  ---------------------------
// 
//
//  Description:	FaultArraySet class use for build array fault
//					information
// 
//  Documents:		FaultArraySet.h
//
//  Creation date  : 01/11/2022
//
//  Author:			Pascal Siegel - Rabah Namar
//
//  Copyright (C) 2022 by in2earth
//
//---------------------------- FaultArraySet.h  ---------------------------
#pragma once

#ifndef faultarrayset_h
#define faultarrayset_h

// local includes
#include <remfracres/REMFracRes.h>
#include <remfracres/FracResClustersFaultTypeEnum.h>
#include <remfracres/FracResGeometryCorrelationTypeEnum.h>
#include <remfracres/FracResDataTypeEnum.h>
#include <remfracres/BoxAABBEvalDistance.h>
#include <remfracres/DistributionPower.h>
#include <numeric/model.h>
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

#include <geode/geosciences/explicit/common.h>
#include <geode/geosciences_io/mesh/common.h>
#include <geode/geosciences_io/model/common.h>



#include <geode/geosciences/explicit/representation/builder/structural_model_builder.h>
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

#include <primitive_mesh/fault_arrays_set_geo.h>
#include <primitive_mesh/fault_array.h>
#include <remfracres/fractureset.h>

namespace remfracres{

/**
 *    FaultArraySet contains the parameters for the simulation
 *    of a Boolean model and the generated Fault Array
 * 
 *    @author Pascal Siegel - Rabah Namar - Olivier Jaquet 
 *    @date 2022
 */
 
 class FaultArraySet
{

// variables

public:

	/**
	 *  array type of fracture set if single unknown
	*/

	enum cluster_type {no_cluster, dip_slip_normal, dip_slip_reverse, strike_slip};

	/**
	 *  array type of fracture set if single unknown
	*/

	enum array_type { unkown, bedding, anastomosing, echelon, relay, corridor, fault_zone, cluster};


	const double ndvalue_ = -9999.0;

	/**
	 *  Distribution type
	*/

	enum distribution_type { uniform, triangulated, normal, lognormal,power};
	/**
	 *      name of FractureSet  
	 */

	std::string name_;
	/**
	 *      index of FractureSet
	 */
	int fault_array_set_index_;

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

	/**
	 *      name of property used for
	 *		the determination of the fracture density
	 */


	std::string damage_zone_hanging_property_name_;

	std::string damage_zone_footwall_property_name_;


	/**
	 *      name of property used for
	 *		the determination of the fracture density
	 */


	std::string damage_zone_hanging_distribution_name_;
	std::string damage_zone_footwall_distribution_name_;


	/**
	 *      name of density used for
	 *		the determination of the array density
	 */

	std::string density_distribution_name_;

	/**
	 *
	 * 		density min et max
	 *
	*/

	std::vector< double > density_min_max_;
	/**
	 *
	 * 		minimum space between array distribution name
	 *
	*/

	std::string minimum_space_distribution_name_;

	/**
	 *
	 * 		minimum space between array property name
	 *
	*/

	std::string minimum_space_property_name_;



	/**
	 *      name of property used for  
	 *		constraining of the fracture 
	 */
	
	std::string constraint_property_name_;

	/**
	 *      fracture density , mean density for property grid
	 */
	
	double density_;

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

	int faultArray_number_;

	/**
	 *		number of possible generated fractures  
	 */

	int evaluate_faultArray_number_;

	/**
	 *		number of inserted fractures
	 */

	int faultArray_inserted_;

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

	int distribution_minimum_space_type_index_;

	int distribution_damage_zone_hanging_index_;
	int distribution_damage_zone_footwall_index_;

	/**
		 *		true if the fracture is created
		 */

	int density_type_index_;

	float distr_length1_val_;
	float distr_length2_val_;
	float distr_azimut_val_;

	float distr_angle_val_;

	float distr_aperture_val_;
	float distr_minimum_space_val_;
	float distr_damage_zone_hanging_val_;
	float distr_damage_zone_footwall_val_;

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
	std::vector< double > minimum_space_val_;
	std::vector< double > damage_zone_hanging_val_;
	std::vector< double > damage_zone_footwall_val_;

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
	distribution_type minimum_space_type_;
	distribution_type damage_zone_footwall_type_;
	distribution_type damage_zone_hanging_type_;
	std::string length1_property_name_;
	std::string length2_property_name_;
	std::string azimut_property_name_;
	std::string dip_angle_property_name_;
	std::string aperture_property_name_;
	std::string damage_hanging_wall_property_name_;
	std::string damage_foot_wall_property_name_;


	std::vector< double > minmaxmode_values_;
	std::vector< double > minmaxmode_density_values_;
	std::vector< double > minmaxmode_minimum_space_values_;
	std::vector< double > minmaxmode_damage_zone_hanging_values_;
	std::vector< double > minmaxmode_damage_zone_footwall_values_;

	boost::function< double()> randLength1_;
	boost::function< double()> randLength2_;
	boost::function< double()> randAzimut_;
	boost::function< double()> randDip_;
	boost::function< double()> randAperture_;
	boost::function< double()> randDensity_;
	boost::function< double()> randMinimumSpace_;
	boost::function< double()> randDamageHangingZone_;
	boost::function< double()> randDamageFootwallZone_;

	std::shared_ptr< fractures_intersect::fault_arrays_set_geo  > fault_array_geo_set_;

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

	/**
	 *      volume of the grid where the generation take place
	 */

	float volume_;

	/**
	 *		area of the grid where the generation take place
	 */

	float area_;

 	/**
 	 *		list of  fractures set element
 	 */
     std::vector< FractureSet > fracturesSet_list_;

 	/**
 	 *      fault array set number
 	 */

 	int nb_fracture_sets_;

	/**
 	 *      fault array set number
 	 */

 	int nb_fault_array_sets_;


 	FracResClustersFaultTypeEnum fracResClusterType_;

 	bool usePreviousFracture_;

 	array_type typeFamilyArray_;

 	cluster_type typeCluster_;

 	FracResGeometryCorrelationTypeEnum heightCorrelationType_;
 	double heightCorrelationValue_;
 	bool heightCorrelationVariable_;
 	std::string heightCorrelationProperty_name_;

 	FracResGeometryCorrelationTypeEnum apertureCorrelationType_;
 	double apertureCorrelationValue_;
 	bool apertureCorrelationVariable_;
 	std::string apertureCorrelationProperty_name_;

 	FracResDataTypeEnum faultCoreType_;
 	double faultCoreValue_;
 	std::string faultCoreDistribution_name_;
 	std::string faultCoreProperty_name_;

	double rate_reactive_;

	bool is_active_;


	int model_geo_cont_index_;

	std::string output_prefixe_name_;
	double total_simulation_time_;
	double total_writting_output_time_;

	double length_max_;
	DistributionPower randDistPowerLength1_;
	DistributionPower randDistPowerLength2_;
	DistributionPower randDistPowerAzimuth_;
	DistributionPower randDistPowerDipAngle_;
	DistributionPower randDistPowerAperture_;
	DistributionPower randDistPowerDensity_;
	DistributionPower randDistPowerMinimumSpace_;
	DistributionPower randDistPowerHangingWall_;
	DistributionPower randDistPowerFootWall_;

public:
	
	
	/**
	 *      number of activable fracture (fracture not yet used in 
	 *      the karstification phases)   
	 */
	
	int faultArray_activable_;
	
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



	/**
	 *      maximun fracture density in property grid
	 */
	
	float max_density_;

	/**
	 *		index number of each property for include data property values
	 */ 

     int index_leng1_, index_leng2_	, index_orient_ ,index_dip_ ,index_aper_ ;


     fractures_intersect::model model_fractures_;

     fractures_intersect::Cboite_rect* bbox_;


     /**
      * Density law type
      * 0 : POO   nb_fractures
      * 1 : P32   density (m3) * volume
      * 2 : P31   density (m2/m3) / surface_moyenne * volume;
      */

     int density_law_type_;

     /**
	  * minimum_space_between_array_law_type_
	  * 0 : constant
	  * 1 : distribution name
	  * 2 : property name;
	  */

	 int minimum_space_between_array_law_type_;

	 double minimum_space_between_array_law_value_;

	 std::string minimum_space_between_array_distribution_name_;

	 std::vector< double >  minimum_space_between_array_min_max_mode_;

	 std::string minimum_space_between_array_property_name_;


     int nb_failles_;

     std::vector< double > cells_density_;
// functions

public: 

	/**
	 *      Default constructor 
	 */

     FaultArraySet();

	/**
	 *		Copy constructor
	 * 
	 *      @param from The value to copy to this object
	 */

     FaultArraySet(const FaultArraySet& from);

	/**
	 *		Destructor 
	 */
        
	~FaultArraySet();

	/** 
     *      Assignment operator 
	 *
     *      @param from the value to assign to this object 
     * 
     *      @return A reference to this object 
     */
  
	FaultArraySet& operator=(const FaultArraySet& from);

	/** 
     *      copy function 
	 *
     *      @param from the value to assign to this object 
     */
  
	void copy(const FaultArraySet& from);
	
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
	 *      fractures_generation: generation of the fractures 
	 */
	
	void fault_array_generation(const geode::StructuralModel &model, std::shared_ptr< geode::AABBTree3D >& tree,

    std::shared_ptr< BoxAABBEvalDistance3D >& disteval,  std::vector<CellId> cells,fractures_intersect::Ccontrainte_Geo& geo_contraint,bool estim_volume = false);


	/**
	 *      get_fractures_number: get the number of fractures 
	 *
	 *		@return  The number of fractures
	 */
	
	int get_fault_array_number();
	

	/**
	 *		set the number of fractures
	 * 
	 *		@param n_frac the number of fractures
	 */

	void set_fault_array_number(int n_frac);


	/**
	 *      get_fractures_number: get the evaluate number of fractures 
	 *
	 *		@return  The number of fractures
	 */
	
	int get_evaluate_fault_array_number();
	


	/**
	 *      get_fractures_activables: get the number of activable fractures
	 *
	 *		@return The number of activable fractures
	 */
	
	int get_fault_array_activables();
	
	/**
	 *      desactivate_fractures: desactivate a given number of fractures
	 *                       
	 *		@param n_fracture: number of fracture to desactivate 
	 *
     *      @return the number of deactivated fractures 
	 */

	int  desactivate_fault_array(int n_fracture);
	 
	/**
	 *      dactivate_fractures: activate a given number of fractures
	 *                       
	 *		@param n_fracture: number of fracture to activate 
	 *
     *      @return the number of activated fractures 
	 */
	
	int activate_fault_array(int n_fracture);


    


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
    void initialize_density_P31_from_structural_model(const geode::StructuralModel& model,std::vector<CellId>& cells);


      /**
       * @fn void initialize_volume_from_global_grid_of_structural_model(const geode::StructuralModel&)
       * @brief
       *
       * @param model
       */
      void initialize_volume_from_global_grid_of_structural_model(const geode::StructuralModel& model);
      /**
       * @fn float evaluate_fault_array_number(const geode::StructuralModel&, bool)
       * @brief
       *
       * @param model
       * @param density_variable_in_grid
       * @return
       */

      float evaluate_fault_array_number(const geode::StructuralModel& model, std::vector<CellId> cells,bool density_variable_in_grid);
      /**
       * @fn void evaluate_final_fracture_number(boost::uniform_01<boost::mt19937>&)
       * @brief
       *
       * @param center_uni
       */

       void evaluate_final_fault_array_number( boost::uniform_01<boost::mt19937>& center_uni);
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
	   * @fn void initialize_aperture_distribution(boost::mt19937&)
	   * @brief
	   *
	   * @param center_gen
	   */
	  void initialize_density_distribution(boost::mt19937 &center_gen);

	  /**
	   * @fn void initialize_fracture_set_distribution(boost::mt19937&)
	   * @brief
	   *
	   * @param center_gen
	   */
	  void initialize_fault_array_set_distribution(boost::mt19937& center_gen);

	  /**
	   * @fn double update_data_from_distribution(boost::mt19937&)
	   * @brief
	   *
	   * @param center_gen
	   * @return
	   */

	  double update_data_from_distribution(boost::mt19937 &center_gen);

	  bool add_fault_array(geode::Point3D p0_xyz, geode::Point3D p0_uv , double l1 , double l2, double orientation,
	      						   double dip, double aperture, bool check_intersect ,float* constraint_data);

	/**
	*      generate_fault array_density: estimate the fault array density
	*      for the grid
	*/

	void generate_fault_array_density (bool density_variable_in_grid,double density_val,const geode::StructuralModel &model,std::vector< remfracres::CellId >& cells,double lfrac_mean, int density_type);
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


	void initialize_geometrical_fractures_data() ;

	void init_surface_geometry( int indice,const geode::TriangulatedSurface3D& surface , geode::TriangulatedSurfaceBuilder3D& builder );

	 void initialize_minimum_space_distribution(boost::mt19937 &center_gen) ;

	 void fault_corridor_generation(const geode::StructuralModel &model, std::shared_ptr< geode::AABBTree3D >& tree,

			    std::shared_ptr< BoxAABBEvalDistance3D >& disteval, std::vector<CellId> cells, fractures_intersect::Ccontrainte_Geo& geo_contraint,bool estim_volume );

	 void fault_zone_generation(const geode::StructuralModel &model, std::shared_ptr< geode::AABBTree3D >& tree,

			    std::shared_ptr< BoxAABBEvalDistance3D >& disteval, std::vector<CellId> cells,fractures_intersect::Ccontrainte_Geo& geo_contraint, bool estim_volume );

	 void initialize_damage_zone_hanging_distribution(boost::mt19937& center_gen);
	 void initialize_damage_zone_footwall_distribution(boost::mt19937& center_gen);

	 void update_values_of_array_index_from_property( remfracres::CellId id_cell, const geode::StructuralModel &model,int ifrac);

	 int calculation_array_from_property(const geode::StructuralModel &model, std::shared_ptr< geode::AABBTree3D >& tree,

			    std::shared_ptr< BoxAABBEvalDistance3D >& disteval, std::vector<CellId> cells,fractures_intersect::Ccontrainte_Geo &geo_contraint,bool estim_volume, boost::uniform_01<boost::mt19937>& center_uni);

	 int  calculation_array_fault_zone_from_property(const geode::StructuralModel &model, std::shared_ptr< geode::AABBTree3D >& tree,

	  			    std::shared_ptr< BoxAABBEvalDistance3D >& disteval, std::vector<CellId> cells,fractures_intersect::Ccontrainte_Geo& geo_contraint,bool estim_volume, boost::uniform_01<boost::mt19937>& center_uni);

	 int  calculation_array_corridor_from_property(const geode::StructuralModel &model, std::shared_ptr< geode::AABBTree3D >& tree,

	  			    std::shared_ptr< BoxAABBEvalDistance3D >& disteval, std::vector<CellId> cells,fractures_intersect::Ccontrainte_Geo& geo_contraint,bool estim_volume, boost::uniform_01<boost::mt19937>& center_uni);

	  int fault_array_fractures_evaluation(const geode::StructuralModel &model, std::shared_ptr< geode::AABBTree3D >& tree, std::shared_ptr< BoxAABBEvalDistance3D >& disteval, std::vector<CellId> cells,bool estim_volume );

};
// inline methods


 }
#endif  // fractureset_h
