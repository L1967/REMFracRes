//----------------------------  RemFracResFileManagement.h  ---------------------------
// 
//
//  Description:	RemFracResFileManagement class use for File FracRes management read
//					information
// 
//  Documents:		RemFracResFileManagement.h
//
//  Creation date  : 02/03/2023
//
//  Author:			Rabah Namar
//
//  Copyright (C) 2022 by in2earth
//
//---------------------------- RemFracResFileManagement.h  ---------------------------
#pragma once

#ifndef RemFracResFileManagement_h
#define RemFracResFileManagement_h

// local includes
#include <remfracres/FracResStage.h>
#include <remfracres/FracResUnit.h>

//#include <remfracres/FaultArraySet.h>
#include <remfracres/fracResFileKeyWord.h>
#include <primitive_mesh/tetra.h>
#include <numeric/model.h>
// standard includes
#include <iomanip>
#include <vector>
#include <list>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <utility>

// boost includes 

// geode include




namespace remfracres{

/**
 *    FractureSet contains the parameters for the simulation 
 *    of a Boolean model and the generated fractures
 * 
 *    @author Pascal Siegel - Rabah Namar - Olivier Jaquet 
 *    @date 2023
 */
 
 class RemFracResFileManagement
{

// variables

public:

	/**
	 *      name of RemFracresFileManagement
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


	FracResUnit* fracResUnit_;


 	/**
  	 *      general seed value for intersection and fracture set
  	 */

  	int seed_;
 	/**
  	 *		boolean for outpout dircetory same as input file
  	 */
      bool use_global_outpout_directory_path_;

      /**
      */
      std::string outpout_directory_file_;

      FracResFileKeyWord*  fracresFileKeyWord_;


// functions

public: 

	/**
	 *      Default constructor 
	 */

     RemFracResFileManagement();

	/**
	 *		Copy constructor
	 * 
	 *      @param from The value to copy to this object
	 */

     RemFracResFileManagement(const RemFracResFileManagement& from);

	/**
	 *		Destructor 
	 */
        
	~RemFracResFileManagement();

	/** 
     *      Assignment operator 
	 *
     *      @param from the value to assign to this object 
     * 
     *      @return A reference to this object 
     */
  
	RemFracResFileManagement& operator=(const RemFracResFileManagement& from);

	/** 
     *      copy function 
	 *
     *      @param from the value to assign to this object 
     */
  
	void copy(const RemFracResFileManagement& from);
	

      /**
       * @fn std::string enumCategoryToString(fractures_intersect::FracResIntersectionCategoryEnum)
       * @brief
       *
       * @param cat
       * @return
       */
      std::string enumCategoryToString(fractures_intersect::FracResIntersectionCategoryEnum cat);

      /**
       * @fn std::string enumIntersectionToString(fractures_intersect::FracResIntersectionTypeEnum)
       * @brief
       *
       * @param cat
       * @return
       */

      std::string enumIntersectionToString(fractures_intersect::FracResIntersectionTypeEnum cat);


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

	/**
	* @fn bool openFracres(std::string)
	 * @brief
	 *
	 * @param filename
	 * @return
	 */

	bool openFracres(std::string filename);

	/**
	* @fn void manageScenario(std::fstream&)
	 * @brief
	 *
	 * @param fin
	 */
	void manageScenario(std::fstream& fin);
	/**
		 * @fn void manageRunParameter(std::fstream&)
	 * @brief
	 *
	 * @param fin
	 */
	void manageRunParameter(std::fstream& fin);

	/**
		 * @fn void manageMatrixProp(std::fstream&)
	 * @brief
	 *
	 * @param fin
	 */
	void manageMatrixProp(std::fstream& fin);
	/**
		 * @fn void manageOutputs(std::fstream&)
	 * @brief
	 *
	 * @param fin
	 */
	void manageOutputs(std::fstream& fin);
/**
 * @fn void manageStageHypScenarioKeyVal(std::fstream&)
 * @brief
 *
 * @param fin
 */
    void manageStageHypScenarioKeyVal(std::fstream& fin);
    /**
         * @fn void manageUpscalingRep(std::fstream&)
     * @brief
     *
     * @param fin
     */
	void manageUpscalingRep(std::fstream& fin);
	/**
		 * @fn void manageRepresentationObject(std::fstream&)
	 * @brief
	 *
	 * @param fin
	 */
	void manageRepresentationObject(std::fstream& fin);
	/**
		 * @fn FracResFracture manageFracture(std::fstream&)
	 * @brief
	 *
	 * @param fin
	 * @return
	 */
	FracResFracture manageFracture(std::fstream& fin);
	/**
		 * @fn FracResHypothesis manageHypothesis(std::fstream&)
	 * @brief
	 *
	 * @param fin
	 * @return
	 */

	FracResHypothesis manageHypothesis(std::fstream& fin);
	/**
		 * @fn void manageBeddingStageHypothesis(std::fstream&, FracResHypothesis&)
	 * @brief
	 *
	 * @param fin
	 * @param fracResHypothesis
	 */

	void manageBeddingStageHypothesis(std::fstream& fin, FracResHypothesis& fracResHypothesis );
	/**
		 * @fn FracResBeddingAperture manageBeddingAperture*(std::fstream&)
	 * @brief
	 *
	 * @param fin
	 * @return
	 */
	FracResBeddingAperture* manageBeddingAperture(std::fstream& fin);
	/**
		 * @fn void manageClassicStageHypothesis(std::fstream&, FracResStage&, FracResHypothesis&)
	 * @brief
	 *
	 * @param fin
	 * @param fracresStage
	 * @param fracResHypothesis
	 */
	void manageClassicStageHypothesis(std::fstream& fin, FracResStage& fracresStage, FracResHypothesis& fracResHypothesis );
	/**
		 * @fn void manageFractureProperty(std::fstream&, FracResStage&, FracResHypothesis&)
	 * @brief
	 *
	 * @param fin
	 * @param fracresStage
	 * @param fracResHypothesis
	 */
	void manageFractureProperty(std::fstream& fin, FracResStage& fracresStage, FracResHypothesis& fracResHypothesis);
	/**
		 * @fn fractures_intersect::FracResIntersectionTypeEnum manageIntersection(std::fstream&, int&)
	 * @brief
	 *
	 * @param fin
	 * @param index
	 * @return
	 */
	fractures_intersect::FracResIntersectionTypeEnum  manageIntersection(std::fstream& fin, int& index);
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
	 * @param min_max_mode
	 * @return
	 */
	std::vector< std::string > manageDistribution(std::fstream& fin,int& distribution_parameter_type,int& distribution_array_type, int& distribution_data_type, int& distribution_density_law_type,double& distribution_value,std::vector<double> &min_max_mode);
	/**
		 * @fn std::vector<std::string> manageGeometry(std::fstream&, int&, int&, int&, double&, std::vector<double>&)
	 * @brief
	 *
	 * @param fin
	 * @param distribution_parameter_type
	 * @param distribution_array_type
	 * @param distribution_data_type
	 * @param distribution_value
	 * @param min_max_mode
	 * @return
	 */
	std::vector< std::string >  manageGeometry(std::fstream& fin, int& distribution_parameter_type,int& distribution_array_type, int& distribution_data_type,double& distribution_value,std::vector< double>& min_max_mode);
	/**
		 * @fn void manageBedParallel(std::fstream&, FracResStage&, FracResHypothesis&)
	 * @brief
	 *
	 * @param fin
	 * @param fracresStage
	 * @param fracResHypothesis
	 */
	void manageBedParallel(std::fstream& fin, FracResStage& fracresStage, FracResHypothesis& fracResHypothesis);
	/**
		 * @fn void split(const string&, std::vector<std::pair<int,int>>&)
	 * @brief
	 *
	 * @param chaine
	 * @param elements
	 */
	void split(const string& chaine,std::vector< std::pair<int,int> >& elements);
	/**
		 * @fn void manageTotalDistribution(std::fstream&, FaultArraySet*, FractureSet*)
	 * @brief
	 *
	 * @param fin
	 * @param faultArray
	 * @param faultSingle
	 */
	void manageTotalDistribution(std::fstream& fin, FaultArraySet* faultArray,FractureSet* faultSingle);
	/**
		 * @fn void manageTotalGeometry(std::fstream&, FaultArraySet*, FractureSet*)
	 * @brief
	 *
	 * @param fin
	 * @param faultArray
	 * @param faultSingle
	 */
	void manageTotalGeometry(std::fstream& fin, FaultArraySet* faultArray,FractureSet* faultSingle);
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
	void setGeometryDipAzimuth(int geometry_data_type,double geometry_value, std::vector< std::string >& name,const std::vector<double> &min_max_mode, FractureSet *faultSingle);
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
	void setGeometryDipAngle(int geometry_data_type,double geometry_value, std::vector< std::string >& name,const std::vector<double> &min_max_mode, FractureSet *faultSingle);
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
	void setGeometryHeight(int geometry_data_type,double geometry_value, std::vector< std::string >& name,const std::vector<double> &min_max_mode, FractureSet *faultSingle);
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
	void setGeometryLength(int geometry_data_type,double geometry_value, std::vector< std::string >& name,const std::vector<double> &min_max_mode, FractureSet *faultSingle);
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
	void setGeometryAperture(int geometry_data_type,double geometry_value, std::vector< std::string >& name,const std::vector<double> &min_max_mode, FractureSet *faultSingle);
/**
 * @fn void setGeometryAzimuthArray(int, double, std::vector<std::string>&, const std::vector<double>&, FaultArraySet*)
 * @brief
 *
 * @param geometry_data_type
 * @param geometry_value
 * @param name
 * @param min_max_mode
 * @param faultSingle
 */
	void setGeometryAzimuthArray(int geometry_data_type,double geometry_value, std::vector< std::string >& name,const std::vector<double> &min_max_mode, FaultArraySet *faultSingle);
	/**
		 * @fn void setGeometryAngleArray(int, double, std::vector<std::string>&, const std::vector<double>&, FaultArraySet*)
	 * @brief
	 *
	 * @param geometry_data_type
	 * @param geometry_value
	 * @param name
	 * @param min_max_mode
	 * @param faultSingle
	 */
	void setGeometryAngleArray(int geometry_data_type,double geometry_value, std::vector< std::string >& name,const std::vector<double> &min_max_mode, FaultArraySet *faultSingle);
	/**
		 * @fn void setGeometryHeightArray(int, double, std::vector<std::string>&, const std::vector<double>&, FaultArraySet*)
	 * @brief
	 *
	 * @param geometry_data_type
	 * @param geometry_value
	 * @param name
	 * @param min_max_mode
	 * @param faultSingle
	 */
	void setGeometryHeightArray(int geometry_data_type,double geometry_value, std::vector< std::string >& name,const std::vector<double> &min_max_mode, FaultArraySet *faultSingle);

	 /**
	 	  * @fn void setGeometryLengthArray(int, double, std::vector<std::string>&, const std::vector<double>&, FaultArraySet*)
	  * @brief
	  *
	  * @param geometry_data_type
	  * @param geometry_value
	  * @param name
	  * @param min_max_mode
	  * @param faultSingle
	  */
	void setGeometryLengthArray(int geometry_data_type,double geometry_value, std::vector< std::string >& name,const std::vector<double> &min_max_mode, FaultArraySet *faultSingle);
	/**
		 * @fn void setGeometryApertureArray(int, double, std::vector<std::string>&, const std::vector<double>&, FaultArraySet*)
	 * @brief
	 *
	 * @param geometry_data_type
	 * @param geometry_value
	 * @param name
	 * @param min_max_mode
	 * @param faultSingle
	 */
	void setGeometryApertureArray(int geometry_data_type,double geometry_value, std::vector< std::string >& name,const std::vector<double> &min_max_mode, FaultArraySet *faultSingle);
	/**
		 * @fn void setGeometryDamageHangingArray(int, double, std::vector<std::string>&, const std::vector<double>&, FaultArraySet*)
	 * @brief
	 *
	 * @param geometry_data_type
	 * @param geometry_value
	 * @param name
	 * @param min_max_mode
	 * @param faultSingle
	 */
	void setGeometryDamageHangingArray(int geometry_data_type,double geometry_value, std::vector< std::string >& name,const std::vector<double> &min_max_mode, FaultArraySet *faultSingle);
	/**
		 * @fn void setGeometryDamageFootwallArray(int, double, std::vector<std::string>&, const std::vector<double>&, FaultArraySet*)
	 * @brief
	 *
	 * @param geometry_data_type
	 * @param geometry_value
	 * @param name
	 * @param min_max_mode
	 * @param faultSingle
	 */
	void setGeometryDamageFootwallArray(int geometry_data_type,double geometry_value, std::vector< std::string >& name,const std::vector<double> &min_max_mode, FaultArraySet *faultSingle);
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
	void setGeometrySingleFaultCorridor(int geometry_parameter_type,int geometry_data_type, double geometry_value, std::vector<std::string>& name,const std::vector<double> &min_max_mode, FractureSet *faultSingle);
	/**
		 * @fn void setGeometrySingleFaultZone(int, int, double, std::vector<std::string>&, const std::vector<double>&, FaultArraySet*)
	 * @brief
	 *
	 * @param geometry_parameter_type
	 * @param geometry_data_type
	 * @param geometry_value
	 * @param name
	 * @param min_max_mode
	 * @param faultSingle
	 */

	void setGeometrySingleFaultZone(int geometry_parameter_type,int geometry_data_type, double geometry_value, std::vector<std::string>& name,const std::vector<double> &min_max_mode, FaultArraySet *faultSingle);
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
	void setGeometrySinglefault(int geometry_parameter_type,
			int geometry_data_type, double geometry_value,
			std::vector< std::string > &name, const std::vector<double> &min_max_mode,
			FractureSet *faultSingle);
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
	void setGeometryArray(int geometry_parameter_type,
			int geometry_data_type, double geometry_value,
			std::vector< std::string > &name, const std::vector<double> &min_max_mode,
			FaultArraySet *faultArray);
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
	void setGeometryCorridor(int geometry_parameter_type,
		int geometry_data_type, double geometry_value, std::vector<std::string> &name,
		const std::vector<double> &min_max_mode, FaultArraySet *faultArray);
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
	void setGeometryFaultZone(int geometry_parameter_type,
		int geometry_data_type, double geometry_value, std::vector<std::string> &name,
		const std::vector<double> &min_max_mode, FaultArraySet *faultArray);
/**
 * @fn void setDistributionSingleFault(int, int, double, std::vector<std::string>&, std::vector<double>&, FractureSet*)
 * @brief
 *
 * @param distribution_parameter_type
 * @param distribution_density_law_type
 * @param distribution_value
 * @param name
 * @param min_mode_max
 * @param faultSingle
 */
	void setDistributionSingleFault(int distribution_parameter_type,
			int distribution_density_law_type, double distribution_value,
			std::vector< std::string > &name, std::vector< double >& min_mode_max, FractureSet *faultSingle);
/**
 * @fn void setDistributionArray(int, int, double, std::vector<std::string>&, std::vector<double>&, FaultArraySet*)
 * @brief
 *
 * @param distribution_parameter_type
 * @param distribution_density_law_type
 * @param distribution_value
 * @param name
 * @param min_mode_max
 * @param faultArray
 */
	void setDistributionArray(int distribution_parameter_type,
			int distribution_density_law_type, double distribution_value,
			std::vector< std::string > &name, std::vector< double >& min_mode_max, FaultArraySet *faultArray);
	/**
	* @fn void setDistributionArrayFault(int, int, double, std::vector<std::string>&, std::vector<double>&, FaultArraySet*)
	 * @brief
	 *
	 * @param distribution_density_law_type
	 * @param distribution_data_type
	 * @param distribution_value
	 * @param name
	 * @param min_mode_max
	 * @param faultArray
	 */
	void setDistributionArrayFault(int distribution_density_law_type,int distribution_data_type,
			double distribution_value, std::vector< std::string > &name, std::vector< double >& min_mode_max,
			FaultArraySet *faultArray);
	/**
	* @fn void setDistributionArrayStressZoneFault(int, double, std::vector<std::string>&, std::vector<double>&, FaultArraySet*)
	 * @brief
	 *
	 * @param distribution_density_law_type
	 * @param distribution_value
	 * @param name
	 * @param min_mode_max
	 * @param faultArray
	 */
	void setDistributionArrayStressZoneFault(int distribution_density_law_type,
			double distribution_value, std::vector< std::string > &name, std::vector< double >& min_mode_max,
			FaultArraySet *faultArray);
	/**
	* @fn void setDistributionArrayFaultFamili1(int, int, double, std::vector<std::string>&, std::vector<double>&, FaultArraySet*)
	 * @brief
	 *
	 * @param distribution_density_law_type
	 * @param distribution_data_type
	 * @param distribution_value
	 * @param name
	 * @param min_mode_max
	 * @param faultArray
	 */
	void setDistributionArrayFaultFamili1(int distribution_density_law_type,int distribution_data_type,
			double distribution_value, std::vector< std::string > &name, std::vector< double >& min_mode_max,
			FaultArraySet *faultArray);
/**
 * @fn void setDistributionHangingWall(int, int, double, std::vector<std::string>&, std::vector<double>&, FaultArraySet*)
 * @brief
 *
 * @param distribution_density_law_type
 * @param distribution_data_type
 * @param distribution_value
 * @param name
 * @param min_mode_max
 * @param faultArray
 */
	void setDistributionHangingWall(int distribution_density_law_type,int distribution_data_type,
			double distribution_value, std::vector< std::string > &name, std::vector< double >& min_mode_max,
			FaultArraySet *faultArray);
/**
 * @fn void setDistributionFootWall(int, int, double, std::vector<std::string>&, std::vector<double>&, FaultArraySet*)
 * @brief
 *
 * @param distribution_density_law_type
 * @param distribution_data_type
 * @param distribution_value
 * @param name
 * @param min_mode_max
 * @param faultArray
 */
	void setDistributionFootWall(int distribution_density_law_type,int distribution_data_type,
			double distribution_value, std::vector< std::string > &name, std::vector< double >& min_mode_max,
			FaultArraySet *faultArray);
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
	void setDistributionArray(int distribution_density_law_type,int distribution_array_type,
			int distribution_parameter_type, int distribution_data_type,
			double distribution_value, std::vector< std::string > &name, std::vector< double >& min_mode_max,
			FaultArraySet *faultArray);
	/**
	* @fn std::vector<std::string> splitLine(std::istringstream&, char)
	 * @brief
	 *
	 * @param my_stream
	 * @param delim
	 * @return
	 */
	std::vector<std::string> splitLine(std::istringstream& my_stream, char delim);

	/**
		 * @fn void setDistributionArrayEchelonSpaceFault(int, double, std::vector<std::string>&, std::vector<double>&, FaultArraySet*)
	 * @brief
	 *
	 * @param distribution_data_type
	 * @param distribution_value
	 * @param name
	 * @param min_max_mode
	 * @param faultArray
	 */
	void setDistributionArrayEchelonSpaceFault(int distribution_data_type, double distribution_value,std::vector< std::string >& name,std::vector< double >& min_max_mode,FaultArraySet *faultArray) ;
/**
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
	}else if(value=="FRACTURE_DENSITY"){
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
	}else if(value=="FRACTURE_DENSITY_CORRIDORS"){
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
	}else if (value == "FRACTURE_DENSITY_CORRIDORS"){
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
	}else if(value=="FRACTURE_CLUSTERS_ZONE_FRACTURES_HANGING"){
		cat_index = 5;
	}else if(value=="FRACTURE_CLUSTERS_ZONE_FRACTURES_FOOTWALL"){
		cat_index = 6;
	}else if(value=="FRACTURE_CLUSTERS_ZONE_FRACTURES"){
		cat_index = 7;
	}else if(value=="FRACTURE_CLUSTERS_CORRIDOR_FRACTURES"){
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
	}else if(value=="FAULT_ARRAY_COMMON"){
	   cat_index = 2;
	}else if(value=="FAULT_ARRAY_FAMILY_1"){
		cat_index = 3;
	}else if(value=="FAULt_ARRAY_FAMILY_2"){
		cat_index = 4;
	}else if(value=="FRACTURE_CLUSTERS_ZONE_CLUSTERS"){
		cat_index = 5;
	}else if(value=="FRACTURE_CLUSTERS_ZONE_FRACTURES"){
		cat_index = 6;
	}else if(value=="FRACTURE_CLUSTERS_CORRIDOR_CLUSTERS"){
		cat_index = 7;
	}else if(value=="FRACTURE_CLUSTERS_CORRIDOR_FRACTURES"){
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
};
// inline methods


 }
#endif  // RemFracResFileManagement_h
