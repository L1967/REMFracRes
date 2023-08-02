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

//----------------------------RemFracResFileManagement.cpp  ---------------------------
//
//  RemFracResFileManagement contains the parameters for reading FracRes Sismage File parameter
//
//
//  Author:  Pascal Siegel - Rabah Namar - Olivier Jaquet
//
//	Date: 2002-2022
//
//----------------------------  RemFracResFileManagement.cpp  ---------------------------
// local includes
#include <remfracres/RemFracResFileManagement.h>
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

// lib includes 

// system includes                       

namespace remfracres {

//-----------------------------------------------
// default constructor of the class RemFracResFileManagement
//-----------------------------------------------
/**
 * @fn  RemFracResFileManagement()
 * @brief
 *
 */
RemFracResFileManagement::RemFracResFileManagement() {
	name_ = "simultion_0";
	seed_ = 1337;
	use_global_outpout_directory_path_ = false;
	outpout_directory_file_ = remfracres::data_path;
	manage_bed_interface_ = false;
	fracresFileKeyWord_ = new FracResFileKeyWord();
	fracResUnit_ = new FracResUnit();
	fracResStageList_.clear();

}

//---------------------------------------------
//  copy constructor of the class RemFracResFileManagement
//---------------------------------------------
/**
 * @fn  RemFracResFileManagement(const RemFracResFileManagement&)
 * @brief
 *
 * @param from
 */
RemFracResFileManagement::RemFracResFileManagement(
		const RemFracResFileManagement &from) {
	name_ = from.name_;
	seed_ = from.seed_;

	use_global_outpout_directory_path_ =
			from.use_global_outpout_directory_path_;
	outpout_directory_file_ = from.outpout_directory_file_;

	manage_bed_interface_ = from.manage_bed_interface_;
	fracresFileKeyWord_ = from.fracresFileKeyWord_;
	fracResUnit_ = from.fracResUnit_;
	fracResStageList_.clear();

}

//-------------------------------------------
// destructor of the class RemFracResFileManagement
//-------------------------------------------
/**
 * @fn  ~RemFracResFileManagement()
 * @brief
 *
 */
RemFracResFileManagement::~RemFracResFileManagement() {
	fracResStageList_.clear();
}

//-----------------------------------------------
// assignment operator of the class RemFracResFileManagement
//----------------------------------------------
/**
 * @fn RemFracResFileManagement operator =&(const RemFracResFileManagement&)
 * @brief
 *
 * @param from
 * @return
 */
RemFracResFileManagement& RemFracResFileManagement::operator=(
		const RemFracResFileManagement &from) {
	name_ = from.name_;
	seed_ = from.seed_;

	use_global_outpout_directory_path_ =
			from.use_global_outpout_directory_path_;
	outpout_directory_file_ = from.outpout_directory_file_;

	manage_bed_interface_ = from.manage_bed_interface_;
	fracresFileKeyWord_ = from.fracresFileKeyWord_;
	fracResUnit_ = from.fracResUnit_;
	return *this;
}

//--------------------------------------------------------
//      copy function 
//--------------------------------------------------------
/**
 * @fn void copy(const RemFracResFileManagement&)
 * @brief
 *
 * @param from
 */
void RemFracResFileManagement::copy(const RemFracResFileManagement &from) {

	use_global_outpout_directory_path_ =
			from.use_global_outpout_directory_path_;
	outpout_directory_file_ = from.outpout_directory_file_;
	manage_bed_interface_ = from.manage_bed_interface_;
	fracresFileKeyWord_ = from.fracresFileKeyWord_;
	fracResUnit_ = from.fracResUnit_;
	seed_ = from.seed_;

}

/**
 * @fn std::string enumCategoryToString(fractures_intersect::FracResIntersectionCategoryEnum)
 * @brief
 *
 * @param cat
 * @return
 */

std::string RemFracResFileManagement::enumCategoryToString(
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
		buffer = "FAULTS";
		break;
	case fractures_intersect::FracResIntersectionCategoryEnum::FAULT_ARRAYS:
		buffer = "FAULT_ARRAYS";
		break;
	case fractures_intersect::FracResIntersectionCategoryEnum::JOINTS:
		buffer = "JOINTS";
		break;
	case fractures_intersect::FracResIntersectionCategoryEnum::FRACTURE_CLUSTERS_CORRIDORS_ZONES:
		buffer = "FRACTURE_CLUSTERS_CORRIDORS_ZONES";
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
std::string RemFracResFileManagement::enumIntersectionToString(
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
bool RemFracResFileManagement::openFracres(std::string filename) {

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
std::vector<std::string> RemFracResFileManagement::splitLine(
		std::istringstream &my_stream, char delim) {
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
bool RemFracResFileManagement::manageLine(std::fstream &fin) {
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
		} else if (tokenVector[0]
				== fracresFileKeyWord_->manage_bed_interface) {
			bool b = (tokenVector[1] == "true");
			manage_bed_interface_ = b;

		} else if (tokenVector[0] == fracresFileKeyWord_->unit_begin) {
			manageUnit(fin);
		} else if (tokenVector[0] == fracresFileKeyWord_->stage_begin) {
			FracResStage fracResStage = manageStage(fin);
			fracResStageList_.push_back(fracResStage);
		} else if (tokenVector[0] == fracresFileKeyWord_->scenario_begin) {
			manageScenario(fin);
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
void RemFracResFileManagement::manageUnit(std::fstream &fin) {

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
FracResStage RemFracResFileManagement::manageStage(std::fstream &fin) {
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
FracResHypothesis RemFracResFileManagement::manageRepHypothesis(
		std::fstream &fin) {
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
void RemFracResFileManagement::manageMatrixProp(std::fstream &fin) {
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
void RemFracResFileManagement::manageOutputs(std::fstream &fin) {
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
			outpout_directory_file_ = tokenVector[1];
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
void RemFracResFileManagement::manageScenario(std::fstream &fin) {
	std::string contenu_ligne;
	std::string comment;
	std::string separator;

	while (getline(fin, contenu_ligne)) // tant que l'on peut mettre la ligne dans "contenu"
	{
		std::istringstream my_stream;

		cout << contenu_ligne << endl;  // on l'affiche
		my_stream.str(contenu_ligne);

		std::vector<std::string> tokenVector = splitLine(my_stream, ';');
		if (tokenVector[0] == fracresFileKeyWord_->outputs_density) {

		} else if (tokenVector[0]
				== fracresFileKeyWord_->scenario_is_modeling) {
			bool test = (tokenVector[1] == "true");
		} else if (tokenVector[0] == fracresFileKeyWord_->scenario_is_qc) {

		} else if (tokenVector[0] == fracresFileKeyWord_->scenario_nb_real) {
			int test2 = stoi(tokenVector[1]);
		} else if (tokenVector[0]
				== fracresFileKeyWord_->scenario_rep_hyp_index) {
			int test0 = stoi(tokenVector[1]);
		} else if (tokenVector[0]
				== fracresFileKeyWord_->scenario_nb_estimate_fractures) {

		} else if (tokenVector[0] == fracresFileKeyWord_->scenario_name) {

		} else if (tokenVector[0] == fracresFileKeyWord_->scenario_index) {
			int test8 = stoi(tokenVector[1]);
		} else if (tokenVector[0]
				== fracresFileKeyWord_->scenario_stage_hyp_link) {

		} else if (tokenVector[0]
				== fracresFileKeyWord_->scenario_stages_hyp_map) {

		} else if (tokenVector[0] == fracresFileKeyWord_->scenario_end) {

			break;
		}
		my_stream.clear();
	}

}
/**
 * @fn void manageStageHypScenarioKeyVal(std::fstream&)
 * @brief
 *
 * @param fin
 */
void RemFracResFileManagement::manageStageHypScenarioKeyVal(std::fstream &fin) {
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
void RemFracResFileManagement::manageUpscalingRep(std::fstream &fin) {
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
void RemFracResFileManagement::manageRepresentationObject(std::fstream &fin) {
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
void RemFracResFileManagement::manageRunParameter(std::fstream &fin) {
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
FracResFracture RemFracResFileManagement::manageFracture(std::fstream &fin) {
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
FracResHypothesis RemFracResFileManagement::manageHypothesis(
		std::fstream &fin) {
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
void RemFracResFileManagement::manageBeddingStageHypothesis(std::fstream &fin,
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
			fracResHypothesis.beddingRegionName_ = tokenVector[1];
		} else if (tokenVector[0]
				== fracresFileKeyWord_->hyp_bed_prop_discrete_name) {
			fracResHypothesis.beddingPropertyDiscreteName_ = tokenVector[1];
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
FracResBeddingAperture* RemFracResFileManagement::manageBeddingAperture(
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
			fracResBeddingAperture->discrete_property_name_ = tokenVector[1];
		} else if (tokenVector[0]
				== fracresFileKeyWord_->bed_aperture_property_name) {
			fracResBeddingAperture->discrete_property_name_ = tokenVector[1];
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

void RemFracResFileManagement::split(const string &chaine,
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
void RemFracResFileManagement::manageClassicStageHypothesis(std::fstream &fin,
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
void RemFracResFileManagement::manageFractureProperty(std::fstream &fin,
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

				faultArray->heightCorrelationProperty_name_ = tokenVector[1];
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

				faultArray->apertureCorrelationProperty_name_ = tokenVector[1];
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

				faultArray->faultCoreProperty_name_ = tokenVector[1];
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
void RemFracResFileManagement::setGeometryDipAzimuth(int geometry_data_type,
		double geometry_value, std::vector<std::string> &name,
		const std::vector<double> &min_max_mode, FractureSet *faultSingle) {

	if (std::isnan(geometry_value))
		geometry_value = 1e-6;
	if (geometry_data_type == 0) {
		faultSingle->distribution_azimut_type_index_ = 1;
		faultSingle->distr_azimut_val_ = geometry_value;
	} else if (geometry_data_type == 1) {
		faultSingle->distribution_azimut_type_index_ = 0;
		if (name[1] == "normal") {
			faultSingle->azimuth_type_
					== FractureSet::distribution_type::normal;
		} else if (name[1] == "triangulated") {
			faultSingle->azimuth_type_
					== FractureSet::distribution_type::triangulated;
			faultSingle->minmaxmode_values_[8] = min_max_mode[1];
		} else if (name[1] == "uniform") {
			faultSingle->azimuth_type_
					== FractureSet::distribution_type::uniform;
		} else if (name[1] == "lognormal") {
			faultSingle->azimuth_type_
					== FractureSet::distribution_type::lognormal;
		}

		faultSingle->minmaxmode_values_[6] = min_max_mode[0];
		faultSingle->minmaxmode_values_[7] = min_max_mode[2];
	} else {
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
void RemFracResFileManagement::setGeometryDipAngle(int geometry_data_type,
		double geometry_value, std::vector<std::string> &name,
		const std::vector<double> &min_max_mode, FractureSet *faultSingle) {
	if (std::isnan(geometry_value))
		geometry_value = 90;
	if (geometry_data_type == 0) {
		faultSingle->distribution_angle_type_index_ = 1;
		faultSingle->distr_angle_val_ = geometry_value;
	} else if (geometry_data_type == 1) {
		faultSingle->distribution_angle_type_index_ = 0;

		if (name[1] == "normal") {
			faultSingle->dip_angle_type_ =
					FractureSet::distribution_type::normal;
		} else if (name[1] == "triangulated") {

			faultSingle->dip_angle_type_ =
					FractureSet::distribution_type::triangulated;
			faultSingle->minmaxmode_values_[11] = min_max_mode[1];

		} else if (name[1] == "uniform") {
			faultSingle->dip_angle_type_ =
					FractureSet::distribution_type::uniform;
		} else if (name[1] == "lognormal") {
			faultSingle->dip_angle_type_ =
					FractureSet::distribution_type::lognormal;
		}

		faultSingle->minmaxmode_values_[9] = min_max_mode[0];
		faultSingle->minmaxmode_values_[10] = min_max_mode[2];

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
void RemFracResFileManagement::setGeometryLength(int geometry_data_type,
		double geometry_value, std::vector<std::string> &name,
		const std::vector<double> &min_max_mode, FractureSet *faultSingle) {
	if (std::isnan(geometry_value))
		geometry_value = 1e-6;
	if (geometry_data_type == 0) {
		faultSingle->distribution_length1_type_index_ = 1;
		faultSingle->distr_length1_val_ = geometry_value;
	} else if (geometry_data_type == 1) {
		faultSingle->distribution_length1_type_index_ = 0;

		if (name[1] == "normal") {
			faultSingle->length1_type_ = FractureSet::distribution_type::normal;
		} else if (name[1] == "triangulated") {

			faultSingle->length1_type_ =
					FractureSet::distribution_type::triangulated;
			faultSingle->minmaxmode_values_[2] = min_max_mode[1];

		} else if (name[1] == "uniform") {
			faultSingle->length1_type_ =
					FractureSet::distribution_type::uniform;
		} else if (name[1] == "lognormal") {
			faultSingle->length1_type_ =
					FractureSet::distribution_type::lognormal;
		}

		faultSingle->minmaxmode_values_[0] = min_max_mode[0];
		faultSingle->minmaxmode_values_[1] = min_max_mode[2];

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
void RemFracResFileManagement::setGeometryHeight(int geometry_data_type,
		double geometry_value, std::vector<std::string> &name,
		const std::vector<double> &min_max_mode, FractureSet *faultSingle) {
	if (std::isnan(geometry_value))
		geometry_value = 1e-6;
	if (geometry_data_type == 0) {
		faultSingle->distribution_length2_type_index_ = 1;
		faultSingle->distr_length2_val_ = geometry_value;
	} else if (geometry_data_type == 1) {
		faultSingle->distribution_length2_type_index_ = 0;

		if (name[1] == "normal") {
			faultSingle->length2_type_ = FractureSet::distribution_type::normal;
		} else if (name[1] == "triangulated") {

			faultSingle->length2_type_ =
					FractureSet::distribution_type::triangulated;
			faultSingle->minmaxmode_values_[5] = min_max_mode[1];

		} else if (name[1] == "uniform") {
			faultSingle->length2_type_ =
					FractureSet::distribution_type::uniform;
		} else if (name[1] == "lognormal") {
			faultSingle->length2_type_ =
					FractureSet::distribution_type::lognormal;
		}

		faultSingle->minmaxmode_values_[3] = min_max_mode[0];
		faultSingle->minmaxmode_values_[4] = min_max_mode[2];

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
void RemFracResFileManagement::setGeometryAperture(int geometry_data_type,
		double geometry_value, std::vector<std::string> &name,
		const std::vector<double> &min_max_mode, FractureSet *faultSingle) {
	if (std::isnan(geometry_value))
		geometry_value = 1e-6;
	if (geometry_data_type == 0) {
		faultSingle->distribution_aperture_type_index_ = 1;
		faultSingle->distr_aperture_val_ = geometry_value;
	} else if (geometry_data_type == 1) {
		faultSingle->distribution_aperture_type_index_ = 0;

		if (name[1] == "normal") {
			faultSingle->aperture_type_ =
					FractureSet::distribution_type::normal;
		} else if (name[1] == "triangulated") {

			faultSingle->aperture_type_ =
					FractureSet::distribution_type::triangulated;
			faultSingle->minmaxmode_values_[14] = min_max_mode[1];

		} else if (name[1] == "uniform") {
			faultSingle->aperture_type_ =
					FractureSet::distribution_type::uniform;
		} else if (name[1] == "lognormal") {
			faultSingle->aperture_type_ =
					FractureSet::distribution_type::lognormal;
		}

		faultSingle->minmaxmode_values_[12] = min_max_mode[0];
		faultSingle->minmaxmode_values_[13] = min_max_mode[2];

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
void RemFracResFileManagement::setGeometrySinglefault(
		int geometry_parameter_type, int geometry_data_type,
		double geometry_value, std::vector<std::string> &name,
		const std::vector<double> &min_max_mode, FractureSet *faultSingle) {

	if (geometry_parameter_type == 2) {
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
void RemFracResFileManagement::setGeometrySingleFaultCorridor(
		int geometry_parameter_type, int geometry_data_type,
		double geometry_value, std::vector<std::string> &name,
		const std::vector<double> &min_max_mode, FractureSet *faultSingle) {

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
void RemFracResFileManagement::setGeometrySingleFaultZone(
		int geometry_parameter_type, int geometry_data_type,
		double geometry_value, std::vector<std::string> &name,
		const std::vector<double> &min_max_mode, FaultArraySet *faultArray) {

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
void RemFracResFileManagement::setGeometryAzimuthArray(int geometry_data_type,
		double geometry_value, std::vector<std::string> &name,
		const std::vector<double> &min_max_mode, FaultArraySet *faultArray) {

	if (geometry_data_type == 0) {
		if (std::isnan(geometry_value))
			geometry_value = 180;
		faultArray->distribution_azimut_type_index_ = 1;
		faultArray->distr_azimut_val_ = geometry_value;
	} else if (geometry_data_type == 1) {
		faultArray->distribution_azimut_type_index_ = 0;
		if (name[1] == "normal") {
			faultArray->azimuth_type_ =
					FaultArraySet::distribution_type::normal;
		} else if (name[1] == "triangulated") {
			faultArray->azimuth_type_ =
					FaultArraySet::distribution_type::triangulated;
			faultArray->minmaxmode_values_[8] = min_max_mode[1];
		} else if (name[1] == "uniform") {
			faultArray->azimuth_type_ =
					FaultArraySet::distribution_type::uniform;
		} else if (name[1] == "lognormal") {
			faultArray->azimuth_type_ =
					FaultArraySet::distribution_type::lognormal;
		}

		faultArray->minmaxmode_values_[6] = min_max_mode[0];
		faultArray->minmaxmode_values_[7] = min_max_mode[2];
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
void RemFracResFileManagement::setGeometryAngleArray(int geometry_data_type,
		double geometry_value, std::vector<std::string> &name,
		const std::vector<double> &min_max_mode, FaultArraySet *faultArray) {
	if (std::isnan(geometry_value))
		geometry_value = 45;
	if (geometry_data_type == 0) {
		faultArray->distribution_angle_type_index_ = 1;
		faultArray->distr_angle_val_ = geometry_value;
	} else if (geometry_data_type == 1 || geometry_data_type == 3) {
		faultArray->distribution_angle_type_index_ = 0;

		if (name[1] == "normal") {
			faultArray->dip_angle_type_ =
					FaultArraySet::distribution_type::normal;
		} else if (name[1] == "triangulated") {

			faultArray->dip_angle_type_ =
					FaultArraySet::distribution_type::triangulated;
			faultArray->minmaxmode_values_[11] = min_max_mode[1];

		} else if (name[1] == "uniform") {
			faultArray->dip_angle_type_ =
					FaultArraySet::distribution_type::uniform;
		} else if (name[1] == "lognormal") {
			faultArray->dip_angle_type_ =
					FaultArraySet::distribution_type::lognormal;
		}

		faultArray->minmaxmode_values_[9] = min_max_mode[0];
		faultArray->minmaxmode_values_[10] = min_max_mode[2];

	} else {
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
void RemFracResFileManagement::setGeometryLengthArray(int geometry_data_type,
		double geometry_value, std::vector<std::string> &name,
		const std::vector<double> &min_max_mode, FaultArraySet *faultArray) {
	if (std::isnan(geometry_value))
		geometry_value = 1;
	if (geometry_data_type == 0) {
		faultArray->distribution_length1_type_index_ = 1;
		faultArray->distr_length1_val_ = geometry_value;
	} else if (geometry_data_type == 1) {
		faultArray->distribution_length1_type_index_ = 0;

		if (name[1] == "normal") {
			faultArray->length1_type_ =
					FaultArraySet::distribution_type::normal;
		} else if (name[1] == "triangulated") {

			faultArray->length1_type_ =
					FaultArraySet::distribution_type::triangulated;
			faultArray->minmaxmode_values_[2] = min_max_mode[1];

		} else if (name[1] == "uniform") {
			faultArray->length1_type_ =
					FaultArraySet::distribution_type::uniform;
		} else if (name[1] == "lognormal") {
			faultArray->length1_type_ =
					FaultArraySet::distribution_type::lognormal;
		}

		faultArray->minmaxmode_values_[0] = min_max_mode[0];
		faultArray->minmaxmode_values_[1] = min_max_mode[2];

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
void RemFracResFileManagement::setGeometryHeightArray(int geometry_data_type,
		double geometry_value, std::vector<std::string> &name,
		const std::vector<double> &min_max_mode, FaultArraySet *faultArray) {
	if (std::isnan(geometry_value))
		geometry_value = 1;
	if (geometry_data_type == 0) {
		faultArray->distribution_length2_type_index_ = 1;
		faultArray->distr_length2_val_ = geometry_value;
	} else if (geometry_data_type == 1) {
		faultArray->distribution_length2_type_index_ = 0;

		if (name[1] == "normal") {
			faultArray->length2_type_ =
					FaultArraySet::distribution_type::normal;

		} else if (name[1] == "triangulated") {

			faultArray->length2_type_ =
					FaultArraySet::distribution_type::triangulated;
			faultArray->minmaxmode_values_[5] = min_max_mode[1];

		} else if (name[1] == "uniform") {
			faultArray->length2_type_ =
					FaultArraySet::distribution_type::uniform;
		} else if (name[1] == "lognormal") {
			faultArray->length2_type_ =
					FaultArraySet::distribution_type::lognormal;
		}

		faultArray->minmaxmode_values_[3] = min_max_mode[0];
		faultArray->minmaxmode_values_[4] = min_max_mode[2];

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

void RemFracResFileManagement::setGeometryDamageHangingArray(
		int geometry_data_type, double geometry_value,
		std::vector<std::string> &name, const std::vector<double> &min_max_mode,
		FaultArraySet *faultArray) {
	if (std::isnan(geometry_value))
		geometry_value = 1e-6;
	if (geometry_data_type == 0) {
		faultArray->distribution_damage_zone_hanging_index_ = 1;
		faultArray->distr_damage_zone_hanging_val_ = geometry_value;
	} else if (geometry_data_type == 1) {
		faultArray->distribution_damage_zone_hanging_index_ = 0;

		if (name[1] == "normal") {
			faultArray->damage_zone_hanging_type_ =
					FaultArraySet::distribution_type::normal;
		} else if (name[1] == "triangulated") {

			faultArray->damage_zone_hanging_type_ =
					FaultArraySet::distribution_type::triangulated;
			faultArray->minmaxmode_damage_zone_hanging_values_[1] =
					min_max_mode[1];

		} else if (name[1] == "uniform") {
			faultArray->damage_zone_hanging_type_ =
					FaultArraySet::distribution_type::uniform;
		} else if (name[1] == "lognormal") {
			faultArray->damage_zone_hanging_type_ =
					FaultArraySet::distribution_type::lognormal;
		}

		faultArray->minmaxmode_damage_zone_hanging_values_[0] = min_max_mode[0];
		faultArray->minmaxmode_damage_zone_hanging_values_[2] = min_max_mode[2];

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
void RemFracResFileManagement::setGeometryDamageFootwallArray(
		int geometry_data_type, double geometry_value,
		std::vector<std::string> &name, const std::vector<double> &min_max_mode,
		FaultArraySet *faultArray) {
	if (std::isnan(geometry_value))
		geometry_value = 1e-6;
	if (geometry_data_type == 0) {
		faultArray->distribution_damage_zone_footwall_index_ = 1;
		faultArray->distr_damage_zone_footwall_val_ = geometry_value;
	} else if (geometry_data_type == 1) {
		faultArray->distribution_damage_zone_footwall_index_ = 0;

		if (name[1] == "normal") {
			faultArray->damage_zone_footwall_type_ =
					FaultArraySet::distribution_type::normal;
		} else if (name[1] == "triangulated") {

			faultArray->damage_zone_footwall_type_ =
					FaultArraySet::distribution_type::triangulated;
			faultArray->minmaxmode_damage_zone_footwall_values_[1] =
					min_max_mode[1];

		} else if (name[1] == "uniform") {
			faultArray->damage_zone_footwall_type_ =
					FaultArraySet::distribution_type::uniform;
		} else if (name[1] == "lognormal") {
			faultArray->damage_zone_footwall_type_ =
					FaultArraySet::distribution_type::lognormal;
		}

		faultArray->minmaxmode_damage_zone_footwall_values_[0] =
				min_max_mode[0];
		faultArray->minmaxmode_damage_zone_footwall_values_[2] =
				min_max_mode[2];

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
void RemFracResFileManagement::setGeometryApertureArray(int geometry_data_type,
		double geometry_value, std::vector<std::string> &name,
		const std::vector<double> &min_max_mode, FaultArraySet *faultArray) {
	if (std::isnan(geometry_value))
		geometry_value = 1e-6;
	if (geometry_data_type == 0) {
		faultArray->distribution_aperture_type_index_ = 1;
		faultArray->distr_aperture_val_ = geometry_value;
	} else if (geometry_data_type == 1) {
		faultArray->distribution_aperture_type_index_ = 0;

		if (name[1] == "normal") {
			faultArray->aperture_type_ =
					FaultArraySet::distribution_type::normal;
		} else if (name[1] == "triangulated") {

			faultArray->aperture_type_ =
					FaultArraySet::distribution_type::triangulated;
			faultArray->minmaxmode_values_[14] = min_max_mode[1];

		} else if (name[1] == "uniform") {
			faultArray->aperture_type_ =
					FaultArraySet::distribution_type::uniform;
		} else if (name[1] == "lognormal") {
			faultArray->aperture_type_ =
					FaultArraySet::distribution_type::lognormal;
		}

		faultArray->minmaxmode_values_[12] = min_max_mode[0];
		faultArray->minmaxmode_values_[13] = min_max_mode[2];

	} else {
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
void RemFracResFileManagement::setGeometryArray(int geometry_parameter_type,
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
	} else if (geometry_parameter_type == 6) {
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
void RemFracResFileManagement::setGeometryCorridor(int geometry_parameter_type,
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
	} else if (geometry_parameter_type == 7) {
		faultArray->distribution_length2_type_index_ = 1;
		faultArray->distr_length2_val_ = geometry_value;

	} else if (geometry_parameter_type == 6) {
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

void RemFracResFileManagement::setGeometryFaultZone(int geometry_parameter_type,
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

void RemFracResFileManagement::manageTotalGeometry(std::fstream &fin,
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
void RemFracResFileManagement::setDistributionSingleFault(
		int distribution_parameter_type, int distribution_data_type,
		double distribution_value, std::vector<std::string> &name,
		std::vector<double> &min_max_mode, FractureSet *faultSingle) {

	if (distribution_parameter_type == 3) {
		if (std::isnan(distribution_value))
			distribution_value = 1e-6;

		if (distribution_data_type == 0) {
			faultSingle->density_ = distribution_value;
		} else if (distribution_data_type == 1) {
			faultSingle->density_distribution_name_ = name[0];

			if (name[1] == "normal") {
				faultSingle->density_type_
						== FractureSet::distribution_type::normal;
			} else if (name[1] == "triangulated") {

				faultSingle->density_type_ =
						FractureSet::distribution_type::triangulated;
				faultSingle->minmaxmode_density_values_[2] = min_max_mode[1];

			} else if (name[1] == "uniform") {
				faultSingle->density_type_ =
						FractureSet::distribution_type::uniform;
			} else if (name[1] == "lognormal") {
				faultSingle->density_type_ =
						FractureSet::distribution_type::lognormal;
			}

			faultSingle->minmaxmode_density_values_[0] = min_max_mode[0];
			faultSingle->minmaxmode_density_values_[1] = min_max_mode[2];

		} else {
			faultSingle->density_property_name_ = name[0];
		}
	} else if (distribution_parameter_type == 11) {
		if (std::isnan(distribution_value))
			distribution_value = 0.1;

		if (distribution_data_type == 0) {
			faultSingle->default_stress_zone_max_ = distribution_value;
		} else if (distribution_data_type == 1) {
			faultSingle->width_stress_zone_distribution_name_ = name[0];

			faultSingle->distribution_shadow_zone_type_index_ = 0;
			if (name[1] == "normal") {
				faultSingle->shadow_zone_type_
						== FractureSet::distribution_type::normal;
			} else if (name[1] == "triangulated") {

				faultSingle->shadow_zone_type_ =
						FractureSet::distribution_type::triangulated;
				faultSingle->minmaxmode_shadow_zone_values_[2] =
						min_max_mode[1];

			} else if (name[1] == "uniform") {
				faultSingle->shadow_zone_type_ =
						FractureSet::distribution_type::uniform;
			} else if (name[1] == "lognormal") {
				faultSingle->shadow_zone_type_ =
						FractureSet::distribution_type::lognormal;
			}

			faultSingle->minmaxmode_shadow_zone_values_[0] = min_max_mode[0];
			faultSingle->minmaxmode_shadow_zone_values_[1] = min_max_mode[2];

		} else {
			faultSingle->distribution_shadow_zone_type_index_ = 1;
			faultSingle->width_stress_zone_property_name_ = name[0];
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
void RemFracResFileManagement::setDistributionHangingWall(
		int distribution_parameter_type, int distribution_data_type,
		double distribution_value, std::vector<std::string> &name,
		std::vector<double> &min_max_mode, FaultArraySet *faultArray) {
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
			if (name[1] == "normal") {
				faultArray->damage_zone_hanging_type_
						== FaultArraySet::distribution_type::normal;
			} else if (name[1] == "triangulated") {

				faultArray->damage_zone_hanging_type_ =
						FaultArraySet::distribution_type::triangulated;
				faultArray->minmaxmode_damage_zone_hanging_values_[2] =
						min_max_mode[1];

			} else if (name[1] == "uniform") {
				faultArray->damage_zone_hanging_type_ =
						FaultArraySet::distribution_type::uniform;
			} else if (name[1] == "lognormal") {
				faultArray->damage_zone_hanging_type_ =
						FaultArraySet::distribution_type::lognormal;
			}

			faultArray->minmaxmode_damage_zone_hanging_values_[0] =
					min_max_mode[0];
			faultArray->minmaxmode_damage_zone_hanging_values_[1] =
					min_max_mode[2];

		} else {
			faultArray->distribution_damage_zone_hanging_index_ = 1;
			faultArray->damage_zone_hanging_property_name_ = name[0];
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
void RemFracResFileManagement::setDistributionFootWall(
		int distribution_parameter_type, int distribution_data_type,
		double distribution_value, std::vector<std::string> &name,
		std::vector<double> &min_max_mode, FaultArraySet *faultArray) {
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
			if (name[1] == "normal") {
				faultArray->damage_zone_footwall_type_
						== FaultArraySet::distribution_type::normal;
			} else if (name[1] == "triangulated") {

				faultArray->damage_zone_footwall_type_ =
						FaultArraySet::distribution_type::triangulated;
				faultArray->minmaxmode_damage_zone_footwall_values_[2] =
						min_max_mode[1];

			} else if (name[1] == "uniform") {
				faultArray->damage_zone_footwall_type_ =
						FaultArraySet::distribution_type::uniform;
			} else if (name[1] == "lognormal") {
				faultArray->damage_zone_footwall_type_ =
						FaultArraySet::distribution_type::lognormal;
			}

			faultArray->minmaxmode_damage_zone_footwall_values_[0] =
					min_max_mode[0];
			faultArray->minmaxmode_damage_zone_footwall_values_[1] =
					min_max_mode[2];

		} else {
			faultArray->distribution_damage_zone_footwall_index_ = 1;
			faultArray->damage_zone_footwall_property_name_ = name[0];
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
void RemFracResFileManagement::setDistributionArray(
		int distribution_parameter_type, int distribution_data_type,
		double distribution_value, std::vector<std::string> &name,
		std::vector<double> &min_max_mode, FaultArraySet *faultArray) {
	if (distribution_parameter_type == 5 || distribution_parameter_type == 1
			|| distribution_parameter_type == 0) {
		if (std::isnan(distribution_value))
			distribution_value = 1e-6;

		if (distribution_data_type == 0) {
			faultArray->density_ = distribution_value;
		} else if (distribution_data_type == 1) {
			faultArray->density_distribution_name_ = name[0];

			if (name[1] == "normal") {
				faultArray->density_type_
						== FaultArraySet::distribution_type::normal;
			} else if (name[1] == "triangulated") {

				faultArray->density_type_ =
						FaultArraySet::distribution_type::triangulated;
				faultArray->minmaxmode_density_values_[2] = min_max_mode[1];

			} else if (name[1] == "uniform") {
				faultArray->density_type_ =
						FaultArraySet::distribution_type::uniform;
			} else if (name[1] == "lognormal") {
				faultArray->density_type_ =
						FaultArraySet::distribution_type::lognormal;
			}

			faultArray->minmaxmode_density_values_[0] = min_max_mode[0];
			faultArray->minmaxmode_density_values_[1] = min_max_mode[2];

		} else {
			faultArray->density_property_name_ = name[0];
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
			if (name[1] == "normal") {
				faultArray->minimum_space_type_
						== FaultArraySet::distribution_type::normal;
			} else if (name[1] == "triangulated") {

				faultArray->minimum_space_type_ =
						FaultArraySet::distribution_type::triangulated;
				faultArray->minmaxmode_minimum_space_values_[2] =
						min_max_mode[1];

			} else if (name[1] == "uniform") {
				faultArray->minimum_space_type_ =
						FaultArraySet::distribution_type::uniform;
			} else if (name[1] == "lognormal") {
				faultArray->minimum_space_type_ =
						FaultArraySet::distribution_type::lognormal;
			}

			faultArray->minmaxmode_minimum_space_values_[0] = min_max_mode[0];
			faultArray->minmaxmode_minimum_space_values_[1] = min_max_mode[2];

		} else {
			faultArray->distribution_minimum_space_type_index_ = 1;
			faultArray->minimum_space_between_array_property_name_ = name[0];
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
void RemFracResFileManagement::setDistributionArrayFault(
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
		if (distribution_data_type == 0) {
			if (std::isnan(distribution_value))
				distribution_value = 1e-6;
			faultArray->fracturesSet_list_[0].density_ = distribution_value;

		} else if (distribution_data_type == 1) {
			faultArray->fracturesSet_list_[0].density_distribution_name_ =
					name[0];
			if (name[1] == "normal") {
				faultArray->fracturesSet_list_[0].density_type_
						== FractureSet::distribution_type::normal;
			} else if (name[1] == "triangulated") {

				faultArray->fracturesSet_list_[0].density_type_ =
						FractureSet::distribution_type::triangulated;
				faultArray->fracturesSet_list_[0].minmaxmode_density_values_[2] =
						min_max_mode[1];

			} else if (name[1] == "uniform") {
				faultArray->fracturesSet_list_[0].density_type_ =
						FractureSet::distribution_type::uniform;
			} else if (name[1] == "lognormal") {
				faultArray->fracturesSet_list_[0].density_type_ =
						FractureSet::distribution_type::lognormal;
			}

			faultArray->fracturesSet_list_[0].minmaxmode_density_values_[0] =
					min_max_mode[0];
			faultArray->fracturesSet_list_[0].minmaxmode_density_values_[1] =
					min_max_mode[2];

		} else {
			faultArray->fracturesSet_list_[0].density_property_name_ = name[0];
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
void RemFracResFileManagement::setDistributionArrayStressZoneFault(
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

		if (name[1] == "normal") {
			faultArray->fracturesSet_list_[0].shadow_zone_type_
					== FractureSet::distribution_type::normal;
		} else if (name[1] == "triangulated") {

			faultArray->fracturesSet_list_[0].shadow_zone_type_ =
					FractureSet::distribution_type::triangulated;
			faultArray->fracturesSet_list_[0].minmaxmode_shadow_zone_values_[2] =
					min_max_mode[1];

		} else if (name[1] == "uniform") {
			faultArray->fracturesSet_list_[0].shadow_zone_type_ =
					FractureSet::distribution_type::uniform;
		} else if (name[1] == "lognormal") {
			faultArray->fracturesSet_list_[0].shadow_zone_type_ =
					FractureSet::distribution_type::lognormal;
		}

		faultArray->fracturesSet_list_[0].minmaxmode_shadow_zone_values_[0] =
				min_max_mode[0];
		faultArray->fracturesSet_list_[0].minmaxmode_shadow_zone_values_[1] =
				min_max_mode[2];
	} else {
		faultArray->fracturesSet_list_[0].width_stress_zone_property_name_ =
				name[0];
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
void RemFracResFileManagement::setDistributionArrayEchelonSpaceFault(
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

		if (name[1] == "normal") {
			faultArray->fracturesSet_list_[0].echelon_space_type_
					== FractureSet::distribution_type::normal;
		} else if (name[1] == "triangulated") {

			faultArray->fracturesSet_list_[0].echelon_space_type_ =
					FractureSet::distribution_type::triangulated;
			faultArray->fracturesSet_list_[0].minmaxmode_echelon_space_values_[2] =
					min_max_mode[1];

		} else if (name[1] == "uniform") {
			faultArray->fracturesSet_list_[0].echelon_space_type_ =
					FractureSet::distribution_type::uniform;
		} else if (name[1] == "lognormal") {
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
void RemFracResFileManagement::setDistributionArrayFaultFamili1(
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
			if (name[1] == "normal") {
				faultArray->fracturesSet_list_[0].density_type_
						== FractureSet::distribution_type::normal;
			} else if (name[1] == "triangulated") {

				faultArray->fracturesSet_list_[0].density_type_ =
						FractureSet::distribution_type::triangulated;
				faultArray->fracturesSet_list_[0].minmaxmode_density_values_[2] =
						min_max_mode[1];

			} else if (name[1] == "uniform") {
				faultArray->fracturesSet_list_[0].density_type_ =
						FractureSet::distribution_type::uniform;
			} else if (name[1] == "lognormal") {
				faultArray->fracturesSet_list_[0].density_type_ =
						FractureSet::distribution_type::lognormal;
			}

			faultArray->fracturesSet_list_[0].minmaxmode_density_values_[0] =
					min_max_mode[0];
			faultArray->fracturesSet_list_[0].minmaxmode_density_values_[1] =
					min_max_mode[2];

		} else {
			faultArray->fracturesSet_list_[0].density_property_name_ = name[0];
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
void RemFracResFileManagement::setDistributionArray(
		int distribution_density_law_type, int distribution_array_type,
		int distribution_parameter_type, int distribution_data_type,
		double distribution_value, std::vector<std::string> &name,
		std::vector<double> &min_mode_max, FaultArraySet *faultArray) {
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
void RemFracResFileManagement::manageTotalDistribution(std::fstream &fin,
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

fractures_intersect::FracResIntersectionTypeEnum RemFracResFileManagement::manageIntersection(
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
std::vector<std::string> RemFracResFileManagement::manageDistribution(
		std::fstream &fin, int &distribution_parameter_type,
		int &distribution_array_type, int &distribution_data_type,
		int &distribution_density_law_type, double &distribution_value,
		std::vector<double> &minmax) {
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
std::vector<std::string> RemFracResFileManagement::manageGeometry(
		std::fstream &fin, int &geometry_parameter_type,
		int &geometry_array_type, int &geometry_data_type,
		double &distribution_value, std::vector<double> &min_max_mode) {
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
				name[1] = "normal";
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
				std::vector<std::string> tokenVectorValues = splitLine(
						my_stream3, ':');
				min_max_mode[0] = std::strtod(tokenVectorValues[1].c_str(),
				NULL);
				min_max_mode[1] = std::strtod(tokenVectorValues[3].c_str(),
				NULL);
				min_max_mode[2] = std::strtod(tokenVectorValues[5].c_str(),
				NULL);
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
void RemFracResFileManagement::manageBedParallel(std::fstream &fin,
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
			fracResBeddingProperties->bed_parrallel_discrete_property_name_ =
					tokenVector[1];
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

}

