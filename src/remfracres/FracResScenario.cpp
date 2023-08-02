/*
 * Copyright (c) 2019 - 2022 Geode-solutions
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 */



#include <remfracres/FracResScenario.h>

namespace remfracres
{

FracResScenario::FracResScenario(){

	scenario_file_name_="Scenario_00";
    name_ = "Scenario_0";
    path_result_directory_ = "Result_";
	scenarioIndex_ = 0;
	estimate_Fracture_Number_ =0;
	modeling_ = false;
	numberOfRealisation_ = 1;
	scenario_stage_hypothesis_index_vector_.clear();
	success_list_.clear();
	model_fractures_list_.clear();
	nb_failles_ = 0;
	total_simulation_time_list_.clear();
	total_writting_output_time_list_.clear();
	evaluate_fractures_number_list_.clear();



}

FracResScenario::FracResScenario(const FracResScenario& from){

	name_ = from.name_;
	scenario_file_name_ = from.scenario_file_name_;
	scenarioIndex_ = from.scenarioIndex_;
	estimate_Fracture_Number_ = from.estimate_Fracture_Number_;
	modeling_ = from.modeling_;
	numberOfRealisation_ = from.numberOfRealisation_;
	scenario_stage_hypothesis_index_vector_ = from.scenario_stage_hypothesis_index_vector_;
	path_result_directory_ = from.path_result_directory_;
	success_list_ = from.success_list_;
	model_fractures_list_ = from.model_fractures_list_;
	nb_failles_ = from.nb_failles_;
	total_simulation_time_list_ = from.total_simulation_time_list_;
	total_writting_output_time_list_ = from.total_writting_output_time_list_;
	evaluate_fractures_number_list_ = from.evaluate_fractures_number_list_;

}


FracResScenario& FracResScenario::operator =(const FracResScenario& from){

	name_ = from.name_;
	scenario_file_name_ = from.scenario_file_name_;
	scenarioIndex_ = from.scenarioIndex_;
	estimate_Fracture_Number_ = from.estimate_Fracture_Number_;
	modeling_ = from.modeling_;
	numberOfRealisation_ = from.numberOfRealisation_;
	scenario_stage_hypothesis_index_vector_ = from.scenario_stage_hypothesis_index_vector_;
	path_result_directory_ = from.path_result_directory_;
	success_list_ = from.success_list_;
	model_fractures_list_ = from.model_fractures_list_;
	nb_failles_ = from.nb_failles_;
	total_simulation_time_list_ = from.total_simulation_time_list_;
	total_writting_output_time_list_ = from.total_writting_output_time_list_;
	evaluate_fractures_number_list_ = from.evaluate_fractures_number_list_;

	return *this;
}

void FracResScenario::copy(const FracResScenario& from){

	name_ = from.name_;
	scenario_file_name_ = from.scenario_file_name_;
	scenarioIndex_ = from.scenarioIndex_;
	estimate_Fracture_Number_ = from.estimate_Fracture_Number_;
	modeling_ = from.modeling_;
	numberOfRealisation_ = from.numberOfRealisation_;
	scenario_stage_hypothesis_index_vector_ = from.scenario_stage_hypothesis_index_vector_;
	path_result_directory_ = from.path_result_directory_;
	success_list_ = from.success_list_;
	model_fractures_list_ = from.model_fractures_list_;
	nb_failles_ = from.nb_failles_;
	total_simulation_time_list_ = from.total_simulation_time_list_;
	total_writting_output_time_list_ = from.total_writting_output_time_list_;
	evaluate_fractures_number_list_ = from.evaluate_fractures_number_list_;

}
} // namespace FracResScenario
