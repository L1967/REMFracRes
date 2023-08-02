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



#include <remfracres/FracResBeddingProperties.h>

namespace remfracres
{


FracResBeddingProperties::FracResBeddingProperties()
{
	stage_name_ = "None";

	fracres_hypothesis_name_= "None";

	name_ = "None";

	slope_parameter_name_ = "None";

	slope_data_type_ = FracResDataTypeEnum::CONSTANT;

	slope_unit_name_ = "None";

	slope_data_value_ = 0.0;

	slope_ratio_data_value_ = 1.0;

	bed_parrallel_slope_distribution_name_ = "None";

	bed_parrallel_slope_distribution_type_name_ = "None";

	bed_parrallel_discrete_property_name_ = "None";

	bed_parrallel_slope_prop_name_ = "None";

	region_name_ = "None";

	propertIndex_ = 0;

	is_region_ = false;

	bedParallelPropertyActive_ = false;

	bedding_list_.clear();
	min_mode_max_.resize(3);
	min_mode_max_[0]=0;
	min_mode_max_[1]=1;
	min_mode_max_[2]=0;

   index_fracture_set_geo_begin_ = 0;

   nb_fracture_set_geo_ = 0;

   output_directory_path_ = "";

   output_prefixe_name_ = "";

}

FracResBeddingProperties::FracResBeddingProperties(const FracResBeddingProperties& from)
{
	stage_name_ = from.stage_name_;

	fracres_hypothesis_name_= from.fracres_hypothesis_name_;

	name_ = from.name_;

	slope_parameter_name_ = from.slope_parameter_name_;

	slope_data_type_ = from.slope_data_type_;

	slope_unit_name_ = from.slope_unit_name_;

	slope_data_value_ = from.slope_data_value_;

	slope_ratio_data_value_ = from.slope_ratio_data_value_;

	bed_parrallel_slope_distribution_name_ = from.bed_parrallel_slope_distribution_name_;

	bed_parrallel_slope_distribution_type_name_ = from.bed_parrallel_slope_distribution_type_name_;

	min_mode_max_ = from.min_mode_max_;

	bed_parrallel_discrete_property_name_ = from.bed_parrallel_discrete_property_name_;

	bed_parrallel_slope_prop_name_ = from.bed_parrallel_slope_prop_name_;

	region_name_ = from.region_name_;

	propertIndex_ = from.propertIndex_;

	is_region_ = from.is_region_;

	bedParallelPropertyActive_ = from.bedParallelPropertyActive_;

	bedding_list_ = from.bedding_list_;

	index_fracture_set_geo_begin_ =  from.index_fracture_set_geo_begin_ ;

	nb_fracture_set_geo_ = from.nb_fracture_set_geo_;

	output_directory_path_ = from.output_directory_path_;

	output_prefixe_name_ = from.output_prefixe_name_;

}


FracResBeddingProperties& FracResBeddingProperties::operator =(const FracResBeddingProperties& from){

	stage_name_ = from.stage_name_;

	fracres_hypothesis_name_= from.fracres_hypothesis_name_;

	name_ = from.name_;

	slope_parameter_name_ = from.slope_parameter_name_;

	slope_data_type_ = from.slope_data_type_;

	slope_unit_name_ = from.slope_unit_name_;

	slope_data_value_ = from.slope_data_value_;

	slope_ratio_data_value_ = from.slope_ratio_data_value_;

	bed_parrallel_slope_distribution_name_ = from.bed_parrallel_slope_distribution_name_;

	bed_parrallel_slope_distribution_type_name_ = from.bed_parrallel_slope_distribution_type_name_;

	min_mode_max_ = from.min_mode_max_;

	bed_parrallel_discrete_property_name_ = from.bed_parrallel_discrete_property_name_;

	bed_parrallel_slope_prop_name_ = from.bed_parrallel_slope_prop_name_;

	region_name_ = from.region_name_;

	propertIndex_ = from.propertIndex_;

	is_region_ = from.is_region_;

	bedParallelPropertyActive_ = from.bedParallelPropertyActive_;

	bedding_list_ = from.bedding_list_;

	index_fracture_set_geo_begin_ =  from.index_fracture_set_geo_begin_ ;

	nb_fracture_set_geo_ = from.nb_fracture_set_geo_;

	output_directory_path_ = from.output_directory_path_;

	output_prefixe_name_ = from.output_prefixe_name_;

	return *this;
}

void FracResBeddingProperties::copy(const FracResBeddingProperties& from){

	stage_name_ = from.stage_name_;

	fracres_hypothesis_name_= from.fracres_hypothesis_name_;

	name_ = from.name_;

	slope_parameter_name_ = from.slope_parameter_name_;

	slope_data_type_ = from.slope_data_type_;

	slope_unit_name_ = from.slope_unit_name_;

	slope_data_value_ = from.slope_data_value_;

	slope_ratio_data_value_ = from.slope_ratio_data_value_;

	bed_parrallel_slope_distribution_name_ = from.bed_parrallel_slope_distribution_name_;

	bed_parrallel_slope_distribution_type_name_ = from.bed_parrallel_slope_distribution_type_name_;

	min_mode_max_ = from.min_mode_max_;

	bed_parrallel_discrete_property_name_ = from.bed_parrallel_discrete_property_name_;

	bed_parrallel_slope_prop_name_ = from.bed_parrallel_slope_prop_name_;

	region_name_ = from.region_name_;

	propertIndex_ = from.propertIndex_;

	is_region_ = from.is_region_;

	bedParallelPropertyActive_ = from.bedParallelPropertyActive_;

	bedding_list_ = from.bedding_list_;

	index_fracture_set_geo_begin_ =  from.index_fracture_set_geo_begin_ ;

	nb_fracture_set_geo_ = from.nb_fracture_set_geo_;

	output_directory_path_ = from.output_directory_path_;

	output_prefixe_name_ = from.output_prefixe_name_;
}

} // namespace FracResBeddingProperties
