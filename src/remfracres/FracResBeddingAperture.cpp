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



#include <remfracres/FracResBeddingAperture.h>

namespace remfracres
{


FracResBeddingAperture::FracResBeddingAperture()
{

	aperture_name_="None";

	parameter_name_ = "Aperture";

	name_ = "bed_interface";

	unit_name_ = "None";

	discrete_property_name_ = "None";

	dataType_ = FracResDataTypeEnum::CONSTANT;

	value_ = 0.0;

	distribution_name_ = "None";

	discrete_property_selected_values_.clear();

	index_=0;

	min_mode_max_.resize(3);
	min_mode_max_[0]=0;
	min_mode_max_[1]=1;
	min_mode_max_[2]=0;
	is_active_ = true;
	output_directory_path_ = "";

	output_prefixe_name_ = "";
	index_begin_ =0;
	index_end_ =0;

}

FracResBeddingAperture::FracResBeddingAperture(const FracResBeddingAperture& from)
{

	name_ = from.name_;
	aperture_name_= from.aperture_name_;
	unit_name_ = from.unit_name_;
	discrete_property_name_ = from.discrete_property_name_;
	dataType_ = from.dataType_;
	value_ = from.value_;
	distribution_name_ = from.distribution_name_;
	discrete_property_selected_values_ = from.discrete_property_selected_values_;
	parameter_name_ = from.parameter_name_;
	index_ = from.index_;
	min_mode_max_ = from.min_mode_max_;
	is_active_ = from.is_active_;
	output_directory_path_ = from.output_directory_path_;
	output_prefixe_name_ = from.output_prefixe_name_;
	index_begin_ = from.index_begin_;
	index_end_  = from.index_end_;
}


FracResBeddingAperture& FracResBeddingAperture::operator =(const FracResBeddingAperture& from){

	name_ = from.name_;
	aperture_name_= from.aperture_name_;
	unit_name_ = from.unit_name_;
	discrete_property_name_ = from.discrete_property_name_;
	dataType_ = from.dataType_;
	value_ = from.value_;
	distribution_name_ = from.distribution_name_;
	discrete_property_selected_values_ = from.discrete_property_selected_values_;
	parameter_name_ = from.parameter_name_;
	index_ = from.index_;
	min_mode_max_ = from.min_mode_max_;
	is_active_ = from.is_active_;
	output_directory_path_ = from.output_directory_path_;
	output_prefixe_name_ = from.output_prefixe_name_;
	index_begin_ = from.index_begin_;
	index_end_  = from.index_end_;
	return *this;
}

void FracResBeddingAperture::copy(const FracResBeddingAperture& from){

	name_ = from.name_;
	aperture_name_= from.aperture_name_;
	unit_name_ = from.unit_name_;
	discrete_property_name_ = from.discrete_property_name_;
	dataType_ = from.dataType_;
	value_ = from.value_;
	distribution_name_ = from.distribution_name_;
	discrete_property_selected_values_ = from.discrete_property_selected_values_;
	parameter_name_ = from.parameter_name_;
	index_ = from.index_;
	min_mode_max_ = from.min_mode_max_;
	is_active_ = from.is_active_;
	output_directory_path_ = from.output_directory_path_;
	output_prefixe_name_ = from.output_prefixe_name_;
	index_begin_ = from.index_begin_;
	index_end_  = from.index_end_;
}

} // namespace FracResBeddingAperture
