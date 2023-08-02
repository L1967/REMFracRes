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



#include <remfracres/FracResHypothesis.h>

namespace remfracres
{

FracResHypothesis::FracResHypothesis(){

	stage_name_="stage00";
    name_ = "Hypothesis_0";
	beddingIsRegion_=false;
	beddingRegionName_ = "None";
	beddingPropertyDiscreteName_ = "None";
	hypothesisIndex_ = 0;
    faultArrayList_.clear();
    fractureList_.clear();
    beddingList_.clear();
    bed_parallel_List_.clear();
	type_dfn_.clear();
	index_dfn_.clear();
	faultCorridorList_.clear();
	faultZoneList_.clear();
	evaluate_fractures_number_total_ = 0;
	nb_bed_parallel_active_ = 0;
	nb_bed_interface_active_ = 0;
	nb_single_fault_active_ = 0;
	nb_fault_array_active_ = 0;
	nb_cluster_corridors_active_ = 0;
	nb_cluster_fault_zone_active_ = 0;


}

FracResHypothesis::FracResHypothesis(const FracResHypothesis& from){

	name_ = from.name_;
	stage_name_ = from.stage_name_;
	hypothesisIndex_ = from.hypothesisIndex_;
	beddingIsRegion_ = from.beddingIsRegion_;
	beddingRegionName_ = from.beddingRegionName_;
	beddingPropertyDiscreteName_ = from.beddingPropertyDiscreteName_;
	faultArrayList_ = from.faultArrayList_;
	fractureList_ = from.fractureList_;
	beddingList_ = from.beddingList_;
    bed_parallel_List_= from.bed_parallel_List_;
	type_dfn_ = from.type_dfn_;
	index_dfn_ = from.index_dfn_;
	faultCorridorList_ = from.faultCorridorList_;
	faultZoneList_ = from.faultZoneList_;
	evaluate_fractures_number_total_ = from.evaluate_fractures_number_total_;
	nb_bed_parallel_active_ = from.nb_bed_parallel_active_;
	nb_bed_interface_active_ = from.nb_bed_interface_active_;
	nb_single_fault_active_ = from.nb_single_fault_active_;
	nb_fault_array_active_ = from.nb_fault_array_active_;
	nb_cluster_corridors_active_ = from.nb_cluster_corridors_active_;
	nb_cluster_fault_zone_active_ = from.nb_cluster_fault_zone_active_;
}


FracResHypothesis& FracResHypothesis::operator =(const FracResHypothesis& from){

	name_ = from.name_;
	stage_name_ = from.stage_name_;
	hypothesisIndex_ = from.hypothesisIndex_;
	beddingIsRegion_ = from.beddingIsRegion_;
	beddingRegionName_ = from.beddingRegionName_;
	beddingPropertyDiscreteName_ = from.beddingPropertyDiscreteName_;
	faultArrayList_ = from.faultArrayList_;
	fractureList_ = from.fractureList_;
	beddingList_ = from.beddingList_;
    bed_parallel_List_= from.bed_parallel_List_;
	type_dfn_ = from.type_dfn_;
	index_dfn_ = from.index_dfn_;
	faultCorridorList_ = from.faultCorridorList_;
	faultZoneList_ = from.faultZoneList_;
	evaluate_fractures_number_total_ = from.evaluate_fractures_number_total_;
	nb_bed_parallel_active_ = from.nb_bed_parallel_active_;
	nb_bed_interface_active_ = from.nb_bed_interface_active_;
	nb_single_fault_active_ = from.nb_single_fault_active_;
	nb_fault_array_active_ = from.nb_fault_array_active_;
	nb_cluster_corridors_active_ = from.nb_cluster_corridors_active_;
	nb_cluster_fault_zone_active_ = from.nb_cluster_fault_zone_active_;
	return *this;
}

void FracResHypothesis::copy(const FracResHypothesis& from){

	name_ = from.name_;
	stage_name_ = from.stage_name_;
	hypothesisIndex_ = from.hypothesisIndex_;
	beddingIsRegion_ = from.beddingIsRegion_;
	beddingRegionName_ = from.beddingRegionName_;
	beddingPropertyDiscreteName_ = from.beddingPropertyDiscreteName_;
	faultArrayList_ = from.faultArrayList_;
	fractureList_ = from.fractureList_;
	beddingList_ = from.beddingList_;
    bed_parallel_List_= from.bed_parallel_List_;
	type_dfn_ = from.type_dfn_;
	index_dfn_ = from.index_dfn_;
	faultCorridorList_ = from.faultCorridorList_;
	faultZoneList_ = from.faultZoneList_;
	evaluate_fractures_number_total_ = from.evaluate_fractures_number_total_;
	nb_bed_parallel_active_ = from.nb_bed_parallel_active_;
	nb_bed_interface_active_ = from.nb_bed_interface_active_;
	nb_single_fault_active_ = from.nb_single_fault_active_;
	nb_fault_array_active_ = from.nb_fault_array_active_;
	nb_cluster_corridors_active_ = from.nb_cluster_corridors_active_;
	nb_cluster_fault_zone_active_ = from.nb_cluster_fault_zone_active_;
}

void FracResHypothesis::update_active_dfn_number(){

	nb_bed_parallel_active_ = bed_parallel_List_.size();
	nb_bed_interface_active_ = beddingList_.size();

	nb_single_fault_active_ = 0;
	for( FractureSet* setFrac : fractureList_){
		if(setFrac->is_active_) nb_single_fault_active_++;
	}

	nb_fault_array_active_ = 0;
	for( FaultArraySet*  setArray : faultArrayList_){
		if(setArray->is_active_) nb_fault_array_active_++;
	}

	nb_cluster_corridors_active_ = 0;
	for( FaultArraySet*  setArrayCorridor : faultCorridorList_){
		if(setArrayCorridor->is_active_) nb_cluster_corridors_active_++;
	}

	nb_cluster_fault_zone_active_ = 0;
	for( FaultArraySet*  setArrayFaultZone : faultZoneList_){
		if(setArrayFaultZone->is_active_) nb_cluster_fault_zone_active_++;
	}
}
} // namespace FracResHypothesis
