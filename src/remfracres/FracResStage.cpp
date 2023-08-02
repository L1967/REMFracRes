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



#include <remfracres/FracResStage.h>

namespace remfracres
{


FracResStage::FracResStage()
{

	 name_="stage00";

	stageIndex_=0;

	beddingStage_=false;

	hypotheisisList_.clear();

	fractureList_.clear();
	evaluate_fractures_number_total_per_hypothesis_.clear();
}

FracResStage::FracResStage(const FracResStage& from)
{

	 name_ = from.name_;

	stageIndex_ = from.stageIndex_;

	beddingStage_ = from.beddingStage_;

	hypotheisisList_ = from.hypotheisisList_;

	fractureList_ = from.fractureList_;
	evaluate_fractures_number_total_per_hypothesis_ = from.evaluate_fractures_number_total_per_hypothesis_;

}


FracResStage& FracResStage::operator =(const FracResStage& from){

	 name_ = from.name_;

	stageIndex_ = from.stageIndex_;

	beddingStage_ = from.beddingStage_;

	hypotheisisList_ = from.hypotheisisList_;

	fractureList_ = from.fractureList_;
	evaluate_fractures_number_total_per_hypothesis_ = from.evaluate_fractures_number_total_per_hypothesis_;


	return *this;
}

void FracResStage::copy(const FracResStage& from){

	 name_ = from.name_;

	stageIndex_ = from.stageIndex_;

	beddingStage_ = from.beddingStage_;

	hypotheisisList_ = from.hypotheisisList_;

	fractureList_ = from.fractureList_;
	evaluate_fractures_number_total_per_hypothesis_ = from.evaluate_fractures_number_total_per_hypothesis_;

}

} // namespace FracResStage
