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



#include <remfracres/FracResUnit.h>

namespace remfracres
{

FracResUnit::FracResUnit(){

    name_ = "Unit_0";
	unit_age_= "MA";
	unit_duration_= "KY";
	unit_distance_= "M";
	unit_altitude_= "M";
	unit_production_duration_= "KY";
	unit_erosion_duration_= "KY";
	unit_subsidence_duration_= "M";
	unit_river_sediment_duration_= "YEAR";
	unit_river_sediment_volume_= "0.0";
	unit_temperature_= "°C";
	unit_slope_= "DEGRE";

}

FracResUnit::FracResUnit(const FracResUnit& from){

	 name_ = from.name_;
		unit_age_= from.unit_age_;
		unit_duration_= from.unit_duration_;
		unit_distance_= from.unit_distance_;
		unit_altitude_= from.unit_altitude_;
		unit_production_duration_= from.unit_production_duration_;
		unit_erosion_duration_= from.unit_erosion_duration_;
		unit_subsidence_duration_= from.unit_subsidence_duration_;
		unit_river_sediment_duration_= from.unit_river_sediment_duration_;
		unit_river_sediment_volume_= from.unit_river_sediment_volume_;
		unit_temperature_= from.unit_temperature_;
		unit_slope_= from.unit_slope_;

}

FracResUnit& FracResUnit::operator =(const FracResUnit& from){

	 name_ = from.name_;
	unit_age_= from.unit_age_;
	unit_duration_= from.unit_duration_;
	unit_distance_= from.unit_distance_;
	unit_altitude_= from.unit_altitude_;
	unit_production_duration_= from.unit_production_duration_;
	unit_erosion_duration_= from.unit_erosion_duration_;
	unit_subsidence_duration_= from.unit_subsidence_duration_;
	unit_river_sediment_duration_= from.unit_river_sediment_duration_;
	unit_river_sediment_volume_= from.unit_river_sediment_volume_;
	unit_temperature_= from.unit_temperature_;
	unit_slope_= from.unit_slope_;

	return *this;
}

void FracResUnit::copy(const FracResUnit& from){


	 name_ = from.name_;
	unit_age_= from.unit_age_;
	unit_duration_= from.unit_duration_;
	unit_distance_= from.unit_distance_;
	unit_altitude_= from.unit_altitude_;
	unit_production_duration_= from.unit_production_duration_;
	unit_erosion_duration_= from.unit_erosion_duration_;
	unit_subsidence_duration_= from.unit_subsidence_duration_;
	unit_river_sediment_duration_= from.unit_river_sediment_duration_;
	unit_river_sediment_volume_= from.unit_river_sediment_volume_;
	unit_temperature_= from.unit_temperature_;
	unit_slope_= from.unit_slope_;


}
} // namespace FracResUnit
