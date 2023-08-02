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


#include <fstream>
#include <iostream>
#include <memory>
#include <string>

// lib includes

#include <geode/basic/assert.h>
#include <geode/basic/logger.h>

#include <geode/mesh/core/tetrahedral_solid.h>
#include <geode/geometry/basic_objects/tetrahedron.h>
#include <geode/geometry/mensuration.h>
#include <geode/mesh/helpers/convert_surface_mesh.h>

#include <geode/mesh/core/geode_triangulated_surface.h>
#include <geode/mesh/core/geode_polygonal_surface.h>
#include <geode/mesh/core/polygonal_surface.h>
#include <geode/mesh/core/triangulated_surface.h>
#include <geode/mesh/io/triangulated_surface_input.h>
#include <geode/mesh/io/triangulated_surface_output.h>
#include <geode/mesh/io/polygonal_surface_input.h>
#include <geode/mesh/io/polygonal_surface_output.h>
#include <geode/mesh/builder/geode_triangulated_surface_builder.h>
#include <geode/mesh/builder/geode_polygonal_surface_builder.h>
#include <geode/mesh/builder/triangulated_surface_builder.h>
#include <geode/mesh/builder/polygonal_surface_builder.h>
#include <geode/geosciences/common.h>
#include <geode/io/geosciences/common.h>
#include <geode/io/mesh/common.h>
#include <geode/io/model/common.h>
#include <remfracres/fractureset.h>
int main()
{
	using namespace geode;
	using namespace std;
	using namespace remfracres;

	std::cout << "Hello world";

	try {
		geode::OpenGeodeGeosciencesGeosciences::initialize();
	    geode::OpenGeodeGeosciencesIOGeosciences::initialize();
		geode::OpenGeodeIOModel::initialize();
		geode::OpenGeodeIOMesh::initialize();

		remfracres::FractureSet obj;
		obj.test_All_fractures_boolean();

		Logger::info("TEST SUCCESS");

		return 0;
	} catch (...) {
		return geode::geode_lippincott();
	}

	return 0;
}
