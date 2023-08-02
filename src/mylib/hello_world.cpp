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

#include <mylib/hello_world.h>

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


namespace remfracres
{
    using namespace geode;
    using namespace std;
    using namespace fractures_intersect;

    bool hello_world()
    {
        geode::Logger::info( "Hello Geode World!" );
        return true;
    }
    bool test_fractures_boolean() {

    	geode::OpenGeodeGeosciencesGeosciences::initialize();
        geode::OpenGeodeGeosciencesIOGeosciences::initialize();
    	geode::OpenGeodeIOModel::initialize();
    	geode::OpenGeodeIOMesh::initialize();
    	std::string filename = "LS2.";
    	std::string prefixe_fracture = "FractureSet_";
    	std::string surface_filename =  absl::StrCat(prefixe_fracture,filename);
    	auto surface = geode::TriangulatedSurface3D::create(
    	geode::OpenGeodeTriangulatedSurface3D::impl_name_static() );
    	auto builder = geode::TriangulatedSurfaceBuilder3D::create( *surface );

    	const geode::StructuralModel& model = geode::load_structural_model(absl::StrCat(remfracres::data_path, filename, "lso"));
    	remfracres::FractureSet obj;
    	obj.fracture_set_index_=1;
    	obj.set_fractures_number(2);
    	obj.density_=0.01;
    	for(int i=0; i < 10;i++){
    		obj.length1_.push_back(0.5*(i+1));
    		obj.length2_.push_back(0.1*(i+1));
    		obj.Dip_.push_back(10*(i+1));
    		obj.Orientation_.push_back(180*(i+1));
    		obj.aperture_.push_back(0.5*(i+1));
    	}
    	obj.clean_vertices_and_triangles_array();
    	obj.set_generation_box(model.bounding_box());
    	obj.fractures_generation(model,false);
    	obj.init_surface_geometry(obj.fracture_set_index_,*surface, *builder);

    	const auto output_file_native = absl::StrCat( obj.directory_path_,surface_filename, "vtp");

    	const auto output_file_surface_native = absl::StrCat( obj.directory_path_,surface_filename, "dxf");

    	// auto poly_surface = geode::convert_surface_mesh_into_polygonal_surface( *surface );


    	//geode::save_polygonal_surface(*poly_surface, output_file_surface_native);

    	geode::save_triangulated_surface( *surface, output_file_native );


    	return true;
    }
} // namespace mymodule
