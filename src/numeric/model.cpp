
//---------------------------- model.cpp  ---------------------------
// 
//  Description: numeric model 
//
//  Author: Pascal Siegel 
//
//	Date: 2006
//
//---------------------------- model.cpp  ---------------------------



// local includes 
#include <numeric/model.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <map>
#include <primitive_mesh/fault_array.h>
#include <primitive_mesh/fault_arrays_set_geo.h>

namespace fractures_intersect {

//------------------------------------------
// create the mesh  
//------------------------------------------

model::model()
	{

/*
	Point3D ref_pt;

	std::shared_ptr<Cboite_rect> new_box = std::shared_ptr<Cboite_rect>(new Cboite_rect ());
	new_box->Set_point_ref(ref_pt);
	new_box->Set_dx (100.0);
	new_box->Set_dy (100.0);
	new_box->Set_dz (100.0);
	new_box->Set_ind_materiaux (-1);
	new_box->Set_type_poly(Zone_affinement);
	new_box->Set_taille_element(10);
	new_box->Set_nom("cell_box");
	new_box->Actualise_plans();
	geo_cont_.Add_primitive(new_box);

	std::vector<std::shared_ptr<Cfaille> > tab_fractures_set0;
	std::vector<std::shared_ptr<Cfaille> > tab_fractures_set1;

	tab_fractures_set0.push_back(std::shared_ptr<Cfaille>(new Cfaille(0,Point3D(500.0,600.0,500.0), 0.0, 45.0, 200.0, 200.0, 0.01, 20.0)));
	tab_fractures_set0.push_back(std::shared_ptr<Cfaille>(new Cfaille(0,Point3D(500.0,750.0,500.0), 90.0, 45.0, 400.0, 500.0, 0.05, 10.0)));

	std::shared_ptr<Cfaille> faille_1 = std::shared_ptr<Cfaille>(new Cfaille(0,Point3D(480.0,500.0,500.0), 0.0, 90.0, 1000.0, 1000.0, 0.02, 50.0));

	tab_fractures_set1.push_back(std::shared_ptr<Cfaille>(new Cfaille(0,Point3D(600.0,750.0,500.0), 90.0, 45.0, 400.0, 500.0, 0.05, 10.0)));

	std::vector<FracResIntersectionTypeEnum> tab_type_intersection(7,ABUTTING);

	geo_cont_.add_fracture_set(std::shared_ptr<fracture_set_geo> (new fracture_set_geo(0, FAULTS, 12154, tab_type_intersection,tab_fractures_set0)));
	geo_cont_.add_fracture_set(std::shared_ptr<fracture_set_geo> (new fracture_set_geo(1, FAULTS, 12154, tab_type_intersection,tab_fractures_set1)));


	std::shared_ptr<fracture_set_geo> stress_zone_set(new fracture_set_geo(0, FAULTS, 12154, tab_type_intersection));

	stress_zone_set->add_fracture_and_create_cylinder(faille_1);


	std::shared_ptr<Cfaille> faille = std::shared_ptr<Cfaille>(new Cfaille(0,Point3D(300.0,300.0,500.0), 70.0, 45.0, 500.0, 800.0, 0.01, 20.0));


	stress_zone_set->add_fracture(std::shared_ptr<Cfaille>(new Cfaille(0,Point3D(300.0,300.0,500.0), 70.0, 45.0, 500.0, 800.0, 0.01, 20.0)));

	bool ducon = stress_zone_set->in_all_stress_zone__abutting_3d(faille, 12);


	geo_cont_.intersect_all_fracture_set(10);

	std::vector<Point3D> tab_point;
	std::vector<std::vector<int>> tab_triangle;
	std::vector<property<float>> tab_triangle_face_property;

	geo_cont_.get_fracture_set_surface_mesh(0, tab_point,tab_triangle,tab_triangle_face_property);
	save_meshed_vtk_plan("faille_mesh.vtk", tab_point, tab_triangle, tab_triangle_face_property);

	tab_point.clear();
	tab_triangle.clear();
	tab_triangle_face_property.clear();

	geo_cont_.get_fracture_set_surface_mesh(1, tab_point,tab_triangle,tab_triangle_face_property);
	save_meshed_vtk_plan("faille_mesh_1.vtk", tab_point, tab_triangle, tab_triangle_face_property);

	tab_point.clear();
	tab_triangle.clear();
	tab_triangle_face_property.clear();

	stress_zone_set->get_fracture_set_surface_mesh(tab_point,tab_triangle,tab_triangle_face_property);
	save_meshed_vtk_plan("faille_mesh_cyl.vtk", tab_point, tab_triangle, tab_triangle_face_property);

	tab_point.clear();
	tab_triangle.clear();
	tab_triangle_face_property.clear();


	stress_zone_set->get_fracture_set_protection_zone_surface_mesh(tab_point,tab_triangle, tab_triangle_face_property);
	save_meshed_vtk_plan("cyl_mesh.vtk", tab_point,  tab_triangle, tab_triangle_face_property,true);

	tab_point.clear();
	tab_triangle.clear();
	tab_triangle_face_property.clear();

	std::shared_ptr<fault_arrays_set_geo> fault_array_zone_set(new fault_arrays_set_geo(0,2568));
//	fault_array_zone_set->add_array(std::shared_ptr<fault_array>(new fault_array(0,Point3D(300.0,450.0,500.0), 70.0, 65.0, 500.0, 800.0, 150,10)));

//	fault_array_zone_set->add_array(std::shared_ptr<fault_array>(new fault_array(0,Point3D(300.0,450.0,500.0), 90.0, 156.0, 250, 500.0, 260,10)));

	std::shared_ptr<fault_array> test_array = std::shared_ptr<fault_array>(new fault_array(0,Point3D(100.0,250.0,500.0), 45.0, 45.0, 500.0, 1000.0, 562,10));
	fault_array_zone_set->add_array(test_array);

	test_array->tab_fractures_set_.push_back(std::shared_ptr<fracture_set_geo> (new fracture_set_geo(0, FAULTS, 12154, tab_type_intersection)));

	test_array->tab_fractures_set_[0]->add_fracture(std::shared_ptr<Cfaille>(new Cfaille(0,test_array->generate_point(0.2,0.5,0.5), 90.0, 30.0, 200.0, 300.0, 0.01, 20.0)));
	test_array->tab_fractures_set_[0]->add_fracture(std::shared_ptr<Cfaille>(new Cfaille(0,test_array->generate_point(0.4,0.5,0.5), 90.0, 30.0, 200.0, 300.0, 0.01, 20.0)));
	test_array->tab_fractures_set_[0]->add_fracture(std::shared_ptr<Cfaille>(new Cfaille(0,test_array->generate_point(0.6,0.5,0.5), 90.0, 30.0, 200.0, 300.0, 0.01, 20.0)));
	test_array->tab_fractures_set_[0]->add_fracture(std::shared_ptr<Cfaille>(new Cfaille(0,test_array->generate_point(0.8,0.5,0.5), 90.0, 30.0, 200.0, 300.0, 0.01, 20.0)));


	test_array->tab_fractures_set_.push_back(std::shared_ptr<fracture_set_geo> (new fracture_set_geo(1, FAULTS, 12154, tab_type_intersection)));

	test_array->tab_fractures_set_[1]->add_fracture(std::shared_ptr<Cfaille>(new Cfaille(0,test_array->generate_point(0.5,0.2,0.5), 30.0, 0.0, 200.0, 300.0, 0.01, 20.0)));
	test_array->tab_fractures_set_[1]->add_fracture(std::shared_ptr<Cfaille>(new Cfaille(0,test_array->generate_point(0.5,0.4,0.5), 30.0, 0.0, 200.0, 300.0, 0.01, 20.0)));
	test_array->tab_fractures_set_[1]->add_fracture(std::shared_ptr<Cfaille>(new Cfaille(0,test_array->generate_point(0.5,0.6,0.5), 30.0, 0.0, 200.0, 300.0, 0.01, 20.0)));
	test_array->tab_fractures_set_[1]->add_fracture(std::shared_ptr<Cfaille>(new Cfaille(0,test_array->generate_point(0.5,0.8,0.5), 30.0, 0.0, 200.0, 300.0, 0.01, 20.0)));


	test_array->tab_fractures_set_.push_back(std::shared_ptr<fracture_set_geo> (new fracture_set_geo(2, FAULTS, 12154, tab_type_intersection)));

	test_array->tab_fractures_set_[2]->add_fracture(std::shared_ptr<Cfaille>(new Cfaille(0,test_array->generate_point(0.5,0.5,0.2), 70.0, 0.0, 200.0, 300.0, 0.01, 20.0)));
	test_array->tab_fractures_set_[2]->add_fracture(std::shared_ptr<Cfaille>(new Cfaille(0,test_array->generate_point(0.5,0.5,0.4), 70.0, 0.0, 200.0, 300.0, 0.01, 20.0)));
	test_array->tab_fractures_set_[2]->add_fracture(std::shared_ptr<Cfaille>(new Cfaille(0,test_array->generate_point(0.5,0.5,0.6), 70.0, 0.0, 200.0, 300.0, 0.01, 20.0)));
	test_array->tab_fractures_set_[2]->add_fracture(std::shared_ptr<Cfaille>(new Cfaille(0,test_array->generate_point(0.5,0.5,0.8), 70.0, 0.0, 200.0, 300.0, 0.01, 20.0)));



	fault_array_zone_set->get_array_surface_mesh(tab_point,tab_triangle, tab_triangle_face_property);
	save_meshed_vtk_plan("fault_array.vtk", tab_point,  tab_triangle, tab_triangle_face_property,true);

	tab_point.clear();
	tab_triangle.clear();
	tab_triangle_face_property.clear();

	fault_array_zone_set->get_fracture_set_surface_mesh(tab_point,tab_triangle, tab_triangle_face_property);
	save_meshed_vtk_plan("fault_array_fractrure_set.vtk", tab_point,  tab_triangle, tab_triangle_face_property,true);
*/
}



//----------------------------------------------------------
//	read the fractures data   
//----------------------------------------------------------

bool  model::save_meshed_plan(std::string filename)
{
	std::fstream os(filename.data(),std::ios_base::out);

	if (os.is_open() == false )
	{
		std::cout << "Could not save the mesh plan file" <<  filename << std::endl;
		return false;
	}


	for (size_t i=0;i<geo_cont_.tab_plans.size();i++) // boucle sur les plans
	{
		os<< geo_cont_.tab_plans[i];
	}

	for (size_t i=0;i<geo_cont_.tab_polyedres.size();i++) // boucle sur les plans
	{
		os<< geo_cont_.tab_polyedres[i];
	}


	os<<"END"<<std::endl;

	os.close();
	return true;

}

//-----------------------------------------------------
//save the meshed fractures data
//-----------------------------------------------------

bool save_meshed_vtk_plan(std::string filename, std::vector<Point3D>& tab_point, std::vector<std::vector<int>>& tab_triangle,
		std::vector<property<float>>& tab_triangle_face_property, bool is_binary)
{
	std::fstream os;

	if(is_binary) {
		os.open(filename.data(), std::ios::out | std::ios::binary);
	}
	else {
		os.open(filename.data(),std::ios_base::out);
	}

	if (os.is_open() == false )
	{
		std::cout << "Could not save vtk the mesh plan file" <<  filename << std::endl;
		return false;
	}else{
		std::cout << "Save vtk the mesh plan file" <<  filename << std::endl;
	}

	os<<"# vtk DataFile Version 2.0"<<std::endl;
	os<<"Cube faille"<<std::endl;
	if (is_binary) os<<"BINARY"<<std::endl;
	else os<<"ASCII"<<std::endl;

	os<<"DATASET POLYDATA"<<std::endl;
	os<<"POINTS " <<  tab_point.size() <<" double" <<std::endl;


	for(size_t i=0;i<tab_point.size();i++)
	{
		if (is_binary) {


			double x = tab_point[i][0];
			double y = tab_point[i][1];
		    double z = tab_point[i][2];

			SwapEnd(x);
			SwapEnd(y);
			SwapEnd(z);

			os.write(reinterpret_cast<char*>(&x), sizeof(double));
			os.write(reinterpret_cast<char*>(&y), sizeof(double));
			os.write(reinterpret_cast<char*>(&z), sizeof(double));
		}
		else {
			os<< tab_point[i] <<std::endl;
		}
	}

	if (is_binary) os<<std::endl;

	os<<"TRIANGLE_STRIPS " << tab_triangle.size() << " " << tab_triangle.size()*4<<std::endl;

	for(size_t i=0;i<tab_triangle.size();i++)
	{
		if (is_binary) {

			 int n_size = 3;
			 SwapEnd(n_size);

			 os.write(reinterpret_cast<char*>(&n_size), sizeof(int));

			 int  n1 = tab_triangle[i][0];
			 int  n2 = tab_triangle[i][1];
			 int  n3 = tab_triangle[i][2];

			 SwapEnd(n1);
			 SwapEnd(n2);
			 SwapEnd(n3);

			 os.write(reinterpret_cast<char*>(&n1), sizeof(int));
			 os.write(reinterpret_cast<char*>(&n2), sizeof(int));
			 os.write(reinterpret_cast<char*>(&n3), sizeof(int));

		}
		else {
			os<< "3 " << tab_triangle[i][0] <<" " << tab_triangle[i][1] << " " << tab_triangle[i][2]<<std::endl;
		}
	}

	if (is_binary) os<<std::endl;


	os<<"CELL_DATA " << tab_triangle.size() <<std::endl;
	for (size_t pi = 0; pi < tab_triangle_face_property.size();pi++)
	{
		os<<"SCALARS "<< tab_triangle_face_property[pi].name_ <<" float "<<std::endl;
		os<<"LOOKUP_TABLE default"<<std::endl;
		for(size_t i=0;i<tab_triangle_face_property[pi].size();i++)
		{
			float value = tab_triangle_face_property[pi].value(i);

			if (is_binary) {
				SwapEnd(value);
				os.write(reinterpret_cast<char*>(&value), sizeof(float));
			}
			else os<< value <<" ";

		}
		os<<std::endl;

	}


	os.close();
	return true;
}




}


	







