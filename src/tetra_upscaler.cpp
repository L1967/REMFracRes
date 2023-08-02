// tetra_upscaler.cpp : Defines the entry point for the console application.

#include <primitive_mesh/varglob.h>
#include <primitive_mesh/contrainte_geo.h>
#include <numeric/model.h>
#include <iomanip>
#include <iostream>
#include <fstream>


int main(int argc,  char * argv[])
{

	std::cout<< "start intersection " << std::endl;

	// Definition du model
	fractures_intersect::model* hydro_model = new fractures_intersect::model(); ;

	hydro_model->save_meshed_plan("test_mesh.ts");

	delete hydro_model;

	std::cout<< "end intersection " << std::endl;
	return 0;
}

