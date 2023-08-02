/*
 * tetra.h
 *
 *  Created on: Nov 4, 2014
 *      Author: pascal
 */

#ifndef TETRA_H_
#define TETRA_H_




namespace fractures_intersect {


//---------------------------------------------------
// definition des numero de noeud pour les 6 tetra
// d'un element hexa
//---------------------------------------------------

const int hexa_tetra[6][4]={6,0,2,3, // element tetra for a hexa
						   7,3,0,6,
						   0,6,4,7,
						   6,0,1,2,
						   0,1,5,6,
						   0,5,4,6};
//---------------------------------------------------
// definition des numero de noeud par face
// pour un element tétraèdrique
//---------------------------------------------------

const int tetra_face[4][3]={3,2,1, // tetra face
					       0,2,3,
						   3,1,0,
						   0,1,2};


const int tetra_face_unstructured[4][3]={ 0,2,1, // tetra face
						   0,3,2,
						   0,1,3,
						   1,2,3};

//---------------------------------------------------
// definition voisinage pour chaque tetra
// d'un element hexa -1 i+1 , -2 i-1 , -3 j+1 , -4 j-1 , -5 k+1 , -6 k-1
//---------------------------------------------------

const int hexa_tetra_neigbour[6][4]={-6,-3,1,3, // hexa tetra neighbour
						   0,2,-3,-2,
						   -5,-2,1,5,
						   -6,-1,0,4,
						   -1,5,3,-4,
						   -5,2,4,-4};

//---------------------------------------------------
// dead tetra with common node
//---------------------------------------------------

const int hexa_dead_tetra[4][2]={2,5, // for common node 0 - 4  dead tetra
						   4,-1,
						   0,3,
						   1,-1};





class tetra
{

public :

	tetra* t_neighbour_[4];

	double t_face_data_[4];

	int t_face_index_[4];

	int t_face_index_global_[4];

	bool t_face_visited_[4];

	bool t_face_border_[4];

    int t_face_discontinuity_[4];

	int t_edge_conduit_[6];

	int t_nodes_[4];

	bool dead_;

	int active_index_;
	int cell_hexa_index_;

	int geode_block_index_;

	int geode_polyhedron_id_;

	bool border_;

	double permx_;
	double permy_;
	double permz_;
	double facies_;
	double porosity_;
	double transmissivity_;

public :

	tetra (int node_1, int node_2, int node_3, int node_4, int active_index,int cell_index)
	{
		t_nodes_[0] = node_1;
		t_nodes_[1] = node_2;
		t_nodes_[2] = node_3;
		t_nodes_[3] = node_4;
		dead_ = false;
		for (int i = 0; i < 4; i++) {
			t_neighbour_[i] = 0;
			t_face_data_[i] = 1;
			t_face_index_[i] = -1;
			t_face_visited_[i] = false;
			t_face_index_global_[i] = -1;
			t_face_border_[4]=false;
			t_face_discontinuity_[4]=-1;
		}
		for (int j = 0; j < 6; j++) {
			t_edge_conduit_[j] = -1;
		}

		active_index_ = active_index;
		cell_hexa_index_ = cell_index;
		geode_block_index_=-1;
		geode_polyhedron_id_=-1;
		border_ = false;
	};

	tetra(bool dead, int cell_index)
	{
		dead_ = dead;

		for (int i = 0; i < 4; i++) {
			t_nodes_[i] = -1;
			t_neighbour_[i] = 0;
			t_face_data_[i] = 1;
			t_face_index_[i] = -1;
			t_face_visited_[i] = false;
			t_face_index_global_[i] = -1;
			t_face_border_[4]=false;
			t_face_discontinuity_[4]=-1;
		}
		for (int j = 0; j < 6; j++) {
			t_edge_conduit_[j] = -1;
		}
		active_index_ = -1;
		cell_hexa_index_ = cell_index;
		geode_block_index_=-1;
		geode_polyhedron_id_=-1;
		border_ = false;
	};

	int find_tetra_face_index(tetra* tetra_to_find)
	{
		for (int i = 0; i < 4; i++)
		{
			if (tetra_to_find == t_neighbour_[i]) return i;
		}

		return -1;
	}


	void clear_neigbours()
	{
		for (int i = 0; i < 4; i++) {
			t_neighbour_[i] = 0;
			t_face_index_[i] = -1;
			t_face_visited_[i] = false;
			t_face_index_global_[i] = -1;
			t_face_border_[4]=false;
			t_face_discontinuity_[4]=-1;
		}
		for (int j = 0; j < 6; j++) {
			t_edge_conduit_[j] = -1;
		}
	}

	void clear_face_index()
	{
		for (int i = 0; i < 4; i++) {
			t_face_index_[i] = -1;
			t_face_visited_[i] = false;
			t_face_index_global_[i] = -1;
			t_face_border_[4]=false;
			t_face_discontinuity_[4]=-1;

		}
		for (int j = 0; j < 6; j++) {
			t_edge_conduit_[j] = -1;
		}
	}

	int face_size() {return 4;};

	int tetra_size() {return 4;};

	int node_size() {return 4;};

	int neighbour_size() {return 4;};

};


} /* namespace common_solver */
#endif /* TETRA_H_ */
