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



#include <primitive_mesh/FracResBeddingLocalisation.h>

namespace fractures_intersect
{


FracResBeddingLocalisation::FracResBeddingLocalisation()
{

	index_= -1;

	index_polyedre_1_ = -1;

	index_polyedre_2_= -1;

}

FracResBeddingLocalisation::FracResBeddingLocalisation( int indice, int polyhedron_index1, int polyhedron_index2, geode::uuid id_block )
{

	index_= indice;

	index_polyedre_1_ = polyhedron_index1;

	index_polyedre_2_= polyhedron_index2;
	block_id_ = id_block;

}

FracResBeddingLocalisation::FracResBeddingLocalisation(const FracResBeddingLocalisation& from)
{

	index_= from.index_;

	index_polyedre_1_ = from.index_polyedre_1_;

	index_polyedre_2_= from.index_polyedre_2_;

	block_id_ = from.block_id_;
}


FracResBeddingLocalisation& FracResBeddingLocalisation::operator =(const FracResBeddingLocalisation& from){

	index_= from.index_;

	index_polyedre_1_ = from.index_polyedre_1_;

	index_polyedre_2_= from.index_polyedre_2_;

	block_id_ = from.block_id_;

	return *this;
}

void FracResBeddingLocalisation::copy(const FracResBeddingLocalisation& from){

	index_= from.index_;

	index_polyedre_1_ = from.index_polyedre_1_;

	index_polyedre_2_= from.index_polyedre_2_;

	block_id_ = from.block_id_;
}

} // namespace FracResBeddingLocalisation
