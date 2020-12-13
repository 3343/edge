/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (alex.breuer AT uni-jena.de)
 *
 * @section LICENSE
 * Copyright (c) 2020, Friedrich Schiller University Jena
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 * Grid holding input data.
 **/
#include "Grid.h"

edge_v::io::Grid::Grid( io::Hdf5 const * i_reader ) {
  m_reader = i_reader;
}

void edge_v::io::Grid::free() {
  if( m_data != nullptr ) delete[] m_data;
}

edge_v::io::Grid::~Grid() {
  // free memory if allocated
  free();
}

void edge_v::io::Grid::init( t_idx             i_nPts,
                             double   const (* i_pts)[3],
                             io::Hdf5 const  * i_reader,
                             float           * o_data ) {
  // get the number of x- and y-values
  t_idx l_nx = i_reader->nVas("/x");
  t_idx l_ny = i_reader->nVas("/y");
  t_idx l_nz = i_reader->nVas("/z");

  // check sanity
  EDGE_V_CHECK_GT( l_nz, 0 );
  EDGE_V_CHECK_EQ( l_nz % l_nx, 0 );
  EDGE_V_CHECK_EQ( l_nz % l_ny, 0 );

  // read the coordinates
  double * l_x = new double[ l_nx ];
  double * l_y = new double[ l_ny ];
  float  * l_z = new float[  l_nz ];

  i_reader->get( "/x",
                 l_x );

  i_reader->get( "/y",
                 l_y );

  i_reader->get( "/z",
                 l_z );

  // iterate over all pts and query
  for( t_idx l_pt = 0; l_pt < i_nPts; l_pt++ ) {
    // get grid id
    t_idx l_id0 = std::lower_bound( l_y, l_y + l_ny, i_pts[l_pt][1] ) - l_y;
    t_idx l_id1 = std::lower_bound( l_x, l_x + l_nx, i_pts[l_pt][0] ) - l_x;
    t_idx l_id = l_id0 * l_nx + l_id1;

    // check that the point is within the grid
    EDGE_V_CHECK_GT( i_pts[l_pt][0], l_x[0] );
    EDGE_V_CHECK_GT( i_pts[l_pt][1], l_y[0] );
    EDGE_V_CHECK_LT( l_id1, l_nx );
    EDGE_V_CHECK_LT( l_id0, l_ny );

    o_data[l_pt] = l_z[l_id];
  }

  // free temporary memory
  delete[] l_x;
  delete[] l_y;
  delete[] l_z;
}

void edge_v::io::Grid::init( t_idx           i_nPts,
                             double const (* i_pts)[3] ) {

  // free memory if allocated
  free();

  // read bathymetry data
  m_data = new float[i_nPts];
  init( i_nPts,
        i_pts,
        m_reader,
        m_data );
}