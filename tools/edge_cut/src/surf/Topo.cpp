/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 * @author David Lenz (dlenz AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2018, Regents of the University of California
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
 * Meshing of topographic data.
 **/

#include "Topo.h"
#include "io/logging.hpp"

#include <fstream>


edge_cut::surf::Topo::Topo( std::string const & i_topoFile )
{
  std::ifstream l_ptFile( i_topoFile, std::ios::in );

  // check for NaNs in the input
  std::string l_val;
  while( l_ptFile >> l_val ) {
    EDGE_CHECK_NE( l_val, "NaN" );
  }
  l_ptFile.clear();
  l_ptFile.seekg(0, std::ios::beg);

  std::istream_iterator< TopoPoint > l_fileBegin( l_ptFile );
  std::istream_iterator< TopoPoint > l_fileEnd;

  m_delTria = new Triangulation( l_fileBegin, l_fileEnd );
}


edge_cut::surf::Topo::~Topo()
{
  delete m_delTria;
}

std::ostream &
edge_cut::surf::Topo::writeTriaToOff( std::ostream & io_os ) const
{
  typedef typename Triangulation::Vertex_handle                Vertex_handle;
  typedef typename Triangulation::Finite_vertices_iterator     Vertex_iterator;
  typedef typename Triangulation::Finite_faces_iterator        Face_iterator;

  io_os << "OFF" << std::endl;

  // outputs the number of vertices and faces
  std::size_t const l_nVerts = m_delTria->number_of_vertices();
  std::size_t const l_nFaces = m_delTria->number_of_faces();
  std::size_t const l_nEdges = 0;                          //Assumption

  io_os << l_nVerts << " " << l_nFaces << " " << l_nEdges << std::endl;

  // write the vertices
  std::map<Vertex_handle, std::size_t> l_idVertMap;

  std::size_t l_vertID = 0;
  for( Vertex_iterator l_vit = m_delTria->finite_vertices_begin();
       l_vit != m_delTria->finite_vertices_end();
       ++l_vit, ++l_vertID )
  {
      io_os << *l_vit << std::endl;
      l_idVertMap[ l_vit ] = l_vertID;
  }
  EDGE_CHECK_EQ( l_vertID, l_nVerts );

  // write the vertex indices of each full_cell
  std::size_t l_faceID = 0;
  for( Face_iterator l_fit = m_delTria->finite_faces_begin();
       l_fit != m_delTria->finite_faces_end();
       ++l_fit, ++l_faceID )
  {
      io_os << 3;
      for( int j = 0; j < 3; ++j )
      {
        io_os << ' ' << l_idVertMap[ l_fit->vertex(j) ];
      }
      io_os << std::endl;
  }
  EDGE_CHECK_EQ( l_faceID, l_nFaces );

  return io_os;
}
