/**
 * @file This file is part of EDGE.
 *
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
 * Delete extraneous bits of mesh which extend above topographical surface
 **/

#ifndef EDGE_CUT_BDRY_TRIMMER_H_
#define EDGE_CUT_BDRY_TRIMMER_H_

#include "io/logging.hpp"

namespace edge_cut {
  namespace surf {
    template< class Polyhedron > class BdryTrimmer;
  }
}

template< class Polyhedron >
class edge_cut::surf::BdryTrimmer{
public:
  // Const types are used for the (unchanged) topo mesh
  typedef typename Polyhedron::Vertex_const_handle       Vertex_const;

  // Non-const types are for the boundary mesh
  typedef typename Polyhedron::Vertex_handle             Vertex;
  typedef typename Polyhedron::Vertex_iterator           VertexIt;
  typedef typename Polyhedron::Halfedge_handle           Halfedge;
  typedef typename Polyhedron::Face_handle               Face;

  // We cannot compare vertices and halfedges on two different polyhedral meshes,
  // because they will always have different connectivity information.
  // The Point type has no connectivity data, so it can be used to compare
  // vertices on two different meshes.
  typedef typename Polyhedron::Point                     Point;

  BdryTrimmer( Polyhedron& i_bdry, Polyhedron& i_topo );

  void trim();

private:
  Vertex_const getPrevTopoVert( Vertex_const i_vertTopo );

  Vertex_const getNextTopoVert( Vertex_const i_vertTopo );

  Vertex getBdryVertex( Vertex_const i_v );

  Halfedge getBdryHalfedge( Vertex_const i_v1, Vertex_const i_v2 );

  bool isPosOriented( Halfedge l_halfedge );

  void deleteConnCompsExcept( Face i_rep );

public:
  Polyhedron& m_bdry;

private:
  const Polyhedron& m_topo;
};

#include "BdryTrimmer.inc"

#endif
