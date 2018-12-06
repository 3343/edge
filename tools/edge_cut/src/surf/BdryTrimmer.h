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
class edge_cut::surf::BdryTrimmer {
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

  /**
   * Constructor: Checks validity of polyhedron input
   *
   * @param io_bdry polyhedron for boundary mesh
   * @param i_topo poolyhedron for topography mesh
   **/
  BdryTrimmer( Polyhedron & io_bdry, Polyhedron const & i_topo );

  /**
  * Trims the boundary mesh of all facets which lie above the topography mesh
  *
  * @return None
  **/
  void trim();

private:
  /**
   * Gets the previous vertex along the boundary of topography mesh
   *
   * @param i_vertTopo vertex on topography mesh
   * @return const handle to previous vertex
   **/
  Vertex_const getPrevTopoVert( Vertex_const i_vertTopo ) const;

  /**
   * Gets the next vertex along the boundary of topography mesh
   *
   * @param i_vertTopo vertex on topography mesh
   * @return const handle to next vertex
   **/
  Vertex_const getNextTopoVert( Vertex_const i_vertTopo ) const;

  /**
   * Gets the vertex on the boundary mesh which coincides with the given vertex
   * on the topography mesh
   *
   * @param i_v vertex on the topography mesh to match
   * @return handle to the matching vertex on the boundary mesh
   **/
  Vertex getBdryVertex( Vertex_const i_v ) const;

  /**
   * Gets the halfedge on the boundary mesh with endpoints coinciding with
   * the two given vertices on the topography mesh
   *
   * @param i_v1 topography vertex which is the base point of the desired halfedge
   * @param i_v2 topography vertex which is the incident point of the desired halfedge
   * @param handle to the matching halfedge on the boundary
   **/
  Halfedge getBdryHalfedge( Vertex_const i_v1, Vertex_const i_v2 ) const;

  /**
   * Checks the relative orientations of the topography mesh and the boundary mesh.
   * Given a directed edge on the boundary mesh which also represents an edge on the
   * topography mesh, determines whether that halfedge belongs to the face above
   * the topography or below.
   *
   * @param i_halfedge halfedge on the boundary to test
   * @return True if the halfedge belongs to a face which extends above
   * the topography mesh
   **/
  bool isPosOriented( Halfedge i_halfedge ) const;

  /**
   * Deletes all connected components of m_bdry except for the one containing
   * the face i_rep
   *
   * @param i_rep facet in the boundary mesh which is contained in the connected
   *              component to save
   * @return None
   **/
  void deleteConnCompsExcept( Face i_rep );

public:
  // Boundary mesh
  Polyhedron       & m_bdry;

private:
  // Topography mesh
  Polyhedron const & m_topo;
};

#include "BdryTrimmer.inc"

#endif
