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

template< class Polyhedron >
edge_cut::surf::BdryTrimmer<Polyhedron>::BdryTrimmer( Polyhedron& i_bdry, Polyhedron& i_topo ) : m_bdry(i_bdry), m_topo(i_topo), m_isValid(true)
{
  // TODO add sanity checks on polyhedra: non empty, non closed, has border edges, etc
  m_topo.normalize_border();
  m_bdry.normalize_border();
}

template< class Polyhedron >
typename edge_cut::surf::BdryTrimmer<Polyhedron>::Vertex_const
edge_cut::surf::BdryTrimmer<Polyhedron>::getPrevTopoVert( Vertex_const i_vertTopo )
{
  auto l_edge = i_vertTopo->vertex_begin();
  while( !l_edge->is_border_edge() ){ l_edge++; }

  return l_edge->prev()->vertex();
}

template< class Polyhedron >
typename edge_cut::surf::BdryTrimmer<Polyhedron>::Vertex_const
edge_cut::surf::BdryTrimmer<Polyhedron>::getNextTopoVert( Vertex_const i_vertTopo )
{
  auto l_edge = i_vertTopo->vertex_begin();
  while( ! l_edge->next()->is_border_edge() ){ l_edge++; }

  return l_edge->next()->vertex();
}

template< class Polyhedron >
typename edge_cut::surf::BdryTrimmer<Polyhedron>::Vertex
edge_cut::surf::BdryTrimmer<Polyhedron>::getBdryVertex( Vertex_const i_v )
{
  VertexIt l_vIt = m_bdry.vertices_begin();
  while( l_vIt->point() != i_v->point() ) {
    l_vIt++;
    if ( l_vIt == m_bdry.vertices_end() ) {
      EDGE_LOG_INFO << "Could not find coincident vertex between topography and border mesh - mesh will not be trimmed";
      EDGE_LOG_INFO << "Missing point: " << i_v->point();

      m_isValid = false;
      return l_vIt;
    }
  }
  return l_vIt;
}

template< class Polyhedron >
typename edge_cut::surf::BdryTrimmer<Polyhedron>::Halfedge
edge_cut::surf::BdryTrimmer<Polyhedron>::getBdryHalfedge( Vertex_const i_v1, Vertex_const i_v2 )
{
  Vertex l_vert1 = getBdryVertex( i_v1 );
  Vertex l_vert2 = getBdryVertex( i_v2 );

  auto l_edge = l_vert2->vertex_begin();
  while( l_edge->prev()->vertex()->point() != l_vert1->point() ) { l_edge++; }

  return l_edge;
}

template< class Polyhedron >
bool
edge_cut::surf::BdryTrimmer<Polyhedron>::isPosOriented( Halfedge l_halfedge )
{

  if ( l_halfedge->is_border() ) return true;

  Point v1 = l_halfedge->prev()->vertex()->point();
  Point v2 = l_halfedge->vertex()->point();
  Point v3 = l_halfedge->next()->vertex()->point();

  //                   (v2) .-------------------
  //                      ./
  // [above topo]       . /
  //                  .  /
  //                .   /
  //         (v3) .    / U      [below topo]
  //             \    /
  //            W \  /
  // ______________\/   (v1)
  //
  // Consider vector u = v_2 - v_1   and   w = v_3 - v_1
  // We will compute the perpendicular vector from u into v_3 : if the z
  // component of this vector is positive, then we are "above" the topography.
  // Otherwise, we are below, and have to flip the halfedge to the other side.
  double l_dotProdUW = (v3.x()-v1.x()) * (v2.x()-v1.x()) + (v3.y()-v1.y()) * (v2.y()-v1.y()) + (v3.z()-v1.z()) * (v2.z()-v1.z());
  double l_dotProdUU = (v2.x()-v1.x()) * (v2.x()-v1.x()) + (v2.y()-v1.y()) * (v2.y()-v1.y()) + (v2.z()-v1.z()) * (v2.z()-v1.z());

  // This is the z component of: w - (u*w)/(u*u)u
  double l_z = v3.z()-v1.z() - (l_dotProdUW/l_dotProdUU)*(v2.z()-v1.z());

  if ( l_z > 0 ) return true;
  else return false;
}

template< class Polyhedron >
bool
edge_cut::surf::BdryTrimmer<Polyhedron>::trim()
{
  Vertex_const l_topoVert = m_topo.border_halfedges_begin()->vertex();
  Vertex_const l_topoBaseVert = getPrevTopoVert( l_topoVert );
  Vertex_const l_topoVertNext;

  Halfedge l_rootHalfedge, l_h, l_hNext;
  l_rootHalfedge = getBdryHalfedge( l_topoBaseVert, l_topoVert );

  bool l_posOriented = isPosOriented( l_rootHalfedge );

  if ( l_posOriented )  l_h = l_rootHalfedge;
  else                  l_h = l_rootHalfedge->opposite();

  Vertex_const l_startVert = l_topoVert;
  do {
    l_topoVertNext = getNextTopoVert( l_topoVert );

    while ( l_h->next()->vertex()->point() != l_topoVertNext->point() ) {
      if ( l_h->is_border() ) break;

      if ( l_posOriented )  l_hNext = l_h->next()->opposite();
      else                  l_hNext = l_h->prev()->opposite();

      m_bdry.erase_facet( l_h );

      l_h = l_hNext;
    }

    if ( l_posOriented )  l_h = l_h->next();
    else                  l_h = l_h->prev();

    l_topoVert = l_topoVertNext;
  } while ( l_topoVert->point() != l_startVert->point() );


  // Iterate over all non-border halfedges for the one with the highest elevation
  // This halfedge will be on the "upper" connected component
  // TODO WARNING this assumes that the trim splits the mesh into two components
  // TODO add PMP check that there is more than one CC
  m_bdry.normalize_border();
  Halfedge l_hMax = m_bdry.halfedges_begin();
  for ( auto l_it = m_bdry.halfedges_begin(); l_it != m_bdry.border_halfedges_begin(); l_it++ ) {
    if ( l_it->vertex()->point().z() > l_hMax->vertex()->point().z() )
      l_hMax = l_it;
  }

  m_bdry.erase_connected_component( l_hMax );
  m_bdry.normalize_border();

  return m_isValid;
}
