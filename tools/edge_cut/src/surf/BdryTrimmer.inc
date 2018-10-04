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
edge_cut::surf::BdryTrimmer<Polyhedron>::BdryTrimmer( Polyhedron& i_bdry, Polyhedron& i_topo ) :
  m_bdry(i_bdry),
  m_topo(i_topo)
{
  m_bdry.normalize_border();

  assert( !m_topo.is_empty() );
  assert( !m_topo.is_closed() );
  assert( m_topo.is_pure_triangle() );
  assert( m_topo.normalized_border_is_valid() );

  assert( !m_bdry.is_empty() );
  assert( !m_bdry.is_closed() );
  assert( m_bdry.is_pure_triangle() );
  assert( m_bdry.normalized_border_is_valid() );
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
  for ( ;
        l_vIt != m_bdry.vertices_end();
        l_vIt++ )
  {
    if ( l_vIt->point() == i_v->point() ) break;
  }

  if ( l_vIt == m_bdry.vertices_end() )
  {
    EDGE_LOG_INFO << "Could not find coincident vertex between topography and border mesh - mesh will not be trimmed";
    EDGE_LOG_INFO << "Missing point: " << i_v->point();
    assert( l_vIt != m_bdry.vertices_end() );
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

// Checks the relative orientation of the two meshes. Given a halfedge on the
// boundary mesh which also represents an edge on the topography mesh,
// determines whether that halfedge belongs to the face above the topography or
// below.
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
void
edge_cut::surf::BdryTrimmer<Polyhedron>::deleteConnCompsExcept( Face i_rep )
{
  typedef typename boost::graph_traits<Polyhedron>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<Polyhedron>::vertex_descriptor vert_descriptor;

  std::map<face_descriptor, std::size_t> fi_map;
  std::map<vert_descriptor, std::size_t> vi_map;
  std::size_t fid = 0;
  std::size_t vid = 0;
  BOOST_FOREACH(face_descriptor f, faces(m_bdry))
  {
    fi_map[f]=fid++;
  }
  BOOST_FOREACH(vert_descriptor v, vertices(m_bdry))
  {
    vi_map[v]=vid++;
  }

  std::vector< Face > l_ccsToKeep( 1, i_rep );
  CGAL::Polygon_mesh_processing::keep_connected_components( m_bdry,
                                                            l_ccsToKeep,
                                                            CGAL::parameters::
                                                              face_index_map( boost::make_assoc_property_map( fi_map )).
                                                              vertex_index_map( boost::make_assoc_property_map( vi_map )) );

  m_bdry.normalize_border();

  return;
}


template< class Polyhedron >
void
edge_cut::surf::BdryTrimmer<Polyhedron>::trim()
{
  Vertex_const l_topoVert = m_topo.border_halfedges_begin()->vertex();
  Vertex_const l_topoBaseVert = getPrevTopoVert( l_topoVert );
  Vertex_const l_topoVertNext;

  Halfedge l_rootHalfedge, l_h, l_hNext;
  l_rootHalfedge = getBdryHalfedge( l_topoBaseVert, l_topoVert );

  // We require that the two meshes have a fixed (arbitrarily chosen) relative
  // orientation with each other. In particular, we require that if the border
  // of the topography mesh is traversed in order (as defined by getNextTopoVert),
  // then the corresponding halfedges on the boundary mesh which span subsequent
  // vertices will lie above the topography. See isPosOriented() for more info.
  if( !isPosOriented( l_rootHalfedge ) ) {
    m_bdry.inside_out();
    l_rootHalfedge = getBdryHalfedge( l_topoBaseVert, l_topoVert );
    assert( isPosOriented( l_rootHalfedge ) );
  }

  // Get a representative face for the region below the topography
  // We will use this when deleting connected components later.
  Face l_ccRep = l_rootHalfedge->opposite()->face();

  l_h = l_rootHalfedge;

  bool l_vertDone = false;
  Vertex_const l_startVert = l_topoVert;
  do { // For each vertex on the border of topography
    l_topoVertNext = getNextTopoVert( l_topoVert );

    l_vertDone = false;
    while ( !l_vertDone ) { // Remove each facet adjacent to vert which is above topo

      // Get halfedge handle to next adjacent face above topography, if there
      // is one. Otherwise, get halfedge incident to topoVertNext
      if ( l_h->next()->vertex()->point() == l_topoVertNext->point() ){
        l_hNext = l_h->next();
        l_vertDone = true;
      }
      else if ( !l_h->next_on_vertex()->is_border() ){
        l_hNext = l_h->next_on_vertex();
      }
      else {
        if ( l_h->next_on_vertex()->next()->vertex()->point() == l_topoVertNext->point() ) {
          l_hNext = l_h->next_on_vertex()->next();
          l_vertDone = true;
        }
        else if ( !l_h->next_on_vertex()->next_on_vertex()->is_border() ) {
          l_hNext = l_h->next_on_vertex()->next_on_vertex();
        }
        else {
          // This case should not be possible in my understanding
          assert( false );
        }
      }

      // Delete current facet above topography, if there is one
      if ( !l_h->is_border() ) {
        m_bdry.erase_facet( l_h );
      }
      l_h = l_hNext;
    }

    l_topoVert = l_topoVertNext;
  } while ( l_topoVert->point() != l_startVert->point() );

  deleteConnCompsExcept( l_ccRep );

  return;
}
