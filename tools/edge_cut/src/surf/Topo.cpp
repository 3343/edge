/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2017, Regents of the University of California
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


edge_cut::surf::Topo::Topo( std::string const & i_topoFile, double i_xMin, double i_xMax, double i_yMin, double i_yMax, double i_zMin, double i_zMax ) :
  m_xMin( i_xMin ),
  m_xMax( i_xMax ),
  m_yMin( i_yMin ),
  m_yMax( i_yMax ),
  m_zMin( i_zMin ),
  m_zMax( i_zMax )
{
  std::ifstream l_ptFile( i_topoFile, std::ios::in );
  std::istream_iterator< TopoPoint > l_fileBegin( l_ptFile );
  std::istream_iterator< TopoPoint > l_fileEnd;

  m_delTria = new Triangulation( l_fileBegin, l_fileEnd );

  EDGE_LOG_INFO << "  computed delaunay triangulation:";
  EDGE_LOG_INFO << "    #vertices: " << m_delTria->number_of_vertices();
  EDGE_LOG_INFO << "    #faces:    " << m_delTria->number_of_faces();
}

edge_cut::surf::Topo::~Topo()
{
  delete m_delTria;
}

double edge_cut::surf::Topo::topoDisp( TopoPoint const & i_pt ) const {
  double l_zTopoDisp = 1;

  Triangulation::Face_handle l_faceHa = m_delTria->locate( i_pt );

  // return if no face qualifies
  if( l_faceHa == NULL ) {
    std::cout << "Error: Got NULL value when locating point on topography"
              << std::endl;
    return 0;
  }
  // do the 3D intersection otherwise and check the side w.r.t. to the face
  else {
    // set up the face
    CGAL::Triangle_3< K > l_face(
      l_faceHa->vertex(0)->point(),
      l_faceHa->vertex(1)->point(),
      l_faceHa->vertex(2)->point()
    );

    // set up the vertical ray
    CGAL::Direction_3< K > l_dirPos(0,0,1);
    CGAL::Direction_3< K > l_dirNeg(0,0,-1);
    CGAL::Ray_3< K > l_rayPos( i_pt, l_dirPos );
    CGAL::Ray_3< K > l_rayNeg( i_pt, l_dirNeg );

    // compute the intersection (if available)
    auto l_interPos = CGAL::intersection( l_face, l_rayPos );
    auto l_interNeg = CGAL::intersection( l_face, l_rayNeg );
    if ( l_interPos )
      l_zTopoDisp = i_pt.z() - boost::get<CGAL::Point_3< K > >(&*l_interPos)->z();
    else if ( l_interNeg )
      l_zTopoDisp = i_pt.z() - boost::get<CGAL::Point_3< K > >(&*l_interNeg)->z();
  }

  return l_zTopoDisp;
}


edge_cut::surf::TopoPoint edge_cut::surf::Topo::interpolatePt( double i_x, double i_y ) const {
  TopoPoint l_basePt( i_x, i_y, 0 );
  double l_zDist = topoDisp( l_basePt );

  return TopoPoint( i_x, i_y, -1*l_zDist );
}


bool edge_cut::surf::Topo::interRay( TopoPoint const & i_pt,
                                     bool  i_positive ) const {
  Triangulation::Face_handle l_faceHa;

  // get the possible 2D face
  l_faceHa = m_delTria->locate( i_pt );

  // return if no face qualifies
  if( l_faceHa == NULL ) return false;
  // do the 3D intersection otherwise and check the side w.r.t. to the face
  else {
      // set up the face
      CGAL::Triangle_3< K > l_face(
        l_faceHa->vertex(0)->point(),
        l_faceHa->vertex(1)->point(),
        l_faceHa->vertex(2)->point()
      );

    // set up the vertical ray
    CGAL::Direction_3< K > l_dir(0,0,1);
    if( i_positive == false ) l_dir = -l_dir;
    CGAL::Ray_3< K > l_ray( i_pt, l_dir );

    // compute the intersection (if available)
    auto l_inter = CGAL::intersection( l_face, l_ray );
    if( l_inter ) return true;
  }

  return false;
}

unsigned short edge_cut::surf::Topo::interSeg( TopoPoint const & i_segPt1,
                                               TopoPoint const & i_segPt2,
                                               TopoPoint         o_inters [C_MAX_SURF_INTER] ) const {
  // set up segment
  CGAL::Segment_3< K > l_seg( i_segPt1, i_segPt2 );

  // get intersections in 2D
  Triangulation::Line_face_circulator l_lineWalk = m_delTria->line_walk( i_segPt1, i_segPt2 ), l_lineWalkDone(l_lineWalk);

  unsigned int l_interCount = 0;

  if( l_lineWalk != 0) {
    do {
      // set up the face
      CGAL::Triangle_3< K > l_face(
        (*l_lineWalk).vertex(0)->point(),
        (*l_lineWalk).vertex(1)->point(),
        (*l_lineWalk).vertex(2)->point()
      );

      // compute the intersection (if available)
      auto l_inter = CGAL::intersection( l_face, l_seg );

      // only continue for non-empty intersections
      if( l_inter ) {
        // save intersection (dependent on return type of the intersection)
        if( auto l_pt = boost::get< CGAL::Point_3< K > >( &*l_inter ) ) {
          o_inters[l_interCount] = *l_pt;
        }
        else EDGE_LOG_FATAL << "intersection return type not supported";

        // increase counter
        l_interCount++;
      }
    }
    while(++l_lineWalk != l_lineWalkDone);
  }

  // return number of intersections
  return l_interCount;
}

std::ostream & edge_cut::surf::Topo::writeTriaToOff( std::ostream & os ) const {
  typedef typename Triangulation::Vertex_handle                Vertex_handle;
  typedef typename Triangulation::Finite_vertices_iterator     Vertex_iterator;
  typedef typename Triangulation::Finite_faces_iterator        Face_iterator;

  os << "OFF" << std::endl;

  // outputs the number of vertices and faces
  std::size_t num_verts = m_delTria->number_of_vertices();
  std::size_t num_faces = m_delTria->number_of_faces();
  std::size_t num_edges = 0;                          //Assumption

  os << num_verts << " " << num_faces << " " << num_edges << std::endl;

  // write the vertices
  std::map<Vertex_handle, std::size_t> index_of_vertex;

  std::size_t v_idx = 0;
  for( Vertex_iterator it = m_delTria->finite_vertices_begin(); it != m_delTria->finite_vertices_end(); ++it, ++v_idx )
  {
      os << *it << std::endl;
      index_of_vertex[it] = v_idx;
  }
  CGAL_assertion( v_idx == num_verts );

  // write the vertex indices of each full_cell
  std::size_t f_idx = 0;
  for( Face_iterator it = m_delTria->finite_faces_begin(); it != m_delTria->finite_faces_end(); ++it, ++f_idx )
  {
      os << 3;
      for( int j = 0; j < 3; ++j )
      {
        os << ' ' << index_of_vertex[ it->vertex(j) ];
      }
      os << std::endl;
  }
  CGAL_assertion( f_idx == num_faces );

  return os;
}
