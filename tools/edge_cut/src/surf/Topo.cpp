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

void edge_cut::surf::Topo::computeDelaunay( std::string const & i_topoFile ) {
  // parse the points
  CGAL::Projection_traits_xy_3<
    CGAL::Cartesian<double>
  >::Point l_pt;

  std::ifstream l_ptFile( i_topoFile, std::ios::in );

  while ( l_ptFile >> l_pt ) {
    m_delTria.insert( l_pt );
  }

  EDGE_LOG_INFO << "  computed delaunay triangulation:";
  EDGE_LOG_INFO << "    #vertices: " << m_delTria.number_of_faces();
  EDGE_LOG_INFO << "    #faces:    " << m_delTria.number_of_vertices();
}

edge_cut::surf::Topo::Topo( std::string const & i_topoFile ) {
  computeDelaunay( i_topoFile );
}

bool edge_cut::surf::Topo::interRay( CGAL::Point_3< CGAL::Cartesian<double> > const & i_pt,
                                     bool                                             i_positive ) const {
  CGAL::Delaunay_triangulation_2 <
    CGAL::Projection_traits_xy_3 <
      CGAL::Cartesian<double>
    >
  >::Face_handle l_faceHa;

  // get the possible 2D face 
  l_faceHa = m_delTria.locate( i_pt );

  // return if no face qualifies
  if( l_faceHa == NULL ) return false;
  // do the 3D intersection otherwise and check the side w.r.t. to the face
  else {
      // set up the face
      CGAL::Triangle_3< CGAL::Cartesian<double> > l_face(
        l_faceHa->vertex(0)->point(),
        l_faceHa->vertex(1)->point(),
        l_faceHa->vertex(2)->point()
      );

    // set up the vertical ray
    CGAL::Direction_3< CGAL::Cartesian<double> > l_dir(0,0,1);
    if( i_positive == false ) l_dir = -l_dir;
    CGAL::Ray_3< CGAL::Cartesian<double> > l_ray( i_pt, l_dir );

    // compute the intersection (if available)
    auto l_inter = CGAL::intersection( l_face, l_ray );
    if( l_inter ) return true;
  }

  return false;
}

unsigned short edge_cut::surf::Topo::interSeg( CGAL::Point_3< CGAL::Cartesian<double> > const & i_segPt1,
                                               CGAL::Point_3< CGAL::Cartesian<double> > const & i_segPt2,
                                               CGAL::Point_3< CGAL::Cartesian<double> >         o_inters [C_MAX_SURF_INTER] ) const {
  // set up segment
  CGAL::Segment_3< CGAL::Cartesian<double> > l_seg( i_segPt1, i_segPt2 );

  // get intersections in 2D
  CGAL::Delaunay_triangulation_2 <
    CGAL::Projection_traits_xy_3 <
      CGAL::Cartesian<double>
    >
  >::Line_face_circulator l_lineWalk = m_delTria.line_walk( i_segPt1, i_segPt2 ), l_lineWalkDone(l_lineWalk);

  unsigned int l_interCount = 0;

  if( l_lineWalk != 0) {
    do {
      // set up the face
      CGAL::Triangle_3< CGAL::Cartesian<double> > l_face(
        (*l_lineWalk).vertex(0)->point(),
        (*l_lineWalk).vertex(1)->point(),
        (*l_lineWalk).vertex(2)->point()
      );

      // compute the intersection (if available)
      auto l_inter = CGAL::intersection( l_face, l_seg );

      // only continue for non-empty intersections
      if( l_inter ) {
        // save intersection (dependent on return type of the intersection)
        if( auto l_pt = boost::get< CGAL::Point_3< CGAL::Cartesian<double> > >( &*l_inter ) ) {
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
