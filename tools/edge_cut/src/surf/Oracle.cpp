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
 * Oracle for the surface mesh.
 **/

#include "Oracle.h"
#include "io/logging.hpp"

edge_cut::surf::Oracle::Oracle( double              i_box[5],
                                std::string const & i_topoFile ): m_topo( i_topoFile ) {
  // copy box to member
  for( unsigned short l_bd = 0; l_bd < 5; l_bd++ ) m_box[l_bd] = i_box[l_bd];
}

CGAL::Surface_mesh_default_triangulation_3::Geom_traits::FT edge_cut::surf::Oracle::operator()(
  CGAL::Surface_mesh_default_triangulation_3::Geom_traits::Point_3 i_pt
) const {
  bool l_inside = true;

  // left and right boundary
  l_inside = l_inside && ( i_pt.x() > m_box[0] && i_pt.x() < m_box[1] );
  // front and back
  l_inside = l_inside && ( i_pt.y() > m_box[2] && i_pt.y() < m_box[3] );
  // bottom
  l_inside = l_inside && ( i_pt.z() > m_box[4] );

  // return if already outside
  if( !l_inside ) {
    return CGAL::Surface_mesh_default_triangulation_3::Geom_traits::FT( 0 );
  }
  else {
    // check for intersection
    CGAL::Point_3< K > l_pt( CGAL::to_double( i_pt.x() ),
                                                   CGAL::to_double( i_pt.y() ),
                                                   CGAL::to_double( i_pt.z() ) );
    if( !m_topo.interRay( l_pt ) ) l_inside = false;
  }

  if( !l_inside )
    return CGAL::Surface_mesh_default_triangulation_3::Geom_traits::FT( 0 );
  else
    return CGAL::Surface_mesh_default_triangulation_3::Geom_traits::FT( 1 );
}
