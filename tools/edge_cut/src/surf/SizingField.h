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
 * Class describing how mesh refinement varies throughout the domain
 **/
#ifndef EDGE_CUT_SIZE_FIELD_H

#include <map>
#include "io/logging.hpp"

namespace edge_cut {
  namespace surf {
    template< class Mesh_domain > struct SizingField;
  }
}

// Mesh_domain must be model of CGAL concept MeshDomain_3
template< class Mesh_domain >
struct edge_cut::surf::SizingField {
  typedef typename Mesh_domain::FT       FT;
  typedef typename Mesh_domain::Point_3  Point_3;
  typedef typename Mesh_domain::Index    Index;

  /**
   * Constructor
   *
   * @param i_innerVal  Base value to be scaled
   * @param i_scale     Factor by which value is increased outside of outer radius
   * @param i_center    Center of the circular refinement region
   * @param i_innerRad  Inner radius of the circular refinement region
   * @param i_outerRad  Outer radius of the circular refinement region
   * @param i_layers    Set of (depth, scale) pairs describing how value is scaled at different depths
   **/
  SizingField( FT i_innerVal, FT i_scale, Point_3 i_center, FT i_innerRad, FT i_outerRad, std::map< FT, FT, std::greater<FT> > i_layers );

  /**
   * Calculates the scaled value for the point p
   *
   * @param p Point at which to calculate the associated value
   * @return The scaled value
   **/
  FT operator()( const Point_3& p, const int, const Index& ) const;

  /**
   * Gets the scaling factor corresponding to a particular depth
   *
   * @param i_depth Depth level
   * @return The scaling factor for depth i_z
   **/
  FT getDepthScale( FT i_depth ) const;

  // Base value for the sizing field that will be scaled, depending on the point
  FT m_innerVal;

  // Scale factor for circular-region refinement; m_innerVal will scale by this factor outside of the outer radius
  FT m_scale;

  // Center point of the circular-region refinement
  Point_3 m_center;

  // Inner radius for the circular-region refinement
  FT m_innerRad;

  // Outer radius for the circular-region refinement
  FT m_outerRad;

  // Reciprocal of the difference between the inner and outer radii
  FT m_widthInv;

  // (Depth, Scale) pairs, such that points with z-coord above Depth will be
  // scaled by a factor of Scale. In cases with multiple layers, the Scale associated
  // to the highest Depth takes precedence.
  std::map< FT, FT, std::greater<FT> >  m_layers;
};

#include "SizingField.inc"

#endif
