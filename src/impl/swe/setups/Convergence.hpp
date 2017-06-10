/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2016-2017, Regents of the University of California
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
 * Setup for convergence studies of the shallow water equations.
 **/
#ifndef CONVERGENCE_HPP
#define CONVERGENCE_HPP

#include "constants.hpp"
#include <cassert>
#include "mesh/common.hpp"
#include "dg/QuadraturePoints.h"

namespace edge {
  namespace swe {
    namespace setups {
      class Convergence;
    }
  }
}

class edge::swe::setups::Convergence {
  public:
#if defined PP_T_ELEMENTS_LINE
    /**
     * Sets the specified dam break problem.
     *
     * @param i_cfr concurrent forward run.
     * @param i_nElements number of elements.
     * @param i_connElVe connectivity information from elements to vertices.
     * @param i_vertexChars vertex characteristics.
     * @param i_elementChars element characteristics.
     * @param i_discL location of the left discontinuity.
     * @param i_discR location of the right discontinuity.
     * @param i_hL left water height.
     * @param i_hR right water height.
     * @param o_b will be set to bathymetry.
     * @param o_dofs will be set to dofs.
     **/
    static void setDamBreak1D(       int_cfr          i_cfr,
                                     int_el           i_nElements,
                               const int_el         (*i_connElVe)[C_ENT[T_SDISC.ELEMENT].N_VERTICES],
                               const t_vertexChars   *i_vertexChars,
                               const t_elementChars  *i_elementChars,
                                     real_mesh       i_discL,
                                     real_mesh       i_discR,
                                     real_base       i_hL,
                                     real_base       i_hR,
                                     real_base     (*o_b)[1][1],
                                     real_base     (*o_dofs)[N_QUANTITIES][1][N_CRUNS] ) {
      assert( i_discL < i_discR );

      // iterate over elements
      for( int_el l_el = 0; l_el < i_nElements; l_el++ ) {
        for( int_md l_md = 0; l_md < N_ELEMENT_MODES; l_md++ ) {
          // set momentum and bathymetry straight to zero
          o_dofs[l_el][1][l_md][i_cfr] = 0;
          o_b[l_el][0][l_md] = 0;

          // get the vertex coords
          real_mesh l_veCoords[3][C_ENT[T_SDISC.ELEMENT].N_VERTICES];
          mesh::common< T_SDISC.ELEMENT >::getElVeCoords( l_el, i_connElVe, i_vertexChars, l_veCoords );

          // get gaussian quadrature points and their weigths
          std::vector< real_mesh > l_ptsX, l_ptsY, l_ptsZ, l_weights;

          dg::QuadraturePoints::getQpts( T_SDISC.ELEMENT,
                                         ORDER+1,
                                         l_veCoords,
                                         l_ptsX, l_ptsY, l_ptsZ, l_weights );

          // iterate over gauss points
          for( int_md l_gp = 0; l_gp < l_ptsX.size(); l_gp++ ) {
            if( l_ptsX[l_gp] > i_discL &&
                l_ptsX[l_gp] < i_discR ) o_dofs[l_el][0][l_md][i_cfr] = i_hL * l_weights[l_gp];
            else                         o_dofs[l_el][0][l_md][i_cfr] = i_hR * l_weights[l_gp];
          }
        }
      }
    }
#else
#error element type not supported.
#endif
};

#endif
