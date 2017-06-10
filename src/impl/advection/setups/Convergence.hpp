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
 * Setup for convergence studies of the advection equation.
 **/

#ifndef CONVERGENCE_HPP
#define CONVERGENCE_HPP

#include "constants.hpp"
#include "mesh/common.hpp"
#include "dg/QuadraturePoints.h"
#include <cmath>
#include <cassert>

namespace edge {
  namespace advection {
    namespace setups {
      class Convergence;
    }
  }
}

class edge::advection::setups::Convergence {
  private:

  public:
    /**
     * Gets the initial value for a given point.
     *
     * @param i_pos location of the point.
     * @param i_off offset for the initial data.
     * @param i_a start of the considered interval.
     * @param i_b end of the considered interval.
     * @param i_k number of waves ( 2*PI/k ).
     **/
    static real_base getSineWavesVal1D( real_mesh    i_pos,
                                        real_mesh    i_off,
                                        real_mesh    i_a,
                                        real_mesh    i_b,
                                        unsigned int i_k ) {
      assert( i_b - i_a > TOL.MESH );

      // compute repeat interval length
      real_mesh l_length = ( i_b - i_a ) / i_k;

      // offset gauss point
      real_mesh l_pos = i_pos - i_off;

      // handle possible invalid values due to offset
      if( l_pos < i_a ) l_pos = i_b - (  i_a - l_pos );

      // shift point to first interval
      while( l_pos > i_a ) l_pos -= l_length;
      l_pos += l_length;

      // shift to [0, ..
      l_pos -= i_a;

      // map to [0,2pi]
      l_pos /= l_length;

      l_pos *= 2 * M_PI;

      // return value
      return std::sin( l_pos );
    }

    /**
     * Sets a sine curve with variable frequency.
     *
     * @param i_cfr concurrent forward run.
     * @param i_basis DG-basis.
     * @param i_nElements number of elements.
     * @param i_connElVe connectivity information from elements to vertices.
     * @param i_vertexChars vertex characteristics.
     * @param i_elementChars element characteristics.
     * @param o_dofs will be set to initial DOFs.
     * @param i_aX start of the considered domain in x-direction.
     * @param i_bX end of the considered domain in x-drection.
     * @param i_aY start of the considered domain in y-direction.
     * @param i_bY end of the considered domain in y-direction.
     * @param i_aZ start of the considered domain in z-direction.
     * @param i_bZ end of the considered domain in z-direction.
     * @param i_offsetX offset used for the initial data in x-direction.
     * @param i_offsetY offset used for the initial data in y-direction
     * @param i_offsetZ offset used for the initial data in z-direction.
     * @param i_k number of waves ( 2*PI/k )
     **/
    static void setSineWaves(       int_cfr          i_cfr,
                              const dg::Basis       &i_basis,
                                    int_el           i_nElements,
                              const int_el         (*i_connElVe)[C_ENT[T_SDISC.ELEMENT].N_VERTICES],
                              const t_vertexChars   *i_vertexChars,
                              const t_elementChars  *i_elementChars,
                                    real_base      (*o_dofs)[1][N_ELEMENT_MODES][N_CRUNS],
                                    real_mesh        i_aX = 0,
                                    real_mesh        i_bX = 1,
                                    real_mesh        i_aY = 0,
                                    real_mesh        i_bY = 1,
                                    real_mesh        i_aZ = 0,
                                    real_mesh        i_bZ = 1,
                                    real_mesh        i_offsetX = 0,
                                    real_mesh        i_offsetY = 0,
                                    real_mesh        i_offsetZ = 0,
                                    unsigned int     i_k = 1 ) {
      assert( i_bX > i_aX );
      assert( i_bY > i_aY || N_DIM == 1 );
      assert( i_bZ > i_aZ || N_DIM == 1 || N_DIM == 2 );

      // iterate over elements
      for( int_el l_el = 0; l_el < i_nElements; l_el++ ) {
        // get the vertex coords
        real_mesh l_veCoords[3][C_ENT[T_SDISC.ELEMENT].N_VERTICES];
        mesh::common< T_SDISC.ELEMENT >::getElVeCoords( l_el, i_connElVe, i_vertexChars, l_veCoords );

        // get gaussian quadrature points and their weigths
        std::vector< real_mesh > l_ptsX, l_ptsY, l_ptsZ, l_weights;

        dg::QuadraturePoints::getQpts( T_SDISC.ELEMENT,
                                       ORDER+1,
                                       l_veCoords,
                                       l_ptsX, l_ptsY, l_ptsZ, l_weights );

        // evaluated initial quantities at quad points
        std::vector< real_base > l_q0; l_q0.resize( l_weights.size() );

        // iterate over quadrature points
        for( int_md l_qp = 0; l_qp < l_ptsX.size(); l_qp++ ) {
          // return value
          l_q0[l_qp] = getSineWavesVal1D( l_ptsX[l_qp],
                                          i_offsetX,
                                          i_aX,
                                          i_bX,
                                          i_k );
          if( N_DIM > 1 ) {
            l_q0[l_qp] *= getSineWavesVal1D( l_ptsY[l_qp],
                                             i_offsetY,
                                             i_aY,
                                             i_bY,
                                             i_k );
          }
          if( N_DIM > 2 ) {
            l_q0[l_qp] *= getSineWavesVal1D( l_ptsZ[l_qp],
                                             i_offsetZ,
                                             i_aZ,
                                             i_bZ,
                                             i_k );
          }
        }

        // map to modes
        real_base l_modes[N_ELEMENT_MODES];
        i_basis.qpts2modal( &l_q0[0],
                            ORDER+1,
                            l_modes );

        // scatter to cfr-structure
        for( int_md l_md = 0; l_md < N_ELEMENT_MODES; l_md++ ) {
          o_dofs[l_el][0][l_md][i_cfr] = l_modes[l_md];
        }
      }
    }

    /**
     * Sets constant back ground wave speeds
     *
     * @param i_nElements number of elements.
     * @param o_waveSpeeds will be set to constant wave speeds
     * @param i_cX wave speed in x-direction.
     * @param i_cY wave speed in y-direction (ignored for #dims==1).
     * @param i_cZ wave speed in z-direction (ignored for #dims!=3).
     **/
    static void setConstantSpeed( int_el      i_nElements,
                                  t_bgPars  (*o_waveSpeeds)[1],
                                  real_base   i_cX,
                                  real_base   i_cY = 0,
                                  real_base   i_cZ = 0) {
      for( int_el l_el = 0; l_el < i_nElements; l_el++ ) {
        o_waveSpeeds[l_el][0].a = i_cX;
#if PP_N_DIM > 1
        o_waveSpeeds[l_el][0].b = i_cY;
#endif
#if PP_N_DIM > 2
        o_waveSpeeds[l_el][0].c = i_cZ;
#endif
      }
    }

    /**
     * Gets the error norms for the sine setup.
     *
     * @param i_basis basis.
     * @param i_cfr concurrent forward to get the norms for.
     * @param i_inMap index mapping between computational and mesh data.
     * @param i_elLayout element layout.
     * @param i_connect connectivity information.
     * @param i_elVe elements' adjacent vertices.
     * @param i_vertexChars vertex characteristics.
     * @param i_dofs DOFs.
     * @param o_norms will be set to L1, L2 and Linf-norm.
     * @param i_aX start of the considered domain in x-direction.
     * @param i_bX end of the considered domain in x-drection.
     * @param i_aY start of the considered domain in y-direction.
     * @param i_bY end of the considered domain in y-direction.
     * @param i_aZ start of the considered domain in z-direction.
     * @param i_bZ end of the considered domain in z-direction.
     * @param i_offsetX offset used for the initial data in x-direction.
     * @param i_offsetY offset used for the initial data in y-direction
     * @param i_offsetZ offset used for the initial data in z-direction.
     * @param i_k number of sine curves.
     **/
    static void getSineErrorNorms( const  edge::dg::Basis &i_basis,
                                          int_cfr          i_cfr,
                                   const  t_inMap         *i_inMap,
                                   const  t_enLayout      &i_elLayout,
                                   const  t_connect       &i_connect,
                                          t_vertexChars   *i_vertexChars,
                                          real_base      (*i_dofs)[1][N_ELEMENT_MODES][N_CRUNS],
                                          double           o_norms[3][1],
                                          real_mesh        i_aX = 0,
                                          real_mesh        i_bX = 1,
                                          real_mesh        i_aY = 0,
                                          real_mesh        i_bY = 1,
                                          real_mesh        i_aZ = 0,
                                          real_mesh        i_bZ = 1,
                                          real_mesh        i_offsetX = 0,
                                          real_mesh        i_offsetY = 0,
                                          real_mesh        i_offsetZ = 0,
                                          unsigned int     i_k = 1  ) {
      // reset norms
      o_norms[0][0] = 0; o_norms[1][0] = 0; o_norms[2][0] = 0;

      // iterate over time groups
      for( unsigned short l_tg = 0; l_tg < i_elLayout.timeGroups.size(); l_tg++ ) {
        // iterate over owned elements
        int_el l_lower = i_elLayout.timeGroups[l_tg].inner.first;
        int_el l_upper = l_lower + i_elLayout.timeGroups[l_tg].nEntsOwn;
        for( int_el l_el = l_lower; l_el < l_upper; l_el++ ) {
          // determine first data element associated
          int_el l_uEl = i_inMap->elMeDa[ i_inMap->elDaMe[l_el] ];
          // count element-error only once (copy element might be duplicated)
          if( l_el != l_uEl ) continue;

          // get the vertex coords
          real_mesh l_veCoords[3][C_ENT[T_SDISC.ELEMENT].N_VERTICES];
          mesh::common< T_SDISC.ELEMENT >::getElVeCoords( l_el, i_connect.elVe, i_vertexChars, l_veCoords );

          /**
           * Get exact data.
           **/
          // get gaussian quadrature points and their weigths in physical coords
          std::vector< real_mesh > l_ptsX, l_ptsY, l_ptsZ, l_weights;

          dg::QuadraturePoints::getQpts( T_SDISC.ELEMENT,
                                         ORDER+1,
                                         l_veCoords,
                                         l_ptsX, l_ptsY, l_ptsZ, l_weights );

          // evaluated initial quantities at quad points
          std::vector< real_base > l_q0; l_q0.resize( l_weights.size() );

          // iterate over gauss points
          for( int_md l_qp = 0; l_qp < l_ptsX.size(); l_qp++ ) {
            l_q0[l_qp] = getSineWavesVal1D( l_ptsX[l_qp],
                                            i_offsetX,
                                            i_aX,
                                            i_bX,
                                            i_k );
            if( N_DIM > 1 ) {
              l_q0[l_qp] *= getSineWavesVal1D( l_ptsY[l_qp],
                                               i_offsetY,
                                               i_aY,
                                               i_bY,
                                               i_k );
            }
            if( N_DIM > 2 ) {
              l_q0[l_qp] *= getSineWavesVal1D( l_ptsZ[l_qp],
                                               i_offsetZ,
                                               i_aZ,
                                               i_bZ,
                                               i_k );
            }
          }
          /**
           * Get solution.
           **/
          // get gaussian quadrature points and their weigths for the reference element
          std::vector< real_mesh > l_soPtsX, l_soPtsY, l_soPtsZ, l_soWeights;

          dg::QuadraturePoints::getQpts( T_SDISC.ELEMENT,
                                         ORDER+1,
                                         C_REF_ELEMENT.VE.ENT[T_SDISC.ELEMENT],
                                         l_soPtsX, l_soPtsY, l_soPtsZ, l_soWeights );

          assert( l_soPtsX.size() == l_ptsX.size() );

          // solution at quad points
          std::vector< real_base > l_qS; l_qS.resize( l_soWeights.size() );

          // gather modes
          real_base l_modes[N_ELEMENT_MODES];
          for( int_md l_md = 0; l_md < N_ELEMENT_MODES; l_md++ ) {
            l_modes[l_md] = i_dofs[l_el][0][l_md][i_cfr];
          }

          // iterate over quad points and get the solution
          for( int_md l_qp = 0; l_qp < l_soPtsX.size(); l_qp++ ) {
            l_qS[l_qp] = i_basis.modal2refPtVal( ORDER+1, l_qp, l_modes );
          }

          // add errors
          for( int_md l_qp = 0; l_qp < l_ptsX.size(); l_qp++ ) {
            double l_diff = std::abs( l_q0[l_qp] - l_qS[l_qp] );
            o_norms[0][0] += l_diff *          l_weights[l_qp];
            o_norms[1][0] += l_diff * l_diff * l_weights[l_qp];
            o_norms[2][0] = std::max( o_norms[2][0], l_diff );
          }
        }
      }

      // do the MPI-reductions
#ifdef PP_USE_MPI
      double l_rNorms[3];

      MPI_Allreduce( &o_norms[0][0], &l_rNorms[0], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
      MPI_Allreduce( &o_norms[1][0], &l_rNorms[1], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
      MPI_Allreduce( &o_norms[2][0], &l_rNorms[2], 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
      for( unsigned short l_no = 0; l_no < 3; l_no++ ) o_norms[l_no][0] = l_rNorms[l_no];
#endif

      o_norms[1][0] = std::sqrt(o_norms[1][0]);
    }
};

#endif
