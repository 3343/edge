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
 * Setup for convergence studies of the elastic wave equations.
 **/

#include "constants.hpp"
#include <cmath>
#include "../common.hpp"
#include "mesh/common.hpp"
#include "dg/QuadraturePoints.h"
#include "data/layout.hpp"

namespace edge {
  namespace elastic {
    namespace setups {
      class Convergence;
    }
  }
}

class edge::elastic::setups::Convergence {
  private:
    /**
     * Derices the initial quantities for the plane wave convergence setup.
     * Travelling direction is n=(1,1). p-wave velocity is 2, s-wave velocity is 1. density is 1.
     *
     * @param i_x position in x-direction.
     * @param i_y position in y-direction.
     * @param o_quantities will be set to initial quantities (sigma_xx, sigma_yy, sigma_xy, u, v) at the given point.
     **/
    static void getPlaneWavesPoint2D( real_mesh i_x,
                                      real_mesh i_y,
                                      real_base o_quantities[5] ) {
      // wave number
      real_base l_k = ( (real_base) 2.0 * M_PI ) / (real_base) 25.0;

      // derive scaling
      real_base l_scale = std::sin( i_x * l_k + i_y * l_k );

      // lame parameters
      real_base l_lambda = 2; real_base l_mu = 1; real_base l_rho = 1;

      // compute wave speeds
      real_base l_cp = common::getVelP( l_rho, l_lambda, l_mu );

      real_base l_cs = common::getVelS( l_rho, l_mu );

      // derive 2nd and fifth eigenvector
      real_base l_r2[5]; real_base l_r5[5];
      real_base l_nx = 1 / std::sqrt(2);
      real_base l_ny = 1 / std::sqrt(2);

      l_r2[0] = -2.0 * l_mu * l_nx * l_ny;
      l_r2[1] =  2.0 * l_mu * l_nx * l_ny;
      l_r2[2] =  l_mu * ( l_nx*l_nx - l_ny*l_ny);
      l_r2[3] = -l_ny * l_cs;
      l_r2[4] =  l_nx * l_cs;

      l_r5[0] = l_lambda + 2.0 * l_mu * l_nx * l_nx;
      l_r5[1] = l_lambda + 2.0 * l_mu * l_ny * l_ny;
      l_r5[2] = 2.0 * l_mu * l_nx * l_ny;
      l_r5[3] = -l_nx * l_cp;
      l_r5[4] = -l_ny * l_cp;

      // set result
      for( unsigned int l_q = 0; l_q < 5; l_q++ )  o_quantities[l_q] = ( l_r2[l_q] + l_r5[l_q] ) * l_scale;
    }

    /**
     * Derices the initial quantities for the plane wave convergence setup.
     * Travelling direction is n=(1,1,1). p-wave velocity is 2, s-wave velocity is 1. density is 1.
     *
     * @param i_x position in x-direction.
     * @param i_y position in y-direction.
     * @param i_z position in z-direciton.
     * @param o_quantities will be set to initial quantities (sigma_xx, sigma_yy, sigma_zz, sigma_xy, sigma_xz, sigma_yz, u, v, w) at the given point.
     **/
    static void getPlaneWavesPoint3D( real_mesh i_x,
                                      real_mesh i_y,
                                      real_mesh i_z,
                                      real_base o_quantities[9] ) {
      // wave number
      real_base l_k = ( (real_base) 2.0 * M_PI ) / (real_base) 100.0;

      // derive scaling
      real_base l_scale = std::sin( i_x * l_k + i_y * l_k  + i_z * l_k);

      // derive first and eights eigenvector
      real_base l_r1[9]; real_base l_r8[9];

      l_r1[0] = 4.0 / std::sqrt( 3 );
      l_r1[1] = 4.0 / std::sqrt( 3 );
      l_r1[2] = 4.0 / std::sqrt( 3 );
      l_r1[3] = 1.0 / std::sqrt( 3 );
      l_r1[4] = 1.0 / std::sqrt( 3 );
      l_r1[5] = 1.0 / std::sqrt( 3 );
      l_r1[6] = 1.0;
      l_r1[7] = 1.0;
      l_r1[8] = 1.0;

      // remark: choice of #8 is non-unique since there are two eigenvalues with speed c_s
      l_r8[0] =  2.0 / std::sqrt( 3 );
      l_r8[1] =  0.0;
      l_r8[2] = -2.0 / std::sqrt( 3 );
      l_r8[3] =  1.0 / std::sqrt( 3 );
      l_r8[4] = -1.0 / std::sqrt( 3 );
      l_r8[5] =  0.0;
      l_r8[6] = -1.0;
      l_r8[7] =  0.0;
      l_r8[8] =  1.0;

      // set result
      for( unsigned int l_q = 0; l_q < 9; l_q++ )  o_quantities[l_q] = ( l_r1[l_q] + l_r8[l_q] ) * l_scale;
    }

    /**
     * Gets the initial plane wave quantities at the element's physical quadrature points.
     *
     * @param i_el element which solution at the quad points is queried.
     * @param i_connElVe connectivity info from elements to vertices.
     * @param i_vertexChars vertex characteristics.
     * @param i_offsetX offset in x-direction.
     * @param i_offsetY offset in y-direction.
     * @param i_offsetZ offset in z-direction.
     * @param o_weights will be set to weights of quad points.
     * @param o_q0 will be set to initial solution at quad points.
     */
    static void getPlaneWavesQp(       int_el                     i_el,
                                 const int_el                   (*i_connElVe)[C_ENT[T_SDISC.ELEMENT].N_VERTICES],
                                 const t_vertexChars             *i_vertexChars,
                                       real_mesh                  i_offsetX,
                                       real_mesh                  i_offsetY,
                                       real_mesh                  i_offsetZ,
                                       std::vector< real_mesh >  &o_weights,
                                       std::vector< real_base >  *o_q0 ) {
      PP_INSTR_FUN("get_pwav_qp")

      // get elements vertices
      real_mesh l_veCoords[3][C_ENT[T_SDISC.ELEMENT].N_VERTICES];
      mesh::common< T_SDISC.ELEMENT >::getElVeCoords( i_el, i_connElVe, i_vertexChars, l_veCoords );

      // get gaussian quadrature points and their weigths
      std::vector< real_mesh > l_ptsX, l_ptsY, l_ptsZ;

      dg::QuadraturePoints::getQpts( T_SDISC.ELEMENT,
                                     ORDER+1,
                                     l_veCoords,
                                     l_ptsX, l_ptsY, l_ptsZ, o_weights );

     for( int_qt l_qt = 0; l_qt < N_QUANTITIES; l_qt++ ) o_q0[l_qt].resize( o_weights.size() );

      // iterate over quad points and get solution
      for( unsigned int l_qp = 0; l_qp < l_ptsX.size(); l_qp++ ) {
        real_mesh l_x = l_ptsX[l_qp] + i_offsetX;
        real_mesh l_y = l_ptsY[l_qp] + i_offsetY;
#if PP_N_DIM == 3
        real_mesh l_z = l_ptsZ[l_qp] + i_offsetZ;
#endif

        // all quantities at quad point in AoS
        real_base l_q0tmp[N_QUANTITIES];
#if PP_N_DIM == 2
        getPlaneWavesPoint2D( l_x, l_y, l_q0tmp );
#elif PP_N_DIM == 3
        getPlaneWavesPoint3D( l_x, l_y, l_z, l_q0tmp );
#else
        assert( false );
#endif

        // soa
        for( int_qt l_qt = 0; l_qt < N_QUANTITIES; l_qt++ ) o_q0[l_qt][l_qp] = l_q0tmp[l_qt];
      }
    }

  public:
    /**
     * Sets the given DOFs to plane waves.
     *
     * @param i_cfr concurrent forward run.
     * @param i_basis basis.
     * @param i_nElements number of elements.
     * @param i_connElVe connectivity information from elements to vertices.
     * @param i_vertexChars vertex characteristics.
     * @param i_elementChars element characteristics.
     * @param o_bgPars will be set to background parameters.
     * @param o_dofs will be set to plane waves.
     * @param i_offsetX offset of the waves in x-direction.
     * @param i_offsetY offset of the waves in y-direction.
     * @param i_offsetZ offset of the waves in z-direction (only used for 3D elastics).
     **/
    static void setPlaneWaves(       int_cfr                         i_cfr,
                               const dg::Basis                      &i_basis,
                                     int_el                          i_nElements,
                               const int_el                        (*i_connElVe)[C_ENT[T_SDISC.ELEMENT].N_VERTICES],
                               const t_vertexChars                  *i_vertexChars,
                               const t_elementChars                 *i_elementChars,
                                     t_bgPars                      (*o_bgPars)[1],
                                     real_base                     (*o_dofs)[N_QUANTITIES][N_ELEMENT_MODES][N_CRUNS],
                                     real_mesh                       i_offsetX = 0,
                                     real_mesh                       i_offsetY = 0,
                                     real_mesh                       i_offsetZ = 0 ) {
       PP_INSTR_FUN("set_pwav")
#ifdef PP_USE_OMP
#pragma omp parallel for
#endif
      for( int_el l_el = 0; l_el < i_nElements; l_el++ ) {
        // set background parameters
        o_bgPars[l_el][0].rho = 1.0;
        o_bgPars[l_el][0].lam = 2.0;
        o_bgPars[l_el][0].mu  = 1.0;

        // evaluated initial quantities at quad points
        std::vector< real_mesh > l_weights;
        std::vector< real_base > l_q0[N_QUANTITIES];
        getPlaneWavesQp( l_el,
                         i_connElVe,
                         i_vertexChars,
                         i_offsetX,
                         i_offsetY,
                         i_offsetZ,
                         l_weights,
                         l_q0 );

        for( int_qt l_qt = 0; l_qt < N_QUANTITIES; l_qt++ ) {
          // map to modes
          real_base l_modes[N_ELEMENT_MODES];
          i_basis.qpts2modal( &l_q0[l_qt][0],
                              ORDER+1,
                              l_modes );

          // scatter to cfr-structure
          for( int_md l_md = 0; l_md < N_ELEMENT_MODES; l_md++ ) {
            o_dofs[l_el][l_qt][l_md][i_cfr] = l_modes[l_md];
          }
        }
      }
    }

    /**
     * Gets the error norms for the plance convergence setup.
     *
     * @param i_cfr concurrent forward run which norms are returned.
     * @param i_basis basis.
     * @param i_inMap index mapping between computational and mesh data.
     * @param i_elLayout element layout.
     * @param i_connect connectivity information.
     * @param i_vertexChars vertex characteristics.
     * @param i_elementChars element characteristics.
     * @param i_dofs DOFs.
     * @param o_norms will be set to norms. [*][]: L1, L2, Linf, [][*]: sigma_xx, sigma_yy, (sigma_zz), sigma_xy, (sigma_xz, sigmay_z), u, v, (w); parentheses: 3D.
     * @param i_offsetX x-offset of the integration points' positions in the setup of the convergence benchmark.
     * @param i_offsetY y-offset of the integration points' positions in the setup of the convergence benchmark.
     * @param i_offsetZ z-offset of the integration points' positions in the setup of the convergence benchmark.
     **/
    static void getPlaneErrorNorms(        int_cfr          i_cfr,
                                    const  dg::Basis       &i_basis,
                                    const  t_inMap         *i_inMap,
                                    const  t_enLayout      &i_elLayout,
                                    const  t_connect       &i_connect,
                                    const  t_vertexChars   *i_vertexChars,
                                    const  t_elementChars  *i_elementChars,
                                    const  real_base      (*i_dofs)[N_QUANTITIES][N_ELEMENT_MODES][N_CRUNS],
                                           double           o_norms[3][N_QUANTITIES],
                                           real_mesh        i_offsetX = 0,
                                           real_mesh        i_offsetY = 0,
                                           real_mesh        i_offsetZ = 0 ) {
      // reset norms
      for( unsigned int l_norm = 0; l_norm < 3; l_norm++ ) {
        for( int_qt l_q = 0; l_q < N_QUANTITIES; l_q++ ) {
          o_norms[l_norm][l_q] = 0.0;
        }
      }

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

          /**
           * Get exact data.
           **/
          // evaluated initial quantities at quad points
          std::vector< real_base > l_q0[N_QUANTITIES];
          std::vector< real_mesh > l_weights;
          getPlaneWavesQp( l_el,
                           i_connect.elVe,
                           i_vertexChars,
                           i_offsetX,
                           i_offsetY,
                           i_offsetZ,
                           l_weights,
                           l_q0 );

          /**
           * Get numerical solution
           **/
          // get gaussian quadrature points and their weigths for the reference element
          std::vector< real_mesh > l_soPtsX, l_soPtsY, l_soPtsZ, l_soWeights;

          dg::QuadraturePoints::getQpts( T_SDISC.ELEMENT,
                                         ORDER+1,
                                         C_REF_ELEMENT.VE.ENT[T_SDISC.ELEMENT],
                                         l_soPtsX, l_soPtsY, l_soPtsZ, l_soWeights );

          assert( l_soPtsX.size() == l_weights.size() );

          // solution at quad points
          std::vector< real_base > l_qS[N_QUANTITIES];
          for( int_qt l_qt = 0; l_qt < N_QUANTITIES; l_qt++ ) l_qS[l_qt].resize( l_soWeights.size() );

          // gather modes
          real_base l_modes[N_QUANTITIES][N_ELEMENT_MODES];
          for( int_qt l_qt = 0; l_qt < N_QUANTITIES; l_qt++ ) {
            for( int_md l_md = 0; l_md < N_ELEMENT_MODES; l_md++ ) {
              l_modes[l_qt][l_md] = i_dofs[l_el][l_qt][l_md][i_cfr];
            }
          }

          // iterate over quad points and get the solution
          for( int_qt l_qt = 0; l_qt < N_QUANTITIES; l_qt++ ) {
            for( unsigned int l_qp = 0; l_qp < l_soPtsX.size(); l_qp++ ) {
              l_qS[l_qt][l_qp] = i_basis.modal2refPtVal( ORDER+1, l_qp, l_modes[l_qt] );
            }
          }

          // add errors
          for( unsigned int l_qp = 0; l_qp < l_weights.size(); l_qp++ ) {
            for( unsigned int l_qt = 0; l_qt < N_QUANTITIES; l_qt++ ) {
              double l_diff = std::abs( l_q0[l_qt][l_qp] - l_qS[l_qt][l_qp] );

              o_norms[0][l_qt] += l_diff *          l_weights[l_qp];
              o_norms[1][l_qt] += l_diff * l_diff * l_weights[l_qp];
              o_norms[2][l_qt] = std::max( o_norms[2][l_qt], l_diff );
            }
          }
        }
      }

      // do the MPI-reductions
#ifdef PP_USE_MPI
      double l_rNorms[3][N_QUANTITIES];
      for( unsigned int l_qt = 0; l_qt < N_QUANTITIES; l_qt++ ) {
        MPI_Allreduce( &o_norms[0][l_qt], &l_rNorms[0][l_qt], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
        MPI_Allreduce( &o_norms[1][l_qt], &l_rNorms[1][l_qt], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
        MPI_Allreduce( &o_norms[2][l_qt], &l_rNorms[2][l_qt], 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
        for( unsigned short l_no = 0; l_no < 3; l_no++ ) o_norms[l_no][l_qt] = l_rNorms[l_no][l_qt];
      }
#endif

      for( unsigned int l_qt = 0; l_qt < N_QUANTITIES; l_qt++ ) {
        o_norms[1][l_qt] = std::sqrt(o_norms[1][l_qt]);
      }
    }
};
