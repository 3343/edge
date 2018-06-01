/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *         Alexander Heinecke (alexander.heinecke AT intel.com)
 *
 * @section LICENSE
 * Copyright (c) 2016-2017, Regents of the University of California
 * Copyright (c) 2016, Intel Corporation
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
 * ADER-DG solver for the elastic wave equations.
 **/
#ifndef ADER_DG_HPP
#define ADER_DG_HPP

#include <limits>
#include <cassert>
#include "constants.hpp"
#include "mesh/common.hpp"
#include "linalg/Mappings.hpp"
#include "TimePred.hpp"
#include "VolInt.hpp"
#include "SurfInt.hpp"
#include "io/Receivers.h"
#include "InternalBoundary.hpp"
#include "FrictionLaws.hpp"

#if defined(PP_T_KERNELS_XSMM) || defined(PP_T_KERNELS_XSMM_DENSE_SINGLE)
#include <libxsmm.h>
#endif

namespace edge {
  namespace elastic {
    namespace solvers { 
      class AderDg;
    }
  }
}

class edge::elastic::solvers::AderDg {
  //private:
    /**
     * Sets the given matrix to zero.
     *
     * @param o_mat matrix which will be set to 0.
     **/
    static void zero( real_base o_mat[N_QUANTITIES][N_ELEMENT_MODES][N_CRUNS] ) {
#if __has_builtin(__builtin_assume_aligned)
      // share alignment with compiler
      (void) __builtin_assume_aligned(o_mat, ALIGNMENT.ELEMENT_MODES.PRIVATE);
#endif

      // reset result to zero
      for( int_qt l_qt = 0; l_qt < N_QUANTITIES; l_qt++ ) {
        for( int_md l_md = 0; l_md < N_ELEMENT_MODES; l_md++ ) {
          for( int_cfr l_cfr = 0; l_cfr < N_CRUNS; l_cfr++ ) {
            o_mat[l_qt][l_md][l_cfr] = 0;
          }
        }
      }
    }

    /**
     * Performs the operation C = A.B with private per-run data in C and B.
     *
     * @param i_a matrix A.
     * @param i_b matrix B.
     * @param o_c matrix C, will bet set to A.B.
     **/
    static void matMulB0( const real_base i_a[N_QUANTITIES][N_ELEMENT_MODES][N_CRUNS],
                          const real_base i_b[N_ELEMENT_MODES][N_ELEMENT_MODES],
                                real_base o_c[N_QUANTITIES][N_ELEMENT_MODES][N_CRUNS] ) {
#if __has_builtin(__builtin_assume_aligned)
      // share alignment with compiler
      (void) __builtin_assume_aligned(i_a, ALIGNMENT.ELEMENT_MODES.PRIVATE);
      (void) __builtin_assume_aligned(o_c, ALIGNMENT.ELEMENT_MODES.PRIVATE);
#endif

      // reset result to zero
      for( int_qt l_qt = 0; l_qt < N_QUANTITIES; l_qt++ ) {
        for( int_md l_md = 0; l_md < N_ELEMENT_MODES; l_md++ ) {
          for( int_cfr l_cfr = 0; l_cfr < N_CRUNS; l_cfr++ ) {
            o_c[l_qt][l_md][l_cfr] = 0;
          }
        }
      }

      for( int_qt l_k = 0; l_k < N_ELEMENT_MODES; l_k++ ) {
        for( int_md l_m = 0; l_m < N_QUANTITIES; l_m++ ) {
          for( int_md l_n = 0; l_n < N_ELEMENT_MODES; l_n++ ) {
            for( int_cfr l_cfr = 0; l_cfr < N_CRUNS; l_cfr++ ) {
              o_c[l_m][l_n][l_cfr] += i_a[l_m][l_k][l_cfr] * i_b[l_k][l_n];
            }
          }
        }
      }
    }

    /**
     * Performs the operation C += A.B with private per-run data in C and B.
     *
     * @param i_a matrix A.
     * @param i_b matrix B.
     * @param i_c matrix C.
     **/
    static void matMulB1( const real_base i_a[N_QUANTITIES][N_QUANTITIES],
                          const real_base i_b[N_QUANTITIES][N_ELEMENT_MODES][N_CRUNS],
                                real_base o_c[N_QUANTITIES][N_ELEMENT_MODES][N_CRUNS] ) {
#if __has_builtin(__builtin_assume_aligned)
      // share alignment with compiler
      (void) __builtin_assume_aligned(i_b, ALIGNMENT.ELEMENT_MODES.PRIVATE);
      (void) __builtin_assume_aligned(o_c, ALIGNMENT.ELEMENT_MODES.PRIVATE);
#endif

      for( int_qt l_k = 0; l_k < N_QUANTITIES; l_k++ ) {
        for( int_md l_m = 0; l_m < N_QUANTITIES; l_m++ ) {
          for( int_md l_n = 0; l_n < N_ELEMENT_MODES; l_n++ ) {
            for( int_cfr l_cfr = 0; l_cfr < N_CRUNS; l_cfr++ ) {
              o_c[l_m][l_n][l_cfr] += i_a[l_m][l_k] * i_b[l_k][l_n][l_cfr];
            }
          }
        }
      }
    }

  public:
    /**
     * Gets the Jacobians of the elastic wave equations in 2 dimensions.
     *
     * @param i_rho density rho.
     * @param i_lam Lame parameter lambda.
     * @param i_mu Lame parameter mu.
     * @param o_A will be set to Jacobians.
     **/
    template <typename T>
    static void getJac2D( T i_rho, T i_lam, T i_mu,
                          T o_A[2][5][5] ) {
      // reset to zero
      for( unsigned short l_di = 0; l_di < 2; l_di++ ) {
        for( int_qt l_q1 = 0; l_q1 < 5; l_q1++ ) {
          for( int_qt l_q2 = 0; l_q2 < 5; l_q2++ ) {
            o_A[l_di][l_q1][l_q2] = 0;
          }
        }
      }

     /*
      * Jacobians in q_t = A(x,y) * q_x + B(x,y) * q_y
      *
      * A:
      *    _____0__1_______2_______________3____4
      *  0|     0, 0,      0, -lambda - 2*mu,   0|0
      *  1|     0, 0,      0,        -lambda,   0|1
      *  2|     0, 0,      0,              0, -mu|2
      *  3|-1/rho, 0,      0,              0,   0|3
      *  4|     0, 0, -1/rho,              0,   0|4
      *    -----0--1-------2---------------3----4
      *
      * B:
      *    0_______1_______2____3_______________4
      *  0|0,      0,      0,   0,        -lambda|0
      *  1|0,      0,      0,   0, -lambda - 2*mu|1
      *  2|0,      0,      0, -mu,              0|2
      *  3|0,      0, -1/rho,   0,              0|3
      *  4|0, -1/rho,      0,   0,              0|4
      *    0-------1-------2----3---------------4
      */
      o_A[0][0][3] = -i_lam - T(2) * i_mu;
      o_A[0][1][3] = -i_lam;
      o_A[0][2][4] = -i_mu;
      o_A[0][3][0] = -T(1) / i_rho;
      o_A[0][4][2] = -T(1) / i_rho;

      o_A[1][0][4] = -i_lam;
      o_A[1][1][4] = -i_lam - T(2) * i_mu;
      o_A[1][2][3] = -i_mu;
      o_A[1][3][2] = -T(1) / i_rho;
      o_A[1][4][1] = -T(1) / i_rho;
    }

    /**
     * Gets the Jacobians of the elastic wave equations in 3 dimensions.
     *
     * @param i_rho density rho.
     * @param i_lam Lame parameter lambda.
     * @param i_mu Lame parameter mu.
     * @param o_A will be set to Jacobians.
     **/
    template <typename T>
    static void getJac3D( T i_rho, T i_lam, T i_mu,
                          T o_A[3][9][9] ) {
      // reset to zero
      for( unsigned short l_di = 0; l_di < 3; l_di++ ) {
        for( int_qt l_q1 = 0; l_q1 < 9; l_q1++ ) {
          for( int_qt l_q2 = 0; l_q2 < 9; l_q2++ ) {
            o_A[l_di][l_q1][l_q2] = 0;
          }
        }
      }

      /*
       * Jacobians in q_t = A(x,y,z) * q_x + B(x,y,z) * q_y + C(x,y,z) * q_z
       *
       * A:
       *   _____0__1__2_______3__4_______5_______________6____7____8
       * 0|     0, 0, 0,      0  0,      0, -lambda - 2*mu,   0,   0|0
       * 1|     0, 0, 0,      0, 0,      0,        -lambda,   0,   0|1
       * 2|     0, 0, 0,      0, 0,      0,        -lambda,   0,   0|2
       * 3|     0, 0, 0,      0, 0,      0,              0, -mu,   0|3
       * 4|     0, 0, 0,      0, 0,      0,              0,   0,   0|4
       * 5|     0, 0, 0,      0, 0,      0,              0,   0, -mu|5
       * 6|-1/rho, 0, 0,      0, 0,      0,              0,   0,   0|6
       * 7|     0, 0, 0, -1/rho, 0,      0,              0,   0,   0|7
       * 8|     0, 0, 0,      0, 0, -1/rho,              0,   0,   0|8
       *   -----0--1--2-------3--4-------5---------------6----7----8
       *
       * B:
       *   0_______1__2_______3_______4__5____6_________________7______8
       * 0|0,      0, 0,      0       0, 0,   0,          -lambda,     0|0
       * 1|0,      0, 0,      0,      0, 0,   0,   -lambda - 2*mu,     0|1
       * 2|0,      0, 0,      0,      0, 0,   0,          -lambda,     0|2
       * 3|0,      0, 0,      0,      0, 0, -mu,                0,     0|3
       * 4|0,      0, 0,      0,      0, 0,   0,                0,   -mu|4
       * 5|0,      0, 0,      0,      0, 0,   0,                0,     0|5
       * 6|0,      0, 0, -1/rho,      0, 0,   0,                0,     0|6
       * 7|0, -1/rho, 0,      0,      0, 0,   0,                0,     0|7
       * 8|0,      0, 0,      0, -1/rho, 0,   0,                0,     0|8
       *   0-------1--2-------3-------4--5---------------6------7------8
       *
       * C:
       *   0__1_______2__3_______4_______5____6____7_______________8
       * 0|0, 0,      0, 0       0,      0,   0,   0,        -lambda|0
       * 1|0, 0,      0, 0,      0,      0,   0,   0,        -lambda|1
       * 2|0, 0,      0, 0,      0,      0,   0,   0, -lambda - 2*mu|2
       * 3|0, 0,      0, 0,      0,      0,   0,   0,              0|3
       * 4|0, 0,      0, 0,      0,      0,   0, -mu,              0|4
       * 5|0, 0,      0, 0,      0,      0, -mu,   0,              0|5
       * 6|0, 0,      0, 0,      0, -1/rho,   0,   0,              0|6
       * 7|0, 0,      0, 0, -1/rho,      0,   0,   0,              0|7
       * 8|0, 0, -1/rho, 0,      0,      0,   0,   0,              0|8
       *   0--1-------2--3-------4-------5----6----7---------------8
       */
      o_A[0][0][6] = -i_lam - T(2) * i_mu;
      o_A[0][1][6] = -i_lam;
      o_A[0][2][6] = -i_lam;
      o_A[0][3][7] = -i_mu;
      o_A[0][5][8] = -i_mu;
      o_A[0][6][0] = -T(1) / i_rho;
      o_A[0][7][3] = -T(1) / i_rho;
      o_A[0][8][5] = -T(1) / i_rho;

      o_A[1][0][7] = -i_lam;
      o_A[1][1][7] = -i_lam - T(2) * i_mu;
      o_A[1][2][7] = -i_lam;
      o_A[1][3][6] = -i_mu;
      o_A[1][4][8] = -i_mu;
      o_A[1][6][3] = -T(1) / i_rho;
      o_A[1][7][1] = -T(1) / i_rho;
      o_A[1][8][4] = -T(1) / i_rho;

      o_A[2][0][8] = -i_lam;
      o_A[2][1][8] = -i_lam;
      o_A[2][2][8] = -i_lam - T(2) * i_mu;
      o_A[2][4][7] = -i_mu;
      o_A[2][5][6] = -i_mu;
      o_A[2][6][5] = -T(1) / i_rho;
      o_A[2][7][4] = -T(1) / i_rho;
      o_A[2][8][2] = -T(1) / i_rho;
    }

    /**
     * Gets the Jacobians of the elastic wave equations.
     *
     * @param i_rho density rho.
     * @param i_lam Lame parameter lambda.
     * @param i_mu Lame parameter mu.
     * @param o_A will be set to Jacobians.
     * @param i_nDim number of dimensions (2 or 3).
     **/
    template <typename T>
    static void getJac( T               i_rho,
                        T               i_lam,
                        T               i_mu,
                        T              *o_A,
                        unsigned short  i_nDim = N_DIM ) {
      if( i_nDim == 2 )      getJac2D( i_rho, i_lam, i_mu, (T (*)[5][5]) o_A);
      else if( i_nDim == 3 ) getJac3D( i_rho, i_lam, i_mu, (T (*)[9][9]) o_A );
      else                   EDGE_LOG_FATAL << "dimensions not supported: " << i_nDim;
    }

    /**
     * Sets up the star matrices, which are a linear combination of the Jacobians.
     *
     * @param i_nElements number of elements.
     * @param i_vertexChars vertex characteristics.
     * @parma i_elVe vertices adjacent to the elements.
     * @param i_bgPars background parameters.
     * @param o_starMatrices will be set to star matrices.
     **/
    static void setupStarM(       int_el           i_nElements,
                            const t_vertexChars   *i_vertexChars,
                            const int_el         (*i_elVe)[ C_ENT[T_SDISC.ELEMENT].N_VERTICES ],
                            const t_bgPars       (*i_bgPars)[1],
                                  t_matStar      (*o_starMatrices)[N_DIM] ) {
      PP_INSTR_FUN("star_matrices")

      // iterate over elements
#ifdef PP_USE_OMP
#pragma omp parallel for num_threads( std::max( parallel::g_nThreads-1, 1 ) )
#endif
      for( int_el l_el = 0; l_el < i_nElements; l_el++ ) {
        // derive jacobians of the pdes
        real_base l_A[N_DIM][N_QUANTITIES][N_QUANTITIES];

        for( unsigned int l_di = 0; l_di < N_DIM; l_di++ ) {
#if defined PP_T_KERNELS_XSMM
          for( int_qt l_nz = 0; l_nz < N_MAT_STAR; l_nz++ ) o_starMatrices[l_el][l_di].mat[l_nz] = 0;
#else
          for( int_qt l_qt1 = 0; l_qt1 < N_QUANTITIES; l_qt1++ ) {
            for( int_qt l_qt2 = 0; l_qt2 < N_QUANTITIES; l_qt2++ ) {
              o_starMatrices[l_el][l_di].mat[l_qt1][l_qt2] = 0;
            }
          }
#endif
        }

        EDGE_CHECK( i_bgPars[l_el][0].rho > TOL.SOLVER );

        getJac( i_bgPars[l_el][0].rho, i_bgPars[l_el][0].lam, i_bgPars[l_el][0].mu, l_A[0][0], N_DIM );

        // derive vertex coords
        real_mesh l_veCoords[3][C_ENT[T_SDISC.ELEMENT].N_VERTICES];
        mesh::common< T_SDISC.ELEMENT >::getElVeCoords( l_el, i_elVe, i_vertexChars, l_veCoords );

        // get inverse jacobian
        real_mesh l_jac[N_DIM][N_DIM];
        linalg::Mappings::evalJac( T_SDISC.ELEMENT, l_veCoords[0], l_jac[0] );

        real_mesh l_jacInv[N_DIM][N_DIM];
#if PP_N_DIM == 2
        linalg::Matrix::inv2x2( l_jac, l_jacInv );
#elif PP_N_DIM == 3
        linalg::Matrix::inv3x3( l_jac, l_jacInv );
#else
#error invalid dimension.
#endif

        // set star matrices
        // iterate over reference dimensions
        for( unsigned int l_d1 = 0; l_d1 < N_DIM; l_d1++ ) {
#if defined PP_T_KERNELS_XSMM

#if PP_N_DIM == 2
          o_starMatrices[l_el][l_d1].mat[0] += l_A[0][0][3] * l_jacInv[0][l_d1];
          o_starMatrices[l_el][l_d1].mat[2] += l_A[0][1][3] * l_jacInv[0][l_d1];
          o_starMatrices[l_el][l_d1].mat[5] += l_A[0][2][4] * l_jacInv[0][l_d1];
          o_starMatrices[l_el][l_d1].mat[6] += l_A[0][3][0] * l_jacInv[0][l_d1];
          o_starMatrices[l_el][l_d1].mat[9] += l_A[0][4][2] * l_jacInv[0][l_d1];

          o_starMatrices[l_el][l_d1].mat[1] += l_A[1][0][4] * l_jacInv[1][l_d1];
          o_starMatrices[l_el][l_d1].mat[3] += l_A[1][1][4] * l_jacInv[1][l_d1];
          o_starMatrices[l_el][l_d1].mat[4] += l_A[1][2][3] * l_jacInv[1][l_d1];
          o_starMatrices[l_el][l_d1].mat[7] += l_A[1][3][2] * l_jacInv[1][l_d1];
          o_starMatrices[l_el][l_d1].mat[8] += l_A[1][4][1] * l_jacInv[1][l_d1];
#elif PP_N_DIM == 3
          o_starMatrices[l_el][l_d1].mat[ 0] += l_A[0][0][6] * l_jacInv[0][l_d1];
          o_starMatrices[l_el][l_d1].mat[ 3] += l_A[0][1][6] * l_jacInv[0][l_d1];
          o_starMatrices[l_el][l_d1].mat[ 6] += l_A[0][2][6] * l_jacInv[0][l_d1];
          o_starMatrices[l_el][l_d1].mat[10] += l_A[0][3][7] * l_jacInv[0][l_d1];
          o_starMatrices[l_el][l_d1].mat[14] += l_A[0][5][8] * l_jacInv[0][l_d1];
          o_starMatrices[l_el][l_d1].mat[15] += l_A[0][6][0] * l_jacInv[0][l_d1];
          o_starMatrices[l_el][l_d1].mat[19] += l_A[0][7][3] * l_jacInv[0][l_d1];
          o_starMatrices[l_el][l_d1].mat[23] += l_A[0][8][5] * l_jacInv[0][l_d1];

          o_starMatrices[l_el][l_d1].mat[ 1] += l_A[1][0][7] * l_jacInv[1][l_d1];
          o_starMatrices[l_el][l_d1].mat[ 4] += l_A[1][1][7] * l_jacInv[1][l_d1];
          o_starMatrices[l_el][l_d1].mat[ 7] += l_A[1][2][7] * l_jacInv[1][l_d1];
          o_starMatrices[l_el][l_d1].mat[ 9] += l_A[1][3][6] * l_jacInv[1][l_d1];
          o_starMatrices[l_el][l_d1].mat[12] += l_A[1][4][8] * l_jacInv[1][l_d1];
          o_starMatrices[l_el][l_d1].mat[16] += l_A[1][6][3] * l_jacInv[1][l_d1];
          o_starMatrices[l_el][l_d1].mat[18] += l_A[1][7][1] * l_jacInv[1][l_d1];
          o_starMatrices[l_el][l_d1].mat[22] += l_A[1][8][4] * l_jacInv[1][l_d1];

          o_starMatrices[l_el][l_d1].mat[ 2] += l_A[2][0][8] * l_jacInv[2][l_d1];
          o_starMatrices[l_el][l_d1].mat[ 5] += l_A[2][1][8] * l_jacInv[2][l_d1];
          o_starMatrices[l_el][l_d1].mat[ 8] += l_A[2][2][8] * l_jacInv[2][l_d1];
          o_starMatrices[l_el][l_d1].mat[11] += l_A[2][4][7] * l_jacInv[2][l_d1];
          o_starMatrices[l_el][l_d1].mat[13] += l_A[2][5][6] * l_jacInv[2][l_d1];
          o_starMatrices[l_el][l_d1].mat[17] += l_A[2][6][5] * l_jacInv[2][l_d1];
          o_starMatrices[l_el][l_d1].mat[20] += l_A[2][7][4] * l_jacInv[2][l_d1];
          o_starMatrices[l_el][l_d1].mat[21] += l_A[2][8][2] * l_jacInv[2][l_d1];
#else
#error not defined
#endif

#else
          for( unsigned int l_d2 = 0; l_d2 < N_DIM; l_d2++ ) {
            for( int_qt l_qt1 = 0; l_qt1 < N_QUANTITIES; l_qt1++ ) {
              for( int_qt l_qt2 = 0; l_qt2 < N_QUANTITIES; l_qt2++ ) {
                o_starMatrices[l_el][l_d1].mat[l_qt1][l_qt2] += l_A[l_d2][l_qt1][l_qt2] * l_jacInv[l_d2][l_d1];
               }
             }
           }
#endif
         }
      }
    }

    /**
     * Local step: Cauchy Kowalevski + volume.
     *
     * @param i_first first element considered.
     * @param i_nElements number of elements.
     * @param i_time time of the initial DOFs.
     * @param i_dT time step.
     * @param i_firstSpRp first sparse rupture-element.
     * @param i_firstSpRe first sparse receiver entity.
     * @param i_elFa faces adjacent to the elements.
     * @param i_faChars face characteristics.
     * @param i_elChars element characteristics.
     * @param i_dg const DG data.
     * @param i_starM star matrices.
     * @param i_fluxSolvers flux solvers for the local element's contribution.
     * @param io_dofs DOFs.
     * @param o_tInt will be set to time integrated DOFs.
     * @param o_tRup will be set to DOFs for rupture elements.
     * @param io_recvs will be updated with receiver info.
     * @param i_kernels kernels of XSMM-library for the local step (if enabled).
     *
     * @paramt TL_T_INT_LID integer type of local entity ids.
     * @paramt TL_T_REAL floating point type.
     * @paramt TL_T_MM matrix-matrix multiplication kernels.
     **/
    template < typename TL_T_INT_LID,
               typename TL_T_REAL,
               typename TL_T_MM >
    static void local( TL_T_INT_LID                     i_first,
                       TL_T_INT_LID                     i_nElements,
                       double                           i_time,
                       double                           i_dT,
                       TL_T_INT_LID                     i_firstSpRp,
                       TL_T_INT_LID                     i_firstSpRe,
                       TL_T_INT_LID            const (* i_elFa)[ C_ENT[T_SDISC.ELEMENT].N_FACES ],
                       t_faceChars             const  * i_faChars,
                       t_elementChars          const  * i_elChars,
                       t_dg                    const  & i_dg,
                       t_matStar               const (* i_starM)[N_DIM],
                       t_fluxSolver            const (* i_fluxSolvers)[ C_ENT[T_SDISC.ELEMENT].N_FACES ],
                       TL_T_REAL                     (* io_dofs)[N_QUANTITIES][N_ELEMENT_MODES][N_CRUNS],
                       TL_T_REAL                     (* o_tInt)[N_QUANTITIES][N_ELEMENT_MODES][N_CRUNS],
                       TL_T_REAL                     (* o_tRup)[N_QUANTITIES][N_ELEMENT_MODES][N_CRUNS],
                       edge::io::Receivers            & io_recvs,
                       TL_T_MM const                  & i_mm ) {
#if __has_builtin(__builtin_assume_aligned)
      // share alignment with compiler
      (void) __builtin_assume_aligned(io_dofs, ALIGNMENT.ELEMENT_MODES.PRIVATE);
      (void) __builtin_assume_aligned(o_tInt,  ALIGNMENT.ELEMENT_MODES.PRIVATE);
#endif

      // counter of elements with faces having rupture physics
      TL_T_INT_LID l_elRp = i_firstSpRp;

      // counter for receivers
      unsigned int l_enRe = i_firstSpRe;

      // temporary data structurre for product for two-way mult and receivers
      TL_T_REAL (*l_tmpEl)[N_ELEMENT_MODES][N_CRUNS] = parallel::g_scratchMem->tRes;

      // buffer for derivatives
      TL_T_REAL (*l_derBuffer)[N_QUANTITIES][N_ELEMENT_MODES][N_CRUNS] = parallel::g_scratchMem->dBuf;

      // iterate over all elements
      for( TL_T_INT_LID l_el = i_first; l_el < i_first+i_nElements; l_el++ ) {
        if( (i_elChars[l_el].spType & RUPTURE) != RUPTURE ) {}
        else {
          // save DOFs of the rupture element
          for( unsigned short l_qt = 0; l_qt < N_QUANTITIES; l_qt++ )
            for( unsigned short l_md = 0; l_md < N_ELEMENT_MODES; l_md++ )
              for( unsigned short l_ru = 0; l_ru < N_CRUNS; l_ru++ )
                o_tRup[l_elRp][l_qt][l_md][l_ru] = io_dofs[l_el][l_qt][l_md][l_ru];

          // increase rupture element counter
          l_elRp++;
        }

        /*
         * compute ader time integration
         */
        TimePred< T_SDISC.ELEMENT,
                  N_QUANTITIES,
                  ORDER,
                  ORDER,
                  N_CRUNS >::ck( (TL_T_REAL)  i_dT,
                                              i_dg.mat.stiffT,
                                            &(i_starM[l_el][0].mat), // TODO: fix struct
                                              io_dofs[l_el],
                                              i_mm,
                                              l_tmpEl,
                                              l_derBuffer,
                                              o_tInt[l_el] );

        /*
         * Write receivers (if required)
         */
        if( !( (i_elChars[l_el].spType & RECEIVER) == RECEIVER) ) {} // no receivers in the current element
        else { // we have receivers in the current element
          while( true ) { // iterate of possible multiple receiver-ouput per time step
            double l_rePt = io_recvs.getRecvTimeRel( l_enRe, i_time, i_dT );
            if( !(l_rePt >= 0) ) break;
            else {
              TL_T_REAL l_rePts = l_rePt;
              // eval time prediction at the given point
              TimePred< T_SDISC.ELEMENT,
                        N_QUANTITIES,
                        ORDER,
                        ORDER,
                        N_CRUNS >::evalTimePrediction(  1,
                                                       &l_rePts,
                                                        l_derBuffer,
                          (TL_T_REAL (*)[N_QUANTITIES][N_ELEMENT_MODES][N_CRUNS])l_tmpEl );

              // write this time prediction
              io_recvs.writeRecvAll( l_enRe, l_tmpEl );
            }
          }
          l_enRe++;
        }

        /*
         * compute volume contribution
         */
        VolInt< T_SDISC.ELEMENT,
                N_QUANTITIES,
                ORDER,
                N_CRUNS >::apply(   i_dg.mat.stiff,
                                  & i_starM[l_el][0].mat, // TODO: fix struct
                                    o_tInt[l_el],
                                    i_mm,
                                    io_dofs[l_el],
                                    l_tmpEl );

        /*
         * prefetches for next iteration
         */
        const TL_T_REAL (* l_preDofs)[N_ELEMENT_MODES][N_CRUNS] = nullptr;
        const TL_T_REAL (* l_preTint)[N_ELEMENT_MODES][N_CRUNS] = nullptr;
        if( l_el < i_first+i_nElements-1 ) {
          l_preDofs = io_dofs[l_el+1];
          l_preTint = o_tInt[l_el+1];
        }
        else {
          l_preDofs = io_dofs[l_el];
          l_preTint = o_tInt[l_el];
        }

         /*
          * compute local surface contribution
          */
        // reuse derivative buffer
        TL_T_REAL (*l_tmpFa)[N_QUANTITIES][N_FACE_MODES][N_CRUNS] =
          (TL_T_REAL (*)[N_QUANTITIES][N_FACE_MODES][N_CRUNS]) parallel::g_scratchMem->dBuf;
        // call kernel
        SurfInt< T_SDISC.ELEMENT,
                 N_QUANTITIES,
                 ORDER,
                 ORDER,
                 N_CRUNS >::local( i_dg.mat.fluxL,
                                   i_dg.mat.fluxT,
                                   ( TL_T_REAL (*)[N_QUANTITIES][N_QUANTITIES] )  ( i_fluxSolvers[l_el][0].solver[0] ), // TODO: fix struct
                                   o_tInt[l_el],
                                   i_mm,
                                   io_dofs[l_el],
                                   l_tmpFa,
                                   l_preDofs,
                                   l_preTint );
      }
    }

    /**
     * Solves rupture physics for the given faces.
     *
     * @param i_first first rupture faces.
     * @param i_nFaces number of faces.
     * @param i_firstReSp first sparse id of the receivers.
     * @param i_time current time of the faces.
     * @param i_dT time step of the two adjacent elements.
     * @param i_dg discontinuous galerkin data structures.
     * @param i_starM star matrices.
     * @param i_faElSpRp adjacency information from sparse rupture faces to sparse rupture elements.
     * @param i_spTypesFaSp sparse types of the sparse rupture faces.
     * @param i_frictionGlobal global data of the friction law.
     * @param i_frictionFace face-local data of the friction law.
     * @param i_frictionQuadPoint data of the friction law local to the spatial quadrature points of the faces.
     * @param i_solver solvers used for internal boundaries. [0]: trafo to face-align coords, [1]: middle states [2]: trafo to physical coords and flux (left side) [3]: trafo to physical coords and flux (right side).
     * @param i_tDofs DOFs of the sparse rupture elements.
     * @param o_updates will be set to element-updates resulting from the flux computation. [0]: left side, [1]: right side.
     * @param io_recvsQuad receivers at quadrature points.
     *
     * @paramt TL_T_INT_LID integer type of local entity ids.
     * @paramt TL_T_REAL type used for floating point arithmetic.
     * @paramt TL_T_FRI_GL struct representing global friction data.
     * @paramt TL_T_FRI_FA struct representing face-local friction data.
     * @paramt TL_T_FRI_QP struct representing quad-point-local friction data.
     * @paramt TL_T_INT_SP integer type of the sparse type.
     * @paramt TL_T_RECVQ type of the quadrature receivers.
     **/
    template< typename TL_T_INT_LID,
              typename TL_T_REAL,
              typename TL_T_FRI_GL,
              typename TL_T_FRI_FA,
              typename TL_T_FRI_QP,
              typename TL_T_INT_SP,
              typename TL_T_RECVQ >
    static void rupture( TL_T_INT_LID                     i_first,
                         TL_T_INT_LID                     i_nFaces,
                         TL_T_INT_LID                     i_firstSpRe,
                         TL_T_REAL                        i_time,
                         TL_T_REAL               const    i_dT,
                         t_dg                    const  & i_dg,
                         t_matStar               const (* i_starM)[N_DIM],
                         TL_T_INT_LID            const (* i_faElSpRp)[2],
                         t_InternalBoundaryFace<
                           TL_T_REAL,
                           TL_T_INT_SP
                         >                       const  * i_iBnd,
                         TL_T_FRI_GL                    & i_frictionGlobal,
                         TL_T_FRI_FA                    * i_frictionFace,
                         TL_T_FRI_QP                   (* i_frictionQuadPoint)[N_FACE_QUAD_POINTS],
                         TL_T_REAL               const (* i_solvers)[4][N_QUANTITIES][N_QUANTITIES],
                         TL_T_REAL               const (* i_tDofs)[N_QUANTITIES][N_ELEMENT_MODES][N_CRUNS],
                         TL_T_REAL                     (* o_updates)[2][N_QUANTITIES][N_ELEMENT_MODES][N_CRUNS],
                         TL_T_RECVQ                     & io_recvsQuad
#if defined (PP_T_KERNELS_XSMM) || defined (PP_T_KERNELS_XSMM_DENSE_SINGLE)
#if PP_PRECISION == 64
                      ,const libxsmm_dmmfunction *i_kernels
#else
                      ,const libxsmm_smmfunction *i_kernels
#endif
#endif
                     ) {
      EDGE_LOG_FATAL << "removed quadrature-based rupture physics";
    }

    /**
     * Performs the neighboring updates of the ADER-DG scheme.
     *
     * @param i_first first element considered.
     * @param i_nElements number of elements.
     * @param i_dg constant DG data.
     * @param i_faChars face characteristics.
     * @param i_fluxSolvers flux solvers for the neighboring elements' contribution.
     * @param i_elFa elements' adjacent faces.
     * @param i_elFaEl face-neighboring elements.
     * @param i_faElSpRp adjacency information from sparse rupture faces to sparse rupture elements.
     * @param i_elFaSpRp adjacnecy information from sparse rupture elements to sparse rupture faces.
     * @param i_fIdElFaEl local face ids of face-neighboring elememts.
     * @param i_vIdElFaEl local vertex ids w.r.t. the shared face from the neighboring elements' perspsective.
     * @param i_tInt time integrated degrees of freedom.
     * @param i_updatesSpRp surface updates resulting from rupture physics.
     * @param io_dofs DOFs which will be updated with neighboring elements' contribution.
     * @param i_kernels kernels of XSMM-library for the neighboring step (if enabled).
     *
     * @paramt TL_T_INT_LID integer type of local entity ids.
     * @paramt TL_T_REAL type used for floating point arithmetic.
     * @paramt TL_T_MM type of the matrix-matrix multiplication kernels.
     **/
    template< typename TL_T_INT_LID,
              typename TL_T_REAL,
              typename TL_T_MM >
    static void neigh( TL_T_INT_LID            i_first,
                       TL_T_INT_LID            i_nElements,
                       TL_T_INT_LID            i_firstSpRp,
                       t_dg           const  & i_dg,
                       t_faceChars    const  * i_faChars,
                       t_fluxSolver   const (* i_fluxSolvers)[ C_ENT[T_SDISC.ELEMENT].N_FACES ],
                       TL_T_INT_LID   const (* i_elFa)[C_ENT[T_SDISC.ELEMENT].N_FACES],
                       TL_T_INT_LID   const (* i_elFaEl)[C_ENT[T_SDISC.ELEMENT].N_FACES],
                       TL_T_INT_LID  const  (* i_faElSpRp)[2],
                       TL_T_INT_LID  const  (* i_elFaSpRp)[C_ENT[T_SDISC.ELEMENT].N_FACES],
                       unsigned short const (* i_fIdElFaEl)[C_ENT[T_SDISC.ELEMENT].N_FACES],
                       unsigned short const (* i_vIdElFaEl)[C_ENT[T_SDISC.ELEMENT].N_FACES],
                       TL_T_REAL      const (* i_tInt)[N_QUANTITIES][N_ELEMENT_MODES][N_CRUNS],
                       TL_T_REAL      const (* i_updatesSpRp)[2][N_QUANTITIES][N_ELEMENT_MODES][N_CRUNS],
                       TL_T_REAL            (* io_dofs)[N_QUANTITIES][N_ELEMENT_MODES][N_CRUNS],
                       TL_T_MM const         & i_mm ) {
#if __has_builtin(__builtin_assume_aligned)
      // share alignment with compiler
      (void) __builtin_assume_aligned(i_tInt,  ALIGNMENT.ELEMENT_MODES.PRIVATE);
      (void) __builtin_assume_aligned(io_dofs, ALIGNMENT.ELEMENT_MODES.PRIVATE);
#endif

      // temporary product for three-way mult
        TL_T_REAL (*l_tmpFa)[N_QUANTITIES][N_FACE_MODES][N_CRUNS] =
          (TL_T_REAL (*)[N_QUANTITIES][N_FACE_MODES][N_CRUNS]) parallel::g_scratchMem->dBuf;

      // iterate over elements
      for( TL_T_INT_LID l_el = i_first; l_el < i_first+i_nElements; l_el++ ) {

        // add neighboring contribution
        for( TL_T_INT_LID l_fa = 0; l_fa < C_ENT[T_SDISC.ELEMENT].N_FACES; l_fa++ ) {
          TL_T_INT_LID l_faId = i_elFa[l_el][l_fa];
          TL_T_INT_LID l_ne;

          // determine flux matrix id
          unsigned short l_fId = SurfInt< T_SDISC.ELEMENT,
                                          N_QUANTITIES,
                                          ORDER,
                                          ORDER,
                                          N_CRUNS >::fMatId( i_vIdElFaEl[l_el][l_fa],
                                                             i_fIdElFaEl[l_el][l_fa] );

          if( (i_faChars[l_faId].spType & OUTFLOW) != OUTFLOW ) {
            // derive neighbor
            if( (i_faChars[l_faId].spType & FREE_SURFACE) != FREE_SURFACE )
              l_ne = i_elFaEl[l_el][l_fa];
            else
              l_ne = l_el;

            /*
             * prefetches
             */
            const TL_T_REAL (* l_pre)[N_ELEMENT_MODES][N_CRUNS] = nullptr;
            TL_T_INT_LID l_neUp = std::numeric_limits<TL_T_INT_LID>::max();
            // prefetch for the upcoming surface integration of this element
            if( l_fa < C_ENT[T_SDISC.ELEMENT].N_FACES-1 ) l_neUp = i_elFaEl[l_el][l_fa+1];
            // first surface integration of the next element
            else if( l_el < i_first+i_nElements-1 ) l_neUp = i_elFaEl[l_el+1][0];

            // only proceed with adjacent data if the element exists
            if( l_neUp != std::numeric_limits<TL_T_INT_LID>::max() ) l_pre = i_tInt[l_neUp];
            // next element data in case of boundary conditions
            else if( l_el < i_first+i_nElements-1 )              l_pre = io_dofs[l_el+1];
            // default to element data to avoid performance penality
            else                                                 l_pre = io_dofs[l_el];

            /*
             * solve
             */
            SurfInt< T_SDISC.ELEMENT,
                     N_QUANTITIES,
                     ORDER,
                     ORDER,
                     N_CRUNS >::neigh( ((i_faChars[l_faId].spType & FREE_SURFACE) != FREE_SURFACE ) ? i_dg.mat.fluxN[l_fId] :
                                                                                                      i_dg.mat.fluxL[l_fa],
                                       i_dg.mat.fluxT[l_fa],
                                       ( TL_T_REAL (*)[N_QUANTITIES] )  ( i_fluxSolvers[l_el][l_fa].solver[0] ), // TODO: fix struct
                                       i_tInt[l_ne],
                                       i_mm,
                                       io_dofs[l_el],
                                       l_tmpFa,
                                       l_pre,
                                       l_fa,
                                       ((i_faChars[l_faId].spType & FREE_SURFACE) != FREE_SURFACE ) ? l_fId + C_ENT[T_SDISC.ELEMENT].N_FACES: l_fa );
        }
      }
    }
  }
};

#endif
