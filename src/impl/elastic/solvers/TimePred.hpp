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
 * Time predictions through the ADER scheme for the elastic wave equations.
 **/

#ifndef TIME_PRED_HPP
#define TIME_PRED_HPP

#include "constants.hpp"
#include "linalg/Matrix.h"

namespace edge {
  namespace elastic {
    namespace solvers {
      template <t_entityType TL_TYPE_ELEMENT, unsigned short TL_N_QUANTITIES, unsigned short TL_ORDER, unsigned short TL_N_CRUNS, unsigned short TL_N_POINTS_TIME>
      class TimePred;
    }
  }
}

/**
 * ADER-related functions:
 *   1) Computation of time predictions (time derivatives and time integrated DOFs) through the Cauchy–Kowalevski procedure.
 *   2) Evaluation of time prediction at specific points in time.
 *
 * @paramt TL_TYPE_ELEMENT element type.
 * @paramt TL_N_QUANTITIES number of quantities.
 * @paramt TL_ORDER order of the ADER scheme.
 * @paramt TL_N_CRUNS number of concurrent forward runs (fused simulations)
 * @paramt TL_N_POINTS_TIME number of points in the evaluation of the time predictor.
 **/
template <t_entityType TL_TYPE_ELEMENT, unsigned short TL_N_QUANTITIES, unsigned short TL_ORDER, unsigned short TL_N_CRUNS, unsigned short TL_N_POINTS_TIME=1>
class edge::elastic::solvers::TimePred {
  private:
    // assemble derived template parameters
    //! dimension of the element
    static unsigned short const TL_N_DIM           = C_ENT[TL_TYPE_ELEMENT].N_DIM;
    //! number of element modes
    static unsigned short const TL_N_ELEMENT_MODES = CE_N_ELEMENT_MODES( TL_TYPE_ELEMENT, TL_ORDER );

    /**
     * Sets the given matrix to zero.
     *
     * @param o_mat matrix which will be set to 0.
     *
     * @paramt floating point precision of the matrix.
     **/
    template < typename TL_T_REAL >
    static void zero( TL_T_REAL o_mat[TL_N_QUANTITIES][TL_N_ELEMENT_MODES][TL_N_CRUNS] ) {
      // reset result to zero
      for( int_qt l_qt = 0; l_qt < TL_N_QUANTITIES; l_qt++ ) {
        for( int_md l_md = 0; l_md < TL_N_ELEMENT_MODES; l_md++ ) {
          for( int_cfr l_cfr = 0; l_cfr < TL_N_CRUNS; l_cfr++ ) {
            o_mat[l_qt][l_md][l_cfr] = 0;
          }
        }
      }
    }

  public:
    /**
     * Applies the Cauchy–Kowalevski procedure (vanilla implementation) and computes time derivatives and time integrated DOFs.
     *
     * @param i_dT time step.
     * @param i_stiffT transposed stiffness matrix (multiplied with inverse mass matrix).
     * @param i_star star matrices.
     * @param i_dofs DOFs.
     * @param o_scratch will be used as scratch memory.
     * @param o_der will be set to time derivatives.
     * @param o_tInt will be set to time integrated DOFs.
     **/
    template <typename TL_T_PRECISION>
    static void inline ckVanilla( TL_T_PRECISION       i_dT,
                                  TL_T_PRECISION const i_stiffT[TL_N_DIM][TL_N_ELEMENT_MODES][TL_N_ELEMENT_MODES],
                                  TL_T_PRECISION const i_star[TL_N_DIM][TL_N_QUANTITIES][TL_N_QUANTITIES],
                                  TL_T_PRECISION const i_dofs[TL_N_QUANTITIES][TL_N_ELEMENT_MODES][TL_N_CRUNS],
                                  TL_T_PRECISION       o_scratch[TL_N_QUANTITIES][TL_N_ELEMENT_MODES][TL_N_CRUNS],
                                  TL_T_PRECISION       o_der[TL_ORDER][TL_N_QUANTITIES][TL_N_ELEMENT_MODES][TL_N_CRUNS],
                                  TL_T_PRECISION       o_tInt[TL_N_QUANTITIES][TL_N_ELEMENT_MODES][TL_N_CRUNS] ) {
      // scalar for the time integration
      TL_T_PRECISION l_scalar = i_dT;

      // initialize zero-derivative, reset time integrated dofs
      for( int_qt l_qt = 0; l_qt < TL_N_QUANTITIES; l_qt++ ) {
        for( int_md l_md = 0; l_md < TL_N_ELEMENT_MODES; l_md++ ) {
          for( int_cfr l_cfr = 0; l_cfr < TL_N_CRUNS; l_cfr++ ) {
            o_der[0][l_qt][l_md][l_cfr] = i_dofs[l_qt][l_md][l_cfr];
            o_tInt[l_qt][l_md][l_cfr]   = l_scalar * i_dofs[l_qt][l_md][l_cfr];
          }
        }
      }

      // iterate over time derivatives
      for( unsigned int l_de = 1; l_de < TL_ORDER; l_de++ ) {
        // reset this derivative
        for( int_qt l_qt = 0; l_qt < TL_N_QUANTITIES; l_qt++ )
          for( int_md l_md = 0; l_md < TL_N_ELEMENT_MODES; l_md++ )
            for( int_cfr l_cr = 0; l_cr < TL_N_CRUNS; l_cr++ ) o_der[l_de][l_qt][l_md][l_cr] = 0;

        // compute the derivatives
        for( unsigned short l_di = 0; l_di < TL_N_DIM; l_di++ ) {
          // multiply with transposed stiffness matrices and inverse mass matrix
          linalg::Matrix::matMulB0FusedAC( TL_N_CRUNS,
                                           TL_N_QUANTITIES,
                                           TL_N_ELEMENT_MODES,
                                           TL_N_ELEMENT_MODES,
                                           o_der[l_de-1][0][0],
                                           i_stiffT[l_di][0],
                                           o_scratch[0][0] );

          // multiply with star matrices
          linalg::Matrix::matMulB1FusedBC( TL_N_CRUNS,
                                           TL_N_QUANTITIES,
                                           TL_N_ELEMENT_MODES,
                                           TL_N_QUANTITIES,
                                           i_star[l_di][0],
                                           o_scratch[0][0],
                                           o_der[l_de][0][0] );
        }

        // update scalar
        l_scalar *= -i_dT / (l_de+1);

        // update time integrated dofs
        for( int_qt l_qt = 0; l_qt < TL_N_QUANTITIES; l_qt++ ) {
          for( int_md l_md = 0; l_md < CE_N_ELEMENT_MODES_CK( TL_TYPE_ELEMENT, TL_ORDER, l_de ); l_md++ ) {
            for( int_cfr l_cr = 0; l_cr < TL_N_CRUNS; l_cr++ )
              o_tInt[l_qt][l_md][l_cr] += l_scalar * o_der[l_de][l_qt][l_md][l_cr];
          }
        }
      }
    }

    /**
     * Applies the Cauchy–Kowalevski procedure (fused LIBXSMM version) and computes time derivatives and time integrated DOFs.
     *
     * @param i_dT time step.
     * @param i_stiffT transposed stiffness matrix (multiplied with inverse mass matrix).
     * @param i_star star matrices.
     * @param i_dofs DOFs.
     * @param i_kernels generated LIBXSMM kernels.
     * @param o_scratch will be used as scratch memory.
     * @param o_der will be set to time derivatives.
     * @param o_tInt will be set to time integrated DOFs.
     **/
    template <typename TL_T_PRECISION, typename TL_T_XSMM_KERNEL >
    static void inline ckXsmmFused( TL_T_PRECISION                 i_dT,
                                    TL_T_PRECISION   const * const i_stiffT[(TL_ORDER>1) ? (TL_ORDER-1)*TL_N_DIM : TL_N_DIM],
                                    TL_T_PRECISION   const         i_star[TL_N_DIM][(TL_N_DIM==2) ? 10 : 24],
                                    TL_T_PRECISION   const         i_dofs[TL_N_QUANTITIES][TL_N_ELEMENT_MODES][TL_N_CRUNS],
                                    TL_T_XSMM_KERNEL const *       i_kernels,
                                    TL_T_PRECISION                 o_scratch[TL_N_QUANTITIES][TL_N_ELEMENT_MODES][TL_N_CRUNS],
                                    TL_T_PRECISION                 o_der[TL_ORDER][TL_N_QUANTITIES][TL_N_ELEMENT_MODES][TL_N_CRUNS],
                                    TL_T_PRECISION                 o_tInt[TL_N_QUANTITIES][TL_N_ELEMENT_MODES][TL_N_CRUNS] ) {
      // scalar for the time integration
      TL_T_PRECISION l_scalar = i_dT;

      // initialize zero-derivative, reset time integrated dofs
      for( int_qt l_qt = 0; l_qt < TL_N_QUANTITIES; l_qt++ ) {
        for( int_md l_md = 0; l_md < TL_N_ELEMENT_MODES; l_md++ ) {
#pragma omp simd
          for( int_cfr l_cfr = 0; l_cfr < TL_N_CRUNS; l_cfr++ ) {
            o_der[0][l_qt][l_md][l_cfr] = i_dofs[l_qt][l_md][l_cfr];
            o_tInt[l_qt][l_md][l_cfr]   = l_scalar * i_dofs[l_qt][l_md][l_cfr];
          }
        }
      }

      // iterate over time derivatives
      for( unsigned int l_de = 1; l_de < TL_ORDER; l_de++ ) {
        // reset this derivative
        for( int_qt l_qt = 0; l_qt < TL_N_QUANTITIES; l_qt++ )
          for( int_md l_md = 0; l_md < TL_N_ELEMENT_MODES; l_md++ )
#pragma omp simd
            for( int_cfr l_cr = 0; l_cr < TL_N_CRUNS; l_cr++ ) o_der[l_de][l_qt][l_md][l_cr] = 0;

        // compute the derivatives
        for( unsigned short l_di = 0; l_di < TL_N_DIM; l_di++ ) {
          // reset to zero for basis in non-hierarchical storage
          if( TL_TYPE_ELEMENT == TET4 || TL_TYPE_ELEMENT == TRIA3 ) {}
          else {
            zero( o_scratch );
          }

          // multiply with transposed stiffness matrices and inverse mass matrix
          i_kernels[(l_de-1)*(TL_N_DIM+1)+l_di]( o_der[l_de-1][0][0],
                                                 i_stiffT[(l_de-1)*TL_N_DIM+l_di],
                                                 o_scratch[0][0] );
          // multiply with star matrices
          i_kernels[(l_de-1)*(TL_N_DIM+1)+TL_N_DIM]( i_star[l_di],
                                                     o_scratch[0][0],
                                                     o_der[l_de][0][0] );
        }

        // update scalar
        l_scalar *= -i_dT / (l_de+1);

        // update time integrated dofs
        for( int_qt l_qt = 0; l_qt < TL_N_QUANTITIES; l_qt++ ) {
          for( int_md l_md = 0; l_md < CE_N_ELEMENT_MODES_CK( TL_TYPE_ELEMENT, TL_ORDER, l_de ); l_md++ ) {
#pragma omp simd
            for( int_cfr l_cr = 0; l_cr < TL_N_CRUNS; l_cr++ )
              o_tInt[l_qt][l_md][l_cr] += l_scalar * o_der[l_de][l_qt][l_md][l_cr];
          }
        }
      }
    }

    /**
     * Applies the Cauchy–Kowalevski procedure (single forward run LIBXSMM version) and computes time derivatives and time integrated DOFs.
     *
     * @param i_dT time step.
     * @param i_stiffT transposed stiffness matrix (multiplied with inverse mass matrix).
     * @param i_star star matrices.
     * @param i_dofs DOFs.
     * @param i_kernels LIBXSMM kernels.
     * @param o_scratch will be used as scratch memory.
     * @param o_der will be set to time derivatives.
     * @param o_tInt will be set to time integrated DOFs.
     **/
    template <typename TL_T_PRECISION, typename TL_T_XSMM_KERNEL >
    static void inline ckXsmmSingle( TL_T_PRECISION           i_dT,
                                     TL_T_PRECISION   const   i_stiffT[TL_N_DIM][TL_N_ELEMENT_MODES][TL_N_ELEMENT_MODES],
                                     TL_T_PRECISION   const   i_star[TL_N_DIM][TL_N_QUANTITIES][TL_N_QUANTITIES],
                                     TL_T_PRECISION   const   i_dofs[TL_N_QUANTITIES][TL_N_ELEMENT_MODES][1],
                                     TL_T_XSMM_KERNEL const * i_kernels,
                                     TL_T_PRECISION           o_scratch[TL_N_QUANTITIES][TL_N_ELEMENT_MODES][1],
                                     TL_T_PRECISION           o_der[TL_ORDER][TL_N_QUANTITIES][TL_N_ELEMENT_MODES][1],
                                     TL_T_PRECISION           o_tInt[TL_N_QUANTITIES][TL_N_ELEMENT_MODES][1] ) {
      // scalar for the time integration
      TL_T_PRECISION l_scalar = i_dT;

      // initialize zero-derivative, reset time integrated dofs
      for( int_qt l_qt = 0; l_qt < TL_N_QUANTITIES; l_qt++ ) {
#pragma omp simd
        for( int_md l_md = 0; l_md < TL_N_ELEMENT_MODES; l_md++ ) {
          o_der[0][l_qt][l_md][0] = i_dofs[l_qt][l_md][0];
          o_tInt[l_qt][l_md][0]   = l_scalar * i_dofs[l_qt][l_md][0];
        }
      }

      // iterate over time derivatives
      for( unsigned int l_de = 1; l_de < TL_ORDER; l_de++ ) {
        // reset this derivative
        for( int_qt l_qt = 0; l_qt < TL_N_QUANTITIES; l_qt++ )
#pragma omp simd
          for( int_md l_md = 0; l_md < TL_N_ELEMENT_MODES; l_md++ ) o_der[l_de][l_qt][l_md][0] = 0;

        // compute the derivatives
        for( unsigned short l_di = 0; l_di < TL_N_DIM; l_di++ ) {
          // multiply with transposed stiffness matrices and inverse mass matrix
          i_kernels[(l_de-1)*2]( i_stiffT[l_di][0],
                                 o_der[l_de-1][0][0],
                                 o_scratch[0][0] );
          // multiply with star matrices
          i_kernels[((l_de-1)*2)+1]( o_scratch[0][0],
                                     i_star[l_di][0],
                                     o_der[l_de][0][0] );
        }

        // update scalar
        l_scalar *= -i_dT / (l_de+1);

        // update time integrated dofs
        for( int_qt l_qt = 0; l_qt < TL_N_QUANTITIES; l_qt++ ) {
#pragma omp simd
          for( int_md l_md = 0; l_md < CE_N_ELEMENT_MODES_CK( TL_TYPE_ELEMENT, TL_ORDER, l_de ); l_md++ ) {
            o_tInt[l_qt][l_md][0] += l_scalar * o_der[l_de][l_qt][l_md][0];
          }
        }
      }
    }

    /**
     * Evaluates the time prediction, given by the time derivatives at the given points in time.
     * The points are relative to the time at which the time prediction was obtained.
     * Example:
     *   0    1.5  2.0  2.4 2.9 absolute time
     *   |-----|----x----x---x----------------->
     *         |   0.5  0.9 1.4 relative time (expected as input)
     *       time
     *    prediction
     *
     * @param i_pts relative pts in time.
     * @param i_der time prediction given through the time derivatives.
     * @param o_preDofs will be set to the predicted DOFs at the points in time.
     **/
    template <typename TL_T_PRECISION >
    static void inline evalTimePrediction( TL_T_PRECISION const i_pts[TL_N_POINTS_TIME],
                                           TL_T_PRECISION const i_der[TL_ORDER][TL_N_QUANTITIES][TL_N_ELEMENT_MODES][TL_N_CRUNS],
                                           TL_T_PRECISION       o_preDofs[TL_N_POINTS_TIME][TL_N_QUANTITIES][TL_N_ELEMENT_MODES][TL_N_CRUNS] ) {
      for( unsigned short l_pt = 0; l_pt < TL_N_POINTS_TIME; l_pt++ ) {
        // init dofs
        for( int_qt l_qt = 0; l_qt < TL_N_QUANTITIES; l_qt++ )
          for( int_md l_md = 0; l_md < TL_N_ELEMENT_MODES; l_md++ )
            for( int_cfr l_ru = 0; l_ru < TL_N_CRUNS; l_ru++ ) o_preDofs[l_pt][l_qt][l_md][l_ru] = i_der[0][l_qt][l_md][l_ru];

        // evaluate time derivatives at given point in time
        real_base l_scalar = 1.0;

        // iterate over derivatives
        for( unsigned short l_de = 1; l_de < TL_ORDER; l_de++ ) {
          // update scalar
          l_scalar *= -i_pts[l_pt] / l_de;

          for( int_qt l_qt = 0; l_qt < TL_N_QUANTITIES; l_qt++ ) {
            for( int_md l_md = 0; l_md < CE_N_ELEMENT_MODES_CK( TL_TYPE_ELEMENT, TL_ORDER, l_de ); l_md++ ) {
              for( int_cfr l_ru = 0; l_ru < TL_N_CRUNS; l_ru++ )
                o_preDofs[l_pt][l_qt][l_md][l_ru] += l_scalar * i_der[l_de][l_qt][l_md][l_ru];
            }
          }
        }
      }
    }
};

#endif
