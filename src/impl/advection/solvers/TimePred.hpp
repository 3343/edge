/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2016-2018, Regents of the University of California
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
 * ADER time prediction for the advection equation.
 **/
#ifndef EDGE_ADVECTION_TIME_PRED_HPP
#define EDGE_ADVECTION_TIME_PRED_HPP

#include "constants.hpp"
#include "linalg/Matrix.h"

namespace edge {
  namespace advection {
    namespace solvers {
      template< t_entityType   TL_T_EL,
                unsigned short TL_O_SP,
                unsigned short TL_O_TI,
                unsigned short TL_N_CRS >
      class TimePred;
    }
  }
}

/**
 * ADER time predictions for the advection equation.
 *
 * @param TL_T_EL element typ.
 * @param TL_O_SP order in space.
 * @param TL_O_TI order in time.
 * @param TL_N_CRS number of fused simulations
 **/
template< t_entityType   TL_T_EL,
          unsigned short TL_O_SP,
          unsigned short TL_O_TI,
          unsigned short TL_N_CRS >
class edge::advection::solvers::TimePred {
  private:
    //! number of dimensions
    static unsigned short const TL_N_DIM = C_ENT[TL_T_EL].N_DIM;

    //! number of element modes.
    static unsigned short const TL_N_MDS = CE_N_ELEMENT_MODES( TL_T_EL, TL_O_SP );

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
     *
     * @paramt TL_T_REAL float point precision.
     **/
    template< typename TL_T_REAL >
    static void inline ckVanilla( double          i_dT,
                                  TL_T_REAL const i_stiffT[TL_N_DIM][TL_N_MDS][TL_N_MDS],
                                  TL_T_REAL const i_star[TL_N_DIM],
                                  TL_T_REAL const i_dofs[TL_N_MDS][TL_N_CRS],
                                  TL_T_REAL       o_scratch[TL_N_MDS][TL_N_CRS],
                                  TL_T_REAL       o_der[TL_O_TI][TL_N_MDS][TL_N_CRS],
                                  TL_T_REAL       o_tInt[TL_N_MDS][TL_N_CRS] ) {
      // scalar coefficients in taylor expansion
      TL_T_REAL l_scalar = i_dT;

      // initialize zero derivatives, reset time integrated dofs
      for( unsigned short l_md = 0; l_md < TL_N_MDS; l_md++ ) {
        for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {
          o_der[0][l_md][l_cr] = i_dofs[l_md][l_cr];
          o_tInt[l_md][l_cr]   = l_scalar * i_dofs[l_md][l_cr];
        }
      }

      // iterate over time derivatives
      for( unsigned int l_de = 1; l_de < TL_O_TI; l_de++ ) {
        // reset derivative
        for( unsigned short l_md = 0; l_md < TL_N_MDS; l_md++ )
          for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) o_der[l_de][l_md][l_cr] = 0;

        for( unsigned int l_di = 0; l_di < TL_N_DIM; l_di++ ) {
          // multiply with transposed stiffness matrices and inverse mass matrix
          linalg::Matrix::matMulFusedAC( TL_N_CRS,                    // #fused
                                         1,                           // m
                                         TL_N_MDS,                    // n
                                         TL_N_MDS,                    // k
                                         TL_N_MDS,                    // ldA
                                         TL_N_MDS,                    // ldB
                                         TL_N_MDS,                    // ldC
                                         static_cast<real_base>(0.0), // beta
                                         o_der[l_de-1][0],
                                         i_stiffT[l_di][0],
                                         o_scratch[0] );

          // multiply with star "matrices"
          linalg::Matrix::matMulFusedBC( TL_N_CRS,                    // #fused
                                         1,                           // m
                                         TL_N_MDS,                    // n
                                         1,                           // k
                                         1,                           // ldA
                                         TL_N_MDS,                    // ldB
                                         TL_N_MDS,                    // ldC
                                         static_cast<real_base>(1.0), // beta
                                         i_star+l_di,
                                         o_scratch[0],
                                         o_der[l_de][0] );
        }

        // update scalar
        l_scalar *= -i_dT / (l_de+1);

        // update time integrated dofs
        for( unsigned short l_md = 0; l_md < TL_N_MDS; l_md++ )
          for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) o_tInt[l_md][l_cr] += l_scalar * o_der[l_de][l_md][l_cr];
      }

    }
};

#endif
