/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2019, Alexander Breuer
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
 * Volume integration for the advection equation.
 **/
#ifndef EDGE_ADVECTION_KERNELS_VOL_INT_HPP
#define EDGE_ADVECTION_KERNELS_VOL_INT_HPP

#include "constants.hpp"
#include "linalg/Matrix.h"
#include "dg/VolInt.hpp"

namespace edge {
  namespace advection {
    namespace kernels {
      template< typename       TL_T_REAL,
                t_entityType   TL_T_EL,
                unsigned short TL_O_SP,
                unsigned short TL_N_CRS >
      class VolInt;
    }
  }
}

/**
 * Quadrature-free ADER-DG volume integration for the advection equation.
 *
 * @paramt TL_T_REAL real type.
 * @paramt TL_T_EL element type.
 * @paramt TL_O_SP spatial order.
 * @paramt TL_N_CRS number of fused simulations.
 **/
template< typename       TL_T_REAL,
          t_entityType   TL_T_EL,
          unsigned short TL_O_SP,
          unsigned short TL_N_CRS >
class edge::advection::kernels::VolInt {
  private:
    //! number of dimensions
    static unsigned short const TL_N_DIS = C_ENT[TL_T_EL].N_DIM;

    //! number of element modes
    static unsigned short const TL_N_MDS = CE_N_ELEMENT_MODES( TL_T_EL, TL_O_SP );

    //! pointers to the stiffness matrices
    TL_T_REAL *m_stiff[TL_N_DIS];

  public:
    /**
     * Constructor of the volume integration.
     *
     * @param io_dynMem dynamic memory allocations.
     **/
    VolInt( data::Dynamic & io_dynMem ) {
      // store stiffness matrices
      dg::VolInt< TL_T_EL,
                  TL_O_SP >::storeStiffDense( io_dynMem,
                                              m_stiff );
    }

    /**
     * Applies the volume contribution.
     *
     * @param i_star star matrices.
     * @param i_tDofs time integrated degrees of freedom.
     * @param io_dofs will be updated with the contribution of the volume integral.
     **/
    void apply( TL_T_REAL const i_star[TL_N_DIS],
                TL_T_REAL const i_tDofs[1][TL_N_MDS][TL_N_CRS],
                TL_T_REAL       io_dofs[1][TL_N_MDS][TL_N_CRS] ) {
      // temporary product for two-way mult
      TL_T_REAL l_tmpProd[TL_N_MDS][TL_N_CRS];

      // iterate over dimensions
      for( unsigned short l_di = 0; l_di < TL_N_DIS; l_di++ ) {
        // multiply with stiffness and inverse mass matrix
        linalg::Matrix::matMulFusedAC( TL_N_CRS,       // #fused
                                       1,              // m
                                       TL_N_MDS,       // n
                                       TL_N_MDS,       // k
                                       TL_N_MDS,       // ldA
                                       TL_N_MDS,       // ldB
                                       TL_N_MDS,       // ldC
                                       TL_T_REAL(1.0), // alpha
                                       TL_T_REAL(0.0), // beta
                                       i_tDofs[0][0],  // A
                                       m_stiff[l_di],  // B
                                       l_tmpProd[0] ); // C

        // multiply with star "matrix"
        linalg::Matrix::matMulFusedBC( TL_N_CRS,       // #fused
                                       1,              // m
                                       TL_N_MDS,       // n
                                       1,              // k
                                       TL_N_MDS,       // ldA
                                       TL_N_MDS,       // ldB
                                       TL_N_MDS,       // ldC
                                       TL_T_REAL(1.0), // alpha
                                       TL_T_REAL(1.0), // beta
                                       i_star+l_di,    // A
                                       l_tmpProd[0],   // B
                                       io_dofs[0][0]); // C
      }
    }
};

#endif