/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *         Alexander Heinecke (alexander.heinecke AT intel.com)
 *
 * @section LICENSE
 * Copyright (c) 2019, Alexander Breuer
 * Copyright (c) 2016-2018, Regents of the University of California
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
 * Quadrature-free ADER-DG volume integration for seismic wave propagation using vanilla kernels.
 **/
#ifndef EDGE_SEISMIC_KERNELS_VOL_INT_VANILLA_HPP
#define EDGE_SEISMIC_KERNELS_VOL_INT_VANILLA_HPP

#include "VolInt.hpp"
#include "dg/Basis.h"
#include "data/MmVanilla.hpp"

namespace edge {
  namespace elastic {
    namespace kernels { 
      template< typename       TL_T_REAL,
                t_entityType   TL_T_EL,
                unsigned short TL_O_SP,
                unsigned short TL_N_CRS >
      class VolIntVanilla;
    }
  }
}

/**
 * Quadrature-free ADER-DG volume integration for seismic wave propagation using vanilla kernels.
 *
 * @paramt TL_T_REAL floating point precision.
 * @paramt TL_T_EL element type.
 * @paramt TL_O_SP spatial order.
 * @paramt TL_N_CRS number of fused simulations.
 **/
template< typename       TL_T_REAL,
          t_entityType   TL_T_EL,
          unsigned short TL_O_SP,
          unsigned short TL_N_CRS >
class edge::elastic::kernels::VolIntVanilla: edge::elastic::kernels::VolInt < TL_T_REAL,
                                                                              TL_T_EL,
                                                                              TL_O_SP,
                                                                              TL_N_CRS > {
  private:
    //! dimension of the element
    static unsigned short const TL_N_DIS = C_ENT[TL_T_EL].N_DIM;

    //! number of element modes
    static unsigned short const TL_N_MDS = CE_N_ELEMENT_MODES( TL_T_EL, TL_O_SP );

    //! number of elastic quantities
    static unsigned short const TL_N_QTS_E = (TL_N_DIS == 2) ? 5 : 9;

    //! matrix kernels
    edge::data::MmVanilla< TL_T_REAL > m_mm;

    //! pointers to the (possibly recursive) stiffness matrices
    TL_T_REAL *m_stiff[TL_N_DIS];

    /**
     * Generates the matrix kernels for the stiffness matrices and star matrices.
     **/
    void generateKernels() {
      // multiplication with stiffness matrix and star matrix
      m_mm.add( 0,                                            // group
                TL_N_QTS_E,                                   // m
                TL_N_MDS,                                     // n
                CE_N_ELEMENT_MODES_CK( TL_T_EL, TL_O_SP, 1 ), // k
                TL_N_MDS,                                     // ldA
                TL_N_MDS,                                     // ldB
                TL_N_MDS,                                     // ldC
                static_cast<TL_T_REAL>(1.0),                  // alpha
                static_cast<TL_T_REAL>(0.0),                  // beta
                true,                                         // fused AC
                false,                                        // fused BC
                TL_N_CRS );

      m_mm.add( 0,                           // group
                TL_N_QTS_E,                  // m
                TL_N_MDS,                    // n
                TL_N_QTS_E,                  // k
                TL_N_QTS_E,                  // ldA
                TL_N_MDS,                    // ldB
                TL_N_MDS,                    // ldC
                static_cast<TL_T_REAL>(1.0), // alpha
                static_cast<TL_T_REAL>(1.0), // beta
                false,                       // fused AC
                true,                        // fused BC
                TL_N_CRS );
    }

  public:
    /**
     * Constructor of the vanilla volume integration.
     *
     * @param io_dynMem dynamic memory allocations.
     **/
    VolIntVanilla( data::Dynamic & io_dynMem ) {
      // formulation of the basis in terms of the reference element
      dg::Basis l_basis( TL_T_EL,
                         TL_O_SP );

      // get stiffness matrices
      TL_T_REAL l_stiff[TL_N_DIS][TL_N_MDS][TL_N_MDS];
      l_basis.getStiffMm1Dense( TL_N_MDS,
                                l_stiff[0][0],
                                false );

      // store stiffness matrices dense
      this->storeStiffDense( l_stiff,
                             io_dynMem,
                             m_stiff );


      generateKernels();
    }

    /**
     * Volume contribution using vanilla matrix-matrix multiplication kernels.
     *
     * @param i_jac jacobians.
     * @param i_tDofs time integerated DG-DOFs.
     * @param io_dofs will be updated with local contribution of the element to the volume integral.
     * @param o_scratch will be used as scratch space for the computations.
     **/
    void apply( TL_T_REAL const i_jac[TL_N_DIS][TL_N_QTS_E][TL_N_QTS_E],
                TL_T_REAL const i_tDofs[TL_N_QTS_E][TL_N_MDS][TL_N_CRS],
                TL_T_REAL       io_dofs[TL_N_QTS_E][TL_N_MDS][TL_N_CRS],
                TL_T_REAL       o_scratch[TL_N_QTS_E][TL_N_MDS][TL_N_CRS] ) const {
      // iterate over dimensions
      for( unsigned short l_di = 0; l_di < TL_N_DIS; l_di++ ) {
        // multiply with stiffness and inverse mass matrix
        m_mm.m_kernels[0][0]( i_tDofs[0][0],
                              m_stiff[l_di],
                              o_scratch[0][0] );

        // multiply with star matrix
        m_mm.m_kernels[0][1]( i_jac[l_di][0],
                              o_scratch[0][0],
                              io_dofs[0][0] );
      }
    }
};

#endif