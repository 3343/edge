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
 * Surface integration for the advection equation.
 **/
#ifndef EDGE_ADVECTION_KERNELS_SURF_INT_HPP
#define EDGE_ADVECTION_KERNELS_SURF_INT_HPP

#include "constants.hpp"
#include "linalg/Matrix.h"
#include "dg/SurfInt.hpp"

namespace edge {
  namespace advection {
    namespace kernels {
      template< typename       TL_T_REAL,
                t_entityType   TL_T_EL,
                unsigned short TL_O_SP,
                unsigned short TL_N_CRS >
      class SurfInt;
    }
  }
}

/**
 * Quadrature-free ADER-DG surface integration for the advection equation.
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
class edge::advection::kernels::SurfInt {
  private:
    //! number of dimensions
    static unsigned short const TL_N_DIS = C_ENT[TL_T_EL].N_DIM;

    //! number of faces
    static unsigned short const TL_N_FAS = C_ENT[TL_T_EL].N_FACES;

    //! number of face modes
    static unsigned short const TL_N_MDS_FA = CE_N_ELEMENT_MODES( C_ENT[TL_T_EL].TYPE_FACES, TL_O_SP );

    //! number of element modes
    static unsigned short const TL_N_MDS_EL = CE_N_ELEMENT_MODES( TL_T_EL, TL_O_SP );

    //! number of neigboring contribution flux matrices
    static unsigned short const TL_N_FMNS = CE_N_FLUXN_MATRICES( TL_T_EL );

    //! number of element modes
    static unsigned short const TL_N_MDS = CE_N_ELEMENT_MODES( TL_T_EL, TL_O_SP );

    //! pointers to the local flux matrices
    TL_T_REAL *m_fIntLN[TL_N_FAS+TL_N_FMNS];

    //! pointers to the transposed flux matrices
    TL_T_REAL *m_fIntT[TL_N_FAS];

  public:
    /**
     * Constructor of the surface integration.
     *
     * @param io_dynMem dynamic memory allocations.
     **/
    SurfInt( data::Dynamic & io_dynMem ) {
      // store stiffness matrices
      dg::SurfInt< TL_T_EL,
                   TL_O_SP >::storeFluxDense( io_dynMem,
                                              m_fIntLN,
                                              m_fIntT );
    }

    /**
     * Element-local contribution.
     *
     * @param i_fs flux solvers.
     * @param i_tDofs time integrated degrees of freedom.
     * @param io_dofs will be updated with the contribution of the local surface integral.
     **/
    void local( TL_T_REAL const i_fs[TL_N_FAS],
                TL_T_REAL const i_tDofs[1][TL_N_MDS][TL_N_CRS],
                TL_T_REAL       io_dofs[1][TL_N_MDS][TL_N_CRS] ) {
      for( unsigned short l_fa = 0; l_fa < TL_N_FAS; l_fa++ ) {
        // scratch space for three-way product
        TL_T_REAL l_scratch[2][TL_N_MDS_FA][TL_N_CRS];

        linalg::Matrix::matMulFusedAC( TL_N_CRS,          // #fused
                                       1,                 // m
                                       TL_N_MDS_FA,       // n
                                       TL_N_MDS_EL,       // k
                                       TL_N_MDS_EL,       // ldA
                                       TL_N_MDS_FA,       // ldB
                                       TL_N_MDS_FA,       // ldC
                                       TL_T_REAL(1.0),    // alpha
                                       TL_T_REAL(0.0),    // beta
                                       i_tDofs[0][0],     // A
                                       m_fIntLN[l_fa],    // B
                                       l_scratch[0][0] ); // C

        linalg::Matrix::matMulFusedBC( TL_N_CRS,          // #fused
                                       1,                 // m
                                       TL_N_MDS_FA,       // n
                                       1,                 // k
                                       1,                 // ldA
                                       TL_N_MDS_FA,       // ldB
                                       TL_N_MDS_FA,       // ldC
                                       TL_T_REAL(1.0),    // alpha
                                       TL_T_REAL(0.0),    // beta
                                       i_fs+l_fa,         // A
                                       l_scratch[0][0],   // B
                                       l_scratch[1][0] ); // C

        linalg::Matrix::matMulFusedAC( TL_N_CRS,        // #fused
                                       1,               // m
                                       TL_N_MDS_EL,     // n
                                       TL_N_MDS_FA,     // k
                                       TL_N_MDS_FA,     // ldA
                                       TL_N_MDS_EL,     // ldB
                                       TL_N_MDS_EL,     // ldC
                                       TL_T_REAL(1.0),  // alpha
                                       TL_T_REAL(1.0),  // beta
                                       l_scratch[1][0], // A
                                       m_fIntT[l_fa],   // B
                                       io_dofs[0][0] ); // C
      }
    }

    /**
     * Contribution of a neighboring element to the surface integral.
     *
     * @param i_fa local face.
     * @param i_vId id of the vertex, matching the element's vertex 0, from the perspective of the adjacent element w.r.t. to the reference element.
     * @param i_fId id of the face from the perspective of the adjacent element w.r.t. to the reference element.
     * @param i_fs flux solver.
     * @param i_tDofs time integrated degrees of freedom.
     * @param io_dofs will be update with the neighboring element's contribution.
     **/
    void neigh( unsigned short       i_fa,
                unsigned short       i_vId,
                unsigned short       i_fId,
                TL_T_REAL            i_fs,
                TL_T_REAL      const i_tDofs[1][TL_N_MDS_EL][TL_N_CRS],
                TL_T_REAL            io_dofs[1][TL_N_MDS_EL][TL_N_CRS] ) {
      // derive the id of the neighboring flux matrix
      unsigned short l_fMatId = std::numeric_limits< unsigned short >::max();
      if( i_vId != std::numeric_limits< unsigned short >::max() ) {
        l_fMatId = TL_N_FAS + i_vId * TL_N_FAS;
        l_fMatId += i_fId;
      }
      else {
        l_fMatId = i_fa;
      }

      // scratch space for three-way product
      TL_T_REAL l_scratch[2][TL_N_MDS_FA][TL_N_CRS];

      linalg::Matrix::matMulFusedAC( TL_N_CRS,           // #fused
                                     1,                  // m
                                     TL_N_MDS_FA,        // n
                                     TL_N_MDS_EL,        // k
                                     TL_N_MDS_EL,        // ldA
                                     TL_N_MDS_FA,        // ldB
                                     TL_N_MDS_FA,        // ldC
                                     TL_T_REAL(1.0),     // alpha
                                     TL_T_REAL(0.0),     // beta
                                     i_tDofs[0][0],      // A
                                     m_fIntLN[l_fMatId], // B
                                     l_scratch[0][0] );  // C

      linalg::Matrix::matMulFusedBC( TL_N_CRS,          // #fused
                                     1,                 // m
                                     TL_N_MDS_FA,       // n
                                     1,                 // k
                                     1,                 // ldA
                                     TL_N_MDS_FA,       // ldB
                                     TL_N_MDS_FA,       // ldC
                                     TL_T_REAL(1.0),    // alpha
                                     TL_T_REAL(0.0),    // beta
                                     &i_fs,             // A
                                     l_scratch[0][0],   // B
                                     l_scratch[1][0] ); // C

      linalg::Matrix::matMulFusedAC( TL_N_CRS,        // #fused
                                     1,               // m
                                     TL_N_MDS_EL,     // n
                                     TL_N_MDS_FA,     // k
                                     TL_N_MDS_FA,     // ldA
                                     TL_N_MDS_EL,     // ldB
                                     TL_N_MDS_EL,     // ldC
                                     TL_T_REAL(1.0),  // alpha
                                     TL_T_REAL(1.0),  // beta
                                     l_scratch[1][0], // A
                                     m_fIntT[i_fa],   // B
                                     io_dofs[0][0] ); // C
    }
};

#endif