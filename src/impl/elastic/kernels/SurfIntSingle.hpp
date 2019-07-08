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
 * Optimized quadrature-free ADER-DG surface integration for single seismic forward simulations.
 **/
#ifndef EDGE_SEISMIC_KERNELS_SURF_INT_SINGLE_HPP
#define EDGE_SEISMIC_KERNELS_SURF_INT_SINGLE_HPP

#include "SurfInt.hpp"
#include "dg/Basis.h"
#include "data/MmXsmmSingle.hpp"

namespace edge {
  namespace elastic {
    namespace kernels { 
      template< typename       TL_T_REAL,
                t_entityType   TL_T_EL,
                unsigned short TL_O_SP >
      class SurfIntSingle;
    }
  }
}

/**
 * Optimized quadrature-free ADER-DG surface integration for single seismic forward simulations.
 *
 * @paramt TL_T_REAL floating point precision.
 * @paramt TL_T_EL element type.
 * @paramt TL_O_SP spatial order.
 **/
template< typename       TL_T_REAL,
          t_entityType   TL_T_EL,
          unsigned short TL_O_SP >
class edge::elastic::kernels::SurfIntSingle: edge::elastic::kernels::SurfInt < TL_T_REAL,
                                                                               TL_T_EL,
                                                                               TL_O_SP,
                                                                               1 > {
  private:
    //! number of dimensions
    static unsigned short const TL_N_DIS = C_ENT[TL_T_EL].N_DIM;

    //! number of faces
    static unsigned short const TL_N_FAS = C_ENT[TL_T_EL].N_FACES;

    //! number of DG face modes
    static unsigned short const TL_N_MDS_FA = CE_N_ELEMENT_MODES( C_ENT[TL_T_EL].TYPE_FACES, TL_O_SP );

    //! number of DG element modes
    static unsigned short const TL_N_MDS_EL = CE_N_ELEMENT_MODES( TL_T_EL, TL_O_SP );

    //! number of neigboring contribution flux matrices
    static unsigned short const TL_N_FMNS = CE_N_FLUXN_MATRICES( TL_T_EL );

    //! number of elastic quantities
    static unsigned short const TL_N_QTS_E = (TL_N_DIS == 2) ? 5 : 9;

    //! pointers to the local flux matrices
    TL_T_REAL *m_fIntLN[TL_N_FAS+TL_N_FMNS];

    //! pointers to the transposed flux matrices
    TL_T_REAL *m_fIntT[TL_N_FAS];

    //! matrix kernels
    edge::data::MmXsmmSingle< TL_T_REAL > m_mm;

    /**
     * Generates the matrix kernels for the flux matrices and flux solvers.
     **/
    void generateKernels() {
      // add first flux matrix
      m_mm.add( 0,                           // group
                TL_N_MDS_FA,                 // m
                TL_N_QTS_E,                  // n
                TL_N_MDS_EL,                 // k
                TL_N_MDS_FA,                 // ldA
                TL_N_MDS_EL,                 // ldB
                TL_N_MDS_FA,                 // ldC
                static_cast<real_base>(1.0), // alpha
                static_cast<real_base>(0.0), // beta
                LIBXSMM_GEMM_PREFETCH_AL2BL2_VIA_C_AHEAD );

      // add flux solver
      m_mm.add( 0,                           // group
                TL_N_MDS_FA,                 // m
                TL_N_QTS_E,                  // n
                TL_N_QTS_E,                  // k
                TL_N_MDS_FA,                 // ldA
                TL_N_QTS_E,                  // ldB
                TL_N_MDS_FA,                 // ldC
                static_cast<real_base>(1.0), // alpha
                static_cast<real_base>(0.0), // beta
                LIBXSMM_GEMM_PREFETCH_NONE );

      // add second flux matrix
      m_mm.add( 0,                           // group
                TL_N_MDS_EL,                 // m
                TL_N_QTS_E,                  // n
                TL_N_MDS_FA,                 // k
                TL_N_MDS_EL,                 // ldA
                TL_N_MDS_FA,                 // ldB
                TL_N_MDS_EL,                 // ldC
                static_cast<real_base>(1.0), // alpha
                static_cast<real_base>(1.0), // beta
                LIBXSMM_GEMM_PREFETCH_AL2BL2_VIA_C_AHEAD );
    }

  public:
    /**
     * Constructor of the surface integrations for single forward simulations.
     *
     * @param io_dynMem dynamic memory allocations.
     **/
    SurfIntSingle( data::Dynamic & io_dynMem ) {
      // formulation of the basis in terms of the reference element
      dg::Basis l_basis( TL_T_EL,
                         TL_O_SP );

      // get flux matrices
      TL_T_REAL l_fIntL[TL_N_FAS][TL_N_MDS_EL][TL_N_MDS_FA];
      TL_T_REAL l_fIntN[TL_N_FMNS][TL_N_MDS_EL][TL_N_MDS_FA];
      TL_T_REAL l_fIntT[TL_N_FAS][TL_N_MDS_FA][TL_N_MDS_EL];
      l_basis.getFluxDense( l_fIntL[0][0],
                            l_fIntN[0][0],
                            l_fIntT[0][0] );

      // store flux matrices dense
      this->storeFluxDense( l_fIntL,
                            l_fIntN,
                            l_fIntT,
                            io_dynMem,
                            m_fIntLN,
                            m_fIntT );

      // generate kernels
      generateKernels();
    }

    /**
     * Element local contribution for single forward simulations.
     *
     * @param i_fSol flux solvers.
     * @param i_tDofs time integrated DG-DOFs.
     * @param io_dofs will be updated with local contribution of the element to the surface integral.
     * @param o_scratch will be used as scratch space for the computations.
     * @param i_dofsP DOFs for prefetching (not used).
     * @param i_tDofsP time integrated DOFs for prefetching (not used).
     **/
    void local( TL_T_REAL const i_fSol[TL_N_FAS][TL_N_QTS_E][TL_N_QTS_E],
                TL_T_REAL const i_tDofs[TL_N_QTS_E][TL_N_MDS_EL][1],
                TL_T_REAL       io_dofs[TL_N_QTS_E][TL_N_MDS_EL][1],
                TL_T_REAL       o_scratch[2][TL_N_QTS_E][TL_N_MDS_FA][1],
                TL_T_REAL const i_dofsP[TL_N_QTS_E][TL_N_MDS_EL][1] = nullptr,
                TL_T_REAL const i_tDofsP[TL_N_QTS_E][TL_N_MDS_EL][1] = nullptr ) const {
      // iterate over faces
      for( unsigned short l_fa = 0; l_fa < TL_N_FAS; l_fa++ ) {
        // multiply with first face integration matrix
        m_mm.m_kernels[0][0]( m_fIntLN[l_fa],
                              i_tDofs[0][0],
                              o_scratch[0][0][0],
                              nullptr,
                              i_dofsP[0][0],
                              nullptr );

        // multiply with flux solver
        m_mm.m_kernels[0][1]( o_scratch[0][0][0],
                              i_fSol[l_fa][0],
                              o_scratch[1][0][0] );

        // multiply with second face integration matrix
        m_mm.m_kernels[0][2]( m_fIntT[l_fa],
                              o_scratch[1][0][0],
                              io_dofs[0][0],
                              nullptr,
                              i_tDofsP[0][0],
                              nullptr );
      }
    }

    /**
     * Neighboring contribution of a single adjacent element for single forward simulations.
     *
     * @param i_fa local face.
     * @param i_vId id of the vertex, matching the element's vertex 0, from the perspective of the adjacent element w.r.t. to the reference element.
     * @param i_fId id of the face from the perspective of the adjacent element w.r.t. to the reference element.
     * @param i_fSol flux solvers.
     * @param i_tDofs time integrated DG-DOFs.
     * @param io_dofs will be updated with the contribution of the adjacent element to the surface integral.
     * @param o_scratch will be used as scratch space for the computations.
     * @param i_pre DOFs or tDOFs for prefetching.
     **/
    void neigh( unsigned short       i_fa,
                unsigned short       i_vId,
                unsigned short       i_fId,
                TL_T_REAL      const i_fSol[TL_N_QTS_E][TL_N_QTS_E],
                TL_T_REAL      const i_tDofs[TL_N_QTS_E][TL_N_MDS_EL][1],
                TL_T_REAL            io_dofs[TL_N_QTS_E][TL_N_MDS_EL][1],
                TL_T_REAL            o_scratch[2][TL_N_QTS_E][TL_N_MDS_FA][1],
                TL_T_REAL      const i_pre[TL_N_QTS_E][TL_N_MDS_EL][1] = nullptr ) const {
      // derive the id of the neighboring flux matrix
      unsigned short l_fMatId = std::numeric_limits< unsigned short >::max();
      if( i_vId != std::numeric_limits< unsigned short >::max() ) {
        l_fMatId = TL_N_FAS + this->fMatId( i_vId, i_fId );
      }
      else {
        l_fMatId = i_fa;
      }

      // multiply with first face integration matrix
      m_mm.m_kernels[0][0]( m_fIntLN[l_fMatId],
                            i_tDofs[0][0],
                            o_scratch[0][0][0],
                            nullptr,
                            i_pre[0][0],
                            nullptr );

      // multiply with flux solver
      m_mm.m_kernels[0][1]( o_scratch[0][0][0],
                            i_fSol[0],
                            o_scratch[1][0][0] );

      // multiply with second face integration matrix
      m_mm.m_kernels[0][2]( m_fIntT[i_fa],
                            o_scratch[1][0][0],
                            io_dofs[0][0],
                            nullptr,
                            i_pre[0][0],
                            nullptr );
    }
};

#endif