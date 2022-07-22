/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *         Alexander Heinecke (alexander.heinecke AT intel.com)
 *
 * @section LICENSE
 * Copyright (c) 2020, Friedrich Schiller University Jena
 * Copyright (c) 2019-2020, Alexander Breuer
 * Copyright (c) 2016-2018, Regents of the University of California
 * Copyright (c) 2016-2020, Intel Corporation
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
 * Optimized quadrature-free ADER-DG surface integration for fused seismic forward simulations.
 **/
#ifndef EDGE_SEISMIC_KERNELS_SURF_INT_FUSED_HPP
#define EDGE_SEISMIC_KERNELS_SURF_INT_FUSED_HPP

#include "SurfInt.hpp"
#include "dg/Basis.h"
#include "FakeMats.hpp"

#include "data/MmXsmmFused.hpp"
#include "data/UnaryXsmm.hpp"
#include "data/BinaryXsmm.hpp"
#include "data/TernaryXsmm.hpp"

#define ELTWISE_TPP

#ifdef PP_MMKERNEL_PERF
#include "parallel/global.h"
#endif 

namespace edge {
  namespace seismic {
    namespace kernels { 
      template< typename       TL_T_REAL,
                unsigned short TL_N_RMS,
                t_entityType   TL_T_EL,
                unsigned short TL_O_SP,
                unsigned short TL_N_CRS >
      class SurfIntFused;
    }
  }
}

/**
 * Quadrature-free ADER-DG surface integration for fused seismic forward simulations.
 *
 * @paramt TL_T_REAL floating point precision.
 * @paramt TL_N_RMS number of relaxation mechanisms.
 * @paramt TL_T_EL element type.
 * @paramt TL_O_SP spatial order.
 * @paramt TL_N_CRS number of fused simulations.
 **/
template< typename       TL_T_REAL,
          unsigned short TL_N_RMS,
          t_entityType   TL_T_EL,
          unsigned short TL_O_SP,
          unsigned short TL_N_CRS >
class edge::seismic::kernels::SurfIntFused: public edge::seismic::kernels::SurfInt< TL_T_REAL,
                                                                                    TL_N_RMS,
                                                                                    TL_T_EL,
                                                                                    TL_O_SP,
                                                                                    TL_N_CRS > {
  private:
    //! number of dimensions
    static unsigned short const TL_N_DIS = C_ENT[TL_T_EL].N_DIM;

    //! number of faces
    static unsigned short const TL_N_FAS = C_ENT[TL_T_EL].N_FACES;

    //! half the number of faces
    static unsigned short const TL_N_FAS_DIV2 =  TL_N_FAS / 2;

    //! number of DG face modes
    static unsigned short const TL_N_MDS_FA = CE_N_ELEMENT_MODES( C_ENT[TL_T_EL].TYPE_FACES, TL_O_SP );

    //! number of DG element modes
    static unsigned short const TL_N_MDS_EL = CE_N_ELEMENT_MODES( TL_T_EL, TL_O_SP );

    //! number of neigboring contribution flux matrices
    static unsigned short const TL_N_FMNS = CE_N_FLUXN_MATRICES( TL_T_EL );

    //! number of elastic quantities
    static unsigned short const TL_N_QTS_E = CE_N_QTS_E( TL_N_DIS );

    //! number of quantities per relaxation mechanism
    static unsigned short const TL_N_QTS_M = CE_N_QTS_M( TL_N_DIS );

    //! number of entries in the elastic flux solvers
    static unsigned short const TL_N_ENS_FS_E = CE_N_ENS_FS_E_DE( TL_N_DIS );

    //! number of entries in the anelastic flux solvers
    static unsigned short const TL_N_ENS_FS_A = CE_N_ENS_FS_A_DE( TL_N_DIS );

    //! pointers to the local flux matrices
    TL_T_REAL *m_fIntLN[TL_N_FAS+TL_N_FMNS] = {};

    //! pointers to the transposed flux matrices
    TL_T_REAL *m_fIntT[TL_N_FAS] = {};

    //! matrix kernels
    edge::data::MmXsmmFused< TL_T_REAL > m_mm;

    //! unary kernels
    edge::data::UnaryXsmm< TL_T_REAL > u_unary;

    //! binary kernels
    edge::data::BinaryXsmm< TL_T_REAL > b_binary;

    //! ternary kernels
    edge::data::TernaryXsmm< TL_T_REAL > t_ternary;

    /**
     * Gets the sparse matrix structure of the dense input matrices.
     * 
     * @param i_fIntL local flux matrices.
     * @param i_fIntN neighboring flux matrices.
     * @param i_fIntT transposed flux matrices.
     * @param o_offsets will be set to the offsets (counting non-zero entries) of the sparse matrices.
     * @param o_nonZeros will be set to the raw non-zero entries of the sparse matrices.
     * @param o_mats will be set to the CSC representation of the sparse matrices.
     **/
    static void getCscFlux( TL_T_REAL                const   i_fIntL[TL_N_FAS][TL_N_MDS_EL][TL_N_MDS_FA],
                            TL_T_REAL                const   i_fIntN[TL_N_FMNS][TL_N_MDS_EL][TL_N_MDS_FA],
                            TL_T_REAL                const   i_fIntT[TL_N_FAS][TL_N_MDS_FA][TL_N_MDS_EL],
                            std::vector< size_t >          & o_offsets,
                            std::vector< TL_T_REAL >       & o_nonZeros,
                            std::vector< t_matCsc >        & o_mats ) {
      // determine CSC fill-in strategy
      std::string l_cscFillIn = "none";

      if ( libxsmm_get_target_archid() == LIBXSMM_X86_AVX512_KNM ) {
        l_cscFillIn = "qfma";
      }

      // reset and init output
      o_offsets.resize( 0 );
      o_offsets.push_back( 0 );
      o_nonZeros.resize( 0 );
      o_mats.resize( 0 );

      // local contribution flux matrices
      for( unsigned short l_fl = 0; l_fl < TL_N_FAS; l_fl++ ) {
        t_matCsc l_fluxCsc;

        edge::linalg::Matrix::denseToCsc< TL_T_REAL >( TL_N_MDS_EL,
                                                       TL_N_MDS_FA,
                                                       i_fIntL[l_fl][0],
                                                       l_fluxCsc,
                                                       TOL.BASIS,
                                                       std::numeric_limits< unsigned int >::max(),
                                                       std::numeric_limits< unsigned int >::max(),
                                                       l_cscFillIn );
        o_mats.push_back( l_fluxCsc );

        for( unsigned short l_nz = 0; l_nz < l_fluxCsc.val.size(); l_nz++ ) {
          o_nonZeros.push_back( l_fluxCsc.val[l_nz] );
        }
        o_offsets.push_back( o_offsets.back() + l_fluxCsc.val.size() );
      }

      // neighboring contribution flux matrices
      for( unsigned short l_fn = 0; l_fn < TL_N_FMNS; l_fn++ ) {
        t_matCsc l_fluxCsc;

        edge::linalg::Matrix::denseToCsc< TL_T_REAL >( TL_N_MDS_EL,
                                                       TL_N_MDS_FA,
                                                       i_fIntN[l_fn][0],
                                                       l_fluxCsc,
                                                       TOL.BASIS,
                                                       std::numeric_limits< unsigned int >::max(),
                                                       std::numeric_limits< unsigned int >::max(),
                                                       l_cscFillIn );
        o_mats.push_back( l_fluxCsc );

        for( unsigned short l_nz = 0; l_nz < l_fluxCsc.val.size(); l_nz++ ) {
          o_nonZeros.push_back( l_fluxCsc.val[l_nz] );
        }
        o_offsets.push_back( o_offsets.back() + l_fluxCsc.val.size() );
      }

      // transposed flux matrices
      for( unsigned short l_ft = 0; l_ft < TL_N_FAS; l_ft++ ) {
        t_matCsc l_fluxCsc;

        edge::linalg::Matrix::denseToCsc< TL_T_REAL >( TL_N_MDS_FA,
                                                       TL_N_MDS_EL,
                                                       i_fIntT[l_ft][0],
                                                       l_fluxCsc,
                                                       TOL.BASIS,
                                                       std::numeric_limits< unsigned int >::max(),
                                                       std::numeric_limits< unsigned int >::max(),
                                                       l_cscFillIn );
        o_mats.push_back( l_fluxCsc );

        for( unsigned short l_nz = 0; l_nz < l_fluxCsc.val.size(); l_nz++ ) {
          o_nonZeros.push_back( l_fluxCsc.val[l_nz] );
        }
        o_offsets.push_back( o_offsets.back() + l_fluxCsc.val.size() );
      }
    }

    /**
     * Generates the matrix kernels for the flux matrices and flux solvers.
     *
     * @param i_fIntL local flux matrices.
     * @param i_fIntN neighboring flux matrices.
     * @param i_fIntT transposed flux matrices.
     **/
    void generateKernels( TL_T_REAL const   i_fIntL[TL_N_FAS][TL_N_MDS_EL][TL_N_MDS_FA],
                          TL_T_REAL const   i_fIntN[TL_N_FMNS][TL_N_MDS_EL][TL_N_MDS_FA],
                          TL_T_REAL const   i_fIntT[TL_N_FAS][TL_N_MDS_FA][TL_N_MDS_EL] ) {
      // convert flux matrices to CSC (incl. possible fill-in)
      std::vector< size_t > l_offsets;
      std::vector< TL_T_REAL > l_nonZeros;
      std::vector< t_matCsc > l_fIntCsc;
      getCscFlux( i_fIntL,
                  i_fIntN,
                  i_fIntT,
                  l_offsets,
                  l_nonZeros,
                  l_fIntCsc );

      // local contribution flux matrices
      for( unsigned short l_fl = 0; l_fl < TL_N_FAS; l_fl++ ) {
        m_mm.add(  0,                         // group
                   TL_N_CRS,                  // number of sims
                   false,                     // csc
                  &l_fIntCsc[l_fl].colPtr[0], // column pointer
                  &l_fIntCsc[l_fl].rowIdx[0], // row index
                  &l_fIntCsc[l_fl].val[0],    // values
                   TL_N_QTS_E,                // m
                   TL_N_MDS_FA,               // n
                   TL_N_MDS_EL,               // k
                   TL_N_MDS_EL,               // ldA
                   0,                         // ldB
                   TL_N_MDS_FA,               // ldC
                   TL_T_REAL(1.0),            // alpha
                   TL_T_REAL(0.0),            // beta
                   LIBXSMM_GEMM_PREFETCH_NONE );
      }

      // neighboring contribution flux matrices
      for( unsigned short l_fn = 0; l_fn < TL_N_FMNS; l_fn++ ) {
        unsigned short l_ma = TL_N_FAS + l_fn;

        m_mm.add(  0,                         // group
                   TL_N_CRS,                  // number of sims
                   false,                     // csc
                  &l_fIntCsc[l_ma].colPtr[0], // column pointer
                  &l_fIntCsc[l_ma].rowIdx[0], // row index
                  &l_fIntCsc[l_ma].val[0],    // values
                   TL_N_QTS_E,                // m
                   TL_N_MDS_FA,               // n
                   TL_N_MDS_EL,               // k
                   TL_N_MDS_EL,               // ldA
                   0,                         // ldB
                   TL_N_MDS_FA,               // ldC
                   TL_T_REAL(1.0),            // alpha
                   TL_T_REAL(0.0),            // beta
                   LIBXSMM_GEMM_PREFETCH_NONE );
      }

      // transposed flux matrices for elastic update
      for( unsigned short l_ft = 0; l_ft < TL_N_FAS; l_ft++ ) {
        unsigned short l_ma = TL_N_FAS + TL_N_FMNS + l_ft;

        m_mm.add(  2,                         // group
                   TL_N_CRS,                  // number of sims
                   false,                     // csc
                  &l_fIntCsc[l_ma].colPtr[0], // column pointer
                  &l_fIntCsc[l_ma].rowIdx[0], // row index
                  &l_fIntCsc[l_ma].val[0],    // values
                   TL_N_QTS_E,                // m
                   TL_N_MDS_EL,               // n
                   TL_N_MDS_FA,               // k
                   TL_N_MDS_FA,               // ldA
                   0,                         // ldB
                   TL_N_MDS_EL,               // ldC
                   TL_T_REAL(1.0),            // alpha
                   TL_T_REAL(1.0),            // beta
                   LIBXSMM_GEMM_PREFETCH_NONE );
      }

      // transposed flux matrices for viscoelastic update
      for( unsigned short l_ft = 0; l_ft < TL_N_FAS; l_ft++ ) {
        unsigned short l_ma = TL_N_FAS + TL_N_FMNS + l_ft;

        m_mm.add(  4,                         // group
                   TL_N_CRS,                  // number of sims
                   false,                     // csc
                  &l_fIntCsc[l_ma].colPtr[0], // column pointer
                  &l_fIntCsc[l_ma].rowIdx[0], // row index
                  &l_fIntCsc[l_ma].val[0],    // values
                   TL_N_QTS_M,                // m
                   TL_N_MDS_EL,               // n
                   TL_N_MDS_FA,               // k
                   TL_N_MDS_FA,               // ldA
                   0,                         // ldB
                   TL_N_MDS_EL,               // ldC
                   TL_T_REAL(1.0),            // alpha
                   TL_T_REAL(1.0),            // beta
                   LIBXSMM_GEMM_PREFETCH_NONE );
      }

      // elastic flux solver
      t_matCsr l_fsCsrE;
      FakeMats< TL_N_DIS >::fsCsrE( l_fsCsrE );
      EDGE_CHECK_EQ( l_fsCsrE.val.size(), TL_N_QTS_E*TL_N_QTS_E );

      m_mm.add(  1,                   // group
                 TL_N_CRS,            // number of sims
                 true,                // csr
                 &l_fsCsrE.rowPtr[0], // row pointer
                 &l_fsCsrE.colIdx[0], // column index
                 &l_fsCsrE.val[0],    // values
                 TL_N_QTS_E,          // m
                 TL_N_MDS_FA,         // n
                 TL_N_QTS_E,          // k
                 0,                   // ldA
                 TL_N_MDS_FA,         // ldB
                 TL_N_MDS_FA,         // ldC
                 TL_T_REAL(1.0),      // alpha
                 TL_T_REAL(0.0),      // beta
                 LIBXSMM_GEMM_PREFETCH_BL2_VIA_C );

      // anelastic flux solver
      if( TL_N_RMS > 0 ) {
        t_matCsr l_fsCsrA;
        FakeMats< TL_N_DIS >::fsCsrA( l_fsCsrA );
        EDGE_CHECK_EQ( l_fsCsrA.val.size(), TL_N_QTS_M*TL_N_QTS_E );

        m_mm.add(  3,                   // group
                   TL_N_CRS,            // number of sims
                   true,                // csr
                   &l_fsCsrA.rowPtr[0], // row pointer
                   &l_fsCsrA.colIdx[0], // column index
                   &l_fsCsrA.val[0],    // values
                   TL_N_QTS_M,          // m
                   TL_N_MDS_FA,         // n
                   TL_N_QTS_M,          // k
                   0,                   // ldA
                   TL_N_MDS_FA,         // ldB
                   TL_N_MDS_FA,         // ldC
                   TL_T_REAL(1.0),      // alpha
                   TL_T_REAL(0.0),      // beta
                   LIBXSMM_GEMM_PREFETCH_NONE );

#ifdef ELTWISE_TPP
        // zeroing the anelastic update
        u_unary.add(0, TL_N_CRS, TL_N_MDS_EL * TL_N_QTS_M /* m, n */, LIBXSMM_MELTW_TYPE_UNARY_XOR, LIBXSMM_MELTW_FLAG_UNARY_NONE);
#endif
      }
    }

    /**
     * Stores the flux matrices.
     * 
     * @param i_fIntL local flux matrices.
     * @param i_fIntN neighboring flux matrices.
     * @param i_fIntT transposed flux matrices.
     * @param io_dynMem dynamic memory management, which will be used for the respective allocations.
     * @param o_fIntLN will contain pointers to memory for the local and neighboring flux matrices.
     * @param o_fIntT will contain pointers to memory for the transposed flux matrices.
     **/
    static void storeFluxSparse( TL_T_REAL     const   i_fIntL[TL_N_FAS][TL_N_MDS_EL][TL_N_MDS_FA],
                                 TL_T_REAL     const   i_fIntN[TL_N_FMNS][TL_N_MDS_EL][TL_N_MDS_FA],
                                 TL_T_REAL     const   i_fIntT[TL_N_FAS][TL_N_MDS_FA][TL_N_MDS_EL],
                                 data::Dynamic       & io_dynMem,
                                 TL_T_REAL           * o_fIntLN[TL_N_FAS+TL_N_FMNS],
                                 TL_T_REAL           * o_fIntT[TL_N_FAS] ) {
      // convert flux matrices to CSC (incl. possible fill-in)
      std::vector< size_t > l_offsets;
      std::vector< TL_T_REAL > l_nonZeros;
      std::vector< t_matCsc > l_fIntCsc;
      getCscFlux( i_fIntL,
                  i_fIntN,
                  i_fIntT,
                  l_offsets,
                  l_nonZeros,
                  l_fIntCsc );

      // copy sparse matrices to a permanent data structure
      TL_T_REAL * l_fIntRaw = (TL_T_REAL*) io_dynMem.allocate( l_nonZeros.size() * sizeof(TL_T_REAL),
                                                               4096,
                                                               true );
      for( std::size_t l_en = 0; l_en < l_nonZeros.size(); l_en++ ) {
        l_fIntRaw[l_en] = l_nonZeros[l_en];
      }

      // check that we have offsets for all matrices
      EDGE_CHECK_EQ( l_offsets.size(), TL_N_FAS + TL_N_FMNS + TL_N_FAS + 1 );

      // set the pointers
      for( unsigned short l_ma = 0; l_ma < TL_N_FAS+TL_N_FMNS; l_ma++ ) {
        o_fIntLN[l_ma] = l_fIntRaw + l_offsets[l_ma];
      }
      for( unsigned short l_fa = 0; l_fa < TL_N_FAS; l_fa++ ) {
        o_fIntT[l_fa] = l_fIntRaw + l_offsets[TL_N_FAS+TL_N_FMNS + l_fa];
      }
    }

    inline void execMmKernel( unsigned short         i_group, 
                              unsigned short         i_kernel,
                              TL_T_REAL      const * i_a,
                              TL_T_REAL      const * i_b,
                              TL_T_REAL            * io_c ) {
#ifdef PP_MMKERNEL_PERF
      size_t start = _rdtsc();
#endif
      //m_mm.m_kernels[i_group][i_kernel]( i_a, i_b, io_c );
      m_mm.execute(i_group, i_kernel, i_a, i_b, io_c);
#ifdef PP_MMKERNEL_PERF
      size_t duration = _rdtsc() - start;
      m_mm.m_kernelStats[edge::parallel::g_thread][i_group][i_kernel].invocations++;
      m_mm.m_kernelStats[edge::parallel::g_thread][i_group][i_kernel].cycles += duration; 
#endif
    }

    inline void execMmKernelPf( unsigned short         i_group, 
                                unsigned short         i_kernel,
                                TL_T_REAL      const * i_a,
                                TL_T_REAL      const * i_b,
                                TL_T_REAL            * io_c,
                                TL_T_REAL      const * i_a_pf,
                                TL_T_REAL      const * i_b_pf,
                                TL_T_REAL      const * i_c_pf  ) {
#ifdef PP_MMKERNEL_PERF
      size_t start = _rdtsc();
#endif
      //m_mm.m_kernels[i_group][i_kernel]( i_a, i_b, io_c, i_a_pf, i_b_pf, i_c_pf );
      m_mm.execute(i_group, i_kernel, i_a, i_b, io_c, i_a_pf, i_b_pf, i_c_pf);
#ifdef PP_MMKERNEL_PERF
      size_t duration = _rdtsc() - start;
      m_mm.m_kernelStats[edge::parallel::g_thread][i_group][i_kernel].invocations++;
      m_mm.m_kernelStats[edge::parallel::g_thread][i_group][i_kernel].cycles += duration; 
#endif
    }

  public:
    /**
     * Constructor of the fused surface integration.
     *
     * @param i_rfs relaxation frequencies, use nullptr if TL_N_RMS==0.
     * @param io_dynMem dynamic memory allocations.
     **/
    SurfIntFused( TL_T_REAL     const * i_rfs,
                  data::Dynamic       & io_dynMem ): SurfInt< TL_T_REAL,
                                                              TL_N_RMS,
                                                              TL_T_EL,
                                                              TL_O_SP,
                                                              TL_N_CRS >( i_rfs,
                                                                          io_dynMem ) {
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

      // store flux matrices sparse
      storeFluxSparse( l_fIntL,
                       l_fIntN,
                       l_fIntT,
                       io_dynMem,
                       m_fIntLN,
                       m_fIntT );

      // generate kernels
      generateKernels( l_fIntL,
                       l_fIntN,
                       l_fIntT );
    }

    /**
     * Element local contribution for fused seismic simulations.
     *
     * @param i_fsE elastic flux solvers.
     * @param i_fsA anelastic flux solvers, use nullptr if TL_N_RMS==0.
     * @param i_tDofsE elastic time integerated DG-DOFs.
     * @param io_dofsE will be updated with local elastic contribution of the element to the surface integral.
     * @param io_dofsA will be updated with local anelastic contribution of the element to the surface integral, use nullptr for TL_N_RMS==0.
     * @param o_scratch will be used as scratch space for the computations.
     * @param i_dofsP DOFs for prefetching (not used).
     * @param i_tDofsP time integrated DOFs for prefetching (not used).
     **/
    void local( TL_T_REAL const   i_fsE[TL_N_FAS][TL_N_ENS_FS_E],
                TL_T_REAL const (*i_fsA)[TL_N_ENS_FS_A],
                TL_T_REAL const   i_tDofsE[TL_N_QTS_E][TL_N_MDS_EL][TL_N_CRS],
                TL_T_REAL         io_dofsE[TL_N_QTS_E][TL_N_MDS_EL][TL_N_CRS],
                TL_T_REAL       (*io_dofsA)[TL_N_QTS_M][TL_N_MDS_EL][TL_N_CRS],
                TL_T_REAL         o_scratch[2][TL_N_QTS_E][TL_N_MDS_FA][TL_N_CRS],
                TL_T_REAL const   i_dofsP[TL_N_QTS_E][TL_N_MDS_EL][TL_N_CRS] = nullptr,
                TL_T_REAL const   i_tDofsP[TL_N_QTS_E][TL_N_MDS_EL][TL_N_CRS] = nullptr ) {
      // anelastic update
      TL_T_REAL l_upAn[TL_N_QTS_M][TL_N_MDS_EL][TL_N_CRS];
      if( TL_N_RMS > 0 ) {
#ifdef ELTWISE_TPP
        u_unary.execute (0, 0, (TL_T_REAL*)&l_upAn[0][0][0]);
#else
        for( unsigned short l_qt = 0; l_qt < TL_N_QTS_M; l_qt++ ) {
          for( unsigned short l_md = 0; l_md < TL_N_MDS_EL; l_md++ ) {
#pragma omp simd
            for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {
              l_upAn[l_qt][l_md][l_cr] = 0;
            }
          }
        }
#endif
      }

      // iterate over faces
      for( unsigned short l_fa = 0; l_fa < TL_N_FAS; l_fa++ ) {
        // local flux matrix
        execMmKernel( 0, l_fa, i_tDofsE[0][0],
                               m_fIntLN[l_fa],
                               o_scratch[0][0][0] );

        // flux solver
        execMmKernelPf( 1, 0, i_fsE[l_fa],
                              o_scratch[0][0][0],
                              o_scratch[1][0][0],
                              nullptr,
                              (l_fa < TL_N_FAS_DIV2) ? i_dofsP[0][0] : i_tDofsP[0][0],
                              nullptr );

        // transposed flux matrix
        execMmKernel( 2, l_fa, o_scratch[1][0][0],
                               m_fIntT[l_fa],
                               io_dofsE[0][0] );

        if( TL_N_RMS > 0 ) {
          // anelastic flux solver
          execMmKernel( 3, 0, i_fsA[l_fa],
                              o_scratch[0][0][0],
                              o_scratch[1][0][0] );

          // transposed flux matrix
          execMmKernel( 4, l_fa, o_scratch[1][0][0],
                                 m_fIntT[l_fa],
                                 l_upAn[0][0] );
        }
      }

      // scatter to anelastic DOFs
      if( TL_N_RMS > 0) this->scatterUpdateA( l_upAn, io_dofsA );
    }

    /**
     * Applies the first first face-integration matrix to the elastic DOFs.
     *
     * @param i_fa local face.
     * @param i_vId id of the vertex, matching the element's vertex 0, from the perspective of the adjacent element w.r.t. to the reference element.
     * @param i_fId id of the face from the perspective of the adjacent element w.r.t. to the reference element.
     * @param i_tDofsE elastic time integrated DOFs.
     * @param o_tDofsFiE elastic time integrated DOFs after application of the first face-int matrix.
     * @param i_pre DOFs or tDOFs for prefetching.
     **/
    void neighFluxInt( unsigned short       i_fa,
                       unsigned short       i_vId,
                       unsigned short       i_fId,
                       TL_T_REAL      const i_tDofsE[TL_N_QTS_E][TL_N_MDS_EL][TL_N_CRS],
                       TL_T_REAL            o_tDofsFiE[TL_N_QTS_E][TL_N_MDS_FA][TL_N_CRS],
                       TL_T_REAL      const i_pre[TL_N_QTS_E][TL_N_MDS_EL][TL_N_CRS] = nullptr ) {
      // derive the id of the neighboring flux matrix
      unsigned short l_fMatId = std::numeric_limits< unsigned short >::max();
      if( i_vId != std::numeric_limits< unsigned short >::max() ) {
        l_fMatId = TL_N_FAS + this->fMatId( i_vId, i_fId );
      }
      else {
        l_fMatId = i_fa;
      }

      // local or neighboring flux matrix
      execMmKernel( 0, l_fMatId, i_tDofsE[0][0],
                                 m_fIntLN[l_fMatId],
                                 o_tDofsFiE[0][0] );
    }

    /**
     * Neighboring contribution of a single adjacent element for fused simulations.
     *
     * @param i_fa local face.
     * @param i_vId id of the vertex, matching the element's vertex 0, from the perspective of the adjacent element w.r.t. to the reference element.
     * @param i_fId id of the face from the perspective of the adjacent element w.r.t. to the reference element.
     * @param i_fsE elastic flux solver.
     * @param i_fsA anelastic flux solver
     * @param i_tDofsE elastic time integrated DG-DOFs.
     * @param io_dofsE will be updated with the elastic contribution of the adjacent element to the surface integral.
     * @param io_dofsA will be updated with the unscaled (w.r.t. frequencies) anelastic contribution of the adjacent element tot the surface integral, use nullptr for TL_N_RMS==0.
     * @param o_scratch will be used as scratch space for the computations.
     * @param i_pre DOFs or tDOFs for prefetching.
     **/
    void neigh( unsigned short       i_fa,
                unsigned short       i_vId,
                unsigned short       i_fId,
                TL_T_REAL      const i_fsE[TL_N_ENS_FS_E],
                TL_T_REAL      const i_fsA[TL_N_ENS_FS_A],
                TL_T_REAL      const i_tDofsE[TL_N_QTS_E][TL_N_MDS_EL][TL_N_CRS],
                TL_T_REAL      const i_tDofsFiE[TL_N_QTS_E][TL_N_MDS_FA][TL_N_CRS],
                TL_T_REAL            io_dofsE[TL_N_QTS_E][TL_N_MDS_EL][TL_N_CRS],
                TL_T_REAL            io_dofsA[TL_N_QTS_M][TL_N_MDS_EL][TL_N_CRS],
                TL_T_REAL            o_scratch[2][TL_N_QTS_E][TL_N_MDS_FA][TL_N_CRS],
                TL_T_REAL      const i_pre[TL_N_QTS_E][TL_N_MDS_EL][TL_N_CRS] = nullptr ) {
      // apply first face integration matrix or set pointer to pre-computed data
      TL_T_REAL const * l_tDofsFiE = o_scratch[0][0][0];

      if( i_tDofsFiE == nullptr ) {
        neighFluxInt( i_fa,
                      i_vId,
                      i_fId,
                      i_tDofsE,
                      o_scratch[0],
                      i_pre );
      }
      else {
        l_tDofsFiE = i_tDofsFiE[0][0];
      }

      // flux solver
      execMmKernelPf( 1, 0,                     i_fsE,
                                                l_tDofsFiE,
                                                o_scratch[1][0][0],
                                                nullptr,
                            (TL_T_REAL const *) i_pre,
                                                nullptr );

      // transposed flux matrix
      execMmKernel( 2, i_fa, o_scratch[1][0][0],
                             m_fIntT[i_fa],
                             io_dofsE[0][0] );

      if( TL_N_RMS > 0 ) {
        // anelastic flux solver
        execMmKernel( 3, 0, i_fsA,
                            l_tDofsFiE,
                            o_scratch[1][0][0] );

        // transposed flux matrix
        execMmKernel( 4, i_fa, o_scratch[1][0][0],
                               m_fIntT[i_fa],
                               io_dofsA[0][0] );
      }
    }
};

#endif
