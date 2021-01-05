/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *         Alexander Heinecke (alexander.heinecke AT intel.com)
 *
 * @section LICENSE
 * Copyright (c) 2020, Friedrich Schiller University Jena
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
 * Optimized quadrature-free ADER-DG volume integration for fused seismic simulations.
 **/
#ifndef EDGE_SEISMIC_KERNELS_VOL_INT_FUSED_HPP
#define EDGE_SEISMIC_KERNELS_VOL_INT_FUSED_HPP

#include "VolInt.hpp"
#include "dg/Basis.h"
#include "data/MmXsmmFused.hpp"
#include "FakeMats.hpp"

namespace edge {
  namespace seismic {
    namespace kernels { 
      template< typename       TL_T_REAL,
                unsigned short TL_N_RMS,
                t_entityType   TL_T_EL,
                unsigned short TL_O_SP,
                unsigned short TL_N_CRS >
      class VolIntFused;
    }
  }
}

/**
 * Optimized quadrature-free ADER-DG volume integration for fused seismic simulations.
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
class edge::seismic::kernels::VolIntFused: edge::seismic::kernels::VolInt< TL_T_REAL,
                                                                           TL_N_RMS,
                                                                           TL_T_EL,
                                                                           TL_O_SP,
                                                                           TL_N_CRS > {
  private:
    //! dimension of the element
    static unsigned short const TL_N_DIS = C_ENT[TL_T_EL].N_DIM;

    //! number of element modes
    static unsigned short const TL_N_MDS = CE_N_ELEMENT_MODES( TL_T_EL, TL_O_SP );

    //! number of elastic quantities
    static unsigned short const TL_N_QTS_E = CE_N_QTS_E( TL_N_DIS );

    //! number of quantities per relaxation mechanism
    static unsigned short const TL_N_QTS_M = CE_N_QTS_M( TL_N_DIS );

    //! number of non-zeros in the star matrices
    static unsigned short const TL_N_NZS_STAR_E = CE_N_ENS_STAR_E_SP( TL_N_DIS );

    //! number of non-zeros in the anelastic star matrices
    static unsigned short const TL_N_NZS_STAR_A = CE_N_ENS_STAR_A_SP( TL_N_DIS );

    //! number of non-zeros in the anelastic soruce matrices
    static unsigned short const TL_N_NZS_SRC_A = CE_N_ENS_SRC_A_SP( TL_N_DIS );

    //! matrix kernels
    edge::data::MmXsmmFused< TL_T_REAL > m_mm;

    //! pointers to the stiffness matrices
    TL_T_REAL *m_stiff[TL_N_DIS] = {};

    /**
     * Gets the sparse matrix structure of the dense input matrices.
     * 
     * @param i_stiff stiffness matrices.
     * @param o_maxNzRow will be set to maximum number non-zero rows of the three stiffness matrices.
     * @param o_offsets will be set to the offsets (counting non-zero entries) of the sparse matrices.
     * @param o_nonZeros will be set to the raw non-zero entries of the sparse matrices.
     * @param o_mats will be set to the CSC representation of the sparse matrices.
     **/
    static void getCscStiff( TL_T_REAL               const   i_stiff[TL_N_DIS][TL_N_MDS][TL_N_MDS],
                             unsigned int                  & o_maxNzRow,
                             std::vector< size_t >         & o_offsets,
                             std::vector< TL_T_REAL >      & o_nonZeros,
                             std::vector< t_matCsc >       & o_mats ) {
      // determine CSC fill-in strategy
      std::string l_cscFillIn = "none";

      if ( libxsmm_get_target_archid() == LIBXSMM_X86_AVX512_KNM ) {
        l_cscFillIn = "qfma";
      }

      // reset and init output
      o_maxNzRow = 0;
      o_offsets.resize( 0 );
      o_offsets.push_back( 0 );
      o_nonZeros.resize( 0 );
      o_mats.resize( 0 );

      // derive the non-zeros of the volume integration
      t_matCrd l_stiffCrd[TL_N_DIS];
      for( unsigned short l_di = 0; l_di < TL_N_DIS; l_di++ ) {
        edge::linalg::Matrix::denseToCrd< TL_T_REAL >( TL_N_MDS,
                                                       TL_N_MDS,
                                                       i_stiff[l_di][0],
                                                       l_stiffCrd[l_di],
                                                       TOL.BASIS );
      }

      // nz-blocks
      unsigned int l_nzBl[2][2][2];
      // init with matrix dim
      l_nzBl[0][0][0] = l_nzBl[0][1][0] = 0;
      l_nzBl[0][0][1] = l_nzBl[0][1][1] = TL_N_MDS-1;

      // get max #nz-rows
      for( unsigned short l_di = 0; l_di < N_DIM; l_di++ ) {
        edge::linalg::Matrix::getBlockNz( l_stiffCrd[l_di], l_nzBl[0], l_nzBl[1] );
        o_maxNzRow = std::max( o_maxNzRow, l_nzBl[1][0][1] );
      }

#ifdef PP_T_BASIS_HIERARCHICAL
      // check that size is one "order" less
      EDGE_CHECK_EQ( o_maxNzRow+1, CE_N_ELEMENT_MODES( TL_T_EL, TL_O_SP-1 ) );
#endif

      // convert to CSC
      t_matCsc l_stiffCsc[TL_N_DIS];
      for( unsigned short l_di = 0; l_di < TL_N_DIS; l_di++ ) {
        edge::linalg::Matrix::denseToCsc< TL_T_REAL >( TL_N_MDS,
                                                       TL_N_MDS,
                                                       i_stiff[l_di][0],
                                                       l_stiffCsc[l_di],
                                                       TOL.BASIS,
                                                       std::numeric_limits< unsigned int >::max(),
                                                       std::numeric_limits< unsigned int >::max(),
                                                       l_cscFillIn );
        o_mats.push_back( l_stiffCsc[l_di] );
      }

      // add data for the non-zero entries and sparse offsets
      for( unsigned short l_di = 0; l_di < TL_N_DIS; l_di++ ) {
        for( unsigned int l_nz = 0; l_nz < l_stiffCsc[l_di].val.size(); l_nz++ ) {
          o_nonZeros.push_back( l_stiffCsc[l_di].val[l_nz] );
        }
        o_offsets.push_back( o_offsets.back() + l_stiffCsc[l_di].val.size() );
      }
    }

    /**
     * Generates the matrix kernels for the stiffness matrices and star matrices.
     *
     * @param i_stiff dense stiffness matrices.
     **/
    void generateKernels( TL_T_REAL const i_stiff[TL_N_DIS][TL_N_MDS][TL_N_MDS] ) {
      // convert stiffness matrices to CSC (incl. possible fill-ins)
      unsigned int l_maxNzRow;
      std::vector< size_t > l_offsets;
      std::vector< TL_T_REAL > l_nonZeros;
      std::vector< t_matCsc > l_stiffCsc;
      getCscStiff( i_stiff,
                   l_maxNzRow,
                   l_offsets,
                   l_nonZeros,
                   l_stiffCsc );

      // get csr elastic star matrix
      t_matCsr l_starCsrE;
      FakeMats< TL_N_DIS >::starCsrE( l_starCsrE );
      EDGE_CHECK_EQ( l_starCsrE.val.size(), TL_N_NZS_STAR_E );

      // stiffness matrices
      for( unsigned short l_di = 0; l_di < N_DIM; l_di++ ) {
        m_mm.add(  0,                                                 // group
                   TL_N_CRS,                                          // number of sims
                   false,                                             // csc
                  &l_stiffCsc[l_di].colPtr[0],                        // column pointer
                  &l_stiffCsc[l_di].rowIdx[0],                        // row index
                  &l_stiffCsc[l_di].val[0],                           // values
                   TL_N_QTS_E,                                        // m
                   TL_N_MDS,                                          // n
                   l_maxNzRow+1,                                      // k
                   (TL_N_RMS == 0) ? l_maxNzRow+1 : TL_N_MDS,         // ldA
                   0,                                                 // ldB
                   TL_N_MDS,                                          // ldC
                   TL_T_REAL(1.0),                                    // alpha
                   (TL_N_RMS == 0) ? TL_T_REAL(1.0) : TL_T_REAL(0.0), // beta
                   LIBXSMM_GEMM_PREFETCH_NONE ); // Remark: Star matrix is multiplied first in elastic solver
      }

      // elastic star matrix
      m_mm.add(  1,                                                 // group
                 TL_N_CRS,                                          // number of sims
                 true,                                              // csr
                &l_starCsrE.rowPtr[0],                              // row pointer
                &l_starCsrE.colIdx[0],                              // column index
                &l_starCsrE.val[0],                                 // values
                 TL_N_QTS_E,                                        // m
                 (TL_N_RMS == 0) ? l_maxNzRow+1 : TL_N_MDS,         // n
                 TL_N_QTS_E,                                        // k
                 0,                                                 // ldA
                 TL_N_MDS,                                          // ldB
                 (TL_N_RMS == 0) ? l_maxNzRow+1 : TL_N_MDS,         // ldC
                 TL_T_REAL(1.0),                                    // alpha
                 (TL_N_RMS == 0) ? TL_T_REAL(0.0) : TL_T_REAL(1.0), // beta
                 LIBXSMM_GEMM_PREFETCH_NONE );

      if( TL_N_RMS > 0 ) {
        // get csr anelastic source matrix
        t_matCsr l_starCsrA;
        FakeMats< TL_N_DIS >::starCsrA( l_starCsrA );
        EDGE_CHECK_EQ( l_starCsrA.val.size(), TL_N_NZS_STAR_A );

        // get csr source matrix
        t_matCsr l_srcCsrA;
        FakeMats< TL_N_DIS >::srcCsrA( l_srcCsrA );
        EDGE_CHECK_EQ( l_srcCsrA.val.size(), TL_N_NZS_SRC_A );

        // multiplication with anelastic star matrix
        m_mm.add( 2,                     // group
                  TL_N_CRS,              // number of sims
                  true,                  // csr
                  &l_starCsrA.rowPtr[0], // ptr
                  &l_starCsrA.colIdx[0], // ids
                  &l_starCsrA.val[0],    // val
                  TL_N_QTS_M,            // m
                  TL_N_MDS,              // n
                  TL_N_QTS_E,            // k
                  0,                     // ldA
                  TL_N_MDS,              // ldB
                  TL_N_MDS,              // ldC
                  TL_T_REAL(1.0),        // alpha
                  TL_T_REAL(1.0),        // beta
                  LIBXSMM_GEMM_PREFETCH_NONE );

        // multiplication with anelastic source matrices
        m_mm.add( 2,                    // group
                  TL_N_CRS,             // number of sims
                  true,                 // csr
                  &l_srcCsrA.rowPtr[0], // ptr
                  &l_srcCsrA.colIdx[0], // ids
                  &l_srcCsrA.val[0],    // val
                  TL_N_QTS_M,           // m
                  TL_N_MDS,             // n
                  TL_N_QTS_M,           // k
                  0,                    // ldA
                  TL_N_MDS,             // ldB
                  TL_N_MDS,             // ldC
                  TL_T_REAL(1.0),       // alpha
                  TL_T_REAL(1.0),       // beta
                  LIBXSMM_GEMM_PREFETCH_NONE );
        }
    }

    /**
     * Stores the stiffness matrices.
     * 
     * @param i_stiff dense stiffness matrices.
     * @param io_dynMem dynamic memory management, which will be used for the respective allocations.
     * @param o_stiff will contain pointers to memory for the individual matrices.
     **/
    static void storeStiffSparse( TL_T_REAL     const     i_stiff[TL_N_DIS][TL_N_MDS][TL_N_MDS],
                                  data::Dynamic       &   io_dynMem,
                                  TL_T_REAL           *   o_stiff[TL_N_DIS]  ) {
      // convert stiffness matrices to CSC (incl. possible fill-in)
      unsigned int l_maxNzRow;
      std::vector< size_t > l_offsets;
      std::vector< TL_T_REAL > l_nonZeros;
      std::vector< t_matCsc > l_cscStiff;
      getCscStiff( i_stiff,
                   l_maxNzRow,
                   l_offsets,
                   l_nonZeros,
                   l_cscStiff );

      // copy sparse matrices to a permanent data structure
      TL_T_REAL * l_stiffRaw = (TL_T_REAL*) io_dynMem.allocate( l_nonZeros.size() * sizeof(TL_T_REAL),
                                                                4096,
                                                                true );
      for( std::size_t l_en = 0; l_en < l_nonZeros.size(); l_en++ ) {
        l_stiffRaw[l_en] = l_nonZeros[l_en];
      }

      // check that we have offsets for all matrices
      EDGE_CHECK_EQ( l_offsets.size(), TL_N_DIS + 1 );

      // set the pointers
      for( unsigned int l_di = 0; l_di < N_DIM; l_di++ ) {
        o_stiff[l_di] = l_stiffRaw + l_offsets[l_di];
      }
    }


  public:
    /**
     * Constructor of the optimized volume integration for fused forward simulations.
     *
     * @param i_rfs relaxation frequencies, use nullptr if TL_N_RMS==0.
     * @param io_dynMem dynamic memory allocations.
     **/
    VolIntFused( TL_T_REAL     const * i_rfs,
                 data::Dynamic       & io_dynMem ): VolInt< TL_T_REAL,
                                                            TL_N_RMS,
                                                            TL_T_EL,
                                                            TL_O_SP,
                                                            TL_N_CRS >( i_rfs,
                                                                        io_dynMem ) {
      // formulation of the basis in terms of the reference element
      dg::Basis l_basis( TL_T_EL,
                         TL_O_SP );

      // get stiffness matrices
      TL_T_REAL l_stiff[TL_N_DIS][TL_N_MDS][TL_N_MDS];
      l_basis.getStiffMm1Dense( TL_N_MDS,
                                l_stiff[0][0],
                                false );

      // store stiffness matrices dense
      this->storeStiffSparse( l_stiff,
                              io_dynMem,
                              m_stiff );

      generateKernels( l_stiff );
    }

    /**
     * Optimized volume contribution for fused seismic forward simulations.
     *
     * @param i_starE elastic star matrices.
     * @param i_starA anelastic star matrices, use nullptr if TL_N_RMS==0.
     * @param i_tDofsE time integrated elastic DOFs.
     * @param i_tDofsA time integrated anselastic DOFs.
     * @param io_dofsE will be updated with local elastic contribution of the element to the volume integral.
     * @param io_dofsA will be updated with local anelastic contribution of the element to the volume integral, use nullptr if TL_N_RMS==0.
     * @param o_scratch will be used as scratch space for the computations.
     **/
    void apply( TL_T_REAL const   i_starE[TL_N_DIS][TL_N_NZS_STAR_E],
                TL_T_REAL const (*i_starA)[TL_N_NZS_STAR_A],
                TL_T_REAL const (*i_srcA)[TL_N_NZS_SRC_A],
                TL_T_REAL const   i_tDofsE[TL_N_QTS_E][TL_N_MDS][TL_N_CRS],
                TL_T_REAL const (*i_tDofsA)[TL_N_QTS_M][TL_N_MDS][TL_N_CRS],
                TL_T_REAL         io_dofsE[TL_N_QTS_E][TL_N_MDS][TL_N_CRS],
                TL_T_REAL       (*io_dofsA)[TL_N_QTS_M][TL_N_MDS][TL_N_CRS],
                TL_T_REAL         o_scratch[TL_N_QTS_E][TL_N_MDS][TL_N_CRS] ) const {
      // relaxation frequencies
      TL_T_REAL const *l_rfs = this->m_rfs;

      // buffer for the anelastic part
      TL_T_REAL l_scratch[TL_N_QTS_M][TL_N_MDS][TL_N_CRS];
      if( TL_N_RMS > 0 ) {
        for( unsigned short l_qt = 0; l_qt < TL_N_QTS_M; l_qt++ )
          for( unsigned short l_md = 0; l_md < TL_N_MDS; l_md++ )
#pragma omp simd
            for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) l_scratch[l_qt][l_md][l_cr] = 0;
      }

      // set scratch memory to zero for attenuation, which multiplies by stiffness matrices first
      if( TL_N_RMS > 0 ) {
        for( unsigned short l_qt = 0; l_qt < TL_N_QTS_E; l_qt++ )
          for( unsigned short l_md = 0; l_md < TL_N_MDS; l_md++ )
#pragma omp simd
            for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) o_scratch[l_qt][l_md][l_cr] = 0;
      }

      // iterate over dimensions
      for( unsigned short l_di = 0; l_di < TL_N_DIS; l_di++ ) {
        // elastic: utilize zero-block in star matrix multiplication
        if( TL_N_RMS == 0 ) {
          // multiply with star matrix
          m_mm.m_kernels[1][0]( i_starE[l_di],
                                i_tDofsE[0][0],
                                o_scratch[0][0] );

          // multiply with stiffness and inverse mass matrix
          m_mm.m_kernels[0][l_di]( o_scratch[0][0],
                                   m_stiff[l_di],
                                   io_dofsE[0][0] );
        }
        // viscoelastic: re-use stiffness matrix multiplication
        else {
          // multiply with stiffness and inverse mass matrix
          m_mm.m_kernels[0][l_di]( i_tDofsE[0][0],
                                   m_stiff[l_di],
                                   o_scratch[0][0] );

          // multiply with elastic star matrix
          m_mm.m_kernels[1][0]( i_starE[l_di],
                                o_scratch[0][0],
                                io_dofsE[0][0] );

          // multiply with anelastic star matrices
          m_mm.m_kernels[2][0]( i_starA[l_di],
                                o_scratch[0][0],
                                l_scratch[0][0] );
        }
      }

      for( unsigned short l_rm = 0; l_rm < TL_N_RMS; l_rm++ ) {
        // add contribution of source matrix
        m_mm.m_kernels[2][1]( i_srcA[l_rm],
                              i_tDofsA[l_rm][0][0],
                              io_dofsE[0][0] );

        // multiply with relaxation frequency and add
        for( unsigned short l_qt = 0; l_qt < TL_N_QTS_M; l_qt++ ) {
          for( unsigned short l_md = 0; l_md < TL_N_MDS; l_md++ ) {
            for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {
              io_dofsA[l_rm][l_qt][l_md][l_cr] += l_rfs[l_rm] * ( l_scratch[l_qt][l_md][l_cr] - i_tDofsA[l_rm][l_qt][l_md][l_cr] );
            }
          }
        }
      }
    }
};

#endif