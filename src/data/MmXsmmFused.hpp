/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (alex.breuer AT uni-jena.de)
 *         Alexander Heinecke (alexander.heinecke AT intel.com)
 *
 * @section LICENSE
 * Copyright (c) 2020, Friedrich Schiller University Jena
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
 * Data structures of the fused LIBXSMM, matrix-matrix multiplication kernels.
 **/

#ifndef EDGE_DATA_MM_XSMM_FUSED_HPP
#define EDGE_DATA_MM_XSMM_FUSED_HPP

#include <vector>
#include "constants.hpp"
#include "io/logging.h"
#include "linalg/Matrix.h"
#include "parallel/global.h"

#include "XsmmUtils.hpp"

#include <libxsmm.h>

namespace edge {
  namespace data {
    template< typename TL_T_REAL >
    class MmXsmmFused;

    typedef struct MmXsmmStats {
      size_t invocations;
      size_t cycles;
    } MmXsmmStats;
  }
}

/**
 * Holds LIBXSMM kernels for fused simulations.
 **/
template<typename TL_T_REAL>
class edge::data::MmXsmmFused {
  public:
    //! generated kernels of libxsmm
    std::vector< std::vector< libxsmm_gemmfunction > > m_kernels;

    //! number of flops performed by each libxsmm kernel
    std::vector< std::vector< size_t > > m_kernelFlops;

    //! stats for kernels
    std::vector< std::vector< std::vector< MmXsmmStats > > > m_kernelStats;

    /**
     * @brief Constructor, which limits the LIBXSMM target architecture, if required.
     */
    MmXsmmFused() {
      m_kernelStats.resize(edge::parallel::g_nThreads);
    }

    /**
     * Adds a sparse libxsmm-kernel for the given matrix in CSR- or CSC-format.
     *
     * @param i_group id of the kernel group.
     * @param i_nCrs number of used simulations.
     * @param i_csr true for CSR, false for CSC.
     * @param i_ptr row-pointer of CSR or column pointer of CSC; last element holds the number of non-zero entries.
     * @param i_idx column index of CSR or row index of CSC.
     * @param i_val non-zero values.
     * @param i_m number of rows in column-major A and C.
     * @param i_n number of columns in column-major B and C.
     * @param i_k number of columns/rows in column-major A/B.
     * @param i_ldA leading dimension of column-major A.
     * @param i_ldB leading dimension of column-major B.
     * @param i_ldC leading dimension of column-major C.
     * @param i_alpha alpha parameter (needs to be 1.0 for now).
     * @param i_beta beta parameter (need to be 0.0/1.0 for now).
     * @param i_prefetch prefetch strategy.
     **/
    void add( unsigned short                    i_group,
              unsigned short                    i_nCrs,
              bool                              i_csr,
              unsigned int               const *i_ptr,
              unsigned int               const *i_idx,
              TL_T_REAL                  const *i_val,
              unsigned int                      i_m,
              unsigned int                      i_n,
              unsigned int                      i_k,
              unsigned int                      i_ldA,
              unsigned int                      i_ldB,
              unsigned int                      i_ldC,
              double                            i_alpha,
              double                            i_beta,
              libxsmm_gemm_prefetch_type        i_prefetch ) {
      EDGE_VLOG(1) << "  adding X precision XSMM-kernel #" << m_kernels.size() << " (sparse)"
                   << " M=" << i_m << " N=" << i_n << " K=" << i_k
                   << " ldA=" << i_ldA << " ldB=" << i_ldB << " ldC=" << i_ldC
                   << " alpha=" << i_alpha << " beta=" << i_beta;

      // add kernel groups, if required
      if( i_group >= m_kernels.size() ) {
        m_kernels.resize( i_group+1 );
        m_kernelFlops.resize( i_group+1 );
        for ( int i = 0; i < edge::parallel::g_nThreads; ++i ) {
          m_kernelStats[i].resize( i_group+1 );
        }
      }

      // add description
      /*
      libxsmm_descriptor_blob l_xgemmBlob;
      const libxsmm_gemm_descriptor* l_desc = 0;
      const int l_flags = LIBXSMM_GEMM_FLAGS('N', 'N');
      l_desc = libxsmm_gemm_descriptor_dinit(&l_xgemmBlob, LIBXSMM_GEMM_PRECISION_F64,
        i_m, i_n, i_k, i_ldA, i_ldB, i_ldC, i_alpha, i_beta, l_flags, i_prefetch);

      m_descs[i_group].push_back( l_desc );
      */

      // @FIXME: Do we need LIBXSMM_GEMM_FLAG_USE_XGEMM_ABI?
      int l_flags = LIBXSMM_GEMM_FLAGS('N', 'N');
      if (i_beta == 0.0)
          l_flags |= LIBXSMM_GEMM_FLAG_BETA_0;

      libxsmm_datatype l_dtype_A    = XsmmDtype<TL_T_REAL>();
      libxsmm_datatype l_dtype_B    = XsmmDtype<TL_T_REAL>();
      libxsmm_datatype l_dtype_C    = XsmmDtype<TL_T_REAL>();
      libxsmm_datatype l_dtype_comp = XsmmDtype<TL_T_REAL>();

      libxsmm_gemm_shape l_gemm_shape = libxsmm_create_gemm_shape(i_m, i_n, i_k, i_ldA, i_ldB, i_ldC,
                                                        l_dtype_A, l_dtype_B, l_dtype_C, l_dtype_comp);

      // generate and store function for this kernels
      if( i_csr )
        //m_kernels[i_group].push_back( libxsmm_create_packed_spxgemm_csr( m_descs[i_group].back(), i_nCrs, i_ptr, i_idx, i_val ).dmm );
        m_kernels[i_group].push_back( libxsmm_create_packed_spgemm_csr_v2(l_gemm_shape, l_flags, i_prefetch, i_nCrs,
                                      i_ptr, i_idx, i_val) );
      else
        //m_kernels[i_group].push_back( libxsmm_create_packed_spxgemm_csc( m_descs[i_group].back(), i_nCrs, i_ptr, i_idx, i_val ).dmm );
        m_kernels[i_group].push_back( libxsmm_create_packed_spgemm_csc_v2(l_gemm_shape, l_flags, i_prefetch, i_nCrs,
                                      i_ptr, i_idx, i_val) );

      // check that we generated a kernel
      EDGE_CHECK( m_kernels[i_group].back() != 0 );

      // read flops and store them
      libxsmm_kernel_info l_kinfo;
      libxsmm_get_kernel_info( (const void*)m_kernels[i_group].back(), &l_kinfo );
      m_kernelFlops[i_group].push_back( l_kinfo.nflops );

      // Initalize stats telemetry
      MmXsmmStats l_mystats = { 0 , 0 };
      for ( int i = 0; i < edge::parallel::g_nThreads; ++i ) {
        m_kernelStats[i][i_group].push_back( l_mystats );
      }
    }

    /**
     * Adds a dense libxsmm-kernel.
     *
     * @param i_group id of the kernel group.
     * @param i_nCrs number of used simulations.
     * @param i_m number of rows in column-major A and C.
     * @param i_n number of columns in column-major B and C.
     * @param i_k number of columns/rows in column-major A/B.
     * @param i_ldA leading dimension of column-major A.
     * @param i_ldB leading dimension of column-major B.
     * @param i_ldC leading dimension of column-major C.
     * @param i_alpha alpha parameter (needs to be 1.0 for now).
     * @param i_beta beta parameter (need to be 0.0/1.0 for now).
     * @param i_fusedAC true if matrices A and C are fused.
     * @param i_fusedBC true if matrices B and C are fused.
     * @param i_prefetch prefetch strategy.
     **/
    void add( unsigned short             i_group,
              unsigned short             i_nCrs,
              unsigned int               i_m,
              unsigned int               i_n,
              unsigned int               i_k,
              unsigned int               i_ldA,
              unsigned int               i_ldB,
              unsigned int               i_ldC,
              TL_T_REAL                  i_alpha,
              TL_T_REAL                  i_beta,
              bool                       i_fusedAC,
              bool                       i_fusedBC,
              libxsmm_gemm_prefetch_type i_prefetch ) {
      // check that a fused kernel is requested
      EDGE_CHECK( i_fusedBC || i_fusedAC );

      // add kernel groups, if required
      if( i_group >= m_kernels.size() ) {
        //m_descs.resize( i_group+1 );
        m_kernels.resize( i_group+1 );
        m_kernelFlops.resize( i_group+1 );
        for ( int i = 0; i < edge::parallel::g_nThreads; ++i ) {
          m_kernelStats[i].resize( i_group+1 );
        }
      }

      // add description
      /*
      libxsmm_descriptor_blob l_xgemmBlob;
      const libxsmm_gemm_descriptor* l_desc = 0;
      const int l_flags = LIBXSMM_GEMM_FLAGS('N', 'N') | LIBXSMM_GEMM_FLAG_USE_XGEMM_ABI;
      l_desc = libxsmm_gemm_descriptor_dinit2( &l_xgemmBlob,
                                               LIBXSMM_DATATYPE_F64, LIBXSMM_DATATYPE_F64,
                                               i_m, i_n, i_k,
                                               (i_fusedBC ? 0 : i_ldA), i_ldB, i_ldC,
                                               i_alpha, i_beta,
                                               l_flags,
                                               i_prefetch);

      m_descs[i_group].push_back( l_desc );
      */

      int l_flags = LIBXSMM_GEMM_FLAGS('N', 'N') | LIBXSMM_GEMM_FLAG_USE_XGEMM_ABI;
      if (i_beta == 0.0)
        l_flags |= LIBXSMM_GEMM_FLAG_BETA_0;

      libxsmm_datatype l_dtype_A    = XsmmDtype<TL_T_REAL>();
      libxsmm_datatype l_dtype_B    = XsmmDtype<TL_T_REAL>();
      libxsmm_datatype l_dtype_C    = XsmmDtype<TL_T_REAL>();
      libxsmm_datatype l_dtype_comp = XsmmDtype<TL_T_REAL>();

      libxsmm_gemm_shape l_gemm_shape = libxsmm_create_gemm_shape(i_m, i_n, i_k, (i_fusedBC ? 0 : i_ldA),
                                                                  i_ldB, i_ldC,  l_dtype_A, l_dtype_B, l_dtype_C, l_dtype_comp);

      if( i_fusedBC ) {
        // generate fake CSR-structure
        unsigned int *l_rows = nullptr;
        unsigned int *l_cols = nullptr;
        double       *l_vals = nullptr;
        linalg::Matrix::fakeCsr( i_m, i_n, i_k,
                                 l_rows, l_cols, l_vals );


        // generate and store function for this kernels
        //m_kernels[i_group].push_back( libxsmm_create_packed_spxgemm_csr( m_descs[i_group].back(), i_nCrs, l_rows, l_cols, l_vals ).dmm );
        m_kernels[i_group].push_back( libxsmm_create_packed_spgemm_csr_v2(l_gemm_shape, l_flags, i_prefetch, i_nCrs,
                                      l_rows, l_cols, l_vals) );

        // free memory of fake CSR-structure
        delete[] l_rows; delete[] l_cols; delete[] l_vals;
      }
      else {
        //m_kernels[i_group].push_back( libxsmm_create_packed_xgemm_ac_rm( m_descs[i_group].back(), i_nCrs ).dmm );
        m_kernels[i_group].push_back( libxsmm_create_packed_gemm_ac_rm_v2( l_gemm_shape, l_flags, i_prefetch, i_nCrs) );
      }

      // check that we generated a kernel
      EDGE_CHECK( m_kernels[i_group].back() != 0 );

      // read flops and store them
      libxsmm_kernel_info l_kinfo;
      libxsmm_get_kernel_info( (const void*)m_kernels[i_group].back(), &l_kinfo );
      m_kernelFlops[i_group].push_back( l_kinfo.nflops );

      // Initalize stats telemetry
      MmXsmmStats l_mystats = { 0 , 0 };
      for ( int i = 0; i < edge::parallel::g_nThreads; ++i ) {
        m_kernelStats[i][i_group].push_back( l_mystats );
      }
    }

    void execute( const unsigned short i_group,
                  const unsigned short i_entry,
                  TL_T_REAL      const * i_a,
                  TL_T_REAL      const * i_b,
                  TL_T_REAL            * io_c) const {

      libxsmm_gemm_param gemm_param;
      memset( &gemm_param, 0, sizeof(libxsmm_gemm_param) );

      gemm_param.a.primary    = (void*)i_a;
      gemm_param.b.primary    = (void*)i_b;
      gemm_param.c.primary    = (void*)io_c;

      m_kernels[i_group][i_entry]( &gemm_param );
    }

    void execute( const unsigned short i_group,
                  const unsigned short i_entry,
                  TL_T_REAL      const * i_a,
                  TL_T_REAL      const * i_b,
                  TL_T_REAL            * io_c,
                  TL_T_REAL      const * i_a_pf,
                  TL_T_REAL      const * i_b_pf,
                  TL_T_REAL      const * i_c_pf) const {

      libxsmm_gemm_param gemm_param;
      memset( &gemm_param, 0, sizeof(libxsmm_gemm_param) );

      gemm_param.a.primary    = (void*)i_a;
      gemm_param.b.primary    = (void*)i_b;
      gemm_param.c.primary    = (void*)io_c;

      gemm_param.a.quaternary = (void*)i_a_pf;
      gemm_param.b.quaternary = (void*)i_b_pf;
      gemm_param.c.quaternary = (void*)i_c_pf;

      m_kernels[i_group][i_entry]( &gemm_param );
    }

};
#endif

