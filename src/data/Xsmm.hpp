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
 * Data structures of the XSMM-library.
 **/

#if !defined (PP_T_KERNELS_XSMM) && !defined (PP_T_KERNELS_XSMM_DENSE_SINGLE)
#error error: libxsmm disabled, but compiling data-file.
#endif

#ifndef XSMM_HPP
#define XSMM_HPP

#include "constants.hpp"
#include <libxsmm.h>
#include <io/logging.h>
#include <vector>

namespace edge {
  namespace data {
    class Xsmm;
  }
}

class edge::data::Xsmm {
  //private:
    //! gemm descriptors of libxsmm
    std::vector< libxsmm_gemm_descriptor > m_descs;

  public:
    //! generated kernels of libxsmm
#if PP_PRECISION == 64
    std::vector< libxsmm_dmmfunction > m_kernels;
#else
    std::vector< libxsmm_smmfunction > m_kernels;
#endif

    /**
     * Adds a sparse libxsmm-kernel for the given matrix in CSR-format.
     *
     * @param i_rowPtr row-pointer of CSR; last element holds the number of non-zero entries.
     * @param i_colIdx column index of CSR.
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
    template <typename T>
    void addSparseSoa( const unsigned int              *i_rowPtr,
                       const unsigned int              *i_colIdx,
                       const T                         *i_val,
                             unsigned int               i_m,
                             unsigned int               i_n,
                             unsigned int               i_k,
                             unsigned int               i_ldA,
                             unsigned int               i_ldB,
                             unsigned int               i_ldC,
                             T                          i_alpha,
                             T                          i_beta,
                             libxsmm_gemm_prefetch_type i_prefetch = LIBXSMM_PREFETCH_NONE ) {
      EDGE_VLOG(1) << "  adding XSMM-kernel #" << m_kernels.size() << " (sparse)"
                   << " M=" << i_m << " N=" << i_n << " K=" << i_k
                   << " ldA=" << i_ldA << " ldB=" << i_ldB << " ldC=" << i_ldC
                   << " alpha=" << i_alpha << " beta=" << i_beta;

      // check precision
#if (PP_PRECISION == 32) and defined(__AVX512F__)
      EDGE_CHECK( sizeof(T)== 4 );
      EDGE_CHECK( PP_N_CRUNS == 16 );
#elif (PP_PRECISION == 32) and defined(__AVX2__)
      EDGE_CHECK( sizeof(T)== 4 );
      EDGE_CHECK( PP_N_CRUNS == 8 );
#elif (PP_PRECISION == 64) and defined(__AVX512F__)
      EDGE_CHECK( sizeof(T)== 8 );
      EDGE_CHECK( PP_N_CRUNS == 8  );
#elif (PP_PRECISION == 64) and defined(__AVX2__)
      EDGE_CHECK( sizeof(T)== 8 );
      EDGE_CHECK( PP_N_CRUNS == 4  );
#else
#error precision not supported.
#endif

      // add description
      libxsmm_gemm_descriptor l_desc;
      LIBXSMM_GEMM_DESCRIPTOR( l_desc, (PP_PRECISION == 64) ? LIBXSMM_GEMM_PRECISION_F64 : LIBXSMM_GEMM_PRECISION_F32, 0,
                               i_m, i_n, i_k, i_ldA, i_ldB, i_ldC,
                               i_alpha, i_beta, i_prefetch );

      m_descs.push_back( l_desc );

      // generate and store function for this kernels
#if PP_PRECISION == 64
      m_kernels.push_back( libxsmm_create_xcsr_soa( &m_descs.back(), i_rowPtr, i_colIdx, i_val ).dmm );
#else
      m_kernels.push_back( libxsmm_create_xcsr_soa( &m_descs.back(), i_rowPtr, i_colIdx, i_val ).smm );
#endif

      // check that we generated a kernel
      EDGE_CHECK( m_kernels.back() != 0 );
    }

    /**
     * Adds a libxsmm dense GEMM kernel, LIBXSMM is col-major and so is this call for row-major usage
     * please flip A and B
     *
     * @param i_m number of rows in column-major A and C
     * @param i_n number of columns in column-major B and C
     * @param i_k number of columns/rows in column-major A/B
     * @param i_ldA leading dimension of column-major A
     * @param i_ldB leading dimension of column-major B
     * @param i_ldC leading dimension of column-major C
     * @param i_alpha alpha parameter (needs to be 1.0 for now)
     * @param i_beta beta parameter (need to be 0.0/1.0 for now)
     * @param i_prefetch prefetch strategy.
     **/
    template <typename T>
    void add_gemm( unsigned int               i_m,
                   unsigned int               i_n,
                   unsigned int               i_k,
                   unsigned int               i_ldA,
                   unsigned int               i_ldB,
                   unsigned int               i_ldC,
                   T                          i_alpha,
                   T                          i_beta,
                   libxsmm_gemm_prefetch_type i_prefetch = LIBXSMM_PREFETCH_NONE ) {
      EDGE_VLOG(1) << "  adding XSMM-kernel gemm #" << m_kernels.size()
                   << " M=" << i_m << " N=" << i_n << " K=" << i_k
                   << " ldA=" << i_ldA << " ldB=" << i_ldB << " ldC=" << i_ldC
                   << " alpha=" << i_alpha << " beta=" << i_beta;

      // check fused runs
      EDGE_CHECK( PP_N_CRUNS == 1 );

      // check precision
#if PP_PRECISION == 32
      EDGE_CHECK( sizeof(T)== 4 );

#elif PP_PRECISION == 64
      EDGE_CHECK( sizeof(T)== 8 );
#else
#error precision not supported.
#endif

      // add description
      libxsmm_gemm_descriptor l_desc;
      LIBXSMM_GEMM_DESCRIPTOR( l_desc, (PP_PRECISION == 64) ? LIBXSMM_GEMM_PRECISION_F64 : LIBXSMM_GEMM_PRECISION_F32, 0,
                               i_m, i_n, i_k, i_ldA, i_ldB, i_ldC,
                               i_alpha, i_beta, i_prefetch );
      m_descs.push_back( l_desc );
      
      // generate and store function for this kernels
#if PP_PRECISION == 64
      m_kernels.push_back( libxsmm_xmmdispatch( &m_descs.back() ).dmm );
#else
      m_kernels.push_back( libxsmm_xmmdispatch( &m_descs.back() ).smm );
#endif

      // check that we generated a kernel
      EDGE_CHECK( m_kernels.back() != 0 );
    }   
};

#endif
