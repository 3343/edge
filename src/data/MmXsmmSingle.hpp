/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *         Alexander Heinecke (alexander.heinecke AT intel.com)
 *
 * @section LICENSE
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
 * Data structures of the non-fused LIBXSMM, matrix-matrix multiplication kernels.
 **/

#ifndef EDGE_DATA_MM_XSMM_SINGLE_HPP
#define EDGE_DATA_MM_XSMM_SINGLE_HPP
 
#include <vector>
#include "constants.hpp"
#include "io/logging.h"
 
#include "XsmmUtils.hpp"

#include <libxsmm.h>
 
namespace edge {
  namespace data {
    template< typename TL_T_REAL >
    class MmXsmmSingle;
  }
}
 
/**
 * Holds LIBXSMM gemm kernels for non-fused, single precision simulations.
 **/
template< typename TL_T_REAL >
class edge::data::MmXsmmSingle {
  private:
    //! generated kernels of libxsmm
    std::vector< std::vector< libxsmm_gemmfunction > > m_kernels;
  
  public:

    /**
     * Adds a libxsmm dense GEMM kernel
     * Remark: LIBXSMM is col-major and so is this call,
     * for row-major usage please flip A and B
     *
     * @param i_group id of the kernel group.
     * @param i_m number of rows in column-major A and C
     * @param i_n number of columns in column-major B and C
     * @param i_k number of columns/rows in column-major A/B
     * @param i_ldA leading dimension of column-major A
     * @param i_ldB leading dimension of column-major B
     * @param i_ldC leading dimension of column-major C
     * @param i_alpha alpha parameter (needs to be 1.0 for now)
     * @param i_beta beta parameter (need to be 0.0/1.0 for now)
     * @param i_prefetch prefetching strategy.
     *
     **/
    void add( unsigned short             i_group,
              unsigned int               i_m,
              unsigned int               i_n,
              unsigned int               i_k,
              unsigned int               i_ldA,
              unsigned int               i_ldB,
              unsigned int               i_ldC,
              TL_T_REAL                  i_alpha,
              TL_T_REAL                  i_beta,
              libxsmm_gemm_prefetch_type i_prefetch ) {
      EDGE_VLOG(1) << "  adding, X precision XSMM-kernel gemm #" << m_kernels.size()
                   << " M=" << i_m << " N=" << i_n << " K=" << i_k
                   << " ldA=" << i_ldA << " ldB=" << i_ldB << " ldC=" << i_ldC
                   << " alpha=" << i_alpha << " beta=" << i_beta;
 
      // add kernel groups, if required
      if( i_group >= m_kernels.size() ) {
        m_kernels.resize( i_group+1 );
      }

      int l_flags = LIBXSMM_GEMM_FLAGS('N', 'N') | LIBXSMM_GEMM_FLAG_USE_XGEMM_ABI;
      if (i_beta == 0.0)
        l_flags |= LIBXSMM_GEMM_FLAG_BETA_0;

      libxsmm_datatype l_dtype_A    = XsmmDtype<TL_T_REAL>();
      libxsmm_datatype l_dtype_B    = XsmmDtype<TL_T_REAL>();
      libxsmm_datatype l_dtype_C    = XsmmDtype<TL_T_REAL>();
      libxsmm_datatype l_dtype_comp = XsmmDtype<TL_T_REAL>();

      libxsmm_gemm_shape l_gemm_shape = libxsmm_create_gemm_shape(i_m, i_n, i_k, i_ldA, i_ldB, i_ldC,
                                                                  l_dtype_A, l_dtype_B, l_dtype_C, l_dtype_comp);

      // generate and store function for this kernels
      m_kernels[i_group].push_back( libxsmm_dispatch_gemm_v2(l_gemm_shape, l_flags, i_prefetch) );

      // check that we generated a kernel
      EDGE_CHECK_NE( m_kernels[i_group].back(), 0 );
    }

    void execute( const unsigned short i_group,
                  const unsigned short i_entry,
                  const TL_T_REAL*     i_a,
                  const TL_T_REAL*     i_b,
                        TL_T_REAL*     io_c,
                  const TL_T_REAL*     i_pf_a,
                  const TL_T_REAL*     i_pf_b,
                  const TL_T_REAL*     i_pf_c ) const {
      // @TODO we can optimize by putting this into the constructor (if we avoid data races)
      libxsmm_gemm_param gemm_param;
      memset( &gemm_param, 0, sizeof(libxsmm_gemm_param) );

      gemm_param.a.primary    = (void*)i_a;
      gemm_param.b.primary    = (void*)i_b;
      gemm_param.c.primary    = (void*)io_c;
      gemm_param.a.quaternary = (void*)i_pf_a;
      gemm_param.b.quaternary = (void*)i_pf_b;
      gemm_param.c.quaternary = (void*)i_pf_c;

      m_kernels[i_group][i_entry]( &gemm_param );
    }
};
#endif

