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

#ifndef EDGE_DATA_UNARY_XSMM_HPP
#define EDGE_DATA_UNARY_XSMM_HPP
 
#include <vector>
#include "constants.hpp"
#include "io/logging.h"
#include "XsmmUtils.hpp"
 
#include <libxsmm.h>
 
namespace edge {
  namespace data {
    template< typename TL_T_REAL >
    class UnaryXsmm;
  }
}
 
/**
 * Holds LIBXSMM gemm kernels for non-fused, single precision simulations.
 **/
template<typename TL_T_REAL>
class edge::data::UnaryXsmm {
  private:
    //! generated unary kernels of libxsmm
    std::vector< std::vector< libxsmm_meltwfunction_unary > > u_kernels;

  public:

    /**
     * Adds a libxsmm unary kernel
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
              unsigned int               i_ldi,
              unsigned int               i_ldo,
              libxsmm_meltw_unary_type   type,
              libxsmm_bitfield           flags) {
      EDGE_VLOG(1) << "  adding, X precision XSMM-kernel unary #" << u_kernels.size()
                   << " M=" << i_m << " N=" << i_n
                   << " ldi=" << i_ldi << " ldo=" << i_ldo
                   << " type=" << type << " flags=" << flags;

      // add kernel groups, if required
      if( i_group >= u_kernels.size() ) {
        u_kernels.resize( i_group+1 );
      }

      libxsmm_datatype dtype_in   = XsmmDtype<TL_T_REAL>();
      libxsmm_datatype dtype_out  = XsmmDtype<TL_T_REAL>();
      libxsmm_datatype dtype_comp = XsmmDtype<TL_T_REAL>();

      EDGE_VLOG(1) << "  dtype = " << dtype_in;

      libxsmm_meltw_unary_shape unary_shape = libxsmm_create_meltw_unary_shape(i_m, i_n, i_ldi, i_ldo, dtype_in, dtype_out, dtype_comp);

      // generate and store function for this kernels
      u_kernels[i_group].push_back(libxsmm_dispatch_meltw_unary_v2(type, unary_shape, flags));

      // check that we generated a kernel
      EDGE_CHECK_NE( u_kernels[i_group].back(), 0 );
    }

    void add( unsigned short             i_group,
              unsigned int               i_m,
              unsigned int               i_n,
              libxsmm_meltw_unary_type   type,
              libxsmm_bitfield           flags) {
      add(i_group, i_m, i_n, i_m, i_m, type, flags);
    }

    /* @FIXME: Could have added more of different overloads for various unary kernel calls as in PT extensions */
    void execute( const unsigned short i_group,
                  const unsigned short i_entry,
                        TL_T_REAL*     io_io) const {

      libxsmm_meltw_unary_param unary_param;
      memset( &unary_param, 0, sizeof(libxsmm_meltw_unary_param) );
      unary_param.out.primary   = (void*)io_io;

      u_kernels[i_group][i_entry]( &unary_param );
    }

    void execute( const unsigned short i_group,
                  const unsigned short i_entry,
                  const TL_T_REAL*     i_in,
                        TL_T_REAL*     io_o) const {

      libxsmm_meltw_unary_param unary_param;
      memset( &unary_param, 0, sizeof(libxsmm_meltw_unary_param) );

      unary_param.in.primary    = (void*)i_in;
      unary_param.out.primary   = (void*)io_o;

      u_kernels[i_group][i_entry]( &unary_param );
    }

    void execute( const unsigned short i_group,
                  const unsigned short i_entry,
                  const TL_T_REAL*     i_in,
                        TL_T_REAL*     io_o0,
                        TL_T_REAL*     io_o1) const {

      libxsmm_meltw_unary_param unary_param;
      memset( &unary_param, 0, sizeof(libxsmm_meltw_unary_param) );

      unary_param.in.primary    = (void*)i_in;
      unary_param.out.primary   = (void*)io_o0;
      unary_param.out.secondary = (void*)io_o1;

      u_kernels[i_group][i_entry]( &unary_param );
    }

    void execute( const unsigned short i_group,
                  const unsigned short i_entry,
                  const TL_T_REAL*     i_in0,
                  const TL_T_REAL*     i_in1,
                  const TL_T_REAL*     i_in2,
                  const TL_T_REAL*     i_op0,
                  const TL_T_REAL*     i_op1,
                  const TL_T_REAL*     i_op2,
                        TL_T_REAL*     io_o0,
                        TL_T_REAL*     io_o1) const {

      libxsmm_meltw_unary_param unary_param;
      memset( &unary_param, 0, sizeof(libxsmm_meltw_unary_param) );

      unary_param.in.primary    = (void*)i_in0;
      unary_param.in.secondary  = (void*)i_in1;
      unary_param.in.tertiary   = (void*)i_in2;
      unary_param.op.primary    = (void*)i_op0;
      unary_param.op.secondary  = (void*)i_op1;
      unary_param.op.tertiary   = (void*)i_op2;
      unary_param.out.primary   = (void*)io_o0;
      unary_param.out.secondary = (void*)io_o1;

      u_kernels[i_group][i_entry]( &unary_param );
    }
};
#endif
 
