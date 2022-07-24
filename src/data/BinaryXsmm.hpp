/**
 * @file This file is part of EDGE.
 *
 * @author Kirill Voronin (kirill.voronin AT intel.com)
 *
 * @section LICENSE
 * Copyright (c) 2022, Intel Corporation
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
 * Data structures for binary TPP kernels from LIBXSMM.
 **/

#ifndef EDGE_DATA_BINARY_XSMM_HPP
#define EDGE_DATA_BINARY_XSMM_HPP
 
#include <vector>
#include "constants.hpp"
#include "io/logging.h"
#include "XsmmUtils.hpp"
 
#include <libxsmm.h>
 
namespace edge {
  namespace data {
    template< typename TL_T_REAL >
    class BinaryXsmm;
  }
}
 
/**
 * Holds LIBXSMM binary kernels.
 **/
template<typename TL_T_REAL>
class edge::data::BinaryXsmm {
  private:
    //! generated binary kernels of libxsmm
    std::vector< std::vector< libxsmm_meltwfunction_binary > > b_kernels;

  public:

    /**
     * Adds a libxsmm binary kernel
     * Remark: LIBXSMM is col-major and so is this call,
     * for row-major usage please flip the shape
     *
     * @param i_group id of the kernel group.
     * @param i_m number of rows in column-major inputs/output
     * @param i_n number of columns in column-major inputs/output
     * @param i_ldi0 leading dimension of column-major first input
     * @param i_ldi1 leading dimension of column-major second input
     * @param i_ldo leading dimension of column-major output
     * @param i_type type of the binary operation
     * @param i_flags additional operation modifiers.
     *
     **/

    void add( unsigned short             i_group,
              unsigned int               i_m,
              unsigned int               i_n,
              unsigned int               i_ldi0,
              unsigned int               i_ldi1,
              unsigned int               i_ldo,
              libxsmm_meltw_binary_type  i_type,
              libxsmm_bitfield           i_flags) {
      EDGE_VLOG(1) << "  adding, X precision XSMM-kernel binary #" << b_kernels.size()
                   << " M=" << i_m << " N=" << i_n
                   << " ldi0=" << i_ldi0 << " ldi1=" << i_ldi1 << " ldo=" << i_ldo
                   << " type=" << i_type << " flags=" << i_flags;

      // add kernel groups, if required
      if( i_group >= b_kernels.size() ) {
        b_kernels.resize( i_group+1 );
      }

      libxsmm_datatype l_dtype_in0  = XsmmDtype<TL_T_REAL>();
      libxsmm_datatype l_dtype_in1  = XsmmDtype<TL_T_REAL>();
      libxsmm_datatype l_dtype_out  = XsmmDtype<TL_T_REAL>();
      libxsmm_datatype l_dtype_comp = XsmmDtype<TL_T_REAL>();

      libxsmm_meltw_binary_shape l_binary_shape =
          libxsmm_create_meltw_binary_shape(i_m, i_n, i_ldi0, i_ldi1, i_ldo,
                                                l_dtype_in0, l_dtype_in1,
                                                l_dtype_out, l_dtype_comp);

      // generate and store function for this kernels
      b_kernels[i_group].push_back(libxsmm_dispatch_meltw_binary_v2(i_type,
                                                    l_binary_shape, i_flags));

      // check that we generated a kernel
      EDGE_CHECK_NE( b_kernels[i_group].back(), 0 );
    }

    void add( unsigned short             i_group,
              unsigned int               i_m,
              unsigned int               i_n,
              libxsmm_meltw_binary_type  type,
              libxsmm_bitfield           flags) {
      add(i_group, i_m, i_n, i_m, i_m, i_m, type, flags);
    }

    void execute( const unsigned short i_group,
                  const unsigned short i_entry,
                  const TL_T_REAL*     i_in0,
                  const TL_T_REAL*     i_in1,
                        TL_T_REAL*     io_o) const {

      libxsmm_meltw_binary_param binary_param;
      memset( &binary_param, 0, sizeof(libxsmm_meltw_binary_param) );

      binary_param.in0.primary = (void*)i_in0;
      binary_param.in1.primary = (void*)i_in1;
      binary_param.out.primary = (void*)io_o;

      b_kernels[i_group][i_entry]( &binary_param );
    }

};
#endif
 
