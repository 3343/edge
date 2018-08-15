/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *         Alexander Heinecke (alexander.heinecke AT intel.com)
 *
 * @section LICENSE
 * Copyright (c) 2017-2018, Regents of the University of California
 * Copyright (c) 2018, Intel Corporation
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
 * Data structures of the vanilla matrix-matrix multiplication kernels.
 **/

#ifndef EDGE_DATA_MM_VANILLA_HPP
#define EDGE_DATA_MM_VANILLA_HPP

#include <vector>
#include "constants.hpp"
#include "io/logging.h"

#include "linalg/Matrix.h"

namespace edge {
  namespace data {
    template< typename TL_T_REAL >
    class MmVanilla;
  }
}

/**
 * Holds vanilla kernels for fused and nonfused simulations.
 *
 * @paramt TL_T_REAL floating point precision.
 **/
template< typename TL_T_REAL >
class edge::data::MmVanilla {
  private:
    /**
     * Vanilla matrix kernels which store the BLAS identifiers and
     * offers an overloaded function call operator.
     **/
    class Vanilla {
      private:
        //! number of rows in column-major A and C
        const unsigned int m_m;

        //! number of columns in column-major B and C
        const unsigned int m_n;

        //! number of columns/rows in column-major A/B
        const unsigned int m_k;

        //! leading dimension of column-major A
        const unsigned int m_ldA;

        //! leading dimension of column-major B
        const unsigned int m_ldB;

        //! leading dimension of column-major C
        const unsigned int m_ldC;

        //! beta parameter
        const TL_T_REAL m_beta;

        //! true if matrices A and C carry fused simulations
        const bool m_fusedAC;

        //! true if matrices B and C carry fused simulations
        const bool m_fusedBC;

        //! number of fused runs
        const unsigned short m_nCrs;

      public:
        /**
         * Constructor.
         *
         * @param i_m number of rows in column-major A and C.
         * @param i_n number of columns in column-major B and C.
         * @param i_k number of columns/rows in column-major A/B.
         * @param i_ldA leading dimension of column-major A.
         * @param i_ldB leading dimension of column-major B.
         * @param i_ldC leading dimension of column-major C.
         * @param i_beta beta parameter (needs to be 0 or 1 for now).
         * @param i_fusedAC true if matrices A and C are fused.
         * @param i_fusedBC true if matrices B and C are fused.
         * @param i_nCrs number of fused runs.
         **/
        Vanilla( unsigned int   i_m,
                 unsigned int   i_n,
                 unsigned int   i_k,
                 unsigned int   i_ldA,
                 unsigned int   i_ldB,
                 unsigned int   i_ldC,
                 TL_T_REAL      i_beta,
                 bool           i_fusedAC,
                 bool           i_fusedBC,
                 unsigned short i_nCrs ): m_m( i_m ),
                                          m_n( i_n ),
                                          m_k( i_k ),
                                          m_ldA( i_ldA ),
                                          m_ldB( i_ldB ),
                                          m_ldC( i_ldC ),
                                          m_beta( i_beta ),
                                          m_fusedAC( i_fusedAC ),
                                          m_fusedBC( i_fusedBC ),
                                          m_nCrs( i_nCrs ) {};

          /**
           * Overloads the function call operator with the kernel execution.
           *
           * @param i_a matrix A.
           * @param i_b matrix B.
           * @param io_c matrix C.
           **/
          void operator()( TL_T_REAL const * i_a,
                           TL_T_REAL const * i_b,
                           TL_T_REAL       * io_c ) const {
            if( m_fusedAC ) {
              linalg::Matrix::matMulFusedAC( m_nCrs,
                                             m_m,   m_n,   m_k,
                                             m_ldA, m_ldB, m_ldC,
                                             m_beta,
                                             i_a,   i_b,   io_c );
            }
            else if( m_fusedBC ) {
              linalg::Matrix::matMulFusedBC( m_nCrs,
                                             m_m,   m_n,   m_k,
                                             m_ldA, m_ldB, m_ldC,
                                             m_beta,
                                             i_a,   i_b,   io_c );
            }
            else EDGE_LOG_FATAL << "matrix structure not supported";
          };
    };

  public:
    //! vanilla matrix kernels
    std::vector< std::vector< Vanilla > > m_kernels;

    /**
     * Adds a vanilla kernels (either fused or non-fused).
     * If the given kernel group does not exist, new groups until the given id are created.
     *
     * @param i_group id of the kernel group.
     * @param i_m number of rows in column-major A and C.
     * @param i_n number of columns in column-major B and C.
     * @param i_k number of columns/rows in column-major A/B.
     * @param i_ldA leading dimension of column-major A.
     * @param i_ldB leading dimension of column-major B.
     * @param i_ldC leading dimension of column-major C.
     * @param i_alpha parameter alpha, ignored.
     * @param i_beta parameter beta, needs to be TL_T_REAL(0) or TL_T_REAL(1) for now.
     * @param i_fusedAC true if matrices A and C are fused.
     * @param i_fusedBC true if matrices B and C are fused.
     * @param i_nCfr number of fused simulations.
     **/
    void add( unsigned short i_group,
              unsigned int   i_m,
              unsigned int   i_n,
              unsigned int   i_k,
              unsigned int   i_ldA,
              unsigned int   i_ldB,
              unsigned int   i_ldC,
              TL_T_REAL      i_alpha,
              TL_T_REAL      i_beta,
              bool           i_fusedAC,
              bool           i_fusedBC,
              unsigned short i_nCfr ) {
      // verbose output
      EDGE_VLOG(1) << "  adding vanilla-kernel #" << m_kernels.size() << " (dense)"
                   << " M=" << i_m << " N=" << i_n << " K=" << i_k
                   << " ldA=" << i_ldA << " ldB=" << i_ldB << " ldC=" << i_ldC
                   << " alpha=" << i_alpha << " beta=" << i_beta
                   << " fusedAC=" << i_fusedAC << " fusedBC=" << i_fusedBC
                   << " cfr=" << i_nCfr;

      EDGE_CHECK( i_nCfr == 1 || (i_fusedAC || i_fusedBC) );

      // add kernel groups, if required
      if( i_group >= m_kernels.size() ) m_kernels.resize( i_group+1 );

      // add kernel
      m_kernels[i_group].push_back( Vanilla( i_m, i_n, i_k,
                                             i_ldA, i_ldB, i_ldC,
                                             i_beta,
                                             i_fusedAC,
                                             i_fusedBC,
                                             i_nCfr ) );
    }
};
#endif
