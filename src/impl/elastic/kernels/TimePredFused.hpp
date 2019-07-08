/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *         Alexander Heinecke (alexander.heinecke AT intel.com)
 *
 * @section LICENSE
 * Copyright (c) 2019, Alexander Breuer
 * Copyright (c) 2016-2019, Regents of the University of California
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
 * Time predictions through the ADER scheme for seismic setups with fused forward simulations.
 **/
#ifndef EDGE_SEISMIC_KERNELS_TIME_PRED_FUSED_HPP
#define EDGE_SEISMIC_KERNELS_TIME_PRED_FUSED_HPP

#include "TimePred.hpp"
#include "data/MmXsmmFused.hpp"
#include "impl/elastic/common.hpp"

namespace edge {
  namespace elastic {
    namespace kernels {
      template< typename       TL_T_REAL,
                t_entityType   TL_T_EL,
                unsigned short TL_O_SP,
                unsigned short TL_O_TI,
                unsigned short TL_N_CRS >
      class TimePredFused;
    }
  }
}

/**
 * Optimized ADER time prediction for fused seismic forward simulations.
 *
 * @paramt TL_T_REAL floating point precision.
 * @paramt TL_T_EL element type.
 * @paramt TL_O_SP order in space.
 * @paramt TL_O_TI order in time.
 * @paramt TL_N_CRS number of fused simulations. 
 **/
template< typename       TL_T_REAL,
          t_entityType   TL_T_EL,
          unsigned short TL_O_SP,
          unsigned short TL_O_TI,
          unsigned short TL_N_CRS >
class edge::elastic::kernels::TimePredFused: public edge::elastic::kernels::TimePred < TL_T_REAL,
                                                                                       TL_T_EL,
                                                                                       TL_O_SP,
                                                                                       TL_O_TI,
                                                                                       TL_N_CRS > {
  private:
    //! dimension of the element
    static unsigned short const TL_N_DIS = C_ENT[TL_T_EL].N_DIM;

    //! number of element modes
    static unsigned short const TL_N_MDS = CE_N_ELEMENT_MODES( TL_T_EL, TL_O_SP );

    //! number of elastic quantities
    static unsigned short const TL_N_QTS_E = (TL_N_DIS == 2) ? 5 : 9;

    //! number of non-zeroes in the star matrices
    static unsigned short const TL_N_NZS_STAR_E = (TL_N_DIS==2) ? 10 : 24;

    //! pointers to the (possibly recursive) stiffness matrices
    TL_T_REAL *m_stiffT[CE_MAX(TL_O_SP-1,1)][TL_N_DIS];

    //! matrix kernels
    edge::data::MmXsmmFused< TL_T_REAL > m_mm;

    /**
     * Sets the given matrix to zero.
     *
     * @param o_mat matrix which will be set to 0.
     **/
    static void inline zero( TL_T_REAL o_mat[TL_N_QTS_E][TL_N_MDS][TL_N_CRS] ) {
      // reset result to zero
      for( unsigned short l_qt = 0; l_qt < TL_N_QTS_E; l_qt++ ) {
        for( unsigned short l_md = 0; l_md < TL_N_MDS; l_md++ ) {
          for( unsigned short l_cfr = 0; l_cfr < TL_N_CRS; l_cfr++ ) {
            o_mat[l_qt][l_md][l_cfr] = 0;
          }
        }
      }
    }

    /**
     * Gets the recursive sparse matrix structures of the dense input matrices.
     * 
     * @param i_stiffT transposed stiffness matrices.
     * @param o_maxNzCols will be set to maximum non-zero column of the three stiffness matrices for every recursion.
     * @param o_offsets will be set to the offsets (counting non-zero entries) of the sparse matrices.
     * @param o_nonZeros will be set to the raw non-zero entries of the sparse matrices.
     * @param o_mats will be set to the CSC representation of the sparse matrices.
     **/
    static void getCscStiffT( TL_T_REAL               const   i_stiffT[TL_N_DIS][TL_N_MDS][TL_N_MDS],
                              std::vector< size_t >         & o_maxNzCols,
                              std::vector< size_t >         & o_offsets,
                              std::vector< TL_T_REAL >      & o_nonZeros,
                              std::vector< t_matCsc >       & o_mats ) {
      // determine CSC fill-in strategy
      std::string l_cscFillIn = "none";

      if ( libxsmm_get_target_archid() == LIBXSMM_X86_AVX512_KNM ) {
        l_cscFillIn = "qfma";
      }

      // reset and init output
      o_maxNzCols.resize( 0 );
      o_offsets.resize( 0 );
      o_offsets.push_back( 0 );
      o_nonZeros.resize( 0 );
      o_mats.resize( 0 );

      // derive the non-zeros of the recursive ADER scheme.
      t_matCrd l_stiffTCrd[TL_N_DIS];
      for( unsigned short l_di = 0; l_di < TL_N_DIS; l_di++ ) {
        edge::linalg::Matrix::denseToCrd< TL_T_REAL >( TL_N_MDS,
                                                       TL_N_MDS,
                                                       i_stiffT[l_di][0],
                                                       l_stiffTCrd[l_di],
                                                       TOL.BASIS );
      }

      // nz-blocks
      unsigned int l_nzBl[2][2][2];
      // init with matrix dim
      l_nzBl[0][0][0] = l_nzBl[0][1][0] = 0;
      l_nzBl[0][0][1] = l_nzBl[0][1][1] = TL_N_MDS-1;

      for( unsigned short l_de = 1; l_de < TL_O_SP; l_de++ ) {
        // determine non-zero block for the next iteration
        unsigned int l_maxNzCol = 0;

        for( unsigned short l_di = 0; l_di < TL_N_DIS; l_di++ ) {
          edge::linalg::Matrix::getBlockNz( l_stiffTCrd[l_di], l_nzBl[0], l_nzBl[1] );
          l_maxNzCol = std::max( l_maxNzCol, l_nzBl[1][1][1] );
        }
        o_maxNzCols.push_back( l_maxNzCol );

        t_matCsc l_stiffTCsc[TL_N_DIS];
        for( unsigned short l_di = 0; l_di < TL_N_DIS; l_di++ ) {
          edge::linalg::Matrix::denseToCsc< TL_T_REAL >( TL_N_MDS,
                                                         TL_N_MDS,
                                                         i_stiffT[l_di][0],
                                                         l_stiffTCsc[l_di],
                                                         TOL.BASIS,
                                                         l_nzBl[0][0][1]+1,
                                                         l_maxNzCol+1,
                                                         l_cscFillIn );
          o_mats.push_back( l_stiffTCsc[l_di] );
        }

        // add data for the first CK-stiff matrix or shrinking stiff matrices
        if( l_de == 1 || l_maxNzCol < l_nzBl[0][0][1] ) {
          for( unsigned short l_di = 0; l_di < TL_N_DIS; l_di++ ) {
            for( unsigned int l_nz = 0; l_nz < l_stiffTCsc[l_di].val.size(); l_nz++ ) {
              o_nonZeros.push_back( l_stiffTCsc[l_di].val[l_nz] );
            }
            o_offsets.push_back( o_offsets.back() + l_stiffTCsc[l_di].val.size() );
          }
        }
        else {
          for( unsigned short l_di = 0; l_di < TL_N_DIS; l_di++ ) {
            o_offsets.push_back( o_offsets.back() );
          }
        }

#ifdef PP_T_BASIS_HIERARCHICAL
      // check that size goes down with the number of derivatives
      EDGE_CHECK_EQ( l_maxNzCol+1, CE_N_ELEMENT_MODES( TL_T_EL, TL_O_SP-l_de ) );
#endif

        // reduce relevant rows due to generated zero block
        l_nzBl[0][0][1] = l_maxNzCol;
      }
    }

    /**
     * Generates the matrix kernels for the transposed stiffness matrices and star matrices.
     * 
     * @param i_stiffT dense representation of the transposed matrices.
     **/
    void generateKernels( TL_T_REAL const i_stiffT[TL_N_DIS][TL_N_MDS][TL_N_MDS] ) {
      // convert stiffness matrices to CSC (incl. possible fill-in)
      std::vector< size_t > l_maxNzCols;
      std::vector< size_t > l_offsets;
      std::vector< TL_T_REAL > l_nonZeros;
      std::vector< t_matCsc > l_cscStiffT;
      getCscStiffT( i_stiffT,
                    l_maxNzCols,
                    l_offsets,
                    l_nonZeros,
                    l_cscStiffT );

      // assemble fake star matrix
      TL_T_REAL l_star[TL_N_DIS][TL_N_QTS_E][TL_N_QTS_E];
      edge::elastic::common::getJac( (TL_T_REAL) 1.0,
                                     (TL_T_REAL) 1.0,
                                     (TL_T_REAL) 1.0,
                                                 l_star[0][0],
                                                 TL_N_DIS );
      for( unsigned short l_di = 1; l_di < TL_N_DIS; l_di++ ) {
        for( int_qt l_q1 = 0; l_q1 < TL_N_QTS_E; l_q1++ ) {
          for( int_qt l_q2 = 0; l_q2 < TL_N_QTS_E; l_q2++ ) {
            l_star[0][l_q1][l_q2] =   std::abs(l_star[0][l_q1][l_q2])
                                    + std::abs(l_star[l_di][l_q1][l_q2]);
          }
        }
      }

      // get csr star matrix
      t_matCsr l_starCsr;
      edge::linalg::Matrix::denseToCsr< TL_T_REAL >( TL_N_QTS_E,
                                                     TL_N_QTS_E,
                                                     l_star[0][0],
                                                     l_starCsr,
                                                     TOL.BASIS );

      // derive sparse AoSoA-LIBXSMM kernels
      for( unsigned short l_de = 1; l_de < TL_O_SP; l_de++ ) {
        // generate libxsmm kernel for transposed stiffness matrices
        for( unsigned short l_di = 0; l_di < TL_N_DIS; l_di++ ) {
          m_mm.add( 0,
                    false,
                    &l_cscStiffT[(l_de-1)*TL_N_DIS + l_di].colPtr[0], 
                    &l_cscStiffT[(l_de-1)*TL_N_DIS + l_di].rowIdx[0],
                    &l_cscStiffT[(l_de-1)*TL_N_DIS + l_di].val[0],
                    TL_N_QTS_E,
                    l_maxNzCols[l_de-1]+1,
                    (l_de == 1) ? TL_N_MDS : l_maxNzCols[l_de-2]+1,
                    TL_N_MDS,
                    0,
                    l_maxNzCols[l_de-1]+1,
                    TL_T_REAL(1.0),
                    TL_T_REAL(0.0),
                    LIBXSMM_GEMM_PREFETCH_NONE );
        }
        // generate libxsmm kernel for star matrix
        m_mm.add( 1,
                  true,
                  &l_starCsr.rowPtr[0],
                  &l_starCsr.colIdx[0],
                  &l_starCsr.val[0],
                  TL_N_QTS_E,
                  l_maxNzCols[l_de-1]+1,
                  TL_N_QTS_E,
                  0,
                  l_maxNzCols[l_de-1]+1,
                  TL_N_MDS,
                  TL_T_REAL(1.0),
                  TL_T_REAL(1.0),
                  LIBXSMM_GEMM_PREFETCH_NONE );
      }
    }

    /**
     * Stores the transposed stiffness matrices.
     * 
     * @param i_stiffT dense stiffness matrices.
     * @param io_dynMem dynamic memory management, which will be used for the respective allocations.
     * @param o_stiffT will contain pointers to memory for the individual matrices.
     **/
    static void storeStiffTSparse( TL_T_REAL     const     i_stiffT[TL_N_DIS][TL_N_MDS][TL_N_MDS],
                                   data::Dynamic       &   io_dynMem,
                                   TL_T_REAL           *   o_stiffT[CE_MAX(TL_O_SP-1,1)][TL_N_DIS]  ) {
      // convert stiffness matrices to CSC (incl. possible fill-in)
      std::vector< size_t > l_maxNzCols;
      std::vector< size_t > l_offsets;
      std::vector< TL_T_REAL > l_nonZeros;
      std::vector< t_matCsc > l_cscStiffT;
      getCscStiffT( i_stiffT,
                    l_maxNzCols,
                    l_offsets,
                    l_nonZeros,
                    l_cscStiffT );

      // copy sparse matrices to a permanent data structure
      TL_T_REAL * l_stiffTRaw = (TL_T_REAL*) io_dynMem.allocate( l_nonZeros.size() * sizeof(TL_T_REAL),
                                                                 4096,
                                                                 true );
      for( std::size_t l_en = 0; l_en < l_nonZeros.size(); l_en++ ) {
        l_stiffTRaw[l_en] = l_nonZeros[l_en];
      }

      // check that we have offsets for all matrices
      EDGE_CHECK_EQ( l_offsets.size(), (TL_O_SP-1)*TL_N_DIS+1 );

      // set the pointers
      unsigned short l_mat = 0;
      for( unsigned int l_de = 1; l_de < ORDER; l_de++ ) {
        for( unsigned int l_di = 0; l_di < N_DIM; l_di++ ) {
          o_stiffT[(l_de-1)][l_di] = l_stiffTRaw + l_offsets[l_mat];
          l_mat++;
        }
      }
    }

  public:
    /**
     * Constructor of the fused time prediction.
     *
     * @param io_dynMem dynamic memory allocations.
     **/
    TimePredFused( data::Dynamic & io_dynMem ) {
      // formulation of the basis in terms of the reference element
      dg::Basis l_basis( TL_T_EL,
                         TL_O_SP );

      // get stiffness matrices
      TL_T_REAL l_stiffT[TL_N_DIS][TL_N_MDS][TL_N_MDS];
      l_basis.getStiffMm1Dense( TL_N_MDS,
                                l_stiffT[0][0],
                                true );


      // store stiffness matrices dense
      this->storeStiffTSparse( l_stiffT,
                               io_dynMem,
                               m_stiffT );

      // generate kernels
      generateKernels( l_stiffT );
    };

    /**
     * Applies the Cauchy–Kowalevski procedure (fused LIBXSMM version) and computes time derivatives and time integrated DOFs.
     *
     * @param i_dT time step.
     * @param i_star star matrices.
     * @param i_dofs DOFs.
     * @param o_scratch will be used as scratch memory.
     * @param o_der will be set to time derivatives.
     * @param o_tInt will be set to time integrated DOFs.
     **/
    void ck( TL_T_REAL       i_dT,
             TL_T_REAL const i_star[TL_N_DIS][TL_N_NZS_STAR_E],
             TL_T_REAL const i_dofs[TL_N_QTS_E][TL_N_MDS][TL_N_CRS],
             TL_T_REAL       o_scratch[TL_N_QTS_E][TL_N_MDS][TL_N_CRS],
             TL_T_REAL       o_der[TL_O_TI][TL_N_QTS_E][TL_N_MDS][TL_N_CRS],
             TL_T_REAL       o_tInt[TL_N_QTS_E][TL_N_MDS][TL_N_CRS] ) const {
      // scalar for the time integration
      TL_T_REAL l_scalar = i_dT;

      // initialize zero-derivative, reset time integrated dofs
      for( int_qt l_qt = 0; l_qt < TL_N_QTS_E; l_qt++ ) {
        for( int_md l_md = 0; l_md < TL_N_MDS; l_md++ ) {
#pragma omp simd
          for( int_cfr l_cfr = 0; l_cfr < TL_N_CRS; l_cfr++ ) {
            o_der[0][l_qt][l_md][l_cfr] = i_dofs[l_qt][l_md][l_cfr];
            o_tInt[l_qt][l_md][l_cfr]   = l_scalar * i_dofs[l_qt][l_md][l_cfr];
          }
        }
      }

      // iterate over time derivatives
      for( unsigned int l_de = 1; l_de < TL_O_TI; l_de++ ) {
        // reset this derivative
        for( int_qt l_qt = 0; l_qt < TL_N_QTS_E; l_qt++ )
          for( int_md l_md = 0; l_md < TL_N_MDS; l_md++ )
#pragma omp simd
            for( int_cfr l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) o_der[l_de][l_qt][l_md][l_cr] = 0;

        // compute the derivatives
        for( unsigned short l_di = 0; l_di < TL_N_DIS; l_di++ ) {
          // reset to zero for basis in non-hierarchical storage
          if( TL_T_EL == TET4 || TL_T_EL == TRIA3 ) {}
          else {
            zero( o_scratch );
          }

          // multiply with transposed stiffness matrices and inverse mass matrix
          m_mm.m_kernels[0][(l_de-1)*(TL_N_DIS)+l_di]( o_der[l_de-1][0][0],
                                                       m_stiffT[l_de-1][l_di],
                                                       o_scratch[0][0] );
          // multiply with star matrices
          m_mm.m_kernels[1][l_de-1]( i_star[l_di],
                                     o_scratch[0][0],
                                     o_der[l_de][0][0] );
        }

        // update scalar
        l_scalar *= -i_dT / (l_de+1);

        // update time integrated dofs
        for( int_qt l_qt = 0; l_qt < TL_N_QTS_E; l_qt++ ) {
          for( int_md l_md = 0; l_md < CE_N_ELEMENT_MODES_CK( TL_T_EL, TL_O_SP, l_de ); l_md++ ) {
#pragma omp simd
            for( int_cfr l_cr = 0; l_cr < TL_N_CRS; l_cr++ )
              o_tInt[l_qt][l_md][l_cr] += l_scalar * o_der[l_de][l_qt][l_md][l_cr];
          }
        }
      }
    }
};

#endif