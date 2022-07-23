/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (alex.breuer AT uni-jena.de)
 *         Alexander Heinecke (alexander.heinecke AT intel.com)
 *
 * @section LICENSE
 * Copyright (c) 2021, Friedrich Schiller University Jena
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
 * Time predictions through the ADER scheme for seismic setups with a single forward simulation.
 **/
#ifndef EDGE_SEISMIC_KERNELS_TIME_PRED_SINGLE_HPP
#define EDGE_SEISMIC_KERNELS_TIME_PRED_SINGLE_HPP

#include "TimePred.hpp"
#include "data/MmXsmmSingle.hpp"
#include "data/UnaryXsmm.hpp"
#include "data/BinaryXsmm.hpp"
#include "data/TernaryXsmm.hpp"

#define ELTWISE_TPP

namespace edge {
  namespace seismic {
    namespace kernels {
      template< typename       TL_T_REAL,
                unsigned short TL_N_RMS,
                t_entityType   TL_T_EL,
                unsigned short TL_O_SP,
                unsigned short TL_O_TI >
      class TimePredSingle;
    }
  }
}

/**
 * Optimized ADER time predictions for single forward simulations.
 * 
 * @paramt TL_T_REAL floating point precision.
 * @paramt TL_N_RMS number of relaxation mechanisms.
 * @paramt TL_T_EL element type.
 * @paramt TL_O_SP order in space.
 * @paramt TL_O_TI order in time.
 **/
template< typename       TL_T_REAL,
          unsigned short TL_N_RMS,
          t_entityType   TL_T_EL,
          unsigned short TL_O_SP,
          unsigned short TL_O_TI >
class edge::seismic::kernels::TimePredSingle: public edge::seismic::kernels::TimePred < TL_T_REAL,
                                                                                        TL_N_RMS,
                                                                                        TL_T_EL,
                                                                                        TL_O_SP,
                                                                                        TL_O_TI,
                                                                                        1 > {
  private:
    //! dimension of the element
    static unsigned short const TL_N_DIS = C_ENT[TL_T_EL].N_DIM;

    //! number of element modes
    static unsigned short const TL_N_MDS = CE_N_ELEMENT_MODES( TL_T_EL, TL_O_SP );

    //! number of elastic quantities
    static unsigned short const TL_N_QTS_E = CE_N_QTS_E( TL_N_DIS );

    //! number of quantities per relaxation mechanism
    static unsigned short const TL_N_QTS_M = CE_N_QTS_M( TL_N_DIS );

    //! number of non-zeros in the elastic star matrices
    static unsigned short const TL_N_ENS_STAR_E = CE_N_ENS_STAR_E_DE( TL_N_DIS );

    //! number of non-zeros in the anelastic star matrices
    static unsigned short const TL_N_ENS_STAR_A = CE_N_ENS_STAR_A_DE( TL_N_DIS );

    //! number of non-zeros in the anelastic source matrices
    static unsigned short const TL_N_ENS_SRC_A = CE_N_ENS_SRC_A_DE( TL_N_DIS );

    //! pointers to the (possibly recursive) stiffness matrices
    TL_T_REAL *m_stiffT[CE_MAX(TL_O_TI-1,1)][TL_N_DIS] = {};

    //! matrix kernels
    edge::data::MmXsmmSingle< TL_T_REAL > m_mm;

    //! unary kernels
    edge::data::UnaryXsmm< TL_T_REAL > u_unary;

    //! binary kernels
    edge::data::BinaryXsmm< TL_T_REAL > b_binary;

    //! ternary kernels
    edge::data::TernaryXsmm< TL_T_REAL > t_ternary;

    /**
     * Generates the matrix kernels for the transposed stiffness matrices and star matrices.
     **/
    void generateKernels() {
      // (O-1)*2 kernels for time integration
      for( unsigned int l_de = 1; l_de < TL_O_TI; l_de++ ) {
        // transposed stiffness matrix
        m_mm.add( 0,                                                 // group
                  CE_N_ELEMENT_MODES_CK( TL_T_EL, TL_O_SP, l_de ),   // m
                  TL_N_QTS_E,                                        // n
                  CE_N_ELEMENT_MODES_CK( TL_T_EL, TL_O_SP, l_de-1 ), // k
                  TL_N_MDS,                                          // ldA
                  TL_N_MDS,                                          // ldB
                  TL_N_MDS,                                          // ldC
                  static_cast<TL_T_REAL>(1.0),                       // alpha
                  static_cast<TL_T_REAL>(0.0),                       // beta
                  LIBXSMM_GEMM_PREFETCH_NONE );
      
        // star matrix
        m_mm.add( 1,                                               // group
                  CE_N_ELEMENT_MODES_CK( TL_T_EL, TL_O_SP, l_de ), // m
                  TL_N_QTS_E,                                      // n
                  TL_N_QTS_E,                                      // k
                  TL_N_MDS,                                        // ldA
                  TL_N_QTS_E,                                      // ldB
                  TL_N_MDS,                                        // ldC
                  static_cast<TL_T_REAL>(1.0),                     // alpha
                  static_cast<TL_T_REAL>(1.0),                     // beta
                  LIBXSMM_GEMM_PREFETCH_NONE );

        // anelastic kernels
        if( TL_N_RMS > 0 ) {
          // anelastic star matrix
          m_mm.add( 2,                                            // group
                    CE_N_ELEMENT_MODES_CK( TL_T_EL, TL_O_SP, 1 ), // m
                    TL_N_QTS_M,                                   // n
                    TL_N_DIS,                                     // k
                    TL_N_MDS,                                     // ldA
                    TL_N_DIS,                                     // ldB
                    TL_N_MDS,                                     // ldC
                    static_cast<TL_T_REAL>(1.0),                  // alpha
                    static_cast<TL_T_REAL>(1.0),                  // beta
                    LIBXSMM_GEMM_PREFETCH_NONE );

          // anelastic source matrix
          m_mm.add( 2,                           // group
                    TL_N_MDS,                    // m
                    TL_N_QTS_M,                  // n
                    TL_N_QTS_M,                  // k
                    TL_N_MDS,                    // ldA
                    TL_N_QTS_M,                  // ldB
                    TL_N_MDS,                    // ldC
                    static_cast<TL_T_REAL>(1.0), // alpha
                    static_cast<TL_T_REAL>(1.0), // beta
                    LIBXSMM_GEMM_PREFETCH_NONE );
        }
      }

#ifdef ELTWISE_TPP
      // initialize zero-derivative, reset time integrated dofs
      u_unary.add(0, TL_N_MDS, TL_N_QTS_E /* m, n */, TL_N_MDS, TL_N_MDS /* ldi, ldo */, LIBXSMM_MELTW_TYPE_UNARY_IDENTITY, LIBXSMM_MELTW_FLAG_UNARY_NONE);
      b_binary.add(0, TL_N_MDS, TL_N_QTS_E /* m, n */, TL_N_MDS, TL_N_MDS, TL_N_MDS /* ldi0, ldi1, ldo */,
                    LIBXSMM_MELTW_TYPE_BINARY_MUL, LIBXSMM_MELTW_FLAG_BINARY_BCAST_SCALAR_IN_0);

      // anelastic: init zero-derivative, reset tDofs
      u_unary.add(1,  TL_N_MDS * TL_N_QTS_M, TL_N_RMS /* m, n */, TL_N_MDS * TL_N_QTS_M, TL_N_MDS * TL_N_QTS_M * TL_O_TI /*ldi, ldo */, LIBXSMM_MELTW_TYPE_UNARY_IDENTITY, LIBXSMM_MELTW_FLAG_UNARY_NONE);
      b_binary.add(1, TL_N_MDS * TL_N_QTS_M, TL_N_RMS /* m, n */, LIBXSMM_MELTW_TYPE_BINARY_MUL, LIBXSMM_MELTW_FLAG_BINARY_BCAST_SCALAR_IN_0);

      // zeroing: elastic: reset this derivative
      u_unary.add(2, TL_N_MDS, TL_N_QTS_E /* m, n */, LIBXSMM_MELTW_TYPE_UNARY_XOR, LIBXSMM_MELTW_FLAG_UNARY_NONE);

      // zeroing: buffer for the anelastic computations
      u_unary.add(2, TL_N_MDS, TL_N_QTS_M /* m, n */, LIBXSMM_MELTW_TYPE_UNARY_XOR, LIBXSMM_MELTW_FLAG_UNARY_NONE);

      // multiply with relaxation frequency and add
      // addition
      b_binary.add(3, TL_N_MDS, TL_N_QTS_M /* m, n */, TL_N_MDS, TL_N_MDS, TL_N_MDS /* ldi0, ldi1, ldo */,
                    LIBXSMM_MELTW_TYPE_BINARY_ADD, LIBXSMM_MELTW_FLAG_BINARY_NONE);
      // mult + assign
      b_binary.add(3, TL_N_MDS, TL_N_QTS_M /* m, n */, TL_N_MDS, TL_N_MDS, TL_N_MDS /* ldi0, ldi1, ldo */,
                    LIBXSMM_MELTW_TYPE_BINARY_MUL, LIBXSMM_MELTW_FLAG_BINARY_BCAST_SCALAR_IN_0);
      // accumulation 2
      b_binary.add(3, TL_N_MDS, TL_N_QTS_M /* m, n */, TL_N_MDS, TL_N_MDS, TL_N_MDS /* ldi0, ldi1, ldo */,
                    LIBXSMM_MELTW_TYPE_BINARY_MUL, LIBXSMM_MELTW_FLAG_BINARY_BCAST_SCALAR_IN_0);
      b_binary.add(3, TL_N_MDS, TL_N_QTS_M /* m, n */, TL_N_MDS, TL_N_MDS, TL_N_MDS /* ldi0, ldi1, ldo */,
                    LIBXSMM_MELTW_TYPE_BINARY_ADD, LIBXSMM_MELTW_FLAG_BINARY_NONE);


      // update time integrated dofs
      for( unsigned short l_de = 1; l_de < TL_O_TI; l_de++ ) {
        unsigned short l_nCpMds = (TL_N_RMS == 0) ? CE_N_ELEMENT_MODES_CK( TL_T_EL, TL_O_SP, l_de ) : TL_N_MDS;
        //std::cout << "l_de = " << l_de << " l_nCpMds = " << l_nCpMds << " TL_N_MDS = " << TL_N_MDS << std::endl;
        b_binary.add(4, l_nCpMds, TL_N_QTS_E /* m, n */, TL_N_MDS, TL_N_MDS, TL_N_MDS /* ldi0, ldi1, ldo */,
                      LIBXSMM_MELTW_TYPE_BINARY_MUL, LIBXSMM_MELTW_FLAG_BINARY_BCAST_SCALAR_IN_0);
        b_binary.add(5, l_nCpMds, TL_N_QTS_E /* m, n */, TL_N_MDS, TL_N_MDS, TL_N_MDS /* ldi0, ldi1, ldo */,
                      LIBXSMM_MELTW_TYPE_BINARY_ADD, LIBXSMM_MELTW_FLAG_BINARY_NONE);
      }
#endif
    }

  public:
    /**
     * Constructor of the optimized time prediction for single forward simlations.
     *
     * @param i_rfs relaxation frequencies, use nullptr if TL_N_RMS==0.
     * @param io_dynMem dynamic memory allocations.
     **/
    TimePredSingle( TL_T_REAL     const * i_rfs,
                    data::Dynamic       & io_dynMem ): TimePred< TL_T_REAL,
                                                                 TL_N_RMS,
                                                                 TL_T_EL,
                                                                 TL_O_SP,
                                                                 TL_O_TI,
                                                                 1 >( i_rfs,
                                                                      io_dynMem ) {
      // store stiffness matrices dense
      this->storeStiffTDense( io_dynMem,
                              m_stiffT );

      // generate kernels
      generateKernels();
    };

    /**
     * Applies the Cauchy–Kowalevski procedure (single forward run LIBXSMM version) and computes time derivatives and time integrated DOFs.
     *
     * @param i_dT time step.
     * @param i_starE elastic star matrices.
     * @param i_starA anelastic star matrices, use nullptr if TL_N_RMS==0.
     * @param i_srcA anelastic source matrices, use nullptr if TL_N_RMS==0.
     * @param i_dofsE elastic DOFs.
     * @param i_dofsA anelastic DOFs, use nullptr if TL_N_RMS==0.
     * @param o_scratch will be used as scratch memory.
     * @param o_derE will be set to elastic time derivatives.
     * @param o_derA will be set to anelastic time derivatives (ignored if TL_N_RMS==0, use nullptr).
     * @param o_tIntE will be set to elastic time integrated DOFs.
     * @param o_tIntA will be set to anelastic time integrated DOFS (ignored if TL_N_RMS==0, use nullptr).
     **/
    void ck( TL_T_REAL         i_dT,
             TL_T_REAL const   i_starE[TL_N_DIS][TL_N_ENS_STAR_E],
             TL_T_REAL const (*i_starA)[TL_N_ENS_STAR_A],
             TL_T_REAL const (*i_srcA)[TL_N_ENS_SRC_A],
             TL_T_REAL const   i_dofsE[TL_N_QTS_E][TL_N_MDS][1],
             TL_T_REAL const (*i_dofsA)[TL_N_QTS_M][TL_N_MDS][1],
             TL_T_REAL         o_scratch[TL_N_QTS_E][TL_N_MDS][1],
             TL_T_REAL         o_derE[TL_O_TI][TL_N_QTS_E][TL_N_MDS][1],
             TL_T_REAL       (*o_derA)[TL_O_TI][TL_N_QTS_M][TL_N_MDS][1],
             TL_T_REAL         o_tIntE[TL_N_QTS_E][TL_N_MDS][1],
             TL_T_REAL       (*o_tIntA)[TL_N_QTS_M][TL_N_MDS][1] ) const {
      // relaxation frequencies
      TL_T_REAL const *l_rfs = this->m_rfs;

      // scalar for the time integration
      TL_T_REAL l_scalar = i_dT;

      // initialize zero-derivative, reset time integrated dofs
#ifdef ELTWISE_TPP
      /* 1. o_derE[0][l_qt][l_md][0] = i_dofsE[l_qt][l_md][0] */
      u_unary.execute(0, 0, &i_dofsE[0][0][0], &o_derE[0][0][0][0]);
      /* 2. o_tIntE[l_qt][l_md][0]   = l_scalar * i_dofsE[l_qt][l_md][0] */
      b_binary.execute(0, 0, &l_scalar, &i_dofsE[0][0][0], &o_tIntE[0][0][0]);
#else
      for( unsigned short l_qt = 0; l_qt < TL_N_QTS_E; l_qt++ ) {
#pragma omp simd
        for( unsigned short l_md = 0; l_md < TL_N_MDS; l_md++ ) {
          o_derE[0][l_qt][l_md][0] = i_dofsE[l_qt][l_md][0];
          o_tIntE[l_qt][l_md][0]   = l_scalar * i_dofsE[l_qt][l_md][0];
        }
      }
#endif

      // anelastic: init zero-derivative, reset tDofs
#ifdef ELTWISE_TPP
        /* 1. o_derA[l_rm][0][l_qt][l_md][0] = i_dofsA[l_rm][l_qt][l_md][0]; */
        u_unary.execute(1, 0, &i_dofsA[0][0][0][0], &o_derA[0][0][0][0][0]);
        /* 2. o_tIntA[l_rm][l_qt][l_md][0] = l_scalar * i_dofsA[l_rm][l_qt][l_md][0] */
        b_binary.execute(1, 0, &l_scalar, &i_dofsA[0][0][0][0], &o_tIntA[0][0][0][0]);
#else
      for( unsigned short l_rm = 0; l_rm < TL_N_RMS; l_rm++ ) {
        for( unsigned short l_qt = 0; l_qt < TL_N_QTS_M; l_qt++ ) {
#pragma omp simd
          for( unsigned short l_md = 0; l_md < TL_N_MDS; l_md++ ) {
            o_derA[l_rm][0][l_qt][l_md][0] = i_dofsA[l_rm][l_qt][l_md][0];
            o_tIntA[l_rm][l_qt][l_md][0] = l_scalar * i_dofsA[l_rm][l_qt][l_md][0];
          }
        }
      }
#endif

      // iterate over time derivatives
      for( unsigned short l_de = 1; l_de < TL_O_TI; l_de++ ) {
        // recursive id for the non-zero blocks
        unsigned short l_re = (TL_N_RMS == 0) ? l_de : 1;

        // elastic: reset this derivative
#ifdef ELTWISE_TPP
        u_unary.execute(2, 0, &o_derE[l_de][0][0][0]);
#else
        for( unsigned short l_qt = 0; l_qt < TL_N_QTS_E; l_qt++ )
#pragma omp simd
          for( unsigned short l_md = 0; l_md < TL_N_MDS; l_md++ ) o_derE[l_de][l_qt][l_md][0] = 0;
#endif

        // buffer for the anelastic computations
        TL_T_REAL l_scratch[TL_N_QTS_M][TL_N_MDS];
#ifdef ELTWISE_TPP
        // buffer for relaxation computations
        TL_T_REAL l_scratch2[TL_N_QTS_M][TL_N_MDS];
        // buffer for time integrated dofs
        TL_T_REAL l_scratch3[TL_N_QTS_E][TL_N_MDS];
#endif

        if( TL_N_RMS > 0 ) {
#ifdef ELTWISE_TPP
          u_unary.execute(2, 1, &l_scratch[0][0]);
#else
          for( unsigned short l_qt = 0; l_qt < TL_N_QTS_M; l_qt++ ) {
#pragma omp simd
            for( unsigned short l_md = 0; l_md < TL_N_MDS; l_md++ ) l_scratch[l_qt][l_md] = 0;
          }
#endif
        }

        // compute the derivatives
        for( unsigned short l_di = 0; l_di < TL_N_DIS; l_di++ ) {
          // multiply with transposed stiffness matrices and inverse mass matrix
          m_mm.execute( 0, l_re-1,
                        m_stiffT[l_re-1][l_di],
                        o_derE[l_de-1][0][0],
                        o_scratch[0][0],
                        nullptr,
                        nullptr,
                        nullptr );
          // multiply with star matrices
          m_mm.execute( 1, l_re-1,
                        o_scratch[0][0],
                        i_starE[l_di],
                        o_derE[l_de][0][0],
                        nullptr,
                        nullptr,
                        nullptr );

          if( TL_N_RMS > 0 ) {
            // multiply with anelastic star matrices
            m_mm.execute( 2, 0,
                          o_scratch[TL_N_QTS_M][0],
                          i_starA[l_di],
                          l_scratch[0],
                          nullptr,
                          nullptr,
                          nullptr );
          }
        }

        // update scalar
        l_scalar *= i_dT / (l_de+1);

        for( unsigned short l_rm = 0; l_rm < TL_N_RMS; l_rm++ ) {
          // add contribution of source matrix
          m_mm.execute( 2, 1,
                        o_derA[l_rm][l_de-1][0][0],
                        i_srcA[l_rm],
                        o_derE[l_de][0][0],
                        nullptr,
                        nullptr,
                        nullptr );

          // multiply with relaxation frequency and add
#ifdef ELTWISE_TPP
          /* @TODO: One could use a ternary here (but likely it is not possible right now due
              to the missing support for TERNARY_BCAST flags outside equations) */

          /* 1. l_scratch2[][] = l_scratch[l_qt][l_md] + o_derA[l_rm][l_de-1][l_qt][l_md][0] */
          b_binary.execute(3, 0, &l_scratch[0][0], &o_derA[l_rm][l_de-1][0][0][0], &l_scratch2[0][0]);

          /* 2  o_derA[l_rm][l_de][l_qt][l_md][0] = l_rfs[l_rm] * l_scratch2[l_qt][l_md][0] */
          b_binary.execute(3, 1, &l_rfs[l_rm], &l_scratch2[0][0], &o_derA[l_rm][l_de][0][0][0]);

          /* 3.1 l_scratch2 = l_scalar * o_derA[l_rm][l_de][l_qt][l_md][0] */
          b_binary.execute(3, 2, &l_scalar, &o_derA[l_rm][l_de][0][0][0], &l_scratch2[0][0]);
          /* 3.2 o_tIntA[l_rm][l_qt][l_md][0] += l_scratch2[l_qt][l_md][0] */
          b_binary.execute(3, 3, &o_tIntA[l_rm][0][0][0], &l_scratch2[0][0], &o_tIntA[l_rm][0][0][0]);
#else
          for( unsigned short l_qt = 0; l_qt < TL_N_QTS_M; l_qt++ ) {
#pragma omp simd
            for( unsigned short l_md = 0; l_md < TL_N_MDS; l_md++ ) {
              o_derA[l_rm][l_de][l_qt][l_md][0] = l_rfs[l_rm] * ( l_scratch[l_qt][l_md]+ o_derA[l_rm][l_de-1][l_qt][l_md][0] );
              o_tIntA[l_rm][l_qt][l_md][0] += l_scalar * o_derA[l_rm][l_de][l_qt][l_md][0];
            }
          }
#endif
        }

        // elastic: update time integrated DOFs
#ifdef ELTWISE_TPP
        /* @TODO: One could use a ternary here (but likely it is not possible right now due
            to the missing support for TERNARY_BCAST flags outside equations) */

        /* 1.1 l_scratch3 = l_scalar * o_derE[l_de][l_qt][l_md][0] */
        b_binary.execute(4, l_de-1, &l_scalar, &o_derE[l_de][0][0][0], &l_scratch3[0][0]);

        /* 1.2 o_tIntE[l_qt][l_md][0] += l_scratch3 */
        b_binary.execute(5, l_de-1, &o_tIntE[0][0][0], &l_scratch3[0][0], &o_tIntE[0][0][0]);
#else
        unsigned short l_nCpMds = (TL_N_RMS == 0) ? CE_N_ELEMENT_MODES_CK( TL_T_EL, TL_O_SP, l_de ) : TL_N_MDS;
        for( unsigned short l_qt = 0; l_qt < TL_N_QTS_E; l_qt++ ) {
#pragma omp simd
          for( unsigned short l_md = 0; l_md < l_nCpMds; l_md++ ) {
            o_tIntE[l_qt][l_md][0] += l_scalar * o_derE[l_de][l_qt][l_md][0];
          }
        }
#endif
      }
    }
};

#endif
