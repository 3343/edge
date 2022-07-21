/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (alex.breuer AT uni-jena.de)
 *         Alexander Heinecke (alexander.heinecke AT intel.com)
 *
 * @section LICENSE
 * Copyright (c) 2021, Friedrich Schiller University Jena
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
 * Optimized quadrature-free ADER-DG volume integration for single seismic wave propagation.
 **/
#ifndef EDGE_SEISMIC_KERNELS_VOL_INT_SINGLE_HPP
#define EDGE_SEISMIC_KERNELS_VOL_INT_SINGLE_HPP

#include "VolInt.hpp"
#include "dg/Basis.h"
#include "data/MmXsmmSingle.hpp"
#include "data/UnaryXsmm.hpp"
#include "data/BinaryXsmm.hpp"
#include "data/TernaryXsmm.hpp"

#define ELTWISE_TPP
//#define USE_TERNARY

#ifdef USE_TERNARY
  #error "USE_TERNARY requires (currently missing) support for TERNARY_BCAST flags in LIBXSMM (hence switched off here)"
#endif

namespace edge {
  namespace seismic {
    namespace kernels { 
      template< typename       TL_T_REAL,
                unsigned short TL_N_RMS,
                t_entityType   TL_T_EL,
                unsigned short TL_O_SP >
      class VolIntSingle;
    }
  }
}

/**
 * Optimized quadrature-free ADER-DG volume integration for single seismic forward simulations.
 *
 * @paramt TL_T_REAL floating point precision.
 * @paramt TL_N_RMS number of relaxation mechanisms.
 * @paramt TL_T_EL element type.
 * @paramt TL_O_SP spatial order.
 **/
template< typename       TL_T_REAL,
          unsigned short TL_N_RMS,
          t_entityType   TL_T_EL,
          unsigned short TL_O_SP >
class edge::seismic::kernels::VolIntSingle: edge::seismic::kernels::VolInt < TL_T_REAL,
                                                                             TL_N_RMS,
                                                                             TL_T_EL,
                                                                             TL_O_SP,
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

    //! matrix kernels
    edge::data::MmXsmmSingle< TL_T_REAL > m_mm;

    //! unary kernels
    edge::data::UnaryXsmm< TL_T_REAL > u_unary;

    //! binary kernels
    edge::data::BinaryXsmm< TL_T_REAL > b_binary;

    //! ternary kernels
    edge::data::TernaryXsmm< TL_T_REAL > t_ternary;

    //! pointers to the stiffness matrices
    TL_T_REAL *m_stiff[TL_N_DIS] = {};

    /**
     * Generates the matrix kernels for the stiffness matrices and star matrices.
     **/
    void generateKernels() {
      // stiffness matrix
      m_mm.add( 0,                                            // group
                TL_N_MDS,                                     // m
                TL_N_QTS_E,                                   // n
                CE_N_ELEMENT_MODES_CK( TL_T_EL, TL_O_SP, 1 ), // k
                TL_N_MDS,                                     // ldA
                TL_N_MDS,                                     // ldB
                TL_N_MDS,                                     // ldC
                static_cast<TL_T_REAL>(1.0),                  // alpha
                static_cast<TL_T_REAL>(0.0),                  // beta
                LIBXSMM_GEMM_PREFETCH_NONE );

      // elastic star matrix
      m_mm.add( 0,                           // group
                TL_N_MDS,                    // m
                TL_N_QTS_E,                  // n
                TL_N_QTS_E,                  // k
                TL_N_MDS,                    // ldA
                TL_N_QTS_E,                  // ldB
                TL_N_MDS,                    // ldC
                static_cast<TL_T_REAL>(1.0), // alpha
                static_cast<TL_T_REAL>(1.0), // beta
                LIBXSMM_GEMM_PREFETCH_NONE );

      if( TL_N_RMS > 0 ) {
        // anelastic star matrix
        m_mm.add( 1,                           // group
                  TL_N_MDS,                    // m
                  TL_N_QTS_M,                  // n
                  TL_N_DIS,                    // k
                  TL_N_MDS,                    // ldA
                  TL_N_DIS,                    // ldB
                  TL_N_MDS,                    // ldC
                  static_cast<TL_T_REAL>(1.0), // alpha
                  static_cast<TL_T_REAL>(1.0), // beta
                  LIBXSMM_GEMM_PREFETCH_NONE );

        // anelastic source matrix
        m_mm.add( 1,                           // group
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

#ifdef ELTWISE_TPP
      // zeroing the scratch
      u_unary.add(0, TL_N_MDS, TL_N_QTS_M /* m, n */, TL_N_MDS, TL_N_MDS /* ldi, ldo */, LIBXSMM_MELTW_TYPE_UNARY_XOR, LIBXSMM_MELTW_FLAG_UNARY_NONE);

      // subtraction
      b_binary.add(0, TL_N_MDS, TL_N_QTS_M /* m, n */, TL_N_MDS, TL_N_MDS, TL_N_MDS /* ldi0, ldi1, ldo */,
                    LIBXSMM_MELTW_TYPE_BINARY_SUB, LIBXSMM_MELTW_FLAG_BINARY_NONE);
      // accumulation
#ifdef USE_TERNARY
      /* TERNARY_BCAST is only supported in equations but not in standalone TPPs */
      t_ternary.add(0, TL_N_MDS, TL_N_QTS_M /* m, n */, TL_N_MDS, TL_N_MDS, TL_N_MDS, TL_N_MDS /* ldi0, ldi1, ldi2, ldo */,
                    LIBXSMM_MELTW_TYPE_TERNARY_MULADD, LIBXSMM_MELTW_FLAG_TERNARY_BCAST_SCALAR_IN_1);
#else
      b_binary.add(1, TL_N_MDS, TL_N_QTS_M /* m, n */, TL_N_MDS, TL_N_MDS, TL_N_MDS /* ldi0, ldi1, ldo */,
                    LIBXSMM_MELTW_TYPE_BINARY_MUL, LIBXSMM_MELTW_FLAG_BINARY_BCAST_SCALAR_IN_0);
      b_binary.add(1, TL_N_MDS, TL_N_QTS_M /* m, n */, TL_N_MDS, TL_N_MDS, TL_N_MDS /* ldi0, ldi1, ldo */,
                    LIBXSMM_MELTW_TYPE_BINARY_ADD, LIBXSMM_MELTW_FLAG_BINARY_NONE);
#endif

#endif
    }

  public:
    /**
     * Constructor of the optimized volume integration for single forward simulations.
     *
     * @param i_rfs relaxation frequencies, use nullptr if TL_N_RMS==0.
     * @param io_dynMem dynamic memory allocations.
     **/
    VolIntSingle( TL_T_REAL     const * i_rfs,
                  data::Dynamic       & io_dynMem ): VolInt< TL_T_REAL,
                                                             TL_N_RMS,
                                                             TL_T_EL,
                                                             TL_O_SP,
                                                             1 >( i_rfs,
                                                                  io_dynMem ) {
      // store stiffness matrices dense
      this->storeStiffDense( io_dynMem,
                             m_stiff );

      // generate matrix kernels
      generateKernels();
    }

    /**
     * Optimized volume contribution for single forward simulations.
     *
     * @param i_starE elastic star matrices.
     * @param i_starA anelastic star matrices, use nullptr if TL_N_RMS==0.
     * @param i_tDofsE time integrated elastic DOFs.
     * @param i_tDofsA time integrated anselastic DOFs.
     * @param io_dofsE will be updated with local elastic contribution of the element to the volume integral.
     * @param io_dofsA will be updated with local anelastic contribution of the element to the volume integral, use nullptr if TL_N_RMS==0.
     * @param o_scratch will be used as scratch space for the computations.
     **/
    void apply( TL_T_REAL const   i_starE[TL_N_DIS][TL_N_ENS_STAR_E],
                TL_T_REAL const (*i_starA)[TL_N_ENS_STAR_A],
                TL_T_REAL const (*i_srcA)[TL_N_ENS_SRC_A],
                TL_T_REAL const   i_tDofsE[TL_N_QTS_E][TL_N_MDS][1],
                TL_T_REAL const (*i_tDofsA)[TL_N_QTS_M][TL_N_MDS][1],
                TL_T_REAL         io_dofsE[TL_N_QTS_E][TL_N_MDS][1],
                TL_T_REAL       (*io_dofsA)[TL_N_QTS_M][TL_N_MDS][1],
                TL_T_REAL         o_scratch[TL_N_QTS_E][TL_N_MDS][1] ) const {
      // relaxation frequencies
      TL_T_REAL const *l_rfs = this->m_rfs;

      // buffer for anelastic part
      TL_T_REAL l_scratch[TL_N_QTS_M][TL_N_MDS][1];
      // buffer for relaxation computations
      TL_T_REAL l_scratch2[TL_N_QTS_M][TL_N_MDS][1];

/*
      TL_T_REAL l_dbg1[TL_N_QTS_M][TL_N_MDS][1];
      TL_T_REAL l_dbg2[TL_N_QTS_M][TL_N_MDS][1];
      TL_T_REAL l_dbg3[TL_N_QTS_M][TL_N_MDS][1];
*/

#ifdef ELTWISE_TPP
      u_unary.execute (0, 0, (TL_T_REAL*)&l_scratch[0][0][0]);
#else
      for( unsigned short l_qt = 0; l_qt < TL_N_QTS_M; l_qt++ )
#pragma omp simd
        for( unsigned short l_md = 0; l_md < TL_N_MDS; l_md++ ) l_scratch[l_qt][l_md][0] = 0;
#endif

      // iterate over dimensions
      for( unsigned short l_di = 0; l_di < TL_N_DIS; l_di++ ) {
        // stiffness and inverse mass matrix
        m_mm.execute( 0, 0,
                      m_stiff[l_di],
                      i_tDofsE[0][0],
                      o_scratch[0][0],
                      nullptr,
                      nullptr,
                      nullptr );

        // star matrix
        m_mm.execute( 0, 1, 
                      o_scratch[0][0],
                      i_starE[l_di],
                      io_dofsE[0][0],
                      nullptr,
                      nullptr,
                      nullptr );

        if( TL_N_RMS > 0 ) {
          // anelastic star matrix
          m_mm.execute( 1, 0,
                        o_scratch[TL_N_QTS_M][0],
                        i_starA[l_di],
                        l_scratch[0][0],
                        nullptr,
                        nullptr,
                        nullptr );
        }
      }

      for( unsigned short l_rm = 0; l_rm < TL_N_RMS; l_rm++ ) {
        // add contribution of source matrix
        m_mm.execute( 1, 1,
                      i_tDofsA[l_rm][0][0],
                      i_srcA[l_rm],
                      io_dofsE[0][0],
                      nullptr,
                      nullptr,
                      nullptr );

        // multiply with relaxation frequency and add
//#if 0
#ifdef ELTWISE_TPP
/*
        for( unsigned short l_qt = 0; l_qt < TL_N_QTS_M; l_qt++ )
          for( unsigned short l_md = 0; l_md < TL_N_MDS; l_md++ )
            l_dbg1[l_qt][l_md][0] = io_dofsA[l_rm][l_qt][l_md][0];

        for( unsigned short l_qt = 0; l_qt < TL_N_QTS_M; l_qt++ )
          for( unsigned short l_md = 0; l_md < TL_N_MDS; l_md++ )
            l_dbg2[l_qt][l_md][0] = ( l_scratch[l_qt][l_md][0] - i_tDofsA[l_rm][l_qt][l_md][0] );

        for( unsigned short l_qt = 0; l_qt < TL_N_QTS_M; l_qt++ )
          for( unsigned short l_md = 0; l_md < TL_N_MDS; l_md++ )
            l_dbg3[l_qt][l_md][0] = l_rfs[l_rm] * ( l_scratch[l_qt][l_md][0] - i_tDofsA[l_rm][l_qt][l_md][0] );
*/

        // @TODO: Could be replaced by a single equation (provided no data races occur)
        /* 1. l_scratch2[][] = l_scratch[l_qt][l_md][0] - i_tDofsA[l_rm][l_qt][l_md][0] */
        b_binary.execute(0, 0, &l_scratch[0][0][0], &i_tDofsA[l_rm][0][0][0], &l_scratch2[0][0][0]);

/*
        for( unsigned short l_qt = 0; l_qt < TL_N_QTS_M; l_qt++ )
          for( unsigned short l_md = 0; l_md < TL_N_MDS; l_md++ )
            if (fabs( (float)(l_scratch2[l_qt][l_md][0] - l_dbg2[l_qt][l_md][0])) > 1.0e-5)
              printf("l_scratch2 %p is different from dbg2 [%d][%d] = %15.15f %15.15f \n", (void*)(&l_scratch2[0][0][0]), l_qt, l_md, (float)l_scratch2[l_qt][l_md][0], (float)l_dbg2[l_qt][l_md][0]);
*/

#ifdef USE_TERNARY
        /* This does not work, see comments in generateKernels() */
        /* 2. io_dofsA[l_rm][l_qt][l_md][0] += l_rfs[l_rm] * l_scratch2[l_qt][l_md][0] */
        t_ternary.execute(0, 0, &io_dofsA[l_rm][0][0][0], &l_rfs[l_rm], &l_scratch2[0][0][0], &io_dofsA[l_rm][0][0][0]);
#else
        /* 2.1 l_scratch2[l_qt][l_md][0] *= l_rfs[l_rm] */
        b_binary.execute(1, 0, &l_rfs[l_rm], &l_scratch2[0][0][0], &l_scratch2[0][0][0]);
/*
        for( unsigned short l_qt = 0; l_qt < TL_N_QTS_M; l_qt++ )
          for( unsigned short l_md = 0; l_md < TL_N_MDS; l_md++ )
            if (fabs( (float)(l_scratch2[l_qt][l_md][0] - l_rfs[l_rm] * l_dbg2[l_qt][l_md][0])) > 1.0e-5)
              printf("l_scratch2 %p is different from dbg2 * l_rfs [%d][%d] = %15.15f %15.15f \n", (void*)(&l_scratch2[0][0][0]), l_qt, l_md, (float)l_scratch2[l_qt][l_md][0], (float)l_rfs[l_rm] * (float)l_dbg2[l_qt][l_md][0]);
*/
        /* 2.1 io_dofsA[l_rm][l_qt][l_md][0] += l_scratch2[l_qt][l_md][0] */
        b_binary.execute(1, 1, &io_dofsA[l_rm][0][0][0], &l_scratch2[0][0][0], &io_dofsA[l_rm][0][0][0]);
#endif
/*
        for( unsigned short l_qt = 0; l_qt < TL_N_QTS_M; l_qt++ )
          for( unsigned short l_md = 0; l_md < TL_N_MDS; l_md++ )
            if (fabs( (float)(io_dofsA[l_rm][l_qt][l_md][0] - (l_dbg1[l_qt][l_md][0] + l_dbg3[l_qt][l_md][0]))) > 1.0e-5)
              printf("final updated value %p is different from dbg1 + dbg3 [%d][%d] = %15.15f %15.15f \n", (void*)(&io_dofsA[0][0][0][0]), l_qt, l_md,
                        (float)io_dofsA[l_rm][l_qt][l_md][0], (float)l_dbg1[l_qt][l_md][0] + (float)l_dbg3[l_qt][l_md][0]);
*/
#else
        for( unsigned short l_qt = 0; l_qt < TL_N_QTS_M; l_qt++ ) {
#pragma omp simd
          for( unsigned short l_md = 0; l_md < TL_N_MDS; l_md++ ) {
            io_dofsA[l_rm][l_qt][l_md][0] += l_rfs[l_rm] * ( l_scratch[l_qt][l_md][0] - i_tDofsA[l_rm][l_qt][l_md][0] );
          }
        }
#endif
      }
    }
};

#endif
