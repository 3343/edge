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
#include "data/XsmmUtils.hpp"

#define ELTWISE_TPP

#ifdef ELTWISE_TPP
#  define EQUATION_TPP
#endif

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

#ifdef EQUATION_TPP
    std::vector<libxsmm_matrix_eqn_function>  e_eqns00;
    std::vector<libxsmm_matrix_eqn_function>  e_eqns01;
    std::vector<libxsmm_matrix_eqn_function>  e_eqns1;
#endif

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
      } /* loop over l_de for adding matrix kernels */

#ifdef ELTWISE_TPP
      // initialize zero-derivative, reset time integrated dofs
      u_unary.add(0, TL_N_MDS * TL_N_QTS_E, 1 /* m, n */,  LIBXSMM_MELTW_TYPE_UNARY_IDENTITY, LIBXSMM_MELTW_FLAG_UNARY_NONE);
      b_binary.add(0, TL_N_MDS * TL_N_QTS_E, 1 /* m, n */, LIBXSMM_MELTW_TYPE_BINARY_MUL, LIBXSMM_MELTW_FLAG_BINARY_BCAST_SCALAR_IN_0);

      // anelastic: init zero-derivative, reset tDofs
      u_unary.add(1,  TL_N_MDS * TL_N_QTS_M, TL_N_RMS /* m, n */, TL_N_MDS * TL_N_QTS_M, TL_N_MDS * TL_N_QTS_M * TL_O_TI /*ldi, ldo */, LIBXSMM_MELTW_TYPE_UNARY_IDENTITY, LIBXSMM_MELTW_FLAG_UNARY_NONE);
      b_binary.add(1, TL_N_MDS * TL_N_QTS_M, TL_N_RMS /* m, n */, LIBXSMM_MELTW_TYPE_BINARY_MUL, LIBXSMM_MELTW_FLAG_BINARY_BCAST_SCALAR_IN_0);

      // zeroing: elastic: reset this derivative
      u_unary.add(2, TL_N_MDS * TL_N_QTS_E, 1 /* m, n */, LIBXSMM_MELTW_TYPE_UNARY_XOR, LIBXSMM_MELTW_FLAG_UNARY_NONE);

      // zeroing: buffer for the anelastic computations
      u_unary.add(2, TL_N_MDS * TL_N_QTS_M, 1 /* m, n */, LIBXSMM_MELTW_TYPE_UNARY_XOR, LIBXSMM_MELTW_FLAG_UNARY_NONE);

      // multiply with relaxation frequency and add
      // addition
      b_binary.add(3, TL_N_MDS * TL_N_QTS_M, 1 /* m, n */, LIBXSMM_MELTW_TYPE_BINARY_ADD, LIBXSMM_MELTW_FLAG_BINARY_NONE);
      // mult + assign
      b_binary.add(3, TL_N_MDS * TL_N_QTS_M, 1 /* m, n */, LIBXSMM_MELTW_TYPE_BINARY_MUL, LIBXSMM_MELTW_FLAG_BINARY_BCAST_SCALAR_IN_0);
      // accumulation 2
      b_binary.add(3, TL_N_MDS * TL_N_QTS_M, 1 /* m, n */, LIBXSMM_MELTW_TYPE_BINARY_MUL, LIBXSMM_MELTW_FLAG_BINARY_BCAST_SCALAR_IN_0);
      b_binary.add(3, TL_N_MDS * TL_N_QTS_M, 1 /* m, n */, LIBXSMM_MELTW_TYPE_BINARY_ADD, LIBXSMM_MELTW_FLAG_BINARY_NONE);


      // update time integrated dofs
      for( unsigned short l_de = 1; l_de < TL_O_TI; l_de++ ) {
        unsigned short l_nCpMds = (TL_N_RMS == 0) ? CE_N_ELEMENT_MODES_CK( TL_T_EL, TL_O_SP, l_de ) : TL_N_MDS;
        b_binary.add(4, l_nCpMds, TL_N_QTS_E /* m, n */, TL_N_MDS, TL_N_MDS, TL_N_MDS /* ldi0, ldi1, ldo */,
                      LIBXSMM_MELTW_TYPE_BINARY_MUL, LIBXSMM_MELTW_FLAG_BINARY_BCAST_SCALAR_IN_0);
        b_binary.add(5, l_nCpMds, TL_N_QTS_E /* m, n */, TL_N_MDS, TL_N_MDS, TL_N_MDS /* ldi0, ldi1, ldo */,
                      LIBXSMM_MELTW_TYPE_BINARY_ADD, LIBXSMM_MELTW_FLAG_BINARY_NONE);
      }

#  ifdef EQUATION_TPP
      for( unsigned short l_de = 1; l_de < TL_O_TI; l_de++ ) {
        unsigned short l_nCpMds = (TL_N_RMS == 0) ? CE_N_ELEMENT_MODES_CK( TL_T_EL, TL_O_SP, l_de ) : TL_N_MDS;

        // common part for all equations added below

        libxsmm_datatype dtype      = XsmmDtype<TL_T_REAL>();
        libxsmm_datatype dtype_comp = XsmmDtype<TL_T_REAL>();

        libxsmm_meqn_arg_shape  eqn_out_arg_shape;
        libxsmm_meqn_arg_shape  arg_shape;

        libxsmm_matrix_arg_attributes arg_singular_attr;

        libxsmm_matrix_eqn_arg_metadata arg_metadata;
        libxsmm_matrix_eqn_op_metadata  op_metadata;

        libxsmm_bitfield binary_flags;
        libxsmm_bitfield ternary_flags;

        arg_singular_attr.type = LIBXSMM_MATRIX_ARG_TYPE_SINGULAR;

        // adding to eqns00 (multiply with relaxation frequency and add, part 1)

        libxsmm_blasint my_eqn00 = libxsmm_matrix_eqn_create();         /* o_derA[l_de] = l_rfs[l_rm] * (l_scratch + o_derA[l-de-1]) */

        binary_flags             = LIBXSMM_MELTW_FLAG_BINARY_BCAST_SCALAR_IN_1;
        op_metadata.eqn_idx      = my_eqn00;
        op_metadata.op_arg_pos   = -1;
        libxsmm_matrix_eqn_push_back_binary_op_v2(op_metadata, LIBXSMM_MELTW_TYPE_BINARY_MUL, dtype_comp, binary_flags);

        binary_flags             = LIBXSMM_MELTW_FLAG_BINARY_NONE;
        op_metadata.eqn_idx      = my_eqn00;
        op_metadata.op_arg_pos   = -1;
        libxsmm_matrix_eqn_push_back_binary_op_v2(op_metadata, LIBXSMM_MELTW_TYPE_BINARY_ADD, dtype_comp, binary_flags);

        arg_metadata.eqn_idx     = my_eqn00;
        arg_metadata.in_arg_pos  = 0;
        arg_shape.m    = TL_N_MDS;                                      /* l_scratch, [TL_N_MDS][TL_N_QTS_M] */
        arg_shape.n    = TL_N_QTS_M;
        arg_shape.ld   = TL_N_MDS;
        arg_shape.type = dtype;
        libxsmm_matrix_eqn_push_back_arg_v2(arg_metadata, arg_shape, arg_singular_attr);

        arg_metadata.eqn_idx     = my_eqn00;
        arg_metadata.in_arg_pos  = 1;
        arg_shape.m    = TL_N_MDS;                                      /* o_derA[l_de-1], [TL_N_MDS][TL_N_QTS_M] */
        arg_shape.n    = TL_N_QTS_M;
        arg_shape.ld   = TL_N_MDS;
        arg_shape.type = dtype;
        libxsmm_matrix_eqn_push_back_arg_v2(arg_metadata, arg_shape, arg_singular_attr);

        arg_metadata.eqn_idx     = my_eqn00;
        arg_metadata.in_arg_pos  = 2;
        arg_shape.m    = 1;                                             /* l_rfs[l_rm], [1] */
        arg_shape.n    = 1;
        arg_shape.ld   = 1;
        arg_shape.type = dtype;
        libxsmm_matrix_eqn_push_back_arg_v2(arg_metadata, arg_shape, arg_singular_attr);

        eqn_out_arg_shape.m    = TL_N_MDS;                             /* o_derA[l_de], [TL_N_MDS][TL_N_QTS_M] */
        eqn_out_arg_shape.n    = TL_N_QTS_M;
        eqn_out_arg_shape.ld   = TL_N_MDS;
        eqn_out_arg_shape.type = dtype;

        libxsmm_matrix_eqn_tree_print( my_eqn00 );
        libxsmm_matrix_eqn_rpn_print ( my_eqn00 );
        libxsmm_matrix_eqn_function func00 = libxsmm_dispatch_matrix_eqn_v2( my_eqn00, eqn_out_arg_shape );
        if ( func00 == NULL) {
          fprintf( stderr, "JIT for TPP equation func00 (eqn00) failed. Bailing...!\n");
          exit(-1);
        }
        e_eqns00.push_back(func00);

        // adding to eqns01 (multiply with relaxation frequency and add, part 2)

        libxsmm_blasint my_eqn01 = libxsmm_matrix_eqn_create();         /* o_tintA += l_scalar * o_derA */

        ternary_flags            = LIBXSMM_MELTW_FLAG_TERNARY_BCAST_SCALAR_IN_0;
        op_metadata.eqn_idx      = my_eqn01;
        op_metadata.op_arg_pos   = -1;
        libxsmm_matrix_eqn_push_back_ternary_op_v2(op_metadata, LIBXSMM_MELTW_TYPE_TERNARY_MULADD, dtype_comp, ternary_flags);

        arg_metadata.eqn_idx     = my_eqn01;
        arg_metadata.in_arg_pos  = 0;
        arg_shape.m    = 1;                                             /* l_scalar, [1] */
        arg_shape.n    = 1;
        arg_shape.ld   = 1;
        arg_shape.type = dtype;
        libxsmm_matrix_eqn_push_back_arg_v2(arg_metadata, arg_shape, arg_singular_attr);

        arg_metadata.eqn_idx     = my_eqn01;
        arg_metadata.in_arg_pos  = 1;
        arg_shape.m    = TL_N_MDS;                                      /* o_derA, [TL_N_MDS][TL_N_QTS_M] */
        arg_shape.n    = TL_N_QTS_M;
        arg_shape.ld   = TL_N_MDS;
        arg_shape.type = dtype;
        libxsmm_matrix_eqn_push_back_arg_v2(arg_metadata, arg_shape, arg_singular_attr);

        arg_metadata.eqn_idx     = my_eqn01;
        arg_metadata.in_arg_pos  = 2;
        arg_shape.m    = TL_N_MDS;                                      /* o_tintA, [TL_N_MDS][TL_N_QTS_M] */
        arg_shape.n    = TL_N_QTS_M;
        arg_shape.ld   = TL_N_MDS;
        arg_shape.type = dtype;
        libxsmm_matrix_eqn_push_back_arg_v2(arg_metadata, arg_shape, arg_singular_attr);

        eqn_out_arg_shape.m    = TL_N_MDS;                             /* o_tintA, [TL_N_MDS][TL_N_QTS_M] */
        eqn_out_arg_shape.n    = TL_N_QTS_M;
        eqn_out_arg_shape.ld   = TL_N_MDS;
        eqn_out_arg_shape.type = dtype;

        libxsmm_matrix_eqn_tree_print( my_eqn01 );
        libxsmm_matrix_eqn_rpn_print ( my_eqn01 );
        libxsmm_matrix_eqn_function func01 = libxsmm_dispatch_matrix_eqn_v2( my_eqn01, eqn_out_arg_shape );
        if ( func01 == NULL) {
          fprintf( stderr, "JIT for TPP equation func01 (eqn01) failed. Bailing...!\n");
          exit(-1);
        }
        e_eqns01.push_back(func01);

        // adding to eqns1 (update time integrated DOFs)

        libxsmm_blasint my_eqn1 = libxsmm_matrix_eqn_create();          /* o_tintE += l_scalar * o_derE */

        ternary_flags            = LIBXSMM_MELTW_FLAG_TERNARY_BCAST_SCALAR_IN_0;
        op_metadata.eqn_idx      = my_eqn1;
        op_metadata.op_arg_pos   = -1;
        libxsmm_matrix_eqn_push_back_ternary_op_v2(op_metadata, LIBXSMM_MELTW_TYPE_TERNARY_MULADD, dtype_comp, ternary_flags);

        arg_metadata.eqn_idx     = my_eqn1;
        arg_metadata.in_arg_pos  = 0;
        arg_shape.m    = 1;                                             /* l_scalar, [1]*/
        arg_shape.n    = 1;
        arg_shape.ld   = 1;
        arg_shape.type = dtype;
        libxsmm_matrix_eqn_push_back_arg_v2(arg_metadata, arg_shape, arg_singular_attr);

        arg_metadata.eqn_idx     = my_eqn1;
        arg_metadata.in_arg_pos  = 1;
        arg_shape.m    = l_nCpMds;                                      /* o_derE[l_de], [l_ncpMds*][TL_N_QTS_E] */
        arg_shape.n    = TL_N_QTS_E;
        arg_shape.ld   = TL_N_MDS;
        arg_shape.type = dtype;
        libxsmm_matrix_eqn_push_back_arg_v2(arg_metadata, arg_shape, arg_singular_attr);

        arg_metadata.eqn_idx     = my_eqn1;
        arg_metadata.in_arg_pos  = 2;
        arg_shape.m    = l_nCpMds;                                      /* o_tIntE, [l_ncpMds*][TL_N_QTS_E] */
        arg_shape.n    = TL_N_QTS_E;
        arg_shape.ld   = TL_N_MDS;
        arg_shape.type = dtype;
        libxsmm_matrix_eqn_push_back_arg_v2(arg_metadata, arg_shape, arg_singular_attr);


        eqn_out_arg_shape.m    = l_nCpMds;                             /* o_tIntE, [l_ncpMds*][TL_N_QTS_E] */
        eqn_out_arg_shape.n    = TL_N_QTS_E;
        eqn_out_arg_shape.ld   = TL_N_MDS;
        eqn_out_arg_shape.type = dtype;

        libxsmm_matrix_eqn_tree_print( my_eqn1 );
        libxsmm_matrix_eqn_rpn_print ( my_eqn1 );
        libxsmm_matrix_eqn_function func1 = libxsmm_dispatch_matrix_eqn_v2( my_eqn1, eqn_out_arg_shape );
        if ( func1 == NULL) {
          fprintf( stderr, "JIT for TPP equation func1 (eqn1) failed. Bailing...!\n");
          exit(-1);
        }
        e_eqns1.push_back(func1);
      } /* loop over l_de for TPP equations */
#  endif
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

#if defined(ELTWISE_TPP) and defined(EQUATION_TPP)
      libxsmm_matrix_arg arg_array[3];
      libxsmm_matrix_eqn_param eqn_param;
      memset( &eqn_param, 0, sizeof(eqn_param));
      eqn_param.inputs = arg_array;
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
#if defined(ELTWISE_TPP) and !defined(EQUATION_TPP)
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
#  ifdef EQUATION_TPP
          arg_array[0].primary     = &l_scratch[0][0];               
          arg_array[1].primary     = &o_derA[l_rm][l_de-1][0][0][0];
          arg_array[2].primary     = const_cast<void*>(reinterpret_cast<const void*>(&l_rfs[l_rm]));
          eqn_param.output.primary = &o_derA[l_rm][l_de][0][0][0];
          e_eqns00[l_de-1](&eqn_param);

          arg_array[0].primary     = &l_scalar;
          arg_array[1].primary     = &o_derA[l_rm][l_de][0][0][0];
          arg_array[2].primary     = &o_tIntA[l_rm][0][0][0];
          eqn_param.output.primary = &o_tIntA[l_rm][0][0][0];
          e_eqns01[l_de-1](&eqn_param);
#  else
          /* 1. l_scratch2[][] = l_scratch[l_qt][l_md] + o_derA[l_rm][l_de-1][l_qt][l_md][0] */
          b_binary.execute(3, 0, &l_scratch[0][0], &o_derA[l_rm][l_de-1][0][0][0], &l_scratch2[0][0]);

          /* 2  o_derA[l_rm][l_de][l_qt][l_md][0] = l_rfs[l_rm] * l_scratch2[l_qt][l_md][0] */
          b_binary.execute(3, 1, &l_rfs[l_rm], &l_scratch2[0][0], &o_derA[l_rm][l_de][0][0][0]);

          /* 3.1 l_scratch2 = l_scalar * o_derA[l_rm][l_de][l_qt][l_md][0] */
          b_binary.execute(3, 2, &l_scalar, &o_derA[l_rm][l_de][0][0][0], &l_scratch2[0][0]);
          /* 3.2 o_tIntA[l_rm][l_qt][l_md][0] += l_scratch2[l_qt][l_md][0] */
          b_binary.execute(3, 3, &o_tIntA[l_rm][0][0][0], &l_scratch2[0][0], &o_tIntA[l_rm][0][0][0]);
#  endif

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
#  ifdef EQUATION_TPP
        arg_array[0].primary     = &l_scalar;                                      /* [1] */
        arg_array[1].primary     = &o_derE[l_de][0][0][0];                         /* [l_nCpMds, TL_N_QTS_E] */
        arg_array[2].primary     = &o_tIntE[0][0][0];                              /* [l_nCpMds, TL_N_QTS_E] */
        eqn_param.output.primary = &o_tIntE[0][0][0];                              /* [l_nCpMds, TL_N_QTS_E] */
        e_eqns1[l_de-1](&eqn_param);
#  else
        /* @TODO: One could use a ternary here (but likely it is not possible right now due
            to the missing support for TERNARY_BCAST flags outside equations) */

        /* 1.1 l_scratch3 = l_scalar * o_derE[l_de][l_qt][l_md][0] */
        b_binary.execute(4, l_de-1, &l_scalar, &o_derE[l_de][0][0][0], &l_scratch3[0][0]);

        /* 1.2 o_tIntE[l_qt][l_md][0] += l_scratch3 */
        b_binary.execute(5, l_de-1, &o_tIntE[0][0][0], &l_scratch3[0][0], &o_tIntE[0][0][0]);
#  endif
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
