/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *         Alexander Heinecke (alexander.heinecke AT intel.com)
 *
 * @section LICENSE
 * Copyright (c) 2017, Regents of the University of California
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
 * Generic initialization of the matrix kernels.
 **/

#ifndef EDGE_SEISMIC_KERNELS_MM_HPP
#define EDGE_SEISMIC_KERNELS_MM_HPP

#include "constants.hpp"
#include "io/logging.h"

namespace edge {
  namespace elastic {
    namespace setups {
      class MmKernels;
    }
  }
}

/**
 * Sets up matrix kernels used for seismic simulations.
 **/
class edge::elastic::setups::MmKernels {
  public:
#if defined PP_T_KERNELS_XSMM_DENSE_SINGLE
    /**
    * Adds single (non-fused) LIBXSMM matrix multiplication kernels for the ADER-DG solver.
    * Remark: LIBXSMM is column-major, EDGE is row-major.
    *
    * @param i_tEl entity type.
    * @param i_order order of convergence (space and time).
    * @param i_nQts number of quantities.
    * @param i_nCrs number fused runs.
    * @param io_mm matrix-matrix multiplication kernels to to which the single, non-fused ADER-DG LIBXSMM kernels will be added.
    *
    * @paramt TL_T_REAL floating point precision.
    **/
    template< typename TL_T_REAL >
    static void add( t_entityType                     i_tEl,
                     unsigned short                   i_order,
                     unsigned short                   i_nQts,
                     unsigned short                   i_nCrs,                     
                     data::MmXsmmSingle< TL_T_REAL > &io_mm ) {
      // check for non-fused setup
      EDGE_CHECK( i_nCrs == 1 );

      unsigned short l_nMdsFa = CE_N_ELEMENT_MODES( C_ENT[i_tEl].TYPE_FACES, i_order );
      unsigned short l_nMdsEl = CE_N_ELEMENT_MODES( i_tEl,                   i_order );

      // (O-1)*2 kernels for time integration
      // (multiplication with transposed stiffness matrix and star matrix)
      for( unsigned int l_de = 1; l_de < i_order; l_de++ ) {
        // multiplication with transpose stiffness matrix
        io_mm.add( CE_N_ELEMENT_MODES_CK( i_tEl, i_order, l_de ),    // m
                   i_nQts,                                           // n
                   CE_N_ELEMENT_MODES_CK( i_tEl, i_order, l_de-1 ),  // k
                   l_nMdsEl,                                         // ldA
                   l_nMdsEl,                                         // ldB
                   l_nMdsEl,                                         // ldC
                   static_cast<real_base>(1.0),                      // alpha
                   static_cast<real_base>(0.0),                      // beta
                   LIBXSMM_PREFETCH_NONE );
      
        // multiplication with star matrix
        io_mm.add( CE_N_ELEMENT_MODES_CK( i_tEl, i_order, l_de ), // m
                   i_nQts,                                        // n
                   i_nQts,                                        // k
                   l_nMdsEl,                                      // ldA
                   i_nQts,                                        // ldB
                   l_nMdsEl,                                      // ldC
                   static_cast<real_base>(1.0),                   // alpha
                   static_cast<real_base>(1.0),                   // beta
                   LIBXSMM_PREFETCH_NONE );
      }

      // add two volume integration kernels
      // (multiplication with stiffness matrix and star matrix)
      io_mm.add( l_nMdsEl,                                   // m
                 i_nQts,                                     // n
                 CE_N_ELEMENT_MODES_CK( i_tEl, i_order, 1 ), // k
                 l_nMdsEl,                                   // ldA
                 l_nMdsEl,                                   // ldB
                 l_nMdsEl,                                   // ldC
                 static_cast<real_base>(1.0),                // alpha
                 static_cast<real_base>(0.0),                // beta
                 LIBXSMM_PREFETCH_NONE ); // (ORDER-1)*2

      io_mm.add( l_nMdsEl,                    // m
                 i_nQts,                      // n
                 i_nQts,                      // k
                 l_nMdsEl,                    // ldA
                 i_nQts,                      // ldB
                 l_nMdsEl,                    // ldC
                 static_cast<real_base>(1.0), // alpha
                 static_cast<real_base>(1.0), // beta
                 LIBXSMM_PREFETCH_NONE ); // (ORDER-1)*2+1

      // add first flux matrix
      io_mm.add( l_nMdsFa,                    // m
                 i_nQts,                      // n
                 l_nMdsEl,                    // k
                 l_nMdsFa,                    // ldA
                 l_nMdsEl,                    // ldB
                 l_nMdsFa,                    // ldC
                 static_cast<real_base>(1.0), // alpha
                 static_cast<real_base>(0.0), // beta
                 LIBXSMM_PREFETCH_NONE ); // (ORDER-1)*2+2

      // add flux solver
      io_mm.add( l_nMdsFa,                    // m
                 i_nQts,                      // n
                 i_nQts,                      // k
                 l_nMdsFa,                    // ldA
                 i_nQts,                      // ldB
                 l_nMdsFa,                    // ldC
                 static_cast<real_base>(1.0), // alpha
                 static_cast<real_base>(0.0), // beta
                 LIBXSMM_PREFETCH_NONE ); // (ORDER-1)*2+3

      // add second flux matrix
      io_mm.add( l_nMdsEl,                    // m
                 i_nQts,                      // n
                 l_nMdsFa,                    // k
                 l_nMdsEl,                    // ldA
                 l_nMdsFa,                    // ldB
                 l_nMdsEl,                    // ldC
                 static_cast<real_base>(1.0), // alpha
                 static_cast<real_base>(1.0), // beta
                 LIBXSMM_PREFETCH_NONE ); // (ORDER-1)*2+4
    }
#endif

#if defined PP_T_KERNELS_VANILLA
    /**
    * Adds single and fused matrix multiplication kernels for the ADER-DG solver.
    *
    * @param i_tEl entity type.
    * @param i_order order of convergence (space and time).
    * @param i_nQts number of quantities.
    * @param i_nCrs number of fused simulations.
    * @param io_kernels matrix-matrix multiplication kernels to which to which the ADER-DG vanilla kernels will be added.
    *
    * @paramt TL_T_REAL floating point precision.
    **/
    template< typename TL_T_REAL >
    static void add( t_entityType                  i_tEl,
                     unsigned short                i_order,
                     unsigned short                i_nQts,
                     unsigned short                i_nCrs,
                     data::MmVanilla< TL_T_REAL > &io_mm ) {
      unsigned short l_nMdsFa = CE_N_ELEMENT_MODES( C_ENT[i_tEl].TYPE_FACES, i_order );
      unsigned short l_nMdsEl = CE_N_ELEMENT_MODES( i_tEl, i_order );

      // (O-1)*2 kernels for time integration
      // (multiplication with transposed stiffness matrix and star matrix)
      for( unsigned int l_de = 1; l_de < i_order; l_de++ ) {
        // multiplication with transpose stiffness matrix
        io_mm.add( i_nQts,                                           // m
                   CE_N_ELEMENT_MODES_CK( i_tEl, i_order, l_de-1 ),  // n
                   CE_N_ELEMENT_MODES_CK( i_tEl, i_order, l_de ),    // k
                   l_nMdsEl,                                         // ldA
                   l_nMdsEl,                                         // ldB
                   l_nMdsEl,                                         // ldC
                   static_cast<real_base>(1.0),                      // alpha
                   static_cast<real_base>(0.0),                      // beta
                   true,                                             // fused AC
                   false,                                            // fused BC
                   i_nCrs );
      
        // // multiplication with star matrix
        io_mm.add( i_nQts,                                        // m
                   CE_N_ELEMENT_MODES_CK( i_tEl, i_order, l_de ), // n
                   i_nQts,                                        // k
                   i_nQts,                                        // ldA
                   l_nMdsEl,                                      // ldB
                   l_nMdsEl,                                      // ldC
                   static_cast<real_base>(1.0),                   // alpha
                   static_cast<real_base>(1.0),                   // beta
                   false,                                         // fused AC
                   true,                                          // fused BC
                   i_nCrs );
      }

      // add two volume integration kernels
      // (multiplication with stiffness matrix and star matrix)
      io_mm.add( i_nQts,                                     // m
                 l_nMdsEl,                                   // n
                 CE_N_ELEMENT_MODES_CK( i_tEl, i_order, 1 ), // k
                 l_nMdsEl,                                   // ldA
                 l_nMdsEl,                                   // ldB
                 l_nMdsEl,                                   // ldC
                 static_cast<real_base>(1.0),                // alpha
                 static_cast<real_base>(0.0),                // beta
                 true,                                       // fused AC
                 false,                                      // fused BC
                 i_nCrs ); // (ORDER-1)*2

      io_mm.add( i_nQts,                      // m
                 l_nMdsEl,                    // n
                 i_nQts,                      // k
                 i_nQts,                      // ldA
                 l_nMdsEl,                    // ldB
                 l_nMdsEl,                    // ldC
                 static_cast<real_base>(1.0), // alpha
                 static_cast<real_base>(1.0), // beta
                 false,                       // fused AC
                 true,                        // fused BC
                 i_nCrs ); // (ORDER-1)*2+1

      // add first flux matrix
      io_mm.add( i_nQts,                      // m
                 l_nMdsFa,                    // n
                 l_nMdsEl,                    // k
                 l_nMdsEl,                    // ldA
                 l_nMdsFa,                    // ldB
                 l_nMdsFa,                    // ldC
                 static_cast<real_base>(1.0), // alpha
                 static_cast<real_base>(0.0), // beta
                 true,                        // fused AC
                 false,                       // fused BC
                 i_nCrs);                     // (ORDER-1)*2+2

      // add flux solver
      io_mm.add( i_nQts,                      // m
                 l_nMdsFa,                    // n
                 i_nQts,                      // k
                 i_nQts,                      // ldA
                 l_nMdsFa,                    // ldB
                 l_nMdsFa,                    // ldC
                 static_cast<real_base>(1.0), // alpha
                 static_cast<real_base>(0.0), // beta
                 false,                       // fused AC
                 true,                        // fused BC
                 i_nCrs );                    // (ORDER-1)*2+3

      // add second flux matrix
      io_mm.add( i_nQts,                      // m
                 l_nMdsEl,                    // n
                 l_nMdsFa,                    // k
                 l_nMdsFa,                    // ldA
                 l_nMdsEl,                    // ldB
                 l_nMdsEl,                    // ldC
                 static_cast<real_base>(1.0), // alpha
                 static_cast<real_base>(1.0), // beta
                 true,                        // fused AC
                 false,                       // fused BC
                 i_nCrs );                    // (ORDER-1)*2+4
    }
#endif

};

#endif
