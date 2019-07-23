/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2019, Alexander Breuer
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
 * Gets fake sparse matrices for kernel inits.
 **/
#ifndef EDGE_SEISMIC_KERNELS_FAKE_MATS_HPP
#define EDGE_SEISMIC_KERNELS_FAKE_MATS_HPP

#include "../common.hpp"
#include "../setups/ViscoElasticity.h"

namespace edge {
  namespace seismic {
    namespace kernels {
      template< unsigned short TL_N_DIS >
      class FakeMats;
    }
  }
}

template< unsigned short TL_N_DIS >
class edge::seismic::kernels::FakeMats {
  private:
    //! number of elastic quantities
    static unsigned short const TL_N_QTS_E = CE_N_QTS_E( TL_N_DIS );

    //! number of quantities per relaxation mechanism
    static unsigned short const TL_N_QTS_M = CE_N_QTS_M( TL_N_DIS );

  public:
    /**
     * Gets a fake elastic star matrix in CSR format.
     *
     * @param o_starCsrE will be set to star matrix.
     **/
    static void starCsrE( t_matCsr & o_starCsrE ) {
      // assemble fake star matrices
      double l_starE[TL_N_DIS][TL_N_QTS_E][TL_N_QTS_E];
      edge::seismic::common::getJac( 1.0,
                                     1.0,
                                     1.0,
                                     l_starE[0][0],
                                     TL_N_DIS );
      for( unsigned short l_di = 1; l_di < TL_N_DIS; l_di++ ) {
        for( int_qt l_q1 = 0; l_q1 < TL_N_QTS_E; l_q1++ ) {
          for( int_qt l_q2 = 0; l_q2 < TL_N_QTS_E; l_q2++ ) {
            l_starE[0][l_q1][l_q2] =   std::abs(l_starE[0][l_q1][l_q2])
                                     + std::abs(l_starE[l_di][l_q1][l_q2]);
          }
        }
      }

      // get csr star matrices
      edge::linalg::Matrix::denseToCsr< double >( TL_N_QTS_E,
                                                  TL_N_QTS_E,
                                                  l_starE[0][0],
                                                  o_starCsrE,
                                                  0 );
    }

    /**
     * Gets a fake viscoelastic star matrix in CSR format.
     * In contrast to the dense-version, all columns are considered in codIdx.
     *
     * @param o_starCsrA will be set to star matrix.
     **/
    static void starCsrA( t_matCsr & o_starCsrA ) {
      // assemble fake star matrices
      double l_starA[TL_N_DIS][TL_N_QTS_M * TL_N_DIS];
      double l_jacInv[TL_N_DIS][TL_N_DIS];
      for( unsigned short l_d0 = 0; l_d0 < TL_N_DIS; l_d0++ )
        for( unsigned short l_d1 = 0; l_d1 < TL_N_DIS; l_d1++ ) l_jacInv[l_d0][l_d1] = 1;

      setups::ViscoElasticity::star( l_jacInv, l_starA );

      // convert to csr
      edge::linalg::Matrix::denseToCsr< double >( TL_N_QTS_M,
                                                  TL_N_DIS,
                                                  l_starA[0],
                                                  o_starCsrA,
                                                  0 );

      // add zero columns
      for( unsigned short l_nz = 0; l_nz < o_starCsrA.colIdx.size(); l_nz++ ) {
        o_starCsrA.colIdx[l_nz] += TL_N_QTS_M;
      }
    }

    /**
     * Gets a fake viscoelastic source matrix in CSR format.
     *
     * @param o_srcCsrA will be set to star matrix.
     **/
    static void srcCsrA( t_matCsr & o_srcCsrA ) {
      // assemble fake source matrix
      double l_srcA[1][TL_N_QTS_M*TL_N_QTS_M];
      double l_lameEl[2];
      setups::ViscoElasticity::src( 1,
                                    2.5,
                                    100,
                                    60,
                                    30,
                                    32.4E9,
                                    32.4E9,
                                    l_lameEl[0],
                                    l_lameEl[1],
                                    l_srcA );

      // convert to csr
      edge::linalg::Matrix::denseToCsr< double >( TL_N_QTS_M,
                                                  TL_N_QTS_M,
                                                  l_srcA[0],
                                                  o_srcCsrA,
                                                  0 );
    }

    /**
     * Gets a fake elastic flux solver matrix in CSR format.
     *
     * @param o_fsCsrE will be set to flux solver.
     **/
    static void fsCsrE( t_matCsr & o_fsCsrE ) {
      // assemble fake solver
      double l_fsE[TL_N_QTS_E][TL_N_QTS_E];
      for( unsigned short l_q0 = 0; l_q0 < TL_N_QTS_E; l_q0++ ) {
        for( unsigned short l_q1 = 0; l_q1 < TL_N_QTS_E; l_q1++ ) {
          l_fsE[l_q0][l_q1] = 1;
        }
      }

      // convert to csr
      edge::linalg::Matrix::denseToCsr< double >( TL_N_QTS_E,
                                                  TL_N_QTS_E,
                                                  l_fsE[0],
                                                  o_fsCsrE,
                                                  0 );

    }

    /**
     * Gets a fake anelastic flux solver matrix in CSR format.
     *
     * @param o_fsCsrA will be set to flux solver.
     **/
    static void fsCsrA( t_matCsr & o_fsCsrA ) {
      // assemble fake solver
      double l_fsA[TL_N_QTS_M][TL_N_QTS_E];
      for( unsigned short l_q0 = 0; l_q0 < TL_N_QTS_M; l_q0++ ) {
        for( unsigned short l_q1 = 0; l_q1 < TL_N_QTS_E; l_q1++ ) {
          l_fsA[l_q0][l_q1] = 1;
        }
      }

      // convert to csr
      edge::linalg::Matrix::denseToCsr< double >( TL_N_QTS_M,
                                                  TL_N_QTS_E,
                                                  l_fsA[0],
                                                  o_fsCsrA,
                                                  0 );
    }

};

#endif