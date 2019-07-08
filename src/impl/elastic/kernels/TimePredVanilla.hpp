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
 * Time predictions through the ADER scheme for seismic setups using vanilla kernels.
 **/
#ifndef EDGE_SEISMIC_KERNELS_TIME_PRED_VANILLA_HPP
#define EDGE_SEISMIC_KERNELS_TIME_PRED_VANILLA_HPP

#include "TimePred.hpp"
#include "data/MmVanilla.hpp"

namespace edge {
  namespace elastic {
    namespace kernels {
      template< typename       TL_T_REAL,
                t_entityType   TL_T_EL,
                unsigned short TL_O_SP,
                unsigned short TL_O_TI,
                unsigned short TL_N_CRS >
      class TimePredVanilla;
    }
  }
}

/**
 * Vanilla ADER time predictions for seismic configurations.
 * 
 * @paramt TL_T_REAL floating point precision.
 * @paramt TL_T_EL element type.
 * @paramt TL_O_SP order in space.
 * @paramt TL_O_TI order in time.
 * @paramt TL_N_CRS number of fused forward simulations.
 **/
template< typename       TL_T_REAL,
          t_entityType   TL_T_EL,
          unsigned short TL_O_SP,
          unsigned short TL_O_TI,
          unsigned short TL_N_CRS >
class edge::elastic::kernels::TimePredVanilla: public edge::elastic::kernels::TimePred < TL_T_REAL,
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

    //! matrix kernels
    edge::data::MmVanilla< TL_T_REAL > m_mm;

    //! pointers to the (possibly recursive) stiffness matrices
    TL_T_REAL *m_stiffT[CE_MAX(TL_O_SP-1,1)][TL_N_DIS];

    /**
     * Generates the matrix kernels for the transposed stiffness matrices and star matrices.
     **/
    void generateKernels() {
      // (O-1)*2 kernels for time integration
      for( unsigned short l_de = 1; l_de < TL_O_SP; l_de++ ) {
        // multiplication with transpose stiffness matrix
        m_mm.add( 0,                                                 // group
                  TL_N_QTS_E,                                        // m
                  CE_N_ELEMENT_MODES_CK( TL_T_EL, TL_O_SP, l_de   ), // n
                  CE_N_ELEMENT_MODES_CK( TL_T_EL, TL_O_SP, l_de-1 ), // k
                  TL_N_MDS,                                          // ldA
                  TL_N_MDS,                                          // ldB
                  TL_N_MDS,                                          // ldC
                  static_cast<TL_T_REAL>(1.0),                       // alpha
                  static_cast<TL_T_REAL>(0.0),                       // beta
                  true,                                              // fused AC
                  false,                                             // fused BC
                  TL_N_CRS );

        // multiplication with star matrix
        m_mm.add( 1,                                               // group
                  TL_N_QTS_E,                                      // m
                  CE_N_ELEMENT_MODES_CK( TL_T_EL, TL_O_SP, l_de ), // n
                  TL_N_QTS_E,                                      // k
                  TL_N_QTS_E,                                      // ldA
                  TL_N_MDS,                                        // ldB
                  TL_N_MDS,                                        // ldC
                  static_cast<TL_T_REAL>(1.0),                     // alpha
                  static_cast<TL_T_REAL>(1.0),                     // beta
                  false,                                           // fused AC
                  true,                                            // fused BC
                  TL_N_CRS );
      }
    }

  public:
    /**
     * Constructor of the vanilla time prediction.
     * 
     * @param io_dynMem dynamic memory allocations.
     **/
    TimePredVanilla( data::Dynamic & io_dynMem ) {
      // formulation of the basis in terms of the reference element
      dg::Basis l_basis( TL_T_EL,
                         TL_O_SP );

      // get stiffness matrices
      TL_T_REAL l_stiffT[TL_N_DIS][TL_N_MDS][TL_N_MDS];
      l_basis.getStiffMm1Dense( TL_N_MDS,
                                l_stiffT[0][0],
                                true );


      // store stiffness matrices dense
      this->storeStiffTDense( l_stiffT,
                              io_dynMem,
                              m_stiffT );

      // generate kernels
      generateKernels();
    };

    /**
     * Applies the Cauchy–Kowalevski procedure (vanilla implementation) and computes time derivatives and time integrated DOFs.
     *
     * @param i_dT time step.
     * @param i_star star matrices.
     * @param i_dofs DOFs.
     * @param o_scratch will be used as scratch memory.
     * @param o_der will be set to time derivatives.
     * @param o_tInt will be set to time integrated DOFs.
     **/
    void ck( TL_T_REAL       i_dT,
             TL_T_REAL const i_star[TL_N_DIS][TL_N_QTS_E][TL_N_QTS_E],
             TL_T_REAL const i_dofs[TL_N_QTS_E][TL_N_MDS][TL_N_CRS],
             TL_T_REAL       o_scratch[TL_N_QTS_E][TL_N_MDS][TL_N_CRS],
             TL_T_REAL       o_der[TL_O_TI][TL_N_QTS_E][TL_N_MDS][TL_N_CRS],
             TL_T_REAL       o_tInt[TL_N_QTS_E][TL_N_MDS][TL_N_CRS] ) const {
      // scalar for the time integration
      TL_T_REAL l_scalar = i_dT;

      // initialize zero-derivative, reset time integrated dofs
      for( unsigned short l_qt = 0; l_qt < TL_N_QTS_E; l_qt++ ) {
        for( unsigned short l_md = 0; l_md < TL_N_MDS; l_md++ ) {
          for( unsigned short l_cfr = 0; l_cfr < TL_N_CRS; l_cfr++ ) {
            o_der[0][l_qt][l_md][l_cfr] = i_dofs[l_qt][l_md][l_cfr];
            o_tInt[l_qt][l_md][l_cfr]   = l_scalar * i_dofs[l_qt][l_md][l_cfr];
          }
        }
      }

      // iterate over time derivatives
      for( unsigned short l_de = 1; l_de < TL_O_TI; l_de++ ) {
        // reset this derivative
        for( unsigned short l_qt = 0; l_qt < TL_N_QTS_E; l_qt++ )
          for( unsigned short l_md = 0; l_md < TL_N_MDS; l_md++ )
            for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) o_der[l_de][l_qt][l_md][l_cr] = 0;

        // compute the derivatives
        for( unsigned short l_di = 0; l_di < TL_N_DIS; l_di++ ) {
          // multiply with transposed stiffness matrices and inverse mass matrix
          m_mm.m_kernels[0][l_de-1]( o_der[l_de-1][0][0],
                                     m_stiffT[l_de-1][l_di],
                                     o_scratch[0][0] );
          // multiply with star matrices
          m_mm.m_kernels[1][l_de-1]( i_star[l_di][0],
                                     o_scratch[0][0],
                                     o_der[l_de][0][0] );
        }

        // update scalar
        l_scalar *= -i_dT / (l_de+1);

        // update time integrated dofs
        for( unsigned short l_qt = 0; l_qt < TL_N_QTS_E; l_qt++ ) {
          for( unsigned short l_md = 0; l_md < CE_N_ELEMENT_MODES_CK( TL_T_EL, TL_O_SP, l_de ); l_md++ ) {
            for( int_cfr l_cr = 0; l_cr < TL_N_CRS; l_cr++ )
              o_tInt[l_qt][l_md][l_cr] += l_scalar * o_der[l_de][l_qt][l_md][l_cr];
          }
        }
      }
    }
};

#endif