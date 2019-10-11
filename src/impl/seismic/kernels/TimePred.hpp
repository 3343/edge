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
 * Time predictions through the ADER scheme for seismic setups.
 **/

#ifndef EDGE_SEISMIC_KERNELS_TIME_PRED_HPP
#define EDGE_SEISMIC_KERNELS_TIME_PRED_HPP

#include "constants.hpp"
#include "data/Dynamic.h"
#include "dg/TimePred.hpp"

namespace edge {
  namespace seismic {
    namespace kernels {
      template< typename       TL_T_REAL,
                unsigned short TL_N_RMS,
                t_entityType   TL_T_EL,
                unsigned short TL_O_SP,
                unsigned short TL_O_TI,
                unsigned short TL_N_CRS >
      class TimePred;
    }
  }
}

/**
 * ADER-related functions:
 *   1) Computation of time predictions (time derivatives and time integrated DOFs) through the Cauchy–Kowalevski procedure.
 *   2) Evaluation of time prediction at specific points in time.
 *
 * @paramt TL_T_REAL floating point precision of the solver.
 * @paramt TL_N_RMS number of relaxation mechanisms.
 * @paramt TL_T_EL element type.
 * @paramt TL_O_SP order in space.
 * @paramt TL_O_TI order in time.
 * @paramt TL_N_CRS number of concurrent forward runs (fused simulations)
 **/
template< typename       TL_T_REAL,
          unsigned short TL_N_RMS,
          t_entityType   TL_T_EL,
          unsigned short TL_O_SP,
          unsigned short TL_O_TI,
          unsigned short TL_N_CRS >
class edge::seismic::kernels::TimePred {
  private:
    //! dimension of the element
    static unsigned short const TL_N_DIS = C_ENT[TL_T_EL].N_DIM;

    //! number of element modes
    static unsigned short const TL_N_MDS = CE_N_ELEMENT_MODES( TL_T_EL, TL_O_SP );

    //! number of elastic quantities
    static unsigned short const TL_N_QTS_E = CE_N_QTS_E( TL_N_DIS );

  protected:
    //! relaxation frequencies
    TL_T_REAL * m_rfs = nullptr;

    /**
     * Stores the transposed stiffness matrices as dense.
     * 
     * @param io_dynMem dynamic memory management, which will be used for the respective allocations.
     * @param o_stiffT will contain pointers to memory for the individual matrices.
     **/
    static void storeStiffTDense( data::Dynamic         & io_dynMem,
                                  TL_T_REAL             * o_stiffT[CE_MAX(TL_O_SP-1,1)][TL_N_DIS] ) {
      dg::TimePred< TL_T_EL,
                    TL_O_SP,
                    TL_O_TI >::storeStiffTDense( io_dynMem,
                                                 o_stiffT );
    }

  public:
    /**
     * Constructor of the time prediction.
     *
     * @param i_rfs relaxation frequencies, use nullptr if TL_N_RMS==0.
     * @param io_dynMem dynamic memory allocations.
     **/
    TimePred( TL_T_REAL     const * i_rfs,
              data::Dynamic       & io_dynMem ) {
      if( TL_N_RMS > 0 ) {
        std::size_t l_size = TL_N_RMS * sizeof(TL_T_REAL);
        m_rfs = (TL_T_REAL *) io_dynMem.allocate( l_size );

        for( unsigned short l_rm = 0; l_rm < TL_N_RMS; l_rm++ ) {
          m_rfs[l_rm] = i_rfs[l_rm];
        }
      }
    }

    /**
     * Evaluates the time prediction, given by the time derivatives at the given points in time.
     * The points are relative to the time at which the time prediction was obtained.
     * Example:
     *   0    1.5  2.0  2.4 2.9 absolute time
     *   |-----|----x----x---x----------------->
     *         |   0.5  0.9 1.4 relative time (expected as input)
     *       time
     *    prediction
     *
     * @param i_nPts number of points.
     * @param i_pts relative pts in time.
     * @param i_der time prediction given through the time derivatives.
     * @param o_preDofs will be set to the predicted DOFs at the points in time.
     **/
    static void inline evalTimePrediction( unsigned short   i_nPts,
                                           TL_T_REAL const *i_pts,
                                           TL_T_REAL const  i_der[TL_O_TI][TL_N_QTS_E][TL_N_MDS][TL_N_CRS],
                                           TL_T_REAL       (*o_preDofs)[TL_N_QTS_E][TL_N_MDS][TL_N_CRS] ) {
      for( unsigned short l_pt = 0; l_pt < i_nPts; l_pt++ ) {
        // init dofs
        for( unsigned short l_qt = 0; l_qt < TL_N_QTS_E; l_qt++ )
          for( unsigned short l_md = 0; l_md < TL_N_MDS; l_md++ )
            for( unsigned short l_ru = 0; l_ru < TL_N_CRS; l_ru++ ) o_preDofs[l_pt][l_qt][l_md][l_ru] = i_der[0][l_qt][l_md][l_ru];

        // evaluate time derivatives at given point in time
        TL_T_REAL l_scalar = 1.0;

        // iterate over derivatives
        for( unsigned short l_de = 1; l_de < TL_O_TI; l_de++ ) {
          // update scalar
          l_scalar *= i_pts[l_pt] / l_de;

          for( unsigned short l_qt = 0; l_qt < TL_N_QTS_E; l_qt++ ) {
            for( unsigned short l_md = 0; l_md < CE_N_ELEMENT_MODES_CK( TL_T_EL, TL_O_SP, l_de ); l_md++ ) {
              for( unsigned short l_ru = 0; l_ru < TL_N_CRS; l_ru++ )
                o_preDofs[l_pt][l_qt][l_md][l_ru] += l_scalar * i_der[l_de][l_qt][l_md][l_ru];
            }
          }
        }
      }
    }

    /**
     * Integrates the DOFs over the interval [0, dt].
     *
     * @param i_dt time step.
     * @param i_derE elastic derivatives.
     * @param o_tIntE will be set to elastic time-integrated DOFs.
     **/
    static void integrate( TL_T_REAL       i_dt,
                           TL_T_REAL const i_derE[TL_O_TI][TL_N_QTS_E][TL_N_MDS][TL_N_CRS],
                           TL_T_REAL       o_tIntE[TL_N_QTS_E][TL_N_MDS][TL_N_CRS] ) {
      // scalar for the time integration
      TL_T_REAL l_sca = i_dt;

      // init time integrated dofs
      for( unsigned short l_qt = 0; l_qt < TL_N_QTS_E; l_qt++ )
        for( unsigned short l_md = 0; l_md < TL_N_MDS; l_md++ )
          for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ )
            o_tIntE[l_qt][l_md][l_cr] = l_sca * i_derE[0][l_qt][l_md][l_cr];

      // iterate over time derivatives
      for( unsigned int l_de = 1; l_de < TL_O_TI; l_de++ ) {
        // update scalar
        l_sca *= i_dt / (l_de+1);

        // elastic: update time integrated DOFs
        unsigned short l_nCpMds = (TL_N_RMS == 0) ? CE_N_ELEMENT_MODES_CK( TL_T_EL, TL_O_SP, l_de ) : TL_N_MDS;

        for( unsigned short l_qt = 0; l_qt < TL_N_QTS_E; l_qt++ )
          for( unsigned short l_md = 0; l_md < l_nCpMds; l_md++ )
            for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ )
              o_tIntE[l_qt][l_md][l_cr] += l_sca * i_derE[l_de][l_qt][l_md][l_cr];
      }
    }

};

#endif
