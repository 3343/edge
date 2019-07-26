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
 * Data structures for ADER time predictions.
 **/

#ifndef EDGE_DG_TIME_PRED_HPP
#define EDGE_DG_TIME_PRED_HPP

#include "constants.hpp"
#include "data/Dynamic.h"
#include "Basis.h"

namespace edge {
  namespace dg {
    template< t_entityType   TL_T_EL,
              unsigned short TL_O_SP,
              unsigned short TL_O_TI >
    class TimePred;
  }
}

/**
 * ADER-related functions.
 *
 * @paramt TL_T_EL element type.
 * @paramt TL_O_SP order in space.
 * @paramt TL_O_TI order in time.
 **/
template< t_entityType   TL_T_EL,
          unsigned short TL_O_SP,
          unsigned short TL_O_TI >
class edge::dg::TimePred {
  private:
    //! dimension of the element
    static unsigned short const TL_N_DIS = C_ENT[TL_T_EL].N_DIM;

    //! number of element modes
    static unsigned short const TL_N_MDS = CE_N_ELEMENT_MODES( TL_T_EL, TL_O_SP );

  private:
    /**
     * Stores the transposed stiffness matrices as dense.
     * This includes multiplications with (-1) for kernels with support for alpha==1 only.
     * 
     * @param i_stiffT dense stiffness matrices.
     * @param io_dynMem dynamic memory management, which will be used for the respective allocations.
     * @param o_stiffT will contain pointers to memory for the individual matrices.
     *
     * @paramt TL_T_REAL real type.
     **/
    template< typename TL_T_REAL >
    static void storeStiffTDense( TL_T_REAL     const     i_stiffT[TL_N_DIS][TL_N_MDS][TL_N_MDS],
                                  data::Dynamic         & io_dynMem,
                                  TL_T_REAL             * o_stiffT[CE_MAX(TL_O_TI-1,1)][TL_N_DIS] ) {
      // allocate raw memory for the stiffness matrices
      std::size_t l_size  = TL_N_DIS * std::size_t(TL_N_MDS) * TL_N_MDS;
                  l_size *= sizeof(TL_T_REAL);
      TL_T_REAL* l_stiffTRaw = (TL_T_REAL*) io_dynMem.allocate( l_size,
                                                                4096,
                                                                false,
                                                                true );

      // copy data
      std::size_t l_en = 0;
      for( unsigned short l_di = 0; l_di < TL_N_DIS; l_di++ ) {
        for( unsigned short l_m0 = 0; l_m0 < TL_N_MDS; l_m0++ ) {
          for( unsigned short l_m1 = 0; l_m1 < TL_N_MDS; l_m1++ ) {
            l_stiffTRaw[l_en] = -i_stiffT[l_di][l_m0][l_m1];
            l_en++;
          }
        }
      }

      // assign pointers
      for( unsigned short l_de = 0; l_de < TL_O_TI-1; l_de++ ) {
        for( unsigned short l_di = 0; l_di < TL_N_DIS; l_di++ ) {
          o_stiffT[l_de][l_di] = l_stiffTRaw + l_di * std::size_t(TL_N_MDS) * TL_N_MDS; 
        }
      }
    }

  public:
    /**
     * Stores the transposed stiffness matrices as dense.
     * This includes multiplications with (-1) for kernels with support for alpha==1 only.
     * 
     * @param io_dynMem dynamic memory management, which will be used for the respective allocations.
     * @param o_stiffT will contain pointers to memory for the individual matrices.
     *
     * @paramt TL_T_REAL real type.
     **/
    template< typename TL_T_REAL >
    static void storeStiffTDense( data::Dynamic & io_dynMem,
                                  TL_T_REAL     * o_stiffT[CE_MAX(TL_O_TI-1,1)][TL_N_DIS] ) {
      // formulation of the basis in terms of the reference element
      dg::Basis l_basis( TL_T_EL,
                         TL_O_SP );

      // get stiffness matrices
      TL_T_REAL l_stiffT[TL_N_DIS][TL_N_MDS][TL_N_MDS];
      l_basis.getStiffMm1Dense( TL_N_MDS,
                                l_stiffT[0][0],
                                true );

      // store stiffness matrices dense
      storeStiffTDense( l_stiffT,
                        io_dynMem,
                        o_stiffT );
    }
};

#endif
