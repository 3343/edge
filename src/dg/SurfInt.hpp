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
 * Data structures for ADER-DG surface integrations.
 **/
#ifndef EDGE_DG_SURF_INT_HPP
#define EDGE_DG_SURF_INT_HPP
#include "constants.hpp"
#include "data/Dynamic.h"

namespace edge {
  namespace dg {
      template< t_entityType   TL_T_EL,
                unsigned short TL_O_SP >
      class SurfInt;
  }
}

/**
 * Functions for the setup of surface integration data structures.
 *
 * @paramt TL_T_EL element type.
 * @paramt TL_O_SP spatial order.
 **/
template< t_entityType   TL_T_EL,
          unsigned short TL_O_SP >
class edge::dg::SurfInt {
  private:
    //! number of dimensions
    static unsigned short const TL_N_DIS = C_ENT[TL_T_EL].N_DIM;

    //! number of faces
    static unsigned short const TL_N_FAS = C_ENT[TL_T_EL].N_FACES;

    //! number of DG face modes
    static unsigned short const TL_N_MDS_FA = CE_N_ELEMENT_MODES( C_ENT[TL_T_EL].TYPE_FACES, TL_O_SP );

    //! number of DG element modes
    static unsigned short const TL_N_MDS_EL = CE_N_ELEMENT_MODES( TL_T_EL, TL_O_SP );

    //! number of neigboring contribution flux matrices
    static unsigned short const TL_N_FMNS = CE_N_FLUXN_MATRICES( TL_T_EL );

  private:
    /**
     * Stores the flux matrices as dense.
     * 
     * @param i_fIntL local flux matrices.
     * @param i_fIntN neighboring flux matrices.
     * @param i_fIntT transposed flux matrices.
     * @param io_dynMem dynamic memory management, which will be used for the respective allocations.
     * @param o_fIntLN will contain pointers to memory for the local and neighboring flux matrices.
     * @param o_fIntT will contain pointers to memory for the transposed flux matrices.
     *
     * @paramt TL_T_REAL real type.
     **/
    template< typename TL_T_REAL >
    static void storeFluxDense( TL_T_REAL     const   i_fIntL[TL_N_FAS][TL_N_MDS_EL][TL_N_MDS_FA],
                                TL_T_REAL     const   i_fIntN[TL_N_FMNS][TL_N_MDS_EL][TL_N_MDS_FA],
                                TL_T_REAL     const   i_fIntT[TL_N_FAS][TL_N_MDS_FA][TL_N_MDS_EL],
                                data::Dynamic       & io_dynMem,
                                TL_T_REAL           * o_fIntLN[TL_N_FAS+TL_N_FMNS],
                                TL_T_REAL           * o_fIntT[TL_N_FAS] ) {
      // allocate raw memory for the flux integration matrices
      std::size_t l_size  = TL_N_FAS  * std::size_t(TL_N_MDS_EL) * TL_N_MDS_FA;
                  l_size += TL_N_FMNS * std::size_t(TL_N_MDS_EL) * TL_N_MDS_FA;
                  l_size += TL_N_FAS  * std::size_t(TL_N_MDS_FA) * TL_N_MDS_EL;
                  l_size *= sizeof(TL_T_REAL);
      TL_T_REAL * l_fIntRaw = (TL_T_REAL*) io_dynMem.allocate( l_size,
                                                               4096,
                                                               false,
                                                               true );

      // local
      std::size_t l_en = 0;
      for( unsigned short l_fa = 0; l_fa < TL_N_FAS; l_fa++ ) {
        // set pointer
        o_fIntLN[l_fa] = l_fIntRaw + l_en;

        // set data
        for( unsigned short l_m0 = 0; l_m0 < TL_N_MDS_EL; l_m0++ ) {
          for( unsigned short l_m1 = 0; l_m1 < TL_N_MDS_FA; l_m1++ ) {
            l_fIntRaw[l_en] = i_fIntL[l_fa][l_m0][l_m1];
            l_en++;
          }
        }
      }

      // neighboring
      for( unsigned short l_ne = 0; l_ne < TL_N_FMNS; l_ne++ ) {
        // set pointer
        o_fIntLN[TL_N_FAS + l_ne] = l_fIntRaw + l_en;

        // set data
        for( unsigned short l_m0 = 0; l_m0 < TL_N_MDS_EL; l_m0++ ) {
          for( unsigned short l_m1 = 0; l_m1 < TL_N_MDS_FA; l_m1++ ) {
            l_fIntRaw[l_en] = i_fIntN[l_ne][l_m0][l_m1];
            l_en++;
          }
        }
      }

      // transposed
      for( unsigned short l_fa = 0; l_fa < TL_N_FAS; l_fa++ ) {
        // set pointer
        o_fIntT[l_fa] = l_fIntRaw + l_en;

        // set data
        for( unsigned short l_m0 = 0; l_m0 < TL_N_MDS_FA; l_m0++ ) {
          for( unsigned short l_m1 = 0; l_m1 < TL_N_MDS_EL; l_m1++ ) {
            l_fIntRaw[l_en] = i_fIntT[l_fa][l_m0][l_m1];
            l_en++;
          }
        }
      }
    }

  public:
    /**
     * Stores the flux matrices as dense.
     * 
     * @param io_dynMem dynamic memory management, which will be used for the respective allocations.
     * @param o_fIntLN will contain pointers to memory for the local and neighboring flux matrices.
     * @param o_fIntT will contain pointers to memory for the transposed flux matrices.
     *
     * @paramt TL_T_REAL real type.
     **/
    template< typename TL_T_REAL >
    static void storeFluxDense( data::Dynamic       & io_dynMem,
                                TL_T_REAL           * o_fIntLN[TL_N_FAS+TL_N_FMNS],
                                TL_T_REAL           * o_fIntT[TL_N_FAS] ) {
      // formulation of the basis in terms of the reference element
      dg::Basis l_basis( TL_T_EL,
                         TL_O_SP );

      // get flux matrices
      TL_T_REAL l_fIntL[TL_N_FAS][TL_N_MDS_EL][TL_N_MDS_FA];
      TL_T_REAL l_fIntN[TL_N_FMNS][TL_N_MDS_EL][TL_N_MDS_FA];
      TL_T_REAL l_fIntT[TL_N_FAS][TL_N_MDS_FA][TL_N_MDS_EL];
      l_basis.getFluxDense( l_fIntL[0][0],
                            l_fIntN[0][0],
                            l_fIntT[0][0] );

      // store flux matrices dense
      storeFluxDense( l_fIntL,
                      l_fIntN,
                      l_fIntT,
                      io_dynMem,
                      o_fIntLN,
                      o_fIntT );
    }
};

#endif
