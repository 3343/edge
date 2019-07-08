/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *         Alexander Heinecke (alexander.heinecke AT intel.com)
 *
 * @section LICENSE
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
 * Quadrature-free ADER-DG surface integration for seismic wave propagation.
 **/
#ifndef EDGE_SEISMIC_KERNELS_SURF_INT_HPP
#define EDGE_SEISMIC_KERNELS_SURF_INT_HPP
#include "constants.hpp"
#include "data/Dynamic.h"

namespace edge {
  namespace elastic {
    namespace kernels {
      template< typename       TL_T_REAL,
                t_entityType   TL_T_EL,
                unsigned short TL_O_SP,
                unsigned short TL_N_CRS >
      class SurfInt;
    }
  }
}

/**
 * Quadrature-free ADER-DG surface integration for seismic wave propagation.
 *
 * @paramt TL_T_REAL floating point precision.
 * @paramt TL_T_EL element type.
 * @paramt TL_O_SP spatial order.
 * @paramt TL_N_CRS number of fused simulations.
 **/
template< typename       TL_T_REAL,
          t_entityType   TL_T_EL,
          unsigned short TL_O_SP,
          unsigned short TL_N_CRS >
class edge::elastic::kernels::SurfInt {
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

  protected:
    /**
     * Stores the flux matrices as dense.
     * 
     * @param i_fIntL local flux matrices.
     * @param i_fIntN neighboring flux matrices.
     * @param i_fIntT transposed flux matrices.
     * @param io_dynMem dynamic memory management, which will be used for the respective allocations.
     * @param o_fIntLN will contain pointers to memory for the local and neighboring flux matrices.
     * @param o_fIntT will contain pointers to memory for the transposed flux matrices.
     **/
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

    /**
     * Determines the flux matrix id for neighboring contribution of the quadrature-free face integral.
     *
     * @param i_vIdElFaEl id of the vertex, matching the element's vertex 0, from the perspective of the adjacent element w.r.t. to the reference element.
     * @param i_fIdElFaEl id of the face from the perspective of the adjacent element w.r.t. to the reference element.
     *
     * @return flux matrix id.
     **/
    static unsigned short inline fMatId( unsigned short i_vIdElFaEl,
                                         unsigned short i_fIdElFaEl ) {
      return i_vIdElFaEl*TL_N_FAS+i_fIdElFaEl;
    }
};

#endif
