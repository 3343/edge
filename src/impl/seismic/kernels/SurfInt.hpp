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
#include "dg/SurfInt.hpp"

namespace edge {
  namespace seismic {
    namespace kernels {
      template< typename       TL_T_REAL,
                unsigned short TL_N_RMS,
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
 * @paramt TL_N_RMS number of relaxation mechanisms.
 * @paramt TL_T_EL element type.
 * @paramt TL_O_SP spatial order.
 * @paramt TL_N_CRS number of fused simulations.
 **/
template< typename       TL_T_REAL,
          unsigned short TL_N_RMS,
          t_entityType   TL_T_EL,
          unsigned short TL_O_SP,
          unsigned short TL_N_CRS >
class edge::seismic::kernels::SurfInt {
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

    //! number of quantities per relaxation mechanism
    static unsigned short const TL_N_QTS_M = CE_N_QTS_M( TL_N_DIS );

    //! relaxation frequencies
    TL_T_REAL * m_rfs = nullptr;

  protected:
    /**
     * Stores the flux matrices as dense.
     *
     * @param io_dynMem dynamic memory management, which will be used for the respective allocations.
     * @param o_fIntLN will contain pointers to memory for the local and neighboring flux matrices.
     * @param o_fIntT will contain pointers to memory for the transposed flux matrices.
     **/
    static void storeFluxDense( data::Dynamic       & io_dynMem,
                                TL_T_REAL           * o_fIntLN[TL_N_FAS+TL_N_FMNS],
                                TL_T_REAL           * o_fIntT[TL_N_FAS] ) {
      dg::SurfInt< TL_T_EL,
                   TL_O_SP >::storeFluxDense( io_dynMem,
                                              o_fIntLN,
                                              o_fIntT );
    }

    /**
     * Constructor of the surface integration.
     *
     * @param i_rfs relaxation frequencies, use nullptr if TL_N_RMS==0.
     * @param io_dynMem dynamic memory allocations.
     **/
    SurfInt( TL_T_REAL     const * i_rfs,
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

  public:
    /**
     * Scatters the anelastic update (excluding frequency scaling) to the anelastic DOFs.
     *
     * @param i_update update, which is scattered.
     * @param io_dofsA anelastic DOFs, which are updated.
     **/
    void scatterUpdateA( TL_T_REAL const   i_update[TL_N_QTS_M][TL_N_MDS_EL][TL_N_CRS],
                         TL_T_REAL       (*io_dofsA)[TL_N_QTS_M][TL_N_MDS_EL][TL_N_CRS] ) const {
      for( unsigned short l_rm = 0; l_rm < TL_N_RMS; l_rm++ ) {
        for( unsigned short l_qt = 0; l_qt < TL_N_QTS_M; l_qt++ ) {
          for( unsigned short l_md = 0; l_md < TL_N_MDS_EL; l_md++ ) {
            for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ){
              io_dofsA[l_rm][l_qt][l_md][l_cr] += m_rfs[l_rm] * i_update[l_qt][l_md][l_cr];
            }
          }
        }
      }
    }
};

#endif
