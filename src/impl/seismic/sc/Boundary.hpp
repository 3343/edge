/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2018, Regents of the University of California
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
 * Boundary conditions for the seismic sub-cell solvers.
 **/
#ifndef EDGE_SEISMIC_SC_BOUNDARY_HPP
#define EDGE_SEISMIC_SC_BOUNDARY_HPP

#include "constants.hpp"
#include "io/logging.h"

namespace edge {
  namespace elastic {
    namespace sc {
      template< t_entityType   TL_T_EL,
                unsigned short TL_O_SP,
                unsigned short TL_N_CRS >
      class Boundary;
    }
  }
}

/**
 * Boundary conditions for the sub-cell solvers.
 *
 * @paramt TL_T_EL element type.
 * @paramt TL_O_SP spatial order of the DG-solver.
 * @paramt TL_N_CRS number of fused simulations.
 **/
template< t_entityType   TL_T_EL,
          unsigned short TL_O_SP,
          unsigned short TL_N_CRS >
class edge::elastic::sc::Boundary {
  private:
    //! number of dimensions
    static unsigned short const TL_N_DIS = C_ENT[TL_T_EL].N_DIM;

    //! number of faces per element
    static const unsigned short TL_N_FAS = C_ENT[TL_T_EL].N_FACES;

    //! number of elastic quantities
    static const unsigned short TL_N_QTS = (TL_N_DIS == 2) ? 5 : 9;

    //! number of sub-faces per element face
    static unsigned short const TL_N_SFS = CE_N_SUB_FACES( TL_T_EL, TL_O_SP );

    //! number of sub-cells per element
    static unsigned short const TL_N_SCS = CE_N_SUB_CELLS( TL_T_EL, TL_O_SP );

  public:
    /**
     * Applies boundary conditions for the elastic wave equations to a sub-grid.
     *
     * @param i_scSfSc sub-cells adjacent to sub-cells (sub-faces as bridge).
     * @param i_spTypes sparse types of the DG-faces.
     * @param io_dofsSc sub-grid for which the ghost-sub-cells with boundary conditions are overwritten.
     *
     * @paramt TL_T_SP integral type for sparse ids.
     * @paramt TL_T_REAL floating point precision.
     **/
    template< typename TL_T_SP,
              typename TL_T_REAL >
    static void apply( unsigned short const i_scSfSc[TL_N_SCS + TL_N_FAS * TL_N_SFS][TL_N_FAS],
                       TL_T_SP        const i_spTypes[TL_N_FAS],
                       TL_T_REAL            io_dofsSc[TL_N_QTS][TL_N_SCS+TL_N_FAS*TL_N_SFS][TL_N_CRS] ) {
      // iterate over faces
      for( unsigned short l_fa = 0; l_fa < TL_N_FAS; l_fa++ ) {
        // default: no boundary condition
        if(    (i_spTypes[l_fa] & FREE_SURFACE) != FREE_SURFACE
            && (i_spTypes[l_fa] & OUTFLOW     ) != OUTFLOW ){}
        else if( (i_spTypes[l_fa] & FREE_SURFACE) == FREE_SURFACE ) {
          // set free surface ghost-data
          for( unsigned short l_sf = 0; l_sf < TL_N_SFS; l_sf++ ) {
            unsigned short l_sc = i_scSfSc[TL_N_SCS + TL_N_SFS*l_fa + l_sf][0];
            for( unsigned short l_qt = 0; l_qt < TL_N_QTS; l_qt++ )
              for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ )
                io_dofsSc[l_qt][TL_N_SCS + TL_N_SFS*l_fa + l_sf][l_cr] = io_dofsSc[l_qt][l_sc][l_cr];
          }
        }
        else if( (i_spTypes[l_fa] & OUTFLOW ) == OUTFLOW ) {
          for( unsigned short l_sf = 0; l_sf < TL_N_SFS; l_sf++ )
            for( unsigned short l_qt = 0; l_qt < TL_N_QTS; l_qt++ )
              for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ )
                io_dofsSc[l_qt][TL_N_SCS + TL_N_SFS*l_fa + l_sf][l_cr] = 0;
        }
      }
    }
};

#endif
