/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2017, Regents of the University of California
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
 * Sub-grid operations.
 **/
#ifndef SC_SUBGRID_HPP
#define SC_SUBGRID_HPP

#include "constants.hpp"

namespace edge {
  namespace sc {
    template< t_entityType   TL_T_EL,
              unsigned short TL_O_SP >
    class SubGrid;
  }
}

/**
 * Sub-grid oprations.
 *
 * @paramt TL_T_EL element type.
 * @paramt TL_O_SP spatial order of the DG-scheme.
 **/
template< t_entityType   TL_T_EL,
          unsigned short TL_O_SP >
class edge::sc::SubGrid {
  private:
    //! number of dimensions
    static unsigned short const TL_N_DIS = C_ENT[ TL_T_EL ].N_DIM;

    //! number of vertices per DG-element / sub-cell
    static unsigned short const TL_N_VES = C_ENT[TL_T_EL].N_VERTICES;

    //! number of sub-vertices
    static unsigned short const TL_N_SVS = CE_N_SUB_VERTICES( TL_T_EL, TL_O_SP );

    //! number of sub-cells per element
    static unsigned short const TL_N_SCS  = CE_N_SUB_CELLS( TL_T_EL, TL_O_SP );

  public:
    /**
     * Returns the closest sub-cell to a point.
     *
     * @param i_pt reference coordinates of the point.
     * @param i_scSv sub-vertices adjacent to the sub-cells (no bridge).
     * @param i_svChars characteristics of the sub-vertices.
     *
     * @paramt TL_T_LID integral type of local ids.
     * @paramt TL_T_REAL floating point type.
     * @paramt TL_T_SV_CHARS type of the sub-vertex characteristics, offers coordinates through .coords.
     **/
    template< typename TL_T_LID,
              typename TL_T_REAL,
              typename TL_T_SV_CHARS >
    static TL_T_LID ptSc( TL_T_REAL     const i_pt[TL_N_DIS],
                          TL_T_LID      const i_scSv[TL_N_SCS][TL_N_VES],
                          TL_T_SV_CHARS const i_svChars[TL_N_SVS] ) {
      // sub-cell which is clostest to the point
      TL_T_LID l_scM = std::numeric_limits< TL_T_LID >::max();

      // current minimum, squared distance (over all vertices)
      TL_T_REAL l_dxM = std::numeric_limits< TL_T_REAL >::max();

      // iterate over sub-cells
      for( TL_T_LID l_sc = 0; l_sc < TL_N_SCS; l_sc++ ) {
        // determine combined distance
        TL_T_REAL l_dx = 0;
        for( unsigned short l_ve = 0; l_ve < TL_N_VES; l_ve++ ) {
          TL_T_LID l_sv = i_scSv[l_sc][l_ve];
          for( unsigned short l_di = 0; l_di < TL_N_DIS; l_di++ ) {
            TL_T_REAL l_diff = i_svChars[l_sv].coords[l_di] - i_pt[l_di];
            l_dx += l_diff * l_diff;
          }
        }

        // update minima if possible
        if( l_dx < l_dxM ) {
          l_scM = l_sc;
          l_dxM = l_dx;
        }
      }

      return l_scM;
    }
};

#endif
