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
 * Initializes internal boundary data.
 **/
#ifndef EDGE_SC_IBND_INIT_HPP
#define EDGE_SC_IBND_INIT_HPP

#include "constants.hpp"
#include "InternalBoundary.type"
#include "data/SparseEntities.hpp"

namespace edge {
  namespace sc {
    namespace ibnd {
      template< t_entityType TL_T_EL >
      class Init;
    }
  }
}

/**
 * Initialization of internal boundaries.
 *
 * @paramt TL_T_EL element type.
 **/
template< t_entityType TL_T_EL >
class edge::sc::ibnd::Init {
  private:
    //! number of faces.
    static const unsigned short TL_N_FAS = C_ENT[TL_T_EL].N_FACES;

  public:
    /**
     * Initializes the connectivity info.
     *
     * @param i_nFa number of faces.
     * @param i_nEl number of elements.
     * @param i_spIb sparse type of the internal boudary.
     * @param i_spLi sparse type of limited elements.
     * @param i_faEl elements adjacent to faces (no bridge).
     * @param i_elFa faces adjacent to elements (no bridge).
     * @param i_charsFa face characteristics.
     * @param i_charsEl element characteristics.
     * @param o_conn will be set to connectivity info.
     *
     * @paramt TL_T_LID integral type for local ids.
     * @paramt TL_T_SP sparse type.
     * @paramt TL_T_CHARS_FA characteristics of DG-faces, offers member .spType.
     * @paramt TL_T_CHARS_EL characteristics of DG-elements, offers member .spType.
     **/
    template< typename     TL_T_LID,
              typename     TL_T_SP,
              typename     TL_T_CHARS_FA,
              typename     TL_T_CHARS_EL >
    static void connect( TL_T_LID                      i_nFa,
                         TL_T_LID                      i_nEl,
                         TL_T_SP                       i_spIb,
                         TL_T_SP                       i_spLi,
                         TL_T_LID                    (*i_faEl)[2],
                         TL_T_LID              const (*i_elFa)[ C_ENT[TL_T_EL].N_FACES ],
                         TL_T_CHARS_FA         const  *i_charsFa,
                         TL_T_CHARS_EL         const  *i_charsEl,
                         t_connect< TL_T_LID,
                                    TL_T_EL  >        &o_conn ) {
      // link internal boundary faces and internal boundary elements
      data::SparseEntities::linkSpAdjDst( i_nFa,
                                          2,
                                          i_faEl[0],
                                          i_spIb,
                                          i_spIb,
                                          i_charsFa,
                                          i_charsEl,
                                          o_conn.bfBe[0] );

      // link internal boundary faces and limited elements
      data::SparseEntities::linkSpAdjDst( i_nFa,
                                          2,
                                          i_faEl[0],
                                          i_spIb,
                                          i_spLi,
                                          i_charsFa,
                                          i_charsEl,
                                          o_conn.bfLe[0] );
                                          
      // link internal boundary elements and internal boundary faces
      data::SparseEntities::linkSpAdjSst( i_nEl,
                                          TL_N_FAS,
                                          i_elFa[0],
                                          i_spIb,
                                          i_charsFa,
                                          o_conn.beBf[0] );

      // link internal boundary faces and dense elements
      edge::data::SparseEntities::linkSpDe( i_nEl,
                                            i_spIb,
                                            i_charsEl,
                                            o_conn.beEl );
    }
};

#endif
