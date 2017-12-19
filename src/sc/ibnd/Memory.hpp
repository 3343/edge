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
 * Allocates memory for internal boudary data.
 **/
#ifndef EDGE_SC_IBND_MEMORY_HPP
#define EDGE_SC_IBND_MEMORY_HPP

#include "data/Dynamic.h"

namespace edge {
  namespace sc {
    namespace ibnd {
      class Memory;
    }
  }
}

/**
 * Memory allocations for the sub-cell limited internal boundaries.
 **/
class edge::sc::ibnd::Memory {
  public:
    /**
     * Allocates memory for the sub-cell limited internal boundaries.
     *
     * @param i_nBf number of faces at the internal boundary.
     * @param i_nBe number of elements (faces as bridge) at the internal boundary.
     * @param io_dynMem will be sued for dynamic memory allocations.
     * @param o_iBnd internal boundary for which data will memory will be allocated.
     *
     * @paramt TL_T_LID integral type for local ids.
     * @paramt TL_T_REAL floating point type.
     * @paramt TL_T_SP sparse integral type.
     * @paramt TL_T_EL element type.
     * @paramt TL_N_QTS number of quantities.
     **/
    template< typename       TL_T_LID,
              typename       TL_T_REAL,
              typename       TL_T_SP,
              t_entityType   TL_T_EL,
              unsigned short TL_N_QTS >
    static void alloc( TL_T_LID                        i_nBf,
                       TL_T_LID                        i_nBe,
                       edge::data::Dynamic            &io_dynMem,
                       t_InternalBoundary< TL_T_LID,
                                           TL_T_REAL,
                                           TL_T_SP,
                                           TL_T_EL,
                                           TL_N_QTS > &o_iBnd ) {
      //! number of faces per element
      static unsigned short const TL_N_FAS = C_ENT[TL_T_EL].N_FACES;


      // ibnd fas adjacent to limited els
      std::size_t l_bfLeSize = std::size_t(2) * i_nBf * sizeof(TL_T_LID);
      o_iBnd.connect.bfLe = (TL_T_LID (*)[2]) io_dynMem.allocate( l_bfLeSize );

      // ibnd els adjacent to ibnd fas
      std::size_t l_bfBeSize = std::size_t(2) * i_nBf * sizeof(TL_T_LID);
      o_iBnd.connect.bfBe = (TL_T_LID (*)[2]) io_dynMem.allocate( l_bfBeSize );

      // ibnd els adjacent to ibnd fas
      std::size_t l_beBfSize = std::size_t(TL_N_FAS) * i_nBe * sizeof(TL_T_LID);
      o_iBnd.connect.beBf = (TL_T_LID (*)[TL_N_FAS]) io_dynMem.allocate( l_beBfSize );

      // link between dense els and ibnd els
      std::size_t l_beElSize = i_nBe * sizeof(TL_T_LID);
      o_iBnd.connect.beEl = (TL_T_LID *) io_dynMem.allocate( l_beElSize ); 


      // ibnd fa chars
      std::size_t l_bfCharsSize = i_nBf * sizeof( t_bfChars< TL_T_SP > );
      o_iBnd.bfChars = (t_bfChars< TL_T_SP > *) io_dynMem.allocate( l_bfCharsSize );


      // middle state solvers
      std::size_t l_mssSize = i_nBf * std::size_t(4) * TL_N_QTS * TL_N_QTS * sizeof(TL_T_REAL);
      o_iBnd.mss = (TL_T_REAL (*)[4][TL_N_QTS][TL_N_QTS]) io_dynMem.allocate( l_mssSize );
    }
};

#endif
