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
 * Allocates memory for sub-cell data.
 **/
#ifndef EDGE_SC_MEMORY_HPP
#define EDGE_SC_MEMORY_HPP

#include "constants.hpp"
#include "data/Dynamic.h"

namespace edge {
  namespace sc {
    class Memory;
  }
}

/**
 * Memory allocations for the sub-cell limiter.
 **/
class edge::sc::Memory {
  public:
    /**
     * Allocates memory for the sub-cell limiter.
     *
     * @param i_nLim number of limited DG elements.
     * @param i_nLimPlus number of limited plus DG elements (limited + face neighbors).
     * @param i_nExt number of elements providing extrema.
     * @param io_dynMem will be called for dynamic memory allocation.
     * @param o_lim sub-cell limiter, for which memory is allocated.
     *
     * @paramt TL_T_LID integral type of local ids.
     * @paramt TL_T_REAL floating point type.
     * @paramt TL_T_EL element type.
     * @paramt TL_O_SP spatial order of the DG-scheme.
     * @paramt TL_N_QTS number of quantities.
     * @paramt TL_N_CRS number of fused simulations.
     **/
    template< typename       TL_T_LID,
              typename       TL_T_REAL,
              t_entityType   TL_T_EL,
              unsigned short TL_O_SP,
              unsigned short TL_N_QTS,
              unsigned short TL_N_CRS >
    static void alloc( TL_T_LID              i_nLim,
                       TL_T_LID              i_nLimPlus,
                       TL_T_LID              i_nExt,
                       edge::data::Dynamic  &io_dynMem,
                       t_SubCell< TL_T_LID,
                                  TL_T_REAL,
                                  TL_T_EL,
                                  TL_O_SP,
                                  TL_N_QTS,
                                  TL_N_CRS > &o_lim ) {
      //! number of bytes for based pointer aligned of all DOF-related data structures
      const std::size_t TL_N_ALIGN = CE_MAX( TL_N_CRS * sizeof(TL_T_REAL), std::size_t(64) );

      //! number of sub-cells per element
      const unsigned short TL_N_SCS  = CE_N_SUB_CELLS( TL_T_EL, TL_O_SP );

      //! number of sub-faces per element face
      const unsigned short TL_N_SFS = CE_N_SUB_FACES( TL_T_EL, TL_O_SP );

      //! number of faces
      const unsigned short TL_N_FAS = C_ENT[TL_T_EL].N_FACES;

      // DOFs
      std::size_t l_dofsSize  = TL_N_QTS * (std::size_t) TL_N_SCS * TL_N_CRS;
                  l_dofsSize *= i_nLim * sizeof(TL_T_REAL);
      o_lim.dofs =  (TL_T_REAL (*) [TL_N_QTS][TL_N_SCS][TL_N_CRS]) io_dynMem.allocate( l_dofsSize, TL_N_ALIGN );

      // tDofs (this overestimates the memory requirements, as lp only requires selcted faces)
      std::size_t l_tDofsRawSize  = TL_N_FAS * TL_N_QTS * (std::size_t) TL_N_SFS * TL_N_CRS;
                  l_tDofsRawSize *= i_nLimPlus * sizeof(TL_T_REAL);
      for( unsigned short l_bu = 0; l_bu < 2; l_bu++ )
        o_lim.tDofsRaw[l_bu] = (TL_T_REAL (*) [TL_N_QTS][TL_N_SFS][TL_N_CRS]) io_dynMem.allocate( l_tDofsRawSize, TL_N_ALIGN );

      // size of the tDofs pointer structures
      std::size_t l_tDofsPtrSize  = TL_N_FAS;
                  l_tDofsPtrSize *= i_nLimPlus * sizeof(TL_T_REAL*);
      for( unsigned short l_bu = 0; l_bu < 2; l_bu++ )
        o_lim.tDofs[l_bu] =  (TL_T_REAL (* (*) [TL_N_FAS]) [TL_N_QTS][TL_N_SFS][TL_N_CRS]) io_dynMem.allocate( l_tDofsPtrSize );

      // admissibility of the DG solution
      std::size_t l_admSize =  TL_N_CRS;
                  l_admSize *= i_nLim * sizeof(bool);
      o_lim.adm[0] = (bool (*) [TL_N_CRS]) io_dynMem.allocate( l_admSize );
      o_lim.adm[1] = (bool (*) [TL_N_CRS]) io_dynMem.allocate( l_admSize );
      o_lim.adm[2] = (bool (*) [TL_N_CRS]) io_dynMem.allocate( l_admSize );
      o_lim.adm[3] = (bool (*) [TL_N_CRS]) io_dynMem.allocate( l_admSize );

      // lock for sub-cell solution
      std::size_t l_lockSize =  TL_N_CRS;
                  l_lockSize *= i_nLim * sizeof(bool);
      o_lim.lock = (bool (*) [TL_N_CRS]) io_dynMem.allocate( l_lockSize );

      // #limits since sync
      std::size_t l_limSyncSize  = TL_N_CRS;
                  l_limSyncSize *= i_nLim * sizeof(unsigned int);
      o_lim.limSync = (unsigned int (*) [TL_N_CRS]) io_dynMem.allocate( l_limSyncSize );

      // extrema of the DG / sub-cell solution
      std::size_t l_extSize  = 2 * (std::size_t) TL_N_QTS * TL_N_CRS;
                  l_extSize *= i_nExt * sizeof(TL_T_REAL);
      o_lim.ext[0] = (TL_T_REAL (*)[2][TL_N_QTS][TL_N_CRS]) io_dynMem.allocate( l_extSize, TL_N_ALIGN );
      o_lim.ext[1] = (TL_T_REAL (*)[2][TL_N_QTS][TL_N_CRS]) io_dynMem.allocate( l_extSize,TL_N_ALIGN );

      // link between between dominant limited possibly redundant (no bridge)
      std::size_t l_liDoLiDuSize  = TL_N_FAS;
                  l_liDoLiDuSize *= i_nLim * sizeof(TL_T_LID);
      o_lim.connect.liDoLiDu = (TL_T_LID (*)[TL_N_FAS]) io_dynMem.allocate( l_liDoLiDuSize );

      // link between limited and limited plus (no bridge)
      std::size_t l_liLpSize = i_nLim * sizeof(TL_T_LID);
      o_lim.connect.liLp = (TL_T_LID*) io_dynMem.allocate( l_liLpSize );

      // link between limited plus and dense DG elements (no bridge)
      std::size_t l_lpElSize  = i_nLimPlus * sizeof(TL_T_LID);
      o_lim.connect.lpEl = (TL_T_LID (*)) io_dynMem.allocate( l_lpElSize );

      // link between limited plus and limited elements (no bridge)
      std::size_t l_lpLiSize  = i_nLimPlus * sizeof(TL_T_LID);
      o_lim.connect.lpLi = (TL_T_LID (*)) io_dynMem.allocate( l_lpLiSize );

      // link between limited plus and limited DG elements (faces as bridge)
      std::size_t l_lpFaLpSize  = C_ENT[TL_T_EL].N_FACES;
                  l_lpFaLpSize *= i_nLimPlus * sizeof(TL_T_LID);
      o_lim.connect.lpFaLp = (TL_T_LID (*)[C_ENT[TL_T_EL].N_FACES]) io_dynMem.allocate( l_lpFaLpSize );

      // link between limited elements and extrema
      std::size_t l_liExSize = i_nLim * sizeof(TL_T_LID);
      o_lim.connect.liEx = (TL_T_LID*) io_dynMem.allocate( l_liExSize );

      // link between limited DG elements and those providing extrema (vertices as bridge)
      std::size_t l_liVeExSize = (i_nLim+1) * sizeof(TL_T_LID*);
      o_lim.connect.liVeEx = (TL_T_LID**) io_dynMem.allocate( l_liVeExSize );
    }
};

#endif
