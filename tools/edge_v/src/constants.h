/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 * @author Rajdeep Konwar (rkonwar AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2020, Friedrich Schiller University Jena
 * Copyright (c) 2019-2020, Alexander Breuer
 * Copyright (c) 2017-2018, Regents of the University of California
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
 * This file contains definition of compile time constants for Edge-V.
 **/

#ifndef EDGE_V_CONSTANTS_H
#define EDGE_V_CONSTANTS_H

#include <cstddef>

namespace edge_v {
  // integral type for ids
  typedef std::size_t t_idx;

  // sparse type
  typedef int t_sparseType;

  // entity types
  typedef enum {
    POINT    = 0,
    LINE     = 1,
    QUAD4R   = 2,
    TRIA3    = 3,
    HEX8R    = 4,
    TET4     = 5,
    INVALID  = 99
  } t_entityType;

  /**
   * Gets the number of dimensions for the given entity type.
   *
   * @param i_enTy entity type.
   * @return number of dimensions.
   **/
  constexpr unsigned short CE_N_DIS( t_entityType i_enTy ) {
    return (i_enTy == POINT ) ? 0 :
           (i_enTy == LINE  ) ? 1 :
           (i_enTy == QUAD4R) ? 2 :
           (i_enTy == TRIA3 ) ? 2 :
           (i_enTy == HEX8R ) ? 3 :
           (i_enTy == TET4  ) ? 3 : 0;
  }

  /**
   * Gets the number of vertices for the given entity type.
   *
   * @param i_enTy entity type.
   * @return number of vertices.
   **/
  constexpr unsigned short CE_N_VES( t_entityType i_enTy ) {
    return (i_enTy == POINT ) ? 1 :
           (i_enTy == LINE  ) ? 2 :
           (i_enTy == QUAD4R) ? 4 :
           (i_enTy == TRIA3 ) ? 3 :
           (i_enTy == HEX8R ) ? 8 :
           (i_enTy == TET4  ) ? 4 : 0;
  }

  /**
   * Gets the number of faces for the given entity type.
   *
   * @param i_enTy entity type.
   * @return number of faces.
   **/
  constexpr unsigned short CE_N_FAS( t_entityType i_enTy ) {
    return (i_enTy == POINT ) ? 0 :
           (i_enTy == LINE  ) ? 2 :
           (i_enTy == QUAD4R) ? 4 :
           (i_enTy == TRIA3 ) ? 3 :
           (i_enTy == HEX8R ) ? 6 :
           (i_enTy == TET4  ) ? 4 : 0;
  }

  /**
   * Gets the face type for the given element type.
   *
   * @param i_elTy element type.
   * @return entity type of the faces.
   **/
  constexpr t_entityType CE_T_FA( t_entityType i_enTy ) {
    return (i_enTy == LINE  ) ? POINT :
           (i_enTy == QUAD4R) ? LINE :
           (i_enTy == TRIA3 ) ? LINE :
           (i_enTy == HEX8R ) ? QUAD4R :
           (i_enTy == TET4  ) ? TRIA3 : INVALID;
  }
}

#endif