/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2019-2020, Alexander Breuer
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
 * Generic geometry computations.
 **/
#ifndef EDGE_V_GEOM_GENERIC_H
#define EDGE_V_GEOM_GENERIC_H

#include <cstddef>
#include "../constants.h"

namespace edge_v {
  namespace geom {
    class Generic;
  }
}

class edge_v::geom::Generic {
  public:
    /**
     * Computes the centroid of the given entity.
     *
     * @param i_enTy entity type.
     * @param i_veCrds vertex coordinates.
     * @param o_cen will be set to centroid of the entity.
     **/
    static void centroid( t_entityType         i_enTy,
                          double       const (*i_veCrds)[3],
                          double               o_cen[3] );

    /**
     * Sorts the given data lexicographically.
     * Total number of entries: n0 x n1.
     *
     * @param i_n0 number of entries in the slow dimension.
     * @param i_n1 number of entries in the fast dimension.
     * @param io_data data which is sorted.
     * @param o_mapping will be set to mapping from sorted to original. nullptr if unused.
     **/
    static void sortLex( t_idx   i_n0,
                         t_idx   i_n1,
                         t_idx * io_data,
                         t_idx * o_mapping = nullptr );

    /**
     * Gets the adjacent elements' face ids for the given element-faces.
     *
     * @param i_elTy element type.
     * @param i_nFas number of element-faces to get ids for.
     * @param i_elOff element offset (added to element-ids in queries).
     * @param i_el elements to which the faces belong.
     * @param i_fa local face ids w.r.t. the elements.
     * @param i_elFaEl elements adjacent to elemetns (faces as bridge).
     * @param o_faIdsAd will be set to faces ids w.r.t. to the adjacent elements.
     **/
    static void getFaIdsAd( t_entityType           i_elTy,
                            t_idx                  i_nFas,
                            t_idx                  i_elOff,
                            t_idx    const       * i_el,
                            unsigned short const * i_fa,
                            t_idx    const       * i_elFaEl,
                            unsigned short       * o_faIdsAd );
};

#endif