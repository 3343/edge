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
 * Sparse types in the mesh.
 **/
 
#ifndef SPARSE_TYPES_HPP
#define SPARSE_TYPES_HPP

#include <vector>

namespace edge {
  namespace mesh {
    template< t_entityType TL_T_EN >
    class SparseTypes;
  }
}

/**
 * Sparse types.
 *
 * @paramt TL_T_EN entity type.
 **/
template< t_entityType TL_T_EN >
class edge::mesh::SparseTypes {
  private:
    //! number of vertices
    static unsigned short const TL_N_VE = C_ENT[ TL_T_EN ].N_VERTICES;

  public:
    /**
     * Sets the given sparse types of all entities in the respective domains.
     *
     * @param i_nEnts number of entities.
     * @param i_enVe vertices adjacent to the entities.
     * @param i_spTypes sparse types in the domains.
     * @param i_doms domains where the sparse types of the entitities will be adjusted.
     * @param i_veChars vertex characteristics.
     * @param io_enChars will be updated with the given sparse types (bitwise or) if the respective entities are inside the domain.
     *
     * @paramt TL_T_INT_LID integer type of local ids.
     * @paramt TL_T_INT_SP integer type of the sparse type.
     * @paramt TL_T_DOM type of the domains, offes inside-query.
     * @paramt TL_T_VE_CHARS type of the vertex characteristics, has a coords-array for coordinates.
     * @oaramt TL_T_EN_CHARS type of the entity characteristics, has a member spType for the sparse type.
     **/
    template< typename TL_T_INT_LID,
              typename TL_T_INT_SP,
              typename TL_T_DOM,
              typename TL_T_VE_CHARS,
              typename TL_T_EN_CHARS >
    static void set( TL_T_INT_LID                     i_nEnts,
                     TL_T_INT_LID            const (* i_enVe)[TL_N_VE],
                     TL_T_INT_SP             const  * i_spTypes,
                     std::vector< TL_T_DOM > const  & i_doms,
                     TL_T_VE_CHARS           const  * i_veChars,
                     TL_T_EN_CHARS                  * io_enChars ) {
      // iterate over entities
      for( TL_T_INT_LID l_en = 0; l_en < i_nEnts; l_en++ ) {
        // iterate over domains
        for( std::size_t l_do = 0; l_do < i_doms.size(); l_do++ ) {
          // true if the entire entity is inside the domain
          bool l_in = true;

          // iterate over vertices
          for( unsigned short l_ve = 0; l_ve < TL_N_VE; l_ve++ ) {

            // get the vertex id
            TL_T_INT_LID l_veId = i_enVe[l_en][l_ve];

            if( !i_doms[l_do].inside( i_veChars[l_veId].coords ) ) {
              l_in = false;
              break;
            }
          }

          // set the sparse type accordingly if the entity is inside the domain
          if( l_in ) {
            io_enChars[l_en].spType |= i_spTypes[l_do];
          }
        }
      }
    }

};

#endif
