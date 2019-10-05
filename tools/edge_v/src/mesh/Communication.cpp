/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2019, Alexander Breuer
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
 * Derives communication structures.
 **/
#include "Communication.h"
#include "io/logging.h"

void edge_v::mesh::Communication::getPaTgPairs( t_entityType                           i_elTy,
                                                std::size_t                            i_first,
                                                std::size_t                            i_size,
                                                std::size_t                    const * i_elFaEl,
                                                unsigned short                 const * i_elTg,
                                                std::size_t                    const * i_elPa,
                                                std::set<
                                                std::pair< std::size_t,
                                                            unsigned short > >       & o_pairs ) {
  unsigned short l_nElFas = CE_N_FAS( i_elTy );
  EDGE_V_CHECK_GT( i_size, 0 );

  // reset pairs
  o_pairs.empty();

  // get the region's partition and time group
  unsigned short l_pa = i_elPa[i_first];
  unsigned short l_tg = i_elTg[i_first];

  // iterate over the elements of the region
  for( std::size_t l_el = i_first; l_el < i_first+i_size; l_el++ ) {
    // check that the partition and time group match
    EDGE_V_CHECK_EQ( i_elPa[l_el], l_pa );
    EDGE_V_CHECK_EQ( i_elTg[l_el], l_tg );

    // iterate over the face-adjacent elements
    for( unsigned short l_fa = 0; l_fa < l_nElFas; l_fa++ ) {
      std::size_t l_ad = i_elFaEl[ l_el*l_nElFas + l_fa ];
      if( l_ad != std::numeric_limits< std::size_t >::max() ) {
        // get adjacent partition and time group
        std::size_t l_adPa = i_elPa[l_ad];
        std::size_t l_adTg = i_elTg[l_ad];

        // communication is required, if adjacent elements belongs to a different partition
        if( l_pa != l_adPa ) {
          std::pair< std::size_t, unsigned short > l_pair( {l_adPa, l_adTg} );
          o_pairs.insert( l_pair );
        }
      }
    }
  }
}