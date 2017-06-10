/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2016, Regents of the University of California
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
 * Base discretization into quads or hexes.
 **/

#ifndef BASE_H_
#define BASE_H_

#include "constants.hpp"
#include "io/logging.h"

namespace edge {
  namespace mesh {
    namespace regular {
      class Base;
    }
  }
}

class edge::mesh::regular::Base {
  public:
    /**
     * Gets the base mesh setup in 2D.
     *
     * @param i_nRanks number of ranks.
     * @param i_rank rank for which the base mesh is derived.
     * @param i_periodic if true, periodic boundaries are initialized.
     * @param i_nQuadAll number of quads for all partitions.
     * @param o_nQuad number of base elements in x- and y-direction (local rank).
     * @param o_nPart number of partitions in x- and y-direction.
     * @param o_part partition of local rank in x- and y-direction.
     **/
    static void getSetup2d(       int          i_nRanks,
                                  int          i_rank,
                                  bool         i_periodic,
                            const unsigned int i_nQuadAll[2],
                                  unsigned int o_nQuad[2],
                                  unsigned int o_nPart[2],
                                  unsigned int o_part[2] ) {
      EDGE_CHECK( i_periodic ); // TODO: Add support for non-periodic boundaries

      EDGE_CHECK( i_rank < i_nRanks );

      // derive the per-dimension paritioning
      o_nPart[0] = o_nPart[1] = 1;

      // find biggest (balanced) divisor
      int l_di = 1;
      for(  l_di = 1; l_di < i_nRanks; l_di++ ) {
        if( i_nRanks % l_di == 0 && l_di >= (i_nRanks / l_di) ) break;
      }

      if( i_nQuadAll[0] > i_nQuadAll[1] ) {
        o_nPart[0] = l_di;
        o_nPart[1] = i_nRanks / l_di;
      }
      else {
        o_nPart[0] = i_nRanks / l_di;
        o_nPart[1] = l_di;
      }

      // check that the #parts match
      EDGE_CHECK( o_nPart[0] * o_nPart[1] == (unsigned int) i_nRanks );
      // check that there's elements left in the partitions
      EDGE_CHECK( i_nQuadAll[0] / o_nPart[0] > 1 );
      EDGE_CHECK( i_nQuadAll[1] / o_nPart[1] > 1 );

      // derive local "coords" of the partition
      o_part[0] = i_rank % o_nPart[0];
      o_part[1] = i_rank / o_nPart[0];

      // derive local sizes
      o_nQuad[0] = i_nQuadAll[0] / o_nPart[0];
      o_nQuad[1] = i_nQuadAll[1] / o_nPart[1];

      // distribute remaining elements equally
      if( i_nQuadAll[0]%o_nPart[0] > o_part[0] ) o_nQuad[0]++;
      if( i_nQuadAll[1]%o_nPart[1] > o_part[1] ) o_nQuad[1]++;
    }

    /**
     * Gets the rank associated with a partition.
     *
     * @param i_nPart number of partitions in x-, y- and z-direction.
     * @param i_partX partition of local rank in x-, y- and z-direction.
     **/
    static int getRank( unsigned int i_nPart[3],
                        unsigned int i_partX[3] ) {
      return i_partX[2] * (i_nPart[0]*i_nPart[1]) + i_partX[1] * i_nPart[0] + i_partX[0];
    }

    /**
     * Gets the base mesh setup in 3D.
     *
     * The implementation tries to find a homogeneous number of partitions in all dimensions.
     * If this is not possible the dimension with the largest number of elements gets the
     * highest number of partitions.
     *
     * Example: 2048 partitions, gives a 8*8*32 splitting, the dimension with the most base
     *          hexes gets the 32 partitions; no attempts are made to reach a more balanced
     *          splitting, e.g. 8*16*16.
     *
     * @param i_nRanks number of ranks.
     * @param i_rank rank for which the base mesh is derived.
     * @param i_periodic if true, periodic boundaries are initialized.
     * @param i_nHexAll number of elements in x-, y- and z-direction.
     * @param o_nHex number of base elements in x-, y- and z-direction (local rank).
     * @param o_nPart number of partitions in x-, y- and z-direction.
     * @param o_part partition of local rank in x-, y- and z-direction.
     **/
    static void getSetup3d(       int          i_nRanks,
                                  int          i_rank,
                                  bool         i_periodic,
                            const unsigned int i_nHexAll[3],
                                  unsigned int o_nHex[3],
                                  unsigned int o_nPart[3],
                                  unsigned int o_part[3] ) {
      EDGE_CHECK( i_periodic ); // TODO: Add support for non-periodic boundaries

      EDGE_CHECK( i_rank < i_nRanks );

      // derive the per-dimension paritioning
      o_nPart[0] = o_nPart[1] = o_nPart[2] = 1;

      // find biggest (balanced) divisor
      unsigned int l_di = 1;
      for(  int l_it = 1; l_it < i_nRanks; l_it++ ) {
        if( i_nRanks % (l_it*l_it) == 0 && l_it <= (i_nRanks / (l_it*l_it)) ) l_di = l_it;
      }

      // distribute the two dimensions with divisor partitioning and the remainder
      unsigned int l_rem = i_nRanks / (l_di*l_di);
      bool l_remMax = l_rem > l_di;
      bool l_remSet = false;

      if( i_nHexAll[0] > i_nHexAll[1] && i_nHexAll[0] > i_nHexAll[2] ) {
        o_nPart[0] = std::max( l_rem, l_di );
        if( l_remMax ) l_remSet = true;
      }
      else {
        o_nPart[0] = std::min( l_rem, l_di );
        if( !l_remMax ) l_remSet = true;
      }

      if( l_remSet ) o_nPart[1] = l_di;
      else if( i_nHexAll[1] > i_nHexAll[0] && i_nHexAll[1] > i_nHexAll[2] ) {
        o_nPart[1] = std::max( l_rem, l_di );
        if( l_remMax ) l_remSet = true;
      }
      else {
        o_nPart[1] = std::min( l_rem, l_di );
        if( !l_remMax ) l_remSet = true;
      }

      if( l_remSet ) o_nPart[2] = l_di;
      else           o_nPart[2] = l_rem;

      // check that the #parts match
      EDGE_CHECK( o_nPart[0] * o_nPart[1] * o_nPart[2] == (unsigned int) i_nRanks );
      // check that there's elements left in the partitions
      EDGE_CHECK( i_nHexAll[0] / o_nPart[0] > 1 );
      EDGE_CHECK( i_nHexAll[1] / o_nPart[1] > 1 );
      EDGE_CHECK( i_nHexAll[2] / o_nPart[2] > 1 );

      // derive local "coords" of the partition
      o_part[0] =   i_rank %  o_nPart[0];
      o_part[1] = ( i_rank % (o_nPart[1] * o_nPart[0]) ) / o_nPart[0];
      o_part[2] =   i_rank / (o_nPart[1] * o_nPart[0]);

      // derive local sizes
      o_nHex[0] = i_nHexAll[0] / o_nPart[0];
      o_nHex[1] = i_nHexAll[1] / o_nPart[1];
      o_nHex[2] = i_nHexAll[2] / o_nPart[2];

      // distribute remaining elements equally
      if( i_nHexAll[0]%o_nPart[0] > o_part[0] ) o_nHex[0]++;
      if( i_nHexAll[1]%o_nPart[1] > o_part[1] ) o_nHex[1]++;
      if( i_nHexAll[2]%o_nPart[2] > o_part[2] ) o_nHex[2]++;
    }

    /**
     * Gets the number of vertices in the domain.
     * Remark: Includes vertices of corners.
     *
     * @param i_nHex number of owned, base hexes in x-, y- and z-direction.
     **/
    static int_el getNVeHex( unsigned int i_nHex[3] ) {
       return (i_nHex[0]+3)*(i_nHex[1]+3)*(i_nHex[2]+3);
    }

    /**
     * Gets the vertices associated with a hex.
     *
     * TODO: This scheme breaks with the inner-send-recv ordering; adjust if necessary.
     *
     * @param i_xh x-position of the containing hex (-1 or #hexes if non-owned).
     * @param i_yh x-position of the containing hex (-1 or #hexes if non-owned).
     * @param i_zh x-position of the containing hex (-1 or #hexes if non-owned).
     * @param i_nXHex number of owned, base hexes in x-direction.
     * @param i_nYHex number of owned, base hexes in y-direction.
     * @param i_nZHex number of owned, base hexes in z-direction.
     * @param o_ves will be set to vertices
     **/
     static void getVesHex( int           i_xh,
                            int           i_yh,
                            int           i_zh,
                            unsigned int  i_nXHex,
                            unsigned int  i_nYHex,
                            unsigned int  i_nZHex,
                            int_el        o_ves[8] ) {
       o_ves[0] = (i_zh+1)*(i_nXHex+3)*(i_nYHex+3) +
                  (i_yh+1)*(i_nXHex+3) +
                  (i_xh+1);
       o_ves[1] = o_ves[0]+1;

       o_ves[3] = o_ves[0] + (i_nXHex+3);
       o_ves[2] = o_ves[3]+1;

       o_ves[4] = o_ves[0] + (i_nXHex+3)*(i_nYHex+3);
       o_ves[5] = o_ves[4] + 1;

       o_ves[7] = o_ves[4] + (i_nXHex+3);
       o_ves[6] = o_ves[7] + 1;
     }
};

#endif
