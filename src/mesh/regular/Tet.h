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
 * Global time stepping setup of regular tetrahedral meshes.
 **/

#ifndef TET_H_
#define TET_H_

#include "constants.hpp"
#include "data/layout.hpp"
#include <vector>
#include <limits>
#include "io/logging.h"

namespace edge {
  namespace mesh {
    namespace regular {
      class Tet;
    }
  }
}

class edge::mesh::regular::Tet {
  private:
    //! number of hexahedral elements in x-, y- and z-direction.
    unsigned int m_nHex[3];

    //! corner of lower-left-bottom, owned tet
    double m_corner[3];

    //! hex sizes
    double m_dX[3];

    //! index mappings
    t_inMap m_inMap;

    //! true if the boundaries are periodic
    bool m_periodic = false;

    //! mpi neighbors in all dimensions, left and right
    int m_mpiNe[3][2] = {{0,0},{0,0},{0,0}};

    /**
     * Maps the local vertices of a tetrahedral face to those of the containing hexahedron.
     * Remark: Order is dimension-wise, not by id. This matches our global sorting of the vertices.
     * [*][][]: Element type, [][*][]: face, [][][*]: vertex
     **/
    const unsigned short m_mapFaVe[2][16][3] {
      // type A)
      {
        {0,3,7},
        {0,4,7},
        {1,2,5},
        {2,5,6},
        {0,1,5},
        {0,4,5},
        {3,2,7},
        {2,7,6},
        {0,1,2},
        {0,3,2},
        {4,5,7},
        {5,7,6},
        {0,2,5},
        {0,2,7},
        {0,5,7},
        {2,5,7}
      },
      // type B)
      {
        {0,3,4},
        {3,4,7},
        {1,2,6},
        {1,5,6},
        {0,1,4},
        {1,4,5},
        {3,2,6},
        {3,7,6},
        {0,1,3},
        {1,3,2},
        {4,5,6},
        {4,7,6},
        {1,3,4},
        {1,3,6},
        {1,4,6},
        {3,4,6}
      }
    };

   /**
    * Maps the local vertices of a tet to the corresponding vertices of the hex.
    * Remark: The order is dimension-wise, following our global ordering.
    *
    * [*][][]: type, [][*][]: tet in the hex, [][][*]: vertex
    **/
    const unsigned short m_mapElVe[2][5][4] {
      {
        {0,1,2,5},
        {0,3,2,7},
        {0,2,5,7},
        {0,4,5,7},
        {2,5,7,6}
      },
      {
        {0,1,3,4},
        {1,3,2,6},
        {1,3,4,6},
        {1,4,5,6},
        {3,4,7,6}
      }
    };

    /**
     * Maps one tet-face in a hex to the two adjacent tets (maybe in different hexes).
     * Sorted by global element-id.
     * [*][][]: Element type, [][*][]: face, [][][*]: left/right
     **/
    const unsigned short m_mapFaEl[2][16][2] {
      // type A)
      {
        {1,1}, // 0-3-7
        {3,3}, // 0-4-7
        {0,0}, // 1-2-5
        {4,4}, // 2-5-6
        {1,0}, // 0-1-5
        {4,3}, // 0-4-5
        {1,0}, // 2-3-7
        {4,3}, // 2-6-7
        {3,0}, // 0-1-2
        {4,1}, // 0-2-3
        {3,0}, // 4-5-7
        {4,1}, // 5-6-7
        {0,2}, // 0-2-5
        {1,2}, // 0-2-7
        {2,3}, // 0-5-7
        {2,4}, // 2-5-7
      },
      // type B)
      {
        {0,0}, // 0-3-4
        {4,4}, // 3-4-7
        {1,1}, // 1-2-6
        {3,3}, // 1-5-6
        {1,0}, // 0-1-4
        {4,3}, // 1-4-5
        {1,0}, // 2-3-6
        {4,3}, // 3-6-7
        {3,0}, // 0-1-3
        {4,1}, // 1-2-3
        {3,0}, // 4-5-6
        {4,1}, // 4-6-7
        {0,2}, // 1-3-4
        {1,2}, // 1-3-6
        {2,3}, // 1-4-6
        {2,4}, // 3-4-6
      }
    };

    /**
     * Maps a given tetrahedral face in a hexahedral element to the offset w.r.t. to the hexe's position.
     * [*][][]: face, [][*][]: left/right, [][*]: dimension
     **/
    const short m_mapFaHxOff[16][2][3] {
      { {-1, 0, 0}, { 0, 0, 0} }, // 0
      { {-1, 0, 0}, { 0, 0, 0} }, // 1

      { { 0, 0, 0}, { 1, 0, 0} }, // 2
      { { 0, 0, 0}, { 1, 0, 0} }, // 3

      { { 0,-1, 0}, { 0, 0, 0} }, // 4
      { { 0,-1, 0}, { 0, 0, 0} }, // 5

      { { 0, 0, 0}, { 0, 1, 0} }, // 6
      { { 0, 0, 0}, { 0, 1, 0} }, // 7

      { { 0, 0,-1}, { 0, 0, 0} }, // 8
      { { 0, 0,-1}, { 0, 0, 0} }, // 9

      { { 0, 0, 0}, { 0, 0, 1} }, // 10
      { { 0, 0, 0}, { 0, 0, 1} }, // 11

      { { 0, 0, 0}, { 0, 0, 0} }, // 12
      { { 0, 0, 0}, { 0, 0, 0} }, // 13
      { { 0, 0, 0}, { 0, 0, 0} }, // 14
      { { 0, 0, 0}, { 0, 0, 0} }  // 15
    };

    /**
     * Maps one tet to the adjacent faces.
     * Sorted by global face-id.
     * [*][][]: Element type, [][*][]: element, [][][*]: tet-face
     **/
   const unsigned short m_mapElFa[2][5][4] {
     // type A)
     {
       {  8,  4, 12,  2}, // 0-1-2-5: 0-1-2, 0-1-5, 0-2-5, 1-2-5
       {  9,  0, 13,  6}, // 0-2-3-7: 0-2-3, 0-3-7, 0-2-7, 2-3-7
       { 12, 13, 14, 15}, // 0-2-5-7: 0-2-5, 0-2-7, 0-5-7, 2-5-7
       {  5,  1, 14, 10}, // 0-4-5-7: 0-4-5, 0-4-7, 0-5-7, 4-5-7
       { 15,  3,  7, 11}  // 2-5-6-7: 2-5-7, 2-5-6, 2-6-7, 5-6-7
     },
     // type B)
     {
       {  8,  4,  0, 12}, // 0-1-3-4: 0-1-3, 0-1-4, 0-3-4, 1-3-4
       {  9, 13,  2,  6}, // 1-2-3-6: 1-2-3, 1-3-6, 1-2-6, 2-3-6
       { 12, 13, 14, 15}, // 1-3-4-6: 1-3-4, 1-3-6, 1-4-6, 3-4-6
       {  5, 14,  3, 10}, // 1-4-5-6: 1-4-5, 1-4-6, 1-5-6, 4-5-6
       {  1, 15,  7, 11}  // 3-4-6-7: 3-4-7, 3-4-6, 3-6-7, 4-6-7
     }
   };

    /**
     * Maps one tet to the adjacent tets.
     * Sorted by global face-id.
     * [*][][][]: Element type, [][*][][]: element, [][][*][]: face
     * [][][][0]: tet id, [][][][1]: x-offset, [][][][2]: y-offset, [][][][3]: z-offset 
     **/
   const short m_mapElFaEl[2][5][4][4] {
     // type A)
     {
       {
         {3, 0, 0,-1}, // 0-1-2
         {1, 0,-1, 0}, // 0-1-5
         {2, 0, 0, 0}, // 0-2-5
         {0, 1, 0, 0}, // 1-2-5
       },
       {
         {4, 0, 0,-1}, // 0-2-3
         {1,-1, 0, 0}, // 0-3-7
         {2, 0, 0, 0}, // 0-2-7
         {0, 0, 1, 0}, // 2-3-7
       },
       {
         {0, 0, 0, 0}, // 0-2-5
         {1, 0, 0, 0}, // 0-2-7
         {3, 0, 0, 0}, // 0-5-7
         {4, 0, 0, 0}, // 2-5-7
       },
       {
         {4, 0,-1, 0}, // 0-4-5
         {3,-1, 0, 0}, // 0-4-7
         {2, 0, 0, 0}, // 0-5-7
         {0, 0, 0, 1}, // 4-5-7
       },
       {
         {2, 0, 0, 0}, // 2-5-7
         {4, 1, 0, 0}, // 2-5-6
         {3, 0, 1, 0}, // 2-6-7
         {1, 0, 0, 1}, // 5-6-7
       },
     },
     // type B)
     {

       {
         {3, 0, 0,-1}, // 0-1-3
         {1, 0,-1, 0}, // 0-1-4
         {0,-1, 0, 0}, // 0-3-4
         {2, 0, 0, 0}, // 1-3-4
       },
       {
         {4, 0, 0,-1}, // 1-2-3
         {2, 0, 0, 0}, // 1-3-6
         {1, 1, 0, 0}, // 1-2-6
         {0, 0, 1, 0}, // 2-3-6
       },
       {
         {0, 0, 0, 0}, // 1-3-4
         {1, 0, 0, 0}, // 1-3-6
         {3, 0, 0, 0}, // 1-4-6
         {4, 0, 0, 0}, // 3-4-6
       },
       {
         {4, 0,-1, 0}, // 1-4-5
         {2, 0, 0, 0}, // 1-4-6
         {3, 1, 0, 0}, // 1-5-6
         {0, 0, 0, 1}, // 4-5-6
       },
       {
         {4,-1, 0, 0}, // 3-4-7
         {2, 0, 0, 0}, // 3-4-6
         {3, 0, 1, 0}, // 3-6-7
         {1, 0, 0, 1}, // 4-6-7
       }
     },
   };

    //! divided hex
    typedef struct {
      bool   type;
      int_el ve[ 8];
      int_el fa[16];
      std::vector< int_el > el[5];
    } t_hexDiv;

    //! hexes
    std::vector< t_hexDiv > m_hexes;

    //! number of elements in the mesh (non-redundant data)
    int_el m_nMeshEl = 0;

    //! entity layouts
    t_enLayout m_veLayout;
    t_enLayout m_faLayout;
    t_enLayout m_elLayout;

    //! entity in a hex
    typedef struct {
      bool           type;
      unsigned short en;
      int            xHex;
      int            yHex;
      int            zHex;
    } t_enInHex;

    //! inner elements
    std::vector< t_enInHex > m_elInner;

    //! send elements
    std::vector< t_enInHex > m_elSend[6];

    //! receive elements
    std::vector< t_enInHex > m_elRecv[6];

    //! global ids
    std::vector< int_gid > m_gIdsVe;
    std::vector< int_gid > m_gIdsFa;
    std::vector< int_gid > m_gIdsEl;

    /**
     * Determines if a given neighboring rank qualifies for an MPI-boundary.
     *
     * @param i_rank rank in question.
     * @return true if the rank qualifies, false otherwise.
     **/
     bool isMpiBnd( int i_rank ) const {
       return i_rank >= 0 && i_rank != parallel::g_rank && i_rank < parallel::g_nRanks;
     }

    /**
     * Splitting A):
     *
     * 0: 0-1-2-5
     * 1: 0-2-3-7
     * 2: 0-2-5-7
     * 3: 0-4-5-7
     * 4: 2-5-6-7
     *
     *          7                                  6
     *           x*******************************x
     *           * *                             *
     *          *.   *.                         **
     *         * .      *.                     * *
     *        *. .         *                  *  *
     *       *   .           .*              *   *
     *      *    .              .*          *    *
     *     *  .  .                . *      *     *
     *  4 *      .                   . *  *      *
     *   x*******************************x  5    *
     *   *       .                     * *       *
     *   *   .   .                       *  .    *
     *   *       . 3            *        *    .  *
     *   *       x.......................*.......x 2
     *   *  .   .                        *  .   *
     *   *     .       *             .   *     *
     *   *    .                .         *    *
     *   * . .    *       .              *   *
     *   *  .        .                   *  *
     *   * . *  .                        * *
     *   *.                              **
     *   x*******************************x
     *  0                                 1
     *
     *
     * Splitting B):
     *
     * 0: 0-1-3-4
     * 1: 1-2-3-6
     * 2: 1-3-4-6
     * 3: 1-4-5-6
     * 4: 3-4-6-7
     *
     *          7                                  6
     *           x*******************************x
     *           *                           *   *
     *          *.                     *    .   **
     *         * .                 *           * *
     *        *  .            *       .       *  *
     *       *   .        *                  * . *
     *      *    .   *                      *    *
     *     *    *.             .           *     *
     *  4 *      .                        *   .  *
     *   x*******************************x  5    *
     *   *  .    .      .                *       *
     *   *  *    .   .                   *   .   *
     *   *     . . 3                     *       *
     *   *       x.......................*.......x 2
     *   *      .   *.                   *  .   *
     *   *     .        .                *     *
     *   *    .           *              *    *
     *   *   .                .          * . *
     *   *  .                   *        *  *
     *   * .                        .    * *
     *   *.                           *  **
     *   x*******************************x
     *  0                                 1
     **/

    /**
     * Gets the number of owned vertices (inner+send).
     *
     * @return number of owned vertices.
     **/
    int_el getNVeOwned() const;

    /**
     * Gets the number of inner-vertices.
     *
     * @return number of inner-vertices.
     **/
    int_el getNVeInner() const;

    /**
     * Gets the number of send-vertices
     *
     * @return number of send-vertices.
     **/
    int_el getNVeSend() const;

    /**
     * Gets the number of receive-vertices.
     *
     * @return number of receive-vertices.
     **/
    int_el getNVeRecv() const;

    /**
     * Gets the number of vertices (inner+send+receive)
     *
     * @return number of vertices.
     **/
    int_el getNVe() const;

    /**
     * Gets the number of send-faces.
     * Remark: Only the right, back, top face is considered to be send-faces.
     *
     * @return number of send-faces.
     **/
    int_el getNFaSend() const;

    /**
     * Gets the number of receive-faces.
     * Remark: Only the left, front, bottom face is considered to be receive-faces.
     *
     * @return number of receive-faces.
     **/
    int_el getNFaRecv() const;

    /**
     * Gets the number of inner-faces.
     *
     * @return number of inner-faces.
     **/
    int_el getNFaInner() const;


    /**
     * Gets the number of faces (inner+send+receive).
     **/
    int_el getNFa() const;

    /**
     * Gets the number of send-elements.
     *
     * @return number of send-elements.
     **/
    int_el getNElSend() const;

    /**
     * Gets the number of receive-elements.
     *
     * @return number of receive-elements.
     **/
    int_el getNElRecv() const;

    /**
     * Gets the number of owned elements (inner+send).
     **/
    int_el getNElOwned() const;

    /**
     * Gets the number of inner elements.
     **/
    int_el getNElInner() const;

    /**
     * Gets the number of elements (inner+send+receive).
     **/
    int_el getNEl() const;

    /**
     * Tests if a given tet in a hex is a receive element.
     *
     * @param i_tet tet id within the hex.
     * @parma i_hxPos x-, y- and z-position of the containing hex (-1 or #hexes if non-owned).
     * @return true if the element is receive, false otherwise.
     **/
    bool isRecvEl(       unsigned short i_tet,
                   const int            i_hxPos[3] ) const;

    /**
     * Tests if a given tet in a hex is a send element.
     *
     * @param i_tet tet id within the hex.
     * @parma i_hxPos x-, y- and z-position of the containing hex (-1 or #hexes if non-owned).
     * @param o_send true if send-element in x-direction ([0]), y-direction ([1]), z-direction([2]); false otherwise.
     **/
    void isSendEl(       unsigned short i_tet,
                   const int            i_hxPos[3],
                         bool           o_send[3] ) const;

    /**
     * Test if the given element is an inner element.
     *
     * @param i_tet tet id within the hex.
     * @param i_hxPos x-, y- and z-position of the containing hex (-1 or #hexes if non-owned).
     *
     * @return true if the element is an inner element, false otherwise.
     **/
    bool isInnerEl(       unsigned short i_tet,
                    const int            i_hxPos[3] ) const;

    /**
     * Gets the index of an hex specified through x-, y-, z-position.
     *
     * @param i_xh x-pos.
     * @param i_yh y-pos.
     * @param i_zh z-pos.
     * @param i_nXHex number of owned hexes in x-direction.
     * @param i_nYHex number of owned hexes in y-direction.
     * @parma i_nZHex number of owned hexes in z-direction.
     * @return index.
     **/
    static unsigned int getHxId(          int i_xh,             int i_yh,             int i_zh,
                                 unsigned int i_nXHex, unsigned int i_nYHex, unsigned int i_nZHex ) {
      // compute resulting hex
      return (i_zh+1)*(i_nXHex+2)*(i_nYHex+2) +
             (i_yh+1)*(i_nXHex+2) +
             (i_xh+1);
    }

    /**
     * Gets the index of an hex, but respect periodic boundaries if present.
     *
     * @param i_xh x-pos.
     * @param i_yh y-pos.
     * @param i_zh z-pos.
     * @param i_nXHex number of owned hexes in x-direction.
     * @param i_nYHex number of owned hexes in y-direction.
     * @parma i_nZHex number of owned hexes in z-direction.
     * @return index, numeric_limits<unsigned int>::max() if the hex logically does not exist.
     **/
    unsigned int getHxIdBnd(              int i_xh,             int i_yh,             int i_zh,
                                 unsigned int i_nXHex, unsigned int i_nYHex, unsigned int i_nZHex ) const {
      // compute resulting hex
      // derive normalized hex coords
      int l_xn, l_yn, l_zn;
      l_xn = l_yn = l_zn = std::numeric_limits< int_el >::max();

      if( m_periodic ) {
        if(      i_xh == -1              && m_mpiNe[0][0] == parallel::g_rank ) l_xn = m_nHex[0]-1;
        else if( i_xh == (int) m_nHex[0] && m_mpiNe[0][1] == parallel::g_rank ) l_xn = 0;
        else                                                                    l_xn = i_xh;

        if(      i_yh == -1              && m_mpiNe[1][0] == parallel::g_rank ) l_yn = m_nHex[1]-1;
        else if( i_yh == (int) m_nHex[1] && m_mpiNe[1][1] == parallel::g_rank ) l_yn = 0;
        else                                                                    l_yn = i_yh;

        if(      i_zh == -1              && m_mpiNe[2][0] == parallel::g_rank ) l_zn = m_nHex[2]-1;
        else if( i_zh == (int) m_nHex[2] && m_mpiNe[2][1] == parallel::g_rank ) l_zn = 0;
        else                                                                    l_zn = i_zh;
      }
      else {
        if( i_xh == -1              || i_yh == -1              || i_zh == -1 ||
            i_xh == (int) m_nHex[0] || i_yh == (int) m_nHex[1] || i_zh == (int) m_nHex[2] ) {
          return std::numeric_limits<unsigned int>::max();
        }
        else {
          l_xn = i_xh;
          l_yn = i_yh;
          l_zn = i_zh;
        }
      }

      return getHxId( l_xn, l_yn, l_zn, i_nXHex, i_nYHex, i_nZHex );
    }

    /**
     * Gets the vertices adjacent to the faces.
     *
     * @param o_faVe will be set to vertices adjacent to the faces
     **/
    void getFaVe( int_el (*o_faVe)[3] );

    /**
     * Gets the vertices adjacent to the elements.
     *
     * @param o_elVe will be set to vertices adjacent to the elements.
     **/
    void getElVe( int_el (*o_elVe)[4] );

    /**
     * Gets the elements adjacent to the faces.
     *
     * @param o_faVe will be set to vertices adjacent to the faces
     **/
    void getFaEl( int_el (*o_faEl)[2] );

    /**
     * Gets the faces adjacent to the elements.
     *
     * @param o_elFa will be set to faces adjacent to elements.
     **/
    void getElFa( int_el (*o_elFa)[4] );

    /**
     * Gets the elements adjacent to elements through faces.
     *
     * @param o_elFaEl will be set to elements adjacent to elements through faces.
     **/
    void getElFaEl( int_el (*o_elFaEl)[4] );

    /**
     * Derives the index mappings for vertices and faces.
     *
     * TODO: In the case of vertices potential duplicates for MPI are ignored.
     **/
    void deriveInMapVeAndFa();

    /**
     * Initializes the send faces for the given base hex.
     *
     * @param i_dim mpi-dimension which is initialized.
     * @param i_hxPos pos of the hex in x-, y- and z- direction.
     * @param i_cIds communication ids.
     **/
    void initSendFa(       unsigned short i_dim,
                     const int            i_hxPos[3],
                     const unsigned short i_cIds[6] );

    /**
     * Initializes the receive faces for the given base hex.
     *
     * @param i_dim mpi-dimension which is initialized.
     * @param i_id pos of the hex in x-, y- and z- direction.
     * @param i_cIds communication ids.
     **/
    void initRecvFa(       unsigned short i_dim,
                     const int            i_id[3],
                     const unsigned short i_cIds[6] );

    /**
     * Initializes the send elements for the given base hex.
     *
     * @param i_dim mpi-dimension which is initialized.
     * @param i_id pos of the hex in x-, y- and z- direction.
     * @param i_cIds communication ids.
     **/
    void initSendEl(       unsigned short i_dim,
                     const int            i_id[3],
                     const unsigned short i_cIds[6] );

    /**
     * Initializes the receive elements for the given base hex.
     *
     * @param i_dim mpi-dimension which is initialized.
     * @param i_id pos of the hex in x-, y- and z- direction.
     * @param i_cIds communication ids.
     **/
    void initRecvEl(       unsigned short i_dim,
                     const int            i_id[3],
                     const unsigned short i_cIds[6] );

    /**
     * Initializes the inner-, send- and recv-first position in the entity layout based on the sizes.
     *
     * @param io_enLo entity layout which first positions are set.
     **/
    void initFirstPos( t_enLayout &io_enLo );

    /**
     * Initializes the regular tetrahedral mesh (private).
     *
     * @param i_type if false, splitting A) is assumed for the lower-bottom-left owned hex-base-element, else B).
     * @param i_nHex number of owned, base hexes in x-, y- and z-direction.
     * @param i_mpiNe mpi-rank of left and right x-, y- and z-bnd. <0 if bnd cond.
     * @param i_corner x-, y- and z-coords of bottom-lower-left vertex (adjacent to an owned element).
     * @param i_dX width dx, dy and dz of the hexes.
     **/
    void init(       bool         i_type,
               const unsigned int i_nHex[3],
               const int          i_mpiNe[3][2],
               const double       i_corner[3],
               const double       i_dX[3] );

  public:
    /**
     * Initializes the regualar tetrahedral mesh.
     *
     * @param i_nHex[3] number of owned, base hexes in x-, y- and z-direction.
     * @param i_rank calling rank.
     * @param i_nRanks number of ranks.
     * @param i_corner x-, y- and z-coords of bottom-lower-left vertex (adjacent to an owned element).
     * @param i_dX width dx, dy and dz of the hexes.
     * @param i_periodic true if boundaries are periodic.
     **/
    void init( const unsigned int i_nHex[3],
                     int          i_rank,
                     int          i_nRanks,
               const double       i_corner[3],
               const double       i_dX[3],
                     bool         i_periodic=true );

    /**
     * Gets the vertex-layout.
     *
     * Remark: Matching indices of faVe and elVe connect info, the vertex layout
     *         contains all vertices of the hexahedral base elements.
     *
     * TODO: Only the #vertices and #owned vertices is returned.
     *       No seperation among MPI ranks.
     *
     * @return vertex-layout.
     **/
    t_enLayout getVeLayout() { return m_veLayout; };

    /**
     * Gets the face-layout.
     *
     * @return face-layout.
     **/
    t_enLayout getFaLayout() { return m_faLayout; };

    /**
     * Gets the element-layout.
     *
     * @return element-layout.
     **/
    t_enLayout getElLayout() { return m_elLayout; };

    /**
     * Gets the vertex characteristics for all vertices (including those of ghost elements).
     *
     * Remark: Matching indices of faVe and elVe connect info, the vertex characteristics
     *         are set for all vertices of the hexahedral base elements.
     *
     * @param o_veChars will be set to vertex characteristics (should be large enough to cover all ves of hex-base).
     **/
    void getVeChars( t_vertexChars *o_veChars );

    /**
     * Gets the face characteristics for all faces.
     *
     * @param o_faChars will be set to face characteristics.
     **/
    void getFaChars( t_faceChars *o_faChars );

    /**
     * Gets the element characteristics for all elements.
     *
     * @param o_elChars will be set to element characteristics.
     **/
    void getElChars( t_elementChars *o_elChars );

    /**
     * Gets the global ids of the vertices.
     *
     * TODO: This is a dummy implementation!
     **/
    void getGIdsVe( std::vector< int_gid > &o_gIds );

    /**
     * Gets the global ids of the faces.
     *
     * TODO: This is a dummy implementation!
     **/
    void getGIdsFa( std::vector< int_gid > &o_gIds );

    /**
     * Gets the global ids of the elements.
     *
     * TODO: This is a dummy implementation!
     **/
    void getGIdsEl( std::vector< int_gid > &o_gIds );

    /**
     * Gets the index mapping for mesh-to-data and data-to-mesh.
     *
     * @return will be set to index mapping.
     **/
     const t_inMap* getInMap() const;

    /**
      * Gets the connectivity information for all elements.
      *
      * @param i_veChars vertex chars, will be used to determine the mapping to the reference element.
      * @param i_faChars face chars, will be used to identify boundary faces in consistency checks.
      * @param o_connect will be set to connectity information.
      **/
    void getConnect( const t_vertexChars  *i_veChars,
                     const t_faceChars    *i_faChars,
                           t_connect      &o_connect );
};

#endif
