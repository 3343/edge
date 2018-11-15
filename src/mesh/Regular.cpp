/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2015-2018, Regents of the University of California
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
 * Regular meshes in EDGE.
 **/

#include "Regular.h"
#include "common.hpp"
#include <cassert>
#include <io/logging.h>

edge::mesh::Regular::Regular( Type   i_elementType,
                              int_el i_nX,
                              double i_sizeX,
                              int_el i_nY,
                              double i_sizeY,
                              int_el i_nZ,
                              double i_sizeZ ):
 m_dX(0),
 m_dY(0),
 m_dZ(0) {
  // check for valid requests
  EDGE_CHECK( i_sizeX > 0 );
  EDGE_CHECK( i_sizeY > 0 || i_elementType == Line );
  EDGE_CHECK( i_sizeZ > 0 || i_elementType != Hexahedral );
  EDGE_CHECK( i_nX > 1                             || i_elementType != Line          );
  EDGE_CHECK( ( i_nX > 1 && i_nY > 1 )             || i_elementType != Quadrilateral );
  EDGE_CHECK( ( i_nX > 2 && i_nY > 2 && i_nZ > 2 ) || i_elementType != Hexahedral    );

  // store the data
  m_elementType = i_elementType;
  m_nX    = i_nX;    m_nY = i_nY;       m_nZ = i_nZ;
  m_sizeX = i_sizeX; m_sizeY = i_sizeY; m_sizeZ = i_sizeZ;

  // derive mesh widths
  m_dX = m_sizeX / m_nX;

  if( i_elementType == Quadrilateral ||
      i_elementType == Hexahedral ) {
    m_dY = m_sizeY / m_nY;
  }

  if( i_elementType == Hexahedral ) {
    m_dZ = m_sizeZ / m_nZ;
  }

  // derive number of base elements
  if( i_elementType == Line ) {
    m_nBaseElements = m_nX;
    m_nBaseFaces    = m_nBaseElements;     // periodic boundaries
    m_nBaseVertices = m_nX+1;
  }
  else if( i_elementType == Quadrilateral ) {
    m_nBaseElements = m_nX * m_nY;
    m_nBaseFaces    = m_nBaseElements * 2; // periodic boundaries
    m_nBaseVertices = (m_nX+1)*(m_nY+1);
  }
  else if( i_elementType == Hexahedral ) {
    m_nBaseElements = m_nX * m_nY * m_nZ;
    m_nBaseFaces    = m_nBaseElements * 3; // periodic boundaries
    m_nBaseVertices = (m_nX+1)*(m_nY+1)*(m_nZ+1);
  }

  // derive number of requested elements
  if( i_elementType == Line ||
      i_elementType == Quadrilateral ||
      i_elementType == Hexahedral ) {
    m_nRequestedElements = m_nBaseElements;
    m_nRequestedFaces    = m_nBaseFaces;
    m_nRequestedVertices = m_nBaseVertices;
  }

  // setup dummy mapping
#ifdef PP_USE_MPI
  EDGE_LOG_FATAL << "fix mapping for mpi";
#endif

  m_inMap.elMeDa.resize( m_nRequestedElements ); m_inMap.elDaMe.resize( m_nRequestedElements );
  m_inMap.faMeDa.resize( m_nRequestedFaces    ); m_inMap.faDaMe.resize( m_nRequestedFaces    );
  m_inMap.veMeDa.resize( m_nRequestedVertices ); m_inMap.veDaMe.resize( m_nRequestedVertices );
  for( int_el l_el = 0; l_el < m_nRequestedElements; l_el++ ) m_inMap.elMeDa[l_el] = m_inMap.elDaMe[l_el] = l_el;
  for( int_el l_fa = 0; l_fa < m_nRequestedFaces;    l_fa++ ) m_inMap.faMeDa[l_fa] = m_inMap.faDaMe[l_fa] = l_fa;
  for( int_el l_ve = 0; l_ve < m_nRequestedVertices; l_ve++ ) m_inMap.veMeDa[l_ve] = m_inMap.veDaMe[l_ve] = l_ve;

  m_gIdsVe.resize( m_nRequestedVertices );
  m_gIdsFa.resize( m_nRequestedFaces    );
  m_gIdsEl.resize( m_nRequestedElements );
  for( int_el l_ve = 0; l_ve < m_nRequestedVertices; l_ve++ ) m_gIdsVe[l_ve] = l_ve;
  for( int_el l_fa = 0; l_fa < m_nRequestedFaces;    l_fa++ ) m_gIdsFa[l_fa] = l_fa;
  for( int_el l_el = 0; l_el < m_nRequestedElements; l_el++ ) m_gIdsEl[l_el] = l_el;
}


t_enLayout edge::mesh::Regular::getElLayout() const {
#ifdef PP_USE_MPI
  EDGE_LOG_FATAL << "mpi-implementation missing";
#endif

  t_enLayout l_elLayout;
  l_elLayout.nEnts = getNElements();
  l_elLayout.timeGroups.resize( 1 );
  l_elLayout.timeGroups[0].nEntsOwn = l_elLayout.nEnts;
  l_elLayout.timeGroups[0].nEntsNotOwn = 0;
  l_elLayout.timeGroups[0].inner.first = 0;
  l_elLayout.timeGroups[0].inner.size = l_elLayout.nEnts;

  return l_elLayout;
}

t_enLayout edge::mesh::Regular::getFaLayout() const {
#ifdef PP_USE_MPI
  EDGE_LOG_FATAL << "mpi-implementation missing";
#endif

  t_enLayout l_faLayout;
  l_faLayout.nEnts = getNFaces();
  l_faLayout.timeGroups.resize( 1 );
  l_faLayout.timeGroups[0].nEntsOwn = l_faLayout.nEnts;
  l_faLayout.timeGroups[0].nEntsNotOwn = 0;
  l_faLayout.timeGroups[0].inner.first = 0;
  l_faLayout.timeGroups[0].inner.size = l_faLayout.nEnts;

  return l_faLayout;
}

t_enLayout edge::mesh::Regular::getVeLayout() const {
#ifdef PP_USE_MPI
  EDGE_LOG_FATAL << "mpi-implementation missing";
#endif

  t_enLayout l_veLayout;
  l_veLayout.nEnts = getNVertices();
  l_veLayout.timeGroups.resize( 1 );
  l_veLayout.timeGroups[0].nEntsOwn = l_veLayout.nEnts;
  l_veLayout.timeGroups[0].nEntsNotOwn = 0;
  l_veLayout.timeGroups[0].inner.first = 0;
  l_veLayout.timeGroups[0].inner.size = l_veLayout.nEnts;

  return l_veLayout;
}

const t_inMap* edge::mesh::Regular::getInMap() const {
  return &m_inMap;
}

int_el edge::mesh::Regular::getNElements() const {
  return m_nRequestedElements;
}

int_el edge::mesh::Regular::getNFaces() const {
  return m_nRequestedFaces;
}

int_el edge::mesh::Regular::getNVertices() const {
  return m_nRequestedVertices;
}

void edge::mesh::Regular::getElementFaceNeighbors( int_el  i_elementId,
                                                   int_el *o_neighboringIds ) const {
  // line element
  if( m_elementType == Line ) {
    // left boundary
    if( i_elementId == 0 ) {
      o_neighboringIds[0] = m_nRequestedElements - 1;
      o_neighboringIds[1] = i_elementId + 1;
    }
    // right boundary
    else if( i_elementId == m_nRequestedElements - 1 ) {
      o_neighboringIds[0] = i_elementId - 1;
      o_neighboringIds[1] = 0;
    }
    // default case, not part of any boundary
    else {
      o_neighboringIds[0] = i_elementId - 1;
      o_neighboringIds[1] = i_elementId + 1;
    }
  }

  /*
   * Sorting of elements:
   *   ___________ _____
   *  |     |     |     |
   *  |  3  |  4  |  5  |
   *  |_____|_____|_____|
   *  |     |     |     |
   *  |  0  |  1  |  2  |
   *  |_____|_____|_____|
   *
   * Face Neighbors (counter-clockwise):
   *   _______________
   *  |      /|\      |
   *  |       2       |
   *  |               |
   *  |<- 3      1 -> |
   *  |               |
   *  |       0       |
   *  |______\|/______|
   *
   */
  else if( m_elementType == Quadrilateral ) {
    // compute x- and y-position from element id
    int_el l_posX = i_elementId % m_nX;
    int_el l_posY = i_elementId / m_nX;

    // set default cases
    o_neighboringIds[0] = ( (l_posY - 1) * m_nX ) + l_posX;     // bottom
    o_neighboringIds[1] = (  l_posY      * m_nX ) + l_posX + 1; // right
    o_neighboringIds[2] = ( (l_posY + 1) * m_nX ) + l_posX;     // top
    o_neighboringIds[3] = (  l_posY      * m_nX ) + l_posX - 1; // left

    // left boundary
    if( l_posX == 0 ) {
      o_neighboringIds[3] = (l_posY * m_nX) + m_nX - 1;
    }
    // right boundary
    else if( l_posX == m_nX-1 ) {
      o_neighboringIds[1] = (l_posY * m_nX);
    }

    // bottom boundary
    if( l_posY == 0 ) {
      o_neighboringIds[0] = ( (m_nY-1) * m_nX ) + l_posX;
    }
    // top boundary
    else if( l_posY == m_nY - 1 ) {
      o_neighboringIds[2] = l_posX;
    }
  }
  /*
   * Sorting of elements:
   *       __________________
   *      /  9  / 10  /  11 /|
   *     /_____/_____/_____/1|
   *    /  6  /  7  /  8  /|1|
   *   /____ /_____/_____/ | .
   *  |     |     |     | 8|/|
   *  |  6  |  7  |  8  | .| |
   *  |_____|_____|_____|/ |5.
   *  |     |     |     | 2|/
   *  |  0  |  1  |  2  | .
   *  |_____|_____|_____|/
   *
   * Face Neighbors:
   *   bottom: 0
   *   front:  1
   *   right:  2
   *   back:   3
   *   left:   4
   *   top:    5
   *
   *             x*******************x
   *            **                  **
   *           * *                 * *
   *          *  *          /     *  *
   *         *   *         3     *   *
   *        *    *              *    *
   *       x*******************x     *
   *       * -4  *             *  2- *
   *       *     x************ * ****x
   *       *    *              *    *
   *       *   *      .1       *   *
   *       |  /                *  *
   *  zeta | / eta             * *
   *       |/                  **
   *       x---****************x
   *         xi
   *
   */
  else if( m_elementType == Hexahedral ) {
    // compute x-, y- and z-position from element id
    int_el l_posX = i_elementId % m_nX;

    int_el l_posY = i_elementId % (m_nX*m_nY);
           l_posY = l_posY / m_nX;

    int_el l_posZ = i_elementId / (m_nX*m_nY);

    /*
     * default cases
     */
    o_neighboringIds[0] =   ( l_posZ-1 ) * m_nX * m_nY // bottom
                          +   l_posY     * m_nX
                          +   l_posX;

    o_neighboringIds[1] =     l_posZ     * m_nX * m_nY
                          + ( l_posY-1 ) * m_nX        // front
                          +   l_posX;

    o_neighboringIds[2] =     l_posZ     * m_nX * m_nY
                          +   l_posY     * m_nX
                          +   l_posX+1;                // right

    o_neighboringIds[3] =     l_posZ     * m_nX * m_nY
                          + ( l_posY+1 ) * m_nX        // back
                          +   l_posX;

    o_neighboringIds[4] =     l_posZ     * m_nX * m_nY
                          +   l_posY     * m_nX
                          +   l_posX-1;                // left

    o_neighboringIds[5] =   ( l_posZ+1 ) * m_nX * m_nY // top
                          +   l_posY     * m_nX
                          +   l_posX;

    /*
     * boundaries
     */
    if( l_posZ == 0 ) { // bottom
      o_neighboringIds[0] =   ( m_nZ-1 ) * m_nX * m_nY
                            +   l_posY   * m_nX
                            +   l_posX;
    }
    if( l_posY == 0 ) { // front
      o_neighboringIds[1] =     l_posZ   * m_nX * m_nY
                            + ( m_nY-1 ) * m_nX
                            +   l_posX;
    }
    if( l_posX == m_nX-1 ) { // right
      o_neighboringIds[2] =     l_posZ     * m_nX * m_nY
                            +   l_posY     * m_nX
                            +   0;
    }
    if( l_posY == m_nY-1 ) { // back
      o_neighboringIds[3] =     l_posZ     * m_nX * m_nY
                            +   0          * m_nX
                            +   l_posX;
    }
    if( l_posX == 0 ) { // left
      o_neighboringIds[4] =     l_posZ     * m_nX * m_nY
                            +   l_posY     * m_nX
                            +   m_nX-1;
    }
    if( l_posZ == m_nZ-1 ) { //top
      o_neighboringIds[5] =     0          * m_nX * m_nY
                            +   l_posY     * m_nX
                            +   l_posX;
    }
  }
}

void edge::mesh::Regular::getElementsAdjacentFaces( int_el (*o_elementAdjacentFaces)[C_ENT[T_SDISC.ELEMENT].N_FACES] ) const {
#if defined PP_T_ELEMENTS_LINE
  for( int_el l_element = 0; l_element < m_nRequestedElements; l_element++ ) {
    if( l_element < m_nRequestedElements-1 ) {
      o_elementAdjacentFaces[l_element][0] = l_element;
      o_elementAdjacentFaces[l_element][1] = l_element+1;
    }
    else {
      o_elementAdjacentFaces[l_element][0] = l_element;
      o_elementAdjacentFaces[l_element][1] = 0;
    }
  }
#elif defined PP_T_ELEMENTS_QUAD4R
  /*
   * Sorting of faces:
   *   ___________ _____
   *  |     |     |     |
   *  |6    |8    |10   |
   *  |__7__|__9__|__11_|
   *  |     |     |     |
   *  |0    |2    |4    |
   *  |__1__|__3__|__5__|
   */
  for( int_el l_element = 0; l_element < m_nRequestedElements; l_element++ ) {
    // "element-local" faces
    o_elementAdjacentFaces[l_element][3] = l_element*2; // left
    o_elementAdjacentFaces[l_element][0] = l_element*2 + 1; // bottom

    // "neighboring elements'" faces
    o_elementAdjacentFaces[l_element][1] = (l_element+1)*2; // right
    o_elementAdjacentFaces[l_element][2] = (l_element+m_nX)*2 + 1; // top

    // compute x- and y-position from element id
    int_el l_posX = l_element % m_nX;
    int_el l_posY = l_element / m_nX;

    // boundaries
    if( l_posX  == m_nX-1 ) { // right
      int_el l_neigh = l_posY * m_nX;
      o_elementAdjacentFaces[l_element][1] = l_neigh*2;
    }
    if( l_posY == m_nY-1 ) { // top
      int_el l_neigh = l_posX;
      o_elementAdjacentFaces[l_element][2] = l_neigh*2+1;
    }
  }
#elif defined PP_T_ELEMENTS_HEX8R
  /*
   * Sorting of faces:
   *   1st number: bottom
   *   2nd number: front
   *   3rd number: right
   *
   *               30,31 33,34
   *      27________32____35__
   *     28/.    /.    /.    /|
   *    29/_____/_____/_____/ |
   *     /  .  /  .  /  .  /| |
   *    /_____/_____/_____/ | |
   * 18|    .|21  .|24  .| _|/|15
   * 19|    .|22,23|25,26| /|_|16
   * 20|____ |_____|_____|/ | /17
   *   |  ._ |_ _ _|_ _._| _|/
   *   | .   | .   | .   | .
   *   |.____|.____|.____|/
   *    0,1,2 3,4,5 6,7,8
   *
   */
  for( int_el l_el = 0; l_el < m_nRequestedElements; l_el++ ) {
    // element local faces
    o_elementAdjacentFaces[l_el][0] = l_el*3;     // bottom
    o_elementAdjacentFaces[l_el][1] = l_el*3 + 1; // front
    o_elementAdjacentFaces[l_el][2] = l_el*3 + 2; // right

    // neighboring elements' faces
    o_elementAdjacentFaces[l_el][3] = ( l_el+(m_nX) )      * 3 + 1; // back
    o_elementAdjacentFaces[l_el][4] = ( l_el-1 )           * 3 + 2; // left
    o_elementAdjacentFaces[l_el][5] = ( l_el+(m_nX*m_nY) ) * 3 + 0; // top

    // derive element's position
    int_el l_posX = l_el % m_nX;
    int_el l_posY = l_el % (m_nX*m_nY);
           l_posY = l_posY / m_nX;
    int_el l_posZ = l_el / (m_nX*m_nY);

    // take care of periodic boundaries
    if( l_posX == 0 ) { // left bnd
      int_el l_neigh = l_el+m_nX-1;
      o_elementAdjacentFaces[l_el][4] = l_neigh*3 + 2;
    }
    if( l_posY == m_nY-1 ) { // back bnd
      int_el l_neigh = l_el - m_nX*(m_nY-1);
      o_elementAdjacentFaces[l_el][3] = l_neigh*3 + 1;
    }
    if( l_posZ == m_nZ-1 ) { // top bnd
      int_el l_neigh = l_el - m_nX*m_nY*(m_nZ-1);
      o_elementAdjacentFaces[l_el][5] = l_neigh*3 + 0;
    }
  }
#else
  assert( false );
#endif
}

int_el edge::mesh::Regular::getNelVeEl() const {
  if(      m_elementType == Line          ) return m_nRequestedElements  * 2;
  else if( m_elementType == Quadrilateral ) return m_nRequestedElements  * 8;
  else if( m_elementType == Hexahedral    ) return m_nRequestedElements * 26;

  EDGE_LOG_FATAL << "unsupported element type";
  return std::numeric_limits< int_el >::max();
}

void edge::mesh::Regular::getElVeEl( int_el** o_elVeEl ) const {
  // pointer to raw data
  int_el* l_raw = o_elVeEl[0];

  if( m_elementType == Line ) {
    for( int_el l_el = 0; l_el < m_nRequestedElements; l_el++ ) {
      l_raw[0] = l_el-1;
      l_raw[1] = l_el+1;
      o_elVeEl[l_el] = l_raw;
      l_raw += 2;
    }
    // adjust periodic boundaries
    o_elVeEl[0][0] = m_nRequestedElements-1;
    o_elVeEl[m_nRequestedElements-1][1] = 0;
  }
  else if( m_elementType == Quadrilateral ) {
    for( int_el l_el = 0; l_el < m_nRequestedElements; l_el++ ) {
      int_el l_px = l_el % m_nX;
      int_el l_py = l_el / m_nX;

      // compute offsets for perodic boundaries
      int l_offx[2] = {0,0};
      int l_offy[2] = {0,0};
      if(      l_px == 0      ) l_offx[0] =  (m_nX-1);
      else if( l_px == m_nX-1 ) l_offx[1] = -(m_nX-1);

      if(      l_py == 0      ) l_offy[0] =  m_nX*(m_nY-1);
      else if( l_py == m_nY-1 ) l_offy[1] = -m_nX*(m_nY-1);

      // set vertex-neighbors, starting bottom-left, counter-clockwise
      l_raw[0] = ( (l_py-1) * m_nX ) + l_px - 1 + l_offx[0] + l_offy[0];
      l_raw[1] = ( (l_py-1) * m_nX ) + l_px + 0 + 0         + l_offy[0];
      l_raw[2] = ( (l_py-1) * m_nX ) + l_px + 1 + l_offx[1] + l_offy[0];

      l_raw[3] = ( (l_py+0) * m_nX ) + l_px + 1 + l_offx[1] + 0;
      l_raw[4] = ( (l_py+1) * m_nX ) + l_px + 1 + l_offx[1] + l_offy[1];
      l_raw[5] = ( (l_py+1) * m_nX ) + l_px + 0 + 0         + l_offy[1];
      l_raw[6] = ( (l_py+1) * m_nX ) + l_px - 1 + l_offx[0] + l_offy[1];

      l_raw[7] = ( (l_py+0) * m_nX ) + l_px - 1 + l_offx[0] + 0;

      o_elVeEl[l_el] = l_raw;
      l_raw += 8;
    }

  }
  else if( m_elementType == Hexahedral ) {
    for( int_el l_el = 0; l_el < m_nRequestedElements; l_el++ ) {
      // derive element's position
      int_el l_px = l_el % m_nX;
      int_el l_py = l_el % (m_nX*m_nY);
             l_py = l_py / m_nX;
      int_el l_pz = l_el / (m_nX*m_nY);

      // compute offsets for perodic boundaries
      int l_offx[2] = {0,0};
      int l_offy[2] = {0,0};
      int l_offz[2] = {0,0};
      if(      l_px == 0      ) l_offx[0] =  (m_nX-1);
      else if( l_px == m_nX-1 ) l_offx[1] = -(m_nX-1);

      if(      l_py == 0      ) l_offy[0] =  m_nX*(m_nY-1);
      else if( l_py == m_nY-1 ) l_offy[1] = -m_nX*(m_nY-1);

      if(      l_pz == 0      ) l_offz[0] =  m_nX*m_nY*(m_nZ-1);
      else if( l_pz == m_nZ-1 ) l_offz[1] = -m_nX*m_nY*(m_nZ-1);

      // set vertex neighbors, starting lower-bottom-left
      l_raw[ 0] =  ( (l_pz-1) * m_nX * m_nY ) + ( (l_py-1) * m_nX ) + (l_px-1) + l_offx[0] + l_offy[0] + l_offz[0];
      l_raw[ 1] =  ( (l_pz-1) * m_nX * m_nY ) + ( (l_py-1) * m_nX ) + (l_px  ) + 0         + l_offy[0] + l_offz[0];
      l_raw[ 2] =  ( (l_pz-1) * m_nX * m_nY ) + ( (l_py-1) * m_nX ) + (l_px+1) + l_offx[1] + l_offy[0] + l_offz[0];

      l_raw[ 3] =  ( (l_pz-1) * m_nX * m_nY ) + ( (l_py  ) * m_nX ) + (l_px-1) + l_offx[0] + 0         + l_offz[0];
      l_raw[ 4] =  ( (l_pz-1) * m_nX * m_nY ) + ( (l_py  ) * m_nX ) + (l_px  ) + 0         + 0         + l_offz[0];
      l_raw[ 5] =  ( (l_pz-1) * m_nX * m_nY ) + ( (l_py  ) * m_nX ) + (l_px+1) + l_offx[1] + 0         + l_offz[0];

      l_raw[ 6] =  ( (l_pz-1) * m_nX * m_nY ) + ( (l_py+1) * m_nX ) + (l_px-1) + l_offx[0] + l_offy[1] + l_offz[0];
      l_raw[ 7] =  ( (l_pz-1) * m_nX * m_nY ) + ( (l_py+1) * m_nX ) + (l_px  ) + 0         + l_offy[1] + l_offz[0];
      l_raw[ 8] =  ( (l_pz-1) * m_nX * m_nY ) + ( (l_py+1) * m_nX ) + (l_px+1) + l_offx[1] + l_offy[1] + l_offz[0];


      l_raw[ 9] =  ( (l_pz  ) * m_nX * m_nY ) + ( (l_py-1) * m_nX ) + (l_px-1) + l_offx[0] + l_offy[0] + 0;
      l_raw[10] =  ( (l_pz  ) * m_nX * m_nY ) + ( (l_py-1) * m_nX ) + (l_px  ) + 0         + l_offy[0] + 0;
      l_raw[11] =  ( (l_pz  ) * m_nX * m_nY ) + ( (l_py-1) * m_nX ) + (l_px+1) + l_offx[1] + l_offy[0] + 0;

      l_raw[12] =  ( (l_pz  ) * m_nX * m_nY ) + ( (l_py  ) * m_nX ) + (l_px-1) + l_offx[0] + 0         + 0;

      l_raw[13] =  ( (l_pz  ) * m_nX * m_nY ) + ( (l_py  ) * m_nX ) + (l_px+1) + l_offx[1] + 0         + 0;

      l_raw[14] =  ( (l_pz  ) * m_nX * m_nY ) + ( (l_py+1) * m_nX ) + (l_px-1) + l_offx[0] + l_offy[1] + 0;
      l_raw[15] =  ( (l_pz  ) * m_nX * m_nY ) + ( (l_py+1) * m_nX ) + (l_px  ) + 0         + l_offy[1] + 0;
      l_raw[16] =  ( (l_pz  ) * m_nX * m_nY ) + ( (l_py+1) * m_nX ) + (l_px+1) + l_offx[1] + l_offy[1] + 0;


      l_raw[17] =  ( (l_pz+1) * m_nX * m_nY ) + ( (l_py-1) * m_nX ) + (l_px-1) + l_offx[0] + l_offy[0] + l_offz[1];
      l_raw[18] =  ( (l_pz+1) * m_nX * m_nY ) + ( (l_py-1) * m_nX ) + (l_px  ) + 0         + l_offy[0] + l_offz[1];
      l_raw[19] =  ( (l_pz+1) * m_nX * m_nY ) + ( (l_py-1) * m_nX ) + (l_px+1) + l_offx[1] + l_offy[0] + l_offz[1];

      l_raw[20] =  ( (l_pz+1) * m_nX * m_nY ) + ( (l_py  ) * m_nX ) + (l_px-1) + l_offx[0] + 0         + l_offz[1];
      l_raw[21] =  ( (l_pz+1) * m_nX * m_nY ) + ( (l_py  ) * m_nX ) + (l_px  ) + 0         + 0         + l_offz[1];
      l_raw[22] =  ( (l_pz+1) * m_nX * m_nY ) + ( (l_py  ) * m_nX ) + (l_px+1) + l_offx[1] + 0         + l_offz[1];

      l_raw[23] =  ( (l_pz+1) * m_nX * m_nY ) + ( (l_py+1) * m_nX ) + (l_px-1) + l_offx[0] + l_offy[1] + l_offz[1];
      l_raw[24] =  ( (l_pz+1) * m_nX * m_nY ) + ( (l_py+1) * m_nX ) + (l_px  ) + 0         + l_offy[1] + l_offz[1];
      l_raw[25] =  ( (l_pz+1) * m_nX * m_nY ) + ( (l_py+1) * m_nX ) + (l_px+1) + l_offx[1] + l_offy[1] + l_offz[1];

      o_elVeEl[l_el] = l_raw;
      l_raw += 26;
    }
  }
  else EDGE_LOG_FATAL;

  // set ghost entry
  o_elVeEl[m_nRequestedElements] = l_raw;
}

void edge::mesh::Regular::getElementsFaceNeighbors( int_el (*o_neighboringIds)[C_ENT[T_SDISC.ELEMENT].N_FACES] ) const {
  for( int_el l_element = 0; l_element < m_nRequestedElements; l_element++ ) {
    getElementFaceNeighbors( l_element, o_neighboringIds[l_element] );
  }
}

void edge::mesh::Regular::getFaceAdjacentElements( int_el i_faceId,
                                                   int_el o_neighboringIds[2] ) const {
  // line element
#if defined PP_T_ELEMENTS_LINE
  // left or right boundary
  if( i_faceId == 0 ) {
    o_neighboringIds[0] = m_nRequestedElements - 1;
    o_neighboringIds[1] = 0;
  }
  else {
    o_neighboringIds[0] = i_faceId - 1;
    o_neighboringIds[1] = i_faceId;
  }
#elif defined PP_T_ELEMENTS_QUAD4R
  // derive face type. 0: left/right, 1: bottom/top
  unsigned int l_type = i_faceId % 2;

  // derive right/top element
  o_neighboringIds[1] = (i_faceId - l_type) / 2;

  // derive left/bottom element
  int_el l_neighboringIds[C_ENT[T_SDISC.ELEMENT].N_FACES];

  mesh::Regular::getElementFaceNeighbors( o_neighboringIds[1],
                                          l_neighboringIds );

  if( l_type == 0 ) {
    o_neighboringIds[0] = l_neighboringIds[3];
  }
  else {
    o_neighboringIds[0] = l_neighboringIds[0];
  }
#elif defined PP_T_ELEMENTS_HEX8R
  // derive face type. 0: z, 1: y, 2: x
  unsigned int l_type = i_faceId % 3;

  // top,back,right equals the face's associated element
  o_neighboringIds[1] = (i_faceId - l_type) / (int_el) 3;

  // derive bottom/front/left element
  int_el l_neighboringIds[C_ENT[T_SDISC.ELEMENT].N_FACES];

  mesh::Regular::getElementFaceNeighbors( o_neighboringIds[1],
                                          l_neighboringIds );

  o_neighboringIds[0] = l_neighboringIds[l_type];
#else
  assert(false);
#endif
  // left should'nt be equal to right..
  EDGE_CHECK( o_neighboringIds[0] != o_neighboringIds[1] );

  // store smaller element first
  if( o_neighboringIds[0] > o_neighboringIds[1] ) {
    int_el l_tmpId = o_neighboringIds[0];
    o_neighboringIds[0] = o_neighboringIds[1];
    o_neighboringIds[1] = l_tmpId;
  }
}

void edge::mesh::Regular::getFaceAdjacentVertices( int_el i_fa, int_el o_faVe[C_ENT[T_SDISC.ELEMENT].N_FACE_VERTICES] ) const {
#if defined PP_T_ELEMENTS_LINE
  o_faVe[0] = i_fa;
#elif defined PP_T_ELEMENTS_QUAD4R
  // determine adjacent elements
  int_el l_elAd[2];
  getFaceAdjacentElements( i_fa, l_elAd );

  // get vertices of right element
  int_el l_elVe[4];

  getElementAdjacentVertices( l_elVe,
                              l_elAd[1]%m_nX,
                              l_elAd[1]/m_nX );

  // assign vertices for vertical faces
  if( i_fa%2 == 0 ) {
    o_faVe[0] = l_elVe[0];
    o_faVe[1] = l_elVe[3];
  }
  // assign vertices for horizontal faces
  else {
    o_faVe[0] = l_elVe[0];
    o_faVe[1] = l_elVe[1];
  }
#elif defined PP_T_ELEMENTS_HEX8R
  // derive face type. 0: z, 1: y, 2: x
  unsigned short l_type = i_fa % 3;

  // top,back,right equals the face's associated element
  int_el l_el = (i_fa - l_type) / (int_el) 3;

  // derive element's position
  int_el l_x = l_el % m_nX;
  int_el l_y = l_el % (m_nX*m_nY);
         l_y = l_y / m_nX;
  int_el l_z = l_el / (m_nX*m_nY);

  // derive vertex ids of the element
  int_el l_veId0 =   l_z * (m_nX+1) * (m_nY+1)
                   + l_y * (m_nX+1)
                   + l_x;
  int_el l_veId1 = l_veId0 +            1;
  int_el l_veId2 = l_veId0 + (m_nX+1) + 1;
  int_el l_veId3 = l_veId0 + (m_nX+1) + 0;

  int_el l_veId4 = l_veId0 + (m_nX+1) * (m_nY+1);
  int_el l_veId5 = l_veId1 + (m_nX+1) * (m_nY+1);
  int_el l_veId6 = l_veId2 + (m_nX+1) * (m_nY+1);

  // determine the corresponding vertices based on the face type
  if( l_type == 0 ) { // bottom face of the element
    o_faVe[0] = l_veId0;
    o_faVe[1] = l_veId1;
    o_faVe[2] = l_veId2;
    o_faVe[3] = l_veId3;
  }
  else if( l_type == 1 ) { // front face of the element
    o_faVe[0] = l_veId0;
    o_faVe[1] = l_veId1;
    o_faVe[2] = l_veId5;
    o_faVe[3] = l_veId4;
  }
  else if( l_type == 2 ) { // right face of the element
    o_faVe[0] = l_veId1;
    o_faVe[1] = l_veId5;
    o_faVe[2] = l_veId6;
    o_faVe[3] = l_veId2;
  }
#else
  EDGE_LOG_FATAL;
#endif
}

void edge::mesh::Regular::getFacesAdjacentVertices( int_el (*o_faVe)[C_ENT[T_SDISC.ELEMENT].N_FACE_VERTICES] ) const {
  for( int_el l_fa = 0; l_fa <  m_nRequestedFaces; l_fa++ ) {
    getFaceAdjacentVertices( l_fa, o_faVe[l_fa] );
  }
}

void edge::mesh::Regular::getFacesAdjacentElements( int_el (*o_neighboringIds)[2] ) const {
  for( int_el l_face = 0; l_face < m_nRequestedFaces; l_face++ ) {
    getFaceAdjacentElements( l_face, o_neighboringIds[l_face] );
  }
}

void edge::mesh::Regular::getElementAdjacentVertices( int_el o_elementAdjacentVertices[C_ENT[T_SDISC.ELEMENT].N_VERTICES],
                                                      int_el i_px,
                                                      int_el i_py,
                                                      int_el i_pz ) const {
#if defined PP_T_ELEMENTS_LINE
  /*
   * Sorting of vertices:
   *  |0____|1____|2____|3
   */
  o_elementAdjacentVertices[0] = i_px;
  o_elementAdjacentVertices[1] = i_px+1;
#elif defined PP_T_ELEMENTS_QUAD4R

  /*
   * Sorting of vertices:
   *   8_____9____ 10___11
   *  |     |     |     |
   *  |     |     |     |
   *  |4____|5____|6____|7
   *  |     |     |     |
   *  |     |     |     |
   *  |0____|1____|2____|3
   */
  // derive vertex ids
  int_el l_veId0 =  i_py      * (m_nX+1) + i_px;
  int_el l_veId1 =  i_py      * (m_nX+1) + i_px + 1;
  int_el l_veId2 = (i_py + 1) * (m_nX+1) + i_px;
  int_el l_veId3 = (i_py + 1) * (m_nX+1) + i_px + 1;

  // set result, using counter-clockwise orientation
  // side effect: keep this ordering for visulization purposes (no diagonals) and high-order configs
  o_elementAdjacentVertices[0] = l_veId0;
  o_elementAdjacentVertices[1] = l_veId1;
  o_elementAdjacentVertices[2] = l_veId3;
  o_elementAdjacentVertices[3] = l_veId2;
#elif defined PP_T_ELEMENTS_HEX8R
  /*
   * Sorting of vertices:
   *
   *      32______33____34____35
   *       /.    /.    /.    /|
   *    28/___29/___30/___31/ |
   *     /  .  /  .  /  .  /| |
   * 24 /___25/___26/___27/ | |23
   *   |    .|    .|    .| _|/|
   *   |    .|_ _ _|_ _ _| /|_|11
   * 12|___13|___14|___15|/ | /
   *   |  ._ |_ _ _|_ _._| _|/
   *   | .   | .   | .   | . 7
   *   |.____|.____|.____|/
   *  0      1     2      3
   */
  // derive vertex ids
  int_el l_veId0 =   i_pz * (m_nX+1) * (m_nY+1)
                   + i_py * (m_nX+1)
                   + i_px;
  int_el l_veId1 = l_veId0 +            1;
  int_el l_veId2 = l_veId0 + (m_nX+1) + 1;
  int_el l_veId3 = l_veId0 + (m_nX+1) + 0;

  int_el l_veId4 = l_veId0 + (m_nX+1) * (m_nY+1);
  int_el l_veId5 = l_veId1 + (m_nX+1) * (m_nY+1);
  int_el l_veId6 = l_veId2 + (m_nX+1) * (m_nY+1);
  int_el l_veId7 = l_veId3 + (m_nX+1) * (m_nY+1);

  // set result
  o_elementAdjacentVertices[0] = l_veId0;
  o_elementAdjacentVertices[1] = l_veId1;
  o_elementAdjacentVertices[2] = l_veId2;
  o_elementAdjacentVertices[3] = l_veId3;
  o_elementAdjacentVertices[4] = l_veId4;
  o_elementAdjacentVertices[5] = l_veId5;
  o_elementAdjacentVertices[6] = l_veId6;
  o_elementAdjacentVertices[7] = l_veId7;
#else
  assert( false );
#endif
}

void edge::mesh::Regular::getElementsAdjacentVertices( int_el (*o_elementAdjacentVertices)[C_ENT[T_SDISC.ELEMENT].N_VERTICES] ) const {
#if defined PP_T_ELEMENTS_LINE
  for( int_el l_el = 0; l_el < m_nRequestedElements; l_el++ ) {
    getElementAdjacentVertices(  o_elementAdjacentVertices[l_el],
                                 l_el );
  }
#elif defined PP_T_ELEMENTS_QUAD4R
  // iterate over elements dimension-wise
  for( int_el l_y = 0; l_y < m_nY; l_y++ ) {
    for( int_el l_x = 0; l_x < m_nX; l_x++ ) {
      // derive element id
      int_el l_elId = l_y * m_nX + l_x;

      getElementAdjacentVertices( o_elementAdjacentVertices[l_elId],
                                  l_x, l_y );
    }
  }
#elif defined PP_T_ELEMENTS_HEX8R
  /*
   * Sorting of vertices:
   *
   *      32______33____34____35
   *       /.    /.    /.    /|
   *    28/___29/___30/___31/ |
   *     /  .  /  .  /  .  /| |
   * 24 /___25/___26/___27/ | |23
   *   |    .|    .|    .| _|/|
   *   |    .|_ _ _|_ _ _| /|_|11
   * 12|___13|___14|___15|/ | /
   *   |  ._ |_ _ _|_ _._| _|/
   *   | .   | .   | .   | . 7
   *   |.____|.____|.____|/
   *  0      1     2      3
   */
  // iterate over elements
  for( int_el l_z = 0; l_z < m_nZ; l_z++ ) {
    for( int_el l_y = 0; l_y < m_nY; l_y++ ) {
      for( int_el l_x = 0; l_x < m_nX; l_x++ ) {
        // derive element id
        int_el l_elId = l_z * m_nX * m_nY +
                        l_y * m_nX +
                        l_x;

        getElementAdjacentVertices( o_elementAdjacentVertices[l_elId],
                                    l_x, l_y, l_z );
      }
    }
  }
#else
  assert( false );
#endif
}

void edge::mesh::Regular::getConnect( const t_vertexChars *i_veChars,
                                      const t_faceChars   *i_faChars,
                                            t_connect     &o_connect ) const {
  // get the connectivity info
  getElementsAdjacentVertices( o_connect.elVe   );
  getElementsAdjacentFaces(    o_connect.elFa   );
  getFacesAdjacentVertices(    o_connect.faVe   );
  getFacesAdjacentElements(    o_connect.faEl   );
  //getElVeEl(                   o_connect.elVeEl );
  getElementsFaceNeighbors(    o_connect.elFaEl );

  t_enLayout l_veLayout, l_faLayout, l_elLayout;
  l_veLayout = getVeLayout(); l_faLayout = getFaLayout(); l_elLayout = getElLayout();

  // check the consitency of the infos
  common< T_SDISC.ELEMENT >::checkConsAdj( l_faLayout,
                                           i_faChars,
                                           m_inMap.faMeDa,
                                           m_inMap.faDaMe,
                                           o_connect.elFaEl,
                                           o_connect.elFa,
                                           o_connect.faEl );

  // get the neighboring elements' local face ids
  common< T_SDISC.ELEMENT >::getFIdsElFaEl( l_elLayout, m_gIdsEl, o_connect.elFa, o_connect.elFaEl, o_connect.fIdElFaEl );

  // get the neighboring elements' local vertex ids, which match the faces' first vertex
  common< T_SDISC.ELEMENT >::getVIdsElFaEl( l_elLayout,
                                            m_gIdsVe,
                                            o_connect.elFaEl,
                                            o_connect.elVe,
                                            o_connect.fIdElFaEl,
                                            o_connect.vIdElFaEl );
}

void edge::mesh::Regular::computeElementVolume( int_el     i_element,
                                                real_mesh &o_volume ) const {
  if( m_elementType == Line ) {
    o_volume = m_dX;
  }
  else if( m_elementType == Quadrilateral ) {
    o_volume = m_dX * m_dY;
  }
  else if( m_elementType == Hexahedral ) {
    o_volume = m_dX * m_dY * m_dZ;
  }
  else {
    assert( false );
  }
}

void edge::mesh::Regular::computeElementInDia( int_el     i_element,
                                               real_mesh &o_dia ) const {
  if( m_elementType == Line ) {
    o_dia = m_dX;
  }
  else if( m_elementType == Quadrilateral ) {
    o_dia = std::min( m_dX, m_dY );
  }
  else if( m_elementType == Hexahedral ) {
    o_dia = std::min( m_dX, m_dY );
    o_dia = std::min( o_dia, m_dZ );
  }
  else assert(false);
}

void edge::mesh::Regular::computeFaceArea( int_el    i_face,
                                           real_mesh &o_area ) const {
  if( m_elementType == Line ) {
    // default to 1, which matches contributions in Riemann solvers
    o_area = 1;
  }
  else if( m_elementType == Quadrilateral ) {
    if( i_face % 2 == 0 ) o_area = m_dY; // vertical face
    else                  o_area = m_dX; // horizontal face
  }
  else if( m_elementType == Hexahedral ) {
    unsigned short l_faType = i_face%3;

    if(      l_faType == 0 ) o_area = m_dX * m_dY; // bottom
    else if( l_faType == 1 ) o_area = m_dX * m_dZ; // front
    else if( l_faType == 2 ) o_area = m_dY * m_dZ; // right
    else assert( false );
  }
  else assert(false);
}

void edge::mesh::Regular::getElChars( t_elementChars* o_elementChars ) const {
  for( int_el l_el = 0; l_el < m_nRequestedElements; l_el++ ) {
    computeElementVolume( l_el,
                          o_elementChars[l_el].volume );
    computeElementInDia( l_el,
                         o_elementChars[l_el].inDia );

    o_elementChars[l_el].spType = MESH_TYPE_NONE;
  }
}

void edge::mesh::Regular::computeOutPointNormal( int_el     i_face,
                                                 real_mesh &o_x,
                                                 real_mesh &o_y,
                                                 real_mesh &o_z ) const {
  if( m_elementType == Line ) { 
    o_x = 1; o_y = 0; o_z = 0;
  }
  else if( m_elementType == Quadrilateral ) {
    if( i_face % 2 == 0 ) { // vertical face
      o_x = 1; o_y = 0; o_z = 0;
    }
    else { // horizontal face
      o_x = 0; o_y = 1; o_z = 0;
    }
  }
  else if( m_elementType == Hexahedral ) {
    unsigned short l_faType = i_face % (int_el) 3;

    // all outer-pointing normals point in the direction of the canonical basis
    if( l_faType  == 0 ) { // bottom
      o_x = 0; o_y = 0; o_z = 1;
    }
    else if( l_faType == 1 ) {  // front
      o_x = 0; o_y = 1; o_z = 0;
    }
    else if( l_faType == 2 ) { // right
      o_x = 1; o_y = 0; o_z = 0;
    }
    else assert( false );
  }
  else assert(false);
}

void edge::mesh::Regular::computeTangents( int_el i_face,
                                           real_mesh o_tangent0[3],
                                           real_mesh o_tangent1[3] ) const {
  // reset
  for( unsigned short l_dim = 0; l_dim < 3; l_dim++ ) {
    o_tangent0[l_dim] = o_tangent1[l_dim] = 0.0;
  }

  if( m_elementType == Line ) {
    // no tangents for line elements
  }
  else if( m_elementType == Quadrilateral ) {
    if( i_face % 2 == 0 ) { // vertical face
      o_tangent0[1] = 1.0;
    }
    else { // horizontal face
      o_tangent0[0] = 1.0;
    }
  }
  else if( m_elementType == Hexahedral ) {
    unsigned short l_faType = i_face % (int_el) 3;

    // all outer-pointing normals point into the associated element of the face
    if( l_faType  == 0 ) { // bottom
      o_tangent0[0] = 1.0;
      o_tangent1[1] = 1.0;
    }
    else if( l_faType == 1 ) {  // front
      o_tangent0[0] = 1.0;
      o_tangent1[2] = 1.0;
    }
    else if( l_faType == 2 ) { //right
      o_tangent0[1] = 1.0;
      o_tangent1[2] = 1.0;
    }
    else assert( false );
  }
  else assert(false);
}

void edge::mesh::Regular::getFaChars( t_faceChars* o_faceChars ) const {
  for( int_el l_fa = 0; l_fa < m_nRequestedFaces; l_fa++ ) {
    // compute face characteristics
    computeFaceArea( l_fa,
                     o_faceChars[l_fa].area );
    computeOutPointNormal( l_fa,
                           o_faceChars[l_fa].outNormal[0],
                           o_faceChars[l_fa].outNormal[1],
                           o_faceChars[l_fa].outNormal[2] );

    // derive associated element
    if( m_elementType == Line ){
      if( l_fa == 0 ) o_faceChars[l_fa].outNormal[0] *= -1;
    }
    else if( m_elementType == Quadrilateral ) {
      int_el l_el = l_fa / 2;

      // derive position
      int_el l_posX = l_el % m_nX;
      int_el l_posY = l_el / m_nX;

      int_el l_ty = l_fa % 2;

      // invert if at periodic boundary to ensure normals, pointing from the smaller to the larger element-id
      if(    (l_posX == 0 && l_ty == 0)
          || (l_posY == 0 && l_ty == 1) ) {
        for( unsigned short l_di = 0; l_di < 2; l_di++ )
          o_faceChars[l_fa].outNormal[l_di] *= -1;
      }
    }
    else if( m_elementType == Hexahedral ) {
      int_el l_el = l_fa / 3;

      // derive position
      int_el l_posX = l_el % m_nX;
      int_el l_posY = l_el % (m_nX*m_nY);
            l_posY = l_posY / m_nX;
      int_el l_posZ = l_el / (m_nX*m_nY);

      int_el l_ty = l_fa % 3;

      // invert if at periodic boundary to ensure normals, pointing from the smaller to the larger element-id
      if(    (l_posX == m_nX-1 && l_ty == 2)
          || (l_posY == 0      && l_ty == 1)
          || (l_posZ == 0      && l_ty == 0) ) {
        for( unsigned short l_di = 0; l_di < 3; l_di++ )
          o_faceChars[l_fa].outNormal[l_di] *= -1;
      }
    }
    else EDGE_LOG_FATAL;

    computeTangents( l_fa,
                     o_faceChars[l_fa].tangent0,
                     o_faceChars[l_fa].tangent1 );

    o_faceChars[l_fa].spType = MESH_TYPE_NONE;
  }
}

void edge::mesh::Regular::getVeChars( t_vertexChars* o_vertexChars ) const {
#if defined PP_T_ELEMENTS_LINE
  for( int_el l_x = 0; l_x < m_nX+1; l_x++ ) {
    o_vertexChars[l_x].coords[0] = l_x * m_dX;
    o_vertexChars[l_x].coords[1] = 0;
    o_vertexChars[l_x].coords[2] = 0;

    o_vertexChars[l_x].spType = MESH_TYPE_NONE;
  }
#elif defined PP_T_ELEMENTS_QUAD4R
  // iterate over vertices dimension-wise
  for( int_el l_y = 0; l_y < m_nY+1; l_y++ ) {
    for( int_el l_x = 0; l_x < m_nX+1; l_x++ ) {
      // derive vertex id
      int_el l_veId = l_y*(m_nX+1) + l_x;

      // derive positions
      o_vertexChars[l_veId].coords[0] = (real_mesh) l_x * m_dX;
      o_vertexChars[l_veId].coords[1] = (real_mesh) l_y * m_dY;
      o_vertexChars[l_veId].coords[2] = 0;

      o_vertexChars[l_veId].spType = MESH_TYPE_NONE;
    }
  }
#elif defined PP_T_ELEMENTS_HEX8R
  // iterate over vertices dimension-wise
  for( int_el l_z = 0; l_z < m_nZ+1; l_z++ ) {
    for( int_el l_y = 0; l_y < m_nY+1; l_y++ ) {
      for( int_el l_x = 0; l_x < m_nX+1; l_x++ ) {
        // derive vertex id
        int_el l_veId = l_z * (m_nX+1) * (m_nY+1) +
                        l_y * (m_nX+1)            +
                        l_x;

        // derive positions
        o_vertexChars[l_veId].coords[0] = (real_mesh) l_x * m_dX;
        o_vertexChars[l_veId].coords[1] = (real_mesh) l_y * m_dY;
        o_vertexChars[l_veId].coords[2] = (real_mesh) l_z * m_dZ;

        o_vertexChars[l_veId].spType = MESH_TYPE_NONE;
      }
    }
  }
#else
  assert( false );
#endif
}
