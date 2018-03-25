/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2015-2016, Regents of the University of California
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
 * Unstructured mesh representation.
 **/

#include "Unstructured.h"
#include "common.hpp"
#include "linalg/Geom.hpp"
#include <io/logging.h>
#include <cassert>

edge::mesh::Unstructured::Unstructured(       unsigned int  i_nBndVals,
                                        const int          *i_bndVals,
                                              int           i_periodicVal ):
#ifdef PP_USE_MOAB
 m_moab( N_DIM,
         i_nBndVals,
         i_bndVals,
         i_periodicVal )
#endif
{
  if( i_nBndVals > 0 && i_bndVals == NULL ) {
    EDGE_LOG_FATAL << "non-zero number of bnd vals given, but no values: " << i_nBndVals;
  }

#ifdef PP_USE_MPI
  if( i_periodicVal != std::numeric_limits<int>::max() ) {
    EDGE_LOG_FATAL << "trying to use periodic boundaries with MPI enabled.";
  }
#endif
}

void edge::mesh::Unstructured::printStats() {
  EDGE_LOG_INFO << "check it out: mesh statistics; element: full dim, face: dim-1, vertex: dim-2";
#ifdef PP_USE_MPI
  unsigned int l_nEn[3][3][3];

  // tmp values, 0: entities, 1: mpi-type
  unsigned int l_tmp[3][3];

  // iterate over entities (ve, fa, el)
  for( unsigned short l_en = 0; l_en < 3; l_en++ ) {
    // iterate over mpi-types (inner, send, receive)
    for( unsigned short l_mt = 0; l_mt < 3; l_mt++ ) {
      l_tmp[l_en][l_mt] = 0;
      // iterate over statistics (min, ave, max );
      for( unsigned short l_st = 0; l_st < 3; l_st++ ) {
        l_nEn[l_en][l_mt][l_st] = 0;
      }
    }
  }

  // get entity layouts
  t_enLayout l_veL = getVeLayout();
  t_enLayout l_faL = getFaLayout();
  t_enLayout l_elL = getElLayout();

  for( int_tg l_tg = 0; l_tg < l_veL.timeGroups.size(); l_tg++ ) {
    l_tmp[0][0] += l_veL.timeGroups[l_tg].inner.size;
    l_tmp[0][1] += l_veL.timeGroups[l_tg].nEntsOwn - l_veL.timeGroups[l_tg].inner.size;
    l_tmp[0][2] += l_veL.timeGroups[l_tg].nEntsNotOwn;

    l_tmp[1][0] += l_faL.timeGroups[l_tg].inner.size;
    l_tmp[1][1] += l_faL.timeGroups[l_tg].nEntsOwn - l_faL.timeGroups[l_tg].inner.size;
    l_tmp[1][2] += l_faL.timeGroups[l_tg].nEntsNotOwn;

    l_tmp[2][0] += l_elL.timeGroups[l_tg].inner.size;
    l_tmp[2][1] += l_elL.timeGroups[l_tg].nEntsOwn - l_elL.timeGroups[l_tg].inner.size;
    l_tmp[2][2] += l_elL.timeGroups[l_tg].nEntsNotOwn;
  }

  // get the statistics
  for( unsigned short l_en = 0; l_en < 3; l_en++ ) {
    for( unsigned short l_mt = 0; l_mt < 3; l_mt++ ) {
      int l_err;
      l_err = MPI_Allreduce( &l_tmp[l_en][l_mt], l_nEn[l_en][l_mt]+0, 1, MPI_UNSIGNED, MPI_MIN, MPI_COMM_WORLD );
      EDGE_CHECK_EQ( l_err, MPI_SUCCESS );
      l_err = MPI_Allreduce( &l_tmp[l_en][l_mt], l_nEn[l_en][l_mt]+1, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD );
      EDGE_CHECK_EQ( l_err, MPI_SUCCESS );
      l_err = MPI_Allreduce( &l_tmp[l_en][l_mt], l_nEn[l_en][l_mt]+2, 1, MPI_UNSIGNED, MPI_MAX, MPI_COMM_WORLD );
      EDGE_CHECK_EQ( l_err, MPI_SUCCESS );
    }
  }

  // share our wisdom
  EDGE_LOG_INFO << "  #vertices:";
  EDGE_LOG_INFO << "    inner (min/ave/max): " << l_nEn[0][0][0] << " / "
                                               << l_nEn[0][0][1] / (double) parallel::g_nRanks << " / "
                                               << l_nEn[0][0][2];
  EDGE_LOG_INFO << "    send  (min/ave/max): " << l_nEn[0][1][0] << " / "
                                               << l_nEn[0][1][1] / (double) parallel::g_nRanks << " / "
                                               << l_nEn[0][1][2];
  EDGE_LOG_INFO << "    recv  (min/ave/max): " << l_nEn[0][2][0] << " / "
                                               << l_nEn[0][2][1] / (double) parallel::g_nRanks << " / "
                                               << l_nEn[0][2][2];
  EDGE_LOG_INFO << "  #faces:";
  EDGE_LOG_INFO << "    inner (min/ave/max): " << l_nEn[1][0][0] << " / "
                                               << l_nEn[1][0][1] / (double) parallel::g_nRanks << " / "
                                               << l_nEn[1][0][2];
  EDGE_LOG_INFO << "    send  (min/ave/max): " << l_nEn[1][1][0] << " / "
                                               << l_nEn[1][1][1] / (double) parallel::g_nRanks << " / "
                                               << l_nEn[1][1][2];
  EDGE_LOG_INFO << "    recv  (min/ave/max): " << l_nEn[1][2][0] << " / "
                                               << l_nEn[1][2][1] / (double) parallel::g_nRanks << " / "
                                               << l_nEn[1][2][2];
  EDGE_LOG_INFO << "  #elements:";
  EDGE_LOG_INFO << "    inner (min/ave/max): " << l_nEn[2][0][0] << " / "
                                               << l_nEn[2][0][1] / (double) parallel::g_nRanks << " / "
                                               << l_nEn[2][0][2];
  EDGE_LOG_INFO << "    send  (min/ave/max): " << l_nEn[2][1][0] << " / "
                                               << l_nEn[2][1][1] / (double) parallel::g_nRanks << " / "
                                               << l_nEn[2][1][2];
  EDGE_LOG_INFO << "    recv  (min/ave/max): " << l_nEn[2][2][0] << " / "
                                               << l_nEn[2][2][1] / (double) parallel::g_nRanks << " / "
                                               << l_nEn[2][2][2];
#else
  EDGE_LOG_INFO << "  #vertices: " << getNVertices();
  EDGE_LOG_INFO << "  #faces: "    << getNFaces();
  EDGE_LOG_INFO << "  #elements: " << getNElements();
#endif

#ifdef PP_USE_MOAB
  std::vector<int> l_meshRegions;
  m_moab.getMeshRegions( l_meshRegions );

  if( l_meshRegions.size() > 0 ) {
    EDGE_LOG_INFO << "  interested in mesh regions, hu? here's our data";
    for( unsigned int l_region = 0; l_region < l_meshRegions.size(); l_region++ ) {
      EDGE_LOG_INFO << "    region #" << l_region << ": " << l_meshRegions[l_region];
    }
  }
#endif

}

void edge::mesh::Unstructured::read( const std::string &i_pathToMesh,
                                     const std::string &i_optRead ) {
#ifdef PP_USE_MOAB
  m_moab.read( i_pathToMesh,
               i_optRead );
#else
  assert( false );
#endif

  printStats();
}

void edge::mesh::Unstructured::write( const std::string &i_pathToMesh ) {
#ifdef PP_USE_MOAB
  m_moab.write( i_pathToMesh );
#else
  assert( false );
#endif
}

t_enLayout edge::mesh::Unstructured::getVeLayout() const {
#ifdef PP_USE_MOAB
  return m_moab.getVeLayout();
#else
  assert(false);
#endif
}

t_enLayout edge::mesh::Unstructured::getFaLayout() const {
#ifdef PP_USE_MOAB
  return m_moab.getFaLayout();
#else
  assert(false);
#endif
}

t_enLayout edge::mesh::Unstructured::getElLayout() const {
#ifdef PP_USE_MOAB
  return m_moab.getElLayout();
#else
  assert(false);
#endif
}

const t_inMap* edge::mesh::Unstructured::getInMap() const {
#ifdef PP_USE_MOAB
  return m_moab.getInMap();
#else
  assert(false);
#endif
}

int_el edge::mesh::Unstructured::getNVertices() const {
#ifdef PP_USE_MOAB
  return m_moab.getNVertices();
#else
  assert(false);
#endif
}

int_el edge::mesh::Unstructured::getNFaces() const {
#ifdef PP_USE_MOAB
  return m_moab.getNFaces();
#else
  assert(false);
#endif
}

int_el edge::mesh::Unstructured::getNElements() const {
#ifdef PP_USE_MOAB
  return m_moab.getNElements();
#else
  assert(false);
#endif
}


int_el edge::mesh::Unstructured::getNelVeEl() {
#ifdef PP_USE_MOAB
  return m_moab.getNelVeEl();
#else
  assert(false);
#endif
}

void edge::mesh::Unstructured::getElVe( int_el (*o_elVe)[C_ENT[T_SDISC.ELEMENT].N_VERTICES] ) {
#ifdef PP_USE_MOAB
  m_moab.getElVe( o_elVe );
#else
  assert(false);
#endif
}

void edge::mesh::Unstructured::getElementsAdjacentFaces( int_el (*o_elementsAdjacentFaces)[C_ENT[T_SDISC.ELEMENT].N_FACES] ) {
#ifdef PP_USE_MOAB
  m_moab.getElementsAdjacentFaces( o_elementsAdjacentFaces );
#else
  assert(false);
#endif
}

void edge::mesh::Unstructured::getFacesAdjacentElements( int_el (*o_neighboringIds)[2] ) {
#ifdef PP_USE_MOAB
  m_moab.getFacesAdjacentElements( o_neighboringIds );
#else
  assert(false);
#endif
}

void edge::mesh::Unstructured::getElVeEl( int_el **o_elVeEl ) {
#ifdef PP_USE_MOAB
  m_moab.getElVeEl( o_elVeEl );
#else
  assert(false);
#endif
}

void edge::mesh::Unstructured::getFacesAdjacentVertices( int_el (*o_neighVeIds)[C_ENT[T_SDISC.FACE].N_VERTICES] ) {
#ifdef PP_USE_MOAB
  m_moab.getFacesAdjacentVertices( o_neighVeIds );
#else
  assert(false);
#endif
}

void edge::mesh::Unstructured::getElementsFaceNeighbors( int_el (*o_neighboringIds)[C_ENT[T_SDISC.ELEMENT].N_FACES] ) {
#ifdef PP_USE_MOAB
  m_moab.getElementsFaceNeighbors( o_neighboringIds );
#else
  assert(false);
#endif
}

void edge::mesh::Unstructured::getGIdsVe( std::vector< int_gid > &o_gIds ) const {
#ifdef PP_USE_MOAB
  m_moab.getGIdsVe( o_gIds );
#else
  assert(false);
  return NULL;
#endif
}

void edge::mesh::Unstructured::getGIdsFa( std::vector< int_gid > &o_gIds ) const {
#ifdef PP_USE_MOAB
  m_moab.getGIdsFa( o_gIds );
#else
  assert(false);
  return NULL;
#endif
}

void edge::mesh::Unstructured::getGIdsEl( std::vector< int_gid > &o_gIds ) const {
#ifdef PP_USE_MOAB
  m_moab.getGIdsEl( o_gIds );
#else
  assert(false);
  return NULL;
#endif
}

void edge::mesh::Unstructured::getConnect( const t_vertexChars *i_veChars,
                                           const t_faceChars   *i_faChars,
                                                 t_connect     &o_connect ) {
  // get the connectivity info
  getElVe(                     o_connect.elVe   );
  getElementsAdjacentFaces(    o_connect.elFa   );
  getFacesAdjacentElements(    o_connect.faEl   );
  getFacesAdjacentVertices(    o_connect.faVe   );
  //getElVeEl(                   o_connect.elVeEl );
  getElementsFaceNeighbors(    o_connect.elFaEl );

  // get the element and face layout
  t_enLayout l_faLayout, l_elLayout;
  l_faLayout = getFaLayout();
  l_elLayout = getElLayout();

  // get the index mappings
  t_inMap const * l_inMap = getInMap();

  std::vector< int_gid > l_gIdsVe, l_gIdsFa, l_gIdsEl;
  getGIdsVe( l_gIdsVe ); getGIdsFa( l_gIdsFa ); getGIdsEl( l_gIdsEl );

  // check the constitency of the performed under the hood ordering
  // otherwise the upcoming mapping to the reference element is invalid since it assumes an ascending ordering
  // of elVe and elFa + elFaEL w.r.t. to the vertices of the faces; e.g. for the vertex faces
  // [1, 4, 10] < [2, 4, 10]; [1, 4, 10] < [1, 5, 10]; [1, 4, 10] < [1, 4, 11]
  // [1, 4, 10] > [0, 4, 10]; [1, 4, 10] > [1, 3, 10]; [1, 4, 10] > [1, 4, 9]
  common< T_SDISC.ELEMENT >::checkConsOrd( l_faLayout,
                                           l_elLayout,
                                           l_gIdsVe,
                                           l_gIdsFa,
                                           l_gIdsEl,
                                           o_connect.elVe,
                                           o_connect.elFa,
                                           o_connect.faEl,
                                           o_connect.faVe
#ifdef PP_USE_MPI
                                          ,false
#endif
                                         );

#if defined PP_T_ELEMENTS_TRIA3
  // derive a normalized mapping to the reference element by ensuring counter-clockwise ordering of the vertices and faces
  common< T_SDISC.ELEMENT >::normOrdTria( l_elLayout,
                                          i_veChars,
                                          o_connect.elVe,
                                          o_connect.elFa,
                                          o_connect.elFaEl );
#elif defined PP_T_ELEMENTS_TET4
  // derive a normalized mapping to the reference element by ensuring counter-clockwise ordering of vertices 2-3 based on 0 and 1.
  common< T_SDISC.ELEMENT >::normOrdTet4( l_elLayout,
                                          i_veChars,
                                          o_connect.elVe,
                                          o_connect.elFa,
                                          o_connect.elFaEl );
#else
  assert( false );
#endif

  // check the consistency of the adjacency info
  common< T_SDISC.ELEMENT >::checkConsAdj( l_faLayout,
                                           i_faChars,
                                           l_inMap->faMeDa,
                                           l_inMap->faDaMe,
                                           o_connect.elFaEl,
                                           o_connect.elFa,
                                           o_connect.faEl );

  // get the neighboring elements' local face ids
  common< T_SDISC.ELEMENT >::getFIdsElFaEl( l_elLayout,
                                            l_gIdsEl,
                                            o_connect.elFa,
                                            o_connect.elFaEl,
                                            o_connect.fIdElFaEl );

  // get the neighboring elements' local vertex ids, which match the faces' first vertex
  common< T_SDISC.ELEMENT >::getVIdsElFaEl( l_elLayout,
                                            l_gIdsVe,
                                            o_connect.elFaEl,
                                            o_connect.elVe,
                                            o_connect.fIdElFaEl,
                                            o_connect.vIdElFaEl,
#ifndef PP_USE_MPI
                                            true,
#else
                                            false,
#endif
                                            i_veChars  );
}

void edge::mesh::Unstructured::getVeChars( t_vertexChars *o_veChars ) {
  for( int_el l_ve = 0; l_ve < getNVertices(); l_ve++ ) {
#ifdef PP_USE_MOAB
    m_moab.getVeCoords( l_ve,
                        o_veChars[l_ve].coords[0],
                        o_veChars[l_ve].coords[1],
                        o_veChars[l_ve].coords[2] );

    // get vertex type
    m_moab.getMatVal( 0, l_ve, o_veChars[l_ve].spType );
#else
    assert(false);
#endif
  }
}

void edge::mesh::Unstructured::getFaChars( t_faceChars *o_faChars ) {
  for( int_el l_fa = 0; l_fa < getNFaces(); l_fa++ ) {
    double l_faceVertices[3][C_ENT[T_SDISC.FACE].N_VERTICES];

    // coordinates of point in the left element, corresponding to 0 in the face's adjacency
    double l_normalPoint[3];

#ifdef PP_USE_MOAB
    m_moab.getFaceVerticesCoords( l_fa, l_faceVertices[0], l_faceVertices[1], l_faceVertices[2] );
    bool l_left;
    m_moab.getNormalPointCoords( l_fa, l_left, l_normalPoint[0], l_normalPoint[1], l_normalPoint[2] );

    // compute face characteristics
    o_faChars[l_fa].area = linalg::Geom::volume( T_SDISC.FACE, l_faceVertices[0], N_DIM );
    linalg::Geom::computeOutPtNormal( T_SDISC.FACE,
                                      l_faceVertices[0],
                                      l_normalPoint,
                                      o_faChars[l_fa].outNormal );
    linalg::Geom::computeTangents( T_SDISC.FACE,
                                   l_faceVertices[0],
                                   o_faChars[l_fa].tangent0,
                                   o_faChars[l_fa].tangent1 );

    // change orientation of outer-pointing normal if we derived it with a point on the right side.
    if( l_left == false ) {
      for( unsigned int l_dim = 0; l_dim < 3; l_dim ++ ) {
        o_faChars[l_fa].outNormal[l_dim] = -o_faChars[l_fa].outNormal[l_dim];
      }
    }

    // get face type
    m_moab.getMatVal( N_DIM-1, l_fa, o_faChars[l_fa].spType );
#else
    assert(false);
#endif
  }
}

void edge::mesh::Unstructured::getElChars( t_elementChars *o_elChars ) {
  for( int_el l_el = 0; l_el < getNElements(); l_el++ ) {
    real_mesh l_veCrds[3][C_ENT[T_SDISC.ELEMENT].N_VERTICES];
#ifdef PP_USE_MOAB
    m_moab.getElementVerticesCoords( l_el,
                                     l_veCrds[0],
                                     l_veCrds[1],
                                     l_veCrds[2] );
#else
    EDGE_LOG_FATAL;
#endif

    // compute element characteristics
    o_elChars[l_el].volume = linalg::Geom::volume( T_SDISC.ELEMENT, l_veCrds[0] );
    o_elChars[l_el].inDia  = linalg::Geom::inDia(  T_SDISC.ELEMENT, l_veCrds[0] );

#if defined PP_T_ELEMENTS_TRIA3
    // check the that volume of incircle is smaller than triangle volume
    EDGE_CHECK_LT( o_elChars[l_el].inDia * o_elChars[l_el].inDia * M_PI/4.0, o_elChars[l_el].volume );
#elif defined PP_T_ELEMENTS_TET4
    // check that the volume of the insphere is smaller than the tretrahedral volume
    EDGE_CHECK_LT( o_elChars[l_el].inDia * o_elChars[l_el].inDia * o_elChars[l_el].inDia * M_PI/6.0,
                   o_elChars[l_el].volume );
#endif

    // get element type
    m_moab.getMatVal( N_DIM, l_el, o_elChars[l_el].spType );
  }
}
