/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2016-2018, Regents of the University of California
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
 * Interface to MOAB.
 **/

#include "Moab.h"
#include <cassert>
#include <cmath>
#ifdef PP_USE_MPI
#include <MBParallelConventions.h>
#endif
#include <MBTagConventions.hpp>

edge::mesh::Moab::Moab(       unsigned int  i_dim,
                              unsigned int  i_nBndVals,
                        const int          *i_bndVals,
                              int           i_periodicVal ):
#ifdef PP_USE_MPI
  m_core(),
  m_pcomm( &m_core, MPI_COMM_WORLD ),
#endif
 m_dim(i_dim),
 m_nBndVals(i_nBndVals),
 m_bndVals(i_bndVals)
#ifndef PP_USE_MPI
 ,m_periodicVal(i_periodicVal)
#endif
  {
  m_tagMat = 0;
  m_tagLId = 0;
  m_tagGId = 0;

  EDGE_LOG_INFO << "  initialized MOAB:";
  printMoab();
}

void edge::mesh::Moab::printMoab() {
  // query moab info
  std::string l_apiVersion, l_implVersion;
  m_core.api_version( &l_apiVersion ); 
  m_core.impl_version( &l_implVersion );

  EDGE_LOG_INFO << "    " << l_apiVersion;
  EDGE_LOG_INFO << "    " << l_implVersion;
}

void edge::mesh::Moab::initEn( const std::vector< moab::EntityHandle > &i_ents,
                                     std::vector< int_el >             &o_enMeDa,
                                     std::vector< int_el >             &o_enDaMe ) {
  moab::ErrorCode l_error;

  // id in the mesh
  int_el l_meshId = 0;

  for( std::size_t l_en = 0; l_en < i_ents.size(); l_en++ ) {
    moab::EntityHandle l_enHan = i_ents[l_en];

    // get the current id
    int_el l_curId;
    l_error = m_core.tag_get_data( m_tagLId, &l_enHan, 1, &l_curId ); EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );

    // we haven't touched this
    if( l_curId == m_defaultTagLId ) {
      l_error = m_core.tag_set_data( m_tagLId, &l_enHan, 1, &l_en ); EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );

      o_enDaMe.push_back( l_meshId );
      o_enMeDa.push_back( l_en );
      l_meshId++;
    }
    // this is duplicated data
    else {
      EDGE_CHECK_LT( l_curId, (int_el) o_enDaMe.size() );
      o_enDaMe.push_back( o_enDaMe[l_curId] );
    }
  }
}

void edge::mesh::Moab::init() {
  moab::ErrorCode l_error;

  // create local id tag of entities
  l_error = m_core.tag_get_handle(  "local_id",
                                    1,
                                    moab::MB_TYPE_INTEGER,
                                    m_tagLId,
                                    moab::MB_TAG_CREAT|moab::MB_TAG_EXCL|moab::MB_TAG_DENSE,
                                   &m_defaultTagLId );
  EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );

  // get the size of the tag and check with internal definition
  int l_lIdSize = 0;
  l_error = m_core.tag_get_bytes( m_tagLId, l_lIdSize ); EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );
  CHECK( l_lIdSize == sizeof( int_el ) );

  initEn( m_vertices, m_inMap.veMeDa, m_inMap.veDaMe );
  initEn( m_faces,    m_inMap.faMeDa, m_inMap.faDaMe );
  initEn( m_elements, m_inMap.elMeDa, m_inMap.elDaMe );
}

moab::EntityHandle edge::mesh::Moab::getEnHandle( unsigned short i_nDim,
                                                  int_el         i_en ) {
  // get entity handle
  moab::EntityHandle l_ha = std::numeric_limits< moab::EntityHandle >::max();
  if( i_nDim == 0 )
    l_ha = m_vertices[i_en];
  else if( i_nDim == m_dim-1 )
    l_ha = m_faces[i_en];
  else if( i_nDim == m_dim )
    l_ha = m_elements[i_en];
  else EDGE_LOG_FATAL << "invalid dimension: " << i_nDim << " " << i_en;

  return l_ha;
}

bool edge::mesh::Moab::isBnd( moab::EntityHandle i_ent ) {
  bool l_isBnd = false;

  // iterate over boundary meshsets and check for the entity
  for( std::size_t l_set = 0; l_set < m_bndSets.size(); l_set++ ) {
    l_isBnd = m_core.contains_entities(  m_bndSets[l_set],
                                        &i_ent,
                                         1 );
    if( l_isBnd == true ) break;
  }

#ifdef PP_USE_MPI
  // ghost entities are at MPI-"boundaries"
  if( getEnMpiType(i_ent) == 2 ) l_isBnd = true;
#endif

  return l_isBnd;
}

void edge::mesh::Moab::setPeriodicAdjacenciesDim( const std::vector< int_el >  &i_faces,
                                                  const moab::Range            &i_handles,
                                                  const real_mesh             (*i_freeDims)[2] ) {
  moab::ErrorCode l_error;

  for( std::size_t l_fa = 0; l_fa < i_faces.size(); l_fa++ ) {
    // get local array id
    int_el l_faId = i_faces[l_fa];

    // find the partner
    for( std::size_t l_faP = l_fa+1; l_faP < i_faces.size(); l_faP++ ) {
      // get partner array id
      int_el l_faPId = i_faces[l_faP];

      // check if we have a match
      if( std::abs( i_freeDims[l_faId][0] - i_freeDims[l_faPId][0] ) < TOL.MESH &&
          std::abs( i_freeDims[l_faId][1] - i_freeDims[l_faPId][1] ) < TOL.MESH ) {
        // get the elements on either side
        moab::Range        l_element;
        moab::EntityHandle l_faHandle = i_handles[l_faId];
       
        l_error = m_core.get_adjacencies( &l_faHandle,
                                           1,
                                           m_dim,
                                           true,
                                           l_element,
                                           moab::Interface::UNION );
        EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );
        assert( l_element.size() == 1 );

        moab::Range        l_elementP;
        moab::EntityHandle l_faHandleP = i_handles[l_faPId];
       
        l_error = m_core.get_adjacencies( &l_faHandleP,
                                           1,
                                           m_dim,
                                           true,
                                           l_elementP,
                                           moab::Interface::UNION );
        EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );
        assert( l_elementP.size() == 1 );

        // store the adjacencies between the two elements in a custom lookup table.. urgh.. finally!
        m_periodicFaces.push_back(          l_faHandle    );
        m_periodicFaces.push_back(          l_faHandleP   );

        m_periodicElementsRemote.push_back( l_elementP[0] );
        m_periodicElementsRemote.push_back( l_element[0]  );

        // remove one face from faces
        m_faces.erase( std::remove( m_faces.begin(), m_faces.end(), l_faHandleP ), m_faces.end() );

        // TODO: Needs to be fixed for LTS, this is GTS only
        m_faLayout.nEnts--;
        m_faLayout.timeGroups[0].inner.size--;
        m_faLayout.timeGroups[0].nEntsOwn--;
      }
    }
  }
}

void edge::mesh::Moab::setPeriodicAdjacencies( int i_tagValue ) {
  moab::ErrorCode l_error;

  // get the entity sets with material sets defined
  moab::Range l_sets;
  l_error = m_core.get_entities_by_type_and_tag(  0, moab::MBENTITYSET, &m_tagMat, NULL, 1, l_sets );
  EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );
  assert( l_sets.size() > 0 );

  // find the periodic boundary set
  moab::EntityHandle l_periodicSet = std::numeric_limits< moab::EntityHandle >::max();
  for( std::size_t l_set = 0; l_set < l_sets.size(); l_set++ ) {
    // get value of material set w.r.t. to set
    int l_matVal;
    l_periodicSet = l_sets[l_set];
    l_error = m_core.tag_get_data(  m_tagMat,
                                   &l_periodicSet,
                                    1,
                                   &l_matVal );
    EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );

    // success for periodic value
    if( l_matVal == i_tagValue ) break;

    // inform user if we don't have a periodic boundary
    if( l_set == l_sets.size()-1 ) {
      EDGE_LOG_INFO << "  couldn't find a periodic mesh region " << i_tagValue << ", continuing w/o";
      return;
    }
  }
  EDGE_CHECK_NE( l_periodicSet, std::numeric_limits< moab::EntityHandle >::max() );

  // get boundary faces
  moab::Range l_bndFaces;
  l_error = m_core.get_entities_by_dimension( l_periodicSet, m_dim-1, l_bndFaces );
  EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );

  // assert every face has a partner
  assert( l_bndFaces.size() % 2 == 0 );

  // iterate over the entities and get the mid points per dimension

  // indices of x-, y- and z- bnd faces
  std::vector< int_el > l_bndXId, l_bndYId, l_bndZId;

  // positions of bnd-midpoints w.r.t. to remaining free dimensions
  real_mesh (*l_bndPos)[2] =  (real_mesh(*)[2]) new real_mesh[2*l_bndFaces.size()];

  for( std::size_t l_fa = 0; l_fa < l_bndFaces.size(); l_fa++ ) {
    // init with invalid values
    l_bndPos[l_fa][0] = l_bndPos[l_fa][1] = std::numeric_limits< real_mesh >::max();

    // get vertices of the face
    moab::EntityHandle l_face = l_bndFaces[l_fa];
    moab::Range l_vertices;
    l_error = m_core.get_adjacencies( &l_face, 1, 0, true, l_vertices, moab::Interface::UNION );
    EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );
    // assert we have a enough vertices for a face
    assert( l_vertices.size() > m_dim-1 );

    // reset mid points
    real_mesh l_midX = 0; real_mesh l_midY = 0; real_mesh l_midZ = 0;

    // coordinates of the vertices
    std::vector< real_mesh > l_veCoords[3];

    // iterate over vertices and compute midpoints dimension-wise
    for( std::size_t l_ve = 0; l_ve < l_vertices.size(); l_ve++ ) {
      real_mesh l_coords[3];
      moab::EntityHandle l_vertex = l_vertices[l_ve];
      // get coordinates of vertex
      l_error = m_core.get_coords( &l_vertex, 1, l_coords ); EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );

      // add value to face-mid points
      l_midX += l_coords[0]; l_midY += l_coords[1]; l_midZ += l_coords[2];

      // store vertex coords
      l_veCoords[0].push_back( l_coords[0] );
      l_veCoords[1].push_back( l_coords[1] );
      l_veCoords[2].push_back( l_coords[2] );

      // scale decide which boundary we are checking mid point against last vertex
      if( l_ve == l_vertices.size()-1 ) {
        l_midX /= l_vertices.size(); l_midY /= l_vertices.size(); l_midZ /= l_vertices.size(); 

        real_mesh l_vars[3] = {0, 0, 0};

        // compute variations of vertex coords
        for( unsigned short l_veVar = 1; l_veVar < l_vertices.size(); l_veVar++ ) {
          l_vars[0] += std::abs( l_veCoords[0][l_veVar] - l_veCoords[0][l_veVar-1] );
          l_vars[1] += std::abs( l_veCoords[1][l_veVar] - l_veCoords[1][l_veVar-1] );
          l_vars[2] += std::abs( l_veCoords[2][l_veVar] - l_veCoords[2][l_veVar-1] );
        }

        // determine bnd-dimension
        if(      l_vars[0] < TOL.MESH ) {
          l_bndPos[l_fa][0] = l_midY; l_bndPos[l_fa][1] = l_midZ;
          l_bndXId.push_back( l_fa );
        }
        else if( l_vars[1] < TOL.MESH ) {
          l_bndPos[l_fa][0] = l_midX; l_bndPos[l_fa][1] = l_midZ;
          l_bndYId.push_back( l_fa );
        }
        else if( l_vars[2] < TOL.MESH ) {
          l_bndPos[l_fa][0] = l_midX; l_bndPos[l_fa][1] = l_midY;
          l_bndZId.push_back( l_fa );
        }
        else assert(false);
      }
    }

    // check for valid entries
    EDGE_CHECK( l_bndPos[l_fa][0] != std::numeric_limits< real_mesh >::max() );
    EDGE_CHECK( l_bndPos[l_fa][1] != std::numeric_limits< real_mesh >::max() );
  }

  EDGE_LOG_INFO << "  #periodic faces:";
  EDGE_LOG_INFO << "    x: " << l_bndXId.size();
  EDGE_LOG_INFO << "    y: " << l_bndYId.size();
  EDGE_LOG_INFO << "    z: " << l_bndZId.size();

  // ensure every face has a partner
  assert( l_bndXId.size() % 2 == 0 );
  assert( l_bndYId.size() % 2 == 0 );
  assert( l_bndZId.size() % 2 == 0 );

  setPeriodicAdjacenciesDim( l_bndXId, l_bndFaces, l_bndPos );
  setPeriodicAdjacenciesDim( l_bndYId, l_bndFaces, l_bndPos );
  setPeriodicAdjacenciesDim( l_bndZId, l_bndFaces, l_bndPos );

  delete[] l_bndPos;
}

void edge::mesh::Moab::fixPeriodicTagIds() {
  moab::ErrorCode l_error;

  for( std::size_t l_fa = 0; l_fa < m_periodicFaces.size(); l_fa += 2 ) {
    // get id of face still under control
    int_el l_validId;
    moab::EntityHandle l_validFace = m_periodicFaces[l_fa];
    l_error = m_core.tag_get_data(  m_tagLId,
                                   &l_validFace,
                                    1,
                                   &l_validId );
    EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );

    // copy over
    moab::EntityHandle l_invalidFace = m_periodicFaces[l_fa+1];
    l_error = m_core.tag_set_data(  m_tagLId,
                                   &l_invalidFace,
                                    1,
                                   &l_validId );
    EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );
  }
}

#ifdef PP_USE_MPI
unsigned short edge::mesh::Moab::getEnMpiType( moab::EntityHandle i_en,
                                               int                i_rank ) {
  EDGE_CHECK_NE( i_rank, parallel::g_rank );

  moab::ErrorCode l_err;

  // get the parallel status of the entity
  unsigned char l_ps;
  l_err = m_pcomm.get_pstatus( i_en, l_ps ); EDGE_CHECK_EQ( l_err, moab::MB_SUCCESS );

  // check if the entity is owned
  bool l_own = !(l_ps & PSTATUS_NOT_OWNED);

  // no other flags set means inner-entity
  if( l_own && !(l_ps & (PSTATUS_SHARED | PSTATUS_MULTISHARED | PSTATUS_INTERFACE | PSTATUS_GHOST)) ) {
    // double check that we didn't miss anything
    EDGE_CHECK_EQ( l_ps, 0 );
    return 0;
  }

  /*
   * Check if the entity MPI-status is for real:
   *   MOAB requires us to consider vertices as bridge dim.
   *   Faces as bridge-dim is what we are interested in.
   * -> The entity (if a face itself) or one adjacent face of the entity
   *    has an adjacent local-rank and neighboring-rank element.
   */

  // get the dimension of the entity
  unsigned short l_dim = m_core.dimension_from_handle( i_en );

  // derive adjacent faces
  std::vector< moab::EntityHandle > l_enFa;

  // vertex: most complex since all vertices belonging to ghost elements are valid ghost vertices
  if( l_dim == 0 ) {
    // deterime adjacent elements
    std::vector< moab::EntityHandle > l_veEl;
    l_err = m_core.get_adjacencies( &i_en, 1, m_dim, true, l_veEl ); EDGE_CHECK_EQ( l_err, moab::MB_SUCCESS );
    EDGE_CHECK_GE( l_veEl.size(), 1 );

    // determine faces of these elements
    l_err = m_core.get_adjacencies( &l_veEl[0], l_veEl.size(), m_dim-1, true, l_enFa, moab::Core::UNION );
    EDGE_CHECK_EQ( l_err, moab::MB_SUCCESS );
  }
  // faces
  else if( l_dim == m_dim-1 ) {
    l_enFa.push_back( i_en );
  }
  // element
  else {
    // it has to be a vertex, face or element
    EDGE_CHECK_EQ( l_dim, m_dim );

    l_err = m_core.get_adjacencies( &i_en, 1, m_dim-1, true, l_enFa ); EDGE_CHECK_EQ( l_err, moab::MB_SUCCESS );
  }

  // check that we found an adjacent face
  EDGE_CHECK_GT( l_enFa.size(), 0 ) << "couldn't find adjacent faces " << l_dim;


  // MPI-status of the entity.
  // [0] == true: one of the adjacent faces has two adjacent elements owned by the current rank; bnd-conditions
  //              count as "adjacent element"
  // [1] == true: one of the adjacent faces has one adjacent element owned by the current rank and
  //              one owned by a valid (see input of the function) neighboring rank.
  // [2] == true: one of the adjacent faces has one adjacent element owned by an invalid neighboring rank.
  bool l_mpiStatus[3] = { false, false, false };

  // iterate over the faces
  for( unsigned short l_fa = 0; l_fa < l_enFa.size(); l_fa++ ) {
    // get adjacent elements
    std::vector< moab::EntityHandle > l_enFaEl;
    l_err = m_core.get_adjacencies( &l_enFa[l_fa], 1, m_dim, true, l_enFaEl ); EDGE_CHECK_EQ( l_err, moab::MB_SUCCESS );

    // check that there's at least one adjacent element
    EDGE_CHECK_GE( l_enFaEl.size(), 1 );

    // faces with one adjacent element might describe bnd-conditions and could be inner faces
    if( l_enFaEl.size() == 1 ) {
      // get the parallel status of the one element
      unsigned char l_psAdj;
      l_err = m_pcomm.get_pstatus( l_enFaEl[0], l_psAdj ); EDGE_CHECK_EQ( l_err, moab::MB_SUCCESS );
      // non-relevant if this is ghost element
      if( l_psAdj & PSTATUS_NOT_OWNED ) {
        l_mpiStatus[2] = true;
      }
      // this must be a bnd-face, obviously no ghost elements present
      else l_mpiStatus[0] = true;
    }
    else {
      // get the ranks owning the two elements
      int l_neRanks[2];
      l_err = m_pcomm.get_owner( l_enFaEl[0], l_neRanks[0] ); EDGE_CHECK_EQ( l_err, moab::MB_SUCCESS );
      l_err = m_pcomm.get_owner( l_enFaEl[1], l_neRanks[1] ); EDGE_CHECK_EQ( l_err, moab::MB_SUCCESS );

      // [0]: #local ranks, [2]: #valid remote ranks
      int l_count[2] = {0,0};
      for( unsigned short l_ne = 0; l_ne < 2; l_ne++ ) {
        if( l_neRanks[l_ne] == parallel::g_rank )
          l_count[0]++;
        else if( l_neRanks[l_ne] == i_rank || i_rank == -1 )
          l_count[1]++;
      }

      if(      l_count[0] == 2                    ) l_mpiStatus[0] = true;
      else if( l_count[0] == 1 && l_count[1] == 1 ) l_mpiStatus[1] = true;
      else                                          l_mpiStatus[2] = true;
    }
  }
  // abort if this is not a valid MPI-entity nor an inner-entity
  if( !l_mpiStatus[0] && !l_mpiStatus[1] ) return std::numeric_limits< unsigned short >::max();

  // this is a valid mpi-entity
  else if( l_mpiStatus[1] ) {
    if(  l_ps & PSTATUS_NOT_OWNED ) return 2;
    else                            return 1;
  }
  else if( l_mpiStatus[2] ) {
    return std::numeric_limits< unsigned short >::max();
  }
  // else: inner-entity
  EDGE_CHECK( l_mpiStatus[0] );
  EDGE_CHECK( !(l_ps & PSTATUS_NOT_OWNED) );

  return 0;
}
#endif

void edge::mesh::Moab::setupEnLayout( int_tg                             i_tg,
                                      moab::Range                       &i_ents,
                                      t_enLayout                        &o_enLayout,
                                      std::vector< moab::EntityHandle > &o_ents ) {
#ifdef PP_USE_MPI
  moab::ErrorCode l_error;

  // allocate total ent size with 5% overhead
  o_ents.reserve( i_ents.size() * 1.05 );

  // get inner-elements
  for( moab::Range::const_iterator l_en = i_ents.begin(); l_en != i_ents.end(); l_en++ ) {
    if( getEnMpiType(*l_en) == 0 ) o_ents.push_back( *l_en );
  }
  EDGE_CHECK_GT( o_ents.size(), 0 ); // make sure we have inner entities

  // store the inner info
  o_enLayout.timeGroups[i_tg].inner.first = o_enLayout.nEnts;
  o_enLayout.timeGroups[i_tg].inner.size  = o_ents.size();
  o_enLayout.nEnts                       += o_ents.size();
  o_enLayout.timeGroups[i_tg].nEntsOwn    = o_ents.size();

  // get neighboring ranks (sets are sorted)
  std::set< unsigned int > l_neRanks;
  l_error = m_pcomm.get_comm_procs( l_neRanks ); EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );

  // parallel status for everything MPI-releated
  unsigned char l_psMpi = PSTATUS_GHOST | PSTATUS_SHARED | PSTATUS_MULTISHARED| PSTATUS_INTERFACE;

  // iterate over the neighoring ranks and determine true neighbors:
  // since we are deriving face-adjacency through vertex adjacency, certain
  // ranks might only share elements through vertices which would lead to empty messages..
  for( std::set< unsigned int >::iterator l_nr = l_neRanks.begin(); l_nr != l_neRanks.end(); ) {
    // get the entities shared with this processor
    moab::Range l_enMpi;
    l_error = m_pcomm.filter_pstatus( i_ents, l_psMpi, PSTATUS_OR, *l_nr, &l_enMpi ); EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );

    bool l_realMpi = false;
    // iterate over the shared ents and check for send-/receive-ents
    for( std::size_t l_en = 0; l_en < l_enMpi.size(); l_en++ ) {
      unsigned short l_mpiType = getEnMpiType( l_enMpi[l_en], *l_nr );
      // abort in the case of send-/recv-ent
      if( l_mpiType == 1 || l_mpiType == 2 ) {
        l_realMpi = true;
        break;
      }
    }

    // erase rank if pseudo (vertex neigbhbors only)
    if( !l_realMpi ) l_nr = l_neRanks.erase( l_nr );
    else l_nr++;
  }

  // resize containers of comm. data
  o_enLayout.timeGroups[i_tg].send.resize(    l_neRanks.size() );
  o_enLayout.timeGroups[i_tg].receive.resize( l_neRanks.size() );
  o_enLayout.timeGroups[i_tg].neRanks.resize( l_neRanks.size() );
  o_enLayout.timeGroups[i_tg].neTgs.resize(   l_neRanks.size() );

  // neighbring id
  unsigned int l_nId = 0;

  // iterate over send regions
  for( std::set< unsigned int >::iterator l_nr = l_neRanks.begin(); l_nr != l_neRanks.end(); l_nr++ ) {
    // get the entities shared with this processor
    moab::Range l_enMpi;
    l_error = m_pcomm.filter_pstatus( i_ents, l_psMpi, PSTATUS_OR, *l_nr, &l_enMpi ); EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );

    // derive send entities
    std::vector< moab::EntityHandle > l_enSe;
    for( moab::Range::const_iterator l_en = l_enMpi.begin(); l_en != l_enMpi.end(); l_en++ ) {
      if( getEnMpiType(*l_en, *l_nr) == 1 ) l_enSe.push_back( *l_en );
    }

    // sort the send entities by their global id
    sortGId( l_enSe );

    // add send entities
    for( std::size_t l_en = 0; l_en < l_enSe.size(); l_en++ ) o_ents.push_back( l_enSe[l_en] );

    // store the send info
    o_enLayout.timeGroups[i_tg].send[l_nId].first = o_enLayout.nEnts;
    o_enLayout.timeGroups[i_tg].send[l_nId].size  = l_enSe.size();
    o_enLayout.nEnts                             += l_enSe.size();
    o_enLayout.timeGroups[i_tg].nEntsOwn         += l_enSe.size();
    o_enLayout.timeGroups[i_tg].neRanks[l_nId]    = *l_nr;

    // set neighboring time group
    EDGE_CHECK( i_tg == 0 ); // TODO: Add LTS support for MPI
    o_enLayout.timeGroups[i_tg].neTgs[l_nId] = 0;

    l_nId++;
  }

  // reset neighboring id
  l_nId = 0;
  o_enLayout.timeGroups[i_tg].nEntsNotOwn = 0;

  // iterate over receive regions
  for( std::set< unsigned int >::iterator l_nr = l_neRanks.begin(); l_nr != l_neRanks.end(); l_nr++) {
    // get the entities shared with this processor
    moab::Range l_enMpi;
    l_error = m_pcomm.filter_pstatus( i_ents, l_psMpi, PSTATUS_OR, *l_nr, &l_enMpi ); EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );

    // derive receive entities
    std::vector< moab::EntityHandle > l_enRe;
    for( moab::Range::const_iterator l_en = l_enMpi.begin(); l_en != l_enMpi.end(); l_en++ ) {
      if( getEnMpiType(*l_en, *l_nr) == 2 ) l_enRe.push_back( *l_en );
    }

    // sort the receive entities by their global id
    sortGId( l_enRe );

    // add receive entities
    for( std::size_t l_en = 0; l_en < l_enRe.size(); l_en++ ) o_ents.push_back( l_enRe[l_en] );

    // store the receive info
    o_enLayout.timeGroups[i_tg].receive[l_nId].first = o_enLayout.nEnts;
    o_enLayout.timeGroups[i_tg].receive[l_nId].size  = l_enRe.size();
    o_enLayout.nEnts                                += l_enRe.size();
    o_enLayout.timeGroups[i_tg].nEntsNotOwn         += l_enRe.size();

    l_nId++;
  }
#else
  // non-MPI layout
  for( std::size_t l_en = 0; l_en < i_ents.size(); l_en++ ) o_ents.push_back( i_ents[l_en] );

  o_enLayout.timeGroups[i_tg].inner.first = o_enLayout.nEnts;
  o_enLayout.timeGroups[i_tg].inner.size  = i_ents.size();
  o_enLayout.timeGroups[i_tg].nEntsOwn    = i_ents.size();
  o_enLayout.nEnts += i_ents.size();
#endif
}

void edge::mesh::Moab::setupDataLayout() {
  // setup a single time stepping group for now
  unsigned int l_nTimeGroups = 1;
  m_elLayout.nEnts = m_faLayout.nEnts = m_veLayout.nEnts = 0;
  m_elLayout.timeGroups.resize(l_nTimeGroups);
  m_faLayout.timeGroups.resize(l_nTimeGroups);
  m_veLayout.timeGroups.resize(l_nTimeGroups);

  moab::ErrorCode l_error;

  // get ghost entities
#ifdef PP_USE_MPI
  // higher-dimensional (>0) bridges in MOAB do not work in parallel; use vertices and derive face-adjacency manually
  l_error = m_pcomm.exchange_ghost_cells( m_dim, 0, 1, m_dim-1, true, true ); EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );
#endif

  // elements, faces and vertices
  moab::Range l_elements, l_faces, l_vertices;

  // get elements
  l_error = m_core.get_entities_by_dimension( 0, m_dim, l_elements ); EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );
  EDGE_CHECK( l_elements.size() != 0 );

  // get vertices
  l_error = m_core.get_entities_by_dimension( 0, 0, l_vertices ); EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );
  EDGE_CHECK( l_vertices.size() != 0 );

  // create & get faces
  l_error = m_core.get_adjacencies( l_elements, m_dim-1, true, l_faces, moab::Interface::UNION );
  EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );
  EDGE_CHECK_NE( l_faces.size(), 0 );

  // create global id tag of entities
  l_error = m_core.tag_get_handle(  "global_id",
                                    1,
                                    moab::MB_TYPE_INTEGER,
                                    m_tagGId,
                                    moab::MB_TAG_CREAT|moab::MB_TAG_EXCL|moab::MB_TAG_DENSE,
                                   &m_defaultTagGId );
  EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );

  // check that our gid fits into int
  CHECK_EQ( sizeof(int), sizeof( int_gid ) );

  // get MOAB's global id tag of entities
  moab::Tag l_tagMoabGId;
  l_error = m_core.tag_get_handle( GLOBAL_ID_TAG_NAME,
                                   1,
                                   moab::MB_TYPE_INTEGER,
                                   l_tagMoabGId,
                                   moab::MB_TAG_DENSE|moab::MB_TAG_CREAT );
  EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );

  // get the size of the tag and check with internal definition
  int l_gIdSize = 0;
  l_error = m_core.tag_get_bytes( l_tagMoabGId, l_gIdSize ); EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );
  CHECK( l_gIdSize == sizeof( int_gid ) );

  // MPI: create a consistent view of global id tag w.r.t. ghost elements and vertices
#ifdef PP_USE_MPI
  l_error = m_pcomm.exchange_tags( l_tagMoabGId, l_vertices ); EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );
  l_error = m_pcomm.exchange_tags( l_tagMoabGId, l_elements ); EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );
#endif

  // store mesh's global ids of els and ves before they get overwritten by assign_global_ids
  std::vector< int_gid > l_tmpGId( l_vertices.size() );
  l_error = m_core.tag_get_data( l_tagMoabGId, l_vertices, &l_tmpGId[0] ); EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );
  l_error = m_core.tag_set_data( m_tagGId,     l_vertices, &l_tmpGId[0] ); EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );

  l_tmpGId.resize( l_elements.size() );
  l_error = m_core.tag_get_data( l_tagMoabGId, l_elements, &l_tmpGId[0] ); EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );
  l_error = m_core.tag_set_data( m_tagGId,     l_elements, &l_tmpGId[0] ); EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );

#ifdef PP_USE_MPI
  // assign global ids to the entities
  l_error = m_pcomm.assign_global_ids( 0, m_dim-1 ); EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );
#endif

  // global ids for faces as determined by moab
  l_tmpGId.resize( l_faces.size() );
  l_error = m_core.tag_get_data( l_tagMoabGId, l_faces, &l_tmpGId[0] ); EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );
  l_error = m_core.tag_set_data( m_tagGId,     l_faces, &l_tmpGId[0] ); EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );

  // setup the data layout for the entities
  for( int_tg l_tg = 0; l_tg < l_nTimeGroups; l_tg++ ) {
    setupEnLayout( l_tg, l_vertices, m_veLayout, m_vertices );
    setupEnLayout( l_tg, l_faces,    m_faLayout, m_faces    );
    setupEnLayout( l_tg, l_elements, m_elLayout, m_elements );
  }
}

void edge::mesh::Moab::sortGId( std::vector< moab::EntityHandle > &io_ents ) {
  moab::ErrorCode l_error;

  // nothing todo for empty arrays
  if( io_ents.size() == 0 ) return;

     // get the global ids of the entities
  std::vector< int_gid > l_gIds;
  l_gIds.resize( io_ents.size() );

  l_error = m_core.tag_get_data(  m_tagGId,
                                 &io_ents[0],
                                  io_ents.size(),
                                 &l_gIds[0] );
  EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );

  // sort the entities by their id
  for( std::size_t l_en = 0; l_en < io_ents.size(); l_en++ ) {
   int_el l_minPos = l_en;
   // find position of minimum id of upcoming vertices
   for( std::size_t l_enUp = l_en+1; l_enUp < io_ents.size(); l_enUp++ ) {
     if( l_gIds[l_enUp] < l_gIds[l_minPos] ) {
       l_minPos = l_enUp;
     }
   }
   // exchange the global ids with the upcoming min
   int_gid l_tmpGid     = l_gIds[l_en];
   l_gIds[l_en]     = l_gIds[l_minPos];
   l_gIds[l_minPos] = l_tmpGid;

   // exchange ent with the upoming minimum
   moab::EntityHandle l_tmpEn = io_ents[l_en];
   io_ents[l_en]              = io_ents[l_minPos];
   io_ents[l_minPos]          = l_tmpEn;
  }
}

void edge::mesh::Moab::sortElFaEntities( std::vector< moab::EntityHandle > &io_faHas ) {
  moab::ErrorCode l_error;

  // global ids of the faces' vertices
  int_gid l_faVes[C_ENT[T_SDISC.ELEMENT].N_FACES][C_ENT[T_SDISC.FACE].N_VERTICES];
  assert( C_ENT[T_SDISC.ELEMENT].N_FACES >= io_faHas.size() );

  // iterate over faces and get the sorted ids of the vertices
  for( unsigned int l_fa = 0; l_fa < io_faHas.size(); l_fa++ ) {
    // get the handles of the face's vertices
    moab::Range l_verts;
    moab::EntityHandle l_faHan = io_faHas[l_fa];
    l_error = m_core.get_adjacencies( &l_faHan, 1, 0, true, l_verts, moab::Interface::UNION );
    EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );
    assert( l_verts.size() == C_ENT[T_SDISC.FACE].N_VERTICES );

    // get the global ids of the vertices
    l_error = m_core.tag_get_data( m_tagGId, l_verts, l_faVes[l_fa] ); EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );

    // sort them for this face
    for( unsigned int l_ve = 0; l_ve < C_ENT[T_SDISC.FACE].N_VERTICES; l_ve++ ) {
      unsigned int l_minPos = l_ve;
      for( unsigned int l_veUp = l_ve+1; l_veUp < C_ENT[T_SDISC.FACE].N_VERTICES; l_veUp++ ) {
        if( l_faVes[l_fa][l_minPos] > l_faVes[l_fa][l_veUp] ) {
          l_minPos = l_veUp;
        }
      }
      // exchange values
      int_gid l_veTmp = l_faVes[l_fa][l_ve];
      l_faVes[l_fa][l_ve] = l_faVes[l_fa][l_minPos];
      l_faVes[l_fa][l_minPos] = l_veTmp;
    }
  }

  // lets go, sort the entities
  for( std::size_t l_fa = 0; l_fa < io_faHas.size(); l_fa++ ) {
    unsigned int l_minPos = l_fa;
    for( std::size_t l_faUp = l_fa+1; l_faUp < io_faHas.size(); l_faUp++ ) {
      // iterate over vertices and check if the next face has a smaller combined index
      for( unsigned int l_ve = 0; l_ve < C_ENT[T_SDISC.FACE].N_VERTICES; l_ve++ ) {
        if( l_faVes[l_minPos][l_ve] > l_faVes[l_faUp][l_ve] ) {
          l_minPos = l_faUp;
          break;
        }
        else if ( l_faVes[l_minPos][l_ve] < l_faVes[l_faUp][l_ve] ) {
          break;
        }
      }
    }

    // do the exchange w.r.t to the array of vertices
    for( unsigned int l_ve = 0; l_ve < C_ENT[T_SDISC.FACE].N_VERTICES; l_ve++ ) {
      int_gid l_veTmp = l_faVes[l_fa][l_ve];
      l_faVes[l_fa][l_ve] = l_faVes[l_minPos][l_ve];
      l_faVes[l_minPos][l_ve] = l_veTmp;
    }

    // accomplish our original goal: sorting the entity ids
    moab::EntityHandle l_entTmp = io_faHas[l_fa];
    io_faHas[l_fa] = io_faHas[l_minPos];
    io_faHas[l_minPos] = l_entTmp;
  }
}

void edge::mesh::Moab::read( const std::string &i_pathToMesh,
                             const std::string &i_optRead ) {
  moab::ErrorCode l_error;
  // load the mesh
  EDGE_LOG_INFO << "  loading mesh";
  l_error = m_core.load_file( i_pathToMesh.c_str(), 0, i_optRead.c_str() ); EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );

  EDGE_LOG_INFO << "  processing mesh info";
  // get tags defined by the mesh
  l_error = m_core.tag_get_tags( m_tagsMesh ); EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );

  // get names of tags
  m_tagsNames.resize( m_tagsMesh.size() );
  for( std::size_t l_ta = 0; l_ta < m_tagsMesh.size(); l_ta++ ) {
    l_error = m_core.tag_get_name( m_tagsMesh[l_ta], m_tagsNames[l_ta] ); EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );
  }

  // set material tag
  m_tagMat = m_tagsMesh[0];

  // get name of material tag
  std::string l_matName;
  l_error = m_core.tag_get_name( m_tagMat, l_matName ); EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );

  // double-check that this is the right one
  if( l_matName != "MATERIAL_SET" ) {
    EDGE_LOG_ERROR << "  material tag name is not matching";
  }
  // get the size of the tag and check with internal definition
  int l_matSize = 0;
  l_error = m_core.tag_get_bytes( m_tagMat, l_matSize ); EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );

  // check that the tag has the size of an int
  EDGE_CHECK_EQ( (std::size_t) l_matSize, sizeof( int ) );
  // check that the tag fits in our sparse type
  EDGE_CHECK_LT( (std::size_t) l_matSize, sizeof( int_spType ) );

  // setup our custom data layout
  setupDataLayout();

  // get the boundary mesh sets
  if( m_nBndVals > 0 ) {
    moab::Range l_mSets;

    // get all meshsets
    l_error = m_core.get_entities_by_type( 0, moab::MBENTITYSET, l_mSets );
    EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );

    // get mat-tag data
    std::vector< int > l_matData; l_matData.resize( l_mSets.size() );
    l_error = m_core.tag_get_data( m_tagMat, l_mSets, &l_matData[0] );
    EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );

    for( std::size_t l_se = 0; l_se < l_mSets.size(); l_se++ ) {
      // add to boundary sets if it matches input data
      for( unsigned int l_bv = 0; l_bv < m_nBndVals; l_bv++ ) {
        if( l_matData[l_se] == m_bndVals[l_bv] ) m_bndSets.insert( l_mSets[l_se] );
      }
    }
    EDGE_CHECK( m_bndSets.size() <= m_nBndVals );
  }

#ifndef PP_USE_MPI
  // hack in periodic boundary conditions, this can be removed if MOAB get native support for those..
  EDGE_VLOG(3) << "  setting custom periodic boundaries";
  setPeriodicAdjacencies( m_periodicVal );
#endif

  // initialize
  init();

#ifndef PP_USE_MPI
  // hack in ids for removed periodic faces
  fixPeriodicTagIds();
#endif
}

void edge::mesh::Moab::addPeriodicVertices( moab::EntityHandle  i_ve,
                                            moab::Range        &io_ves ) {

  // return if no periodic faces are available
  if( m_periodicFaces.size() == 0 ) return;

  // vertices to test for
  moab::Range l_test; l_test.insert( i_ve );

  // error code
  moab::ErrorCode l_err;

  // go through the process multiple times to allow for multi-hops crossing all dimensions
  for( unsigned short l_rp = 0; l_rp < m_dim; l_rp++ ) {
    // iterate over test vertices
    for( std::size_t l_te = 0; l_te < l_test.size(); l_te++ ) {
      // test entity
      moab::EntityHandle l_teEn = l_test[l_te];

      // query coordinates of the test vertex
      double l_veCrds[3];
      l_err = m_core.get_coords( &l_teEn, 1, l_veCrds );
      EDGE_CHECK_EQ( l_err, moab::MB_SUCCESS );

      // determine adjacent faces
      moab::Range l_veFa;
      l_err = m_core.get_adjacencies( &l_teEn, 1, m_dim-1, true, l_veFa, moab::Interface::UNION );
      EDGE_CHECK_EQ( l_err, moab::MB_SUCCESS );

      // determine periodic equivalents
      moab::Range l_veFaP;

      for( std::size_t l_fa = 0; l_fa < m_periodicFaces.size(); l_fa++ ) {
        for( std::size_t l_vf = 0; l_vf < l_veFa.size(); l_vf++ ) {
          if( m_periodicFaces[l_fa] == l_veFa[l_vf] ) {
            // add periodic face
            if( l_fa%2 == 0 ) l_veFaP.insert( m_periodicFaces[l_fa+1] );
            else              l_veFaP.insert( m_periodicFaces[l_fa-1] );
          }
        }
      }

      // determine face-vertices, which match the queried one
      moab::Range l_fpVe;
      l_err = m_core.get_adjacencies( l_veFaP, 0, true, l_fpVe, moab::Interface::UNION );
      EDGE_CHECK_EQ( l_err, moab::MB_SUCCESS );

      if( l_fpVe.size() > 0 ) {
        std::vector< double > l_veCrdsP( 3*l_fpVe.size() );
        l_err = m_core.get_coords( l_fpVe, &l_veCrdsP[0] );
        EDGE_CHECK_EQ( l_err, moab::MB_SUCCESS );

        // iterate over vertices
        for( std::size_t l_ve = 0; l_ve < l_fpVe.size(); l_ve++ ) {
          // check for the number of matching dimensions
          unsigned l_nDi = 0;

          for( unsigned short l_di = 0; l_di < m_dim; l_di++ ) {
            if( std::abs( l_veCrds[l_di] - l_veCrdsP[l_ve*3 + l_di] ) < TOL.MESH )
              l_nDi++;
          }

          // check that at most #dims-1 dimensions match
          EDGE_CHECK_LT( l_nDi, m_dim );

          // if only one dim doesnt match, we have our periodic vertex
          if( l_nDi == m_dim-1 ) {
            l_test.insert( l_fpVe[l_ve] );
          }
        }
      }

    }
  }

  // remove the input vertex
  l_test.erase( i_ve );

  // check number of found periodic vertices
  if(      m_dim == 2 ) { EDGE_CHECK_LE( l_test.size(), 3 ); }
  else if( m_dim == 3 ) { EDGE_CHECK_LE( l_test.size(), 7 ); }

  // save results
  io_ves.merge( l_test );
}

void edge::mesh::Moab::addPeriodicElement( moab::EntityHandle  i_face,
                                           moab::Range        &io_elements ) {
  for( std::size_t l_fa = 0; l_fa < m_periodicFaces.size(); l_fa++ ) {
    if( m_periodicFaces[l_fa] == i_face ) {
      io_elements.insert( m_periodicElementsRemote[l_fa] );
      return;
    }
  }
}

void edge::mesh::Moab::addPeriodicElement( moab::EntityHandle                 i_face,
                                           std::vector< moab::EntityHandle > &io_elements ) {
  for( std::size_t l_fa = 0; l_fa < m_periodicFaces.size(); l_fa++ ) {
    if( m_periodicFaces[l_fa] == i_face ) {
      io_elements.push_back( m_periodicElementsRemote[l_fa] );
      return;
    }
  }
}

void edge::mesh::Moab::write( const std::string &i_pathToMesh ) {
  moab::ErrorCode l_error;

  std::string l_pathToMesh = std::to_string(edge::parallel::g_rank)+"_"+i_pathToMesh;

  // write mesh to disk
  EDGE_LOG_INFO_ALL << "writing mesh " << l_pathToMesh;
  l_error = m_core.write_mesh( l_pathToMesh.c_str() );
  EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );
}

void edge::mesh::Moab::getElVeEl( int_el  &o_nElVeEl,
                                  int_el **o_elVeEl ) {
  // init element counter
  o_nElVeEl = 0;

  // raw ptr to adjaceny info
  int_el *l_raw = nullptr;
  if( o_elVeEl != nullptr ) l_raw = o_elVeEl[0];

  // moab error code
  moab::ErrorCode l_err;

  // vertices
  std::vector< moab::EntityHandle > l_ves;

  // elements
  std::vector< moab::EntityHandle > l_els;

  // iterate over elements
  for( std::size_t l_el = 0; l_el < m_elements.size(); l_el++ ) {
    // get entity
    moab::EntityHandle l_elEnt = m_elements[l_el];

    // get adjacent vertices
    l_ves.clear();
    l_err = m_core.get_adjacencies( &l_elEnt, 1, 0, true, l_ves, moab::Interface::UNION );
    EDGE_CHECK_EQ( l_err, moab::MB_SUCCESS ); EDGE_CHECK_GT( l_ves.size(), m_dim );

    // take care of periodic boundaries
    moab::Range l_vesP;
    for( std::size_t l_ve = 0; l_ve < l_ves.size(); l_ve++ ) {
      addPeriodicVertices( l_ves[l_ve], l_vesP );
    }
    // add to vector
    for( std::size_t l_vp = 0; l_vp < l_vesP.size(); l_vp++ ) {
      l_ves.push_back( l_vesP[l_vp] );
    }

    // get elements adjacent to the vertices
    l_els.clear();
    l_err = m_core.get_adjacencies( &l_ves[0], l_ves.size(), m_dim, true, l_els, moab::Interface::UNION );
    EDGE_CHECK_EQ( l_err, moab::MB_SUCCESS ); EDGE_CHECK_GT( l_els.size(), m_dim );

    // sort elements by their global id
    sortGId( l_els );

    // store entries
    for( std::size_t l_ad = 0; l_ad < l_els.size(); l_ad++ ) {
      if( l_els[l_ad] != l_elEnt ) {
        // query local id
        if( l_raw != nullptr ) {
          l_err = m_core.tag_get_data( m_tagLId, &l_els[l_ad], 1, l_raw );
          EDGE_CHECK_EQ( l_err, moab::MB_SUCCESS );

          l_raw++;
        }

        // increase element count
        o_nElVeEl++;
      }
    }

    // assign next pointer
    if( o_elVeEl != nullptr && l_raw != nullptr ) o_elVeEl[l_el+1] = l_raw;
  }
}

void edge::mesh::Moab::getElVeEl( int_el **o_elVeEl ) {
  // number of adjacent elements
  int_el l_nElVeEl;

  // get adjacent elements
  getElVeEl( l_nElVeEl, o_elVeEl );

  // check size
  EDGE_CHECK_EQ( o_elVeEl[ m_elements.size() ] - o_elVeEl[0], l_nElVeEl );
}

int_el edge::mesh::Moab::getNelVeEl() {
  // number of adjacent elements
  int_el l_nElVeEl;

  getElVeEl( l_nElVeEl, nullptr );

  return l_nElVeEl;
}

void edge::mesh::Moab::getElementFaceNeighbors( moab::EntityHandle i_element,
                                                int_el             o_neighboringIds[C_ENT[T_SDISC.ELEMENT].N_FACES] ) {
  moab::ErrorCode l_error;

  // get the faces
  std::vector< moab::EntityHandle > l_faces;
  l_error = m_core.get_adjacencies( &i_element, 1, m_dim-1, true, l_faces, moab::Interface::UNION  );
  EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );
  assert( l_faces.size() == C_ENT[T_SDISC.ELEMENT].N_FACES );

  // sort the faces by the global ids of their vertices
  sortElFaEntities( l_faces );

  // iterate over faces and get the corresponding local ids of the adjacent elements
  for( std::size_t l_fa = 0; l_fa < l_faces.size(); l_fa++ ) {
    moab::Range l_faceNeighbor;

    // get the face neighbor
    moab::EntityHandle l_face = l_faces[l_fa];
    l_error = m_core.get_adjacencies( &l_face, 1, m_dim, true, l_faceNeighbor, moab::Interface::UNION  );
    EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );

    // remove current element
    l_faceNeighbor.erase( i_element );
    assert( l_faceNeighbor.size() <= 1 );

    // try to add periodic neighbors, could also be a non-periodic boundary though
    addPeriodicElement( l_face, l_faceNeighbor );

    if( l_faceNeighbor.size() != 0 ) {
      // get local id of the face neighbor
      l_error = m_core.tag_get_data(  m_tagLId, l_faceNeighbor, &o_neighboringIds[l_fa] );
      EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );
    }
    else {
      o_neighboringIds[l_fa] = std::numeric_limits<int_el>::max();
    }
  }
}

void edge::mesh::Moab::getElementsFaceNeighbors( int_el (*o_neighboringIds)[C_ENT[T_SDISC.ELEMENT].N_FACES] ) {
  // iterate over elements
  for( std::size_t l_el = 0; l_el < m_elements.size(); l_el++ ) {
    // get entity of this element
    moab::EntityHandle l_elEnt = m_elements[l_el];

    // get face neighbors of the element
    getElementFaceNeighbors( l_elEnt, o_neighboringIds[l_el] );
  }
}

void edge::mesh::Moab::getMeshRegions( std::vector<int> &o_values ) {
  moab::ErrorCode l_error;

  // get the entity sets with material tags
  moab::Range l_sets;
  l_error = m_core.get_entities_by_type_and_tag(0, moab::MBENTITYSET, &m_tagMat, NULL, 1, l_sets);
  EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );

  // resize output values to match the #sets
  o_values.resize( l_sets.size() );

  // get material values
  l_error = m_core.tag_get_data( m_tagMat, l_sets, &o_values[0] );
  EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );
}

void edge::mesh::Moab::getElVe( int_el (*o_elVe)[C_ENT[T_SDISC.ELEMENT].N_VERTICES] ) {
  moab::ErrorCode l_error;

  for( std::size_t l_el = 0; l_el < m_elements.size(); l_el++ ) {
    // get the vertex neighbors
    moab::EntityHandle l_elHan = m_elements[l_el];
    std::vector< moab::EntityHandle > l_vertices;
    l_error = m_core.get_adjacencies( &l_elHan, 1, 0, true, l_vertices, moab::Interface::UNION );
    EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );
    // check that we have everybody
    EDGE_CHECK_EQ( l_vertices.size(), C_ENT[T_SDISC.ELEMENT].N_VERTICES );

    // sort vertices by their global id
    sortGId( l_vertices );

    // now get the local ids of the vertices
    l_error = m_core.tag_get_data(  m_tagLId, &l_vertices[0], l_vertices.size(), o_elVe[l_el] );
    EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );
  }
}

void edge::mesh::Moab::getElementAdjacentFaces( moab::EntityHandle i_element,
                                                int_el             o_elementAdjacentFaces[C_ENT[T_SDISC.ELEMENT].N_FACES] ) {
  moab::ErrorCode l_error;
  // get the faces
  std::vector< moab::EntityHandle > l_faces;
  l_error = m_core.get_adjacencies( &i_element, 1, m_dim-1, true, l_faces, moab::Interface::UNION  );
  EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );
  EDGE_CHECK_EQ( l_faces.size(), C_ENT[T_SDISC.ELEMENT].N_FACES );

  // sort the faces ascending by their vertices' global ids
  sortElFaEntities( l_faces );

  // get their local ids
  l_error = m_core.tag_get_data( m_tagLId, &l_faces[0], l_faces.size(), o_elementAdjacentFaces );
  EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );

  // get mpi type of the element
  unsigned short l_mpiType = 0;
#ifdef PP_USE_MPI
  l_mpiType = getEnMpiType( i_element );
#endif

  // make sure that every face is accounted for
  for( unsigned short l_fa = 0; l_fa < C_ENT[T_SDISC.ELEMENT].N_FACES; l_fa++ ) {
    EDGE_CHECK( l_mpiType != 0 ||
                o_elementAdjacentFaces[l_fa] != std::numeric_limits<int_el>::max() ) << i_element;
  }
}

void edge::mesh::Moab::getElementsAdjacentFaces( int_el (*o_elementsAdjacentFaces)[C_ENT[T_SDISC.ELEMENT].N_FACES] ) {
  // iterate over elements
  for( std::size_t l_el = 0; l_el < m_elements.size(); l_el++ ) {
    // get entity of this element
    moab::EntityHandle l_elEnt = m_elements[l_el];

    // get faces neighboring the element
    getElementAdjacentFaces( l_elEnt,
                             o_elementsAdjacentFaces[l_el] );
  }
}

void edge::mesh::Moab::getFacesAdjacentElements( int_el (*o_neighboringIds)[2] ) {
  moab::ErrorCode l_error;

  // iterate over faces
  for( std::size_t l_fa = 0; l_fa < m_faces.size(); l_fa++ ) {
    // get entity of this face
    moab::EntityHandle l_faHan = m_faces[l_fa];

    // get id of the face (precaution)
    int_el l_localId;
    l_error = m_core.tag_get_data( m_tagLId, &l_faHan, 1, &l_localId ); EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );

    // get parallel status for MPI-"elements" which might be duplicated.
    unsigned char l_pstatus = 0;
#ifdef PP_USE_MPI
    l_error = m_pcomm.get_pstatus( l_faHan, l_pstatus ); EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );
#endif

    // check it matches the local id or is send-/receive-ent
    assert( l_localId == (int_el) l_fa || l_pstatus != 0 );

    // get adjacent elements
    std::vector< moab::EntityHandle > l_elements;
    l_error = m_core.get_adjacencies( &l_faHan, 1, m_dim, true, l_elements, moab::Interface::UNION );
    EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );
    assert( l_elements.size() > 0 && l_elements.size() <= 2 );

    // try to find second face neighbor for periodic faces
    addPeriodicElement( l_faHan, l_elements );
    assert( l_elements.size() <= 2 );

    // set 2nd neighbor invalid for non-existing neighbors
    if( l_elements.size() == 1 ) o_neighboringIds[l_fa][1] = std::numeric_limits<int_el>::max();

    // get local ids of neighboring elements
    l_error = m_core.tag_get_data( m_tagLId, &l_elements[0], l_elements.size(), o_neighboringIds[l_fa] );
    EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );

    // get global ids of neighboring elements
    int_gid l_gId[2];
    l_error = m_core.tag_get_data( m_tagGId, &l_elements[0], l_elements.size(), l_gId ); EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );


    // enforce an ascending ordering
    if( l_elements.size() == 2 && l_gId[0] > l_gId[1] ) {
      int_el l_elTmp = o_neighboringIds[l_fa][1];
      o_neighboringIds[l_fa][1] = o_neighboringIds[l_fa][0];
      o_neighboringIds[l_fa][0] = l_elTmp;
    }
  }
}

void edge::mesh::Moab::getFacesAdjacentVertices( int_el (*o_neighVeIds)[C_ENT[T_SDISC.FACE].N_VERTICES] )  {
  moab::ErrorCode l_error;

  // iterate over faces
  for( std::size_t l_fa = 0; l_fa < m_faces.size(); l_fa++ ) {
    // get entity of this face
    moab::EntityHandle l_faHan = m_faces[l_fa];

    // get id of the face (pre caution)
    int_el l_localId;
    l_error = m_core.tag_get_data( m_tagLId, &l_faHan, 1, &l_localId ); EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );

    // get parallel status for MPI-"elements" which might be duplicated.
    unsigned char l_pstatus = 0;
#ifdef PP_USE_MPI
    l_error = m_pcomm.get_pstatus( l_faHan, l_pstatus ); EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );
#endif
    assert( l_localId == (int_el) l_fa || l_pstatus != 0 );

    // get neighboring vertices
    std::vector< moab::EntityHandle > l_verts;
    l_error = m_core.get_adjacencies( &l_faHan, 1, 0, true, l_verts, moab::Interface::UNION );
    EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );

    // check that we have our verts
    assert( l_verts.size() == C_ENT[T_SDISC.FACE].N_VERTICES );

    sortGId( l_verts );

    // now get the local ids of the vertices
    l_error = m_core.tag_get_data(  m_tagLId,
                                   &l_verts[0],
                                    C_ENT[T_SDISC.FACE].N_VERTICES,
                                    o_neighVeIds[l_fa] );
    EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );
  }
}

void edge::mesh::Moab::getElementVerticesCoords( int_el i_element,
                                                 double o_x[C_ENT[T_SDISC.ELEMENT].N_VERTICES],
                                                 double o_y[C_ENT[T_SDISC.ELEMENT].N_VERTICES],
                                                 double o_z[C_ENT[T_SDISC.ELEMENT].N_VERTICES] ) {
  moab::ErrorCode l_error;

  // get the vertices
  moab::EntityHandle l_element = m_elements[i_element];
  moab::Range l_vertices;
  l_error = m_core.get_adjacencies( &l_element, 1, 0, true, l_vertices, moab::Interface::UNION );
  EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );
  assert( l_vertices.size() == C_ENT[T_SDISC.ELEMENT].N_VERTICES );

  // get the coordinates
  l_error = m_core.get_coords( l_vertices, o_x, o_y, o_z ); EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );
}

void edge::mesh::Moab::getGIdsVe( std::vector< int_gid > &o_gIds ) const {
  moab::ErrorCode l_error;
  o_gIds.resize( m_vertices.size() );

  // get the global ids
  l_error = m_core.tag_get_data(  m_tagGId, &m_vertices[0], m_vertices.size(), &o_gIds[0] );
  EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );
}

void edge::mesh::Moab::getGIdsFa( std::vector< int_gid > &o_gIds ) const {
  moab::ErrorCode l_error;
  o_gIds.resize( m_faces.size() );

  // get the global ids
  l_error = m_core.tag_get_data(  m_tagGId, &m_faces[0], m_faces.size(), &o_gIds[0] );
  EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );
}

void edge::mesh::Moab::getGIdsEl( std::vector< int_gid > &o_gIds ) const {
  moab::ErrorCode l_error;
  o_gIds.resize( m_elements.size() );

  // get the global ids
  l_error = m_core.tag_get_data(  m_tagGId, &m_elements[0], m_elements.size(), &o_gIds[0] );
  EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );
}

void edge::mesh::Moab::syncTags( std::vector< unsigned short > const &i_tids,
                                 unsigned short                       i_nDim  ) {
#ifdef PP_USE_MPI
  moab::ErrorCode l_err;

  // get MOAB tag-handles from local ids
  std::vector< moab::Tag > l_tags;
  for( std::size_t l_ta = 0; l_ta < i_tids.size(); l_ta++ )
    l_tags.push_back( m_tagsMesh[ i_tids[l_ta] ] );

  moab::Range l_enMpi;
  // get handles of given dimension
  l_err = m_core.get_entities_by_dimension( 0, i_nDim, l_enMpi ); EDGE_CHECK_EQ( l_err, moab::MB_SUCCESS ) << l_err;

  // filter handles by MPI-status
  unsigned char l_psMpi = PSTATUS_GHOST | PSTATUS_SHARED | PSTATUS_MULTISHARED | PSTATUS_INTERFACE;
  l_err = m_pcomm.filter_pstatus( l_enMpi, l_psMpi, PSTATUS_OR, -1 ); EDGE_CHECK_EQ( l_err, moab::MB_SUCCESS ) << l_err;

  // exchange tags
  l_err = m_pcomm.exchange_tags( l_tags, l_tags, l_enMpi ); EDGE_CHECK_EQ( l_err, moab::MB_SUCCESS ) << l_err;
#endif
}

void edge::mesh::Moab::getTagsNames( std::vector< std::string > &o_tagsNames ) {
  // copy values
  o_tagsNames = m_tagsNames;
}

int edge::mesh::Moab::getTagBytes( unsigned short i_tid ) {
  moab::ErrorCode l_err;

  int l_bytes = std::numeric_limits< int >::max();

  l_err = m_core.tag_get_bytes( m_tagsMesh[i_tid], l_bytes ); EDGE_CHECK_EQ( l_err, moab::MB_SUCCESS );

  return l_bytes;
}

void edge::mesh::Moab::getTagData( unsigned short  i_tid,
                                   unsigned short  i_nDim,
                                   int_el          i_en,
                                   void           *o_val ) {
  // error code
  moab::ErrorCode l_err;

  // get entity handle
  moab::EntityHandle l_ha = getEnHandle( i_nDim, i_en );

  // query for value
  l_err = m_core.tag_get_data(  m_tagsMesh[i_tid],
                               &l_ha,
                                1,
                                o_val ); EDGE_CHECK_EQ( l_err, moab::MB_SUCCESS );
}

void edge::mesh::Moab::getMatVal( unsigned short  i_dim,
                                  int_el          i_ent,
                                  int_spType     &o_val ) {
  moab::ErrorCode l_error;

  // get entity handle
  moab::EntityHandle l_ha = getEnHandle( i_dim, i_ent );

  // get the entity sets with material tags
  moab::Range l_sets;
  l_error = m_core.get_entities_by_type_and_tag( 0, moab::MBENTITYSET, &m_tagMat, NULL, 1, l_sets );
  EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );

  // check if the entity is part of one of sets
  unsigned int l_conSet = std::numeric_limits<unsigned int>::max();
  for( std::size_t l_set = 0; l_set < l_sets.size(); l_set++ ) {
    if( m_core.contains_entities( l_sets[l_set], &l_ha, 1 ) ) {
      l_conSet = l_set;
      break;
    }
  }

  // default case: no material value define, set default
  if( l_conSet == std::numeric_limits<unsigned int>::max() ) {
    o_val = MESH_TYPE_NONE;
  }
  // query moab for mat val
  else {
    moab::EntityHandle l_setHa = l_sets[l_conSet];

    int l_val;
    l_error = m_core.tag_get_data( m_tagMat, &l_setHa, 1, &l_val );
    EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );
    o_val = l_val; // change type of value
  }

  // negative values are not allowed since the last half of bits is reserved for apps.
  EDGE_CHECK_GE( o_val, 0 ) << "Only non-negative values are allowed for sparse types as mesh-input";
}

void edge::mesh::Moab::getVeCoords( int_el  i_vertex,
                                    double &o_x,
                                    double &o_y,
                                    double &o_z ) {
  moab::ErrorCode l_error;

  moab::EntityHandle l_veHan = m_vertices[i_vertex];
  double l_coords[3];
  l_error = m_core.get_coords( &l_veHan, 1, l_coords ); EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );

  o_x = l_coords[0]; o_y = l_coords[1]; o_z = l_coords[2];
}

void edge::mesh::Moab::getFaceVerticesCoords( int_el i_face,
                                              double o_x[C_ENT[T_SDISC.FACE].N_VERTICES],
                                              double o_y[C_ENT[T_SDISC.FACE].N_VERTICES],
                                              double o_z[C_ENT[T_SDISC.FACE].N_VERTICES] ) {
  moab::ErrorCode l_error;

  //get the vertices
  moab::EntityHandle l_face = m_faces[i_face];
  moab::Range l_vertices;
  l_error = m_core.get_adjacencies( &l_face, 1, 0, true, l_vertices, moab::Interface::UNION );
  EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );
  assert( l_vertices.size() == C_ENT[T_SDISC.FACE].N_VERTICES );

  // get the coordinates
  l_error = m_core.get_coords( l_vertices, o_x, o_y, o_z );
  EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );
}

void edge::mesh::Moab::getNormalPointCoords( int_el  i_face,
                                             bool   &o_left,
                                             double &o_x,
                                             double &o_y,
                                             double &o_z ) {
#if defined PP_T_ELEMENTS_TRIA3 || defined PP_T_ELEMENTS_TET4
  moab::ErrorCode l_error;

  // get adjacent elements
  std::vector< moab::EntityHandle > l_elements;

  EDGE_CHECK( i_face < (int_el) m_faces.size() );

  moab::EntityHandle l_faHan = m_faces[i_face];
  l_error = m_core.get_adjacencies( &l_faHan, 1, m_dim, true, l_elements, moab::Interface::UNION );
  EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );
  assert( l_elements.size() >= 1 );

  addPeriodicElement( m_faces[i_face], l_elements );

  // face has to be at the boundary or define two neighbors
  EDGE_CHECK( l_elements.size() == 2 || (isBnd(l_faHan) && l_elements.size() == 1) )
  << "before reporting this, double-check that your mesh defines all necessary boundary-conditions";

  // get global ids and sort ascending
  int_gid l_gIds[2];
  l_error = m_core.tag_get_data( m_tagGId, &l_elements[0], l_elements.size(), l_gIds );
  EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );

  if( l_elements.size() == 2 ) o_left = l_gIds[0] < l_gIds[1];
  else                         o_left = true; // boundary elements are left by default


  // handle of the first element
  moab::EntityHandle l_elFirst = l_elements[0];

  // get vertices of first element
  moab::Range l_veFirst;
  l_error = m_core.get_adjacencies( &l_elFirst, 1, 0, true, l_veFirst, moab::Interface::UNION );
  EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );
  assert( l_veFirst.size() == C_ENT[T_SDISC.ELEMENT].N_VERTICES );

  // get the vertices of the face
  moab::Range l_veFace;
  l_error = m_core.get_adjacencies( &l_faHan, 1, 0, true, l_veFace, moab::Interface::UNION );
  EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );
  assert( l_veFace.size() == C_ENT[T_SDISC.FACE].N_VERTICES );

  // get remaining vertex
  l_veFirst -= l_veFace;
  assert( l_veFirst.size() == 1);

  // get the coordinates and out of here!
  l_error = m_core.get_coords( l_veFirst, &o_x, &o_y, &o_z );
  EDGE_CHECK_EQ( l_error, moab::MB_SUCCESS );
#else
  assert( false );
#endif
}
