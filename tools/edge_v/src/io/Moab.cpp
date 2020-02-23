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
 * MOAB interface.
 **/
#include "Moab.h"
#include <moab/ReorderTool.hpp>
#include <iomanip>

moab::EntityType edge_v::io::Moab::getMoabType( t_entityType i_enTy ) {
  moab::EntityType l_ty = moab::MBVERTEX;

  if( i_enTy == LINE ) {
    l_ty = moab::MBEDGE;
  }
  else if( i_enTy == QUAD4R ) {
    l_ty = moab::MBQUAD;
  }
  else if( i_enTy == TRIA3 ) {
    l_ty = moab::MBTRI;
  }
  else if( i_enTy == HEX8R ) {
    l_ty = moab::MBHEX;
  }
  else if( i_enTy == TET4 ) {
    l_ty = moab::MBTET;
  }
  else EDGE_V_CHECK_EQ( i_enTy, POINT );

  return l_ty;
}

void edge_v::io::Moab::getMemStats( unsigned long long & o_totalStorage,
                                    unsigned long long & o_totalAmortized,
                                    unsigned long long & o_entityStorage,
                                    unsigned long long & o_entityAmortized,
                                    unsigned long long & o_adjacencyStorage,
                                    unsigned long long & o_adjacencyAmortized,
                                    unsigned long long & o_tagStorage,
                                    unsigned long long & o_tagAmortized ) const {
  m_moab->estimated_memory_use( 0,
                                0,
                                &o_totalStorage,
                                &o_totalAmortized,
                                &o_entityStorage,
                                &o_entityAmortized,
                                &o_adjacencyStorage,
                                &o_adjacencyAmortized,
                                0,
                                0,
                                &o_tagStorage,
                                &o_tagAmortized );
}

edge_v::io::Moab::Moab( std::string const & i_pathToMesh ) {
  // init moab
  m_moab = new moab::Core;
  moab::ErrorCode l_err = m_moab->load_file( i_pathToMesh.c_str() );
  EDGE_V_CHECK_EQ( l_err, moab::MB_SUCCESS );

  m_root = m_moab->get_root_set();

  // get number of dimensions (or what MOAB claims)
  int l_nDis = std::numeric_limits< int >::max();
  l_err = m_moab->get_dimension( l_nDis );
  EDGE_V_CHECK_EQ( l_err, moab::MB_SUCCESS );

  // get elements and number of dimensions
  std::vector< moab::EntityHandle > l_els;
  for( unsigned short l_di = l_nDis; l_di > 0; l_di-- ) {
    l_els.resize(0);
    l_err = m_moab->get_entities_by_dimension( m_root,
                                               l_di,
                                               l_els );
    EDGE_V_CHECK_EQ( l_err, moab::MB_SUCCESS );

    // abort if we found the element
    if( l_els.size() > 0 ) {
      l_nDis = l_di;
      break;
    }
  }

  // get vertices
  moab::Range l_ves;
  l_err = m_moab->get_entities_by_dimension( m_root,
                                             0,
                                             l_ves );
  EDGE_V_CHECK_EQ( l_err, moab::MB_SUCCESS );

  // get faces
  moab::Range l_fas;
  l_err = m_moab->get_adjacencies( l_ves,
                                   l_nDis-1,
                                   true,
                                   l_fas,
                                   moab::Interface::UNION );
  EDGE_V_CHECK_EQ( l_err, moab::MB_SUCCESS );

  // get the number of entities
  int l_nEns; // moab forces int-queries
  m_moab->get_number_entities_by_type( m_root,
                                       getMoabType( POINT ),
                                       l_nEns );
  m_nPoint = l_nEns;

  m_moab->get_number_entities_by_type( m_root,
                                       getMoabType( LINE ),
                                       l_nEns );
  m_nLine = l_nEns;

  m_moab->get_number_entities_by_type( m_root,
                                       getMoabType( QUAD4R ),
                                       l_nEns );
  m_nQuad4r = l_nEns;

  m_moab->get_number_entities_by_type( m_root,
                                       getMoabType( TRIA3 ),
                                       l_nEns );
  m_nTria3 = l_nEns;

  m_moab->get_number_entities_by_type( m_root,
                                       getMoabType( HEX8R ),
                                       l_nEns );
  m_nHex8r = l_nEns;

  m_moab->get_number_entities_by_type( m_root,
                                       getMoabType( TET4 ),
                                       l_nEns );
  m_nTet4 = l_nEns;
}

edge_v::io::Moab::~Moab() {
  delete m_moab;
}

void edge_v::io::Moab::printStats() const {
  EDGE_V_LOG_INFO << "moab stats:";
  EDGE_V_LOG_INFO << "  #point:  " << m_nPoint;
  EDGE_V_LOG_INFO << "  #line:   " << m_nLine;
  EDGE_V_LOG_INFO << "  #quad4r: " << m_nQuad4r;
  EDGE_V_LOG_INFO << "  #tria3:  " << m_nTria3;
  EDGE_V_LOG_INFO << "  #hex8r:  " << m_nHex8r;
  EDGE_V_LOG_INFO << "  #tet4:   " << m_nTet4;
  std::vector< std::string > l_tagNames;
  getTagNames( l_tagNames );
  if( l_tagNames.size() > 1 ) {
    EDGE_V_LOG_INFO << "  defined tags:";
    for( std::size_t l_ta = 0; l_ta < l_tagNames.size(); l_ta++ ) {
      EDGE_V_LOG_INFO << "    #" << l_ta << ": " << l_tagNames[l_ta];
    }
  }

  // get memory statistics
  unsigned long long l_totalStorage,     l_totalAmortized,
                     l_entityStorage,    l_entityAmortized,
                     l_adjacencyStorage, l_adjacencyAmortized,
                     l_tagStorage,       l_tagAmortized;

  getMemStats( l_totalStorage,
               l_totalAmortized,
               l_entityStorage,
               l_entityAmortized,
               l_adjacencyStorage,
               l_adjacencyAmortized,
               l_tagStorage,
               l_tagAmortized );
  EDGE_V_LOG_INFO << "  memory (used / allocated in GiB):";
  EDGE_V_LOG_INFO << "    total: "     << std::setprecision(3)
                                       << byteToGib(l_totalStorage) << " / "
                                       << byteToGib(l_totalAmortized);
  EDGE_V_LOG_INFO << "    entity: "    << std::setprecision(3)
                                       << byteToGib(l_entityStorage) << " / "
                                       << byteToGib(l_entityAmortized);
  EDGE_V_LOG_INFO << "    adjacency: " << std::setprecision(3)
                                       << byteToGib(l_adjacencyStorage) << " / "
                                       << byteToGib(l_adjacencyAmortized);
  EDGE_V_LOG_INFO << "    tag: "       << std::setprecision(3)
                                       << byteToGib(l_tagStorage) << " / "
                                       << byteToGib(l_tagAmortized);
}

edge_v::t_entityType edge_v::io::Moab::getElType() const {
  // check entity counts for mesh type
  if( m_nTet4 > 0 ) {
    EDGE_V_CHECK_EQ( m_nHex8r,  0 );
    EDGE_V_CHECK_EQ( m_nQuad4r, 0 );

    EDGE_V_CHECK_GT( m_nTria3,  0 );

    return TET4;
  }
  else if( m_nHex8r > 0 ) {
    EDGE_V_CHECK_EQ( m_nTria3, 0 );

    EDGE_V_CHECK_GT( m_nQuad4r, 0 );

    return HEX8R;
  }
  else if( m_nTria3 > 0 ) {
    EDGE_V_CHECK_EQ( m_nQuad4r, 0 );

    return TRIA3;
  }
  else if( m_nQuad4r > 0 ) {
    return QUAD4R;
  }
  else {
    return LINE;
  }
}

std::size_t edge_v::io::Moab::nEnsByDis( unsigned short i_nDis ) const {
  int l_nEns = std::numeric_limits< int >::max();
  moab::ErrorCode l_err = m_moab->get_number_entities_by_dimension( m_root,
                                                                    i_nDis,
                                                                    l_nEns );
  EDGE_V_CHECK_EQ( l_err, moab::MB_SUCCESS );

  return l_nEns;
}

std::size_t edge_v::io::Moab::nEnsByType( t_entityType i_enTy ) const {
  std::size_t l_nEns = 0;

  if(      i_enTy == POINT  ) l_nEns = m_nPoint;
  else if( i_enTy == LINE   ) l_nEns = m_nLine;
  else if( i_enTy == QUAD4R ) l_nEns = m_nQuad4r;
  else if( i_enTy == TRIA3  ) l_nEns = m_nTria3;
  else if( i_enTy == HEX8R  ) l_nEns = m_nHex8r;
  else if( i_enTy == TET4   ) l_nEns = m_nTet4;
  else EDGE_V_LOG_FATAL;

  return l_nEns;
}

void edge_v::io::Moab::getVeCrds( double(*o_veCrds)[3] ) const {
  // get the vertices in the mesh
  std::vector< moab::EntityHandle > l_ves;
  moab::ErrorCode l_err = m_moab->get_entities_by_dimension( m_root,
                                                             0,
                                                             l_ves );
  EDGE_V_CHECK_EQ( l_err, moab::MB_SUCCESS );

  // get their coords
  l_err = m_moab->get_coords( & l_ves[0],
                                l_ves.size(),
                                o_veCrds[0] );
  EDGE_V_CHECK_EQ( l_err, moab::MB_SUCCESS );
}

void edge_v::io::Moab::getEnVe( t_entityType   i_enTy,
                                std::size_t  * o_enVe ) const {
  moab::EntityType l_ty = getMoabType( i_enTy );

  // get mapping from entities to vertices
  std::vector< moab::EntityHandle > l_enCo;
  moab::ErrorCode l_err = m_moab->get_connectivity_by_type( l_ty,
                                                            l_enCo );
  EDGE_V_CHECK_EQ( l_err, moab::MB_SUCCESS );

  // translate to ids
  for( std::size_t l_ve = 0; l_ve < l_enCo.size(); l_ve++ ) {
    o_enVe[l_ve] = m_moab->id_from_handle( l_enCo[l_ve] ) - 1;
  }
}

void edge_v::io::Moab::getFaEl( t_entityType   i_elTy,
                                std::size_t  * o_faEl ) const {
  unsigned short l_nDis = CE_N_DIS( i_elTy );
  std::size_t l_nFas = nEnsByType( CE_T_FA(i_elTy) );
  moab::EntityType l_faTy = getMoabType( CE_T_FA(i_elTy) );

  // init
  for( std::size_t l_fa = 0; l_fa < l_nFas; l_fa++ ) {
    for( std::size_t l_sd = 0; l_sd < 2; l_sd++ ) {
      o_faEl[l_fa*2 + l_sd] = std::numeric_limits< std::size_t >::max();
    }
  }

  // get faces
  std::vector< moab::EntityHandle > l_fas;
  moab::ErrorCode l_err = m_moab->get_entities_by_type( m_root,
                                                        l_faTy,
                                                        l_fas );
  EDGE_V_CHECK_EQ( l_err, moab::MB_SUCCESS );
  EDGE_V_CHECK_EQ( l_fas.size(), l_nFas );

  for( std::size_t l_fa = 0; l_fa < l_nFas; l_fa++ ) {
    // get the face-adjacent elements
    std::vector< moab::EntityHandle > l_els;

    l_err = m_moab->get_adjacencies( &l_fas[l_fa],
                                     1,
                                     l_nDis,
                                     true,
                                     l_els,
                                     moab::Interface::UNION );
    EDGE_V_CHECK_EQ( l_err, moab::MB_SUCCESS );
    EDGE_V_CHECK_GT( l_els.size(), 0 );
    EDGE_V_CHECK_LE( l_els.size(), 2 );

    // assign ids
    for( std::size_t l_si = 0; l_si < l_els.size(); l_si++ ) {
      o_faEl[l_fa*2 + l_si] = m_moab->id_from_handle( l_els[l_si] ) - 1;
    }
    EDGE_V_CHECK_LT( o_faEl[l_fa*2 + 0], o_faEl[l_fa*2 + 1] );
  }
}

void edge_v::io::Moab::getElFa( t_entityType   i_elTy,
                                std::size_t  * o_elFa ) const {
  unsigned short l_nDis = CE_N_DIS( i_elTy );
  unsigned short l_nElFas = CE_N_FAS( i_elTy );
  std::size_t l_nEls = nEnsByType( i_elTy );

  // get elements
  moab::EntityType l_elTy = getMoabType( i_elTy );
  std::vector< moab::EntityHandle > l_els;
  moab::ErrorCode l_err = m_moab->get_entities_by_type( m_root,
                                                        l_elTy,
                                                        l_els );
  EDGE_V_CHECK_EQ( l_err, moab::MB_SUCCESS );
  EDGE_V_CHECK_EQ( l_els.size(), l_nEls );

  for( std::size_t l_el = 0; l_el < l_nEls; l_el++ ) {
    // get the element's adjacent faces
    std::vector< moab::EntityHandle > l_fas;

    l_err = m_moab->get_adjacencies( &l_els[l_el],
                                     1,
                                     l_nDis-1,
                                     true,
                                     l_fas,
                                     moab::Interface::UNION );
    EDGE_V_CHECK_EQ( l_err, moab::MB_SUCCESS );
    EDGE_V_CHECK_EQ( l_fas.size(), l_nElFas );

    // assign ids
    for( unsigned short l_fa = 0; l_fa < l_nElFas; l_fa++ ) {
      o_elFa[l_el*l_nElFas + l_fa] = m_moab->id_from_handle( l_fas[l_fa] ) - 1;
    }
  }
}

void edge_v::io::Moab::getElFaEl( t_entityType   i_elTy,
                                  std::size_t  * o_elFaEl ) const {
#ifdef PP_USE_MPI
EDGE_V_LOG_FATAL; // TODO: MOAB had isssues with non-vertex bridges in MPI-settings
#endif
  unsigned short l_nDis = CE_N_DIS( i_elTy );
  unsigned short l_nElFas = CE_N_FAS( i_elTy );
  moab::EntityType l_ty = getMoabType( i_elTy );

  // get elements
  std::vector< moab::EntityHandle > l_els;
  moab::ErrorCode l_err = m_moab->get_entities_by_type( m_root,
                                                        l_ty,
                                                        l_els );
  EDGE_V_CHECK_EQ( l_err, moab::MB_SUCCESS );

  // iterate over all entities, determine faces and respective adjacent elements
  for( std::size_t l_el = 0; l_el < l_els.size(); l_el++ ) {
    // element's faces
    std::vector< moab::EntityHandle > l_fas;
    l_err = m_moab->get_adjacencies( &l_els[l_el],
                                     1,
                                     l_nDis-1,
                                     true,
                                     l_fas,
                                     moab::Interface::UNION );
    EDGE_V_CHECK_EQ( l_err, moab::MB_SUCCESS );
    EDGE_V_CHECK_EQ( l_fas.size(), l_nElFas );

    // face adjacent elements
    for( unsigned short l_fa = 0; l_fa < l_nElFas; l_fa++ ) {
      std::vector< moab::EntityHandle > l_elFaEl;
      l_err = m_moab->get_adjacencies( &l_fas[l_fa],
                                       1,
                                       l_nDis,
                                       true,
                                       l_elFaEl,
                                       moab::Interface::UNION );
      EDGE_V_CHECK_EQ( l_err, moab::MB_SUCCESS );

      // ensure at least one adjacent element (size includes element itself)
      EDGE_V_CHECK( l_elFaEl.size() == 1 || l_elFaEl.size() == 2 );

      // init adjacency info
      o_elFaEl[l_el * l_nElFas + l_fa] = std::numeric_limits< std::size_t >::max();

      // store info obtained from MOAB
      for( std::size_t l_ad = 0; l_ad < l_elFaEl.size(); l_ad++ ) {
        if( l_elFaEl[l_ad] != l_els[l_el] ) {
          o_elFaEl[l_el * l_nElFas + l_fa] = m_moab->id_from_handle( l_elFaEl[l_ad] ) - 1;
        }
      }
    }
  }
}

void edge_v::io::Moab::getEnMsInter( std::vector< moab::EntityHandle > const & i_ens,
                                     std::vector< moab::EntityHandle > const & i_ms,
                                     std::vector< moab::EntityHandle >       & o_ms ) {
  o_ms.resize(0);
  for( std::size_t l_id = 0; l_id < i_ms.size(); l_id++ ) {
    // create temporary meshset from the entities
    moab::EntityHandle l_msEn;
    moab::ErrorCode l_err = m_moab->create_meshset( moab::MESHSET_SET,
                                                    l_msEn );
    EDGE_V_CHECK_EQ( l_err, moab::MB_SUCCESS );

    l_err = m_moab->add_entities( l_msEn,
                                  i_ens.data(),
                                  i_ens.size() );
    EDGE_V_CHECK_EQ( l_err, moab::MB_SUCCESS );

    // intersect with the input
    l_err = m_moab->intersect_meshset( l_msEn,
                                       i_ms[l_id] );
    EDGE_V_CHECK_EQ( l_err, moab::MB_SUCCESS );

    // store result
    o_ms.push_back( l_msEn );
  }
}

void edge_v::io::Moab::setTags( std::vector< moab::Tag >          const & i_tags,
                                std::vector< moab::EntityHandle > const & i_ens,
                                std::vector< moab::EntityHandle >       & io_ens ) {
  EDGE_V_CHECK_EQ( i_ens.size(), io_ens.size() );

  // iterate over tags
  for( std::size_t l_ta = 0; l_ta < i_tags.size(); l_ta++ ) {
    // get size in bytes
    int l_size = 0;
    moab::ErrorCode l_err = m_moab->tag_get_bytes( i_tags[l_ta],
                                                   l_size );
    EDGE_V_CHECK_EQ( l_err, moab::MB_SUCCESS );
    EDGE_V_CHECK_GT( l_size, 0 );

    int l_length = 0;
    l_err = m_moab->tag_get_length( i_tags[l_ta],
                                    l_length );
    EDGE_V_CHECK_EQ( l_err, moab::MB_SUCCESS );
    EDGE_V_CHECK_GT( l_length, 0 );

    // store tag data in buffer
    unsigned char *l_buf = new unsigned char[ std::size_t(l_size) * l_length * i_ens.size() ];
    l_err = m_moab->tag_get_data( i_tags[l_ta],
                                  i_ens.data(),
                                  i_ens.size(),
                                  l_buf );
    EDGE_V_CHECK_EQ( l_err, moab::MB_SUCCESS );

    // set tag data
    l_err = m_moab->tag_set_data( i_tags[l_ta],
                                  io_ens.data(),
                                  io_ens.size(),
                                  l_buf );
    EDGE_V_CHECK_EQ( l_err, moab::MB_SUCCESS );

    // free buffer
    delete[] l_buf;
  }
}

void edge_v::io::Moab::writeMesh( std::string const & i_pathToMesh ) {
  moab::ErrorCode l_err = m_moab->write_file( i_pathToMesh.c_str() );
  EDGE_V_CHECK_EQ( l_err, moab::MB_SUCCESS );
}

void edge_v::io::Moab::writeMesh( std::size_t         i_first,
                                  std::size_t         i_nEls,
                                  std::string const & i_pathToMesh ) {
  // get the elements
  std::vector< moab::EntityHandle > l_els;
  moab::ErrorCode l_err = m_moab->get_entities_by_type( m_root,
                                                        getMoabType( getElType() ),
                                                        l_els );
  EDGE_V_CHECK_EQ( l_err, moab::MB_SUCCESS );

  // only continue with given subset
  l_els = std::vector< moab::EntityHandle >( &l_els[i_first],
                                             &l_els[i_first+i_nEls] );

  // add requested elements to mesh set
  moab::EntityHandle l_msEl;
  l_err = m_moab->create_meshset( moab::MESHSET_SET,
                                  l_msEl );
  EDGE_V_CHECK_EQ( l_err, moab::MB_SUCCESS );

  l_err = m_moab->add_entities( l_msEl,
                                l_els.data(),
                                l_els.size() );
  EDGE_V_CHECK_EQ( l_err, moab::MB_SUCCESS );

  // preserve material tag
  moab::Tag l_tagMat;
  l_err = m_moab->tag_get_handle( "MATERIAL_SET",
                                  l_tagMat );
  EDGE_V_CHECK_EQ( l_err, moab::MB_SUCCESS );

  // get all sets which define the material tag
  moab::Range l_msMatRa;
  l_err = m_moab->get_entities_by_type_and_tag( m_root,
                                                moab::MBENTITYSET,
                                                &l_tagMat,
                                                nullptr,
                                                1,
                                                l_msMatRa );
  EDGE_V_CHECK_EQ( l_err, moab::MB_SUCCESS );

  // convert to vector
  std::vector< moab::EntityHandle > l_msMat;
  for( std::size_t l_ms = 0; l_ms < l_msMatRa.size(); l_ms++ ) l_msMat.push_back( l_msMatRa[l_ms] );

  // output mesh sets
  std::vector< moab::EntityHandle > l_msOut;
  l_msOut.push_back( l_msEl );

  unsigned short l_nDis = CE_N_DIS( getElType() );
  for( unsigned short l_di = 0; l_di < l_nDis; l_di++ ) {
    std::vector< moab::EntityHandle > l_ens;
    l_err = m_moab->get_adjacencies( l_els.data(),
                                     l_els.size(),
                                     l_di,
                                     true,
                                     l_ens,
                                     moab::Interface::UNION );
    EDGE_V_CHECK_EQ( l_err, moab::MB_SUCCESS );

    // intersect with material sets
    std::vector< moab::EntityHandle > l_msInter;
    getEnMsInter( l_ens,
                  l_msMat,
                  l_msInter );

    // copy matrial tag values
    std::vector< moab::Tag > l_tags;
    l_tags.push_back( l_tagMat );
    setTags( l_tags,
             l_msMat,
             l_msInter );

    l_msOut.insert( l_msOut.end(),
                    l_msInter.begin(),
                    l_msInter.end());
  }

  // write the meshsets and their entities
  l_err = m_moab->write_file( i_pathToMesh.c_str(),
                              0,
                              0,
                              l_msOut.data(),
                              l_msOut.size() );
  EDGE_V_CHECK_EQ( l_err, moab::MB_SUCCESS );

  // remove the created output meshsets
  l_err = m_moab->delete_entities( l_msOut.data(),
                                   l_msOut.size() );
  EDGE_V_CHECK_EQ( l_err, moab::MB_SUCCESS );
}

void edge_v::io::Moab::setEnData( t_entityType        i_enTy,
                                  moab::DataType      i_daTy,
                                  std::string const & i_tagName,
                                  void        const * i_data ) {
  moab::EntityType l_ty = getMoabType(  i_enTy );

  // get the entities by type
  std::vector< moab::EntityHandle > l_ens;
  moab::ErrorCode l_err = m_moab->get_entities_by_type( 0,
                                                        l_ty,
                                                        l_ens );
  EDGE_V_CHECK_EQ( l_err, moab::MB_SUCCESS );

  // create the tag
  moab::Tag l_tag;
  l_err = m_moab->tag_get_handle( i_tagName.c_str(),
                                  1,
                                  i_daTy,
                                  l_tag,
                                  moab::MB_TAG_CREAT|moab::MB_TAG_DENSE );
  EDGE_V_CHECK_EQ( l_err, moab::MB_SUCCESS );

  // store the data
  l_err = m_moab->tag_set_data( l_tag,
                                &l_ens[0],
                                l_ens.size(),
                                i_data );
  EDGE_V_CHECK_EQ( l_err, moab::MB_SUCCESS );
}

void edge_v::io::Moab::setEnData( t_entityType        i_enTy,
                                  std::string const & i_tagName,
                                  float       const * i_data ) {
  // convert to double
  int l_nEns;
  moab::ErrorCode l_err = m_moab->get_number_entities_by_type( m_root,
                                                               getMoabType( i_enTy ),
                                                               l_nEns );
  EDGE_V_CHECK_EQ( l_err, moab::MB_SUCCESS );

  double *l_data = new double[l_nEns];
  convert( l_nEns,
           i_data,
           l_data );

  setEnData( i_enTy,
             moab::MB_TYPE_DOUBLE,
             i_tagName,
             l_data );

  delete[] l_data;
}

void edge_v::io::Moab::setEnData( t_entityType        i_enTy,
                                  std::string const & i_tagName,
                                  double      const * i_data ) {
  setEnData( i_enTy,
             moab::MB_TYPE_DOUBLE,
             i_tagName,
             i_data );
}

void edge_v::io::Moab::setEnData( t_entityType           i_enTy,
                                  std::string    const & i_tagName,
                                  unsigned short const * i_data ) {
  // convert to int
  int l_nEns;
  moab::ErrorCode l_err = m_moab->get_number_entities_by_type( m_root,
                                                               getMoabType( i_enTy ),
                                                               l_nEns );
  EDGE_V_CHECK_EQ( l_err, moab::MB_SUCCESS );

  int *l_data = new int[l_nEns];
  convert( l_nEns,
           i_data,
           l_data );

  setEnData( i_enTy,
             moab::MB_TYPE_INTEGER,
             i_tagName,
             l_data );

  // free memory
  delete[] l_data;
}

void edge_v::io::Moab::setEnData( t_entityType         i_enTy,
                                  std::string  const & i_tagName,
                                  std::size_t  const * i_data ) {
  // convert to int
  int l_nEns;
  moab::ErrorCode l_err = m_moab->get_number_entities_by_type( m_root,
                                                               getMoabType( i_enTy ),
                                                               l_nEns );
  EDGE_V_CHECK_EQ( l_err, moab::MB_SUCCESS );

  int *l_data = new int[l_nEns];
  convert( l_nEns,
           i_data,
           l_data );

  setEnData( i_enTy,
             moab::MB_TYPE_INTEGER,
             i_tagName,
             l_data );

  // free memory
  delete[] l_data;
}

void edge_v::io::Moab::getEnData( t_entityType           i_enTy,
                                  std::string    const & i_tagName,
                                  unsigned short       * o_data ) const {
  // get entities and allocate temporary array
  std::vector< moab::EntityHandle > l_ens;
  moab::ErrorCode l_err = m_moab->get_entities_by_type( m_root,
                                                        getMoabType( i_enTy ),
                                                        l_ens );
  EDGE_V_CHECK_EQ( l_err, moab::MB_SUCCESS );

  int *l_data = new int[l_ens.size()];

  // query mesh for tag
  moab::Tag l_tag;
  l_err = m_moab->tag_get_handle( i_tagName.c_str(),
                                  l_tag );
  EDGE_V_CHECK_EQ( l_err, moab::MB_SUCCESS );

  // get data
  m_moab->tag_get_data( l_tag,
                        l_ens.data(),
                        l_ens.size(),
                        l_data );

  // convert
  convert( l_ens.size(),
           l_data,
           o_data );

  delete[] l_data;
}

void edge_v::io::Moab::getEnDataFromSet( t_entityType         i_enTy,
                                         std::string  const & i_tagName,
                                         int                * o_data ) const {
  // query mesh for tag
  moab::Tag l_tag;
  moab::ErrorCode l_err = m_moab->tag_get_handle( i_tagName.c_str(),
                                                  l_tag );
  EDGE_V_CHECK_EQ( l_err, moab::MB_SUCCESS );

  // get all sets which define the tag
  moab::Range l_sets;
  l_err = m_moab->get_entities_by_type_and_tag( m_root,
                                                moab::MBENTITYSET,
                                                &l_tag,
                                                nullptr,
                                                1,
                                                l_sets );
  EDGE_V_CHECK_EQ( l_err, moab::MB_SUCCESS );

  // get entities
  std::vector< moab::EntityHandle > l_ens;
  l_err = m_moab->get_entities_by_type( m_root,
                                        getMoabType( i_enTy ),
                                                     l_ens );
  EDGE_V_CHECK_EQ( l_err, moab::MB_SUCCESS );

  // lambda which finds the cotaining mesh set for an entity if any
  auto l_findMs = [ this, l_sets, l_ens ]( std::size_t i_en ) {
    std::size_t l_coSe = std::numeric_limits< std::size_t >::max();
    for( std::size_t l_se = 0; l_se < l_sets.size(); l_se++ ) {
      if( m_moab->contains_entities( l_sets[l_se], l_ens.data()+i_en, 1 ) ) {
        l_coSe = l_se;
        break;
      }
    }
    return l_coSe;
  };

  // iterate over the entities
  for( std::size_t l_en = 0; l_en < l_ens.size(); l_en++ ) {
    // determine containg mesh set for the entity
    std::size_t l_coSe = l_findMs(l_en);

    // get tag-data if part of a mesh set
    if( l_coSe < std::numeric_limits< std::size_t >::max() ) {
      EDGE_V_CHECK_LT( l_coSe, l_sets.size() );

      moab::EntityHandle l_setHa = l_sets[l_coSe];
      l_err = m_moab->tag_get_data( l_tag,
                                    &l_setHa,
                                    1,
                                    o_data+l_en );
      EDGE_V_CHECK_EQ( l_err, moab::MB_SUCCESS );
    }
    else {
      o_data[l_en] = std::numeric_limits< int >::max();
    }
  }
}

void edge_v::io::Moab::getTagNames( std::vector< std::string > & o_tagNames ) const {
  // clear output vector
  o_tagNames.clear();

  // obtain tag handles
  std::vector< moab::Tag > l_tags;
  moab::ErrorCode l_err = m_moab->tag_get_tags( l_tags );
  EDGE_V_CHECK_EQ( l_err, moab::MB_SUCCESS );

  // get their names
  for( std::size_t l_ta = 0; l_ta < l_tags.size(); l_ta++ ) {
    std::string l_name;
    m_moab->tag_get_name( l_tags[l_ta],
                          l_name );
    EDGE_V_CHECK_EQ( l_err, moab::MB_SUCCESS );

    o_tagNames.push_back( l_name );
  }
}

bool edge_v::io::Moab::deleteTag( std::string const & i_tagName ) {
  // get tag and exit if it doesn't exist
  moab::Tag l_tag;
  moab::ErrorCode l_err = m_moab->tag_get_handle( i_tagName.c_str(),
                                                  l_tag );
  EDGE_V_CHECK(    l_err == moab::MB_SUCCESS
                || l_err == moab::MB_TAG_NOT_FOUND );
  if( l_err == moab::MB_TAG_NOT_FOUND ) return false;

  // delete tag
  m_moab->tag_delete( l_tag );
  return true;
}

void edge_v::io::Moab::reorder( t_entityType        i_enTy,
                                std::string const & i_tagName,
                                int                 i_skipValue ) {
  // query mesh for tag, which defines the priorities
  moab::Tag l_tagPrio;
  moab::ErrorCode l_err = m_moab->tag_get_handle( i_tagName.c_str(),
                                                  l_tagPrio );
  EDGE_V_CHECK_EQ( l_err, moab::MB_SUCCESS );

  // sanity check of the type
  moab::DataType l_daTy;
  l_err = m_moab->tag_get_data_type( l_tagPrio,
                                     l_daTy );
  EDGE_V_CHECK_EQ( l_err, moab::MB_SUCCESS );
  EDGE_V_CHECK_EQ( l_daTy, moab::DataType::MB_TYPE_INTEGER );

  // use reorder tool to get the permutation of the entities
  moab::Tag l_tagPerm;
  moab::EntityHandle l_zero = 0;
  l_err = m_moab->tag_get_handle( 0,
                                  1,
                                  moab::DataType::MB_TYPE_HANDLE,
                                  l_tagPerm,
                                  moab::MB_TAG_DENSE|moab::MB_TAG_CREAT|moab::MB_TAG_EXCL,
                                  &l_zero );
  EDGE_V_CHECK_EQ( l_err, moab::MB_SUCCESS );

  moab::ReorderTool l_tool( (moab::Core*) m_moab );
  l_err = l_tool.handle_order_from_int_tag( getMoabType(i_enTy),
                                            CE_N_VES(i_enTy),
                                            l_tagPrio,
                                            i_skipValue,
                                            l_tagPerm );
  EDGE_V_CHECK_EQ( l_err, moab::MB_SUCCESS );

  // perform the actual reordering
  l_err = l_tool.reorder_entities( l_tagPerm );
  EDGE_V_CHECK_EQ( l_err, moab::MB_SUCCESS );

  // delete the permutation tag
  l_err = m_moab->tag_delete( l_tagPerm );
  EDGE_V_CHECK_EQ( l_err, moab::MB_SUCCESS );
}
