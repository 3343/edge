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
 * MOAB interface.
 **/
#include "Moab.h"
#include <moab/ReorderTool.hpp>

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

edge_v::io::Moab::Moab( std::string const & i_pathToMesh ) {
  // init moab
  m_moab = new moab::Core;
  moab::ErrorCode l_err = m_moab->load_file( i_pathToMesh.c_str() );
  EDGE_V_CHECK_EQ( l_err, moab::MB_SUCCESS );

  m_root = m_moab->get_root_set();

  // get mesh dimension
  int l_nDis = 0;
  l_err = m_moab->get_dimension( l_nDis );
  EDGE_V_CHECK_EQ( l_err, moab::MB_SUCCESS );
  EDGE_V_CHECK_GT( l_nDis, 0 );

  // get elements
  moab::Range l_els;
  l_err = m_moab->get_entities_by_dimension( m_root,
                                             l_nDis,
                                             l_els );
  EDGE_V_CHECK_EQ( l_err, moab::MB_SUCCESS );

  // create all other entities
  for( unsigned short l_di = 1; l_di < l_nDis; l_di++ ) {
    moab::Range l_ens;
    l_err = m_moab->get_adjacencies( l_els, l_di, true, l_ens, moab::Interface::UNION );
    EDGE_V_CHECK_EQ( l_err, moab::MB_SUCCESS );
  }

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
}

edge_v::t_entityType edge_v::io::Moab::getElType() const {
  EDGE_V_CHECK_GT( m_nLine, 0 );

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

void edge_v::io::Moab::getEnFaEn( t_entityType   i_enTy,
                                  std::size_t  * o_enFaEn ) const {
#ifdef PP_USE_MPI
EDGE_V_LOG_FATAL; // TODO: MOAB had isssues with non-vertex bridges in MPI-settings
#endif
  unsigned short l_nDis = CE_N_DIS( i_enTy );
  unsigned short l_nEnFas = CE_N_FAS( i_enTy );
  moab::EntityType l_ty = getMoabType( i_enTy );

  // get entities of the input type
  std::vector< moab::EntityHandle > l_ens;
  moab::ErrorCode l_err = m_moab->get_entities_by_type( m_root,
                                                        l_ty,
                                                        l_ens );
  EDGE_V_CHECK_EQ( l_err, moab::MB_SUCCESS );

  // iterate over all entities, determine faces and respective adjacent elements
  for( std::size_t l_en = 0; l_en < l_ens.size(); l_en++ ) {
    // element's faces
    std::vector< moab::EntityHandle > l_fas;
    l_err = m_moab->get_adjacencies( &l_ens[l_en],
                                     1,
                                     l_nDis-1,
                                     true,
                                     l_fas,
                                     moab::Interface::UNION );
    EDGE_V_CHECK_EQ( l_err, moab::MB_SUCCESS );

    // face adjacent elements
    std::vector< moab::EntityHandle > l_enFaEn;
    l_err = m_moab->get_adjacencies( &l_fas[0],
                                     l_fas.size(),
                                     l_nDis,
                                     true,
                                     l_enFaEn,
                                     moab::Interface::UNION );
    EDGE_V_CHECK_EQ( l_err, moab::MB_SUCCESS );

    // ensure at least one adjacent element (size includes element itself)
    EDGE_V_CHECK_GT( l_enFaEn.size(), 1 );

    // init adjacency info
    for( unsigned short l_fa = 0; l_fa < l_nEnFas; l_fa++ ) {
      o_enFaEn[l_en * l_nEnFas + l_fa] = std::numeric_limits< std::size_t >::max();
    }

    // store info obtained from MOAB
    unsigned short l_fa = 0;
    for( std::size_t l_ad = 0; l_ad < l_enFaEn.size(); l_ad++ ) {
      if( l_enFaEn[l_ad] != l_ens[l_en] ) {
        o_enFaEn[l_en * l_nEnFas + l_fa] = m_moab->id_from_handle( l_enFaEn[l_ad] ) - 1;
        l_fa++;  
      }
    }
  }
}

void edge_v::io::Moab::writeMesh( std::string const & i_pathToMesh ) {
  moab::ErrorCode l_err = m_moab->write_file( i_pathToMesh.c_str() );
  EDGE_V_CHECK_EQ( l_err, moab::MB_SUCCESS );
}

void edge_v::io::Moab::setGlobalData( moab::DataType         i_daTy,
                                      std::string    const & i_tagName,
                                      std::size_t            i_nValues,
                                      void           const * i_data ) {
  // create the tag
  moab::Tag l_tag;
  moab::ErrorCode l_err = m_moab->tag_get_handle( i_tagName.c_str(),
                                                  i_nValues,
                                                  i_daTy,
                                                  l_tag,
                                                  moab::MB_TAG_CREAT|moab::MB_TAG_SPARSE );
  EDGE_V_CHECK_EQ( l_err, moab::MB_SUCCESS );

  // store the data
  l_err = m_moab->tag_set_data(  l_tag,
                                &m_root,
                                 1,
                                 i_data );
  EDGE_V_CHECK_EQ( l_err, moab::MB_SUCCESS );
}

void edge_v::io::Moab::setGlobalData( std::string const & i_tagName,
                                      std::size_t         i_nValues,
                                      double      const * i_data ) {
  setGlobalData( moab::DataType::MB_TYPE_DOUBLE,
                 i_tagName,
                 i_nValues,
                 i_data );
}

std::size_t edge_v::io::Moab::getGlobalDataSize( std::string const & i_tagName ) const {
  // get the tag
  moab::Tag l_tag;
  moab::ErrorCode l_err = m_moab->tag_get_handle( i_tagName.c_str(),
                                                  l_tag );
  EDGE_V_CHECK_EQ( l_err, moab::MB_SUCCESS );

  // get the size
  int l_size = 0;
  l_err = m_moab->tag_get_length( l_tag,
                                  l_size );
  EDGE_V_CHECK_EQ( l_err, moab::MB_SUCCESS );

  return l_size;
}

void edge_v::io::Moab::getGlobalData( std::string const & i_tagName,
                                      double            * o_data ) const {
  // get the tag
  moab::Tag l_tag;
  moab::ErrorCode l_err = m_moab->tag_get_handle( i_tagName.c_str(),
                                                  l_tag );
  EDGE_V_CHECK_EQ( l_err, moab::MB_SUCCESS );

  // check the type
  moab::DataType l_daTy;
  l_err = m_moab->tag_get_data_type( l_tag,
                                     l_daTy );
  EDGE_V_CHECK_EQ( l_err, moab::MB_SUCCESS );
  EDGE_V_CHECK_EQ( l_daTy, moab::DataType::MB_TYPE_DOUBLE );

  // get data
  l_err = m_moab->tag_get_data(  l_tag,
                                &m_root,
                                 1,
                                 o_data );
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
#ifdef PP_USE_OMP
#pragma omp parallel for
#endif
  for( int l_en = 0; l_en < l_nEns; l_en++ )
    l_data[l_en] = i_data[l_en];

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
#ifdef PP_USE_OMP
#pragma omp parallel for
#endif
  for( std::size_t l_en = 0; l_en < l_ens.size(); l_en++ ) {
    o_data[l_en] = l_data[l_en];
  }

  delete[] l_data;
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