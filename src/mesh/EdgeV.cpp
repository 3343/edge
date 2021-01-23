/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2020, Friedrich Schiller University Jena
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
 * Mesh-interface using EDGE-V.
 **/
#include "io/logging.h"
#include "parallel/global.h"

#include "EdgeV.h"
#include "../data/EntityLayout.h"


void edge::mesh::EdgeV::setElLayout( unsigned short         i_nTgs,
                                     std::size_t    const * i_nTgElsIn,
                                     std::size_t    const * i_nTgElsSe,
                                     t_enLayout           & o_elLay ) {
  // assign sizes
  o_elLay.timeGroups.resize( 0 );
  o_elLay.timeGroups.resize( i_nTgs );

  for( unsigned short l_tg = 0; l_tg < i_nTgs; l_tg++ ) {
    o_elLay.timeGroups[l_tg].inner.size = i_nTgElsIn[l_tg];

    // TODO: dummy entries for send and receive, will be removed
    o_elLay.timeGroups[l_tg].send.resize( 1 );
    o_elLay.timeGroups[l_tg].send[0].size = i_nTgElsSe[l_tg];
    o_elLay.timeGroups[l_tg].receive.resize( 1 );
    o_elLay.timeGroups[l_tg].receive[0].size = 0;
  }

  // derive the rest
  data::EntityLayout::sizesToLayout( o_elLay );
}


void edge::mesh::EdgeV::setLtsTypes( std::size_t            i_nEls,
                                     unsigned short         i_nElFas,
                                     std::size_t    const * i_elFaEl,
                                     unsigned short         i_nTgs,
                                     std::size_t    const * i_nTgElsIn,
                                     std::size_t    const * i_nTgElsSe,
                                     unsigned short const * i_sendFa,
                                     std::size_t    const * i_sendEl,
                                     std::size_t    const * i_commStruct,
                                     long long              i_elEq,
                                     long long              i_elLt,
                                     long long              i_elGt,
                                     long long      const * i_adEq,
                                     long long      const * i_adLt,
                                     long long      const * i_adGt,
                                     long long            * o_spTys ) {
  // define lambda which finds the time group of an element
  auto l_getTg = [ i_nTgs, i_nTgElsIn, i_nTgElsSe ]( std::size_t i_el ) {
    std::size_t l_first = 0;

    // check inner elements
    for( unsigned short l_tg = 0; l_tg < i_nTgs; l_tg++ ) {
      if( i_el >= l_first && i_el < l_first + i_nTgElsIn[l_tg] ) {
        return l_tg;
      }
      l_first += i_nTgElsIn[l_tg];
    }

    // check send elements
    for( unsigned short l_tg = 0; l_tg < i_nTgs; l_tg++ ) {
      if( i_el >= l_first && i_el < l_first + i_nTgElsSe[l_tg] ) {
        return l_tg;
      }
      l_first += i_nTgElsSe[l_tg];
    }

    // fail if we didn't find a time group
    EDGE_LOG_FATAL;
    return std::numeric_limits< unsigned short >::max();
  };

  // attempts to find the time group of an MPI adjacent element
  auto l_getTgComm = [ i_sendFa, i_sendEl, i_commStruct ]( std::size_t    i_el,
                                                           unsigned short i_fa ) {
    unsigned short l_tg = std::numeric_limits< unsigned short >::max();

    if( i_sendFa == nullptr ) {
      return l_tg;
    }

    std::size_t l_nChs = i_commStruct[0];
    std::size_t l_first = 0;
    for( std::size_t l_ch = 0; l_ch < l_nChs; l_ch++ ) {
      std::size_t l_tgAd  = i_commStruct[1 + l_ch*4 + 2];
      std::size_t l_nSeRe = i_commStruct[1 + l_ch*4 + 3];

      for( std::size_t l_co = l_first; l_co < l_first+l_nSeRe; l_co++ ) {
        if( i_fa == i_sendFa[l_co] && i_el == i_sendEl[l_co] ) {
          l_tg = l_tgAd;
          return l_tg;
        }
      }
      l_first += l_nSeRe;
    }

    return l_tg;
  };

  // iterate over the elements
  for( std::size_t l_el = 0; l_el < i_nEls; l_el++ ) {
    unsigned short l_tgEl = l_getTg( l_el );

    bool l_elEq = false;
    bool l_elLt = false;
    bool l_elGt = false;

    // iterate over face-neighbors
    for( unsigned short l_fa = 0; l_fa < i_nElFas; l_fa++ ) {
      std::size_t l_ad = i_elFaEl[l_el*i_nElFas + l_fa];

      unsigned short l_tgAd = std::numeric_limits< unsigned short >::max();

      if( l_ad == std::numeric_limits< std::size_t >::max() ) {
        // try to derive time group for MPI boundaries
        l_tgAd = l_getTgComm( l_el,
                              l_fa );

        // set GTS at domain boundaries
        if( l_tgAd == std::numeric_limits< unsigned short >::max() ) {
          l_tgAd = l_tgEl;
        }
      }
      else {
        l_tgAd = l_getTg( l_ad );
      }

      // set the LTS relations w.r.t. this face-adjacent element
      if( l_tgEl == l_tgAd ) {
        o_spTys[l_el] |= i_adEq[l_fa];
        l_elEq = true;
      }
      else if( l_tgEl < l_tgAd ) {
        EDGE_CHECK_EQ( l_tgEl+1, l_tgAd );
        o_spTys[l_el] |= i_adLt[l_fa];
        l_elLt = true;
      }
      else if( l_tgEl > l_tgAd ) {
        EDGE_CHECK_EQ( l_tgEl, l_tgAd+1 );
        o_spTys[l_el] |= i_adGt[l_fa];
        l_elGt = true;
      }
      else EDGE_LOG_FATAL;
    }

    // set element's LTS configuration
    if( l_elEq ) {
      o_spTys[l_el] |= i_elEq;
    }
    if( l_elLt ) {
      o_spTys[l_el] |= i_elLt;
    }
    if( l_elGt ) {
      o_spTys[l_el] |= i_elGt;
    }
  }
}

edge::mesh::EdgeV::EdgeV( std::string const & i_pathToMesh,
                          std::string const & i_pathToSupplement,
                          int                 i_periodic ) {
  m_gmsh = new edge_v::io::Gmsh();
  m_gmsh->open( i_pathToMesh );
  m_gmsh->readMesh();
  m_mesh = new edge_v::mesh::Mesh( *m_gmsh,
                                   i_periodic );
  m_hdf = new edge_v::io::Hdf5( i_pathToSupplement );

  // check if the mesh is EDGE-V annotated for LTS
  unsigned short l_nLtsTags = 0;
  if( m_hdf->exists( "/edge_v/n_time_group_elements_inner" ) ) l_nLtsTags++;
  if( m_hdf->exists( "/edge_v/n_time_group_elements_send"  ) ) l_nLtsTags++;
  if( m_hdf->exists( "/edge_v/relative_time_steps"         ) ) l_nLtsTags++;
  EDGE_CHECK( l_nLtsTags == 0 || l_nLtsTags == 3 );

  // get the LTS info
  m_nTgs = 1;
  if( l_nLtsTags > 0 ) {
    m_nTgs = m_hdf->nVas( "/edge_v/n_time_group_elements_inner" );
  }

  m_nTgElsIn = new std::size_t[ m_nTgs ];
  m_nTgElsSe = new std::size_t[ m_nTgs ];
  if( l_nLtsTags > 0 ) {
    m_hdf->get( "/edge_v/n_time_group_elements_inner",
                m_nTgElsIn );
    m_hdf->get( "/edge_v/n_time_group_elements_send",
                m_nTgElsSe );
  }
  else {
    m_nTgElsIn[0] = m_mesh->nEls();
    m_nTgElsSe[0] = 0;
  }

  // set number of inner and send elements
  m_nElsIn = m_nElsSe = 0;
  for( unsigned short l_tg = 0; l_tg < m_nTgs; l_tg++ ) {
    m_nElsIn += m_nTgElsIn[l_tg];
    m_nElsSe += m_nTgElsSe[l_tg];
  }

  // allocate memory for relative time steps and init with GTS
  m_relDt = new double[m_nTgs+1];
  m_relDt[0] = 1;
  m_relDt[1] = std::numeric_limits< double >::max();

  if( m_nTgs > 1 ) {
    // get relative time steps and check for rate-2
    m_hdf->get( "/edge_v/relative_time_steps",
                m_relDt );

    for( unsigned short l_tg = 0; l_tg < m_nTgs-1; l_tg++ ) {
      double l_rate = m_relDt[l_tg+1] / m_relDt[l_tg];
      EDGE_CHECK_LT( std::abs(l_rate-2.0), 1E-5 );
    }
  }

  setElLayout( m_nTgs,
               m_nTgElsIn,
               m_nTgElsSe,
               m_elLay );

  // check if the mesh is annotated for MPI
  unsigned short l_nMpiTags = 0;
  if( m_hdf->exists( "/edge_v/communication_structure" ) ) l_nMpiTags++;
  if( m_hdf->exists( "/edge_v/send_el"                 ) ) l_nMpiTags++;
  if( m_hdf->exists( "/edge_v/send_fa"                 ) ) l_nMpiTags++;
  if( m_hdf->exists( "/edge_v/recv_el"                 ) ) l_nMpiTags++;
  if( m_hdf->exists( "/edge_v/recv_fa"                 ) ) l_nMpiTags++;
  if( m_hdf->exists( "/edge_v/send_vertex_ids"         ) ) l_nMpiTags++;
  if( m_hdf->exists( "/edge_v/send_face_ids"           ) ) l_nMpiTags++;

  // check number of mpi tags
  if( parallel::g_nRanks > 1 ) {
    EDGE_CHECK_EQ( l_nMpiTags, 7 );
  }
  else {
    EDGE_CHECK_EQ( l_nMpiTags, 0 );
  }

  if( l_nMpiTags == 7 ) {
    // get communication structure
    std::size_t l_commSize = m_hdf->nVas( "/edge_v/communication_structure" );
    EDGE_CHECK_EQ( l_commSize%4, 1 );
    m_commStruct = new std::size_t[ l_commSize ];

    m_hdf->get( "/edge_v/communication_structure",
                m_commStruct );

    // check ranks in comm structure
    std::size_t l_nChs = m_commStruct[0];
    EDGE_CHECK_GT( l_nChs, 0 );
    for( std::size_t l_ch = 0; l_ch < l_nChs; l_ch++ ) {
      std::size_t l_raAd = m_commStruct[1 + l_ch*4 + 1];
      EDGE_CHECK_NE( parallel::g_rank,   (int) l_raAd );
      EDGE_CHECK_GT( parallel::g_nRanks, (int) l_raAd );
    }

    // get communicating faces
                m_nCommElFa  = m_hdf->nVas( "/edge_v/send_fa" );
    std::size_t l_nSendEls   = m_hdf->nVas( "/edge_v/send_el" );
    std::size_t l_nRecvFas   = m_hdf->nVas( "/edge_v/recv_fa" );
    std::size_t l_nRecvEls   = m_hdf->nVas( "/edge_v/recv_el" );
    std::size_t l_nSeVeIdsAd = m_hdf->nVas( "/edge_v/send_vertex_ids" );
    std::size_t l_nSeFaIdsAd = m_hdf->nVas( "/edge_v/send_face_ids" );

    EDGE_CHECK_EQ( m_nCommElFa, l_nSendEls );
    EDGE_CHECK_EQ( m_nCommElFa, l_nRecvFas );
    EDGE_CHECK_EQ( m_nCommElFa, l_nRecvEls );
    EDGE_CHECK_EQ( m_nCommElFa, l_nSeVeIdsAd );
    EDGE_CHECK_EQ( m_nCommElFa, l_nSeFaIdsAd );

    m_sendFa      = new unsigned short[ m_nCommElFa ];
    m_sendEl      = new std::size_t[ m_nCommElFa ];
    m_recvFa      = new unsigned short[ m_nCommElFa ];
    m_recvEl      = new std::size_t[ m_nCommElFa ];
    m_sendVeIdsAd = new unsigned short[ m_nCommElFa ];
    m_sendFaIdsAd = new unsigned short[ m_nCommElFa ];

    m_hdf->get( "/edge_v/send_fa",
                m_sendFa );
    m_hdf->get( "/edge_v/send_el",
                m_sendEl );
    m_hdf->get( "/edge_v/recv_fa",
                m_recvFa );
    m_hdf->get( "/edge_v/recv_el",
                m_recvEl );
    m_hdf->get( "/edge_v/send_vertex_ids",
                m_sendVeIdsAd );
    m_hdf->get( "/edge_v/send_face_ids",
                m_sendFaIdsAd );
  }
}

edge::mesh::EdgeV::~EdgeV() {
  if( m_gmsh        != nullptr ) delete   m_gmsh;
  if( m_mesh        != nullptr ) delete   m_mesh;
  if( m_hdf         != nullptr ) delete   m_hdf;
  if( m_nTgElsIn    != nullptr ) delete[] m_nTgElsIn;
  if( m_nTgElsSe    != nullptr ) delete[] m_nTgElsSe;
  if( m_relDt       != nullptr ) delete[] m_relDt;
  if( m_commStruct  != nullptr ) delete[] m_commStruct;
  if( m_sendFa      != nullptr ) delete[] m_sendFa;
  if( m_sendEl      != nullptr ) delete[] m_sendEl;
  if( m_recvFa      != nullptr ) delete[] m_recvFa;
  if( m_recvEl      != nullptr ) delete[] m_recvEl;
  if( m_sendVeIdsAd != nullptr ) delete[] m_sendVeIdsAd;
  if( m_sendFaIdsAd != nullptr ) delete[] m_sendFaIdsAd;
}

void edge::mesh::EdgeV::setLtsTypes( t_elementChars * io_elChars ) const {
  // allocate and init temporary array
  std::size_t l_nEls = m_mesh->nEls();
  long long *l_spTys = new long long[ l_nEls ];
  for( std::size_t l_el = 0; l_el < l_nEls; l_el++ ) l_spTys[l_el] = 0;

  long long l_adEq[6];
  long long l_adLt[6];
  long long l_adGt[6];
  for( unsigned short l_fa = 0; l_fa < 6; l_fa++ ) {
    l_adEq[l_fa] = C_LTS_AD[l_fa][AD_EQ];
    l_adLt[l_fa] = C_LTS_AD[l_fa][AD_LT];
    l_adGt[l_fa] = C_LTS_AD[l_fa][AD_GT];
  }

  // set the LTS types
  edge_v::t_entityType l_elTy = m_mesh->getTypeEl();
  unsigned short l_nElFas = edge_v::CE_N_FAS( l_elTy );

  setLtsTypes( l_nEls,
               l_nElFas,
               m_mesh->getElFaEl(),
               m_nTgs,
               m_nTgElsIn,
               m_nTgElsSe,
               m_sendFa,
               m_sendEl,
               m_commStruct,
               C_LTS_EL[EL_INT_EQ],
               C_LTS_EL[EL_INT_LT],
               C_LTS_EL[EL_INT_GT],
               l_adEq,
               l_adLt,
               l_adGt,
               l_spTys );

  // copy over sparse types
  EDGE_CHECK_EQ( sizeof(long long), sizeof( io_elChars->spType ) );
#ifdef PP_USE_OMP
#pragma omp parallel for
#endif
  for( std::size_t l_el = 0; l_el < l_nEls; l_el++ ) {
    io_elChars[l_el].spType |= l_spTys[l_el];
  }

  delete[] l_spTys;
}

void edge::mesh::EdgeV::setSpTypes( t_vertexChars  * io_veChars,
                                    t_faceChars    * io_faChars,
                                    t_elementChars * io_elChars ) const {
  // vertices
  edge_v::t_sparseType const * l_spTypeVe = m_mesh->getSpTypeVe();

#ifdef PP_USE_OMP
#pragma omp parallel for
#endif
  for( std::size_t l_ve = 0; l_ve < nVes(); l_ve++ ) {
    io_veChars[l_ve].spType |= l_spTypeVe[l_ve];
  }

  // faces
  edge_v::t_sparseType const * l_spTypeFa = m_mesh->getSpTypeFa();

#ifdef PP_USE_OMP
#pragma omp parallel for
#endif
  for( std::size_t l_fa = 0; l_fa < nFas(); l_fa++ ) {
    io_faChars[l_fa].spType |= l_spTypeFa[l_fa];
  }

  // elements
  edge_v::t_sparseType const * l_spTypeEl = m_mesh->getSpTypeEl();

#ifdef PP_USE_OMP
#pragma omp parallel for
#endif
  for( std::size_t l_el = 0; l_el < nEls(); l_el++ ) {
    io_elChars[l_el].spType |= l_spTypeEl[l_el];
  }
}

void edge::mesh::EdgeV::setSeVeFaIdsAd( unsigned short * o_veIdsAd,
                                        unsigned short * o_faIdsAd ) const {
  edge_v::t_entityType l_elTy = m_mesh->getTypeEl();
  unsigned short l_nElFas = edge_v::CE_N_FAS( l_elTy );

  // set ids for inner entities
#ifdef PP_USE_OMP
#pragma omp parallel for
#endif
  for( std::size_t l_el = 0; l_el < nEls(); l_el++ ) {
    for( unsigned short l_fa = 0; l_fa < l_nElFas; l_fa++ ) {
      o_veIdsAd[l_el*l_nElFas + l_fa] = std::numeric_limits< unsigned short >::max();
      o_faIdsAd[l_el*l_nElFas + l_fa] = std::numeric_limits< unsigned short >::max();

      if( getElFaEl()[l_el*l_nElFas + l_fa] != std::numeric_limits< std::size_t >::max() ) {
        m_mesh->getVeIdsAd( 1,
                            0,
                            &l_el,
                            &l_fa,
                            o_veIdsAd+(l_el*l_nElFas + l_fa) );

        m_mesh->getFaIdsAd( 1,
                            0,
                            &l_el,
                            &l_fa,
                            o_faIdsAd+(l_el*l_nElFas + l_fa) );
      }
    }
  }

  // set ids for send entities
  if( m_sendVeIdsAd != nullptr ) {
    for( std::size_t l_co = 0; l_co < m_nCommElFa; l_co++ ) {
      std::size_t l_el = m_sendEl[l_co];
      unsigned short l_fa = m_sendFa[l_co];

      o_veIdsAd[l_el*l_nElFas + l_fa] = m_sendVeIdsAd[l_co];
      o_faIdsAd[l_el*l_nElFas + l_fa] = m_sendFaIdsAd[l_co];
    }
  }
}

std::size_t edge::mesh::EdgeV::nVas( std::string const & i_name ) const {
  // abort if the data set does not exists
  if( m_hdf->exists(i_name) == false ) return std::numeric_limits< std::size_t >::max();

  return m_hdf->nVas(i_name);
}

void edge::mesh::EdgeV::getData( std::string const & i_name,
                                 float             * o_data ) const {
  m_hdf->get( i_name,
              o_data );
}