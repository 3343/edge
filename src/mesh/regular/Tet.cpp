/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2016-2017, Regents of the University of California
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

#include "Tet.h"
#include "Base.h"
#include "../common.hpp"
#include "linalg/Geom.hpp"

int_el edge::mesh::regular::Tet::getNVeOwned() const {
  return (m_nHex[0]+1) * (m_nHex[1]+1) * (m_nHex[2]+1);
}

int_el edge::mesh::regular::Tet::getNVeInner() const {
  // get owned vertices first
  int_el l_nVe = getNVeOwned();

  // remove send-vertices
  if( isMpiBnd(m_mpiNe[0][0]) ) l_nVe -= (m_nHex[1]+1) * (m_nHex[2]+1);
  if( isMpiBnd(m_mpiNe[0][1]) ) l_nVe -= (m_nHex[1]+1) * (m_nHex[2]+1);

  if( isMpiBnd(m_mpiNe[1][0]) ) l_nVe -= (m_nHex[0]+1) * (m_nHex[2]+1);
  if( isMpiBnd(m_mpiNe[1][1]) ) l_nVe -= (m_nHex[0]+1) * (m_nHex[2]+1);

  if( isMpiBnd(m_mpiNe[2][0]) ) l_nVe -= (m_nHex[0]+1) * (m_nHex[1]+1);
  if( isMpiBnd(m_mpiNe[2][0]) ) l_nVe -= (m_nHex[0]+1) * (m_nHex[1]+1);

  return l_nVe;
}

int_el edge::mesh::regular::Tet::getNVeSend() const {
  int_el l_nVe  = getNVeOwned();
         l_nVe -= getNVeInner();

  return l_nVe;
}

int_el edge::mesh::regular::Tet::getNVeRecv() const {
  EDGE_LOG_FATAL << "missing implementation.";
  return 0;
}

int_el edge::mesh::regular::Tet::getNVe() const {
  int_el l_nVe  = getNVeOwned();
         l_nVe += getNVeRecv();

  return l_nVe;
}

int_el edge::mesh::regular::Tet::getNFaSend() const {
  int_el l_nFa = 0;

  if( isMpiBnd(m_mpiNe[0][1]) ) l_nFa += m_nHex[1] * m_nHex[2] * 2;
  if( isMpiBnd(m_mpiNe[1][1]) ) l_nFa += m_nHex[0] * m_nHex[2] * 2;
  if( isMpiBnd(m_mpiNe[2][1]) ) l_nFa += m_nHex[0] * m_nHex[1] * 2;

  return l_nFa;
}

int_el edge::mesh::regular::Tet::getNFaRecv() const {
  int_el l_nFa = 0;

  if( isMpiBnd(m_mpiNe[0][0]) ) l_nFa += m_nHex[1] * m_nHex[2] * 2;
  if( isMpiBnd(m_mpiNe[1][0]) ) l_nFa += m_nHex[0] * m_nHex[2] * 2;
  if( isMpiBnd(m_mpiNe[2][0]) ) l_nFa += m_nHex[0] * m_nHex[1] * 2;

  return l_nFa;
}

int_el edge::mesh::regular::Tet::getNFaInner() const {
  int_el l_nFa = 0;

  // add faces aligned with hex faces
  l_nFa +=  m_nHex[0]    *  m_nHex[1]    * (m_nHex[2]+1) * 2;
  l_nFa +=  m_nHex[0]    * (m_nHex[1]+1) *  m_nHex[2]    * 2;
  l_nFa += (m_nHex[0]+1) *  m_nHex[1]    *  m_nHex[2]    * 2;

  // add inner faces
  l_nFa += m_nHex[0] * m_nHex[1] * m_nHex[2] * 4;

  // check that there's either two or no MPI at periodic boundaries
  if( m_periodic ) {
    for( unsigned short l_di = 0; l_di < 3; l_di++ )
      EDGE_CHECK( m_mpiNe[l_di][0] ==m_mpiNe[l_di][1] );
  }

  // periodic, non-MPI faces are only counted once
  if( m_periodic && isMpiBnd(m_mpiNe[0][0]) == false ) l_nFa -= m_nHex[1] * m_nHex[2] * 2;
  if( m_periodic && isMpiBnd(m_mpiNe[1][0]) == false ) l_nFa -= m_nHex[0] * m_nHex[2] * 2;
  if( m_periodic && isMpiBnd(m_mpiNe[2][0]) == false ) l_nFa -= m_nHex[0] * m_nHex[1] * 2;

  // remove send and receive faces
  l_nFa -= getNFaSend();

  l_nFa -= getNFaRecv();

  return l_nFa;
}

int_el edge::mesh::regular::Tet::getNFa() const {
  return m_faLayout.nEnts;
}

int_el edge::mesh::regular::Tet::getNElSend() const {
  int_el l_nEl = 0;

  if( isMpiBnd(m_mpiNe[0][0]) ) l_nEl += m_nHex[1]*m_nHex[2]*2;
  if( isMpiBnd(m_mpiNe[0][1]) ) l_nEl += m_nHex[1]*m_nHex[2]*2;

  if( isMpiBnd(m_mpiNe[1][0]) ) l_nEl += m_nHex[0]*m_nHex[2]*2;
  if( isMpiBnd(m_mpiNe[1][1]) ) l_nEl += m_nHex[0]*m_nHex[2]*2;

  if( isMpiBnd(m_mpiNe[2][0]) ) l_nEl += m_nHex[0]*m_nHex[1]*2;
  if( isMpiBnd(m_mpiNe[2][1]) ) l_nEl += m_nHex[0]*m_nHex[1]*2;

  return l_nEl;
}

int_el edge::mesh::regular::Tet::getNElRecv() const {
  // equivalent to send elements
  return getNElSend();
}

int_el edge::mesh::regular::Tet::getNElOwned() const {
  EDGE_CHECK( m_elLayout.timeGroups.size() == 1 );
  return  m_elLayout.timeGroups[0].nEntsOwn;
}

int_el edge::mesh::regular::Tet::getNElInner() const {
  EDGE_CHECK( m_elLayout.timeGroups.size() == 1 );
  return  m_elLayout.timeGroups[0].inner.size;
}

int_el edge::mesh::regular::Tet::getNEl() const {
  return  m_elLayout.nEnts;
}

bool edge::mesh::regular::Tet::isRecvEl(       unsigned short i_tet,
                                         const int            i_hxPos[3] ) const {
  unsigned short l_recv = 0;
  for( unsigned short l_di = 0; l_di < 3; l_di ++ )
    if( i_hxPos[l_di] == -1 || i_hxPos[l_di] == (int) m_nHex[l_di] ) l_recv++;

  // no receive elements for owned hexes or corners
  if( l_recv == 0 || l_recv > 1 ) return false;

  unsigned int l_hx = getHxId( i_hxPos[0], i_hxPos[1], i_hxPos[2],
                               m_nHex[0],  m_nHex[1],  m_nHex[2] );
  bool l_type = m_hexes[l_hx].type;

  // y-z MPI-bnd
  if( i_hxPos[0]== -1 ) {
    if(      l_type==false && (i_tet==0 || i_tet==4) ) return isMpiBnd(m_mpiNe[0][0]);
    else if( l_type==true  && (i_tet==1 || i_tet==3) ) return isMpiBnd(m_mpiNe[0][0]);
  }
  else if( i_hxPos[0] == (int) m_nHex[0] ) {
    if(      l_type==false && (i_tet==1 || i_tet==3) ) return isMpiBnd(m_mpiNe[0][1]);
    else if( l_type==true  && (i_tet==0 || i_tet==4) ) return isMpiBnd(m_mpiNe[0][1]);
  }

  // x-z MPI-bnd
  if( i_hxPos[1] == -1 ) {
    if(      l_type==false && (i_tet==1 || i_tet==4) ) return isMpiBnd(m_mpiNe[1][0]);
    else if( l_type==true  && (i_tet==1 || i_tet==4) ) return isMpiBnd(m_mpiNe[1][0]);
  }
  else if( i_hxPos[1] == (int) m_nHex[1] ) {
    if(      l_type==false && (i_tet==0 || i_tet==3) ) return isMpiBnd(m_mpiNe[1][1]);
    else if( l_type==true  && (i_tet==0 || i_tet==3) ) return isMpiBnd(m_mpiNe[1][1]);
  }

  // x-y MPI-bnd
  if( i_hxPos[2] == -1 ) {
    if(      l_type==false && (i_tet==3 || i_tet==4) ) return isMpiBnd(m_mpiNe[2][0]);
    else if( l_type==true  && (i_tet==3 || i_tet==4) ) return isMpiBnd(m_mpiNe[2][0]);
  }
  else if( i_hxPos[2] == (int) m_nHex[2] ) {
    if(      l_type==false && (i_tet==0 || i_tet==1) ) return isMpiBnd(m_mpiNe[2][1]);
    else if( l_type==true  && (i_tet==0 || i_tet==1) ) return isMpiBnd(m_mpiNe[2][1]);
  }

  return false;
}

void edge::mesh::regular::Tet::isSendEl(       unsigned short i_tet,
                                         const int            i_hxPos[3],
                                               bool           o_send[3] ) const {
  // reset output info
  o_send[0] = o_send[1] = o_send[2] = false;

  // no send elements if not owned
  for( unsigned int l_di = 0; l_di < 3; l_di++ )
    if( i_hxPos[l_di] == -1 || i_hxPos[l_di] == (int) m_nHex[l_di] ) return;

  // get hex type
  unsigned int l_hx = getHxId( i_hxPos[0], i_hxPos[1], i_hxPos[2],
                               m_nHex[0],  m_nHex[1],  m_nHex[2] );
  bool l_type = m_hexes[l_hx].type;

  // y-z MPI-bnd
  if( i_hxPos[0] == 0 ) {
    if(      l_type==false && (i_tet==1 || i_tet==3) ) o_send[0] = isMpiBnd(m_mpiNe[0][0]);
    else if( l_type==true  && (i_tet==0 || i_tet==4) ) o_send[0] = isMpiBnd(m_mpiNe[0][0]);
  }
  else if( i_hxPos[0] == (int) m_nHex[0]-1 ) {
    if(      l_type==false && (i_tet==0 || i_tet==4) ) o_send[0] = isMpiBnd(m_mpiNe[0][1]);
    else if( l_type==true  && (i_tet==1 || i_tet==3) ) o_send[0] = isMpiBnd(m_mpiNe[0][1]);
  }

  // x-z MPI-bnd
  if( i_hxPos[1] == 0 ) {
    if(      l_type==false && (i_tet==0 || i_tet==3) ) o_send[1] = isMpiBnd(m_mpiNe[1][0]);
    else if( l_type==true  && (i_tet==0 || i_tet==3) ) o_send[1] = isMpiBnd(m_mpiNe[1][0]);
  }
  else if( i_hxPos[1] == (int) m_nHex[1]-1 ) {
    if(      l_type==false && (i_tet==1 || i_tet==4) ) o_send[1] = isMpiBnd(m_mpiNe[1][1]);
    else if( l_type==true  && (i_tet==1 || i_tet==4) ) o_send[1] = isMpiBnd(m_mpiNe[1][1]);
  }

  // x-y MPI-bnd
  if( i_hxPos[2] == 0 ) {
    if(      l_type==false && (i_tet==0 || i_tet==1) ) o_send[2] = isMpiBnd(m_mpiNe[2][0]);
    else if( l_type==true  && (i_tet==0 || i_tet==1) ) o_send[2] = isMpiBnd(m_mpiNe[2][0]);
  }
  else if( i_hxPos[2] == (int) m_nHex[2]-1 ) {
    if(      l_type==false && (i_tet==3 || i_tet==4) ) o_send[2] = isMpiBnd(m_mpiNe[2][1]);
    else if( l_type==true  && (i_tet==3 || i_tet==4) ) o_send[2] = isMpiBnd(m_mpiNe[2][1]);
  }
}

bool edge::mesh::regular::Tet::isInnerEl(       unsigned short i_tet,
                                          const int            i_hxPos[3] ) const {
  // rule out trivial cases
  if( i_hxPos[0] > 0 && i_hxPos[0] < (int) m_nHex[0]-1 &&
      i_hxPos[1] > 0 && i_hxPos[1] < (int) m_nHex[1]-1 &&
      i_hxPos[2] > 0 && i_hxPos[2] < (int) m_nHex[2]-1 ) return true;
  // if not check for send elements
  else if( i_hxPos[0] >= 0 && i_hxPos[0] < (int) m_nHex[0] &&
           i_hxPos[1] >= 0 && i_hxPos[1] < (int) m_nHex[1] &&
           i_hxPos[2] >= 0 && i_hxPos[2] < (int) m_nHex[2] ) {
    bool l_send[3];
    isSendEl( i_tet, i_hxPos, l_send );
    return !( l_send[0] || l_send[1] || l_send[2] );
  }

  return false;
}

void edge::mesh::regular::Tet::deriveInMapVeAndFa() {
  m_inMap.veMeDa.resize( m_veLayout.nEnts );
  m_inMap.veDaMe.resize( m_veLayout.nEnts );
  // setup dummy info for vertices
  for( int_el l_ve = 0; l_ve < m_veLayout.nEnts; l_ve++ ) {
    m_inMap.veMeDa[l_ve] = m_inMap.veDaMe[l_ve] = l_ve;
  }

  // faces aren't duplicated at all
  m_inMap.faMeDa.resize( m_faLayout.nEnts );
  m_inMap.faDaMe.resize( m_faLayout.nEnts );
  for( int_el l_fa = 0; l_fa < m_faLayout.nEnts; l_fa++ ) {
    m_inMap.faMeDa[l_fa] = m_inMap.faDaMe[l_fa] = l_fa;
  }
}

void edge::mesh::regular::Tet::initSendFa(       unsigned short i_dim,
                                           const int            i_hxPos[3],
                                           const unsigned short i_cIds[6] ) {
  // get hex index
  unsigned int l_hx = getHxId( i_hxPos[0], i_hxPos[1], i_hxPos[2], m_nHex[0], m_nHex[1], m_nHex[2] );

  std::function< void( unsigned short, unsigned short )  > l_updateBnd = [&]( unsigned short i_fa,
                                                                              unsigned short i_cId ) {
    // set the face
    m_hexes[l_hx].fa[i_fa] = m_faLayout.nEnts;

    // update the counter
    m_faLayout.nEnts++;
    m_faLayout.timeGroups[0].nEntsOwn++;
    m_faLayout.timeGroups[0].send[i_cId].size++;
  };

  /**
   * Face setup
   **/
  // check if the hex is owned
  bool l_own = true;
  for( unsigned short l_di = 0; l_di < 3; l_di++ ) {
    l_own = l_own && ( i_hxPos[l_di] >= 0 && i_hxPos[l_di] < (int) m_nHex[l_di] );
  }

  // only consider owned hexes for faces
  if( l_own ) {
    // right MPI-bnd
    if( i_dim==0 && i_hxPos[0] == (int) m_nHex[0]-1 && isMpiBnd(m_mpiNe[0][1]) ) {
      for( unsigned short l_fa =  2; l_fa <  4; l_fa++ ) l_updateBnd( l_fa, i_cIds[1] );
    }
    // back MPI-bnd
    if( i_dim==1 && i_hxPos[1] == (int) m_nHex[1]-1 && isMpiBnd(m_mpiNe[1][1]) ) {
      for( unsigned short l_fa =  6; l_fa <  8; l_fa++ ) l_updateBnd( l_fa, i_cIds[3] );
    }
    // top MPI-bnd
    if( i_dim==2 && i_hxPos[2] == (int) m_nHex[2]-1 && isMpiBnd(m_mpiNe[2][1]) ) {
      for( unsigned short l_fa = 10; l_fa < 12; l_fa++ ) l_updateBnd( l_fa, i_cIds[5] );
    }
  }
}

void edge::mesh::regular::Tet::initSendEl(       unsigned short i_dim,
                                           const int            i_id[3],
                                           const unsigned short i_cIds[6] ) {
  // get hex index
  unsigned int l_hx = getHxId( i_id[0], i_id[1], i_id[2], m_nHex[0], m_nHex[1], m_nHex[2] );

  std::function< void( unsigned short, unsigned short )  > l_updateBnd = [&]( unsigned short i_te,
                                                                              unsigned short i_cId ) {
    m_hexes[l_hx].el[i_te].push_back( m_elLayout.nEnts );

    if( m_hexes[l_hx].el[i_te].size() == 1 ) {
      m_inMap.elMeDa.push_back( m_elLayout.nEnts );
      m_inMap.elDaMe.push_back( m_nMeshEl        );
      m_nMeshEl++;
    }
    else {
      int_el l_daId = m_hexes[l_hx].el[i_te][0];
      m_inMap.elDaMe.push_back( m_inMap.elDaMe[l_daId] );
    }

    m_elLayout.nEnts++;
    m_elLayout.timeGroups[0].nEntsOwn++;
    m_elLayout.timeGroups[0].send[i_cId].size++;
  };

  /**
   * Element setup
   **/
  for( unsigned short l_te = 0; l_te < 5; l_te++ ) {
    bool l_se[3];
    isSendEl( l_te, i_id, l_se );

    // x mpi bnds
    if(      i_dim==0 && i_id[0] == 0                 && l_se[0] ) l_updateBnd( l_te, i_cIds[0] );
    else if( i_dim==0 && i_id[0] == (int) m_nHex[0]-1 && l_se[0] ) l_updateBnd( l_te, i_cIds[1] );

    // y mpi bnds
    if(      i_dim==1 && i_id[1] == 0                 && l_se[1] ) l_updateBnd( l_te, i_cIds[2] );
    else if( i_dim==1 && i_id[1] == (int) m_nHex[1]-1 && l_se[1] ) l_updateBnd( l_te, i_cIds[3] );

    // z mpi bnds
    if(      i_dim==2 && i_id[2] == 0                 && l_se[2] ) l_updateBnd( l_te, i_cIds[4] );
    else if( i_dim==2 && i_id[2] == (int) m_nHex[2]-1 && l_se[2] ) l_updateBnd( l_te, i_cIds[5] );
  }
}

void edge::mesh::regular::Tet::initRecvFa(       unsigned short i_dim,
                                           const int            i_id[3],
                                           const unsigned short i_cIds[6] ) {
  // get hex index
  unsigned int l_hx = getHxId( i_id[0], i_id[1], i_id[2], m_nHex[0], m_nHex[1], m_nHex[2] );

  std::function< void( unsigned short, unsigned short )  > l_updateBnd = [&]( unsigned short i_fa,
                                                                              unsigned short i_cId ) {
        m_hexes[l_hx].fa[i_fa] = m_faLayout.nEnts;

        m_faLayout.nEnts++;
        m_faLayout.timeGroups[0].nEntsNotOwn++;
        m_faLayout.timeGroups[0].receive[i_cId].size++;
  };

  /**
   * Face setup
   **/
  // check if the hex is owned
  bool l_own = true;
  for( unsigned short l_di = 0; l_di < 3; l_di++ ) {
    l_own = l_own && ( i_id[l_di] >= 0 && i_id[l_di] < (int) m_nHex[l_di] );
  }

  // only consider owned hexes for faces
  if( l_own ) {
    // left MPI-bnd
    if( i_dim==0 && i_id[0] == 0 && isMpiBnd(m_mpiNe[0][0]) ) {
      for( unsigned short l_fa = 0; l_fa <  2; l_fa++ ) l_updateBnd( l_fa, i_cIds[0] );
    }
    // front MPI-bnd
    if( i_dim==1 && i_id[1] == 0 && isMpiBnd(m_mpiNe[1][0]) ) {
      for( unsigned short l_fa = 4; l_fa <  6; l_fa++ ) l_updateBnd( l_fa, i_cIds[2] );
    }
    // bottom MPI-bnd
    if( i_dim==2 && i_id[2] == 0 && isMpiBnd(m_mpiNe[2][0]) ) {
      for( unsigned short l_fa = 8; l_fa < 10; l_fa++ ) l_updateBnd( l_fa, i_cIds[4] );
    }
  }
}

void edge::mesh::regular::Tet::initRecvEl(       unsigned short i_dim,
                                           const int            i_id[3],
                                           const unsigned short i_cIds[6] ) {
  // get hex index
  unsigned int l_hx = getHxId( i_id[0], i_id[1], i_id[2], m_nHex[0], m_nHex[1], m_nHex[2] );

  std::function< void( unsigned short, unsigned short )  > l_updateBnd = [&]( unsigned short i_te,
                                                                              unsigned short i_cId ) {
      m_hexes[l_hx].el[i_te].push_back( m_elLayout.nEnts );

      EDGE_CHECK( m_hexes[l_hx].el[i_te].size() == 1 );

      m_inMap.elMeDa.push_back( m_elLayout.nEnts );
      m_inMap.elDaMe.push_back( m_nMeshEl        );
      m_nMeshEl++;

      m_elLayout.nEnts++;
      m_elLayout.timeGroups[0].nEntsNotOwn++;
      m_elLayout.timeGroups[0].receive[i_cId].size++;
  };

  /**
   * Element setup
   **/
  for( unsigned short l_te = 0; l_te < 5; l_te++ ) {
    bool l_re = isRecvEl( l_te, i_id );

    if(      i_dim==0 && i_id[0] == -1              && l_re ) l_updateBnd( l_te, i_cIds[0] );
    else if( i_dim==0 && i_id[0] == (int) m_nHex[0] && l_re ) l_updateBnd( l_te, i_cIds[1] );

    if(      i_dim==1 && i_id[1] == -1              && l_re ) l_updateBnd( l_te, i_cIds[2] );
    else if( i_dim==1 && i_id[1] == (int) m_nHex[1] && l_re ) l_updateBnd( l_te, i_cIds[3] );

    if(      i_dim==2 && i_id[2] == -1              && l_re ) l_updateBnd( l_te, i_cIds[4] );
    else if( i_dim==2 && i_id[2] == (int) m_nHex[2] && l_re ) l_updateBnd( l_te, i_cIds[5] );
  }

}

void edge::mesh::regular::Tet::initFirstPos( t_enLayout &io_enLo ) {
  int_el l_elPos = 0;
  for( unsigned short l_tg = 0; l_tg < io_enLo.timeGroups.size(); l_tg++ ) {
    io_enLo.timeGroups[l_tg].inner.first = l_elPos;
    l_elPos += io_enLo.timeGroups[l_tg].inner.size;

    for( unsigned int l_se = 0; l_se < io_enLo.timeGroups[l_tg].send.size(); l_se++ ) {
      io_enLo.timeGroups[l_tg].send[l_se].first = l_elPos;
      l_elPos += io_enLo.timeGroups[l_tg].send[l_se].size;
    }

    for( unsigned int l_re = 0; l_re < io_enLo.timeGroups[l_tg].receive.size(); l_re++ ) {
      io_enLo.timeGroups[l_tg].receive[l_re].first = l_elPos;
      l_elPos += io_enLo.timeGroups[l_tg].receive[l_re].size;
    }
  }
}

void edge::mesh::regular::Tet::init(       bool         i_type,
                                     const unsigned int i_nHex[3],
                                     const int          i_mpiNe[3][2],
                                     const double       i_corner[3],
                                     const double       i_dX[3] ) {
  for( unsigned short l_di = 0; l_di < 3; l_di++ ) {
    m_nHex[l_di] = i_nHex[l_di];
    for( unsigned short l_si = 0; l_si < 2; l_si++ ) {
      m_mpiNe[l_di][l_si] = i_mpiNe[l_di][l_si];
    }
  }
  m_nMeshEl = 0;

  EDGE_CHECK( i_corner != NULL );
  EDGE_CHECK( i_dX     != NULL );

  for( unsigned short l_di = 0; l_di < 3; l_di++ ) {
    m_corner[l_di] = i_corner[l_di];
    m_dX[l_di]     = i_dX[l_di];
  }

  // setup comm ids
  unsigned short l_cIds[6];
  short l_tmpId = -1;

  for( unsigned short l_di = 0; l_di < 3; l_di ++ ) {
    if( isMpiBnd(m_mpiNe[l_di][0]) ) l_tmpId++;
    l_cIds[l_di*2+0] = l_tmpId;
    if( isMpiBnd(m_mpiNe[l_di][1]) && m_mpiNe[l_di][1] != m_mpiNe[l_di][0] ) l_tmpId++;
    l_cIds[l_di*2+1] = l_tmpId;
  }

  // init data layouts
  std::function< void( t_enLayout& ) > l_initEnLo =
    [&]( t_enLayout &o_lo ) {
      o_lo.nEnts = 0;
      o_lo.timeGroups.clear();
      o_lo.timeGroups.resize( 1 );
      o_lo.timeGroups[0].nEntsOwn    = 0;
      o_lo.timeGroups[0].nEntsNotOwn = 0;
      o_lo.timeGroups[0].inner.first = 0;
      o_lo.timeGroups[0].inner.size  = 0;

      for( unsigned short l_di = 0; l_di < 3; l_di ++ ) {
         if( isMpiBnd(m_mpiNe[l_di][0]) ) {
           o_lo.timeGroups[0].neRanks.push_back( m_mpiNe[l_di][0] );
           o_lo.timeGroups[0].neTgs.push_back( 0 );
         }
         if( isMpiBnd(m_mpiNe[l_di][1]) && m_mpiNe[l_di][1] != m_mpiNe[l_di][0] ) {
           o_lo.timeGroups[0].neRanks.push_back( m_mpiNe[l_di][1] );
           o_lo.timeGroups[0].neTgs.push_back( 0 );
         }
      }

      o_lo.timeGroups[0].send.resize(    o_lo.timeGroups[0].neRanks.size() );
      o_lo.timeGroups[0].receive.resize( o_lo.timeGroups[0].neRanks.size() );

      for( unsigned short l_ne = 0; l_ne < o_lo.timeGroups[0].neRanks.size(); l_ne++ ) {
        o_lo.timeGroups[0].send[l_ne].first    = 0;
        o_lo.timeGroups[0].send[l_ne].size     = 0;
        o_lo.timeGroups[0].receive[l_ne].first = 0;
        o_lo.timeGroups[0].receive[l_ne].size  = 0;
      }
  };

  l_initEnLo( m_veLayout );
  l_initEnLo( m_faLayout );
  l_initEnLo( m_elLayout );

  EDGE_CHECK( m_faLayout.timeGroups[0].send.size() == m_faLayout.timeGroups[0].receive.size() );
  EDGE_CHECK( m_elLayout.timeGroups[0].send.size() == m_faLayout.timeGroups[0].send.size() );
  EDGE_CHECK( m_elLayout.timeGroups[0].send.size() == m_elLayout.timeGroups[0].receive.size() );

  m_veLayout.nEnts = Base::getNVeHex( m_nHex );
  EDGE_CHECK( m_veLayout.timeGroups.size() == 1 );
  m_veLayout.timeGroups[0].nEntsOwn = getNVeOwned();
  // TODO: Add neighboring info if required.

  // init hexes
  m_hexes.clear();
  m_hexes.resize( (m_nHex[0]+2)*(m_nHex[1]+2)*(m_nHex[2]+2) );
  for( std::size_t l_hx = 0; l_hx < m_hexes.size(); l_hx++ ) {
    m_hexes[l_hx].type = false;
    for( unsigned short l_ve = 0; l_ve <  8; l_ve++ ) m_hexes[l_hx].ve[l_ve] = std::numeric_limits<int_el>::max();
    for( unsigned short l_fa = 0; l_fa < 16; l_fa++ ) m_hexes[l_hx].fa[l_fa] = std::numeric_limits<int_el>::max();
  }

  std::function< void( unsigned int, unsigned short )  > l_updateFaInfo = [&]( unsigned int   i_hx,
                                                                               unsigned short i_fa ) {
    m_hexes[i_hx].fa[i_fa] = m_faLayout.nEnts;

    m_faLayout.nEnts++;
    m_faLayout.timeGroups[0].inner.size++;
    m_faLayout.timeGroups[0].nEntsOwn++;
  };

  // iterate over all hexes and set up type and inner ents, and vertices
  for( int l_zh = -1; l_zh < (int) m_nHex[2]+1; l_zh++ ) {
    for( int l_yh = -1; l_yh < (int) m_nHex[1]+1; l_yh++ ) {
      for( int l_xh = -1; l_xh < (int) m_nHex[0]+1; l_xh++ ) {
        // get hex index
        unsigned int l_hx = getHxId( l_xh, l_yh, l_zh, m_nHex[0], m_nHex[1], m_nHex[2] );
        EDGE_CHECK( l_hx < m_hexes.size() );

        m_hexes[l_hx].type = i_type;
        if( std::abs(l_xh)%2 == 1 ) m_hexes[l_hx].type = !m_hexes[l_hx].type;
        if( std::abs(l_yh)%2 == 1 ) m_hexes[l_hx].type = !m_hexes[l_hx].type;
        if( std::abs(l_zh)%2 == 1 ) m_hexes[l_hx].type = !m_hexes[l_hx].type;

        /*
         * Face setup
         */
        // check if the hex is owned
        bool l_xOwn = l_xh >= 0 && l_xh < (int) m_nHex[0];
        bool l_yOwn = l_yh >= 0 && l_yh < (int) m_nHex[1];
        bool l_zOwn = l_zh >= 0 && l_zh < (int) m_nHex[2];

        // compute neighboring indices
        unsigned int l_hxLe = getHxIdBnd( l_xh-1, l_yh,   l_zh,   m_nHex[0], m_nHex[1], m_nHex[2] );
        unsigned int l_hxFr = getHxIdBnd( l_xh,   l_yh-1, l_zh,   m_nHex[0], m_nHex[1], m_nHex[2] );
        unsigned int l_hxBo = getHxIdBnd( l_xh,   l_yh,   l_zh-1, m_nHex[0], m_nHex[1], m_nHex[2] );
        unsigned int l_hxRi = getHxIdBnd( l_xh+1, l_yh,   l_zh,   m_nHex[0], m_nHex[1], m_nHex[2] );
        unsigned int l_hxBa = getHxIdBnd( l_xh,   l_yh+1, l_zh,   m_nHex[0], m_nHex[1], m_nHex[2] );
        unsigned int l_hxTo = getHxIdBnd( l_xh,   l_yh,   l_zh+1, m_nHex[0], m_nHex[1], m_nHex[2] );

        // only consider owned hexes for faces
        if( l_xOwn && l_yOwn && l_zOwn ) {
          // iterate over left faces
          for( unsigned short l_fa = 0; l_fa < 2; l_fa++ ) {
            // left-bnd, non-mpi faces are inner faces and have to be created
            if( l_xh == 0 && !isMpiBnd(m_mpiNe[0][0]) ) {
              l_updateFaInfo( l_hx, l_fa );
            }
            // all other faces are part of the left-neighboring hex
            else if( l_xh > 0 ) {
              m_hexes[l_hx].fa[l_fa] = m_hexes[l_hxLe].fa[l_fa+2];
            }
          }

          // iterate over right faces
          for( unsigned short l_fa = 2; l_fa < 4; l_fa++ ) {
            // create right-inner faces
            if(    (l_xh == (int) m_nHex[0]-1 && !isMpiBnd(m_mpiNe[0][1]) && !m_periodic)
                ||  l_xh <  (int) m_nHex[0]-1 ) l_updateFaInfo( l_hx, l_fa );
            else if( m_periodic && !isMpiBnd(m_mpiNe[0][1]) ) m_hexes[l_hx].fa[l_fa] = m_hexes[l_hxRi].fa[l_fa-2];
          }

          // iterate over front faces
          for( unsigned short l_fa = 4; l_fa < 6; l_fa++ ) {
            // front-bd, non-mpi faces are inner and have to be created
            if( l_yh == 0 && !isMpiBnd(m_mpiNe[1][0]) ) l_updateFaInfo( l_hx, l_fa );
            // all other faces are part of the front-neighbor
            else if( l_yh > 0 ) {
              m_hexes[l_hx].fa[l_fa] = m_hexes[l_hxFr].fa[l_fa+2];
            }
          }

          // iterate over back faces
          for( unsigned short l_fa = 6; l_fa < 8; l_fa++ ) {
            // inner back-faces belong to this hex
            if(    (l_yh == (int) m_nHex[1]-1 && !isMpiBnd(m_mpiNe[1][1]) && !m_periodic)
                ||  l_yh <  (int) m_nHex[1]-1 ) l_updateFaInfo( l_hx, l_fa );
            else if( m_periodic && !isMpiBnd(m_mpiNe[1][1]) ) m_hexes[l_hx].fa[l_fa] = m_hexes[l_hxBa].fa[l_fa-2];
          }

          // iterate over bottom faces
          for( unsigned short l_fa = 8; l_fa < 10; l_fa++ ) {
            // bottom-bnd, non-mpi faces are inner and created
            if( l_zh == 0 && !isMpiBnd(m_mpiNe[2][0]) ) l_updateFaInfo( l_hx, l_fa );
            // all other faces are part of the neighboring hex
            else if( l_zh > 0 ) {
              m_hexes[l_hx].fa[l_fa] = m_hexes[l_hxBo].fa[l_fa+2];
            }
          }

          // iterate over top faces
          for( unsigned short l_fa = 10; l_fa < 12; l_fa++ ) {
            // top inner-faces are always created
            if(    (l_zh == (int) m_nHex[2]-1 && !isMpiBnd(m_mpiNe[2][1]) && !m_periodic)
                ||  l_zh <  (int) m_nHex[2]-1 ) l_updateFaInfo( l_hx, l_fa );
            else if( m_periodic && !isMpiBnd(m_mpiNe[2][1]) ) m_hexes[l_hx].fa[l_fa] = m_hexes[l_hxTo].fa[l_fa-2];
          }

          // iterate over faces inside the hex and create them, at least no shares here..
          for( unsigned int short l_fa = 12; l_fa < 16; l_fa++ ) l_updateFaInfo( l_hx, l_fa );
        }

        /*
         * Element setup
         */
        // iterate over tets
        for( unsigned short l_te = 0; l_te < 5; l_te++ ) {
          int l_hxPos[3];
          l_hxPos[0] = l_xh; l_hxPos[1] = l_yh; l_hxPos[2] = l_zh;
          bool l_in = isInnerEl( l_te, l_hxPos );
          if( l_in ) {
            m_hexes[l_hx].el[l_te].push_back( m_elLayout.nEnts );

            // inner elmements are unique
            EDGE_CHECK( m_elLayout.nEnts == m_nMeshEl );

            m_inMap.elMeDa.push_back( m_elLayout.nEnts );
            m_inMap.elDaMe.push_back( m_nMeshEl        );
            m_nMeshEl++;

            m_elLayout.nEnts++;
            m_elLayout.timeGroups[0].inner.size++;
            m_elLayout.timeGroups[0].nEntsOwn++;
          }
        }

        /**
         * Vertex setup
         **/
        Base::getVesHex( l_xh,
                         l_yh,
                         l_zh,
                         m_nHex[0],
                         m_nHex[1],
                         m_nHex[2],
                         m_hexes[l_hx].ve );
      }
    }
  }

  // set up communication, layer by layer
  int l_ids[3];

  // left
  for( int l_zh = 0; l_zh < (int) m_nHex[2]; l_zh++ ) {
    for( int l_yh = 0; l_yh < (int) m_nHex[1]; l_yh++ ) {
      l_ids[0] = 0; l_ids[1] = l_yh; l_ids[2] = l_zh;
      initSendEl( 0, l_ids, l_cIds );
      initRecvFa( 0, l_ids, l_cIds );
    }
  }
  // right
  for( int l_zh = 0; l_zh < (int) m_nHex[2]; l_zh++ ) {
    for( int l_yh = 0; l_yh < (int) m_nHex[1]; l_yh++ ) {
      l_ids[0] = m_nHex[0]-1; l_ids[1] = l_yh; l_ids[2] = l_zh;
      initSendEl( 0, l_ids, l_cIds );
      initSendFa( 0, l_ids, l_cIds );
    }
  }
  // front
  for( int l_zh = 0; l_zh < (int) m_nHex[2]; l_zh++ ) {
    for( int l_xh = 0; l_xh < (int) m_nHex[0]; l_xh++ ) {
      l_ids[0] = l_xh; l_ids[1] = 0; l_ids[2] = l_zh;
      initSendEl( 1, l_ids, l_cIds );
      initRecvFa( 1, l_ids, l_cIds );
    }
  }
  // back
  for( int l_zh = 0; l_zh < (int) m_nHex[2]; l_zh++ ) {
    for( int l_xh = 0; l_xh < (int) m_nHex[0]; l_xh++ ) {
      l_ids[0] = l_xh; l_ids[1] = m_nHex[1]-1; l_ids[2] = l_zh;
      initSendEl( 1, l_ids, l_cIds );
      initSendFa( 1, l_ids, l_cIds );
    }
  }
  // bottom
  for( int l_yh = 0; l_yh < (int) m_nHex[1]; l_yh++ ) {
    for( int l_xh = 0; l_xh < (int) m_nHex[0]; l_xh++ ) {
      l_ids[0] = l_xh; l_ids[1] = l_yh; l_ids[2] = 0;
      initSendEl( 2, l_ids, l_cIds );
      initRecvFa( 2, l_ids, l_cIds );
    }
  }
  // top
  for( int l_yh = 0; l_yh < (int) m_nHex[1]; l_yh++ ) {
    for( int l_xh = 0; l_xh < (int) m_nHex[0]; l_xh++ ) {
      l_ids[0] = l_xh; l_ids[1] = l_yh; l_ids[2] = m_nHex[2]-1;
      initSendEl( 2, l_ids, l_cIds );
      initSendFa( 2, l_ids, l_cIds );
    }
  }

  /*
   * iterate over tets and set up receive elements
   *
   * If two MPI-ranks share more than one domain-face, we have to reorder the messages to match:
   *
   *  ****************************************
   *  *  9  10  11  12  13   14  15  16  17  *
   *  *                                      *
   *  *       Default for rank 1's           *
   *  *       receive elements               *
   *  *                                      *
   *  *  0   1   2   3   4   5   6   7   8   *
   *  ****************************************
   *  *  9  10  11  12  13   14  15  16  17  *
   *  *                                      *
   *  *       Default for rank 0's           *
   *  *       send elements                  *
   *  *                                      *
   *  *  0   1   2   3   4   5   6   7   8   *
   *  ****************************************
   *
   * --------->
   *
   *  ****************************************
   *  *  0   1   2   3   4   5   6   7   8   *
   *  *                                      *
   *  *       Fixed receive elements         *
   *  *       for rank 1                     *
   *  *                                      *
   *  *  9  10  11  12  13   14  15  16  17  *
   *  ****************************************
   *  *  9  10  11  12  13   14  15  16  17  *
   *  *                                      *
   *  *       Default for rank 0's           *
   *  *       send elements                  *
   *  *                                      *
   *  *  0   1   2   3   4   5   6   7   8   *
   *  ****************************************
   *
   * This is why we change the order of the layers w.r.t. to send ents
   * for two boundaries going to the same rank.
   */
  // set up x-ids
  int l_idsRe[2][3];
  if( m_mpiNe[0][0] != m_mpiNe[0][1] ) {
    l_idsRe[0][0] = -1;
    l_idsRe[1][0] = m_nHex[0];
  }
  else {
    l_idsRe[0][0] = m_nHex[0];
    l_idsRe[1][0] = -1;
  }

  // left-right
  for( unsigned short l_si = 0; l_si < 2; l_si++ ) {
    for( int l_zh = 0; l_zh < (int) m_nHex[2]; l_zh++ ) {
      for( int l_yh = 0; l_yh < (int) m_nHex[1]; l_yh++ ) {
        l_idsRe[l_si][1] = l_yh; l_idsRe[l_si][2] = l_zh;
        initRecvEl( 0, l_idsRe[l_si], l_cIds );
      }
    }
  }

  // set up y-ids
  if( m_mpiNe[1][0] != m_mpiNe[1][1] ) {
    l_idsRe[0][1] = -1;
    l_idsRe[1][1] = m_nHex[1];
  }
  else {
    l_idsRe[0][1] = m_nHex[1];
    l_idsRe[1][1] = -1;
  }

  // front-back
  for( unsigned short l_si = 0; l_si < 2; l_si++ ) {
    for( int l_zh = 0; l_zh < (int) m_nHex[2]; l_zh++ ) {
      for( int l_xh = 0; l_xh < (int) m_nHex[0]; l_xh++ ) {
        l_idsRe[l_si][0] = l_xh; l_idsRe[l_si][2] = l_zh;
        initRecvEl( 1, l_idsRe[l_si], l_cIds );
      }
    }
  }

  // set up z-dis
  if( m_mpiNe[2][0] != m_mpiNe[2][1] ) {
    l_idsRe[0][2] = -1;
    l_idsRe[1][2] = m_nHex[2];
  }
  else {
    l_idsRe[0][2] = m_nHex[2];
    l_idsRe[1][2] = -1;
  }

  // bottom-top
  for( unsigned short l_si = 0; l_si < 2; l_si++ ) {
    for( int l_yh = 0; l_yh < (int) m_nHex[1]; l_yh++ ) {
      for( int l_xh = 0; l_xh < (int) m_nHex[0]; l_xh++ ) {
        l_idsRe[l_si][0] = l_xh; l_idsRe[l_si][1] = l_yh;
        initRecvEl( 2, l_idsRe[l_si], l_cIds );
      }
    }
  }

  // set up face info for ghost elements
  if( m_periodic ) {
    for( int l_zh = 0; l_zh < (int) m_nHex[2]; l_zh++ ) {
      for( int l_yh = 0; l_yh < (int) m_nHex[1]; l_yh++ ) {
        // left
        if( m_mpiNe[0][0] != parallel::g_rank ) {
          unsigned int l_hx   = getHxIdBnd(  0, l_yh, l_zh, m_nHex[0], m_nHex[1], m_nHex[2] );
          unsigned int l_hxNe = getHxIdBnd( -1, l_yh, l_zh, m_nHex[0], m_nHex[1], m_nHex[2] );

          for( unsigned short l_fa = 2; l_fa < 4; l_fa++ )
            m_hexes[l_hxNe].fa[l_fa] = m_hexes[l_hx].fa[l_fa-2];
        }
        // right
        if( m_mpiNe[0][1] != parallel::g_rank ) {
          unsigned int l_hx   = getHxIdBnd(  m_nHex[0]-1, l_yh, l_zh, m_nHex[0], m_nHex[1], m_nHex[2] );
          unsigned int l_hxNe = getHxIdBnd(  m_nHex[0]  , l_yh, l_zh, m_nHex[0], m_nHex[1], m_nHex[2] );
          for( unsigned short l_fa = 0; l_fa < 2; l_fa++ )
            m_hexes[l_hxNe].fa[l_fa] = m_hexes[l_hx].fa[l_fa+2];
        }
      }
    }

  for( int l_zh = 0; l_zh < (int) m_nHex[2]; l_zh++ ) {
    for( int l_xh = 0; l_xh < (int) m_nHex[0]; l_xh++ ) {
        // front
        if( m_mpiNe[1][0] != parallel::g_rank ) {
          unsigned int l_hx   = getHxIdBnd(  l_xh,  0, l_zh, m_nHex[0], m_nHex[1], m_nHex[2] );
          unsigned int l_hxNe = getHxIdBnd(  l_xh, -1, l_zh, m_nHex[0], m_nHex[1], m_nHex[2] );

          for( unsigned short l_fa = 6; l_fa < 8; l_fa++ )
            m_hexes[l_hxNe].fa[l_fa] = m_hexes[l_hx].fa[l_fa-2];
        }
        // back
        if( m_mpiNe[1][1] != parallel::g_rank ) {
          unsigned int l_hx   = getHxIdBnd(  l_xh, m_nHex[1]-1, l_zh, m_nHex[0], m_nHex[1], m_nHex[2] );
          unsigned int l_hxNe = getHxIdBnd(  l_xh, m_nHex[1],   l_zh, m_nHex[0], m_nHex[1], m_nHex[2] );
          for( unsigned short l_fa = 4; l_fa < 6; l_fa++ )
            m_hexes[l_hxNe].fa[l_fa] = m_hexes[l_hx].fa[l_fa+2];
        }
      }
    }

    for( int l_yh = 0; l_yh < (int) m_nHex[1]; l_yh++ ) {
      for( int l_xh = 0; l_xh < (int) m_nHex[0]; l_xh++ ) {
        // bottom
        if( m_mpiNe[2][0] != parallel::g_rank ) {
          unsigned int l_hx   = getHxIdBnd( l_xh, l_yh,  0, m_nHex[0], m_nHex[1], m_nHex[2] );
          unsigned int l_hxNe = getHxIdBnd( l_xh, l_yh, -1, m_nHex[0], m_nHex[1], m_nHex[2] );

          for( unsigned short l_fa = 10; l_fa < 12; l_fa++ )
            m_hexes[l_hxNe].fa[l_fa] = m_hexes[l_hx].fa[l_fa-2];
        }
        // top
        if( m_mpiNe[2][1] != parallel::g_rank ) {
          unsigned int l_hx   = getHxIdBnd( l_xh, l_yh,  m_nHex[2]-1, m_nHex[0], m_nHex[1], m_nHex[2] );
          unsigned int l_hxNe = getHxIdBnd( l_xh, l_yh,  m_nHex[2],   m_nHex[0], m_nHex[1], m_nHex[2] );
          for( unsigned short l_fa = 8; l_fa < 10; l_fa++ )
            m_hexes[l_hxNe].fa[l_fa] = m_hexes[l_hx].fa[l_fa+2];
        }
      }
    }
  }

  // derive index mappings
  deriveInMapVeAndFa();

  // get global ids
  getGIdsVe( m_gIdsVe ); getGIdsFa( m_gIdsFa ); getGIdsEl( m_gIdsEl );

  // intialize 'first'-positions
  initFirstPos( m_veLayout ); initFirstPos( m_faLayout ); initFirstPos( m_elLayout );
}

void edge::mesh::regular::Tet::init( const unsigned int i_nHex[3],
                                                    int i_rank,
                                                    int i_nRanks,
                                     const double       i_corner[3],
                                     const double       i_dX[3],
                                           bool         i_periodic ) {
  EDGE_CHECK( i_periodic == true ); // TODO: Add support for non-periodic boundaries
  m_periodic = i_periodic;

  // periodic setups must have alternating types
  if( m_periodic ) {
    EDGE_CHECK( i_nHex[0]%2 == 0 );
    EDGE_CHECK( i_nHex[1]%2 == 0 );
    EDGE_CHECK( i_nHex[2]%2 == 0 );
  }

  // get a base disc
  unsigned int l_nHexLoc[3], l_nPart[3], l_part[3];

  Base::getSetup3d( i_nRanks,
                    i_rank,
                    m_periodic,
                    i_nHex,
                    l_nHexLoc,
                    l_nPart,
                    l_part );

  EDGE_LOG_INFO << "  our per-dimension mpi-setup (x,y,z): "
                << l_nPart[0] << " " << l_nPart[1] << " " << l_nPart[2];

  // iterate over previous ranks and compute the element offset
  unsigned int l_offset[3] = {0,0,0};

  for( unsigned short l_di = 0; l_di < 3; l_di++ ) {
    for( unsigned int l_pr = 0; l_pr < l_part[l_di]; l_pr++ ) {
      unsigned int l_partPr[3];
      l_partPr[0] = l_part[0]; l_partPr[1] = l_part[1]; l_partPr[2] = l_part[2];
      l_partPr[l_di] = l_pr;

      int l_prRank = Base::getRank( l_nPart, l_partPr );

      unsigned int l_nHexPr[3], l_nPartPr[3];
      Base::getSetup3d( i_nRanks,
                        l_prRank,
                        m_periodic, // TODO: fix periodic
                        i_nHex,
                        l_nHexPr,
                        l_nPartPr,
                        l_partPr );
      l_offset[l_di] += l_nHexPr[l_di];
    }
  }

  // derive local corner of the rank's domain
  double l_coLoc[3];
  for( unsigned short l_di = 0; l_di < 3; l_di++ ) {
    l_coLoc[l_di] = i_corner[l_di] + l_offset[l_di] * i_dX[l_di];
  }

  // setup mpi neighbors
  int l_mpiNe[3][2];

  for( unsigned short l_di = 0; l_di < 3; l_di++ ) {
    unsigned int l_partNe[3];

    // 'left' neighbor
    if( l_part[l_di] > 0 ) {
      l_partNe[0] = l_part[0]; l_partNe[1] = l_part[1]; l_partNe[2] = l_part[2];
      l_partNe[l_di] -= 1;

      l_mpiNe[l_di][0] = Base::getRank( l_nPart, l_partNe );
    }
    else {
      if( !m_periodic ) l_mpiNe[l_di][0] = -1;
      else {
        l_partNe[0] = l_part[0]; l_partNe[1] = l_part[1]; l_partNe[2] = l_part[2];
        l_partNe[l_di] = l_nPart[l_di]-1;

        l_mpiNe[l_di][0] = Base::getRank( l_nPart, l_partNe );
      }
    }

    // 'right' neighbor
    if( l_part[l_di] < l_nPart[l_di]-1 ) {
      l_partNe[0] = l_part[0]; l_partNe[1] = l_part[1]; l_partNe[2] = l_part[2];
      l_partNe[l_di] += 1;

      l_mpiNe[l_di][1] = Base::getRank( l_nPart, l_partNe );
    }
    else {
      if( !m_periodic ) l_mpiNe[l_di][1] = -1;
      else {
        l_partNe[0] = l_part[0]; l_partNe[1] = l_part[1]; l_partNe[2] = l_part[2];
        l_partNe[l_di] = 0;

        l_mpiNe[l_di][1] = Base::getRank( l_nPart, l_partNe );
      }
    }
  }

  // derive type
  bool l_type = (l_offset[0] + l_offset[1] + l_offset[2])%2;

  // do the initialization
  init( l_type,
        l_nHexLoc,
        l_mpiNe,
        l_coLoc,
        i_dX );
}

void edge::mesh::regular::Tet::getElVe( int_el (*o_elVe)[4] ) {
  // iterate over hexes
  for( std::size_t l_hx = 0; l_hx < m_hexes.size(); l_hx++ ) {
    // iterate over tets
    for( unsigned short l_te = 0; l_te < 5; l_te++ ) {
      // iterate over ids
      for( unsigned short l_el = 0; l_el < m_hexes[l_hx].el[l_te].size(); l_el++ ) {
        for( unsigned short l_ve = 0; l_ve < 4; l_ve++ ) {
          unsigned short l_locId = m_mapElVe[m_hexes[l_hx].type][l_te][l_ve];

          o_elVe[ m_hexes[l_hx].el[l_te][l_el] ][l_ve] = m_hexes[l_hx].ve[l_locId];
        }
      }
    }
  }
}

void edge::mesh::regular::Tet::getFaVe( int_el (*o_faVe)[3] ) {
  // iterate over owned hexes having faces
  for( int l_zh = 0; l_zh < (int) m_nHex[2]; l_zh++ ) {
    for( int l_yh = 0; l_yh < (int) m_nHex[1]; l_yh++ ) {
      for( int l_xh = 0; l_xh < (int) m_nHex[0]; l_xh++ ) {
        unsigned int l_hx = getHxId( l_xh, l_yh, l_zh, m_nHex[0], m_nHex[1], m_nHex[2] );

        // iterate over hex faces
        for( unsigned short l_fa = 0; l_fa < 16; l_fa++ ) {
          // iterate over face vertices
          for( unsigned short l_ve = 0; l_ve < 3; l_ve++ ) {
            unsigned short l_veId = m_mapFaVe[m_hexes[l_hx].type][l_fa][l_ve];

            o_faVe[ m_hexes[l_hx].fa[l_fa] ][l_ve] = m_hexes[l_hx].ve[l_veId];
          }
        }
      }
    }
  }
}

void edge::mesh::regular::Tet::getFaEl( int_el (*o_faEl)[2] ) {
  // maps faces to tets (with sorted global ids of the element ids)
  std::function< void( unsigned short, const int i_hxPos[3], int_el[2] )  > l_map =
    [&](       unsigned short i_fa,
         const int            i_hxPos[3],
               int_el         o_faEls[2] ) {
      // get id of owning hex
      unsigned int l_hx = getHxId( i_hxPos[0], i_hxPos[1], i_hxPos[2], m_nHex[0], m_nHex[1], m_nHex[2] );

      // get offsets
      unsigned int l_hxOff[2][3];
      for( unsigned short l_di = 0; l_di < 3; l_di++ ) {
        l_hxOff[0][l_di] = m_mapFaHxOff[i_fa][0][l_di];
        l_hxOff[1][l_di] = m_mapFaHxOff[i_fa][1][l_di];
      }

      // get indices of neighboring hexes
      unsigned int l_hxNe[2];
      for( unsigned short l_si = 0; l_si < 2; l_si++ ) {
        l_hxNe[l_si] = getHxIdBnd( i_hxPos[0]+l_hxOff[l_si][0], i_hxPos[1]+l_hxOff[l_si][1], i_hxPos[2]+l_hxOff[l_si][2],
                                   m_nHex[0],                   m_nHex[1],                   m_nHex[2] );
      }

      // at least one hex must be valid
      EDGE_CHECK( l_hxNe[0] != std::numeric_limits< unsigned int >::max() ||
                  l_hxNe[1] != std::numeric_limits< unsigned int >::max() );

      // get the local ids of the adjacent elements
      unsigned short l_el[2];
      l_el[0] = m_mapFaEl[ m_hexes[l_hx].type ][i_fa][0];
      l_el[1] = m_mapFaEl[ m_hexes[l_hx].type ][i_fa][1];

      // assign the elements
      if( l_hxNe[0] == std::numeric_limits<unsigned int>::max() ) { // boundary condition for hx0

        o_faEls[0] = m_hexes[l_hxNe[1]].el[ l_el[1] ][0];
        o_faEls[1] = std::numeric_limits<int_el>::max();
      }
      else if( l_hxNe[1] == std::numeric_limits<unsigned int>::max() ) { // boundary condition for hx1
        EDGE_CHECK( m_hexes[l_hxNe[0]].el[ l_el[1] ].size() > 0 );
        o_faEls[0] = m_hexes[l_hxNe[0]].el[ l_el[0] ][0];
        o_faEls[1] = std::numeric_limits<int_el>::max();
      }
      else { // non-bnd
        EDGE_CHECK( m_hexes[l_hxNe[0]].el[ l_el[0] ].size() > 0 );
        EDGE_CHECK( m_hexes[l_hxNe[1]].el[ l_el[1] ].size() > 0 );
        EDGE_CHECK( m_hexes[l_hxNe[0]].el[ l_el[0] ][0] != m_hexes[l_hxNe[1]].el[ l_el[1] ][0] );

        int_el l_gIds[2];
        l_gIds[0] = m_gIdsEl[ m_hexes[l_hxNe[0]].el[ l_el[0] ][0] ];
        l_gIds[1] = m_gIdsEl[ m_hexes[l_hxNe[1]].el[ l_el[1] ][0] ];

        // same element id can not happen
        EDGE_CHECK( l_gIds[0] != l_gIds[1] );

        if( l_gIds[0] < l_gIds[1] ) {
          o_faEls[0] = m_hexes[l_hxNe[0]].el[ l_el[0] ][0];
          o_faEls[1] = m_hexes[l_hxNe[1]].el[ l_el[1] ][0];
        }
        else {
          o_faEls[0] = m_hexes[l_hxNe[1]].el[ l_el[1] ][0];
          o_faEls[1] = m_hexes[l_hxNe[0]].el[ l_el[0] ][0];
        }
      }
    };

  // iterate over owned hexes having faces
  int l_hxPos[3] = {0,0,0};
  for( l_hxPos[2] = 0; l_hxPos[2] < (int) m_nHex[2]; l_hxPos[2]++ ) {
    for( l_hxPos[1] = 0; l_hxPos[1] < (int) m_nHex[1]; l_hxPos[1]++ ) {
      for( l_hxPos[0] = 0; l_hxPos[0] < (int) m_nHex[0]; l_hxPos[0]++ ) {

        for( unsigned short l_fa = 0; l_fa < 16; l_fa++ ) {
          unsigned int l_hx = getHxId( l_hxPos[0], l_hxPos[1], l_hxPos[2], m_nHex[0], m_nHex[1], m_nHex[2] );

          l_map( l_fa,
                 l_hxPos,
                 o_faEl[ m_hexes[l_hx].fa[l_fa] ] );
        }
      }
    }
  }
}

void edge::mesh::regular::Tet::getElFa( int_el (*o_elFa)[4] ) {
  // maps tets to faces (with sorted global ids of the element ids)
  std::function< void( unsigned short, int i_xh, int i_yh, int i_zh, int_el[4] )  > l_map =
    [&]( unsigned short i_te,
         int            i_xh,
         int            i_yh,
         int            i_zh,
         int_el         o_elFas[4] ) {
      // get id of owning hex
      unsigned int l_hx = getHxId( i_xh, i_yh, i_zh, m_nHex[0], m_nHex[1], m_nHex[2] );

      // iterate over faces
      for( unsigned short l_fa = 0; l_fa < 4; l_fa++ ) {
        bool l_type = m_hexes[l_hx].type;
        unsigned short l_faId = m_mapElFa[l_type][i_te][l_fa];
        o_elFas[l_fa] = m_hexes[l_hx].fa[l_faId];
      }
    };

  // iterate over all hexes
  for( int l_zh = -1; l_zh < (int) m_nHex[2]+1; l_zh++ ) {
    for( int l_yh = -1; l_yh < (int) m_nHex[1]+1; l_yh++ ) {
      for( int l_xh = -1; l_xh < (int) m_nHex[0]+1; l_xh++ ) {
        unsigned int l_hx = getHxId( l_xh, l_yh, l_zh, m_nHex[0], m_nHex[1], m_nHex[2] );

        // iterate over possible tets
        for( unsigned short l_te = 0; l_te < 5; l_te++ ) {
          // iterate over existing tets
          for( unsigned short l_ex = 0; l_ex < m_hexes[l_hx].el[l_te].size(); l_ex++ ) {
            l_map( l_te,
                   l_xh, l_yh, l_zh,
                   o_elFa[ m_hexes[l_hx].el[l_te][l_ex] ] );
          }
        }
      }
    }
  }
}

void edge::mesh::regular::Tet::getElFaEl( int_el (*o_elFaEl)[4] ) {
  // maps tets to face-adjacent tets
  std::function< void( unsigned short, int i_xh, int i_yh, int i_zh, int_el[4] )  > l_map =
    [&]( unsigned short i_te,
         int            i_xh,
         int            i_yh,
         int            i_zh,
         int_el         o_elFaElLoc[4] ) {
      // get id of owning hex
      unsigned int l_hx = getHxId( i_xh, i_yh, i_zh, m_nHex[0], m_nHex[1], m_nHex[2] );

      // iterate over faces
      for( unsigned short l_fa = 0; l_fa < 4; l_fa++ ) {
        bool l_ty = m_hexes[l_hx].type;

        // get the mapping info
        short l_mapElFaEl[4];
        for( unsigned short l_i = 0; l_i < 4; l_i++ ) l_mapElFaEl[l_i] = m_mapElFaEl[l_ty][i_te][l_fa][l_i];

        // derive the hex id containing the adjacent tet
        unsigned int l_hxN = getHxIdBnd( i_xh+l_mapElFaEl[1],
                                         i_yh+l_mapElFaEl[2],
                                         i_zh+l_mapElFaEl[3],
                                         m_nHex[0], m_nHex[1], m_nHex[2] );

        // check alternating hex types or same element
        EDGE_CHECK(  l_hxN == std::numeric_limits< unsigned int >::max() ||
                     m_hexes[l_hx].type != m_hexes[l_hxN].type ||
                     l_hx               == l_hxN );

        // check that the tet exists for owned elements
        EDGE_CHECK( l_hxN == std::numeric_limits< unsigned int >::max() ||
                    m_hexes[l_hxN].el[ l_mapElFaEl[0] ].size() > 0 ||
                    i_xh==-1 || i_yh==-1 || i_zh==-1 ||
                    i_xh==(int)m_nHex[0] || i_yh==(int)m_nHex[1] || i_zh==(int)m_nHex[2] );

        // assigned the info
        if( l_hxN != std::numeric_limits< unsigned int >::max() &&
            m_hexes[l_hxN].el[ l_mapElFaEl[0] ].size() > 0 ) {
          o_elFaElLoc[l_fa] = m_hexes[l_hxN].el[ l_mapElFaEl[0] ][0];
        }
        else o_elFaElLoc[l_fa] = std::numeric_limits< int_el >::max();
      }
    };

  // iterate over all hexes
  for( int l_zh = -1; l_zh < (int) m_nHex[2]+1; l_zh++ ) {
    for( int l_yh = -1; l_yh < (int) m_nHex[1]+1; l_yh++ ) {
      for( int l_xh = -1; l_xh < (int) m_nHex[0]+1; l_xh++ ) {
        unsigned int l_hx = getHxId( l_xh, l_yh, l_zh, m_nHex[0], m_nHex[1], m_nHex[2] );

        // iterate over possible tets
        for( unsigned short l_te = 0; l_te < 5; l_te++ ) {
          // iterate over existing tets
          for( unsigned short l_ex = 0; l_ex < m_hexes[l_hx].el[l_te].size(); l_ex++ ) {
            l_map( l_te,
                   l_xh, l_yh, l_zh,
                   o_elFaEl[ m_hexes[l_hx].el[l_te][l_ex] ] );
          }
        }
      }
    }
  }
}

void edge::mesh::regular::Tet::getVeChars( t_vertexChars *o_veChars ) {
  // iterate over hexes
  for( int l_zh = -1; l_zh < (int) m_nHex[2]+1; l_zh++ ) {
    for( int l_yh = -1; l_yh < (int) m_nHex[1]+1; l_yh++ ) {
      for( int l_xh = -1; l_xh < (int) m_nHex[0]+1; l_xh++ ) {
        unsigned int l_hx = getHxId( l_xh, l_yh, l_zh, m_nHex[0], m_nHex[1], m_nHex[2] );

        for( unsigned short l_ve = 0; l_ve < 8; l_ve++ ) {
          double l_crds[3];
          l_crds[0] = m_corner[0];
          l_crds[1] = m_corner[1];
          l_crds[2] = m_corner[2];

          l_crds[0] += l_xh * m_dX[0];
          l_crds[1] += l_yh * m_dX[1];
          l_crds[2] += l_zh * m_dX[2];

          if( l_ve == 1 || l_ve == 2 || l_ve == 5 || l_ve == 6 ) l_crds[0] += m_dX[0];
          if( l_ve == 2 || l_ve == 3 || l_ve == 6 || l_ve == 7 ) l_crds[1] += m_dX[1];
          if( l_ve == 4 || l_ve == 5 || l_ve == 6 || l_ve == 7 ) l_crds[2] += m_dX[2];

          o_veChars[ m_hexes[l_hx].ve[l_ve] ].coords[0] = l_crds[0];
          o_veChars[ m_hexes[l_hx].ve[l_ve] ].coords[1] = l_crds[1];
          o_veChars[ m_hexes[l_hx].ve[l_ve] ].coords[2] = l_crds[2];

          // set dummy vertex type
          o_veChars[ m_hexes[l_hx].ve[l_ve] ].spType = MESH_TYPE_NONE;
        }
      }
    }
  }
}

void edge::mesh::regular::Tet::getFaChars( t_faceChars *o_faChars ) {
  // get the vertices adjacent to the faces
  int_el l_nFa = getNFa();
  int_el (*l_faVe)[3] = ( int_el(*)[3] ) new int_el[3 * l_nFa];
  getFaVe( l_faVe );

  // get vertices coords
  int_el l_nVe = m_veLayout.nEnts;
  t_vertexChars *l_veChars = new t_vertexChars[l_nVe];
  getVeChars( l_veChars );

  // iterate over faces
  for( int_el l_fa = 0; l_fa < l_nFa; l_fa++ ) {
    // get the vertices
    real_mesh l_ves[3][3];
    for( unsigned short l_ve = 0; l_ve < 3; l_ve++ ) {
      for( unsigned short l_di = 0; l_di < 3; l_di++ ) {
        l_ves[l_di][l_ve] = l_veChars[ l_faVe[l_fa][l_ve] ].coords[l_di];
      }
    }

    // assign faces' area
    o_faChars[l_fa].area = edge::linalg::Geom::volume( TRIA3, l_ves[0], 3 );

    // set dummy face type
    o_faChars[l_fa].spType = MESH_TYPE_NONE;

    // assign tangents
    linalg::Geom::computeTangents( TRIA3,
                                   l_ves[0],
                                   o_faChars[l_fa].tangent0,
                                   o_faChars[l_fa].tangent1 );
  }

  // get the elements
  int_el l_nEl = m_elLayout.nEnts;
  int_el (*l_elFa)[4] = ( int_el(*)[4] ) new int_el[4 * l_nEl];
  getElFa( l_elFa );

  int_el (*l_elFaEl)[4] = ( int_el(*)[4] ) new int_el[4 * l_nEl];
  getElFaEl( l_elFaEl );

  int_el (*l_elVe)[4] = ( int_el(*)[4] ) new int_el[4 * l_nEl];
  getElVe( l_elVe );

  // iterate over tets
  for( int_el l_el = 0; l_el < l_nEl ; l_el++ ) {
    // iterate over faces
    for( unsigned short l_fa = 0; l_fa < 4; l_fa++ ) {
      int_el l_faId = l_elFa[l_el][l_fa];

      // ignore ghost faces
      if( l_faId == std::numeric_limits<int_el>::max() ) continue;

      // get face coords
      real_mesh l_ves[3][3];
      for( unsigned short l_ve = 0; l_ve < 3; l_ve++ ) {
        for( unsigned short l_di = 0; l_di < 3; l_di++ ) {
          int_el l_veId = l_faVe[l_faId][l_ve];
          l_ves[l_di][l_ve] = l_veChars[l_veId].coords[l_di];
        }
      }

      // get normal point of tet
      real_mesh l_np[3];
      int_el l_npId = l_elVe[l_el][3-l_fa];
      l_np[0] = l_veChars[l_npId].coords[0];
      l_np[1] = l_veChars[l_npId].coords[1];
      l_np[2] = l_veChars[l_npId].coords[2];

      linalg::Geom::computeOutPtNormal( TRIA3,
                                        l_ves[0],
                                        l_np,
                                        o_faChars[l_faId].outNormal );

      // invert if out normal points inside the tet
      if( l_elFaEl[l_el][l_fa] != std::numeric_limits< int_el >::max() &&
          m_gIdsEl[l_elFaEl[l_el][l_fa]] < m_gIdsEl[l_el] ) {
        for( unsigned short l_di = 0; l_di < 3; l_di++ ) o_faChars[l_faId].outNormal[l_di] *= -1;
      }
    }
  }

  delete[] l_faVe;
  delete[] l_veChars;
  delete[] l_elFa;
  delete[] l_elVe;
  delete[] l_elFaEl;
}

void edge::mesh::regular::Tet::getElChars( t_elementChars *o_elChars ) {
  // get vertices coords
  int_el l_nVe = m_veLayout.nEnts;
  t_vertexChars *l_veChars = new t_vertexChars[l_nVe];
  getVeChars( l_veChars );

  // get the elements
  int_el l_nEl = m_elLayout.nEnts;
  int_el (*l_elVe)[4] = ( int_el(*)[4] ) new int_el[4 * l_nEl];
  getElVe( l_elVe );

  // iterate over tets
  for( int_el l_el = 0; l_el < l_nEl ; l_el++ ) {
    // get vetex coords
    real_mesh l_ves[3][4];
    for( unsigned short l_ve = 0; l_ve < 4; l_ve++ ) {
      for( unsigned short l_di = 0; l_di < 3; l_di++ ) {
        int_el l_veId = l_elVe[l_el][l_ve];
        l_ves[l_di][l_ve] = l_veChars[l_veId].coords[l_di];
      }
    }

    // compute volume
    o_elChars[l_el].volume = edge::linalg::Geom::volume( TET4, l_ves[0] );
    // compute in-sphere diameter
    o_elChars[l_el].inDia = edge::linalg::Geom::inDia(  TET4, l_ves[0] );

    // set dummy element type
    o_elChars[l_el].spType = MESH_TYPE_NONE;
  }

  delete[] l_elVe;
  delete[] l_veChars;
}

void edge::mesh::regular::Tet::getGIdsVe( std::vector< int_gid > &o_gIds ) {
  o_gIds.resize( m_veLayout.nEnts );
  // set info
  for( int_el l_gi = 0; l_gi < m_veLayout.nEnts; l_gi++ ) {
    o_gIds[l_gi] = l_gi;
  }
}

void edge::mesh::regular::Tet::getGIdsFa( std::vector< int_gid > &o_gIds ) {
  o_gIds.resize( m_faLayout.nEnts );
  // set info
  for( int_el l_gi = 0; l_gi < m_faLayout.nEnts; l_gi++ ) {
    o_gIds[l_gi] = l_gi;
  }
}

void edge::mesh::regular::Tet::getGIdsEl( std::vector< int_gid > &o_gIds ) {
  o_gIds.resize( m_elLayout.nEnts );
  for( std::size_t l_id = 0; l_id < o_gIds.size(); l_id++ ) o_gIds[l_id] = std::numeric_limits<int_gid>::max();

  int_gid l_gTmp = 0;

  // store global ids dimension wise
  // TODO: For a real global id, we would have to consider the offsets by the neighbors
  for( int l_zh = -1; l_zh < (int) m_nHex[2]+1; l_zh++ ) {
    for( int l_yh = -1; l_yh < (int) m_nHex[1]+1; l_yh++ ) {
      for( int l_xh = -1; l_xh < (int) m_nHex[0]+1; l_xh++ ) {
        unsigned int l_hx = getHxId( l_xh, l_yh, l_zh, m_nHex[0], m_nHex[1], m_nHex[2] );
        EDGE_CHECK( l_hx < m_hexes.size() );

        for( unsigned short l_te = 0; l_te < 5; l_te++ ) {
          for( unsigned short l_el = 0; l_el < m_hexes[l_hx].el[l_te].size(); l_el++ ) {
            int_el l_elId = m_hexes[l_hx].el[l_te][l_el];
            o_gIds[l_elId] = l_gTmp;
          }
          l_gTmp++;
        }
      }
    }
  }

  // check that every element has an assigned global id
  for( std::size_t l_id = 0; l_id < o_gIds.size(); l_id++ )
    EDGE_CHECK( o_gIds[l_id] != std::numeric_limits<int_gid>::max() );
}

const t_inMap* edge::mesh::regular::Tet::getInMap() const {
  return &m_inMap;
}

//////////////////////////////////////////////////////////////////////////////////////
//// TODO: This is very close to the Unstructured impl, reconsider implementation ////
//////////////////////////////////////////////////////////////////////////////////////
void edge::mesh::regular::Tet::getConnect( const t_vertexChars  *i_veChars,
                                           const t_faceChars    *i_faChars,
                                                 t_connect      &o_connect ) {
  // get info
  getFaVe(   o_connect.faVe   );
  getElVe(   o_connect.elVe   );
  getFaEl(   o_connect.faEl   );
  getElFa(   o_connect.elFa   );
  getElFaEl( o_connect.elFaEl );

  // check consistency of our ordering
  common< TET4 >::checkConsOrd( m_faLayout,
                                m_elLayout,
                                m_gIdsVe,
                                m_gIdsFa,
                                m_gIdsEl,
                                o_connect.elVe,
                                o_connect.elFa,
                                o_connect.faEl,
                                o_connect.faVe,
                                m_periodic
                              );


  // normalize vertex ordering
  common< TET4 >::normOrdTet4( m_elLayout,
                               i_veChars,
                               o_connect.elVe,
                               o_connect.elFa,
                               o_connect.elFaEl );

  // get the index mappings
  t_inMap const * l_inMap = getInMap();

  // check consistency of the adjacency infos
  common< TET4 >::checkConsAdj( m_faLayout,
                                i_faChars,
                                l_inMap->faMeDa,
                                l_inMap->faDaMe,
                                o_connect.elFaEl,
                                o_connect.elFa,
                                o_connect.faEl );

  // get the neighboring elements' local face ids
  common< TET4 >::getFIdsElFaEl( m_elLayout,
                                 m_gIdsEl,
                                 o_connect.elFa,
                                 o_connect.elFaEl,
                                 o_connect.fIdElFaEl );

  // get the neighboring elements' local vertex ids, which match the faces' first vertex
  common< TET4 >::getVIdsElFaEl( m_elLayout,
                                 m_gIdsEl,
                                 o_connect.elFaEl,
                                 o_connect.elVe,
                                 o_connect.fIdElFaEl,
                                 o_connect.vIdElFaEl,
                                 true, // TODO: periodic enforced here
                                 i_veChars  );
}
