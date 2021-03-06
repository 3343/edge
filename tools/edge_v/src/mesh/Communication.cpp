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
 * Derives communication structures.
 **/
#include "Communication.h"
#include "../io/logging.h"

bool edge_v::mesh::Communication::isComm( t_entityType         i_elTy,
                                          t_idx                i_el,
                                          t_idx        const * i_elFaEl,
                                          t_idx        const * i_elPa ) {
  unsigned short l_nElFas = CE_N_FAS( i_elTy );
  bool l_isComm = false;

  // element's partition
  t_idx l_elPa = i_elPa[i_el];

  for( unsigned short l_fa = 0; l_fa < l_nElFas; l_fa++ ) {
    t_idx l_ad = i_elFaEl[ i_el*l_nElFas + l_fa ];
    if( l_ad != std::numeric_limits< t_idx >::max() ) {
      // neighboring element's partition
      t_idx l_adPa = i_elPa[l_ad];
      if( l_elPa != l_adPa ) l_isComm = true;
    }
  }
  return l_isComm;
};

void edge_v::mesh::Communication::getPaElComm( t_entityType         i_elTy,
                                               t_idx                i_nEls,
                                               t_idx        const * i_elFaEl,
                                               t_idx        const * i_elPa,
                                               t_idx              * o_first,
                                               t_idx              * o_size ) {
  // partition of the communication region
  t_idx l_regPa = std::numeric_limits< t_idx >::max();

  for( t_idx l_el = 0; l_el < i_nEls; l_el++ ) {
    bool l_commEl = isComm( i_elTy,
                            l_el,
                            i_elFaEl,
                            i_elPa );

    t_idx l_elPa = i_elPa[l_el];

    // init a new comm region if a new comm region shows up
    if( l_regPa != l_elPa && l_commEl == true  ) {
      // sanity check on the element order
      EDGE_V_CHECK( l_regPa == std::numeric_limits< t_idx >::max() ||
                    l_elPa  == l_regPa+1 );

      o_first[l_elPa] = l_el;
      l_regPa = l_elPa;
      o_size[l_elPa] = 1;
    }
    // increase size of the comm region for every element inside
    else if( (l_regPa == l_elPa) && l_commEl ) {
      o_size[l_elPa] += 1;
    }
  }
}

void edge_v::mesh::Communication::getPaTgPairs( t_entityType                          i_elTy,
                                                t_idx                                 i_first,
                                                t_idx                                 i_size,
                                                t_idx                         const * i_elFaEl,
                                                unsigned short                const * i_elTg,
                                                t_idx                         const * i_elPa,
                                                std::set<
                                                std::pair< t_idx,
                                                           unsigned short > >       & o_pairs ) {
  unsigned short l_nElFas = CE_N_FAS( i_elTy );
  EDGE_V_CHECK_GT( i_size, 0 );

  // reset pairs
  o_pairs.empty();

  // get the region's partition and time group
  unsigned short l_pa = i_elPa[i_first];
  unsigned short l_tg = i_elTg[i_first];

  // iterate over the elements of the region
  for( t_idx l_el = i_first; l_el < i_first+i_size; l_el++ ) {
    // check that the partition and time group match
    EDGE_V_CHECK_EQ( i_elPa[l_el], l_pa );
    EDGE_V_CHECK_EQ( i_elTg[l_el], l_tg );

    // iterate over the face-adjacent elements
    for( unsigned short l_fa = 0; l_fa < l_nElFas; l_fa++ ) {
      t_idx l_ad = i_elFaEl[ l_el*l_nElFas + l_fa ];
      if( l_ad != std::numeric_limits< t_idx >::max() ) {
        // get adjacent partition and time group
        t_idx l_adPa = i_elPa[l_ad];
        t_idx l_adTg = i_elTg[l_ad];

        // communication is required, if adjacent elements belongs to a different partition
        if( l_pa != l_adPa ) {
          std::pair< t_idx, unsigned short > l_pair( {l_adPa, l_adTg} );
          o_pairs.insert( l_pair );
        }
      }
    }
  }
}

void edge_v::mesh::Communication::getMsgsSend( t_entityType                   i_elTy,
                                               t_idx                          i_first,
                                               t_idx                          i_size,
                                               t_idx                  const * i_elFaEl,
                                               t_idx                  const * i_elPa,
                                               unsigned short         const * i_elTg,
                                               std::vector< Message >       & o_send ) {
  // get the communication partition-time group pairs
  std::set< std::pair< t_idx,
                       unsigned short > > l_pairs;

  getPaTgPairs( i_elTy,
                i_first,
                i_size,
                i_elFaEl,
                i_elTg,
                i_elPa,
                l_pairs );

  // init the output
  o_send.resize(0);
  o_send.resize( l_pairs.size() );

  // store partitions and time groups
  t_idx l_msg = 0;
  for( auto l_it = l_pairs.begin(); l_it != l_pairs.end(); l_it++ ) {
    o_send[l_msg].pa = std::get<0>(*l_it);
    o_send[l_msg].tg = std::get<1>(*l_it);
    l_msg++;
  }

  // get the communication region's partition and time region
  t_idx l_paReg = i_elPa[i_first];

  unsigned short l_nElFas = CE_N_FAS( i_elTy );

  // lambda, which inserts an element-face pair into the send messages
  auto l_insert = [   i_elFaEl,
                      l_nElFas,
                    & o_send    ]( t_idx          i_el,
                                   unsigned short i_fa,
                                   t_idx          i_pa,
                                   unsigned short i_tg ) mutable {
    for( t_idx l_se = 0; l_se < o_send.size(); l_se++ ) {
      if( o_send[l_se].pa == i_pa && o_send[l_se].tg == i_tg ) {
        o_send[l_se].el.push_back( i_el );
        o_send[l_se].fa.push_back( i_fa );

        t_idx l_elAd = i_elFaEl[i_el*l_nElFas + i_fa];
        unsigned short l_faAd = std::numeric_limits< unsigned short >::max();
        for( unsigned short l_fa = 0; l_fa < l_nElFas; l_fa++ ) {
          if( i_elFaEl[l_elAd*l_nElFas + l_fa] == i_el ) l_faAd = l_fa;
        }
        EDGE_V_CHECK_NE( l_faAd, std::numeric_limits< unsigned short >::max() );
        o_send[l_se].elAd.push_back( l_elAd );
        o_send[l_se].faAd.push_back( l_faAd );

        break;
      }
      EDGE_V_CHECK_NE( l_se, o_send.size()-1 );
    }
  };

  // assemble the messages
  for( t_idx l_el = i_first; l_el < i_first+i_size; l_el++ ) {
    for( unsigned short l_fa = 0; l_fa < l_nElFas; l_fa++ ) {
      t_idx l_ad = i_elFaEl[l_el*l_nElFas + l_fa];

      if( l_ad != std::numeric_limits< t_idx >::max() ) {
        t_idx l_paAd = i_elPa[l_ad];
        unsigned short l_tgAd = i_elTg[l_ad];

        // insert the element-face pair into a message if the adjacent partition is different from that of the region
        if( l_paAd != l_paReg ) {
          l_insert( l_el,
                    l_fa,
                    l_paAd,
                    l_tgAd );
        }
      }
    }
  }
}

void edge_v::mesh::Communication::getStruct( t_entityType                   i_elTy,
                                             t_idx                          i_nEls,
                                             t_idx                  const * i_elFaEl,
                                             t_idx                  const * i_elPa,
                                             unsigned short         const * i_elTg,
                                             std::vector< Partition >     & o_struct ) {
  t_idx l_nPas = i_elPa[i_nEls-1]+1;
  t_idx *l_firstComm = new t_idx[ l_nPas ];
  t_idx *l_sizeComm  = new t_idx[ l_nPas ];

  for( t_idx l_pa = 0; l_pa < l_nPas; l_pa++ ) {
    l_firstComm[l_pa] = 0;
    l_sizeComm[l_pa] = 0;
  }

  // get communicating elements
  getPaElComm( i_elTy,
               i_nEls,
               i_elFaEl,
               i_elPa,
               l_firstComm,
               l_sizeComm );

  // lambda which gets the time groups of a communication region
  auto l_getTgs = [ i_elTg ]( t_idx                           i_first,
                              t_idx                           i_size,
                              std::vector< unsigned short > & o_tgs,
                              std::vector< t_idx >          & o_firstTg,
                              std::vector< t_idx >          & o_sizeTg ) {
    o_tgs.resize(0);
    o_firstTg.resize(0);
    o_sizeTg.resize(0);
    if( i_size == 0 ) return;

    unsigned short l_tg = std::numeric_limits< unsigned short >::max();
    for( t_idx l_el = i_first; l_el < i_first+i_size; l_el++ ) {
      unsigned short l_elTg = i_elTg[l_el];
      if( l_elTg != l_tg ) {
        o_tgs.push_back( l_elTg );
        o_firstTg.push_back( l_el );
        o_sizeTg.push_back( 1 );
        l_tg = l_elTg;
      }
      else {
        o_sizeTg.back()++;
      }
    }
  };

  // init comm structures
  o_struct.resize( 0 );
  o_struct.resize( l_nPas );

  // get the send parts
  for( unsigned short l_pa = 0; l_pa < l_nPas; l_pa++ ) {
    // assemble the time groups
    std::vector< t_idx > l_firstTg;
    std::vector< t_idx > l_sizeTg;
    std::vector< unsigned short > l_tgs;

    l_getTgs( l_firstComm[l_pa],
              l_sizeComm[l_pa],
              l_tgs,
              l_firstTg,
              l_sizeTg );

    o_struct[l_pa].tr.resize( l_tgs.size() );

    for( t_idx l_tg = 0; l_tg < l_tgs.size(); l_tg++ ) {
      o_struct[l_pa].tr[l_tg].tg = l_tgs[l_tg];

      getMsgsSend( i_elTy,
                   l_firstTg[l_tg],
                   l_sizeTg[l_tg],
                   i_elFaEl,
                   i_elPa,
                   i_elTg,
                   o_struct[l_pa].tr[l_tg].send );
    }
  }

  // mirror remote send for recv
  for( t_idx l_pa = 0; l_pa < l_nPas; l_pa++ ) {
    for( t_idx l_t0 = 0; l_t0 < o_struct[l_pa].tr.size(); l_t0++ ) {
      unsigned short l_tg = o_struct[l_pa].tr[l_t0].tg;

      t_idx l_nSend = o_struct[l_pa].tr[l_t0].send.size();
      o_struct[l_pa].tr[l_t0].recv.resize( l_nSend );

      // iterate over the send messages
      for( t_idx l_s0 = 0; l_s0 < l_nSend; l_s0++ ) {
        unsigned short l_tgAd = o_struct[l_pa].tr[l_t0].send[l_s0].tg;
        t_idx l_paAd = o_struct[l_pa].tr[l_t0].send[l_s0].pa;

        // find the matching send of the remote partition which reflects our recv
        for( t_idx l_t1 = 0; l_t1 < o_struct[l_paAd].tr.size(); l_t1++ ) {
          if( o_struct[l_paAd].tr[l_t1].tg == l_tgAd ) {
            for( t_idx l_s1 = 0; l_s1 < o_struct[l_paAd].tr[l_t1].send.size(); l_s1++ ) {
              if(    o_struct[l_paAd].tr[l_t1].send[l_s1].tg == l_tg
                  && o_struct[l_paAd].tr[l_t1].send[l_s1].pa == l_pa ) {
                // check for matching sizes
                EDGE_V_CHECK_EQ( o_struct[l_pa  ].tr[l_t0].send[l_s0].el.size(),
                                 o_struct[l_paAd].tr[l_t1].send[l_s1].el.size() );
                EDGE_V_CHECK_EQ( o_struct[l_pa  ].tr[l_t0].send[l_s0].fa.size(),
                                 o_struct[l_paAd].tr[l_t1].send[l_s1].fa.size() );

                o_struct[l_pa].tr[l_t0].recv[l_s0] = o_struct[l_paAd].tr[l_t1].send[l_s1];
                // reverse order of local and remote
                std::vector< t_idx > l_tmpEl = o_struct[l_pa].tr[l_t0].recv[l_s0].el;
                o_struct[l_pa].tr[l_t0].recv[l_s0].el = o_struct[l_pa].tr[l_t0].recv[l_s0].elAd;
                o_struct[l_pa].tr[l_t0].recv[l_s0].elAd = l_tmpEl;

                std::vector< unsigned short > l_tmpFa = o_struct[l_pa].tr[l_t0].recv[l_s0].fa;
                o_struct[l_pa].tr[l_t0].recv[l_s0].fa = o_struct[l_pa].tr[l_t0].recv[l_s0].faAd;
                o_struct[l_pa].tr[l_t0].recv[l_s0].faAd = l_tmpFa;

                o_struct[l_pa].tr[l_t0].recv[l_s0].pa = l_paAd;
                o_struct[l_pa].tr[l_t0].recv[l_s0].tg = l_tgAd;
                break;
              }
              EDGE_V_CHECK_NE( l_s1+1, o_struct[l_paAd].tr[l_t1].send.size() );
            }
            break;
          }
        }
      }
    }
  }

  delete[] l_firstComm;
  delete[] l_sizeComm;
}

void edge_v::mesh::Communication::setChOff( std::vector< Partition > const & i_struct,
                                            t_idx                          * o_chOff ) {
  o_chOff[0] = 0;
  for( t_idx l_pa = 0; l_pa < i_struct.size(); l_pa++ ) {
    o_chOff[l_pa+1] = o_chOff[l_pa];

    for( unsigned short l_tg = 0; l_tg < i_struct[l_pa].tr.size(); l_tg++ ) {
      o_chOff[l_pa+1] += i_struct[l_pa].tr[l_tg].send.size();
    }
  }
}

void edge_v::mesh::Communication::setChs( std::vector< Partition > const & i_struct,
                                          t_idx                    const * i_chOff,
                                          t_idx                          * o_chs ) {
  for( t_idx l_pa = 0; l_pa < i_struct.size(); l_pa++ ) {
    t_idx l_off = i_chOff[l_pa]*4 + l_pa;
    o_chs[l_off] = i_chOff[l_pa+1] - i_chOff[l_pa];
    l_off++;

    for( unsigned short l_tg = 0; l_tg < i_struct[l_pa].tr.size(); l_tg++ ) {
      for( unsigned short l_se = 0; l_se < i_struct[l_pa].tr[l_tg].send.size(); l_se++ ) {
        o_chs[l_off+0] = i_struct[l_pa].tr[l_tg].tg;
        o_chs[l_off+1] = i_struct[l_pa].tr[l_tg].send[l_se].pa;
        o_chs[l_off+2] = i_struct[l_pa].tr[l_tg].send[l_se].tg;
        o_chs[l_off+3] = i_struct[l_pa].tr[l_tg].send[l_se].el.size();
        l_off += 4;
      }
    }
  }
}

void edge_v::mesh::Communication::nElsInSe( t_entityType           i_elTy,
                                            t_idx                  i_nPas,
                                            unsigned short         i_nTgs,
                                            t_idx                  i_nEls,
                                            t_idx          const * i_elFaEl,
                                            t_idx          const * i_elPa,
                                            unsigned short const * i_elTg,
                                            t_idx                * o_nElsIn,
                                            t_idx                * o_nElsSe ) {
  // init
  for( t_idx l_pa = 0; l_pa < i_nPas; l_pa++ ) {
    for( unsigned short l_tg = 0; l_tg < i_nTgs; l_tg++ ) {
      o_nElsIn[l_pa*i_nTgs + l_tg] = 0;
      o_nElsSe[l_pa*i_nTgs + l_tg] = 0;
    }
  }

  // assign counts
  for( t_idx l_el = 0; l_el < i_nEls; l_el++ ) {
    bool l_commEl = isComm( i_elTy,
                            l_el,
                            i_elFaEl,
                            i_elPa );
    t_idx l_pa = i_elPa[l_el];
    unsigned short l_tg = i_elTg[l_el];

    if( !l_commEl ) {
      o_nElsIn[l_pa*i_nTgs + l_tg]++;
    }
    else {
      o_nElsSe[l_pa*i_nTgs + l_tg]++;
    }
  }
}

void edge_v::mesh::Communication::setSeReElFaOff( std::vector< Partition > const & i_struct,
                                                  t_idx                          * o_off ) {
  o_off[0] = 0;

  for( t_idx l_pa = 0; l_pa < i_struct.size(); l_pa++ ) {
    o_off[l_pa+1] = o_off[l_pa];

    for( t_idx l_tg = 0; l_tg < i_struct[l_pa].tr.size(); l_tg++ ) {
      for( t_idx l_ms = 0; l_ms < i_struct[l_pa].tr[l_tg].send.size(); l_ms++ ) {
        o_off[l_pa+1] += i_struct[l_pa].tr[l_tg].send[l_ms].el.size();
      }
    }
  }
}

void edge_v::mesh::Communication::setSeReElFa( t_idx                    const * i_nPaEls,
                                               std::vector< Partition > const & i_struct,
                                               unsigned short                 * o_sendFa,
                                               t_idx                          * o_sendEl,
                                               unsigned short                 * o_recvFa,
                                               t_idx                          * o_recvEl ) {
  t_idx l_id = 0;
  t_idx l_elFirst = 0;

  for( t_idx l_pa = 0; l_pa < i_struct.size(); l_pa++ ) {
    for( t_idx l_tg = 0; l_tg < i_struct[l_pa].tr.size(); l_tg++ ) {
      for( t_idx l_ms = 0; l_ms < i_struct[l_pa].tr[l_tg].send.size(); l_ms++ ) {
        for( t_idx l_en = 0; l_en < i_struct[l_pa].tr[l_tg].send[l_ms].fa.size(); l_en++ ) {
          o_sendFa[l_id] = i_struct[l_pa].tr[l_tg].send[l_ms].fa[l_en];
          o_sendEl[l_id] = i_struct[l_pa].tr[l_tg].send[l_ms].el[l_en] - l_elFirst;
          o_recvFa[l_id] = i_struct[l_pa].tr[l_tg].recv[l_ms].fa[l_en];
          o_recvEl[l_id] = i_struct[l_pa].tr[l_tg].recv[l_ms].el[l_en] - l_elFirst;
          l_id++;
        }
      }
    }

    l_elFirst += i_nPaEls[l_pa];
  }
}

edge_v::mesh::Communication::Communication( unsigned short         i_nTgs,
                                            t_entityType           i_elTy,
                                            t_idx                  i_nEls,
                                            t_idx          const * i_elFaEl,
                                            t_idx                  i_nPas,
                                            t_idx          const * i_nPaEls,
                                            unsigned short const * i_elTg ) {
  EDGE_V_CHECK_GT( i_nPas, 0 );
  m_nPas = i_nPas;
  m_nTgs = i_nTgs;

  // init temporary data structure for the parts
  t_idx * l_elPa = new t_idx[ i_nEls ];
  t_idx l_paFirst = 0;

  for( t_idx l_pa = 0; l_pa < i_nPas; l_pa++ ) {
    t_idx l_paSize = i_nPaEls[l_pa];
    EDGE_V_CHECK_LE( l_paFirst+l_paSize, i_nEls );
    for( t_idx l_el = l_paFirst; l_el < l_paFirst+l_paSize; l_el++ ) {
      l_elPa[l_el] = l_pa;
    }

    l_paFirst += i_nPaEls[l_pa];
  }

  // derive communication structure
  getStruct( i_elTy,
             i_nEls,
             i_elFaEl,
             l_elPa,
             i_elTg,
             m_struct );

  // get channels offsets
  m_chOff = new t_idx[i_nPas+1];
  setChOff( m_struct,
            m_chOff );

  // assign simplified comm structure
  m_chs = new t_idx[ i_nPas + m_chOff[i_nPas]*4 ];
  setChs( m_struct,
          m_chOff,
          m_chs );

  m_nElsIn = new t_idx[ i_nPas*m_nTgs ];
  m_nElsSe = new t_idx[ i_nPas*m_nTgs ];

  nElsInSe( i_elTy,
            m_nPas,
            m_nTgs,
            i_nEls,
            i_elFaEl,
            l_elPa,
            i_elTg,
            m_nElsIn,
            m_nElsSe );

  delete[] l_elPa;

  m_sendRecvOff = new t_idx[ i_nPas+1 ];
  setSeReElFaOff( m_struct,
                  m_sendRecvOff );

  t_idx l_nSeRe = m_sendRecvOff[ i_nPas ];
  m_sendFa = new unsigned short[ l_nSeRe ];
  m_sendEl = new t_idx[ l_nSeRe ];
  m_recvFa = new unsigned short[ l_nSeRe ];
  m_recvEl = new t_idx[ l_nSeRe ];
  setSeReElFa( i_nPaEls,
               m_struct,
               m_sendFa,
               m_sendEl,
               m_recvFa,
               m_recvEl );
}

edge_v::mesh::Communication::~Communication() {
  delete[] m_chOff;
  delete[] m_chs;
  delete[] m_nElsIn;
  delete[] m_nElsSe;
  delete[] m_sendRecvOff;
  delete[] m_sendFa;
  delete[] m_sendEl;
  delete[] m_recvFa;
  delete[] m_recvEl;
}

edge_v::t_idx const * edge_v::mesh::Communication::getStruct( t_idx i_pa ) const {
  return m_chs+m_chOff[i_pa]*4 + i_pa;
}