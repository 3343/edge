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
 * Derives communication structures.
 **/
#include "Communication.h"
#include "io/logging.h"

bool edge_v::mesh::Communication::isComm( t_entityType         i_elTy,
                                          std::size_t          i_el,
                                          std::size_t  const * i_elFaEl,
                                          std::size_t  const * i_elPa ) {
  unsigned short l_nElFas = CE_N_FAS( i_elTy );
  bool l_isComm = false;

  // element's partition
  std::size_t l_elPa = i_elPa[i_el];

  for( unsigned short l_fa = 0; l_fa < l_nElFas; l_fa++ ) {
    std::size_t l_ad = i_elFaEl[ i_el*l_nElFas + l_fa ];
    if( l_ad != std::numeric_limits< std::size_t >::max() ) {
      // neighboring element's partition
      std::size_t l_adPa = i_elPa[l_ad];
      if( l_elPa != l_adPa ) l_isComm = true;
    }
  }
  return l_isComm;
};

void edge_v::mesh::Communication::getPaElComm( t_entityType        i_elTy,
                                               std::size_t         i_nEls,
                                               std::size_t const * i_elFaEl,
                                               std::size_t const * i_elPa,
                                               std::size_t       * o_first,
                                               std::size_t       * o_size ) {
  // partition of the communication region
  std::size_t l_regPa = std::numeric_limits< std::size_t >::max();

  for( std::size_t l_el = 0; l_el < i_nEls; l_el++ ) {
    bool l_commEl = isComm( i_elTy,
                            l_el,
                            i_elFaEl,
                            i_elPa );

    std::size_t l_elPa = i_elPa[l_el];

    // init a new comm region if a new comm region shows up
    if( l_regPa != l_elPa && l_commEl == true  ) {
      // sanity check on the element order
      EDGE_V_CHECK( l_regPa == std::numeric_limits< std::size_t >::max() ||
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

void edge_v::mesh::Communication::getPaTgPairs( t_entityType                           i_elTy,
                                                std::size_t                            i_first,
                                                std::size_t                            i_size,
                                                std::size_t                    const * i_elFaEl,
                                                unsigned short                 const * i_elTg,
                                                std::size_t                    const * i_elPa,
                                                std::set<
                                                std::pair< std::size_t,
                                                            unsigned short > >       & o_pairs ) {
  unsigned short l_nElFas = CE_N_FAS( i_elTy );
  EDGE_V_CHECK_GT( i_size, 0 );

  // reset pairs
  o_pairs.empty();

  // get the region's partition and time group
  unsigned short l_pa = i_elPa[i_first];
  unsigned short l_tg = i_elTg[i_first];

  // iterate over the elements of the region
  for( std::size_t l_el = i_first; l_el < i_first+i_size; l_el++ ) {
    // check that the partition and time group match
    EDGE_V_CHECK_EQ( i_elPa[l_el], l_pa );
    EDGE_V_CHECK_EQ( i_elTg[l_el], l_tg );

    // iterate over the face-adjacent elements
    for( unsigned short l_fa = 0; l_fa < l_nElFas; l_fa++ ) {
      std::size_t l_ad = i_elFaEl[ l_el*l_nElFas + l_fa ];
      if( l_ad != std::numeric_limits< std::size_t >::max() ) {
        // get adjacent partition and time group
        std::size_t l_adPa = i_elPa[l_ad];
        std::size_t l_adTg = i_elTg[l_ad];

        // communication is required, if adjacent elements belongs to a different partition
        if( l_pa != l_adPa ) {
          std::pair< std::size_t, unsigned short > l_pair( {l_adPa, l_adTg} );
          o_pairs.insert( l_pair );
        }
      }
    }
  }
}

void edge_v::mesh::Communication::getMsgsSend( t_entityType                   i_elTy,
                                               std::size_t                    i_first,
                                               std::size_t                    i_size,
                                               std::size_t            const * i_elFaEl,
                                               std::size_t            const * i_elPa,
                                               unsigned short         const * i_elTg,
                                               std::vector< Message >       & o_send ) {
  // get the communication partition-time group pairs
  std::set< std::pair< std::size_t,
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
  std::size_t l_msg = 0;
  for( auto l_it = l_pairs.begin(); l_it != l_pairs.end(); l_it++ ) {
    o_send[l_msg].pa = std::get<0>(*l_it);
    o_send[l_msg].tg = std::get<1>(*l_it);
    l_msg++;
  }

  // get the communication region's partition and time region
  std::size_t l_paReg = i_elPa[i_first];

  unsigned short l_nElFas = CE_N_FAS( i_elTy );

  // lambda, which inserts an element-face pair into the send messages
  auto l_insert = [   i_elFaEl,
                      l_nElFas,
                    & o_send    ]( std::size_t    i_el,
                                   unsigned short i_fa,
                                   std::size_t    i_pa,
                                   unsigned short i_tg ) mutable {
    for( std::size_t l_se = 0; l_se < o_send.size(); l_se++ ) {
      if( o_send[l_se].pa == i_pa && o_send[l_se].tg == i_tg ) {
        o_send[l_se].el.push_back( i_el );
        o_send[l_se].fa.push_back( i_fa );

        std::size_t l_elAd = i_elFaEl[i_el*l_nElFas + i_fa];
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
  for( std::size_t l_el = i_first; l_el < i_first+i_size; l_el++ ) {
    for( unsigned short l_fa = 0; l_fa < l_nElFas; l_fa++ ) {
      std::size_t l_ad = i_elFaEl[l_el*l_nElFas + l_fa];

      if( l_ad != std::numeric_limits< std::size_t >::max() ) {
        std::size_t l_paAd = i_elPa[l_ad];
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
                                             std::size_t                    i_nEls,
                                             std::size_t            const * i_elFaEl,
                                             std::size_t            const * i_elPa,
                                             unsigned short         const * i_elTg,
                                             std::vector< Partition >     & o_struct ) {
  std::size_t l_nPas = i_elPa[i_nEls-1]+1;
  std::size_t *l_firstComm = new std::size_t[ l_nPas ];
  std::size_t *l_sizeComm  = new std::size_t[ l_nPas ];

  // get communicating elements
  getPaElComm( i_elTy,
               i_nEls,
               i_elFaEl,
               i_elPa,
               l_firstComm,
               l_sizeComm );

  // lambda which gets the time groups of a communication region
  auto l_getTgs = [ i_elTg ]( std::size_t                     i_first,
                              std::size_t                     i_size,
                              std::vector< unsigned short > & o_tgs,
                              std::vector< std::size_t >    & o_firstTg,
                              std::vector< std::size_t >    & o_sizeTg ) {
    o_tgs.resize(0);
    o_firstTg.resize(0);
    o_sizeTg.resize(0);
    if( i_size == 0 ) return;

    unsigned short l_tg = std::numeric_limits< unsigned short >::max();
    for( std::size_t l_el = i_first; l_el < i_first+i_size; l_el++ ) {
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
    std::vector< std::size_t > l_firstTg;
    std::vector< std::size_t > l_sizeTg;
    std::vector< unsigned short > l_tgs;

    l_getTgs( l_firstComm[l_pa],
              l_sizeComm[l_pa],
              l_tgs,
              l_firstTg,
              l_sizeTg );

    o_struct[l_pa].tr.resize( l_tgs.size() );

    for( std::size_t l_tg = 0; l_tg < l_tgs.size(); l_tg++ ) {
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
  for( std::size_t l_pa = 0; l_pa < l_nPas; l_pa++ ) {
    for( std::size_t l_t0 = 0; l_t0 < o_struct[l_pa].tr.size(); l_t0++ ) {
      std::size_t l_nSend = o_struct[l_pa].tr[l_t0].send.size();
      o_struct[l_pa].tr[l_t0].recv.resize( l_nSend );

      // iterate over the send messages
      for( std::size_t l_s0 = 0; l_s0 < l_nSend; l_s0++ ) {
        unsigned short l_tgAd = o_struct[l_pa].tr[l_t0].send[l_s0].tg;
        std::size_t l_paAd = o_struct[l_pa].tr[l_t0].send[l_s0].pa;

        // find the matching send of the remote partition which reflects our recv
        for( std::size_t l_t1 = 0; l_t1 < o_struct[l_paAd].tr.size(); l_t1++ ) {
          if( o_struct[l_paAd].tr[l_t1].tg == l_tgAd ) {
            for( std::size_t l_s1 = 0; l_s1 < o_struct[l_paAd].tr[l_t1].send.size(); l_s1++ ) {
              if( o_struct[l_paAd].tr[l_t1].send[l_s1].pa == l_pa ) {
                o_struct[l_pa].tr[l_t0].recv[l_s0] = o_struct[l_paAd].tr[l_t1].send[l_s1];
                // reverse order of local and remote
                std::vector< std::size_t > l_tmpEl = o_struct[l_pa].tr[l_t0].recv[l_s0].el;
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
                                            std::size_t                    * o_chOff ) {
  o_chOff[0] = 0;
  for( std::size_t l_pa = 0; l_pa < i_struct.size(); l_pa++ ) {
    o_chOff[l_pa+1] = o_chOff[l_pa];

    for( unsigned short l_tg = 0; l_tg < i_struct[l_pa].tr.size(); l_tg++ ) {
      o_chOff[l_pa+1] += i_struct[l_pa].tr[l_tg].send.size();
    }
  }
}

void edge_v::mesh::Communication::setChs( std::vector< Partition > const & i_struct,
                                          std::size_t              const * i_chOff,
                                          std::size_t                    * o_chs ) {
  for( std::size_t l_pa = 0; l_pa < i_struct.size(); l_pa++ ) {
    std::size_t l_off = i_chOff[l_pa]*4 + l_pa;
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

edge_v::mesh::Communication::Communication( t_entityType           i_elTy,
                                            std::size_t            i_nEls,
                                            std::size_t    const * i_elFaEl,
                                            std::size_t            i_nPas,
                                            std::size_t    const * i_nPaEls,
                                            unsigned short const * i_elTg ) {
  EDGE_V_CHECK_GT( i_nPas, 0 );

  // init temporary data structure for the parts
  std::size_t * l_elPa = new std::size_t[ i_nEls ];
  std::size_t l_paFirst = 0;

  for( std::size_t l_pa = 0; l_pa < i_nPas; l_pa++ ) {
    std::size_t l_paSize = i_nPaEls[l_pa];
    EDGE_V_CHECK_LE( l_paFirst+l_paSize, i_nEls );
    for( std::size_t l_el = l_paFirst; l_el < l_paFirst+l_paSize; l_el++ ) {
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

  delete[] l_elPa;

  // get channels offsets
  m_chOff = new std::size_t[i_nPas+1];
  setChOff( m_struct,
            m_chOff );

  // assign simplified comm structure
  m_chs = new std::size_t[ i_nPas + m_chOff[i_nPas]*4 ];
  setChs( m_struct,
          m_chOff,
          m_chs );
}

edge_v::mesh::Communication::~Communication() {
  delete[] m_chOff;
  delete[] m_chs;
}

std::size_t const * edge_v::mesh::Communication::getStruct( std::size_t i_pa ) const {
  return m_chs+m_chOff[i_pa]*4 + i_pa;
}