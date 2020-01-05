/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2019-2020, Alexander Breuer
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
 * Basic interface for distributed memory parallelizations.
 **/
#include "Distributed.h"

void edge::parallel::Distributed::init( unsigned short         i_nTgs,
                                        unsigned short         i_nElFas,
                                        std::size_t            i_nEls,
                                        std::size_t            i_nByFa,
                                        std::size_t    const * i_commStruct,
                                        unsigned short const * i_sendFa,
                                        std::size_t    const * i_sendEl,
                                        unsigned short const * i_recvFa,
                                        std::size_t    const * i_recvEl,
                                        data::Dynamic        & io_dynMem ) {
  m_nTgs = i_nTgs;

  // derive the number of communication channels and communicating faces
  if( i_commStruct != nullptr ) {
    m_nChs = i_commStruct[0];
  }
  m_nSeRe = (std::size_t *) io_dynMem.allocate( m_nChs * sizeof(std::size_t) );

  std::size_t l_sizeSend = 0;
  std::size_t l_sizeRecv = 0;
  for( std::size_t l_ch = 0; l_ch < m_nChs; l_ch++ ) {
    std::size_t l_tg    = i_commStruct[1 + l_ch*4 + 0];
    std::size_t l_tgAd  = i_commStruct[1 + l_ch*4 + 2];
    std::size_t l_nSeRe = i_commStruct[1 + l_ch*4 + 3];
    m_nSeRe[l_ch] = l_nSeRe;

    l_sizeSend += (l_tg > l_tgAd) ? 2 * l_nSeRe * i_nByFa : l_nSeRe * i_nByFa;
    l_sizeRecv += (l_tg < l_tgAd) ? 2 * l_nSeRe * i_nByFa : l_nSeRe * i_nByFa;
  }

  // allocate send and receive buffer
  l_sizeSend *= sizeof(unsigned char);
  l_sizeRecv *= sizeof(unsigned char);
  m_sendBufferSize = l_sizeSend;
  m_sendBuffer = (unsigned char*) io_dynMem.allocate( l_sizeSend );
  m_recvBufferSize = l_sizeRecv;
  m_recvBuffer = (unsigned char*) io_dynMem.allocate( l_sizeRecv );

  // allocate pointer data structure
  m_sendPtrs = (unsigned char**) io_dynMem.allocate( i_nEls * i_nElFas * sizeof(unsigned char*) );
  m_recvPtrs = (unsigned char**) io_dynMem.allocate( i_nEls * i_nElFas * sizeof(unsigned char*) );

  // init with null pointers (no communication)
#ifdef PP_USE_OMP
#pragma omp parallel for
#endif
  for( std::size_t l_el = 0; l_el < i_nEls; l_el++ ) {
    for( unsigned short l_fa = 0; l_fa < i_nElFas; l_fa++ ) {
      m_sendPtrs[ l_el*i_nElFas + l_fa ] = nullptr;
      m_recvPtrs[ l_el*i_nElFas + l_fa ] = nullptr;
    }
  }

  // init message structure and communication data
  m_sendMsgs = (t_msg*) io_dynMem.allocate( m_nChs * sizeof(t_msg) );
  m_recvMsgs = (t_msg*) io_dynMem.allocate( m_nChs * sizeof(t_msg) );

  std::size_t l_offSend = 0;
  std::size_t l_offRecv = 0;
  std::size_t l_first = 0;
  for( std::size_t l_ch = 0; l_ch < m_nChs; l_ch++ ) {
    // unpack channel: local time group, adjacent rank, adjacent time group, #comm faces
    std::size_t l_tg    = i_commStruct[1 + l_ch*4 + 0];
    std::size_t l_raAd  = i_commStruct[1 + l_ch*4 + 1];
    std::size_t l_tgAd  = i_commStruct[1 + l_ch*4 + 2];
    std::size_t l_nSeRe = i_commStruct[1 + l_ch*4 + 3];

    l_sizeSend = (l_tg > l_tgAd) ? l_nSeRe * i_nByFa * 2 : l_nSeRe * i_nByFa;
    l_sizeRecv = (l_tg < l_tgAd) ? l_nSeRe * i_nByFa * 2 : l_nSeRe * i_nByFa;

    // assign
    m_sendMsgs[l_ch].lt      = (l_tg < l_tgAd);
    m_sendMsgs[l_ch].tg      = l_tg;
    m_sendMsgs[l_ch].rank    = l_raAd;
    m_sendMsgs[l_ch].size    = l_sizeSend;
    m_sendMsgs[l_ch].offL    = l_offSend;

    m_recvMsgs[l_ch].lt      = (l_tg < l_tgAd);
    m_recvMsgs[l_ch].tg      = l_tg;
    m_recvMsgs[l_ch].rank    = l_raAd;
    m_recvMsgs[l_ch].size    = l_sizeRecv;
    m_recvMsgs[l_ch].offL    = l_offRecv;

    for( std::size_t l_co = 0; l_co < l_nSeRe; l_co++ ) {
      std::size_t l_seFa = i_sendFa[l_first+l_co];
      std::size_t l_seEl = i_sendEl[l_first+l_co];
      std::size_t l_reFa = i_recvFa[l_first+l_co];
      std::size_t l_reEl = i_recvEl[l_first+l_co];

      m_sendPtrs[ l_seEl*i_nElFas + l_seFa ] = m_sendBuffer+l_offSend;
      m_recvPtrs[ l_reEl*i_nElFas + l_reFa ] = m_recvBuffer+l_offRecv;

      l_offSend += (l_tg > l_tgAd) ? i_nByFa * 2 : i_nByFa;
      l_offRecv += (l_tg < l_tgAd) ? i_nByFa * 2 : i_nByFa;
    }
    l_first += l_nSeRe;
  }
}