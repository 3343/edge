/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2017, Regents of the University of California
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
 * Data layout of entities.
 **/

#include "io/logging.h"
#include "EntityLayout.h"

void edge::data::EntityLayout::sizesToLayout( t_enLayout &io_enLayout ) {
  // init the missing parameters
  io_enLayout.nEnts = 0;
  for( std::size_t l_tg = 0; l_tg < io_enLayout.timeGroups.size(); l_tg++ ) {
    io_enLayout.timeGroups[l_tg].nEntsOwn    = 0;
    io_enLayout.timeGroups[l_tg].nEntsNotOwn = 0;
    io_enLayout.timeGroups[l_tg].inner.first = 0;

    EDGE_CHECK( io_enLayout.timeGroups[l_tg].send.size() == io_enLayout.timeGroups[l_tg].receive.size() );

    for( unsigned int l_nr = 0; l_nr < io_enLayout.timeGroups[l_tg].send.size(); l_nr++ ) {
      io_enLayout.timeGroups[l_tg].send[l_nr].first    = 0;
      io_enLayout.timeGroups[l_tg].receive[l_nr].first = 0;
    }
  }

  // iterate over the time groups
  for( std::size_t l_tg = 0; l_tg < io_enLayout.timeGroups.size(); l_tg++ ) {
    // first element of inner entities based on previous time groups
    if( l_tg > 0 ) {
      io_enLayout.timeGroups[l_tg].inner.first += io_enLayout.timeGroups[l_tg-1].inner.first;
      io_enLayout.timeGroups[l_tg].inner.first += io_enLayout.timeGroups[l_tg-1].nEntsOwn;
      io_enLayout.timeGroups[l_tg].inner.first += io_enLayout.timeGroups[l_tg-1].nEntsNotOwn;
    }

    // add inner elements to owned elements
    io_enLayout.timeGroups[l_tg].nEntsOwn += io_enLayout.timeGroups[l_tg].inner.size;

    // iterate over send entities
    for( unsigned int l_nr = 0; l_nr < io_enLayout.timeGroups[l_tg].send.size(); l_nr++ ) {
      // first region is based on inner elements
      if( l_nr == 0 ) {
        io_enLayout.timeGroups[l_tg].send[0].first += io_enLayout.timeGroups[l_tg].inner.first;
        io_enLayout.timeGroups[l_tg].send[0].first += io_enLayout.timeGroups[l_tg].inner.size;
      }
      // all other regions rely on the previous one
      else {
        io_enLayout.timeGroups[l_tg].send[l_nr].first += io_enLayout.timeGroups[l_tg].send[l_nr-1].first;
        io_enLayout.timeGroups[l_tg].send[l_nr].first += io_enLayout.timeGroups[l_tg].send[l_nr-1].size;
      }

      // add send region to owned elements
      io_enLayout.timeGroups[l_tg].nEntsOwn += io_enLayout.timeGroups[l_tg].send[l_nr].size;
    }

    // iterate over receive entities
    for( unsigned int l_nr = 0; l_nr < io_enLayout.timeGroups[l_tg].receive.size(); l_nr++ ) {
      // first region is based on inner and send elements
      if( l_nr == 0 ) {
        io_enLayout.timeGroups[l_tg].receive[0].first += io_enLayout.timeGroups[l_tg].inner.first;
        io_enLayout.timeGroups[l_tg].receive[0].first += io_enLayout.timeGroups[l_tg].nEntsOwn;
      }
      // all other regions use the previous receive region
      else {
        io_enLayout.timeGroups[l_tg].receive[l_nr].first += io_enLayout.timeGroups[l_tg].receive[l_nr-1].first;
        io_enLayout.timeGroups[l_tg].receive[l_nr].first += io_enLayout.timeGroups[l_tg].receive[l_nr-1].size;
      }

      // add receive region to non-owned elements
      io_enLayout.timeGroups[l_tg].nEntsNotOwn += io_enLayout.timeGroups[l_tg].receive[l_nr].size;
    }

    // update total number of elements
    io_enLayout.nEnts += io_enLayout.timeGroups[l_tg].nEntsOwn;
    io_enLayout.nEnts += io_enLayout.timeGroups[l_tg].nEntsNotOwn;
  }
}
