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
 * Derives mesh partitions.
 **/
#include "Partition.h"
#include <metis.h>

void edge_v::mesh::Partition::getElPr( edge_v::t_entityType         i_elTy,
                                       std::size_t                  i_nEls,
                                       std::size_t          const * i_elFaEl,
                                       std::size_t          const * i_elPa,
                                       unsigned short       const * i_elTg,
                                       std::size_t                * o_elPr ) {
  // determine number of time groups
  unsigned short l_nTgs = 0;
  for( std::size_t l_el = 0; l_el < i_nEls; l_el++ ) {
    l_nTgs = std::max( l_nTgs, i_elTg[l_el] );
  }
  l_nTgs++;

  // lambda, which determines if an element is a send-element
  auto l_send = [ i_elTy, i_elFaEl, i_elPa ]( std::size_t i_el ) {
    unsigned short l_nElFas = CE_N_FAS( i_elTy );
    bool l_isSend = false;

    // element's partition
    std::size_t l_elPa = i_elPa[i_el];

    for( unsigned short l_fa = 0; l_fa < l_nElFas; l_fa++ ) {
      std::size_t l_ad = i_elFaEl[ i_el*l_nElFas + l_fa ];
      if( l_ad != std::numeric_limits< std::size_t >::max() ) {
        // neighboring element's partition
        std::size_t l_adPa = i_elPa[l_ad];
        if( l_elPa != l_adPa ) l_isSend = true;
      }
    }
    return l_isSend;
  };

#ifdef PP_USE_OMP
#pragma omp parallel for
#endif
  for( std::size_t l_el = 0; l_el < i_nEls; l_el++ ) {
    // init priority based on element's rank
    o_elPr[l_el] = i_elPa[l_el]*2*l_nTgs;

    // add offset of send elements
    if( l_send(l_el) ) o_elPr[l_el] += l_nTgs;

    // add time group priority
    o_elPr[l_el] += i_elTg[l_el];
  }
}

edge_v::mesh::Partition::~Partition() {
  // free memory
  if( m_elPa != nullptr ) {
    delete[] m_elPa;
  }
}

void edge_v::mesh::Partition::kWay( std::size_t            i_nParts,
                                    unsigned short const * i_elTg,
                                    unsigned short         i_nCuts ) {
  // get info from mesh
  edge_v::t_entityType l_elTy = m_mesh.getTypeEl();
  unsigned short l_nElFas = CE_N_FAS( l_elTy );
  std::size_t l_nEls = m_mesh.nEls();
  std::size_t const * l_elFaEl = m_mesh.getElFaEl();

  // assemble the adjacency info for the metis-call
  idx_t * l_xadj = new idx_t[ l_nEls+1 ];
  idx_t * l_adjncy = new idx_t[ l_nEls*l_nElFas ];

  std::size_t l_adId = 0;
  l_xadj[0] = 0;
  for( std::size_t l_el = 0; l_el < l_nEls; l_el++ ) {
    l_xadj[l_el+1] = l_xadj[l_el];

    for( std::size_t l_fa = 0; l_fa < l_nElFas; l_fa++ ) {
      std::size_t l_ad = l_elFaEl[l_el*l_nElFas + l_fa];

      if( l_ad != std::numeric_limits< std::size_t >::max() ) {
        l_adjncy[l_adId] = l_ad;
        l_adId++;
        l_xadj[l_el+1]++;
      }
    }

    // sort by element id
    std::sort( l_adjncy+l_xadj[l_el], l_adjncy+l_xadj[l_el+1] );
  }

  // assemble vertex and edge weights
  idx_t * l_vwgt = nullptr;
  idx_t * l_adjwgt = nullptr;
  if( i_elTg != nullptr ) {
    l_vwgt = new idx_t[ l_nEls ];
    l_adjwgt = new idx_t[ l_xadj[l_nEls] ];

    unsigned short l_tgMax = 0;
    for( std::size_t l_el = 0; l_el < l_nEls; l_el++ ) {
      l_tgMax = std::max( l_tgMax, i_elTg[l_el] );
      l_vwgt[l_el] = 1;
    }

    for( idx_t l_ad = 0; l_ad < l_xadj[l_nEls]; l_ad++ )
      l_adjwgt[l_ad] = 1;

    std::size_t l_adId = 0;
    for( std::size_t l_el = 0; l_el < l_nEls; l_el++ ) {
      // set vertex weight
      for( unsigned short l_tg = i_elTg[l_el]; l_tg < l_tgMax; l_tg++ ) {
        l_vwgt[l_el] *= 2;
      }

      for( std::size_t l_fa = 0; l_fa < l_nElFas; l_fa++ ) {
        std::size_t l_ad = l_elFaEl[l_el*l_nElFas + l_fa];

        if( l_ad != std::numeric_limits< std::size_t >::max() ) {
          // larger time group elements have to sent twice the amount
          // -> comm volume is given by frequency of min time group
          unsigned short l_minTg = std::min( i_elTg[l_el], i_elTg[l_ad] );

          // set edge weight
          for( unsigned short l_tg = l_minTg; l_tg < l_tgMax; l_tg++ ) {
            l_adjwgt[l_adId] *= 2;
          }

          l_adId++;
        }
      }
    }
  }

  // set remaining metis parameters
  idx_t l_nvtxs = l_nEls;
  idx_t l_ncon = 1;
  idx_t l_objVal = 0;
  idx_t * l_elPa = new idx_t[ l_nEls ];
  idx_t l_nParts = i_nParts;
  idx_t l_opts[METIS_NOPTIONS];
  int l_err = METIS_SetDefaultOptions(l_opts);
  EDGE_V_CHECK_EQ( l_err, METIS_OK );
  l_opts[METIS_OPTION_CONTIG] = 1;
  l_opts[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_VOL;
  l_opts[METIS_OPTION_NCUTS] = i_nCuts;

  // call metis
  l_err = METIS_PartGraphKway( &l_nvtxs,
                               &l_ncon,
                                l_xadj,
                                l_adjncy,
                                l_vwgt,
                                NULL,
                                l_adjwgt,
                                &l_nParts,
                                NULL,
                                NULL,
                                l_opts,
                                &l_objVal,
                                l_elPa );
  EDGE_V_CHECK_EQ( l_err, METIS_OK );

  // free intermediate adjacency memory
  if( i_elTg != nullptr ) {
    delete[] l_vwgt;
    delete[] l_adjwgt;
  }
  delete[] l_adjncy;
  delete[] l_xadj;

  // allocate memory for the results
  if( m_elPa == nullptr ) {
    m_elPa = new std::size_t[ l_nEls ];
  }

  // store results
#ifdef PP_USE_OMP
#pragma omp parallel for
#endif
  for( std::size_t l_el = 0; l_el < l_nEls; l_el++ ) {
    m_elPa[l_el] = l_elPa[l_el];
  }

  // free intermediate partition storage
  delete[] l_elPa;
}