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
 * Time step groups.
 **/
#include "Groups.h"

#include "io/logging.h"

void edge_v::time::Groups::getLoads( std::size_t            i_nEls,
                                     double         const * i_tsCfl,
                                     double         const * i_tsGroups,
                                     unsigned short const * i_elTg,
                                     double               & o_gts,
                                     double               & o_ltsGrouped,
                                     double               & o_ltsPerElement ) {
  // gts is scaled by 1
  o_gts = i_nEls;

  o_ltsGrouped = 0;
  o_ltsPerElement = 0;

  for( std::size_t l_el = 0; l_el < i_nEls; l_el++ ) {
    unsigned short l_tg = i_elTg[l_el];
    o_ltsGrouped += 1.0 / i_tsGroups[l_tg];
    o_ltsPerElement += 1.0 / i_tsCfl[l_el];
  }
}

std::size_t edge_v::time::Groups::normalizeElTgs( t_entityType           i_elTy,
                                                  std::size_t            i_nEls,
                                                  std::size_t    const * i_elFaEl,
                                                  unsigned short       * io_elTg ) {
  std::size_t l_normAll = 0;

  unsigned short l_nElFas = CE_N_FAS( i_elTy );

  while( true ) {
    std::size_t l_normIt = 0;

    for( std::size_t l_el = 0; l_el < i_nEls; l_el++ ) {
      unsigned short l_tgEl = io_elTg[l_el];

      for( unsigned short l_fa = 0; l_fa < l_nElFas; l_fa++ )  {
        std::size_t l_ad = i_elFaEl[l_el * l_nElFas + l_fa];

        // ignore faces at the boundary
        if( l_ad < std::numeric_limits< std::size_t >::max() ) {
          unsigned short l_tgAd = io_elTg[l_ad];
          // lower time group if required
          if( l_tgEl > l_tgAd+1 ) {
            io_elTg[l_el] = l_tgAd+1;
            l_normIt++;
            l_normAll++;
          }
        }
      }
    }

    if( l_normIt == 0 ) break;
  }

  return l_normAll;
}


void edge_v::time::Groups::setElTg( std::size_t            i_nEls,
                                    unsigned short         i_nGroups,
                                    double         const * i_tsGroups,
                                    double         const * i_tsCfl,
                                    unsigned short       * o_elTg ) {
#ifdef PP_USE_OMP
#pragma omp parallel for
#endif
  for( std::size_t l_el = 0; l_el < i_nEls; l_el++ ) {
    for( unsigned short l_tg = 0; l_tg < i_nGroups; l_tg++ ) {
      if( i_tsCfl[l_el] < i_tsGroups[l_tg+1] ) {
        o_elTg[l_el] = l_tg;
        break;
      }
    }
  }
}

void edge_v::time::Groups::nGroupEls( std::size_t            i_nEls,
                                      unsigned short         i_nGroups,
                                      unsigned short const * i_elTg,
                                      std::size_t          * o_nGroupEls ) {
  for( unsigned short l_tg = 0; l_tg < i_nGroups; l_tg++ )
    o_nGroupEls[l_tg] = 0;

  for( std::size_t l_el = 0; l_el < i_nEls; l_el++ ) {
    unsigned short l_tg = i_elTg[l_el];
    o_nGroupEls[l_tg]++;
  }
}


edge_v::time::Groups::Groups( t_entityType           i_elTy,
                              std::size_t            i_nEls,
                              std::size_t    const * i_elFaEl,
                              unsigned short         i_nRates,
                              double         const * i_rates,
                              double                 i_funDt,
                              double         const * i_ts ) {
  // check for valid rates and minimum relative time step
  for( unsigned short l_ra = 0; l_ra < i_nRates; l_ra++ ) {
    EDGE_V_CHECK_GT( i_rates[l_ra], 1.0 );
    EDGE_V_CHECK_LE( i_rates[l_ra], 2.0 );
  }
  EDGE_V_CHECK_GT( i_funDt, 0 );
  EDGE_V_CHECK_LE( i_funDt, 1 );

  m_nEls = i_nEls;

  // allocate memory for the groups and init
  m_nGroups = i_nRates+1;
  m_tsIntervals = new double[ m_nGroups+1 ];
  m_tsIntervals[0] = i_funDt;
  m_tsIntervals[i_nRates+1] = std::numeric_limits< double >::max();
  for( unsigned short l_ra = 0; l_ra < i_nRates; l_ra++ ) {
    m_tsIntervals[l_ra+1] = m_tsIntervals[l_ra] * i_rates[l_ra];
  }

  // allocate memory for the element-to-group assignment and init
  m_elTg = new unsigned short[i_nEls];
  setElTg( i_nEls,
           m_nGroups,
           m_tsIntervals,
           i_ts,
           m_elTg );

  // store loads excluding normalization
  getLoads( m_nEls,
            i_ts,
            m_tsIntervals,
            m_elTg,
            m_loads[0],
            m_loads[1],
            m_loads[3] );

  // normalize time group assignments
  normalizeElTgs( i_elTy,
                  i_nEls,
                  i_elFaEl,
                  m_elTg );

  // derive number of elements per time group
  m_nGroupEls = new std::size_t[m_nGroups];
  nGroupEls( i_nEls,
             m_nGroups,
             m_elTg,
             m_nGroupEls );

  // store loads including normalization
  getLoads( m_nEls,
            i_ts,
            m_tsIntervals,
            m_elTg,
            m_loads[0],
            m_loads[2],
            m_loads[3] );
}

edge_v::time::Groups::~Groups() {
  // free memory
  delete[] m_tsIntervals;
  delete[] m_elTg;
  delete[] m_nGroupEls;
}

void edge_v::time::Groups::printStats() const {
  EDGE_V_LOG_INFO << "time step histogram (group / #elements / range):";
  for( unsigned short l_tg = 0; l_tg < m_nGroups; l_tg++ ) {
    EDGE_V_LOG_INFO << "  " << l_tg
                            << " " << m_nGroupEls[l_tg]
                            << " [" << m_tsIntervals[l_tg] << ", " << m_tsIntervals[l_tg+1] << "[";
  }

  EDGE_V_LOG_INFO << "theoretical LTS speedups:";
  EDGE_V_LOG_INFO << "  grouped              over GTS:                  " << m_loads[0] / m_loads[1];
  EDGE_V_LOG_INFO << "  grouped (normalized) over GTS:                  " << m_loads[0] / m_loads[2];
  EDGE_V_LOG_INFO << "  per-element          over GTS:                  " << m_loads[0] / m_loads[3];
  EDGE_V_LOG_INFO << "  per-element          over grouped:              " << m_loads[1] / m_loads[3];
  EDGE_V_LOG_INFO << "  per-element          over grouped (normalized): " << m_loads[2] / m_loads[3];
}