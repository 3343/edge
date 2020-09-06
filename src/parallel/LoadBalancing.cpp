/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2020, Friedrich Schiller University Jena
 * Copyright (c) 2018-2020, Regents of the University of California
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
 * Dynamic load balancing of work regions.
 **/
#include <map>
#ifdef PP_USE_MPI
#include "mpi_wrapper.inc"
#endif

#include "LoadBalancing.h"

void edge::parallel::LoadBalancing::resolveSpEn( unsigned short i_id ) {
  // work packages of the region
  WrkPkgLb *l_wps = m_wrkRgns[i_id].wrkPkgs.data();

  // iterate over the defined sparse types
  for( unsigned short l_ty = 0; l_ty < m_wrkRgns[i_id].firstSp.size(); l_ty++ ) {
    // current relative position in the regions' sparse entities
    std::size_t l_sp = 0;

    // iterate over the work packages
    for( unsigned short l_wp = 0; l_wp < m_wrkRgns[i_id].wrkPkgs.size(); l_wp++ ) { 
      // assign first sparse id
      l_wps[l_wp].firstSp[l_ty] = m_wrkRgns[i_id].firstSp[l_ty] + l_sp;

      // count number of sparse entities
      while( l_sp < m_wrkRgns[i_id].spDe[l_ty].size() ) {
        if( m_wrkRgns[i_id].spDe[l_ty][l_sp] < l_wps[l_wp].first+l_wps[l_wp].size ) {
          l_sp++;
        }
        else break;
      }
    }

    // check that nobody is left behind
    EDGE_CHECK_EQ( l_sp, m_wrkRgns[i_id].spDe[l_ty].size() );
  }
}

void edge::parallel::LoadBalancing::init( unsigned int i_nWrks) {
  m_nWrks = i_nWrks;
}

void edge::parallel::LoadBalancing::balanceWrkRgn( unsigned short i_id ) {
  // work packages of the region
  WrkPkgLb *l_wps = m_wrkRgns[i_id].wrkPkgs.data();

  // gather elapsed times
  m_wrkRgns[i_id].elaMin = std::numeric_limits< double >::max();
  m_wrkRgns[i_id].elaMax = std::numeric_limits< double >::lowest();
  m_wrkRgns[i_id].elaSum = 0;

  std::vector< double > l_elapsed;
  l_elapsed.reserve( m_nWrks );
  for( unsigned int l_wo = 0; l_wo < m_nWrks; l_wo++ ) {
    // store the data and reset the timer
    l_elapsed.push_back( l_wps[l_wo].timer.elapsed() );
    l_wps[l_wo].timer.reset();

    // do the reductions
    m_wrkRgns[i_id].elaMin  = std::min( m_wrkRgns[i_id].elaMin, l_elapsed.back() );
    m_wrkRgns[i_id].elaMax  = std::max( m_wrkRgns[i_id].elaMax, l_elapsed.back() );
    m_wrkRgns[i_id].elaSum += l_elapsed.back();
  }

  // derive imbalance
  double l_imbalance  = (m_wrkRgns[i_id].elaMax - m_wrkRgns[i_id].elaMin);
         l_imbalance /= std::max( m_wrkRgns[i_id].elaSum / m_nWrks, m_zeroTime );

  // fill in pseudo-data if any of the elapsed times is non-positive
  if( m_wrkRgns[i_id].elaMin < m_zeroTime ) {
    for( unsigned int l_wo = 0; l_wo < m_nWrks; l_wo++ ) {
      l_elapsed[l_wo] = 3343;
      l_wps[l_wo].size = 1;
    }
  }
  // continue with current balancing if treshold is not reached
  else if( l_imbalance < m_maxImbalance ) return;

  // derive entity throughput per second
  std::vector< double > l_through;
  l_through.reserve( m_nWrks );
  double l_throughSum = 0;
  for( unsigned int l_wo = 0; l_wo < m_nWrks; l_wo++ ) {
    l_through.push_back( l_wps[l_wo].size / l_elapsed[l_wo] );
    l_throughSum += l_through.back();
  }

  // distribute the work by relative throughput
  std::size_t l_dist = 0;
  for( unsigned int l_wo = 0; l_wo < m_nWrks; l_wo++ ) {
    // continue for undefined
    if( l_throughSum == 0.0  ) {
      l_wps[l_wo].size = 0;
      continue;
    }

    // resize
    l_wps[l_wo].size = ( l_through[l_wo] / l_throughSum ) * m_wrkRgns[i_id].size;
    l_dist += l_wps[l_wo].size;
  }

  // remove entries, if too many have been distributed
  while( l_dist > m_wrkRgns[i_id].size ) {
    for( unsigned int l_wo = 0; l_wo < m_nWrks; l_wo++ ) {
      if( l_wps[l_wo].size > 0 ) {
        l_wps[l_wo].size--;
        l_dist--;
      }
      if( l_dist == m_wrkRgns[i_id].size ) break;
    }
  }

  // add entries, if too little have been distributed
  while( l_dist < m_wrkRgns[i_id].size ) {
    for( unsigned int l_wo = 0; l_wo < m_nWrks; l_wo++ ) {
      l_wps[l_wo].size++;
      l_dist++;
      if( l_dist == m_wrkRgns[i_id].size ) break;
    }
  }

  // derive the first entries
  l_wps[0].first = m_wrkRgns[i_id].first;
  for( unsigned int l_wo = 1; l_wo < m_nWrks; l_wo++ ) {
    l_wps[l_wo].first = l_wps[l_wo-1].first + l_wps[l_wo-1].size;
  }

  // double check that we got all entries
  EDGE_CHECK_EQ( l_wps[m_nWrks-1].first+l_wps[m_nWrks-1].size,
                 m_wrkRgns[i_id].first+m_wrkRgns[i_id].size );

  // assign sparse first ids
  resolveSpEn( i_id );
}

void edge::parallel::LoadBalancing::balance() {
  for( std::size_t l_rg = 0; l_rg < m_wrkRgns.size(); l_rg++ ) balanceWrkRgn( l_rg );
  m_nBalanced++;
}

void edge::parallel::LoadBalancing::startClock( unsigned short i_wrkRgn,
                                                unsigned short i_thread ) {
  m_wrkRgns[i_wrkRgn].wrkPkgs[i_thread].timer.start();
}

void edge::parallel::LoadBalancing::stopClock( unsigned short i_wrkRgn,
                                               unsigned short i_thread ) {
  m_wrkRgns[i_wrkRgn].wrkPkgs[i_thread].timer.end();
}

void edge::parallel::LoadBalancing::print() {
  EDGE_LOG_INFO << "printing statistics on data, used for load balancing #" << m_nBalanced;

  // collect stats
  std::vector< double > l_min, l_minG;
  std::vector< double > l_ave, l_aveG;
  std::vector< double > l_max, l_maxG;
  std::vector< double > l_imb, l_imbG;

  for( std::size_t l_wr = 0; l_wr < m_wrkRgns.size(); l_wr++ ) {
    l_min.push_back( m_wrkRgns[l_wr].elaMin );
    l_max.push_back( m_wrkRgns[l_wr].elaMax );
    l_ave.push_back( m_wrkRgns[l_wr].elaSum );
    double l_im = 0;
    if( l_min.back() > m_zeroTime )
      l_im = ( l_max.back()-l_min.back() ) / ( l_ave.back() / m_nWrks );
    l_imb.push_back( l_im );
  }

  l_minG.resize( l_min.size() );
  l_aveG.resize( l_ave.size() );
  l_maxG.resize( l_max.size() );
  l_imbG.resize( l_imb.size() );

  // derive global stats
#ifdef PP_USE_MPI
  MPI_Allreduce( l_min.data(), l_minG.data(), l_min.size(), MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );
  MPI_Allreduce( l_ave.data(), l_aveG.data(), l_ave.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
  MPI_Allreduce( l_max.data(), l_maxG.data(), l_max.size(), MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
  MPI_Allreduce( l_imb.data(), l_imbG.data(), l_imb.size(), MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
#else
  for( unsigned short l_wr = 0; l_wr < m_wrkRgns.size(); l_wr++ ) {
    l_minG[l_wr] = l_min[l_wr];
    l_aveG[l_wr] = l_ave[l_wr];
    l_maxG[l_wr] = l_max[l_wr];
    l_imbG[l_wr] = l_imb[l_wr];
  }
#endif

  // compute average
  for( unsigned short l_wr = 0; l_wr < m_wrkRgns.size(); l_wr++ ) {
    l_aveG[l_wr] /= m_nWrks * g_nRanks;
  }

  for( unsigned short l_wr = 0; l_wr < m_wrkRgns.size(); l_wr++ ) {
    EDGE_LOG_INFO << "  work region at position #" << l_wr << ":";
    EDGE_LOG_INFO << "    min time of any worker:       " << l_minG[l_wr] << "s";
    EDGE_LOG_INFO << "    ave time of all workers:      " << l_aveG[l_wr] << "s";
    EDGE_LOG_INFO << "    max time of any worker:       " << l_maxG[l_wr] << "s";
    EDGE_LOG_INFO << "    max imbalance over all ranks: " << l_imbG[l_wr];
  } 
}