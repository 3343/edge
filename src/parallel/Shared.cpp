/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2015-2018, Regents of the University of California
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
 * Shared memory parallelization.
 **/
#include "Shared.h"
#include "io/logging.h"
#include <atomic>

#ifdef PP_USE_OMP
#include <omp.h>
#endif

void edge::parallel::Shared::init( unsigned int i_nWrks ) {
#ifdef PP_USE_OMP
#pragma omp parallel
{
  g_thread = omp_get_thread_num();
  std::string l_threadStr = std::to_string(g_thread);
  l_threadStr.copy( g_threadStr, 10 );

  // do the initialization of shared vars only on thread 0
  if( g_thread == 0 ) {
    g_nThreads = omp_get_num_threads();

    // assign worker threads
    if( i_nWrks == 0 ) {
#ifdef PP_USE_MPI
      m_nWrks = std::max( g_nThreads - 1, 1);
#else
      m_nWrks = g_nThreads;
#endif
    }
    else m_nWrks = i_nWrks;

    // check for at least one worker
    EDGE_CHECK( m_nWrks > 0 );
  }

  // check that no other threads wrote in private thread num
  // no guarantee that this check will discover this however
  EDGE_CHECK( g_thread == omp_get_thread_num() );
}
#else
  g_thread   = 0;
  g_nThreads = 1;
  m_nWrks    = 1;
#endif

  // init the load balancing
  m_balancing.init( m_nWrks );
}

void edge::parallel::Shared::print() {
  EDGE_LOG_INFO << "shared memory setup:";
#ifdef PP_USE_OMP
  EDGE_LOG_INFO << "  omp-version: "
#if _OPENMP==200505
                << "2.5";
#elif _OPENMP==200805
                << "3.0";
#elif _OPENMP==201107
                << "3.1";
#elif _OPENMP==201307
                << "4.0";
#elif _OPENMP==201511
                << "4.5";
#else
                << "unknown, raw is " << _OPENMP;
#endif

#endif
  EDGE_LOG_INFO << "  #threads: " << g_nThreads;
  EDGE_LOG_INFO << "  #workers: " << m_nWrks;
}

bool edge::parallel::Shared::isCommLead() {
  if( m_nWrks == g_nThreads ) return g_thread==0;
  else                        return g_thread==m_nWrks;
}

bool edge::parallel::Shared::isWrk() {
  return g_thread < m_nWrks;
}

bool edge::parallel::Shared::isSched() {
  EDGE_CHECK( m_nWrks > 0 );
  // if every thread is a worker, the last one schedules
  if( m_nWrks == g_nThreads ) return g_thread == m_nWrks-1;
  // if scheduling and communication threads are available, the first one of these schedules work
  else                        return g_thread == m_nWrks;
}

bool edge::parallel::Shared::isComm() {
  EDGE_CHECK( m_nWrks > 0 );
  // if every thread is a worker, every thread communicates
  if( m_nWrks == g_nThreads ) return true;
  // if scheduling and comm threads are available, every thread but the workers communicates
  else                        return g_thread >= m_nWrks;
}

std::size_t edge::parallel::Shared::getWrkRgn( unsigned int i_id ) {
  std::size_t l_rg;
  for( l_rg = 0; l_rg < m_wrkRgns.size(); l_rg++ ) {
    if( m_wrkRgns[l_rg].id == i_id ) break;
    EDGE_CHECK( l_rg != m_wrkRgns.size()-1 );
  }
  return l_rg;
}

bool edge::parallel::Shared::getWrkTd( int_tg         & o_tg,
                                       unsigned short & o_step,
                                       unsigned int   & o_id,
                                       int_el         & o_first,
                                       int_el         & o_size,
                                       int_el         * o_firstSp ) {
  EDGE_CHECK( g_thread < m_nWrks );

  // invalid everything by default
  o_tg        = std::numeric_limits< int_tg >::max();
  o_step      = std::numeric_limits< unsigned short >::max();
  o_id        = std::numeric_limits< unsigned int >::max();
  o_first     = std::numeric_limits< int_el       >::max();
  o_size      = std::numeric_limits< int_el       >::max();

  // iterate over work region
  for( std::size_t l_rg = 0; l_rg < m_wrkRgns.size(); l_rg++ ) {
    volatile WrkPkg* l_wps = &m_wrkRgns[l_rg].wrkPkgs[0];

    if( l_wps[g_thread].status == RDY ) {
      o_tg    = m_wrkRgns[l_rg].tg;
      o_step  = m_wrkRgns[l_rg].step;
      o_id    = m_wrkRgns[l_rg].id;

      m_balancing.getWrkTd( l_rg,
                            g_thread,
                            o_first,
                            o_size,
                            o_firstSp );

      // flush for a consistent view
#ifdef PP_USE_OMP
#pragma omp flush
#endif

      return true;
    }
  }

  // nothing found return false
  return false;
}

void edge::parallel::Shared::setStatusAll( t_status     i_status,
                                           unsigned int i_id ) {
  // flush for a consistent view
#ifdef PP_USE_OMP
#pragma omp flush
#endif

  std::size_t l_rg = getWrkRgn( i_id );

  volatile WrkPkg* l_wps = &m_wrkRgns[l_rg].wrkPkgs[0];

  // iterate over all workers and set status
  for( int l_td = 0; l_td < m_nWrks; l_td++ ) {
    if( i_status == RDY ) {
      EDGE_CHECK( l_wps[l_td].status == FIN || l_wps[l_td].status == WAI );
    }
    else {
      EDGE_LOG_FATAL << "status change not allowed";
    }

    // set status
    l_wps[l_td].status = i_status;
  }
}

void edge::parallel::Shared::resetStatus( t_status i_status ) {
  // flush for a consistent view
#ifdef PP_USE_OMP
#pragma omp flush
#endif

  // iterate over all regions
  for( unsigned short l_rg = 0; l_rg < m_wrkRgns.size(); l_rg++ ) {
    volatile WrkPkg* l_wps = &m_wrkRgns[l_rg].wrkPkgs[0];

    // iterate over all workers
    for( int l_td = 0; l_td < m_nWrks; l_td++ ) {
      // set status
      l_wps[l_td].status = i_status;
    }
  }
}

void edge::parallel::Shared::setStatusTd(  t_status     i_status,
                                           unsigned int i_id ) {
  // flush for a consistent view
#ifdef PP_USE_OMP
#pragma omp flush
#endif

  // check that the calling thread is a worker
  EDGE_CHECK( g_thread < m_nWrks );

  std::size_t l_rg = getWrkRgn( i_id );

  volatile t_status *l_st = &m_wrkRgns[l_rg].wrkPkgs[g_thread].status;

  // check that the previous status matches
  if( i_status == IPR ) {
    EDGE_CHECK( *l_st == RDY );
    // start the timer for this work package, now having status "in progress"
    m_balancing.startClock( l_rg, g_thread );
  }
  else if( i_status == FIN ) {
    EDGE_CHECK( *l_st == IPR );
    // stop the timer for this work package, now having status "finished"
    m_balancing.stopClock( l_rg, g_thread );
  }
  else {
    EDGE_LOG_FATAL << "previous status not matching: " << i_status;
  }

  // assign
  *l_st = i_status;
}

bool edge::parallel::Shared::getStatusAll( t_status     i_status,
                                           unsigned int i_id ) {
  // find the correct work region
  std::size_t l_rg = getWrkRgn( i_id );

  volatile WrkPkg* l_wps = &m_wrkRgns[l_rg].wrkPkgs[0];

  // iterate over all workers and check if we are finished
  for( int l_td = 0; l_td < m_nWrks; l_td++ ) {
    if(l_wps[l_td].status != i_status) return false;
  }

#ifdef PP_USE_OMP
#pragma omp flush
#endif

  return true;
}

void edge::parallel::Shared::balance() {
  m_balancing.balance();
  if( EDGE_VLOG_IS_ON(2) )
    m_balancing.print();
}