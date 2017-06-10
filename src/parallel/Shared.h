/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2015-2016, Regents of the University of California
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
#ifndef SHARED_H_
#define SHARED_H_

#include <vector>
#include "data/layout.hpp"
#include "data/SparseEntities.hpp"
#include "parallel/global.h"
#include "io/logging.h"

namespace edge {
  namespace parallel {
    class Shared;
  }
}

class edge::parallel::Shared {
  private:
    //! number of workers
    int m_nWrks;

  public:
    // per-thread status of a work package
    typedef enum {
      RDY, // ready
      IPR, // in progress
      FIN, // finished
      WAI  // waiting to be scheduled
    } t_status;

   private:
    // definition of a reoccurring a work package
    struct WrkPkg {
      //! status
      t_status status;

      //! entities covered by this work package
      t_timeRegion ents;

      //! sparse entities covered by the work package; dimension corresponds to the sparse types given at init
      std::vector< t_timeRegion> spEn;
    };

    // work region containing work packages of all threads for this region
    struct WrkRgn {
      //! id of the region
      unsigned int id;

      //! step in the control flow
      unsigned short step;

      //! time group
      int_tg tg;

      //! priority of the region
      int prio;

      //! sparse types
      std::vector< int_spType > spTypes;

      //! work packages of the threads
      std::vector< WrkPkg > wrkPkgs;
    };

    //! work regions present in the simulation, sorted by priority (descending).
    std::vector< WrkRgn > m_wrkRgns;

    /**
     * Gets the work region for the given id.
     *
     * @param i_id id for which the work region is derived.
     **/
    std::size_t getWrkRgn( unsigned int i_id );

  public:
    /**
     * Prints the shared memory config.
     **/
    void print();

    /**
     * Constructor.
     *
     * Remark: This should be called outside of the omp-parallel region.
     **/
    Shared(){};

    /**
     * Destructor.
     *
     * Remark: This should be called outside of the omp-parallel region.
     **/
     ~Shared(){};

    /**
     * Initialization of the shared memory parallelization.
     * Worker threads will be assigned to threads, 0..(nWrks-1).
     *
     * Remark: This should be called outside of the omp-parallel region.
     *
     * @param i_nWrk snumber of of worker-threads. If 0 the class decides.
     **/
    void init( unsigned int i_nWrks = 0 );

    /**
     * Determine if the thread is the lead of the communication threads
     *
     * @return true if the thread is the lead of the communication threads, false otherwise.
     **/
    bool isCommLead();

    /**
     * Determines if the calling thread is a worker.
     *
     * @return true if the thread is a worker, false otherwise.
     **/
    bool isWrk();

    /**
     * Determines if the calling thread schedules work and communication tasks.
     *
     * @return true if the thread schedules work and comm tasks, false otherwise.
     **/
    bool isSched();

    /**
     * Determines if the calling thread communicates.
     *
     * @return true if the thread communicates, false otherwise.
     **/
    bool isComm();

    /**
     * Registers a work region in the shared memory parallelization.
     *
     * The number of entities and first entity (sparse counting) in the work region for
     * all given sparse types are stored for all threads.
     *
     * @param i_tg time group.
     * @param i_step step in the control flow.
     * @param i_id id of the region.
     * @param i_first first dense-entity covered.
     * @param i_size number of dense-entities.
     * @param i_prio priority of the work region.
     * @param i_nSpTypes number of sparse types.
     * @param i_spType bitmasks considered for sparse-entities.
     * @param i_enChars entity characteristics.
     **/
    template <typename T = t_vertexChars>
    void regWrkRgn( int_tg          i_tg,
                    unsigned short  i_step,
                    unsigned int    i_id,
                    int_el          i_first,
                    int_el          i_size,
                    int             i_prio=0,
                    unsigned short  i_nSpTypes=0,
                    int_spType     *i_spType=NULL,
                    const T        *i_enChars=NULL ) {
    // wait for all threads since memory modifications before might result in inconsitent results
#ifdef PP_USE_OMP
#pragma omp barrier
#endif
    if( g_thread == 0 ) {
      // create a local work region
      WrkRgn l_wrkRgn;
      l_wrkRgn.wrkPkgs.resize( m_nWrks );
      l_wrkRgn.tg   = i_tg;
      l_wrkRgn.step = i_step;
      l_wrkRgn.id   = i_id;
      l_wrkRgn.prio = i_prio;

      // derive shared size per worker
      int_el l_size = i_size / m_nWrks;

      // distribute the equal contribution
      for( int l_td = 0; l_td < m_nWrks; l_td++ ) {
        l_wrkRgn.wrkPkgs[l_td].status    = WAI;
        l_wrkRgn.wrkPkgs[l_td].ents.size = l_size;
      }

      // derive remainder which doesn't match the number of workers
      l_size = i_size % m_nWrks;

      // distribute remainder round-robin
      int_el l_tdRr = m_wrkRgns.size() % m_nWrks;
      for( int_el l_rm = l_size; l_rm > 0; l_rm-- ) {
        l_wrkRgn.wrkPkgs[l_tdRr].ents.size++;
        l_tdRr++;
        l_tdRr = l_tdRr % m_nWrks;
      }

      // set first positions
      l_wrkRgn.wrkPkgs[0].ents.first = i_first;
      for( int l_td = 1; l_td < m_nWrks; l_td++ ) {
        l_wrkRgn.wrkPkgs[l_td].ents.first = l_wrkRgn.wrkPkgs[l_td-1].ents.first+
                                            l_wrkRgn.wrkPkgs[l_td-1].ents.size;
      }

      // determine the first sparse entry in the step for every type
      for( unsigned short l_st = 0; l_st < i_nSpTypes; l_st++ ) {
        l_wrkRgn.spTypes.push_back( i_spType[l_st] );
      }

      std::vector< int_el > l_wpSizes(m_nWrks);
      for( int l_td = 0; l_td < m_nWrks; l_td++ ) {
        l_wpSizes[l_td] = l_wrkRgn.wrkPkgs[l_td].ents.size;
        l_wrkRgn.wrkPkgs[l_td].spEn.resize(i_nSpTypes);
      }

      std::vector< int_el > l_wpFirstSp(m_nWrks);
      std::vector< int_el > l_wpSizesSp(m_nWrks);
      for( unsigned short l_st = 0; l_st < i_nSpTypes; l_st++ ) {
        data::SparseEntities::subRgnsSpId( i_first,
                                           m_nWrks,
                                           l_wpSizes.data(),
                                           i_spType[l_st],
                                           i_enChars,
                                           l_wpFirstSp.data(),
                                           l_wpSizesSp.data() );

        for( int l_td = 0; l_td < m_nWrks; l_td++ ) {
          l_wrkRgn.wrkPkgs[l_td].spEn[l_st].first = l_wpFirstSp[l_td];
          l_wrkRgn.wrkPkgs[l_td].spEn[l_st].size  = l_wpSizesSp[l_td];
        }
      }

      // put the work region in the right spot based on priority
      std::size_t l_ps = m_wrkRgns.size();

      for( std::size_t l_rg = 0; l_rg < m_wrkRgns.size(); l_rg++ ) {
        if( m_wrkRgns[l_rg].prio < l_wrkRgn.prio ) {
          l_ps = l_rg;
          break;
        }
      }

      m_wrkRgns.insert( m_wrkRgns.begin()+l_ps, l_wrkRgn );
    }

// sync memory view
#ifdef PP_USE_OMP
#pragma omp barrier
#endif
    }

    /**
     * Gets work for the calling thread.
     *
     * @param o_tg time group.
     * @param o_step step in the computational scheme.
     * @param o_id id of the work region; numeric_limits<unsigned int>::max() if none available.
     * @param o_first first entity to work on; numeric_limits<int_el>::max() if none available.
     * @param o_size number of entities to work on; numeric_limits<int_el>::max() if none available.
     * @param o_spEn sparse-entities, one range for every defined sparse type.
     * @return true if work is available; false otherwise.
     **/
    bool getWrkTd( int_tg          &o_tg,
                   unsigned short  &o_step,
                   unsigned int    &o_id,
                   int_el          &o_first,
                   int_el          &o_size,
                   t_timeRegion   **o_spEn=nullptr );

    /**
     * Sets the given status of the work package in the respective region for the calling thread.
     *
     * @param i_st status which is set.
     * @param i_id id of the work region.
     **/
    void setStatusTd( t_status     i_status,
                      unsigned int i_id );

    /**
     * Checks if the status of all workers matches for the region.
     *
     * @param i_status status to check.
     * @param i_id id of the region.
     * @return true if all workers' statuses matches, false otherwise.
     **/
    bool getStatusAll( t_status     i_status,
                       unsigned int i_id );

    /**
     * Sets the status of the region for all workers.
     *
     * @param i_status status to set.
     * @param i_id id of the region.
     **/
    void setStatusAll( t_status     i_status,
                       unsigned int i_id );

    /**
     * Resets the status for all regions for all workers to the given value.
     *
     * @param i_status status to reset to.
     **/
    void resetStatus( t_status i_status );
};

#endif
