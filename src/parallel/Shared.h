/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (alex.breuer AT uni-jena.de)
 *
 * @section LICENSE
 * Copyright (c) 2021, Friedrich Schiller University Jena
 * Copyright (c) 2019, Alexander Breuer
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
#ifndef EDGE_PARALLEL_SHARED_H
#define EDGE_PARALLEL_SHARED_H

#include <cstdint>
#include <vector>
#include "data/SparseEntities.hpp"
#include "data/EntityLayout.type"
#include "parallel/global.h"
#include "LoadBalancing.h"
#include "io/logging.h"

namespace edge {
  namespace parallel {
    class Shared;
  }
}

class edge::parallel::Shared {
  public:
    //! number of workers
    int m_nWrks;

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

      //! 64byte padding for separate signaling cache lines
      uint64_t padding[8];
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

    //! dynamic load balancing
    LoadBalancing m_balancing;

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
     * If separateWrks is set, thread 0 will be used for exclusively scheduling and MPI.
     *
     * Remark: This should be called outside of the omp-parallel region.
     *
     * @param i_separateWrks if true, the workers will be isolated from the scheduling/comm thread
     **/
    void init( bool i_separateWrks = true );

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
    void regWrkRgn( unsigned short   i_tg,
                    unsigned short   i_step,
                    unsigned int     i_id,
                    std::size_t      i_first,
                    std::size_t      i_size,
                    int              i_prio=0,
                    unsigned short   i_nSpTypes=0,
                    int_spType     * i_spType = nullptr,
                    const T        * i_enChars = nullptr ) {
    // wait for all threads since memory modifications before might result in inconsistent results
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

      // init to wait
      for( int l_td = 0; l_td < m_nWrks; l_td++ ) {
        l_wrkRgn.wrkPkgs[l_td].status = WAI;
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

      // register the work region in the dynamic load balancing
      m_balancing.regWrkRgn( l_ps,
                             i_first,
                             i_size,
                             i_nSpTypes,
                             i_spType,
                             i_enChars );
    }

// sync memory view
#ifdef PP_USE_OMP
#pragma omp barrier
#endif
    }

    /**
     * Gets work for the calling thread.
     * If work is available, OMP-flush is called.
     *
     * @param o_tg time group.
     * @param o_step step in the computational scheme.
     * @param o_id id of the work region; numeric_limits<unsigned int>::max() if none available.
     * @param o_first first entity to work on; numeric_limits<int_el>::max() if none available.
     * @param o_size number of entities to work on; numeric_limits<int_el>::max() if none available.
     * @param o_spEn sparse-entities, one range for every defined sparse type.
     *
     * @return true if work is available; false otherwise.
     **/
    bool getWrkTd( unsigned short & o_tg,
                   unsigned short & o_step,
                   unsigned int   & o_id,
                   int_el         & o_first,
                   int_el         & o_size,
                   int_el         * o_spEn );

    /**
     * Sets the given status of the work package in the respective region for the calling thread.
     * Additionally, OMP-flush is called.
     *
     * @param i_st status which is set.
     * @param i_id id of the work region.
     **/
    void setStatusTd( t_status     i_status,
                      unsigned int i_id );

    /**
     * Checks if the status of all workers matches for the region.
     * If this is true, OMP-flush is called.
     *
     * @param i_status status to check.
     * @param i_id id of the region.
     * @return true if all workers' statuses matches, false otherwise.
     **/
    bool getStatusAll( t_status     i_status,
                       unsigned int i_id );

    /**
     * Sets the status of the region for all workers.
     * Additionally, OMP-flush is called.
     *
     * @param i_status status to set.
     * @param i_id id of the region.
     **/
    void setStatusAll( t_status     i_status,
                       unsigned int i_id );

    /**
     * Resets the status for all regions for all workers to the given value.
     * Additionally, OMP-flush is called.
     *
     * @param i_status status to reset to.
     **/
    void resetStatus( t_status i_status );

    /**
     * @brief Balances the work packages.
     */
    void balance();

    /**
     * @brief Performs NUMA-aware zero-initialization of the given array through first-touch.
     *        Should be called within an OpenMP-parallel region from all threads.
     *        The inits of the workers are done as OpenMP-critical (one by one).
     *        After all init an OpenMP-barrier is called.
     *
     * @param i_nWrks number of workers.
     * @param i_nEns number of entries in the given array, which will be split equally among the available workers. Remainders will be added to the first workers.
     * @param o_arr array, which will be initialized.
     *
     * @paramt TL_T_EN type of the array entries.
     */
    template< typename TL_T_EN >
    static void numaInit( unsigned int        i_nWrks,
                          std::size_t const   i_nEns,
                          TL_T_EN           * o_arr ) {
      // return early, if this is not a worker.
      if( g_worker < 0 ) {
        // wait for other threads
#ifdef PP_USE_OMP
#pragma omp barrier
#endif

        return;
      }

      // split among the workers
      std::size_t l_split = i_nEns / std::size_t(i_nWrks);
      std::size_t l_rem   = i_nEns % std::size_t(i_nWrks);

      // derive start position of the calling worker in the array
      TL_T_EN * l_arr = o_arr;
      l_arr += l_split * (std::size_t) g_worker;
      if( l_rem > 0 ) l_arr += std::min( l_rem, (std::size_t) g_worker );

      // derive number of entries under control of the worker
      std::size_t l_nEns = l_split;
      if( (std::size_t) g_worker < l_rem ) l_nEns++;

      // perform NUMA-aware init
#ifdef PP_USE_OMP
#pragma omp critical
#endif
      for( std::size_t l_en = 0; l_en < l_nEns; l_en++ ) l_arr[l_en] = 0;

      // wait for other threads
#ifdef PP_USE_OMP
#pragma omp barrier
#endif
    }

    /**
     * @brief Performs NUMA-aware zero-initialization of the given array through first-touch.
     *        Should be called within an OpenMP-parallel region from all threads.
     *        For every region, the inits of the workers are done as OpenMP-critical (one by one).
     *        After all inits in a region, an OpenMP-barrier is called.
     *
     * @param i_nTgs number of time groups.
     * @param i_nTgEnsIn number of inner entities per time group.
     * @param i_nTgEnsSe number of send entities per time group.
     * @param i_nWrks number of workers.
     * @param i_nVasPerEn number of values per entity
     * @param o_arr array, which will be initialized.
     *
     * @paramt TL_T_VA type of the values.
     */
    template< typename TL_T_LID,
              typename TL_T_VA >
    static void numaInit( unsigned short         i_nTgs,
                          TL_T_LID       const * i_nTgEnsIn,
                          TL_T_LID       const * i_nTgEnsSe,
                          unsigned int           i_nWrks,
                          std::size_t    const   i_nVasPerEn,
                          TL_T_VA              * o_arr ) {
      // offset
      std::size_t l_off = 0;

      // inner
      for( unsigned short l_tg = 0; l_tg < i_nTgs; l_tg++ ) {
        std::size_t l_nVas = i_nTgEnsIn[l_tg];
        l_nVas *= i_nVasPerEn;

        edge::parallel::Shared::numaInit( i_nWrks,
                                          l_nVas,
                                          o_arr+l_off );
        l_off += l_nVas;
      }

      // send
      for( unsigned short l_tg = 0; l_tg < i_nTgs; l_tg++ ) {
        std::size_t l_nVas = i_nTgEnsSe[l_tg];
        l_nVas *= i_nVasPerEn;

        edge::parallel::Shared::numaInit( i_nWrks,
                                          l_nVas,
                                          o_arr+l_off );
        l_off += l_nVas;
      }
    }
};

#endif
