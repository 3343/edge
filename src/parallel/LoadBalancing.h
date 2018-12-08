/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2018, Regents of the University of California
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
#ifndef EDGE_PARALLEL_LOADBALANCING_H
#define EDGE_PARALLEL_LOADBALANCING_H

#include <vector>
#include "constants.hpp"
#include "data/SparseEntities.hpp"
#include "monitor/Timer.hpp"
#include "io/logging.h"

namespace edge {
  namespace parallel {
    class LoadBalancing;
  }
}

class edge::parallel::LoadBalancing {
  private:
    //! time tolerance, if the measured time falls below it for a thread and work package, the load is equally distributed
    const double m_zeroTime;

    //! maximum allowed imbalance
    const double m_maxImbalance;

    //! number of performed balancing steps
    std::size_t m_nBalanced;

    //! number of workers
    unsigned int m_nWrks;

    //! load balancing definition of a work package
    struct WrkPkg {
      //! timer for the work package
      edge::monitor::Timer timer;

      //! first entity of the work package
      std::size_t first;

      //! size of the work package
      std::size_t size;

      //! first sparse entities of the work package
      std::vector< size_t > firstSp;
    };

    //! load balancing definition of a work region
    struct WrkRgn {
      //! first entity
      std::size_t first;

      //! number of entities
      std::size_t size;

      //! work packages of the region
      std::vector< WrkPkg > wrkPkgs;

      //! first sparse entity id
      std::vector< std::size_t > firstSp;

      //! sparse-dense link for the given work region
      std::vector< std::vector < std::size_t > > spDe;

      //! minimum elapsed time of previous balancing (not the current one)
      double elaMin;

      //! summed elapsed time of previous balancing (not the current one)
      double elaSum;

      //! maximum elapsed time of previous balancing (not the current one)
      double elaMax;
    };

    //! work regions present in the simulation
    std::vector< WrkRgn > m_wrkRgns;

    /**
     * @brief Resolves the sparse entities, based on a computed balancing.
     * 
     * @param i_id work region for which the sparse entities are resolved.
     */
    void resolveSpEn( unsigned short i_id );

    /**
     * @brief (Re-)Balances the work region with the given id into work packages.
     *        The balancing is based on internal performance time monitoring of the workers.
     *        If any of the internal clocks is non-positive, an equal work distribution is returned.
     *        This is the case, if called for the first time for a work region.
     *        This steps also resets the internal clocks of the work packages to zero.
     *
     * @param i_id id of the work region, which gets balanced.
     */
    void balanceWrkRgn( unsigned short i_id );

  public:
    /**
     * @brief Constructor.
     *        If ts is the time of the slowest worker and that of the fastest, the imbalance is defined as:
     *        (ts - tf) / ave
     *
     * @param i_zeroTime if any of the workers falls below the zero time in a work region, work is equally distributed.
     * @param i_maxImbalance maximum allowed imbalance for a work region.
     */
    LoadBalancing( double i_zeroTime=1E-3,
                   double i_maxImbalance=2.5E-2 ): m_zeroTime(     i_zeroTime     ),
                                                   m_maxImbalance( i_maxImbalance ),
                                                   m_nBalanced(0){}; 

    /**
     * @brief Initializes the dynamic load balancing.
     * 
     * @param i_nWrks number of workers.
     */
    void init( unsigned int i_nWrks );

    /**
     * @brief (Re-)Balances all work regions.
     *        The balancing is based on internal performance time monitoring of the workers.
     *        If any of the internal clocks is non-positive for a work region, an equal work distribution is returned for the work region.
     *        This is the case, if called for the first time for a work region.
     *        This steps also resets the internal clocks of the work packages to zero.
     */
    void balance();

    /**
     * @brief Registers a work region in the load balancing.
     *
     * @param i_pos position of the work region.
     * @param i_first first dense-entity covered.
     * @param i_size number of dense-entities.
     * @param i_nSpTypes number of sparse types.
     * @param i_spType bitmasks considered for the sparse-entities.
     * @param i_enChars entity characteristics.
     *
     * @paramt TL_T_SP sparse type.
     * @paramt TL_T_EN_CHARS type of the entity characteristics.
     **/
    template< typename TL_T_SP = int,
              typename TL_T_EN_CHARS = t_vertexChars >
    void regWrkRgn( unsigned short        i_pos,
                    std::size_t           i_first,
                    std::size_t           i_size,
                    unsigned short        i_nSpTypes=0,
                    TL_T_SP             * i_spType=NULL,
                    const TL_T_EN_CHARS * i_enChars=NULL ) {
      // create new work region
      WrkRgn l_wrkRgn;
      l_wrkRgn.wrkPkgs.resize( m_nWrks );

      // init
      for( unsigned short l_wo = 0; l_wo < m_nWrks; l_wo++ ) {
        l_wrkRgn.wrkPkgs[l_wo].size  = std::numeric_limits< std::size_t >::max();
        l_wrkRgn.wrkPkgs[l_wo].first = std::numeric_limits< std::size_t >::max();
        l_wrkRgn.wrkPkgs[l_wo].firstSp.resize( i_nSpTypes );
        l_wrkRgn.wrkPkgs[l_wo].elaMin = 0;
        l_wrkRgn.wrkPkgs[l_wo].elaSum = 0;
        l_wrkRgn.wrkPkgs[l_wo].elaMax = 0;
      }

      // set boundaries of the region
      l_wrkRgn.first = i_first;
      l_wrkRgn.size = i_size;

      // derive generic sparse-dense link
      l_wrkRgn.spDe.resize( i_nSpTypes );
      for( unsigned short l_ty = 0; l_ty < i_nSpTypes; l_ty++ ) {
        // first sparse entity of the region
        std::size_t l_firstSp = data::SparseEntities::nSp( i_first,
                                                           i_spType[l_ty],
                                                           i_enChars );
        l_wrkRgn.firstSp.push_back( l_firstSp );

        // number of sparse entities in the work region
        std::size_t l_nSp = data::SparseEntities::nSp( i_size,
                                                       i_spType[l_ty],
                                                       i_enChars+i_first );

        l_wrkRgn.spDe[l_ty].resize( l_nSp );
        // sparse-dense link (relative to the region)
        data::SparseEntities::linkSpDe( i_size,
                                        i_spType[l_ty],
                                        i_enChars+i_first,
                                        l_wrkRgn.spDe[l_ty].data() );

        // add offset
        for( std::size_t l_sp = 0; l_sp < l_nSp; l_sp++ ) l_wrkRgn.spDe[l_ty][l_sp] += i_first;
      }

      // insert the work region
      EDGE_CHECK_LE( i_pos, m_wrkRgns.size() );
      m_wrkRgns.insert( m_wrkRgns.begin()+i_pos, l_wrkRgn );

      // perform initial balancing
      balanceWrkRgn( i_pos );
    }

    /**
     * @brief Gets work for the given worker in the specified region.
     * 
     * @param i_wrkRgn id of the work region.
     * @param i_worker id of the worker.
     * @param o_first will be set to first entity, covered by the worker.
     * @param o_size will be set to number of entities, covered by the worker.
     * @param o_firstSp will be set to first sparse entities, covered by the worker (if any).
     *
     * @paramt TL_T_LID type of the local ids.
     */
    template< typename TL_T_LID >
    void getWrkTd( unsigned short   i_wrkRgn,
                   unsigned short   i_worker,
                   TL_T_LID       & o_first,
                   TL_T_LID       & o_size,
                   TL_T_LID       * o_firstSp ) {
      o_first = m_wrkRgns[i_wrkRgn].wrkPkgs[i_worker].first;
      o_size  = m_wrkRgns[i_wrkRgn].wrkPkgs[i_worker].size;

      unsigned short l_nSpTy = m_wrkRgns[i_wrkRgn].wrkPkgs[i_worker].firstSp.size();
      for( unsigned short l_ty = 0; l_ty < l_nSpTy; l_ty++ ) {
        o_firstSp[l_ty] = m_wrkRgns[i_wrkRgn].wrkPkgs[i_worker].firstSp[l_ty];
      }
    }

    /**
     * @brief Starts the time monitoring for the given work package.
     * 
     * @param i_wrkRgn id of the work region
     * @param i_worker id of the worker. 
     */
    void startClock( unsigned short i_wrkRgn,
                     unsigned short i_worker );

    /**
     * @brief Stops the clock for the given work package.
     *        The elapsed time between start and stop is used for the load balancing.
     * 
     * @param i_wrkRgn id of the work region
     * @param i_worker id of the worker. 
     */
    void stopClock( unsigned short i_wrkRgn,
                    unsigned short i_worker );

    /**
     * @brief Prints summarized information on the performed load balancing.
     */
    void print();
};

#endif