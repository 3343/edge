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
 * LTS cluster of EDGE with static, data-independent time step characterisitcs.
 **/

#ifndef EDGE_TIME_TIME_GROUP_STATIC_H
#define EDGE_TIME_TIME_GROUP_STATIC_H

#include "data/Internal.hpp"
#include "data/EntityLayout.type"
#include "io/Receivers.h"
#include "io/ReceiversSf.hpp"

namespace edge {
  namespace time {
    class TimeGroupStatic;
  }
}

class edge::time::TimeGroupStatic {
  //private:
    //! rate of the cluster
//    const int_ts m_rate;

    //! fundamental time step multiple of the cluster
    const int_ts m_funMult;

    //! number of performend updates since last synchronization
    volatile int_ts m_updatesSync;

    //! total number of performed updates
    volatile int_ts m_updatesPer;

    //! number of required updates until synchronization
    int_ts m_updatesReq;

    //! general time step for all but last update
    double m_dTgen;

    //! final time step
    double m_dTfin;

    //! time step of the current update
    volatile double m_dT;

    //! covered simulation time
    double m_covSimTime;

    //! elements/faces under control of this cluster
    data::Internal &m_internal;

    /**
     * Sets the time step for the current update.
     **/
    void setDt();

  public:
    /**
     * Constructor.
     *
     * @param i_rate local rate of this cluster.
     * @param i_funMult global rate with respect to the fundamental time step.
     * @param i_internal data under control of this cluster.
     **/
    TimeGroupStatic( int_ts          i_rate,
                     int_ts          i_funMult,
                     data::Internal &i_internal );

    /**
     * Sets up the cluster for iterations until the given synchronization point.
     *
     * @param i_dTfun fundamental time step.
     * @param i_time time to advance forward in time.
     **/
    void setUp( double i_dTfun,
                double i_time );

    /**
     * Updates the time step info.
     **/
    void updateTsInfo();

    /**
     * Computes a step of the cluster.
     *
     * @param i_step id of the step which is computed.
     * @param i_first first entity to which this applies.
     * @param i_size number of entities involved in this step and for the calling worker thread.
     * @param i_enSp sparse entities.
     * @param i_recvs receiver output (sensitive to high output frequencies).
     * @param io_recvsSf receivers at sub-faces.
     **/

    void computeStep( unsigned short                              i_step,
                      int_el                                      i_first,
                      int_el                                      i_size,
                      int_el                              const * i_enSp,
                      io::Receivers                             & io_recvs,
                      io::ReceiversSf< real_base,
                                       T_SDISC.ELEMENT,
                                       ORDER,
                                       N_CRUNS >                & io_recvsSf );

    /**
     * Prepare the limiter for synchronization by resetting the ids for admissibility and extrema.
     **/
    void limSync();

    /**
     * Gets the number of updates the time group performed since the last synchronization.
     *
     * @return number of updates since last synchronization
     **/
    int_ts getUpdatesSync() { return m_updatesSync; };

    /**
     * Gets the number of updates the time group performed.
     *
     * @return number of performed updates.
     **/
    int_ts getUpdatesPer() { return m_updatesPer; }

    /**
     * Gets the covered simulation time.
     *
     * @return covered simulation time.
     **/
    double getCovSimTime() { return m_covSimTime; };

    /**
     * Checks if the cluster is performing its last time step
     *
     * @return true if the cluster is performing its last time step. false otherwise
     **/
    bool lastTimeStep() const {
      return m_updatesReq == 1;
    }

    /**
     * Checks if the cluster reached the synchronization time.
     *
     * @return true if the cluster reached the synchronization time. false otherwise.
     **/
    bool finished() const {
      return m_updatesReq == 0;
    }

    /**
     * Gets the unique identifier for the data of the admissibility.
     *
     * @param i_adm local id of the admissibility (previous, canditate, limited #1, limited #2).
     **/
    std::uintptr_t getAdmDataId( unsigned short i_adm ) {
      EDGE_CHECK_LT( i_adm, 4 );
      return reinterpret_cast< std::intptr_t >( m_internal.m_globalShared2[0].adm[i_adm] );
    }

    /**
     * Gets the unique identifier for the data of the SC Dofs.
     *
     * @param i_do local id of the SC DOFs (previous/current).
     **/
    std::uintptr_t getDofsScDataId( unsigned short i_do ) {
      EDGE_CHECK_LT( i_do, 2 );
      return reinterpret_cast< std::intptr_t >( m_internal.m_globalShared2[0].tDofs[i_do] );
    }

    /**
     * Gets the unique identifier for the data of the extrema.
     *
     * @param i_ex local id of the extrema (previous/current).
     **/
    std::uintptr_t getExDataId( unsigned short i_ex ) {
      EDGE_CHECK_LT( i_ex, 2 );
      return reinterpret_cast< std::intptr_t >( m_internal.m_globalShared2[0].ext[i_ex] );
    }
};

#endif
