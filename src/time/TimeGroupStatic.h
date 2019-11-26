/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
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
 * LTS cluster of EDGE with static, DOF-independent time steps.
 **/

#ifndef EDGE_TIME_TIME_GROUP_STATIC_H
#define EDGE_TIME_TIME_GROUP_STATIC_H

#include "data/Internal.hpp"
#include "data/EntityLayout.type"
#include "io/Receivers.h"

namespace edge {
  namespace time {
    class TimeGroupStatic;
  }
}

class edge::time::TimeGroupStatic {
  private:
    //! multiplier of the fundamental time step which defines the time step of this group
    unsigned short m_funMul;

    //! divisor of the largest LTS group's time step which results in this group's time step
    unsigned short m_maxDiv;

    //! number of performend time steps since last synchronization
    volatile std::size_t m_nTsSync;

    //! total number of performed time steps
    volatile std::size_t m_nTsPer;

    //! number of required full time steps before synchronization phase
    std::size_t m_nTsReqFull;

    //! number of performed time predictions (inner, send) since last sync
    std::size_t m_nTimePredSync[2];

    //! number of performed DOF updates (inner, send) since last sync
    std::size_t m_nDofUpSync[2];

    //! full time step of the group (without synchronization adjustments)
    double m_dtFull;

    //! time step used during the synchronization phase
    double m_dtSync;

    //! time step of the current update
    volatile double m_dt;

    //! covered simulation time
    double m_covSimTime;

    //! internal state
    data::Internal & m_internal;

    //! send pointers
    unsigned char ** m_sendPtrs = nullptr;

    //! receive pointers
    unsigned char ** m_recvPtrs = nullptr;

    /**
     * Sets the time step for the current update.
     **/
    void setDt();

  public:
    /**
     * Constructor.
     *
     * @param i_nTgs number of time groups.
     * @param i_tgId id of this time group.
     * @param i_internal internal data.
     * @param i_sendPtrs send pointers.
     * @param i_recvPtrs receive pointers.
     **/
    TimeGroupStatic( unsigned short   i_nTgs,
                     unsigned short   i_tgId,
                     data::Internal & i_internal,
                     unsigned char ** i_sendPtrs,
                     unsigned char ** i_recvPtrs );

    /**
     * Sets up the cluster for iterations until the given synchronization point.
     *
     * @param i_dtFun fundamental time step.
     * @param i_time time to advance forward in time.
     **/
    void setUp( double i_dtFun,
                double i_time );

    /**
     * Increases the time prediction counter for the inner elements.
     **/
    void updateTimePredInner(){ m_nTimePredSync[0]++; };

    /**
     * Increases the time prediction counter for the send elements.
     **/
    void updateTimePredSend(){ m_nTimePredSync[1]++; };

    /**
     * Increases the DOF update counter for the innner elements.
     **/
    void updateDofUpInner(){ m_nDofUpSync[0]++; };

    /**
     * Increases the DOF update counter for the send elements.
     **/
    void updateDofUpSend(){ m_nDofUpSync[1]++; };

    /**
     * Gets the number of time predicitons since sync for inner elements.
     *
     * @return number of time predictions.
     **/
    std::size_t nTimePredInner(){ return m_nTimePredSync[0]; }

    /**
     * Gets the number of time predicitons since sync for send elements.
     *
     * @return number of time predictions.
     **/
    std::size_t nTimePredSend(){ return m_nTimePredSync[1]; }

    /**
     * Gets the number of DOF updates since sync for inner elements.
     *
     * @return number of DOF updates.
     **/
    std::size_t nDofUpInner(){ return m_nDofUpSync[0]; }

    /**
     * Gets the number of DOF updates since sync for send elements.
     *
     * @return number of DOF updates.
     **/
    std::size_t nDofUpSend(){ return m_nDofUpSync[1]; }

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
     **/
    void computeStep( unsigned short         i_step,
                      int_el                 i_first,
                      int_el                 i_size,
                      int_el         const * i_enSp,
                      io::Receivers        & io_recvs );

    /**
     * Gets the number of updates the time group performed since the last synchronization.
     *
     * @return number of updates since last synchronization
     **/
    std::size_t getUpdatesSync() { return m_nTsSync; };

    /**
     * Gets the number of updates the time group performed.
     *
     * @return number of performed updates.
     **/
    std::size_t getUpdatesPer() { return m_nTsPer; }

    /**
     * Gets the covered simulation time.
     *
     * @return covered simulation time.
     **/
    double getCovSimTime() { return m_covSimTime; };

    /**
     * Checks if the time group is performing its last time step
     *
     * @return true if the group is performing its last time step. false otherwise
     **/
    bool lastTimeStep() const {
      return (m_nTsSync == m_nTsReqFull + m_maxDiv - 1 );
    }

    /**
     * Checks if the time group reached the synchronization time.
     *
     * @return true if the group reached the synchronization time. false otherwise.
     **/
    bool finished() const {
      return (m_nTsSync == m_nTsReqFull + m_maxDiv );
    }
};

#endif
