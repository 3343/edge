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
 * Management of the time stepping.
 **/

#ifndef EDGE_TIME_MANAGER_H
#define EDGE_TIME_MANAGER_H

#include "parallel/Shared.h"
#include "constants.hpp"
#include "io/Receivers.h"
#include "TimeGroupStatic.h"
#include <vector>

namespace edge {
  namespace time {
    class Manager;
  }
}

class edge::time::Manager {
  private:
    //! fundamental time step
    const double m_dTfun;

    //! shared memory parallelization
    parallel::Shared & m_shared;

    //! mpi parallelization
    parallel::Mpi & m_mpi;

    //! receiver output
    io::Receivers & m_recvs;

    //! clusters under control of the time manager
    std::vector< TimeGroupStatic* > m_timeGroups;

    //! control flow of the scheme
    unsigned short (*m_cflow)[N_ENTRIES_CONTROL_FLOW];

    //! true if the manager reached the desired synchronization point
    volatile bool m_finished = false;

    /**
     * Returns true if the time predictions of the neighboring smaller and large time group (on this rank) are available for an update.
     *
     * @param i_tg time group for which the information is requested.
     * @return true if available, false if not.
     **/
    bool getTimePredAvailable( unsigned short i_tg );

    /**
     * Returns true if this cluster's time predicitions have been consumed (on this rank) by neighboring clusters.
     *
     * @param i_tg time group for which the information is requested.
     * @return true if consumed, false if not.
     **/
    bool getTimePredConsumed( unsigned short i_tg );

    /**
     * Runs scheduling tasks.
     **/
    void schedule();

    /**
     * Runs communication tasks.
     **/
    void communicate();

    /**
     * Performs computations.
     **/
    void compute();

  public:
    /**
     * Constructor of the time stepping management.
     *
     * @param i_dt fundamental time step.
     * @param i_shared shared memory parallelization.
     * @param i_mpi mpi parallelization.
     * @param i_timeGroups time groups
     * @param i_recvs modal receivers.
     **/
    Manager( double                               i_dt,
             parallel::Shared                   & i_shared,
             parallel::Mpi                      & i_mpi,
             std::vector< TimeGroupStatic >     & i_timeGroups,
             io::Receivers                      & i_recvs );

    /**
     * Destructor.
     **/
    ~Manager();

    /**
     * Advances in time for the given time.
     *
     * @param i_time time to advance forward in time.
     **/
    void simulate( double i_time );
};

#endif
