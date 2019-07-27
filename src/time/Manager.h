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
 * Management of the time stepping.
 **/

#ifndef EDGE_TIME_MANAGER_H
#define EDGE_TIME_MANAGER_H

#include "parallel/Shared.h"
#include "constants.hpp"
#include "io/Receivers.h"
#include "io/ReceiversSf.hpp"
#include "TimeGroupStatic.h"
#include <vector>

namespace edge {
  namespace time {
    class Manager;
  }
}

class edge::time::Manager {
  //private:
    //! fundamental time step
    const double m_dTfun;

    //! shared memory parallelization
    parallel::Shared &m_shared;

    //! mpi parallelization
    parallel::Mpi &m_mpi;

    //! receiver output
    io::Receivers &m_recvs;

    //! receiver output at sub-faces
    io::ReceiversSf<real_base, T_SDISC.ELEMENT, ORDER, N_CRUNS> &m_recvsSf;

    //! clusters under control of the time manager
    std::vector< TimeGroupStatic* > m_timeGroups;

    //! control flow of the scheme
    unsigned short m_cflow[N_ENTRIES_CONTROL_FLOW];

    //! true if the manager reached the desired synchronization point
    volatile bool m_finished = false;

    //! scheduling loop
    void schedule();

    //! communication loop
    void communicate();

    //! computational loop
    void compute();

  public:
    /**
     * Constructor of time step management.
     *
     * @param i_dT fundamental time step.
     * @param i_shared shared memory parallelization.
     * @param i_mpi mpi parallelization.
     * @param i_recvs modal receivers.
     * @param i_recvsSf receivers at sub-faces.
     **/
    Manager(       double                              i_dT,
                   parallel::Shared                    &i_shared,
                   parallel::Mpi                       &i_mpi,
                   io::Receivers                       &i_recvs,
                   io::ReceiversSf< real_base,
                                    T_SDISC.ELEMENT,
                                    ORDER,
                                    N_CRUNS >        &i_recvsSf ):
     m_dTfun(i_dT), m_shared(i_shared), m_mpi(i_mpi), m_recvs(i_recvs), m_recvsSf(i_recvsSf){};

    /**
     * Adds a time group to the time manager.
     *
     * @param i_timeGroup time group to add.
     **/
    void add( TimeGroupStatic *i_timeGroup );

    /**
     * Advances in time for the given time.
     *
     * @param i_time time to advance forward in time.
     **/
    void simulate( double i_time );
};

#endif
