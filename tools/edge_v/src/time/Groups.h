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
#ifndef EDGE_V_TIME_GROUPS_H
#define EDGE_V_TIME_GROUPS_H

#include <cstddef>
#include "../constants.h"

namespace edge_v {
  namespace time {
    class Groups;
  }
}

/**
 * Derivation of time step groups.
 * Regularizes local time stepping.
 **/
class edge_v::time::Groups {
  private:
    //! number of of elements
    std::size_t m_nEls = 0;

    //! number of time groups
    unsigned short m_nGroups = 0;

    //! time steps of the time groups
    double * m_tsIntervals = nullptr;

    //! time groups of the elements
    unsigned short * m_elTg;

    //! number of elements in each time groups
    std::size_t * m_nGroupEls;

    //! theoretical loads of the different time stepping schemes: GTS, grouped LTS before normalization, grouped LTS after normalization, per-element LTS
    double m_loads[4] = {};

    /**
     * Gets the theoretical loads of global time stepping, grouped local time stepping and per-element time stepping.
     *
     * @param i_nEls number of elements.
     * @param i_tsCfl cfl-conditioned time steps of the elements.
     * @param i_tsGroups time steps of the groups.
     * @param i_elTg time groups of the elements.
     * @param o_gts will be set to the load per GTS time step.
     * @param o_ltsGrouped will be set to the summed load per element update, when considering the time groups.
     * @param o_ltsPerElement will be set to the summed load per element update, when assuming that every element operates its CFL-time step.
     **/
    static void getLoads( std::size_t            i_nEls,
                          double         const * i_tsCfl,
                          double         const * i_tsGroups,
                          unsigned short const * i_elTg,
                          double               & o_gts,
                          double               & o_ltsGrouped,
                          double               & o_ltsPerElement );

    /**
     * Normalizes the element time groups, by enforcing single time group jumps between face-adjacent elements.
     *
     * @param i_elTy element type.
     * @param i_nEls number of elements.
     * @param i_elFaEl elements adjacent to elements.
     * @param io_elTg time groups of the elements, which will be updated.
     *
     * @return number of normalized elements.
     **/
    static std::size_t normalizeElTgs( t_entityType           i_elTy,
                                       std::size_t            i_nEls,
                                       std::size_t    const * i_elFaEl,
                                       unsigned short       * io_elTg );

    /**
     * Sets the elements' time groups.
     *
     * @param i_nEls number of elements.
     * @param i_nGroups number of groups.
     * @param i_tsGroups time steps of the groups.
     * @param i_tsCfl CFL-conditioned time steps of the elements.
     * @param o_elTg will be set to the time groups of the elements.
     **/
    static void setElTg( std::size_t            i_nEls,
                         unsigned short         i_nGroups,
                         double         const * i_tsGroups,
                         double         const * i_tsCfl,
                         unsigned short       * o_elTg );

    /**
     * Sets the number of elements per group.
     *
     * @param i_nEls number of elements.
     * @param i_nGroups number of groups.
     * @param i_elTg time groups of the elements.
     * @param o_nGroupEls will be set to the number of elements for every time group.
     **/
    static void nGroupEls( std::size_t            i_nEls,
                           unsigned short         i_nGroups,
                           unsigned short const * i_elTg,
                           std::size_t          * o_nGroupEls );

  public:
    /**
     * Constructor.
     * The term rate refers to the ratio a time step group has to the next one.
     * For example, setting the rates of (1.7, 2.0, 1.9) will assign groups as follows:
     *
     * [1dt, ..., 1.7dt, ..., 3.4, ..., 6.46, ..., infty )
     *     group 1  |  group 2 |  group 3 |  group 4
     *
     * If no rates are given, we get global time stepping.
     *
     * @param i_elTy entity type of the elements.
     * @param i_nEls number of elements.
     * @param i_elFaEl elements adjacent to elements (faces as bridge).
     * @param i_nRates number of rates.
     * @param i_rates rates.
     * @param i_ts CFL-conditioned and normalized (minimum is 1) time steps of the elements.
     **/
    Groups( t_entityType           i_elTy,
            std::size_t            i_nEls,
            std::size_t    const * i_elFaEl,
            unsigned short         i_nRates,
            double         const * i_rates,
            double         const * i_ts );

    /**
     * Destructor.
     **/
    ~Groups();

    /**
     * Prints the stats of the grouping.
     **/
    void printStats() const;

    /**
     * Gets the number of time groups.
     *
     * @return number of time groups.
     **/
    unsigned short nGroups() const { return m_nGroups; }

    /**
     * Gets the time step intervals of the groups.
     * For example for rates (1.5, 2.0, 2.0), the following values are returned:
     *   1, 1.5, 3.0, 6.0, std::numeric_limits< double >::max()
     *
     * @return time step intervals.
     **/
    double const * getTsIntervals() const { return m_tsIntervals; }

    /**
     * Gets the elements' time groups.
     *
     * @return elements' time groups.
     **/
    unsigned short const * getElTg() const { return m_elTg; }
};

#endif