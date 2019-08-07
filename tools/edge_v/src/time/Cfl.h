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
 * Time steps, based on stability constraints.
 **/
#ifndef EDGE_V_TIME_CFL_H
#define EDGE_V_TIME_CFL_H

#include "models/Model.hpp"
#include "constants.h"

namespace edge_v {
  namespace time {
    class Cfl;
  }
}

/**
 * Time step related data and function.
 * Everything is stability-related (CFL-condition) and element-local.
 * No grouping, normalization etc. is performed.
 **/
class edge_v::time::Cfl {
  private:
    //! number of elements
    std::size_t m_nEls = 0;

    //! minimum absolute time step
    double m_tsAbsMin = 0;

    //! normalized time steps of the elements
    double * m_ts = nullptr;

    /**
     * Sets the time steps of the elements.
     * Time step means diameter divided by maximum wave speed in every element.
     * Thus, any constant scaling coefficients are ignored.
     *
     * @param i_elTy element type.
     * @param i_nVes number of vertices.
     * @param i_nEls number of elements.
     * @param i_elVe vertices adjacent to the elements.
     * @param i_veCrds vertex coordinates.
     * @param i_inDia length, incircle or insphere diameters of the elements.
     * @param io_velMod velocity model.
     * @param o_tsAbsMin will be set to absolute minimum time step.
     * @param o_ts will be set to the normalized time steps (divided by absolute minimum) of the elements.
     **/
    static void setTimeSteps( t_entityType  const    i_elTy,
                              std::size_t            i_nVes,
                              std::size_t            i_nEls,
                              std::size_t   const  * i_elVe,
                              double        const (* i_veCrds)[3],
                              double        const  * i_inDia,
                              models::Model        & io_velMod,
                              double               & o_tsAbsMin,
                              double               * o_ts );

  public:
    /**
     * Constructor.
     *
     * @param i_elTy element type.
     * @param i_nVes number of vertices.
     * @param i_nEls number of elements.
     * @param i_elVe vertices adjacent to the elements.
     * @param i_veCrds vertex coordinates.
     * @param i_inDia length, incircle or insphere diameters of the elements.
     * @param io_velMod velocity model.
     **/
    Cfl( t_entityType  const    i_elTy,
         std::size_t            i_nVes,
         std::size_t            i_nEls,
         std::size_t   const  * i_elVe,
         double        const (* i_veCrds)[3],
         double        const  * i_inDia,
         models::Model        & io_velMod );

    /**
     * Destructor.
     **/
    ~Cfl();

    /**
     * Prints the statistices of the CFL-based time steps.
     **/
    void printStats() const;

    /**
     * Gets the normalized (divided by absolute minimum) time steps.
     *
     * @return normalized time steps.
     **/
    double const * getTimeSteps() const { return m_ts; }
};

#endif