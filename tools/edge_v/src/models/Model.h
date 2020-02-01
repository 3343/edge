/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2019-2020, Alexander Breuer
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
 * Base class of a velocity model.
 **/
#ifndef EDGE_V_MODELS_MODEL_H
#define EDGE_V_MODELS_MODEL_H

#include <cstddef>

namespace edge_v {
  namespace models {
    class Model;
  }
}

/**
 * Base velocity model.
 **/
class edge_v::models::Model {
  public:
    /**
     * Virtual destructor for base class.
     **/
    virtual ~Model(){};

    /**
     * Inits the velocity model at the given points.
     *
     * @param i_nPts number of points.
     * @param i_pts coordinates of the points.
     **/
    virtual void init( std::size_t          i_nPts,
                       double      const (* i_pts)[3] ) = 0;

    /**
     * Frees point-related data (allocated and initialized in the init call).
     **/
    virtual void free() = 0;

    /**
     * Gets the minimum wave speed at a point.
     *
     * @param i_pt point at which the minimum wave speed is derived.
     **/
    virtual double getMinSpeed( std::size_t i_pt ) const = 0;

    /**
     * Gets the maximum wave speed at a point.
     *
     * @param i_pt point at which the maximum wave speed is derived.
     **/
    virtual double getMaxSpeed( std::size_t i_pt ) const = 0;

    /**
     * Averages the given velocities at the vertices for the elements.
     *
     * @param i_nElVes number of vertices for every element.
     * @parma i_nEls number of elements.
     * @param i_elVe vertices adjacent to the elements.
     * @param i_velVe velocities at the vertices, which are averaged.
     * @param o_velEl will we set to the averaged velocities of the elements.
     **/
    void getElAve( unsigned short         i_nElVes,
                   std::size_t            i_nEls,
                   std::size_t    const * i_elVe,
                   float          const * i_velVe,
                   float                * o_velEl );
};

#endif
