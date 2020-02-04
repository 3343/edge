/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (breuer AT mytum.de)
 *
 * @section LICENSE
 * Copyright (c) 2020, Alexander Breuer
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
 * Velocity-based mesh refinement.
 **/
#ifndef EDGE_V_MESH_REFINEMENT_H
#define EDGE_V_MESH_REFINEMENT_H

#include <string>
#include "models/Model.h"

namespace edge_v {
  namespace mesh {
    class Refinement;
  }
}

/**
 * Mesh refinement.
 **/
class edge_v::mesh::Refinement {
  private:
    //! evaluated mesh refinement at vertices
    float * m_refVe = nullptr;

    //! evaluated mesh refinement at elements
    float * m_refEl = nullptr;

    /**
     * Frees the memory.
     **/
    void free();

  public:
    /**
     * Destructor.
     **/
    ~Refinement();

    /**
     * Inits the mesh refinement by computing the target refinement for all given elements.
     * The mesh refinement at the vertices is then set to the minimum of all adjacent elements.
     *
     * Inputs for the expression are:
     *   x: x-coordinate (averaged over vertices).
     *   y: y-coordinate (averaged over vertices).
     *   z: z-coordinate (averaged over vertices).
     *   maximum_wave_speed_ratio: Given the per-vertex maximum wave speeds,
     *                             this is the max of these divided by the respective min.
     *
     * @param i_nVes number of vertices.
     * @param i_nEls number of elements.
     * @param i_nElVes number of vertices per element.
     * @param i_elVe vertices adjacent to the elements.
     * @param i_veCrds coordinates of the vertices.
     * @param i_refExpr expression used for the refinement, which defines 'frequency' and 'edges_per_wave_length'.
     * @param i_velMod velocity model.
     **/
    void init( std::size_t             i_nVes,
               std::size_t             i_nEls,
               unsigned short          i_nElVes,
               std::size_t    const  * i_elVe,
               double         const (* i_veCrds)[3],
               std::string    const  & i_refExpr,
               models::Model  const  & i_velMod );

    /**
     * Gets the target lengths at the vertices.
     *
     * @return target lengths.
     **/
    float const * getTargetLengthsVe(){ return m_refVe; }

    /**
     * Gets the target lengths of the elements.
     *
     * @return target lengths.
     **/
    float const * getTargetLengthsEl(){ return m_refEl; }
};

#endif