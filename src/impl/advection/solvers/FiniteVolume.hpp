/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2019, Alexander Breuer
 * Copyright (c) 2016-2018, Regents of the University of California
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
 * Finite Volume solver for the advection equation: f(q) = a * q || f(q) = a*q + b*q
 **/
#ifndef EDGE_ADVECTION_FINITE_VOLUME_HPP
#define EDGE_ADVECTION_FINITE_VOLUME_HPP

#include <limits>
#include <cassert>
#include "constants.hpp"

namespace edge {
  namespace advection {
    namespace solvers { 
      class FiniteVolume;
    }
  }
}

class edge::advection::solvers::FiniteVolume {
  public:
    /**
     * Computes the (trivial) time integrated DOFs.
     *
     * @param i_first first element considered.
     * @param i_nElement number of elements.
     * @param i_dT time step ("integrated interval").
     * @param i_dofs DOFs.
     * @parma o_tInt will be set to time integrated DOFs.
     **/
    static void timeInt(       int_el      i_first,
                               int_el      i_nElements,
                               double      i_dT,
                         const real_base (*i_dofs)[1][N_ELEMENT_MODES][N_CRUNS],
                               real_base (**o_tInt)[N_ELEMENT_MODES][N_CRUNS] ) {
      // compute time integated DOFs
      for( int_el l_element = i_first; l_element < i_first+i_nElements; l_element++ ) {
        for( int_md l_mode = 0; l_mode < N_ELEMENT_MODES; l_mode++ ) {
          for( int_cfr l_run = 0; l_run < N_CRUNS; l_run++ ) {
            o_tInt[l_element][0][l_mode][l_run] = i_dT * i_dofs[l_element][0][l_mode][l_run];
          }
        }
      }
    }

    /**
     * Performs a single time step of the advection solver.
     *
     * @param i_first first element considered.
     * @param i_nElement number of elements.
     * @param i_elFaEl ids of face-neighboring elements.
     * @param i_fluxSolvers flux solvers.
     * @param i_tInt time integrated DOFs.
     * @param io_dofs DOFs.
     **/
    static void update(       int_el             i_first,
                              int_el             i_nElements,
                        const int_el           (*i_elFaEl)[C_ENT[T_SDISC.ELEMENT].N_FACES],
                        const t_elementShared2 (*i_fluxSolvers)[ C_ENT[T_SDISC.ELEMENT].N_FACES*2 ],
                              real_base        (**i_tInt)[1][N_CRUNS],
                              real_base        (*io_dofs)[1][1][N_CRUNS] ) {
      // iterate over elements
      for( int_el l_el = i_first; l_el < i_first+i_nElements; l_el++ ) {
          // fluxes
          for( int_el l_fa = 0; l_fa < C_ENT[T_SDISC.ELEMENT].N_FACES; l_fa++ ) {
            // derive neighbor
            int_el l_neigh = i_elFaEl[l_el][l_fa];

            for( int_cfr l_run = 0; l_run < N_CRUNS; l_run++ ) {
              // local cont
              io_dofs[l_el][0][0][l_run] += i_fluxSolvers[l_el][l_fa] * i_tInt[l_el][0][0][l_run];

              // neighboring cont
              io_dofs[l_el][0][0][l_run] += i_fluxSolvers[l_el][C_ENT[T_SDISC.ELEMENT].N_FACES+l_fa] * i_tInt[l_neigh][0][0][l_run];
            }
          }
      }
    }
};

#endif
