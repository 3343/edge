/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2016-2017, Regents of the University of California
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
 * Setup for convergence studies of the advection equation.
 **/

#ifndef CONVERGENCE_HPP
#define CONVERGENCE_HPP

#include "constants.hpp"

namespace edge {
  namespace advection {
    namespace setups {
      class Convergence;
    }
  }
}

class edge::advection::setups::Convergence {
  private:

  public:
    /**
     * Sets constant back ground wave speeds
     *
     * @param i_nElements number of elements.
     * @param o_waveSpeeds will be set to constant wave speeds
     * @param i_cX wave speed in x-direction.
     * @param i_cY wave speed in y-direction (ignored for #dims==1).
     * @param i_cZ wave speed in z-direction (ignored for #dims!=3).
     **/
    static void setConstantSpeed( int_el      i_nElements,
                                  t_bgPars  (*o_waveSpeeds)[1],
                                  real_base   i_cX,
                                  real_base   i_cY = 0,
                                  real_base   i_cZ = 0) {
      for( int_el l_el = 0; l_el < i_nElements; l_el++ ) {
        o_waveSpeeds[l_el][0].a = i_cX;
#if PP_N_DIM > 1
        o_waveSpeeds[l_el][0].b = i_cY;
#endif
#if PP_N_DIM > 2
        o_waveSpeeds[l_el][0].c = i_cZ;
#endif
      }
    }
};

#endif
