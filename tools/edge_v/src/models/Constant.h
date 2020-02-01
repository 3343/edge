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
 * Constant velocity model.
 **/
#ifndef EDGE_V_MODELS_CONSTANT_H
#define EDGE_V_MODELS_CONSTANT_H

#include "Model.h"

namespace edge_v {
  namespace models {
    class Constant;
  }
}

/**
 * Constant velocity model.
 **/
class edge_v::models::Constant: public Model {
  private:
    //! contant velocity
    double m_vel;

  public:
    /**
     * Constructor.
     *
     * @param i_vel velocity of the constant model. 
     **/
    Constant( double i_vel );

    /**
     * Constant model: Dummy function.
     *
     * @param i_nPts ignored.
     * @param i_pts ignored.
     **/
    void init( std::size_t          i_nPts,
               double      const (* i_pts)[3] );

    /**
     * Constant model: Dummy function
     **/
    void free();

    /**
     * Gets the constant wave speed.
     *
     * @param i_pt ignored.
     **/
    double getMinSpeed( std::size_t i_pt ) const;

    /**
     * Gets the constant wave speed.
     *
     * @param i_pt ignored.
     **/
    double getMaxSpeed( std::size_t i_pt ) const;
};

#endif