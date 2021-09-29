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
 * Expression-based velocity model.
 **/
#ifndef EDGE_V_MODELS_SEISMIC_EXPRESSION_H
#define EDGE_V_MODELS_SEISMIC_EXPRESSION_H

#include "../Model.h"
#include <string>

namespace edge_v {
  namespace models {
    namespace seismic {
      class Expression;
    }
  }
}

/**
 * Expression-based seismic velocity model.
 **/
class edge_v::models::seismic::Expression: public Model {
  private:
    //! expression which is evaluated in the velocity model
    std::string m_expr;

    //! p-wave velocities
    double * m_vp = nullptr;

    //! s-wave velocities
    double * m_vs = nullptr;

    //! quality factor qp
    double * m_qp = nullptr;

    //! quality factor qs
    double * m_qs = nullptr;

  public:
    /**
     * Constructor.
     *
     * @param i_expr relevant expression.
     **/
    Expression( std::string const & i_expr );

    /**
     * Destructor.
     **/
    ~Expression(); 

    /**
     * Inits the velocity model by evaluating the expressions at the given points.
     *
     * @param i_nPts number of points.
     * @param i_pts coordinates of the points.
     **/
    void init( t_idx           i_nPts,
               double const (* i_pts)[3] );

    /**
     * Frees point-related data (allocated and initialized in the init call).
     **/
    void free();

    /**
     * Gets the minimum wave speed at a point.
     *
     * @param i_pt point at which the minimum wave speed is derived.
     **/
    double getMinSpeed( t_idx i_pt ) const;

    /**
     * Gets the maximum wave speed at a point.
     *
     * @param i_pt point at which the maximum wave speed is derived.
     **/
    double getMaxSpeed( t_idx i_pt ) const;
};

#endif