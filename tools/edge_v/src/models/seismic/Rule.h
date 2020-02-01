/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2020, Alexander Breuer
 * Copyright (c) 2018, Regents of the University of California
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
 * Velocity rules, typically applied in the upper layers.
 **/
#ifndef EDGE_V_VEL_RULE_H
#define EDGE_V_VEL_RULE_H

#include "io/ExprTk.h"
#include <cassert>
#include <string>

namespace edge_v {
  namespace models {
    namespace seismic {
      class Rule;
    }
  }
}

/**
 * Rules for seismic velocity models.
 **/
class edge_v::models::seismic::Rule {
  private:
    //! p-wave velocity in the expression
    double m_vp = 0;

    //! s-wave velocity in the expression
    double m_vs = 0;

    //! density rho in the expression
    double m_rho = 0;

    //! expression which is evaluated
    io::ExprTk m_exprTk;

  public:
    /**
     * Constructor which initializes the expression interface.
     *
     * @param i_rule velocity rule.
     **/
    Rule( std::string const & i_rule );

    /**
     * Applies the velocity rule.
     * 
     * @param io_vp p-wave velocity.
     * @param io_vs s-wave velocity.
     * @param io_rho density.
     **/
    void apply( float & io_vp,
                float & io_vs,
                float & io_rho );
};

#endif