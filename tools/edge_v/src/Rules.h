/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
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
#ifndef EDGE_V_VEL_RULES_H
#define EDGE_V_VEL_RULES_H

#include <string>

namespace edge_v {
  namespace vel {
    class Rules;
  }
}

/**
 * @brief Runtime configuration.
 */
class edge_v::vel::Rules {
  private:
    /**
     * @brief TPV34 benchmark.
     * 
     * @param io_vp p-wave velocity.
     * @param io_vs s-wave velocity.
     * @param io_rho density.
     */
    static void tpv34( double &io_vp,
                       double &io_vs,
                       double &io_rho );

  public:
    /**
     * @brief Applies the velocity rule.
     * 
     * @param i_rule string representation of the rule, which gets applied.
     * @param io_vp p-wave velocity.
     * @param io_vs s-wave velocity.
     * @param io_rho density.
     */
    static void apply( std::string &i_rule,
                       double      &io_vp,
                       double      &io_vs,
                       double      &io_rho );
};

#endif