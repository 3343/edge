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
#include "Rule.h"
#include "../../io/logging.h"

edge_v::models::seismic::Rule::Rule( std::string const & i_rule ) {
  std::string l_vpName = "vp";
  std::string l_vsName = "vs";
  std::string l_rhoName = "rho";

  m_exprTk.addVar( l_vpName,
                   m_vp );
  m_exprTk.addVar( l_vsName,
                   m_vs );
  m_exprTk.addVar( l_rhoName,
                   m_rho );

  std::string l_rule = i_rule;
  if( l_rule == "" ) l_rule = "0==0;";

  m_exprTk.compile( l_rule );
}

void edge_v::models::seismic::Rule::apply( float & io_vp,
                                           float & io_vs,
                                           float & io_rho ) {
  // check for valid input
  EDGE_V_CHECK_GT( io_vp, 0 );
  EDGE_V_CHECK_GT( io_vs, 0 );
  EDGE_V_CHECK_GT( io_rho, 0 );

  // copy to local double
  m_vp = io_vp;
  m_vs = io_vs;
  m_rho = io_rho;

  // eval expression
  m_exprTk.eval();

  // copy back to float
  io_vp = m_vp;
  io_vs = m_vs;
  io_rho = m_rho;
}