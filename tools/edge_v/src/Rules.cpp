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
#include "Rules.h"
#include <cassert>

void edge_v::vel::Rules::tpv34( float &io_vp,
                                float &io_vs,
                                float &io_rho ) {
  if( io_vp <= 2984.0 || io_vs <= 1400.0 ) {
    io_vp =  2984;
    io_vs =  1400;
    io_rho = 2220.34;
  }
}

void edge_v::vel::Rules::highf2018( float &io_vp,
                                    float &io_vs,
                                    float &io_rho ) {
  // limit vs to 500 m/s
  if( io_vs < 500 ) {
    //  derive vp/vs ratio
    float l_sca = io_vp / io_vs;
    // adjust velocities
    io_vs = 500;
    io_vp = 500*l_sca;
  }

  // limit vp to 1500 m/s
  if( io_vp < 1500 )
    io_vp = 1500;

  // limit vp/vs ratio to 3
  float l_sca = io_vp / io_vs;
  if( l_sca > 3 ) {
    io_vp = 3*io_vs;
  }

  // adjust lambda
  l_sca = io_vp / io_vs;
  if( l_sca < float(1.45) )
    io_vp = float(1.45) * io_vs;
}

void edge_v::vel::Rules::apply( std::string &i_rule,
                                float       &io_vp,
                                float       &io_vs,
                                float       &io_rho ) {
  // check for valid input
  assert( io_vp  > 0 );
  assert( io_vs  > 0 );
  assert( io_rho > 0 );

  if( i_rule == "tpv34" ) {
    tpv34( io_vp,
           io_vs,
           io_rho );
  }
  else if( i_rule == "highf2018" ) {
    highf2018( io_vp,
               io_vs,
               io_rho );
  }
}