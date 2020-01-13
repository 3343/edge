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
#include "Refinement.h"
#include "io/logging.h"
#include "io/ExprTk.h"

void edge_v::mesh::Refinement::free() {
  if( m_ref != nullptr ) delete[] m_ref;
}

edge_v::mesh::Refinement::~Refinement() {
  free();
}

void edge_v::mesh::Refinement::init( std::size_t            i_nPts,
                                     double        const (* i_pts)[3],
                                     std::string   const  & i_refExpr,
                                     models::Model const  & i_velMod ) {
  // allocate output memory
  free();
  m_ref = new float[i_nPts];

  // variables
  double l_crds[3] = { std::numeric_limits< double >::max(),
                       std::numeric_limits< double >::max(),
                       std::numeric_limits< double >::max() };
  double l_elspwl = std::numeric_limits< double >::max();
  double l_freq = std::numeric_limits< double >::max();

  // names of the variables
  std::string l_crdsName[3] = {"x", "y", "z"};
  std::string l_freqName = "frequency";
  std::string l_elspwlName = "elements_per_wave_length";

  // set up expression-tk interface
  io::ExprTk m_exprTk;
  for( unsigned short l_di = 0; l_di < 3; l_di++ ) {
    m_exprTk.addVar( l_crdsName[l_di],
                     l_crds[l_di] );
  }
  m_exprTk.addVar( l_elspwlName,
                   l_elspwl );
  m_exprTk.addVar( l_freqName,
                   l_freq );
  m_exprTk.compile( i_refExpr );

  // derive refinemnt
  for( std::size_t l_pt = 0; l_pt < i_nPts; l_pt++ ) {
    // copy coords over
    for( unsigned short l_di = 0; l_di < 3; l_di++ ) {
      l_crds[l_di] = i_pts[l_pt][l_di];
    }

    // eval expression at the point
    m_exprTk.eval();
    EDGE_V_CHECK_GT( l_freq, 0 );
    EDGE_V_CHECK_GT( l_elspwl, 0 );

    // query velocity model
    double l_maxWs = i_velMod.getMaxSpeed( l_pt );

    // set target refinement value
    m_ref[l_pt] = l_maxWs / l_freq;
    m_ref[l_pt] *= l_elspwl;
  }
}