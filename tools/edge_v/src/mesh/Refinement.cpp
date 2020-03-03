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
#include "../io/logging.h"
#include "../io/ExprTk.h"

void edge_v::mesh::Refinement::free() {
  if( m_refVe != nullptr ) delete[] m_refVe;
  if( m_refEl != nullptr ) delete[] m_refEl;
}

edge_v::mesh::Refinement::~Refinement() {
  free();
}

void edge_v::mesh::Refinement::init( t_idx                   i_nVes,
                                     t_idx                   i_nEls,
                                     unsigned short          i_nElVes,
                                     t_idx          const  * i_elVe,
                                     double         const (* i_veCrds)[3],
                                     std::string    const  & i_refExpr,
                                     models::Model  const  & i_velMod ) {
  // allocate output memory
  free();
  m_refVe = new float[i_nVes];
  m_refEl = new float[i_nEls];

#ifdef PP_USE_OMP
#pragma omp parallel
{
#endif
  // variables
  double l_crds[3] = { std::numeric_limits< double >::max(),
                       std::numeric_limits< double >::max(),
                       std::numeric_limits< double >::max() };
  double l_elspwl = std::numeric_limits< double >::max();
  double l_freq = std::numeric_limits< double >::max();
  double l_maxWsRatio = std::numeric_limits< double >::max();

  // names of the variables
  std::string l_crdsName[3] = {"x", "y", "z"};
  std::string l_freqName = "frequency";
  std::string l_edspwlName = "edges_per_wave_length";
  std::string l_maxWsRatioName = "maximum_wave_speed_ratio";

  // set up expression-tk interface
  io::ExprTk l_exprTk;
  for( unsigned short l_di = 0; l_di < 3; l_di++ ) {
    l_exprTk.addVar( l_crdsName[l_di],
                     l_crds[l_di] );
  }
  l_exprTk.addVar( l_maxWsRatioName,
                   l_maxWsRatio );
  l_exprTk.addVar( l_edspwlName,
                   l_elspwl );
  l_exprTk.addVar( l_freqName,
                   l_freq );
  l_exprTk.compile( i_refExpr );

  // derive refinement
#ifdef PP_USE_OMP
#pragma omp for
#endif
  for( t_idx l_el = 0; l_el < i_nEls; l_el++ ) {
    // derive averaged coordinates and wave speed ratio
    for( unsigned short l_di = 0; l_di < 3; l_di++ ) l_crds[l_di] = 0;
    double l_wsMin = std::numeric_limits< double >::max();
    double l_wsMax = std::numeric_limits< double >::lowest();

    for( unsigned short l_ve = 0; l_ve < i_nElVes; l_ve++ ) {
      t_idx l_veId = i_elVe[l_el * i_nElVes + l_ve];
      for( unsigned short l_di = 0; l_di < 3; l_di++ ) {
        l_crds[l_di] += i_veCrds[l_veId][l_di] / i_nElVes;
      }

      l_wsMin = std::min( l_wsMin, i_velMod.getMinSpeed( l_veId ) );
      l_wsMax = std::max( l_wsMax, i_velMod.getMinSpeed( l_veId ) );
    }

    l_maxWsRatio = l_wsMax / l_wsMin;

    // eval expression
    l_exprTk.eval();
    EDGE_V_CHECK_GT( l_freq, 0 );
    EDGE_V_CHECK_GT( l_elspwl, 0 );

    // store result
    m_refEl[l_el] = 1.0;
    m_refEl[l_el] /= l_freq * l_elspwl;
  }

#ifdef PP_USE_OMP
}
#endif
  for( t_idx l_ve = 0; l_ve < i_nVes; l_ve++ ) m_refVe[l_ve] = std::numeric_limits< float >::max();

  // propagate refinement to vertices
  for( t_idx l_el = 0; l_el < i_nEls; l_el++ ) {
    float l_wsAve = 0;
    for( unsigned short l_ve = 0; l_ve < i_nElVes; l_ve++ ) {
      t_idx l_veId = i_elVe[l_el * i_nElVes + l_ve];

      float l_ws = i_velMod.getMinSpeed( l_veId );
      l_wsAve += l_ws;

      float l_refVe = l_ws * m_refEl[l_el];
      m_refVe[l_veId] = std::min( m_refVe[l_veId], l_refVe );
    }

    l_wsAve /= i_nElVes;
    m_refEl[l_el] *= l_wsAve;
  }
}
