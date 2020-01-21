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

void edge_v::mesh::Refinement::init( std::size_t             i_nVes,
                                     std::size_t             i_nEls,
                                     unsigned short          i_nElVes,
                                     std::size_t    const  * i_elVe,
                                     double         const (* i_veCrds)[3],
                                     std::string    const  & i_refExpr,
                                     models::Model  const  & i_velMod ) {
  double *l_scale = new double[i_nEls];

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
  std::string l_elspwlName = "elements_per_wave_length";
  std::string l_maxWsRatioName = "maximum_wave_speed_ratio";

  // set up expression-tk interface
  io::ExprTk m_exprTk;
  for( unsigned short l_di = 0; l_di < 3; l_di++ ) {
    m_exprTk.addVar( l_crdsName[l_di],
                     l_crds[l_di] );
  }
  m_exprTk.addVar( l_maxWsRatioName,
                   l_maxWsRatio );

  m_exprTk.addVar( l_elspwlName,
                   l_elspwl );
  m_exprTk.addVar( l_freqName,
                   l_freq );
  m_exprTk.compile( i_refExpr );

  // derive refinement
#ifdef PP_USE_OMP
#pragma omp for
#endif
  for( std::size_t l_el = 0; l_el < i_nEls; l_el++ ) {
    // derive averaged coordinates and wave speed ratio
    for( unsigned short l_di = 0; l_di < 3; l_di++ ) l_crds[l_di] = 0;
    double l_wsMin = std::numeric_limits< double >::max();
    double l_wsMax = std::numeric_limits< double >::lowest();

    for( unsigned short l_ve = 0; l_ve < i_nElVes; l_ve++ ) {
      std::size_t l_veId = i_elVe[l_el * i_nElVes + l_ve];
      for( unsigned short l_di = 0; l_di < 3; l_di++ ) {
        l_crds[l_di] += i_veCrds[l_veId][l_di] / i_nElVes;
      }

      l_wsMin = std::min( l_wsMin, i_velMod.getMaxSpeed( l_veId ) );
      l_wsMax = std::max( l_wsMax, i_velMod.getMaxSpeed( l_veId ) );
    }

    l_maxWsRatio = l_wsMax / l_wsMin;

    // eval expression
    m_exprTk.eval();
    EDGE_V_CHECK_GT( l_freq, 0 );
    EDGE_V_CHECK_GT( l_elspwl, 0 );

    // store result
    l_scale[l_el] = 1.0;
    l_scale[l_el] /= l_freq * l_elspwl;
  }

#ifdef PP_USE_OMP
}
#endif

  // allocate output memory
  free();
  m_ref = new float[i_nVes];

  for( std::size_t l_ve = 0; l_ve < i_nVes; l_ve++ ) m_ref[l_ve] = std::numeric_limits< double >::max();

  // propagate refinement to vertices
  for( std::size_t l_el = 0; l_el < i_nEls; l_el++ ) {
    for( unsigned short l_ve = 0; l_ve < i_nElVes; l_ve++ ) {
      std::size_t l_veId = i_elVe[l_el * i_nElVes + l_ve];

      double l_ref = i_velMod.getMaxSpeed( l_veId );
             l_ref *= l_scale[l_el];
      m_ref[l_veId] = std::min( m_ref[l_veId], float(l_ref) );
    }
  }

  // free temporary memory
  delete[] l_scale;
}