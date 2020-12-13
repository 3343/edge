/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (alex.breuer AT uni-jena.de)
 *
 * @section LICENSE
 * Copyright (c) 2020, Friedrich Schiller University Jena
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
 * Expression-based velocity model using an input grid.
 **/
#include "GridExpression.h"
#include "../io/ExprTk.h"

edge_v::models::GridExpression::GridExpression( io::Grid    const * i_grid,
                                                std::string const & i_expr ) {
  m_grid = i_grid;

  if( i_expr != "" ) m_expr.str = i_expr;
}

void edge_v::models::GridExpression::free() {
  if( m_speedsMin != nullptr ) delete[] m_speedsMin;
  if( m_speedsMax != nullptr ) delete[] m_speedsMax;
}

edge_v::models::GridExpression::~GridExpression() {
  free();
}

void edge_v::models::GridExpression::init( t_idx           i_nPts,
                                           double const (* i_pts)[3] ) {

  // free memory if allocated
  free();

  // allocate memory for the speeds
  m_speedsMin = new double[i_nPts];
  m_speedsMax = new double[i_nPts];

#ifdef PP_USE_OMP
#pragma omp parallel
#endif
  {
    // set up expression
    io::ExprTk l_expr;
    double l_crds[3] = {0, 0, 0};
    double l_data = 0;
    double l_speeds[2] = {0, 0};

    l_expr.addVar( m_expr.xName,
                   l_crds[0] );
    l_expr.addVar( m_expr.yName,
                   l_crds[1] );
    l_expr.addVar( m_expr.zName,
                   l_crds[2] );
    l_expr.addVar( m_expr.dataName,
                   l_data );
    l_expr.addVar( m_expr.minSpeedName,
                   l_speeds[0] );
    l_expr.addVar( m_expr.maxSpeedName,
                   l_speeds[1] );

    l_expr.compile( m_expr.str );

    // eval expression at points
#ifdef PP_USE_OMP
#pragma omp for
#endif
    for( t_idx l_pt = 0; l_pt < i_nPts; l_pt++ ) {
      for( unsigned short l_di = 0; l_di < 3; l_di++ ) {
        l_crds[l_di] = i_pts[l_pt][l_di];
      }
      l_data = m_grid->getData()[l_pt];
      l_expr.eval();
      m_speedsMin[l_pt] = l_speeds[0];
      m_speedsMax[l_pt] = l_speeds[1];
    }

  }
}

double edge_v::models::GridExpression::getMinSpeed( t_idx i_pt ) const {
  return m_speedsMin[i_pt];
}

double edge_v::models::GridExpression::getMaxSpeed( t_idx i_pt ) const {
  return m_speedsMax[i_pt];
}