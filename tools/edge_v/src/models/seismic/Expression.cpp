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
#include "Expression.h"
#include "io/logging.h"
#include "io/ExprTk.h"

edge_v::models::seismic::Expression::Expression( std::string const & i_expr ) {
  m_expr = i_expr;
}

edge_v::models::seismic::Expression::~Expression() {
  free();
}

void edge_v::models::seismic::Expression::init( std::size_t          i_nPts,
                                                double      const (* i_pts)[3] ) {
  // variables in the velocity model
  double l_crds[3] = { std::numeric_limits< double >::max(),
                       std::numeric_limits< double >::max(),
                       std::numeric_limits< double >::max() };
  double l_vp = std::numeric_limits< double >::max();
  double l_vs = std::numeric_limits< double >::max();
  double l_qp = std::numeric_limits< double >::max();
  double l_qs = std::numeric_limits< double >::max();

  // names of the variables
  std::string l_crdsName[3] = {"x", "y", "z"};
  std::string l_vpName = "vp";
  std::string l_vsName = "vs";
  std::string l_qpName = "qp";
  std::string l_qsName = "qs";

  // set up expression-tk interface
  io::ExprTk l_expr;

  for( unsigned short l_di = 0; l_di < 3; l_di++ ) {
    l_expr.addVar( l_crdsName[l_di],
                   l_crds[l_di] );
  }
  l_expr.addVar( l_vpName,
                 l_vp );
  l_expr.addVar( l_vsName,
                 l_vs );
  l_expr.addVar( l_qpName,
                 l_qp );
  l_expr.addVar( l_qsName,
                 l_qs );

  l_expr.compile( m_expr );

  // allocate memory
  EDGE_V_CHECK_EQ( m_vp, nullptr );
  m_vp = new double[ i_nPts ];
  EDGE_V_CHECK_EQ( m_vs, nullptr );
  m_vs = new double[ i_nPts ];
  EDGE_V_CHECK_EQ( m_qp, nullptr );
  m_qp = new double[ i_nPts ];
  EDGE_V_CHECK_EQ( m_qs, nullptr );
  m_qs = new double[ i_nPts ];

  for( std::size_t l_pt = 0; l_pt < i_nPts; l_pt++ ) {
    // copy coords over
    for( unsigned short l_di = 0; l_di < 3; l_di++ ) {
      l_crds[l_di] = i_pts[l_pt][l_di];
    }

    // eval expression at the point
    l_expr.eval();

    // store results
    m_vp[l_pt] = l_vp;
    m_vs[l_pt] = l_vs;
    m_qp[l_pt] = l_qp;
    m_qs[l_pt] = l_qs;
  }
}

void edge_v::models::seismic::Expression::free() {
  if( m_vp != nullptr ) delete[] m_vp;
  if( m_vs != nullptr ) delete[] m_vs;
  if( m_qp != nullptr ) delete[] m_qp;
  if( m_qs != nullptr ) delete[] m_qs;
}

double edge_v::models::seismic::Expression::getMinSpeed( std::size_t i_pt ) const {
  return m_vs[i_pt];
}

double edge_v::models::seismic::Expression::getMaxSpeed( std::size_t i_pt ) const {
  return m_vp[i_pt];
}