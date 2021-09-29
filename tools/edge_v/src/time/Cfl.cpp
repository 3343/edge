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
 * Time steps, based on stability constraints.
 **/
#include "Cfl.h"

#include "../io/logging.h"

edge_v::time::Cfl::Cfl( t_entityType  const    i_elTy,
                        t_idx                  i_nVes,
                        t_idx                  i_nEls,
                        t_idx         const  * i_elVe,
                        double        const (* i_veCrds)[3],
                        double        const  * i_inDia,
                        bool                   i_velModEl,
                        models::Model        & io_velMod ) {
  m_nEls = i_nEls;

  // allocate memory
  t_idx l_size = m_nEls;
  m_ts = new double[l_size];

  // set the time steps
  setTimeSteps( i_elTy,
                i_nVes,
                i_nEls,
                i_elVe,
                i_veCrds,
                i_inDia,
                i_velModEl,
                io_velMod,
                m_tsAbsMin,
                m_ts );
}

edge_v::time::Cfl::~Cfl() {
  delete[] m_ts;
}

void edge_v::time::Cfl::printStats() const {
  EDGE_V_LOG_INFO << "printing CFL-based time step stats (absolute / relative):";

  // collect min, mean, max
  double l_mean = 0;
  double l_max = std::numeric_limits< double >::lowest();

  for( t_idx l_el = 0; l_el < m_nEls; l_el++ ) {
    l_mean += m_ts[l_el];
    l_max = std::max( m_ts[l_el], l_max );
  }
  l_mean /= m_nEls;

  EDGE_V_LOG_INFO << "  min:  " << m_tsAbsMin        << " / 1.0";
  EDGE_V_LOG_INFO << "  mean: " << l_mean*m_tsAbsMin << " / " << l_mean;
  EDGE_V_LOG_INFO << "  max:  " << l_max*m_tsAbsMin  << " / " << l_max;
}

void edge_v::time::Cfl::setTimeSteps( t_entityType  const    i_elTy,
                                      t_idx                  i_nVes,
                                      t_idx                  i_nEls,
                                      t_idx         const  * i_elVe,
                                      double        const (* i_veCrds)[3],
                                      double        const  * i_inDia,
                                      bool                   i_velModEl,
                                      models::Model        & io_velMod,
                                      double               & o_tsAbsMin,
                                      double               * o_ts ) {
  unsigned short l_nElVes = CE_N_VES( i_elTy );

  o_tsAbsMin = std::numeric_limits< double >::max();

#ifdef PP_USE_OMP
#pragma omp parallel for reduction(min:o_tsAbsMin)
#endif
  // compute absolute time steps
  for( t_idx l_el = 0; l_el < i_nEls; l_el++ ) {
    double l_cMean = 0;
    if( i_velModEl ) {
      l_cMean = io_velMod.getMaxSpeed( l_el );
    }
    else {
      // iterate over vertices and determine mean wave speed
      for( unsigned short l_ve = 0; l_ve < l_nElVes; l_ve++ ) {
        // get the id of the vertex
        t_idx l_veId = i_elVe[ l_el*l_nElVes + l_ve ];

        // add to velocity
        l_cMean += io_velMod.getMaxSpeed( l_veId );
      }
      l_cMean /= l_nElVes;
    }
    EDGE_V_CHECK_GT( l_cMean, 0 );

    // compute time step
    o_ts[l_el] = i_inDia[l_el] / l_cMean;

    // update minimum
    o_tsAbsMin = std::min( o_tsAbsMin, o_ts[l_el] );
  }
  EDGE_V_CHECK_GT( o_tsAbsMin, 0 );

  // normalize
  double l_minInv = 1 / o_tsAbsMin;

#ifdef PP_USE_OMP
#pragma omp parallel for
#endif
  for( t_idx l_el = 0; l_el < i_nEls; l_el++ ) {
    o_ts[l_el] *= l_minInv;
  }
}