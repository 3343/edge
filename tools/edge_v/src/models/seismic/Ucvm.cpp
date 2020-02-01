/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2019-2020, Alexander Breuer
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
 * Ucvm velocity model.
 **/
#include "Ucvm.h"

edge_v::models::seismic::Ucvm::Ucvm( io::Ucvm          & i_ucvmReader,
                                     double      const   i_trafoSrc[3][3],
                                     std::string const & i_projSrc,
                                     std::string const & i_projDes,
                                     std::string const & i_ucvmType ): m_ucvmReader( i_ucvmReader ) {
  // store info
  for( unsigned short l_d0 = 0; l_d0 < 3; l_d0++ )
    for( unsigned short l_d1 = 0; l_d1 < 3; l_d1++ )
      m_trafoSrc[l_d0][l_d1] = i_trafoSrc[l_d0][l_d1];

  m_projSrc = i_projSrc;
  m_projDes = i_projDes;
  m_ucvmType = i_ucvmType;
}

void edge_v::models::seismic::Ucvm::free() {
  if( m_velP != nullptr ) delete[] m_velP;
  if( m_velS != nullptr ) delete[] m_velS;
  if( m_rho  != nullptr ) delete[] m_rho;
}

edge_v::models::seismic::Ucvm::~Ucvm() {
  // free memory if allocated
  free();
}

void edge_v::models::seismic::Ucvm::init( std::size_t          i_nPts,
                                          double      const (* i_pts)[3] ) {
  // free memory if allocated
  free();

  // allocate memory
  m_velP = new float[i_nPts];
  m_velS = new float[i_nPts];
  m_rho  = new float[i_nPts];

  m_ucvmReader.getVels( i_nPts,
                        m_trafoSrc,
                        m_projSrc,
                        m_projDes,
                        m_ucvmType,
                        i_pts,
                        m_velP,
                        m_velS,
                        m_rho );
}

double edge_v::models::seismic::Ucvm::getMinSpeed( std::size_t i_pt ) const {
  return std::abs( m_velS[i_pt] );
}

double edge_v::models::seismic::Ucvm::getMaxSpeed( std::size_t i_pt ) const {
  return std::abs( m_velP[i_pt] );
}