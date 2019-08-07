/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2019, Alexander Breuer
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
 * Mesh functions and storage.
 **/
#include "Mesh.h"

#include "../io/Moab.h"
#include "Geom.h"
#include "io/logging.h"

void edge_v::mesh::Mesh::getEnVeCrds( t_entityType          i_enTy,
                                      std::size_t  const  * i_enVe,
                                      double       const (* i_veCrds)[3],
                                      double             (* o_enVeCrds)[3] ) {
  // number of vertices for the entity type
  unsigned short l_nEnVes = CE_N_VES( i_enTy );

  // gather coordinates
  for( unsigned short l_ve = 0; l_ve < l_nEnVes; l_ve++ ) {
    std::size_t l_veId = i_enVe[l_ve];

    for( unsigned short l_di = 0; l_di < 3; l_di++ ) {
      o_enVeCrds[l_ve][l_di] = i_veCrds[l_veId][l_di];
    }
  }
}

void edge_v::mesh::Mesh::setInDiameter( t_entityType          i_enTy,
                                        std::size_t           i_nEns,
                                        std::size_t  const  * i_enVe,
                                        double       const (* i_veCrds)[3],
                                        double              * o_inDia ) {
  // get the number of vertices
  unsigned short l_nEnVes = CE_N_VES( i_enTy );

  // check the size of the stack-memory below
  EDGE_V_CHECK_LE( l_nEnVes, 8 );

#ifdef PP_USE_OMP
#pragma omp parallel for
#endif
  for( std::size_t l_en = 0; l_en < i_nEns; l_en++ ) {
    // get vertex coordinates
    double l_veCrds[8][3] = {};
    getEnVeCrds( i_enTy,
                 i_enVe + (l_nEnVes * l_en),
                 i_veCrds,
                 l_veCrds );

    // compute diameter
    o_inDia[l_en] = Geom::inDiameter( i_enTy,
                                      l_veCrds );
  }
}

edge_v::mesh::Mesh::Mesh( edge_v::io::Moab const & i_moab ) {
  // get the element type of the mesh
  m_elTy = i_moab.getElType();
  t_entityType l_faTy = CE_T_FA( m_elTy );

  // derive element properties
  unsigned short l_nElVes = CE_N_VES( m_elTy );
  unsigned short l_nElFas = CE_N_FAS( m_elTy );

  // query mesh
  m_nVes = i_moab.nEnsByType( POINT  );
  m_nFas = i_moab.nEnsByType( l_faTy );
  m_nEls = i_moab.nEnsByType( m_elTy );

  // allocate memory
  std::size_t l_size  = m_nVes * 3;
  m_veCrds = (double (*)[3]) new double[ l_size ];

  l_size  = m_nEls * l_nElVes;
  m_elVe = new std::size_t[ l_size ];

  l_size  = m_nEls;
  m_inDiaEl = new double[ l_size ];

  l_size = m_nEls * l_nElFas;
  m_elFaEl = new std::size_t[ l_size ];

  // query moab
  i_moab.getVeCrds( m_veCrds );
  i_moab.getEnVe( m_elTy, m_elVe );
  i_moab.getEnFaEn( m_elTy, m_elFaEl );

  // compute mesh properties
  setInDiameter( m_elTy,
                 m_nEls,
                 m_elVe,
                 m_veCrds,
                 m_inDiaEl );
}

edge_v::mesh::Mesh::~Mesh() {
  delete[] m_veCrds;
  delete[] m_elVe;
  delete[] m_elFaEl;
  delete[] m_inDiaEl;
}

void edge_v::mesh::Mesh::printStats() const {
  EDGE_V_LOG_INFO << "mesh stats:";
  EDGE_V_LOG_INFO << "  #vertices: " << m_nVes;
  EDGE_V_LOG_INFO << "  #faces:    " << m_nFas;
  EDGE_V_LOG_INFO << "  #elements: " << m_nEls;
}