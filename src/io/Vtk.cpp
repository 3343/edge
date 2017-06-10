/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2016-2017, Regents of the University of California
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
 * VTK output.
 **/

#include "Vtk.h"
#include "data/common.hpp"
#include "submodules/include/visit_writer.h"

void edge::io::Vtk::init(       int_el                 i_nVe,
                          const std::vector< int_el > &i_elPrint,
                          const t_vertexChars         *i_veChars,
                          const int_el               (*i_elVe)[C_ENT[T_SDISC.ELEMENT].N_VERTICES] ) {
  // allocate buffers for output matching the format of the visit_writer
  m_coordsVe       = (float*)  data::common::allocate( sizeof(float)  * i_nVe * 3                                            );
  m_connElVe       = (int*)    data::common::allocate( sizeof(int)    * i_elPrint.size() * C_ENT[T_SDISC.ELEMENT].N_VERTICES );
  m_dofs           = (float*)  data::common::allocate( sizeof(float)  * i_elPrint.size() * N_QUANTITIES * N_CRUNS            );
  m_dofsPtrs       = (float**) data::common::allocate( sizeof(float*)                    * N_QUANTITIES * N_CRUNS            );

  // set the visit element type dependent on our build config
  if(       T_SDISC.ELEMENT == LINE   ) m_visitElType = VISIT_LINE;
  else if ( T_SDISC.ELEMENT == TRIA3  ) m_visitElType = VISIT_TRIANGLE;
  else if ( T_SDISC.ELEMENT == QUAD4R ) m_visitElType = VISIT_QUAD;
  else if ( T_SDISC.ELEMENT == HEX8R  ) m_visitElType = VISIT_HEXAHEDRON;
  else if ( T_SDISC.ELEMENT == TET4   ) m_visitElType = VISIT_TETRA;
  else {
   EDGE_LOG_FATAL << "missing element type " << T_SDISC.ELEMENT;
  }

  // set up variable names and dofs ptrs
  for( int_cfr l_run = 0; l_run < N_CRUNS; l_run++ ) {
    for( int_md l_q = 0; l_q < N_QUANTITIES; l_q++ ) {
      m_varNames[   l_run*N_QUANTITIES+l_q] =   "crun_" + std::to_string( (unsigned long long) l_run)
                                              + "_var_" + std::to_string( (unsigned long long) l_q);
      m_varNamesC[  l_run*N_QUANTITIES+l_q] = m_varNames[l_run*N_QUANTITIES+l_q].c_str();
      m_dofsPtrs[l_run*N_QUANTITIES+l_q]    = m_dofs+(i_elPrint.size()*(N_QUANTITIES*l_run + l_q));
    }
  }

  // setup vertices coords
  for( int_el l_ve = 0; l_ve < i_nVe; l_ve++ ) {
    for( int l_dim = 0; l_dim < 3; l_dim++ ) {
      m_coordsVe[l_ve*3+l_dim] = i_veChars[l_ve].coords[l_dim];
    }
  }

  // setup connectivity
  for( int_el l_el = 0; l_el < (int_el) i_elPrint.size(); l_el++ ) {
    int_el l_elId = i_elPrint[l_el];
    for( int_md l_ve = 0; l_ve < C_ENT[T_SDISC.ELEMENT].N_VERTICES; l_ve++ ) {
      m_connElVe[ l_el * C_ENT[T_SDISC.ELEMENT].N_VERTICES + l_ve ] = i_elVe[l_elId][l_ve];
    }
  }
}

edge::io::Vtk::~Vtk() {
  if( m_initialized ) {
     data::common::release( m_coordsVe );
     data::common::release( m_connElVe );
     data::common::release( m_dofs     );
     data::common::release( m_dofsPtrs );
  }
}

void edge::io::Vtk::write( const std::string           &i_outFile,
                                 bool                   i_binary,
                                 int_el                 i_nVe,
                           const std::vector< int_el > &i_elPrint,
                           const t_vertexChars         *i_veChars,
                           const int_el               (*i_elVe)[C_ENT[T_SDISC.ELEMENT].N_VERTICES],
                           const real_base            (*i_dofs)[N_QUANTITIES][N_ELEMENT_MODES][N_CRUNS] ) {
  // do the init if not already accomplished
  if( !m_initialized ) {
    init( i_nVe,
          i_elPrint,
          i_veChars,
          i_elVe );

    m_initialized = true;
  }

  // reorder the DOFs and fill the buffer
  for( int_el l_el = 0; l_el < (int_el) i_elPrint.size(); l_el++ ) {
    int_el l_elId = i_elPrint[l_el];
    for( int_md l_q = 0; l_q < N_QUANTITIES; l_q++ ) {
      for( int_cfr l_crun = 0; l_crun < N_CRUNS; l_crun++ ) {
        m_dofs[i_elPrint.size()*(N_QUANTITIES*l_crun + l_q) + l_el] = i_dofs[l_elId][l_q][0][l_crun];
      }
    }
  }

  // write the data, now..
  edge_write_unstructured_mesh( i_outFile.c_str(),
                                i_binary,
                                i_nVe,
                                m_coordsVe,
                                i_elPrint.size(),
                                m_visitElType,
                                m_connElVe,
                                m_varNamesC,
                                m_dofsPtrs );
}
