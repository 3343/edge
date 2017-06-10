/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2016, Regents of the University of California
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
 * Wave field writer.
 **/

#include "WaveField.h"

#include <cassert>
#include <io/logging.h>
#include <sys/types.h>
#include <cstring>
#include "io/FileSystem.hpp"
#include "monitor/instrument.hpp"

edge::io::WaveField::WaveField(       std::string      i_type,
                                      std::string      i_outFile,
                                const t_enLayout      &i_elLayout,
                                const t_inMap         *i_inMap,
                                const t_vertexChars   *i_veChars,
                                const int_el         (*i_elVe)[C_ENT[T_SDISC.ELEMENT].N_VERTICES],
                                const real_base      (*i_dofs)[N_QUANTITIES][N_ELEMENT_MODES][N_CRUNS] ):

  m_veChars(i_veChars), m_elVe(i_elVe), m_dofs(i_dofs) {
//  if(      i_type == "netcdf"     ) m_type = netcdf;
  if( i_type == "vtk_ascii"  )      m_type = vtkAscii;
  else if( i_type == "vtk_binary" ) m_type = vtkBinary;
  else                              m_type = none;

  // create new directory only for non-empty paths
  if( m_type != none ) {
    EDGE_LOG_INFO << "setting up wave field output";
    FileSystem::createDir( i_outFile );
  }

  m_nVe = i_inMap->veDaMe.size();

  m_writeStep  = 0;
  m_outFile    = i_outFile;

  // derive the print elements
  for( int_tg l_tg = 0; l_tg < i_elLayout.timeGroups.size(); l_tg++ ) {
    int_el l_first = i_elLayout.timeGroups[l_tg].inner.first;
    int_el l_size  = i_elLayout.timeGroups[l_tg].nEntsOwn;
    // iterate over owned elements of this time group
    for( int_el l_el = l_first; l_el < l_first+l_size; l_el++ ) {
      // get unique element
      int_el l_elUn = i_inMap->elDaMe[ l_el   ];
             l_elUn = i_inMap->elMeDa[ l_elUn ];

      m_elPrint.push_back( l_elUn );
    }
  }
  // remove duplicates
  std::sort( m_elPrint.begin(), m_elPrint.end() );
  m_elPrint.erase( unique( m_elPrint.begin(), m_elPrint.end() ), m_elPrint.end() );
}

void edge::io::WaveField::write( double i_time ) {
  PP_INSTR_FUN("write_wf")

//  if( m_type == netcdf ) writeNetcdf( i_time, i_dofs );
  if ( m_type == vtkAscii || m_type == vtkBinary ) {
    // create file name
    std::string l_outFile = m_outFile;
    l_outFile += "_" + parallel::g_rankStr + "_" + std::to_string((unsigned long long) m_writeStep) + ".vtk";

    // write output
    m_vtk.write( l_outFile,
                 m_type==vtkBinary,
                 m_nVe,
                 m_elPrint,
                 m_veChars,
                 m_elVe,
                 m_dofs );
  }

  m_writeStep++;
}
