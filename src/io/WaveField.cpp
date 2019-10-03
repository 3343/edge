/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2019, Alexander Breuer
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
 * Wave field writer.
 **/
#include <algorithm>
#include <cassert>
#include <sys/types.h>
#include <cstring>

#include "WaveField.h"
#include "logging.h"
#include "FileSystem.hpp"
#include "monitor/instrument.hpp"

edge::io::WaveField::WaveField( std::string             i_type,
                                std::string             i_outFile,
                                std::size_t             i_nVes,
                                t_enLayout     const  & i_elLayout,
                                t_vertexChars  const  * i_veChars,
                                t_elementChars const  * i_elChars,
                                std::size_t    const (* i_elVe)[C_ENT[T_SDISC.ELEMENT].N_VERTICES],
                                real_base      const (* i_dofs)[N_QUANTITIES][N_ELEMENT_MODES][N_CRUNS],
                                int_spType              i_spType ):
 m_veChars(i_veChars),
 m_elChars(i_elChars),
 m_elVe(i_elVe),
 m_dofs(i_dofs),
 m_nVes(i_nVes) {
  // derive the print elements
  for( int_tg l_tg = 0; l_tg < i_elLayout.timeGroups.size(); l_tg++ ) {
    std::size_t l_first = i_elLayout.timeGroups[l_tg].inner.first;
    std::size_t l_size  = i_elLayout.timeGroups[l_tg].nEntsOwn;
    // iterate over owned elements of this time group
    for( std::size_t l_el = l_first; l_el < l_first+l_size; l_el++ ) {
      // only add if element has the desired sparse type
      if(     i_spType == std::numeric_limits< int_spType >::max()
          || (i_elChars[l_el].spType & i_spType) == i_spType )
        m_elPrint.push_back( l_el );
    }
  }
  // remove duplicates
  std::sort( m_elPrint.begin(), m_elPrint.end() );
  m_elPrint.erase( unique( m_elPrint.begin(), m_elPrint.end() ), m_elPrint.end() );

  /*
   * derive sparse ids of limited elements.
   */
  // last print element
  std::size_t l_lastPrint = (m_elPrint.size() == 0) ? 0 : m_elPrint.back()+1;

  // current print element
  std::size_t l_ep = 0;

  // iterate over dense elements
  for( std::size_t l_el = 0; l_el < l_lastPrint; l_el++ ) {
    // ignore if not print
    if( l_el != m_elPrint[l_ep] ) continue;

    // increase counter for print elements
    l_ep++;

    if( l_ep >= m_elPrint.size() ) break;
  }

  if(      i_type == "vtk_ascii"  ) m_type = vtkAscii;
  else if( i_type == "vtk_binary" ) m_type = vtkBinary;
  else                              m_type = none;

  // create new directory
  if( m_type != none ) {
    EDGE_LOG_INFO << "setting up wave field output";
    std::string l_dir, l_file;
    FileSystem::splitPathLast( i_outFile, l_dir, l_file );

    if( m_elPrint.size() > 0 ) {
      std::string l_dirCreate = l_dir + "/" + std::to_string(parallel::g_rank);
      FileSystem::createDir( l_dirCreate );

      l_dir = l_dir + "/" + std::to_string(parallel::g_rank) + '/';
      m_outFile = l_dir + l_file;
    }
  }

  m_writeStep = 0;
}

void edge::io::WaveField::write( double i_time ) {
  PP_INSTR_FUN("write_wf")

  // create file name
  std::string l_outFile = m_outFile;
  l_outFile += "_" + parallel::g_rankStr + "_" + std::to_string((unsigned long long) m_writeStep) + ".vtk";

  // write output
  if( m_elPrint.size() > 0 )
    m_vtk.write( l_outFile,
                  m_type==vtkBinary,
                  m_nVes,
                  m_elPrint,
                  m_veChars,
                  m_elVe,
                  m_dofs );

  m_writeStep++;
}
