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
 * Writes background meshes in msh4 format.
 **/
#include "BgMeshMsh4.h"
#include "logging.h"

void edge_v::io::BgMeshMsh4::write( t_entityType               i_elTy,
                                    t_idx                      i_nVes,
                                    t_idx                      i_nEls,
                                    t_idx             const  * i_elVe,
                                    double            const (* i_veCrds)[3],
                                    float             const  * i_targetLengths,
                                    std::ostream             & io_stream ) {
  // write header
  io_stream << "$MeshFormat\n";
  io_stream << "4.1 0 8\n";
  io_stream << "$EndMeshFormat\n";

  // write vertices
  unsigned short l_nDis = CE_N_DIS(i_elTy);
  io_stream << "$Nodes\n";
  io_stream << "1 " << i_nVes << " 0 " << i_nVes-1 << "\n";
  io_stream << l_nDis << " 1 " << 0 << " " << i_nVes << "\n";
  for( t_idx l_ve = 0; l_ve < i_nVes; l_ve++ ) {
    io_stream << l_ve << "\n";
  }
  for( t_idx l_ve = 0; l_ve < i_nVes; l_ve++ ) {
    io_stream << i_veCrds[l_ve][0] << " "
              << i_veCrds[l_ve][1] << " "
              << i_veCrds[l_ve][2] << "\n";
  }
  io_stream << "$EndNodes\n";

  // write elements
  unsigned short l_nElVes = CE_N_VES(i_elTy);
  unsigned short l_elTyGmsh = 0;
  if(      i_elTy == TRIA3 ) l_elTyGmsh = 2;
  else if( i_elTy == TET4  ) l_elTyGmsh = 4;
  else EDGE_V_LOG_FATAL;

  io_stream << "$Elements\n";
  io_stream << "1 " << i_nEls << " 0 " << i_nEls-1 << "\n";
  io_stream << l_nDis << " 1 " << l_elTyGmsh << " " << i_nEls << "\n";
  for( t_idx l_el = 0; l_el < i_nEls; l_el++ ) {
    io_stream << l_el;
    for( unsigned short l_ve = 0; l_ve < l_nElVes; l_ve++ ) {
      io_stream << " " << i_elVe[l_el*l_nElVes + l_ve];
    }
    io_stream << "\n";
  }
  io_stream << "$EndElements\n";

  // write the values of the background mesh
  io_stream << "$NodeData\n";
  io_stream << "1\n";
  io_stream << "\"EDGE-V Background Mesh (git-version: " << PP_EDGE_VERSION << ")\"\n";
  io_stream << "1\n";
  io_stream << "0.0\n";
  io_stream << "3\n";
  io_stream << "0\n";
  io_stream << "1\n";
  io_stream << i_nVes << "\n";
  for( t_idx l_ve = 0; l_ve < i_nVes; l_ve++ ) {
    io_stream << l_ve << " " << i_targetLengths[l_ve] << "\n";
  }
  io_stream << "$EndNodeData\n";
}