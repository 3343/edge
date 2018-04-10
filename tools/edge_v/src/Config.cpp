/**
 * @file This file is part of EDGE.
 *
 * @author Junyi Qiu (juq005 AT ucsd.edu)
 * @author Rajdeep Konwar (rkonwar AT ucsd.edu)
 * @author David Lenz (dlenz AT ucsd.edu)
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2017-2018, Regents of the University of California
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
 * Runtime configuration.
 **/
#include "Config.h"
#include <iostream>
#include <fstream>
#include <cassert>

void edge_v::io::Config::vecStringToDouble( char         i_sep,
                                            std::string &i_string,
                                            double*      o_val ) {
  std::string l_right = i_string;

  unsigned short l_pos = 0;
  while( l_right.size() > 0 ) {
    std::size_t l_split = l_right.find_first_of(i_sep);
    if( l_split < l_right.size() ) {
      std::string l_left =  l_right.substr( 0, l_split );
                  l_right = l_right.substr( l_split+1  );

      o_val[l_pos] = atof(l_left.c_str());

      l_pos++;
    }
    else {
      o_val[l_pos] = atof(l_right.c_str());
      break;
    }
  }
}

edge_v::io::Config::Config( const std::string &i_pathToFile ) {
  clock_t l_t = clock();

  m_antnCfgFn     = i_pathToFile;

  m_minVp         = 0.0;
  m_minVs         = 0.0;
  m_minVs2        = 0.0;
  m_maxVpVsRatio  = 0.0;
  m_elmtsPerWave  = 0.0;

  for( unsigned short l_d1 = 0; l_d1 < 3; l_d1++ ) {
    for( unsigned short l_d2 = 0; l_d2 < 3; l_d2++ )
      m_trafo[l_d1][l_d2] = 0;

    m_trafo[l_d1][l_d1] = 1.0;
  }

  m_tetRefinement = 0;

  std::cout << "Reading Config File: " << i_pathToFile << "... " << std::flush;

  std::ifstream l_mshFs( i_pathToFile.c_str(), std::ios::in );
  if( !l_mshFs.is_open() ) {
    std::cout << "Failed." << std::endl;
    std::cerr << "Error: cannot open the config file." << std::endl;
    exit( EXIT_FAILURE );
  }

  std::string l_lineBuf;
  while( getline( l_mshFs, l_lineBuf ) ) {
    size_t l_i = -1, l_j;

    while( (++l_i < l_lineBuf.length()) && (l_lineBuf[l_i] == ' ') );

    if( (l_i >= l_lineBuf.length()) || (l_lineBuf[l_i] == '#') )
      continue;

    l_j = l_i - 1;

    while( (++l_j < l_lineBuf.length()) && (l_lineBuf[l_j] != '=') );

    if( l_j >= l_lineBuf.length() )
      continue;

    std::string l_varName   = l_lineBuf.substr( l_i, l_j - l_i );
    std::string l_varValue  = l_lineBuf.substr( l_j + 1 );

    if( l_varName.compare( "ucvm_config" ) == 0 )
      m_ucvmCfgFn     = l_varValue;
    else if( l_varName.compare( "ucvm_model_list" ) == 0 )
      m_ucvmModelList = l_varValue;
    else if( l_varName.compare( "ucvm_cmode" ) == 0 ) {
      if( l_varValue.compare( "UCVM_COORD_GEO_DEPTH" ) == 0 )
        m_ucvmCmode   = UCVM_COORD_GEO_DEPTH;
      else if( l_varValue.compare( "UCVM_COORD_GEO_ELEV" ) == 0 )
        m_ucvmCmode   = UCVM_COORD_GEO_ELEV;
      else std::cerr << "unknown UCVM_COORD_GEO mode";
    }
    else if( l_varName.compare( "ucvm_type"        ) == 0 ) m_ucvmType     = atoi( l_varValue.c_str() );
    else if( l_varName.compare( "vel_rule"         ) == 0 ) m_velRule      = l_varValue;
    else if( l_varName.compare( "min_vp"           ) == 0 ) m_minVp        = atof( l_varValue.c_str() );
    else if( l_varName.compare( "min_vs"           ) == 0 ) m_minVs        = atof( l_varValue.c_str() );
    else if( l_varName.compare( "min_vs2"          ) == 0 ) m_minVs2       = atof( l_varValue.c_str() );
    else if( l_varName.compare( "max_vp_vs_ratio"  ) == 0 ) m_maxVpVsRatio = atof( l_varValue.c_str() );
    else if( l_varName.compare( "elmts_per_wave"   ) == 0 ) m_elmtsPerWave = atof( l_varValue.c_str() );
    else if( l_varName.compare( "trafo_x"          ) == 0 ) vecStringToDouble( ' ',  l_varValue, m_trafo[0] );
    else if( l_varName.compare( "trafo_y"          ) == 0 ) vecStringToDouble( ' ',  l_varValue, m_trafo[1] );
    else if( l_varName.compare( "trafo_z"          ) == 0 ) vecStringToDouble( ' ',  l_varValue, m_trafo[2] );
    else if( l_varName.compare( "proj_mesh"        ) == 0 ) m_projMesh     = l_varValue;
    else if( l_varName.compare( "proj_vel"         ) == 0 ) m_projVel      = l_varValue;
    else if( l_varName.compare( "mesh_file"        ) == 0 ) m_meshFn       = l_varValue;
    else if( l_varName.compare( "node_vm_file"     ) == 0 ) m_vmNodeFn     = l_varValue;
    else if( l_varName.compare( "elmt_vm_file"     ) == 0 ) m_vmElmtFn     = l_varValue;
    else if( l_varName.compare( "h5m_file"         ) == 0 ) m_h5mFn        = l_varValue;
    else if( l_varName.compare( "pos_file"         ) == 0 ) m_posFn        = l_varValue;
    else if( l_varName.compare( "fault_input_file" ) == 0 ) m_faultInputFns.push_back( l_varValue );
    else if( l_varName.compare( "tet_refinement"   ) == 0 ) m_tetRefinement = std::stoi( l_varValue );
    else std::cout << "\nUnknown setting (" << l_varName << "). Ignored." << std::endl;
  }

  l_mshFs.close();

  std::cout << "Done! ";
  l_t = clock() - l_t;
  std::cout << "(" << (float) l_t / CLOCKS_PER_SEC << "s)" << std::endl;
}
