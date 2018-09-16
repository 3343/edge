/**
 * @file This file is part of EDGE.
 *
 * @author Rajdeep Konwar (rkonwar AT ucsd.edu)
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
 * This is the main file of EDGE-V.
 **/

#include "io/Config.h"
#include "Rules.h"
#include "io/Moab.hpp"
#include "io/Ucvm.hpp"
#include "io/GmshView.hpp"
#include <iostream>

int main( int i_argc, char **i_argv ) {
  // check input arguments
  if( i_argc != 3 ) {
    std::cerr << "Usage: " << i_argv[0] << " -f edge_v.conf" << std::endl;
    return EXIT_FAILURE;
  } else if( (i_argv == nullptr) || (i_argv[1][0] != '-') || (i_argv[1][1] != 'f') ) {
    std::cerr << "Usage: " << i_argv[0] << " -f edge_v.conf" << std::endl;
    return EXIT_FAILURE;
  }

  // start time
  clock_t l_tp = clock();

  // parse config
  std::string l_configFile = std::string( i_argv[2] );
  edge_v::io::Config l_config( l_configFile );

  // init UCVM
  std::cout << "initializing UCVM" << std::endl;
  edge_v::io::Ucvm l_ucvm( l_config.m_ucvmCfgFn,
                           l_config.m_ucvmModelList,
                           l_config.m_ucvmCmode );

  // init MOAB mesh interface
  std::cout << "reading mesh: " << l_config.m_meshFn << std::endl;
  edge_v::io::Moab< int > l_moab( l_config.m_meshFn );

  // get number of entities
  int l_nEns[4];
  for( unsigned short l_di = 0; l_di < 4; l_di++ )
    l_nEns[l_di] = l_moab.nEnsByDis( l_di );

  std::cout << "entities by number of dimensions: " << std::endl << "  ";
  for( unsigned short l_di = 0; l_di < 4; l_di++ ) {
    std::cout << l_di << ": " << l_nEns[l_di];
    if( l_di < 3 ) std::cout << ", ";
  }
  std::cout << std::endl;

  // get the coordinates of the vertices
  double (*l_veCrds)[3] = new double[ l_nEns[0] ][3];
  l_moab.getVeCrds( l_veCrds );

  // get connectivity
  std::vector< int > l_elVe;
  l_elVe.resize( l_nEns[3]*4 );
  l_moab.getEnVe( "tet4",
                  &l_elVe[0] );

  // assemble characteristic lengths and store Gmsh view, if requested
  if( l_config.m_posFn != "" ) {
    // get velocity model
    float (*l_veVps)  = new float[ l_nEns[0] ];
    float (*l_veVss)  = new float[ l_nEns[0] ];
    float (*l_veRhos) = new float[ l_nEns[0] ];

    // get velocities from UCVM
    l_ucvm.getVels( l_nEns[0],
                    l_config.m_trafo,
                    l_config.m_projMesh,
                    l_config.m_projVel,
                    l_config.m_ucvmType,
                    l_veCrds,
                    l_veVps,
                    l_veVss,
                    l_veRhos );

    // apply velocity rules
    for( int l_ve = 0; l_ve < l_nEns[0]; l_ve++ ) {
      // massage velocities according to given rule
      edge_v::vel::Rules::apply( l_config.m_velRule,
                                 l_veVps[l_ve],
                                 l_veVss[l_ve],
                                 l_veRhos[l_ve] );
    }

    // assemble characteristic lengths
    float (*l_cls)  = new float[ l_nEns[0] ];

    // scale vs to get char lengths
    float l_sca = float(1.0) / l_config.m_elmtsPerWave;
    for( int l_ve = 0; l_ve < l_nEns[0]; l_ve++ ) {
      l_cls[l_ve] = l_veVss[l_ve] * l_sca;
    }

    // write to disk
    edge_v::io::GmshView::write(  l_config.m_posFn,
                                  "tet4",
                                  l_elVe.size()/4,
                                 &l_elVe[0],
                                  l_veCrds,
                                  l_cls );

    // free memory
    delete[] l_cls;
    delete[] l_veVps;
    delete[] l_veVss;
    delete[] l_veRhos;
  }

  // assemble elements' Lame parameters and store annotated mesh, if requested
  if( l_config.m_annoFn != "" ) {
    std::cout << "annotating mesh with velocity model" << std::endl;

    // barycenters
    double (* l_elBars)[3] = new double[ l_nEns[3] ][3];

    // compute bary centers
    for( int l_el = 0; l_el < l_nEns[3]; l_el++ ) {
      for( unsigned short l_di = 0; l_di < 3; l_di++ ) {
        l_elBars[l_el][l_di] = 0;
        for( unsigned short l_ve = 0; l_ve < 4; l_ve++ ) {
          int l_veId = l_elVe[l_el*4+l_ve];
          l_elBars[l_el][l_di] += l_veCrds[l_veId][l_di];
        }
        l_elBars[l_el][l_di] *= 0.25;
      }
    }

    // query UCVM for velocity model at barycenters
    double (*l_elVps) = new double[ l_nEns[3] ];
    double (*l_elVss) = new double[ l_nEns[3] ];
    double (*l_elRhos) = new double[ l_nEns[3] ];

    l_ucvm.getVels( l_nEns[3],
                    l_config.m_trafo,
                    l_config.m_projMesh,
                    l_config.m_projVel,
                    l_config.m_ucvmType,
                    l_elBars,
                    l_elVps,
                    l_elVss,
                    l_elRhos );

    // apply velocity rules
    for( int l_el = 0; l_el < l_nEns[3]; l_el++ ) {
      // massage velocities according to given rule
      edge_v::vel::Rules::apply( l_config.m_velRule,
                                 l_elVps[l_el],
                                 l_elVss[l_el],
                                 l_elRhos[l_el] );
    }

    // derive elements' Lame parameters
    double (*l_elLams) = new double[ l_nEns[3] ];
    double (*l_elMus)  = new double[ l_nEns[3] ];

    for( int l_el = 0; l_el < l_nEns[3]; l_el++ ) {
      l_elMus[l_el] = l_elVss[l_el] * l_elVss[l_el] * l_elRhos[l_el];
      l_elLams[l_el] = l_elVps[l_el] * l_elVps[l_el] * l_elRhos[l_el];
      l_elLams[l_el] -= 2.0 * l_elMus[l_el];
    }

    // store data in MOAB
    l_moab.setEnData( "tet4",
                      "LAMBDA",
                      l_elLams );

    l_moab.setEnData( "tet4",
                      "MU",
                      l_elMus );

    l_moab.setEnData( "tet4",
                      "RHO",
                      l_elRhos );

    std::cout << "writing annotated mesh: " << l_config.m_annoFn << std::endl;
    l_moab.writeMesh( l_config.m_annoFn );

    // free memory
    delete[] l_elBars;
    delete[] l_elLams;
    delete[] l_elMus;
    delete[] l_elRhos;
  }

  // free memory
  delete[] l_veCrds;

  // finish time
  l_tp = clock() - l_tp;
  std::cout << "Time taken: " << (float) l_tp / CLOCKS_PER_SEC << "s\n";

  return EXIT_SUCCESS;
}