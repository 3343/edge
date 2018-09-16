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

#include "vm_utility.h"
#include "Config.h"
#include "Rules.h"
#include "io/Moab.hpp"
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
  edge_v::vm::Utility::ucvmInit( l_config );

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

  // UCVM interface
  ucvm_point_t *l_ucvmPts = new ucvm_point_t[l_nEns[0]];
  ucvm_data_t *l_ucvmData = new ucvm_data_t[ l_nEns[0]];

  // get the coordinates of the vertices
  double (*l_veCrds)[3] = new double[ l_nEns[0] ][3];
  l_moab.getVeCrds( l_veCrds );

  // get connectivity
  std::vector< int > l_elVe;
  l_elVe.resize( l_nEns[3]*4 );
  l_moab.getEnVe( "tet4",
                  &l_elVe[0] );

  // for all vertices in mesh, populate ucvm points with coords
  for( int l_ve = 0; l_ve < l_nEns[0]; l_ve++ ) {
    // apply trafo
    double l_tmp[3] = {0,0,0};
    for( unsigned short l_d1 = 0; l_d1 < 3; l_d1++ )
      for( unsigned short l_d2 = 0; l_d2 < 3; l_d2++ )
        l_tmp[l_d1] += l_config.m_trafo[l_d1][l_d2] * l_veCrds[l_ve][l_d2];
    for( unsigned short l_di = 0; l_di < 3; l_di++ )
      l_ucvmPts[l_ve].coord[l_di] = l_tmp[l_di];
  }

  // proj.4 projections
  projPJ l_pjSrc  = pj_init_plus( l_config.m_projMesh.c_str() );
  projPJ l_pjDest = pj_init_plus( l_config.m_projVel.c_str()  );

  // proj.4 trafo: UTM -> long,lat,elev
  double *l_xPtr  = &(l_ucvmPts[0].coord[0]);
  double *l_yPtr  = &(l_ucvmPts[0].coord[1]);
  double *l_zPtr  = &(l_ucvmPts[0].coord[2]);
  int l_pntOffset = sizeof( ucvm_point_t ) / sizeof( double );
  pj_transform( l_pjSrc, l_pjDest, l_nEns[0], l_pntOffset,
                l_xPtr, l_yPtr, l_zPtr );

  // apply rad to degree
  for( int_v l_nid = 0; l_nid < l_nEns[0]; l_nid++ ) {
    l_ucvmPts[l_nid].coord[0] *= RAD_TO_DEG;
    l_ucvmPts[l_nid].coord[1] *= RAD_TO_DEG;
  }

  clock_t l_c = clock();

  // UCVM Query on ucvm points to get ucvm data
  std::cout << "UCVM Query... ";
  std::cout.flush();

  int l_ucvmStatus  = ucvm_query( l_nEns[0], l_ucvmPts, l_ucvmData );
  if( l_ucvmStatus ) {
    std::cout << "Failed!" << std::endl;
    std::cerr << "Error: cannot complete UCVM query." << std::endl;
  }

  std::cout << "Done! ";
  l_c = clock() - l_c;
  std::cout << "(" << (float) l_c / CLOCKS_PER_SEC << "s)" << std::endl;

  // get velocity model
  float (*l_vps)  = new float[ l_nEns[0] ];
  float (*l_vss)  = new float[ l_nEns[0] ];
  float (*l_rhos) = new float[ l_nEns[0] ];

  for( int l_ve = 0; l_ve < l_nEns[0]; l_ve++ ) {
    // get pointer to use properties
    ucvm_prop_t *l_propPtr;
    if( l_config.m_ucvmType == "gtl" )
      l_propPtr = &(l_ucvmData[l_ve].gtl);
    else if( l_config.m_ucvmType == "crust" )
      l_propPtr = &(l_ucvmData[l_ve].crust);
    else
      l_propPtr = &(l_ucvmData[l_ve].cmb);

    // store velocities
    l_vps[l_ve] = l_propPtr->vp;
    l_vss[l_ve] = l_propPtr->vs;
    l_rhos[l_ve] = l_propPtr->rho;

    // massage velocities according to given rule
    edge_v::vel::Rules::apply( l_config.m_velRule,
                               l_vps[l_ve],
                               l_vss[l_ve],
                               l_rhos[l_ve] );
  }

  // assemble characteristic lengths, if necessary
  if( l_config.m_posFn != "" ) {
    // assemble characteristic lengths
    float (*l_cls)  = new float[ l_nEns[0] ];

    // scale vs to get char lengths
    float l_sca = float(1.0) / l_config.m_elmtsPerWave;
    for( int l_ve = 0; l_ve < l_nEns[0]; l_ve++ ) {
      l_cls[l_ve] = l_vss[l_ve] * l_sca;
    }

    // write to disk
    edge_v::io::GmshView::write(  l_config.m_posFn,
                                  "tet4",
                                  l_elVe.size()/4,
                                 &l_elVe[0],
                                  l_veCrds,
                                  l_cls );

    delete[] l_cls;
  }

  // free memory
  delete[] l_veCrds;
  delete[] l_ucvmPts;
  delete[] l_ucvmData;
  delete[] l_vps;
  delete[] l_vss;
  delete[] l_rhos;

  // finish time
  l_tp = clock() - l_tp;
  std::cout << "Time taken: " << (float) l_tp / CLOCKS_PER_SEC << "s\n";

  return EXIT_SUCCESS;
}