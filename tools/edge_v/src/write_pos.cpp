/**
 * @file This file is part of EDGE.
 *
 * @author Rajdeep Konwar (rkonwar AT ucsd.edu)
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
 * This is the main file of write_pos.
 **/

#include "vm_utility.h"
#include "Config.h"

int main( int i_argc, char **i_argv ) {
  //! Check input arguments
  if( i_argc != 3 ) {
    std::cerr << "Usage: " << i_argv[0] << " -f config_file.log" << std::endl;
    return EXIT_FAILURE;
  } else if( (i_argv == nullptr) || (i_argv[1][0] != '-') || (i_argv[1][1] != 'f') ) {
    std::cerr << "Usage: " << i_argv[0] << " -f config_file.log" << std::endl;
    return EXIT_FAILURE;
  }

  //! Start time
  clock_t l_tp = clock();

  //! pos configuration file
  std::string l_configFile = std::string( i_argv[2] );
  edge_v::io::Config l_posCfg( l_configFile );
  edge_v::vm::Utility::ucvmInit( l_posCfg );

  //! MOAB mesh object
  moab_mesh l_msh;
  edge_v::vm::Utility::meshInit( l_msh, l_posCfg );

  //! UCVM pointers
  ucvm_point_t *l_ucvmPts = new ucvm_point_t[l_msh.m_numNodes];
  ucvm_data_t *l_ucvmData = new ucvm_data_t[ l_msh.m_numNodes];

  //! MOAB variables
  moab::EntityID                    l_elmtId,   l_nodeId;
  moab::EntityHandle                l_elmtHndl, l_nodeHndl;
  moab::ErrorCode                   l_rval;
  std::vector< moab::EntityHandle > l_vertices;

  int_v l_eid, l_nid; //! Loop counters

  //! For all nodes in mesh, populate ucvm points with coords
  for( l_nid = 0; l_nid < l_msh.m_numNodes; l_nid++ ) {
    l_nodeId  = l_nid + 1;

    //! Get node handle from its id
    l_rval    = l_msh.m_intf->handle_from_id( moab::MBVERTEX, l_nodeId, l_nodeHndl );
    assert( l_rval == moab::MB_SUCCESS );

    //! Get coordinates of node and populate in ucvm pts
    l_rval    = l_msh.m_intf->get_coords( &l_nodeHndl, 1, l_ucvmPts[l_nid].coord );
    assert( l_rval == moab::MB_SUCCESS );
  }

  //! Proj4 variables
  std::string l_pjInitParams = "+proj=tmerc +units=m +axis=enu +no_defs \
                                +datum=WGS84 +k=0.9996 +lon_0=-117.916 +lat_0=33.933";
  projPJ l_pjSrc  = pj_init_plus( l_pjInitParams.c_str() );
  projPJ l_pjDest = pj_init_plus( "+proj=latlong +datum=WGS84" );

  //! Proj4 Transform: UTM -> long,lat,elev
  double *l_xPtr  = &(l_ucvmPts[0].coord[0]);
  double *l_yPtr  = &(l_ucvmPts[0].coord[1]);
  double *l_zPtr  = &(l_ucvmPts[0].coord[2]);
  int l_pntOffset = sizeof( ucvm_point_t ) / sizeof( double );
  pj_transform( l_pjSrc, l_pjDest, l_msh.m_numNodes, l_pntOffset,
                l_xPtr, l_yPtr, l_zPtr );

  //! Apply Rad to Degree
  for( l_nid = 0; l_nid < l_msh.m_numNodes; l_nid++ ) {
    l_ucvmPts[l_nid].coord[0] *= RAD_TO_DEG;
    l_ucvmPts[l_nid].coord[1] *= RAD_TO_DEG;
  }

  clock_t l_c = clock();

  //! UCVM Query on ucvm points to get ucvm data
  std::cout << "UCVM Query... ";
  std::cout.flush();

  int l_ucvmStatus  = ucvm_query( l_msh.m_numNodes, l_ucvmPts, l_ucvmData );
  if( l_ucvmStatus ) {
    std::cout << "Failed!" << std::endl;
    std::cerr << "Error: cannot complete UCVM query." << std::endl;
  }

  std::cout << "Done! ";
  l_c = clock() - l_c;
  std::cout << "(" << (float) l_c / CLOCKS_PER_SEC << "s)" << std::endl;

  //! pos model
  edge_v::vm::Utility::posModel l_posModel;
  edge_v::vm::Utility::posInit( l_posModel, l_msh );

  ucvm_prop_t *l_propPtr;   //! UCVM prop pointer
  real l_vs;

  for( l_eid = 0; l_eid < l_msh.m_numElmts; l_eid++ ) {
    //! Element id starts index at 1
    l_elmtId  = l_eid + 1;

    //! Get element handle from id
    if( ELMTTYPE == 3 )
      l_rval  = l_msh.m_intf->handle_from_id( moab::MBTRI, l_elmtId, l_elmtHndl );
    else if( ELMTTYPE == 4 )
      l_rval  = l_msh.m_intf->handle_from_id( moab::MBTET, l_elmtId, l_elmtHndl );
    assert( l_rval == moab::MB_SUCCESS );
    l_vertices.clear();

    //! Get handles to vertices of element
    l_rval  = l_msh.m_intf->get_adjacencies( &l_elmtHndl, 1, 0, false, l_vertices );
    assert( l_rval == moab::MB_SUCCESS );
    assert( l_vertices.size() == ELMTTYPE );

    for( unsigned int l_vid = 0; l_vid < ELMTTYPE; l_vid++ ) {
      //! Each vertex corresponds to a node, get node id from handle
      l_nodeId  = l_msh.m_intf->id_from_handle( l_vertices[l_vid] );
      l_nid     = l_nodeId - 1;   //! Node index starts from 0

      //! Get coordinates of node and repopulate ucvm pts
      l_rval    = l_msh.m_intf->get_coords( &l_vertices[l_vid], 1,
                                          l_ucvmPts[l_nid].coord );
      assert( l_rval == moab::MB_SUCCESS );

      //! Store x,y,z coordinates
      l_posModel.m_posList[l_eid].m_xyzPts[l_vid].m_x = l_ucvmPts[l_nid].coord[0];
      l_posModel.m_posList[l_eid].m_xyzPts[l_vid].m_y = l_ucvmPts[l_nid].coord[1];
      l_posModel.m_posList[l_eid].m_xyzPts[l_vid].m_z = l_ucvmPts[l_nid].coord[2];

      //! cmb: combination crustal+gtl
      l_propPtr = &(l_ucvmData[l_nid].cmb);
      l_vs      = l_propPtr->vs;

      //! Cutoff velocity (to avoid very low vs to avoid large elmts per wavelength)
      if( l_vs < l_posCfg.m_minVs )
        l_vs  = l_posCfg.m_minVs;

      //! Rescale elmt-size according to no. of elmts required per wavelength
      if( l_posCfg.m_elmtsPerWave > 1.0 )
        l_vs  /= l_posCfg.m_elmtsPerWave;

      l_posModel.m_posList[l_eid].m_vs[l_vid] = l_vs;
    }
  }

  //! write pos file
  edge_v::vm::Utility::writePos( l_posModel, l_posCfg, l_msh );

  edge_v::vm::Utility::posFinalize( l_posModel );

  //! Memory deallocation
  delete[] l_ucvmPts;
  delete[] l_ucvmData;

  //! Finish time
  l_tp = clock() - l_tp;
  std::cout << "Time taken: " << (float) l_tp / CLOCKS_PER_SEC << "s\n";

  return EXIT_SUCCESS;
}
