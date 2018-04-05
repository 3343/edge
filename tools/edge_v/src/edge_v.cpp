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
 * This is the main file of EDGE-V.
 **/

#include "vm_utility.h"
#include "FaultModel.h"

int main( int i_argc, char **i_argv ) {
  if( i_argc != 3 ) {
    std::cerr << "Usage: " << i_argv[0] << " -f config_file.log" << std::endl;
    exit( EXIT_FAILURE );
  } else if( (i_argv == nullptr) || (i_argv[1][0] != '-') || (i_argv[1][1] != 'f') ) {
    std::cerr << "Usage: " << i_argv[0] << " -f config_file.log" << std::endl;
    exit( EXIT_FAILURE );
  }

  //! Start time
  clock_t l_te = clock();

  std::string l_configFile = std::string( i_argv[2] );
  edge_v::io::Config l_aCfg( l_configFile );
  edge_v::vm::Utility::ucvmInit( l_aCfg );

  moab_mesh l_msh;
  edge_v::vm::Utility::meshInit( l_msh, l_aCfg );

  //! Phase-1: Nodes
  edge_v::vm::Utility::vmodel l_vModelNodes;
  edge_v::vm::Utility::vmNodeInit( l_vModelNodes, l_msh );

  ucvm_point_t *l_ucvmPoints  = new ucvm_point_t[l_msh.m_numNodes];
  ucvm_data_t *l_ucvmProps    = new ucvm_data_t[ l_msh.m_numNodes];

  edge_v::vm::Utility::worker_reg l_wrkRg;
  edge_v::vm::Utility::workerInit( l_wrkRg,
                                   l_msh.m_numNodes,
                                   l_aCfg.m_projMesh,
                                   l_aCfg.m_projVel );

  moab::EntityID      l_pEntId;
  moab::EntityHandle  l_pHandle;
  moab::ErrorCode     l_rval;

  int_v l_pid;
  const int_v l_pOfs = l_wrkRg.m_workerTid * l_wrkRg.m_workSize;

  for( l_pid = 0; l_pid < l_wrkRg.m_numPrvt; l_pid++ ) {
    l_pEntId  = l_pid + l_pOfs + 1;
    l_rval    = l_msh.m_intf->handle_from_id( moab::MBVERTEX, l_pEntId, l_pHandle );
    assert( l_rval == moab::MB_SUCCESS );

    l_rval    = l_msh.m_intf->get_coords( &l_pHandle, 1,
                                          l_ucvmPoints[l_pid+l_pOfs].coord );
    assert( l_rval == moab::MB_SUCCESS );
  }

  //! Proj4 Transform: UTM->Long,Lat,Elv
  double *l_xPtr  = &(l_ucvmPoints[l_pOfs].coord[0]);
  double *l_yPtr  = &(l_ucvmPoints[l_pOfs].coord[1]);
  double *l_zPtr  = &(l_ucvmPoints[l_pOfs].coord[2]);
  int l_pntOfs    = sizeof( ucvm_point_t ) / sizeof( double );
  pj_transform( l_wrkRg.m_pjUtm, l_wrkRg.m_pjGeo, l_wrkRg.m_numPrvt, l_pntOfs,
                l_xPtr, l_yPtr, l_zPtr );

  //! Apply Rad to Degree
  for( l_pid = 0; l_pid < l_wrkRg.m_numPrvt; l_pid++ ) {
    l_ucvmPoints[l_pid+l_pOfs].coord[0] *= RAD_TO_DEG;
    l_ucvmPoints[l_pid+l_pOfs].coord[1] *= RAD_TO_DEG;

    //! (Raj): Depth (m) does not need transformation from rad to deg
    // ucvmPoints[l_pid+l_pOfs].m_coord[2] *= RAD_TO_DEG;
  }

  //! UCVM Query
  clock_t l_c = clock();

  std::cout << "UCVM Query... ";
  std::cout.flush();

  int l_ucvmStatus = ucvm_query( l_msh.m_numNodes, l_ucvmPoints, l_ucvmProps );
  if( l_ucvmStatus != 0 ) {
    std::cout << "Failed." << std::endl;
    std::cerr << "Warning: cannot complete UCVM query, won't annotate velocities." << std::endl;
  }
  else {
    l_c = clock() - l_c;
    std::cout << "(" << (float) l_c / CLOCKS_PER_SEC << "s)" << std::endl;

    ucvm_prop_t *l_propPtr;

    //! Move to VM Nodes Array
    for( l_pid = 0; l_pid < l_wrkRg.m_numPrvt; l_pid++ ) {
      switch( l_aCfg.m_ucvmType ) {
        case 0:   //! Crustal
          l_propPtr = &(l_ucvmProps[l_pid+l_pOfs].crust);
          break;
        case 1:   //! Geotechnical layer
          l_propPtr = &(l_ucvmProps[l_pid+l_pOfs].gtl);
          break;
        case 2:   //! Combination
          l_propPtr = &(l_ucvmProps[l_pid+l_pOfs].cmb);
          break;
        default:
          l_propPtr = &(l_ucvmProps[l_pid+l_pOfs].cmb);
      }

      l_vModelNodes.m_vmList[l_pid+l_pOfs].m_data[0] = l_propPtr->vp;
      l_vModelNodes.m_vmList[l_pid+l_pOfs].m_data[1] = l_propPtr->vs;
      l_vModelNodes.m_vmList[l_pid+l_pOfs].m_data[2] = l_propPtr->rho;
    }

    delete[] l_ucvmPoints;
    delete[] l_ucvmProps;

    edge_v::vm::Utility::writeVMNodes( l_vModelNodes, l_aCfg, l_msh );

    //! Phase-2: Elements
    edge_v::vm::Utility::vmodel l_vModelElmts;
    edge_v::vm::Utility::vmElmtInit( l_vModelElmts, l_msh );

    edge_v::vm::Utility::workerInit( l_wrkRg,
                                      l_msh.m_numElmts,
                                     l_aCfg.m_projMesh,
                                     l_aCfg.m_projVel );

    edge_v::vm::Utility::vmodel * const l_pVModelNodes = &l_vModelNodes;

    //! Averaging and Clipping
    real l_vp, l_vs, l_rho, l_vpVsRatio, l_lam, l_mu;
    moab::EntityID                    l_eEntId;
    moab::EntityHandle                l_eHandle;
    std::vector< moab::EntityHandle > l_eVertices;

    const int_v l_eOfs = l_wrkRg.m_workerTid * l_wrkRg.m_workSize;
    unsigned int l_pIdx;

    for( int_v l_eid = 0; l_eid < l_wrkRg.m_numPrvt; l_eid++ ) {
      l_eEntId  = l_eid + l_eOfs + 1;
      l_rval    = l_msh.m_intf->handle_from_id( moab::MBTET, l_eEntId, l_eHandle );
      assert( l_rval == moab::MB_SUCCESS );
      l_eVertices.clear();

      l_rval = l_msh.m_intf->get_adjacencies( &l_eHandle, 1, 0, false, l_eVertices );
      assert( l_rval == moab::MB_SUCCESS );
      assert( l_eVertices.size() == ELMTTYPE );

      l_vp  = 0.0;
      l_vs  = 0.0;
      l_rho = 0.0;

      for( int_v l_vId = 0; l_vId < 4; l_vId++ ) {
        moab::EntityID l_pEntId = l_msh.m_intf->id_from_handle( l_eVertices[l_vId] );
        l_pIdx = l_pEntId - 1;

        l_vp  += l_pVModelNodes->m_vmList[l_pIdx].m_data[0];
        l_vs  += l_pVModelNodes->m_vmList[l_pIdx].m_data[1];
        l_rho += l_pVModelNodes->m_vmList[l_pIdx].m_data[2];
      }

      l_vp  /= 4.0;
      l_vs  /= 4.0;
      l_rho /= 4.0;

      //! Second Stop
      if( l_vs < l_aCfg.m_minVs ) {
        l_vpVsRatio = l_vp / l_vs;
        l_vs        = l_aCfg.m_minVs;
        l_vp        = l_aCfg.m_minVs * l_vpVsRatio;
      }

      //! Third Stop
      l_mu = l_rho * l_vs * l_vs;
      if( l_vp > (l_vs * l_aCfg.m_maxVpVsRatio) ) {
        l_lam = l_rho * l_vs * l_vs * l_aCfg.m_maxVpVsRatio * l_aCfg.m_maxVpVsRatio
                - 2.0 * l_mu;
      } else {
        l_lam = l_rho * l_vp * l_vp - 2.0 * l_mu;
      }

      if( l_lam < 0 ) {
        if( l_vs < l_aCfg.m_minVs ) {
          l_vp = 2.45 * l_vs;
        } else if( l_vs < l_aCfg.m_minVs2 ) {
          l_vp = 2 * l_vs;
        } else {
          l_vp = 1.87 * l_vs;
        }

        l_lam = l_rho * l_vp * l_vp;
      }

      l_vModelElmts.m_vmList[l_eid+l_eOfs].m_data[0] = l_lam;
      l_vModelElmts.m_vmList[l_eid+l_eOfs].m_data[1] = l_mu;
      l_vModelElmts.m_vmList[l_eid+l_eOfs].m_data[2] = l_rho;
    } //! Tet Elements Loop Over

    edge_v::vm::Utility::vmNodeFinalize( l_vModelNodes );

    edge_v::vm::Utility::writeVMElmts( l_vModelElmts, l_aCfg, l_msh );

    edge_v::vm::Utility::writeVMTags(  l_vModelElmts, l_aCfg, l_msh );

    edge_v::vm::Utility::vmElmtFinalize( l_vModelElmts );
  }

  faultAntn( l_aCfg, l_msh );
  if( l_ucvmStatus != 0 ) {
    l_msh.m_intf->write_mesh( l_aCfg.m_h5mFn.c_str() );
  }

  edge_v::vm::Utility::meshFinalize( l_msh );

  //! Finish time
  l_te = clock() - l_te;
  std::cout << "Time taken: " << (float) l_te / CLOCKS_PER_SEC << "s\n";

  return EXIT_SUCCESS;
}
