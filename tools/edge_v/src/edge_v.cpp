/**
 * @file This file is part of EDGE.
 *
 * @author Junyi Qiu (juq005 AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2017-2019, Regents of the University of California
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
 * This is the main file of Edge-V.
 **/

#include "vm_utility.h"

int main( int argc, char **argv ) {
  if( argc != 3 ) {
    std::cerr << "Usage : " << argv[0] << " -f config_file.log" << std::endl;
    exit( -1 );
  } else if( (argv == NULL) || (argv[1][0] != '-') || (argv[1][1] != 'f') ) {
    std::cerr << "Usage : " << argv[0] << " -f config_file.log" << std::endl;
    exit( -1 );
  }

  antn_cfg aCfg;
  std::string configFile = std::string( argv[2] );
  antnInit( aCfg, configFile );
  ucvmInit( aCfg );

  moab_mesh mMsh;
  meshInit( mMsh, aCfg );


  //! Phase-1: Nodes

  vmodel vModelNodes;
  vmNodeInit( vModelNodes, mMsh );

  ucvm_point_t *ucvmPoints  = new ucvm_point_t[mMsh.num_nodes];
  ucvm_data_t *ucvmProps    = new ucvm_data_t[mMsh.num_nodes];

  #pragma omp parallel
  {
    worker_reg wrkRg;
    workerInit( wrkRg, mMsh.num_nodes );

    moab::EntityID      pEntId;
    moab::EntityHandle  pHandle;
    moab::ErrorCode     rval;

    const int_v p_ofs = wrkRg.worker_tid * wrkRg.work_size;
    for( int_v pid = 0; pid < wrkRg.num_prvt; pid++ ) {
      pEntId = pid + p_ofs + 1;
      rval = mMsh.intf->handle_from_id( moab::MBVERTEX, pEntId, pHandle );
      assert( rval == moab::MB_SUCCESS );

      rval = mMsh.intf->get_coords( &pHandle, 1, ucvmPoints[pid+p_ofs].coord );
      assert( rval == moab::MB_SUCCESS );
    }

    //! Proj4 Transform: UTM->Long,Lat,Elv
    double *xPtr  = &(ucvmPoints[p_ofs].coord[0]);
    double *yPtr  = &(ucvmPoints[p_ofs].coord[1]);
    double *zPtr  = &(ucvmPoints[p_ofs].coord[2]);
    int pntOfs    = sizeof( ucvm_point_t ) / sizeof( double );
    int pjstatus  = pj_transform( wrkRg.pj_utm, wrkRg.pj_geo, wrkRg.num_prvt,
                                  pntOfs, xPtr, yPtr, zPtr );

    //! Apply Rad to Degree
    for( int_v pid = 0; pid < wrkRg.num_prvt; pid++ ) {
      ucvmPoints[pid+p_ofs].coord[0] *= RAD_TO_DEG;
      ucvmPoints[pid+p_ofs].coord[1] *= RAD_TO_DEG;
      //! (Raj): Depth (m) does not need transformation from rad to deg
      // ucvmPoints[pid+p_ofs].coord[2] *= RAD_TO_DEG;
    }

    //! UCVM Query
#pragma omp master
{
    std::cout << "UCVM Query ... ";
    std::cout.flush();
    int ucvmStatus = ucvm_query( mMsh.num_nodes, ucvmPoints, ucvmProps );
    if( ucvmStatus != 0 ) {
      std::cout << "Failed." << std::endl;
      std::cerr << "Error: cannot complete UCVM query." << std::endl;
    }
    std::cout << "Done!" << std::endl;
}
#pragma omp barrier

    ucvm_prop_t *propPtr;

    //! Move to VM Nodes Array
    for( int_v pid = 0; pid < wrkRg.num_prvt; pid++ ) {
      switch( aCfg.ucvm_type ) {
        case 0:   //! Crustal
          propPtr = &(ucvmProps[pid+p_ofs].crust);
          break;
        case 1:   //! Geotechnical layer
          propPtr = &(ucvmProps[pid+p_ofs].gtl);
          break;
        case 2:   //! Combination
          propPtr = &(ucvmProps[pid+p_ofs].cmb);
          break;
        default:
          propPtr = &(ucvmProps[pid+p_ofs].cmb);
      }

      vModelNodes.vm_list[pid+p_ofs].data[0] = propPtr->vp;
      vModelNodes.vm_list[pid+p_ofs].data[1] = propPtr->vs;
      vModelNodes.vm_list[pid+p_ofs].data[2] = propPtr->rho;
    }
  } //! Exit Parallel Region

  delete[] ucvmPoints;
  delete[] ucvmProps;

  writeVMNodes( vModelNodes, aCfg, mMsh );


  //! Phase-2: Elements

  vmodel vModelElmts;
  vmElmtInit( vModelElmts, mMsh );

  #pragma omp parallel
  {
    worker_reg wrkRg;
    workerInit( wrkRg, mMsh.num_elmts );

    const real c_min_vs           = aCfg.min_vs;
    const real c_min_vs2          = aCfg.min_vs2;
    const real c_max_vp_vs_ratio  = aCfg.max_vp_vs_ratio;

    vmodel * const pVModelNodes = &vModelNodes;

    //! Averaging and Clipping
    real l_vp, l_vs, l_rho, l_vp_vs_ratio, l_lam, l_mu;
    moab::EntityID                    eEntId;
    moab::EntityHandle                eHandle;
    std::vector< moab::EntityHandle > eVertices;
    moab::ErrorCode                   rval;

    const int_v e_ofs = wrkRg.worker_tid * wrkRg.work_size;
    for( int_v eid = 0; eid < wrkRg.num_prvt; eid++ ) {
      eEntId = eid + e_ofs + 1;
      rval = mMsh.intf->handle_from_id( moab::MBTET, eEntId, eHandle );
      assert( rval == moab::MB_SUCCESS );
      eVertices.clear();

      rval = mMsh.intf->get_adjacencies( &eHandle, 1, 0, false, eVertices );
      assert( rval == moab::MB_SUCCESS );
      assert( eVertices.size() == ELMTTYPE );

      l_vp  = 0;
      l_vs  = 0;
      l_rho = 0;

      for( int vid = 0; vid < 4; vid++ ) {
        moab::EntityID pEntId = mMsh.intf->id_from_handle( eVertices[vid] );
        unsigned int pidx = pEntId - 1;

        l_vp  += pVModelNodes->vm_list[pidx].data[0];
        l_vs  += pVModelNodes->vm_list[pidx].data[1];
        l_rho += pVModelNodes->vm_list[pidx].data[2];
      }

      l_vp  /= 4.0;
      l_vs  /= 4.0;
      l_rho /= 4.0;

      //! Second Stop
      if( l_vs < c_min_vs ) {
        l_vp_vs_ratio = l_vp / l_vs;
        l_vs          = aCfg.min_vs;
        l_vp          = aCfg.min_vs * l_vp_vs_ratio;
      }

      //! Third Stop
      l_mu = l_rho * l_vs * l_vs;
      if( l_vp > l_vs * c_max_vp_vs_ratio ) {
        l_lam = l_rho * l_vs * l_vs * c_max_vp_vs_ratio * c_max_vp_vs_ratio - 2
                * l_mu;
      } else {
        l_lam = l_rho * l_vp * l_vp - 2 * l_mu;
      }

      if( l_lam < 0 ) {
        if( l_vs < c_min_vs ) {
          l_vp = 2.45 * l_vs;
        } else if( l_vs * c_min_vs2 ) {
          l_vp = 2 * l_vs;
        } else {
          l_vp = 1.87 * l_vs;
        }

        l_lam = l_rho * l_vp * l_vp;
      }

      vModelElmts.vm_list[eid+e_ofs].data[0] = l_lam;
      vModelElmts.vm_list[eid+e_ofs].data[1] = l_mu;
      vModelElmts.vm_list[eid+e_ofs].data[2] = l_rho;
    } //! Tet Elements Loop Over
  } //! Exit Parallel Region

  vmNodeFinalize( vModelNodes );

  fault_antn( mMsh, aCfg );

  writeVMElmts( vModelElmts, aCfg, mMsh );
  writeVMTags(  vModelElmts, aCfg, mMsh );

  vmElmtFinalize( vModelElmts );

  meshFinalize( mMsh );

  return 0;
}
