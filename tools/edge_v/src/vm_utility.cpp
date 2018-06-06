/**
 * @file This file is part of EDGE.
 *
 * @author Junyi Qiu (juq005 AT ucsd.edu)
 * @author Rajdeep Konwar (rkonwar AT ucsd.edu)
 * @author David Lenz (dlenz AT ucsd.edu)
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
 * This file contains the utility routines of Edge-V.
 **/

#include <fstream>
#include <omp.h>
#include <sstream>
#include <moab/Core.hpp>

#include "vm_utility.h"

void edge_v::vm::Utility::lamePar( double  i_vp,
                                   double  i_vs,
                                   double  i_rho,
                                   double &o_lam,
                                   double &o_mu ) {
  o_mu   = i_vs * i_vs * i_rho;
  o_lam  = i_vp * i_vp * i_rho;
  o_lam -= 2*o_mu;
}

int edge_v::vm::Utility::ucvmInit( const io::Config &i_cfg ) {
  clock_t l_t = clock();
  std::cout << "Initializing UCVM... " << std::flush;

  if( ucvm_init( i_cfg.m_ucvmCfgFn.c_str() ) != UCVM_CODE_SUCCESS ) {
    std::cerr << "failed initializing UCVM" << std::endl;
    return 1;
  };

  if( ucvm_add_model_list( i_cfg.m_ucvmModelList.c_str() ) != UCVM_CODE_SUCCESS ) {
    std::cerr << "failed setting model: " << i_cfg.m_ucvmModelList << std::endl;
    return 1;
  }

  if( ucvm_setparam( UCVM_PARAM_QUERY_MODE, i_cfg.m_ucvmCmode ) != UCVM_CODE_SUCCESS  ) {
    std::cerr << "failed setting param query mode: " << i_cfg.m_ucvmModelList << std::endl;
    return 1;
  }

  std::cout << "Done! ";
  l_t = clock() - l_t;
  std::cout << "(" << (float) l_t / CLOCKS_PER_SEC << "s)" << std::endl;

  return 0;
}

int edge_v::vm::Utility::meshInit( moab_mesh &i_mesh,
              io::Config  &i_cfg ) {
  clock_t l_t = clock();

  if( i_mesh.m_intf != nullptr ) {
    std::cout << "Failed." << std::endl;
    std::cerr << "Error: The mesh has been initialized." << std::endl;
    exit( EXIT_FAILURE );
  }

  i_mesh.m_intf           = new moab::Core;
  moab::Interface *l_face = i_mesh.m_intf;

  std::cout << "Reading mesh file: " << i_cfg.m_meshFn << "... " << std::endl;

  //! Load the mesh from msh file
  moab::ErrorCode l_rval = l_face->load_mesh( i_cfg.m_meshFn.c_str() );
  if( l_rval != moab::MB_SUCCESS ) {
    std::cout << "Failed." << std::endl;
    std::cerr << "Error: cannot open the mesh file." << std::endl;
    exit( EXIT_FAILURE );
  }
  assert( l_rval == moab::MB_SUCCESS );

  //! Get Nodes && Tets Statistics
  int         l_numNodes, l_numElmts;
  moab::Range l_verts,    l_elems;

  l_rval = l_face->get_number_entities_by_type( 0, moab::MBVERTEX, l_numNodes );
  assert( l_rval == moab::MB_SUCCESS );

  if( ELMTTYPE == 3 )
    l_rval = l_face->get_number_entities_by_type( 0, moab::MBTRI, l_numElmts );
  else if( ELMTTYPE == 4 )
    l_rval = l_face->get_number_entities_by_type( 0, moab::MBTET, l_numElmts );
  assert( l_rval == moab::MB_SUCCESS );

  l_rval = l_face->get_entities_by_type( 0, moab::MBVERTEX, l_verts );
  assert( l_rval == moab::MB_SUCCESS );
  assert( l_verts.size() == (size_t) l_numNodes );

  if( ELMTTYPE == 3 )
    l_rval = l_face->get_entities_by_type( 0, moab::MBTRI, l_elems );
  else if( ELMTTYPE == 4 )
    l_rval = l_face->get_entities_by_type( 0, moab::MBTET, l_elems );
  assert( l_rval == moab::MB_SUCCESS );
  assert( l_elems.size() == (size_t) l_numElmts );

  i_mesh.m_numNodes = l_numNodes;
  i_mesh.m_numElmts = l_numElmts;
  std::cout << " | Number of vertices is " << i_mesh.m_numNodes << std::endl;
  std::cout << " | Number of elements is " << i_mesh.m_numElmts << std::endl;

  std::cout << "Done! ";
  l_t = clock() - l_t;
  std::cout << "(" << (float) l_t / CLOCKS_PER_SEC << "s)" << std::endl;

  return 0;
}

int edge_v::vm::Utility::meshFinalize( moab_mesh &i_mesh ) {
  if( i_mesh.m_intf != nullptr ) {
    delete i_mesh.m_intf;
  } else {
    std::cout << "Failed." << std::endl;
    std::cerr << "Error: mesh object corrupted.";
    exit( EXIT_FAILURE );
  }

  return 0;
}

int edge_v::vm::Utility::vmNodeInit(       vmodel    &i_vmNodes,
                                     const moab_mesh &i_mesh ) {
  i_vmNodes.m_vmList = new vm_datum[i_mesh.m_numNodes];

  return 0;
}

int edge_v::vm::Utility::vmNodeFinalize( vmodel &i_vmNodes ) {
  if( i_vmNodes.m_vmList != nullptr )
    delete[] i_vmNodes.m_vmList;

  return 0;
}

int edge_v::vm::Utility::vmElmtInit(       vmodel    &i_vmElmts,
                const moab_mesh &i_mesh ) {
  i_vmElmts.m_vmList = new vm_datum[i_mesh.m_numElmts];

  return 0;
}

int edge_v::vm::Utility::vmElmtFinalize( vmodel &i_vmElmts ) {
  if( i_vmElmts.m_vmList != nullptr )
    delete[] i_vmElmts.m_vmList;

  return 0;
}

int edge_v::vm::Utility::workerInit(       worker_reg  &i_wrkRg,
                                     const int_v       &i_totalNum,
                                     const std::string &i_projMesh,
                                     const std::string &i_projVel   ) {
  int_v l_tid                 = omp_get_thread_num();
  int_v l_numThrds            = omp_get_num_threads();
  i_wrkRg.m_workerTid         = l_tid;

  int_v l_workSize            = (i_totalNum + l_numThrds - 1) / l_numThrds;
  int_v l_nPrivate            = std::min( l_workSize * (l_tid + 1), i_totalNum )
                                - l_workSize * l_tid;
  i_wrkRg.m_workSize          = l_workSize;
  i_wrkRg.m_numPrvt           = l_nPrivate;

  i_wrkRg.m_pjUtm             = pj_init_plus( i_projMesh.c_str() );
  i_wrkRg.m_pjGeo             = pj_init_plus( i_projVel.c_str() );

  return 0;
}

int edge_v::vm::Utility::writeVMNodes(       vmodel    &i_vmNodes,
                                       const io::Config  &i_cfg,
                                       const moab_mesh &i_mesh ) {
  clock_t l_t = clock();

  std::cout << "Writing velocity model: " << i_cfg.m_vmNodeFn << "... "
            << std::flush;

  std::ofstream l_vmNodeFs( i_cfg.m_vmNodeFn.c_str(), std::ios::out );
  if( !l_vmNodeFs.is_open() ) {
    std::cout << "Failed" << std::endl;
    std::cerr << "Error: cannot generate the velocity model for nodes."
              << std::endl;
    exit( EXIT_FAILURE );
  }

  int_v l_bufferSize = 10000;
  std::stringstream l_ss;

  //! Write down the headers
  l_ss << "$UcvmModel" << std::endl;
  l_ss << i_cfg.m_ucvmModelList << std::endl;
  l_ss << "$EndUcvmModel" << std::endl;

  //! Write down the velocity model data
  int_v l_numNodes = i_mesh.m_numNodes;

  l_ss << "$NodesVelocityModel" << std::endl;
  l_ss << l_numNodes << std::endl;

  real l_data0, l_data1, l_data2;
  for( int_v l_pid = 0; l_pid < l_numNodes; l_pid++ ) {
    if( l_pid % l_bufferSize == 0 ) {
      //! Write buffer to file
      l_vmNodeFs << l_ss.rdbuf();
      l_vmNodeFs.flush();

      //! Reset buffer
      l_ss.str( std::string() );
      l_ss.clear();
    }

    l_data0 = i_vmNodes.m_vmList[l_pid].m_data[0];
    l_data1 = i_vmNodes.m_vmList[l_pid].m_data[1];
    l_data2 = i_vmNodes.m_vmList[l_pid].m_data[2];

    l_ss << (l_pid + 1) << " " << l_data0 << " " << l_data1 << " "
         << l_data2 << std::endl;
  }

  l_ss << "$NodesVelocityModel" << std::endl;

  //! Final buffer write to file
  l_vmNodeFs << l_ss.rdbuf();
  l_vmNodeFs.flush();
  l_vmNodeFs.close();

  std::cout << "Done! ";
  l_t = clock() - l_t;
  std::cout << "(" << (float) l_t / CLOCKS_PER_SEC << "s)" << std::endl;

  return 0;
}

int edge_v::vm::Utility::writeVMElmts(       vmodel    &i_vmElmts,
                                       const io::Config  &i_cfg,
                                       const moab_mesh &i_mesh ) {
  clock_t l_t = clock();

  std::cout << "Writing velocity model: " << i_cfg.m_vmElmtFn << "... "
            << std::flush;

  std::ofstream l_vmElmtFs( i_cfg.m_vmElmtFn.c_str(), std::ios::out );
  if( !l_vmElmtFs.is_open() ) {
    std::cout << "Failed" << std::endl;
    std::cerr << "Error: cannot generate the velocity model for elements."
              << std::endl;
    exit( EXIT_FAILURE );
  }

  int_v l_bufferSize = 10000;
  std::stringstream l_ss;

  //! Write down the headers
  l_ss << "$UcvmModel" << std::endl;
  l_ss << i_cfg.m_ucvmModelList << std::endl;
  l_ss << "$EndUcvmModel" << std::endl;

  //! Write down the velocity model data
  int_v l_numElmts = i_mesh.m_numElmts;

  l_ss << "$ElementsVelocityModel" << std::endl;
  l_ss << l_numElmts << std::endl;

  real l_data0, l_data1, l_data2;
  for( int_v l_eid = 0; l_eid < l_numElmts; l_eid++ ) {
    if( l_eid % l_bufferSize == 0 ) {
      //! Write buffer to file
      l_vmElmtFs << l_ss.rdbuf();
      l_vmElmtFs.flush();

      //! Reset buffer
      l_ss.str( std::string() );
      l_ss.clear();
    }

    l_data0 = i_vmElmts.m_vmList[l_eid].m_data[0];
    l_data1 = i_vmElmts.m_vmList[l_eid].m_data[1];
    l_data2 = i_vmElmts.m_vmList[l_eid].m_data[2];

    l_ss << (l_eid + 1) << " " << l_data0 << " " << l_data1 << " "
         << l_data2 << std::endl;
  }

  l_ss << "$ElementsVelocityModel" << std::endl;

  //! Final buffer write to file
  l_vmElmtFs << l_ss.rdbuf();
  l_vmElmtFs.flush();
  l_vmElmtFs.close();

  std::cout << "Done! ";
  l_t = clock() - l_t;
  std::cout << "(" << (float) l_t / CLOCKS_PER_SEC << "s)" << std::endl;

  return 0;
}

int edge_v::vm::Utility::writeVMTags(       vmodel    &i_vmElmts,
                                      const io::Config  &i_cfg,
                                      const moab_mesh &i_mesh ) {
  clock_t l_t = clock();

  moab::Interface *l_face = i_mesh.m_intf;
  if( l_face == nullptr ) {
    std::cout << "Failed." << std::endl;
    std::cerr << "Error: cannot operate the mesh." << std::endl;
    exit( EXIT_FAILURE );
  }

  std::cout << "Writing annotated mesh file: " << i_cfg.m_h5mFn << "... "
            << std::flush;

  moab::Range l_elems;
  moab::ErrorCode l_rval = l_face->get_entities_by_type( 0, moab::MBTET, l_elems );
  assert( l_rval == moab::MB_SUCCESS );

  //! Create Tags
  moab::Tag l_tagLambda, l_tagMu, l_tagRho;
  l_rval = l_face->tag_get_handle( "LAMBDA", 1, moab::MB_TYPE_DOUBLE, l_tagLambda,
                                   moab::MB_TAG_CREAT | moab::MB_TAG_DENSE );
  assert( l_rval == moab::MB_SUCCESS );

  l_rval = l_face->tag_get_handle( "MU", 1, moab::MB_TYPE_DOUBLE, l_tagMu,
                                   moab::MB_TAG_CREAT | moab::MB_TAG_DENSE );
  assert( l_rval == moab::MB_SUCCESS );

  l_rval = l_face->tag_get_handle( "RHO", 1, moab::MB_TYPE_DOUBLE, l_tagRho,
                                   moab::MB_TAG_CREAT | moab::MB_TAG_DENSE );
  assert( l_rval == moab::MB_SUCCESS );

  //! Set Tags
  int_v l_eid;
  int_v l_numElmts = i_mesh.m_numElmts;
  real  l_lambda, l_mu, l_rho;

  for( moab::Range::iterator l_rit = l_elems.begin(); l_rit != l_elems.end(); ++l_rit ) {
    moab::EntityID l_entId = l_face->id_from_handle( *l_rit );

    l_eid = l_entId - 1;
    assert( l_eid < l_numElmts );

    l_lambda  = i_vmElmts.m_vmList[l_eid].m_data[0];
    l_mu      = i_vmElmts.m_vmList[l_eid].m_data[1];
    l_rho     = i_vmElmts.m_vmList[l_eid].m_data[2];


    l_face->tag_set_data( l_tagLambda, &(*l_rit), 1, &l_lambda );
    assert( l_rval == moab::MB_SUCCESS );

    l_rval = l_face->tag_set_data( l_tagMu, &(*l_rit), 1, &l_mu );
    assert( l_rval == moab::MB_SUCCESS );

    l_rval = l_face->tag_set_data( l_tagRho, &(*l_rit), 1, &l_rho );
    assert( l_rval == moab::MB_SUCCESS );
  }

  //! Write to H5M File
  l_rval = l_face->write_file( i_cfg.m_h5mFn.c_str(), "H5M" );
  assert( l_rval == moab::MB_SUCCESS );

  std::cout << "Done! ";
  l_t = clock() - l_t;
  std::cout << "(" << (float) l_t / CLOCKS_PER_SEC << "s)" << std::endl;

  return 0;
}

int edge_v::vm::Utility::posInit(       posModel  &i_posModel,
                                  const moab_mesh &i_msh ) {
  i_posModel.m_posList = new pos_datum[i_msh.m_numElmts];

  return 0;
}

int edge_v::vm::Utility::posFinalize( posModel &i_posModel ) {
  if( i_posModel.m_posList != nullptr )
    delete[] i_posModel.m_posList;

  return 0;
}

int edge_v::vm::Utility::writePos( const posModel  &i_posModel,
                                   const io::Config  &i_cfg,
                                   const moab_mesh &i_msh ) {
  clock_t l_t = clock();

  std::cout << "Writing pos file: " << i_cfg.m_posFn << "... " << std::flush;

  std::ofstream l_outPos( i_cfg.m_posFn.c_str(), std::ios::out );
  if( !l_outPos.is_open() ) {
    std::cout << "Failed" << std::endl;
    std::cerr << "Error: cannot generate the pos file." << std::endl;
    exit( EXIT_FAILURE );
  }

  int_v l_bufferSize = 10000;
  std::stringstream l_ss;

  l_ss << "/*********************************************************************"
       << std::endl << " *" << std::endl << " *  EDGE generated pos file"
       << std::endl << " *" << std::endl
       << " *********************************************************************/"
       << std::endl << std::endl;

  l_ss << "View \"EDGEpos\" {" << std::endl;

  for( int_v l_eid = 0; l_eid < i_msh.m_numElmts; l_eid++ ) {
    if( l_eid % l_bufferSize == 0 ) {
      //! Write buffer to file
      l_outPos << l_ss.rdbuf();
      l_outPos.flush();

      //! Reset buffer
      l_ss.str( std::string() );
      l_ss.clear();
    }

    if( ELMTTYPE == 3 )
      l_ss << "ST(";
    else if( ELMTTYPE == 4 )
      l_ss << "SS(";

    for( unsigned int l_vid = 0; l_vid < ELMTTYPE; l_vid++ ) {
      l_ss << i_posModel.m_posList[l_eid].m_xyzPts[l_vid].m_x << ","
           << i_posModel.m_posList[l_eid].m_xyzPts[l_vid].m_y << ","
           << i_posModel.m_posList[l_eid].m_xyzPts[l_vid].m_z;

      if( l_vid != ELMTTYPE - 1 )
        l_ss << ",";
    }

    l_ss << "){";

    for( unsigned int l_vid = 0; l_vid < ELMTTYPE; l_vid++ ) {
      l_ss << i_posModel.m_posList[l_eid].m_vs[l_vid];

      if( l_vid != ELMTTYPE - 1 )
        l_ss << ",";
    }

    l_ss << "};" << std::endl;
  }

  l_ss << "};" << std::endl;

  //! Final buffer write to file
  l_outPos << l_ss.rdbuf();
  l_outPos.flush();
  l_outPos.close();

  std::cout << "Done! ";
  l_t = clock() - l_t;
  std::cout << "(" << (float) l_t / CLOCKS_PER_SEC << "s)" << std::endl;

  return 0;
}