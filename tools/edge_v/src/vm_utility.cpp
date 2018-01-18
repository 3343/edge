/**
 * @file This file is part of EDGE.
 *
 * @author Junyi Qiu (juq005 AT ucsd.edu)
 * @author Rajdeep Konwar (rkonwar AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2017, Regents of the University of California
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

#include "vm_utility.h"

int antnInit( antn_cfg & a_cfg, const std::string &cfg_f ) {
  a_cfg.antn_cfg_fn = cfg_f;

  a_cfg.min_vp          = 0;
  a_cfg.min_vs          = 0;
  a_cfg.min_vs2         = 0;
  a_cfg.max_vp_vs_ratio = 0;

  a_cfg.hypoc.lon = 0;
  a_cfg.hypoc.lat = 0;

  std::cout << "Reading Annotation Config File: " << cfg_f << " ... ";
  std::cout.flush();

  std::ifstream iMshFs( cfg_f.c_str(), std::ios::in );
  if( !iMshFs.is_open() ) {
    std::cout << "Failed." << std::endl;
    std::cerr << "Error: cannot open the config file." << std::endl;
    exit( -1 );
  }

  std::string lineBuf;
  while( getline( iMshFs, lineBuf ) ) {
    int i = -1, j;

    while( (++i < lineBuf.length()) && (lineBuf[i] == ' ') );

    if( (i >= lineBuf.length()) || (lineBuf[i] == '#') )
      continue;

    j = i - 1;

    while( (++j < lineBuf.length()) && (lineBuf[j] != '=') );

    if( j >= lineBuf.length() )
      continue;

    std::string varName   = lineBuf.substr( i, j - i );
    std::string varValue  = lineBuf.substr( j + 1 );

    if( varName.compare( "ucvm_config" ) == 0 )
      a_cfg.ucvm_cfg_fn     = varValue;
    else if( varName.compare( "ucvm_model_list" ) == 0 )
      a_cfg.ucvm_model_list = varValue;
    else if( varName.compare( "mesh_file" ) == 0 )
      a_cfg.mesh_fn         = varValue;
    else if( varName.compare( "node_vm_file" ) == 0 )
      a_cfg.vm_node_fn      = varValue;
    else if( varName.compare( "elmt_vm_file" ) == 0 )
      a_cfg.vm_elmt_fn      = varValue;
    else if( varName.compare( "h5m_file" ) == 0 )
      a_cfg.h5m_fn          = varValue;
    else
      std::cout << "\nUnknown setting (" << varName << "). Ignored." << std::endl;
  }

  a_cfg.ucvm_cmode  = UCVMCMODE;
  a_cfg.ucvm_type   = UCVMTYPE;

  a_cfg.min_vp = (a_cfg.min_vp == 0) ? MIN_VP : a_cfg.min_vp;
  a_cfg.min_vs = (a_cfg.min_vs == 0) ? MIN_VS : a_cfg.min_vs;
  a_cfg.min_vs2 = (a_cfg.min_vs2 == 0) ? MIN_VS2 : a_cfg.min_vs2;
  a_cfg.max_vp_vs_ratio = (a_cfg.max_vp_vs_ratio == 0) ? \
                          MAX_VP_VS_RATIO : a_cfg.max_vp_vs_ratio;

  a_cfg.hypoc.lon = CENTERICLON;
  a_cfg.hypoc.lat = CENTERICLAT;

  a_cfg.elmt_type = ELMTTYPE;

  std::cout << "Done!" << std::endl;

  return 0;
}

int ucvmInit( const antn_cfg &a_cfg ) {
  ucvm_init( a_cfg.ucvm_cfg_fn.c_str() );
  ucvm_add_model_list( a_cfg.ucvm_model_list.c_str() );
  ucvm_setparam( UCVM_PARAM_QUERY_MODE, a_cfg.ucvm_cmode );

  return 0;
}

int meshInit( moab_mesh &m_mesh, antn_cfg &a_cfg ) {
  if( m_mesh.intf != nullptr ) {
    std::cout << "Failed." << std::endl;
    std::cerr << "Error: The mesh has been initialized." << std::endl;
    exit( -1 );
  }

  m_mesh.intf = new moab::Core;
  moab::Interface *iface = m_mesh.intf;

  std::cout << "Reading Mesh File: " << a_cfg.mesh_fn << " ... " << std::endl;

  //! Load the mesh from msh file
  moab::ErrorCode rval = iface->load_mesh( a_cfg.mesh_fn.c_str() );
  if( rval != moab::MB_SUCCESS ) {
    std::cout << "Failed." << std::endl;
    std::cerr << "Error: cannot open the mesh file." << std::endl;
    exit( -1 );
  }
  assert( rval == moab::MB_SUCCESS );

  //! Get Nodes && Tets Statistics
  int numNodes, numElmts;
  moab::Range verts, elems;

  rval = iface->get_number_entities_by_type( 0, moab::MBVERTEX, numNodes );
  assert( rval == moab::MB_SUCCESS );

  rval = iface->get_number_entities_by_type( 0, moab::MBTET, numElmts );
  assert( rval == moab::MB_SUCCESS );

  rval = iface->get_entities_by_type( 0, moab::MBVERTEX, verts );
  assert( rval == moab::MB_SUCCESS );
  assert( verts.size() == numNodes );

  rval = iface->get_entities_by_type( 0, moab::MBTET, elems );
  assert( rval == moab::MB_SUCCESS );
  assert( elems.size() == numElmts );

  m_mesh.num_nodes = numNodes;
  m_mesh.num_elmts = numElmts;
  std::cout << " | Number of vertices is " << m_mesh.num_nodes << std::endl;
  std::cout << " | Number of elements is " << m_mesh.num_elmts << std::endl;

  std::cout << "Done!" << std::endl;

  return 0;
}

int meshFinalize( moab_mesh &m_mesh ) {
  if( m_mesh.intf != nullptr ) {
    delete m_mesh.intf;
  } else {
    std::cout << "Failed." << std::endl;
    std::cerr << "Error: mesh object corrupted.";
    exit( -1 );
  }

  return 0;
}

int vmNodeInit( vmodel &vm_nodes, const moab_mesh &m_mesh ) {
  vm_nodes.vm_list = new vm_datum[m_mesh.num_nodes];

  return 0;
}

int vmNodeFinalize( vmodel &vm_nodes ) {
  if( vm_nodes.vm_list != nullptr )
    delete[] vm_nodes.vm_list;

  return 0;
}

int vmElmtInit( vmodel &vm_elmts, const moab_mesh &m_mesh ) {
  vm_elmts.vm_list = new vm_datum[m_mesh.num_elmts];

  return 0;
}

int vmElmtFinalize( vmodel &vm_elmts ) {
  if( vm_elmts.vm_list != nullptr )
    delete[] vm_elmts.vm_list;

  return 0;
}

int workerInit( worker_reg &wrk_rg, int_v totalNum ) {
  int_v tid         = omp_get_thread_num();
  int_v numThrds    = omp_get_num_threads();
  wrk_rg.worker_tid = tid;

  int_v workSize    = (totalNum + numThrds - 1) / numThrds;
  int_v nPrivate    = std::min( workSize * (tid + 1), totalNum ) - workSize * tid;
  wrk_rg.work_size  = workSize;
  wrk_rg.num_prvt   = nPrivate;

  std::string pjInitParams = "+proj=tmerc +units=m +axis=enu +no_defs \
                          +datum=WGS84 +k=0.9996 +lon_0=-117.916 +lat_0=33.933";
  wrk_rg.pj_utm = pj_init_plus( pjInitParams.c_str() );
  wrk_rg.pj_geo = pj_init_plus( "+proj=latlong +datum=WGS84" );

  return 0;
}

int writeVMNodes( vmodel &vm_nodes, const antn_cfg &a_cfg,
                  const moab_mesh &m_mesh ) {
  std::cout << "Write Velocity Model: " << a_cfg.vm_node_fn << " ... ";
  std::cout.flush();

  std::ofstream oVmNodeFs( a_cfg.vm_node_fn.c_str(), std::ios::out );
  if( !oVmNodeFs.is_open() ) {
    std::cout << "Failed" << std::endl;
    std::cerr << "Error: cannot generate the velocity model for nodes."
              << std::endl;
    exit( -1 );
  }
  
  //! Write down the headers
  oVmNodeFs << "$UcvmModel\n";
  oVmNodeFs << a_cfg.ucvm_model_list << std::endl;
  oVmNodeFs << "$EndUcvmModel\n";

  //! Write down the velocity model data
  int_v numNodes = m_mesh.num_nodes;

  oVmNodeFs << "$NodesVelocityModel\n";
  oVmNodeFs << numNodes << std::endl;
  for( int_v pid = 0; pid < numNodes; pid++ ) {
    real data0 = vm_nodes.vm_list[pid].data[0];
    real data1 = vm_nodes.vm_list[pid].data[1];
    real data2 = vm_nodes.vm_list[pid].data[2];

    oVmNodeFs << (pid + 1) << " " << data0 << " " << data1 << " " << data2
              << std::endl;
    oVmNodeFs.flush();
  }
  oVmNodeFs << "$NodesVelocityModel\n";

  oVmNodeFs.close();
  std::cout << "Done!" << std::endl;

  return 0;
}

int writeVMElmts( vmodel &vm_elmts, const antn_cfg &a_cfg,
                  const moab_mesh &m_mesh ) {
  std::cout << "Writing Velocity Model: " << a_cfg.vm_elmt_fn << " ... ";
  std::cout.flush();

  std::ofstream oVmElmtFs( a_cfg.vm_elmt_fn.c_str(), std::ios::out );
  if( !oVmElmtFs.is_open() ) {
    std::cout << "Failed" << std::endl;
    std::cerr << "Error: cannot generate the velocity model for elements."
              << std::endl;
    exit( -1 );
  }

  //! Write down the headers
  oVmElmtFs << "$UcvmModel\n";
  oVmElmtFs << a_cfg.ucvm_model_list << std::endl;
  oVmElmtFs << "$EndUcvmModel\n";

  //! Write down the velocity model data
  int_v numElmts = m_mesh.num_elmts;

  oVmElmtFs << "$ElementsVelocityModel\n";
  oVmElmtFs << numElmts << std::endl;
  for( int_v eid = 0; eid < numElmts; eid++ ) {
    real data0 = vm_elmts.vm_list[eid].data[0];
    real data1 = vm_elmts.vm_list[eid].data[1];
    real data2 = vm_elmts.vm_list[eid].data[2];

    oVmElmtFs << (eid + 1) << " " << data0 << " " << data1 << " " << data2
              << std::endl;
    oVmElmtFs.flush();
  }
  oVmElmtFs << "$ElementsVelocityModel\n";

  oVmElmtFs.close();
  std::cout << "Done!" << std::endl;

  return 0;
}

int writeVMTags( vmodel &vm_elmts, const antn_cfg &a_cfg,
                 const moab_mesh &m_mesh ) {
  moab::Interface *iface = m_mesh.intf;
  if( iface == nullptr ) {
    std::cout << "Failed." << std::endl;
    std::cerr << "Error: cannot operate the mesh." << std::endl;
    exit( -1 );
  }

  std::cout << "Writing Annotated Mesh File: " << a_cfg.h5m_fn << " ... ";
  std::cout.flush();

  moab::Range elems;
  moab::ErrorCode rval = iface->get_entities_by_type( 0, moab::MBTET, elems );
  assert( rval == moab::MB_SUCCESS );

  //! Create Tags
  moab::Tag tag_lambda, tag_mu, tag_rho;
  rval = iface->tag_get_handle( "LAMBDA", 4, moab::MB_TYPE_OPAQUE, tag_lambda,
                    moab::MB_TAG_CREAT|moab::MB_TAG_DENSE|moab::MB_TAG_BYTES );
  assert( rval == moab::MB_SUCCESS );

  rval = iface->tag_get_handle( "MU", 4, moab::MB_TYPE_OPAQUE, tag_mu,
                    moab::MB_TAG_CREAT|moab::MB_TAG_DENSE|moab::MB_TAG_BYTES );
  assert( rval == moab::MB_SUCCESS );

  rval = iface->tag_get_handle( "RHO", 4, moab::MB_TYPE_OPAQUE, tag_rho,
                    moab::MB_TAG_CREAT|moab::MB_TAG_DENSE|moab::MB_TAG_BYTES );
  assert( rval == moab::MB_SUCCESS );

  //! Set Tags
  int_v eid;
  int_v numElmts = m_mesh.num_elmts;
  real l_lambda, l_mu, l_rho;

  for( moab::Range::iterator rit = elems.begin(); rit != elems.end(); ++rit ) {
    moab::EntityID entId = iface->id_from_handle( *rit );

    eid = entId - 1;
    assert( eid < numElmts );

    l_lambda  = vm_elmts.vm_list[eid].data[0];
    l_mu      = vm_elmts.vm_list[eid].data[1];
    l_rho     = vm_elmts.vm_list[eid].data[2];


    iface->tag_set_data( tag_lambda, &(*rit), 1, &l_lambda );
    assert( rval == moab::MB_SUCCESS );

    rval = iface->tag_set_data( tag_mu, &(*rit), 1, &l_mu );
    assert( rval == moab::MB_SUCCESS );

    rval = iface->tag_set_data( tag_rho, &(*rit), 1, &l_rho );
    assert( rval == moab::MB_SUCCESS );
  }

  //! Write to H5M File
  rval = iface->write_file( a_cfg.h5m_fn.c_str(), "H5M" );
  assert( rval == moab::MB_SUCCESS );

  std::cout << "Done!" << std::endl;

  return 0;
}
