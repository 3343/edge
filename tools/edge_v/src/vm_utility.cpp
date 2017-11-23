/**
 * @file This file is part of EDGE.
 *
 * @author Junyi Qiu (juq005 AT ucsd.edu)
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

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <assert.h>
#include <omp.h>

#include "vm_utility.h"

extern "C" {
#include "ucvm.h"
}
#include "proj_api.h"
#include "moab/Core.hpp"

using namespace moab;
using namespace std;

#define UCVMCMODE UCVM_COORD_GEO_ELEV

#define MIN_VP 1500.0
#define MIN_VS 500.0
#define MIN_VS2 1200.0
#define MAX_VP_VS_RATIO 3.0

#define CENTERICLON -117.916
#define CENTERICLAT 33.933

#define ELMTTYPE 4 //TODO

#define min(X, Y) (((X) < (Y)) ? (X) : (Y))
#define max(X, Y) (((X) < (Y)) ? (X) : (Y))



int antnInit(antn_cfg & a_cfg, const string &cfg_f) {
  a_cfg.antn_cfg_fn = cfg_f;

  a_cfg.min_vp = 0;
  a_cfg.min_vs = 0;
  a_cfg.min_vs2 = 0;
  a_cfg.max_vp_vs_ratio = 0;

  a_cfg.hypoc.lon = 0;
  a_cfg.hypoc.lat = 0;

  cout << "Reading Annotation Config File: " << cfg_f << " ... ";
  cout.flush();

  ifstream iMshFs(cfg_f.c_str(), ios::in);
  if (iMshFs == NULL) {
    cout << "Failed." << endl;
    cerr << "Error: cannot open the annotation config file." << endl;
    exit(-1);
  }

  string lineBuf;
  while (getline(iMshFs, lineBuf)) {
    int i = -1, j;
    while ((++i < lineBuf.length()) && (lineBuf[i] == ' ')) ;
    if ((i >= lineBuf.length()) || (lineBuf[i] == '#')) continue;

    j = i-1;
    while ((++j < lineBuf.length()) && (lineBuf[j] != '=')) ;
    if (j >= lineBuf.length()) continue;

    string varName = lineBuf.substr(i, j-i);
    string varValue = lineBuf.substr(j+1);
    if (varName.compare("ucvm_config") == 0) a_cfg.ucvm_cfg_fn = varValue;
    else if (varName.compare("ucvm_model_list") == 0) a_cfg.ucvm_model_list = varValue;
    else if (varName.compare("mesh_file") == 0) a_cfg.mesh_fn = varValue;
    else if (varName.compare("node_vm_file") == 0) a_cfg.vm_node_fn = varValue;
    else if (varName.compare("elmt_vm_file") == 0) a_cfg.vm_elmt_fn = varValue;
    else if (varName.compare("h5m_file") == 0) a_cfg.h5m_fn = varValue;
    else 
      cout << "\nUnknown setting (" << varName << "). Ignored." << endl;
  }


  a_cfg.ucvm_cmode=UCVMCMODE;
  a_cfg.ucvm_type = 2; //TODO

  a_cfg.min_vp = (a_cfg.min_vp == 0) ? MIN_VP : a_cfg.min_vp;
  a_cfg.min_vs = (a_cfg.min_vs == 0) ? MIN_VS : a_cfg.min_vs;
  a_cfg.min_vs2 = (a_cfg.min_vs2 == 0) ? MIN_VS2 : a_cfg.min_vs2;
  a_cfg.max_vp_vs_ratio = (a_cfg.max_vp_vs_ratio == 0) ? MAX_VP_VS_RATIO : a_cfg.max_vp_vs_ratio;

  a_cfg.hypoc.lon = CENTERICLON;
  a_cfg.hypoc.lat = CENTERICLAT;

  a_cfg.elmt_type = ELMTTYPE;


  cout << "Done!" << endl;

  return 0;
}

int ucvmInit(const antn_cfg &a_cfg) {
  ucvm_init(a_cfg.ucvm_cfg_fn.c_str());
  ucvm_add_model_list(a_cfg.ucvm_model_list.c_str());
  ucvm_setparam(UCVM_PARAM_QUERY_MODE, a_cfg.ucvm_cmode);

  return 0;
}

int meshInit(moab_mesh &m_mesh, antn_cfg &a_cfg){
  if (m_mesh.intf != NULL) {
    cout << "Failed." << endl;
    cerr << "Error: The mesh has been initialized." << endl;
    exit(-1);
  }

  m_mesh.intf = new Core;
  Interface *iface = m_mesh.intf;

  cout << "Reading Mesh File: " << a_cfg.mesh_fn << " ... " << endl;

  // Load the mesh from msh file
  ErrorCode rval = iface->load_mesh( a_cfg.mesh_fn.c_str() );
  if (rval != MB_SUCCESS) {
    cout << "Failed." << endl;
    cerr << "Error: cannot open the mesh file." << endl;
    exit(-1);
  }
  assert(rval == MB_SUCCESS);

  // Get Nodes && Tets Statistics
  int numNodes, numElmts;
  Range verts, elems;
  rval = iface->get_number_entities_by_type(0, MBVERTEX, numNodes);
  assert(rval == MB_SUCCESS);
  rval = iface->get_number_entities_by_type(0, MBTET, numElmts);
  assert(rval == MB_SUCCESS);
  rval = iface->get_entities_by_type(0, MBVERTEX, verts);
  assert(rval == MB_SUCCESS);
  assert(verts.size() == numNodes);
  rval = iface->get_entities_by_type(0, MBTET, elems);
  assert(rval == MB_SUCCESS);
  assert(elems.size() == numElmts);

  m_mesh.num_nodes = numNodes;
  m_mesh.num_elmts = numElmts;
  cout << " | Number of vertices is " << m_mesh.num_nodes << endl;
  cout << " | Number of elements is " << m_mesh.num_elmts << endl;

  cout << "Done!" << endl;
  
  return 0;
}

int meshFinalize(moab_mesh &m_mesh){
  if (m_mesh.intf != NULL) {
    delete m_mesh.intf;
  } else {
    cout << "Failed." << endl;
    cerr << "Error: mesh object corrupted.";
    exit(-1);
  }
  
  return 0;
}

int vmNodeInit(vmodel &vm_nodes, const moab_mesh &m_mesh) {
  vm_nodes.vm_list = new vm_datum[m_mesh.num_nodes];

  return 0;
}

int vmNodeFinalize(vmodel &vm_nodes) {
  if (vm_nodes.vm_list != NULL)
    delete[] vm_nodes.vm_list;

  return 0;
}

int vmElmtInit(vmodel &vm_elmts, const moab_mesh &m_mesh) {
  vm_elmts.vm_list = new vm_datum[m_mesh.num_elmts];

  return 0;
}

int vmElmtFinalize(vmodel &vm_elmts) {
  if (vm_elmts.vm_list != NULL)
    delete[] vm_elmts.vm_list;

  return 0;
}

int workerInit(worker_reg &wrk_rg, unsigned int totalNum) {
  unsigned int tid = omp_get_thread_num();
  unsigned int numThrds = omp_get_num_threads();
  wrk_rg.worker_tid = tid;

  unsigned int workSize = (totalNum + numThrds - 1) / numThrds;
  unsigned int nPrivate = min(workSize*(tid+1), totalNum) - workSize*tid;
  wrk_rg.work_size = workSize;
  wrk_rg.num_prvt = nPrivate;

  string pjInitParams = "+proj=tmerc +units=m +axis=enu +no_defs +datum=WGS84 +k=0.9996 +lon_0=-117.916 +lat_0=33.933";
  wrk_rg.pj_utm = pj_init_plus(pjInitParams.c_str());
  wrk_rg.pj_geo = pj_init_plus("+proj=latlong +datum=WGS84");

  return 0;
}

int writeVMNodes(vmodel &vm_nodes, const antn_cfg &a_cfg, const moab_mesh &m_mesh) {
  cout << "Write Velocity Model: " << a_cfg.vm_node_fn << " ... ";
  cout.flush();

  ofstream oVmNodeFs(a_cfg.vm_node_fn.c_str(), ios::out);
  if (oVmNodeFs == NULL) {
    cout << "Failed" << endl;
    cerr << "Error: cannot generate the velocity model for nodes." << endl;
    exit(-1);
  }
  
  // Write down the headers
  oVmNodeFs << "$UcvmModel\n";
  oVmNodeFs << a_cfg.ucvm_model_list << endl;
  oVmNodeFs << "$EndUcvmModel\n";

  // Write down the velocity model data
  unsigned int numNodes = m_mesh.num_nodes;

  oVmNodeFs << "$NodesVelocityModel\n";
  oVmNodeFs << numNodes << endl;
  for (int pid = 0; pid < numNodes; pid++) {
    double data0 = vm_nodes.vm_list[pid].data[0];
    double data1 = vm_nodes.vm_list[pid].data[1];
    double data2 = vm_nodes.vm_list[pid].data[2];

    oVmNodeFs << (pid+1) << " " << data0 << " " << data1 << " " << data2 << endl;
    oVmNodeFs.flush();
  } 
  oVmNodeFs << "$NodesVelocityModel\n";

  oVmNodeFs.close();
  cout << "Done!" << endl;

  return 0;
}

int writeVMElmts(vmodel &vm_elmts, const antn_cfg &a_cfg, const moab_mesh &m_mesh) {
  cout << "Writing Velocity Model: " << a_cfg.vm_elmt_fn << " ... ";
  cout.flush();

  ofstream oVmElmtFs(a_cfg.vm_elmt_fn.c_str(), ios::out);
  if (oVmElmtFs == NULL) {
    cout << "Failed" << endl;
    cerr << "Error: cannot generate the velocity model for elements." << endl;
    exit(-1);
  }
  
  // Write down the headers
  oVmElmtFs << "$UcvmModel\n";
  oVmElmtFs << a_cfg.ucvm_model_list << endl;
  oVmElmtFs << "$EndUcvmModel\n";

  // Write down the velocity model data
  unsigned int numElmts = m_mesh.num_elmts;

  oVmElmtFs << "$ElementsVelocityModel\n";
  oVmElmtFs << numElmts << endl;
  for (int eid = 0; eid < numElmts; eid++) {
    double data0 = vm_elmts.vm_list[eid].data[0];
    double data1 = vm_elmts.vm_list[eid].data[1];
    double data2 = vm_elmts.vm_list[eid].data[2];

    oVmElmtFs << (eid+1) << " " << data0 << " " << data1 << " " << data2 << endl;
    oVmElmtFs.flush();
  } 
  oVmElmtFs << "$ElementsVelocityModel\n";

  oVmElmtFs.close();
  cout << "Done!" << endl;

  return 0;
}

int writeVMTags(vmodel &vm_elmts, const antn_cfg &a_cfg, const moab_mesh &m_mesh) {
  Interface *iface = m_mesh.intf;
  if (iface == NULL) {
    cout << "Failed." << endl;
    cerr << "Error: cannot operate the mesh." << endl;
    exit(-1);
  }

  cout << "Writing Annotated Mesh File: " << a_cfg.h5m_fn << " ... ";
  cout.flush();

  Range elems;
  ErrorCode rval = iface->get_entities_by_type(0, MBTET, elems);
  assert(rval == MB_SUCCESS);

  // Create Tags
  Tag tag_lambda, tag_mu, tag_rho;
  rval = iface->tag_get_handle("LAMBDA", 4, MB_TYPE_OPAQUE, tag_lambda, MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES);
  assert(rval == MB_SUCCESS);
  rval = iface->tag_get_handle("MU", 4, MB_TYPE_OPAQUE, tag_mu, MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES);
  assert(rval == MB_SUCCESS);
  rval = iface->tag_get_handle("RHO", 4, MB_TYPE_OPAQUE, tag_rho, MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES);
  assert(rval == MB_SUCCESS);

  // Set Tags
  unsigned int numElmts = m_mesh.num_elmts;
  for (Range::iterator rit = elems.begin(); rit != elems.end(); rit++) {
    EntityID entId = iface->id_from_handle(*rit);
    
    unsigned long long eid = entId - 1;
    assert(eid < int(numElmts));

    float l_lambda, l_mu, l_rho;
    l_lambda = vm_elmts.vm_list[eid].data[0];
    l_mu = vm_elmts.vm_list[eid].data[1];
    l_rho = vm_elmts.vm_list[eid].data[2];


    iface->tag_set_data(tag_lambda, &(*rit), 1, &l_lambda);
    assert(rval == MB_SUCCESS);
    rval = iface->tag_set_data(tag_mu, &(*rit), 1, &l_mu);
    assert(rval == MB_SUCCESS);
    rval = iface->tag_set_data(tag_rho, &(*rit), 1, &l_rho);
    assert(rval == MB_SUCCESS);
  }

  // Write to H5M File
  rval = iface->write_file(a_cfg.h5m_fn.c_str(), "H5M");
  assert(rval == MB_SUCCESS);

  cout << "Done!" << endl;
}










