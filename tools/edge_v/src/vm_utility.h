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
 * This is the header file for the utility routines of Edge-V.
 **/

#include <iostream>
#include <string>

extern "C" {
#include "ucvm.h"
}

#include "proj_api.h"

#include "moab/Core.hpp"

using namespace std;

#define EDGE_V_OUT (cout << "EDGE-V INFO: ")
#define EDGE_V_ERR (cerr << "EDGE-V ERR : ")


// *** Coordination System ***

typedef struct xyz_point_t {
  double x;
  double y;
  double z;
} xyz_point_t;


typedef struct geo_point_t {
  double lon;
  double lat;
  double dep;
} geo_point_t;

// **************************


// *** Configuration File ***

typedef struct antn_cfg {
  string antn_cfg_fn;

  string ucvm_cfg_fn;
  string ucvm_model_list;
  ucvm_ctype_t ucvm_cmode;
  int ucvm_type;

  double min_vp;
  double min_vs;
  double min_vs2;
  double max_vp_vs_ratio;

  string mesh_fn;
  geo_point_t hypoc;

  unsigned int elmt_type;

  string vm_node_fn;
  string vm_elmt_fn;
  string h5m_fn;

  unsigned int parallel_mode;
  unsigned int num_worker;

} antn_cfg;


int antnInit(antn_cfg &, const string &);
int ucvmInit(const antn_cfg &);

// // ***************************


// // *** Coordination System ***

// int xyz2geo(const xyz_point_t &, geo_point_t &);
// int geo2xyz(const geo_point_t &, xyz_point_t &);

// // ***************************


// *** Parallel Module ***

typedef struct worker_reg
{
  unsigned int worker_tid;

  projPJ pj_utm;
  projPJ pj_geo;

  unsigned int vm_type;
  unsigned int work_size;
  unsigned int num_prvt;


} worker_reg;

int workerInit(worker_reg &, unsigned int);

// int pjUtmInit(projPJ **);
// int pjUtmFinalize(projPJ *);

// int pjGeoInit(projPJ **);
// int pjGeoFinalize(projPJ *);

// ***********************


// *** Mesh Module ***

typedef struct elmt_t
{
  unsigned int vertex[4];
} elmt;

// typedef struct msh_node
// {
//   xyz_point_t *node_list;
// } msh_node;

// typedef struct msh_elmt
// {
//   elmt_t *elmt_list;
// } msh_elmt;

typedef struct moab_mesh
{
  moab::Interface *intf = NULL;
  int num_nodes;
  int num_elmts;
} moab_mesh;


// int meshNodeInit(msh_node &, antn_cfg &);
// int meshNodeFinalize(msh_node &);

// int meshElmtInit(msh_elmt &, antn_cfg &);
// int meshElmtFinalize (msh_elmt &);

int meshInit(moab_mesh &, antn_cfg &);
int meshFinalize(moab_mesh &);

// *******************


// *** Velocity Model Module ***

typedef struct vm_datum
{
  double data[3];
} vm_datum;

typedef struct vmodel
{
  vm_datum *vm_list;
} vmodel;

int vmNodeInit(vmodel &, const moab_mesh &);
int vmNodeFinalize(vmodel &);
int vmElmtInit(vmodel &, const moab_mesh &);
int vmElmtFinalize(vmodel &);

int writeVMNodes(vmodel &, const antn_cfg &, const moab_mesh &);
int writeVMElmts(vmodel &, const antn_cfg &, const moab_mesh &);
int writeVMTags(vmodel &, const antn_cfg &, const moab_mesh &);


