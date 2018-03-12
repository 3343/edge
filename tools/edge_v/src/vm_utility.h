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
 * This is the header file for the utility routines of Edge-V.
 **/

#ifndef VM_UTILITY_H
#define VM_UTILITY_H

#include <iostream>
#include <iomanip>
#include <cassert>
#include <string>
#include <cmath>

extern "C" {
#include "ucvm.h"
}

#include "proj_api.h"
#include "moab/Core.hpp"
#include "vm_constants.h"


// *** Coordinate System ***

typedef struct xyz_point_t {
  real x;
  real y;
  real z;
} xyz_point_t;

typedef struct geo_point_t {
  real lon;
  real lat;
  real dep;
} geo_point_t;

// **************************


// *** Configuration File ***

typedef struct antn_cfg {
  std::string antn_cfg_fn;

  std::string ucvm_cfg_fn;
  std::string ucvm_model_list;
  std::vector< std::string > fault_input_fns;
  ucvm_ctype_t ucvm_cmode;
  int ucvm_type;

  real min_vp;
  real min_vs;
  real min_vs2;
  real max_vp_vs_ratio;

  std::string mesh_fn;
  geo_point_t hypoc;

  unsigned int elmt_type;

  std::string vm_node_fn;
  std::string vm_elmt_fn;
  std::string h5m_fn;

  unsigned int parallel_mode;
  unsigned int num_worker;
} antn_cfg;

int antnInit( antn_cfg &, const std::string & );

int ucvmInit( const antn_cfg & );

// // ***************************


// // *** Coordination System ***

// int xyz2geo( const xyz_point_t &, geo_point_t & );
// int geo2xyz( const geo_point_t &, xyz_point_t & );

// ***********************


// *** Parallel Module ***

typedef struct worker_reg {
  int_v worker_tid;

  projPJ pj_utm;
  projPJ pj_geo;

  int_v work_size;
  int_v num_prvt;
} worker_reg;

int workerInit( worker_reg &, int_v );

// int pjUtmInit(projPJ **);
// int pjUtmFinalize(projPJ *);

// int pjGeoInit(projPJ **);
// int pjGeoFinalize(projPJ *);

// ***********************


// *** Mesh Module ***

typedef struct elmt_t {
  unsigned int vertex[ELMTTYPE];
} elmt;

// typedef struct msh_node
// {
//   xyz_point_t *node_list;
// } msh_node;

// typedef struct msh_elmt
// {
//   elmt_t *elmt_list;
// } msh_elmt;

typedef struct moab_mesh {
  moab::Interface *intf = nullptr;
  int_v num_nodes;
  int_v num_elmts;
} moab_mesh;


// int meshNodeInit(msh_node &, antn_cfg &);
// int meshNodeFinalize(msh_node &);

// int meshElmtInit(msh_elmt &, antn_cfg &);
// int meshElmtFinalize (msh_elmt &);

int meshInit( moab_mesh &, antn_cfg & );
int meshFinalize( moab_mesh & );

// *******************


// *** Velocity Model Module ***

typedef struct vm_datum {
  real data[3];
} vm_datum;

typedef struct vmodel {
  vm_datum *vm_list;
} vmodel;


int vmNodeInit( vmodel &, const moab_mesh & );
int vmElmtInit( vmodel &, const moab_mesh & );

int vmNodeFinalize( vmodel & );
int vmElmtFinalize( vmodel & );

int writeVMNodes( vmodel &, const antn_cfg &, const moab_mesh & );
int writeVMElmts( vmodel &, const antn_cfg &, const moab_mesh & );
int writeVMTags(  vmodel &, const antn_cfg &, const moab_mesh & );


//*****************************************
// *** Fault Model ***

class Triangle{
public:
  xyz_point_t m_v1, m_v2, m_v3;

  Triangle() = default;
  Triangle( xyz_point_t i_v1, xyz_point_t i_v2, xyz_point_t i_v3 ) :
    m_v1( i_v1 ),
    m_v2( i_v2 ),
    m_v3( i_v3 )
  { }
  Triangle( double* );

  xyz_point_t centroid();
  xyz_point_t bary_to_phy( xyz_point_t );
  std::vector< Triangle > subdivide( unsigned int );
};

class FDatum{
public:
  // NOTE it's important that these are type double (not typedef-ed)
  // because they are used in conjuction with the moab MB_DOUBLE type
  std::vector< double > m_sFric;      // static coefficient of friction
  std::vector< double > m_nStress;    // initial normal stress
  std::vector< double > m_sStress;    // initial shear stress
};

class FModel{
public:
  // Number of cells in x,y directions
  int m_Nx, m_Ny;

  // Min and Max physical coordinates
  double m_xMin, m_xMax, m_yMin, m_yMax;

  // Reciprocal of the cell width in x and y directions
  double m_xScaleInv, m_yScaleInv;

  // Angle fault makes with pos x-axis, increasing toward pos z-axis
  double m_faultAngle;

  // File names to read from; Number of fused runs is inferred from here
  std::vector< std::string > m_filenames;

  // Main data array for the model, one FDatum per grid point,
  // datum at (nx, ny) corresponds to m_faultData[ nx + ny*(N_x+1) ]
  FDatum *m_faultData;

  FModel( std::vector< std::string > );
  ~FModel();

  // No copy or copy assignment constructors
  FModel( const FModel& ) = delete;
  FModel& operator=( const FModel& ) = delete;

  void populate();
  xyz_point_t proj_to_plane( xyz_point_t );
  int get_nearest_nx( xyz_point_t );
  int get_nearest_ny( xyz_point_t );
  FDatum get_datum_pt( xyz_point_t );
  FDatum get_datum_tri( Triangle );
};

moab::Range get_fault_faces( moab_mesh & );
int fault_antn( moab_mesh &, const antn_cfg & );



// *******************


#endif //! VM_UTILITY_H
