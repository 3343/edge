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
 * This is the header file for the utility routines of Edge-V.
 **/

#ifndef VM_UTILITY_H
#define VM_UTILITY_H

#include <iostream>
#include <iomanip>
#include <cassert>
#include <string>
#include <cmath>
#include <map>

extern "C" {
#include "ucvm.h"
}

#include "proj_api.h"
#include "moab/Core.hpp"
#include "vm_constants.h"


// *** Coordinate System ***
typedef struct xyz_point_t {
  real m_x;
  real m_y;
  real m_z;
} xyz_point_t;

typedef struct geo_point_t {
  real m_lon;
  real m_lat;
  real m_dep;
} geo_point_t;
// *************************


// *** Configuration File ***
typedef struct antn_cfg {
  std::string   m_antnCfgFn;

  std::string   m_ucvmCfgFn;
  std::string   m_ucvmModelList;
  ucvm_ctype_t  m_ucvmCmode;
  int           m_ucvmType;

  real          m_minVp;
  real          m_minVs;
  real          m_minVs2;
  real          m_maxVpVsRatio;
  real          m_elmtsPerWave;

  std::string   m_meshFn;
  geo_point_t   m_hypoc;

  int           m_tetRefinement;
  std::vector< std::string > m_faultInputFns;
  std::string   m_vmNodeFn;
  std::string   m_vmElmtFn;
  std::string   m_h5mFn;
  std::string   m_posFn;

  unsigned int  m_parallelMode;
  unsigned int  m_numWorker;
} antn_cfg;

int antnInit(       antn_cfg    &,
              const std::string & );

int ucvmInit( const antn_cfg & );
// ***********************


// *** Parallel Module ***
typedef struct worker_reg {
  int_v   m_workerTid;

  projPJ  m_pjUtm;
  projPJ  m_pjGeo;

  int_v   m_workSize;
  int_v   m_numPrvt;
} worker_reg;

int workerInit(       worker_reg &,
                const int_v      & );
// ***********************


// *** Mesh Module ***
typedef struct elmt_t {
  unsigned int m_vertex[ELMTTYPE];
} elmt;

typedef struct moab_mesh {
  moab::Interface *m_intf = nullptr;
  int_v            m_numNodes;
  int_v            m_numElmts;
} moab_mesh;

int meshInit( moab_mesh &,
              antn_cfg  & );
int meshFinalize( moab_mesh & );
// *******************


// *** Velocity Model Module ***
typedef struct vm_datum {
  real m_data[3];
} vm_datum;

typedef struct vmodel {
  vm_datum *m_vmList;
} vmodel;

int vmNodeInit(       vmodel    &,
                const moab_mesh & );
int vmElmtInit(       vmodel    &,
                const moab_mesh & );

int vmNodeFinalize( vmodel & );
int vmElmtFinalize( vmodel & );

int writeVMNodes(       vmodel    &,
                  const antn_cfg  &,
                  const moab_mesh & );
int writeVMElmts(       vmodel    &,
                  const antn_cfg  &,
                  const moab_mesh & );
int writeVMTags(        vmodel    &,
                  const antn_cfg  &,
                  const moab_mesh & );
// *******************


// *** pos related ***
typedef struct pos_datum {
  xyz_point_t m_xyzPts[ELMTTYPE];
  real        m_vs[ELMTTYPE];
} pos_datum;

typedef struct posModel {
  pos_datum *m_posList;
} posModel;

int posInit(       posModel  &,
             const moab_mesh & );
int posFinalize( posModel & );
int writePos( const posModel  &,
              const antn_cfg  &,
              const moab_mesh & );
// *******************


// *** Fault Model ***
class Triangle {
public:
  xyz_point_t m_v1, m_v2, m_v3;

  Triangle() = default;
  Triangle( const xyz_point_t &i_v1,
            const xyz_point_t &i_v2,
            const xyz_point_t &i_v3 ) : m_v1( i_v1 ), m_v2( i_v2 ), m_v3( i_v3 ) {}
  Triangle( const real * );

  xyz_point_t centroid();
  xyz_point_t baryToPhy( const xyz_point_t & );
  std::vector< Triangle > subdivide( const int & );
};

typedef std::map< std::string, std::vector< real > > FDatum;

class FModel {
public:
  //! Number of cells in x,y directions
  int m_Nx, m_Ny;

  //! Min and Max physical coordinates
  real m_xMin, m_xMax, m_yMin, m_yMax;

  //! Reciprocal of the cell width in x and y directions
  real m_xScaleInv, m_yScaleInv;

  //! Angle fault makes with pos x-axis, increasing toward pos z-axis
  real m_faultAngle;

  // List of quantities being described by the input file
  std::vector< std::string > m_qtyNames;

  //! File names to read from; Number of fused runs is inferred from here
  std::vector< std::string > m_filenames;

  //! Main data array for the model, one FDatum per grid point,
  //! datum at (nx, ny) corresponds to m_faultData[ nx + ny*(N_x+1) ]
  FDatum *m_faultData;

  FModel( const std::vector< std::string > & );
  ~FModel();

  //! No copy constructor or copy assignment
  FModel( const FModel & ) = delete;
  FModel & operator=( const FModel & ) = delete;

  void populate();
  xyz_point_t projToPlane( const xyz_point_t & );
  int getNearestNx( const xyz_point_t & );
  int getNearestNy( const xyz_point_t & );
  FDatum getDatumPt( const xyz_point_t & );
  FDatum getDatumTri( Triangle & );
};

moab::Range getFaultFaces( moab_mesh & );
int faultAntn( const  antn_cfg &,
                      moab_mesh & );
// *******************

#endif //! VM_UTILITY_H
