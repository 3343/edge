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
 * This is the header file for the velocity-related utility of EDGE-V.
 **/

#ifndef EDGE_VM_UTILITY_H
#define EDGE_VM_UTILITY_H

#include <iostream>
#include <iomanip>
#include <cassert>
#include <string>
#include <cmath>

extern "C" {
#include "ucvm.h"
}

#include "proj_api.h"
#include "moab/Interface.hpp"
#include "constants.h"
#include "Config.h"
#include "Geometry.h"

namespace edge_v {
  namespace vm {
    class Utility;
  }
}

class edge_v::vm::Utility {
  public:
    static int ucvmInit( const io::Config & );

    // *** Parallel Module ***
    typedef struct worker_reg {
      int_v   m_workerTid;

      projPJ  m_pjUtm;
      projPJ  m_pjGeo;

      int_v   m_workSize;
      int_v   m_numPrvt;
    } worker_reg;

    static int workerInit(       worker_reg &,
                    const int_v      & );

    // *** Mesh Module ***
    static int meshInit( moab_mesh &,
                  io::Config & );
    static int meshFinalize( moab_mesh & );


    // *** Velocity Model Module ***
    typedef struct vm_datum {
      real m_data[3];
    } vm_datum;

    typedef struct vmodel {
      vm_datum *m_vmList;
    } vmodel;

    static int vmNodeInit(       vmodel    &,
                    const moab_mesh & );
    static int vmElmtInit(       vmodel    &,
                    const moab_mesh & );

    static int vmNodeFinalize( vmodel & );
    static int vmElmtFinalize( vmodel & );

    static int writeVMNodes(       vmodel    &,
                      const io::Config  &,
                      const moab_mesh & );
    static int writeVMElmts(       vmodel    &,
                      const io::Config  &,
                      const moab_mesh & );
    static int writeVMTags(        vmodel    &,
                      const io::Config  &,
                      const moab_mesh & );

    // *** pos related ***
    typedef struct pos_datum {
      xyz_point_t m_xyzPts[ELMTTYPE];
      real        m_vs[ELMTTYPE];
    } pos_datum;

    typedef struct posModel {
      pos_datum *m_posList;
    } posModel;

    static int posInit(       posModel  &,
                       const moab_mesh & );
    static int posFinalize( posModel & );
    static int writePos( const posModel  &,
                         const io::Config  &,
                         const moab_mesh & );
};

#endif
