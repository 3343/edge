/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2020, Alexander Breuer
 * Copyright (c) 2018, Regents of the University of California
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
 * Interface to UCVM.
 **/
#ifndef EDGE_V_IO_UCVM_H
#define EDGE_V_IO_UCVM_H

#include <string>

namespace edge_v {
  namespace io {
    class Ucvm;
  }
}

/**
 * Interface to UCVM.
 */
class edge_v::io::Ucvm {
  public:
    /**
     * Constructor which initializes the UCVM interface.
     * 
     * @param i_config path to the UCVM config.
     * @param i_models UCVM models.
     * @param i_crdMode UCVM coordinate mode.
     **/
    Ucvm( std::string const & i_config,
          std::string const & i_models,
          std::string const & i_crdMode );

    /**
     * Get the velocities from UCVM.
     * Given a point p, the following is done:
     *   1) All coordinates (x, y, z) of p are used to obtain p* = i_trafoSrc(p).
     *   2) x and y of p* are transformed by proj according to the given projections (ignoring z).
     *   3) The obtained x and y coordinates plus the original z of p* are used to query UCVM.
     * 
     * @param i_nPts number of points.
     * @param i_trafoSrc trafo, applied to the source coordinates before the proj.4 projection.
     * @param i_projSrc proj.4 string of the source projection.
     * @param i_projDes proh.4 string of the destination projection (used to query UCVM).
     * @param i_ucvmType ucvm type: gtl, crust or cmb.
     * @param i_pts points for which the vp, vs and rho are queried.
     * @param o_vps will be set to p-wave velocities of the points.
     * @param o_vss will be set to s-wave velocities of the points.
     * @param o_rhos will be set to densities of the points.
     **/
    void getVels( std::size_t          i_nPts,
                  double      const    i_trafoSrc[3][3],
                  std::string const  & i_projSrc,
                  std::string const  & i_projDes,
                  std::string const  & i_ucvmType,
                  double            (* i_pts)[3],
                  float              * o_vps,
                  float              * o_vss,
                  float              * o_rhos );
};

#endif