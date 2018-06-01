/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2016-2018, Regents of the University of California
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
 * Common functions for elastics.
 **/
#ifndef ELASTIC_COMMON_HPP
#define ELASTIC_COMMON_HPP

#include <cassert>
#include "constants.hpp"
#include "monitor/instrument.hpp"

namespace edge {
  namespace elastic {
    class common;
  }
}

class edge::elastic::common {
  public:
    /**
     * Gets the p-wave velocity for the given parameters.
     *
     * @param i_rho density rho.
     * @param i_lam lame parameter lambda
     * @param i_mu lame parameter mu
     * @return p-wave velocity
     **/
    static real_base getVelP( real_base i_rho,
                              real_base i_lam,
                              real_base i_mu ) {
      real_base l_cp = (i_lam + (real_base) 2.0 * i_mu)/i_rho;
      assert( l_cp > 0 );
      l_cp = std::sqrt( l_cp );

      return l_cp;
    }

    /**
     * Gets the s-wave velocity for the given parameters.
     *
     * @param i_rho density rho.
     * @param i_lam lame parameter lambda.
     * @param i_mu lame parameter mu.
     * @return p-wave velocity.
     **/
    static real_base getVelS( real_base i_rho,
                              real_base i_mu ) {
      real_base l_s = i_mu / i_rho;
      assert( l_s > 0 );
      l_s = std::sqrt( l_s );

      return l_s;
    }

    /**
     * Gets the time step according to the CFL-condition.
     *
     * @param i_rho density rho.
     * @param i_lam lame parameter lambda.
     * @param i_mu lame parameter mu.
     * @param i_dIns maximum insphere diameter.
     * @param i_scale scaled time step, e.g. influence of CFL number or order of convergence.
     **/
    static double getTimeStepCFL( real_base i_rho,
                                  real_base i_lam,
                                  real_base i_mu,
                                  real_base i_dIns,
                                  double    i_scale = 1.0 ) {
      double l_dt  = i_dIns;
             l_dt *= i_scale;
             l_dt /= getVelP( i_rho, i_lam, i_mu );
             l_dt /= 2.0*ORDER-1;

      assert( l_dt > 0 );

      return l_dt;
    }

    /**
     * Gets the time step statistics according to CFL.
     *
     * @param i_nElements number of elements.
     * @param i_elementChars element characteristics.
     * @param i_bgPars background parameters.
     * @param o_minDt minimum time step.
     * @param o_aveDt average time step.
     * @param o_maxDt maximum time step.
     **/
    static void getTimeStepStatsCFL(       int_el           i_nElements,
                                     const t_elementChars (*i_elementChars),
                                     const t_bgPars       (*i_bgPars)[1],
                                           double          &o_minDt,
                                           double          &o_aveDt,
                                           double          &o_maxDt ) {
      PP_INSTR_FUN("cfl_stats")

      // initialize time steps
      o_minDt = std::numeric_limits<double>::max();
      o_maxDt = 0;
      o_aveDt = 0;

      for( int_el l_el = 0; l_el < i_nElements; l_el++ ) {
        double l_dt = getTimeStepCFL( i_bgPars[l_el][0].rho,
                                      i_bgPars[l_el][0].lam,
                                      i_bgPars[l_el][0].mu,
                                      i_elementChars[l_el].inDia,
                                      SCALE_CFL );

        o_minDt  = std::min( o_minDt, l_dt );
        o_maxDt  = std::max( o_maxDt, l_dt );
        o_aveDt += l_dt / i_nElements;
      }
    }

    /**
     * Sets up the 3D transformation of quantities, from a local, face-aligned coordinate system to physical xyz-coords.
     * Remark: n is face's the normal, s and t the two tangents in xyz. All of which have length one.
     *
     * @param i_nx x-coordinate of the normal.
     * @param i_ny y-coordinate of the normal.
     * @param i_nz z-coordinate of the normal.
     * @param i_sx x-coordinate of the first tangent.
     * @param i_sy y-coordinate of the first tangent.
     * @param i_sz z-coordinate of the first tangent.
     * @param i_tx x-coordinate of the second tangent.
     * @param i_ty y-coordinate of the second tangent.
     * @param i_tz z-coordinate of the second tangent.
     * @param o_t will be set to the transformation matrix.
     **/
    template <typename TL_T_MESH, typename TL_T_TRAFO>
    static void setupTrafo3d( TL_T_MESH i_nx, TL_T_MESH i_ny, TL_T_MESH i_nz,
                              TL_T_MESH i_sx, TL_T_MESH i_sy, TL_T_MESH i_sz,
                              TL_T_MESH i_tx, TL_T_MESH i_ty, TL_T_MESH i_tz,
                              TL_T_TRAFO o_t[9][9] ) {
#include "impl/elastic/generated/Trafo3D.inc"
    }

    /**
     * Sets up the 3D transformation of quantities, from physical coordinates to a local, face-aligned coordinate system.
     * Remark: n is face's the normal, s and t the two tangents in xyz. All of which have length one.
     *
     * @param i_nx x-coordinate of the normal.
     * @param i_ny y-coordinate of the normal.
     * @param i_nz z-coordinate of the normal.
     * @param i_sx x-coordinate of the first tangent.
     * @param i_sy y-coordinate of the first tangent.
     * @param i_sz z-coordinate of the first tangent.
     * @param i_tx x-coordinate of the second tangent.
     * @param i_ty y-coordinate of the second tangent.
     * @param i_tz z-coordinate of the second tangent.
     * @param o_tm1 will be set to the transformation matrix.
     **/
    template <typename TL_T_MESH, typename TL_T_TRAFO>
    static void setupTrafoInv3d( TL_T_MESH i_nx, TL_T_MESH i_ny, TL_T_MESH i_nz,
                                 TL_T_MESH i_sx, TL_T_MESH i_sy, TL_T_MESH i_sz,
                                 TL_T_MESH i_tx, TL_T_MESH i_ty, TL_T_MESH i_tz,
                                 TL_T_TRAFO o_tm1[9][9] ) {
#include "impl/elastic/generated/TrafoInv3D.inc"
    }

    /**
     * Sets up the 2D transformation of quantities, from a local, face-aligned coordinate system to physical xy-coords.
     * Remark: n is face's the normal, which has length one.
     *
     * @param i_nx x-coordinate of the normal.
     * @param i_ny y-coordinate of the normal.
     * @param o_t will be set to the transformation matrix.
     **/
    template <typename TL_T_MESH, typename TL_T_TRAFO>
    static void setupTrafo2d( TL_T_MESH i_nx, TL_T_MESH i_ny,
                              TL_T_TRAFO o_t[5][5] ) {
#include "impl/elastic/generated/TrafoElastic2D.inc"
    }

    /**
     * Sets up the 3D transformation of quantities, from physical coordinates to a local, face-aligned coordinate system.
     * Remark: n is face's the normal, which has length one.
     *
     * @param i_nx x-coordinate of the normal.
     * @param i_ny y-coordinate of the normal.
     * @param o_tm1 will be set to the transformation matrix.
     **/
    template <typename TL_T_MESH, typename TL_T_TRAFO>
    static void setupTrafoInv2d( TL_T_MESH i_nx, TL_T_MESH i_ny,
                                 TL_T_TRAFO o_tm1[5][5] ) {
#include "impl/elastic/generated/TrafoInvElastic2D.inc"
    }
};

#endif
