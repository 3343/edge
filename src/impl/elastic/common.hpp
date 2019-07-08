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
 * Common functions for seismic simulations.
 **/
#ifndef EDGE_SEISMIC_COMMON_HPP
#define EDGE_SEISMIC_COMMON_HPP

#include <cassert>
#include <cmath>
#include "constants.hpp"
#include "monitor/instrument.hpp"
#include "io/logging.h"

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
     * Gets the Jacobians of the elastic wave equations in 2 dimensions.
     *
     * @param i_rho density rho.
     * @param i_lam Lame parameter lambda.
     * @param i_mu Lame parameter mu.
     * @param o_A will be set to Jacobians.
     **/
    template <typename T>
    static void getJac2D( T i_rho, T i_lam, T i_mu,
                          T o_A[2][5][5] ) {
      // reset to zero
      for( unsigned short l_di = 0; l_di < 2; l_di++ ) {
        for( int_qt l_q1 = 0; l_q1 < 5; l_q1++ ) {
          for( int_qt l_q2 = 0; l_q2 < 5; l_q2++ ) {
            o_A[l_di][l_q1][l_q2] = 0;
          }
        }
      }

     /*
      * Jacobians in q_t = A(x,y) * q_x + B(x,y) * q_y
      *
      * A:
      *    _____0__1_______2_______________3____4
      *  0|     0, 0,      0, -lambda - 2*mu,   0|0
      *  1|     0, 0,      0,        -lambda,   0|1
      *  2|     0, 0,      0,              0, -mu|2
      *  3|-1/rho, 0,      0,              0,   0|3
      *  4|     0, 0, -1/rho,              0,   0|4
      *    -----0--1-------2---------------3----4
      *
      * B:
      *    0_______1_______2____3_______________4
      *  0|0,      0,      0,   0,        -lambda|0
      *  1|0,      0,      0,   0, -lambda - 2*mu|1
      *  2|0,      0,      0, -mu,              0|2
      *  3|0,      0, -1/rho,   0,              0|3
      *  4|0, -1/rho,      0,   0,              0|4
      *    0-------1-------2----3---------------4
      */
      o_A[0][0][3] = -i_lam - T(2) * i_mu;
      o_A[0][1][3] = -i_lam;
      o_A[0][2][4] = -i_mu;
      o_A[0][3][0] = -T(1) / i_rho;
      o_A[0][4][2] = -T(1) / i_rho;

      o_A[1][0][4] = -i_lam;
      o_A[1][1][4] = -i_lam - T(2) * i_mu;
      o_A[1][2][3] = -i_mu;
      o_A[1][3][2] = -T(1) / i_rho;
      o_A[1][4][1] = -T(1) / i_rho;
    }

    /**
     * Gets the Jacobians of the elastic wave equations in 3 dimensions.
     *
     * @param i_rho density rho.
     * @param i_lam Lame parameter lambda.
     * @param i_mu Lame parameter mu.
     * @param o_A will be set to Jacobians.
     **/
    template <typename T>
    static void getJac3D( T i_rho, T i_lam, T i_mu,
                          T o_A[3][9][9] ) {
      // reset to zero
      for( unsigned short l_di = 0; l_di < 3; l_di++ ) {
        for( int_qt l_q1 = 0; l_q1 < 9; l_q1++ ) {
          for( int_qt l_q2 = 0; l_q2 < 9; l_q2++ ) {
            o_A[l_di][l_q1][l_q2] = 0;
          }
        }
      }

      /*
       * Jacobians in q_t = A(x,y,z) * q_x + B(x,y,z) * q_y + C(x,y,z) * q_z
       *
       * A:
       *   _____0__1__2_______3__4_______5_______________6____7____8
       * 0|     0, 0, 0,      0  0,      0, -lambda - 2*mu,   0,   0|0
       * 1|     0, 0, 0,      0, 0,      0,        -lambda,   0,   0|1
       * 2|     0, 0, 0,      0, 0,      0,        -lambda,   0,   0|2
       * 3|     0, 0, 0,      0, 0,      0,              0, -mu,   0|3
       * 4|     0, 0, 0,      0, 0,      0,              0,   0,   0|4
       * 5|     0, 0, 0,      0, 0,      0,              0,   0, -mu|5
       * 6|-1/rho, 0, 0,      0, 0,      0,              0,   0,   0|6
       * 7|     0, 0, 0, -1/rho, 0,      0,              0,   0,   0|7
       * 8|     0, 0, 0,      0, 0, -1/rho,              0,   0,   0|8
       *   -----0--1--2-------3--4-------5---------------6----7----8
       *
       * B:
       *   0_______1__2_______3_______4__5____6_________________7______8
       * 0|0,      0, 0,      0       0, 0,   0,          -lambda,     0|0
       * 1|0,      0, 0,      0,      0, 0,   0,   -lambda - 2*mu,     0|1
       * 2|0,      0, 0,      0,      0, 0,   0,          -lambda,     0|2
       * 3|0,      0, 0,      0,      0, 0, -mu,                0,     0|3
       * 4|0,      0, 0,      0,      0, 0,   0,                0,   -mu|4
       * 5|0,      0, 0,      0,      0, 0,   0,                0,     0|5
       * 6|0,      0, 0, -1/rho,      0, 0,   0,                0,     0|6
       * 7|0, -1/rho, 0,      0,      0, 0,   0,                0,     0|7
       * 8|0,      0, 0,      0, -1/rho, 0,   0,                0,     0|8
       *   0-------1--2-------3-------4--5---------------6------7------8
       *
       * C:
       *   0__1_______2__3_______4_______5____6____7_______________8
       * 0|0, 0,      0, 0       0,      0,   0,   0,        -lambda|0
       * 1|0, 0,      0, 0,      0,      0,   0,   0,        -lambda|1
       * 2|0, 0,      0, 0,      0,      0,   0,   0, -lambda - 2*mu|2
       * 3|0, 0,      0, 0,      0,      0,   0,   0,              0|3
       * 4|0, 0,      0, 0,      0,      0,   0, -mu,              0|4
       * 5|0, 0,      0, 0,      0,      0, -mu,   0,              0|5
       * 6|0, 0,      0, 0,      0, -1/rho,   0,   0,              0|6
       * 7|0, 0,      0, 0, -1/rho,      0,   0,   0,              0|7
       * 8|0, 0, -1/rho, 0,      0,      0,   0,   0,              0|8
       *   0--1-------2--3-------4-------5----6----7---------------8
       */
      o_A[0][0][6] = -i_lam - T(2) * i_mu;
      o_A[0][1][6] = -i_lam;
      o_A[0][2][6] = -i_lam;
      o_A[0][3][7] = -i_mu;
      o_A[0][5][8] = -i_mu;
      o_A[0][6][0] = -T(1) / i_rho;
      o_A[0][7][3] = -T(1) / i_rho;
      o_A[0][8][5] = -T(1) / i_rho;

      o_A[1][0][7] = -i_lam;
      o_A[1][1][7] = -i_lam - T(2) * i_mu;
      o_A[1][2][7] = -i_lam;
      o_A[1][3][6] = -i_mu;
      o_A[1][4][8] = -i_mu;
      o_A[1][6][3] = -T(1) / i_rho;
      o_A[1][7][1] = -T(1) / i_rho;
      o_A[1][8][4] = -T(1) / i_rho;

      o_A[2][0][8] = -i_lam;
      o_A[2][1][8] = -i_lam;
      o_A[2][2][8] = -i_lam - T(2) * i_mu;
      o_A[2][4][7] = -i_mu;
      o_A[2][5][6] = -i_mu;
      o_A[2][6][5] = -T(1) / i_rho;
      o_A[2][7][4] = -T(1) / i_rho;
      o_A[2][8][2] = -T(1) / i_rho;
    }

    /**
     * Gets the Jacobians of the elastic wave equations.
     *
     * @param i_rho density rho.
     * @param i_lam Lame parameter lambda.
     * @param i_mu Lame parameter mu.
     * @param o_A will be set to Jacobians.
     * @param i_nDim number of dimensions (2 or 3).
     **/
    template< typename TL_T_REAL >
    static void getJac( TL_T_REAL       i_rho,
                        TL_T_REAL       i_lam,
                        TL_T_REAL       i_mu,
                        TL_T_REAL     * o_A,
                        unsigned short  i_nDim ) {
      if( i_nDim == 2 )      getJac2D( i_rho, i_lam, i_mu, (TL_T_REAL (*)[5][5]) o_A );
      else if( i_nDim == 3 ) getJac3D( i_rho, i_lam, i_mu, (TL_T_REAL (*)[9][9]) o_A );
      else                   EDGE_LOG_FATAL << "dimensions not supported: " << i_nDim;
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
