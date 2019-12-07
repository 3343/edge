/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2019, Alexander Breuer
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
 * Common functionality of the solvers for the advection equations.
 **/

#ifndef EDGE_ADVECTION_SOLVERS_COMMON_HPP
#define EDGE_ADVECTION_SOLVERS_COMMON_HPP

#include "io/logging.h"

namespace edge {
  namespace advection {
    namespace solvers {
      class common;
    }
  }
}

class edge::advection::solvers::common {
  //private:
    /**
     * Computes the time step in a given element.
     *
     * @param i_inDia insphere diameter of the element.
     * @param i_bgPars background parameters of an element.
     * @param i_cfl cfl number.
     **/
    static double computeTimeStep( double           i_inDia,
                                   t_bgPars const & i_bgPars,
                                   double           i_cfl  = 0.25 ) {
       // get abs advection speed
#if PP_N_DIM == 1
       double l_lambda = std::abs( i_bgPars.a );
#endif
#if PP_N_DIM == 2
       double l_lambda  = i_bgPars.a * i_bgPars.a;
              l_lambda += i_bgPars.b * i_bgPars.b;
              l_lambda  = std::sqrt(l_lambda);
#endif
#if PP_N_DIM == 3
       double l_lambda  = i_bgPars.a * i_bgPars.a;
              l_lambda += i_bgPars.b * i_bgPars.b;
              l_lambda += i_bgPars.c * i_bgPars.c;
              l_lambda  = std::sqrt(l_lambda);
#endif

       // compute time step
       double l_dt  = i_inDia / l_lambda;
              l_dt *= i_cfl;

       // scale with order
       l_dt /= 2.0*ORDER-1;

       return l_dt;
    }

  public:
    /**
     * Gets the time step statitics over the elements.
     *
     * @param i_nEls number of elements.
     * @param i_elChars element characteristics.
     * @param i_bgPars background parameters of an element.
     * @param o_minDt set to minimum occurring time step.
     * @param o_aveDt set to average occuring time step.
     * @param o_maxDt set to maximum occuring time step.
     **/
    static void getTimeStepStatistics( std::size_t             i_nEls,
                                       t_elementChars const  * i_elChars,
                                       t_bgPars       const (* i_bgPars)[1],
                                       double                & o_minDt,
                                       double                & o_aveDt,
                                       double                & o_maxDt ) {
      // intialize statistics
      o_minDt = std::numeric_limits< double >::max();
      o_aveDt = 0;
      o_maxDt = 0;

      for( std::size_t l_el = 0; l_el < i_nEls; l_el++ ) {
        // compute time step
        double l_dT = computeTimeStep( i_elChars[l_el].inDia,
                                       i_bgPars[l_el][0],
                                       SCALE_CFL );

        // add element to stats
        o_minDt  = std::min( o_minDt, l_dT );
        o_aveDt += l_dT;
        o_maxDt  = std::max( o_maxDt, l_dT );
      }

      // average
      o_aveDt /= i_nEls;
    }

    /**
     * Sets up the solvers.
     * We store one solver for the elements own contribution in pos [0, N_FACES] and the neigh contribution in [N_FACES+1, 2*N_FACES].
     *
     * @param i_nElsIn number of inner elements.
     * @param i_nElsSe number of send elements.
     * @param i_nFas number of faces in the mesh.
     * @param i_nCommElFa nuber of communicating faces.
     * @param i_recvFa local ids of the faces w.r.t. their position in the elements.
     * @param i_recvEl ids of the respective element for the communicating faces.
     * @param i_elVe vertices adjacent to the elements.
     * @param i_faEl elements adjacent to the faces.
     * @param i_elFa faces adjacent to the elements.
     * @param i_elFaEl adjacent elements (faces as bridge).
     * @param i_veChars vertex characteristics.
     * @param i_faChars mesh characteristics of the faces.
     * @param i_bgPars backgrorund parameters (velocities).
     * @param i_bgParsRe background parameter of the receive element, matching the el-fa comm strucuture.
     * @param o_fs will be set to the flux solvers.
     *
     * @paramt TL_T_REAL real type.
     **/
    template< typename TL_T_REAL >
    static void setupSolvers( std::size_t             i_nElsIn,
                              std::size_t             i_nElsSe,
                              std::size_t             i_nFas,
                              std::size_t             i_nCommElFa,
                              unsigned short const  * i_recvFa,
                              std::size_t    const  * i_recvEl,
                              std::size_t    const (* i_elVe)[ C_ENT[T_SDISC.ELEMENT].N_VERTICES ],
                              std::size_t    const (* i_faEl)[2],
                              std::size_t    const (* i_elFa)[C_ENT[T_SDISC.ELEMENT].N_FACES],
                              std::size_t    const (* i_elFaEl)[C_ENT[T_SDISC.ELEMENT].N_FACES],
                              t_vertexChars  const  * i_veChars,
                              t_faceChars    const  * i_faChars,
                              t_bgPars       const (* i_bgPars)[1],
                              t_bgPars       const  * i_bgParsRe,
                              TL_T_REAL            (* o_fs)[C_ENT[T_SDISC.ELEMENT].N_FACES*2] ) {
      std::size_t l_nEls = i_nElsIn + i_nElsSe;

      // iterate over elements and reset flux solvers
      for( std::size_t l_el = 0; l_el < l_nEls; l_el++ ) {
        for( int_md l_fa = 0; l_fa < C_ENT[T_SDISC.ELEMENT].N_FACES*2; l_fa++ ) {
          o_fs[l_el][l_fa] = 0;
        }
      }

      // iterate over faces
      for( std::size_t l_fa = 0; l_fa < i_nFas; l_fa++ ) {
        // get left and right element
        std::size_t l_elL = i_faEl[l_fa][0];
        std::size_t l_elR = i_faEl[l_fa][1];

        // determine existance of elements (boundary conditions)
        bool l_exL = (l_elL < l_nEls);
        bool l_exR = (l_elR < l_nEls);

        // find ids of local contribution flux solvers for the elements
        unsigned short l_fLocIdL = std::numeric_limits<unsigned short>::max();
        unsigned short l_fLocIdR = std::numeric_limits<unsigned short>::max();

        for( unsigned short l_fi = 0; l_fi < C_ENT[T_SDISC.ELEMENT].N_FACES; l_fi++ ) {
          if( l_exL && i_elFa[l_elL][l_fi] == l_fa ) {
            // face should only by found once
            EDGE_CHECK_EQ( l_fLocIdL, std::numeric_limits< unsigned short >::max() );
            l_fLocIdL = l_fi;
          }
          if( l_exR && i_elFa[l_elR][l_fi] == l_fa ) {
            // face should only by found once
            EDGE_CHECK_EQ( l_fLocIdR, std::numeric_limits <unsigned short >::max() );
            l_fLocIdR = l_fi;
          }
        }

        // derive advection speeds in normal direction
        TL_T_REAL l_advNormalL = std::numeric_limits<TL_T_REAL>::max();
        TL_T_REAL l_advNormalR = std::numeric_limits<TL_T_REAL>::max();

        // project left speed on normal
        if( l_exL ) {
          l_advNormalL  = i_bgPars[l_elL][0].a * i_faChars[l_fa].outNormal[0];
#if PP_N_DIM > 1
          l_advNormalL += i_bgPars[l_elL][0].b * i_faChars[l_fa].outNormal[1];
#endif
#if PP_N_DIM > 2
          l_advNormalL += i_bgPars[l_elL][0].c * i_faChars[l_fa].outNormal[2];
#endif
        }

        // check for receive element if "right" doesn't exist
        if( !l_exR && l_elL >= i_nElsIn ) {
          for( std::size_t l_co = 0; l_co < i_nCommElFa; l_co++ ) {
            if( i_recvEl[l_co] == l_elL && i_recvFa[l_co] == l_fLocIdL ) {
              l_advNormalR  = i_bgParsRe[l_co].a * i_faChars[l_fa].outNormal[0];
#if PP_N_DIM > 1
              l_advNormalR += i_bgParsRe[l_co].b * i_faChars[l_fa].outNormal[1];
#endif
#if PP_N_DIM > 2
              l_advNormalR += i_bgParsRe[l_co].c * i_faChars[l_fa].outNormal[2];
#endif
              break;
            }
          }
        }
        else if( l_exR ) {
          l_advNormalR  = i_bgPars[l_elR][0].a * i_faChars[l_fa].outNormal[0];
#if PP_N_DIM > 1
          l_advNormalR += i_bgPars[l_elR][0].b * i_faChars[l_fa].outNormal[1];
#endif
#if PP_N_DIM > 2
          l_advNormalR += i_bgPars[l_elR][0].c * i_faChars[l_fa].outNormal[2];
#endif
        }

        bool l_periodic = (i_faChars[l_fa].spType & PERIODIC) == PERIODIC;
        EDGE_CHECK( !l_exL || l_fLocIdL != std::numeric_limits< unsigned short >::max() || l_periodic ) << l_fa;
        EDGE_CHECK( !l_exR || l_fLocIdR != std::numeric_limits< unsigned short >::max() || l_periodic ) << l_fa;

        // disable elements without face ids (allowed for periodic faces, where elFa only points to one of them)
        l_exL = l_exL && l_fLocIdL != std::numeric_limits< unsigned short >::max();
        l_exR = l_exL && l_fLocIdR != std::numeric_limits< unsigned short >::max();

        // ids of neigh contribution flux solvers for the elements
        unsigned short l_fNegIdL = l_exL ? l_fLocIdL + C_ENT[T_SDISC.ELEMENT].N_FACES :
                                   std::numeric_limits< unsigned short >::max();
        unsigned short l_fNegIdR = l_exR ? l_fLocIdR + C_ENT[T_SDISC.ELEMENT].N_FACES :
                                   std::numeric_limits< unsigned short >::max();

        /*
         * setup solvers
         *
         * n is the normal -> positive advection speeds go L->R
         *                 -> negative advection speeds go R->L
         */
         // check for normal direction of left element's wave speeds
         if( l_advNormalL > 0 ) {
           if( l_exL ) o_fs[l_elL][l_fLocIdL] = -l_advNormalL;
           if( l_exR ) o_fs[l_elR][l_fNegIdR] = l_advNormalL;
         }

         // check for normal direction of right element's wave speeds
         if( l_advNormalR < 0 ) {
           if( l_exL ) o_fs[l_elL][l_fNegIdL] = -l_advNormalR;
           if( l_exR ) o_fs[l_elR][l_fLocIdR] = l_advNormalR;
         }

        // derive vertex coords
        double l_veCoordsL[N_DIM][C_ENT[T_SDISC.ELEMENT].N_VERTICES];
        double l_veCoordsR[N_DIM][C_ENT[T_SDISC.ELEMENT].N_VERTICES];

        for( unsigned int l_ve = 0; l_ve < C_ENT[T_SDISC.ELEMENT].N_VERTICES; l_ve++ ) {
          std::size_t l_veIdL = std::numeric_limits<std::size_t>::max();
          std::size_t l_veIdR = std::numeric_limits<std::size_t>::max();
          if( l_exL == true ) l_veIdL = i_elVe[l_elL][l_ve];
          if( l_exR == true ) l_veIdR = i_elVe[l_elR][l_ve];

          for( unsigned int l_dim = 0; l_dim < N_DIM; l_dim++ ) {
            if( l_exL == true ) l_veCoordsL[l_dim][l_ve] = i_veChars[l_veIdL].coords[l_dim];
            if( l_exR == true ) l_veCoordsR[l_dim][l_ve] = i_veChars[l_veIdR].coords[l_dim];
          }
        }

        // compute determinant of the mapping's jacobian
        double l_jacL[N_DIM][N_DIM], l_jacR[N_DIM][N_DIM];
        for( unsigned short l_dm1 = 0; l_dm1 < N_DIM; l_dm1++ )
          for( unsigned short l_dm2 = 0; l_dm2 < N_DIM; l_dm2++ )
            l_jacL[l_dm1][l_dm2] = l_jacR[l_dm1][l_dm2] = std::numeric_limits<double>::max();
        double l_jL = std::numeric_limits<double>::max();
        double l_jR = std::numeric_limits<double>::max();
        if( l_exL == true ) linalg::Mappings::evalJac( T_SDISC.ELEMENT, l_veCoordsL[0], l_jacL[0] );
        if( l_exR == true ) linalg::Mappings::evalJac( T_SDISC.ELEMENT, l_veCoordsR[0], l_jacR[0] );
#if PP_N_DIM == 1
        if( l_exL == true ) l_jL = l_jacL[0][0];
        if( l_exR == true ) l_jR = l_jacR[0][0];
#elif PP_N_DIM == 2
        if( l_exL == true ) l_jL = linalg::Matrix::det( l_jacL );
        if( l_exR == true ) l_jR = linalg::Matrix::det( l_jacR );
#elif PP_N_DIM == 3
        if( l_exL == true ) l_jL = linalg::Matrix::det( l_jacL );
        if( l_exR == true ) l_jR = linalg::Matrix::det( l_jacR );
#else
#error number of dimensions not supported.
#endif

        if( T_SDISC.ELEMENT == TET4 && ORDER == 1 ) {
          // Remark: For tets the determinant is not equal to the volume of the tets.
          //         We apply the effect of the inverse mass matrix right here.
          //         Effectively this just scales the flux solver by the volume of the tets.
          if( l_exL == true ) l_jL /= (double) 6;
          if( l_exR == true ) l_jR /= (double) 6;
        }

        // ensure positive determinants
        EDGE_CHECK( l_exL == false || l_jL > 0 );  EDGE_CHECK( l_exR == false || l_jR > 0 );
        EDGE_CHECK_GT( i_faChars[l_fa].area, 0 );

        // scale the flux solvers by dets of volume integration
        if( l_exL == true ) o_fs[l_elL][l_fLocIdL] /= l_jL;
        if( l_exR == true ) o_fs[l_elR][l_fLocIdR] /= l_jR;
        if( l_exL == true ) o_fs[l_elL][l_fNegIdL] /= l_jL;
        if( l_exR == true ) o_fs[l_elR][l_fNegIdR] /= l_jR;

        // scale by the face surfaces
        if( l_exL == true ) o_fs[l_elL][l_fLocIdL] *= i_faChars[l_fa].area;
        if( l_exR == true ) o_fs[l_elR][l_fLocIdR] *= i_faChars[l_fa].area;
        if( l_exL == true ) o_fs[l_elL][l_fNegIdL] *= i_faChars[l_fa].area;
        if( l_exR == true ) o_fs[l_elR][l_fNegIdR] *= i_faChars[l_fa].area;

#if defined PP_T_ELEMENTS_TET4
        // the det of the surface-jac is twice the area of the triangle
        // TODO: Replace this formulation with Jacobi-determinants.
        if( l_exL == true ) o_fs[l_elL][l_fLocIdL] *= 2;
        if( l_exR == true ) o_fs[l_elR][l_fLocIdR] *= 2;
        if( l_exL == true ) o_fs[l_elL][l_fNegIdL] *= 2;
        if( l_exR == true ) o_fs[l_elR][l_fNegIdR] *= 2;
#endif
      }
    }
};

#endif