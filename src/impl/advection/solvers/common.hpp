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
 * Common functionality of the solvers for the advection equations.
 **/

#ifndef EDGE_ADVECTION_SOLVERS_COMMON_HPP
#define EDGE_ADVECTION_SOLVERS_COMMON_HPP

#include <cassert>

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
     * @param i_inSphereDiameter insphere diameter of the element.
     * @param i_bgPars background parameters of an element.
     * @param i_cfl cfl number.
     **/
    static double computeTimeStep(       double    i_inSphereDiameter,
                                   const t_bgPars &i_bgPars,
                                         double    i_cfl  = 0.25 ) {
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
       double l_dT  = i_inSphereDiameter / l_lambda;
              l_dT *= i_cfl;

       // scale with order
       l_dT /= 2.0*ORDER-1;

       return l_dT;
    }

  public:
    /**
     * Gets the time step statitics over the elements.
     *
     * @param i_nElements number of elements.
     * @param i_elementChars element characteristics.
     * @param i_bgPars background parameters of an element.
     * @param o_minDt set to minimum occurring time step.
     * @param o_aveDt set to average occuring time step.
     * @param o_maxDt set to maximum occuring time step.
     **/
    static void getTimeStepStatistics(       int_el           i_nElements,
                                       const t_elementChars  *i_elementChars,
                                       const t_bgPars       (*i_bgPars)[1],
                                             double          &o_minDt,
                                             double          &o_aveDt,
                                             double          &o_maxDt ) {
      // intialize statistics
      o_minDt = std::numeric_limits< double >::max();
      o_aveDt = 0;
      o_maxDt = 0;

      for( int_el l_element = 0; l_element < i_nElements; l_element++ ) {
        // compute time step
        double l_dT = computeTimeStep( i_elementChars[l_element].inDia,
                                       i_bgPars[l_element][0],
                                       SCALE_CFL );

        // add element to stats
        o_minDt  = std::min( o_minDt, l_dT );
        o_aveDt += l_dT;
        o_maxDt  = std::max( o_maxDt, l_dT );
      }

      // average
      o_aveDt /= i_nElements;
    }

    /**
     * Sets up the solvers.
     * We store one solver for the elements own contribution in pos [0, N_FACES] and the neigh contribution in [N_FACES+1, 2*N_FACES].
     *
     * @param i_nElements number of elements in the mesh.
     * @param i_nFaces number of faces in the mesh.
     * @param i_faMeDa mapping of face ids: mesh-to-data.
     * @param i_faDaMe mapping of face ids: data-to-mesh.
     * @param i_elMeDa mapping of element ids: mesh-to-data.
     * @param i_elDaMe mapping of element ids: data-to-mesh.
     * @param i_elVe vertices adjacent to the elements.
     * @param i_faceAdjacentElements elements adjacent to the faces.
     * @param i_elementAdjacentFaces faces adjacent to the elements.
     * @param i_vertexChars vertex characteristics.
     * @param i_faceChars mesh characteristics of the faces.
     * @param i_bgPars backgrorund parameters (velocities).
     * @param o_fluxSolvers will be set to the flux solvers.
     **/
    static void setupSolvers(       int_el                 i_nElements,
                                    int_el                 i_nFaces,
                              const std::vector< int_el > &i_faMeDa,
                              const std::vector< int_el > &i_faDaMe,
                              const std::vector< int_el > &i_elMeDa,
                              const std::vector< int_el > &i_elDaMe,
                              const int_el               (*i_elVe)[ C_ENT[T_SDISC.ELEMENT].N_VERTICES ],
                              const int_el               (*i_faceAdjacentElements)[2],
                              const int_el               (*i_elementAdjacentFaces)[C_ENT[T_SDISC.ELEMENT].N_FACES],
                              const t_vertexChars         *i_vertexChars,
                              const t_faceChars           *i_faceChars,
                              const t_bgPars             (*i_bgPars)[1],
                                    real_base            (*o_fluxSolvers)[C_ENT[T_SDISC.ELEMENT].N_FACES*2] ) {
      // iterate over elements and reset flux solvers
      for( int_el l_el = 0; l_el < i_nElements; l_el++ ) {
        for( int_md l_fa = 0; l_fa < C_ENT[T_SDISC.ELEMENT].N_FACES*2; l_fa++ ) {
          o_fluxSolvers[l_el][l_fa] = 0;
        }
      }

      // iterate over faces
      for( int_el l_fa = 0; l_fa < i_nFaces; l_fa++ ) {
        // check that the face are unique
        EDGE_CHECK( l_fa == i_faMeDa[ i_faDaMe[l_fa] ] );

        // get left and right element
        int_el l_elL = i_faceAdjacentElements[l_fa][0];
        int_el l_elR = i_faceAdjacentElements[l_fa][1];

        // determine existance of elements (boundary conditions)
        bool l_exL = (l_elL < i_nElements);
        bool l_exR = (l_elR < i_nElements);

        // compute advection speeds
        real_base l_advSpeedL  = std::numeric_limits<real_base>::max();
        real_base l_advSpeedR  = std::numeric_limits<real_base>::max();
        real_base l_advNormalL = std::numeric_limits<real_base>::max();
        real_base l_advNormalR = std::numeric_limits<real_base>::max();

        if( l_exL ) {
          l_advSpeedL  = i_bgPars[l_elL][0].a * i_bgPars[l_elL][0].a;
#if PP_N_DIM > 1
          l_advSpeedL += i_bgPars[l_elL][0].b * i_bgPars[l_elL][0].b;
#endif
#if PP_N_DIM > 2
          l_advSpeedL += i_bgPars[l_elL][0].c * i_bgPars[l_elL][0].c;
#endif
          l_advSpeedL  = std::sqrt( l_advSpeedL );
        }

        if( l_exR ) {
          l_advSpeedR  = i_bgPars[l_elR][0].a * i_bgPars[l_elR][0].a;
#if PP_N_DIM > 1
          l_advSpeedR += i_bgPars[l_elR][0].b * i_bgPars[l_elR][0].b;
#endif
#if PP_N_DIM > 2
          l_advSpeedR += i_bgPars[l_elR][0].c * i_bgPars[l_elR][0].c;
#endif
          l_advSpeedR = std::sqrt( l_advSpeedR );
        }

        // project advection on normal
        if( l_exL ) {
          l_advNormalL  = i_bgPars[l_elL][0].a * i_faceChars[l_fa].outNormal[0];
#if PP_N_DIM > 1
          l_advNormalL += i_bgPars[l_elL][0].b * i_faceChars[l_fa].outNormal[1];
#endif
#if PP_N_DIM > 2
          l_advNormalL += i_bgPars[l_elL][0].c * i_faceChars[l_fa].outNormal[2];
#endif
        }

        if( l_exR ) {
          l_advNormalR  = i_bgPars[l_elR][0].a * i_faceChars[l_fa].outNormal[0];
#if PP_N_DIM > 1
          l_advNormalR += i_bgPars[l_elR][0].b * i_faceChars[l_fa].outNormal[1];
#endif
#if PP_N_DIM > 2
          l_advNormalR += i_bgPars[l_elR][0].c * i_faceChars[l_fa].outNormal[2];
#endif
        }

        // check that our directed normal speeds are lower than the dimensionless advection speed
        assert( l_exL == false || l_advSpeedL + TOL.SOLVER > std::abs(l_advNormalL) );
        assert( l_exR == false || l_advSpeedR + TOL.SOLVER > std::abs(l_advNormalR) );

        // find ids of local contribution flux solvers for the elements
        unsigned short l_fLocIdL = std::numeric_limits<unsigned short>::max();
        unsigned short l_fLocIdR = std::numeric_limits<unsigned short>::max();

        for( unsigned int l_elFa = 0; l_elFa < C_ENT[T_SDISC.ELEMENT].N_FACES; l_elFa++ ) {
          if( l_exL == true && i_elementAdjacentFaces[l_elL][l_elFa] == l_fa ) {
            assert( l_fLocIdL == std::numeric_limits<unsigned short>::max() ); // this must be unique
            l_fLocIdL = l_elFa; // id of local contribution solver
          }
          if( l_exR == true && i_elementAdjacentFaces[l_elR][l_elFa] == l_fa ) {
            assert( l_fLocIdR == std::numeric_limits<unsigned short>::max() ); // this must be unique
            l_fLocIdR = l_elFa; // id of local contribution solver
          }
        }
        assert( l_exL == false || l_fLocIdL != std::numeric_limits<unsigned short>::max() );
        assert( l_exR == false || l_fLocIdR != std::numeric_limits<unsigned short>::max() );

        // find ids of neigh contribution flux solvers for the elements
        unsigned short l_fNegIdL = l_fLocIdL + C_ENT[T_SDISC.ELEMENT].N_FACES;
        unsigned short l_fNegIdR = l_fLocIdR + C_ENT[T_SDISC.ELEMENT].N_FACES;

        /*
         * setup solvers of elements
         *
         * n is the normal -> positive advection speeds go L->R
         *                 -> negative advection speeds go R->L
         */
         // check for normal direction of left element's wave speeds
         if( l_exR == false || (l_exL == true && l_advNormalL > 0) ) {
           assert( l_exL == false || std::abs( o_fluxSolvers[l_elL][l_fLocIdL] ) < TOL.SOLVER );
           if( l_exL ) o_fluxSolvers[l_elL][l_fLocIdL] -= l_advNormalL;

           assert( l_exR == false || std::abs( o_fluxSolvers[l_elR][l_fNegIdR] ) < TOL.SOLVER );
           if( l_exR ) o_fluxSolvers[l_elR][l_fNegIdR] += l_advNormalL;

           if( l_exR == false &&  l_advNormalL < 0 ) {
             o_fluxSolvers[l_elL][l_fLocIdL] *= -1;
           }
         }

         // check for normal direction of right element's wave speeds
         if( l_exL == false || (l_exR == true && l_advNormalR < 0) ) {
           assert( l_exR == false || std::abs( o_fluxSolvers[l_elR][l_fLocIdR] ) < TOL.SOLVER );
           if( l_exR ) o_fluxSolvers[l_elR][l_fLocIdR] += l_advNormalR;

           assert( l_exL == false || std::abs( o_fluxSolvers[l_elL][l_fNegIdL] ) < TOL.SOLVER );
           if( l_exL ) o_fluxSolvers[l_elL][l_fNegIdL] -= l_advNormalR;

           if( l_exL == false && l_advNormalR > 0) {
             o_fluxSolvers[l_elR][l_fNegIdR] *= -1;
           }
         }

        // derive vertex coords
        real_mesh l_veCoordsL[N_DIM][C_ENT[T_SDISC.ELEMENT].N_VERTICES];
        real_mesh l_veCoordsR[N_DIM][C_ENT[T_SDISC.ELEMENT].N_VERTICES];

        for( unsigned int l_ve = 0; l_ve < C_ENT[T_SDISC.ELEMENT].N_VERTICES; l_ve++ ) {
          int_el l_veIdL = std::numeric_limits<int_el>::max();
          int_el l_veIdR = std::numeric_limits<int_el>::max();
          if( l_exL == true ) l_veIdL = i_elVe[l_elL][l_ve];
          if( l_exR == true ) l_veIdR = i_elVe[l_elR][l_ve];

          for( unsigned int l_dim = 0; l_dim < N_DIM; l_dim++ ) {
            if( l_exL == true ) l_veCoordsL[l_dim][l_ve] = i_vertexChars[l_veIdL].coords[l_dim];
            if( l_exR == true ) l_veCoordsR[l_dim][l_ve] = i_vertexChars[l_veIdR].coords[l_dim];
          }
        }

        // compute determinant of the mapping's jacobian
        real_mesh l_jacL[N_DIM][N_DIM], l_jacR[N_DIM][N_DIM];
        for( unsigned short l_dm1 = 0; l_dm1 < N_DIM; l_dm1++ )
          for( unsigned short l_dm2 = 0; l_dm2 < N_DIM; l_dm2++ )
            l_jacL[l_dm1][l_dm2] = l_jacR[l_dm1][l_dm2] = std::numeric_limits<real_mesh>::max();
        real_mesh l_jL = std::numeric_limits<real_mesh>::max();
        real_mesh l_jR = std::numeric_limits<real_mesh>::max();
        if( l_exL == true ) linalg::Mappings::evalJac( T_SDISC.ELEMENT, l_veCoordsL[0], l_jacL[0] );
        if( l_exR == true ) linalg::Mappings::evalJac( T_SDISC.ELEMENT, l_veCoordsR[0], l_jacR[0] );
#if PP_N_DIM == 1
        if( l_exL == true ) l_jL = l_jacL[0][0];
        if( l_exR == true ) l_jR = l_jacR[0][0];
#elif PP_N_DIM == 2
        if( l_exL == true ) l_jL = linalg::Matrix::det2x2( l_jacL );
        if( l_exR == true ) l_jR = linalg::Matrix::det2x2( l_jacR );
#elif PP_N_DIM == 3
        if( l_exL == true ) l_jL = linalg::Matrix::det3x3( l_jacL );
        if( l_exR == true ) l_jR = linalg::Matrix::det3x3( l_jacR );
#else
#error number of dimensions not supported.
#endif

        if( T_SDISC.ELEMENT == TET4 && ORDER == 1 ) {
          // Remark: For tets the determinant is not equal to the volume of the tets.
          //         We apply the effect of the inverse mass matrix right here.
          //         Effectively this just scales the flux solver by the volume of the tets.
          if( l_exL == true ) l_jL /= (real_mesh) 6;
          if( l_exR == true ) l_jR /= (real_mesh) 6;
        }

        // ensure positive determinants
        assert( l_exL == false || l_jL > 0 );  assert( l_exR == false || l_jR > 0 );
        assert( i_faceChars[l_fa].area > 0 );

        // scale the flux solvers by dets of volume integration
        if( l_exL == true ) o_fluxSolvers[l_elL][l_fLocIdL] /= l_jL;
        if( l_exR == true ) o_fluxSolvers[l_elR][l_fLocIdR] /= l_jR;
        if( l_exL == true ) o_fluxSolvers[l_elL][l_fNegIdL] /= l_jL;
        if( l_exR == true ) o_fluxSolvers[l_elR][l_fNegIdR] /= l_jR;

        // scale by the face surfaces
        if( l_exL == true ) o_fluxSolvers[l_elL][l_fLocIdL] *= i_faceChars[l_fa].area;
        if( l_exR == true ) o_fluxSolvers[l_elR][l_fLocIdR] *= i_faceChars[l_fa].area;
        if( l_exL == true ) o_fluxSolvers[l_elL][l_fNegIdL] *= i_faceChars[l_fa].area;
        if( l_exR == true ) o_fluxSolvers[l_elR][l_fNegIdR] *= i_faceChars[l_fa].area;

#if defined PP_T_ELEMENTS_TET4
        // the det of the surface-jac is twice the area of the triangle
        // TODO: Replace this formulation with Jacobi-determinants.
        if( l_exL == true ) o_fluxSolvers[l_elL][l_fLocIdL] *= 2;
        if( l_exR == true ) o_fluxSolvers[l_elR][l_fLocIdR] *= 2;
        if( l_exL == true ) o_fluxSolvers[l_elL][l_fNegIdL] *= 2;
        if( l_exR == true ) o_fluxSolvers[l_elR][l_fNegIdR] *= 2;
#endif
      }

      // take care of duplicated elements
      for( int_el l_el = 0; l_el < i_nElements; l_el++ ) {
        // get dominant id
        int_el l_elDo = i_elMeDa[ i_elDaMe[l_el] ];

        for( int_md l_fa = 0; l_fa < C_ENT[T_SDISC.ELEMENT].N_FACES*2; l_fa++ ) {
          o_fluxSolvers[l_el][l_fa] = o_fluxSolvers[l_elDo][l_fa];
        }
      }
    }
};

#endif
