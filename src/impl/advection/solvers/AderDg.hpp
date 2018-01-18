/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2016-2017, Regents of the University of California
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
 * ADER-DG solver for linear equations (e.g. advection equation: f(q) = a * q || f(q) = a*q + b*q or elastics ).
 **/
#ifndef ADER_DG_HPP
#define ADER_DG_HPP

#include <limits>
#include <cassert>
#include "constants.hpp"
#include "mesh/common.hpp"
#include "linalg/Matrix.h"
#include "linalg/Mappings.hpp"
#include "TimePred.hpp"

namespace edge {
  namespace advection {
    namespace solvers {
      class AderDg;
    }
  }
}

class edge::advection::solvers::AderDg {
  public:
    /**
     * Sets up the star matrices, which are a linear combination of the Jacobians.
     *
     * @param i_nElements number of elements.
     * @param i_vertexChars vertex characteristics.
     * @param i_elVe vertices adjacent to the elements.
     * @param i_bgPars background parameters.
     * @param o_starMatrices will be set to star matrices.
     **/
    static void setupStarM(       int_el           i_nElements,
                            const t_vertexChars   *i_vertexChars,
                            const int_el         (*i_elVe)[ C_ENT[T_SDISC.ELEMENT].N_VERTICES ],
                            const t_bgPars       (*i_bgPars)[1],
                                  real_base      (*o_starMatrices)[N_DIM] ) {
      // iterate over elements
      for( int_el l_el = 0; l_el < i_nElements; l_el++ ) {
        // derive vertex coords
        real_mesh l_veCoords[3][C_ENT[T_SDISC.ELEMENT].N_VERTICES];
        mesh::common< T_SDISC.ELEMENT >::getElVeCoords( l_el, i_elVe, i_vertexChars, l_veCoords );

        // get inverse jacobian
        real_mesh l_jac[N_DIM][N_DIM];
        linalg::Mappings::evalJac( T_SDISC.ELEMENT, l_veCoords[0], l_jac[0] );

        real_mesh l_jacInv[N_DIM][N_DIM];
#if PP_N_DIM == 1
        l_jacInv[0][0] = 1 / l_jac[0][0];
#elif PP_N_DIM == 2
        linalg::Matrix::inv2x2( l_jac, l_jacInv );
#elif PP_N_DIM == 3
        linalg::Matrix::inv3x3( l_jac, l_jacInv );
#else
#error invalid dimension.
#endif

        // set star matrices
        // iterate over reference dimensions
        for( unsigned int l_dim = 0; l_dim < N_DIM; l_dim++ ) {
          o_starMatrices[l_el][l_dim]  = i_bgPars[l_el][0].a * l_jacInv[0][l_dim];
#if PP_N_DIM > 1
          o_starMatrices[l_el][l_dim] += i_bgPars[l_el][0].b * l_jacInv[1][l_dim];
#endif
#if PP_N_DIM > 2
          o_starMatrices[l_el][l_dim] += i_bgPars[l_el][0].c * l_jacInv[2][l_dim];
#endif
         }
      }
    }

    /**
     * Local step: Cauchy Kowalewski + volume.
     *
     * @param i_first first element considered.
     * @param i_nElement number of elements.
     * @param i_dT time step.
     * @param i_dg constant DG data.
     * @param i_starM star matrices.
     * @param i_fluxSolvers flux solvers.
     * @param io_dofs DOFs.
     * @parma o_tInt will be set to time integrated DOFs.
     **/
    static void local(       int_el       i_first,
                             int_el       i_nElements,
                             real_base    i_dT,
                       const t_dg        &i_dg,
                       const real_base (*i_starM)[N_DIM],
                       const real_base (*i_fluxSolvers)[ C_ENT[T_SDISC.ELEMENT].N_FACES*2 ],
                             real_base (*io_dofs)[1][N_ELEMENT_MODES][N_CRUNS],
                             real_base (*o_tInt)[1][N_ELEMENT_MODES][N_CRUNS] ) {
      // iterate over all elements
#ifdef PP_USE_OMP
#pragma omp parallel for
#endif
      for( int_el l_el = i_first; l_el < i_first+i_nElements; l_el++ ) {
        // temporary product for two-way mult
        real_base l_tmpProd[N_ELEMENT_MODES][N_CRUNS];

        /*
         * compute ader time integration
         */
        // buffer for derivatives
        real_base l_derBuffer[ORDER][N_ELEMENT_MODES][N_CRUNS];

#if defined PP_T_KERNELS_VANILLA
        TimePred<
          T_SDISC.ELEMENT,
          ORDER,
          ORDER,
          N_CRUNS >::ckVanilla( i_dT,
                                i_dg.mat.stiffT,
                                i_starM[l_el],
                                io_dofs[l_el][0],
                                l_tmpProd,
                                l_derBuffer,
                                o_tInt[l_el][0] );
#else
#error kernels not implemented.
#endif

        /*
         * compute volume contribution
         */
         for( unsigned int l_di = 0; l_di < N_DIM; l_di++ ) {
#if defined PP_T_KERNELS_VANILLA
            linalg::Matrix::matMulFusedAC( N_CRUNS,                     // #fused
                                           1,                           // m
                                           N_ELEMENT_MODES,             // n
                                           N_ELEMENT_MODES,             // k
                                           N_ELEMENT_MODES,             // ldA
                                           N_ELEMENT_MODES,             // ldB
                                           N_ELEMENT_MODES,             // ldC
                                           static_cast<real_base>(0.0), // beta
                                           o_tInt[l_el][0][0],          // A
                                           i_dg.mat.stiff[l_di][0],     // B
                                           l_tmpProd[0]  );             // C

            linalg::Matrix::matMulFusedBC( N_CRUNS,                     // #fused
                                           1,                           // m
                                           N_ELEMENT_MODES,             // n
                                           1,                           // k
                                           1,                           // ldA
                                           N_ELEMENT_MODES,             // ldB
                                           N_ELEMENT_MODES,             // ldC
                                           static_cast<real_base>(1.0), // beta
                                           i_starM[l_el]+l_di,          // A
                                           l_tmpProd[0],                // B
                                           io_dofs[l_el][0][0] );       // C
#else
           EDGE_LOG_FATAL << "not implemented;"
#endif
         }

         /*
          * compute local surface contribution
          */
         for( unsigned int l_fa = 0; l_fa < C_ENT[T_SDISC.ELEMENT].N_FACES; l_fa++ ) {
#if defined PP_T_KERNELS_VANILLA
           // scratch space for three-way product
           real_base l_scratch[2][N_FACE_MODES][N_CRUNS];

           linalg::Matrix::matMulFusedAC( N_CRUNS,                     // #fused
                                          1,                           // m
                                          N_FACE_MODES,                // n
                                          N_ELEMENT_MODES,             // k
                                          N_ELEMENT_MODES,             // ldA
                                          N_FACE_MODES,                // ldB
                                          N_FACE_MODES,                // ldC
                                          static_cast<real_base>(0.0), // beta
                                          o_tInt[l_el][0][0],          // A
                                          i_dg.mat.fluxL[l_fa][0],     // B
                                          l_scratch[0][0]  );          // C

           linalg::Matrix::matMulFusedBC( N_CRUNS,                     // #fused
                                          1,                           // m
                                          N_FACE_MODES,                // n
                                          1,                           // k
                                          1,                           // ldA
                                          N_FACE_MODES,                // ldB
                                          N_FACE_MODES,                // ldC
                                          static_cast<real_base>(0.0), // beta
                                          i_fluxSolvers[l_el]+l_fa,    // A
                                          l_scratch[0][0],             // B
                                          l_scratch[1][0] );           // C

           linalg::Matrix::matMulFusedAC( N_CRUNS,                     // #fused
                                          1,                           // m
                                          N_ELEMENT_MODES,             // n
                                          N_FACE_MODES,                // k
                                          N_FACE_MODES,                // ldA
                                          N_ELEMENT_MODES,             // ldB
                                          N_ELEMENT_MODES,             // ldC
                                          static_cast<real_base>(1.0), // beta
                                          l_scratch[1][0],             // A
                                          i_dg.mat.fluxT[l_fa][0],     // B
                                          io_dofs[l_el][0][0]  );      // C
#else
           EDGE_LOG_FATAL << "not implemented;"
#endif
         }
      }
    }

    /**
     * Performs the neighboring updates of the ADER-DG scheme.
     *
     * @param i_first first element considered.
     * @param i_nElements number of elements.
     * @param i_dT time step.
     * @param i_dg constant DG data.
     * @parma i_faChars face characteristics.
     * @param i_fluxSolvers flux solvers.
     * @param i_elFa elements' adjacent faces.
     * @param i_elFaEl elements' adjacent elements (through faces).
     * @param i_fIdElFaEl local face ids of face-neighboring elememts.
     * @param i_vIdElFaEl local vertex ids w.r.t. the shared face from the neighboring elements' perspective.
     * @param i_tInt time integrated degrees of freedom.
     * @param io_dofs DOFs which will be updated with neighboring elements' contribution.
     **/
    static void neigh(       int_el           i_first,
                             int_el           i_nElements,
                             real_base        i_dT,
                       const t_dg            &i_dg,
                       const t_faceChars     *i_faChars,
                       const real_base      (*i_fluxSolvers)[ C_ENT[T_SDISC.ELEMENT].N_FACES*2 ],
                       const int_el         (*i_elFa)[C_ENT[T_SDISC.ELEMENT].N_FACES],
                       const int_el         (*i_elFaEl)[C_ENT[T_SDISC.ELEMENT].N_FACES],
                       const unsigned short (*i_fIdElFaEl)[C_ENT[T_SDISC.ELEMENT].N_FACES],
                       const unsigned short (*i_vIdElFaEl)[C_ENT[T_SDISC.ELEMENT].N_FACES],
                       const real_base      (*i_tInt)[1][N_ELEMENT_MODES][N_CRUNS],
                             real_base      (*io_dofs)[1][N_ELEMENT_MODES][N_CRUNS] ) {
      // iterate over elements
#ifdef PP_USE_OMP
#pragma omp parallel for
#endif
      for( int_el l_el = i_first; l_el < i_first+i_nElements; l_el++ ) {
        // add neighboring contribution
        for( int_el l_fa = 0; l_fa < C_ENT[T_SDISC.ELEMENT].N_FACES; l_fa++ ) {
          int_el l_faId = i_elFa[l_el][l_fa];
          int_el l_ne;

          if( (i_faChars[l_faId].spType & OUTFLOW) != OUTFLOW ) {
            // derive neighbor
            l_ne = i_elFaEl[l_el][l_fa];
          }
          // outflow boundary conditions
          else {
            l_ne = l_el;
          }

          unsigned short l_fId =  i_vIdElFaEl[l_el][l_fa] * C_ENT[T_SDISC.ELEMENT].N_FACES;
                         l_fId += i_fIdElFaEl[l_el][l_fa];

          // scratch space for three-way product
          real_base l_scratch[2][N_FACE_MODES][N_CRUNS];

          linalg::Matrix::matMulFusedAC( N_CRUNS,                                                                          // #fused
                                         1,                                                                                // m
                                         N_FACE_MODES,                                                                     // n
                                         N_ELEMENT_MODES,                                                                  // k
                                         N_ELEMENT_MODES,                                                                  // ldA
                                         N_FACE_MODES,                                                                     // ldB
                                         N_FACE_MODES,                                                                     // ldC
                                         static_cast<real_base>(0.0),                                                      // beta
                                         i_tInt[l_ne][0][0],                                                               // A
                                         ( (i_faChars[l_faId].spType & OUTFLOW) != OUTFLOW ) ? i_dg.mat.fluxN[l_fId][0] :
                                                                                               i_dg.mat.fluxL[l_fa][0],    // B
                                         l_scratch[0][0]  );                                                               // C

          linalg::Matrix::matMulFusedBC( N_CRUNS,                           // #fused
                                         1,                                 // m
                                         N_FACE_MODES,                      // n
                                         1,                                 // k
                                         1,                                 // ldA
                                         N_FACE_MODES,                      // ldB
                                         N_FACE_MODES,                      // ldC
                                         static_cast<real_base>(0.0),       // beta
                                         i_fluxSolvers[l_el] +
                                           C_ENT[T_SDISC.ELEMENT].N_FACES +
                                           l_fa,                            // A
                                         l_scratch[0][0],                   // B
                                         l_scratch[1][0] );                 // C

          linalg::Matrix::matMulFusedAC( N_CRUNS,                           // #fused
                                         1,                                 // m
                                         N_ELEMENT_MODES,                   // n
                                         N_FACE_MODES,                      // k
                                         N_FACE_MODES,                      // ldA
                                         N_ELEMENT_MODES,                   // ldB
                                         N_ELEMENT_MODES,                   // ldC
                                         static_cast<real_base>(1.0),       // beta
                                         l_scratch[1][0],                   // A
                                         i_dg.mat.fluxT[l_fa][0],           // B
                                         io_dofs[l_el][0][0]  );            // C
        }
      }
    }
};

#endif
