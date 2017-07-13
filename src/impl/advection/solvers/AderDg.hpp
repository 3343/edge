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
  //private:
    /**
     * Performs the operation C = A.B with private per-run data in C and B.
     *
     * @param i_a vector A.
     * @param i_b matrix B.
     * @param o_c vector C, will bet set to A.B.
     **/
    static void matMulB0( const real_base i_a[N_ELEMENT_MODES][N_CRUNS],
                          const real_base i_b[N_ELEMENT_MODES][N_ELEMENT_MODES],
                                real_base o_c[N_ELEMENT_MODES][N_CRUNS] ) {
#if __has_builtin(__builtin_assume_aligned)
      // share alignment with compiler
      (void) __builtin_assume_aligned(i_a, ALIGNMENT.ELEMENT_MODES.PRIVATE);
      (void) __builtin_assume_aligned(o_c, ALIGNMENT.ELEMENT_MODES.PRIVATE);
#endif

      // reset result to zero
      for( int_md l_md = 0; l_md < N_ELEMENT_MODES; l_md++ ) {
        for( int_cfr l_cfr = 0; l_cfr < N_CRUNS; l_cfr++ ) {
          o_c[l_md][l_cfr] = 0;
        }
      }

      for( unsigned int l_i = 0; l_i < N_ELEMENT_MODES; l_i++ ) {
        for( unsigned int l_j = 0; l_j < N_ELEMENT_MODES; l_j++ ) {
          for( int_cfr l_cfr = 0; l_cfr < N_CRUNS; l_cfr++ ) {
            o_c[l_j][l_cfr] += i_a[l_i][l_cfr] * i_b[l_i][l_j];
          }
        }
      }
    }

    /**
     * Performs the operation C += A.B with private per-run data in C and B.
     *
     * @param i_a scalar A.
     * @param i_b vector B.
     * @param i_c vector C.
     **/
    static void matMulB1( const real_base i_a,
                          const real_base i_b[N_ELEMENT_MODES][N_CRUNS],
                                real_base o_c[N_ELEMENT_MODES][N_CRUNS] ) {
#if __has_builtin(__builtin_assume_aligned)
      // share alignment with compiler
      (void) __builtin_assume_aligned(i_b, ALIGNMENT.ELEMENT_MODES.PRIVATE);
      (void) __builtin_assume_aligned(o_c, ALIGNMENT.ELEMENT_MODES.PRIVATE);
#endif

      for( int_md l_md = 0; l_md < N_ELEMENT_MODES; l_md++ ) {
        for( int_cfr l_cfr = 0; l_cfr < N_CRUNS; l_cfr++ ) {
          o_c[l_md][l_cfr] += i_a * i_b[l_md][l_cfr];
        }
      }
    }

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
                             double       i_dT,
                       const t_dg        &i_dg,
                       const real_base (*i_starM)[N_DIM],
                       const real_base (*i_fluxSolvers)[ C_ENT[T_SDISC.ELEMENT].N_FACES*2 ],
                             real_base (*io_dofs)[1][N_ELEMENT_MODES][N_CRUNS],
                             real_base (*o_tInt)[1][N_ELEMENT_MODES][N_CRUNS] ) {
#if __has_builtin(__builtin_assume_aligned)
      // share alignment with compiler
      (void) __builtin_assume_aligned(io_dofs, ALIGNMENT.ELEMENT_MODES.PRIVATE);
      (void) __builtin_assume_aligned(o_tInt,  ALIGNMENT.ELEMENT_MODES.PRIVATE);
#endif

      // iterate over all elements
#ifdef PP_USE_OMP
#pragma omp parallel for
#endif
      for( int_el l_el = i_first; l_el < i_first+i_nElements; l_el++ ) {
        // temporary product for two-way mult
        real_base l_tmpProd[N_ELEMENT_MODES][N_CRUNS]  __attribute__ ((aligned (ALIGNMENT.BASE.STACK)));

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
         for( unsigned int l_dim = 0; l_dim < N_DIM; l_dim++ ) {
#if defined PP_T_KERNELS_VANILLA
           // multiply with stiffness and inverse mass matrix
           matMulB0( o_tInt[l_el][0], i_dg.mat.stiff[l_dim], l_tmpProd );

           // multiply with star "matrix"
           matMulB1( i_starM[l_el][l_dim], l_tmpProd, io_dofs[l_el][0] );
#else
           EDGE_LOG_FATAL << "not implemented;"
#endif
         }

         /*
          * compute local surface contribution
          */
         for( unsigned int l_fa = 0; l_fa < C_ENT[T_SDISC.ELEMENT].N_FACES; l_fa++ ) {
#if defined PP_T_KERNELS_VANILLA
           // multiply with flux matrix
           matMulB0( o_tInt[l_el][0], i_dg.mat.flux[l_fa], l_tmpProd );

           // multiply with flux solver
           matMulB1( i_fluxSolvers[l_el][l_fa], l_tmpProd, io_dofs[l_el][0] );
#else
           EDGE_LOG_FATAL << "not implemented;"
#endif
         }
      }
    }

    /**
     * Quadrature-free computation of the surface intergral (neighboring contribution) for the given face.
     *
     * @param i_fa element local id of the face.
     * @param i_fIdElFaEl local face id of the face-neighboring element.
     * @param i_vIdElFaEl local vertex id w.r.t. the shared face from the neighboring element's perspective.
     * @param i_dgMat DG matrices.
     * @param i_faChars characteristics of the face.
     * @param i_fluxSolvers flux solvers of the face.
     * @param i_tIntNe time integrated DOFs of the adjacent element.
     * @param io_dofs will be updated with neighboring contribution of the face to the surface integral.
     */
    static void surfIntNeQuadFree( unsigned short        i_fa,
                                   unsigned short        i_fIdElFaEl,
                                   unsigned short        i_vIdElFaEl,
                                   t_dgMat        const &i_dgMat,
                                   t_faceChars    const &i_faChars,
                                   real_base      const  i_fluxSolvers[ C_ENT[T_SDISC.ELEMENT].N_FACES*2 ],
                                   real_base      const  i_tIntNe[N_ELEMENT_MODES][N_CRUNS],
                                   real_base             io_dofs[N_ELEMENT_MODES][N_CRUNS] ) {
#if __has_builtin(__builtin_assume_aligned)
      // share alignment with compiler
      (void) __builtin_assume_aligned(i_tIntNe, ALIGNMENT.ELEMENT_MODES.PRIVATE);
      (void) __builtin_assume_aligned(io_dofs,  ALIGNMENT.ELEMENT_MODES.PRIVATE);
#endif

      // temporary product for two-way mult
      real_base (*l_tmpProd)[N_CRUNS] = parallel::g_scratchMem->tRes;

      // id of the flux matrix
      unsigned short l_fId;

      if( (i_faChars.spType & OUTFLOW) != OUTFLOW ) {
        /*
         * derive id of flux matrix
         */
  #if defined PP_T_ELEMENTS_HEX8R
        // shortcut for rectangular, 8-node hexes, having only six flux matrices
        l_fId  = C_ENT[T_SDISC.ELEMENT].N_FACES + i_fa;
  #else
        // only jump over vertex combinations for 3D elements
        unsigned short l_vertexJump = std::max( (N_DIM / 3) * C_ENT[T_SDISC.FACE].N_VERTICES, 1);

        // jump over local flux matrices
        l_fId  = C_ENT[T_SDISC.ELEMENT].N_FACES;
        // jump over local face
        l_fId += i_fa * C_ENT[T_SDISC.ELEMENT].N_FACES * l_vertexJump;

        // jump over neighboring face
        l_fId += i_fIdElFaEl * l_vertexJump;

        // jump over vertices
        l_fId += i_vIdElFaEl;
  #endif
      }
      // outflow boundary conditions
      else {
        l_fId = i_fa;
      }

      /*
       * solve
       */
  #if defined PP_T_KERNELS_VANILLA
      // multiply with flux matrix
      matMulB0( i_tIntNe, i_dgMat.flux[l_fId], l_tmpProd );

      // multiply with flux solver
      matMulB1( i_fluxSolvers[ C_ENT[T_SDISC.ELEMENT].N_FACES + i_fa ], l_tmpProd, io_dofs );
  #else
      EDGE_LOG_FATAL << "not implemented;"
  #endif
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
                             double           i_dT,
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

          // compute quadrature free neighboring contribution to the face's surface integral
          surfIntNeQuadFree( l_fa,
                             i_fIdElFaEl[l_el][l_fa],
                             i_vIdElFaEl[l_el][l_fa],
                             i_dg.mat,
                             i_faChars[l_faId],
                             i_fluxSolvers[l_el],
                             i_tInt[l_ne][0],
                             io_dofs[l_el][0] );
        }
      }
    }
};

#endif
