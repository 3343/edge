/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
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
 * Support for internal boundary conditions through quadrature rules.
 **/

#ifndef INTERNAL_BOUNDARY_HPP
#define INTERNAL_BOUNDARY_HPP

#include "InternalBoundary.type"
#include "constants.hpp"
#include "io/logging.h"
#include "dg/QuadratureEval.hpp"
#include "linalg/Matrix.h"
#include "TimePred.hpp"
#include "common.hpp"

namespace edge {
  namespace elastic {
    namespace solvers {
      template< t_entityType TL_T_EL, unsigned short TL_N_QU, unsigned short TL_O_SP, unsigned short TL_N_CRUNS, unsigned short TL_O_TI >
      class InternalBoundary;

      template< t_entityType TL_T_EL >
      class InternalBoundaryTypes;

#ifdef __INTEL_COMPILER
      template< t_entityType TL_T_EL, unsigned short TL_N_DIM=N_DIM >
#else
      template< t_entityType TL_T_EL, unsigned short TL_N_DIM=C_ENT[TL_T_EL].N_DIM >
#endif
      class InternalBoundarySolvers;

      template< t_entityType TL_T_EL >
      class InternalBoundarySolvers< TL_T_EL, 2 >;

      template< t_entityType TL_T_EL >
      class InternalBoundarySolvers< TL_T_EL, 3 >;
    }
  }
}

/**
 * Internal boundary conditions through quadrature rules.
 *
 * @paramt TL_T_EL element type.
 * @paramt TL_N_QU number of quantities.
 * @paramt TL_O_SP order of the used quadrature in space.
 * @paramt TL_N_CRUNS number of concurrent forward runs (fused simulations).
 * @paramt TL_O_TI order of the used quadrature in time.
 **/
template <t_entityType TL_T_EL, unsigned short TL_N_QU, unsigned short TL_O_SP, unsigned short TL_N_CRUNS, unsigned short TL_O_TI=TL_O_SP>
class edge::elastic::solvers::InternalBoundary {
  private:
    // assemble derived template parameters
    //! dimension of the element
    static unsigned short const TL_N_DIM              = C_ENT[TL_T_EL].N_DIM;
    //! number of element modes
    static unsigned short const TL_N_ELEMENT_MODES    = CE_N_ELEMENT_MODES( TL_T_EL, TL_O_SP );
    //! number of options for the layout of quad points
    static unsigned short const TL_N_FACE_QUAD_OPTS   = (CE_N_FACE_VERTEX_OPTS(TL_T_EL)+1) *
                                                         C_ENT[TL_T_EL].N_FACES;
    //! number of quadrature points per face
    static unsigned short const TL_N_FACE_QUAD_POINTS = CE_N_FACE_QUAD_POINTS( TL_T_EL, TL_O_SP );

  public:
    /**
     * Dummy solver perturbating nothing.
     *
     * @param i_ms middle state.
     * @param o_msL will be set to middle state.
     * @param o_msR will be set to middle state.
     *
     * @paramt TL_T_REAL type of floating point arithmetic.
     * @paramt TL_T_FA_DATA face data.
     **/
    class DummySolv {
      public:
        /**
         * Dummy perturbations.
         *
         **/
        template< typename TL_T_REAL, typename TL_T_FA_DATA >
        static void inline perturb( unsigned short,
                                    TL_T_REAL,
                                    TL_T_REAL            i_ms[TL_N_QU][TL_N_CRUNS],
                                    TL_T_FA_DATA const *,
                                    TL_T_REAL            o_msL[TL_N_QU][TL_N_CRUNS],
                                    TL_T_REAL            o_msR[TL_N_QU][TL_N_CRUNS] ) {
          for( unsigned short l_qt = 0; l_qt < TL_N_QU; l_qt++ ) {
            for( unsigned short l_ru = 0; l_ru < TL_N_CRUNS; l_ru++ ) {
              o_msL[l_qt][l_ru] = i_ms[l_qt][l_ru];
              o_msR[l_qt][l_ru] = i_ms[l_qt][l_ru];
            }
          }
        }
    };

    /**
     * Evaluates the internal boundary condition in space at a face of the given element type.
     *
     * Tria3-example:
     *
     *
     * 1) Illustatration of an element in physical coordinates
     *
     *   x: face-local quadrature points
     *
     *   vertices to faces: f0: 0-1, f1: 1-2, f2: 2-0
     *
     *              *
     *             *2* 1  *
     *            *   x         *
     *           *     *           0   *
     *          *0  L   *     R     *
     *             *     x        *
     *                *   *     *
     *                   *1*2 *
     *                      * 
     *   In this example the left element's face f1 matches the right elements face f2.
     *   Vertex 2 of the right element lies on the first vertex of the left elements face f1.
     *
     * 2) Illustration of the quadrature point-local Riemann problem.
     *
     *   Q^L: left-side middle state     Q^R: right-side middle state
     *   c^L_s: left going s-wave        c^L_s right going s-wave
     *   c^L_p: left going p-wave        c^R_p right going p-wave
     *
     *    c^L_p       c^L_s         /|\     c^R_s            c^R_p
     *     *            *          L | R    *                 *
     *         *          *          |     *              *
     *             *        *     Q^L|Q^R *           *
     *                 *      *      |   *        *
     *                     *    *    |t *     *
     *                         *  *  | *  *
     *                   ___________*|*___________\
     *                             (0,0)       x  /
     *
     *   Remark: In general Q^L != Q^R due to stationary waves at x=0.
     *           However the non-equal quantities don't contribute to the fluxes
     *           due to non-zero zero wave speeds of the respective eigenvectors.
     *
     * 3) After solving the Riemann problem, the flux function is applied to the middle states
     *    and the repspective DOF-update stemming from the side obtained through multiplication
     *    with the test functions and quadrature.
     *
     * @paramt TL_T_REAL precision of the evaluation.
     * @param i_faIdL local face id of the left element.
     * @param i_faIdR local face id of the right element.
     * @param i_vIdR local id of the right element's vertex lying on the left element's first face-vertex.
     * @param i_massI diagonal of the inverse mass matrix (orthogonal basis is assumed).
     * @param i_weightsFaces weights of the face's quadrature point.
     * @param i_basisFaces evaluated basis at the quad points.
     *                     [*][][]: options of the quad point layout,
     *                     [][*][]: quad points of the option,
     *                     [][][*]: evaluated basis functions per quad point.
     * @param i_tm1 transformation matrix from physical coordinates for face-aligned coordinates.
     * @param i_solMsJumpL solver for the single jump from the left element's quantities to the middle state.
     * @param i_solMsFluxL flux solver using (probably perturbed) middle states for the left element.
     * @param i_solMsFluxR flux solver using (probably perturbed) middle states for the right element.
     * @param i_dofsL modal DOFs of the left element.
     * @param i_dofsR modal DOFs of the right element.
     * @param o_surfUpdateL will be set to left-going surface update of this part of the internal boundary.
     * @param o_surfUpdateR will be set to right-going surface update of this part of the internal boundary.
     * @param i_dt associated "time step" of this evaluation, might be used internally to compute slip from the slip rate, for example.
     * @param io_faData data used in the pertubation of the middle states.
     *
     * @paramt TL_T_REAL floating point type.
     * @paramt TL_T_MS_SOLV middle state "solver", offers member function .perturb.
     * @paramt TL_T_FA_DATA data passed to middle state solver.
     **/
    template< typename TL_T_REAL,
              typename TL_T_MS_SOLV = DummySolv,
              typename TL_T_FA_DATA = void >
    static void evalSpace( unsigned short        i_faIdL,
                           unsigned short        i_faIdR,
                           unsigned short        i_vIdR,
                           TL_T_REAL      const  i_massI[TL_N_ELEMENT_MODES],
                           TL_T_REAL      const  i_weightsFaces[TL_N_FACE_QUAD_POINTS],
                           TL_T_REAL      const  i_basisFaces[TL_N_FACE_QUAD_OPTS][TL_N_FACE_QUAD_POINTS][TL_N_ELEMENT_MODES],
                           TL_T_REAL      const  i_tm1[TL_N_QU][TL_N_QU],
                           TL_T_REAL      const  i_solMsJumpL[TL_N_QU][TL_N_QU],
                           TL_T_REAL      const  i_solMsFluxL[TL_N_QU][TL_N_QU],
                           TL_T_REAL      const  i_solMsFluxR[TL_N_QU][TL_N_QU],
                           TL_T_REAL      const  i_dofsL[TL_N_QU][TL_N_ELEMENT_MODES][TL_N_CRUNS],
                           TL_T_REAL      const  i_dofsR[TL_N_QU][TL_N_ELEMENT_MODES][TL_N_CRUNS],
                           TL_T_REAL             o_surfUpdateL[TL_N_QU][TL_N_ELEMENT_MODES][TL_N_CRUNS],
                           TL_T_REAL             o_surfUpdateR[TL_N_QU][TL_N_ELEMENT_MODES][TL_N_CRUNS],
                           TL_T_REAL             i_dt = 0,
                           TL_T_FA_DATA         *io_faData = nullptr ) {
      // temporary storage for the middle states
      TL_T_REAL l_msTmp[2][TL_N_QU][TL_N_ELEMENT_MODES][TL_N_CRUNS];
      for( int_qt l_qt = 0; l_qt < TL_N_QU; l_qt++ ) {
        for( int_md l_md = 0; l_md < TL_N_ELEMENT_MODES; l_md++ ) {
          for( int_cfr l_ru = 0; l_ru < TL_N_CRUNS; l_ru++ ) {
            l_msTmp[0][l_qt][l_md][l_ru] = 0;
            l_msTmp[1][l_qt][l_md][l_ru] = 0;
          }
        }
      }

      // rotate the DOFs from physical coordinates to face-aligned coords
      // remark: the back-rotation to physical coordinates is part of the the flux solver
      TL_T_REAL l_dofs[2][TL_N_QU][TL_N_ELEMENT_MODES][TL_N_CRUNS];
      linalg::Matrix::matMulB0FusedBC( TL_N_CRUNS,
                                       TL_N_QU, TL_N_ELEMENT_MODES, TL_N_QU,
                                       i_tm1[0], i_dofsL[0][0], l_dofs[0][0][0] );
      linalg::Matrix::matMulB0FusedBC( TL_N_CRUNS,
                                       TL_N_QU, TL_N_ELEMENT_MODES, TL_N_QU,
                                       i_tm1[0], i_dofsR[0][0], l_dofs[1][0][0] );

      // derive face quad pos of right element
      unsigned short l_posR  = C_ENT[TL_T_EL].N_FACES;
      l_posR                += i_faIdR * CE_N_FACE_VERTEX_OPTS(TL_T_EL);
      l_posR                += i_vIdR;

      // iterate over the quad points in space
      for( int_md l_qp = 0; l_qp < CE_N_FACE_QUAD_POINTS( TL_T_EL, TL_O_SP ); l_qp++ ) {
        // temporary values at the quad points
        TL_T_REAL l_qEv[2][TL_N_QU][TL_N_CRUNS];
        // jump in quantities
        TL_T_REAL l_qJump[TL_N_QU][TL_N_CRUNS];

        // eval left and right elements' DOFs at quad points
        for( int_qt l_qt = 0; l_qt < TL_N_QU; l_qt++ ) {
          dg::QuadratureEval<TL_T_EL, TL_O_SP, TL_N_CRUNS>::evalBasis( i_basisFaces[i_faIdL][l_qp],
                                                                       l_dofs[0][l_qt],
                                                                       l_qEv[0][l_qt] );

          dg::QuadratureEval<TL_T_EL, TL_O_SP, TL_N_CRUNS>::evalBasis( i_basisFaces[l_posR][l_qp],
                                                                       l_dofs[1][l_qt],
                                                                       l_qEv[1][l_qt] );

          // compute the jump in quantities
          for( int_cfr l_ru = 0; l_ru < TL_N_CRUNS; l_ru++ ) {
            l_qJump[l_qt][l_ru] = l_qEv[1][l_qt][l_ru] - l_qEv[0][l_qt][l_ru];
          }
        }
        // jump over waves with negative speeds from the left to get the left-side middle state
        linalg::Matrix::matMulB1FusedBC( TL_N_CRUNS,
                                         TL_N_QU, 1, TL_N_QU,
                                         i_solMsJumpL[0],
                                         l_qJump[0],
                                         l_qEv[0][0] );

         // perturb if necessary
         TL_T_REAL l_ms[2][TL_N_QU][TL_N_CRUNS];
         TL_T_MS_SOLV::perturb( l_qp,
                                i_dt,
                                l_qEv[0],
                                io_faData,
                                l_ms[0],
                                l_ms[1] );

         // compute the contribution of this quad point
         // remark: due to linear fluxes, the flux computation is applied at the very end.
         for( int_qt l_qt = 0; l_qt < TL_N_QU; l_qt++ ) {
           for( int_md l_md = 0; l_md < TL_N_ELEMENT_MODES; l_md++ ) {
             // precompute weights
             TL_T_REAL l_weightL = i_weightsFaces[l_qp] *              // weight of the face
                                   i_basisFaces[i_faIdL][l_qp][l_md] * // test function
                                   i_massI[l_md];                      // inverse mass matrix

             TL_T_REAL l_weightR = i_weightsFaces[l_qp] *              // weight of the face
                                   i_basisFaces[l_posR][l_qp][l_md] *  // test function
                                   i_massI[l_md];                      // inverse mass matrix

             // add contribution
             for( int_cfr l_ru = 0; l_ru < TL_N_CRUNS; l_ru++ ) {
               // left-going fluxes are subtracted
               l_msTmp[0][l_qt][l_md][l_ru] -= l_weightL * l_ms[0][l_qt][l_ru];
               // right-going fluxes are added
               l_msTmp[1][l_qt][l_md][l_ru] += l_weightR * l_ms[1][l_qt][l_ru];
             }
           }
         }
      }

      // compute fluxes and rotate DOFs back to physical coordinate system
      linalg::Matrix::matMulB0FusedBC( TL_N_CRUNS,
                                       TL_N_QU, TL_N_ELEMENT_MODES, TL_N_QU,
                                       i_solMsFluxL[0],
                                       l_msTmp[0][0][0],
                                       o_surfUpdateL[0][0] );

      linalg::Matrix::matMulB0FusedBC( TL_N_CRUNS,
                                       TL_N_QU, TL_N_ELEMENT_MODES, TL_N_QU,
                                       i_solMsFluxR[0],
                                       l_msTmp[1][0][0],
                                       o_surfUpdateR[0][0] );
    }

    /**
     * Evaluates the internal boundary condition in space and time at a face of the given element type.
     *
     * @paramt TL_T_REAL precision of the evaluation.
     * @param i_faIdL local face id of the left element.
     * @param i_faIdR local face id of the right element.
     * @param i_veIdR local id of the right element's vertex lying on the left element's first face-vertex.
     * @param i_massI diagonal of the inverse mass matrix (orthogonal basis is assumed).
     * @param i_dT time step.
     * @param i_ptsLine quadrature points for the unit line element [0,1].
     * @param i_weightsLine quadrature weights for the unit line element [0,1].
     * @param i_weightsFaces weights of the face's quadrature point.
     * @param i_basisFaces evaluated basis at the quad points.
     *                     [*][][]: options of the quad point layout,
     *                     [][*][]: quad points of the option,
     *                     [][][*]: evaluated basis functions per quad point.
     * @param i_tm1 transformation matrix from physical coordinates for face-aligned coordinates.
     * @param i_solMsJumpL solver for the single jump from the left element's quantities to the middle state.
     * @param i_solMsFluxL flux solver using (probably perturbed) middle states for the left element.
     * @param i_solMsFluxR flux solver using (probably perturbed) middle states for the right element.
     * @param i_tDersL modal time derivatives of the left element (time prediction).
     * @param i_tDersR modal time derivatives of the right element (time prediction).
     * @param o_scratch will be used as scratch memory.
     * @param o_surfUpdateL will be set to left-going surface update of this part of the internal boundary.
     * @param o_surfUpdateR will be set to right-going surface update of this part of the internal boundary.
     *
     * @paramt TL_T_REAL floating point type.
     * @paramt TL_T_MS_SOLV middle state "solver", offers member functions .perturb.
     * @paramt TL_T_FA_DATA data passed to middle state solver.
     **/
    template< typename TL_T_REAL,
              typename TL_T_MS_SOLV = DummySolv,
              typename TL_T_FA_DATA = void >
    static void evalSpaceTime( unsigned short       i_faIdL,
                               unsigned short       i_faIdR,
                               unsigned short       i_veIdR,
                               TL_T_REAL      const i_massI[TL_N_ELEMENT_MODES],
                               TL_T_REAL      const i_dT,
                               TL_T_REAL      const i_ptsLine[ TL_O_TI ],
                               TL_T_REAL      const i_weightsLine[ TL_O_TI ],
                               TL_T_REAL      const i_weightsFaces[TL_N_FACE_QUAD_POINTS],
                               TL_T_REAL      const i_basisFaces[TL_N_FACE_QUAD_OPTS][TL_N_FACE_QUAD_POINTS][TL_N_ELEMENT_MODES],
                               TL_T_REAL      const i_tm1[TL_N_QU][TL_N_QU],
                               TL_T_REAL      const i_solMsJumpL[TL_N_QU][TL_N_QU],
                               TL_T_REAL      const i_solMsFluxL[TL_N_QU][TL_N_QU],
                               TL_T_REAL      const i_solMsFluxR[TL_N_QU][TL_N_QU],
                               TL_T_REAL      const i_tDersL[TL_O_SP][TL_N_QU][TL_N_ELEMENT_MODES][TL_N_CRUNS],
                               TL_T_REAL      const i_tDersR[TL_O_SP][TL_N_QU][TL_N_ELEMENT_MODES][TL_N_CRUNS],
                               TL_T_REAL            o_scratch[4][TL_N_QU][TL_N_ELEMENT_MODES][TL_N_CRUNS],
                               TL_T_REAL            o_surfUpdateL[TL_N_QU][TL_N_ELEMENT_MODES][TL_N_CRUNS],
                               TL_T_REAL            o_surfUpdateR[TL_N_QU][TL_N_ELEMENT_MODES][TL_N_CRUNS],
                               TL_T_FA_DATA        *io_faData = nullptr ) {
      // reset updates
      for( int_qt l_qt = 0; l_qt < TL_N_QU; l_qt++ ) {
        for( int_md l_md = 0; l_md < TL_N_ELEMENT_MODES; l_md++ ) {
          for( int_cfr l_ru = 0; l_ru < TL_N_CRUNS; l_ru++ ) {
            o_surfUpdateL[l_qt][l_md][l_ru] = 0;
            o_surfUpdateR[l_qt][l_md][l_ru] = 0;
          }
        }
      }

      // assign pointers to scratch memory
      TL_T_REAL (*l_dofsL  )[TL_N_ELEMENT_MODES][TL_N_CRUNS] = o_scratch[0];
      TL_T_REAL (*l_dofsR  )[TL_N_ELEMENT_MODES][TL_N_CRUNS] = o_scratch[1];
      TL_T_REAL (*l_surfUpL)[TL_N_ELEMENT_MODES][TL_N_CRUNS] = o_scratch[2];
      TL_T_REAL (*l_surfUpR)[TL_N_ELEMENT_MODES][TL_N_CRUNS] = o_scratch[3];

      // iterate over quad points in time
      for( unsigned short l_qp = 0; l_qp < TL_O_TI; l_qp++ ) {
        TL_T_REAL l_ptTime[1];
        l_ptTime[0] = i_ptsLine[l_qp] * i_dT;
        // evaluate the time prediction at the quad point in time
        TimePred< TL_T_EL,
                  TL_N_QU,
                  TL_O_SP,
                  TL_N_CRUNS, 1 >::evalTimePrediction(                       l_ptTime,
                                                                             i_tDersL,
                    (TL_T_REAL (*)[TL_N_QU][TL_N_ELEMENT_MODES][TL_N_CRUNS]) l_dofsL );

        TimePred< TL_T_EL,
                  TL_N_QU,
                  TL_O_SP,
                  TL_N_CRUNS, 1 >::evalTimePrediction(                       l_ptTime,
                                                                             i_tDersR,
                    (TL_T_REAL (*)[TL_N_QU][TL_N_ELEMENT_MODES][TL_N_CRUNS]) l_dofsR );

        // get the "time step" of the spatial integration (used for internal middle state pertubations)
        TL_T_REAL l_dtPt = ( TL_O_TI == 1      ) ? 1                                   : // [0,          qp,            1]
                           ( l_qp == 0         ) ? i_ptsLine[l_qp]                     : // [0, qp, x, x,   [...],      1]
                           ( l_qp < TL_O_TI-1  ) ? i_ptsLine[l_qp] - i_ptsLine[l_qp-1] : // [0, x, [...], qp, x, [...], 1]
                                                   1               - i_ptsLine[l_qp-1];  // [0, x,     [...]     x, qp, 1]
        l_dtPt *= i_dT;

        // perform quadrature in space
        evalSpace<
          TL_T_REAL,
          TL_T_MS_SOLV,
          TL_T_FA_DATA >( i_faIdL, i_faIdR, i_veIdR,
                          i_massI, i_weightsFaces, i_basisFaces,
                          i_tm1,
                          i_solMsJumpL,
                          i_solMsFluxL, i_solMsFluxR,
                          l_dofsL, l_dofsR,
                          l_surfUpL, l_surfUpR,
                          l_dtPt,
                          io_faData );

        // add the the contribution of this temporal quad point to the update
        TL_T_REAL l_scale = i_weightsLine[l_qp] *  i_dT; // scale weights (based on [0,1])
        for( int_qt l_qt = 0; l_qt < TL_N_QU; l_qt++ ) {
          for( int_md l_md = 0; l_md < TL_N_ELEMENT_MODES; l_md++ ) {
            for( int_cfr l_ru = 0; l_ru < TL_N_CRUNS; l_ru++ ) {
              o_surfUpdateL[l_qt][l_md][l_ru] += l_surfUpL[l_qt][l_md][l_ru] * l_scale;
              o_surfUpdateR[l_qt][l_md][l_ru] += l_surfUpR[l_qt][l_md][l_ru] * l_scale;
            }
          }
        }
      }

    }
};

template< t_entityType TL_T_EL >
class edge::elastic::solvers::InternalBoundaryTypes {
  private:
    static unsigned short const TL_N_EL_FA = C_ENT[T_SDISC.ELEMENT].N_FACES;

  public:
    /**
     * Initializes the data of an internal boundary faces.
     *
     * @param i_nFaDe number of dense faces.
     * @param i_spType sparse type of the internal boundary.
     * @param i_charsFa characteristics of the dense faces.
     * @param i_faEl elements adjacent to the dense faces.
     * @param i_elFa dense faces adjacent to the dense elements.
     * @param i_vIdElFaEl vertex ids of the shared face with respect to the dense element's adjacent dense elements.
     *
     * @paramt TL_T_INT_LID integer type of local ids.
     * @paramt TL_T_REAL real type used in arithmetic operations.
     * @paramt TL_T_INT_SP type of the sparse type.
     * @paramt TL_T_CHARS_FA struct of the face characteristics. provides .spType member for comparison with the sparse type.
     **/
    template< typename TL_T_INT_LID,
              typename TL_T_REAL,
              typename TL_T_INT_SP,
              typename TL_T_CHARS_FA >
    static void initFaces( TL_T_INT_LID               i_nFaDe,
                           TL_T_INT_SP                i_spType,
                           TL_T_CHARS_FA  const     * i_charsFa,
                           TL_T_INT_LID   const    (* i_faEl)[2],
                           TL_T_INT_LID   const    (* i_elFa)[ TL_N_EL_FA ],
                           unsigned short const    (* i_vIdElFaEl)[ TL_N_EL_FA ],
                           t_InternalBoundaryFace<
                             TL_T_REAL,
                             TL_T_INT_SP
                           >                        * o_intFa ) {
      // id of the sparse internal boundary faces
      TL_T_INT_LID l_spId = 0;

      // iterate over dense faces
      for( TL_T_INT_LID l_fa = 0; l_fa < i_nFaDe; l_fa++ ) {
        // check if this is a face of the internal boundary
        if( (i_charsFa[l_fa].spType & i_spType) != i_spType ) continue;

        // set the sparse type
        o_intFa[l_spId].spType = (TL_T_INT_SP) i_charsFa[l_fa].spType;

        // iterate over adjacent elements
        for( unsigned short l_sd = 0; l_sd < 2; l_sd++ ) {
          TL_T_INT_LID l_el = i_faEl[l_fa][l_sd];
          EDGE_CHECK( l_el != std::numeric_limits< TL_T_INT_LID >::max() );

          // find the local face and vertex id
          for( unsigned short l_fe = 0; l_fe < TL_N_EL_FA; l_fe++ ) {
            if( i_elFa[l_el][l_fe] == l_fa ) {
              o_intFa[l_spId].fIdFaEl[l_sd] = l_fe;

              // the left element holds the right elements vertex id
              if( l_sd == 0 ) o_intFa[l_spId].vIdFaElR = i_vIdElFaEl[l_el][l_fe];
              break;
            }
            // check that we found every thing
            EDGE_CHECK( l_fe != TL_N_EL_FA-1 );
          }
        }

        l_spId++;
      }
    }

    /**
     * Manipulate the LTS-types to match the internal boundary requirements:
     *   If internal boundary face is at an MPI-boundary, the adajcent elements store the
     *   plain DOFs in the tDOFs, from which the quad points in the internal boundary solver
     *   are assembled.
     *   Additionally, all other surface intergration, which used this data, has to perform
     *   an additional time integration.
     *
     * @param i_faLayout layout of the faces.
     * @param i_spType sparse type of the internal boundary.
     * @param i_faEl elements adjacent to the faces.
     * @param io_faChars face characteristics which will be updated accordingly.
     * @param io_elChars element characteristiscs which wil be updated accordingly.
     **/
    template< typename TL_T_INT_SP,
              typename TL_T_INT_LID,
              typename TL_T_CHARS_FA,
              typename TL_T_CHARS_EL >
    static void initMpi( t_enLayout   const  &i_faLayout,
                         TL_T_INT_SP          i_spType,
                         TL_T_INT_LID const (*i_faEl),
                         TL_T_CHARS_FA        io_faChars,
                         TL_T_CHARS_EL        io_elChars ) {
      // iterate over the time groups
      for( std::size_t l_tg = 0; l_tg < i_faLayout.timeGroups.size(); l_tg++ ) {
        // determine first send-face
        TL_T_INT_LID l_first  = i_faLayout.timeGroups[l_tg].inner.first;
                     l_first += i_faLayout.timeGroups[l_tg].inner.size;
        // determine number of send/receive faces
        TL_T_INT_LID l_size  = i_faLayout.timeGroups[l_tg].nEntsOwn;
                     l_size += i_faLayout.timeGroups[l_tg].nEntsNotOwn;
                     l_size -= i_faLayout.timeGroups[l_tg].inner.size;

        // iterate over send and receive faces
        for( TL_T_INT_LID l_fa = l_first; l_fa < l_first+l_size; l_fa++ ) {
          // check if the face is part of the internal boundary
          if( (io_faChars.spType & i_spType) == i_spType ) {
            // determine adjacent elements
            TL_T_INT_LID l_el[2] = i_faEl[l_fa];

            // check for extising adjacent elements
            EDGE_CHECK( l_el[0] < std::numeric_limits< TL_T_INT_LID >::max() );
            EDGE_CHECK( l_el[1] < std::numeric_limits< TL_T_INT_LID >::max() );

            // set the elements' LTS flags for plain DOFs
            for( unsigned short l_sd = 0; l_sd < 2; l_sd++ ) {
              io_elChars[l_el[l_sd]].spType |= C_LTS_EL[EL_DOFS];
            }
          }
        }
      }
    }
};

/**
 * Initializes 2D solvers used for internal boundaries.
 *
 * @paramt TL_T_EL two-dimensional element type.
 **/
template< t_entityType TL_T_EL >
class edge::elastic::solvers::InternalBoundarySolvers<TL_T_EL, 2> {
  //! #vertices of the elements
  static unsigned short const TL_N_EL_VE = C_ENT[TL_T_EL].N_VERTICES;

#ifndef __INTEL_COMPILER
  static_assert( C_ENT[TL_T_EL].N_DIM==2, "template type and dimensions don't match" );
#endif

  private:
    /**
     * Initializes the middle state solver for the 2D elastic wave equations in velocity-stress formulation.
     * The result of the solver S the combines the jump over all left-going waves.
     * S can be used directly to compute the middle state:
     *   q^l_m = q_l + S * (q_r - q_l)
     *
     * Remark: It is assumed that the normal points from left to right.
     *
     * @paramt TL_T_REAL precision of the derived middle state solver.
     * @param i_lamL lame parameter lambda of the left element.
     * @param i_lamR lame parameter lambda of the right element.
     * @param i_muL lame parameter mu of the left element.
     * @param i_muR lame parameter mu of the right element.
     * @param i_rhoL density rho of the left element.
     * @param i_rhoR density rho of the right element.
     * @param o_msJump will be set to solver computing the single jump to the middle state from the left.
     **/
    template <typename TL_T_REAL>
    static void solvMsJumpL( TL_T_REAL i_lamL, TL_T_REAL i_lamR,
                             TL_T_REAL i_muL,  TL_T_REAL i_muR,
                             TL_T_REAL i_rhoL, TL_T_REAL i_rhoR,
                             TL_T_REAL o_msJump[5][5] ) {
#include "../generated/MiddleStateJumpL2D.inc"
    }

    /**
     * Initializes the middle state solver for the 2D elastic wave equations in velocity-stress formulation.
     * The result of the solver S the combines the jump over all right-going waves.
     * S can be used directly to compute the middle state:
     *   q^l_m = q_r + S * (q_r - q_l)
     *
     * Remark: It is assumed that the normal points from left to right.
     *
     * @paramt TL_T_REAL precision of the derived middle state solver.
     * @param i_lamL lame parameter lambda of the left element.
     * @param i_lamR lame parameter lambda of the right element.
     * @param i_muL lame parameter mu of the left element.
     * @param i_muR lame parameter mu of the right element.
     * @param i_rhoL density rho of the left element.
     * @param i_rhoR density rho of the right element.
     * @param o_msJump will be set to solver computing the single jump to the middle state from the right.
     **/
    template <typename TL_T_REAL>
    static void solvMsJumpR( TL_T_REAL i_lamL, TL_T_REAL i_lamR,
                             TL_T_REAL i_muL,  TL_T_REAL i_muR,
                             TL_T_REAL i_rhoL, TL_T_REAL i_rhoR,
                             TL_T_REAL o_msJump[5][5] ) {
#include "../generated/MiddleStateJumpR2D.inc"
    }

    /**
     * Initializes the flux solver (using middle states) for the 2D elastic wave equations in velocity-stress formulation.
     *
     * Remark 1: It is assumed that the normal points from left to right.
     *
     * Remark 2: The normal is supposed to have unit length: (i_nx*i_nx+i_ny*i_ny==1) 
     *
     * @paramt TL_T_MESH precision of mesh coordinates.
     * @paramt TL_T_SOLVER precision of the solver.
     * @param i_nx x-component of face's normal.
     * @param i_ny y-component of face's normal.
     * @param i_lam lame parameter lambda.
     * @param i_mu lame parameter mu.
     * @param i_rho density rho.
     * @param o_msFlux will be set to solver for fluxes from computed middle states.
     **/
    template <typename TL_T_MESH, typename TL_T_SOLVER>
    static void solvMsFlux( TL_T_MESH i_nx, TL_T_MESH i_ny,
                            TL_T_SOLVER i_lam, TL_T_SOLVER i_mu, TL_T_SOLVER i_rho,
                            TL_T_SOLVER o_msFlux[5][5] ) {
#include "../generated/MiddleStateFlux2D.inc"
    }

  public:
    /**
     * Initializes the internal boundary solvers for the sparse faces of the given faces.
     *
     * @param i_nFa number of dense faces. 
     * @param i_spType sparse type.
     * @param i_faEl adjacency from faces to elements.
     * @param i_elVe adjacency from elements to vertices.
     * @param i_veChars vertex characteristics.
     * @param i_faChars face characteristics.
     * @param i_bgPars background parameters.
     * @param o_solverSp will be set so solver of the sparse faces. [0]: trafo to face-aligned coords, [1]: middle state, [2]: flux left, [3]: flux right.
     * @param i_bndCrds boundary coordinate system. [0]: vector pointing in the direction "left" to "right" side. If the nullptr is given, the outward-ponting normals of the mesh are used in the setup. These dependent on the global mesh indices and point from the lower to the upper index. Thus, the face-normal coordinate system in incosistent from face to face and might use either side of the internal boundary for the normal. In contrast if a vector is given, all normals are ensured to point into the normalDir-direction.
     *
     * @paramt TL_T_INT_LID integer type of local entities.
     * @paramt TL_T_INT_SP integer type for the sparse type.
     * @paramt TL_T_VE_CHARS struct of the vertex characteristics (defines .coords[2+]).
     * @paramt TL_T_FA_CHARS struct of the face characteristics (defines .outNormal, .area).
     * @paramt TL_T_BG_PARS struct of the background parameters (defines .lam, .mu and .rho).
     * @paramt TL_T_REAL_MESH floating point precision of mesh-related data.
     * @paramt TL_T_REAL_COMP floating point precision of computational data.
     **/
    template< typename TL_T_INT_LID,
              typename TL_T_INT_SP,
              typename TL_T_VE_CHARS,
              typename TL_T_FA_CHARS,
              typename TL_T_BG_PARS,
              typename TL_T_REAL_MESH,
              typename TL_T_REAL_COMP >
    static void init( TL_T_INT_LID           i_nFa,
                      TL_T_INT_SP            i_spType,
                      TL_T_INT_LID   const (*i_faEl)[2],
                      TL_T_INT_LID   const (*i_elVe)[TL_N_EL_VE],
                      TL_T_VE_CHARS  const  *i_veChars,
                      TL_T_FA_CHARS  const  *i_faChars,
                      TL_T_BG_PARS   const  *i_bgPars,
                      TL_T_REAL_COMP       (*o_solversSp)[4][5][5],
                      TL_T_REAL_MESH         i_bndCrds[2][2] = nullptr ) {
      TL_T_INT_LID l_spId = 0;

      // iterate over dense faces
      for( TL_T_INT_LID l_fa = 0; l_fa < i_nFa; l_fa++ ) {
        // only continue for rupture elements
        if( (i_faChars[l_fa].spType & i_spType) == i_spType ) {
          // get elements
          TL_T_INT_LID l_el[2];
          l_el[0] = i_faEl[l_fa][0];
          l_el[1] = i_faEl[l_fa][1];
          EDGE_CHECK( l_el[0] != l_el[1] );

          // set dummy values for boundary conditions
          l_el[0] = (l_el[0] != std::numeric_limits< TL_T_INT_LID >::max() ) ? l_el[0] : l_el[1];
          l_el[1] = (l_el[1] != std::numeric_limits< TL_T_INT_LID >::max() ) ? l_el[1] : l_el[0];

          // determine if the face points into the preferred direction
          bool l_prefDir = true;

          if( i_bndCrds != nullptr ) {
            TL_T_REAL_MESH l_sProd = linalg::Geom::sprod2( i_bndCrds[0], i_faChars[l_fa].outNormal );
            EDGE_CHECK( std::abs(l_sProd) > TOL.MESH );

            l_prefDir = (l_sProd > 0);
          }

          // direction of the face
          TL_T_REAL_MESH l_dir[2];
          l_dir[0] = i_faChars[l_fa].outNormal[0];
          l_dir[1] = i_faChars[l_fa].outNormal[1];

          // change direction if necessary
          if( l_prefDir == false ) {
            l_dir[0] *= -1;
            l_dir[1] *= -1;
          }

          edge::elastic::common::setupTrafoInv2d( l_dir[0],
                                                  l_dir[1],
                                                  o_solversSp[l_spId][0] );

          if( l_prefDir == true ) {
            solvMsJumpL( i_bgPars[l_el[0]].lam, i_bgPars[l_el[1]].lam,
                         i_bgPars[l_el[0]].mu,  i_bgPars[l_el[1]].mu,
                         i_bgPars[l_el[0]].rho, i_bgPars[l_el[1]].rho,
                         o_solversSp[l_spId][1] );
          } else {
            // the position of left and right changes if the face's normal is not
            // in the preferred direction.
            // remark:    logically (mesh perspective) we still use the left element
            //            and its respective material properties
            //         -> "left material parameters used as right"
            solvMsJumpR( i_bgPars[l_el[1]].lam, i_bgPars[l_el[0]].lam,
                         i_bgPars[l_el[1]].mu,  i_bgPars[l_el[0]].mu,
                         i_bgPars[l_el[1]].rho, i_bgPars[l_el[0]].rho,
                         o_solversSp[l_spId][1] );

            // change the sign of the solver, which is multuplied to the jump in quantities Rh-Lh:
            // Rh-Lh = -(R - L) if Rh=L and Lh=R
            for( unsigned short l_q1 = 0; l_q1 < 5; l_q1++ ) {
              for( unsigned short l_q2 = 0; l_q2 < 5; l_q2++ ) {
                o_solversSp[l_spId][1][l_q1][l_q2] *= -1;
              }
            }
          }

          solvMsFlux( l_dir[0],
                      l_dir[1],
                      i_bgPars[l_el[0]].lam,
                      i_bgPars[l_el[0]].mu,
                      i_bgPars[l_el[0]].rho,
                      o_solversSp[l_spId][2] );

          solvMsFlux( l_dir[0],
                      l_dir[1],
                      i_bgPars[l_el[1]].lam,
                      i_bgPars[l_el[1]].mu,
                      i_bgPars[l_el[1]].rho,
                      o_solversSp[l_spId][3] );

          // scale the solvers by scalars
          // TODO: This is redundant to the quad-free implementation and needs work..

          for( unsigned short l_sd = 0; l_sd < 2; l_sd++ ) {
            // vertex coordinates
            TL_T_REAL_MESH l_veCrds[2][TL_N_EL_VE];

            for( unsigned short l_ve = 0; l_ve < TL_N_EL_VE; l_ve++ ) {
              TL_T_INT_LID l_veId = i_elVe[ l_el[l_sd] ][ l_ve ];

              for( unsigned short l_di = 0; l_di < 2; l_di++ ) {
                l_veCrds[l_di][l_ve] = i_veChars[l_veId].coords[l_di];
              }
            }

            // get the jacobian
            TL_T_REAL_MESH l_jac[2][2];
            linalg::Mappings::evalJac( TL_T_EL, l_veCrds[0], l_jac[0] );

            // get determinant
            TL_T_REAL_MESH l_det = linalg::Matrix::det2x2( l_jac );

            // get scaling factor
            TL_T_REAL_MESH l_sca = i_faChars[l_fa].area / l_det;

            // adjust sign of flux if the direction is inversed
            if( l_prefDir == false ) l_sca *= -1;

            // scale the flux solver
            for( unsigned short l_q1 = 0; l_q1 < 5; l_q1++ ) {
              for( unsigned short l_q2 = 0; l_q2 < 5; l_q2++ ) {
                o_solversSp[l_spId][2+l_sd][l_q1][l_q2] *= l_sca;
              }
            }
          }
          l_spId++;
        }
      }
    }
};

/**
 * Initializes 3D solvers used for internal boundaries.
 *
 * @paramt TL_T_EL three-dimensional element type.
 **/
template< t_entityType TL_T_EL >
class edge::elastic::solvers::InternalBoundarySolvers<TL_T_EL, 3> {
  //! #vertices of the elements
  static unsigned short const TL_N_EL_VE = C_ENT[TL_T_EL].N_VERTICES;

#ifndef __INTEL_COMPILER
  static_assert( C_ENT[TL_T_EL].N_DIM==3, "template type and dimensions don't match" );
#endif

  private:
    /**
     * Initializes the middle state solver for the 3D elastic wave equations in velocity-stress formulation.
     * The result of the solver S the combines the jump over all left-going waves.
     * S can be used directly to compute the middle state:
     *   q^l_m = q_l + S * (q_r - q_l)
     *
     * Remark: It is assumed that the normal points from left to right.
     *
     * @paramt TL_T_REAL precision of the derived middle state solver.
     * @param i_lamL lame parameter lambda of the left element.
     * @param i_lamR lame parameter lambda of the right element.
     * @param i_muL lame parameter mu of the left element.
     * @param i_muR lame parameter mu of the right element.
     * @param i_rhoL density rho of the left element.
     * @param i_rhoR density rho of the right element.
     * @param o_msJump will be set to solver computing the single jump to the middle state from the left.
     **/
    template <typename TL_T_REAL>
    static void solvMsJumpL( TL_T_REAL i_lamL, TL_T_REAL i_lamR,
                             TL_T_REAL i_muL,  TL_T_REAL i_muR,
                             TL_T_REAL i_rhoL, TL_T_REAL i_rhoR,
                             TL_T_REAL o_msJump[9][9] ) {
#include "../generated/MiddleStateJumpL3D.inc"
    }

    /**
     * Initializes the middle state solver for the 3D elastic wave equations in velocity-stress formulation.
     * The result of the solver S the combines the jump over all right-going waves.
     * S can be used directly to compute the middle state:
     *   q^l_m = q_r + S * (q_r - q_l)
     *
     * Remark: It is assumed that the normal points from left to right.
     *
     * @paramt TL_T_REAL precision of the derived middle state solver.
     * @param i_lamL lame parameter lambda of the left element.
     * @param i_lamR lame parameter lambda of the right element.
     * @param i_muL lame parameter mu of the left element.
     * @param i_muR lame parameter mu of the right element.
     * @param i_rhoL density rho of the left element.
     * @param i_rhoR density rho of the right element.
     * @param o_msJump will be set to solver computing the single jump to the middle state from the left.
     **/
    template <typename TL_T_REAL>
    static void solvMsJumpR( TL_T_REAL i_lamL, TL_T_REAL i_lamR,
                             TL_T_REAL i_muL,  TL_T_REAL i_muR,
                             TL_T_REAL i_rhoL, TL_T_REAL i_rhoR,
                             TL_T_REAL o_msJump[9][9] ) {
#include "../generated/MiddleStateJumpR3D.inc"
    }


    /**
     * Initializes the flux solver (using middle states) for the 3D elastic wave equations in velocity-stress formulation.
     *
     * Remark 1: It is assumed that the normal points from left to right (see wave strength solver).
     *
     * Remark 2: The basis (normal + 2 tangents) is supposed to be orthonormal.
     *
     * @paramt TL_T_MESH precision of mesh coordinates.
     * @paramt TL_T_SOLVER precision of the solver.
     * @param i_nx x-component of face's normal.
     * @param i_ny y-component of face's normal.
     * @param i_nz z-component of face's normal.
     * @param i_sx x-component of face's first tangent.
     * @param i_sy y-component of face's first tangent.
     * @param i_sz z-component of face's first tangent.
     * @param i_tx x-component of face's second tangent.
     * @param i_ty y-component of face's second tangent.
     * @param i_tz z-component of face's second tangent.
     * @param i_lam lame parameter lambda.
     * @param i_mu lame parameter mu.
     * @param i_rho density rhot.
     * @param o_msFlux will be set to solver for fluxes from computed middle states.
     **/
    template <typename TL_T_MESH, typename TL_T_SOLVER>
    static void solvMsFlux( TL_T_MESH i_nx, TL_T_MESH i_ny, TL_T_MESH i_nz,
                            TL_T_MESH i_sx, TL_T_MESH i_sy, TL_T_MESH i_sz,
                            TL_T_MESH i_tx, TL_T_MESH i_ty, TL_T_MESH i_tz,
                            TL_T_SOLVER i_lam, TL_T_SOLVER i_mu, TL_T_SOLVER i_rho,
                            TL_T_SOLVER o_msFlux[9][9] ) {
#include "../generated/MiddleStateFlux3D.inc"
    }

  public:
    /**
     * Initializes the internal boundary solvers for the sparse faces of the given faces.
     *
     * @param i_nFa number of dense faces. 
     * @param i_spType sparse type.
     * @param i_faEl adjacency from faces to elements.
     * @param i_elVe adjacency from elements to vertices.
     * @param i_veChars vertex characteristics.
     * @param i_faChars face characteristics.
     * @param i_bgPars background parameters.
     * @param o_solverSp will be set so solver of the sparse faces. [0]: trafo to face-aligned coords, [1]: middle state, [2]: flux left, [3]: flux right.
     *
     * @paramt TL_T_INT_LID integer type of local entities.
     * @paramt TL_T_INT_SP integer type for the sparse type.
     * @paramt TL_T_VE_CHARS struct of the vertex characteristics (defines .coords[3]).
     * @paramt TL_T_FA_CHARS struct of the face characteristics (defines .outNormal, .tangent0, .tangent1 and .area).
     * @paramt TL_T_BG_PARS struct of the background parameters (defines .lam, .mu and .rho).
     * @paramt TL_T_REAL_MESH floating point precision of mesh-related data.
     * @paramt TL_T_REAL_COMP floating point precision of computational data.
     **/
    template< typename TL_T_INT_LID,
              typename TL_T_INT_SP,
              typename TL_T_VE_CHARS,
              typename TL_T_FA_CHARS,
              typename TL_T_BG_PARS,
              typename TL_T_REAL_MESH,
              typename TL_T_REAL_COMP >
    static void init( TL_T_INT_LID           i_nFa,
                      TL_T_INT_SP            i_spType,
                      TL_T_INT_LID   const (*i_faEl)[2],
                      TL_T_INT_LID   const (*i_elVe)[TL_N_EL_VE],
                      TL_T_VE_CHARS  const  *i_veChars,
                      TL_T_FA_CHARS  const  *i_faChars,
                      TL_T_BG_PARS   const  *i_bgPars,
                      TL_T_REAL_COMP       (*o_solversSp)[4][9][9],
                      TL_T_REAL_MESH         i_bndCrds[3][3] = nullptr ) {
      TL_T_INT_LID l_spId = 0;

      // check the input boundary coordinate system
      if( i_bndCrds != nullptr ) {
        TL_T_REAL_MESH l_norm;
        l_norm = linalg::Geom::norm3( i_bndCrds[0] );
        EDGE_CHECK( std::abs(l_norm-1.0) < TOL.MESH ) << l_norm;
        l_norm = linalg::Geom::norm3( i_bndCrds[1] );
        EDGE_CHECK( std::abs(l_norm-1.0) < TOL.MESH ) << l_norm;
        l_norm = linalg::Geom::norm3( i_bndCrds[2] );
        EDGE_CHECK( std::abs(l_norm-1.0) < TOL.MESH ) << l_norm;
      }

      // iterate over dense faces
      for( TL_T_INT_LID l_fa = 0; l_fa < i_nFa; l_fa++ ) {
        // only continue for rupture elements
        if( (i_faChars[l_fa].spType & i_spType) == i_spType ) {
          // init face-local coordinate system
          TL_T_REAL_MESH l_faCrds[3][3];
          for( unsigned short l_di = 0; l_di < 3; l_di++ ) {
            l_faCrds[0][l_di] = i_faChars[l_fa].outNormal[l_di];
            l_faCrds[1][l_di] = i_faChars[l_fa].tangent0[l_di];
            l_faCrds[2][l_di] = i_faChars[l_fa].tangent1[l_di];
          }

          // true if the given normal points into the preferred direction
          // otherwise we have to change the definition of left and right in the
          // internal boundary solvers since the given face normal points from the left
          // to the right element
          bool l_prefDir = true;

          // adjust the face coordinate system, if a boundary coordinate system is given
          if( i_bndCrds != nullptr ) {
            // check if the normal points in the same direction
            TL_T_REAL_MESH l_sProd = linalg::Geom::sprod3( i_bndCrds[0], l_faCrds[0] );
            EDGE_CHECK( std::abs(l_sProd) > TOL.MESH );

            l_prefDir = (l_sProd > 0);

            // adjust the face normal if this is not the preferred direction
            if( l_prefDir == false ) {
              for( unsigned short l_di = 0; l_di < 3; l_di++ ) l_faCrds[0][l_di] *= -1;
            }

            // determine the rotation matrix, which brings us from the normal of the boundary's
            // coordinate system to the (probably adjusted) face-normal
            TL_T_REAL_MESH l_rm[3][3];

            linalg::GeomT<3>::rotMat( i_bndCrds[0], l_faCrds[0], l_rm );

            // apply the rotation to the two remaining basis vectors to obtain the tangents
            linalg::Matrix::matMulB0( 3, 1, 3,
                                      l_rm[0], i_bndCrds[1], l_faCrds[1] );

            linalg::Matrix::matMulB0( 3, 1, 3,
                                      l_rm[0], i_bndCrds[2], l_faCrds[2] );

            // double check that we have a valid face-local coordinate system
            l_sProd = linalg::Geom::sprod3( l_faCrds[0], l_faCrds[1] );
            EDGE_CHECK( std::abs(l_sProd) < TOL.MESH );
            l_sProd = linalg::Geom::sprod3( l_faCrds[0], l_faCrds[2] );
            EDGE_CHECK( std::abs(l_sProd) < TOL.MESH );
            l_sProd = linalg::Geom::sprod3( l_faCrds[1], l_faCrds[2] );
            EDGE_CHECK( std::abs(l_sProd) < TOL.MESH );

            TL_T_REAL_MESH l_norm;
            l_norm = linalg::Geom::norm3( l_faCrds[0] );
            EDGE_CHECK( std::abs(l_norm-1.0) < TOL.MESH ) << l_norm;
            l_norm = linalg::Geom::norm3( l_faCrds[1] );
            EDGE_CHECK( std::abs(l_norm-1.0) < TOL.MESH ) << l_norm;
            l_norm = linalg::Geom::norm3( l_faCrds[2] );
            EDGE_CHECK( std::abs(l_norm-1.0) < TOL.MESH ) << l_norm;
          }

          // get elements
          TL_T_INT_LID l_el[2];
          l_el[0] = i_faEl[l_fa][0];
          l_el[1] = i_faEl[l_fa][1];
          EDGE_CHECK( l_el[0] != l_el[1] );

          // set dummy values for boundary conditions
          l_el[0] = (l_el[0] != std::numeric_limits< TL_T_INT_LID >::max() ) ? l_el[0] : l_el[1];
          l_el[1] = (l_el[1] != std::numeric_limits< TL_T_INT_LID >::max() ) ? l_el[1] : l_el[0];

          edge::elastic::common::setupTrafoInv3d( l_faCrds[0][0],
                                                  l_faCrds[0][1],
                                                  l_faCrds[0][2],
                                                  l_faCrds[1][0],
                                                  l_faCrds[1][1],
                                                  l_faCrds[1][2],
                                                  l_faCrds[2][0],
                                                  l_faCrds[2][1],
                                                  l_faCrds[2][2],
                                                  o_solversSp[l_spId][0] );

          if( l_prefDir == true ) {
            solvMsJumpL( i_bgPars[l_el[0]].lam, i_bgPars[l_el[1]].lam,
                         i_bgPars[l_el[0]].mu,  i_bgPars[l_el[1]].mu,
                         i_bgPars[l_el[0]].rho, i_bgPars[l_el[1]].rho,
                         o_solversSp[l_spId][1] );
          }
          else {
            // the position of left and right changes if the face's normal is not
            // in the preferred direction.
            // remark:    logically (mesh perspective) we still use the left element
            //            and its respective material properties
            //         -> "left material parameters used as right"
            solvMsJumpR( i_bgPars[l_el[1]].lam, i_bgPars[l_el[0]].lam,
                         i_bgPars[l_el[1]].mu,  i_bgPars[l_el[0]].mu,
                         i_bgPars[l_el[1]].rho, i_bgPars[l_el[0]].rho,
                         o_solversSp[l_spId][1] );

            // change the sign of the solver, which is multiplied to the jump in quantities Rh-Lh:
            // Rh-Lh = -(R - L) if Rh=L and Lh=R
            for( unsigned short l_q1 = 0; l_q1 < 9; l_q1++ ) {
              for( unsigned short l_q2 = 0; l_q2 < 9; l_q2++ ) {
                o_solversSp[l_spId][1][l_q1][l_q2] *= -1;
              }
            }
          }

          solvMsFlux( l_faCrds[0][0],
                      l_faCrds[0][1],
                      l_faCrds[0][2],
                      l_faCrds[1][0],
                      l_faCrds[1][1],
                      l_faCrds[1][2],
                      l_faCrds[2][0],
                      l_faCrds[2][1],
                      l_faCrds[2][2],
                      i_bgPars[l_el[0]].lam,
                      i_bgPars[l_el[0]].mu,
                      i_bgPars[l_el[0]].rho,
                      o_solversSp[l_spId][2] );

          solvMsFlux( l_faCrds[0][0],
                      l_faCrds[0][1],
                      l_faCrds[0][2],
                      l_faCrds[1][0],
                      l_faCrds[1][1],
                      l_faCrds[1][2],
                      l_faCrds[2][0],
                      l_faCrds[2][1],
                      l_faCrds[2][2],
                      i_bgPars[l_el[1]].lam,
                      i_bgPars[l_el[1]].mu,
                      i_bgPars[l_el[1]].rho,
                      o_solversSp[l_spId][3] );

          // scale the solvers by scalars
          // TODO: This is redundant to the quad-free implementation and needs work..
          for( unsigned short l_sd = 0; l_sd < 2; l_sd++ ) {
            // vertex coordinates
            TL_T_REAL_MESH l_veCrds[3][TL_N_EL_VE];

            for( unsigned short l_ve = 0; l_ve < TL_N_EL_VE; l_ve++ ) {
              TL_T_INT_LID l_veId = i_elVe[ l_el[l_sd] ][ l_ve ];

              for( unsigned short l_di = 0; l_di < 3; l_di++ ) {
                l_veCrds[l_di][l_ve] = i_veChars[l_veId].coords[l_di];
              }
            }

            // get the jacobian
            TL_T_REAL_MESH l_jac[3][3];
            linalg::Mappings::evalJac( TL_T_EL, l_veCrds[0], l_jac[0] );

            // get determinant
            TL_T_REAL_MESH l_det = linalg::Matrix::det3x3( l_jac );

            // get scaling factor
            TL_T_REAL_MESH l_sca = i_faChars[l_fa].area / l_det;

            // adjust sign of flux if the direction is inversed
            if( l_prefDir == false ) l_sca *= -1;

            // scale with 2 to account for 0.5 volume of triangular reference element
            if( TL_T_EL == TET4 ) l_sca *= 2;

            // scale the flux solver
            for( unsigned short l_q1 = 0; l_q1 < 9; l_q1++ ) {
              for( unsigned short l_q2 = 0; l_q2 < 9; l_q2++ ) {
                o_solversSp[l_spId][2+l_sd][l_q1][l_q2] *= l_sca;
              }
            }
          }
          l_spId++;
        }
      }
    }
};

#endif
