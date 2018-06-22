/**
 * @file This file is part of EDGE.
 *
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
 * Support for seismic internal boundary conditions through sub-cell limiting.
 **/

#ifndef EDGE_SEISMIC_INTERNAL_BOUNDARY_HPP
#define EDGE_SEISMIC_INTERNAL_BOUNDARY_HPP

#include "constants.hpp"
#include "io/logging.h"
#include "linalg/Matrix.h"
#include "linalg/Geom.hpp"
#include "common.hpp"

namespace edge {
  namespace elastic {
    namespace solvers {
      template< t_entityType   TL_T_EL,
                unsigned short TL_N_QTS,
                unsigned short TL_O_SP,
                unsigned short TL_N_CRS >
      class InternalBoundary;

      template< t_entityType TL_T_EL >
      class InternalBoundaryTypes;

#ifdef __INTEL_COMPILER
      template< t_entityType   TL_T_EL,
                unsigned short TL_O_SP,
                unsigned short TL_N_DIM=N_DIM >
#else
      template< t_entityType   TL_T_EL,
                unsigned short TL_O_SP,
                unsigned short TL_N_DIM=C_ENT[TL_T_EL].N_DIM >
#endif
      class InternalBoundarySolvers;

      template< t_entityType   TL_T_EL,
                unsigned short TL_O_SP >
      class InternalBoundarySolvers< TL_T_EL, TL_O_SP, 2 >;

      template< t_entityType   TL_T_EL,
                unsigned short TL_O_SP >
      class InternalBoundarySolvers< TL_T_EL, TL_O_SP, 3 >;
    }
  }
}

/**
 * Internal boundary conditions through sub-cell limiting.
 *
 * @paramt TL_T_EL element type.
 * @paramt TL_N_QTS number of quantities.
 * @paramt TL_O_SP order of the used DG method in space.
 * @paramt TL_N_CRS number of concurrent forward runs (fused simulations).
 **/
template< t_entityType   TL_T_EL,
          unsigned short TL_N_QTS,
          unsigned short TL_O_SP,
          unsigned short TL_N_CRS >
class edge::elastic::solvers::InternalBoundary {
  private:
    //! number of sub-faces
    static unsigned short const TL_N_SFS = CE_N_SUB_FACES( TL_T_EL, TL_O_SP );

    //! id of the matrix kernel group
    static unsigned short const MM_GR = static_cast< unsigned short >( t_mm::SUB_CELL );

  public:
    /**
     * Dummy solver perturbating nothing.
     **/
    class DummySolv {
      public:
        /**
         * Dummy perturbations.
         *
         * @param i_ms middle state which is not perturbed.
         * @param o_msL will be set to input middle state.
         * @param o_msR will be set to input middle state.
         * @param o_per will be set to false.
         *
         * @paramt TL_T_REAL floating point precision.
         * @paramt TL_T_FA_DATA type of face data (unused).
         **/
        template< typename TL_T_REAL, typename TL_T_FA_DATA >
        static void inline perturb( unsigned short,
                                    TL_T_REAL,
                                    TL_T_REAL            i_ms[TL_N_QTS][TL_N_CRS],
                                    TL_T_FA_DATA const *,
                                    TL_T_REAL            o_msL[TL_N_QTS][TL_N_CRS],
                                    TL_T_REAL            o_msR[TL_N_QTS][TL_N_CRS],
                                    TL_T_REAL            o_per[TL_N_CRS] ) {
          for( unsigned short l_qt = 0; l_qt < TL_N_QTS; l_qt++ ) {
            for( unsigned short l_ru = 0; l_ru < TL_N_CRS; l_ru++ ) {
              o_msL[l_qt][l_ru] = i_ms[l_qt][l_ru];
              o_msR[l_qt][l_ru] = i_ms[l_qt][l_ru];
              o_per[l_ru] = false;
            }
          }
        }
    };

    /**
     * Evaluates the internal boundary condition for all sub-faces at a DG-face of the given element type.
     *
     * Tria3-example:
     *
     *
     * 1) Illustatration of an element in physical coordinates
     *
     *   ***xxx***: DG-face sub-divided into three sub-faces
     *
     *   vertices to faces: f0: 0-1, f1: 1-2, f2: 2-0
     *
     *              *
     *             *2* 1  *
     *            *   *         *
     *           *     x           0  *
     *          *0  L   x     R     *
     *             *     x        *
     *                *   *     *
     *                   *1*2 *
     *                      * 
     *   In this example the left element's face f1 matches the right elements face f2.
     *   Vertex 2 of the right element lies on the first vertex of the left elements face f1.
     *
     * 2) Illustration of the sub-face-local Riemann problem.
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
     *    and the repspective sub-cell net-updates are obtained for the two sides.
     *
     * @param i_tm1 transformation matrix from physical coordinates for face-aligned coordinates.
     * @param i_solMsJumpL solver for the single jump from the left element's quantities to the middle state.
     * @param i_solMsFluxL flux solver using (probably perturbed) middle states for the left element.
     * @param i_solMsFluxR flux solver using (probably perturbed) middle states for the right element.
     * @param i_dofsL DOFs of the left element's face-adjacent sub-cells.
     * @param i_dofsR DOFs of the right element's face-adjacent sub-cells.
     * @param o_netUpsL will be set to left-going sub-cell net-updates.
     * @param o_netUpsL will be set to right-going sub-cell net-updates.
     * @param o_per will be set to true if middle states were perturbed, false otherwise.
     * @param i_dt time step, used for scaling the net-updates; might be used internally to compute slip from the slip rate, for example.
     * @param io_faData data used in the pertubation of the middle states.
     *
     * @paramt TL_T_REAL floating point type.
     * @paramt TL_T_MS_SOLV middle state "solver", offers member function .perturb.
     * @paramt TL_T_FA_DATA data passed to middle state solver.
     **/
    template< typename TL_T_REAL,
              typename TL_T_MM,
              typename TL_T_MS_SOLV = DummySolv,
              typename TL_T_FA_DATA = void >
    static void netUpdates( TL_T_REAL      const  i_tm1[TL_N_QTS][TL_N_QTS],
                            TL_T_REAL      const  i_solMsJumpL[TL_N_QTS][TL_N_QTS],
                            TL_T_REAL      const  i_solMsFluxL[TL_N_QTS][TL_N_QTS],
                            TL_T_REAL      const  i_solMsFluxR[TL_N_QTS][TL_N_QTS],
                            TL_T_REAL      const  i_dofsL[TL_N_QTS][TL_N_SFS][TL_N_CRS],
                            TL_T_REAL      const  i_dofsR[TL_N_QTS][TL_N_SFS][TL_N_CRS],
                            TL_T_REAL             o_netUpsL[TL_N_QTS][TL_N_SFS][TL_N_CRS],
                            TL_T_REAL             o_netUpsR[TL_N_QTS][TL_N_SFS][TL_N_CRS],
                            bool                  o_per[TL_N_SFS][TL_N_CRS],
                            TL_T_MM        const &i_mm,
                            TL_T_REAL             i_dt = 0,
                            TL_T_FA_DATA         *io_faData = nullptr
) {
      // temporary storage for the middle states
      TL_T_REAL l_msTmp[2][TL_N_QTS][TL_N_SFS][TL_N_CRS];
      for( unsigned short l_qt = 0; l_qt < TL_N_QTS; l_qt++ ) {
        for( unsigned short l_sf = 0; l_sf < TL_N_SFS; l_sf++ ) {
          for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {
            l_msTmp[0][l_qt][l_sf][l_cr] = 0;
            l_msTmp[1][l_qt][l_sf][l_cr] = 0;
          }
        }
      }

      // rotate the DOFs from physical coordinates to face-aligned coords
      // remark: the back-rotation to physical coordinates is part of the the flux solver
      TL_T_REAL l_dofs[2][TL_N_QTS][TL_N_SFS][TL_N_CRS];
#if defined(PP_T_KERNELS_XSMM_DENSE_SINGLE)
      i_mm.m_kernels[MM_GR][4]( i_dofsL[0][0], i_tm1[0], l_dofs[0][0][0] );
      i_mm.m_kernels[MM_GR][4]( i_dofsR[0][0], i_tm1[0], l_dofs[1][0][0] );
#else
      i_mm.m_kernels[MM_GR][4]( i_tm1[0], i_dofsL[0][0], l_dofs[0][0][0] );
      i_mm.m_kernels[MM_GR][4]( i_tm1[0], i_dofsR[0][0], l_dofs[1][0][0] );
#endif

      // iterate over sub-faces
      for( unsigned short l_sf = 0; l_sf < TL_N_SFS; l_sf++ ) {
        // temporary values at the sub-faces
        TL_T_REAL l_qVal[2][TL_N_QTS][TL_N_CRS];
        // jump in quantities
        TL_T_REAL l_qJump[TL_N_QTS][TL_N_CRS];

        // copy left and right elements' DOFs at sub-face (SoA -> AoS), compute jump
        for( unsigned short l_qt = 0; l_qt < TL_N_QTS; l_qt++ ) {
          // compute the jump in quantities
          for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {
            l_qVal[0][l_qt][l_cr] = l_dofs[0][l_qt][l_sf][l_cr];
            l_qVal[1][l_qt][l_cr] = l_dofs[1][l_qt][l_sf][l_cr];

            l_qJump[l_qt][l_cr] = l_qVal[1][l_qt][l_cr] - l_qVal[0][l_qt][l_cr];
          }
        }

        // jump over waves with negative speeds from the left to get the left-side middle state
        linalg::Matrix::matMulFusedBC( TL_N_CRS,
                                       TL_N_QTS, 1, TL_N_QTS,
                                       TL_N_QTS, 1, 1,
                                       static_cast<TL_T_REAL>(1.0),
                                       i_solMsJumpL[0],
                                       l_qJump[0],
                                       l_qVal[0][0] );

        // perturb if necessary
        TL_T_REAL l_ms[2][TL_N_QTS][TL_N_CRS];
        TL_T_MS_SOLV::perturb( l_sf,
                               i_dt,
                               l_qVal[0],
                               io_faData,
                               l_ms[0],
                               l_ms[1],
                               o_per[l_sf] );

        // scale and save middle states (AoS -> SoA)
        for( unsigned short l_qt = 0; l_qt < TL_N_QTS; l_qt++ ) {
          for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {
            // left-going fluxes are subtracted
            l_msTmp[0][l_qt][l_sf][l_cr] = -l_ms[0][l_qt][l_cr] * i_dt;
            // right-going fluxes are added
            l_msTmp[1][l_qt][l_sf][l_cr] =  l_ms[1][l_qt][l_cr] * i_dt;
          }
        }
      }

      // compute fluxes and rotate DOFs back to physical coordinate system
#if defined(PP_T_KERNELS_XSMM_DENSE_SINGLE)
      i_mm.m_kernels[MM_GR][4]( l_msTmp[0][0][0], i_solMsFluxL[0], o_netUpsL[0][0] );
      i_mm.m_kernels[MM_GR][4]( l_msTmp[1][0][0], i_solMsFluxR[0], o_netUpsR[0][0] );
#else
      i_mm.m_kernels[MM_GR][4]( i_solMsFluxL[0], l_msTmp[0][0][0], o_netUpsL[0][0] );
      i_mm.m_kernels[MM_GR][4]( i_solMsFluxR[0], l_msTmp[1][0][0], o_netUpsR[0][0] );
#endif
    }
};

template< t_entityType TL_T_EL >
class edge::elastic::solvers::InternalBoundaryTypes {
  private:
    //! number of vertices per element
    static unsigned short const TL_N_EL_VE = C_ENT[T_SDISC.ELEMENT].N_VERTICES;
    //! number of faces per element
    static unsigned short const TL_N_EL_FA = C_ENT[T_SDISC.ELEMENT].N_FACES;


  public:
    /**
     * Initializes the data of the internal boundary faces.
     *
     * @param i_nFa number of faces.
     * @param i_nDe number of elements.
     * @param i_spType sparse type of faces at the internal boundary.
     * @param i_charsFa characteristics of the faces.
     * @param i_charsEl characteristics of the elements.
     * @param i_faEl elements adjacent to the faces.
     * @param i_elFa faces adjacent to the elements.
     * @param i_vIdElFaEl vertex ids of the shared face with respect to the dense element's adjacent dense elements.
     *
     * @paramt TL_T_LID integer type of local ids.
     * @paramt TL_T_SP the sparse type of the internal boundary.
     * @paramt TL_T_CHARS_FA struct of the face characteristics. provides .spType.
     * @paramt TL_T_CHARS_EL struct of element chars, provides .spType.
     **/
    template< typename TL_T_LID,
              typename TL_T_SP,
              typename TL_T_CHARS_FA,
              typename TL_T_CHARS_EL >
    static void initFaces( TL_T_LID                   i_nFa,
                           TL_T_LID                   i_nEl,
                           TL_T_SP                    i_spType,
                           TL_T_CHARS_FA  const     * i_charsFa,
                           TL_T_CHARS_EL  const     * i_charsEl,
                           TL_T_LID       const    (* i_faEl)[2],
                           TL_T_LID       const    (* i_elFa)[ TL_N_EL_FA ],
                           unsigned short const (* i_vIdElFaEl)[ TL_N_EL_FA ],
                           edge::sc::ibnd::t_bfChars<
                             TL_T_SP
                           >                        * o_bfChars ) {
      // id of the sparse internal boundary faces
      TL_T_LID l_bf = 0;

      // iterate over faces
      for( TL_T_LID l_fa = 0; l_fa < i_nFa; l_fa++ ) {
        // check if this is a face of the internal boundary
        if( (i_charsFa[l_fa].spType & i_spType) != i_spType ) continue;

        // set the sparse type
        o_bfChars[l_bf].spType = (TL_T_SP) i_charsFa[l_fa].spType;

        o_bfChars[l_bf].vIdFaElR = std::numeric_limits< unsigned short >::max();

        // iterate over adjacent elements
        for( unsigned short l_sd = 0; l_sd < 2; l_sd++ ) {
          o_bfChars[l_bf].fIdBfEl[l_sd] = std::numeric_limits< unsigned short >::max();

          TL_T_LID l_el = i_faEl[l_fa][l_sd];
          EDGE_CHECK( l_el != std::numeric_limits< TL_T_LID >::max() );

          // find the local face and vertex id
          for( unsigned short l_fe = 0; l_fe < TL_N_EL_FA; l_fe++ ) {
            if( i_elFa[l_el][l_fe] == l_fa ) {
              o_bfChars[l_bf].fIdBfEl[l_sd] = l_fe;

              EDGE_CHECK(    o_bfChars[l_bf].vIdFaElR == std::numeric_limits< unsigned short >::max()
                          || i_vIdElFaEl[l_el][l_fe]  == std::numeric_limits< unsigned short >::max()
                          || o_bfChars[l_bf].vIdFaElR == i_vIdElFaEl[l_el][l_fe] );

              // compute minimum, as one of the elements might be recv
              o_bfChars[l_bf].vIdFaElR = std::min( o_bfChars[l_bf].vIdFaElR,
                                                   i_vIdElFaEl[l_el][l_fe] );
              break;
            }
            // check that we found every thing
            EDGE_CHECK_NE( l_fe, TL_N_EL_FA-1 );
          }
        }
        for( unsigned short l_sd = 0; l_sd < 2; l_sd++ )
          EDGE_CHECK_LT( o_bfChars[l_bf].fIdBfEl[l_sd], TL_N_EL_FA );
        EDGE_CHECK_LT( o_bfChars[l_bf].vIdFaElR , TL_N_EL_VE );

        l_bf++;
      }
    }
};

/**
 * Initializes 2D solvers used for internal boundaries.
 *
 * @paramt TL_T_EL two-dimensional element type.
 * @paramt TL_O_SP order of the used DG method in space.
 **/
template< t_entityType   TL_T_EL,
          unsigned short TL_O_SP >
class edge::elastic::solvers::InternalBoundarySolvers< TL_T_EL, TL_O_SP, 2 > {
  //! #vertices of the elements
  static unsigned short const TL_N_EL_VE = C_ENT[TL_T_EL].N_VERTICES;

  //! number of sub-faces
  static unsigned short const TL_N_SFS = CE_N_SUB_FACES( TL_T_EL, TL_O_SP );

  //! number of sub-cells
  static unsigned short const TL_N_SCS = CE_N_SUB_CELLS( TL_T_EL, TL_O_SP );

#if !defined(__INTEL_COMPILER) && !defined(__COVERITY__)
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
     * @param i_faChars face characteristics.
     * @param i_elChars element characteristics.
     * @param i_bgPars background parameters.
     * @param o_solverSp will be set so solver of the sparse faces. [0]: trafo to face-aligned coords, [1]: middle state, [2]: flux left, [3]: flux right.
     * @param i_bndCrds boundary coordinate system. [0]: vector pointing in the direction "left" to "right" side. If the nullptr is given, the outward-ponting normals of the mesh are used in the setup. These dependent on the global mesh indices and point from the lower to the upper index. Thus, the face-normal coordinate system in incosistent from face to face and might use either side of the internal boundary for the normal. In contrast if a vector is given, all normals are ensured to point into the normalDir-direction.
     *
     * @paramt TL_T_INT_LID integer type of local entities.
     * @paramt TL_T_INT_SP integer type for the sparse type.
     * @paramt TL_T_FA_CHARS struct of the face characteristics (defines .outNormal, .area).
     * @paramt TL_T_EL_CHARS struct of the element characteristics (defines .volume).
     * @paramt TL_T_BG_PARS struct of the background parameters (defines .lam, .mu and .rho).
     * @paramt TL_T_REAL_MESH floating point precision of mesh-related data.
     * @paramt TL_T_REAL_COMP floating point precision of computational data.
     **/
    template< typename TL_T_INT_LID,
              typename TL_T_INT_SP,
              typename TL_T_FA_CHARS,
              typename TL_T_EL_CHARS,
              typename TL_T_BG_PARS,
              typename TL_T_REAL_MESH,
              typename TL_T_REAL_COMP >
    static void init( TL_T_INT_LID           i_nFa,
                      TL_T_INT_SP            i_spType,
                      TL_T_INT_LID   const (*i_faEl)[2],
                      TL_T_INT_LID   const (*i_elVe)[TL_N_EL_VE],
                      TL_T_FA_CHARS  const  *i_faChars,
                      TL_T_EL_CHARS  const  *i_elChars,
                      TL_T_BG_PARS   const  *i_bgPars,
                      TL_T_REAL_COMP       (*o_solversSp)[4][5][5],
                      TL_T_REAL_MESH         i_bndCrds[2][2] = nullptr ) {
      TL_T_INT_LID l_spId = 0;

      // iterate over dense faces
      for( TL_T_INT_LID l_fa = 0; l_fa < i_nFa; l_fa++ ) {
        // only continue for elements of the internal boundary elements
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

            // change the sign of the solver, which is multiplied to the jump in quantities Rh-Lh:
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
          for( unsigned short l_sd = 0; l_sd < 2; l_sd++ ) {
            // get scaling factor
            TL_T_REAL_MESH l_sca  = i_faChars[ l_fa       ].area   / TL_N_SFS;
                           l_sca /= i_elChars[ l_el[l_sd] ].volume / TL_N_SCS;

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
 * @paramt TL_O_SP order of the used DG method in space.
 **/
template< t_entityType   TL_T_EL,
          unsigned short TL_O_SP >
class edge::elastic::solvers::InternalBoundarySolvers< TL_T_EL, TL_O_SP, 3 > {
  //! #vertices of the elements
  static unsigned short const TL_N_EL_VE = C_ENT[TL_T_EL].N_VERTICES;

  //! number of sub-faces
  static unsigned short const TL_N_SFS = CE_N_SUB_FACES( TL_T_EL, TL_O_SP );

  //! number of sub-cells
  static unsigned short const TL_N_SCS = CE_N_SUB_CELLS( TL_T_EL, TL_O_SP );

#if !defined(__INTEL_COMPILER) && !defined(__COVERITY__)
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
     * @param i_faChars face characteristics.
     * @param i_elChars element characteristics.
     * @param i_bgPars background parameters.
     * @param o_solverSp will be set so solver of the sparse faces. [0]: trafo to face-aligned coords, [1]: middle state, [2]: flux left, [3]: flux right.
     *
     * @paramt TL_T_INT_LID integer type of local entities.
     * @paramt TL_T_INT_SP integer type for the sparse type.
     * @paramt TL_T_FA_CHARS struct of the face characteristics (defines .outNormal, .tangent0, .tangent1 and .area).
     * @paramt TL_T_EL_CHARS struct of the element characteristics (defines .volume).
     * @paramt TL_T_BG_PARS struct of the background parameters (defines .lam, .mu and .rho).
     * @paramt TL_T_REAL_MESH floating point precision of mesh-related data.
     * @paramt TL_T_REAL_COMP floating point precision of computational data.
     **/
    template< typename TL_T_INT_LID,
              typename TL_T_INT_SP,
              typename TL_T_FA_CHARS,
              typename TL_T_EL_CHARS,
              typename TL_T_BG_PARS,
              typename TL_T_REAL_MESH,
              typename TL_T_REAL_COMP >
    static void init( TL_T_INT_LID           i_nFa,
                      TL_T_INT_SP            i_spType,
                      TL_T_INT_LID   const (*i_faEl)[2],
                      TL_T_INT_LID   const (*i_elVe)[TL_N_EL_VE],
                      TL_T_FA_CHARS  const  *i_faChars,
                      TL_T_EL_CHARS  const  *i_elChars,
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
          for( unsigned short l_sd = 0; l_sd < 2; l_sd++ ) {
            // get scaling factor
            TL_T_REAL_MESH l_sca  = i_faChars[ l_fa       ].area   / TL_N_SFS;
                           l_sca /= i_elChars[ l_el[l_sd] ].volume / TL_N_SCS;

            // adjust sign of flux if the direction is inversed
            if( l_prefDir == false ) l_sca *= -1;

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
