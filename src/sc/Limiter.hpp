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
 * Subcell limiter.
 **/
#ifndef EDGE_SC_LIMITER_HPP
#define EDGE_SC_LIMITER_HPP

#include "constants.hpp"
#include "io/logging.h"
#include "linalg/Matrix.h"
#include "Detections.hpp"
#include "ibnd/SuperCell.hpp"

namespace edge {
  namespace sc {
    template< t_entityType   TL_T_EL,
              unsigned short TL_O_SP,
              unsigned short TL_N_QTS,
              unsigned short TL_N_CRS,
              typename       TL_T_BND_SC >
    class Limiter;
  }
}

/**
 * Subcell limiter.
 *
 * @paramt TL_T_EL element type.
 * @paramt TL_O_SP spatial order of the DG-solver.
 * @paramt TL_N_QTS number of quantities.
 * @paramt TL_N_CRS number of fused forward runs.
 * @paramt TL_T_BND_SC type of the boundary conditions.
 **/
template< t_entityType   TL_T_EL,
          unsigned short TL_O_SP,
          unsigned short TL_N_QTS,
          unsigned short TL_N_CRS,
          typename       TL_T_BND_SC >
class edge::sc::Limiter {
  private:
    //! number of dimensions
    static unsigned short const TL_N_DIS = C_ENT[TL_T_EL].N_DIM;

    //! number of vertices
    static unsigned short const TL_N_VES = C_ENT[TL_T_EL].N_VERTICES;

    //! number of faces
    static unsigned short const TL_N_FAS = C_ENT[TL_T_EL].N_FACES;

    //! number of DG face modes
    static unsigned short const TL_N_MDS_FA = CE_N_ELEMENT_MODES( C_ENT[TL_T_EL].TYPE_FACES, TL_O_SP );

    //! number of DG element modes
    static unsigned short const TL_N_MDS_EL = CE_N_ELEMENT_MODES( TL_T_EL, TL_O_SP );

    //! number of DG flux matrices
    static unsigned short const TL_N_FLUXN = CE_N_FLUXN_MATRICES( TL_T_EL );

    //! number of sub-edges per element edge
    static unsigned short const TL_N_SES = CE_N_SUB_EDGES( TL_N_DIS, TL_O_SP );

    //! number of sub-faces per element face
    static unsigned short const TL_N_SFS = CE_N_SUB_FACES( TL_T_EL, TL_O_SP );

    //! number of subcells
    static unsigned short const TL_N_SCS = CE_N_SUB_CELLS( TL_T_EL, TL_O_SP );

    /**
     * Rolls back the DG surface integration for an DG element's face and replaces it with the sub-cell integral.
     * The face integral is replaced only if the candidate solution of the adjacent element is not admissible.
     *   Remark: This function assumes a limited plus element and an adjacent (face as bridge), limited element.
     *
     * @param i_dt time step.
     * @param i_lp id of the limited plus element.
     * @param i_fa respective face of the limited plus element.
     * @param i_admAdP admissibility of the adjacent element's DG solution for the previous solution.
     * @param i_admAdC admissibility of the adjacent element's DG solution for the candidate solution.
     * @param i_fIntL local flux matrix.
     * @param i_fIntN neighboring flux matrix.
     * @param i_fIntT transposed flux matrix.
     * @param i_scatterSf scatter operator for the element's sub-cells adjacent to the face in the reference element.
     * @param i_scatterSfAd scatter operator for the adjacent element's sub-cells adjacent to the face w.r.t. the reference element.
     * @param i_sfInt sfInt operator, which computes the DG face integral from the sub-cell fluxes.
     * @param i_fluxSolver flux solver for the element's local contribution to the face integral.
     * @param i_fluxSolverAd flux solver for the adjacent element's contribution to the face integral.
     * @param i_tDofsDgP time integrated DOFs of the DG element's previous solution.
     * @param i_tDofsDgAdP time integrated DOFs of the adjacent element's previous solution
     * @param i_tDofsScAdP previous sub-cell limited solution of the adjacent element for sub-cells adjacent to the element's face (if available).
     * @param i_dofsDgP DOFs of the element's previous DG solution.
     * @param io_dofsDg current candidate DG solution of the element, for which the DG surface integral is replaced with the sub-cell surface integral.
     * @param i_mm matrix multiplication kernels.
     * @param i_solvSc sub-cell solver.
     *
     * @paramt TL_T_LID ingral type of local ids.
     * @paramt TL_T_REAL floating point type.
     * @paramt TL_T_MM type of the matrix-matrix multiplication kernels.
     * @paramt TL_T_SOLV_SC type of the sub-cell solver.
     **/
   template< typename TL_T_LID,
             typename TL_T_REAL,
             typename TL_T_MM,
             typename TL_T_SOLV_SC >
   static void surfIntRb( TL_T_REAL           i_dt,
                          TL_T_LID            i_lp,
                          unsigned short      i_fa,
                          bool         const  i_admAdC[TL_N_CRS],
                          unsigned short const i_scDgAd[TL_N_SFS],
                          TL_T_REAL    const  i_fIntL[TL_N_MDS_EL][TL_N_MDS_FA],
                          TL_T_REAL    const  i_fIntN[TL_N_MDS_EL][TL_N_MDS_FA],
                          TL_T_REAL    const  i_fIntT[TL_N_MDS_FA][TL_N_MDS_EL],
                          TL_T_REAL    const  i_sfInt[TL_N_SFS][TL_N_MDS_EL],
                          TL_T_REAL    const  i_fluxSolver[TL_N_QTS][TL_N_QTS],
                          TL_T_REAL    const  i_fluxSolverAd[TL_N_QTS][TL_N_QTS],
                          TL_T_REAL    const  i_tDofsDgP[  TL_N_QTS][TL_N_MDS_EL][TL_N_CRS],
                          TL_T_REAL    const  i_tDofsDgAdP[TL_N_QTS][TL_N_MDS_EL][TL_N_CRS],
                          TL_T_REAL    const  i_tDofsScP[TL_N_QTS][TL_N_SFS][TL_N_CRS],
                          TL_T_REAL    const  i_tDofsScAdP[TL_N_QTS][TL_N_SFS][TL_N_CRS],
                          TL_T_REAL           io_dofsDg[   TL_N_QTS][TL_N_MDS_EL][TL_N_CRS],
                          TL_T_MM      const &i_mm,
                          TL_T_SOLV_SC const &i_solvSc ) {
      // determine if a rollback is required: a single fused sim of the adjacent limited element is not admissible
      bool l_rb = false;
      for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {
        l_rb = l_rb || (i_admAdC[l_cr] == false );
      }

      // only continue if a rollback is required
      if( l_rb ) {
        // sub-cell solution on both sides of the DG element's face
        TL_T_REAL l_scFa[2][TL_N_QTS][TL_N_SFS][TL_N_CRS];
        for( unsigned short l_qt = 0; l_qt < TL_N_QTS; l_qt++ ) {
          for( unsigned short l_sf = 0; l_sf < TL_N_SFS; l_sf++ ) {
            unsigned short l_sfRe = i_scDgAd[l_sf];

            for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {
              l_scFa[0][l_qt][l_sf][l_cr] = i_tDofsScP[  l_qt][l_sfRe][l_cr];
              l_scFa[1][l_qt][l_sf][l_cr] = i_tDofsScAdP[l_qt][l_sf  ][l_cr];
            }
          }
        }


        // compute sub-cell FV fluxes between the two DG elements
        TL_T_REAL l_fluxes[TL_N_QTS][TL_N_SFS][TL_N_CRS];
        i_solvSc.nuFaSf( i_dt,
                         i_lp,
                         i_fa,
                         l_scFa,
                         l_fluxes );

        // scale fluxes with number of sub-faces, as this is part of the surface integral
        for( unsigned short l_qt = 0; l_qt < TL_N_QTS; l_qt++ ) {
          for( unsigned short l_sf = 0; l_sf < TL_N_SFS; l_sf++ ) {
            for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {
              l_fluxes[l_qt][l_sf][l_cr] *= TL_N_SFS;
              l_fluxes[l_qt][l_sf][l_cr] /= TL_N_SCS;

              // scale fluxes for trias and tets (DG is scaled by det. of Jacobian, not the volume)
              if( TL_T_EL == TRIA3 ) l_fluxes[l_qt][l_sf][l_cr] *= TL_T_REAL(0.5);
              if( TL_T_EL == TET4  ) l_fluxes[l_qt][l_sf][l_cr] /= TL_T_REAL(3.0);
            }
          }
        }

        // face integral of the sub-cell fluxes
        TL_T_REAL l_sIntSc[TL_N_QTS][TL_N_MDS_EL][TL_N_CRS];

        // integrate fluxes over DG element's boundary and project back to DG space
        edge::sc::Kernels<
          TL_T_EL,
          TL_O_SP,
          TL_N_QTS,
          TL_N_CRS >::sfIntVanilla( l_fluxes,
                                    i_sfInt,
                                    l_sIntSc );

        // face integral of the DG fluxes
        TL_T_REAL l_sIntDg[TL_N_QTS][TL_N_MDS_EL][TL_N_CRS];
        for( unsigned short l_qt = 0; l_qt < TL_N_QTS; l_qt++ )
          for( unsigned short l_md = 0; l_md < TL_N_MDS_EL; l_md++ )
            for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ )
              l_sIntDg[l_qt][l_md][l_cr] = 0;

        TL_T_REAL l_scratch[2][TL_N_QTS][TL_N_MDS_FA][TL_N_CRS];

        // local
        elastic::solvers::SurfInt< TL_T_EL,
                                   TL_N_QTS,
                                   TL_O_SP,
                                   TL_O_SP,
                                   TL_N_CRS >::neigh( i_fIntL,
                                                      i_fIntT,
                                                      i_fluxSolver,
                                                      i_tDofsDgP,
                                                      i_mm,
                                                      l_sIntDg,
                                                      l_scratch,
                                                      i_tDofsDgAdP );

        // adjacent
        elastic::solvers::SurfInt< TL_T_EL,
                                   TL_N_QTS,
                                   TL_O_SP,
                                   TL_O_SP,
                                   TL_N_CRS >::neigh( i_fIntN,
                                                      i_fIntT,
                                                      i_fluxSolverAd,
                                                      i_tDofsDgAdP,
                                                      i_mm,
                                                      l_sIntDg,
                                                      l_scratch,
                                                      io_dofsDg );

        // update the DOFs accordingly
        for( unsigned short l_qt = 0; l_qt < TL_N_QTS; l_qt++ ) {
          for( unsigned short l_md = 0; l_md < TL_N_MDS_EL; l_md++ ) {
            for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {
              if( i_admAdC[l_cr] == false ) {
                io_dofsDg[l_qt][l_md][l_cr] -= l_sIntDg[l_qt][l_md][l_cr];
                io_dofsDg[l_qt][l_md][l_cr] += l_sIntSc[l_qt][l_md][l_cr];
              }
            }
          }
        }

      }
   }

    /**
     * Computes the limited solution of the given element.
     *
     * @param i_dt time step.
     * @param i_lp id of the limited plus element.
     * @param i_scSfSc sub-cells adjacent to sub-cells (sub-faces as bridge).
     * @param i_scTySf types of sub-faces adjacent to sub-cells.
     * @param i_iBndSt stencil for super-cells at internal boundaries.
     * @param i_admP admissibility of the element's previous solution.
     * @param i_scatter sub-cell scatter operator.
     * @param i_netUpSc net-updates if provided for the faces.
     * @param i_dofsDgP DOFs of the previous DG solution.
     * @param i_dofsScAdP sub-cell solution of the adjacent sub-cells (DG-faces as bridge).
     * @param io_dofsSc will be updated with the limited sub-cell solution.
     * @param i_solvSc sub-cell solver.
     *
     * @paramt TL_T_LID integral type of local ids.
     * @paramt TL_T_REAL floating point type.
     * @paramt TL_T_SP sparse type.
     * @paramt TL_T_SOLV_SC type of the sub-cell solver.
     **/
    template< typename TL_T_LID,
              typename TL_T_REAL,
              typename TL_T_SP,
              typename TL_T_SOLV_SC >
    static void limit( TL_T_REAL                      i_dt,
                       TL_T_LID                       i_lp,
                       unsigned short const           i_scSfSc[TL_N_SCS + TL_N_FAS * TL_N_SFS][TL_N_FAS],
                       unsigned short const           i_scTySf[TL_N_SCS][TL_N_FAS],
                       sc::ibnd::t_SuperCell<
                         TL_T_EL,
                         TL_O_SP >            const  &i_iBndSt,
                       TL_T_SP                        i_spTypes[TL_N_FAS],
                       bool                   const   i_admP[TL_N_CRS],
                       TL_T_REAL              const   i_scatter[TL_N_MDS_EL][TL_N_SCS],
                       TL_T_REAL              const (*i_netUpSc[TL_N_FAS])[TL_N_SFS][TL_N_CRS],
                       TL_T_REAL              const   i_dofsDgP[TL_N_QTS][TL_N_MDS_EL][TL_N_CRS],
                       TL_T_REAL              const (*i_dofsScAdP[TL_N_FAS])[TL_N_SFS][TL_N_CRS],
                       TL_T_REAL                      io_dofsSc[TL_N_QTS][TL_N_SCS][TL_N_CRS],
                       TL_T_SOLV_SC           const  &i_solvSc ) {
      // TODO: Use scratch memory
      TL_T_REAL l_subCell[TL_N_QTS][TL_N_SCS][TL_N_CRS];

      // derive element's sub-cell solution
      edge::sc::Kernels<
        TL_T_EL,
        TL_O_SP,
        TL_N_QTS,
        TL_N_CRS >::scatterReplaceVanilla( i_dofsDgP,
                                           i_scatter,
                                           io_dofsSc,
                                           i_admP,
                                           l_subCell );

      // assemble sub-grid, TODO: replace with strides in scatter ops
      TL_T_REAL l_sg[TL_N_QTS][TL_N_SCS + TL_N_FAS * TL_N_SFS][TL_N_CRS];
      for( unsigned short l_qt = 0; l_qt < TL_N_QTS; l_qt++ )
        for( unsigned short l_sc = 0; l_sc < TL_N_SCS; l_sc++ )
          for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ )
            l_sg[l_qt][l_sc][l_cr] = l_subCell[l_qt][l_sc][l_cr];

      for( unsigned short l_fa = 0; l_fa < TL_N_FAS; l_fa++ ) {
        for( unsigned short l_qt = 0; l_qt < TL_N_QTS; l_qt++ ) {
          for( unsigned short l_sf = 0; l_sf < TL_N_SFS; l_sf++ ) {
            for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {
              // default case: flux computation
              if( i_netUpSc[l_fa] == nullptr ) {
                if( i_dofsScAdP[l_fa] != nullptr ) {
                  l_sg[l_qt][TL_N_SCS + l_fa*TL_N_SFS + l_sf][l_cr] = i_dofsScAdP[l_fa][l_qt][l_sf][l_cr];
                }
              }
              // net-updates in ghost sub-cells
              else {
                l_sg[l_qt][TL_N_SCS + l_fa*TL_N_SFS + l_sf][l_cr] = i_netUpSc[l_fa][l_qt][l_sf][l_cr];
              }
            }
          }
        }
      }

      TL_T_BND_SC::apply( i_scSfSc,
                          i_spTypes,
                          l_sg );

      bool l_netUps[TL_N_FAS];
      for( unsigned short l_fa = 0; l_fa < TL_N_FAS; l_fa++ ) {
        if( i_netUpSc[l_fa] == nullptr )
          l_netUps[l_fa] = false;
        else l_netUps[l_fa] = true;
      }
      i_solvSc.tsSc( i_dt,
                     i_lp,
                     i_scSfSc,
                     i_scTySf,
                     l_netUps,
                     l_sg,
                     io_dofsSc );

      // apply super-cell stencils, if required
      for( unsigned short l_fa = 0; l_fa < TL_N_FAS; l_fa++ ) {
        if( l_netUps[l_fa] ) {
          ibnd::SuperCell< TL_T_EL,
                           TL_O_SP,
                           TL_N_QTS,
                           TL_N_CRS >::applyFa( i_iBndSt.col[l_fa],
                                                io_dofsSc );
        }
      }

    }

  public:
    /**
     * Applies the a-posteriori sub-cell limiter.
     *
     * @param i_first first limited plus element (limited + adjacent elements through faces).
     * @param i_nLps number of limited plus elements.
     * @param i_dt time step.
     * @param i_fIntL local flux integration matrices.
     * @param i_fIntN neighboring flux integration matrices.
     * @param i_fIntT transposed flux integration mattrices.
     * @param i_scOps sub-cell operators.
     * @param i_scConn sub-cell connectivity.
     * @param i_iBndSt super-cell stencils at internal boundaries.
     * @param i_charsFa face characteristics.
     * @param i_fsDg DG flux solvers for local contributions.
     * @param i_fsDgAd DG flux solvers for contributions of adjacent elements.
     * @param i_elFa DG faces adjacent to elements (no bridge).
     * @param i_elFaEl DG elements adjacent to elements (faces as bridge).
     * @param i_vIdElFaEl vertex ids of adjacent DG-elements (faces as bridge).
     * @param i_fIdElFaEl face ids of adjacent DG-elements (faces as bridge).
     * @param i_admP admissibility of the previous solution.
     * @param i_admC admissibility of the candidate solution.
     * @param o_admL0 will be set to admissibility of the limited solution.
     * @param o_admL1 will be set to admissibility of the limtied solution.
     * @param io_lock if set to true, the solution is locked in sub-cell space; will be reset to false.
     * @param io_limSync counter for number of limits since synchronization. Will be updated if appropiate.
     * @param i_tDofsDg temporary DOFs of the DG-solution.
     * @param io_dofsDg will be updated with DG-projection of the limited solution.
     * @param i_tDofsScP previous sub-cell solution for sub-cells adjacent to DG-faces.
     * @param i_tDofsScL will be set to sub-cell tDofs of limited solution.
     * @param io_dofsSc will be updated with sub-cell solution.
     * @param i_extP extrema of the previous solution.
     * @param o_extL will be set to extrema of the limited solution.
     * @param i_mm matrix multiplication kernels.
     * @param i_solvSc sub-cell solver.
     *
     * @paramt TL_T_LID integral type of local ids.
     * @paramt TL_T_REAL floating point type.
     * @paramt TL_T_CHARS_FA face characteristics, offering member .spType.
     * @paramt TL_T_MM type of the matrix-matrix multiplication kernels.
     * @paramt TL_T_SOLV_SC type of the sub-cell solver.
     **/
    template< typename TL_T_LID,
              typename TL_T_REAL,
              typename TL_T_CHARS_FA,
              typename TL_T_MM,
              typename TL_T_SOLV_SC >
    static void aPost( TL_T_LID                                        i_first,
                       TL_T_LID                                        i_nLps,
                       TL_T_REAL                                       i_dt,
                       TL_T_REAL              const                    i_fIntL[TL_N_FAS][TL_N_MDS_EL][TL_N_MDS_FA],
                       TL_T_REAL              const                    i_fIntN[TL_N_FAS][TL_N_MDS_EL][TL_N_MDS_FA],
                       TL_T_REAL              const                    i_fIntT[TL_N_FAS][TL_N_MDS_FA][TL_N_MDS_EL],
                       sc::t_ops<
                         TL_T_REAL,
                         TL_T_EL,
                         TL_O_SP >            const                   &i_scOps,
                       sc::t_connect<
                         TL_T_LID,
                         TL_T_EL,
                         TL_O_SP >            const                   &i_scConn,
                       sc::ibnd::t_SuperCell<
                         TL_T_EL,
                         TL_O_SP >            const                   &i_iBndSt,
                       TL_T_CHARS_FA          const                  (*i_charsFa),
                       TL_T_REAL              const                  (*i_fsDg)[TL_N_FAS][TL_N_QTS][TL_N_QTS],
                       TL_T_REAL              const                  (*i_fsDgAd)[TL_N_FAS][TL_N_QTS][TL_N_QTS],
                       TL_T_LID               const                  (*i_elFa)[TL_N_FAS],
                       TL_T_LID               const                  (*i_elFaEl)[TL_N_FAS],
                       unsigned short         const                  (*i_vIdElFaEl)[TL_N_FAS],
                       unsigned short         const                  (*i_fIdElFaEl)[TL_N_FAS],
                       bool                   const                  (*i_admP)[TL_N_CRS],
                       bool                                          (*i_admC)[TL_N_CRS], // TODO: change back to const (currenlty modified for iBnds)
                       bool                                          (*o_admL0)[TL_N_CRS],
                       bool                                          (*o_admL1)[TL_N_CRS],
                       bool                                          (*io_lock)[TL_N_CRS],
                       unsigned int                                  (*io_limSync)[TL_N_CRS],
#ifndef __INTEL_COMPILER
                       TL_T_REAL              const (* const * const   i_tDofsDg[2])[TL_N_MDS_EL][TL_N_CRS],
#else
                       TL_T_REAL                                    (**i_tDofsDg[2])[TL_N_MDS_EL][TL_N_CRS],
#endif
                       TL_T_REAL                                     (*io_dofsDg)[TL_N_QTS][TL_N_MDS_EL][TL_N_CRS],
                       TL_T_REAL                                  (* (*i_tDofsScP) [TL_N_FAS])[TL_N_QTS][TL_N_SFS][TL_N_CRS],
                       TL_T_REAL                                  (* (*o_tDofsScL) [TL_N_FAS])[TL_N_QTS][TL_N_SFS][TL_N_CRS],
                       TL_T_REAL                                     (*io_dofsSc)[TL_N_QTS][TL_N_SCS][TL_N_CRS],
                       TL_T_REAL                                     (*i_extP)[2][TL_N_QTS][TL_N_CRS],
                       TL_T_REAL                                     (*o_extL)[2][TL_N_QTS][TL_N_CRS],
                       TL_T_MM                const                   &i_mm,
                       TL_T_SOLV_SC           const                   &i_solvSc ) {
      // iterate over limited plus elements
      for( TL_T_LID l_lp = i_first; l_lp < i_first+i_nLps; l_lp++ ) {
        // dense id
        TL_T_LID l_el = i_scConn.lpEl[l_lp];

        // limited id
        TL_T_LID l_li = i_scConn.lpLi[l_lp];

        // derive the adjacent limited elements (faces as bridge)
        TL_T_LID l_lpFaLi[TL_N_FAS];
        for( unsigned short l_fa = 0; l_fa < TL_N_FAS; l_fa++ ) {
          if( i_scConn.lpFaLp[l_lp][l_fa] < std::numeric_limits< TL_T_LID >::max() )
            l_lpFaLi[l_fa] = i_scConn.lpLi[ i_scConn.lpFaLp[l_lp][l_fa] ];
          else
            l_lpFaLi[l_fa] = std::numeric_limits< TL_T_LID >::max();
        }


        // trigger sub-cell solver on both sides of internal boundary faces
        // TODO: Intermediate solution.
        // TODO: Consider using default-solver w/o net-updates if the fault is not slipping.
        // TODO: Consider doing the sync in the internal boundary solver (requires neighboring updates before rupture computations).
        for( unsigned short l_fa = 0; l_fa < TL_N_FAS; l_fa++ ) {
          TL_T_LID l_adFa = i_elFa[l_el][l_fa];
          if( ( i_charsFa[l_adFa].spType & t_spTypeElastic::RUPTURE ) == t_spTypeElastic::RUPTURE ) {
            for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {
              i_admC[l_li][l_cr] = i_admC[l_li][l_cr] && i_admC[ l_lpFaLi[l_fa] ][l_cr];
            }
          }
        }

        // iterate over faces and replace DG surface integral with sub-cell surface integral if required
        for( unsigned short l_fa = 0; l_fa < TL_N_FAS; l_fa++ ) {
          TL_T_LID l_liAd = l_lpFaLi[l_fa];
          TL_T_LID l_lpAd = i_scConn.lpFaLp[l_lp][l_fa];

          // only continue if the face is adjacent to a limited element
          if( l_liAd != std::numeric_limits< TL_T_LID >::max() ) {
            TL_T_LID l_elAd = i_elFaEl[l_el][l_fa];

            // id of the face of adjacent element w.r.t. reference element
            unsigned short l_fId = i_fIdElFaEl[l_el][l_fa];

            unsigned short l_fMatId  = i_vIdElFaEl[l_el][l_fa] * TL_N_FAS;
                           l_fMatId += i_fIdElFaEl[l_el][l_fa];

            // rollback DG surface integration of the face and replace with sub-cell integral
            surfIntRb( i_dt,
                       l_lp,
                       l_fa,
                       i_admC[l_liAd],
                       i_scConn.scDgAd[ i_vIdElFaEl[l_el][l_fa] ],
                       i_fIntL[l_fa],
                       i_fIntN[l_fMatId],
                       i_fIntT[l_fa],
                       i_scOps.sfInt[l_fa],
                       i_fsDg[l_el][l_fa],
                       i_fsDgAd[l_el][l_fa],
                       i_tDofsDg[0][l_el],
                       i_tDofsDg[0][l_elAd],
                     *(i_tDofsScP[l_lp][l_fa]),
                     *(i_tDofsScP[l_lpAd][l_fId]),
                       io_dofsDg[l_el],
                       i_mm,
                       i_solvSc );
          }
        }

        /*
         * perform limiting if required
         */
        if( l_li != std::numeric_limits< TL_T_LID >::max() ) {
          // derive if limiting is requ. and init admissibility of the limited solution
          bool l_lim = false;
          for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {
            l_lim = l_lim || (i_admC[l_li][l_cr] == false);
            o_admL0[l_li][l_cr] = i_admC[l_li][l_cr];
            o_admL1[l_li][l_cr] = i_admC[l_li][l_cr];
          }

          if( l_lim ) {
            // net-updates for internal boundaries
            TL_T_REAL const (*l_netUpSc[TL_N_FAS])[TL_N_SFS][TL_N_CRS];
            // sub-cell solution of adjacent element
            TL_T_REAL const (*l_dofsScAdP[TL_N_FAS])[TL_N_SFS][TL_N_CRS];

            int_spType l_spTypes[TL_N_FAS];

            // assemble info of adjacent (face as bridge) elements
            for( unsigned short l_fa = 0; l_fa < TL_N_FAS; l_fa++ ) {
              unsigned short l_fId = i_fIdElFaEl[l_el][l_fa];
              TL_T_LID l_adFa = i_elFa[l_el][l_fa];
              TL_T_LID l_adLp = i_scConn.lpFaLp[l_lp][l_fa];

              l_spTypes[l_fa] = i_charsFa[l_adFa].spType;

              // TODO: replace with generic ibnd sp-type
              if( ( l_spTypes[l_fa] & t_spTypeElastic::RUPTURE ) == t_spTypeElastic::RUPTURE ) {
                l_netUpSc[l_fa] = *(i_tDofsScP[l_lp][l_fa]);
              }
              else {
                l_netUpSc[l_fa] = nullptr;
              }

              // set adjacent admissibility and SC dofs
              if( l_adLp < std::numeric_limits< TL_T_LID >::max() ) {
                l_dofsScAdP[l_fa] = *(i_tDofsScP[l_adLp][l_fId]);
              }
              else {
                l_dofsScAdP[l_fa] = nullptr;
              }
            }

            // compute limited sub-cell solution
            limit( i_dt,
                   l_lp,
                   i_scConn.scSfSc,
                   i_scConn.scTySf,
                   i_iBndSt,
                   l_spTypes,
                   i_admP[l_li],
                   i_scOps.scatter,
                   l_netUpSc,
                   i_tDofsDg[1][l_el],
                   l_dofsScAdP,
                   io_dofsSc[l_li],
                   i_solvSc );

            // DG DOFs of the projected sub-cell solution
            TL_T_REAL l_dofsDg[TL_N_QTS][TL_N_MDS_EL][TL_N_CRS];
            // project the DG solution back
            edge::sc::Kernels< TL_T_EL,
                               TL_O_SP,
                               TL_N_QTS,
                               TL_N_CRS >::gatherVanilla( io_dofsSc[l_li],
                                                          i_scOps.gather,
                                                          l_dofsDg );
            // copy only relevant DOFs back
            for( unsigned short l_qt = 0; l_qt < TL_N_QTS; l_qt++ ) {
              for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {
                if( i_admC[l_li][l_cr] == false ) {
                  for( unsigned short l_md = 0; l_md < TL_N_MDS_EL; l_md++ ) {
                    io_dofsDg[l_el][l_qt][l_md][l_cr] = l_dofsDg[l_qt][l_md][l_cr];
                  }
                }
              }
            }

            // update admissibility based on DG-projection of the sub-cell solution
            bool l_admL[TL_N_CRS];
            TL_T_REAL l_scDgL[TL_N_QTS][TL_N_SCS][TL_N_CRS];
            TL_T_REAL l_exDgL[2][TL_N_QTS][TL_N_CRS];

            edge::sc::Kernels< TL_T_EL,
                               TL_O_SP,
                               TL_N_QTS,
                               TL_N_CRS >::dgExtremaVanilla(  io_dofsDg[l_el],
                                                              i_scOps.scatter,
                                                              l_scDgL,
                                                              l_exDgL[0],
                                                              l_exDgL[1] );

            edge::sc::Detections< TL_T_EL,
                                  TL_N_QTS,
                                  TL_N_CRS >::dmpFa( i_extP[l_lp],
                                                     i_extP,
                                                     l_exDgL,
                                                     i_scConn.lpFaLp[l_lp],
                                                     l_admL );

            for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {
              // lock solution in sub-cell space if required
              l_admL[l_cr] = l_admL[l_cr] && !io_lock[l_li][l_cr];

              // set admissibility of the limited solution
              if( i_admC[l_li][l_cr] == false ) {
                o_admL0[l_li][l_cr] = l_admL[l_cr];
                o_admL1[l_li][l_cr] = l_admL[l_cr];
              }
            }

            // update the extrema
            TL_T_REAL l_extSc[2][TL_N_QTS][TL_N_CRS];
            edge::sc::Kernels< TL_T_EL,
                               TL_O_SP,
                               TL_N_QTS,
                               TL_N_CRS >::scExtrema(  io_dofsSc[l_li],
                                                       l_extSc[0],
                                                       l_extSc[1] );

            for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {
              if( i_admC[l_li][l_cr] == false ) {
                // increase counter for number of limited solutions since synchronization
                io_limSync[l_li][l_cr]++;

                if( l_admL[l_cr] == false ) {
                  o_extL[l_lp][0][0][l_cr] = l_extSc[0][0][l_cr];
                  o_extL[l_lp][1][0][l_cr] = l_extSc[1][0][l_cr];
                }
                else {
                  o_extL[l_lp][0][0][l_cr] = l_exDgL[0][0][l_cr];
                  o_extL[l_lp][1][0][l_cr] = l_exDgL[1][0][l_cr];
                }
              }
            }


          }
        }


        // everything is updated: the final solution is either DG or SC: update sub-cell tDofs for next time step
        for( unsigned short l_fa = 0 ; l_fa < TL_N_FAS; l_fa++ ) {
          if( o_tDofsScL[l_lp][l_fa] != nullptr ) {
            unsigned short l_fMatId  = i_vIdElFaEl[l_el][l_fa] * TL_N_FAS;
                           l_fMatId += l_fa;
            edge::sc::Kernels<
              TL_T_EL,
              TL_O_SP,
              TL_N_QTS,
              TL_N_CRS >::scatterFaVanilla(   io_dofsDg[l_el],
                                              i_scOps.scatterSurf[TL_N_FAS+l_fMatId],
                                            *(o_tDofsScL[l_lp][l_fa]) );
          }
        }
        if( l_li != std::numeric_limits< TL_T_LID >::max() ) {
          // reset lock
          for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ )
            io_lock[l_li][l_cr] = false;

          edge::sc::Kernels< TL_T_EL,
                             TL_O_SP,
                             TL_N_QTS,
                             TL_N_CRS >::gatherSurfDofs( io_dofsSc[l_li],
                                                         i_scConn.faSfSc,
                                                         i_scConn.scDgAd,
                                                         i_vIdElFaEl[l_el],
                                                         o_tDofsScL[l_lp],
                                                         o_admL0[l_li] );
        }

      }
    }
};

#endif
