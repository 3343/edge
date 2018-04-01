/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
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
 * Upwind sub-cell solver.
 **/
#ifndef EDGE_SEISMIC_SC_UPWIND_HPP
#define EDGE_SEISMIC_SC_UPWIND_HPP

#include "constants.hpp"
#include "io/logging.h"
#include "data/Dynamic.h"
#include "../solvers/common.hpp"
#include "linalg/Matrix.h"
#include "linalg/Geom.hpp"

namespace edge {
  namespace elastic {
    namespace sc {
      template< unsigned short N_DIS >
      class UpWindSolver;

      template<>
      class UpWindSolver< 2 >;

      template<>
      class UpWindSolver< 3 >;

      template< typename       TL_T_REAL,
                t_entityType   TL_T_EL,
                unsigned short TL_O_SP,
                unsigned short TL_N_CRS >
      class UpWind;
    }
  }
}

/**
 * Sets up the face-local matrices for the two-dimensional upwind flux solver.
 **/
template<>
class edge::elastic::sc::UpWindSolver< 2 > {
  public:
    /**
     * Initializes the upwind flux solvers.
     *
     * @param i_mpL left element's material parameters: lambda, mu, rho.
     * @param i_mpR right element's material parameters: lambda, mu, rho.
     * @param i_n face-normal, pointing from the left to the right element.
     * @param i_t0 first face-tangent.
     * @param i_t1 second face-tangent, normal to the normal and t0.
     * @param o_upL will be set to solver for left element update [0]: own contribution, [1]: neighboring contribution.
     * @param o_upR will be set to solver for right element update [0]: own contribution, [1]: neighboring contribution.
     * @param i_vis scaling for the viscosity.
     * @param i_fs if true, waves 0 and 2 of the right solver are mirrored for free-surface boundaries.
     *
     * @paramt TL_T_REAL floating point precision.
     **/
    template< typename TL_T_REAL >
    static void lr( TL_T_REAL const i_mpL[3],
                    TL_T_REAL const i_mpR[3],
                    TL_T_REAL const i_n[2],
                    TL_T_REAL const *,
                    TL_T_REAL const *,
                    TL_T_REAL       o_upL[2][5][5],
                    TL_T_REAL       o_upR[2][5][5],
                    TL_T_REAL       i_vis,
                    bool            i_fs ) {

      solvers::common::setupSolver2d( i_mpL[2], i_mpR[2],
                                      i_mpL[0], i_mpR[0],
                                      i_mpL[1], i_mpR[1],
                                      i_n[0], i_n[1], 0,
                                      o_upL[0],
                                      o_upL[1],
                                      i_fs );

      solvers::common::setupSolver2d( i_mpL[2], i_mpR[2],
                                      i_mpL[0], i_mpR[0],
                                      i_mpL[1], i_mpR[1],
                                      -i_n[0], -i_n[1], 0,
                                      o_upR[0],
                                      o_upR[1] );

      // determine maximum signal speed of both sides
      TL_T_REAL l_mssEl =           elastic::common::getVelP( i_mpL[2], i_mpL[0], i_mpL[1] );
                l_mssEl = std::max( elastic::common::getVelP( i_mpR[2], i_mpR[0], i_mpR[1] ),
                                    l_mssEl );

      // get the trafos
      TL_T_REAL l_t[5][5];
      TL_T_REAL l_tm1[5][5];

      elastic::common::setupTrafo2d(  i_n[0],  i_n[1],
                                      l_t );

      elastic::common::setupTrafoInv2d( i_n[0],  i_n[1],
                                        l_tm1 );

      // derive scaling for artificial viscosity for zero-speed wave
      TL_T_REAL l_entFix[5][5];
      for( unsigned short l_ro = 0; l_ro < 5; l_ro++ )
        for( unsigned short l_co = 0; l_co < 5; l_co++ )
          l_entFix[l_ro][l_co] = 0;
      l_entFix[1][1]  = l_mssEl;
      l_entFix[1][1] *= i_vis;

      // do the matrix mults
      TL_T_REAL l_tmp[5][5];
      linalg::Matrix::matMulB0( 5, 5, 5,
                                l_t[0], l_entFix[0], l_tmp[0] );
      linalg::Matrix::matMulB0( 5, 5, 5,
                                l_tmp[0], l_tm1[0], l_entFix[0] );

      // add viscosity to solvers
      for( unsigned short l_ro = 0; l_ro < 5; l_ro++ ) {
        for( unsigned short l_co = 0; l_co < 5; l_co++ ) {
          o_upL[0][l_ro][l_co] += l_entFix[l_ro][l_co];
          o_upR[0][l_ro][l_co] += l_entFix[l_ro][l_co];

          o_upL[1][l_ro][l_co] -= l_entFix[l_ro][l_co];
          o_upR[1][l_ro][l_co] -= l_entFix[l_ro][l_co];
        }
      }

      // scale with -1 since we add the fluxes
      for( unsigned l_co = 0; l_co < 2; l_co++ ) {
        for( unsigned short l_q1 = 0; l_q1 < 5; l_q1++ ) {
          for( unsigned short l_q2 = 0; l_q2 < 5; l_q2++ ) {
            o_upL[l_co][l_q1][l_q2] *= (TL_T_REAL) -1.0;
            o_upR[l_co][l_q1][l_q2] *= (TL_T_REAL) -1.0;
          }
        }
      }

    }
};

/**
 * Sets up the face-local matrices for the three-dimensional upwind flux solver.
 **/
template<>
class edge::elastic::sc::UpWindSolver< 3 > {
  public:
    /**
     * Initializes the upwind flux solvers.
     *
     * @param i_mpL left element's material parameters: lambda, mu, rho.
     * @param i_mpR right element's material parameters: lambda, mu, rho.
     * @param i_n face-normal, pointing from the left to the right element.
     * @param i_t0 first face-tangent.
     * @param i_t1 second face-tangent, normal to the normal and t0.
     * @param o_upL will be set to solver for left element update [0]: own contribution, [1]: neighboring contribution.
     * @param o_upR will be set to solver for right element update [0]: own contribution, [1]: neighboring contribution.
     * @param i_vis scaling for the viscosity.
     * @param i_fs if true, waves 0, 3 and 5 of the left solver are mirrored for free-surface boundaries.
     * 
     * @paramt TL_T_REAL floating point precision.
     **/
    template< typename TL_T_REAL >
    static void lr( TL_T_REAL const i_mpL[3],
                    TL_T_REAL const i_mpR[3],
                    TL_T_REAL const i_n[3],
                    TL_T_REAL const i_t0[3],
                    TL_T_REAL const i_t1[3],
                    TL_T_REAL       o_upL[2][9][9],
                    TL_T_REAL       o_upR[2][9][9],
                    TL_T_REAL       i_vis,
                    bool            i_fs ) {
      solvers::common::setupSolver3d( i_mpL[2], i_mpR[2],
                                      i_mpL[0], i_mpR[0],
                                      i_mpL[1], i_mpR[1],
                                      i_n[0], i_n[1], i_n[2],
                                      i_t0[0], i_t0[1], i_t0[2],
                                      i_t1[0], i_t1[1], i_t1[2],
                                      o_upL[0],
                                      o_upL[1],
                                      i_fs );

      solvers::common::setupSolver3d( i_mpL[2], i_mpR[2],
                                      i_mpL[0], i_mpR[0],
                                      i_mpL[1], i_mpR[1],
                                      -i_n[0], -i_n[1], -i_n[2],
                                      i_t0[0], i_t0[1], i_t0[2],
                                      i_t1[0], i_t1[1], i_t1[2],
                                      o_upR[0],
                                      o_upR[1] );

      // determine maximum signal speed of both sides
      TL_T_REAL l_mssEl =           elastic::common::getVelP( i_mpL[2], i_mpL[0], i_mpL[1] );
                l_mssEl = std::max( elastic::common::getVelP( i_mpR[2], i_mpR[0], i_mpR[1] ),
                                    l_mssEl );

      // get the trafos
      TL_T_REAL l_t[9][9];
      TL_T_REAL l_tm1[9][9];

      elastic::common::setupTrafo3d(  i_n[0],  i_n[1],  i_n[2],
                                      i_t0[0], i_t0[1], i_t0[2],
                                      i_t1[0], i_t1[1], i_t1[2],
                                      l_t );

      elastic::common::setupTrafoInv3d( i_n[0],  i_n[1],  i_n[2],
                                        i_t0[0], i_t0[1], i_t0[2],
                                        i_t1[0], i_t1[1], i_t1[2],
                                        l_tm1 );

      // derive scaling for artificial viscosity for zero-speed waves
      TL_T_REAL l_entFix[9][9];
      for( unsigned short l_ro = 0; l_ro < 9; l_ro++ )
        for( unsigned short l_co = 0; l_co < 9; l_co++ )
          l_entFix[l_ro][l_co] = 0;
      l_entFix[1][1]  = l_mssEl;
      l_entFix[1][1] *= i_vis;
      l_entFix[2][2]  = l_entFix[1][1];

      // do the matrix mults
      TL_T_REAL l_tmp[9][9];
      linalg::Matrix::matMulB0( 9, 9, 9,
                                l_t[0], l_entFix[0], l_tmp[0] );
      linalg::Matrix::matMulB0( 9, 9, 9,
                                l_tmp[0], l_tm1[0], l_entFix[0] );

      // add viscosity to solvers
      for( unsigned short l_ro = 0; l_ro < 9; l_ro++ ) {
        for( unsigned short l_co = 0; l_co < 9; l_co++ ) {
          o_upL[0][l_ro][l_co] += l_entFix[l_ro][l_co];
          o_upR[0][l_ro][l_co] += l_entFix[l_ro][l_co];

          o_upL[1][l_ro][l_co] -= l_entFix[l_ro][l_co];
          o_upR[1][l_ro][l_co] -= l_entFix[l_ro][l_co];
        }
      }

      // scale with -1, since we add the fluxes
      for( unsigned l_co = 0; l_co < 2; l_co++ ) {
        for( unsigned short l_q1 = 0; l_q1 < 9; l_q1++ ) {
          for( unsigned short l_q2 = 0; l_q2 < 9; l_q2++ ) {
            o_upL[l_co][l_q1][l_q2] *= (TL_T_REAL) -1.0;
            o_upR[l_co][l_q1][l_q2] *= (TL_T_REAL) -1.0;
          }
        }
      }

    }
};

/**
 * Upwind sub-cell solver.
 *
 * @paramt TL_T_REAL floating point precision.
 * @paramt TL_T_EL element type.
 * @paramt TL_O_SP spatial order of the DG-solver.
 * @paramt TL_N_CRS number of fused simulations.
 **/
template< typename       TL_T_REAL,
          t_entityType   TL_T_EL,
          unsigned short TL_O_SP,
          unsigned short TL_N_CRS >
class edge::elastic::sc::UpWind {
  private:
    //! number of dimensions
    static const unsigned short TL_N_DIS = C_ENT[TL_T_EL].N_DIM;

    //! number of vertices per element
    static const unsigned short TL_N_VES = C_ENT[TL_T_EL].N_VERTICES;

    //! number of faces per element
    static const unsigned short TL_N_FAS = C_ENT[TL_T_EL].N_FACES;

    //! number of elastic quantities
    static const unsigned short TL_N_QTS = (TL_N_DIS == 2) ? 5 : 9;

    //! number of sub-faces per element face
    static unsigned short const TL_N_SFS = CE_N_SUB_FACES( TL_T_EL, TL_O_SP );

    //! number of sub-cells per element
    static unsigned short const TL_N_SCS = CE_N_SUB_CELLS( TL_T_EL, TL_O_SP );

    //! number of sub-face types
    static unsigned short const TL_N_TYSF = (TL_T_EL != TET4) ? TL_N_FAS : 6;

    // upwind solvers.
    // 0 - TL_N_TYSF-1: homogeneous (incorporates the element's material parameters), left update
    // TL_N_TYSF - 2*TL_N_TYSF-1: homogeneous, right update
    // 2*TL_N_TYSF-1 - 2*TL_N_TYSF-1+TL_N_FAS-1: heterogeneous (incorporates material parameters of the adjacent elements)
    // [][0]: local contribution
    // [][1]: neighboring contribution
    TL_T_REAL (* m_up)[TL_N_TYSF*2+TL_N_FAS][2][TL_N_QTS][TL_N_QTS];

  private:
    /**
     * Derives the upwinds solver's ids and side, based on the sub-cell's face-type.
     *
     * @param i_scTySf sub-cell's face type, DG-surf left: 0 - TL_N_FAS-1, DG-surf right: TL_N_FAS - 2*TL_N_FAS-1, Inner left: 2*TL_N_FAS - 2*TL_N_FAS+TL_N_TYSF-1, Inner right: 2*TL_N_FAS+TL_N_TYSF - 2*TL_N_FAS+2*TL_N_TYSF.
     * @return id of the pair of upwind solvers.
     **/
    static unsigned short upId( unsigned short i_scTySf ) {
      // DG surface or inner sub-face
      bool l_dgs = (i_scTySf < 2*TL_N_FAS) ? true : false;

      // contribution of the sub-cell
      unsigned short l_upId = (l_dgs) ?  i_scTySf : // implicitly assuming that there is no "right" for DG
                                        (i_scTySf - 2*TL_N_FAS);

      // switch between homogeneous and heterogeneous
      l_upId += (l_dgs) ? 2*TL_N_TYSF : 0;

      return l_upId;
    }

  public:
    /**
     * Allocates memory for the upwind solver.
     *
     * @param i_nLimPlus number of limited plus DG elements (limited + face neighbors).
     * @param io_dynMem will be called for dynamic memory allocation.
     *
     * @paramt TL_T_LID integral type of local ids.
     **/
    template< typename TL_T_LID >
    void alloc( TL_T_LID              i_nLimPlus,
                edge::data::Dynamic  &io_dynMem ) {
      // size of the homogenous and heterogeneous central flux solver per limited plus element
      std::size_t l_upSize = (std::size_t) (2*TL_N_TYSF+TL_N_FAS) * 2 * TL_N_QTS * TL_N_QTS;
                  l_upSize *= i_nLimPlus * sizeof(TL_T_REAL);
      m_up = (TL_T_REAL (*) [2*TL_N_TYSF+TL_N_FAS][2][TL_N_QTS][TL_N_QTS]) io_dynMem.allocate( l_upSize );
    }

    /**
     * Initializes the upwind sub-cell solver.
     *
     * @param i_firstLp first limited plus element.
     * @param i_sizeLp number of limited plus elements, which are initialized.
     * @param i_lpEl connectivity: limited plus element -> dense element.
     * @param i_gId global ids of the dense elements.
     * @param i_elVe connectivity: dense elemenet -> dense vertex.
     * @param i_elFa connectivity: dense element -> dense face.
     * @param i_elFaEl dense element -> dense element (faces as bridge).
     * @param i_charsVe characteristics of the vertices.
     * @param i_charsFa characteristics of the faces.
     * @param i_charsEl characteristics of the elements.
     * @param i_matPars material parameters.
     * @param i_vis amount of viscosity (0.5 is LLF w.r.t. zero-waves).
     *
     * @paramt TL_T_LID integral type for local ids.
     * @paramt TL_T_GID integral type for gloval ids.
     * @paramt TL_T_CHARS_VE characteristics of the vertices, providing member .coords.
     * @paramt TL_T_CHARS_FA characteristics of the faces, providing members .outNormal, .tangent0, tangent1, .area.
     * @paramt TL_T_CHARS_EL characteristics of the elements, providing member .volume.
     * @paramt TL_T_MAT_PARS material parameters, providing members .rho, .lam, and .mu (Lame parameters).
     **/
    template< typename TL_T_LID,
              typename TL_T_GID,
              typename TL_T_CHARS_VE,
              typename TL_T_CHARS_FA,
              typename TL_T_CHARS_EL,
              typename TL_T_MAT_PARS >
    void init( TL_T_LID                  i_firstLp,
               TL_T_LID                  i_sizeLp,
               TL_T_LID          const (*i_lpEl),
               TL_T_GID          const (*i_gId),
               TL_T_LID          const (*i_elVe)[TL_N_VES],
               TL_T_LID          const (*i_elFa)[TL_N_FAS],
               TL_T_LID          const (*i_elFaEl)[TL_N_FAS],
               TL_T_CHARS_VE     const (*i_charsVe),
               TL_T_CHARS_FA     const (*i_charsFa),
               TL_T_CHARS_EL     const (*i_charsEl),
               TL_T_MAT_PARS     const (*i_matPars),
               TL_T_REAL                 i_visDrHom = TL_T_REAL(0.025),
               TL_T_REAL                 i_visDrHet = TL_T_REAL(0.05) ) {
      // iterate over limited plus elements
      for( TL_T_LID l_lp = i_firstLp; l_lp < i_firstLp+i_sizeLp; l_lp++ ) {
        // derive dense id
        TL_T_LID l_el = i_lpEl[l_lp];

        // get material parameters of the element
        TL_T_REAL l_mpEl[3];
        l_mpEl[0] = i_matPars[l_el].lam;
        l_mpEl[1] = i_matPars[l_el].mu;
        l_mpEl[2] = i_matPars[l_el].rho;


        // assemble normals, tangents and face areas
        TL_T_REAL l_n[TL_N_TYSF][TL_N_DIS];
        TL_T_REAL l_t0[TL_N_TYSF][TL_N_DIS];
        TL_T_REAL l_t1[TL_N_TYSF][TL_N_DIS];
        TL_T_REAL l_ar[TL_N_TYSF];

        // determine DR viscosity
        TL_T_REAL l_drVisHom  = (i_charsEl[l_el].spType & RUPTURE) == RUPTURE;
                  l_drVisHom *= i_visDrHom;
        TL_T_REAL l_drVisHet[TL_N_FAS];
        for( unsigned short l_fa = 0; l_fa < TL_N_FAS; l_fa++ ) {
          l_drVisHet[l_fa] = 0;

          TL_T_LID l_elAd = i_elFaEl[l_el][l_fa];
          bool l_ruEl =     l_elAd != std::numeric_limits< TL_T_LID >::max()
                        && (i_charsEl[l_el].spType & RUPTURE) == RUPTURE;

          bool l_ruAd =     l_elAd != std::numeric_limits< TL_T_LID >::max()
                        && (i_charsEl[l_elAd].spType & RUPTURE) == RUPTURE;

          // add heterogeneous viscosity only for rupture elements and sub-cells adjacent to rupture elements
          if( l_ruAd || l_ruEl ) {
            l_drVisHet[l_fa] = i_visDrHet;
          }
        }

        for( unsigned short l_ty = 0; l_ty < TL_N_TYSF; l_ty++ ) {
          // default: sub-faces are parallel to element faces
          if( l_ty < TL_N_FAS ) {
            // get id of the the adjacent element
            TL_T_LID l_elAd = i_elFaEl[l_el][l_ty];

            // get face's id
            TL_T_LID l_faId = i_elFa[l_el][l_ty];

            for( unsigned short l_di = 0; l_di < TL_N_DIS; l_di++ ) {
              l_n[l_ty][ l_di] = i_charsFa[l_faId].outNormal[l_di];
              l_t0[l_ty][l_di] = i_charsFa[l_faId].tangent0[ l_di];
              l_t1[l_ty][l_di] = i_charsFa[l_faId].tangent1[ l_di];
            }

            // change direction of normal, if the local element is right (order of face -> element is global)
            if(    l_elAd != std::numeric_limits< TL_T_LID >::max()
                && i_gId[l_el] > i_gId[l_elAd] ) {
              for( unsigned short l_di = 0; l_di < TL_N_DIS; l_di++ )
                l_n[l_ty][l_di] *= -1;
            }

            l_ar[l_ty] = i_charsFa[l_faId].area;
          }
          // special handling for limiting-specific introduced sub-faces
          else {
            EDGE_CHECK_EQ( TL_T_EL, TET4 );

            // assemble tetrahedral vertices
            TL_T_REAL l_veCrds[3][4];
            for( unsigned short l_di = 0; l_di < 3; l_di++ )
              for( unsigned short l_ve = 0; l_ve < 4; l_ve++ )
                l_veCrds[l_di][l_ve] = i_charsVe[ i_elVe[l_el][l_ve] ].coords[l_di];

            // assemble normal and tangents
            linalg::Geom::sfAdd( TL_T_EL,
                                 (l_ty == TL_N_FAS) ? 0 : 1,
                                 l_veCrds[0],
                                 l_n[l_ty],
                                 l_t0[l_ty],
                                 l_t1[l_ty],
                                 l_ar[l_ty] );
          }
        }

        // compute unscaled upwind solvers
        for( unsigned short l_ty = 0; l_ty < TL_N_TYSF; l_ty++ ) {
          // default: sub-faces are parallel to element faces
          if( l_ty < TL_N_FAS ) {
            // get id of the the adjacent element
            TL_T_LID l_elAd = i_elFaEl[l_el][l_ty];

            // use own element, if at a boudary
            if( l_elAd == std::numeric_limits< TL_T_LID >::max() ) l_elAd = l_el;

            // get material parameters of the adjacent element
            TL_T_REAL l_mpElAd[3];
            l_mpElAd[0] = i_matPars[l_elAd].lam;
            l_mpElAd[1] = i_matPars[l_elAd].mu;
            l_mpElAd[2] = i_matPars[l_elAd].rho;

            // dummy solver for right element
            TL_T_REAL l_upDummy[2][TL_N_QTS][TL_N_QTS];

            // compute heterogeneous upwind solvers
            // free surface boundary conditions for the second solver in the case of undefined adjacency.
            // valid for outflow boundaries, since we use zero-valued ghost cells in that case.
            UpWindSolver< TL_N_DIS >::lr( l_mpEl,
                                          l_mpElAd,
                                          l_n[l_ty],
                                          l_t0[l_ty],
                                          l_t1[l_ty],
                                          m_up[l_lp][2*TL_N_TYSF+l_ty],
                                          l_upDummy,
                                          l_drVisHet[l_ty],
                                          (i_elFaEl[l_el][l_ty] == std::numeric_limits< TL_T_LID >::max())
                                        );
          }

          // compute homogeneous upwind solvers
          UpWindSolver< TL_N_DIS >::lr( l_mpEl,
                                        l_mpEl,
                                        l_n[l_ty],
                                        l_t0[l_ty],
                                        l_t1[l_ty],
                                        m_up[l_lp][l_ty],
                                        m_up[l_lp][TL_N_TYSF+l_ty],
                                        l_drVisHom,
                                        false );
        }


        // scale flux solvers and viscosity term by sub-cells' volume and  sub-faces' area
        for( unsigned short l_ty = 0; l_ty < TL_N_TYSF; l_ty++ ) {
          TL_T_REAL l_sca  = l_ar[l_ty] / TL_N_SFS;
                    l_sca *= TL_N_SCS / i_charsEl[l_el].volume;

          // upwind contribution
          for( unsigned short l_co = 0; l_co < 2; l_co++ ) {
            for( unsigned short l_q1 = 0; l_q1 < TL_N_QTS; l_q1++ ) {
              for( unsigned short l_q2 = 0; l_q2 < TL_N_QTS; l_q2++ ) {
                // homogeneous
                m_up[l_lp][              l_ty][l_co][l_q1][l_q2] *= l_sca;
                m_up[l_lp][    TL_N_TYSF+l_ty][l_co][l_q1][l_q2] *= l_sca;
                // heterogeneous
                if( l_ty < TL_N_FAS )
                  m_up[l_lp][2*TL_N_TYSF+l_ty][l_co][l_q1][l_q2] *= l_sca;
              }
            }
          }
        }

      }
    }

    /**
     * Performs a sub-cell time step.
     * 
     * Remark on net-updates:
     *   Optionally net-updates might be provided for sub-faces at the DG-faces.
     *   Internal dynamic rupture boundaries, providing net-updates in the case
     *   of a sliding fault, are an important use case for this feature.
     *
     *   If the flag for a DG-face is set to true, sub-cells adjacent to the DG
     *   face do not compute fluxes for the sub-faces discretizing the DG-face.
     *   Instead the net-updates are assumed to be stored in the corresponding
     *   ghost-sub-cells. These net-updates are added directly to the DOFs of
     *   the sub-cell. No scaling or processing in any way is performed. All
     *   other flux computations for the inner non-DG sub-faces are done through
     *   via the upwind solver (incl. LLF viscosity).
     *
     * @param i_dt time step.
     * @param i_lp id of the corresponding limited-plus element.
     * @param i_scSfSc sub-cells adjacent to sub-cells (sub-faces as bridge).
     * @param i_scTySf types of the the sub-cell's faces.
     * @param i_netUpSc optional switch applying net-updates for the sub-cells at the DG-faces.
     * @param i_dofsSc sub-cell DOFs at the current time step. 0 - N_SCS: sub-cells of the element, N_SCS - N_SCS + N_FAS*N_SFS-1: adjacent elements' sub-cell at DG-boundary.
     * @param o_dofsSc will be set to sub-cell DOFs at the next time step. 0 - N_SCS: sub-cells of the element.
     *
     * @paramt integral type for local ids.
     **/
    template< typename TL_T_LID >
    void tsSc( TL_T_REAL            i_dt,
               TL_T_LID             i_lp,
               unsigned short const i_scSfSc[TL_N_SCS + TL_N_FAS * TL_N_SFS][TL_N_FAS],
               unsigned short const i_scTySf[TL_N_SCS][TL_N_FAS],
               bool           const i_netUpSc[TL_N_FAS],
               TL_T_REAL            io_dofsSc[TL_N_QTS][TL_N_SCS+TL_N_FAS*TL_N_SFS][TL_N_CRS],
               TL_T_REAL            o_dofsSc[TL_N_QTS][TL_N_SCS][TL_N_CRS] ) const {
      // TODO: Currently, we compute the fluxes for each sub-cell.
      //       This could be reduced to one flux-computation per sub-face

      // iterate over the sub-cells
      for( unsigned short l_sc = 0; l_sc < TL_N_SCS; l_sc++ ) {
        // init solution
        for( unsigned short l_qt = 0; l_qt < TL_N_QTS; l_qt++ )
          for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ )
            o_dofsSc[l_qt][l_sc][l_cr] = io_dofsSc[l_qt][l_sc][l_cr];

        // iterate over faces of the sub-cell and compute net-updates
        for( unsigned short l_fa = 0; l_fa < TL_N_FAS; l_fa++ ) {
          // derive adjacent sub-cell
          TL_T_LID l_scAd = i_scSfSc[l_sc][l_fa];

          // derive type of the sub-face
          unsigned short l_ty = i_scTySf[l_sc][l_fa];

          // default case: flux computation
          if( l_ty > TL_N_FAS-1 || i_netUpSc[ l_ty%TL_N_FAS ] == false ) {
            // choose solver pair
            unsigned short l_upId = upId( l_ty );

            // iterate over target quantities
            for( unsigned short l_q1 = 0; l_q1 < TL_N_QTS; l_q1++ ) {
              // iterate over coupled quantities
              for( unsigned short l_q2  = 0; l_q2 < TL_N_QTS; l_q2++ ) {
                // iterate over fused runs
                for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {
                  // own flux contribution
                  o_dofsSc[l_q1][l_sc][l_cr] +=   io_dofsSc[l_q2][l_sc][l_cr]
                                                * m_up[i_lp][ l_upId ][0][l_q1][l_q2]
                                                * i_dt;

                  // neighboring flux contribution
                  o_dofsSc[l_q1][l_sc][l_cr] +=   io_dofsSc[l_q2][l_scAd][l_cr]
                                                * m_up[i_lp][ l_upId ][1][l_q1][l_q2]
                                                * i_dt;
                }
              }
            }
          }
          // net-updates are provided for the DG sub-face
          else {
            for( unsigned short l_q1 = 0; l_q1 < TL_N_QTS; l_q1++ )
              for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ )
                o_dofsSc[l_q1][l_sc][l_cr] += io_dofsSc[l_q1][l_scAd][l_cr];
          } 
        }
      }
    }

    /**
     * Computes the net-updates for the left element's sub-cells at a DG-face.
     * Remark:
     *  The net-updates are scaled by the area of the sub-faces and the inverse of the
     *  sub-cells' volume.
     *
     * @param i_dt time step.
     * @param i_lp id of the corresponding limited plus element.
     * @param i_fa respective face of the limited plus element.
     * @param i_dofsSc sub-cells DOFs adjacent to the DG-face. 0: limited plus (left), 1: right
     * @param o_netUps will be set to the net-updates for the limited plus element's sub-cells at the DG-face.
     **/
    template< typename TL_T_LID >
    void nuFaSf( TL_T_REAL      i_dt,
                 TL_T_LID       i_lp,
                 unsigned short i_fa,
                 TL_T_REAL      i_dofsSc[2][TL_N_QTS][TL_N_SFS][TL_N_CRS],
                 TL_T_REAL      o_netUps[TL_N_QTS][TL_N_SFS][TL_N_CRS] ) const {
      // iterate over sub-faces at DG-boundary
      for( unsigned short l_sf = 0; l_sf < TL_N_SFS; l_sf++ ) {

        // iterate over target quantities
        for( unsigned short l_q1 = 0; l_q1 < TL_N_QTS; l_q1++ ) {
          // init net-updates
          for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {
            // own sub-cell's viscosity contribution
            o_netUps[l_q1][l_sf][l_cr] = 0;
          }

          // iterate over coupled quantities
          for( unsigned short l_q2  = 0; l_q2 < TL_N_QTS; l_q2++ ) {
            // iterate over fused runs
            for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {
              // own flux contribution
              o_netUps[l_q1][l_sf][l_cr] +=   i_dofsSc[0][l_q2][l_sf][l_cr]
                                            * m_up[i_lp][2*TL_N_TYSF+i_fa][0][l_q1][l_q2]
                                            * i_dt;

              // neighboring flux contribution
              o_netUps[l_q1][l_sf][l_cr] +=   i_dofsSc[1][l_q2][l_sf][l_cr]
                                            * m_up[i_lp][2*TL_N_TYSF+i_fa][1][l_q1][l_q2]
                                            * i_dt;
            }
          }
        }

      }
    }
};

#endif
