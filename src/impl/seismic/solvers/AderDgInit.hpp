/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2019, Alexander Breuer
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
 * Initialization of the ADER-DG data structures.
 **/
#ifndef EDGE_SEISMIC_SOLVERS_ADER_DG_INIT_HPP
#define EDGE_SEISMIC_SOLVERS_ADER_DG_INIT_HPP

#include "../setups/Elasticity.h"
#include "../setups/ViscoElasticity.h"
#include "mesh/common.hpp"

namespace edge {
  namespace seismic {
    namespace solvers {
      template< t_entityType TL_T_EL,
                bool         TL_MATS_SP >
      class AderDgInit;
    }
  }
}

/**
 * Init of ADER-DG data.
 *
 * @paramt TL_T_EL element type.
 * @paramt TL_MATS_SP true if the matrices are initialized as sparse.
 **/
template< t_entityType TL_T_EL,
          bool         TL_MATS_SP >
class edge::seismic::solvers::AderDgInit {
  private:
    //! number of dimensions
    static unsigned short const TL_N_DIS = C_ENT[TL_T_EL].N_DIM;

    //! number of vertices per element
    static unsigned short const TL_N_VES_EL = C_ENT[TL_T_EL].N_VERTICES;

    //! number of faces
    static unsigned short const TL_N_FAS = C_ENT[TL_T_EL].N_FACES;

    //! number of entries in the elastic star matrices
    static unsigned short const TL_N_ENS_STAR_E = (TL_MATS_SP) ? CE_N_ENS_STAR_E_SP( TL_N_DIS )
                                                               : CE_N_ENS_STAR_E_DE( TL_N_DIS );

    //! number of entries in the elastic flux solvers
    static unsigned short const TL_N_ENS_FS_E = CE_N_ENS_FS_E_DE( TL_N_DIS );

    //! number of entries in the anelastic source matrices
    static unsigned short const TL_N_ENS_SRC_A = (TL_MATS_SP) ? CE_N_ENS_SRC_A_SP( TL_N_DIS )
                                                              : CE_N_ENS_SRC_A_DE( TL_N_DIS );

    //! number of entries in the anelastic star matrices
    static unsigned short const TL_N_ENS_STAR_A = (TL_MATS_SP) ? CE_N_ENS_STAR_A_SP( TL_N_DIS )
                                                               : CE_N_ENS_STAR_A_DE( TL_N_DIS );

    //! number of entries in the anelastic flux solvers
    static unsigned short const TL_N_ENS_FS_A = CE_N_ENS_FS_A_DE( TL_N_DIS );

  public:
    /**
     * Initializes the anelastic sources matrices.
     *
     * @param i_nEls number of elements.
     * @param i_nRms number of relaxation mechanisms.
     * @param i_freqCen center frequency.
     * @param i_freqRat frequency ratio.
     * @param i_bgPars background parameters, lame parameters lambda and mu will be overwritten with elastic ones.
     * @param o_srcA will be set to anelastic source matrices.
     *
     * @paramt TL_T_LID local integral type.
     * @paramt TL_T_REAL floating point type.
     **/
    template< typename TL_T_LID,
              typename TL_T_REAL >
    static void initSrcA( TL_T_LID          i_nEls,
                          unsigned short    i_nRms,
                          double            i_freqCen,
                          double            i_freqRat,
                          t_bgPars        * io_bgPars,
                          TL_T_REAL      (* o_srcA)[TL_N_ENS_SRC_A] ) {
      double l_lameE[2] = { std::numeric_limits< double >::max(), std::numeric_limits< double >::max() };

      for( TL_T_LID l_el = 0; l_el < i_nEls; l_el++ ) {
        seismic::setups::ViscoElasticity::src( i_nRms,
                                                i_freqCen,
                                                i_freqRat,
                                                io_bgPars[l_el].qp,
                                                io_bgPars[l_el].qs,
                                                io_bgPars[l_el].lam,
                                                io_bgPars[l_el].mu,
                                                l_lameE[0],
                                                l_lameE[1],
                                                (o_srcA+l_el*std::size_t(i_nRms)) );

        // write unrelaxed lame parameters
        io_bgPars[l_el].lam = l_lameE[0];
        io_bgPars[l_el].mu  = l_lameE[1];
      }
    }

    /**
     * Initializes the star matrices.
     *
     * @param i_nEls number of elements.
     * @param i_elVe vertices adjacent to elements.
     * @param i_veChars vertex characteristics.
     * @param i_bgPars background parameters.
     * @param o_starE will be set to elastic star matrices.
     * @param o_starA will be set to anelastic star matrices.
     **/
    template< typename TL_T_LID,
              typename TL_T_REAL >
    static void initStar( TL_T_LID               i_nEls,
                          TL_T_LID      const (* i_elVe)[TL_N_VES_EL],
                          t_vertexChars const  * i_veChars,
                          t_bgPars      const  * i_bgPars,
                          TL_T_REAL           (* o_starE)[TL_N_DIS][TL_N_ENS_STAR_E],
                          TL_T_REAL           (* o_starA)[TL_N_DIS][TL_N_ENS_STAR_A] ) {
      for( TL_T_LID l_el = 0; l_el < i_nEls; l_el++ ) {
        // derive vertex coords
        double l_veCrds[TL_N_DIS][TL_N_VES_EL];
        mesh::common< TL_T_EL >::getElVeCrds( l_el,
                                              i_elVe,
                                              i_veChars,
                                              l_veCrds );

        // jacobian
        double l_jac[TL_N_DIS][TL_N_DIS];
        double l_jacInv[TL_N_DIS][TL_N_DIS];
        linalg::Mappings::evalJac( TL_T_EL, l_veCrds[0], l_jac[0] );

        // get inverse jacobian
        linalg::Matrix::inv( l_jac, l_jacInv );

        // init elastic star matrices
        double l_starE[TL_N_DIS][TL_N_ENS_STAR_E];
        seismic::setups::Elasticity::star( i_bgPars[l_el].rho,
                                           i_bgPars[l_el].lam,
                                           i_bgPars[l_el].mu,
                                           l_jacInv,
                                           l_starE );
        for( unsigned short l_di = 0; l_di < TL_N_DIS; l_di++ )
          for( unsigned short l_en = 0; l_en < TL_N_ENS_STAR_E; l_en++ )
            o_starE[l_el][l_di][l_en] = l_starE[l_di][l_en];

        // init anelastic star matrices
        if( o_starA != nullptr ) {
          double l_starA[TL_N_DIS][TL_N_ENS_STAR_A];
          seismic::setups::ViscoElasticity::star( l_jacInv,
                                                  l_starA );
          for( unsigned short l_di = 0; l_di < TL_N_DIS; l_di++ )
            for( unsigned short l_en = 0; l_en < TL_N_ENS_STAR_A; l_en++ )
              o_starA[l_el][l_di][l_en] = l_starA[l_di][l_en];
        }
      }
    }

    /**
     * Sets up the two-dimensional flux solvers for a single face.
     *
     * @param i_rhoL density of the left element.
     * @param i_rhoL density of the right element.
     * @param i_lamL Lame parameter lambda of the left element.
     * @param i_lamR Lame parameter lambda of the right element.
     * @param i_muL Lame parameter mu of the left element.
     * @param i_muR Lame parameter mu of the right element.
     * @param i_nx x-component of the face-normal.
     * @param i_ny y-component of the face-normal.
     * @param o_fsAl will be set to anelastic flux solver, applied to the left element's DOFs.
     * @param o_fsAr will be set to anelastic flux solver, applied to the right element's DOFs.
     * @param i_freeSurface true if free surface boundary conditions are applied to the right element.
     *
     * @paramt TL_T_REAL floating point type.
     **/
    template< typename TL_T_REAL >
    static void setUpFs( double    i_rhoL,
                         double    i_rhoR,
                         double    i_lamL,
                         double    i_lamR,
                         double    i_muL,
                         double    i_muR,
                         double    i_nx,
                         double    i_ny,
                         TL_T_REAL o_fsEl[5*5],
                         TL_T_REAL o_fsEr[5*5],
                         TL_T_REAL o_fsAl[3*5],
                         TL_T_REAL o_fsAr[3*5],
                         bool      i_freeSurface ) {
      // intermediate matrices
      double l_tE[5][5];
      double l_tA[3][3];
      double l_tm1[5][5];

      // compute trafos
      seismic::common::setupTrafo2d( i_nx, i_ny, l_tE );

      // extract "anelastic" back-rotation (stresses only)
      for( unsigned short l_q0 = 0; l_q0 < 3; l_q0++ )
        for( unsigned short l_q1 = 0; l_q1 < 3; l_q1++ )
          l_tA[l_q0][l_q1] = l_tE[l_q0][l_q1];

      seismic::common::setupTrafoInv2d( i_nx, i_ny, l_tm1 );

      // compute elastic flux solvers
      double l_fsE[2][5*5];
      setups::Elasticity::fs( i_rhoL,
                              i_rhoR,
                              i_lamL,
                              i_lamR,
                              i_muL,
                              i_muR,
                              l_tE,
                              l_tm1,
                              l_fsE[0],
                              l_fsE[1],
                              i_freeSurface );
      // store results
      for( unsigned short l_en = 0; l_en < 5*5; l_en++ ) {
        o_fsEl[l_en] = l_fsE[0][l_en];
        o_fsEr[l_en] = l_fsE[1][l_en];
      }

      if( o_fsAl != nullptr && o_fsAr != nullptr ) {
        // compute anelastic flux solvers
        double l_fsA[2][3*5];
        setups::ViscoElasticity::fs( i_rhoL,
                                     i_rhoR,
                                     i_lamL,
                                     i_lamR,
                                     i_muL,
                                     i_muR,
                                     l_tA,
                                     l_tm1,
                                     l_fsA[0],
                                     l_fsA[1],
                                     i_freeSurface );

        // store results
        for( unsigned short l_en = 0; l_en < 3*5; l_en++ ) {
          o_fsAl[l_en] = l_fsA[0][l_en];
          o_fsAr[l_en] = l_fsA[1][l_en];
        }
      }
    }

    /**
     * Sets up the three-dimensional flux solvers for a single face.
     *
     * @param i_rhoL density of the left element.
     * @param i_rhoL density of the right element.
     * @param i_lamL Lame parameter lambda of the left element.
     * @param i_lamR Lame parameter lambda of the right element.
     * @param i_muL Lame parameter mu of the left element.
     * @param i_muR Lame parameter mu of the right element.
     * @param i_nx x-component of the face-normal.
     * @param i_ny y-component of the face-normal.
     * @param i_nz z-component of the face-normal.
     * @param i_sx x-component of the first face-tangent.
     * @param i_sy y-component of the first face-tangent.
     * @param i_sz z-component of the first face-tangent.
     * @param i_tx x-component of the second face-tangent.
     * @param i_ty y-component of the second face-tangent.
     * @param i_tz z-component of the second face-tangent.
     * @param o_fsAl will be set to anelastic flux solver, applied to the left element's DOFs.
     * @param o_fsAr will be set to anelastic flux solver, applied to the right element's DOFs.
     * @param i_freeSurface true if free surface boundary conditions are applied to the right element.
     *
     * @paramt TL_T_REAL floating point type.
     **/
    template< typename TL_T_REAL >
    static void setUpFs( double    i_rhoL,
                         double    i_rhoR,
                         double    i_lamL,
                         double    i_lamR,
                         double    i_muL,
                         double    i_muR,
                         double    i_nx,
                         double    i_ny,
                         double    i_nz,
                         double    i_sx,
                         double    i_sy,
                         double    i_sz,
                         double    i_tx,
                         double    i_ty,
                         double    i_tz,
                         TL_T_REAL o_fsEl[9*9],
                         TL_T_REAL o_fsEr[9*9],
                         TL_T_REAL o_fsAl[6*9],
                         TL_T_REAL o_fsAr[6*9],
                         bool      i_freeSurface ) {
      // intermediate matrices
      double l_tE[9][9];
      double l_tA[6][6];
      double l_tm1[9][9];

      // compute trafos
      seismic::common::setupTrafo3d( i_nx, i_ny, i_nz,
                                     i_sx, i_sy, i_sz,
                                     i_tx, i_ty, i_tz,
                                     l_tE );

      // extract "anelastic" back-rotation (stresses only)
      for( unsigned short l_q0 = 0; l_q0 < 6; l_q0++ )
        for( unsigned short l_q1 = 0; l_q1 < 6; l_q1++ )
          l_tA[l_q0][l_q1] = l_tE[l_q0][l_q1];

      seismic::common::setupTrafoInv3d( i_nx, i_ny, i_nz,
                                        i_sx, i_sy, i_sz,
                                        i_tx, i_ty, i_tz,
                                        l_tm1 );

      // compute elastic flux solvers
      double l_fsE[2][9*9];
      setups::Elasticity::fs( i_rhoL,
                              i_rhoR,
                              i_lamL,
                              i_lamR,
                              i_muL,
                              i_muR,
                              l_tE,
                              l_tm1,
                              l_fsE[0],
                              l_fsE[1],
                              i_freeSurface );
      // store results
      for( unsigned short l_en = 0; l_en < 9*9; l_en++ ) {
        o_fsEl[l_en] = l_fsE[0][l_en];
        o_fsEr[l_en] = l_fsE[1][l_en];
      }

      if( o_fsAl != nullptr && o_fsAr != nullptr ) {
        // compute anelastic flux solvers
        double l_fsA[2][6*9];
        setups::ViscoElasticity::fs( i_rhoL,
                                     i_rhoR,
                                     i_lamL,
                                     i_lamR,
                                     i_muL,
                                     i_muR,
                                     l_tA,
                                     l_tm1,
                                     l_fsA[0],
                                     l_fsA[1],
                                     i_freeSurface );

        // store results
        for( unsigned short l_en = 0; l_en < 6*9; l_en++ ) {
          o_fsAl[l_en] = l_fsA[0][l_en];
          o_fsAr[l_en] = l_fsA[1][l_en];
        }
      }
    }

    /**
     * Initializes the flux solvers.
     * 
     * @param i_nEls number of elements.
     * @param i_nFas number of faces.
     * @param i_faEl elements adjacent to faces.
     * @param i_elVe vertices adjacent to elements.
     * @param i_elFaEl faces adjacent to elements.
     * @param i_veChars vertex characteristics.
     * @param i_faChars face characteristics.
     * @param i_elChars element characteristics.
     * @param i_bgPars background parameters.
     * @param o_fsA will be set to the anelastic flux solvers for the own and neighboring contributions.
     *
     * @paramt TL_T_LID local integral type.
     * @paramt TL_T_REAL floating point type.
     **/
    template< typename TL_T_LID,
              typename TL_T_REAL >
    static void initFs( TL_T_LID                i_nEls,
                        TL_T_LID                i_nFas,
                        TL_T_LID       const (* i_faEl)[2],
                        TL_T_LID       const (* i_elVe)[TL_N_VES_EL],
                        TL_T_LID       const (* i_elFaEl)[TL_N_FAS],
                        t_vertexChars  const  * i_veChars,
                        t_faceChars    const  * i_faChars,
                        t_elementChars const  * i_elChars,
                        t_bgPars       const  * i_bgPars,
                        TL_T_REAL            (* o_fsE[2])[TL_N_FAS][TL_N_ENS_FS_E],
                        TL_T_REAL            (* o_fsA[2])[TL_N_FAS][TL_N_ENS_FS_A] ) {
      PP_INSTR_FUN("flux_solvers")

      // init invalid to avoid silent errors
#ifdef PP_USE_OMP
#pragma omp parallel for
#endif
      for( TL_T_LID l_el = 0; l_el < i_nEls; l_el++ ) {
        for( TL_T_LID l_fa = 0; l_fa < TL_N_FAS; l_fa++ ) {
          for( unsigned short l_sd = 0; l_sd < 2; l_sd++ ) {
            // init elastic flux solvers
            for( unsigned short l_en = 0; l_en < TL_N_ENS_FS_E; l_en++ ) {
              o_fsE[l_sd][l_el][l_fa][l_en] = std::numeric_limits< TL_T_REAL >::max();
            }
            // init anelastic flux solvers
            if( o_fsA[0] != nullptr && o_fsA[1] != nullptr ) {
              for( unsigned short l_en = 0; l_en < TL_N_ENS_FS_A; l_en++ ) {
                o_fsA[l_sd][l_el][l_fa][l_en] = std::numeric_limits< TL_T_REAL >::max();
              }
            }
          }
        }
      }

      // do the init per face
#ifdef PP_USE_OMP
#pragma omp parallel for
#endif
      for( TL_T_LID l_fa = 0; l_fa < i_nFas; l_fa++ ) {
        // get left and right element
        TL_T_LID l_elL = i_faEl[l_fa][0];
        TL_T_LID l_elR = i_faEl[l_fa][1];

        // check if the neighbors exist, if not: bnd-conditions
        bool l_exL = l_elL != std::numeric_limits< TL_T_LID >::max();
        bool l_exR = l_elR != std::numeric_limits< TL_T_LID >::max();

        // find ids of neighboring flux solvers in the elements
        unsigned short l_fIdL = std::numeric_limits< unsigned short >::max();
        unsigned short l_fIdR = std::numeric_limits< unsigned short >::max();
        for( unsigned short l_fi = 0; l_fi < TL_N_FAS; l_fi++ ) {
          if( l_exL && i_elFaEl[l_elL][l_fi] == l_elR ) {
            // face should only by found once
            EDGE_CHECK_EQ( l_fIdL, std::numeric_limits< unsigned short >::max() );
            l_fIdL = l_fi;
          }
          if( l_exR && i_elFaEl[l_elR][l_fi] == l_elL ) {
            // face should only by found once
            EDGE_CHECK_EQ( l_fIdR, std::numeric_limits <unsigned short >::max() );
            l_fIdR = l_fi;
          }
        }
        EDGE_CHECK( !l_exL || l_fIdL != std::numeric_limits< unsigned short >::max() ) << l_fa;
        EDGE_CHECK( !l_exR || l_fIdR != std::numeric_limits< unsigned short >::max() ) << l_fa;

        // collect lame parameters
        TL_T_REAL l_rhoL, l_rhoR, l_lamL, l_lamR, l_muL, l_muR;
        if( l_exL ) {
          l_rhoL = i_bgPars[l_elL].rho; l_lamL = i_bgPars[l_elL].lam; l_muL = i_bgPars[l_elL].mu;
        }
        else {
          EDGE_CHECK( l_exR );
          // mirror right elements paramters for non-existing left element
          l_rhoL = i_bgPars[l_elR].rho; l_lamL = i_bgPars[l_elR].lam; l_muL = i_bgPars[l_elR].mu;
        }
        if( l_exR ) {
          l_rhoR = i_bgPars[l_elR].rho; l_lamR = i_bgPars[l_elR].lam; l_muR = i_bgPars[l_elR].mu;
        }
        else {
          EDGE_CHECK( l_exL );
          // mirror left elements paramters for non-existing right element
          l_rhoR = i_bgPars[l_elL].rho; l_lamR = i_bgPars[l_elL].lam; l_muR = i_bgPars[l_elL].mu;
        }

        // compute solvers for the left element
        if( l_exL ) {
          if( TL_N_DIS == 2 ) {
            setUpFs( l_rhoL, l_rhoR,
                     l_lamL, l_lamR,
                     l_muL,  l_muR,
                     i_faChars[l_fa].outNormal[0],
                     i_faChars[l_fa].outNormal[1],
                     o_fsE[0][l_elL][l_fIdL],
                     o_fsE[1][l_elL][l_fIdL],
                     (o_fsA[0] != nullptr) ? o_fsA[0][l_elL][l_fIdL] : nullptr,
                     (o_fsA[1] != nullptr) ? o_fsA[1][l_elL][l_fIdL] : nullptr,
                     (i_faChars[l_fa].spType & FREE_SURFACE) == FREE_SURFACE );
          }
          else if( TL_N_DIS == 3 ) {
            setUpFs( l_rhoL, l_rhoR,
                     l_lamL, l_lamR,
                     l_muL,  l_muR,
                     i_faChars[l_fa].outNormal[0],
                     i_faChars[l_fa].outNormal[1],
                     i_faChars[l_fa].outNormal[2],
                     i_faChars[l_fa].tangent0[0],
                     i_faChars[l_fa].tangent0[1],
                     i_faChars[l_fa].tangent0[2],
                     i_faChars[l_fa].tangent1[0],
                     i_faChars[l_fa].tangent1[1],
                     i_faChars[l_fa].tangent1[2],
                     o_fsE[0][l_elL][l_fIdL],
                     o_fsE[1][l_elL][l_fIdL],
                     (o_fsA[0] != nullptr) ? o_fsA[0][l_elL][l_fIdL] : nullptr,
                     (o_fsA[1] != nullptr) ? o_fsA[1][l_elL][l_fIdL] : nullptr,
                     (i_faChars[l_fa].spType & FREE_SURFACE) == FREE_SURFACE );
          }
        }

        // compute solvers for the right element
        if( l_exR ) {
          // boundary conditions have the element on the left-side per definition
          EDGE_CHECK_NE( (i_faChars[l_fa].spType & FREE_SURFACE), FREE_SURFACE );

          if( TL_N_DIS == 2 ) {
            setUpFs(  l_rhoR, l_rhoL,
                      l_lamR, l_lamL,
                      l_muR,  l_muL,
                     -i_faChars[l_fa].outNormal[0],
                     -i_faChars[l_fa].outNormal[1],
                      o_fsE[0][l_elR][l_fIdR],
                      o_fsE[1][l_elR][l_fIdR],
                      (o_fsA[0] != nullptr) ? o_fsA[0][l_elR][l_fIdR] : nullptr,
                      (o_fsA[1] != nullptr) ? o_fsA[1][l_elR][l_fIdR] : nullptr,
                      false );
          }
          else if( TL_N_DIS == 3 ) {
            setUpFs(  l_rhoR, l_rhoL,
                      l_lamR, l_lamL,
                      l_muR,  l_muL,
                     -i_faChars[l_fa].outNormal[0],
                     -i_faChars[l_fa].outNormal[1],
                     -i_faChars[l_fa].outNormal[2],
                      i_faChars[l_fa].tangent0[0],
                      i_faChars[l_fa].tangent0[1],
                      i_faChars[l_fa].tangent0[2],
                      i_faChars[l_fa].tangent1[0],
                      i_faChars[l_fa].tangent1[1],
                      i_faChars[l_fa].tangent1[2],
                      o_fsE[0][l_elR][l_fIdR],
                      o_fsE[1][l_elR][l_fIdR],
                      (o_fsA[0] != nullptr) ? o_fsA[0][l_elR][l_fIdR] : nullptr,
                      (o_fsA[1] != nullptr) ? o_fsA[1][l_elR][l_fIdR] : nullptr,
                      false );
          }
        }

        // derive vertex coords
        TL_T_REAL l_veCrds[2][TL_N_DIS][TL_N_VES_EL];
        if( l_exL ) {
          mesh::common< TL_T_EL >::getElVeCrds( l_elL,
                                                i_elVe,
                                                i_veChars,
                                                l_veCrds[0] );
        }
        if( l_exR ) {
          mesh::common< TL_T_EL >::getElVeCrds( l_elR,
                                                i_elVe,
                                                i_veChars,
                                                l_veCrds[1] );
        }

        // compute determinant of the mapping's jacobian
        TL_T_REAL l_jDet[2] = { std::numeric_limits< TL_T_REAL >::max(), std::numeric_limits< TL_T_REAL >::max() };

        TL_T_REAL l_jac[TL_N_DIS][TL_N_DIS];
        if( l_exL ) {
          linalg::Mappings::evalJac( TL_T_EL, l_veCrds[0][0], l_jac[0] );
          l_jDet[0] = linalg::Matrix::det( l_jac );
        }
        if( l_exR ) {
          linalg::Mappings::evalJac( TL_T_EL, l_veCrds[1][0], l_jac[0] );
          l_jDet[1] = linalg::Matrix::det( l_jac );
        }

        // ensure positive determinants
        EDGE_CHECK( !l_exL || l_jDet[0] > 0 );
        EDGE_CHECK( !l_exR || l_jDet[1] > 0 );
        EDGE_CHECK( i_faChars[l_fa].area > 0 );

        // compute scalar scaling of the solvers
        TL_T_REAL l_sca[2] = { std::numeric_limits< TL_T_REAL >::max(), std::numeric_limits< TL_T_REAL >::max() };
        if( l_exL ) {
          l_sca[0] = -i_faChars[l_fa].area / l_jDet[0];
          if( TL_T_EL == TET4 ) l_sca[0] *= 2;
        }
        if( l_exR ) {
          l_sca[1] = -i_faChars[l_fa].area / l_jDet[1];
          if( TL_T_EL == TET4 ) l_sca[1] *= 2;
        }

        // scale solvers and double-check, that we didn't mess up
        for( unsigned short l_en = 0; l_en < TL_N_ENS_FS_E; l_en++ ) {
            EDGE_CHECK( !l_exL || o_fsE[0][l_elL][l_fIdL][l_en] != std::numeric_limits< TL_T_REAL >::max() );
            EDGE_CHECK( !l_exR || o_fsE[0][l_elR][l_fIdR][l_en] != std::numeric_limits< TL_T_REAL >::max() );
            EDGE_CHECK( !l_exL || o_fsE[1][l_elL][l_fIdL][l_en] != std::numeric_limits< TL_T_REAL >::max() );
            EDGE_CHECK( !l_exR || o_fsE[1][l_elR][l_fIdR][l_en] != std::numeric_limits< TL_T_REAL >::max() );

            // scale left elements' solvers
            if( l_exL ) o_fsE[0][l_elL][l_fIdL][l_en] *= l_sca[0];
            if( l_exL ) o_fsE[1][l_elL][l_fIdL][l_en] *= l_sca[0];

            // scale right elements' solvers
            if( l_exR ) o_fsE[0][l_elR][l_fIdR][l_en] *= l_sca[1];
            if( l_exR ) o_fsE[1][l_elR][l_fIdR][l_en] *= l_sca[1];
        }

        if( o_fsA[0] != nullptr && o_fsA[1] != nullptr ) {
          for( unsigned short l_en = 0; l_en < TL_N_ENS_FS_A; l_en++ ) {
            EDGE_CHECK( !l_exL || o_fsA[0][l_elL][l_fIdL][l_en] != std::numeric_limits< TL_T_REAL >::max() );
            EDGE_CHECK( !l_exR || o_fsA[0][l_elR][l_fIdR][l_en] != std::numeric_limits< TL_T_REAL >::max() );
            EDGE_CHECK( !l_exL || o_fsA[1][l_elL][l_fIdL][l_en] != std::numeric_limits< TL_T_REAL >::max() );
            EDGE_CHECK( !l_exR || o_fsA[1][l_elR][l_fIdR][l_en] != std::numeric_limits< TL_T_REAL >::max() );

            // scale left elements' solvers
            if( l_exL ) o_fsA[0][l_elL][l_fIdL][l_en] *= l_sca[0];
            if( l_exL ) o_fsA[1][l_elL][l_fIdL][l_en] *= l_sca[0];

            // scale right elements' solvers
            if( l_exR ) o_fsA[0][l_elR][l_fIdR][l_en] *= l_sca[1];
            if( l_exR ) o_fsA[1][l_elR][l_fIdR][l_en] *= l_sca[1];
          }
        }
      }
    }

};

#endif