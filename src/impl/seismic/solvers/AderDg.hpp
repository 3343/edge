/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *         Alexander Heinecke (alexander.heinecke AT intel.com)
 *
 * @section LICENSE
 * Copyright (c) 2019, Alexander Breuer
 * Copyright (c) 2016-2018, Regents of the University of California
 * Copyright (c) 2016, Intel Corporation
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
 * ADER-DG solver for the elastic wave equations.
 **/
#ifndef EDGE_SEISMIC_SOLVERS_ADER_DG_HPP
#define EDGE_SEISMIC_SOLVERS_ADER_DG_HPP

#include <limits>
#include "constants.hpp"
#include "mesh/common.hpp"
#include "impl/seismic/common.hpp"
#include "linalg/Matrix.h"
#include "linalg/Mappings.hpp"
#include "../kernels/Kernels.hpp"
#include "io/Receivers.h"
#include "InternalBoundary.hpp"
#include "FrictionLaws.hpp"
#include "sc/Kernels.hpp"
#include "sc/Detections.hpp"
#include "AderDgInit.hpp"

namespace edge {
  namespace seismic {
    namespace solvers { 
      template< typename       TL_T_REAL,
                unsigned short TL_N_RMS,
                t_entityType   TL_T_EL,
                unsigned short TL_O_SP,
                unsigned short TL_O_TI,
                unsigned short TL_N_CRS,
                bool           TL_MATS_SP >
      class AderDg;
    }
  }
}

/**
 * ADER-DG solver for the elastic wave equations, split into local and neighboring updates.
 *
 * @paramt TL_T_REAL floating point precision.
 * @paramt TL_N_RMS number of relaxation mechanisms.
 * @paramt TL_T_EL element type.
 * @paramt TL_O_SP spatial order.
 * @paramt TL_O_TI temporal order.
 * @paramt TL_N_CRS number of fused simulations.
 * @paramt TL_MATS_SPARSE true if the element-local matrices are initialized as sparse.
 **/
template< typename       TL_T_REAL,
          unsigned short TL_N_RMS,
          t_entityType   TL_T_EL,
          unsigned short TL_O_SP,
          unsigned short TL_O_TI,
          unsigned short TL_N_CRS,
          bool           TL_MATS_SP >
class edge::seismic::solvers::AderDg {
  private:
    //! number of dimensions
    static unsigned short const TL_N_DIS = C_ENT[TL_T_EL].N_DIM;

    //! number of vertices per face
    static unsigned short const TL_N_VES_FA = C_ENT[TL_T_EL].N_FACE_VERTICES;

    //! number of vertices per element
    static unsigned short const TL_N_VES_EL = C_ENT[TL_T_EL].N_VERTICES;

    //! number of faces
    static unsigned short const TL_N_FAS = C_ENT[TL_T_EL].N_FACES;

    //! number of DG modes
    static unsigned short const TL_N_MDS = CE_N_ELEMENT_MODES( TL_T_EL, TL_O_SP );

    //! number of elastic quantities
    static unsigned short const TL_N_QTS_E = CE_N_QTS_E( TL_N_DIS );

    //! number of quantities per relaxation mechanism
    static unsigned short const TL_N_QTS_M = CE_N_QTS_M( TL_N_DIS );

    //! number of sub-faces per DG-face
    static unsigned short const TL_N_SFS = CE_N_SUB_FACES( TL_T_EL, TL_O_SP );

    //! number of subcells per DG-element
    static unsigned short const TL_N_SCS = CE_N_SUB_CELLS( TL_T_EL, TL_O_SP );

    //! anelastic source matrices
    static unsigned short const TL_N_ENS_SRC_A = (TL_MATS_SP) ? CE_N_ENS_SRC_A_SP( TL_N_DIS )
                                                              : CE_N_ENS_SRC_A_DE( TL_N_DIS );
    TL_T_REAL (*m_srcA)[TL_N_ENS_SRC_A] = nullptr;

    //! anelastic star matrices
    static unsigned short const TL_N_ENS_STAR_A = (TL_MATS_SP) ? CE_N_ENS_STAR_A_SP( TL_N_DIS )
                                                               : CE_N_ENS_STAR_A_DE( TL_N_DIS );
    TL_T_REAL (*m_starA)[TL_N_DIS][TL_N_ENS_STAR_A] = nullptr;

    //! anelastic flux solvers
    static unsigned short const TL_N_ENS_FS_A = CE_N_ENS_FS_A_DE( TL_N_DIS );
    TL_T_REAL (*m_fsA[2])[TL_N_FAS][TL_N_ENS_FS_A] = { nullptr, nullptr };

    //! kernels
    kernels::Kernels< TL_T_REAL,
                      TL_N_RMS,
                      TL_T_EL,
                      TL_O_SP,
                      TL_O_TI,
                      TL_N_CRS > * m_kernels;

    /**
     * Allocates the constant data of the ADER-DG solver.
     *
     * @param i_nEls number of elements.
     * @parma i_align alignment of the individual arrays.
     * @param io_dynMem dynmic memory allocations.
     **/
    void alloc( std::size_t     i_nEls,
                std::size_t     i_align,
                data::Dynamic & io_dynMem ) {
      // size of the allocs in byte
      std::size_t l_size = std::numeric_limits< size_t >::max();

      if( TL_N_RMS > 0 ) {
        // anelastic source matrices
        l_size = i_nEls * std::size_t(TL_N_RMS) * TL_N_ENS_SRC_A ;
        l_size *= sizeof(TL_T_REAL);

        m_srcA = ( TL_T_REAL (*) [TL_N_ENS_SRC_A ] ) io_dynMem.allocate( l_size,
                                                                         i_align,
                                                                         false,
                                                                         true );

        // anelastic star matrices
        l_size = i_nEls * std::size_t(TL_N_DIS) * TL_N_ENS_STAR_A;
        l_size *= sizeof(TL_T_REAL);

        m_starA = ( TL_T_REAL (*) [TL_N_DIS][TL_N_ENS_STAR_A] ) io_dynMem.allocate( l_size,
                                                                                    i_align,
                                                                                    false,
                                                                                    true );

        // anelastic flux solvers
        l_size = i_nEls * std::size_t(TL_N_FAS) * TL_N_ENS_FS_A;
        l_size *= sizeof(TL_T_REAL);

        for( unsigned short l_sd = 0; l_sd < 2; l_sd++ ) {
          m_fsA[l_sd] = ( TL_T_REAL (*) [TL_N_FAS][TL_N_ENS_FS_A] ) io_dynMem.allocate( l_size,
                                                                                        i_align,
                                                                                        false,
                                                                                        true );
        }
      }
    }

  public:
    /**
     * Initializes the ADER-DG solver.
     *
     * @param i_nEls number of elements.
     * @param i_nFas number of faces.
     * @param i_faEl elements adjacent to faces.
     * @param i_elVe vertices adjacent to elements.
     * @param i_elFa faces adjacent to elements.
     * @param i_elMeDa mesh to data mapping for elements.
     * @param i_elDaMe data to mesh mapping for elements.
     * @param i_veChars vertex characteristics.
     * @param i_faChars face characteristics.
     * @param i_elChars elements characteristics.
     * @param io_bgPars background parameters (phase lambda and mu will be replaced with elastic lambda an mu if TL_N_RMS>0).
     * @param i_freqCen central frequency for attenuation.
     * @param i_freqRat frequency ratio between upper and lower frequencies for attenuation.
     * @param io_dynMem dynamic memory management.
     *
     * @paramt TL_T_LID integral type of local ids.
     */
    template< typename TL_T_LID >
    AderDg(TL_T_LID                i_nEls,
           TL_T_LID                i_nFas,
           TL_T_LID       const (* i_faEl)[2],
           TL_T_LID       const (* i_elVe)[TL_N_VES_EL],
           TL_T_LID       const (* i_elFa)[TL_N_FAS],
           TL_T_LID       const  * i_elMeDa,
           TL_T_LID       const  * i_elDaMe,
           t_vertexChars  const  * i_veChars,
           t_faceChars    const  * i_faChars,
           t_elementChars const  * i_elChars,
           t_bgPars              * io_bgPars,
           double                  i_freqCen,
           double                  i_freqRat,
           data::Dynamic         & io_dynMem ) {
      // alloc and init kernels
      TL_T_REAL *l_rfs = nullptr;
      if( TL_N_RMS > 0 ) {
        // get relaxation frequencies
        l_rfs = new TL_T_REAL[CE_MAX(2*TL_N_RMS-1, 1)];
        setups::ViscoElasticity::frequencies( TL_N_RMS,
                                              i_freqCen,
                                              i_freqRat,
                                              l_rfs );
        for( unsigned short l_rm = 0; l_rm < TL_N_RMS; l_rm++ ) {
          l_rfs[l_rm] = l_rfs[l_rm*2];
        }
      }
      m_kernels = new kernels::Kernels< TL_T_REAL,
                                        TL_N_RMS,
                                        TL_T_EL,
                                        TL_O_SP,
                                        TL_O_TI,
                                        TL_N_CRS >( l_rfs, io_dynMem );
      if( TL_N_RMS > 0 ) {
        delete[] l_rfs;
      }

      // allocate constant data
      alloc( i_nEls,
             ALIGNMENT.BASE.HEAP,
             io_dynMem );

      // init constant data
      if( TL_N_RMS > 0 ) {
        // init sources and replace lame parameters with elastic ones
        AderDgInit< TL_T_EL,
                    TL_MATS_SP >::initSrcA( i_nEls,
                                            TL_N_RMS,
                                            i_freqCen,
                                            i_freqRat,
                                            io_bgPars,
                                            m_srcA );

        // init star matrices
        AderDgInit< TL_T_EL,
                    TL_MATS_SP >::initStar( i_nEls,
                                            i_elVe,
                                            i_veChars,
                                            io_bgPars,
                                            m_starA );

        // init flux solvers
        AderDgInit< TL_T_EL,
                    TL_MATS_SP >::initFs( i_nEls,
                                          i_nFas,
                                          i_faEl,
                                          i_elVe,
                                          i_elFa,
                                          i_elMeDa,
                                          i_elDaMe,
                                          i_veChars,
                                          i_faChars,
                                          i_elChars,
                                          io_bgPars,
                                          m_fsA );
      }
    }

    /**
     * Destructor.
     **/
    ~AderDg() {
      delete m_kernels;
    }

    /**
     * Sets up the star matrices, which are a linear combination of the Jacobians.
     *
     * @param i_nElements number of elements.
     * @param i_vertexChars vertex characteristics.
     * @parma i_elVe vertices adjacent to the elements.
     * @param i_bgPars background parameters.
     * @param o_starMatrices will be set to star matrices.
     **/
    static void setupStarM(       int_el           i_nElements,
                            const t_vertexChars   *i_vertexChars,
                            const int_el         (*i_elVe)[TL_N_VES_EL],
                            const t_bgPars       (*i_bgPars)[1],
                                  t_matStar      (*o_starMatrices)[TL_N_DIS] ) {
      PP_INSTR_FUN("star_matrices")

      // iterate over elements
#ifdef PP_USE_OMP
#pragma omp parallel for
#endif
      for( int_el l_el = 0; l_el < i_nElements; l_el++ ) {
        // derive jacobians of the pdes
        real_base l_A[TL_N_DIS][TL_N_QTS_E][TL_N_QTS_E];

        for( unsigned int l_di = 0; l_di < TL_N_DIS; l_di++ ) {
#if defined PP_T_KERNELS_XSMM
          for( int_qt l_nz = 0; l_nz < N_MAT_STAR; l_nz++ ) o_starMatrices[l_el][l_di].mat[l_nz] = 0;
#else
          for( int_qt l_qt1 = 0; l_qt1 < TL_N_QTS_E; l_qt1++ ) {
            for( int_qt l_qt2 = 0; l_qt2 < TL_N_QTS_E; l_qt2++ ) {
              o_starMatrices[l_el][l_di].mat[l_qt1][l_qt2] = 0;
            }
          }
#endif
        }

        EDGE_CHECK( i_bgPars[l_el][0].rho > TOL.SOLVER );

        elastic::common::getJac( i_bgPars[l_el][0].rho,
                                 i_bgPars[l_el][0].lam,
                                 i_bgPars[l_el][0].mu,
                                 l_A[0][0],
                                 TL_N_DIS );

        // derive vertex coords
        real_mesh l_veCoords[TL_N_DIS][TL_N_VES_EL];
        mesh::common< TL_T_EL >::getElVeCrds( l_el, i_elVe, i_vertexChars, l_veCoords );

        // get inverse jacobian
        real_mesh l_jac[TL_N_DIS][TL_N_DIS];
        linalg::Mappings::evalJac( TL_T_EL, l_veCoords[0], l_jac[0] );

        real_mesh l_jacInv[TL_N_DIS][TL_N_DIS];
        linalg::Matrix::inv( l_jac, l_jacInv );

        // set star matrices
        // iterate over reference dimensions
        for( unsigned int l_d1 = 0; l_d1 < TL_N_DIS; l_d1++ ) {
#if defined PP_T_KERNELS_XSMM

#if PP_N_DIM == 2
          o_starMatrices[l_el][l_d1].mat[0] += l_A[0][0][3] * l_jacInv[0][l_d1];
          o_starMatrices[l_el][l_d1].mat[2] += l_A[0][1][3] * l_jacInv[0][l_d1];
          o_starMatrices[l_el][l_d1].mat[5] += l_A[0][2][4] * l_jacInv[0][l_d1];
          o_starMatrices[l_el][l_d1].mat[6] += l_A[0][3][0] * l_jacInv[0][l_d1];
          o_starMatrices[l_el][l_d1].mat[9] += l_A[0][4][2] * l_jacInv[0][l_d1];

          o_starMatrices[l_el][l_d1].mat[1] += l_A[1][0][4] * l_jacInv[1][l_d1];
          o_starMatrices[l_el][l_d1].mat[3] += l_A[1][1][4] * l_jacInv[1][l_d1];
          o_starMatrices[l_el][l_d1].mat[4] += l_A[1][2][3] * l_jacInv[1][l_d1];
          o_starMatrices[l_el][l_d1].mat[7] += l_A[1][3][2] * l_jacInv[1][l_d1];
          o_starMatrices[l_el][l_d1].mat[8] += l_A[1][4][1] * l_jacInv[1][l_d1];
#elif PP_N_DIM == 3
          o_starMatrices[l_el][l_d1].mat[ 0] += l_A[0][0][6] * l_jacInv[0][l_d1];
          o_starMatrices[l_el][l_d1].mat[ 3] += l_A[0][1][6] * l_jacInv[0][l_d1];
          o_starMatrices[l_el][l_d1].mat[ 6] += l_A[0][2][6] * l_jacInv[0][l_d1];
          o_starMatrices[l_el][l_d1].mat[10] += l_A[0][3][7] * l_jacInv[0][l_d1];
          o_starMatrices[l_el][l_d1].mat[14] += l_A[0][5][8] * l_jacInv[0][l_d1];
          o_starMatrices[l_el][l_d1].mat[15] += l_A[0][6][0] * l_jacInv[0][l_d1];
          o_starMatrices[l_el][l_d1].mat[19] += l_A[0][7][3] * l_jacInv[0][l_d1];
          o_starMatrices[l_el][l_d1].mat[23] += l_A[0][8][5] * l_jacInv[0][l_d1];

          o_starMatrices[l_el][l_d1].mat[ 1] += l_A[1][0][7] * l_jacInv[1][l_d1];
          o_starMatrices[l_el][l_d1].mat[ 4] += l_A[1][1][7] * l_jacInv[1][l_d1];
          o_starMatrices[l_el][l_d1].mat[ 7] += l_A[1][2][7] * l_jacInv[1][l_d1];
          o_starMatrices[l_el][l_d1].mat[ 9] += l_A[1][3][6] * l_jacInv[1][l_d1];
          o_starMatrices[l_el][l_d1].mat[12] += l_A[1][4][8] * l_jacInv[1][l_d1];
          o_starMatrices[l_el][l_d1].mat[16] += l_A[1][6][3] * l_jacInv[1][l_d1];
          o_starMatrices[l_el][l_d1].mat[18] += l_A[1][7][1] * l_jacInv[1][l_d1];
          o_starMatrices[l_el][l_d1].mat[22] += l_A[1][8][4] * l_jacInv[1][l_d1];

          o_starMatrices[l_el][l_d1].mat[ 2] += l_A[2][0][8] * l_jacInv[2][l_d1];
          o_starMatrices[l_el][l_d1].mat[ 5] += l_A[2][1][8] * l_jacInv[2][l_d1];
          o_starMatrices[l_el][l_d1].mat[ 8] += l_A[2][2][8] * l_jacInv[2][l_d1];
          o_starMatrices[l_el][l_d1].mat[11] += l_A[2][4][7] * l_jacInv[2][l_d1];
          o_starMatrices[l_el][l_d1].mat[13] += l_A[2][5][6] * l_jacInv[2][l_d1];
          o_starMatrices[l_el][l_d1].mat[17] += l_A[2][6][5] * l_jacInv[2][l_d1];
          o_starMatrices[l_el][l_d1].mat[20] += l_A[2][7][4] * l_jacInv[2][l_d1];
          o_starMatrices[l_el][l_d1].mat[21] += l_A[2][8][2] * l_jacInv[2][l_d1];
#else
#error not defined
#endif

#else
          for( unsigned int l_d2 = 0; l_d2 < TL_N_DIS; l_d2++ ) {
            for( int_qt l_qt1 = 0; l_qt1 < TL_N_QTS_E; l_qt1++ ) {
              for( int_qt l_qt2 = 0; l_qt2 < TL_N_QTS_E; l_qt2++ ) {
                o_starMatrices[l_el][l_d1].mat[l_qt1][l_qt2] += l_A[l_d2][l_qt1][l_qt2] * l_jacInv[l_d2][l_d1];
               }
             }
           }
#endif
         }
      }
    }

    /**
     * Local step: ADER + volume + local surface.
     *
     * @param i_first first element considered.
     * @param i_nElements number of elements.
     * @param i_time time of the initial DOFs.
     * @param i_dt time step.
     * @param i_firstSpRe first sparse receiver entity.
     * @param i_elChars element characteristics.
     * @param i_starM star matrices.
     * @param i_fluxSolvers flux solvers for the local element's contribution.
     * @param io_dofsE elastic DOFs.
     * @param io_dofsA anelastic DOFs.
     * @param o_tDofsDg will be set to temporary DOFs of the DG solution, [0]: time integrated, [1]: DOFs of previous time step (if required).
     * @param io_recvs will be updated with receiver info.
     *
     * @paramt TL_T_LID integer type of local entity ids.
     **/
    template < typename TL_T_LID >
    void local( TL_T_LID                             i_first,
                TL_T_LID                             i_nElements,
                double                               i_time,
                double                               i_dt,
                TL_T_LID                             i_firstSpRe,
                t_elementChars              const  * i_elChars,
                t_matStar                   const (* i_starM)[TL_N_DIS],
                t_fluxSolver                const (* i_fluxSolvers)[TL_N_FAS],
                TL_T_REAL                         (* io_dofsE)[TL_N_QTS_E][TL_N_MDS][TL_N_CRS],
                TL_T_REAL                         (* io_dofsA)[TL_N_MDS][TL_N_CRS],
                TL_T_REAL        (* const * const    o_tDofsDg[2])[TL_N_MDS][TL_N_CRS],
                edge::io::Receivers                & io_recvs ) const {
      // counter for receivers
      unsigned int l_enRe = i_firstSpRe;

      // temporary data structurre for product for two-way mult and receivers
      TL_T_REAL (*l_tmp)[TL_N_MDS][TL_N_CRS] = parallel::g_scratchMem->tRes[0];

      // buffer for derivatives
      TL_T_REAL (*l_derBuffer)[TL_N_QTS_E][TL_N_MDS][TL_N_CRS] = parallel::g_scratchMem->dBuf;

      // iterate over all elements
      for( TL_T_LID l_el = i_first; l_el < i_first+i_nElements; l_el++ ) {
        // store DOFs where required
        if( (i_elChars[l_el].spType & C_LTS_EL[EL_DOFS]) != C_LTS_EL[EL_DOFS] ) {}
        else {
          for( unsigned short l_qt = 0; l_qt < TL_N_QTS_E; l_qt++ )
            for( unsigned short l_md = 0; l_md < TL_N_MDS; l_md++ )
              for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ )
                o_tDofsDg[1][l_el][l_qt][l_md][l_cr] = io_dofsE[l_el][l_qt][l_md][l_cr];
        }

        // pointer to anelastic dofs
        TL_T_REAL (*l_dofsA)[TL_N_QTS_M][TL_N_MDS][TL_N_CRS] =
          (TL_T_REAL (*) [TL_N_QTS_M][TL_N_MDS][TL_N_CRS]) (io_dofsA+l_el*std::size_t(TL_N_RMS)*std::size_t(TL_N_QTS_M));

        TL_T_REAL l_tDofsA[CE_MAX(int(TL_N_RMS),1)][TL_N_QTS_M][TL_N_MDS][TL_N_CRS];
        TL_T_REAL l_derA[CE_MAX(int(TL_N_RMS),1)][TL_O_SP][TL_N_QTS_M][TL_N_MDS][TL_N_CRS];

        // compute ADER time integration
        m_kernels->m_time.ck( (TL_T_REAL)  i_dt,
                                         &(i_starM[l_el][0].mat), // TODO: fix struct
                                           m_starA[l_el],
                                           m_srcA+l_el*std::size_t(TL_N_RMS),
                                           io_dofsE[l_el],
                                           l_dofsA,
                                           l_tmp,
                                           l_derBuffer,
                                           l_derA,
                                           o_tDofsDg[0][l_el],
                                           l_tDofsA );

        // write receivers (if required)
        if( !( (i_elChars[l_el].spType & RECEIVER) == RECEIVER) ) {} // no receivers in the current element
        else { // we have receivers in the current element
          while( true ) { // iterate of possible multiple receiver-ouput per time step
            double l_rePt = io_recvs.getRecvTimeRel( l_enRe, i_time, i_dt );
            if( !(l_rePt >= 0) ) break;
            else {
              TL_T_REAL l_rePts = l_rePt;
              // eval time prediction at the given point
              m_kernels->m_time.evalTimePrediction(  1,
                                                   & l_rePts,
                                                     l_derBuffer,
                          (TL_T_REAL (*)[N_QUANTITIES][N_ELEMENT_MODES][N_CRUNS])l_tmp );

              // write this time prediction
              io_recvs.writeRecvAll( l_enRe, l_tmp );
            }
          }
          l_enRe++;
        }

        // compute volume integral
        m_kernels->m_volInt.apply( & i_starM[l_el][0].mat, // TODO: fix struct
                                     m_starA[l_el],
                                     m_srcA+l_el*std::size_t(TL_N_RMS),
                                     o_tDofsDg[0][l_el],
                                     l_tDofsA,
                                     io_dofsE[l_el],
                                     l_dofsA,
                                     l_tmp );

        // prefetches for next iteration
        const TL_T_REAL (* l_preDofs)[TL_N_MDS][TL_N_CRS] = nullptr;
        const TL_T_REAL (* l_preTint)[TL_N_MDS][TL_N_CRS] = nullptr;

        if( l_el < i_first+i_nElements-1 ) {
          l_preDofs = io_dofsE[l_el+1];
          l_preTint = o_tDofsDg[0][l_el+1];
        }
        else {
          l_preDofs = io_dofsE[l_el];
          l_preTint = o_tDofsDg[0][l_el];
        }

         /*
          * compute local surface contribution
          */
        // reuse derivative buffer
        TL_T_REAL (*l_tmpFa)[N_QUANTITIES][N_FACE_MODES][N_CRUNS] = parallel::g_scratchMem->tResSurf;
        // call kernel
        m_kernels->m_surfInt.local( ( TL_T_REAL (*)[TL_N_QTS_E][TL_N_QTS_E] ) ( i_fluxSolvers[l_el][0].solver[0] ), // TODO: fix struct
                                    m_fsA[0][l_el],
                                    o_tDofsDg[0][l_el],
                                    io_dofsE[l_el],
                                    l_dofsA,
                                    l_tmpFa,
                                    l_preDofs,
                                    l_preTint );
      }
    }

    /**
     * Solves rupture physics for the given faces.
     *
     * @param i_first first rupture faces.
     * @param i_nBf number of rupture faces at the internal boundary.
     * @param i_firstSpRe first sparse id of the receivers.
     * @param i_time current time of the faces.
     * @param i_dt time step of the two adjacent elements.
     * @param i_scDgAd adjacency of sub-cells at the faces of face-adjacent elements.
     * @parma i_liDoLiDu possibly duplicated limited elements, adjacent to an the dominant one.
     * @param i_starM star matrices.
     * @param i_iBnd internal boundary data.
     * @param i_frictionGlobal global data of the friction law.
     * @param i_frictionFa face-local data of the friction law.
     * @param i_frictionSf data of the friction law local to the sub-faces of the DG-faces.
     * @param io_tDofs sub-cell tDOFs, which will be updated with net-updates.
     * @param io_admC admissiblity of the candidate solution, will be set to false if any of the DG-faces are rupture faces and the fault failed.
     * @param io_lock lock for the sub-cell solution. Will be set to true for all active rupture elements.
     * @param io_recvsSf receivers at sub-faces.
     * @param i_mm libxsmm kernels
     *
     * @paramt TL_T_LID integral type of local entity ids.
     * @paramt TL_T_REAL type used for floating point arithmetic.
     * @paramt TL_T_FRI_GL struct representing global friction data.
     * @paramt TL_T_FRI_FA struct representing face-local friction data.
     * @paramt TL_T_FRI_SF struct representing sub-face-local friction data.
     * @paramt TL_T_SP integer type of the sparse type.
     * @paramt TL_T_RECV_SF type of the sub-face receivers.
     **/
    template< typename TL_T_LID,
              typename TL_T_FRI_GL,
              typename TL_T_FRI_FA,
              typename TL_T_FRI_SF,
              typename TL_T_SP,
              typename TL_T_RECV_SF,
              typename TL_T_MM >
    static void rupture( TL_T_LID                                       i_first,
                         TL_T_LID                                       i_nBf,
                         TL_T_LID                                       i_firstSpRe,
                         TL_T_REAL                                      i_time,
                         TL_T_REAL                                      i_dt,
                         unsigned short                      const      i_scDgAd[TL_N_VES_FA][TL_N_SFS],
                         TL_T_LID                            const   (* i_liDoLiDu)[TL_N_FAS],
                         TL_T_LID                            const   (* i_liLp),
                         t_matStar                           const   (* i_starM)[TL_N_DIS],
                         edge::sc::ibnd::t_InternalBoundary<
                           TL_T_LID,
                           TL_T_REAL,
                           TL_T_SP,
                           TL_T_EL,
                           TL_N_QTS_E >                      const     & i_iBnd,
                         TL_T_FRI_GL                                   & i_frictionGlobal,
                         TL_T_FRI_FA                                   * i_frictionFa,
                         TL_T_FRI_SF                                  (* i_frictionSf)[TL_N_SFS],
                         TL_T_REAL                                  (* (*io_tDofs) [TL_N_FAS])[TL_N_QTS_E][TL_N_SFS][TL_N_CRS],
                         bool                                         (* io_admC)[TL_N_CRS],
                         bool                                         (* io_lock)[TL_N_CRS],
                         TL_T_RECV_SF                                  & io_recvsSf,
                         TL_T_MM                              const    & i_mm
                         ) {
      // store sparse receiver id
      TL_T_LID l_faRe = i_firstSpRe;

      // struct for the pertubation of the middle states
      struct {
        TL_T_FRI_GL  *gl;
        TL_T_FRI_FA  *fa;
        TL_T_FRI_SF (*sf)[TL_N_SFS];
      } l_faData;
      l_faData.gl = &i_frictionGlobal;
      l_faData.fa =  i_frictionFa+i_first;
      l_faData.sf =  i_frictionSf+i_first;

      // iterate over internal boundary faces
      for( TL_T_LID l_bf = i_first; l_bf < i_first+i_nBf; l_bf++ ) {
        // get limited elements
        TL_T_LID const *l_li = i_iBnd.connect.bfLe[l_bf];

        // get corresponding limited plus elements
        TL_T_LID l_lp[2];
        for( unsigned short l_sd = 0; l_sd < 2; l_sd++ ) l_lp[l_sd] = i_liLp[ l_li[l_sd] ];

        // get local face ids and vertex id
        unsigned short const *l_fIdBfEl = i_iBnd.bfChars[l_bf].fIdBfEl;
        unsigned short const l_vIdFaEl = i_iBnd.bfChars[l_bf].vIdFaElR;

        // get sub-cell solution at DG-faces
        TL_T_REAL (*l_dofsSc)[TL_N_QTS_E][TL_N_SFS][TL_N_CRS] = parallel::g_scratchMem->scFa;
        for( unsigned short l_qt = 0; l_qt < TL_N_QTS_E; l_qt++ ) {
          for( unsigned short l_sf = 0; l_sf < TL_N_SFS; l_sf++ ) {
            unsigned short l_sfRe = i_scDgAd[l_vIdFaEl][l_sf];
#pragma omp simd
            for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {
              l_dofsSc[0][l_qt][l_sf][l_cr] = (*io_tDofs[ l_lp[0] ][ l_fIdBfEl[0] ])[l_qt][l_sfRe][l_cr];
              l_dofsSc[1][l_qt][l_sf][l_cr] = (*io_tDofs[ l_lp[1] ][ l_fIdBfEl[1] ])[l_qt][l_sf  ][l_cr];
            }
          }
        }

        // rupture indicator for the sub-faces
        bool l_rup[TL_N_SFS][TL_N_CRS];

        // net-updates
        TL_T_REAL (*l_netUps)[TL_N_QTS_E][TL_N_SFS][TL_N_CRS] = parallel::g_scratchMem->scFa+2;

        // call the rupture solver
        edge::elastic::solvers::InternalBoundary<
          TL_T_EL,
          TL_N_QTS_E,
          TL_O_SP,
          TL_N_CRS >::template netUpdates<
            TL_T_REAL, TL_T_MM, 
            edge::elastic::solvers::FrictionLaws< TL_N_DIS, TL_N_CRS >
          > (  i_iBnd.mss[l_bf][0],
               i_iBnd.mss[l_bf][1],
               i_iBnd.mss[l_bf][2],
               i_iBnd.mss[l_bf][3],
               l_dofsSc[0],
               l_dofsSc[1],
               l_netUps[0],
               l_netUps[1],
               l_rup,
               i_mm,
               i_dt,
              &l_faData );
 
        // store net-updates
        for( unsigned short l_qt = 0; l_qt < TL_N_QTS_E; l_qt++ ) {
          for( unsigned short l_sf = 0; l_sf < TL_N_SFS; l_sf++ ) {
            // sub-cells are ordered based on the left element; reorder for right element
            unsigned short l_sfRe = i_scDgAd[l_vIdFaEl][l_sf];

#pragma omp simd
            for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {
               (*io_tDofs[ l_lp[0] ][ l_fIdBfEl[0] ])[l_qt][l_sf][l_cr] = l_netUps[0][l_qt][l_sf][l_cr];
               (*io_tDofs[ l_lp[1] ][ l_fIdBfEl[1] ])[l_qt][l_sf][l_cr] = l_netUps[1][l_qt][l_sfRe][l_cr];
            }
          }
        }

        // update admissibility
        for( unsigned short l_sf = 0; l_sf < TL_N_SFS; l_sf++ ) {
#pragma omp simd
          for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {
            // enforce sub-cell solver for DR
            io_admC[ l_li[0] ][l_cr] = false;
            io_admC[ l_li[1] ][l_cr] = false;

            io_lock[ l_li[0] ][l_cr] = true;
            io_lock[ l_li[1] ][l_cr] = true;
          }
        }

        // synchronize possible duplicates
        for( unsigned short l_sd = 0; l_sd < 2; l_sd++ ) {
          for( unsigned short l_fa = 0; l_fa < TL_N_FAS; l_fa++ ) {
            TL_T_LID l_liDu = i_liDoLiDu[ l_li[l_sd] ][l_fa];

            if( l_liDu == std::numeric_limits< TL_T_LID >::max() ) break;
            else {
              // copy over admissibility and lock
              for( unsigned short l_sf = 0; l_sf < TL_N_SFS; l_sf++ ) {
#pragma omp simd
                for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {
                  io_admC[ l_liDu ][l_cr] = false;
                  io_lock[ l_liDu ][l_cr] = true;
                }
              }

              TL_T_LID l_lpDu = i_liLp[ l_liDu ];

              // copy over data
              for( unsigned short l_qt = 0; l_qt < TL_N_QTS_E; l_qt++ ) {
                for( unsigned short l_sf = 0; l_sf < TL_N_SFS; l_sf++ ) {
                  // sub-cells are ordered based on the left element; reorder for right element
                  unsigned short l_sfRe = (l_sd == 0) ? l_sf : i_scDgAd[l_vIdFaEl][l_sf];

#pragma omp simd
                  for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ )
                    (*io_tDofs[ l_lpDu ][ l_fIdBfEl[l_sd] ])[l_qt][l_sf][l_cr] = l_netUps[l_sd][l_qt][l_sfRe][l_cr];
                }
              }
            }
          }
        }

        // check if this face requires receiver output
        if( (i_iBnd.bfChars[l_bf].spType & RECEIVER) != RECEIVER ){}
        else {
          // check if the receiver requires output
          if( io_recvsSf.getRecvTimeRel( l_faRe, i_time, i_dt ) >= -TOL.TIME ) {
            // gather receiver data, TODO: outsource
            TL_T_REAL l_buff[ (TL_N_DIS-1)*3 ][TL_N_SFS][TL_N_CRS];

            for( unsigned short l_sf = 0; l_sf < TL_N_SFS; l_sf++ ) {
              for( unsigned short l_di = 0; l_di < TL_N_DIS-1; l_di++ ) {
#pragma omp simd
                for( unsigned short l_ru = 0; l_ru < TL_N_CRS; l_ru++ ) {
                  l_buff[             0+l_di][l_sf][l_ru] = i_frictionSf[l_bf][l_sf].tr[l_di][l_ru];
                  l_buff[  (TL_N_DIS-1)+l_di][l_sf][l_ru] = i_frictionSf[l_bf][l_sf].sr[l_di][l_ru];
                  l_buff[2*(TL_N_DIS-1)+l_di][l_sf][l_ru] = i_frictionSf[l_bf][l_sf].dd[l_di][l_ru];
                }
              }
            }

            // write the receiver info
            io_recvsSf.writeRecvAll( i_time, i_dt, l_faRe, l_buff );
          }

          l_faRe++;
        }

        // update pointers
        l_faData.fa++;
        l_faData.sf++;
      }
    }

    /**
     * Performs the neighboring updates of the ADER-DG scheme.
     *
     * @param i_first first element considered.
     * @param i_nElements number of elements.
     * @param i_firstLi first limited element.
     * @param i_firstLp first limited plus element.
     * @param i_firstEx first element computing extrema.
     * @param i_scatter scatter opterators (DG -> sub-cells).
     * @param i_faChars face characteristics.
     * @param i_elChars element characteristics.
     * @param i_fluxSolvers flux solvers for the neighboring elements' contribution.
     * @param i_lpFaLp limited plus elements adjacent to limited plus (faces as bridge).
     * @param i_elFa elements' adjacent faces.
     * @param i_elFaEl face-neighboring elements.
     * @param i_fIdElFaEl local face ids of face-neighboring elememts.
     * @param i_vIdElFaEl local vertex ids w.r.t. the shared face from the neighboring elements' perspsective.
     * @param i_tDofsDg temporarary DG DOFs ([0]: time integrated, [1]: DOFs of previous time step).
     * @param io_dofs DOFs which will be updated with neighboring elements' contribution.
     * @param io_admC will be updated with the admissibility of the candidate solution.
     * @param i_extP extreme of the previous solution.
     * @param o_extC will be set to extreme of the candidate solution.
     * @param i_mm matrix-matrix multiplication kernels.
     *
     * @paramt TL_T_LID integer type of local entity ids.
     * @paramt TL_T_MM type of the matrix-matrix multiplication kernels.
     **/
    template< typename TL_T_LID,
              typename TL_T_MM >
    void neigh( TL_T_LID                              i_first,
                TL_T_LID                              i_nElements,
                TL_T_LID                              i_firstLi,
                TL_T_LID                              i_firstLp,
                TL_T_LID                              i_firstEx,
                TL_T_REAL      const                  i_scatter[TL_N_MDS][TL_N_SCS],
                t_faceChars    const                * i_faChars,
                t_elementChars const                * i_elChars,
                t_fluxSolver   const               (* i_fluxSolvers)[TL_N_FAS],
                unsigned short const                  i_faSfSc[TL_N_FAS][TL_N_SFS],
                unsigned short const                  i_scDgAd[TL_N_VES_FA][TL_N_SFS],
                TL_T_LID       const               (* i_lpFaLp)[TL_N_FAS],
                TL_T_LID       const               (* i_elFa)[TL_N_FAS],
                TL_T_LID       const               (* i_elFaEl)[TL_N_FAS],
                unsigned short const               (* i_fIdElFaEl)[TL_N_FAS],
                unsigned short const               (* i_vIdElFaEl)[TL_N_FAS],
                TL_T_REAL            (* const * const i_tDofsDg[2])[TL_N_MDS][TL_N_CRS],
                TL_T_REAL                          (* io_dofsE)[TL_N_QTS_E][TL_N_MDS][TL_N_CRS],
                TL_T_REAL                          (* io_dofsA)[TL_N_MDS][TL_N_CRS],
                bool                               (* io_admC)[TL_N_CRS],
                TL_T_REAL      const               (* i_extP)[2][TL_N_QTS_E][TL_N_CRS],
                TL_T_REAL                          (* o_extC)[2][TL_N_QTS_E][TL_N_CRS],
                TL_T_REAL                        (* (*o_tDofsSc) [TL_N_FAS])[TL_N_QTS_E][TL_N_SFS][TL_N_CRS],
                TL_T_MM        const                & i_mm ) const {
      // counter for elements computing extrema
      TL_T_LID l_ex = i_firstEx;

      // counter for limited plus elements
      TL_T_LID l_lp = i_firstLp;

      // counter for limited elements
      TL_T_LID l_li = i_firstLi;

      // temporary product for three-way mult
      TL_T_REAL (*l_tmpFa)[N_QUANTITIES][N_FACE_MODES][N_CRUNS] = parallel::g_scratchMem->tResSurf;

      // iterate over elements
      for( TL_T_LID l_el = i_first; l_el < i_first+i_nElements; l_el++ ) {
        // anelastic updates (excluding frequency scaling)
        TL_T_REAL l_upA[TL_N_QTS_M][TL_N_MDS][TL_N_CRS];
        if( TL_N_RMS > 0) {
          for( unsigned short l_qt = 0; l_qt < TL_N_QTS_M; l_qt++ )
            for( unsigned short l_md = 0; l_md < TL_N_MDS; l_md++ )
              for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ )
                l_upA[l_qt][l_md][l_cr] = 0;
        }

        // add neighboring contribution
        for( TL_T_LID l_fa = 0; l_fa < TL_N_FAS; l_fa++ ) {
          TL_T_LID l_faId = i_elFa[l_el][l_fa];
          TL_T_LID l_ne;

          if( (i_faChars[l_faId].spType & OUTFLOW) != OUTFLOW ) {
            // derive neighbor
            if( (i_faChars[l_faId].spType & FREE_SURFACE) != FREE_SURFACE )
              l_ne = i_elFaEl[l_el][l_fa];
            else
              l_ne = l_el;

            /*
             * prefetches
             */
            const TL_T_REAL (* l_pre)[TL_N_MDS][TL_N_CRS] = nullptr;
            TL_T_LID l_neUp = std::numeric_limits<TL_T_LID>::max();
            // prefetch for the upcoming surface integration of this element
            if( l_fa < TL_N_FAS-1 ) l_neUp = i_elFaEl[l_el][l_fa+1];
            // first surface integration of the next element
            else if( l_el < i_first+i_nElements-1 ) l_neUp = i_elFaEl[l_el+1][0];

            // only proceed with adjacent data if the element exists
            if( l_neUp != std::numeric_limits<TL_T_LID>::max() ) l_pre = i_tDofsDg[0][l_neUp];
            // next element data in case of boundary conditions
            else if( l_el < i_first+i_nElements-1 )              l_pre = io_dofsE[l_el+1];
            // default to element data to avoid performance penality
            else                                                 l_pre = io_dofsE[l_el];

            /*
             * solve
             */
            // default are free-surface boundaries
            unsigned short l_vId = std::numeric_limits< unsigned short >::max();
            unsigned short l_fId = std::numeric_limits< unsigned short >::max();
            // switch to mesh-ids if not at the free surface
            if( (i_faChars[l_faId].spType & FREE_SURFACE) != FREE_SURFACE ) {
              l_vId = i_vIdElFaEl[l_el][l_fa];
              l_fId = i_fIdElFaEl[l_el][l_fa];
            }

            m_kernels->m_surfInt.neigh( l_fa,
                                        l_vId,
                                        l_fId,
                                        ( TL_T_REAL (*)[TL_N_QTS_E] )  ( i_fluxSolvers[l_el][l_fa].solver[0] ), // TODO: fix struct
                                        m_fsA[1][l_el][l_fa],
                                        i_tDofsDg[0][l_ne],
                                        io_dofsE[l_el],
                                        l_upA,
                                        l_tmpFa,
                                        l_pre );
          }
        }

        // update anelastic DOFs
        if( TL_N_RMS > 0 ) {
          TL_T_REAL (*l_dofsA)[TL_N_QTS_M][TL_N_MDS][TL_N_CRS] =
            (TL_T_REAL (*) [TL_N_QTS_M][TL_N_MDS][TL_N_CRS]) (io_dofsA+l_el*std::size_t(TL_N_RMS)*std::size_t(TL_N_QTS_M));

          m_kernels->m_surfInt.scatterUpdateA( l_upA, l_dofsA );
        }

        // compute extrema (if required)
        if( (i_elChars[l_el].spType & EXTREMA) != EXTREMA ) {}
        else {
          //! TODO: Use dedicated scratch memory for this
          TL_T_REAL (*l_sg)[TL_N_SCS][TL_N_CRS] = parallel::g_scratchMem->sg;

          // compute DG extrema
          edge::sc::Kernels< TL_T_EL,
                             TL_O_SP,
                             TL_N_QTS_E,
                             TL_N_CRS >::dgExtrema(  i_mm,
                                                     io_dofsE[l_el],
                                                     i_scatter,
                                                     l_sg,
                                                     o_extC[l_ex][0],
                                                     o_extC[l_ex][1] );

          // store the surface sub-cells
          if( ( i_elChars[l_el].spType & LIMIT_PLUS ) == LIMIT_PLUS ) {
            // set admissibility
            if( ( i_elChars[l_el].spType & LIMIT ) == LIMIT ) {
              if( ( i_elChars[l_el].spType & RUPTURE ) != RUPTURE ) {
                bool l_adm[TL_N_CRS];
                edge::sc::Detections< TL_T_EL,
                                      TL_N_QTS_E,
                                      TL_N_CRS >::dmpFa( i_extP[l_ex],
                                                         i_extP,
                                                         o_extC[l_ex],
                                                         i_lpFaLp[l_lp],
                                                         l_adm );

                // update admissibility
#pragma omp simd
                for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {
                  // only write "false" to memory, not "true" (avoids conflicts with rupture-admissibility in shared memory parallelization)
                  if( l_adm[l_cr] == false )
                    io_admC[l_li][l_cr] = false;
                }
              }

              // increase counter of limited elements
              l_li++;
            }
            // increase counter of limited plus elements
            l_lp++;
          }

          // increase counter of extrema elements
          l_ex++;
        }
      }
    }
};

#endif
