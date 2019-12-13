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

    //! number of DG face modes
    static unsigned short const TL_N_MDS_FA = CE_N_ELEMENT_MODES( C_ENT[TL_T_EL].TYPE_FACES, TL_O_SP );

    //! number of DG element modes
    static unsigned short const TL_N_MDS_EL = CE_N_ELEMENT_MODES( TL_T_EL, TL_O_SP );

    //! number of elastic quantities
    static unsigned short const TL_N_QTS_E = CE_N_QTS_E( TL_N_DIS );

    //! number of quantities per relaxation mechanism
    static unsigned short const TL_N_QTS_M = CE_N_QTS_M( TL_N_DIS );

    //! number of sub-faces per DG-face
    static unsigned short const TL_N_SFS = CE_N_SUB_FACES( TL_T_EL, TL_O_SP );

    //! number of subcells per DG-element
    static unsigned short const TL_N_SCS = CE_N_SUB_CELLS( TL_T_EL, TL_O_SP );

    //! elastic star matrices
    static unsigned short const TL_N_ENS_STAR_E = (TL_MATS_SP) ? CE_N_ENS_STAR_E_SP( TL_N_DIS )
                                                               : CE_N_ENS_STAR_E_DE( TL_N_DIS );
    TL_T_REAL (*m_starE)[TL_N_DIS][TL_N_ENS_STAR_E] = nullptr;

    //! elastic flux solvers
    static unsigned short const TL_N_ENS_FS_E = CE_N_ENS_FS_E_DE( TL_N_DIS );
    TL_T_REAL (*m_fsE[2])[TL_N_FAS][TL_N_ENS_FS_E] = { nullptr, nullptr };

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

      // elastic star matrices
      l_size = i_nEls * std::size_t(TL_N_DIS) * TL_N_ENS_STAR_E;
      l_size *= sizeof(TL_T_REAL);

      m_starE = ( TL_T_REAL (*) [TL_N_DIS][TL_N_ENS_STAR_E] ) io_dynMem.allocate( l_size,
                                                                                  i_align,
                                                                                  false,
                                                                                  true );

      // elastic flux solvers
      l_size = i_nEls * std::size_t(TL_N_FAS) * TL_N_ENS_FS_E;
      l_size *= sizeof(TL_T_REAL);

      for( unsigned short l_sd = 0; l_sd < 2; l_sd++ ) {
        m_fsE[l_sd] = ( TL_T_REAL (*) [TL_N_FAS][TL_N_ENS_FS_E] ) io_dynMem.allocate( l_size,
                                                                                      i_align,
                                                                                      false,
                                                                                      true );
      }

      // anelastic part
      if( TL_N_RMS > 0 ) {
        // anelastic source matrices
        l_size = i_nEls * std::size_t(TL_N_RMS) * TL_N_ENS_SRC_A;
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
     * @param i_nElsIn number of inner elements.
     * @param i_nElsSe number of send elements.
     * @param i_nFas number of faces.
     * @param i_nCommElFa number of communicating element-face pairs.
     * @param i_recvFa local ids of the receiving face within the elements.
     * @param i_recvEl ids of the receiving elements.
     * @param i_faEl elements adjacent to faces.
     * @param i_elVe vertices adjacent to elements.
     * @param i_elFa elements adjacent to elements through faces as bridge.
     * @param i_veChars vertex characteristics.
     * @param i_faChars face characteristics.
     * @param i_elChars elements characteristics.
     * @param io_bgPars background parameters (phase lambda and mu will be replaced with elastic lambda an mu if TL_N_RMS>0).
     * @param i_bgParsRe background parameters of the receive elements.
     * @param i_freqCen central frequency for attenuation.
     * @param i_freqRat frequency ratio between upper and lower frequencies for attenuation.
     * @param io_dynMem dynamic memory management.
     *
     * @paramt TL_T_LID integral type of local ids.
     */
    template< typename TL_T_LID >
    AderDg( TL_T_LID                i_nElsIn,
            TL_T_LID                i_nElsSe,
            TL_T_LID                i_nFas,
            TL_T_LID                i_nCommElFa,
            unsigned short const  * i_recvFa,
            TL_T_LID       const  * i_recvEl,
            TL_T_LID       const (* i_faEl)[2],
            TL_T_LID       const (* i_elVe)[TL_N_VES_EL],
            TL_T_LID       const (* i_elFa)[TL_N_FAS],
            t_vertexChars  const  * i_veChars,
            t_faceChars    const  * i_faChars,
            t_elementChars const  * i_elChars,
            t_bgPars              * io_bgPars,
            t_bgPars       const  * i_bgParsRe,
            double                  i_freqCen,
            double                  i_freqRat,
            data::Dynamic         & io_dynMem ) {
      // total number of elements
      std::size_t l_nEls = i_nElsIn + i_nElsSe;

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
      alloc( l_nEls,
             ALIGNMENT.BASE.HEAP,
             io_dynMem );

      // init anelastic source matrices and compute elastic Lame parameters in viscoelastic settings
      if( TL_N_RMS > 0 ) {
        AderDgInit< TL_T_EL,
                    TL_MATS_SP >::initSrcA( l_nEls,
                                            TL_N_RMS,
                                            i_freqCen,
                                            i_freqRat,
                                            io_bgPars,
                                            m_srcA );
      }

      // init star matrices
      AderDgInit< TL_T_EL,
                  TL_MATS_SP >::initStar( l_nEls,
                                          i_elVe,
                                          i_veChars,
                                          io_bgPars,
                                          m_starE,
                                          m_starA );

      // init flux solvers
      AderDgInit< TL_T_EL,
                  TL_MATS_SP >::initFs( i_nElsIn,
                                        i_nElsSe,
                                        i_nFas,
                                        i_nCommElFa,
                                        i_recvFa,
                                        i_recvEl,
                                        i_faEl,
                                        i_elVe,
                                        i_elFa,
                                        i_veChars,
                                        i_faChars,
                                        i_elChars,
                                        io_bgPars,
                                        i_bgParsRe,
                                        m_fsE,
                                        m_fsA );
    }

    /**
     * Destructor.
     **/
    ~AderDg() {
      delete m_kernels;
    }

    /**
     * Local step: ADER + volume + local surface.
     *
     * @param i_first first element considered.
     * @param i_nEls number of elements.
     * @param i_firstTs true if this is the first time step of every rate-2 ts pair.
     * @param i_time time of the initial DOFs.
     * @param i_dt time step.
     * @param i_firstSpRe first sparse receiver entity.
     * @param i_elChars element characteristics.
     * @param io_dofsE elastic DOFs.
     * @param io_dofsA anelastic DOFs.
     * @param o_tDofs time integrated DG DOFs (==, <, >) which will be set or updated.
     * @param io_recvs will be updated with receiver info.
     *
     * @paramt TL_T_LID integer type of local entity ids.
     **/
    template < typename TL_T_LID >
    void local( TL_T_LID                             i_first,
                TL_T_LID                             i_nEls,
                bool                                 i_firstTs,
                double                               i_time,
                double                               i_dt,
                TL_T_LID                             i_firstSpRe,
                t_elementChars              const  * i_elChars,
                unsigned short const              (* i_vIdElFaEl)[TL_N_FAS],
                TL_T_REAL                         (* io_dofsE)[TL_N_QTS_E][TL_N_MDS_EL][TL_N_CRS],
                TL_T_REAL                         (* io_dofsA)[TL_N_MDS_EL][TL_N_CRS],
                TL_T_REAL        (* const * const    o_tDofs[3])[TL_N_MDS_EL][TL_N_CRS],
                TL_T_REAL        (* const * const    o_sendDofs)[TL_N_MDS_FA][TL_N_CRS],
                edge::io::Receivers                & io_recvs ) const {
      // counter for receivers
      unsigned int l_enRe = i_firstSpRe;

      // temporary data structurre for product for two-way mult and receivers
      TL_T_REAL (*l_tmp)[TL_N_MDS_EL][TL_N_CRS] = parallel::g_scratchMem->tRes[0];

      // buffer for derivatives
      TL_T_REAL (*l_derBuffer)[TL_N_QTS_E][TL_N_MDS_EL][TL_N_CRS] = parallel::g_scratchMem->dBuf;

      // iterate over all elements
      for( TL_T_LID l_el = i_first; l_el < i_first+i_nEls; l_el++ ) {
        // pointer to anelastic dofs
        TL_T_REAL (*l_dofsA)[TL_N_QTS_M][TL_N_MDS_EL][TL_N_CRS] =
          (TL_T_REAL (*) [TL_N_QTS_M][TL_N_MDS_EL][TL_N_CRS]) (io_dofsA+l_el*std::size_t(TL_N_RMS)*std::size_t(TL_N_QTS_M));

        TL_T_REAL l_tDofsA[CE_MAX(int(TL_N_RMS),1)][TL_N_QTS_M][TL_N_MDS_EL][TL_N_CRS];
        TL_T_REAL l_derA[CE_MAX(int(TL_N_RMS),1)][TL_O_SP][TL_N_QTS_M][TL_N_MDS_EL][TL_N_CRS];

        // compute ADER time integration
        m_kernels->m_time.ck( i_dt,
                              m_starE[l_el],
                              (TL_N_RMS > 0) ? m_starA[l_el] : nullptr,
                              m_srcA+l_el*std::size_t(TL_N_RMS),
                              io_dofsE[l_el],
                              l_dofsA,
                              l_tmp,
                              l_derBuffer,
                              l_derA,
                              o_tDofs[0][l_el],
                              l_tDofsA );

        // update summed time integrated elastic DOFs, if an adjacent element has a larger time step
        if( (i_elChars[l_el].spType & C_LTS_EL[EL_INT_LT]) == C_LTS_EL[EL_INT_LT] ) {
          // reset, if required
          if( i_firstTs ) {
            for( unsigned short l_qt = 0; l_qt < TL_N_QTS_E; l_qt++ )
              for( unsigned short l_md = 0; l_md < TL_N_MDS_EL; l_md++ )
                for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ )
                  o_tDofs[1][l_el][l_qt][l_md][l_cr] = 0;
          }

          // add tDofs of this time step
          for( unsigned short l_qt = 0; l_qt < TL_N_QTS_E; l_qt++ )
            for( unsigned short l_md = 0; l_md < TL_N_MDS_EL; l_md++ )
              for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ )
                o_tDofs[1][l_el][l_qt][l_md][l_cr] += o_tDofs[0][l_el][l_qt][l_md][l_cr];
        }

        // compute [0, 0.5dt] time integrated DOFs, if an adjacent element has a smaller time step
        if( (i_elChars[l_el].spType & C_LTS_EL[EL_INT_GT]) == C_LTS_EL[EL_INT_GT] ) {
          m_kernels->m_time.integrate( TL_T_REAL(0.5*i_dt),
                                       l_derBuffer,
                                       o_tDofs[2][l_el] );
        }

#ifdef PP_USE_MPI
        for( unsigned short l_fa = 0; l_fa < TL_N_FAS; l_fa++ ) {
          unsigned short l_vId = i_vIdElFaEl[l_el][l_fa];

          if( o_sendDofs[l_el*TL_N_FAS + l_fa] != nullptr ) {
            // gts
            if( (i_elChars[l_el].spType & C_LTS_AD[l_fa][AD_EQ]) == C_LTS_AD[l_fa][AD_EQ] ) {
              m_kernels->m_surfInt.neighFluxInt( std::numeric_limits< unsigned short >::max(),
                                                 l_vId,
                                                 l_fa,
                                                 o_tDofs[0][l_el],
                                                 o_sendDofs[l_el*TL_N_FAS + l_fa] );
            }
            // less than
            else if( (i_elChars[l_el].spType & C_LTS_AD[l_fa][AD_LT]) == C_LTS_AD[l_fa][AD_LT] ) {
              if( !i_firstTs ) {
                m_kernels->m_surfInt.neighFluxInt( std::numeric_limits< unsigned short >::max(),
                                                   l_vId,
                                                   l_fa,
                                                   o_tDofs[1][l_el],
                                                   o_sendDofs[l_el*TL_N_FAS + l_fa] );
              }
            }
            // greater than
            else {
              m_kernels->m_surfInt.neighFluxInt( std::numeric_limits< unsigned short >::max(),
                                                 l_vId,
                                                 l_fa,
                                                 o_tDofs[2][l_el],
                                                 o_sendDofs[l_el*TL_N_FAS + l_fa] );

              m_kernels->m_surfInt.neighFluxInt( std::numeric_limits< unsigned short >::max(),
                                                 l_vId,
                                                 l_fa,
                                                 o_tDofs[0][l_el],
                                                 o_sendDofs[l_el*TL_N_FAS + l_fa]+TL_N_QTS_E );

              // move second integral from [0, dt] to [1/2dt, dt]
              for( unsigned short l_qt = 0; l_qt < TL_N_QTS_E; l_qt++ ) {
                for( unsigned short l_md = 0; l_md < TL_N_MDS_FA; l_md++ ) {
                  for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {
                    o_sendDofs[l_el*TL_N_FAS + l_fa][TL_N_QTS_E+l_qt][l_md][l_cr] -= o_sendDofs[l_el*TL_N_FAS + l_fa][l_qt][l_md][l_cr];
                  }
                }
              }
            }
          }
        }
#endif

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
        m_kernels->m_volInt.apply( m_starE[l_el],
                                   (TL_N_RMS > 0) ? m_starA[l_el] : nullptr,
                                   m_srcA+l_el*std::size_t(TL_N_RMS),
                                   o_tDofs[0][l_el],
                                   l_tDofsA,
                                   io_dofsE[l_el],
                                   l_dofsA,
                                   l_tmp );

        // prefetches for next iteration
        const TL_T_REAL (* l_preDofs)[TL_N_MDS_EL][TL_N_CRS] = nullptr;
        const TL_T_REAL (* l_preTint)[TL_N_MDS_EL][TL_N_CRS] = nullptr;

        if( l_el < i_first+i_nEls-1 ) {
          l_preDofs = io_dofsE[l_el+1];
          l_preTint = o_tDofs[0][l_el+1];
        }
        else {
          l_preDofs = io_dofsE[l_el];
          l_preTint = o_tDofs[0][l_el];
        }

         /*
          * compute local surface contribution
          */
        // reuse derivative buffer
        TL_T_REAL (*l_tmpFa)[N_QUANTITIES][N_FACE_MODES][N_CRUNS] = parallel::g_scratchMem->tResSurf;
        // call kernel
        m_kernels->m_surfInt.local( m_fsE[0][l_el],
                                    (TL_N_RMS > 0) ? m_fsA[0][l_el] : nullptr,
                                    o_tDofs[0][l_el],
                                    io_dofsE[l_el],
                                    l_dofsA,
                                    l_tmpFa,
                                    l_preDofs,
                                    l_preTint );
      }
    }

    /**
     * Performs the neighboring updates of the ADER-DG scheme.
     *
     * @param i_first first element considered.
     * @param i_nEls number of elements.
     * @param i_firstTs true if this is the first time step of every rate-2 ts pair.
     * @param i_faChars face characteristics.
     * @param i_elChars element characteristics.
     * @param i_elFa elements' adjacent faces.
     * @param i_elFaEl face-neighboring elements.
     * @param i_fIdElFaEl local face ids of face-neighboring elememts.
     * @param i_vIdElFaEl local vertex ids w.r.t. the shared face from the neighboring elements' perspsective.
     * @param i_tDofs time integrated DG DOFs (==, <, >).
     * @param io_dofs DOFs which will be updated with neighboring elements' contribution.
     *
     * @paramt TL_T_LID integer type of local entity ids.
     **/
    template< typename TL_T_LID >
    void neigh( TL_T_LID                              i_first,
                TL_T_LID                              i_nEls,
                bool                                  i_firstTs,
                t_faceChars    const                * i_faChars,
                t_elementChars const                * i_elChars,
                TL_T_LID       const               (* i_elFa)[TL_N_FAS],
                TL_T_LID       const               (* i_elFaEl)[TL_N_FAS],
                unsigned short const               (* i_fIdElFaEl)[TL_N_FAS],
                unsigned short const               (* i_vIdElFaEl)[TL_N_FAS],
                TL_T_REAL            (* const * const i_tDofs[3])[TL_N_MDS_EL][TL_N_CRS],
                TL_T_REAL                          (* io_dofsE)[TL_N_QTS_E][TL_N_MDS_EL][TL_N_CRS],
                TL_T_REAL                          (* io_dofsA)[TL_N_MDS_EL][TL_N_CRS],
                TL_T_REAL      const (* const * const i_recvDofs)[TL_N_MDS_FA][TL_N_CRS] ) const {
      // temporary product for three-way mult
      TL_T_REAL (*l_tmpFa)[N_QUANTITIES][N_FACE_MODES][N_CRUNS] = parallel::g_scratchMem->tResSurf;

      // iterate over elements
      for( TL_T_LID l_el = i_first; l_el < i_first+i_nEls; l_el++ ) {
        // anelastic updates (excluding frequency scaling)
        TL_T_REAL l_upA[TL_N_QTS_M][TL_N_MDS_EL][TL_N_CRS];
        if( TL_N_RMS > 0) {
          for( unsigned short l_qt = 0; l_qt < TL_N_QTS_M; l_qt++ )
            for( unsigned short l_md = 0; l_md < TL_N_MDS_EL; l_md++ )
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
            const TL_T_REAL (* l_pre)[TL_N_MDS_EL][TL_N_CRS] = nullptr;
            TL_T_LID l_neUp = std::numeric_limits<TL_T_LID>::max();
            // prefetch for the upcoming surface integration of this element
            if( l_fa < TL_N_FAS-1 ) l_neUp = i_elFaEl[l_el][l_fa+1];
            // first surface integration of the next element
            else if( l_el < i_first+i_nEls-1 ) l_neUp = i_elFaEl[l_el+1][0];

            // only proceed with adjacent data if the element exists
            if( l_neUp != std::numeric_limits<TL_T_LID>::max() ) l_pre = i_tDofs[0][l_neUp];
            // next element data in case of boundary conditions
            else if( l_el < i_first+i_nEls-1 )              l_pre = io_dofsE[l_el+1];
            // default to element data to avoid performance penality
            else                                                 l_pre = io_dofsE[l_el];

            // assemble the neighboring time integrated DOFs
            TL_T_REAL l_tDofs[TL_N_QTS_E][TL_N_MDS_EL][TL_N_CRS];
            TL_T_REAL const (*l_tDofsFiE)[TL_N_MDS_FA][TL_N_CRS] = nullptr;

            if( i_recvDofs == nullptr || i_recvDofs[l_el*TL_N_FAS + l_fa] == nullptr ) {
              for( unsigned short l_qt = 0; l_qt < TL_N_QTS_E; l_qt++ ) {
                for( unsigned short l_md = 0; l_md < TL_N_MDS_EL; l_md++ ) {
                  for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {
                    // element and face-adjacent one have an equal time step
                    if( (i_elChars[l_el].spType & C_LTS_AD[l_fa][AD_EQ]) == C_LTS_AD[l_fa][AD_EQ] ) {
                      l_tDofs[l_qt][l_md][l_cr] = i_tDofs[0][l_ne][l_qt][l_md][l_cr];
                    }
                    // element has a greater time step than the face-adjacent one
                    else if( (i_elChars[l_el].spType & C_LTS_AD[l_fa][AD_GT]) == C_LTS_AD[l_fa][AD_GT] ) {
                      l_tDofs[l_qt][l_md][l_cr] = i_tDofs[1][l_ne][l_qt][l_md][l_cr];
                    }
                    // element has a time step less than the adjacent one
                    else {
                      if( i_firstTs )
                        l_tDofs[l_qt][l_md][l_cr] = i_tDofs[2][l_ne][l_qt][l_md][l_cr];
                      else
                        l_tDofs[l_qt][l_md][l_cr] = i_tDofs[0][l_ne][l_qt][l_md][l_cr] - i_tDofs[2][l_ne][l_qt][l_md][l_cr];
                    }
                  }
                }
              }
            }
#ifdef PP_USE_MPI
            else {
              // derive offset
              std::size_t l_off = 0;
              if( (i_elChars[l_el].spType & C_LTS_AD[l_fa][AD_LT]) == C_LTS_AD[l_fa][AD_LT] ) {
                l_off = (i_firstTs) ? 0 : TL_N_QTS_E;
              }

              l_tDofsFiE = i_recvDofs[l_el*TL_N_FAS + l_fa]+l_off;
            }
#else
            else EDGE_LOG_FATAL;
#endif
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
                                        m_fsE[1][l_el][l_fa],
                                        (TL_N_RMS > 0) ? m_fsA[1][l_el][l_fa] : nullptr,
                                        l_tDofs,
                                        l_tDofsFiE,
                                        io_dofsE[l_el],
                                        l_upA,
                                        l_tmpFa,
                                        l_pre );
          }
        }

        // update anelastic DOFs
        if( TL_N_RMS > 0 ) {
          TL_T_REAL (*l_dofsA)[TL_N_QTS_M][TL_N_MDS_EL][TL_N_CRS] =
            (TL_T_REAL (*) [TL_N_QTS_M][TL_N_MDS_EL][TL_N_CRS]) (io_dofsA+l_el*std::size_t(TL_N_RMS)*std::size_t(TL_N_QTS_M));

          m_kernels->m_surfInt.scatterUpdateA( l_upA, l_dofsA );
        }
      }
    }
};

#endif
