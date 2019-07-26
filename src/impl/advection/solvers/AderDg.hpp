/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2019, Alexander Breuer
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
 * ADER-DG solver for the advection equation.
 **/
#ifndef EDGE_ADVECTION_ADER_DG_HPP
#define EDGE_ADVECTION_ADER_DG_HPP

#include "constants.hpp"
#include "mesh/common.hpp"
#include "linalg/Matrix.h"
#include "linalg/Mappings.hpp"
#include "linalg/Matrix.h"
#include "../kernels/TimePred.hpp"
#include "../kernels/VolInt.hpp"
#include "../kernels/SurfInt.hpp"
#include "sc/Kernels.hpp"
#include "sc/Detections.hpp"

namespace edge {
  namespace advection {
    namespace solvers {
      template< typename       TL_T_REAL,
                t_entityType   TL_T_EL,
                unsigned short TL_O_SP,
                unsigned short TL_O_TI,
                unsigned short TL_N_CRS >
      class AderDg;
    }
  }
}

/**
 * ADER-DG solver for the advection equation split into local and neighboring updates.
 *
 * @paramt TL_T_REAL real type.
 * @paramt TL_T_EL element type.
 * @paramt TL_O_SP spatial order.
 * @paramt TL_O_TI temporal order.
 * @paramt TL_N_CRS number of fused simulations.
 **/
template< typename       TL_T_REAL,
          t_entityType   TL_T_EL,
          unsigned short TL_O_SP,
          unsigned short TL_O_TI,
          unsigned short TL_N_CRS >
class edge::advection::solvers::AderDg {
  private:
    //! number of dimensions
    static unsigned short const TL_N_DIS = C_ENT[TL_T_EL].N_DIM;

    //! number of vertices
    static unsigned short const TL_N_VES = C_ENT[TL_T_EL].N_VERTICES;

    //! number of faces
    static unsigned short const TL_N_FAS = C_ENT[TL_T_EL].N_FACES;

    //! number of quantities
    static unsigned short const TL_N_QTS = 1;

    //! number of DG modes
    static unsigned short const TL_N_MDS = CE_N_ELEMENT_MODES( TL_T_EL, TL_O_SP );

    //! time prediction kernel
    kernels::TimePred< TL_T_REAL,
                       TL_T_EL,
                       TL_O_SP,
                       TL_O_TI,
                       TL_N_CRS > * m_time;

    //! volume kernel
    kernels::VolInt< TL_T_REAL,
                     TL_T_EL,
                     TL_O_SP,
                     TL_N_CRS > * m_volInt;

    //! surface kernel
    kernels::SurfInt< TL_T_REAL,
                      TL_T_EL,
                      TL_O_SP,
                      TL_N_CRS > * m_surfInt;

  public:
    /**
     * Constructor.
     *
     * @param io_dynMem dynamic memory management, which will be used for the respective allocations.
     **/
    AderDg( data::Dynamic & io_dynMem ) {
      // init kernels
      m_time = new kernels::TimePred< TL_T_REAL,
                                      TL_T_EL,
                                      TL_O_SP,
                                      TL_O_TI,
                                      TL_N_CRS >( io_dynMem );

      m_volInt = new kernels::VolInt< TL_T_REAL,
                                      TL_T_EL,
                                      TL_O_SP,
                                      TL_N_CRS >( io_dynMem );

      m_surfInt = new kernels::SurfInt< TL_T_REAL,
                                        TL_T_EL,
                                        TL_O_SP,
                                        TL_N_CRS >( io_dynMem );
    }

    /**
     *  Destructor
     **/
    ~AderDg() {
      // free memory
      delete m_time;
      delete m_volInt;
      delete m_surfInt;
    }

    /**
     * Sets up the star matrices, which are a linear combination of the Jacobians.
     *
     * @param i_nEls number of elements.
     * @param i_vertexChars vertex characteristics.
     * @param i_elVe vertices adjacent to the elements.
     * @param i_bgPars background parameters.
     * @param o_starMatrices will be set to star matrices.
     **/
    static void setupStarM(       int_el           i_nEls,
                            const t_vertexChars   *i_vertexChars,
                            const int_el         (*i_elVe)[TL_N_VES],
                            const t_bgPars       (*i_bgPars)[1],
                                  TL_T_REAL      (*o_starMatrices)[TL_N_DIS] ) {
      // iterate over elements
      for( int_el l_el = 0; l_el < i_nEls; l_el++ ) {
        // derive vertex coords
        double l_veCoords[TL_N_DIS][TL_N_VES];
        mesh::common< TL_T_EL >::getElVeCrds( l_el, i_elVe, i_vertexChars, l_veCoords );

        // get inverse jacobian
        double l_jac[TL_N_DIS][TL_N_DIS];
        linalg::Mappings::evalJac( TL_T_EL, l_veCoords[0], l_jac[0] );

        double l_jacInv[TL_N_DIS][TL_N_DIS];

// TODO: Fix preprocessor usage
#if PP_N_DIM == 1
        l_jacInv[0][0] = 1 / l_jac[0][0];
#else
        linalg::Matrix::inv( l_jac, l_jacInv );
#endif

        // set star matrices
        // iterate over reference dimensions
        for( unsigned int l_dim = 0; l_dim < TL_N_DIS; l_dim++ ) {
          o_starMatrices[l_el][l_dim]  = i_bgPars[l_el][0].a * l_jacInv[0][l_dim];

// TODO: Fix preprocessor usage
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
     * Local step: Time predictions + volume + local surface.
     *
     * @param i_first first element considered.
     * @param i_nEls number of elements.
     * @param i_dt time step.
     * @param i_elChars element characteristics.
     * @param i_starM star matrices.
     * @param i_fluxSolvers flux solvers.
     * @param io_dofsDg DG DOFs which will be updates with the elements' local contributions.
     * @param o_tDofsDg temporary DG DOFs (buffers or DOFs and buffers), which will be updated.
     *
     * @paramt TL_T_LID integral type of local ids.
     * @paramt TL_T_CHARS_EL element characteristics, offering sparse type .spType.
     **/
    template< typename TL_T_LID,
              typename TL_T_CHARS_EL >
    void local( TL_T_LID                                i_first,
                TL_T_LID                                i_nEls,
                double                                  i_dt,
                TL_T_CHARS_EL  const                   *i_elChars,
                TL_T_REAL      const                  (*i_starM)[TL_N_DIS],
                TL_T_REAL      const                  (*i_fluxSolvers)[TL_N_FAS*2],
                TL_T_REAL                             (*io_dofsDg)[TL_N_QTS][TL_N_MDS][TL_N_CRS],
                TL_T_REAL            (* const * const   o_tDofsDg[2])[TL_N_MDS][TL_N_CRS] ) {
      // iterate over all elements
      for( TL_T_LID l_el = i_first; l_el < i_first+i_nEls; l_el++ ) {
        // compute ader time prediction
        TL_T_REAL l_derBuffer[TL_O_TI][TL_N_MDS][TL_N_CRS];
        m_time->ck( i_dt,
                    i_starM[l_el],
                    io_dofsDg[l_el][0],
                    l_derBuffer,
                    o_tDofsDg[0][l_el][0] );

        // compute volume contribution
        m_volInt->apply( i_starM[l_el],
                         o_tDofsDg[0][l_el],
                         io_dofsDg[l_el] );

        // compute local surface contribution
        m_surfInt->local( i_fluxSolvers[l_el],
                          o_tDofsDg[0][l_el],
                          io_dofsDg[l_el] );
      }
    }

    /**
     * Performs the neighboring updates of the ADER-DG scheme.
     *
     * @param i_first first element considered.
     * @param i_nEls number of elements.
     * @param i_faChars face characteristics.
     * @param i_elChars element characteristics.
     * @param i_fluxSolvers flux solvers.
     * @param i_elFa elements' adjacent faces (no bridge).
     * @param i_elFaEl elements' adjacent elements (faces as bridge).
     * @param i_fIdElFaEl local face ids of face-neighboring elements.
     * @param i_vIdElFaEl local vertex ids w.r.t. the shared face from the neighboring elements' perspective.
     * @param i_tDofs temporary DOFs 1) buffers or 2) DOFs and buffers, which will be used for the elements' updates.
     * @param io_dofs DOFs which will be updated with neighboring elements' contribution.
     *
     * @paramt TL_T_LID integral type of local ids.
     * @paramt TL_T_CHARS_FA face characteristics, offering sparse type .spType.
     * @paramt TL_T_CHARS_EL element characteristics, offering sparse type .spType.
     **/
    template< typename TL_T_LID,
              typename TL_T_CHARS_FA,
              typename TL_T_CHARS_EL >
    void neigh( TL_T_LID                                i_first,
                TL_T_LID                                i_nEls,
                TL_T_CHARS_FA  const                   *i_faChars,
                TL_T_CHARS_EL  const                   *i_elChars,
                TL_T_REAL      const                  (*i_fluxSolvers)[TL_N_FAS*2],
                TL_T_LID       const                  (*i_elFa)[TL_N_FAS],
                TL_T_LID       const                  (*i_elFaEl)[TL_N_FAS],
                unsigned short const                  (*i_fIdElFaEl)[TL_N_FAS],
                unsigned short const                  (*i_vIdElFaEl)[TL_N_FAS],
                TL_T_REAL            (* const * const   i_tDofs[2])[TL_N_MDS][TL_N_CRS],
                TL_T_REAL                             (*io_dofs)[TL_N_QTS][TL_N_MDS][TL_N_CRS] ) {
      // iterate over elements
      for( TL_T_LID l_el = i_first; l_el < i_first+i_nEls; l_el++ ) {
        for( unsigned short l_fa = 0; l_fa < TL_N_FAS; l_fa++ ) {
          TL_T_LID l_faId = i_elFa[l_el][l_fa];
          TL_T_LID l_ne;

          if( (i_faChars[l_faId].spType & OUTFLOW) != OUTFLOW ) {
            // derive neighbor
            l_ne = i_elFaEl[l_el][l_fa];
          }
          // outflow boundary conditions
          else {
            l_ne = l_el;
          }

          // default are outflow boundaries
          unsigned short l_vId = std::numeric_limits< unsigned short >::max();
          unsigned short l_fId = std::numeric_limits< unsigned short >::max();
          // switch to mesh-ids if not outflow
          if( (i_faChars[l_faId].spType & OUTFLOW) != OUTFLOW ) {
            l_vId = i_vIdElFaEl[l_el][l_fa];
            l_fId = i_fIdElFaEl[l_el][l_fa];
          }

          m_surfInt->neigh( l_fa,
                            l_vId,
                            l_fId,
                            i_fluxSolvers[l_el][TL_N_FAS+l_fa],
                            i_tDofs[0][l_ne],
                            io_dofs[l_el] );
      }
    }
  }
};

#endif
