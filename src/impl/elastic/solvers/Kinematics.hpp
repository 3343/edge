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
 * Support for kinematic sources.
 **/
#ifndef KINEMATICS_SOLVER_HPP
#define KINEMATICS_SOLVER_HPP

#include "io/logging.h"
#include "Kinematics.type"
#include "linalg/Series.hpp"
#include "constants.hpp"

namespace edge {
  namespace elastic {
    namespace solvers {
      template< t_entityType   TL_T_EL,
                unsigned short TL_N_QTS,
                unsigned short TL_O_SP,
                unsigned short TL_N_FRUNS,
                unsigned short TL_N_FSRCS
              >
      class Kinematics;
    }
  }
}

/**
 * Solver for kinematic sources.
 *
 * We assume one, and only one kinematic source per run.
 * Since no restrictions on the corresponding points sources are enforced,
 * this is only an implementation-specific consideration rather than a limitation.
 *
 * Kinematics are able to operate in two modes.
 *
 *   1) TL_N_FRUNS equals TL_N_FSRCS
 *      All kinematic sources share the same location w.r.t. their corresponding points sources.
 *      In this case the number of fused sources equals the number of fused runs.
 *      This implies that the onset time of the points sources and the number of samples match.
 *      We can enforce this in preprocessing (if the locations match). However, if the time intervals
 *      are not very similar, case 2) might be a better match.
 *
 *   2) TL_N_FSRCS equals 1.
 *      The position of the corresponding point sources is not shared among the fused simulations.
 *
 * @paramt TL_T_EL element type.
 * @paramt TL_N_QTS number of quantities.
 * @paramt TL_O_SP order of the used solver.
 * @paramt TL_N_FRUNS number of fused forward runs.
 * @paramt TL_N_FSRCS number of fused sources.
 **/
template< t_entityType   TL_T_EL,
          unsigned short TL_N_QTS,
          unsigned short TL_O_SP,
          unsigned short TL_N_FRUNS,
          unsigned short TL_N_FSRCS
        >
class edge::elastic::solvers::Kinematics {
  private:
    //! number of dimensions
    static unsigned short const TL_N_DIM = C_ENT[TL_T_EL].N_DIM;
    //! number fo element modes
    static unsigned short const TL_N_MDS = CE_N_ELEMENT_MODES( TL_T_EL, TL_O_SP );
    //! number of independent, non-fused kinematic source descriptions
    static unsigned short const TL_N_IND_SRCS = (TL_N_FSRCS == 1) ? TL_N_FRUNS : 1;
    //! number of values in the stress tensors
    static unsigned short const TL_N_STRESS = (TL_N_DIM==2) ? 3 : 6;

  public:
    /**
     * Apply kinematic sources via dirac delta distribution.
     *
     *  The layout of i_elSpSo gives the first possible point source in a source-element.
     *  "possible" means that these are only source terms of the element, if
     *  the number of source in this element is > 0. The number of sources in an
     *  source-element is >0, if the next source-element's first source terms differs.
     *  The last entry is a ghost element, which holds the total number of sources.
     *
     *  Example with two independent source terms:
     *    el   0   1   2   3   4   5
     *    so |0 0|2 0|2 0|2 3|5 3|6 6|
     *
     *   * Two point sources of the first kinematic source description are in element 0 (first: 0).
     *   * No point sources are in element 1.
     *   * Three point sources of the second kinematic source description are in element 2 (first 0).
     *   * Three point sources of the first kinematic source description are in element 3 (first: 2).
     *   * One point source of the 1st and three point sources of the 2nd desc are in element 4 (first: 5/3).
     *
     *
     * @param i_first first sparse source-element to which sources are applied.
     * @param i_size number of sparse source elements.
     * @param i_t1 start time at which the sources are applied.
     * @param i_t2 end time untile which the sources are applied.
     * @param i_solver solvers for the kinematic source teroms.
     * @param i_elSpSo first id of the first possible souce term in a asparseelement.
     * @param io_dofs will be updated with source-contributions.
     *
     * @paramt TL_T_INT_LID intgral type of local ids.
     * @paramt TL_T_REAL floating point type.
     **/
    template< typename       TL_T_INT_LID,
              typename       TL_T_REAL >
    static void applyDirac( TL_T_INT_LID                    i_first,
                            TL_T_INT_LID                    i_size,
                            TL_T_REAL                       i_t1,
                            TL_T_REAL                       i_t2,
                            TL_T_INT_LID            const (*i_elSpSo)[TL_N_IND_SRCS],
                            solvers::t_Kinematics<
                              TL_N_DIM,
                              TL_N_MDS,
                              TL_N_FSRCS,
                              TL_T_REAL,
                              TL_T_INT_LID >        const   i_solvers[TL_N_IND_SRCS],
                            TL_T_REAL                     (*io_dofs)[TL_N_QTS][TL_N_MDS][TL_N_FRUNS] ) {
      // iterate over sparse source elements
      for( TL_T_INT_LID l_el = i_first; l_el < i_first+i_size; l_el++ ) {
        // iterate over kinematic sources
        for( TL_T_INT_LID l_ki = 0; l_ki < TL_N_IND_SRCS; l_ki++ ) {
          // iterate over sources in this element
          for( TL_T_INT_LID l_so = i_elSpSo[l_el][l_ki]; l_so < i_elSpSo[l_el+1][l_ki]; l_so++ ) {
            // iterate over slip-directions
            for( unsigned short l_sd = 0; l_sd < TL_N_DIM; l_sd++ ) {

              // only continue for active slip
              if( i_solvers[l_ki].aSlip[l_sd] ) {
                // compute slip
                TL_T_REAL l_slip[TL_N_FSRCS];

                linalg::Series< 
                  TL_N_FSRCS
                >::integrate( i_solvers[l_ki].dt[l_so],
                              i_solvers[l_ki].onSet[l_so],
                              i_solvers[l_ki].first[l_sd][l_so+1]-i_solvers[l_ki].first[l_sd][l_so],
                              i_solvers[l_ki].sr[l_sd]+i_solvers[l_ki].first[l_sd][l_so],
                              i_t1,
                              i_t2,
                              l_slip,
                              (TL_T_REAL) 0.0 );

                // iterate over stresses
                for( unsigned short l_qt = 0; l_qt < TL_N_STRESS; l_qt++ ) {
                  TL_T_REAL l_slipS[TL_N_FSRCS];
                  // scale slip with moment tensor coefficients
                  for( unsigned short l_fs = 0; l_fs < TL_N_FSRCS; l_fs++ ) {
                    l_slipS[l_fs] = l_slip[l_fs] * i_solvers[l_ki].sSca[l_sd][l_so][l_qt][l_fs];
                  }

                  // apply source contribution of the slip
                  for( unsigned short l_md = 0; l_md < TL_N_MDS; l_md++ ) {
                    for( unsigned short l_fs = 0; l_fs < TL_N_FSRCS; l_fs++ ) {
                      // derive id of fused run
                      unsigned short l_ru = l_ki * TL_N_FSRCS + l_fs;
                      // dense element id
                      TL_T_INT_LID l_elDe = i_solvers[l_ki].soElDe[l_so];

                      io_dofs[l_elDe][l_qt][l_md][l_ru] += l_slipS[l_fs] *
                                                           i_solvers[l_ki].bEval[l_so][l_md];
                    }
                  }
                }
              }

            }
          }
        }
      }
    }
};

#endif
