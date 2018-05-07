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
 * Initialization of kinematic sources.
 **/

#ifndef KINEMATICS_SETUP_HPP
#define KINEMATICS_SETUP_HPP

#include "data/Dynamic.h"
#include "../solvers/Kinematics.type"
#include "constants.hpp"
#include "data/layout.hpp"
#include "mesh/common.hpp"
#include "dg/Basis.h"
#include "linalg/Mappings.hpp"
#include "linalg/Matrix.h"
#include <algorithm>

namespace edge {
  namespace elastic {
    namespace setups {
      template< typename       TL_T_INT_LID,
                typename       TL_T_REAL_MESH,
                typename       TL_T_REAL_COMP,
                t_entityType   TL_T_EL,
                unsigned short TL_O_SP,
                unsigned short TL_N_FRUNS,
                unsigned short TL_N_FSRCS
              >
      class KinematicsInit;
    }
  }
}

/**
 * Initialization of kinematics sources.
 *
 * We assume one, and only one kinematic source per run.
 * Since no restrictions on the corresponding points sources are enforced,
 * this is only an implementation-specific consideration rather than a limitation.
 *
 * Kinematics are able to operate in two modes.
 *
 *   1) TL_N_FSRCS equals 1.
 *      Each of the fused runs has a completely independent, non-fused kinematic source description.
 *
 *   2) TL_N_FRUNS equals TL_N_FSRCS
 *      All kinematic sources are identical except for the slip-rates and the coefficients of the slip rates (encoding slip-direction, subfault area, etc.).
 *
 * @paramt TL_T_LID integral type of local ids.
 * @paramt TL_T_REAL_MESH type for floating point values for mesh-related data.
 * @paramt TL_T_REAL_COMP type for floating point values for computational data.
 * @paramt TL_T_EL element type.
 * @paramt TL_O_SP order in space.
 * @paramt TL_N_FRUNS number of fused forward runs.
 * @paramt TL_N_FSRCS number of fused sources.
 **/
template< typename       TL_T_INT_LID,
          typename       TL_T_REAL_MESH,
          typename       TL_T_REAL_COMP,
          t_entityType   TL_T_EL,
          unsigned short TL_O_SP,
          unsigned short TL_N_FRUNS,
          unsigned short TL_N_FSRCS
        >
class edge::elastic::setups::KinematicsInit {
  private:
    // require that the number of fused srcs either matches the number of fused runs or is 1
    static_assert( TL_N_FRUNS == TL_N_FSRCS || TL_N_FSRCS == 1,
                   "#fused sources does not match source requirements" );

    //! number of dimensions.
    static unsigned short const TL_N_DIM      = C_ENT[TL_T_EL].N_DIM;
    //! number of values in the stress tensors
    static unsigned short const TL_N_STRESS = (TL_N_DIM==2) ? 3 : 6;
    //! number fo element modes.
    static unsigned short const TL_N_EL_MODES = CE_N_ELEMENT_MODES( TL_T_EL, TL_O_SP );
    //! number of independent, non-fused kinematic source descriptions
    static unsigned short const TL_N_IND_SRCS = (TL_N_FSRCS == 1) ? TL_N_FRUNS : 1;
    //! number of vertices
    static unsigned short const TL_N_VE       = C_ENT[TL_T_EL].N_VERTICES;

    /**
     * Determine the point sources local to the given domain.
     *
     * @param i_nSrcs number of point sources for the independent kinematic source configurations.
     * @param i_srcCrds coordinates of the point sources [*][][]: src config [][*][]: point source [][][*]: dim.
     * @param i_elLayout element layout.
     * @param i_elVe vertices adjacent to the elements.
     * @param i_veChars vertex characteristics.
     * @param i_gId global ids of the elements.
     * @param o_srcEl will be set to ids of elements having a point src. If an element has more than one point source, it will appear multiple times. [*][] src config, [][*]: local element ids.
     * @param o_srcId will be set to source ids of the active source local to the given domain.
     *
     * @paramt TL_T_INT_GID integer type of global ids.
     **/
    template< typename TL_T_INT_GID >
    static void getLocal( TL_T_INT_LID                const           i_nSrcs[TL_N_IND_SRCS],
#ifdef __INTEL_COMPILER
                          double                            (*        i_srcCrds[TL_N_IND_SRCS])[TL_N_DIM],
#else
                          double                      const (* const  i_srcCrds[TL_N_IND_SRCS])[TL_N_DIM],
#endif
                          t_enLayout                  const          &i_elLayout,
                          TL_T_INT_LID                const         (*i_elVe)[TL_N_VE],
                          t_vertexChars               const          *i_veChars,
                          TL_T_INT_GID                const          *i_gIds,
                          std::vector< TL_T_INT_LID >                 o_srcEl[TL_N_IND_SRCS],
                          std::vector< TL_T_INT_LID >                 o_srcId[TL_N_IND_SRCS] ) {
      PP_INSTR_FUN("get_local")

      // iterate over independent source configurations
      for( unsigned short l_is = 0; l_is < TL_N_IND_SRCS; l_is++ ) {
        // determine elements with minimum global ids holding the sources
        std::vector< TL_T_INT_LID > l_lIds( i_nSrcs[l_is] );
        std::vector< TL_T_INT_GID > l_gIds( i_nSrcs[l_is] );
        std::vector< int_tg       > l_tgs(  i_nSrcs[l_is] );
#ifdef PP_USE_OMP
#pragma omp parallel for
#endif
        for( TL_T_INT_LID l_so = 0; l_so < i_nSrcs[l_is]; l_so++ ) {
          // get source coordinates
          double l_crds[TL_N_DIM];
          for( unsigned short l_di = 0; l_di < TL_N_DIM; l_di++ ) l_crds[l_di] = i_srcCrds[l_is][l_so][l_di];

          // get the info of the miminum id element holding the source
          bool l_owned = mesh::common< TL_T_EL >::findMinGid( l_crds,
                                                              i_elLayout,
                                                              i_elVe[0],
                                                              i_veChars,
                                                              i_gIds,
                                                              l_lIds[l_so], l_gIds[l_so], l_tgs[l_so] );
          // set invalid if not owned
          if( l_owned ){}
          else {
            l_lIds[l_so] = std::numeric_limits< TL_T_INT_LID >::max();
            l_gIds[l_so] = std::numeric_limits< TL_T_INT_GID >::max();
            l_tgs[ l_so] = std::numeric_limits< int_tg       >::max();
          }
        }

        // reset size
        o_srcEl[l_is].resize(0); o_srcId[l_is].resize(0);
        // reserve memory for sources
        o_srcEl[l_is].reserve(i_nSrcs[l_is]); o_srcId[l_is].reserve(i_nSrcs[l_is]);

        // setup the source layout
        for( TL_T_INT_LID l_so = 0; l_so < i_nSrcs[l_is]; l_so++ ) {
          // get local id, global id and time group of the source
          TL_T_INT_LID l_lId = l_lIds[l_so];
          TL_T_INT_GID l_gId = l_gIds[l_so];
          int_tg       l_tg  = l_tgs[l_so];

          // only considered owned sources
          if( l_lId == std::numeric_limits< TL_T_INT_LID >::max() ) continue;

          // check if this is an inner src
          TL_T_INT_LID l_inFirst = i_elLayout.timeGroups[l_tg].inner.first;
          TL_T_INT_LID l_inSize  = i_elLayout.timeGroups[l_tg].inner.size;

          // check if this an inner-src
          bool l_innerSrc = l_lId < l_inFirst + l_inSize;

          // inner-srcs are only appliad once
          if( l_innerSrc ) {
            o_srcEl[l_is].push_back( l_lId );
            o_srcId[l_is].push_back( l_so );
          }
          // send-srcs are applied to every send-element
          else {
            // iterate over neighboring ranks
            for( unsigned int l_nr = 0; l_nr < i_elLayout.timeGroups[l_tg].send.size(); l_nr++ ) {
              TL_T_INT_LID l_sendFirst = i_elLayout.timeGroups[l_tg].send[l_nr].first;
              TL_T_INT_LID l_sendSize  = i_elLayout.timeGroups[l_tg].send[l_nr].size;

              for( TL_T_INT_LID l_en = l_sendFirst; l_en < l_sendFirst+l_sendSize; l_en++ ) {
                if( i_gIds[l_en] == l_gId ) {
                  o_srcEl[l_is].push_back( l_en );
                  o_srcId[l_is].push_back( l_so );
                }
              }
            }
          }

        }
      }
    }

    /**
     * Updates the sparse types of the elements.
     *
     * @param i_srcEl source elements.
     * @param i_spTypeSrcs sparse type of the sources.
     * @param io_elChars will be updates with sparse types of sources if element has a source.
     *
     * @paramt intgral type of sparse types.
     **/
    template< typename TL_T_INT_SP >
    static void updateSpTypes( std::vector< TL_T_INT_LID >  i_srcEl[TL_N_IND_SRCS],
                               TL_T_INT_SP                  i_spTypeSrcs,
                               t_elementChars              *io_elChars ) {
      PP_INSTR_FUN("update_sp_types")

      // iterate over independent src configs
      for( unsigned short l_is = 0; l_is < TL_N_IND_SRCS; l_is++ ) {
        TL_T_INT_LID l_nSrcs = i_srcEl[l_is].size();

        // iterate over sources and set sparse type
#ifdef PP_USE_OMP
#pragma omp parallel for
#endif
        for( TL_T_INT_LID l_so = 0; l_so < l_nSrcs; l_so++ ) {
          io_elChars[ i_srcEl[l_is][l_so] ].spType |= i_spTypeSrcs;
        }
      }
    }

    /**
     * Prepares the mapping from sparse source elements to point sources.
     *
     * @param i_nEl number of dense elements.
     * @param i_spTypeSrcs sparse type of the sources.
     * @param i_elChars elements characteristics with sparse types set for source terms.
     * @param i_srcElP sparse source elements, sorted by their dense ids. If an elements holds more than one entry, multplie entries occur accordingly.
     * @param io_mem dynamic memory allocations.
     * @param o_elSpSo will be set to the mapping.
     *
     * @paramt integral type of sparse types.
     **/
    template< typename TL_T_INT_SP >
    static void prepareElSpSo( TL_T_INT_LID                         i_nEl,
                               TL_T_INT_SP                          i_spTypeSrcs,
                               t_elementChars              const   *i_elChars,
                               std::vector< TL_T_INT_LID > const    i_srcElP[TL_N_IND_SRCS],
                               data::Dynamic                       &io_mem,
                               TL_T_INT_LID                      (*&o_elSpSo)[TL_N_IND_SRCS] ) {
      // derive the number of source elements
      TL_T_INT_LID l_nElSrc = 0;

      // iterate over dense elements and determine size
      for( TL_T_INT_LID l_el = 0; l_el < i_nEl; l_el++ ) {
        if( (i_elChars[l_el].spType & i_spTypeSrcs) == i_spTypeSrcs ) l_nElSrc++;
      }

      // allocate memory for the mapping
      o_elSpSo =  (TL_T_INT_LID (*)[TL_N_IND_SRCS]) io_mem.allocate(
                    sizeof(TL_T_INT_LID) * TL_N_IND_SRCS * (l_nElSrc+1) );

      // init the first entries
      for( unsigned short l_kId = 0; l_kId < TL_N_IND_SRCS; l_kId++ ) {
        o_elSpSo[0][l_kId] = 0;
      }

      // sparse source element
      TL_T_INT_LID l_spElSrc = 0;

      // current point source index of the kinematic source descriptions
      TL_T_INT_LID l_ptSrc[TL_N_IND_SRCS];
      for( unsigned short l_kId = 0; l_kId < TL_N_IND_SRCS; l_kId++ ) l_ptSrc[l_kId] = 0;

      // iterate over dense elements and set mapping
      for( TL_T_INT_LID l_el = 0; l_el < i_nEl; l_el++ ) {
        if( (i_elChars[l_el].spType & i_spTypeSrcs) == i_spTypeSrcs ) {
          for( unsigned short l_kId = 0; l_kId < TL_N_IND_SRCS; l_kId++ ) {
            // set mapping
            o_elSpSo[l_spElSrc][l_kId] = l_ptSrc[l_kId];

            // update the index
            for( std::size_t l_so = l_ptSrc[l_kId]; l_so < i_srcElP[l_kId].size(); l_so++ ) {
              // check that we are ascending
              EDGE_CHECK_LE( l_el, i_srcElP[l_kId][l_so] );

              // increase index until we are in the next element
              if( i_srcElP[l_kId][l_so] == l_el ) l_ptSrc[l_kId]++;
              else                                break;
            }
          }

          // increase sparse source counter
          l_spElSrc++;
        }
      }

      // set ghost entries
      for( unsigned short l_kId = 0; l_kId < TL_N_IND_SRCS; l_kId++ ) {
        o_elSpSo[l_nElSrc][l_kId] = l_ptSrc[l_kId];
      }
    }

    /**
     * Prepares the local sources by allocating memory and by setting appropiate ids and boundaries.
     *
     * @param i_nSrcsG number of global point sources.
     * @param i_nSplsG number of global slip rate samples per dimension.
     * @param i_srcIdP dense source ids ordered by the memory layout.
     * @param i_srcElP dense element ids ordered by the memory layout.
     * @param i_firstG first global slip-rate entries of the sources per dimensions.
     * @param io_mem will be used for dynamic memory allocations.
     * @param o_solvers solvers which will be prepared.
     **/
    static void prepareSolvers( TL_T_INT_LID                                         i_nSrcsG,
                                std::size_t                           const          i_nSplsG[TL_N_DIM],
                                std::vector< TL_T_INT_LID >           const &        i_srcIdP,
                                std::vector< TL_T_INT_LID >           const &        i_srcElP,
                                unsigned int                          const * const  i_firstG[TL_N_DIM],
                                data::Dynamic                                       &io_mem,
                                solvers::t_Kinematics< TL_N_DIM,
                                                       TL_N_EL_MODES,
                                                       TL_N_FSRCS,
                                                       TL_T_REAL_COMP,
                                                       TL_T_INT_LID >               &o_solvers ) {
      PP_INSTR_FUN("prepare_solvers")

      // allocate memory for dense element ids
      o_solvers.soElDe = (TL_T_INT_LID*) io_mem.allocate( sizeof(TL_T_INT_LID) * i_srcElP.size() );
      // copy over dense element ids
      for( std::size_t l_so = 0; l_so < i_srcElP.size(); l_so++ )
        o_solvers.soElDe[l_so] = i_srcElP[l_so];

      // allocate memory for onset times
      o_solvers.onSet = (TL_T_REAL_COMP*) io_mem.allocate( sizeof(TL_T_REAL_COMP) * i_srcIdP.size() );

      // allocate memory for delta t of samplings
      o_solvers.dt = (TL_T_REAL_COMP*) io_mem.allocate( sizeof(TL_T_REAL_COMP) * i_srcIdP.size() );

      // allocate memory for evaluated basis functions
      o_solvers.bEval = (TL_T_REAL_COMP (*)[TL_N_EL_MODES]) io_mem.allocate(
                         sizeof(TL_T_REAL_COMP) * TL_N_EL_MODES * i_srcIdP.size() );

      // iterate over dimensions
      for( unsigned short l_di = 0; l_di < TL_N_DIM; l_di++ ) {
        // number of sources
        o_solvers.nSrcs = i_srcIdP.size();

        // only proceed if the slip-direction is active
        if( i_nSplsG[l_di] > 0 ) {
          // set active slip in this direction
          o_solvers.aSlip[l_di] = true;

          // allocate memory for position of first slip-rate sample for every source + ghost
          o_solvers.first[l_di] = (TL_T_INT_LID*) io_mem.allocate(   sizeof(TL_T_INT_LID)
                                                                   * (i_srcIdP.size()+1) );

          // allocate memory for scaling in the comptutations of the moment tensor
          o_solvers.sSca[l_di] = (TL_T_REAL_COMP (*)[TL_N_STRESS][TL_N_FSRCS] )
                                   io_mem.allocate(   sizeof(TL_T_REAL_COMP)
                                                    * TL_N_FSRCS*TL_N_DIM*2*i_srcIdP.size() );

          // local number of slip-rate samples
          std::size_t l_nSlpL = 0;

          // determine local number of slip-rate samples and set local first-entries
          o_solvers.first[l_di][0] = 0;
          for( std::size_t l_so = 1; l_so < i_srcIdP.size()+1; l_so++ ) {
            // init with previous source's offset
            o_solvers.first[l_di][l_so] = o_solvers.first[l_di][l_so-1];

            // global id of local source
            TL_T_INT_LID l_id;
            if( l_so < i_srcIdP.size() ) l_id = i_srcIdP[l_so];
            else                         l_id = i_srcIdP[l_so-1]+1;
            // global if of previous local source
            TL_T_INT_LID l_idP = i_srcIdP[l_so-1];

            if( l_id < i_nSrcsG ) {
              // #samples of previous source
              TL_T_INT_LID l_nSplsP = i_firstG[l_di][l_idP+1] - i_firstG[l_di][l_idP];

              o_solvers.first[l_di][l_so] += l_nSplsP;
              l_nSlpL += l_nSplsP;
            }
            else {
              // check for ghost-entry
              EDGE_CHECK_EQ( l_id, i_nSrcsG );

              // special handing for ghost-entry
              o_solvers.first[l_di][l_so] += i_nSplsG[l_di] - i_firstG[l_di][l_id-1];
              l_nSlpL += i_nSplsG[l_di] - i_firstG[l_di][l_id-1];
            }
          }

          // allocate memory for slip-rates dependent on the slip-rate sample
          o_solvers.sr[l_di] = (TL_T_REAL_COMP (*)[TL_N_FSRCS] )
                                 io_mem.allocate(   sizeof(TL_T_REAL_COMP)
                                                  * TL_N_FSRCS*l_nSlpL );
        }
        else o_solvers.aSlip[l_di] = false;
      }
    }

    /**
     * Sets up the point sources:
     *   1) Onset times
     *   2) dt of the pt-src samples.
     *   3) Slip-rates of the sources in the solver.
     *   4) Precomputed coefficients in the solvers for the entries in the moment tensor.
     *
     * Remark: The source specifications are accessed ascending for linear access
     *         when reading from disk. Random access happens when writing the coefficients
     *         of the moment tensors.
     *
     * @param i_kId id of the kinematic source in the source specifications.
     * @param i_srcs source specifications.
     * @param i_srcId ordered ids of the sources in the source specification.
     * @param i_srcEl dense element-ids w.r.t. the (ordered) sources.
     * @param i_permI inverse permutation. For example, index 0 with value 3 means: local source 0 goes to position 3, where local source 0 has the third smallest source id.
     * @param i_elVe vertices adjacent to elements.
     * @param i_veChars vertex characteristics, which are used to access the vertex coordinates.
     * @param i_elChars element characteristics, which are used for Lame parameter lambda and as fallback for mu if the value is below tolerance TOL.SOLVER.
     * @param i_bgPars background paramters. Lame parameters mu and lambda might be used if source reader does not provide them.
     * @parma i_massI inverse mass matrix.
     * @param io_solver will be updated with moment-tensor coefficients and slip-rates.
     *
     * @paramt TL_T_IN_SRC type of the source reader.
     **/
    template< typename TL_T_IN_SRC >
    static void setPtSrcs( unsigned short                                i_kId,
                           TL_T_IN_SRC                                  &i_srcs,
                           std::vector< TL_T_INT_LID >           const  &i_srcId,
                           std::vector< TL_T_INT_LID >           const  &i_srcEl,
                           std::vector< TL_T_INT_LID >           const  &i_permI,
                           TL_T_INT_LID                          const (*i_elVe)[TL_N_VE],
                           t_vertexChars                         const  *i_veChars,
                           t_elementChars                        const  *i_elChars,
                           t_bgPars                              const  *i_bgPars,
                           TL_T_REAL_COMP                        const   i_massI[TL_N_EL_MODES],
                           solvers::t_Kinematics< TL_N_DIM,
                                                  TL_N_EL_MODES,
                                                  TL_N_FSRCS,
                                                  TL_T_REAL_COMP,
                                                  TL_T_INT_LID >       &io_solvers ) {
      PP_INSTR_FUN("set_pt_srcs")

      // iterate over active sources
      for( std::size_t l_so = 0; l_so < i_srcId.size(); l_so++ ) {
        // get onset time
        io_solvers.onSet[ i_permI[l_so] ] = i_srcs.getOnSet( i_kId,
                                                             i_srcId[l_so] );
        // get dt
        io_solvers.dt[ i_permI[l_so] ] = i_srcs.getDt( i_kId,
                                                       i_srcId[l_so] );

        // get slip rates
        for( unsigned short l_sd = 0; l_sd < TL_N_DIM; l_sd++ ) {
          if( io_solvers.aSlip[l_sd] ) {
            // determine first entry
            TL_T_INT_LID l_first = io_solvers.first[l_sd][ i_permI[l_so] ];

            i_srcs.getSrs( i_kId,
                           i_srcId[l_so],
                           l_sd,
                           io_solvers.sr[l_sd][l_first] );

          }
        }

        /*
         * determine slip-rate coeffiecients
         */
        // slip direction
        double l_sds[TL_N_DIM][TL_N_DIM];
        i_srcs.getSds( i_kId,
                       i_srcId[l_so],
                       l_sds );

        // lame parameter mu
        TL_T_REAL_COMP l_mu = i_srcs.getMu( i_kId,
                                            i_srcId[l_so] );

        // overwrite mu by value of background velocity model if not provided
        if( l_mu < -TOL.SOLVER ) {
          l_mu = i_bgPars[ i_srcEl[l_so] ].mu;
        }

        // lame paramter lamda
        TL_T_REAL_COMP l_lam = i_bgPars[ i_srcEl[l_so] ].lam;

        // area
        TL_T_REAL_COMP l_A = i_srcs.getA( i_kId,
                                          i_srcId[l_so] );

        /* Slip coefficients:
         *   delta u * [ lambda * l_k * n_k * delta_{ij} + mu * ( l_i * n_j + l_j * n_i )
         *
         * delta u: slip cofficient is the term multiplied with it.
         * lambda: Lame parameter lambda
         * delta_{ij}: Kronecker delta
         * mu: Lame parameter mu
         * l: unit vector in slip-direction  l = ( l_1, l_2 ) or l = ( l_1, l_2, l_3 )
         * n: unit normal n = (n_1, n_2) or n = (n_1, n_2, n_3)
         *
         * Reference:
         *   Source Mechanisms of Earthquakes, 2014
         *   Theory and Practice
         *   Udias, Madariaga, Buforn
         *   Eq. (5.26)
         */
        for( unsigned short l_sd = 0; l_sd < TL_N_DIM; l_sd++ ) {
          if( io_solvers.aSlip[l_sd] ) {
            // init normal stresses
            for( unsigned short l_di = 0; l_di < TL_N_DIM; l_di++ ) {
              io_solvers.sSca[l_sd][ i_permI[l_so] ][l_di][0] = 0;
            }

            // scalar product with normal
            TL_T_REAL_COMP l_spN = 0;
            for( unsigned short l_di = 0; l_di < TL_N_DIM; l_di++ )
              l_spN += l_sds[l_sd][l_di]*l_sds[0][l_di];

            // contribution if slip is in normal direction
            for( unsigned short l_di = 0; l_di < TL_N_DIM; l_di++ ) {
              io_solvers.sSca[l_sd][ i_permI[l_so] ][l_di][0]
              += l_A * l_lam * l_spN;
            }

            // normal stresses
            for( unsigned short l_di = 0; l_di < TL_N_DIM; l_di++ ) {
              io_solvers.sSca[l_sd][ i_permI[l_so] ][l_di][0]
              += l_A * 2.0 * l_mu * l_sds[l_sd][l_di] * l_sds[0][l_di];
            }

            // moment tensor indices
            unsigned short l_mId[3][2] = { {0,1}, {1,2}, {0,2} };

            // shear stresses
            unsigned short l_nShear = TL_N_STRESS - TL_N_DIM;
            for( unsigned short l_di = 0; l_di < l_nShear; l_di++ ) {
              unsigned short l_mId1 = l_mId[l_di][0];
              unsigned short l_mId2 = l_mId[l_di][1];

              io_solvers.sSca[l_sd][ i_permI[l_so] ][l_di+TL_N_DIM][0]
              = l_A * l_mu * l_sds[l_sd][l_mId1] * l_sds[0][l_mId2] +
                l_A * l_mu * l_sds[l_sd][l_mId2] * l_sds[0][l_mId1];
            }

          }
        }

        /**
         * Evalute basis at the position of the point source and scale with inverse Jacobian.
         **/
        // get elements' ve-coords
        TL_T_REAL_MESH l_veCrds[3][TL_N_VE];
        mesh::common< TL_T_EL >::getElVeCoords( i_srcEl[l_so], i_elVe, i_veChars, l_veCrds );

        // get source coordinates
        TL_T_REAL_MESH l_srcCrds[TL_N_DIM];
        i_srcs.getSrcCrds( i_kId,
                           i_srcId[l_so],
                           l_srcCrds );

        // get position in ref element
        TL_T_REAL_MESH l_ref[3] = {0,0,0};
        edge::linalg::Mappings::phyToRef( TL_T_EL,
                                          l_veCrds[0],
                                          l_srcCrds,
                                          l_ref );
        // check for reasonable coords
        for( unsigned short l_di = 0; l_di < TL_N_DIM; l_di++ ) {
          EDGE_CHECK_GT( l_ref[l_di], -TOL.MESH );
          EDGE_CHECK_LT( l_ref[l_di], 1+TOL.MESH );
        }

        // setup the evaluated basis
        for( unsigned short l_md = 0; l_md < TL_N_EL_MODES; l_md++ ) {
          dg::Basis::evalBasis( l_md,
                                TL_T_EL,
                                io_solvers.bEval[ i_permI[l_so] ][l_md],
                                l_ref[0], l_ref[1], l_ref[2] );
        }

        // multiply with inverse mass matrix
        for( unsigned short l_md = 0; l_md < TL_N_EL_MODES; l_md++ ) {
          io_solvers.bEval[ i_permI[l_so] ][l_md] *= i_massI[l_md];
        }

        // divide by jacobian's det
        TL_T_REAL_MESH l_jac[TL_N_DIM][TL_N_DIM];
        edge::linalg::Mappings::evalJac( TL_T_EL, l_veCrds[0], l_jac[0] );

        TL_T_REAL_MESH l_jacDet;
        l_jacDet = edge::linalg::Matrix::det( TL_N_DIM, l_jac[0] );

        for( unsigned short l_md = 0; l_md < TL_N_EL_MODES; l_md++ )
          io_solvers.bEval[ i_permI[l_so] ][l_md] /= l_jacDet;
      }
    }

    /**
     * Derives the source permutations, which order the sources by the elements' ids.
     *
     * @param i_elId element ids, as given by the order of the point sources.
     * @param o_per will be set to the permutation. For example, index 0 with value 3 means: local source 3 goes to position 0, where 3 holds the smallest element id.
     * @param o_perI will be set to the inverse permutation. For example, index 0 with value 3 means: local source 0 goes to position 3, where local source 0 has the third smallest source id.
     **/
    static void getSrcPerms( std::vector< TL_T_INT_LID > const &i_elId,
                             std::vector< TL_T_INT_LID >       &o_per,
                             std::vector< TL_T_INT_LID >       &o_perI ) {
      PP_INSTR_FUN("get_src_perms")

      o_per.resize(  i_elId.size() );
      o_perI.resize( i_elId.size() );
      for( std::size_t l_el = 0; l_el < i_elId.size(); l_el++ ) {
        o_per[l_el] = l_el;
        o_perI[l_el] = l_el;
      }
      // determine permutation
      std::sort( o_per.begin(), o_per.end(),
                 [&](const int& a, const int& b) { return i_elId[a] < i_elId[b]; }
               );
      // determine inverse permutation
      std::sort( o_perI.begin(), o_perI.end(),
                 [&](const int& a, const int& b) { return o_per[a] < o_per[b]; }
               );
    }

  public:
    /**
     * Initializes the kinematic source solvers.
     *
     * @param i_elLayout dense element layout.
     * @param i_elVe vertices adjacent to the elements.
     * @param i_spTypeSrcs sparse type (bit-mask) of used for source elements.
     * @param i_veChars vertex characteristics.
     * @param i_gIds global ids of the elements.
     * @param i_srcs source reader.
     * @param i_massI inverse mass matrix.
     * @param i_bgPars background paramters. Lame parameters mu and lambda might be used if source reader does not provide them.
     * @param io_elChars element characteristics which will be update by the sparse source type if the element contains one or more point-sources.
     * @param io_mem will be used for dynamic memory allocations.
     * @param o_elSpSo will be set to first source ids of the sparse source elements.
     * @param o_solver will be set to kinematic source solvers.
     *
     * @paramt TL_T_INT_GID integral type of global ids.
     * @paramt TL_T_INT_SP integral type of sparse types.
     * @paramt TL_T_IN_SRC type of the source reader.
     **/
    template< typename TL_T_INT_GID,
              typename TL_T_INT_SP,
              typename TL_T_IN_SRC >
    static void solvers( t_enLayout                            const   &i_elLayout,
                         TL_T_INT_LID                          const  (*i_elVe)[TL_N_VE],
                         TL_T_INT_SP                                    i_spTypeSrcs,
                         t_vertexChars                         const   *i_veChars,
                         TL_T_INT_GID                          const   *i_gIds,
                         TL_T_IN_SRC                                   &i_srcs,
                         TL_T_REAL_COMP                        const    i_massI[TL_N_EL_MODES],
                         t_bgPars                              const   *i_bgPars,
                         t_elementChars                                *io_elChars,
                         data::Dynamic                                 &io_mem,
                         TL_T_INT_LID                                (*&o_elSpSo)[TL_N_IND_SRCS],
                         solvers::t_Kinematics< TL_N_DIM,
                                                TL_N_EL_MODES,
                                                TL_N_FSRCS,
                                                TL_T_REAL_COMP,
                                                TL_T_INT_LID >         o_solvers[TL_N_IND_SRCS] ) {
      PP_INSTR_FUN("solvers")

      // check that we have enough kinematic sources
      EDGE_CHECK_LE( TL_N_IND_SRCS, i_srcs.nIn() );

      // get the number of global sources
      std::vector< TL_T_INT_LID > l_nSrcsG( TL_N_IND_SRCS );
      for( unsigned short l_is = 0; l_is < TL_N_IND_SRCS; l_is++ ) {
        std::size_t l_tmp = i_srcs.nSrcsG( l_is );

        // make sure that our local ids are large enough to fit all sources.
        // if this fails, our implementation has bigger problems than the integer types..
        EDGE_CHECK_LT( l_tmp, std::numeric_limits< TL_T_INT_LID >::max() );

        l_nSrcsG[l_is] = l_tmp;
      }

      // query source coordinates
      double (*l_srcCrds[TL_N_IND_SRCS])[TL_N_DIM];
      for( unsigned short l_is = 0; l_is < TL_N_IND_SRCS; l_is++ ) {
        l_srcCrds[l_is] = new double[ l_nSrcsG[l_is] ][TL_N_DIM];

        i_srcs.getSrcCrds( l_is, l_srcCrds[l_is], 0, l_nSrcsG[l_is] );
      }

      // determine sources local to the MPI-domain
      std::vector< TL_T_INT_LID > l_srcEl[TL_N_IND_SRCS]; // dense element ids
      std::vector< TL_T_INT_LID > l_srcId[TL_N_IND_SRCS]; // dense source ids

      getLocal( &l_nSrcsG[0],
                 l_srcCrds,
                 i_elLayout,
                 i_elVe,
                 i_veChars,
                 i_gIds,
                 l_srcEl,
                 l_srcId );

      // determine permutations
      std::vector< TL_T_INT_LID > l_perm[TL_N_IND_SRCS];
      std::vector< TL_T_INT_LID > l_permI[TL_N_IND_SRCS];
      // permutations of source ids and source elements w.r.t. to the memory layout
      std::vector< TL_T_INT_LID > l_srcIdP[TL_N_IND_SRCS];
      std::vector< TL_T_INT_LID > l_srcElP[TL_N_IND_SRCS];

      for( unsigned short l_is = 0; l_is < TL_N_IND_SRCS; l_is++ ) {
        getSrcPerms( l_srcEl[l_is], l_perm[l_is], l_permI[l_is] );
        // apply forward permutation
        l_srcIdP[l_is].resize( l_perm[l_is].size() );
        l_srcElP[l_is].resize( l_perm[l_is].size() );
        for( std::size_t l_so = 0; l_so < l_srcId[l_is].size(); l_so++ ) {
          TL_T_INT_LID l_soP = l_perm[l_is][l_so];
          l_srcIdP[l_is][l_so] = l_srcId[l_is][l_soP];
          l_srcElP[l_is][l_so] = l_srcEl[l_is][l_soP];
        }
      }

      // update the sparse types
      updateSpTypes( l_srcEl,
                     i_spTypeSrcs,
                     io_elChars );

      // prepare mapping from elements to sources
      prepareElSpSo( i_elLayout.nEnts,
                     i_spTypeSrcs,
                     io_elChars,
                     l_srcElP,
                     io_mem,
                     o_elSpSo );

      for( unsigned short l_is = 0; l_is < TL_N_IND_SRCS; l_is++ ) {
        // global number of slip-rate samples per dimension
        std::size_t l_nSplsG[TL_N_DIM];

        // temporary memory holding all first global slip-rate entries of the sources
        unsigned int *l_firstG[TL_N_DIM];

        // get first-entries for all global sources
        for( unsigned short l_di = 0; l_di < TL_N_DIM; l_di++ ) {
          l_nSplsG[l_di] = i_srcs.nSplsG( l_is, l_di );

          l_firstG[l_di] = new unsigned int[ l_nSrcsG[l_is]+1 ];
          i_srcs.getOffSetsG( l_is, l_di, l_firstG[l_di] );
        }

        // prepare memory and sizes
        prepareSolvers( l_nSrcsG[l_is],
                        l_nSplsG,
                        l_srcIdP[l_is],
                        l_srcElP[l_is],
                        l_firstG,
                        io_mem,
                        o_solvers[l_is] );

        // set moment tensor coefficients and slip rates
        setPtSrcs( l_is,
                   i_srcs,
                   l_srcId[l_is],
                   l_srcEl[l_is],
                   l_permI[l_is],
                   i_elVe,
                   i_veChars,
                   io_elChars,
                   i_bgPars,
                   i_massI,
                   o_solvers[l_is] );

        // free temporary offsets
        for( unsigned short l_di = 0; l_di < TL_N_DIM; l_di++ ) {
         delete[] l_firstG[l_di];
        }
     }

      // free memory
      for( unsigned short l_is = 0; l_is < TL_N_IND_SRCS; l_is++ ) {
        delete[] l_srcCrds[l_is];
      }
    }
};

#endif
