/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2016-2017, Regents of the University of California
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
 * Finite Volume solver for elastic wave equations.
 **/
#ifndef FINITE_VOLUME_HPP
#define FINITE_VOLUME_HPP

#include "constants.hpp"
#include "io/Receivers.h"
#include "InternalBoundary.hpp"
#include "FrictionLaws.hpp"

namespace edge {
  namespace elastic {
    namespace solvers {
      class FiniteVolume;
    }
  }
}

class edge::elastic::solvers::FiniteVolume {
  //private:
    /**
     * Updates the DOFs with the contribution of a flux solver.
     *
     * @param i_a flux solver.
     * @param i_b dofs.
     * @param io_c new dofs, updated with the contribution.
     **/
    static void computeFluxSolver( const real_base  i_a[N_QUANTITIES][N_QUANTITIES],
                                   const real_base  i_b[N_QUANTITIES][N_ELEMENT_MODES][N_CRUNS],
                                         real_base io_c[N_QUANTITIES][N_ELEMENT_MODES][N_CRUNS] ) {
#if __has_builtin(__builtin_assume_aligned)
      // share alignment with compiler
      (void) __builtin_assume_aligned(i_b, ALIGNMENT.ELEMENT_MODES.PRIVATE);
      (void) __builtin_assume_aligned(io_c, ALIGNMENT.ELEMENT_MODES.PRIVATE);
#endif

      // iterate over rows of the flux solver
      for( unsigned int l_row = 0; l_row < N_QUANTITIES; l_row++ ) {
        // iterate over columns of the flux solver
        for( unsigned int l_col = 0; l_col < N_QUANTITIES; l_col++ ) {
          // iteare over concurrent runs
#pragma omp simd
          for( int_cfr l_crun = 0; l_crun < N_CRUNS; l_crun++ ) {
            io_c[l_row][0][l_crun] += i_a[l_row][l_col] * i_b[l_col][0][l_crun];
          }
        }
      }
    }

  public:
    /**
     * Compute (trivial) time prediction.
     *
     * @param i_first first element cosidered.
     * @param i_nElements number of elements.
     * @param i_time time of the DOFs.
     * @param i_dT time step ("integrated interval")
     * @param i_firstSpRp first sparse rupture-element.
     * @param i_firstSpRc first sparse receiver-element.
     * @param i_elChars element characteristics.
     * @param i_dofs DOFs.
     * @param o_tInt will be set to time integrated DOFs.
     * @param o_tRup will be set to DOFs ("derivatives") for rupture elements.
     * @param io_recvs will be updated with receiver info.
     *
     * @paramt TL_T_INT_LID integer type of local entity ids.
     * @paramt TL_T_REAL floating point type.
     **/
    template < typename TL_T_INT_LID,
               typename TL_T_REAL >
    static void timePred(       TL_T_INT_LID          i_first,
                                TL_T_INT_LID          i_nElements,
                                double                i_time,
                                double                i_dT,
                                TL_T_INT_LID          i_firstSpRp,
                                TL_T_INT_LID          i_firstSpRc,
                                t_elementChars       *i_elChars,
                          const TL_T_REAL           (*i_dofs)[N_QUANTITIES][1][N_CRUNS],
                                TL_T_REAL           (*o_tInt)[N_QUANTITIES][1][N_CRUNS],
                                TL_T_REAL           (*o_tRup)[N_QUANTITIES][1][N_CRUNS],
                                edge::io::Receivers  &io_recvs ) {
#if __has_builtin(__builtin_assume_aligned)
      // share alignment with compiler
      (void) __builtin_assume_aligned(i_dofs, ALIGNMENT.ELEMENT_MODES.PRIVATE);
      (void) __builtin_assume_aligned(o_tInt, ALIGNMENT.ELEMENT_MODES.PRIVATE);
#endif

      // counter of elements with faces having rupture physics
      TL_T_INT_LID l_elRp = i_firstSpRp;

      // counter of elements with receivers
      TL_T_INT_LID l_elRe = i_firstSpRc;

      // iterate over elements
      for( TL_T_INT_LID l_el = i_first; l_el < i_first+i_nElements; l_el++ ) {
        // update receivers if required
        if( !( (i_elChars[l_el].spType & RECEIVER) == RECEIVER) ) {}
        else {
          while( true ) {
            double l_reDt = io_recvs.getRecvTimeRel( l_elRe, i_time, i_dT );
            if( !(l_reDt >= 0) ) break;
            else                 io_recvs.writeRecvAll( l_elRe, i_dofs[l_el] );
          }
          l_elRe++;
        }

        // compute time integrated DOFs
        for( int_qt l_qt = 0; l_qt < N_QUANTITIES; l_qt++ ) {
#pragma omp simd
          for( int_cfr l_ru = 0; l_ru < N_CRUNS; l_ru++ ) {
            o_tInt[l_el][l_qt][0][l_ru] = i_dT * i_dofs[l_el][l_qt][0][l_ru];
          }
        }

        // copy dofs to sparse "derivative" of rupture solver (if required)
        if( (i_elChars[l_el].spType & RUPTURE) != RUPTURE ) {}
        else {
          for( int_qt l_qt = 0; l_qt < N_QUANTITIES; l_qt++ ) {
            for( int_cfr l_ru = 0; l_ru < N_CRUNS; l_ru++ ) {
              o_tRup[l_elRp][l_qt][0][l_ru] = i_dofs[l_el][l_qt][0][l_ru];
            }
          }

          // update sparse rupture-element counter
          l_elRp++;
        }
      }
    }

    /**
     * Solves rupture physics for the given faces.
     *
     * @param i_first first rupture faces.
     * @param i_nFaces number of faces.
     * @param i_firstReSp first sparse id of the receivers.
     * @param i_time current time of the faces.
     * @param i_dT time step of the two adjacent elements.
     * @param i_faElSpRp adjacency information from sparse rupture faces to sparse rupture elements.
     * @param i_spTypesFaSp sparse types of the sparse rupture faces.
     * @param i_frictionGlobal global data of the friction law.
     * @param i_frictionFace face-local data of the friction law.
     * @param i_frictionQuadPoint data of the friction law local to the spatial quadrature points of the faces.
     * @param i_solver solvers used for internal boundaries. [0]: trafo to face-align coords, [1]: middle states [2]: trafo to physical coords and flux (left side) [3]: trafo to physical coords and flux (right side).
     * @param i_tDofs DOFs of the sparse rupture elements.
     * @param o_updates will be set to element-updates resulting from the flux computation. [0]: left side, [1]: right side.
     * @param io_recvsQuad receivers at quadrature points.
     *
     * @paramt TL_T_INT_LID integer type of local entity ids.
     * @paramt TL_T_REAL type used for floating point arithmetic.
     * @paramt TL_T_FRI_GL struct representing global friction data.
     * @paramt TL_T_FRI_FA struct representing face-local friction data.
     * @paramt TL_T_FRI_QP struct representing quad-point-local friction data.
     * @paramt TL_T_INT_SP integer type of the sparse type.
     * @paramt TL_T_RECVQ type of the quadrature receivers.
     **/
    template< typename TL_T_INT_LID,
              typename TL_T_REAL,
              typename TL_T_FRI_GL,
              typename TL_T_FRI_FA,
              typename TL_T_FRI_QP,
              typename TL_T_INT_SP,
              typename TL_T_RECVQ >
    static void rupture( TL_T_INT_LID                     i_first,
                         TL_T_INT_LID                     i_nFaces,
                         TL_T_INT_LID                     i_firstReSp,
                         TL_T_REAL                        i_time,
                         TL_T_REAL               const    i_dT,
                         TL_T_INT_LID            const (* i_faElSpRp)[2],
                         t_InternalBoundaryFace<
                           TL_T_REAL,
                           TL_T_INT_SP
                         >                       const  * i_iBnd,
                         TL_T_FRI_GL                    & i_frictionGlobal,
                         TL_T_FRI_FA                    * i_frictionFace,
                         TL_T_FRI_QP                   (* i_frictionQuadPoint)[1],
                         TL_T_REAL               const (* i_solvers)[4][N_QUANTITIES][N_QUANTITIES],
                         TL_T_REAL               const (* i_tDofs)[N_QUANTITIES][1][N_CRUNS],
                         TL_T_REAL                     (* o_updates)[2][N_QUANTITIES][1][N_CRUNS],
                         TL_T_RECVQ                     & io_recvsQuad ) {
      // store sparse receiver id
      TL_T_INT_LID l_faRe = i_firstReSp;

      // set up "dummy" values for high-order discretization
      TL_T_REAL l_massI[1]            = {1 / *C_REF_ELEMENT.VOL.ENT[T_SDISC.ELEMENT] };
      TL_T_REAL l_ptsLine[1]          = {0.5};
      TL_T_REAL l_weightsLine[1]      = {1};
      TL_T_REAL l_weightsFaces[1]     = {1};
      unsigned short const l_quadOpts =   (CE_N_FACE_VERTEX_OPTS(T_SDISC.ELEMENT)+1)
                                        *  C_ENT[T_SDISC.ELEMENT].N_FACES;
      TL_T_REAL l_basisFaces[l_quadOpts][1][1];
      for( unsigned short l_fa = 0; l_fa < l_quadOpts; l_fa++ ) {
        l_basisFaces[l_fa][0][0] = 1;
      }

      // scratch memory
      TL_T_REAL l_scratch[4][N_QUANTITIES][1][N_CRUNS];

      // struct for the pertubation of the middle states
      struct {
        TL_T_FRI_GL  *gl;
        TL_T_FRI_FA  *fa;
        TL_T_FRI_QP (*qp)[1];
      } l_faData;
      l_faData.gl = &i_frictionGlobal;
      l_faData.fa =  i_frictionFace+i_first;
      l_faData.qp =  i_frictionQuadPoint+i_first;

      // iterate over the rupture faces
      for( TL_T_INT_LID l_fa = i_first; l_fa < i_first+i_nFaces; l_fa++ ) {
        // get sparse rupture elements at the left and right side of the face
        TL_T_INT_LID l_spL = i_faElSpRp[l_fa][0];
        TL_T_INT_LID l_spR = i_faElSpRp[l_fa][1];

        // call internal boundary solver
        InternalBoundary<
          T_SDISC.ELEMENT,
          N_QUANTITIES,
          1,
          N_CRUNS >::template evalSpaceTime<
            TL_T_REAL,
            FrictionLaws< N_DIM, N_CRUNS >
          >(  0, 0, 0, // irrelevant for FV
              l_massI,
              i_dT,
              l_ptsLine,
              l_weightsLine,
              l_weightsFaces,
              l_basisFaces,
              i_solvers[l_fa][0],
              i_solvers[l_fa][1],
              i_solvers[l_fa][2],
              i_solvers[l_fa][3],
              i_tDofs+l_spL, // derivative (simply the DOF) "computation"
              i_tDofs+l_spR, // derivative (simply the DOF) "computation"
              l_scratch,
              o_updates[l_fa][0],
              o_updates[l_fa][1],
             &l_faData );

        // check if this face requires receiver output
        if( (i_iBnd[l_fa].spType & RECEIVER) != RECEIVER ){}
        else {
          // check if the receiver requires output
          if( io_recvsQuad.getRecvTimeRel( l_faRe, i_time, i_dT ) >= -TOL.TIME ) {
            // gather receiver data, TODO: outsource
            TL_T_REAL l_buff[ (N_DIM-1)*3 ][N_FACE_QUAD_POINTS][N_CRUNS];

            for( unsigned short l_qp = 0; l_qp < N_FACE_QUAD_POINTS; l_qp++ ) {
              for( unsigned short l_di = 0; l_di < N_DIM-1; l_di++ ) {
                for( unsigned short l_ru = 0; l_ru < N_CRUNS; l_ru++ ) {
                  l_buff[          0+l_di][l_qp][l_ru] = i_frictionQuadPoint[l_fa][l_qp].tr[l_di][l_ru];
                  l_buff[  (N_DIM-1)+l_di][l_qp][l_ru] = i_frictionQuadPoint[l_fa][l_qp].sr[l_di][l_ru];
                  l_buff[2*(N_DIM-1)+l_di][l_qp][l_ru] = i_frictionQuadPoint[l_fa][l_qp].dd[l_di][l_ru];
                }
              }
            }

            // write the receiver info (potentially multiple time if time step is larger than receiver frequency)
            io_recvsQuad.writeRecvAll( i_time, i_dT, l_faRe, l_buff );
          }

          l_faRe++;
        }

        // update pointers
        l_faData.fa++;
        l_faData.qp++;
      }
    }

    /**
     * Updates the DOFs using the flux solvers.
     *
     * @param i_first first element.
     * @param i_nElements number of elements.
     * @param i_firstSpRp first sparse rupture element.
     * @param i_elFaEl face-adjacent elements.
     * @param i_elFa faces adjacent to elements.
     * @param i_faElSpRp adjacency information from sparse rupture faces to sparse rupture elements.
     * @param i_elFaSpRp adjacnecy information from sparse rupture elements to sparse rupture faces.
     * @param i_fluxSolversOwn flux solvers for the elements' own contribution.
     * @param i_fluxSolversNeigh flux solvers for the neighboring elements' contribution.
     * @param i_tInt time integrated DOFs.
     * @param i_updatesSpRp updates resulting from the numerical flux at rupture faces.
     * @param io_dofs will updated by the flux contributions.
     *
     * @paramt TL_T_INT_LID integer type of local entity ids.
     * @paramt TL_T_REAL type used for floating point arithmetic.
     **/
    template< typename TL_T_INT_LID,
              typename TL_T_REAL >
    static void update( TL_T_INT_LID         i_first,
                        TL_T_INT_LID         i_nElements,
                        TL_T_INT_LID         i_firstSpRp,
                        TL_T_INT_LID const (*i_elFaEl)[C_ENT[T_SDISC.ELEMENT].N_FACES],
                        TL_T_INT_LID const (*i_elFa)[C_ENT[T_SDISC.ELEMENT].N_FACES],
                        TL_T_INT_LID const (*i_faElSpRp)[2],
                        TL_T_INT_LID const (*i_elFaSpRp)[C_ENT[T_SDISC.ELEMENT].N_FACES],
                        t_faceChars  const (*i_faChars),
                        t_fluxSolver const (*i_fluxSolversOwn)[C_ENT[T_SDISC.ELEMENT].N_FACES],
                        t_fluxSolver const (*i_fluxSolversNeigh)[C_ENT[T_SDISC.ELEMENT].N_FACES],
                        TL_T_REAL    const (*i_tInt)[N_QUANTITIES][1][N_CRUNS],
                        TL_T_REAL    const (*i_updatesSpRp)[2][N_QUANTITIES][1][N_CRUNS],
                        TL_T_REAL          (*io_dofs)[N_QUANTITIES][1][N_CRUNS] ) {
#if __has_builtin(__builtin_assume_aligned)
      // share alignment with compiler
      (void) __builtin_assume_aligned(i_tInt, ALIGNMENT.ELEMENT_MODES.PRIVATE);
      (void) __builtin_assume_aligned(io_dofs, ALIGNMENT.ELEMENT_MODES.PRIVATE);
#endif

      // counter of elements with faces having rupture physics
      TL_T_INT_LID l_elRp = i_firstSpRp;

      // iterate over elements
      for( TL_T_INT_LID l_el = i_first; l_el < i_first+i_nElements; l_el++ ) {
        // will be set to true if one of the faces enforces rupture physics
        bool l_rp = false;

        // iterate over faces
        for( unsigned int l_fa = 0; l_fa < C_ENT[T_SDISC.ELEMENT].N_FACES; l_fa++ ) {
          TL_T_INT_LID l_faId = i_elFa[l_el][l_fa];

          // call flux solver for the element's own contribution (if not a rupture face)
          if( (i_faChars[l_faId].spType & RUPTURE) != RUPTURE ) {
            computeFluxSolver( i_fluxSolversOwn[l_el][l_fa].solver,
                               i_tInt[l_el],
                               io_dofs[l_el] );
          }

          // call flux solver for neighboring elements' contributions (if not a outflow boundary or rupture face)
          if(    ( (i_faChars[l_faId].spType & OUTFLOW) != OUTFLOW )
              && ( (i_faChars[l_faId].spType & RUPTURE) != RUPTURE ) ){
            // determine neighbor
            TL_T_INT_LID l_neigh;

            if( (i_faChars[l_faId].spType & FREE_SURFACE) != FREE_SURFACE ) l_neigh = i_elFaEl[l_el][l_fa];
            else                                                            l_neigh = l_el;

            // call flux solver for the neighboring element's contribution
            computeFluxSolver( i_fluxSolversNeigh[l_el][l_fa].solver,
                               i_tInt[l_neigh],
                               io_dofs[l_el] );
          }
          // apply updates of the flux computation directly if this is a rupture face
          else if( ( (i_faChars[l_faId].spType & RUPTURE) == RUPTURE ) ) {
            // get the sparse index of the rupture face
            TL_T_INT_LID l_faIdRp = i_elFaSpRp[l_elRp][l_fa];

            // determine if this is the left or right side element
            TL_T_INT_LID l_elL = i_faElSpRp[l_faIdRp][0];
            TL_T_INT_LID l_sd  = (l_elRp == l_elL) ? 0 : 1;

            // update the DOFs
            for( unsigned short l_qt = 0; l_qt < N_QUANTITIES; l_qt++ ) {
              for( unsigned short l_ru = 0; l_ru < N_CRUNS; l_ru++ ) {
                io_dofs[l_el][l_qt][0][l_ru] += i_updatesSpRp[l_faIdRp][l_sd][l_qt][0][l_ru];
              }
            }

            // remember this rupture face to increase the rupture element counter
            l_rp = true;
          }
        }

        // increase the sparse counter for the rupture elements
        if( l_rp == false ){}
        else l_elRp++;
      }
    }
};

#endif
