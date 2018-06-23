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
 * Super-cells (collection of sub-cells).
 **/
#ifndef EDGE_SC_SUPERCELL_HPP
#define EDGE_SC_SUPERCELL_HPP

#include "constants.hpp"
#include "io/logging.h"
#include <vector>
#include "SuperCell.type"

namespace edge {
  namespace sc {
    namespace ibnd {
      template< t_entityType   TL_T_EL,
                unsigned short TL_O_SP >
      class SuperCellInit;

      template< t_entityType   TL_T_EL,
                unsigned short TL_O_SP,
                unsigned short TL_N_QTS,
                unsigned short TL_N_CRS >
      class SuperCell;
    }
  }
}

/**
 * Initialization of the super-cells.
 *
 * @paramt TL_T_EL element type.
 * @paramt TL_O_SP spatial order of the DG-scheme.
 **/
template< t_entityType   TL_T_EL,
          unsigned short TL_O_SP >
class edge::sc::ibnd::SuperCellInit {
  private:
    //! number of vertices per DG-face / sub-face
    static unsigned short const TL_N_VES_FA = C_ENT[TL_T_EL].N_FACE_VERTICES;

    //! number of faces
    static unsigned short const TL_N_FAS = C_ENT[TL_T_EL].N_FACES;

    //! number of sub-faces
    static unsigned short const TL_N_SFS = CE_N_SUB_FACES( TL_T_EL, TL_O_SP );

    //! number of sub-cells per element
    static unsigned short const TL_N_SCS  = CE_N_SUB_CELLS( TL_T_EL, TL_O_SP );

    /**
     * Assembles a super-cell stencil.
     *   Output is sorted by the id of the limited plus elements, but not by the sub-cell ids.
     * 
     * @param i_scSrc sub-cell of the source limited plus element, which is part of the stencil.
     * @param i_scDgAd sub-cell adjacency when moving from one element to the next (faces as bridge).
     * @param i_scSfSc sub-cell adjacency within an element (faces as bridge).
     * @param i_vIdElFaEl vertex ids of the element faces.
     * @param i_fIdElFaEl face ids of the element faces.
     * @param i_lpEl limited plus to element (no bridge).
     * @param i_lpFaLp limited plus adjacent to limited plus (faces as bridge).
     * @param i_nLp number of limited plus elements, possibly involved in the stencil.
     * @param i_lps ids of the limited plus elements, possibly involved in the stencil.
     * @param i_iBndFas local face ids of the limited plus elements, which build the internal boundary.
     * @param o_suLp will be set to limited plus elements, which build the stencil (1 lp per stencil component).
     * @param o_scs will be set to local sub-cell ids for of the stencil
     *
     * @paramt TL_T_LID integral type of local ids.
     **/
    template< typename TL_T_LID >
    static void stAssemble( unsigned short                        i_scSrc,
                            unsigned short                const   i_scDgAd[TL_N_VES_FA][TL_N_SFS],
                            unsigned short                const   i_scSfSc[TL_N_SCS + TL_N_FAS * TL_N_SFS][TL_N_FAS],
                            unsigned short                const (*i_vIdElFaEl)[TL_N_FAS],
                            unsigned short                const (*i_fIdElFaEl)[TL_N_FAS],
                            TL_T_LID                      const  *i_lpEl,
                            TL_T_LID                      const (*i_lpFaLp)[TL_N_FAS],
                            unsigned short                        i_nLp,
                            TL_T_LID                      const  *i_lps,
                            unsigned short                const  *i_iBndFas,
                            std::vector< TL_T_LID       >        &o_suLp,
                            std::vector< unsigned short >        &o_scs ) {
      EDGE_CHECK_GT( i_nLp, 0 );

      // init output
      o_suLp.resize(0);
      o_scs.resize(0);

      // involved super-elements
      std::vector< TL_T_LID > l_suLp;
      // involved sub-cells
      std::vector< unsigned short > l_scs;
      // involved boundary faces
      std::vector< unsigned short > l_iBndFas;

      // init
      l_suLp.push_back( i_lps[0] );
      l_scs.push_back( i_scSrc );
      l_iBndFas.push_back( i_iBndFas[0] );

      // first entry of the previous level
      unsigned short l_firstLvl = 0;

      while( l_firstLvl != l_scs.size() ){
         unsigned short l_lvlSize = l_scs.size() - l_firstLvl;

        // add everything at the internal boundary, which is local
        for( unsigned short l_sc = l_firstLvl; l_sc < l_firstLvl+l_lvlSize; l_sc++ ) {
          unsigned short l_scId   = l_scs[l_sc];
          unsigned short l_lp     = l_suLp[l_sc];
          unsigned short l_iBndFa = l_iBndFas[l_sc];

          // iterate over adjacent sub-cells
          for( unsigned short l_s1 = 0; l_s1 < TL_N_FAS; l_s1++ ) {
            unsigned short l_scAd = i_scSfSc[l_scId][l_s1];

            // ignore ghost sub-cels
            if( l_scAd >= TL_N_SCS ) continue;

            // check if the adjacent sub-cell is at the internal boundary
            bool l_bnd = false;
            for( unsigned short l_s2 = 0; l_s2 < TL_N_FAS; l_s2++ ) {
              if(    i_scSfSc[l_scAd][l_s2] >= TL_N_SCS +  l_iBndFa    * TL_N_SFS
                  && i_scSfSc[l_scAd][l_s2]  < TL_N_SCS + (l_iBndFa+1) * TL_N_SFS ) l_bnd = true;
            }

            // add the adjacent sub-cell to the stencil (if not already there)
            if( l_bnd ) {
              // determine, if already part of the stencil
              bool l_present = false;
              for( unsigned short l_cu = l_firstLvl; l_cu < l_firstLvl+l_lvlSize; l_cu++ ) {
                if(    l_scs[l_cu]  == l_scAd
                    && l_suLp[l_cu] == l_lp ) l_present = true;
              }

              // add if not already present
              if( !l_present ) {
                l_suLp.push_back( l_lp );
                l_scs.push_back( l_scAd );
                l_iBndFas.push_back( l_iBndFa );
              }
            }

          }
        }

        // next level-start, since we are adding sub-cells of adjacent elements now
        unsigned short l_firstLvlNext = l_scs.size();

        // iterate over current levels sub-cells and try to find new, DG-adjacent sub-cells
        for( unsigned short l_s1 = l_firstLvl; l_s1 < l_firstLvlNext; l_s1++ ) {
          unsigned short l_scId   = l_scs[l_s1];
          unsigned short l_lp     = l_suLp[l_s1];

          // iterate over faces of the limited plus element and see if we can reach a new sub-cell
          for( unsigned short l_f1 = 0; l_f1 < TL_N_FAS; l_f1++ ) {
            TL_T_LID l_lpAd = i_lpFaLp[l_lp][l_f1];

            // determine internal boundary face and check if this limited plus element is part of our search range
            unsigned short l_iBndFaAd = std::numeric_limits< unsigned short >::max();
            for( unsigned short l_sr = 0; l_sr < i_nLp; l_sr++ )
              if( l_lpAd == i_lps[l_sr] ) l_iBndFaAd = i_iBndFas[l_sr];

            // abort if the face is not part of our search range
            if( l_iBndFaAd == std::numeric_limits< unsigned short >::max() ) continue;

            // destination face
            unsigned short l_faDe = l_f1;
            // vertex id
            unsigned short l_vId = i_vIdElFaEl[ i_lpEl[l_lp] ][ l_faDe ];
            // src face w.r.t. to next element
            unsigned short l_faSrc = i_fIdElFaEl[ i_lpEl[l_lp] ][ l_faDe ];

            // iterate over sub-faces of the DG-face
            for( unsigned short l_sf = 0; l_sf < TL_N_SFS; l_sf++ ) {
              unsigned short l_sfSc = i_scSfSc[ TL_N_SCS + l_faDe * TL_N_SFS + l_sf ][0];

              // continue if our sub-face is part of the DG-face
              if( l_sfSc == l_scId ) {
                // perform sub-face reordering
                unsigned short l_sfRe = i_scDgAd[ l_vId ][l_sf];

                // determine remote sub-cell
                unsigned short l_scRe = i_scSfSc[ TL_N_SCS + l_faSrc * TL_N_SFS + l_sfRe ][0];

                // determine, if already part of the stencil
                bool l_present = false;
                for( unsigned short l_cu = 0; l_cu < l_scs.size(); l_cu++ ) {
                  if(    l_scs[l_cu]  == l_scRe
                      && l_suLp[l_cu] == l_lpAd ) l_present = true;
                }

                // add if not present
                if( !l_present ) {
                  l_suLp.push_back( l_lpAd );
                  l_scs.push_back( l_scRe );
                  l_iBndFas.push_back( l_iBndFaAd );
                }
              }
            }
          }

        }

        // store the counter for the current level
        l_firstLvl = l_firstLvlNext;
      }

      // perform bloody inefficient reordering by limited plus ids (shouldn't be too large anyways..)
      TL_T_LID l_lpMin = l_suLp[0];
      TL_T_LID l_lpMax = l_suLp[0];
      for( unsigned short l_lp = 1; l_lp < l_suLp.size(); l_lp++ ) {
        l_lpMin = std::min( l_suLp[l_lp], l_lpMin );
        l_lpMax = std::max( l_suLp[l_lp], l_lpMax );
      }

      while( true ) {
        // perform reordering
        for( unsigned short l_sc = 0; l_sc < l_suLp.size(); l_sc++ ) {
          if( l_suLp[l_sc] == l_lpMin ) {
             o_suLp.push_back( l_suLp[l_sc] );
             o_scs.push_back( l_scs[l_sc] );
          }
        }

        // abort if min is max
        if( l_lpMin == l_lpMax ) break;

        // get new min
        TL_T_LID l_newMin =  l_lpMax;
        for( unsigned short l_lp = 0; l_lp < l_suLp.size(); l_lp++ ) {
          if(    l_suLp[l_lp] < l_newMin
              && l_suLp[l_lp] > l_lpMin ) l_newMin = l_suLp[l_lp];
        }
        l_lpMin = l_newMin;
      }

    }

  public:
    /**
     * Initializes the faces stencils for internal boundaries.
     *
     * @param i_scSfSc sub-cells adjacent to sub-cells (sub-faces as bridge).
     * @param i_scTySf type of the sub-cells faces.
     **/
    static void faSt( unsigned short const    i_scDgAd[TL_N_VES_FA][TL_N_SFS],
                      unsigned short const    i_scSfSc[TL_N_SCS + TL_N_FAS * TL_N_SFS][TL_N_FAS],
                      t_SuperCell< TL_T_EL,
                                   TL_O_SP > &o_su ) {
      // abort for non-tets
      if( TL_T_EL != TET4 ) return;

      // iterate over faces
      for( unsigned short l_fa = 0; l_fa < TL_N_FAS; l_fa++ ) {
        // counter for stencils
        unsigned short l_nSt = 0;

        // iterate over sub-faces at the faces
        for( unsigned short l_sf = 0; l_sf < TL_N_SFS; l_sf++ ) {
          unsigned short l_sc = i_scSfSc[TL_N_SCS + l_fa * TL_N_SFS + l_sf ][0];

          // dummy limited plus elements of the stencil
          std::vector< short > l_suLp;

          // sub-cells of the stencil
          std::vector< unsigned short > l_scs;

          // dummy data for the query
          short l_lp = 0;
          short l_lpFaLp[TL_N_FAS];
          for( unsigned short l_en = 0; l_en < TL_N_FAS; l_en++ )
            l_lpFaLp[l_en] = -1;

          // call stencil assembly with dummy arguments
          stAssemble(                                   l_sc,
                                                        i_scDgAd,
                                                        i_scSfSc,
                      (unsigned short (*) [TL_N_FAS] )  nullptr,
                      (unsigned short (*) [TL_N_FAS] )  nullptr,
                      (short *)                         nullptr,
                                                       &l_lpFaLp,
                      (unsigned short)                  1,
                                                       &l_lp,
                                                       &l_fa,
                                                        l_suLp,
                                                        l_scs );

          // check if this is a new stencil
          if( l_scs.size() > 1 && l_scs[0] < l_scs[1] ) {
            EDGE_CHECK_LT( l_nSt, o_su.col[l_fa].TL_N_ST );
            // store the stencil
            o_su.col[l_fa].sc[l_nSt][0] = l_scs[0];
            o_su.col[l_fa].sc[l_nSt][1] = l_scs[1];
            l_nSt++;
          }

        }
      }
    }
};

/**
 * Application of the super-cell stencil.
 *
 * @paramt TL_T_EL element type.
 * @paramt TL_O_SP spatial order of the DG-scheme.
 * @paramt TL_N_QTS number of quantities..
 * @paramt TL_N_CRS number of fused simulations.
 **/
template< t_entityType   TL_T_EL,
          unsigned short TL_O_SP,
          unsigned short TL_N_QTS,
          unsigned short TL_N_CRS >
class edge::sc::ibnd::SuperCell {
  private:
    //! number of sub-cells per element
    static unsigned short const TL_N_SCS  = CE_N_SUB_CELLS( TL_T_EL, TL_O_SP );

  public:
    /**
     * Applies the face stencil.
     *
     * @param i_co stencil collection for the face.
     * @param io_dofsSc sub-cell DOFs to which the stencils are applied.
     *
     * @paramt TL_T_REAL floating point type.
     **/
    template< typename TL_T_REAL >
    static void applyFa( t_stCol< TL_T_EL,
                                  TL_O_SP > const  &i_co,
                         TL_T_REAL                  io_dofsSc[TL_N_QTS][TL_N_SCS][TL_N_CRS] ) {
     // only apply for tets
     if( TL_T_EL == TET4 ) {
       // iterate of stencil
       for( unsigned short l_st = 0; l_st < i_co.TL_N_ST; l_st++ ) {
         for( unsigned short l_qt = 0; l_qt < TL_N_QTS; l_qt++ ) {
           // gather average
           TL_T_REAL l_ave[TL_N_CRS];

           for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {
             l_ave[l_cr]  = TL_T_REAL(0.5) * io_dofsSc[l_qt][ i_co.sc[l_st][0] ][l_cr];
             l_ave[l_cr] += TL_T_REAL(0.5) * io_dofsSc[l_qt][ i_co.sc[l_st][1] ][l_cr];
           }

           // scatter average
           for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {
              io_dofsSc[l_qt][ i_co.sc[l_st][0] ][l_cr] = l_ave[l_cr];
              io_dofsSc[l_qt][ i_co.sc[l_st][1] ][l_cr] = l_ave[l_cr];
           }
         }
       }
     }
    }
};

#endif