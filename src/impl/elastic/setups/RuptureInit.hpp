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
 * Initialization of rupture physics.
 **/

#ifndef EDGE_SEISMIC_RUPTURE_INIT_HPP
#define EDGE_SEISMIC_RUPTURE_INIT_HPP

#include "io/logging.h"
#include "linalg/Geom.hpp"
#include "data/Dynamic.h"

namespace edge {
  namespace elastic {
    namespace setups {
      template< t_entityType   TL_T_EL,
                unsigned short TL_O_SP,
                unsigned short TL_N_CRS >
      class RuptureInit;
    }
  }
}

/**
 * Initialization of rupture physics.
 *
 * @paramt TL_T_EL element type (e.g. TET4) for which the rupture physics are initialized.
 * @paramt TL_T_O_SP convergence rate in space.
 * @paramt TL_N_CRS number of fused runs.
 **/
template< t_entityType   TL_T_EL,
          unsigned short TL_O_SP,
          unsigned short TL_N_CRS >
class edge::elastic::setups::RuptureInit {
  private:
    //! number of dimensions
    static unsigned short const TL_N_DIS = C_ENT[TL_T_EL].N_DIM;

    //! entity type of the faces
    static t_entityType const TL_T_FA = C_ENT[TL_T_EL].TYPE_FACES;

    //! #vertices of the faces
    static unsigned short const TL_N_VES_FA = C_ENT[TL_T_FA].N_VERTICES;

    //! #vertices of the elements
    static unsigned short const TL_N_VES_EL = C_ENT[TL_T_EL].N_VERTICES;

    //! #faces of the elements
    static unsigned short const TL_N_FAS_EL = C_ENT[TL_T_EL].N_FACES;

    //! number of sub-faces per DG-face
    static unsigned short const TL_N_SFS = CE_N_SUB_FACES( TL_T_EL, TL_O_SP );

    // dummy fault data
    typedef struct { double sn0[TL_N_SFS][TL_N_CRS];
                     double ss0[TL_N_DIS-1][TL_N_SFS][TL_N_CRS]; } t_faultDataDummy;

  public:
    /**
     * Allocates memory for linear slip weakening.
     *
     * @param i_nRf number of rupture faces.
     * @param io_dynMem dynamic memory allocator.
     * @param o_lsw linear slip weakening data.
     *
     * @paramt TL_T_LID integral type for local ids.
     * @paramt TL_T_REAL floating point type.
     **/
    template< typename TL_T_LID,
              typename TL_T_REAL >
    static void alloc( TL_T_LID                              i_nRf,
                       data::Dynamic                        &io_dynMem,
                       solvers::t_LinSlipWeak< TL_T_REAL,
                                               TL_T_EL,
                                               TL_O_SP,
                                               TL_N_CRS > &o_lsw ) {
      // allocate face-local data
      std::size_t l_faSize  = i_nRf;
                  l_faSize *= sizeof( solvers::t_LinSlipWeakFace< TL_T_REAL > );
      o_lsw.fa = (solvers::t_LinSlipWeakFace< TL_T_REAL > *) io_dynMem.allocate( l_faSize );

      // allocate sub-face local data
      std::size_t l_sfSize  = std::size_t(i_nRf) * CE_N_SUB_FACES(TL_T_EL, TL_O_SP);
                  l_sfSize *= sizeof( solvers::t_LinSlipWeakSubFace< TL_T_REAL,
                                                                     TL_N_DIS,
                                                                     TL_N_CRS > );
      o_lsw.sf = ( solvers::t_LinSlipWeakSubFace< TL_T_REAL,
                                                  TL_N_DIS,
                                                  TL_N_CRS > (*) [CE_N_SUB_FACES(TL_T_EL, TL_O_SP)] ) io_dynMem.allocate( l_sfSize );
    }
    /////////////////////////////////////////////////
    // TODO: Add support for per-sub-face-point setup! //
    /////////////////////////////////////////////////

    /**
     * Initializes rupture physics using the linear slip weakening law.
     *
     * @param i_nFaces number of dense faces.
     * @param i_spType sparse type of the rupture faces.
     * @param i_scDgAd sub-cell reordering based on DG connectivity.
     * @param i_faVe vertices adjacent to the dense faces.
     * @param i_faEl elements adjacent to the dense faces.
     * @param i_elVe vertices adjacent to the elements.
     * @param i_elFa faces adjacent to the faces.
     * @param i_veChars vertex characteristics.
     * @param i_faChars face characteristics.
     * @param i_bgPars background parameters of the dense elements.
     * @param i_faultCrdSys fault coordinate systems. [0]: vector pointing from the plus to the minus side of the fault, [1]: along-strike direction, [2] (if 3D): along-dip direction
     * @param i_lswPars parameters of the linear slip weakening law.
     * @param i_doms domains where initial stress values at the fault are defined.
     * @param i_stressInit initial normal and shear stress values in the domains.
     * @param o_lsw linear slip weakening data.
     * @param i_faultData input fault data. Only used if not nullptr. If both domains and fault data are given, domains have priority.
     *
     * @paramt TL_T_LID integer type of per-rank local ids.
     * @paramt TL_T_INT_SP integer type of the sparse type.
     * @paramt TL_T_VE_CHARS struct summarizing the vertex characteristics.
     * @paramt TL_T_FA_CHARS struct summarizing the face characteristics.
     * @paramt TL_T_BG_PARS struct summarizing the background parameters.
     * @paramt TL_T_REAL_MESH type used in floating point arithmetic for mesh-related data.
     * @paramt TL_T_REAL_COMP type used in floating point arithmetic for computational data.
     * @paramt TL_T_DOM class implementing the concept of a domain.
     * @paramt TL_T_FAULT_DATA type of the fault data, offering members .ss0 for shear stresses and .sn0 for normal stress.
     **/
    template< typename TL_T_LID,
              typename TL_T_INT_SP,
              typename TL_T_VE_CHARS,
              typename TL_T_FA_CHARS,
              typename TL_T_BG_PARS,
              typename TL_T_REAL_MESH,
              typename TL_T_REAL_COMP,
              typename TL_T_DOM,
              typename TL_T_FAULT_DATA = t_faultDataDummy >
    static void linSlipWeak( TL_T_LID                                              const    i_nFaces,
                             TL_T_INT_SP                                           const    i_spType,
                             unsigned short                                        const    i_scDgAd[TL_N_VES_FA][TL_N_SFS],
                             TL_T_LID                                              const (* i_faVe)[TL_N_VES_FA],
                             TL_T_LID                                              const (* i_faEl)[2],
                             TL_T_LID                                              const (* i_elVe)[TL_N_VES_EL],
                             TL_T_LID                                              const (* i_elFa)[TL_N_FAS_EL],
                             TL_T_VE_CHARS                                         const  * i_veChars,
                             TL_T_FA_CHARS                                         const  * i_faChars,
                             TL_T_BG_PARS                                          const  * i_bgPars,
                             TL_T_REAL_MESH                                        const    i_faultCrdSys[TL_N_DIS][TL_N_DIS],
                             TL_T_REAL_COMP                                        const    i_lswPars[TL_N_CRS][3],
                             std::vector< TL_T_DOM >                               const    i_doms[TL_N_CRS],
                             std::vector< std::array< TL_T_REAL_COMP, TL_N_DIS > > const    i_stressInit[TL_N_CRS],
                             solvers::t_LinSlipWeak< TL_T_REAL_COMP,
                                                     TL_T_EL,
                                                     TL_O_SP,
                                                     TL_N_CRS >                           & o_lsw,
                             TL_T_FAULT_DATA                                        const * i_faultData = nullptr ) {
      // iterate over runs and setup friction laws
      for( unsigned short l_ru = 0; l_ru < TL_N_CRS; l_ru++ ) {
        EDGE_CHECK( i_lswPars[l_ru][0] > TOL.SOLVER );
        o_lsw.gl.mus[l_ru] = i_lswPars[l_ru][0];
        EDGE_CHECK( i_lswPars[l_ru][1] > TOL.SOLVER );
        o_lsw.gl.mud[l_ru] = i_lswPars[l_ru][1];
        EDGE_CHECK( i_lswPars[l_ru][2] > TOL.SOLVER );
        o_lsw.gl.dcInv[l_ru] = (TL_T_REAL_COMP) 1 / i_lswPars[l_ru][2];
      } 

      // iterate over dense faces
      TL_T_LID l_spId = 0;

      for( TL_T_LID l_fa = 0; l_fa < i_nFaces; l_fa++ ) {
        // ignore all non-rupture faces
        if( (i_faChars[l_fa].spType & i_spType) != i_spType ) continue;

        // adjacent elements id
        TL_T_LID l_elL = i_faEl[l_fa][0];
        TL_T_LID l_elR = i_faEl[l_fa][1];

        // set up face data
        TL_T_REAL_MESH l_sProd = (TL_N_DIS == 2) ? linalg::Geom::sprod2( i_faultCrdSys[0], i_faChars[l_fa].outNormal ) :
                                                   linalg::Geom::sprod3( i_faultCrdSys[0], i_faChars[l_fa].outNormal );
        EDGE_CHECK( std::abs(l_sProd) > TOL.MESH );

        // set minus-plus distinction
        o_lsw.fa[l_spId].lEqM = (l_sProd < 0);

        // compute shear modulus divided by Lame parmeter mu
        TL_T_REAL_COMP l_csDivMuL  = std::sqrt( (TL_T_REAL_COMP) 1 / (i_bgPars[l_elL].mu * i_bgPars[l_elL].rho) );
        TL_T_REAL_COMP l_csDivMuR  = std::sqrt( (TL_T_REAL_COMP) 1 / (i_bgPars[l_elR].mu * i_bgPars[l_elR].rho) );
        o_lsw.fa[l_spId].csDmuM = o_lsw.fa[l_spId].lEqM  ? l_csDivMuL : l_csDivMuR;
        o_lsw.fa[l_spId].csDmuP = o_lsw.fa[l_spId].lEqM  ? l_csDivMuR : l_csDivMuL;

        // derive pt of the face (might be outside for degenerated elements)
        // TODO: use quad points for derivation of location
        TL_T_REAL_MESH l_pt[TL_N_DIS];
        for( unsigned short l_di = 0; l_di < TL_N_DIS; l_di++ ) {
          l_pt[l_di] = 0;
          for( unsigned short l_ve = 0; l_ve < TL_N_VES_FA; l_ve++ ) {
            l_pt[l_di] += i_veChars[ i_faVe[l_fa][l_ve] ].coords[l_di];
          }
          l_pt[l_di] /= TL_N_VES_FA;
        }

        // determine the local rupture face-id of the left element
        unsigned short l_faL = std::numeric_limits< unsigned short >::max();
        for( unsigned short l_fl = 0; l_fl < TL_N_FAS_EL; l_fl++ ) {
          if( i_elFa[l_elL][l_fl] == l_fa ) l_faL = l_fl;
        }
        EDGE_CHECK_NE( l_faL, std::numeric_limits< unsigned short >::max() );

        // left element's local vertex ids, associated with the faces
        unsigned short const *l_elVeFa = C_REF_ELEMENT.FA_VE_CC.ENT[TL_T_EL];
        l_elVeFa += l_faL * TL_N_VES_FA;

        // check that the first element-face vertex matches the first face vertex for 3D elements
        EDGE_CHECK( (TL_N_DIS != 3) ||
                    (i_elVe[l_elL][ l_elVeFa[0] ] == i_faVe[l_fa][0]) );

        // determine if we have reorder the sub-faces
        bool l_reorder = false;
        if( TL_N_DIS == 2) {
          // 2D elements are reorderd if the first vertex has large id than the second one
          l_reorder = (i_elVe[l_elL][l_elVeFa[1]] < i_elVe[l_elL][l_elVeFa[0]]);
        }
        else {
          for( unsigned short l_v0 = 2; l_v0 < TL_N_VES_FA; l_v0++ ) {
            // 3D elements are reorder if the second id is not the second smallest one (by construction 0 holds the smallest id)
            if( i_elVe[l_elL][l_elVeFa[1]] > i_elVe[l_elL][l_elVeFa[l_v0]] )
              l_reorder = true;
          }
        }

        // iterate over fused runs
        for( unsigned short l_ru = 0; l_ru < TL_N_CRS; l_ru++ ) {
          // iterate over sub-faces
          for( unsigned short l_sf = 0; l_sf < TL_N_SFS; l_sf++ ) {
            // set initial values invalid
            o_lsw.sf[l_spId][l_sf].sn0[l_ru] = std::numeric_limits< TL_T_REAL_COMP >::max();
            o_lsw.sf[l_spId][l_sf].co0[l_ru] = 0;
            o_lsw.sf[l_spId][l_sf].st[l_ru]  = 0;
            o_lsw.sf[l_spId][l_sf].ss0A[l_ru] = 0;
            for( unsigned short l_di = 0; l_di < TL_N_DIS-1; l_di++ ) {
              o_lsw.sf[l_spId][l_sf].ss0[l_di][l_ru] = std::numeric_limits< TL_T_REAL_COMP >::max();
            }
            o_lsw.sf[l_spId][l_sf].muf[l_ru] = std::numeric_limits< TL_T_REAL_COMP >::max();
            for( unsigned short l_di = 0; l_di < TL_N_DIS-1; l_di++ ) {
              o_lsw.sf[l_spId][l_sf].dd[l_di][l_ru] = 0;
              o_lsw.sf[l_spId][l_sf].sr[l_di][l_ru] = 0;
              o_lsw.sf[l_spId][l_sf].tr[l_di][l_ru] = 0;
            }

            // data-based input
            if( i_faultData != nullptr ) {
              unsigned short l_sfRe = l_reorder ? i_scDgAd[0][l_sf] : l_sf;
              // assign values
              o_lsw.sf[l_spId][l_sf].sn0[l_ru] = i_faultData[l_spId].sn0[l_sfRe][l_ru];
              o_lsw.sf[l_spId][l_sf].co0[l_ru] = i_faultData[l_spId].co0[l_sfRe][l_ru];
              for( unsigned short l_di = 0; l_di < TL_N_DIS-1; l_di++ )
                o_lsw.sf[l_spId][l_sf].ss0[l_di][l_ru] = i_faultData[l_spId].ss0[l_di][l_sfRe][l_ru];
            }

            // config-based input iterate over domains
            for( std::size_t l_do = 0; l_do < i_doms[l_ru].size(); l_do++ ) {
              // check if the point is inside
              if( i_doms[l_ru][l_do].inside( l_pt ) ) {
                // set intial stress values
                o_lsw.sf[l_spId][l_sf].sn0[l_ru] = i_stressInit[l_ru][l_do][0];
                for( unsigned short l_di = 0; l_di < TL_N_DIS-1; l_di++ ) {
                  o_lsw.sf[l_spId][l_sf].ss0[0+l_di][l_ru] = i_stressInit[l_ru][l_do][1+l_di];
                }
                // init friction coefficient with static friction coefficient
                o_lsw.sf[l_spId][l_sf].muf[l_ru] = i_lswPars[l_ru][0];

                // abort after the first found region
                break;
              }
            }
          }
        }
        // increase rupture-face counter
        l_spId++;
      }

      // assemble names of quantities
      std::string l_str;
      for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {
        l_str = "sn0_" + std::to_string(l_cr);
        strcpy( o_lsw.sfQtNa[0*TL_N_CRS + l_cr], l_str.c_str() );

        l_str = "co0_" + std::to_string(l_cr);
        strcpy( o_lsw.sfQtNa[1*TL_N_CRS + l_cr], l_str.c_str() );

        l_str = "ss0A_" + std::to_string(l_cr);
        strcpy( o_lsw.sfQtNa[(2+TL_N_DIS-1)*TL_N_CRS + l_cr], l_str.c_str() );

        l_str = "muf_"  + std::to_string(l_cr);
        strcpy( o_lsw.sfQtNa[(2+TL_N_DIS-1+1)*TL_N_CRS + l_cr], l_str.c_str() );

        l_str = "st_"  + std::to_string(l_cr);
        strcpy( o_lsw.sfQtNa[(2+TL_N_DIS-1+2)*TL_N_CRS + l_cr], l_str.c_str() );

        for( unsigned short l_di = 0; l_di < TL_N_DIS-1; l_di++ ) {
          l_str = "ss0_" + std::to_string(l_di) + "_" + std::to_string(l_cr);;
          strcpy( o_lsw.sfQtNa[(2)*TL_N_CRS + l_di*TL_N_CRS + l_cr], l_str.c_str() );

          l_str = "tr_"  + std::to_string(l_di) + "_" + std::to_string(l_cr);
          strcpy( o_lsw.sfQtNa[(2+TL_N_DIS-1+3)*TL_N_CRS + l_di*TL_N_CRS + l_cr], l_str.c_str() );

          l_str = "sr_"  + std::to_string(l_di) + "_" + std::to_string(l_cr);
          strcpy( o_lsw.sfQtNa[(2+TL_N_DIS-1+3+(TL_N_DIS-1)  )*TL_N_CRS + l_di*TL_N_CRS + l_cr], l_str.c_str() );

          l_str ="dd_"  + std::to_string(l_di) + "_" + std::to_string(l_cr);
          strcpy( o_lsw.sfQtNa[(2+TL_N_DIS-1+3+2*(TL_N_DIS-1))*TL_N_CRS + l_di*TL_N_CRS + l_cr], l_str.c_str() );
        }
      }
      // set pointers
      for( unsigned short l_en = 0; l_en < 5*TL_N_CRS + 4 * (TL_N_DIS-1) * TL_N_CRS; l_en++ ) {
        o_lsw.sfQtNaPtr[l_en] = o_lsw.sfQtNa[l_en];
      }
    }
};
#endif
