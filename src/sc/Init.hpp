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
 * Initializes sub-cell data.
 **/
#ifndef SC_INIT_HPP
#define SC_INIT_HPP

#include "constants.hpp"
#include "io/logging.h"

/**
 * EDGEpre's preprocessed data.
 **/
namespace edge {
  namespace pre {
    namespace sc {
      //! sub-vertices adjacent to sub-cells
      extern const unsigned short *g_scsvRaw;
      extern const std::size_t     g_scsvSize;

      //! sub-cells adjacent to sub-cells (sub-faces as bridge)
      extern const unsigned short *g_scsfscRaw;
      extern const std::size_t     g_scsfscSize;

      //! types of the sub-cells' sub-faces
      extern const unsigned short *g_sctysfRaw;
      extern const std::size_t     g_sctysfSize;

      //! reference coordinates of the sub-vertices
      extern const double      *g_svcrdsRaw;
      extern const std::size_t  g_svcrdsSize;

      //! scatter matrix
      extern const double      *g_scatterRaw;
      extern const std::size_t  g_scatterSize;

      //! scatter matrix for sub-cells at the DG-element's faces
      extern const double      *g_scattersurfRaw;
      extern const std::size_t  g_scattersurfSize;

      //! gather matrix
      extern const double      *g_gatherRaw;
      extern const std::size_t  g_gatherSize;

      //! surface integration matrix for the sub-faces at the DG-surface
      extern const double      *g_sfintRaw;
      extern const std::size_t  g_sfintSize;
    }
  }
}

namespace edge {
  namespace sc {
    template< t_entityType   TL_T_EL,
              unsigned short TL_O_SP >
    class Init;
  }
}

/**
 * Subcell resolution in EDGE.
 *
 * @paramt TL_T_EL element type.
 * @paramt TL_O_SP order of the used solver.
 **/
template< t_entityType   TL_T_EL,
          unsigned short TL_O_SP >
class edge::sc::Init {
  private:
    //! number of dimensions
    static unsigned short const TL_N_DIS = C_ENT[ TL_T_EL ].N_DIM;

    //! number of vertices per element / sub-cell
    static unsigned short const TL_N_VES = C_ENT[ TL_T_EL ].N_VERTICES;

    //! number of faces per element / sub-cell
    static unsigned short const TL_N_FAS = C_ENT[ TL_T_EL ].N_FACES;

    //! number of sub-vertices
    static unsigned short const TL_N_SVS = CE_N_SUB_VERTICES( TL_T_EL, TL_O_SP );

    //! number of sub-faces
    static unsigned short const TL_N_SFS = CE_N_SUB_FACES( TL_T_EL, TL_O_SP );

    //! number of sub-cells
    static unsigned short const TL_N_SCS = CE_N_SUB_CELLS( TL_T_EL, TL_O_SP );

    //! number of modes
    static unsigned short const TL_N_MDS = CE_N_ELEMENT_MODES( TL_T_EL, TL_O_SP );

    /**
     * Gets the send sub-cells adjacent to DG-faces through sub-grid connectivity.
     *
     * @param i_scSfSc sub-cells adjacent to sub-cells (sub-faces as bridge).
     * @param o_faSfSc will be set to sub-cells adjacent to DG-faces (sub-faces as "bridge").
     **/
    template< typename TL_T_LID >
    static void getFaSfSc( TL_T_LID const i_scSfSc[ TL_N_SCS + TL_N_FAS * TL_N_SFS ][ TL_N_FAS ],
                           TL_T_LID       o_faSfSc[TL_N_FAS][TL_N_SFS] ) {
      // check for type of raw scSfSc
      static_assert(
        std::is_same< decltype( pre::sc::g_scsfscRaw ),
                      const unsigned short* >::value,
        "g_scsfscRaw is assumed to be const unsigned short*" );

      // iterate over DG surface
      for( unsigned short l_f1 = 0; l_f1 < TL_N_FAS; l_f1++ ) {
        // iterate over sub-faces of the DG-face
        for( unsigned short l_sf = 0; l_sf < TL_N_SFS; l_sf++ ) {
          // determine sub-cell
          TL_T_LID l_sc = i_scSfSc[ TL_N_SCS + l_f1 * TL_N_SFS + l_sf ][0];
          o_faSfSc[l_f1][l_sf] = l_sc;
          // check for a valid sub-cell id
          EDGE_CHECK_NE( o_faSfSc[l_f1][l_sf],
                         std::numeric_limits< unsigned short >::max() );

          // check that the rest points to nirvana
          for( unsigned short l_f2 = 1; l_f2 < TL_N_FAS; l_f2++ ) {
            EDGE_CHECK_EQ( i_scSfSc[ TL_N_SCS + l_f1 * TL_N_SFS + l_sf ][l_f2],
                           std::numeric_limits< unsigned short >::max() );
          }
        }
      }
    }

  public:
    /**
     * Initializes the sub-grid connectivity info.
     *
     * @param o_conn will be set to connectivity info.
     * @paramt integral type for local ids.
     **/
    template< typename TL_T_LID >
    static void connect( t_connect< TL_T_LID,
                                    TL_T_EL,
                                    TL_O_SP  > &o_conn ) {
      // size checks
      std::size_t l_size;

      // total number of subcells (including receive sub-cells)
      unsigned short l_nScs = TL_N_SCS + TL_N_FAS * TL_N_SFS;

      // raw pointer
      unsigned short const * l_ptr;

      /*
       * scSv
       */
      // check sizes
      l_size  = l_nScs;
      l_size *= TL_N_VES;
      EDGE_CHECK_EQ( edge::pre::sc::g_scsvSize,
                     l_size );

      // assign info
      l_ptr = edge::pre::sc::g_scsvRaw;

      for( unsigned short l_sc = 0; l_sc < l_nScs; l_sc++ ) {
        for( unsigned short l_ve = 0; l_ve < TL_N_VES; l_ve++ ) {
          o_conn.scSv[l_sc][l_ve] = *l_ptr;
          l_ptr++;
        }
      }

      /*
       * scSfSc
       */
      // check sizes
      l_size = l_nScs;
      l_size *= TL_N_FAS;
      EDGE_CHECK_EQ( edge::pre::sc::g_scsfscSize,
                     l_size );

      // assign info
      l_ptr = edge::pre::sc::g_scsfscRaw;

      for( unsigned short l_sc = 0; l_sc < l_nScs; l_sc++ ) {
        for( unsigned short l_fa = 0; l_fa < TL_N_FAS; l_fa++ ) {
          o_conn.scSfSc[l_sc][l_fa] = *l_ptr;
          l_ptr++;
        }
      }

      /*
       * faSfSc
       */
      getFaSfSc( o_conn.scSfSc,
                 o_conn.faSfSc );

      /*
       * scTySf
       */
      // check sizes
      l_size =  TL_N_SCS;
      l_size *= TL_N_FAS;
      EDGE_CHECK_EQ( edge::pre::sc::g_sctysfSize,
                     l_size );

      // assign info
      l_ptr = edge::pre::sc::g_sctysfRaw;

      for( unsigned short l_sc = 0; l_sc < TL_N_SCS; l_sc++ ) {
        for( unsigned short l_fa = 0; l_fa < TL_N_FAS; l_fa++ ) {
          o_conn.scTySf[l_sc][l_fa] = *l_ptr;
          l_ptr++;
        }
      }
    }

    /**
     * Initializes the sub-vertex characteristics.
     *
     * @param o_svChars will be set to connectivity info.
     * @paramt TL_T_REAL floating point tupe.
     **/
    template < typename TL_T_REAL >
    static void svChars( t_svChars< TL_T_REAL,
                                    TL_T_EL > o_svChars[TL_N_SVS] ) {
      // size checks
      std::size_t l_size;

      // raw pointer
      double const * l_ptr;

      /*
       * coords
       */
      // check sizes
      l_size = TL_N_SVS * TL_N_DIS;
      EDGE_CHECK_EQ( edge::pre::sc::g_svcrdsSize,
                     l_size );

      // assign info
      l_ptr = edge::pre::sc::g_svcrdsRaw;

      for( unsigned short l_sv = 0; l_sv < TL_N_SVS; l_sv++ ) {
        for( unsigned short l_di = 0; l_di < TL_N_DIS; l_di++ ) {
          o_svChars[l_sv].coords[l_di] = *l_ptr;
          l_ptr++;
        }
      }
    }

    /**
     * Initializes the operators of the sub-cell limiter.
     **/
    template < typename TL_T_REAL >
    static void ops( t_ops< TL_T_REAL,
                            TL_T_EL,
                            TL_O_SP > &o_ops ) {
      // size checks
      std::size_t l_size;

      // raw pointer
      double const * l_ptr;

      /**
       * scatter
       **/
      // check sizes
      l_size = TL_N_MDS;
      l_size *= TL_N_SCS;
      EDGE_CHECK_EQ( edge::pre::sc::g_scatterSize,
                     l_size );

      // assign info
      l_ptr = edge::pre::sc::g_scatterRaw;

      for( unsigned short l_md = 0; l_md < TL_N_MDS; l_md++ ) {
        for( unsigned short l_sc = 0; l_sc < TL_N_SCS; l_sc++ ) {
          o_ops.scatter[l_md][l_sc] = *l_ptr;
          l_ptr++;
        }
      }


      /**
       * scatter surf
       **/
      l_size  = TL_N_MDS;
      l_size *= TL_N_FAS*2;
      l_size *= TL_N_SFS;
      EDGE_CHECK_EQ( edge::pre::sc::g_scattersurfSize,
                     l_size );

      // assign info
      l_ptr = edge::pre::sc::g_scattersurfRaw;

      for( unsigned short l_fa = 0; l_fa < TL_N_FAS*2; l_fa++ ) {
        for( unsigned short l_md = 0; l_md < TL_N_MDS; l_md++ ) {
          for( unsigned short l_sf = 0; l_sf < TL_N_SFS; l_sf++ ) {
            o_ops.scatterSurf[l_fa][l_md][l_sf] = *l_ptr;
            l_ptr++;
          }
        }
      }

      /**
       * gather
       **/
       // check sizes
       l_size = TL_N_SCS;
       l_size *= TL_N_MDS;
       EDGE_CHECK_EQ( edge::pre::sc::g_gatherSize,
                      l_size );

       // assign info
       l_ptr = edge::pre::sc::g_gatherRaw;

       for( unsigned short l_sc = 0; l_sc < TL_N_SCS; l_sc++ ) {
         for( unsigned short l_md = 0; l_md < TL_N_MDS; l_md++ ) {
           o_ops.gather[l_sc][l_md] = *l_ptr;
           l_ptr++;
         }
       }

       /**
        * sfint
        **/
        // check size
        l_size = TL_N_FAS * TL_N_SFS;
        l_size *= TL_N_MDS;
        EDGE_CHECK_EQ( edge::pre::sc::g_sfintSize,
                       l_size );

        // assign info
        l_ptr = edge::pre::sc::g_sfintRaw;

        for( unsigned short l_fa = 0; l_fa < TL_N_FAS; l_fa++ ) {
          for( unsigned short l_sf = 0; l_sf < TL_N_SFS; l_sf++ ) {
            for( unsigned short l_md = 0; l_md < TL_N_MDS; l_md++ ) {
              o_ops.sfInt[l_fa][l_sf][l_md] = *l_ptr;
              l_ptr++;
            }
          }
        }
    }
};

#endif
