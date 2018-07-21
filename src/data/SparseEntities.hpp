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
 * Support for sparse entities.
 **/

#ifndef EDGE_DATA_SPARSE_ENTITIES_HPP
#define EDGE_DATA_SPARSE_ENTITIES_HPP

#include <limits>
#include "parallel/Mpi.h"
#include "io/logging.h"

#include "EntityLayout.h"
namespace edge {
  namespace data {
    class SparseEntities;
  }
}

/**
 * Support for sparse entities.
 **/
class edge::data::SparseEntities {
  private:
    /**
     * Initializes a partial sparse layout by 1) setting entity related sizes to zero,
     * 2) copying dimensions of the data structures, and 3) copying MPI-neighboring ranks and time groups.
     *
     * @param i_deLayout dense layout.
     * @param o_spLayout will be initialized based on the dense layout.
     *
     * @paramt TL_T_LAYOUT type of the layouts.
     **/
    template< typename TL_T_LAYOUT >
    static void initPartLayout( TL_T_LAYOUT const & i_deLayout,
                                TL_T_LAYOUT       & o_spLayout ) {
      o_spLayout.timeGroups.resize( i_deLayout.timeGroups.size() );
      for( std::size_t l_tg = 0; l_tg < o_spLayout.timeGroups.size(); l_tg++ ) {
        o_spLayout.timeGroups[l_tg].inner.size  = 0;

        // check matching size of MPI-data
        EDGE_CHECK_EQ( i_deLayout.timeGroups[l_tg].neRanks.size(),
                       i_deLayout.timeGroups[l_tg].send.size() );
        EDGE_CHECK_EQ( i_deLayout.timeGroups[l_tg].neRanks.size(),
                       i_deLayout.timeGroups[l_tg].receive.size() );
        EDGE_CHECK_EQ( i_deLayout.timeGroups[l_tg].neRanks.size(),
                       i_deLayout.timeGroups[l_tg].neTgs.size() );

        // resize
        o_spLayout.timeGroups[l_tg].send.resize(    i_deLayout.timeGroups[l_tg].send.size()    );
        o_spLayout.timeGroups[l_tg].receive.resize( i_deLayout.timeGroups[l_tg].receive.size() );

        o_spLayout.timeGroups[l_tg].neRanks.resize( i_deLayout.timeGroups[l_tg].neRanks.size() );
        o_spLayout.timeGroups[l_tg].neTgs.resize(   i_deLayout.timeGroups[l_tg].neTgs.size() );

        //init MPI-data
        for( unsigned int l_nr = 0; l_nr < o_spLayout.timeGroups[l_tg].send.size(); l_nr++ ) {
          o_spLayout.timeGroups[l_tg].send[l_nr].size     = 0;
          o_spLayout.timeGroups[l_tg].receive[l_nr].size  = 0;

          // copy shared data to sparse layout
          o_spLayout.timeGroups[l_tg].neRanks[l_nr] = i_deLayout.timeGroups[l_tg].neRanks[l_nr];
          o_spLayout.timeGroups[l_tg].neTgs[l_nr]   = i_deLayout.timeGroups[l_tg].neTgs[l_nr];
        }
      }
    }

  public:
    /**
     * Determines the number of sparse entities for the given dense entities.
     *
     * @param i_nen number of dense entities.
     * @param i_spType value of the sparse type, which is used (as a bit mask) to count the sparse entities.
     * @param i_chars entity characteristics corresponding to the dense layout (includes the sparse type).
     * @return number of sparse entities of the given type in the dense layout.
     *
     * @paramt TL_T_INT_LID integer for local ids.
     * @paramt TL_T_INT_SP integer type of the sparse type.
     * @paramt TL_T_EN_CHARS struct of the entity characteristics.
     **/
    template< typename TL_T_INT_LID, typename TL_T_INT_SP, typename TL_T_EN_CHARS >
    static TL_T_INT_LID nSp( TL_T_INT_LID          i_nEn,
                             TL_T_INT_SP           i_spType,
                             TL_T_EN_CHARS const * i_chars ) {
      TL_T_INT_LID l_nSp = 0;

      for( TL_T_INT_LID l_en = 0; l_en < i_nEn; l_en++ )
        if( (i_chars[l_en].spType & i_spType) == i_spType ) l_nSp++;

      return l_nSp;
    }

    /**
     * Extracts a sparse layout based from a dense layout based on the given sparse type.
     *
     * Remark: The sparse type is handle as a bit flag.
     *         An entity in the dense layout is assumed to be in the sparse layout if all corresponding bits are set.
     *
     * @param i_spType value of the sparse type, which is used (as a bit mask) to derive entities of the sparse layout.
     * @param i_chars entity characteristics corresponding to the dense layout (includes the sparse type).
     * @param i_deLayout dense entity layout.
     * @param o_spLayout will be set to sparse entity layout.
     *
     * @paramt TL_T_INT_SP integer type of the sparse type.
     * @paramt TL_T_EN_CHARS struct of the entity characteristics.
     **/
    template <typename TL_T_INT_SP, typename TL_T_EN_CHARS>
    static void denseToSparse( TL_T_INT_SP               i_spType,
                               TL_T_EN_CHARS     const * i_chars,
                               t_enLayout        const & i_deLayout,
                               t_enLayout              & o_spLayout ) {
      // init the partial layout
      initPartLayout( i_deLayout, o_spLayout );

      // iterate over the dense layout and assign sizes of the groups
      for( std::size_t l_tg = 0; l_tg < i_deLayout.timeGroups.size(); l_tg++ ) {
        // iterate over dense inner entities
        for( int_el l_de = i_deLayout.timeGroups[l_tg].inner.first; l_de < i_deLayout.timeGroups[l_tg].inner.first + i_deLayout.timeGroups[l_tg].inner.size; l_de++ ) {
          if( (i_chars[l_de].spType & i_spType) == i_spType ) o_spLayout.timeGroups[l_tg].inner.size++;
        }

        // iter over dense send entities
        for( unsigned int l_nr = 0; l_nr < o_spLayout.timeGroups[l_tg].send.size(); l_nr++ ) {
          for( int_el l_de = i_deLayout.timeGroups[l_tg].send[l_nr].first; l_de < i_deLayout.timeGroups[l_tg].send[l_nr].first + i_deLayout.timeGroups[l_tg].send[l_nr].size; l_de++ ) {
            if( (i_chars[l_de].spType & i_spType) == i_spType ) o_spLayout.timeGroups[l_tg].send[l_nr].size++;
          }
        }

        // iter over dense receive entities
        for( unsigned int l_nr = 0; l_nr < o_spLayout.timeGroups[l_tg].receive.size(); l_nr++ ) {
          for( int_el l_de = i_deLayout.timeGroups[l_tg].receive[l_nr].first; l_de < i_deLayout.timeGroups[l_tg].receive[l_nr].first + i_deLayout.timeGroups[l_tg].receive[l_nr].size; l_de++ ) {
            if( (i_chars[l_de].spType & i_spType) == i_spType ) o_spLayout.timeGroups[l_tg].receive[l_nr].size++;
          }
        }
      }

      // complete the layout be setting remaining derivable quantities
      edge::data::EntityLayout::sizesToLayout( o_spLayout );
    }

    /**
     * Extracts a sparse layout for the entities based on adjacent entities with sparse types.
     *
     * Example:     *********  el0
     *             *  *     *
     *        el3 *   1*0 2 *
     *           *2   spf3 1*
     *          *  0      * *
     *         *****spf1****
     *       * *   2     *
     *    *0   *       * el2
     *  *     2*0   spf2
     *    *1   *  1*
     *  el1  * * *
     *         *
     *
     * Sparse adjacencies:
     *
     *   el0-0: spf3
     *   el0-1: ----
     *   el0-2: ----
     *
     *   el1-0: ----
     *   el1-1: ----
     *   el1-2: ----
     *
     *   el2-0: ----
     *   el2-1: spf2
     *   el2-2: spf1
     *
     *   el3-0: spf1
     *   el3-1: spf3
     *   el3-2: ----
     *
     * @param i_nAdjPerEn number entities adjacent to each of the dense entities.
     * @param i_enEn adjacent entities (nAdjPerEn) for each of the dense entities. Assumed is a flat array, meaning: [de*i_nAdjPerEn] gives the de's dense entity and [de*i_nAdjPerEn + ae] the ae's adjacent entity.
     * @param i_spType sparse type identitying the sparse entities through matching bits.
     * @param i_charsAdj characteristics of the adjecent entities (having a member .spType).
     * @param i_deLayout dense entity layout.
     * @param o_spLayout will be set to sparse layout of entities having one or more adjacent entities with the given sparse type.
     *
     * @paramt TL_T_INT_LID integer type of local ids.
     * @paramt TL_T_INT_SP integer type of the sparse type.
     * @paramt TL_T_ADJ_CHARS struct of the adjacent entities' characteristics. Offers a member .spType for comparison with i_spType.
     * @paramt TL_T_LAYOUT struct represnting the entity layout.
     **/
    template< typename TL_T_INT_LID,
              typename TL_T_INT_SP,
              typename TL_T_ADJ_CHARS,
              typename TL_T_LAYOUT >
    static void denseToSparseAdj( unsigned short           i_nAdjPerEn,
                                  TL_T_INT_LID     const * i_enEn,
                                  TL_T_INT_SP              i_spType,
                                  TL_T_ADJ_CHARS * const   i_charsAdj,
                                  TL_T_LAYOUT      const & i_deLayout,
                                  TL_T_LAYOUT            & o_spLayout ) {
      // init the partial layout
      initPartLayout( i_deLayout, o_spLayout );

      // iterate over the dense layout and assign sizes of the groups
      for( std::size_t l_tg = 0; l_tg < i_deLayout.timeGroups.size(); l_tg++ ) {
        // iterate over dense inner entities
        for( TL_T_INT_LID l_de = i_deLayout.timeGroups[l_tg].inner.first; l_de < i_deLayout.timeGroups[l_tg].inner.first + i_deLayout.timeGroups[l_tg].inner.size; l_de++ ) {
          // iterate over adjacent entities
          for( unsigned short l_ae = 0; l_ae < i_nAdjPerEn; l_ae++ ) {
            TL_T_INT_LID l_aeId = i_enEn[ l_de * i_nAdjPerEn + l_ae ];
            if( (i_charsAdj[l_aeId].spType & i_spType) == i_spType ) {
              o_spLayout.timeGroups[l_tg].inner.size++;
              break;
            }
          }
        }

#ifdef PP_USE_MPI
        // assemble data of dense send entities
        for( unsigned int l_nr = 0; l_nr < o_spLayout.timeGroups[l_tg].send.size(); l_nr++ ) {
          for( int_el l_de = i_deLayout.timeGroups[l_tg].send[l_nr].first; l_de < i_deLayout.timeGroups[l_tg].send[l_nr].first + i_deLayout.timeGroups[l_tg].send[l_nr].size; l_de++ ) {
            for( unsigned short l_ae = 0; l_ae < i_nAdjPerEn; l_ae++ ) {
              TL_T_INT_LID l_aeId = i_enEn[ l_de * i_nAdjPerEn + l_ae ];
              if( (i_charsAdj[l_aeId].spType & i_spType) == i_spType ) {
                o_spLayout.timeGroups[l_tg].send[l_nr].size++;
                break;
              }
            }
          }
        }

        // allocate send and receiver buffers
        TL_T_INT_LID l_nRgns = o_spLayout.timeGroups[l_tg].neRanks.size();
        if( l_nRgns == 0 ) continue; // continue with next time group for MPI with 1 rank

        std::vector< TL_T_INT_LID > l_buffSend( l_nRgns );
        std::vector< TL_T_INT_LID > l_buffRecv( l_nRgns );
        std::vector< MPI_Request  > l_reqSend( l_nRgns );
        std::vector< MPI_Request  > l_reqRecv( l_nRgns );

        // copy send data to buffer
        for( unsigned int l_nr = 0; l_nr < o_spLayout.timeGroups[l_tg].send.size(); l_nr++ ) {
          l_buffSend[l_nr] = o_spLayout.timeGroups[l_tg].send[l_nr].size;
        }

        // fix MPI-tags below as the default 0 might lead to colliding messages
        // if more than one message goes out to the same rank (multiple time groups)
        EDGE_CHECK_EQ( l_tg, 0 );

        // exchange the data
        parallel::Mpi::iSendTgRg( (unsigned char *) &l_buffSend[0],
                                  sizeof(TL_T_INT_LID),
                                  l_nRgns,
                                  &o_spLayout.timeGroups[l_tg].neRanks[0],
                                  &l_reqSend[0] );
        parallel::Mpi::iRecvTgRg( (unsigned char *) &l_buffRecv[0],
                                  sizeof(TL_T_INT_LID),
                                  l_nRgns,
                                  &o_spLayout.timeGroups[l_tg].neRanks[0],
                                  &l_reqRecv[0] );

        parallel::Mpi::waitAll( l_nRgns,
                                &l_reqSend[0] );
        parallel::Mpi::waitAll( l_nRgns,
                                &l_reqRecv[0] );

        // set receive layout
        for( unsigned int l_nr = 0; l_nr < o_spLayout.timeGroups[l_tg].receive.size(); l_nr++ ) {
          o_spLayout.timeGroups[l_tg].receive[l_nr].size = l_buffRecv[l_nr];
        }
#endif
      }

      // complete the layout be setting remaining derivable quantities
      edge::data::EntityLayout::sizesToLayout( o_spLayout );
    }

    /**
     * Determines the first sparse ids and sizes of a sequence of sub-regions.
     *
     * Example (ex1: example 1, ex2: example 2, sr*: subregion):
     *   dense ids:  [  0  1  2  3  4  5  6  7  8  9  10  11  12  13  14  15  16  17  18  19 ]
     *   sparse ids: [        0     1  2        3          4       5   6           7       8 ]
     *
     *   region ex1:                     [3          |4      |5      |6              ]
     *                                        sr0        sr1    sr2           sr3
     *
     *   region ex2:        [0                 |3                        |7              ]
     *                               sr0                   sr1                  sr2
     *
     * @param i_rgnFirst first dense entry of the region.
     * @param i_nSubRgns number of sub-regions.
     * @param i_subRgnSizes sizes of the sub-regions.
     * @param i_spType sparse type used for the bit comparisons.
     * @param i_deChars dense entity characteristics hosting the respective sparse type.
     * @param o_subRgnsFirstSp will be set to first sparse id of every subregion.
     * @param o_subRgnsSizeSp will be set to sparse size of every subregion.
     *
     * @paramt TL_T_INT_LID integer type of per-rank local ids.
     * @paramt TL_T_INT_SP integer type of the sparse type.
     * @paramt TL_T_DE_CHAR dense entity characteristics with an accessible class member .spType
     **/
    template< typename TL_T_INT_LID,
              typename TL_T_INT_SP,
              typename TL_T_DE_CHARS >
    static void subRgnsSpId( TL_T_INT_LID         i_rgnFirst,
                             unsigned int         i_nSubRgns,
                             TL_T_INT_LID const * i_subRgnSizes,
                             TL_T_INT_SP          i_spType,
                             TL_T_DE_CHARS        i_deChars,
                             TL_T_INT_LID       * o_subRgnsFirstSp,
                             TL_T_INT_LID       * o_subRgnsSizeSp ) {
      // determine the first sparse entity of the entire region
      TL_T_INT_LID l_first = 0;
      for( TL_T_INT_LID l_en = 0; l_en < i_rgnFirst; l_en++ ) {
        if( (i_deChars[l_en].spType & i_spType) == i_spType ) l_first++;
      }

      // update the first ids of the sub-regions
      TL_T_INT_LID l_en = i_rgnFirst;
      for( unsigned int l_sr = 0; l_sr < i_nSubRgns; l_sr++ ) {
        o_subRgnsFirstSp[l_sr] = l_first;

        TL_T_INT_LID l_up = l_en +  i_subRgnSizes[l_sr];
        o_subRgnsSizeSp[l_sr] = 0;

        for( ; l_en < l_up; l_en++ ) {
          if( (i_deChars[l_en].spType & i_spType) == i_spType ) {
            l_first++;
            o_subRgnsSizeSp[l_sr]++;
          }
        }
      }
    }

    /**
     * Inherits sparse information within the entities.
     *
     * Example:
     *   "parent" sparse type: 01
     *   "child"  sparse type: 10
     *
     *   Input:      Output:
     *
     *   en | bits   en | bits
     *   0  | 01     0  | 11
     *   1  | 11     1  | 11
     *   2  | 01     2  | 11
     *   3  | 10     3  | 10
     *   4  | 01     4  | 11
     *   5  | 00     5  | 00
     **/
    template< typename TL_T_INT_LID,
              typename TL_T_INT_SP,
              typename TL_T_CHARS >
    static void inherit( TL_T_INT_LID  i_nEn,
                         TL_T_INT_SP   i_spTypeIn,
                         TL_T_INT_SP   i_spTypeOut,
                         TL_T_CHARS   *io_chars ) {
      // iterate over dense entities
      for( TL_T_INT_LID l_en = 0; l_en < i_nEn; l_en++ ) {
        if( (io_chars[l_en].spType & i_spTypeIn) == i_spTypeIn ) {
          io_chars[l_en].spType |= i_spTypeOut;
        }
      }
    }

    /**
     * Syncs the entity characteristics, of duplicated entities.
     *
     * @param i_nEn number of entities.
     * @param i_enDaMe mapping of entity ids: data-to-mesh.
     * @param i_enMeDa mapping of element ids: mesh-to-data.
     *
     * @paramt TL_T_LID integral type of local ids.
     * @paramt TL_T_CHARS entity characteristics, offering member .spType.
     **/
    template< typename TL_T_LID,
              typename TL_T_CHARS >
    static void syncDupl( TL_T_LID          i_nEn,
                          TL_T_LID   const *i_enDaMe,
                          TL_T_LID   const *i_enMeDa,
                          TL_T_CHARS       *io_chars ) {
      // iterate over entities
#ifdef PP_USE_OMP
#pragma omp parallel for
#endif
      for( TL_T_LID l_en = 0; l_en < i_nEn; l_en++ ) {
        TL_T_LID l_enDo = i_enDaMe[ l_en   ];
                 l_enDo = i_enMeDa[ l_enDo ];
        io_chars[l_en].spType = io_chars[l_enDo].spType;
      }
    }

    /**
     * Propagates sparse information to adjacent entities.
     *
     * Example (i_spTypeIn == i_spTypeOut):
     *
     *   Input bits  Input bits  Adjacency of    Result (critical bit is
     *   of en0:     of en1:     entity 0 to 1:  the second (x): xy:
     *
     *   en0 | bits  en0 | bits  en0 | en1       en1
     *   0   | 00    0   | 00    0   | 0-2       0 | 00
     *   1   | 11    1   | 01    1   | 1-2       1 | 11
     *   2   | 00    2   | 00    2   | 3-1       2 | 10
     *   3   | 11    3   | 01    3   | 4-2       3 | 01
     *   4   | 11    4   | 01    4   | 1-4       4 | 11
     *   5   | 00                5   | 2-3
     *
     * @param i_nEn0 number of en0-entities.
     * @param i_nEn0PerEn1 number of en0-entities adjacent to a single en1-entity.
     * @param i_en0En1 adjacency information from en0 to en1. Assumed is a flat array, meaning: id0*i_nEn0PerEn1 + id1] gives the entity (en1) adjacent to id0 (en0) at position id1.
     * @param i_spTypeIn sparse type, which triggers propagation.
     * @param i_spTypeOut sparse type which is set in elements adjacent to a triggered element.
     * @param i_charsEn0 characteristics of entity type en0 (has to provide the member .spType).
     * @param i_charsEn1 characteristics of entity type en1 (has to provide the member .spType) which will be updated via bit-wise | if an adjacent en0 entity has the respective sparse type.
     * @param i_enDaMe mapping of entity ids: data-to-mesh.
     * @param i_enMeDa mapping of element ids: mesh-to-data.
     *
     * @paramt TL_T_LID integer type of local ids.
     * @paramt TL_T_INT_SP integer type of the sparse type.
     * @paramt TL_T_EN0_CHARS struct of the en0 entities' characteristics. Offers a member .spType for comparison with i_spType.
     * @paramt TL_T_EN1_CHARS struct of the en1 entities' characteristics. Offers a member .spType for comparison with i_spType.
     **/
    template< typename TL_T_LID,
              typename TL_T_INT_SP,
              typename TL_T_EN0_CHARS,
              typename TL_T_EN1_CHARS >
    static void propAdj( TL_T_LID               i_nEn0,
                         unsigned short         i_nEn0PerEn1,
                         TL_T_LID       const * i_en0En1,
                         TL_T_INT_SP            i_spTypeIn,
                         TL_T_INT_SP            i_spTypeOut,
                         TL_T_EN0_CHARS const * i_charsEn0,
                         TL_T_EN1_CHARS       * o_charsEn1,
                         TL_T_LID       const * i_enDaMe = nullptr,
                         TL_T_LID       const * i_enMeDa = nullptr ) {
      // maximum number of modified output entities
      TL_T_LID l_nEn1 = 0;

      // iterate over en0
      for( TL_T_LID l_en = 0; l_en < i_nEn0; l_en++ ) {
        // check for sparse type
        if( (i_charsEn0[l_en].spType & i_spTypeIn) == i_spTypeIn ) {
          // iterate over adjacent entities
          for( unsigned short l_ae = 0; l_ae < i_nEn0PerEn1; l_ae++ ) {
            TL_T_LID l_aeId = i_en0En1[ l_en*i_nEn0PerEn1 + l_ae ];

            // propagate the info (if entity exists)
            if( l_aeId != std::numeric_limits< TL_T_LID >::max() ) {
              l_nEn1 = std::max( l_nEn1, l_aeId+1 );
              o_charsEn1[l_aeId].spType |= i_spTypeOut;
            }
          }
        }
      }

      // synchronize output characteristics
      if( i_enDaMe != nullptr ) {
        EDGE_CHECK_NE( i_enMeDa, nullptr );
        syncDupl( l_nEn1,
                  i_enDaMe,
                  i_enMeDa,
                  o_charsEn1 );
      }
    }

    /**
     * Propagates sparse information to adjacent entities.
     *   Array of pointers (i_en0en1) equivalent of the single array function.
     *   Allows for dynamic number of adjacent entities, e.g. for vertives.
     *
     * @param i_nEn0 number of en0-entities.
     * @param i_en0En1 adjacency information from en0 to en1. For entity en0_n, the expression i_en0En1[en0_n+1] - i_en0En1[en0_n] is assumed to give the number of adjacent en1 entities.
     * @param i_spTypeIn sparse type, which triggers propagation.
     * @param i_spTypeOut sparse type which is set in elements adjacent to a triggered element.
     * @param i_charsEn0 characteristics of entity type en0 (has to provide the member .spType).
     * @param i_charsEn1 characteristics of entity type en1 (has to provide the member .spType) which will be updated via bit-wise | if an adjacent en0 entity has the respective sparse type.
     * @param i_enDaMe mapping of entity ids: data-to-mesh.
     * @param i_enMeDa mapping of element ids: mesh-to-data.
     *
     * @paramt TL_T_LID integer type of local ids.
     * @paramt TL_T_INT_SP integer type of the sparse type.
     * @paramt TL_T_EN0_CHARS struct of the en0 entities' characteristics. Offers a member .spType for comparison with i_spType.
     * @paramt TL_T_EN1_CHARS struct of the en1 entities' characteristics. Offers a member .spType for comparison with i_spType.
     **/
    template< typename TL_T_LID,
              typename TL_T_INT_SP,
              typename TL_T_EN0_CHARS,
              typename TL_T_EN1_CHARS >
    static void propAdj( TL_T_LID                   i_nEn0,
                         TL_T_LID   const * const * i_en0En1,
                         TL_T_INT_SP                    i_spTypeIn,
                         TL_T_INT_SP                    i_spTypeOut,
                         TL_T_EN0_CHARS         const * i_charsEn0,
                         TL_T_EN1_CHARS               * o_charsEn1,
                         TL_T_LID               const * i_enDaMe = nullptr,
                         TL_T_LID               const * i_enMeDa = nullptr ) {
      // maximum number of modified output entities
      TL_T_LID l_nEn1 = 0;

      // iterate over en0
      for( TL_T_LID l_en = 0; l_en < i_nEn0; l_en++ ) {
        // check for sparse type
        if( (i_charsEn0[l_en].spType & i_spTypeIn) == i_spTypeIn ) {
          // iterate over adjacent entities
          for( unsigned short l_ae = 0; l_ae < i_en0En1[l_en+1]-i_en0En1[l_en]; l_ae++ ) {
            TL_T_LID l_aeId = i_en0En1[l_en][l_ae];

            // propagate the info (if entity exists)
            if( l_aeId != std::numeric_limits< TL_T_LID >::max() ) {
              o_charsEn1[l_aeId].spType |= i_spTypeOut;
              l_nEn1 = std::max( l_nEn1, l_aeId+1 );
            }
          }
        }
      }

      // synchronize output characteristics
      if( i_enDaMe != nullptr ) {
        EDGE_CHECK_NE( i_enMeDa, nullptr );
        syncDupl( l_nEn1,
                  i_enDaMe,
                  i_enMeDa,
                  o_charsEn1 );
      }
    }

    /**
     * Links sparse entities to their dense counter-part.
     *
     * Example:
     *
     *   Input (implicit):    Output
     *
     *    SpId | DeId         SpId | DeId
     *         | 0               0 | 2
     *         | 1               1 | 4
     *       0 | 2               2 | 5
     *         | 3               3 | 7
     *       1 | 4
     *       2 | 5
     *         | 6
     *       3 | 7
     *         | 8
     *
     * @param i_nDe number of dense entities.
     * @param i_spType sparse type for which the link is determined.
     * @param i_chars entity characteristics with a sparse type.
     * @param o_spDe will be set to link between the sparse and dense elements.
     *
     * @paramt TL_T_INT_LID integer type of local ids.
     * @paramt TL_T_INT_SP integer type of the sparse type.
     * @paramt TL_T_CHARS struct of the dense entities' characteristics. Offers a member .spType for comparison with i_spType.
     **/
    template< typename TL_T_INT_LID,
              typename TL_T_INT_SP,
              typename TL_T_CHARS >
    static void linkSpDe( TL_T_INT_LID       i_nDe,
                          TL_T_INT_SP        i_spType,
                          TL_T_CHARS const * i_chars,
                          TL_T_INT_LID     * o_spDe ) {
      // init sparse id
      TL_T_INT_LID l_sp = 0;

      // iterate over dense entities
      for( TL_T_INT_LID l_de = 0; l_de < i_nDe; l_de++ ) {
        if( (i_chars[l_de].spType & i_spType) == i_spType ) {
          o_spDe[l_sp] = l_de;
          l_sp++;
        }
      }
    }

    /**
     * Links sparse entities to their sparse counter-part.
     *
     * Example:
     *
     *   Input (implicit):    Output
     *
     *    SpId0 | SpId1 | DeId         SpId0 | SpId1
     *          | 0     | 0                0 | 1
     *          |       | 1                1 | -
     *        0 | 1     | 2                2 | -
     *          |       | 3                3 | 2
     *        1 |       | 4
     *        2 |       | 5            undefined values (- in the table above)
     *          |       | 6            are set to std::numeric_limits< TL_T_LID >::max()
     *        3 | 2     | 7
     *          |       | 8
     *
     * @param i_nDe number dense entities.
     * @param i_spType0 first sparse type for which the link is determined.
     * @param i_spType1 second sparse type for which the link is determined.
     * @param i_chars entity characteristics with a sparse type.
     * @param o_sp0Sp1 will be set to link between the sparse elements.
     *
     * @paramt TL_T_LID integer type of local ids.
     * @paramt TL_T_SP integer type of the sparse type.
     * @paramt TL_T_CHARS struct of the dense entities' characteristics. Offers a member .spType for comparison with i_spType.
     **/
    template< typename TL_T_LID,
              typename TL_T_SP,
              typename TL_T_CHARS >
    static void linkSpSp( TL_T_LID           i_nDe,
                          TL_T_SP            i_spType0,
                          TL_T_SP            i_spType1,
                          TL_T_CHARS const * i_chars,
                          TL_T_LID         * o_sp0Sp1 ) {
      // init sparse ids
      TL_T_LID l_sp0 = 0;
      TL_T_LID l_sp1 = 0;

      // iterate over dense entities
      for( TL_T_LID l_de = 0; l_de < i_nDe; l_de++ ) {
        if( (i_chars[l_de].spType & i_spType0) == i_spType0 ) {
          if( (i_chars[l_de].spType & i_spType1) == i_spType1 )
            o_sp0Sp1[l_sp0] = l_sp1;
          else
            o_sp0Sp1[l_sp0] = std::numeric_limits< TL_T_LID >::max();

          l_sp0++;
        }

        if( (i_chars[l_de].spType & i_spType1) == i_spType1 )
          l_sp1++;
      }
    }

    /**
     * Links sparse entities based on adjacency information (single sparse type).
     *
     * Single vs. double sparse type:
     *   The case of single sparse type, only considers the to-entities having the sparse type.
     *   The case of double sparse type, assumes that both side have sparse types defined and uses
     *   the from-sparse ids for indexing.
     *
     * Example (single sparse type):
     *
     *   Dense adjacency      faces with        dense illustration  compressed sparse info
     *   information of       sparse ids        of derived          as given as output
     *   elements connected   (implicit info):  adjacency info:     of this function:
     *   to faces:
     *
     *   el | fa              fa | spId         el | spId           elSpId | faSpId
     *   0  | 1-4             0  | -            0  | x-x            0      | 2-1
     *   1  | 5-3             1  | -            1  | 2-1            1      | 3-x
     *   2  | 9-7             2  | 0            2  | 3-x            2      | x-1
     *   3  | 8-3             3  | 1            3  | x-1            3      | 0-x
     *   4  | 2-1             4  | -            4  | 0-x            4      | 1-x
     *   5  | 4-6             5  | 2            5  | x-x            5      | 3-0
     *   6  | 3-8             6  | -            6  | 1-x
     *   7  | 0-4             7  | -            7  | x-x
     *   8  | 9-2             8  | -            8  | 3-0
     *   9  | 4-1             9  | 3            9  | x-x
     *   10 | 4-1                               10 | x-x
     *
     * Remark: In the case only a subset is linked per entity, invalid ids
     *         std::numeric_limits<  TL_T_INT_LID >::max() are set.
     *         Example: Invalid are all faces with an x in the rightmost table.
     *
     * @param i_nEn number of dense entities which having adjacent entities.
     * @param i_nAdjPerEn number of adjacent entities (to) for each of the dense entities (from).
     * @param i_enEn adjacency information.nAdjPerEn. Adjacency information with std::numeric_limits< TL_T_INT_LID > is ignored (assuming boundary conditions).
     * @param i_spType sparse type used for the bit comparisons.
     * @param i_charsAdj characteristics of the adjecent entities (having a member .spType).
     * @param o_spLink will be set to the linked sparse ids of the adjacent entities. Assumed is a flat array, meaning: [de*i_nAdjPerEn] gives the de's dense entity and [de*i_nAdjPerEn + ae] the ae's adjacent entity.
     *
     * @paramt TL_T_INT_LID integer type of local ids.
     * @paramt TL_T_INT_SP integer type of the sparse type.
     * @paramt TL_T_ADJ_CHARS struct of the adjacent entities' characteristics. Offers a member .spType for comparison with i_spType.
     **/
    template< typename TL_T_INT_LID,
              typename TL_T_INT_SP,
              typename TL_T_ADJ_CHARS >
    static void linkSpAdjSst( TL_T_INT_LID           i_nEn,
                              unsigned short         i_nAdjPerEn,
                              TL_T_INT_LID   const * i_enEn,
                              TL_T_INT_SP            i_spType,
                              TL_T_ADJ_CHARS const * i_charsAdj,
                              TL_T_INT_LID         * o_spLink ) {
      // nothing to do if no entities are given
      if( i_nEn == 0 ) return;
      EDGE_CHECK( i_nAdjPerEn > 0 );

      // get number of adjacent entities to consider
      TL_T_INT_LID l_nAdjEn = 0;
      for( TL_T_INT_LID l_de = 0; l_de < i_nEn; l_de++ ) {
        for( unsigned short l_ae = 0; l_ae < i_nAdjPerEn; l_ae++ ) {
          TL_T_INT_LID l_aeId = i_enEn[ l_de * i_nAdjPerEn + l_ae ];

          // ignore undefined adjacencies
          if( l_aeId == std::numeric_limits< TL_T_INT_LID >::max() ) continue;
          l_nAdjEn = std::max( l_nAdjEn, l_aeId );
        }
      }

      // increase by one: size, not index
      l_nAdjEn++;

      // assemble lookup for the sparse id
      std::vector< TL_T_INT_LID > l_adjDeToSp;
      l_adjDeToSp.resize( l_nAdjEn );
      TL_T_INT_LID l_spId = 0;
      for( TL_T_INT_LID l_de = 0; l_de < l_nAdjEn; l_de++ ) {
        if( ( i_charsAdj[l_de].spType & i_spType ) == i_spType ) {
          l_adjDeToSp[l_de] = l_spId;
          l_spId++;
        }
        else l_adjDeToSp[l_de] = std::numeric_limits< TL_T_INT_LID >::max();
      }

      // set up compressed sparse info
      l_spId = 0;
      for( TL_T_INT_LID l_de = 0; l_de < i_nEn; l_de++ ) {
        // check if this an entity with adjacent entities having the desired sparse tag
        bool l_sp = false;
        for( unsigned short l_ae = 0; l_ae < i_nAdjPerEn; l_ae++ ) {
          TL_T_INT_LID l_aeId = i_enEn[ l_de * i_nAdjPerEn + l_ae ];

          // ignore undefined adjacencies
          if( l_aeId == std::numeric_limits< TL_T_INT_LID >::max() ) continue;

          if( ( i_charsAdj[l_aeId ].spType & i_spType ) == i_spType ) l_sp = true;
        }

        // link the sparse entities
        if( l_sp ) {
          for( unsigned short l_ae = 0; l_ae < i_nAdjPerEn; l_ae++ ) {
            TL_T_INT_LID l_aeId = i_enEn[ l_de * i_nAdjPerEn + l_ae ];

            // ignore undefined adjacencies
            if( l_aeId == std::numeric_limits< TL_T_INT_LID >::max() ) continue;

            if( ( i_charsAdj[l_aeId].spType & i_spType ) == i_spType ) o_spLink[ l_spId*i_nAdjPerEn + l_ae ] = l_adjDeToSp[l_aeId];
            else                                                       o_spLink[ l_spId*i_nAdjPerEn + l_ae ] = std::numeric_limits< TL_T_INT_LID >::max();
          }
          l_spId++;
        }
      }
    }

    /**
     * Links sparse entities based on adjacency information (double sparse type).
     *   This implementation assumes a flat array for the adjacency information (i_enEn) and
     *   a fixed number of to-entities for every from entity (i_nAdjPerEn).
     *
     * Single vs. double sparse type:
     *   The case of single sparse type, only considers the to-entities having the sparse type.
     *   The case of double sparse type, assumes that both side have sparse types defined and uses
     *   the from-sparse ids for indexing.
     *
     *   An example use case for double sparse types are internal boundaries (faces),
     *   which propagate their sparse info to adjacent elements, "internal boundary elements".
     *   Now, if a mapping faToEl from faces to elements is assumed as input, we are only interested
     *   in faces who were originally part of the internal boundary, not all faces adjacent to
     *   "internal boundary elements".
     *
     * Example (double sparse type, flat)
     *
     *   Dense adjacency          elements with     dense illustration  compressed sparse info
     *   information of           sparse ids        of derived          as given as output
     *   faces connected          (implicit info):  adjacency info:     of this function:
     *   to elements:
     *
     *   fa | el  | spId           el | spId         fa | spId           faSpId | elSpId
     *   0  | 1-4 |                0  | -            0  | x-x            0      | 2-1
     *   1  | 5-3 |  0             1  | -            1  | 2-1            1      | 3-x
     *   2  | 9-7 |  1             2  | 0            2  | 3-x            2      | x-1
     *   3  | 8-3 |  2             3  | 1            3  | x-1            3      | 0-x
     *   4  | 2-1 |  3             4  | -            4  | 0-x            4      | 1-x
     *   5  | 2-6 | ***            5  | 2            5  | x-x            5      | 3-0
     *   6  | 3-8 |  4             6  | -            6  | 1-x
     *   7  | 0-4 |                7  | -            7  | x-x
     *   8  | 9-2 |  5             8  | -            8  | 3-0
     *   9  | 4-1 |                9  | 3            9  | x-x
     *   10 | 4-9 | ***                              10 | x-x
     *
     *   Faces marked with spId *** are adjacent to "internal boundary elements", but not
     *   returned in the compressed sparse info, since the don't have the sparse flag set.
     *
     *
     * Remark: In the case only a subset is linked per entity, invalid ids
     *         std::numeric_limits<  TL_T_INT_LID >::max() are set.
     *         Example: Invalid are all faces with an x in the rightmost table.
     *
     * @param i_nEn number of dense entities which having adjacent entities.
     * @param i_nAdjPerEn number of adjacent entities (to) for each of the dense entities (from).
     * @param i_enEn adjacency information. Adjacency information with std::numeric_limits< TL_T_INT_LID >::max being ignored (assuming boundary conditions).
     * @param i_spTypeFrom sparse type used for the bit comparisons of the from-entities.
     * @param i_spTypeTo sparse type used for the bit comparisons of the to-entities.
     * @param i_charsFrom characteristics of the from entities (having a member .spType).
     * @param i_charsTo characteristics of the adjecent entities (having a member .spType).
     * @param o_spLink will be set to the linked sparse ids of the adjacent entities. Assumed is a flat array, meaning: [de*i_nAdjPerEn] gives the de's dense entity and [de*i_nAdjPerEn + ae] the ae's adjacent entity.
     *
     * @paramt TL_T_INT_LID integer type of local ids.
     * @paramt TL_T_INT_SP integer type of the sparse type.
     * @paramt TL_T_FROM_CHARS struct of the from entities' characteristics. Offers a member .spType for comparison with i_spType.
     * @paramt TL_T_TO_CHARS struct of the adjacent entities' characteristics. Offers a member .spType for comparison with i_spType.
     **/
    template< typename TL_T_INT_LID,
              typename TL_T_INT_SP,
              typename TL_T_CHARS_FROM,
              typename TL_T_CHARS_TO >
    static void linkSpAdjDst( TL_T_INT_LID            i_nEn,
                              unsigned short          i_nAdjPerEn,
                              TL_T_INT_LID    const * i_enEn,
                              TL_T_INT_SP             i_spTypeFrom,
                              TL_T_INT_SP             i_spTypeTo,
                              TL_T_CHARS_FROM const * i_charsFrom,
                              TL_T_CHARS_TO   const * i_charsTo,
                              TL_T_INT_LID          * o_spLink ) {
      // nothing to do if no entities are given
      if( i_nEn == 0 ) return;
      EDGE_CHECK( i_nAdjPerEn > 0 );

      // get number of adjacent entities to consider
      TL_T_INT_LID l_nAdjEn = 0;
      for( TL_T_INT_LID l_de = 0; l_de < i_nEn; l_de++ ) {
        for( unsigned short l_ae = 0; l_ae < i_nAdjPerEn; l_ae++ ) {
          TL_T_INT_LID l_aeId = i_enEn[ l_de * i_nAdjPerEn + l_ae ];

          // ignore undefined adjacencies
          if( l_aeId == std::numeric_limits< TL_T_INT_LID >::max() ) continue;

          l_nAdjEn = std::max( l_nAdjEn, l_aeId );
        }
      }

      // increase by one: size, not index
      l_nAdjEn++;

      // assemble lookup for the sparse id
      std::vector< TL_T_INT_LID > l_adjDeToSp;
      l_adjDeToSp.resize( l_nAdjEn );
      TL_T_INT_LID l_spId = 0;
      for( TL_T_INT_LID l_de = 0; l_de < l_nAdjEn; l_de++ ) {
        if( ( i_charsTo[l_de].spType & i_spTypeTo ) == i_spTypeTo ) {
          l_adjDeToSp[l_de] = l_spId;
          l_spId++;
        }
        else l_adjDeToSp[l_de] = std::numeric_limits< TL_T_INT_LID >::max();
      }

      // set up compressed sparse info
      l_spId = 0;
      for( TL_T_INT_LID l_de = 0; l_de < i_nEn; l_de++ ) {
        // link the sparse entities
        if( (i_charsFrom[l_de].spType & i_spTypeFrom) == i_spTypeFrom ) {
          for( unsigned short l_ae = 0; l_ae < i_nAdjPerEn; l_ae++ ) {
            TL_T_INT_LID l_aeId = i_enEn[ l_de * i_nAdjPerEn + l_ae ];

            // ignore undefined adjacencies
            if( l_aeId == std::numeric_limits< TL_T_INT_LID >::max() ) {
              o_spLink[ l_spId*i_nAdjPerEn + l_ae ] = std::numeric_limits< TL_T_INT_LID >::max();
              continue;
            }

            if( ( i_charsTo[l_aeId].spType & i_spTypeTo ) == i_spTypeTo ) o_spLink[ l_spId*i_nAdjPerEn + l_ae ] = l_adjDeToSp[l_aeId];
            else                                                          o_spLink[ l_spId*i_nAdjPerEn + l_ae ] = std::numeric_limits< TL_T_INT_LID >::max();
          }

          l_spId++;
        }
      }
    }

    /**
     * Determines the sizes of the link for sparse entities based on adjacency information (double sparse type).
     *   This implementation assumes dynamic storage for the adjacency information (i_enEn).
     *   E.g., for vertices as bridge between elements, a dynamic number of adjacent to-elements per from-entity is allowed.
     *
     * Single vs. double sparse type:
     *   The case of single sparse type, only considers the to-entities having the sparse type.
     *   The case of double sparse type, assumes that both side have sparse types defined and uses
     *   the from-sparse ids for indexing.
     *
     * Example (double sparse type, dynamic)
     *
     *   Dense adjacency                elements with      dense illustration  compressed sparse info.
     *   information of                 sparse ids         of derived          the number of entries of both columns
     *   vertices connected             (implicit info):   adjacency info:     as output of this function:
     *   to elements:
     *
     *   ve | el        | spId           el | spId         ve | spId           veSpId | elSpId
     *   0  | 1-4       |  -             0  | -            0  | x-x            0      | 2-1
     *   1  | 8-5-0-3   |  0             1  | -            1  | x-2-x-1        1      | -
     *   2  | 8-7-4     |  1             2  | 0            2  | x-x-x          2      | 1
     *   3  | 8-3-6-7-8 |  2             3  | 1            3  | x-1-x-x-x      3      | 0
     *   4  | 2-1       |  3             4  | -            4  | 0-x            4      | 1
     *   5  | 2-6-5-3   |  -             5  | 2            5  | x-x-x-x        5      | 3-0-2-1
     *   6  | 3-8       |  4             6  | -            6  | 1-x
     *   7  | 0-4       |  -             7  | -            7  | x-x
     *   8  | 9-2-5-3-7 |  5             8  | -            8  | 3-0-2-1-x
     *   9  | 4-1-9-8   |  -             9  | 3            9  | x-x-x-x
     *   10 | 4-9       |  -                               10 | x-x
     *
     * @param i_nEn number of dense entities having adjacent entities.
     * @param i_enEn adjacency information. A ghost cell at position i_nEn is assumed, giving the number of adjacent to-entities for the last from-entity (i_nEn-1).
     * @param i_spTypeFrom sparse type used for the bit comparisons of the from-entities.
     * @param i_spTypeTo sparse type used for the bit comparisons of the to-entities.
     * @param i_charsFrom characteristics of the from-entities (having a member .spType).
     * @param i_charsTo characteristics of the adjecent to-entities (having a member .spType).
     * @param o_nSpLinkRaw number of link items in terms of adjacenct to-entities.
     * @param o_nSpLinkPtr number of link items in terms of the number of pointers, excludes ghost entry.
     *
     * @paramt TL_T_INT_LID integer type of local ids.
     * @paramt TL_T_INT_SP integer type of the sparse type.
     * @paramt TL_T_CHARS_FROM struct of the from entities' characteristics. Offers a member .spType for comparison with i_spType.
     * @paramt TL_T_CHARS_TO struct of the adjacent entities' characteristics. Offers a member .spType for comparison with i_spType.
     **/
    template< typename TL_T_INT_LID,
              typename TL_T_INT_SP,
              typename TL_T_CHARS_FROM,
              typename TL_T_CHARS_TO >
    static void nLinkSpAdjDst( TL_T_INT_LID                    i_nEn,
                               TL_T_INT_LID    const * const * i_enEn,
                               TL_T_INT_SP                     i_spTypeFrom,
                               TL_T_INT_SP                     i_spTypeTo,
                               TL_T_CHARS_FROM         const * i_charsFrom,
                               TL_T_CHARS_TO           const * i_charsTo,
                               TL_T_INT_LID                  & o_nSpLinkRaw,
                               TL_T_INT_LID                  & o_nSpLinkPtr ) {
      // init output
      o_nSpLinkRaw = o_nSpLinkPtr = 0;

      // get number of adjacent entities to consider
      for( TL_T_INT_LID l_de = 0; l_de < i_nEn; l_de++ ) {
        // ignore from-entities that do not match the from sparse type
        if( (i_charsFrom[l_de].spType & i_spTypeFrom) != i_spTypeFrom ) continue;

        // increase counter for pointers
        o_nSpLinkPtr++;

        // iterate over adjacent entities
        for( unsigned short l_ae = 0; l_ae < i_enEn[l_de+1]-i_enEn[l_de]; l_ae++ ) {
          TL_T_INT_LID l_aeId = i_enEn[l_de][l_ae];

          // require valid ids
          EDGE_CHECK( l_aeId != std::numeric_limits< TL_T_INT_LID >::max() );

          // increase counter for raw data
          if( (i_charsTo[l_aeId].spType & i_spTypeTo) == i_spTypeTo ) o_nSpLinkRaw++;
        }
      }
    }

    /**
     * Determines the link for sparse entities based on adjacency information (double sparse type).
     *   This implementation assumes dynamic storage for the adjacency information (i_enEn).
     *   E.g., for vertices as bridge between elements, a dynamic number of adjacent to-elements per from-entity is allowed.
     *
     * See nLinkSpAdjDst for an example. In contrast this function outputs the actual link and not the number of entries.
     *
     * @param i_nEn number of dense entities having adjacent entities.
     * @param i_enEn adjacency information. A ghost cell at position i_nEn is assumed, giving the number of adjacent to-entities for the last from-entity (i_nEn-1).
     * @param i_spTypeFrom sparse type used for the bit comparisons of the from-entities.
     * @param i_spTypeTo sparse type used for the bit comparisons of the to-entities.
     * @param i_charsFrom characteristics of the from-entities (having a member .spType).
     * @param i_charsTo characteristics of the adjecent to-entities (having a member .spType).
     * @param o_spLinkRaw raw entries of the sparse link.
     * @param o_spLinkPtr pointers to entries for every sparse from-entity of the link.
     *
     * @paramt TL_T_INT_LID integer type of local ids.
     * @paramt TL_T_INT_SP integer type of the sparse type.
     * @paramt TL_T_CHARS_FROM struct of the from entities' characteristics. Offers a member .spType for comparison with i_spType.
     * @paramt TL_T_CHARS_TO struct of the adjacent entities' characteristics. Offers a member .spType for comparison with i_spType.
     **/
    template< typename TL_T_INT_LID,
              typename TL_T_INT_SP,
              typename TL_T_CHARS_FROM,
              typename TL_T_CHARS_TO >
    static void linkSpAdjDst( TL_T_INT_LID                    i_nEn,
                              TL_T_INT_LID    const * const * i_enEn,
                              TL_T_INT_SP                     i_spTypeFrom,
                              TL_T_INT_SP                     i_spTypeTo,
                              TL_T_CHARS_FROM         const * i_charsFrom,
                              TL_T_CHARS_TO           const * i_charsTo,
                              TL_T_INT_LID                  * o_spLinkRaw,
                              TL_T_INT_LID                 ** o_spLinkPtr ) {
      // init output
      o_spLinkPtr[0] = o_spLinkRaw;

      // get number of adjacent entities to consider
      TL_T_INT_LID l_nAdjEn = 0;
      for( TL_T_INT_LID l_de = 0; l_de < i_nEn; l_de++ ) {
        for( unsigned short l_ae = 0; l_ae < i_enEn[l_de+1]-i_enEn[l_de]; l_ae++ ) {
          TL_T_INT_LID l_aeId = i_enEn[l_de][l_ae];
          l_nAdjEn = std::max( l_nAdjEn, l_aeId );
        }
      }

      // increase by one: size, not index
      l_nAdjEn++;

      // assemble lookup for the sparse to-id
      std::vector< TL_T_INT_LID > l_adjDeToSp;
      l_adjDeToSp.resize( l_nAdjEn );

      TL_T_INT_LID l_spId = 0;
      for( TL_T_INT_LID l_de = 0; l_de < l_nAdjEn; l_de++ ) {
        if( (i_charsTo[l_de].spType & i_spTypeTo) == i_spTypeTo ) {
          l_adjDeToSp[l_de] = l_spId;
          l_spId++;
        }
        else l_adjDeToSp[l_de] = std::numeric_limits< TL_T_INT_LID >::max();
      }


      // set up compressed sparse info
      l_spId = 0;
      TL_T_INT_LID * l_raw = o_spLinkRaw;

      for( TL_T_INT_LID l_de = 0; l_de < i_nEn; l_de++ ) {
        // link the sparse entities
        if( (i_charsFrom[l_de].spType & i_spTypeFrom) == i_spTypeFrom ) {
          // set raw pointer
          o_spLinkPtr[l_spId] = l_raw;

          for( unsigned short l_ae = 0; l_ae < i_enEn[l_de+1]-i_enEn[l_de]; l_ae++ ) {
            TL_T_INT_LID l_aeId = i_enEn[l_de][l_ae];

            if( (i_charsTo[l_aeId].spType & i_spTypeTo) == i_spTypeTo ) {
              *l_raw = l_adjDeToSp[l_aeId];
              l_raw++;
            }
          }

          l_spId++;
        }
      }

      // set ghost entry
      o_spLinkPtr[l_spId] = l_raw;
    }
};
#endif
