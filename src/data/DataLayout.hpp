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
 * Special purpose data layouts.
 **/
#ifndef EDGE_DATA_DATA_LAYOUT_HPP
#define EDGE_DATA_DATA_LAYOUT_HPP

#include "io/logging.h"
#include "EntityLayout.type"

namespace edge {
  namespace data {
    class DataLayout;
  }
}

/**
 * Special purpose data layouts.
 **/
class edge::data::DataLayout {
  public:
    /**
     * Assembles a pointered data structure, adjacency-dependent data.
     * For example, if an sparse element provides data for face-adjacent sparse elements,
     * the corresponding face-pointers would only be set, if the adjacent element has the corresponding sparse type.
     * 
     * Additionally, the assembly respects communication-related store, not requiring artificial MPI-buffers.
     * For example, if we communicate data of faces, all of these entries for a communication region would
     * be stored linearly in memory.
     *
     * Here, entity duplication is assumed.
     * Write something about redundancy for send -> inner.
     * but not for send -> ghost.
     *
     * @param i_nAdPerEn number of adjacent entities for each entity.
     * @param i_enAdEn adjacent entities.
     * @param i_lay entity layout.
     * @param i_raw location of the raw data.
     * @param o_ptrs will be set pointers into the raw data.
     * @param o_send will be set to pointers into the raw data, for use within MPI-send. Size is the number of MPI regions+1 for the derivation of messages sizes.
     * @param o_recv will be set to pointers into the raw data, for use within MPI-recv. SIze is the number of MPI regions+1 for message size derivation.
     * @param i_enSpEn additional sparse "adjacency" (no bridge) for determining if the adjacent entity should be considered.
     *
     * @paramt TL_T_INT_LID integral type for the local ids.
     * @paramt TL_T_RAW type of the raw data.
     * @paramt TL_T_PTR type of the pointered data, given as array of [i_enAdEn] pointers to TL_T_RAW.
     **/
    template< typename TL_T_INT_LID,
              typename TL_T_RAW,
              typename TL_T_PTR >
    static void adj( unsigned short          i_nAdPerEn,
                     TL_T_INT_LID    const * i_enAdEn,
                     t_enLayout      const & i_lay,
                     TL_T_RAW              * i_raw,
                     TL_T_PTR              * o_ptrs,
                     std::vector<
                       std::vector<
                         unsigned char *
                       >
                      >                    & o_send,
                     std::vector<
                       std::vector<
                         unsigned char *
                       >
                      >                    & o_recv,
                      TL_T_INT_LID const   * i_enSpEn = nullptr ) {
    // ptr to raw data
    TL_T_RAW * l_raw = i_raw;

    // init all pointers invalid
    for( TL_T_INT_LID l_en = 0; l_en < i_lay.nEnts; l_en++ )
      for( TL_T_INT_LID l_ad = 0; l_ad < i_nAdPerEn; l_ad++ )
        o_ptrs[l_en][l_ad] = nullptr;

    // iterate over time groups
    for( std::size_t l_tg = 0; l_tg < i_lay.timeGroups.size(); l_tg++ ) {
      // add empty communication data
      o_send.resize( o_send.size()+1 );
      o_recv.resize( o_recv.size()+1 );

      // iterate over inner entities
      TL_T_INT_LID l_inFirst = i_lay.timeGroups[l_tg].inner.first;
      TL_T_INT_LID l_inSize  = i_lay.timeGroups[l_tg].inner.size;
      for( TL_T_INT_LID l_in = l_inFirst; l_in < l_inFirst+l_inSize; l_in++ ) {
        // iterate over adjacent entities
        for( unsigned short l_ad = 0; l_ad < i_nAdPerEn; l_ad++ ) {
          TL_T_INT_LID l_enAd = i_enAdEn[ l_in * i_nAdPerEn + l_ad ];

          // ignore if that doesn't match our sparse type
          if( l_enAd == std::numeric_limits< TL_T_INT_LID >::max() ) continue;
          if(    i_enSpEn != nullptr
              && i_enSpEn[l_in] == std::numeric_limits< TL_T_INT_LID >::max()
              && i_enSpEn[l_enAd] == std::numeric_limits< TL_T_INT_LID >::max() ) continue;

          // check if the adjacency is defined
          EDGE_CHECK_EQ( o_ptrs[l_in][l_ad], nullptr );
          // assign data
          o_ptrs[l_in][l_ad] = l_raw;
          // move on
          l_raw++;
        }
      }

      // iterate over send regions and set everything which goes to owned entities
      for( std::size_t l_sr = 0; l_sr < i_lay.timeGroups[l_tg].send.size(); l_sr++ ) {
        // iterate over send entities
        TL_T_INT_LID l_seFirst = i_lay.timeGroups[l_tg].send[l_sr].first;
        TL_T_INT_LID l_seSize = i_lay.timeGroups[l_tg].send[l_sr].size;
        for( TL_T_INT_LID l_se = l_seFirst; l_se < l_seFirst+l_seSize; l_se++ ) {
          // iterate over adjacent entities
          for( unsigned short l_ad = 0; l_ad < i_nAdPerEn; l_ad++ ) {
            // determine adjacent entity
            TL_T_INT_LID l_enAd = i_enAdEn[ l_se * i_nAdPerEn + l_ad ];

            // ignore if that doesn't match our sparse type
            if( l_enAd == std::numeric_limits< TL_T_INT_LID >::max() ) continue;
            if(    i_enSpEn != nullptr
                && i_enSpEn[l_se] == std::numeric_limits< TL_T_INT_LID >::max()
                && i_enSpEn[l_enAd] == std::numeric_limits< TL_T_INT_LID >::max() ) continue;

            // check if data for this adjacency should be added
            bool l_add = false;

            // iterate over time groups
            for( std::size_t l_at = 0; l_at < i_lay.timeGroups.size(); l_at++ ) {
              // get first and size of owned entities
              TL_T_INT_LID l_atFirst = i_lay.timeGroups[l_at].inner.first;
              TL_T_INT_LID l_atSize  = i_lay.timeGroups[l_at].nEntsOwn;

              // check if the adjacent entity is part of this
              if(    l_enAd >= l_atFirst
                  && l_enAd  < l_atFirst+l_atSize ) l_add = true;
              break;
            }

            if( l_add ) {
              EDGE_CHECK_EQ( o_ptrs[l_se][l_ad], nullptr );
              // assign data
              o_ptrs[l_se][l_ad] = l_raw;
              // move on
              l_raw++;
            }
          }
        }
      }

      // iterate over send regions and set everything, which goes to corresponding ghost-entities
      for( std::size_t l_sr = 0; l_sr < i_lay.timeGroups[l_tg].send.size(); l_sr++ ) {
        // add communication data for the send-region
        o_send.back().push_back( (unsigned char *) l_raw );

        TL_T_INT_LID l_seFirst = i_lay.timeGroups[l_tg].send[l_sr].first;
        TL_T_INT_LID l_seSize  = i_lay.timeGroups[l_tg].send[l_sr].size;

        // corresponding receive region
        TL_T_INT_LID l_reFirst = i_lay.timeGroups[l_tg].receive[l_sr].first;
        TL_T_INT_LID l_reSize  = i_lay.timeGroups[l_tg].receive[l_sr].size;

        // iterate over send entities
        for( TL_T_INT_LID l_se = l_seFirst; l_se < l_seFirst+l_seSize; l_se++ ) {
          // iterate over adjacent entities
          for( unsigned short l_ad = 0; l_ad < i_nAdPerEn; l_ad++ ) {
            // determine adjacent entity
            TL_T_INT_LID l_enAd = i_enAdEn[ l_se * i_nAdPerEn + l_ad ];

            // ignore if that doesn't match our sparse type
            if( l_enAd == std::numeric_limits< TL_T_INT_LID >::max() ) continue;
            if(    i_enSpEn != nullptr
                && i_enSpEn[l_se] == std::numeric_limits< TL_T_INT_LID >::max()
                && i_enSpEn[l_enAd] == std::numeric_limits< TL_T_INT_LID >::max() ) continue;

            // add, if and only if this connects us to the respective receive region
            if(    l_enAd >= l_reFirst
                && l_enAd  < l_reFirst+l_reSize ) {
              EDGE_CHECK_EQ( o_ptrs[l_se][l_ad], nullptr );
              // assign data
              o_ptrs[l_se][l_ad] = l_raw;
              // move on
              l_raw++;
            }
          }
        }
      }

      // add final pointer for size computations of messages
      o_send.back().push_back( (unsigned char *) l_raw );

      // iterate over recv regions and set everything, which goes to corresponding send-entities
      for( std::size_t l_rr = 0; l_rr < i_lay.timeGroups[l_tg].receive.size(); l_rr++ ) {
        // add communication data for the recv-region
        o_recv.back().push_back( (unsigned char *) l_raw );

        TL_T_INT_LID l_reFirst = i_lay.timeGroups[l_tg].receive[l_rr].first;
        TL_T_INT_LID l_reSize  = i_lay.timeGroups[l_tg].receive[l_rr].size;

        // corresponding send region
        TL_T_INT_LID l_seFirst = i_lay.timeGroups[l_tg].send[l_rr].first;
        TL_T_INT_LID l_seSize  = i_lay.timeGroups[l_tg].send[l_rr].size;

        // iterate over receive entities
        for( TL_T_INT_LID l_re = l_reFirst; l_re < l_reFirst+l_reSize; l_re++ ) {
          // iterate over adjacent entities
          for( unsigned short l_ad = 0; l_ad < i_nAdPerEn; l_ad++ ) {
            // determine adjacent entity
            TL_T_INT_LID l_enAd = i_enAdEn[ l_re * i_nAdPerEn + l_ad ];

            // ignore if that doesn't match our sparse type
            if( l_enAd == std::numeric_limits< TL_T_INT_LID >::max() ) continue;
            if(    i_enSpEn != nullptr
                && i_enSpEn[l_re] == std::numeric_limits< TL_T_INT_LID >::max()
                && i_enSpEn[l_enAd] == std::numeric_limits< TL_T_INT_LID >::max() ) continue;

            // add, if and only if this connects us to the respective send region
            if(    l_enAd >= l_seFirst
                && l_enAd  < l_seFirst+l_seSize ) {
              EDGE_CHECK_EQ( o_ptrs[l_re][l_ad], nullptr );

              // assign data
              o_ptrs[l_re][l_ad] = l_raw;
              // move on
              l_raw++;
            }
          }
        }
      }

      // add final pointer for size computations of messages
      o_recv.back().push_back( (unsigned char *) l_raw );

    }
  }
};

#endif
