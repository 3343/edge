/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2017-2019, Regents of the University of California
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
 * Output of receivers at sub-faces.
 **/

#ifndef EDGE_IO_RECEIVERS_SF_HPP
#define EDGE_IO_RECEIVERS_SF_HPP

#include "io/logging.h"
#include "Receivers.h"
#include <string>
#include "data/EntityLayout.type"
#include "linalg/Mappings.hpp"
#include "linalg/Geom.hpp"
#include "FileSystem.hpp"
#include <limits>

namespace edge {
  namespace io {
    template< typename       TL_T_REAL,
              t_entityType   TL_T_EL,
              unsigned short TL_O_SP,
              unsigned short TL_N_CRS >
    class ReceiversSf;
  }
}

/**
 * Sub-face version of the receivers.
 *
 * @paramt TL_T_REAL real type used for arithmetic operations.
 * @paramt TL_T_EL element type.
 * @paramt TL_T_O_SP order of the DG-scheme in space.
 * @paramt TL_N_CRS number of fused simulations.
 **/
template< typename       TL_T_REAL,
          t_entityType   TL_T_EL,
          unsigned short TL_O_SP,
          unsigned short TL_N_CRS >
class edge::io::ReceiversSf: public Receivers {
  private:
    //! dimension of the element
    static unsigned short const TL_N_DIS = C_ENT[TL_T_EL].N_DIM;

    //! number of sub-vertices
    static unsigned short const TL_N_SVS = CE_N_SUB_VERTICES( TL_T_EL, TL_O_SP );

    //! number of sub-cells per element
    static unsigned short const TL_N_SCS  = CE_N_SUB_CELLS( TL_T_EL, TL_O_SP );

    //! number of sub-faces per element face
    static unsigned short const TL_N_SFS = CE_N_SUB_FACES( TL_T_EL, TL_O_SP );

    //! number of vertices per faces
    static unsigned short const TL_N_VES_FA = C_ENT[TL_T_EL].N_FACE_VERTICES;

    //! number of vertices per element
    static unsigned short const TL_N_VES_EL = C_ENT[TL_T_EL].N_VERTICES;

    //! number of faces
    static unsigned short const TL_N_FAS = C_ENT[TL_T_EL].N_FACES;

    struct RecvSf {
      //! sub-face of the receiver
      unsigned short sf;
      //! coordinates of the receiver
      TL_T_REAL crds[TL_N_DIS];
    };

    //! properties specific to receivers at sub-faces
    std::vector< RecvSf > m_recvsSf;

  public:
    /**
     * Intitializes receivers at sub-faces.
     *
     * The receivers are based on the nearest found sub-face.
     * with the given sparse type.
     *
     * To avoid MPI-related issues,
     *   1) only owned elements' faces are considered.
     *   2) this face's left element has to match the considered owned element.
     *
     * @param i_nRecvs number of receivers.
     * @param i_spType sparse type of the face's considered for receiver output.
     * @param i_nQts number of quantities per receiver (type: TL_T_REAL).
     * @param i_outDir output directory for the csv-files.
     * @param i_recvNames names of the receivers.
     * @param i_recvCrds coordinates of the receivers.
     * @param i_freq frequency of the receiver output.
     * @param i_faLayout entity layout of the faces.
     * @param i_elLayout entity layout of the elements.
     * @param i_scSv sub-vertices, adjacent to sub-cells.
     * @param i_faEl elements adjacent to the faces.
     * @param i_elVe vertices adjacent to the elements.
     * @param i_elFa faces adjacent to the elements.
     * @param i_svChars sub-vertex characteristics.
     * @param i_veChars vertex characteristics.
     * @param io_faChars, face characteristics, member .spType will be updated with the bitmask RECEIVER, if a receiver is located at this face.
     * @param i_bufferSize size of the receiver-buffer.
     * @param i_time initial time of the receivers.
     *
     * @paramt TL_T_INT_SP integer type of the sparse type.
     * @paramt TL_T_LAYOUT structure of the entity layout.
     * @paramt TL_T_INT_LID integer type of local entity ids.
     * @paramt TL_T_CHARS_SV sub-vertex type, offering member .coords.
     * @paramt TL_T_CHARS_VE structure of the vertex characteristics. provides .coords member for the vertices' coordinates.
     * @paramt TL_T_CHARS_FA structure of the face characteristics. provides .spType member for access to the faces' sparse types.

     **/
    template< typename TL_T_INT_SP,
              typename TL_T_LAYOUT,
              typename TL_T_INT_LID,
              typename TL_T_CHARS_VE,
              typename TL_T_CHARS_FA,
              typename TL_T_CHARS_SC >
    void init( unsigned int           i_nRecvs,
               TL_T_INT_SP            i_spType,
               unsigned short         i_nQts,
               std::string    const  &i_outDir,
               std::string    const (*i_recvNames),
               double         const (*i_recvCrds)[TL_N_DIS],
               double                 i_freq,
               TL_T_LAYOUT    const  &i_faLayout,
               TL_T_LAYOUT    const  &i_elLayout,
               unsigned short const   i_scSv[ TL_N_SCS + TL_N_FAS * TL_N_SFS ][ TL_N_VES_EL ],
               TL_T_INT_LID   const (*i_faEl)[2],
               TL_T_INT_LID   const (*i_elVe)[TL_N_VES_EL],
               TL_T_INT_LID   const (*i_elFa)[TL_N_FAS],
               TL_T_CHARS_SC          i_svChars[TL_N_SVS],
               TL_T_CHARS_VE  const  *i_veChars,
               TL_T_CHARS_FA         *io_faChars,
               unsigned int           i_bufferSize=100,
               double                 i_time=0  ) {
      // minimal distance to all receivers
      std::vector< double > l_minDist(i_nRecvs);
      // id of face having the minimum distance to the receivers
      std::vector< TL_T_INT_LID > l_minFa(i_nRecvs);
      // id of sub-face point having the minimum distance to the receivers
      std::vector< unsigned short > l_minSf(i_nRecvs);
      // sub-face coordinates
      std::vector< double > l_sfCrds[TL_N_DIS];
      for( unsigned short l_di = 0; l_di < TL_N_DIS; l_di++ ) l_sfCrds[l_di].resize( i_nRecvs );

      // init minimum data structures
      for( unsigned int l_re = 0; l_re < i_nRecvs; l_re++ ) {
        l_minDist[l_re] = std::numeric_limits< double >::max();
        l_minFa[l_re]   = std::numeric_limits< TL_T_INT_LID   >::max();
        l_minSf[l_re]   = std::numeric_limits< unsigned short >::max();

        for( unsigned short l_di = 0; l_di < TL_N_DIS; l_di++ ) {
          l_sfCrds[l_di][l_re] = std::numeric_limits< double >::max();
        }
      }


      TL_T_INT_LID l_first = 0;

      // iterate over all time groups
      for( std::size_t l_tg = 0; l_tg < i_elLayout.timeGroups.size(); l_tg++ ) {
        // iterate over all owned entities
        TL_T_INT_LID l_size = i_elLayout.timeGroups[l_tg].nEntsOwn;

        for( TL_T_INT_LID l_el = l_first; l_el < l_first+l_size; l_el++ ) {
          // iterate over the element's faces
          for( unsigned short l_fa = 0; l_fa < TL_N_FAS; l_fa++ ) {
            // determine the dense id of the face
            TL_T_INT_LID l_faId = i_elFa[l_el][l_fa];

            if( l_faId < i_faLayout.nEnts                          && // valid face
                (io_faChars[l_faId].spType & i_spType) == i_spType && // continue if sparse type matches
                 i_faEl[l_faId][0] == l_el ) {                        // continue if "left" == trivial mapping

              // iterate over sub-faces
              for( unsigned short l_sf = 0; l_sf < TL_N_SFS; l_sf++ ) {
                // assemble vertex coordinates
                double l_veCrds[TL_N_DIS][TL_N_VES_EL];
                for( unsigned short l_ve = 0; l_ve < TL_N_VES_EL; l_ve++ ) {
                  TL_T_INT_LID l_veId = i_elVe[l_el][l_ve];

                  for( unsigned short l_di = 0; l_di < TL_N_DIS; l_di++ ) {
                    l_veCrds[l_di][l_ve] = i_veChars[l_veId].coords[l_di];
                  }
                }

                // derive mesh coordinates of sub-face by averaging
                double l_meshCrds[TL_N_DIS];

                // TODO: generalize
                double l_refCrds[3] = {0, 0, 0}; // TODO: work around for non-template
                for( unsigned short l_sv = 0; l_sv < TL_N_VES_FA; l_sv++ ) {
                  // get id of the sub-vertex
                  unsigned short l_svId = i_scSv[TL_N_SCS + l_fa * TL_N_SFS + l_sf][l_sv];

                  for( unsigned short l_di = 0; l_di < TL_N_DIS; l_di++ ) {
                    // add contribution of the sub-face
                    l_refCrds[l_di] += i_svChars[l_svId].coords[l_di] * ( double(1) / TL_N_VES_FA);
                  }
                }


                linalg::Mappings::refToPhy( TL_T_EL,
                                            l_veCrds[0],
                                            l_refCrds,
                                            l_meshCrds );

                // iterate over receivers and determine the distance
                for( unsigned int l_re = 0; l_re < i_nRecvs; l_re++ ) {
                  double l_dist = linalg::GeomT<TL_N_DIS>::norm( l_meshCrds,
                                                                 i_recvCrds[l_re] );

                  // store the info if this a new minimum
                  if( l_minDist[l_re] > l_dist ) {
                    l_minDist[l_re] = l_dist;
                    l_minFa[l_re]   = l_faId;
                    l_minSf[l_re]   = l_sf;

                    for( unsigned short l_di = 0; l_di < TL_N_DIS; l_di++ ) {
                      l_sfCrds[l_di][l_re] = l_meshCrds[l_di];
                    }
                  }
                }
              }
            }
          }
        }

        // update first element
        l_first += l_size;
        l_first += i_elLayout.timeGroups[l_tg].nEntsNotOwn;
    }

    // determine if we hold the sub-face with the minimum distance to the receiver
    std::vector< unsigned short > l_recvOwn( i_nRecvs );
    parallel::Mpi::min( i_nRecvs,
                        l_minDist.data(),
                        l_recvOwn.data() );

    /*
     * set up the receivers data structures
     */
    // buffer size
    m_buffSize = i_bufferSize;

    // #vars
    m_nQts = i_nQts;
    // frequency
    m_freq = i_freq;

    TL_T_INT_LID l_spId = 0;

    // number of sparse receiver entities
    TL_T_INT_LID l_spRe = 0;

    l_first = 0;

    for( std::size_t l_tg = 0; l_tg < i_faLayout.timeGroups.size(); l_tg++ ) {
      TL_T_INT_LID l_size  = i_faLayout.timeGroups[l_tg].nEntsOwn;
                   l_size += i_faLayout.timeGroups[l_tg].nEntsNotOwn;

      TL_T_INT_LID l_spIdTg = 0;

      for( TL_T_INT_LID l_fa = l_first; l_fa < l_first+l_size; l_fa++ ) {
        if( (io_faChars[l_fa].spType & i_spType) == i_spType ) {
          bool l_reFa = false;

          // iterate over receivers
          for( unsigned int l_re = 0; l_re < i_nRecvs; l_re++ ) {
            if(    l_recvOwn[l_re] == 1
                && l_minFa[l_re] == l_fa ) {
              // add this receiver
              m_recvs.resize( m_recvs.size() + 1 );
              m_recvsSf.resize( m_recvsSf.size() + 1 );

              // update face type
              io_faChars[l_fa].spType |= RECEIVER;

              // init receiver data
              m_recvs.back().nBuff = 0;
              m_recvs.back().buffer.resize( i_bufferSize*i_nQts*TL_N_CRS );
              m_recvs.back().buffTime.resize( m_buffSize );
              m_recvs.back().time  = i_time;
              m_recvs.back().tg    = l_tg;
              m_recvs.back().en    = l_spId;
              m_recvs.back().enTg  = l_spIdTg;
              std::string l_dir    = i_outDir + "/" + std::to_string(parallel::g_rank);
              m_recvs.back().path  = l_dir + "/" + i_recvNames[l_re]+".csv";

              m_recvsSf.back().sf = l_minSf[l_re];
              for( unsigned short l_di = 0; l_di < TL_N_DIS; l_di++ )
                m_recvsSf.back().crds[l_di] = l_sfCrds[l_di][l_re];

              // set sparse entitity ids of the receivers
              if( l_reFa == false ) {
                m_spEnToRecv.resize( l_spRe+1 );
                m_spEnToRecv[l_spRe] = m_recvs.size()-1;
                l_reFa = true;
              }
            }
          }
          if( l_reFa == true ) l_spRe++;

          l_spId++;
          l_spIdTg++;
        }
      }
      l_first += l_size;
    }

    // create directories and touch output
    if( i_nRecvs > 0 ) {
      std::string l_dirCreate = i_outDir + "/" + std::to_string(parallel::g_nRanks);
      FileSystem::createDir( l_dirCreate );

      touchOutput();
    }
  }

  /**
   * Writes the output for all receivers associated with the sparse receiver entity.
   * "All" receivers: A single sparse receiver entity might contain more than one
   *                  sub-face and thus might have multiple receivers per entity.
   *
   * @param i_time time of the receiver output.
   * @param i_dt time step, the maximum of the receiver's frequency and the time step will be used to determine the next output point.
   * @param i_reSp sparse entity.
   * @param i_data data for the receiver output. [*][][]: quantities, [][*][]: sub-faces, all of them per entity, [][][*]: fused forward simulations.
   *
   * @paramt TL_T_INT_SP integer type of the sparse entities.
   **/
  template< typename TL_T_INT_SP >
  void writeRecvAll( double             i_time,
                     double             i_dt,
                     TL_T_INT_SP        i_reSp,
                     TL_T_REAL   const (*i_data)[TL_N_SFS][TL_N_CRS] ) {
    // get first matching receiver
    std::size_t l_first = m_spEnToRecv[i_reSp];
    EDGE_CHECK_LT( l_first, m_recvs.size() );
    int_el l_firstEn = m_recvs[l_first].en;

    // iterate over all receivers corresponding to this receiver entity
    for( std::size_t l_re = l_first; l_re < m_recvs.size(); l_re++ ) {
      // receiver matches the sparse entity
      if( m_recvs[l_re].en == l_firstEn ) {
        // check our bufffer isn't overflowing
        EDGE_CHECK_LT( m_recvs[l_re].nBuff, m_buffSize );

        // set quantities in buffer
        for( unsigned short l_qt = 0; l_qt < m_nQts; l_qt++ ) {
          for( unsigned short l_ru = 0; l_ru < TL_N_CRS; l_ru++ ) {
            unsigned int l_pos  = m_recvs[l_re].nBuff*m_nQts*TL_N_CRS;
                         l_pos += l_qt*TL_N_CRS;
                         l_pos += l_ru;
            m_recvs[l_re].buffer[l_pos] = i_data[l_qt][m_recvsSf[l_re].sf][l_ru];
          }
        }

        // set time in buffer
        m_recvs[l_re].buffTime[ m_recvs[l_re].nBuff ] = i_time;

        // update receiver stats
        m_recvs[l_re].time = i_time + std::max( m_freq, i_dt );

        m_recvs[l_re].nBuff++;
      }
      // receiver doesn't match, abort (ascending order)
      else break;
    }
  }
};
#endif
