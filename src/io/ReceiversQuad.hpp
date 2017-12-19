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
 * Output of receivers at quadrature points.
 **/

#ifndef RECEIVERS_QUAD_HPP
#define RECEIVERS_QUAD_HPP

#include "io/logging.h"
#include "Receivers.h"
#include <string>
#include "data/layout.hpp"
#include "linalg/Mappings.hpp"
#include "linalg/Geom.hpp"
#include <limits>

namespace edge {
  namespace io {
    template< typename TL_T_REAL, t_entityType TL_T_EL, unsigned short TL_O_QUAD, unsigned short TL_N_CRUNS >
    class ReceiversQuad;
  }
}

/**
 * Quadrature version of the receivers.
 *
 * @paramt TL_T_REAL real type used for arithmetic operations.
 * @paramt TL_T_EL element type.
 * @paramt TL_T_O_QUAD order of the quadrature.
 * @paramt TL_N_CRUNS number of fused simulations.
 **/
template< typename TL_T_REAL, t_entityType TL_T_EL, unsigned short TL_O_QUAD, unsigned short TL_N_CRUNS >
class edge::io::ReceiversQuad: public Receivers {
  private:
    //! dimension of the element
    static unsigned short const TL_N_DIM = C_ENT[TL_T_EL].N_DIM;

    //! number of vertices per element
    static unsigned short const TL_N_EL_VE = C_ENT[TL_T_EL].N_VERTICES;

    //! number of vertex options
    static unsigned short const TL_N_FA_VEOPS = CE_N_FACE_VERTEX_OPTS(TL_T_EL);

    //! number of quad points per face
    static unsigned short const TL_N_FA_QPTS = CE_N_FACE_QUAD_POINTS( TL_T_EL, TL_O_QUAD );

    //! number of faces
    static unsigned short const TL_N_FA = C_ENT[TL_T_EL].N_FACES;

    struct RecvQuad {
      //! quadrature point of the receiver
      unsigned short qp;
      //! coordinates of the receiver
      TL_T_REAL crds[TL_N_DIM];
    };

    //! properties specific to receivers at quadrature points
    std::vector< RecvQuad > m_recvsQuad;

  public:
    /**
     * Prints statistics of the receivers.
     **/
    void print() {
      for( std::size_t l_re = 0; l_re < m_recvsQuad.size(); l_re++ ) {
        // get location
        std::string l_loc = "";
        for( unsigned short l_di = 0; l_di < TL_N_DIM; l_di++ ) {
          l_loc += " ";
          l_loc += std::to_string( m_recvsQuad[l_re].crds[l_di] );
        }

        EDGE_LOG_INFO_ALL << "    receiver " << m_recvs[l_re].path
                          << " is located at" <<  l_loc;
      }
    }

    /**
     * Intitializes receivers at quadrature points of faces.
     *
     * The receivers are based on the nearest found quadrature point at faces
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
     * @param i_qPts location of the the element face's quadrature points. Since only left elements are considered, no vertex options have to be considered. [*][][]: faces in the reference element. [][*][]: quadrature points per face, [][][*]: reference dimensions.
     * @param i_faLayout entity layout of the faces.
     * @param i_elLayout entity layout of the elements.
     * @param i_faEl elements adjacent to the faces.
     * @param i_elVe vertices adjacent to the elements.
     * @param i_elFa faces adjacent to the elements.
     * @param i_veChars vertex characteristics.
     * @param io_faChars, face characteristics, member .spType will be updated with the bitmask RECEIVER, if a receiver is located at this face.
     * @param i_bufferSize size of the receiver-buffer.
     * @param i_time initial time of the receivers.
     *
     * @paramt TL_T_INT_SP integer type of the sparse type.
     * @paramt TL_T_LAYOUT structure of the entity layout.
     * @paramt TL_T_INT_LID integer type of local entity ids.
     * @paramt TL_T_CHARS_VE structure of the vertex characteristics. provides .coords member for the vertices' coordinates.
     * @paramt TL_T_CHARS_FA structure of the face characteristics. provides .spType member for access to the faces' sparse types.
     * @paramt TL_T_REAL_MESH precision of mesh-releated data
     **/
    template< typename TL_T_INT_SP,
              typename TL_T_LAYOUT,
              typename TL_T_INT_LID,
              typename TL_T_CHARS_VE,
              typename TL_T_CHARS_FA,
              typename TL_T_REAL_MESH >
    void init( unsigned int           i_nRecvs,
               TL_T_INT_SP            i_spType,
               unsigned short         i_nQts,
               std::string    const  &i_outDir,
               std::string    const (*i_recvNames),
               TL_T_REAL_MESH const (*i_recvCrds)[TL_N_DIM],
               double                i_freq,
               TL_T_REAL_MESH const   i_qPts[TL_N_FA][TL_N_FA_QPTS][TL_N_DIM],
               TL_T_LAYOUT    const  &i_faLayout,
               TL_T_LAYOUT    const  &i_elLayout,
               TL_T_INT_LID   const (*i_faEl)[2],
               TL_T_INT_LID   const (*i_elVe)[TL_N_EL_VE],
               TL_T_INT_LID   const (*i_elFa)[TL_N_FA],
               TL_T_CHARS_VE  const  *i_veChars,
               TL_T_CHARS_FA         *io_faChars,
               unsigned int           i_bufferSize=100,
               double                 i_time=0  ) {
      // minimal distance to all receivers
      std::vector< TL_T_REAL_MESH > l_minDist(i_nRecvs);
      // id of face having the minimum distance to the receivers
      std::vector< TL_T_INT_LID   > l_minFa(i_nRecvs);
      // id of quad point having the minimum distance to the receivers
      std::vector< unsigned short > l_minQp(i_nRecvs);
      // quad points coordinates
      std::vector< TL_T_REAL_MESH > l_qpCrds[TL_N_DIM];
      for( unsigned short l_di = 0; l_di < TL_N_DIM; l_di++ ) l_qpCrds[l_di].resize( i_nRecvs );

      // init minimum data structures
      for( unsigned int l_re = 0; l_re < i_nRecvs; l_re++ ) {
        l_minDist[l_re] = std::numeric_limits< TL_T_REAL_MESH >::max();
        l_minFa[l_re]   = std::numeric_limits< TL_T_INT_LID   >::max();
        l_minQp[l_re]   = std::numeric_limits< unsigned short >::max();

        for( unsigned short l_di = 0; l_di < TL_N_DIM; l_di++ ) {
          l_qpCrds[l_di][l_re] = std::numeric_limits< TL_T_REAL_MESH>::max();
        }
      }


      TL_T_INT_LID l_first = 0;

      // iterate over all time groups
      for( std::size_t l_tg = 0; l_tg < i_elLayout.timeGroups.size(); l_tg++ ) {
        // iterate over all owned entities
        TL_T_INT_LID l_size = i_elLayout.timeGroups[l_tg].nEntsOwn;

        for( TL_T_INT_LID l_el = l_first; l_el < l_first+l_size; l_el++ ) {
          // iterate over the element's faces
          for( unsigned short l_fa = 0; l_fa < TL_N_FA; l_fa++ ) {
            // determine the dense id of the face
            TL_T_INT_LID l_faId = i_elFa[l_el][l_fa];

            if( l_faId < i_faLayout.nEnts                          && // valid face
                (io_faChars[l_faId].spType & i_spType) == i_spType && // continue if sparse type matches
                 i_faEl[l_faId][0] == l_el ) {                        // continue if "left" == trivial mapping

              // iterate over the quad points
              for( unsigned short l_qp = 0; l_qp < TL_N_FA_QPTS; l_qp++ ) {
                // assemble vertex coordinates
                TL_T_REAL_MESH l_veCrds[TL_N_DIM][TL_N_EL_VE];
                for( unsigned short l_ve = 0; l_ve < TL_N_EL_VE; l_ve++ ) {
                  TL_T_INT_LID l_veId = i_elVe[l_el][l_ve];

                  for( unsigned short l_di = 0; l_di < TL_N_DIM; l_di++ ) {
                    l_veCrds[l_di][l_ve] = i_veChars[l_veId].coords[l_di];
                  }
                }

                // derive mesh coorinates of the quad point
                TL_T_REAL_MESH l_meshCrds[TL_N_DIM];

                // TODO: generalize
                TL_T_REAL_MESH l_refCrds[3] = {0, 0, 0}; // TODO: work around for non-template
                for( unsigned short l_di = 0; l_di < TL_N_DIM; l_di++ ) {
                  l_refCrds[l_di] = i_qPts[l_fa][l_qp][l_di];
                }


                linalg::Mappings::refToPhy( TL_T_EL,
                                            l_veCrds[0],
                                            l_refCrds,
                                            l_meshCrds );

                // iterate over receivers and determine the distance
                for( unsigned int l_re = 0; l_re < i_nRecvs; l_re++ ) {
                  TL_T_REAL_MESH l_dist = linalg::GeomT<TL_N_DIM>::norm( l_meshCrds,
                                                                         i_recvCrds[l_re] );

                  // store the info if this a new minimum
                  if( l_minDist[l_re] > l_dist ) {
                    l_minDist[l_re] = l_dist;
                    l_minFa[l_re]   = l_faId;
                    l_minQp[l_re]   = l_qp;

                    for( unsigned short l_di = 0; l_di < TL_N_DIM; l_di++ ) {
                      l_qpCrds[l_di][l_re] = l_meshCrds[l_di];
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

    // TODO: Eliminate MPI-duplicates here
    EDGE_CHECK( i_elLayout.timeGroups[0].neRanks.size() == 0 );

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
            if( l_minFa[l_re] == l_fa ) {
              // add this receiver
              m_recvs.resize( m_recvs.size() + 1 );
              m_recvsQuad.resize( m_recvsQuad.size() + 1 );

              // update face type
              io_faChars[l_fa].spType |= i_spType;

              // init receiver data
              m_recvs.back().nBuff = 0;
              m_recvs.back().buffer.resize( i_bufferSize*i_nQts*TL_N_CRUNS );
              m_recvs.back().buffTime.resize( m_buffSize );
              m_recvs.back().time  = i_time;
              m_recvs.back().tg    = l_tg;
              m_recvs.back().en    = l_spId;
              m_recvs.back().enTg  = l_spIdTg;
              m_recvs.back().path  = i_outDir+"/"+i_recvNames[l_re]+".csv";

              m_recvsQuad.back().qp = l_minQp[l_re];
              for( unsigned short l_di = 0; l_di < TL_N_DIM; l_di++ )
                m_recvsQuad.back().crds[l_di] = l_qpCrds[l_di][l_re];

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

    // touch output
    if( i_nRecvs > 0 ) touchOutput( i_outDir );
  }

  /**
   * Writes the output for all receivers associated with the sparse receiver entity.
   * "All" receivers: A single sparse receiver entity might contain more than one
   *                  quadrature point and thus might have multiple receivers per entity.
   *
   * @param i_time time of the receiver output.
   * @param i_dt time step, the maximum of the receiver's frequency and the time step will be used to determine the next output point.
   * @param i_reSp sparse entity.
   * @param i_data data for the receiver output. [*][][]: quantities, [][*][]: quadrature points, all of them per entity, [][][*]: fused forward simulations.
   *
   * @paramt TL_T_INT_SP integer type of the sparse entities.
   **/
  template< typename TL_T_INT_SP >
  void writeRecvAll( double             i_time,
                     double             i_dt,
                     TL_T_INT_SP        i_reSp,
                     TL_T_REAL   const (*i_data)[TL_N_FA_QPTS][TL_N_CRUNS] ) {
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
          for( unsigned short l_ru = 0; l_ru < TL_N_CRUNS; l_ru++ ) {
            unsigned int l_pos  = m_recvs[l_re].nBuff*m_nQts*TL_N_CRUNS;
                         l_pos += l_qt*TL_N_CRUNS;
                         l_pos += l_ru;
            m_recvs[l_re].buffer[l_pos] = i_data[l_qt][m_recvsQuad[l_re].qp][l_ru];
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
