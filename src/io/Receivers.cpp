/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2016-2019, Regents of the University of California
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
 * Output of receivers.
 **/
#include "data/SparseEntities.hpp"
#include "Receivers.h"
#include "FileSystem.hpp"
#include "linalg/Geom.hpp"
#include "dg/Basis.h"
#include "io/logging.h"
#include <set>
#include <fstream>
#include <sstream>

void edge::io::Receivers::init( t_entityType             i_enType,
                                unsigned short           i_nTgs,
                                std::size_t    const   * i_nTgEnsIn,
                                std::size_t    const   * i_nTgEnsSe,
                                unsigned int             i_nRecvs,
                                std::string    const   & i_outDir,
                                std::string    const ( * i_recvNames),
                                double         const ( * i_recvCrds)[3],
                                double                   i_freq,
                                std::size_t    const ( * i_enVe),
                                t_vertexChars  const   * i_veChars,
                                unsigned int             i_bufferSize,
                                double                   i_time ) {
  // init class-wide vars
  m_freq = i_freq;
  m_buffSize = i_bufferSize;
  m_nQts = N_QUANTITIES;

  std::size_t l_nEns = 0;
  for( unsigned short l_tg = 0; l_tg < i_nTgs; l_tg++ )
    l_nEns += i_nTgEnsIn[l_tg] + i_nTgEnsSe[l_tg];

  // number of vertices
  unsigned short l_nVe = C_ENT[i_enType].N_VERTICES;

  // derive the dense ids of the receivers
  int_el *l_deIds = new int_el[i_nRecvs];
  edge::data::SparseEntities::ptToEn( i_enType,
                                      int_el(i_nRecvs),
                                      i_recvCrds,
                                      l_nEns,
                                      i_enVe,
                                      i_veChars,
                                      l_deIds );

  // lambda which inits the receivers for an entity if any
  auto l_addEn = [ i_enVe,
                   i_enType,
                   i_nRecvs,
                   i_veChars,
                   i_recvCrds,
                   i_time,
                   i_outDir,
                   i_recvNames,
                   l_deIds,
                   l_nVe,
                   this ]( unsigned short i_tg,
                           std::size_t    i_en,
                           std::size_t    i_enTg ) {
    // iterate over the receivers
    for( unsigned int l_re = 0; l_re < i_nRecvs; l_re++ ) {
      // only continue if the receiver is owned and matches the element
      if( i_en == l_deIds[l_re] ) {
        // get the vertices of the entity
        real_mesh l_tmpVe[ 3 * 8];
        for( unsigned short l_ve = 0; l_ve < l_nVe; l_ve++ ) {
          int_el l_veId = i_enVe[i_en*l_nVe+l_ve];

          for( unsigned short l_di = 0; l_di < 3; l_di++ ) {
            l_tmpVe[l_di*l_nVe + l_ve] = i_veChars[l_veId].coords[l_di];
          }
        }

        // init the receiver coordinates
        real_mesh l_recvCrds[3];
        for( unsigned short l_di = 0; l_di < 3; l_di++ )
          l_recvCrds[l_di] = i_recvCrds[l_re][l_di];

        // project receiver to the element's surface if required
        edge::linalg::Geom::closestPoint( i_enType,
                                          l_tmpVe,
                                          l_recvCrds );

        // set info
        if( m_recvs.size() > 0 && m_recvs.back().en < i_en ) {
          m_spEnToRecv.push_back( m_recvs.size() );
        }
        else if( m_recvs.size() == 0 ) m_spEnToRecv.push_back( 0 );

        // add a new receiver
        m_recvs.resize( m_recvs.size()+1 );
        m_recvs.back().id = l_re;
        for( unsigned short l_di = 0; l_di < 3; l_di++ )
          m_recvs.back().coords[l_di] = l_recvCrds[l_di];
        m_recvs.back().buffer.resize( N_QUANTITIES*N_CRUNS*m_buffSize );
        m_recvs.back().buffTime.resize( m_buffSize );
        m_recvs.back().nBuff = 0;
        m_recvs.back().time  = i_time;
        m_recvs.back().tg    = i_tg;
        m_recvs.back().en    = i_en;
        m_recvs.back().enTg  = i_enTg;
        std::string l_dir = i_outDir + "/" + std::to_string(parallel::g_rank);
        m_recvs.back().path  = l_dir + "/" + i_recvNames[l_re]+".csv";

        // determine the location in reference coordinates
        real_mesh l_ref[3] = {0,0,0};
        linalg::Mappings::phyToRef( i_enType, l_tmpVe, l_recvCrds, l_ref );

        // check for reasonable coords
        for( unsigned short l_di = 0; l_di < C_ENT[i_enType].N_DIM; l_di++ ) {
          EDGE_CHECK_GT( l_ref[l_di], -TOL.MESH );
          EDGE_CHECK_LT( l_ref[l_di], 1+TOL.MESH );
        }

        // evaluate the basis at the given locations
        for( int_md l_md = 0; l_md < N_ELEMENT_MODES; l_md++ ) {
          m_recvs.back().evaBasis[l_md] = dg::Basis::evalBasis( l_md,
                                                                T_SDISC.ELEMENT,
                                                                l_ref[0],
                                                                l_ref[1],
                                                                l_ref[2] );
        }
      }
    }
  };

  std::size_t l_off = 0;
  // inner
  for( unsigned short l_tg = 0; l_tg < i_nTgs; l_tg++ ) {
    for( std::size_t l_en = 0; l_en < i_nTgEnsIn[l_tg]; l_en++ ) {
      l_addEn( l_tg,
               l_en + l_off,
               l_en );
    }
    l_off += i_nTgEnsIn[l_tg];
  }
  // send
  for( unsigned short l_tg = 0; l_tg < i_nTgs; l_tg++ ) {
    for( std::size_t l_en = 0; l_en < i_nTgEnsSe[l_tg]; l_en++ ) {
      l_addEn( l_tg,
               l_en + l_off,
               l_en );
    }
    l_off += i_nTgEnsSe[l_tg];
  }

  // create directories and touch output
  if( m_recvs.size() > 0 ) {
    std::string l_dirCreate = i_outDir + "/" + std::to_string(parallel::g_rank);
    FileSystem::createDir( l_dirCreate );

    touchOutput( i_recvNames,
                 i_recvCrds );
  }

  // free memory
  delete[] l_deIds;
}

void edge::io::Receivers::touchOutput( std::string const (*i_recvNames),
                                       real_mesh   const (*i_recvCrds)[3] ) {
  // shared header
  std::string l_headerSh  = "# EDGE\n";
              l_headerSh += "# code version: " + std::string(PP_EDGE_VERSION) + "\n";
              l_headerSh += "# build date / time: " + std::string(__DATE__) + " / " + std::string(__TIME__) + "\n";

  // define column names
  std::string l_colNames = "time";
  for( int_qt l_qt = 0; l_qt < m_nQts; l_qt++ ) {
    for( int_cfr l_cr = 0; l_cr < N_CRUNS; l_cr++ ) {
      l_colNames += ",Q" + std::to_string(l_qt) + "_C" + std::to_string(l_cr);
    }
  }
  l_colNames += '\n';

  // write header
  for( std::size_t l_re = 0; l_re < m_recvs.size(); l_re++ ) {
    std::ofstream l_file;
    l_file.open( m_recvs[l_re].path );

    if( l_file.is_open() ) {
      unsigned int l_id = m_recvs[l_re].id;

      l_file << l_headerSh;

      if( i_recvNames != nullptr && i_recvCrds != nullptr ) {
        l_file << "# receiver name: " << i_recvNames[l_id] << "\n";
        l_file << "# specified coordinates: " << i_recvCrds[l_id][0] << " "
                                              << i_recvCrds[l_id][1] << " "
                                              << i_recvCrds[l_id][2] << "\n";
        l_file << "# projected coordinates: " << m_recvs[l_re].coords[0] << " "
                                              << m_recvs[l_re].coords[1] << " "
                                              << m_recvs[l_re].coords[2] << "\n";
      }

      l_file << l_colNames;
    }
    else EDGE_LOG_FATAL;
  }
}

void edge::io::Receivers::getEnRecv( std::vector< int_el > & o_en ) {
  o_en.resize( m_spEnToRecv.size() );
  for( std::size_t l_en = 0; l_en < m_spEnToRecv.size(); l_en++ ) {
    std::size_t l_recvId = m_spEnToRecv[l_en];

    EDGE_CHECK_LT( l_recvId, m_recvs.size() );

    o_en[l_en] = m_recvs[l_recvId].en;
  }
}

double edge::io::Receivers::getRecvTimeRel( int_el i_spEn,
                                            double i_time,
                                            double i_dt ) {
  unsigned int l_re = m_spEnToRecv[i_spEn];
  if( m_recvs[l_re].time > i_time-TOL.TIME && m_recvs[l_re].time < i_time+i_dt-TOL.TIME ) {
    // determine relative dt
    double l_dt = m_recvs[l_re].time - i_time;
    EDGE_CHECK_GT( l_dt, -TOL.TIME );
    return std::abs( l_dt );
  }

  // return an invalid value if no output is expected
  return -std::numeric_limits< double >::max();
}

void edge::io::Receivers::writeRecvAll(       int_el    i_spEn,
                                        const real_base i_dofs[N_QUANTITIES][N_ELEMENT_MODES][N_CRUNS] ) {
  EDGE_CHECK_LT( (std::size_t) i_spEn, m_spEnToRecv.size() );

  // get id and entity of first receiver
  std::size_t l_first = m_spEnToRecv[i_spEn];
  EDGE_CHECK_LT( l_first, m_recvs.size() );
  int_el l_firstEn = m_recvs[l_first].en;

  for( std::size_t l_re = l_first; l_re < m_recvs.size(); l_re++ ) {
    if( m_recvs[l_re].en == l_firstEn ) {
      // check our bufffer isn't overflowing
      EDGE_CHECK_LT( m_recvs[l_re].nBuff, m_buffSize );

      // iterate over quantities
      for( int_qt l_qt = 0; l_qt < N_QUANTITIES; l_qt++ ) {
        // reset buffer
        for( int_cfr l_cr = 0; l_cr < N_CRUNS; l_cr++ ) {
          m_recvs[l_re].buffer[ m_recvs[l_re].nBuff*N_QUANTITIES*N_CRUNS + l_qt*N_CRUNS + l_cr ] = 0;
        }

        // iterate over modes
        for( int_md l_md = 0; l_md < N_ELEMENT_MODES; l_md++ ) {
          // eval the DOFS and store values
          for( int_cfr l_cr = 0; l_cr < N_CRUNS; l_cr++ ) {
            m_recvs[l_re].buffer[ m_recvs[l_re].nBuff*N_QUANTITIES*N_CRUNS + l_qt*N_CRUNS + l_cr ] +=
              m_recvs[l_re].evaBasis[l_md] * i_dofs[l_qt][l_md][l_cr];
          }
        }
      }
      // set time
      m_recvs[l_re].buffTime[ m_recvs[l_re].nBuff ] = m_recvs[l_re].time;

      // update receiver stats
      m_recvs[l_re].time += m_freq;
      m_recvs[l_re].nBuff++;
    }
    // abort if all of sparse entities' associated receivers have been written
    else break;
  }
}

void edge::io::Receivers::flush( unsigned int i_re ) {
  std::ofstream l_file;
  if( m_recvs[i_re].nBuff > 0 ) {
    l_file.open( m_recvs[i_re].path, std::ios_base::app );

    if( l_file.is_open() ) {
      // stream buffer
      std::ostringstream l_stream;

      // assemble output stream
      for( unsigned int l_bu = 0; l_bu < m_recvs[i_re].nBuff; l_bu++ ) {
        // write time info
        l_stream << std::to_string( m_recvs[i_re].buffTime[l_bu] );
        // write recv values
        for( unsigned int l_va = 0; l_va < m_nQts*N_CRUNS; l_va++ ) {
          l_stream << "," << std::scientific << m_recvs[i_re].buffer[ l_bu*m_nQts*N_CRUNS+l_va ];
        }
        l_stream << "\n";
      }

      // write stream to file
      l_file << l_stream.str();
    }
    else EDGE_LOG_FATAL << "could not open the recv-file: " << m_recvs[i_re].path;
  }
  m_recvs[i_re].nBuff = 0;
}

void edge::io::Receivers::flushAll() {
  // iterate over all receivers
  for( std::size_t l_re = 0; l_re < m_recvs.size(); l_re++ ) flush( l_re );
}

void edge::io::Receivers::flushIf( unsigned int i_tresh ) {
  for( std::size_t l_re = 0; l_re < m_recvs.size(); l_re++ ) {
    if( m_buffSize - m_recvs[l_re].nBuff < i_tresh ) flush( l_re );
  }
}
