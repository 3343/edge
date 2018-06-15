/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2016, Regents of the University of California
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
#include "parallel/Mpi.h"
#include "Receivers.h"
#include "FileSystem.hpp"
#include "linalg/Geom.hpp"
#include "dg/Basis.h"
#include "io/logging.h"
#include <set>
#include <fstream>
#include <sstream>

void edge::io::Receivers::print() {
  // rank's local receivers
  unsigned int l_nRecvsL = m_recvs.size();
  // global receivers
  unsigned int l_nRecvsG;
#ifdef PP_USE_MPI
  int l_err = MPI_Allreduce( &l_nRecvsL, &l_nRecvsG, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD );
  EDGE_CHECK_EQ( l_err, MPI_SUCCESS );
#else
  l_nRecvsG = l_nRecvsL;
#endif
  EDGE_LOG_INFO << "  we have " << l_nRecvsG << " active receivers in the domain" << ", buffer size: " << m_buffSize;
}

void edge::io::Receivers::init(       t_entityType    i_enType,
                                      unsigned int    i_nRecvs,
                                const std::string    &i_outDir,
                                const std::string   (*i_recvNames),
                                const real_mesh     (*i_recvCrds)[3],
                                      double          i_freq,
                                const t_enLayout     &i_enLayout,
                                const int_el        (*i_enVe),
                                const t_vertexChars  *i_veChars,
                                      unsigned int    i_bufferSize,
                                      double          i_time ) {
  // init class-wide vars
  m_freq = i_freq;
  m_buffSize = i_bufferSize;
  m_nQts = N_QUANTITIES;

  unsigned short l_nVe = C_ENT[i_enType].N_VERTICES;

  // set of receivers to find
  std::set< unsigned int > l_rSet;
  for( unsigned int l_rc = 0; l_rc < i_nRecvs; l_rc++ ) l_rSet.insert( l_rc );

  // iterate over the time groups
  int_el l_first = 0;

  for( int_tg l_tg = 0; l_tg < i_enLayout.timeGroups.size(); l_tg++ ) {
    int_el l_size  = i_enLayout.timeGroups[l_tg].nEntsOwn;

    // iterate over the owned entities
    for( int_el l_en = l_first; l_en < l_first+l_size; l_en++ ) {
      // buffer ves
      EDGE_CHECK_LE( l_nVe, 8 );
      real_mesh l_tmpVe[ 3*8 ];

      // get the vertices
      for( unsigned short l_ve = 0; l_ve < l_nVe; l_ve++ ) {
        int_el l_veId = i_enVe[l_en*l_nVe+l_ve];

        for( unsigned short l_di = 0; l_di < 3; l_di++ ) {
          l_tmpVe[l_di*l_nVe + l_ve] = i_veChars[l_veId].coords[l_di];
        }
      }

      // iterate over the non-found receivers
      for( std::set< unsigned int >::iterator l_rc = l_rSet.begin(); l_rc != l_rSet.end(); ) {
        // check if the receiver is inside
        if( linalg::Geom::inside( i_enType, l_tmpVe, i_recvCrds[*l_rc] ) != 0 ) {
          // set info
          if( m_recvs.size() > 0 && m_recvs.back().en < l_en ) {
            m_spEnToRecv.push_back( m_recvs.size() );
          }
          else if( m_recvs.size() == 0 ) m_spEnToRecv.push_back( 0 );

          // add a new receiver
          m_recvs.resize( m_recvs.size()+1 );

          m_recvs.back().buffer.resize( N_QUANTITIES*N_CRUNS*m_buffSize );
          m_recvs.back().buffTime.resize( m_buffSize );
          m_recvs.back().nBuff = 0;
          m_recvs.back().time  = i_time;
          m_recvs.back().tg    = l_tg;
          m_recvs.back().en    = l_en;
          m_recvs.back().enTg  = l_en-l_first;
          m_recvs.back().path  = i_outDir+"/"+i_recvNames[*l_rc]+".csv";

          // determine the location in reference coordinates
          real_mesh l_ref[3] = {0,0,0};
          linalg::Mappings::phyToRef( i_enType, l_tmpVe, i_recvCrds[*l_rc], l_ref );

          // check for reasonable coords
          for( unsigned short l_di = 0; l_di < C_ENT[i_enType].N_DIM; l_di++ ) {
            EDGE_CHECK_GT( l_ref[l_di], -TOL.MESH );
            EDGE_CHECK_LT( l_ref[l_di], 1+TOL.MESH );
          }

          // evaluate the basis at the given locations
          for( int_md l_md = 0; l_md < N_ELEMENT_MODES; l_md++ ) {
            dg::Basis::evalBasis( l_md, T_SDISC.ELEMENT, m_recvs.back().evaBasis[l_md], l_ref[0], l_ref[1], l_ref[2] );
          }

          // delete this receiver from the list of non-found receivers
          l_rc = l_rSet.erase(l_rc);
        }
        else l_rc++;
      }
    }
    l_first += i_enLayout.timeGroups[l_tg].nEntsOwn +
               i_enLayout.timeGroups[l_tg].nEntsNotOwn;
  }

  // touch output
  if( i_nRecvs > 0 ) touchOutput( i_outDir );
}

void edge::io::Receivers::touchOutput( const std::string &i_outDir ) {
  // create ouput-directories if not present
  if( parallel::g_rank == 0 )
    FileSystem::createDir( i_outDir );
#ifdef PP_USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  // define a header
  std::string l_header = "time";
  for( int_qt l_qt = 0; l_qt < m_nQts; l_qt++ ) {
    for( int_cfr l_cr = 0; l_cr < N_CRUNS; l_cr++ ) {
      l_header += ",Q" + std::to_string(l_qt) + "_C" + std::to_string(l_cr);
    }
  }
  l_header += '\n';

  // write header
  for( std::size_t l_re = 0; l_re < m_recvs.size(); l_re++ ) {
    std::ofstream l_file;
    l_file.open( m_recvs[l_re].path );

    if( l_file.is_open() ) {
      l_file << l_header;
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
