/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2016-2018, Regents of the University of California
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
 * Output for error norms.
 **/

#include "ErrorNorms.h"
#include <io/logging.h>
#include "submodules/pugixml/src/pugixml.hpp"

void edge::io::ErrorNorms::print( const double i_errorNorms[3][N_QUANTITIES][N_CRUNS] ) {
  EDGE_LOG_INFO << "fasten your seat belts, error norms cming next.";
  for( unsigned int l_norm = 0; l_norm < 3; l_norm++ ) {
    for( unsigned int l_q = 0; l_q < N_QUANTITIES; l_q++ ) {
      for( int_cfr l_cfr = 0; l_cfr < N_CRUNS; l_cfr++ ) {
        if( l_norm == 0 ) {
          EDGE_LOG_INFO << "  L1   error, quantity #"
                        << l_q << ", cfr #" << l_cfr << ": " << i_errorNorms[0][l_q][l_cfr];
        }
        else if( l_norm == 1 ) {
          EDGE_LOG_INFO << "  L2   error, quantity #"
                        << l_q << ", cfr #" << l_cfr << ": " << i_errorNorms[1][l_q][l_cfr];
        }
        else if( l_norm == 2 ) {
          EDGE_LOG_INFO << "  Linf error, quantity #"
                        << l_q << ", cfr #" << l_cfr << ": " << i_errorNorms[2][l_q][l_cfr];
        }
      }
    }
  }
}

void edge::io::ErrorNorms::writeXml( const double i_errorNorms[3][N_QUANTITIES][N_CRUNS] ) {
  if( parallel::g_rank == 0 ) {
    pugi::xml_document l_doc;

    pugi::xml_node l_normsNd = l_doc.append_child( "error_norms" );

    // iterate over norms
    for( unsigned int l_norm = 0; l_norm < 3; l_norm++ ) {
      pugi::xml_node l_normNd;

      // add norm nodes to xml
      if( l_norm == 0 )
        l_normNd = l_normsNd.append_child( "l1" );
      else if( l_norm == 1 )
        l_normNd = l_normsNd.append_child( "l2" );
      else if( l_norm == 2 )
        l_normNd = l_normsNd.append_child( "linf" );

      for( unsigned int l_q = 0; l_q < N_QUANTITIES; l_q++ ) {
        pugi::xml_node l_qNd;

        // add quantity to xml
        l_qNd = l_normNd.append_child( "q" );

        // iterate over runs
        for( int_cfr l_cfr = 0; l_cfr < N_CRUNS; l_cfr++ ) {
          pugi::xml_node l_cfrNd;

          l_cfrNd = l_qNd.append_child( "cfr" );
          l_cfrNd.text() = i_errorNorms[l_norm][l_q][l_cfr];
        }
      }
    }

    if( !l_doc.save_file(m_file.c_str()) ) EDGE_LOG_ERROR << "failed to write error norms to xml: " << m_file;
  }
}

void edge::io::ErrorNorms::write( const double i_errorNorms[3][N_QUANTITIES][N_CRUNS] ) {
  if( m_outType == sout || m_outType == sout_file ) print( i_errorNorms );
  if( m_outType == file || m_outType == sout_file ) writeXml( i_errorNorms );
}
