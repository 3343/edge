/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2016-2017, Regents of the University of California
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
 * Runtime configuration of elastics.
 **/

#include "Config.h"
#include "io/Config.h"
#include "io/logging.h"
#include <string>

void edge::elastic::io::Config::print() {
  EDGE_LOG_INFO << "  printing implementation-specific config for elastics (if any)";

  // print info about the velocity model
  if( m_velDoms.size() > 0 ) {
    EDGE_LOG_INFO << "    found " << m_velDoms.size() << " velmodel-domains in the config: ";
    for( std::size_t l_do = 0; l_do < m_velDoms.size(); l_do++ ) {
      EDGE_LOG_INFO << "      domain #" << l_do << ":";
      EDGE_LOG_INFO << "        rho: " << m_velVals[l_do][0] << ", "
                    <<         "lam: " << m_velVals[l_do][1] << ", "
                    <<         "mu:  " << m_velVals[l_do][1] << ", ";
      std::vector< std::string > l_doStrs = m_velDoms[l_do].toString();
      for( std::size_t l_ob = 0; l_ob < l_doStrs.size(); l_ob++ ) {
        EDGE_LOG_INFO << "        object #" << l_ob << ": " << l_doStrs[l_ob];
      }
    }
  }

  // print info about kinematic sources terms
  if( m_kinSrcs.size() > 0 ) {
    EDGE_LOG_INFO << "    there are " << m_kinSrcs.size() << " kinematic source files in your config:";
    for( unsigned short l_ki = 0; l_ki < m_kinSrcs.size(); l_ki++ )
      EDGE_LOG_INFO << "      #" << l_ki << ": " << m_kinSrcs[l_ki];
#ifndef PP_HAS_NETCDF
    EDGE_LOG_FATAL << "netCDF is required for kinematic sources descriptions and not supported by your build.";
#endif
  }

  // print rupture info
  if( m_frictionLaw == "lsw" ) {
    EDGE_LOG_INFO << "    the config has spontaneous rupture setups:";
    EDGE_LOG_INFO << "      fault coordinate system:";
    for( unsigned short l_cd = 0; l_cd < N_DIM; l_cd++ ) {
      std::string l_dir = (l_cd == 0) ? "normal" :
                          (l_cd == 1) ? "strike" :
                                        "dip";

       EDGE_LOG_INFO << "        " << l_dir << ": "
                     <<  m_faultCrds[l_cd][0] << " "
                     <<  m_faultCrds[l_cd][1] << " "
                     << ((N_DIM > 2) ?  std::to_string(m_faultCrds[l_cd][2]) : "");
    }

    for( int_cfr l_ru = 0; l_ru < N_CRUNS; l_ru++ ) {
      EDGE_LOG_INFO << "      rupture setup #" << l_ru << ":";
      EDGE_LOG_INFO << "        mus: " << m_lsw[l_ru][0] << ", "
                    <<         "mud: " << m_lsw[l_ru][1] << ", "
                    <<          "dc: " << m_lsw[l_ru][2];
       for( std::size_t l_do = 0; l_do < m_rupDoms[l_ru].size(); l_do++ ) {
          EDGE_LOG_INFO << "        domain #" << l_do << ":";
          EDGE_LOG_INFO << "          sn0: " << m_stressInit[l_ru][l_do][0] << ", "
                        <<           "ss0: " << m_stressInit[l_ru][l_do][1] << " (normal-strike)"
                        << ((N_DIM==3) ? ", " + std::to_string(m_stressInit[l_ru][l_do][2]) + " (normal-dip)"
                                       : "" ) ;
          std::vector< std::string > l_doStrs = m_rupDoms[l_ru][l_do].toString();
          for( std::size_t l_ob = 0; l_ob < l_doStrs.size(); l_ob++ ) {
            EDGE_LOG_INFO << "          object #" << l_ob << ": " << l_doStrs[l_ob];
          }
       }
    }
  }
}

edge::elastic::io::Config::Config( const pugi::xml_document &i_xml ) {
  pugi::xml_node l_setups = i_xml.child("edge").child("cfr").child("setups");

  /*
   * read kinematic source info if available
   */
  for( pugi::xml_node l_ki = l_setups.child("kinematic_sources").child("file");
       l_ki;
       l_ki = l_ki.next_sibling("file") ) {
    m_kinSrcs.push_back( l_ki.text().as_string() );
  }

  /*
   * read velocity model, if available
   */
  pugi::xml_node l_velMod = i_xml.child("edge").child("cfr").child("velocity_model");

  for( pugi::xml_node l_do = l_velMod.child("domain"); l_do; l_do = l_do.next_sibling("domain") ) {
    // add this domain
    m_velDoms.resize( m_velDoms.size() + 1 );

    // initial velocity values of this domain
    std::vector< real_base > l_vals;

    // names of the values
    std::vector< std::vector< std::string > > l_valsN(3);
    l_valsN[0].push_back( "rho" );
    l_valsN[1].push_back( "lambda" );
    l_valsN[2].push_back( "mu" );

    // parse the domain
    edge::io::ConfigDoms< N_DIM >::parse( l_do,
                                          l_valsN,
                                          edge::io::ConfigDoms< N_DIM >::F64,
                                          l_vals,
                                          m_velDoms.back() );
    EDGE_CHECK_EQ( l_valsN.size(), 3 );

    // copy over velocity values
    m_velVals.resize( m_velVals.size() + 1 );
    for( unsigned short l_va = 0; l_va < 3; l_va++ ) {
      m_velVals.back()[l_va] = l_vals[l_va];
    }
  }

  /*
   * read rupture setups, if available
   */
  int_cfr l_ru = 0;

  // set fruction law to LSW by default
  if( l_setups.child("rupture") ) m_frictionLaw = "lsw";

  // set minus-to-plus vector
  for( unsigned short l_d0 = 0; l_d0 < N_DIM; l_d0++ ) {
    std::string l_d = (l_d0 == 0) ? "normal" :
                      (l_d0 == 1) ? "strike" :
                                    "dip";

    for( unsigned short l_d1 = 0; l_d1 < N_DIM; l_d1++ ) {
      std::string l_c = (l_d1 == 0) ? "x" :
                        (l_d1 == 1) ? "y" :
                                      "z";

      m_faultCrds[l_d0][l_d1] = l_setups.child("fault_coordinates").child(l_d.c_str()).child(l_c.c_str()).text().as_double();
    }
  }

  // iterate over rupture setups
  for( pugi::xml_node l_rs = l_setups.child("rupture"); l_rs; l_rs = l_rs.next_sibling("rupture") ) {
    // read friction parameters
    pugi::xml_node l_par = l_rs.child("friction_law");
    m_lsw[l_ru][0] = l_par.child("mus").text().as_double();
    m_lsw[l_ru][1] = l_par.child("mud").text().as_double();
    m_lsw[l_ru][2] = l_par.child("dc").text().as_double();

    // read initial values
    pugi::xml_node l_sInit = l_rs.child("stress_init");
    for( pugi::xml_node l_do = l_sInit.child("domain"); l_do; l_do = l_do.next_sibling("domain") ) {
      // add this domain
      m_rupDoms[l_ru].resize( m_rupDoms[l_ru].size() + 1 );

      // initial stress values of this domain
      std::vector< real_base > l_stressInit;

      // names of the data
      std::vector< std::vector< std::string > > l_dataN(2);
      l_dataN[0].push_back( "values" ); l_dataN[0].push_back( "sn0" );
      l_dataN[1].push_back( "values" ); l_dataN[1].push_back( "ss0" ); l_dataN[1].push_back( "normal_strike" );
      if( N_DIM == 3 ) {
        l_dataN.resize( 3 );
        l_dataN[2].push_back( "values" ); l_dataN[2].push_back( "ss0" ); l_dataN[2].push_back( "normal_dip" );
      }

      edge::io::ConfigDoms< N_DIM >::parse( l_do,
                                            l_dataN,
                                            edge::io::ConfigDoms< N_DIM >::F64,
                                            l_stressInit,
                                            m_rupDoms[l_ru].back() );
      EDGE_CHECK_EQ( l_stressInit.size(), 1 + (N_DIM-1) );

      // copy over initial stress values
      m_stressInit[l_ru].resize( m_stressInit[l_ru].size() + 1 );
      for( unsigned short l_st = 0; l_st < 1 + (N_DIM-1); l_st++ ) {
        m_stressInit[l_ru].back()[l_st] = l_stressInit[l_st];
      }
    }

    l_ru++;
    // abort if more setups than runs are available
    if( l_ru >= N_CRUNS ) break;
  }

  // print config
  print();
}
