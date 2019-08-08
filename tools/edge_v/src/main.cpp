/**
 * @file This file is part of EDGE.
 *
 * @author Rajdeep Konwar (rkonwar AT ucsd.edu)
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2019, Alexander Breuer
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
 * Entry point of standalone EDGE-V.
 **/
#include "io/Config.h"
#include "io/Csv.h"
#include "models/Constant.h"
#include "mesh/Mesh.h"
#include "time/Cfl.h"
#include "time/Groups.h"

#include "io/logging.h"
#ifdef PP_USE_EASYLOGGING
INITIALIZE_EASYLOGGINGPP
#endif

int main( int i_argc, char *i_argv[] ) {
#ifdef PP_USE_EASYLOGGING
  START_EASYLOGGINGPP( i_argc, i_argv );
#endif

  EDGE_V_LOG_INFO << "VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV";
  EDGE_V_LOG_INFO << "VVVVVVVVVVVVVV   VVVVVVVVVVVVVV            VVVVVVVVVVVVVVV  VVVVVVVVVVVVVV";
  EDGE_V_LOG_INFO << "VVVVVVVVVVVVVV   VVVVVVVVVVVVVVV         VVVVVVVVVVVVVVVV   VVVVVVVVVVVVVV";
  EDGE_V_LOG_INFO << "VVVVV            VVVVV       VVVVV      VVVVVV                       VVVVV";
  EDGE_V_LOG_INFO << "VVVVV            VVVVV        VVVVV    VVVVV                         VVVVV";
  EDGE_V_LOG_INFO << "VVVVVVVVVVVVV    VVVVV         VVVVV  VVVVV                  VVVVVVVVVVVVV";
  EDGE_V_LOG_INFO << "VVVVVVVVVVVVV    VVVVV         VVVVV  VVVVV      VVVVVVVVV   VVVVVVVVVVVVV";
  EDGE_V_LOG_INFO << "VVVVV            VVVVV         VVVVV  VVVVV      VVVVVVVVV           VVVVV";
  EDGE_V_LOG_INFO << "VVVVV            VVVVV        VVVVV    VVVVV        VVVVVV           VVVVV";
  EDGE_V_LOG_INFO << "VVVVV            VVVVV       VVVVV      VVVVV       VVVVV            VVVVV";
  EDGE_V_LOG_INFO << "VVVVVVVVVVVVVVV  VVVVVVVVVVVVVVV         VVVVVVVVVVVVVVV   VVVVVVVVVVVVVVV";
  EDGE_V_LOG_INFO << "VVVVVVVVVVVVVVV  VVVVVVVVVVVVVV           VVVVVVVVVVVVV    VVVVVVVVVVVVVVV";
  EDGE_V_LOG_INFO << "VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV";
  EDGE_V_LOG_INFO << "";
  EDGE_V_LOG_INFO << "              EDGE-V is available from: https://dial3343.org";
  EDGE_V_LOG_INFO << "";

  // parse argumenets
  std::vector< std::string > l_args( i_argv, i_argv + i_argc );
  if( l_args.size() != 3 || (l_args[1] != "-x" && l_args[1] != "--xml") ) {
    EDGE_V_LOG_INFO << "USAGE: " << l_args[0] << " [options]";
    EDGE_V_LOG_INFO;
    EDGE_V_LOG_INFO << "  Options:";
    EDGE_V_LOG_INFO << "   --help, -h Print usage and exit.";
    EDGE_V_LOG_INFO << "   --xml,  -x Location of the XML configuration.";

    for( std::size_t l_ar = 0; l_ar < l_args.size(); l_ar++ ) {
      if( l_args[l_ar] == "-h" || l_args[l_ar] == "--help" ) {
        return EXIT_SUCCESS;
      }
    }

    return EXIT_FAILURE;
  }

  EDGE_V_LOG_INFO << "parsing xml config";
  std::string l_xml = l_args[2];
  edge_v::io::Config l_config( l_xml );

  EDGE_V_LOG_INFO << "initializing MOAB";
  std::string l_meshIn = l_config.getMeshIn();
  edge_v::io::Moab l_moab( l_meshIn );
  l_moab.printStats();

  // initializing and setting mesh data
  EDGE_V_LOG_INFO << "initializing mesh interface";
  edge_v::mesh::Mesh l_mesh( l_moab );
  l_mesh.printStats();

  EDGE_V_LOG_INFO << "computing CFL time steps";
  edge_v::models::Constant l_veMod( 2 );
  edge_v::time::Cfl l_cfl( l_mesh.getElTy(),
                           l_mesh.nVes(),
                           l_mesh.nEls(),
                           l_mesh.getElVe(),
                           l_mesh.getVeCrds(),
                           l_mesh.getInDiaEl(),
                           l_veMod );
  l_cfl.printStats();

  if( l_config.getTsOut() != "" ) {
    EDGE_V_LOG_INFO << "writing time steps";
    std::string l_colNames = "ts_cfl";
    double const * l_data[1] = { l_cfl.getTimeSteps() };
    edge_v::io::Csv::write( l_config.getTsOut(),
                            1,
                            l_mesh.nEls(),
                            &l_colNames,
                            5,
                            l_data );
  }

  EDGE_V_LOG_INFO << "computing time step groups";
  unsigned short l_nTsGroups = l_config.nTsGroups();
  double *l_rates = new double[l_nTsGroups];
  for( unsigned short l_tg = 0; l_tg < l_nTsGroups; l_tg++ ) {
    l_rates[l_tg] = 2.0;
  }

  double l_funDt = l_config.getFunDt();
  // search for fundamental dt, if not specified
  double l_speedUp = 0;
  if( l_funDt == 0 ) {
    double l_dt = 1.0;
    while( l_dt > 0.5 ) {
      edge_v::time::Groups l_tsGroups(  l_mesh.getElTy(),
                                        l_mesh.nEls(),
                                        l_mesh.getElFaEl(),
                                        l_nTsGroups,
                                        l_rates,
                                        l_dt,
                                        l_cfl.getTimeSteps() );
      if( l_tsGroups.getSpeedUp() > l_speedUp ) {
        l_speedUp = l_tsGroups.getSpeedUp();
        l_funDt = l_dt;
        EDGE_V_LOG_INFO << "  found new best fundamental dt / speedup: " << l_funDt << " / " << l_speedUp;
      }
      l_dt -= 0.01;
    }
  }

  edge_v::time::Groups l_tsGroups(  l_mesh.getElTy(),
                                    l_mesh.nEls(),
                                    l_mesh.getElFaEl(),
                                    l_nTsGroups,
                                    l_rates,
                                    l_funDt,
                                    l_cfl.getTimeSteps() );
  delete[] l_rates;
  l_tsGroups.printStats();

  EDGE_V_LOG_INFO << "storing elements' time step groups";
  std::string l_tagElTg = "edge_v_element_time_groups";
  l_moab.deleteTag( l_tagElTg );
  l_moab.setEnData( l_mesh.getElTy(),
                    l_tagElTg,
                    l_tsGroups.getElTg() );

  EDGE_V_LOG_INFO << "storing relative time steps of the groups";
  std::string l_tagRelTs = "edge_v_relative_time_steps";
  l_moab.deleteTag( l_tagRelTs );
  l_moab.setGlobalData( l_tagRelTs,
                        l_tsGroups.nGroups()+1,
                        l_tsGroups.getTsIntervals() );

  EDGE_V_LOG_INFO << "reordering by time step groups";
  l_moab.reorder( l_mesh.getElTy(),
                  l_tagElTg );

  if( l_config.getMeshOut() != "" ) {
    EDGE_V_LOG_INFO << "writing mesh";
    l_moab.writeMesh( l_config.getMeshOut() );
  }

  return EXIT_SUCCESS;
}