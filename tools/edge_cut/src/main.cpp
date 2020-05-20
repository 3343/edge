/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (breuer AT mytum.de)
 * @author David Lenz (dlenz AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2020, Alexander Breuer
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
 * This is the main file of EDGEcut.
 **/
#include "io/logging.hpp"
INITIALIZE_EASYLOGGINGPP
#include "io/OptionParser.h"
#include "io/Config.h"
#include "mesh/Extrude.h"

int main( int i_argc, char *i_argv[] ) {
  EDGE_CUT_LOG_INFO << "##########################################################################";
  EDGE_CUT_LOG_INFO << "##############   ##############            ###############  ##############";
  EDGE_CUT_LOG_INFO << "##############   ###############         ################   ##############";
  EDGE_CUT_LOG_INFO << "#####            #####       #####      ######                       #####";
  EDGE_CUT_LOG_INFO << "#####            #####        #####    #####                         #####";
  EDGE_CUT_LOG_INFO << "#############    #####         #####  #####                  #############";
  EDGE_CUT_LOG_INFO << "#############    #####         #####  #####      #########   #############";
  EDGE_CUT_LOG_INFO << "#####            #####         #####  #####      #########           #####";
  EDGE_CUT_LOG_INFO << "#####            #####        #####    #####        ######           #####";
  EDGE_CUT_LOG_INFO << "#####            #####       #####      #####       #####            #####";
  EDGE_CUT_LOG_INFO << "###############  ###############         ###############   ###############";
  EDGE_CUT_LOG_INFO << "###############  ##############           #############    ###############";
  EDGE_CUT_LOG_INFO << "#######################################################################cut";
  EDGE_CUT_LOG_INFO << "";
  EDGE_CUT_LOG_INFO << "              EDGEcut is available from: https://dial3343.org";
  EDGE_CUT_LOG_INFO << "";

  edge_cut::io::OptionParser l_options( i_argc, i_argv );

  EDGE_CUT_LOG_INFO << "EDGE version: " << PP_EDGE_VERSION;
  EDGE_CUT_LOG_INFO << "parsing runtime config";
  edge_cut::io::Config l_config( l_options.getXml() );
  l_config.print();

  // read the input
  EDGE_CUT_LOG_INFO << "reading the input grid";
  std::ifstream l_read( l_config.getGridIn(),
                        std::ios::in );

  std::vector< double > l_veCrds;
  double l_crd;
  while( l_read >> l_crd ) {
    l_veCrds.push_back( l_crd );
  }
  EDGE_CUT_CHECK_EQ( l_veCrds.size()%3, 0 );

  EDGE_CUT_LOG_INFO << "generating extruded surfaces";
  edge_cut::mesh::Extrude l_ext( l_veCrds.size()/3,
                                 (double (*)[3]) l_veCrds.data(),
                                 l_config.getExtZ(),
                                 l_config.getExtLvls() );

  // open file streams
  std::ofstream l_left(   l_config.getMeshOutLeft(),
                          std::ios::out );
  std::ofstream l_right(  l_config.getMeshOutRight(),
                          std::ios::out );
  std::ofstream l_front(  l_config.getMeshOutFront(),
                          std::ios::out );
  std::ofstream l_back(   l_config.getMeshOutBack(),
                          std::ios::out );
  std::ofstream l_bottom( l_config.getMeshOutBottom(),
                          std::ios::out );
  std::ofstream l_top(    l_config.getMeshOutTop(),
                          std::ios::out );

  // write
  EDGE_CUT_LOG_INFO << "writing output meshes";
  l_ext.writeOff( l_left,
                  l_right,
                  l_front,
                  l_back,
                  l_bottom,
                  l_top );

  EDGE_CUT_LOG_INFO << "thank you for using EDGEcut!";

  return EXIT_SUCCESS;
}