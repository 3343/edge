/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2015-2016, Regents of the University of California
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
 * Parser for command line options.
 **/

#include "OptionParser.h"
#include <vector>
#include "logging.h"

edge::io::OptionParser::OptionParser( int i_argc, char ** i_argv ) {
  /*
   * parse command line arguments
   */
  enum  optionIndex { UNKNOWN,
                      HELP,
                      XML,
                      VERBOSE };
  const option::Descriptor l_usage[] = {
    { UNKNOWN,             0, "",  ""    ,    Arg::Unknown,      "USAGE: ./edge [options]\n\n Options:" },
    { HELP,                0, "h", "help",    option::Arg::None, "  --help, -h  \tPrint usage and exit." },
    { XML,                 0, "x", "xml" ,    Arg::NonEmpty,     "  --xml, -x  \tLocation of the XML configuration." },
    { VERBOSE,             0, "v", "verbose", option::Arg::None, "  --vebose, -v \tActivates maximum logging verbosity." },
    { 0,0,0,0,0,0 }
  };

  // parse easylogging args first
  START_EASYLOGGINGPP( i_argc, i_argv);

  // ignore program name
  i_argc -= ( i_argc>0 );
  i_argv += ( i_argc>0 );

  option::Stats l_stats( l_usage, i_argc, i_argv );

  std::vector<option::Option> l_options;
                              l_options.resize( l_stats.options_max );
  std::vector<option::Option> l_buffer;
                              l_buffer.resize(  l_stats.buffer_max  );

  option::Parser l_parse(  l_usage,
                           i_argc,
                           i_argv,
                          &l_options[0],
                          &l_buffer[0]);

  if( l_parse.error() ) {
    exit( EXIT_FAILURE );
  }

  // exit for help or no options
  if( l_options[HELP] || i_argc == 0 || l_parse.optionsCount() == 0 ) {
    option::printUsage( std::cout, l_usage);
    exit( EXIT_SUCCESS );
  }

  m_xmlPath = l_options[XML].arg;

  EDGE_LOG_INFO << "  xml: " << m_xmlPath;
}
