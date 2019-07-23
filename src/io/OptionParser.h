/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2015-2018, Regents of the University of California
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

#ifndef EDGE_IO_OPTION_PARSER_H_
#define EDGE_IO_OPTION_PARSER_H_

#include <string>

#pragma GCC diagnostic push
#if !defined(__clang__)
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
#include "submodules/optionparser/optionparser.h"
#pragma GCC diagnostic pop

#include "logging.h"

namespace edge {
  namespace io {
    class OptionParser;
  }
}

class edge::io::OptionParser {
  //private:
    struct Arg: public option::Arg {
      static option::ArgStatus Unknown(const option::Option& i_option,
                                             bool            i_msg) {
        if (i_msg) EDGE_LOG_ERROR << "Unknown option '" << i_option.name << "'\n";
        return option::ARG_ILLEGAL;
      }

      static option::ArgStatus NonEmpty( const option::Option& i_option,
                                               bool            i_msg ) {
        if (i_option.arg != 0 && i_option.arg[0] != 0)
          return option::ARG_OK;

        if (i_msg) EDGE_LOG_ERROR << "Option '" << i_option.name << "' requires a non-empty argument\n";
        return option::ARG_ILLEGAL;
      }
    };

    std::string m_xmlPath;

  public:
    /**
     * Constructor of the option parser.
     **/
    OptionParser( int i_argc, char **i_argv );

    /**
     * Gets the path to the XML-file.
     *
     * @return path to xml-file.
     **/
    std::string getXmlPath(){ return m_xmlPath; };
};

#endif
