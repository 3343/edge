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
 * Output for error norms.
 **/

#ifndef EDGE_IO_ERROR_NORMS_H
#define EDGE_IO_ERROR_NORMS_H

#include <string>
#include "constants.hpp"

namespace edge {
  namespace io {
    class ErrorNorms;
  }
}

class edge::io::ErrorNorms {
  private:
    /**
     * Prints the error norms to out stream.
     *
     * @param i_errorNorms error norms. [0][][]: L1, [1][][]: L2, [2][][]: Linf; [][*][]: quantity; [][][*]: cfr
     **/
    void print( const double i_errorNorms[N_CRUNS][3][N_QUANTITIES] );

    /**
     * Writes the error norms to disk.
     *
     * @param i_errorNorms error norms. [0][][]: L1, [1][][]: L2, [2][][]: Linf; [][*][]: quantity; [][][*]: cfr
     **/
    void writeXml( const double i_errorNorms[N_CRUNS][3][N_QUANTITIES] );

  public:
    enum OutType {
      none,
      sout,
      file,
      sout_file
    };

    // output type
    OutType m_outType;

    // path to file
    std::string m_file;

    /**
     * Constructor
     *
     * @param i_outType output type.
     * @param i_file path to file for output.
     **/
    ErrorNorms( OutType i_outType, std::string i_file="" ): m_outType(i_outType), m_file(i_file) {};

    /**
     * Constructor
     *
     * @param i_outType output type.
     * @param i_file path to file for output.
     **/
    ErrorNorms( std::string i_outType, std::string i_file = "" ): m_file(i_file) {
      if( i_outType == "sout" )           m_outType = sout;
      else if( i_outType == "file" )      m_outType = file;
      else if( i_outType == "sout_file" ) m_outType = sout_file; 
      else                                m_outType = none;
    }

    /**
     * Checks if error-norm output is enabled.
     *
     * @return true if enabled, false otherwise.
     **/
    bool outEnabled() { return !(m_outType==none); };

    /**
     * Writes the error norms.
     *
     * @param i_errorNorms error norms. [0][][]: L1, [1][][]: L2, [2][][]: Linf; [][*][]: quantity; [][][*]: cfr
     **/
    void write( const double i_errorNorms[N_CRUNS][3][N_QUANTITIES] );
};

#endif
