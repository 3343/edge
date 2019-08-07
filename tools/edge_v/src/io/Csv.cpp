/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2019, Alexander Breuer
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
 * Writes data as CSV-files.
 **/
#include "Csv.h"

#include <fstream>
#include <iomanip> 

void edge_v::io::Csv::write( std::string    const         & i_csv,
                             unsigned short                 i_nCols,
                             std::size_t                    i_nRows,
                             std::string    const         * i_cols,
                             unsigned short                 i_precision,
                             double         const * const * i_data ) {
  // assemble header
  std::string l_header  = "# EDGE-V\n";
              l_header += "# code version: " + std::string(PP_EDGE_VERSION) + "\n";
              l_header += "# build date / time: " + std::string(__DATE__) + " / " + std::string(__TIME__) + "\n";
  for( unsigned short l_co = 0; l_co < i_nCols; l_co++ ) {
    if( l_co > 0 ) l_header += ",";
    l_header += i_cols[l_co];
  }
  l_header += '\n';

  // open the file
  std::ofstream l_csv;
  l_csv.open( i_csv );

  if( l_csv.is_open() ) {
    // write header
    l_csv << l_header;

    // set precision and write data
    l_csv << std::setprecision( i_precision );
    for( std::size_t l_ro = 0; l_ro < i_nRows; l_ro++ ) {
      for( unsigned short l_co = 0; l_co < i_nCols; l_co++ ) {
        if( l_co > 0 ) l_csv << ",";
        l_csv << i_data[l_co][l_ro];
      }
      l_csv << "\n";
    }
  }
}