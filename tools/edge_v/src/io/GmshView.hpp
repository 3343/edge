/**
 * @file This file is part of EDGE.
 *
 * @author Rajdeep Konwar (rkonwar AT ucsd.edu)
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
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
 * Writes a GMSH view, based on the given coordinates and values.
 **/
#ifndef EDGEV_IO_GMSHVIEW_HPP
#define EDGEV_IO_GMSHVIEW_HPP

#include <string>
#include <fstream>
#include <sstream>

namespace edge_v {
  namespace io {
    class GmshView;
  }
}

/**
 * @brief Inteface to Gmsh's view.
 */
class edge_v::io::GmshView {
  public:
    /**
     * @brief Writes the Gmsh-view.
     * 
     * @param i_path path of the output file.
     * @param i_elType element type: 'tria3' or 'tet4'
     * @param i_nEls number of elements.
     * @param i_elVe vertices adjacent to the elements.
     * @param i_veCrds coordinates of the vertices.
     * @param i_values values, which are assigned to the nodes.
     *
     * @paramt TL_T_REAL floating point precision of the values.
     */
    template< typename TL_T_LID,
              typename TL_T_REAL >
    static void write( std::string const  & i_path,
                       std::string const  & i_elType,
                       std::size_t          i_nEls,
                       TL_T_LID          (* i_elVe),
                       double      const (* i_veCrds)[3],
                       TL_T_REAL   const  * i_values ) {
      // open file stream
      std::ofstream l_stream( i_path,
                              std::ios::out );

      // check that the stream is point
      if( !l_stream.is_open() ) {
        std::cerr << "Error: Failed generating output file for Gmsh-view." << std::endl;
        exit( EXIT_FAILURE );
      }

      unsigned int l_bufferSize = 10000;
      std::stringstream l_ss;

      // derive number of vertices
      unsigned short l_nVes = 0;
      if( i_elType == "tria3" ) {
        l_nVes = 3;
      }
      else if( i_elType == "tet4" ) {
        l_nVes = 4;
      }
      else assert( false );

      l_ss << "/******************************"  << std::endl;
      l_ss << " * EDGE-V generated Gmsh-view *"  << std::endl;
      l_ss << " ******************************/" << std::endl;
      l_ss << "View \"EDGE-V\" {" << std::endl;

      // iterate over elements
      for( std::size_t l_el = 0; l_el < i_nEls; l_el++ ) {
        // write to disk, if required
        if( l_el % l_bufferSize == 0 ) {
          l_stream << l_ss.rdbuf();
          l_stream.flush();

          l_ss.str( std::string() );
          l_ss.clear();
        }

        // vertices of the element
        TL_T_LID *l_ves = i_elVe + l_el * l_nVes;

        // open a element type point
        if( i_elType == "tria3" )
          l_ss << "ST(";
        else if( i_elType == "tet4" )
          l_ss << "SS("; 

        // write the vertices' coordinates
        for( unsigned short l_ve = 0; l_ve < l_nVes; l_ve++ ) {
          l_ss << i_veCrds[ l_ves[l_ve] ][0] << ","
               << i_veCrds[ l_ves[l_ve] ][1] << ","
               << i_veCrds[ l_ves[l_ve] ][2];
          if( l_ve != l_nVes-1) l_ss << ",";
        }

        // write the vertices values
        l_ss << "){";
        for( unsigned short l_ve = 0; l_ve < l_nVes; l_ve++ ) {
              l_ss << i_values[ l_ves[l_ve] ];
              if( l_ve != l_nVes-1) l_ss << ",";
        }
        l_ss << "};" << std::endl;
      }

      l_ss << "};" << std::endl;

      // final buffer write to file
      l_stream << l_ss.rdbuf();
      l_stream.flush();
      l_stream.close();
    }
};

#endif