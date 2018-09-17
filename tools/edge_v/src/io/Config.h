/**
 * @file This file is part of EDGE.
 *
 * @author Junyi Qiu (juq005 AT ucsd.edu)
 * @author Rajdeep Konwar (rkonwar AT ucsd.edu)
 * @author David Lenz (dlenz AT ucsd.edu)
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
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
 * Runtime configuration.
 **/
#ifndef EDGE_V_IO_CONFIG_H
#define EDGE_V_IO_CONFIG_H

#include <string>
extern "C" {
#include "ucvm.h"
}
#include "constants.h"
#include "Geometry.h"

namespace edge_v {
  namespace io {
    class Config;
  }
}

/**
 * @brief Runtime configuration.
 */
class edge_v::io::Config {
  /**
   * @brief Converts the given string vector to a real array.
   *
   * @param i_sep seperator of the string vector.
   * @param i_string string, which is converted.
   * @param o_val will be set to converted vector.
   *
   * @paramt TL_T_REAL floating point type.
   */
  template< typename TL_T_REAL >
  void vecStringToReal( char         i_sep,
                        std::string &i_string,
                        TL_T_REAL*   o_val ) {
    std::string l_right = i_string;

    unsigned short l_pos = 0;
    while( l_right.size() > 0 ) {
      std::size_t l_split = l_right.find_first_of(i_sep);
      if( l_split < l_right.size() ) {
        std::string l_left =  l_right.substr( 0, l_split );
                    l_right = l_right.substr( l_split+1  );

        o_val[l_pos] = atof(l_left.c_str());

        l_pos++;
      }
      else {
        o_val[l_pos] = atof(l_right.c_str());
        break;
      }
    }
  }

  public:
    std::string                m_antnCfgFn;

    std::string                m_ucvmCfgFn;
    std::string                m_ucvmModelList;
    std::string                m_ucvmCmode;
    std::string                m_ucvmType;

    std::string                m_velRule;
    std::string                m_meshFn;

    //! center of the refinement in the xy-plane
    float                      m_refCenter[2] = {0, 0};

    //! inner radius and out radius for the linear transition between the two char lengths
    float                      m_refRadii[2] = {0, std::numeric_limits<float>::max() };

    //! characteristic lengths, 0: inside the inner radius, 1: outside the outer radius
    float                      m_refCls[2] = {1, 1};

    //! trafo, applied to the nodes before querying the UCVM
    real                       m_trafo[3][3];

    //! projection used for generation of the mesh
    std::string                m_projMesh;
    //! projection of the velocity model
    std::string                m_projVel;

    std::vector< std::string > m_faultInputFns;
    std::string                m_annoFn;
    std::string                m_posFn;

    int                        m_tetRefinement = 0;

    /**
     * @brief Initializes the configuration.
     * 
     * @param i_pathToFile path to configuration file.
     */
    Config( const std::string &i_pathToFile );
};

#endif