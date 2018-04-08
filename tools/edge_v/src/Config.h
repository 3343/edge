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
   * @brief Converts the given string vector to a double array.
   *
   * @param i_sep seperator of the string vector.
   * @param i_string string, which is converted.
   * @param o_val will be set to converted vector.
   */
  void vecStringToDouble( char         i_sep,
                          std::string &i_string,
                          double*      o_val );

  public:
    std::string                m_antnCfgFn;

    std::string                m_ucvmCfgFn;
    std::string                m_ucvmModelList;
    ucvm_ctype_t               m_ucvmCmode;
    int                        m_ucvmType;

    real                       m_minVp;
    real                       m_minVs;
    real                       m_minVs2;
    real                       m_maxVpVsRatio;
    real                       m_elmtsPerWave;

    std::string                m_meshFn;

    //! trafo, applied to the nodes before querying the UCVM
    real                       m_trafo[3][3];

    //! projection used for generation of the mesh
    std::string                m_projMesh;
    //! projection of the velocity model
    std::string                m_projVel;

    int                        m_tetRefinement;
    std::vector< std::string > m_faultInputFns;
    std::string                m_vmNodeFn;
    std::string                m_vmElmtFn;
    std::string                m_h5mFn;
    std::string                m_posFn;

    unsigned int               m_parallelMode;
    unsigned int               m_numWorker;

    /**
     * @brief Initializes the configuration.
     * 
     * @param i_pathToFile path to configuration file.
     */
    Config( const std::string &i_pathToFile );
};

#endif