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
 * EDGE-V config.
 **/
#ifndef EDGE_V_IO_CONFIG_H
#define EDGE_V_IO_CONFIG_H

#include <vector>
#include <string>

namespace edge_v {
  namespace io {
    class Config;
  }
}

/**
 * EDGE-V config.
 **/
class edge_v::io::Config {
  private:
    //! path to the input mesh
    std::string m_meshIn = "";

    //! path to the output mesh
    std::string m_meshOut = "";

    //! path to the output-csv for the time steps
    std::string m_tsOut = "";

    //! number of time step groups
    unsigned short m_nTsGroups = 0;

    //! relative minimum time step
    double m_funDt = 0;

  public:
    /**
     * Constructor.
     *
     * @param i_xml xml file, which is parsed.
     **/
    Config( std::string & i_xml );

    /**
     * Gets the input mesh.
     *
     * @return input mesh.
     **/
    std::string const & getMeshIn() const { return m_meshIn; }

    /**
     * Gets the output mesh.
     *
     * @return output mesh.
     **/
    std::string const & getMeshOut() const { return m_meshOut; }

    /**
     * Gets the the number of time step groups.
     *
     * @return number of time step groups.
     **/
    unsigned short nTsGroups() const { return m_nTsGroups; }

    /**
     * Gets the relative fundamental time step.
     *
     * @return relative fundamental time step.
     **/
    double getFunDt() const { return m_funDt; }

    /**
     * Gets the output file for the time steps of the elements.
     *
     * @return output file for time steps.
     **/
    std::string const & getTsOut() const { return m_tsOut; }
};

#endif