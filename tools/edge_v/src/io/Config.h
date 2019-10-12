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

    //! base name for the mesh output by partition
    std::string m_meshOutPa[2] = {"", ""};

    std::string m_seismicExpr = "";

    //! writes element annotations if true
    bool m_writeElAn = false;

    //! periodic boundary conditions
    bool m_periodic = false;

    //! number of partitions to derive
    std::size_t m_nPartitions = 0;

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
     * Gets the base name of the by partition mesh output.
     *
     * @return base name of the mesh output.
     **/
    std::string const & getMeshOutPaBase() const { return m_meshOutPa[0]; }

    /**
     * Gets the file extension of the by partition mesh output.
     *
     * @return file extension the mesh output.
     **/
    std::string const & getMeshOutPaExt() const { return m_meshOutPa[1]; }

    /**
     * Gets the expression string for a seismic expression-based velocity model.
     *
     * @return expression string.
     **/
    std::string const & getVelModSeismicExpr() const { return m_seismicExpr; }

    /**
     * Gets the configuration of periodic boundary conditions.
     *
     * @return true if the mesh has periodic boundaries, false if not.
     **/
    bool getPeriodic() const { return m_periodic; }

    /**
     * Gets the configuration for element annotations in the output mesh.
     *
     * @return true if the mesh's elements are annotated.
     **/
    bool getWriteElAn() const { return m_writeElAn; }

    /**
     * Gets the number of partitions in the output mesh.
     *
     * @return number of partitions.
     **/
    std::size_t nPartitions() const { return m_nPartitions; }

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