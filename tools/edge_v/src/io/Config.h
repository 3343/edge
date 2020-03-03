/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2019-2020, Alexander Breuer
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
#include <limits>

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

    //! base name, mesh extension, meta data extension
    std::string m_meshOut[2] = {"", ""};

    //! seismic expression
    std::string m_seismicExpr = "";

    //! ucvm parameters
    struct {
      double trafoSrc[3][3] = { {1, 0, 0},
                                {0, 1, 0},
                                {0, 0, 1} };
      std::string projSrc = "";
      std::string projDes = "";

      std::string models = "";
      std::string modelType = "";
      std::string crdMode = "";
      std::string rule = "";
      bool lowerToSurf = false;
    } m_ucvm;

    //! mesh refinement parameters
    struct {
      std::string expr = "";
      std::string out = "";
    } m_ref;

    //! writes element annotations if true
    bool m_writeElAn = false;

    //! periodic boundary conditions
    int m_periodic = std::numeric_limits< int >::max();

    //! number of partitions to derive
    std::size_t m_nPartitions = 1;

    //! path to the output-csv for the time steps
    std::string m_tsOut = "";

    //! number of time step groups
    unsigned short m_nTsGroups = 1;

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
     * Gets the base name of the mesh output.
     *
     * @return base name of the mesh output.
     **/
    std::string const & getMeshOutBase() const { return m_meshOut[0]; }

    /**
     * Gets the file extension of the by mesh output.
     *
     * @return file extension the mesh output.
     **/
    std::string const & getMeshOutExt() const { return m_meshOut[1]; }

    /**
     * Gets the expression string for a seismic expression-based velocity model.
     *
     * @return expression string.
     **/
    std::string const & getVelModSeismicExpr() const { return m_seismicExpr; }

    /**
     * Gets the initial transformation matrix, applied before the projections and UCVM-queries.
     *
     * @return transformation matrix.
     **/
    double const (* getVelModUcvmTrafoSrc() const) [3] { return m_ucvm.trafoSrc; }

    /**
     * Gets the proj source projection applied (in conjunction with the destination projection) before the UCVM-queries.
     *
     * @return source projection.
     **/
    std::string const & getVelModUcvmProjSrc() const { return m_ucvm.projSrc; }

    /**
     * Gets the proj destination projection applied (in conjunction with the source projection) before the UCVM-queries.
     *
     * @return destination projection.
     **/
    std::string const & getVelModUcvmProjDes() const { return m_ucvm.projDes; }

    /**
     * Gets the UCVM models.
     *
     * @return UCVM models.
     **/
    std::string const & getVelModUcvmModels() const { return m_ucvm.models; }

    /**
     * Gets the UCVM model type.
     *
     * @return UCVM model type.
     **/
    std::string const & getVelModUcvmModelType() const { return m_ucvm.modelType; }

    /**
     * Gets the UCVM coordinate mode.
     *
     * @return UCVM coordinate mode.
     **/
    std::string const & getVelModUcvmCrdMode() const { return m_ucvm.crdMode; }

    /**
     * Gets the normalization rule, applied to the UCVM output.
     *
     * @return normalization rule.
     **/
    std::string const & getVelModUcvmRule() const { return m_ucvm.rule; }

    /**
     * Gets the lower-to-surface option for the UCVM queries.
     *
     * @return lower-to-surface option.
     **/
    bool const & getLowerToSurf() const { return m_ucvm.lowerToSurf; }

    /**
     * Gets the expression of the mesh refinement.
     *
     * @return expression for mesh refinement.
     **/
    std::string const & getRefExpr() const { return m_ref.expr; }

    /**
     * Gets the output path for the mesh refinement.
     *
     * @return path to mesh refinement output.
     **/
    std::string const & getRefOut() const { return m_ref.out; }

    /**
     * Gets the configuration of periodic boundary conditions.
     *
     * @return true if the mesh has periodic boundaries, false if not.
     **/
    int getPeriodic() const { return m_periodic; }

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