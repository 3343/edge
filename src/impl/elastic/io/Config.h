/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2019, Alexander Breuer
 * Copyright (c) 2016-2018, Regents of the University of California
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
 * Runtime configuration of elastics.
 **/
#ifndef EDGE_SEISMIC_CONFIG_H
#define EDGE_SEISMIC_CONFIG_H

#include "constants.hpp"
#include <vector>
#include <array>
#include "submodules/include/pugixml.hpp"
#include "linalg/HalfSpace.hpp"
#include "linalg/Domain.hpp"

namespace edge {
  namespace elastic {
    namespace io {
      class Config;
    }
  }
}

class edge::elastic::io::Config {
    /**
     * Prints the configuration.
     **/
    void print();

  public:
    //! files holding point source descriptions
    std::vector< std::string > m_ptSrcs;

    //! domains of the velocity model
    std::vector< linalg::Domain< real_mesh, N_DIM, edge::linalg::HalfSpace > > m_velDoms;

    //! values in the boxed velocity model
    std::vector< std::array< real_base, 5 > > m_velVals;

    // attenuation: central frequency and frequency ratio
    double m_attFreqs[2];

    //! friction law
    std::string m_frictionLaw = "";

    //! fault coordinate system: [0]: normal, [1]: along-strike, [2]: along-dip
    real_mesh m_faultCrds[N_DIM][N_DIM];

    //! linear slip weakening parameters
    real_base m_lsw[N_CRUNS][3];

    //! domains for rupture initialization
    std::vector< linalg::Domain< real_mesh, N_DIM, edge::linalg::HalfSpace > > m_rupDoms[N_CRUNS];

    //! initial values in the rupture domains
    std::vector< std::array< real_base, 1 + (N_DIM-1) > > m_stressInit[N_CRUNS];

    /**
     * Constructor of the config.
     *
     * @param i_xmlPath path to xml file.
     **/
    Config( const pugi::xml_document &i_xml );
};

#endif
