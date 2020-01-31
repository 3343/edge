/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (breuer AT mytum.de)
 * @author David Lenz (dlenz AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2020, Alexander Breuer
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
 * Configuration class for EDGEcut.
 **/
#ifndef EDGE_CUT_CONFIG_H
#define EDGE_CUT_CONFIG_H

#include <string>
#include "logging.hpp"

namespace edge_cut {
  namespace io {
    class Config;
  }
}

class edge_cut::io::Config {
  private:
    //! input grid
    std::string m_gridIn;

    //! left output mesh
    std::string m_meshOutLeft;

    //! right output mesh
    std::string m_meshOutRight;

    //! front output mesh
    std::string m_meshOutFront;

    //! back output mesh
    std::string m_meshOutBack;

    //! bottom output mesh
    std::string m_meshOutBottom;

    //! top output mesh
    std::string m_meshOutTop;

    //! target z of the extrude
    double m_extrudeZ = 0;

  public:
    /**
     * Constructor.
     *
     * @param i_xmlPath path to the XML config.
     **/
    Config( std::string i_xmlPath );

    /**
     * Prints the config.
     **/
    void print();

    /**
     * Gets the input grid.
     *
     * @return path to the input grid.
     **/
    std::string getGridIn() { return m_gridIn; }

    /**
     * Gets the left output mesh.
     *
     * @return path to the left output mesh.
     **/ 
    std::string getMeshOutLeft() { return m_meshOutLeft; }

    /**
     * Gets the right output mesh.
     *
     * @return path to the right output mesh.
     **/ 
    std::string getMeshOutRight() { return m_meshOutRight; }

    /**
     * Gets the front output mesh.
     *
     * @return path to the front output mesh.
     **/ 
    std::string getMeshOutFront() { return m_meshOutFront; }

    /**
     * Gets the back output mesh.
     *
     * @return path to the back output mesh.
     **/ 
    std::string getMeshOutBack() { return m_meshOutBack; }

    /**
     * Gets the bottom output mesh.
     *
     * @return path to the bottom output mesh.
     **/ 
    std::string getMeshOutBottom() { return m_meshOutBottom; }

    /**
     * Gets the top output mesh.
     *
     * @return path to the top output mesh.
     **/ 
    std::string getMeshOutTop() { return m_meshOutTop; }

    /**
     * Gets the target z of the extrude operation.
     *
     * @return target z.
     **/
    double getExtZ() { return m_extrudeZ; }
};

#endif
