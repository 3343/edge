/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2019-2020, Alexander Breuer
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
 * Ucvm velocity model.
 **/
#ifndef EDGE_V_MODELS_UCVM_H
#define EDGE_V_MODELS_UCVM_H

#include "io/Ucvm.h"
#include "../Model.h"

namespace edge_v {
  namespace models {
    namespace seismic {
      class Ucvm;
    }
  }
}

/**
 * EDGE-V velocity model, using UCVM.
 **/
class edge_v::models::seismic::Ucvm: public Model {
  private:
    //! number of points
    std::size_t m_nPts = 0;

    //! p-wave velocities
    float *m_velP = nullptr;

    //! s-wave velocities
    float *m_velS = nullptr;

    //! densities
    float *m_rho = nullptr;

    //! ucvm reader
    io::Ucvm & m_ucvmReader;

    //! trafo applied before calling proj
    double m_trafoSrc[3][3] = {0};

    //! proj source trafo
    std::string m_projSrc = "";

    //! proj destination trafo
    std::string m_projDes = "";

    //! ucvm type
    std::string m_ucvmType = "";

    /**
     * Frees the classe's memory.
     **/
    void free();

  public:
    /**
     * Constructor.
     *
     * @param i_ucvmReader UCVM reader.
     * @param i_trafoSrc trafo, applied to the source coordinates before the proj.4 projection.
     * @param i_projSrc proj.4 string of the source projection.
     * @param i_projDes proh.4 string of the destination projection (used to query UCVM).
     * @param i_ucvmType ucvm type: gtl, crust or cmb.
     */
    Ucvm( io::Ucvm          & i_ucvmReader,
          double      const   i_trafoSrc[3][3],
          std::string const & i_projSrc,
          std::string const & i_projDes,
          std::string const & i_ucvmType );

    /**
     * Destructor.
     **/
    ~Ucvm();

   /**
     * Inits the velocity model at the given points.
     *
     * @param i_nPts number of points.
     * @param i_pts coordinates of the points.
     **/
    void init( std::size_t          i_nPts,
               double      const (* i_pts)[3] );

    /**
     * Gets the minimum wave speed at a point.
     *
     * @param i_pt point at which the minimum wave speed is derived.
     **/
    double getMinSpeed( std::size_t i_pt ) const;

    /**
     * Gets the maximum wave speed at a point.
     *
     * @param i_pt point at which the maximum wave speed is derived.
     **/
    double getMaxSpeed( std::size_t i_pt ) const;

    /**
     * Gets the p-wave velocities.
     *
     * @return p-wave velocities.
     **/
    float const * getVelP() { return m_velP; }

    /**
     * Gets the s-wave velocities.
     *
     * @return s-wave velocities.
     **/
    float const * getVelS() { return m_velS; }

    /**
     * Gets the material densities.
     *
     * @return densities.
     **/
    float const * getRho() { return m_rho; }
};

#endif