/**
 * @file This file is part of EDGE.
 *
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
 * Fault model.
 **/
#ifndef EDGE_V_FAULT_MODEL_H
#define EDGE_V_FAULT_MODEL_H

#include <map>
#include <vector>
#include "moab/Interface.hpp"
#include "constants.h"
#include "io/Config.h"
#include "Geometry.h"

namespace edge_v {
  namespace fault {
    class Model;
  }
}

// *** Fault Model ***
class Triangle {
public:
  xyz_point_t m_v1, m_v2, m_v3;

  Triangle() = default;
  Triangle( const xyz_point_t &i_v1,
            const xyz_point_t &i_v2,
            const xyz_point_t &i_v3 ) : m_v1( i_v1 ), m_v2( i_v2 ), m_v3( i_v3 ) {}
  Triangle( const real * );

  xyz_point_t centroid();
  xyz_point_t baryToPhy( const xyz_point_t & );
  std::vector< Triangle > subdivide( const int & );
};

typedef std::map< std::string, std::vector< real > > FDatum;

class edge_v::fault::Model {
  public:
    //! Number of cells in x,y directions
    int m_Nx, m_Ny;

    //! Min and Max physical coordinates
    real m_xMin, m_xMax, m_yMin, m_yMax;

    //! Reciprocal of the cell width in x and y directions
    real m_xScaleInv, m_yScaleInv;

    //! Angle fault makes with pos x-axis, increasing toward pos z-axis
    real m_faultAngle;

    // List of quantities being described by the input file
    std::vector< std::string > m_qtyNames;

    //! File names to read from; Number of fused runs is inferred from here
    std::vector< std::string > m_filenames;

    //! Main data array for the model, one FDatum per grid point,
    //! datum at (nx, ny) corresponds to m_faultData[ nx + ny*(N_x+1) ]
    FDatum *m_faultData;

    Model( const std::vector< std::string > & );
    ~Model();

    //! No copy constructor or copy assignment
    Model( const Model & ) = delete;
    Model & operator=( const Model & ) = delete;

    void populate();
    xyz_point_t projToPlane( const xyz_point_t & );
    int getNearestNx( const xyz_point_t & );
    int getNearestNy( const xyz_point_t & );
    FDatum getDatumPt( const xyz_point_t & );
    FDatum getDatumTri( Triangle & );
};

moab::Range getFaultFaces( moab_mesh & );
int faultAntn( const edge_v::io::Config & i_cfg,
                     moab_mesh & );

#endif