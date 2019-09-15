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
 * Mesh-interface using EDGE-V.
 **/
#ifndef EDGE_MESH_EDGEV_H
#define EDGE_MESH_EDGEV_H

#include "../data/EntityLayout.type"
#include <edge_v/edge_v.h>

namespace edge {
  namespace mesh {
    class EdgeV;
  }
}

/**
 * Mesh-related functions and data.
 **/
class edge::mesh::EdgeV {
  private:
    //! EDGE-V moab interface
    edge_v::io::Moab m_moab;

    //! EDGE-V mesh interface
    edge_v::mesh::Mesh m_mesh;

    //! relative time steps, first is fundamental
    double *m_relDt = nullptr;

    //! number of time groups
    unsigned short m_nTgs = 1;

    //! number of elements per time group
    std::size_t * m_nTgEls = nullptr;

    //! element layout
    t_enLayout m_elLay;

    /**
     * Sets the data layout of the elements.
     *
     * @param i_nTgs number of time groups.
     * @param i_nTgEls number of elements in the time groups.
     * @param o_elLay will be set to the element layout.
     **/
    static void setElLayout( unsigned short         i_nTgs,
                             std::size_t    const * m_nTgEls,
                             t_enLayout           & o_elLay );

  public:
    /**
     * Sets LTS types (bitmasks) of the elements.
     *
     * @param i_nEls number of elements.
     * @param i_nElFas number of faces per element.
     * @param i_elFaEl face-adjacent elements.
     * @param i_nTgs number of time groups.
     * @param i_nTgsEl number of elements in the time groups.
     * @param i_elEq type if any of the face-adjacent elements has an equal relationship.
     * @param i_elLt type if any of the face-adjacent elements has an less-than relationship.
     * @param i_elGt type if any of the face-adjacent elements has an greather-than relation ship.
     * @param i_adEq types (one per face) if the adjacent element is in equal relationship (tsEl == tsEl).
     * @param i_adLt types (one per face) if the adjacent element is in less-than relationship (tsEl < tsAd).
     * @param i_adGt types (one per face) if the adjacent element is in greater-than relation ship (tsEl > tsAd).
     * @param o_spTys element sparse types which will be adjusted.
     **/
    static void setLtsTypes( std::size_t            i_nEls,
                             unsigned short         i_nElFas,
                             std::size_t    const * i_elFaEl,
                             unsigned short         i_nTgs,
                             std::size_t    const * i_nTgEls,
                             long long              i_elEq,
                             long long              i_elLt,
                             long long              i_elGt,
                             long long      const * i_adEq,
                             long long      const * i_adLt,
                             long long      const * i_adGt,
                             long long            * o_spTys );
    /**
     * Constructor.
     *
     * @param i_pathToMesh path to the mesh file.
     * @parma i_periodic if true, periodic boundaries as enforced.
     **/
    EdgeV( std::string const & i_pathToMesh,
           bool                i_periodic );

    /**
     * Destructor.
     **/
    ~EdgeV();

    /**
     * Gets the relative time steps.
     *
     * @return relative time steps.
     **/
    double const * getRelDt() const { return m_relDt; }

    /**
     * Gets the number of time groups.
     *
     * @return number of time groups.
     **/
    unsigned short nTgs() const { return m_nTgs; }

    /**
     * Gets the number of elements for the time groups.
     *
     * @return number of elements for the time groups.
     **/
    std::size_t const * nTgEls() const { return m_nTgEls; }

    /**
     * Gets the element layout.
     *
     * @return element layout.
     **/
    t_enLayout getElLayout() const { return m_elLay; }

    /**
     * Sets LTS type in the sparse types of the elements.
     *
     * @param io_elChars element characteristics whose sparse type is adjusted
     **/
    void setLtsTypes( t_elementChars * io_elChars ) const;
};

#endif