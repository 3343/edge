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

    //! number of inner elements per time group
    std::size_t * m_nTgElsIn = nullptr;

    //! number of send elements per time group
    std::size_t * m_nTgElsSe = nullptr;

    //! element layout
    t_enLayout m_elLay;

    //! communication structure
    std::size_t * m_commStruct = nullptr;

    //! number of communicating element-face pairs
    std::size_t m_nCommElFa = 0;

    //! send faces
    unsigned short * m_sendFa = nullptr;

    //! send elements
    std::size_t    * m_sendEl = nullptr;

    //! recv faces
    unsigned short * m_recvFa = nullptr;

    //! recv elements
    std::size_t    * m_recvEl = nullptr;

    //! vertex ids of the send elements' faces pointing to ghost elements
    unsigned short * m_sendVeIdsAd = nullptr;

    //! face ids of the send elements' faces pointing to ghost elements
    unsigned short * m_sendFaIdsAd = nullptr;

    /**
     * Sets the data layout of the elements.
     *
     * @param i_nTgs number of time groups.
     * @param i_nTgElsIn number of inner elements in the time groups.
     * @param i_nTgElsSe number of send elements in the time groups.
     * @param o_elLay will be set to the element layout.
     **/
    static void setElLayout( unsigned short         i_nTgs,
                             std::size_t    const * i_nTgElsIn,
                             std::size_t    const * i_nTgElsSe,
                             t_enLayout           & o_elLay );

  public:
    /**
     * Sets LTS types (bitmasks) of the elements.
     *
     * @param i_nEls number of elements.
     * @param i_nElFas number of faces per element.
     * @param i_elFaEl face-adjacent elements.
     * @param i_nTgs number of time groups.
     * @param i_nTgsElIn number of inner elements in the time groups.
     * @param i_nTgsElSe number of send elements in the time groups.
     * @param i_sendFa local face ids of the send element-face pairs.
     * @param i_sendEl element ids of the send element-face pairs.
     * @param i_commStruct communication structure.
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
                             std::size_t    const * i_nTgElsIn,
                             std::size_t    const * i_nTgElsSe,
                             unsigned short const * i_sendFa,
                             std::size_t    const * i_sendEl,
                             std::size_t    const * i_commStruct,
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
     * Gets the number of elements for the time groups for the inner part.
     *
     * @return number of elements for the time groups.
     **/
    std::size_t const * nTgElsIn() const { return m_nTgElsIn; }

    /**
     * Gets the number of elements for the time groups for the send part.
     *
     * @return number of elements for the time groups.
     **/
    std::size_t const * nTgElsSe() const { return m_nTgElsSe; }

    /**
     * Gets the number of vertices.
     *
     * @return number of vertices.
     **/
    std::size_t nVes() const { return m_mesh.nVes(); }

    /**
     * Gets the number of faces.
     *
     * @return number of faces.
     **/
    std::size_t nFas() const { return m_mesh.nFas(); }

    /**
     * Gets the number of elements.
     *
     * @return number of elements.
     **/
    std::size_t nEls() const { return m_mesh.nEls(); }

    /**
     * Sets LTS type in the sparse types of the elements.
     *
     * @param io_elChars element characteristics whose sparse type is adjusted.
     **/
    void setLtsTypes( t_elementChars * io_elChars ) const;

    /**
     * Sets the sparse types based on mesh-definitions.
     *
     * @param io_veChars vertex characteristics whose sparse type is adjusted.
     * @param io_faChars face characteristics whose sparse type is adjusted.
     * @param io_elChars elemenet characteristics whose sparse type is adjusted.
     **/
    void setSpTypes( t_vertexChars  * io_veChars,
                     t_faceChars    * io_faChars,
                     t_elementChars * io_elChars ) const;

    /**
     * Sets the adjacent vertex- and face-ids for respective ghost elements.
     *
     * @param io_veIdsAd will be updated with the vertex ids of the adjacent ghost elements.
     * @param io_faIdsAd will be updated with the face ids of the adjacent ghost elements.
     **/
    void setSeVeFaIdsAd( unsigned short * io_veIdsAd,
                         unsigned short * io_faIdsAd ) const;

    /**
     * Gets the element layout.
     *
     * @return element layout.
     **/
    t_enLayout getElLayout() const { return m_elLay; }

    /**
     * Gets the vertices of the faces.
     *
     * @return faVe info.
     **/
    std::size_t const * getFaVe() const { return m_mesh.getFaVe(); }

    /**
     * Gets the elements adjacent to the faces.
     *
     * @return faEl info.
     **/
    std::size_t const * getFaEl() const { return m_mesh.getFaEl(); }

    /**
     * Gets the vertices adjacent to the elements.
     *
     * @return elVe info.
     **/
    std::size_t const * getElVe() const { return m_mesh.getElVe(); }

    /**
     * Gets the faces adjacent to the elements.
     *
     * @return elFa info.
     **/
    std::size_t const * getElFa() const { return m_mesh.getElFa(); }

    /**
     * Gets the face-adjacent elements.
     *
     * @return elFaEl info.
     **/
    std::size_t const * getElFaEl() const { return m_mesh.getElFaEl(); }

    /**
     * Gets the vertex coordinates.
     *
     * @return vertex coordinates.
     **/
    double const (* getVeCrds() )[3]{ return m_mesh.getVeCrds(); }

    /**
     * Gets the normals of the faces.
     *
     * @return normals of the faces.
     **/
    double const (* getNormalsFa() )[3] { return m_mesh.getNormalsFa(); }

    /**
     * Gets the tangents of the faces.
     *
     * @return tangents of the faces.
     **/
    double const (* getTangentsFa() )[2][3] { return m_mesh.getTangentsFa(); }

    /**
     * Gets the areas of the faces.
     *
     * @return areas.
     **/
    double const * getAreasFa() { return m_mesh.getAreasFa(); }

    /**
     * Gets the volumes of the elements.
     *
     * @return volumes.
     **/
    double const * getVolumesEl() { return m_mesh.getVolumesEl(); }

    /**
     * Gets the in-diameters of the elements.
     *
     * @return in-diameters.
     **/
    double const * getInDiasEl() { return m_mesh.getInDiasEl(); }

    /**
     * Gets the number of communicating element-face pairs.
     *
     * @return number of communicating element-face pairs.
     **/
    std::size_t nCommElFa() { return m_nCommElFa; }

    /**
     * Gets the communication structure.
     *
     * @return comm structure
     **/
    std::size_t const * getCommStruct() { return m_commStruct; }

    /**
     * Gets the local, per-element ids of the sending faces.
     *
     * @return face ids.
     **/
    unsigned short const * getSendFa() { return m_sendFa; }

    /**
     * Gets the ids of the sending elements.
     *
     * @return element ids.
     **/
    std::size_t const * getSendEl() { return m_sendEl; }

    /**
     * Gets the local, per-element ids of the receiving faces.
     *
     * @return face ids.
     **/
    unsigned short const * getRecvFa() { return m_recvFa; }

    /**
     * Gets the ids of the receiving elements.
     *
     * @return element ids.
     **/
    std::size_t const * getRecvEl() { return m_recvEl; }
};

#endif