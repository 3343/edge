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
 * Mesh functions and storage.
 **/
#ifndef EDGE_V_MESH_MESH_H
#define EDGE_V_MESH_MESH_H

#include "../constants.h"

#include "../io/Moab.h"

namespace edge_v {
  namespace mesh {
    class Mesh;
  }
}

/**
 * Mesh-related functions and data.
 **/
class edge_v::mesh::Mesh {
  private:
    //! element type
    t_entityType m_elTy;

    //! number of vertices in the mesh
    std::size_t m_nVes;

    //! number of faces in the mesh
    std::size_t m_nFas;

    //! number of elements in the mesh
    std::size_t m_nEls;

    //! vertex coordinates
    double (*m_veCrds)[3];

    //! vertices adjacent to the faces
    std::size_t *m_faVe;

    //! elements adjacent to the faces
    std::size_t *m_faEl;

    //! faces adjacent to the elements
    std::size_t *m_elFa;

    //! vertices adjacent to the elements
    std::size_t *m_elVe;

    //! elements adjacent to elements, faces as bridge
    std::size_t *m_elFaEl;

    //! indiameter of the elements
    double *m_inDiasEl;

    //! area/length of the faces
    double * m_volFa = nullptr;

    //! volumes of the elements
    double * m_volEl = nullptr;

    //! normals of the faces
    double (* m_normals)[3] = nullptr;

    //! tangents of the faces
    double (* m_tangents)[2][3] = nullptr;

    /**
     * Returns an entry in the second array, which is not in the first one.
     *
     * @param i_sizeFirst size of the first array.
     * @param i_sizeSecond size of the second array.
     * @param i_first first array.
     * @param i_second second array.
     * @return first found entry present in the second but not in the first one. numeric_limits< std::size_t >::max() if none.
     **/
    static std::size_t getAddEntry( std::size_t   i_sizeFirst,
                                    std::size_t   i_sizeSecond,
                                    std::size_t * i_first,
                                    std::size_t * i_second );

    /**
     * Gathers the vertex coordinates of the given entity.
     *
     * @param i_enTy entity type.
     * @param i_enVe vertices adjacent to this entity.
     * @param i_veCrd vertex coordinates in the mesh.
     * @param o_enVeCrds will be set to the coordinates of the entity's vertices.
     **/
    static void getEnVeCrds( t_entityType          i_enTy,
                             std::size_t  const  * i_enVe,
                             double       const (* i_veCrds)[3],
                             double             (* o_enVeCrds)[3] );

    /**
     * Computes the lengths (1d), incircle (2d) or insphere (3d) diameters.
     *
     * @param i_enTy entity type.
     * @param i_nEns number of entities.
     * @param i_enVe vertices adjacent to the entities.
     * @param i_veCrds coordinates of the vertices.
     * @param o_inDia will be set to the entities' diameters.
     **/
    static void setInDiameter( t_entityType         i_enTy,
                               std::size_t          i_nEns,
                               std::size_t const  * i_enVe,
                               double      const (* i_veCrds)[3],
                               double             * o_inDia );

    /**
     * Inserts periodic boundaries into faEl and elFaEl by comparing the vertex coordinates of the faces.
     * The function assumes that the periodic boundaries are coordinate aligned.
     *
     * @param i_elTy element type.
     * @param i_nFas number of faces.
     * @param i_faVe vertices of the faces.
     * @param i_veCrds vertex coordinates.
     * @param io_faEl elements adjacent to the faces.
     * @param i_elFa faces adjacent to the elements.
     * @param io_elFaEl elements adjacent to elements (faces as bridge).
     * @param o_pFasGt will be set to periodic faces originally only adjacent to the element with the greater id.
     **/
    static void setPeriodicBnds( t_entityType                        i_elTy,
                                 std::size_t                         i_nFas,
                                 std::size_t                const  * i_faVe,
                                 double                     const (* i_veCrds)[3],
                                 std::size_t                       * io_faEl,
                                 std::size_t                const  * i_elFa,
                                 std::size_t                       * io_elFaEl,
                                 std::vector< std::size_t >        & o_pFasGt );

    /**
     * Normalizes the order of the given adjacency information.
     *   1) faVe: vertices with smaller id first for each face.
     *   2) faEl: elements with smaller id first for each face.
     *   3) elVe: smallest vertex first, then ensuring counterclockwise ordering.
     *   4) elFa: following the element's vertex ordering w.r.t. face-assignments in the reference element.
     *   5) elFaEl: same as elFa.
     *
     * @param i_elTy element type.
     * @param i_nFas number of faces.
     * @param i_nEls number of elements.
     * @param i_veCrds coordinates of the element's vertices.
     * @param i_faVe vertices of the faces.
     * @param i_faEl elements adjacent to the faces.
     * @param io_elVe vertices adjacent to the element (ordered ascending by the ids).
     * @param io_elFa faces adjacent to the element (ordered ascending by the vertex ids of the faces).
     * @param io_elFaEl elements adjacent to the elements (ordered ascending by the vertex ids of the faces).
     **/
    static void normOrder( t_entityType         i_elTy,
                           std::size_t          i_nFas,
                           std::size_t          i_nEls,
                           double      const (* i_veCrds)[3],
                           std::size_t        * io_faVe,
                           std::size_t        * io_faEl,
                           std::size_t        * io_elVe,
                           std::size_t        * io_elFa,
                           std::size_t        * io_elFaEl );

  public:
    /**
     * Constructor.
     *
     * @param i_moab moab interface.
     * @param i_periodic if true search insert periodic adjacency info for faces at the boundary.
     **/
    Mesh( io::Moab const & i_moab,
          bool             i_periodic = false );

    /**
     * Destructor
     **/
    ~Mesh();

    /**
     * Prints the mesh statistics.
     **/
    void printStats() const;

    /**
     * Gets the vertex type of the mesh.
     *
     * @return vertex type.
     **/
    t_entityType getTypeVe() const{ return edge_v::POINT; }

    /**
     * Gets the face type of the mesh.
     *
     * @return face type.
     **/
    t_entityType getTypeFa() const{ return CE_T_FA(m_elTy); }

    /**
     * Gets the element type of the mesh.
     *
     * @return element type.
     **/
    t_entityType getTypeEl() const { return m_elTy; }

    /**
     * Gets the number of vertices.
     *
     * @return number of vertices.
     **/
    std::size_t nVes() const { return m_nVes; }

    /**
     * Gets the number of faces.
     *
     * @return number of faces.
     **/
    std::size_t nFas() const { return m_nFas; }

    /**
     * Gets the number of elements.
     *
     * @return number of elements in the mesh.
     **/
    std::size_t nEls() const { return m_nEls; }

    /**
     * Gets the vertices adjacent to the elements.
     *
     * @return connectivity info.
     **/
    std::size_t const * getElVe() const { return m_elVe; }

    /**
     * Gets the vertices adjacent to the faces.
     *
     * @return faVe info.
     **/
    std::size_t const * getFaVe() const { return m_faVe; }

    /**
     * Gets the elements adjacent to the faces.
     *
     * @return faEl info.
     **/
    std::size_t const * getFaEl() const { return m_faEl; }

    /**
     * Gets the faces adjacent to the elements.
     *
     * @return elFa info.
     **/
    std::size_t const * getElFa() const { return m_elFa; }

    /**
     * Gets the elements adjacent to the elements (faces as bridge).
     *
     * @return elFaEl info.
     **/
    std::size_t const * getElFaEl() const { return m_elFaEl; }

    /**
     * Gets the vertex coordinates.
     *
     * @return vertex coordinates.
     **/
    double const (* getVeCrds() const)[3] { return m_veCrds; }

    /**
     * Gets the area/length of the faces.
     *
     * @return area/length.
     **/
    double const * getAreasFa();

    /**
     * Gets the volumes of the elements.
     *
     * @return volumes.
     **/
    double const * getVolumesEl();

    /**
     * Gets the normals of the faces.
     *
     * @return normals of the faces.
     **/
    double const (* getNormalsFa() )[3];

    /**
     * Gets the tangents of the faces.
     *
     * @return tangents of the faces.
     **/
    double const (* getTangentsFa() )[2][3];

    /**
     * Gets the length (1d), incircle (2d) or insphere (3d) diameters of the elements.
     *
     * @return element diameter.
     **/
    double const * getInDiasEl() const { return m_inDiasEl; }
};

#endif