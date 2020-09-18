/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2020, Friedrich Schiller University Jena
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
 * Mesh functions and storage.
 **/
#ifndef EDGE_V_MESH_MESH_H
#define EDGE_V_MESH_MESH_H

#include "../constants.h"
#include "../io/Gmsh.h"
#include <limits>

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
    t_idx m_nVes;

    //! number of faces in the mesh
    t_idx m_nFas;

    //! number of elements in the mesh
    t_idx m_nEls;

    //! vertex coordinates
    double (*m_veCrds)[3] = nullptr;

    //! vertices adjacent to the faces
    t_idx * m_faVe = nullptr;

    //! elements adjacent to the faces
    t_idx * m_faEl = nullptr;

    //! faces adjacent to the elements
    t_idx * m_elFa = nullptr;

    //! vertices adjacent to the elements
    t_idx * m_elVe = nullptr;

    //! elements adjacent to elements, faces as bridge
    t_idx * m_elFaEl = nullptr;

    //! sparse type of the vertices
    t_sparseType * m_spTypeVe = nullptr;

    //! sparse type of the faces
    t_sparseType * m_spTypeFa = nullptr;

    //! sparse type of the elements
    t_sparseType * m_spTypeEl = nullptr;

    //! indiameter of the elements
    double * m_inDiasEl = nullptr;

    //! area/length of the faces
    double * m_volFa = nullptr;

    //! volumes of the elements
    double * m_volEl = nullptr;

    //! normals of the faces
    double (* m_normals)[3] = nullptr;

    //! tangents of the faces
    double (* m_tangents)[2][3] = nullptr;

   /**
    * Sorts the given data lexicographically.
    * Total number of entries: n0 x n1.
    *
    * @param i_n0 number of entries in the slow dimension.
    * @param i_n1 number of entries in the fast dimension.
    * @param io_data data which is sorted.
    * @param o_mapping will be set to mapping from sorted to original. nullptr if unused.
    **/
   static void sortLex( t_idx   i_n0,
                        t_idx   i_n1,
                        t_idx * io_data,
                        t_idx * o_mapping = nullptr );

    /**
     * Gets the element-to-element adjacency (faces as bridge).
     *
     * @param i_elTy element type.
     * @param i_nEls number of elements.
     * @param i_faEl elements adjacent to the faces.
     * @param i_elFa faces adjacent to the elements.
     * @param o_elFaEl will be set to element-to-element adjacency.
     **/
    static void getElFaEl( t_entityType         i_elTy,
                           t_idx                i_nEls,
                           t_idx        const * i_faEl,
                           t_idx        const * i_elFa,
                           t_idx              * o_elFaEl );

    /**
     * Returns an entry in the second array, which is not in the first one.
     *
     * @param i_sizeFirst size of the first array.
     * @param i_sizeSecond size of the second array.
     * @param i_first first array.
     * @param i_second second array.
     * @return first found entry present in the second but not in the first one. numeric_limits< t_idx >::max() if none.
     **/
    static t_idx getAddEntry( t_idx         i_sizeFirst,
                              t_idx         i_sizeSecond,
                              t_idx const * i_first,
                              t_idx const * i_second );

    /**
     * Gathers the vertex coordinates of the given entity.
     *
     * @param i_enTy entity type.
     * @param i_enVe vertices adjacent to this entity.
     * @param i_veCrd vertex coordinates in the mesh.
     * @param o_enVeCrds will be set to the coordinates of the entity's vertices.
     **/
    static void getEnVeCrds( t_entityType          i_enTy,
                             t_idx        const  * i_enVe,
                             double       const (* i_veCrds)[3],
                             double             (* o_enVeCrds)[3] );

    /**
     * Adds the given type to the sparse entities.
     *
     * Remark:
     *   All per-entity vertices are assumed to be sorted ascending.
     *   All entity arrays are assumed to be sorted lexicographically.
     *
     * @param i_enTy entity type.
     * @param i_nEnsDe number of dense entities.
     * @param i_nEnsSp number of sparse entities.
     * @param i_enVeDe vertex ids of the dense entities.
     * @param i_enVeSp vertex ids of the sparse entities.
     * @param i_spTypeAdd sparse type which will be added to the entities' sparse types.
     * @param io_spType will be updated by spTypeAdd if defined for the entity.
     **/
    static void addSparseTypeEn( t_entityType         i_enTy,
                                 t_idx                i_nEnsDe,
                                 t_idx                i_nEnsSp,
                                 t_idx        const * i_enVeDe,
                                 t_idx        const * i_enVeSp,
                                 t_sparseType         i_spTypeAdd,
                                 t_sparseType       * io_spType );

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
                               t_idx                i_nEns,
                               t_idx       const  * i_enVe,
                               double      const (* i_veCrds)[3],
                               double             * o_inDia );

    /**
     * Inserts periodic boundaries into faEl and elFaEl by comparing the vertex coordinates of the faces.
     * The function assumes that the periodic boundaries are coordinate aligned.
     *
     * @param i_elTy element type.
     * @param i_nFas number of faces.
     * @param i_peBndTy periodic boundary type.
     * @param i_faBndTys boundary types of the faces.
     * @param i_faVe vertices of the faces.
     * @param i_veCrds vertex coordinates.
     * @param io_faEl elements adjacent to the faces.
     * @param i_elFa faces adjacent to the elements.
     * @param io_elFaEl elements adjacent to elements (faces as bridge).
     * @param o_pFasGt will be set to periodic faces originally only adjacent to the element with the greater id.
     **/
    static void setPeriodicBnds( t_entityType                  i_elTy,
                                 t_idx                         i_nFas,
                                 int                           i_peBndTy,
                                 int                  const  * i_faBndTys,
                                 t_idx                const  * i_faVe,
                                 double               const (* i_veCrds)[3],
                                 t_idx                       * io_faEl,
                                 t_idx                const  * i_elFa,
                                 t_idx                       * io_elFaEl,
                                 std::vector< t_idx >        & o_pFasGt );

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
    static void normOrder( t_entityType          i_elTy,
                           t_idx                 i_nFas,
                           t_idx                 i_nEls,
                           double       const (* i_veCrds)[3],
                           t_idx               * io_faVe,
                           t_idx               * io_faEl,
                           t_idx               * io_elVe,
                           t_idx               * io_elFa,
                           t_idx               * io_elFaEl );

  public:
    /**
     * Constructor.
     *
     * @param i_gmsh gmsh interface.
     * @param i_periodic if not max insert periodic adjacency info for faces at the boundary if available.
     **/
    Mesh( io::Gmsh const & i_gmsh,
          int              i_periodic = std::numeric_limits< int >::max() );

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
     * Gets the number of vertices for an element.
     *
     * @return number of vertices for an element.
     **/
    unsigned short nElVes() const { return CE_N_VES( m_elTy ); }

    /**
     * Gets the number of vertices.
     *
     * @return number of vertices.
     **/
    t_idx nVes() const { return m_nVes; }

    /**
     * Gets the number of faces.
     *
     * @return number of faces.
     **/
    t_idx nFas() const { return m_nFas; }

    /**
     * Gets the number of elements.
     *
     * @return number of elements in the mesh.
     **/
    t_idx nEls() const { return m_nEls; }

    /**
     * Gets the vertices adjacent to the elements.
     *
     * @return connectivity info.
     **/
    t_idx const * getElVe() const { return m_elVe; }

    /**
     * Gets the vertices adjacent to the faces.
     *
     * @return faVe info.
     **/
    t_idx const * getFaVe() const { return m_faVe; }

    /**
     * Gets the elements adjacent to the faces.
     *
     * @return faEl info.
     **/
    t_idx const * getFaEl() const { return m_faEl; }

    /**
     * Gets the faces adjacent to the elements.
     *
     * @return elFa info.
     **/
    t_idx const * getElFa() const { return m_elFa; }

    /**
     * Gets the elements adjacent to the elements (faces as bridge).
     *
     * @return elFaEl info.
     **/
    t_idx const * getElFaEl() const { return m_elFaEl; }

    /**
     * Gets the sparse types of the vertices.
     *
     * @return vertex sparse types.
     **/
    t_sparseType const * getSpTypeVe() const { return m_spTypeVe; }

    /**
     * Gets the sparse types of the faces.
     *
     * @return face sparse types.
     **/
    t_sparseType const * getSpTypeFa() const { return m_spTypeFa; }

    /**
     * Gets the sparse types of the elements.
     *
     * @return element sparse types.
     **/
    t_sparseType const * getSpTypeEl() const { return m_spTypeEl; }

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

    /**
     * Gets the adjacent elements' face ids for the given element-faces.
     *
     * @param i_nFas number of element-faces to get ids for.
     * @param i_elOff element offset (added to element-ids in queries).
     * @param i_el elements to which the faces belong.
     * @param i_fa local face ids w.r.t. the elements.
     * @param o_faIdsAd will be set to faces ids w.r.t. to the adjacent elements.
     **/
    void getFaIdsAd( t_idx                  i_nFas,
                     t_idx                  i_elOff,
                     t_idx          const * i_el,
                     unsigned short const * i_fa,
                     unsigned short       * o_faIdsAd ) const;

    /**
     * Gets the adjacent elements' vertex ids for the given element-faces.
     * This is the local id of the vertex, which matches the faces' first vertex.
     *
     * @param i_nFas number of element-faces to get ids for.
     * @param i_elOff element offset (added to element-ids in queries).
     * @param i_el elements to which the faces belong.
     * @param i_fa local face ids w.r.t. the elements.
     * @param o_veIdsAd will be set to vertex ids w.r.t. to the adjacent elements.
     **/
    void getVeIdsAd( t_idx                  i_nFas,
                     t_idx                  i_elOff,
                     t_idx          const * i_el,
                     unsigned short const * i_fa,
                     unsigned short       * o_veIdsAd ) const;
};

#endif