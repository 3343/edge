/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2015-2018, Regents of the University of California
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
 * Regular meshes in EDGE.
 **/

#ifndef EDGE_MESH_REGULAR_H_
#define EDGE_MESH_REGULAR_H_

#include "constants.hpp"
#include "data/EntityLayout.type"

namespace edge {
  namespace mesh {
    class Regular;
  }
}

/**
 * Regular meshes based on a hexahedral base.
 * Other element types (tetrahedrons, etc.) are derived by splitting the hexahedrons.
 * 2D elements are obtained by using collapsing the z-dimension, 1D elements by the y- and z-dimension.
 * Meshes are created on the fly, no file I/O is required.
 * TODO: periodic boundaries are assumed for the time being.
 **/
class edge::mesh::Regular {
  public:
    // Requested element type
    enum Type {
      Line,
      Quadrilateral,
      Hexahedral
    };

  private:
    //! type of the elements
    Type m_elementType;

    //! number of base hexahedrons in x-, y- and z-direction
    int_el m_nX, m_nY, m_nZ;

    //! size of the local domain in x-, y- and z-direction
    double m_sizeX, m_sizeY, m_sizeZ;

    //! base width in x-, y- and z_direction
    double m_dX, m_dY, m_dZ;

    //! number of base elements
    int_el m_nBaseElements;

    //! number of base faces
    int_el m_nBaseFaces;

    //! number of base vertices
    int_el m_nBaseVertices;

    //! number of requested elements
    int_el m_nRequestedElements;

    //! number of requested faces
    int_el m_nRequestedFaces;

    //! number of requested vertices
    int_el m_nRequestedVertices;

    //! mapping of indices from mesh to data and data to mesh
    t_inMap m_inMap;

    //! global ids of the vertices
    std::vector< int_gid > m_gIdsVe;

    //! global ids of the faces
    std::vector< int_gid > m_gIdsFa;

    //! global ids of the elements
    std::vector< int_gid > m_gIdsEl;

    /**
     * Gets the ids of the elements adjacent to the face.
     *
     * @param i_faceId id of the faces.
     * @param o_neighboringIds will be set to ids of the adjacent elements.
     **/
    void getFaceAdjacentElements( int_el i_faceId,
                                  int_el o_neighboringIds[2] ) const;

    /**
     * Gets the ids of the face neighboring elements.
     *
     * @param i_elementId id of the element.
     * @param o_neighboringIds will be set to ids of the face neighbors.
     **/
    void getElementFaceNeighbors( int_el  i_elementId,
                                  int_el *o_neighboringIds ) const;


    /**
     * Gets the ids of vertex neighboring elements.
     * The first entry of the array of pointers is expected to hold the raw data for all info.
     * Given nEl elements, an additional pointer (nEl+1) will be set for consistent size computations.
     *
     * @param o_elVeEl will be set to pointers pointing to the data of the respective element.
     **/
    void getElVeEl( int_el** o_elVeEl ) const;

    /**
     * Initializes the given array with the face neighbor ids for all elements.
     *
     * @param o_neighboringIds will be set to ids of all elements' face neighbors.
     **/
    void getElementsFaceNeighbors( int_el (*o_neighboringIds)[C_ENT[T_SDISC.ELEMENT].N_FACES] ) const;

    /**
     * Initializes the given array with ids of adjacent faces for all elements.
     *
     * @param o_elementAdjacentFaces will be set to ids of all elements' adjacent faces.
     **/
    void getElementsAdjacentFaces( int_el (*o_elementAdjacentFaces)[C_ENT[T_SDISC.ELEMENT].N_FACES] ) const;

    /**
     * Gets the adjacent vertices for the given faces.
     *
     * @param i_fa considered face.
     * @param o_faVe will be set to adjacent vertices of the faces.
     **/
    void getFaceAdjacentVertices( int_el i_fa,
                                  int_el o_faVe[C_ENT[T_SDISC.ELEMENT].N_FACE_VERTICES] ) const;

    /**
     * Initializes the given array with the adjacent vertices for all faces.
     *
     * @param o_faVe will be set to adjacent vertices of the faces.
     **/
    void getFacesAdjacentVertices( int_el (*o_faVe)[C_ENT[T_SDISC.ELEMENT].N_FACE_VERTICES] ) const;

    /**
     * Initializes the given array with the adjacent elements for all faces.
     *
     * @param o_neighboringIds will be set to neighboring element ids for all faces.
     **/
    void getFacesAdjacentElements( int_el (*o_neighboringIds)[2] ) const;

    /**
     * Gets the vertices adjacent to a single element.
     *
     * @param o_elementAdjacentVertices will be set to vertices adjacent to the element.
     * @param i_px element's position in x-direction.
     * @param i_py element's position in y-direction.
     * @param i_pz element's position in z-direction.
     **/
    void getElementAdjacentVertices( int_el o_elementAdjacentVertices[C_ENT[T_SDISC.ELEMENT].N_VERTICES],
                                     int_el i_px,
                                     int_el i_py = 0,
                                     int_el i_pz = 0 ) const;

    /**
     * Gets the vertices adjacent to the elements.
     * Remark: Vertces on periodic boundaries are included "redundantly".
     *
     * @param o_elementAdjacentVertices will be set to vertices adjacent to the elements.
     **/
    void getElementsAdjacentVertices( int_el (*o_elementAdjacentVertices)[C_ENT[T_SDISC.ELEMENT].N_VERTICES] ) const;

    /**
     * Computes the volume of an element.
     *
     * @param i_element element which volume is computed.
     * @param o_volume volume of the element.
     **/
    void computeElementVolume( int_el     i_element,
                               real_mesh &o_volume ) const;

    /**
     * Compute the insphere diameter of an element.
     *
     * @param i_element element which diameter is computed.
     * @param o_inDia will be set to insphere diameter of the element.
     **/
    void computeElementInDia( int_el     i_element,
                              real_mesh &o_dia ) const;

    /**
     * Computes the area of a face.
     *
     * @param i_face face which area is computed.
     * @param o_area will be set to area of the face.
     **/
    void computeFaceArea( int_el     i_face,
                          real_mesh &o_area ) const;

    /**
     * Computes the outer pointing normal of a face.
     * Remark: Our outer pointing normals always point from face-neighbor 0 to face-neighbor 1 in o_neighboringIds of get FacesAdjacentElements.
     *
     * @param i_face face which outer pointing normal is computed.
     * @param o_x will be set to x-coordinate of outer pointing normal.
     * @param o_y will be set to y-coordinate of outer pointing normal.
     * @param o_z will be set to z-coordinate of outer pointing normal.
     **/
    void computeOutPointNormal( int_el     i_face,
                                real_mesh &o_x,
                                real_mesh &o_y,
                                real_mesh &o_z ) const;

    /**
     * Computes the tangents for a face.
     *
     * @param i_face face which tangents are computed.
     * @param o_tangent0 will be set to first tangent.
     * @param o_tangent1 will be set to second tangent.
     **/
    void computeTangents( int_el    i_face,
                          real_mesh o_tangent0[3],
                          real_mesh o_tangent1[3] ) const;

  public:
    /**
     * Sets up a new regular mesh.
     *
     * @param i_elementType type of the requested elements.
     * @param i_nX number of base elements in x-direction.
     * @param i_sizeX size of local domain in x-direction.
     * @param i_nY number of base elements in y-direction.
     * @param i_sizeY size of local domain in y-direction
     * @param i_nZ number of base elements in z-direction.
     * @param i_sizeZ size of local domain in z-direction
     **/
    Regular( Type i_elementType,
             int_el i_nX,
             double i_sizeX,
             int_el i_ny    = 0,
             double i_sizeY = 0.0,
             int_el i_nz    = 0,
             double i_sizeZ = 0.0 );

    /**
     * Gets the global ids of the vertices.
     * TODO: Fix implementation, this copies all data!
     **/
    void getGIdsVe( std::vector< int_gid > &o_gIds ) { o_gIds = m_gIdsVe; };

    /**
     * Gets the global ids of the faces.
     * TODO: Fix implementation, this copies all data!
     **/
    void getGIdsFa( std::vector< int_gid > &o_gIds ) { o_gIds = m_gIdsFa; };

    /**
     * Gets the global ids of the elements.
     * TODO: Fix implementation, this copies all data!
     **/
    void getGIdsEl( std::vector< int_gid > &o_gIds ) { o_gIds = m_gIdsEl; };

    /**
     * Gets the data layout of the elements (N_DIM).
     *
     * @return data layout of the elements.
     **/
    t_enLayout getElLayout() const;

    /**
     * Gets the data layout of the faces (N_DIM-1).
     *
     * @return data layout of the faces.
     **/
    t_enLayout getFaLayout() const;

    /**
     * Gets the data layout of the vertices (0).
     *
     * @return data layout of the vertices.
     **/
    t_enLayout getVeLayout() const;

    /**
     * Gets the index mapping for mesh-to-data and data-to-mesh.
     *
     * @return will be set to index mapping.
     **/
    const t_inMap* getInMap() const;

    /**
     * Gets the number of elements
     *
     * @return number of elements.
     **/
    int_el getNElements() const;

    /**
     * Gets the number of faces
     *
     * @return number of faces.
     **/
    int_el getNFaces() const;

    /**
     * Gets the number of vertices.
     *
     * @return number of vertices.
     **/
    int_el getNVertices() const;

    /**
     * Gets the total number of entries for elements adjacent through vertices.
     *
     * @return total number of entries.
     **/
    int_el getNelVeEl() const;

    /**
     * Gets the connectivity information.
     *  elVe:   element -> vertices
     *  elFa:   element -> faces
     *  faEl:   face    -> elements
     *  elFaEl: element -> faces -> elements (excluding starting element)
     *
     * Remark: Faces on periodic boundaries show up only once (in no particular order).
     *         Vertices on periodic boundaries show up reduantly matching the mesh.
     *
     * @param i_veChars vertex characteritics.
     * @param i_faChars face characteristics.
     * @param o_connect will be set to connectity information
     **/
    void getConnect( const t_vertexChars *i_veChars,
                     const t_faceChars   *faChars,
                           t_connect     &o_connect ) const;

    /**
     * Gets the elements' characteristics.
     *
     * @param o_elementChars will be written to the elements' characteristics.
     **/
    void getElChars( t_elementChars* o_elementChars ) const;

    /**
     * Gets the faces' characteristics.
     *
     * @param o_faceChars will be written to the faces' charactersitiscs.
     **/
    void getFaChars( t_faceChars* o_faceChars ) const;

    /**
     * Gets the vertices' characteristics.
     *
     * @param o_vertexChars will be written to vertex characteristics.
     **/
    void getVeChars( t_vertexChars* o_vertexChars ) const;
};

#endif
