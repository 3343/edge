/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
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
 * Interface to MOAB.
 **/

#ifndef EDGE_MESH_MOAB_H
#define EDGE_MESH_MOAB_H

#ifdef PP_USE_MPI
#include <moab/ParallelComm.hpp>
#endif
#include "constants.hpp"
#include "data/EntityLayout.type"
#include <string>
#include <cstdlib>
#include <moab/Core.hpp>
#include "io/logging.h"

namespace edge {
  namespace mesh {
    class Moab;
  }
}

class edge::mesh::Moab {
  //private:
    //! moab core
    moab::Core m_core;

#ifdef PP_USE_MPI
    //! parallel communications in moab
    moab::ParallelComm m_pcomm;
#endif

    //! dimension of the mesh
    const unsigned int m_dim;

    //! number of values for boundary regions
    const unsigned int m_nBndVals;

    //! boundary regions
    const int *m_bndVals;

    //! boundary sets of the mesh
    moab::Range m_bndSets;

    //! element layout
    t_enLayout m_elLayout;

    //! face layout
    t_enLayout m_faLayout;

    //! vertex layout
    t_enLayout m_veLayout;

    //! vertices (including duplicated and ghost vertices)
    std::vector< moab::EntityHandle > m_vertices;

    //! faces (including duplicated and ghost faces)
    std::vector< moab::EntityHandle > m_faces;

    //! elements (including duplicated ghost elements)
    std::vector< moab::EntityHandle > m_elements;

    //! mapping of indices from mesh to data and data to mesh
    t_inMap m_inMap;

    //! tags in the mesh
    std::vector< moab::Tag > m_tagsMesh;

    //! names of the tags
    std::vector< std::string > m_tagsNames;

    //! material tag, as defined in MATERIAL_SET of MOAB.
    // While we use this to define true material regions w.r.t to elements, this tag is also used to define boundary conditions on faces.
    moab::Tag m_tagMat;

    //! local id of ents w.r.t. the current partition. This id is equivalent to the sorting/storage in memory.
    moab::Tag m_tagLId;

    //! global if of ents
    moab::Tag m_tagGId;

    //! default value for the local ids
    const int_el m_defaultTagLId = std::numeric_limits<int_el>::max();

    //! default value for the global ids
    const int_el m_defaultTagGId = std::numeric_limits<int_gid>::max();

#ifndef PP_USE_MPI
    //! periodic tag value
    const int m_periodicVal;
#endif

    //! handles of the periodic faces, the matching handle of the element, not in MOAB's adjacencies is in periodic elements remote vector.
    std::vector< moab::EntityHandle > m_periodicFaces;

    //! handles of the remote periodic element, the matching handle of the face, not in MOAB's adjacencies is in the periodic face vector
    std::vector< moab::EntityHandle > m_periodicElementsRemote;

    /**
     * Gets the entity handle for the given entity.
     *
     * @param i_nDim number of dimensions of the entity.
     * @param i_en local id of the entity.
     **/
    moab::EntityHandle getEnHandle( unsigned short i_nDim,
                                    int_el         i_en );

    /**
     * Determines if the given entity is at the boundary.
     * Supported are faces only.
     * For faces (dim-1) this is true if the face is a boundary face.
     *
     * @param entity handle of the given entity.
     * @return true if at a boundary, false otherwise.
     **/
    bool isBnd( moab::EntityHandle i_ent );

    /**
     * Adds the periodic vertices (if existing), which are equivalent to the given vertex.
     *
     * @param i_ve vertex for which the periodic equivalents are determined.
     * @param io_ves will be updated with ids of the periodic vertices.
     **/
    void addPeriodicVertices( moab::EntityHandle  i_ve,
                              moab::Range        &io_ves );

    /**
     * Adds the periodic element (if existing), which is adjacent to the given face.
     *
     * @param i_face entity handle of the face.
     * @param io_elements entity handle of the range to which the element will be added to.
     **/
    void addPeriodicElement( moab::EntityHandle  i_face,
                             moab::Range        &io_elements );

    /**
     * Adds the periodic element (if existing), which is adjacent to the given face.
     *
     * @param i_face entity handle of the face.
     * @param io_elements vector to which the element will be added to (at the end).
     **/
    void addPeriodicElement( moab::EntityHandle                 i_face,
                             std::vector< moab::EntityHandle > &io_elements );

    /**
     * Sets the periodic adjacencies for a single dimension.
     *
     * @param i_faces ids of faces in the current dimension. These are equivalent to the array indices in i_handles and i_freeDims.
     * @param i_handles handles of all periodic boundary faces, not only the current dimension.
     * @param i_freeDims free dimensions of the boundary face.
     **/
    void setPeriodicAdjacenciesDim( const std::vector< int_el >  &i_faces,
                                    const moab::Range            &i_handles,
                                    const real_mesh             (*i_freeDims)[2] );

    /**
     * Sets the adjacencies for periodic boundaries and removes one of the redundant faces from local storage.
     * Remark: This only works for periodic boundaries aligned to the mesh dimensions.
     *
     * @param i_tagValue is the value defined for the faces in the mesh, which belong to periodic boundaries.
     **/
    void setPeriodicAdjacencies( int i_tagValue );

    /**
     * Fixes the id tags for faces, which have been deleted from m_faces but might be queried from the core.
     * Remark: This function assumes that every second face in the periodic face list was removed from m_faces
     *         in the call of setPeriodicAdjacencies.
     **/
    void fixPeriodicTagIds(); 

#ifdef PP_USE_MPI
    /**
     * Gets the MPI-type of an entity.
     *
     * @param i_en entity handle of the entity in question.
     * @param i_rank if -1 all MPI-neighbors are considered for the neighboring ranks w.r.t. to send- and recv-entities.
     *               otherwise the non-local rank has to match the given rank.
     * @return 0 if inner-entity,
     *         1 if send-entity,
     *         2 if recv-entity,
     *         std::numeric_limits< unsigned short >::max() if none of the three.
     **/
    unsigned short getEnMpiType( moab::EntityHandle i_en,
                                 int                i_rank=-1 );
#endif

    /**
     * Sets up the data layout with respect to a single entity and time group.
     *
     * Remark: It's assumed that nents of the enLayout is set appropiately.
     *
     * @param i_tg time group for which the data layout is set.
     * @param i_ents input entities (including non-owned) as given by MOAB.
     * @param o_enLayout will be set to data layout of the entity.
     * @param o_ents will be set to a vector of entitiy handle following the data layout.
     **/
    void setupEnLayout( int_tg                             i_tg,
                        moab::Range                       &i_ents,
                        t_enLayout                        &o_enLayout,
                        std::vector< moab::EntityHandle > &o_ents );

    /**
     * Sets up the data layout (m_elements, m_faces, m_vertices).
     *
     * Our layout consists of a strict ordering of elements, faces and vertices according to their time step
     * (elements and faces associated to a single time group only) and MPI-functionality.
     *
     * First we sort our entities by their time step groups. Every group consists of inner-entities, send-entities and
     * receive entities. Inner-entities are independent from neighboring partitions within one time step of the
     * corresponding time step group. Send-entities are owned by the partition and operate on data which is required
     * by neighboring partitions. Receive-entities are owned by neighboring partitions and contain data which is required
     * by copy-elements. Further the order of the send- and receive-entities is by the neighboring rank's id. If a send-
     * entity is required by more than one time step group-neighbor, it is duplicated. Send- and receive- entities
     * with respect to a a single neighboring rank are sorted by the global ids of the entities. This is to ensure
     * a consistent ordering of the send-entities w.r.t. to the remote receive-ents and vice versa.
     * For vertices and faces, we consider only consider interface-data to be shared.
     **/
    void setupDataLayout();

    /**
     * Sorts the given array by the global ids of the entities.
     *
     * @param io_ents array of entities which will be replaced with sorted entities.
     **/
    void sortGId( std::vector< moab::EntityHandle > &io_ents );

    /**
     * Sorts the given array of entities in ascending order with respect to the global id of the faces' vertices.
     *
     * @param io_faHas array of faces which will be replaced with sorted faces.
     **/
    void sortElFaEntities( std::vector< moab::EntityHandle > &io_faHas );

    /**
     * Gets the elements adjacent to elements through vertices.
     * Remark: Order with respect to vertices is ascending.
     *
     * @param o_nElVeEl will be set to number of adjacent elements.
     * @param o_elVeEl will be set to elements adjacent to elements through vertices; optional, nullptr ignores the argument.
     **/
    void getElVeEl( int_el  &o_nElVeEl,
                    int_el **o_elVeEl );

    /**
     * Gets the ids of the face neighboring elements.
     *
     * Remark: Ordering w.r.t. to the global ids of the vertices of the shared face is guaranteed to be ascending.
     *         This includes periodic boundary faces since the MOAB core is queried for this operation.
     *
     * @param i_element handle of the element.
     * @param o_neighboringIds will be set to ids of the face neighbors. std::numeric_limits<int_el>::max() is returned for non-existing neighbors.
     **/
    void getElementFaceNeighbors( moab::EntityHandle i_element,
                                  int_el             o_neighboringIds[C_ENT[T_SDISC.ELEMENT].N_FACES] );

    /**
     * Gets the ids of adjacent faces for a single element.
     * Remark: Ordering of faces is ascending w.r.t. to the vertices.
     *         Periodic faces is internal vertices defined in the MOAB core.
     *
     * @param i_element handle of the element.
     * @param o_elementAdjacentFaces will be set to ids of neighboring faces.
     **/
    void getElementAdjacentFaces( moab::EntityHandle i_element,
                                  int_el             o_elementAdjacentFaces[C_ENT[T_SDISC.ELEMENT].N_FACES] );

    /**
     * Initializes an entity by setting the local id-tag and the mesh-to-data and data-to-mesh mapping.
     *
     * @param i_ents entities which tages will be set.
     * @param o_enMeDa will be set to mesh-to-data mapping.
     * @param o_enDaMe will be set to data-to-mesh mapping.
     **/
    void initEn( const std::vector< moab::EntityHandle > &i_ents,
                       std::vector< int_el >             &o_enMeDa,
                       std::vector< int_el >             &o_enDaMe );

    /**
     * Initializes the MOAB-interface by creating the local id-tag and setting its values.
     * Remark: This should be called after we created the data layout.
     **/
    void init();

    /**
     * Prints MOAB information.
     **/
    void printMoab();

  public:
    /**
     * Sets up a new MOAB interface.
     *
     * @param i_dim dimensions of the interface
     * @param i_nBndVals number of boundary values in the mesh.
     * @param i_bndVals bnd regions in the mesh. This is whenever a face has no 2nd neighbor and bnd conditions apply.
     * @param i_periodic tag value of periodic boundaries (if any are present).
     **/
    Moab(       unsigned int  i_dim,
                unsigned int  i_nBndVals    = 0,
          const int          *i_bndVals     = NULL,
                int           i_periodicVal = std::numeric_limits<int>::max() );

    /**
     * Reads a given mesh.
     *
     * @param i_pathToMesh path to mesh.
     * @param i_optRead read options forwarded to MOAB.
     **/
    void read( const std::string &i_pathToMesh,
               const std::string &i_optRead );

    /**
     * Writes the mesh, including all information added during processing.
     *
     * @param i_pathToMesh path to mesh.
     **/
    void write( const std::string &i_pathToMesh );

    /**
     * Gets the data layout of the vertices (0).
     *
     * @return data layout of the vertices.
     **/
    t_enLayout getVeLayout() const { return m_veLayout; };

    /**
     * Gets the data layout of the faces (N_DIM-1).
     *
     * @return data layout of the faces.
     **/
    t_enLayout getFaLayout() const { return m_faLayout; };

    /**
     * Gets the data layout of the elements (N_DIM).
     *
     * @return data layout of the elements.
     **/
    t_enLayout getElLayout() const { return m_elLayout; };

    /**
     * Gets the mapping of indices from mesh to data and data to mesh.
     *
     * @return index map.
     **/
    const t_inMap* getInMap() const { return &m_inMap; };

    /**
     * Gets the number of vertices.
     *
     * Remark: This includes receive-entities and duplicated send-entities.
     *
     * @return number of vertices.
     **/
    int_el getNVertices() const { return m_vertices.size(); };

    /**
     * Gets the number of faces.
     *
     * Remark: This includes receive-entities and duplicated send-entities.
     *
     * @return number of faces.
     **/
    int_el getNFaces() const { return m_faces.size(); };

    /**
     * Gets the number of elements
     *
     * Remark: This includes receive-entities and duplicated send-entities.
     *
     * @return number of elements.
     **/
    int_el getNElements() const { return m_elements.size(); };

    /**
     * Gets the regions defined by the mesh.
     * This includes all entitites (e.g. boundary conditions for faces).
     *
     * @param o_values vector to which the values will be written.
     **/
     void getMeshRegions( std::vector<int> &o_values );

    /**
     * Gets the local ids of all adjacent vertices for all elements.
     *
     * Remark 1: This includes "redundant" vertices of periodic boundaries.
     * Remark 2: Ordering of vertices' ids is ascending.
     *
     * @param o_elementAdjacentVertices vertices adjacent to the elements.
     **/
    void getElVe( int_el (*o_elVe)[C_ENT[T_SDISC.ELEMENT].N_VERTICES] );

    /**
     * Gets the local ids of adjacent faces for all elements.
     *
     * Remark: Order with respect to the faces' vertices is ascending.
     *
     * @param o_elementsAdjacentFaces will be set to ids of all elements' adjacent faces.
     **/
    void getElementsAdjacentFaces( int_el (*o_elementsAdjacentFaces)[C_ENT[T_SDISC.ELEMENT].N_FACES] );


    /**
     * Gets the elements adjacent to elements through vertices.
     * Remark: Order with respect to vertices is ascending.
     *
     * @param o_elVeEl will be set to element adjacent to elements through vertices.
     **/
    void getElVeEl( int_el **o_elVeEl );

    /**
     * Gets the number of elements adjacent to the elements through vertices.
     *
     * @return number of adjacent elements (vertices as bridge).
     **/
    int_el getNelVeEl();

    /**
     * Gets the local ids of all face neighbors of all elements.
     *
     * Remark: Order with respect to the faces' vertices is ascending.
     *
     * @o_neighboringIds array to which face-neighboring elements' ids will be written.
     **/
    void getElementsFaceNeighbors( int_el (*o_neighboringIds)[C_ENT[T_SDISC.ELEMENT].N_FACES] );

    /**
     * Gets the local ids of adjacent elements for all faces.
     *
     * Remark 1: If only one neighbor exists, the second element [][1] is set to: std::numeric_limits<int_el>::max().
     * Remark 2: The ordering of the adj. elements is guaranteed to be ascending.
     *
     * @param o_neighboringIds will be set to neighboring element ids for all faces.
     **/
    void getFacesAdjacentElements( int_el (*o_neighboringIds)[2] );

    /**
     * Gets the local ids of adjacent vertices for all faces.
     *
     * Remark: Ordering of vertices is ascending.
     *
     * @param o_neighVeIds will be set to ids of the faces' adjacent vertices.
     **/
    void getFacesAdjacentVertices( int_el (*o_neighVeIds)[C_ENT[T_SDISC.FACE].N_VERTICES] );

    /**
     * Gets the global ids of the vertices.
     *
     * @param o_gIds will be set to reference to the global vertex ids.
     **/
    void getGIdsVe( std::vector< int_gid > &o_gIds ) const;

    /**
     * Gets the global ids of the faces.
     *
     * @param o_gIds will be set to reference to the global face ids.
     **/
    void getGIdsFa( std::vector< int_gid > &o_gIds ) const;

    /**
     * Gets the global ids of the elements.
     *
     * @param o_gIds will bet set to reference to the global element ids.
     **/
    void getGIdsEl( std::vector< int_gid > &o_gIds ) const;

    /**
     * Gets the the names of the tags in the mesh.
     * The positions in the vector are identical to local ids for querying the tags.
     *
     * @param o_tagsNames will be set to names of tags.
     **/
    void getTagsNames( std::vector< std::string > &o_tagsNames );

    /**
     * Synchronizes the given tags for shared entitities in distributed memory settings.
     *
     * @param i_tids local ids of the tags.
     * @param i_nDim number of entity dimensions. 0: vertex, m_dim-1: face, m_dim: element.
     **/
    void syncTags( std::vector< unsigned short > const &i_tids,
                   unsigned short                       i_nDim  );

    /**
     *  Gets the number of bytes for a tag.
     *
     *  @param i_tid local id of the tag.
     *  @return size of the tag in bytes.
     **/
    int getTagBytes( unsigned short i_tid );

    /**
     * Gets tag data for the given entity.
     *
     * @param i_tid local id of the tag.
     * @param i_nDim number of entity dimensions. 0: vertex, m_dim-1: face, m_dim: element.
     * @param i_en local id of the entity.
     * @param o_val will be set to the value of the tag.
     **/
    void getTagData( unsigned short  i_tid,
                     unsigned short  i_nDim,
                     int_el          i_en,
                     void           *o_val );

    /**
     * Gets the material tag value for the given entity.
     *
     * @param i_dim dimensions of the entity: 0: vertex, m_dim-1: face, m_dim: element
     * @param i_ent local id of the entity.
     * @param o_val will be set to the value of the material tag.
     **/
    void getMatVal( unsigned short  i_dim,
                    int_el          i_ent,
                    int_spType     &o_val );

    /**
     * Gets the coordinates of a vertex.
     *
     * @param i_vertex local id of the vertex.
     * @param o_x position in x-direction.
     * @param o_y position in y-direction.
     * @param o_z position in z-direction.
     **/
    void getVeCoords( int_el  i_vertex,
                      double &o_x,
                      double &o_y,
                      double &o_z );

    /**
     * Gets the coordinates of the face's vertices.
     *
     * @param i_face local id of the face.
     * @param o_x position of vertices in x-direction.
     * @param o_y position of vertices in y-direction.
     * @param o_z position of vertices in z-direction.
     **/
    void getFaceVerticesCoords( int_el i_face,
                                double o_x[C_ENT[T_SDISC.FACE].N_VERTICES],
                                double o_y[C_ENT[T_SDISC.FACE].N_VERTICES],
                                double o_z[C_ENT[T_SDISC.FACE].N_VERTICES] );

    /**
     * Gets the coordinates of the element's vertices.
     *
     * @param i_element local id of the element.
     * @param o_x position of vertices in x-direction.
     * @param o_y position of vertices in y-direction.
     * @param o_z position of vertices in z-direction.
     **/
    void getElementVerticesCoords( int_el    i_element,
                                   double o_x[C_ENT[T_SDISC.ELEMENT].N_VERTICES],
                                   double o_y[C_ENT[T_SDISC.ELEMENT].N_VERTICES],
                                   double o_z[C_ENT[T_SDISC.ELEMENT].N_VERTICES] );

    /**
     * Gets the coordinates of a normal point with respect to the face.
     * This is the remaining vertex of the respective element for triangles and tets.
     *
     * Remark 1: We return true or false to decide on which side we are.
     *           We can not simply return the coords of the left elements additional point due to periodic boundaries.
     *
     * Remark 2: Boundary elements don't have a neighbor and are on the "left" side always.
     *
     * @param i_face id of the face.
     * @param o_left will be set to true if the element whose point is returned has smaller id than the second element.
     * @param o_x will be set to x coordinate.
     * @param o_y will be set to y coordinate.
     * @param o_z will be set to z coordinate.
     **/
    void getNormalPointCoords( int_el  i_face,
                               bool   &o_left,
                               double &o_x,
                               double &o_y,
                               double &o_z );
};

#endif
