/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2016, Regents of the University of California
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
 * Unstructured mesh representation.
 **/
#ifndef UNSTRUCTURED_H_
#define UNSTRUCTURED_H_

#include "constants.hpp"
#include <string>

#ifdef PP_USE_MOAB
#include "Moab.h"
#endif

namespace edge {
  namespace mesh {
    class Unstructured;
  }
}

class edge::mesh::Unstructured {
  //private:

#ifdef PP_USE_MOAB
    Moab m_moab;
#endif

    /**
     * Prints the mesh statistics.
     **/
    void printStats();

    /**
     * Initializes the given array with the adjacenet vertices of all elements.
     * Remark 1: This includes "redundant" vertices on periodic faces, which is different from the face handling (only considered once).
     * Remark 2: Resulting order of the per-element vertices is guaranteed to be ascending.
     *           E.g. if an element has vertices 21, 42, 90, 20, the resulting array will be [20, 21, 42, 90].
     *
     * @param o_elVe will be set to ids of all elements' adjacent vertices.
     **/
    void getElVe( int_el (*o_elVe)[C_ENT[T_SDISC.ELEMENT].N_VERTICES] );

    /**
     * Initializes the given array with ids of adjacent faces of all elements.
     * Remark 1: Periodic faces are removed in this process. Two elements neighboring through a periodic face
     *           will return the same face. While this is natural, it is different from the mesh representation
     *           and from the handling of vertices, which contain mesh-specified entries independent of periodic boundaries.
     * Remark 2: Resulting order of faces is ascending w.r.t. to their vertex ids.
     *           For periodic boundaries the vertices within the original non-modified mesh is considered.
     *
     * @param o_elFa will be set to ids of all elements' adjacent faces.
     **/
    void getElementsAdjacentFaces( int_el (*o_elFa)[C_ENT[T_SDISC.ELEMENT].N_FACES] );

    /**
     * Initializes the given array with the adjacent elements of all faces.
     * Remark: Resulting order of the elements is guaranteed to be ascending.
     *
     * @param o_neighboringIds will be set to neighboring element ids of all faces.
     **/
    void getFacesAdjacentElements( int_el (*o_neighboringIds)[2] );

    /**
     * Initializes the given array with the adjacent vertices of all faces.
     * Remark: Resulting order of vertices is ascending.
     *
     * @param o_neighVeIds will be set to the ids of all faces' adjacent vertices.
     **/
    void getFacesAdjacentVertices( int_el (*o_neighVeIds)[C_ENT[T_SDISC.FACE].N_VERTICES] );


    /**
     * Determines the element-to-element adjacency through vertices as bridge.
     * Remark: It is expected that o_elVeEl is an array of pointers of size #nel+1 (incl. ghost element) for consistent size computations.
     *         Futher, o_elVeEl[0] is expected to point to raw memory, large enough to hold all adjacent elemenets of all elements.
     *
     * @param o_elVeEl will be set to elements adjacent to elements (vertices as bridge). [*][]: ptr to elements, [][*]: adjacent elements.
     **/
    void getElVeEl( int_el **o_elVeEl );

    /**
     * Initializes the given array with the face neighbors' ids of all elements.
     * Remark: Ordering w.r.t. to the vertices of the shared face is guaranteed to be ascending.
     *         "Redundant" faces are also respected internally and can not be reproduced from face adjacent vertices,
     *          only one of the redundant faces is returned by the respective calls.
     *
     * @param o_neighboringIds will be set to ids of all elements' face neighbors.
     **/
    void getElementsFaceNeighbors( int_el (*o_neighboringIds)[C_ENT[T_SDISC.ELEMENT].N_FACES] );

    /**
     * Computes the volume of an element.
     *
     * @param i_vertices vertices of the element.
     * @param o_volume will bet set to the volume of the element.
     **/
    void computeElementVolume( const double     i_vertices[3][C_ENT[T_SDISC.ELEMENT].N_VERTICES],
                                     real_mesh &o_volume ) const;

    /**
     * Computes the insphere diameter of an element.
     *
     * @param i_vertices vertices of an element.
     * @param o_dia will be set to insphere diameter.
     **/
    void computeElementInDiameter( double     i_vertices[3][C_ENT[T_SDISC.ELEMENT].N_VERTICES],
                                   real_mesh &o_dia ) const;

  public:
    /**
     * Sets up a new unstructured mesh.
     *
     * @param i_nBndVals number of boundary values.
     * @param i_bndVals values of the boundary regions.
     * @param i_periodic tag value of periodic boundaries (if any are present).
     **/
    Unstructured(       unsigned int  i_nBndVals    = 0,
                  const int          *i_bndVals     = NULL,
                        int           i_periodicVal = std::numeric_limits<int>::max() );

    /**
     * Reads a given unstructured mesh.
     *
     * @param i_pathToMesh path to mesh.
     * @param i_optRead read options forwarded to the meshing interface.
     **/
    void read( const std::string &i_pathToMesh,
               const std::string &i_optRead );

    /**
     * Writes the processed, unstructured mesh.
     *
     * @param i_pathToMesh
     **/
    void write( const std::string &i_pathToMesh );

    /////////////////////////////////////////////////////////
    // TODO: Fix assignment of global IDs, this is a copy. //
    /////////////////////////////////////////////////////////
    /**
     * Gets the global ids of the vertices.
     **/
    void getGIdsVe( std::vector< int_gid > &o_gIds ) const;

    /**
     * Gets the global ids of the faces.
     **/
    void getGIdsFa( std::vector< int_gid > &o_gIds ) const;

    /**
     * Gets the global ids of the elements.
     **/
    void getGIdsEl( std::vector< int_gid > &o_gIds ) const;

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
     * Gets the number of vertices in the mesh (N_DIM-2).
     *
     * @return number of elements.
     **/
    int_el getNVertices() const;

    /**
     * Gets the number of faces in the mesh (N_DIM-1).
     *
     * @return number of faces.
     **/
    int_el getNFaces() const;

    /**
     * Gets the number of elements in the mesh (N_DIM).
     *
     * @return number of elements.
     **/
    int_el getNElements() const;

    /**
     * Gets the number of elements adjacent to all elements (vertices as bridge).
     *
     * @return number of adjacent elements.
     **/
    int_el getNelVeEl();

    /**
     * Gets the connectivity information.
     *  elVe:   element -> vertices
     *    triangles are guaranteed to have a counter-clockwise ordering of the vertices.
     *    tets are guaranteed to have a counter-clockwise order of vertices
     *    2,3,4 w.r.t. to the normal pointing towards vertex 1.
     *  elFa:   element -> faces
     *  faEl:   face    -> elements
     *  elFaEl: element -> faces -> elements (excluding starting element)
     *
     * Remark: Faces on periodic boundaries show up only once (in no particular order).
     *         Vertices on periodic boundaries show up reduantly matching the mesh.
     *
     * @param i_veChars vertex chars, will be used to determine the mapping to the reference element.
     * @param i_faChars face chars, will be used to identify boundary faces in consistency checks.
     * @param o_connect will be set to connectity information.
     **/
    void getConnect( const t_vertexChars *i_veChars,
                     const t_faceChars   *i_faChars,
                           t_connect     &o_connect );

    /**
     * Gets the vertex characteristics.
     *
     * @param o_veChars will be set to the vertices' characteristics.
     **/
    void getVeChars( t_vertexChars *o_veChars );

    /**
     * Gets the characteristics of the elements.
     *
     * @param o_elChars will be set to the elements' characteristics
     **/
    void getElChars( t_elementChars *o_elChars );

    /**
     * Gets the characteristics of the faces.
     *
     * @param o_faChars will be set to the faces' characteristics.
     **/
    void getFaChars( t_faceChars *o_faChars );

    /**
     * Gets the given parameters for all dense entities (including duplicates).
     *
     * @param i_nDim number of dimensions of the dense entities. 0: vertices, NDIM-1: faces, NDIM: elemennts, where NDIM is the number of dimensions for the simulation.
     * @param i_nPars number of parameters, which are queried.
     * @param i_names names of the paramters.
     * @param o_pars output to which the parameters are written as array of struct. Fastest dimension are parameters.
     * @return error code. 0 if succecessful, 1 if parameter is missing, 2 if size of template parameter does not match mesh specifications (size in byte).
     *
     * @paramt TL_T_PAR type of the parameter.
     **/
    template< typename TL_T_PAR >
    bool getParsDe( unsigned short   i_nDim,
                    unsigned short   i_nPars,
                    std::string    * i_names,
                    TL_T_PAR       * o_pars ) {
      // get the names of available tags
      std::vector< std::string > l_tagsNames;
#ifdef PP_USE_MOAB
      m_moab.getTagsNames( l_tagsNames );
#else
      EDGE_LOG_FATAL;
#endif

      // determine local ids of the parameters
      std::vector< unsigned short > l_idsPar;
      for( unsigned short l_pa = 0; l_pa < i_nPars; l_pa++ ) {
        // iterate over mesh tag-names and find matches
        for( std::size_t l_pm = 0; l_pm < l_tagsNames.size(); l_pm++ ) {
          // store local id if this is our parameter
          if( i_names[l_pa] == l_tagsNames[l_pm] ) {
            l_idsPar.push_back( l_pm );
            break;
          }

          // abort execution if the parameter does not exist in the mesh
          if( l_pm == l_tagsNames.size()-1 ) {
            return 1;
          }
        }
      }

      // check that the template parameters size matches the mesh's specs
#ifdef PP_USE_MOAB
      for( unsigned short l_pa = 0; l_pa < i_nPars; l_pa++ ) {
        int l_nParBytes = m_moab.getTagBytes( l_idsPar[l_pa] );

        EDGE_CHECK_EQ(l_nParBytes, sizeof( TL_T_PAR ) )
          << "size of mesh parameter " << l_tagsNames[ l_idsPar[l_pa] ]
          << " and template don't match; "
          << "mesh: " << l_nParBytes
          << ", template: " << sizeof( TL_T_PAR );
      }
#else
      EDGE_LOG_FATAL;
#endif

      // synchronize tags
#ifdef PP_USE_MOAB
      m_moab.syncTags( l_idsPar, N_DIM );
#else
      EDGE_LOG_FATAL;
#endif

      // determine number of entities
      int_el l_nEn = ( i_nDim == 0       ) ? getNVertices() :
                     ( i_nDim == N_DIM-1 ) ? getNFaces() :
                                             getNElements();
      // iterate over entities
      for( int_el l_en = 0; l_en < l_nEn; l_en++ ) {
        // iterate over parameters
        for( unsigned short l_pa = 0; l_pa < i_nPars; l_pa++ ) {
#ifdef PP_USE_MOAB
          // get data from MOAB
          m_moab.getTagData( l_idsPar[l_pa],
                             i_nDim,
                             l_en,
                             o_pars+(i_nPars*l_en)+l_pa );
#else
          EDGE_LOG_FATAL;
#endif
        }
      }

      return 0;

    }
};

#endif
