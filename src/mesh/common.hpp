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
 * Common functions for the mesh
 **/
#ifndef EDGE_MESH_COMMON_HPP
#define EDGE_MESH_COMMON_HPP

#include "linalg/Geom.hpp"
#include "io/logging.h"
#include "monitor/instrument.hpp"
#include "constants.hpp"
#include "data/EntityLayout.type"
#include "linalg/Matrix.h"
#include <cassert>
#include <cmath>
#include <limits>

namespace edge {
  namespace mesh {
    template < t_entityType TL_T_EL >
    class common;
  }
}

template < t_entityType TL_T_EL >
class edge::mesh::common {
  private:
    //! number of dimensions
    static const unsigned short TL_N_DIS = C_ENT[TL_T_EL].N_DIM;

    //! face type
    static const t_entityType TL_T_FA = C_ENT[TL_T_EL].TYPE_FACES;

    //! number of face vertices
    static unsigned short const TL_N_FA_VES = C_ENT[TL_T_FA].N_VERTICES;

    //! number of element vertices
    static unsigned short const TL_N_EL_VES = C_ENT[TL_T_EL].N_VERTICES;

    //! number of element faces
    static unsigned short const TL_N_EL_FAS = C_ENT[TL_T_EL].N_FACES;

    /**
     * Sets global vertex ids to dummy values for use in periodic boundaries.
     *
     *  Remark: In the case of periodic boundaries we have unique, non-shared vertices for the elements.
     *          Therefore we can not compate vertex ids to find corresponding node mappings of faces.
     *          Thus we go through the geometric location (periodic boundaries are aligned to either x-, y- or z-coord).
     *
     * @param i_gIdsVe global ids of the vertices.
     * @param i_elIdLoc id of the local element.
     * @param i_elIdNeg if of the neighboring element.
     * @param i_localFace local face of the element (0-4).
     * @param i_elVe elements' vertices.
     * @param i_veChars vertex characteristics.
     * @param o_vIds will be set to dummy global ids of vertices containing the local elements vertex. all but the dominant (w.r.t. to the vertex ordering) entries are set to std::numeric_limits<int_gid>::max().
     **/
#if defined PP_T_ELEMENTS_TET4
    static void setPeriodicVesTet4( const std::vector< int_gid > &i_gIdsVe,
                                           int_el                 i_elIdLoc,
                                           int_el                 i_elIdNeg,
                                           unsigned short         i_localFace,
                                     const int_el               (*i_elVe)[TL_N_EL_VES],
                                     const t_vertexChars         *i_veChars,
                                           int_gid                o_vIds[4] ) {
      // remark: periodic faces have outer-pointing normals pointing in opposite directions

      // get local face's vertex coords
      real_mesh l_locFaCoords[3][3];
      real_mesh l_normalPt[3];
      for( unsigned short l_di = 0; l_di < 3; l_di++ )
        l_normalPt[l_di] = std::numeric_limits< real_mesh >::max();

      for( unsigned int l_dim = 0; l_dim < 3; l_dim++ ) {
        if( i_localFace == 0 ) {
          l_locFaCoords[l_dim][0] = i_veChars[ i_elVe[i_elIdLoc][0] ].coords[l_dim];
          l_locFaCoords[l_dim][1] = i_veChars[ i_elVe[i_elIdLoc][2] ].coords[l_dim];
          l_locFaCoords[l_dim][2] = i_veChars[ i_elVe[i_elIdLoc][1] ].coords[l_dim];
          l_normalPt[l_dim] = i_veChars[ i_elVe[i_elIdLoc][3] ].coords[l_dim];
        }
        else if( i_localFace == 1 ) {
          l_locFaCoords[l_dim][0] = i_veChars[ i_elVe[i_elIdLoc][0] ].coords[l_dim];
          l_locFaCoords[l_dim][1] = i_veChars[ i_elVe[i_elIdLoc][1] ].coords[l_dim];
          l_locFaCoords[l_dim][2] = i_veChars[ i_elVe[i_elIdLoc][3] ].coords[l_dim];
          l_normalPt[l_dim] = i_veChars[ i_elVe[i_elIdLoc][2] ].coords[l_dim];
        }
        else if( i_localFace == 2 ) {
          l_locFaCoords[l_dim][0] = i_veChars[ i_elVe[i_elIdLoc][0] ].coords[l_dim];
          l_locFaCoords[l_dim][1] = i_veChars[ i_elVe[i_elIdLoc][3] ].coords[l_dim];
          l_locFaCoords[l_dim][2] = i_veChars[ i_elVe[i_elIdLoc][2] ].coords[l_dim];
          l_normalPt[l_dim] = i_veChars[ i_elVe[i_elIdLoc][1] ].coords[l_dim];
        }
        else if( i_localFace == 3 ) {
          l_locFaCoords[l_dim][0] = i_veChars[ i_elVe[i_elIdLoc][1] ].coords[l_dim];
          l_locFaCoords[l_dim][1] = i_veChars[ i_elVe[i_elIdLoc][2] ].coords[l_dim];
          l_locFaCoords[l_dim][2] = i_veChars[ i_elVe[i_elIdLoc][3] ].coords[l_dim];
          l_normalPt[l_dim] = i_veChars[ i_elVe[i_elIdLoc][0] ].coords[l_dim];
        }
        else EDGE_LOG_FATAL << i_localFace;
      }

      // get outer point normal of local face
      real_mesh l_locOutNormal[3];
      for( unsigned short l_di = 0; l_di < 3; l_di++ )
        l_locOutNormal[l_di] = std::numeric_limits< real_mesh >::max();

      linalg::Geom::computeOutPtNormal( TRIA3,
                                        l_locFaCoords[0],
                                        l_normalPt,
                                        l_locOutNormal );

      // iterate over the neighbor's faces and find the matching face
      for( unsigned short l_faNe = 0; l_faNe < 4; l_faNe++ ) {
        // get neighbor face's vertex coords
        real_mesh l_neFaCoords[3][3];
        real_mesh l_neNormalPt[3];
        for( unsigned int l_dim = 0; l_dim < 3; l_dim++ ) {
          if( l_faNe == 0 ) {
            l_neFaCoords[l_dim][0] = i_veChars[ i_elVe[i_elIdNeg][0] ].coords[l_dim];
            l_neFaCoords[l_dim][1] = i_veChars[ i_elVe[i_elIdNeg][2] ].coords[l_dim];
            l_neFaCoords[l_dim][2] = i_veChars[ i_elVe[i_elIdNeg][1] ].coords[l_dim];
            l_neNormalPt[l_dim]    = i_veChars[ i_elVe[i_elIdNeg][3] ].coords[l_dim];
          }
          else if( l_faNe == 1 ) {
            l_neFaCoords[l_dim][0] = i_veChars[ i_elVe[i_elIdNeg][0] ].coords[l_dim];
            l_neFaCoords[l_dim][1] = i_veChars[ i_elVe[i_elIdNeg][1] ].coords[l_dim];
            l_neFaCoords[l_dim][2] = i_veChars[ i_elVe[i_elIdNeg][3] ].coords[l_dim];
            l_neNormalPt[l_dim]    = i_veChars[ i_elVe[i_elIdNeg][2] ].coords[l_dim];
          }
          else if( l_faNe == 2 ) {
            l_neFaCoords[l_dim][0] = i_veChars[ i_elVe[i_elIdNeg][0] ].coords[l_dim];
            l_neFaCoords[l_dim][1] = i_veChars[ i_elVe[i_elIdNeg][3] ].coords[l_dim];
            l_neFaCoords[l_dim][2] = i_veChars[ i_elVe[i_elIdNeg][2] ].coords[l_dim];
            l_neNormalPt[l_dim]    = i_veChars[ i_elVe[i_elIdNeg][1] ].coords[l_dim];
          }
          else if( l_faNe == 3 ) {
            l_neFaCoords[l_dim][0] = i_veChars[ i_elVe[i_elIdNeg][1] ].coords[l_dim];
            l_neFaCoords[l_dim][1] = i_veChars[ i_elVe[i_elIdNeg][2] ].coords[l_dim];
            l_neFaCoords[l_dim][2] = i_veChars[ i_elVe[i_elIdNeg][3] ].coords[l_dim];
            l_neNormalPt[l_dim]    = i_veChars[ i_elVe[i_elIdNeg][0] ].coords[l_dim];
          }
        }

       // get outer point normal of neighboring face
       real_mesh l_neOutNormal[3];
       linalg::Geom::computeOutPtNormal( TRIA3,
                                         l_neFaCoords[0],
                                         l_neNormalPt,
                                         l_neOutNormal );

       // check if we have opposite direction out normals
       if(    std::abs( l_locOutNormal[0] + l_neOutNormal[0] ) < TOL.MESH
           && std::abs( l_locOutNormal[1] + l_neOutNormal[1] ) < TOL.MESH
           && std::abs( l_locOutNormal[2] + l_neOutNormal[2] ) < TOL.MESH ) {
         // reset neighbor
         for( unsigned short l_ve = 0; l_ve < 4; l_ve++ ) {
           o_vIds[l_ve] = std::numeric_limits<int_gid>::max();
         }

         // assign dominant vertex
         int_el l_domTmp;
         if( i_localFace <= 2 ) l_domTmp = i_elVe[i_elIdLoc][0];
         else                   l_domTmp = i_elVe[i_elIdLoc][1];

         // get local pos of dominant vertex
         for( unsigned short l_ve = 0; l_ve < 3; l_ve++ ) {
           unsigned int short l_nTmpMatches = 0;
           for( unsigned int l_dim = 0; l_dim < 3; l_dim++ ) {
             if( std::abs( l_neFaCoords[l_dim][l_ve] - i_veChars[l_domTmp].coords[l_dim] ) < TOL.MESH ) {
               l_nTmpMatches++;
             }
           }
           assert( l_nTmpMatches != 3 );

           // the 'dominant' vertex has two share two dimensions with the neighbor
           if( l_nTmpMatches == 2 ) {
             // assign the vertex
             if( l_faNe == 0 ) {
               if(      l_ve == 0 ) o_vIds[0] = i_gIdsVe[ l_domTmp ];
               else if( l_ve == 1 ) o_vIds[2] = i_gIdsVe[ l_domTmp ];
               else if( l_ve == 2 ) o_vIds[1] = i_gIdsVe[ l_domTmp ];
             }
             else if( l_faNe == 1 ) {
               if(      l_ve == 0 ) o_vIds[0] = i_gIdsVe[ l_domTmp ];
               else if( l_ve == 1 ) o_vIds[1] = i_gIdsVe[ l_domTmp ];
               else if( l_ve == 2 ) o_vIds[3] = i_gIdsVe[ l_domTmp ];
             }
             else if( l_faNe == 2 ) {
               if(      l_ve == 0 ) o_vIds[0] = i_gIdsVe[ l_domTmp ];
               else if( l_ve == 1 ) o_vIds[3] = i_gIdsVe[ l_domTmp ];
               else if( l_ve == 2 ) o_vIds[2] = i_gIdsVe[ l_domTmp ];
             }
             else if( l_faNe == 3 ) {
               if(      l_ve == 0 ) o_vIds[1] = i_gIdsVe[ l_domTmp ];
               else if( l_ve == 1 ) o_vIds[2] = i_gIdsVe[ l_domTmp ];
               else if( l_ve == 2 ) o_vIds[3] = i_gIdsVe[ l_domTmp ];
             }
             else assert( false );
           }
         }
       }
     }
   }
#endif

  public:
    /**
     * Computes the barycenter of an element.
     *
     * @param i_vertices vertices of the element.
     * @param o_baryCtr will be set to barycenter of the element.
     **/
    static void computeBaryCtr( const real_mesh i_vertices[3][TL_N_EL_VES],
                                      real_mesh o_baryCtr[3] ) {
#if defined PP_T_ELEMENTS_LINE || defined PP_T_ELEMENTS_QUAD4R || defined PP_T_ELEMENTS_TRIA3 || defined PP_T_ELEMENTS_HEX8R
      // reset baryctr;
      o_baryCtr[0] = o_baryCtr[1] = o_baryCtr[2] = 0;

      // derive  bary ctr
      for( int_md l_ve = 0; l_ve < TL_N_EL_VES; l_ve++ ) {
        o_baryCtr[0] += i_vertices[0][l_ve] / TL_N_EL_VES;
        o_baryCtr[1] += i_vertices[1][l_ve] / TL_N_EL_VES;
        o_baryCtr[2] += i_vertices[2][l_ve] / TL_N_EL_VES;
      }
#else
    assert(false);
#endif
    }

    /**
     * Normalizes the connectivity information of triangular meshes.
     * Remark 1: After the normalization the connectivity is not necessarily ordered w.r.t. to vertex ids.
     * Remark 2: We expect sorted input of all connectitvity information.
     *           No additional checks or changes are performed for this.
     * Remark 3: The mapping to the reference elements is unique by the choice of the first vertex
                 and a counter-clockwise orientation.
     *           Local id of adjacent faces and elements is uniquely given by the mapping.
     * Remark 4: The term "clockwise" is determined by the determinant of the vertices' matrix.
     * Remark 5: Faces follow our counter-clockwise storage of the vertices.
     *
     *
     * Example for the mapping of the initial, vertex-ordered face assignment to an oriantation-based assignment:
     *
     *                 15                             15
     *                 *                              *
     *               *    *             map         *    *
     *             * 0    1  *         ---->      * 0     2 *
     *           *      2       *               *       1      *
     *       35 ******************** 80     35 ******************** 80
     *
     * Example for a valid mapping (0->1, 1->2, 2-> 0 is counter-clockwise):
     *
     *                0                   2 *
     *                *              Ref.   * *
     *             *      *         ---->   *   *
     *           *             *            *     *
     *       1 ******************** 2     0 ********* 1
     *
     * In this case we only check that the first face has vertices 0 and 1, the second face has vertices 1 and 2 and the
     * third face has vertices 2 and 0. No changes are performed.
     *
     * Example for an invalid mapping (0->1, 1->2, 2-> 0 is clockwise):
     *
     *                0                   1 *
     *                *              Ref.   * *
     *             *      *         ---->   *   *
     *           *             *            *     *
     *       2 ******************** 1     0 ********* 2
     *
     * In this case we exchange the local position of vertices 1 and 2.
     * Additional we change the ordering of the faces to match this.
     *
     * @param i_elLayout data layout of the elements.
     * @param i_veChars characteristics of vertices.
     * @param io_elVe vertices adjacent to the elements, sorted for every element by the vertex ids.
     * @param io_elFa faces adjacent to the elements, sorted for every element by the faces' vertex ids.
     * @param io_elFaEl elements adjacent to the elements, sorted for every (local) element by the vertex ids of the shared face.
     **/
    static void normOrdTria( const t_enLayout     &i_elLayout,
                             const t_vertexChars  *i_veChars,
                                   int_el        (*io_elVe)[3],
                                   int_el        (*io_elFa)[3],
                                   int_el        (*io_elFaEl)[3] ) {
      // iterate over all elements
      for( int_el l_el = 0; l_el < i_elLayout.nEnts; l_el++ ) {
        /*
         * reorder face-information in ascending order of vertex ids:
         * id0--fa0-->id1--fa1-->id2--fa2-->
         */
        int_el l_faTmp = io_elFa[l_el][2];
        io_elFa[l_el][2] = io_elFa[l_el][1];
        io_elFa[l_el][1] = l_faTmp;

        int_el l_elTmp = io_elFaEl[l_el][2];
        io_elFaEl[l_el][2] = io_elFaEl[l_el][1];
        io_elFaEl[l_el][1] = l_elTmp;

        /*
         * Too easy.. let's work on the counter-clockwise ordering.
         */

        // get the coords of the thre verts
        real_mesh l_veCoords[2][3];

        for( unsigned int l_ve = 0; l_ve < 3; l_ve++ ) {
          int_el l_veId = io_elVe[l_el][l_ve];

          // set the coords, we ignore the third dimension
          l_veCoords[0][l_ve] = i_veChars[l_veId].coords[0];
          l_veCoords[1][l_ve] = i_veChars[l_veId].coords[1];
          assert( std::abs( i_veChars[l_veId].coords[2] ) < TOL.MESH );
        }

        // vector pointing from vertex 0 to 1
        real_mesh l_0to1[2];
        l_0to1[0] = l_veCoords[0][1] - l_veCoords[0][0];
        l_0to1[1] = l_veCoords[1][1] - l_veCoords[1][0];

        // vector pointing from vertex 2 to 0
        real_mesh l_2to0[2];
        l_2to0[0] = l_veCoords[0][2] - l_veCoords[0][0];
        l_2to0[1] = l_veCoords[1][2] - l_veCoords[1][0];


        // compute the determinant
        real_mesh l_det  = l_0to1[0] * l_2to0[1]; // a*d
                  l_det -= l_2to0[0] * l_0to1[1]; // -b*c

        // assert our vectors are linear independent
        assert( std::abs( l_det ) > TOL.MESH );

        // negative determinant -> clockwise -> change pos of 2nd and 3rd vertex
        if( l_det < 0 ) {
          int_el l_veTmp = io_elVe[l_el][1];
          io_elVe[l_el][1] = io_elVe[l_el][2];
          io_elVe[l_el][2] = l_veTmp;

          /*
           * change position of face-information accordingly
           *
           *              v0                                v0
           *               *                                *
           *          f2 *     * f0      ---->        f0  *     * f2
           *           *    f1    *                    *      f1    *
           *      v2 ****************** v1         v1 ****************** v2
           *
           * Swapping positions of vertices v1 and v2 means that we have to swap positions
           * of faces f0 and f2 also.
           *
           */
          l_faTmp = io_elFa[l_el][0];
          io_elFa[l_el][0] = io_elFa[l_el][2];
          io_elFa[l_el][2] = l_faTmp;

          l_elTmp = io_elFaEl[l_el][0];
          io_elFaEl[l_el][0] = io_elFaEl[l_el][2];
          io_elFaEl[l_el][2] = l_elTmp;
        }
      }
    }

    /**
     * Enforce a consistent ordering for tetrahedral elements.
     *
     * Given the tetrahedral vertices v0-v3. We keep the assignment of vertices v0 and v1.
     * Then we consider outer pointing normal of face v1-v2-v3.
     * Looking from the outside to this face (opposite direction of the normal),
     * we enforce counter-clockwise storage of vertices v2 and v3 w.r.t. to face 3.
     *
     * While exchanging v2 and v3 leaves the nodes (not the ordering) of face 2 (0-3-2)
     * and face 3 (1-2-3) untouched, we have to change the exchange the positions of face 0 and 1:
     *   face 0 (0-2-1) gets face 1 (0-1-3).
     *   face 1 (0-1-3) gets face 0 (0-2-2).
     *
     *
     * Example:
     *
     *                 v2                                          v3
     *                 *                                           *
     *               * . *                                       * . *
     *             *    .  *                                   *    .  *
     *           *       .   *                               *       .   *
     *         *          .    *                           *          .    *
     *       *  .            .   *              ->       *  .            .   *
     *  v1 *                    .  *                v1 *                    .  *
     *            *                . *                        *                . *
     *                    *           .*                              *           .*
     *                            *      *                                    *      *
     *                                     *                                            *
     *                                       v3                                          v2
     *
     * @param i_elLayout data layout of the elements.
     * @param i_veChars characteristics of vertices.
     * @param io_elVe vertices adjacent to the elements, sorted for every element by the vertex ids.
     * @param io_elFa faces adjacent to the elements, sorted for every element by the faces' vertex ids.
     * @param io_elFaEl elements adjacent to the elements, sorted for every (local) element by the vertex ids of the shared face.
     **/
    static void normOrdTet4( const t_enLayout     &i_elLayout,
                             const t_vertexChars  *i_veChars,
                                   int_el        (*io_elVe)[4],
                                   int_el        (*io_elFa)[4],
                                   int_el        (*io_elFaEl)[4] ) {
      // iterate over all elements
      for( int_el l_el = 0; l_el < i_elLayout.nEnts; l_el++ ) {
        // get the local ids of the vertices
        int_el l_vIds[4];
        for( unsigned short l_ve = 0; l_ve < 4; l_ve++ ) {
          l_vIds[l_ve] = io_elVe[l_el][l_ve];
        }

        // get the vertex coords of the tet
        real_mesh l_veCoords[3][4];

        for( unsigned short l_dim = 0; l_dim < 3; l_dim++ ) {
          for( unsigned short l_ve = 0; l_ve < 4; l_ve++ ) {
            l_veCoords[l_dim][l_ve] = i_veChars[l_vIds[l_ve]].coords[l_dim];
          }
        }

        // derive the three vectors pointing from 0 to face 3's vertices
        real_mesh l_directed[3][3];

        for( unsigned short l_dim = 0; l_dim < 3; l_dim++ ) {
          l_directed[l_dim][0] = l_veCoords[l_dim][1] - l_veCoords[l_dim][0];
          l_directed[l_dim][1] = l_veCoords[l_dim][2] - l_veCoords[l_dim][0];
          l_directed[l_dim][2] = l_veCoords[l_dim][3] - l_veCoords[l_dim][0];
        }

        // compute determinant
        real_mesh l_det = linalg::Matrix::det( l_directed );

        // assert non-planar vertices
        assert( std::abs( l_det ) > TOL.MESH );

        // exchange vertices 2,3 and faces 0,1 if we are clockwise
        if( l_det < 0 ) {
          // exchange nodes
          io_elVe[l_el][2] = l_vIds[3];
          io_elVe[l_el][3] = l_vIds[2];

          // exchange faces
          int_el l_tmp;
          l_tmp = io_elFa[l_el][0];
          io_elFa[l_el][0] = io_elFa[l_el][1];
          io_elFa[l_el][1] = l_tmp;

          l_tmp = io_elFaEl[l_el][0];
          io_elFaEl[l_el][0] = io_elFaEl[l_el][1];
          io_elFaEl[l_el][1] = l_tmp;
        }
      }
    }

    /**
     * Checks the consistency of the mesh's adjacency information.
     * Are our face-neighboring relations consistent with element neighbors, especially in terms of local locations.
     * Example:
     *
     *            * * * * * * * Two neighboring elements a and b, they share edge x.
     *          *  x         *  a's local face going to b is 0.
     *        *     x   a   *   b's local face going to c is 2.
     *      *     b  x     *
     *         *      x   *     This means:
     *             *   x *      elementFaceNeighbors[a][0] == b
     *                * *       elementFaceNeighbors[b][2] == a
     *                          elementAdjacentFaces[a][0] == x
     *                          elementAdjacentFaces[b][2] == x
     *                          facesAdjacentElements[0] == a && facesAdjacentElements[1] == b
     *                                                        || (depending on the normal)
     *                          facesAdjacentElements[1] == b && facesAdjacentElements[0] == a
     *
     * @param i_faLayout data layout of the faces.
     * @param i_faceChars characteristics of the faces.
     * @param i_faMeDa mapping of faces: mesh-to-data.
     * @param i_faDaMe mapping of faces: data-to-mesh.
     * @param i_elFaEl ids of elements' face-neighbors.
     * @param i_elFa ids of elements' neighboring faces.
     * @param i_faEl ids of faces' adjacent elements.
     **/
    static void checkConsAdj( const t_enLayout            &i_faLayout,
                              const t_faceChars           *i_faChars,
                              const std::vector< int_el > &i_faMeDa,
                              const std::vector< int_el > &i_faDaMe,
                              const int_el               (*i_elFaEl)[TL_N_EL_FAS],
                              const int_el               (*i_elFa)[TL_N_EL_FAS],
                              const int_el               (*i_faEl)[2] ) {
      // iterate over faces
      for( int_el l_fa = 0; l_fa < i_faLayout.nEnts; l_fa++ ) {
        // continue for face with special handling
        if( (i_faChars[l_fa].spType & MESH_TYPE_NONE) != MESH_TYPE_NONE ) continue;

        // get left and right of this face
        int_el l_elL = i_faEl[l_fa][0];
        int_el l_elR = i_faEl[l_fa][1];

        if( l_elL == l_elR ) {
          EDGE_LOG_FATAL << "same element on both sides of a face, strange: "
                         << l_fa << " " << l_elL << " " << l_elR;
        }

        // find the local position of the face in both elements
        int_el l_faL = std::numeric_limits<int_el>::max();
        int_el l_faR = std::numeric_limits<int_el>::max();

        // get unique id
        int_el l_faUId = i_faMeDa[ i_faDaMe[l_fa] ];

        for( unsigned int l_locFa = 0; l_locFa < TL_N_EL_FAS; l_locFa++ ) {
          if( i_elFa[l_elL][l_locFa] == l_faUId ) l_faL = l_locFa;
          if( i_elFa[l_elR][l_locFa] == l_faUId ) l_faR = l_locFa;
        }

        EDGE_CHECK( l_faL != std::numeric_limits<int_el>::max() );
        EDGE_CHECK( l_faR != std::numeric_limits<int_el>::max() );

        // check that the face neighbor array matches our expectations
        if( i_elFaEl[l_elL][l_faL] != l_elR ) {
          EDGE_LOG_FATAL << "houston, this looks like a bug: "
                         << l_elL << " " << l_faL << " " << i_elFaEl[l_elL][l_faL] << " " << l_elR;
        }
        if( i_elFaEl[l_elR][l_faR] != l_elL ) {
          EDGE_LOG_FATAL << "oops, not good: "
                         << l_elR << " " << l_faR << " " << i_elFaEl[l_elR][l_faR] << " " << l_elL;
        }
      }
    }

    /**
     * Checks the consistency of the ordering.
     *
     * faEl: elements adjacent to faces have to be in ascending order.
     * faVe: vertices adjacent to faces have to be in ascending order.
     *
     * elVe: vertices adjacent to elements have to be in ascending order.
     * elFa + faVe: faces of elements have to be in ascending order w.r.t. to their vertices.
     *
     * @param i_faLayout face layout.
     * @param i_elLayout element layout.
     * @param i_gIdsVe global ids of the vertices.
     * @aram  i_gIdsEl global ids of the elements.
     * @param i_nElements number of elements.
     * @param i_elVe elements' vertices.
     * @param i_elFa elements' faces.
     * @param i_faEl faces' elements.
     * @param i_faVe faces' vertices
     * @param i_periodic true if periodic boundaries exist, limits checks of face-ordering to elements in possession of all their faces' vertices.
     **/
    static void checkConsOrd( const t_enLayout             &i_faLayout,
                              const t_enLayout             &i_elLayout,
                              const std::vector< int_gid > &i_gIdsVe,
                              const std::vector< int_gid > &i_gIdsFa,
                              const std::vector< int_gid > &i_gIdsEl,
                              const int_el                (*i_elVe)[TL_N_EL_VES],
                              const int_el                (*i_elFa)[TL_N_EL_FAS],
                              const int_el                (*i_faEl)[2],
                              const int_el                (*i_faVe)[TL_N_FA_VES],
                                    bool                    i_periodic = true  ) {
       // check size of global ids
       EDGE_CHECK( i_faLayout.nEnts == (int_el) i_gIdsFa.size() );
       EDGE_CHECK( i_elLayout.nEnts == (int_el) i_gIdsEl.size() );

       // iterate over faces and check ascending order of elements and vertices
       for( int_el l_fa = 0; l_fa < i_faLayout.nEnts; l_fa++ ) {
         // check that the first entity exists
         EDGE_CHECK( i_faEl[l_fa][0] < i_elLayout.nEnts );

         // continue for valid adjacent only (bnd conditions)
         if( i_faEl[l_fa][1] != std::numeric_limits< int_el >::max() ) {
           // get global element ids
           int_gid l_gIdsEl[2];
           l_gIdsEl[0] = i_gIdsEl[ i_faEl[l_fa][0] ];
           l_gIdsEl[1] = i_gIdsEl[ i_faEl[l_fa][1] ];

           // check ascending order
           EDGE_CHECK( l_gIdsEl[1] > l_gIdsEl[0] );
         }

         // check order of adjacent verts
         for( unsigned int l_ve = 1; l_ve < TL_N_FA_VES; l_ve++ ) {
           int_gid l_gIdsVe[2];
           l_gIdsVe[0] =  i_gIdsVe[ i_faVe[l_fa][l_ve-1] ];
           l_gIdsVe[1] =  i_gIdsVe[ i_faVe[l_fa][l_ve  ] ];

           EDGE_CHECK( l_gIdsVe[1] > l_gIdsVe[0] );
         }
       }

       // iterate over elements and check order of vertices and faces
       for( int_el l_el = 0; l_el < i_elLayout.nEnts; l_el++ ) {
         /*
          * check ascending order of vertices
          */
         for( unsigned l_ve = 1; l_ve < TL_N_EL_VES; l_ve++ ) {
           // get global ids of the vertices
           int_gid l_gIdsVe[2];
           l_gIdsVe[0] =  i_gIdsVe[ i_elVe[l_el][l_ve-1] ];
           l_gIdsVe[1] =  i_gIdsVe[ i_elVe[l_el][l_ve  ] ];

           // check if ascending, otherwise mission failed.
           EDGE_CHECK( l_gIdsVe[1] > l_gIdsVe[0] );
         }
       }

       /*
        * check ascending order of owned elements' faces.
        */
       for( int_tg l_tg = 0; l_tg < i_elLayout.timeGroups.size(); l_tg++ ) {
         int_el l_first = i_elLayout.timeGroups[l_tg].inner.first;
         int_el l_size = i_elLayout.timeGroups[l_tg].nEntsOwn;

         for( int_el l_el = l_first; l_el < l_first+l_size; l_el++ ) {
           // get the vertices of the adjacent faces
           int_el l_faVes[TL_N_EL_FAS][TL_N_FA_VES];

           for( unsigned int l_fa = 0; l_fa < TL_N_EL_FAS; l_fa++ ) {
             int_el l_faId = i_elFa[l_el][l_fa];
             for( unsigned int l_ve = 0; l_ve < TL_N_FA_VES; l_ve++ ) {
               l_faVes[l_fa][l_ve] = i_faVe[l_faId][l_ve];
             }
           }

           // check if the faces are ascending
           bool l_asc = true;
           for( unsigned int l_fa = 1; l_fa < TL_N_EL_FAS; l_fa++ ) {
             for( unsigned int l_ve = 0; l_ve < TL_N_FA_VES; l_ve++ ) {
               // get global ids of the vertices
               int_gid l_gIdsVe[2];
               l_gIdsVe[0] =  i_gIdsVe[ l_faVes[l_fa-1][l_ve] ];
               l_gIdsVe[1] =  i_gIdsVe[ l_faVes[l_fa  ][l_ve] ];

               // continue with next face if we are ascending
               if( l_gIdsVe[0] < l_gIdsVe[1] ) {
                 break;
               }
               // store for now if we are not ascending
               else if( l_gIdsVe[0] > l_gIdsVe[1] ) {
                 l_asc = false;
                 break;
               }
               // equal: continue searching
               else if( l_ve == TL_N_FA_VES - 1) {
                 EDGE_LOG_FATAL << "two faces of an element have identical vertices "
                                << l_el << " " << l_fa << " " << i_elFa[l_el][l_fa];
               }
             }
           }

           // are all the faces' vertices covered by the element's vertices?
           // -> we can identify an issue only if we are not on a periodic boundary
           bool l_covered = true;
           for( unsigned short l_fa = 0; l_fa < TL_N_EL_FAS; l_fa++ ) {
             for( unsigned short l_faVe = 0; l_faVe < TL_N_FA_VES; l_faVe++ ) {
               for( unsigned short l_elVe = 0; l_elVe < TL_N_EL_VES; l_elVe++ ) {
                 // continue with next vertex of the faces if its part of the elemnt
                 if( l_faVes[l_fa][l_faVe] == i_elVe[l_el][l_elVe] ) {
                   break;
                 }
                 // this face's vertex is not part of the element
                 else if( l_elVe == TL_N_EL_VES-1 ) {
                   l_covered = false;
                 }
               }
             }
           }

           if( i_periodic == false && l_covered == false ) {
             EDGE_LOG_FATAL << "not all vertices of the adjacent faces are covered by the element: " << l_el;
           }
           else if( l_asc == false && l_covered == true ) {
             EDGE_LOG_FATAL << "faces of the element not in ascending order, this should not happen: " << l_el;
           }
         }
       }
    }

    /**
     * Gets the coordinates of an elements vertices based on connectivity information.
     *
     * @param i_el id of the entity.
     * @param i_elVe entities' vertices.
     * @param i_veChars characteristics of the vertices.
     * @param o_veCrds will be set to coordinates of the vertices.
     *
     * @paramt TL_T_LID integral type of the local ids.
     * @paramt TL_T_REAL floating point type.
     * @paramt TL_T_VE_CHARS type of the vertex characteristics, offering member .coords.
     **/
    template< typename TL_T_LID,
              typename TL_T_REAL,
              typename TL_T_VE_CHARS >
    static void getElVeCrds(       TL_T_LID        i_el,
                             const TL_T_LID      (*i_elVe)[TL_N_EL_VES],
                             const TL_T_VE_CHARS  *i_veChars,
                                   TL_T_REAL       o_veCrds[TL_N_DIS][TL_N_EL_VES] ) {
      // get elements vertices
      for( unsigned short l_ve = 0; l_ve < TL_N_EL_VES; l_ve++ ) {
        TL_T_LID l_veId       = i_elVe[i_el][l_ve];

        for( unsigned short l_di = 0; l_di < TL_N_DIS; l_di++ ) {
          o_veCrds[l_di][l_ve] = i_veChars[l_veId].coords[l_di];
        }
      }
    }

    /**
     * Gets the local face ids from the adjacent elements' perspective.
     *
     * Remark: This covers adjacent, owned elements only.
     *
     * Example:
     *
     * Element L 2nd local face neighbor is element N. Element N neighbors L via face 1.
     * Therefore o_fIdEleFaEl[L][2] == 1.
     *       **********
     *     *   * 1 N *
     *    * L 2 *   *
     *   *       * *
     *  ***********
     *
     * @param i_elLayout data layout of the elements.
     * @param i_gIdsEl global ids of the elements.
     * @param i_elFa elements' adjacent faces.
     * @param i_elFaEl elements' adjacent elememnts connected through faces.
     * @param i_faceChars face characteristics.
     * @param o_fIdElFaEl will bet set to local face ids of face-neighboring elememts.
     **/
    static void getFIdsElFaEl( const t_enLayout             &i_elLayout,
                               const std::vector< int_gid > &i_gIdsEl,
                               const int_el                (*i_elFa)[TL_N_EL_FAS],
                               const int_el                (*i_elFaEl)[TL_N_EL_FAS],
                                     unsigned short        (*o_fIdElFaEl)[TL_N_EL_FAS] ) {
      // iterate over own entities
      for( int_tg l_tg = 0; l_tg < i_elLayout.timeGroups.size(); l_tg++ ) {
        int_el l_first = i_elLayout.timeGroups[l_tg].inner.first;
        int_el l_size  = i_elLayout.timeGroups[l_tg].nEntsOwn;

        for( int_el l_el = l_first; l_el < l_first+l_size; l_el++ ) {
          // global element id
          assert( l_el < (int_el) i_gIdsEl.size() );
          int_gid l_gIdEl = i_gIdsEl[l_el];

          for( unsigned int l_fa = 0; l_fa < TL_N_EL_FAS; l_fa++ ) {
            // assert an existing face
            EDGE_CHECK( i_elFa[l_el][l_fa] != std::numeric_limits<int_el>::max() );

            o_fIdElFaEl[l_el][l_fa] = std::numeric_limits<unsigned short>::max();

            // continue for boundary conditions
            if( i_elFaEl[l_el][l_fa] > i_elLayout.nEnts ) {
              o_fIdElFaEl[l_el][l_fa] = std::numeric_limits<unsigned short>::max();
              continue;
            }

            // local id of the neighbor
            int_el l_lIdElNe = i_elFaEl[l_el][l_fa];

            // assert an existing, adjacent element
            EDGE_CHECK( l_lIdElNe < (int_el) i_gIdsEl.size() );

            // find the local face mapping back to our element
            for( unsigned int l_nFa = 0; l_nFa < TL_N_EL_FAS; l_nFa++ ) {
              // continue for receive-elements pointing further out
              if( i_elFaEl[l_lIdElNe][l_nFa] == std::numeric_limits< int_el >::max() ) {
                continue;
              }

              EDGE_CHECK( i_elFaEl[l_lIdElNe][l_nFa] < (int_el) i_gIdsEl.size() );
              if( i_gIdsEl[ i_elFaEl[l_lIdElNe][l_nFa] ] == l_gIdEl ) {
                // this must be unique
                EDGE_CHECK( o_fIdElFaEl[l_el][l_fa] == std::numeric_limits<unsigned short>::max() );

                // let's store the info
                o_fIdElFaEl[l_el][l_fa] = l_nFa;
              }
            }

            // check that we found the corresponding face-id
            EDGE_CHECK( o_fIdElFaEl[l_el][l_fa] != std::numeric_limits<unsigned short>::max() );
          }
        }
      }
    }

    /**
     * Gets the adjancent elements' vertex information.
     *
     * @param i_elLayout data layout of the elements.
     * @param i_gIdsVe global ids of the vertices.
     * @param i_elFaEl elements' face neighboring elements.
     * @param i_elVe elements vertices.
     * @param i_fIdElFaEl local face ids of face-neighboring elememts.
     * @param o_vIdElFaEl will be set set to corresponding local vertex id w.r.t. the shared face from the neighboring elements' perspsective.
     * @param i_periodic if true we try to go through vertex coords for periodic faces.
     * @param i_veChars vertex characteristics. used when we assume periodic boundaries (having duplicated vertices).
     **/
    static void getVIdsElFaEl( const t_enLayout             &i_elLayout,
                               const std::vector< int_gid > &i_gIdsVe,
                               const int_el                (*i_elFaEl)[TL_N_EL_FAS],
                               const int_el                (*i_elVe)[TL_N_EL_VES],
                               const unsigned short        (*i_fIdElFaEl)[TL_N_EL_FAS],
                                     unsigned short        (*o_vIdElFaEl)[TL_N_EL_FAS],
                                     bool                    i_periodic=false,
                               const t_vertexChars          *i_veChars=NULL ) {
      // init
#ifdef PP_USE_OMP
#pragma omp parallel for
#endif
      for( int_el l_el = 0; l_el < i_elLayout.nEnts; l_el++ )
        for( unsigned short l_fa = 0; l_fa < TL_N_EL_FAS; l_fa++ )
          o_vIdElFaEl[l_el][l_fa] = std::numeric_limits< unsigned short >::max();

      for( int_tg l_tg = 0; l_tg < i_elLayout.timeGroups.size(); l_tg++ ) {
        int_el l_first = i_elLayout.timeGroups[l_tg].inner.first;
        int_el l_size  = i_elLayout.timeGroups[l_tg].nEntsOwn;

        for( int_el l_el = l_first; l_el < l_first+l_size; l_el++ ) {
#if defined PP_T_MESH_REGULAR || defined PP_T_ELEMENTS_TRIA3
          for( unsigned short l_fa = 0; l_fa < TL_N_EL_FAS; l_fa++ ) {
            // no rotation possible.
            o_vIdElFaEl[l_el][l_fa] = 0;
          }
#elif defined PP_T_ELEMENTS_TET4
          // iterate over element's faces
          for( unsigned short l_fa = 0; l_fa < TL_N_EL_FAS; l_fa++ ) {
            // neighboring element
            int_el l_ne = i_elFaEl[l_el][l_fa];

            // set invalid and continue for boundary conditions
            if( i_elFaEl[l_el][l_fa] == std::numeric_limits<int_el>::max() ) {
              o_vIdElFaEl[l_el][l_fa] = std::numeric_limits<unsigned short>::max();
              continue;
            }

            // check if we are at a periodic boundary
            unsigned int l_nSharedVes = 0;
            for( unsigned l_ve1 = 0; l_ve1 < 4; l_ve1++ ) {
              for( unsigned int l_ve2 = 0; l_ve2 < 4; l_ve2++ ) {
                if( i_elVe[l_el][l_ve1] == i_elVe[l_ne][l_ve2] ) l_nSharedVes++;
              }
            }
            // this might fail for very small meshes, where all tets are connected to each other (l_nSharedVes)
            // we keep the check as this is unrealistic
            EDGE_CHECK( l_nSharedVes == 3 || l_nSharedVes == 0 );

            // global ids of neighboring vertices
            int_gid l_neVe[4];

            // default, we found a neighbor
            if( l_nSharedVes == 3 ) {
              for( unsigned short l_ve = 0; l_ve < 4; l_ve++ ) l_neVe[l_ve] = i_gIdsVe[ i_elVe[l_ne][l_ve] ];
            }
            // either abort or try to save things for periodic boundaries
            else {
              if( i_periodic == false ) EDGE_LOG_FATAL << "couldn't derive vertex ids, assuming a non-periodic mesh";

              // go through geometric properties and set dummy values to the 'neighboring' vertices
              setPeriodicVesTet4( i_gIdsVe, l_el, l_ne, l_fa, i_elVe, i_veChars, l_neVe );
            }

            // dominant vertex, which dictates the vertex id
            int_gid l_domVe;

            // in the case of faces 0-2 we are interested in the position of local vertex 0
            if( l_fa <= 2 ) l_domVe = i_gIdsVe[ i_elVe[l_el][0] ];
            // face 3 requires the position of local vertex 1 w.r.t. the neighboring face
            else            l_domVe = i_gIdsVe[ i_elVe[l_el][1] ];

            // find the position of the vertex in the other element
            unsigned short l_locNe = 4;

            for( unsigned int l_ve2 = 0; l_ve2 < 4; l_ve2++ ) {
              if( l_neVe[l_ve2] == l_domVe ) {
                EDGE_CHECK_EQ( l_locNe, 4 );
                l_locNe = l_ve2;
              }
            }
            assert( l_locNe != 4 );

            /*
             * decide, depending on the neighboring face, what vertex combination we have
             *
             * Example:
             *
             *      face 0                   face 1
             *         3                       2
             *         *                         *
             *       *    *         neighbors    *  *
             *     *        *                    *     *
             * 0  ************* 1              0 ********** 3 <-- dominant 0 goes here
             *
             * -> We have vertex combi 1 out of possible 0-2.
             */
             if( i_fIdElFaEl[l_el][l_fa] == 0 ) {
               if(      l_locNe == 0 ) o_vIdElFaEl[l_el][l_fa] = 0;
               else if( l_locNe == 2 ) o_vIdElFaEl[l_el][l_fa] = 1;
               else if( l_locNe == 1 ) o_vIdElFaEl[l_el][l_fa] = 2;
               else assert( false );
             }
             else if( i_fIdElFaEl[l_el][l_fa] == 1 ) {
               if(      l_locNe == 0 ) o_vIdElFaEl[l_el][l_fa] = 0;
               else if( l_locNe == 1 ) o_vIdElFaEl[l_el][l_fa] = 1;
               else if( l_locNe == 3 ) o_vIdElFaEl[l_el][l_fa] = 2;
               else assert( false );
             }
             else if( i_fIdElFaEl[l_el][l_fa] == 2 ) {
               if(      l_locNe == 0 ) o_vIdElFaEl[l_el][l_fa] = 0;
               else if( l_locNe == 3 ) o_vIdElFaEl[l_el][l_fa] = 1;
               else if( l_locNe == 2 ) o_vIdElFaEl[l_el][l_fa] = 2;
               else assert( false );
             }
             else if( i_fIdElFaEl[l_el][l_fa] == 3 ) {
               if(      l_locNe == 1 ) o_vIdElFaEl[l_el][l_fa] = 0;
               else if( l_locNe == 2 ) o_vIdElFaEl[l_el][l_fa] = 1;
               else if( l_locNe == 3 ) o_vIdElFaEl[l_el][l_fa] = 2;
               else assert( false );
             }
             else assert( false );
          }
#else
         assert( false );
#endif
        }
      }
    }

    /**
     * Prints a summary of the present neighboring relations.
     *
     * The index is given by:
     *     local_face_id-neighboring_face_id-neighboring_vertex_id
     *
     * @param i_nElements number of elements.
     * @param i_faIdElFaEl face relations.
     * @param i_veIdElFaEl vertex relations.
     **/
    static void printNeighRel(  const t_enLayout     &i_elLayout,
                                const unsigned short *i_fIdElFaEl,
                                const unsigned short *i_vIdElFaEl ) {
      t_entityType l_elType = TL_T_EL;
      t_entityType l_faType = TL_T_FA;

      // jump for the neighboring face (1 if not 3D)
      unsigned short l_veJump = 1;
      if( N_DIM == 3 ) {
        l_veJump = C_ENT[l_faType].N_VERTICES;
      }

      // compute number of options
      unsigned short l_nOpts = C_ENT[l_elType].N_FACES*C_ENT[l_elType].N_FACES*l_veJump;

      // reset options
      unsigned long *l_optsLo = new unsigned long[l_nOpts];
      for( unsigned short l_opt = 0; l_opt < l_nOpts; l_opt++ ) l_optsLo[l_opt] = 0;

      for( int_tg l_tg = 0; l_tg < i_elLayout.timeGroups.size(); l_tg++ ) {
        int_el l_first = i_elLayout.timeGroups[l_tg].inner.first;
        int_el l_size  = i_elLayout.timeGroups[l_tg].nEntsOwn;

        for( int_el l_el = l_first; l_el < l_first+l_size; l_el++ ) {
          for( unsigned short l_fa = 0; l_fa < C_ENT[l_elType].N_FACES; l_fa++ ) {
            if( i_fIdElFaEl[l_el*C_ENT[l_elType].N_FACES+l_fa] < C_ENT[l_elType].N_FACES ) {
              unsigned short l_locOpt = l_fa * C_ENT[l_elType].N_FACES*l_veJump
                                        +
                                        i_fIdElFaEl[l_el*C_ENT[l_elType].N_FACES+l_fa] * l_veJump
                                        +
                                        i_vIdElFaEl[l_el*C_ENT[l_elType].N_FACES+l_fa];
              l_optsLo[l_locOpt]++;
            }
          }
        }
      }

      unsigned long *l_optsGl;
#ifdef PP_USE_MPI
      l_optsGl = new unsigned long[l_nOpts];
      // compute sum over all ranks
      MPI_Allreduce( l_optsLo, l_optsGl, l_nOpts, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD );
#else
      l_optsGl = l_optsLo;
#endif

      for( unsigned short l_opt = 0; l_opt < l_nOpts; l_opt++ ) {
        EDGE_LOG_INFO << "    id " << l_opt / (C_ENT[l_elType].N_FACES*l_veJump)
                      << "-" << ( l_opt%(C_ENT[l_elType].N_FACES*l_veJump) ) / l_veJump
                      << "-" << l_opt%l_veJump
                      << ": " << l_optsGl[l_opt];
      }

      delete[] l_optsLo;
#ifdef PP_USE_MPI
      delete[] l_optsGl;
#endif
    }
};
#endif
