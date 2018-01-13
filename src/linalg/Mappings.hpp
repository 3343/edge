/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2016-2017, Regents of the University of California
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
 * Coordinate transformations and mappings.
 **/

#ifndef MAPPINGS_HPP
#define MAPPINGS_HPP

#include <cmath>
#include "io/logging.h"
#include "constants.hpp"
#include "Shape.hpp"

namespace edge {
  namespace linalg {
    class Mappings;
  }
}

class edge::linalg::Mappings {
  public:
    /**
     * Computes the distance between two points.
     *
     * @param i_pts coords of points.
     * @return distance.
     **/
    static real_mesh distPts( const real_mesh i_pts[3][2] ) {
      real_mesh l_dist  = ( i_pts[0][1] - i_pts[0][0] ) * ( i_pts[0][1] - i_pts[0][0] );
                l_dist += ( i_pts[1][1] - i_pts[1][0] ) * ( i_pts[1][1] - i_pts[1][0] );
                l_dist += ( i_pts[2][1] - i_pts[2][0] ) * ( i_pts[2][1] - i_pts[2][0] );
                l_dist  = std::sqrt( l_dist );

      return l_dist;
    }

    /**
     * Derives the coordinates of a point in a line element w.r.t. to another element.
     *
     * 2D example (derives coords of out from in):
     *
     *               x from[*][1]
     *              *
     *      1st    o in[*]
     *    element *             x to[*][1]
     *  /|\      *               *            2nd
     *   |      *                 *          element
     * y |     x from[*][0]        o out[*]
     *   |_____\           to[*][1] x
     *      x  /
     *
     * @param i_veCoordsFrom coordinates of 1st elements vertices.
     * @param i_veCoordsTo coordinates of 2nd elements vertices.
     * @param i_coordsIn coordinates of 1st elements point, which is mapped to the second element.
     * @param o_coordsOut will be set to coordinates of input point w.r.t. to the second element.
     **/
    static void map2dLine( const real_mesh i_veCoordsFrom[3][2],
                           const real_mesh i_veCoordsTo[3][2],
                           const real_mesh i_coordsIn[3],
                                 real_mesh o_coordsOut[3] ) {
      // compute distance of 1st element's first vertex and input
      real_mesh l_distIn  = ( i_coordsIn[0] - i_veCoordsFrom[0][0] ) * ( i_coordsIn[0] - i_veCoordsFrom[0][0] );
                l_distIn += ( i_coordsIn[1] - i_veCoordsFrom[1][0] ) * ( i_coordsIn[1] - i_veCoordsFrom[1][0] );
                l_distIn += ( i_coordsIn[2] - i_veCoordsFrom[2][0] ) * ( i_coordsIn[2] - i_veCoordsFrom[2][0] );
                l_distIn  = std::sqrt( l_distIn );

      // compute distance of 1st element's first and second vertex
      real_mesh l_distVe1  = distPts( i_veCoordsFrom );

      // derive output by starting at 2nd element's first vertex
      o_coordsOut[0] = i_veCoordsTo[0][0];
      o_coordsOut[1] = i_veCoordsTo[1][0];
      o_coordsOut[2] = i_veCoordsTo[2][0];

      // move towards the 2nd vertex by the given length
      real_mesh l_scale = l_distIn / l_distVe1;
      o_coordsOut[0] += l_scale * ( i_veCoordsTo[0][1] - i_veCoordsTo[0][0] );
      o_coordsOut[1] += l_scale * ( i_veCoordsTo[1][1] - i_veCoordsTo[1][0] );
      o_coordsOut[2] += l_scale * ( i_veCoordsTo[2][1] - i_veCoordsTo[2][0] );
    }

   /*
    * Maps the point of one quad to the corresponding point of another quad in 3d.
    * This implementation is limited to the reference-faces of HEX8R elements.
    *
    * The unqiue face neighbor of HEX8R elements have identical, local coordinates systems:
    *   Example:
    *     face 1 (0-1-5-4, eta=0) is mapped to face 3 (3-2-6-7, eta=1).
    *     This means that both have a local counter-clockwise storage.
    *     Face coords of face 1 are simply shifted by eta=1 to match face 3.
    *     Face coords of face 3 are simply shifted by eta=-1 to match face 1.
    *
    *    4               5                        7               6
    *      *************                            *************
    *      *           *                            *           *
    *      |           *            ---->           |           *
    * zeta |           *                       zeta |           *
    *      |           *                            |           *
    *      x---*********                            x---*********
    *    0   xi          1                        3   xi          2
    *
    * @param i_fa1 id of the local face.
    * @param i_fa2 id of the neighboring face.
    * @param i_coordsIn reference coordinates of the first, local quad.
    * @param o_coordsOut will be set to coords w.r.t. to the second, neighboring quad.
    *
    */
    static void map3dQuad4(        unsigned int i_fa1,
                                   unsigned int i_fa2,
                             const real_mesh i_coordsIn[3],
                                   real_mesh o_coordsOut[3] ) {
      // initialize output
      o_coordsOut[0] = i_coordsIn[0];
      o_coordsOut[1] = i_coordsIn[1];
      o_coordsOut[2] = i_coordsIn[2];

      if(       i_fa1 == 0 ) {
        EDGE_CHECK_EQ( i_fa2, 5 );
        o_coordsOut[2] += 1;
      }
      else if(  i_fa1 == 1 ) {
        EDGE_CHECK_EQ( i_fa2, 3 );
        o_coordsOut[1] += 1;
      }
      else if(  i_fa1 == 2 ) {
        EDGE_CHECK_EQ( i_fa2, 4 );
        o_coordsOut[0] -= 1;
      }
      else if(  i_fa1 == 3 ) {
        EDGE_CHECK_EQ( i_fa2, 1 );
        o_coordsOut[1] -= 1;
      }
      else if(  i_fa1 == 4 ) {
        EDGE_CHECK_EQ( i_fa2, 2 );
        o_coordsOut[0] += 1;
      }
      else if(  i_fa1 == 5 ) {
        EDGE_CHECK_EQ( i_fa2, 0 );
        o_coordsOut[2] -= 1;
      }
      else EDGE_LOG_FATAL;

      // make sure we are still in bounds
      EDGE_CHECK_GT( o_coordsOut[0] + TOL.MESH, 0 );
      EDGE_CHECK_GT( o_coordsOut[1] + TOL.MESH, 0 );
      EDGE_CHECK_GT( o_coordsOut[2] + TOL.MESH, 0 );
      EDGE_CHECK_LT( o_coordsOut[0] - TOL.MESH, 1 );
      EDGE_CHECK_LT( o_coordsOut[1] - TOL.MESH, 1 );
      EDGE_CHECK_LT( o_coordsOut[2] - TOL.MESH, 1 );
    }

    /**
     * Evaluates the jacobian for the given element type with the given coordinates at the given location (optional).
     *
     * @param i_entType entity type.
     * @param i_veCoords coordinates of the entity.
     * @param o_jac will be se to evaluated jacobian.
     * @param i_xi xi position at which the jacobian is evaluated; only relevant for non-constant jacobians.
     * @param i_eta eta position at which the jacobian is evaluated; only relevant for non-constant jacobians.
     * @param i_zeta zeta position at which the jacobian is evaluated; only relevant for non-constant jacobians.
     **/
    static void evalJac(       t_entityType  i_entType,
                         const real_mesh    *i_veCoords,
                               real_mesh    *o_jac,
                               real_mesh     i_xi   = 0,
                               real_mesh     i_eta  = 0,
                               real_mesh     i_zeta = 0 ) {
      // derivatives of the shape functions, we use 3 dims and 8 nodes as current upper bound
      real_mesh l_evalSfDer[3*8];
      for( unsigned short l_en = 0; l_en < 3*8; l_en++ ) l_evalSfDer[l_en] = std::numeric_limits< real_mesh >::max();
      EDGE_CHECK_LE( C_ENT[i_entType].N_VERTICES, 8 );

      if( i_entType == LINE ) {
        Shape::derLine( l_evalSfDer );
      }
      else if( i_entType == QUAD4R ) {
        Shape::derQuad4( 0, 0, (real_mesh (*)[4]) l_evalSfDer );
      }
      else if( i_entType == TRIA3 ) {
        Shape::derTria3( (real_mesh (*)[3]) l_evalSfDer );
      }
      else if( i_entType == HEX8R ) {
        Shape::derHex8Lin( 0, 0, 0, (real_mesh (*)[8]) l_evalSfDer );
      }
      else if( i_entType == TET4 ) {
        Shape::derTet4( (real_mesh(*)[4]) l_evalSfDer );
      }
      else EDGE_LOG_FATAL;

      // get #dims and vertices
      unsigned short l_nDims = C_ENT[i_entType].N_DIM;
      unsigned short l_nVes  = C_ENT[i_entType].N_VERTICES;

      // compute the jacobian
      for( unsigned int l_dim1 = 0; l_dim1 < l_nDims; l_dim1++ ) {
        for( unsigned int l_dim2 = 0; l_dim2 < l_nDims; l_dim2++ ) {
          o_jac[l_dim1*l_nDims + l_dim2] = 0;
          for( unsigned int l_pt = 0; l_pt < l_nVes; l_pt++ ) {
            o_jac[l_dim1*l_nDims + l_dim2] += i_veCoords[l_dim2*l_nVes+l_pt] * l_evalSfDer[l_dim1*l_nVes+l_pt];
          }
        }
      }
    }

    /**
     * Maps the point in reference coordinates xi to physical coords
     * w.r.t. to the vertices of the line element.
     *
     * @param i_xi coord in xi-direction.
     * @param i_veCoords coords of the line element.
     * @param o_pt will be set to resulting point in physical coordinates.
     *
     * @paramt TL_T_REAL floating point precision.
     **/
    template< typename TL_T_REAL >
    static void mapRefToPhyLine( TL_T_REAL       i_xi,
                                 TL_T_REAL const i_veCoords[1][2],
                                 TL_T_REAL       o_pt[1] ) {
      o_pt[0]  = i_veCoords[0][1] - i_veCoords[0][0];
      o_pt[0] *= i_xi;
      o_pt[0] += i_veCoords[0][0];
    }

    /**
     * Maps the point in reference coordinates xi and eta to physical coords
     * w.r.t. to the vertices of the quadrilateral.
     *
     * @param i_xi coord in xi-direction.
     * @param i_eta coord in eta-direction.
     * @param i_veCoords coords of the quadrilateral.
     * @param o_pt will be set to resulting point in physical coordinates.
     *
     * @paramt TL_T_REAL floating point precision.
     **/
    template< typename TL_T_REAL >
    static void mapRefToPhyQuad4( TL_T_REAL       i_xi,
                                  TL_T_REAL       i_eta,
                                  TL_T_REAL const i_veCoords[2][4],
                                  TL_T_REAL       o_pt[2] ) {
      // get the evaluated shape functions
      TL_T_REAL l_evalSf[4];
      Shape::quad4( i_xi, i_eta, l_evalSf );

      // perform the mapping
      for( unsigned int l_dim = 0; l_dim < 2; l_dim++ ) {
        o_pt[l_dim] = 0;
        for( unsigned int l_pt = 0; l_pt < 4; l_pt++ ) {
          o_pt[l_dim] += i_veCoords[l_dim][l_pt] * l_evalSf[l_pt];
        }
      }
    }

    /**
     * Maps the point in reference coordinates xi and eta to physical coords
     * w.r.t. to the vertices of the triangle.
     *
     * @param i_xi coord in xi-direction.
     * @param i_eta coord in eta-direction.
     * @param i_veCoords coords of the triangle.
     * @param o_pt will be set to resulting point in physical coordinates.
     *
     * @paramt TL_T_REAL floating point precision.
     **/
    template< typename TL_T_REAL >
    static void mapRefToPhyTria3( TL_T_REAL       i_xi,
                                  TL_T_REAL       i_eta,
                                  TL_T_REAL const i_veCoords[2][3],
                                  TL_T_REAL       o_pt[2] ) {
      // get the evaluated shape functions
      TL_T_REAL l_evalSf[3];
      Shape::tria3( i_xi, i_eta, l_evalSf );

      // perform the mapping
      for( unsigned int l_dim = 0; l_dim < 2; l_dim++ ) {
        o_pt[l_dim] = 0;
        for( unsigned int l_pt = 0; l_pt < 3; l_pt++ ) {
          o_pt[l_dim] += i_veCoords[l_dim][l_pt] * l_evalSf[l_pt];
        }
      }
    }

    /**
     * Maps a point in reference coordinates xi, eta and zeta of a linear 8-node hexahedral element
     * to physical coords w.r.t. to the vertices of the hex.
     *
     * @param i_xi coord in xi-direction.
     * @param i_eta coord in eta-direction.
     * @param i_zeta coord in zeta-direction
     * @param i_veCoords coords of the hexahedron.
     * @param o_pt will be set to resulting point in physical coordinates.
     *
     * @paramt TL_T_REAL floating point precision.
     **/
    template< typename TL_T_REAL >
    static void mapRefToPhyHex8r( TL_T_REAL       i_xi,
                                  TL_T_REAL       i_eta,
                                  TL_T_REAL       i_zeta,
                                  TL_T_REAL const i_veCoords[3][8],
                                  TL_T_REAL       o_pt[3] ) {
      // get the evaluated shape functions
      TL_T_REAL l_evalSf[8];
      Shape::hex8Lin( i_xi, i_eta, i_zeta, l_evalSf );

      // perform the mapping
      for( unsigned int l_dim = 0; l_dim < 3; l_dim++ ) {
        o_pt[l_dim] = 0;
        for( unsigned int l_pt = 0; l_pt < 8; l_pt++ ) {
          o_pt[l_dim] += i_veCoords[l_dim][l_pt] * l_evalSf[l_pt];
        }
      }
    }

    /**
     * Maps the point in reference coordinates xi, eta and zeta to physical coords
     * w.r.t. to the vertices of the tetrahedron.
     *
     * @param i_xi coord in xi-direction.
     * @param i_eta coord in eta-direction.
     * @param i_zeta coord in zeta-direction
     * @param i_veCoords coords of the tetrahedron.
     * @param o_pt will be set to resulting point in physical coordinates.
     *
     * @paramt TL_T_REAL floating point precision.
     **/
    template< typename TL_T_REAL >
    static void mapRefToPhyTet4( TL_T_REAL       i_xi,
                                 TL_T_REAL       i_eta,
                                 TL_T_REAL       i_zeta,
                                 TL_T_REAL const i_veCoords[3][4],
                                 TL_T_REAL       o_pt[3] ) {
      // get the evaluated shape functions
      TL_T_REAL l_evalSf[4];
      Shape::tet4( i_xi, i_eta, i_zeta, l_evalSf );

      // perform the mapping
      for( unsigned int l_dim = 0; l_dim < 3; l_dim++ ) {
        o_pt[l_dim] = 0;
        for( unsigned int l_pt = 0; l_pt < 4; l_pt++ ) {
          o_pt[l_dim] += i_veCoords[l_dim][l_pt] * l_evalSf[l_pt];
        }
      }
    }

    /**
     * Maps the given point in reference coordinates to physical coordinates.
     *
     * @param i_enType entity type.
     * @param i_ves vertices of the entity in physical coordinates (vertices are the fastest dimension).
     * @param i_xi point in reference coords.
     * @param o_x will be set point in physical coords.
     *
     * @paramt TL_T_REAL floating point precision.
     **/
    template< typename TL_T_REAL >
    static void refToPhy( t_entityType        i_enType,
                          TL_T_REAL    const *i_ves,
                          TL_T_REAL    const *i_xi,
                          TL_T_REAL          *o_x ) {
     if(      i_enType == LINE   ) mapRefToPhyLine( i_xi[0],
                                                    (TL_T_REAL const (*)[2]) i_ves,
                                                    o_x );
     else if( i_enType == QUAD4R ) mapRefToPhyQuad4( i_xi[0], i_xi[1],
                                                     (TL_T_REAL const (*)[4]) i_ves,
                                                     o_x );
     else if( i_enType == TRIA3  ) mapRefToPhyTria3( i_xi[0], i_xi[1],
                                                     (TL_T_REAL const (*)[3]) i_ves,
                                                     o_x );
     else if( i_enType == HEX8R  ) mapRefToPhyHex8r( i_xi[0], i_xi[1], i_xi[2],
                                                     (TL_T_REAL const (*)[8]) i_ves,
                                                     o_x );
     else if( i_enType == TET4   ) mapRefToPhyTet4( i_xi[0], i_xi[1], i_xi[2],
                                                    (TL_T_REAL (*)[4]) i_ves,
                                                    o_x  );
     else                          EDGE_LOG_FATAL;
    }

    /**
     * Maps the given point of a tet4-element in physical coordinates to reference coordinates.
     *
     * Reference: Advanced Finite Element Methods (ASEN 6367) - Spring 2013 
     *            Department of Aerospace Engineering Sciences University of Colorado at Boulder
     *            Chapter 9: The Linear Tetrahedron
     *
     * @param i_ves vertices of the tet in physical coordinates ([*][]: dim, [][*]: vertex).
     * @param i_x point in physical coords.
     * @parma o_xi will be set to point in reference coordinates.
     **/
    template <typename T>
    static void phyToRefTet4( const T i_ves[3][4],
                              const T i_x[3],
                                    T o_xi[3] ) {
      // determine inverse matrix of the trafo from ref-to-physical coords

      // generic inversion might be problematic, we go through explicit formulas
#if 0
      const T *l_ves = i_ves[0];
      T l_inv[16];
#include "generated/ptr_inv_4x4.inc"
#endif

      T l_m[4][4];
      //    x2 (y3 z4 - y4 z3) + x3 (y4 z2 - y2 z4) + x4 (y2 z3 − y3 z2)
      // -> 01  12 23 - 13 22  + 02  13 21 - 11 23  + 03  11 22 - 12 21
      l_m[0][0] =   i_ves[0][1]*(i_ves[1][2]*i_ves[2][3] - i_ves[1][3]*i_ves[2][2])
                  + i_ves[0][2]*(i_ves[1][3]*i_ves[2][1] - i_ves[1][1]*i_ves[2][3])
                  + i_ves[0][3]*(i_ves[1][1]*i_ves[2][2] - i_ves[1][2]*i_ves[2][1]);

      //    x1 (y4 z3 - y3 z4) + x3 (y1 z4 - y4 z1) + x4 (y3 z1 - y1 z3)
      // -> 00  13 22   12 23    02  10 23   13 20    03  12 20   10 22
      l_m[1][0] =   i_ves[0][0]*(i_ves[1][3]*i_ves[2][2] - i_ves[1][2]*i_ves[2][3])
                  + i_ves[0][2]*(i_ves[1][0]*i_ves[2][3] - i_ves[1][3]*i_ves[2][0])
                  + i_ves[0][3]*(i_ves[1][2]*i_ves[2][0] - i_ves[1][0]*i_ves[2][2]);

      //    x1 (y2 z4 - y4 z2) + x2 (y4 z1 - y1 z4) + x4 (y1 z2 - y2 z1)
      // -> 00  11 23   13 21    01  13 20   10 23    03  10 21   11 20
      l_m[2][0] =   i_ves[0][0]*(i_ves[1][1]*i_ves[2][3] - i_ves[1][3]*i_ves[2][1])
                  + i_ves[0][1]*(i_ves[1][3]*i_ves[2][0] - i_ves[1][0]*i_ves[2][3])
                  + i_ves[0][3]*(i_ves[1][0]*i_ves[2][1] - i_ves[1][1]*i_ves[2][0]);

      //    x1 (y3 z2 - y2 z3) + x2 (y1 z3 - y3 z1) + x3 (y2 z1 - y1 z2)
      // -> 00  12 21   11 22    01  10 22   12 20    02  11 20   10 21
      l_m[3][0] =   i_ves[0][0]*(i_ves[1][2]*i_ves[2][1] - i_ves[1][1]*i_ves[2][2])
                  + i_ves[0][1]*(i_ves[1][0]*i_ves[2][2] - i_ves[1][2]*i_ves[2][0])
                  + i_ves[0][2]*(i_ves[1][1]*i_ves[2][0] - i_ves[1][0]*i_ves[2][1]);

      //    y 42 z 32 - y 32 z 42
      // -> 1 31 2 21   1 21 2 31
      l_m[0][1] =  (i_ves[1][3]-i_ves[1][1]) * (i_ves[2][2]-i_ves[2][1])
                  -(i_ves[1][2]-i_ves[1][1]) * (i_ves[2][3]-i_ves[2][1]);
      //    x 32 z 42 - x 42 z 32
      // -> 0 21 2 31   0 31 2 21
      l_m[0][2] =  (i_ves[0][2]-i_ves[0][1]) * (i_ves[2][3]-i_ves[2][1])
                  -(i_ves[0][3]-i_ves[0][1]) * (i_ves[2][2]-i_ves[2][1]);
      //    x 42 y 32 - x 32 y 42
      // -> 0 31 1 21   0 21 1 31
      l_m[0][3] =  (i_ves[0][3]-i_ves[0][1]) * (i_ves[1][2]-i_ves[1][1])
                  -(i_ves[0][2]-i_ves[0][1]) * (i_ves[1][3]-i_ves[1][1]);

      //    y 31 z 43 - y 34 z 13
      // -> 1 20 2 32   1 23 2 02
      l_m[1][1] =  (i_ves[1][2]-i_ves[1][0]) * (i_ves[2][3]-i_ves[2][2])
                  -(i_ves[1][2]-i_ves[1][3]) * (i_ves[2][0]-i_ves[2][2]);
      //    x 43 z 31 - x 13 z 34
      // -> 0 32 2 20   0 02 2 23
      l_m[1][2] =  (i_ves[0][3]-i_ves[0][2]) * (i_ves[2][2]-i_ves[2][0])
                  -(i_ves[0][0]-i_ves[0][2]) * (i_ves[2][2]-i_ves[2][3]);
      //    x 31 y 43 - x 34 y 13
      // -> 0 20 1 32   0 23 1 02
      l_m[1][3] =  (i_ves[0][2]-i_ves[0][0]) * (i_ves[1][3]-i_ves[1][2])
                  -(i_ves[0][2]-i_ves[0][3]) * (i_ves[1][0]-i_ves[1][2]);

      //    y 24 z 14 - y 14 z 24
      // -> 1 13 2 03   1 03 2 13
      l_m[2][1] =  (i_ves[1][1]-i_ves[1][3]) * (i_ves[2][0]-i_ves[2][3])
                  -(i_ves[1][0]-i_ves[1][3]) * (i_ves[2][1]-i_ves[2][3]);
      //    x 14 z 24 - x 24 z 14
      // -> 0 03 2 13   0 13 2 03
      l_m[2][2] =  (i_ves[0][0]-i_ves[0][3]) * (i_ves[2][1]-i_ves[2][3])
                  -(i_ves[0][1]-i_ves[0][3]) * (i_ves[2][0]-i_ves[2][3]);
      //    x 24 y 14 - x 14 y 24
      // -> 0 13 1 03   0 03 1 13
      l_m[2][3] =  (i_ves[0][1]-i_ves[0][3]) * (i_ves[1][0]-i_ves[1][3])
                  -(i_ves[0][0]-i_ves[0][3]) * (i_ves[1][1]-i_ves[1][3]);

      //    y 13 z 21 - y 12 z 31
      // -> 1 02 2 10   1 01 2 20
      l_m[3][1] =  (i_ves[1][0]-i_ves[1][2]) * (i_ves[2][1]-i_ves[2][0])
                  -(i_ves[1][0]-i_ves[1][1]) * (i_ves[2][2]-i_ves[2][0]);
      //    x 21 z 13 - x 31 z 12
      // -> 0 10 2 02   0 20 2 01
      l_m[3][2] =  (i_ves[0][1]-i_ves[0][0]) * (i_ves[2][0]-i_ves[2][2])
                  -(i_ves[0][2]-i_ves[0][0]) * (i_ves[2][0]-i_ves[2][1]);
      //    x 13 y 21 - x 12 y 31
      // -> 0 02 1 10   0 01 1 20
      l_m[3][3] =  (i_ves[0][0]-i_ves[0][2]) * (i_ves[1][1]-i_ves[1][0])
                  -(i_ves[0][0]-i_ves[0][1]) * (i_ves[1][2]-i_ves[1][0]);

      // divide by 6x the volume
      T l_div = l_m[0][0] + l_m[1][0] + l_m[2][0] + l_m[3][0];
      l_div = T(1) / l_div;

      for( unsigned short l_ro = 0; l_ro < 4; l_ro++ )
        for( unsigned short l_co = 0; l_co < 4; l_co++ )
          l_m[l_ro][l_co] *= l_div;

      T *l_inv = l_m[0];

      // derive the evaluated shape functions 1-3 (not 0) which directly translate to the reference coordinates
      for( unsigned short l_ro = 1; l_ro < 4; l_ro++ ) {
        o_xi[l_ro-1] = l_inv[l_ro*4+0];
        for( unsigned short l_ve = 1; l_ve < 4; l_ve++ )
          o_xi[l_ro-1] += l_inv[l_ro*4+l_ve] * i_x[l_ve-1];
      }
    }

    /**
     * Maps the given point of a tria3-element in physical coordinates to reference coordinates.
     *
     * Reference: An arbitrary high-order discontinuous Galerkin method for elastic
     *            waves on unstructured meshes – I. The two-dimensional isotropic
     *            case with external source terms
     *            Kaeser, Dumbser
     *
     * @param i_ves vertices of the triangle in physical coordinates ([*][]: dim, [][*]: vertex).
     * @param i_x point in physical coords.
     * @parma o_xi will be set to point in reference coordinates.
     **/
    template <typename T>
    static void phyToRefTria3( const T i_ves[2][3],
                               const T i_x[2],
                                     T o_xi[2] ) {
      //    x 3 y 1 - x 1 y 3
      // -> 0 2 1 0 - 0 0 1 2
      o_xi[0]  = i_ves[0][2] * i_ves[1][0] - i_ves[0][0] * i_ves[1][2];
      //    y 3 - y 1 + x 1 - x 3
      // -> 1 2 - 1 0 + 0 0 - 0 2
      o_xi[0] += i_x[0] * ( i_ves[1][2] - i_ves[1][0] ) + i_x[1] * ( i_ves[0][0] - i_ves[0][2] );

      //    x 1 y 2 - x 2 y 1
      //    0 0 1 1 - 0 1 1 0
      o_xi[1]  = i_ves[0][0] * i_ves[1][1] - i_ves[0][1] * i_ves[1][0];
      //    y 1 - y 2 + x 2 - x 1
      //    1 0 - 1 1 + 0 1 - 0 0
      o_xi[1] += i_x[0] * ( i_ves[1][0] - i_ves[1][1] ) + i_x[1] * ( i_ves[0][1] - i_ves[0][0] );

      //    x 2 - x 1  y 3 - y 1
      // -> 0 1 - 0 0  1 2 - 1 0
      T l_det  = (i_ves[0][1] - i_ves[0][0]) * (i_ves[1][2] - i_ves[1][0]);
      //    x 3 - x 1  y 2 - y 1
      // -> 0 2 - 0 0  1 1 - 1 0
        l_det -= (i_ves[0][2] - i_ves[0][0]) * (i_ves[1][1] - i_ves[1][0]);

      EDGE_CHECK_GT( l_det, TOL.MESH );
      o_xi[0] /= l_det;
      o_xi[1] /= l_det;
    }

    /**
     * Maps the given point of a quad4r-element in physical coordinates to reference coordinates.
     *
     * @param i_ves vertices of the quad in physical coordinates ([*][]: dim, [][*]: vertex).
     * @param i_x point in physical coords.
     * @parma o_xi will be set to point in reference coordinates.
     **/
    template <typename T>
    static void phyToRefQuad4r( const T i_ves[2][4],
                                const T i_x[2],
                                      T o_xi[2] ) {
      T l_dx = i_ves[0][1] - i_ves[0][0];
      EDGE_CHECK_GT( l_dx, TOL.MESH );
      T l_dy = i_ves[1][3] - i_ves[1][0];
      EDGE_CHECK_GT( l_dy, TOL.MESH );

      // shift
      o_xi[0] = i_x[0] - i_ves[0][0];
      o_xi[1] = i_x[1] - i_ves[1][0];
      // scale
      o_xi[0] /= l_dx;
      o_xi[1] /= l_dy;
    }

    /**
     * Maps the given point of a hex8r-element in physical coordinates to reference coordinates.
     *
     * @param i_ves vertices of the hex8r element in physical coordinates ([*][]: dim, [][*]: vertex).
     * @param i_x point in physical coords.
     * @param o_xi will be set to point in reference coordinates.
     *
     * @paramt TL_T_REAL floating point precision.
     **/
    template <typename TL_T_REAL>
    static void phyToRefHex8r( TL_T_REAL const i_ves[3][8],
                               TL_T_REAL const i_x[3],
                               TL_T_REAL       o_xi[3] ) {
      TL_T_REAL l_dx = i_ves[0][1] - i_ves[0][0];
      EDGE_CHECK_GT( l_dx, TOL.MESH );
      TL_T_REAL l_dy = i_ves[1][3] - i_ves[1][0];
      EDGE_CHECK_GT( l_dy, TOL.MESH );
      TL_T_REAL l_dz = i_ves[2][4] - i_ves[2][0];
      EDGE_CHECK_GT( l_dz, TOL.MESH );

      // shift
      o_xi[0] = i_x[0] - i_ves[0][0];
      o_xi[1] = i_x[1] - i_ves[1][0];
      o_xi[2] = i_x[2] - i_ves[2][0];
      // scale
      o_xi[0] /= l_dx;
      o_xi[1] /= l_dy;
      o_xi[2] /= l_dz;
    }

    /**
     * Maps the given point in physical coordinates to reference coordinates.
     *
     * @param i_enType entity type.
     * @param i_ves vertices of the entity in physical coordinates (vertices are the fastest dimension).
     * @param i_x point in physical coords.
     * @parma o_xi will be set to point in reference coordinates.
     **/
    template <typename T>
    static void phyToRef( t_entityType  i_enType,
                          const T      *i_ves,
                          const T      *i_x,
                                T      *o_xi ) {
     if(      i_enType == QUAD4R ) phyToRefQuad4r((T(*)[4]) i_ves, i_x, o_xi );
     else if( i_enType == TRIA3  ) phyToRefTria3( (T(*)[3]) i_ves, i_x, o_xi );
     else if( i_enType == HEX8R  ) phyToRefHex8r( (T(*)[8]) i_ves, i_x, o_xi );
     else if( i_enType == TET4   ) phyToRefTet4(  (T(*)[4]) i_ves, i_x, o_xi );
     else                          EDGE_LOG_FATAL;
    }

    /**
     * Gets the vertex coordinates of a given local face in the reference element.
     *
     * @param i_entType entity type of the elemenet.
     * @param i_face face id.
     * @param o_coords will be set to coords of the faces, [*][]: dimension, [][*]: vertex
     **/
    static void getFaceCoords( const t_entityType     i_entType,
                                     unsigned short   i_face,
                                     real_mesh       *o_coords ) {
      // get the vertex coords of the face
      for( unsigned int l_dim = 0; l_dim < 3; l_dim++ ) {
        for( unsigned int l_ve = 0; l_ve < C_ENT[i_entType].N_FACE_VERTICES; l_ve++ ) {
          // derive id of vertex coord
          unsigned int l_id  = l_dim       * C_ENT[i_entType].N_FACES * C_ENT[i_entType].N_FACE_VERTICES;
                       l_id += i_face                                 * C_ENT[i_entType].N_FACE_VERTICES;
                       l_id += l_ve;

          o_coords[ l_dim * C_ENT[i_entType].N_FACE_VERTICES + l_ve ] = C_REF_ELEMENT.FA.ENT[i_entType][l_id];
        }
      }
    }

    /**
     * Maps the volume coordinates of one element's face w.r.t.
     * to the ref. element to those of another, neighboring element.
     *
     * Quadrilateral example:
     *
     * Physical:           | Mapping of pt ">*<"   | Mapping of pt ">*<"
     *                     | to ref element for a: | to ref element for b:
     *         a2 b1    b0 | 3            2        | 3            2
     * a3  *************   |   **********          |   **********
     *     *     *     *   |   *        *          |   *       >*<
     *     *    >*<    *   |   *        *          |   *        *
     * a0  *************   |   *       >*<         |   *        *
     *         a1 b2    b3 |   **********          |   **********
     *                     | 0            1        | 0            1
     *
     * Therefore (xi_a, eta_a) cooresponds to (xi_b, eta_b) = (xi_a, 1-eta_a).
     *
     * @param i_entType entity type of the elements.
     * @param i_localFace id of the shared face in the ref. element w.r.t. to the local element.
     * @param i_neighFace id of the shared face in the ref. element w.r.t. to the neigh element.
     * @param i_vertexMap mapping of the vertices w.r.t. to the shared face.
     * @param i_xi 1st element's volume coords in xi-direction.
     * @param i_eta 1st element's volume coords in eta-direction.
     * @param i_zeta 1st element's volume coords in zeta-direction.
     * @param o_xi will be set to 2nd element's volume coords in xi-direction.
     * @param o_eta will bet set to 2nd element's volume coords in eta-direction.
     * @param o_zeta will be set to 2nd element's volume coords in zeta direcation.
     **/
    static void mapFaceCoords( t_entityType    i_entType,
                               unsigned short  i_localFace,
                               unsigned short  i_neighFace,
                               unsigned short  i_vertexMap,
                               real_mesh       i_xi,
                               real_mesh       i_eta,
                               real_mesh       i_zeta,
                               real_mesh      &o_xi,
                               real_mesh      &o_eta,
                               real_mesh      &o_zeta ) {
      // in and output of the mapping
      real_mesh l_coordsIn[3], l_coordsOut[3];
      l_coordsIn[0] = i_xi; l_coordsIn[1] = i_eta; l_coordsIn[2] = i_zeta;

      // rectangular, 8-node hexahedralelements.
      if( i_entType == HEX8R ) {
        EDGE_CHECK_EQ( i_vertexMap, 0 );

        // do the mapping
        map3dQuad4( i_localFace,
                    i_neighFace,
                    l_coordsIn,
                    l_coordsOut );
      }
      // we don't have HEX8R elements
      else {
        // four is sufficient to store all verts for faces of our element types -> hexes
        real_mesh l_veCoordsLoc[3*4];
        real_mesh l_veCoordsNei[3*4];

        // get the face vertices' coords
        getFaceCoords( i_entType,
                       i_localFace,
                       l_veCoordsLoc );

        getFaceCoords( i_entType,
                       i_neighFace,
                       l_veCoordsNei );

        // perform reordering based on vertex mapping for 3D
        EDGE_CHECK_LT( C_ENT[i_entType].N_DIM, 3 );

        // change orientations
        for( unsigned int l_dim = 0; l_dim < 3; l_dim++ ) {
          for( unsigned int l_ve = 0; l_ve < C_ENT[i_entType].N_FACE_VERTICES/2; l_ve++ ) {
            unsigned int l_oldId = l_dim*C_ENT[i_entType].N_FACE_VERTICES + l_ve;
            unsigned int l_newId = l_dim*C_ENT[i_entType].N_FACE_VERTICES + C_ENT[i_entType].N_FACE_VERTICES - 1 - l_ve;

            real_mesh l_tmp = l_veCoordsNei[l_oldId];
            l_veCoordsNei[l_oldId] = l_veCoordsNei[l_newId];
            l_veCoordsNei[l_newId] = l_tmp;
          }
        }

        // do the mapping
        if( i_entType == LINE ) {
          // mapping of points is trivial
          l_coordsOut[0] = l_veCoordsNei[0];
          l_coordsOut[1] = l_veCoordsNei[1];
          l_coordsOut[2] = l_veCoordsNei[2];
        }
        else if( i_entType == QUAD4R || i_entType == TRIA3 ) {
          map2dLine( (real_mesh (*)[2]) l_veCoordsLoc,
                     (real_mesh (*)[2]) l_veCoordsNei,
                                        l_coordsIn,
                                        l_coordsOut );
        }
        else EDGE_LOG_FATAL;
      }

      // return the result
      o_xi   = l_coordsOut[0];
      o_eta  = l_coordsOut[1];
      o_zeta = l_coordsOut[2];
    }

    /**
     * Derives the volume coordinates based on the regular, 4-node quadrilateral's face-coordinates.
     *
     * @param i_fa local face id of the quad.
     * @param i_faPoint face-local coordinates of the point.
     * @param i_volPoint will be set to volume coordinates of the point.
     **/
    template <typename T>
    static void faToVolRefQuad4r( unsigned short i_fa,
                                  T              i_faPoint,
                                  T              o_volPoint[2] ) {
      EDGE_CHECK( i_fa < 4 );

      for( unsigned short l_di = 0; l_di < 2; l_di++ ) {
        // set direction of face
        o_volPoint[l_di]  = C_REF_ELEMENT.FA.QUAD[l_di][i_fa][1];
        o_volPoint[l_di] -= C_REF_ELEMENT.FA.QUAD[l_di][i_fa][0];
        // scale with face parameter
        o_volPoint[l_di] *= i_faPoint;
        // offset
        o_volPoint[l_di] += C_REF_ELEMENT.FA.QUAD[l_di][i_fa][0];
      }
    }

   /**
     * Derives the volume coordinates based on the 3-node triangle's face-coordinates.
     *
     * @param i_fa local face id of the tria.
     * @param i_faPoint face-local coordinates of the point.
     * @param i_volPoint will be set to volume coordinates of the point.
     **/
    template <typename T>
    static void faToVolRefTria3( unsigned short i_fa,
                                 T              i_faPoint,
                                 T              o_volPoint[2] ) {
      EDGE_CHECK( i_fa < 3 );

      for( unsigned short l_di = 0; l_di < 2; l_di++ ) {
        // set direction of face
        o_volPoint[l_di]  = C_REF_ELEMENT.FA.TRIA[l_di][i_fa][1];
        o_volPoint[l_di] -= C_REF_ELEMENT.FA.TRIA[l_di][i_fa][0];
        // scale with face parameter
        o_volPoint[l_di] *= i_faPoint;
        // offset
        o_volPoint[l_di] += C_REF_ELEMENT.FA.TRIA[l_di][i_fa][0];
      }
    }

   /**
     * Derives the volume coordinates based on the rectangular, 8-node hexahedron's face-coordinates.
     *
     * @param i_fa local face id of the hex.
     * @param i_faPoint face-local coordinates of the point.
     * @param i_volPoint will be set to volume coordinates of the point.
     **/
    template <typename T>
    static void faToVolRefHex8r(       unsigned short i_fa,
                                 const T              i_faPoint[2],
                                       T              o_volPoint[3] ) {
      EDGE_CHECK( i_fa < 6 );

      for( unsigned short l_di = 0; l_di < 3; l_di++ ) {
        // contrib of first face-dim
        T l_tmp  = C_REF_ELEMENT.FA.HEX[l_di][i_fa][1];
        l_tmp   -= C_REF_ELEMENT.FA.HEX[l_di][i_fa][0];
        l_tmp   *= i_faPoint[0];
        o_volPoint[l_di] = l_tmp;

        // contrib of second face-dim
        l_tmp    = C_REF_ELEMENT.FA.HEX[l_di][i_fa][3];
        l_tmp   -= C_REF_ELEMENT.FA.HEX[l_di][i_fa][0];
        l_tmp   *= i_faPoint[1];
        o_volPoint[l_di] += l_tmp;

        // shift by origin
        o_volPoint[l_di] += C_REF_ELEMENT.FA.HEX[l_di][i_fa][0];
      }
    }

    /**
     * Derives the volume coordinates based on the 4-node tetrahedron's face-coordinates.
     *
     * @param i_fa local face id of the tet.
     * @param i_faPoint face-local coordinates of the point.
     * @param i_volPoint will be set to volume coordinates of the point.
     **/
    template <typename T>
    static void faToVolRefTet4(       unsigned short i_fa,
                                const T              i_faPoint[2],
                                      T              o_volPoint[3] ) {
      EDGE_CHECK( i_fa < 4 );

      for( unsigned short l_di = 0; l_di < 3; l_di++ ) {
        // contrib of first face-dim
        T l_tmp  = C_REF_ELEMENT.FA.TET[l_di][i_fa][1];
        l_tmp   -= C_REF_ELEMENT.FA.TET[l_di][i_fa][0];
        l_tmp   *= i_faPoint[0];
        o_volPoint[l_di] = l_tmp;

        // contrib of second face-dim
        l_tmp    = C_REF_ELEMENT.FA.TET[l_di][i_fa][2];
        l_tmp   -= C_REF_ELEMENT.FA.TET[l_di][i_fa][0];
        l_tmp   *= i_faPoint[1];
        o_volPoint[l_di] += l_tmp;

        // shift by origin
        o_volPoint[l_di] += C_REF_ELEMENT.FA.TET[l_di][i_fa][0];
      }
    }


   /**
     * Derives the volume coordinates based on the element's face-coordinates.
     *
     * @param i_fa local face id.
     * @param i_elType element type.
     * @param i_faPoint face-local coordinates of the point.
     * @param i_volPoint will be set to volume coordinates of the point.
     **/
    template <typename T>
    static void faToVolRef(       unsigned short  i_fa,
                                  t_entityType    i_elType,
                            const T              *i_faPoint,
                                  T              *o_volPoint ) {
      if( i_elType == LINE ) {
        EDGE_CHECK( i_fa < 2 );
        o_volPoint[0] = C_REF_ELEMENT.FA.LINE[0][i_fa][0];
      }
      else if( i_elType == QUAD4R ) {
        faToVolRefQuad4r( i_fa, i_faPoint[0], o_volPoint );
      }
      else if( i_elType == TRIA3 ) {
        faToVolRefTria3( i_fa, i_faPoint[0], o_volPoint );
      }
      else if( i_elType == HEX8R ) {
        faToVolRefHex8r( i_fa, i_faPoint, o_volPoint );
      }
      else if( i_elType == TET4 ) {
        faToVolRefTet4( i_fa, i_faPoint, o_volPoint );
      }
      else EDGE_LOG_FATAL << "element type " << i_elType << " not supported";
    }

    /**
     * Expresses a point on a tet's face in terms of a neighboring tet's face-coordinates.
     *
     * @param i_ve vertex of the neighboring face assigned to local face's vertex 0.
     * @param i_faPtL point's face-coordinates w.r.t. the local tet.
     * @param o_faPtN will be set to face-coordinates w.r.t. the neighboring tet.
     *
     * Local Coordinates:       Case 1:       Case 2:      Case 3:
     *
     *           2                1            2            0
     *           |\               |\           |\           |\
     *           | \              | \          | \          | \
     *           |  \             |  \         |  \         |  \
     *          u|   \           i|   \        |   \u      t|   \c
     *          a|2  1\          h|2  1\       |2  1\a     a|2  1\h
     *          t|     \         c|     \      |     \t    u|     \i
     *  /|\      |      \         |      \     |      \     |      \
     *   |__\    |___0___\        |___0___\    |___0___\    |___0___\
     *      /    0  chi   1       0  tau   2   1  chi   0   2        1
     *
     *
     *           | 1 |             | 0 |  Case 1        | 0 |         | 1 |
     * chi_new * |   | + tau_new * |   |    =     chi * |   | + tau * |   |
     *           | 0 |             | 1 |                | 1 |         | 0 |
     *
     *                                    Case 2        |-1 |         |-1 |
     *                                      =     chi * |   | + tau * |   |
     *                                                  | 0 |         | 1 |
     *
     *                                    Case 3        | 1 |         | 0 |
     *                                      =     chi * |   | + tau * | 0 |
     *                                                  |-1 |         |-1 |
     **/
    static void faLocToFaNeiTet4( unsigned short       i_ve,
                                  real_mesh      const i_faPtL[2],
                                  real_mesh            o_faPtN[2] ) {
      if(      i_ve == 0 ) {
        o_faPtN[0] = i_faPtL[1];
        o_faPtN[1] = i_faPtL[0];
      }
      else if( i_ve == 1 ) {
        o_faPtN[0] = 1 - i_faPtL[0] - i_faPtL[1];
        o_faPtN[1] = i_faPtL[1];
      }
      else if( i_ve == 2 ) {
        o_faPtN[0] = i_faPtL[0];
        o_faPtN[1] = 1 - i_faPtL[0] - i_faPtL[1];
      }
      else EDGE_LOG_FATAL;
    }

    /**
     * Expresses a point on an element's face in terms of a neighboring element's face-coordinates.
     *
     * @param i_elType element type.
     * @param i_faPtL point's face-coordinates w.r.t. the local element.
     * @param o_faPtN will be set to face-coordinates w.r.t. the neighboring element.
     * @param i_ve vertex of of the neighboring face assigned to local face's vertex 0 (if relevant)
     **/
    static void faLocToFaNei( t_entityType          i_elType,
                              real_mesh      const *i_faPtL,
                              real_mesh            *o_faPtN,
                              unsigned short        i_ve = std::numeric_limits< unsigned short >::max() ) {
      EDGE_CHECK( i_elType != POINT );

      if( C_ENT[i_elType].N_DIM < 3 ) {
        o_faPtN[0] = 1 - i_faPtL[0];
      }
      else if( i_elType == HEX8R ) {
        o_faPtN[0] = i_faPtL[1];
        o_faPtN[1] = i_faPtL[0];
      }
      else if( i_elType == TET4 ) {
        faLocToFaNeiTet4( i_ve, i_faPtL, o_faPtN );
      }
    }


};

#endif
