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
 * Geometry functions.
 **/

#ifndef EDGE_LINALG_GEOM_HPP
#define EDGE_LINALG_GEOM_HPP

#include <cassert>
#include <cmath>
#include "submodules/wykobi/wykobi.hpp"
#include "constants.hpp"
#include "io/logging.h"
#include "monitor/instrument.hpp"
#include "Matrix.h"
#include "Mappings.hpp"

namespace edge {
  namespace linalg {
    class Geom;

    template< unsigned short TL_N_DIS >
    class GeomT;

    template< unsigned short TL_N_DIS >
    class GeomTs;

    template<>
    class GeomTs<2>;

    template<>
    class GeomTs<3>;
  }
}

class edge::linalg::Geom {
  public:
    /**
     * Computes the scalar product of the two 1-entry "vectors" a and b.
     *
     * @param i_a vector a.
     * @param i_b vector b.
     * @return scalar product.
     *
     * @paramt T floating point type.
     **/
    template <typename T>
    static T sprod1( const T i_a[1],
                     const T i_b[1] ) {
      return i_a[0] * i_b[0];
    }

    /**
     * Computes the scalar product of the two 2-entry vectors a and b.
     *
     * @param i_a vector a.
     * @param i_b vector b.
     * @return scalar product.
     *
     * @paramt T floating point type.
     **/
    template <typename T>
    static T sprod2( const T i_a[2],
                     const T i_b[2] ) {
      return i_a[0] * i_b[0] + i_a[1] * i_b[1];
    }

    /**
     * Computes the scalar product of the two 3-entry vectors a and b.
     *
     * @param i_a vector a.
     * @param i_b vector b.
     * @return scalar product.
     *
     * @paramt T floating point type.
     **/
    template <typename T>
    static T sprod3( const T i_a[3],
                     const T i_b[3] ) {
      return i_a[0] * i_b[0] + i_a[1] * i_b[1] + i_a[2] * i_b[2];
    }

    /**
     * Computes the norm of the two 2-entry vectors vec.
     *
     * @param i_vec vector.
     * @return norm
     *
     * @paramt T floating point type.
     **/
    template <typename T>
    static T norm2( const T i_vec[2] ) {
      return std::sqrt( sprod2( i_vec, i_vec) );
    }

    /**
     * Normalizes the given 1D "vector" to unit length.
     *
     * @param io_vec vector which will be normalized.
     *
     * @paramt T floating point type.
     **/
    template <typename T>
    static void normalize1( T io_vec[1] ) {
      if(      io_vec[0] > 0 ) io_vec[0] =  1;
      else if( io_vec[0] < 0 ) io_vec[0] = -1;
      else EDGE_LOG_FATAL << "can't normalize zero vector";
    }

    /**
     * Normalizes the 2D given vector to unit length.
     *
     * @param io_vec vector which will be normalized.
     *
     * @paramt T floating point type.
     **/
    template <typename T>
    static void normalize2( T io_vec[2] ) {
      // determine norm
      T l_norm = norm2( io_vec );
      EDGE_CHECK( l_norm > TOL.MESH );
      for( unsigned short l_di = 0; l_di < 2; l_di++ ) {
        io_vec[l_di] /= l_norm;
      }
    }

    /**
     * Normalizes the 3D given vector to unit length.
     *
     * @param io_vec vector which will be normalized.
     *
     * @paramt T floating point type.
     **/
    template <typename T>
    static void normalize3( T io_vec[3] ) {
      // determine norm
      T l_norm = norm3( io_vec );
      for( unsigned short l_di = 0; l_di < 3; l_di++ ) {
        io_vec[l_di] /= l_norm;
      }
    }

    /**
     * Computes the norm of the two 3-entry vectors vec.
     *
     * @param i_vec vector.
     * @return norm
     *
     * @paramt T floating point type.
     **/
    template <typename T>
    static T norm3( const T i_vec[3] ) {
      return std::sqrt( sprod3(i_vec, i_vec) );
    }

    /**
     * Computes the cross product between of the two vectors a and b in 3D.
     *
     * @param i_a vector a.
     * @param i_b vector b.
     * @param o_cross resulting cross product.
     *
     * @paramt T floating point type.
     **/
    template <typename T>
    static void cross3( const T i_a[3],
                        const T i_b[3],
                              T o_cross[3] ) {
      o_cross[0] = i_a[1]*i_b[2] - i_a[2]*i_b[1];
      o_cross[1] = i_a[2]*i_b[0] - i_a[0]*i_b[2];
      o_cross[2] = i_a[0]*i_b[1] - i_a[1]*i_b[0];
    }

    /**
     * Computes the insphere diameter of the given entity.
     *
     * @param i_entType type of the entity.
     * @param i_vertices vertices of the entity which diameter is computed, [*][]: dimensions, [][*]: vertices.
     * @param i_nDims number of considered dimensions (see i_vertices also). Default is 4 which selects the lowest possible number of dimensions.
     *
     * @paramt T floating point type.
     **/
    template <typename T>
    static T inDia(       t_entityType    i_entType,
                    const T              *i_vertices,
                          unsigned short  i_nDim = 4 ) {
      // diameter which is computed
      T l_dia = 0;

      // assign dimensions
      unsigned short l_nDim = C_ENT[i_entType].N_DIM;
      if( i_nDim < 4 ) l_nDim = i_nDim;

      /*
       * Tria3, dim 2
       */
      if( i_entType == TRIA3 && l_nDim == 2 ) {
        T l_faceAreas[3];

        // compute areas of the faces
        for( unsigned short l_fa = 0; l_fa < 3; l_fa++ ) {
          T l_faVes[2][2];
          for( unsigned short l_di = 0; l_di < 2; l_di++ )
            for( unsigned short l_ve = 0; l_ve < 2; l_ve++ )
              l_faVes[l_di][l_ve] = i_vertices[l_di*3+(l_fa+l_ve)%3];
          l_faceAreas[l_fa] = volume( LINE, l_faVes[0], l_nDim );
        }

        T l_volume = volume( TRIA3, i_vertices, l_nDim );
        l_dia  = (T) 4 * l_volume;
        l_dia /= l_faceAreas[0] + l_faceAreas[1] + l_faceAreas[2];
      }

      /*
       * Tet4, dim 3
       *
       * Reference: John Burkardt
       *            The Inscribed Sphere of a Tetrahedron
       *            Computational Geometry Lab: TETRAHEDRONS
       *            2010
       */
      else if( i_entType == TET4 && l_nDim == 3 ) {
        // define lambda for computation of normal vector (faces) length
        std::function< T( const T[3], const T[3], const T[3]) > l_faNorm =
          []( const T i_a[3], const T i_b[3], const T i_c[3] ) {
          T l_ba[3], l_ca[3];
          for( unsigned short l_di = 0; l_di < 3; l_di++ ) {
            l_ba[l_di] = i_b[l_di] - i_a[l_di];
            l_ca[l_di] = i_c[l_di] - i_a[l_di];
          }
          T l_cross[3];
          cross3( l_ba, l_ca, l_cross );

          return norm3( l_cross );
        };

        // convert soa to aos
        T l_veAos[4*3];
        for( unsigned int l_ve = 0; l_ve < 4; l_ve++ )
          for( unsigned int l_dim = 0; l_dim < 3; l_dim++ )
            l_veAos[l_ve*3 + l_dim] = i_vertices[l_dim*4 + l_ve];

        l_dia += l_faNorm( l_veAos+3*0, l_veAos+3*1, l_veAos+3*2 );
        l_dia += l_faNorm( l_veAos+3*0, l_veAos+3*1, l_veAos+3*3 );
        l_dia += l_faNorm( l_veAos+3*0, l_veAos+3*2, l_veAos+3*3 );
        l_dia += l_faNorm( l_veAos+3*1, l_veAos+3*2, l_veAos+3*3 );
        l_dia =  volume( TET4, i_vertices )*6 / l_dia;
        l_dia *= 2;
      }
      else EDGE_LOG_FATAL << "unknown setting";

      EDGE_CHECK( l_dia > 0 );

      return l_dia;
    }

    /**
     * Computes the volume of the given entity.
     *
     * @param i_entType type of the entity.
     * @param i_vertices vertices of the entity which volume is computed, [*][]: dimensions, [][*]: vertices.
     * @param i_nDims number of considered dimensions (see i_vertices also). Default is 4 which selects the lowest possible number of dimensions.
     *
     * @paramt T floating point type.
     **/
    template <typename T>
    static T volume(       t_entityType    i_entType,
                     const T              *i_vertices,
                           unsigned short  i_nDim = 4 ) {
      // volume which is computed
      T l_volume = 0;

      // assign dimensions
      unsigned short l_nDim = C_ENT[i_entType].N_DIM;
      if( i_nDim < 4 ) l_nDim = i_nDim;

      if( i_entType == LINE && l_nDim == 2 ) {
        T l_x = i_vertices[0] - i_vertices[1];
        T l_y = i_vertices[2] - i_vertices[3];

        l_volume = std::sqrt( l_x * l_x + l_y * l_y );
      }
      else if( i_entType == TRIA3 && l_nDim == 2 ) {
        // http://mathworld.wolfram.com/TriangleArea.html
        l_volume = -i_vertices[0*3+1]*i_vertices[1*3+0] +
                    i_vertices[0*3+2]*i_vertices[1*3+0] +
                    i_vertices[0*3+0]*i_vertices[1*3+1] +
                   -i_vertices[0*3+2]*i_vertices[1*3+1] +
                   -i_vertices[0*3+0]*i_vertices[1*3+2] +
                    i_vertices[0*3+1]*i_vertices[1*3+2];
         l_volume = std::abs( (T) 0.5 * l_volume );
      }
      else if( i_entType == TRIA3 && l_nDim == 3 ) {
        // get edges
        T l_a[3], l_b[3];
        l_a[0] = i_vertices[1] - i_vertices[0];
        l_a[1] = i_vertices[4] - i_vertices[3];
        l_a[2] = i_vertices[7] - i_vertices[6];

        l_b[0] = i_vertices[2] - i_vertices[0];
        l_b[1] = i_vertices[5] - i_vertices[3];
        l_b[2] = i_vertices[8] - i_vertices[6];

        // compute cross product
        T l_cross[3];
        cross3( l_a, l_b, l_cross );

        // compute volume
        l_volume = norm3(l_cross) / 2;
      }
      else if( i_entType == TET4 ) {
        EDGE_CHECK( l_nDim >= 3 );

        // determinant is six times the volume of the det
        T l_jac[3][3];
        Mappings::evalJac( TET4, i_vertices, l_jac[0] );
        T l_det = linalg::Matrix::det3x3( l_jac );

        // det might be negative since we call this in general before enforcing counterclockwise tets
        l_volume = std::abs(l_det) / (T) 6;
      }
      else EDGE_LOG_FATAL << "unknown setting: " << i_entType << " " << l_nDim;

      return l_volume;
    }

    /**
     * Computes the tangents for a entity.
     *
     * @param i_entType type of the entity which tangens are computed.
     * @param i_vertices vertices of the entity which tangents are computed, [*][]: dimensions, [][*]: vertices.
     * @param o_tangent0 will be set to first tangent.
     * @param o_tangent1 will be set to second tangent.
     *
     * @paramt TL_T_REAL floating point type.
     **/
    template < typename TL_T_REAL >
    static void computeTangents(       t_entityType         i_entType,
                                 const TL_T_REAL    * const i_vertices,
                                       TL_T_REAL            o_tangent0[3],
                                       TL_T_REAL            o_tangent1[3] ) {
      unsigned short l_nVe = C_ENT[i_entType].N_VERTICES;

      if( i_entType == LINE ) {
        o_tangent0[0] = i_vertices[0*l_nVe + 1] - i_vertices[0*l_nVe + 0];
        o_tangent0[1] = i_vertices[1*l_nVe + 1] - i_vertices[1*l_nVe + 0];
        o_tangent0[2] = i_vertices[2*l_nVe + 1] - i_vertices[2*l_nVe + 0];

        o_tangent1[0] = 0;
        o_tangent1[1] = 0;
        o_tangent1[2] = 0;   
      }
      else if( i_entType == TRIA3 ) {
        /*
         * we can derive the second (normal) tangent from the first tangent
         *
         *            x <- v2
         * v1 -> x   *
         *       *  *
         *       * *
         *       x******x <- n: that's the vector we want
         *
         *  <v1, n> = 0
         *  <v1, v2 + alpha * v1> = 0
         *  <v1, v2> + <v1, alpha * v1> = 0
         *  => alpha = - <v1,v2> / <v1,v1>
         *  => n = v2 + alpha * v1
         */
        TL_T_REAL l_v1[3], l_v2[3], l_n[3];

        for( unsigned int l_dim = 0; l_dim < 3; l_dim++ ) {
          l_v1[l_dim] = i_vertices[l_dim*l_nVe + 1] - i_vertices[l_dim*l_nVe + 0];
          l_v2[l_dim] = i_vertices[l_dim*l_nVe + 2] - i_vertices[l_dim*l_nVe + 0];
        }

        // get scalar products
        TL_T_REAL l_s1 = 0;
        TL_T_REAL l_s2 = 0;

        for( unsigned int l_dim = 0; l_dim < 3; l_dim++ ) {
          l_s1 += l_v1[l_dim] * l_v2[l_dim];
          l_s2 += l_v1[l_dim] * l_v1[l_dim];
        }

        assert( std::abs(l_s2) > TOL.MESH );

        TL_T_REAL l_alpha = -l_s1 / l_s2;

        // set the normal
        for( unsigned int l_dim = 0; l_dim < 3; l_dim++ ) {
          l_n[l_dim] = l_v2[l_dim] + l_alpha * l_v1[l_dim]; 
        }

        // get length of tangents
        TL_T_REAL l_lens[2] = {0, 0};
        for( unsigned int l_dim = 0; l_dim < 3; l_dim++ ) {
          l_lens[0] += l_v1[l_dim] * l_v1[l_dim];
          l_lens[1] += l_n[l_dim]  * l_n[l_dim];
        }
        l_lens[0] = std::sqrt( l_lens[0] );
        l_lens[1] = std::sqrt( l_lens[1] );

        assert( l_lens[0] > TOL.MESH );
        assert( l_lens[1] > TOL.MESH );

        // set tangents
        for( unsigned int l_dim = 0; l_dim < 3; l_dim++ ) {
          o_tangent0[l_dim] = l_v1[l_dim] / l_lens[0];
          o_tangent1[l_dim] = l_n[l_dim]  / l_lens[1];
        }
      }
      else assert( false );
    }

    /**
     * Computes the directed, outer pointing normal to a entity.
     *
     * @param i_entType entity type.
     * @param i_vertices vertices of the entity. [*][]: dimensions, [][*]: vertices.
     * @param i_normalPoint additional point lying on the non-normal side of the entity.
     * @param o_outNormal will bet set to the outer pointing normal. origin: (0,0,0), length: 1.0
     *
     * @paramt TL_T_REAL floating point type.
     **/
    template < typename TL_T_REAL >
    static void computeOutPtNormal(       t_entityType  i_entType,
                                    const TL_T_REAL    *i_vertices,
                                    const TL_T_REAL     i_normalPoint[3],
                                          TL_T_REAL     o_outNormal[3] ) {
      unsigned short l_nVe = C_ENT[i_entType].N_VERTICES;

      if( i_entType == LINE ) {
        // derive tangent
        TL_T_REAL l_tX = i_vertices[0*l_nVe + 0] - i_vertices[0*l_nVe + 1];
        TL_T_REAL l_tY = i_vertices[1*l_nVe + 0] - i_vertices[1*l_nVe + 1];

        // get length of tangent (norm)
        TL_T_REAL l_length = std::sqrt(l_tX * l_tX + l_tY * l_tY);

        // scale vector (unit length)
        l_tX /= l_length; l_tY /= l_length;

        // compute normal vector
        TL_T_REAL l_nX = l_tY;
        TL_T_REAL l_nY = -l_tX;

        // check that this is a normal (dot product == 0)
        EDGE_CHECK( std::abs(l_nX * l_tX + l_nY * l_tY) < TOL.MESH );

        // check for length one
        EDGE_CHECK( std::abs(std::sqrt( l_nX * l_nX + l_nY * l_nY ) - 1) < TOL.MESH );

        /*
         * compute dot product of vector pointing from a vertex to the normal point and the normal
         *
         *            x * * * * * o <-- normal point
         *          *  *     .   o
         *        *     *  .  <-----example normal
         *      x     1  *     o <-- vector pointing from a vertex
         *         *      *   o
         *             *   * o
         *                * x
         *  translates to:
         *  |        o
         *  |       o \
         *  |      o   \
         *  |     o   a . <--- projection (a=90deg) by the dot product > 0:  
         *  |    o    .        we are on the wrong side here.
         *  |   o   .
         *  |  o  .
         *  | o .
         *  |o.
         *  |_______
         */
        TL_T_REAL l_dP =   l_nX * (i_normalPoint[0] - i_vertices[0*l_nVe + 0])
                         + l_nY * (i_normalPoint[1] - i_vertices[1*l_nVe + 0]);

        // if the dot product is positive, the angle is below 90deg:
        // we want the normal to point in the other direction; therefore we have to change the sign
        if( l_dP > 0 ) {
          l_nX = -l_nX;
          l_nY = -l_nY;
        }

        // time to return the values and keep going
        o_outNormal[0] = l_nX;
        o_outNormal[1] = l_nY;
        o_outNormal[2] = 0;
      }
      else if( i_entType == TRIA3 ) {
        // compute tangents
        TL_T_REAL l_t1[3], l_t2[3];
        computeTangents( T_SDISC.FACE, i_vertices, l_t1, l_t2 );

        // compute normal via cross product
        TL_T_REAL l_nX = ( l_t1[1] * l_t2[2] ) - ( l_t1[2] * l_t2[1] );
        TL_T_REAL l_nY = ( l_t1[2] * l_t2[0] ) - ( l_t1[0] * l_t2[2] );
        TL_T_REAL l_nZ = ( l_t1[0] * l_t2[1] ) - ( l_t1[1] * l_t2[0] );

        // compute length
        TL_T_REAL l_nLength = std::sqrt( l_nX*l_nX + l_nY*l_nY + l_nZ*l_nZ );
        EDGE_CHECK( l_nLength > TOL.MESH );

        // normalize
        l_nX /= l_nLength; l_nY /= l_nLength; l_nZ /= l_nLength;

        // compute dot product with normal point for right orientation
        TL_T_REAL l_dP =   l_nX * (i_normalPoint[0] - i_vertices[0*l_nVe + 0])
                         + l_nY * (i_normalPoint[1] - i_vertices[1*l_nVe + 0])
                         + l_nZ * (i_normalPoint[2] - i_vertices[2*l_nVe + 0]);
        EDGE_CHECK( std::abs(l_dP) > TOL.MESH );

        // if the dot product is positive, the angle is below 90deg:
        // we want the normal to point in the other direction; therefore we have to change the sign
        if( l_dP > 0 ) {
          l_nX = -l_nX;
          l_nY = -l_nY;
          l_nZ = -l_nZ;
        }

        // set result
        o_outNormal[0] = l_nX;
        o_outNormal[1] = l_nY;
        o_outNormal[2] = l_nZ;
      }
      else EDGE_LOG_FATAL;
    }

    /**
     * Computes the sub-face normal and sub-face tangents for the two additional
     * face types, introduced in a tetrahedral sub-grid.
     * In the reference element, the normals of these face-types are given as:
     *   0) (sqrt(2), sqrt(2), 0)
     *   1) (sqrt(2), 0, sqrt(2))
     *
     * @param i_et entity type.
     * @param i_ty additional face type for which the normal and tangets are requested.
     * @param i_veCrds vertex coordinates [*][] dimension, [][*] vertex.
     * @param o_n will be set to normal of face type.
     * @param o_t0 will be set to first tangent of additional face type.
     * @param o_t1 will be set to second tangent of additional face type.
     * @param o_a will be set to area of the face.
     **/
    template< typename TL_T_REAL >
    static void sfAdd( t_entityType           i_et,
                       unsigned short         i_ty,
                       TL_T_REAL      const * i_veCrds,
                       TL_T_REAL              o_n[3],
                       TL_T_REAL              o_t0[3],
                       TL_T_REAL              o_t1[3],
                       TL_T_REAL            & o_a ) {
      EDGE_CHECK_EQ( i_et, TET4 );
      EDGE_CHECK_LE( i_ty, 1 );

      // by definition "left" is the side of vertex 0
      TL_T_REAL l_np[3];
      for( unsigned short l_di = 0; l_di < 3; l_di++ )
        l_np[l_di] = i_veCrds[ l_di*4 + 0 ];

      // assemble pseudo face
      TL_T_REAL l_faVes[3][3];

      if( i_ty == 0 ) {
        for( unsigned short l_di = 0; l_di < 3; l_di++ ) {
          l_faVes[l_di][0] = i_veCrds[l_di*4 + 1];
          l_faVes[l_di][1] = i_veCrds[l_di*4 + 2];

          l_faVes[l_di][2]  = i_veCrds[l_di*4 + 2];
          l_faVes[l_di][2] += i_veCrds[l_di*4 + 3];
          l_faVes[l_di][2] -= i_veCrds[l_di*4 + 0];
        }
      }
      else {
        for( unsigned short l_di = 0; l_di < 3; l_di++ ) {
          l_faVes[l_di][0]  = i_veCrds[l_di*4 + 1];
          l_faVes[l_di][1]  = i_veCrds[l_di*4 + 2];
          l_faVes[l_di][1] += i_veCrds[l_di*4 + 1];
          l_faVes[l_di][1] -= i_veCrds[l_di*4 + 0];

          l_faVes[l_di][2]  = i_veCrds[l_di*4 + 2];
          l_faVes[l_di][2] += i_veCrds[l_di*4 + 3];
          l_faVes[l_di][2] -= i_veCrds[l_di*4 + 0];
        }
      }

      linalg::Geom::computeOutPtNormal( C_ENT[i_et].TYPE_FACES,
                                        l_faVes[0],
                                        l_np,
                                        o_n );

      linalg::Geom::computeTangents( C_ENT[i_et].TYPE_FACES,
                                     l_faVes[0],
                                     o_t0,
                                     o_t1 );

      o_a = linalg::Geom::volume( C_ENT[i_et].TYPE_FACES,
                                  l_faVes[0],
                                  3 );
    }

    /**
     * Derives if a given point is inside or outside an entity.
     *   Remark: For 3D ents we assume counter-clockwise storage of the face's vertices.
     *
     * @param i_entType entity type.
     * @param i_vertices coords of the entities vertices. [*][]: dimensions, [][*]: vertices.
     * @param i_pt coordinates of the point in question.
     * @return 1 if the point is inside the entity, 2 if on a face, 0 otherwise.
     *
     * @paramt TL_T_REAL floating point type.
     **/
    template< typename TL_T_REAL>
    static unsigned short inside(       t_entityType  i_entType,
                                  const TL_T_REAL    *i_vertices,
                                  const TL_T_REAL    *i_pt ) {
      PP_INSTR_FUN("inside")

      bool l_inside  = true;
      bool l_outside = false;
      bool l_face    = false;

      // define directed vector from entities's face-verts to pt
      assert(  C_ENT[i_entType].N_FACE_VERTICES < 4 ||
              (C_ENT[i_entType].N_FACE_VERTICES == 4 && i_entType == HEX8R) );
      TL_T_REAL l_directed[3*4];
      for( unsigned short l_dv = 0; l_dv < 3*4; l_dv++ ) l_directed[l_dv] = std::numeric_limits< TL_T_REAL >::max();

      unsigned short l_nDims     = C_ENT[i_entType].N_DIM;
      unsigned short l_nEnVes    = C_ENT[i_entType].N_VERTICES;
      unsigned short l_nFaVes    = C_ENT[i_entType].N_FACE_VERTICES;
      // only consider 3 verts per face for hexes
      unsigned short l_nFaVesRed = std::min( (unsigned short) 3, l_nFaVes );

      unsigned short l_nFas      = C_ENT[i_entType].N_FACES;

      // dermine surrounding box of the entity
      TL_T_REAL l_box[3][2];
      for( unsigned short l_di = 0; l_di < 3; l_di++ ) {
        l_box[l_di][0] = std::numeric_limits< TL_T_REAL >::max();
        l_box[l_di][1] = std::numeric_limits< TL_T_REAL >::lowest();
      }
      for( unsigned short l_di = 0; l_di < l_nDims; l_di++ ) {
        for( unsigned short l_ve = 0; l_ve < l_nEnVes; l_ve++ ) {
          l_box[l_di][0] = std::min( l_box[l_di][0], i_vertices[l_di*l_nEnVes + l_ve] );
          l_box[l_di][1] = std::max( l_box[l_di][1], i_vertices[l_di*l_nEnVes + l_ve] );
        }
      }
      // shortcut if point is not in the surrounding box
      for( unsigned short l_di = 0; l_di < l_nDims; l_di++ ) {
        if( l_box[l_di][0] > i_pt[l_di]+TOL.LINALG ||
            l_box[l_di][1] < i_pt[l_di]-TOL.LINALG ) return 0;
      }

      // line elements only require the point in question to be in between the vertices
      if( C_ENT[i_entType].N_FACE_VERTICES == 1 ) {
        l_inside = (i_pt[0] > i_vertices[0]+TOL.LINALG) &&
                   (i_pt[0] < i_vertices[1]-TOL.LINALG);

        l_face =           std::abs(i_pt[0] - i_vertices[0]) < TOL.LINALG;
        l_face = l_face || std::abs(i_pt[0] - i_vertices[1]) < TOL.LINALG;
      }
      // all other elements require a positive determinants of the faces
      else {
        for( unsigned short l_fa = 0; l_fa < l_nFas; l_fa++ ) {
          // derive directed vectors for this face
          for( unsigned short l_ve = 0; l_ve < l_nFaVesRed; l_ve++ ) {
            unsigned short l_veId = C_REF_ELEMENT.FA_VE_CC.ENT[i_entType][l_fa*l_nFaVes + l_ve];
            for( unsigned short l_dm = 0; l_dm < l_nDims; l_dm++ ) {
              l_directed[l_dm*l_nFaVesRed + l_ve] = i_vertices[l_dm*l_nEnVes + l_veId] - i_pt[l_dm];
            }
          }

          TL_T_REAL l_det = 0;
          if( l_nFaVesRed == 2 ) {
            l_det = Matrix::det2x2( (TL_T_REAL (*)[2]) l_directed );
          }
          else if( l_nFaVesRed == 3 ) {
            l_det = Matrix::det3x3( (TL_T_REAL (*)[3]) l_directed );
          }
          else EDGE_LOG_FATAL;
          l_inside  = l_inside  && ( l_det >  TOL.LINALG );
          l_outside = l_outside || ( l_det < -TOL.LINALG );
          l_face    = l_face    ||   std::abs( l_det ) < TOL.LINALG;
        }
      }

      assert( !(l_inside == true && l_face == true) );

      unsigned short l_result = 0;
      if( l_inside             ) l_result = 1;
      if( l_face && !l_outside ) l_result = 2;

      return l_result;
    }

    /**
     * @brief Projects a point outside of an entity to the surface of the entity.
     *        If the point is inside the entity, nothing is done.
     *
     * @param i_enType entity type.
     * @param i_veCrds vertex coordinates of the entity.
     * @param io_pt will be used as input and updated accordingly.
     */
    template< typename TL_T_REAL >
    static void closestPoint( t_entityType         i_enType,
                              TL_T_REAL     const *i_veCrds,
                              TL_T_REAL           *io_pt ) {
      unsigned short l_nVes = C_ENT[i_enType].N_VERTICES;

      if( i_enType == LINE || i_enType == QUAD4R || i_enType == HEX8R ) {
        // number of dimensions
        unsigned short l_nDis = C_ENT[i_enType].N_DIM;

        // derive minimum and maximum vertex coordinates
        TL_T_REAL l_minMax[2][3];
        for( unsigned short l_di = 0; l_di < 3; l_di++ ) {
          l_minMax[0][l_di] = std::numeric_limits< TL_T_REAL >::max();
          l_minMax[1][l_di] = std::numeric_limits< TL_T_REAL >::lowest();
        }

        for( unsigned short l_ve = 0; l_ve < l_nVes; l_ve++ ) {
          for( unsigned short l_di = 0; l_di < l_nDis; l_di++ ) {
            l_minMax[0][l_di] = std::min( l_minMax[0][l_di],
                                          i_veCrds[l_di*l_nVes+l_ve] );
            l_minMax[1][l_di] = std::max( l_minMax[1][l_di],
                                          i_veCrds[l_di*l_nVes+l_ve] );
          }
        }

        // crop by cube
        for( unsigned short l_di  = 0; l_di < l_nDis; l_di++ ) {
          io_pt[l_di] = std::min( io_pt[l_di], l_minMax[1][l_di] );
          io_pt[l_di] = std::max( io_pt[l_di], l_minMax[0][l_di] );
        }
      }
      else if( i_enType == TRIA3 ) {
        // init triangle
        wykobi::triangle< TL_T_REAL, 2 > l_tria;
        for( unsigned short l_di = 0; l_di < 2; l_di++ )
          for( unsigned short l_ve = 0; l_ve < 3; l_ve++ )
            l_tria[l_ve][l_di] = i_veCrds[l_di*l_nVes+l_ve];

        // compute closest point
        wykobi::point2d< TL_T_REAL > l_pt;
        l_pt = wykobi::closest_point_on_triangle_from_point( l_tria,
                                                             io_pt[0], io_pt[1] );

        // save result
        for( unsigned short l_di = 0; l_di < 2; l_di++ )
          io_pt[l_di] = l_pt[l_di];
      }
      else if( i_enType == TET4 ) {
        // query faces if the point is not inside
        if( inside( i_enType, i_veCrds, io_pt ) == 0 ) {
          // minimum distance to point
          TL_T_REAL l_distMin = std::numeric_limits< TL_T_REAL >::max();

          wykobi::point3d< TL_T_REAL > l_ptIn;
          for( unsigned short l_di = 0; l_di < 3; l_di++ )
            l_ptIn[l_di] = io_pt[l_di];

          // iterate over faces
          for( unsigned short l_fa = 0; l_fa < 4; l_fa++ ) {
            // assemble triangle
            wykobi::triangle< TL_T_REAL, 3 > l_tria;
            for( unsigned short l_ve = 0; l_ve < 3; l_ve++ ) {
              unsigned short l_veId = C_REF_ELEMENT.FA_VE_CC.TET[l_fa][l_ve];

              for( unsigned short l_di = 0; l_di < 3; l_di++ )
                l_tria[l_ve][l_di] = i_veCrds[l_di*l_nVes+l_veId];
            }

            // compute closest point
            wykobi::point3d< TL_T_REAL > l_pt;
            l_pt = wykobi::closest_point_on_triangle_from_point( l_tria,
                                                                 l_ptIn );

            // update point for a new minimum
            TL_T_REAL l_dist = wykobi::distance( l_ptIn, l_pt );
            if( l_dist < l_distMin ) {
              l_distMin = l_dist;
              for( unsigned short l_di = 0; l_di < 3; l_di++ )
                io_pt[l_di] = l_pt[l_di];
            }
          }
        }
      }
      else EDGE_LOG_FATAL;
    }
};

template< unsigned short TL_N_DIS >
class edge::linalg::GeomT {
  public:
    /**
     * @brief Computes the scalar product of the two vectors v1 and v2.
     *
     * @param i_v1 first vector.
     * @param i_v2 second vecto.
     * @return scalar product of the two vectors.
     *
     * @tparam TL_T_REAL floating point type.
     */
    template< typename TL_T_REAL >
    static TL_T_REAL sprod( TL_T_REAL const i_v1[TL_N_DIS],
                            TL_T_REAL const i_v2[TL_N_DIS] ) {
      TL_T_REAL l_sp = 0;
      for( unsigned short l_di = 0; l_di < TL_N_DIS; l_di++ )
        l_sp += i_v1[l_di] * i_v2[l_di];

      return l_sp;
    }

    /**
     * @brief Computes the L2 norm of the given vector
     *
     * @param i_vec vector for which the norm is computed.
     * @return L2 norm of the given vector.
     *
     * @paramt floating point type.
     */
    template< typename TL_T_REAL >
    static TL_T_REAL normL2( TL_T_REAL const i_vec[TL_N_DIS] ) {
      TL_T_REAL l_l2 = 0;

      for( unsigned short l_di = 0; l_di < TL_N_DIS; l_di++ )
        l_l2 += i_vec[l_di] * i_vec[l_di];
      l_l2 = std::sqrt( l_l2 );

      return l_l2;
    }

    /**
     * Computes the norm of the vector given by the two points.
     *
     * @param i_pt0 first point.
     * @param i_pt1 second point.
     *
     * @paramt TL_T_REAL floating point type.
     **/
    template< typename TL_T_REAL  >
    static TL_T_REAL norm( TL_T_REAL const i_pt0[TL_N_DIS],
                           TL_T_REAL const i_pt1[TL_N_DIS] ) {

      TL_T_REAL l_length = 0;

      for( unsigned short l_di = 0; l_di < TL_N_DIS; l_di++ ) {
        TL_T_REAL l_diff = i_pt1[l_di] - i_pt0[l_di];
        l_length += l_diff * l_diff;
      }

      l_length = std::sqrt( l_length );
      return l_length;
    }
};

template<>
class edge::linalg::GeomTs< 2 > {
  public:
    /**
     * @brief Computes the angle between two vectors.
     *
     * @param i_v0 first vector.
     * @param i_v1 second vector.
     * @return angle between the two vectors.
     *
     * @paramt TL_T_REAL floating point type.
     */
    template< typename TL_T_REAL >
    static TL_T_REAL angle( TL_T_REAL const i_v0[2],
                            TL_T_REAL const i_v1[2] ) {
      // derive angle in degree
      TL_T_REAL l_an  = wykobi::robust_cartesian_angle( i_v1[0], i_v1[1] );
                l_an -= wykobi::robust_cartesian_angle( i_v0[0], i_v0[1] );
      if( l_an < 0 ) l_an += TL_T_REAL( 360 );

      // scale to radians
      TL_T_REAL l_sca = (2.0*M_PI) / 360.0;
      l_an *= l_sca;

      return l_an;
    }
};

template<>
class edge::linalg::GeomTs< 3 > {
  public:
    /**
     * Computes the norm of the vector given by the two points.
     *
     * @param i_pt0 first point.
     * @param i_pt1 second point.
     * @return norm of the vector.
     *
     * @paramt TL_T_REAL floating point type.
     **/
    template< typename TL_T_REAL  >
    static TL_T_REAL norm( TL_T_REAL const i_pt0[3],
                           TL_T_REAL const i_pt1[3] ) {

      TL_T_REAL l_diff[3];

      for( unsigned short l_di = 0; l_di < 3; l_di++ ) {
        l_diff[l_di] = i_pt1[l_di] - i_pt0[l_di];
      }

      return Geom::norm3( l_diff );
    }

    /**
     * Derives the matrix which rotates a unit vector in v0-direction to a unit vector in v1-direction.
     * Remark: v0 and v1 are not allowed to be orthogonal (zero scalar product).
     *
     * @param i_v0 vector in v0-direction.
     * @param i_v1 vector in v1-direction.
     * @param o_rm will be set to rotation matrix.
     *
     * @paramt TL_T_REAL precision of the arithmetic operations.
     **/
    template< typename TL_T_REAL  >
    static void rotMat( TL_T_REAL const i_v0[3],
                        TL_T_REAL const i_v1[3],
                        TL_T_REAL       o_rm[3][3] ) {
      // copy vectors to local storage
      TL_T_REAL l_v0[3], l_v1[3];
      for( unsigned short l_di = 0; l_di< 3; l_di++ ) {
        l_v0[l_di] = i_v0[l_di];
        l_v1[l_di] = i_v1[l_di];
      }

      // normalize the two vectors
      Geom::normalize3( l_v0 );
      Geom::normalize3( l_v1 );

      // get the cross product
      TL_T_REAL l_cross[3];
      Geom::cross3( l_v0, l_v1, l_cross );

      // assemble skew symmetric cross product matrix
      TL_T_REAL l_sm[3][3];
      l_sm[0][0] = 0;
      l_sm[0][1] = -l_cross[2];
      l_sm[0][2] =  l_cross[1];

      l_sm[1][0] =  l_cross[2];
      l_sm[1][1] = 0;
      l_sm[1][2] = -l_cross[0];

      l_sm[2][0] = -l_cross[1];
      l_sm[2][1] =  l_cross[0];
      l_sm[2][2] = 0;

      // get scaling
      TL_T_REAL l_scale = Geom::sprod3( l_v0, l_v1 );
      EDGE_CHECK( std::abs(l_scale + (TL_T_REAL) 1.0) > TOL.MESH );
      l_scale = (TL_T_REAL) 1.0 / ( (TL_T_REAL) 1.0 + l_scale );

      // assemble rotation matrix
      Matrix::matMulB0( 3, 3, 3,
                        l_sm[0], l_sm[0], o_rm[0] );

      for( unsigned short l_d0 = 0; l_d0 < 3; l_d0++ ) {
        for( unsigned short l_d1 = 0; l_d1 < 3; l_d1++ ) {
          o_rm[l_d0][l_d1] *= l_scale;
          o_rm[l_d0][l_d1] += l_sm[l_d0][l_d1];
        }
      }
      o_rm[0][0] += (TL_T_REAL) 1.0;
      o_rm[1][1] += (TL_T_REAL) 1.0;
      o_rm[2][2] += (TL_T_REAL) 1.0;
    }
};

#endif
