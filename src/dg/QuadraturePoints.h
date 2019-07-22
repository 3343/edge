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
 * Computes quadrature points.
 **/
#ifndef EDGE_DG_QUADRATURE_POINTS_H
#define EDGE_DG_QUADRATURE_POINTS_H

#include "constants.hpp"
#include "linalg/Mappings.hpp"
#include "linalg/Geom.hpp"
#include <vector>
#include <cmath>
#include <cassert>
#include "monitor/instrument.hpp"
#include "submodules/FastGL/fastgl.h"

namespace edge {
  namespace dg {
    class QuadraturePoints;
  }
}

class edge::dg::QuadraturePoints {
  private:
    /**
     * Gets the Gauss-Legendre quadrature points and weights.
     *
     * @param i_nPts number of requested points.
     * @param o_pts will be set to points in [0,1].
     * @param o_weights will be set to weights.
     *
     * @paramt TL_T_REAL_MESH floating point type of the points.
     * @paramt TL_T_REAL_COMP floating point type of the weights.
     **/
    template< typename TL_T_REAL_MESH, typename TL_T_REAL_COMP >
    static void gaussLegendre( unsigned short   i_nPts,
                               TL_T_REAL_MESH  *o_pts,
                               TL_T_REAL_COMP  *o_weights  ) {
      EDGE_CHECK_GT( i_nPts, 0 );

      // temporary storage of a quadrature pair
      fastgl::QuadPair l_tmp;

      // iterate over nodes
      for( unsigned short l_qp = 0; l_qp < i_nPts; l_qp++ ) {
        // get quad-pair in [-1,1]
        l_tmp = fastgl::GLPair( i_nPts, i_nPts-l_qp );

        // copy the result
        o_pts[l_qp]     = l_tmp.x();
        o_weights[l_qp] = l_tmp.weight;

        // scale the result
        o_pts[l_qp] += (TL_T_REAL_MESH) 1.0;
        o_pts[l_qp] *= (TL_T_REAL_MESH) 0.5;

        o_weights[l_qp] *= 0.5;
      }
    }

  public:
    /**
     * Gets the Gaussian quadrature points for a line element.
     * For now only the x-dimension is considered, don't use if vertex coordinates in y- or z- are non-zero.
     * Remark: The weights are already scaled by the determinant of the transformation.
     *
     * @param i_order requested order of the integration rule.
     * @param i_veCoords coordinates of the vertices.
     * @param o_ptsX will be set to x-coords of the quadrature points.
     * @param o_ptsY will be set to y-coords of the quadrature points.
     * @param o_ptsZ will be set to z-coords of the quadrature points.
     * @param o_ptsWeights will be set the weights of quadrature points.
'    *
     * @paramt TL_T_REAL_MESH floating point type of the points.
     * @paramt TL_T_REAL_COMP floating point type of the weights.
     **/
    template< typename TL_T_REAL_MESH, typename TL_T_REAL_COMP >
    static void getQptsLine(       unsigned int                   i_order,
                             const TL_T_REAL_MESH                 i_veCoords[3][2],
                                   std::vector< TL_T_REAL_MESH > &o_ptsX,
                                   std::vector< TL_T_REAL_MESH > &o_ptsY,
                                   std::vector< TL_T_REAL_MESH > &o_ptsZ,
                                   std::vector< TL_T_REAL_COMP > &o_ptsWeights ) {
      // resize the vectors
      o_ptsX.resize( i_order ); o_ptsY.resize( i_order ); o_ptsZ.resize( i_order );
      o_ptsWeights.resize( i_order );

      // quad points in one dimension
      std::vector< TL_T_REAL_MESH > l_pts1D;
      l_pts1D.resize( i_order );
      gaussLegendre( i_order, &l_pts1D[0], &o_ptsWeights[0] );

      // reference line, remark this is local to this function only
      TL_T_REAL_MESH l_refLine[3][2] = { {0,1}, {0,0}, {0,0} };

      // derive scaling for the weights, which equals the length of the line
      TL_T_REAL_COMP l_scale = linalg::Mappings::distPts( i_veCoords );

      // map the quad points to the requested 3D line
      for( unsigned int l_qp = 0; l_qp < i_order; l_qp++ ) {
        // coordinates of quad point w.r.t. to ref line
        TL_T_REAL_MESH l_coordsRef[3] = {0, 0, 0};
        l_coordsRef[0] = l_pts1D[l_qp];

        // coordinates of quad point w.r.t. to requrested line
        TL_T_REAL_MESH l_coordsReq[3];

        linalg::Mappings::map2dLine( l_refLine,
                                     i_veCoords,
                                     l_coordsRef,
                                     l_coordsReq );

        o_ptsX[l_qp] = l_coordsReq[0];
        o_ptsY[l_qp] = l_coordsReq[1];
        o_ptsZ[l_qp] = l_coordsReq[2];

        o_ptsWeights[l_qp] *= l_scale;
      }
    }

    /**
     * Gets the gaussian intergration points in phyiscal coords for quads.
     * Remark: The weights already include the determinant of the transformation to physical space.
     *         (-> Integration of the constant basis gives the volume of the element in physical space).
     * Remark: We assume that the quad is aligned to two of the tree possible dimensions.
     *         By doing so we take a shortcut not going over shape functions.
     *
     * @param i_order order of the 1D gaussian integration rule.
     * @param i_veCoords physical coords of the quad's vertices. A counter-clockwise order of the vertices is assumed.
     * @param o_ptsX will be set to x-coords of the quadrature points.
     * @param o_ptsY will be set to y-coords of the quadrature points.
     * @param o_ptsZ will be set to z-coords of the quadrature points.
     * @param o_ptsWeights will be set the weights of quadrature points.
     *
     * @paramt TL_T_REAL_MESH floating point type of the points.
     * @paramt TL_T_REAL_COMP floating point type of the weights.
     **/
    template < typename TL_T_REAL_MESH, typename TL_T_REAL_COMP >
    static void getQptsQuad4(       unsigned int                   i_order,
                              const TL_T_REAL_MESH                 i_veCoords[3][4],
                                    std::vector< TL_T_REAL_MESH > &o_ptsX,
                                    std::vector< TL_T_REAL_MESH > &o_ptsY,
                                    std::vector< TL_T_REAL_MESH > &o_ptsZ,
                                    std::vector< TL_T_REAL_COMP > &o_ptsWeight ) {
      // get 1D quad points
      std::vector< TL_T_REAL_MESH > l_pts1D;
      std::vector< TL_T_REAL_COMP > l_weights1D;
      l_pts1D.resize( i_order ); l_weights1D.resize( i_order );
      gaussLegendre( i_order, &l_pts1D[0], &l_weights1D[0] );

      // get the min and max-positions of the quad
      TL_T_REAL_MESH l_minMax[3][2];
      // init
      for( unsigned short l_dim = 0; l_dim < 3; l_dim++ ) {
        l_minMax[l_dim][0] = l_minMax[l_dim][1] = i_veCoords[l_dim][0];
      }

      for( unsigned short l_dim = 0; l_dim < 3; l_dim++ ) {
        for( unsigned short l_pt = 1; l_pt < 4; l_pt++ ) {
          l_minMax[l_dim][0] = std::min( l_minMax[l_dim][0], i_veCoords[l_dim][l_pt] );
          l_minMax[l_dim][1] = std::max( l_minMax[l_dim][1], i_veCoords[l_dim][l_pt] );
        }
      }

      // find the free dimensions
      unsigned int l_dims[2];
      unsigned int l_localDim = 0;
      unsigned int l_fixedDim = std::numeric_limits<unsigned int>::max();
      for( unsigned short l_dim = 0; l_dim < 3; l_dim++ ) {
        if( l_minMax[l_dim][1] - l_minMax[l_dim][0] > TOL.MESH ) {
          l_dims[l_localDim] = l_dim;
          l_localDim++;
        }
        else l_fixedDim = l_dim;
      }
      assert( l_localDim == 2);

      // get the integration points
      for( unsigned int l_p1 = 0; l_p1 < i_order; l_p1++ ) {
        for( unsigned int l_p2 = 0; l_p2 < i_order; l_p2++ ) {
          // set free dimensions
          if( l_dims[0] == 0 ) {
            o_ptsX.push_back( (l_minMax[0][1] - l_minMax[0][0]) * l_pts1D[l_p2] + l_minMax[0][0] );
          }
          else if( l_dims[0] == 1 ) {
            o_ptsY.push_back( (l_minMax[1][1] - l_minMax[1][0]) * l_pts1D[l_p2] + l_minMax[1][0] );
          }
          else assert( false );

          if( l_dims[1] == 1 ) {
            o_ptsY.push_back( (l_minMax[1][1] - l_minMax[1][0]) * l_pts1D[l_p1] + l_minMax[1][0] );
          }
          else if( l_dims[1] == 2 ) {
            o_ptsZ.push_back( (l_minMax[2][1] - l_minMax[2][0]) * l_pts1D[l_p1] + l_minMax[2][0] );
          }
          else assert( false );

          // set fixed dimension
          if(      l_fixedDim == 0 ) o_ptsX.push_back( l_minMax[0][0] );
          else if( l_fixedDim == 1 ) o_ptsY.push_back( l_minMax[1][0] );
          else if( l_fixedDim == 2 ) o_ptsZ.push_back( l_minMax[2][0] );
          else assert( false );

          // determine determinant, which is simply the area of the face
          TL_T_REAL_MESH l_absDet =   ( l_minMax[ l_dims[0] ][ 1 ] - l_minMax[ l_dims[0] ][ 0 ] )
                                    * ( l_minMax[ l_dims[1] ][ 1 ] - l_minMax[ l_dims[1] ][ 0 ] );
          assert( l_absDet > TOL.MESH );

          o_ptsWeight.push_back( l_weights1D[l_p1] * l_weights1D[l_p2] * l_absDet );
        }
      }
    }

    /**
     * Gets quadrature points and weights of the given tet3-element in physical space via collapsed coords.
     *
     * Remark 1: The weights already include the determinant of the transformation to physical space.
     *          (-> Integration of the constant basis gives the volume of the element in physical space).
     *
     * Remark 2: We derive the quadpoints by a trafo to the collapsed coordinate system; this artificially
     *           increases the number of quad points, plus better options, e.g. being symmetric, might exist.
     *           However, we are able to go towards arbitrary oders (in theory: machine precision..).
     *
     * Remark 3: We completely unroll the quadrature. More efficient schemes,
     *           utilizing the tensor structure, are possible.
     *
     * Reference: Karniadakis, George, and Spencer Sherwin.
     *            Spectral/hp element methods for computational fluid dynamics.
     *            Oxford University Press, 2013.
     *            Chapter 4.1.1.2 -- modified to match our xi1, xi2 > 0, xi1^2 + xi^2 = 1 reference triangle.
     *
     * @param i_order order of the Gaussian integration.
     * @param i_veCoords physical coords of the triangles' vertices.
     * @param o_ptsX will be set to x-coords of the gaussian integration points.
     * @param o_ptsY will be set to y-coords of the gaussian integration points.
     * @param o_ptsZ will be set to z-coords of the gaussian integration points.
     * @param o_weights will be set to the weights of teh gaussion integration points.
     *
     * @paramt TL_T_REAL_MESH floating point type of the points.
     * @paramt TL_T_REAL_COMP floating point type of the weights.
     **/
    template < typename TL_T_REAL_MESH, typename TL_T_REAL_COMP >
    static void getQPtsTria3(       unsigned short                 i_order,
                              const TL_T_REAL_MESH                 i_veCoords[3][3],
                                    std::vector< TL_T_REAL_MESH > &o_ptsX,
                                    std::vector< TL_T_REAL_MESH > &o_ptsY,
                                    std::vector< TL_T_REAL_MESH > &o_ptsZ,
                                    std::vector< TL_T_REAL_COMP > &o_weights ) {
      // get 1D quad points
      std::vector< TL_T_REAL_MESH > l_pts1D;
      std::vector< TL_T_REAL_COMP > l_weights1D;
      l_pts1D.resize( i_order ); l_weights1D.resize( i_order );
      gaussLegendre( i_order, &l_pts1D[0], &l_weights1D[0] );

      // scale to [-1, 1], the respective mappings to [0, 1] are part of the collapsed coord handling
      for( unsigned short l_qp = 0; l_qp < i_order; l_qp++ ) {
        l_pts1D[l_qp] *= (TL_T_REAL_MESH) 2.0;
        l_pts1D[l_qp] -= (TL_T_REAL_MESH) 1.0;
        l_weights1D[l_qp] *= (TL_T_REAL_COMP) 2.0;
      }


      // get number of unrolled, 2D quadpoints
      unsigned short l_nQpts = i_order*i_order;

      // create local and resize output vectors
      std::vector< TL_T_REAL_MESH > l_ptsXi, l_ptsEta, l_weights;
      l_ptsXi.resize(    l_nQpts );
      l_ptsEta.resize(   l_nQpts );
      l_weights.resize(  l_nQpts );

      o_ptsX.resize(     l_nQpts );
      o_ptsY.resize(     l_nQpts );
      o_ptsZ.resize(     l_nQpts );
      o_weights.resize(  l_nQpts );

      // iterate dimension-wise over quadpoints
      for( unsigned short l_dim1 = 0; l_dim1 < i_order; l_dim1++ ) {
        for( unsigned short l_dim2 = 0; l_dim2 < i_order; l_dim2++ ) {
          // derive index of the quad point
          unsigned short l_qIn = l_dim1*i_order + l_dim2;

          // get pts in tet
          l_ptsXi[l_qIn]  = l_pts1D[l_dim1] + 1;
          l_ptsXi[l_qIn] *= l_pts1D[l_dim2] - 1;
          l_ptsXi[l_qIn] /= -4;

          l_ptsEta[l_qIn]  = l_pts1D[l_dim2] + 1;
          l_ptsEta[l_qIn] /= 2;

          l_weights[l_qIn] = l_weights1D[l_dim1] * l_weights1D[l_dim2];

          // scale weights with jacobian of the trafo
          l_weights[l_qIn] *= (l_pts1D[l_dim2]-1);
          l_weights[l_qIn] /= -8;
        }
      }

      for( unsigned int l_qp = 0; l_qp < l_nQpts; l_qp++ ) {
        // get quadrature point in physical coords
        TL_T_REAL_MESH l_pt[2];
        linalg::Mappings::mapRefToPhyTria3( l_ptsXi[l_qp],
                                            l_ptsEta[l_qp],
                                            i_veCoords,
                                            l_pt );

        o_ptsX[l_qp] = l_pt[0];
        o_ptsY[l_qp] = l_pt[1];
        o_ptsZ[l_qp] = 0; // TODO: Add third dimension

        // determine weights
        TL_T_REAL_MESH l_jac[2][2];
        linalg::Mappings::evalJac( TRIA3, i_veCoords[0], l_jac[0] );
        TL_T_REAL_MESH l_det = linalg::Matrix::det( l_jac );
        assert( l_det > TOL.MESH );

        o_weights[l_qp] = l_weights[l_qp] * l_det;
      }
    }

    /**
     * Gets the gaussian intergration points in phyiscal coords for rectangular, 8-points hexes.
     *
     * Remark: The weights already include the determinant of the transformation to physical space.
     *         (-> Integration of the constant basis gives the volume of the element in physical space).
     *
     * @param i_order order of the 1D gaussian integration rule.
     * @param i_veCoords physical coords of the quad's vertices.
     * @param o_ptsX will be set to x-coords of the quadrature points.
     * @param o_ptsY will be set to y-coords of the quadrature points.
     * @param o_ptsZ will be set to z-coords of the quadrature points.
     * @param o_ptsWeights will be set the weights of quadrature points.
     *
     * @paramt TL_T_REAL_MESH floating point type of the points.
     * @paramt TL_T_REAL_COMP floating point type of the weights.
     **/
    template < typename TL_T_REAL_MESH, typename TL_T_REAL_COMP >
    static void getQPtsHex8r(       unsigned int                   i_order,
                              const TL_T_REAL_MESH                 i_veCoords[3][8],
                                    std::vector< TL_T_REAL_MESH > &o_ptsX,
                                    std::vector< TL_T_REAL_MESH > &o_ptsY,
                                    std::vector< TL_T_REAL_MESH > &o_ptsZ,
                                    std::vector< TL_T_REAL_COMP > &o_ptsWeight ) {
      // get 1D quad points
      std::vector< TL_T_REAL_MESH > l_pts1D;
      std::vector< TL_T_REAL_COMP > l_weights1D;
      l_pts1D.resize( i_order ); l_weights1D.resize( i_order );
      gaussLegendre( i_order, &l_pts1D[0], &l_weights1D[0] );

      for( unsigned int l_p1 = 0; l_p1 < i_order; l_p1++ ) {
        for( unsigned int l_p2 = 0; l_p2 < i_order; l_p2++ ) {
          for( unsigned int l_p3 = 0; l_p3 < i_order; l_p3++ ) {
            // get quad point in physical coords
            TL_T_REAL_MESH l_pt[3];
            linalg::Mappings::mapRefToPhyHex8r( l_pts1D[l_p3],
                                                l_pts1D[l_p2],
                                                l_pts1D[l_p1],
                                                i_veCoords,
                                                l_pt );
            o_ptsX.push_back( l_pt[0] );
            o_ptsY.push_back( l_pt[1] );
            o_ptsZ.push_back( l_pt[2] );

            // determine weights
            TL_T_REAL_MESH l_jac[3][3];
            linalg::Mappings::evalJac( HEX8R, i_veCoords[0], l_jac[0] );
            TL_T_REAL_MESH l_det = linalg::Matrix::det( l_jac );
            assert( l_det > TOL.MESH );

            o_ptsWeight.push_back( l_weights1D[l_p1] * l_weights1D[l_p2] * l_weights1D[l_p3] * l_det );
          }
        }
      }
    }

    /**
     * Gets the quadrature points and weights of a tetrahedron in physical space.
     *
     * Remark 1: The weights already include the determinant of the transformation to physical space.
     *          (-> Integration of the constant basis gives the volume of the element in physical space).
     *
     * Remark 2: We derive the quadpoints by a trafo to the collapsed coordinate system; this artificially
     *           increases the number of quad points, plus better options, e.g. being symmetric, might exist.
     *           However, we are able to go towards arbitrary oders (in theory: machine precision..).
     *
     * Remark 3: We completely unroll the quadrature. More efficient schemes,
     *           utilizing the tensor structure, are possible.
     *
     * Reference: Karniadakis, George, and Spencer Sherwin.
     *            Spectral/hp element methods for computational fluid dynamics.
     *            Oxford University Press, 2013.
     *            Chapter 4.1.1.2 -- modified to match our xi1, xi2, xi3 > 0, xi1^2 + xi^2 + xi^3 = 1
     *                               reference tetrahedron.
     *
     * @param i_order order of the quadrature rule.
     * @param i_veCoords physical coords of the tetrahedrons's vertices.
     * @param o_ptsX will be set to x-coords of the quadraturepoints.
     * @param o_ptsY will be set to y-coords of the quadrature points.
     * @param o_ptsZ will be set to z-coords of the quadrature points.
     * @param o_weights will be set to the weights of the quadrature points.
     *
     * @paramt TL_T_REAL_MESH floating point type of the points.
     * @paramt TL_T_REAL_COMP floating point type of the weights.
     **/
    template< typename TL_T_REAL_MESH, typename TL_T_REAL_COMP >
    static void getQPtsTet4(       unsigned int                   i_order,
                             const TL_T_REAL_MESH                 i_veCoords[3][4],
                                   std::vector< TL_T_REAL_MESH > &o_ptsX,
                                   std::vector< TL_T_REAL_MESH > &o_ptsY,
                                   std::vector< TL_T_REAL_MESH > &o_ptsZ,
                                   std::vector< TL_T_REAL_COMP > &o_weights ) {
      // get 1D quad points
      std::vector< TL_T_REAL_MESH > l_pts1D;
      std::vector< TL_T_REAL_COMP > l_weights1D;
      l_pts1D.resize( i_order ); l_weights1D.resize( i_order );
      gaussLegendre( i_order, &l_pts1D[0], &l_weights1D[0] );

      // scale to [-1, 1], the respective mappings to [0, 1] are part of the collapsed coord handling
      for( unsigned short l_qp = 0; l_qp < i_order; l_qp++ ) {
        l_pts1D[l_qp] *= (TL_T_REAL_MESH) 2.0;
        l_pts1D[l_qp] -= (TL_T_REAL_MESH) 1.0;
        l_weights1D[l_qp] *= (TL_T_REAL_COMP) 2.0;
      }

      // get unrolled, 3D quad points via collapsed coordinate system
      unsigned short l_nQpts = i_order * i_order * i_order;

      // resize output vectors
      o_ptsX.resize(l_nQpts);
      o_ptsY.resize(l_nQpts);
      o_ptsZ.resize(l_nQpts);
      o_weights.resize(l_nQpts);

      // get quad points over reference tetrahedron
      std::vector< TL_T_REAL_MESH > l_ptsXi, l_ptsEta, l_ptsZeta;
      std::vector< TL_T_REAL_COMP > l_wgts;
      l_ptsXi.resize(   l_nQpts );
      l_ptsEta.resize(  l_nQpts );
      l_ptsZeta.resize( l_nQpts );
      l_wgts.resize(    l_nQpts );

      // iterate dimension-wise over quadpoints
      for( unsigned short l_dim1 = 0; l_dim1 < i_order; l_dim1++ ) {
        for( unsigned short l_dim2 = 0; l_dim2 < i_order; l_dim2++ ) {
          for( unsigned short l_dim3 = 0; l_dim3 < i_order; l_dim3++ ) {
            // derive index of the quad point
            unsigned short l_qIn = l_dim1*i_order*i_order + l_dim2*i_order + l_dim3;

            // get pts in tet
            l_ptsXi[l_qIn]    = l_pts1D[l_dim1] + 1;
            l_ptsXi[l_qIn]   *= l_pts1D[l_dim2] - 1;
            l_ptsXi[l_qIn]   *= l_pts1D[l_dim3] - 1;
            l_ptsXi[l_qIn]   /= 8;

            l_ptsEta[l_qIn]   = l_pts1D[l_dim2] + 1;
            l_ptsEta[l_qIn]  *= l_pts1D[l_dim3] - 1;
            l_ptsEta[l_qIn]  /= -4;

            l_ptsZeta[l_qIn]  = l_pts1D[l_dim3] + 1;
            l_ptsZeta[l_qIn] /= 2;

            l_wgts[l_qIn] = l_weights1D[l_dim1] * l_weights1D[l_dim2] * l_weights1D[l_dim3];

            // scale weights with jacobian of the trafo
            l_wgts[l_qIn] *= (l_pts1D[l_dim2]-1);
            l_wgts[l_qIn] *= (l_pts1D[l_dim3]-1);
            l_wgts[l_qIn] *= (l_pts1D[l_dim3]-1);
            l_wgts[l_qIn] /= -64;
          }
        }
      }

      for( unsigned int l_qp = 0; l_qp < l_nQpts; l_qp++ ) {
        // get quadrature point in physical coords
        TL_T_REAL_MESH l_pt[3];
        linalg::Mappings::mapRefToPhyTet4( l_ptsXi[l_qp],
                                           l_ptsEta[l_qp],
                                           l_ptsZeta[l_qp],
                                           i_veCoords,
                                           l_pt );

        o_ptsX[l_qp] = l_pt[0];
        o_ptsY[l_qp] = l_pt[1];
        o_ptsZ[l_qp] = l_pt[2];

        // determine adjusted weights (physical coords)
        TL_T_REAL_MESH l_jac[3][3];
        linalg::Mappings::evalJac( TET4, i_veCoords[0], l_jac[0] );
        TL_T_REAL_MESH l_det = linalg::Matrix::det( l_jac );
        assert( l_det > TOL.MESH );

        o_weights[l_qp] = l_wgts[l_qp] * l_det;
      }
    }

  public:
    /**
     * Gets the quadrature points for the given entity type, order of the quadrature rule and vertices.
     * Remark: The weights already include the determinant of the transformation to physical space.
     *         (-> Integration of the constant basis gives the volume of the element in physical space).
     *
     * @param i_entType requested entity type.
     * @param i_order requested order of the integration rule.
     * @param i_veCoords coordinates of the vertices.
     * @param o_ptsX will be set to x-coords of the quadrature points.
     * @param o_ptsY will be set to y-coords of the quadrature points.
     * @param o_ptsZ will be set to z-coords of the quadrature points.
     * @param o_weights will be set the weights of quadrature points.
     *
     * @paramt TL_T_REAL_MESH floating point type of the points.
     * @paramt TL_T_REAL_COMP floating point type of the weights.
     **/
    template < typename TL_T_REAL_MESH, typename TL_T_REAL_COMP >
    static void getQpts(       t_entityType                   i_entType,
                               unsigned int                   i_order,
                         const TL_T_REAL_MESH                *i_veCoords,
                               std::vector< TL_T_REAL_MESH > &o_ptsX,
                               std::vector< TL_T_REAL_MESH > &o_ptsY,
                               std::vector< TL_T_REAL_MESH > &o_ptsZ,
                               std::vector< TL_T_REAL_COMP > &o_weights ) {
      PP_INSTR_FUN("get_qps")

      if( i_entType == LINE ) {
        getQptsLine(                         i_order,
                     (TL_T_REAL_MESH (*)[2]) i_veCoords,
                                             o_ptsX,
                                             o_ptsY,
                                             o_ptsZ,
                                             o_weights );
      }
      else if( i_entType == QUAD4R ) {
        getQptsQuad4(                         i_order,
                      (TL_T_REAL_MESH (*)[4]) i_veCoords,
                                              o_ptsX,
                                              o_ptsY,
                                              o_ptsZ,
                                              o_weights );
      }
      else if( i_entType == TRIA3 ) {
        getQPtsTria3(                         i_order,
                      (TL_T_REAL_MESH (*)[3]) i_veCoords,
                                              o_ptsX,
                                              o_ptsY,
                                              o_ptsZ,
                                              o_weights );
      }
      else if( i_entType == HEX8R ) {
        getQPtsHex8r(                         i_order,
                      (TL_T_REAL_MESH (*)[8]) i_veCoords,
                                              o_ptsX,
                                              o_ptsY,
                                              o_ptsZ,
                                              o_weights );
      }
      else if( i_entType == TET4 ) {
        getQPtsTet4(                         i_order,
                     (TL_T_REAL_MESH (*)[4]) i_veCoords,
                                             o_ptsX,
                                             o_ptsY,
                                             o_ptsZ,
                                             o_weights );
      }
      else assert(false);
    }

    /**
     * Gets the quadrature points for the given entity type, order of the quadrature rule and vertices.
     * Remark: The weights already include the determinant of the transformation to physical space.
     *         (-> Integration of the constant basis gives the volume of the element in physical space).
     *
     * @param i_entType requested entity type.
     * @param i_order requested order of the integration rule.
     * @param i_veCoords coordinates of the vertices.
     * @param o_ptsX will be set to x-coords of the quadrature points.
     * @param o_ptsY will be set to y-coords of the quadrature points.
     * @param o_ptsZ will be set to z-coords of the quadrature points.
     * @param o_weights will be set the weights of quadrature points.
     *
     * @paramt TL_T_REAL_MESH floating point type of the points.
     * @paramt TL_T_REAL_COMP floating point type of the weights.
     **/
    template< int nVeCoords, typename TL_T_REAL_MESH, typename TL_T_REAL_COMP >
    static void getQpts(       t_entityType                   i_entType,
                               unsigned int                   i_order,
                         const TL_T_REAL_MESH                 i_veCoords[3][nVeCoords],
                               std::vector< TL_T_REAL_MESH > &o_ptsX,
                               std::vector< TL_T_REAL_MESH > &o_ptsY,
                               std::vector< TL_T_REAL_MESH > &o_ptsZ,
                               std::vector< TL_T_REAL_COMP > &o_weights ) {
      // get back some type-safety at runtime
      assert( nVeCoords == C_ENT[i_entType].N_VERTICES );

      // call 1D-array version
      getQpts(                   i_entType,
                                 i_order,
               (TL_T_REAL_MESH*) i_veCoords,
                                 o_ptsX,
                                 o_ptsY,
                                 o_ptsZ,
                                 o_weights );
    }
};
#endif
