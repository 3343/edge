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
 * Shape functions.
 **/

#ifndef SHAPE_HPP
#define SHAPE_HPP

#include <cmath>
#include "constants.hpp"

namespace edge {
  namespace linalg {
    class Shape;
  }
}

class edge::linalg::Shape {
  public:
    /*
     * Evaluates the nodal shape functions for line elements at the given point.
     *
     * Node ordering for our ref. element [0,1]:
     *
     *  0*********1
     *
     * @param i_xi xi-coord.
     * @param i_eta eta-coord.
     * @param o_evalSf will be set to evaluated shape functions.
     *
     * @paramt TL_T_REAL floating point type.
     **/
    template< typename TL_T_REAL >
    static void line( TL_T_REAL i_xi,
                      TL_T_REAL o_evalSf[2] ) {
      o_evalSf[0] = 1 - i_xi;
      o_evalSf[1] = i_xi;
    }

    /*
     * Evaluates the derivatives of the nodal shape functions for line elements at the given point.
     *
     * Node ordering for our ref. element [0,1]:
     *
     *  0*********1
     *
     * @param i_xi xi-coord.
     * @param i_eta eta-coord.
     * @param o_evalSf will be set to evaluated derivatives of the shape functions.
     *
     * @paramt TL_T_REAL floating point type.
     **/
    template< typename TL_T_REAL >
    static void derLine( TL_T_REAL o_evalSf[2] ) {
      o_evalSf[0] = -1;
      o_evalSf[1] =  1;
    }

    /*
     * Evaluates the nodal shape functions of 4-point quadrilaterals at the given point.
     *
     * Node ordering for our ref element [0,1]^2:
     *
     *   3           2
     *     *********
     *     *       *
     * eta *       *
     *     *       *
     *     *********
     *   0     xi    1

     * Reference: Diersch, FEFLOW: Finite Element Modeling of Flow, Mass and Heat Transport
     *            in Porous and Fractured Media, Springer 2014
     *            (H.46) -- modified to match [0,1]^3 reference element
     *
     * @param i_xi xi-coord.
     * @param i_eta eta-coord.
     * @param o_evalSf will be set to evaluated shape functions.
     *        If called for vertex v with xi_v, eta_v, o_evalSf[v] = 1, 0 elsewhere.
     *
     * @paramt TL_T_REAL floating point type.
     **/
    template< typename TL_T_REAL >
    static void quad4( TL_T_REAL i_xi,
                       TL_T_REAL i_eta,
                       TL_T_REAL o_evalSf[4] ) {
      o_evalSf[0] = ( 1.0 - i_xi ) * ( 1.0 - i_eta );
      o_evalSf[1] =   i_xi   *       ( 1.0 - i_eta );
      o_evalSf[2] =   i_xi   *         i_eta;
      o_evalSf[3] = ( 1.0 - i_xi ) *   i_eta;
    }

    /**
     * Evaluates the derivatives of the shape function of 4-point quadrilaterals at the given point.
     *
     * Node ordering for our ref element [0,1]^2: see quad4.
     *
     * @param i_xi xi-coord.
     * @param i_eta eta-coord.
     * @param o_evalDerSf evaluated derivatives of shape functions, [0][*]: d / dxi, [1][*]: d / deta.
     *
     * @paramt TL_T_REAL floating point type.
     **/
    template< typename TL_T_REAL >
    static void derQuad4( TL_T_REAL i_xi,
                          TL_T_REAL i_eta,
                          TL_T_REAL o_evalSfDer[2][4] ) {
      o_evalSfDer[0][0] = -( 1.0 - i_eta );
      o_evalSfDer[0][1] =  ( 1.0 - i_eta );
      o_evalSfDer[0][2] =  i_eta;
      o_evalSfDer[0][3] = -i_eta;

      o_evalSfDer[1][0] = -( 1.0 - i_xi );
      o_evalSfDer[1][1] = -i_xi;
      o_evalSfDer[1][2] =  i_xi;
      o_evalSfDer[1][3] =  ( 1.0 - i_xi );
    }

    /*
     * Evaluates the nodal shape functions of 3-point triangular elements at the given point.
     *
     * Node ordering for our ref element [0,1]x[0,1-xi]:
     *
     *   2
     *     *
     *     * *
     * eta *   *
     *     *     *
     *     *********
     *   0     xi    1
     *
     * Reference: Diersch, FEFLOW: Finite Element Modeling of Flow, Mass and Heat Transport
     *            in Porous and Fractured Media, Springer 2014
     *            (H.11)
     *
     * @param i_xi xi-coord.
     * @param i_eta eta-coord.
     * @param o_evalSf will be set to evaluated shape functions.
     *        If called for vertex v with xi_v, eta_v, o_evalSf[v] = 1, 0 elsewhere.
     *
     * @paramt TL_T_REAL floating point type.
     **/
    template< typename TL_T_REAL >
    static void tria3( TL_T_REAL i_xi,
                       TL_T_REAL i_eta,
                       TL_T_REAL o_evalSf[3] ) {
      o_evalSf[0] = TL_T_REAL(1) - i_xi - i_eta;
      o_evalSf[1] = i_xi;
      o_evalSf[2] = i_eta;
    }

    /**
     * Evaluates the derivatives of the shape functions of 3-point triangular elements at the given point.
     *
     * Node ordering for our ref element [0,1]x[0,1-xi]: see tria3.
     *
     * @param o_evalDerSf evaluated derivatives of shape functions, [0][*]: d / dxi, [1][*]: d / deta.
     *
     * @paramt TL_T_REAL floating point type.
     **/
    template< typename TL_T_REAL >
    static void derTria3( TL_T_REAL o_evalSfDer[2][3] ) {
      o_evalSfDer[0][0] = -1;
      o_evalSfDer[0][1] =  1;
      o_evalSfDer[0][2] =  0;

      o_evalSfDer[1][0] = -1;
      o_evalSfDer[1][1] =  0;
      o_evalSfDer[1][2] =  1;
    }

    /**
     * Evaluates the nodal shape functions of linear, 8-point hexahedral elements at the given point.
     *
     * Node ordering for our ref element [0,1]^3:
     *
     *           7 x*******************x 6
     *            **                  **
     *           * *                 * *
     *          *  *                *  *
     *         *   *               *   *
     *        *    *              *    *
     *     4 x*******************x 5   *
     *       *     *             *     *
     *       *   3 x************ * ****x 2
     *       *    *              *    *
     *       *   *               *   *
     *       |  /                *  *
     *  zeta | / eta             * *
     *       |/                  **
     *       x---****************x
     *     0   xi                 1
     *
     * Reference: Diersch, FEFLOW: Finite Element Modeling of Flow, Mass and Heat Transport
     *            in Porous and Fractured Media, Springer 2014
     *            (H.54) -- modified to match [0,1]^3 reference element
     *
     * @param i_xi xi-coord.
     * @param i_eta eta-coord.
     * @param i_zeta zeta-coord.
     * @param o_evalSf evaluated shape functions at the given coordinates.
     *
     * @paramt TL_T_REAL floating point type.
     **/
    template< typename TL_T_REAL >
    static void hex8Lin( TL_T_REAL i_xi,
                         TL_T_REAL i_eta,
                         TL_T_REAL i_zeta,
                         TL_T_REAL o_evalSf[8] ) {
      o_evalSf[0] = (1 - i_xi) * (1 - i_eta) * (1 - i_zeta);
      o_evalSf[1] = i_xi       * (1 - i_eta) * (1 - i_zeta);
      o_evalSf[2] = i_xi       * i_eta       * (1 - i_zeta);
      o_evalSf[3] = (1 - i_xi) * i_eta       * (1 - i_zeta);

      o_evalSf[4] = (1 - i_xi) * (1 - i_eta) * i_zeta;
      o_evalSf[5] = i_xi       * (1 - i_eta) * i_zeta;
      o_evalSf[6] = i_xi       * i_eta       * i_zeta;
      o_evalSf[7] = (1 - i_xi) * i_eta       * i_zeta;
    }

    /*
     * Evaluates the derivatives of the nodal shape functions for linear,
     * 8-point hexahedral elements at the given point.
     *
     * Node ordering: see hex8Lin.
     *
     * @param i_xi xi-coord.
     * @param i_eta eta-coord.
     * @param i_zeta zeta-coord.
     * @param o_evalDerSf evaluated derivatives of shape functions, [0][*]: d / dxi, [1][*]: d / deta, [2][*] d / dzeta
     *
     * @paramt TL_T_REAL floating point type.
     **/
    template< typename TL_T_REAL >
    static void derHex8Lin( TL_T_REAL i_xi,
                            TL_T_REAL i_eta,
                            TL_T_REAL i_zeta,
                            TL_T_REAL o_evalSf[3][8] ) {
      o_evalSf[0][0] = -1 * (1 - i_eta) * (1 - i_zeta);
      o_evalSf[0][1] =  1 * (1 - i_eta) * (1 - i_zeta);
      o_evalSf[0][2] =  1 * i_eta       * (1 - i_zeta);
      o_evalSf[0][3] = -1 * i_eta       * (1 - i_zeta);

      o_evalSf[0][4] = -1 * (1 - i_eta) * i_zeta;
      o_evalSf[0][5] =  1 * (1 - i_eta) * i_zeta;
      o_evalSf[0][6] =  1 * i_eta       * i_zeta;
      o_evalSf[0][7] = -1 * i_eta       * i_zeta;

      o_evalSf[1][0] = (1 - i_xi) * -1 * (1 - i_zeta);
      o_evalSf[1][1] = i_xi       * -1 * (1 - i_zeta);
      o_evalSf[1][2] = i_xi       *  1 * (1 - i_zeta);
      o_evalSf[1][3] = (1 - i_xi) *  1 * (1 - i_zeta);

      o_evalSf[1][4] = (1 - i_xi) * -1 * i_zeta;
      o_evalSf[1][5] = i_xi       * -1 * i_zeta;
      o_evalSf[1][6] = i_xi       *  1 * i_zeta;
      o_evalSf[1][7] = (1 - i_xi) *  1 * i_zeta;

      o_evalSf[2][0] = (1 - i_xi) * (1 - i_eta) * -1;
      o_evalSf[2][1] = i_xi       * (1 - i_eta) * -1;
      o_evalSf[2][2] = i_xi       * i_eta       * -1;
      o_evalSf[2][3] = (1 - i_xi) * i_eta       * -1;

      o_evalSf[2][4] = (1 - i_xi) * (1 - i_eta) * 1;
      o_evalSf[2][5] = i_xi       * (1 - i_eta) * 1;
      o_evalSf[2][6] = i_xi       * i_eta       * 1;
      o_evalSf[2][7] = (1 - i_xi) * i_eta       * 1;
    }

    /**
     * Evaluates the nodal shape functions of 4-point tetrahedrons at the given point.
     *
     * Node ordering for our ref-element:
     * 0 =< xi1
     * 0 =< xi2
     * 0 =< xi3
     * xi1 + xi2 + xi3 <= 1
     *
     *                 zeta
     *                  3: x
     *                   *
     *                  *    *
     *                 *   *
     *                *         *
     *               *
     *              *       *      *
     *             *
     *            *                    *
     *           * origin(0): x
     *          *           *     *       *
     *         *        *              *
     *        *      *                       x
     *       *    *               *          2: eta
     *      *  *     *
     * 1: x
     * xi
     *
     * Reference: Diersch, FEFLOW: Finite Element Modeling of Flow, Mass and Heat Transport
     *            in Porous and Fractured Media, Springer 2014
     *            (H.27)
     *
     * @param i_xi xi-coord.
     * @param i_eta eta-coord.
     * @param i_zeta zeta-coord.
     * @param o_evalSf evaluated shape functions at the given coordinate.
     *
     * @paramt TL_T_REAL floating point type.
     **/
    template< typename TL_T_REAL >
    static void tet4( TL_T_REAL i_xi,
                      TL_T_REAL i_eta,
                      TL_T_REAL i_zeta,
                      TL_T_REAL o_evalSf[4] ) {
      o_evalSf[0] = 1 - i_xi - i_eta - i_zeta;
      o_evalSf[1] = i_xi;
      o_evalSf[2] = i_eta;
      o_evalSf[3] = i_zeta;
    }

    /**
     * Evaluates the derivatives of the shape functions of 4-point tetrahedral elements.
     *
     * Node ordering for our ref-element: see tet4.
     *
     * @param o_evalDerSf evaluated derivatives of shape functions, [0][*]: d / dxi, [1][*]: d / deta, [2][*] d / dzeta
     *
     * @paramt TL_T_REAL floating point type.
     **/
    template< typename TL_T_REAL >
    static void derTet4( TL_T_REAL o_evalSfDer[3][4] ) {
      o_evalSfDer[0][0] = -1;
      o_evalSfDer[0][1] =  1;
      o_evalSfDer[0][2] =  0;
      o_evalSfDer[0][3] =  0;

      o_evalSfDer[1][0] = -1;
      o_evalSfDer[1][1] =  0;
      o_evalSfDer[1][2] =  1;
      o_evalSfDer[1][3] =  0;

      o_evalSfDer[2][0] = -1;
      o_evalSfDer[2][1] =  0;
      o_evalSfDer[2][2] =  0;
      o_evalSfDer[2][3] =  1;
    }
};

#endif
