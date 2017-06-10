/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2017, Regents of the University of California
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
 * Support for evals through quadrature rules.
 **/

#ifndef QUADRATURE_EVAL_HPP
#define QUADRATURE_EVAL_HPP

#include "constants.hpp"
#include "io/logging.h"
#include "QuadraturePoints.h"
#include "linalg/Mappings.hpp"
#include "dg/Basis.h"

namespace edge {
  namespace dg {
    template <unsigned short TL_ORDER>
    class QuadratureEvalBase;

    template <t_entityType TL_TYPE_ELEMENT, unsigned short TL_ORDER, unsigned short TL_N_CRUNS=1>
    class QuadratureEval;

    template <unsigned short TL_ORDER, unsigned short TL_N_CRUNS>
    class QuadratureEval< LINE, TL_ORDER, TL_N_CRUNS >;
  }
}

template <unsigned short TL_ORDER>
class edge::dg::QuadratureEvalBase {
  public:
    /**
     * Gets the 1D quad points and weights for the line reference element.
     *
     * @param o_pts will be set to location of quad points.
     * @param o_weights will be set to weights of quad points.
     *
     * @paramt TL_T_REAL_MESH floating point type of the points.
     * @paramt TL_T_REAL_COMP floating point type of the weights.
     **/
    template< typename TL_T_REAL_MESH, typename TL_T_REAL_COMP >
    static void line(  TL_T_REAL_MESH o_pts[ CE_N_ELEMENT_MODES( LINE, TL_ORDER ) ],
                       TL_T_REAL_COMP o_weights[ CE_N_ELEMENT_MODES( LINE, TL_ORDER ) ] ) {
      // get 1D quadrature points and weights
      std::vector< TL_T_REAL_MESH > l_qps[3];
      std::vector< TL_T_REAL_COMP > l_weights;

      dg::QuadraturePoints::getQpts( LINE,
                                     TL_ORDER,
                                     C_REF_ELEMENT.VE.ENT[LINE],
                                     l_qps[0], l_qps[1], l_qps[2], l_weights );

      EDGE_CHECK( l_weights.size() == CE_N_ELEMENT_MODES( LINE, TL_ORDER ) );

      for( unsigned short l_qp = 0; l_qp < CE_N_ELEMENT_MODES( LINE, TL_ORDER ); l_qp++ ) {
         o_pts[l_qp] = l_qps[0][l_qp];
         o_weights[l_qp] = l_weights[l_qp];
      }
    }
};

template <t_entityType TL_TYPE_ELEMENT, unsigned short TL_ORDER, unsigned short TL_N_CRUNS>
class edge::dg::QuadratureEval: public QuadratureEvalBase<TL_ORDER> {
  public:
    /**
     * Gets the face's quadrature points.
     *
     * Remark 1: Collapsed coordinates for tetrahedral elements are used,
     *           resulting in more quadrature points than necessary for the given order of convergence.
     *
     * Remark 2: The weights are based on the face type's reference element.
     *           The sum gives the volume of this element.
     *           The face type is given implicitly by the template parameter TL_TYPE_ELEMENT.
     *
     * @param o_pts will be set to quadrature points. [*][][]: flux options, [][*][]: quad point [][][*]: dim.
     * @param o_weights will be set to weights of the quad points.
     * @param o_basisEval will be set to evaluated basis functiosn (volume/element) at the quad points.
     *
     * @paramt TL_T_REAL_MESH floating point type of the points.
     * @paramt TL_T_REAL_COMP floating point type of the basis.
     **/
#ifndef __INTEL_COMPILER
    template < typename TL_T_REAL_MESH, typename TL_T_REAL_COMP >
    static void faces( TL_T_REAL_MESH o_pts[ (CE_N_FACE_VERTEX_OPTS(TL_TYPE_ELEMENT)+1) * C_ENT[TL_TYPE_ELEMENT].N_FACES ]
                                           [  CE_N_FACE_QUAD_POINTS( TL_TYPE_ELEMENT, TL_ORDER )                         ]
                                           [  C_ENT[TL_TYPE_ELEMENT].N_DIM                                               ],
                       TL_T_REAL_COMP o_weights[ CE_N_FACE_QUAD_POINTS( TL_TYPE_ELEMENT, TL_ORDER ) ],
                       TL_T_REAL_COMP o_basisEval[ (CE_N_FACE_VERTEX_OPTS(TL_TYPE_ELEMENT)+1) * C_ENT[TL_TYPE_ELEMENT].N_FACES ]
                                                 [  CE_N_FACE_QUAD_POINTS( TL_TYPE_ELEMENT, TL_ORDER )                         ]
                                                 [  CE_N_ELEMENT_MODES( TL_TYPE_ELEMENT, TL_ORDER )                            ] ) {
#else
    // TODO: Hard-coded work around for Intel compiler.
    //       Fails to compile the template above; icpc version 17.0.3 (gcc version 4.8.5 compatibility)
    static void faces( real_mesh o_pts[ (CE_N_FACE_VERTEX_OPTS(TL_TYPE_ELEMENT)+1) * C_ENT[TL_TYPE_ELEMENT].N_FACES ]
                                      [  CE_N_FACE_QUAD_POINTS( TL_TYPE_ELEMENT, TL_ORDER )                         ]
                                      [  C_ENT[TL_TYPE_ELEMENT].N_DIM                                               ],
                       real_base o_weights[ CE_N_FACE_QUAD_POINTS( TL_TYPE_ELEMENT, TL_ORDER ) ],
                       real_base o_basisEval[ (CE_N_FACE_VERTEX_OPTS(TL_TYPE_ELEMENT)+1) * C_ENT[TL_TYPE_ELEMENT].N_FACES ]
                                            [  CE_N_FACE_QUAD_POINTS( TL_TYPE_ELEMENT, TL_ORDER )                         ]
                                            [  CE_N_ELEMENT_MODES( TL_TYPE_ELEMENT, TL_ORDER )                            ] ) {
      typedef real_mesh TL_T_REAL_MESH;
      typedef real_base TL_T_REAL_COMP;
#endif
      // get quad points in face local coords
      std::vector< TL_T_REAL_MESH > l_faPtsChi[3];
      std::vector< TL_T_REAL_COMP > l_weights;

      dg::QuadraturePoints::getQpts( C_ENT[TL_TYPE_ELEMENT].TYPE_FACES,
                                     TL_ORDER,
                                     C_REF_ELEMENT.VE.ENT[ C_ENT[TL_TYPE_ELEMENT].TYPE_FACES ],
                                     l_faPtsChi[0], l_faPtsChi[1], l_faPtsChi[2], l_weights );

      // double-check number of faces' quad points
      EDGE_CHECK( l_weights.size() == CE_N_FACE_QUAD_POINTS( TL_TYPE_ELEMENT, TL_ORDER ) );

      // iterate over faces and set volume cords of quad points
      for( unsigned short l_fa = 0; l_fa < C_ENT[TL_TYPE_ELEMENT].N_FACES; l_fa++ ) {
        for( unsigned int l_qp = 0; l_qp < CE_N_FACE_QUAD_POINTS( TL_TYPE_ELEMENT, TL_ORDER ); l_qp++ ) {
          TL_T_REAL_MESH l_faPt[ C_ENT[ TL_TYPE_ELEMENT ].N_DIM-1 ];

          for( unsigned short l_di = 0; l_di < C_ENT[ TL_TYPE_ELEMENT ].N_DIM-1; l_di++ ) {
            l_faPt[l_di] = l_faPtsChi[l_di][l_qp];
          }

          // set local quad points (volume coordinates)
          linalg::Mappings::faToVolRef( l_fa, TL_TYPE_ELEMENT, l_faPt, o_pts[l_fa][l_qp] );

          // set "neighboring" quad points
          for( unsigned short l_ve = 0; l_ve < CE_N_FACE_VERTEX_OPTS( TL_TYPE_ELEMENT ); l_ve++ ) {
            // local face coords
            TL_T_REAL_MESH l_faPtL[ C_ENT[ TL_TYPE_ELEMENT ].N_DIM-1  ];
            linalg::Mappings::faLocToFaNei( TL_TYPE_ELEMENT, l_faPt, l_faPtL, l_ve );

            // set quad points (volume coordinates)
            unsigned short l_pos =  C_ENT[TL_TYPE_ELEMENT].N_FACES;
                           l_pos += l_fa * CE_N_FACE_VERTEX_OPTS( TL_TYPE_ELEMENT );
                           l_pos += l_ve;
            linalg::Mappings::faToVolRef( l_fa, TL_TYPE_ELEMENT, l_faPtL, o_pts[l_pos][l_qp] );
          }
        }
      }

      // set quadrature weights
      for( unsigned int l_qp = 0; l_qp < CE_N_FACE_QUAD_POINTS( TL_TYPE_ELEMENT, TL_ORDER ); l_qp++ ) {
        o_weights[l_qp] = l_weights[l_qp];
      };

      // evaluate basis at face quad points
      for( unsigned short l_op = 0; l_op < (CE_N_FACE_VERTEX_OPTS(TL_TYPE_ELEMENT)+1) * C_ENT[TL_TYPE_ELEMENT].N_FACES; l_op++ ) {
        for( unsigned short l_qp = 0; l_qp < CE_N_FACE_QUAD_POINTS( TL_TYPE_ELEMENT, TL_ORDER ); l_qp++ ) {
          for( unsigned int   l_md = 0;  l_md < CE_N_ELEMENT_MODES( TL_TYPE_ELEMENT, TL_ORDER ); l_md++ ) {
            dg::Basis::evalBasis( l_md,
                                  TL_TYPE_ELEMENT,
                                  o_basisEval[l_op][l_qp][l_md],
                                  o_pts[l_op][l_qp],
                                  -1, TL_ORDER );
          }
        }
      }
    }

    /**
     * Computes the point-wise solution by using the evaluated basis functions at the point.
     *
     * @param i_basisEval evaluated basis functions at the point.
     * @param i_dofs modal degrees of freedom.
     * @param o_quadEval will be set to evaluated solution at the implicitly given point.
     *
     * @paramt TL_T_REAL precision of the arithmetic operations.
     **/
    template< typename TL_T_REAL >
    static void inline evalBasis( TL_T_REAL      const i_basisEval[ CE_N_ELEMENT_MODES( TL_TYPE_ELEMENT, TL_ORDER ) ],
                                  TL_T_REAL      const i_dofs[      CE_N_ELEMENT_MODES( TL_TYPE_ELEMENT, TL_ORDER ) ][ TL_N_CRUNS ],
                                  TL_T_REAL            o_quadEval[  TL_N_CRUNS                                      ]   ) {
      // set initial values
      for( int_cfr l_ru = 0; l_ru < TL_N_CRUNS; l_ru++ ) {
        o_quadEval[l_ru] = i_dofs[0][l_ru] * i_basisEval[0];
      }

      // iterate over the modes and complete eval
      for( int_md l_md = 1; l_md < CE_N_ELEMENT_MODES( TL_TYPE_ELEMENT, TL_ORDER ); l_md++ ) {
        for( int_cfr l_ru = 0; l_ru < TL_N_CRUNS; l_ru++ ) {
          o_quadEval[l_ru] += i_dofs[l_md][l_ru] * i_basisEval[l_md];
        }
      }
    }
};

template <unsigned short TL_ORDER, unsigned short TL_N_CRUNS>
class edge::dg::QuadratureEval< LINE, TL_ORDER, TL_N_CRUNS >: public QuadratureEvalBase<TL_ORDER> {
  public:
    /**
     * Dummy implementation for one-dimensional face quad points of the line element.
     *
     * @param o_pts will be set to quadrature points. [*][][]: flux options, [][*][]: quad point [][][*]: dim.
     * @param o_weights will be set to weights of the quad points.
     * @param o_basisEval will be set to evaluated basis functiosn (volume/element) at the quad points.
     *
     * @paramt TL_T_REAL_MESH floating point type of the points.
     * @paramt TL_T_REAL_COMP floating point type of the basis.
     **/
    template< typename TL_T_REAL_MESH, typename TL_T_REAL_COMP >
    static void faces( TL_T_REAL_MESH o_pts[4][1][1],
                       TL_T_REAL_COMP o_weights[1],
                       TL_T_REAL_COMP o_basisEval[4][1][TL_ORDER] ) {
      o_pts[0][0][0] = C_REF_ELEMENT.VE.ENT[ LINE ][0];
      o_pts[1][0][0] = C_REF_ELEMENT.VE.ENT[ LINE ][1];
      o_pts[2][0][0] = C_REF_ELEMENT.VE.ENT[ LINE ][1];
      o_pts[3][0][0] = C_REF_ELEMENT.VE.ENT[ LINE ][0];
      o_weights[0] = 1;

      // evaluate basis
      for( unsigned short l_op = 0; l_op < 4; l_op++ ) {
          for( unsigned short l_md = 0; l_md < TL_ORDER; l_md++ ) {
            dg::Basis::evalBasis( l_md,
                                  LINE,
                                  o_basisEval[l_op][0][l_md],
                                  o_pts[l_op][0] );
          }
      }
    }
};
#endif
