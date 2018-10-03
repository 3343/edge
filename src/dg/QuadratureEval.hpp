/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2017-2018, Regents of the University of California
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

#ifndef EDGE_DG_QUADRATURE_EVAL_HPP
#define EDGE_DG_QUADRATURE_EVAL_HPP

#include "constants.hpp"
#include "io/logging.h"
#include "QuadraturePoints.h"
#include "linalg/Mappings.hpp"
#include "dg/Basis.h"

namespace edge {
  namespace dg {
    template< unsigned short TL_O_SP >
    class QuadratureEvalBase;

    template< t_entityType   TL_T_EL,
              unsigned short TL_O_SP,
              unsigned short TL_N_CRS=1 >
    class QuadratureEval;

    template< unsigned short TL_O_SP,
              unsigned short TL_N_CRS >
    class QuadratureEval< LINE, TL_O_SP, TL_N_CRS >;
  }
}

template< unsigned short TL_O_SP >
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
    static void line(  TL_T_REAL_MESH o_pts[ CE_N_ELEMENT_MODES( LINE, TL_O_SP ) ],
                       TL_T_REAL_COMP o_weights[ CE_N_ELEMENT_MODES( LINE, TL_O_SP ) ] ) {
      // get 1D quadrature points and weights
      std::vector< TL_T_REAL_MESH > l_qps[3];
      std::vector< TL_T_REAL_COMP > l_weights;

      dg::QuadraturePoints::getQpts( LINE,
                                     TL_O_SP,
                                     C_REF_ELEMENT.VE.ENT[LINE],
                                     l_qps[0], l_qps[1], l_qps[2], l_weights );

      EDGE_CHECK( l_weights.size() == CE_N_ELEMENT_MODES( LINE, TL_O_SP ) );

      for( unsigned short l_qp = 0; l_qp < CE_N_ELEMENT_MODES( LINE, TL_O_SP ); l_qp++ ) {
         o_pts[l_qp] = l_qps[0][l_qp];
         o_weights[l_qp] = l_weights[l_qp];
      }
    }
};

template< t_entityType   TL_T_EL,
          unsigned short TL_O_SP,
          unsigned short TL_N_CRS >
class edge::dg::QuadratureEval: public QuadratureEvalBase< TL_O_SP > {
  private:
    //! number of dimensions
    static const unsigned short TL_N_DIS = C_ENT[TL_T_EL].N_DIM;
    //! number of vertex options
    static const unsigned short TL_N_VE_OPTS = CE_N_FACE_VERTEX_OPTS(TL_T_EL);
    //! number of quadrature points per face
    static const unsigned short TL_N_QPS_FA = CE_N_FACE_QUAD_POINTS( TL_T_EL, TL_O_SP );
    //! number of quadrature points per element
    static const unsigned short TL_N_QPS_EL = CE_N_QUAD_POINTS( TL_T_EL, TL_O_SP );
    //! number of modes per element
    static const unsigned short TL_N_MDS_EL = CE_N_ELEMENT_MODES( TL_T_EL, TL_O_SP );

  public:
    /**
     * Gets the faces' quadrature points.
     *
     * Remark 1: Collapsed coordinates for tetrahedral elements are used,
     *           resulting in more quadrature points than necessary for the given order of convergence.
     *
     * Remark 2: The weights are based on the face type's reference element.
     *           The sum gives the volume of this element.
     *           The face type is given implicitly by the template parameter TL_T_EL.
     *
     * @param o_pts will be set to quadrature points. [*][][]: flux options, [][*][]: quad point [][][*]: dim.
     * @param o_weights will be set to weights of the quad points.
     * @param o_basisEval will be set to evaluated basis functiosn (volume/element) at the quad points.
     *
     * @paramt TL_T_REAL_MESH floating point type of the points.
     * @paramt TL_T_REAL_COMP floating point type of the basis.
     **/
#ifndef __INTEL_COMPILER
    template< typename TL_T_REAL_MESH,
              typename TL_T_REAL_COMP >
    static void faces( TL_T_REAL_MESH o_pts[ (TL_N_VE_OPTS+1) * C_ENT[TL_T_EL].N_FACES ]
                                           [  TL_N_QPS_FA                              ]
                                           [  TL_N_DIS                                 ],
                       TL_T_REAL_COMP o_weights[ TL_N_QPS_FA ],
                       TL_T_REAL_COMP o_basisEval[ (TL_N_VE_OPTS+1) * C_ENT[TL_T_EL].N_FACES ]
                                                 [  TL_N_QPS_FA                        ]
                                                 [  TL_N_MDS_EL                        ] ) {
#else
    // TODO: Hard-coded work around for Intel compiler.
    //       Fails to compile the template above; icpc version 17.0.3 (gcc version 4.8.5 compatibility)
    static void faces( real_mesh o_pts[ (TL_N_VE_OPTS+1) * C_ENT[TL_T_EL].N_FACES ]
                                      [  TL_N_QPS_FA                              ]
                                      [  TL_N_DIS                     ],
                       real_base o_weights[ TL_N_QPS_FA ],
                       real_base o_basisEval[ (TL_N_VE_OPTS+1) * C_ENT[TL_T_EL].N_FACES ]
                                            [  TL_N_QPS_FA                              ]
                                            [  TL_N_MDS_EL                              ] ) {
      typedef real_mesh TL_T_REAL_MESH;
      typedef real_base TL_T_REAL_COMP;
#endif
      // get quad points in face local coords
      std::vector< TL_T_REAL_MESH > l_faPtsChi[3];
      std::vector< TL_T_REAL_COMP > l_weights;

      dg::QuadraturePoints::getQpts( C_ENT[TL_T_EL].TYPE_FACES,
                                     TL_O_SP,
                                     C_REF_ELEMENT.VE.ENT[C_ENT[TL_T_EL].TYPE_FACES],
                                     l_faPtsChi[0], l_faPtsChi[1], l_faPtsChi[2], l_weights );

      // double-check number of faces' quad points
      EDGE_CHECK( l_weights.size() == TL_N_QPS_FA );

      // iterate over faces and set volume cords of quad points
      for( unsigned short l_fa = 0; l_fa < C_ENT[TL_T_EL].N_FACES; l_fa++ ) {
        for( unsigned int l_qp = 0; l_qp < TL_N_QPS_FA; l_qp++ ) {
          TL_T_REAL_MESH l_faPt[ C_ENT[ TL_T_EL ].N_DIM-1 ];

          for( unsigned short l_di = 0; l_di < C_ENT[ TL_T_EL ].N_DIM-1; l_di++ ) {
            l_faPt[l_di] = l_faPtsChi[l_di][l_qp];
          }

          // set local quad points (volume coordinates)
          linalg::Mappings::faToVolRef( l_fa, TL_T_EL, l_faPt, o_pts[l_fa][l_qp] );

          // set "neighboring" quad points
          for( unsigned short l_ve = 0; l_ve < CE_N_FACE_VERTEX_OPTS( TL_T_EL ); l_ve++ ) {
            // local face coords
            TL_T_REAL_MESH l_faPtL[ C_ENT[ TL_T_EL ].N_DIM-1  ];
            linalg::Mappings::faLocToFaNei( TL_T_EL, l_faPt, l_faPtL, l_ve );

            // set quad points (volume coordinates)
            unsigned short l_pos =  C_ENT[TL_T_EL].N_FACES;
                           l_pos += l_fa * CE_N_FACE_VERTEX_OPTS( TL_T_EL );
                           l_pos += l_ve;
            linalg::Mappings::faToVolRef( l_fa, TL_T_EL, l_faPtL, o_pts[l_pos][l_qp] );
          }
        }
      }

      // set quadrature weights
      for( unsigned int l_qp = 0; l_qp < TL_N_QPS_FA; l_qp++ ) {
        o_weights[l_qp] = l_weights[l_qp];
      };

      // evaluate basis at face quad points
      for( unsigned short l_op = 0; l_op < (TL_N_VE_OPTS+1) * C_ENT[TL_T_EL].N_FACES; l_op++ ) {
        for( unsigned short l_qp = 0; l_qp < TL_N_QPS_FA; l_qp++ ) {
          for( unsigned int   l_md = 0;  l_md < TL_N_MDS_EL; l_md++ ) {
            dg::Basis::evalBasis( l_md,
                                  TL_T_EL,
                                  o_basisEval[l_op][l_qp][l_md],
                                  o_pts[l_op][l_qp],
                                  -1, TL_O_SP );
          }
        }
      }
    }

    /**
     * @brief Derives quadrature points, weights and evaluated basis functions w.r.t. the reference element.
     * 
     * @param o_pts will be set to the coordinates of the quadrature points.
     * @param o_wes will be set to the weights of the quadrature points.
     * @param o_eva will be set to the evaluated basis functions at the quadrature points.
     *
     * @paramt TL_T_REAL floating precision.
     */
    template< typename TL_T_REAL >
    static void element( TL_T_REAL o_pts[TL_N_QPS_EL][TL_N_DIS],
                         TL_T_REAL o_wes[TL_N_QPS_EL],
                         TL_T_REAL o_eva[TL_N_QPS_EL][TL_N_MDS_EL] ) {
      // get the quad points and weights
      std::vector< double > l_pts[3];
      std::vector< double > l_wes;

      dg::QuadraturePoints::getQpts( TL_T_EL,
                                     TL_O_SP,
                                     C_REF_ELEMENT.VE.ENT[TL_T_EL],
                                     l_pts[0], l_pts[1], l_pts[2],
                                     l_wes );

      EDGE_CHECK_EQ( l_pts[0].size(), TL_N_QPS_EL );
      EDGE_CHECK_EQ( l_wes.size(), TL_N_QPS_EL );

      // store data
      for( unsigned short l_qp = 0; l_qp < TL_N_QPS_EL; l_qp++ ) {
        for( unsigned short l_di = 0; l_di < TL_N_DIS; l_di++ ) {
          o_pts[l_qp][l_di] = l_pts[l_di][l_qp];
        }
        o_wes[l_qp] = l_wes[l_qp];
      }

      // evaluate the faces at the quad points
      for( unsigned short l_qp = 0; l_qp < TL_N_QPS_EL; l_qp++ ) {
        for( unsigned short l_md = 0;  l_md < TL_N_MDS_EL; l_md++ ) {
          dg::Basis::evalBasis( l_md,
                                TL_T_EL,
                                o_eva[l_qp][l_md],
                                o_pts[l_qp],
                                -1,
                                TL_O_SP );
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
    static void inline evalBasis( TL_T_REAL const i_basisEval[TL_N_MDS_EL],
                                  TL_T_REAL const i_dofs[TL_N_MDS_EL][TL_N_CRS],
                                  TL_T_REAL       o_quadEval[TL_N_CRS] ) {
      // set initial values
      for( int_cfr l_ru = 0; l_ru < TL_N_CRS; l_ru++ ) {
        o_quadEval[l_ru] = i_dofs[0][l_ru] * i_basisEval[0];
      }

      // iterate over the modes and complete eval
      for( int_md l_md = 1; l_md < TL_N_MDS_EL; l_md++ ) {
        for( int_cfr l_ru = 0; l_ru < TL_N_CRS; l_ru++ ) {
          o_quadEval[l_ru] += i_dofs[l_md][l_ru] * i_basisEval[l_md];
        }
      }
    }
};

template< unsigned short TL_O_SP,
          unsigned short TL_N_CRS >
class edge::dg::QuadratureEval< LINE, TL_O_SP, TL_N_CRS >: public QuadratureEvalBase<TL_O_SP> {
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
                       TL_T_REAL_COMP o_basisEval[4][1][TL_O_SP] ) {
      o_pts[0][0][0] = C_REF_ELEMENT.VE.ENT[ LINE ][0];
      o_pts[1][0][0] = C_REF_ELEMENT.VE.ENT[ LINE ][1];
      o_pts[2][0][0] = C_REF_ELEMENT.VE.ENT[ LINE ][1];
      o_pts[3][0][0] = C_REF_ELEMENT.VE.ENT[ LINE ][0];
      o_weights[0] = 1;

      // evaluate basis
      for( unsigned short l_op = 0; l_op < 4; l_op++ ) {
          for( unsigned short l_md = 0; l_md < TL_O_SP; l_md++ ) {
            dg::Basis::evalBasis( l_md,
                                  LINE,
                                  o_basisEval[l_op][0][l_md],
                                  o_pts[l_op][0] );
          }
      }
    }
};
#endif
