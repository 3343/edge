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
 * Initial values used for the Degrees Of Freedom.
 **/
#ifndef EDGE_INITIAL_DOFS_HPP
#define EDGE_INITIAL_DOFS_HPP

#include "constants.hpp"
#include "io/logging.h"
#include "data/Expression.hpp"
#include "mesh/common.hpp"
#include "dg/QuadraturePoints.h"
#include "sc/SubGrid.hpp"

namespace edge {
  namespace setups {
    template< t_entityType   TL_T_EL,
              unsigned short TL_O_SP,
              unsigned short TL_N_QTS,
              unsigned short TL_N_CRS >
    class InitialDofs;
  }
}

/**
 * Initial values used for the Degrees Of Freedom.
 *
 * @paramt TL_T_EL element type.
 * @paramt TL_O_SP order of the used solver.
 * @paramt TL_N_QTS number of quantities.
 * @paramt TL_N_CRS number of fused runs.
 **/
template< t_entityType   TL_T_EL,
          unsigned short TL_O_SP,
          unsigned short TL_N_QTS,
          unsigned short TL_N_CRS >
class edge::setups::InitialDofs {
  private:
    //! number of dimensions
    static unsigned short const TL_N_DIMS = C_ENT[ TL_T_EL ].N_DIM;

    //! number of vertices
    static unsigned short const TL_N_VES = C_ENT[ TL_T_EL ].N_VERTICES;

    //! number of sub-cells
    static unsigned short const TL_N_SCS = CE_N_SUB_CELLS( TL_T_EL, TL_O_SP );

    //! number fo element modes.
    static unsigned short const TL_N_MDS = CE_N_ELEMENT_MODES( TL_T_EL, TL_O_SP );

    //! number of quad points for order+1
    static unsigned int const TL_N_QPS1 = CE_N_QUAD_POINTS( TL_T_EL, TL_O_SP+1 );

    /**
     * Determines the quadrature points and weights for an element.
     *
     * @param i_order order of the quadrature rule.
     * @param i_el element's id.
     * @param i_elVe vertices adjacent to the elements.
     * @param i_veChars vertex characteristics.
     * @param o_pts will be set to quadrature points.
     * @param o_wes will be set to quadrature weights.
     *
     * @paramt TL_T_LID integral type of local ids.
     * @paramt TL_T_REAL floating point type.
     * @paramt TL_T_VE_CHARS vertex characteristics, offering coordinates through .coords.
     **/
    template< typename TL_T_LID,
              typename TL_T_REAL,
              typename TL_T_VE_CHARS >
    static void qps( unsigned short                   i_order,
                     TL_T_LID                         i_el,
                     TL_T_LID                 const (*i_elVe)[TL_N_VES],
                     TL_T_VE_CHARS            const  *i_veChars,
                     std::vector< TL_T_REAL >         o_pts[3],
                     std::vector< TL_T_REAL >        &o_wes ) {
      // get the element's vertices
      TL_T_REAL l_veCrds[3][TL_N_VES];

      // init with zero, TODO: 3 dims are work-around to ensure compability with quad-pt computation
      for( unsigned short l_di = 0; l_di < 3; l_di++ )
        for( unsigned short l_ve = 0; l_ve < TL_N_VES; l_ve++ ) l_veCrds[l_di][l_ve] = 0;

      mesh::common<
        TL_T_EL
      >::getElVeCrds( i_el,
                      i_elVe,
                      i_veChars,
                      l_veCrds );

      // get quad points and weights
      dg::QuadraturePoints::getQpts( TL_T_EL,
                                     i_order,
                                     l_veCrds,
                                     o_pts[0], o_pts[1], o_pts[2], o_wes );
    }

    /**
     * 1) Binds the coordinates as input, and the quantities as output to the expressions.
     * 2) Compiles the expressions.
     *
     * @param i_exprStrs expression strings.
     * @param i_crds memory location of coordinates which is used as input.
     * @param o_qts memory location of quantities which is used as output.
     * @param io_exprs expressions to which the memory locations are bound.
     *
     * @paramt TL_T_REAL floating point type.
     **/
    template< typename TL_T_REAL >
    static void bc( std::string                         const i_exprStrs[TL_N_CRS],
                    TL_T_REAL                                 i_crds[TL_N_DIMS],
                    TL_T_REAL                                 o_qts[TL_N_CRS][TL_N_QTS],
                    edge::data::Expression< TL_T_REAL >       io_exprs[TL_N_CRS] ) {
      // bind variables and compile expressions
      for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {
        io_exprs[l_cr].bindCrds( i_crds, TL_N_DIMS );
        io_exprs[l_cr].bind( "q", o_qts[l_cr], TL_N_QTS );

        io_exprs[l_cr].compile( i_exprStrs[l_cr] );
      }
    }

  public:
    /**
     * Initializes the DG-DOFs of the elements according to the given expressions.
     * We use numerical quadrature is used to query the expressions point-wise before
     * projecting the initial solution to DG-space. 
     *
     * @param i_first id of the first element.
     * @param i_size number of elements.
     * @param i_exprStrs expressiosn strings. symbols x, y, z for the coordinates are provided, symboles q[0], .., q[TL_N_QTS-1] for the quantities read.
     * @param i_basis DG basis.
     * @param i_elVe vertices adjacents to elements.
     * @param i_veChars vertex characteristics.
     * @param o_dofs will be set to DOFs.
     *
     * @paramt TL_T_LID integral type for local ids.
     * @paramt TL_T_REAL floating point type.
     * @paramt TL_T_VE_CHARS type of the vertex characteristics with member .coords giving the vertices' coordinates.
     **/
    template< typename TL_T_LID,
              typename TL_T_REAL,
              typename TL_T_VE_CHARS >
    static void dg( TL_T_LID              i_first,
                    TL_T_LID              i_size,
                    std::string   const   i_exprStrs[TL_N_CRS],
                    dg::Basis     const  &i_basis,
                    TL_T_LID            (*i_elVe)[TL_N_VES],
                    TL_T_VE_CHARS const  *i_veChars,
                    TL_T_REAL           (*o_dofs)[TL_N_QTS][TL_N_MDS][TL_N_CRS] ) {
      // iterate over DG-elements
#ifdef PP_USE_OMP
#pragma omp parallel
#endif
      {
        // coordinates
        TL_T_REAL l_crds[TL_N_DIMS];

        // quantities
        TL_T_REAL l_qts[TL_N_CRS][TL_N_QTS];

        // expressions
        edge::data::Expression< TL_T_REAL > l_exprs[TL_N_CRS];

        // bind and compile expressions
        bc( i_exprStrs, l_crds, l_qts, l_exprs );

#ifdef PP_USE_OMP
#pragma omp for
#endif
        for( TL_T_LID l_el = i_first; l_el < i_first+i_size; l_el++ ) {
          // get quad-points
          std::vector< TL_T_REAL > l_pts[3], l_wes;
          qps( TL_O_SP+1,
              l_el,
              i_elVe,
              i_veChars,
              l_pts,
              l_wes );
          EDGE_CHECK( l_wes.size() == TL_N_QPS1 ); // check compability of work-around

          // solution at quad points
          TL_T_REAL l_q0[TL_N_CRS][TL_N_QTS][TL_N_QPS1];

          // get solution at quad points
          for( unsigned short l_qp = 0; l_qp < TL_N_QPS1; l_qp++ ) {
            // copy coordinates
            for( unsigned short l_di = 0; l_di < TL_N_DIMS; l_di++ )
              l_crds[l_di] = l_pts[l_di][l_qp];

            // eval expressions and store results
            for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {
              // reset quantities
              for( unsigned short l_qt = 0; l_qt < TL_N_QTS; l_qt++ )
                l_qts[l_cr][l_qt]  = 0;

              l_exprs[l_cr].eval();

              for( unsigned short l_qt = 0; l_qt < TL_N_QTS; l_qt++ )
                l_q0[l_cr][l_qt][l_qp] = l_qts[l_cr][l_qt];
            }
          }

          // project to modes
          for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {
            for( unsigned short l_qt = 0; l_qt < TL_N_QTS; l_qt++ ) {
              TL_T_REAL l_mds[TL_N_MDS];

              i_basis.qpts2modal( l_q0[l_cr][l_qt],
                                  TL_O_SP+1,
                                  l_mds );

              // store modes
              for( unsigned short l_md = 0; l_md < TL_N_MDS; l_md++ ) {
                o_dofs[l_el][l_qt][l_md][l_cr] = l_mds[l_md];
              }

            }
          }

        }
      }
    }

    /**
     * Initializes the sub-cell-DOFs of limited limited elements according to the given expressions.
     * The expressions are evaluated at the vertices of the sub-cells and the sub-cell values are
     * set to accordingly to the average of all adjacent sub-vertices.
     *
     * @param i_first first DG-element.
     * @param i_size number of DG elements.
     * @param i_liFirst first limited element.
     * @param i_spType sparse type indicating a limited element.
     * @param i_exprStrs expressions strings which are evaluated at sub-vertices of limited elements.
     * @param i_elVe vertices adjacent to DG-elements.
     * @param i_scSv sub-vertices adjacent to sub-cells.
     * @param i_veChars characteristics of DG-vertices (prodividing their spatial coordinates).
     * @param i_svChars characteristics of sub-vertices (providing their spatial coordinates).
     * @param i_elChars characteristics of the DG-elements (providing their sparse type).
     * @param o_dofsSc sub-cell DOFs, which will be initialized.
     *
     * @paramt TL_T_LID integral type of local ids.
     * @paramt TL_T_REAL floating point type.
     * @paramt TL_T_SP sparse type.
     * @paramt TL_T_VE_CHARS vertex characteristics, offering coordinates through .coords.
     * @paramt TL_T_SV_CHARS sub-vertex characteristics, offering coordinates through .coords.
     * @paramt TL_T_EL_chars element characteristics, offering sparse type through .spType.
     **/
    template< typename TL_T_LID,
              typename TL_T_REAL,
              typename TL_T_SP,
              typename TL_T_VE_CHARS,
              typename TL_T_SV_CHARS,
              typename TL_T_EL_CHARS >
    static void sc( TL_T_LID              i_first,
                    TL_T_LID              i_size,
                    TL_T_LID              i_liFirst,
                    TL_T_SP               i_spType,
                    std::string   const   i_exprStrs[TL_N_CRS],
                    TL_T_LID      const (*i_elVe)[TL_N_VES],
                    TL_T_LID      const (*i_scSv)[TL_N_VES],
                    TL_T_VE_CHARS const  *i_veChars,
                    TL_T_SV_CHARS const  *i_svChars,
                    TL_T_EL_CHARS const  *i_elChars,
                    TL_T_REAL           (*o_dofsSc)[TL_N_QTS][TL_N_SCS][TL_N_CRS] ) {
      // coordinates
      TL_T_REAL l_crds[TL_N_DIMS];

      // quantities
      TL_T_REAL l_qts[TL_N_CRS][TL_N_QTS];

      // expressions
      edge::data::Expression< TL_T_REAL > l_exprs[TL_N_CRS];

      // bind and compile expressions
      bc( i_exprStrs, l_crds, l_qts, l_exprs );

      // init limited elements
      TL_T_LID l_li = i_liFirst;

      // averaging scalar
      TL_T_REAL l_sca = 1; l_sca /= TL_N_VES;

      // iterate over DG-elements
      for( TL_T_LID l_el = i_first; l_el < i_first+i_size; l_el++ ) {
        // check if DG-element is limited
        if( (i_elChars[l_el].spType & i_spType) == i_spType ) {
          // get vertex coordinates of the DG-element
          TL_T_REAL l_veCrds[TL_N_DIMS][TL_N_VES];
          mesh::common<
            TL_T_EL
          >::getElVeCrds( l_el,
                          i_elVe,
                          i_veChars,
                          l_veCrds );

          // iterate over sub-cells
          for( unsigned short l_sc = 0; l_sc < TL_N_SCS; l_sc++ ) {
            // init sc-dofs
            for( unsigned short l_qt = 0; l_qt < TL_N_QTS; l_qt++ ) {
              for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {
                o_dofsSc[l_li][l_qt][l_sc][l_cr] = 0;
              }
            }

            // iterate over sub-vertices of the sub-cell
            for( unsigned short l_ve = 0; l_ve < TL_N_VES; l_ve++ ) {
              // get sub-vertex id
              TL_T_LID l_sv = i_scSv[l_sc][l_ve];

              // get physical coordinates of the vertex
              linalg::Mappings::refToPhy( TL_T_EL,
                                          l_veCrds[0],
                                          i_svChars[l_sv].coords,
                                          l_crds );

              // eval expressions of all forward runs and add results
              for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {
                // reset quantities
                for( unsigned short l_qt = 0; l_qt < TL_N_QTS; l_qt++ )
                  l_qts[l_cr][l_qt] = 0;

                l_exprs[l_cr].eval();

                for( unsigned short l_qt = 0; l_qt < TL_N_QTS; l_qt++ )
                  o_dofsSc[l_li][l_qt][l_sc][l_cr] += l_qts[l_cr][l_qt] * l_sca;
              }

            }
          }

          l_li++;
        }
      }

    }

    /**
     * Computes the L1, L2, and Linf error of the given, possibly limited solution.
     *
     * @param i_first firse DG-element.
     * @param i_size number of DG-elements.
     * @param i_liFirst first limited DG-element.
     * @param i_spType sparse type of limited elements.
     * @param i_exprStrs expressions string, encoding the reference solution.
     * @param i_basis DG basis.
     * @param i_elVe vertices adjacent to DG-elements (no bridge).
     * @param i_scSv sub-vertices adjacent to sub-cells (no bridge).
     * @param i_veChars vertex characteristics.
     * @param i_svChars vertex characteristics of sub-cell-vertices.
     * @param i_elChars element characteristics.
     * @param i_adm admissibility of limited elements.
     * @param i_dofsDg DOFs of the DG-solution.
     * @param i_dofsSc DOFs of the sub-cell solution.
     * @param o_l1 will be set to L1 error.
     * @param o_l2p2 will be set to squared (to the power of 2) L2 error.
     * @param o_lInf will be set to Linf error.
     *
     * @paramt TL_T_LID integral type of local ids.
     * @paramt TL_T_REAL floating point type.
     * @paramt TL_T_SP sparse type.
     * @paramt TL_T_VE_CHARS vertex characteristics, offering coordinates through .coords.
     * @paramt TL_T_SV_CHARS sub-vertex characteristics, offering coordinates through .coords.
     * @paramt TL_T_EL_chars element characteristics, offering sparse type through .spType.
     **/
    template< typename TL_T_LID,
              typename TL_T_REAL,
              typename TL_T_SP,
              typename TL_T_VE_CHARS,
              typename TL_T_SV_CHARS,
              typename TL_T_EL_CHARS >
    static void err( TL_T_LID               i_first,
                     TL_T_LID               i_size,
                     TL_T_LID               i_liFirst,
                     TL_T_SP                i_spType,
                     std::string    const   i_exprStrs[TL_N_CRS],
                     dg::Basis      const  &i_basis,
                     TL_T_LID       const (*i_elVe)[TL_N_VES],
                     unsigned short const (*i_scSv)[TL_N_VES],
                     TL_T_VE_CHARS  const  *i_veChars,
                     TL_T_SV_CHARS  const  *i_svChars,
                     TL_T_EL_CHARS  const  *i_elChars,
                     bool           const (*i_adm)[TL_N_CRS],
                     TL_T_REAL      const (*i_dofsDg)[TL_N_QTS][TL_N_MDS][TL_N_CRS],
                     TL_T_REAL      const (*i_dofsSc)[TL_N_QTS][TL_N_SCS][TL_N_CRS],
                     double                 o_l1[TL_N_QTS][TL_N_CRS],
                     double                 o_l2p2[TL_N_QTS][TL_N_CRS],
                     double                 o_lInf[TL_N_QTS][TL_N_CRS] ) {
      // init error norms
      for( unsigned short l_qt = 0; l_qt < TL_N_QTS; l_qt++ )
        for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ )
          o_l1[l_qt][l_cr] = o_l2p2[l_qt][l_cr] = o_lInf[l_qt][l_cr] = 0;

      // limited DG-element id
      TL_T_LID l_li = i_liFirst;

      // coordinates and quantities of the expression
      double l_cE[TL_N_DIMS];
      double l_qE[TL_N_CRS][TL_N_QTS];

      // expressions
      edge::data::Expression< double > l_exprs[TL_N_CRS];

      // bind and compile expressions
      bc( i_exprStrs, l_cE, l_qE, l_exprs );

      // get quad-points in reference coordinates
      std::vector< double > l_ptsR[3], l_wesR;
      dg::QuadraturePoints::getQpts( TL_T_EL,
                                     TL_O_SP+1,
                                     C_REF_ELEMENT.VE.ENT[TL_T_EL],
                                     l_ptsR[0], l_ptsR[1], l_ptsR[2], l_wesR );
      EDGE_CHECK( l_wesR.size() == TL_N_QPS1 ); // check compability of work-around

      // determine sub-cells closest to the quad points
      TL_T_LID l_qpSc[TL_N_QPS1];
      for( unsigned short l_qp = 0; l_qp < TL_N_QPS1; l_qp++ ) {
        // TODO: work-around for SoA
        double l_pt[TL_N_DIMS];
        for( unsigned short l_di = 0; l_di < TL_N_DIMS; l_di++ ) l_pt[l_di] = l_ptsR[l_di][l_qp];

        // set sub-cell
        l_qpSc[l_qp] = sc::SubGrid< TL_T_EL, TL_O_SP >::ptSc( l_pt, i_scSv, i_svChars );
      }


      // iterate over DG-elements
      for( TL_T_LID l_el = i_first; l_el < i_first+i_size; l_el++ ) {
        // get quad-points in physical coordinates
        std::vector< double > l_ptsP[3], l_wesP;
        qps( TL_O_SP+1,
             l_el,
             i_elVe,
             i_veChars,
             l_ptsP,
             l_wesP );
        // check compability of work-around
        EDGE_CHECK_EQ( l_wesP.size(), TL_N_QPS1 );

        // reference solution
        double l_qR[TL_N_CRS][TL_N_QTS][TL_N_QPS1];
        // compute ref solution at quad points
        for( unsigned short l_qp = 0; l_qp < TL_N_QPS1; l_qp++ ) {
          // copy over coordinates
          for( unsigned short l_di = 0; l_di < TL_N_DIMS; l_di++ )
            l_cE[l_di] = l_ptsP[l_di][l_qp];

          // get ref solution
          for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {
            // eval expression
            l_exprs[l_cr].eval();
            // save results
            for( unsigned short l_qt = 0; l_qt < TL_N_QTS; l_qt++ )
              l_qR[l_cr][l_qt][l_qp] = l_qE[l_cr][l_qt];
          }
        }

        // iterate over fused sims, get numerical solution, and update errors
        for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {
          // numerical solution at quad points
          TL_T_REAL l_qN[TL_N_QTS][TL_N_QPS1];

          // limited element
          if( (i_elChars[l_el].spType & i_spType) == i_spType &&
               i_adm[l_li][l_cr] == false ) {
            // iterate over quad points
            for( unsigned short l_qp = 0; l_qp < TL_N_QPS1; l_qp++ ) {
              // id of corresponding sub-cell
              TL_T_LID l_sc = l_qpSc[l_qp];

              // set numerical solution
              for( unsigned short l_qt = 0; l_qt < TL_N_QTS; l_qt++ )
                l_qN[l_qt][l_qp] = i_dofsSc[l_li][l_qt][l_sc][l_cr];
            }
          }
          else {
            // iterate over quantities
            for( unsigned short l_qt = 0; l_qt < TL_N_QTS; l_qt++ ) {
              // gather modes
              TL_T_REAL l_mds[TL_N_MDS];
              for( unsigned short l_md = 0; l_md < TL_N_MDS; l_md++ )
                l_mds[l_md] = i_dofsDg[l_el][l_qt][l_md][l_cr];

              // set numerical solution
              for( unsigned short l_qp = 0; l_qp < TL_N_QPS1; l_qp++ )
                l_qN[l_qt][l_qp] = i_basis.modal2refPtVal( TL_O_SP+1, l_qp, l_mds );
            }
          }

          // compute errors
          for( unsigned short l_qt = 0; l_qt < TL_N_QTS; l_qt++ ) {
            for( unsigned short l_qp = 0; l_qp < TL_N_QPS1; l_qp++ ) {
              double l_diff = std::abs( l_qN[l_qt][l_qp] - l_qR[l_cr][l_qt][l_qp] );
              o_l1[l_qt][l_cr]   += l_diff *          l_wesP[l_qp];
              o_l2p2[l_qt][l_cr] += l_diff * l_diff * l_wesP[l_qp];
              o_lInf[l_qt][l_cr]  = std::max( o_lInf[l_qt][l_cr], l_diff );
            }
          }
        }

        // increase limited DG-element id if required
        if( (i_elChars[l_el].spType & i_spType) == i_spType ) l_li++;
      }
    }
};

#endif
