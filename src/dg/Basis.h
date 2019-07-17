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
 * Basis of an element.
 **/

#ifndef BASIS_HPP
#define BASIS_HPP

#include <cassert>
#include <vector>
#include <array>
#include "constants.hpp"
#include "linalg/Matrix.h"

namespace edge {
  namespace pre {
    namespace dg {
      // mass matrix
      extern double const * g_massRaw;
      extern std::size_t const g_massSize;

      // stiffness matrices, premultiplied by the inverse mass matrix
      extern double const * g_stiffVRaw;
      extern std::size_t const g_stiffVSize;

      // transposed stiffness matrices, premultiplied by the inverse mass matrix (after transpose)
      extern double const * g_stiffTRaw;
      extern std::size_t const g_stiffTSize;

      // local contribution flux matrices
      extern double const * g_fluxLRaw;
      extern std::size_t const g_fluxLSize;

      // neighboring contribution flux matrices
      extern double const * g_fluxNRaw;
      extern std::size_t const g_fluxNSize;

      // "transposed" flux matrices (fa -> el basis + inverse mass)
      extern double const * g_fluxTRaw;
      extern std::size_t const g_fluxTSize;
    }
  }
}

namespace edge {
  namespace dg {
    class Basis;
  }
}

class edge::dg::Basis {
  private:
    // Evaluated basis for different orders at a specific quad point.
    struct EvalBasis {
      //! values of the basis; [*]: mode
      std::vector< real_base > val;
      //! values of the partial derivatives; [*][]: dim, [][*]: mode
      std::array< std::vector< real_base >,  3 > valD;
    };

    // Quad points and evaluated basis functions at the quad points of the reference elements for different order.
    struct EvalQpRefElOrder {
      //! reference dimension xi1
      std::vector< real_mesh > xi1;
      //! reference dimension xi2
      std::vector< real_mesh > xi2;
      //! reference dimensions xi2
      std::vector< real_mesh > xi3;
      //! quadrature weigts
      std::vector< real_mesh > weights;
      //! evaluated basis at quad points dependent on the order; [*][]: quad poiints, [][*]: basis of given order
      std::vector< std::vector< EvalBasis > > basis;
    };

    // different evaluated orders of quad rules
    std::vector< EvalQpRefElOrder > m_qpEval;

    //! mass matrix
    t_matCrd m_mass;

    //! flux matrices
    std::vector< t_matCrd > m_flux;

    //! entity for which the basis is defined
    const t_entityType m_entType;

    //! order of convergence
    const unsigned short m_order;

    //! number of basis functions
    unsigned short m_nBaseFuncs;

    /**
     * Initializes the precomputed quad point positions and their basis values.
     *
     * @param i_order maximum order (quadrature & basis) for all ref-elements.
     **/
    void initEvalQpRefEl( unsigned short i_order );

    /**
     * Inits the mass matrix using numerical integration.
     **/
    void initMassMatrix();

    /**
     * Evaluates a basis functions for the line reference element.
     *
     * @param i_b id of the basis function which get evaluated. First basis is 0.
     * @param i_xi location in the reference element, where the basis is evaluated.
     * @param o_val will be set to the value of the basis function at the given point.
     **/
    static void evalBasisLine( unsigned int  i_b,
                               real_mesh     i_xi,
                               real_base    &o_val );

    /**
     * Evaluates the derivative of a basis function for the line reference element.
     *
     * @param i_b id of the basis function which derivative gets evaluated. First id is 0.
     * @param i_xi location in the reference element, where the derivative of the basis function is evaluated.
     * @param o_valDxi will bet set the value of the derivative of the basis function at the given point.
     **/
    static void evalBasisDerLine( unsigned int  i_b,
                                  real_mesh     i_xi,
                                  real_base    &o_valDxi );

    /**
     * Evaluates a basis function for the quadrilateral reference element by going through the tensor structure.
     *
     * @param i_b id of the basis function which gets evaluated. First basis is 0.
     * @param i_xi xi-coord in the reference element.
     * @param i_eta eta-coord in the reference element.
     * @param o_val will be set to the value of the basis function at the given point.
     * @param i_der index for derivative evaluation:
     *              0               : der in first (xi),
     *              1               : der in second (eta),
     *              everything else : no derivatives in both 1D functions
     * @param i_order order of the basis (to be removed when hierachical).
     **/
    static void evalBasisQuad( unsigned int  i_b,
                               real_mesh     i_xi,
                               real_mesh     i_eta,
                               real_base    &o_val,
                               int           i_der = -1,
                               unsigned int  i_order = PP_ORDER );

    /**
     * Evaluates a basis function for the triangular reference element.
     *
     * @param i_b id of the basis function which gets evaluated. First basis is 0.
     * @param i_xi xi-coord in the reference element.
     * @param i_eta eta-coord in the reference element.
     * @param o_val will be set to the value of the basis function at the given point.
     * @param i_der index for derivative evaluation:
     *              0               : evaluate der of basis w.r.t. to xi,
     *              1               : evaluate der of basis w.r.t. to eta,
     *              everything else : no derivative at all.
     **/
    static void evalBasisTria( unsigned int  i_b,
                               real_mesh     i_xi,
                               real_mesh     i_eta,
                               real_base    &o_val,
                               int           i_der = -1 );

    /**
     * Evaluates a basis function for the hexahedral reference element by going through the tensor structure.
     *
     * @param i_b id of the basis function which gets evaluated. First basis is 0.
     * @param i_xi xi-coord in the reference element.
     * @param i_eta eta-coord in the reference element.
     * @param i_zeta zeta-coord in the reference element.
     * @param o_val will be set to the value of the basis function at the given point.
     * @param i_der index for derivative evaluation:
     *              0               : der in first (xi),
     *              1               : der in second (eta),
     *              2               : der in third (zeta),
     *              everything else : no derivatives in both 1D functions
     * @param i_order order of the basis (to be removed when hierachical).
     **/
    static void evalBasisHex( unsigned int i_b,
                              real_mesh    i_xi,
                              real_mesh    i_eta,
                              real_mesh    i_zeta,
                              real_base   &o_val,
                              int          i_der = -1,
                              unsigned int i_order = PP_ORDER );

    /**
     * Evaluates a basis function for the tetrahedral reference element.
     *
     * @param i_b id of the basis function which gets evaluated. First basis is 0.
     * @param i_xi xi-coord in the reference element.
     * @param i_eta eta-coord in the reference element.
     * @param i_zeta zeta-coord in the reference element.
     * @param o_val will be set to the value of the basis function at the given point.
     * @param i_der index for derivative evaluation:
     *              0               : der in first (xi),
     *              1               : der in second (eta),
     *              2               : der in third (zeta)
     *              everything else : no derivatives in both 1D functions,
     **/
    static void evalBasisTet( unsigned int  i_b,
                              real_mesh     i_xi_1,
                              real_mesh     i_xi_2,
                              real_mesh     i_xi_3,
                              real_base    &o_val,
                              int           i_der = -1 );
  public:
    /**
     * Constructur.
     *
     * @param i_entityType entity type of the basis.
     * @param i_order maximum oder of the basis (such that sub-matrices of a hierarchical basis can be queried).
     **/
    Basis( t_entityType   i_entityType,
           unsigned short i_order );

    /**
     * Evaluates a basis function for different reference elements.
     *
     * @param i_b id of the basis function which gets evaluated. First basis is 0.
     * @param i_entType type of the reference element which basis is evaluated.
     * @param o_val will be set to the value of the basis function at the given point.
     * @param i_xi xi-coord in the reference element.
     * @param i_eta eta-coord in the reference element.
     * @param i_zeta zeta-coord in the reference element.
     * @param i_der index for derivative evaluation:
     *              0               : evaluate der of basis w.r.t. to xi,
     *              1               : evaluate der of basis w.r.t. to eta,
     *              2               : evaluate der of basis w.r.t. to zeta,
     *              everything else : no derivative at all.
     * @param i_order order of the basis (to be removed when hierachical).
     **/
    static void evalBasis( unsigned int  i_b,
                           t_entityType  i_entType,
                           real_base    &o_val,
                           real_mesh     i_xi,
                           real_mesh     i_eta  = 0,
                           real_mesh     i_zeta = 0,
                           int           i_der  = -1,
                           unsigned int i_order = PP_ORDER );

    /**
     * Evaluates a basis function for different reference elements.
     *
     * @param i_b id of the basis function which gets evaluated. First basis is 0.
     * @param i_entType type of the reference element for which basis is evaluated.
     * @param o_val will be set to the value of the basis function at the given point.
     * @param i_pt pt where the basis is evaluated.
     * @param i_der index for derivative evaluation:
     *              0               : evaluate der of basis w.r.t. to xi,
     *              1               : evaluate der of basis w.r.t. to eta,
     *              2               : evaluate der of basis w.r.t. to zeta,
     *              everything else : no derivative at all.
     * @param i_order order of the basis (to be removed when hierachical).
     **/
    static void evalBasis( unsigned int        i_b,
                           t_entityType        i_enType,
                           real_base          &o_val,
                           real_mesh*   const  i_pt,
                           int                 i_der  = -1,
                           unsigned int        i_order = PP_ORDER );

    /**
     * Prints the mass matrix, stiffness matrix/matrices and flux matrices.
     **/
    void print() const;

    /**
     * Maps evaluated quadrature points to modal representation.
     *
     * @param i_evalF evaluated function at quadrature points.
     * @param i_orderQr order of the quadrature rule.
     * @param o_modes will bet set to modes.
     **/
    void qpts2modal( const real_base    *i_evalF,
                           unsigned int  i_orderQr,
                           real_base    *o_modes ) const;

    /**
     * Evaluates the modal representation at a given point of the reference element (everything else gets extrapolated).
     *
     * @param i_pt point of the ref. element.
     * @param i_modes modes (coefficients) of the basis.
     **/
    real_base modal2ptval( const real_mesh  i_pt[3],
                           const real_base *i_modes ) const;

    /**
     * Evaluates the modal representation at a given point of the reference element.
     * In contrast to modal2ptval the point matches a (precomputed) quadrature point.
     *
     * @param i_orderQr order of the quadrature rule (implicitly changes quad points).
     * @param i_quadPoi id of the requested quad point evaluation.
     * @param i_modes modes.
     **/
    real_base modal2refPtVal(       unsigned short  i_orderQr,
                                    unsigned int    i_quadPoi,
                              const real_base      *i_modes ) const;

    /**
     * Gets the inverse mass matrix.
     * Format is dense (including zero entries).
     * Storage is row-major or column-major w.r.t. to the i-th basis function phi[i] in:
     * Int[Ref. element] ( phi[row] * phi[col] )
     *
     * @param i_nModes number of modes.
     * @param o_matrix will be set to inverse mass matrix in dense storage: [nmodes][nmodes].
     * @param i_rowMajor true if row major storage is requested.
     **/
    void getMassInvDense( int_md     i_nModes,
                          real_base *o_matrix,
                          bool       i_rowMajor=true ) const;

    /**
     * Gets the stiffness matrices multiplied with the inverse mass matrix.
     * Format is dense (including zero entries).
     * Storage is row-major w.r.t. to the i-th basis function phi[i] in:
     * Int[Ref. element] ( phi[row] * phi[col]_x )
     *
     * @param i_nModes number of modes.
     * @param o_matrices will be set to stiffness matrices in dense storage: [ndim][nmodes][nmodes].
     * @param i_tStiff true if the stiffness matrix should be transposed before multiplication with the inverse mass matrix.
     *
     * @paramt TL_T_REAL floating point precision.
     **/
    template< typename TL_T_REAL >
    void getStiffMm1Dense( int_md     i_nModes,
                           TL_T_REAL *o_matrices,
                           bool       i_tStiff ) const {
      if( i_tStiff == false ) {
        // check that the size matches
        EDGE_CHECK_EQ( pre::dg::g_stiffVSize, i_nModes*i_nModes*C_ENT[m_entType].N_DIM );

        // set values
        for( std::size_t l_va = 0; l_va < pre::dg::g_stiffVSize; l_va++ )
          o_matrices[l_va] = pre::dg::g_stiffVRaw[l_va];
      }
      else {
        // check that size matches
        EDGE_CHECK_EQ( pre::dg::g_stiffTSize, i_nModes*i_nModes*C_ENT[m_entType].N_DIM );

        // set values
        for( std::size_t l_va = 0; l_va < pre::dg::g_stiffTSize; l_va++ )
          o_matrices[l_va] = pre::dg::g_stiffTRaw[l_va];
      }
    }

    /**
     * Gets the flux matrices multiplied with the inverse mass matrix.
     * Format is dense (including zero entries).
     * Order of storage for neighboring contrib. flux matrices is 1) vertices, 2) adjacent face.
     *
     * @param o_fluxL will be set to local contribution flux matrices.
     * @param o_fluxN will be set to neighboring contribution flux matrices.
     * @param i_fluxT will be set to 2D->3D basis change, premultiplied by the inverse mass matrix.
     **/
    template< typename TL_T_REAL >
    void getFluxDense( TL_T_REAL *o_fluxL,
                       TL_T_REAL *o_fluxN,
                       TL_T_REAL *o_fluxT ) const {
      // check for reasonable sizes
      EDGE_CHECK_EQ( pre::dg::g_fluxLSize,
                    C_ENT[T_SDISC.ELEMENT].N_FACES * N_ELEMENT_MODES * N_FACE_MODES  );


      EDGE_CHECK_EQ( pre::dg::g_fluxNSize,
                    N_FLUXN_MATRICES * N_ELEMENT_MODES * N_FACE_MODES  );

      EDGE_CHECK_EQ( pre::dg::g_fluxTSize,
                    C_ENT[T_SDISC.ELEMENT].N_FACES * N_ELEMENT_MODES * N_FACE_MODES  );

      // assign local
      for( std::size_t l_va = 0; l_va < pre::dg::g_fluxLSize; l_va++ )
        o_fluxL[l_va] = pre::dg::g_fluxLRaw[l_va];

      // assign neighboring
      for( std::size_t l_va = 0; l_va < pre::dg::g_fluxNSize; l_va++ )
        o_fluxN[l_va] = pre::dg::g_fluxNRaw[l_va];

      // assign transposed
      for( std::size_t l_va = 0; l_va < pre::dg::g_fluxTSize; l_va++ )
        o_fluxT[l_va] = pre::dg::g_fluxTRaw[l_va];
    }
};

#endif
