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
 * Kernels for sub-cell resolution.
 **/
#ifndef SC_KERNELS_HPP
#define SC_KERNELS_HPP

#include "constants.hpp"
#include "linalg/Matrix.h"

namespace edge {
  namespace sc {
    template< t_entityType   TL_T_EL,
              unsigned short TL_O_SP,
              unsigned short TL_N_QTS,
              unsigned short TL_N_CRS >
    class Kernels;
  }
}

/**
 * Kernels for sub-cell resolution.
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
class edge::sc::Kernels {
  private:
    //! number of sub-faces per element face
    static unsigned short const TL_N_SFS = CE_N_SUB_FACES( TL_T_EL, TL_O_SP );

    //! number of sub-cells
    static unsigned short const TL_N_SCS = CE_N_SUB_CELLS( TL_T_EL, TL_O_SP );

    //! number of modes
    static unsigned short const TL_N_MDS = CE_N_ELEMENT_MODES( TL_T_EL, TL_O_SP );

  public:
    /**
     * Applies the scatter matrix using vanilla kernels.
     * Scatter: Projects DG solution to sub-cell solution.
     *
     * @param i_dofsDg DOFs of DG solution.
     * @param i_mat scatter matrix.
     * @param o_scDofs will be set to degrees of freedom of the sub-cell solution.
     *
     * @paramt TL_T_REAL floating point type.
     **/
    template < typename TL_T_REAL >
    static void scatterVanilla( TL_T_REAL const i_dofsDg[TL_N_QTS][TL_N_MDS][TL_N_CRS],
                                TL_T_REAL const i_mat[TL_N_MDS][TL_N_SCS],
                                TL_T_REAL       o_scDofs[TL_N_QTS][TL_N_SCS][TL_N_CRS] ) {
      // project element's DG solution to sub-cells
      linalg::Matrix::matMulB0FusedAC( TL_N_CRS,
                                       TL_N_QTS, TL_N_SCS, TL_N_MDS, // m, n, k
                                       TL_N_MDS, TL_N_SCS, TL_N_SCS, // ldA, ldB, ldC
                                       i_dofsDg[0][0],               // A
                                       i_mat[0],                     // B
                                       o_scDofs[0][0] );             // C
    }

    /**
     * Applies the scatter matrix for sub-cells at the element's face using vanilla kernels.
     * Scatter: Projects DG solution to sub-cell solution.
     *
     * @param i_dofsDg DOFs of DG solution.
     * @param i_mat scatter matrix.
     * @param o_scDofs will be set to degrees of freedom of the sub-cell solution.
     *
     * @paramt TL_T_REAL floating point type.
     **/
    template < typename TL_T_REAL >
    static void scatterFaVanilla( TL_T_REAL const i_dofsDg[TL_N_QTS][TL_N_MDS][TL_N_CRS],
                                  TL_T_REAL const i_mat[TL_N_MDS][TL_N_SFS],
                                  TL_T_REAL       o_scDofs[TL_N_QTS][TL_N_SFS][TL_N_CRS] ) {
      // project element's DG solution to sub-cells
      linalg::Matrix::matMulB0FusedAC( TL_N_CRS,
                                       TL_N_QTS, TL_N_SFS, TL_N_MDS, // m, n, k
                                       TL_N_MDS, TL_N_SFS, TL_N_SFS, // ldA, ldB, ldC
                                       i_dofsDg[0][0],               // A
                                       i_mat[0],                     // B
                                       o_scDofs[0][0] );             // C
    }

    /**
     * 1) Projects DG solution to sub-cell solution.
     * 2) Replaces projected sub-cell solution with stored sub-cell solution
     *    if admissibility flag is false.
     *
     * @param i_dofsDg DOFs of DG solution.
     * @param i_mat scatter matrix.
     * @param i_dofsSc stored sub-cell solution.
     * @param i_admP admissibility of the DG-solution w.r.t. the previous step.
     * @param o_dofsSc will be set to degrees of freedom of the sub-cell solution.
     *
     * @paramt TL_T_REAL floating point type.
     **/
    template < typename TL_T_REAL >
    static void scatterReplaceVanilla( TL_T_REAL const i_dofsDg[TL_N_QTS][TL_N_MDS][TL_N_CRS],
                                       TL_T_REAL const i_mat[TL_N_MDS][TL_N_SCS],
                                       TL_T_REAL const i_dofsSc[TL_N_QTS][TL_N_SCS][TL_N_CRS],
                                       bool      const i_admP[TL_N_CRS],
                                       TL_T_REAL       o_dofsSc[TL_N_QTS][TL_N_SCS][TL_N_CRS] ) {
      // project DG to sub-cell solution
      scatterVanilla( i_dofsDg,
                      i_mat,
                      o_dofsSc );

      // overwrite project sub-cell solution with stored sub-cell solution
      for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {
        if( i_admP[l_cr] == false ) {
          for( unsigned short l_qt = 0; l_qt < TL_N_QTS; l_qt++ ) {
            for( unsigned short l_sc = 0; l_sc < TL_N_SCS; l_sc++ ) {
              o_dofsSc[l_qt][l_sc][l_cr] = i_dofsSc[l_qt][l_sc][l_cr];
            }
          }
        }
      }
    }

    /**
     * 1) Projects DG solution to sub-cell solution at faces of the DG element.
     * 2) Replaces project sub-cell solution with stored sub-cell solution
     *    if admissibility flag is given and false.
     *
     * @param i_dofsDg DOFs of DG solution.
     * @param i_mat scatter matrix.
     * @param i_dofsSc stored sub-cell solution (if avaliable).
     * @param i_admP admissibility of the DG-solution w.r.t. the previous step. If a nullptr is given, admissibility is assumed for all fused simulations.
     * @param o_dofsSc will be set to degrees of freedom of the sub-cell solution.
     *
     * @paramt TL_T_REAL floating point type.
     **/
    template < typename TL_T_REAL >
    static void scatterReplaceFaVanilla( TL_T_REAL const i_dofsDg[TL_N_QTS][TL_N_MDS][TL_N_CRS],
                                         TL_T_REAL const i_mat[TL_N_MDS][TL_N_SFS],
                                         TL_T_REAL const i_dofsSc[TL_N_QTS][TL_N_SFS][TL_N_CRS],
                                         bool      const i_admP[TL_N_CRS],
                                         TL_T_REAL       o_dofsSc[TL_N_QTS][TL_N_SFS][TL_N_CRS] ) {
      // project DG to sub-cell solution
      scatterFaVanilla( i_dofsDg,
                        i_mat,
                        o_dofsSc );

      // overwrite project sub-cell solution with stored sub-cell solution
      if( i_admP != nullptr ) {
        for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {
          if( i_admP[l_cr] == false ) {
            for( unsigned short l_qt = 0; l_qt < TL_N_QTS; l_qt++ ) {
              for( unsigned short l_sf = 0; l_sf < TL_N_SFS; l_sf++ ) {
                o_dofsSc[l_qt][l_sf][l_cr] = i_dofsSc[l_qt][l_sf][l_cr];
              }
            }
          }
        }
      }
    }

    /**
     * Applies the gather matrix using vanilla kernels.
     * Gather: Projects sub-cell solution to DG-solution.
     *
     * @param i_scDofs DOFs of sub-cell solution.
     * @param i_mat gather matrix.
     * @param o_dgDofs will be set to degrees of freedom of DG solution.
     *
     * @paramt TL_T_REAL floating point type.
     **/
    template < typename TL_T_REAL >
    static void gatherVanilla( TL_T_REAL const i_scDofs[TL_N_QTS][TL_N_SCS][TL_N_CRS],
                               TL_T_REAL const i_mat[TL_N_SCS][TL_N_MDS],
                               TL_T_REAL       o_dgDofs[TL_N_QTS][TL_N_MDS][TL_N_CRS] ) {
      // project element's sub-cell solution to DG solution
      linalg::Matrix::matMulB0FusedAC( TL_N_CRS,
                                       TL_N_QTS, TL_N_MDS, TL_N_SCS, // m, n, k
                                       TL_N_SCS, TL_N_MDS, TL_N_MDS, // ldA, ldB, ldC
                                       i_scDofs[0][0],               // A
                                       i_mat[0],                     // B
                                       o_dgDofs[0][0] );             // C
    }

    /**
     * Applies the sfInt operator using vanilla kernels.
     * sfInt: Integration of piecewise-constant fluxes of a DG-face.
     *
     * @param i_scFluxes fluxes of the sub-cell solution at the face.
     * @param i_mat sfInt matrix.
     * @param o_int will be set to surface integral.
     *
     * @paramt TL_T_REAL floating point type.
     **/
    template< typename TL_T_REAL >
    static void sfIntVanilla( TL_T_REAL const i_scFluxes[TL_N_QTS][TL_N_SFS][TL_N_CRS],
                              TL_T_REAL const i_mat[TL_N_SCS][TL_N_MDS],
                              TL_T_REAL       o_int[TL_N_QTS][TL_N_MDS][TL_N_CRS] ) {
      // project element's sub-cell solution at faces to DG solution
      linalg::Matrix::matMulB0FusedAC( TL_N_CRS,
                                       TL_N_QTS, TL_N_MDS, TL_N_SFS, // m, n, k
                                       TL_N_SFS, TL_N_MDS, TL_N_MDS, // ldA, ldB, ldC
                                       i_scFluxes[0][0],             // A
                                       i_mat[0],                     // B
                                       o_int[0][0] );                // C
    }

    /**
     * Computes the extrema of a sub-cell solution.
     *
     * @param i_dofsSc DOFs of the sub-cell solution.
     * @param o_min will be set to minima of the sub-cell solution.
     * @param o_max will be set to maxima of the sub-cell solution.
     *
     * @paramt TL_T_REAL floating point type.
     **/
    template < typename TL_T_REAL >
    static void scExtrema( TL_T_REAL const i_dofsSc[TL_N_QTS][TL_N_SCS][TL_N_CRS],
                           TL_T_REAL       o_min[TL_N_QTS][TL_N_CRS],
                           TL_T_REAL       o_max[TL_N_QTS][TL_N_CRS] ) {
      // init extrema
      for( unsigned short l_qt = 0; l_qt < TL_N_QTS; l_qt++ ) {
        for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {
          o_min[l_qt][l_cr] = std::numeric_limits< TL_T_REAL >::max();
          o_max[l_qt][l_cr] = std::numeric_limits< TL_T_REAL >::lowest();
        }
      }

      // determine extrema
      for( unsigned short l_sc = 0; l_sc < TL_N_SCS; l_sc++ ) {
        for( unsigned short l_qt = 0; l_qt < TL_N_QTS; l_qt++ ) {
          for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {
            o_min[l_qt][l_cr] = std::min( o_min[l_qt][l_cr], i_dofsSc[l_qt][l_sc][l_cr] );
            o_max[l_qt][l_cr] = std::max( o_max[l_qt][l_cr], i_dofsSc[l_qt][l_sc][l_cr] );
          }
        }
      }
    }

    /**
     * Computes the extrema of a DG solution.
     *
     * @param i_dofsDg degrees of freedom of the DG solution.
     * @param i_matScatter scatter matrix.
     * @param o_subcell will be set to the sub-cell solution.
     * @param o_min will be set to minima of the corresponding sub-cell solution.
     * @param o_max will be set to maxima of the corresponding sub-cell solution.
     *
     * @paramt TL_T_REAL floating point type.
     **/
    template < typename TL_T_REAL >
    static void dgExtremaVanilla( TL_T_REAL const i_dofsDg[TL_N_QTS][TL_N_MDS][TL_N_CRS],
                                  TL_T_REAL const i_matScatter[TL_N_MDS][TL_N_SCS],
                                  TL_T_REAL       o_subcell[TL_N_QTS][TL_N_SCS][TL_N_CRS],
                                  TL_T_REAL       o_min[TL_N_QTS][TL_N_CRS],
                                  TL_T_REAL       o_max[TL_N_QTS][TL_N_CRS] ) {
      // determine sub-cell solution
      scatterVanilla( i_dofsDg,
                      i_matScatter,
                      o_subcell );

      // determine extrema
      scExtrema( o_subcell,
                 o_min,
                 o_max );
    }
};

#endif
