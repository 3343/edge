/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *         Alexander Heinecke (alexander.heinecke AT intel.com)
 *
 * @section LICENSE
 * Copyright (c) 2017-2018, Regents of the University of California
 * Copyright (c) 2018, Intel Corporation
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
#ifndef EDGE_SC_KERNELS_HPP
#define EDGE_SC_KERNELS_HPP

#include "constants.hpp"

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
    //! number of faces
    static unsigned short const TL_N_FAS = C_ENT[TL_T_EL].N_FACES;

    //! number of vertices per face
    static unsigned short const TL_N_VES_FA = C_ENT[TL_T_EL].N_FACE_VERTICES;

    //! number of sub-faces per element face
    static unsigned short const TL_N_SFS = CE_N_SUB_FACES( TL_T_EL, TL_O_SP );

    //! number of sub-cells
    static unsigned short const TL_N_SCS = CE_N_SUB_CELLS( TL_T_EL, TL_O_SP );

    //! number of modes
    static unsigned short const TL_N_MDS = CE_N_ELEMENT_MODES( TL_T_EL, TL_O_SP );

    //! id of the matrix kernel group
    static unsigned short const MM_GR = static_cast< unsigned short >( t_mm::SUB_CELL );

  public:
    /**
     * Applies the scatter matrix.
     * Scatter: Projects DG solution to sub-cell solution.
     *
     * @param i_mm Fused dense JITed matrix kernels 
     * @param i_dofsDg DOFs of DG solution.
     * @param i_mat scatter matrix.
     * @param o_scDofs will be set to degrees of freedom of the sub-cell solution.
     *
     * @paramt TL_T_REAL floating point type.
     * @paramt TL_T_MM type of the matrix-matrix multiplication kernels.
     **/
    template< typename TL_T_REAL,
              typename TL_T_MM >
    static void scatter( TL_T_MM   const & i_mm,
                         TL_T_REAL const   i_dofsDg[TL_N_QTS][TL_N_MDS][TL_N_CRS],
                         TL_T_REAL const   i_mat[TL_N_MDS][TL_N_SCS],
                         TL_T_REAL         o_scDofs[TL_N_QTS][TL_N_SCS][TL_N_CRS] ) {
#if defined(PP_T_KERNELS_XSMM_DENSE_SINGLE)
      i_mm.m_kernels[MM_GR][0]( i_mat[0], i_dofsDg[0][0], o_scDofs[0][0] );
#else
      i_mm.m_kernels[MM_GR][0]( i_dofsDg[0][0], i_mat[0], o_scDofs[0][0] );
#endif
    }

    /**
     * Applies the scatter matrix for sub-cells at the element's face.
     * Scatter: Projects DG solution to sub-cell solution.
     *
     * @param i_dofsDg DOFs of DG solution.
     * @param i_mat scatter matrix.
     * @param o_scDofs will be set to degrees of freedom of the sub-cell solution.
     *
     * @paramt TL_T_REAL floating point type.
     * @paramt TL_T_MM type of the matrix-matrix multiplication kernels.
     **/
    template< typename TL_T_REAL,
              typename TL_T_MM >
    static void scatterFa( TL_T_MM   const & i_mm,
                           TL_T_REAL const   i_dofsDg[TL_N_QTS][TL_N_MDS][TL_N_CRS],
                           TL_T_REAL const   i_mat[TL_N_MDS][TL_N_SFS],
                           TL_T_REAL         o_scDofs[TL_N_QTS][TL_N_SFS][TL_N_CRS] ) {
#if defined(PP_T_KERNELS_XSMM_DENSE_SINGLE)
      i_mm.m_kernels[MM_GR][1]( i_mat[0], i_dofsDg[0][0], o_scDofs[0][0] );
#else
      i_mm.m_kernels[MM_GR][1]( i_dofsDg[0][0], i_mat[0], o_scDofs[0][0] );
#endif
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
     * @paramt TL_T_MM type of the matrix-matrix multiplication kernels.
     **/
    template< typename TL_T_REAL,
              typename TL_T_MM >
    static void scatterReplace( TL_T_MM   const & i_mm,
                                TL_T_REAL const   i_dofsDg[TL_N_QTS][TL_N_MDS][TL_N_CRS],
                                TL_T_REAL const   i_mat[TL_N_MDS][TL_N_SCS],
                                TL_T_REAL const   i_dofsSc[TL_N_QTS][TL_N_SCS][TL_N_CRS],
                                bool      const   i_admP[TL_N_CRS],
                                TL_T_REAL         o_dofsSc[TL_N_QTS][TL_N_SCS][TL_N_CRS] ) {
      // project DG to sub-cell solution
      scatter( i_mm,
               i_dofsDg,
               i_mat,
               o_dofsSc );

      // overwrite project sub-cell solution with stored sub-cell solution
      for( unsigned short l_qt = 0; l_qt < TL_N_QTS; l_qt++ ) {
        for( unsigned short l_sc = 0; l_sc < TL_N_SCS; l_sc++ ) {
#pragma omp simd 
          for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {  
              o_dofsSc[l_qt][l_sc][l_cr] = ( i_admP[l_cr] == false ) ? i_dofsSc[l_qt][l_sc][l_cr] : o_dofsSc[l_qt][l_sc][l_cr];
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
     * @paramt TL_T_MM type of the matrix-matrix multiplication kernels.
     **/
    template< typename TL_T_REAL,
              typename TL_T_MM >
    static void scatterReplaceFa( TL_T_MM   const & i_mm,
                                  TL_T_REAL const   i_dofsDg[TL_N_QTS][TL_N_MDS][TL_N_CRS],
                                  TL_T_REAL const   i_mat[TL_N_MDS][TL_N_SFS],
                                  TL_T_REAL const   i_dofsSc[TL_N_QTS][TL_N_SFS][TL_N_CRS],
                                  bool      const   i_admP[TL_N_CRS],
                                  TL_T_REAL         o_dofsSc[TL_N_QTS][TL_N_SFS][TL_N_CRS] ) {
      // project DG to sub-cell solution
      scatterFa( i_mm,
                 i_dofsDg,
                 i_mat,
                 o_dofsSc );

      // overwrite project sub-cell solution with stored sub-cell solution
      if( i_admP != nullptr ) {
        for( unsigned short l_qt = 0; l_qt < TL_N_QTS; l_qt++ ) {
          for( unsigned short l_sf = 0; l_sf < TL_N_SFS; l_sf++ ) {
#pragma omp simd
            for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {
              o_dofsSc[l_qt][l_sf][l_cr] = ( i_admP[l_cr] == false ) ? i_dofsSc[l_qt][l_sf][l_cr] : o_dofsSc[l_qt][l_sf][l_cr];
            }
          }
        }
      }
    }

    /**
     * Applies the gather matrix.
     * Gather: Projects sub-cell solution to DG-solution.
     *
     * @param i_scDofs DOFs of sub-cell solution.
     * @param i_mat gather matrix.
     * @param o_dgDofs will be set to degrees of freedom of DG solution.
     *
     * @paramt TL_T_REAL floating point type.
     * @paramt TL_T_MM type of the matrix-matrix multiplication kernels.
     **/
    template< typename TL_T_REAL,
              typename TL_T_MM >
    static void gather( TL_T_MM   const &i_mm,
                        TL_T_REAL const i_scDofs[TL_N_QTS][TL_N_SCS][TL_N_CRS],
                        TL_T_REAL const i_mat[TL_N_SCS][TL_N_MDS],
                        TL_T_REAL       o_dgDofs[TL_N_QTS][TL_N_MDS][TL_N_CRS] ) {
      // project element's sub-cell solution to DG solution
#if defined(PP_T_KERNELS_XSMM_DENSE_SINGLE)
      i_mm.m_kernels[MM_GR][2]( i_mat[0], i_scDofs[0][0], o_dgDofs[0][0] );
#else
      i_mm.m_kernels[MM_GR][2]( i_scDofs[0][0], i_mat[0], o_dgDofs[0][0] );
#endif
    }

    /**
     * Gathers the DOFs of sub-cells at the surface of the element and stores them.
     * The sub-cells are stored from the view of the face-adjacent element and account for possible face-vertex combinations.
     *
     * @param i_scDofs sub-cell DOFs of the element.
     * @param i_faSfSc sub-cells adjacent to the sub-faces of a face.
     * @param i_scDgAd sub-cell reordering, based for adjacent element.
     * @param i_vIdElFaEl vertex combis.
     * @param o_scSurfDofs will be set to the surface DOFs, unless a nullptr is set.
     * @param i_adm if not nullptr, only fused non-admissible fused runs are updated
     *
     * @paramt TL_T_REAL floating point precision.
     **/
    template< typename TL_T_REAL >
    static void gatherSurfDofs( TL_T_REAL      const    i_scDofs[TL_N_QTS][TL_N_SCS][TL_N_CRS],
                                unsigned short const    i_faSfSc[TL_N_FAS][TL_N_SFS],
                                unsigned short const    i_scDgAd[TL_N_VES_FA][TL_N_SFS],
                                unsigned short const    i_vIdElFaEl[TL_N_FAS],
                                TL_T_REAL            (* o_scSurfDofs [TL_N_FAS])[TL_N_QTS][TL_N_SFS][TL_N_CRS],
                                bool           const    i_adm[TL_N_CRS] ) {
      // iterate over faces
      for( unsigned short l_fa = 0; l_fa < TL_N_FAS; l_fa++ ) {
        // only continue if the DOFs are required
        if( o_scSurfDofs[l_fa] != nullptr ) {
          // get vertex combination
          unsigned short l_vId = i_vIdElFaEl[l_fa];

          // iterate over quantities
          for( unsigned short l_qt = 0; l_qt < TL_N_QTS; l_qt++ ) {
            // iterate over sub-faces
            for( unsigned short l_sf = 0; l_sf < TL_N_SFS; l_sf++ ) {
              // determine sub-face from the view of the adjacent element
              unsigned short l_sfRe = i_scDgAd[l_vId][l_sf];

              // get id of copy sub-cell (faces as bridge)
              unsigned short l_sc = i_faSfSc[l_fa][l_sfRe];

              // copy over the sub-cell DOFs
#pragma omp simd
              for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {
                (*o_scSurfDofs[l_fa])[l_qt][l_sf][l_cr] = ( i_adm[l_cr] == false ) ? i_scDofs[l_qt][l_sc][l_cr] : (*o_scSurfDofs[l_fa])[l_qt][l_sf][l_cr];
              }
            }
          }

        }
      }
    }

    /**
     * Applies the sfInt operator.
     * sfInt: Integration of piecewise-constant fluxes of a DG-face.
     *
     * @param i_scFluxes fluxes of the sub-cell solution at the face.
     * @param i_mat sfInt matrix.
     * @param o_int will be set to surface integral.
     *
     * @paramt TL_T_REAL floating point type.
     * @paramt TL_T_MM type of the matrix-matrix multiplication kernels.
     **/
    template< typename TL_T_REAL,
              typename TL_T_MM >
    static void sfInt( TL_T_MM   const &i_mm,
                       TL_T_REAL const i_scFluxes[TL_N_QTS][TL_N_SFS][TL_N_CRS],
                       TL_T_REAL const i_mat[TL_N_SCS][TL_N_MDS],
                       TL_T_REAL       o_int[TL_N_QTS][TL_N_MDS][TL_N_CRS] ) {
      // project element's sub-cell solution at faces to DG solution
#if defined(PP_T_KERNELS_XSMM_DENSE_SINGLE)
      i_mm.m_kernels[MM_GR][3]( i_mat[0], i_scFluxes[0][0], o_int[0][0] );
#else
      i_mm.m_kernels[MM_GR][3]( i_scFluxes[0][0], i_mat[0], o_int[0][0] );
#endif
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
    template< typename TL_T_REAL >
    static void scExtrema( TL_T_REAL const i_dofsSc[TL_N_QTS][TL_N_SCS][TL_N_CRS],
                           TL_T_REAL       o_min[TL_N_QTS][TL_N_CRS],
                           TL_T_REAL       o_max[TL_N_QTS][TL_N_CRS] ) {
      // init extrema
      for( unsigned short l_qt = 0; l_qt < TL_N_QTS; l_qt++ ) {
#pragma omp simd
        for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ ) {
          o_min[l_qt][l_cr] = std::numeric_limits< TL_T_REAL >::max();
          o_max[l_qt][l_cr] = std::numeric_limits< TL_T_REAL >::lowest();
        }
      }

      // determine extrema
      for( unsigned short l_sc = 0; l_sc < TL_N_SCS; l_sc++ ) {
        for( unsigned short l_qt = 0; l_qt < TL_N_QTS; l_qt++ ) {
#pragma omp simd
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
     * @paramt TL_T_MM type of the matrix-matrix multiplication kernels.
     **/
    template< typename TL_T_REAL,
              typename TL_T_MM >
    static void dgExtrema( TL_T_MM   const &i_mm,
                           TL_T_REAL const i_dofsDg[TL_N_QTS][TL_N_MDS][TL_N_CRS],
                           TL_T_REAL const i_matScatter[TL_N_MDS][TL_N_SCS],
                           TL_T_REAL       o_subcell[TL_N_QTS][TL_N_SCS][TL_N_CRS],
                           TL_T_REAL       o_min[TL_N_QTS][TL_N_CRS],
                           TL_T_REAL       o_max[TL_N_QTS][TL_N_CRS] ) {
      // determine sub-cell solution
      scatter( i_mm,
               i_dofsDg,
               i_matScatter,
               o_subcell );

      // determine extrema
      scExtrema( o_subcell,
                 o_min,
                 o_max );
    }
};

#endif
