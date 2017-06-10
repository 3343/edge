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
 * Support for viscosity.
 **/

#ifndef VISCOSITY_HPP
#define VISCOSITY_HPP

#include "../common.hpp"
#include "linalg/Matrix.h"

namespace edge {
  namespace elastic {
    namespace solvers {
      template< unsigned short TL_N_DIM >
      class Viscosity;

      template<>
      class Viscosity< 2 >;

      template<>
      class Viscosity< 3 >;
    }
  }
}

/**
 * Viscosity part of the two-dimensional flux solvers.
 **/
template<>
class edge::elastic::solvers::Viscosity< 2 > {
  private:
    /**
     * Harten's entropy fix in two dimensions.
     *
     * Reference: Finite Volume Methods for Hyperbolic Problems
     *            Randy LeVeque, 2002
     *            Formula: (15.53) & (15.54) page 326
     *
     * Remark: The paramter delta in literature is already set to the maximum p-wave speed
     *         of both adjacent elements.
     *
     * @param i_lamL Lame parameter lambda of the left element.
     * @param i_lamR Lame parameter lambda of the right element.
     * @param i_muL Lame parameter mu of the left element.
     * @param i_muR Lame parameter mu of the right element.
     * @param i_rhoL density rho of the left element.
     * @param i_rhoR density rho of the right element.
     * @param i_scale scaling of the parameter delta (set to c_p if scale is 1).
     * @param o_entFixHarten will be set to flux solver carrying the contribution of Harten's entropy fix in face-aligned coordinates.
     *
     * @paramt TL_T_REAL floating point precision.
     **/
    template< typename TL_T_REAL >
    static void entFixHarten( TL_T_REAL i_lamL,
                              TL_T_REAL i_lamR,
                              TL_T_REAL i_muL,
                              TL_T_REAL i_muR,
                              TL_T_REAL i_rhoL,
                              TL_T_REAL i_rhoR,
                              TL_T_REAL i_scale,
                              TL_T_REAL o_entFixHarten[5][5] ) {
#include "../generated/EntFixHarten2D.inc"
    }

  public:
    /**
     * Harten's entropy fix in two dimensions.
     *
     * Reference: Finite Volume Methods for Hyperbolic Problems
     *            Randy LeVeque, 2002
     *            Formula: (15.53) & (15.54) page 326
     *
     * Remark: The paramter delta in literature is already set to the maximum p-wave speed
     *         of both adjacent elements.
     *
     * @param i_lamL Lame parameter lambda of the left element.
     * @param i_lamR Lame parameter lambda of the right element.
     * @param i_muL Lame parameter mu of the left element.
     * @param i_muR Lame parameter mu of the right element.
     * @param i_rhoL density rho of the left element.
     * @param i_rhoR density rho of the right element.
     * @param i_scale scaling of the parameter delta (set to c_p is scale is 1).
     * @param i_n outer point normal of the face (length 1).
     * @param o_entFixHarten will be set to flux solver carrying the contribution of Harten's entropy fix, including the corresponding rotation to and from face-aligned coordinates.
     *
     * @paramt TL_T_REAL_MESH floating point precision of mesh-related data.
     * @paramt TL_T_REAL_COMP floating point precision of computational data.
     **/
    template< typename TL_T_REAL_MESH, typename TL_T_REAL_COMP >
    static void entFixHarten( TL_T_REAL_COMP       i_lamL,
                              TL_T_REAL_COMP       i_lamR,
                              TL_T_REAL_COMP       i_muL,
                              TL_T_REAL_COMP       i_muR,
                              TL_T_REAL_COMP       i_rhoL,
                              TL_T_REAL_COMP       i_rhoR,
                              TL_T_REAL_COMP       i_scale,
                              TL_T_REAL_MESH const i_n[2],
                              TL_T_REAL_COMP       o_entFixHarten[5][5] ) {
      // get the trafos
      TL_T_REAL_MESH l_t[5][5];
      TL_T_REAL_MESH l_tm1[5][5];

      elastic::common::setupTrafo2d(    i_n[0], i_n[1], l_t );
      elastic::common::setupTrafoInv2d( i_n[0], i_n[1], l_tm1 );

      // get the face-aligned llf entropy fix
      TL_T_REAL_COMP l_entFix[5][5];
      entFixHarten( i_lamL, i_lamR,
                    i_muL,  i_muR,
                    i_rhoL, i_rhoR,
                    i_scale,
                    l_entFix );

      // do the matrix mults
      TL_T_REAL_COMP l_tmp[5][5];
      linalg::Matrix::matMulB0( 5, 5, 5,
                                l_t[0], l_entFix[0], l_tmp[0] );
      linalg::Matrix::matMulB0( 5, 5, 5,
                                l_tmp[0], l_tm1[0], o_entFixHarten[0] );
    }

};

/**
 * Viscosity part of the three-dimensional flux solvers.
 *
 **/
template<>
class edge::elastic::solvers::Viscosity< 3 > {
  private:
    /**
     * Harten's entropy fix in three dimensions.
     *
     * Reference: Finite Volume Methods for Hyperbolic Problems
     *            Randy LeVeque, 2002
     *            Formula: (15.53) & (15.54) page 326
     *
     * Remark: The paramter delta in literature is already set to the maximum
     *         s-wave or p-wave speed p-wave speed of both adjacent elements.
     *
     * @param i_lamL Lame parameter lambda of the left element.
     * @param i_lamR Lame parameter lambda of the right element.
     * @param i_muL Lame parameter mu of the left element.
     * @param i_muR Lame parameter mu of the right element.
     * @param i_rhoL density rho of the left element.
     * @param i_rhoR density rho of the right element.
     * @param i_scale scaling of the parameter delta (set to c_p is scale is 1).
     * @param o_entFixHarten will be set to flux solver carrying the contribution of Harten's entropy fix in face-aligned coordinates.
     *
     * @paramt TL_T_REAL floating point precision.
     **/
    template< typename TL_T_REAL >
    static void entFixHarten( TL_T_REAL i_lamL,
                              TL_T_REAL i_lamR,
                              TL_T_REAL i_muL,
                              TL_T_REAL i_muR,
                              TL_T_REAL i_rhoL,
                              TL_T_REAL i_rhoR,
                              TL_T_REAL i_scale,
                              TL_T_REAL o_entFixHarten[9][9] ) {
#include "../generated/EntFixHarten3D.inc"
    }

  public:
    /**
     * Harten's entropy fix in three dimensions.
     *
     * Reference: Finite Volume Methods for Hyperbolic Problems
     *            Randy LeVeque, 2002
     *            Formula: (15.53) & (15.54) page 326
     *
     * Remark: The paramter delta in literature is already set to the maximum p-wave speed
     *         of both adjacent elements.
     *
     * @param i_lamL Lame parameter lambda of the left element.
     * @param i_lamR Lame parameter lambda of the right element.
     * @param i_muL Lame parameter mu of the left element.
     * @param i_muR Lame parameter mu of the right element.
     * @param i_rhoL density rho of the left element.
     * @param i_rhoR density rho of the right element.
     * @param i_scale scaling of the parameter delta (set to c_p/c_s if scale is 1).
     * @param i_n outer point normal of the face (length 1).
     * @param i_s first tangent of the face (normal to normal and second tangent, length 1).
     * @parma i_t second tangent of the facce (normal to nromal and first tangent, length 1).
     * @param o_entFixHarten will be set to flux solver carrying the contribution of Harten's entropy fix, including the corresponding rotation to and from face-aligned coordinates.
     *
     * @paramt TL_T_REAL_MESH floating point precision of mesh-related data.
     * @paramt TL_T_REAL_COMP floating point precision of computational data.
     **/
    template< typename TL_T_REAL_MESH, typename TL_T_REAL_COMP >
    static void entFixHarten( TL_T_REAL_COMP       i_lamL,
                              TL_T_REAL_COMP       i_lamR,
                              TL_T_REAL_COMP       i_muL,
                              TL_T_REAL_COMP       i_muR,
                              TL_T_REAL_COMP       i_rhoL,
                              TL_T_REAL_COMP       i_rhoR,
                              TL_T_REAL_COMP       i_scale,
                              TL_T_REAL_MESH const i_n[3],
                              TL_T_REAL_MESH const i_s[3],
                              TL_T_REAL_MESH const i_t[3],
                              TL_T_REAL_COMP       o_entFixHarten[9][9] ) {
      // get the trafos
      TL_T_REAL_MESH l_t[9][9];
      TL_T_REAL_MESH l_tm1[9][9];

      elastic::common::setupTrafo3d(    i_n[0], i_n[1], i_n[2],
                                        i_s[0], i_s[1], i_s[2],
                                        i_t[0], i_t[1], i_t[2],
                                        l_t );
      elastic::common::setupTrafoInv3d( i_n[0], i_n[1], i_n[2],
                                        i_s[0], i_s[1], i_s[2],
                                        i_t[0], i_t[1], i_t[2],
                                        l_tm1 );

      // get the face-aligned llf entropy fix
      TL_T_REAL_COMP l_entFix[9][9];
      entFixHarten( i_lamL, i_lamR,
                    i_muL,  i_muR,
                    i_rhoL, i_rhoR,
                    i_scale,
                    l_entFix );

      // do the matrix mults
      TL_T_REAL_COMP l_tmp[9][9];
      linalg::Matrix::matMulB0( 9, 9, 9,
                                l_t[0], l_entFix[0], l_tmp[0] );
      linalg::Matrix::matMulB0( 9, 9, 9,
                                l_tmp[0], l_tm1[0], o_entFixHarten[0] );
    }

};

#endif
