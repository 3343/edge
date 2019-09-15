/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *         Alexander Heinecke (alexander.heinecke AT intel.com)
 *
 * @section LICENSE
 * Copyright (c) 2016-2017, Regents of the University of California
 * Copyright (c) 2016, Intel Corporation
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
 * Compile time constants for ADER-DG simulations.
 **/

#if defined PP_T_ELEMENTS_TRIA3
#define PP_T_BASIS_HIERARCHICAL
#elif defined PP_T_ELEMENTS_TET4
#define PP_T_BASIS_HIERARCHICAL
#endif

#define PP_N_GLOBAL_SHARED_1 1

// derive number of flux matrices
const unsigned short N_FLUXN_MATRICES = CE_N_FLUXN_MATRICES( T_SDISC.ELEMENT );
const unsigned short N_FACE_QUAD_POINTS = CE_N_FACE_QUAD_POINTS( T_SDISC.ELEMENT, ORDER );

const unsigned short N_FACE_MODES = CE_N_ELEMENT_MODES( T_SDISC.FACE, ORDER );

// human-readable access of LTS-types
typedef enum {
  EL_INT_EQ = 0, // GTS tInt "="
  EL_INT_LT = 1, // LTS tInt "<"
  EL_INT_GT = 2  // LTS tInt ">"
} t_ltsEl;
const int_spType C_LTS_EL[3] = { // data offered by the element
  //   xxxxxxxxxxxxxxxoxxxxxxxxxxxxxxxx <- time bits
  //        56|     48|     40|     32|     24|     16|      8|      0|
  // 0b0000000000000000000000000000000100000000000000000000000000000000 // GTS: tInt
  4294967296,
  // 0b0000000000000000000000000000001000000000000000000000000000000000 // LTS "<"
  8589934592,
  // 0b0000000000000000000000000000010000000000000000000000000000000000 // LTS ">"
  17179869184
};

// the following ADER-LTS types are only defined for 6 faces
static_assert( C_ENT[T_SDISC.ELEMENT].N_FACES <= 6, "more than 6 face not supported" );
typedef enum {
  AD_EQ = 0, // adjacent element has an equal time step
  AD_LT = 1, // element has a time step less than the adjacent one ("< relation")
  AD_GT = 2  // element has time step greater than the adjacent one ("> relation")
} t_ltsAd;
const int_spType C_LTS_AD[6][3] = { // LTS relations with adjacent elements
  {
    // elFaEl #0 or faEl #0
    //   xxxxxxxxxxxxxxxoxxxxxxxxxxxxxxxx <- time bits
    //        56|     48|     40|     32|     24|     16|      8|      0|
    // 0b0000000000000000000000000001000000000000000000000000000000000000 // == relation
    68719476736,
    // 0b0000000000000000000000000010000000000000000000000000000000000000 // < relation
    137438953472,
    // 0b0000000000000000000000000100000000000000000000000000000000000000 // > relation
    274877906944
  },
  {
    // elFaEl #1 or faEl #1
    //   xxxxxxxxxxxxxxxoxxxxxxxxxxxxxxxx <- time bits
    //        56|     48|     40|     32|     24|     16|      8|      0|
    // 0b0000000000000000000000001000000000000000000000000000000000000000 // == relation
    549755813888,
    // 0b0000000000000000000000010000000000000000000000000000000000000000 // < relation
    1099511627776,
    // 0b0000000000000000000000100000000000000000000000000000000000000000 // > relation
    2199023255552
  },
  {
    // elFaEl #2
    //   xxxxxxxxxxxxxxxoxxxxxxxxxxxxxxxx <- time bits
    //        56|     48|     40|     32|     24|     16|      8|      0|
    // 0b0000000000000000000001000000000000000000000000000000000000000000 // == relation
    4398046511104,
    // 0b0000000000000000000010000000000000000000000000000000000000000000 // < relation
    8796093022208,
    // 0b0000000000000000000100000000000000000000000000000000000000000000 // > relation
    17592186044416
  },
  {
    // elFaEl #3
    //   xxxxxxxxxxxxxxxoxxxxxxxxxxxxxxxx <- time bits
    //        56|     48|     40|     32|     24|     16|      8|      0|
    // 0b0000000000000000001000000000000000000000000000000000000000000000 // == relation
    35184372088832,
    // 0b0000000000000000010000000000000000000000000000000000000000000000 // < relation
    70368744177664,
    // 0b0000000000000000100000000000000000000000000000000000000000000000 // > relation
    140737488355328
  },
  {
    // elFaEl #4
    //   xxxxxxxxxxxxxxxoxxxxxxxxxxxxxxxx <- time bits
    //        56|     48|     40|     32|     24|     16|      8|      0|
    // 0b0000000000000001000000000000000000000000000000000000000000000000 // == relation
    281474976710656,
    // 0b0000000000000010000000000000000000000000000000000000000000000000 // < relation
    562949953421312,
    // 0b0000000000000100000000000000000000000000000000000000000000000000 // > relation
    1125899906842624
  },
  {
    // elFaEl #5
    //   xxxxxxxxxxxxxxxoxxxxxxxxxxxxxxxx <- time bits
    //        56|     48|     40|     32|     24|     16|      8|      0|
    // 0b0000000000001000000000000000000000000000000000000000000000000000 // == relation
    2251799813685248,
    // 0b0000000000010000000000000000000000000000000000000000000000000000 // < relation
    4503599627370496,
    // 0b0000000000100000000000000000000000000000000000000000000000000000 // > relation
    9007199254740992
  }
};

typedef struct {
  //! inverse mass matrix
  real_base massI[N_ELEMENT_MODES];
} t_dgMat;

typedef struct {
  //! 1D quadrature points
  real_mesh ptsLine[ CE_N_ELEMENT_MODES( LINE, ORDER ) ];

  //! weights of 1D quadrature points
  real_base weightsLine[ CE_N_ELEMENT_MODES( LINE, ORDER ) ];

  //! quadrature points in reference coordinates at the face's quad points
  real_mesh ptsFaces[ (CE_N_FACE_VERTEX_OPTS(T_SDISC.ELEMENT)+1) * C_ENT[T_SDISC.ELEMENT].N_FACES ]
                    [ CE_N_FACE_QUAD_POINTS( T_SDISC.ELEMENT, ORDER )                             ]
                    [ C_ENT[T_SDISC.ELEMENT].N_DIM                                                ];

  //! weights of the quadrature points (sum is equal to the surface of the face's ref element)
  real_base weightsFaces[ CE_N_FACE_QUAD_POINTS( T_SDISC.ELEMENT, ORDER ) ];

  //! evaluated basis functions at the face's quad points
  real_base basisFaces[ (CE_N_FACE_VERTEX_OPTS(T_SDISC.ELEMENT)+1) * C_ENT[T_SDISC.ELEMENT].N_FACES ]
                      [ CE_N_FACE_QUAD_POINTS( T_SDISC.ELEMENT, ORDER )                             ]
                      [ CE_N_ELEMENT_MODES( T_SDISC.ELEMENT, ORDER )                                ];
} t_dgQuadEval;

typedef struct {
  t_dgMat mat;
  t_dgQuadEval quadEval;
} t_dg;

typedef t_dg t_globalShared1;
