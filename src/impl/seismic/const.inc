/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2019, Alexander Breuer
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
 * Compile time constants seismic simulations.
 **/

/**
 * Gets the number of elastic quantities.
 *
 * @param i_nDis number of dimensions.
 * @return number of quantities.
 **/
constexpr unsigned short CE_N_QTS_E( unsigned short i_nDis ) {
  return (i_nDis == 2) ? 5 : 9;
}

/**
 * Gets the number of quantities per relaxation mechanism.
 *
 * @param i_nDis number of dimensions.
 * @return number of quantities.
 **/
constexpr unsigned short CE_N_QTS_M( unsigned short i_nDis ) {
  return (i_nDis == 2) ? 3 : 6;
}

/**
 * Gets the number of entries in a dense elastic star matrix.
 *
 * @param i_nDis number of dimensions.
 * @return number of entries.
 **/
constexpr unsigned short CE_N_ENS_STAR_E_DE( unsigned short i_nDis ) {
  return CE_N_QTS_E(i_nDis) * CE_N_QTS_E(i_nDis);
}

/**
 * Gets the number of non-zeros in a sparse elastic star matrix.
 *
 * @param i_nDis number of dimensions.
 * @return number of non-zeros.
 **/
constexpr unsigned short CE_N_ENS_STAR_E_SP( unsigned short i_nDis ) {
  return (i_nDis == 2) ? 10 : 24;
}

/**
 * Gets the number of entries in a dense anelastic star matrix.
 *
 * @param i_nDis number of dimensions.
 * @return number of entries.
 **/
constexpr unsigned short CE_N_ENS_STAR_A_DE( unsigned short i_nDis ) {
  return i_nDis * CE_N_QTS_M( i_nDis );
}

/**
 * Gets the number of non-zeros in a sparse anelastic star matrix.
 *
 * @param i_nDis number of dimensions.
 * @return number of non-zeros.
 **/
constexpr unsigned short CE_N_ENS_STAR_A_SP( unsigned short i_nDis ) {
  return (i_nDis == 2) ? 4 : 9;
}

/**
 * Gets the number of entries in a dense anelastic source matrix.
 *
 * @param i_nDis number of dimensions.
 * @return number of entries.
 **/
constexpr unsigned short CE_N_ENS_SRC_A_DE( unsigned short i_nDis ) {
  return CE_N_QTS_M( i_nDis ) * CE_N_QTS_M( i_nDis );
}

/**
 * Gets the number of non-zeros in a sparse anelastic source matrix.
 *
 * @param i_nDis number of dimensions.
 * @return number of non-zeros.
 **/
constexpr unsigned short CE_N_ENS_SRC_A_SP( unsigned short i_nDis ) {
  return (i_nDis == 2) ? 5 : 12;
}

/**
 * Gets the number of entries in a dense elastic flux solver.
 *
 * @param i_nDis number of dimensions.
 * @return number of entries.
 **/
constexpr unsigned short CE_N_ENS_FS_E_DE( unsigned short i_nDis ) {
  return CE_N_QTS_E( i_nDis ) * CE_N_QTS_E( i_nDis );
}

/**
 * Gets the number of entries in a dense anelastic flux solver.
 *
 * @param i_nDis number of dimensions.
 * @return number of entries.
 **/
constexpr unsigned short CE_N_ENS_FS_A_DE( unsigned short i_nDis ) {
  return CE_N_QTS_E( i_nDis ) * CE_N_QTS_M( i_nDis );
}

// number of entries in the control flow
const unsigned short N_ENTRIES_CONTROL_FLOW=9;

// number of relaxation mechanisms
const unsigned short N_RELAXATION_MECHANISMS=PP_N_RELAXATION_MECHANISMS;

#if PP_ORDER > 9
#error not supported
#endif
// CFL-number
const double SCALE_CFL_ORDER[9] = { 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9 };

const double SCALE_CFL=(1.0/N_DIM)*SCALE_CFL_ORDER[PP_ORDER-1];

/*
 * Bitmasks for elastics.
 **/
static_assert( sizeof(int_spType) >= 4, "size of entity type too small" );
typedef enum: int_spType {
  FREE_SURFACE = 101,
  OUTFLOW      = 105,
  PERIODIC     = 106,
  RUPTURE      = 201,
  SOURCE       = 65536 // 0b0000000000000000000000000000000000000000000000010000000000000000
} t_spTypeElastic;

// elastic quantities
const unsigned short N_QUANTITIES = CE_N_QTS_E(PP_N_DIM);

// number of anelastic quantities
const unsigned short N_QUANTITIES_A = N_RELAXATION_MECHANISMS * CE_N_QTS_M(PP_N_DIM);

typedef struct {
  // density rho
  real_base rho;
  // lame parameter lambda
  real_base lam;
  // lame parameter mu
  real_base mu;
  // q-factor qp
  real_base qp;
  // q-factor qs
  real_base qs;
} t_bgPars;

/**
 * Scratch memory (per thread)
 **/
#define PP_SCRATCH_MEMORY
#define PP_SCRATCH_MEMORY_HBW
#define PP_SCRATCH_MEMORY_HUGE
struct scratchMem {
  // temporary results
  real_base tRes[2][N_QUANTITIES][N_ELEMENT_MODES][N_CRUNS] __attribute__ ((aligned (ALIGNMENT.ELEMENT_MODES.PRIVATE)));
  // derivative buffer
  real_base dBuf[ORDER][N_QUANTITIES][N_ELEMENT_MODES][N_CRUNS] __attribute__ ((aligned (ALIGNMENT.ELEMENT_MODES.PRIVATE)));
  // scratch memory for the surface integration
  real_base tResSurf[2][N_QUANTITIES][N_FACE_MODES][N_CRUNS] __attribute__ ((aligned (ALIGNMENT.FACE_MODES.PRIVATE)));
};
typedef scratchMem t_scratchMem;