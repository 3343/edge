/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2016-2017, Regents of the University of California
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
 * Compile time constants for the advection equation.
 **/

// advection equation performs two steps per time step
// 0) time prediction
// 1) update
// 2) limit
const unsigned short N_STEPS_PER_UPDATE=3;
const unsigned short N_ENTRIES_CONTROL_FLOW=7;

#if PP_ORDER > 9
#error not supported
#endif
// CFL-numbers based on (experimental) 1D-stability
const double SCALE_CFL_ORDER[9] = { 0.45, 0.9, 0.8, 0.7, 0.6, 0.5, 0.45, 0.4, 0.4};

#if defined PP_T_ELEMENTS_LINE
const double SCALE_CFL=SCALE_CFL_ORDER[PP_ORDER-1];
#else
const double SCALE_CFL=0.5*SCALE_CFL_ORDER[PP_ORDER-1];
#endif

const unsigned short N_QUANTITIES=1;

/**
 * Scratch memory (per thread)
 **/
#define PP_SCRATCH_MEMORY
#define PP_SCRATCH_MEMORY_HBW
struct scratchMem {
  // scratch memory for the surface integration
  real_base tResSurf[2][N_QUANTITIES][N_FACE_MODES][N_CRUNS] __attribute__ ((aligned (ALIGNMENT.FACE_MODES.PRIVATE)));
  // temporary results
  real_base tRes[2][N_ELEMENT_MODES][N_CRUNS] __attribute__ ((aligned (ALIGNMENT.ELEMENT_MODES.PRIVATE)));
  // sub-cells at the DG faces
  real_base scFa[4][N_QUANTITIES][CE_N_SUB_FACES(T_SDISC.ELEMENT, ORDER)][N_CRUNS] __attribute__ ((aligned (ALIGNMENT.CRUNS)));
  // sub-grid
  real_base sg[N_QUANTITIES][CE_N_SUB_CELLS(T_SDISC.ELEMENT,ORDER)][N_CRUNS] __attribute__ ((aligned (ALIGNMENT.CRUNS)));
  // sub-grid with ghost cells
  real_base sgGh[N_QUANTITIES][   CE_N_SUB_CELLS(T_SDISC.ELEMENT, ORDER)
                                + C_ENT[T_SDISC.ELEMENT].N_FACES *  CE_N_SUB_FACES(T_SDISC.ELEMENT, ORDER)
                              ][N_CRUNS] __attribute__ ((aligned (ALIGNMENT.CRUNS)));
};
typedef scratchMem t_scratchMem;


typedef enum {
  OUTFLOW = 105
} t_bndConds;

/*
 * Advection equation.
 * Varying coefficients (space-depenent scalar a, b, c) are allowed:
 *
 * q_t + f(q)_x + g(q)_y + h(q)_z = 0
 *
 * with
 *   f(q) = a(x)*q
 *   g(q) = b(x)*q
 *   h(q) = c(x)*q
 */
typedef struct {
  // advection speed in x-direction
  real_base a;
#if defined PP_T_ELEMENTS_QUAD4R || defined PP_T_ELEMENTS_TRIA3 || defined PP_T_ELEMENTS_HEX8R || defined PP_T_ELEMENTS_TET4
  // advection speed in y-direction
  real_base b;
#endif
#if defined PP_T_ELEMENTS_HEX8R || defined PP_T_ELEMENTS_TET4
  // advection speed in z-direction
  real_base c;
#endif
} t_bgPars;

// DOFs
#define PP_N_ELEMENT_MODE_PRIVATE_1 1
typedef real_base t_elementModePrivate1;

// flex data holding a time buffer or a time buffer and element-DOFs per entry
#define PP_N_GLOBAL_SHARED_3 3
typedef real_base (**t_globalShared3)[N_ELEMENT_MODES][N_CRUNS];

// parameters a, b, c
#define PP_N_ELEMENT_SHARED_1 1
typedef t_bgPars t_elementShared1;

// star matrices in xi-, eta- and zeta-direction
#define PP_N_ELEMENT_SHARED_2 N_DIM
typedef real_base t_elementShared2;

// flux solvers
#define PP_N_ELEMENT_SHARED_3 (C_ENT[T_SDISC.ELEMENT].N_FACES*2)
typedef real_base t_elementShared3;
