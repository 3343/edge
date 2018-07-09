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
 * Compile time constants for the shallow water equations.
 **/

// shallow water equations perform a single update
const unsigned short N_STEPS_PER_UPDATE=2;
const unsigned short N_ENTRIES_CONTROL_FLOW=2;

#if PP_ORDER > 1
#error only fv for swe.
#endif

typedef real_base t_elementModePrivate1;
typedef real_base t_elementModeShared1;
typedef real_base t_faceModePrivate1;


// bathymetry
#define PP_N_ELEMENT_MODE_SHARED_1 1
const unsigned short N_QUANTITIES=N_DIM+1;

// DOFs are h and hu
#define PP_N_ELEMENT_MODE_PRIVATE_1 N_QUANTITIES

// left- and right-going net-updates for h and hu
#define PP_N_FACE_MODE_PRIVATE_1 4

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