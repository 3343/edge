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
 * Compile time constants for the elatic wave equations.
 **/

#include "sc/ibnd/InternalBoundary.type"

// elastics perform three steps per time step
// 0) time prediction + local cont
// 1) neighboring cont
// 2) sources
const unsigned short N_STEPS_PER_UPDATE=3;
#ifdef PP_T_EQUATIONS_ELASTIC_RUPTURE
const unsigned short N_ENTRIES_CONTROL_FLOW=10;
#else
// TODO: Add limiter
const unsigned short N_ENTRIES_CONTROL_FLOW=8;
#endif

const real_base C_SCALE_ENTROPY_FIX_HARTEN = 0.05;

#if PP_ORDER > 9
#error not supported
#endif
// CFL-numbers based on (experimental) 1D-stability of advection equation
const double SCALE_CFL_ORDER[9] = { 0.45, 0.9, 0.8, 0.7, 0.6, 0.5, 0.45, 0.4, 0.4};

const double SCALE_CFL=0.5*SCALE_CFL_ORDER[PP_ORDER-1];

/*
 * Bitsmasks for elastics.
 **/
static_assert( sizeof(int_spType) >= 4, "size of entity type too small" );
typedef enum: int_spType {
  FREE_SURFACE = 101,  // TODO: change to bitmask
  OUTFLOW      = 105,  // TODO: change to bitmask
  RUPTURE      = 201,  // TODO: change to bitmask
  SOURCE       = 65536 // 0b0000000000000000000000000000000000000000000000010000000000000000
} t_spTypeElastic;

#if PP_N_DIM == 2
// two normal stress components (sigma_xx, sigma_yy), one shear stress component (sigma_xy), two particle velocites (u, v)
const unsigned short N_QUANTITIES = 5;
#elif PP_N_DIM == 3
// three normal stress components (signa_xx, sigma_yy, sigmag_zz), three shear stress (sigma_xy, sigma_yz, sigma_xz), three particle velocities (u,v,w)
const unsigned short N_QUANTITIES = 9;
#else
#error elastic constants for other dimensions than 2 and 3 undefined.
#endif

typedef struct {
  // density rho
  real_base rho;
  // lame parameter lambda
  real_base lam;
  // lame parameter mu
  real_base mu;
} t_bgPars;

typedef struct {
  real_base solver[N_QUANTITIES][N_QUANTITIES];
} t_fluxSolver;


/**
 * Scratch memory (per thread)
 **/
#define PP_SCRATCH_MEMORY
#define PP_SCRATCH_MEMORY_HBW
struct scratchMem {
  // temporary results
  real_base tRes[N_QUANTITIES][N_ELEMENT_MODES][N_CRUNS] __attribute__ ((aligned (ALIGNMENT.ELEMENT_MODES.PRIVATE)));
  // derivative buffer
  real_base dBuf[ORDER][N_QUANTITIES][N_ELEMENT_MODES][N_CRUNS] __attribute__ ((aligned (ALIGNMENT.ELEMENT_MODES.PRIVATE)));
};
typedef scratchMem t_scratchMem;

#include "src/impl/elastic/solvers/Kinematics.type"
#include "src/impl/elastic/solvers/FrictionLaws.type"
#include "src/impl/elastic/solvers/InternalBoundary.type"