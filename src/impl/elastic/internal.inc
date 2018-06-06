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
 * Internal data of seismic simulations.
 **/

/*
 * First global data encodes the setup
 */
#define PP_N_GLOBAL_PRIVATE_1 1
struct setupChars {
  enum InitialDofs { zero, planeWaves } initialDofs;
};
typedef setupChars t_globalPrivate1;

/*
 * Source terms
 */
// maximum number of independent sources if given by the number of fused runs
#define PP_N_GLOBAL_SHARED_3 PP_N_CRUNS
typedef edge::elastic::solvers::t_Kinematics< N_DIM,
                                              N_ELEMENT_MODES,
                                              1,
                                              real_base,
                                              int_el > t_globalShared3;

// positions of first point source associated to the sparse source elements
#define PP_N_GLOBAL_SHARED_4 1
typedef int_el (*t_globalShared4)[N_CRUNS];

/*
 * Rupture physics
 */
#define PP_N_GLOBAL_SHARED_5 1
typedef edge::elastic::solvers::t_LinSlipWeak< real_base,
                                               T_SDISC.ELEMENT,
                                               ORDER,
                                               N_CRUNS> t_globalShared5;
// DOFs of the rupture elements
#define PP_N_ELEMENT_SPARSE_SHARED_3 1
typedef real_base t_elementSparseShared3[N_QUANTITIES][N_ELEMENT_MODES][N_CRUNS];

/*
 * First shared element data are the background parameters.
 */
#define PP_N_ELEMENT_SHARED_1 1
typedef t_bgPars t_elementShared1;

/*
 * Second shared element data are the flux solvers for the element's own contribution.
 */
#define PP_N_ELEMENT_SHARED_2 C_ENT[T_SDISC.ELEMENT].N_FACES
typedef t_fluxSolver t_elementShared2;

/*
 * Third shared element data are the flux solvers for the neighboring elements' contribution.
 */
#define PP_N_ELEMENT_SHARED_3 C_ENT[T_SDISC.ELEMENT].N_FACES
typedef t_fluxSolver t_elementShared3;

// setup number of entries in the star matrices
#if defined PP_T_KERNELS_XSMM
const unsigned short N_MAT_STAR = (N_DIM==2) ? 10 : 24;
#endif

/*
 * Fourth shared element data are the star matrices.
 */
typedef struct {
#if defined PP_T_KERNELS_XSMM
  real_base mat[N_MAT_STAR];
#else
  real_base mat[N_QUANTITIES][N_QUANTITIES];
#endif
} t_matStar;

#define PP_N_ELEMENT_SHARED_4 N_DIM
typedef t_matStar t_elementShared4;

/*
 * First private element mode data are the first set of dofs.
 */
#define PP_N_ELEMENT_MODE_PRIVATE_1 N_QUANTITIES
typedef real_base t_elementModePrivate1;

// flex data holding a time buffer or a time buffer and element-DOFs per entry
#define PP_N_GLOBAL_SHARED_6 2
typedef real_base (**t_globalShared6)[N_ELEMENT_MODES][N_CRUNS];

#define PP_N_GLOBAL_SHARED_7 1
typedef edge::sc::ibnd::t_InternalBoundary< int_el,
                                            real_base,
                                            t_spTypeElastic,
                                            T_SDISC.ELEMENT,
                                            N_QUANTITIES > t_globalShared7;

//#include "impl/elastic/sc/Llf.hpp"
#include "impl/elastic/sc/UpWind.hpp"
#define PP_N_GLOBAL_SHARED_8 1
typedef edge::elastic::sc::UpWind< real_base,
                                   T_SDISC.ELEMENT,
                                   ORDER,
                                   N_CRUNS > t_globalShared8;