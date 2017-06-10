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
typedef real_base t_elementShared1;
typedef real_base t_faceModePrivate1;

#if defined PP_T_ELEMENTS_LINE

const unsigned short N_QUANTITIES=2;

// DOFs are h and hu
#define PP_N_ELEMENT_MODE_PRIVATE_1 N_QUANTITIES

// FV for now
#define PP_N_ELEMENT_MODES 1

// bathymetry
#define PP_N_ELEMENT_MODE_SHARED_1 1

// dx
#define PP_N_ELEMENT_SHARED_1 1

// net-update are stored for 1 mode per face
const unsigned short N_FACE_MODES=1;

// left- and right-going net-updates for h and hu
#define PP_N_FACE_MODE_PRIVATE_1 4


#elif defined PP_T_ELEMENTS_QUAD4R

const unsigned short N_QUANTITIES=3;

#define PP_N_ELEMENT_MODE_PRIVATE_1 N_QUANTITIES
#define PP_N_ELEMENT_MODES 1
#define PP_N_ELEMENT_MODE_SHARED_1 1
#define PP_N_ELEMENT_SHARED_1 1
const unsigned short N_FACE_MODES=1;
#define PP_N_FACE_MODE_PRIVATE_1 4

#else
#error shallow water constants for element type undefined.
#endif
