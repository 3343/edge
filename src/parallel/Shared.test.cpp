/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2018, Regents of the University of California
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
 * Unit tests of the shared memory implementation.
 **/
#include <catch.hpp>

#ifdef PP_USE_OMP
#include <omp.h>
#endif

#define private public
#include "Shared.h"
#undef private

TEST_CASE( "NUMA-aware initialization without separate workers.", "[numaInitNoSep]" ) {

  // our test arrays, which are have to be set to zero
  float l_arr1[73*31];
  for( unsigned int l_en = 0; l_en < 73*31; l_en++ ) l_arr1[l_en] = float(1);
  float l_arr2[3] = {1, 1, 1};

  // open OMP-region if required
#ifdef PP_USE_OMP
#pragma omp parallel
#endif
  {

#ifdef PP_USE_OMP
    edge::parallel::Shared l_shared;
    l_shared.init( false );
#else
    edge::parallel::g_thread   = 0;
    edge::parallel::g_nThreads = 1;
    edge::parallel::
#endif

    l_shared.numaInit( 73*31, l_arr1 );
    l_shared.numaInit( 3,     l_arr2 );
  }

  // check the result
  for( unsigned int l_en = 0; l_en < 73*31; l_en++ )  REQUIRE( l_arr1[l_en] == float(0) );
  for( unsigned int l_en = 0; l_en <     3; l_en++ )  REQUIRE( l_arr2[l_en] == float(0) );
}

TEST_CASE( "NUMA-aware initialization with separate workers.", "[numaInitSep]" ) {

  // our test arrays, which are have to be set to zero
  float l_arr1[73*31];
  for( unsigned int l_en = 0; l_en < 73*31; l_en++ ) l_arr1[l_en] = float(1);
  float l_arr2[3] = {1, 1, 1};

  // open OMP-region if required
#ifdef PP_USE_OMP
#pragma omp parallel
#endif
  {

#ifdef PP_USE_OMP
    edge::parallel::Shared l_shared;
    l_shared.init( true );
#else
    edge::parallel::g_thread   = 0;
    edge::parallel::g_nThreads = 1;
    edge::parallel::
#endif

    l_shared.numaInit( 73*31, l_arr1 );
    l_shared.numaInit( 3,     l_arr2 );
  }

  // check the result
  for( unsigned int l_en = 0; l_en < 73*31; l_en++ )  REQUIRE( l_arr1[l_en] == float(0) );
  for( unsigned int l_en = 0; l_en <     3; l_en++ )  REQUIRE( l_arr2[l_en] == float(0) );
}