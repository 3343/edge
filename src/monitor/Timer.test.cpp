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
 * Tests the timer.
 **/
#include <catch.hpp>
#include "io/logging.h"
#define private public
#include "Timer.hpp"
#undef private

TEST_CASE( "Monitor: Timer.", "[timer][monitor]" ) {
  // create a timer
  edge::monitor::Timer l_timer;

  // check query of wall-clock time
  REQUIRE( l_timer.getWtime() > 0 );

  // start the timer
  l_timer.start();

  // keep the CPU busy
  unsigned int l_state = 0;
  for( unsigned int l_it = 0; l_it < 10000; l_it++ ) l_state+=1.0;
  REQUIRE( l_state == 10000 );

  // stop time and check that time elapsed
  l_timer.end();
  REQUIRE( l_timer.elapsed() > 0 );

  // reset time and check for zero seconds elapsed
  l_timer.reset();
  REQUIRE( l_timer.elapsed() == 0.0 );

  // start the timer
  l_timer.start();

  // do some work
  l_state = 0;
  for( unsigned int l_it = 0; l_it < 10000; l_it++ ) l_state+=1.0;
  REQUIRE( l_state == 10000 );

  // stop timer and get time
  l_timer.end();
  double l_elapsed = l_timer.elapsed();

  l_timer.start();

  // do some work
  l_state = 0;
  for( unsigned int l_it = 0; l_it < 10000; l_it++ ) l_state+=1.0;
  REQUIRE( l_state == 10000 );

  // stop timer and get time
  l_timer.end();

  // check that the elapsed time increased
  REQUIRE( l_timer.elapsed() > l_elapsed );
}