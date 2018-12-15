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
 * Sets processor options.
 **/
#include "Cpu.h"
#include "xmmintrin.h"
#include "pmmintrin.h"
#include "io/logging.h"

bool edge::setups::Cpu::getFlushToZero() {
#ifdef __SSE__
  unsigned int l_ftz = _MM_GET_FLUSH_ZERO_MODE();

  if( l_ftz == _MM_FLUSH_ZERO_ON ) {
    return true;
  }
  else if( l_ftz != _MM_FLUSH_ZERO_OFF ) {
    EDGE_LOG_FATAL;
  }
#endif

  return false;
}

void edge::setups::Cpu::setFlushToZero( bool i_on ) {
#ifdef PP_USE_OMP
#pragma omp parallel
#endif
{

#ifdef __SSE__
  if( i_on )
    _MM_SET_FLUSH_ZERO_MODE( _MM_FLUSH_ZERO_ON );
  else
    _MM_SET_FLUSH_ZERO_MODE( _MM_FLUSH_ZERO_OFF );
#endif

}
}

bool edge::setups::Cpu::getDenormalsAreZero() {
#ifdef __SSE__
  unsigned int l_ftz = _MM_GET_DENORMALS_ZERO_MODE();

  if( l_ftz == _MM_DENORMALS_ZERO_ON ) {
    return true;
  }
  else if( l_ftz != _MM_DENORMALS_ZERO_OFF ) {
    EDGE_LOG_FATAL;
  }
#endif

  return false;
}

void edge::setups::Cpu::setDenormalsAreZero( bool i_on ) {
#ifdef PP_USE_OMP
#pragma omp parallel
#endif
{

#ifdef __SSE__
  if( i_on )
    _MM_SET_DENORMALS_ZERO_MODE( _MM_DENORMALS_ZERO_ON );
  else
    _MM_SET_DENORMALS_ZERO_MODE( _MM_DENORMALS_ZERO_OFF );
#endif

}
}