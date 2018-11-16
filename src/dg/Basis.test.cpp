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
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONsTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 * Unit tests of the DG basis.
 **/
#include <catch.hpp>
#define private public
#include "Basis.h"
#undef private

TEST_CASE( "Tests the precomputation of the basis values and quad points.", "[basis][evalQpRefEl]" ) {
#ifndef PP_T_ELEMENTS_TET4
  // return if not compiled for tets
  return;
#endif

  edge::dg::Basis l_basis( TET4, PP_ORDER );

  for( unsigned short l_pq = 0; l_pq < 2+1; l_pq++ ) {
    REQUIRE( l_basis.m_qpEval[l_pq].xi1.size()     == (l_pq+1) * (l_pq+1) * (l_pq+1) );
    REQUIRE( l_basis.m_qpEval[l_pq].xi2.size()     == (l_pq+1) * (l_pq+1) * (l_pq+1) );
    REQUIRE( l_basis.m_qpEval[l_pq].xi3.size()     == (l_pq+1) * (l_pq+1) * (l_pq+1) );
    REQUIRE( l_basis.m_qpEval[l_pq].weights.size() == (l_pq+1) * (l_pq+1) * (l_pq+1) );
    REQUIRE( l_basis.m_qpEval[l_pq].basis.size()   == (l_pq+1) * (l_pq+1) * (l_pq+1) );
    REQUIRE( l_basis.m_qpEval[l_pq].basis.size()   == (l_pq+1) * (l_pq+1) * (l_pq+1) );

    for( int l_lo = 0; l_lo < (l_pq+1) * (l_pq+1) * (l_pq+1); l_lo++ ) {
      for( unsigned short l_pb = 0; l_pb < 2+1; l_pb++ ) {
        unsigned short l_nModes = std::numeric_limits< unsigned short >::max();
        if(      l_pb == 0 ) l_nModes = 1;
        else if( l_pb == 1 ) l_nModes = 4;
        else if( l_pb == 2 ) l_nModes = 10;
        else if( l_pb == 3 ) l_nModes = 20;
        else if( l_pb == 4 ) l_nModes = 35;

        REQUIRE( l_basis.m_qpEval[l_pq].basis[l_lo][l_pb].val[0]         == Approx(1) );
        REQUIRE( l_basis.m_qpEval[l_pq].basis[l_lo][l_pb].valD[0].size() == l_nModes  );
        REQUIRE( l_basis.m_qpEval[l_pq].basis[l_lo][l_pb].valD[1].size() == l_nModes  );
        REQUIRE( l_basis.m_qpEval[l_pq].basis[l_lo][l_pb].valD[2].size() == l_nModes  );
      }
    }
  }
}
