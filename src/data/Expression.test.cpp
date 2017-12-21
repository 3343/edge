/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2017, Regents of the University of California
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
 * Unit tests for EDGE's expressions.
 **/
#include <catch.hpp>
#define private public
#include "Expression.hpp"
#undef private

TEST_CASE( "Expression: Scalar.", "[expression][scalar]" ) {
  edge::data::Expression< float > l_expr1;

  std::string l_sym1  = "x";
  std::string l_sym2  = "y";
  std::string l_sym3  = "z";

  float l_in1[2];
  float l_out1;

  std::string l_exprStr1 = "z := x - (3 * y)";


  l_expr1.bind( l_sym1, l_in1[ 0] );
  l_expr1.bind( l_sym2, l_in1[ 1] );
  l_expr1.bind( l_sym3, l_out1    );
  l_expr1.compile( l_exprStr1 );

  // check compiled expression
  l_in1[0] = 2;
  l_in1[1] = 3;
  l_expr1.eval();
  REQUIRE( l_out1 == Approx( -7.0 ) );


  l_in1[0] = -2;
  l_in1[1] = 8;
  l_expr1.eval();
  REQUIRE( l_out1 == Approx( -26 ) );
}

TEST_CASE( "Expression: Array.", "[expression][array]" ) {
  std::string l_sym1  = "x";
  std::string l_sym2  = "y";
  std::string l_sym3  = "q";

  double l_in1[2];
  double l_out1[6];

  std::string l_exprStr1 = "q[0] := x + y;\
                            q[1] := x * y;\
                            q[2] := x - y;\
                            q[3] := x;\
                            q[4] := y;\
                            q[5] := x / y;\
                           ";

  edge::data::Expression< double > l_expr1;
  l_expr1.bind( l_sym1, l_in1[0] );
  l_expr1.bind( l_sym2, l_in1[1] );
  l_expr1.bind( l_sym3, l_out1, 6 );
  l_expr1.compile( l_exprStr1 );

  // check compiled expressions
  l_in1[0] =  2.5;
  l_in1[1] = -3.0;
  l_expr1.eval();

  REQUIRE( l_out1[0] == Approx( -0.5 ) );
  REQUIRE( l_out1[1] == Approx( -7.5 ) );
  REQUIRE( l_out1[2] == Approx(  5.5 ) );
  REQUIRE( l_out1[3] == Approx(  2.5 ) );
  REQUIRE( l_out1[4] == Approx( -3.0 ) );
  REQUIRE( l_out1[5] == Approx( -2.5 / 3.0 ) );
}

TEST_CASE( "Expression: Bind coordinates.", "[expression][bindCrds]" ) {
  std::string l_sym1  = "someVar";

  double l_in1[3];
  double l_out1;

  edge::data::Expression< double > l_expr1;
  l_expr1.bindCrds( l_in1, 3 );
  l_expr1.bind( l_sym1, l_out1 );

  std::string l_exprStr1 = "someVar := cos(x) * y + z;";
  l_expr1.compile( l_exprStr1 );

  l_in1[0] =  0.0;
  l_in1[1] =  2.0;
  l_in1[2] = -4.0;
  l_expr1.eval();

  REQUIRE( l_out1 == Approx(-2.0) );
}
