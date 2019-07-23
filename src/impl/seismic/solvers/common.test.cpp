/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2016-2018, Regents of the University of California
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
 * Unit tests for common solvers of the elastic wave equations.
 **/

#include <catch.hpp>
#define private public
#include "common.hpp"
#undef private

TEST_CASE( "Setup of the two-dimensional, elastic solvers.", "[solversCommon][elastic2d]" ) {
  double l_fL[5][5];
  double l_fR[5][5];

  edge::seismic::solvers::common::setupSolver2d( 1.0, 1.0,
                                                 2.0, 2.0,
                                                 1.0, 1.0,
                                                 -0.2, std::sqrt(1-0.04),
                                                 l_fL,
                                                 l_fR );

  // check some values
  REQUIRE( l_fL[0][0] == Approx( 0.0592)              );
  REQUIRE( l_fL[1][0] == Approx( 0.0008)              );
  REQUIRE( l_fL[2][0] == Approx(-0.094060406122874)   );
  REQUIRE( l_fL[3][0] == Approx( 0.1)                 );
  REQUIRE( l_fL[4][0] == Approx( 0.0)                 );

  REQUIRE( l_fL[0][1] == Approx( 0.4608)              );
  REQUIRE( l_fL[1][1] == Approx( 0.9792)              );
  REQUIRE( l_fL[2][1] == Approx(-0.00391918358845308) );
  REQUIRE( l_fL[3][1] == Approx( 0.0)                 );
  REQUIRE( l_fL[4][1] == Approx(-0.489897948556636)   );

  REQUIRE( l_fR[0][2] == Approx( 0.384079991668402)   );
  REQUIRE( l_fR[1][2] == Approx( 0.20379754659956)    );
  REQUIRE( l_fR[2][2] == Approx(-0.4616)              );
  REQUIRE( l_fR[3][2] == Approx(-0.489897948556636)   );
  REQUIRE( l_fR[4][2] == Approx( 0.1)                 );

  REQUIRE( l_fR[0][3] == Approx( 0.4)                 );
  REQUIRE( l_fR[1][3] == Approx( 0.2)                 );
  REQUIRE( l_fR[2][3] == Approx(-0.489897948556636)   );
  REQUIRE( l_fR[3][3] == Approx(-0.52)                );
  REQUIRE( l_fR[4][3] == Approx( 0.0979795897113271)  );

  REQUIRE( l_fL[0][4] == Approx(-0.979795897113271)   );
  REQUIRE( l_fL[1][4] == Approx(-1.95959179422654)    );
  REQUIRE( l_fL[2][4] == Approx( 0.1)                 );
  REQUIRE( l_fL[3][4] == Approx(-0.0979795897113271)  );
  REQUIRE( l_fL[4][4] == Approx( 0.98)                );

  edge::seismic::solvers::common::setupSolver2d( 1.0, 1.0,
                                                 2.0, 2.0,
                                                 1.0, 1.0,
                                                 0.5, -std::sqrt(0.75),
                                                 l_fL,
                                                 l_fR );

  // check some values
  REQUIRE( l_fL[0][0] == Approx( 0.34375)            );
  REQUIRE( l_fL[1][0] == Approx( 0.03125)            );
  REQUIRE( l_fL[2][0] == Approx(-0.162379763209582)  );
  REQUIRE( l_fL[3][0] == Approx(-0.25)               );
  REQUIRE( l_fL[4][0] == Approx( 0.0)                );

  REQUIRE( l_fL[0][1] == Approx( 0.28125)            );
  REQUIRE( l_fL[1][1] == Approx( 0.84375)            );
  REQUIRE( l_fL[2][1] == Approx(-0.0541265877365274) );
  REQUIRE( l_fL[3][1] == Approx( 0.0)                );
  REQUIRE( l_fL[4][1] == Approx( 0.433012701892219)  );

  REQUIRE( l_fR[0][2] == Approx( 0.757772228311384)  );
  REQUIRE( l_fR[1][2] == Approx( 0.541265877365274)  );
  REQUIRE( l_fR[2][2] == Approx(-0.3125)             );
  REQUIRE( l_fR[3][2] == Approx( 0.433012701892219)  );
  REQUIRE( l_fR[4][2] == Approx(-0.25)               );

  REQUIRE( l_fL[0][3] == Approx(-1.0)                );
  REQUIRE( l_fL[1][3] == Approx(-0.5)                );
  REQUIRE( l_fL[2][3] == Approx( 0.433012701892219)  );
  REQUIRE( l_fL[3][3] == Approx( 0.625)              );
  REQUIRE( l_fL[4][3] == Approx(-0.21650635094611)   );

  REQUIRE( l_fL[0][4] == Approx( 0.866025403784439)  );
  REQUIRE( l_fL[1][4] == Approx( 1.73205080756888)   );
  REQUIRE( l_fL[2][4] == Approx(-0.25)               );
  REQUIRE( l_fL[3][4] == Approx(-0.21650635094611)   );
  REQUIRE( l_fL[4][4] == Approx( 0.875)              );

  edge::seismic::solvers::common::setupSolver2d( 1.0, 1.0,
                                                 2.0, 2.0,
                                                 1.0, 1.0,
                                                 1/std::sqrt(2), 1/std::sqrt(2),
                                                 l_fL,
                                                 l_fR );

  // check some values
  REQUIRE( l_fL[0][0] == Approx( 0.625)               );
  REQUIRE( l_fL[1][0] == Approx( 0.125)               );
  REQUIRE( l_fL[2][0] == Approx( 0.125)               );
  REQUIRE( l_fL[3][0] == Approx(-0.25 * std::sqrt(2)) );
  REQUIRE( l_fL[4][0] == Approx( 0.0)                 );

  REQUIRE( l_fL[0][1] == Approx( 0.125)               );
  REQUIRE( l_fL[1][1] == Approx( 0.625)               );
  REQUIRE( l_fL[2][1] == Approx( 0.125)               );
  REQUIRE( l_fL[3][1] == Approx( 0.0)                 );
  REQUIRE( l_fL[4][1] == Approx(-0.25 * std::sqrt(2)) );

  REQUIRE( l_fR[0][2] == Approx(-0.75)                );
  REQUIRE( l_fR[1][2] == Approx(-0.75)                );
  REQUIRE( l_fR[2][2] == Approx(-0.25)                );
  REQUIRE( l_fR[3][2] == Approx(-0.25 * std::sqrt(2)) );
  REQUIRE( l_fR[4][2] == Approx(-0.25 * std::sqrt(2)) );

  REQUIRE( l_fL[0][3] == Approx(-std::sqrt(2))        );
  REQUIRE( l_fL[1][3] == Approx(-0.5  * std::sqrt(2)) );
  REQUIRE( l_fL[2][3] == Approx(-0.25 * std::sqrt(2)) );
  REQUIRE( l_fL[3][3] == Approx( 0.75)                );
  REQUIRE( l_fL[4][3] == Approx( 0.25)                );

  REQUIRE( l_fL[0][4] == Approx(-0.5 *  std::sqrt(2)) );
  REQUIRE( l_fL[1][4] == Approx(-std::sqrt(2))        );
  REQUIRE( l_fL[2][4] == Approx(-0.25 * std::sqrt(2)) );
  REQUIRE( l_fL[3][4] == Approx( 0.25)                );
  REQUIRE( l_fL[4][4] == Approx( 0.75)                );
}
