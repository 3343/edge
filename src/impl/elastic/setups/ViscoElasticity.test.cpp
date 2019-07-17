/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2019, Alexander Breuer
 * Copyright (c) 2019, Regents of the University of California
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
 * Unit tests for viscoelastic attenuation.
 **/

#include <catch.hpp>
#define private public
#include "ViscoElasticity.h"
#undef private

/*
 * The solutions of the test cases in this section can be
 * reproduced with the following python code:
 *
 *
    import sympy
    sympy.init_printing(use_unicode=True)

    l_fc = 2
    l_fr = 10

    l_omin = ( 2 * sympy.pi * l_fc )  / sympy.sqrt(l_fr)

    # compute relaxation frequencies
    l_freqs = []
    l_freqs.append( l_omin )
    for l_k in [1,2,3,4]:
      l_freqs.append( sympy.exp( sympy.ln(l_omin) + ( l_k*sympy.ln(l_fr) ) / (2* (3 - 1) ) ) )

    for l_k in l_freqs:
      print( sympy.N(l_k) )
  #
  # 3.97383530631844
  # 7.06658950411180
  # 12.5663706143592
  # 22.3465181224331
  # 39.7383530631844
  #
    # assemble matrices for least squares solve
    l_mats = []

    for l_qf in [2000*0.1, 2000*0.1*1.5]:
      l_mats.append( sympy.zeros(5,3) )

      for l_l in range(3):
        for l_k in range(5):
          l_mats[-1][l_k,l_l] = l_freqs[2*l_l] * l_freqs[l_k] + l_freqs[2*l_l]**2 * (1/l_qf)
          l_mats[-1][l_k,l_l] = l_mats[-1][l_k,l_l] / ( l_freqs[l_k]**2 + l_freqs[2*l_l]**2 )

    # assemble right hand sides
    l_rhs = []
    for l_qf in [2000*0.1, 2000*0.1*1.5]:
      l_rhs.append( sympy.ones(5,1) * (1/l_qf) )
      
    # perform the least square solves
    l_fit = []
    for l_f in range(2):
      l_fit.append( l_mats[l_f].solve_least_squares( l_rhs[l_f] ) )

    # compute solution
    l_fit = [ sympy.N( l_fit[0] ), sympy.N( l_fit[1] ) ]
    print( l_fit[0] )
    print( l_fit[1] )
  #
  # Matrix([[0.00788212871452061], [0.000633529971773939], [0.00797631601780113]])
  # Matrix([[0.00527973491702664], [0.000423653588369432], [0.00532171233332960]])
  #
    l_lam = 20.8E9
    l_mu = 10.4E9

    import math

    l_yp = l_fit[0]
    l_ys = l_fit[1]

    l_omegaRef = 2 * math.pi * l_fc

    l_theta1kappa = 1
    l_theta1mu = 1
    l_theta2kappa = 0
    l_theta2mu = 0
    for l_l in range(3):
      l_div = 1
      l_div /= 1 + ( l_omegaRef / sympy.N(l_freqs[l_l*2]) )**2

      l_theta1kappa -= l_yp[l_l] * l_div
      l_theta2kappa += l_yp[l_l] * (l_omegaRef / sympy.N(l_freqs[l_l*2])) * l_div

      l_theta1mu -= l_ys[l_l] * l_div
      l_theta2mu += l_ys[l_l] * (l_omegaRef / sympy.N(l_freqs[l_l*2])) * l_div
        
    l_rkappa = math.sqrt( l_theta1kappa**2 + l_theta2kappa**2 )
    l_rmu = math.sqrt( l_theta1mu**2 + l_theta2mu**2 )

    l_muUnrelaxed = l_mu * (l_rmu + l_theta1mu) / (2*l_rmu**2)

    l_kappa = (l_lam + 2*l_mu) * (l_rkappa + l_theta1kappa) / (2*l_rkappa**2)
    l_lamUnrelaxed = l_kappa - 2*l_muUnrelaxed

    print( 'elastic lam, mu' )
    print( l_lamUnrelaxed, l_muUnrelaxed )
  # elastic lam, mu
  # 21031265249.9435 10457744636.3789

    # compute lame coefficients
    l_ylam = ( 1 + (2*l_muUnrelaxed) / l_lamUnrelaxed ) * l_yp - (2*l_muUnrelaxed / l_lamUnrelaxed) * l_ys
    l_ymu = l_ys
    print( 'lame coeffs' )
    print( l_ylam )
    print( l_ymu )
  #
  # lame coeffs
  # Matrix([[0.0104701964749766], [0.000842250997044352], [0.0106163062526320]])
  # Matrix([[0.00527973491702664], [0.000423653588369427], [0.00532171233332960]])
  #
    # compute source matrix entries
    print( 'source matrix entries' )
    print( -l_lamUnrelaxed * l_ylam - 2 * l_muUnrelaxed * l_ymu )
    print( -l_lamUnrelaxed * l_ylam )
    print( - 2 * l_muUnrelaxed * l_ymu )
  #
  # source matrix entries
  # Matrix([[-330629718.304332], [-26574526.2087753], [-334580569.994195]])
  # Matrix([[-220201479.284256], [-17713604.1258692], [-223274352.773737]])
  # Matrix([[-110428239.020076], [-8860922.08290611], [-111306217.220458]])
  */

TEST_CASE( "Derivation of relaxation frequencies.", "[ViscoElasticity][frequencies]" ) {
  double l_freqs1[5];

  edge::seismic::setups::ViscoElasticity::frequencies( 3,
                                                       2,
                                                       10,
                                                       l_freqs1 );

  REQUIRE( l_freqs1[0] == Approx(3.97383530631844) );
  REQUIRE( l_freqs1[1] == Approx(7.06658950411180) );
  REQUIRE( l_freqs1[2] == Approx(12.5663706143592) );
  REQUIRE( l_freqs1[3] == Approx(22.3465181224331) );
  REQUIRE( l_freqs1[4] == Approx(39.7383530631844) );
}

TEST_CASE( "Derivation of anelastic coefficients, p/s-version", "[ViscoElasticity][anelasticCoeffPs]" ) {
  double l_freqs1[5];
  edge::seismic::setups::ViscoElasticity::frequencies( 3,
                                                       2,
                                                       10,
                                                       l_freqs1 );

  double l_vs1 = 2000;
  double l_qs1 = l_vs1 * 0.1;
  double l_qp1 = l_qs1 * 1.5;

  double l_fit1[3];

  edge::seismic::setups::ViscoElasticity::anelasticCoeffPs( 3,
                                                            l_qs1,
                                                            l_freqs1,
                                                            l_fit1 );

  REQUIRE( l_fit1[0] == Approx(0.00788212871452061)  );
  REQUIRE( l_fit1[1] == Approx(0.000633529971773939) );
  REQUIRE( l_fit1[2] == Approx(0.00797631601780113)  );

  edge::seismic::setups::ViscoElasticity::anelasticCoeffPs( 3,
                                                            l_qp1,
                                                            l_freqs1,
                                                            l_fit1 );

  REQUIRE( l_fit1[0] == Approx(0.00527973491702664)  );
  REQUIRE( l_fit1[1] == Approx(0.000423653588369432) );
  REQUIRE( l_fit1[2] == Approx(0.00532171233332960)  );
}

TEST_CASE( "Derivation of elastic Lame parameters.", "[ViscoElasticity][elasticLame]" ) {
  double l_lamElastic = std::numeric_limits< double >::max();
  double l_muElastic = std::numeric_limits< double >::max();

  double l_coeffsKappa[3] = { 0.00788212871452061,
                              0.000633529971773939,
                              0.00797631601780113 };

  double l_coeffsMu[3] = { 0.00527973491702664,
                           0.000423653588369432,
                           0.00532171233332960 };

  double l_freqRef = 2 * M_PI * 2;

  double l_freqs[5] = { 3.97383530631844,
                        7.06658950411180,
                        12.5663706143592,
                        22.3465181224331,
                        39.7383530631844 };

  edge::seismic::setups::ViscoElasticity::elasticLame( 3,
                                                       l_freqRef,
                                                       l_freqs,
                                                       l_coeffsKappa,
                                                       l_coeffsMu,
                                                       20.8E9,
                                                       10.4E9,
                                                       l_lamElastic,
                                                       l_muElastic );
  REQUIRE( l_lamElastic == Approx(21031265249.9435) );
  REQUIRE( l_muElastic == Approx(10457744636.3789) );
}

TEST_CASE( "Derivation of anelastic coefficients, lame-version", "[ViscoElasticity][anelasticCoeffLame]" ) {
  double l_coeffs1[2*3];

  double l_lameE[2] = {-1, -1};

  edge::seismic::setups::ViscoElasticity::anelasticCoeffLame( 3,
                                                              2,
                                                              10,
                                                              2000*0.1,
                                                              2000*0.1*1.5,
                                                              20.8E9,
                                                              10.4E9,
                                                              l_lameE[0],
                                                              l_lameE[1],
                                                              l_coeffs1 );
  REQUIRE( l_lameE[0] == Approx(21031265249.9435) );
  REQUIRE( l_lameE[1] == Approx(10457744636.3789) );

  REQUIRE( l_coeffs1[0+0] == Approx(0.0104701964749766)   );
  REQUIRE( l_coeffs1[0+1] == Approx(0.000842250997044352) );
  REQUIRE( l_coeffs1[0+2] == Approx(0.0106163062526320)   );

  REQUIRE( l_coeffs1[3+0] == Approx(0.00527973491702664)  );
  REQUIRE( l_coeffs1[3+1] == Approx(0.000423653588369432) );
  REQUIRE( l_coeffs1[3+2] == Approx(0.00532171233332960)  );
}

TEST_CASE( "Derivation of the source matrix in 3D", "[ViscoElasticity][src]" ) {
  double l_srcElastic[3][6*6];
  double l_lameE[2] = {-1, -1};

  edge::seismic::setups::ViscoElasticity::src( 3,
                                               2,
                                               10,
                                               2000*0.1,
                                               2000*0.1*1.5,
                                               20.8E9,
                                               10.4E9,
                                               l_lameE[0],
                                               l_lameE[1],
                                               l_srcElastic );

  REQUIRE( l_lameE[0] == Approx(21031265249.9435) );
  REQUIRE( l_lameE[1] == Approx(10457744636.3789) );

  REQUIRE( l_srcElastic[0][ 0] == Approx(-330629718.304332) );
  REQUIRE( l_srcElastic[0][ 1] == Approx(-220201479.284256) );
  REQUIRE( l_srcElastic[0][ 2] == Approx(-220201479.284256) );
  REQUIRE( l_srcElastic[0][ 3] == Approx(0                ) );
  REQUIRE( l_srcElastic[0][ 4] == Approx(0                ) );
  REQUIRE( l_srcElastic[0][ 5] == Approx(0                ) );

  REQUIRE( l_srcElastic[0][ 6] == Approx(-220201479.284256) );
  REQUIRE( l_srcElastic[0][ 7] == Approx(-330629718.304332) );
  REQUIRE( l_srcElastic[0][ 8] == Approx(-220201479.284256) );
  REQUIRE( l_srcElastic[0][ 9] == Approx(0                ) );
  REQUIRE( l_srcElastic[0][10] == Approx(0                ) );
  REQUIRE( l_srcElastic[0][11] == Approx(0                ) );

  REQUIRE( l_srcElastic[0][12] == Approx(-220201479.284256) );
  REQUIRE( l_srcElastic[0][13] == Approx(-220201479.284256) );
  REQUIRE( l_srcElastic[0][14] == Approx(-330629718.304332) );
  REQUIRE( l_srcElastic[0][15] == Approx(0                ) );
  REQUIRE( l_srcElastic[0][16] == Approx(0                ) );
  REQUIRE( l_srcElastic[0][17] == Approx(0                ) );

  REQUIRE( l_srcElastic[0][18] == Approx(0                ) );
  REQUIRE( l_srcElastic[0][19] == Approx(0                ) );
  REQUIRE( l_srcElastic[0][20] == Approx(0                ) );
  REQUIRE( l_srcElastic[0][21] == Approx(-110428239.020076) );
  REQUIRE( l_srcElastic[0][22] == Approx(0                ) );
  REQUIRE( l_srcElastic[0][23] == Approx(0                ) );

  REQUIRE( l_srcElastic[0][24] == Approx(0                ) );
  REQUIRE( l_srcElastic[0][25] == Approx(0                ) );
  REQUIRE( l_srcElastic[0][26] == Approx(0                ) );
  REQUIRE( l_srcElastic[0][27] == Approx(0                ) );
  REQUIRE( l_srcElastic[0][28] == Approx(-110428239.020076) );
  REQUIRE( l_srcElastic[0][29] == Approx(0                ) );

  REQUIRE( l_srcElastic[0][30] == Approx(0                ) );
  REQUIRE( l_srcElastic[0][31] == Approx(0                ) );
  REQUIRE( l_srcElastic[0][32] == Approx(0                ) );
  REQUIRE( l_srcElastic[0][33] == Approx(0                ) );
  REQUIRE( l_srcElastic[0][34] == Approx(0                ) );
  REQUIRE( l_srcElastic[0][35] == Approx(-110428239.020076) );


  REQUIRE( l_srcElastic[1][ 0] == Approx(-26574526.2087753) );
  REQUIRE( l_srcElastic[1][ 1] == Approx(-17713604.1258692) );
  REQUIRE( l_srcElastic[1][ 2] == Approx(-17713604.1258692) );
  REQUIRE( l_srcElastic[1][ 3] == Approx(0                ) );
  REQUIRE( l_srcElastic[1][ 4] == Approx(0                ) );
  REQUIRE( l_srcElastic[1][ 5] == Approx(0                ) );

  REQUIRE( l_srcElastic[1][ 6] == Approx(-17713604.1258692) );
  REQUIRE( l_srcElastic[1][ 7] == Approx(-26574526.2087753) );
  REQUIRE( l_srcElastic[1][ 8] == Approx(-17713604.1258692) );
  REQUIRE( l_srcElastic[1][ 9] == Approx(0                ) );
  REQUIRE( l_srcElastic[1][10] == Approx(0                ) );
  REQUIRE( l_srcElastic[1][11] == Approx(0                ) );

  REQUIRE( l_srcElastic[1][12] == Approx(-17713604.1258692) );
  REQUIRE( l_srcElastic[1][13] == Approx(-17713604.1258692) );
  REQUIRE( l_srcElastic[1][14] == Approx(-26574526.2087753) );
  REQUIRE( l_srcElastic[1][15] == Approx(0                ) );
  REQUIRE( l_srcElastic[1][16] == Approx(0                ) );
  REQUIRE( l_srcElastic[1][17] == Approx(0                ) );

  REQUIRE( l_srcElastic[1][18] == Approx(0                ) );
  REQUIRE( l_srcElastic[1][19] == Approx(0                ) );
  REQUIRE( l_srcElastic[1][20] == Approx(0                ) );
  REQUIRE( l_srcElastic[1][21] == Approx(-8860922.08290611) );
  REQUIRE( l_srcElastic[1][22] == Approx(0                ) );
  REQUIRE( l_srcElastic[1][23] == Approx(0                ) );

  REQUIRE( l_srcElastic[1][24] == Approx(0                ) );
  REQUIRE( l_srcElastic[1][25] == Approx(0                ) );
  REQUIRE( l_srcElastic[1][26] == Approx(0                ) );
  REQUIRE( l_srcElastic[1][27] == Approx(0                ) );
  REQUIRE( l_srcElastic[1][28] == Approx(-8860922.08290611) );
  REQUIRE( l_srcElastic[1][29] == Approx(0                ) );

  REQUIRE( l_srcElastic[1][30] == Approx(0                ) );
  REQUIRE( l_srcElastic[1][31] == Approx(0                ) );
  REQUIRE( l_srcElastic[1][32] == Approx(0                ) );
  REQUIRE( l_srcElastic[1][33] == Approx(0                ) );
  REQUIRE( l_srcElastic[1][34] == Approx(0                ) );
  REQUIRE( l_srcElastic[1][35] == Approx(-8860922.08290611) );


  REQUIRE( l_srcElastic[2][ 0] == Approx(-334580569.994195) );
  REQUIRE( l_srcElastic[2][ 1] == Approx(-223274352.773737) );
  REQUIRE( l_srcElastic[2][ 2] == Approx(-223274352.773737) );
  REQUIRE( l_srcElastic[2][ 3] == Approx(0                ) );
  REQUIRE( l_srcElastic[2][ 4] == Approx(0                ) );
  REQUIRE( l_srcElastic[2][ 5] == Approx(0                ) );

  REQUIRE( l_srcElastic[2][ 6] == Approx(-223274352.773737) );
  REQUIRE( l_srcElastic[2][ 7] == Approx(-334580569.994195) );
  REQUIRE( l_srcElastic[2][ 8] == Approx(-223274352.773737) );
  REQUIRE( l_srcElastic[2][ 9] == Approx(0                ) );
  REQUIRE( l_srcElastic[2][10] == Approx(0                ) );
  REQUIRE( l_srcElastic[2][11] == Approx(0                ) );

  REQUIRE( l_srcElastic[2][12] == Approx(-223274352.773737) );
  REQUIRE( l_srcElastic[2][13] == Approx(-223274352.773737) );
  REQUIRE( l_srcElastic[2][14] == Approx(-334580569.994195) );
  REQUIRE( l_srcElastic[2][15] == Approx(0                ) );
  REQUIRE( l_srcElastic[2][16] == Approx(0                ) );
  REQUIRE( l_srcElastic[2][17] == Approx(0                ) );

  REQUIRE( l_srcElastic[2][18] == Approx(0                ) );
  REQUIRE( l_srcElastic[2][19] == Approx(0                ) );
  REQUIRE( l_srcElastic[2][20] == Approx(0                ) );
  REQUIRE( l_srcElastic[2][21] == Approx(-111306217.220458) );
  REQUIRE( l_srcElastic[2][22] == Approx(0                ) );
  REQUIRE( l_srcElastic[2][23] == Approx(0                ) );

  REQUIRE( l_srcElastic[2][24] == Approx(0                ) );
  REQUIRE( l_srcElastic[2][25] == Approx(0                ) );
  REQUIRE( l_srcElastic[2][26] == Approx(0                ) );
  REQUIRE( l_srcElastic[2][27] == Approx(0                ) );
  REQUIRE( l_srcElastic[2][28] == Approx(-111306217.220458) );
  REQUIRE( l_srcElastic[2][29] == Approx(0                ) );

  REQUIRE( l_srcElastic[2][30] == Approx(0                ) );
  REQUIRE( l_srcElastic[2][31] == Approx(0                ) );
  REQUIRE( l_srcElastic[2][32] == Approx(0                ) );
  REQUIRE( l_srcElastic[2][33] == Approx(0                ) );
  REQUIRE( l_srcElastic[2][34] == Approx(0                ) );
  REQUIRE( l_srcElastic[2][35] == Approx(-111306217.220458) );
}

TEST_CASE( "Derivation of the viscoelastic part of the star matrix in 3D", "[ViscoElasticity][starMatrix3d]" ) {
  float l_jacInv[3][3] = { { 4, 8, 0 },
                           { 2, 8, 4 },
                           { 8, 0, 4 } };
  float l_starA[3][6 * 3];

  edge::seismic::setups::ViscoElasticity::star( l_jacInv,
                                                l_starA );

  REQUIRE( l_starA[0][0*3 + 6 - 6] == Approx( -1.0 * 4 ) );
  REQUIRE( l_starA[0][3*3 + 7 - 6] == Approx( -0.5 * 4 ) );
  REQUIRE( l_starA[0][5*3 + 8 - 6] == Approx( -0.5 * 4 ) );

  REQUIRE( l_starA[0][1*3 + 7 - 6] == Approx( -1.0 * 2 ) );
  REQUIRE( l_starA[0][3*3 + 6 - 6] == Approx( -0.5 * 2 ) );
  REQUIRE( l_starA[0][4*3 + 8 - 6] == Approx( -0.5 * 2 ) );

  REQUIRE( l_starA[0][2*3 + 8 - 6] == Approx( -1.0 * 8 ) );
  REQUIRE( l_starA[0][4*3 + 7 - 6] == Approx( -0.5 * 8 ) );
  REQUIRE( l_starA[0][5*3 + 6 - 6] == Approx( -0.5 * 8 ) );
}

TEST_CASE( "Derivation of the viscoelastic part of the flux solvers in 3D", "[ViscoElasticity][fluxSolver3d]" ) {
  double l_rhoL1 =  1.2E3;
  double l_rhoR1 =  1.1E3;
  double l_lamL1 = 20.8E9;
  double l_lamR1 = 19.4E9;
  double l_muL1  = 10.4E9;
  double l_muR1  = 13.1E9;

  float l_fsMidL1[6][9];
  float l_fsMidR1[6][9];

  edge::seismic::setups::ViscoElasticity::fsMid( l_rhoL1, l_rhoR1,
                                                 l_lamL1, l_lamR1,
                                                 l_muL1,  l_muR1,
                                                 l_fsMidL1,
                                                 l_fsMidR1 );

  REQUIRE( l_fsMidL1[0][0] == Approx( 7.06824616167993E-8) );
  REQUIRE( l_fsMidL1[0][1] == 0                            );
  REQUIRE( l_fsMidL1[0][2] == 0                            );
  REQUIRE( l_fsMidL1[0][3] == 0                            );
  REQUIRE( l_fsMidL1[0][4] == 0                            );
  REQUIRE( l_fsMidL1[0][5] == 0                            );
  REQUIRE( l_fsMidL1[0][6] == Approx(-0.499400478754375  ) );
  REQUIRE( l_fsMidL1[0][7] == 0                            );
  REQUIRE( l_fsMidL1[0][8] == 0                            );

  for( unsigned short l_co = 0; l_co < 9; l_co++ ) {
    REQUIRE( l_fsMidL1[1][l_co] == 0 );
    REQUIRE( l_fsMidL1[2][l_co] == 0 );
  }

  REQUIRE( l_fsMidL1[3][0] == 0                            );
  REQUIRE( l_fsMidL1[3][1] == 0                            );
  REQUIRE( l_fsMidL1[3][2] == 0                            );
  REQUIRE( l_fsMidL1[3][3] == Approx( 6.82244126138142E-8) );
  REQUIRE( l_fsMidL1[3][4] == 0                            );
  REQUIRE( l_fsMidL1[3][5] == 0                            );
  REQUIRE( l_fsMidL1[3][6] == 0                            );
  REQUIRE( l_fsMidL1[3][7] == Approx(-0.241016678980355)   );
  REQUIRE( l_fsMidL1[3][8] == 0                            );

  for( unsigned short l_co = 0; l_co < 9; l_co++ ) {
    REQUIRE( l_fsMidL1[4][l_co] == 0 );
  }

  REQUIRE( l_fsMidL1[5][0] == 0                            );
  REQUIRE( l_fsMidL1[5][1] == 0                            );
  REQUIRE( l_fsMidL1[5][2] == 0                            );
  REQUIRE( l_fsMidL1[5][3] == 0                            );
  REQUIRE( l_fsMidL1[5][4] == 0                            );
  REQUIRE( l_fsMidL1[5][5] == Approx(6.82244126138142E-8)  );
  REQUIRE( l_fsMidL1[5][6] == 0                            );
  REQUIRE( l_fsMidL1[5][7] == 0                            );
  REQUIRE( l_fsMidL1[5][8] == Approx(-0.241016678980355)   );



  REQUIRE( l_fsMidR1[0][0] == Approx(-7.06824616167993E-8) );
  REQUIRE( l_fsMidR1[0][1] == 0                            );
  REQUIRE( l_fsMidR1[0][2] == 0                            );
  REQUIRE( l_fsMidR1[0][3] == 0                            );
  REQUIRE( l_fsMidR1[0][4] == 0                            );
  REQUIRE( l_fsMidR1[0][5] == 0                            );
  REQUIRE( l_fsMidR1[0][6] == Approx(-0.500599521245625  ) );
  REQUIRE( l_fsMidR1[0][7] == 0                            );
  REQUIRE( l_fsMidR1[0][8] == 0                            );

  for( unsigned short l_co = 0; l_co < 9; l_co++ ) {
    REQUIRE( l_fsMidL1[1][l_co] == 0 );
    REQUIRE( l_fsMidL1[2][l_co] == 0 );
  }

  REQUIRE( l_fsMidR1[3][0] == 0                            );
  REQUIRE( l_fsMidR1[3][1] == 0                            );
  REQUIRE( l_fsMidR1[3][2] == 0                            );
  REQUIRE( l_fsMidR1[3][3] == Approx(-6.82244126138142E-8) );
  REQUIRE( l_fsMidR1[3][4] == 0                            );
  REQUIRE( l_fsMidR1[3][5] == 0                            );
  REQUIRE( l_fsMidR1[3][6] == 0                            );
  REQUIRE( l_fsMidR1[3][7] == Approx(-0.258983321019645)   );
  REQUIRE( l_fsMidR1[3][8] == 0                            );

  for( unsigned short l_co = 0; l_co < 9; l_co++ ) {
    REQUIRE( l_fsMidR1[4][l_co] == 0 );
  }

  REQUIRE( l_fsMidR1[5][0] == 0                            );
  REQUIRE( l_fsMidR1[5][1] == 0                            );
  REQUIRE( l_fsMidR1[5][2] == 0                            );
  REQUIRE( l_fsMidR1[5][3] == 0                            );
  REQUIRE( l_fsMidR1[5][4] == 0                            );
  REQUIRE( l_fsMidR1[5][5] == Approx(-6.82244126138142E-8)  );
  REQUIRE( l_fsMidR1[5][6] == 0                            );
  REQUIRE( l_fsMidR1[5][7] == 0                            );
  REQUIRE( l_fsMidR1[5][8] == Approx(-0.258983321019645)   );
}