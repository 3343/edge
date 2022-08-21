/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2019, Alexander Breuer
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
 * Reference solutions for the elastic data generation.
 **/
// python code of the reference solution
/*
  import sympy
  import numpy as np
  sympy.init_printing(use_unicode=True)

  def printCStyle( i_mat ):
    print( '{' )
    for l_ro in i_mat:
      print( '  {', end="" )
      for l_co in l_ro:
        print( l_co, end="," )
      print( '},' )
    print( '},' )

  # Lame parameters + density
  la = sympy.symbols('lambda', real=True)
  mu = sympy.symbols('mu', real=True)
  rh = sympy.symbols('rho', real=True)

  # elastic Jacobian matrix
  Ae = sympy.Matrix(
  [
    [ 0, 0,  0, -(la+2*mu), 0],
    [ 0, 0,  0, -la, 0],
    
    [ 0, 0,  0, 0, -mu],

    [ -1/rh, 0, 0,  0, 0],
    [ 0, 0,  -1/rh, 0, 0],
  ]
  )

  Be = sympy.Matrix(
  [
    [ 0, 0,  0, 0, -la],
    [ 0, 0,  0, 0, -(la+2*mu)],
    
    [ 0, 0,  0, -mu, 0],

    [ 0, 0, -1/rh,  0, 0],
    [ 0, -1/rh,  0, 0, 0],
  ]
  )

  AeNum = Ae.subs( rh, 2600 ).subs( la, 20.8E9 ).subs( mu, 9.9E9 )
  BeNum = Be.subs( rh, 2600 ).subs( la, 20.8E9 ).subs( mu, 9.9E9 )

  AeStar = 0.4 * AeNum + 0.5 * BeNum
  BeStar = 0.3 * AeNum + 0.1 * BeNum
  print('AeStar')
  printCStyle( np.array(AeStar) )
  print('BeStar')
  printCStyle( np.array(BeStar) )


  Ae = sympy.Matrix([
    [0, 0, 0, 0, 0, 0, -(la+2*mu), 0, 0],
    [0, 0, 0, 0, 0, 0, -la, 0, 0],
    [0, 0, 0, 0, 0, 0, -la, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, -mu, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, -mu],
    [-1/rh, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, -1/rh, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, -1/rh, 0, 0, 0]])

  Be = sympy.Matrix([
    [0, 0, 0, 0, 0, 0, 0, -la, 0],
    [0, 0, 0, 0, 0, 0, 0, -(la+2*mu), 0],
    [0, 0, 0, 0, 0, 0, 0, -la, 0],
    [0, 0, 0, 0, 0, 0, -mu, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, -mu],
    [0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, -1/rh, 0, 0, 0, 0, 0],
    [0, -1/rh, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, -1/rh, 0, 0, 0, 0]])

  Ce = sympy.Matrix([
    [0, 0, 0, 0, 0, 0, 0, 0, -la],
    [0, 0, 0, 0, 0, 0, 0, 0, -la],
    [0, 0, 0, 0, 0, 0, 0, 0, -(la+2*mu)],
    [0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, -mu, 0],
    [0, 0, 0, 0, 0, 0, -mu, 0, 0],
    [0, 0, 0, 0, 0, -1/rh, 0, 0, 0],
    [0, 0, 0, 0, -1/rh, 0, 0, 0, 0],
    [0, 0, -1/rh, 0, 0, 0, 0, 0, 0]])

  AeNum = Ae.subs( rh, 2600 ).subs( la, 20.8E9 ).subs( mu, 9.9E9 )
  BeNum = Be.subs( rh, 2600 ).subs( la, 20.8E9 ).subs( mu, 9.9E9 )
  CeNum = Ce.subs( rh, 2600 ).subs( la, 20.8E9 ).subs( mu, 9.9E9 )

  AeStar = 0.4 * AeNum + 0.5 * BeNum + 0.2 * CeNum
  BeStar = 0.3 * AeNum + 0.1 * BeNum + 0.4 * CeNum
  CeStar = 0.5 * AeNum + 0.3 * BeNum + 0.2 * CeNum

  print('AeStar')
  printCStyle( np.array(AeStar) )
  print('BeStar')
  printCStyle( np.array(BeStar) )
  print('CeStar')
  printCStyle( np.array(CeStar) )
*/
#pragma GCC diagnostic push
#if !defined(__clang__)
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#pragma GCC diagnostic ignored "-Wunused-variable"
#endif

double l_starE2d[2][5][5] =
{
  {
    {0,0,0,-16240000000.0000,-10400000000.0000,},
    {0,0,0,-8320000000.00000,-20300000000.0000,},
    {0,0,0,-4950000000.00000,-3960000000.00000,},
    {-0.000153846153846154,0,-0.000192307692307692,0,0,},
    {0,-0.000192307692307692,-0.000153846153846154,0,0,},
  },
  {
    {0,0,0,-12180000000.0000,-2080000000.00000,},
    {0,0,0,-6240000000.00000,-4060000000.00000,},
    {0,0,0,-990000000.000000,-2970000000.00000,},
    {-0.000115384615384615,0,-3.84615384615385e-5,0,0,},
    {0,-3.84615384615385e-5,-0.000115384615384615,0,0,},
  }
};


double l_starE3d[3][9][9] =
{
  {
    {0,0,0,0,0,0,-16240000000.0000,-10400000000.0000,-4160000000.00000,},
    {0,0,0,0,0,0,-8320000000.00000,-20300000000.0000,-4160000000.00000,},
    {0,0,0,0,0,0,-8320000000.00000,-10400000000.0000,-8120000000.00000,},
    {0,0,0,0,0,0,-4950000000.00000,-3960000000.00000,0,},
    {0,0,0,0,0,0,0,-1980000000.00000,-4950000000.00000,},
    {0,0,0,0,0,0,-1980000000.00000,0,-3960000000.00000,},
    {-0.000153846153846154,0,0,-0.000192307692307692,0,-7.69230769230769e-5,0,0,0,},
    {0,-0.000192307692307692,0,-0.000153846153846154,-7.69230769230769e-5,0,0,0,0,},
    {0,0,-7.69230769230769e-5,0,-0.000192307692307692,-0.000153846153846154,0,0,0,},
  },
  {
    {0,0,0,0,0,0,-12180000000.0000,-2080000000.00000,-8320000000.00000,},
    {0,0,0,0,0,0,-6240000000.00000,-4060000000.00000,-8320000000.00000,},
    {0,0,0,0,0,0,-6240000000.00000,-2080000000.00000,-16240000000.0000,},
    {0,0,0,0,0,0,-990000000.000000,-2970000000.00000,0,},
    {0,0,0,0,0,0,0,-3960000000.00000,-990000000.000000,},
    {0,0,0,0,0,0,-3960000000.00000,0,-2970000000.00000,},
    {-0.000115384615384615,0,0,-3.84615384615385e-5,0,-0.000153846153846154,0,0,0,},
    {0,-3.84615384615385e-5,0,-0.000115384615384615,-0.000153846153846154,0,0,0,0,},
    {0,0,-0.000153846153846154,0,-3.84615384615385e-5,-0.000115384615384615,0,0,0,},
  },
  {
    {0,0,0,0,0,0,-20300000000.0000,-6240000000.00000,-4160000000.00000,},
    {0,0,0,0,0,0,-10400000000.0000,-12180000000.0000,-4160000000.00000,},
    {0,0,0,0,0,0,-10400000000.0000,-6240000000.00000,-8120000000.00000,},
    {0,0,0,0,0,0,-2970000000.00000,-4950000000.00000,0,},
    {0,0,0,0,0,0,0,-1980000000.00000,-2970000000.00000,},
    {0,0,0,0,0,0,-1980000000.00000,0,-4950000000.00000,},
    {-0.000192307692307692,0,0,-0.000115384615384615,0,-7.69230769230769e-5,0,0,0,},
    {0,-0.000115384615384615,0,-0.000192307692307692,-7.69230769230769e-5,0,0,0,0,},
    {0,0,-7.69230769230769e-5,0,-0.000115384615384615,-0.000192307692307692,0,0,0,},
  }
};

#pragma GCC diagnostic pop