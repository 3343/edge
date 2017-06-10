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
 * Generated tetrahedral basis.
 **/

assert( b < 286 );
if( b == 0 ) {
  val = 1 ;
}
if( b == 1 ) {
  val = 4*xi_3 - 1 ;
}
if( b == 2 ) {
  val = 3*xi_2 + xi_3 - 1 ;
}
if( b == 3 ) {
  val = 2*xi_1 + xi_2 + xi_3 - 1 ;
}
if( b == 4 ) {
  val = 15*std::pow(xi_3, 2) - 10*xi_3 + 1 ;
}
if( b == 5 ) {
  val = (6*xi_3 - 1)*(3*xi_2 + xi_3 - 1) ;
}
if( b == 6 ) {
  val = 10*std::pow(xi_2, 2) + 8*xi_2*xi_3 - 8*xi_2 + std::pow(xi_3, 2) - 2*xi_3 + 1 ;
}
if( b == 7 ) {
  val = (6*xi_3 - 1)*(2*xi_1 + xi_2 + xi_3 - 1) ;
}
if( b == 8 ) {
  val = (5*xi_2 + xi_3 - 1)*(2*xi_1 + xi_2 + xi_3 - 1) ;
}
if( b == 9 ) {
  val = -1.0L/2.0L*std::pow(xi_2 + xi_3 - 1, 2) + (3.0L/2.0L)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) ;
}
if( b == 10 ) {
  val = 56*std::pow(xi_3, 3) - 63*std::pow(xi_3, 2) + 18*xi_3 - 1 ;
}
if( b == 11 ) {
  val = (3*xi_2 + xi_3 - 1)*(14*xi_3 + 7*std::pow(2*xi_3 - 1, 2) - 6) ;
}
if( b == 12 ) {
  val = -(8*xi_3 - 1)*(4*xi_2*std::pow(xi_3 - 1, 2) + 3*std::pow(xi_3 - 1, 3) - 5*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2))/(2*xi_3 - 2) ;
}
if( b == 13 ) {
  val = 35*std::pow(xi_2, 3) + 45*std::pow(xi_2, 2)*xi_3 - 45*std::pow(xi_2, 2) + 15*xi_2*std::pow(xi_3, 2) - 30*xi_2*xi_3 + 15*xi_2 + std::pow(xi_3, 3) - 3*std::pow(xi_3, 2) + 3*xi_3 - 1 ;
}
if( b == 14 ) {
  val = (14*xi_3 + 7*std::pow(2*xi_3 - 1, 2) - 6)*(2*xi_1 + xi_2 + xi_3 - 1) ;
}
if( b == 15 ) {
  val = (8*xi_3 - 1)*(5*xi_2 + xi_3 - 1)*(2*xi_1 + xi_2 + xi_3 - 1) ;
}
if( b == 16 ) {
  val = -(36*xi_2*std::pow(xi_3 - 1, 2) + 17*std::pow(xi_3 - 1, 3) - 21*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2))*(2*xi_1 + xi_2 + xi_3 - 1)/(4*xi_3 - 4) ;
}
if( b == 17 ) {
  val = (1.0L/2.0L)*(8*xi_3 - 1)*(-std::pow(xi_2 + xi_3 - 1, 2) + 3*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2)) ;
}
if( b == 18 ) {
  val = (1.0L/2.0L)*(-std::pow(xi_2 + xi_3 - 1, 2) + 3*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2))*(7*xi_2 + xi_3 - 1) ;
}
if( b == 19 ) {
  val = -(6*xi_1*std::pow(xi_2 + xi_3 - 1, 3) + 3*std::pow(xi_2 + xi_3 - 1, 4) - 5*(xi_2 + xi_3 - 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3))/(2*xi_2 + 2*xi_3 - 2) ;
}
if( b == 20 ) {
  val = 210*std::pow(xi_3, 4) - 336*std::pow(xi_3, 3) + 168*std::pow(xi_3, 2) - 28*xi_3 + 1 ;
}
if( b == 21 ) {
  val = (3*xi_2 + xi_3 - 1)*(6*xi_3 + 15*std::pow(2*xi_3 - 1, 3) + 18*std::pow(2*xi_3 - 1, 2) - 4) ;
}
if( b == 22 ) {
  val = -(108*xi_3 + 45*std::pow(2*xi_3 - 1, 2) - 41)*(4*xi_2*std::pow(xi_3 - 1, 2) + 3*std::pow(xi_3 - 1, 3) - 5*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2))/(8*xi_3 - 8) ;
}
if( b == 23 ) {
  val = -(10*xi_3 - 1)*(30*xi_2*std::pow(xi_3 - 1, 3) + 12*std::pow(xi_3 - 1, 4) + 15*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 2) - 35*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 3))/(8*xi_3 - 8) ;
}
if( b == 24 ) {
  val = 3*xi_2*std::pow(xi_3 - 1, 3) + (15.0L/8.0L)*std::pow(xi_3 - 1, 4) - 21.0L/4.0L*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 2) - 7.0L/2.0L*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 3) + (63.0L/8.0L)*std::pow(2*xi_2 + xi_3 - 1, 4) ;
}
if( b == 25 ) {
  val = (2*xi_1 + xi_2 + xi_3 - 1)*(6*xi_3 + 15*std::pow(2*xi_3 - 1, 3) + 18*std::pow(2*xi_3 - 1, 2) - 4) ;
}
if( b == 26 ) {
  val = (1.0L/4.0L)*(5*xi_2 + xi_3 - 1)*(108*xi_3 + 45*std::pow(2*xi_3 - 1, 2) - 41)*(2*xi_1 + xi_2 + xi_3 - 1) ;
}
if( b == 27 ) {
  val = -(10*xi_3 - 1)*(36*xi_2*std::pow(xi_3 - 1, 2) + 17*std::pow(xi_3 - 1, 3) - 21*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2))*(2*xi_1 + xi_2 + xi_3 - 1)/(4*xi_3 - 4) ;
}
if( b == 28 ) {
  val = (1.0L/2.0L)*(2*std::pow(xi_3 - 1, 3) - 21*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2) + 21*std::pow(2*xi_2 + xi_3 - 1, 3))*(2*xi_1 + xi_2 + xi_3 - 1) ;
}
if( b == 29 ) {
  val = (1.0L/8.0L)*(-std::pow(xi_2 + xi_3 - 1, 2) + 3*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2))*(108*xi_3 + 45*std::pow(2*xi_3 - 1, 2) - 41) ;
}
if( b == 30 ) {
  val = (1.0L/2.0L)*(10*xi_3 - 1)*(-std::pow(xi_2 + xi_3 - 1, 2) + 3*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2))*(7*xi_2 + xi_3 - 1) ;
}
if( b == 31 ) {
  val = (std::pow(xi_2 + xi_3 - 1, 2) - 3*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2))*(20*xi_2*std::pow(xi_3 - 1, 2) + 8*std::pow(xi_3 - 1, 3) - 9*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2))/(2*xi_3 - 2) ;
}
if( b == 32 ) {
  val = -(10*xi_3 - 1)*(6*xi_1*std::pow(xi_2 + xi_3 - 1, 3) + 3*std::pow(xi_2 + xi_3 - 1, 4) - 5*(xi_2 + xi_3 - 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3))/(2*xi_2 + 2*xi_3 - 2) ;
}
if( b == 33 ) {
  val = (1.0L/2.0L)*(9*xi_2 + xi_3 - 1)*(6*xi_1*std::pow(xi_2 + xi_3 - 1, 3) + 5*(-xi_2 - xi_3 + 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3) + 3*std::pow(xi_2 + xi_3 - 1, 4))/(-xi_2 - xi_3 + 1) ;
}
if( b == 34 ) {
  val = (3.0L/8.0L)*std::pow(xi_2 + xi_3 - 1, 4) - 15.0L/4.0L*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + (35.0L/8.0L)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 4) ;
}
if( b == 35 ) {
  val = 792*std::pow(xi_3, 5) - 1650*std::pow(xi_3, 4) + 1200*std::pow(xi_3, 3) - 360*std::pow(xi_3, 2) + 40*xi_3 - 1 ;
}
if( b == 36 ) {
  val = (1.0L/16.0L)*(3*xi_2 + xi_3 - 1)*(-216*xi_3 + 495*std::pow(2*xi_3 - 1, 4) + 660*std::pow(2*xi_3 - 1, 3) + 90*std::pow(2*xi_3 - 1, 2) + 91) ;
}
if( b == 37 ) {
  val = -(4*xi_2*std::pow(xi_3 - 1, 2) + 3*std::pow(xi_3 - 1, 3) - 5*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2))*(120*xi_3 + 110*std::pow(2*xi_3 - 1, 3) + 165*std::pow(2*xi_3 - 1, 2) - 59)/(8*xi_3 - 8) ;
}
if( b == 38 ) {
  val = -(88*xi_3 + 33*std::pow(2*xi_3 - 1, 2) - 31)*(30*xi_2*std::pow(xi_3 - 1, 3) + 12*std::pow(xi_3 - 1, 4) + 15*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 2) - 35*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 3))/(16*xi_3 - 16) ;
}
if( b == 39 ) {
  val = (1.0L/8.0L)*(12*xi_3 - 1)*(24*xi_2*std::pow(xi_3 - 1, 3) + 15*std::pow(xi_3 - 1, 4) - 42*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 2) - 28*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 3) + 63*std::pow(2*xi_2 + xi_3 - 1, 4)) ;
}
if( b == 40 ) {
  val = (1.0L/16.0L)*(70*xi_2*std::pow(xi_3 - 1, 5) - 105*std::pow(-xi_3 + 1, 2)*std::pow(-2*xi_2 - xi_3 + 1, 4) + 231*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 5) + 30*std::pow(xi_3 - 1, 6) + 70*std::pow(xi_3 - 1, 4)*std::pow(-2*xi_2 - xi_3 + 1, 2) + 210*std::pow(xi_3 - 1, 3)*std::pow(-2*xi_2 - xi_3 + 1, 3))/(xi_3 - 1) ;
}
if( b == 41 ) {
  val = (1.0L/16.0L)*(2*xi_1 + xi_2 + xi_3 - 1)*(-216*xi_3 + 495*std::pow(2*xi_3 - 1, 4) + 660*std::pow(2*xi_3 - 1, 3) + 90*std::pow(2*xi_3 - 1, 2) + 91) ;
}
if( b == 42 ) {
  val = (1.0L/4.0L)*(5*xi_2 + xi_3 - 1)*(2*xi_1 + xi_2 + xi_3 - 1)*(120*xi_3 + 110*std::pow(2*xi_3 - 1, 3) + 165*std::pow(2*xi_3 - 1, 2) - 59) ;
}
if( b == 43 ) {
  val = -(88*xi_3 + 33*std::pow(2*xi_3 - 1, 2) - 31)*(36*xi_2*std::pow(xi_3 - 1, 2) + 17*std::pow(xi_3 - 1, 3) - 21*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2))*(2*xi_1 + xi_2 + xi_3 - 1)/(8*xi_3 - 8) ;
}
if( b == 44 ) {
  val = (1.0L/2.0L)*(12*xi_3 - 1)*(2*std::pow(xi_3 - 1, 3) - 21*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2) + 21*std::pow(2*xi_2 + xi_3 - 1, 3))*(2*xi_1 + xi_2 + xi_3 - 1) ;
}
if( b == 45 ) {
  val = (1.0L/8.0L)*(2*xi_1 + xi_2 + xi_3 - 1)*(88*xi_2*std::pow(xi_3 - 1, 3) + 41*std::pow(xi_3 - 1, 4) - 18*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 2) - 180*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 3) + 165*std::pow(2*xi_2 + xi_3 - 1, 4)) ;
}
if( b == 46 ) {
  val = (1.0L/8.0L)*(-std::pow(xi_2 + xi_3 - 1, 2) + 3*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2))*(120*xi_3 + 110*std::pow(2*xi_3 - 1, 3) + 165*std::pow(2*xi_3 - 1, 2) - 59) ;
}
if( b == 47 ) {
  val = (1.0L/4.0L)*(-std::pow(xi_2 + xi_3 - 1, 2) + 3*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2))*(7*xi_2 + xi_3 - 1)*(88*xi_3 + 33*std::pow(2*xi_3 - 1, 2) - 31) ;
}
if( b == 48 ) {
  val = (12*xi_3 - 1)*(std::pow(xi_2 + xi_3 - 1, 2) - 3*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2))*(20*xi_2*std::pow(xi_3 - 1, 2) + 8*std::pow(xi_3 - 1, 3) - 9*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2))/(2*xi_3 - 2) ;
}
if( b == 49 ) {
  val = (1.0L/16.0L)*(-std::pow(xi_2 + xi_3 - 1, 2) + 3*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2))*(126*xi_2*std::pow(xi_3 - 1, 3) - 225*std::pow(-xi_3 + 1, 2)*std::pow(-2*xi_2 - xi_3 + 1, 2) + 165*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 3) + 68*std::pow(xi_3 - 1, 4))/(xi_3 - 1) ;
}
if( b == 50 ) {
  val = -(88*xi_3 + 33*std::pow(2*xi_3 - 1, 2) - 31)*(6*xi_1*std::pow(xi_2 + xi_3 - 1, 3) + 3*std::pow(xi_2 + xi_3 - 1, 4) - 5*(xi_2 + xi_3 - 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3))/(4*xi_2 + 4*xi_3 - 4) ;
}
if( b == 51 ) {
  val = (1.0L/2.0L)*(12*xi_3 - 1)*(9*xi_2 + xi_3 - 1)*(6*xi_1*std::pow(xi_2 + xi_3 - 1, 3) + 5*(-xi_2 - xi_3 + 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3) + 3*std::pow(xi_2 + xi_3 - 1, 4))/(-xi_2 - xi_3 + 1) ;
}
if( b == 52 ) {
  val = (1.0L/8.0L)*(6*xi_1*std::pow(xi_2 + xi_3 - 1, 3) + 5*(-xi_2 - xi_3 + 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3) + 3*std::pow(xi_2 + xi_3 - 1, 4))*(140*xi_2*std::pow(xi_3 - 1, 2) + 55*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 2) + 51*std::pow(xi_3 - 1, 3))/((xi_3 - 1)*(xi_2 + xi_3 - 1)) ;
}
if( b == 53 ) {
  val = (1.0L/8.0L)*(12*xi_3 - 1)*(3*std::pow(xi_2 + xi_3 - 1, 4) - 30*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + 35*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 4)) ;
}
if( b == 54 ) {
  val = (1.0L/8.0L)*(11*xi_2 + xi_3 - 1)*(3*std::pow(xi_2 + xi_3 - 1, 4) - 30*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + 35*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 4)) ;
}
if( b == 55 ) {
  val = (1.0L/8.0L)*(30*xi_1*std::pow(xi_2 + xi_3 - 1, 5) + 63*(-xi_2 - xi_3 + 1)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 5) + 15*std::pow(xi_2 + xi_3 - 1, 6) + 70*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 3))/(xi_2 + xi_3 - 1) ;
}
if( b == 56 ) {
  val = 3003*std::pow(xi_3, 6) - 7722*std::pow(xi_3, 5) + 7425*std::pow(xi_3, 4) - 3300*std::pow(xi_3, 3) + 675*std::pow(xi_3, 2) - 54*xi_3 + 1 ;
}
if( b == 57 ) {
  val = (1.0L/16.0L)*(3*xi_2 + xi_3 - 1)*(-190*xi_3 + 1001*std::pow(2*xi_3 - 1, 5) + 1430*std::pow(2*xi_3 - 1, 4) + 110*std::pow(2*xi_3 - 1, 3) - 440*std::pow(2*xi_3 - 1, 2) + 105) ;
}
if( b == 58 ) {
  val = -(4*xi_2*std::pow(xi_3 - 1, 2) + 3*std::pow(xi_3 - 1, 3) - 5*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2))*(-88*xi_3 + 1001*std::pow(2*xi_3 - 1, 4) + 1716*std::pow(2*xi_3 - 1, 3) + 726*std::pow(2*xi_3 - 1, 2) + 5)/(32*xi_3 - 32) ;
}
if( b == 59 ) {
  val = -(150*xi_3 + 91*std::pow(2*xi_3 - 1, 3) + 156*std::pow(2*xi_3 - 1, 2) - 67)*(30*xi_2*std::pow(xi_3 - 1, 3) + 12*std::pow(xi_3 - 1, 4) + 15*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 2) - 35*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 3))/(16*xi_3 - 16) ;
}
if( b == 60 ) {
  val = (1.0L/32.0L)*(260*xi_3 + 91*std::pow(2*xi_3 - 1, 2) - 87)*(24*xi_2*std::pow(xi_3 - 1, 3) + 15*std::pow(xi_3 - 1, 4) - 42*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 2) - 28*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 3) + 63*std::pow(2*xi_2 + xi_3 - 1, 4)) ;
}
if( b == 61 ) {
  val = (1.0L/16.0L)*(14*xi_3 - 1)*(70*xi_2*std::pow(xi_3 - 1, 5) - 105*std::pow(-xi_3 + 1, 2)*std::pow(-2*xi_2 - xi_3 + 1, 4) + 231*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 5) + 30*std::pow(xi_3 - 1, 6) + 70*std::pow(xi_3 - 1, 4)*std::pow(-2*xi_2 - xi_3 + 1, 2) + 210*std::pow(xi_3 - 1, 3)*std::pow(-2*xi_2 - xi_3 + 1, 3))/(xi_3 - 1) ;
}
if( b == 62 ) {
  val = (1.0L/16.0L)*(-60*xi_2*std::pow(xi_3 - 1, 6) - 135*std::pow(-xi_3 + 1, 5)*std::pow(-2*xi_2 - xi_3 + 1, 2) + 198*std::pow(-xi_3 + 1, 2)*std::pow(-2*xi_2 - xi_3 + 1, 5) - 429*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 6) - 35*std::pow(xi_3 - 1, 7) - 180*std::pow(xi_3 - 1, 4)*std::pow(-2*xi_2 - xi_3 + 1, 3) - 495*std::pow(xi_3 - 1, 3)*std::pow(-2*xi_2 - xi_3 + 1, 4))/(xi_3 - 1) ;
}
if( b == 63 ) {
  val = (1.0L/16.0L)*(2*xi_1 + xi_2 + xi_3 - 1)*(-190*xi_3 + 1001*std::pow(2*xi_3 - 1, 5) + 1430*std::pow(2*xi_3 - 1, 4) + 110*std::pow(2*xi_3 - 1, 3) - 440*std::pow(2*xi_3 - 1, 2) + 105) ;
}
if( b == 64 ) {
  val = (1.0L/16.0L)*(5*xi_2 + xi_3 - 1)*(2*xi_1 + xi_2 + xi_3 - 1)*(-88*xi_3 + 1001*std::pow(2*xi_3 - 1, 4) + 1716*std::pow(2*xi_3 - 1, 3) + 726*std::pow(2*xi_3 - 1, 2) + 5) ;
}
if( b == 65 ) {
  val = -(36*xi_2*std::pow(xi_3 - 1, 2) + 17*std::pow(xi_3 - 1, 3) - 21*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2))*(2*xi_1 + xi_2 + xi_3 - 1)*(150*xi_3 + 91*std::pow(2*xi_3 - 1, 3) + 156*std::pow(2*xi_3 - 1, 2) - 67)/(8*xi_3 - 8) ;
}
if( b == 66 ) {
  val = (1.0L/8.0L)*(260*xi_3 + 91*std::pow(2*xi_3 - 1, 2) - 87)*(2*std::pow(xi_3 - 1, 3) - 21*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2) + 21*std::pow(2*xi_2 + xi_3 - 1, 3))*(2*xi_1 + xi_2 + xi_3 - 1) ;
}
if( b == 67 ) {
  val = (1.0L/8.0L)*(14*xi_3 - 1)*(2*xi_1 + xi_2 + xi_3 - 1)*(88*xi_2*std::pow(xi_3 - 1, 3) + 41*std::pow(xi_3 - 1, 4) - 18*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 2) - 180*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 3) + 165*std::pow(2*xi_2 + xi_3 - 1, 4)) ;
}
if( b == 68 ) {
  val = (1.0L/32.0L)*(2*xi_1 + xi_2 + xi_3 - 1)*(-90*xi_2*std::pow(xi_3 - 1, 5) - 70*std::pow(-xi_3 + 1, 6) - 1485*std::pow(-xi_3 + 1, 2)*std::pow(-2*xi_2 - xi_3 + 1, 4) + 1287*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 5) + 630*std::pow(xi_3 - 1, 4)*std::pow(-2*xi_2 - xi_3 + 1, 2) + 330*std::pow(xi_3 - 1, 3)*std::pow(-2*xi_2 - xi_3 + 1, 3))/(xi_3 - 1) ;
}
if( b == 69 ) {
  val = (1.0L/32.0L)*(-std::pow(xi_2 + xi_3 - 1, 2) + 3*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2))*(-88*xi_3 + 1001*std::pow(2*xi_3 - 1, 4) + 1716*std::pow(2*xi_3 - 1, 3) + 726*std::pow(2*xi_3 - 1, 2) + 5) ;
}
if( b == 70 ) {
  val = (1.0L/4.0L)*(-std::pow(xi_2 + xi_3 - 1, 2) + 3*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2))*(7*xi_2 + xi_3 - 1)*(150*xi_3 + 91*std::pow(2*xi_3 - 1, 3) + 156*std::pow(2*xi_3 - 1, 2) - 67) ;
}
if( b == 71 ) {
  val = (std::pow(xi_2 + xi_3 - 1, 2) - 3*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2))*(260*xi_3 + 91*std::pow(2*xi_3 - 1, 2) - 87)*(20*xi_2*std::pow(xi_3 - 1, 2) + 8*std::pow(xi_3 - 1, 3) - 9*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2))/(8*xi_3 - 8) ;
}
if( b == 72 ) {
  val = (1.0L/16.0L)*(14*xi_3 - 1)*(-std::pow(xi_2 + xi_3 - 1, 2) + 3*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2))*(126*xi_2*std::pow(xi_3 - 1, 3) - 225*std::pow(-xi_3 + 1, 2)*std::pow(-2*xi_2 - xi_3 + 1, 2) + 165*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 3) + 68*std::pow(xi_3 - 1, 4))/(xi_3 - 1) ;
}
if( b == 73 ) {
  val = (1.0L/32.0L)*(-std::pow(xi_2 + xi_3 - 1, 2) + 3*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2))*(200*xi_2*std::pow(xi_3 - 1, 3) + 71*std::pow(xi_3 - 1, 4) + 330*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 2) - 1100*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 3) + 715*std::pow(2*xi_2 + xi_3 - 1, 4)) ;
}
if( b == 74 ) {
  val = -(6*xi_1*std::pow(xi_2 + xi_3 - 1, 3) + 3*std::pow(xi_2 + xi_3 - 1, 4) - 5*(xi_2 + xi_3 - 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3))*(150*xi_3 + 91*std::pow(2*xi_3 - 1, 3) + 156*std::pow(2*xi_3 - 1, 2) - 67)/(4*xi_2 + 4*xi_3 - 4) ;
}
if( b == 75 ) {
  val = (1.0L/8.0L)*(9*xi_2 + xi_3 - 1)*(260*xi_3 + 91*std::pow(2*xi_3 - 1, 2) - 87)*(6*xi_1*std::pow(xi_2 + xi_3 - 1, 3) + 5*(-xi_2 - xi_3 + 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3) + 3*std::pow(xi_2 + xi_3 - 1, 4))/(-xi_2 - xi_3 + 1) ;
}
if( b == 76 ) {
  val = (1.0L/8.0L)*(14*xi_3 - 1)*(6*xi_1*std::pow(xi_2 + xi_3 - 1, 3) + 5*(-xi_2 - xi_3 + 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3) + 3*std::pow(xi_2 + xi_3 - 1, 4))*(140*xi_2*std::pow(xi_3 - 1, 2) + 55*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 2) + 51*std::pow(xi_3 - 1, 3))/((xi_3 - 1)*(xi_2 + xi_3 - 1)) ;
}
if( b == 77 ) {
  val = -(6*xi_1*std::pow(xi_2 + xi_3 - 1, 3) + 3*std::pow(xi_2 + xi_3 - 1, 4) - 5*(xi_2 + xi_3 - 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3))*(198*xi_2*std::pow(xi_3 - 1, 2) + 92*std::pow(xi_3 - 1, 3) - 231*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2) + 143*std::pow(2*xi_2 + xi_3 - 1, 3))/(8*xi_2 + 8*xi_3 - 8) ;
}
if( b == 78 ) {
  val = (1.0L/32.0L)*(260*xi_3 + 91*std::pow(2*xi_3 - 1, 2) - 87)*(3*std::pow(xi_2 + xi_3 - 1, 4) - 30*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + 35*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 4)) ;
}
if( b == 79 ) {
  val = (1.0L/8.0L)*(14*xi_3 - 1)*(11*xi_2 + xi_3 - 1)*(3*std::pow(xi_2 + xi_3 - 1, 4) - 30*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + 35*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 4)) ;
}
if( b == 80 ) {
  val = -(108*xi_2*std::pow(xi_3 - 1, 2) + 37*std::pow(xi_3 - 1, 3) - 39*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2))*(3*std::pow(xi_2 + xi_3 - 1, 4) - 30*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + 35*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 4))/(16*xi_3 - 16) ;
}
if( b == 81 ) {
  val = (1.0L/8.0L)*(14*xi_3 - 1)*(30*xi_1*std::pow(xi_2 + xi_3 - 1, 5) + 63*(-xi_2 - xi_3 + 1)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 5) + 15*std::pow(xi_2 + xi_3 - 1, 6) + 70*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 3))/(xi_2 + xi_3 - 1) ;
}
if( b == 82 ) {
  val = (1.0L/8.0L)*(13*xi_2 + xi_3 - 1)*(30*xi_1*std::pow(xi_2 + xi_3 - 1, 5) + 63*(-xi_2 - xi_3 + 1)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 5) + 15*std::pow(xi_2 + xi_3 - 1, 6) + 70*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 3))/(xi_2 + xi_3 - 1) ;
}
if( b == 83 ) {
  val = -5.0L/16.0L*std::pow(xi_2 + xi_3 - 1, 6) + (105.0L/16.0L)*std::pow(xi_2 + xi_3 - 1, 4)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) - 315.0L/16.0L*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 4) + (231.0L/16.0L)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 6) ;
}
if( b == 84 ) {
  val = 11440*std::pow(xi_3, 7) - 35035*std::pow(xi_3, 6) + 42042*std::pow(xi_3, 5) - 25025*std::pow(xi_3, 4) + 7700*std::pow(xi_3, 3) - 1155*std::pow(xi_3, 2) + 70*xi_3 - 1 ;
}
if( b == 85 ) {
  val = (1.0L/16.0L)*(3*xi_2 + xi_3 - 1)*(198*xi_3 + 2002*std::pow(2*xi_3 - 1, 6) + 3003*std::pow(2*xi_3 - 1, 5) - 1430*std::pow(2*xi_3 - 1, 3) - 330*std::pow(2*xi_3 - 1, 2) - 83) ;
}
if( b == 86 ) {
  val = -(4*xi_2*std::pow(xi_3 - 1, 2) + 3*std::pow(xi_3 - 1, 3) - 5*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2))*(-600*xi_3 + 2184*std::pow(2*xi_3 - 1, 5) + 4095*std::pow(2*xi_3 - 1, 4) + 1820*std::pow(2*xi_3 - 1, 3) - 390*std::pow(2*xi_3 - 1, 2) + 283)/(32*xi_3 - 32) ;
}
if( b == 87 ) {
  val = -(30*xi_2*std::pow(xi_3 - 1, 3) + 12*std::pow(xi_3 - 1, 4) + 15*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 2) - 35*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 3))*(156*xi_3 + 455*std::pow(2*xi_3 - 1, 4) + 910*std::pow(2*xi_3 - 1, 3) + 546*std::pow(2*xi_3 - 1, 2) - 87)/(32*xi_3 - 32) ;
}
if( b == 88 ) {
  val = (1.0L/32.0L)*(588*xi_3 + 280*std::pow(2*xi_3 - 1, 3) + 525*std::pow(2*xi_3 - 1, 2) - 249)*(24*xi_2*std::pow(xi_3 - 1, 3) + 15*std::pow(xi_3 - 1, 4) - 42*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 2) - 28*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 3) + 63*std::pow(2*xi_2 + xi_3 - 1, 4)) ;
}
if( b == 89 ) {
  val = (1.0L/16.0L)*(90*xi_3 + 30*std::pow(2*xi_3 - 1, 2) - 29)*(70*xi_2*std::pow(xi_3 - 1, 5) - 105*std::pow(-xi_3 + 1, 2)*std::pow(-2*xi_2 - xi_3 + 1, 4) + 231*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 5) + 30*std::pow(xi_3 - 1, 6) + 70*std::pow(xi_3 - 1, 4)*std::pow(-2*xi_2 - xi_3 + 1, 2) + 210*std::pow(xi_3 - 1, 3)*std::pow(-2*xi_2 - xi_3 + 1, 3))/(xi_3 - 1) ;
}
if( b == 90 ) {
  val = -(16*xi_3 - 1)*(60*xi_2*std::pow(xi_3 - 1, 6) + 35*std::pow(xi_3 - 1, 7) - 135*std::pow(xi_3 - 1, 5)*std::pow(2*xi_2 + xi_3 - 1, 2) - 180*std::pow(xi_3 - 1, 4)*std::pow(2*xi_2 + xi_3 - 1, 3) + 495*std::pow(xi_3 - 1, 3)*std::pow(2*xi_2 + xi_3 - 1, 4) + 198*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 5) - 429*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 6))/(16*xi_3 - 16) ;
}
if( b == 91 ) {
  val = (1.0L/128.0L)*(-630*xi_2*std::pow(xi_3 - 1, 7) - 3465*std::pow(-xi_3 + 1, 5)*std::pow(2*xi_2 + xi_3 - 1, 3) + 3465*std::pow(-xi_3 + 1, 4)*std::pow(2*xi_2 + xi_3 - 1, 4) - 6435*(-xi_3 + 1)*std::pow(2*xi_2 + xi_3 - 1, 7) - 280*std::pow(xi_3 - 1, 8) - 945*std::pow(xi_3 - 1, 6)*std::pow(2*xi_2 + xi_3 - 1, 2) - 9009*std::pow(xi_3 - 1, 3)*std::pow(2*xi_2 + xi_3 - 1, 5) - 3003*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 6))/(xi_3 - 1) ;
}
if( b == 92 ) {
  val = (1.0L/16.0L)*(2*xi_1 + xi_2 + xi_3 - 1)*(198*xi_3 + 2002*std::pow(2*xi_3 - 1, 6) + 3003*std::pow(2*xi_3 - 1, 5) - 1430*std::pow(2*xi_3 - 1, 3) - 330*std::pow(2*xi_3 - 1, 2) - 83) ;
}
if( b == 93 ) {
  val = (1.0L/16.0L)*(5*xi_2 + xi_3 - 1)*(2*xi_1 + xi_2 + xi_3 - 1)*(-600*xi_3 + 2184*std::pow(2*xi_3 - 1, 5) + 4095*std::pow(2*xi_3 - 1, 4) + 1820*std::pow(2*xi_3 - 1, 3) - 390*std::pow(2*xi_3 - 1, 2) + 283) ;
}
if( b == 94 ) {
  val = -(36*xi_2*std::pow(xi_3 - 1, 2) + 17*std::pow(xi_3 - 1, 3) - 21*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2))*(2*xi_1 + xi_2 + xi_3 - 1)*(156*xi_3 + 455*std::pow(2*xi_3 - 1, 4) + 910*std::pow(2*xi_3 - 1, 3) + 546*std::pow(2*xi_3 - 1, 2) - 87)/(16*xi_3 - 16) ;
}
if( b == 95 ) {
  val = (1.0L/8.0L)*(2*std::pow(xi_3 - 1, 3) - 21*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2) + 21*std::pow(2*xi_2 + xi_3 - 1, 3))*(2*xi_1 + xi_2 + xi_3 - 1)*(588*xi_3 + 280*std::pow(2*xi_3 - 1, 3) + 525*std::pow(2*xi_3 - 1, 2) - 249) ;
}
if( b == 96 ) {
  val = (1.0L/8.0L)*(90*xi_3 + 30*std::pow(2*xi_3 - 1, 2) - 29)*(2*xi_1 + xi_2 + xi_3 - 1)*(88*xi_2*std::pow(xi_3 - 1, 3) + 41*std::pow(xi_3 - 1, 4) - 18*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 2) - 180*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 3) + 165*std::pow(2*xi_2 + xi_3 - 1, 4)) ;
}
if( b == 97 ) {
  val = (1.0L/32.0L)*(16*xi_3 - 1)*(2*xi_1 + xi_2 + xi_3 - 1)*(-90*xi_2*std::pow(xi_3 - 1, 5) - 70*std::pow(-xi_3 + 1, 6) - 1485*std::pow(-xi_3 + 1, 2)*std::pow(-2*xi_2 - xi_3 + 1, 4) + 1287*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 5) + 630*std::pow(xi_3 - 1, 4)*std::pow(-2*xi_2 - xi_3 + 1, 2) + 330*std::pow(xi_3 - 1, 3)*std::pow(-2*xi_2 - xi_3 + 1, 3))/(xi_3 - 1) ;
}
if( b == 98 ) {
  val = -(2*xi_1 + xi_2 + xi_3 - 1)*(780*xi_2*std::pow(xi_3 - 1, 6) + 365*std::pow(xi_3 - 1, 7) + 165*std::pow(xi_3 - 1, 5)*std::pow(2*xi_2 + xi_3 - 1, 2) - 3740*std::pow(xi_3 - 1, 4)*std::pow(2*xi_2 + xi_3 - 1, 3) + 2145*std::pow(xi_3 - 1, 3)*std::pow(2*xi_2 + xi_3 - 1, 4) + 6006*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 5) - 5005*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 6))/(64*xi_3 - 64) ;
}
if( b == 99 ) {
  val = (1.0L/32.0L)*(-std::pow(xi_2 + xi_3 - 1, 2) + 3*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2))*(-600*xi_3 + 2184*std::pow(2*xi_3 - 1, 5) + 4095*std::pow(2*xi_3 - 1, 4) + 1820*std::pow(2*xi_3 - 1, 3) - 390*std::pow(2*xi_3 - 1, 2) + 283) ;
}
if( b == 100 ) {
  val = (1.0L/8.0L)*(-std::pow(xi_2 + xi_3 - 1, 2) + 3*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2))*(7*xi_2 + xi_3 - 1)*(156*xi_3 + 455*std::pow(2*xi_3 - 1, 4) + 910*std::pow(2*xi_3 - 1, 3) + 546*std::pow(2*xi_3 - 1, 2) - 87) ;
}
if( b == 101 ) {
  val = (std::pow(xi_2 + xi_3 - 1, 2) - 3*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2))*(20*xi_2*std::pow(xi_3 - 1, 2) + 8*std::pow(xi_3 - 1, 3) - 9*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2))*(588*xi_3 + 280*std::pow(2*xi_3 - 1, 3) + 525*std::pow(2*xi_3 - 1, 2) - 249)/(8*xi_3 - 8) ;
}
if( b == 102 ) {
  val = (1.0L/16.0L)*(-std::pow(xi_2 + xi_3 - 1, 2) + 3*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2))*(90*xi_3 + 30*std::pow(2*xi_3 - 1, 2) - 29)*(126*xi_2*std::pow(xi_3 - 1, 3) - 225*std::pow(-xi_3 + 1, 2)*std::pow(-2*xi_2 - xi_3 + 1, 2) + 165*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 3) + 68*std::pow(xi_3 - 1, 4))/(xi_3 - 1) ;
}
if( b == 103 ) {
  val = (1.0L/32.0L)*(16*xi_3 - 1)*(-std::pow(xi_2 + xi_3 - 1, 2) + 3*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2))*(200*xi_2*std::pow(xi_3 - 1, 3) + 71*std::pow(xi_3 - 1, 4) + 330*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 2) - 1100*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 3) + 715*std::pow(2*xi_2 + xi_3 - 1, 4)) ;
}
if( b == 104 ) {
  val = (std::pow(xi_2 + xi_3 - 1, 2) - 3*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2))*(770*xi_2*std::pow(xi_3 - 1, 5) + 386*std::pow(xi_3 - 1, 6) - 990*std::pow(xi_3 - 1, 4)*std::pow(2*xi_2 + xi_3 - 1, 2) - 1430*std::pow(xi_3 - 1, 3)*std::pow(2*xi_2 + xi_3 - 1, 3) + 5005*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 4) - 3003*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 5))/(64*xi_3 - 64) ;
}
if( b == 105 ) {
  val = -(6*xi_1*std::pow(xi_2 + xi_3 - 1, 3) + 3*std::pow(xi_2 + xi_3 - 1, 4) - 5*(xi_2 + xi_3 - 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3))*(156*xi_3 + 455*std::pow(2*xi_3 - 1, 4) + 910*std::pow(2*xi_3 - 1, 3) + 546*std::pow(2*xi_3 - 1, 2) - 87)/(8*xi_2 + 8*xi_3 - 8) ;
}
if( b == 106 ) {
  val = (1.0L/8.0L)*(9*xi_2 + xi_3 - 1)*(6*xi_1*std::pow(xi_2 + xi_3 - 1, 3) + 5*(-xi_2 - xi_3 + 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3) + 3*std::pow(xi_2 + xi_3 - 1, 4))*(588*xi_3 + 280*std::pow(2*xi_3 - 1, 3) + 525*std::pow(2*xi_3 - 1, 2) - 249)/(-xi_2 - xi_3 + 1) ;
}
if( b == 107 ) {
  val = (1.0L/8.0L)*(90*xi_3 + 30*std::pow(2*xi_3 - 1, 2) - 29)*(6*xi_1*std::pow(xi_2 + xi_3 - 1, 3) + 5*(-xi_2 - xi_3 + 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3) + 3*std::pow(xi_2 + xi_3 - 1, 4))*(140*xi_2*std::pow(xi_3 - 1, 2) + 55*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 2) + 51*std::pow(xi_3 - 1, 3))/((xi_3 - 1)*(xi_2 + xi_3 - 1)) ;
}
if( b == 108 ) {
  val = -(16*xi_3 - 1)*(6*xi_1*std::pow(xi_2 + xi_3 - 1, 3) + 3*std::pow(xi_2 + xi_3 - 1, 4) - 5*(xi_2 + xi_3 - 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3))*(198*xi_2*std::pow(xi_3 - 1, 2) + 92*std::pow(xi_3 - 1, 3) - 231*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2) + 143*std::pow(2*xi_2 + xi_3 - 1, 3))/(8*xi_2 + 8*xi_3 - 8) ;
}
if( b == 109 ) {
  val = (6*xi_1*std::pow(xi_2 + xi_3 - 1, 3) + 3*std::pow(xi_2 + xi_3 - 1, 4) - 5*(xi_2 + xi_3 - 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3))*(168*xi_2*std::pow(xi_3 - 1, 3) + 127*std::pow(xi_3 - 1, 4) - 1326*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 2) + 2548*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 3) - 1365*std::pow(2*xi_2 + xi_3 - 1, 4))/(32*xi_2 + 32*xi_3 - 32) ;
}
if( b == 110 ) {
  val = (1.0L/32.0L)*(3*std::pow(xi_2 + xi_3 - 1, 4) - 30*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + 35*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 4))*(588*xi_3 + 280*std::pow(2*xi_3 - 1, 3) + 525*std::pow(2*xi_3 - 1, 2) - 249) ;
}
if( b == 111 ) {
  val = (1.0L/8.0L)*(11*xi_2 + xi_3 - 1)*(90*xi_3 + 30*std::pow(2*xi_3 - 1, 2) - 29)*(3*std::pow(xi_2 + xi_3 - 1, 4) - 30*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + 35*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 4)) ;
}
if( b == 112 ) {
  val = -(16*xi_3 - 1)*(108*xi_2*std::pow(xi_3 - 1, 2) + 37*std::pow(xi_3 - 1, 3) - 39*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2))*(3*std::pow(xi_2 + xi_3 - 1, 4) - 30*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + 35*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 4))/(16*xi_3 - 16) ;
}
if( b == 113 ) {
  val = (1.0L/64.0L)*(3*std::pow(xi_2 + xi_3 - 1, 4) - 30*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + 35*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 4))*(858*xi_2*std::pow(xi_3 - 1, 3) - 819*std::pow(-xi_3 + 1, 2)*std::pow(-2*xi_2 - xi_3 + 1, 2) + 455*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 3) + 372*std::pow(xi_3 - 1, 4))/(xi_3 - 1) ;
}
if( b == 114 ) {
  val = (1.0L/8.0L)*(90*xi_3 + 30*std::pow(2*xi_3 - 1, 2) - 29)*(30*xi_1*std::pow(xi_2 + xi_3 - 1, 5) + 63*(-xi_2 - xi_3 + 1)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 5) + 15*std::pow(xi_2 + xi_3 - 1, 6) + 70*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 3))/(xi_2 + xi_3 - 1) ;
}
if( b == 115 ) {
  val = (1.0L/8.0L)*(16*xi_3 - 1)*(13*xi_2 + xi_3 - 1)*(30*xi_1*std::pow(xi_2 + xi_3 - 1, 5) + 63*(-xi_2 - xi_3 + 1)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 5) + 15*std::pow(xi_2 + xi_3 - 1, 6) + 70*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 3))/(xi_2 + xi_3 - 1) ;
}
if( b == 116 ) {
  val = -1.0L/32.0L*(308*xi_2*std::pow(xi_3 - 1, 2) + 101*std::pow(xi_3 - 1, 3) - 105*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2))*(30*xi_1*std::pow(xi_2 + xi_3 - 1, 5) + 15*std::pow(xi_2 + xi_3 - 1, 6) - 70*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3) + 63*(xi_2 + xi_3 - 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 5))/((xi_3 - 1)*(xi_2 + xi_3 - 1)) ;
}
if( b == 117 ) {
  val = (1.0L/16.0L)*(16*xi_3 - 1)*(-5*std::pow(xi_2 + xi_3 - 1, 6) + 105*std::pow(xi_2 + xi_3 - 1, 4)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) - 315*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 4) + 231*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 6)) ;
}
if( b == 118 ) {
  val = (1.0L/16.0L)*(15*xi_2 + xi_3 - 1)*(-5*std::pow(xi_2 + xi_3 - 1, 6) + 105*std::pow(xi_2 + xi_3 - 1, 4)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) - 315*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 4) + 231*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 6)) ;
}
if( b == 119 ) {
  val = -(70*xi_1*std::pow(xi_2 + xi_3 - 1, 7) + 35*std::pow(xi_2 + xi_3 - 1, 8) - 315*std::pow(xi_2 + xi_3 - 1, 5)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3) + 693*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 5) - 429*(xi_2 + xi_3 - 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 7))/(16*xi_2 + 16*xi_3 - 16) ;
}
if( b == 120 ) {
  val = 43758*std::pow(xi_3, 8) - 155584*std::pow(xi_3, 7) + 224224*std::pow(xi_3, 6) - 168168*std::pow(xi_3, 5) + 70070*std::pow(xi_3, 4) - 16016*std::pow(xi_3, 3) + 1848*std::pow(xi_3, 2) - 88*xi_3 + 1 ;
}
if( b == 121 ) {
  val = (1.0L/16.0L)*(3*xi_2 + xi_3 - 1)*(252*xi_3 + 3978*std::pow(2*xi_3 - 1, 7) + 6188*std::pow(2*xi_3 - 1, 6) - 546*std::pow(2*xi_3 - 1, 5) - 4095*std::pow(2*xi_3 - 1, 4) - 910*std::pow(2*xi_3 - 1, 3) + 546*std::pow(2*xi_3 - 1, 2) - 133) ;
}
if( b == 122 ) {
  val = -(4*xi_2*std::pow(xi_3 - 1, 2) + 3*std::pow(xi_3 - 1, 3) - 5*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2))*(-156*xi_3 + 4641*std::pow(2*xi_3 - 1, 6) + 9282*std::pow(2*xi_3 - 1, 5) + 4095*std::pow(2*xi_3 - 1, 4) - 1820*std::pow(2*xi_3 - 1, 3) - 1365*std::pow(2*xi_3 - 1, 2) + 107)/(32*xi_3 - 32) ;
}
if( b == 123 ) {
  val = -(30*xi_2*std::pow(xi_3 - 1, 3) + 12*std::pow(xi_3 - 1, 4) + 15*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 2) - 35*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 3))*(-210*xi_3 + 1071*std::pow(2*xi_3 - 1, 5) + 2380*std::pow(2*xi_3 - 1, 4) + 1610*std::pow(2*xi_3 - 1, 3) + 210*std::pow(2*xi_3 - 1, 2) + 87)/(32*xi_3 - 32) ;
}
if( b == 124 ) {
  val = (1.0L/32.0L)*(600*xi_3 + 765*std::pow(2*xi_3 - 1, 4) + 1700*std::pow(2*xi_3 - 1, 3) + 1230*std::pow(2*xi_3 - 1, 2) - 291)*(24*xi_2*std::pow(xi_3 - 1, 3) + 15*std::pow(xi_3 - 1, 4) - 42*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 2) - 28*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 3) + 63*std::pow(2*xi_2 + xi_3 - 1, 4)) ;
}
if( b == 125 ) {
  val = (1.0L/16.0L)*(252*xi_3 + 102*std::pow(2*xi_3 - 1, 3) + 204*std::pow(2*xi_3 - 1, 2) - 103)*(70*xi_2*std::pow(xi_3 - 1, 5) - 105*std::pow(-xi_3 + 1, 2)*std::pow(-2*xi_2 - xi_3 + 1, 4) + 231*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 5) + 30*std::pow(xi_3 - 1, 6) + 70*std::pow(xi_3 - 1, 4)*std::pow(-2*xi_2 - xi_3 + 1, 2) + 210*std::pow(xi_3 - 1, 3)*std::pow(-2*xi_2 - xi_3 + 1, 3))/(xi_3 - 1) ;
}
if( b == 126 ) {
  val = -(476*xi_3 + 153*std::pow(2*xi_3 - 1, 2) - 149)*(60*xi_2*std::pow(xi_3 - 1, 6) + 35*std::pow(xi_3 - 1, 7) - 135*std::pow(xi_3 - 1, 5)*std::pow(2*xi_2 + xi_3 - 1, 2) - 180*std::pow(xi_3 - 1, 4)*std::pow(2*xi_2 + xi_3 - 1, 3) + 495*std::pow(xi_3 - 1, 3)*std::pow(2*xi_2 + xi_3 - 1, 4) + 198*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 5) - 429*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 6))/(64*xi_3 - 64) ;
}
if( b == 127 ) {
  val = -(18*xi_3 - 1)*(630*xi_2*std::pow(xi_3 - 1, 7) + 280*std::pow(xi_3 - 1, 8) + 945*std::pow(xi_3 - 1, 6)*std::pow(2*xi_2 + xi_3 - 1, 2) - 3465*std::pow(xi_3 - 1, 5)*std::pow(2*xi_2 + xi_3 - 1, 3) - 3465*std::pow(xi_3 - 1, 4)*std::pow(2*xi_2 + xi_3 - 1, 4) + 9009*std::pow(xi_3 - 1, 3)*std::pow(2*xi_2 + xi_3 - 1, 5) + 3003*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 6) - 6435*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 7))/(128*xi_3 - 128) ;
}
if( b == 128 ) {
  val = (35.0L/8.0L)*xi_2*std::pow(xi_3 - 1, 7) + (385.0L/32.0L)*(-xi_3 + 1)*std::pow(xi_3 - 1, 5)*std::pow(2*xi_2 + xi_3 - 1, 2) + (385.0L/16.0L)*(-xi_3 + 1)*std::pow(xi_3 - 1, 4)*std::pow(2*xi_2 + xi_3 - 1, 3) + (5005.0L/32.0L)*(-xi_3 + 1)*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 6) + (715.0L/16.0L)*(-xi_3 + 1)*std::pow(2*xi_2 + xi_3 - 1, 7) + (315.0L/128.0L)*std::pow(xi_3 - 1, 8) + (5005.0L/64.0L)*std::pow(xi_3 - 1, 4)*std::pow(2*xi_2 + xi_3 - 1, 4) + (1001.0L/16.0L)*std::pow(xi_3 - 1, 3)*std::pow(2*xi_2 + xi_3 - 1, 5) + (12155.0L/128.0L)*std::pow(2*xi_2 + xi_3 - 1, 8) ;
}
if( b == 129 ) {
  val = (1.0L/16.0L)*(2*xi_1 + xi_2 + xi_3 - 1)*(252*xi_3 + 3978*std::pow(2*xi_3 - 1, 7) + 6188*std::pow(2*xi_3 - 1, 6) - 546*std::pow(2*xi_3 - 1, 5) - 4095*std::pow(2*xi_3 - 1, 4) - 910*std::pow(2*xi_3 - 1, 3) + 546*std::pow(2*xi_3 - 1, 2) - 133) ;
}
if( b == 130 ) {
  val = (1.0L/16.0L)*(5*xi_2 + xi_3 - 1)*(2*xi_1 + xi_2 + xi_3 - 1)*(-156*xi_3 + 4641*std::pow(2*xi_3 - 1, 6) + 9282*std::pow(2*xi_3 - 1, 5) + 4095*std::pow(2*xi_3 - 1, 4) - 1820*std::pow(2*xi_3 - 1, 3) - 1365*std::pow(2*xi_3 - 1, 2) + 107) ;
}
if( b == 131 ) {
  val = -(36*xi_2*std::pow(xi_3 - 1, 2) + 17*std::pow(xi_3 - 1, 3) - 21*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2))*(2*xi_1 + xi_2 + xi_3 - 1)*(-210*xi_3 + 1071*std::pow(2*xi_3 - 1, 5) + 2380*std::pow(2*xi_3 - 1, 4) + 1610*std::pow(2*xi_3 - 1, 3) + 210*std::pow(2*xi_3 - 1, 2) + 87)/(16*xi_3 - 16) ;
}
if( b == 132 ) {
  val = (1.0L/8.0L)*(2*std::pow(xi_3 - 1, 3) - 21*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2) + 21*std::pow(2*xi_2 + xi_3 - 1, 3))*(2*xi_1 + xi_2 + xi_3 - 1)*(600*xi_3 + 765*std::pow(2*xi_3 - 1, 4) + 1700*std::pow(2*xi_3 - 1, 3) + 1230*std::pow(2*xi_3 - 1, 2) - 291) ;
}
if( b == 133 ) {
  val = (1.0L/8.0L)*(2*xi_1 + xi_2 + xi_3 - 1)*(252*xi_3 + 102*std::pow(2*xi_3 - 1, 3) + 204*std::pow(2*xi_3 - 1, 2) - 103)*(88*xi_2*std::pow(xi_3 - 1, 3) + 41*std::pow(xi_3 - 1, 4) - 18*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 2) - 180*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 3) + 165*std::pow(2*xi_2 + xi_3 - 1, 4)) ;
}
if( b == 134 ) {
  val = (1.0L/128.0L)*(476*xi_3 + 153*std::pow(2*xi_3 - 1, 2) - 149)*(2*xi_1 + xi_2 + xi_3 - 1)*(-90*xi_2*std::pow(xi_3 - 1, 5) - 70*std::pow(-xi_3 + 1, 6) - 1485*std::pow(-xi_3 + 1, 2)*std::pow(-2*xi_2 - xi_3 + 1, 4) + 1287*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 5) + 630*std::pow(xi_3 - 1, 4)*std::pow(-2*xi_2 - xi_3 + 1, 2) + 330*std::pow(xi_3 - 1, 3)*std::pow(-2*xi_2 - xi_3 + 1, 3))/(xi_3 - 1) ;
}
if( b == 135 ) {
  val = -(18*xi_3 - 1)*(2*xi_1 + xi_2 + xi_3 - 1)*(780*xi_2*std::pow(xi_3 - 1, 6) + 365*std::pow(xi_3 - 1, 7) + 165*std::pow(xi_3 - 1, 5)*std::pow(2*xi_2 + xi_3 - 1, 2) - 3740*std::pow(xi_3 - 1, 4)*std::pow(2*xi_2 + xi_3 - 1, 3) + 2145*std::pow(xi_3 - 1, 3)*std::pow(2*xi_2 + xi_3 - 1, 4) + 6006*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 5) - 5005*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 6))/(64*xi_3 - 64) ;
}
if( b == 136 ) {
  val = (1.0L/32.0L)*(2*xi_1 + xi_2 + xi_3 - 1)*(154*xi_2*std::pow(xi_3 - 1, 7) - 924*std::pow(-xi_3 + 1, 6)*std::pow(-2*xi_2 - xi_3 + 1, 2) - 6006*std::pow(-xi_3 + 1, 2)*std::pow(-2*xi_2 - xi_3 + 1, 6) + 4862*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 7) + 98*std::pow(xi_3 - 1, 8) + 5005*std::pow(xi_3 - 1, 4)*std::pow(-2*xi_2 - xi_3 + 1, 4) + 3003*std::pow(xi_3 - 1, 3)*std::pow(-2*xi_2 - xi_3 + 1, 5))/(xi_3 - 1) ;
}
if( b == 137 ) {
  val = (1.0L/32.0L)*(-std::pow(xi_2 + xi_3 - 1, 2) + 3*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2))*(-156*xi_3 + 4641*std::pow(2*xi_3 - 1, 6) + 9282*std::pow(2*xi_3 - 1, 5) + 4095*std::pow(2*xi_3 - 1, 4) - 1820*std::pow(2*xi_3 - 1, 3) - 1365*std::pow(2*xi_3 - 1, 2) + 107) ;
}
if( b == 138 ) {
  val = (1.0L/8.0L)*(-std::pow(xi_2 + xi_3 - 1, 2) + 3*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2))*(7*xi_2 + xi_3 - 1)*(-210*xi_3 + 1071*std::pow(2*xi_3 - 1, 5) + 2380*std::pow(2*xi_3 - 1, 4) + 1610*std::pow(2*xi_3 - 1, 3) + 210*std::pow(2*xi_3 - 1, 2) + 87) ;
}
if( b == 139 ) {
  val = (std::pow(xi_2 + xi_3 - 1, 2) - 3*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2))*(20*xi_2*std::pow(xi_3 - 1, 2) + 8*std::pow(xi_3 - 1, 3) - 9*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2))*(600*xi_3 + 765*std::pow(2*xi_3 - 1, 4) + 1700*std::pow(2*xi_3 - 1, 3) + 1230*std::pow(2*xi_3 - 1, 2) - 291)/(8*xi_3 - 8) ;
}
if( b == 140 ) {
  val = (1.0L/16.0L)*(-std::pow(xi_2 + xi_3 - 1, 2) + 3*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2))*(252*xi_3 + 102*std::pow(2*xi_3 - 1, 3) + 204*std::pow(2*xi_3 - 1, 2) - 103)*(126*xi_2*std::pow(xi_3 - 1, 3) - 225*std::pow(-xi_3 + 1, 2)*std::pow(-2*xi_2 - xi_3 + 1, 2) + 165*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 3) + 68*std::pow(xi_3 - 1, 4))/(xi_3 - 1) ;
}
if( b == 141 ) {
  val = (1.0L/128.0L)*(-std::pow(xi_2 + xi_3 - 1, 2) + 3*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2))*(476*xi_3 + 153*std::pow(2*xi_3 - 1, 2) - 149)*(200*xi_2*std::pow(xi_3 - 1, 3) + 71*std::pow(xi_3 - 1, 4) + 330*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 2) - 1100*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 3) + 715*std::pow(2*xi_2 + xi_3 - 1, 4)) ;
}
if( b == 142 ) {
  val = (18*xi_3 - 1)*(std::pow(xi_2 + xi_3 - 1, 2) - 3*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2))*(770*xi_2*std::pow(xi_3 - 1, 5) + 386*std::pow(xi_3 - 1, 6) - 990*std::pow(xi_3 - 1, 4)*std::pow(2*xi_2 + xi_3 - 1, 2) - 1430*std::pow(xi_3 - 1, 3)*std::pow(2*xi_2 + xi_3 - 1, 3) + 5005*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 4) - 3003*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 5))/(64*xi_3 - 64) ;
}
if( b == 143 ) {
  val = (1.0L/32.0L)*(std::pow(xi_2 + xi_3 - 1, 2) - 3*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2))*(96*xi_2*std::pow(xi_3 - 1, 5) + 23*std::pow(xi_3 - 1, 6) + 780*std::pow(xi_3 - 1, 4)*std::pow(2*xi_2 + xi_3 - 1, 2) - 1820*std::pow(xi_3 - 1, 3)*std::pow(2*xi_2 + xi_3 - 1, 3) - 1365*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 4) + 5460*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 5) - 3094*std::pow(2*xi_2 + xi_3 - 1, 6)) ;
}
if( b == 144 ) {
  val = -(6*xi_1*std::pow(xi_2 + xi_3 - 1, 3) + 3*std::pow(xi_2 + xi_3 - 1, 4) - 5*(xi_2 + xi_3 - 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3))*(-210*xi_3 + 1071*std::pow(2*xi_3 - 1, 5) + 2380*std::pow(2*xi_3 - 1, 4) + 1610*std::pow(2*xi_3 - 1, 3) + 210*std::pow(2*xi_3 - 1, 2) + 87)/(8*xi_2 + 8*xi_3 - 8) ;
}
if( b == 145 ) {
  val = (1.0L/8.0L)*(9*xi_2 + xi_3 - 1)*(6*xi_1*std::pow(xi_2 + xi_3 - 1, 3) + 5*(-xi_2 - xi_3 + 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3) + 3*std::pow(xi_2 + xi_3 - 1, 4))*(600*xi_3 + 765*std::pow(2*xi_3 - 1, 4) + 1700*std::pow(2*xi_3 - 1, 3) + 1230*std::pow(2*xi_3 - 1, 2) - 291)/(-xi_2 - xi_3 + 1) ;
}
if( b == 146 ) {
  val = (1.0L/8.0L)*(6*xi_1*std::pow(xi_2 + xi_3 - 1, 3) + 5*(-xi_2 - xi_3 + 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3) + 3*std::pow(xi_2 + xi_3 - 1, 4))*(140*xi_2*std::pow(xi_3 - 1, 2) + 55*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 2) + 51*std::pow(xi_3 - 1, 3))*(252*xi_3 + 102*std::pow(2*xi_3 - 1, 3) + 204*std::pow(2*xi_3 - 1, 2) - 103)/((xi_3 - 1)*(xi_2 + xi_3 - 1)) ;
}
if( b == 147 ) {
  val = -(476*xi_3 + 153*std::pow(2*xi_3 - 1, 2) - 149)*(6*xi_1*std::pow(xi_2 + xi_3 - 1, 3) + 3*std::pow(xi_2 + xi_3 - 1, 4) - 5*(xi_2 + xi_3 - 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3))*(198*xi_2*std::pow(xi_3 - 1, 2) + 92*std::pow(xi_3 - 1, 3) - 231*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2) + 143*std::pow(2*xi_2 + xi_3 - 1, 3))/(32*xi_2 + 32*xi_3 - 32) ;
}
if( b == 148 ) {
  val = (1.0L/32.0L)*(18*xi_3 - 1)*(6*xi_1*std::pow(xi_2 + xi_3 - 1, 3) + 5*(-xi_2 - xi_3 + 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3) + 3*std::pow(xi_2 + xi_3 - 1, 4))*(-168*xi_2*std::pow(xi_3 - 1, 3) - 127*std::pow(xi_3 - 1, 4) + 1326*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 2) - 2548*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 3) + 1365*std::pow(2*xi_2 + xi_3 - 1, 4))/(-xi_2 - xi_3 + 1) ;
}
if( b == 149 ) {
  val = (1.0L/16.0L)*(6*xi_1*std::pow(xi_2 + xi_3 - 1, 3) + 5*(-xi_2 - xi_3 + 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3) + 3*std::pow(xi_2 + xi_3 - 1, 4))*(390*xi_2*std::pow(xi_3 - 1, 5) + 1820*std::pow(-xi_3 + 1, 3)*std::pow(2*xi_2 + xi_3 - 1, 3) + 1547*(-xi_3 + 1)*std::pow(2*xi_2 + xi_3 - 1, 5) + 174*std::pow(xi_3 - 1, 6) + 3185*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 4))/((xi_3 - 1)*(xi_2 + xi_3 - 1)) ;
}
if( b == 150 ) {
  val = (1.0L/32.0L)*(3*std::pow(xi_2 + xi_3 - 1, 4) - 30*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + 35*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 4))*(600*xi_3 + 765*std::pow(2*xi_3 - 1, 4) + 1700*std::pow(2*xi_3 - 1, 3) + 1230*std::pow(2*xi_3 - 1, 2) - 291) ;
}
if( b == 151 ) {
  val = (1.0L/8.0L)*(11*xi_2 + xi_3 - 1)*(3*std::pow(xi_2 + xi_3 - 1, 4) - 30*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + 35*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 4))*(252*xi_3 + 102*std::pow(2*xi_3 - 1, 3) + 204*std::pow(2*xi_3 - 1, 2) - 103) ;
}
if( b == 152 ) {
  val = -(476*xi_3 + 153*std::pow(2*xi_3 - 1, 2) - 149)*(108*xi_2*std::pow(xi_3 - 1, 2) + 37*std::pow(xi_3 - 1, 3) - 39*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2))*(3*std::pow(xi_2 + xi_3 - 1, 4) - 30*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + 35*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 4))/(64*xi_3 - 64) ;
}
if( b == 153 ) {
  val = (1.0L/64.0L)*(18*xi_3 - 1)*(3*std::pow(xi_2 + xi_3 - 1, 4) - 30*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + 35*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 4))*(858*xi_2*std::pow(xi_3 - 1, 3) - 819*std::pow(-xi_3 + 1, 2)*std::pow(-2*xi_2 - xi_3 + 1, 2) + 455*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 3) + 372*std::pow(xi_3 - 1, 4))/(xi_3 - 1) ;
}
if( b == 154 ) {
  val = (1.0L/32.0L)*(3*std::pow(xi_2 + xi_3 - 1, 4) - 30*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + 35*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 4))*(-336*xi_2*std::pow(xi_3 - 1, 3) - 171*std::pow(xi_3 - 1, 4) + 840*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 2) - 1260*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 3) + 595*std::pow(2*xi_2 + xi_3 - 1, 4)) ;
}
if( b == 155 ) {
  val = (1.0L/8.0L)*(252*xi_3 + 102*std::pow(2*xi_3 - 1, 3) + 204*std::pow(2*xi_3 - 1, 2) - 103)*(30*xi_1*std::pow(xi_2 + xi_3 - 1, 5) + 63*(-xi_2 - xi_3 + 1)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 5) + 15*std::pow(xi_2 + xi_3 - 1, 6) + 70*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 3))/(xi_2 + xi_3 - 1) ;
}
if( b == 156 ) {
  val = (1.0L/32.0L)*(13*xi_2 + xi_3 - 1)*(476*xi_3 + 153*std::pow(2*xi_3 - 1, 2) - 149)*(30*xi_1*std::pow(xi_2 + xi_3 - 1, 5) + 63*(-xi_2 - xi_3 + 1)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 5) + 15*std::pow(xi_2 + xi_3 - 1, 6) + 70*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 3))/(xi_2 + xi_3 - 1) ;
}
if( b == 157 ) {
  val = -1.0L/32.0L*(18*xi_3 - 1)*(308*xi_2*std::pow(xi_3 - 1, 2) + 101*std::pow(xi_3 - 1, 3) - 105*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2))*(30*xi_1*std::pow(xi_2 + xi_3 - 1, 5) + 15*std::pow(xi_2 + xi_3 - 1, 6) - 70*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3) + 63*(xi_2 + xi_3 - 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 5))/((xi_3 - 1)*(xi_2 + xi_3 - 1)) ;
}
if( b == 158 ) {
  val = (1.0L/8.0L)*(30*xi_1*std::pow(xi_2 + xi_3 - 1, 5) + 63*(-xi_2 - xi_3 + 1)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 5) + 15*std::pow(xi_2 + xi_3 - 1, 6) + 70*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 3))*(195*xi_2*std::pow(xi_3 - 1, 2) + 81*std::pow(xi_3 - 1, 3) - 165*(xi_3 - 1)*std::pow(-2*xi_2 - xi_3 + 1, 2) - 85*std::pow(-2*xi_2 - xi_3 + 1, 3))/(xi_2 + xi_3 - 1) ;
}
if( b == 159 ) {
  val = (1.0L/64.0L)*(476*xi_3 + 153*std::pow(2*xi_3 - 1, 2) - 149)*(-5*std::pow(xi_2 + xi_3 - 1, 6) + 105*std::pow(xi_2 + xi_3 - 1, 4)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) - 315*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 4) + 231*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 6)) ;
}
if( b == 160 ) {
  val = (1.0L/16.0L)*(18*xi_3 - 1)*(15*xi_2 + xi_3 - 1)*(-5*std::pow(xi_2 + xi_3 - 1, 6) + 105*std::pow(xi_2 + xi_3 - 1, 4)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) - 315*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 4) + 231*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 6)) ;
}
if( b == 161 ) {
  val = (104*xi_2*std::pow(xi_3 - 1, 2) + 33*std::pow(xi_3 - 1, 3) - 34*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2))*(5*std::pow(xi_2 + xi_3 - 1, 6) - 105*std::pow(xi_2 + xi_3 - 1, 4)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + 315*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 4) - 231*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 6))/(16*xi_3 - 16) ;
}
if( b == 162 ) {
  val = -(18*xi_3 - 1)*(70*xi_1*std::pow(xi_2 + xi_3 - 1, 7) + 35*std::pow(xi_2 + xi_3 - 1, 8) - 315*std::pow(xi_2 + xi_3 - 1, 5)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3) + 693*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 5) - 429*(xi_2 + xi_3 - 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 7))/(16*xi_2 + 16*xi_3 - 16) ;
}
if( b == 163 ) {
  val = -(17*xi_2 + xi_3 - 1)*(70*xi_1*std::pow(xi_2 + xi_3 - 1, 7) + 35*std::pow(xi_2 + xi_3 - 1, 8) - 315*std::pow(xi_2 + xi_3 - 1, 5)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3) + 693*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 5) - 429*(xi_2 + xi_3 - 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 7))/(16*xi_2 + 16*xi_3 - 16) ;
}
if( b == 164 ) {
  val = (35.0L/128.0L)*std::pow(xi_2 + xi_3 - 1, 8) - 315.0L/32.0L*std::pow(xi_2 + xi_3 - 1, 6)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + (3465.0L/64.0L)*std::pow(xi_2 + xi_3 - 1, 4)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 4) - 3003.0L/32.0L*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 6) + (6435.0L/128.0L)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 8) ;
}
if( b == 165 ) {
  val = 167960*std::pow(xi_3, 9) - 680238*std::pow(xi_3, 8) + 1145664*std::pow(xi_3, 7) - 1039584*std::pow(xi_3, 6) + 550368*std::pow(xi_3, 5) - 171990*std::pow(xi_3, 4) + 30576*std::pow(xi_3, 3) - 2808*std::pow(xi_3, 2) + 108*xi_3 - 1 ;
}
if( b == 166 ) {
  val = (1.0L/128.0L)*(3*xi_2 + xi_3 - 1)*(-1456*xi_3 + 62985*std::pow(2*xi_3 - 1, 8) + 100776*std::pow(2*xi_3 - 1, 7) - 18564*std::pow(2*xi_3 - 1, 6) - 86632*std::pow(2*xi_3 - 1, 5) - 17290*std::pow(2*xi_3 - 1, 4) + 18200*std::pow(2*xi_3 - 1, 3) + 4732*std::pow(2*xi_3 - 1, 2) + 609) ;
}
if( b == 167 ) {
  val = -(4*xi_2*std::pow(xi_3 - 1, 2) + 3*std::pow(xi_3 - 1, 3) - 5*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2))*(616*xi_3 + 9690*std::pow(2*xi_3 - 1, 7) + 20349*std::pow(2*xi_3 - 1, 6) + 8568*std::pow(2*xi_3 - 1, 5) - 6545*std::pow(2*xi_3 - 1, 4) - 4830*std::pow(2*xi_3 - 1, 3) - 105*std::pow(2*xi_3 - 1, 2) - 287)/(32*xi_3 - 32) ;
}
if( b == 168 ) {
  val = -(30*xi_2*std::pow(xi_3 - 1, 3) + 12*std::pow(xi_3 - 1, 4) + 15*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 2) - 35*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 3))*(-600*xi_3 + 4845*std::pow(2*xi_3 - 1, 6) + 11628*std::pow(2*xi_3 - 1, 5) + 8415*std::pow(2*xi_3 - 1, 4) + 680*std::pow(2*xi_3 - 1, 3) - 1245*std::pow(2*xi_3 - 1, 2) + 301)/(64*xi_3 - 64) ;
}
if( b == 169 ) {
  val = (1.0L/32.0L)*(24*xi_2*std::pow(xi_3 - 1, 3) + 15*std::pow(xi_3 - 1, 4) - 42*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 2) - 28*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 3) + 63*std::pow(2*xi_2 + xi_3 - 1, 4))*(-20*xi_3 + 1938*std::pow(2*xi_3 - 1, 5) + 4845*std::pow(2*xi_3 - 1, 4) + 4080*std::pow(2*xi_3 - 1, 3) + 1190*std::pow(2*xi_3 - 1, 2) - 21) ;
}
if( b == 170 ) {
  val = (1.0L/256.0L)*(5848*xi_3 + 4845*std::pow(2*xi_3 - 1, 4) + 11628*std::pow(2*xi_3 - 1, 3) + 9486*std::pow(2*xi_3 - 1, 2) - 2687)*(70*xi_2*std::pow(xi_3 - 1, 5) - 105*std::pow(-xi_3 + 1, 2)*std::pow(-2*xi_2 - xi_3 + 1, 4) + 231*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 5) + 30*std::pow(xi_3 - 1, 6) + 70*std::pow(xi_3 - 1, 4)*std::pow(-2*xi_2 - xi_3 + 1, 2) + 210*std::pow(xi_3 - 1, 3)*std::pow(-2*xi_2 - xi_3 + 1, 3))/(xi_3 - 1) ;
}
if( b == 171 ) {
  val = -(1584*xi_3 + 570*std::pow(2*xi_3 - 1, 3) + 1197*std::pow(2*xi_3 - 1, 2) - 631)*(60*xi_2*std::pow(xi_3 - 1, 6) + 35*std::pow(xi_3 - 1, 7) - 135*std::pow(xi_3 - 1, 5)*std::pow(2*xi_2 + xi_3 - 1, 2) - 180*std::pow(xi_3 - 1, 4)*std::pow(2*xi_2 + xi_3 - 1, 3) + 495*std::pow(xi_3 - 1, 3)*std::pow(2*xi_2 + xi_3 - 1, 4) + 198*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 5) - 429*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 6))/(64*xi_3 - 64) ;
}
if( b == 172 ) {
  val = -(304*xi_3 + 95*std::pow(2*xi_3 - 1, 2) - 93)*(630*xi_2*std::pow(xi_3 - 1, 7) + 280*std::pow(xi_3 - 1, 8) + 945*std::pow(xi_3 - 1, 6)*std::pow(2*xi_2 + xi_3 - 1, 2) - 3465*std::pow(xi_3 - 1, 5)*std::pow(2*xi_2 + xi_3 - 1, 3) - 3465*std::pow(xi_3 - 1, 4)*std::pow(2*xi_2 + xi_3 - 1, 4) + 9009*std::pow(xi_3 - 1, 3)*std::pow(2*xi_2 + xi_3 - 1, 5) + 3003*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 6) - 6435*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 7))/(256*xi_3 - 256) ;
}
if( b == 173 ) {
  val = (1.0L/128.0L)*(20*xi_3 - 1)*(560*xi_2*std::pow(xi_3 - 1, 7) - 1540*std::pow(-xi_3 + 1, 6)*std::pow(2*xi_2 + xi_3 - 1, 2) + 3080*std::pow(-xi_3 + 1, 5)*std::pow(2*xi_2 + xi_3 - 1, 3) - 20020*std::pow(-xi_3 + 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 6) + 5720*(-xi_3 + 1)*std::pow(2*xi_2 + xi_3 - 1, 7) + 315*std::pow(xi_3 - 1, 8) + 10010*std::pow(xi_3 - 1, 4)*std::pow(2*xi_2 + xi_3 - 1, 4) + 8008*std::pow(xi_3 - 1, 3)*std::pow(2*xi_2 + xi_3 - 1, 5) + 12155*std::pow(2*xi_2 + xi_3 - 1, 8)) ;
}
if( b == 174 ) {
  val = (1.0L/256.0L)*(1386*xi_2*std::pow(xi_3 - 1, 9) - 18018*std::pow(-xi_3 + 1, 6)*std::pow(-2*xi_2 - xi_3 + 1, 4) + 54054*std::pow(-xi_3 + 1, 5)*std::pow(-2*xi_2 - xi_3 + 1, 5) - 21879*std::pow(-xi_3 + 1, 2)*std::pow(-2*xi_2 - xi_3 + 1, 8) + 46189*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 9) + 630*std::pow(xi_3 - 1, 10) + 2772*std::pow(xi_3 - 1, 8)*std::pow(-2*xi_2 - xi_3 + 1, 2) + 12012*std::pow(xi_3 - 1, 7)*std::pow(-2*xi_2 - xi_3 + 1, 3) + 36036*std::pow(xi_3 - 1, 4)*std::pow(-2*xi_2 - xi_3 + 1, 6) + 87516*std::pow(xi_3 - 1, 3)*std::pow(-2*xi_2 - xi_3 + 1, 7))/(xi_3 - 1) ;
}
if( b == 175 ) {
  val = (1.0L/128.0L)*(2*xi_1 + xi_2 + xi_3 - 1)*(-1456*xi_3 + 62985*std::pow(2*xi_3 - 1, 8) + 100776*std::pow(2*xi_3 - 1, 7) - 18564*std::pow(2*xi_3 - 1, 6) - 86632*std::pow(2*xi_3 - 1, 5) - 17290*std::pow(2*xi_3 - 1, 4) + 18200*std::pow(2*xi_3 - 1, 3) + 4732*std::pow(2*xi_3 - 1, 2) + 609) ;
}
if( b == 176 ) {
  val = (1.0L/16.0L)*(5*xi_2 + xi_3 - 1)*(2*xi_1 + xi_2 + xi_3 - 1)*(616*xi_3 + 9690*std::pow(2*xi_3 - 1, 7) + 20349*std::pow(2*xi_3 - 1, 6) + 8568*std::pow(2*xi_3 - 1, 5) - 6545*std::pow(2*xi_3 - 1, 4) - 4830*std::pow(2*xi_3 - 1, 3) - 105*std::pow(2*xi_3 - 1, 2) - 287) ;
}
if( b == 177 ) {
  val = -(36*xi_2*std::pow(xi_3 - 1, 2) + 17*std::pow(xi_3 - 1, 3) - 21*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2))*(2*xi_1 + xi_2 + xi_3 - 1)*(-600*xi_3 + 4845*std::pow(2*xi_3 - 1, 6) + 11628*std::pow(2*xi_3 - 1, 5) + 8415*std::pow(2*xi_3 - 1, 4) + 680*std::pow(2*xi_3 - 1, 3) - 1245*std::pow(2*xi_3 - 1, 2) + 301)/(32*xi_3 - 32) ;
}
if( b == 178 ) {
  val = (1.0L/8.0L)*(2*std::pow(xi_3 - 1, 3) - 21*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2) + 21*std::pow(2*xi_2 + xi_3 - 1, 3))*(2*xi_1 + xi_2 + xi_3 - 1)*(-20*xi_3 + 1938*std::pow(2*xi_3 - 1, 5) + 4845*std::pow(2*xi_3 - 1, 4) + 4080*std::pow(2*xi_3 - 1, 3) + 1190*std::pow(2*xi_3 - 1, 2) - 21) ;
}
if( b == 179 ) {
  val = (1.0L/128.0L)*(2*xi_1 + xi_2 + xi_3 - 1)*(5848*xi_3 + 4845*std::pow(2*xi_3 - 1, 4) + 11628*std::pow(2*xi_3 - 1, 3) + 9486*std::pow(2*xi_3 - 1, 2) - 2687)*(88*xi_2*std::pow(xi_3 - 1, 3) + 41*std::pow(xi_3 - 1, 4) - 18*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 2) - 180*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 3) + 165*std::pow(2*xi_2 + xi_3 - 1, 4)) ;
}
if( b == 180 ) {
  val = (1.0L/128.0L)*(2*xi_1 + xi_2 + xi_3 - 1)*(1584*xi_3 + 570*std::pow(2*xi_3 - 1, 3) + 1197*std::pow(2*xi_3 - 1, 2) - 631)*(-90*xi_2*std::pow(xi_3 - 1, 5) - 70*std::pow(-xi_3 + 1, 6) - 1485*std::pow(-xi_3 + 1, 2)*std::pow(-2*xi_2 - xi_3 + 1, 4) + 1287*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 5) + 630*std::pow(xi_3 - 1, 4)*std::pow(-2*xi_2 - xi_3 + 1, 2) + 330*std::pow(xi_3 - 1, 3)*std::pow(-2*xi_2 - xi_3 + 1, 3))/(xi_3 - 1) ;
}
if( b == 181 ) {
  val = -(304*xi_3 + 95*std::pow(2*xi_3 - 1, 2) - 93)*(2*xi_1 + xi_2 + xi_3 - 1)*(780*xi_2*std::pow(xi_3 - 1, 6) + 365*std::pow(xi_3 - 1, 7) + 165*std::pow(xi_3 - 1, 5)*std::pow(2*xi_2 + xi_3 - 1, 2) - 3740*std::pow(xi_3 - 1, 4)*std::pow(2*xi_2 + xi_3 - 1, 3) + 2145*std::pow(xi_3 - 1, 3)*std::pow(2*xi_2 + xi_3 - 1, 4) + 6006*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 5) - 5005*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 6))/(128*xi_3 - 128) ;
}
if( b == 182 ) {
  val = (1.0L/32.0L)*(20*xi_3 - 1)*(2*xi_1 + xi_2 + xi_3 - 1)*(154*xi_2*std::pow(xi_3 - 1, 7) - 924*std::pow(-xi_3 + 1, 6)*std::pow(-2*xi_2 - xi_3 + 1, 2) - 6006*std::pow(-xi_3 + 1, 2)*std::pow(-2*xi_2 - xi_3 + 1, 6) + 4862*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 7) + 98*std::pow(xi_3 - 1, 8) + 5005*std::pow(xi_3 - 1, 4)*std::pow(-2*xi_2 - xi_3 + 1, 4) + 3003*std::pow(xi_3 - 1, 3)*std::pow(-2*xi_2 - xi_3 + 1, 5))/(xi_3 - 1) ;
}
if( b == 183 ) {
  val = (1.0L/128.0L)*(2*xi_1 + xi_2 + xi_3 - 1)*(1680*xi_2*std::pow(xi_3 - 1, 7) + 13832*std::pow(-xi_3 + 1, 5)*std::pow(2*xi_2 + xi_3 - 1, 3) - 30940*std::pow(-xi_3 + 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 6) + 47736*(-xi_3 + 1)*std::pow(2*xi_2 + xi_3 - 1, 7) + 791*std::pow(xi_3 - 1, 8) + 1092*std::pow(xi_3 - 1, 6)*std::pow(2*xi_2 + xi_3 - 1, 2) + 2730*std::pow(xi_3 - 1, 4)*std::pow(2*xi_2 + xi_3 - 1, 4) + 50232*std::pow(xi_3 - 1, 3)*std::pow(2*xi_2 + xi_3 - 1, 5) + 37791*std::pow(2*xi_2 + xi_3 - 1, 8)) ;
}
if( b == 184 ) {
  val = (1.0L/32.0L)*(-std::pow(xi_2 + xi_3 - 1, 2) + 3*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2))*(616*xi_3 + 9690*std::pow(2*xi_3 - 1, 7) + 20349*std::pow(2*xi_3 - 1, 6) + 8568*std::pow(2*xi_3 - 1, 5) - 6545*std::pow(2*xi_3 - 1, 4) - 4830*std::pow(2*xi_3 - 1, 3) - 105*std::pow(2*xi_3 - 1, 2) - 287) ;
}
if( b == 185 ) {
  val = (1.0L/16.0L)*(-std::pow(xi_2 + xi_3 - 1, 2) + 3*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2))*(7*xi_2 + xi_3 - 1)*(-600*xi_3 + 4845*std::pow(2*xi_3 - 1, 6) + 11628*std::pow(2*xi_3 - 1, 5) + 8415*std::pow(2*xi_3 - 1, 4) + 680*std::pow(2*xi_3 - 1, 3) - 1245*std::pow(2*xi_3 - 1, 2) + 301) ;
}
if( b == 186 ) {
  val = (std::pow(xi_2 + xi_3 - 1, 2) - 3*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2))*(20*xi_2*std::pow(xi_3 - 1, 2) + 8*std::pow(xi_3 - 1, 3) - 9*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2))*(-20*xi_3 + 1938*std::pow(2*xi_3 - 1, 5) + 4845*std::pow(2*xi_3 - 1, 4) + 4080*std::pow(2*xi_3 - 1, 3) + 1190*std::pow(2*xi_3 - 1, 2) - 21)/(8*xi_3 - 8) ;
}
if( b == 187 ) {
  val = (1.0L/256.0L)*(-std::pow(xi_2 + xi_3 - 1, 2) + 3*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2))*(126*xi_2*std::pow(xi_3 - 1, 3) - 225*std::pow(-xi_3 + 1, 2)*std::pow(-2*xi_2 - xi_3 + 1, 2) + 165*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 3) + 68*std::pow(xi_3 - 1, 4))*(5848*xi_3 + 4845*std::pow(2*xi_3 - 1, 4) + 11628*std::pow(2*xi_3 - 1, 3) + 9486*std::pow(2*xi_3 - 1, 2) - 2687)/(xi_3 - 1) ;
}
if( b == 188 ) {
  val = (1.0L/128.0L)*(-std::pow(xi_2 + xi_3 - 1, 2) + 3*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2))*(1584*xi_3 + 570*std::pow(2*xi_3 - 1, 3) + 1197*std::pow(2*xi_3 - 1, 2) - 631)*(200*xi_2*std::pow(xi_3 - 1, 3) + 71*std::pow(xi_3 - 1, 4) + 330*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 2) - 1100*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 3) + 715*std::pow(2*xi_2 + xi_3 - 1, 4)) ;
}
if( b == 189 ) {
  val = (std::pow(xi_2 + xi_3 - 1, 2) - 3*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2))*(304*xi_3 + 95*std::pow(2*xi_3 - 1, 2) - 93)*(770*xi_2*std::pow(xi_3 - 1, 5) + 386*std::pow(xi_3 - 1, 6) - 990*std::pow(xi_3 - 1, 4)*std::pow(2*xi_2 + xi_3 - 1, 2) - 1430*std::pow(xi_3 - 1, 3)*std::pow(2*xi_2 + xi_3 - 1, 3) + 5005*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 4) - 3003*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 5))/(128*xi_3 - 128) ;
}
if( b == 190 ) {
  val = (1.0L/32.0L)*(20*xi_3 - 1)*(std::pow(xi_2 + xi_3 - 1, 2) - 3*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2))*(96*xi_2*std::pow(xi_3 - 1, 5) + 23*std::pow(xi_3 - 1, 6) + 780*std::pow(xi_3 - 1, 4)*std::pow(2*xi_2 + xi_3 - 1, 2) - 1820*std::pow(xi_3 - 1, 3)*std::pow(2*xi_2 + xi_3 - 1, 3) - 1365*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 4) + 5460*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 5) - 3094*std::pow(2*xi_2 + xi_3 - 1, 6)) ;
}
if( b == 191 ) {
  val = (1.0L/64.0L)*(-std::pow(xi_2 + xi_3 - 1, 2) + 3*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2))*(910*xi_2*std::pow(xi_3 - 1, 7) - 819*std::pow(-xi_3 + 1, 6)*std::pow(-2*xi_2 - xi_3 + 1, 2) + 4641*std::pow(-xi_3 + 1, 3)*std::pow(-2*xi_2 - xi_3 + 1, 5) - 23205*std::pow(-xi_3 + 1, 2)*std::pow(-2*xi_2 - xi_3 + 1, 6) + 12597*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 7) + 448*std::pow(xi_3 - 1, 8) + 5005*std::pow(xi_3 - 1, 5)*std::pow(-2*xi_2 - xi_3 + 1, 3) + 11375*std::pow(xi_3 - 1, 4)*std::pow(-2*xi_2 - xi_3 + 1, 4))/(xi_3 - 1) ;
}
if( b == 192 ) {
  val = -(6*xi_1*std::pow(xi_2 + xi_3 - 1, 3) + 3*std::pow(xi_2 + xi_3 - 1, 4) - 5*(xi_2 + xi_3 - 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3))*(-600*xi_3 + 4845*std::pow(2*xi_3 - 1, 6) + 11628*std::pow(2*xi_3 - 1, 5) + 8415*std::pow(2*xi_3 - 1, 4) + 680*std::pow(2*xi_3 - 1, 3) - 1245*std::pow(2*xi_3 - 1, 2) + 301)/(16*xi_2 + 16*xi_3 - 16) ;
}
if( b == 193 ) {
  val = (1.0L/8.0L)*(9*xi_2 + xi_3 - 1)*(6*xi_1*std::pow(xi_2 + xi_3 - 1, 3) + 5*(-xi_2 - xi_3 + 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3) + 3*std::pow(xi_2 + xi_3 - 1, 4))*(-20*xi_3 + 1938*std::pow(2*xi_3 - 1, 5) + 4845*std::pow(2*xi_3 - 1, 4) + 4080*std::pow(2*xi_3 - 1, 3) + 1190*std::pow(2*xi_3 - 1, 2) - 21)/(-xi_2 - xi_3 + 1) ;
}
if( b == 194 ) {
  val = (1.0L/128.0L)*(6*xi_1*std::pow(xi_2 + xi_3 - 1, 3) + 5*(-xi_2 - xi_3 + 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3) + 3*std::pow(xi_2 + xi_3 - 1, 4))*(140*xi_2*std::pow(xi_3 - 1, 2) + 55*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 2) + 51*std::pow(xi_3 - 1, 3))*(5848*xi_3 + 4845*std::pow(2*xi_3 - 1, 4) + 11628*std::pow(2*xi_3 - 1, 3) + 9486*std::pow(2*xi_3 - 1, 2) - 2687)/((xi_3 - 1)*(xi_2 + xi_3 - 1)) ;
}
if( b == 195 ) {
  val = -(6*xi_1*std::pow(xi_2 + xi_3 - 1, 3) + 3*std::pow(xi_2 + xi_3 - 1, 4) - 5*(xi_2 + xi_3 - 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3))*(1584*xi_3 + 570*std::pow(2*xi_3 - 1, 3) + 1197*std::pow(2*xi_3 - 1, 2) - 631)*(198*xi_2*std::pow(xi_3 - 1, 2) + 92*std::pow(xi_3 - 1, 3) - 231*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2) + 143*std::pow(2*xi_2 + xi_3 - 1, 3))/(32*xi_2 + 32*xi_3 - 32) ;
}
if( b == 196 ) {
  val = (304*xi_3 + 95*std::pow(2*xi_3 - 1, 2) - 93)*(6*xi_1*std::pow(xi_2 + xi_3 - 1, 3) + 3*std::pow(xi_2 + xi_3 - 1, 4) - 5*(xi_2 + xi_3 - 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3))*(168*xi_2*std::pow(xi_3 - 1, 3) + 127*std::pow(xi_3 - 1, 4) - 1326*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 2) + 2548*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 3) - 1365*std::pow(2*xi_2 + xi_3 - 1, 4))/(64*xi_2 + 64*xi_3 - 64) ;
}
if( b == 197 ) {
  val = (1.0L/16.0L)*(20*xi_3 - 1)*(6*xi_1*std::pow(xi_2 + xi_3 - 1, 3) + 5*(-xi_2 - xi_3 + 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3) + 3*std::pow(xi_2 + xi_3 - 1, 4))*(390*xi_2*std::pow(xi_3 - 1, 5) + 1820*std::pow(-xi_3 + 1, 3)*std::pow(2*xi_2 + xi_3 - 1, 3) + 1547*(-xi_3 + 1)*std::pow(2*xi_2 + xi_3 - 1, 5) + 174*std::pow(xi_3 - 1, 6) + 3185*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 4))/((xi_3 - 1)*(xi_2 + xi_3 - 1)) ;
}
if( b == 198 ) {
  val = -(6*xi_1*std::pow(xi_2 + xi_3 - 1, 3) + 3*std::pow(xi_2 + xi_3 - 1, 4) - 5*(xi_2 + xi_3 - 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3))*(588*xi_2*std::pow(xi_3 - 1, 5) + 317*std::pow(xi_3 - 1, 6) - 1995*std::pow(xi_3 - 1, 4)*std::pow(2*xi_2 + xi_3 - 1, 2) + 980*std::pow(xi_3 - 1, 3)*std::pow(2*xi_2 + xi_3 - 1, 3) + 8925*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 4) - 14994*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 5) + 6783*std::pow(2*xi_2 + xi_3 - 1, 6))/(32*xi_2 + 32*xi_3 - 32) ;
}
if( b == 199 ) {
  val = (1.0L/32.0L)*(3*std::pow(xi_2 + xi_3 - 1, 4) - 30*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + 35*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 4))*(-20*xi_3 + 1938*std::pow(2*xi_3 - 1, 5) + 4845*std::pow(2*xi_3 - 1, 4) + 4080*std::pow(2*xi_3 - 1, 3) + 1190*std::pow(2*xi_3 - 1, 2) - 21) ;
}
if( b == 200 ) {
  val = (1.0L/128.0L)*(11*xi_2 + xi_3 - 1)*(3*std::pow(xi_2 + xi_3 - 1, 4) - 30*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + 35*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 4))*(5848*xi_3 + 4845*std::pow(2*xi_3 - 1, 4) + 11628*std::pow(2*xi_3 - 1, 3) + 9486*std::pow(2*xi_3 - 1, 2) - 2687) ;
}
if( b == 201 ) {
  val = -(108*xi_2*std::pow(xi_3 - 1, 2) + 37*std::pow(xi_3 - 1, 3) - 39*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2))*(3*std::pow(xi_2 + xi_3 - 1, 4) - 30*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + 35*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 4))*(1584*xi_3 + 570*std::pow(2*xi_3 - 1, 3) + 1197*std::pow(2*xi_3 - 1, 2) - 631)/(64*xi_3 - 64) ;
}
if( b == 202 ) {
  val = (1.0L/128.0L)*(304*xi_3 + 95*std::pow(2*xi_3 - 1, 2) - 93)*(3*std::pow(xi_2 + xi_3 - 1, 4) - 30*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + 35*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 4))*(858*xi_2*std::pow(xi_3 - 1, 3) - 819*std::pow(-xi_3 + 1, 2)*std::pow(-2*xi_2 - xi_3 + 1, 2) + 455*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 3) + 372*std::pow(xi_3 - 1, 4))/(xi_3 - 1) ;
}
if( b == 203 ) {
  val = (1.0L/32.0L)*(20*xi_3 - 1)*(3*std::pow(xi_2 + xi_3 - 1, 4) - 30*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + 35*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 4))*(-336*xi_2*std::pow(xi_3 - 1, 3) - 171*std::pow(xi_3 - 1, 4) + 840*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 2) - 1260*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 3) + 595*std::pow(2*xi_2 + xi_3 - 1, 4)) ;
}
if( b == 204 ) {
  val = -(3*std::pow(xi_2 + xi_3 - 1, 4) - 30*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + 35*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 4))*(330*xi_2*std::pow(xi_3 - 1, 5) + 114*std::pow(xi_3 - 1, 6) + 1170*std::pow(xi_3 - 1, 4)*std::pow(2*xi_2 + xi_3 - 1, 2) - 5270*std::pow(xi_3 - 1, 3)*std::pow(2*xi_2 + xi_3 - 1, 3) + 6885*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 4) - 2907*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 5))/(64*xi_3 - 64) ;
}
if( b == 205 ) {
  val = (1.0L/128.0L)*(30*xi_1*std::pow(xi_2 + xi_3 - 1, 5) + 63*(-xi_2 - xi_3 + 1)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 5) + 15*std::pow(xi_2 + xi_3 - 1, 6) + 70*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 3))*(5848*xi_3 + 4845*std::pow(2*xi_3 - 1, 4) + 11628*std::pow(2*xi_3 - 1, 3) + 9486*std::pow(2*xi_3 - 1, 2) - 2687)/(xi_2 + xi_3 - 1) ;
}
if( b == 206 ) {
  val = (1.0L/32.0L)*(13*xi_2 + xi_3 - 1)*(1584*xi_3 + 570*std::pow(2*xi_3 - 1, 3) + 1197*std::pow(2*xi_3 - 1, 2) - 631)*(30*xi_1*std::pow(xi_2 + xi_3 - 1, 5) + 63*(-xi_2 - xi_3 + 1)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 5) + 15*std::pow(xi_2 + xi_3 - 1, 6) + 70*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 3))/(xi_2 + xi_3 - 1) ;
}
if( b == 207 ) {
  val = -1.0L/64.0L*(304*xi_3 + 95*std::pow(2*xi_3 - 1, 2) - 93)*(308*xi_2*std::pow(xi_3 - 1, 2) + 101*std::pow(xi_3 - 1, 3) - 105*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2))*(30*xi_1*std::pow(xi_2 + xi_3 - 1, 5) + 15*std::pow(xi_2 + xi_3 - 1, 6) - 70*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3) + 63*(xi_2 + xi_3 - 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 5))/((xi_3 - 1)*(xi_2 + xi_3 - 1)) ;
}
if( b == 208 ) {
  val = (1.0L/8.0L)*(20*xi_3 - 1)*(30*xi_1*std::pow(xi_2 + xi_3 - 1, 5) + 63*(-xi_2 - xi_3 + 1)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 5) + 15*std::pow(xi_2 + xi_3 - 1, 6) + 70*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 3))*(195*xi_2*std::pow(xi_3 - 1, 2) + 81*std::pow(xi_3 - 1, 3) - 165*(xi_3 - 1)*std::pow(-2*xi_2 - xi_3 + 1, 2) - 85*std::pow(-2*xi_2 - xi_3 + 1, 3))/(xi_2 + xi_3 - 1) ;
}
if( b == 209 ) {
  val = (1.0L/32.0L)*(30*xi_1*std::pow(xi_2 + xi_3 - 1, 5) + 63*(-xi_2 - xi_3 + 1)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 5) + 15*std::pow(xi_2 + xi_3 - 1, 6) + 70*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 3))*(-968*xi_2*std::pow(xi_3 - 1, 3) - 455*std::pow(xi_3 - 1, 4) + 1734*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 2) - 2244*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 3) + 969*std::pow(2*xi_2 + xi_3 - 1, 4))/(xi_2 + xi_3 - 1) ;
}
if( b == 210 ) {
  val = (1.0L/64.0L)*(1584*xi_3 + 570*std::pow(2*xi_3 - 1, 3) + 1197*std::pow(2*xi_3 - 1, 2) - 631)*(-5*std::pow(xi_2 + xi_3 - 1, 6) + 105*std::pow(xi_2 + xi_3 - 1, 4)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) - 315*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 4) + 231*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 6)) ;
}
if( b == 211 ) {
  val = (1.0L/32.0L)*(15*xi_2 + xi_3 - 1)*(304*xi_3 + 95*std::pow(2*xi_3 - 1, 2) - 93)*(-5*std::pow(xi_2 + xi_3 - 1, 6) + 105*std::pow(xi_2 + xi_3 - 1, 4)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) - 315*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 4) + 231*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 6)) ;
}
if( b == 212 ) {
  val = (20*xi_3 - 1)*(104*xi_2*std::pow(xi_3 - 1, 2) + 33*std::pow(xi_3 - 1, 3) - 34*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2))*(5*std::pow(xi_2 + xi_3 - 1, 6) - 105*std::pow(xi_2 + xi_3 - 1, 4)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + 315*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 4) - 231*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 6))/(16*xi_3 - 16) ;
}
if( b == 213 ) {
  val = (1.0L/128.0L)*(2550*xi_2*std::pow(xi_3 - 1, 3) - 1989*std::pow(-xi_3 + 1, 2)*std::pow(-2*xi_2 - xi_3 + 1, 2) + 969*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 3) + 1028*std::pow(xi_3 - 1, 4))*(-5*std::pow(xi_2 + xi_3 - 1, 6) + 105*std::pow(xi_2 + xi_3 - 1, 4)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) - 315*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 4) + 231*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 6))/(xi_3 - 1) ;
}
if( b == 214 ) {
  val = -(304*xi_3 + 95*std::pow(2*xi_3 - 1, 2) - 93)*(70*xi_1*std::pow(xi_2 + xi_3 - 1, 7) + 35*std::pow(xi_2 + xi_3 - 1, 8) - 315*std::pow(xi_2 + xi_3 - 1, 5)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3) + 693*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 5) - 429*(xi_2 + xi_3 - 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 7))/(32*xi_2 + 32*xi_3 - 32) ;
}
if( b == 215 ) {
  val = -(20*xi_3 - 1)*(17*xi_2 + xi_3 - 1)*(70*xi_1*std::pow(xi_2 + xi_3 - 1, 7) + 35*std::pow(xi_2 + xi_3 - 1, 8) - 315*std::pow(xi_2 + xi_3 - 1, 5)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3) + 693*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 5) - 429*(xi_2 + xi_3 - 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 7))/(16*xi_2 + 16*xi_3 - 16) ;
}
if( b == 216 ) {
  val = (1.0L/64.0L)*(540*xi_2*std::pow(xi_3 - 1, 2) + 171*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 2) + 167*std::pow(xi_3 - 1, 3))*(70*xi_1*std::pow(xi_2 + xi_3 - 1, 7) + 315*std::pow(-xi_2 - xi_3 + 1, 5)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3) + 429*(-xi_2 - xi_3 + 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 7) + 35*std::pow(xi_2 + xi_3 - 1, 8) + 693*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 5))/((xi_3 - 1)*(xi_2 + xi_3 - 1)) ;
}
if( b == 217 ) {
  val = (1.0L/128.0L)*(20*xi_3 - 1)*(35*std::pow(xi_2 + xi_3 - 1, 8) - 1260*std::pow(xi_2 + xi_3 - 1, 6)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + 6930*std::pow(xi_2 + xi_3 - 1, 4)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 4) - 12012*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 6) + 6435*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 8)) ;
}
if( b == 218 ) {
  val = (1.0L/128.0L)*(19*xi_2 + xi_3 - 1)*(35*std::pow(xi_2 + xi_3 - 1, 8) - 1260*std::pow(xi_2 + xi_3 - 1, 6)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + 6930*std::pow(xi_2 + xi_3 - 1, 4)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 4) - 12012*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 6) + 6435*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 8)) ;
}
if( b == 219 ) {
  val = (1.0L/128.0L)*(630*xi_1*std::pow(xi_2 + xi_3 - 1, 9) + 18018*std::pow(-xi_2 - xi_3 + 1, 5)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 5) + 12155*(-xi_2 - xi_3 + 1)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 9) + 315*std::pow(xi_2 + xi_3 - 1, 10) + 4620*std::pow(xi_2 + xi_3 - 1, 7)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 3) + 25740*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 7))/(xi_2 + xi_3 - 1) ;
}
if( b == 220 ) {
  val = 646646*std::pow(xi_3, 10) - 2939300*std::pow(xi_3, 9) + 5668650*std::pow(xi_3, 8) - 6046560*std::pow(xi_3, 7) + 3898440*std::pow(xi_3, 6) - 1559376*std::pow(xi_3, 5) + 382200*std::pow(xi_3, 4) - 54600*std::pow(xi_3, 3) + 4095*std::pow(xi_3, 2) - 130*xi_3 + 1 ;
}
if( b == 221 ) {
  val = -1.0L/128.0L*(3*xi_2 + xi_3 - 1)*(2394*xi_3 - 124355*std::pow(2*xi_3 - 1, 9) - 203490*std::pow(2*xi_3 - 1, 8) + 58140*std::pow(2*xi_3 - 1, 7) + 217056*std::pow(2*xi_3 - 1, 6) + 36414*std::pow(2*xi_3 - 1, 5) - 64260*std::pow(2*xi_3 - 1, 4) - 17220*std::pow(2*xi_3 - 1, 3) + 5040*std::pow(2*xi_3 - 1, 2) - 1239) ;
}
if( b == 222 ) {
  val = -(4*xi_2*std::pow(xi_3 - 1, 2) + 3*std::pow(xi_3 - 1, 3) - 5*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2))*(2640*xi_3 + 159885*std::pow(2*xi_3 - 1, 8) + 348840*std::pow(2*xi_3 - 1, 7) + 135660*std::pow(2*xi_3 - 1, 6) - 162792*std::pow(2*xi_3 - 1, 5) - 117810*std::pow(2*xi_3 - 1, 4) + 4760*std::pow(2*xi_3 - 1, 3) + 14700*std::pow(2*xi_3 - 1, 2) - 1499)/(256*xi_3 - 256) ;
}
if( b == 223 ) {
  val = -(30*xi_2*std::pow(xi_3 - 1, 3) + 12*std::pow(xi_3 - 1, 4) + 15*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 2) - 35*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 3))*(182*xi_3 + 10659*std::pow(2*xi_3 - 1, 7) + 27132*std::pow(2*xi_3 - 1, 6) + 20349*std::pow(2*xi_3 - 1, 5) - 5355*std::pow(2*xi_3 - 1, 3) - 1428*std::pow(2*xi_3 - 1, 2) - 59)/(64*xi_3 - 64) ;
}
if( b == 224 ) {
  val = (1.0L/512.0L)*(24*xi_2*std::pow(xi_3 - 1, 3) + 15*std::pow(xi_3 - 1, 4) - 42*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 2) - 28*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 3) + 63*std::pow(2*xi_2 + xi_3 - 1, 4))*(-10812*xi_3 + 74613*std::pow(2*xi_3 - 1, 6) + 203490*std::pow(2*xi_3 - 1, 5) + 188955*std::pow(2*xi_3 - 1, 4) + 58140*std::pow(2*xi_3 - 1, 3) - 6885*std::pow(2*xi_3 - 1, 2) + 5011) ;
}
if( b == 225 ) {
  val = (1.0L/256.0L)*(3150*xi_3 + 13167*std::pow(2*xi_3 - 1, 5) + 35910*std::pow(2*xi_3 - 1, 4) + 34770*std::pow(2*xi_3 - 1, 3) + 13680*std::pow(2*xi_3 - 1, 2) - 1669)*(70*xi_2*std::pow(xi_3 - 1, 5) - 105*std::pow(-xi_3 + 1, 2)*std::pow(-2*xi_2 - xi_3 + 1, 4) + 231*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 5) + 30*std::pow(xi_3 - 1, 6) + 70*std::pow(xi_3 - 1, 4)*std::pow(-2*xi_2 - xi_3 + 1, 2) + 210*std::pow(xi_3 - 1, 3)*std::pow(-2*xi_2 - xi_3 + 1, 3))/(xi_3 - 1) ;
}
if( b == 226 ) {
  val = -(11704*xi_3 + 7315*std::pow(2*xi_3 - 1, 4) + 18620*std::pow(2*xi_3 - 1, 3) + 16530*std::pow(2*xi_3 - 1, 2) - 5209)*(60*xi_2*std::pow(xi_3 - 1, 6) + 35*std::pow(xi_3 - 1, 7) - 135*std::pow(xi_3 - 1, 5)*std::pow(2*xi_2 + xi_3 - 1, 2) - 180*std::pow(xi_3 - 1, 4)*std::pow(2*xi_2 + xi_3 - 1, 3) + 495*std::pow(xi_3 - 1, 3)*std::pow(2*xi_2 + xi_3 - 1, 4) + 198*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 5) - 429*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 6))/(256*xi_3 - 256) ;
}
if( b == 227 ) {
  val = -(1170*xi_3 + 385*std::pow(2*xi_3 - 1, 3) + 840*std::pow(2*xi_3 - 1, 2) - 457)*(630*xi_2*std::pow(xi_3 - 1, 7) + 280*std::pow(xi_3 - 1, 8) + 945*std::pow(xi_3 - 1, 6)*std::pow(2*xi_2 + xi_3 - 1, 2) - 3465*std::pow(xi_3 - 1, 5)*std::pow(2*xi_2 + xi_3 - 1, 3) - 3465*std::pow(xi_3 - 1, 4)*std::pow(2*xi_2 + xi_3 - 1, 4) + 9009*std::pow(xi_3 - 1, 3)*std::pow(2*xi_2 + xi_3 - 1, 5) + 3003*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 6) - 6435*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 7))/(256*xi_3 - 256) ;
}
if( b == 228 ) {
  val = (1.0L/512.0L)*(756*xi_3 + 231*std::pow(2*xi_3 - 1, 2) - 227)*(560*xi_2*std::pow(xi_3 - 1, 7) - 1540*std::pow(-xi_3 + 1, 6)*std::pow(2*xi_2 + xi_3 - 1, 2) + 3080*std::pow(-xi_3 + 1, 5)*std::pow(2*xi_2 + xi_3 - 1, 3) - 20020*std::pow(-xi_3 + 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 6) + 5720*(-xi_3 + 1)*std::pow(2*xi_2 + xi_3 - 1, 7) + 315*std::pow(xi_3 - 1, 8) + 10010*std::pow(xi_3 - 1, 4)*std::pow(2*xi_2 + xi_3 - 1, 4) + 8008*std::pow(xi_3 - 1, 3)*std::pow(2*xi_2 + xi_3 - 1, 5) + 12155*std::pow(2*xi_2 + xi_3 - 1, 8)) ;
}
if( b == 229 ) {
  val = (1.0L/256.0L)*(22*xi_3 - 1)*(1386*xi_2*std::pow(xi_3 - 1, 9) - 18018*std::pow(-xi_3 + 1, 6)*std::pow(-2*xi_2 - xi_3 + 1, 4) + 54054*std::pow(-xi_3 + 1, 5)*std::pow(-2*xi_2 - xi_3 + 1, 5) - 21879*std::pow(-xi_3 + 1, 2)*std::pow(-2*xi_2 - xi_3 + 1, 8) + 46189*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 9) + 630*std::pow(xi_3 - 1, 10) + 2772*std::pow(xi_3 - 1, 8)*std::pow(-2*xi_2 - xi_3 + 1, 2) + 12012*std::pow(xi_3 - 1, 7)*std::pow(-2*xi_2 - xi_3 + 1, 3) + 36036*std::pow(xi_3 - 1, 4)*std::pow(-2*xi_2 - xi_3 + 1, 6) + 87516*std::pow(xi_3 - 1, 3)*std::pow(-2*xi_2 - xi_3 + 1, 7))/(xi_3 - 1) ;
}
if( b == 230 ) {
  val = (1.0L/256.0L)*(-1260*xi_2*std::pow(xi_3 - 1, 10) - 4095*std::pow(-xi_3 + 1, 9)*std::pow(-2*xi_2 - xi_3 + 1, 2) + 49140*std::pow(-xi_3 + 1, 6)*std::pow(-2*xi_2 - xi_3 + 1, 5) - 139230*std::pow(-xi_3 + 1, 5)*std::pow(-2*xi_2 - xi_3 + 1, 6) + 41990*std::pow(-xi_3 + 1, 2)*std::pow(-2*xi_2 - xi_3 + 1, 9) - 88179*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 10) - 693*std::pow(xi_3 - 1, 11) - 10920*std::pow(xi_3 - 1, 8)*std::pow(-2*xi_2 - xi_3 + 1, 3) - 40950*std::pow(xi_3 - 1, 7)*std::pow(-2*xi_2 - xi_3 + 1, 4) - 79560*std::pow(xi_3 - 1, 4)*std::pow(-2*xi_2 - xi_3 + 1, 7) - 188955*std::pow(xi_3 - 1, 3)*std::pow(-2*xi_2 - xi_3 + 1, 8))/(xi_3 - 1) ;
}
if( b == 231 ) {
  val = -1.0L/128.0L*(2*xi_1 + xi_2 + xi_3 - 1)*(2394*xi_3 - 124355*std::pow(2*xi_3 - 1, 9) - 203490*std::pow(2*xi_3 - 1, 8) + 58140*std::pow(2*xi_3 - 1, 7) + 217056*std::pow(2*xi_3 - 1, 6) + 36414*std::pow(2*xi_3 - 1, 5) - 64260*std::pow(2*xi_3 - 1, 4) - 17220*std::pow(2*xi_3 - 1, 3) + 5040*std::pow(2*xi_3 - 1, 2) - 1239) ;
}
if( b == 232 ) {
  val = (1.0L/128.0L)*(5*xi_2 + xi_3 - 1)*(2*xi_1 + xi_2 + xi_3 - 1)*(2640*xi_3 + 159885*std::pow(2*xi_3 - 1, 8) + 348840*std::pow(2*xi_3 - 1, 7) + 135660*std::pow(2*xi_3 - 1, 6) - 162792*std::pow(2*xi_3 - 1, 5) - 117810*std::pow(2*xi_3 - 1, 4) + 4760*std::pow(2*xi_3 - 1, 3) + 14700*std::pow(2*xi_3 - 1, 2) - 1499) ;
}
if( b == 233 ) {
  val = -(36*xi_2*std::pow(xi_3 - 1, 2) + 17*std::pow(xi_3 - 1, 3) - 21*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2))*(2*xi_1 + xi_2 + xi_3 - 1)*(182*xi_3 + 10659*std::pow(2*xi_3 - 1, 7) + 27132*std::pow(2*xi_3 - 1, 6) + 20349*std::pow(2*xi_3 - 1, 5) - 5355*std::pow(2*xi_3 - 1, 3) - 1428*std::pow(2*xi_3 - 1, 2) - 59)/(32*xi_3 - 32) ;
}
if( b == 234 ) {
  val = (1.0L/128.0L)*(2*std::pow(xi_3 - 1, 3) - 21*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2) + 21*std::pow(2*xi_2 + xi_3 - 1, 3))*(2*xi_1 + xi_2 + xi_3 - 1)*(-10812*xi_3 + 74613*std::pow(2*xi_3 - 1, 6) + 203490*std::pow(2*xi_3 - 1, 5) + 188955*std::pow(2*xi_3 - 1, 4) + 58140*std::pow(2*xi_3 - 1, 3) - 6885*std::pow(2*xi_3 - 1, 2) + 5011) ;
}
if( b == 235 ) {
  val = (1.0L/128.0L)*(2*xi_1 + xi_2 + xi_3 - 1)*(88*xi_2*std::pow(xi_3 - 1, 3) + 41*std::pow(xi_3 - 1, 4) - 18*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 2) - 180*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 3) + 165*std::pow(2*xi_2 + xi_3 - 1, 4))*(3150*xi_3 + 13167*std::pow(2*xi_3 - 1, 5) + 35910*std::pow(2*xi_3 - 1, 4) + 34770*std::pow(2*xi_3 - 1, 3) + 13680*std::pow(2*xi_3 - 1, 2) - 1669) ;
}
if( b == 236 ) {
  val = (1.0L/512.0L)*(2*xi_1 + xi_2 + xi_3 - 1)*(11704*xi_3 + 7315*std::pow(2*xi_3 - 1, 4) + 18620*std::pow(2*xi_3 - 1, 3) + 16530*std::pow(2*xi_3 - 1, 2) - 5209)*(-90*xi_2*std::pow(xi_3 - 1, 5) - 70*std::pow(-xi_3 + 1, 6) - 1485*std::pow(-xi_3 + 1, 2)*std::pow(-2*xi_2 - xi_3 + 1, 4) + 1287*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 5) + 630*std::pow(xi_3 - 1, 4)*std::pow(-2*xi_2 - xi_3 + 1, 2) + 330*std::pow(xi_3 - 1, 3)*std::pow(-2*xi_2 - xi_3 + 1, 3))/(xi_3 - 1) ;
}
if( b == 237 ) {
  val = -(2*xi_1 + xi_2 + xi_3 - 1)*(1170*xi_3 + 385*std::pow(2*xi_3 - 1, 3) + 840*std::pow(2*xi_3 - 1, 2) - 457)*(780*xi_2*std::pow(xi_3 - 1, 6) + 365*std::pow(xi_3 - 1, 7) + 165*std::pow(xi_3 - 1, 5)*std::pow(2*xi_2 + xi_3 - 1, 2) - 3740*std::pow(xi_3 - 1, 4)*std::pow(2*xi_2 + xi_3 - 1, 3) + 2145*std::pow(xi_3 - 1, 3)*std::pow(2*xi_2 + xi_3 - 1, 4) + 6006*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 5) - 5005*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 6))/(128*xi_3 - 128) ;
}
if( b == 238 ) {
  val = (1.0L/128.0L)*(756*xi_3 + 231*std::pow(2*xi_3 - 1, 2) - 227)*(2*xi_1 + xi_2 + xi_3 - 1)*(154*xi_2*std::pow(xi_3 - 1, 7) - 924*std::pow(-xi_3 + 1, 6)*std::pow(-2*xi_2 - xi_3 + 1, 2) - 6006*std::pow(-xi_3 + 1, 2)*std::pow(-2*xi_2 - xi_3 + 1, 6) + 4862*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 7) + 98*std::pow(xi_3 - 1, 8) + 5005*std::pow(xi_3 - 1, 4)*std::pow(-2*xi_2 - xi_3 + 1, 4) + 3003*std::pow(xi_3 - 1, 3)*std::pow(-2*xi_2 - xi_3 + 1, 5))/(xi_3 - 1) ;
}
if( b == 239 ) {
  val = (1.0L/128.0L)*(22*xi_3 - 1)*(2*xi_1 + xi_2 + xi_3 - 1)*(1680*xi_2*std::pow(xi_3 - 1, 7) + 13832*std::pow(-xi_3 + 1, 5)*std::pow(2*xi_2 + xi_3 - 1, 3) - 30940*std::pow(-xi_3 + 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 6) + 47736*(-xi_3 + 1)*std::pow(2*xi_2 + xi_3 - 1, 7) + 791*std::pow(xi_3 - 1, 8) + 1092*std::pow(xi_3 - 1, 6)*std::pow(2*xi_2 + xi_3 - 1, 2) + 2730*std::pow(xi_3 - 1, 4)*std::pow(2*xi_2 + xi_3 - 1, 4) + 50232*std::pow(xi_3 - 1, 3)*std::pow(2*xi_2 + xi_3 - 1, 5) + 37791*std::pow(2*xi_2 + xi_3 - 1, 8)) ;
}
if( b == 240 ) {
  val = -(2*xi_1 + xi_2 + xi_3 - 1)*(1638*xi_2*std::pow(xi_3 - 1, 9) + 966*std::pow(xi_3 - 1, 10) - 9828*std::pow(xi_3 - 1, 8)*std::pow(2*xi_2 + xi_3 - 1, 2) - 5460*std::pow(xi_3 - 1, 7)*std::pow(2*xi_2 + xi_3 - 1, 3) + 90090*std::pow(xi_3 - 1, 6)*std::pow(2*xi_2 + xi_3 - 1, 4) - 27846*std::pow(xi_3 - 1, 5)*std::pow(2*xi_2 + xi_3 - 1, 5) - 241332*std::pow(xi_3 - 1, 4)*std::pow(2*xi_2 + xi_3 - 1, 6) + 151164*std::pow(xi_3 - 1, 3)*std::pow(2*xi_2 + xi_3 - 1, 7) + 188955*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 8) - 146965*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 9))/(256*xi_3 - 256) ;
}
if( b == 241 ) {
  val = (1.0L/256.0L)*(-std::pow(xi_2 + xi_3 - 1, 2) + 3*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2))*(2640*xi_3 + 159885*std::pow(2*xi_3 - 1, 8) + 348840*std::pow(2*xi_3 - 1, 7) + 135660*std::pow(2*xi_3 - 1, 6) - 162792*std::pow(2*xi_3 - 1, 5) - 117810*std::pow(2*xi_3 - 1, 4) + 4760*std::pow(2*xi_3 - 1, 3) + 14700*std::pow(2*xi_3 - 1, 2) - 1499) ;
}
if( b == 242 ) {
  val = (1.0L/16.0L)*(-std::pow(xi_2 + xi_3 - 1, 2) + 3*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2))*(7*xi_2 + xi_3 - 1)*(182*xi_3 + 10659*std::pow(2*xi_3 - 1, 7) + 27132*std::pow(2*xi_3 - 1, 6) + 20349*std::pow(2*xi_3 - 1, 5) - 5355*std::pow(2*xi_3 - 1, 3) - 1428*std::pow(2*xi_3 - 1, 2) - 59) ;
}
if( b == 243 ) {
  val = (std::pow(xi_2 + xi_3 - 1, 2) - 3*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2))*(20*xi_2*std::pow(xi_3 - 1, 2) + 8*std::pow(xi_3 - 1, 3) - 9*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2))*(-10812*xi_3 + 74613*std::pow(2*xi_3 - 1, 6) + 203490*std::pow(2*xi_3 - 1, 5) + 188955*std::pow(2*xi_3 - 1, 4) + 58140*std::pow(2*xi_3 - 1, 3) - 6885*std::pow(2*xi_3 - 1, 2) + 5011)/(128*xi_3 - 128) ;
}
if( b == 244 ) {
  val = (1.0L/256.0L)*(-std::pow(xi_2 + xi_3 - 1, 2) + 3*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2))*(126*xi_2*std::pow(xi_3 - 1, 3) - 225*std::pow(-xi_3 + 1, 2)*std::pow(-2*xi_2 - xi_3 + 1, 2) + 165*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 3) + 68*std::pow(xi_3 - 1, 4))*(3150*xi_3 + 13167*std::pow(2*xi_3 - 1, 5) + 35910*std::pow(2*xi_3 - 1, 4) + 34770*std::pow(2*xi_3 - 1, 3) + 13680*std::pow(2*xi_3 - 1, 2) - 1669)/(xi_3 - 1) ;
}
if( b == 245 ) {
  val = (1.0L/512.0L)*(-std::pow(xi_2 + xi_3 - 1, 2) + 3*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2))*(11704*xi_3 + 7315*std::pow(2*xi_3 - 1, 4) + 18620*std::pow(2*xi_3 - 1, 3) + 16530*std::pow(2*xi_3 - 1, 2) - 5209)*(200*xi_2*std::pow(xi_3 - 1, 3) + 71*std::pow(xi_3 - 1, 4) + 330*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 2) - 1100*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 3) + 715*std::pow(2*xi_2 + xi_3 - 1, 4)) ;
}
if( b == 246 ) {
  val = (std::pow(xi_2 + xi_3 - 1, 2) - 3*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2))*(1170*xi_3 + 385*std::pow(2*xi_3 - 1, 3) + 840*std::pow(2*xi_3 - 1, 2) - 457)*(770*xi_2*std::pow(xi_3 - 1, 5) + 386*std::pow(xi_3 - 1, 6) - 990*std::pow(xi_3 - 1, 4)*std::pow(2*xi_2 + xi_3 - 1, 2) - 1430*std::pow(xi_3 - 1, 3)*std::pow(2*xi_2 + xi_3 - 1, 3) + 5005*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 4) - 3003*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 5))/(128*xi_3 - 128) ;
}
if( b == 247 ) {
  val = (1.0L/128.0L)*(std::pow(xi_2 + xi_3 - 1, 2) - 3*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2))*(756*xi_3 + 231*std::pow(2*xi_3 - 1, 2) - 227)*(96*xi_2*std::pow(xi_3 - 1, 5) + 23*std::pow(xi_3 - 1, 6) + 780*std::pow(xi_3 - 1, 4)*std::pow(2*xi_2 + xi_3 - 1, 2) - 1820*std::pow(xi_3 - 1, 3)*std::pow(2*xi_2 + xi_3 - 1, 3) - 1365*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 4) + 5460*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 5) - 3094*std::pow(2*xi_2 + xi_3 - 1, 6)) ;
}
if( b == 248 ) {
  val = (1.0L/64.0L)*(22*xi_3 - 1)*(-std::pow(xi_2 + xi_3 - 1, 2) + 3*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2))*(910*xi_2*std::pow(xi_3 - 1, 7) - 819*std::pow(-xi_3 + 1, 6)*std::pow(-2*xi_2 - xi_3 + 1, 2) + 4641*std::pow(-xi_3 + 1, 3)*std::pow(-2*xi_2 - xi_3 + 1, 5) - 23205*std::pow(-xi_3 + 1, 2)*std::pow(-2*xi_2 - xi_3 + 1, 6) + 12597*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 7) + 448*std::pow(xi_3 - 1, 8) + 5005*std::pow(xi_3 - 1, 5)*std::pow(-2*xi_2 - xi_3 + 1, 3) + 11375*std::pow(xi_3 - 1, 4)*std::pow(-2*xi_2 - xi_3 + 1, 4))/(xi_3 - 1) ;
}
if( b == 249 ) {
  val = (1.0L/256.0L)*(-std::pow(xi_2 + xi_3 - 1, 2) + 3*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2))*(112*xi_2*std::pow(xi_3 - 1, 7) - 119*std::pow(-xi_3 + 1, 8) + 17080*std::pow(-xi_3 + 1, 5)*std::pow(2*xi_2 + xi_3 - 1, 3) - 55930*std::pow(-xi_3 + 1, 4)*std::pow(2*xi_2 + xi_3 - 1, 4) + 193800*(-xi_3 + 1)*std::pow(2*xi_2 + xi_3 - 1, 7) + 9660*std::pow(xi_3 - 1, 6)*std::pow(2*xi_2 + xi_3 - 1, 2) + 128520*std::pow(xi_3 - 1, 3)*std::pow(2*xi_2 + xi_3 - 1, 5) + 27132*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 6) + 101745*std::pow(2*xi_2 + xi_3 - 1, 8)) ;
}
if( b == 250 ) {
  val = -(6*xi_1*std::pow(xi_2 + xi_3 - 1, 3) + 3*std::pow(xi_2 + xi_3 - 1, 4) - 5*(xi_2 + xi_3 - 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3))*(182*xi_3 + 10659*std::pow(2*xi_3 - 1, 7) + 27132*std::pow(2*xi_3 - 1, 6) + 20349*std::pow(2*xi_3 - 1, 5) - 5355*std::pow(2*xi_3 - 1, 3) - 1428*std::pow(2*xi_3 - 1, 2) - 59)/(16*xi_2 + 16*xi_3 - 16) ;
}
if( b == 251 ) {
  val = (1.0L/128.0L)*(9*xi_2 + xi_3 - 1)*(6*xi_1*std::pow(xi_2 + xi_3 - 1, 3) + 5*(-xi_2 - xi_3 + 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3) + 3*std::pow(xi_2 + xi_3 - 1, 4))*(-10812*xi_3 + 74613*std::pow(2*xi_3 - 1, 6) + 203490*std::pow(2*xi_3 - 1, 5) + 188955*std::pow(2*xi_3 - 1, 4) + 58140*std::pow(2*xi_3 - 1, 3) - 6885*std::pow(2*xi_3 - 1, 2) + 5011)/(-xi_2 - xi_3 + 1) ;
}
if( b == 252 ) {
  val = (1.0L/128.0L)*(6*xi_1*std::pow(xi_2 + xi_3 - 1, 3) + 5*(-xi_2 - xi_3 + 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3) + 3*std::pow(xi_2 + xi_3 - 1, 4))*(140*xi_2*std::pow(xi_3 - 1, 2) + 55*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 2) + 51*std::pow(xi_3 - 1, 3))*(3150*xi_3 + 13167*std::pow(2*xi_3 - 1, 5) + 35910*std::pow(2*xi_3 - 1, 4) + 34770*std::pow(2*xi_3 - 1, 3) + 13680*std::pow(2*xi_3 - 1, 2) - 1669)/((xi_3 - 1)*(xi_2 + xi_3 - 1)) ;
}
if( b == 253 ) {
  val = -(6*xi_1*std::pow(xi_2 + xi_3 - 1, 3) + 3*std::pow(xi_2 + xi_3 - 1, 4) - 5*(xi_2 + xi_3 - 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3))*(198*xi_2*std::pow(xi_3 - 1, 2) + 92*std::pow(xi_3 - 1, 3) - 231*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2) + 143*std::pow(2*xi_2 + xi_3 - 1, 3))*(11704*xi_3 + 7315*std::pow(2*xi_3 - 1, 4) + 18620*std::pow(2*xi_3 - 1, 3) + 16530*std::pow(2*xi_3 - 1, 2) - 5209)/(128*xi_2 + 128*xi_3 - 128) ;
}
if( b == 254 ) {
  val = (6*xi_1*std::pow(xi_2 + xi_3 - 1, 3) + 3*std::pow(xi_2 + xi_3 - 1, 4) - 5*(xi_2 + xi_3 - 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3))*(1170*xi_3 + 385*std::pow(2*xi_3 - 1, 3) + 840*std::pow(2*xi_3 - 1, 2) - 457)*(168*xi_2*std::pow(xi_3 - 1, 3) + 127*std::pow(xi_3 - 1, 4) - 1326*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 2) + 2548*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 3) - 1365*std::pow(2*xi_2 + xi_3 - 1, 4))/(64*xi_2 + 64*xi_3 - 64) ;
}
if( b == 255 ) {
  val = (1.0L/64.0L)*(756*xi_3 + 231*std::pow(2*xi_3 - 1, 2) - 227)*(6*xi_1*std::pow(xi_2 + xi_3 - 1, 3) + 5*(-xi_2 - xi_3 + 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3) + 3*std::pow(xi_2 + xi_3 - 1, 4))*(390*xi_2*std::pow(xi_3 - 1, 5) + 1820*std::pow(-xi_3 + 1, 3)*std::pow(2*xi_2 + xi_3 - 1, 3) + 1547*(-xi_3 + 1)*std::pow(2*xi_2 + xi_3 - 1, 5) + 174*std::pow(xi_3 - 1, 6) + 3185*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 4))/((xi_3 - 1)*(xi_2 + xi_3 - 1)) ;
}
if( b == 256 ) {
  val = -(22*xi_3 - 1)*(6*xi_1*std::pow(xi_2 + xi_3 - 1, 3) + 3*std::pow(xi_2 + xi_3 - 1, 4) - 5*(xi_2 + xi_3 - 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3))*(588*xi_2*std::pow(xi_3 - 1, 5) + 317*std::pow(xi_3 - 1, 6) - 1995*std::pow(xi_3 - 1, 4)*std::pow(2*xi_2 + xi_3 - 1, 2) + 980*std::pow(xi_3 - 1, 3)*std::pow(2*xi_2 + xi_3 - 1, 3) + 8925*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 4) - 14994*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 5) + 6783*std::pow(2*xi_2 + xi_3 - 1, 6))/(32*xi_2 + 32*xi_3 - 32) ;
}
if( b == 257 ) {
  val = -1.0L/32.0L*(6*xi_1*std::pow(xi_2 + xi_3 - 1, 3) + 3*std::pow(xi_2 + xi_3 - 1, 4) - 5*(xi_2 + xi_3 - 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3))*(630*xi_2*std::pow(xi_3 - 1, 7) + 272*std::pow(xi_3 - 1, 8) + 1155*std::pow(xi_3 - 1, 6)*std::pow(2*xi_2 + xi_3 - 1, 2) - 7735*std::pow(xi_3 - 1, 5)*std::pow(2*xi_2 + xi_3 - 1, 3) + 5355*std::pow(xi_3 - 1, 4)*std::pow(2*xi_2 + xi_3 - 1, 4) + 20349*std::pow(xi_3 - 1, 3)*std::pow(2*xi_2 + xi_3 - 1, 5) - 33915*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 6) + 14535*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 7))/((xi_3 - 1)*(xi_2 + xi_3 - 1)) ;
}
if( b == 258 ) {
  val = (1.0L/512.0L)*(3*std::pow(xi_2 + xi_3 - 1, 4) - 30*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + 35*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 4))*(-10812*xi_3 + 74613*std::pow(2*xi_3 - 1, 6) + 203490*std::pow(2*xi_3 - 1, 5) + 188955*std::pow(2*xi_3 - 1, 4) + 58140*std::pow(2*xi_3 - 1, 3) - 6885*std::pow(2*xi_3 - 1, 2) + 5011) ;
}
if( b == 259 ) {
  val = (1.0L/128.0L)*(11*xi_2 + xi_3 - 1)*(3*std::pow(xi_2 + xi_3 - 1, 4) - 30*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + 35*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 4))*(3150*xi_3 + 13167*std::pow(2*xi_3 - 1, 5) + 35910*std::pow(2*xi_3 - 1, 4) + 34770*std::pow(2*xi_3 - 1, 3) + 13680*std::pow(2*xi_3 - 1, 2) - 1669) ;
}
if( b == 260 ) {
  val = -(108*xi_2*std::pow(xi_3 - 1, 2) + 37*std::pow(xi_3 - 1, 3) - 39*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2))*(3*std::pow(xi_2 + xi_3 - 1, 4) - 30*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + 35*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 4))*(11704*xi_3 + 7315*std::pow(2*xi_3 - 1, 4) + 18620*std::pow(2*xi_3 - 1, 3) + 16530*std::pow(2*xi_3 - 1, 2) - 5209)/(256*xi_3 - 256) ;
}
if( b == 261 ) {
  val = (1.0L/128.0L)*(3*std::pow(xi_2 + xi_3 - 1, 4) - 30*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + 35*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 4))*(1170*xi_3 + 385*std::pow(2*xi_3 - 1, 3) + 840*std::pow(2*xi_3 - 1, 2) - 457)*(858*xi_2*std::pow(xi_3 - 1, 3) - 819*std::pow(-xi_3 + 1, 2)*std::pow(-2*xi_2 - xi_3 + 1, 2) + 455*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 3) + 372*std::pow(xi_3 - 1, 4))/(xi_3 - 1) ;
}
if( b == 262 ) {
  val = (1.0L/128.0L)*(756*xi_3 + 231*std::pow(2*xi_3 - 1, 2) - 227)*(3*std::pow(xi_2 + xi_3 - 1, 4) - 30*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + 35*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 4))*(-336*xi_2*std::pow(xi_3 - 1, 3) - 171*std::pow(xi_3 - 1, 4) + 840*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 2) - 1260*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 3) + 595*std::pow(2*xi_2 + xi_3 - 1, 4)) ;
}
if( b == 263 ) {
  val = -(22*xi_3 - 1)*(3*std::pow(xi_2 + xi_3 - 1, 4) - 30*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + 35*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 4))*(330*xi_2*std::pow(xi_3 - 1, 5) + 114*std::pow(xi_3 - 1, 6) + 1170*std::pow(xi_3 - 1, 4)*std::pow(2*xi_2 + xi_3 - 1, 2) - 5270*std::pow(xi_3 - 1, 3)*std::pow(2*xi_2 + xi_3 - 1, 3) + 6885*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 4) - 2907*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 5))/(64*xi_3 - 64) ;
}
if( b == 264 ) {
  val = (1.0L/64.0L)*(3*std::pow(xi_2 + xi_3 - 1, 4) - 30*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + 35*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 4))*(972*xi_2*std::pow(xi_3 - 1, 5) + 467*std::pow(xi_3 - 1, 6) - 1275*std::pow(xi_3 - 1, 4)*std::pow(2*xi_2 + xi_3 - 1, 2) - 3060*std::pow(xi_3 - 1, 3)*std::pow(2*xi_2 + xi_3 - 1, 3) + 14535*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 4) - 17442*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 5) + 6783*std::pow(2*xi_2 + xi_3 - 1, 6)) ;
}
if( b == 265 ) {
  val = (1.0L/128.0L)*(30*xi_1*std::pow(xi_2 + xi_3 - 1, 5) + 63*(-xi_2 - xi_3 + 1)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 5) + 15*std::pow(xi_2 + xi_3 - 1, 6) + 70*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 3))*(3150*xi_3 + 13167*std::pow(2*xi_3 - 1, 5) + 35910*std::pow(2*xi_3 - 1, 4) + 34770*std::pow(2*xi_3 - 1, 3) + 13680*std::pow(2*xi_3 - 1, 2) - 1669)/(xi_2 + xi_3 - 1) ;
}
if( b == 266 ) {
  val = (1.0L/128.0L)*(13*xi_2 + xi_3 - 1)*(30*xi_1*std::pow(xi_2 + xi_3 - 1, 5) + 63*(-xi_2 - xi_3 + 1)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 5) + 15*std::pow(xi_2 + xi_3 - 1, 6) + 70*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 3))*(11704*xi_3 + 7315*std::pow(2*xi_3 - 1, 4) + 18620*std::pow(2*xi_3 - 1, 3) + 16530*std::pow(2*xi_3 - 1, 2) - 5209)/(xi_2 + xi_3 - 1) ;
}
if( b == 267 ) {
  val = -1.0L/64.0L*(308*xi_2*std::pow(xi_3 - 1, 2) + 101*std::pow(xi_3 - 1, 3) - 105*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2))*(1170*xi_3 + 385*std::pow(2*xi_3 - 1, 3) + 840*std::pow(2*xi_3 - 1, 2) - 457)*(30*xi_1*std::pow(xi_2 + xi_3 - 1, 5) + 15*std::pow(xi_2 + xi_3 - 1, 6) - 70*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3) + 63*(xi_2 + xi_3 - 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 5))/((xi_3 - 1)*(xi_2 + xi_3 - 1)) ;
}
if( b == 268 ) {
  val = (1.0L/32.0L)*(756*xi_3 + 231*std::pow(2*xi_3 - 1, 2) - 227)*(30*xi_1*std::pow(xi_2 + xi_3 - 1, 5) + 63*(-xi_2 - xi_3 + 1)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 5) + 15*std::pow(xi_2 + xi_3 - 1, 6) + 70*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 3))*(195*xi_2*std::pow(xi_3 - 1, 2) + 81*std::pow(xi_3 - 1, 3) - 165*(xi_3 - 1)*std::pow(-2*xi_2 - xi_3 + 1, 2) - 85*std::pow(-2*xi_2 - xi_3 + 1, 3))/(xi_2 + xi_3 - 1) ;
}
if( b == 269 ) {
  val = (1.0L/32.0L)*(22*xi_3 - 1)*(30*xi_1*std::pow(xi_2 + xi_3 - 1, 5) + 63*(-xi_2 - xi_3 + 1)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 5) + 15*std::pow(xi_2 + xi_3 - 1, 6) + 70*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 3))*(-968*xi_2*std::pow(xi_3 - 1, 3) - 455*std::pow(xi_3 - 1, 4) + 1734*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 2) - 2244*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 3) + 969*std::pow(2*xi_2 + xi_3 - 1, 4))/(xi_2 + xi_3 - 1) ;
}
if( b == 270 ) {
  val = (1.0L/256.0L)*(30*xi_1*std::pow(xi_2 + xi_3 - 1, 5) + 15*std::pow(xi_2 + xi_3 - 1, 6) - 70*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3) + 63*(xi_2 + xi_3 - 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 5))*(2210*xi_2*std::pow(xi_3 - 1, 5) + 1358*std::pow(xi_3 - 1, 6) - 16830*std::pow(xi_3 - 1, 4)*std::pow(2*xi_2 + xi_3 - 1, 2) + 48450*std::pow(xi_3 - 1, 3)*std::pow(2*xi_2 + xi_3 - 1, 3) - 53295*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 4) + 20349*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 5))/((xi_3 - 1)*(xi_2 + xi_3 - 1)) ;
}
if( b == 271 ) {
  val = (1.0L/256.0L)*(-5*std::pow(xi_2 + xi_3 - 1, 6) + 105*std::pow(xi_2 + xi_3 - 1, 4)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) - 315*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 4) + 231*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 6))*(11704*xi_3 + 7315*std::pow(2*xi_3 - 1, 4) + 18620*std::pow(2*xi_3 - 1, 3) + 16530*std::pow(2*xi_3 - 1, 2) - 5209) ;
}
if( b == 272 ) {
  val = (1.0L/32.0L)*(15*xi_2 + xi_3 - 1)*(1170*xi_3 + 385*std::pow(2*xi_3 - 1, 3) + 840*std::pow(2*xi_3 - 1, 2) - 457)*(-5*std::pow(xi_2 + xi_3 - 1, 6) + 105*std::pow(xi_2 + xi_3 - 1, 4)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) - 315*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 4) + 231*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 6)) ;
}
if( b == 273 ) {
  val = (756*xi_3 + 231*std::pow(2*xi_3 - 1, 2) - 227)*(104*xi_2*std::pow(xi_3 - 1, 2) + 33*std::pow(xi_3 - 1, 3) - 34*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2))*(5*std::pow(xi_2 + xi_3 - 1, 6) - 105*std::pow(xi_2 + xi_3 - 1, 4)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + 315*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 4) - 231*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 6))/(64*xi_3 - 64) ;
}
if( b == 274 ) {
  val = (1.0L/128.0L)*(22*xi_3 - 1)*(2550*xi_2*std::pow(xi_3 - 1, 3) - 1989*std::pow(-xi_3 + 1, 2)*std::pow(-2*xi_2 - xi_3 + 1, 2) + 969*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 3) + 1028*std::pow(xi_3 - 1, 4))*(-5*std::pow(xi_2 + xi_3 - 1, 6) + 105*std::pow(xi_2 + xi_3 - 1, 4)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) - 315*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 4) + 231*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 6))/(xi_3 - 1) ;
}
if( b == 275 ) {
  val = (1.0L/256.0L)*(-5*std::pow(xi_2 + xi_3 - 1, 6) + 105*std::pow(xi_2 + xi_3 - 1, 4)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) - 315*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 4) + 231*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 6))*(-8424*xi_2*std::pow(xi_3 - 1, 3) - 3803*std::pow(xi_3 - 1, 4) + 12654*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 2) - 14820*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 3) + 5985*std::pow(2*xi_2 + xi_3 - 1, 4)) ;
}
if( b == 276 ) {
  val = -(1170*xi_3 + 385*std::pow(2*xi_3 - 1, 3) + 840*std::pow(2*xi_3 - 1, 2) - 457)*(70*xi_1*std::pow(xi_2 + xi_3 - 1, 7) + 35*std::pow(xi_2 + xi_3 - 1, 8) - 315*std::pow(xi_2 + xi_3 - 1, 5)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3) + 693*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 5) - 429*(xi_2 + xi_3 - 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 7))/(32*xi_2 + 32*xi_3 - 32) ;
}
if( b == 277 ) {
  val = -(17*xi_2 + xi_3 - 1)*(756*xi_3 + 231*std::pow(2*xi_3 - 1, 2) - 227)*(70*xi_1*std::pow(xi_2 + xi_3 - 1, 7) + 35*std::pow(xi_2 + xi_3 - 1, 8) - 315*std::pow(xi_2 + xi_3 - 1, 5)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3) + 693*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 5) - 429*(xi_2 + xi_3 - 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 7))/(64*xi_2 + 64*xi_3 - 64) ;
}
if( b == 278 ) {
  val = (1.0L/64.0L)*(22*xi_3 - 1)*(540*xi_2*std::pow(xi_3 - 1, 2) + 171*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 2) + 167*std::pow(xi_3 - 1, 3))*(70*xi_1*std::pow(xi_2 + xi_3 - 1, 7) + 315*std::pow(-xi_2 - xi_3 + 1, 5)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3) + 429*(-xi_2 - xi_3 + 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 7) + 35*std::pow(xi_2 + xi_3 - 1, 8) + 693*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 5))/((xi_3 - 1)*(xi_2 + xi_3 - 1)) ;
}
if( b == 279 ) {
  val = -(1938*xi_2*std::pow(xi_3 - 1, 2) + 764*std::pow(xi_3 - 1, 3) - 1425*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2) + 665*std::pow(2*xi_2 + xi_3 - 1, 3))*(70*xi_1*std::pow(xi_2 + xi_3 - 1, 7) + 35*std::pow(xi_2 + xi_3 - 1, 8) - 315*std::pow(xi_2 + xi_3 - 1, 5)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3) + 693*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 5) - 429*(xi_2 + xi_3 - 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 7))/(64*xi_2 + 64*xi_3 - 64) ;
}
if( b == 280 ) {
  val = (1.0L/512.0L)*(756*xi_3 + 231*std::pow(2*xi_3 - 1, 2) - 227)*(35*std::pow(xi_2 + xi_3 - 1, 8) - 1260*std::pow(xi_2 + xi_3 - 1, 6)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + 6930*std::pow(xi_2 + xi_3 - 1, 4)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 4) - 12012*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 6) + 6435*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 8)) ;
}
if( b == 281 ) {
  val = (1.0L/128.0L)*(22*xi_3 - 1)*(19*xi_2 + xi_3 - 1)*(35*std::pow(xi_2 + xi_3 - 1, 8) - 1260*std::pow(xi_2 + xi_3 - 1, 6)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + 6930*std::pow(xi_2 + xi_3 - 1, 4)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 4) - 12012*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 6) + 6435*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 8)) ;
}
if( b == 282 ) {
  val = -(340*xi_2*std::pow(xi_3 - 1, 2) + 103*std::pow(xi_3 - 1, 3) - 105*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2))*(35*std::pow(xi_2 + xi_3 - 1, 8) - 1260*std::pow(xi_2 + xi_3 - 1, 6)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + 6930*std::pow(xi_2 + xi_3 - 1, 4)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 4) - 12012*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 6) + 6435*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 8))/(256*xi_3 - 256) ;
}
if( b == 283 ) {
  val = (1.0L/128.0L)*(22*xi_3 - 1)*(630*xi_1*std::pow(xi_2 + xi_3 - 1, 9) + 18018*std::pow(-xi_2 - xi_3 + 1, 5)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 5) + 12155*(-xi_2 - xi_3 + 1)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 9) + 315*std::pow(xi_2 + xi_3 - 1, 10) + 4620*std::pow(xi_2 + xi_3 - 1, 7)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 3) + 25740*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 7))/(xi_2 + xi_3 - 1) ;
}
if( b == 284 ) {
  val = (1.0L/128.0L)*(21*xi_2 + xi_3 - 1)*(630*xi_1*std::pow(xi_2 + xi_3 - 1, 9) + 18018*std::pow(-xi_2 - xi_3 + 1, 5)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 5) + 12155*(-xi_2 - xi_3 + 1)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 9) + 315*std::pow(xi_2 + xi_3 - 1, 10) + 4620*std::pow(xi_2 + xi_3 - 1, 7)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 3) + 25740*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 7))/(xi_2 + xi_3 - 1) ;
}
if( b == 285 ) {
  val = -63.0L/256.0L*std::pow(xi_2 + xi_3 - 1, 10) + (3465.0L/256.0L)*std::pow(xi_2 + xi_3 - 1, 8)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) - 15015.0L/128.0L*std::pow(xi_2 + xi_3 - 1, 6)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 4) + (45045.0L/128.0L)*std::pow(xi_2 + xi_3 - 1, 4)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 6) - 109395.0L/256.0L*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 8) + (46189.0L/256.0L)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 10) ;
}
