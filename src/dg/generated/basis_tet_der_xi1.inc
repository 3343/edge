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
 * Generated derivatives of tetrahedral basis w.r.t. to xi1.
 **/

assert( b < 286 );
if( b == 0 ) {
  valDxi1 = 0 ;
}
if( b == 1 ) {
  valDxi1 = 0 ;
}
if( b == 2 ) {
  valDxi1 = 0 ;
}
if( b == 3 ) {
  valDxi1 = 2 ;
}
if( b == 4 ) {
  valDxi1 = 0 ;
}
if( b == 5 ) {
  valDxi1 = 0 ;
}
if( b == 6 ) {
  valDxi1 = 0 ;
}
if( b == 7 ) {
  valDxi1 = 12*xi_3 - 2 ;
}
if( b == 8 ) {
  valDxi1 = 10*xi_2 + 2*xi_3 - 2 ;
}
if( b == 9 ) {
  valDxi1 = 12*xi_1 + 6*xi_2 + 6*xi_3 - 6 ;
}
if( b == 10 ) {
  valDxi1 = 0 ;
}
if( b == 11 ) {
  valDxi1 = 0 ;
}
if( b == 12 ) {
  valDxi1 = 0 ;
}
if( b == 13 ) {
  valDxi1 = 0 ;
}
if( b == 14 ) {
  valDxi1 = 28*xi_3 + 14*std::pow(2*xi_3 - 1, 2) - 12 ;
}
if( b == 15 ) {
  valDxi1 = 2*(8*xi_3 - 1)*(5*xi_2 + xi_3 - 1) ;
}
if( b == 16 ) {
  valDxi1 = -2*(36*xi_2*std::pow(xi_3 - 1, 2) + 17*std::pow(xi_3 - 1, 3) - 21*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2))/(4*xi_3 - 4) ;
}
if( b == 17 ) {
  valDxi1 = (1.0L/2.0L)*(8*xi_3 - 1)*(24*xi_1 + 12*xi_2 + 12*xi_3 - 12) ;
}
if( b == 18 ) {
  valDxi1 = (1.0L/2.0L)*(7*xi_2 + xi_3 - 1)*(24*xi_1 + 12*xi_2 + 12*xi_3 - 12) ;
}
if( b == 19 ) {
  valDxi1 = -(6*std::pow(xi_2 + xi_3 - 1, 3) - 30*(xi_2 + xi_3 - 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2))/(2*xi_2 + 2*xi_3 - 2) ;
}
if( b == 20 ) {
  valDxi1 = 0 ;
}
if( b == 21 ) {
  valDxi1 = 0 ;
}
if( b == 22 ) {
  valDxi1 = 0 ;
}
if( b == 23 ) {
  valDxi1 = 0 ;
}
if( b == 24 ) {
  valDxi1 = 0 ;
}
if( b == 25 ) {
  valDxi1 = 12*xi_3 + 30*std::pow(2*xi_3 - 1, 3) + 36*std::pow(2*xi_3 - 1, 2) - 8 ;
}
if( b == 26 ) {
  valDxi1 = (1.0L/2.0L)*(5*xi_2 + xi_3 - 1)*(108*xi_3 + 45*std::pow(2*xi_3 - 1, 2) - 41) ;
}
if( b == 27 ) {
  valDxi1 = -2*(10*xi_3 - 1)*(36*xi_2*std::pow(xi_3 - 1, 2) + 17*std::pow(xi_3 - 1, 3) - 21*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2))/(4*xi_3 - 4) ;
}
if( b == 28 ) {
  valDxi1 = 2*std::pow(xi_3 - 1, 3) - 21*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2) + 21*std::pow(2*xi_2 + xi_3 - 1, 3) ;
}
if( b == 29 ) {
  valDxi1 = (1.0L/8.0L)*(108*xi_3 + 45*std::pow(2*xi_3 - 1, 2) - 41)*(24*xi_1 + 12*xi_2 + 12*xi_3 - 12) ;
}
if( b == 30 ) {
  valDxi1 = (1.0L/2.0L)*(10*xi_3 - 1)*(7*xi_2 + xi_3 - 1)*(24*xi_1 + 12*xi_2 + 12*xi_3 - 12) ;
}
if( b == 31 ) {
  valDxi1 = (20*xi_2*std::pow(xi_3 - 1, 2) + 8*std::pow(xi_3 - 1, 3) - 9*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2))*(-24*xi_1 - 12*xi_2 - 12*xi_3 + 12)/(2*xi_3 - 2) ;
}
if( b == 32 ) {
  valDxi1 = -(10*xi_3 - 1)*(6*std::pow(xi_2 + xi_3 - 1, 3) - 30*(xi_2 + xi_3 - 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2))/(2*xi_2 + 2*xi_3 - 2) ;
}
if( b == 33 ) {
  valDxi1 = (1.0L/2.0L)*(30*(-xi_2 - xi_3 + 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + 6*std::pow(xi_2 + xi_3 - 1, 3))*(9*xi_2 + xi_3 - 1)/(-xi_2 - xi_3 + 1) ;
}
if( b == 34 ) {
  valDxi1 = -15.0L/4.0L*std::pow(xi_2 + xi_3 - 1, 2)*(8*xi_1 + 4*xi_2 + 4*xi_3 - 4) + 35*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3) ;
}
if( b == 35 ) {
  valDxi1 = 0 ;
}
if( b == 36 ) {
  valDxi1 = 0 ;
}
if( b == 37 ) {
  valDxi1 = 0 ;
}
if( b == 38 ) {
  valDxi1 = 0 ;
}
if( b == 39 ) {
  valDxi1 = 0 ;
}
if( b == 40 ) {
  valDxi1 = 0 ;
}
if( b == 41 ) {
  valDxi1 = -27*xi_3 + (495.0L/8.0L)*std::pow(2*xi_3 - 1, 4) + (165.0L/2.0L)*std::pow(2*xi_3 - 1, 3) + (45.0L/4.0L)*std::pow(2*xi_3 - 1, 2) + 91.0L/8.0L ;
}
if( b == 42 ) {
  valDxi1 = (1.0L/2.0L)*(5*xi_2 + xi_3 - 1)*(120*xi_3 + 110*std::pow(2*xi_3 - 1, 3) + 165*std::pow(2*xi_3 - 1, 2) - 59) ;
}
if( b == 43 ) {
  valDxi1 = -2*(88*xi_3 + 33*std::pow(2*xi_3 - 1, 2) - 31)*(36*xi_2*std::pow(xi_3 - 1, 2) + 17*std::pow(xi_3 - 1, 3) - 21*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2))/(8*xi_3 - 8) ;
}
if( b == 44 ) {
  valDxi1 = (12*xi_3 - 1)*(2*std::pow(xi_3 - 1, 3) - 21*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2) + 21*std::pow(2*xi_2 + xi_3 - 1, 3)) ;
}
if( b == 45 ) {
  valDxi1 = 22*xi_2*std::pow(xi_3 - 1, 3) + (41.0L/4.0L)*std::pow(xi_3 - 1, 4) - 9.0L/2.0L*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 2) - 45*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 3) + (165.0L/4.0L)*std::pow(2*xi_2 + xi_3 - 1, 4) ;
}
if( b == 46 ) {
  valDxi1 = (1.0L/8.0L)*(24*xi_1 + 12*xi_2 + 12*xi_3 - 12)*(120*xi_3 + 110*std::pow(2*xi_3 - 1, 3) + 165*std::pow(2*xi_3 - 1, 2) - 59) ;
}
if( b == 47 ) {
  valDxi1 = (1.0L/4.0L)*(7*xi_2 + xi_3 - 1)*(88*xi_3 + 33*std::pow(2*xi_3 - 1, 2) - 31)*(24*xi_1 + 12*xi_2 + 12*xi_3 - 12) ;
}
if( b == 48 ) {
  valDxi1 = (12*xi_3 - 1)*(20*xi_2*std::pow(xi_3 - 1, 2) + 8*std::pow(xi_3 - 1, 3) - 9*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2))*(-24*xi_1 - 12*xi_2 - 12*xi_3 + 12)/(2*xi_3 - 2) ;
}
if( b == 49 ) {
  valDxi1 = (1.0L/16.0L)*(24*xi_1 + 12*xi_2 + 12*xi_3 - 12)*(126*xi_2*std::pow(xi_3 - 1, 3) - 225*std::pow(-xi_3 + 1, 2)*std::pow(-2*xi_2 - xi_3 + 1, 2) + 165*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 3) + 68*std::pow(xi_3 - 1, 4))/(xi_3 - 1) ;
}
if( b == 50 ) {
  valDxi1 = -(6*std::pow(xi_2 + xi_3 - 1, 3) - 30*(xi_2 + xi_3 - 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2))*(88*xi_3 + 33*std::pow(2*xi_3 - 1, 2) - 31)/(4*xi_2 + 4*xi_3 - 4) ;
}
if( b == 51 ) {
  valDxi1 = (1.0L/2.0L)*(12*xi_3 - 1)*(30*(-xi_2 - xi_3 + 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + 6*std::pow(xi_2 + xi_3 - 1, 3))*(9*xi_2 + xi_3 - 1)/(-xi_2 - xi_3 + 1) ;
}
if( b == 52 ) {
  valDxi1 = (1.0L/8.0L)*(30*(-xi_2 - xi_3 + 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + 6*std::pow(xi_2 + xi_3 - 1, 3))*(140*xi_2*std::pow(xi_3 - 1, 2) + 55*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 2) + 51*std::pow(xi_3 - 1, 3))/((xi_3 - 1)*(xi_2 + xi_3 - 1)) ;
}
if( b == 53 ) {
  valDxi1 = (1.0L/8.0L)*(12*xi_3 - 1)*(-30*std::pow(xi_2 + xi_3 - 1, 2)*(8*xi_1 + 4*xi_2 + 4*xi_3 - 4) + 280*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3)) ;
}
if( b == 54 ) {
  valDxi1 = (1.0L/8.0L)*(-30*std::pow(xi_2 + xi_3 - 1, 2)*(8*xi_1 + 4*xi_2 + 4*xi_3 - 4) + 280*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3))*(11*xi_2 + xi_3 - 1) ;
}
if( b == 55 ) {
  valDxi1 = (1.0L/8.0L)*(-630*(-xi_2 - xi_3 + 1)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 4) + 30*std::pow(xi_2 + xi_3 - 1, 5) - 420*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 2))/(xi_2 + xi_3 - 1) ;
}
if( b == 56 ) {
  valDxi1 = 0 ;
}
if( b == 57 ) {
  valDxi1 = 0 ;
}
if( b == 58 ) {
  valDxi1 = 0 ;
}
if( b == 59 ) {
  valDxi1 = 0 ;
}
if( b == 60 ) {
  valDxi1 = 0 ;
}
if( b == 61 ) {
  valDxi1 = 0 ;
}
if( b == 62 ) {
  valDxi1 = 0 ;
}
if( b == 63 ) {
  valDxi1 = -95.0L/4.0L*xi_3 + (1001.0L/8.0L)*std::pow(2*xi_3 - 1, 5) + (715.0L/4.0L)*std::pow(2*xi_3 - 1, 4) + (55.0L/4.0L)*std::pow(2*xi_3 - 1, 3) - 55*std::pow(2*xi_3 - 1, 2) + 105.0L/8.0L ;
}
if( b == 64 ) {
  valDxi1 = (1.0L/8.0L)*(5*xi_2 + xi_3 - 1)*(-88*xi_3 + 1001*std::pow(2*xi_3 - 1, 4) + 1716*std::pow(2*xi_3 - 1, 3) + 726*std::pow(2*xi_3 - 1, 2) + 5) ;
}
if( b == 65 ) {
  valDxi1 = -2*(36*xi_2*std::pow(xi_3 - 1, 2) + 17*std::pow(xi_3 - 1, 3) - 21*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2))*(150*xi_3 + 91*std::pow(2*xi_3 - 1, 3) + 156*std::pow(2*xi_3 - 1, 2) - 67)/(8*xi_3 - 8) ;
}
if( b == 66 ) {
  valDxi1 = (1.0L/4.0L)*(260*xi_3 + 91*std::pow(2*xi_3 - 1, 2) - 87)*(2*std::pow(xi_3 - 1, 3) - 21*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2) + 21*std::pow(2*xi_2 + xi_3 - 1, 3)) ;
}
if( b == 67 ) {
  valDxi1 = (1.0L/4.0L)*(14*xi_3 - 1)*(88*xi_2*std::pow(xi_3 - 1, 3) + 41*std::pow(xi_3 - 1, 4) - 18*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 2) - 180*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 3) + 165*std::pow(2*xi_2 + xi_3 - 1, 4)) ;
}
if( b == 68 ) {
  valDxi1 = (1.0L/16.0L)*(-90*xi_2*std::pow(xi_3 - 1, 5) - 70*std::pow(-xi_3 + 1, 6) - 1485*std::pow(-xi_3 + 1, 2)*std::pow(-2*xi_2 - xi_3 + 1, 4) + 1287*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 5) + 630*std::pow(xi_3 - 1, 4)*std::pow(-2*xi_2 - xi_3 + 1, 2) + 330*std::pow(xi_3 - 1, 3)*std::pow(-2*xi_2 - xi_3 + 1, 3))/(xi_3 - 1) ;
}
if( b == 69 ) {
  valDxi1 = (1.0L/32.0L)*(24*xi_1 + 12*xi_2 + 12*xi_3 - 12)*(-88*xi_3 + 1001*std::pow(2*xi_3 - 1, 4) + 1716*std::pow(2*xi_3 - 1, 3) + 726*std::pow(2*xi_3 - 1, 2) + 5) ;
}
if( b == 70 ) {
  valDxi1 = (1.0L/4.0L)*(7*xi_2 + xi_3 - 1)*(24*xi_1 + 12*xi_2 + 12*xi_3 - 12)*(150*xi_3 + 91*std::pow(2*xi_3 - 1, 3) + 156*std::pow(2*xi_3 - 1, 2) - 67) ;
}
if( b == 71 ) {
  valDxi1 = (260*xi_3 + 91*std::pow(2*xi_3 - 1, 2) - 87)*(20*xi_2*std::pow(xi_3 - 1, 2) + 8*std::pow(xi_3 - 1, 3) - 9*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2))*(-24*xi_1 - 12*xi_2 - 12*xi_3 + 12)/(8*xi_3 - 8) ;
}
if( b == 72 ) {
  valDxi1 = (1.0L/16.0L)*(14*xi_3 - 1)*(24*xi_1 + 12*xi_2 + 12*xi_3 - 12)*(126*xi_2*std::pow(xi_3 - 1, 3) - 225*std::pow(-xi_3 + 1, 2)*std::pow(-2*xi_2 - xi_3 + 1, 2) + 165*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 3) + 68*std::pow(xi_3 - 1, 4))/(xi_3 - 1) ;
}
if( b == 73 ) {
  valDxi1 = (1.0L/32.0L)*(24*xi_1 + 12*xi_2 + 12*xi_3 - 12)*(200*xi_2*std::pow(xi_3 - 1, 3) + 71*std::pow(xi_3 - 1, 4) + 330*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 2) - 1100*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 3) + 715*std::pow(2*xi_2 + xi_3 - 1, 4)) ;
}
if( b == 74 ) {
  valDxi1 = -(6*std::pow(xi_2 + xi_3 - 1, 3) - 30*(xi_2 + xi_3 - 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2))*(150*xi_3 + 91*std::pow(2*xi_3 - 1, 3) + 156*std::pow(2*xi_3 - 1, 2) - 67)/(4*xi_2 + 4*xi_3 - 4) ;
}
if( b == 75 ) {
  valDxi1 = (1.0L/8.0L)*(30*(-xi_2 - xi_3 + 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + 6*std::pow(xi_2 + xi_3 - 1, 3))*(9*xi_2 + xi_3 - 1)*(260*xi_3 + 91*std::pow(2*xi_3 - 1, 2) - 87)/(-xi_2 - xi_3 + 1) ;
}
if( b == 76 ) {
  valDxi1 = (1.0L/8.0L)*(14*xi_3 - 1)*(30*(-xi_2 - xi_3 + 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + 6*std::pow(xi_2 + xi_3 - 1, 3))*(140*xi_2*std::pow(xi_3 - 1, 2) + 55*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 2) + 51*std::pow(xi_3 - 1, 3))/((xi_3 - 1)*(xi_2 + xi_3 - 1)) ;
}
if( b == 77 ) {
  valDxi1 = -(6*std::pow(xi_2 + xi_3 - 1, 3) - 30*(xi_2 + xi_3 - 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2))*(198*xi_2*std::pow(xi_3 - 1, 2) + 92*std::pow(xi_3 - 1, 3) - 231*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2) + 143*std::pow(2*xi_2 + xi_3 - 1, 3))/(8*xi_2 + 8*xi_3 - 8) ;
}
if( b == 78 ) {
  valDxi1 = (1.0L/32.0L)*(-30*std::pow(xi_2 + xi_3 - 1, 2)*(8*xi_1 + 4*xi_2 + 4*xi_3 - 4) + 280*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3))*(260*xi_3 + 91*std::pow(2*xi_3 - 1, 2) - 87) ;
}
if( b == 79 ) {
  valDxi1 = (1.0L/8.0L)*(14*xi_3 - 1)*(-30*std::pow(xi_2 + xi_3 - 1, 2)*(8*xi_1 + 4*xi_2 + 4*xi_3 - 4) + 280*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3))*(11*xi_2 + xi_3 - 1) ;
}
if( b == 80 ) {
  valDxi1 = -(-30*std::pow(xi_2 + xi_3 - 1, 2)*(8*xi_1 + 4*xi_2 + 4*xi_3 - 4) + 280*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3))*(108*xi_2*std::pow(xi_3 - 1, 2) + 37*std::pow(xi_3 - 1, 3) - 39*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2))/(16*xi_3 - 16) ;
}
if( b == 81 ) {
  valDxi1 = (1.0L/8.0L)*(14*xi_3 - 1)*(-630*(-xi_2 - xi_3 + 1)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 4) + 30*std::pow(xi_2 + xi_3 - 1, 5) - 420*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 2))/(xi_2 + xi_3 - 1) ;
}
if( b == 82 ) {
  valDxi1 = (1.0L/8.0L)*(13*xi_2 + xi_3 - 1)*(-630*(-xi_2 - xi_3 + 1)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 4) + 30*std::pow(xi_2 + xi_3 - 1, 5) - 420*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 2))/(xi_2 + xi_3 - 1) ;
}
if( b == 83 ) {
  valDxi1 = (105.0L/16.0L)*std::pow(xi_2 + xi_3 - 1, 4)*(8*xi_1 + 4*xi_2 + 4*xi_3 - 4) - 315.0L/2.0L*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3) + (693.0L/4.0L)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 5) ;
}
if( b == 84 ) {
  valDxi1 = 0 ;
}
if( b == 85 ) {
  valDxi1 = 0 ;
}
if( b == 86 ) {
  valDxi1 = 0 ;
}
if( b == 87 ) {
  valDxi1 = 0 ;
}
if( b == 88 ) {
  valDxi1 = 0 ;
}
if( b == 89 ) {
  valDxi1 = 0 ;
}
if( b == 90 ) {
  valDxi1 = 0 ;
}
if( b == 91 ) {
  valDxi1 = 0 ;
}
if( b == 92 ) {
  valDxi1 = (99.0L/4.0L)*xi_3 + (1001.0L/4.0L)*std::pow(2*xi_3 - 1, 6) + (3003.0L/8.0L)*std::pow(2*xi_3 - 1, 5) - 715.0L/4.0L*std::pow(2*xi_3 - 1, 3) - 165.0L/4.0L*std::pow(2*xi_3 - 1, 2) - 83.0L/8.0L ;
}
if( b == 93 ) {
  valDxi1 = (1.0L/8.0L)*(5*xi_2 + xi_3 - 1)*(-600*xi_3 + 2184*std::pow(2*xi_3 - 1, 5) + 4095*std::pow(2*xi_3 - 1, 4) + 1820*std::pow(2*xi_3 - 1, 3) - 390*std::pow(2*xi_3 - 1, 2) + 283) ;
}
if( b == 94 ) {
  valDxi1 = -2*(36*xi_2*std::pow(xi_3 - 1, 2) + 17*std::pow(xi_3 - 1, 3) - 21*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2))*(156*xi_3 + 455*std::pow(2*xi_3 - 1, 4) + 910*std::pow(2*xi_3 - 1, 3) + 546*std::pow(2*xi_3 - 1, 2) - 87)/(16*xi_3 - 16) ;
}
if( b == 95 ) {
  valDxi1 = (1.0L/4.0L)*(2*std::pow(xi_3 - 1, 3) - 21*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2) + 21*std::pow(2*xi_2 + xi_3 - 1, 3))*(588*xi_3 + 280*std::pow(2*xi_3 - 1, 3) + 525*std::pow(2*xi_3 - 1, 2) - 249) ;
}
if( b == 96 ) {
  valDxi1 = (1.0L/4.0L)*(90*xi_3 + 30*std::pow(2*xi_3 - 1, 2) - 29)*(88*xi_2*std::pow(xi_3 - 1, 3) + 41*std::pow(xi_3 - 1, 4) - 18*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 2) - 180*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 3) + 165*std::pow(2*xi_2 + xi_3 - 1, 4)) ;
}
if( b == 97 ) {
  valDxi1 = (1.0L/16.0L)*(16*xi_3 - 1)*(-90*xi_2*std::pow(xi_3 - 1, 5) - 70*std::pow(-xi_3 + 1, 6) - 1485*std::pow(-xi_3 + 1, 2)*std::pow(-2*xi_2 - xi_3 + 1, 4) + 1287*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 5) + 630*std::pow(xi_3 - 1, 4)*std::pow(-2*xi_2 - xi_3 + 1, 2) + 330*std::pow(xi_3 - 1, 3)*std::pow(-2*xi_2 - xi_3 + 1, 3))/(xi_3 - 1) ;
}
if( b == 98 ) {
  valDxi1 = -2*(780*xi_2*std::pow(xi_3 - 1, 6) + 365*std::pow(xi_3 - 1, 7) + 165*std::pow(xi_3 - 1, 5)*std::pow(2*xi_2 + xi_3 - 1, 2) - 3740*std::pow(xi_3 - 1, 4)*std::pow(2*xi_2 + xi_3 - 1, 3) + 2145*std::pow(xi_3 - 1, 3)*std::pow(2*xi_2 + xi_3 - 1, 4) + 6006*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 5) - 5005*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 6))/(64*xi_3 - 64) ;
}
if( b == 99 ) {
  valDxi1 = (1.0L/32.0L)*(24*xi_1 + 12*xi_2 + 12*xi_3 - 12)*(-600*xi_3 + 2184*std::pow(2*xi_3 - 1, 5) + 4095*std::pow(2*xi_3 - 1, 4) + 1820*std::pow(2*xi_3 - 1, 3) - 390*std::pow(2*xi_3 - 1, 2) + 283) ;
}
if( b == 100 ) {
  valDxi1 = (1.0L/8.0L)*(7*xi_2 + xi_3 - 1)*(24*xi_1 + 12*xi_2 + 12*xi_3 - 12)*(156*xi_3 + 455*std::pow(2*xi_3 - 1, 4) + 910*std::pow(2*xi_3 - 1, 3) + 546*std::pow(2*xi_3 - 1, 2) - 87) ;
}
if( b == 101 ) {
  valDxi1 = (20*xi_2*std::pow(xi_3 - 1, 2) + 8*std::pow(xi_3 - 1, 3) - 9*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2))*(-24*xi_1 - 12*xi_2 - 12*xi_3 + 12)*(588*xi_3 + 280*std::pow(2*xi_3 - 1, 3) + 525*std::pow(2*xi_3 - 1, 2) - 249)/(8*xi_3 - 8) ;
}
if( b == 102 ) {
  valDxi1 = (1.0L/16.0L)*(90*xi_3 + 30*std::pow(2*xi_3 - 1, 2) - 29)*(24*xi_1 + 12*xi_2 + 12*xi_3 - 12)*(126*xi_2*std::pow(xi_3 - 1, 3) - 225*std::pow(-xi_3 + 1, 2)*std::pow(-2*xi_2 - xi_3 + 1, 2) + 165*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 3) + 68*std::pow(xi_3 - 1, 4))/(xi_3 - 1) ;
}
if( b == 103 ) {
  valDxi1 = (1.0L/32.0L)*(16*xi_3 - 1)*(24*xi_1 + 12*xi_2 + 12*xi_3 - 12)*(200*xi_2*std::pow(xi_3 - 1, 3) + 71*std::pow(xi_3 - 1, 4) + 330*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 2) - 1100*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 3) + 715*std::pow(2*xi_2 + xi_3 - 1, 4)) ;
}
if( b == 104 ) {
  valDxi1 = (-24*xi_1 - 12*xi_2 - 12*xi_3 + 12)*(770*xi_2*std::pow(xi_3 - 1, 5) + 386*std::pow(xi_3 - 1, 6) - 990*std::pow(xi_3 - 1, 4)*std::pow(2*xi_2 + xi_3 - 1, 2) - 1430*std::pow(xi_3 - 1, 3)*std::pow(2*xi_2 + xi_3 - 1, 3) + 5005*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 4) - 3003*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 5))/(64*xi_3 - 64) ;
}
if( b == 105 ) {
  valDxi1 = -(6*std::pow(xi_2 + xi_3 - 1, 3) - 30*(xi_2 + xi_3 - 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2))*(156*xi_3 + 455*std::pow(2*xi_3 - 1, 4) + 910*std::pow(2*xi_3 - 1, 3) + 546*std::pow(2*xi_3 - 1, 2) - 87)/(8*xi_2 + 8*xi_3 - 8) ;
}
if( b == 106 ) {
  valDxi1 = (1.0L/8.0L)*(30*(-xi_2 - xi_3 + 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + 6*std::pow(xi_2 + xi_3 - 1, 3))*(9*xi_2 + xi_3 - 1)*(588*xi_3 + 280*std::pow(2*xi_3 - 1, 3) + 525*std::pow(2*xi_3 - 1, 2) - 249)/(-xi_2 - xi_3 + 1) ;
}
if( b == 107 ) {
  valDxi1 = (1.0L/8.0L)*(30*(-xi_2 - xi_3 + 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + 6*std::pow(xi_2 + xi_3 - 1, 3))*(90*xi_3 + 30*std::pow(2*xi_3 - 1, 2) - 29)*(140*xi_2*std::pow(xi_3 - 1, 2) + 55*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 2) + 51*std::pow(xi_3 - 1, 3))/((xi_3 - 1)*(xi_2 + xi_3 - 1)) ;
}
if( b == 108 ) {
  valDxi1 = -(16*xi_3 - 1)*(6*std::pow(xi_2 + xi_3 - 1, 3) - 30*(xi_2 + xi_3 - 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2))*(198*xi_2*std::pow(xi_3 - 1, 2) + 92*std::pow(xi_3 - 1, 3) - 231*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2) + 143*std::pow(2*xi_2 + xi_3 - 1, 3))/(8*xi_2 + 8*xi_3 - 8) ;
}
if( b == 109 ) {
  valDxi1 = (6*std::pow(xi_2 + xi_3 - 1, 3) - 30*(xi_2 + xi_3 - 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2))*(168*xi_2*std::pow(xi_3 - 1, 3) + 127*std::pow(xi_3 - 1, 4) - 1326*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 2) + 2548*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 3) - 1365*std::pow(2*xi_2 + xi_3 - 1, 4))/(32*xi_2 + 32*xi_3 - 32) ;
}
if( b == 110 ) {
  valDxi1 = (1.0L/32.0L)*(-30*std::pow(xi_2 + xi_3 - 1, 2)*(8*xi_1 + 4*xi_2 + 4*xi_3 - 4) + 280*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3))*(588*xi_3 + 280*std::pow(2*xi_3 - 1, 3) + 525*std::pow(2*xi_3 - 1, 2) - 249) ;
}
if( b == 111 ) {
  valDxi1 = (1.0L/8.0L)*(-30*std::pow(xi_2 + xi_3 - 1, 2)*(8*xi_1 + 4*xi_2 + 4*xi_3 - 4) + 280*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3))*(11*xi_2 + xi_3 - 1)*(90*xi_3 + 30*std::pow(2*xi_3 - 1, 2) - 29) ;
}
if( b == 112 ) {
  valDxi1 = -(16*xi_3 - 1)*(-30*std::pow(xi_2 + xi_3 - 1, 2)*(8*xi_1 + 4*xi_2 + 4*xi_3 - 4) + 280*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3))*(108*xi_2*std::pow(xi_3 - 1, 2) + 37*std::pow(xi_3 - 1, 3) - 39*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2))/(16*xi_3 - 16) ;
}
if( b == 113 ) {
  valDxi1 = (1.0L/64.0L)*(-30*std::pow(xi_2 + xi_3 - 1, 2)*(8*xi_1 + 4*xi_2 + 4*xi_3 - 4) + 280*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3))*(858*xi_2*std::pow(xi_3 - 1, 3) - 819*std::pow(-xi_3 + 1, 2)*std::pow(-2*xi_2 - xi_3 + 1, 2) + 455*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 3) + 372*std::pow(xi_3 - 1, 4))/(xi_3 - 1) ;
}
if( b == 114 ) {
  valDxi1 = (1.0L/8.0L)*(90*xi_3 + 30*std::pow(2*xi_3 - 1, 2) - 29)*(-630*(-xi_2 - xi_3 + 1)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 4) + 30*std::pow(xi_2 + xi_3 - 1, 5) - 420*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 2))/(xi_2 + xi_3 - 1) ;
}
if( b == 115 ) {
  valDxi1 = (1.0L/8.0L)*(16*xi_3 - 1)*(13*xi_2 + xi_3 - 1)*(-630*(-xi_2 - xi_3 + 1)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 4) + 30*std::pow(xi_2 + xi_3 - 1, 5) - 420*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 2))/(xi_2 + xi_3 - 1) ;
}
if( b == 116 ) {
  valDxi1 = -1.0L/32.0L*(308*xi_2*std::pow(xi_3 - 1, 2) + 101*std::pow(xi_3 - 1, 3) - 105*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2))*(30*std::pow(xi_2 + xi_3 - 1, 5) - 420*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + 630*(xi_2 + xi_3 - 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 4))/((xi_3 - 1)*(xi_2 + xi_3 - 1)) ;
}
if( b == 117 ) {
  valDxi1 = (1.0L/16.0L)*(16*xi_3 - 1)*(105*std::pow(xi_2 + xi_3 - 1, 4)*(8*xi_1 + 4*xi_2 + 4*xi_3 - 4) - 2520*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3) + 2772*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 5)) ;
}
if( b == 118 ) {
  valDxi1 = (1.0L/16.0L)*(15*xi_2 + xi_3 - 1)*(105*std::pow(xi_2 + xi_3 - 1, 4)*(8*xi_1 + 4*xi_2 + 4*xi_3 - 4) - 2520*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3) + 2772*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 5)) ;
}
if( b == 119 ) {
  valDxi1 = -(70*std::pow(xi_2 + xi_3 - 1, 7) - 1890*std::pow(xi_2 + xi_3 - 1, 5)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + 6930*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 4) - 6006*(xi_2 + xi_3 - 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 6))/(16*xi_2 + 16*xi_3 - 16) ;
}
if( b == 120 ) {
  valDxi1 = 0 ;
}
if( b == 121 ) {
  valDxi1 = 0 ;
}
if( b == 122 ) {
  valDxi1 = 0 ;
}
if( b == 123 ) {
  valDxi1 = 0 ;
}
if( b == 124 ) {
  valDxi1 = 0 ;
}
if( b == 125 ) {
  valDxi1 = 0 ;
}
if( b == 126 ) {
  valDxi1 = 0 ;
}
if( b == 127 ) {
  valDxi1 = 0 ;
}
if( b == 128 ) {
  valDxi1 = 0 ;
}
if( b == 129 ) {
  valDxi1 = (63.0L/2.0L)*xi_3 + (1989.0L/4.0L)*std::pow(2*xi_3 - 1, 7) + (1547.0L/2.0L)*std::pow(2*xi_3 - 1, 6) - 273.0L/4.0L*std::pow(2*xi_3 - 1, 5) - 4095.0L/8.0L*std::pow(2*xi_3 - 1, 4) - 455.0L/4.0L*std::pow(2*xi_3 - 1, 3) + (273.0L/4.0L)*std::pow(2*xi_3 - 1, 2) - 133.0L/8.0L ;
}
if( b == 130 ) {
  valDxi1 = (1.0L/8.0L)*(5*xi_2 + xi_3 - 1)*(-156*xi_3 + 4641*std::pow(2*xi_3 - 1, 6) + 9282*std::pow(2*xi_3 - 1, 5) + 4095*std::pow(2*xi_3 - 1, 4) - 1820*std::pow(2*xi_3 - 1, 3) - 1365*std::pow(2*xi_3 - 1, 2) + 107) ;
}
if( b == 131 ) {
  valDxi1 = -2*(36*xi_2*std::pow(xi_3 - 1, 2) + 17*std::pow(xi_3 - 1, 3) - 21*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2))*(-210*xi_3 + 1071*std::pow(2*xi_3 - 1, 5) + 2380*std::pow(2*xi_3 - 1, 4) + 1610*std::pow(2*xi_3 - 1, 3) + 210*std::pow(2*xi_3 - 1, 2) + 87)/(16*xi_3 - 16) ;
}
if( b == 132 ) {
  valDxi1 = (1.0L/4.0L)*(2*std::pow(xi_3 - 1, 3) - 21*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2) + 21*std::pow(2*xi_2 + xi_3 - 1, 3))*(600*xi_3 + 765*std::pow(2*xi_3 - 1, 4) + 1700*std::pow(2*xi_3 - 1, 3) + 1230*std::pow(2*xi_3 - 1, 2) - 291) ;
}
if( b == 133 ) {
  valDxi1 = (1.0L/4.0L)*(252*xi_3 + 102*std::pow(2*xi_3 - 1, 3) + 204*std::pow(2*xi_3 - 1, 2) - 103)*(88*xi_2*std::pow(xi_3 - 1, 3) + 41*std::pow(xi_3 - 1, 4) - 18*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 2) - 180*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 3) + 165*std::pow(2*xi_2 + xi_3 - 1, 4)) ;
}
if( b == 134 ) {
  valDxi1 = (1.0L/64.0L)*(476*xi_3 + 153*std::pow(2*xi_3 - 1, 2) - 149)*(-90*xi_2*std::pow(xi_3 - 1, 5) - 70*std::pow(-xi_3 + 1, 6) - 1485*std::pow(-xi_3 + 1, 2)*std::pow(-2*xi_2 - xi_3 + 1, 4) + 1287*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 5) + 630*std::pow(xi_3 - 1, 4)*std::pow(-2*xi_2 - xi_3 + 1, 2) + 330*std::pow(xi_3 - 1, 3)*std::pow(-2*xi_2 - xi_3 + 1, 3))/(xi_3 - 1) ;
}
if( b == 135 ) {
  valDxi1 = -2*(18*xi_3 - 1)*(780*xi_2*std::pow(xi_3 - 1, 6) + 365*std::pow(xi_3 - 1, 7) + 165*std::pow(xi_3 - 1, 5)*std::pow(2*xi_2 + xi_3 - 1, 2) - 3740*std::pow(xi_3 - 1, 4)*std::pow(2*xi_2 + xi_3 - 1, 3) + 2145*std::pow(xi_3 - 1, 3)*std::pow(2*xi_2 + xi_3 - 1, 4) + 6006*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 5) - 5005*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 6))/(64*xi_3 - 64) ;
}
if( b == 136 ) {
  valDxi1 = (1.0L/16.0L)*(154*xi_2*std::pow(xi_3 - 1, 7) - 924*std::pow(-xi_3 + 1, 6)*std::pow(-2*xi_2 - xi_3 + 1, 2) - 6006*std::pow(-xi_3 + 1, 2)*std::pow(-2*xi_2 - xi_3 + 1, 6) + 4862*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 7) + 98*std::pow(xi_3 - 1, 8) + 5005*std::pow(xi_3 - 1, 4)*std::pow(-2*xi_2 - xi_3 + 1, 4) + 3003*std::pow(xi_3 - 1, 3)*std::pow(-2*xi_2 - xi_3 + 1, 5))/(xi_3 - 1) ;
}
if( b == 137 ) {
  valDxi1 = (1.0L/32.0L)*(24*xi_1 + 12*xi_2 + 12*xi_3 - 12)*(-156*xi_3 + 4641*std::pow(2*xi_3 - 1, 6) + 9282*std::pow(2*xi_3 - 1, 5) + 4095*std::pow(2*xi_3 - 1, 4) - 1820*std::pow(2*xi_3 - 1, 3) - 1365*std::pow(2*xi_3 - 1, 2) + 107) ;
}
if( b == 138 ) {
  valDxi1 = (1.0L/8.0L)*(7*xi_2 + xi_3 - 1)*(24*xi_1 + 12*xi_2 + 12*xi_3 - 12)*(-210*xi_3 + 1071*std::pow(2*xi_3 - 1, 5) + 2380*std::pow(2*xi_3 - 1, 4) + 1610*std::pow(2*xi_3 - 1, 3) + 210*std::pow(2*xi_3 - 1, 2) + 87) ;
}
if( b == 139 ) {
  valDxi1 = (20*xi_2*std::pow(xi_3 - 1, 2) + 8*std::pow(xi_3 - 1, 3) - 9*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2))*(-24*xi_1 - 12*xi_2 - 12*xi_3 + 12)*(600*xi_3 + 765*std::pow(2*xi_3 - 1, 4) + 1700*std::pow(2*xi_3 - 1, 3) + 1230*std::pow(2*xi_3 - 1, 2) - 291)/(8*xi_3 - 8) ;
}
if( b == 140 ) {
  valDxi1 = (1.0L/16.0L)*(24*xi_1 + 12*xi_2 + 12*xi_3 - 12)*(252*xi_3 + 102*std::pow(2*xi_3 - 1, 3) + 204*std::pow(2*xi_3 - 1, 2) - 103)*(126*xi_2*std::pow(xi_3 - 1, 3) - 225*std::pow(-xi_3 + 1, 2)*std::pow(-2*xi_2 - xi_3 + 1, 2) + 165*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 3) + 68*std::pow(xi_3 - 1, 4))/(xi_3 - 1) ;
}
if( b == 141 ) {
  valDxi1 = (1.0L/128.0L)*(476*xi_3 + 153*std::pow(2*xi_3 - 1, 2) - 149)*(24*xi_1 + 12*xi_2 + 12*xi_3 - 12)*(200*xi_2*std::pow(xi_3 - 1, 3) + 71*std::pow(xi_3 - 1, 4) + 330*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 2) - 1100*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 3) + 715*std::pow(2*xi_2 + xi_3 - 1, 4)) ;
}
if( b == 142 ) {
  valDxi1 = (18*xi_3 - 1)*(-24*xi_1 - 12*xi_2 - 12*xi_3 + 12)*(770*xi_2*std::pow(xi_3 - 1, 5) + 386*std::pow(xi_3 - 1, 6) - 990*std::pow(xi_3 - 1, 4)*std::pow(2*xi_2 + xi_3 - 1, 2) - 1430*std::pow(xi_3 - 1, 3)*std::pow(2*xi_2 + xi_3 - 1, 3) + 5005*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 4) - 3003*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 5))/(64*xi_3 - 64) ;
}
if( b == 143 ) {
  valDxi1 = (1.0L/32.0L)*(-24*xi_1 - 12*xi_2 - 12*xi_3 + 12)*(96*xi_2*std::pow(xi_3 - 1, 5) + 23*std::pow(xi_3 - 1, 6) + 780*std::pow(xi_3 - 1, 4)*std::pow(2*xi_2 + xi_3 - 1, 2) - 1820*std::pow(xi_3 - 1, 3)*std::pow(2*xi_2 + xi_3 - 1, 3) - 1365*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 4) + 5460*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 5) - 3094*std::pow(2*xi_2 + xi_3 - 1, 6)) ;
}
if( b == 144 ) {
  valDxi1 = -(6*std::pow(xi_2 + xi_3 - 1, 3) - 30*(xi_2 + xi_3 - 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2))*(-210*xi_3 + 1071*std::pow(2*xi_3 - 1, 5) + 2380*std::pow(2*xi_3 - 1, 4) + 1610*std::pow(2*xi_3 - 1, 3) + 210*std::pow(2*xi_3 - 1, 2) + 87)/(8*xi_2 + 8*xi_3 - 8) ;
}
if( b == 145 ) {
  valDxi1 = (1.0L/8.0L)*(30*(-xi_2 - xi_3 + 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + 6*std::pow(xi_2 + xi_3 - 1, 3))*(9*xi_2 + xi_3 - 1)*(600*xi_3 + 765*std::pow(2*xi_3 - 1, 4) + 1700*std::pow(2*xi_3 - 1, 3) + 1230*std::pow(2*xi_3 - 1, 2) - 291)/(-xi_2 - xi_3 + 1) ;
}
if( b == 146 ) {
  valDxi1 = (1.0L/8.0L)*(30*(-xi_2 - xi_3 + 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + 6*std::pow(xi_2 + xi_3 - 1, 3))*(140*xi_2*std::pow(xi_3 - 1, 2) + 55*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 2) + 51*std::pow(xi_3 - 1, 3))*(252*xi_3 + 102*std::pow(2*xi_3 - 1, 3) + 204*std::pow(2*xi_3 - 1, 2) - 103)/((xi_3 - 1)*(xi_2 + xi_3 - 1)) ;
}
if( b == 147 ) {
  valDxi1 = -(6*std::pow(xi_2 + xi_3 - 1, 3) - 30*(xi_2 + xi_3 - 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2))*(476*xi_3 + 153*std::pow(2*xi_3 - 1, 2) - 149)*(198*xi_2*std::pow(xi_3 - 1, 2) + 92*std::pow(xi_3 - 1, 3) - 231*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2) + 143*std::pow(2*xi_2 + xi_3 - 1, 3))/(32*xi_2 + 32*xi_3 - 32) ;
}
if( b == 148 ) {
  valDxi1 = (1.0L/32.0L)*(18*xi_3 - 1)*(30*(-xi_2 - xi_3 + 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + 6*std::pow(xi_2 + xi_3 - 1, 3))*(-168*xi_2*std::pow(xi_3 - 1, 3) - 127*std::pow(xi_3 - 1, 4) + 1326*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 2) - 2548*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 3) + 1365*std::pow(2*xi_2 + xi_3 - 1, 4))/(-xi_2 - xi_3 + 1) ;
}
if( b == 149 ) {
  valDxi1 = (1.0L/16.0L)*(30*(-xi_2 - xi_3 + 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + 6*std::pow(xi_2 + xi_3 - 1, 3))*(390*xi_2*std::pow(xi_3 - 1, 5) + 1820*std::pow(-xi_3 + 1, 3)*std::pow(2*xi_2 + xi_3 - 1, 3) + 1547*(-xi_3 + 1)*std::pow(2*xi_2 + xi_3 - 1, 5) + 174*std::pow(xi_3 - 1, 6) + 3185*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 4))/((xi_3 - 1)*(xi_2 + xi_3 - 1)) ;
}
if( b == 150 ) {
  valDxi1 = (1.0L/32.0L)*(-30*std::pow(xi_2 + xi_3 - 1, 2)*(8*xi_1 + 4*xi_2 + 4*xi_3 - 4) + 280*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3))*(600*xi_3 + 765*std::pow(2*xi_3 - 1, 4) + 1700*std::pow(2*xi_3 - 1, 3) + 1230*std::pow(2*xi_3 - 1, 2) - 291) ;
}
if( b == 151 ) {
  valDxi1 = (1.0L/8.0L)*(-30*std::pow(xi_2 + xi_3 - 1, 2)*(8*xi_1 + 4*xi_2 + 4*xi_3 - 4) + 280*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3))*(11*xi_2 + xi_3 - 1)*(252*xi_3 + 102*std::pow(2*xi_3 - 1, 3) + 204*std::pow(2*xi_3 - 1, 2) - 103) ;
}
if( b == 152 ) {
  valDxi1 = -(-30*std::pow(xi_2 + xi_3 - 1, 2)*(8*xi_1 + 4*xi_2 + 4*xi_3 - 4) + 280*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3))*(476*xi_3 + 153*std::pow(2*xi_3 - 1, 2) - 149)*(108*xi_2*std::pow(xi_3 - 1, 2) + 37*std::pow(xi_3 - 1, 3) - 39*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2))/(64*xi_3 - 64) ;
}
if( b == 153 ) {
  valDxi1 = (1.0L/64.0L)*(18*xi_3 - 1)*(-30*std::pow(xi_2 + xi_3 - 1, 2)*(8*xi_1 + 4*xi_2 + 4*xi_3 - 4) + 280*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3))*(858*xi_2*std::pow(xi_3 - 1, 3) - 819*std::pow(-xi_3 + 1, 2)*std::pow(-2*xi_2 - xi_3 + 1, 2) + 455*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 3) + 372*std::pow(xi_3 - 1, 4))/(xi_3 - 1) ;
}
if( b == 154 ) {
  valDxi1 = (1.0L/32.0L)*(-30*std::pow(xi_2 + xi_3 - 1, 2)*(8*xi_1 + 4*xi_2 + 4*xi_3 - 4) + 280*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3))*(-336*xi_2*std::pow(xi_3 - 1, 3) - 171*std::pow(xi_3 - 1, 4) + 840*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 2) - 1260*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 3) + 595*std::pow(2*xi_2 + xi_3 - 1, 4)) ;
}
if( b == 155 ) {
  valDxi1 = (1.0L/8.0L)*(-630*(-xi_2 - xi_3 + 1)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 4) + 30*std::pow(xi_2 + xi_3 - 1, 5) - 420*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 2))*(252*xi_3 + 102*std::pow(2*xi_3 - 1, 3) + 204*std::pow(2*xi_3 - 1, 2) - 103)/(xi_2 + xi_3 - 1) ;
}
if( b == 156 ) {
  valDxi1 = (1.0L/32.0L)*(13*xi_2 + xi_3 - 1)*(476*xi_3 + 153*std::pow(2*xi_3 - 1, 2) - 149)*(-630*(-xi_2 - xi_3 + 1)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 4) + 30*std::pow(xi_2 + xi_3 - 1, 5) - 420*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 2))/(xi_2 + xi_3 - 1) ;
}
if( b == 157 ) {
  valDxi1 = -1.0L/32.0L*(18*xi_3 - 1)*(308*xi_2*std::pow(xi_3 - 1, 2) + 101*std::pow(xi_3 - 1, 3) - 105*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2))*(30*std::pow(xi_2 + xi_3 - 1, 5) - 420*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + 630*(xi_2 + xi_3 - 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 4))/((xi_3 - 1)*(xi_2 + xi_3 - 1)) ;
}
if( b == 158 ) {
  valDxi1 = (1.0L/8.0L)*(-630*(-xi_2 - xi_3 + 1)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 4) + 30*std::pow(xi_2 + xi_3 - 1, 5) - 420*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 2))*(195*xi_2*std::pow(xi_3 - 1, 2) + 81*std::pow(xi_3 - 1, 3) - 165*(xi_3 - 1)*std::pow(-2*xi_2 - xi_3 + 1, 2) - 85*std::pow(-2*xi_2 - xi_3 + 1, 3))/(xi_2 + xi_3 - 1) ;
}
if( b == 159 ) {
  valDxi1 = (1.0L/64.0L)*(476*xi_3 + 153*std::pow(2*xi_3 - 1, 2) - 149)*(105*std::pow(xi_2 + xi_3 - 1, 4)*(8*xi_1 + 4*xi_2 + 4*xi_3 - 4) - 2520*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3) + 2772*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 5)) ;
}
if( b == 160 ) {
  valDxi1 = (1.0L/16.0L)*(18*xi_3 - 1)*(15*xi_2 + xi_3 - 1)*(105*std::pow(xi_2 + xi_3 - 1, 4)*(8*xi_1 + 4*xi_2 + 4*xi_3 - 4) - 2520*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3) + 2772*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 5)) ;
}
if( b == 161 ) {
  valDxi1 = (104*xi_2*std::pow(xi_3 - 1, 2) + 33*std::pow(xi_3 - 1, 3) - 34*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2))*(-105*std::pow(xi_2 + xi_3 - 1, 4)*(8*xi_1 + 4*xi_2 + 4*xi_3 - 4) + 2520*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3) - 2772*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 5))/(16*xi_3 - 16) ;
}
if( b == 162 ) {
  valDxi1 = -(18*xi_3 - 1)*(70*std::pow(xi_2 + xi_3 - 1, 7) - 1890*std::pow(xi_2 + xi_3 - 1, 5)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + 6930*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 4) - 6006*(xi_2 + xi_3 - 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 6))/(16*xi_2 + 16*xi_3 - 16) ;
}
if( b == 163 ) {
  valDxi1 = -(17*xi_2 + xi_3 - 1)*(70*std::pow(xi_2 + xi_3 - 1, 7) - 1890*std::pow(xi_2 + xi_3 - 1, 5)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + 6930*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 4) - 6006*(xi_2 + xi_3 - 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 6))/(16*xi_2 + 16*xi_3 - 16) ;
}
if( b == 164 ) {
  valDxi1 = -315.0L/32.0L*std::pow(xi_2 + xi_3 - 1, 6)*(8*xi_1 + 4*xi_2 + 4*xi_3 - 4) + (3465.0L/8.0L)*std::pow(xi_2 + xi_3 - 1, 4)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3) - 9009.0L/8.0L*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 5) + (6435.0L/8.0L)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 7) ;
}
if( b == 165 ) {
  valDxi1 = 0 ;
}
if( b == 166 ) {
  valDxi1 = 0 ;
}
if( b == 167 ) {
  valDxi1 = 0 ;
}
if( b == 168 ) {
  valDxi1 = 0 ;
}
if( b == 169 ) {
  valDxi1 = 0 ;
}
if( b == 170 ) {
  valDxi1 = 0 ;
}
if( b == 171 ) {
  valDxi1 = 0 ;
}
if( b == 172 ) {
  valDxi1 = 0 ;
}
if( b == 173 ) {
  valDxi1 = 0 ;
}
if( b == 174 ) {
  valDxi1 = 0 ;
}
if( b == 175 ) {
  valDxi1 = -91.0L/4.0L*xi_3 + (62985.0L/64.0L)*std::pow(2*xi_3 - 1, 8) + (12597.0L/8.0L)*std::pow(2*xi_3 - 1, 7) - 4641.0L/16.0L*std::pow(2*xi_3 - 1, 6) - 10829.0L/8.0L*std::pow(2*xi_3 - 1, 5) - 8645.0L/32.0L*std::pow(2*xi_3 - 1, 4) + (2275.0L/8.0L)*std::pow(2*xi_3 - 1, 3) + (1183.0L/16.0L)*std::pow(2*xi_3 - 1, 2) + 609.0L/64.0L ;
}
if( b == 176 ) {
  valDxi1 = (1.0L/8.0L)*(5*xi_2 + xi_3 - 1)*(616*xi_3 + 9690*std::pow(2*xi_3 - 1, 7) + 20349*std::pow(2*xi_3 - 1, 6) + 8568*std::pow(2*xi_3 - 1, 5) - 6545*std::pow(2*xi_3 - 1, 4) - 4830*std::pow(2*xi_3 - 1, 3) - 105*std::pow(2*xi_3 - 1, 2) - 287) ;
}
if( b == 177 ) {
  valDxi1 = -2*(36*xi_2*std::pow(xi_3 - 1, 2) + 17*std::pow(xi_3 - 1, 3) - 21*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2))*(-600*xi_3 + 4845*std::pow(2*xi_3 - 1, 6) + 11628*std::pow(2*xi_3 - 1, 5) + 8415*std::pow(2*xi_3 - 1, 4) + 680*std::pow(2*xi_3 - 1, 3) - 1245*std::pow(2*xi_3 - 1, 2) + 301)/(32*xi_3 - 32) ;
}
if( b == 178 ) {
  valDxi1 = (1.0L/4.0L)*(2*std::pow(xi_3 - 1, 3) - 21*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2) + 21*std::pow(2*xi_2 + xi_3 - 1, 3))*(-20*xi_3 + 1938*std::pow(2*xi_3 - 1, 5) + 4845*std::pow(2*xi_3 - 1, 4) + 4080*std::pow(2*xi_3 - 1, 3) + 1190*std::pow(2*xi_3 - 1, 2) - 21) ;
}
if( b == 179 ) {
  valDxi1 = (1.0L/64.0L)*(5848*xi_3 + 4845*std::pow(2*xi_3 - 1, 4) + 11628*std::pow(2*xi_3 - 1, 3) + 9486*std::pow(2*xi_3 - 1, 2) - 2687)*(88*xi_2*std::pow(xi_3 - 1, 3) + 41*std::pow(xi_3 - 1, 4) - 18*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 2) - 180*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 3) + 165*std::pow(2*xi_2 + xi_3 - 1, 4)) ;
}
if( b == 180 ) {
  valDxi1 = (1.0L/64.0L)*(1584*xi_3 + 570*std::pow(2*xi_3 - 1, 3) + 1197*std::pow(2*xi_3 - 1, 2) - 631)*(-90*xi_2*std::pow(xi_3 - 1, 5) - 70*std::pow(-xi_3 + 1, 6) - 1485*std::pow(-xi_3 + 1, 2)*std::pow(-2*xi_2 - xi_3 + 1, 4) + 1287*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 5) + 630*std::pow(xi_3 - 1, 4)*std::pow(-2*xi_2 - xi_3 + 1, 2) + 330*std::pow(xi_3 - 1, 3)*std::pow(-2*xi_2 - xi_3 + 1, 3))/(xi_3 - 1) ;
}
if( b == 181 ) {
  valDxi1 = -2*(304*xi_3 + 95*std::pow(2*xi_3 - 1, 2) - 93)*(780*xi_2*std::pow(xi_3 - 1, 6) + 365*std::pow(xi_3 - 1, 7) + 165*std::pow(xi_3 - 1, 5)*std::pow(2*xi_2 + xi_3 - 1, 2) - 3740*std::pow(xi_3 - 1, 4)*std::pow(2*xi_2 + xi_3 - 1, 3) + 2145*std::pow(xi_3 - 1, 3)*std::pow(2*xi_2 + xi_3 - 1, 4) + 6006*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 5) - 5005*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 6))/(128*xi_3 - 128) ;
}
if( b == 182 ) {
  valDxi1 = (1.0L/16.0L)*(20*xi_3 - 1)*(154*xi_2*std::pow(xi_3 - 1, 7) - 924*std::pow(-xi_3 + 1, 6)*std::pow(-2*xi_2 - xi_3 + 1, 2) - 6006*std::pow(-xi_3 + 1, 2)*std::pow(-2*xi_2 - xi_3 + 1, 6) + 4862*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 7) + 98*std::pow(xi_3 - 1, 8) + 5005*std::pow(xi_3 - 1, 4)*std::pow(-2*xi_2 - xi_3 + 1, 4) + 3003*std::pow(xi_3 - 1, 3)*std::pow(-2*xi_2 - xi_3 + 1, 5))/(xi_3 - 1) ;
}
if( b == 183 ) {
  valDxi1 = (105.0L/4.0L)*xi_2*std::pow(xi_3 - 1, 7) + (1729.0L/8.0L)*std::pow(-xi_3 + 1, 5)*std::pow(2*xi_2 + xi_3 - 1, 3) - 7735.0L/16.0L*std::pow(-xi_3 + 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 6) + (5967.0L/8.0L)*(-xi_3 + 1)*std::pow(2*xi_2 + xi_3 - 1, 7) + (791.0L/64.0L)*std::pow(xi_3 - 1, 8) + (273.0L/16.0L)*std::pow(xi_3 - 1, 6)*std::pow(2*xi_2 + xi_3 - 1, 2) + (1365.0L/32.0L)*std::pow(xi_3 - 1, 4)*std::pow(2*xi_2 + xi_3 - 1, 4) + (6279.0L/8.0L)*std::pow(xi_3 - 1, 3)*std::pow(2*xi_2 + xi_3 - 1, 5) + (37791.0L/64.0L)*std::pow(2*xi_2 + xi_3 - 1, 8) ;
}
if( b == 184 ) {
  valDxi1 = (1.0L/32.0L)*(24*xi_1 + 12*xi_2 + 12*xi_3 - 12)*(616*xi_3 + 9690*std::pow(2*xi_3 - 1, 7) + 20349*std::pow(2*xi_3 - 1, 6) + 8568*std::pow(2*xi_3 - 1, 5) - 6545*std::pow(2*xi_3 - 1, 4) - 4830*std::pow(2*xi_3 - 1, 3) - 105*std::pow(2*xi_3 - 1, 2) - 287) ;
}
if( b == 185 ) {
  valDxi1 = (1.0L/16.0L)*(7*xi_2 + xi_3 - 1)*(24*xi_1 + 12*xi_2 + 12*xi_3 - 12)*(-600*xi_3 + 4845*std::pow(2*xi_3 - 1, 6) + 11628*std::pow(2*xi_3 - 1, 5) + 8415*std::pow(2*xi_3 - 1, 4) + 680*std::pow(2*xi_3 - 1, 3) - 1245*std::pow(2*xi_3 - 1, 2) + 301) ;
}
if( b == 186 ) {
  valDxi1 = (20*xi_2*std::pow(xi_3 - 1, 2) + 8*std::pow(xi_3 - 1, 3) - 9*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2))*(-24*xi_1 - 12*xi_2 - 12*xi_3 + 12)*(-20*xi_3 + 1938*std::pow(2*xi_3 - 1, 5) + 4845*std::pow(2*xi_3 - 1, 4) + 4080*std::pow(2*xi_3 - 1, 3) + 1190*std::pow(2*xi_3 - 1, 2) - 21)/(8*xi_3 - 8) ;
}
if( b == 187 ) {
  valDxi1 = (1.0L/256.0L)*(24*xi_1 + 12*xi_2 + 12*xi_3 - 12)*(126*xi_2*std::pow(xi_3 - 1, 3) - 225*std::pow(-xi_3 + 1, 2)*std::pow(-2*xi_2 - xi_3 + 1, 2) + 165*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 3) + 68*std::pow(xi_3 - 1, 4))*(5848*xi_3 + 4845*std::pow(2*xi_3 - 1, 4) + 11628*std::pow(2*xi_3 - 1, 3) + 9486*std::pow(2*xi_3 - 1, 2) - 2687)/(xi_3 - 1) ;
}
if( b == 188 ) {
  valDxi1 = (1.0L/128.0L)*(24*xi_1 + 12*xi_2 + 12*xi_3 - 12)*(1584*xi_3 + 570*std::pow(2*xi_3 - 1, 3) + 1197*std::pow(2*xi_3 - 1, 2) - 631)*(200*xi_2*std::pow(xi_3 - 1, 3) + 71*std::pow(xi_3 - 1, 4) + 330*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 2) - 1100*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 3) + 715*std::pow(2*xi_2 + xi_3 - 1, 4)) ;
}
if( b == 189 ) {
  valDxi1 = (304*xi_3 + 95*std::pow(2*xi_3 - 1, 2) - 93)*(-24*xi_1 - 12*xi_2 - 12*xi_3 + 12)*(770*xi_2*std::pow(xi_3 - 1, 5) + 386*std::pow(xi_3 - 1, 6) - 990*std::pow(xi_3 - 1, 4)*std::pow(2*xi_2 + xi_3 - 1, 2) - 1430*std::pow(xi_3 - 1, 3)*std::pow(2*xi_2 + xi_3 - 1, 3) + 5005*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 4) - 3003*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 5))/(128*xi_3 - 128) ;
}
if( b == 190 ) {
  valDxi1 = (1.0L/32.0L)*(20*xi_3 - 1)*(-24*xi_1 - 12*xi_2 - 12*xi_3 + 12)*(96*xi_2*std::pow(xi_3 - 1, 5) + 23*std::pow(xi_3 - 1, 6) + 780*std::pow(xi_3 - 1, 4)*std::pow(2*xi_2 + xi_3 - 1, 2) - 1820*std::pow(xi_3 - 1, 3)*std::pow(2*xi_2 + xi_3 - 1, 3) - 1365*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 4) + 5460*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 5) - 3094*std::pow(2*xi_2 + xi_3 - 1, 6)) ;
}
if( b == 191 ) {
  valDxi1 = (1.0L/64.0L)*(24*xi_1 + 12*xi_2 + 12*xi_3 - 12)*(910*xi_2*std::pow(xi_3 - 1, 7) - 819*std::pow(-xi_3 + 1, 6)*std::pow(-2*xi_2 - xi_3 + 1, 2) + 4641*std::pow(-xi_3 + 1, 3)*std::pow(-2*xi_2 - xi_3 + 1, 5) - 23205*std::pow(-xi_3 + 1, 2)*std::pow(-2*xi_2 - xi_3 + 1, 6) + 12597*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 7) + 448*std::pow(xi_3 - 1, 8) + 5005*std::pow(xi_3 - 1, 5)*std::pow(-2*xi_2 - xi_3 + 1, 3) + 11375*std::pow(xi_3 - 1, 4)*std::pow(-2*xi_2 - xi_3 + 1, 4))/(xi_3 - 1) ;
}
if( b == 192 ) {
  valDxi1 = -(6*std::pow(xi_2 + xi_3 - 1, 3) - 30*(xi_2 + xi_3 - 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2))*(-600*xi_3 + 4845*std::pow(2*xi_3 - 1, 6) + 11628*std::pow(2*xi_3 - 1, 5) + 8415*std::pow(2*xi_3 - 1, 4) + 680*std::pow(2*xi_3 - 1, 3) - 1245*std::pow(2*xi_3 - 1, 2) + 301)/(16*xi_2 + 16*xi_3 - 16) ;
}
if( b == 193 ) {
  valDxi1 = (1.0L/8.0L)*(30*(-xi_2 - xi_3 + 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + 6*std::pow(xi_2 + xi_3 - 1, 3))*(9*xi_2 + xi_3 - 1)*(-20*xi_3 + 1938*std::pow(2*xi_3 - 1, 5) + 4845*std::pow(2*xi_3 - 1, 4) + 4080*std::pow(2*xi_3 - 1, 3) + 1190*std::pow(2*xi_3 - 1, 2) - 21)/(-xi_2 - xi_3 + 1) ;
}
if( b == 194 ) {
  valDxi1 = (1.0L/128.0L)*(30*(-xi_2 - xi_3 + 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + 6*std::pow(xi_2 + xi_3 - 1, 3))*(140*xi_2*std::pow(xi_3 - 1, 2) + 55*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 2) + 51*std::pow(xi_3 - 1, 3))*(5848*xi_3 + 4845*std::pow(2*xi_3 - 1, 4) + 11628*std::pow(2*xi_3 - 1, 3) + 9486*std::pow(2*xi_3 - 1, 2) - 2687)/((xi_3 - 1)*(xi_2 + xi_3 - 1)) ;
}
if( b == 195 ) {
  valDxi1 = -(6*std::pow(xi_2 + xi_3 - 1, 3) - 30*(xi_2 + xi_3 - 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2))*(1584*xi_3 + 570*std::pow(2*xi_3 - 1, 3) + 1197*std::pow(2*xi_3 - 1, 2) - 631)*(198*xi_2*std::pow(xi_3 - 1, 2) + 92*std::pow(xi_3 - 1, 3) - 231*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2) + 143*std::pow(2*xi_2 + xi_3 - 1, 3))/(32*xi_2 + 32*xi_3 - 32) ;
}
if( b == 196 ) {
  valDxi1 = (6*std::pow(xi_2 + xi_3 - 1, 3) - 30*(xi_2 + xi_3 - 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2))*(304*xi_3 + 95*std::pow(2*xi_3 - 1, 2) - 93)*(168*xi_2*std::pow(xi_3 - 1, 3) + 127*std::pow(xi_3 - 1, 4) - 1326*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 2) + 2548*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 3) - 1365*std::pow(2*xi_2 + xi_3 - 1, 4))/(64*xi_2 + 64*xi_3 - 64) ;
}
if( b == 197 ) {
  valDxi1 = (1.0L/16.0L)*(20*xi_3 - 1)*(30*(-xi_2 - xi_3 + 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + 6*std::pow(xi_2 + xi_3 - 1, 3))*(390*xi_2*std::pow(xi_3 - 1, 5) + 1820*std::pow(-xi_3 + 1, 3)*std::pow(2*xi_2 + xi_3 - 1, 3) + 1547*(-xi_3 + 1)*std::pow(2*xi_2 + xi_3 - 1, 5) + 174*std::pow(xi_3 - 1, 6) + 3185*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 4))/((xi_3 - 1)*(xi_2 + xi_3 - 1)) ;
}
if( b == 198 ) {
  valDxi1 = -(6*std::pow(xi_2 + xi_3 - 1, 3) - 30*(xi_2 + xi_3 - 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2))*(588*xi_2*std::pow(xi_3 - 1, 5) + 317*std::pow(xi_3 - 1, 6) - 1995*std::pow(xi_3 - 1, 4)*std::pow(2*xi_2 + xi_3 - 1, 2) + 980*std::pow(xi_3 - 1, 3)*std::pow(2*xi_2 + xi_3 - 1, 3) + 8925*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 4) - 14994*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 5) + 6783*std::pow(2*xi_2 + xi_3 - 1, 6))/(32*xi_2 + 32*xi_3 - 32) ;
}
if( b == 199 ) {
  valDxi1 = (1.0L/32.0L)*(-30*std::pow(xi_2 + xi_3 - 1, 2)*(8*xi_1 + 4*xi_2 + 4*xi_3 - 4) + 280*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3))*(-20*xi_3 + 1938*std::pow(2*xi_3 - 1, 5) + 4845*std::pow(2*xi_3 - 1, 4) + 4080*std::pow(2*xi_3 - 1, 3) + 1190*std::pow(2*xi_3 - 1, 2) - 21) ;
}
if( b == 200 ) {
  valDxi1 = (1.0L/128.0L)*(-30*std::pow(xi_2 + xi_3 - 1, 2)*(8*xi_1 + 4*xi_2 + 4*xi_3 - 4) + 280*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3))*(11*xi_2 + xi_3 - 1)*(5848*xi_3 + 4845*std::pow(2*xi_3 - 1, 4) + 11628*std::pow(2*xi_3 - 1, 3) + 9486*std::pow(2*xi_3 - 1, 2) - 2687) ;
}
if( b == 201 ) {
  valDxi1 = -(-30*std::pow(xi_2 + xi_3 - 1, 2)*(8*xi_1 + 4*xi_2 + 4*xi_3 - 4) + 280*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3))*(108*xi_2*std::pow(xi_3 - 1, 2) + 37*std::pow(xi_3 - 1, 3) - 39*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2))*(1584*xi_3 + 570*std::pow(2*xi_3 - 1, 3) + 1197*std::pow(2*xi_3 - 1, 2) - 631)/(64*xi_3 - 64) ;
}
if( b == 202 ) {
  valDxi1 = (1.0L/128.0L)*(-30*std::pow(xi_2 + xi_3 - 1, 2)*(8*xi_1 + 4*xi_2 + 4*xi_3 - 4) + 280*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3))*(304*xi_3 + 95*std::pow(2*xi_3 - 1, 2) - 93)*(858*xi_2*std::pow(xi_3 - 1, 3) - 819*std::pow(-xi_3 + 1, 2)*std::pow(-2*xi_2 - xi_3 + 1, 2) + 455*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 3) + 372*std::pow(xi_3 - 1, 4))/(xi_3 - 1) ;
}
if( b == 203 ) {
  valDxi1 = (1.0L/32.0L)*(20*xi_3 - 1)*(-30*std::pow(xi_2 + xi_3 - 1, 2)*(8*xi_1 + 4*xi_2 + 4*xi_3 - 4) + 280*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3))*(-336*xi_2*std::pow(xi_3 - 1, 3) - 171*std::pow(xi_3 - 1, 4) + 840*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 2) - 1260*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 3) + 595*std::pow(2*xi_2 + xi_3 - 1, 4)) ;
}
if( b == 204 ) {
  valDxi1 = -(-30*std::pow(xi_2 + xi_3 - 1, 2)*(8*xi_1 + 4*xi_2 + 4*xi_3 - 4) + 280*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3))*(330*xi_2*std::pow(xi_3 - 1, 5) + 114*std::pow(xi_3 - 1, 6) + 1170*std::pow(xi_3 - 1, 4)*std::pow(2*xi_2 + xi_3 - 1, 2) - 5270*std::pow(xi_3 - 1, 3)*std::pow(2*xi_2 + xi_3 - 1, 3) + 6885*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 4) - 2907*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 5))/(64*xi_3 - 64) ;
}
if( b == 205 ) {
  valDxi1 = (1.0L/128.0L)*(-630*(-xi_2 - xi_3 + 1)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 4) + 30*std::pow(xi_2 + xi_3 - 1, 5) - 420*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 2))*(5848*xi_3 + 4845*std::pow(2*xi_3 - 1, 4) + 11628*std::pow(2*xi_3 - 1, 3) + 9486*std::pow(2*xi_3 - 1, 2) - 2687)/(xi_2 + xi_3 - 1) ;
}
if( b == 206 ) {
  valDxi1 = (1.0L/32.0L)*(13*xi_2 + xi_3 - 1)*(-630*(-xi_2 - xi_3 + 1)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 4) + 30*std::pow(xi_2 + xi_3 - 1, 5) - 420*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 2))*(1584*xi_3 + 570*std::pow(2*xi_3 - 1, 3) + 1197*std::pow(2*xi_3 - 1, 2) - 631)/(xi_2 + xi_3 - 1) ;
}
if( b == 207 ) {
  valDxi1 = -1.0L/64.0L*(304*xi_3 + 95*std::pow(2*xi_3 - 1, 2) - 93)*(308*xi_2*std::pow(xi_3 - 1, 2) + 101*std::pow(xi_3 - 1, 3) - 105*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2))*(30*std::pow(xi_2 + xi_3 - 1, 5) - 420*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + 630*(xi_2 + xi_3 - 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 4))/((xi_3 - 1)*(xi_2 + xi_3 - 1)) ;
}
if( b == 208 ) {
  valDxi1 = (1.0L/8.0L)*(20*xi_3 - 1)*(-630*(-xi_2 - xi_3 + 1)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 4) + 30*std::pow(xi_2 + xi_3 - 1, 5) - 420*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 2))*(195*xi_2*std::pow(xi_3 - 1, 2) + 81*std::pow(xi_3 - 1, 3) - 165*(xi_3 - 1)*std::pow(-2*xi_2 - xi_3 + 1, 2) - 85*std::pow(-2*xi_2 - xi_3 + 1, 3))/(xi_2 + xi_3 - 1) ;
}
if( b == 209 ) {
  valDxi1 = (1.0L/32.0L)*(-630*(-xi_2 - xi_3 + 1)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 4) + 30*std::pow(xi_2 + xi_3 - 1, 5) - 420*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 2))*(-968*xi_2*std::pow(xi_3 - 1, 3) - 455*std::pow(xi_3 - 1, 4) + 1734*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 2) - 2244*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 3) + 969*std::pow(2*xi_2 + xi_3 - 1, 4))/(xi_2 + xi_3 - 1) ;
}
if( b == 210 ) {
  valDxi1 = (1.0L/64.0L)*(105*std::pow(xi_2 + xi_3 - 1, 4)*(8*xi_1 + 4*xi_2 + 4*xi_3 - 4) - 2520*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3) + 2772*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 5))*(1584*xi_3 + 570*std::pow(2*xi_3 - 1, 3) + 1197*std::pow(2*xi_3 - 1, 2) - 631) ;
}
if( b == 211 ) {
  valDxi1 = (1.0L/32.0L)*(15*xi_2 + xi_3 - 1)*(304*xi_3 + 95*std::pow(2*xi_3 - 1, 2) - 93)*(105*std::pow(xi_2 + xi_3 - 1, 4)*(8*xi_1 + 4*xi_2 + 4*xi_3 - 4) - 2520*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3) + 2772*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 5)) ;
}
if( b == 212 ) {
  valDxi1 = (20*xi_3 - 1)*(104*xi_2*std::pow(xi_3 - 1, 2) + 33*std::pow(xi_3 - 1, 3) - 34*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2))*(-105*std::pow(xi_2 + xi_3 - 1, 4)*(8*xi_1 + 4*xi_2 + 4*xi_3 - 4) + 2520*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3) - 2772*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 5))/(16*xi_3 - 16) ;
}
if( b == 213 ) {
  valDxi1 = (1.0L/128.0L)*(105*std::pow(xi_2 + xi_3 - 1, 4)*(8*xi_1 + 4*xi_2 + 4*xi_3 - 4) - 2520*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3) + 2772*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 5))*(2550*xi_2*std::pow(xi_3 - 1, 3) - 1989*std::pow(-xi_3 + 1, 2)*std::pow(-2*xi_2 - xi_3 + 1, 2) + 969*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 3) + 1028*std::pow(xi_3 - 1, 4))/(xi_3 - 1) ;
}
if( b == 214 ) {
  valDxi1 = -(304*xi_3 + 95*std::pow(2*xi_3 - 1, 2) - 93)*(70*std::pow(xi_2 + xi_3 - 1, 7) - 1890*std::pow(xi_2 + xi_3 - 1, 5)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + 6930*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 4) - 6006*(xi_2 + xi_3 - 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 6))/(32*xi_2 + 32*xi_3 - 32) ;
}
if( b == 215 ) {
  valDxi1 = -(20*xi_3 - 1)*(17*xi_2 + xi_3 - 1)*(70*std::pow(xi_2 + xi_3 - 1, 7) - 1890*std::pow(xi_2 + xi_3 - 1, 5)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + 6930*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 4) - 6006*(xi_2 + xi_3 - 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 6))/(16*xi_2 + 16*xi_3 - 16) ;
}
if( b == 216 ) {
  valDxi1 = (1.0L/64.0L)*(540*xi_2*std::pow(xi_3 - 1, 2) + 171*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 2) + 167*std::pow(xi_3 - 1, 3))*(1890*std::pow(-xi_2 - xi_3 + 1, 5)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + 6006*(-xi_2 - xi_3 + 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 6) + 70*std::pow(xi_2 + xi_3 - 1, 7) + 6930*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 4))/((xi_3 - 1)*(xi_2 + xi_3 - 1)) ;
}
if( b == 217 ) {
  valDxi1 = (1.0L/128.0L)*(20*xi_3 - 1)*(-1260*std::pow(xi_2 + xi_3 - 1, 6)*(8*xi_1 + 4*xi_2 + 4*xi_3 - 4) + 55440*std::pow(xi_2 + xi_3 - 1, 4)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3) - 144144*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 5) + 102960*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 7)) ;
}
if( b == 218 ) {
  valDxi1 = (1.0L/128.0L)*(19*xi_2 + xi_3 - 1)*(-1260*std::pow(xi_2 + xi_3 - 1, 6)*(8*xi_1 + 4*xi_2 + 4*xi_3 - 4) + 55440*std::pow(xi_2 + xi_3 - 1, 4)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3) - 144144*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 5) + 102960*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 7)) ;
}
if( b == 219 ) {
  valDxi1 = (1.0L/128.0L)*(-180180*std::pow(-xi_2 - xi_3 + 1, 5)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 4) - 218790*(-xi_2 - xi_3 + 1)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 8) + 630*std::pow(xi_2 + xi_3 - 1, 9) - 27720*std::pow(xi_2 + xi_3 - 1, 7)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 2) - 360360*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 6))/(xi_2 + xi_3 - 1) ;
}
if( b == 220 ) {
  valDxi1 = 0 ;
}
if( b == 221 ) {
  valDxi1 = 0 ;
}
if( b == 222 ) {
  valDxi1 = 0 ;
}
if( b == 223 ) {
  valDxi1 = 0 ;
}
if( b == 224 ) {
  valDxi1 = 0 ;
}
if( b == 225 ) {
  valDxi1 = 0 ;
}
if( b == 226 ) {
  valDxi1 = 0 ;
}
if( b == 227 ) {
  valDxi1 = 0 ;
}
if( b == 228 ) {
  valDxi1 = 0 ;
}
if( b == 229 ) {
  valDxi1 = 0 ;
}
if( b == 230 ) {
  valDxi1 = 0 ;
}
if( b == 231 ) {
  valDxi1 = -1197.0L/32.0L*xi_3 + (124355.0L/64.0L)*std::pow(2*xi_3 - 1, 9) + (101745.0L/32.0L)*std::pow(2*xi_3 - 1, 8) - 14535.0L/16.0L*std::pow(2*xi_3 - 1, 7) - 6783.0L/2.0L*std::pow(2*xi_3 - 1, 6) - 18207.0L/32.0L*std::pow(2*xi_3 - 1, 5) + (16065.0L/16.0L)*std::pow(2*xi_3 - 1, 4) + (4305.0L/16.0L)*std::pow(2*xi_3 - 1, 3) - 315.0L/4.0L*std::pow(2*xi_3 - 1, 2) + 1239.0L/64.0L ;
}
if( b == 232 ) {
  valDxi1 = (1.0L/64.0L)*(5*xi_2 + xi_3 - 1)*(2640*xi_3 + 159885*std::pow(2*xi_3 - 1, 8) + 348840*std::pow(2*xi_3 - 1, 7) + 135660*std::pow(2*xi_3 - 1, 6) - 162792*std::pow(2*xi_3 - 1, 5) - 117810*std::pow(2*xi_3 - 1, 4) + 4760*std::pow(2*xi_3 - 1, 3) + 14700*std::pow(2*xi_3 - 1, 2) - 1499) ;
}
if( b == 233 ) {
  valDxi1 = -2*(36*xi_2*std::pow(xi_3 - 1, 2) + 17*std::pow(xi_3 - 1, 3) - 21*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2))*(182*xi_3 + 10659*std::pow(2*xi_3 - 1, 7) + 27132*std::pow(2*xi_3 - 1, 6) + 20349*std::pow(2*xi_3 - 1, 5) - 5355*std::pow(2*xi_3 - 1, 3) - 1428*std::pow(2*xi_3 - 1, 2) - 59)/(32*xi_3 - 32) ;
}
if( b == 234 ) {
  valDxi1 = (1.0L/64.0L)*(2*std::pow(xi_3 - 1, 3) - 21*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2) + 21*std::pow(2*xi_2 + xi_3 - 1, 3))*(-10812*xi_3 + 74613*std::pow(2*xi_3 - 1, 6) + 203490*std::pow(2*xi_3 - 1, 5) + 188955*std::pow(2*xi_3 - 1, 4) + 58140*std::pow(2*xi_3 - 1, 3) - 6885*std::pow(2*xi_3 - 1, 2) + 5011) ;
}
if( b == 235 ) {
  valDxi1 = (1.0L/64.0L)*(88*xi_2*std::pow(xi_3 - 1, 3) + 41*std::pow(xi_3 - 1, 4) - 18*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 2) - 180*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 3) + 165*std::pow(2*xi_2 + xi_3 - 1, 4))*(3150*xi_3 + 13167*std::pow(2*xi_3 - 1, 5) + 35910*std::pow(2*xi_3 - 1, 4) + 34770*std::pow(2*xi_3 - 1, 3) + 13680*std::pow(2*xi_3 - 1, 2) - 1669) ;
}
if( b == 236 ) {
  valDxi1 = (1.0L/256.0L)*(11704*xi_3 + 7315*std::pow(2*xi_3 - 1, 4) + 18620*std::pow(2*xi_3 - 1, 3) + 16530*std::pow(2*xi_3 - 1, 2) - 5209)*(-90*xi_2*std::pow(xi_3 - 1, 5) - 70*std::pow(-xi_3 + 1, 6) - 1485*std::pow(-xi_3 + 1, 2)*std::pow(-2*xi_2 - xi_3 + 1, 4) + 1287*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 5) + 630*std::pow(xi_3 - 1, 4)*std::pow(-2*xi_2 - xi_3 + 1, 2) + 330*std::pow(xi_3 - 1, 3)*std::pow(-2*xi_2 - xi_3 + 1, 3))/(xi_3 - 1) ;
}
if( b == 237 ) {
  valDxi1 = -2*(1170*xi_3 + 385*std::pow(2*xi_3 - 1, 3) + 840*std::pow(2*xi_3 - 1, 2) - 457)*(780*xi_2*std::pow(xi_3 - 1, 6) + 365*std::pow(xi_3 - 1, 7) + 165*std::pow(xi_3 - 1, 5)*std::pow(2*xi_2 + xi_3 - 1, 2) - 3740*std::pow(xi_3 - 1, 4)*std::pow(2*xi_2 + xi_3 - 1, 3) + 2145*std::pow(xi_3 - 1, 3)*std::pow(2*xi_2 + xi_3 - 1, 4) + 6006*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 5) - 5005*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 6))/(128*xi_3 - 128) ;
}
if( b == 238 ) {
  valDxi1 = (1.0L/64.0L)*(756*xi_3 + 231*std::pow(2*xi_3 - 1, 2) - 227)*(154*xi_2*std::pow(xi_3 - 1, 7) - 924*std::pow(-xi_3 + 1, 6)*std::pow(-2*xi_2 - xi_3 + 1, 2) - 6006*std::pow(-xi_3 + 1, 2)*std::pow(-2*xi_2 - xi_3 + 1, 6) + 4862*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 7) + 98*std::pow(xi_3 - 1, 8) + 5005*std::pow(xi_3 - 1, 4)*std::pow(-2*xi_2 - xi_3 + 1, 4) + 3003*std::pow(xi_3 - 1, 3)*std::pow(-2*xi_2 - xi_3 + 1, 5))/(xi_3 - 1) ;
}
if( b == 239 ) {
  valDxi1 = (1.0L/64.0L)*(22*xi_3 - 1)*(1680*xi_2*std::pow(xi_3 - 1, 7) + 13832*std::pow(-xi_3 + 1, 5)*std::pow(2*xi_2 + xi_3 - 1, 3) - 30940*std::pow(-xi_3 + 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 6) + 47736*(-xi_3 + 1)*std::pow(2*xi_2 + xi_3 - 1, 7) + 791*std::pow(xi_3 - 1, 8) + 1092*std::pow(xi_3 - 1, 6)*std::pow(2*xi_2 + xi_3 - 1, 2) + 2730*std::pow(xi_3 - 1, 4)*std::pow(2*xi_2 + xi_3 - 1, 4) + 50232*std::pow(xi_3 - 1, 3)*std::pow(2*xi_2 + xi_3 - 1, 5) + 37791*std::pow(2*xi_2 + xi_3 - 1, 8)) ;
}
if( b == 240 ) {
  valDxi1 = -2*(1638*xi_2*std::pow(xi_3 - 1, 9) + 966*std::pow(xi_3 - 1, 10) - 9828*std::pow(xi_3 - 1, 8)*std::pow(2*xi_2 + xi_3 - 1, 2) - 5460*std::pow(xi_3 - 1, 7)*std::pow(2*xi_2 + xi_3 - 1, 3) + 90090*std::pow(xi_3 - 1, 6)*std::pow(2*xi_2 + xi_3 - 1, 4) - 27846*std::pow(xi_3 - 1, 5)*std::pow(2*xi_2 + xi_3 - 1, 5) - 241332*std::pow(xi_3 - 1, 4)*std::pow(2*xi_2 + xi_3 - 1, 6) + 151164*std::pow(xi_3 - 1, 3)*std::pow(2*xi_2 + xi_3 - 1, 7) + 188955*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 8) - 146965*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 9))/(256*xi_3 - 256) ;
}
if( b == 241 ) {
  valDxi1 = (1.0L/256.0L)*(24*xi_1 + 12*xi_2 + 12*xi_3 - 12)*(2640*xi_3 + 159885*std::pow(2*xi_3 - 1, 8) + 348840*std::pow(2*xi_3 - 1, 7) + 135660*std::pow(2*xi_3 - 1, 6) - 162792*std::pow(2*xi_3 - 1, 5) - 117810*std::pow(2*xi_3 - 1, 4) + 4760*std::pow(2*xi_3 - 1, 3) + 14700*std::pow(2*xi_3 - 1, 2) - 1499) ;
}
if( b == 242 ) {
  valDxi1 = (1.0L/16.0L)*(7*xi_2 + xi_3 - 1)*(24*xi_1 + 12*xi_2 + 12*xi_3 - 12)*(182*xi_3 + 10659*std::pow(2*xi_3 - 1, 7) + 27132*std::pow(2*xi_3 - 1, 6) + 20349*std::pow(2*xi_3 - 1, 5) - 5355*std::pow(2*xi_3 - 1, 3) - 1428*std::pow(2*xi_3 - 1, 2) - 59) ;
}
if( b == 243 ) {
  valDxi1 = (20*xi_2*std::pow(xi_3 - 1, 2) + 8*std::pow(xi_3 - 1, 3) - 9*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2))*(-24*xi_1 - 12*xi_2 - 12*xi_3 + 12)*(-10812*xi_3 + 74613*std::pow(2*xi_3 - 1, 6) + 203490*std::pow(2*xi_3 - 1, 5) + 188955*std::pow(2*xi_3 - 1, 4) + 58140*std::pow(2*xi_3 - 1, 3) - 6885*std::pow(2*xi_3 - 1, 2) + 5011)/(128*xi_3 - 128) ;
}
if( b == 244 ) {
  valDxi1 = (1.0L/256.0L)*(24*xi_1 + 12*xi_2 + 12*xi_3 - 12)*(126*xi_2*std::pow(xi_3 - 1, 3) - 225*std::pow(-xi_3 + 1, 2)*std::pow(-2*xi_2 - xi_3 + 1, 2) + 165*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 3) + 68*std::pow(xi_3 - 1, 4))*(3150*xi_3 + 13167*std::pow(2*xi_3 - 1, 5) + 35910*std::pow(2*xi_3 - 1, 4) + 34770*std::pow(2*xi_3 - 1, 3) + 13680*std::pow(2*xi_3 - 1, 2) - 1669)/(xi_3 - 1) ;
}
if( b == 245 ) {
  valDxi1 = (1.0L/512.0L)*(24*xi_1 + 12*xi_2 + 12*xi_3 - 12)*(11704*xi_3 + 7315*std::pow(2*xi_3 - 1, 4) + 18620*std::pow(2*xi_3 - 1, 3) + 16530*std::pow(2*xi_3 - 1, 2) - 5209)*(200*xi_2*std::pow(xi_3 - 1, 3) + 71*std::pow(xi_3 - 1, 4) + 330*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 2) - 1100*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 3) + 715*std::pow(2*xi_2 + xi_3 - 1, 4)) ;
}
if( b == 246 ) {
  valDxi1 = (-24*xi_1 - 12*xi_2 - 12*xi_3 + 12)*(1170*xi_3 + 385*std::pow(2*xi_3 - 1, 3) + 840*std::pow(2*xi_3 - 1, 2) - 457)*(770*xi_2*std::pow(xi_3 - 1, 5) + 386*std::pow(xi_3 - 1, 6) - 990*std::pow(xi_3 - 1, 4)*std::pow(2*xi_2 + xi_3 - 1, 2) - 1430*std::pow(xi_3 - 1, 3)*std::pow(2*xi_2 + xi_3 - 1, 3) + 5005*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 4) - 3003*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 5))/(128*xi_3 - 128) ;
}
if( b == 247 ) {
  valDxi1 = (1.0L/128.0L)*(756*xi_3 + 231*std::pow(2*xi_3 - 1, 2) - 227)*(-24*xi_1 - 12*xi_2 - 12*xi_3 + 12)*(96*xi_2*std::pow(xi_3 - 1, 5) + 23*std::pow(xi_3 - 1, 6) + 780*std::pow(xi_3 - 1, 4)*std::pow(2*xi_2 + xi_3 - 1, 2) - 1820*std::pow(xi_3 - 1, 3)*std::pow(2*xi_2 + xi_3 - 1, 3) - 1365*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 4) + 5460*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 5) - 3094*std::pow(2*xi_2 + xi_3 - 1, 6)) ;
}
if( b == 248 ) {
  valDxi1 = (1.0L/64.0L)*(22*xi_3 - 1)*(24*xi_1 + 12*xi_2 + 12*xi_3 - 12)*(910*xi_2*std::pow(xi_3 - 1, 7) - 819*std::pow(-xi_3 + 1, 6)*std::pow(-2*xi_2 - xi_3 + 1, 2) + 4641*std::pow(-xi_3 + 1, 3)*std::pow(-2*xi_2 - xi_3 + 1, 5) - 23205*std::pow(-xi_3 + 1, 2)*std::pow(-2*xi_2 - xi_3 + 1, 6) + 12597*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 7) + 448*std::pow(xi_3 - 1, 8) + 5005*std::pow(xi_3 - 1, 5)*std::pow(-2*xi_2 - xi_3 + 1, 3) + 11375*std::pow(xi_3 - 1, 4)*std::pow(-2*xi_2 - xi_3 + 1, 4))/(xi_3 - 1) ;
}
if( b == 249 ) {
  valDxi1 = (1.0L/256.0L)*(24*xi_1 + 12*xi_2 + 12*xi_3 - 12)*(112*xi_2*std::pow(xi_3 - 1, 7) - 119*std::pow(-xi_3 + 1, 8) + 17080*std::pow(-xi_3 + 1, 5)*std::pow(2*xi_2 + xi_3 - 1, 3) - 55930*std::pow(-xi_3 + 1, 4)*std::pow(2*xi_2 + xi_3 - 1, 4) + 193800*(-xi_3 + 1)*std::pow(2*xi_2 + xi_3 - 1, 7) + 9660*std::pow(xi_3 - 1, 6)*std::pow(2*xi_2 + xi_3 - 1, 2) + 128520*std::pow(xi_3 - 1, 3)*std::pow(2*xi_2 + xi_3 - 1, 5) + 27132*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 6) + 101745*std::pow(2*xi_2 + xi_3 - 1, 8)) ;
}
if( b == 250 ) {
  valDxi1 = -(6*std::pow(xi_2 + xi_3 - 1, 3) - 30*(xi_2 + xi_3 - 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2))*(182*xi_3 + 10659*std::pow(2*xi_3 - 1, 7) + 27132*std::pow(2*xi_3 - 1, 6) + 20349*std::pow(2*xi_3 - 1, 5) - 5355*std::pow(2*xi_3 - 1, 3) - 1428*std::pow(2*xi_3 - 1, 2) - 59)/(16*xi_2 + 16*xi_3 - 16) ;
}
if( b == 251 ) {
  valDxi1 = (1.0L/128.0L)*(30*(-xi_2 - xi_3 + 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + 6*std::pow(xi_2 + xi_3 - 1, 3))*(9*xi_2 + xi_3 - 1)*(-10812*xi_3 + 74613*std::pow(2*xi_3 - 1, 6) + 203490*std::pow(2*xi_3 - 1, 5) + 188955*std::pow(2*xi_3 - 1, 4) + 58140*std::pow(2*xi_3 - 1, 3) - 6885*std::pow(2*xi_3 - 1, 2) + 5011)/(-xi_2 - xi_3 + 1) ;
}
if( b == 252 ) {
  valDxi1 = (1.0L/128.0L)*(30*(-xi_2 - xi_3 + 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + 6*std::pow(xi_2 + xi_3 - 1, 3))*(140*xi_2*std::pow(xi_3 - 1, 2) + 55*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 2) + 51*std::pow(xi_3 - 1, 3))*(3150*xi_3 + 13167*std::pow(2*xi_3 - 1, 5) + 35910*std::pow(2*xi_3 - 1, 4) + 34770*std::pow(2*xi_3 - 1, 3) + 13680*std::pow(2*xi_3 - 1, 2) - 1669)/((xi_3 - 1)*(xi_2 + xi_3 - 1)) ;
}
if( b == 253 ) {
  valDxi1 = -(6*std::pow(xi_2 + xi_3 - 1, 3) - 30*(xi_2 + xi_3 - 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2))*(198*xi_2*std::pow(xi_3 - 1, 2) + 92*std::pow(xi_3 - 1, 3) - 231*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2) + 143*std::pow(2*xi_2 + xi_3 - 1, 3))*(11704*xi_3 + 7315*std::pow(2*xi_3 - 1, 4) + 18620*std::pow(2*xi_3 - 1, 3) + 16530*std::pow(2*xi_3 - 1, 2) - 5209)/(128*xi_2 + 128*xi_3 - 128) ;
}
if( b == 254 ) {
  valDxi1 = (6*std::pow(xi_2 + xi_3 - 1, 3) - 30*(xi_2 + xi_3 - 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2))*(1170*xi_3 + 385*std::pow(2*xi_3 - 1, 3) + 840*std::pow(2*xi_3 - 1, 2) - 457)*(168*xi_2*std::pow(xi_3 - 1, 3) + 127*std::pow(xi_3 - 1, 4) - 1326*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 2) + 2548*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 3) - 1365*std::pow(2*xi_2 + xi_3 - 1, 4))/(64*xi_2 + 64*xi_3 - 64) ;
}
if( b == 255 ) {
  valDxi1 = (1.0L/64.0L)*(30*(-xi_2 - xi_3 + 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + 6*std::pow(xi_2 + xi_3 - 1, 3))*(756*xi_3 + 231*std::pow(2*xi_3 - 1, 2) - 227)*(390*xi_2*std::pow(xi_3 - 1, 5) + 1820*std::pow(-xi_3 + 1, 3)*std::pow(2*xi_2 + xi_3 - 1, 3) + 1547*(-xi_3 + 1)*std::pow(2*xi_2 + xi_3 - 1, 5) + 174*std::pow(xi_3 - 1, 6) + 3185*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 4))/((xi_3 - 1)*(xi_2 + xi_3 - 1)) ;
}
if( b == 256 ) {
  valDxi1 = -(22*xi_3 - 1)*(6*std::pow(xi_2 + xi_3 - 1, 3) - 30*(xi_2 + xi_3 - 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2))*(588*xi_2*std::pow(xi_3 - 1, 5) + 317*std::pow(xi_3 - 1, 6) - 1995*std::pow(xi_3 - 1, 4)*std::pow(2*xi_2 + xi_3 - 1, 2) + 980*std::pow(xi_3 - 1, 3)*std::pow(2*xi_2 + xi_3 - 1, 3) + 8925*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 4) - 14994*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 5) + 6783*std::pow(2*xi_2 + xi_3 - 1, 6))/(32*xi_2 + 32*xi_3 - 32) ;
}
if( b == 257 ) {
  valDxi1 = -1.0L/32.0L*(6*std::pow(xi_2 + xi_3 - 1, 3) - 30*(xi_2 + xi_3 - 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2))*(630*xi_2*std::pow(xi_3 - 1, 7) + 272*std::pow(xi_3 - 1, 8) + 1155*std::pow(xi_3 - 1, 6)*std::pow(2*xi_2 + xi_3 - 1, 2) - 7735*std::pow(xi_3 - 1, 5)*std::pow(2*xi_2 + xi_3 - 1, 3) + 5355*std::pow(xi_3 - 1, 4)*std::pow(2*xi_2 + xi_3 - 1, 4) + 20349*std::pow(xi_3 - 1, 3)*std::pow(2*xi_2 + xi_3 - 1, 5) - 33915*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 6) + 14535*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 7))/((xi_3 - 1)*(xi_2 + xi_3 - 1)) ;
}
if( b == 258 ) {
  valDxi1 = (1.0L/512.0L)*(-30*std::pow(xi_2 + xi_3 - 1, 2)*(8*xi_1 + 4*xi_2 + 4*xi_3 - 4) + 280*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3))*(-10812*xi_3 + 74613*std::pow(2*xi_3 - 1, 6) + 203490*std::pow(2*xi_3 - 1, 5) + 188955*std::pow(2*xi_3 - 1, 4) + 58140*std::pow(2*xi_3 - 1, 3) - 6885*std::pow(2*xi_3 - 1, 2) + 5011) ;
}
if( b == 259 ) {
  valDxi1 = (1.0L/128.0L)*(-30*std::pow(xi_2 + xi_3 - 1, 2)*(8*xi_1 + 4*xi_2 + 4*xi_3 - 4) + 280*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3))*(11*xi_2 + xi_3 - 1)*(3150*xi_3 + 13167*std::pow(2*xi_3 - 1, 5) + 35910*std::pow(2*xi_3 - 1, 4) + 34770*std::pow(2*xi_3 - 1, 3) + 13680*std::pow(2*xi_3 - 1, 2) - 1669) ;
}
if( b == 260 ) {
  valDxi1 = -(-30*std::pow(xi_2 + xi_3 - 1, 2)*(8*xi_1 + 4*xi_2 + 4*xi_3 - 4) + 280*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3))*(108*xi_2*std::pow(xi_3 - 1, 2) + 37*std::pow(xi_3 - 1, 3) - 39*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2))*(11704*xi_3 + 7315*std::pow(2*xi_3 - 1, 4) + 18620*std::pow(2*xi_3 - 1, 3) + 16530*std::pow(2*xi_3 - 1, 2) - 5209)/(256*xi_3 - 256) ;
}
if( b == 261 ) {
  valDxi1 = (1.0L/128.0L)*(-30*std::pow(xi_2 + xi_3 - 1, 2)*(8*xi_1 + 4*xi_2 + 4*xi_3 - 4) + 280*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3))*(1170*xi_3 + 385*std::pow(2*xi_3 - 1, 3) + 840*std::pow(2*xi_3 - 1, 2) - 457)*(858*xi_2*std::pow(xi_3 - 1, 3) - 819*std::pow(-xi_3 + 1, 2)*std::pow(-2*xi_2 - xi_3 + 1, 2) + 455*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 3) + 372*std::pow(xi_3 - 1, 4))/(xi_3 - 1) ;
}
if( b == 262 ) {
  valDxi1 = (1.0L/128.0L)*(-30*std::pow(xi_2 + xi_3 - 1, 2)*(8*xi_1 + 4*xi_2 + 4*xi_3 - 4) + 280*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3))*(756*xi_3 + 231*std::pow(2*xi_3 - 1, 2) - 227)*(-336*xi_2*std::pow(xi_3 - 1, 3) - 171*std::pow(xi_3 - 1, 4) + 840*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 2) - 1260*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 3) + 595*std::pow(2*xi_2 + xi_3 - 1, 4)) ;
}
if( b == 263 ) {
  valDxi1 = -(22*xi_3 - 1)*(-30*std::pow(xi_2 + xi_3 - 1, 2)*(8*xi_1 + 4*xi_2 + 4*xi_3 - 4) + 280*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3))*(330*xi_2*std::pow(xi_3 - 1, 5) + 114*std::pow(xi_3 - 1, 6) + 1170*std::pow(xi_3 - 1, 4)*std::pow(2*xi_2 + xi_3 - 1, 2) - 5270*std::pow(xi_3 - 1, 3)*std::pow(2*xi_2 + xi_3 - 1, 3) + 6885*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 4) - 2907*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 5))/(64*xi_3 - 64) ;
}
if( b == 264 ) {
  valDxi1 = (1.0L/64.0L)*(-30*std::pow(xi_2 + xi_3 - 1, 2)*(8*xi_1 + 4*xi_2 + 4*xi_3 - 4) + 280*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3))*(972*xi_2*std::pow(xi_3 - 1, 5) + 467*std::pow(xi_3 - 1, 6) - 1275*std::pow(xi_3 - 1, 4)*std::pow(2*xi_2 + xi_3 - 1, 2) - 3060*std::pow(xi_3 - 1, 3)*std::pow(2*xi_2 + xi_3 - 1, 3) + 14535*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 4) - 17442*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 5) + 6783*std::pow(2*xi_2 + xi_3 - 1, 6)) ;
}
if( b == 265 ) {
  valDxi1 = (1.0L/128.0L)*(-630*(-xi_2 - xi_3 + 1)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 4) + 30*std::pow(xi_2 + xi_3 - 1, 5) - 420*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 2))*(3150*xi_3 + 13167*std::pow(2*xi_3 - 1, 5) + 35910*std::pow(2*xi_3 - 1, 4) + 34770*std::pow(2*xi_3 - 1, 3) + 13680*std::pow(2*xi_3 - 1, 2) - 1669)/(xi_2 + xi_3 - 1) ;
}
if( b == 266 ) {
  valDxi1 = (1.0L/128.0L)*(13*xi_2 + xi_3 - 1)*(-630*(-xi_2 - xi_3 + 1)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 4) + 30*std::pow(xi_2 + xi_3 - 1, 5) - 420*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 2))*(11704*xi_3 + 7315*std::pow(2*xi_3 - 1, 4) + 18620*std::pow(2*xi_3 - 1, 3) + 16530*std::pow(2*xi_3 - 1, 2) - 5209)/(xi_2 + xi_3 - 1) ;
}
if( b == 267 ) {
  valDxi1 = -1.0L/64.0L*(308*xi_2*std::pow(xi_3 - 1, 2) + 101*std::pow(xi_3 - 1, 3) - 105*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2))*(30*std::pow(xi_2 + xi_3 - 1, 5) - 420*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + 630*(xi_2 + xi_3 - 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 4))*(1170*xi_3 + 385*std::pow(2*xi_3 - 1, 3) + 840*std::pow(2*xi_3 - 1, 2) - 457)/((xi_3 - 1)*(xi_2 + xi_3 - 1)) ;
}
if( b == 268 ) {
  valDxi1 = (1.0L/32.0L)*(756*xi_3 + 231*std::pow(2*xi_3 - 1, 2) - 227)*(-630*(-xi_2 - xi_3 + 1)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 4) + 30*std::pow(xi_2 + xi_3 - 1, 5) - 420*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 2))*(195*xi_2*std::pow(xi_3 - 1, 2) + 81*std::pow(xi_3 - 1, 3) - 165*(xi_3 - 1)*std::pow(-2*xi_2 - xi_3 + 1, 2) - 85*std::pow(-2*xi_2 - xi_3 + 1, 3))/(xi_2 + xi_3 - 1) ;
}
if( b == 269 ) {
  valDxi1 = (1.0L/32.0L)*(22*xi_3 - 1)*(-630*(-xi_2 - xi_3 + 1)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 4) + 30*std::pow(xi_2 + xi_3 - 1, 5) - 420*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 2))*(-968*xi_2*std::pow(xi_3 - 1, 3) - 455*std::pow(xi_3 - 1, 4) + 1734*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 2) - 2244*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 3) + 969*std::pow(2*xi_2 + xi_3 - 1, 4))/(xi_2 + xi_3 - 1) ;
}
if( b == 270 ) {
  valDxi1 = (1.0L/256.0L)*(30*std::pow(xi_2 + xi_3 - 1, 5) - 420*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + 630*(xi_2 + xi_3 - 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 4))*(2210*xi_2*std::pow(xi_3 - 1, 5) + 1358*std::pow(xi_3 - 1, 6) - 16830*std::pow(xi_3 - 1, 4)*std::pow(2*xi_2 + xi_3 - 1, 2) + 48450*std::pow(xi_3 - 1, 3)*std::pow(2*xi_2 + xi_3 - 1, 3) - 53295*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 4) + 20349*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 5))/((xi_3 - 1)*(xi_2 + xi_3 - 1)) ;
}
if( b == 271 ) {
  valDxi1 = (1.0L/256.0L)*(105*std::pow(xi_2 + xi_3 - 1, 4)*(8*xi_1 + 4*xi_2 + 4*xi_3 - 4) - 2520*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3) + 2772*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 5))*(11704*xi_3 + 7315*std::pow(2*xi_3 - 1, 4) + 18620*std::pow(2*xi_3 - 1, 3) + 16530*std::pow(2*xi_3 - 1, 2) - 5209) ;
}
if( b == 272 ) {
  valDxi1 = (1.0L/32.0L)*(15*xi_2 + xi_3 - 1)*(105*std::pow(xi_2 + xi_3 - 1, 4)*(8*xi_1 + 4*xi_2 + 4*xi_3 - 4) - 2520*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3) + 2772*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 5))*(1170*xi_3 + 385*std::pow(2*xi_3 - 1, 3) + 840*std::pow(2*xi_3 - 1, 2) - 457) ;
}
if( b == 273 ) {
  valDxi1 = (756*xi_3 + 231*std::pow(2*xi_3 - 1, 2) - 227)*(104*xi_2*std::pow(xi_3 - 1, 2) + 33*std::pow(xi_3 - 1, 3) - 34*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2))*(-105*std::pow(xi_2 + xi_3 - 1, 4)*(8*xi_1 + 4*xi_2 + 4*xi_3 - 4) + 2520*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3) - 2772*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 5))/(64*xi_3 - 64) ;
}
if( b == 274 ) {
  valDxi1 = (1.0L/128.0L)*(22*xi_3 - 1)*(105*std::pow(xi_2 + xi_3 - 1, 4)*(8*xi_1 + 4*xi_2 + 4*xi_3 - 4) - 2520*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3) + 2772*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 5))*(2550*xi_2*std::pow(xi_3 - 1, 3) - 1989*std::pow(-xi_3 + 1, 2)*std::pow(-2*xi_2 - xi_3 + 1, 2) + 969*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 3) + 1028*std::pow(xi_3 - 1, 4))/(xi_3 - 1) ;
}
if( b == 275 ) {
  valDxi1 = (1.0L/256.0L)*(105*std::pow(xi_2 + xi_3 - 1, 4)*(8*xi_1 + 4*xi_2 + 4*xi_3 - 4) - 2520*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3) + 2772*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 5))*(-8424*xi_2*std::pow(xi_3 - 1, 3) - 3803*std::pow(xi_3 - 1, 4) + 12654*std::pow(xi_3 - 1, 2)*std::pow(2*xi_2 + xi_3 - 1, 2) - 14820*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 3) + 5985*std::pow(2*xi_2 + xi_3 - 1, 4)) ;
}
if( b == 276 ) {
  valDxi1 = -(1170*xi_3 + 385*std::pow(2*xi_3 - 1, 3) + 840*std::pow(2*xi_3 - 1, 2) - 457)*(70*std::pow(xi_2 + xi_3 - 1, 7) - 1890*std::pow(xi_2 + xi_3 - 1, 5)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + 6930*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 4) - 6006*(xi_2 + xi_3 - 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 6))/(32*xi_2 + 32*xi_3 - 32) ;
}
if( b == 277 ) {
  valDxi1 = -(17*xi_2 + xi_3 - 1)*(756*xi_3 + 231*std::pow(2*xi_3 - 1, 2) - 227)*(70*std::pow(xi_2 + xi_3 - 1, 7) - 1890*std::pow(xi_2 + xi_3 - 1, 5)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + 6930*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 4) - 6006*(xi_2 + xi_3 - 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 6))/(64*xi_2 + 64*xi_3 - 64) ;
}
if( b == 278 ) {
  valDxi1 = (1.0L/64.0L)*(22*xi_3 - 1)*(540*xi_2*std::pow(xi_3 - 1, 2) + 171*(-xi_3 + 1)*std::pow(-2*xi_2 - xi_3 + 1, 2) + 167*std::pow(xi_3 - 1, 3))*(1890*std::pow(-xi_2 - xi_3 + 1, 5)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + 6006*(-xi_2 - xi_3 + 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 6) + 70*std::pow(xi_2 + xi_3 - 1, 7) + 6930*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 4))/((xi_3 - 1)*(xi_2 + xi_3 - 1)) ;
}
if( b == 279 ) {
  valDxi1 = -(1938*xi_2*std::pow(xi_3 - 1, 2) + 764*std::pow(xi_3 - 1, 3) - 1425*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2) + 665*std::pow(2*xi_2 + xi_3 - 1, 3))*(70*std::pow(xi_2 + xi_3 - 1, 7) - 1890*std::pow(xi_2 + xi_3 - 1, 5)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 2) + 6930*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 4) - 6006*(xi_2 + xi_3 - 1)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 6))/(64*xi_2 + 64*xi_3 - 64) ;
}
if( b == 280 ) {
  valDxi1 = (1.0L/512.0L)*(756*xi_3 + 231*std::pow(2*xi_3 - 1, 2) - 227)*(-1260*std::pow(xi_2 + xi_3 - 1, 6)*(8*xi_1 + 4*xi_2 + 4*xi_3 - 4) + 55440*std::pow(xi_2 + xi_3 - 1, 4)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3) - 144144*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 5) + 102960*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 7)) ;
}
if( b == 281 ) {
  valDxi1 = (1.0L/128.0L)*(22*xi_3 - 1)*(19*xi_2 + xi_3 - 1)*(-1260*std::pow(xi_2 + xi_3 - 1, 6)*(8*xi_1 + 4*xi_2 + 4*xi_3 - 4) + 55440*std::pow(xi_2 + xi_3 - 1, 4)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3) - 144144*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 5) + 102960*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 7)) ;
}
if( b == 282 ) {
  valDxi1 = -(340*xi_2*std::pow(xi_3 - 1, 2) + 103*std::pow(xi_3 - 1, 3) - 105*(xi_3 - 1)*std::pow(2*xi_2 + xi_3 - 1, 2))*(-1260*std::pow(xi_2 + xi_3 - 1, 6)*(8*xi_1 + 4*xi_2 + 4*xi_3 - 4) + 55440*std::pow(xi_2 + xi_3 - 1, 4)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3) - 144144*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 5) + 102960*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 7))/(256*xi_3 - 256) ;
}
if( b == 283 ) {
  valDxi1 = (1.0L/128.0L)*(22*xi_3 - 1)*(-180180*std::pow(-xi_2 - xi_3 + 1, 5)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 4) - 218790*(-xi_2 - xi_3 + 1)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 8) + 630*std::pow(xi_2 + xi_3 - 1, 9) - 27720*std::pow(xi_2 + xi_3 - 1, 7)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 2) - 360360*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 6))/(xi_2 + xi_3 - 1) ;
}
if( b == 284 ) {
  valDxi1 = (1.0L/128.0L)*(21*xi_2 + xi_3 - 1)*(-180180*std::pow(-xi_2 - xi_3 + 1, 5)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 4) - 218790*(-xi_2 - xi_3 + 1)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 8) + 630*std::pow(xi_2 + xi_3 - 1, 9) - 27720*std::pow(xi_2 + xi_3 - 1, 7)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 2) - 360360*std::pow(xi_2 + xi_3 - 1, 3)*std::pow(-2*xi_1 - xi_2 - xi_3 + 1, 6))/(xi_2 + xi_3 - 1) ;
}
if( b == 285 ) {
  valDxi1 = (3465.0L/256.0L)*std::pow(xi_2 + xi_3 - 1, 8)*(8*xi_1 + 4*xi_2 + 4*xi_3 - 4) - 15015.0L/16.0L*std::pow(xi_2 + xi_3 - 1, 6)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 3) + (135135.0L/32.0L)*std::pow(xi_2 + xi_3 - 1, 4)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 5) - 109395.0L/16.0L*std::pow(xi_2 + xi_3 - 1, 2)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 7) + (230945.0L/64.0L)*std::pow(2*xi_1 + xi_2 + xi_3 - 1, 9) ;
}
