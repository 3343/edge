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
 * Generated derivatives of triangular basis w.r.t. to xi2.
 **/

assert( b < 55 );
if( b == 0 ) {
  valDxi2 = 0 ;
}
if( b == 1 ) {
  valDxi2 = 3 ;
}
if( b == 2 ) {
  valDxi2 = 1 ;
}
if( b == 3 ) {
  valDxi2 = 20*xi_2 - 8 ;
}
if( b == 4 ) {
  valDxi2 = 10*xi_1 + 10*xi_2 - 6 ;
}
if( b == 5 ) {
  valDxi2 = 6*xi_1 + 2*xi_2 - 2 ;
}
if( b == 6 ) {
  valDxi2 = 105*std::pow(xi_2, 2) - 90*xi_2 + 15 ;
}
if( b == 7 ) {
  valDxi2 = 9*xi_2 + (21.0L/4.0L)*std::pow(2*xi_2 - 1, 2) + (1.0L/4.0L)*(168*xi_2 - 48)*(2*xi_1 + xi_2 - 1) - 17.0L/4.0L ;
}
if( b == 8 ) {
  valDxi2 = -7.0L/2.0L*std::pow(xi_2 - 1, 2) + (1.0L/2.0L)*(7*xi_2 - 1)*(12*xi_1 + 4*xi_2 - 4) + (21.0L/2.0L)*std::pow(2*xi_1 + xi_2 - 1, 2) ;
}
if( b == 9 ) {
  valDxi2 = -(18*xi_1*std::pow(xi_2 - 1, 2) + 12*std::pow(xi_2 - 1, 3) - 15*(xi_2 - 1)*std::pow(2*xi_1 + xi_2 - 1, 2) - 5*std::pow(2*xi_1 + xi_2 - 1, 3))/(2*xi_2 - 2) + 2*(6*xi_1*std::pow(xi_2 - 1, 3) + 3*std::pow(xi_2 - 1, 4) - 5*(xi_2 - 1)*std::pow(2*xi_1 + xi_2 - 1, 3))/std::pow(2*xi_2 - 2, 2) ;
}
if( b == 10 ) {
  valDxi2 = 504*std::pow(xi_2, 3) - 672*std::pow(xi_2, 2) + 252*xi_2 - 24 ;
}
if( b == 11 ) {
  valDxi2 = (21.0L/2.0L)*std::pow(2*xi_2 - 1, 3) + (21.0L/2.0L)*std::pow(2*xi_2 - 1, 2) + (1.0L/2.0L)*(2*xi_1 + xi_2 - 1)*(168*xi_2 + 126*std::pow(2*xi_2 - 1, 2) - 84) - 1 ;
}
if( b == 12 ) {
  valDxi2 = (1.0L/2.0L)*(72*xi_2 - 16)*(-std::pow(xi_2 - 1, 2) + 3*std::pow(2*xi_1 + xi_2 - 1, 2)) + (1.0L/2.0L)*(12*xi_1 + 4*xi_2 - 4)*(20*xi_2 + 9*std::pow(2*xi_2 - 1, 2) - 8) ;
}
if( b == 13 ) {
  valDxi2 = -(9*xi_2 - 1)*(18*xi_1*std::pow(xi_2 - 1, 2) + 12*std::pow(xi_2 - 1, 3) - 15*(xi_2 - 1)*std::pow(2*xi_1 + xi_2 - 1, 2) - 5*std::pow(2*xi_1 + xi_2 - 1, 3))/(2*xi_2 - 2) - 9*(6*xi_1*std::pow(xi_2 - 1, 3) + 3*std::pow(xi_2 - 1, 4) - 5*(xi_2 - 1)*std::pow(2*xi_1 + xi_2 - 1, 3))/(2*xi_2 - 2) + 2*(9*xi_2 - 1)*(6*xi_1*std::pow(xi_2 - 1, 3) + 3*std::pow(xi_2 - 1, 4) - 5*(xi_2 - 1)*std::pow(2*xi_1 + xi_2 - 1, 3))/std::pow(2*xi_2 - 2, 2) ;
}
if( b == 14 ) {
  valDxi2 = (3.0L/2.0L)*std::pow(xi_2 - 1, 3) - 15.0L/4.0L*std::pow(xi_2 - 1, 2)*(4*xi_1 + 2*xi_2 - 2) - 15.0L/4.0L*(2*xi_2 - 2)*std::pow(2*xi_1 + xi_2 - 1, 2) + (35.0L/2.0L)*std::pow(2*xi_1 + xi_2 - 1, 3) ;
}
if( b == 15 ) {
  valDxi2 = 2310*std::pow(xi_2, 4) - 4200*std::pow(xi_2, 3) + 2520*std::pow(xi_2, 2) - 560*xi_2 + 35 ;
}
if( b == 16 ) {
  valDxi2 = -11*xi_2 + (165.0L/8.0L)*std::pow(2*xi_2 - 1, 4) + (45.0L/2.0L)*std::pow(2*xi_2 - 1, 3) - 9.0L/4.0L*std::pow(2*xi_2 - 1, 2) + (1.0L/8.0L)*(2*xi_1 + xi_2 - 1)*(-144*xi_2 + 1320*std::pow(2*xi_2 - 1, 3) + 1080*std::pow(2*xi_2 - 1, 2) - 16) + 41.0L/8.0L ;
}
if( b == 17 ) {
  valDxi2 = (1.0L/16.0L)*(-std::pow(xi_2 - 1, 2) + 3*std::pow(2*xi_1 + xi_2 - 1, 2))*(1800*xi_2 + 990*std::pow(2*xi_2 - 1, 2) - 774) + (1.0L/16.0L)*(12*xi_1 + 4*xi_2 - 4)*(126*xi_2 + 165*std::pow(2*xi_2 - 1, 3) + 225*std::pow(2*xi_2 - 1, 2) - 68) ;
}
if( b == 18 ) {
  valDxi2 = -(440*xi_2 - 80)*(6*xi_1*std::pow(xi_2 - 1, 3) + 3*std::pow(xi_2 - 1, 4) - 5*(xi_2 - 1)*std::pow(2*xi_1 + xi_2 - 1, 3))/(8*xi_2 - 8) - (140*xi_2 + 55*std::pow(2*xi_2 - 1, 2) - 51)*(18*xi_1*std::pow(xi_2 - 1, 2) + 12*std::pow(xi_2 - 1, 3) - 15*(xi_2 - 1)*std::pow(2*xi_1 + xi_2 - 1, 2) - 5*std::pow(2*xi_1 + xi_2 - 1, 3))/(8*xi_2 - 8) + 8*(140*xi_2 + 55*std::pow(2*xi_2 - 1, 2) - 51)*(6*xi_1*std::pow(xi_2 - 1, 3) + 3*std::pow(xi_2 - 1, 4) - 5*(xi_2 - 1)*std::pow(2*xi_1 + xi_2 - 1, 3))/std::pow(8*xi_2 - 8, 2) ;
}
if( b == 19 ) {
  valDxi2 = (33.0L/8.0L)*std::pow(xi_2 - 1, 4) - 165.0L/4.0L*std::pow(xi_2 - 1, 2)*std::pow(2*xi_1 + xi_2 - 1, 2) + (1.0L/8.0L)*(11*xi_2 - 1)*(12*std::pow(xi_2 - 1, 3) - 30*std::pow(xi_2 - 1, 2)*(4*xi_1 + 2*xi_2 - 2) - 30*(2*xi_2 - 2)*std::pow(2*xi_1 + xi_2 - 1, 2) + 140*std::pow(2*xi_1 + xi_2 - 1, 3)) + (385.0L/8.0L)*std::pow(2*xi_1 + xi_2 - 1, 4) ;
}
if( b == 20 ) {
  valDxi2 = (1.0L/8.0L)*(150*xi_1*std::pow(xi_2 - 1, 4) - 315*(-xi_2 + 1)*std::pow(-2*xi_1 - xi_2 + 1, 4) + 90*std::pow(xi_2 - 1, 5) - 210*std::pow(xi_2 - 1, 3)*std::pow(-2*xi_1 - xi_2 + 1, 2) + 210*std::pow(xi_2 - 1, 2)*std::pow(-2*xi_1 - xi_2 + 1, 3) - 63*std::pow(-2*xi_1 - xi_2 + 1, 5))/(xi_2 - 1) - 1.0L/8.0L*(30*xi_1*std::pow(xi_2 - 1, 5) + 63*(-xi_2 + 1)*std::pow(-2*xi_1 - xi_2 + 1, 5) + 15*std::pow(xi_2 - 1, 6) + 70*std::pow(xi_2 - 1, 3)*std::pow(-2*xi_1 - xi_2 + 1, 3))/std::pow(xi_2 - 1, 2) ;
}
if( b == 21 ) {
  valDxi2 = 10296*std::pow(xi_2, 5) - 23760*std::pow(xi_2, 4) + 19800*std::pow(xi_2, 3) - 7200*std::pow(xi_2, 2) + 1080*xi_2 - 48 ;
}
if( b == 22 ) {
  valDxi2 = -45.0L/16.0L*xi_2 + (1287.0L/32.0L)*std::pow(2*xi_2 - 1, 5) + (1485.0L/32.0L)*std::pow(2*xi_2 - 1, 4) - 165.0L/16.0L*std::pow(2*xi_2 - 1, 3) - 315.0L/16.0L*std::pow(2*xi_2 - 1, 2) - 1.0L/32.0L*(2*xi_1 + xi_2 - 1)*(5040*xi_2 - 12870*std::pow(2*xi_2 - 1, 4) - 11880*std::pow(2*xi_2 - 1, 3) + 1980*std::pow(2*xi_2 - 1, 2) - 2430) + 35.0L/16.0L ;
}
if( b == 23 ) {
  valDxi2 = (1.0L/32.0L)*(-std::pow(xi_2 - 1, 2) + 3*std::pow(2*xi_1 + xi_2 - 1, 2))*(2640*xi_2 + 5720*std::pow(2*xi_2 - 1, 3) + 6600*std::pow(2*xi_2 - 1, 2) - 1520) + (1.0L/32.0L)*(12*xi_1 + 4*xi_2 - 4)*(-200*xi_2 + 715*std::pow(2*xi_2 - 1, 4) + 1100*std::pow(2*xi_2 - 1, 3) + 330*std::pow(2*xi_2 - 1, 2) + 71) ;
}
if( b == 24 ) {
  valDxi2 = -(1848*xi_2 + 858*std::pow(2*xi_2 - 1, 2) - 726)*(6*xi_1*std::pow(xi_2 - 1, 3) + 3*std::pow(xi_2 - 1, 4) - 5*(xi_2 - 1)*std::pow(2*xi_1 + xi_2 - 1, 3))/(8*xi_2 - 8) - (198*xi_2 + 143*std::pow(2*xi_2 - 1, 3) + 231*std::pow(2*xi_2 - 1, 2) - 92)*(18*xi_1*std::pow(xi_2 - 1, 2) + 12*std::pow(xi_2 - 1, 3) - 15*(xi_2 - 1)*std::pow(2*xi_1 + xi_2 - 1, 2) - 5*std::pow(2*xi_1 + xi_2 - 1, 3))/(8*xi_2 - 8) + 8*(6*xi_1*std::pow(xi_2 - 1, 3) + 3*std::pow(xi_2 - 1, 4) - 5*(xi_2 - 1)*std::pow(2*xi_1 + xi_2 - 1, 3))*(198*xi_2 + 143*std::pow(2*xi_2 - 1, 3) + 231*std::pow(2*xi_2 - 1, 2) - 92)/std::pow(8*xi_2 - 8, 2) ;
}
if( b == 25 ) {
  valDxi2 = (1.0L/16.0L)*(312*xi_2 - 48)*(3*std::pow(xi_2 - 1, 4) - 30*std::pow(xi_2 - 1, 2)*std::pow(2*xi_1 + xi_2 - 1, 2) + 35*std::pow(2*xi_1 + xi_2 - 1, 4)) + (1.0L/16.0L)*(108*xi_2 + 39*std::pow(2*xi_2 - 1, 2) - 37)*(12*std::pow(xi_2 - 1, 3) - 30*std::pow(xi_2 - 1, 2)*(4*xi_1 + 2*xi_2 - 2) - 30*(2*xi_2 - 2)*std::pow(2*xi_1 + xi_2 - 1, 2) + 140*std::pow(2*xi_1 + xi_2 - 1, 3)) ;
}
if( b == 26 ) {
  valDxi2 = (1.0L/8.0L)*(13*xi_2 - 1)*(150*xi_1*std::pow(xi_2 - 1, 4) - 315*(-xi_2 + 1)*std::pow(-2*xi_1 - xi_2 + 1, 4) + 90*std::pow(xi_2 - 1, 5) - 210*std::pow(xi_2 - 1, 3)*std::pow(-2*xi_1 - xi_2 + 1, 2) + 210*std::pow(xi_2 - 1, 2)*std::pow(-2*xi_1 - xi_2 + 1, 3) - 63*std::pow(-2*xi_1 - xi_2 + 1, 5))/(xi_2 - 1) + (13.0L/8.0L)*(30*xi_1*std::pow(xi_2 - 1, 5) + 63*(-xi_2 + 1)*std::pow(-2*xi_1 - xi_2 + 1, 5) + 15*std::pow(xi_2 - 1, 6) + 70*std::pow(xi_2 - 1, 3)*std::pow(-2*xi_1 - xi_2 + 1, 3))/(xi_2 - 1) - 1.0L/8.0L*(13*xi_2 - 1)*(30*xi_1*std::pow(xi_2 - 1, 5) + 63*(-xi_2 + 1)*std::pow(-2*xi_1 - xi_2 + 1, 5) + 15*std::pow(xi_2 - 1, 6) + 70*std::pow(xi_2 - 1, 3)*std::pow(-2*xi_1 - xi_2 + 1, 3))/std::pow(xi_2 - 1, 2) ;
}
if( b == 27 ) {
  valDxi2 = -15.0L/8.0L*std::pow(xi_2 - 1, 5) + (105.0L/16.0L)*std::pow(xi_2 - 1, 4)*(4*xi_1 + 2*xi_2 - 2) + (105.0L/4.0L)*std::pow(xi_2 - 1, 3)*std::pow(2*xi_1 + xi_2 - 1, 2) - 315.0L/4.0L*std::pow(xi_2 - 1, 2)*std::pow(2*xi_1 + xi_2 - 1, 3) - 315.0L/16.0L*(2*xi_2 - 2)*std::pow(2*xi_1 + xi_2 - 1, 4) + (693.0L/8.0L)*std::pow(2*xi_1 + xi_2 - 1, 5) ;
}
if( b == 28 ) {
  valDxi2 = 45045*std::pow(xi_2, 6) - 126126*std::pow(xi_2, 5) + 135135*std::pow(xi_2, 4) - 69300*std::pow(xi_2, 3) + 17325*std::pow(xi_2, 2) - 1890*xi_2 + 63 ;
}
if( b == 29 ) {
  valDxi2 = (195.0L/16.0L)*xi_2 + (5005.0L/64.0L)*std::pow(2*xi_2 - 1, 6) + (3003.0L/32.0L)*std::pow(2*xi_2 - 1, 5) - 2145.0L/64.0L*std::pow(2*xi_2 - 1, 4) - 935.0L/16.0L*std::pow(2*xi_2 - 1, 3) - 165.0L/64.0L*std::pow(2*xi_2 - 1, 2) - 1.0L/64.0L*(2*xi_1 + xi_2 - 1)*(1320*xi_2 - 60060*std::pow(2*xi_2 - 1, 5) - 60060*std::pow(2*xi_2 - 1, 4) + 17160*std::pow(2*xi_2 - 1, 3) + 22440*std::pow(2*xi_2 - 1, 2) - 1440) - 365.0L/64.0L ;
}
if( b == 30 ) {
  valDxi2 = (1.0L/64.0L)*(-std::pow(xi_2 - 1, 2) + 3*std::pow(2*xi_1 + xi_2 - 1, 2))*(-7920*xi_2 + 30030*std::pow(2*xi_2 - 1, 4) + 40040*std::pow(2*xi_2 - 1, 3) + 8580*std::pow(2*xi_2 - 1, 2) + 3190) + (1.0L/64.0L)*(12*xi_1 + 4*xi_2 - 4)*(-770*xi_2 + 3003*std::pow(2*xi_2 - 1, 5) + 5005*std::pow(2*xi_2 - 1, 4) + 1430*std::pow(2*xi_2 - 1, 3) - 990*std::pow(2*xi_2 - 1, 2) + 386) ;
}
if( b == 31 ) {
  valDxi2 = -(6*xi_1*std::pow(xi_2 - 1, 3) + 3*std::pow(xi_2 - 1, 4) - 5*(xi_2 - 1)*std::pow(2*xi_1 + xi_2 - 1, 3))*(10608*xi_2 + 10920*std::pow(2*xi_2 - 1, 3) + 15288*std::pow(2*xi_2 - 1, 2) - 5136)/(32*xi_2 - 32) - (18*xi_1*std::pow(xi_2 - 1, 2) + 12*std::pow(xi_2 - 1, 3) - 15*(xi_2 - 1)*std::pow(2*xi_1 + xi_2 - 1, 2) - 5*std::pow(2*xi_1 + xi_2 - 1, 3))*(168*xi_2 + 1365*std::pow(2*xi_2 - 1, 4) + 2548*std::pow(2*xi_2 - 1, 3) + 1326*std::pow(2*xi_2 - 1, 2) - 127)/(32*xi_2 - 32) + 32*(6*xi_1*std::pow(xi_2 - 1, 3) + 3*std::pow(xi_2 - 1, 4) - 5*(xi_2 - 1)*std::pow(2*xi_1 + xi_2 - 1, 3))*(168*xi_2 + 1365*std::pow(2*xi_2 - 1, 4) + 2548*std::pow(2*xi_2 - 1, 3) + 1326*std::pow(2*xi_2 - 1, 2) - 127)/std::pow(32*xi_2 - 32, 2) ;
}
if( b == 32 ) {
  valDxi2 = (1.0L/64.0L)*(6552*xi_2 + 2730*std::pow(2*xi_2 - 1, 2) - 2418)*(3*std::pow(xi_2 - 1, 4) - 30*std::pow(xi_2 - 1, 2)*std::pow(2*xi_1 + xi_2 - 1, 2) + 35*std::pow(2*xi_1 + xi_2 - 1, 4)) + (1.0L/64.0L)*(858*xi_2 + 455*std::pow(2*xi_2 - 1, 3) + 819*std::pow(2*xi_2 - 1, 2) - 372)*(12*std::pow(xi_2 - 1, 3) - 30*std::pow(xi_2 - 1, 2)*(4*xi_1 + 2*xi_2 - 2) - 30*(2*xi_2 - 2)*std::pow(2*xi_1 + xi_2 - 1, 2) + 140*std::pow(2*xi_1 + xi_2 - 1, 3)) ;
}
if( b == 33 ) {
  valDxi2 = (1.0L/32.0L)*(840*xi_2 - 112)*(30*xi_1*std::pow(xi_2 - 1, 5) + 63*(-xi_2 + 1)*std::pow(-2*xi_1 - xi_2 + 1, 5) + 15*std::pow(xi_2 - 1, 6) + 70*std::pow(xi_2 - 1, 3)*std::pow(-2*xi_1 - xi_2 + 1, 3))/(xi_2 - 1) + (1.0L/32.0L)*(308*xi_2 + 105*std::pow(2*xi_2 - 1, 2) - 101)*(150*xi_1*std::pow(xi_2 - 1, 4) - 315*(-xi_2 + 1)*std::pow(-2*xi_1 - xi_2 + 1, 4) + 90*std::pow(xi_2 - 1, 5) - 210*std::pow(xi_2 - 1, 3)*std::pow(-2*xi_1 - xi_2 + 1, 2) + 210*std::pow(xi_2 - 1, 2)*std::pow(-2*xi_1 - xi_2 + 1, 3) - 63*std::pow(-2*xi_1 - xi_2 + 1, 5))/(xi_2 - 1) - 1.0L/32.0L*(308*xi_2 + 105*std::pow(2*xi_2 - 1, 2) - 101)*(30*xi_1*std::pow(xi_2 - 1, 5) + 63*(-xi_2 + 1)*std::pow(-2*xi_1 - xi_2 + 1, 5) + 15*std::pow(xi_2 - 1, 6) + 70*std::pow(xi_2 - 1, 3)*std::pow(-2*xi_1 - xi_2 + 1, 3))/std::pow(xi_2 - 1, 2) ;
}
if( b == 34 ) {
  valDxi2 = -75.0L/16.0L*std::pow(xi_2 - 1, 6) + (1575.0L/16.0L)*std::pow(xi_2 - 1, 4)*std::pow(2*xi_1 + xi_2 - 1, 2) - 4725.0L/16.0L*std::pow(xi_2 - 1, 2)*std::pow(2*xi_1 + xi_2 - 1, 4) + (1.0L/16.0L)*(15*xi_2 - 1)*(-30*std::pow(xi_2 - 1, 5) + 105*std::pow(xi_2 - 1, 4)*(4*xi_1 + 2*xi_2 - 2) + 420*std::pow(xi_2 - 1, 3)*std::pow(2*xi_1 + xi_2 - 1, 2) - 1260*std::pow(xi_2 - 1, 2)*std::pow(2*xi_1 + xi_2 - 1, 3) - 315*(2*xi_2 - 2)*std::pow(2*xi_1 + xi_2 - 1, 4) + 1386*std::pow(2*xi_1 + xi_2 - 1, 5)) + (3465.0L/16.0L)*std::pow(2*xi_1 + xi_2 - 1, 6) ;
}
if( b == 35 ) {
  valDxi2 = (1.0L/16.0L)*(-490*xi_1*std::pow(xi_2 - 1, 6) - 3003*(-xi_2 + 1)*std::pow(2*xi_1 + xi_2 - 1, 6) - 280*std::pow(xi_2 - 1, 7) + 945*std::pow(xi_2 - 1, 5)*std::pow(2*xi_1 + xi_2 - 1, 2) + 1575*std::pow(xi_2 - 1, 4)*std::pow(2*xi_1 + xi_2 - 1, 3) - 3465*std::pow(xi_2 - 1, 3)*std::pow(2*xi_1 + xi_2 - 1, 4) - 2079*std::pow(xi_2 - 1, 2)*std::pow(2*xi_1 + xi_2 - 1, 5) + 429*std::pow(2*xi_1 + xi_2 - 1, 7))/(xi_2 - 1) - 1.0L/16.0L*(-70*xi_1*std::pow(xi_2 - 1, 7) - 429*(-xi_2 + 1)*std::pow(2*xi_1 + xi_2 - 1, 7) - 35*std::pow(xi_2 - 1, 8) + 315*std::pow(xi_2 - 1, 5)*std::pow(2*xi_1 + xi_2 - 1, 3) - 693*std::pow(xi_2 - 1, 3)*std::pow(2*xi_1 + xi_2 - 1, 5))/std::pow(xi_2 - 1, 2) ;
}
if( b == 36 ) {
  valDxi2 = 194480*std::pow(xi_2, 7) - 640640*std::pow(xi_2, 6) + 840840*std::pow(xi_2, 5) - 560560*std::pow(xi_2, 4) + 200200*std::pow(xi_2, 3) - 36960*std::pow(xi_2, 2) + 3080*xi_2 - 80 ;
}
if( b == 37 ) {
  valDxi2 = (77.0L/16.0L)*xi_2 + (2431.0L/16.0L)*std::pow(2*xi_2 - 1, 7) + (3003.0L/16.0L)*std::pow(2*xi_2 - 1, 6) - 3003.0L/32.0L*std::pow(2*xi_2 - 1, 5) - 5005.0L/32.0L*std::pow(2*xi_2 - 1, 4) + (231.0L/8.0L)*std::pow(2*xi_2 - 1, 2) + (1.0L/32.0L)*(2*xi_1 + xi_2 - 1)*(7392*xi_2 + 68068*std::pow(2*xi_2 - 1, 6) + 72072*std::pow(2*xi_2 - 1, 5) - 30030*std::pow(2*xi_2 - 1, 4) - 40040*std::pow(2*xi_2 - 1, 3) - 3542) - 49.0L/16.0L ;
}
if( b == 38 ) {
  valDxi2 = (1.0L/32.0L)*(-std::pow(xi_2 - 1, 2) + 3*std::pow(2*xi_1 + xi_2 - 1, 2))*(-6240*xi_2 + 37128*std::pow(2*xi_2 - 1, 5) + 54600*std::pow(2*xi_2 - 1, 4) + 10920*std::pow(2*xi_2 - 1, 3) - 10920*std::pow(2*xi_2 - 1, 2) + 3216) + (1.0L/32.0L)*(12*xi_1 + 4*xi_2 - 4)*(96*xi_2 + 3094*std::pow(2*xi_2 - 1, 6) + 5460*std::pow(2*xi_2 - 1, 5) + 1365*std::pow(2*xi_2 - 1, 4) - 1820*std::pow(2*xi_2 - 1, 3) - 780*std::pow(2*xi_2 - 1, 2) - 23) ;
}
if( b == 39 ) {
  valDxi2 = -(6*xi_1*std::pow(xi_2 - 1, 3) + 3*std::pow(xi_2 - 1, 4) - 5*(xi_2 - 1)*std::pow(2*xi_1 + xi_2 - 1, 3))*(15470*std::pow(2*xi_2 - 1, 4) + 25480*std::pow(2*xi_2 - 1, 3) + 10920*std::pow(2*xi_2 - 1, 2) - 390)/(16*xi_2 - 16) - (18*xi_1*std::pow(xi_2 - 1, 2) + 12*std::pow(xi_2 - 1, 3) - 15*(xi_2 - 1)*std::pow(2*xi_1 + xi_2 - 1, 2) - 5*std::pow(2*xi_1 + xi_2 - 1, 3))*(-390*xi_2 + 1547*std::pow(2*xi_2 - 1, 5) + 3185*std::pow(2*xi_2 - 1, 4) + 1820*std::pow(2*xi_2 - 1, 3) + 174)/(16*xi_2 - 16) + 16*(6*xi_1*std::pow(xi_2 - 1, 3) + 3*std::pow(xi_2 - 1, 4) - 5*(xi_2 - 1)*std::pow(2*xi_1 + xi_2 - 1, 3))*(-390*xi_2 + 1547*std::pow(2*xi_2 - 1, 5) + 3185*std::pow(2*xi_2 - 1, 4) + 1820*std::pow(2*xi_2 - 1, 3) + 174)/std::pow(16*xi_2 - 16, 2) ;
}
if( b == 40 ) {
  valDxi2 = (1.0L/32.0L)*(3*std::pow(xi_2 - 1, 4) - 30*std::pow(xi_2 - 1, 2)*std::pow(2*xi_1 + xi_2 - 1, 2) + 35*std::pow(2*xi_1 + xi_2 - 1, 4))*(6720*xi_2 + 4760*std::pow(2*xi_2 - 1, 3) + 7560*std::pow(2*xi_2 - 1, 2) - 3024) + (1.0L/32.0L)*(12*std::pow(xi_2 - 1, 3) - 30*std::pow(xi_2 - 1, 2)*(4*xi_1 + 2*xi_2 - 2) - 30*(2*xi_2 - 2)*std::pow(2*xi_1 + xi_2 - 1, 2) + 140*std::pow(2*xi_1 + xi_2 - 1, 3))*(336*xi_2 + 595*std::pow(2*xi_2 - 1, 4) + 1260*std::pow(2*xi_2 - 1, 3) + 840*std::pow(2*xi_2 - 1, 2) - 171) ;
}
if( b == 41 ) {
  valDxi2 = (1.0L/8.0L)*(1320*xi_2 + 510*std::pow(2*xi_2 - 1, 2) - 465)*(30*xi_1*std::pow(xi_2 - 1, 5) + 63*(-xi_2 + 1)*std::pow(-2*xi_1 - xi_2 + 1, 5) + 15*std::pow(xi_2 - 1, 6) + 70*std::pow(xi_2 - 1, 3)*std::pow(-2*xi_1 - xi_2 + 1, 3))/(xi_2 - 1) + (1.0L/8.0L)*(195*xi_2 + 85*std::pow(2*xi_2 - 1, 3) + 165*std::pow(2*xi_2 - 1, 2) - 81)*(150*xi_1*std::pow(xi_2 - 1, 4) - 315*(-xi_2 + 1)*std::pow(-2*xi_1 - xi_2 + 1, 4) + 90*std::pow(xi_2 - 1, 5) - 210*std::pow(xi_2 - 1, 3)*std::pow(-2*xi_1 - xi_2 + 1, 2) + 210*std::pow(xi_2 - 1, 2)*std::pow(-2*xi_1 - xi_2 + 1, 3) - 63*std::pow(-2*xi_1 - xi_2 + 1, 5))/(xi_2 - 1) - 1.0L/8.0L*(195*xi_2 + 85*std::pow(2*xi_2 - 1, 3) + 165*std::pow(2*xi_2 - 1, 2) - 81)*(30*xi_1*std::pow(xi_2 - 1, 5) + 63*(-xi_2 + 1)*std::pow(-2*xi_1 - xi_2 + 1, 5) + 15*std::pow(xi_2 - 1, 6) + 70*std::pow(xi_2 - 1, 3)*std::pow(-2*xi_1 - xi_2 + 1, 3))/std::pow(xi_2 - 1, 2) ;
}
if( b == 42 ) {
  valDxi2 = (1.0L/16.0L)*(272*xi_2 - 32)*(-5*std::pow(xi_2 - 1, 6) + 105*std::pow(xi_2 - 1, 4)*std::pow(2*xi_1 + xi_2 - 1, 2) - 315*std::pow(xi_2 - 1, 2)*std::pow(2*xi_1 + xi_2 - 1, 4) + 231*std::pow(2*xi_1 + xi_2 - 1, 6)) + (1.0L/16.0L)*(104*xi_2 + 34*std::pow(2*xi_2 - 1, 2) - 33)*(-30*std::pow(xi_2 - 1, 5) + 105*std::pow(xi_2 - 1, 4)*(4*xi_1 + 2*xi_2 - 2) + 420*std::pow(xi_2 - 1, 3)*std::pow(2*xi_1 + xi_2 - 1, 2) - 1260*std::pow(xi_2 - 1, 2)*std::pow(2*xi_1 + xi_2 - 1, 3) - 315*(2*xi_2 - 2)*std::pow(2*xi_1 + xi_2 - 1, 4) + 1386*std::pow(2*xi_1 + xi_2 - 1, 5)) ;
}
if( b == 43 ) {
  valDxi2 = -(17*xi_2 - 1)*(490*xi_1*std::pow(xi_2 - 1, 6) + 280*std::pow(xi_2 - 1, 7) - 945*std::pow(xi_2 - 1, 5)*std::pow(2*xi_1 + xi_2 - 1, 2) - 1575*std::pow(xi_2 - 1, 4)*std::pow(2*xi_1 + xi_2 - 1, 3) + 3465*std::pow(xi_2 - 1, 3)*std::pow(2*xi_1 + xi_2 - 1, 4) + 2079*std::pow(xi_2 - 1, 2)*std::pow(2*xi_1 + xi_2 - 1, 5) - 3003*(xi_2 - 1)*std::pow(2*xi_1 + xi_2 - 1, 6) - 429*std::pow(2*xi_1 + xi_2 - 1, 7))/(16*xi_2 - 16) - 17*(70*xi_1*std::pow(xi_2 - 1, 7) + 35*std::pow(xi_2 - 1, 8) - 315*std::pow(xi_2 - 1, 5)*std::pow(2*xi_1 + xi_2 - 1, 3) + 693*std::pow(xi_2 - 1, 3)*std::pow(2*xi_1 + xi_2 - 1, 5) - 429*(xi_2 - 1)*std::pow(2*xi_1 + xi_2 - 1, 7))/(16*xi_2 - 16) + 16*(17*xi_2 - 1)*(70*xi_1*std::pow(xi_2 - 1, 7) + 35*std::pow(xi_2 - 1, 8) - 315*std::pow(xi_2 - 1, 5)*std::pow(2*xi_1 + xi_2 - 1, 3) + 693*std::pow(xi_2 - 1, 3)*std::pow(2*xi_1 + xi_2 - 1, 5) - 429*(xi_2 - 1)*std::pow(2*xi_1 + xi_2 - 1, 7))/std::pow(16*xi_2 - 16, 2) ;
}
if( b == 44 ) {
  valDxi2 = (35.0L/16.0L)*std::pow(xi_2 - 1, 7) - 315.0L/32.0L*std::pow(xi_2 - 1, 6)*(4*xi_1 + 2*xi_2 - 2) - 945.0L/16.0L*std::pow(xi_2 - 1, 5)*std::pow(2*xi_1 + xi_2 - 1, 2) + (3465.0L/16.0L)*std::pow(xi_2 - 1, 4)*std::pow(2*xi_1 + xi_2 - 1, 3) + (3465.0L/16.0L)*std::pow(xi_2 - 1, 3)*std::pow(2*xi_1 + xi_2 - 1, 4) - 9009.0L/16.0L*std::pow(xi_2 - 1, 2)*std::pow(2*xi_1 + xi_2 - 1, 5) - 3003.0L/32.0L*(2*xi_2 - 2)*std::pow(2*xi_1 + xi_2 - 1, 6) + (6435.0L/16.0L)*std::pow(2*xi_1 + xi_2 - 1, 7) ;
}
if( b == 45 ) {
  valDxi2 = 831402*std::pow(xi_2, 8) - 3150576*std::pow(xi_2, 7) + 4900896*std::pow(xi_2, 6) - 4036032*std::pow(xi_2, 5) + 1891890*std::pow(xi_2, 4) - 504504*std::pow(xi_2, 3) + 72072*std::pow(xi_2, 2) - 4752*xi_2 + 99 ;
}
if( b == 46 ) {
  valDxi2 = -105.0L/8.0L*xi_2 + (37791.0L/128.0L)*std::pow(2*xi_2 - 1, 8) + (5967.0L/16.0L)*std::pow(2*xi_2 - 1, 7) - 7735.0L/32.0L*std::pow(2*xi_2 - 1, 6) - 6279.0L/16.0L*std::pow(2*xi_2 - 1, 5) + (1365.0L/64.0L)*std::pow(2*xi_2 - 1, 4) + (1729.0L/16.0L)*std::pow(2*xi_2 - 1, 3) + (273.0L/32.0L)*std::pow(2*xi_2 - 1, 2) + (1.0L/128.0L)*(2*xi_1 + xi_2 - 1)*(8736*xi_2 + 604656*std::pow(2*xi_2 - 1, 7) + 668304*std::pow(2*xi_2 - 1, 6) - 371280*std::pow(2*xi_2 - 1, 5) - 502320*std::pow(2*xi_2 - 1, 4) + 21840*std::pow(2*xi_2 - 1, 3) + 82992*std::pow(2*xi_2 - 1, 2) - 6048) + 791.0L/128.0L ;
}
if( b == 47 ) {
  valDxi2 = (1.0L/64.0L)*(-std::pow(xi_2 - 1, 2) + 3*std::pow(2*xi_1 + xi_2 - 1, 2))*(6552*xi_2 + 176358*std::pow(2*xi_2 - 1, 6) + 278460*std::pow(2*xi_2 - 1, 5) + 46410*std::pow(2*xi_2 - 1, 4) - 91000*std::pow(2*xi_2 - 1, 3) - 30030*std::pow(2*xi_2 - 1, 2) - 2366) + (1.0L/64.0L)*(12*xi_1 + 4*xi_2 - 4)*(910*xi_2 + 12597*std::pow(2*xi_2 - 1, 7) + 23205*std::pow(2*xi_2 - 1, 6) + 4641*std::pow(2*xi_2 - 1, 5) - 11375*std::pow(2*xi_2 - 1, 4) - 5005*std::pow(2*xi_2 - 1, 3) + 819*std::pow(2*xi_2 - 1, 2) - 448) ;
}
if( b == 48 ) {
  valDxi2 = -(6*xi_1*std::pow(xi_2 - 1, 3) + 3*std::pow(xi_2 - 1, 4) - 5*(xi_2 - 1)*std::pow(2*xi_1 + xi_2 - 1, 3))*(-15960*xi_2 + 81396*std::pow(2*xi_2 - 1, 5) + 149940*std::pow(2*xi_2 - 1, 4) + 71400*std::pow(2*xi_2 - 1, 3) - 5880*std::pow(2*xi_2 - 1, 2) + 7392)/(32*xi_2 - 32) - (18*xi_1*std::pow(xi_2 - 1, 2) + 12*std::pow(xi_2 - 1, 3) - 15*(xi_2 - 1)*std::pow(2*xi_1 + xi_2 - 1, 2) - 5*std::pow(2*xi_1 + xi_2 - 1, 3))*(-588*xi_2 + 6783*std::pow(2*xi_2 - 1, 6) + 14994*std::pow(2*xi_2 - 1, 5) + 8925*std::pow(2*xi_2 - 1, 4) - 980*std::pow(2*xi_2 - 1, 3) - 1995*std::pow(2*xi_2 - 1, 2) + 317)/(32*xi_2 - 32) + 32*(6*xi_1*std::pow(xi_2 - 1, 3) + 3*std::pow(xi_2 - 1, 4) - 5*(xi_2 - 1)*std::pow(2*xi_1 + xi_2 - 1, 3))*(-588*xi_2 + 6783*std::pow(2*xi_2 - 1, 6) + 14994*std::pow(2*xi_2 - 1, 5) + 8925*std::pow(2*xi_2 - 1, 4) - 980*std::pow(2*xi_2 - 1, 3) - 1995*std::pow(2*xi_2 - 1, 2) + 317)/std::pow(32*xi_2 - 32, 2) ;
}
if( b == 49 ) {
  valDxi2 = (1.0L/64.0L)*(3*std::pow(xi_2 - 1, 4) - 30*std::pow(xi_2 - 1, 2)*std::pow(2*xi_1 + xi_2 - 1, 2) + 35*std::pow(2*xi_1 + xi_2 - 1, 4))*(9360*xi_2 + 29070*std::pow(2*xi_2 - 1, 4) + 55080*std::pow(2*xi_2 - 1, 3) + 31620*std::pow(2*xi_2 - 1, 2) - 5010) + (1.0L/64.0L)*(12*std::pow(xi_2 - 1, 3) - 30*std::pow(xi_2 - 1, 2)*(4*xi_1 + 2*xi_2 - 2) - 30*(2*xi_2 - 2)*std::pow(2*xi_1 + xi_2 - 1, 2) + 140*std::pow(2*xi_1 + xi_2 - 1, 3))*(-330*xi_2 + 2907*std::pow(2*xi_2 - 1, 5) + 6885*std::pow(2*xi_2 - 1, 4) + 5270*std::pow(2*xi_2 - 1, 3) + 1170*std::pow(2*xi_2 - 1, 2) + 114) ;
}
if( b == 50 ) {
  valDxi2 = (1.0L/32.0L)*(13872*xi_2 + 7752*std::pow(2*xi_2 - 1, 3) + 13464*std::pow(2*xi_2 - 1, 2) - 5968)*(30*xi_1*std::pow(xi_2 - 1, 5) + 63*(-xi_2 + 1)*std::pow(-2*xi_1 - xi_2 + 1, 5) + 15*std::pow(xi_2 - 1, 6) + 70*std::pow(xi_2 - 1, 3)*std::pow(-2*xi_1 - xi_2 + 1, 3))/(xi_2 - 1) + (1.0L/32.0L)*(968*xi_2 + 969*std::pow(2*xi_2 - 1, 4) + 2244*std::pow(2*xi_2 - 1, 3) + 1734*std::pow(2*xi_2 - 1, 2) - 455)*(150*xi_1*std::pow(xi_2 - 1, 4) - 315*(-xi_2 + 1)*std::pow(-2*xi_1 - xi_2 + 1, 4) + 90*std::pow(xi_2 - 1, 5) - 210*std::pow(xi_2 - 1, 3)*std::pow(-2*xi_1 - xi_2 + 1, 2) + 210*std::pow(xi_2 - 1, 2)*std::pow(-2*xi_1 - xi_2 + 1, 3) - 63*std::pow(-2*xi_1 - xi_2 + 1, 5))/(xi_2 - 1) - 1.0L/32.0L*(30*xi_1*std::pow(xi_2 - 1, 5) + 63*(-xi_2 + 1)*std::pow(-2*xi_1 - xi_2 + 1, 5) + 15*std::pow(xi_2 - 1, 6) + 70*std::pow(xi_2 - 1, 3)*std::pow(-2*xi_1 - xi_2 + 1, 3))*(968*xi_2 + 969*std::pow(2*xi_2 - 1, 4) + 2244*std::pow(2*xi_2 - 1, 3) + 1734*std::pow(2*xi_2 - 1, 2) - 455)/std::pow(xi_2 - 1, 2) ;
}
if( b == 51 ) {
  valDxi2 = (1.0L/128.0L)*(15912*xi_2 + 5814*std::pow(2*xi_2 - 1, 2) - 5406)*(-5*std::pow(xi_2 - 1, 6) + 105*std::pow(xi_2 - 1, 4)*std::pow(2*xi_1 + xi_2 - 1, 2) - 315*std::pow(xi_2 - 1, 2)*std::pow(2*xi_1 + xi_2 - 1, 4) + 231*std::pow(2*xi_1 + xi_2 - 1, 6)) + (1.0L/128.0L)*(2550*xi_2 + 969*std::pow(2*xi_2 - 1, 3) + 1989*std::pow(2*xi_2 - 1, 2) - 1028)*(-30*std::pow(xi_2 - 1, 5) + 105*std::pow(xi_2 - 1, 4)*(4*xi_1 + 2*xi_2 - 2) + 420*std::pow(xi_2 - 1, 3)*std::pow(2*xi_1 + xi_2 - 1, 2) - 1260*std::pow(xi_2 - 1, 2)*std::pow(2*xi_1 + xi_2 - 1, 3) - 315*(2*xi_2 - 2)*std::pow(2*xi_1 + xi_2 - 1, 4) + 1386*std::pow(2*xi_1 + xi_2 - 1, 5)) ;
}
if( b == 52 ) {
  valDxi2 = -(1368*xi_2 - 144)*(70*xi_1*std::pow(xi_2 - 1, 7) + 35*std::pow(xi_2 - 1, 8) - 315*std::pow(xi_2 - 1, 5)*std::pow(2*xi_1 + xi_2 - 1, 3) + 693*std::pow(xi_2 - 1, 3)*std::pow(2*xi_1 + xi_2 - 1, 5) - 429*(xi_2 - 1)*std::pow(2*xi_1 + xi_2 - 1, 7))/(64*xi_2 - 64) - (540*xi_2 + 171*std::pow(2*xi_2 - 1, 2) - 167)*(490*xi_1*std::pow(xi_2 - 1, 6) + 280*std::pow(xi_2 - 1, 7) - 945*std::pow(xi_2 - 1, 5)*std::pow(2*xi_1 + xi_2 - 1, 2) - 1575*std::pow(xi_2 - 1, 4)*std::pow(2*xi_1 + xi_2 - 1, 3) + 3465*std::pow(xi_2 - 1, 3)*std::pow(2*xi_1 + xi_2 - 1, 4) + 2079*std::pow(xi_2 - 1, 2)*std::pow(2*xi_1 + xi_2 - 1, 5) - 3003*(xi_2 - 1)*std::pow(2*xi_1 + xi_2 - 1, 6) - 429*std::pow(2*xi_1 + xi_2 - 1, 7))/(64*xi_2 - 64) + 64*(540*xi_2 + 171*std::pow(2*xi_2 - 1, 2) - 167)*(70*xi_1*std::pow(xi_2 - 1, 7) + 35*std::pow(xi_2 - 1, 8) - 315*std::pow(xi_2 - 1, 5)*std::pow(2*xi_1 + xi_2 - 1, 3) + 693*std::pow(xi_2 - 1, 3)*std::pow(2*xi_1 + xi_2 - 1, 5) - 429*(xi_2 - 1)*std::pow(2*xi_1 + xi_2 - 1, 7))/std::pow(64*xi_2 - 64, 2) ;
}
if( b == 53 ) {
  valDxi2 = (665.0L/128.0L)*std::pow(xi_2 - 1, 8) - 5985.0L/32.0L*std::pow(xi_2 - 1, 6)*std::pow(2*xi_1 + xi_2 - 1, 2) + (65835.0L/64.0L)*std::pow(xi_2 - 1, 4)*std::pow(2*xi_1 + xi_2 - 1, 4) - 57057.0L/32.0L*std::pow(xi_2 - 1, 2)*std::pow(2*xi_1 + xi_2 - 1, 6) + (1.0L/128.0L)*(19*xi_2 - 1)*(280*std::pow(xi_2 - 1, 7) - 1260*std::pow(xi_2 - 1, 6)*(4*xi_1 + 2*xi_2 - 2) - 7560*std::pow(xi_2 - 1, 5)*std::pow(2*xi_1 + xi_2 - 1, 2) + 27720*std::pow(xi_2 - 1, 4)*std::pow(2*xi_1 + xi_2 - 1, 3) + 27720*std::pow(xi_2 - 1, 3)*std::pow(2*xi_1 + xi_2 - 1, 4) - 72072*std::pow(xi_2 - 1, 2)*std::pow(2*xi_1 + xi_2 - 1, 5) - 12012*(2*xi_2 - 2)*std::pow(2*xi_1 + xi_2 - 1, 6) + 51480*std::pow(2*xi_1 + xi_2 - 1, 7)) + (122265.0L/128.0L)*std::pow(2*xi_1 + xi_2 - 1, 8) ;
}
if( b == 54 ) {
  valDxi2 = (1.0L/128.0L)*(5670*xi_1*std::pow(xi_2 - 1, 8) - 90090*std::pow(-xi_2 + 1, 5)*std::pow(-2*xi_1 - xi_2 + 1, 4) - 90090*std::pow(-xi_2 + 1, 4)*std::pow(-2*xi_1 - xi_2 + 1, 5) - 109395*(-xi_2 + 1)*std::pow(-2*xi_1 - xi_2 + 1, 8) + 3150*std::pow(xi_2 - 1, 9) - 13860*std::pow(xi_2 - 1, 7)*std::pow(-2*xi_1 - xi_2 + 1, 2) + 32340*std::pow(xi_2 - 1, 6)*std::pow(-2*xi_1 - xi_2 + 1, 3) - 180180*std::pow(xi_2 - 1, 3)*std::pow(-2*xi_1 - xi_2 + 1, 6) + 77220*std::pow(xi_2 - 1, 2)*std::pow(-2*xi_1 - xi_2 + 1, 7) - 12155*std::pow(-2*xi_1 - xi_2 + 1, 9))/(xi_2 - 1) - 1.0L/128.0L*(630*xi_1*std::pow(xi_2 - 1, 9) + 18018*std::pow(-xi_2 + 1, 5)*std::pow(-2*xi_1 - xi_2 + 1, 5) + 12155*(-xi_2 + 1)*std::pow(-2*xi_1 - xi_2 + 1, 9) + 315*std::pow(xi_2 - 1, 10) + 4620*std::pow(xi_2 - 1, 7)*std::pow(-2*xi_1 - xi_2 + 1, 3) + 25740*std::pow(xi_2 - 1, 3)*std::pow(-2*xi_1 - xi_2 + 1, 7))/std::pow(xi_2 - 1, 2) ;
}
