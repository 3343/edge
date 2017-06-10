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
 * Generated triangular basis.
 **/

assert( b < 55 );
if( b == 0 ) {
  val = 1 ;
}
if( b == 1 ) {
  val = 3*xi_2 - 1 ;
}
if( b == 2 ) {
  val = 2*xi_1 + xi_2 - 1 ;
}
if( b == 3 ) {
  val = 10*std::pow(xi_2, 2) - 8*xi_2 + 1 ;
}
if( b == 4 ) {
  val = (5*xi_2 - 1)*(2*xi_1 + xi_2 - 1) ;
}
if( b == 5 ) {
  val = -1.0L/2.0L*std::pow(xi_2 - 1, 2) + (3.0L/2.0L)*std::pow(2*xi_1 + xi_2 - 1, 2) ;
}
if( b == 6 ) {
  val = 35*std::pow(xi_2, 3) - 45*std::pow(xi_2, 2) + 15*xi_2 - 1 ;
}
if( b == 7 ) {
  val = (1.0L/4.0L)*(2*xi_1 + xi_2 - 1)*(36*xi_2 + 21*std::pow(2*xi_2 - 1, 2) - 17) ;
}
if( b == 8 ) {
  val = (1.0L/2.0L)*(7*xi_2 - 1)*(-std::pow(xi_2 - 1, 2) + 3*std::pow(2*xi_1 + xi_2 - 1, 2)) ;
}
if( b == 9 ) {
  val = -(6*xi_1*std::pow(xi_2 - 1, 3) + 3*std::pow(xi_2 - 1, 4) - 5*(xi_2 - 1)*std::pow(2*xi_1 + xi_2 - 1, 3))/(2*xi_2 - 2) ;
}
if( b == 10 ) {
  val = 126*std::pow(xi_2, 4) - 224*std::pow(xi_2, 3) + 126*std::pow(xi_2, 2) - 24*xi_2 + 1 ;
}
if( b == 11 ) {
  val = (1.0L/2.0L)*(2*xi_1 + xi_2 - 1)*(21*std::pow(2*xi_2 - 1, 3) + 21*std::pow(2*xi_2 - 1, 2) - 2) ;
}
if( b == 12 ) {
  val = (1.0L/2.0L)*(-std::pow(xi_2 - 1, 2) + 3*std::pow(2*xi_1 + xi_2 - 1, 2))*(20*xi_2 + 9*std::pow(2*xi_2 - 1, 2) - 8) ;
}
if( b == 13 ) {
  val = -(9*xi_2 - 1)*(6*xi_1*std::pow(xi_2 - 1, 3) + 3*std::pow(xi_2 - 1, 4) - 5*(xi_2 - 1)*std::pow(2*xi_1 + xi_2 - 1, 3))/(2*xi_2 - 2) ;
}
if( b == 14 ) {
  val = (3.0L/8.0L)*std::pow(xi_2 - 1, 4) - 15.0L/4.0L*std::pow(xi_2 - 1, 2)*std::pow(2*xi_1 + xi_2 - 1, 2) + (35.0L/8.0L)*std::pow(2*xi_1 + xi_2 - 1, 4) ;
}
if( b == 15 ) {
  val = 462*std::pow(xi_2, 5) - 1050*std::pow(xi_2, 4) + 840*std::pow(xi_2, 3) - 280*std::pow(xi_2, 2) + 35*xi_2 - 1 ;
}
if( b == 16 ) {
  val = (1.0L/8.0L)*(2*xi_1 + xi_2 - 1)*(-88*xi_2 + 165*std::pow(2*xi_2 - 1, 4) + 180*std::pow(2*xi_2 - 1, 3) - 18*std::pow(2*xi_2 - 1, 2) + 41) ;
}
if( b == 17 ) {
  val = (1.0L/16.0L)*(-std::pow(xi_2 - 1, 2) + 3*std::pow(2*xi_1 + xi_2 - 1, 2))*(126*xi_2 + 165*std::pow(2*xi_2 - 1, 3) + 225*std::pow(2*xi_2 - 1, 2) - 68) ;
}
if( b == 18 ) {
  val = -(140*xi_2 + 55*std::pow(2*xi_2 - 1, 2) - 51)*(6*xi_1*std::pow(xi_2 - 1, 3) + 3*std::pow(xi_2 - 1, 4) - 5*(xi_2 - 1)*std::pow(2*xi_1 + xi_2 - 1, 3))/(8*xi_2 - 8) ;
}
if( b == 19 ) {
  val = (1.0L/8.0L)*(11*xi_2 - 1)*(3*std::pow(xi_2 - 1, 4) - 30*std::pow(xi_2 - 1, 2)*std::pow(2*xi_1 + xi_2 - 1, 2) + 35*std::pow(2*xi_1 + xi_2 - 1, 4)) ;
}
if( b == 20 ) {
  val = (1.0L/8.0L)*(30*xi_1*std::pow(xi_2 - 1, 5) + 63*(-xi_2 + 1)*std::pow(-2*xi_1 - xi_2 + 1, 5) + 15*std::pow(xi_2 - 1, 6) + 70*std::pow(xi_2 - 1, 3)*std::pow(-2*xi_1 - xi_2 + 1, 3))/(xi_2 - 1) ;
}
if( b == 21 ) {
  val = 1716*std::pow(xi_2, 6) - 4752*std::pow(xi_2, 5) + 4950*std::pow(xi_2, 4) - 2400*std::pow(xi_2, 3) + 540*std::pow(xi_2, 2) - 48*xi_2 + 1 ;
}
if( b == 22 ) {
  val = -1.0L/32.0L*(2*xi_1 + xi_2 - 1)*(90*xi_2 - 1287*std::pow(2*xi_2 - 1, 5) - 1485*std::pow(2*xi_2 - 1, 4) + 330*std::pow(2*xi_2 - 1, 3) + 630*std::pow(2*xi_2 - 1, 2) - 70) ;
}
if( b == 23 ) {
  val = (1.0L/32.0L)*(-std::pow(xi_2 - 1, 2) + 3*std::pow(2*xi_1 + xi_2 - 1, 2))*(-200*xi_2 + 715*std::pow(2*xi_2 - 1, 4) + 1100*std::pow(2*xi_2 - 1, 3) + 330*std::pow(2*xi_2 - 1, 2) + 71) ;
}
if( b == 24 ) {
  val = -(6*xi_1*std::pow(xi_2 - 1, 3) + 3*std::pow(xi_2 - 1, 4) - 5*(xi_2 - 1)*std::pow(2*xi_1 + xi_2 - 1, 3))*(198*xi_2 + 143*std::pow(2*xi_2 - 1, 3) + 231*std::pow(2*xi_2 - 1, 2) - 92)/(8*xi_2 - 8) ;
}
if( b == 25 ) {
  val = (1.0L/16.0L)*(108*xi_2 + 39*std::pow(2*xi_2 - 1, 2) - 37)*(3*std::pow(xi_2 - 1, 4) - 30*std::pow(xi_2 - 1, 2)*std::pow(2*xi_1 + xi_2 - 1, 2) + 35*std::pow(2*xi_1 + xi_2 - 1, 4)) ;
}
if( b == 26 ) {
  val = (1.0L/8.0L)*(13*xi_2 - 1)*(30*xi_1*std::pow(xi_2 - 1, 5) + 63*(-xi_2 + 1)*std::pow(-2*xi_1 - xi_2 + 1, 5) + 15*std::pow(xi_2 - 1, 6) + 70*std::pow(xi_2 - 1, 3)*std::pow(-2*xi_1 - xi_2 + 1, 3))/(xi_2 - 1) ;
}
if( b == 27 ) {
  val = -5.0L/16.0L*std::pow(xi_2 - 1, 6) + (105.0L/16.0L)*std::pow(xi_2 - 1, 4)*std::pow(2*xi_1 + xi_2 - 1, 2) - 315.0L/16.0L*std::pow(xi_2 - 1, 2)*std::pow(2*xi_1 + xi_2 - 1, 4) + (231.0L/16.0L)*std::pow(2*xi_1 + xi_2 - 1, 6) ;
}
if( b == 28 ) {
  val = 6435*std::pow(xi_2, 7) - 21021*std::pow(xi_2, 6) + 27027*std::pow(xi_2, 5) - 17325*std::pow(xi_2, 4) + 5775*std::pow(xi_2, 3) - 945*std::pow(xi_2, 2) + 63*xi_2 - 1 ;
}
if( b == 29 ) {
  val = -1.0L/64.0L*(2*xi_1 + xi_2 - 1)*(-780*xi_2 - 5005*std::pow(2*xi_2 - 1, 6) - 6006*std::pow(2*xi_2 - 1, 5) + 2145*std::pow(2*xi_2 - 1, 4) + 3740*std::pow(2*xi_2 - 1, 3) + 165*std::pow(2*xi_2 - 1, 2) + 365) ;
}
if( b == 30 ) {
  val = (1.0L/64.0L)*(-std::pow(xi_2 - 1, 2) + 3*std::pow(2*xi_1 + xi_2 - 1, 2))*(-770*xi_2 + 3003*std::pow(2*xi_2 - 1, 5) + 5005*std::pow(2*xi_2 - 1, 4) + 1430*std::pow(2*xi_2 - 1, 3) - 990*std::pow(2*xi_2 - 1, 2) + 386) ;
}
if( b == 31 ) {
  val = -(6*xi_1*std::pow(xi_2 - 1, 3) + 3*std::pow(xi_2 - 1, 4) - 5*(xi_2 - 1)*std::pow(2*xi_1 + xi_2 - 1, 3))*(168*xi_2 + 1365*std::pow(2*xi_2 - 1, 4) + 2548*std::pow(2*xi_2 - 1, 3) + 1326*std::pow(2*xi_2 - 1, 2) - 127)/(32*xi_2 - 32) ;
}
if( b == 32 ) {
  val = (1.0L/64.0L)*(3*std::pow(xi_2 - 1, 4) - 30*std::pow(xi_2 - 1, 2)*std::pow(2*xi_1 + xi_2 - 1, 2) + 35*std::pow(2*xi_1 + xi_2 - 1, 4))*(858*xi_2 + 455*std::pow(2*xi_2 - 1, 3) + 819*std::pow(2*xi_2 - 1, 2) - 372) ;
}
if( b == 33 ) {
  val = (1.0L/32.0L)*(308*xi_2 + 105*std::pow(2*xi_2 - 1, 2) - 101)*(30*xi_1*std::pow(xi_2 - 1, 5) + 63*(-xi_2 + 1)*std::pow(-2*xi_1 - xi_2 + 1, 5) + 15*std::pow(xi_2 - 1, 6) + 70*std::pow(xi_2 - 1, 3)*std::pow(-2*xi_1 - xi_2 + 1, 3))/(xi_2 - 1) ;
}
if( b == 34 ) {
  val = (1.0L/16.0L)*(15*xi_2 - 1)*(-5*std::pow(xi_2 - 1, 6) + 105*std::pow(xi_2 - 1, 4)*std::pow(2*xi_1 + xi_2 - 1, 2) - 315*std::pow(xi_2 - 1, 2)*std::pow(2*xi_1 + xi_2 - 1, 4) + 231*std::pow(2*xi_1 + xi_2 - 1, 6)) ;
}
if( b == 35 ) {
  val = (1.0L/16.0L)*(-70*xi_1*std::pow(xi_2 - 1, 7) - 429*(-xi_2 + 1)*std::pow(2*xi_1 + xi_2 - 1, 7) - 35*std::pow(xi_2 - 1, 8) + 315*std::pow(xi_2 - 1, 5)*std::pow(2*xi_1 + xi_2 - 1, 3) - 693*std::pow(xi_2 - 1, 3)*std::pow(2*xi_1 + xi_2 - 1, 5))/(xi_2 - 1) ;
}
if( b == 36 ) {
  val = 24310*std::pow(xi_2, 8) - 91520*std::pow(xi_2, 7) + 140140*std::pow(xi_2, 6) - 112112*std::pow(xi_2, 5) + 50050*std::pow(xi_2, 4) - 12320*std::pow(xi_2, 3) + 1540*std::pow(xi_2, 2) - 80*xi_2 + 1 ;
}
if( b == 37 ) {
  val = (1.0L/32.0L)*(2*xi_1 + xi_2 - 1)*(154*xi_2 + 4862*std::pow(2*xi_2 - 1, 7) + 6006*std::pow(2*xi_2 - 1, 6) - 3003*std::pow(2*xi_2 - 1, 5) - 5005*std::pow(2*xi_2 - 1, 4) + 924*std::pow(2*xi_2 - 1, 2) - 98) ;
}
if( b == 38 ) {
  val = (1.0L/32.0L)*(-std::pow(xi_2 - 1, 2) + 3*std::pow(2*xi_1 + xi_2 - 1, 2))*(96*xi_2 + 3094*std::pow(2*xi_2 - 1, 6) + 5460*std::pow(2*xi_2 - 1, 5) + 1365*std::pow(2*xi_2 - 1, 4) - 1820*std::pow(2*xi_2 - 1, 3) - 780*std::pow(2*xi_2 - 1, 2) - 23) ;
}
if( b == 39 ) {
  val = -(6*xi_1*std::pow(xi_2 - 1, 3) + 3*std::pow(xi_2 - 1, 4) - 5*(xi_2 - 1)*std::pow(2*xi_1 + xi_2 - 1, 3))*(-390*xi_2 + 1547*std::pow(2*xi_2 - 1, 5) + 3185*std::pow(2*xi_2 - 1, 4) + 1820*std::pow(2*xi_2 - 1, 3) + 174)/(16*xi_2 - 16) ;
}
if( b == 40 ) {
  val = (1.0L/32.0L)*(3*std::pow(xi_2 - 1, 4) - 30*std::pow(xi_2 - 1, 2)*std::pow(2*xi_1 + xi_2 - 1, 2) + 35*std::pow(2*xi_1 + xi_2 - 1, 4))*(336*xi_2 + 595*std::pow(2*xi_2 - 1, 4) + 1260*std::pow(2*xi_2 - 1, 3) + 840*std::pow(2*xi_2 - 1, 2) - 171) ;
}
if( b == 41 ) {
  val = (1.0L/8.0L)*(195*xi_2 + 85*std::pow(2*xi_2 - 1, 3) + 165*std::pow(2*xi_2 - 1, 2) - 81)*(30*xi_1*std::pow(xi_2 - 1, 5) + 63*(-xi_2 + 1)*std::pow(-2*xi_1 - xi_2 + 1, 5) + 15*std::pow(xi_2 - 1, 6) + 70*std::pow(xi_2 - 1, 3)*std::pow(-2*xi_1 - xi_2 + 1, 3))/(xi_2 - 1) ;
}
if( b == 42 ) {
  val = (1.0L/16.0L)*(104*xi_2 + 34*std::pow(2*xi_2 - 1, 2) - 33)*(-5*std::pow(xi_2 - 1, 6) + 105*std::pow(xi_2 - 1, 4)*std::pow(2*xi_1 + xi_2 - 1, 2) - 315*std::pow(xi_2 - 1, 2)*std::pow(2*xi_1 + xi_2 - 1, 4) + 231*std::pow(2*xi_1 + xi_2 - 1, 6)) ;
}
if( b == 43 ) {
  val = -(17*xi_2 - 1)*(70*xi_1*std::pow(xi_2 - 1, 7) + 35*std::pow(xi_2 - 1, 8) - 315*std::pow(xi_2 - 1, 5)*std::pow(2*xi_1 + xi_2 - 1, 3) + 693*std::pow(xi_2 - 1, 3)*std::pow(2*xi_1 + xi_2 - 1, 5) - 429*(xi_2 - 1)*std::pow(2*xi_1 + xi_2 - 1, 7))/(16*xi_2 - 16) ;
}
if( b == 44 ) {
  val = (35.0L/128.0L)*std::pow(xi_2 - 1, 8) - 315.0L/32.0L*std::pow(xi_2 - 1, 6)*std::pow(2*xi_1 + xi_2 - 1, 2) + (3465.0L/64.0L)*std::pow(xi_2 - 1, 4)*std::pow(2*xi_1 + xi_2 - 1, 4) - 3003.0L/32.0L*std::pow(xi_2 - 1, 2)*std::pow(2*xi_1 + xi_2 - 1, 6) + (6435.0L/128.0L)*std::pow(2*xi_1 + xi_2 - 1, 8) ;
}
if( b == 45 ) {
  val = 92378*std::pow(xi_2, 9) - 393822*std::pow(xi_2, 8) + 700128*std::pow(xi_2, 7) - 672672*std::pow(xi_2, 6) + 378378*std::pow(xi_2, 5) - 126126*std::pow(xi_2, 4) + 24024*std::pow(xi_2, 3) - 2376*std::pow(xi_2, 2) + 99*xi_2 - 1 ;
}
if( b == 46 ) {
  val = (1.0L/128.0L)*(2*xi_1 + xi_2 - 1)*(-1680*xi_2 + 37791*std::pow(2*xi_2 - 1, 8) + 47736*std::pow(2*xi_2 - 1, 7) - 30940*std::pow(2*xi_2 - 1, 6) - 50232*std::pow(2*xi_2 - 1, 5) + 2730*std::pow(2*xi_2 - 1, 4) + 13832*std::pow(2*xi_2 - 1, 3) + 1092*std::pow(2*xi_2 - 1, 2) + 791) ;
}
if( b == 47 ) {
  val = (1.0L/64.0L)*(-std::pow(xi_2 - 1, 2) + 3*std::pow(2*xi_1 + xi_2 - 1, 2))*(910*xi_2 + 12597*std::pow(2*xi_2 - 1, 7) + 23205*std::pow(2*xi_2 - 1, 6) + 4641*std::pow(2*xi_2 - 1, 5) - 11375*std::pow(2*xi_2 - 1, 4) - 5005*std::pow(2*xi_2 - 1, 3) + 819*std::pow(2*xi_2 - 1, 2) - 448) ;
}
if( b == 48 ) {
  val = -(6*xi_1*std::pow(xi_2 - 1, 3) + 3*std::pow(xi_2 - 1, 4) - 5*(xi_2 - 1)*std::pow(2*xi_1 + xi_2 - 1, 3))*(-588*xi_2 + 6783*std::pow(2*xi_2 - 1, 6) + 14994*std::pow(2*xi_2 - 1, 5) + 8925*std::pow(2*xi_2 - 1, 4) - 980*std::pow(2*xi_2 - 1, 3) - 1995*std::pow(2*xi_2 - 1, 2) + 317)/(32*xi_2 - 32) ;
}
if( b == 49 ) {
  val = (1.0L/64.0L)*(3*std::pow(xi_2 - 1, 4) - 30*std::pow(xi_2 - 1, 2)*std::pow(2*xi_1 + xi_2 - 1, 2) + 35*std::pow(2*xi_1 + xi_2 - 1, 4))*(-330*xi_2 + 2907*std::pow(2*xi_2 - 1, 5) + 6885*std::pow(2*xi_2 - 1, 4) + 5270*std::pow(2*xi_2 - 1, 3) + 1170*std::pow(2*xi_2 - 1, 2) + 114) ;
}
if( b == 50 ) {
  val = (1.0L/32.0L)*(30*xi_1*std::pow(xi_2 - 1, 5) + 63*(-xi_2 + 1)*std::pow(-2*xi_1 - xi_2 + 1, 5) + 15*std::pow(xi_2 - 1, 6) + 70*std::pow(xi_2 - 1, 3)*std::pow(-2*xi_1 - xi_2 + 1, 3))*(968*xi_2 + 969*std::pow(2*xi_2 - 1, 4) + 2244*std::pow(2*xi_2 - 1, 3) + 1734*std::pow(2*xi_2 - 1, 2) - 455)/(xi_2 - 1) ;
}
if( b == 51 ) {
  val = (1.0L/128.0L)*(2550*xi_2 + 969*std::pow(2*xi_2 - 1, 3) + 1989*std::pow(2*xi_2 - 1, 2) - 1028)*(-5*std::pow(xi_2 - 1, 6) + 105*std::pow(xi_2 - 1, 4)*std::pow(2*xi_1 + xi_2 - 1, 2) - 315*std::pow(xi_2 - 1, 2)*std::pow(2*xi_1 + xi_2 - 1, 4) + 231*std::pow(2*xi_1 + xi_2 - 1, 6)) ;
}
if( b == 52 ) {
  val = -(540*xi_2 + 171*std::pow(2*xi_2 - 1, 2) - 167)*(70*xi_1*std::pow(xi_2 - 1, 7) + 35*std::pow(xi_2 - 1, 8) - 315*std::pow(xi_2 - 1, 5)*std::pow(2*xi_1 + xi_2 - 1, 3) + 693*std::pow(xi_2 - 1, 3)*std::pow(2*xi_1 + xi_2 - 1, 5) - 429*(xi_2 - 1)*std::pow(2*xi_1 + xi_2 - 1, 7))/(64*xi_2 - 64) ;
}
if( b == 53 ) {
  val = (1.0L/128.0L)*(19*xi_2 - 1)*(35*std::pow(xi_2 - 1, 8) - 1260*std::pow(xi_2 - 1, 6)*std::pow(2*xi_1 + xi_2 - 1, 2) + 6930*std::pow(xi_2 - 1, 4)*std::pow(2*xi_1 + xi_2 - 1, 4) - 12012*std::pow(xi_2 - 1, 2)*std::pow(2*xi_1 + xi_2 - 1, 6) + 6435*std::pow(2*xi_1 + xi_2 - 1, 8)) ;
}
if( b == 54 ) {
  val = (1.0L/128.0L)*(630*xi_1*std::pow(xi_2 - 1, 9) + 18018*std::pow(-xi_2 + 1, 5)*std::pow(-2*xi_1 - xi_2 + 1, 5) + 12155*(-xi_2 + 1)*std::pow(-2*xi_1 - xi_2 + 1, 9) + 315*std::pow(xi_2 - 1, 10) + 4620*std::pow(xi_2 - 1, 7)*std::pow(-2*xi_1 - xi_2 + 1, 3) + 25740*std::pow(xi_2 - 1, 3)*std::pow(-2*xi_1 - xi_2 + 1, 7))/(xi_2 - 1) ;
}
