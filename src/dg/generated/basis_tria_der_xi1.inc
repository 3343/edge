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
 * Generated derivatives of triangular basis w.r.t. to xi1.
 **/

assert( b < 55 );
if( b == 0 ) {
  valDxi1 = 0 ;
}
if( b == 1 ) {
  valDxi1 = 0 ;
}
if( b == 2 ) {
  valDxi1 = 2 ;
}
if( b == 3 ) {
  valDxi1 = 0 ;
}
if( b == 4 ) {
  valDxi1 = 10*xi_2 - 2 ;
}
if( b == 5 ) {
  valDxi1 = 12*xi_1 + 6*xi_2 - 6 ;
}
if( b == 6 ) {
  valDxi1 = 0 ;
}
if( b == 7 ) {
  valDxi1 = 18*xi_2 + (21.0L/2.0L)*std::pow(2*xi_2 - 1, 2) - 17.0L/2.0L ;
}
if( b == 8 ) {
  valDxi1 = (1.0L/2.0L)*(7*xi_2 - 1)*(24*xi_1 + 12*xi_2 - 12) ;
}
if( b == 9 ) {
  valDxi1 = -(6*std::pow(xi_2 - 1, 3) - 30*(xi_2 - 1)*std::pow(2*xi_1 + xi_2 - 1, 2))/(2*xi_2 - 2) ;
}
if( b == 10 ) {
  valDxi1 = 0 ;
}
if( b == 11 ) {
  valDxi1 = 21*std::pow(2*xi_2 - 1, 3) + 21*std::pow(2*xi_2 - 1, 2) - 2 ;
}
if( b == 12 ) {
  valDxi1 = (1.0L/2.0L)*(24*xi_1 + 12*xi_2 - 12)*(20*xi_2 + 9*std::pow(2*xi_2 - 1, 2) - 8) ;
}
if( b == 13 ) {
  valDxi1 = -(9*xi_2 - 1)*(6*std::pow(xi_2 - 1, 3) - 30*(xi_2 - 1)*std::pow(2*xi_1 + xi_2 - 1, 2))/(2*xi_2 - 2) ;
}
if( b == 14 ) {
  valDxi1 = -15.0L/4.0L*std::pow(xi_2 - 1, 2)*(8*xi_1 + 4*xi_2 - 4) + 35*std::pow(2*xi_1 + xi_2 - 1, 3) ;
}
if( b == 15 ) {
  valDxi1 = 0 ;
}
if( b == 16 ) {
  valDxi1 = -22*xi_2 + (165.0L/4.0L)*std::pow(2*xi_2 - 1, 4) + 45*std::pow(2*xi_2 - 1, 3) - 9.0L/2.0L*std::pow(2*xi_2 - 1, 2) + 41.0L/4.0L ;
}
if( b == 17 ) {
  valDxi1 = (1.0L/16.0L)*(24*xi_1 + 12*xi_2 - 12)*(126*xi_2 + 165*std::pow(2*xi_2 - 1, 3) + 225*std::pow(2*xi_2 - 1, 2) - 68) ;
}
if( b == 18 ) {
  valDxi1 = -(6*std::pow(xi_2 - 1, 3) - 30*(xi_2 - 1)*std::pow(2*xi_1 + xi_2 - 1, 2))*(140*xi_2 + 55*std::pow(2*xi_2 - 1, 2) - 51)/(8*xi_2 - 8) ;
}
if( b == 19 ) {
  valDxi1 = (1.0L/8.0L)*(11*xi_2 - 1)*(-30*std::pow(xi_2 - 1, 2)*(8*xi_1 + 4*xi_2 - 4) + 280*std::pow(2*xi_1 + xi_2 - 1, 3)) ;
}
if( b == 20 ) {
  valDxi1 = (1.0L/8.0L)*(-630*(-xi_2 + 1)*std::pow(-2*xi_1 - xi_2 + 1, 4) + 30*std::pow(xi_2 - 1, 5) - 420*std::pow(xi_2 - 1, 3)*std::pow(-2*xi_1 - xi_2 + 1, 2))/(xi_2 - 1) ;
}
if( b == 21 ) {
  valDxi1 = 0 ;
}
if( b == 22 ) {
  valDxi1 = -45.0L/8.0L*xi_2 + (1287.0L/16.0L)*std::pow(2*xi_2 - 1, 5) + (1485.0L/16.0L)*std::pow(2*xi_2 - 1, 4) - 165.0L/8.0L*std::pow(2*xi_2 - 1, 3) - 315.0L/8.0L*std::pow(2*xi_2 - 1, 2) + 35.0L/8.0L ;
}
if( b == 23 ) {
  valDxi1 = (1.0L/32.0L)*(24*xi_1 + 12*xi_2 - 12)*(-200*xi_2 + 715*std::pow(2*xi_2 - 1, 4) + 1100*std::pow(2*xi_2 - 1, 3) + 330*std::pow(2*xi_2 - 1, 2) + 71) ;
}
if( b == 24 ) {
  valDxi1 = -(6*std::pow(xi_2 - 1, 3) - 30*(xi_2 - 1)*std::pow(2*xi_1 + xi_2 - 1, 2))*(198*xi_2 + 143*std::pow(2*xi_2 - 1, 3) + 231*std::pow(2*xi_2 - 1, 2) - 92)/(8*xi_2 - 8) ;
}
if( b == 25 ) {
  valDxi1 = (1.0L/16.0L)*(-30*std::pow(xi_2 - 1, 2)*(8*xi_1 + 4*xi_2 - 4) + 280*std::pow(2*xi_1 + xi_2 - 1, 3))*(108*xi_2 + 39*std::pow(2*xi_2 - 1, 2) - 37) ;
}
if( b == 26 ) {
  valDxi1 = (1.0L/8.0L)*(13*xi_2 - 1)*(-630*(-xi_2 + 1)*std::pow(-2*xi_1 - xi_2 + 1, 4) + 30*std::pow(xi_2 - 1, 5) - 420*std::pow(xi_2 - 1, 3)*std::pow(-2*xi_1 - xi_2 + 1, 2))/(xi_2 - 1) ;
}
if( b == 27 ) {
  valDxi1 = (105.0L/16.0L)*std::pow(xi_2 - 1, 4)*(8*xi_1 + 4*xi_2 - 4) - 315.0L/2.0L*std::pow(xi_2 - 1, 2)*std::pow(2*xi_1 + xi_2 - 1, 3) + (693.0L/4.0L)*std::pow(2*xi_1 + xi_2 - 1, 5) ;
}
if( b == 28 ) {
  valDxi1 = 0 ;
}
if( b == 29 ) {
  valDxi1 = (195.0L/8.0L)*xi_2 + (5005.0L/32.0L)*std::pow(2*xi_2 - 1, 6) + (3003.0L/16.0L)*std::pow(2*xi_2 - 1, 5) - 2145.0L/32.0L*std::pow(2*xi_2 - 1, 4) - 935.0L/8.0L*std::pow(2*xi_2 - 1, 3) - 165.0L/32.0L*std::pow(2*xi_2 - 1, 2) - 365.0L/32.0L ;
}
if( b == 30 ) {
  valDxi1 = (1.0L/64.0L)*(24*xi_1 + 12*xi_2 - 12)*(-770*xi_2 + 3003*std::pow(2*xi_2 - 1, 5) + 5005*std::pow(2*xi_2 - 1, 4) + 1430*std::pow(2*xi_2 - 1, 3) - 990*std::pow(2*xi_2 - 1, 2) + 386) ;
}
if( b == 31 ) {
  valDxi1 = -(6*std::pow(xi_2 - 1, 3) - 30*(xi_2 - 1)*std::pow(2*xi_1 + xi_2 - 1, 2))*(168*xi_2 + 1365*std::pow(2*xi_2 - 1, 4) + 2548*std::pow(2*xi_2 - 1, 3) + 1326*std::pow(2*xi_2 - 1, 2) - 127)/(32*xi_2 - 32) ;
}
if( b == 32 ) {
  valDxi1 = (1.0L/64.0L)*(-30*std::pow(xi_2 - 1, 2)*(8*xi_1 + 4*xi_2 - 4) + 280*std::pow(2*xi_1 + xi_2 - 1, 3))*(858*xi_2 + 455*std::pow(2*xi_2 - 1, 3) + 819*std::pow(2*xi_2 - 1, 2) - 372) ;
}
if( b == 33 ) {
  valDxi1 = (1.0L/32.0L)*(308*xi_2 + 105*std::pow(2*xi_2 - 1, 2) - 101)*(-630*(-xi_2 + 1)*std::pow(-2*xi_1 - xi_2 + 1, 4) + 30*std::pow(xi_2 - 1, 5) - 420*std::pow(xi_2 - 1, 3)*std::pow(-2*xi_1 - xi_2 + 1, 2))/(xi_2 - 1) ;
}
if( b == 34 ) {
  valDxi1 = (1.0L/16.0L)*(15*xi_2 - 1)*(105*std::pow(xi_2 - 1, 4)*(8*xi_1 + 4*xi_2 - 4) - 2520*std::pow(xi_2 - 1, 2)*std::pow(2*xi_1 + xi_2 - 1, 3) + 2772*std::pow(2*xi_1 + xi_2 - 1, 5)) ;
}
if( b == 35 ) {
  valDxi1 = (1.0L/16.0L)*(-6006*(-xi_2 + 1)*std::pow(2*xi_1 + xi_2 - 1, 6) - 70*std::pow(xi_2 - 1, 7) + 1890*std::pow(xi_2 - 1, 5)*std::pow(2*xi_1 + xi_2 - 1, 2) - 6930*std::pow(xi_2 - 1, 3)*std::pow(2*xi_1 + xi_2 - 1, 4))/(xi_2 - 1) ;
}
if( b == 36 ) {
  valDxi1 = 0 ;
}
if( b == 37 ) {
  valDxi1 = (77.0L/8.0L)*xi_2 + (2431.0L/8.0L)*std::pow(2*xi_2 - 1, 7) + (3003.0L/8.0L)*std::pow(2*xi_2 - 1, 6) - 3003.0L/16.0L*std::pow(2*xi_2 - 1, 5) - 5005.0L/16.0L*std::pow(2*xi_2 - 1, 4) + (231.0L/4.0L)*std::pow(2*xi_2 - 1, 2) - 49.0L/8.0L ;
}
if( b == 38 ) {
  valDxi1 = (1.0L/32.0L)*(24*xi_1 + 12*xi_2 - 12)*(96*xi_2 + 3094*std::pow(2*xi_2 - 1, 6) + 5460*std::pow(2*xi_2 - 1, 5) + 1365*std::pow(2*xi_2 - 1, 4) - 1820*std::pow(2*xi_2 - 1, 3) - 780*std::pow(2*xi_2 - 1, 2) - 23) ;
}
if( b == 39 ) {
  valDxi1 = -(6*std::pow(xi_2 - 1, 3) - 30*(xi_2 - 1)*std::pow(2*xi_1 + xi_2 - 1, 2))*(-390*xi_2 + 1547*std::pow(2*xi_2 - 1, 5) + 3185*std::pow(2*xi_2 - 1, 4) + 1820*std::pow(2*xi_2 - 1, 3) + 174)/(16*xi_2 - 16) ;
}
if( b == 40 ) {
  valDxi1 = (1.0L/32.0L)*(-30*std::pow(xi_2 - 1, 2)*(8*xi_1 + 4*xi_2 - 4) + 280*std::pow(2*xi_1 + xi_2 - 1, 3))*(336*xi_2 + 595*std::pow(2*xi_2 - 1, 4) + 1260*std::pow(2*xi_2 - 1, 3) + 840*std::pow(2*xi_2 - 1, 2) - 171) ;
}
if( b == 41 ) {
  valDxi1 = (1.0L/8.0L)*(-630*(-xi_2 + 1)*std::pow(-2*xi_1 - xi_2 + 1, 4) + 30*std::pow(xi_2 - 1, 5) - 420*std::pow(xi_2 - 1, 3)*std::pow(-2*xi_1 - xi_2 + 1, 2))*(195*xi_2 + 85*std::pow(2*xi_2 - 1, 3) + 165*std::pow(2*xi_2 - 1, 2) - 81)/(xi_2 - 1) ;
}
if( b == 42 ) {
  valDxi1 = (1.0L/16.0L)*(104*xi_2 + 34*std::pow(2*xi_2 - 1, 2) - 33)*(105*std::pow(xi_2 - 1, 4)*(8*xi_1 + 4*xi_2 - 4) - 2520*std::pow(xi_2 - 1, 2)*std::pow(2*xi_1 + xi_2 - 1, 3) + 2772*std::pow(2*xi_1 + xi_2 - 1, 5)) ;
}
if( b == 43 ) {
  valDxi1 = -(17*xi_2 - 1)*(70*std::pow(xi_2 - 1, 7) - 1890*std::pow(xi_2 - 1, 5)*std::pow(2*xi_1 + xi_2 - 1, 2) + 6930*std::pow(xi_2 - 1, 3)*std::pow(2*xi_1 + xi_2 - 1, 4) - 6006*(xi_2 - 1)*std::pow(2*xi_1 + xi_2 - 1, 6))/(16*xi_2 - 16) ;
}
if( b == 44 ) {
  valDxi1 = -315.0L/32.0L*std::pow(xi_2 - 1, 6)*(8*xi_1 + 4*xi_2 - 4) + (3465.0L/8.0L)*std::pow(xi_2 - 1, 4)*std::pow(2*xi_1 + xi_2 - 1, 3) - 9009.0L/8.0L*std::pow(xi_2 - 1, 2)*std::pow(2*xi_1 + xi_2 - 1, 5) + (6435.0L/8.0L)*std::pow(2*xi_1 + xi_2 - 1, 7) ;
}
if( b == 45 ) {
  valDxi1 = 0 ;
}
if( b == 46 ) {
  valDxi1 = -105.0L/4.0L*xi_2 + (37791.0L/64.0L)*std::pow(2*xi_2 - 1, 8) + (5967.0L/8.0L)*std::pow(2*xi_2 - 1, 7) - 7735.0L/16.0L*std::pow(2*xi_2 - 1, 6) - 6279.0L/8.0L*std::pow(2*xi_2 - 1, 5) + (1365.0L/32.0L)*std::pow(2*xi_2 - 1, 4) + (1729.0L/8.0L)*std::pow(2*xi_2 - 1, 3) + (273.0L/16.0L)*std::pow(2*xi_2 - 1, 2) + 791.0L/64.0L ;
}
if( b == 47 ) {
  valDxi1 = (1.0L/64.0L)*(24*xi_1 + 12*xi_2 - 12)*(910*xi_2 + 12597*std::pow(2*xi_2 - 1, 7) + 23205*std::pow(2*xi_2 - 1, 6) + 4641*std::pow(2*xi_2 - 1, 5) - 11375*std::pow(2*xi_2 - 1, 4) - 5005*std::pow(2*xi_2 - 1, 3) + 819*std::pow(2*xi_2 - 1, 2) - 448) ;
}
if( b == 48 ) {
  valDxi1 = -(6*std::pow(xi_2 - 1, 3) - 30*(xi_2 - 1)*std::pow(2*xi_1 + xi_2 - 1, 2))*(-588*xi_2 + 6783*std::pow(2*xi_2 - 1, 6) + 14994*std::pow(2*xi_2 - 1, 5) + 8925*std::pow(2*xi_2 - 1, 4) - 980*std::pow(2*xi_2 - 1, 3) - 1995*std::pow(2*xi_2 - 1, 2) + 317)/(32*xi_2 - 32) ;
}
if( b == 49 ) {
  valDxi1 = (1.0L/64.0L)*(-30*std::pow(xi_2 - 1, 2)*(8*xi_1 + 4*xi_2 - 4) + 280*std::pow(2*xi_1 + xi_2 - 1, 3))*(-330*xi_2 + 2907*std::pow(2*xi_2 - 1, 5) + 6885*std::pow(2*xi_2 - 1, 4) + 5270*std::pow(2*xi_2 - 1, 3) + 1170*std::pow(2*xi_2 - 1, 2) + 114) ;
}
if( b == 50 ) {
  valDxi1 = (1.0L/32.0L)*(-630*(-xi_2 + 1)*std::pow(-2*xi_1 - xi_2 + 1, 4) + 30*std::pow(xi_2 - 1, 5) - 420*std::pow(xi_2 - 1, 3)*std::pow(-2*xi_1 - xi_2 + 1, 2))*(968*xi_2 + 969*std::pow(2*xi_2 - 1, 4) + 2244*std::pow(2*xi_2 - 1, 3) + 1734*std::pow(2*xi_2 - 1, 2) - 455)/(xi_2 - 1) ;
}
if( b == 51 ) {
  valDxi1 = (1.0L/128.0L)*(105*std::pow(xi_2 - 1, 4)*(8*xi_1 + 4*xi_2 - 4) - 2520*std::pow(xi_2 - 1, 2)*std::pow(2*xi_1 + xi_2 - 1, 3) + 2772*std::pow(2*xi_1 + xi_2 - 1, 5))*(2550*xi_2 + 969*std::pow(2*xi_2 - 1, 3) + 1989*std::pow(2*xi_2 - 1, 2) - 1028) ;
}
if( b == 52 ) {
  valDxi1 = -(540*xi_2 + 171*std::pow(2*xi_2 - 1, 2) - 167)*(70*std::pow(xi_2 - 1, 7) - 1890*std::pow(xi_2 - 1, 5)*std::pow(2*xi_1 + xi_2 - 1, 2) + 6930*std::pow(xi_2 - 1, 3)*std::pow(2*xi_1 + xi_2 - 1, 4) - 6006*(xi_2 - 1)*std::pow(2*xi_1 + xi_2 - 1, 6))/(64*xi_2 - 64) ;
}
if( b == 53 ) {
  valDxi1 = (1.0L/128.0L)*(19*xi_2 - 1)*(-1260*std::pow(xi_2 - 1, 6)*(8*xi_1 + 4*xi_2 - 4) + 55440*std::pow(xi_2 - 1, 4)*std::pow(2*xi_1 + xi_2 - 1, 3) - 144144*std::pow(xi_2 - 1, 2)*std::pow(2*xi_1 + xi_2 - 1, 5) + 102960*std::pow(2*xi_1 + xi_2 - 1, 7)) ;
}
if( b == 54 ) {
  valDxi1 = (1.0L/128.0L)*(-180180*std::pow(-xi_2 + 1, 5)*std::pow(-2*xi_1 - xi_2 + 1, 4) - 218790*(-xi_2 + 1)*std::pow(-2*xi_1 - xi_2 + 1, 8) + 630*std::pow(xi_2 - 1, 9) - 27720*std::pow(xi_2 - 1, 7)*std::pow(-2*xi_1 - xi_2 + 1, 2) - 360360*std::pow(xi_2 - 1, 3)*std::pow(-2*xi_1 - xi_2 + 1, 6))/(xi_2 - 1) ;
}
