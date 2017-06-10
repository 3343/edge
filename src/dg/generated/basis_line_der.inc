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
 * Generated derivatives of 1D basis.
 **/

assert( b < 10 );
if( b == 0 ) {
  valDxi = 0 ;
}
if( b == 1 ) {
  valDxi = 2 ;
}
if( b == 2 ) {
  valDxi = 12*xi - 6 ;
}
if( b == 3 ) {
  valDxi = 15*std::pow(2*xi - 1, 2) - 3 ;
}
if( b == 4 ) {
  valDxi = -30*xi + 35*std::pow(2*xi - 1, 3) + 15 ;
}
if( b == 5 ) {
  valDxi = (315.0L/4.0L)*std::pow(2*xi - 1, 4) - 105.0L/2.0L*std::pow(2*xi - 1, 2) + 15.0L/4.0L ;
}
if( b == 6 ) {
  valDxi = (105.0L/2.0L)*xi + (693.0L/4.0L)*std::pow(2*xi - 1, 5) - 315.0L/2.0L*std::pow(2*xi - 1, 3) - 105.0L/4.0L ;
}
if( b == 7 ) {
  valDxi = (3003.0L/8.0L)*std::pow(2*xi - 1, 6) - 3465.0L/8.0L*std::pow(2*xi - 1, 4) + (945.0L/8.0L)*std::pow(2*xi - 1, 2) - 35.0L/8.0L ;
}
if( b == 8 ) {
  valDxi = -315.0L/4.0L*xi + (6435.0L/8.0L)*std::pow(2*xi - 1, 7) - 9009.0L/8.0L*std::pow(2*xi - 1, 5) + (3465.0L/8.0L)*std::pow(2*xi - 1, 3) + 315.0L/8.0L ;
}
if( b == 9 ) {
  valDxi = (109395.0L/64.0L)*std::pow(2*xi - 1, 8) - 45045.0L/16.0L*std::pow(2*xi - 1, 6) + (45045.0L/32.0L)*std::pow(2*xi - 1, 4) - 3465.0L/16.0L*std::pow(2*xi - 1, 2) + 315.0L/64.0L ;
}
