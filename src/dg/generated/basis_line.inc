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
 * Generated 1D basis.
 **/

assert( b < 10 );
if( b == 0 ) {
  val    = 1 ;
}
if( b == 1 ) {
  val    = 2*xi - 1 ;
}
if( b == 2 ) {
  val    = (3.0L/2.0L)*std::pow(2*xi - 1, 2) - 1.0L/2.0L ;
}
if( b == 3 ) {
  val    = -3*xi + (5.0L/2.0L)*std::pow(2*xi - 1, 3) + 3.0L/2.0L ;
}
if( b == 4 ) {
  val    = (35.0L/8.0L)*std::pow(2*xi - 1, 4) - 15.0L/4.0L*std::pow(2*xi - 1, 2) + 3.0L/8.0L ;
}
if( b == 5 ) {
  val    = (15.0L/4.0L)*xi + (63.0L/8.0L)*std::pow(2*xi - 1, 5) - 35.0L/4.0L*std::pow(2*xi - 1, 3) - 15.0L/8.0L ;
}
if( b == 6 ) {
  val    = (231.0L/16.0L)*std::pow(2*xi - 1, 6) - 315.0L/16.0L*std::pow(2*xi - 1, 4) + (105.0L/16.0L)*std::pow(2*xi - 1, 2) - 5.0L/16.0L ;
}
if( b == 7 ) {
  val    = -35.0L/8.0L*xi + (429.0L/16.0L)*std::pow(2*xi - 1, 7) - 693.0L/16.0L*std::pow(2*xi - 1, 5) + (315.0L/16.0L)*std::pow(2*xi - 1, 3) + 35.0L/16.0L ;
}
if( b == 8 ) {
  val    = (6435.0L/128.0L)*std::pow(2*xi - 1, 8) - 3003.0L/32.0L*std::pow(2*xi - 1, 6) + (3465.0L/64.0L)*std::pow(2*xi - 1, 4) - 315.0L/32.0L*std::pow(2*xi - 1, 2) + 35.0L/128.0L ;
}
if( b == 9 ) {
  val    = (315.0L/64.0L)*xi + (12155.0L/128.0L)*std::pow(2*xi - 1, 9) - 6435.0L/32.0L*std::pow(2*xi - 1, 7) + (9009.0L/64.0L)*std::pow(2*xi - 1, 5) - 1155.0L/32.0L*std::pow(2*xi - 1, 3) - 315.0L/128.0L ;
}
