/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2018, Regents of the University of California
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
 * Unit test for super-cells.
 **/
#include <catch.hpp>
#define private public
#include "SuperCell.hpp"
#undef private

unsigned short l_scSfScTet1[63][4] = {
  { 2,4,1,3 },
  { 6,7,14,0 },
  { 0,8,7,15 },
  { 25,0,21,23 },
  { 0,26,22,23 },
  { 27,36,45,20 },
  { 28,1,20,10 },
  { 32,38,1,2 },
  { 33,24,2,12 },
  { 35,40,24,54 },
  { 29,6,50,21 },
  { 30,25,21,13 },
  { 34,8,25,56 },
  { 31,11,53,58 },
  { 1,37,20,16 },
  { 24,39,2,18 },
  { 14,41,47,22 },
  { 26,42,22,19 },
  { 15,43,26,59 },
  { 17,44,49,62 },
  { 46,6,5,14 },
  { 51,11,10,3 },
  { 48,4,16,17 },
  { 3,4,52,61 },
  { 8,9,15,55 },
  { 11,12,3,57 },
  { 4,18,17,60 },
  { 5,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 6,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 10,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 11,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 13,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 7,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 8,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 12,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 9,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 5,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 14,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 7,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 15,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 9,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 16,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 17,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 18,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 19,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 5,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 20,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 16,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 22,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 19,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 10,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 21,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 23,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 13,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 9,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 24,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 12,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 25,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 13,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 18,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 26,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 23,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 19,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() }
};

unsigned short l_scSfScTet2[225][4] = {
  { 3,25,1,10 },
  { 52,53,76,0 },
  { 6,28,4,13 },
  { 0,54,53,77 },
  { 54,55,77,2 },
  { 8,31,7,16 },
  { 2,56,55,78 },
  { 56,57,78,5 },
  { 5,58,57,79 },
  { 14,34,11,19 },
  { 11,0,97,101 },
  { 61,62,10,9 },
  { 17,38,15,21 },
  { 15,2,14,37 },
  { 9,63,62,13 },
  { 63,64,13,12 },
  { 113,5,17,118 },
  { 12,65,64,16 },
  { 22,41,20,23 },
  { 20,9,98,103 },
  { 68,69,19,18 },
  { 114,12,22,120 },
  { 18,70,69,21 },
  { 115,18,99,105 },
  { 29,44,26,35 },
  { 0,26,100,101 },
  { 25,82,81,24 },
  { 32,46,30,39 },
  { 2,30,29,37 },
  { 24,28,82,83 },
  { 28,84,83,27 },
  { 5,116,32,118 },
  { 27,31,84,85 },
  { 40,48,36,42 },
  { 9,36,102,103 },
  { 36,24,102,107 },
  { 34,37,35,33 },
  { 13,28,36,40 },
  { 12,117,40,120 },
  { 117,27,40,123 },
  { 33,38,37,39 },
  { 18,119,104,105 },
  { 119,33,104,109 },
  { 47,50,45,49 },
  { 24,45,106,107 },
  { 44,89,88,43 },
  { 27,121,47,123 },
  { 43,46,89,90 },
  { 33,122,108,109 },
  { 122,43,108,111 },
  { 43,124,110,111 },
  { 125,150,175,96 },
  { 126,1,96,60 },
  { 134,152,1,3 },
  { 135,4,3,62 },
  { 141,154,4,6 },
  { 142,7,6,64 },
  { 146,156,7,8 },
  { 147,112,8,66 },
  { 149,158,112,200 },
  { 127,52,184,97 },
  { 128,11,97,67 },
  { 136,54,11,14 },
  { 137,15,14,69 },
  { 143,56,15,17 },
  { 144,113,17,71 },
  { 148,58,113,202 },
  { 129,61,191,98 },
  { 130,20,98,72 },
  { 138,63,20,22 },
  { 139,114,22,74 },
  { 145,65,114,204 },
  { 131,68,196,99 },
  { 132,115,99,75 },
  { 140,70,115,206 },
  { 133,73,199,208 },
  { 1,151,96,80 },
  { 4,153,3,82 },
  { 7,155,6,84 },
  { 112,157,8,86 },
  { 76,159,177,100 },
  { 26,160,100,87 },
  { 77,161,26,29 },
  { 30,162,29,89 },
  { 78,163,30,32 },
  { 116,164,32,91 },
  { 79,165,116,209 },
  { 81,166,179,106 },
  { 45,167,106,92 },
  { 83,168,45,47 },
  { 121,169,47,94 },
  { 85,170,121,216 },
  { 88,171,181,110 },
  { 124,172,110,95 },
  { 90,173,124,221 },
  { 93,174,183,224 },
  { 176,52,51,76 },
  { 185,61,60,10 },
  { 192,68,67,19 },
  { 197,73,72,23 },
  { 178,25,80,81 },
  { 10,25,186,102 },
  { 187,34,101,35 },
  { 19,34,193,104 },
  { 194,41,103,42 },
  { 23,41,198,215 },
  { 180,44,87,88 },
  { 35,44,188,108 },
  { 189,48,107,49 },
  { 42,48,195,220 },
  { 182,50,92,93 },
  { 49,50,190,223 },
  { 58,59,79,201 },
  { 65,66,16,203 },
  { 70,71,21,205 },
  { 73,74,23,207 },
  { 31,86,85,210 },
  { 38,118,39,212 },
  { 16,31,117,211 },
  { 41,120,42,214 },
  { 21,38,119,213 },
  { 46,91,90,217 },
  { 48,123,49,219 },
  { 39,46,122,218 },
  { 50,94,93,222 },
  { 51,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 52,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 60,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 61,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 67,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 68,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 72,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 73,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 75,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 53,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 54,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 62,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 63,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 69,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 70,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 74,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 55,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 56,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 64,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 65,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 71,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 57,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 58,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 66,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 59,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 51,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 76,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 53,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 77,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 55,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 78,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 57,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 79,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 59,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 80,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 81,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 82,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 83,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 84,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 85,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 86,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 87,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 88,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 89,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 90,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 91,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 92,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 93,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 94,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 95,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 51,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 96,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 80,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 100,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 87,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 106,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 92,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 110,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 95,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 60,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 97,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 101,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 102,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 107,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 108,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 111,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 67,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 98,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 103,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 104,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 109,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 72,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 99,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 105,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 75,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 59,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 112,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 66,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 113,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 71,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 114,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 74,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 115,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 75,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 86,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 116,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 118,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 117,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 120,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 119,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 105,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 91,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 121,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 123,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 122,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 109,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 94,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 124,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 111,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() },
  { 95,std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max(),std::numeric_limits< unsigned short >::max() } };

unsigned short const l_scDgAdTet1[3][9] = {
    {  0,1,5,6,8,2,3,7,4,  },
    {  4,3,2,1,0,7,6,5,8,  },
    {  8,6,7,3,4,5,1,2,0,  },
};

unsigned short const l_scDgAdTet2[3][25] = {
  { 0,1,9,10,16,17,21,22,24,2,3,11,12,18,19,23,4,5,13,14,20,6,7,15,8 },
  { 8,7,6,5,4,3,2,1,0,15,14,13,12,11,10,9,20,19,18,17,16,23,22,21,24 },
  { 24,22,23,19,20,14,15,7,8,21,17,18,12,13,5,6,16,10,11,3,4,9,1,2,0 },
};

TEST_CASE( "SuperCells: Default stencils.", "[superCell][tet4]" ) {
  edge::sc::ibnd::t_SuperCell< TET4,
                               2 > l_su1;

  edge::sc::ibnd::SuperCellInit< TET4,
                                 2 >::faSt( l_scDgAdTet1,
                                            l_scSfScTet1,
                                            l_su1 );

  // check bottom-face
  REQUIRE( l_su1.col[0].sc[0][0] ==  6 );
  REQUIRE( l_su1.col[0].sc[0][1] == 10 );

  REQUIRE( l_su1.col[0].sc[1][0] == 11 );
  REQUIRE( l_su1.col[0].sc[1][1] == 13 );

  REQUIRE( l_su1.col[0].sc[2][0] ==  8 );
  REQUIRE( l_su1.col[0].sc[2][1] == 12 );

  // check front-face
  REQUIRE( l_su1.col[1].sc[0][0] == 14 );
  REQUIRE( l_su1.col[1].sc[0][1] == 16 );

  REQUIRE( l_su1.col[1].sc[1][0] == 15 );
  REQUIRE( l_su1.col[1].sc[1][1] == 18 );

  REQUIRE( l_su1.col[1].sc[2][0] == 17 );
  REQUIRE( l_su1.col[1].sc[2][1] == 19 );

  edge::sc::ibnd::t_SuperCell< TET4,
                               3 > l_su2;

  edge::sc::ibnd::SuperCellInit< TET4,
                                 3 >::faSt( l_scDgAdTet2,
                                            l_scSfScTet2,
                                            l_su2 );

  // check diagonal face (third order)
  REQUIRE( l_su2.col[3].sc[0][0] ==  59 );
  REQUIRE( l_su2.col[3].sc[0][1] == 112 );

  REQUIRE( l_su2.col[3].sc[1][0] ==  66 );
  REQUIRE( l_su2.col[3].sc[1][1] == 113 );

  REQUIRE( l_su2.col[3].sc[2][0] ==  71 );
  REQUIRE( l_su2.col[3].sc[2][1] == 114 );

  REQUIRE( l_su2.col[3].sc[3][0] ==  74 );
  REQUIRE( l_su2.col[3].sc[3][1] == 115 );

  REQUIRE( l_su2.col[3].sc[4][0] ==  86 );
  REQUIRE( l_su2.col[3].sc[4][1] == 116 );

  REQUIRE( l_su2.col[3].sc[5][0] == 117 );
  REQUIRE( l_su2.col[3].sc[5][1] == 118 );

  REQUIRE( l_su2.col[3].sc[6][0] == 119 );
  REQUIRE( l_su2.col[3].sc[6][1] == 120 );

  REQUIRE( l_su2.col[3].sc[7][0] ==  91 );
  REQUIRE( l_su2.col[3].sc[7][1] == 121 );

  REQUIRE( l_su2.col[3].sc[8][0] == 122 );
  REQUIRE( l_su2.col[3].sc[8][1] == 123 );

  REQUIRE( l_su2.col[3].sc[9][0] ==  94 );
  REQUIRE( l_su2.col[3].sc[9][1] == 124 );
}