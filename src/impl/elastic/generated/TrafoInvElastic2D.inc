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
 * Generated inverse transformation of quantities, from physical coordinates to a local, face-aligned coordinate system in 2D.
 **/

o_tm1[0][0]=1 - std::pow(i_ny,2);
o_tm1[0][1]=std::pow(i_ny,2);
o_tm1[0][2]=2*i_nx*i_ny;
o_tm1[0][3]=0;
o_tm1[0][4]=0;
o_tm1[1][0]=std::pow(i_ny,2);
o_tm1[1][1]=1 - std::pow(i_ny,2);
o_tm1[1][2]=-2*i_nx*i_ny;
o_tm1[1][3]=0;
o_tm1[1][4]=0;
o_tm1[2][0]=-(i_nx*i_ny);
o_tm1[2][1]=i_nx*i_ny;
o_tm1[2][2]=1 - 2*std::pow(i_ny,2);
o_tm1[2][3]=0;
o_tm1[2][4]=0;
o_tm1[3][0]=0;
o_tm1[3][1]=0;
o_tm1[3][2]=0;
o_tm1[3][3]=i_nx;
o_tm1[3][4]=i_ny;
o_tm1[4][0]=0;
o_tm1[4][1]=0;
o_tm1[4][2]=0;
o_tm1[4][3]=-i_ny;
o_tm1[4][4]=i_nx;
