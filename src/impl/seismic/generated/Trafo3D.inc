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
 * Generated transformation of quantities, from a local, face-aligned coordinate system to physical xyz-coords in 3D.
 **/

o_t[0][0]=std::pow(i_nx,2);
o_t[0][1]=std::pow(i_sx,2);
o_t[0][2]=std::pow(i_tx,2);
o_t[0][3]=2*i_nx*i_sx;
o_t[0][4]=2*i_sx*i_tx;
o_t[0][5]=2*i_nx*i_tx;
o_t[0][6]=0;
o_t[0][7]=0;
o_t[0][8]=0;
o_t[1][0]=std::pow(i_ny,2);
o_t[1][1]=std::pow(i_sy,2);
o_t[1][2]=std::pow(i_ty,2);
o_t[1][3]=2*i_ny*i_sy;
o_t[1][4]=2*i_sy*i_ty;
o_t[1][5]=2*i_ny*i_ty;
o_t[1][6]=0;
o_t[1][7]=0;
o_t[1][8]=0;
o_t[2][0]=std::pow(i_nz,2);
o_t[2][1]=std::pow(i_sz,2);
o_t[2][2]=std::pow(i_tz,2);
o_t[2][3]=2*i_nz*i_sz;
o_t[2][4]=2*i_sz*i_tz;
o_t[2][5]=2*i_nz*i_tz;
o_t[2][6]=0;
o_t[2][7]=0;
o_t[2][8]=0;
o_t[3][0]=i_nx*i_ny;
o_t[3][1]=i_sx*i_sy;
o_t[3][2]=i_tx*i_ty;
o_t[3][3]=i_ny*i_sx + i_nx*i_sy;
o_t[3][4]=i_sy*i_tx + i_sx*i_ty;
o_t[3][5]=i_ny*i_tx + i_nx*i_ty;
o_t[3][6]=0;
o_t[3][7]=0;
o_t[3][8]=0;
o_t[4][0]=i_ny*i_nz;
o_t[4][1]=i_sy*i_sz;
o_t[4][2]=i_ty*i_tz;
o_t[4][3]=i_nz*i_sy + i_ny*i_sz;
o_t[4][4]=i_sz*i_ty + i_sy*i_tz;
o_t[4][5]=i_nz*i_ty + i_ny*i_tz;
o_t[4][6]=0;
o_t[4][7]=0;
o_t[4][8]=0;
o_t[5][0]=i_nx*i_nz;
o_t[5][1]=i_sx*i_sz;
o_t[5][2]=i_tx*i_tz;
o_t[5][3]=i_nz*i_sx + i_nx*i_sz;
o_t[5][4]=i_sz*i_tx + i_sx*i_tz;
o_t[5][5]=i_nz*i_tx + i_nx*i_tz;
o_t[5][6]=0;
o_t[5][7]=0;
o_t[5][8]=0;
o_t[6][0]=0;
o_t[6][1]=0;
o_t[6][2]=0;
o_t[6][3]=0;
o_t[6][4]=0;
o_t[6][5]=0;
o_t[6][6]=i_nx;
o_t[6][7]=i_sx;
o_t[6][8]=i_tx;
o_t[7][0]=0;
o_t[7][1]=0;
o_t[7][2]=0;
o_t[7][3]=0;
o_t[7][4]=0;
o_t[7][5]=0;
o_t[7][6]=i_ny;
o_t[7][7]=i_sy;
o_t[7][8]=i_ty;
o_t[8][0]=0;
o_t[8][1]=0;
o_t[8][2]=0;
o_t[8][3]=0;
o_t[8][4]=0;
o_t[8][5]=0;
o_t[8][6]=i_nz;
o_t[8][7]=i_sz;
o_t[8][8]=i_tz;
