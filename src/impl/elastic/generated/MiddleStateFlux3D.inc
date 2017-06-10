/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2017, Regents of the University of California
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
 * Generated 3D solver flux computation based on derived middle states.
 **/
o_msFlux[0][0]=0;
o_msFlux[0][1]=0;
o_msFlux[0][2]=0;
o_msFlux[0][3]=0;
o_msFlux[0][4]=0;
o_msFlux[0][5]=0;
o_msFlux[0][6]=-2*i_mu*std::pow(i_nx,2) - i_lam*(std::pow(i_nx,2) + std::pow(i_sx,2) + std::pow(i_tx,2));
o_msFlux[0][7]=-2*i_mu*i_nx*i_sx;
o_msFlux[0][8]=-2*i_mu*i_nx*i_tx;
o_msFlux[1][0]=0;
o_msFlux[1][1]=0;
o_msFlux[1][2]=0;
o_msFlux[1][3]=0;
o_msFlux[1][4]=0;
o_msFlux[1][5]=0;
o_msFlux[1][6]=-2*i_mu*std::pow(i_ny,2) - i_lam*(std::pow(i_ny,2) + std::pow(i_sy,2) + std::pow(i_ty,2));
o_msFlux[1][7]=-2*i_mu*i_ny*i_sy;
o_msFlux[1][8]=-2*i_mu*i_ny*i_ty;
o_msFlux[2][0]=0;
o_msFlux[2][1]=0;
o_msFlux[2][2]=0;
o_msFlux[2][3]=0;
o_msFlux[2][4]=0;
o_msFlux[2][5]=0;
o_msFlux[2][6]=-2*i_mu*std::pow(i_nz,2) - i_lam*(std::pow(i_nz,2) + std::pow(i_sz,2) + std::pow(i_tz,2));
o_msFlux[2][7]=-2*i_mu*i_nz*i_sz;
o_msFlux[2][8]=-2*i_mu*i_nz*i_tz;
o_msFlux[3][0]=0;
o_msFlux[3][1]=0;
o_msFlux[3][2]=0;
o_msFlux[3][3]=0;
o_msFlux[3][4]=0;
o_msFlux[3][5]=0;
o_msFlux[3][6]=-2*i_mu*i_nx*i_ny - i_lam*(i_nx*i_ny + i_sx*i_sy + i_tx*i_ty);
o_msFlux[3][7]=-(i_mu*(i_ny*i_sx + i_nx*i_sy));
o_msFlux[3][8]=-(i_mu*(i_ny*i_tx + i_nx*i_ty));
o_msFlux[4][0]=0;
o_msFlux[4][1]=0;
o_msFlux[4][2]=0;
o_msFlux[4][3]=0;
o_msFlux[4][4]=0;
o_msFlux[4][5]=0;
o_msFlux[4][6]=-2*i_mu*i_ny*i_nz - i_lam*(i_ny*i_nz + i_sy*i_sz + i_ty*i_tz);
o_msFlux[4][7]=-(i_mu*(i_nz*i_sy + i_ny*i_sz));
o_msFlux[4][8]=-(i_mu*(i_nz*i_ty + i_ny*i_tz));
o_msFlux[5][0]=0;
o_msFlux[5][1]=0;
o_msFlux[5][2]=0;
o_msFlux[5][3]=0;
o_msFlux[5][4]=0;
o_msFlux[5][5]=0;
o_msFlux[5][6]=-2*i_mu*i_nx*i_nz - i_lam*(i_nx*i_nz + i_sx*i_sz + i_tx*i_tz);
o_msFlux[5][7]=-(i_mu*(i_nz*i_sx + i_nx*i_sz));
o_msFlux[5][8]=-(i_mu*(i_nz*i_tx + i_nx*i_tz));
o_msFlux[6][0]=-(i_nx/i_rho);
o_msFlux[6][1]=0;
o_msFlux[6][2]=0;
o_msFlux[6][3]=-(i_sx/i_rho);
o_msFlux[6][4]=0;
o_msFlux[6][5]=-(i_tx/i_rho);
o_msFlux[6][6]=0;
o_msFlux[6][7]=0;
o_msFlux[6][8]=0;
o_msFlux[7][0]=-(i_ny/i_rho);
o_msFlux[7][1]=0;
o_msFlux[7][2]=0;
o_msFlux[7][3]=-(i_sy/i_rho);
o_msFlux[7][4]=0;
o_msFlux[7][5]=-(i_ty/i_rho);
o_msFlux[7][6]=0;
o_msFlux[7][7]=0;
o_msFlux[7][8]=0;
o_msFlux[8][0]=-(i_nz/i_rho);
o_msFlux[8][1]=0;
o_msFlux[8][2]=0;
o_msFlux[8][3]=-(i_sz/i_rho);
o_msFlux[8][4]=0;
o_msFlux[8][5]=-(i_tz/i_rho);
o_msFlux[8][6]=0;
o_msFlux[8][7]=0;
o_msFlux[8][8]=0;
