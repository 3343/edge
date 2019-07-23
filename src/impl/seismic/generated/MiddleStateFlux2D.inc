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
 * Generated 2D solver flux computation based on derived middle states.
 **/
o_msFlux[0][0]=0;
o_msFlux[0][1]=0;
o_msFlux[0][2]=0;
o_msFlux[0][3]=-i_lam + 2*i_mu*(-1 + std::pow(i_ny,2));
o_msFlux[0][4]=2*i_mu*i_nx*i_ny;
o_msFlux[1][0]=0;
o_msFlux[1][1]=0;
o_msFlux[1][2]=0;
o_msFlux[1][3]=-i_lam - 2*i_mu*std::pow(i_ny,2);
o_msFlux[1][4]=-2*i_mu*i_nx*i_ny;
o_msFlux[2][0]=0;
o_msFlux[2][1]=0;
o_msFlux[2][2]=0;
o_msFlux[2][3]=-2*i_mu*i_nx*i_ny;
o_msFlux[2][4]=i_mu*(-1 + 2*std::pow(i_ny,2));
o_msFlux[3][0]=-(i_nx/i_rho);
o_msFlux[3][1]=0;
o_msFlux[3][2]=i_ny/i_rho;
o_msFlux[3][3]=0;
o_msFlux[3][4]=0;
o_msFlux[4][0]=-(i_ny/i_rho);
o_msFlux[4][1]=0;
o_msFlux[4][2]=-(i_nx/i_rho);
o_msFlux[4][3]=0;
o_msFlux[4][4]=0;
