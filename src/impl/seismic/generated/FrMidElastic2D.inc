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
 * Generated right side's contribution to the Riemann solution in the face-aligned coordinate system.
 **/

o_frMid[0][0]=-((std::pow((i_lamL + 2*i_muL)/i_rhoL,1.5)*i_rhoL*std::sqrt((i_lamR + 2*i_muR)/i_rhoR))/(i_lamR*std::sqrt((i_lamL + 2*i_muL)/i_rhoL) + 2*i_muR*std::sqrt((i_lamL + 2*i_muL)/i_rhoL) + (i_lamL + 2*i_muL)*std::sqrt((i_lamR + 2*i_muR)/i_rhoR)));
o_frMid[0][1]=0;
o_frMid[0][2]=0;
o_frMid[0][3]=-(((i_lamR + 2*i_muR)*std::pow((i_lamL + 2*i_muL)/i_rhoL,1.5)*i_rhoL)/(i_lamR*std::sqrt((i_lamL + 2*i_muL)/i_rhoL) + 2*i_muR*std::sqrt((i_lamL + 2*i_muL)/i_rhoL) + (i_lamL + 2*i_muL)*std::sqrt((i_lamR + 2*i_muR)/i_rhoR)));
o_frMid[0][4]=0;
o_frMid[1][0]=-((i_lamL*std::sqrt((i_lamL + 2*i_muL)/i_rhoL)*std::sqrt((i_lamR + 2*i_muR)/i_rhoR))/(i_lamR*std::sqrt((i_lamL + 2*i_muL)/i_rhoL) + 2*i_muR*std::sqrt((i_lamL + 2*i_muL)/i_rhoL) + (i_lamL + 2*i_muL)*std::sqrt((i_lamR + 2*i_muR)/i_rhoR)));
o_frMid[1][1]=0;
o_frMid[1][2]=0;
o_frMid[1][3]=-((i_lamL*(i_lamR + 2*i_muR)*std::sqrt((i_lamL + 2*i_muL)/i_rhoL))/(i_lamR*std::sqrt((i_lamL + 2*i_muL)/i_rhoL) + 2*i_muR*std::sqrt((i_lamL + 2*i_muL)/i_rhoL) + (i_lamL + 2*i_muL)*std::sqrt((i_lamR + 2*i_muR)/i_rhoR)));
o_frMid[1][4]=0;
o_frMid[2][0]=0;
o_frMid[2][1]=0;
o_frMid[2][2]=-((i_muL*std::sqrt(i_muL/i_rhoL)*std::sqrt(i_muR/i_rhoR))/(i_muR*std::sqrt(i_muL/i_rhoL) + i_muL*std::sqrt(i_muR/i_rhoR)));
o_frMid[2][3]=0;
o_frMid[2][4]=-((i_muL*i_muR*std::sqrt(i_muL/i_rhoL))/(i_muR*std::sqrt(i_muL/i_rhoL) + i_muL*std::sqrt(i_muR/i_rhoR)));
o_frMid[3][0]=-(((i_lamL + 2*i_muL)*std::sqrt((i_lamR + 2*i_muR)/i_rhoR))/(i_rhoL*(i_lamR*std::sqrt((i_lamL + 2*i_muL)/i_rhoL) + 2*i_muR*std::sqrt((i_lamL + 2*i_muL)/i_rhoL) + (i_lamL + 2*i_muL)*std::sqrt((i_lamR + 2*i_muR)/i_rhoR))));
o_frMid[3][1]=0;
o_frMid[3][2]=0;
o_frMid[3][3]=-(((i_lamL + 2*i_muL)*(i_lamR + 2*i_muR))/(i_rhoL*(i_lamR*std::sqrt((i_lamL + 2*i_muL)/i_rhoL) + 2*i_muR*std::sqrt((i_lamL + 2*i_muL)/i_rhoL) + (i_lamL + 2*i_muL)*std::sqrt((i_lamR + 2*i_muR)/i_rhoR))));
o_frMid[3][4]=0;
o_frMid[4][0]=0;
o_frMid[4][1]=0;
o_frMid[4][2]=-((i_muL*std::sqrt(i_muR/i_rhoR))/(i_muR*std::sqrt(i_muL/i_rhoL)*i_rhoL + i_muL*i_rhoL*std::sqrt(i_muR/i_rhoR)));
o_frMid[4][3]=0;
o_frMid[4][4]=-((i_muL*i_muR)/(i_muR*std::sqrt(i_muL/i_rhoL)*i_rhoL + i_muL*i_rhoL*std::sqrt(i_muR/i_rhoR)));
