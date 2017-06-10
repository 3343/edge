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
 * Generated 3D solver for jump to middle state from the right side.
 **/
o_msJump[0][0]=-(((i_lamR + 2*i_muR)*std::sqrt((i_lamL + 2*i_muL)/i_rhoL))/((i_lamR + 2*i_muR)*std::sqrt((i_lamL + 2*i_muL)/i_rhoL) + (i_lamL + 2*i_muL)*std::sqrt((i_lamR + 2*i_muR)/i_rhoR)));
o_msJump[0][1]=0;
o_msJump[0][2]=0;
o_msJump[0][3]=0;
o_msJump[0][4]=0;
o_msJump[0][5]=0;
o_msJump[0][6]=((i_lamL + 2*i_muL)*(i_lamR + 2*i_muR))/((i_lamR + 2*i_muR)*std::sqrt((i_lamL + 2*i_muL)/i_rhoL) + (i_lamL + 2*i_muL)*std::sqrt((i_lamR + 2*i_muR)/i_rhoR));
o_msJump[0][7]=0;
o_msJump[0][8]=0;
o_msJump[1][0]=-((i_lamR*std::sqrt((i_lamL + 2*i_muL)/i_rhoL))/((i_lamR + 2*i_muR)*std::sqrt((i_lamL + 2*i_muL)/i_rhoL) + (i_lamL + 2*i_muL)*std::sqrt((i_lamR + 2*i_muR)/i_rhoR)));
o_msJump[1][1]=0;
o_msJump[1][2]=0;
o_msJump[1][3]=0;
o_msJump[1][4]=0;
o_msJump[1][5]=0;
o_msJump[1][6]=(i_lamR*(i_lamL + 2*i_muL))/((i_lamR + 2*i_muR)*std::sqrt((i_lamL + 2*i_muL)/i_rhoL) + (i_lamL + 2*i_muL)*std::sqrt((i_lamR + 2*i_muR)/i_rhoR));
o_msJump[1][7]=0;
o_msJump[1][8]=0;
o_msJump[2][0]=-((i_lamR*std::sqrt((i_lamL + 2*i_muL)/i_rhoL))/((i_lamR + 2*i_muR)*std::sqrt((i_lamL + 2*i_muL)/i_rhoL) + (i_lamL + 2*i_muL)*std::sqrt((i_lamR + 2*i_muR)/i_rhoR)));
o_msJump[2][1]=0;
o_msJump[2][2]=0;
o_msJump[2][3]=0;
o_msJump[2][4]=0;
o_msJump[2][5]=0;
o_msJump[2][6]=(i_lamR*(i_lamL + 2*i_muL))/((i_lamR + 2*i_muR)*std::sqrt((i_lamL + 2*i_muL)/i_rhoL) + (i_lamL + 2*i_muL)*std::sqrt((i_lamR + 2*i_muR)/i_rhoR));
o_msJump[2][7]=0;
o_msJump[2][8]=0;
o_msJump[3][0]=0;
o_msJump[3][1]=0;
o_msJump[3][2]=0;
o_msJump[3][3]=-((i_muR*std::sqrt(i_muL/i_rhoL))/(i_muR*std::sqrt(i_muL/i_rhoL) + i_muL*std::sqrt(i_muR/i_rhoR)));
o_msJump[3][4]=0;
o_msJump[3][5]=0;
o_msJump[3][6]=0;
o_msJump[3][7]=(i_muL*i_muR)/(i_muR*std::sqrt(i_muL/i_rhoL) + i_muL*std::sqrt(i_muR/i_rhoR));
o_msJump[3][8]=0;
o_msJump[4][0]=0;
o_msJump[4][1]=0;
o_msJump[4][2]=0;
o_msJump[4][3]=0;
o_msJump[4][4]=0;
o_msJump[4][5]=0;
o_msJump[4][6]=0;
o_msJump[4][7]=0;
o_msJump[4][8]=0;
o_msJump[5][0]=0;
o_msJump[5][1]=0;
o_msJump[5][2]=0;
o_msJump[5][3]=0;
o_msJump[5][4]=0;
o_msJump[5][5]=-((i_muR*std::sqrt(i_muL/i_rhoL))/(i_muR*std::sqrt(i_muL/i_rhoL) + i_muL*std::sqrt(i_muR/i_rhoR)));
o_msJump[5][6]=0;
o_msJump[5][7]=0;
o_msJump[5][8]=(i_muL*i_muR)/(i_muR*std::sqrt(i_muL/i_rhoL) + i_muL*std::sqrt(i_muR/i_rhoR));
o_msJump[6][0]=(std::sqrt((i_lamL + 2*i_muL)/i_rhoL)*std::sqrt((i_lamR + 2*i_muR)/i_rhoR))/((i_lamR + 2*i_muR)*std::sqrt((i_lamL + 2*i_muL)/i_rhoL) + (i_lamL + 2*i_muL)*std::sqrt((i_lamR + 2*i_muR)/i_rhoR));
o_msJump[6][1]=0;
o_msJump[6][2]=0;
o_msJump[6][3]=0;
o_msJump[6][4]=0;
o_msJump[6][5]=0;
o_msJump[6][6]=-(((i_lamL + 2*i_muL)*std::sqrt((i_lamR + 2*i_muR)/i_rhoR))/((i_lamR + 2*i_muR)*std::sqrt((i_lamL + 2*i_muL)/i_rhoL) + (i_lamL + 2*i_muL)*std::sqrt((i_lamR + 2*i_muR)/i_rhoR)));
o_msJump[6][7]=0;
o_msJump[6][8]=0;
o_msJump[7][0]=0;
o_msJump[7][1]=0;
o_msJump[7][2]=0;
o_msJump[7][3]=(std::sqrt(i_muL/i_rhoL)*std::sqrt(i_muR/i_rhoR))/(i_muR*std::sqrt(i_muL/i_rhoL) + i_muL*std::sqrt(i_muR/i_rhoR));
o_msJump[7][4]=0;
o_msJump[7][5]=0;
o_msJump[7][6]=0;
o_msJump[7][7]=-((i_muL*std::sqrt(i_muR/i_rhoR))/(i_muR*std::sqrt(i_muL/i_rhoL) + i_muL*std::sqrt(i_muR/i_rhoR)));
o_msJump[7][8]=0;
o_msJump[8][0]=0;
o_msJump[8][1]=0;
o_msJump[8][2]=0;
o_msJump[8][3]=0;
o_msJump[8][4]=0;
o_msJump[8][5]=(std::sqrt(i_muL/i_rhoL)*std::sqrt(i_muR/i_rhoR))/(i_muR*std::sqrt(i_muL/i_rhoL) + i_muL*std::sqrt(i_muR/i_rhoR));
o_msJump[8][6]=0;
o_msJump[8][7]=0;
o_msJump[8][8]=-((i_muL*std::sqrt(i_muR/i_rhoR))/(i_muR*std::sqrt(i_muL/i_rhoL) + i_muL*std::sqrt(i_muR/i_rhoR)));
