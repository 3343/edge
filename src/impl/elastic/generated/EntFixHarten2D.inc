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
 * Generated Harten's entropy fix in 2D.
 * Adds viscosity to the zero-wave speed component of the Riemann solution.
 * For i_scale=1.0, the wave speed is scale with c_p.
 **/
o_entFixHarten[0][0]=0.;
o_entFixHarten[0][1]=0.;
o_entFixHarten[0][2]=0.;
o_entFixHarten[0][3]=0.;
o_entFixHarten[0][4]=0.;
o_entFixHarten[1][0]=(-0.25*(i_lamR*std::sqrt((i_lamL + 2*i_muL)/i_rhoL) + i_lamL*std::sqrt((i_lamR + 2*i_muR)/i_rhoR))*i_scale*std::max(std::sqrt((i_lamL + 2*i_muL)/i_rhoL),std::sqrt((i_lamR + 2*i_muR)/i_rhoR)))/(i_lamR*std::sqrt((i_lamL + 2*i_muL)/i_rhoL) + 2*i_muR*std::sqrt((i_lamL + 2*i_muL)/i_rhoL) + (i_lamL + 2*i_muL)*std::sqrt((i_lamR + 2*i_muR)/i_rhoR));
o_entFixHarten[1][1]=0.25*i_scale*std::max(std::sqrt((i_lamL + 2*i_muL)/i_rhoL),std::sqrt((i_lamR + 2*i_muR)/i_rhoR));
o_entFixHarten[1][2]=0.;
o_entFixHarten[1][3]=(0.25*(2*i_lamR*i_muL - 2*i_lamL*i_muR)*i_scale*std::max(std::sqrt((i_lamL + 2*i_muL)/i_rhoL),std::sqrt((i_lamR + 2*i_muR)/i_rhoR)))/(i_lamR*std::sqrt((i_lamL + 2*i_muL)/i_rhoL) + 2*i_muR*std::sqrt((i_lamL + 2*i_muL)/i_rhoL) + (i_lamL + 2*i_muL)*std::sqrt((i_lamR + 2*i_muR)/i_rhoR));
o_entFixHarten[1][4]=0.;
o_entFixHarten[2][0]=0.;
o_entFixHarten[2][1]=0.;
o_entFixHarten[2][2]=0.;
o_entFixHarten[2][3]=0.;
o_entFixHarten[2][4]=0.;
o_entFixHarten[3][0]=0.;
o_entFixHarten[3][1]=0.;
o_entFixHarten[3][2]=0.;
o_entFixHarten[3][3]=0.;
o_entFixHarten[3][4]=0.;
o_entFixHarten[4][0]=0.;
o_entFixHarten[4][1]=0.;
o_entFixHarten[4][2]=0.;
o_entFixHarten[4][3]=0.;
o_entFixHarten[4][4]=0.;
