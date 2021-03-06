/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2019, Alexander Breuer
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
 * Third order reference solution for the time prediction unit tests.
 **/

#pragma GCC diagnostic push
#if !defined(__clang__)
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#endif

float l_utDt = 0.017;

float l_utStar[3] = { 0.7876563493948178,0.7019650033591671,0.7065805011608832 };
float l_utDofs[10] = { 0.31502566345387806,0.6686942615102741,0.5240257048836324,0.4324515724140714,0.600432130975785,0.8142123084940156,0.39985121899704656,0.1543684265191141,0.5628388222969362,0.026977345180984713, };
float l_utDers[3][10] = {
  { 0.31502566345387806,0.6686942615102741,0.5240257048836324,0.4324515724140714,0.600432130975785,0.8142123084940156,0.39985121899704656,0.1543684265191141,0.5628388222969362,0.026977345180984713, },
  { -5.6880746896948065,-7.067985114141075,-8.862840982715426,-3.2690887552582133,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0, },
  { 54.6574223802016,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, }
};
float l_utTint[10] = { 0.004578264805414016,0.010346478596681275,0.007627756461019372,0.006879293405904402,0.010207346226588344,0.013841609244398267,0.006797470722949792,0.0026242632508249398,0.009568259979047915,0.00045861486807674014, };

float l_utTint12[10] = { 0.002477830856112612,0.005428570260588983,0.004134048361010281,0.0035577425342359042,0.005103673113294172,0.006920804622199133,0.003398735361474896,0.0013121316254124699,0.004784129989523957,0.00022930743403837007, };
float l_utTint13[10] = { 0.0016954778314748013,0.0036757870542256216,0.002827181158562542,0.0023980718742203145,0.003402448742196115,0.004613869748132755,0.0022658235743165975,0.0008747544169416465,0.0031894199930159717,0.00015287162269224673, };

#pragma GCC diagnostic pop