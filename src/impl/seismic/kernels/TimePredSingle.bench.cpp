/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (alex.breuer AT uni-jena.de)
 *
 * @section LICENSE
 * Copyright (c) 2022, Alexander Breuer
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
 * Performance benchmarks of the optimized time prediction for single forward simulations.
 **/
#define CATCH_CONFIG_ENABLE_BENCHMARKING
#include <catch.hpp>
#include "TimePredSingle.hpp"

TEST_CASE( "Optimized elastic ADER time prediction for single forward simulations (FP32).", "[elastic][TimePredSingleFp32]" ) {
  // set up kernel
  edge::data::Dynamic l_dynMem;
  edge::seismic::kernels::TimePredSingle< float,
                                          0,
                                          TET4,
                                          ORDER,
                                          ORDER > l_pred( nullptr, l_dynMem );

  // set up dummy data structures
  float l_scratch[9][N_ELEMENT_MODES][1]     = {0};
  float l_ders[ORDER][9][N_ELEMENT_MODES][1] = {0};
  float l_dofs[9][N_ELEMENT_MODES][1]        = {0};
  float l_tDofs[9][N_ELEMENT_MODES][1]       = {0};
  float l_star[3][81]                        = {0};

  // benchmark time prediction
  BENCHMARK("elastic TimePredSingle FP32") {
    l_pred.ck( 0.017,
               l_star,
               nullptr,
               nullptr,
               l_dofs,
               nullptr,
               l_scratch,
               l_ders,
               nullptr,
               l_tDofs,
               nullptr );
  };
}


TEST_CASE( "Optimized viscoelastic ADER time prediction for single forward simulations (FP32).", "[visco][TimePredSinglFp32]" ) {
  // set up dummy data structures
  float l_rfs[PP_N_RELAXATION_MECHANISMS]                                 = {0};
  float l_srcA[PP_N_RELAXATION_MECHANISMS][36]                            = {0};
  float l_scratch[9][N_ELEMENT_MODES][1]                                  = {0};
  float l_dersE[ORDER][9][N_ELEMENT_MODES][1]                             = {0};
  float l_dofsE[9][N_ELEMENT_MODES][1]                                    = {0};
  float l_dofsA[PP_N_RELAXATION_MECHANISMS][6][N_ELEMENT_MODES][1]        = {0};
  float l_tDofsE[9][N_ELEMENT_MODES][1]                                   = {0};
  float l_tDofsA[PP_N_RELAXATION_MECHANISMS][6][N_ELEMENT_MODES][1]       = {0};
  float l_dersA[PP_N_RELAXATION_MECHANISMS][ORDER][6][N_ELEMENT_MODES][1] = {0};
  float l_starA[3][18]                                                    = {0};
  float l_starE[3][81]                                                    = {0};

  // set up kernel
  edge::data::Dynamic l_dynMem;
  edge::seismic::kernels::TimePredSingle< float,
                                          PP_N_RELAXATION_MECHANISMS,
                                          TET4,
                                          ORDER,
                                          ORDER > l_pred( l_rfs,
                                                          l_dynMem );

  // benchmark time prediction
  BENCHMARK("visco TimePredSingle FP32") {
    l_pred.ck( 0.017,
               l_starE,
               l_starA,
               l_srcA,
               l_dofsE,
               l_dofsA,
               l_scratch,
               l_dersE,
               l_dersA,
               l_tDofsE,
               l_tDofsA );
  };
}
