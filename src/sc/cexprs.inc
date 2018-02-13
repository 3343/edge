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
 * Constant expressions for the sub-cell limiter.
 **/

/**
 * Derives the number of sub-vertices used for sub-cell limiting
 *
 * @param i_et element type.
 * @param i_os order of the DG-scheme in space.
 **/
constexpr unsigned short CE_N_SUB_VERTICES( t_entityType i_et,
                                            unsigned     i_os ) {
  return (i_et == LINE   ) ?    2 * i_os :
         (i_et == QUAD4R ) ?   (2 * i_os) * (2 * i_os) :
         (i_et == TRIA3  ) ? ( (2 * i_os) * (2 * i_os + 1 ) ) / 2 :
         (i_et == HEX8R  ) ?   (2 * i_os) * (2 * i_os) * (2 * i_os) :
         (i_et == TET4   ) ? ( (2 * i_os) * (2 * i_os + 1 ) * (2 * i_os + 2 ) ) / 6 :
                               (unsigned short) -1;
}

/**
 * Derives the number of sub-edges used for sub-cell limiting.
 *
 * @param i_dim number of dimensions.
 * @param i_os order of the DG-scheme in space.
 * @return number of sub-edges.
 **/
constexpr unsigned short CE_N_SUB_EDGES( unsigned short i_dim, unsigned i_os ) {
  return (i_dim == 1) ? 1 : 2 * (i_os - 1) + 1;
}

/**
 * Derives the number of sub-faces used for sub-cell limiting.
 *
 * @param i_et element type.
 * @param i_os order of the DG-scheme in space.
 * @return number of sub-faces.
 **/
constexpr unsigned short CE_N_SUB_FACES( t_entityType i_et,
                                         unsigned     i_os ) {
  return C_ENT[i_et].N_DIM == 1 ?  1 : (
           C_ENT[i_et].N_DIM == 2 ?  CE_N_SUB_EDGES(2, i_os) :
             CE_N_SUB_EDGES(3, i_os) *  CE_N_SUB_EDGES(3 , i_os)
         );
}

/**
 * Derives the number of sub-cells used for sub-cell limiting.
 *
 * @param i_en element type.
 * @param i_os order of the DG-scheme in space.
 * @return number of sub-cells.
 **/
constexpr unsigned short CE_N_SUB_CELLS( t_entityType i_et,
                                         unsigned     i_os ) {
  return C_ENT[i_et].N_DIM == 1 ?  2 * (i_os - 1) + 1 : (
           C_ENT[i_et].N_DIM == 2 ?  CE_N_SUB_EDGES(2, i_os) *  CE_N_SUB_EDGES(2, i_os) :
             CE_N_SUB_EDGES(3, i_os) *  CE_N_SUB_EDGES(3, i_os) *  CE_N_SUB_EDGES(3, i_os)
         );
}
