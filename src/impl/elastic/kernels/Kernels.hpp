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
 * Time, volume and surface kernels, based on the given build configuration.
 **/

#ifndef EDGE_SEISMIC_KERNELS_KERNELS_HPP
#define EDGE_SEISMIC_KERNELS_KERNELS_HPP

#include "constants.hpp"
#include "data/Dynamic.h"

#if defined(PP_T_KERNELS_VANILLA)
#include "TimePredVanilla.hpp"
#include "VolIntVanilla.hpp"
#include "SurfIntVanilla.hpp"
#elif defined(PP_T_KERNELS_XSMM_DENSE_SINGLE)
#include "TimePredSingle.hpp"
#include "VolIntSingle.hpp"
#include "SurfIntSingle.hpp"
#elif defined(PP_T_KERNELS_XSMM)
#include "TimePredFused.hpp"
#include "VolIntFused.hpp"
#include "SurfIntFused.hpp"
#else
#error kernels not supported
#endif


namespace edge {
  namespace elastic {
    namespace kernels {
      template< typename       TL_T_REAL,
                t_entityType   TL_T_EL,
                unsigned short TL_O_SP,
                unsigned short TL_O_TI,
                unsigned short TL_N_CRS >
      class Kernels;
    }
  }
}

/**
 * Time, volume and surface kernels.
 *
 * @paramt TL_T_REAL floating point precision.
 * @paramt TL_T_EL element type.
 * @paramt TL_O_SP order in space.
 * @paramt TL_O_TI order in time.
 * @paramt TL_N_CRS number of fused simulations. 
 **/
template< typename       TL_T_REAL,
          t_entityType   TL_T_EL,
          unsigned short TL_O_SP,
          unsigned short TL_O_TI,
          unsigned short TL_N_CRS >
class edge::elastic::kernels::Kernels {
  public:
#if defined(PP_T_KERNELS_VANILLA)
    TimePredVanilla< TL_T_REAL,
                     TL_T_EL,
                     TL_O_SP,
                     TL_O_TI,
                     TL_N_CRS > m_time;
    VolIntVanilla< TL_T_REAL,
                   TL_T_EL,
                   TL_O_SP,
                   TL_N_CRS > m_volInt;
    SurfIntVanilla< TL_T_REAL,
                    TL_T_EL,
                    TL_O_SP,
                    TL_N_CRS > m_surfInt;
#elif defined(PP_T_KERNELS_XSMM_DENSE_SINGLE)
    static_assert( TL_N_CRS == 1, "trying to build single kernels in fused setting" );
    TimePredSingle< TL_T_REAL,
                    TL_T_EL,
                    TL_O_SP,
                    TL_O_TI > m_time;
    VolIntSingle< TL_T_REAL,
                  TL_T_EL,
                  TL_O_SP> m_volInt;
    SurfIntSingle< TL_T_REAL,
                   TL_T_EL,
                   TL_O_SP > m_surfInt;
#elif defined(PP_T_KERNELS_XSMM)
static_assert( TL_N_CRS != 1, "trying to build fused kernels in single setting" );
    TimePredFused< TL_T_REAL,
                   TL_T_EL,
                   TL_O_SP,
                   TL_O_TI,
                   TL_N_CRS > m_time;
    VolIntFused< TL_T_REAL,
                 TL_T_EL,
                 TL_O_SP,
                 TL_N_CRS > m_volInt;
    SurfIntFused< TL_T_REAL,
                  TL_T_EL,
                  TL_O_SP,
                  TL_N_CRS > m_surfInt;
#else
#error kernels not supported
#endif

    /**
     * Constructor, which initializes the kernels.
     *
     * @param io_dynMem dynamic memory allocations.
     **/
    Kernels( data::Dynamic & io_dynMem ): m_time( io_dynMem ),
                                          m_volInt( io_dynMem ),
                                          m_surfInt( io_dynMem ) {};
};

#endif