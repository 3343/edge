/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *         Alexander Heinecke (alexander.heinecke AT intel.com)
 *
 * @section LICENSE
 * Copyright (c) 2016-2019, Regents of the University of California
 * Copyright (c) 2016, Intel Corporation
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
 * Setup for ADER-DG simulations.
 **/
// set up quadrature points
edge::dg::QuadratureEval<T_SDISC.ELEMENT, ORDER, ORDER>::line( l_internal.m_globalShared1[0].quadEval.ptsLine,
                                                               l_internal.m_globalShared1[0].quadEval.weightsLine );

edge::dg::QuadratureEval<T_SDISC.ELEMENT, ORDER, ORDER>::faces( l_internal.m_globalShared1[0].quadEval.ptsFaces,
                                                                l_internal.m_globalShared1[0].quadEval.weightsFaces,
                                                                l_internal.m_globalShared1[0].quadEval.basisFaces );

// get inverse mass matrix
{
  real_base l_massI[N_ELEMENT_MODES][N_ELEMENT_MODES];
  l_basis.getMassInvDense( N_ELEMENT_MODES, (real_base *) l_massI, true );
  for( int_md l_md1 = 0; l_md1 < N_ELEMENT_MODES; l_md1++ ) {
    for( int_md l_md2 = 0; l_md2 < N_ELEMENT_MODES; l_md2++ ) {
      if( l_md1 != l_md2 ) {
        EDGE_CHECK( std::abs(l_massI[l_md1][l_md2]) < TOL.BASIS );
      }
      else {
        l_internal.m_globalShared1[0].mat.massI[l_md1] = l_massI[l_md1][l_md2];
      }
    }
  }
}

#if defined PP_T_KERNELS_VANILLA
l_basis.getStiffMm1Dense( N_ELEMENT_MODES, (real_base *) l_internal.m_globalShared1[0].mat.stiffT, true  );
l_basis.getStiffMm1Dense( N_ELEMENT_MODES, (real_base *) l_internal.m_globalShared1[0].mat.stiff,  false );
l_basis.getFluxDense( (real_base *) l_internal.m_globalShared1[0].mat.fluxL,
                      (real_base *) l_internal.m_globalShared1[0].mat.fluxN,
(real_base *) l_internal.m_globalShared1[0].mat.fluxT );
#endif