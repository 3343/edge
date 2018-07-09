/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2016-2018, Regents of the University of California
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
 * Steps for the shallow water equations.
 **/
if( i_step == 0 ) {
  edge::swe::solvers::FiniteVolume<
    T_SDISC.ELEMENT,
    N_CRUNS
  >::netUpdates( i_first,
                 i_size,
                 m_internal.m_connect.faEl,
                 m_internal.m_faceChars,
                 m_internal.m_elementModePrivate1,
                 m_internal.m_elementModeShared1,
                 m_internal.m_faceModePrivate1 );
}
else if( i_step == 1 ) {
  edge::swe::solvers::FiniteVolume<
    T_SDISC.ELEMENT,
    N_CRUNS
  >::update( i_first,
             i_size,
             real_base(m_dT),
             m_internal.m_connect.faEl,
             m_internal.m_connect.elFa,
             m_internal.m_faceChars,
             m_internal.m_elementChars,
             m_internal.m_faceModePrivate1,
             m_internal.m_elementModePrivate1 );
}
else EDGE_LOG_FATAL << "step not supported";
