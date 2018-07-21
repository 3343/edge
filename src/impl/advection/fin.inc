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
 * Finalize for the advection equation.
 **/
edge::io::ErrorNorms l_errorWriter( l_config.m_errorNormsType,
                                    l_config.m_errorNormsFile );

if( l_errorWriter.outEnabled() ) {
  double l_normsTmp[3][1][N_CRUNS];
  edge::setups::InitialDofs<
    T_SDISC.ELEMENT,
    ORDER,
    N_QUANTITIES,
    N_CRUNS >::err( 0,
                    l_internal.m_nElements,
                    0,
                    t_enTypeShared::LIMIT,
                    l_config.m_refValsExprStrs,
                    l_basis,
                    l_internal.m_connect.elVe,
                    l_internal.m_globalShared2[0].connect.scSv,
                    l_internal.m_vertexChars,
                    l_internal.m_globalShared2[0].svChars,
                    l_internal.m_elementChars,
                    l_internal.m_globalShared2[0].adm[0], // dummy
                    l_internal.m_elementModePrivate1,
                    l_internal.m_globalShared2[0].dofs,
                    l_normsTmp[0], l_normsTmp[1], l_normsTmp[2] );

  // perform reductions
  double l_norms[3][1][N_CRUNS];
#ifdef PP_USE_MPI
  MPI_Allreduce( &l_normsTmp[0][0][0], &l_norms[0][0][0], N_CRUNS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
  MPI_Allreduce( &l_normsTmp[1][0][0], &l_norms[1][0][0], N_CRUNS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
  MPI_Allreduce( &l_normsTmp[2][0][0], &l_norms[2][0][0], N_CRUNS, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
#else
  for( unsigned short l_cr = 0; l_cr < N_CRUNS; l_cr++ ) {
    l_norms[0][0][l_cr] = l_normsTmp[0][0][l_cr];
    l_norms[1][0][l_cr] = l_normsTmp[1][0][l_cr];
    l_norms[2][0][l_cr] = l_normsTmp[2][0][l_cr];
  }
#endif

  // compute l2
  for( unsigned short l_cr = 0; l_cr < N_CRUNS; l_cr++ )
    l_norms[1][0][l_cr] = std::sqrt( l_norms[1][0][l_cr] );

  // write
  l_errorWriter.write( l_norms );
}
