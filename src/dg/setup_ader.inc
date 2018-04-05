/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *         Alexander Heinecke (alexander.heinecke AT intel.com)
 *
 * @section LICENSE
 * Copyright (c) 2016-2017, Regents of the University of California
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
edge::dg::QuadratureEval<T_SDISC.ELEMENT, ORDER>::line( l_internal.m_globalShared1[0].quadEval.ptsLine,
                                                        l_internal.m_globalShared1[0].quadEval.weightsLine );

edge::dg::QuadratureEval<T_SDISC.ELEMENT, ORDER>::faces( l_internal.m_globalShared1[0].quadEval.ptsFaces,
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

#if defined PP_T_KERNELS_VANILLA || defined PP_T_KERNELS_XSMM_DENSE_SINGLE
l_basis.getStiffMm1Dense( N_ELEMENT_MODES, (real_base *) l_internal.m_globalShared1[0].mat.stiffT, true  );
l_basis.getStiffMm1Dense( N_ELEMENT_MODES, (real_base *) l_internal.m_globalShared1[0].mat.stiff,  false );
l_basis.getFluxDense( (real_base *) l_internal.m_globalShared1[0].mat.fluxL,
                      (real_base *) l_internal.m_globalShared1[0].mat.fluxN,
                      (real_base *) l_internal.m_globalShared1[0].mat.fluxT );
#elif defined PP_T_KERNELS_XSMM
{
  // determine CSC fill-in strategy
  std::string l_cscFillIn = "none";
  if ( libxsmm_get_target_archid() == LIBXSMM_X86_AVX512_KNM ) {
    l_cscFillIn = "qfma";
  }
  EDGE_LOG_INFO << "  using the following CSC fill-in strategy for DG matrices: " << l_cscFillIn;

  // get matrices as dense
  real_base l_stiffTDe[N_DIM][N_ELEMENT_MODES][N_ELEMENT_MODES];
  real_base l_stiffDe[N_DIM][N_ELEMENT_MODES][N_ELEMENT_MODES];
  real_base l_fluxLDe[ C_ENT[T_SDISC.ELEMENT].N_FACES ][N_ELEMENT_MODES][N_FACE_MODES];
  real_base l_fluxNDe[ N_FLUXN_MATRICES               ][N_ELEMENT_MODES][N_FACE_MODES];
  real_base l_fluxTDe[ C_ENT[T_SDISC.ELEMENT].N_FACES ][N_FACE_MODES][N_ELEMENT_MODES];
  l_basis.getStiffMm1Dense( N_ELEMENT_MODES, (real_base *) l_stiffTDe, true  );
  l_basis.getStiffMm1Dense( N_ELEMENT_MODES, (real_base *) l_stiffDe,  false );
  l_basis.getFluxDense( (real_base *) l_fluxLDe,
                        (real_base *) l_fluxNDe,
                        (real_base *) l_fluxTDe );

  // non-zero entries of the matices
  std::vector< real_base > *l_globalMat = new std::vector< real_base >();

  // offsets of the matrices
  std::vector< unsigned int > l_offsets;
  l_offsets.push_back( 0 );

  /*
   * 1) Cauchy-Kovalewski
   */
  t_matCrd l_stiffTCrd[N_DIM];
  for( unsigned short l_di = 0; l_di < N_DIM; l_di++ ) {
    edge::linalg::Matrix::denseToCrd< real_base >( N_ELEMENT_MODES, N_ELEMENT_MODES,
                                                   l_stiffTDe[l_di][0], l_stiffTCrd[l_di], TOL.BASIS );
  }

  // nz-blocks
  unsigned int l_nzBl[2][2][2];
  // init with matrix dim
  l_nzBl[0][0][0] = l_nzBl[0][1][0] = 0;
  l_nzBl[0][0][1] = l_nzBl[0][1][1] = N_ELEMENT_MODES-1;

  for( unsigned int l_de = 1; l_de < ORDER; l_de++ ) {
    // determine non-zero block for the next iteration
    unsigned int l_maxNzCol = 0;

    for( unsigned short l_di = 0; l_di < N_DIM; l_di++ ) {
      edge::linalg::Matrix::getBlockNz( l_stiffTCrd[l_di], l_nzBl[0], l_nzBl[1] );
      l_maxNzCol = std::max( l_maxNzCol, l_nzBl[1][1][1] );
    }

    t_matCsc l_stiffTCsc[N_DIM];
    for( unsigned short l_di = 0; l_di < N_DIM; l_di++ ) {
      edge::linalg::Matrix::denseToCsc< real_base >( N_ELEMENT_MODES, N_ELEMENT_MODES,
                                                     l_stiffTDe[l_di][0], l_stiffTCsc[l_di], TOL.BASIS,
                                                     l_nzBl[0][0][1]+1, l_maxNzCol+1,
                                                     l_cscFillIn );
    }

    // add data for the first CK-stiff matrix or shrinking stiff matrices
    if( l_de == 1 || l_maxNzCol < l_nzBl[0][0][1] ) {
      for( unsigned short l_di = 0; l_di < N_DIM; l_di++ ) {
        for( unsigned int l_nz = 0; l_nz < l_stiffTCsc[l_di].val.size(); l_nz++ ) {
          l_globalMat->push_back( l_stiffTCsc[l_di].val[l_nz] );
        }
        l_offsets.push_back( l_offsets.back() + l_stiffTCsc[l_di].val.size() );
      }
    }
    else {
      for( unsigned short l_di = 0; l_di < N_DIM; l_di++ ) {
        l_offsets.push_back( l_offsets.back() );
      }
    }

#ifdef PP_T_BASIS_HIERARCHICAL
   // check that size goes down with the number of derivatives
   EDGE_CHECK_EQ( l_maxNzCol+1, CE_N_ELEMENT_MODES( T_SDISC.ELEMENT, ORDER-l_de ) );
#endif

    // reduce relevant rows due to generated zero block
    l_nzBl[0][0][1] = l_maxNzCol;
  }

  /*
   * 2) Volume
   */
  for( unsigned short l_di = 0; l_di < N_DIM; l_di++ ) {
    t_matCsc l_stiffCsc;

    edge::linalg::Matrix::denseToCsc< real_base >( N_ELEMENT_MODES, N_ELEMENT_MODES,
                                                   l_stiffDe[l_di][0], l_stiffCsc, TOL.BASIS,
                                                   std::numeric_limits< unsigned int >::max(), std::numeric_limits< unsigned int >::max(),
                                                   l_cscFillIn );

    for( unsigned short l_nz = 0; l_nz < l_stiffCsc.val.size(); l_nz++ ) {
      l_globalMat->push_back( l_stiffCsc.val[l_nz] );
    }
    l_offsets.push_back( l_offsets.back() + l_stiffCsc.val.size() );
  }

  /*
   * 3) Surface
   */
  // local contribution flux matrices
  for( unsigned short l_fl = 0; l_fl < C_ENT[T_SDISC.ELEMENT].N_FACES; l_fl++ ) {
    t_matCsc l_fluxCsc;

    edge::linalg::Matrix::denseToCsc< real_base >( N_ELEMENT_MODES, N_FACE_MODES,
                                                   l_fluxLDe[l_fl][0], l_fluxCsc, TOL.BASIS,
                                                   std::numeric_limits< unsigned int >::max(), std::numeric_limits< unsigned int >::max(),
                                                   l_cscFillIn );

    for( unsigned short l_nz = 0; l_nz < l_fluxCsc.val.size(); l_nz++ ) {
      l_globalMat->push_back( l_fluxCsc.val[l_nz] );
    }
    l_offsets.push_back( l_offsets.back() + l_fluxCsc.val.size() );
  }

  // neighboring contribution flux matrices
  for( unsigned short l_fn = 0; l_fn < N_FLUXN_MATRICES; l_fn++ ) {
    t_matCsc l_fluxCsc;

    edge::linalg::Matrix::denseToCsc< real_base >( N_ELEMENT_MODES, N_FACE_MODES,
                                                   l_fluxNDe[l_fn][0], l_fluxCsc, TOL.BASIS,
                                                   std::numeric_limits< unsigned int >::max(), std::numeric_limits< unsigned int >::max(),
                                                   l_cscFillIn );

    for( unsigned short l_nz = 0; l_nz < l_fluxCsc.val.size(); l_nz++ ) {
      l_globalMat->push_back( l_fluxCsc.val[l_nz] );
    }
    l_offsets.push_back( l_offsets.back() + l_fluxCsc.val.size() );
  }

  // transposed flux matrices
  for( unsigned short l_ft = 0; l_ft < C_ENT[T_SDISC.ELEMENT].N_FACES; l_ft++ ) {
    t_matCsc l_fluxCsc;

    edge::linalg::Matrix::denseToCsc< real_base >( N_FACE_MODES, N_ELEMENT_MODES,
                                                   l_fluxTDe[l_ft][0], l_fluxCsc, TOL.BASIS,
                                                   std::numeric_limits< unsigned int >::max(), std::numeric_limits< unsigned int >::max(),
                                                   l_cscFillIn );

    for( unsigned short l_nz = 0; l_nz < l_fluxCsc.val.size(); l_nz++ ) {
      l_globalMat->push_back( l_fluxCsc.val[l_nz] );
    }
    l_offsets.push_back( l_offsets.back() + l_fluxCsc.val.size() );
  }

  // set the pointers
  unsigned short l_mat = 0;
  for( unsigned int l_de = 1; l_de < ORDER; l_de++ ) {
    for( unsigned int l_di = 0; l_di < N_DIM; l_di++ ) {
      l_internal.m_globalShared1[0].mat.stiffT[(l_de-1)*N_DIM+l_di] = &(*l_globalMat)[0] + l_offsets[l_mat];
      l_mat++;
    }
  }

  for( unsigned int l_di = 0; l_di < N_DIM; l_di++ ) {
    l_internal.m_globalShared1[0].mat.stiff[l_di] = &(*l_globalMat)[0] + l_offsets[l_mat];
    l_mat++;
  }

  for( unsigned short l_fl = 0; l_fl < C_ENT[T_SDISC.ELEMENT].N_FACES; l_fl++ ) {
    l_internal.m_globalShared1[0].mat.fluxL[l_fl] = &(*l_globalMat)[0] + l_offsets[l_mat];
    l_mat++;
  }
  for( unsigned short l_fn = 0; l_fn < N_FLUXN_MATRICES; l_fn++ ) {
    l_internal.m_globalShared1[0].mat.fluxN[l_fn] = &(*l_globalMat)[0] + l_offsets[l_mat];
    l_mat++;
  }
  for( unsigned short l_ft = 0; l_ft < C_ENT[T_SDISC.ELEMENT].N_FACES; l_ft++ ) {
    l_internal.m_globalShared1[0].mat.fluxT[l_ft] = &(*l_globalMat)[0] + l_offsets[l_mat];
    l_mat++;
  }
}
#else
#error undefined kernels
#endif
