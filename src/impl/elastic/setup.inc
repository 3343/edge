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
 * Setup for elastics.
 **/
#ifdef PP_USE_MPI
  // init mpi layout
  EDGE_CHECK( l_enLayouts[2].timeGroups.size() == 1 );

  l_mpi.initLayout( l_enLayouts[2],
                    l_internal.m_elementModePrivate2[0][0][0],
                    N_QUANTITIES*N_ELEMENT_MODES*N_CRUNS*sizeof(t_elementModePrivate2),
                    0,
                    1,
                    100 );
#endif

// parse config specific to elastics
edge::elastic::io::Config l_elasticConf( l_config.m_doc );

// setup initial DOFs
EDGE_LOG_INFO << "  setting up material parameters and initial DOFs";
for( int_cfr l_run = 0; l_run < N_CRUNS; l_run++ ) {
  if( l_config.m_setups[l_run] == "plane_waves" ) {
    l_internal.m_globalPrivate1[0][l_run].initialDofs = setupChars::InitialDofs::planeWaves;
  }
  else l_internal.m_globalPrivate1[0][l_run].initialDofs = setupChars::InitialDofs::zero;
}

// get the initial setup
PP_INSTR_REG_DEF(dofsMat)
PP_INSTR_REG_BEG(dofsMat,"dofs_mat")

// query velocity model from mesh
unsigned short l_vmMesh = std::numeric_limits< unsigned short >::max();
EDGE_LOG_INFO << "    querying mesh for velocity model";

#if defined PP_T_MESH_UNSTRUCTURED
// allocate memory for velocity model
double (*l_velMod)[3] = (double (*)[3]) new double[ l_internal.m_nElements * 3 ];
std::string l_bgPars[3] = { "LAMBDA",
                            "MU",
                            "RHO" };

l_vmMesh = l_mesh.getParsDe( N_DIM,
                             3,
                             l_bgPars,
                             l_velMod[0] );

if( l_vmMesh == 0 ) {
#ifdef PP_USE_OMP
#pragma omp parallel for
#endif
  for( int_el l_el = 0; l_el < l_internal.m_nElements; l_el++ ) {
    // check for valid lame parameters
    EDGE_CHECK_GT( l_velMod[l_el][0], 0 ) << l_el;
    EDGE_CHECK_GT( l_velMod[l_el][1], 0 ) << l_el;
    EDGE_CHECK_GT( l_velMod[l_el][2], 0 ) << l_el;

    l_internal.m_elementShared1[l_el][0].lam  = l_velMod[l_el][0];
    l_internal.m_elementShared1[l_el][0].mu   = l_velMod[l_el][1];
    l_internal.m_elementShared1[l_el][0].rho  = l_velMod[l_el][2];
  }
}

// free memory of velocity model
delete[] l_velMod;
#endif

if( l_vmMesh == 0 ) {
  EDGE_LOG_INFO << "    successfully obtained velocity model from mesh";
}
else {
  EDGE_LOG_INFO << "    failed obtaining velocity model from mesh, continuining w/o "
                << l_vmMesh;
}

// perform NUMA-aware zero-init of DOFs and tDOFs
#ifdef PP_USE_OMP
#pragma omp parallel for num_threads( std::max( edge::parallel::g_nThreads-1, 1 ) )
#endif
for( int_el l_el = 0; l_el < l_internal.m_nElements; l_el++ ) {
  for( int_qt l_qt = 0; l_qt < N_QUANTITIES; l_qt++ ) {
    for( int_md l_md = 0; l_md < N_ELEMENT_MODES; l_md++ ) {
      for( int_cfr l_ru = 0; l_ru < N_CRUNS; l_ru++ ) {
        l_internal.m_elementModePrivate1[l_el][l_qt][l_md][l_ru] = 0;
        l_internal.m_elementModePrivate2[l_el][l_qt][l_md][l_ru] = 0;
      }
    }
  }
}


EDGE_LOG_INFO << "    setting initial DOFs and velocity model based on user-provided config (if available)";
for( int_cfr l_run = 0; l_run < N_CRUNS; l_run++ ) {
  if( l_internal.m_globalPrivate1[0][l_run].initialDofs == setupChars::InitialDofs::planeWaves ) {
    // setup dofs
    edge::elastic::setups::Convergence::setPlaneWaves(  l_run,
                                                        l_basis,
                                                        l_internal.m_nElements,
                                                        l_internal.m_connect.elVe,
                                                        l_internal.m_vertexChars,
                                                        l_internal.m_elementChars,
                                                        l_internal.m_elementShared1,
                                                        l_internal.m_elementModePrivate1,
                                                        -50+l_run*5,
                                                        -50+l_run*5,
                                                        -50+l_run*5 );
  }
  else if( l_internal.m_globalPrivate1[0][l_run].initialDofs == setupChars::InitialDofs::zero ) {
    // setup material parameters
    for( int_el l_el = 0; l_el < l_internal.m_nElements; l_el++ ) {
      // compute vertices average position in all dimension
      real_mesh l_ave[N_DIM];
      for( unsigned short l_dm = 0; l_dm < N_DIM; l_dm++ ) {
        l_ave[l_dm] = 0;
        for( unsigned int l_ve = 0; l_ve < C_ENT[T_SDISC.ELEMENT].N_VERTICES; l_ve++ ) {
          int_el l_veId = l_internal.m_connect.elVe[l_el][l_ve];
          l_ave[l_dm] += l_internal.m_vertexChars[l_veId].coords[l_dm];
        }
        l_ave[l_dm] /= C_ENT[T_SDISC.ELEMENT].N_VERTICES;
      }

      // check if velocity model is provided in user-config
      EDGE_CHECK( (l_vmMesh == 0) || (l_elasticConf.m_velDoms.size() > 0) )
        << "couldn't find a velocity model in the mesh or user config, aborting";

      // find the matching domain in the velocity model
      for( std::size_t l_do = 0; l_do < l_elasticConf.m_velDoms.size(); l_do++ ) {
        if( l_elasticConf.m_velDoms[l_do].inside(l_ave) ) {
          l_internal.m_elementShared1[l_el][0].rho = l_elasticConf.m_velVals[l_do][0];
          l_internal.m_elementShared1[l_el][0].lam = l_elasticConf.m_velVals[l_do][1];
          l_internal.m_elementShared1[l_el][0].mu  = l_elasticConf.m_velVals[l_do][2];
          break;
        }

        // abort if no matching velocity domain is present
        EDGE_CHECK( (l_vmMesh == 0) ||
                    (l_do != l_elasticConf.m_velDoms.size()-1 ) )
          << "here is the troublesome point: "
          << l_ave[0] << " " << l_ave[1] << " " << ( (N_DIM > 2) ? std::to_string(l_ave[2]) : "" );
      }
    }
  }
  else EDGE_LOG_FATAL;
}
PP_INSTR_REG_END(dofsMat)

// setup kinematic sources sources
if( l_elasticConf.m_kinSrcs.size() > 0 ) {
  PP_INSTR_REG_DEF(kin)
  PP_INSTR_REG_BEG(kin,"kin_srcs")
#ifdef PP_HAS_NETCDF
  EDGE_LOG_INFO << "  initializing kinematic sources";
  // setup reader
  edge::elastic::io::Nrf< N_DIM> l_srcReader;
  for( unsigned short l_ki = 0; l_ki < l_elasticConf.m_kinSrcs.size(); l_ki++ ) {
    l_srcReader.init( l_elasticConf.m_kinSrcs[l_ki] );
  }
  // print info of kinematic src descriptions
  std::vector< std::string > l_srcInfo = l_srcReader.toString();
  for( unsigned short l_ki = 0; l_ki < l_srcInfo.size(); l_ki++ ) {
    EDGE_LOG_INFO << "    " << l_srcInfo[l_ki];
  }

  // init sources
  edge::elastic::setups::KinematicsInit<
    int_el,
    real_mesh,
    real_base,
    T_SDISC.ELEMENT,
    ORDER,
    N_CRUNS,
    1 >::solvers(  l_enLayouts[2],
                   l_internal.m_connect.elVe,
                   t_spTypeElastic::SOURCE,
                   l_internal.m_vertexChars,
                  &l_gIdsEl[0],
                   l_srcReader,
                   l_internal.m_globalShared1[0].mat.massI,
                   l_internal.m_elementShared1[0],
                   l_internal.m_elementChars,
                   l_dynMem,
                   l_internal.m_globalShared4[0],
                   l_internal.m_globalShared3 );
#endif
  PP_INSTR_REG_END(kin)
}
// get layout of sparse source elements
l_enLayouts.resize( l_enLayouts.size() + 1 );
edge::data::SparseEntities::denseToSparse( t_spTypeElastic::SOURCE,
                                           l_internal.m_elementChars,
                                           l_enLayouts[2],
                                           l_enLayouts.back() );
unsigned short l_srcLayout = l_enLayouts.size()-1;

// set up dynamic rupture face layout
l_enLayouts.resize( l_enLayouts.size() + 1 );
unsigned short l_rupLayoutFa = l_enLayouts.size()-1;
edge::data::SparseEntities::denseToSparse( t_spTypeElastic::RUPTURE,
                                           l_internal.m_faceChars,
                                           l_enLayouts[1],
                                           l_enLayouts[l_rupLayoutFa] );

// set up dynamic rupture element layout
l_enLayouts.resize( l_enLayouts.size() + 1 );
unsigned short l_rupLayoutEl = l_enLayouts.size()-1;
edge::data::SparseEntities::denseToSparseAdj( C_ENT[T_SDISC.ELEMENT].N_FACES,
                                              l_internal.m_connect.elFa[0],
                                              t_spTypeElastic::RUPTURE,
                                              l_internal.m_faceChars,
                                              l_enLayouts[2],
                                              l_enLayouts[l_rupLayoutEl] );

// propagate sparse info of rupture faces to elements
edge::data::SparseEntities::propAdj( l_enLayouts[1].nEnts,
                                     2,
                                     l_internal.m_connect.faEl[0],
                                     t_spTypeElastic::RUPTURE,
                                     l_internal.m_faceChars,
                                     l_internal.m_elementChars );

// setup the GEMM kernels for ADER-DG
#if defined PP_T_KERNELS_XSMM
{
  EDGE_LOG_INFO << "  setting up libxsmm, sparse";

  // get matrices as dense
  real_base l_stiffTDense[N_DIM][N_ELEMENT_MODES][N_ELEMENT_MODES];
  real_base l_stiffDense[N_DIM][N_ELEMENT_MODES][N_ELEMENT_MODES];
  real_base l_fluxDense[N_FLUX_MATRICES][N_ELEMENT_MODES][N_ELEMENT_MODES];
  l_basis.getStiffMm1Dense( N_ELEMENT_MODES, (real_base *) l_stiffTDense, true  );
  l_basis.getStiffMm1Dense( N_ELEMENT_MODES, (real_base *) l_stiffDense,  false );
  l_basis.getFluxMm1Dense(  N_ELEMENT_MODES, (real_base *) l_fluxDense,   true  );

  real_base l_star[N_DIM][N_QUANTITIES][N_QUANTITIES];
  edge::elastic::solvers::AderDg::getJac( (real_base) 1.0,
                                          (real_base) 1.0,
                                          (real_base) 1.0,
                                                      l_star[0][0],
                                                      N_DIM );
  real_base l_fSolv[N_QUANTITIES][N_QUANTITIES];
  for( unsigned short l_di = 1; l_di < N_DIM; l_di++ ) {
    for( int_qt l_q1 = 0; l_q1 < N_QUANTITIES; l_q1++ ) {
      for( int_qt l_q2 = 0; l_q2 < N_QUANTITIES; l_q2++ ) {
        l_star[0][l_q1][l_q2] =   std::abs(l_star[0][l_q1][l_q2])
                                + std::abs(l_star[l_di][l_q1][l_q2]);
        l_fSolv[l_q1][l_q2] = 1;
      }
    }
  }

  /*
   * Derive sparse AoSoA-LIBXSMM kernels.
   *
   * 1) Cauchy Kovalewski
   */
  // get sparse, transposed stiffness matrices
  t_matCrd l_stiffTCrd[N_DIM];
  for( unsigned short l_di = 0; l_di < N_DIM; l_di++ ) {
    edge::linalg::Matrix::denseToCrd< real_base >( N_ELEMENT_MODES, N_ELEMENT_MODES,
                                                   l_stiffTDense[l_di][0], l_stiffTCrd[l_di], TOL.BASIS );
  }

  // get csr star matrix
  t_matCsr l_starCsr;
  edge::linalg::Matrix::denseToCsr< real_base >( N_QUANTITIES, N_QUANTITIES,
                                                 l_star[0][0],  l_starCsr,  TOL.BASIS );
  EDGE_CHECK( l_starCsr.val.size() == N_MAT_STAR );

  // exploit potential zero-block generation in recursive CK

  // nz-blocks
  unsigned int l_nzBl[2][2][2];
  // init with matrix dim
  l_nzBl[0][0][0] = l_nzBl[0][1][0] = 0;
  l_nzBl[0][0][1] = l_nzBl[0][1][1] = N_ELEMENT_MODES-1;

  // iterate over derivatives (recusive calls)
  for( unsigned short l_de = 1; l_de < ORDER; l_de++ ) {
    // determine non-zero block in the next iteration
    unsigned int l_maxNzCol = 0;
    for( unsigned short l_di = 0; l_di < N_DIM; l_di++ ) {
      edge::linalg::Matrix::getBlockNz( l_stiffTCrd[l_di], l_nzBl[0], l_nzBl[1] );
      l_maxNzCol = std::max( l_maxNzCol, l_nzBl[1][1][1] );
    }

    // generate libxsmm kernel for transposed stiffness matrices
    for( unsigned short l_di = 0; l_di < N_DIM; l_di++ ) {
      t_matCsr l_stiffTCsr;
      edge::linalg::Matrix::denseToCsr< real_base >( N_ELEMENT_MODES, N_ELEMENT_MODES,
                                                     l_stiffTDense[l_di][0], l_stiffTCsr, TOL.BASIS,
                                                     l_nzBl[0][0][1]+1, l_maxNzCol+1 );

      l_internal.m_mm.add( &l_stiffTCsr.rowPtr[0],  &l_stiffTCsr.colIdx[0], &l_stiffTCsr.val[0],
                            N_QUANTITIES, l_maxNzCol+1, l_nzBl[0][0][1]+1,
                            N_ELEMENT_MODES, 0, l_maxNzCol+1,
                            real_base(1.0), real_base(0.0),
                            LIBXSMM_PREFETCH_NONE );
    }
    // generate libxsmm kernel for star matrix
    l_internal.m_mm.add( &l_starCsr.rowPtr[0],  &l_starCsr.colIdx[0], &l_starCsr.val[0],
                          N_QUANTITIES, l_maxNzCol+1, N_QUANTITIES,
                          0, l_maxNzCol+1, N_ELEMENT_MODES,
                          real_base(1.0), real_base(1.0),
                          LIBXSMM_PREFETCH_NONE );

#ifdef PP_T_BASIS_HIERARCHICAL
   // check that size goes down with the number of derivatives
   EDGE_CHECK_EQ( l_maxNzCol+1, CE_N_ELEMENT_MODES( T_SDISC.ELEMENT, ORDER-l_de ) );
#endif

    // reduce relevant rows due to generated zero block
    l_nzBl[0][0][1] = l_maxNzCol;
  }

  /*
   * 2) add volume kernels
   */
  // get stiffness matrices
  t_matCrd l_stiffCrd[3];
  for( unsigned short l_di = 0; l_di < N_DIM; l_di++ ) {
    edge::linalg::Matrix::denseToCrd< real_base >( N_ELEMENT_MODES, N_ELEMENT_MODES,
                                                   l_stiffDense[l_di][0], l_stiffCrd[l_di], TOL.BASIS );
  }

  // reset zero-blocks
  l_nzBl[0][0][0] = l_nzBl[0][1][0] = 0;
  l_nzBl[0][0][1] = l_nzBl[0][1][1] = N_ELEMENT_MODES-1;
  unsigned int l_maxNzRow = 0;
  // get max #nz-rows
  for( unsigned short l_di = 0; l_di < N_DIM; l_di++ ) {
    edge::linalg::Matrix::getBlockNz( l_stiffCrd[l_di], l_nzBl[0], l_nzBl[1] );
    l_maxNzRow = std::max( l_maxNzRow, l_nzBl[1][0][1] );
  }

#ifdef PP_T_BASIS_HIERARCHICAL
   // check that size is one "order" less
   EDGE_CHECK_EQ( l_maxNzRow+1, CE_N_ELEMENT_MODES( T_SDISC.ELEMENT, ORDER-1 ) );
#endif

  for( unsigned short l_di = 0; l_di < N_DIM; l_di++ ) {
    t_matCsr l_stiffCsr;
    edge::linalg::Matrix::denseToCsr< real_base >( N_ELEMENT_MODES, N_ELEMENT_MODES,
                                                   l_stiffDense[l_di][0],  l_stiffCsr, TOL.BASIS );

    l_internal.m_mm.add( &l_stiffCsr.rowPtr[0], &l_stiffCsr.colIdx[0], &l_stiffCsr.val[0],
                          N_QUANTITIES, N_ELEMENT_MODES, l_maxNzRow+1,
                          l_maxNzRow+1, 0, N_ELEMENT_MODES,
                          real_base(1.0), real_base(1.0),
                          LIBXSMM_PREFETCH_NONE ); // Remark: Star matrix is multiplied first
  }

  // star matrix
  l_internal.m_mm.add( &l_starCsr.rowPtr[0], &l_starCsr.colIdx[0], &l_starCsr.val[0],
                        N_QUANTITIES, l_maxNzRow+1, N_QUANTITIES,
                        0, N_ELEMENT_MODES, l_maxNzRow+1,
                        real_base(1.0), real_base(0.0),
                        LIBXSMM_PREFETCH_NONE );

  /*
   * 3) surface kernels
   */
  // flux matrices
  for( unsigned short l_fm = 0; l_fm < N_FLUX_MATRICES; l_fm++ ) {
    t_matCsr l_fluxCsr;

    edge::linalg::Matrix::denseToCsr< real_base >( N_ELEMENT_MODES, N_ELEMENT_MODES,
                                                   l_fluxDense[l_fm][0], l_fluxCsr, TOL.BASIS );

    l_internal.m_mm.add( &l_fluxCsr.rowPtr[0],  &l_fluxCsr.colIdx[0], &l_fluxCsr.val[0],
                          N_QUANTITIES, N_ELEMENT_MODES, N_ELEMENT_MODES,
                          N_ELEMENT_MODES, 0, N_ELEMENT_MODES,
                          real_base(1.0), real_base(0.0),
                          LIBXSMM_PREFETCH_NONE );
  }

  // flux solver with and without prefetches
  t_matCsr l_fSolvCsr;
  edge::linalg::Matrix::denseToCsr< real_base >( N_QUANTITIES, N_QUANTITIES,
                                                 l_fSolv[0], l_fSolvCsr, TOL.BASIS );
  EDGE_CHECK( l_fSolvCsr.val.size() == N_QUANTITIES*N_QUANTITIES );

  l_internal.m_mm.add( &l_fSolvCsr.rowPtr[0],  &l_fSolvCsr.colIdx[0], &l_fSolvCsr.val[0],
                        N_QUANTITIES, N_ELEMENT_MODES, N_QUANTITIES,
                        0, N_ELEMENT_MODES, N_ELEMENT_MODES,
                        real_base(1.0), real_base(1.0),
                        LIBXSMM_PREFETCH_NONE );

  l_internal.m_mm.add( &l_fSolvCsr.rowPtr[0],  &l_fSolvCsr.colIdx[0], &l_fSolvCsr.val[0],
                        N_QUANTITIES, N_ELEMENT_MODES, N_QUANTITIES,
                        0, N_ELEMENT_MODES, N_ELEMENT_MODES,
                        real_base(1.0), real_base(1.0),
                        LIBXSMM_PREFETCH_BL2_VIA_C );
}
#endif

#ifndef PP_T_KERNELS_XSMM
edge::elastic::setups::MmKernels::add( T_SDISC.ELEMENT,
                                       ORDER,
                                       N_QUANTITIES,
                                       N_CRUNS,
                                       l_internal.m_mm );
#endif

// set up fault receivers
if( l_elasticConf.m_frictionLaw != "" &&
    l_config.m_recvCrds[1].size() > 0 ) {
  EDGE_LOG_INFO << "  setting up fault receivers";

  // TODO: fix dimension incompability of recv implementations
  std::vector< std::array< real_mesh, N_DIM > > l_crds( l_config.m_recvCrds[1].size() );
  for( std::size_t l_re = 0; l_re < l_config.m_recvCrds[1].size(); l_re++ ) {
    for( unsigned short l_di = 0; l_di < N_DIM; l_di++ )
      l_crds[l_re][l_di] = l_config.m_recvCrds[1][l_re][l_di];
  }

  l_recvsQuad.init(                     l_config.m_recvCrds[1].size(),
                                        t_spTypeElastic::RUPTURE,
                                        (N_DIM-1)*3, // TODO: hardcoded to linear slip weakening
                                        l_config.m_recvPath[1],
                                       &l_config.m_recvNames[1][0],
                (real_mesh (*)[N_DIM]) &l_crds[0][0], // TODO: dimension
                                        l_config.m_recvFreq[1],
                                        l_internal.m_globalShared1[0].quadEval.ptsFaces,
                                        l_enLayouts[1],
                                        l_enLayouts[2],
                                        l_internal.m_connect.faEl,
                                        l_internal.m_connect.elVe,
                                        l_internal.m_connect.elFa,
                                        l_internal.m_vertexChars,
                                        l_internal.m_faceChars );
   l_recvsQuad.print();
}

// init sparse, internal data structures
l_internal.initSparse( 0, l_enLayouts[l_rupLayoutFa].nEnts, 0,
                       0, l_enLayouts[l_rupLayoutFa].nEnts, l_enLayouts[l_rupLayoutEl].nEnts,
                       0, l_enLayouts[l_rupLayoutFa].nEnts, l_enLayouts[l_rupLayoutEl].nEnts,
                       0, l_enLayouts[l_rupLayoutFa].nEnts, 0,
                       0, l_enLayouts[l_rupLayoutFa].nEnts, 0,
                       0, l_enLayouts[l_rupLayoutFa].nEnts, 0  );

if( l_elasticConf.m_frictionLaw != "" ) {
  EDGE_LOG_INFO << "  setting up rupture physics";
  // link sparse rupture faces and sparse rupture elements
  edge::data::SparseEntities::linkSpAdj( l_enLayouts[1].nEnts,
                                         2,
                                         l_internal.m_connect.faEl[0],
                                         t_spTypeElastic::RUPTURE,
                                         l_internal.m_faceChars,
                                         l_internal.m_elementChars,
                                         l_internal.m_faceSparseShared3[0] );

  // link sparse rupture elements and sparse faces
  edge::data::SparseEntities::linkSpAdj( l_enLayouts[2].nEnts,
                                         C_ENT[T_SDISC.ELEMENT].N_FACES,
                                         l_internal.m_connect.elFa[0],
                                         t_spTypeElastic::RUPTURE,
                                         l_internal.m_faceChars,
                                         l_internal.m_elementSparseShared2[0][0] );

  edge::elastic::solvers::InternalBoundaryTypes<
    T_SDISC.ELEMENT
  >::initFaces( l_enLayouts[1].nEnts,
                RUPTURE,
                l_internal.m_faceChars,
                l_internal.m_connect.faEl,
                l_internal.m_connect.elFa,
                l_internal.m_connect.vIdElFaEl,
                l_internal.m_faceSparseShared6[0] );

  // set up rupture solvers
  edge::elastic::solvers::InternalBoundarySolvers<
    T_SDISC.ELEMENT
  >::init( l_enLayouts[1].nEnts,
           t_spTypeElastic::RUPTURE,
           l_internal.m_connect.faEl,
           l_internal.m_connect.elVe,
           l_internal.m_vertexChars,
           l_internal.m_faceChars,
           l_internal.m_elementShared1[0],
           l_internal.m_faceSparseShared5,
           l_elasticConf.m_faultCrds );

  // init rupture properties
  edge::elastic::setups::RuptureInit<
    T_SDISC.ELEMENT,
    ORDER,
    N_CRUNS >::linSlipWeak(
      l_enLayouts[1].nEnts,
      t_spTypeElastic::RUPTURE,
      l_internal.m_connect.faVe,
      l_internal.m_connect.faEl,
      l_internal.m_vertexChars,
      l_internal.m_faceChars,
      l_internal.m_elementShared1[0],
      l_elasticConf.m_faultCrds,
      l_elasticConf.m_lsw,
      l_elasticConf.m_rupDoms,
      l_elasticConf.m_stressInit,
      l_internal.m_globalShared2[0],
      l_internal.m_faceSparseShared1[0],
      l_internal.m_faceSparseShared2[0] );
}

#if PP_ORDER > 1
// setup star matrices
edge::elastic::solvers::AderDg::setupStarM( l_internal.m_nElements,
                                            l_internal.m_vertexChars,
                                            l_internal.m_connect.elVe,
                                            l_internal.m_elementShared1,
                                            l_internal.m_elementShared4 );
#endif

// setup solvers
edge::elastic::solvers::common::setupSolvers( l_internal.m_nElements,
                                              l_internal.m_nFaces,
                                              l_mesh.getInMap()->elMeDa,
                                              l_mesh.getInMap()->elDaMe,
                                              l_internal.m_connect.elVe,
                                              l_internal.m_connect.faEl,
                                              l_internal.m_connect.elFa,
                                              l_internal.m_vertexChars,
                                              l_internal.m_faceChars,
                                              l_internal.m_elementChars,
                                              l_internal.m_elementShared1,
                                              l_internal.m_elementShared2,
                                              l_internal.m_elementShared3 );

edge::elastic::common::getTimeStepStatsCFL( l_internal.m_nElements,
                                            l_internal.m_elementChars,
                                            l_internal.m_elementShared1,
                                            l_dT[0], l_dT[1], l_dT[2] );

// setup shared memory parallelization
for( int_tg l_tg = 0; l_tg < l_enLayouts[2].timeGroups.size(); l_tg++ ) {
  int_spType l_spType[3] = { RECEIVER, SOURCE, RUPTURE };

  // local inner-elements
  l_shared.regWrkRgn( l_tg,
                      0,
                      l_tg * N_ENTRIES_CONTROL_FLOW + 0,
                      l_enLayouts[2].timeGroups[l_tg].inner.first,
                      l_enLayouts[2].timeGroups[l_tg].inner.size,
                      l_tg,
                      3, l_spType, l_internal.m_elementChars );

  // local send-elements
  l_shared.regWrkRgn( l_tg,
                      0,
                      l_tg * N_ENTRIES_CONTROL_FLOW + 1,
                      l_enLayouts[2].timeGroups[l_tg].inner.first+
                      l_enLayouts[2].timeGroups[l_tg].inner.size,
                      l_enLayouts[2].timeGroups[l_tg].nEntsOwn-
                      l_enLayouts[2].timeGroups[l_tg].inner.size,
                      l_enLayouts[2].timeGroups.size() + l_tg,
                      3, l_spType, l_internal.m_elementChars );

  // neigh inner-elements
  l_shared.regWrkRgn( l_tg,
                      1,
                      l_tg * N_ENTRIES_CONTROL_FLOW + 3,
                      l_enLayouts[2].timeGroups[l_tg].inner.first,
                      l_enLayouts[2].timeGroups[l_tg].inner.size,
                      l_tg,
                      1, l_spType+2, l_internal.m_elementChars );

  // neigh send-elements
  l_shared.regWrkRgn( l_tg,
                      1,
                      l_tg * N_ENTRIES_CONTROL_FLOW + 4,
                      l_enLayouts[2].timeGroups[l_tg].inner.first+
                      l_enLayouts[2].timeGroups[l_tg].inner.size,
                      l_enLayouts[2].timeGroups[l_tg].nEntsOwn-
                      l_enLayouts[2].timeGroups[l_tg].inner.size,
                      l_enLayouts[2].timeGroups.size() + l_tg,
                      1, l_spType+2, l_internal.m_elementChars );

  if( l_elasticConf.m_frictionLaw == "" ) {
    // src inner-elements
    l_shared.regWrkRgn( l_tg,
                        2,
                        l_tg * N_ENTRIES_CONTROL_FLOW + 6,
                        l_enLayouts[l_srcLayout].timeGroups[l_tg].inner.first,
                        l_enLayouts[l_srcLayout].timeGroups[l_tg].inner.size,
                        l_tg );

    // src send-elements
    l_shared.regWrkRgn( l_tg,
                        2,
                        l_tg * N_ENTRIES_CONTROL_FLOW + 7,
                        l_enLayouts[l_srcLayout].timeGroups[l_tg].inner.first+
                        l_enLayouts[l_srcLayout].timeGroups[l_tg].inner.size,
                        l_enLayouts[l_srcLayout].timeGroups[l_tg].nEntsOwn-
                        l_enLayouts[l_srcLayout].timeGroups[l_tg].inner.size,
                        l_enLayouts[l_srcLayout].timeGroups.size() + l_tg );
  }
  else {
    l_spType[0] = { RECEIVER };

    // rupture physics inner-faces
    l_shared.regWrkRgn( l_tg,
                        2,
                        l_tg * N_ENTRIES_CONTROL_FLOW + 6,
                        l_enLayouts[l_rupLayoutFa].timeGroups[l_tg].inner.first,
                        l_enLayouts[l_rupLayoutFa].timeGroups[l_tg].inner.size,
                        l_tg,
                        1, l_spType, l_internal.m_faceSparseShared6[0] );

    // rupture physics send- and receive-faces
    l_shared.regWrkRgn( l_tg,
                        2,
                        l_tg * N_ENTRIES_CONTROL_FLOW + 7,
                        l_enLayouts[l_rupLayoutFa].timeGroups[l_tg].inner.first+
                        l_enLayouts[l_rupLayoutFa].timeGroups[l_tg].inner.size,
                        l_enLayouts[l_rupLayoutFa].timeGroups[l_tg].nEntsOwn-
                        l_enLayouts[l_rupLayoutFa].timeGroups[l_tg].inner.size+
                        l_enLayouts[l_rupLayoutFa].timeGroups[l_tg].nEntsNotOwn,
                        l_enLayouts[l_rupLayoutFa].timeGroups.size() + l_tg,
                        1, l_spType, l_internal.m_faceSparseShared6[0] );
  }
}
