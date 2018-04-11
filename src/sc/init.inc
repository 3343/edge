/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2017-2018, Regents of the University of California
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
 * Performs the init of sub-cell data, once limited, limited plus and and extrema sparse types are derived.
 **/
#if defined(PP_T_KERNELS_XSMM)
// setup gather/scatter kernels
{
  unsigned int l_n_qts = N_QUANTITIES;
  unsigned int l_n_mds = N_ELEMENT_MODES;
  unsigned int l_n_scs = CE_N_SUB_CELLS(T_SDISC.ELEMENT,ORDER); 
  unsigned int l_n_sfs = CE_N_SUB_FACES(T_SDISC.ELEMENT,ORDER);

  // scatter
  l_internal.m_mm.addScAC( l_n_qts, l_n_scs, l_n_mds, // m, n, k
                           l_n_mds, l_n_scs, l_n_scs, // ldA, ldB, lcC
                           static_cast<real_base>(1.0), static_cast<real_base>(0.0) /*, LIBXSMM_GEMM_PREFETCH_NONE*/ );

  // scatter fa
  l_internal.m_mm.addScAC( l_n_qts, l_n_sfs, l_n_mds, // m, n, k
                           l_n_mds, l_n_sfs, l_n_sfs, // ldA, ldB, lcC
                           static_cast<real_base>(1.0), static_cast<real_base>(0.0) /*, LIBXSMM_GEMM_PREFETCH_NONE*/ );

  // gather
  l_internal.m_mm.addScAC( l_n_qts, l_n_mds, l_n_scs, // m, n, k
                           l_n_scs, l_n_mds, l_n_mds, // ldA, ldB, lcC
                           static_cast<real_base>(1.0), static_cast<real_base>(0.0) /*, LIBXSMM_GEMM_PREFETCH_NONE*/ );

  // sfInt (gather fa)
  l_internal.m_mm.addScAC( l_n_qts, l_n_mds, l_n_sfs, // m, n, k
                           l_n_sfs, l_n_mds, l_n_mds, // ldA, ldB, lcC
                           static_cast<real_base>(1.0), static_cast<real_base>(0.0) /*, LIBXSMM_GEMM_PREFETCH_NONE*/ );

  // internal boundary transform
  l_internal.m_mm.addScBC( l_n_qts, l_n_sfs, l_n_qts,
                           l_n_qts, l_n_sfs, l_n_sfs,
                           static_cast<real_base>(1.0), static_cast<real_base>(0.0) /*, LIBXSMM_GEMM_PREFETCH_NONE*/ ); 
}
#endif
unsigned short l_limLayout;
unsigned short l_limPlusLayout;
unsigned short l_extLayout ;
{
  // derive an entity layout for the limited elements
  l_enLayouts.resize( l_enLayouts.size() + 1 );
  edge::data::SparseEntities::denseToSparse( t_enTypeShared::LIMIT,
                                             l_internal.m_elementChars,
                                             l_enLayouts[2],
                                             l_enLayouts.back() );
  l_limLayout = l_enLayouts.size()-1;

  // derive an entity layout for the limited plus elements
  l_enLayouts.resize( l_enLayouts.size() + 1 );
  edge::data::SparseEntities::denseToSparse( t_enTypeShared::LIMIT_PLUS,
                                             l_internal.m_elementChars,
                                             l_enLayouts[2],
                                             l_enLayouts.back() );
  l_limPlusLayout = l_enLayouts.size()-1;

  // derive an layout for extrema
  l_enLayouts.resize( l_enLayouts.size() + 1 );
  edge::data::SparseEntities::denseToSparse( t_enTypeShared::EXTREMA,
                                             l_internal.m_elementChars,
                                             l_enLayouts[2],
                                             l_enLayouts.back() );
  l_extLayout = l_enLayouts.size()-1;

  // get sub-cell specific entity counts and perform sanity checks
  int_el l_nLi = l_enLayouts[l_limLayout    ].nEnts;
  int_el l_nLp = l_enLayouts[l_limPlusLayout].nEnts;
  int_el l_nEx = l_enLayouts[l_extLayout    ].nEnts;

  EDGE_CHECK_LE( l_nLi, l_nLp );
  EDGE_CHECK_LE( l_nLp, l_nEx );

  // allocate sub-cell memory
  edge::sc::Memory::alloc( l_nLi,
                           l_nLp,
                           l_nEx,
                           l_dynMem,
                           l_internal.m_globalShared2[0] );

  // init sub-cell memory
  edge::sc::Init< T_SDISC.ELEMENT,
                  ORDER,
                  N_QUANTITIES,
                  N_CRUNS >::data( l_nLi,
                                   l_enLayouts[l_limPlusLayout].nEnts,
                                   l_nEx,
                                   l_internal.m_globalShared2[0].dofs,
                                   l_internal.m_globalShared2[0].tDofsRaw,
                                   l_internal.m_globalShared2[0].ext,
                                   l_internal.m_globalShared2[0].adm,
                                   l_internal.m_globalShared2[0].lock,
                                   l_internal.m_globalShared2[0].limSync );

  // link limited elements to limited plus elements
  edge::data::SparseEntities::linkSpSp( l_enLayouts[2].nEnts,
                                        t_enTypeShared::LIMIT,
                                        t_enTypeShared::LIMIT_PLUS,
                                        l_internal.m_elementChars,
                                        l_internal.m_globalShared2[0].connect.liLp );

  // link limited plus elements to dense elements (no bridge)
  edge::data::SparseEntities::linkSpDe( l_enLayouts[2].nEnts,
                                        t_enTypeShared::LIMIT_PLUS,
                                        l_internal.m_elementChars,
                                        l_internal.m_globalShared2[0].connect.lpEl );

  // link limited plus elements to limited elements (no bridge)
  edge::data::SparseEntities::linkSpSp( l_enLayouts[2].nEnts,
                                        t_enTypeShared::LIMIT_PLUS,
                                        t_enTypeShared::LIMIT,
                                        l_internal.m_elementChars,
                                        l_internal.m_globalShared2[0].connect.lpLi );

  // link limited plus elements to limited plus elements (faces as bridge)
  edge::data::SparseEntities::linkSpAdjDst( l_enLayouts[2].nEnts,
                                            C_ENT[T_SDISC.ELEMENT].N_FACES,
                                            l_internal.m_connect.elFaEl[0],
                                            t_enTypeShared::LIMIT_PLUS,
                                            t_enTypeShared::LIMIT_PLUS,
                                            l_internal.m_elementChars,
                                            l_internal.m_elementChars,
                                            l_internal.m_globalShared2[0].connect.lpFaLp[0] );

  // link limited elements to extrema (no bridge)
  edge::data::SparseEntities::linkSpSp( l_enLayouts[2].nEnts,
                                        t_enTypeShared::LIMIT,
                                        t_enTypeShared::EXTREMA,
                                        l_internal.m_elementChars,
                                        l_internal.m_globalShared2[0].connect.liEx );

  // setup pointer data structure for sub-cell tDofs (dual buffer)
  std::vector< std::vector< unsigned char * > > l_tDofsScMpi[2][2];
  for( unsigned short l_bu = 0; l_bu < 2; l_bu++ )
    edge::data::DataLayout::adj( C_ENT[T_SDISC.ELEMENT].N_FACES,
                                 l_internal.m_globalShared2[0].connect.lpFaLp[0],
                                 l_enLayouts[l_limPlusLayout],
                                 l_internal.m_globalShared2[0].tDofsRaw[l_bu],
                                 l_internal.m_globalShared2[0].tDofs[l_bu],
                                 l_tDofsScMpi[l_bu][0],
                                 l_tDofsScMpi[l_bu][1],
                                 l_internal.m_globalShared2[0].connect.lpLi );

  // assemble dominant, limited
  edge::sc::Init< T_SDISC.ELEMENT,
                  ORDER,
                  N_QUANTITIES,
                  N_CRUNS >::liDoLiDu(  l_enLayouts[l_limLayout],
                                        l_internal.m_globalShared2[0].connect.liLp,
                                        l_internal.m_globalShared2[0].connect.lpEl,
                                       &l_mesh.getInMap()->elDaMe[0],
                                       &l_mesh.getInMap()->elMeDa[0],
                                        l_internal.m_globalShared2[0].connect.liDoLiDu  );

#ifdef PP_USE_MPI
  // register MPI layouts
  EDGE_CHECK_EQ( l_enLayouts[l_limPlusLayout].timeGroups.size(), 1 );

  // candidate admissibility in even time steps since sync
  l_mpi.addDefault( l_enLayouts[l_limLayout],
                    l_internal.m_globalShared2[0].adm[1],
                    N_CRUNS*sizeof(bool),
                    0,
                    1,
                    reinterpret_cast< std::intptr_t >( l_internal.m_globalShared2[0].adm[1] ) );

  // candidate admissibility in odd time steps since sync
  l_mpi.addDefault( l_enLayouts[l_limLayout],
                    l_internal.m_globalShared2[0].adm[3],
                    N_CRUNS*sizeof(bool),
                    0,
                    1,
                    reinterpret_cast< std::intptr_t >( l_internal.m_globalShared2[0].adm[3] ) );

  for( unsigned short l_bu = 0; l_bu < 2; l_bu++ ) {
    l_mpi.addCustom( l_enLayouts[l_limPlusLayout],
                     l_tDofsScMpi[l_bu][0],
                     l_tDofsScMpi[l_bu][1],
                     0,
                     1,
                     reinterpret_cast< std::intptr_t >( l_internal.m_globalShared2[0].tDofs[l_bu] ) );

    l_mpi.addDefault( l_enLayouts[l_extLayout],
                      l_internal.m_globalShared2[0].ext[l_bu],
                      2*N_QUANTITIES*N_CRUNS*sizeof(real_base),
                      0,
                      1,
                      reinterpret_cast< std::intptr_t >( l_internal.m_globalShared2[0].ext[l_bu] ) );
  }
#endif
}
