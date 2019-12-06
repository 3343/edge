/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *         Alexander Heinecke (alexander.heinecke AT intel.com)
 *
 * @section LICENSE
 * Copyright (c) 2019, Alexander Breuer
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
 * Setup for seismic simulations.
 **/

// set up sparse types
{
  // propagate sparse info of free-surface boundary to elements
  edge::data::SparseEntities::propAdj( l_edgeV.nFas(),
                                       2,
                                       l_internal.m_connect.faEl[0],
                                       t_spTypeElastic::FREE_SURFACE,
                                       t_spTypeElastic::FREE_SURFACE,
                                       l_internal.m_faceChars,
                                       l_internal.m_elementChars );
}

// allocate flex data for time buffers and DOFs
{
  int_spType  l_spTypes[2] = { C_LTS_EL[t_ltsEl::EL_INT_LT], C_LTS_EL[t_ltsEl::EL_INT_GT] };
  std::size_t l_spSizes[2] = { N_QUANTITIES,                 N_QUANTITIES                 };

  real_base (**l_raw)[N_ELEMENT_MODES][N_CRUNS];

  l_raw = l_dynMem.flex<
    real_base [N_ELEMENT_MODES][N_CRUNS] >( l_edgeV.nEls(),
                                            2,
                                            N_QUANTITIES,
                                            l_spTypes,
                                            l_spSizes,
                                            l_internal.m_elementChars,
                                            ALIGNMENT.BASE.HEAP,
                                            true,
                                            true );

  l_internal.m_globalShared6[0] = l_raw;
  l_internal.m_globalShared6[1] = l_raw +   l_edgeV.nEls();
  l_internal.m_globalShared6[2] = l_raw + 2*l_edgeV.nEls();
}

l_mpiRemix.init( l_edgeV.nTgs(),
                 C_ENT[T_SDISC.ELEMENT].N_FACES,
                 l_edgeV.nEls(),
                 N_QUANTITIES*N_ELEMENT_MODES*N_CRUNS*sizeof(real_base),
                 l_edgeV.getCommStruct(),
                 l_edgeV.getSendFa(),
                 l_edgeV.getSendEl(),
                 l_edgeV.getRecvFa(),
                 l_edgeV.getRecvEl(),
                 l_dynMem );

// parse config specific to elastics
edge::seismic::io::Config l_seismicConf( l_config.m_doc );

// setup initial DOFs
EDGE_LOG_INFO << "  setting up material parameters and initial DOFs";

// get the initial setup
PP_INSTR_REG_DEF(dofsMat)
PP_INSTR_REG_BEG(dofsMat,"dofs_mat")

// query velocity model from mesh
unsigned short l_vmMesh = std::numeric_limits< unsigned short >::max();
EDGE_LOG_INFO << "    querying mesh for velocity model";

// TODO: add velocity model queries to EDGE-V interface

if( l_vmMesh == 0 ) {
  EDGE_LOG_INFO << "    successfully obtained velocity model from mesh";
}
else {
  EDGE_LOG_INFO << "    failed obtaining velocity model from mesh, continuining w/o "
                << l_vmMesh;
}

// perform NUMA-aware zero-init of DOFs and tDOFs
#ifdef PP_USE_OMP
#pragma omp parallel
#endif
{
  // elastic DOFs
  edge::parallel::Shared::numaInit( l_edgeV.nTgs(),
                                    l_edgeV.nTgElsIn(),
                                    l_edgeV.nTgElsSe(),
                                    l_shared.m_nWrks,
                                    N_QUANTITIES*N_ELEMENT_MODES*N_CRUNS,
                                    l_internal.m_elementModePrivate1[0][0][0] );
  // anelastic DOFs
#if (PP_N_RELAXATION_MECHANISMS > 0)
  edge::parallel::Shared::numaInit( l_edgeV.nTgs(),
                                    l_edgeV.nTgElsIn(),
                                    l_edgeV.nTgElsSe(),
                                    l_shared.m_nWrks,
                                    N_QUANTITIES_A*N_ELEMENT_MODES*N_CRUNS,
                                    l_internal.m_elementModePrivate2[0][0][0] );
#endif

  // tDOFs
  edge::parallel::Shared::numaInit( l_edgeV.nTgs(),
                                    l_edgeV.nTgElsIn(),
                                    l_edgeV.nTgElsSe(),
                                    l_shared.m_nWrks,
                                    N_QUANTITIES*N_ELEMENT_MODES*N_CRUNS,
                                    l_internal.m_globalShared6[0][0][0][0] );
}

// init LTS DOFs
#ifdef PP_USE_OMP
#pragma omp parallel for
#endif
for( int_el l_el = 0; l_el < l_edgeV.nEls(); l_el++ ) {
  for( int_qt l_qt = 0; l_qt < N_QUANTITIES; l_qt++ ) {
    for( int_md l_md = 0; l_md < N_ELEMENT_MODES; l_md++ ) {
      for( int_cfr l_ru = 0; l_ru < N_CRUNS; l_ru++ ) {
        if( l_internal.m_globalShared6[1][l_el] != nullptr ) {
          l_internal.m_globalShared6[1][l_el][l_qt][l_md][l_ru] = 0;
        }
        if( l_internal.m_globalShared6[2][l_el] != nullptr ) {
          l_internal.m_globalShared6[2][l_el][l_qt][l_md][l_ru] = 0;
        }
      }
    }
  }
}

EDGE_LOG_INFO << "    setting initial DOFs and velocity model based on user-provided config (if available)";
// initialize DOFs
{
  // DG
  edge::setups::InitialDofs<
    T_SDISC.ELEMENT,
    ORDER,
    N_QUANTITIES,
    N_CRUNS >::dg( std::size_t(0),
                   l_edgeV.nEls(),
                   l_config.m_initValsExprStrs,
                   l_basis,
                   l_internal.m_connect.elVe,
                   l_internal.m_vertexChars,
                   l_internal.m_elementModePrivate1 );
}

// set up material parameters
{
  for( int_el l_el = 0; l_el < l_edgeV.nEls(); l_el++ ) {
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
    EDGE_CHECK( (l_vmMesh == 0) || (l_seismicConf.m_velDoms.size() > 0) )
      << "couldn't find a velocity model in the mesh or user config, aborting";

    // find the matching domain in the velocity model
    for( std::size_t l_do = 0; l_do < l_seismicConf.m_velDoms.size(); l_do++ ) {
      if( l_seismicConf.m_velDoms[l_do].inside(l_ave) ) {
        l_internal.m_elementShared1[l_el][0].rho = l_seismicConf.m_velVals[l_do][0];
        l_internal.m_elementShared1[l_el][0].lam = l_seismicConf.m_velVals[l_do][1];
        l_internal.m_elementShared1[l_el][0].mu  = l_seismicConf.m_velVals[l_do][2];
        if( N_RELAXATION_MECHANISMS > 0 ) {
          l_internal.m_elementShared1[l_el][0].qp = l_seismicConf.m_velVals[l_do][3];
          l_internal.m_elementShared1[l_el][0].qs = l_seismicConf.m_velVals[l_do][4];
        }
        break;
      }

      // abort if no matching velocity domain is present
      EDGE_CHECK( (l_vmMesh == 0) ||
                  (l_do != l_seismicConf.m_velDoms.size()-1 ) )
        << "here is the troublesome point: "
        << l_ave[0] << " " << l_ave[1] << " " << ( (N_DIM > 2) ? std::to_string(l_ave[2]) : "" );
    }
  }
}
PP_INSTR_REG_END(dofsMat)

// initialize ADER-DG solver and determine elastic material parameters
edge::seismic::solvers::AderDg<
 real_base,
 N_RELAXATION_MECHANISMS,
 T_SDISC.ELEMENT,
 ORDER,
 ORDER,
 N_CRUNS,
 MM_KERNELS_SPARSE > l_aderDg( l_edgeV.nEls(),
                               l_edgeV.nFas(),
                               l_internal.m_connect.faEl,
                               l_internal.m_connect.elVe,
                               l_internal.m_connect.elFa,
                               l_internal.m_vertexChars,
                               l_internal.m_faceChars,
                               l_internal.m_elementChars,
                               (t_bgPars*) l_internal.m_elementShared1,
                               l_seismicConf.m_attFreqs[0],
                               l_seismicConf.m_attFreqs[1],
                               l_dynMem );
l_internal.m_globalShared4[0] = &l_aderDg;

// setup point sources
if( l_seismicConf.m_ptSrcs.size() > 0 ) {
  PP_INSTR_REG_DEF(ptsrcs)
  PP_INSTR_REG_BEG(ptsrcs,"pt_srcs")
#ifdef PP_HAS_HDF5
  EDGE_LOG_INFO << "  initializing point sources";
  bool l_psInit = l_internal.m_globalShared3[0].init( l_seismicConf.m_ptSrcs.data(),
                                                      t_spTypeElastic::SOURCE,
                                                      l_edgeV.nEls(),
                                                      l_internal.m_connect.elVe,
                                                      l_internal.m_vertexChars,
                                                      l_internal.m_globalShared1[0].mat.massI,
                                                      l_internal.m_elementChars,
                                                      l_dynMem );
  if( l_psInit ) { EDGE_LOG_INFO << "  point sources have been initialized";                   }
  else {           EDGE_LOG_INFO << "  failed initializing the point sources, continuing w/o"; }
#else
  EDGE_LOG_INFO << "  no HDF5 support, skipping point source initialization";
#endif
  PP_INSTR_REG_END(ptsrcs)
}

// get layout of sparse source elements
std::size_t *l_nTgElsInSrc = new std::size_t[ l_edgeV.nTgs() ];
std::size_t *l_nTgElsSeSrc = new std::size_t[ l_edgeV.nTgs() ];
edge::data::SparseEntities::denseToSparse( l_edgeV.nTgs(),
                                           t_spTypeElastic::SOURCE,
                                           l_internal.m_elementChars,
                                           l_edgeV.nTgElsIn(),
                                           l_edgeV.nTgElsSe(),
                                           l_nTgElsInSrc,
                                           l_nTgElsSeSrc );

// get time step stats
edge::seismic::common::getTimeStepStatsCFL( l_edgeV.nEls(),
                                            l_internal.m_elementChars,
                                            l_internal.m_elementShared1,
                                            l_dT[0], l_dT[1], l_dT[2] );

// setup shared memory parallelization
std::size_t l_firstDe = 0;
std::size_t l_firstSrc = 0;

for( std::size_t l_tg = 0; l_tg < l_edgeV.nTgs(); l_tg++ ) {
  int_spType l_spType[1] = { RECEIVER };

  // local, inner
  l_shared.regWrkRgn( l_tg,
                      0,
                      l_tg * N_ENTRIES_CONTROL_FLOW + 0,
                      l_firstDe,
                      l_edgeV.nTgElsIn()[l_tg],
                      l_edgeV.nTgs() - l_tg,
                      1, l_spType, l_internal.m_elementChars );

  // src, inner
  l_shared.regWrkRgn( l_tg,
                      1,
                      l_tg * N_ENTRIES_CONTROL_FLOW + 2,
                      l_firstSrc,
                      l_nTgElsInSrc[l_tg],
                      l_edgeV.nTgs() - l_tg );

  // neigh, inner
  l_shared.regWrkRgn( l_tg,
                      2,
                      l_tg * N_ENTRIES_CONTROL_FLOW + 4,
                      l_firstDe,
                      l_edgeV.nTgElsIn()[l_tg],
                      l_edgeV.nTgs() - l_tg );

  l_firstDe  += l_edgeV.nTgElsIn()[l_tg];
  l_firstSrc += l_nTgElsInSrc[l_tg];
}

for( std::size_t l_tg = 0; l_tg < l_edgeV.nTgs(); l_tg++ ) {
  int_spType l_spType[1] = { RECEIVER };

  // local, send
  l_shared.regWrkRgn( l_tg,
                      0,
                      l_tg * N_ENTRIES_CONTROL_FLOW + 1,
                      l_firstDe,
                      l_edgeV.nTgElsSe()[l_tg],
                      l_edgeV.nTgs()*2 - l_tg,
                      1, l_spType, l_internal.m_elementChars );

  // src, send
  l_shared.regWrkRgn( l_tg,
                      1,
                      l_tg * N_ENTRIES_CONTROL_FLOW + 3,
                      l_firstSrc,
                      l_nTgElsSeSrc[l_tg],
                      l_edgeV.nTgs()*2 - l_tg );

  // neigh, send
  l_shared.regWrkRgn( l_tg,
                      2,
                      l_tg * N_ENTRIES_CONTROL_FLOW + 5,
                      l_firstDe,
                      l_edgeV.nTgElsSe()[l_tg],
                      l_edgeV.nTgs()*2 - l_tg );

  l_firstDe  += l_edgeV.nTgElsSe()[l_tg];
  l_firstSrc += l_nTgElsSeSrc[l_tg];
}

delete[] l_nTgElsInSrc;
delete[] l_nTgElsSeSrc;