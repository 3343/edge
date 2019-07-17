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
 * Setup for elastics.
 **/

// set up sparse types
{
  // propagate sparse info of rupture faces to elements
  edge::data::SparseEntities::propAdj(  l_enLayouts[1].nEnts,
                                        2,
                                        l_internal.m_connect.faEl[0],
                                        t_spTypeElastic::RUPTURE,
                                        t_spTypeElastic::RUPTURE,
                                        l_internal.m_faceChars,
                                        l_internal.m_elementChars,
                                       &l_mesh.getInMap()->elDaMe[0],
                                       &l_mesh.getInMap()->elMeDa[0] );

#ifdef PP_USE_MPI
  edge::parallel::Mpi::syncSpTypes( l_enLayouts[2],
                                    l_internal.m_elementChars );
#endif

  // inherit limited from dynamic rupture elements
  edge::data::SparseEntities::inherit(              l_internal.m_nElements,
                                       (int_spType) t_spTypeElastic::RUPTURE,
                                       (int_spType) t_enTypeShared::LIMIT,
                                                    l_internal.m_elementChars );

  // propagate sparse info of free-surface boundary to elements
  edge::data::SparseEntities::propAdj( l_enLayouts[1].nEnts,
                                       2,
                                       l_internal.m_connect.faEl[0],
                                       t_spTypeElastic::FREE_SURFACE,
                                       t_spTypeElastic::FREE_SURFACE,
                                       l_internal.m_faceChars,
                                       l_internal.m_elementChars,
                                       &l_mesh.getInMap()->elDaMe[0],
                                       &l_mesh.getInMap()->elMeDa[0] );

#ifdef PP_USE_MPI
  edge::parallel::Mpi::syncSpTypes( l_enLayouts[2],
                                    l_internal.m_elementChars );
#endif

  // perform sweeps over face-neighbors to extend limited region
  for( unsigned short l_sw = 0; l_sw < 3; l_sw++ ) {
    // inherit limit plus from limited elements (superset)
    edge::data::SparseEntities::inherit(              l_internal.m_nElements,
                                         (int_spType) t_enTypeShared::LIMIT,
                                         (int_spType) t_enTypeShared::LIMIT_PLUS,
                                                      l_internal.m_elementChars );

    // propagate limit through faces
    edge::data::SparseEntities::propAdj(               l_internal.m_nElements,
                                                       C_ENT[T_SDISC.ELEMENT].N_FACES,
                                                       l_internal.m_connect.elFaEl[0],
                                         (int_spType)  t_enTypeShared::LIMIT_PLUS,
                                         (int_spType)  t_enTypeShared::LIMIT,
                                                       l_internal.m_elementChars,
                                                       l_internal.m_elementChars,
                                                      &l_mesh.getInMap()->elDaMe[0],
                                                      &l_mesh.getInMap()->elMeDa[0] );

#ifdef PP_USE_MPI
    edge::parallel::Mpi::syncSpTypes( l_enLayouts[2],
                                      l_internal.m_elementChars );
#endif
  }

  // limit plus are the limited elements and their face neighbors
  edge::data::SparseEntities::propAdj(               l_internal.m_nElements,
                                                     C_ENT[T_SDISC.ELEMENT].N_FACES,
                                                     l_internal.m_connect.elFaEl[0],
                                       (int_spType)  t_enTypeShared::LIMIT,
                                       (int_spType)  t_enTypeShared::LIMIT_PLUS,
                                                     l_internal.m_elementChars,
                                                     l_internal.m_elementChars,
                                                    &l_mesh.getInMap()->elDaMe[0],
                                                    &l_mesh.getInMap()->elMeDa[0] );
#ifdef PP_USE_MPI
  edge::parallel::Mpi::syncSpTypes( l_enLayouts[2],
                                    l_internal.m_elementChars );
#endif

  // store DG DOFS for all limit plus elements
  edge::data::SparseEntities::inherit(              l_internal.m_nElements,
                                       (int_spType) t_enTypeShared::LIMIT_PLUS,
                                       (int_spType) C_LTS_EL[t_ltsEl::EL_DOFS],
                                                    l_internal.m_elementChars );

  // store extrema for all limit plus elements
  edge::data::SparseEntities::inherit(              l_internal.m_nElements,
                                       (int_spType) t_enTypeShared::LIMIT_PLUS,
                                       (int_spType) t_enTypeShared::EXTREMA,
                                                    l_internal.m_elementChars );
}

// allocate flex data for time buffers and DOFs
{
  int_spType  l_spTypes[1] = { C_LTS_EL[t_ltsEl::EL_DOFS] };
  std::size_t l_spSizes[1] = { N_QUANTITIES };

  real_base (**l_raw)[N_ELEMENT_MODES][N_CRUNS];

  l_raw = l_dynMem.flex<
    real_base [N_ELEMENT_MODES][N_CRUNS] >( l_internal.m_nElements,
                                            1,
                                            N_QUANTITIES,
                                            l_spTypes,
                                            l_spSizes,
                                            l_internal.m_elementChars,
                                            ALIGNMENT.BASE.HEAP,
                                            true,
                                            true );

  l_internal.m_globalShared6[0] = l_raw;
  l_internal.m_globalShared6[1] = l_raw + l_internal.m_nElements;
}

#ifdef PP_USE_MPI
  // init mpi layout
  EDGE_CHECK( l_enLayouts[2].timeGroups.size() == 1 );

  l_mpi.addDefault( l_enLayouts[2],
                    l_internal.m_globalShared6[0][0][0][0],
                    N_QUANTITIES*N_ELEMENT_MODES*N_CRUNS*sizeof(real_base),
                    0,
                    1 );
#endif

// init data of limiter
#include "sc/init.inc"

// parse config specific to elastics
edge::elastic::io::Config l_elasticConf( l_config.m_doc );

// setup initial DOFs
EDGE_LOG_INFO << "  setting up material parameters and initial DOFs";

// get the initial setup
PP_INSTR_REG_DEF(dofsMat)
PP_INSTR_REG_BEG(dofsMat,"dofs_mat")

// query velocity model from mesh
unsigned short l_vmMesh = std::numeric_limits< unsigned short >::max();
EDGE_LOG_INFO << "    querying mesh for velocity model";

#if defined PP_T_MESH_UNSTRUCTURED
// allocate memory for velocity model
double (*l_velMod)[5] = (double (*)[5]) new double[ l_internal.m_nElements * 5 ];
std::string l_bgPars[5] = { "LAMBDA",
                            "MU",
                            "RHO",
                            "QP",
                            "QS" };

l_vmMesh = l_mesh.getParsDe( N_DIM,
                             (N_RELAXATION_MECHANISMS == 0) ? 3 : 5,
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
    // check for valid q-factors
    if( N_RELAXATION_MECHANISMS > 0 ) {
      EDGE_CHECK_GT( l_velMod[l_el][3], 0 ) << l_el;
      EDGE_CHECK_GT( l_velMod[l_el][4], 0 ) << l_el;
    }

    l_internal.m_elementShared1[l_el][0].lam  = l_velMod[l_el][0];
    l_internal.m_elementShared1[l_el][0].mu   = l_velMod[l_el][1];
    l_internal.m_elementShared1[l_el][0].rho  = l_velMod[l_el][2];
    if( N_RELAXATION_MECHANISMS > 0 ) {
      l_internal.m_elementShared1[l_el][0].qp = l_velMod[l_el][3];
      l_internal.m_elementShared1[l_el][0].qs = l_velMod[l_el][4];
    }
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
#pragma omp parallel
#endif
{
  // star matrix
#if defined PP_T_KERNELS_XSMM
  edge::parallel::Shared::numaInit( l_enLayouts[2],
                                    l_shared.m_nWrks,
                                    N_MAT_STAR,
                                    l_internal.m_elementShared4[0][0].mat );
#else
  edge::parallel::Shared::numaInit( l_enLayouts[2],
                                    l_shared.m_nWrks,
                                    N_QUANTITIES*N_QUANTITIES,
                                    l_internal.m_elementShared4[0][0].mat[0] );
#endif

  // flux solver own
  edge::parallel::Shared::numaInit( l_enLayouts[2],
                                    l_shared.m_nWrks,
                                    N_QUANTITIES*N_QUANTITIES,
                                    l_internal.m_elementShared2[0][0].solver[0] );
  // flux solver neigh
  edge::parallel::Shared::numaInit( l_enLayouts[2],
                                    l_shared.m_nWrks,
                                    N_QUANTITIES*N_QUANTITIES,
                                    l_internal.m_elementShared3[0][0].solver[0] );

  // elastic DOFs
  edge::parallel::Shared::numaInit( l_enLayouts[2],
                                    l_shared.m_nWrks,
                                    N_QUANTITIES*N_ELEMENT_MODES*N_CRUNS,
                                    l_internal.m_elementModePrivate1[0][0][0] );
  // anelastic DOFs
#if (PP_N_RELAXATION_MECHANISMS > 0)
  edge::parallel::Shared::numaInit( l_enLayouts[2],
                                    l_shared.m_nWrks,
                                    N_QUANTITIES_A*N_ELEMENT_MODES*N_CRUNS,
                                    l_internal.m_elementModePrivate2[0][0][0] );
#endif

  // tDOFs
  edge::parallel::Shared::numaInit( l_enLayouts[2],
                                    l_shared.m_nWrks,
                                    N_QUANTITIES*N_ELEMENT_MODES*N_CRUNS,
                                    l_internal.m_globalShared6[0][0][0][0] );
}

// init buffer DOFs
// TODO: Add NUMA init for limiter and LTS
#ifdef PP_USE_OMP
#pragma omp parallel for
#endif
for( int_el l_el = 0; l_el < l_internal.m_nElements; l_el++ ) {
  for( int_qt l_qt = 0; l_qt < N_QUANTITIES; l_qt++ ) {
    for( int_md l_md = 0; l_md < N_ELEMENT_MODES; l_md++ ) {
      for( int_cfr l_ru = 0; l_ru < N_CRUNS; l_ru++ ) {
        if( l_internal.m_globalShared6[1][l_el] != nullptr )
          l_internal.m_globalShared6[1][l_el][l_qt][l_md][l_ru] = 0;
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
    N_CRUNS >::dg( 0,
                   l_internal.m_nElements,
                   l_config.m_initValsExprStrs,
                   l_basis,
                   l_internal.m_connect.elVe,
                   l_internal.m_vertexChars,
                   l_internal.m_elementModePrivate1 );
}

// set up material parameters
{
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
        if( N_RELAXATION_MECHANISMS > 0 ) {
          l_internal.m_elementShared1[l_el][0].qp = l_elasticConf.m_velVals[l_do][3];
          l_internal.m_elementShared1[l_el][0].qs = l_elasticConf.m_velVals[l_do][4];
        }
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
PP_INSTR_REG_END(dofsMat)

// initialize ADER-DG solver and determine elastic material parameters
edge::seismic::solvers::AderDg<
 real_base,
 N_RELAXATION_MECHANISMS,
 T_SDISC.ELEMENT,
 ORDER,
 ORDER,
 N_CRUNS,
 MM_KERNELS_SPARSE > l_aderDg( l_internal.m_nElements,
                               l_internal.m_nFaces,
                               l_internal.m_connect.faEl,
                               l_internal.m_connect.elVe,
                               l_internal.m_connect.elFa,
                               l_mesh.getInMap()->elMeDa.data(),
                               l_mesh.getInMap()->elDaMe.data(),
                               l_internal.m_vertexChars,
                               l_internal.m_faceChars,
                               l_internal.m_elementChars,
                               (t_bgPars*) l_internal.m_elementShared1,
                               l_elasticConf.m_attFreqs[0],
                               l_elasticConf.m_attFreqs[1],
                               l_dynMem );
l_internal.m_globalShared4[0] = &l_aderDg;

// setup point sources
if( l_elasticConf.m_ptSrcs.size() > 0 ) {
  PP_INSTR_REG_DEF(ptsrcs)
  PP_INSTR_REG_BEG(ptsrcs,"pt_srcs")
#ifdef PP_HAS_HDF5
  EDGE_LOG_INFO << "  initializing point sources";
  bool l_psInit = l_internal.m_globalShared3[0].init( l_elasticConf.m_ptSrcs.data(),
                                                      &l_gIdsEl[0],
                                                      t_spTypeElastic::SOURCE,
                                                      l_enLayouts[2],
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

  l_recvsSf.init(                     l_config.m_recvCrds[1].size(),
                                      t_spTypeElastic::RUPTURE,
                                      (N_DIM-1)*3, // TODO: hardcoded to linear slip weakening
                                      l_config.m_recvPath[1],
                                      &l_config.m_recvNames[1][0],
              (real_mesh (*)[N_DIM]) &l_crds[0][0], // TODO: dimension
                                      l_config.m_recvFreq[1],
                                      l_enLayouts[1],
                                      l_enLayouts[2],
                                      l_internal.m_globalShared2[0].connect.scSv,
                                      l_internal.m_connect.faEl,
                                      l_internal.m_connect.elVe,
                                      l_internal.m_connect.elFa,
                                      l_internal.m_globalShared2[0].svChars,
                                      l_internal.m_vertexChars,
                                      l_internal.m_faceChars );
}

// init sparse, internal data structures
l_internal.initSparse( 0, l_enLayouts[l_rupLayoutFa].nEnts, 0,
                       0, l_enLayouts[l_rupLayoutFa].nEnts, l_enLayouts[l_rupLayoutEl].nEnts,
                       0, l_enLayouts[l_rupLayoutFa].nEnts, l_enLayouts[l_rupLayoutEl].nEnts,
                       0, l_enLayouts[l_rupLayoutFa].nEnts, 0,
                       0, l_enLayouts[l_rupLayoutFa].nEnts, 0,
                       0, l_enLayouts[l_rupLayoutFa].nEnts, 0  );

if( l_elasticConf.m_frictionLaw != "" ) {
  // allocate memory for the internal boundary
  edge::sc::ibnd::Memory::alloc( l_enLayouts[l_rupLayoutFa].nEnts,
                                 l_enLayouts[l_rupLayoutEl].nEnts,
                                 l_dynMem,
                                 l_internal.m_globalShared7[0] );

  // init connectivity of the internal boundary
  edge::sc::ibnd::Init< T_SDISC.ELEMENT >::connect( l_enLayouts[1].nEnts,
                                                    l_enLayouts[2].nEnts,
                                       (int_spType) t_spTypeElastic::RUPTURE,
                                       (int_spType) t_enTypeShared::LIMIT,
                                                    l_internal.m_connect.faEl,
                                                    l_internal.m_connect.elFa,
                                                    l_internal.m_faceChars,
                                                    l_internal.m_elementChars,
                                                    l_internal.m_globalShared7[0].connect );

  // set up ibnd face chars
  edge::elastic::solvers::InternalBoundaryTypes<
    T_SDISC.ELEMENT
  >::initFaces( l_enLayouts[1].nEnts,
                l_enLayouts[2].nEnts,
                RUPTURE,
                l_internal.m_faceChars,
                l_internal.m_elementChars,
                l_internal.m_connect.faEl,
                l_internal.m_connect.elFa,
                l_internal.m_connect.vIdElFaEl,
                l_internal.m_globalShared7[0].bfChars );

  // set up rupture solvers
  edge::elastic::solvers::InternalBoundarySolvers<
    T_SDISC.ELEMENT,
    ORDER
  >::init( l_enLayouts[1].nEnts,
           t_spTypeElastic::RUPTURE,
           l_internal.m_connect.faEl,
           l_internal.m_connect.elVe,
           l_internal.m_faceChars,
           l_internal.m_elementChars,
           l_internal.m_elementShared1[0],
           l_internal.m_globalShared7[0].mss,
           l_elasticConf.m_faultCrds );

  // query mesh for fault data
  unsigned short l_faultMesh = std::numeric_limits< unsigned short >::max();
  EDGE_LOG_INFO << "    querying mesh for fault data";

  struct t_FaultData{ double sn0[CE_N_SUB_FACES(T_SDISC.ELEMENT, ORDER)][N_CRUNS];
                      double co0[CE_N_SUB_FACES(T_SDISC.ELEMENT, ORDER)][N_CRUNS];
                      double ss0[N_DIM-1][CE_N_SUB_FACES(T_SDISC.ELEMENT, ORDER)][N_CRUNS];
                    } *l_faultData;

  l_faultData = new t_FaultData[ l_enLayouts[l_rupLayoutFa].nEnts ];
  std::string l_faultPars[4] = { "N_STRESS",
                                 "COHESION",
                                 "S_STRESS_STK",
                                 "S_STRESS_DIP" };

  l_faultMesh = l_mesh.getParsDe(                                   N_DIM-1,
                                                                    4,
                                                                    l_faultPars,
    ( double (*)[CE_N_SUB_FACES(T_SDISC.ELEMENT, ORDER)][N_CRUNS] ) l_faultData,
                                                                    RUPTURE,
                                                                    l_internal.m_faceChars );

  if(      l_faultMesh == 0 ){ EDGE_LOG_INFO << "      succesfully queried mesh for fault data"; }
  else if( l_faultMesh == 1 ){ EDGE_LOG_INFO << "      couldn't find all fault data parameters mesh, continuoing w/o"; }
  else if( l_faultMesh == 2 ){ EDGE_LOG_INFO << "      a fault data paramter doesn't match requirements, continuoing w/o"; }
  else{ EDGE_LOG_FATAL; };

  // set up rupture physics
  EDGE_LOG_INFO << "  setting up rupture physics";
   edge::elastic::setups::RuptureInit<
    T_SDISC.ELEMENT,
    ORDER,
    N_CRUNS >::alloc(
      l_enLayouts[1].nEnts,
      l_dynMem,
      l_internal.m_globalShared5[0] );

  edge::elastic::setups::RuptureInit<
    T_SDISC.ELEMENT,
    ORDER,
    N_CRUNS >::linSlipWeak( l_enLayouts[1].nEnts,
                            t_spTypeElastic::RUPTURE,
                            l_internal.m_globalShared2[0].connect.scDgAd,
                            l_internal.m_connect.faVe,
                            l_internal.m_connect.faEl,
                            l_internal.m_connect.elVe,
                            l_internal.m_connect.elFa,
                            l_internal.m_vertexChars,
                            l_internal.m_faceChars,
                            l_internal.m_elementShared1[0],
                            l_elasticConf.m_faultCrds,
                            l_elasticConf.m_lsw,
                            l_elasticConf.m_rupDoms,
                            l_elasticConf.m_stressInit,
                            l_internal.m_globalShared5[0],
       (l_faultMesh == 0) ? l_faultData : nullptr );

  delete[] l_faultData;
}

// setup star matrices
edge::seismic::solvers::AderDg< real_base,
                                N_RELAXATION_MECHANISMS,
                                T_SDISC.ELEMENT,
                                ORDER,
                                ORDER,
                                N_CRUNS,
                                MM_KERNELS_SPARSE
                              >::setupStarM( l_internal.m_nElements,
                                             l_internal.m_vertexChars,
                                             l_internal.m_connect.elVe,
                                             l_internal.m_elementShared1,
                                             l_internal.m_elementShared4 );

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

// set up sub-cell local solvers
EDGE_CHECK_EQ( l_enLayouts[l_limPlusLayout].timeGroups.size(), 1 );
l_internal.m_globalShared8[0].alloc( l_enLayouts[l_limPlusLayout].timeGroups[0].nEntsOwn,
                                     l_dynMem );
l_internal.m_globalShared8[0].init(   0,
                                      l_enLayouts[l_limPlusLayout].timeGroups[0].nEntsOwn,
                                      l_internal.m_globalShared2[0].connect.lpEl,
                                     &l_gIdsEl[0],
                                      l_internal.m_connect.elVe,
                                      l_internal.m_connect.elFa,
                                      l_internal.m_connect.elFaEl,
                                      l_internal.m_vertexChars,
                                      l_internal.m_faceChars,
                                      l_internal.m_elementChars,
                                      l_internal.m_elementShared1[0] );

// get time step stats
edge::elastic::common::getTimeStepStatsCFL( l_internal.m_nElements,
                                            l_internal.m_elementChars,
                                            l_internal.m_elementShared1,
                                            l_dT[0], l_dT[1], l_dT[2] );

// setup shared memory parallelization
for( int_tg l_tg = 0; l_tg < l_enLayouts[2].timeGroups.size(); l_tg++ ) {
  int_spType l_spType[5] = { RECEIVER, SOURCE, LIMIT, LIMIT_PLUS, EXTREMA };

  // local inner-elements
  l_shared.regWrkRgn( l_tg,
                      0,
                      l_tg * N_ENTRIES_CONTROL_FLOW + 0,
                      l_enLayouts[2].timeGroups[l_tg].inner.first,
                      l_enLayouts[2].timeGroups[l_tg].inner.size,
                      l_tg,
                      1, l_spType, l_internal.m_elementChars );

  // local send-elements
  l_shared.regWrkRgn( l_tg,
                      0,
                      l_tg * N_ENTRIES_CONTROL_FLOW + 1,
                      l_enLayouts[2].timeGroups[l_tg].inner.first+
                      l_enLayouts[2].timeGroups[l_tg].inner.size,
                      l_enLayouts[2].timeGroups[l_tg].nEntsOwn-
                      l_enLayouts[2].timeGroups[l_tg].inner.size,
                      l_enLayouts[2].timeGroups.size() + l_tg,
                      1, l_spType, l_internal.m_elementChars );

  // neigh inner-elements
  l_shared.regWrkRgn( l_tg,
                      2,
                      l_tg * N_ENTRIES_CONTROL_FLOW + 5,
                      l_enLayouts[2].timeGroups[l_tg].inner.first,
                      l_enLayouts[2].timeGroups[l_tg].inner.size,
                      l_tg,
                      3, l_spType+2, l_internal.m_elementChars );

  // neigh send-elements
  l_shared.regWrkRgn( l_tg,
                      2,
                      l_tg * N_ENTRIES_CONTROL_FLOW + 6,
                      l_enLayouts[2].timeGroups[l_tg].inner.first+
                      l_enLayouts[2].timeGroups[l_tg].inner.size,
                      l_enLayouts[2].timeGroups[l_tg].nEntsOwn-
                      l_enLayouts[2].timeGroups[l_tg].inner.size,
                      l_enLayouts[2].timeGroups.size() + l_tg,
                      3, l_spType+2, l_internal.m_elementChars );

    // point sources inner-elements
    l_shared.regWrkRgn( l_tg,
                        3,
                        l_tg * N_ENTRIES_CONTROL_FLOW + 9,
                        l_enLayouts[l_srcLayout].timeGroups[l_tg].inner.first,
                        l_enLayouts[l_srcLayout].timeGroups[l_tg].inner.size,
                        l_tg );

    // point sources send-elements
    l_shared.regWrkRgn( l_tg,
                        3,
                        l_tg * N_ENTRIES_CONTROL_FLOW + 10,
                        l_enLayouts[l_srcLayout].timeGroups[l_tg].inner.first+
                        l_enLayouts[l_srcLayout].timeGroups[l_tg].inner.size,
                        l_enLayouts[l_srcLayout].timeGroups[l_tg].nEntsOwn-
                        l_enLayouts[l_srcLayout].timeGroups[l_tg].inner.size,
                        l_enLayouts[l_srcLayout].timeGroups.size() + l_tg );

    l_spType[0] = { RECEIVER };
    // rupture physics inner-faces
    l_shared.regWrkRgn( l_tg,
                        1,
                        l_tg * N_ENTRIES_CONTROL_FLOW + 3,
                        l_enLayouts[l_rupLayoutFa].timeGroups[l_tg].inner.first,
                        l_enLayouts[l_rupLayoutFa].timeGroups[l_tg].inner.size,
                        l_tg,
                        1, l_spType, l_internal.m_globalShared7[0].bfChars );

    // rupture physics send- and receive-faces
    l_shared.regWrkRgn( l_tg,
                        1,
                        l_tg * N_ENTRIES_CONTROL_FLOW + 4,
                        l_enLayouts[l_rupLayoutFa].timeGroups[l_tg].inner.first+
                        l_enLayouts[l_rupLayoutFa].timeGroups[l_tg].inner.size,
                        l_enLayouts[l_rupLayoutFa].timeGroups[l_tg].nEntsOwn-
                        l_enLayouts[l_rupLayoutFa].timeGroups[l_tg].inner.size+
                        l_enLayouts[l_rupLayoutFa].timeGroups[l_tg].nEntsNotOwn,
                        l_enLayouts[l_rupLayoutFa].timeGroups.size() + l_tg,
                        1, l_spType, l_internal.m_globalShared7[0].bfChars );

  // limit, inner-elements
  l_shared.regWrkRgn( l_tg,
                      4,
                      l_tg * N_ENTRIES_CONTROL_FLOW + 11,
                      l_enLayouts[l_limPlusLayout].timeGroups[l_tg].inner.first,
                      l_enLayouts[l_limPlusLayout].timeGroups[l_tg].inner.size,
                      l_tg );

  // limit, send-elements
  l_shared.regWrkRgn( l_tg,
                      4,
                      l_tg * N_ENTRIES_CONTROL_FLOW + 12,
                      l_enLayouts[l_limPlusLayout].timeGroups[l_tg].inner.first+
                      l_enLayouts[l_limPlusLayout].timeGroups[l_tg].inner.size,
                      l_enLayouts[l_limPlusLayout].timeGroups[l_tg].nEntsOwn-
                      l_enLayouts[l_limPlusLayout].timeGroups[l_tg].inner.size,
                      l_enLayouts[l_limPlusLayout].timeGroups.size() + l_tg );
}

// TODO: work on initialization
// init admissibity to true
int_el l_li = 0;
for( int_el l_el = 0; l_el < l_internal.m_nElements; l_el++ ) {
  if( (l_internal.m_elementChars[l_el].spType & LIMIT) == LIMIT ) {
    for( unsigned short l_cr = 0; l_cr < N_CRUNS; l_cr++ ) {
      l_internal.m_globalShared2[0].adm[0][l_li][l_cr] = true;
      l_internal.m_globalShared2[0].adm[1][l_li][l_cr] = true;
      l_internal.m_globalShared2[0].adm[2][l_li][l_cr] = true;

      l_internal.m_globalShared2[0].lock[l_li][l_cr] = false;
    }
    l_li++;
  }
}

// setup writer for rupture boundaries
edge::io::InternalBoundary< int_el,
                            T_SDISC.ELEMENT,
                            ORDER > l_rupWriter( l_config.m_iBndFile,
                                                 l_config.m_iBndType=="vtk_binary" );

l_rupWriter.alloc( l_enLayouts[l_rupLayoutFa].nEnts,
                   5*N_CRUNS + 4 * (N_DIM-1) * N_CRUNS,
                   l_dynMem );

l_rupWriter.init( l_enLayouts[l_rupLayoutFa].nEnts,
                  l_internal.m_globalShared2[0].connect.scSv,
                  l_internal.m_globalShared7[0].connect.bfBe,
                  l_internal.m_globalShared7[0].connect.beEl,
                  l_internal.m_connect.elVe,
                  l_internal.m_globalShared2[0].svChars,
                  l_internal.m_vertexChars,
                  l_internal.m_globalShared7[0].bfChars );
