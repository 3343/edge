/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2016-2017, Regents of the University of California
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
 * Setup for the advection equation.
 **/

// set up sparse types
{
  // inherit LTS DOFs, limit plus and extrema, which are supersets of limited elements
  int_spType l_limitSup[3] = { t_ltsEl::EL_DOFS, t_enTypeShared::LIMIT_PLUS, t_enTypeShared::EXTREMA };
  for( unsigned short l_ih = 0; l_ih < 3; l_ih++ ) {
    edge::data::SparseEntities::inherit(              l_internal.m_nElements,
                                         (int_spType) t_enTypeShared::LIMIT,
                                                      l_limitSup[l_ih],
                                                      l_internal.m_elementChars );
  }

  // require DOFs for all elements adjacent to possibly limited elements
  edge::data::SparseEntities::propAdj(              l_internal.m_nElements,
                                                    C_ENT[T_SDISC.ELEMENT].N_FACES,
                                                    l_internal.m_connect.elFaEl[0],
                                       (int_spType) t_enTypeShared::LIMIT,
                                       (int_spType) C_LTS_EL[t_ltsEl::EL_DOFS],
                                                    l_internal.m_elementChars,
                                                    l_internal.m_elementChars );
  // limit plus are the limited elements and their face neighbors
  edge::data::SparseEntities::propAdj(              l_internal.m_nElements,
                                                    C_ENT[T_SDISC.ELEMENT].N_FACES,
                                                    l_internal.m_connect.elFaEl[0],
                                       (int_spType) t_enTypeShared::LIMIT,
                                       (int_spType) t_enTypeShared::LIMIT_PLUS,
                                                    l_internal.m_elementChars,
                                                    l_internal.m_elementChars );
  // store extrema for all limited elements and those adjacent through vertices
  edge::data::SparseEntities::propAdj(              l_internal.m_nElements,
                                                    l_internal.m_connect.elVeEl,
                                       (int_spType) t_enTypeShared::LIMIT,
                                       (int_spType) t_enTypeShared::EXTREMA,
                                                    l_internal.m_elementChars,
                                                    l_internal.m_elementChars );
}

// allocate flex data for time buffers and DOFs
{
  int_spType  l_spTypes[1] = { C_LTS_EL[t_ltsEl::EL_DOFS] };
  std::size_t l_spSizes[1] = { 1 };

  real_base (**l_raw)[N_ELEMENT_MODES][N_CRUNS];

  l_raw = l_dynMem.flex<
    real_base [N_ELEMENT_MODES][N_CRUNS] >( l_internal.m_nElements,
                                            1,
                                            1,
                                            l_spTypes,
                                            l_spSizes,
                                            l_internal.m_elementChars );

  l_internal.m_globalShared2[0] = l_raw;
  l_internal.m_globalShared2[1] = l_raw + l_internal.m_nElements;
}


// init data of limiter
// TODO: outsource to avoid redundancy
unsigned short l_limPlusLayout;
{
  // get number of limited elements
  int_el l_nLim = edge::data::SparseEntities::nSp( l_enLayouts[2].nEnts,
                                                   t_enTypeShared::LIMIT,
                                                   l_internal.m_elementChars );

  // derive an entity layout for the limited plus elements
  l_enLayouts.resize( l_enLayouts.size() + 1 );
  edge::data::SparseEntities::denseToSparse( t_enTypeShared::LIMIT_PLUS,
                                             l_internal.m_elementChars,
                                             l_enLayouts[2],
                                             l_enLayouts.back() );
  l_limPlusLayout = l_enLayouts.size()-1;

  // get number of sub-cell extrema
  int_el l_nScExt = edge::data::SparseEntities::nSp( l_enLayouts[2].nEnts,
                                                     t_enTypeShared::EXTREMA,
                                                     l_internal.m_elementChars );

  // allocate sub-cell memory
  edge::sc::Memory< int_el,
                    real_base,
                    T_SDISC.ELEMENT,
                    ORDER,
                    N_QUANTITIES,
                    N_CRUNS >::alloc( l_nLim,
                                      l_enLayouts[l_limPlusLayout].nEnts,
                                      l_nScExt,
                                      l_dynMem,
                                      l_internal.m_globalShared3[0] );

  // link limited plus elements to dense elements (no bridge)
  edge::data::SparseEntities::linkSpDe( l_enLayouts[2].nEnts,
                                        t_enTypeShared::LIMIT_PLUS,
                                        l_internal.m_elementChars,
                                        l_internal.m_globalShared3[0].connect.lpEl );

  // link limited plus elements to limited elements (no bridge)
  edge::data::SparseEntities::linkSpSp( l_enLayouts[2].nEnts,
                                        t_enTypeShared::LIMIT_PLUS,
                                        t_enTypeShared::LIMIT,
                                        l_internal.m_elementChars,
                                        l_internal.m_globalShared3[0].connect.lpLi );

  // link limited plus elements to limited plus elements (faces as bridge)
  edge::data::SparseEntities::linkSpAdjDst( l_enLayouts[2].nEnts,
                                            C_ENT[T_SDISC.ELEMENT].N_FACES,
                                            l_internal.m_connect.elFaEl[0],
                                            t_enTypeShared::LIMIT_PLUS,
                                            t_enTypeShared::LIMIT_PLUS,
                                            l_internal.m_elementChars,
                                            l_internal.m_elementChars,
                                            l_internal.m_globalShared3[0].connect.lpFaLp[0] );

  // link limited elements to extrema (no bridge)
  edge::data::SparseEntities::linkSpSp( l_enLayouts[2].nEnts,
                                        t_enTypeShared::LIMIT,
                                        t_enTypeShared::EXTREMA,
                                        l_internal.m_elementChars,
                                        l_internal.m_globalShared3[0].connect.liEx );

  // link limited elements to sub-cell extrema (vertices as bridge)
  int_el l_nLinkRaw;
  int_el l_nLinkPtr;

  // get number of raw entries and pointers
  edge::data::SparseEntities::nLinkSpAdjDst( l_enLayouts[2].nEnts,
                                             l_internal.m_connect.elVeEl,
                                             t_enTypeShared::LIMIT,
                                             t_enTypeShared::EXTREMA,
                                             l_internal.m_elementChars,
                                             l_internal.m_elementChars,
                                             l_nLinkRaw,
                                             l_nLinkPtr );

  // check that the number of pointers match our expectations
  EDGE_CHECK_EQ( l_nLinkPtr, l_nLim );

  // allocate data for the raw entries
  int_el *l_raw = (int_el *) l_dynMem.allocate( l_nLinkRaw * sizeof(int_el) );

  // derive link
  edge::data::SparseEntities::linkSpAdjDst( l_enLayouts[2].nEnts,
                                            l_internal.m_connect.elVeEl,
                                            t_enTypeShared::LIMIT,
                                            t_enTypeShared::EXTREMA,
                                            l_internal.m_elementChars,
                                            l_internal.m_elementChars,
                                            l_raw,
                                            l_internal.m_globalShared3[0].connect.liVeEx );
}

#ifdef PP_USE_MPI
  // init mpi layout
  EDGE_CHECK( l_enLayouts[2].timeGroups.size() == 1 );

  ///////////////////////////////////////////////
  // TODO: fix MPI support of sub-cell limiter //
  ///////////////////////////////////////////////
  l_mpi.initLayout( l_enLayouts[2],
                    l_internal.m_elementModePrivate2[0][0][0],
                    N_ELEMENT_MODES*N_CRUNS*sizeof(t_elementModePrivate2),
                    0,
                    1 );
#endif

// setup shared memory parallelization
for( int_tg l_tg = 0; l_tg < l_enLayouts[2].timeGroups.size(); l_tg++ ) {
  int_spType l_spTypes[2] = { LIMIT, EXTREMA };

  // local, inner-elements
  l_shared.regWrkRgn( l_tg,
                      0,
                      l_tg * N_ENTRIES_CONTROL_FLOW + 0,
                      l_enLayouts[2].timeGroups[l_tg].inner.first,
                      l_enLayouts[2].timeGroups[l_tg].inner.size,
                      l_tg,
                      1, l_spTypes, l_internal.m_elementChars );

  // local, send-elements
  l_shared.regWrkRgn( l_tg,
                      0,
                      l_tg * N_ENTRIES_CONTROL_FLOW + 1,
                      l_enLayouts[2].timeGroups[l_tg].inner.first+
                      l_enLayouts[2].timeGroups[l_tg].inner.size,
                      l_enLayouts[2].timeGroups[l_tg].nEntsOwn-
                      l_enLayouts[2].timeGroups[l_tg].inner.size,
                      l_enLayouts[2].timeGroups.size() + l_tg,
                      1, l_spTypes, l_internal.m_elementChars );

  // neigh, inner-elements
  l_shared.regWrkRgn( l_tg,
                      1,
                      l_tg * N_ENTRIES_CONTROL_FLOW + 3,
                      l_enLayouts[2].timeGroups[l_tg].inner.first,
                      l_enLayouts[2].timeGroups[l_tg].inner.size,
                      l_tg,
                      2, l_spTypes, l_internal.m_elementChars );

  // neigh, send-elements
  l_shared.regWrkRgn( l_tg,
                      1,
                      l_tg * N_ENTRIES_CONTROL_FLOW + 4,
                      l_enLayouts[2].timeGroups[l_tg].inner.first+
                      l_enLayouts[2].timeGroups[l_tg].inner.size,
                      l_enLayouts[2].timeGroups[l_tg].nEntsOwn-
                      l_enLayouts[2].timeGroups[l_tg].inner.size,
                      l_enLayouts[2].timeGroups.size() + l_tg,
                      2, l_spTypes, l_internal.m_elementChars );

  // limit, inner-elements
  l_shared.regWrkRgn( l_tg,
                      2,
                      l_tg * N_ENTRIES_CONTROL_FLOW + 6,
                      l_enLayouts[l_limPlusLayout].timeGroups[l_tg].inner.first,
                      l_enLayouts[l_limPlusLayout].timeGroups[l_tg].inner.size,
                      l_tg );

  // limit, send-elements
  l_shared.regWrkRgn( l_tg,
                      2,
                      l_tg * N_ENTRIES_CONTROL_FLOW + 7,
                      l_enLayouts[l_limPlusLayout].timeGroups[l_tg].inner.first+
                      l_enLayouts[l_limPlusLayout].timeGroups[l_tg].inner.size,
                      l_enLayouts[l_limPlusLayout].timeGroups[l_tg].nEntsOwn-
                      l_enLayouts[l_limPlusLayout].timeGroups[l_tg].inner.size,
                      l_enLayouts[l_limPlusLayout].timeGroups.size() + l_tg );
}

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

  // SC
  edge::setups::InitialDofs<
    T_SDISC.ELEMENT,
    ORDER,
    N_QUANTITIES,
    N_CRUNS >::sc( 0,
                   l_internal.m_nElements,
                   0,
                   t_enTypeShared::LIMIT,
                   l_config.m_initValsExprStrs,
                   l_internal.m_connect.elVe,
                   l_internal.m_globalShared3[0].connect.scSv,
                   l_internal.m_vertexChars,
                   l_internal.m_globalShared3[0].svChars,
                   l_internal.m_elementChars,
                   l_internal.m_globalShared3[0].dofs );
}

// init sub-cell extrema and admissibility
{
  int_el l_spLim = 0;
  int_el l_spExt = 0;

  for( int_el l_el = 0; l_el < l_internal.m_nElements; l_el++ ) {
    //! TODO: Use dedicated scratch memory for this
    real_base l_scratch[1][CE_N_SUB_CELLS( T_SDISC.ELEMENT, ORDER )][N_CRUNS];

    if( (l_internal.m_elementChars[l_el].spType & EXTREMA) == EXTREMA ) {
      // init extrema with DG-extrema
      edge::sc::Kernels< T_SDISC.ELEMENT,
                         ORDER,
                         1,
                         N_CRUNS >::dgExtremaVanilla( l_internal.m_elementModePrivate1[l_el],
                                                      l_internal.m_globalShared3[0].ops.scatter,
                                                      l_scratch,
                                                      l_internal.m_globalShared3[0].ext[0][l_spExt][0],
                                                      l_internal.m_globalShared3[0].ext[0][l_spExt][1] );

      // copy DG-extrema as candidate solution
      for( unsigned short l_ex = 0; l_ex < 2; l_ex++ )
        for( unsigned short l_qt = 0; l_qt < N_QUANTITIES; l_qt++ )
          for( unsigned short l_cr = 0; l_cr < N_CRUNS; l_cr++ )
            l_internal.m_globalShared3[0].ext[1][l_spExt][l_ex][l_qt][l_cr] =
              l_internal.m_globalShared3[0].ext[0][l_spExt][l_ex][l_qt][l_cr];

      if( (l_internal.m_elementChars[l_el].spType & LIMIT) == LIMIT ) {
        // overwrite extrema of limited elements with sub-cell extrema
        edge::sc::Kernels< T_SDISC.ELEMENT,
                           ORDER,
                           1,
                           N_CRUNS >::scExtrema(  l_internal.m_globalShared3[0].dofs[l_spLim],
                                                  l_internal.m_globalShared3[0].ext[0][l_spExt][0],
                                                  l_internal.m_globalShared3[0].ext[0][l_spExt][1] );

        l_spLim++;
      }

      l_spExt++;
    }
    // every limited element is required to compute extrema
    else EDGE_CHECK_NE( (l_internal.m_elementChars[l_el].spType & LIMIT), LIMIT );
  }

  int_el l_li = 0;
  for( int_el l_el = 0; l_el < l_internal.m_nElements; l_el++ ) {
    if( (l_internal.m_elementChars[l_el].spType & LIMIT) == LIMIT ) {
      int_el l_ex = l_internal.m_globalShared3[0].connect.liEx[l_li];

      edge::sc::Detections< 1,
                            N_CRUNS >::dmp( l_internal.m_globalShared3[0].ext[0][l_ex],
                                            l_internal.m_globalShared3[0].ext[0],
                                            l_internal.m_globalShared3[0].ext[1][l_ex],
                                            l_internal.m_globalShared3[0].connect.liVeEx[l_li+1]-
                                            l_internal.m_globalShared3[0].connect.liVeEx[l_li],
                                            l_internal.m_globalShared3[0].connect.liVeEx[l_li],
                                            l_internal.m_globalShared3[0].adm[0][l_li] );

      l_li++;
    }
  }


}

// initialize wave speeds
edge::advection::setups::Convergence::setConstantSpeed( l_internal.m_nElements,
                                                        l_internal.m_elementShared1,
                                                        1,
                                                        0.5,  // ignored for #dim==1
                                                        -1 ); // ignored for #dim!=3

#if PP_ORDER > 1
// set up star "matrices"
edge::advection::solvers::AderDg<
    T_SDISC.ELEMENT,
    ORDER,
    ORDER,
    N_CRUNS >::setupStarM( l_internal.m_nElements,
                           l_internal.m_vertexChars,
                           l_internal.m_connect.elVe,
                           l_internal.m_elementShared1,
                           l_internal.m_elementShared2 );
#endif

// get time step stats
edge::advection::solvers::common::getTimeStepStatistics( l_internal.m_nElements,
                                                         l_internal.m_elementChars,
                                                         l_internal.m_elementShared1,
                                                         l_dT[0],
                                                         l_dT[1],
                                                         l_dT[2] );

// set up flux solvers
edge::advection::solvers::common::setupSolvers( l_internal.m_nElements,
                                                l_internal.m_nFaces,
                                                l_mesh.getInMap()->faMeDa,
                                                l_mesh.getInMap()->faDaMe,
                                                l_mesh.getInMap()->elMeDa,
                                                l_mesh.getInMap()->elDaMe,
                                                l_internal.m_connect.elVe,
                                                l_internal.m_connect.faEl,
                                                l_internal.m_connect.elFa,
                                                l_internal.m_vertexChars,
                                                l_internal.m_faceChars,
                                                l_internal.m_elementShared1,
                                                l_internal.m_elementShared3 );
