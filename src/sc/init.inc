/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2017, Regents of the University of California
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
  EDGE_CHECK_LE( l_nLim, l_enLayouts[l_limPlusLayout].nEnts );

  // get number of sub-cell extrema
  int_el l_nScExt = edge::data::SparseEntities::nSp( l_enLayouts[2].nEnts,
                                                     t_enTypeShared::EXTREMA,
                                                     l_internal.m_elementChars );
  EDGE_CHECK_LE( l_enLayouts[l_limPlusLayout].nEnts, l_nScExt );

  // allocate sub-cell memory
  edge::sc::Memory::alloc( l_nLim,
                           l_enLayouts[l_limPlusLayout].nEnts,
                           l_nScExt,
                           l_dynMem,
                           l_internal.m_globalShared2[0] );

  // init sub-cell memory
  edge::sc::Init< T_SDISC.ELEMENT,
                  ORDER,
                  N_QUANTITIES,
                  N_CRUNS >::data( l_nLim,
                                   l_nScExt,
                                   l_internal.m_globalShared2[0].dofs,
                                   l_internal.m_globalShared2[0].surf,
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
}
