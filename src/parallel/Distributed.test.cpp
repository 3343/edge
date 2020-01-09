/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2019-2020, Alexander Breuer
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
 * Tests the distributed memory parallelization interface.
 **/
#include <catch.hpp>
#define protected public
#include "DistributedDummy.hpp"
#undef protected

std::size_t l_commStruct[9] = {2, 0, 2, 4, 15, 3, 1, 1, 5};

unsigned short l_sendFa[20] = { 0, 3, 2, 0, 1,
                                2, 3, 3, 3, 1,
                                0, 2, 1, 2, 2,
                                3, 2, 0, 1, 1 };

std::size_t l_sendEl[20] =    { 4, 3, 5, 9, 1,
                                2, 4, 9, 2, 5,
                                0, 7, 2, 8, 3,
                                7, 4, 5, 9, 3 };

unsigned short l_recvFa[20] = { 3, 3, 2, 0, 1,
                                2, 3, 1, 2, 2,
                                0, 2, 3, 2, 3,
                                2, 2, 0, 1, 1 };

std::size_t l_recvEl[20] =    { 4, 3, 5, 9, 3,
                                6, 0, 6, 5, 2,
                                0, 9, 7, 8, 3,
                                7, 4, 2, 9, 4 };

TEST_CASE( "Tests the setup of distributed memory-related data structures.", "[Distributed][data]" ) {
  edge::data::Dynamic l_dynMem;

  edge::parallel::DistributedDummy l_dist( 0, nullptr );
  l_dist.init( 7,
               4,
               200,
               17,
               l_commStruct,
               l_sendFa,
               l_sendEl,
               l_recvFa,
               l_recvEl,
               l_dynMem );

  REQUIRE( l_dist.m_nChs == 2 );

  // check sizes of buffers (tg0 -> tg4: 15, tg3 -> tg1: 5 )
  REQUIRE( l_dist.m_sendBufferSize == (15 + 5*2)* // faces w.r.t. to time stepping
                                      17 );       // size per comm struct
  REQUIRE( l_dist.m_recvBufferSize == (15*2 + 5)* // faces w.r.t. to time stepping
                                      17 );       // size per comm struct

  // check offset of pointers
  REQUIRE( l_dist.m_sendPtrs[1] == l_dist.m_sendPtrs[0] + 200*4 );
  REQUIRE( l_dist.m_recvPtrs[1] == l_dist.m_recvPtrs[0] + 200*4 );

  // check some send pointers
  REQUIRE( l_dist.m_sendPtrs[0][0*4 + 0] == l_dist.m_sendBuffers+10*17 );
  REQUIRE( l_dist.m_sendPtrs[0][0*4 + 1] == nullptr                    );
  REQUIRE( l_dist.m_sendPtrs[0][0*4 + 2] == nullptr                    );
  REQUIRE( l_dist.m_sendPtrs[0][0*4 + 3] == nullptr                    );

  REQUIRE( l_dist.m_sendPtrs[0][1*4 + 0] == nullptr                    );
  REQUIRE( l_dist.m_sendPtrs[0][1*4 + 1] == l_dist.m_sendBuffers+ 4*17 );
  REQUIRE( l_dist.m_sendPtrs[0][1*4 + 2] == nullptr                    );
  REQUIRE( l_dist.m_sendPtrs[0][1*4 + 3] == nullptr                    );

  REQUIRE( l_dist.m_sendPtrs[0][2*4 + 0] == nullptr                    );
  REQUIRE( l_dist.m_sendPtrs[0][2*4 + 1] == l_dist.m_sendBuffers+12*17 );
  REQUIRE( l_dist.m_sendPtrs[0][2*4 + 2] == l_dist.m_sendBuffers+ 5*17 );
  REQUIRE( l_dist.m_sendPtrs[0][2*4 + 3] == l_dist.m_sendBuffers+ 8*17 );
  // [...]
  REQUIRE( l_dist.m_sendPtrs[1][0*4 + 0] == l_dist.m_sendBuffers + l_dist.m_sendBufferSize + 10*17 );
  REQUIRE( l_dist.m_sendPtrs[1][0*4 + 1] == nullptr                                                );
  REQUIRE( l_dist.m_sendPtrs[1][0*4 + 2] == nullptr                                                );
  REQUIRE( l_dist.m_sendPtrs[1][0*4 + 3] == nullptr                                                );

  REQUIRE( l_dist.m_sendPtrs[1][1*4 + 0] == nullptr                                                );
  REQUIRE( l_dist.m_sendPtrs[1][1*4 + 1] == l_dist.m_sendBuffers + l_dist.m_sendBufferSize +  4*17 );
  REQUIRE( l_dist.m_sendPtrs[1][1*4 + 2] == nullptr                                                );
  REQUIRE( l_dist.m_sendPtrs[1][1*4 + 3] == nullptr                                                );

  REQUIRE( l_dist.m_sendPtrs[1][2*4 + 0] == nullptr                      );
  REQUIRE( l_dist.m_sendPtrs[1][2*4 + 1] == l_dist.m_sendBuffers + l_dist.m_sendBufferSize + 12*17 );
  REQUIRE( l_dist.m_sendPtrs[1][2*4 + 2] == l_dist.m_sendBuffers + l_dist.m_sendBufferSize +  5*17 );
  REQUIRE( l_dist.m_sendPtrs[1][2*4 + 3] == l_dist.m_sendBuffers + l_dist.m_sendBufferSize +  8*17 );
  // [...]

  // check some recv pointers
  REQUIRE( l_dist.m_recvPtrs[0][0*4 + 0] == l_dist.m_recvBuffers+10*17*2                             );
  REQUIRE( l_dist.m_recvPtrs[0][0*4 + 1] == nullptr                                                  );
  REQUIRE( l_dist.m_recvPtrs[0][0*4 + 2] == nullptr                                                  );
  REQUIRE( l_dist.m_recvPtrs[0][0*4 + 3] == l_dist.m_recvBuffers+ 6*17*2                             );

  REQUIRE( l_dist.m_recvPtrs[0][1*4 + 0] == nullptr                                                  );
  REQUIRE( l_dist.m_recvPtrs[0][1*4 + 1] == nullptr                                                  );
  REQUIRE( l_dist.m_recvPtrs[0][1*4 + 2] == nullptr                                                  );
  REQUIRE( l_dist.m_recvPtrs[0][1*4 + 3] == nullptr                                                  );
  // [...]
  REQUIRE( l_dist.m_recvPtrs[1][0*4 + 0] == l_dist.m_recvBuffers + 10*17*2                           );
  REQUIRE( l_dist.m_recvPtrs[1][0*4 + 1] == nullptr                                                  );
  REQUIRE( l_dist.m_recvPtrs[1][0*4 + 2] == nullptr                                                  );
  REQUIRE( l_dist.m_recvPtrs[1][0*4 + 3] == l_dist.m_recvBuffers +  6*17*2 );

  REQUIRE( l_dist.m_recvPtrs[1][1*4 + 0] == nullptr                                                  );
  REQUIRE( l_dist.m_recvPtrs[1][1*4 + 1] == nullptr                                                  );
  REQUIRE( l_dist.m_recvPtrs[1][1*4 + 2] == nullptr                                                  );
  REQUIRE( l_dist.m_recvPtrs[1][1*4 + 3] == nullptr                                                  );
  // [...]
  REQUIRE( l_dist.m_recvPtrs[2][0*4 + 0] == l_dist.m_recvBuffers + l_dist.m_recvBufferSize + 10*17*2 );
  REQUIRE( l_dist.m_recvPtrs[2][0*4 + 1] == nullptr                                                  );
  REQUIRE( l_dist.m_recvPtrs[2][0*4 + 2] == nullptr                                                  );
  REQUIRE( l_dist.m_recvPtrs[2][0*4 + 3] == l_dist.m_recvBuffers + l_dist.m_recvBufferSize +  6*17*2 );

  REQUIRE( l_dist.m_recvPtrs[2][1*4 + 0] == nullptr                                                  );
  REQUIRE( l_dist.m_recvPtrs[2][1*4 + 1] == nullptr                                                  );
  REQUIRE( l_dist.m_recvPtrs[2][1*4 + 2] == nullptr                                                  );
  REQUIRE( l_dist.m_recvPtrs[2][1*4 + 3] == nullptr                                                  );
  // [...]
  REQUIRE( l_dist.m_recvPtrs[3][0*4 + 0] == l_dist.m_recvBuffers + l_dist.m_recvBufferSize + 10*17*2 );
  REQUIRE( l_dist.m_recvPtrs[3][0*4 + 1] == nullptr                                                  );
  REQUIRE( l_dist.m_recvPtrs[3][0*4 + 2] == nullptr                                                  );
  REQUIRE( l_dist.m_recvPtrs[3][0*4 + 3] == l_dist.m_recvBuffers + l_dist.m_recvBufferSize +  6*17*2 );

  REQUIRE( l_dist.m_recvPtrs[3][1*4 + 0] == nullptr                                                  );
  REQUIRE( l_dist.m_recvPtrs[3][1*4 + 1] == nullptr                                                  );
  REQUIRE( l_dist.m_recvPtrs[3][1*4 + 2] == nullptr                                                  );
  REQUIRE( l_dist.m_recvPtrs[3][1*4 + 3] == nullptr                                                  );

  REQUIRE( l_dist.m_sendMsgs[0].tgL  == 0       );
  REQUIRE( l_dist.m_sendMsgs[0].tgR  == 4       );
  REQUIRE( l_dist.m_sendMsgs[0].rank == 2       );
  REQUIRE( l_dist.m_sendMsgs[0].size == 15*17   );
  REQUIRE( l_dist.m_sendMsgs[0].offL == 0       );

  REQUIRE( l_dist.m_sendMsgs[1].tgL  == 3       );
  REQUIRE( l_dist.m_sendMsgs[1].tgR  == 1       );
  REQUIRE( l_dist.m_sendMsgs[1].rank == 1       );
  REQUIRE( l_dist.m_sendMsgs[1].size == 5*17*2  );
  REQUIRE( l_dist.m_sendMsgs[1].offL == 15*17   );

  REQUIRE( l_dist.m_recvMsgs[0].tgL  == 0       );
  REQUIRE( l_dist.m_recvMsgs[0].tgR  == 4       );
  REQUIRE( l_dist.m_recvMsgs[0].rank == 2       );
  REQUIRE( l_dist.m_recvMsgs[0].size == 15*17*2 );
  REQUIRE( l_dist.m_recvMsgs[0].offL == 0       );

  REQUIRE( l_dist.m_recvMsgs[1].tgL  == 3       );
  REQUIRE( l_dist.m_recvMsgs[1].tgR  == 1       );
  REQUIRE( l_dist.m_recvMsgs[1].rank == 1       );
  REQUIRE( l_dist.m_recvMsgs[1].size == 5*17    );
  REQUIRE( l_dist.m_recvMsgs[1].offL == 15*17*2 );
}

TEST_CASE( "Tests the derivation of send comm buffer ids for double-buffered schemes.", "[Distributed][sendCommBuffers2]" ) {
  edge::data::Dynamic l_dynMem;
  edge::parallel::DistributedDummy l_dist( 0, nullptr );
  l_dist.init( 7,
               4,
               200,
               17,
               l_commStruct,
               l_sendFa,
               l_sendEl,
               l_recvFa,
               l_recvEl,
               l_dynMem );

  unsigned short l_cbL, l_cbLtR, l_cbGeR;

  l_dist.m_nSendsSync[0] = 0;
  l_dist.sendCommBuffers2( 0,
                           l_cbL,
                           l_cbLtR,
                           l_cbGeR );
  REQUIRE( l_cbL == 0 );
  REQUIRE( l_cbLtR == std::numeric_limits< unsigned short >::max() );
  REQUIRE( l_cbGeR == 0 );

  l_dist.m_nSendsSync[0] = 1;
  l_dist.sendCommBuffers2( 0,
                           l_cbL,
                           l_cbLtR,
                           l_cbGeR );
  REQUIRE( l_cbL == 1 );
  REQUIRE( l_cbLtR == 0 );
  REQUIRE( l_cbGeR == 1 );

  l_dist.m_nSendsSync[0] = 2;
  l_dist.sendCommBuffers2( 0,
                           l_cbL,
                           l_cbLtR,
                           l_cbGeR );
  REQUIRE( l_cbL == 0 );
  REQUIRE( l_cbLtR == std::numeric_limits< unsigned short >::max() );
  REQUIRE( l_cbGeR == 0 );

  l_dist.m_nSendsSync[0] = 3;
  l_dist.sendCommBuffers2( 0,
                           l_cbL,
                           l_cbLtR,
                           l_cbGeR );
  REQUIRE( l_cbL == 1 );
  REQUIRE( l_cbLtR == 1 );
  REQUIRE( l_cbGeR == 1 );

  l_dist.m_nSendsSync[0] = 4;
  l_dist.sendCommBuffers2( 0,
                           l_cbL,
                           l_cbLtR,
                           l_cbGeR );
  REQUIRE( l_cbL == 0 );
  REQUIRE( l_cbLtR == std::numeric_limits< unsigned short >::max() );
  REQUIRE( l_cbGeR == 0 );
}

TEST_CASE( "Tests the derivation of recv comm buffer ids for double-buffered schemes.", "[Distributed][recvCommBuffers2]" ) {
  edge::data::Dynamic l_dynMem;
  edge::parallel::DistributedDummy l_dist( 0, nullptr );
  l_dist.init( 7,
               4,
               200,
               17,
               l_commStruct,
               l_sendFa,
               l_sendEl,
               l_recvFa,
               l_recvEl,
               l_dynMem );

  unsigned short l_cbLtL, l_cbGeL;

  l_dist.m_nRecvsSync[0] = 0;
  l_dist.recvCommBuffers2( 0,
                           l_cbLtL,
                           l_cbGeL );
  REQUIRE( l_cbLtL == 0 );
  REQUIRE( l_cbGeL == 0 );

  l_dist.m_nRecvsSync[0] = 1;
  l_dist.recvCommBuffers2( 0,
                           l_cbLtL,
                           l_cbGeL );
  REQUIRE( l_cbLtL == 0 );
  REQUIRE( l_cbGeL == 1 );

  l_dist.m_nRecvsSync[0] = 2;
  l_dist.recvCommBuffers2( 0,
                           l_cbLtL,
                           l_cbGeL );
  REQUIRE( l_cbLtL == 1 );
  REQUIRE( l_cbGeL == 0 );

  l_dist.m_nRecvsSync[0] = 3;
  l_dist.recvCommBuffers2( 0,
                           l_cbLtL,
                           l_cbGeL );
  REQUIRE( l_cbLtL == 1 );
  REQUIRE( l_cbGeL == 1 );

  l_dist.m_nRecvsSync[0] = 4;
  l_dist.recvCommBuffers2( 0,
                           l_cbLtL,
                           l_cbGeL );
  REQUIRE( l_cbLtL == 0 );
  REQUIRE( l_cbGeL == 0 );
}