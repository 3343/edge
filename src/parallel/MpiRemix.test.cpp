/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2019, Alexander Breuer
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
 * Tests the MPI interface.
 **/
#include <catch.hpp>
#define private public
#include "MpiRemix.h"
#undef private

#ifdef PP_USE_MPI
TEST_CASE( "Tests the setup of MPI-related data structures.", "[MpiRemix][data]" ) {
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

  edge::data::Dynamic l_dynMem;

  edge::parallel::MpiRemix l_mpi( 0,
                                  nullptr );
  l_mpi.init( 7,
              4,
              200,
              17,
              l_commStruct,
              l_sendFa,
              l_sendEl,
              l_recvFa,
              l_recvEl,
              l_dynMem );

  REQUIRE( l_mpi.m_nChs == 2 );

  // check some send pointers
  REQUIRE( l_mpi.m_sendPtrs[0*4 + 0] == l_mpi.m_sendBuffer+10*17  );
  REQUIRE( l_mpi.m_sendPtrs[0*4 + 1] == nullptr                   );
  REQUIRE( l_mpi.m_sendPtrs[0*4 + 2] == nullptr                   );
  REQUIRE( l_mpi.m_sendPtrs[0*4 + 3] == nullptr                   );

  REQUIRE( l_mpi.m_sendPtrs[1*4 + 0] == nullptr                   );
  REQUIRE( l_mpi.m_sendPtrs[1*4 + 1] == l_mpi.m_sendBuffer+ 4*17  );
  REQUIRE( l_mpi.m_sendPtrs[1*4 + 2] == nullptr                   );
  REQUIRE( l_mpi.m_sendPtrs[1*4 + 3] == nullptr                   );

  REQUIRE( l_mpi.m_sendPtrs[2*4 + 0] == nullptr                   );
  REQUIRE( l_mpi.m_sendPtrs[2*4 + 1] == l_mpi.m_sendBuffer+12*17  );
  REQUIRE( l_mpi.m_sendPtrs[2*4 + 2] == l_mpi.m_sendBuffer+ 5*17  );
  REQUIRE( l_mpi.m_sendPtrs[2*4 + 3] == l_mpi.m_sendBuffer+ 8*17  );
  // [...]

  // check some recv pointers
  REQUIRE( l_mpi.m_recvPtrs[0*4 + 0] == l_mpi.m_recvBuffer+10*17*2  );
  REQUIRE( l_mpi.m_recvPtrs[0*4 + 1] == nullptr                     );
  REQUIRE( l_mpi.m_recvPtrs[0*4 + 2] == nullptr                     );
  REQUIRE( l_mpi.m_recvPtrs[0*4 + 3] == l_mpi.m_recvBuffer+ 6*17*2  );

  REQUIRE( l_mpi.m_recvPtrs[1*4 + 0] == nullptr                     );
  REQUIRE( l_mpi.m_recvPtrs[1*4 + 1] == nullptr                     );
  REQUIRE( l_mpi.m_recvPtrs[1*4 + 2] == nullptr                     );
  REQUIRE( l_mpi.m_recvPtrs[1*4 + 3] == nullptr                     );

  REQUIRE( l_mpi.m_sendMsgs[0].lt   == true                         );
  REQUIRE( l_mpi.m_sendMsgs[0].tg   == 0                            );
  REQUIRE( l_mpi.m_sendMsgs[0].rank == 2                            );
  REQUIRE( l_mpi.m_sendMsgs[0].size == 15*17                        );
  REQUIRE( l_mpi.m_sendMsgs[0].tag  == 4                            );
  REQUIRE( l_mpi.m_sendMsgs[0].ptr  == l_mpi.m_sendBuffer + 0       );

  REQUIRE( l_mpi.m_sendMsgs[1].lt   == false                        );
  REQUIRE( l_mpi.m_sendMsgs[1].tg   == 3                            );
  REQUIRE( l_mpi.m_sendMsgs[1].rank == 1                            );
  REQUIRE( l_mpi.m_sendMsgs[1].size == 5*17*2                       );
  REQUIRE( l_mpi.m_sendMsgs[1].tag  == 3*7+1                        );
  REQUIRE( l_mpi.m_sendMsgs[1].ptr  == l_mpi.m_sendBuffer + 15*17   );

  REQUIRE( l_mpi.m_recvMsgs[0].lt   == true                         );
  REQUIRE( l_mpi.m_recvMsgs[0].tg   == 0                            );
  REQUIRE( l_mpi.m_recvMsgs[0].rank == 2                            );
  REQUIRE( l_mpi.m_recvMsgs[0].size == 15*17*2                      );
  REQUIRE( l_mpi.m_recvMsgs[0].tag  == 4*7                          );
  REQUIRE( l_mpi.m_recvMsgs[0].ptr  == l_mpi.m_recvBuffer + 0       );

  REQUIRE( l_mpi.m_recvMsgs[1].lt   == false                        );
  REQUIRE( l_mpi.m_recvMsgs[1].tg   == 3                            );
  REQUIRE( l_mpi.m_recvMsgs[1].rank == 1                            );
  REQUIRE( l_mpi.m_recvMsgs[1].size == 5*17                         );
  REQUIRE( l_mpi.m_recvMsgs[1].tag  == 1*7+3                        );
  REQUIRE( l_mpi.m_recvMsgs[1].ptr  == l_mpi.m_recvBuffer + 15*17*2 );
}
#endif