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
 * Setup for the receivers.
 **/

// receiver-output
edge::io::Receivers l_receivers;
edge::io::ReceiversSf<real_base, T_SDISC.ELEMENT, ORDER, N_CRUNS> l_recvsSf; // implementation-dependent #quantities, setup is in impl

// setup receivers
if( l_config.m_recvCrds[0].size() > 0 ) {
  EDGE_LOG_INFO << "searching for receivers, those outside will be projected to the mesh boundaries..";

  // init receivers and print info
  l_receivers.init(                     T_SDISC.ELEMENT,
                                        l_config.m_recvCrds[0].size(),
                                        l_config.m_recvPath[0],
                                       &l_config.m_recvNames[0][0],
                    (real_mesh (*)[3]) &l_config.m_recvCrds[0][0][0],
                                        l_config.m_recvFreq[0],
                                        l_enLayouts[2],
                                        l_internal.m_connect.elVe[0],
                                        l_internal.m_vertexChars );
   l_receivers.print( l_config.m_recvCrds[0].size() );

   // get dense-entities with receivers
   std::vector< int_el > l_enRecv;
   l_receivers.getEnRecv( l_enRecv );

   // iterate over entities with receivers and enable the receiver-flag in their sparse type
   for( std::size_t l_re = 0; l_re < l_enRecv.size(); l_re++ ) {
     l_internal.m_elementChars[l_enRecv[l_re]].spType |= RECEIVER;
   }
}
