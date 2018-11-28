/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2015-2018, Regents of the University of California
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
 * Logging interface.
 **/

#ifndef EDGE_IO_LOGGING_H_
#define EDGE_IO_LOGGING_H_

#include <iostream>
#include "parallel/global.h"
#ifdef PP_USE_MPI
#include "parallel/mpi_wrapper.inc"
#endif

#ifdef PP_USE_EASYLOGGING

// enable thread-safety in omp-configs
#ifdef PP_USE_OMP
#define ELPP_THREAD_SAFE
#endif

#define ELPP_NO_DEFAULT_LOG_FILE

// silence compilers for easylogging
#pragma GCC system_header
#include <submodules/include/easylogging++.h>

#define EDGE_LOG_INFO_ALL         LOG(INFO)
#define EDGE_LOG_INFO             LOG_IF(edge::parallel::g_rank==0&&edge::parallel::g_thread==0,INFO)
#define EDGE_LOG_TRACE            LOG(TRACE)
#define EDGE_LOG_DEBUG            LOG(DEBUG)
#define EDGE_LOG_FATAL            LOG(FATAL)
#define EDGE_LOG_ERROR            LOG(ERROR)
#define EDGE_LOG_WARNING          LOG(WARNING)
#define EDGE_LOG_VERBOSE          LOG(VERBOSE)
#define EDGE_VLOG_IS_ON(str)      VLOG_IS_ON(str)
#define EDGE_VLOG_ALL(str)        VLOG(str)
#define EDGE_VLOG(str)            VLOG_IF(edge::parallel::g_rank==0,str)

#define EDGE_CHECK(str)           CHECK(str)
#define EDGE_CHECK_EQ(str1, str2) CHECK_EQ(str1, str2)
#define EDGE_CHECK_NE(str1, str2) CHECK_NE(str1, str2)
#define EDGE_CHECK_LT(str1, str2) CHECK_LT(str1, str2)
#define EDGE_CHECK_GT(str1, str2) CHECK_GT(str1, str2)
#define EDGE_CHECK_LE(str1, str2) CHECK_LE(str1, str2)
#define EDGE_CHECK_GE(str1, str2) CHECK_GE(str1, str2)

#else

#include <cassert>
#include <iostream>

#define EDGE_LOG_INFO_ALL    std::cout << ""
#define EDGE_LOG_INFO        std::cout << ""
#define EDGE_LOG_TRACE       std::cout << ""
#define EDGE_LOG_DEBUG       std::cout << ""
#define EDGE_LOG_FATAL       std::cerr << ""
#define EDGE_LOG_ERROR       std::cerr << ""
#define EDGE_LOG_WARNING     std::cerr << ""
#define EDGE_LOG_VERBOSE     std::cout << ""
#define EDGE_VLOG_IS_ON(str) std::cout << ""
#define EDGE_VLOG_ALL(str)   std::cout << ""
#define EDGE_VLOG(str)       std::cout << ""

namespace edge {
  namespace io {
    class NullStream {
      public:
        NullStream() {}

        template<typename T>
        NullStream& operator<<( T const& ) {
          return *this;
        }
    };
  }
}

#define EDGE_CHECK(str)           assert( str          ); edge::io::NullStream() << ""
#define EDGE_CHECK_EQ(str1, str2) assert( str1 == str2 ); edge::io::NullStream() << ""
#define EDGE_CHECK_NE(str1, str2) assert( str1 != str2 ); edge::io::NullStream() << ""
#define EDGE_CHECK_LT(str1, str2) assert( str1 <  str2 ); edge::io::NullStream() << ""
#define EDGE_CHECK_GT(str1, str2) assert( str1 > str2  ); edge::io::NullStream() << ""
#define EDGE_CHECK_LE(str1, str2) assert( str1 <= str2 ); edge::io::NullStream() << ""
#define EDGE_CHECK_GE(str1, str2) assert( str1 >= str2 ); edge::io::NullStream() << ""
#endif

namespace edge {
  namespace io {
    namespace logging {
      /**
       * Configures the logging interface.
       **/
      void config();
    }
  }
}
#endif
