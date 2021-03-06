/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2020, Alexander Breuer
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
 * Logging interface.
 **/

#ifndef EDGE_CUT_LOGGING_HPP
#define EDGE_CUT_LOGGING_HPP

// enable thread-safety in omp-configs
#ifdef PP_USE_OMP
#define ELPP_THREAD_SAFE
#endif

// silence compilers for easylogging
#pragma GCC system_header
#include <easylogging++.h>

#define EDGE_CUT_LOG_INFO             LOG(INFO)
#define EDGE_CUT_LOG_TRACE            LOG(TRACE)
#define EDGE_CUT_LOG_DEBUG            LOG(DEBUG)
#define EDGE_CUT_LOG_FATAL            LOG(FATAL)
#define EDGE_CUT_LOG_ERROR            LOG(ERROR)
#define EDGE_CUT_LOG_WARNING          LOG(WARNING)
#define EDGE_CUT_LOG_VERBOSE          LOG(VERBOSE)
#define EDGE_CUT_VLOG_IS_ON(str)      VLOG_IS_ON(str)
#define EDGE_CUT_VLOG_ALL(str)        VLOG(str)
#define EDGE_CUT_VLOG(str)            VLOG(str)
#define EDGE_CUT_CHECK(str)           CHECK(str)
#define EDGE_CUT_CHECK_EQ(str1, str2) CHECK_EQ(str1, str2)
#define EDGE_CUT_CHECK_NE(str1, str2) CHECK_NE(str1, str2)
#define EDGE_CUT_CHECK_LT(str1, str2) CHECK_LT(str1, str2)
#define EDGE_CUT_CHECK_GT(str1, str2) CHECK_GT(str1, str2)
#define EDGE_CUT_CHECK_LE(str1, str2) CHECK_LE(str1, str2)
#define EDGE_CUT_CHECK_GE(str1, str2) CHECK_GE(str1, str2)

#endif
