/**
 * @file This file is part of EDGE.
 *
 * @author Kirill Voronin (kirill.voronin AT intel.com)
 *
 * @section LICENSE
 * Copyright (c) 2022, Intel Corporation
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
 * Helpers for eltwise TPP (from LIBXSMM) data structures.
 **/

#ifndef EDGE_DATA_XSMM_UTILS_HPP
#define EDGE_DATA_XSMM_UTILS_HPP
 
#include <vector>
#include "constants.hpp"
#include "io/logging.h"
 
#include <libxsmm.h>

template <typename T>
inline libxsmm_datatype XsmmDtype();
template <>
inline libxsmm_datatype XsmmDtype<int64_t>() {
  return LIBXSMM_DATATYPE_I64;
}
template <>
inline libxsmm_datatype XsmmDtype<int32_t>() {
  return LIBXSMM_DATATYPE_I32;
}
template <>
inline libxsmm_datatype XsmmDtype<float>() {
  return LIBXSMM_DATATYPE_F32;
}
template <>
inline libxsmm_datatype XsmmDtype<double>() {
  return LIBXSMM_DATATYPE_F64;
}
/*
template <>
inline libxsmm_datatype XsmmDtype<bfloat16>() {
  return LIBXSMM_DATATYPE_BF16;
}
template <>
inline libxsmm_datatype XsmmDtype<half>() {
  return LIBXSMM_DATATYPE_F16;
}
*/

#endif