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
 * Dynamic memory allocations, bookeeping and memory release on destruction.
 **/

#ifndef DYNAMIC_HPP
#define DYNAMIC_HPP

#include "common.hpp"

namespace edge {
  namespace data {
    class Dynamic;
  }
}

class edge::data::Dynamic {
  //! allocated memory
  std::vector< void* > m_mem;
  //! memory type of allocated memory
  std::vector< bool  > m_hbw;

  public:
   /**
    * Destructor, which frees all allocated memory.
    **/
  ~Dynamic();

    /**
     * Allocates aligned memory of the given size.
     * The memory will be released automatically on destruction of the object.
     * -> Do not free manually.
     *
     * @param i_size size in bytes.
     * @param i_alignment alignment of the base pointer.
     * @param i_hbw if true, function allocates high bandwidth memory if available.
     * @return pointer to memory.
     **/
  void* allocate( size_t i_size,
                  size_t i_alignment=64,
                  bool   i_hbw=false );
};
#endif
