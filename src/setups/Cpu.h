/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2018, Regents of the University of California
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
 * Sets processor options.
 **/
#ifndef EDGE_SETUPS_CPU_H
#define EDGE_SETUPS_CPU_H

namespace edge {
  namespace setups {
    class Cpu;
  }
}

/**
 * CPU options.
 **/
class edge::setups::Cpu {
  public:
    /**
     * @brief Get the flush-to-zero flag of the CPU.
     * 
     * @return true if enabled.
     * @return false if not.
     */
    static bool getFlushToZero();

    /**
     * @brief Sets the flush-to-zero flag of the CPU.
     *
     * @param i_on ftz will be set to on, if true, and to off, if false.
     */
    static void setFlushToZero( bool i_on );

    /**
     * @brief Gets the denormals-are-zero flag of the CPU.
     * 
     * @return true if enabled.
     * @return false if not.
     */
    static bool getDenormalsAreZero();

    /**
     * @brief Sets the denormals-are-zero flag of the CPU.
     *
     * @param i_on daz wil be set to on, if true, and to off if false.
     */
    static void setDenormalsAreZero( bool i_on );
};

#endif