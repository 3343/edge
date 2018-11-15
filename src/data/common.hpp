/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2015, Regents of the University of California
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
 * Shared functions of the data.
 **/
#ifndef DATA_COMMON_HPP
#define DATA_COMMON_HPP

#include "parallel/Mpi.h"
#include "constants.hpp"
#include <fstream>
#include <io/logging.h>

#ifdef PP_USE_MEMKIND
#include <hbwmalloc.h>
#endif

#ifdef PP_USE_NUMA
#include <numa.h>
#endif

namespace edge {
  namespace data {
    class common;
  }
}

class edge::data::common {
  private:
#ifdef PP_USE_NUMA
    /**
     * Checks if the NUMA-lib is available.
     **/
    static void checkNumaAvail() {
      int l_numaAvail = numa_available();
      EDGE_CHECK_NE( l_numaAvail, -1 );
    }

    /**
     * Gets the memory sizes on the node.
     *
     * @param o_mem will be set to size of memory in bytes for every numa node.
     **/
    template <typename T>
    static void getNumaMemSizes( std::vector<T>& o_mem ) {
      checkNumaAvail();

      // get #numa nodes
      int l_nNodes = numa_max_node()+1;
      o_mem.resize( l_nNodes );

      // get memory sizes
      for( std::size_t l_nd = 0; l_nd < o_mem.size(); l_nd++ ) {
        o_mem[l_nd] = numa_node_size( l_nd, NULL );
      }
    }

    /**
     * Gets the Preferred NUMA node.
     *
     * @return preferred NUMA node.
     **/
    static int getNumaPreferred() {
      checkNumaAvail();

      return numa_preferred();
    }
#endif

  public:
    /**
     * Prints memory statistics.
     **/
    static void printMemStats() {
      unsigned long l_mem[3] = {0,0,0}; // 0: total, 1: free: 2: available

      // open meminfo
      std::ifstream l_memInfo("/proc/meminfo");

      // parse file
      if( l_memInfo.is_open() ) {
        std::string l_line;
        while( l_memInfo >> l_line ) {
          if( l_line == "MemTotal:" )
            l_memInfo >> l_mem[0];
          else if( l_line == "MemFree:" )
            l_memInfo >> l_mem[1];
          else if( l_line == "MemAvailable:" ) {
            l_memInfo >> l_mem[2];
            break;
          }
        }

        // stats
        double l_gib = 1024 * 1024 * 1024 / 1000.0;
        double l_stats[3][3];
#ifdef PP_USE_MPI
        unsigned long l_tmp[3][3];
        int l_err;
        l_err = MPI_Allreduce( l_mem, l_tmp[0], 3, MPI_UNSIGNED_LONG, MPI_MIN, MPI_COMM_WORLD );
        EDGE_CHECK_EQ( l_err, MPI_SUCCESS );
        l_err = MPI_Allreduce( l_mem, l_tmp[1], 3, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD );
        EDGE_CHECK_EQ( l_err, MPI_SUCCESS );
        l_err = MPI_Allreduce( l_mem, l_tmp[2], 3, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD );
        EDGE_CHECK_EQ( l_err, MPI_SUCCESS );

        for( unsigned short l_ty = 0; l_ty < 3; l_ty++ ) {
          l_stats[0][l_ty] = l_tmp[0][l_ty] / l_gib;
          l_stats[1][l_ty] = l_tmp[1][l_ty] / ( parallel::g_nRanks*l_gib);
          l_stats[2][l_ty] = l_tmp[2][l_ty] / l_gib;
        }

#else
        for( unsigned short l_st = 0; l_st < 3; l_st++ )
          for( unsigned short l_ty = 0; l_ty < 3; l_ty++ )
            l_stats[l_st][l_ty] = l_mem[l_ty] / l_gib;
#endif

        // print
        EDGE_LOG_INFO << "memory statistics:";
        EDGE_LOG_INFO << "  total (min/ave/max): "
                      << l_stats[0][0] << " / " << l_stats[1][0] << " / " << l_stats[2][0] << " GiB";
        EDGE_LOG_INFO << "  free  (min/ave/max): "
                      << l_stats[0][1] << " / " << l_stats[1][1] << " / " << l_stats[2][1] << " GiB";
        EDGE_LOG_INFO << "  avail (min/ave/max): "
                      << l_stats[0][2] << " / " << l_stats[1][2] << " / " << l_stats[2][2] << " GiB";
      }
    }

    /**
     * Prints the sizes of the NUMA nodes.
     **/
    static void printNumaSizes() {
#ifdef PP_USE_NUMA
      EDGE_LOG_INFO << "numa node sizes (*preferred):";

      // get memory sizes of numa nodes
      std::vector< long > l_mem;
      getNumaMemSizes( l_mem );

      // get mem stats
      double l_gib = 1024 * 1024 * 1024;

      // print
      for( std::size_t l_nd = 0; l_nd < l_mem.size(); l_nd++ ) {
        EDGE_LOG_INFO << "  #" << l_nd << (getNumaPreferred() == (int) l_nd ? "*": "") << ": "
                      << l_mem[l_nd] / l_gib << " GiB";
      }
#endif
    }

    /**
     * Allocates aligned memory of the given size.
     *
     * @param i_size size in bytes.
     * @param i_alignment alignment of the base pointer.
     * @param i_hbw if true, function allocates high bandwidth memory if available.
     * @return pointer to memory.
     **/
    static void* allocate( size_t i_size,
                           size_t i_alignment=64,
                           bool   i_hbw=false ) {
      if( i_size > 0 ) {
        void* l_ptrBuffer;
        bool l_err;
#ifdef PP_USE_MEMKIND
        if( i_hbw ) {
          l_err = (hbw_posix_memalign( &l_ptrBuffer, i_alignment, i_size ) != 0);
        }
        else {
          l_err = (posix_memalign( &l_ptrBuffer, i_alignment, i_size ) != 0);
        }
#else
        l_err = (posix_memalign( &l_ptrBuffer, i_alignment, i_size ) != 0);
#endif
        EDGE_CHECK(!l_err) << "The malloc failed (bytes: " << i_size << ", alignment: " << i_alignment << ").";

        return l_ptrBuffer;
      }
      else return nullptr;
    }

    /**
     * Releases memory.
     *
     * @param io_memory pointer to memory.
     * @param i_hbw true if high bandwidth memory.
     **/
    static void release( void *i_memory,
                         bool  i_hbw=false ) {
#ifdef PP_USE_MEMKIND
      if( i_hbw ) hbw_free(i_memory);
      else free(i_memory);
#else
      free(i_memory);
#endif
    }
};

#endif
