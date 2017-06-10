/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2016, Regents of the University of California
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
 * Interactions with the file system.
 **/

#ifndef FILE_SYSTEM_HPP_
#define FILE_SYSTEM_HPP_

#include "parallel/Mpi.h"
#include <sys/stat.h>
#include "io/logging.h"

namespace edge {
  namespace io {
    class FileSystem;
  }
}

class edge::io::FileSystem {
  public:
    /**
     * Splits the given path into the file name and directory.
     *
     * @param i_path path which gets split.
     * @param o_dir will be set to directory.
     * @param o_file will be set to file.
     **/
    static void splitPath( const std::string &i_path,
                                 std::string &o_dir,
                                 std::string &o_file ) {
      std::size_t l_split = i_path.find_last_of("/\\");
      o_dir  = i_path.substr( 0, l_split );
      o_file = i_path.substr( l_split+1  );
    }

    /**
     * Creates the given directory if not existing already.
     *
     * @param i_dir directory which is created.
     **/
    static void createDir( const std::string &i_dir ) {
      std::string l_dir, l_file;

      // get base directory
      splitPath( i_dir, l_dir, l_file );

      struct stat l_stat;
      //  check if directory exists
      if( stat( l_dir.c_str(), &l_stat ) != 0 ) {
        EDGE_LOG_INFO << "  can't access your base dir: " << l_dir;

        // try to continue by creating the directory
        if( parallel::g_rank == 0 ) {
          int l_error = mkdir( l_dir.c_str(), S_IRWXU | S_IRWXG );
          if( l_error == 0 ) {
            EDGE_LOG_INFO << "  we are still in the game, sucessfully created: " << l_dir;
          }
          else {
            EDGE_LOG_FATAL << "  creating " << l_dir << " failed: " << std::strerror( l_error );
          }
        }
      }
#ifdef PP_USE_MPI
      // wait for all MPI-ranks
      MPI_Barrier( MPI_COMM_WORLD );
#endif
    }
};
#endif
