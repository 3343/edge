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
 * Interactions with the file system.
 **/

#ifndef EDGE_FILE_SYSTEM_HPP_
#define EDGE_FILE_SYSTEM_HPP_

#include <sys/stat.h>
#include "io/logging.h"
#include <unistd.h>
#include <cstring>

namespace edge {
  namespace io {
    class FileSystem;
  }
}

class edge::io::FileSystem {
  public:
    /**
     * Splits the given path into the string before the first '/' and that after.
     *
     * @param i_path path which gets split.
     * @param o_dir will be set to directory.
     * @param o_file will be set to file.
     **/
    static void splitPathFirst( const std::string &i_path,
                                      std::string &o_dir,
                                      std::string &o_file ) {
      std::size_t l_split = i_path.find_first_of("/\\");
      if( l_split != std::string::npos ) {
        o_dir  = i_path.substr( 0, l_split );
        o_file = i_path.substr( l_split+1  );
      }
      else {
        o_dir = "";
        o_file = i_path;
      }
    }

    /**
     * Splits the given path into the string before the last '/' and that after.
     *
     * @param i_path path which gets split.
     * @param o_dir will be set to directory.
     * @param o_file will be set to file.
     **/
    static void splitPathLast( const std::string &i_path,
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
      std::string l_dirCreate;
      std::string l_dir0;
      std::string l_dir1;
      splitPathFirst( i_dir+"/", l_dir0, l_dir1 );

      // add slash, if this is system root
      if( l_dir0 == "" ) l_dir0 = "/";

      struct stat l_stat;

      // create until all sub-directories are there
      while( l_dir0 != "" ) {
        if( l_dir0 != "/" ) l_dir0 += "/";
        l_dirCreate += l_dir0;

        //  check if directory exists
        if( stat( l_dirCreate.c_str(), &l_stat ) != 0 ) {
          // try to continue by creating the directory
          int l_error = mkdir( l_dirCreate.c_str(), S_IRWXU | S_IRWXG );

          // abort on error (other than exist, as this might be called from multiple ranks)
          if( l_error != 0 && errno != EEXIST ) {
            EDGE_LOG_FATAL << "  creating " << l_dirCreate << " failed: " << std::strerror( errno );
          }

          // wait for the file system
          fsync( l_error );
        }

        std::string l_dirTmp = l_dir1;
        splitPathFirst( l_dirTmp, l_dir0, l_dir1 );
      }
    }
};
#endif
