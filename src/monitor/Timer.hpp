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
 * Access to timing functions
 **/
#ifndef TIMER_HPP
#define TIMER_HPP

#include <sys/time.h>

namespace edge {
  namespace monitor {
    class Timer;
  }
}

class edge::monitor::Timer {
  private:
    //! start time
    double l_start;

    //! end time
    double l_end;

    /**
     * Gets the current wall clock time.
     *
     * @return wallclock time in seconds.
     **/
    double getWtime() const {
      struct timeval l_time;

      if( gettimeofday( &l_time, NULL ) ) {
        return 0;
      }

      return (double) l_time.tv_sec + (double) l_time.tv_usec * 0.000001;      
    }

  public:
    Timer(): l_start(0), l_end(0){}

    /**
     * Starts the timer.
     **/
    void start() {
      l_start = getWtime();
    }

    /**
     * Ends the time.
     **/
    void end() {
      l_end = getWtime();
    }

    /**
     * Returns the elapsed time between start and end (seconds).
     *
     * @return elapsed time.
     **/
    double elapsed() const {
      return l_end - l_start;
    }
};

#endif
