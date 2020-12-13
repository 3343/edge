/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (alex.breuer AT uni-jena.de)
 *
 * @section LICENSE
 * Copyright (c) 2020, Friedrich Schiller University Jena
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
 * Grid holding input data.
 **/
#ifndef EDGE_V_IO_GRID_H
#define EDGE_V_IO_GRID_H

#include "Hdf5.h"

namespace edge_v {
  namespace io {
    class Grid;
  }
}

/**
 * Grid hoding input data.
 **/
class edge_v::io::Grid {
  private:
    //! data
    float * m_data = nullptr;

    //! hdf5-reader
    io::Hdf5 const * m_reader = nullptr;

    /**
     * Frees the memory.
     **/
    void free();

    /**
     * Inits the data-array at the given points from the given reader.
     *
     * @param i_nPts number of points.
     * @param i_pts coordinates of the points. The third coordinate of every point is ignored.
     * @param i_reader reader which is used to get the data.
     * @param o_data data-array which will be initialized.
     **/
    static void init( t_idx             i_nPts,
                      double   const (* i_pts)[3],
                      io::Hdf5 const  * i_reader,
                      float           * o_data );

  public:
    /**
     * Constructor.
     *
     * @param i_reader hdf5 reader.
     */
    Grid( io::Hdf5 const * i_reader );

    /**
     * Destructor.
     **/
    ~Grid();

    /**
     * Inits the internal data structures at the given points.
     *
     * @param i_nPts number of points.
     * @param i_pts coordinates of the points. The third coordinate of every point is ignored.
     **/
    void init( t_idx           i_nPts,
               double const (* i_pts)[3] );

    /**
     * Gets the data.
     *
     * @return data.
     **/
    float const * getData() const { return m_data; }
};

#endif