/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2016-2017, Regents of the University of California
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
 * VTK output.
 **/

#ifndef EDGE_IO_VTK_H
#define EDGE_IO_VTK_H

#include "constants.hpp"
#include <string>
#include <vector>

namespace edge {
  namespace io {
    class Vtk;
  }
}

class edge::io::Vtk {
  //private:
    // true if the vtk output is initialized (allocation and setup of reocurring arrays)
    bool   m_initialized;

    //! coords of the vertices associated with the elements.
    float  *m_coordsVe;

    //! connectivity information of elements to vertices.
    int    *m_connElVe;

    //! admissibility; TODO: stored as float since this is all our vtk-writer supports.
    float *m_adm;

    //! 1st order dofs in single precision, storage is element as ld, then quantities, then cruns (slowest dim).
    float  *m_dofs;

    //! pointers to the stride-1 element regions.
    float **m_ptrs;

    //! element type used in visit_writer lib.
    int m_visitElType;

    //! var names in the vtk output
    std::string  m_varNames[N_CRUNS*(N_QUANTITIES+1)];

    //! c-pts to vars names
    char const * m_varNamesC[N_CRUNS*(N_QUANTITIES+1)];

    /**
     * Initializes Vtk including the respective memory allocations.
     *
     * @param i_nVe number of vertices.
     * @param i_elPrint print elements.
     * @param i_veChars vertex characteristics.
     * @param i_elVe vertices adjacent to the elements.
     **/
    void init(       int_el                 i_nVe,
               const std::vector< int_el > &i_elPrint,
               const t_vertexChars         *i_veChars,
               const int_el               (*i_elVe)[C_ENT[T_SDISC.ELEMENT].N_VERTICES] );

  public:
    /**
     * Constructs a new Vtk interface.
     * Vtk uses a slightly modifided version of the visit_writer library.
     * Vtk usage internal storage matching the internal of visit_writer: Even though it is single precision,
     * the overhead >50% of the 1st order DOFs requirements. Future implementations might work with a buffer to reduce this.
     * Remark: The respective data structures are initialized in the first call of write([...]).
     *         -> Constructing a Vtk writer has almost no overhead.
     **/
    Vtk(): m_initialized(false){};

    /**
     * Destructs Vtk (including mem releases).
     **/
    ~Vtk();

    /**
     * Interface to visit_writer output.
     * If this function is called for the first time, the data structures of Vtk are allocated and initialized.
     *
     * @param i_outFile file to which the output is written.
     * @param i_binary true for binary output.
     * @param i_nVe number of vertices.
     * @param i_elPrint print elements.
     * @param i_lePrint sparse ids of limited print elements.
     * @param i_veChars vertex characteristics.
     * @param i_elVe ids of the elements' adjacent vertices.
     * @param i_dofs DOFs.
     * @param i_adm admissibility of the DG solution.
     **/
    void write( const std::string           &i_outFile,
                      bool                   i_binary,
                      int_el                 i_nVe,
                const std::vector< int_el > &i_elPrint,
                const std::vector< int_el > &i_lePrint,
                const t_vertexChars         *i_veChars,
                const int_el               (*i_elVe)[C_ENT[T_SDISC.ELEMENT].N_VERTICES],
                const real_base            (*i_dofs)[N_QUANTITIES][N_ELEMENT_MODES][N_CRUNS],
                const bool                 (*i_adm)[N_CRUNS]=nullptr );
};

#endif
