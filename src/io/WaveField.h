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
 * Wave field writer.
 **/

#ifndef EDGE_IO_WAVE_FIELD_H
#define EDGE_IO_WAVE_FIELD_H

#include "Vtk.h"

#include <string>
#include <limits>
#include "constants.hpp"
#include "data/EntityLayout.type"

namespace edge {
  namespace io {
    class WaveField;
  }
}

class edge::io::WaveField {
  //private:
    enum Type{ none,
               vtkAscii,
               vtkBinary };

    //! vtk interfaces
    Vtk m_vtk;

    //! type of the output
    Type m_type;

    //! numeber of written steps
    size_t m_writeStep;

    //! path to output file
    std::string m_outFile;

    //! vertex chars
    const t_vertexChars *m_veChars;

    //! element chars
    const t_elementChars *m_elChars;

    //! vertices adjacent to the elements
    const int_el (*m_elVe)[C_ENT[T_SDISC.ELEMENT].N_VERTICES];

    //! degrees of freedom
    const real_base (*m_dofs)[N_QUANTITIES][N_ELEMENT_MODES][N_CRUNS];

    //! number of vertices
    int_el m_nVe;

    //! print elements (unqiue owned elements)
    std::vector< int_el > m_elPrint;

    //! sparse ids of limited print element (std::limit<int_el>::max() if not a limited element)
    std::vector< int_el > m_liPrint;

  public:
    /**
     * Constructor of the DoF writer.
     *
     * @param i_type type of the output
     * @param i_outFile output file.
     * @param i_elLayout entity layout of the elements.
     * @param i_inMap index mapping
     * @param i_veChars characteristics of the vertices.
     * @param i_elChars characteristics of the elements.
     * @param i_elVe vertices adjacent to the elements.
     * @param i_dofs location of degrees of freedom, which will get written in corresponding calls.
     * @param i_spType sparse type for elements, which are printed. If numeric_limits<>::max(), all elements are printed.
     **/
    WaveField(       std::string      i_type,
                     std::string      i_outFile,
               const t_enLayout      &i_elLayout,
               const t_inMap         *i_inMap,
               const t_vertexChars   *i_veChars,
               const t_elementChars  *i_elChars,
               const int_el         (*i_elVe)[C_ENT[T_SDISC.ELEMENT].N_VERTICES],
               const real_base      (*i_dofs)[N_QUANTITIES][N_ELEMENT_MODES][N_CRUNS],
                     int_spType       i_spType = std::numeric_limits< int_spType >::max() );

    /**
     * Writes the given dofs.
     *
     * @param i_time time of this snapshot
     * @param i_limSync optional number of times the elements were limited since the last sync.
     **/
    void write( double         i_time,
                unsigned int (*i_limSync)[N_CRUNS] = nullptr );
};

#endif
