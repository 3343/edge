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
 * Mesh stats.
 **/
#ifndef EDGE_STATS_HPP
#define EDGE_STATS_HPP

#include <vtkCellTypes.h>
#include <vtkExtractEdges.h>
#include <vtkLine.h>

namespace edge {
  namespace mesh_probe {
    namespace info {
      class Stats;
    }
  }
}

class edge::mesh_probe::info::Stats {
  private:
    //! number of edges
    vtkIdType m_nEdges;

    //! number of cells
    vtkIdType m_nCells;

    //! edge lengths (min, ave, max)
    double m_edgeLengths[3];

    //! entity types in the mesh
    std::map< vtkIdType, vtkIdType > m_types;

    /**
     * Inits the stats.
     *
     * @param i_ds data set.
     **/
    void init( vtkDataSet* i_ds ) {
      m_nCells = i_ds->GetNumberOfCells();

      // iterate over cells
      for( vtkIdType l_ce = 0; l_ce < m_nCells; l_ce++) {
        // count type
        m_types[ i_ds->GetCellType(l_ce) ]++;
      }

      // get edges
      vtkSmartPointer< vtkExtractEdges > l_edges = vtkSmartPointer< vtkExtractEdges >::New();
      l_edges->SetInputData( i_ds );
      l_edges->Update();

      m_nEdges = l_edges->GetOutput()->GetLines()->GetNumberOfCells();

      // init edge lengths
      m_edgeLengths[0] = std::numeric_limits< double >::max();
      m_edgeLengths[1] = 0;
      m_edgeLengths[2] = std::numeric_limits< double >::lowest();

      // iterate over edges
      for( vtkIdType l_ed = 0; l_ed < m_nEdges; l_ed++ ) {
        vtkSmartPointer< vtkLine > l_edge = vtkLine::SafeDownCast( l_edges->GetOutput()->GetCell(l_ed) );

        double l_length = std::sqrt( l_edge->GetLength2() );

        m_edgeLengths[0]  = std::min( m_edgeLengths[0], l_length );
        m_edgeLengths[1] += l_length / l_edges->GetOutput()->GetLines()->GetNumberOfCells();
        m_edgeLengths[2]  = std::max( m_edgeLengths[2], l_length );
      }
    }

  public:
    /**
     * Constructur parsing the data set for stats.
     *
     * @param i_ds data set to parse.
     **/
    Stats( vtkDataSet* i_ds ) {
      init( i_ds );
    }

    /**
     * Prints the number of entities for each type.
     **/
    void printTypes() {
      std::map< vtkIdType, vtkIdType >::const_iterator l_it = m_types.begin();

      while( l_it != m_types.end() ) {
        std::cout << vtkCellTypes::GetClassNameFromTypeId(l_it->first)
                  << "," << l_it->second << std::endl;
        l_it++;
      }
      std::cout << "edges," << m_nEdges << std::endl;
    }

    /**
     * Prints statistics (min,ave,max) for edge lengths.
     **/
    void printEdgeLengths() {
      std::cout << "edge_length_min," << m_edgeLengths[0] << std::endl;
      std::cout << "edge_length_ave," << m_edgeLengths[1] << std::endl;
      std::cout << "edge_length_max," << m_edgeLengths[2] << std::endl;
    }
};

#endif
