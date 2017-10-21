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
 * XML mesh operators.
 **/
#ifndef EDGE_MESH_XML_HPP
#define EDGE_MESH_XML_HPP

#include <vtkDataSet.h>
#include <vtkSmartPointer.h>
#include <vtksys/SystemTools.hxx>
#include <vtkXMLReader.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLStructuredGridReader.h>
#include <vtkXMLRectilinearGridReader.h>
#include <vtkXMLHyperOctreeReader.h>
#include <vtkXMLCompositeDataReader.h>
#include <vtkXMLStructuredGridReader.h>
#include <vtkXMLImageDataReader.h>
#include <vtkDataSetReader.h>

#include <vtkUnstructuredGrid.h>
#include <vtkRectilinearGrid.h>
#include <vtkHyperOctree.h>
#include <vtkImageData.h>
#include <vtkPolyData.h>
#include <vtkStructuredGrid.h>

namespace edge {
  namespace mesh_probe {
    namespace io {
      class MeshXml;
    }
  }
}

class edge::mesh_probe::io::MeshXml {
  private:
    /**
     * Reads a given XML mesh file.
     *
     * @param i_path path to XML mesh file.
     * @return vtk data set.
     *
     * @paramt TL_T_READER reader which is called.
     **/
    template< class TL_T_READER >
    static vtkDataSet* read( std::string &i_path ) {
      // create reader
      vtkSmartPointer< TL_T_READER > l_read = vtkSmartPointer< TL_T_READER >::New();

      // read file
      l_read->SetFileName( i_path.c_str() );
      l_read->Update();
      l_read->GetOutput()->Register( l_read );

      return vtkDataSet::SafeDownCast( l_read->GetOutput() );
    }

  public:
    /**
     * Reads a given XML mesh file.
     *
     * @param i_path path to XML mesh file.
     * @return vtk data set.
     **/
    static vtkDataSet* read( std::string &i_path ) {
      // get extension file extension
      std::string l_ext = vtksys::SystemTools::GetFilenameLastExtension( i_path );

      // pointer to data set
      vtkDataSet* l_dataSet = nullptr;

      // call appropiate reader based on file extension
      if(      l_ext == ".vtu" ) l_dataSet = read< vtkXMLUnstructuredGridReader >( i_path );
      else if( l_ext == ".vtp" ) l_dataSet = read< vtkXMLPolyDataReader >( i_path );
      else if( l_ext == ".vts" ) l_dataSet = read< vtkXMLStructuredGridReader >( i_path );
      else if( l_ext == ".vtr" ) l_dataSet = read< vtkXMLRectilinearGridReader >( i_path );
      else if( l_ext == ".vti" ) l_dataSet = read< vtkXMLImageDataReader >( i_path );
      else if( l_ext == ".vto" ) l_dataSet = read< vtkXMLHyperOctreeReader >( i_path );
      else if( l_ext == ".vtk" ) l_dataSet = read< vtkDataSetReader >( i_path );
      else std::cerr << "Unknown file type" << std::endl;

      return l_dataSet;
    }
};

#endif
