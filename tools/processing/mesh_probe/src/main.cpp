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
 * Entry point.
 **/

#include "io/MeshXml.hpp"
#include "info/Stats.hpp"

#include <vtkSmartPointer.h>
#include <vtkDataSet.h>
#include <vtkSphere.h>
#include <vtkBox.h>
#include <vtkExtractGeometry.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <cassert>

int main( int i_argc, char *i_argv[] ) {
  if( i_argc < 2) {
    std::cerr << "Usage: " << i_argv[0] << " MESH_FILE [--out OUT_FILE] [--EXT_TYPE EXT_SPECS EXT_SPECS ..]" << std::endl;
    std::cerr << " EXT_TYPE options:" << std::endl;
    std::cerr << "   1) --sphere center_x center_y center_z radius" << std::endl;
    std::cerr << "   2) --box min_x max_x min_y max_y min_z max_z" << std::endl;
    return EXIT_SUCCESS;
  }

  // get path as string
  std::string l_path = std::string( i_argv[1] );

  // read file
  vtkDataSet* l_mesh = edge::mesh_probe::io::MeshXml::read( l_path );

  // working data set
  vtkDataSet* l_work = nullptr;

  // init out-file
  std::string l_outFile("");

  // optional args, 0: out, 1: extract sphere, 2: extract box
  int l_optArgs[3] = {-1, -1, -1};

  if( i_argc > 2 ) {
    for( unsigned short l_ar = 2; l_ar < i_argc; l_ar++ ) {
      if( std::string( i_argv[l_ar] ) == "--out" )
        l_optArgs[0] = l_ar;
      else if( std::string( i_argv[l_ar] ) == "--sphere" )
        l_optArgs[1] = l_ar;
      else if( std::string( i_argv[l_ar] ) == "--box" )
        l_optArgs[2] = l_ar;
    }

    // extract region if requested
    if( l_optArgs[1] != -1 ) {
      assert( i_argc >= 7 );

      // first entry describing the sphere
      int l_first = l_optArgs[1]+1;

      vtkSmartPointer< vtkSphere > l_sphere = vtkSmartPointer< vtkSphere >::New();

      double l_center[3];
      l_center[0] = atof( i_argv[l_first+0] );
      l_center[1] = atof( i_argv[l_first+1] );
      l_center[2] = atof( i_argv[l_first+2] );

      double l_radius = atof( i_argv[l_first+3] );

      l_sphere->SetCenter( l_center[0], l_center[1], l_center[2] );
      l_sphere->SetRadius( l_radius ); 

      vtkSmartPointer< vtkExtractGeometry > l_ext = vtkExtractGeometry::New();
      l_ext->SetImplicitFunction( l_sphere );
      l_ext->SetInputData( l_mesh );
      l_ext->Update();

      l_work = l_ext->GetOutput();
    }
    else if( l_optArgs[2] != -1 ) {
      assert( i_argc >= 9 );

      // first entry describing the box
      int l_first = l_optArgs[2]+1;

      vtkSmartPointer< vtkBox > l_box = vtkSmartPointer< vtkBox >::New();

      double l_boxBnds[6];
      for( unsigned short l_en = 0; l_en < 6; l_en++ ) {
        l_boxBnds[l_en] = atof( i_argv[l_first+l_en] );
      }

      l_box->SetBounds( l_boxBnds );

      vtkSmartPointer< vtkExtractGeometry > l_ext = vtkExtractGeometry::New();
      l_ext->SetImplicitFunction( l_box );
      l_ext->SetInputData( l_mesh );
      l_ext->Update();

      l_work = l_ext->GetOutput();
    }
    else l_work = l_mesh;
  }
  else {
    l_work = l_mesh;
  }

  // generate stats
  edge::mesh_probe::info::Stats l_stats( l_work );

  // print stats
  l_stats.printTypes();
  l_stats.printEdgeLengths();

  // write file if required
  if( l_optArgs[0] != -1 ) {
    vtkSmartPointer< vtkXMLUnstructuredGridWriter > l_writer = vtkSmartPointer< vtkXMLUnstructuredGridWriter >::New();
    l_writer->SetFileName( i_argv[ l_optArgs[0] + 1] );
    l_writer->SetInputData( l_work );
    l_writer->Write();
  }

  // destroy data sets
  if( l_work != l_mesh ) l_mesh->Delete();
  l_work->Delete();

  return EXIT_SUCCESS;
}
