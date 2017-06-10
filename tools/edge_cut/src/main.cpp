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
 * This is the main file of EDGEcut.
 **/

#include "io/logging.hpp"
INITIALIZE_EASYLOGGINGPP

#include <CGAL/Surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>


#include "surf/Oracle.h"
#include "surf/Topo.h"
#include "surf/LayeredCriteria.h"


int main( int i_argc, char *i_argv[] ) {
  
  EDGE_LOG_INFO << "##########################################################################";
  EDGE_LOG_INFO << "##############   ##############            ###############  ##############";
  EDGE_LOG_INFO << "##############   ###############         ################   ##############";
  EDGE_LOG_INFO << "#####            #####       #####      ######                       #####";
  EDGE_LOG_INFO << "#####            #####        #####    #####                         #####";
  EDGE_LOG_INFO << "#############    #####         #####  #####                  #############";
  EDGE_LOG_INFO << "#############    #####         #####  #####      #########   #############";
  EDGE_LOG_INFO << "#####            #####         #####  #####      #########           #####";
  EDGE_LOG_INFO << "#####            #####        #####    #####        ######           #####";
  EDGE_LOG_INFO << "#####            #####       #####      #####       #####            #####";
  EDGE_LOG_INFO << "###############  ###############         ###############   ###############";
  EDGE_LOG_INFO << "###############  ##############           #############    ###############";
  EDGE_LOG_INFO << "#######################################################################cut";
  EDGE_LOG_INFO << "";

  EDGE_LOG_INFO << "ready to go..";
  // delaunay triangulation
  CGAL::Surface_mesh_default_triangulation_3 l_delTria; 
  // complex
  CGAL::Complex_2_in_triangulation_3<
    CGAL::Surface_mesh_default_triangulation_3
  > l_compl(l_delTria);


  // bounding sphere
  CGAL::Surface_mesh_default_triangulation_3::Geom_traits::Point_3  l_bndSphereCenter(400000,3750000,-25000);
  CGAL::Surface_mesh_default_triangulation_3::Geom_traits::FT       l_bndSphereRadiusSquared = 250000.0*250000.0;
  CGAL::Surface_mesh_default_triangulation_3::Geom_traits::Sphere_3 l_bndSphere( l_bndSphereCenter,
                                                                                 l_bndSphereRadiusSquared);

  // surround box:     x1      x2    y1      y2      z1
  // 330000/500000/3700000/3850000
  double l_box[5] = {335000, 495000, 3705000, 3845000, -25000};
  
  // construct oracle
  EDGE_LOG_INFO << "constructing oracle for implicit surface meshing";
  edge_cut::surf::Oracle l_oracle( l_box, "data/map_proj.xyz" );

  EDGE_LOG_INFO << "implicit surface constructor";
  
  CGAL::Implicit_surface_3<
    CGAL::Surface_mesh_default_triangulation_3::Geom_traits,
    edge_cut::surf::Oracle
    > l_implSurf( l_oracle, l_bndSphere, 1e-5 ); 

  EDGE_LOG_INFO << "mesh criteria constructor";
  
  edge_cut::surf::LayeredCriteria<
    CGAL::Surface_mesh_default_triangulation_3
    > l_criteria( 30., 75., 75., {-500, -1000, -3000, -6000, -11000, -21000, -31000}, {4.000, 5.333, 6.222, 8.000, 8.111, 8.444, 10.000});

  
  // derive surface mesh
  EDGE_LOG_INFO << "deriving 3D surface mesh..";
  CGAL::make_surface_mesh( l_compl, l_implSurf, l_criteria, CGAL::Non_manifold_tag() );
  EDGE_LOG_INFO << "surface mesh stats:";
  EDGE_LOG_INFO << "  #vertices: " << l_delTria.number_of_vertices();
  EDGE_LOG_INFO << "  #faces:    " << l_delTria.number_of_facets(); // TODO: fix, this should be edges

  // write surface mesh
  EDGE_LOG_INFO << "writing surface mesh";
  std::ofstream l_out("out.off");
  CGAL::output_surface_facets_to_off( l_out, l_compl );

  EDGE_LOG_INFO << "thank you for using EDGEcut!";
}
