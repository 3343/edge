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

// CGAL mesher debug flags:
#define CGAL_MESH_3_PROFILING 1           // Output basic size and time duration statistics for mesher
#define CGAL_MESH_3_VERBOSE 1             // Show progress of meshing routine, print extra mesh quality statistics from optimizer
// #define CGAL_MESH_3_PROTECTION_DEBUG 1    // Very verbose output of routine preserving 1D features

#include "io/logging.hpp"
INITIALIZE_EASYLOGGINGPP
#include "surf/meshUtils.h"
#include "surf/Topo.h"
#include "implicit_functions.h"

using namespace edge_cut::surf;

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


  // Construct Implicit Domain
  Function l_hSpace1(&posHSpaceDisp);
  Function l_hSpace2(&negHSpaceDisp);
  Function l_uTopo(&belowTopo);
  Function_vector l_halfSpaces;
  l_halfSpaces.push_back( l_hSpace1 );
  l_halfSpaces.push_back( l_hSpace2 );
  l_halfSpaces.push_back( l_uTopo );

  std::vector<std::string> l_vps;
  l_vps.push_back("+--");
  l_vps.push_back("-+-");

  K::Point_3  l_bndSphereCenter( 415000, 3775000, -10000 );
  K::FT       l_bndSphereRadiusSquared = 150000.0*150000.0;
  K::Sphere_3 l_bndSphere( l_bndSphereCenter, l_bndSphereRadiusSquared );

  Mesh_domain l_domain( Function_wrapper( l_halfSpaces, l_vps ), l_bndSphere, 1e-3 );


  // Compute the 1D Features we want to preserve
  // This removes "rounded edges" at the boundary of the domain and creates a
  // shart fault/topography interface
  std::stringstream   topoStream;
  Surf_mesh           topoSurface;
  g_topo.writeTriaToOff( topoStream );
  if (!topoStream || !(topoStream >> topoSurface) || topoSurface.is_empty()
             || !CGAL::is_triangle_mesh(topoSurface)) {
    std::cerr << "Input stream for topography triangulation is not valid." << std::endl;
    return 1;
  }
  topoStream.str( std::string() );
  topoStream.clear();

  std::list< Polyline_type > features;
  features = build1DFeatures( topoSurface, X_MIN, X_MAX, Y_MIN, Y_MAX, Y0, Z_MIN );
  l_domain.add_features( features.begin(), features.end() );


  // Meshing Criteria
  //    Edge Size - Max distance between protecting balls in 1D feature preservation
  //    Facet Size - Radius of Surface Delaunay ball
  //    Facet Distance - Max distance between center for Surface Delaunay Ball and circumcenter of triangle face (surface approximation parameter)
  //    Facet Angle - Max angle of a surface triangle face [ MUST BE <= 30 ]
  //    Cell Size - Tetrahedral circumradius
  //    Cell Radius/Edge Ratio - Max ratio of circumradius to shortest edge (mesh quality parameter) [ MUST BE >= 2 ]
  //
  // WARNING the mesher cannot detect and remove "slivers" - optimizers must be
  //         used to remove them and improve mesh quality
  // NOTE A "Surface Delaunay Ball" is the 3D ball circumscribing a triangle surface facet which also
  //      has its center on the theoretical surface to be meshed. The center of the Surface
  //      Delaunay Ball will not coincide with the triangle circumcenter when the surface mesh
  //      is a poor approximation to the theoretical surface
  std::vector< K::FT > l_layers = {-500, -1000, -3000, -6000, -11000, -21000, -31000};
  std::vector< K::FT > l_scales = {1.000, 2.333, 4.222, 8.000, 8.111, 8.444, 10.000};

  K::FT l_edgeLengthBase, l_facetSizeBase, l_facetApproxBase, l_cellSizeBase;
  K::FT l_angleBound, l_reRatio;
  std::cout << "Base Edge Size: ";
  std::cin >> l_edgeLengthBase;
  std::cout << "Base Facet Size: ";
  std::cin >> l_facetSizeBase;
  std::cout << "Base Facet Distance: ";
  std::cin >> l_facetApproxBase;
  std::cout << "Facet Angle: ";
  std::cin >> l_angleBound;
  std::cout << "Base Cell Size: ";
  std::cin >> l_cellSizeBase;
  std::cout << "Radius-Edge Ratio: ";
  std::cin >> l_reRatio;

  DepthSizingField l_edgeCrit( l_edgeLengthBase, l_layers, l_scales );
  DepthSizingField l_facetCrit( l_facetSizeBase, l_layers, l_scales );
  DepthSizingField l_cellCrit( l_cellSizeBase, l_layers, l_scales );
  DepthSizingField l_approxCrit( l_facetApproxBase, l_layers, l_scales );

  Mesh_criteria   l_criteria( CGAL::parameters::edge_size = l_edgeCrit,
                              CGAL::parameters::facet_size = l_facetCrit,
                              CGAL::parameters::facet_distance = l_approxCrit,
                              CGAL::parameters::facet_angle = l_angleBound,
                              CGAL::parameters::cell_size = l_cellCrit,
                              CGAL::parameters::cell_radius_edge_ratio = l_reRatio );


  // Mesh generation
  // NOTE the optimizers have more options than specified here
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(  l_domain, l_criteria,
                CGAL::parameters::lloyd(    CGAL::parameters::time_limit = 60 ),
                CGAL::parameters::odt(      CGAL::parameters::time_limit = 60 ),
                CGAL::parameters::perturb(  CGAL::parameters::time_limit = 60 ),
                CGAL::parameters::exude(    CGAL::parameters::time_limit = 60 ) );


  // Output
  EDGE_LOG_INFO << "writing surface mesh";
  std::ofstream off_file("o_multImplDomains.off");
  std::ofstream bdry_file("o_bdry.off");
  c3t3.output_facets_in_complex_to_off( off_file );
  c3t3.output_boundary_to_off( bdry_file );

  EDGE_LOG_INFO << "thank you for using EDGEcut!";
}
