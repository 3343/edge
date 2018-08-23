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


  // Compute the 1D Features we want to preserve
  // This removes "rounded edges" at the boundary of the domain and creates a
  // shart fault/topography interface
  std::stringstream   topoStream;
  Polyhedron          topoSurface;
  g_topo.writeTriaToOff( topoStream );
  if (!topoStream || !(topoStream >> topoSurface) || topoSurface.is_empty()
             || !CGAL::is_triangle_mesh(topoSurface)) {
    std::cerr << "Input stream for topography triangulation is not valid." << std::endl;
    return 1;
  }
  topoStream.str( std::string() );
  topoStream.clear();
  std::cout << "Is empty: " << topoSurface.is_empty() << std::endl;
  std::cout << "Is closed: " << topoSurface.is_closed() << std::endl;
  std::cout << "Is pure triangle: " << topoSurface.is_pure_triangle() << std::endl;
  std::cout << "Num verts: " << topoSurface.size_of_vertices() << std::endl;
  std::cout << "Num edges: " << topoSurface.size_of_halfedges() << std::endl;

  Poly_slicer   topoSlicer( topoSurface );
  Polyline_type faultRidges, yMinRidges, yMaxRidges, xMinRidges, xMaxRidges;
  std::list< Polyline_type > features;
  K::Plane_3    yMinPlane  = K::Plane_3( 0, 1, 0, -1 * Y_MIN );
  K::Plane_3    yMaxPlane  = K::Plane_3( 0, 1, 0, -1 * Y_MAX );
  K::Plane_3    xMinPlane  = K::Plane_3( 1, 0, 0, -1 * X_MIN );
  K::Plane_3    xMaxPlane  = K::Plane_3( 1, 0, 0, -1 * X_MAX );

  yMinRidges  = edge_cut::surf::topoIntersect( topoSlicer, yMinPlane );
  yMaxRidges  = edge_cut::surf::topoIntersect( topoSlicer, yMaxPlane );
  xMinRidges  = edge_cut::surf::topoIntersect( topoSlicer, xMinPlane );
  xMaxRidges  = edge_cut::surf::topoIntersect( topoSlicer, xMaxPlane );
  edge_cut::surf::orderPolyline( yMinRidges, 0 );
  edge_cut::surf::orderPolyline( yMaxRidges, 0 );
  edge_cut::surf::orderPolyline( xMinRidges, 1 );
  edge_cut::surf::orderPolyline( xMaxRidges, 1 );
  features.push_back( yMinRidges );
  features.push_back( yMaxRidges );
  features.push_back( xMinRidges );
  features.push_back( xMaxRidges );


  Polyhedron l_freeSurfBdry;
  edge_cut::surf::makeFreeSurfBdry( l_freeSurfBdry );

  // The "re-meshing" form of make-mesh only works when the polyhedral domain
  // is constructed from a vector of polyhedral surfaces (as far as I can tell )
  std::vector< Polyhedron* > topoSurfVector(1, &topoSurface);
  std::vector< Polyhedron* > fsbVector = { &l_freeSurfBdry };


  // Create a polyhedral domain, with only one polyhedron,
  // and no "bounding polyhedron", so the volumetric part of the domain will be
  // empty.
  Mesh_domain l_domain(topoSurfVector.begin(), topoSurfVector.end());
  Mesh_domain fsbDomain( fsbVector.begin(), fsbVector.end() );

  // Add features computed above
  fsbDomain.add_features( features.begin(), features.end() );
  l_domain.add_features( features.begin(), features.end() );
  // l_domain.detect_features(); //includes detection of borders

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
  Mesh_criteria l_criteria( CGAL::parameters::edge_size = 100,
                            CGAL::parameters::facet_angle = 25,
                            CGAL::parameters::facet_size = 600,
                            CGAL::parameters::facet_distance = 15);


  // Mesh generation
  // NOTE the optimizers have more options than specified here
  C3t3 topoComplex = CGAL::make_mesh_3<C3t3>(  l_domain, l_criteria,
                CGAL::parameters::lloyd(    CGAL::parameters::time_limit = 60 ),
                CGAL::parameters::odt(      CGAL::parameters::time_limit = 60 ),
                CGAL::parameters::perturb(  CGAL::parameters::time_limit = 60 ),
                CGAL::parameters::exude(    CGAL::parameters::time_limit = 60 ) );

  C3t3 fsbComplex = CGAL::make_mesh_3<C3t3>( fsbDomain, l_criteria,
                CGAL::parameters::lloyd(    CGAL::parameters::time_limit = 60 ),
                CGAL::parameters::odt(      CGAL::parameters::time_limit = 60 ),
                CGAL::parameters::perturb(  CGAL::parameters::time_limit = 60 ),
                CGAL::parameters::exude(    CGAL::parameters::time_limit = 60 ) );

  EDGE_LOG_INFO << "trimming free surface boundary mesh";
  std::stringstream tempFsbStream, tempTopoStream;
  Polyhedron fsbPoly, topoPoly;
  fsbComplex.output_facets_in_complex_to_off( tempFsbStream );
  tempFsbStream >> fsbPoly;
  tempFsbStream.str( std::string() );
  topoComplex.output_facets_in_complex_to_off( tempTopoStream );
  std::cout << "topoComplexPoly char count = " << tempTopoStream.gcount() << std::endl;
  tempTopoStream >> topoPoly;
  tempTopoStream.str( std::string() );
  topoPoly.normalize_border();

  std::cout << topoSurface.size_of_border_halfedges() << " " << topoPoly.size_of_border_halfedges() << std::endl;
  edge_cut::surf::trimFSB( topoPoly, fsbPoly );



  // Output
  EDGE_LOG_INFO << "writing surface mesh";
  std::ofstream off_file("o_parkfieldTopo.off");
  std::ofstream bdry_file("o_parkfieldBdry.off");
  topoComplex.output_facets_in_complex_to_off( off_file );
  // fsbComplex.output_boundary_to_off( bdry_file );
  bdry_file << fsbPoly;

  EDGE_LOG_INFO << "thank you for using EDGEcut!";
}
