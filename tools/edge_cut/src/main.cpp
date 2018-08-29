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
// #define CGAL_MESH_3_VERBOSE 1             // Show progress of meshing routine, print extra mesh quality statistics from optimizer
// #define CGAL_MESH_3_PROTECTION_DEBUG 1    // Very verbose output of routine preserving 1D features

#include "io/logging.hpp"
INITIALIZE_EASYLOGGINGPP
#include "surf/meshUtils.h"
#include "surf/Topo.h"
#include "surf/BdryTrimmer.h"
// #include "implicit_functions.h"


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

  const double l_xMin =  -24950;
  const double l_xMax =   24950;
  const double l_yMin =  -24950;
  const double l_yMax =   24950;
  const double l_zMin =  -10000;
  const double l_zMax =   80000;

  EDGE_LOG_INFO << "Constructing Delaunay triangulation for topography...";
  const std::string l_topoFile = "./data/topo_alex-100.xyz";
  edge_cut::surf::Topo l_topo( l_topoFile, l_xMin, l_xMax, l_yMin, l_yMax, l_zMin, l_zMax );

  EDGE_LOG_INFO << "Constructing polyhedral surface model of topography...";
  std::stringstream   topoStream;
  Polyhedron          topoSurface;
  l_topo.writeTriaToOff( topoStream );
  if (!topoStream || !(topoStream >> topoSurface) || topoSurface.is_empty()
             || !CGAL::is_triangle_mesh( topoSurface ) ) {
    std::cerr << "Input stream for topography triangulation is not valid." << std::endl;
    return 1;
  }
  topoStream.str( std::string() );

  EDGE_LOG_INFO << "Constructing polyhedral surface model of domain boundary...";
  Polyhedron l_bdryPoly;
  edge_cut::surf::makeFreeSurfBdry( l_bdryPoly, l_topo );

  EDGE_LOG_INFO << "Computing topography-boundary intersection...";
  Poly_slicer   topoSlicer( topoSurface );
  Polyline_type faultRidges, yMinRidges, yMaxRidges, xMinRidges, xMaxRidges;
  std::list< Polyline_type > features;
  K::Plane_3    yMinPlane  = K::Plane_3( 0, 1, 0, -1 * l_yMin );
  K::Plane_3    yMaxPlane  = K::Plane_3( 0, 1, 0, -1 * l_yMax );
  K::Plane_3    xMinPlane  = K::Plane_3( 1, 0, 0, -1 * l_xMin );
  K::Plane_3    xMaxPlane  = K::Plane_3( 1, 0, 0, -1 * l_xMax );

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


  EDGE_LOG_INFO << "Re-meshing polyhedral surfaces according to provided criteria";
  // The "re-meshing" form of make-mesh only works when the polyhedral domain
  // is constructed from a vector of polyhedral surfaces (as far as I can tell )
  std::vector< Polyhedron* > topoSurfVector(1, &topoSurface);
  std::vector< Polyhedron* > l_bdryVector = { &l_bdryPoly };

  // Create a polyhedral domain, with only one polyhedron,
  // and no "bounding polyhedron", so the volumetric part of the domain will be
  // empty.
  Mesh_domain l_domain( topoSurfVector.begin(), topoSurfVector.end() );
  Mesh_domain l_bdryDomain( l_bdryVector.begin(), l_bdryVector.end() );

  // Add features computed above
  l_bdryDomain.add_features( features.begin(), features.end() );
  l_bdryDomain.detect_features();
  l_domain.add_features( features.begin(), features.end() );
  l_domain.detect_features(); //includes detection of borders

  // Meshing Criteria
  //    Edge Size - Max distance between protecting balls in 1D feature preservation
  //    Facet Size - Radius of Surface Delaunay ball
  //    Facet Distance - Max distance between center for Surface Delaunay Ball and circumcenter of triangle face (surface approximation parameter)
  //    Facet Angle - Max angle of a surface triangle face [ MUST BE <= 30 ]
  //    Cell Size - Tetrahedral circumradius
  //    Cell Radius/Edge Ratio - Max ratio of circumradius to shortest edge (mesh quality parameter) [ MUST BE >= 2 ]
  //
  // NOTE the mesher cannot detect and remove "slivers" - optimizers must be
  //      used to remove them and improve mesh quality
  // NOTE A "Surface Delaunay Ball" is the 3D ball circumscribing a triangle surface facet which also
  //      has its center on the theoretical surface to be meshed. The center of the Surface
  //      Delaunay Ball will not coincide with the triangle circumcenter when the surface mesh
  //      is a poor approximation to the theoretical surface
  K::FT l_scale = 10;
  K::Point_3 l_center = CGAL::ORIGIN;
  K::FT l_innerRefineRad = 10000;
  K::FT l_outerRefineRad = 17500;

  K::FT l_edgeLengthBase = 100;
  K::FT l_facetSizeBase = 60;
  K::FT l_facetApproxBase = 80;
  K::FT l_angleBound = 25;


  // NOTE meshes must share common edge refinement criteria in order for borders
  //      to coincide. (recall edge criteria only affects specified 1D features)
  SizingField l_edgeCrit( l_edgeLengthBase, l_scale, l_center, l_innerRefineRad, l_outerRefineRad );

  SizingField l_topoFacetCrit( l_facetSizeBase, l_scale, l_center, l_innerRefineRad, l_outerRefineRad );
  SizingField l_bdryFacetCrit( l_facetSizeBase, l_scale, l_center, 0, 0 );
  SizingField l_topoApproxCrit( l_facetApproxBase, l_scale, l_center, l_innerRefineRad, l_outerRefineRad );
  SizingField l_bdryApproxCrit( l_facetApproxBase, l_scale, l_center, 0, 0 );

  Mesh_criteria   l_topoCriteria( CGAL::parameters::edge_size = l_edgeCrit,
                                  CGAL::parameters::facet_size = l_topoFacetCrit,
                                  CGAL::parameters::facet_distance = l_topoApproxCrit,
                                  CGAL::parameters::facet_angle = l_angleBound         );
  Mesh_criteria   l_bdryCriteria( CGAL::parameters::edge_size = l_edgeCrit,
                                  CGAL::parameters::facet_size = l_bdryFacetCrit,
                                  CGAL::parameters::facet_distance = l_bdryApproxCrit,
                                  CGAL::parameters::facet_angle = l_angleBound         );

  // Mesh generation
  // NOTE the optimizers have more options than specified here
  C3t3 topoComplex = CGAL::make_mesh_3<C3t3>(  l_domain, l_topoCriteria,
                CGAL::parameters::lloyd(    CGAL::parameters::time_limit = 60 ),
                CGAL::parameters::odt(      CGAL::parameters::time_limit = 60 ),
                CGAL::parameters::perturb(  CGAL::parameters::time_limit = 60 ),
                CGAL::parameters::exude(    CGAL::parameters::time_limit = 60 ) );
  C3t3 bdryComplex = CGAL::make_mesh_3<C3t3>( l_bdryDomain, l_bdryCriteria,
                CGAL::parameters::lloyd(    CGAL::parameters::time_limit = 60 ),
                CGAL::parameters::odt(      CGAL::parameters::time_limit = 60 ),
                CGAL::parameters::perturb(  CGAL::parameters::time_limit = 60 ),
                CGAL::parameters::exude(    CGAL::parameters::time_limit = 60 ) );


  EDGE_LOG_INFO << "Trimming boundary mesh...";
  Polyhedron l_bdryPolyTrimmed, topoPoly;
  edge_cut::surf::c3t3ToPolyhedron( topoComplex, topoPoly );
  edge_cut::surf::c3t3ToPolyhedron( bdryComplex, l_bdryPolyTrimmed );
  edge_cut::surf::BdryTrimmer< Polyhedron > l_trimmer( l_bdryPolyTrimmed, topoPoly );

  if ( !l_trimmer.trim() ) {
    EDGE_LOG_INFO << "Error: Unable to trim boundary mesh.";
    return 1;
  }

  // edge_cut::surf::trimFSB( topoPoly, l_bdryPolyTrimmed );



  // Output
  EDGE_LOG_INFO << "Writing meshes...";
  std::ofstream off_file("o_parkfieldTopo.off");
  std::ofstream bdry_file("o_parkfieldBdry.off");
  bdry_file << l_bdryPolyTrimmed;
  topoComplex.output_facets_in_complex_to_off( off_file );
  EDGE_LOG_INFO << "Done.";

  EDGE_LOG_INFO << "Thank you for using EDGEcut!";
}
