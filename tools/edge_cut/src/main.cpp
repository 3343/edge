/**
 * @file This file is part of EDGE.
 *
 * @author David Lenz (dlenz AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2018, Regents of the University of California
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

// Debug options
// #include <gperftools/heap-profiler.h>      // For memory profiling
#define CGAL_MESH_3_PROFILING 1            // Output basic size and time duration statistics for mesher
#define CGAL_MESH_3_VERBOSE 1              // Show progress of meshing routine, print extra mesh quality statistics from optimizer
// #define CGAL_MESH_3_PROTECTION_DEBUG 1     // Very verbose output of routine preserving 1D features

#include "io/logging.hpp"
INITIALIZE_EASYLOGGINGPP
#include "../../../submodules/pugixml/src/pugixml.hpp"
#include "surf/meshUtils.h"
#include "surf/BdryTrimmer.h"

using namespace edge_cut::surf;

int main( int i_argc, char *i_argv[] ) {
  // Check for valid arguments
  std::string l_xmlPath = "";
  if ( i_argc > 2 ) {
    EDGE_LOG_INFO << "Encountered more than one command line option -- exiting";
    return -1;
  }  else if ( i_argc == 1 ) {
    EDGE_LOG_INFO << "Missing command line option (xml config path) -- exiting";
    return -1;
  } else {
    l_xmlPath = i_argv[1];
  }

  // Check for valid config
  pugi::xml_document doc;
  if (!doc.load_file( l_xmlPath.c_str() )) {
    EDGE_LOG_INFO << "Could not open xml config file -- exiting";
    return -1;
  }


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

  const double l_xMin = std::stod( doc.child("bbox").child_value("xMin") );
  const double l_xMax = std::stod( doc.child("bbox").child_value("xMax") );
  const double l_yMin = std::stod( doc.child("bbox").child_value("yMin") );
  const double l_yMax = std::stod( doc.child("bbox").child_value("yMax") );
  const double l_zMin = std::stod( doc.child("bbox").child_value("zMin") );
  const double l_zMax = std::stod( doc.child("bbox").child_value("zMax") );
  const std::string l_topoIn  = doc.child("io").child_value("topo_in");
  const std::string l_topoOut = doc.child("io").child_value("topo_out");
  const std::string l_bdryOut = doc.child("io").child_value("bdry_out");

  double const l_bBox[] = { l_xMin, l_xMax, l_yMin, l_yMax, l_zMin, l_zMax };

  // Create polyhedral meshes of topography (with whatever sampling we were provided),
  // and of domain boundary
  Polyhedron l_topoPoly, l_bdryPoly;
  edge_cut::surf::topoPolyMeshFromXYZ( l_topoPoly, l_topoIn );
  edge_cut::surf::makeBdry( l_bdryPoly, l_bBox );

  // Create polyhedral domains for re-meshing
  // The "re-meshing" form of make-mesh only works when the polyhedral domain
  // is constructed from a vector of polyhedral surfaces (as far as I can tell)
  std::vector< Polyhedron* > l_topoVector( 1, &l_topoPoly );
  std::vector< Polyhedron* > l_bdryVector( 1, &l_bdryPoly );
  Mesh_domain l_topoDomain( l_topoVector.begin(), l_topoVector.end() );
  Mesh_domain l_bdryDomain( l_bdryVector.begin(), l_bdryVector.end() );

  // Preserve the intersection between topography and boundary meshes.
  // This ensures that the two meshes coincide on their boundaries.
  std::list< Polyline_type > l_intersectFeatures = edge_cut::surf::getIntersectionFeatures( l_topoPoly, l_bBox );
  l_bdryDomain.add_features( l_intersectFeatures.begin(), l_intersectFeatures.end() );
  l_topoDomain.add_features( l_intersectFeatures.begin(), l_intersectFeatures.end() );

  // Free memory, a copy is stored in the Mesh_domain object (no move semantics in CGAL yet)
  l_topoPoly.clear();
  l_bdryPoly.clear();

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
  K::FT l_scale = std::stod( doc.child("region").child_value("scale") );
  K::FT l_innerRefineRad = std::stod( doc.child("region").child_value("inner_rad") );
  K::FT l_outerRefineRad = std::stod( doc.child("region").child_value("outer_rad") );
  K::Point_3 l_center = K::Point_3( std::stod( doc.child("region").child("center").child_value("x") ),
                                    std::stod( doc.child("region").child("center").child_value("x") ),
                                    std::stod( doc.child("region").child("center").child_value("x") ) );

  K::FT l_edgeLengthBase  = std::stod( doc.child("refine").child_value("edge") );
  K::FT l_facetSizeBase   = std::stod( doc.child("refine").child_value("facet") );
  K::FT l_facetApproxBase = std::stod( doc.child("refine").child_value("approx") );
  K::FT l_angleBound      = std::stod( doc.child("refine").child_value("angle") );

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
  EDGE_LOG_INFO << "Re-meshing polyhedral surfaces according to provided criteria";
  C3t3 topoComplex = CGAL::make_mesh_3<C3t3>( l_topoDomain,
                                              l_topoCriteria,
                                              CGAL::parameters::lloyd(    CGAL::parameters::time_limit = 60 ),
                                              CGAL::parameters::odt(      CGAL::parameters::time_limit = 60 ),
                                              CGAL::parameters::perturb(  CGAL::parameters::time_limit = 60 ),
                                              CGAL::parameters::exude(    CGAL::parameters::time_limit = 60 ) );

  C3t3 bdryComplex = CGAL::make_mesh_3<C3t3>( l_bdryDomain,
                                              l_bdryCriteria,
                                              CGAL::parameters::lloyd(    CGAL::parameters::time_limit = 60 ),
                                              CGAL::parameters::odt(      CGAL::parameters::time_limit = 60 ),
                                              CGAL::parameters::perturb(  CGAL::parameters::time_limit = 60 ),
                                              CGAL::parameters::exude(    CGAL::parameters::time_limit = 60 ) );

  // Trim bits of boundary mesh which extend above topography
  EDGE_LOG_INFO << "Trimming boundary mesh...";
  Polyhedron l_topoPolyMeshed, l_bdryPolyMeshed;
  edge_cut::surf::c3t3ToPolyhedron( topoComplex, l_topoPolyMeshed );
  edge_cut::surf::c3t3ToPolyhedron( bdryComplex, l_bdryPolyMeshed );
  edge_cut::surf::BdryTrimmer< Polyhedron > l_trimmer( l_bdryPolyMeshed, l_topoPolyMeshed );

  l_trimmer.trim();

  // Output
  EDGE_LOG_INFO << "Writing meshes...";
  std::ofstream topo_file( l_topoOut );
  std::ofstream bdry_file( l_bdryOut );

  bdry_file << l_bdryPolyMeshed;
  topoComplex.output_facets_in_complex_to_off( topo_file );
  EDGE_LOG_INFO << "Done.";
  EDGE_LOG_INFO << "Thank you for using EDGEcut!";
}
