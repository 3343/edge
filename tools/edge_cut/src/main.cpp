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

// Output basic size and time duration statistics for mesher
//#define CGAL_MESH_3_PROFILING 1

// Show progress of meshing routine, print extra mesh quality statistics from optimizer
#define CGAL_MESH_3_VERBOSE 1

// Very verbose output of routine preserving 1D features
//#define CGAL_MESH_3_PROTECTION_DEBUG 1

#include "surf/meshUtils.h"
#include "surf/Topo.h"
#include "implicit_functions.h"

using namespace edge_cut::surf;
using namespace CGAL::parameters;

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
  // Get 1D features above fault plane
  // This reduces "divots" where the two halfspaces meet
  std::stringstream   topoStream;
  Surf_mesh           topoSurface;

  // insert corner points to delaunay triangulation
  // TODO check if this is still necessary for polyline consistency
  g_topo.m_delTria.insert( g_topo.interpolatePt( X_MIN, Y_MIN ) );
  g_topo.m_delTria.insert( g_topo.interpolatePt( X_MIN, Y_MAX ) );
  g_topo.m_delTria.insert( g_topo.interpolatePt( X_MAX, Y_MIN ) );
  g_topo.m_delTria.insert( g_topo.interpolatePt( X_MAX, Y_MAX ) );
  g_topo.m_delTria.insert( g_topo.interpolatePt( X_MIN, Y0 ) );
  g_topo.m_delTria.insert( g_topo.interpolatePt( X_MAX, Y0 ) );


  g_topo.writeTriaToOff( topoStream );
  if (!topoStream || !(topoStream >> topoSurface) || topoSurface.is_empty()
             || !CGAL::is_triangle_mesh(topoSurface)) {
    std::cerr << "Input stream for topography triangulation is not valid." << std::endl;
    return 1;
  }
  topoStream.str( std::string() );
  topoStream.clear();

  Polyline_type faultRidges, yMinRidges, yMaxRidges, xMinRidges, xMaxRidges;
  Polyline_type xMinRidges0, xMinRidges1, xMaxRidges0, xMaxRidges1;
  std::list< Polyline_type > features;
  Poly_slicer   topoSlicer( topoSurface );
  K::Plane_3    faultPlane = K::Plane_3( 0, 1, 0, -1 * Y0 );
  K::Plane_3    yMinPlane  = K::Plane_3( 0, 1, 0, -1 * Y_MIN );
  K::Plane_3    yMaxPlane  = K::Plane_3( 0, 1, 0, -1 * Y_MAX );
  K::Plane_3    xMinPlane  = K::Plane_3( 1, 0, 0, -1 * X_MIN );
  K::Plane_3    xMaxPlane  = K::Plane_3( 1, 0, 0, -1 * X_MAX );


  // Intersect fault plane with topography
  faultRidges = topoIntersect( topoSlicer, faultPlane );
  yMinRidges  = topoIntersect( topoSlicer, yMinPlane );
  yMaxRidges  = topoIntersect( topoSlicer, yMaxPlane );
  xMinRidges  = topoIntersect( topoSlicer, xMinPlane );
  xMaxRidges  = topoIntersect( topoSlicer, xMaxPlane );
  orderPolyline( faultRidges, 0 );
  orderPolyline( yMinRidges, 0 );
  orderPolyline( yMaxRidges, 0 );
  orderPolyline( xMinRidges, 1 );
  orderPolyline( xMaxRidges, 1 );


  // Trim polylines to computational domain and add corner vertices.
  // This is to avoid the case where there already exists a vertex close to one
  // of the corners of the domain. It is not clear what an acceptable tolerance
  // is so tol=10 has been chosen arbitrarily.
  double const tol = 10.;
  Polyline_type l_tempLine;
  l_tempLine.clear();
  for ( std::size_t idx = 0; idx < yMinRidges.size(); idx++ ) {
    if ( yMinRidges[ idx ].x() > X_MIN + tol && yMinRidges[ idx ].x() < X_MAX - tol ) {
      l_tempLine.push_back( yMinRidges[ idx ] );
    }
  }
  l_tempLine.insert( l_tempLine.begin(), g_topo.interpolatePt( X_MIN, Y_MIN ) );
  l_tempLine.push_back( g_topo.interpolatePt( X_MAX, Y_MIN ) );
  for ( std::size_t idx = 1; idx < l_tempLine.size(); idx++ ) {
    assert( l_tempLine[ idx ].x() > l_tempLine[ idx-1 ].x() );
  }
  yMinRidges = l_tempLine;

  l_tempLine.clear();
  for ( std::size_t idx = 0; idx < faultRidges.size(); idx++ ) {
    if ( faultRidges[ idx ].x() > X_MIN + tol && faultRidges[ idx ].x() < X_MAX - tol )
      l_tempLine.push_back( faultRidges[ idx ] );
  }
  l_tempLine.insert( l_tempLine.begin(), g_topo.interpolatePt( X_MIN, Y0 ) );
  l_tempLine.push_back( g_topo.interpolatePt( X_MAX, Y0 ) );
  for ( std::size_t idx = 1; idx < l_tempLine.size(); idx++ ) {
    assert( l_tempLine[ idx ].x() > l_tempLine[ idx-1 ].x() );
  }
  faultRidges = l_tempLine;

  l_tempLine.clear();
  for ( std::size_t idx = 0; idx < yMaxRidges.size(); idx++ ) {
    if ( yMaxRidges[ idx ].x() > X_MIN + tol && yMaxRidges[ idx ].x() < X_MAX - tol )
      l_tempLine.push_back( yMaxRidges[ idx ] );
  }
  l_tempLine.insert( l_tempLine.begin(), g_topo.interpolatePt( X_MIN, Y_MAX ) );
  l_tempLine.push_back( g_topo.interpolatePt( X_MAX, Y_MAX ) );
  for ( std::size_t idx = 1; idx < l_tempLine.size(); idx++ ) {
    assert( l_tempLine[ idx ].x() > l_tempLine[ idx-1 ].x() );
  }
  yMaxRidges = l_tempLine;


  xMinRidges0.clear();
  xMinRidges1.clear();
  for ( std::size_t idx = 0; idx < xMinRidges.size(); idx++ ) {
    if ( xMinRidges[ idx ].y() > Y_MIN + tol && xMinRidges[ idx ].y() < Y0 - tol )
      xMinRidges0.push_back( xMinRidges[ idx ] );
    if ( xMinRidges[ idx ].y() > Y0 + tol && xMinRidges[ idx ].y() < Y_MAX - tol )
      xMinRidges1.push_back( xMinRidges[ idx ] );
  }
  xMinRidges0.insert( xMinRidges0.begin(), g_topo.interpolatePt( X_MIN, Y_MIN ) );
  xMinRidges0.push_back( g_topo.interpolatePt( X_MIN, Y0 ) );
  xMinRidges1.insert( xMinRidges1.begin(), g_topo.interpolatePt( X_MIN, Y0 ) );
  xMinRidges1.push_back( g_topo.interpolatePt( X_MIN, Y_MAX ) );

  xMaxRidges0.clear();
  xMaxRidges1.clear();
  for ( std::size_t idx = 0; idx < xMaxRidges.size(); idx++ ) {
    if ( xMaxRidges[ idx ].y() > Y_MIN + tol && xMaxRidges[ idx ].y() < Y0 - tol )
      xMaxRidges0.push_back( xMaxRidges[ idx ] );
    if ( xMaxRidges[ idx ].y() > Y0 + tol && xMaxRidges[ idx ].y() < Y_MAX - tol )
      xMaxRidges1.push_back( xMaxRidges[ idx ] );
  }
  xMaxRidges0.insert( xMaxRidges0.begin(), g_topo.interpolatePt( X_MAX, Y_MIN ) );
  xMaxRidges0.push_back( g_topo.interpolatePt( X_MAX, Y0 ) );
  xMaxRidges1.insert( xMaxRidges1.begin(), g_topo.interpolatePt( X_MAX, Y0 ) );
  xMaxRidges1.push_back( g_topo.interpolatePt( X_MAX, Y_MAX ) );

  features.push_back( faultRidges );
  features.push_back( yMinRidges );
  features.push_back( yMaxRidges );
  features.push_back( xMinRidges0 );
  features.push_back( xMinRidges1 );
  features.push_back( xMaxRidges0 );
  features.push_back( xMaxRidges1 );

  // Vertical Edges
  features.push_back( { K::Point_3( yMinRidges.front().x(),
                                    yMinRidges.front().y(),
                                    yMinRidges.front().z()  ),
                        K::Point_3( yMinRidges.front().x(),
                                    yMinRidges.front().y(),
                                    Z_MIN                   ) } );
  features.push_back( { K::Point_3( yMinRidges.back().x(),
                                    yMinRidges.back().y(),
                                    yMinRidges.back().z()   ),
                        K::Point_3( yMinRidges.back().x(),
                                    yMinRidges.back().y(),
                                    Z_MIN                   ) } );
  features.push_back( { K::Point_3( yMaxRidges.front().x(),
                                    yMaxRidges.front().y(),
                                    yMaxRidges.front().z()  ),
                        K::Point_3( yMaxRidges.front().x(),
                                    yMaxRidges.front().y(),
                                    Z_MIN                   ) } );
  features.push_back( { K::Point_3( yMaxRidges.back().x(),
                                    yMaxRidges.back().y(),
                                    yMaxRidges.back().z()   ),
                        K::Point_3( yMaxRidges.back().x(),
                                    yMaxRidges.back().y(),
                                    Z_MIN                   ) } );
  features.push_back( { faultRidges.front(),
                        K::Point_3( faultRidges.front().x(),
                                    faultRidges.front().y(),
                                    Z_MIN                   ) } );
  features.push_back( { faultRidges.back(),
                        K::Point_3( faultRidges.back().x(),
                                    faultRidges.back().y(),
                                    Z_MIN                   ) } );

  // Bottom of computational domain
  features.push_back( { K::Point_3( yMinRidges.front().x(),
                                    yMinRidges.front().y(),
                                    Z_MIN                   ),
                        K::Point_3( faultRidges.front().x(),
                                    faultRidges.front().y(),
                                    Z_MIN                   ) } );
  features.push_back( { K::Point_3( faultRidges.front().x(),
                                    faultRidges.front().y(),
                                    Z_MIN                   ),
                        K::Point_3( yMaxRidges.front().x(),
                                    yMaxRidges.front().y(),
                                    Z_MIN                   ) } );
  features.push_back( { K::Point_3( yMaxRidges.front().x(),
                                    yMaxRidges.front().y(),
                                    Z_MIN                   ),
                        K::Point_3( yMaxRidges.back().x(),
                                    yMaxRidges.back().y(),
                                    Z_MIN                   ) } );
  features.push_back( { K::Point_3( yMaxRidges.back().x(),
                                    yMaxRidges.back().y(),
                                    Z_MIN                   ),
                        K::Point_3( faultRidges.back().x(),
                                    faultRidges.back().y(),
                                    Z_MIN                   ) } );
  features.push_back( { K::Point_3( faultRidges.back().x(),
                                    faultRidges.back().y(),
                                    Z_MIN                   ),
                        K::Point_3( yMinRidges.back().x(),
                                    yMinRidges.back().y(),
                                    Z_MIN                   ) } );
  features.push_back( { K::Point_3( yMinRidges.back().x(),
                                    yMinRidges.back().y(),
                                    Z_MIN                   ),
                        K::Point_3( yMinRidges.front().x(),
                                    yMinRidges.front().y(),
                                    Z_MIN                   ) } );
  features.push_back( { K::Point_3( faultRidges.front().x(),
                                    faultRidges.front().y(),
                                    Z_MIN                   ),
                        K::Point_3( faultRidges.back().x(),
                                    faultRidges.back().y(),
                                    Z_MIN                   ) } );


  // Construct Implicit Domain
  Function_vector l_halfSpaces;
  Function l_hSpace1(&posHSpaceDisp);
  Function l_hSpace2(&negHSpaceDisp);
  Function l_uTopo(&belowTopo);
  l_halfSpaces.push_back( l_hSpace1 );
  l_halfSpaces.push_back( l_hSpace2 );
  l_halfSpaces.push_back( l_uTopo );

  std::vector<std::string> vps;
  vps.push_back("+--");
  vps.push_back("-+-");

  K::Point_3  l_bndSphereCenter( 415000, 3775000, -10000 );
  K::FT       l_bndSphereRadiusSquared = 150000.0*150000.0;
  K::Sphere_3 l_bndSphere( l_bndSphereCenter, l_bndSphereRadiusSquared );

  Mesh_domain l_domain( Function_wrapper( l_halfSpaces, vps ), l_bndSphere, 1e-3 );


  // Criteria -- (300, 25, 300, 200, 3, 500) are "reasonable" to start with
  std::vector< K::FT > l_layers = {-500, -1000, -3000, -6000, -11000, -21000, -31000};
  std::vector< K::FT > l_scales = {1.000, 2.333, 4.222, 8.000, 8.111, 8.444, 10.000};

  K::FT l_edgeLengthBase, l_facetSizeBase, l_facetApproxBase, l_cellSizeBase;
  K::FT l_angleBound, l_reRatio;
  std::cout << "base edge length: ";
  std::cin >> l_edgeLengthBase;
  std::cout << "base facet size: ";
  std::cin >> l_facetSizeBase;
  std::cout << "base facet approx: ";
  std::cin >> l_facetApproxBase;
  std::cout << "base cell size: ";
  std::cin >> l_cellSizeBase;
  std::cout << "angle bound: ";
  std::cin >> l_angleBound;
  std::cout << "radius-edge ratio: ";
  std::cin >> l_reRatio;

  DepthSizingField l_edgeCrit( l_edgeLengthBase, l_layers, l_scales );
  DepthSizingField l_facetCrit( l_facetSizeBase, l_layers, l_scales );
  DepthSizingField l_cellCrit( l_cellSizeBase, l_layers, l_scales );
  DepthSizingField l_approxCrit( l_facetApproxBase, l_layers, l_scales );

  Mesh_criteria   l_criteria( CGAL::parameters::edge_size = l_edgeCrit,
                            CGAL::parameters::facet_size = l_facetCrit,
                            CGAL::parameters::cell_size = l_cellCrit,
                            CGAL::parameters::facet_distance = l_approxCrit,
                            CGAL::parameters::facet_angle = l_angleBound,
                            CGAL::parameters::cell_radius_edge_ratio = l_reRatio );

  l_domain.add_features( features.begin(), features.end() );

  // Mesh generation
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(  l_domain,
                                        l_criteria,
                                        CGAL::parameters::no_exude(),
                                        CGAL::parameters::no_perturb() );
  // Perturbation (maximum cpu time: 60s, targeted dihedral angle: default)
  CGAL::perturb_mesh_3( c3t3,
                        l_domain,
                        CGAL::parameters::time_limit = 60 );
  // Exudation
  CGAL::exude_mesh_3( c3t3,
                      CGAL::parameters::time_limit = 180 );


  // Output
  EDGE_LOG_INFO << "writing surface mesh";
  std::ofstream off_file("o_multImplDomains.off");
  std::ofstream bdry_file("o_bdry.off");
  c3t3.output_facets_in_complex_to_off( off_file );
  c3t3.output_boundary_to_off( bdry_file );

  EDGE_LOG_INFO << "thank you for using EDGEcut!";
}
