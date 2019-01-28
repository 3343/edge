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
 * Utilities for managing mesh data structures in EDGEcut
 **/

#include <CGAL/IO/OFF_reader.h>  // Must include here to avoid multiple defn errors (CGAL bug as of v4.13)
#include "surf/meshUtils.h"

edge_cut::surf::Polyline_type
edge_cut::surf::topoIntersect(  Poly_slicer const & i_slicer,
                                K::Plane_3  const & i_plane   )
{
  Polylines     l_polylines;
  Polyline_type l_ridges;

  i_slicer( i_plane, std::back_inserter( l_polylines ) );

  if( l_polylines.size() == 1 ) {
    l_ridges = *( l_polylines.begin() );
  }
  else {
    EDGE_LOG_ERROR << "Polyline slicer did not return precisely one intersection!";
    EDGE_LOG_FATAL << "  Computed " << l_polylines.size() << " intersections";
  }
  return l_ridges;
}


bool
edge_cut::surf::checkMonotonic( Polyline_type const & i_p,
                                unsigned int          i_n,
                                bool                  i_inc )
{
  if ( i_n > 2 ) {
    EDGE_LOG_FATAL << "checkMonotonic: Input index is greater than 2";
    return false;
  }

  for ( std::size_t l_idx = 1; l_idx < i_p.size(); l_idx++ ) {
    if ( i_inc ) {
      // Check monotonically increasing
      if ( i_p[ l_idx - 1 ][ i_n ] >= i_p[ l_idx ][ i_n ] ) return false;
    }
    else {
      // Check monotonically decreasing
      if ( i_p[ l_idx - 1 ][ i_n ] <= i_p[ l_idx ][ i_n ] ) return false;
    }
  }

  return true;
}


void
edge_cut::surf::orderPolyline(  Polyline_type & io_polyline,
                                unsigned int    i_n         )
{
  // i_p is already monotonically increasing, so return
  if ( checkMonotonic( io_polyline, i_n, true) ) {
    return;
  }
  // i_p is monotonically decreasing, so reverse order of all points
  else if ( checkMonotonic( io_polyline, i_n, false ) ) {
    Polyline_type l_temp = io_polyline;
    for ( std::size_t l_idx = 0; l_idx < io_polyline.size(); l_idx++ ) {
      io_polyline[ l_idx ] = l_temp[ io_polyline.size() - 1 - l_idx ];
    }
  }
  // i_p is not monotonic in the (i_n+1)th coordinate, record an error
  else {
    EDGE_LOG_ERROR << "orderPolyline: encountered a non-monotonic polyline - printing polyline:";
    for ( auto const & l_pt : io_polyline )
      EDGE_LOG_ERROR << "  " << l_pt;
    EDGE_LOG_FATAL << "Cannot continue with unordered polyline features";
  }

  return;
}


std::list< edge_cut::surf::Polyline_type >
edge_cut::surf::getIntersectionFeatures(  Polyhedron  const & i_topoSurface,
                                          double      const * i_bBox        )
{
  std::list< Polyline_type > l_features;

  Poly_slicer   l_topoSlicer( i_topoSurface );
  Polyline_type l_yMinRidges, l_yMaxRidges, l_xMinRidges, l_xMaxRidges;

  K::Plane_3    l_yMinPlane  = K::Plane_3( 0, 1, 0, -1 * i_bBox[2] );
  K::Plane_3    l_yMaxPlane  = K::Plane_3( 0, 1, 0, -1 * i_bBox[3] );
  K::Plane_3    l_xMinPlane  = K::Plane_3( 1, 0, 0, -1 * i_bBox[0] );
  K::Plane_3    l_xMaxPlane  = K::Plane_3( 1, 0, 0, -1 * i_bBox[1] );

  l_yMinRidges  = edge_cut::surf::topoIntersect( l_topoSlicer, l_yMinPlane );
  l_yMaxRidges  = edge_cut::surf::topoIntersect( l_topoSlicer, l_yMaxPlane );
  l_xMinRidges  = edge_cut::surf::topoIntersect( l_topoSlicer, l_xMinPlane );
  l_xMaxRidges  = edge_cut::surf::topoIntersect( l_topoSlicer, l_xMaxPlane );
  edge_cut::surf::orderPolyline( l_yMinRidges, 0 );
  edge_cut::surf::orderPolyline( l_yMaxRidges, 0 );
  edge_cut::surf::orderPolyline( l_xMinRidges, 1 );
  edge_cut::surf::orderPolyline( l_xMaxRidges, 1 );
  l_features.push_back( l_yMinRidges );
  l_features.push_back( l_yMaxRidges );
  l_features.push_back( l_xMinRidges );
  l_features.push_back( l_xMaxRidges );

  K::Point_3 l_p000( i_bBox[0], i_bBox[2], i_bBox[4] );
  K::Point_3 l_p010( i_bBox[0], i_bBox[3], i_bBox[4] );
  K::Point_3 l_p100( i_bBox[1], i_bBox[2], i_bBox[4] );
  K::Point_3 l_p110( i_bBox[1], i_bBox[3], i_bBox[4] );

  // Straight edges of bounding box
  l_features.push_back( { l_p000, l_p010 } );
  l_features.push_back( { l_p000, l_p100 } );
  l_features.push_back( { l_p010, l_p110 } );
  l_features.push_back( { l_p100, l_p110 } );
  l_features.push_back( { l_p000, l_xMinRidges.front() } );
  l_features.push_back( { l_p010, l_xMinRidges.back() } );
  l_features.push_back( { l_p100, l_xMaxRidges.front() } );
  l_features.push_back( { l_p110, l_xMaxRidges.back() } );

  return l_features;
}


edge_cut::surf::Polyhedron&
edge_cut::surf::makeBdry( Polyhedron        & o_bdry,
                          double      const * i_bBox  )
{
  if( !o_bdry.is_empty() ) {
    EDGE_LOG_WARNING << "Overwriting existing polyhedron in makeBdry";
  }
  o_bdry.clear();

  K::Point_3 p0, p1, p2, p3, p4, p5, p6, p7;
  std::vector< K::Point_3 > l_points;

  l_points.emplace_back( i_bBox[0], i_bBox[2], i_bBox[4] );
  l_points.emplace_back( i_bBox[1], i_bBox[2], i_bBox[4] );
  l_points.emplace_back( i_bBox[0], i_bBox[3], i_bBox[4] );
  l_points.emplace_back( i_bBox[1], i_bBox[3], i_bBox[4] );
  l_points.emplace_back( i_bBox[0], i_bBox[2], i_bBox[5] );
  l_points.emplace_back( i_bBox[1], i_bBox[2], i_bBox[5] );
  l_points.emplace_back( i_bBox[0], i_bBox[3], i_bBox[5] );
  l_points.emplace_back( i_bBox[1], i_bBox[3], i_bBox[5] );

  std::vector< std::vector< std::size_t > > l_polygons = {  { 0, 1, 2 },
                                                            { 1, 2, 3 },
                                                            { 0, 1, 4 },
                                                            { 1, 4, 5 },
                                                            { 1, 3, 5 },
                                                            { 3, 5, 7 },
                                                            { 0, 2, 4 },
                                                            { 2, 4, 6 },
                                                            { 2, 3, 6 },
                                                            { 3, 6, 7 } };

  CGAL::Polygon_mesh_processing::orient_polygon_soup( l_points, l_polygons );
  CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(  l_points,
                                                                l_polygons,
                                                                o_bdry      );
  o_bdry.normalize_border();

  return o_bdry;
}


void
edge_cut::surf::c3t3ToPolyhedron( C3t3        const & i_c3t3,
                                  Polyhedron        & o_polyhedron )
{
  if( !o_polyhedron.is_empty() ) {
    EDGE_LOG_WARNING << "Overwriting existing polyhedron in c3t3ToPolyhedron";
  }
  o_polyhedron.clear();

  CGAL::facets_in_complex_3_to_triangle_mesh< C3t3, Polyhedron >( i_c3t3, o_polyhedron );
  o_polyhedron.normalize_border();

  EDGE_CHECK( !o_polyhedron.is_empty() );
  EDGE_CHECK( o_polyhedron.is_pure_triangle() );
  EDGE_CHECK( o_polyhedron.normalized_border_is_valid() );

  return;
}


void
edge_cut::surf::c3t3ToSurfMesh( C3t3      const & i_c3t3,
                                Surf_mesh       & o_surfMesh )
{
  if( !o_surfMesh.is_empty() ) {
    EDGE_LOG_WARNING << "Overwriting existing surface mesh in c3t3ToSurfMesh";
  }
  o_surfMesh.clear();

  CGAL::facets_in_complex_3_to_triangle_mesh< C3t3, Surf_mesh >( i_c3t3, o_surfMesh );

  EDGE_CHECK( !o_surfMesh.is_empty() );

  return;
}


edge_cut::surf::Polyhedron&
edge_cut::surf::topoPolyMeshFromXYZ(  Polyhedron        & o_topoPolyMesh,
                                      std::string const & i_topoFile )
{
  if( !o_topoPolyMesh.is_empty() ) {
    EDGE_LOG_WARNING << "Overwriting existing polyhedron in topoPolyMeshFromXYZ";
  }
  o_topoPolyMesh.clear();

  std::stringstream     l_topoStream;
  edge_cut::surf::Topo  l_topo( i_topoFile );

  l_topo.writeTriaToOff( l_topoStream );

  if (  !l_topoStream ||
        !(l_topoStream >> o_topoPolyMesh) ||
        o_topoPolyMesh.is_empty() ||
        !CGAL::is_triangle_mesh( o_topoPolyMesh ) ) {
    EDGE_LOG_FATAL << "Input stream for topography triangulation is not valid.";
  }
  l_topoStream.str( std::string() );

  return o_topoPolyMesh;
}
