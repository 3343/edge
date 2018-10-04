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
 * Utilities for managing mesh data structures in edge_cut
 **/

#include <CGAL/IO/OFF_reader.h>  // Must include here to avoid multiple defn errors (CGAL bug as of v4.13)
#include "surf/meshUtils.h"
#include "io/logging.hpp"

edge_cut::surf::K::FT
edge_cut::surf::SizingField::operator()(  Point_3 const & p,
                                          int     const,
                                          Index   const &   ) const
{
  FT l_distance = std::sqrt( std::pow( p.x()-m_center.x(), 2 ) + std::pow( p.y()-m_center.y(), 2 ) );

  if ( l_distance <= m_innerRad )
    return m_innerVal;
  else if ( l_distance >= m_outerRad )
    return m_scale * m_innerVal;
  else {
    // Control should never reach here if inner and outer radii are equal
    assert( m_outerRad != m_innerRad );
    return ( 1 + ( l_distance - m_innerRad ) * ( m_scale - 1) * m_widthInv ) * m_innerVal;
  }
}


edge_cut::surf::Polyline_type
edge_cut::surf::topoIntersect(  Poly_slicer & i_slicer,
                                K::Plane_3    i_plane   )
{
  Polylines     polylines;
  Polyline_type ridges;

  i_slicer( i_plane, std::back_inserter( polylines ) );

  if( polylines.size() == 1 ) {
    ridges = *(polylines.begin());
  }
  else {
    EDGE_LOG_ERROR << "Error: Polyline slicer did not return precisely one intersection!";
  }
  return ridges;
}


bool
edge_cut::surf::checkMonotonic( Polyline_type & i_p,
                                unsigned int    i_n,
                                bool            i_inc )
{
  if ( i_n > 2 ) {
    EDGE_LOG_ERROR << "checkMonotonic:";
    EDGE_LOG_ERROR << "  Input index is greater than 2";
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
edge_cut::surf::orderPolyline(  Polyline_type & i_p,
                                unsigned int    i_n )
{
  // i_p is already monotonically increasing, so return
  if ( checkMonotonic( i_p, i_n, true) )
    return;
  // i_p is monotonically decreasing, so reverse order of all points
  else if ( checkMonotonic( i_p, i_n, false ) ) {
    Polyline_type l_temp = i_p;
    for ( std::size_t l_idx = 0; l_idx < i_p.size(); l_idx++ ) {
      i_p[ l_idx ] = l_temp[ i_p.size() - 1 - l_idx ];
    }
  }
  // i_p is not monotonic in the (i_n+1)th coordinate, record an error
  else {
    EDGE_LOG_ERROR << "orderPolyline:";
    EDGE_LOG_ERROR << "  encountered a non-monotonic polyline - printing polyline:";
    for ( auto const & l_pt : i_p )
      EDGE_LOG_ERROR << "  " << l_pt;
  }

  return;
}


std::list< edge_cut::surf::Polyline_type >
edge_cut::surf::getIntersectionFeatures(  Polyhedron  const & i_topoSurface,
                                          double      const * i_bBox        )
{
  std::list< Polyline_type > features;

  EDGE_LOG_INFO << "Computing topography-boundary intersection...";
  Poly_slicer   topoSlicer( i_topoSurface );
  Polyline_type faultRidges, yMinRidges, yMaxRidges, xMinRidges, xMaxRidges;
  K::Plane_3    yMinPlane  = K::Plane_3( 0, 1, 0, -1 * i_bBox[2] );
  K::Plane_3    yMaxPlane  = K::Plane_3( 0, 1, 0, -1 * i_bBox[3] );
  K::Plane_3    xMinPlane  = K::Plane_3( 1, 0, 0, -1 * i_bBox[0] );
  K::Plane_3    xMaxPlane  = K::Plane_3( 1, 0, 0, -1 * i_bBox[1] );

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

  K::Point_3 p000( i_bBox[0], i_bBox[2], i_bBox[4] );
  K::Point_3 p010( i_bBox[0], i_bBox[3], i_bBox[4] );
  K::Point_3 p100( i_bBox[1], i_bBox[2], i_bBox[4] );
  K::Point_3 p110( i_bBox[1], i_bBox[3], i_bBox[4] );

  // Straight edges of bounding box
  features.push_back( { p000, p010 } );
  features.push_back( { p000, p100 } );
  features.push_back( { p010, p110 } );
  features.push_back( { p100, p110 } );
  features.push_back( { p000, xMinRidges.front() } );
  features.push_back( { p010, xMinRidges.back() } );
  features.push_back( { p100, xMaxRidges.front() } );
  features.push_back( { p110, xMaxRidges.back() } );

  return features;
}


edge_cut::surf::Polyhedron&
edge_cut::surf::makeBdry( Polyhedron        & io_bdry,
                          double      const * i_bBox  )
{
  EDGE_LOG_INFO << "Constructing polyhedral surface model of domain boundary...";
  if( !io_bdry.is_empty() )
    EDGE_LOG_WARNING << "Overwriting existing polyhedron in makeBdry";
  io_bdry.clear();

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
                                                                io_bdry      );
  io_bdry.normalize_border();

  return io_bdry;
}


void
edge_cut::surf::c3t3ToPolyhedron( C3t3        const & c3t3,
                                  Polyhedron        & polyhedron )
{
  if( !polyhedron.is_empty() )
    EDGE_LOG_WARNING << "Overwriting existing polyhedron in c3t3ToPolyhedron";
  polyhedron.clear();

  std::stringstream sstream;
  std::vector<K::Point_3> points;
  std::vector< std::vector<std::size_t> > polygons;

  c3t3.output_facets_in_complex_to_off( sstream );

  if (!CGAL::read_OFF( sstream, points, polygons)) {
    EDGE_LOG_ERROR << "c3t3ToPolyhedron: " << std::endl;
    EDGE_LOG_ERROR << "Error parsing the OFF stream " << std::endl;
    return;
  }

  CGAL::Polygon_mesh_processing::orient_polygon_soup( points, polygons );

  if( !CGAL::Polygon_mesh_processing::is_polygon_soup_a_polygon_mesh( polygons ) ) {
    EDGE_LOG_ERROR << "c3t3ToPolyhedron: " << std::endl;
    EDGE_LOG_ERROR << "Polygon soup is not a polygon mesh" << std::endl;
    return;
  }

  CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh( points, polygons, polyhedron );
  polyhedron.normalize_border();

  return;
}


void
edge_cut::surf::c3t3ToSurfMesh( C3t3      const & i_c3t3,
                                Surf_mesh       & io_surfMesh )
{
  if( !io_surfMesh.is_empty() )
    EDGE_LOG_WARNING << "Overwriting existing surface mesh in c3t3ToSurfMesh";
  io_surfMesh.clear();

  std::stringstream sstream;
  std::vector<K::Point_3> points;
  std::vector< std::vector<std::size_t> > polygons;

  i_c3t3.output_facets_in_complex_to_off( sstream );

  if (!CGAL::read_OFF( sstream, points, polygons)) {
    EDGE_LOG_ERROR << "c3t3ToSurfMesh: " << std::endl;
    EDGE_LOG_ERROR << "Error parsing the OFF stream " << std::endl;
    return;
  }

  CGAL::Polygon_mesh_processing::orient_polygon_soup( points, polygons );

  if( !CGAL::Polygon_mesh_processing::is_polygon_soup_a_polygon_mesh( polygons ) ) {
    EDGE_LOG_ERROR << "c3t3ToSurfMesh: " << std::endl;
    EDGE_LOG_ERROR << "Polygon soup is not a polygon mesh" << std::endl;
    return;
  }

  CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh( points, polygons, io_surfMesh );

  return;
}


edge_cut::surf::Polyhedron&
edge_cut::surf::topoPolyMeshFromXYZ(  Polyhedron& io_topoPolyMesh,
                                      std::string const & i_topoFile )
{
  EDGE_LOG_INFO << "Constructing Delaunay triangulation for topography...";
  edge_cut::surf::Topo l_topo( i_topoFile );

  EDGE_LOG_INFO << "Constructing polyhedral surface model of topography...";
  std::stringstream   topoStream;
  l_topo.writeTriaToOff( topoStream );

  if (!topoStream || !(topoStream >> io_topoPolyMesh) || io_topoPolyMesh.is_empty()
             || !CGAL::is_triangle_mesh( io_topoPolyMesh ) ) {
    std::cerr << "Input stream for topography triangulation is not valid." << std::endl;
  }
  topoStream.str( std::string() );

  return io_topoPolyMesh;
}
