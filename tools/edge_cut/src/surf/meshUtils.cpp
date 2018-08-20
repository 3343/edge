
#include "surf/meshUtils.h"

edge_cut::surf::K::FT
edge_cut::surf::DepthSizingField::operator()( const Point_3& p, const int, const Index& ) const
{
  // Check if point is above first layer
  if ( p.z() > m_layers[0] )
    return m_baseValue;

  // Test whether point is above increasingly deep layers
  for ( std::size_t l_layer = 1; l_layer < m_layers.size(); l_layer++ ) {
    if ( p.z() > m_layers[ l_layer ] )
      return m_baseValue * m_scales[ l_layer - 1 ];
  }

  // Point must be below deepest layer
  return m_baseValue * m_scales[ m_layers.size() - 1 ];
}


edge_cut::surf::Polyline_type
edge_cut::surf::topoIntersect( Poly_slicer& i_slicer, K::Plane_3 i_plane ) {
  Polylines     polylines;
  Polyline_type ridges;

  i_slicer( i_plane, std::back_inserter( polylines ) );
  std::cout << "The plane " << i_plane.a() << "x + " << i_plane.b() << "y + "
            << i_plane.c() << "z = " << -1*i_plane.d() << " intersects the topography at "
            << polylines.size() << " polylines." << std::endl;
  if( polylines.size() == 1 ) {
    ridges = *(polylines.begin());
  } else {
    std::cout << "Error: Polyline slicer did not return precisely one intersection!" << std::endl;
  }
  return ridges;
}


bool
edge_cut::surf::checkMonotonic( Polyline_type & i_p, unsigned int i_n, bool i_inc ) {
  if ( i_n > 2 ) {
    std::cout << "Error: Index out of bounds in checkMonotonic" << std::endl;
    return false;
  }
  for ( std::size_t l_idx = 1; l_idx < i_p.size(); l_idx++ ) {
    if ( i_inc ) {    // Check monotonically increasing
      if ( i_n == 0 )
      {
        if ( i_p[ l_idx - 1 ].x() >= i_p[ l_idx ].x() ) return false;
      }
      else if ( i_n == 1 )
      {
        if ( i_p[ l_idx - 1 ].y() >= i_p[ l_idx ].y() ) return false;
      }
      else if ( i_n == 2 )
      {
        if ( i_p[ l_idx - 1 ].z() >= i_p[ l_idx ].z() ) return false;
      }
    }
    else              // Check monotonically decreasing
    {
      if ( i_n == 0 )
      {
        if ( i_p[ l_idx - 1 ].x() <= i_p[ l_idx ].x() ) return false;
      }
      else if ( i_n == 1 )
      {
        if ( i_p[ l_idx - 1 ].y() <= i_p[ l_idx ].y() ) return false;
      }
      else if ( i_n == 2 )
      {
        if ( i_p[ l_idx - 1 ].z() <= i_p[ l_idx ].z() ) return false;
      }
    }
  }
  return true;
}


// Assumes that p is either in increasing or decreasing order in the (n+1)th
// component. Modifies p to ensure that it is increasing in (n+1)th component
void
edge_cut::surf::orderPolyline( Polyline_type & i_p, unsigned int i_n ) {
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
    std::cout << "Error: orderPolyline encountered a non-monotonic polyline." << std::endl;
    std::cout << "Printing polyline..." << std::endl;
    for ( auto const & l_pt : i_p )
      std::cout << l_pt << std::endl;

    std::cout << "Error: orderPolyline encountered a non-monotonic polyline." << std::endl;
  }

  return;
}


std::list< edge_cut::surf::Polyline_type >
edge_cut::surf::build1DFeatures(  Surf_mesh&  i_topoSurface,
                                  double      i_xMin,
                                  double      i_xMax,
                                  double      i_yMin,
                                  double      i_yMax,
                                  double      i_y0,
                                  double      i_zMin             )
{
  std::list< Polyline_type > features;
  Polyline_type faultRidges, yMinRidges, yMaxRidges, xMinRidges, xMaxRidges;
  Polyline_type xMinRidges0, xMinRidges1, xMaxRidges0, xMaxRidges1;

  Poly_slicer   topoSlicer( i_topoSurface );
  K::Plane_3    faultPlane = K::Plane_3( 0, 1, 0, -1 * i_y0 );
  K::Plane_3    yMinPlane  = K::Plane_3( 0, 1, 0, -1 * i_yMin );
  K::Plane_3    yMaxPlane  = K::Plane_3( 0, 1, 0, -1 * i_yMax );
  K::Plane_3    xMinPlane  = K::Plane_3( 1, 0, 0, -1 * i_xMin );
  K::Plane_3    xMaxPlane  = K::Plane_3( 1, 0, 0, -1 * i_xMax );


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

  // TODO create function to remove duplicate code
  //
  // Trim polylines to computational domain and add corner vertices.
  // This is to avoid the case where there already exists a vertex close to one
  // of the corners of the domain. It is not clear what an acceptable tolerance
  // is so tol=10 has been chosen arbitrarily.
  double const tol = 10.;
  Polyline_type l_tempLine;
  l_tempLine.clear();
  for ( std::size_t idx = 0; idx < yMinRidges.size(); idx++ ) {
    if ( yMinRidges[ idx ].x() > i_xMin + tol && yMinRidges[ idx ].x() < i_xMax - tol ) {
      l_tempLine.push_back( yMinRidges[ idx ] );
    }
  }
  l_tempLine.insert( l_tempLine.begin(), g_topo.interpolatePt( i_xMin, i_yMin ) );
  l_tempLine.push_back( g_topo.interpolatePt( i_xMax, i_yMin ) );
  for ( std::size_t idx = 1; idx < l_tempLine.size(); idx++ ) {
    assert( l_tempLine[ idx ].x() > l_tempLine[ idx-1 ].x() );
  }
  yMinRidges = l_tempLine;

  l_tempLine.clear();
  for ( std::size_t idx = 0; idx < faultRidges.size(); idx++ ) {
    if ( faultRidges[ idx ].x() > i_xMin + tol && faultRidges[ idx ].x() < i_xMax - tol )
      l_tempLine.push_back( faultRidges[ idx ] );
  }
  l_tempLine.insert( l_tempLine.begin(), g_topo.interpolatePt( i_xMin, Y0 ) );
  l_tempLine.push_back( g_topo.interpolatePt( i_xMax, Y0 ) );
  for ( std::size_t idx = 1; idx < l_tempLine.size(); idx++ ) {
    assert( l_tempLine[ idx ].x() > l_tempLine[ idx-1 ].x() );
  }
  faultRidges = l_tempLine;

  l_tempLine.clear();
  for ( std::size_t idx = 0; idx < yMaxRidges.size(); idx++ ) {
    if ( yMaxRidges[ idx ].x() > i_xMin + tol && yMaxRidges[ idx ].x() < i_xMax - tol )
      l_tempLine.push_back( yMaxRidges[ idx ] );
  }
  l_tempLine.insert( l_tempLine.begin(), g_topo.interpolatePt( i_xMin, i_yMax ) );
  l_tempLine.push_back( g_topo.interpolatePt( i_xMax, i_yMax ) );
  for ( std::size_t idx = 1; idx < l_tempLine.size(); idx++ ) {
    assert( l_tempLine[ idx ].x() > l_tempLine[ idx-1 ].x() );
  }
  yMaxRidges = l_tempLine;


  xMinRidges0.clear();
  xMinRidges1.clear();
  for ( std::size_t idx = 0; idx < xMinRidges.size(); idx++ ) {
    if ( xMinRidges[ idx ].y() > i_yMin + tol && xMinRidges[ idx ].y() < i_y0 - tol )
      xMinRidges0.push_back( xMinRidges[ idx ] );
    if ( xMinRidges[ idx ].y() > i_y0 + tol && xMinRidges[ idx ].y() < i_yMax - tol )
      xMinRidges1.push_back( xMinRidges[ idx ] );
  }
  xMinRidges0.insert( xMinRidges0.begin(), g_topo.interpolatePt( i_xMin, i_yMin ) );
  xMinRidges0.push_back( g_topo.interpolatePt( i_xMin, i_y0 ) );
  xMinRidges1.insert( xMinRidges1.begin(), g_topo.interpolatePt( i_xMin, i_y0 ) );
  xMinRidges1.push_back( g_topo.interpolatePt( i_xMin, i_yMax ) );

  xMaxRidges0.clear();
  xMaxRidges1.clear();
  for ( std::size_t idx = 0; idx < xMaxRidges.size(); idx++ ) {
    if ( xMaxRidges[ idx ].y() > i_yMin + tol && xMaxRidges[ idx ].y() < i_y0 - tol )
      xMaxRidges0.push_back( xMaxRidges[ idx ] );
    if ( xMaxRidges[ idx ].y() > i_y0 + tol && xMaxRidges[ idx ].y() < i_yMax - tol )
      xMaxRidges1.push_back( xMaxRidges[ idx ] );
  }
  xMaxRidges0.insert( xMaxRidges0.begin(), g_topo.interpolatePt( i_xMax, i_yMin ) );
  xMaxRidges0.push_back( g_topo.interpolatePt( i_xMax, i_y0 ) );
  xMaxRidges1.insert( xMaxRidges1.begin(), g_topo.interpolatePt( i_xMax, i_y0 ) );
  xMaxRidges1.push_back( g_topo.interpolatePt( i_xMax, i_yMax ) );

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
                                    i_zMin                   ) } );
  features.push_back( { K::Point_3( yMinRidges.back().x(),
                                    yMinRidges.back().y(),
                                    yMinRidges.back().z()   ),
                        K::Point_3( yMinRidges.back().x(),
                                    yMinRidges.back().y(),
                                    i_zMin                   ) } );
  features.push_back( { K::Point_3( yMaxRidges.front().x(),
                                    yMaxRidges.front().y(),
                                    yMaxRidges.front().z()  ),
                        K::Point_3( yMaxRidges.front().x(),
                                    yMaxRidges.front().y(),
                                    i_zMin                   ) } );
  features.push_back( { K::Point_3( yMaxRidges.back().x(),
                                    yMaxRidges.back().y(),
                                    yMaxRidges.back().z()   ),
                        K::Point_3( yMaxRidges.back().x(),
                                    yMaxRidges.back().y(),
                                    i_zMin                   ) } );
  features.push_back( { faultRidges.front(),
                        K::Point_3( faultRidges.front().x(),
                                    faultRidges.front().y(),
                                    i_zMin                   ) } );
  features.push_back( { faultRidges.back(),
                        K::Point_3( faultRidges.back().x(),
                                    faultRidges.back().y(),
                                    i_zMin                   ) } );

  // Bottom of computational domain
  features.push_back( { K::Point_3( yMinRidges.front().x(),
                                    yMinRidges.front().y(),
                                    i_zMin                   ),
                        K::Point_3( faultRidges.front().x(),
                                    faultRidges.front().y(),
                                    i_zMin                   ) } );
  features.push_back( { K::Point_3( faultRidges.front().x(),
                                    faultRidges.front().y(),
                                    i_zMin                   ),
                        K::Point_3( yMaxRidges.front().x(),
                                    yMaxRidges.front().y(),
                                    i_zMin                   ) } );
  features.push_back( { K::Point_3( yMaxRidges.front().x(),
                                    yMaxRidges.front().y(),
                                    i_zMin                   ),
                        K::Point_3( yMaxRidges.back().x(),
                                    yMaxRidges.back().y(),
                                    i_zMin                   ) } );
  features.push_back( { K::Point_3( yMaxRidges.back().x(),
                                    yMaxRidges.back().y(),
                                    i_zMin                   ),
                        K::Point_3( faultRidges.back().x(),
                                    faultRidges.back().y(),
                                    i_zMin                   ) } );
  features.push_back( { K::Point_3( faultRidges.back().x(),
                                    faultRidges.back().y(),
                                    i_zMin                   ),
                        K::Point_3( yMinRidges.back().x(),
                                    yMinRidges.back().y(),
                                    i_zMin                   ) } );
  features.push_back( { K::Point_3( yMinRidges.back().x(),
                                    yMinRidges.back().y(),
                                    i_zMin                   ),
                        K::Point_3( yMinRidges.front().x(),
                                    yMinRidges.front().y(),
                                    i_zMin                   ) } );
  features.push_back( { K::Point_3( faultRidges.front().x(),
                                    faultRidges.front().y(),
                                    i_zMin                   ),
                        K::Point_3( faultRidges.back().x(),
                                    faultRidges.back().y(),
                                    i_zMin                   ) } );

  return features;
}


edge_cut::surf::Polyhedron& edge_cut::surf::makeFreeSurfBdry( Polyhedron& io_freeSurfBdry )
{
  K::Point_3 p0, p1, p2, p3, p4, p5, p6, p7;
  std::vector< K::Point_3 > l_points;

  l_points.emplace_back( X_MIN, Y_MIN, Z_MIN );
  l_points.emplace_back( X_MAX, Y_MIN, Z_MIN );
  l_points.emplace_back( X_MIN, Y_MAX, Z_MIN );
  l_points.emplace_back( X_MAX, Y_MAX, Z_MIN );
  l_points.emplace_back( X_MIN, Y_MIN, Z_MAX );
  l_points.emplace_back( X_MAX, Y_MIN, Z_MAX );
  l_points.emplace_back( X_MIN, Y_MAX, Z_MAX );
  l_points.emplace_back( X_MAX, Y_MAX, Z_MAX );

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
                                                                io_freeSurfBdry );

  return io_freeSurfBdry;
}
