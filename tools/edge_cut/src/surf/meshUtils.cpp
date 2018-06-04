
#include "surf/meshUtils.h"

// Sizing field
// The vector of depth values, "m_layers," must be in decreasing order
//
//     baseValue
// ----------------- layers[0]
//     scales[0]
// ----------------- layers[1]
//     scales[1]
// ----------------- layers[2]
//       ...
//
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
edge_cut::surf::topoIntersect( Poly_slicer& slicer, K::Plane_3 plane ) {
  Polylines     polylines;
  Polyline_type ridges;

  slicer( plane, std::back_inserter( polylines ) );
  std::cout << "The plane " << plane.a() << "x + " << plane.b() << "y + "
            << plane.c() << "z = " << -1*plane.d() << " intersects the topography at "
            << polylines.size() << " polylines." << std::endl;
  if( polylines.size() == 1 ) {
    ridges = *(polylines.begin());
  } else {
    std::cout << "Error: Polyline slicer did not return precisely one intersection!" << std::endl;
  }
  return ridges;
}


// Returns True if polyline is increasing or decreasing in the (i_n+1)th
// component. If i_inc is True, checks if increasing. Otherwise, checks if
// decreasing.
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

// Surf_mesh
// edge_cut::surf::makePlaneSeg( K::Point_3 i_p1)
