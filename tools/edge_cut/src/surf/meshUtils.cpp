#include <CGAL/IO/OFF_reader.h>  // We get linker error if this included in header. TODO
#include "surf/meshUtils.h"
#include "io/logging.hpp"

// TODO check out gmsh approach to linear scaled fields
edge_cut::surf::K::FT
edge_cut::surf::SizingField::operator()( const Point_3& p, const int, const Index& ) const
{
  FT l_distance = std::sqrt( std::pow( p.x()-m_center.x(), 2 ) + std::pow( p.y()-m_center.y(), 2 ) ); // TODO change to 2D distance
  // Check if point is above first layer
  if ( l_distance < m_innerRad )
    return m_innerVal;
  else if ( l_distance >= m_outerRad )
    return m_scale * m_innerVal;
  else {
    // Control should never reach here if inner and outer radii are equal
    assert( m_outerRad != m_innerRad );
    return ( 1 + ( l_distance - m_innerRad ) * ( m_scale - 1) / ( m_outerRad - m_innerRad ) ) * m_innerVal;
  }
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


// std::list< edge_cut::surf::Polyline_type >
// edge_cut::surf::build1DFeatures(  Surf_mesh&  i_topoSurface,
//                                   double      i_xMin,
//                                   double      i_xMax,
//                                   double      i_yMin,
//                                   double      i_yMax,
//                                   double      i_y0,
//                                   double      i_zMin             )
// {
//   std::list< Polyline_type > features;
//   Polyline_type faultRidges, yMinRidges, yMaxRidges, xMinRidges, xMaxRidges;
//   Polyline_type xMinRidges0, xMinRidges1, xMaxRidges0, xMaxRidges1;
//
//   Poly_slicer   topoSlicer( i_topoSurface );
//   K::Plane_3    faultPlane = K::Plane_3( 0, 1, 0, -1 * i_y0 );
//   K::Plane_3    yMinPlane  = K::Plane_3( 0, 1, 0, -1 * i_yMin );
//   K::Plane_3    yMaxPlane  = K::Plane_3( 0, 1, 0, -1 * i_yMax );
//   K::Plane_3    xMinPlane  = K::Plane_3( 1, 0, 0, -1 * i_xMin );
//   K::Plane_3    xMaxPlane  = K::Plane_3( 1, 0, 0, -1 * i_xMax );
//
//
//   // Intersect fault plane with topography
//   faultRidges = topoIntersect( topoSlicer, faultPlane );
//   yMinRidges  = topoIntersect( topoSlicer, yMinPlane );
//   yMaxRidges  = topoIntersect( topoSlicer, yMaxPlane );
//   xMinRidges  = topoIntersect( topoSlicer, xMinPlane );
//   xMaxRidges  = topoIntersect( topoSlicer, xMaxPlane );
//   orderPolyline( faultRidges, 0 );
//   orderPolyline( yMinRidges, 0 );
//   orderPolyline( yMaxRidges, 0 );
//   orderPolyline( xMinRidges, 1 );
//   orderPolyline( xMaxRidges, 1 );
//
//   // TODO create function to remove duplicate code
//   //
//   // Trim polylines to computational domain and add corner vertices.
//   // This is to avoid the case where there already exists a vertex close to one
//   // of the corners of the domain. It is not clear what an acceptable tolerance
//   // is so tol=10 has been chosen arbitrarily.
//   double const tol = 10.;
//   Polyline_type l_tempLine;
//   l_tempLine.clear();
//   for ( std::size_t idx = 0; idx < yMinRidges.size(); idx++ ) {
//     if ( yMinRidges[ idx ].x() > i_xMin + tol && yMinRidges[ idx ].x() < i_xMax - tol ) {
//       l_tempLine.push_back( yMinRidges[ idx ] );
//     }
//   }
//   l_tempLine.insert( l_tempLine.begin(), g_topo.interpolatePt( i_xMin, i_yMin ) );
//   l_tempLine.push_back( g_topo.interpolatePt( i_xMax, i_yMin ) );
//   for ( std::size_t idx = 1; idx < l_tempLine.size(); idx++ ) {
//     assert( l_tempLine[ idx ].x() > l_tempLine[ idx-1 ].x() );
//   }
//   yMinRidges = l_tempLine;
//
//   l_tempLine.clear();
//   for ( std::size_t idx = 0; idx < faultRidges.size(); idx++ ) {
//     if ( faultRidges[ idx ].x() > i_xMin + tol && faultRidges[ idx ].x() < i_xMax - tol )
//       l_tempLine.push_back( faultRidges[ idx ] );
//   }
//   l_tempLine.insert( l_tempLine.begin(), g_topo.interpolatePt( i_xMin, Y0 ) );
//   l_tempLine.push_back( g_topo.interpolatePt( i_xMax, Y0 ) );
//   for ( std::size_t idx = 1; idx < l_tempLine.size(); idx++ ) {
//     assert( l_tempLine[ idx ].x() > l_tempLine[ idx-1 ].x() );
//   }
//   faultRidges = l_tempLine;
//
//   l_tempLine.clear();
//   for ( std::size_t idx = 0; idx < yMaxRidges.size(); idx++ ) {
//     if ( yMaxRidges[ idx ].x() > i_xMin + tol && yMaxRidges[ idx ].x() < i_xMax - tol )
//       l_tempLine.push_back( yMaxRidges[ idx ] );
//   }
//   l_tempLine.insert( l_tempLine.begin(), g_topo.interpolatePt( i_xMin, i_yMax ) );
//   l_tempLine.push_back( g_topo.interpolatePt( i_xMax, i_yMax ) );
//   for ( std::size_t idx = 1; idx < l_tempLine.size(); idx++ ) {
//     assert( l_tempLine[ idx ].x() > l_tempLine[ idx-1 ].x() );
//   }
//   yMaxRidges = l_tempLine;
//
//
//   xMinRidges0.clear();
//   xMinRidges1.clear();
//   for ( std::size_t idx = 0; idx < xMinRidges.size(); idx++ ) {
//     if ( xMinRidges[ idx ].y() > i_yMin + tol && xMinRidges[ idx ].y() < i_y0 - tol )
//       xMinRidges0.push_back( xMinRidges[ idx ] );
//     if ( xMinRidges[ idx ].y() > i_y0 + tol && xMinRidges[ idx ].y() < i_yMax - tol )
//       xMinRidges1.push_back( xMinRidges[ idx ] );
//   }
//   xMinRidges0.insert( xMinRidges0.begin(), g_topo.interpolatePt( i_xMin, i_yMin ) );
//   xMinRidges0.push_back( g_topo.interpolatePt( i_xMin, i_y0 ) );
//   xMinRidges1.insert( xMinRidges1.begin(), g_topo.interpolatePt( i_xMin, i_y0 ) );
//   xMinRidges1.push_back( g_topo.interpolatePt( i_xMin, i_yMax ) );
//
//   xMaxRidges0.clear();
//   xMaxRidges1.clear();
//   for ( std::size_t idx = 0; idx < xMaxRidges.size(); idx++ ) {
//     if ( xMaxRidges[ idx ].y() > i_yMin + tol && xMaxRidges[ idx ].y() < i_y0 - tol )
//       xMaxRidges0.push_back( xMaxRidges[ idx ] );
//     if ( xMaxRidges[ idx ].y() > i_y0 + tol && xMaxRidges[ idx ].y() < i_yMax - tol )
//       xMaxRidges1.push_back( xMaxRidges[ idx ] );
//   }
//   xMaxRidges0.insert( xMaxRidges0.begin(), g_topo.interpolatePt( i_xMax, i_yMin ) );
//   xMaxRidges0.push_back( g_topo.interpolatePt( i_xMax, i_y0 ) );
//   xMaxRidges1.insert( xMaxRidges1.begin(), g_topo.interpolatePt( i_xMax, i_y0 ) );
//   xMaxRidges1.push_back( g_topo.interpolatePt( i_xMax, i_yMax ) );
//
//   features.push_back( faultRidges );
//   features.push_back( yMinRidges );
//   features.push_back( yMaxRidges );
//   features.push_back( xMinRidges0 );
//   features.push_back( xMinRidges1 );
//   features.push_back( xMaxRidges0 );
//   features.push_back( xMaxRidges1 );
//
//   // Vertical Edges
//   features.push_back( { K::Point_3( yMinRidges.front().x(),
//                                     yMinRidges.front().y(),
//                                     yMinRidges.front().z()  ),
//                         K::Point_3( yMinRidges.front().x(),
//                                     yMinRidges.front().y(),
//                                     i_zMin                   ) } );
//   features.push_back( { K::Point_3( yMinRidges.back().x(),
//                                     yMinRidges.back().y(),
//                                     yMinRidges.back().z()   ),
//                         K::Point_3( yMinRidges.back().x(),
//                                     yMinRidges.back().y(),
//                                     i_zMin                   ) } );
//   features.push_back( { K::Point_3( yMaxRidges.front().x(),
//                                     yMaxRidges.front().y(),
//                                     yMaxRidges.front().z()  ),
//                         K::Point_3( yMaxRidges.front().x(),
//                                     yMaxRidges.front().y(),
//                                     i_zMin                   ) } );
//   features.push_back( { K::Point_3( yMaxRidges.back().x(),
//                                     yMaxRidges.back().y(),
//                                     yMaxRidges.back().z()   ),
//                         K::Point_3( yMaxRidges.back().x(),
//                                     yMaxRidges.back().y(),
//                                     i_zMin                   ) } );
//   features.push_back( { faultRidges.front(),
//                         K::Point_3( faultRidges.front().x(),
//                                     faultRidges.front().y(),
//                                     i_zMin                   ) } );
//   features.push_back( { faultRidges.back(),
//                         K::Point_3( faultRidges.back().x(),
//                                     faultRidges.back().y(),
//                                     i_zMin                   ) } );
//
//   // Bottom of computational domain
//   features.push_back( { K::Point_3( yMinRidges.front().x(),
//                                     yMinRidges.front().y(),
//                                     i_zMin                   ),
//                         K::Point_3( faultRidges.front().x(),
//                                     faultRidges.front().y(),
//                                     i_zMin                   ) } );
//   features.push_back( { K::Point_3( faultRidges.front().x(),
//                                     faultRidges.front().y(),
//                                     i_zMin                   ),
//                         K::Point_3( yMaxRidges.front().x(),
//                                     yMaxRidges.front().y(),
//                                     i_zMin                   ) } );
//   features.push_back( { K::Point_3( yMaxRidges.front().x(),
//                                     yMaxRidges.front().y(),
//                                     i_zMin                   ),
//                         K::Point_3( yMaxRidges.back().x(),
//                                     yMaxRidges.back().y(),
//                                     i_zMin                   ) } );
//   features.push_back( { K::Point_3( yMaxRidges.back().x(),
//                                     yMaxRidges.back().y(),
//                                     i_zMin                   ),
//                         K::Point_3( faultRidges.back().x(),
//                                     faultRidges.back().y(),
//                                     i_zMin                   ) } );
//   features.push_back( { K::Point_3( faultRidges.back().x(),
//                                     faultRidges.back().y(),
//                                     i_zMin                   ),
//                         K::Point_3( yMinRidges.back().x(),
//                                     yMinRidges.back().y(),
//                                     i_zMin                   ) } );
//   features.push_back( { K::Point_3( yMinRidges.back().x(),
//                                     yMinRidges.back().y(),
//                                     i_zMin                   ),
//                         K::Point_3( yMinRidges.front().x(),
//                                     yMinRidges.front().y(),
//                                     i_zMin                   ) } );
//   features.push_back( { K::Point_3( faultRidges.front().x(),
//                                     faultRidges.front().y(),
//                                     i_zMin                   ),
//                         K::Point_3( faultRidges.back().x(),
//                                     faultRidges.back().y(),
//                                     i_zMin                   ) } );
//
//   return features;
// }


edge_cut::surf::Polyhedron& edge_cut::surf::makeFreeSurfBdry( Polyhedron        & io_bdry,
                                                              Topo        const & i_topo  )
{
  K::Point_3 p0, p1, p2, p3, p4, p5, p6, p7;
  std::vector< K::Point_3 > l_points;

  l_points.emplace_back( i_topo.m_xMin, i_topo.m_yMin, i_topo.m_zMin );
  l_points.emplace_back( i_topo.m_xMax, i_topo.m_yMin, i_topo.m_zMin );
  l_points.emplace_back( i_topo.m_xMin, i_topo.m_yMax, i_topo.m_zMin );
  l_points.emplace_back( i_topo.m_xMax, i_topo.m_yMax, i_topo.m_zMin );
  l_points.emplace_back( i_topo.m_xMin, i_topo.m_yMin, i_topo.m_zMax );
  l_points.emplace_back( i_topo.m_xMax, i_topo.m_yMin, i_topo.m_zMax );
  l_points.emplace_back( i_topo.m_xMin, i_topo.m_yMax, i_topo.m_zMax );
  l_points.emplace_back( i_topo.m_xMax, i_topo.m_yMax, i_topo.m_zMax );

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

// TODO WARNING will not work if there are not at least two facets above
// the topography on the free surface boundary
void edge_cut::surf::c3t3ToPolyhedron(  C3t3        const & c3t3,
                                        Polyhedron        & polyhedron )
{
  std::stringstream sstream;
  std::vector<K::Point_3> points;
  std::vector< std::vector<std::size_t> > polygons;

  CGAL::internal::output_facets_in_complex_to_off<C3t3>( c3t3, sstream ); // TODO change to c3t3 member function

  if (!CGAL::read_OFF( sstream, points, polygons))
  {
    std::cerr << "c3t3ToPolyhedron: " << std::endl;
    std::cerr << "Error parsing the OFF file " << std::endl;
    return;
  }
  CGAL::Polygon_mesh_processing::orient_polygon_soup( points, polygons );
  CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh( points, polygons, polyhedron );
  polyhedron.normalize_border();

  return;
}


edge_cut::surf::Polyhedron& edge_cut::surf::trimFSB(  Polyhedron const & i_topo,
                                                      Polyhedron       & io_fsb )
{
  // Const types are used for the (unchanged) topo mesh
  typedef Polyhedron::Vertex_const_handle       Vertex_const;
  typedef Polyhedron::Halfedge_const_handle     Halfedge_const;

  // Non-const types are for the boundary mesh
  typedef Polyhedron::Vertex_handle             Vertex;
  typedef Polyhedron::Vertex_iterator           VertexIt;
  typedef Polyhedron::Halfedge_handle           Halfedge;
  typedef Polyhedron::Halfedge_const_iterator   HalfedgeIt;

  // We cannot compare vertices and halfedges on two different polyhedral meshes,
  // because they will always have different connectivity information.
  // The Point type has no connectivity data, so it can be used to compare
  // vertices on two different meshes.
  typedef Polyhedron::Point                     Point;


  Vertex l_vertFsb, l_vertPrevFsb;
  Vertex_const l_vertTopo;

  // Get a handle to a vertex/halfedge pair on the border of the topo mesh
  // This is where we will start our cut of the surface
  io_fsb.normalize_border();
  if ( !i_topo.normalized_border_is_valid() ) {
    EDGE_LOG_INFO << "Attempted to trim border mesh that was not normalized - mesh will not be trimmed";
    return io_fsb;
  }
  HalfedgeIt l_borderIt = i_topo.border_halfedges_begin();
  l_vertTopo = l_borderIt->vertex();

  // Get a handle within the fsb mesh to the same vertex
  VertexIt l_vIt = io_fsb.vertices_begin();
  while( l_vIt->point() != l_vertTopo->point() ) {
    l_vIt++;
    if ( l_vIt == io_fsb.vertices_end() ) {
      EDGE_LOG_INFO << "Could not find coincident vertex between topography and border mesh - mesh will not be trimmed";
      return io_fsb;
    }
  }
  l_vertFsb = l_vIt;

  // Get a handle within fsb mesh to the "previous" border vertex when circulating the edge of topography
  l_vIt = io_fsb.vertices_begin();
  Vertex_const l_vertPrevTopo;
  auto l_edge = l_vertTopo->vertex_begin();
  while( !l_edge->is_border_edge() ){ l_edge++; }
  l_vertPrevTopo = l_edge->prev()->vertex();

  while( l_vIt->point() != l_vertPrevTopo->point() ) { l_vIt++; }
  l_vertPrevFsb = l_vIt;


  // Get the first halfedge on FSB, incident to vertTopo, and ABOVE the topography
  bool orderedPos;
  Halfedge l_h, l_hNext;
  auto l_edgeCirc = l_vertFsb->vertex_begin();
  while( l_edgeCirc->prev()->vertex()->point() != l_vertPrevFsb->point() ) l_edgeCirc++;
  Halfedge l_rootHalfedge = l_edgeCirc;

  // TODO WARNING Need better test here: this could potentially fail to find
  // the "above" component if topography is steep
  if( l_rootHalfedge->next()->vertex()->point().z() > l_rootHalfedge->opposite()->next()->vertex()->point().z() ){
    orderedPos = true;
    l_h = l_rootHalfedge;
  } else {
    orderedPos = false;
    l_h = l_rootHalfedge->opposite();
  }

  // For each topography vertex, delete all adjacent facets on the fsb mesh
  // which are adjacent to the vertex and above the topography
  Halfedge_const l_borderHalfedge = l_borderIt;
  Halfedge_const l_nextBorderHalfedge;
  Point l_startPoint = l_borderHalfedge->vertex()->point();  //TODO using points here in case halfedge handles get invalidated/changed
  do {
    // Compute the next halfedge on the border of the topography
    auto l_borderEdgeCirc = l_borderHalfedge->vertex()->vertex_begin();
    while( ! l_borderEdgeCirc->next()->is_border_edge() ){ l_borderEdgeCirc++; }
    l_nextBorderHalfedge = l_borderEdgeCirc->next();

    // Erase all but the last facet adjacent to the vertex which are above topography
    // The last facet is also adjacent to the next border vertex, and will be
    // deleted at the start of the next loop
    while( l_h->next()->vertex()->point() != l_nextBorderHalfedge->vertex()->point() ){

      // If halfedge is a border halfedge, then there are no more facets to delete TODO
      if ( l_h->is_border() ) break;

      if( orderedPos )  l_hNext = l_h->next()->opposite();
      else              l_hNext = l_h->prev()->opposite();

      io_fsb.erase_facet( l_h );
      l_h = l_hNext;
    }

    if( orderedPos ) l_h = l_h->next();
    else             l_h = l_h->prev();

    l_borderHalfedge = l_nextBorderHalfedge;
  } while ( l_borderHalfedge->vertex()->point() != l_startPoint );

  // Iterate over all non-border halfedges for the one with the highest elevation
  // This halfedge will be on the "upper" connected component
  // TODO this assumes that the trim splits the mesh into two components
  io_fsb.normalize_border();
  Halfedge l_hMax = io_fsb.halfedges_begin();
  for ( auto l_it = io_fsb.halfedges_begin(); l_it != io_fsb.border_halfedges_begin(); l_it++ ) {
    if ( l_it->vertex()->point().z() > l_hMax->vertex()->point().z() )
      l_hMax = l_it;
  }

  io_fsb.erase_connected_component( l_hMax );
  io_fsb.normalize_border();

  return io_fsb;
}
