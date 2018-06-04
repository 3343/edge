#define CGAL_MESH_3_PROTECTION_DEBUG 1
#define CGAL_MESH_3_VERBOSE 1
#define CGAL_MESH_3_PROFILING 1
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Implicit_mesh_domain_3.h>
#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/Mesh_domain_with_polyline_features_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Triangulation_data_structure.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <boost/optional/optional_io.hpp>


// IO
#include <CGAL/IO/Polyhedron_iostream.h>
// Ouput
#include <CGAL/Mesh_3/Dump_c3t3.h>
// Read 1D features from input file
#include "read_polylines.h"
#include "surf/Topo.h"
#include "implicit_functions.h"
#include "io/logging.hpp"
INITIALIZE_EASYLOGGINGPP


// Sphere Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::FT FT;
typedef K::Point_3 Point;
typedef FT_to_point_function_wrapper<K::FT, K::Point_3> Function;
typedef CGAL::Implicit_multi_domain_to_labeling_function_wrapper<Function>
                                                        Function_wrapper;
typedef Function_wrapper::Function_vector Function_vector;
typedef CGAL::Labeled_mesh_domain_3<Function_wrapper, K> Implicit_domain;
// typedef FT (Function)(const Point&);
// typedef CGAL::Implicit_mesh_domain_3<Function,K> Implicit_domain;


// Polyhedral Domain
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef CGAL::Polyhedral_mesh_domain_3<Polyhedron, K> Polyhedron_domain;


class Hybrid_domain
{
  const Implicit_domain& implicit_domain;
  const Polyhedron_domain& polyhedron_domain;
public:
  Hybrid_domain(const Implicit_domain& implicit_domain,
                const Polyhedron_domain& polyhedron_domain)
    : implicit_domain(implicit_domain)
    , polyhedron_domain(polyhedron_domain)
  {}
  // types required by the `MeshDomain_3` concept
  typedef int Surface_patch_index;
  typedef int Subdomain_index;
  typedef int Index;
  typedef K R;
  typedef K::Point_3 Point_3;
  typedef CGAL::cpp11::tuple<Point_3, Index, int> Intersection;
  CGAL::Bbox_3 bbox() const {
    return implicit_domain.bbox() + polyhedron_domain.bbox();
  }
  struct Construct_initial_points
  {
    Construct_initial_points(const Hybrid_domain& domain)
      : r_domain_(domain) {}
    template<class OutputIterator>
    OutputIterator operator()(OutputIterator pts, const int n = 20) const
    {
      //construct initial points on implicit domain
      typedef Implicit_domain::Index Implicit_Index;
      std::vector<std::pair<Point_3,
                            Implicit_Index> > implicit_points_vector;
      Implicit_domain::Construct_initial_points cstr_implicit_initial_points =
        r_domain_.implicit_domain.construct_initial_points_object();
      cstr_implicit_initial_points(std::back_inserter(implicit_points_vector),
                                   n/2);
      for(std::size_t i = 0, end = implicit_points_vector.size(); i<end; ++i) {
        *pts++ = std::make_pair(implicit_points_vector[i].first, 2);
      }
      //construct initial points on polyhedral domain
      typedef Polyhedron_domain::Index Polyhedron_Index;
      std::vector<std::pair<Point_3,
                            Polyhedron_Index> > polyhedron_points_vector;
      Polyhedron_domain::Construct_initial_points cstr_polyhedron_initial_points =
        r_domain_.polyhedron_domain.construct_initial_points_object();
      cstr_polyhedron_initial_points(std::back_inserter(polyhedron_points_vector),
                                     n/2);
      for(std::size_t i = 0, end = polyhedron_points_vector.size(); i<end; ++i) {
        *pts++ = std::make_pair(polyhedron_points_vector[i].first, 1);
      }
      return pts;
    }
  private:
    const Hybrid_domain& r_domain_;
  }; // end Construct_initial_points_object
  Construct_initial_points construct_initial_points_object() const
  {
    return Construct_initial_points(*this);
  }
  struct Is_in_domain
  {
    Is_in_domain(const Hybrid_domain& domain) : r_domain_(domain) {}
    boost::optional<Subdomain_index> operator()(const K::Point_3& p) const
    {
      boost::optional<Subdomain_index> subdomain_index =
        r_domain_.implicit_domain.is_in_domain_object()(p);
      if(subdomain_index){
        std::cout << "IMPLICIT " << subdomain_index << std::endl;
        return 2;
      }
      else {
        std::cout << p << std::endl;
        subdomain_index = r_domain_.polyhedron_domain.is_in_domain_object()(p);
        if(subdomain_index)
          std::cout << "POLY " << subdomain_index << std::endl;
        return r_domain_.polyhedron_domain.is_in_domain_object()(p);
      }
      // WARNING this assumes only two impl subdomains! FIXME
      // index = 0  <-- not contained
      // index = 1  <-- contained in impl subdomain 1
      // index = 2  <-- contained in impl subdomain 2
      // index = 3  <-- contained in polydedral domain
      // if(subdomain_index){
      //   std::cout << "impl subdomain index  " << subdomain_index << std::endl;
      //   return subdomain_index;
      // }
      // else{
      //   subdomain_index = r_domain_.polyhedron_domain.is_in_domain_object()(p);
      //   std::cout << "poly subdomain index  " << subdomain_index << std::endl;
      //   if( subdomain_index ) return 3;
      //   else return 0;
      // }
    }
  private:
    const Hybrid_domain& r_domain_;
  };
  Is_in_domain is_in_domain_object() const { return Is_in_domain(*this); }
  struct Construct_intersection
  {
    Construct_intersection(const Hybrid_domain& domain)
      : r_domain_(domain) {}
    template <typename Query>
    Intersection operator()(const Query& query) const
    {
      using boost::get;
      //intersection with implicit domain
      Implicit_domain::Intersection implicit_inter =
        r_domain_.implicit_domain.construct_intersection_object()(query);
      //if found, return it
      if(get<2>(implicit_inter) != 0) {
        return Intersection(get<0>(implicit_inter), 2, get<2>(implicit_inter));
      }
      //intersection with polyhedral domain
      Polyhedron_domain::Intersection polyhedron_inter =
        r_domain_.polyhedron_domain.construct_intersection_object()(query);
      //if found, return it
      if(get<2>(polyhedron_inter) != 0) {
        const Point_3 inter_point = get<0>(polyhedron_inter);
        if(!r_domain_.implicit_domain.is_in_domain_object()(inter_point)) {
          return Intersection(inter_point, 1, get<2>(polyhedron_inter));
        }
      }
      //no intersection found
      return Intersection();
    }
  private:
    const Hybrid_domain& r_domain_;
  }; // end Construct_intersection
  Construct_intersection construct_intersection_object() const
  {
    return Construct_intersection(*this);
  }
  //Index types converters
  Index index_from_surface_patch_index(const Surface_patch_index& index) const
  { return index; }
  Index index_from_subdomain_index(const Subdomain_index& index) const
  { return index; }
  Surface_patch_index surface_patch_index(const Index& index) const
  { return index; }
  Subdomain_index subdomain_index(const Index& index) const
  { return index; }
}; // end class Hybrid_domain

typedef CGAL::Mesh_domain_with_polyline_features_3<Hybrid_domain> Domain;
// Triangulation
typedef CGAL::Mesh_triangulation_3<Domain, K>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;
// Criteria
typedef CGAL::Mesh_criteria_3<Tr>     Mesh_criteria;
typedef Mesh_criteria::Edge_criteria  Edge_criteria;
typedef Mesh_criteria::Facet_criteria Facet_criteria;
typedef Mesh_criteria::Cell_criteria  Cell_criteria;
// Function
FT sphere_centered_at_111 (const Point& p)
{
  const FT dx=p.x()-1;
  const FT dy=p.y()-1;
  const FT dz=p.z()-1;

  return dx*dx+dy*dy+dz*dz-1;
}


// This method is fully templated but probably only works for 2d triangulations
// It's adapted from Triangulation_data_structure::write_full_cells
// and from Triangulation_data_structure::operator<<
// Assumes there are no dangling edges in triangulation, only vertices and faces
// template<class Vb, class Fcb>
template<class Tr>
std::ostream &
write_tria_to_off(  std::ostream & os,
                    Tr             triangulation )
                    // CGAL::Triangulation_data_structure_2<Vb, Fcb> tds )
{

  typedef typename Tr::Triangulation_data_structure TDS;
  typedef typename Tr::Vertex_handle                Vertex_handle;
  typedef typename Tr::Finite_vertices_iterator     Vertex_iterator;

  TDS tds = triangulation.tds();

  os << "OFF" << std::endl;

  // outputs the number of vertices and faces
  std::size_t num_verts = triangulation.number_of_vertices();
  std::size_t num_faces = triangulation.number_of_faces();
  std::size_t num_edges = 0;                          //Assumption

  os << num_verts << " " << num_faces << " " << num_edges << std::endl;

  // write the vertices
  std::map<Vertex_handle, int> index_of_vertex;

  // require signed index to handle -1 as index for inifinite vert
  int idx = 0;
  // Vertex_handle inf_v;
  for( Vertex_iterator it = triangulation.finite_vertices_begin(); it != triangulation.finite_vertices_end(); ++it, ++idx )
  {
      os << *it << std::endl;
      index_of_vertex[it] = idx;
  }
  CGAL_assertion( (std::size_t) idx == num_verts );

    const int cur_dim = 2;
    // write the vertex indices of each full_cell
    int i = 0;
    for( auto it = triangulation.finite_faces_begin(); it != triangulation.finite_faces_end(); ++it, ++i )
    {
        os << 3;
        for( int j = 0; j <= cur_dim; ++j )
        {
          assert( index_of_vertex[it->vertex(j)] != -1 ); // assert vertex is reasonable
          os << ' ' << index_of_vertex[it->vertex(j)];
        }
        os << std::endl;
    }

    CGAL_assertion( (std::size_t) i == num_faces );

    return os;
}

// A modifier creating a triangle with the incremental builder.
template <class HDS>
class Build_poly : public CGAL::Modifier_base<HDS> {
public:
  typedef typename HDS::Vertex       PolyVertex;
  typedef typename PolyVertex::Point PolyPoint;
    std::vector< PolyPoint >    m_verts;
    std::vector< unsigned int > m_faces;
    Build_poly( std::vector< PolyPoint > i_verts, std::vector< unsigned int > i_faces ) :
      m_verts( i_verts ),
      m_faces( i_faces )
    { }
    void operator()( HDS& hds) {
        // Postcondition: hds is a valid polyhedral surface.
        CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);
        B.begin_surface( m_verts.size(), m_faces.size() / 3 );

        for( size_t j = 0; j < m_verts.size(); j++ ){
          B.add_vertex( m_verts[ j ] );
        }

        for( size_t i=0; i < m_faces.size(); i+=3 ){
          B.begin_facet();
          B.add_vertex_to_facet( m_faces[ i ] );
          B.add_vertex_to_facet( m_faces[ i+1 ] );
          B.add_vertex_to_facet( m_faces[ i+2 ] );
          B.end_facet();
        }

        B.end_surface();
    }
};

template<class Polyhedron>
Polyhedron build_bbox( Polyhedron& P, double xMin, double xMax, double yMin, double yMax, double zMin, double zMax )
{
  typedef typename Polyhedron::HalfedgeDS   HDS;
  typedef typename HDS::Vertex::Point       PolyPoint;
  PolyPoint x0y0z0( xMin, yMin, zMin );
  PolyPoint x0y0z1( xMin, yMin, zMax );
  PolyPoint x0y1z0( xMin, yMax, zMin );
  PolyPoint x0y1z1( xMin, yMax, zMax );
  PolyPoint x1y0z0( xMax, yMin, zMin );
  PolyPoint x1y0z1( xMax, yMin, zMax );
  PolyPoint x1y1z0( xMax, yMax, zMin );
  PolyPoint x1y1z1( xMax, yMax, zMax );

  std::vector< PolyPoint > verts;
  std::vector< unsigned int > faces;

  verts.push_back( x0y0z0 );
  verts.push_back( x0y0z1 );
  verts.push_back( x0y1z0 );
  verts.push_back( x0y1z1 );
  verts.push_back( x1y0z0 );
  verts.push_back( x1y0z1 );
  verts.push_back( x1y1z0 );
  verts.push_back( x1y1z1 );

  faces = { 0, 2, 6,
            4, 0, 6,
            3, 1, 7,
            7, 1, 5,
            0, 4, 5,
            5, 1, 0,
            6, 2, 7,
            7, 2, 3,
            2, 0, 3,
            3, 0, 1,
            6, 7, 4,
            7, 5, 4 };

  Build_poly<HDS> bbox( verts, faces );

  P.delegate( bbox );
  // Build_triangle<HDS> xyMin1( x0y0z0, x0y1z0, x1y1z0 );
  // Build_triangle<HDS> xyMin2( x1y0z0, x0y0z0, x1y1z0 );
  // Build_triangle<HDS> xyMax1( x0y0z1, x0y1z1, x1y1z1 );
  // Build_triangle<HDS> xyMax2( x0y0z1, x1y1z1, x1y0z1 );
  //
  // std::cout << "XZ";
  // Build_triangle<HDS> xzMin1( x0y0z0, x1y0z0, x1y0z1 );
  // Build_triangle<HDS> xzMin2( x0y0z0, x1y0z1, x0y0z1 );
  // Build_triangle<HDS> xzMax1( x1y1z0 ,x0y1z0, x1y1z1 );
  // Build_triangle<HDS> xzMax2( x1y1z1, x0y1z0, x0y1z1 );
  //
  // std::cout << "YZ" << std::endl;
  // Build_triangle<HDS> yzMin1( x0y1z0, x0y0z0, x0y1z1 );
  // Build_triangle<HDS> yzMin2( x0y1z1, x0y0z0, x0y0z1 );
  // Build_triangle<HDS> yzMax1( x1y1z0, x1y0z0, x1y1z1 );
  // Build_triangle<HDS> yzMax2( x1y0z1, x1y1z1, x1y0z0 );
  // std::cout << "done building tris" << std::endl;
  //
  // P.delegate( xyMin1 ); std::cout << "a";
  // P.delegate( xyMin2 ); std::cout << "b";
  // P.delegate( xyMax1 ); std::cout << "c";
  // P.delegate( xyMax2 ); std::cout << "d";
  // P.delegate( xzMin1 ); std::cout << "e";
  // P.delegate( xzMin2 ); std::cout << "f";
  // P.delegate( xzMax1 ); std::cout << "g";
  // P.delegate( xzMax2 ); std::cout << "h";
  // P.delegate( yzMin1 ); std::cout << "i";
  // P.delegate( yzMin2 ); std::cout << "j";
  // P.delegate( yzMax1 ); std::cout << "k";
  // P.delegate( yzMax2 );

  return P;
}

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;
int main()
{
  // const char* fname = "cube.off";
  // // Create input polyhedron
  // Polyhedron polyhedron;
  // std::ifstream input(fname);
  // input >> polyhedron;
  // if(input.bad()){
  //   std::cerr << "Error: Cannot read file " <<  fname << std::endl;
  //   return EXIT_FAILURE;
  // }
  // input.close();

  // Constants defining bounding box
  const double X_MIN =  330000;
  const double X_MAX =  500000;
  const double Y_MIN =  -10000;
  const double Y_MAX =   25000;
  const double Z_MIN = 3700000;
  const double Z_MAX = 3850000;

  EDGE_LOG_INFO << "Building bounding box...";
  Polyhedron poly_bbox;
  build_bbox( poly_bbox, X_MIN, X_MAX, Y_MIN, Y_MAX, Z_MIN, Z_MAX );
  EDGE_LOG_INFO << "Done.";

  if ( !poly_bbox.is_valid() )
    EDGE_LOG_INFO << "WARNING: bounding box is not valid";

  EDGE_LOG_INFO << "Printing bounding box...";
  EDGE_LOG_INFO << poly_bbox;


  edge_cut::surf::Topo l_topo("/home/david/work/edge/workspace/edge_cut/data/map_proj.xyz");
  std::ofstream outstream("topoTria.off");
  write_tria_to_off( outstream, l_topo.m_delTria );
  outstream.close();

  std::cout << "wrote topography" << std::endl;
  std::cin.get();
  // Use scan_OFF instead of >> to allow for verbose outputs
  // see CGAL/IO/scan_OFF.h for more info
  Polyhedron polyhedron;
  std::ifstream instream("topoTria.off");
  EDGE_LOG_INFO << "Building topography polyhedron from .off file...";
  CGAL::scan_OFF( instream, polyhedron, true );
  // instream >> polyhedron;
  instream.close();

  // Polyhedra domain
  EDGE_LOG_INFO << "Constructing polyhedral domain...";
  Polyhedron_domain polyhedron_domain(polyhedron, poly_bbox);
  EDGE_LOG_INFO << "Done.";

  // Implicit domain
  // - the first argument is the function pointer
  // - the second argument is a bounding sphere of the domain
  // (Warning: Sphere_3 constructor uses square radius !)
  Function_vector l_halfSpaces;
  Function l_hSpace1(&posHSpaceDisp);
  Function l_hSpace2(&negHSpaceDisp);
  l_halfSpaces.push_back( l_hSpace1 );
  l_halfSpaces.push_back( l_hSpace2 );
  //
  K::Point_3  l_bndSphereCenter(415000, 10000, 3775000);
  K::FT       l_bndSphereRadiusSquared = 150000.0*150000.0;
  K::Sphere_3 l_bndSphere( l_bndSphereCenter, l_bndSphereRadiusSquared);
  //
  // std::cout << "Constructing mesh domain" << std::endl;

  std::vector<std::string> vps;
  vps.push_back("+-");
  vps.push_back("-+");

  Implicit_domain impl_domain( Function_wrapper( l_halfSpaces, vps ), l_bndSphere, 1e-3);
  // Implicit_domain impl_domain( l_hSpace1, l_bndSphere );
  // Implicit_domain impl_domain(sphere_centered_at_111,
                                // K::Sphere_3(K::Point_3(1, 1, 1), 2.));

  // Full domain
  Domain domain(impl_domain, polyhedron_domain);
  // Domain domain(hSpaces_domain, polyhedron_domain);

  // Read polyline features
  // NOTE: polylines must intersect at a shared vertex or be disjoint!
  const char* lines_fname = "hybrid_topo.polylines.txt";
  std::vector<std::vector<Point> > featured_curves;
  if (!read_polylines(lines_fname, featured_curves))
  { // see file "read_polylines.h"
    std::cerr << "Error: Cannot read file " << lines_fname << std::endl;
    return EXIT_FAILURE;
  }
  // Add features for protection
  domain.add_features(featured_curves.begin(), featured_curves.end());

  // Criteria
  double l_edge_length, l_angle_bound, l_facet_size, l_facet_app, l_re_ratio, l_cell_size;

  std::cout << "edge length: ";
  std::cin >> l_edge_length;
  // std::cout << std::endl;

  std::cout << "angle bound: ";
  std::cin >> l_angle_bound;
  // std::cout << std::endl;

  std::cout << "face size: ";
  std::cin >> l_facet_size;
  // std::cout << std::endl;

  std::cout << "facet approx: ";
  std::cin >> l_facet_app;
  // std::cout << std::endl;

  std::cout << "radius-edge ratio: ";
  std::cin >> l_re_ratio;
  // std::cout << std::endl;

  std::cout << "cell size: ";
  std::cin >> l_cell_size;
  // std::cout << std::endl;

  Edge_criteria edge_criteria( l_edge_length );
  Facet_criteria facet_criteria( l_angle_bound, l_facet_size, l_facet_app ); // angle, size, approximation
  Cell_criteria cell_criteria( l_re_ratio, l_cell_size ); // radius-edge ratio, size
  Mesh_criteria criteria(edge_criteria, facet_criteria, cell_criteria);

  // Mesh generation (without optimization)
  EDGE_LOG_INFO << "Starting up mesher...";
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria,
                                      no_perturb(), no_exude());
  EDGE_LOG_INFO << "Done!";
  EDGE_LOG_INFO << "Just as a reminder, refinement parameters were:";
  EDGE_LOG_INFO << "Edge Length:        "  << l_edge_length;
  EDGE_LOG_INFO << "Angle Bound:        "  << l_angle_bound;
  EDGE_LOG_INFO << "Facet Size:         "  << l_facet_size;
  EDGE_LOG_INFO << "Facet Approximation: " << l_facet_app;
  EDGE_LOG_INFO << "Radius-Edge Ratio:  "  << l_re_ratio;
  EDGE_LOG_INFO << "Cell Size:          "  << l_cell_size;

  // Output
  // dump_c3t3(c3t3, "out");
  std::string out_file_name = "hybrid.off";
  EDGE_LOG_INFO << "Writing mesh to " << out_file_name << " ...";

  std::ofstream off_file( out_file_name );
  c3t3.output_facets_in_complex_to_off( off_file );
  EDGE_LOG_INFO << "Done writing mesh.";

  return EXIT_SUCCESS;
}
