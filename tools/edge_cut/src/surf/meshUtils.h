#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_slicer.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/internal/Mesh_3/Boundary_of_subdomain_of_complex_3_in_triangulation_3_to_off.h>
#include <boost/optional/optional_io.hpp> // TODO Check if this is still needed


#include "implicit_functions.h" // Definition of Function_wrapper

namespace edge_cut{
  namespace surf{
    // Typedefs for Domain
    typedef CGAL::Exact_predicates_inexact_constructions_kernel                 K;
    typedef CGAL::Mesh_polyhedron_3<K>::type                                    Polyhedron;
    typedef CGAL::Polyhedral_mesh_domain_with_features_3<K>                     Poly_domain;
    typedef CGAL::Mesh_domain_with_polyline_features_3<Poly_domain>             Mesh_domain;

    // Typedefs for Triangulation
    typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type                       Mesh_tria;
    typedef CGAL::Mesh_complex_3_in_triangulation_3<Mesh_tria>                  C3t3;

    // Typedefs for Mesh Criteria
    typedef CGAL::Mesh_criteria_3<Mesh_tria>                                    Mesh_criteria;
    typedef Mesh_criteria::Edge_criteria                                        Edge_criteria;
    typedef Mesh_criteria::Facet_criteria                                       Facet_criteria;
    typedef Mesh_criteria::Cell_criteria                                        Cell_criteria;

    // Typedefs for 1D Features
    typedef CGAL::Surface_mesh<K::Point_3>                                      Surf_mesh;
    typedef CGAL::Polygon_mesh_slicer<Polyhedron, K>                            Poly_slicer;
    typedef std::vector<K::Point_3>                                             Polyline_type;
    typedef std::list< Polyline_type >                                          Polylines;
    typedef std::vector< std::size_t >                                          Polygon;


    /**
     * Sizing field
     * The vector of depth values, "m_layers," must be in decreasing order
     *
     *     baseValue
     * ----------------- layers[0]
     *     scales[0]
     * ----------------- layers[1]
     *     scales[1]
     * ----------------- layers[2]
     *       ...
     **/
    struct SizingField {
      typedef K::FT FT;
      typedef K::Point_3 Point_3;
      typedef Mesh_domain::Index Index;

      FT                m_innerVal;
      FT                m_scale;
      Point_3           m_center;
      FT                m_innerRad;
      FT                m_outerRad;

      SizingField( FT i_innerVal, FT i_scale, Point_3 i_center, FT i_innerRad, FT i_outerRad ) :
        m_innerVal( i_innerVal ),
        m_scale( i_scale ),
        m_center( i_center ),
        m_innerRad( i_innerRad ),
        m_outerRad( i_outerRad )
      { }

      FT operator()(const Point_3& p, const int, const Index&) const;
    };

    /**
     * Computes the intersection between the polygonal surface stored in i_slicer and
     * the plane i_plane, which should be a collection of piecewise-linear curves.
     * Prints an error and returns empty if the intersection contains two or
     * more distinct curves
     *
     * @param i_slicer A representation of a polygonal surface that computes intersections
     * @param i_plane The plane to be intersected with
     * @return The single piecewise-linear intersection, or nothing if more than one
     * intersection is found
     **/
    Polyline_type topoIntersect( Poly_slicer& i_slicer, K::Plane_3 i_plane );

    /**
     * Checks monotonicity of a polyline in the (i_n)th coordinate
     *
     * @param i_p Polyline to be tested
     * @param i_n Index of the cartesian coordinate used to order vertices
     * @param i_inc Flag - true if testing for increasing, false for decreasing
     **/
    bool checkMonotonic( Polyline_type & i_p, unsigned int i_n, bool i_inc );

    /**
     * Reorders i_p so that it is increasing with respect to the (i_n)th coordinate.
     *
     * @pre i_p must be monotonically increasing or decreasing
     *
     * @param i_p The polyline to be ordered
     * @param i_n Index of the cartesian coordinate used to order vertices
     * @return Nothing
     **/
    void orderPolyline( Polyline_type & i_p, unsigned int i_n );

    /**
     * Construct polylines representing the bounding box of the computational region
     * as well as any intersections between topography, the fault, and the bounding box
     *
     * @param i_topoSurface Representation of the topography
     * @param i_xMin Minimum x coordinate of bounding box
     * @param i_xMax Maximum x coordinate of bounding box
     * @param i_yMin Minimum y coordinate of bounding box
     * @param i_yMax Maximum y coordinate of bounding box
     * @param i_zMin Minimum z coordinate of bounding box
     * @param i_xMax Maximum z coordinate of bounding box
     * @return Collection of 1D features
     **/
    // std::list< Polyline_type > build1DFeatures( Surf_mesh&  i_topoSurface,
    //                                             double      i_xMin,
    //                                             double      i_xMax,
    //                                             double      i_yMin,
    //                                             double      i_yMax,
    //                                             double      i_zMin,
    //                                             double      i_zMax         );

    /**
     * Creates a basic polyhedral model of the free surface boundaries of
     * a three-dimensional box with topography.
     **/
    Polyhedron& makeFreeSurfBdry ( Polyhedron        & io_bdry,
                                   Topo        const & i_topo  );

    Polyhedron& trimFSB (          Polyhedron const & i_topo,
                                   Polyhedron       & io_fsb );

    void c3t3ToPolyhedron(         C3t3       const & c3t3,
                                   Polyhedron       & polyhedron );


    // TODO Currently unneeded, should be removed
    // template perameter "Triangulation" must be model of CGAL::Triangulation_2
    template< class Triangulation >
    std::ostream & writeTria2ToOff( Triangulation& tr, std::ostream & os ) {
      typedef typename Triangulation::Vertex_handle                Vertex_handle;
      typedef typename Triangulation::Finite_vertices_iterator     Vertex_iterator;
      typedef typename Triangulation::Finite_faces_iterator        Face_iterator;

      os << "OFF" << std::endl;

      // outputs the number of vertices and faces
      std::size_t num_verts = tr.number_of_vertices();
      std::size_t num_faces = tr.number_of_faces();
      std::size_t num_edges = 0;                          //Assumption

      os << num_verts << " " << num_faces << " " << num_edges << std::endl;

      // write the vertices
      std::map<Vertex_handle, std::size_t> index_of_vertex;

      std::size_t v_idx = 0;
      for( Vertex_iterator it = tr.finite_vertices_begin(); it != tr.finite_vertices_end(); ++it, ++v_idx )
      {
          os << *it << std::endl;
          index_of_vertex[it] = v_idx;
      }
      CGAL_assertion( v_idx == num_verts );

      // write the vertex indices of each full_cell
      std::size_t f_idx = 0;
      for( Face_iterator it = tr.finite_faces_begin(); it != tr.finite_faces_end(); ++it, ++f_idx )
      {
          os << 3;
          for( int j = 0; j < 3; ++j )
          {
            os << ' ' << index_of_vertex[ it->vertex(j) ];
          }
          os << std::endl;
      }
      CGAL_assertion( f_idx == num_faces );

      return os;
    }



    // template perameter "Triangulation" must be model of CGAL::Triangulation_3
    template< class Triangulation >
    std::ostream & writeTria3ToOff( Triangulation& tr, std::ostream & os ) {
      typedef typename Triangulation::Vertex_handle                Vertex_handle;
      typedef typename Triangulation::Finite_vertices_iterator     Vertex_iterator;
      typedef typename Triangulation::Finite_facets_iterator        Face_iterator;

      os << "OFF" << std::endl;

      // outputs the number of vertices and faces
      std::size_t num_verts = tr.number_of_vertices();
      std::size_t num_faces = tr.number_of_finite_facets();
      std::size_t num_edges = 0;                          //Assumption

      os << num_verts << " " << num_faces << " " << num_edges << std::endl;

      // write the vertices
      std::map<Vertex_handle, std::size_t> index_of_vertex;

      std::size_t v_idx = 0;
      for( Vertex_iterator it = tr.finite_vertices_begin(); it != tr.finite_vertices_end(); ++it, ++v_idx )
      {
          os << it->point().x() << ' '
             << it->point().y() << ' '
             << it->point().z() << std::endl;
          std::cout << it->point().x() << ' '
                    << it->point().y() << ' '
                    << it->point().z() << std::endl; // TODO
          index_of_vertex[it] = v_idx;
      }
      CGAL_assertion( v_idx == num_verts );

      // write the vertex indices of each full_cell
      std::size_t f_idx = 0;
      for( Face_iterator it = tr.finite_facets_begin(); it != tr.finite_facets_end(); ++it, ++f_idx )
      {
          os << 3;
          // it is a (Cell_handle, int) pair, where the cell_handle points to
          // cell incident to the facet, and the int is the local index (relative
          // to the cell) for the vertex not on the facet
          // Hence it->second is the only vertex on the cell which is not on the facet
          for( int j = 0; j < 4; ++j )
          {
            // it->first is a handle to the cell that contains the facet
            if ( j != it->second ){
              os << ' ' << index_of_vertex[ it->first->vertex(j) ];
              std::cout << ' ' << index_of_vertex[ it->first->vertex(j) ]; //  TODO
            }
          }
          os << std::endl;
          std::cout << std::endl; // TODO
      }
      CGAL_assertion( f_idx == num_faces );
      std::cout << "End writeTria3ToOff" << std::endl; //TODO
      return os;
    }
  }
}
