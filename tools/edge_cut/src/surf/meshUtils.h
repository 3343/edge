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
    struct DepthSizingField {
      typedef K::FT FT;
      typedef K::Point_3 Point_3;
      typedef Mesh_domain::Index Index;

      FT                m_baseValue;
      std::vector< FT > m_layers;
      std::vector< FT > m_scales;

      DepthSizingField( FT i_baseValue, std::vector< FT > i_layers, std::vector< FT > i_scales ) :
        m_baseValue( i_baseValue ),
        m_layers( i_layers ),
        m_scales( i_scales )
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
    Polyhedron& makeFreeSurfBdry( Polyhedron& i_topo );

  }
}
