#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Implicit_to_labeling_function_wrapper.h>
#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/Mesh_domain_with_polyline_features_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_slicer.h>
#include <boost/optional/optional_io.hpp> // TODO Check if this is still needed

#include "implicit_functions.h" // Definition of Function_wrapper

namespace edge_cut{
  namespace surf{
    // Typedefs for Domain
    typedef CGAL::Exact_predicates_inexact_constructions_kernel                 K;
    typedef FT_to_point_function_wrapper<K::FT, K::Point_3>                     Function;
    typedef CGAL::Implicit_multi_domain_to_labeling_function_wrapper<Function>  Function_wrapper;
    typedef Function_wrapper::Function_vector                                   Function_vector;
    typedef CGAL::Labeled_mesh_domain_3<Function_wrapper, K>                    Labeled_domain;
    typedef CGAL::Mesh_domain_with_polyline_features_3<Labeled_domain>          Mesh_domain;

    // Typedefs for Triangulation
    typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type       Mesh_tria;
    typedef CGAL::Mesh_complex_3_in_triangulation_3<Mesh_tria>  C3t3;

    // Typedefs for Mesh Criteria
    typedef CGAL::Mesh_criteria_3<Mesh_tria>        Mesh_criteria;
    typedef Mesh_criteria::Edge_criteria            Edge_criteria;
    typedef Mesh_criteria::Facet_criteria           Facet_criteria;
    typedef Mesh_criteria::Cell_criteria            Cell_criteria;

    // Typedefs for 1D Features
    typedef CGAL::Surface_mesh<K::Point_3>                    Surf_mesh;
    typedef CGAL::Polygon_mesh_slicer<Surf_mesh, K>           Poly_slicer;
    typedef std::vector<K::Point_3>                           Polyline_type;
    typedef std::list< Polyline_type >                        Polylines;
    typedef std::vector< std::size_t >                        Polygon;


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


    Polyline_type topoIntersect( Poly_slicer& slicer, K::Plane_3 plane );
    bool checkMonotonic( Polyline_type & i_p, unsigned int i_n, bool i_inc );
    void orderPolyline( Polyline_type & i_p, unsigned int i_n );

  }
}
