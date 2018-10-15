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
#ifndef EDGE_CUT_MESH_UTILS_H_
#define EDGE_CUT_MESH_UTILS_H_

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
#include <CGAL/Polygon_mesh_processing/connected_components.h>

#include <boost/foreach.hpp>
#include <boost/property_map/property_map.hpp>

#include "io/logging.hpp"
#include "surf/Topo.h"

namespace edge_cut{
  namespace surf{
    // Typedefs for Domain
    typedef CGAL::Exact_predicates_inexact_constructions_kernel                 K;
    typedef CGAL::Surface_mesh<K::Point_3>                                      Surf_mesh;
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
    typedef CGAL::Polygon_mesh_slicer<Polyhedron, K>                            Poly_slicer;
    typedef std::vector<K::Point_3>                                             Polyline_type;
    typedef std::list< Polyline_type >                                          Polylines;
    typedef std::vector< std::size_t >                                          Polygon;


    struct SizingField {
      typedef K::FT FT;
      typedef K::Point_3 Point_3;
      typedef Mesh_domain::Index Index;

      FT                m_innerVal;
      FT                m_scale;
      Point_3           m_center;
      FT                m_innerRad;
      FT                m_outerRad;
      FT                m_widthInv;

      SizingField( FT i_innerVal, FT i_scale, Point_3 i_center, FT i_innerRad, FT i_outerRad ) :
        m_innerVal( i_innerVal ),
        m_scale( i_scale ),
        m_center( i_center ),
        m_innerRad( i_innerRad ),
        m_outerRad( i_outerRad )
      {
        if ( m_innerRad != m_outerRad )
          m_widthInv = 1 / ( m_outerRad - m_innerRad );
        else
          m_widthInv = 0;
      }

      FT operator()( const Point_3& p, const int, const Index& ) const;
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
     * Computes the intersection of the four bounding planes with the topography
     * and generates line segments spanning the bounding box
     *
     * @param i_topoSurface the topographical surface for which we are computing the boundary
     * @param i_bBox array of doubles specifying the bounds in the x,y,z coordinates
     *
     * @return A list of polylines (AKA vectors of points) which contain all the
     *         points on the boundary of the topography mesh (points may
     *         appear more than once), as well as polylines for each edge of the
     *         bounding box (below the topography)
     **/
    std::list< Polyline_type > getIntersectionFeatures( const Polyhedron& i_topoSurface, double const * i_bBox );


    /**
     * Creates a basic polyhedral model of the boundaries of a three-
     * dimensional box. Overwrites input polyhedron.
     *
     * @param io_bdry reference to polyhedron which will contain boundary
     * @param i_bBox  pointer to array of 6 double defining the extent of the box
     *                array is ordered as { xMin, xMax, yMin, yMax, zMin, zMax }
     *
     * @return A reference to the polyhedron with the constructed boundary mesh
     **/
    Polyhedron& makeBdry ( Polyhedron        & io_bdry,
                           double     const  * i_bBox  );


    /**
     * Converts the triangulation in a CGAL::Complex_3_in_Triangulation_3 to a
     * CGAL::Polyhedron_3
     *
     * @param c3t3 the complex-in-triangulation to convert
     * @param polyhedron a reference to the polyhedron to be created (overwrites any existing data)
     *
     * @return None
     **/
    void c3t3ToPolyhedron(  C3t3       const & c3t3,
                            Polyhedron       & polyhedron );

    /**
     * Converts the triangulation in a CGAL::Complex_3_in_Triangulation_3 to a
     * CGAL::Surface_mesh
     *
     * @param c3t3 reference to the complex-in-triangulation to convert
     * @param io_surfMesh a reference to the surfaceMesh to be created (overwrites any existing data)
     *
     * @return None
     **/
    void c3t3ToSurfMesh(  C3t3      const & i_c3t3,
                          Surf_mesh       & io_surfMesh );


    /**
     * @param io_topoPolyMesh Reference to a CGAL::Polyhedron_3 which will
     *        contain the polyhedral mesh. This Polyhedron_3 must be empty when
     *        passed, as the function will be appending to it.
     * @param i_topoFile ASCII file containing point-set data to be triangulated
     *
     * @return A reference to the constructed topography mesh
     **/
     Polyhedron& topoPolyMeshFromXYZ( Polyhedron& io_topoPolyMesh, std::string const & i_topoFile );
  }
}

#endif
