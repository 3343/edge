/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
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
 * Meshing of topographic data.
 **/

#ifndef EDGE_CUT_TOPO_H_
#define EDGE_CUT_TOPO_H_

#include <string>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Triangulation_hierarchy_2.h>

namespace edge_cut {
  namespace surf {
    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
    typedef CGAL::Projection_traits_xy_3<K>                     Gt;
    typedef CGAL::Triangulation_vertex_base_2<Gt>               Vbb;
    typedef CGAL::Triangulation_hierarchy_vertex_base_2<Vbb>    Vb;
    typedef CGAL::Triangulation_face_base_2<Gt>                 Fb;
    typedef CGAL::Triangulation_data_structure_2<Vb,Fb>         Tds;

    typedef CGAL::Delaunay_triangulation_2<Gt,Tds>              Dt;
    typedef CGAL::Triangulation_hierarchy_2<Dt>                 Triangulation;
    typedef Triangulation::Point                                TopoPoint;

    class Topo;
  }
}

const unsigned short C_MAX_SURF_INTER = 10;

class edge_cut::surf::Topo {
public:
    // 2.5D delaunay triangulation of the topographic data
    Triangulation* m_delTria;


    /**
     * Constructor: Computes the 2.5D Delaunay triangulation from the given topographic data.
     *
     * @param i_topoFile location of the topographic data.
     **/
    Topo( std::string const & i_topoFile );

    /**
     * Destructor: Deletes memory allocated to Delaunay Triangulation
     **/
    ~Topo();

    // Delete copy constructor and copy assignment
    Topo( const Topo& ) = delete;
    Topo& operator=( const Topo& ) = delete;

    /**
     * Writes the topo mesh to the provided stream in .OFF format
     *
     * @param i_os ostream to be written to
     * @return the ostream after writing
     **/
    std::ostream & writeTriaToOff( std::ostream & io_os ) const;
};

#endif
