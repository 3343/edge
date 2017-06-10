/**
 * @file This file is part of EDGE.
 *
 * @author Josh Tobin (tobinrj AT tcd.ie)
 *
 * @section LICENSE
 * Copyright (c) 2017, Regents of the University of California
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
 * Meshing criteria that supports multple layers.
 **/

#ifndef LAYEREDCRITERIA_H_
#define LAYEREDCRITERIA_H_

#include <vector>

#include <CGAL/Surface_mesher/Standard_criteria.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Robust_circumcenter_traits_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Surface_mesh_vertex_base_3.h>
#include <CGAL/Surface_mesh_cell_base_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_with_circumcenter_3.h>
#include <CGAL/tags.h>


namespace edge_cut {
  namespace surf {
    template<class Tr> class LayeredCriteria;
  }
}
  
template <class Tr>
class edge_cut::surf::LayeredCriteria
{
  typedef CGAL::Surface_mesher::Refine_criterion<Tr> Criterion;
  typedef CGAL::Surface_mesher::Standard_criteria<Criterion> Criteria;

public:
  typedef Tr Triangulation;
  typedef typename Tr::Geom_traits::FT FT;

  typedef typename Criteria::Quality Quality;
  typedef typename Tr::Facet Facet;
  typedef typename Tr::Point Point;

  /**
   * Constructor: Creates a default mesh criteria for each layer of the input.  If the number of layers in the depth
   * vector is different from the number of layers in the scale vector, the smaller of the two dimensions is used
   * and the other layers are ignored.
   *
   * @param i_angle_bound        the min angle allowed in a facet.
   * @param i_radius_bound_base  the base max size of radii of Delauney ball.
   * @param i_distance_bound_base  the base max size of distance between circumcenter of mesh facet and center of Delauney ball.
   * @param i_depths  vector of depths of layer interfaces, in decreasing order.
   * @param i_scale   the scale to apply to each layer.
   **/ 
  LayeredCriteria(const FT i_angle_bound,
	          const FT i_radius_bound_base,
	          const FT i_distance_bound_base,
		  std::vector<double> i_depths,
		  std::vector<FT> i_scale
    )
    : m_depths(i_depths)
  {    
    for(unsigned int l_layer=0; l_layer<m_depths.size() && l_layer<i_scale.size(); l_layer++)
    {
      m_criteria.push_back(new CGAL::Surface_mesh_default_criteria_3<Tr>(i_angle_bound, i_radius_bound_base * i_scale[l_layer], i_distance_bound_base * i_scale[l_layer]));
    }
  }

  /**
   * Destructor: Delete the criteria that were allocated in the constructor.
   **/ 
  
  ~LayeredCriteria()
  {
    for(auto l_criteriaPt : m_criteria)
    {
      delete l_criteriaPt;
    }
  }

  /**
   * Constructor: Computes the quality of a given facet.
   *
   * @param i_face    the face whose quality is to be judged.
   * @param o_quality the quality of the inputted face.
   **/   
  bool is_bad (const Facet& i_face, Quality& o_quality) const
  {
    Point l_point1 = i_face.first->vertex((i_face.second+1)&3)->point();
    Point l_point2 = i_face.first->vertex((i_face.second+2)&3)->point();
    Point l_point3 = i_face.first->vertex((i_face.second+3)&3)->point();


    auto l_avg_depth = l_point1.z() + l_point2.z() + l_point3.z();
    l_avg_depth /= 3.;

    for(unsigned int l_layer=0; l_layer<m_depths.size(); l_layer++)
    {
      if(l_avg_depth > m_depths[l_layer])
      {
	return m_criteria[l_layer]->is_bad(i_face, o_quality);
      }
    }

    // This point is below all of the given interface layers, so use the bottom-most layer
    if(m_depths.size() > 0)
      return m_criteria[m_depths.size()-1]->is_bad(i_face,o_quality);

    // If there are no layers at all, any mesh is fine.  
    return false;
  }

private:
  std::vector<CGAL::Surface_mesh_default_criteria_3<Tr>*> m_criteria;
  std::vector<double> m_depths;
};

#endif
