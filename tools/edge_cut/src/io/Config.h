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
 * Configuration class for EDGEcut.
 **/
#ifndef EDGE_CUT_CONFIG_H
#define EDGE_CUT_CONFIG_H

#include <string>
#include <CGAL/make_mesh_3.h>
#include "logging.hpp"
#include "../../../submodules/pugixml/src/pugixml.hpp"

namespace edge_cut {
  namespace io {
    template< class K > class Config;
  }
}

template< class K >
class edge_cut::io::Config {
public:
  Config( std::string i_xmlPath );

  void printConfig();

  // Path to XML file containing the runtime config
  std::string m_xmlPath;

  // XML document containing configuration
  pugi::xml_document m_doc;

  // Coordinates of bounding box
  double m_bBox[6];

  // Input file with topography data
  std::string m_topoIn;

  // Output file for topography surface mesh
  std::string m_topoOut;

  // Input file with triangular mesh model of boundary
  std::string m_bdryIn;

  // Output file for domain boundary surface mesh
  std::string m_bdryOut;

  // Depth layers and scaling factors for depth-based refinement
  std::map< typename K::FT, typename K::FT, std::greater<typename K::FT> > m_layers;

  // Radius for region of maximum refinement
  typename K::FT m_innerRad;

  // Radius outside of which there is minimum refinement
  typename K::FT m_outerRad;

  // Scale factor of refinement level from innerRad -> outerRad
  typename K::FT m_scale;

  // Center of refinement regions
  typename K::Point_3 m_center;

  // Target edge length (for 1D features ONLY) for mesher
  typename K::FT m_edgeBase;

  // Target facet size for mesher
  typename K::FT m_facetSizeBase;

  // Target surface approximation for mesher
  typename K::FT m_facetApproxBase;

  // Minimum allowable facet angle for mesher
  typename K::FT m_angleBound;

  // Time limit for each mesh optimization step (0 == no time limit)
  double m_optimizeTime;

  // CGAL object containing options for lloyd smoothing in mesher
  CGAL::parameters::internal::Lloyd_options   m_lloydOpts;

  // CGAL object containing options for ODT smoothing in mesher
  CGAL::parameters::internal::Odt_options     m_odtOpts;

  // CGAL object containing options for sliver perturbation in mesher
  CGAL::parameters::internal::Perturb_options m_perturbOpts;

  // CGAL object containing options for sliver exudation in mesher
  CGAL::parameters::internal::Exude_options   m_exudeOpts;
};

#include "Config.inc"

#endif
