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
template< class K >
edge_cut::io::Config< K >::Config( std::string i_xmlPath )
{
  // Check for valid config
  pugi::xml_parse_result l_parseResult = m_doc.load_file( i_xmlPath.c_str() );

  // Print any error info
  if( !l_parseResult ) {
    EDGE_LOG_ERROR << "XML [" << i_xmlPath << "] parsed with errors, attr value: ["
                   << m_doc.child("node").attribute("attr").value() << "]";
    EDGE_LOG_FATAL << "Description: " << l_parseResult.description()
                   << " (error at ["   << i_xmlPath << " " << l_parseResult.offset << "]";
  }

  pugi::xml_node l_root = m_doc.child("edge_cut");

  // Computational Domain
  m_bBox[0] = l_root.child("bounding_box").child("xMin").text().as_double();
  m_bBox[1] = l_root.child("bounding_box").child("xMax").text().as_double();
  m_bBox[2] = l_root.child("bounding_box").child("yMin").text().as_double();
  m_bBox[3] = l_root.child("bounding_box").child("yMax").text().as_double();
  m_bBox[4] = l_root.child("bounding_box").child("zMin").text().as_double();
  m_bBox[5] = l_root.child("bounding_box").child("zMax").text().as_double();

  // IO Targets
  m_topoIn  = l_root.child("io").child("topography").child("in").text().as_string();
  m_topoOut = l_root.child("io").child("topography").child("out").text().as_string();
  m_bdryIn  = l_root.child("io").child("boundary").child("in").text().as_string();
  m_bdryOut = l_root.child("io").child("boundary").child("out").text().as_string();

  // Refinement Region Definition
  m_scale    = l_root.child("region").child("scale").text().as_double();
  m_innerRad = l_root.child("region").child("inner_radius").text().as_double();
  m_outerRad = l_root.child("region").child("outer_radius").text().as_double();
  m_center   = typename K::Point_3( l_root.child("region").child("center").child("x").text().as_double(),
                                    l_root.child("region").child("center").child("y").text().as_double(),
                                    l_root.child("region").child("center").child("z").text().as_double() );

  // Meshing Criteria
  m_edgeBase        = l_root.child("refine").child("edge").text().as_double();
  m_facetSizeBase   = l_root.child("refine").child("facet").text().as_double();
  m_facetApproxBase = l_root.child("refine").child("approx").text().as_double();
  m_angleBound      = l_root.child("refine").child("angle").text().as_double();
}
