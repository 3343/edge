/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2020, Friedrich Schiller University Jena
 * Copyright (c) 2019-2020, Alexander Breuer
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
 * EDGE-V config.
 **/
#include "Config.h"

#include "logging.h"
#define PUGIXML_HEADER_ONLY
#include <pugixml.hpp>
#undef PUGIXML_HEADER_ONLY

edge_v::io::Config::Config( std::string & i_xml ) {
  // parse the XML-config
  pugi::xml_document l_doc;
  l_doc.load_file( i_xml.c_str() );

  // read input and ouput mesh
  pugi::xml_node l_mesh = l_doc.child("edge_v").child("mesh");
  m_meshIn   = l_mesh.child("files").child("in").text().as_string();
  m_meshOut[0] = l_mesh.child("files").child("out").child("base").text().as_string();
  m_meshOut[1] = l_mesh.child("files").child("out").child("extension").text().as_string();
  m_writeElAn = l_mesh.child("write_element_annotations").text().as_bool();
  m_periodic = l_mesh.child("periodic").text().as_int();
  m_nPartitions = l_mesh.child("n_partitions").text().as_ullong();
  m_nPartitions = std::max( m_nPartitions, std::size_t(1) );

  // read velocity model
  pugi::xml_node l_velMod = l_doc.child("edge_v").child("velocity_model");
  m_seismicExpr = l_velMod.child("seismic_expression").text().as_string();

  // ucvm
  pugi::xml_node l_ucvm = l_velMod.child("ucvm");
  if( l_ucvm.child("source_transformation") ) {
    m_ucvm.trafoSrc[0][0] = l_ucvm.child("source_transformation").child("x_0").text().as_double();
    m_ucvm.trafoSrc[0][1] = l_ucvm.child("source_transformation").child("x_1").text().as_double();
    m_ucvm.trafoSrc[0][2] = l_ucvm.child("source_transformation").child("x_2").text().as_double();

    m_ucvm.trafoSrc[1][0] = l_ucvm.child("source_transformation").child("y_0").text().as_double();
    m_ucvm.trafoSrc[1][1] = l_ucvm.child("source_transformation").child("y_1").text().as_double();
    m_ucvm.trafoSrc[1][2] = l_ucvm.child("source_transformation").child("y_2").text().as_double();

    m_ucvm.trafoSrc[2][0] = l_ucvm.child("source_transformation").child("z_0").text().as_double();
    m_ucvm.trafoSrc[2][1] = l_ucvm.child("source_transformation").child("z_1").text().as_double();
    m_ucvm.trafoSrc[2][2] = l_ucvm.child("source_transformation").child("z_2").text().as_double();
  }
  m_ucvm.projSrc = l_ucvm.child("projections").child("source").text().as_string();
  m_ucvm.projDes = l_ucvm.child("projections").child("destination").text().as_string();

  m_ucvm.models = l_ucvm.child("models").text().as_string();
  m_ucvm.modelType = l_ucvm.child("model_type").text().as_string();
  m_ucvm.crdMode = l_ucvm.child("coordinate_mode").text().as_string();
  m_ucvm.rule = l_ucvm.child("normalization_rule").text().as_string();
  m_ucvm.lowerToSurf = l_ucvm.child("lower_to_surface").text().as_bool();

  // tsunami models
  pugi::xml_node l_tsunami = l_velMod.child("tsunami");
  m_tsunami.dryTol = l_tsunami.child("dry_tolerance").text().as_double();
  m_tsunami.bath = l_tsunami.child("bathymetry").text().as_string();
  for( pugi::xml_node l_disp = l_tsunami.child("displacements"); l_disp; l_disp = l_disp.next_sibling("displacements") ) {
    m_tsunami.disp.push_back( l_disp.text().as_string() );
  }
  m_tsunami.expr = l_tsunami.child("expression").text().as_string();

  // mesh refinement
  pugi::xml_node l_ref = l_doc.child("edge_v").child("refinement");
  m_ref.expr = l_ref.child("expression").text().as_string();
  m_ref.out = l_ref.child("out").text().as_string();

  // time info
  pugi::xml_node l_time = l_doc.child("edge_v").child("time");
  m_nTsGroups = std::max( l_time.child("n_groups").text().as_uint(), 1u );
  m_funDt     = l_time.child("fundamental_time_step").text().as_double();
  m_tsOut     = l_time.child("files").child("out").child("time_steps").text().as_string();
}