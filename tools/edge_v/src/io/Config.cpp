/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2019, Alexander Breuer
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

#include "io/logging.h"
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
  m_meshOut  = l_mesh.child("files").child("out").text().as_string();
  m_meshOutPa[0] = l_mesh.child("files").child("out_by_partition").child("base").text().as_string();
  m_meshOutPa[1] = l_mesh.child("files").child("out_by_partition").child("extension").text().as_string();
  m_writeElAn = l_mesh.child("write_element_annotations").text().as_bool();
  m_periodic = l_mesh.child("periodic").text().as_bool();
  m_nPartitions = l_mesh.child("n_partitions").text().as_ullong();

  // read velocity model
  pugi::xml_node l_velMod = l_doc.child("edge_v").child("velocity_model");
  m_seismicExpr = l_velMod.child("seismic_expression").text().as_string();

  // time info
  pugi::xml_node l_time = l_doc.child("edge_v").child("time");
  m_nTsGroups = l_time.child("n_groups").text().as_uint();
  m_funDt     = l_time.child("fundamental_time_step").text().as_double();
  m_tsOut     = l_time.child("files").child("out_time_steps").text().as_string();
}