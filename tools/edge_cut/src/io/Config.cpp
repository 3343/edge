/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (breuer AT mytum.de)
 * @author David Lenz (dlenz AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2020, Alexander Breuer
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
#include "Config.h"
#include "submodules/pugixml/src/pugixml.hpp"

edge_cut::io::Config::Config( std::string i_xmlPath ) {
  // parse config
  pugi::xml_document l_doc;
  pugi::xml_parse_result l_parseResult = l_doc.load_file( i_xmlPath.c_str() );

  // print errors
  if( !l_parseResult ) {
    EDGE_CUT_LOG_ERROR << "XML [" << i_xmlPath << "] parsed with errors, attr value: ["
                       << l_doc.child("node").attribute("attr").value() << "]";
    EDGE_CUT_LOG_FATAL << "Description: " << l_parseResult.description()
                       << " (error at ["   << i_xmlPath << " " << l_parseResult.offset << "]";
  }

  // check for EDGEcut root
  pugi::xml_node l_root = l_doc.child("edge_cut");
  if( l_root == pugi::xml_node() ) {
    EDGE_CUT_LOG_FATAL << "XML configuration does not contain <edge_cut> tree";
  }

  // read values
  m_gridIn = l_root.child("extrude").child("in").text().as_string();

  pugi::xml_node l_files = l_root.child("mesh").child("files");
  m_meshOutLeft   = l_files.child("out").child("left").text().as_string();
  m_meshOutRight  = l_files.child("out").child("right").text().as_string();
  m_meshOutFront  = l_files.child("out").child("front").text().as_string();
  m_meshOutBack   = l_files.child("out").child("back").text().as_string();
  m_meshOutBottom = l_files.child("out").child("bottom").text().as_string();
  m_meshOutTop    = l_files.child("out").child("top").text().as_string();

  m_extrudeZ = l_root.child("extrude").child("target_z").text().as_double();
}


void edge_cut::io::Config::print() {
  EDGE_CUT_LOG_INFO << "printing runtime configuration:";
  EDGE_CUT_LOG_INFO << "  input grid: " << m_gridIn;
  EDGE_CUT_LOG_INFO << "  output meshes:";
  EDGE_CUT_LOG_INFO << "    left:   " << m_meshOutLeft;
  EDGE_CUT_LOG_INFO << "    right:  " << m_meshOutRight;
  EDGE_CUT_LOG_INFO << "    front:  " << m_meshOutFront;
  EDGE_CUT_LOG_INFO << "    back:   " << m_meshOutBack;
  EDGE_CUT_LOG_INFO << "    bottom: " << m_meshOutBottom;
  EDGE_CUT_LOG_INFO << "    top:    " << m_meshOutTop;
}
