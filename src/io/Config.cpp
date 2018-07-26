/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2016-2018, Regents of the University of California
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
 * Runtime configuration.
 **/

#include "Config.h"
#include "logging.h"
#include <string>
#include <sstream>
#ifdef PP_USE_MEMKIND
#include <memkind.h>
#endif

void edge::io::Config::printMultiLine( std::string i_pre,
                                       std::string i_str ) {
  if( i_str.length() > 0 ) {
    std::stringstream l_ss( i_str.c_str() );
    std::string l_line;
    while( std::getline( l_ss, l_line, '\n' ) ) {
      EDGE_LOG_INFO << i_pre << l_line;
    }
  }
}

void edge::io::Config::printBuild( pugi::xml_node i_build ) {
  // convert build to string
#if   defined __clang__
  std::string l_compiler = "clang";
#elif defined __ICC
  std::string l_compiler = "intel";
              l_compiler += " ";
              l_compiler += std::to_string( (long long) __INTEL_COMPILER );
              l_compiler += ".";
              l_compiler += std::to_string( (long long) __INTEL_COMPILER_UPDATE );
#elif defined __GNUC__
  std::string l_compiler = "gcc";
#else
  std::string l_compiler = "unknown";
#endif

#if defined __AVX512F__
  std::string l_instSet = "avx512f";
#elif defined __AVX2__
  std::string l_instSet = "avx2";
#elif defined __KNC__
  std::string l_instSet = "knc";
#elif defined __AVX__
  std::string l_instSet = "avx";
#elif defined __SSE2__
  std::string l_instSet = "sse2";
#else
  std::string l_instSet = "unknown";
#endif

  std::string l_elementType;
  if(       T_SDISC.ELEMENT == LINE   ) l_elementType = "line";
  else if ( T_SDISC.ELEMENT == QUAD4R ) l_elementType = "quad4r";
  else if ( T_SDISC.ELEMENT == TRIA3  ) l_elementType = "tria3";
  else if(  T_SDISC.ELEMENT == HEX8R  ) l_elementType = "hex8r";
  else if(  T_SDISC.ELEMENT == TET4   ) l_elementType = "tet4";
  else {
   EDGE_LOG_FATAL << "not good, more element types are required! " << T_SDISC.ELEMENT;
  }

#if defined PP_T_MESH_REGULAR
  std::string l_meshType = "regular";
#elif defined PP_T_MESH_UNSTRUCTURED
  std::string l_meshType = "unstructured";
#else
#error "that's not right, strange mesh type.."
#endif

#if defined PP_T_EQUATIONS_ADVECTION
  std::string l_equations = "advection";
#elif defined PP_T_EQUATIONS_ELASTIC
  std::string l_equations = "elastic";
#elif defined PP_T_EQUATIONS_SWE
  std::string l_equations = "swe";
#else
#error "somebody should implement another set of equations here"
#endif

#if defined PP_T_KERNELS_XSMM_DENSE_SINGLE || defined PP_T_KERNELS_XSMM
  std::string l_xsmm = "yes";
#else
  std::string l_xsmm = "no";
#endif

  EDGE_LOG_INFO << "sharing the build config, note that every change requires a recompile";
  EDGE_LOG_INFO << "  running without a recompile will just use the settings of the";
  EDGE_LOG_INFO << "  last compile and ignore changed runtime parameters:";
  EDGE_LOG_INFO << "    date / time of the build:          " << __DATE__      << " / " << __TIME__;
  EDGE_LOG_INFO << "    compiler name / version:           " << l_compiler    << " / " << __VERSION__;
  EDGE_LOG_INFO << "    inst. set build / arch runtime:    " << l_instSet     << " / " << i_build.child("arch").text().as_string();
  EDGE_LOG_INFO << "    cfr (build / runtime):             " << PP_N_CRUNS    << " / " << i_build.child("cfr").text().as_string();
  EDGE_LOG_INFO << "    equations       (build / runtime): " << l_equations   << " / " << i_build.child("equations").text().as_string();
  EDGE_LOG_INFO << "    element_type    (build / runtime): " << l_elementType << " / " << i_build.child("element_type").text().as_string();
  EDGE_LOG_INFO << "    order           (build / runtime): " << PP_ORDER      << " / " << i_build.child("order").text().as_string();
  EDGE_LOG_INFO << "    precision       (build / runtime): " << PP_PRECISION  << " / " << i_build.child("precision").text().as_string();
  EDGE_LOG_INFO << "    mesh type build / moab runtime:    " << l_meshType    << " / " << i_build.child("moab").text().as_string();
  EDGE_LOG_INFO << "    xsmm            (build / runtime): " << l_xsmm        << " / " << i_build.child("xsmm").text().as_string();
  EDGE_LOG_INFO << "some derived parameters: ";
  EDGE_LOG_INFO << "  #element modes: " << N_ELEMENT_MODES;
#if defined PP_USE_MEMKIND
  EDGE_LOG_INFO << "  memkind version: " << memkind_get_version();
#endif
  EDGE_LOG_INFO << "  memory alignment (bytes): ";
  EDGE_LOG_INFO << "    stack base pointers:   " << ALIGNMENT.BASE.STACK;
  EDGE_LOG_INFO << "    heap base pointers:    " << ALIGNMENT.BASE.HEAP;
  EDGE_LOG_INFO << "    cruns:                 " << ALIGNMENT.CRUNS;
  EDGE_LOG_INFO << "    private element modes: " << ALIGNMENT.ELEMENT_MODES.PRIVATE;
  EDGE_LOG_INFO << "    shared element modes:  " << ALIGNMENT.ELEMENT_MODES.SHARED;
}

void edge::io::Config::printConfig() {
  EDGE_LOG_INFO << "printing the runtime config, just for you!";

  EDGE_LOG_INFO << "  here's the mesh:";
#ifdef PP_T_MESH_REGULAR
  EDGE_LOG_INFO << "    n_elements: ";
  EDGE_LOG_INFO << "      x: " << m_nElementsX;
  EDGE_LOG_INFO << "      y: " << m_nElementsY;
  EDGE_LOG_INFO << "      z: " << m_nElementsZ;

  EDGE_LOG_INFO << "    size: ";
  EDGE_LOG_INFO << "      x: " << m_sizeX;
  EDGE_LOG_INFO << "      y: " << m_sizeY;
  EDGE_LOG_INFO << "      z: " << m_sizeZ;
#elif defined PP_T_MESH_UNSTRUCTURED
  EDGE_LOG_INFO << "    read options: " << m_meshOptRead;
  EDGE_LOG_INFO << "    files:";
  EDGE_LOG_INFO << "      in: " << m_meshFileIn;
  EDGE_LOG_INFO << "      out: " << m_meshFileOut;
  EDGE_LOG_INFO << "    boundary:";
  if( m_periodic != std::numeric_limits<int>::max() ) {
    EDGE_LOG_INFO << "      periodic: " << m_periodic;
  }
  for( unsigned int l_bn = 0; l_bn < m_bndConId.size(); l_bn++ ) {
    EDGE_LOG_INFO << "      " << m_bndConName[l_bn] << ": " << m_bndConId[l_bn];
  }
#endif

  // print info about sparse type domains
  for( unsigned short l_et = 0; l_et < 3; l_et++ ) {
    if( m_spTypesDoms[l_et].size() > 0 ) {
      EDGE_LOG_INFO << "  found " << m_spTypesDoms[l_et].size() << " domains with sparse types for the "
                    << (l_et == 0 ? "vertices" : (l_et == 1 ? "faces" : "elements") )
                    << ":";
      for( std::size_t l_do = 0; l_do < m_spTypesDoms[l_et].size(); l_do++ ) {
        EDGE_LOG_INFO << "    domain #" << l_do << ": " << m_spTypesVals[l_et][l_do];
        std::vector< std::string > l_doStrs = m_spTypesDoms[l_et][l_do].toString();
        for( std::size_t l_ob = 0; l_ob < l_doStrs.size(); l_ob++ ) {
          EDGE_LOG_INFO << "      object #" << l_ob << ": " << l_doStrs[l_ob];
        }
      }
    }
  }

  EDGE_LOG_INFO << "  more? get ready for the setups:";
  for( unsigned short l_cr = 0; l_cr < N_CRUNS; l_cr++ ) {
    std::string l_pre = "    initial_values #" + std::to_string(l_cr) + ": ";
    printMultiLine( l_pre,
                    m_initValsExprStrs[l_cr].c_str() );
  }
  EDGE_LOG_INFO << "  alright, we are sharing parameters again:";
  EDGE_LOG_INFO << "    end_time: " << m_endTime;

  if( m_waveFieldType != "" ) {
    EDGE_LOG_INFO << "  wave_field:";
    EDGE_LOG_INFO << "    type: " << m_waveFieldType;
    if( m_waveFieldSpType != std::numeric_limits< int_spType >::max() )
      EDGE_LOG_INFO << "    sparse_type: " << m_waveFieldSpType;
    EDGE_LOG_INFO << "    file: " << m_waveFieldFile;
    EDGE_LOG_INFO << "    int: "  << m_waveFieldInt;
  }

  if( m_iBndType != "" ) {
    EDGE_LOG_INFO << "  internal_boundary:";
    EDGE_LOG_INFO << "    type: " << m_iBndType;
    EDGE_LOG_INFO << "    file: " << m_iBndFile;
    EDGE_LOG_INFO << "    int: "  << m_iBndInt;
  }

  // iterate over receiver types. 0: element-modal, 1: face-quad
  for( unsigned short l_rt = 0; l_rt < 2; l_rt++ ) {
    std::string l_type;
    if( l_rt == 0 ) l_type = "element modal";
    else            l_type = "face quad";

    if( m_recvNames[l_rt].size() > 0 ) {
      EDGE_LOG_INFO << "  we have " << m_recvNames[l_rt].size() << " " << l_type << " receivers in the config: ";
      EDGE_LOG_INFO << "    sampling frequency: " << m_recvFreq[l_rt];
      EDGE_LOG_INFO << "    path to out-directory: "<< m_recvPath[l_rt];
      EDGE_LOG_INFO << "    and here they are, our receivers: ";
      for( std::size_t l_re = 0; l_re < m_recvNames[l_rt].size(); l_re++ ) {
        EDGE_LOG_INFO << "      #" << l_re << ":";
        EDGE_LOG_INFO << "        name: " << m_recvNames[l_rt][l_re];
#if PP_N_DIM == 1
        EDGE_LOG_INFO << "        x: " << m_recvCrds[l_rt][l_re][0];
#elif PP_N_DIM == 2
        EDGE_LOG_INFO << "        x, y: " << m_recvCrds[l_rt][l_re][0] << ", " << m_recvCrds[l_rt][l_re][1];
#else
        EDGE_LOG_INFO << "        x, y, z: " << m_recvCrds[l_rt][l_re][0] << ", "
                                             << m_recvCrds[l_rt][l_re][1] << ", "
                                             << m_recvCrds[l_rt][l_re][2];
#endif
      }
    }
  }

  if( m_errorNormsType != "" ) {
    EDGE_LOG_INFO << "  error_norms:";
    for( unsigned short l_cr = 0; l_cr < N_CRUNS; l_cr++ ) {
      std::string l_pre = "    reference_values #" + std::to_string(l_cr) + ": ";
      printMultiLine( l_pre,
                      m_refValsExprStrs[l_cr].c_str() );
    }
    EDGE_LOG_INFO << "    type: " << m_errorNormsType;
    EDGE_LOG_INFO << "    file: " << m_errorNormsFile;
  }
}

edge::io::Config::Config( std::string i_xmlPath ):
 m_periodic(std::numeric_limits<int>::max()) {
  pugi::xml_parse_result l_parseResult = m_doc.load_file( i_xmlPath.c_str() );

  // inform user about errors in xml-file
  if(! l_parseResult ) {
    EDGE_LOG_ERROR << "XML [" << i_xmlPath << "] parsed with errors, attr value: ["
                   << m_doc.child("node").attribute("attr").value() << "]";
    EDGE_LOG_FATAL << "Description: " << l_parseResult.description()
                   << " (error at ["   << i_xmlPath << " " << l_parseResult.offset << "]";
  }

  /*
   * print the build
   */
  pugi::xml_node l_build = m_doc.child("edge").child("build");
  printBuild( l_build );

  /*
   * read mesh parameters
   */
  pugi::xml_node l_mesh = m_doc.child("edge").child("cfr").child("mesh");

#ifdef PUGIXML_HAS_LONG_LONG
  m_nElementsX = l_mesh.child("n_elements").child("x").text().as_ullong();
  m_nElementsY = l_mesh.child("n_elements").child("y").text().as_ullong();
  m_nElementsZ = l_mesh.child("n_elements").child("z").text().as_ullong();
#else
  m_nElementsX = l_mesh.child("n_elements").child("x").text().as_uint();
  m_nElementsY = l_mesh.child("n_elements").child("y").text().as_uint();
  m_nElementsZ = l_mesh.child("n_elements").child("z").text().as_uint();
#endif

  m_sizeX = l_mesh.child("size").child("x").text().as_double();
  m_sizeY = l_mesh.child("size").child("y").text().as_double();
  m_sizeZ = l_mesh.child("size").child("z").text().as_double();

  m_meshOptRead = l_mesh.child("options").child("read").text().as_string();
  // set default read options for MPI if not specified
  if( m_meshOptRead == "" ) {
#ifdef PP_USE_MPI
    m_meshOptRead  = "PARALLEL=BCAST_DELETE;PARALLEL_RESOLVE_SHARED_ENTS;PARTITION=PARALLEL_PARTITION;";
#endif
  }

  m_meshFileIn = l_mesh.child("files").child("in").text().as_string();
  m_meshFileOut = l_mesh.child("files").child("out").text().as_string();

  // set periodic boundary value if present
  if( l_mesh.child("boundary").find_child(
       []( pugi::xml_node i_node ){ return std::string(i_node.name()) == "periodic";} ) ) {
    m_periodic = l_mesh.child("boundary").child("periodic").text().as_int();
  }

  // read boundary conditions (non-periodic)
  for (pugi::xml_node l_bd = l_mesh.child("boundary").first_child(); l_bd; l_bd = l_bd.next_sibling() ) {
    std::string l_bdName = l_bd.name();
    if( l_bdName != "periodic" ) {
      m_bndConId.push_back(   l_bd.text().as_int() );
      m_bndConName.push_back( l_bd.name()          );
    }
  }

  // read sparse types
  pugi::xml_node l_spTypes = l_mesh.child( "sparse_types" );

  // iterate over domain of vertices, faces and elements
  for( unsigned short l_et = 0; l_et < 3; l_et++ ) {
    std::string l_enName = ( l_et == 0 ) ? "vertices" :
                           ( l_et == 1 ) ? "faces" :
                                           "elements";

    // get the node for this entity type
    pugi::xml_node l_etSp = l_spTypes.child( l_enName.c_str() );

    // iterate over domains
    for( pugi::xml_node l_do = l_etSp.child("domain"); l_do; l_do = l_do.next_sibling("domain") ) {
      // resize arrays
      m_spTypesDoms[l_et].resize( m_spTypesDoms[l_et].size()+1 );

      // names of the data
      std::vector< std::vector< std::string > > l_dataN(1);
      l_dataN[0].push_back( "type" );

      // parse the domain-node
      ConfigDoms< N_DIM >::parse( l_do,
                                  l_dataN,
                                  ConfigDoms< N_DIM >::INT,
                                  m_spTypesVals[l_et],
                                  m_spTypesDoms[l_et].back() );

      // check for exactly one type per domain
      EDGE_CHECK_EQ( m_spTypesVals[l_et].size(), m_spTypesDoms[l_et].size() );
    }
  }

  /*
   * read setups
   */
  pugi::xml_node l_setups = m_doc.child("edge").child("cfr").child("setups");

  // read initial DOFs
  pugi::xml_node l_initVals = l_setups.child("initial_values");
  unsigned short l_crx = 0;
  for( pugi::xml_node l_id = l_initVals.first_child(); l_id; l_id = l_id.next_sibling() ) {
    if( l_crx < N_CRUNS ) m_initValsExprStrs[l_crx] = l_id.text().as_string();
    l_crx++;
  }

  // read parameters shared among setups
  m_endTime = l_setups.child("end_time").text().as_double();

  /*
   * read output
   */
  pugi::xml_node l_output = m_doc.child("edge").child("cfr").child("output");
  m_waveFieldType = l_output.child("wave_field").child("type").text().as_string();
  if( m_waveFieldType != "" ) {
    m_waveFieldFile = l_output.child("wave_field").child("file").text().as_string();
    m_waveFieldInt  = l_output.child("wave_field").child("int").text().as_double();

    if( l_output.child("wave_field").find_child([]( pugi::xml_node i_node ){ return std::string(i_node.name()) == "sparse_type";}) )
      m_waveFieldSpType = l_output.child("wave_field").child("sparse_type").text().as_uint();
  }
  EDGE_CHECK_GT( m_waveFieldInt, TOL.TIME );

  m_iBndType = l_output.child("internal_boundary").child("type").text().as_string();
  if( m_iBndType != "" ) {
    m_iBndFile = l_output.child("internal_boundary").child("file").text().as_string();
    m_iBndInt  = l_output.child("internal_boundary").child("int").text().as_double();
  }

  // iterate over receiver types. 0: element-modal, 1: face-quad
  for( unsigned short l_rt = 0; l_rt < 2; l_rt++ ) {
    std::string l_type;
    if( l_rt == 0 ) l_type = "receivers";
    else            l_type = "receivers_quad";

    for( pugi::xml_node l_recv = l_output.child(l_type.c_str()).child("receiver"); l_recv; l_recv = l_recv.next_sibling("receiver") ) {
      m_recvCrds[l_rt].resize( m_recvCrds[l_rt].size()+1 );

      m_recvCrds[l_rt].back()[0] = l_recv.child("coords").child("x").text().as_double();
      m_recvCrds[l_rt].back()[1] = l_recv.child("coords").child("y").text().as_double();
      m_recvCrds[l_rt].back()[2] = l_recv.child("coords").child("z").text().as_double();

      std::string l_recvName = l_recv.child("name").text().as_string();
      if( l_recvName != "" ) {
        m_recvNames[l_rt].push_back( l_recvName );
      }
      else m_recvNames[l_rt].push_back( std::to_string(m_recvNames[l_rt].size()) );
    }
    if( m_recvCrds[l_rt].size() > 0 ) {
      double l_freq = l_output.child(l_type.c_str()).child("freq").text().as_double();

      if( l_freq > TOL.TIME ) {
        m_recvFreq[l_rt] = l_freq;
      }
      else m_recvFreq[l_rt] = -std::numeric_limits< double >::max();
    }
    else m_recvFreq[l_rt] = -std::numeric_limits< double >::max();
    m_recvPath[l_rt] = l_output.child(l_type.c_str()).child("path_to_dir").text().as_string();
    // clear invalid input
    if( m_recvFreq[l_rt] < TOL.TIME || m_recvPath[l_rt] == "" ) {
      m_recvCrds[l_rt].clear();
      m_recvNames[l_rt].clear();
    }
  }

  /*
   * error computations
   */
  pugi::xml_node l_refVals = l_output.child("error_norms").child("reference_values");
  l_crx = 0;
  for( pugi::xml_node l_id = l_refVals.first_child(); l_id; l_id = l_id.next_sibling() ) {
    if( l_crx < N_CRUNS ) m_refValsExprStrs[l_crx] = l_id.text().as_string();
    l_crx++;
  }

  m_errorNormsType = l_output.child("error_norms").child("type").text().as_string();
  m_errorNormsFile = l_output.child("error_norms").child("file").text().as_string();

  // print config
  printConfig();
}
