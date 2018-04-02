/**
 * @file This file is part of EDGE.
 *
 * @author Junyi Qiu (juq005 AT ucsd.edu)
 * @author Rajdeep Konwar (rkonwar AT ucsd.edu)
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
 * This file contains the utility routines of Edge-V.
 **/

#include <fstream>
#include <omp.h>

#include "vm_utility.h"

int antnInit( antn_cfg & a_cfg, const std::string &cfg_f ) {
  a_cfg.antn_cfg_fn = cfg_f;

  a_cfg.min_vp          = 0;
  a_cfg.min_vs          = 0;
  a_cfg.min_vs2         = 0;
  a_cfg.max_vp_vs_ratio = 0;

  a_cfg.hypoc.lon = 0;
  a_cfg.hypoc.lat = 0;

  std::cout << "Reading Annotation Config File: " << cfg_f << " ... ";
  std::cout.flush();

  std::ifstream iMshFs( cfg_f.c_str(), std::ios::in );
  if( !iMshFs.is_open() ) {
    std::cout << "Failed." << std::endl;
    std::cerr << "Error: cannot open the config file." << std::endl;
    exit( -1 );
  }

  std::string lineBuf;
  while( getline( iMshFs, lineBuf ) ) {
    int i = -1, j;

    while( (++i < lineBuf.length()) && (lineBuf[i] == ' ') );

    if( (i >= lineBuf.length()) || (lineBuf[i] == '#') )
      continue;

    j = i - 1;

    while( (++j < lineBuf.length()) && (lineBuf[j] != '=') );

    if( j >= lineBuf.length() )
      continue;

    std::string varName   = lineBuf.substr( i, j - i );
    std::string varValue  = lineBuf.substr( j + 1 );

    if( varName.compare( "ucvm_config" ) == 0 )
      a_cfg.ucvm_cfg_fn     = varValue;
    else if( varName.compare( "ucvm_model_list" ) == 0 )
      a_cfg.ucvm_model_list = varValue;
    else if( varName.compare( "fault_input_file" ) == 0 )
      a_cfg.fault_input_fns.push_back( varValue );
    else if( varName.compare( "mesh_file" ) == 0 )
      a_cfg.mesh_fn         = varValue;
    else if( varName.compare( "node_vm_file" ) == 0 )
      a_cfg.vm_node_fn      = varValue;
    else if( varName.compare( "elmt_vm_file" ) == 0 )
      a_cfg.vm_elmt_fn      = varValue;
    else if( varName.compare( "h5m_file" ) == 0 )
      a_cfg.h5m_fn          = varValue;
    else
      std::cout << "\nUnknown setting (" << varName << "). Ignored." << std::endl;
  }

  a_cfg.ucvm_cmode  = UCVMCMODE;
  a_cfg.ucvm_type   = UCVMTYPE;

  a_cfg.min_vp = (a_cfg.min_vp == 0) ? MIN_VP : a_cfg.min_vp;
  a_cfg.min_vs = (a_cfg.min_vs == 0) ? MIN_VS : a_cfg.min_vs;
  a_cfg.min_vs2 = (a_cfg.min_vs2 == 0) ? MIN_VS2 : a_cfg.min_vs2;
  a_cfg.max_vp_vs_ratio = (a_cfg.max_vp_vs_ratio == 0) ? \
                          MAX_VP_VS_RATIO : a_cfg.max_vp_vs_ratio;

  a_cfg.hypoc.lon = CENTERICLON;
  a_cfg.hypoc.lat = CENTERICLAT;

  a_cfg.elmt_type = ELMTTYPE;

  std::cout << "Done!" << std::endl;

  return 0;
}

int ucvmInit( const antn_cfg &a_cfg ) {
  ucvm_init( a_cfg.ucvm_cfg_fn.c_str() );
  ucvm_add_model_list( a_cfg.ucvm_model_list.c_str() );
  ucvm_setparam( UCVM_PARAM_QUERY_MODE, a_cfg.ucvm_cmode );

  return 0;
}

int meshInit( moab_mesh &m_mesh, antn_cfg &a_cfg ) {
  if( m_mesh.intf != nullptr ) {
    std::cout << "Failed." << std::endl;
    std::cerr << "Error: The mesh has been initialized." << std::endl;
    exit( -1 );
  }

  m_mesh.intf = new moab::Core;
  moab::Interface *iface = m_mesh.intf;

  std::cout << "Reading Mesh File: " << a_cfg.mesh_fn << " ... " << std::endl;

  //! Load the mesh from msh file
  moab::ErrorCode rval = iface->load_mesh( a_cfg.mesh_fn.c_str() );
  if( rval != moab::MB_SUCCESS ) {
    std::cout << "Failed." << std::endl;
    std::cerr << "Error: cannot open the mesh file." << std::endl;
    exit( -1 );
  }
  assert( rval == moab::MB_SUCCESS );

  //! Get Nodes && Tets Statistics
  int numNodes, numElmts;
  moab::Range verts, elems;

  rval = iface->get_number_entities_by_type( 0, moab::MBVERTEX, numNodes );
  assert( rval == moab::MB_SUCCESS );

  rval = iface->get_number_entities_by_type( 0, moab::MBTET, numElmts );
  assert( rval == moab::MB_SUCCESS );

  rval = iface->get_entities_by_type( 0, moab::MBVERTEX, verts );
  assert( rval == moab::MB_SUCCESS );
  assert( verts.size() == numNodes );

  rval = iface->get_entities_by_type( 0, moab::MBTET, elems );
  assert( rval == moab::MB_SUCCESS );
  assert( elems.size() == numElmts );

  m_mesh.num_nodes = numNodes;
  m_mesh.num_elmts = numElmts;
  std::cout << " | Number of vertices is " << m_mesh.num_nodes << std::endl;
  std::cout << " | Number of elements is " << m_mesh.num_elmts << std::endl;

  std::cout << "Done!" << std::endl;

  return 0;
}

int meshFinalize( moab_mesh &m_mesh ) {
  if( m_mesh.intf != nullptr ) {
    delete m_mesh.intf;
  } else {
    std::cout << "Failed." << std::endl;
    std::cerr << "Error: mesh object corrupted.";
    exit( -1 );
  }

  return 0;
}

int vmNodeInit( vmodel &vm_nodes, const moab_mesh &m_mesh ) {
  vm_nodes.vm_list = new vm_datum[m_mesh.num_nodes];

  return 0;
}

int vmNodeFinalize( vmodel &vm_nodes ) {
  if( vm_nodes.vm_list != nullptr )
    delete[] vm_nodes.vm_list;

  return 0;
}

int vmElmtInit( vmodel &vm_elmts, const moab_mesh &m_mesh ) {
  vm_elmts.vm_list = new vm_datum[m_mesh.num_elmts];

  return 0;
}

int vmElmtFinalize( vmodel &vm_elmts ) {
  if( vm_elmts.vm_list != nullptr )
    delete[] vm_elmts.vm_list;

  return 0;
}

int workerInit( worker_reg &wrk_rg, int_v totalNum ) {
  int_v tid         = omp_get_thread_num();
  int_v numThrds    = omp_get_num_threads();
  wrk_rg.worker_tid = tid;

  int_v workSize    = (totalNum + numThrds - 1) / numThrds;
  int_v nPrivate    = std::min( workSize * (tid + 1), totalNum ) - workSize * tid;
  wrk_rg.work_size  = workSize;
  wrk_rg.num_prvt   = nPrivate;

  std::string pjInitParams = "+proj=tmerc +units=m +axis=enu +no_defs \
                          +datum=WGS84 +k=0.9996 +lon_0=-117.916 +lat_0=33.933";
  wrk_rg.pj_utm = pj_init_plus( pjInitParams.c_str() );
  wrk_rg.pj_geo = pj_init_plus( "+proj=latlong +datum=WGS84" );

  return 0;
}

int writeVMNodes( vmodel &vm_nodes, const antn_cfg &a_cfg,
                  const moab_mesh &m_mesh ) {
  std::cout << "Write Velocity Model: " << a_cfg.vm_node_fn << " ... ";
  std::cout.flush();

  std::ofstream oVmNodeFs( a_cfg.vm_node_fn.c_str(), std::ios::out );
  if( !oVmNodeFs.is_open() ) {
    std::cout << "Failed" << std::endl;
    std::cerr << "Error: cannot generate the velocity model for nodes."
              << std::endl;
    exit( -1 );
  }

  //! Write down the headers
  oVmNodeFs << "$UcvmModel\n";
  oVmNodeFs << a_cfg.ucvm_model_list << std::endl;
  oVmNodeFs << "$EndUcvmModel\n";

  //! Write down the velocity model data
  int_v numNodes = m_mesh.num_nodes;

  oVmNodeFs << "$NodesVelocityModel\n";
  oVmNodeFs << numNodes << std::endl;
  for( int_v pid = 0; pid < numNodes; pid++ ) {
    real data0 = vm_nodes.vm_list[pid].data[0];
    real data1 = vm_nodes.vm_list[pid].data[1];
    real data2 = vm_nodes.vm_list[pid].data[2];

    oVmNodeFs << (pid + 1) << " " << data0 << " " << data1 << " " << data2
              << std::endl;
    oVmNodeFs.flush();
  }
  oVmNodeFs << "$NodesVelocityModel\n";

  oVmNodeFs.close();
  std::cout << "Done!" << std::endl;

  return 0;
}

int writeVMElmts( vmodel &vm_elmts, const antn_cfg &a_cfg,
                  const moab_mesh &m_mesh ) {
  std::cout << "Writing Velocity Model: " << a_cfg.vm_elmt_fn << " ... ";
  std::cout.flush();

  std::ofstream oVmElmtFs( a_cfg.vm_elmt_fn.c_str(), std::ios::out );
  if( !oVmElmtFs.is_open() ) {
    std::cout << "Failed" << std::endl;
    std::cerr << "Error: cannot generate the velocity model for elements."
              << std::endl;
    exit( -1 );
  }

  //! Write down the headers
  oVmElmtFs << "$UcvmModel\n";
  oVmElmtFs << a_cfg.ucvm_model_list << std::endl;
  oVmElmtFs << "$EndUcvmModel\n";

  //! Write down the velocity model data
  int_v numElmts = m_mesh.num_elmts;

  oVmElmtFs << "$ElementsVelocityModel\n";
  oVmElmtFs << numElmts << std::endl;
  for( int_v eid = 0; eid < numElmts; eid++ ) {
    real data0 = vm_elmts.vm_list[eid].data[0];
    real data1 = vm_elmts.vm_list[eid].data[1];
    real data2 = vm_elmts.vm_list[eid].data[2];

    oVmElmtFs << (eid + 1) << " " << data0 << " " << data1 << " " << data2
              << std::endl;
    oVmElmtFs.flush();
  }
  oVmElmtFs << "$ElementsVelocityModel\n";

  oVmElmtFs.close();
  std::cout << "Done!" << std::endl;

  return 0;
}

int writeVMTags( vmodel &vm_elmts, const antn_cfg &a_cfg,
                 const moab_mesh &m_mesh ) {
  moab::Interface *iface = m_mesh.intf;
  if( iface == nullptr ) {
    std::cout << "Failed." << std::endl;
    std::cerr << "Error: cannot operate the mesh." << std::endl;
    exit( -1 );
  }

  std::cout << "Writing Annotated Mesh File: " << a_cfg.h5m_fn << " ... ";
  std::cout.flush();

  moab::Range elems;
  moab::ErrorCode rval = iface->get_entities_by_type( 0, moab::MBTET, elems );
  assert( rval == moab::MB_SUCCESS );

  //! Create Tags
  moab::Tag tag_lambda, tag_mu, tag_rho;
  rval = iface->tag_get_handle( "LAMBDA", 4, moab::MB_TYPE_OPAQUE, tag_lambda,
                    moab::MB_TAG_CREAT|moab::MB_TAG_DENSE|moab::MB_TAG_BYTES );
  assert( rval == moab::MB_SUCCESS );

  rval = iface->tag_get_handle( "MU", 4, moab::MB_TYPE_OPAQUE, tag_mu,
                    moab::MB_TAG_CREAT|moab::MB_TAG_DENSE|moab::MB_TAG_BYTES );
  assert( rval == moab::MB_SUCCESS );

  rval = iface->tag_get_handle( "RHO", 4, moab::MB_TYPE_OPAQUE, tag_rho,
                    moab::MB_TAG_CREAT|moab::MB_TAG_DENSE|moab::MB_TAG_BYTES );
  assert( rval == moab::MB_SUCCESS );

  //! Set Tags
  int_v eid;
  int_v numElmts = m_mesh.num_elmts;
  real l_lambda, l_mu, l_rho;

  for( moab::Range::iterator rit = elems.begin(); rit != elems.end(); ++rit ) {
    moab::EntityID entId = iface->id_from_handle( *rit );

    eid = entId - 1;
    assert( eid < numElmts );

    l_lambda  = vm_elmts.vm_list[eid].data[0];
    l_mu      = vm_elmts.vm_list[eid].data[1];
    l_rho     = vm_elmts.vm_list[eid].data[2];


    iface->tag_set_data( tag_lambda, &(*rit), 1, &l_lambda );
    assert( rval == moab::MB_SUCCESS );

    rval = iface->tag_set_data( tag_mu, &(*rit), 1, &l_mu );
    assert( rval == moab::MB_SUCCESS );

    rval = iface->tag_set_data( tag_rho, &(*rit), 1, &l_rho );
    assert( rval == moab::MB_SUCCESS );
  }

  //! Write to H5M File
  rval = iface->write_file( a_cfg.h5m_fn.c_str(), "H5M" );
  assert( rval == moab::MB_SUCCESS );

  std::cout << "Done!" << std::endl;

  return 0;
}



//
// Class definitions and helper methods for fault annotation
//
// A planar fault is assumed, with initial data given at a series of
// grid points.  We assume that the y-axis is coplanar with the fault
// plane, but the x and z axes need not be.
// Based on the SCEC/USGSG TPV35 Benchmark. For more information see:
// http://scecdata.usc.edu/cvws/download/TPV35_Description_v05.pdf
// http://scecdata.usc.edu/cvws/tpv35docs.html


// ************************************
// * Begin Triangle class definitions *
// ************************************
Triangle::Triangle( double* i_coords ){
  if( i_coords == nullptr ){
    std::cerr << "Error: Attempted to dereference null ptr in Triangle ctor"
              << std::endl;
    exit( -1 );
  }
  m_v1.x = i_coords[0];
  m_v1.y = i_coords[1];
  m_v1.z = i_coords[2];
  m_v2.x = i_coords[3];
  m_v2.y = i_coords[4];
  m_v2.z = i_coords[5];
  m_v3.x = i_coords[6];
  m_v3.y = i_coords[7];
  m_v3.z = i_coords[8];
}


xyz_point_t Triangle::centroid(){
  xyz_point_t l_centroid;
  l_centroid.x = ( m_v1.x + m_v2.x + m_v3.x ) / 3;
  l_centroid.y = ( m_v1.y + m_v2.y + m_v3.y ) / 3;
  l_centroid.z = ( m_v1.z + m_v2.z + m_v3.z ) / 3;

  return l_centroid;
}


xyz_point_t Triangle::bary_to_phy( xyz_point_t i_b ){
  xyz_point_t l_phys;
  l_phys.x = i_b.x * m_v1.x + i_b.y * m_v2.x + i_b.z * m_v3.x;
  l_phys.y = i_b.x * m_v1.y + i_b.y * m_v2.y + i_b.z * m_v3.y;
  l_phys.z = i_b.x * m_v1.z + i_b.y * m_v2.z + i_b.z * m_v3.z;

  return l_phys;
}


// refine is the number of new subdivisions made to each side of triangle
// EX: refine=0 --> No refinement;    refine=2 --> 9 subtriangles
std::vector< Triangle > Triangle::subdivide( unsigned int i_refine ){
  Triangle                l_subtri;
  std::vector< Triangle > l_subtriangles;
  xyz_point_t             l_bary1, l_bary2, l_bary3, l_bary4;

  double                  h = 1.0 / (i_refine + 1);

//  Reference:
//        v3
//        * *                         b3 * * * * * b4
//   ^    *   *                       *  *          *
//   |    *     *                     *    *        *
//   |    * * * * *                   *      *      *
//   | i  * *     * *                 *        *    *
//   |    *   *   *   *               *          *  *
//        *     * *     *             b1 * * * * * b2
//        v1* * * * * * * v2
//              j
//          ------->
//
// We order the subtriangles increasing from v1 to v2
// and then upward toward v3. Barycentric coordinates
// are used to map from the reference element to physical
// space.
  for( int i = 0; i <= i_refine; i++ ){
    for( int j = 0; j <= i_refine - i; j++ ){
      l_bary1.x = 1 - j*h - i*h;
      l_bary1.y = j*h;
      l_bary1.z = i*h;

      l_bary2.x = 1 - (j+1)*h - i*h;
      l_bary2.y = (j+1)*h;
      l_bary2.z = i*h;

      l_bary3.x = 1 - j*h - (i+1)*h;
      l_bary3.y = j*h;
      l_bary3.z = (i+1)*h;

      l_bary4.x = 1 - (j+1)*h + (i+1)*h;
      l_bary4.y = (j+1)*h;
      l_bary4.z = (i+1)*h;

      l_subtri.m_v1 = bary_to_phy( l_bary1 );
      l_subtri.m_v2 = bary_to_phy( l_bary2 );
      l_subtri.m_v3 = bary_to_phy( l_bary3 );
      l_subtriangles.push_back( l_subtri );

      if( j + i < i_refine ){
        l_subtri.m_v1 = bary_to_phy( l_bary4 );
        l_subtri.m_v2 = bary_to_phy( l_bary3 );
        l_subtri.m_v3 = bary_to_phy( l_bary2 );
        l_subtriangles.push_back( l_subtri );
      }
    }
  }
  return l_subtriangles;
}
// ************************************
// ** End Triangle class definitions **
// ************************************


// ************************************
// ** Begin FModel class definitions **
//*************************************
FModel::FModel( std::vector< std::string > i_fInputFns ){
  // Open Fault Input File
  m_filenames = i_fInputFns;
  std::cout << "Reading Fault Input File header from: " << m_filenames[0]
            << " ... " << std::flush;

  std::ifstream l_faultIfs( m_filenames[0], std::ios::in );
  if( !l_faultIfs.is_open() ){
    std::cout << "Failed." << std::endl;
    std::cerr << "Error: Cannot open the fault input file." << std::endl;
    exit( -1 );
  }

  // Read the header
  std::string l_lineBuf;
  getline( l_faultIfs, l_lineBuf );
  std::istringstream l_lineStream( l_lineBuf );

  l_lineStream >> m_Nx;
  l_lineStream >> m_Ny;
  l_lineStream >> m_xMin;
  l_lineStream >> m_xMax;
  l_lineStream >> m_yMin;
  l_lineStream >> m_yMax;

  l_faultIfs.close();
  std::cout << "Done!" << std::endl;

  size_t l_modelSize = (m_Nx + 1) * (m_Ny + 1);
  m_faultData = new FDatum[ l_modelSize ];

  // Finally, compute the inverse of the width of the grid cells for future
  // convenience, and set fault angle
  m_faultAngle = 0 * M_PI/180.;
  m_xScaleInv = m_Nx / (m_xMax - m_xMin);
  m_yScaleInv = m_Ny / (m_yMax - m_yMin);

  std::cout << " | Number of cells (x-axis): " << m_Nx << std::endl;
  std::cout << " | Number of cells (y-axis): " << m_Ny << std::endl;
  std::cout << " | Min x coordinate:         " << m_xMin << std::endl;
  std::cout << " | Max x coordinate:         " << m_xMax << std::endl;
  std::cout << " | Min y coordinate:         " << m_yMin << std::endl;
  std::cout << " | Max y coordinate:         " << m_yMax << std::endl;
  std::cout << " | Fault Angle:              " << m_faultAngle << std::endl;

  std::cout << "Fault model parameter setup completed." << std::endl;
}


FModel::~FModel(){
  if( m_faultData != nullptr ){
    delete[] m_faultData;
  }
  else{
    std::cerr << "Warning: Attempted to delete non-existant Fault Model Data"
              << std::endl;
  }
}


void FModel::populate(){
  int                 l_nx, l_ny;
  double              l_x, l_y;
  double              l_sFric, l_sStress, l_nStress;
  std::string         l_lineBuf;
  std::ifstream       l_faultIfs;
  std::istringstream  l_lineStream;

  std::cout << "Populating fault model from: " << std::endl;
  for( const auto& l_filename : m_filenames ){
    std::cout << "    > " << l_filename << std::endl;

    l_faultIfs.open( l_filename, std::ios::in );
    if( !l_faultIfs.is_open() ){
      std::cout << "Failed." << std::endl;
      std::cerr << "Error: Cannot open " << l_filename << std::endl;
      exit( -1 );
    }

    getline( l_faultIfs, l_lineBuf );   // Ignore the first line (header)
    while( getline( l_faultIfs, l_lineBuf ) ){
      l_lineStream.clear();
      l_lineStream.str( l_lineBuf );

      l_lineStream >> l_nx;    // integer x coord, (indexed from 0)
      l_lineStream >> l_ny;    // integer y coord, (indexed from 0)
      l_lineStream >> l_x;     // physical x coord
      l_lineStream >> l_y;     // physical y coord

      l_lineStream >> l_sFric;
      l_lineStream >> l_sStress;

      m_faultData[ l_nx + l_ny * (m_Nx+1) ].m_sFric.push_back( l_sFric );
      // m_faultData[ l_nx + l_ny * (m_Nx+1) ].m_nStress.push_back( m_sStress );
      m_faultData[ l_nx + l_ny * (m_Nx+1) ].m_nStress.push_back( 60. );
      m_faultData[ l_nx + l_ny * (m_Nx+1) ].m_sStress.push_back( l_sStress );
    }

    l_faultIfs.close();
  }
  std::cout << "...Done!" << std::endl;
}


// Assumes fault plane is rotated about axis x=0,z=0 only
// (i.e. no variation in y coord)
xyz_point_t FModel::proj_to_plane( xyz_point_t i_pt ){
  xyz_point_t l_proj;

  double s = sin( m_faultAngle );
  double c = cos( m_faultAngle );
  l_proj.x = i_pt.x * c * c + i_pt.z * s * c;
  l_proj.y = i_pt.y;
  l_proj.z = i_pt.x * s * c + i_pt.z * s * s;

  return l_proj;
}


int FModel::get_nearest_nx( xyz_point_t i_pt ){
  int l_nxNearest = std::round( m_xScaleInv * ( i_pt.x - m_xMin ) );

  if( l_nxNearest < 0 )
    l_nxNearest = 0;
  else if( l_nxNearest > m_Nx )
    l_nxNearest = m_Nx;

  return l_nxNearest;
}

int FModel::get_nearest_ny( xyz_point_t i_pt ){
  int l_nyNearest = std::round( m_yScaleInv * ( i_pt.y - m_yMin ) );

  if( l_nyNearest < 0 )
    l_nyNearest = 0;
  else if( l_nyNearest > m_Ny )
    l_nyNearest = m_Ny;

  return l_nyNearest;
}


FDatum FModel::get_datum_pt( xyz_point_t i_pt ){
  int         l_nx, l_ny;
  xyz_point_t l_planePt;

  l_planePt = proj_to_plane( i_pt );
  l_nx = get_nearest_nx( l_planePt );
  l_ny = get_nearest_ny( l_planePt );

  return m_faultData[ l_nx + l_ny * (m_Nx+1) ];
}


FDatum FModel::get_datum_tri( Triangle i_tri ){
  xyz_point_t l_centroid = i_tri.centroid();
  FDatum      l_datum = get_datum_pt( l_centroid );

  return l_datum;
}
// ************************************
// *** End FModel class definitions ***
// ************************************


moab::Range get_fault_faces( moab_mesh &i_mesh ){
  moab::Interface *iface = i_mesh.intf;
  moab::ErrorCode l_rval;

  moab::Range               l_faultFaces;
  moab::Range               l_mSets, l_tempRange;
  std::vector< moab::Tag >  l_tagsMesh;
  moab::Tag                 l_tagMat;
  std::string               l_tagMatName;

  // Get Material-Type tag...
  l_rval = iface->tag_get_tags( l_tagsMesh );
  assert( l_rval == moab::MB_SUCCESS );
  l_tagMat = l_tagsMesh[0];

  // ...and double check that we have the right tag
  l_rval = iface->tag_get_name( l_tagMat, l_tagMatName );
  assert( l_rval == moab::MB_SUCCESS );
  assert( l_tagMatName == "MATERIAL_SET" );

  // Get all entity sets
  l_rval = iface->get_entities_by_type( 0, moab::MBENTITYSET, l_mSets );
  assert( l_rval == moab::MB_SUCCESS );

  // Read off the material type for each entity set
  std::vector< int > l_matData;
  l_matData.resize( l_mSets.size() );
  iface->tag_get_data( l_tagMat, l_mSets, &l_matData[0] );
  assert( l_rval == moab::MB_SUCCESS );

  // Check if entity set has material type 201 (rupture)
  for( size_t l_md = 0; l_md < l_matData.size(); l_md++ ){
    if( l_matData[l_md] == 201 ){
      l_tempRange.clear();
      l_rval = iface->get_entities_by_handle( l_mSets[l_md], l_tempRange );
      l_faultFaces.merge( l_tempRange );
    }
  }

  return l_faultFaces;
}


// Create MOAB tag and tag each fault entity with array of doubles.
// Each tag is an array representing the different initial
// rupture parameters at each subtriangle and each fused run
// See diagram:
//
//  x------------------------------------------------------------x
//  |                         Fault Face                         |
//  x-----------------------------x------------------------------x
//  |         subtriangle1        |          subtriangle2        |
//  x---------x---------x---------x------------------------------x
//  | n_stress| s_stress| s_fric  |
//  x---------x---------x---------x
//  |cfr |cfr |
//  x----x----x
int fault_antn( moab_mesh &i_mesh, const antn_cfg &i_antnCfg ){
  if( i_mesh.intf == nullptr ){
    std::cout << "Failed." << std::endl;
    std::cerr << "Error: Attempted to annotate an uninitialized mesh (fault data)."
              << std::endl;
    exit( -1 );
  }

  // Set parameters for constructing and querying fault model
  const int l_nCfr = i_antnCfg.fault_input_fns.size();// Number of fused runs
  const int l_refine = 0;                             // Subtri. refinement lvl
  const int l_nSubTri = pow( l_refine+1, 2 );         // Subtris. per tet face

  // Set up fault_model
  FModel fault_model( i_antnCfg.fault_input_fns );
  fault_model.populate();

  // Get moab Range of all faces on fault
  moab::ErrorCode         l_rval;
  moab::Interface*        iface = i_mesh.intf;
  moab::Range             l_faultFaces = get_fault_faces( i_mesh );
  const int               l_numFaultFaces = l_faultFaces.size();
  const int               l_fdSize = l_nCfr * 3;
  const int               l_tagSize = l_nSubTri * l_fdSize;

  // Declare some more data structures to be used in the mesh annotation
  Triangle                l_tri;
  std::vector< Triangle > l_subtriangles;
  moab::Range             l_faceVerts;
  double                  l_faceVertsCoords[9];
  moab::Tag               l_tagFData;
  FDatum                  l_fd;

  double*                 l_entData = new double[ l_tagSize ];

  // Get tag to hold fault parameters, creating it if it doesn't exist
  l_rval = iface->tag_get_handle( "FAULT_DATA", l_tagSize, moab::MB_TYPE_DOUBLE,
                                l_tagFData,
                                moab::MB_TAG_CREAT|moab::MB_TAG_DENSE );
  assert( l_rval == moab::MB_SUCCESS );


  // Loop over all fault faces, annotating mesh as we go
  int l_count = 0;
  std::cout << "Annotating mesh with fault stresses ... " << std::endl;
  std::cout << l_count << " of " << l_numFaultFaces << " fault entities annotated."
            << std::flush;

  for( const auto& l_faceHandle : l_faultFaces ){
    l_faceVerts.clear();
    l_rval = iface->get_adjacencies( &l_faceHandle, 1, 0, false, l_faceVerts );
    assert( l_rval == moab::MB_SUCCESS );

    // Ensure that vertices in triangle are ordered according to increasing id
    // Can compare entity handles since entity type is the same for all verts
    assert( l_faceVerts[0] < l_faceVerts[1] );
    assert( l_faceVerts[1] < l_faceVerts[2] );
    l_rval = iface->get_coords( l_faceVerts, l_faceVertsCoords );
    assert( l_rval == moab::MB_SUCCESS );

    // Build a triangle from the vertices of triangle entity, and subdivide
    l_tri = Triangle( l_faceVertsCoords );
    l_subtriangles = l_tri.subdivide( l_refine );
    assert( l_nSubTri == l_subtriangles.size() );

    // Loop over each subtriangle and get approximate fault data
    for( int l_st = 0; l_st < l_nSubTri; l_st++ ){
      l_fd = fault_model.get_datum_tri( l_subtriangles[ l_st ] );

      for( int l_cfr = 0; l_cfr < l_nCfr; l_cfr++ ){
        l_entData[ l_st * l_fdSize + l_cfr ]            = l_fd.m_nStress[l_cfr];
        l_entData[ l_st * l_fdSize + l_nCfr + l_cfr ]   = l_fd.m_sStress[l_cfr];
        l_entData[ l_st * l_fdSize + 2*l_nCfr + l_cfr ] = l_fd.m_sFric[l_cfr];
      }
    }

    // Tag the mesh with the array of stresses
    l_rval = iface->tag_set_data( l_tagFData, &l_faceHandle, 1, l_entData );
    assert( l_rval == moab::MB_SUCCESS);

    // Print out the progress of the mesh annotation
    l_count++;
    std::cout << "\r" << std::setw(70) << "\r" << std::flush;
    std::cout << l_count << " of " << l_numFaultFaces << " completed"
              << std::flush;
  }
  std::cout << std::endl;
  std::cout << "Done!" << std::endl;

  return 0;
}
