/**
 * @file This file is part of EDGE.
 *
 * @author David Lenz (dlenz AT ucsd.edu)
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2017-2018, Regents of the University of California
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
 *
 *  Class definitions and helper methods for fault annotation
 *
 *  A planar fault is assumed, with initial data given at a series of
 *  grid points.  We assume that the y-axis is coplanar with the fault
 *  plane, but the x and z axes need not be.
 *  Based on the SCEC/USGSG TPV35 Benchmark. For more information see:
 *  http://scecdata.usc.edu/cvws/download/TPV35_Description_v05.pdf
 *  http://scecdata.usc.edu/cvws/tpv35docs.html
 **/
#include <iostream>
#include <fstream>
#include <cmath>
#include <cassert>
#include "FaultModel.h"


Triangle::Triangle( const real *i_coords ) {
  if( i_coords == nullptr ) {
    std::cerr << "Error: Attempted to dereference null ptr in Triangle ctor"
              << std::endl;
    exit( EXIT_FAILURE );
  }

  m_v1.m_x = i_coords[0];
  m_v1.m_y = i_coords[1];
  m_v1.m_z = i_coords[2];
  m_v2.m_x = i_coords[3];
  m_v2.m_y = i_coords[4];
  m_v2.m_z = i_coords[5];
  m_v3.m_x = i_coords[6];
  m_v3.m_y = i_coords[7];
  m_v3.m_z = i_coords[8];
}

xyz_point_t Triangle::centroid() {
  xyz_point_t l_centroid;
  l_centroid.m_x = ( m_v1.m_x + m_v2.m_x + m_v3.m_x ) / 3;
  l_centroid.m_y = ( m_v1.m_y + m_v2.m_y + m_v3.m_y ) / 3;
  l_centroid.m_z = ( m_v1.m_z + m_v2.m_z + m_v3.m_z ) / 3;

  return l_centroid;
}

xyz_point_t Triangle::baryToPhy( const xyz_point_t &i_b ) {
  xyz_point_t l_phys;
  l_phys.m_x = i_b.m_x * m_v1.m_x + i_b.m_y * m_v2.m_x + i_b.m_z * m_v3.m_x;
  l_phys.m_y = i_b.m_x * m_v1.m_y + i_b.m_y * m_v2.m_y + i_b.m_z * m_v3.m_y;
  l_phys.m_z = i_b.m_x * m_v1.m_z + i_b.m_y * m_v2.m_z + i_b.m_z * m_v3.m_z;

  return l_phys;
}

/** refine is the number of new subdivisions made to each side of triangle
 *  EX: refine=0 --> No refinement;    refine=2 --> 9 subtriangles
 **/
std::vector< Triangle > Triangle::subdivide( const int &i_refine ) {
  Triangle                l_subtri;
  std::vector< Triangle > l_subtriangles;
  xyz_point_t             l_bary1, l_bary2, l_bary3, l_bary4;

  double                  l_h = 1.0 / (double) i_refine;

/**  Reference:
 *        v3
 *        * *                         b3 * * * * * b4
 *   ^    *   *                       *  *          *
 *   |    *     *                     *    *        *
 *   |    * * * * *                   *      *      *
 *   | i  * *     * *                 *        *    *
 *   |    *   *   *   *               *          *  *
 *        *     * *     *             b1 * * * * * b2
 *        v1* * * * * * * v2
 *              j
 *          ------->
 *
 *  We order the subtriangles increasing from v1 to v2
 *  and then upward toward v3. Barycentric coordinates
 *  are used to map from the reference element to physical
 *  space.
 **/
  for( int l_i = 0; l_i < i_refine; l_i++ ) {
    for( int l_j = 0; l_j < i_refine - l_i; l_j++ ) {
      l_bary1.m_x = 1 - l_j * l_h - l_i * l_h;
      l_bary1.m_y = l_j * l_h;
      l_bary1.m_z = l_i * l_h;

      l_bary2.m_x = 1 - (l_j + 1) * l_h - l_i * l_h;
      l_bary2.m_y = (l_j + 1) * l_h;
      l_bary2.m_z = l_i * l_h;

      l_bary3.m_x = 1 - l_j * l_h - (l_i + 1) * l_h;
      l_bary3.m_y = l_j * l_h;
      l_bary3.m_z = (l_i + 1) * l_h;

      l_bary4.m_x = 1 - (l_j + 1) * l_h - (l_i + 1) * l_h;
      l_bary4.m_y = (l_j + 1) * l_h;
      l_bary4.m_z = (l_i + 1) * l_h;

      l_subtri.m_v1 = baryToPhy( l_bary1 );
      l_subtri.m_v2 = baryToPhy( l_bary2 );
      l_subtri.m_v3 = baryToPhy( l_bary3 );
      l_subtriangles.push_back( l_subtri );

      if( l_j + l_i < i_refine-1 ) {
        l_subtri.m_v1 = baryToPhy( l_bary4 );
        l_subtri.m_v2 = baryToPhy( l_bary3 );
        l_subtri.m_v3 = baryToPhy( l_bary2 );
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
edge_v::fault::Model::Model( const std::vector< std::string >& i_fInputFns ){
  // Open Fault Input File
  m_filenames = i_fInputFns;

  std::cout << "Initializing fault model from: " << std::endl;
  for( const auto& l_filename : m_filenames ){
    std::cout << "    > " << l_filename << std::endl;
  }

  std::ifstream l_faultIfs( m_filenames[0], std::ios::in );
  if( !l_faultIfs.is_open() ){
    std::cout << "Failed." << std::endl;
    std::cerr << "Error: Cannot open " << m_filenames[0] << std::endl;
    exit( EXIT_FAILURE );
  }

  // Read the header
  std::string l_lineBuf;
  std::string l_wordBuf;
  getline( l_faultIfs, l_lineBuf );
  std::istringstream l_lineStream( l_lineBuf );

  l_lineStream >> m_Nx;
  l_lineStream >> m_Ny;
  l_lineStream >> m_xMin;
  l_lineStream >> m_xMax;
  l_lineStream >> m_yMin;
  l_lineStream >> m_yMax;

  while( l_lineStream >> l_wordBuf){
    m_qtyNames.push_back( l_wordBuf );
  }

  size_t l_modelSize = (m_Nx + 1) * (m_Ny + 1);
  m_faultData = new FDatum[ l_modelSize ];

  // Finally, compute the inverse of the width of the grid cells for future
  // convenience, and set fault angle
  real FAULT_ANGLE = 0;
  m_faultAngle = FAULT_ANGLE * M_PI/180.;
  m_xScaleInv = m_Nx / (m_xMax - m_xMin);
  m_yScaleInv = m_Ny / (m_yMax - m_yMin);

  l_faultIfs.close();

  populate();

  std::cout << "...Done!" << std::flush;
}


edge_v::fault::Model::~Model(){
  if( m_faultData != nullptr ){
    delete[] m_faultData;
  }
  else{
    std::cerr << "Warning: Attempted to delete non-existant Fault Model Data"
              << std::endl;
  }
}


void edge_v::fault::Model::populate(){
  int                 l_nx, l_ny;
  real                l_x, l_y, l_z;
  real                l_numBuf;
  std::string         l_lineBuf;
  std::ifstream       l_faultIfs;
  std::istringstream  l_lineStream;

  for( const auto& l_filename : m_filenames ){
    l_faultIfs.open( l_filename, std::ios::in );
    if( !l_faultIfs.is_open() ){
      std::cout << "Failed." << std::endl;
      std::cerr << "Error: Cannot open " << l_filename << std::endl;
      exit( EXIT_FAILURE );
    }

    getline( l_faultIfs, l_lineBuf );   // Ignore the first line (header)
    while( getline( l_faultIfs, l_lineBuf ) ){
      l_lineStream.clear();
      l_lineStream.str( l_lineBuf );

      l_lineStream >> l_nx;    // integer x coord, (indexed from 0)
      l_lineStream >> l_ny;    // integer y coord, (indexed from 0)
      l_lineStream >> l_x;     // physical x coord
      l_lineStream >> l_y;     // physical y coord
      l_lineStream >> l_z;     // physical z coord

      for( const auto& l_qtyName : m_qtyNames ){
        l_lineStream >> l_numBuf;
        m_faultData[ l_nx + l_ny * (m_Nx+1) ][ l_qtyName ].push_back( l_numBuf );
      }
    }

    l_faultIfs.close();
  }
}

//! Assumes fault plane is rotated about axis x=0,z=0 only
//! (i.e. no variation in y coord)
xyz_point_t edge_v::fault::Model::projToPlane( const xyz_point_t &i_pt ) {
  xyz_point_t l_proj;

  real l_s = sin( m_faultAngle );
  real l_c = cos( m_faultAngle );
  l_proj.m_x = i_pt.m_x * l_c * l_c + i_pt.m_z * l_s * l_c;
  l_proj.m_y = i_pt.m_y;
  l_proj.m_z = i_pt.m_x * l_s * l_c + i_pt.m_z * l_s * l_s;

  return l_proj;
}

int edge_v::fault::Model::getNearestNx( const xyz_point_t &i_pt ) {
  int l_nxNearest = std::round( m_xScaleInv * ( i_pt.m_x - m_xMin ) );

  if( l_nxNearest < 0 )
    l_nxNearest = 0;
  else if( l_nxNearest > m_Nx )
    l_nxNearest = m_Nx;

  return l_nxNearest;
}

int edge_v::fault::Model::getNearestNy( const xyz_point_t &i_pt ) {
  int l_nyNearest = std::round( m_yScaleInv * ( i_pt.m_y - m_yMin ) );

  if( l_nyNearest < 0 )
    l_nyNearest = 0;
  else if( l_nyNearest > m_Ny )
    l_nyNearest = m_Ny;

  return l_nyNearest;
}

FDatum edge_v::fault::Model::getDatumPt( const xyz_point_t &i_pt ) {
  int         l_nx, l_ny;
  xyz_point_t l_planePt;

  l_planePt = projToPlane( i_pt );
  l_nx      = getNearestNx( l_planePt );
  l_ny      = getNearestNy( l_planePt );

  return m_faultData[l_nx+l_ny*(m_Nx+1)];
}

FDatum edge_v::fault::Model::getDatumTri( Triangle &i_tri ) {
  xyz_point_t l_centroid = i_tri.centroid();
  FDatum      l_datum    = getDatumPt( l_centroid );

  return l_datum;
}
// ************************************
// *** End FModel class definitions ***
// ************************************


moab::Range getFaultFaces( moab_mesh &i_mesh ) {
  moab::Interface *l_face = i_mesh.m_intf;
  moab::ErrorCode l_rval;

  moab::Range               l_faultFaces;
  moab::Range               l_mSets, l_tempRange;
  std::vector< moab::Tag >  l_tagsMesh;
  moab::Tag                 l_tagMat;
  std::string               l_tagMatName;

  //! Get Material-Type tag...
  l_rval = l_face->tag_get_tags( l_tagsMesh );
  assert( l_rval == moab::MB_SUCCESS );
  l_tagMat = l_tagsMesh[0];

  //! ...and double check that we have the right tag
  l_rval = l_face->tag_get_name( l_tagMat, l_tagMatName );
  assert( l_rval == moab::MB_SUCCESS );
  assert( l_tagMatName == "MATERIAL_SET" );

  //! Get all entity sets
  l_rval = l_face->get_entities_by_type( 0, moab::MBENTITYSET, l_mSets );
  assert( l_rval == moab::MB_SUCCESS );

  //! Read off the material type for each entity set
  std::vector< int > l_matData;
  l_matData.resize( l_mSets.size() );
  l_face->tag_get_data( l_tagMat, l_mSets, &l_matData[0] );
  assert( l_rval == moab::MB_SUCCESS );

  //! Check if entity set has material type 201 (rupture)
  for( size_t l_md = 0; l_md < l_matData.size(); l_md++ ) {
    if( l_matData[l_md] == 201 ) {
      l_tempRange.clear();
      l_rval = l_face->get_entities_by_handle( l_mSets[l_md], l_tempRange );
      l_faultFaces.merge( l_tempRange );
    }
  }

  return l_faultFaces;
}


// Create MOAB tags and tag each fault entity with array of doubles.
// Each tag is an array representing a different initial
// rupture parameter at each subtriangle and each fused run
// See diagram:
//
//  x------------------------------------------------------------x
//  |                         Fault Face                         |
//  x---------------------x------------------x-------------------x
//  |      subtri 1       |     subtri 2     |     subtri 3      |
//  x------x------x-------x------------------x-------------------x
//  | cfr1 | cfr2 | cfr3  |
//  x------x------x-------x
int faultAntn( const edge_v::io::Config & i_cfg,
                     moab_mesh          & i_mesh     )
{
  clock_t l_t = clock();

  if( i_cfg.m_faultInputFns.empty() ){
    std::cout << "No fault input files... skipping fault annotation"
              << std::endl;
    return 0;
  }
  if( i_mesh.m_intf == nullptr ){
    std::cout << "Failed." << std::endl;
    std::cerr << "Error: Attempted to annotate an uninitialized mesh (fault data)."
              << std::endl;
    exit( EXIT_FAILURE );
  }

  // Set parameters for constructing and querying fault model
  const unsigned short l_nCfr    = i_cfg.m_faultInputFns.size(); // Number of fused runs
  const unsigned short l_refine  = 2*i_cfg.m_tetRefinement+1;        // Subtri. refinement lvl
  const unsigned int   l_nSubTri = pow( l_refine, 2 );             // Subtris. per tet face
  const unsigned int   l_tagSize = l_nSubTri * l_nCfr;

  // Set up fault_model
  edge_v::fault::Model fault_model( i_cfg.m_faultInputFns );

  l_t = clock() - l_t;
  std::cout << " (" << (float) l_t / CLOCKS_PER_SEC << "s)" << std::endl;
  l_t = clock();

  // Create maps of tags and data arrays
  moab::ErrorCode                     l_rval;
  moab::Interface*                    l_iface = i_mesh.m_intf;
  std::map< std::string, real* >    l_qtyData;
  std::map< std::string, moab::Tag >  l_tags;

  for( const auto& l_qtyName : fault_model.m_qtyNames ){
    l_qtyData[ l_qtyName ] = new real[ l_tagSize ];

    l_rval = l_iface->tag_get_handle( l_qtyName.c_str(),
                                      l_tagSize,
                                      moab::MB_TYPE_DOUBLE,
                                      l_tags[ l_qtyName ],
                                      moab::MB_TAG_CREAT|moab::MB_TAG_SPARSE );
    assert( l_rval == moab::MB_SUCCESS );
  }

  // Loop over all fault faces, annotating mesh as we go
  Triangle                  l_tri;
  std::vector< Triangle >   l_subtriangles;
  moab::Range               l_faceVerts;
  real                      l_faceVertsCoords[9];
  FDatum                    l_fd;

  moab::Range l_faultFaces = getFaultFaces( i_mesh );

  std::cout << "Annotating mesh with fault stresses... " << std::flush;
  for( const auto& l_faceHandle : l_faultFaces ){
    l_faceVerts.clear();
    l_rval = l_iface->get_adjacencies( &l_faceHandle, 1, 0, false, l_faceVerts );
    assert( l_rval == moab::MB_SUCCESS );

    // Ensure that vertices in triangle are ordered according to increasing id
    // Can compare entity handles since entity type is the same for all verts
    assert( l_faceVerts[0] < l_faceVerts[1] );
    assert( l_faceVerts[1] < l_faceVerts[2] );
    l_rval = l_iface->get_coords( l_faceVerts, l_faceVertsCoords );
    assert( l_rval == moab::MB_SUCCESS );

    // Build a triangle from the vertices of triangle entity, and subdivide
    l_tri = Triangle( l_faceVertsCoords );
    l_subtriangles = l_tri.subdivide( l_refine );
    assert( l_nSubTri == l_subtriangles.size() );

    // Loop over each subtriangle and get approximate fault data
    for( unsigned int l_st = 0; l_st < l_nSubTri; l_st++ ){
      l_fd = fault_model.getDatumTri( l_subtriangles[ l_st ] );

      for( unsigned short l_cfr = 0; l_cfr < l_nCfr; l_cfr++ ){
        for( const auto& l_qtyName : fault_model.m_qtyNames ){
          l_qtyData[ l_qtyName ][ l_st * l_nCfr + l_cfr ] = l_fd[ l_qtyName ][ l_cfr ];
        }
      }
    }

    // Tag the mesh with each quantity
    for( const auto& l_qtyName : fault_model.m_qtyNames ){
      l_rval = l_iface->tag_set_data( l_tags[ l_qtyName ],
                                      &l_faceHandle,
                                      1,
                                      l_qtyData[ l_qtyName ] );
      assert( l_rval == moab::MB_SUCCESS);
    }
  }

  std::cout << "Done!" << std::flush;

  // Clean Up
  for( const auto& l_qtyName : fault_model.m_qtyNames ){
    delete[] l_qtyData[ l_qtyName ];
  }

  l_t = clock() - l_t;
  std::cout << " (" << (float) l_t / CLOCKS_PER_SEC << "s)" << std::endl;

  return 0;
}
