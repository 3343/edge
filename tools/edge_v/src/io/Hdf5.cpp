/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (breuer AT mytum.de)
 *
 * @section LICENSE
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
 * HDF-5 interface.
 **/
#include "Hdf5.h"

edge_v::io::Hdf5::Hdf5( std::string const & i_path,
                        bool                i_readOnly,
                        std::string const & i_group ) {
  // silence error printing
  herr_t l_err = H5Eset_auto( H5P_DEFAULT,
                              NULL,
                              NULL );
  EDGE_V_CHECK_GE( l_err, 0 );

  // create/open the file
  if( i_readOnly ) {
    m_fileId = H5Fopen( i_path.c_str(),
                        H5F_ACC_RDONLY,
                        H5P_DEFAULT );
  }
  else {
    m_fileId = H5Fopen( i_path.c_str(),
                        H5F_ACC_RDWR,
                        H5P_DEFAULT );
  }

  m_groupStr = i_group;
}

edge_v::io::Hdf5::~Hdf5() {
  herr_t l_err = H5Fclose( m_fileId );
  EDGE_V_CHECK_GE( l_err, 0 );
}

bool edge_v::io::Hdf5::exists( std::string const & i_path ) const {
  // open group
  hid_t l_group = openCreateGroup();

  htri_t l_ex = H5Lexists( l_group,
                           i_path.c_str(),
                           H5P_DEFAULT );
  EDGE_V_CHECK_GE( l_ex, 0 );

  // close group
  herr_t l_err = H5Gclose( l_group );
  EDGE_V_CHECK_GE( l_err, 0 );

  return l_ex > 0;
}

std::size_t edge_v::io::Hdf5::nVas( std::string const & i_name ) const {
  // open group
  hid_t l_group = openCreateGroup();

  // open dataset
  hid_t l_dset = H5Dopen2( l_group,
                           i_name.c_str(),
                           H5P_DEFAULT );
  EDGE_V_CHECK_GE( l_dset, 0 );

  // get space
  hid_t l_space = H5Dget_space( l_dset );
  EDGE_V_CHECK_GE( l_space, 0 );

  hsize_t l_nEns;
  int l_nDims = H5Sget_simple_extent_dims( l_space,
                                           &l_nEns,
                                           NULL );
  EDGE_V_CHECK_EQ( l_nDims, 1 );

  // close set, space and group
  herr_t l_err = H5Dclose( l_dset );
  EDGE_V_CHECK_GE( l_err, 0 );
  l_err = H5Sclose( l_space );
  EDGE_V_CHECK_GE( l_err, 0 );
  l_err = H5Gclose( l_group );
  EDGE_V_CHECK_GE( l_err, 0 );

  return l_nEns;
}

hid_t edge_v::io::Hdf5::openCreateGroup() const {
  // try to open the group
  hid_t l_id = H5Gopen( m_fileId,
                        m_groupStr.c_str(),
                        H5P_DEFAULT );

  // create the group if it doesn't exist
  if( l_id < 0 ) {
    l_id = H5Gcreate( m_fileId,
                      m_groupStr.c_str(),
                      H5P_DEFAULT,
                      H5P_DEFAULT,
                      H5P_DEFAULT );
  }

  return l_id;
}

void edge_v::io::Hdf5::get( std::string const & i_name,
                            hid_t               i_memType,
                            void              * o_data ) const {
  // open group
  hid_t l_group = openCreateGroup();

  // open dataset
  hid_t l_dset = H5Dopen2( l_group,
                           i_name.c_str(),
                           H5P_DEFAULT );
  EDGE_V_CHECK_GE( l_dset, 0 );

  // read data
  herr_t l_err = H5Dread( l_dset,
                          i_memType,
                          H5S_ALL,
                          H5S_ALL,
                          H5P_DEFAULT,
                          o_data );
  EDGE_V_CHECK_GE( l_err, 0 );
  
  // close set and group
  l_err = H5Dclose( l_dset );
  EDGE_V_CHECK_GE( l_err, 0 );
  l_err = H5Gclose( l_group );
  EDGE_V_CHECK_GE( l_err, 0 );
}

void edge_v::io::Hdf5::set( std::string const & i_name,
                            std::size_t         i_nValues,
                            void        const * i_data,
                            hid_t               i_memType,
                            hid_t               i_fileType ) const {
  // create dataspace
  hsize_t l_nVas = i_nValues;
  hid_t l_dspace = H5Screate_simple( 1,
                                     &l_nVas,
                                     NULL );
  hid_t l_group = openCreateGroup();

  // create dataset
  hid_t l_dset = H5Dcreate( l_group,
                            i_name.c_str(),
                            i_fileType ,
                            l_dspace,
                            H5P_DEFAULT,
                            H5P_DEFAULT,
                            H5P_DEFAULT );

  // write data
  herr_t l_err = H5Dwrite( l_dset,
                           i_memType,
                           H5S_ALL,
                           H5S_ALL,
                           H5P_DEFAULT,
                           i_data );
  EDGE_V_CHECK_GE( l_err, 0 );

  // close set, space and group
  l_err = H5Dclose( l_dset );
  EDGE_V_CHECK_GE( l_err, 0 );
  l_err = H5Sclose( l_dspace );
  EDGE_V_CHECK_GE( l_err, 0 );
  l_err = H5Gclose( l_group );
  EDGE_V_CHECK_GE( l_err, 0 );
}

void edge_v::io::Hdf5::set( std::string    const & i_name,
                            std::size_t            i_nValues,
                            unsigned short const * i_data ) const {
  set( i_name,
       i_nValues,
       i_data,
       H5T_NATIVE_USHORT,
       H5T_STD_U16LE );
}

void edge_v::io::Hdf5::set( std::string const & i_name,
                            std::size_t         i_nValues,
                            std::size_t const * i_data ) const {
  // convert to unsigned long
  unsigned long *l_data = new unsigned long[i_nValues];
  convert( i_nValues,
           i_data,
           l_data );

  // store
  set( i_name,
       i_nValues,
       l_data,
       H5T_NATIVE_ULONG,
       H5T_STD_I64LE );

  // free memory
  delete[] l_data;
}

void edge_v::io::Hdf5::set( std::string    const & i_name,
                            std::size_t            i_nValues,
                            double         const * i_data ) const {
  set( i_name,
       i_nValues,
       i_data,
       H5T_NATIVE_DOUBLE,
       H5T_IEEE_F64LE );
}

void edge_v::io::Hdf5::get( std::string const & i_name,
                            unsigned short    * o_data ) const {
  get( i_name,
       H5T_NATIVE_USHORT,
       o_data );
}

void edge_v::io::Hdf5::get( std::string const & i_name,
                            std::size_t       * o_data ) const {
  // alloc memory
  std::size_t l_nVas = nVas( i_name );
  unsigned long *l_data = new unsigned long[l_nVas];

  get( i_name,
       H5T_NATIVE_ULONG,
       l_data );

  // convert
  convert( l_nVas,
           l_data,
           o_data );

  // free memory
  delete[] l_data;
}

void edge_v::io::Hdf5::get( std::string const & i_name,
                            double            * o_data ) const {
  get( i_name,
       H5T_NATIVE_DOUBLE,
       o_data );
}