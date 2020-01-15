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
#ifndef EDGE_V_IO_HDF5_H
#define EDGE_V_IO_HDF5_H

#include <hdf5.h>
#include <string>
#include "io/logging.h"

namespace edge_v {
  namespace io {
    class Hdf5;
  }
}

/**
 * Hdf5 interface.
 **/
class edge_v::io::Hdf5 {
  private:
    //! file id
    hid_t m_fileId = 0;

    //! string representation of the group
    std::string m_groupStr;

    /**
     * Converts the given array.
     *
     * @param i_nValues number of values.
     * @param i_data data which is converted.
     * @param o_data will be set to output data.
     *
     * @paramt TL_T_IN type of the input data.
     * @paramt TL_T_OUT type of the output data.
     **/
    template< typename TL_T_IN,
              typename TL_T_OUT >
    static void convert( std::size_t         i_nValues,
                         TL_T_IN     const * i_data,
                         TL_T_OUT          * o_data ) {
#ifdef PP_USE_OMP
#pragma omp parallel for
#endif
      for( std::size_t l_va = 0; l_va < i_nValues; l_va++ ) {
#if !defined(__clang__) && !defined(__INTEL_COMPILER)
        // check bounds
        if( !std::is_signed<TL_T_OUT>() ) {
          EDGE_V_CHECK_GE( i_data[l_va], std::numeric_limits< TL_T_OUT >::lowest() );
        }
        EDGE_V_CHECK_LE( i_data[l_va], std::numeric_limits< TL_T_OUT >::max()    );
#endif

        o_data[l_va] = i_data[l_va];
      }
    }

    /**
     * Opens this object group if it exists.
     * Otherwise: create and open.
     *
     * @return hid_t identifier of the group.
     **/
    hid_t openCreateGroup() const;

    /**
     * Sets the data to a dataset.
     *
     * @param i_name name of the dataset.
     * @param i_nValues number of values which are set.
     * @param i_data data which is is set.
     * @param i_memType type in memory.
     * @param i_fileType type in the HDF5 file.
     **/
    void set( std::string const & i_name,
              std::size_t         i_nValues,
              void        const * i_data,
              hid_t               i_memType,
              hid_t               i_fileType ) const;

    /**
     * Gets the data from a dataset.
     *
     * @param i_name name of the dataset.
     * @param i_memType type in memory.
     * @param o_data array to which the data is copied.
     **/
    void get( std::string const & i_name,
              hid_t               i_memType,
              void              * o_data ) const;

  public:
    /**
     * Constructor.
     *
     * @param i_path path to the HDF5 file which is opened.
     * @param i_readOnly if true the file is opened read only, r+w otherwise.
     * @param i_group group under which the data is stored.
     **/
    Hdf5( std::string const & i_path,
          bool                i_readOnly = true,
          std::string const & i_group = "/edge_v" );

    /**
     * Destructor
     **/
    ~Hdf5();

    /**
     * Tests if the given link in the interface's group exists.
     *
     * @param i_path path of the link.
     * @return true if it exists, false if not.
     **/
    bool exists( std::string const & i_path ) const;

    /**
     * Gets the number of values for the given data set.
     *
     * @param i_name name of the dataset.
     * @return number of values.
     **/
    std::size_t nVas( std::string const & i_name ) const;

    /**
     * Sets the data to a dataset.
     *
     * @param i_name name of the dataset.
     * @param i_nValues number of values which are set.
     * @param i_data data which is is set.
     **/
    void set( std::string    const & i_name,
              std::size_t            i_nValues,
              unsigned short const * i_data ) const;

    /**
     * Sets the data to a dataset.
     *
     * @param i_name name of the dataset.
     * @param i_nValues number of values which are set.
     * @param i_data data which is is set.
     **/
    void set( std::string const & i_name,
              std::size_t         i_nValues,
              std::size_t const * i_data ) const;

    /**
     * Sets the data to a dataset.
     *
     * @param i_name name of the dataset.
     * @param i_nValues number of values which are set.
     * @param i_data data which is is set.
     **/
    void set( std::string const & i_name,
              std::size_t         i_nValues,
              float       const * i_data ) const;

    /**
     * Sets the data to a dataset.
     *
     * @param i_name name of the dataset.
     * @param i_nValues number of values which are set.
     * @param i_data data which is is set.
     **/
    void set( std::string const & i_name,
              std::size_t         i_nValues,
              double      const * i_data ) const;

    /**
     * Gets the data from a dataset.
     *
     * @param i_name name of the dataset.
     * @param o_data array to which the data is written.
     **/
    void get( std::string    const & i_name,
              unsigned short       * o_data ) const;

    /**
     * Gets the data from a dataset.
     *
     * @param i_name name of the dataset.
     * @param o_data array to which the data is written.
     **/
    void get( std::string const & i_name,
              std::size_t       * o_data ) const;

    /**
     * Gets the data from a dataset.
     *
     * @param i_name name of the dataset.
     * @param o_data array to which the data is written.
     **/
    void get( std::string const & i_name,
              float             * o_data ) const;

    /**
     * Gets the data from a dataset.
     *
     * @param i_name name of the dataset.
     * @param o_data array to which the data is written.
     **/
    void get( std::string const & i_name,
              double            * o_data ) const;
};

#endif