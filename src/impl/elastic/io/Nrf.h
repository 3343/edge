/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
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
 * Reader for the NetCDF Rupture Format.
 **/

#ifndef NRF_H
#define NRF_H

#include <string> 
#include <vector>
#include <netcdf.h>
#include "io/logging.h"

namespace edge {
  namespace elastic {
    namespace io {
      template< unsigned short TL_N_DIM >
      class Nrf;
    }
  }
}

/**
 * Reader for the netCDF rupture format.
 *
 * Remark: The reader buffers the number of source terms, as passed to the constructor.
 *         If data of a source beyond the scope of the buffer is accessed,
 *         the data of the following sources according to the buffer size.
 *
 * @paramt number of dimensions.
 **/
template< unsigned short TL_N_DIM >
class edge::elastic::io::Nrf {
  private:
    //! over allocation for storing slip-rate values
    double m_srOvAll;

    //! buffer size of receivers kept in memory.
    const std::size_t m_bSize;

    // NRF subfault data type
    typedef  struct {
      //! first time at which the subfault is active.
      double tInit;

      //! time step of the subfault.
      double dt;

      //! Lame parameter mu.
      double mu;

      //! area (3D simulation) or length (2D simulation) of the subfault.
      double area;

      //! directions, [0][*]: fault-normal, [1][*]: first along-fault, [2][*]: second along-fault
      double dir[TL_N_DIM][TL_N_DIM];
    } t_nrfSub;

    //! NRF wrapper, where only the given number of sources are stored at once. 
    typedef struct {
      //! path of the file.
      std::string path;

      //! netCDF id of the file.
      int ncId;

      //! global number of sources.
      std::size_t nSrcsG;

      //! global number of samples in fault-normal (0), first along-fault (1), and second along-fault (2) direction.
      std::size_t nSplsG[TL_N_DIM];

      //! global id of first of the bSize sources currently stored.
      std::size_t gIdFb;

      //! id of the netCDF variable storing the coordinates of the sources
      int ncVarCrd;

      //! coordinates of the sources.
      double (*crds)[TL_N_DIM];

      //! id of the netCDF variable storing the subfaults
      int ncVarSub;

      //! subfaults described by the point sources in the buffer.
      t_nrfSub *sub;

      //! id of the netCDF variable storing the offsets
      int ncVarSrOff;

      //! true if the wrapper holds valid slip-rates for the buffer
      bool validSr;

      //! slip-rate offset w.r.t. to the buffer+1 (not the absolute offset as in the NRF file format).
      unsigned int (*srOff)[TL_N_DIM];

      //! maximum number of samples covered by allocated memory slip-rates
      std::size_t nSrMax[TL_N_DIM];

      //! id of the netCDF variable storing slip-rates in normal direction
      int ncVarSr[TL_N_DIM];

      //! slip-rates in normal, first along-fault, second along-fault
      double *sr[TL_N_DIM];
    } t_nrf;

    //! NRF wrappers.
    std::vector< t_nrf* > m_nrf;

    /**
     * Updates the source buffer for a given NRT data file.
     *
     * @param i_nrtId class-local id of the NRT data file.
     * @param i_gIdFb global id of the first source of the updated buffer.
     * @param i_sr also update the slip rates, otherwise only the meta-data is read.
     **/
    void updateBuf( unsigned short i_nrtId,
                    size_t         i_gIdFb,
                    bool           i_sr=true ) {
      // error code
      int l_err;

      // hyperslab for two-dimensional netCDF arrays
      std::size_t l_hypSlab[2][2];

      // get the pointer to the NRF data file
      t_nrf *l_nrf = m_nrf[i_nrtId];

      // store global id
      l_nrf->gIdFb = i_gIdFb;

      // buffer size
      EDGE_CHECK_LT( i_gIdFb, m_nrf[i_nrtId]->nSrcsG );
      std::size_t l_bSize   = std::min( m_bSize,   m_nrf[i_nrtId]->nSrcsG - i_gIdFb );
      std::size_t l_bSizeP1 = std::min( m_bSize+1, m_nrf[i_nrtId]->nSrcsG - i_gIdFb );

      // read source positions
      l_err = nc_get_vara(  m_nrf[i_nrtId]->ncId,
                            m_nrf[i_nrtId]->ncVarCrd,
                           &i_gIdFb,
                           &l_bSize,
                            m_nrf[i_nrtId]->crds );
      EDGE_CHECK_EQ( l_err, NC_NOERR ) << nc_strerror( l_err );

      // read subfaults
      l_err = nc_get_vara(  m_nrf[i_nrtId]->ncId,
                            m_nrf[i_nrtId]->ncVarSub,
                           &i_gIdFb,
                           &l_bSize,
                            m_nrf[i_nrtId]->sub );
      EDGE_CHECK_EQ( l_err, NC_NOERR ) << nc_strerror( l_err );

      // read sliprate offsets
      l_hypSlab[0][0] = i_gIdFb;
      l_hypSlab[0][1] = 0;
      l_hypSlab[1][0] = l_bSizeP1;
      l_hypSlab[1][1] = TL_N_DIM;

      l_err = nc_get_vara_uint(  m_nrf[i_nrtId]->ncId,
                                 m_nrf[i_nrtId]->ncVarSrOff,
                                 l_hypSlab[0],
                                 l_hypSlab[1],
                                 m_nrf[i_nrtId]->srOff[0] );
      EDGE_CHECK_EQ( l_err, NC_NOERR ) << nc_strerror( l_err );

      // reorder offsets, 0: normal, 1: first along-fault, 2: second along-fault
#ifdef PP_USE_OMP
#pragma omp parallel for
#endif
      for( std::size_t l_bu = 0; l_bu < l_bSizeP1; l_bu++ ) {
        unsigned int l_off[TL_N_DIM];
        for( unsigned short l_di = 0; l_di < TL_N_DIM; l_di++ )
          l_off[(l_di+1)%TL_N_DIM] = m_nrf[i_nrtId]->srOff[l_bu][l_di];
        for( unsigned short l_di = 0; l_di < TL_N_DIM; l_di++ )
          m_nrf[i_nrtId]->srOff[l_bu][l_di] = l_off[l_di];
      }

      // reorder directions
#ifdef PP_USE_OMP
#pragma omp parallel for
#endif
      for( std::size_t l_bu = 0; l_bu < l_bSize; l_bu++ ) {
        double l_dir[TL_N_DIM][TL_N_DIM];
        for( unsigned short l_d1 = 0; l_d1 < TL_N_DIM; l_d1++ ) {
          for( unsigned short l_d2 = 0; l_d2 < TL_N_DIM; l_d2++ ) {
            l_dir[(l_d1+1)%TL_N_DIM][l_d2] = m_nrf[i_nrtId]->sub[l_bu].dir[l_d1][l_d2];
          }
        }
        for( unsigned short l_d1 = 0; l_d1 < TL_N_DIM; l_d1++ ) {
          for( unsigned short l_d2 = 0; l_d2 < TL_N_DIM; l_d2++ ) {
            m_nrf[i_nrtId]->sub[l_bu].dir[l_d1][l_d2] = l_dir[l_d1][l_d2];
          }
        }
      }

      // set ghost entry for offset if this is the end of the sources
      if( l_bSize == l_bSizeP1 ) {
       for( unsigned short l_di = 0; l_di < TL_N_DIM; l_di++ )
         m_nrf[i_nrtId]->srOff[l_bSize][l_di] = m_nrf[i_nrtId]->nSplsG[l_di];
      }

      // set slip rates invalid and abort if not instructed other wise
      if( !i_sr ) {
        m_nrf[i_nrtId]->validSr = false;
        return;
      }

      // determine number of slip rate samples
      std::size_t l_nSr[TL_N_DIM];
      for( unsigned short l_di = 0; l_di < TL_N_DIM; l_di++ ) {
        l_nSr[l_di] =  m_nrf[i_nrtId]->srOff[l_bSize][l_di] - m_nrf[i_nrtId]->srOff[0][l_di];
      }

      // resize slip-rate buffers if current size is not sufficient
      for( unsigned short l_di = 0; l_di < TL_N_DIM; l_di++ ) {
        if( l_nSr[l_di] > m_nrf[i_nrtId]->nSrMax[l_di] ) {
          // free memory
          if( m_nrf[i_nrtId]->nSrMax[l_di] > 0 ) delete[] m_nrf[i_nrtId]->sr[l_di];

          // allocate new memory, increase accordingly to avoid too many allocs
          if( l_nSr[l_di] > 0 ) m_nrf[i_nrtId]->sr[l_di] = new double[ (std::size_t) (l_nSr[l_di]*m_srOvAll) ];

          // store new memory size
          m_nrf[i_nrtId]->nSrMax[l_di] = l_nSr[l_di]*m_srOvAll;
        }
      }

      // read slip-rates
      for( unsigned short l_di = 0; l_di < TL_N_DIM; l_di++ ) {
        if( l_nSr[l_di] > 0 ) {
          l_hypSlab[0][0] = m_nrf[i_nrtId]->srOff[0][l_di];
          l_hypSlab[0][1] = l_nSr[l_di];

          l_err = nc_get_vara_double( m_nrf[i_nrtId]->ncId,
                                      m_nrf[i_nrtId]->ncVarSr[l_di],
                                      l_hypSlab[0] + 0,
                                      l_hypSlab[0] + 1,
                                      m_nrf[i_nrtId]->sr[l_di] );
          EDGE_CHECK_EQ( l_err, NC_NOERR ) << nc_strerror( l_err );
        }
      }

      // convert offsets to local offsets
    #ifdef PP_USE_OMP
    #pragma omp parallel for
    #endif
      for( std::size_t l_bu = 1; l_bu < l_bSize+1; l_bu++ ) {
        for( unsigned short l_di = 0; l_di < TL_N_DIM; l_di++ ) {
          m_nrf[i_nrtId]->srOff[l_bu][l_di] -= m_nrf[i_nrtId]->srOff[0][l_di];
        }
      }

      // reset offset of first entires
      for( unsigned short l_di = 0; l_di < TL_N_DIM; l_di++ ) {
        m_nrf[i_nrtId]->srOff[0][l_di] = 0; 
      }
    }

    /**
     * Updates the buffer for all NRT data files.
     *
     * @param i_gIdFb global id of first source of the updated buffer.
     * @param i_sr also update the slip rates, otherwise only the meta-data is read.
     **/
    void updateBufAll( size_t i_gIdFb,
                       bool   i_sr=true ) {
      // iterate over all NRF files and update their buffers
      for( unsigned short l_fl = 0; l_fl < m_nrf.size(); l_fl++ ) {
        updateBuf( l_fl, i_gIdFb, i_sr );
      }
    }

    /**
     * Gets the local id of the point source in the buffer and updates the buffer
     * The buffer is updated, if the point source is not covered.
     *
     * @param i_kId id of the kinematic source description.
     * @param i_ptSrcG global id of the point source.
     **/
    std::size_t bufferId( unsigned short i_kId,
                          std::size_t    i_ptSrcG ) {
      // update buffer, if the current point source is not covered
      if(    m_nrf[i_kId]->gIdFb           <= i_ptSrcG
          && m_nrf[i_kId]->gIdFb + m_bSize  >  i_ptSrcG ) {}
      else updateBuf( i_kId, i_ptSrcG );

      // location of the source in the buffer
      return i_ptSrcG - m_nrf[i_kId]->gIdFb;
    }

  public:
    /**
     * Constructor which creates the NRF-reader.
     *
     * @param i_bize buffer size: number of receivers in memory during initialization.
     **/
    Nrf( unsigned int i_bSize=10000 ): m_srOvAll(1.25), m_bSize(i_bSize){};

    /**
     * Init the NRF-reader for the given file.
     * Might be called multiple times for different files.
     *
     * @param i_path path to the NRF-file.
     * @return id of the file within the NRF-reader.
     **/
    unsigned short init( std::string const &i_path ) {
      // error code
      int l_err;

      // dim ids
      int l_dId;

      // add a new netCDF file
      m_nrf.resize( m_nrf.size()+1 );

      // alocate memory for the NRF wrapper
      m_nrf.back() = new t_nrf;

      // allocate memory for the source coordinates of the sources in the buffer
      m_nrf.back()->crds = new double[ m_bSize ][ TL_N_DIM ];

      // allocate memory for the subfaults of the sources in the buffer
      m_nrf.back()->sub = new t_nrfSub[ m_bSize ];

      // allocate memory for the offsets
      m_nrf.back()->srOff = new unsigned int[ m_bSize+1 ][ TL_N_DIM ];

      // init number of slip-rate samples with zero
      for( unsigned short l_di = 0; l_di < TL_N_DIM; l_di++ ) {
        m_nrf.back()->nSrMax[l_di] = 0;
      }

      // save path
      m_nrf.back()->path = i_path;

      // open the file
      l_err = nc_open( i_path.c_str(), NC_NOWRITE, &(m_nrf.back()->ncId) );
      EDGE_CHECK_EQ( l_err, NC_NOERR ) << nc_strerror( l_err );

      // get the number of sources
      l_err = nc_inq_dimid( m_nrf.back()->ncId, "source", &l_dId );
      EDGE_CHECK_EQ( l_err, NC_NOERR ) << nc_strerror( l_err );
      l_err = nc_inq_dimlen( m_nrf.back()->ncId, l_dId, &(m_nrf.back()->nSrcsG) );
      EDGE_CHECK_EQ( l_err, NC_NOERR ) << nc_strerror( l_err );

      // get the number of samples in the different directions
      for( unsigned short l_di = 0; l_di < TL_N_DIM; l_di++ ) {
        std::string l_name;
        if(      l_di == 0 ) l_name = "sample"+std::to_string(TL_N_DIM);
        else if( l_di == 1 ) l_name = "sample1";
        else if( l_di == 2 ) l_name = "sample2";
        else EDGE_LOG_FATAL;

        l_err = nc_inq_dimid( m_nrf.back()->ncId, l_name.c_str(), &l_dId );
        EDGE_CHECK_EQ( l_err, NC_NOERR ) << nc_strerror( l_err );
        l_err = nc_inq_dimlen( m_nrf.back()->ncId, l_dId, (m_nrf.back()->nSplsG)+l_di );
        EDGE_CHECK_EQ( l_err, NC_NOERR ) << nc_strerror( l_err );
      }

      // get the netCDF id of the variables
      l_err = nc_inq_varid( m_nrf.back()->ncId, "centres", &(m_nrf.back()->ncVarCrd) );
      EDGE_CHECK_EQ( l_err, NC_NOERR ) << nc_strerror( l_err );

      l_err = nc_inq_varid( m_nrf.back()->ncId, "subfaults", &(m_nrf.back()->ncVarSub) );
      EDGE_CHECK_EQ( l_err, NC_NOERR ) << nc_strerror( l_err );

      l_err = nc_inq_varid( m_nrf.back()->ncId, "sroffsets", &(m_nrf.back()->ncVarSrOff) );
      EDGE_CHECK_EQ( l_err, NC_NOERR ) << nc_strerror( l_err );

      // get the slip rate variables
      for( unsigned short l_di = 0; l_di < TL_N_DIM; l_di++ ) {
        std::string l_name;
        if(      l_di == 0 ) l_name = "sliprates"+std::to_string(TL_N_DIM);
        else if( l_di == 1 ) l_name = "sliprates1";
        else if( l_di == 2 ) l_name = "sliprates2";
        else EDGE_LOG_FATAL;

        l_err = nc_inq_varid( m_nrf.back()->ncId, l_name.c_str(), &(m_nrf.back()->ncVarSr[l_di]) );
        EDGE_CHECK_EQ( l_err, NC_NOERR ) << nc_strerror( l_err );
      }

      // init the buffer
      updateBufAll( 0 );

      return m_nrf.size();
    }

    /**
     * Generates a vector of strings (one per file) containing information on the NRF data of the reader.
     *
     * @return vector with strings. one string per NRF.
     **/
    std::vector< std::string > toString() const {
      // vector of strings which gets returned
      std::vector< std::string > l_re;

      // assemble info
      for( unsigned short l_fl = 0; l_fl < m_nrf.size(); l_fl++ ) {
        std::string l_str =   "NRF#" + std::to_string(l_fl) + " "
                                + "- #srcs: " + std::to_string(m_nrf[l_fl]->nSrcsG)
                                + ", #samples n: " + std::to_string(m_nrf[l_fl]->nSplsG[0]);
        l_str += (TL_N_DIM > 1) ? ", #samples s: " + std::to_string(m_nrf[l_fl]->nSplsG[1]) : "";
        l_str += (TL_N_DIM > 2) ? ", #samples t: " + std::to_string(m_nrf[l_fl]->nSplsG[2]) : "";

        l_re.push_back( l_str );
      }

      return l_re;
    }

    /**
     * Gets the number of inputs which haven been initialized.
     *
     * @return number of inputs.
     *
     **/
    unsigned short nIn() const {
      return m_nrf.size();
    }

    /**
     * Gets the global number of sources for the given input.
     * 
     * @param i_kId id of the kinematic source for which the number of sources is queried.
     * @return number of sources for the given input.
     **/
    std::size_t nSrcsG( unsigned short i_kId ) const {
      EDGE_CHECK_LT( i_kId, m_nrf.size() );

      return m_nrf[i_kId]->nSrcsG;
    }

    /**
     * Gets the global number of samples.
     *
     * @param i_kId id of the kinematic source for which the number of samples is queried.
     * @param i_dir slip-direction for which the number of samples is queried.
     **/
    std::size_t nSplsG( unsigned short i_i_kId,
                        unsigned short i_dir ) const {
      EDGE_CHECK_LT( i_i_kId,  m_nrf.size() );
      EDGE_CHECK_LT( i_dir, TL_N_DIM );

      return m_nrf[i_i_kId]->nSplsG[i_dir];
    }

    /**
     * Gets the coordinates for the source in the given input.
     *
     * @param i_kId id of the kinematic source for which the source coordinates are queried.
     * @param o_crds memory to whic the source coordinates are written to, [*][]: source, [][*] dimension.
     * @param i_first first source for which the source coordinates are returned.
     * @param i_size number of sources for which the source coordinates are returned, if 0 all sources after the first are returned.
     **/
    void getSrcCrds( unsigned short   i_kId,
                     double         (*o_crds)[TL_N_DIM],
                     std::size_t      i_first=0,
                     std::size_t      i_size=0 ) const {
      // check that we are inbound
      EDGE_CHECK_LT( i_kId, m_nrf.size() );
      std::size_t l_size = m_nrf[i_kId]->nSrcsG;
      if( i_size > 0 ) l_size = i_size;
      EDGE_CHECK_LT( i_first, m_nrf[i_kId]->nSrcsG );
      EDGE_CHECK_LT( i_first+l_size-1, m_nrf[i_kId]->nSrcsG );

      // get the coordinates
      int l_err = nc_get_vara(  m_nrf[i_kId]->ncId,
                                m_nrf[i_kId]->ncVarCrd,
                               &i_first,
                               &l_size,
                                o_crds[0] );
      EDGE_CHECK_EQ( l_err, NC_NOERR ) << nc_strerror( l_err );
    }

    /**
     * Gets the global slip-rate offsets for the given kinematic source description.
     *
     * If no size (i_size equals to 0) the offsets for all sources are written to the output-array.
     * Additional, the last entry is set to a ghost entry, containing the total number of slip rates.
     * Thus if n_s is the number of sources, n_s+1 entries are written. This behavior can be requested
     * explicitly by setting i_first to 0 and i_size to n_s+1.
     *
     * @param i_kId id of the kinematic source description.
     * @param i_dir slip-direction for which the offsets are obtained, 0: normal, 1: first along-fault, 2: second along-fault
     * @param o_srOff will be set to slip-rate offsets for the given kinematic source and direction.
     * @param i_first first point source for which the offset is queried.
     * @param i_size number of point sources for which the offset is queried.
     **/
    void getOffSetsG( unsigned short  i_kId,
                      unsigned short  i_dir,
                      unsigned int   *o_srOff,
                      std::size_t     i_first=0,
                      std::size_t     i_size=0 ) const {
      // check that we are inbound
      EDGE_CHECK_LT( i_kId, m_nrf.size() );
      EDGE_CHECK_LT( i_dir, TL_N_DIM );

      // default to all sources n_s+1 (including total number of sliprates)
      std::size_t l_size = m_nrf[i_kId]->nSrcsG;
      if( l_size > 0 ) l_size++;

      // use user-input for size if given
      if( i_size > 0 ) {
        l_size = i_size;
      }
      // fall-back to all sources
      else {
        EDGE_CHECK( i_first < l_size );
        l_size -= i_first;
      }

      // adjust netCDF-direction to match EDGE's convention, always storing the normal first
      unsigned short l_dir = (i_dir + TL_N_DIM-1)%TL_N_DIM;

      // hyperslab for two-dimensional netCDF arrays
      std::size_t l_hypSlab[2][2];

      l_hypSlab[0][0] = i_first;
      l_hypSlab[0][1] = l_dir;
      l_hypSlab[1][0] = l_size;
      l_hypSlab[1][1] = 1;

      // read data
      int l_err = nc_get_vara_uint(  m_nrf[i_kId]->ncId,
                                     m_nrf[i_kId]->ncVarSrOff,
                                     l_hypSlab[0],
                                     l_hypSlab[1],
                                     o_srOff );
      EDGE_CHECK_EQ( l_err, NC_NOERR ) << nc_strerror( l_err );
    }

    /**
     * Gets the onset time of the given point source.
     *
     * @param i_kId id of the kinematic source description.
     * @param i_ptSrcG global id of the point source.
     * @return onset time.
     **/
    double getOnSet( unsigned short i_kId,
                     std::size_t    i_ptSrcG ) {
      std::size_t l_bId = bufferId( i_kId, i_ptSrcG );
      return m_nrf[i_kId]->sub[l_bId].tInit;
    }

    /**
     * Gets the temporal distance of two source subsequent samples for the given point source.
     *
     * @param i_kId id of the kinematic source description.
     * @param i_ptSrcG global id of the point source.
     * @return temporal distance.
     **/
    double getDt( unsigned short i_kId,
                  std::size_t    i_ptSrcG ) {
      std::size_t l_bId = bufferId( i_kId, i_ptSrcG );
      return m_nrf[i_kId]->sub[l_bId].dt;
    }

    /**
     * Gets the Lame parameter mu for the given point source.
     *
     * @param i_kId id of the kinematic source description.
     * @param i_ptSrcG global id of the point source.
     * @return Lame parameter mu.
     **/
    double getMu( unsigned short i_kId,
                  std::size_t    i_ptSrcG ) {
      std::size_t l_bId = bufferId( i_kId, i_ptSrcG );
      return m_nrf[i_kId]->sub[l_bId].mu;
    }

    /**
     * Gets the area (3D simulations) or length (2D simulations) of the subfault for the given point source.
     *
     * @param i_kId id of the kinematic source description.
     * @param i_ptSrcG global id of the point source.
     * @return area/length.
     **/
    double getA( unsigned short i_kId,
                 std::size_t    i_ptSrcG ) {
      std::size_t l_bId = bufferId( i_kId, i_ptSrcG );
      return m_nrf[i_kId]->sub[l_bId].area;
    }

    /**
     * Gets the slip-direction for the subfault associated to the given point source
     *
     * @param i_kId id of the kinematic source description.
     * @param i_ptSrcG global id of the point source.
     * @param o_sds will be set to directions. [*][]: direction, [][*]: entries of the direction.
     **/
    void getSds( unsigned short i_kId,
                 std::size_t    i_ptSrcG,
                 double         o_sds[TL_N_DIM][TL_N_DIM] ) {
      std::size_t l_bId = bufferId( i_kId, i_ptSrcG );
      for( unsigned short l_d1 = 0; l_d1 < TL_N_DIM; l_d1++ )
        for( unsigned short l_d2 = 0; l_d2 < TL_N_DIM; l_d2++ )
          o_sds[l_d1][l_d2] = m_nrf[i_kId]->sub[l_bId].dir[l_d1][l_d2];
    }

    /**
     * Gets the slip rates for the given point source.
     * TODO: enable fused sources
     *
     * @param i_kId id of the kinematic source description.
     * @param i_ptSrcG global id of the point source.
     * @param i_sd slip direction.
     * @param o_srs will be set to slip-rates.
     *
     * @paramt TL_T_REAL floating-point precision of the requested slip-rates.
     **/
    template < typename TL_T_REAL >
    void getSrs( unsigned short  i_kId,
                 std::size_t     i_ptSrcG,
                 unsigned short  i_sd,
                 TL_T_REAL      *o_srs ) {
      std::size_t l_bId = bufferId( i_kId, i_ptSrcG );

      std::size_t l_first = m_nrf[i_kId]->srOff[l_bId][i_sd];
      std::size_t l_size  = m_nrf[i_kId]->srOff[l_bId+1][i_sd] - l_first;

      // check that the size matches the total number of samples
      EDGE_CHECK_LE( l_size, m_nrf[i_kId]->nSplsG[i_sd] );
      for( std::size_t l_sr = 0; l_sr < l_size; l_sr++ )
        o_srs[l_sr] = m_nrf[i_kId]->sr[i_sd][l_sr+l_first];
    }

    /**
     * Gets the coordinates of a source.
     *
     * @param i_kId id of the kinematic source description.
     * @param i_ptSrcG global id of the point source.
     * @param o_crds will be set to coordiantes.
     **/
    void getSrcCrds( unsigned short i_kId,
                     std::size_t    i_ptSrcG,
                     double         o_crds[TL_N_DIM] ) {
      std::size_t l_bId = bufferId( i_kId, i_ptSrcG );

      for( unsigned short l_di = 0; l_di < TL_N_DIM; l_di++ )
        o_crds[l_di] = m_nrf[i_kId]->crds[l_bId][l_di];
    }


    /**
     * Destructor which closes all netCDF file streams and deallocated the occupied memory.
     **/
    ~Nrf() {
      int l_err;

      for( unsigned short l_fl = 0; l_fl < m_nrf.size(); l_fl++ ) {
        // close the file
        l_err = nc_close( m_nrf[l_fl]->ncId );
        EDGE_CHECK_EQ( l_err, NC_NOERR ) << nc_strerror( l_err );

        // free memory of the source position
        delete[] m_nrf[l_fl]->crds;

        // free memory of the subfaults of the sources in the buffer
        delete[] m_nrf[l_fl]->sub;

        // free memory of the offsets
        delete[] m_nrf[l_fl]->srOff;

        // free memory of the slip rates
        for( unsigned short l_di = 0; l_di < TL_N_DIM; l_di++ )
          if( m_nrf[l_fl]->nSrMax[l_di] > 0 ) delete[] m_nrf[l_fl]->sr[l_di];

        // free the wrapper
        delete m_nrf[l_fl];
      }
    }
};
#endif
