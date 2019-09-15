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
#ifndef EDGE_IO_CONFIG_H
#define EDGE_IO_CONFIG_H

#include <array>
#include <string>
#include <vector>
#include "constants.hpp"
#include "linalg/Domain.hpp"
#include "linalg/HalfSpace.hpp"
#include <pugixml.hpp>

namespace edge {
  namespace io {
    template< unsigned short TL_N_DIM >
    class ConfigDoms;

    class Config;
  }
}

/**
 * Handling of domain-objects in XML-configs.
 *
 * @paramt TL_N_DIM number of dimensions.
 **/
template< unsigned short TL_N_DIM >
class edge::io::ConfigDoms {
  public:
    //! types of values
    enum ValT {
      INT,
      UINT,
      F32,
      F64,
      BOOL,
      LONGLONG,
      ULONGLONG
    };

  private:
    /**
     * Converts a given datum, given as XML-node, to a native C-type.
     *
     * @param i_nd XML-node.
     * @param i_datT type of the datum (used for explicit conversion).
     * @return converted datum.
     * @parmat TL_T_DAT type of the datum (used for implicit conversion, if applicable).
     **/
    template< typename TL_T_DAT >
    static TL_T_DAT convert( pugi::xml_node const &i_nd,
                             ValT                  i_datT) {
      if(      i_datT == INT        ) return i_nd.text().as_int();
      else if( i_datT == UINT       ) return i_nd.text().as_uint();
      else if( i_datT == F32        ) return i_nd.text().as_float();
      else if( i_datT == F64        ) return i_nd.text().as_double();
      else if( i_datT == BOOL       ) return i_nd.text().as_bool();
      else if( i_datT == LONGLONG   ) return i_nd.text().as_llong();
      else if( i_datT == ULONGLONG  ) return i_nd.text().as_ullong();
      else {
        EDGE_LOG_FATAL;
        return 0;
      }
    }

  public:
    /**
     * Parses an xml-tree describing a domain.
     *
     * @param i_domNd xml-desc of the domain.
     * @param i_datN names of the data which will be read. [*][]: datum, [][*]: tree-hierarchy in the domain.
     * @param i_datT type of the data (used for explicit conversion).
     * @param io_dat data of the domain will be appended. For every given name (i_datN) a datum is added. If the type doesn't match the given template parameter, implicit conversions will be used.
     * @param o_dom will be set to domain.
     **/
    template< typename TL_T_REAL_DOM, typename TL_T_DAT >
    static void parse( pugi::xml_node                                const &i_domNd,
                       std::vector< std::vector< std::string > >     const &i_datN,
                       ValT                                                 i_datT,
                       std::vector< TL_T_DAT >                             &io_dat,
                       linalg::Domain<
                         TL_T_REAL_DOM,
                         TL_N_DIM,
                         linalg::HalfSpace
                       >                                                   &o_dom ) {

      // clear the domain
      o_dom.clear();

      // read the geometry
      for( pugi::xml_node l_hs = i_domNd.child("half_space"); l_hs; l_hs = l_hs.next_sibling("half_space") ) {
        // read origin and normal
        TL_T_REAL_DOM l_origin[TL_N_DIM];
        TL_T_REAL_DOM l_normal[TL_N_DIM];

        for( unsigned short l_di = 0; l_di < TL_N_DIM; l_di++ ) {
          l_origin[l_di] = (l_di==0) ? l_hs.child("origin").child("x").text().as_double() :
                           (l_di==1) ? l_hs.child("origin").child("y").text().as_double() :
                                       l_hs.child("origin").child("z").text().as_double();

          l_normal[l_di] = (l_di==0) ? l_hs.child("normal").child("x").text().as_double() :
                           (l_di==1) ? l_hs.child("normal").child("y").text().as_double() :
                                       l_hs.child("normal").child("z").text().as_double();
        }

        // create half-space
        linalg::HalfSpace< TL_T_REAL_DOM, TL_N_DIM> l_hsTmp( l_origin, l_normal );

        // add half-space to domain
        o_dom.add( l_hsTmp );
      }

      // read the data
      for( unsigned short l_va = 0; l_va < i_datN.size(); l_va++ )  {
        pugi::xml_node l_nd = i_domNd;

        for( unsigned short l_ch = 0; l_ch < i_datN[l_va].size(); l_ch++ ) {
          l_nd = l_nd.child( i_datN[l_va][l_ch].c_str() );
        }

        io_dat.push_back( convert< TL_T_DAT >(l_nd, i_datT) );
      }
    }
};

class edge::io::Config {
  private:
    /**
     * Prints the given string line by line.
     *
     * @param i_pre string which is added to the front of every line.
     * @param i_str string which is printed line by line.
     **/
    void printMultiLine( std::string i_pre,
                         std::string i_str );

    /**
     * Print build.
     * Runtime parameters don't have an influence here.
     * This is double checking to avoid changed build configs without a recompile.
     *
     * @param i_build build config which will be printed.
     **/
    void printBuild( pugi::xml_node i_build );

    /**
     * Prints the configuration.
     **/
    void printConfig();

  public:
    //! xml-document containing the config
    pugi::xml_document m_doc;

    /*
     * Mesh parameters
     */
    //! elements in x-dimension (regular only)
    int_el m_nElementsX;

    //! elements in y-dimension (regular only)
    int_el m_nElementsY;

    //! elements in z-dimension (regular only)
    int_el m_nElementsZ;

    //! size in x-dimension (regular only)
    real_mesh m_sizeX;

    //! size in y-dimension (regular only)
    real_mesh m_sizeY;

    //! size in z-dimension (regular only)
    real_mesh m_sizeZ;

    //! id of periodic boundaries (unstructured only)
    int m_periodic;

    //! id of the boundary conditions
    std::vector< int         > m_bndConId;

    //! names of the boundary conditions
    std::vector< std::string > m_bndConName;

    //! mesh read options
    std::string m_meshOptRead;

    //! mesh input file
    std::string m_meshFileIn;

    //! mesh output file
    std::string m_meshFileOut;

    /*
     * Simulation parameters
     */
    //! expression strings for the intial DOFs
    std::string m_initValsExprStrs[N_CRUNS];

    //! expressions strings for the reference DOFs
    std::string m_refValsExprStrs[N_CRUNS];

    //! end time of the simulations
    double m_endTime;

    //! sparse type of the wave field output
    int_spType m_waveFieldSpType = std::numeric_limits< int_spType >::max();

    //! plot type of the wave field output
    std::string m_waveFieldType;

    //! name of the output file for the wave field
    std::string m_waveFieldFile = "";

    //! interval of wave field output (max/2 to prevent inf when used in comparisons)
    double m_waveFieldInt =  std::numeric_limits< double >::max()/2;

    //! maximum synchronization interval (if sync point is reached otherwise before, this is ignored)
    double m_syncMaxInt = std::numeric_limits< double >::max()/2;

    //! type of the internal boundary output
    std::string m_iBndType;

    //! name of the output file for the internal boundary
    std::string m_iBndFile = "";

    //! interval of internal boundary output  (max/2 to prevent inf when used in comparisons)
    double m_iBndInt =  std::numeric_limits< double >::max()/2;

    //! type of the error norms
    std::string m_errorNormsType;

    //! file for xml output of the norms
    std::string m_errorNormsFile;

    //! receiver coordinates
    std::vector< std::array< real_mesh, 3 > > m_recvCrds[2];

    //! receiver names
    std::vector< std::string > m_recvNames[2];

    //! receivers frequency
    double m_recvFreq[2];

    //! path to receiver directory
    std::string m_recvPath[2];

    //! domains for sparse entity types, [0]: vertices, [1]: faces, [2]: elements
    std::vector< linalg::Domain< real_mesh, N_DIM, edge::linalg::HalfSpace > > m_spTypesDoms[3];

    //! values of the sparse entity types in the domains, [0]: vertices, [1]: faces, [2]: elements
    std::vector< int_spType > m_spTypesVals[3];

    /**
     * Constructor of the config.
     *
     * @param i_xmlPath path to xml file.
     **/
    Config( std::string i_xmlPath );
};

#endif
