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
 * Description of a domain.
 **/

#ifndef DOMAIN_HPP
#define DOMAIN_HPP

#include <cstring>

namespace edge {
  namespace linalg {
    template< typename TL_T_REAL, unsigned short TL_N_DIM, template<typename, unsigned short> class TL_T_GEO_OBJ >
    class Domain;
  }
}

/**
 * Spatial domain.
 *
 * @paramt TL_T_REAL precision of the computations.
 * @paramt TL_T_GEO_OBJ type of the geometric objects build the domain.
 * @paramt TL_N_DIM dimension of the domain.
 **/
template< typename                                 TL_T_REAL,
          unsigned short                           TL_N_DIM,
          template<typename, unsigned short> class TL_T_GEO_OBJ >
class edge::linalg::Domain {
  private:
    std::vector< TL_T_GEO_OBJ< TL_T_REAL, TL_N_DIM > > m_geoObjs;

  public:
    /**
     * Constructor.
     **/
    Domain(){};

    /**
     * Clears the domain by removing all geometric objects.
     **/
    void clear() { m_geoObjs.clear(); };

    /**
     * Generates a vector of sting containing the domain's objects info.
     *
     * @return vector containing a string describing each of the domain's objects.
     **/
    std::vector< std::string > toString() {
      std::vector< std::string > l_strs;

      for( std::size_t l_ob = 0; l_ob < m_geoObjs.size(); l_ob++ ) {
        l_strs.push_back( m_geoObjs[l_ob].toString() );
      }
      return l_strs;
    }

    /**
     * Adds a geometric object to the domain.
     *
     * @param i_obj object which will be added.
     **/
    void add( TL_T_GEO_OBJ< TL_T_REAL, TL_N_DIM > const & i_obj ) {
      m_geoObjs.push_back( i_obj );
    }

    /**
     * Determines if the given point is inside the domain, given by the intersection of its geometric objects.
     *
     * Remark: Per definition domains without geometric object always return true;
     *
     * @param i_pt point.
     * @return true if the point is inside all geometric objects, false otherwise.
     *
     * @paramt TL_T_REAL_FUN precision of the point as passed to the function.
     **/
    template < typename TL_T_REAL_FUN >
    bool inside( TL_T_REAL_FUN const i_pt[TL_N_DIM] ) const {
      // convert type
      TL_T_REAL l_pt[TL_N_DIM];
      for( unsigned short l_di = 0; l_di < TL_N_DIM; l_di++ ) l_pt[l_di] = i_pt[l_di];

      for( std::size_t l_ob = 0; l_ob < m_geoObjs.size(); l_ob++ ) {
        if( m_geoObjs[l_ob].inside( l_pt ) == false ) return false;
      }
      return true;
    }
};

#endif
