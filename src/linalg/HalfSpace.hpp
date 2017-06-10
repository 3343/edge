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
 * Description of a half-space with inside operator.
 **/

#ifndef HALF_SPACE_HPP
#define HALF_SPACE_HPP

#include "io/logging.h"
#include <string>
#include "Geom.hpp"

namespace edge {
  namespace linalg {
    template <typename TL_T_REAL, unsigned short TL_N_DIM>
    class HalfSpace;
  }
}

/**
 * Half-space.
 *
 * @paramt TL_N_DIM dimension of the domain.
 * @paramt TL_T_REAL precision of the computations.
 **/
template <typename TL_T_REAL, unsigned short TL_N_DIM>
class edge::linalg::HalfSpace {
  private:
    TL_T_REAL m_origin[TL_N_DIM];
    TL_T_REAL m_normal[TL_N_DIM];
  public:
    /**
     * Sets the origin of the half-space.
     *
     * @param i_origin origin which will be set.
     **/
    void setOrigin( TL_T_REAL i_origin[TL_N_DIM] ) {
      for( unsigned short l_di = 0; l_di < TL_N_DIM; l_di++) {
        m_origin[l_di] = i_origin[l_di];
      }
    }

    /**
     * Sets the normal of the half-space.
     *
     * @param i_normal normal which will be set.
     **/
    void setNormal( TL_T_REAL i_normal[TL_N_DIM] ) {
      for( unsigned short l_di = 0; l_di < TL_N_DIM; l_di++) {
        m_normal[l_di] = i_normal[l_di];
      }

      // normalize normal
      if( TL_N_DIM == 2 ) {
        Geom::normalize2( m_normal );
      }
      else if( TL_N_DIM == 3 ) {
        Geom::normalize3( m_normal );
      }
      else EDGE_LOG_FATAL << "#dimensions not supported";
    }

    /**
     * Constructor of a half-space.
     *
     * @param i_origin origin/offset of the half-space.
     * @param i_normal normal of the half-space, direction is used for "inside".
     **/
    HalfSpace( TL_T_REAL i_origin[TL_N_DIM],
               TL_T_REAL i_normal[TL_N_DIM] ) {
      setOrigin( i_origin );
      setNormal( i_normal );
    }

    /**
     * Copy constructor.
     **/
    HalfSpace( HalfSpace const & i_hs ) {
      for( unsigned short l_di = 0; l_di < TL_N_DIM; l_di++) {
        m_origin[l_di] = i_hs.m_origin[l_di];
        m_normal[l_di] = i_hs.m_normal[l_di];
      }
    }

    /**
     * Generates a string with the half-space's info.
     *
     * @return string containing the origin and normal.
     **/
    std::string toString() {
      std::string l_str = "hs_o:";
      for( unsigned short l_di = 0; l_di < TL_N_DIM; l_di++ ) {
        l_str += " " + std::to_string( m_origin[l_di] );
      }
      l_str += ", hs_n:";
      for( unsigned short l_di = 0; l_di < TL_N_DIM; l_di++ ) {
        l_str += " " + std::to_string( m_normal[l_di] );
      }
      return l_str;
    }

    /**
     * Decides whether the given point is inside or outside of the half-space.
     *
     * 2D example: (n: normal, o: origin)
     *
     *   left part is in  * right part is out
     *                   *
     *                  *
     *                 *
     *                *
     *     n         *
     *         n    *
     *             o
     *            *
     *           *
     *          *
     *         *
     *        *
     *
     *   hyperplane itself: everything close than TL_ZERO_TOL is inside. 
     *
     * @param i_pt point which is considered to inside or outside.
     * @param i_tol zero tolerance used for derivation of in/out from the dot product. If negative the hyperplance is included, if positive it is excluded.
     **/
    bool inside( TL_T_REAL const i_pt[TL_N_DIM],
                 TL_T_REAL       i_tol = -1E-6 ) const {
      // determine directions from origin
      TL_T_REAL l_dir[TL_N_DIM];

      for( unsigned short l_di = 0; l_di < TL_N_DIM; l_di++ ) {
         l_dir[l_di] = i_pt[l_di] - m_origin[l_di];
      }

      // determine dot product
      TL_T_REAL l_dot = 0;
      if( TL_N_DIM == 2 ) {
        l_dot = Geom::sprod2( m_normal, l_dir );
      }
      else if( TL_N_DIM == 3 ) {
        l_dot = Geom::sprod3( m_normal, l_dir );
      }
      else EDGE_LOG_FATAL << "#dimensions not supported";
      if( l_dot > i_tol ) return true;
      return false;
    }
};

#endif
