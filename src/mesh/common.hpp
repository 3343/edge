/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2021, Friedrich Schiller University Jena
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
 * Common functions for the mesh
 **/
#ifndef EDGE_MESH_COMMON_HPP
#define EDGE_MESH_COMMON_HPP

#include "linalg/Geom.hpp"
#include "io/logging.h"
#include "monitor/instrument.hpp"
#include "constants.hpp"
#include "data/EntityLayout.type"
#include "linalg/Matrix.h"
#include <cassert>
#include <cmath>
#include <limits>

namespace edge {
  namespace mesh {
    template < t_entityType TL_T_EL >
    class common;
  }
}

template < t_entityType TL_T_EL >
class edge::mesh::common {
  private:
    //! number of dimensions
    static const unsigned short TL_N_DIS = C_ENT[TL_T_EL].N_DIM;

    //! face type
    static const t_entityType TL_T_FA = C_ENT[TL_T_EL].TYPE_FACES;

    //! number of face vertices
    static unsigned short const TL_N_FA_VES = C_ENT[TL_T_FA].N_VERTICES;

    //! number of element vertices
    static unsigned short const TL_N_EL_VES = C_ENT[TL_T_EL].N_VERTICES;

    //! number of element faces
    static unsigned short const TL_N_EL_FAS = C_ENT[TL_T_EL].N_FACES;

  public:
    /**
     * Gets the coordinates of an elements vertices based on connectivity information.
     *
     * @param i_el id of the entity.
     * @param i_elVe entities' vertices.
     * @param i_veChars characteristics of the vertices.
     * @param o_veCrds will be set to coordinates of the vertices.
     *
     * @paramt TL_T_LID integral type of the local ids.
     * @paramt TL_T_REAL floating point type.
     * @paramt TL_T_VE_CHARS type of the vertex characteristics, offering member .coords.
     **/
    template< typename TL_T_LID,
              typename TL_T_REAL,
              typename TL_T_VE_CHARS >
    static void getElVeCrds(       TL_T_LID        i_el,
                             const TL_T_LID      (*i_elVe)[TL_N_EL_VES],
                             const TL_T_VE_CHARS  *i_veChars,
                                   TL_T_REAL       o_veCrds[TL_N_DIS][TL_N_EL_VES] ) {
      // get elements vertices
      for( unsigned short l_ve = 0; l_ve < TL_N_EL_VES; l_ve++ ) {
        TL_T_LID l_veId       = i_elVe[i_el][l_ve];

        for( unsigned short l_di = 0; l_di < TL_N_DIS; l_di++ ) {
          o_veCrds[l_di][l_ve] = i_veChars[l_veId].coords[l_di];
        }
      }
    }

    /**
     * Prints a summary of the present neighboring relations.
     *
     * The index is given by:
     *     local_face_id-neighboring_face_id-neighboring_vertex_id
     *
     * @param i_nEls number of elements.
     * @param i_faIdElFaEl face relations.
     * @param i_veIdElFaEl vertex relations.
     **/
    static void printNeighRel(  std::size_t            i_nEls,
                                const unsigned short * i_fIdElFaEl,
                                const unsigned short * i_vIdElFaEl ) {
      t_entityType l_elType = TL_T_EL;
      t_entityType l_faType = TL_T_FA;

      // jump for the neighboring face (1 if not 3D)
      unsigned short l_veJump = 1;
      if( N_DIM == 3 ) {
        l_veJump = C_ENT[l_faType].N_VERTICES;
      }

      // compute number of options
      unsigned short l_nOpts = C_ENT[l_elType].N_FACES*C_ENT[l_elType].N_FACES*l_veJump;

      // reset options
      unsigned long *l_optsLo = new unsigned long[l_nOpts];
      for( unsigned short l_opt = 0; l_opt < l_nOpts; l_opt++ ) l_optsLo[l_opt] = 0;

      for( int_el l_el = 0; l_el < i_nEls; l_el++ ) {
        for( unsigned short l_fa = 0; l_fa < C_ENT[l_elType].N_FACES; l_fa++ ) {
          if( i_fIdElFaEl[l_el*C_ENT[l_elType].N_FACES+l_fa] < C_ENT[l_elType].N_FACES ) {
            unsigned short l_locOpt = l_fa * C_ENT[l_elType].N_FACES*l_veJump
                                      +
                                      i_fIdElFaEl[l_el*C_ENT[l_elType].N_FACES+l_fa] * l_veJump
                                      +
                                      i_vIdElFaEl[l_el*C_ENT[l_elType].N_FACES+l_fa];
            l_optsLo[l_locOpt]++;
          }
        }
      }

      unsigned long *l_optsGl;
#ifdef PP_USE_MPI
      l_optsGl = new unsigned long[l_nOpts];
      // compute sum over all ranks
      MPI_Allreduce( l_optsLo, l_optsGl, l_nOpts, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD );
#else
      l_optsGl = l_optsLo;
#endif

      for( unsigned short l_opt = 0; l_opt < l_nOpts; l_opt++ ) {
        EDGE_LOG_INFO << "    id " << l_opt / (C_ENT[l_elType].N_FACES*l_veJump)
                      << "-" << ( l_opt%(C_ENT[l_elType].N_FACES*l_veJump) ) / l_veJump
                      << "-" << l_opt%l_veJump
                      << ": " << l_optsGl[l_opt];
      }

      delete[] l_optsLo;
#ifdef PP_USE_MPI
      delete[] l_optsGl;
#endif
    }
};
#endif
