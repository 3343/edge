/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2018, Regents of the University of California
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
 * Interface to MOAB.
 **/
#ifndef EDGEV_IO_MOAB_H
#define EDGEV_IO_MOAB_H

#include <moab/Interface.hpp>
#include <moab/Core.hpp>

#include <string>

namespace edge_v {
  namespace io {
    template< typename TL_T_GID >
    class Moab;
  }
}

/**
 * @brief Interface to MOAB>
 * 
 * @tparam TL_T_GID integral type of global ids (only changes returned values, but not MOAB's internal storage).
 */
template< typename TL_T_GID >
class edge_v::io::Moab {
  private:
    //! moab interface
    moab::Interface *m_moab;

  public:
    /**
     * @brief Constructs a new Moab object.
     * 
     * @param i_pathToMesh path to the mesh.
     */
    Moab( std::string i_pathToMesh ) {
      m_moab = new moab::Core;

      moab::ErrorCode l_err = m_moab->load_file( i_pathToMesh.c_str() );
      assert( l_err == moab::MB_SUCCESS );
    }

    /**
     * @brief Destroys the Moab object.
     */
    ~Moab() {
      delete m_moab;
    }

    /**
     * @brief Gets the number of entities by the number dimensions.
     * 
     * @param i_nDis number of dimensions.
     *
     * @return number of entities. 
     */
    TL_T_GID nEnsByDis( unsigned short i_nDis ) {
      TL_T_GID l_nEns = std::numeric_limits< TL_T_GID >::max();
      moab::ErrorCode l_err = m_moab->get_number_entities_by_dimension( 0,
                                                                        i_nDis,
                                                                        l_nEns );
      assert( l_err == moab::MB_SUCCESS );

      return l_nEns;
    }

    /**
     * @brief Gets the coordinates of the (ordered) vertices.
     *
     * @param o_veCrds will be set to coordinates of the vertices (SoA).
     */
    void getVeCrds( double(*o_veCrds)[3] ) {
      // get the vertices in the mesh
      std::vector< moab::EntityHandle > l_ves;
      moab::ErrorCode l_err = m_moab->get_entities_by_dimension( 0,
                                                                 0,
                                                                 l_ves );
      assert( l_err == moab::MB_SUCCESS );

      l_err = m_moab->get_coords( &l_ves[0],
                                   l_ves.size(),
                                   o_veCrds[0] );
      assert( l_err == moab::MB_SUCCESS );
    }

    /**
     * @brief Gets the connectivity information for the given entity type.
     * 
     * @param i_enTy entity type, either tet4 or tria3.
     * @param o_enVe will be set vertex ids (starting at 0) adjacent to the entities.
     */
    void getEnVe( std::string  i_enTy,
                  TL_T_GID    *o_enVe ) {
      moab::EntityType l_ty;
      if( i_enTy == "tria3" ) {
        l_ty = moab::MBTRI;
      }
      else if( i_enTy  == "tet4" ) {
        l_ty = moab::MBTET;
      }

      // get mapping from entities to vertices
      std::vector< moab::EntityHandle > l_enCo;
      moab::ErrorCode l_err = m_moab->get_connectivity_by_type( l_ty,
                                                                l_enCo );
      assert( l_err == moab::MB_SUCCESS );

      // translate to ids
      for( std::size_t l_ve = 0; l_ve < l_enCo.size(); l_ve++ ) {
        o_enVe[l_ve] = m_moab->id_from_handle( l_enCo[l_ve] ) - 1;
      }
    }
};

#endif