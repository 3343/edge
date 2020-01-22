/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2019-2020, Alexander Breuer
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
#ifndef EDGE_V_IO_MOAB_H
#define EDGE_V_IO_MOAB_H

#include <moab/Interface.hpp>
#include <moab/Core.hpp>
#include <string>
#include "../constants.h"
#include "logging.h"

namespace edge_v {
  namespace io {
    class Moab;
  }
}

/**
 * MOAB interface.
 **/
class edge_v::io::Moab {
  private:
    //! moab interface
    moab::Interface *m_moab;

    //! root set
    moab::EntityHandle m_root;

    //! number points
    std::size_t m_nPoint;

    //! number of line elements
    std::size_t m_nLine;

    //! number quad4r elements
    std::size_t m_nQuad4r;

    //! number of tria3 elements
    std::size_t m_nTria3;

    //! number of hex8r elements
    std::size_t m_nHex8r;

    //! number of tet4 elements
    std::size_t m_nTet4;

    /**
     * Converts the given type to a native MOAB type.
     *
     * @param i_enTy EDGE-V entity type.
     * @return corresponding MOAB type.
     **/
    static moab::EntityType getMoabType( t_entityType i_enTy );

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
     * Sets the given data in MOAB.
     *
     * @param i_enTy entity type to which this data belongs.
     * @param i_daTy used MOAB data type.
     * @param i_tagName tag name.
     * @param i_data data, which will be stored.
     **/
    void setEnData( t_entityType        i_enTy,
                    moab::DataType      i_daTy,
                    std::string const & i_tagName,
                    void        const * i_data );

    /**
     * Gets the intersection of the entities and the mesh sets.
     * Returned are newly created mesh sets for every intesection with the input mesh sets.
     *
     * @param i_ens entities which are intersected.
     * @param i_ms input mesh sets.
     * @param o_ms output mesh sets.
     **/
    void getEnMsInter( std::vector< moab::EntityHandle > const & i_ens,
                       std::vector< moab::EntityHandle > const & i_ms,
                       std::vector< moab::EntityHandle >       & o_ms );

    /**
     * Sets the tags of the output entities, based on the input entities.
     *
     * @param i_tags tags which are copied.
     * @param i_ens input entities from which the tags are read.
     * @param o_ens output entities for which the tags are set.
     **/
    void setTags( std::vector< moab::Tag >          const & i_tags,
                  std::vector< moab::EntityHandle > const & i_ens,
                  std::vector< moab::EntityHandle >       & io_ens );

  public:
    /**
     * Constructs a new Moab object.
     * 
     * @param i_pathToMesh path to the mesh.
     **/
    Moab( std::string const & i_pathToMesh );

    /**
     * Destroys the Moab object.
     **/
    ~Moab();

    /**
     * Prints mesh statistics.
     **/
    void printStats() const;

    /**
     * Gets the entity type of the mesh's elements.
     *
     * @return element type of the mesh.
     **/
    t_entityType getElType() const;

    /**
     * Gets the number of entities by the number dimensions.
     * 
     * @param i_nDis number of dimensions.
     * @return number of entities.
     **/
    std::size_t nEnsByDis( unsigned short i_nDis ) const;

    /**
     * Gets the number of entities by entity type.
     * 
     * @param i_enTy entity type.
     * @return number of entities.
     **/
    std::size_t nEnsByType( t_entityType i_enTy ) const;

    /**
     * Gets the coordinates of the (ordered) vertices.
     *
     * @param o_veCrds will be set to coordinates of the vertices.
     **/
    void getVeCrds( double(*o_veCrds)[3] ) const;

    /**
     * Gets the connectivity information for the given entity type.
     * The vertices are returned according to the entities canonical ordering.
     *   See: Tim Tautges, Canonical Numbering Systems for Finite-Element Codes.  
     * 
     * @param i_enTy entity type.
     * @param o_enVe will be set vertex ids (starting at 0) adjacent to the entities.
     **/
    void getEnVe( t_entityType   i_enTy,
                  std::size_t  * o_enVe ) const;

    /**
     * Gets the elements adjacent to the faces.
     *
     * @param i_elTy element type to which the faces belong.
     * @param o_faEl will be set to the elements, adjacent to the faces.
     **/
    void getFaEl( t_entityType   i_elTy,
                  std::size_t  * o_faEl ) const;

    /**
     * Gets the faces adjacent to the elements.
     *
     * @param i_elTy element type.
     * @param o_elFa will be set to adjacency information.
     **/
    void getElFa( t_entityType   i_elTy,
                  std::size_t  * o_elFa ) const;

    /**
     * Gets the element-to-element adjacency (faces as bridge).
     *
     * @param i_elTy element type.
     * @param o_elFaEl will be set to adjacency information.
     **/
    void getElFaEl( t_entityType   i_elTy,
                    std::size_t  * o_elFaEl ) const;

    /**
     * Sets the given data in MOAB (stored as native double).
     * The data is stored in MOAB's native double type for portability.
     *
     * @param i_enTy entity type to which this data belongs.
     * @param i_tagName tag name.
     * @param i_data data, which will be stored.
     **/
    void setEnData( t_entityType        i_enTy,
                    std::string const & i_tagName,
                    float       const * i_data );

    /**
     * Sets the given data in MOAB (stored as native double).
     * The data is stored in MOAB's native double type for portability.
     *
     * @param i_enTy entity type to which this data belongs.
     * @param i_tagName tag name.
     * @param i_data data, which will be stored.
     **/
    void setEnData( t_entityType        i_enTy,
                    std::string const & i_tagName,
                    double      const * i_data );

    /**
     * Sets the given data in MOAB.
     * The data is stored in MOAB's native int for portability.
     *
     * @param i_enTy entity type to which this data belongs.
     * @param i_tagName tag name.
     * @param i_data data, which will be stored.
     **/
    void setEnData( t_entityType           i_enTy,
                    std::string    const & i_tagName,
                    unsigned short const * i_data );

    /**
     * Sets the given data in MOAB.
     * The data is stored in MOAB's native int for portability.
     *
     * @param i_enTy entity type to which this data belongs.
     * @param i_tagName tag name.
     * @param i_data data, which will be stored.
     **/
    void setEnData( t_entityType         i_enTy,
                    std::string  const & i_tagName,
                    std::size_t  const * i_data );

    /**
     * Gets the data for the given entity type (stored as native int).
     *
     * @param i_enTy entity type to which this data belongs.
     * @param i_tagName tag name.
     * @param o_data output array.
     **/
    void getEnData( t_entityType           i_enTy,
                    std::string    const & i_tagName,
                    unsigned short       * o_data ) const;

    /**
     * Gets the data for the given entity type from the containing set.
     * If no containing meshset is found, std::numeric_limits< int >::max() is set for the respective entity.
     *
     * @param i_enTy entity type to which this data belongs.
     * @param i_tagName tag name.
     * @param o_data output array.
     **/
    void getEnDataFromSet( t_entityType           i_enTy,
                           std::string    const & i_tagName,
                           int                  * o_data ) const;

    /**
     * Gets the names of the tags.
     *
     * @param o_tagNames will be set to the names of the tags.
     **/
    void getTagNames( std::vector< std::string > & o_tagNames ) const;

    /**
     * Deletes the given tag, if it exists.
     *
     * @param i_tagName name of the tag which will be deleted.
     * @return true if the tag was delete, false if not (didn't exist).
     **/
    bool deleteTag( std::string const & i_tagName );

    /**
     * Reorders the entities by the given int-tag (ascending).
     *
     * @param i_enTy type of the entities, which are reordered.
     * @param i_tagName name of the tag.
     * @param i_skipValue tag with this value are skipped in the reordering.
     **/
    void reorder( t_entityType        i_enTy,
                  std::string const & i_tagName,
                  int                 i_skipValue = std::numeric_limits< int >::max() );

    /**
     * Writes the MOAB database to the given file.
     *
     * @param i_pathToMesh path to the mesh.
     **/
    void writeMesh( std::string const & i_pathToMesh );

    /**
     * Writes the given elements of the MOAB database to the file.
     *
     * @param i_first first element to write.
     * @param i_nEls number of elements which are written.
     * @param i_pathToMesh path to the mesh.
     **/
    void writeMesh( std::size_t         i_first,
                    std::size_t         i_nEls,
                    std::string const & i_pathToMesh );
};

#endif