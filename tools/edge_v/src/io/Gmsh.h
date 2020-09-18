/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (breuer AT mytum.de)
 *
 * @section LICENSE
 * Copyright (c) 2020, Alexander Breuer
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
 * Gmsh interface.
 **/
#ifndef EDGE_V_IO_GMSH_H
#define EDGE_V_IO_GMSH_H

#include "../constants.h"
#include <string>
#include <vector>

namespace edge_v {
  namespace io {
    class Gmsh;
  }
}

/**
 * Gmsh interface.
 **/
class edge_v::io::Gmsh {
  private:
    //! coordinates of the vertices
    std::vector< double > m_veCrds;

    //! tags of the vertices
    std::vector< std::size_t > m_veTags;

    //! tags of the elements
    std::vector< std::size_t > m_elTags;

    //! tags of the elements' vertices
    std::vector< std::size_t > m_elVeTags;

    //! tags of the element-faces' vertices
    std::vector< std::size_t > m_elFaVeTags;

    //! physical groups of the faces
    std::vector< int > m_physicalGroupsFa;

    //! vertices of the faces belonging to the physical groups
    std::vector< std::vector< std::size_t > > m_faVeTagsPhysical;

    /**
     * Converts the given type to a native Gmsh type.
     *
     * @param i_enTy EDGE-V entity type.
     * @return corresponding Gmsh type.
     **/
    static int getGmshType( t_entityType i_enTy );

    /**
     * Converts the given Gmsh type to an entity type of EDGE-V.
     *
     * @param i_gmshType type as defined in Gmsh.
     * @return corresponding EDGE-V entity type.
     **/
    static t_entityType getEntityType( int i_gmshType );

    /**
     * Gets the id of the specified values in the given vector of sorted values.
     *
     * @param i_value value which is searched for.
     * @param i_sortedValues vector of sorted values.
     * @return id of the value in the vector.
     **/
    static t_idx getId( std::size_t                        i_value,
                        std::vector< std::size_t > const & i_sortedValues );

  public:
    /**
     * Constructor which initializes the Gmsh-interface.
     **/
    Gmsh();

    /**
     * Destructor which finalizes the Gmsh-interface.
     **/
    ~Gmsh();

    /**
     * Sets the value of the Gmsh-variable.
     *
     * @param i_name name of the Gmsh-variable.
     * @param i_value value which is set.
     **/
    void setNumber( std::string const & i_name,
                    double              i_value ) const;

    /**
     * Opens the given file, formats as supported by Gmsh.
     *
     * @param i_pathToFile path to the file which is opened.
     **/
    void open( std::string const & i_pathToFile );

    /**
     * Reads the mesh information and stores it in internal data structures.
     **/
    void readMesh();

    /**
     * Writes the given file, formats as supported by Gmsh.
     *
     * @param i_pathToFile path to the file which is written.
     **/
    void write( std::string const & i_pathToFile );

    /**
     * Gets the entity type of the mesh's elements.
     *
     * @return element type of the mesh.
     **/
    t_entityType getElType() const;

    /**
     * Gets the number of vertices in the mesh.
     *
     * @return number of vertices.
     **/
    t_idx nVes() const;

    /**
     * Gets the number of elements in the mesh.
     *
     * @return number of elements.
     **/
    t_idx nEls() const;

    /**
     * Gets the number of physical face-groups.
     *
     * @return number of physical face-groups.
     **/
    t_idx nPhysicalGroupsFa() const;

    /**
     * Gets the tags of the physical face groups.
     *
     * @return tags of the physical face groups.
     **/
    int const * getPhysicalGroupsFa() const;

    /**
     * Gets the number of faces in the given physical group.
     *
     * @param i_physicalGroupFa tag of the physical group.
     * @return number of faces in the given physical group.
     **/
    t_idx nFas( int i_physicalGroupFa ) const;

    /**
     * Gets the coordinates of the vertices.
     *
     * @return coordinates of the vertices.
     **/
    void getVeCrds( double (*o_veCrds)[3] ) const;

    /**
     * Gets the vertices of the faces with the given physical tag.
     *
     * @param i_physicalGroupFa physical tag which is defined on the faces.
     * @param o_faVe will be set to the vertices of the faces.
     **/
    void getFaVe( int     i_physicalGroupFa,
                  t_idx * o_faVe ) const;

    /**
     * Gets the ids of the vertices adjacent to the elements.
     *   Slow dimension: elements.
     *   Fast dimension: each element's vertices.
     *
     * @param o_elVe will be set to the ids of the elements' vertices.
     **/
    void getElVe( t_idx * o_elVe ) const;

    /**
     * Gets the ids of the vertices adjacent to the elements' faces.
     *   Slowest dimension: elements.
     *   2nd slowest dimension: faces.
     *   Fast dimension: vertices of each face.
     *
     * @param o_elFaVe will be set to the ids of the vertices of the elements' faces.
     **/
    void getElFaVe( t_idx * o_elFaVe ) const;

    /**
     * Reorders the elements based on the given priorities.
     *
     * @param i_priorities priorities of the elements (lower is higher, sorted first).
     **/
    void reorder( t_idx const * i_priorities );

    /**
     * Partitions the mesh.
     *   This assumed linear storage of the partitioned elements,
     *   i.e., that the elements of the 1st partition are stored first,
     *   those of the 2nd next.
     *
     * @param i_nPas number of partitions.
     * @param i_nPaEls number of elements in each of the partitions.
     **/
    void partition( t_idx         i_nPas,
                    t_idx const * i_nPaEls ) const;

    /**
     * Writes the given data for the elements.
     *
     * @param i_name name of the data/view.
     * @param i_elData data which is set.
     * @param i_pathToFile path to the file to where the data is written.
     **/
    void writeElData( std::string           const & i_name,
                      std::vector< double > const & i_elData,
                      std::string           const & i_pathToFile  ) const;

    /**
     * Writes the given data for the elements.
     *
     * @param i_name name of the data/view.
     * @param i_elData data which is set.
     * @param i_pathToFile path to the file to where the data is written.
     **/
    template< typename T >
    void writeElData( std::string const & i_name,
                      T           const * i_elData,
                      std::string const & i_pathToFile  ) const {
      // convert data to match gmsh's expected format
      std::vector< double > l_elData;
      l_elData.resize( m_elTags.size() );
#ifdef PP_USE_OMP
#pragma omp parallel for
#endif
      for( std::size_t l_el = 0; l_el < m_elTags.size(); l_el++ ) {
        l_elData[l_el] = i_elData[l_el];
      }

      writeElData( i_name,
                   l_elData,
                   i_pathToFile );
    }
};

#endif