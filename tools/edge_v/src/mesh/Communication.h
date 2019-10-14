/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2019, Alexander Breuer
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
 * Derives communication structures.
 **/
#ifndef EDGE_V_MESH_COMMUNICATION_H
#define EDGE_V_MESH_COMMUNICATION_H

#include "../constants.h"
#include <cstddef>
#include <set>
#include <vector>

namespace edge_v {
  namespace mesh {
    class Communication;
  }
}

/**
 * Communication structures.
 **/
class edge_v::mesh::Communication {
  private:
    //! message of a single time group of a single partition
    struct Message {
      //! remote partition
      std::size_t pa;
      //! remote time group
      unsigned short tg;
      //! communicating elements of the time group
      std::vector< std::size_t > el;
      //! communicating element-faces of the time group
      std::vector< unsigned short > fa;
    };

    //! communication structure of a time region of a single partition
    struct TimeRegion {
      //! time group
      unsigned short tg;
      //! outgoing messages
      std::vector< Message > send;
      //! incoming messages
      std::vector< Message > recv;
    };

    //! communication structure of a single partition, composed of communicating time regions
    struct Partition {
      //! time regions of the partition
      std::vector< TimeRegion > tr;
    };

    //! global communication structure, composed of communicating partitions
    std::vector< Partition > m_struct;

    //! offsets for the channels of the partitions
    std::size_t * m_chOff = nullptr;

    //! communication channels of the partitions
    std::size_t * m_chs = nullptr;

    /**
     * Determines if an element is communicating.
     *
     * @param i_elTy type of the element.
     * @param i_el id of the element.
     * @param i_elFaEl elements adjacent to elements (faces as bridge).
     * @param i_elPa partitions of the elements.
     * @return true if comm, false if not.
     **/
    static bool isComm( t_entityType         i_elTy,
                        std::size_t          i_el,
                        std::size_t  const * i_elFaEl,
                        std::size_t  const * i_elPa );

    /**
     * Gets the partitions' elements which are communicating with elements of other partitions.
     *
     * @param i_elTy element type.
     * @param i_nEls number of elements.
     * @param i_elFaEl elements adjacent to elements, faces as bridge.
     * @param i_elPa partitions of the elements.
     * @param o_first will be set to first elements of the comm regions.
     * @param o_size will be set to the number of elements of the comm regions.
     **/
    static void getPaElComm( t_entityType        i_elTy,
                             std::size_t         i_nEls,
                             std::size_t const * i_elFaEl,
                             std::size_t const * i_elPa,
                             std::size_t       * o_first,
                             std::size_t       * o_size );
    /**
     * Gets the partition-time group pairs in the given comm-region.
     *
     * @param i_elTy type of the elements.
     * @param i_first first element of the region.
     * @param i_size size of the region.
     * @param i_elFaEl elements adjacent to elements (faces as bridge).
     * @param i_elTg time groups of the elements.
     * @param i_elPa partitions of the elements.
     * @param o_pairs will be set to the partition-time group pairs (partition, time group).
     **/
    static void getPaTgPairs( t_entityType                           i_elTy,
                              std::size_t                            i_first,
                              std::size_t                            i_size,
                              std::size_t                    const * i_elFaEl,
                              unsigned short                 const * i_elTg,
                              std::size_t                    const * i_elPa,
                              std::set<
                               std::pair< std::size_t,
                                          unsigned short > >       & o_pairs  );

    /**
     * Gets the send messages for a single communication region.
     *
     * @param i_elTy element type.
     * @param i_first first element of the communication region.
     * @param i_size number of elements in the communication region.
     * @param i_elFaEl elements adjacent to elements, faces as bridge.
     * @param i_elPa partitions of the elements.
     * @param i_elTg time groups of the elements.
     * @param o_send will be set to send messags of the communication region.
     **/
    static void getMsgsSend( t_entityType                   i_elTy,
                             std::size_t                    i_first,
                             std::size_t                    i_size,
                             std::size_t            const * i_elFaEl,
                             std::size_t            const * i_elPa,
                             unsigned short         const * i_elTg,
                             std::vector< Message >       & o_send );

    /**
     * Gets the global communication structure.
     * The input data is assumed to be sorted by partition and time group.
     *
     * @param i_elTy element type.
     * @param i_nEls number of elements.
     * @param i_elFaEl elements adjacent to elements (faces as bridge).
     * @param i_elPa partitions of the elements.
     * @param i_elTg time groups of the elements.
     * @param o_struct will be set to global communication structure.
     **/
    static void getStruct( t_entityType                   i_elTy,
                           std::size_t                    i_nEls,
                           std::size_t            const * i_elFaEl,
                           std::size_t            const * i_elPa,
                           unsigned short         const * i_elTg,
                           std::vector< Partition >     & o_struct );

    /**
     * Sets the channel offsets.
     *
     * @param i_struct global communication structure.
     * @param o_chOff will be set to the per-partition offsets of the channels.
     **/
    static void setChOff( std::vector< Partition > const & i_struct,
                          std::size_t                    * o_chOff );

    /**
     * Sets the structure for the communication channels.
     *
     * @param i_struct global communication structure.
     * @param i_chOff per-partition offsets of the channels.
     * @param o_chs will be set to the communication channels.
     **/
    static void setChs( std::vector< Partition > const & i_struct,
                        std::size_t              const * i_chOff,
                        std::size_t                    * o_chs );

  public:
    /**
     * Constructor which initializes the communication structures.
     *
     * @param i_elTy element type.
     * @param i_nEls number of elements.
     * @param i_elFaEl elements adjacent to elements (faces as bridge).
     * @param i_nPas number of partitions.
     * @param i_nPaEls number of elements per partition.
     * @param i_elTg time groups of the elements.
     **/
    Communication( t_entityType           i_elTy,
                   std::size_t            i_nEls,
                   std::size_t    const * i_elFaEl,
                   std::size_t            i_nPas,
                   std::size_t    const * i_nPaEls,
                   unsigned short const * i_elTg );

    /**
     * Destructor
     **/
    ~Communication();

    /**
     * Gets the communication structure for a given partition.
     * The first entry is the number of communication channels.
     * Next, for each channel, the following 4D structure follows:
     *   First, the local time group.
     *   Next, the remote partition.
     *   Third, the remote time group
     *   Last the number of faces to which this applies.
     *
     * For example [2, 0, 2, 4, 15, 1, 1, 3, 20] means:
     *   Two channels with structure: [0, 2, 4, 15] and [1, 1, 3, 20].
     *
     *   The first channel belongs to local time group 0 and communicates with remote time group 4 of partition 2.
     *   The number of faces to which this applies is 15.
     *
     *   The second channel belongs to local time group 1 and communicates with the third time group on partition 1.
     *   The number of faces is 20.
     *
     * @return number of messages for the partitions.
     **/
    std::size_t const * getStruct( std::size_t i_pa ) const;
};

#endif