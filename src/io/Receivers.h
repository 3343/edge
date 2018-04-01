/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2016, Regents of the University of California
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
 * Output of receivers.
 **/

#ifndef EDGE_IO_RECEIVERS_H
#define EDGE_IO_RECEIVERS_H

#include "constants.hpp"
#include "data/EntityLayout.type"
#include <string>

namespace edge {
  namespace io {
    class Receivers;
  }
}

class edge::io::Receivers {
  protected:
    struct Recv {
      //! number of buffered values
      unsigned int nBuff;
      //! buffer
      std::vector< real_base > buffer;
      //! buffered times
      std::vector< real_base > buffTime;
      //! time of the receiver
      double time;
      //! time group of the receiver
      int_tg tg;
      //! entity of the receiver
      int_el en;
      //! entity of the receiver w.r.t. to the start of the time group
      int_el enTg;
      //! evaluated basis at the receiver's location
      real_base evaBasis[N_ELEMENT_MODES];
      //! path to the receiver's file
      std::string path;
    };
    // receiver under control; entity ids are ascending
    std::vector< Recv > m_recvs;

    //! mapping from sparse entities to the (first) receiver
    std::vector< std::size_t > m_spEnToRecv;

    //! number of quantities
    unsigned short m_nQts = 0;

    //! size of the buffer
    unsigned int m_buffSize;

    //! sampling frequency of the receivers
    double m_freq;

    /**
     * Touches the output for the first time and writes the headers.
     *
     * @param i_outDir output directory which gets created if it does not exist.
     **/
    void touchOutput( const std::string &i_outDir );

    /**
     * Flushes a receiver to disk.
     *
     * @param i_recv receiver which gets flushed.
     **/
    void flush( unsigned int i_recv );

    /**
     * Flushes all receivers to disk.
     **/
    void flushAll();
  public:
    /**
     * Destructor which flushes everything to disk.
     **/
    ~Receivers() { flushAll(); };

    /**
     * Prints statistics of the receivers.
     **/
    void print();

    /**
     * Initialzes the receiver output.
     *   TODO: The current implementation is limited to elements.
     *
     * @param i_enType entity type in which receivers are searched for.
     * @param i_nRecvs number of receivers.
     * @param i_recvNames name of the receivers.
     * @param i_recvCrds coordinates of the receivers.
     * @param i_freq sampling frequency of the receivers.
     * @param i_enLayout layout of the entities in which receivers are search for.
     * @param i_enVe vertices adjacent to the entities.
     * @param i_veChars vertex chars.
     * @param i_bufferSize size of the internal receiver buffer before data gets written to disk.
     * @param i_time time of the first receiver outout
     **/
    void init(       t_entityType     i_enType,
                     unsigned int     i_nRecvs,
               const std::string     &i_outDir,
               const std::string    (*i_recvNames),
               const real_mesh      (*i_recvCrds)[3],
                     double           i_freq,
               const t_enLayout      &i_enLayout,
               const int_el         (*i_enVe),
               const t_vertexChars   *i_veChars,
                     unsigned int     i_bufferSize=100,
                     double           i_time=0 );

    /**
     * Gets the dense-entities with receivers in them.
     * The position in the array is equivalent to the id of the element's associated receiver(s).
     *
     * @param o_en will be set to entities with receivers.
     **/
    void getEnRecv( std::vector< int_el > & o_en );

    /**
     * Gets the relative time (w.r.t. i_time) at which the receiver(s) of the given sparse entity expects the next output.
     *
     * @param i_spEn sparse-entity.
     * @param i_time considered start time.
     * @param i_dt considered delta t, speciying the lengtrh.
     *
     * @return relative dt w.r.t. i_time; if no output is expected a negative value is returned.
     **/
    double getRecvTimeRel( int_el i_spEn,
                           double i_time,
                           double i_dt );

    /**
     * Writes the receiver(s) for the give sparse entity.
     * A single sparse entity can point to multiple receivers.
     *
     * Remark: It is the callers responsibility to ensure that the DOFs are evaluated at the correct time.
     *
     * @param i_spEn sparse entity.
     * @param i_dofs degrees of freedom, which get evaluated in space.
     **/
    void writeRecvAll(        int_el   i_spEn,
                       const real_base i_dofs[N_QUANTITIES][N_ELEMENT_MODES][N_CRUNS] );

    /**
     * Flushes receiver's buffers to disk if the remaining size in the buffer if below the treshold.
     *
     * @param i_treshold, which triggers writers if the buffers remaining entries is below. 
     **/
    virtual void flushIf( unsigned int i_tresh=50 );
};

#endif
