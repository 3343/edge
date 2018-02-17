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
 * Writer for internal sub-face boundaries.
 **/
#ifndef EDGE_IO_INTERNAL_BOUNDARY_HPP
#define EDGE_IO_INTERNAL_BOUNDARY_HPP

#include "constants.hpp"
#include "io/logging.h"
#include "data/Dynamic.h"
#include "FileSystem.hpp"
#include "linalg/Mappings.hpp"
#include "submodules/include/visit_writer.h"

namespace edge {
  namespace io {
    template< typename       TL_T_LID,
              t_entityType   TL_T_EL,
              unsigned short TL_O_SP >
    class InternalBoundary;
  }
}

/**
 * Internal boundaries.
 *
 * @paramt TL_T_LID integral type for local ids.
 * @paramt TL_T_EL element type.
 * @paramt TL_O_SP spatial order of the DG-scheme.
 **/
template< typename       TL_T_LID,
          t_entityType   TL_T_EL,
          unsigned short TL_O_SP >
class edge::io::InternalBoundary {
  private:
    //! number of dims
    static unsigned short const TL_N_DIS = C_ENT[TL_T_EL].N_DIM;

    //! number of vertices per face
    static unsigned short const TL_N_FA_VES = C_ENT[TL_T_EL].N_FACE_VERTICES;

    //! number of vertices per element
    static unsigned short const TL_N_EL_VES = C_ENT[TL_T_EL].N_VERTICES;

    //! number of faces
    static unsigned short const TL_N_FAS = C_ENT[TL_T_EL].N_FACES;

    //! number of sub-vertices per element face
    static unsigned short const TL_N_FA_SVS = CE_N_SUB_VERTICES( C_ENT[TL_T_EL].TYPE_FACES, TL_O_SP );

    //! number of sub-faces per element face
    static unsigned short const TL_N_SFS = CE_N_SUB_FACES( TL_T_EL, TL_O_SP );

    //! number of sub-cells per element
    static unsigned short const TL_N_SCS  = CE_N_SUB_CELLS( TL_T_EL, TL_O_SP );

    //! vertices of the sub-faces
    float (* m_svCrds)[TL_N_FA_SVS][3];

    //! mapping from elements to sub-vertices
    TL_T_LID (* m_bfSfSv)[TL_N_SFS][TL_N_FA_VES];

    //! data buffer
    float *m_buffer;

    //! pointers to buffered data
    float **m_bPtrs;

    //! output base path
    std::string const m_outPath;

    //! number of steps written to disk
    unsigned int m_writeStep = 0;

    //! binary or ascii output
    bool const m_binary;

    //! visit element type
    int m_visitElType;

    /**
     * Derives the face-loca "sub-grid" in reference coordinats.
     *
     * @param i_scSv sub-vertices adjacent to the sub-cells (no bridge).
     * @param o_faSfSvL will be set to sub-vertices adjacent to sub-faces of a face as local sub-grid ids.
     * @param o_faSvR will be set to sub-vertex ids of the new "sub-grid" in terms of the original volume sub-grid.
     **/
    void faSg( unsigned short const i_scSv[ TL_N_SCS + TL_N_FAS * TL_N_SFS ][ TL_N_EL_VES ],
               unsigned short       o_faSfSvL[TL_N_FAS][TL_N_SFS][TL_N_FA_VES],
               unsigned short       o_faSvR[TL_N_FAS][TL_N_FA_SVS] ) {
      for( unsigned short l_fa = 0; l_fa < TL_N_FAS; l_fa++ ) {
        // init with invalid values
        for( unsigned short l_sv = 0; l_sv < TL_N_FA_SVS; l_sv++ )
          o_faSvR[l_fa][l_sv] = std::numeric_limits< unsigned short >::max();

          // iterate over sub-faces
          for( unsigned short l_sf = 0; l_sf < TL_N_SFS; l_sf++ ) {
          // corresponding receive sub-cell
          unsigned short l_scRecv = TL_N_SCS + l_fa * TL_N_SFS + l_sf;

          // assign the sub-vertex ids
          for( unsigned short l_sv = 0; l_sv < TL_N_FA_VES; l_sv++ ) {
            for( unsigned short l_cu = 0; l_cu < TL_N_FA_SVS; l_cu++ ) {
              // assign if we reached the current maximum
              if( o_faSvR[l_fa][l_cu] == i_scSv[l_scRecv][l_sv] ) {
                // store and break
                o_faSfSvL[l_fa][l_sf][l_sv] = l_cu;
                break;
              }
              else if( o_faSvR[l_fa][l_cu] == std::numeric_limits< unsigned short >::max() ) {
                // store and break
                o_faSvR[l_fa][l_cu] = i_scSv[l_scRecv][l_sv];
                o_faSfSvL[l_fa][l_sf][l_sv] = l_cu;
                break;
              }
            }
          }
        }
      }
    }

    /**
     * Copies the given data to the output data.
     *
     * @param i_first first boundary face.
     * @param i_size number of boundary faces.
     * @param i_nQts number of quantities, which are written.
     * @param i_stride stride from one sub-face to the next.
     * @param i_data data which is copied.
     **/
    template< typename TL_T_REAL >
    void copy( TL_T_LID        i_first,
               TL_T_LID        i_size,
               unsigned short  i_nQts,
               unsigned short  i_stride,
               TL_T_REAL      *i_data ) {
      // check for valid input
      EDGE_CHECK_GE( i_first,  0 );
      EDGE_CHECK_GE( i_size,   0 );
      EDGE_CHECK_GT( i_nQts,   0 );
      EDGE_CHECK_GE( i_stride, 1 );

      // iterate over boundary faces
#ifdef PP_USE_OMP
#pragma omp parallel for
#endif
      for( TL_T_LID l_bf = i_first; l_bf < i_first+i_size; l_bf++ ) {
        for( unsigned short l_sf = 0; l_sf < TL_N_SFS; l_sf++ ) {
          // copy the data to the buffer
          for( unsigned short l_qt = 0; l_qt < i_nQts; l_qt++ ) {
            // buffer id
            std::size_t l_bId  = l_qt * std::size_t(i_size) * TL_N_SFS;    // jump over quantities
                        l_bId += (l_bf - i_first) * std::size_t(TL_N_SFS); // jump over boundary faces
                        l_bId += l_sf;                                     // jump over sub-faces

            // data id
            std::size_t l_dId  = l_bf * std::size_t(i_stride) * TL_N_SFS; // jump over boundary faces
                        l_dId += l_sf * i_stride;                         // jump over sub-faces
                        l_dId += l_qt;                                    // jump over quantities

            m_buffer[l_bId] = i_data[l_dId];
          }
        }
      }

      // set pointers to stride-1 data
      for( unsigned short l_qt = 0; l_qt < i_nQts; l_qt++ )
        m_bPtrs[l_qt] = m_buffer + l_qt * std::size_t(i_size) * TL_N_SFS;
    }

  public:
    /**
     * Constructor.
     *
     * @param i_outPath base path to where the output is written.
     **/
    InternalBoundary( std::string &i_outPath,
                      bool i_binary=1 ): m_outPath(i_outPath), m_binary(i_binary) {
      // init directory if given
      FileSystem::createDir( m_outPath );

      // set plotting type
      if(       C_ENT[TL_T_EL].TYPE_FACES == LINE   ) m_visitElType = VISIT_LINE;
      else if ( C_ENT[TL_T_EL].TYPE_FACES == TRIA3  ) m_visitElType = VISIT_TRIANGLE;
      else if ( C_ENT[TL_T_EL].TYPE_FACES == QUAD4R ) m_visitElType = VISIT_QUAD;
      else if ( C_ENT[TL_T_EL].TYPE_FACES == HEX8R  ) m_visitElType = VISIT_HEXAHEDRON;
      else if ( C_ENT[TL_T_EL].TYPE_FACES == TET4   ) m_visitElType = VISIT_TETRA;
      else EDGE_LOG_FATAL << "missing element type " << C_ENT[TL_T_EL].TYPE_FACES;
    }

    /**
     * Allocates memory for the internal boundary writer.
     *   Returns silently if the number of boundary faces or floats is zero.
     *
     * @param i_nBf number of boundary faces.
     * @param i_nMaxQts maximum number of quantities, which the boundary writer can hold for each sub-face.
     * @param io_dynMem dynamic memory allocator.
     **/
    void alloc( TL_T_LID             i_nBf,
                unsigned short       i_nMaxQts,
                edge::data::Dynamic &io_dynMem ) {
      if( i_nBf     == 0 ) return;
      if( i_nMaxQts == 0 ) return;

      // vertices (stores vertices at element edges twice)
      std::size_t l_veCrdsSize  = i_nBf * std::size_t(TL_N_FA_SVS);
                  l_veCrdsSize *= std::size_t(3) * sizeof( float );
      m_svCrds = (float (*)[TL_N_FA_SVS][3]) io_dynMem.allocate( l_veCrdsSize );

      // mapping
      std::size_t l_faSfSvSize  = i_nBf * std::size_t(TL_N_SFS);
                  l_faSfSvSize *= std::size_t( TL_N_FA_VES ) * sizeof( TL_T_LID );
      m_bfSfSv = (TL_T_LID (*)[TL_N_SFS][TL_N_FA_VES]) io_dynMem.allocate( l_faSfSvSize );

      // data buffer
      std::size_t l_bufferSize  = i_nBf * std::size_t(TL_N_SFS);
                  l_bufferSize *= std::size_t(i_nMaxQts) * sizeof( TL_T_LID );
      m_buffer = (float*) io_dynMem.allocate( l_bufferSize );

      // pointers to buffer
      std::size_t l_bPtrs  = i_nMaxQts;
                  l_bPtrs *= sizeof(float*);
      m_bPtrs = (float**) io_dynMem.allocate( l_bPtrs );
    }

    /**
     * Initializes the internal boundary plotter.
     *
     * @param i_nBf number of boundary faces.
     * @param i_scSv sub-vertices adjacent to subcells (no bridge).
     * @param i_bfBe boundary elements adjacent to boundary faces (no bridge).
     * @param i_beEl dense elements ids of boundary elements.
     * @param i_elVe vertices adjacent to elements (no bridge)
     * @param i_charsSv characteristics of the sub-vertices.
     * @param i_charsVe characteristics of the vertices.
     * @param i_charsBf characteristics of the internal boundary faces.
     *
     * @paramt TL_T_CHARS_SV sub-vertex characteristics with member .coords.
     * @paramt TL_T_CHARS_VE vertex characteristics with member .coords.
     * @paramt TL_T_CHARS_BF internal boundary characteristics with member .fIdBfEl.
     **/
    template< typename TL_T_CHARS_SV,
              typename TL_T_CHARS_VE,
              typename TL_T_CHARS_BF >
    void init( TL_T_LID                    i_nBf,
               unsigned short      const   i_scSv[ TL_N_SCS + TL_N_FAS * TL_N_SFS ][ TL_N_EL_VES ],
               TL_T_LID            const (*i_bfBe)[2],
               TL_T_LID            const  *i_beEl,
               TL_T_LID            const (*i_elVe)[TL_N_EL_VES],
               TL_T_CHARS_SV       const  *i_charsSv,
               TL_T_CHARS_VE       const  *i_charsVe,
               TL_T_CHARS_BF       const  *i_charsBf ) {
      // get face "subgrid"
      unsigned short l_faSfSvL[TL_N_FAS][TL_N_SFS][TL_N_FA_VES];
      unsigned short l_faSvR[TL_N_FAS][TL_N_FA_SVS];
      faSg( i_scSv,
            l_faSfSvL,
            l_faSvR );

      // iterate over boundary faces
#ifdef PP_USE_OMP
#pragma omp parallel for
#endif
      for( TL_T_LID l_bf = 0; l_bf < i_nBf; l_bf++ ) {
        // left element's dense id
        TL_T_LID l_el = i_beEl[ i_bfBe[l_bf][0] ];

        // local face id of left element
        unsigned short l_fId = i_charsBf[l_bf].fIdBfEl[0];

        // get coordinates of left element's vertices
        float l_veCrds[TL_N_DIS][TL_N_EL_VES];

        for( unsigned short l_ve = 0; l_ve < TL_N_EL_VES; l_ve++ ) {
          TL_T_LID l_veId = i_elVe[l_el][l_ve];

          for( unsigned short l_di = 0; l_di < TL_N_DIS; l_di++ )
            l_veCrds[l_di][l_ve] = i_charsVe[l_veId].coords[l_di];
        }

        // get coordinates of face "subgrid"
        float l_sgCrds[TL_N_FA_SVS][TL_N_DIS];

        for( unsigned short l_sv = 0; l_sv < TL_N_FA_SVS; l_sv++ )
          for( unsigned short l_di = 0; l_di < TL_N_DIS; l_di++ )
            l_sgCrds[l_sv][l_di] = i_charsSv[ l_faSvR[l_fId][l_sv] ].coords[l_di];

        // determine physical coords of the subgrid
        for( unsigned short l_sv = 0; l_sv < TL_N_FA_SVS; l_sv++ ) {
          // init higher dims with zero
          for( unsigned short l_di = 0; l_di < 3; l_di++ )
            m_svCrds[l_bf][l_sv][l_di] = 0;

          // compute the mapping
          edge::linalg::Mappings::refToPhy( TL_T_EL,
                                            l_veCrds[0],
                                            l_sgCrds[l_sv],
                                            m_svCrds[l_bf][l_sv] );
        }
        // save vertex ids for this face subgrid
        for( unsigned short l_sf = 0; l_sf < TL_N_SFS; l_sf++ ) {
          for( unsigned short l_ve = 0; l_ve < TL_N_FA_VES; l_ve++ ) {
            m_bfSfSv[l_bf][l_sf][l_ve]  = l_bf * TL_N_FA_SVS;
            m_bfSfSv[l_bf][l_sf][l_ve] += l_faSfSvL[l_fId][l_sf][l_ve];
          }
        }
      }
    }

    /**
     * Writes the given internal boundary data to disk.
     *
     * @param i_first first boundary face.
     * @param i_size number of boundary faces.
     * @param i_nQts number of quantities, which are written.
     * @param i_stride stride from one sub-face to the next.
     * @param i_namesQts names of the quantities.
     * @param i_data data which is written to disk.
     **/
    template< typename TL_T_REAL >
    void write( TL_T_LID         i_first,
                TL_T_LID         i_size,
                unsigned short   i_nQts,
                unsigned short   i_stride,
                char const     **i_namesQts,
                TL_T_REAL*       i_data ) {
      // abort if nothing to write
      if( i_size == 0 ) return;

      // check for a positive stride
      EDGE_CHECK_NE( i_stride, 0 );
      EDGE_CHECK_GE( i_stride, i_nQts );

      // copy data to buffer
      copy( i_first,
            i_size,
            i_nQts,
            i_stride,
            i_data );

      // assemble file name
      std::string l_file = m_outPath;
      l_file += "_" + parallel::g_rankStr;
      l_file += "_" + std::to_string((unsigned long long) m_writeStep) + ".vtk";

      int l_nPts = i_size * TL_N_FA_SVS;
      int l_nCells = i_size * TL_N_SFS;

      edge_write_unstructured_mesh( l_file.c_str(),
                                    m_binary,
                                    l_nPts,
                                    m_svCrds[i_first][0],
                                    i_nQts,
                                    l_nCells,
                                    m_visitElType,
                                    m_bfSfSv[i_first][0],
                                    i_namesQts,
                                    m_bPtrs );

      // increase write counter
      m_writeStep++;
    }
};

#endif
