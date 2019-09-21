/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *         Alexander Heinecke (alexander.heinecke AT intel.com)
 *
 * @section LICENSE
 * Copyright (c) 2019, Alexander Breuer
 * Copyright (c) 2016-2018, Regents of the University of California
 * Copyright (c) 2016, Intel Corporation
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
 * Internal data structures in EDGE.
 **/

#ifndef EDGE_DATA_INTERNAL_HPP
#define EDGE_DATA_INTERNAL_HPP

#include <functional>

#include "../internal.hpp"
#include "common.hpp"
#include "io/logging.h"
#include "parallel/Shared.h"

#if defined PP_T_KERNELS_VANILLA
#include "data/MmVanilla.hpp"
#elif defined PP_T_KERNELS_XSMM_DENSE_SINGLE
#include "data/MmXsmmSingle.hpp"
#else
#include "data/MmXsmmFused.hpp"
#endif

namespace edge {
  namespace data {
    class Internal;
  }
}

class edge::data::Internal {
  public:
    /**
     * Memory types.
     * TODO: Extend with more data structures.
     **/
    struct {
      struct {
#ifdef PP_ELEMENT_PRIVATE_1_HBW
        bool elementPrivate1 = true;
#else
        bool elementPrivate1 = false;
#endif
#ifdef PP_ELEMENT_PRIVATE_2_HBW
        bool elementPrivate2 = true;
#else
        bool elementPrivate2 = false;
#endif
#ifdef PP_ELEMENT_PRIVATE_3_HBW
        bool elementPrivate3 = true;
#else
        bool elementPrivate3 = false;
#endif

#ifdef PP_ELEMENT_MODE_PRIVATE_1_HBW
        bool elementModePrivate1 = true;
#else
        bool elementModePrivate1 = false;
#endif
#ifdef PP_ELEMENT_MODE_PRIVATE_2_HBW
        bool elementModePrivate2 = true;
#else
        bool elementModePrivate2 = false;
#endif
#ifdef PP_ELEMENT_MODE_PRIVATE_3_HBW
        bool elementModePrivate3 = true;
#else
        bool elementModePrivate3 = false;
#endif
#ifdef PP_SCRATCH_MEMORY_HBW
        bool scratchMem = true;
#else
        bool scratchMem = false;
#endif
#ifdef PP_ELEMENT_SHARED_1_HBW
        bool elementShared1 = true;
#else
        bool elementShared1 = false;
#endif
#ifdef PP_ELEMENT_SHARED_2_HBW
        bool elementShared2 = true;
#else
        bool elementShared2 = false;
#endif
#ifdef PP_ELEMENT_SHARED_3_HBW
        bool elementShared3 = true;
#else
        bool elementShared3 = false;
#endif
#ifdef PP_ELEMENT_SHARED_4_HBW
        bool elementShared4 = true;
#else
        bool elementShared4 = false;
#endif
      } hbw;
      struct {
#ifdef PP_ELEMENT_MODE_PRIVATE_1_HUGE
        bool elementModePrivate1 = true;
#else
        bool elementModePrivate1 = false;
#endif
#ifdef PP_ELEMENT_MODE_PRIVATE_2_HUGE
        bool elementModePrivate2 = true;
#else
        bool elementModePrivate2 = false;
#endif
#ifdef PP_ELEMENT_MODE_PRIVATE_3_HUGE
        bool elementModePrivate3 = true;
#else
        bool elementModePrivate3 = false;
#endif

#ifdef PP_ELEMENT_SHARED_1_HUGE
        bool elementShared1 = true;
#else
        bool elementShared1 = false;
#endif
#ifdef PP_ELEMENT_SHARED_2_HUGE
        bool elementShared2 = true;
#else
        bool elementShared2 = false;
#endif
#ifdef PP_ELEMENT_SHARED_3_HUGE
        bool elementShared3 = true;
#else
        bool elementShared3 = false;
#endif
#ifdef PP_ELEMENT_SHARED_4_HUGE
        bool elementShared4 = true;
#else
        bool elementShared4 = false;
#endif

#ifdef PP_SCRATCH_MEMORY_HUGE
        bool scratchMem = true;
#else
        bool scratchMem = false;
#endif
      } huge;
    } m_memTypes;

    /**
     * Meta data.
     **/
    bool m_initScratch;
    bool m_initDense;
    bool m_initSparse;

    //! number of vertices
    int_el m_nVertices;

    //! number of faces
    int_el m_nFaces;

    //! number of elements
    int_el m_nElements;

    //! number of sparse vertices #1
    int_el m_nVeSp1;
    //! number of sparse vertices #2
    int_el m_nVeSp2;
    //! number of sparse vertices #3
    int_el m_nVeSp3;
    //! number of sparse vertices #4
    int_el m_nVeSp4;
    //! number of sparse vertices #5
    int_el m_nVeSp5;
    //! number of sparse vertices #6
    int_el m_nVeSp6;

    //! number of sparse faces #1
    int_el m_nFaSp1;
    //! number of sparse faces #2
    int_el m_nFaSp2;
    //! number of sparse faces #3
    int_el m_nFaSp3;
    //! number of sparse faces #4
    int_el m_nFaSp4;
    //! number of sparse faces #5
    int_el m_nFaSp5;
    //! number of sparse faces #6
    int_el m_nFaSp6;

    //! number of sparse elements #1
    int_el m_nElSp1;
    //! number of sparse elements #2
    int_el m_nElSp2;
    //! number of sparse elements #3
    int_el m_nElSp3;
    //! number of sparse elements #4
    int_el m_nElSp4;
    //! number of sparse elements #5
    int_el m_nElSp5;
    //! number of sparse elements #6
    int_el m_nElSp6;

    // shared, private data
#ifdef PP_N_GLOBAL_PRIVATE_1
    t_globalPrivate1 m_globalPrivate1[PP_N_GLOBAL_PRIVATE_1][N_CRUNS] __attribute__ ((aligned (ALIGNMENT.BASE.STACK)));
#endif
#ifdef PP_N_GLOBAL_PRIVATE_2
    t_globalPrivate2 m_globalPrivate2[PP_N_GLOBAL_PRIVATE_2][N_CRUNS] __attribute__ ((aligned (ALIGNMENT.BASE.STACK)));
#endif
#ifdef PP_N_GLOBAL_PRIVATE_3
    t_globalPrivate3 m_globalPrivate3[PP_N_GLOBAL_PRIVATE_3][N_CRUNS] __attribute__ ((aligned (ALIGNMENT.BASE.STACK)));
#endif

    // shared, global data
#ifdef PP_N_GLOBAL_SHARED_1
    t_globalShared1 m_globalShared1[PP_N_GLOBAL_SHARED_1] __attribute__ ((aligned (ALIGNMENT.BASE.STACK)));
#endif
#ifdef PP_N_GLOBAL_SHARED_2
    t_globalShared2 m_globalShared2[PP_N_GLOBAL_SHARED_2] __attribute__ ((aligned (ALIGNMENT.BASE.STACK)));
#endif
#ifdef PP_N_GLOBAL_SHARED_3
    t_globalShared3 m_globalShared3[PP_N_GLOBAL_SHARED_3] __attribute__ ((aligned (ALIGNMENT.BASE.STACK)));
#endif
#ifdef PP_N_GLOBAL_SHARED_4
    t_globalShared4 m_globalShared4[PP_N_GLOBAL_SHARED_4] __attribute__ ((aligned (ALIGNMENT.BASE.STACK)));
#endif
#ifdef PP_N_GLOBAL_SHARED_5
    t_globalShared5 m_globalShared5[PP_N_GLOBAL_SHARED_5] __attribute__ ((aligned (ALIGNMENT.BASE.STACK)));
#endif
#ifdef PP_N_GLOBAL_SHARED_6
    t_globalShared6 m_globalShared6[PP_N_GLOBAL_SHARED_6] __attribute__ ((aligned (ALIGNMENT.BASE.STACK)));
#endif
#ifdef PP_N_GLOBAL_SHARED_7
    t_globalShared7 m_globalShared7[PP_N_GLOBAL_SHARED_7] __attribute__ ((aligned (ALIGNMENT.BASE.STACK)));
#endif
#ifdef PP_N_GLOBAL_SHARED_8
    t_globalShared8 m_globalShared8[PP_N_GLOBAL_SHARED_8] __attribute__ ((aligned (ALIGNMENT.BASE.STACK)));
#endif

    // sparse vertex data for every concurrent forward run
#ifdef PP_N_VERTEX_SPARSE_PRIVATE_1
    t_vertexSparsePrivate1 (*m_vertexSparsePrivate1)[PP_N_VERTEX_SPARSE_PRIVATE_1][N_CRUNS];
#endif
#ifdef PP_N_VERTEX_SPARSE_PRIVATE_2
    t_vertexSparsePrivate2 (*m_vertexSparsePrivate2)[PP_N_VERTEX_SPARSE_PRIVATE_2][N_CRUNS];
#endif
#ifdef PP_N_VERTEX_SPARSE_PRIVATE_3
    t_vertexSparsePrivate3 (*m_vertexSparsePrivate3)[PP_N_VERTEX_SPARSE_PRIVATE_3][N_CRUNS];
#endif
#ifdef PP_N_VERTEX_SPARSE_PRIVATE_4
    t_vertexSparsePrivate4 (*m_vertexSparsePrivate4)[PP_N_VERTEX_SPARSE_PRIVATE_4][N_CRUNS];
#endif
#ifdef PP_N_VERTEX_SPARSE_PRIVATE_5
    t_vertexSparsePrivate5 (*m_vertexSparsePrivate5)[PP_N_VERTEX_SPARSE_PRIVATE_5][N_CRUNS];
#endif
    // shared vertex data shared among the concurrent forward runs
#ifdef PP_N_VERTEX_SPARSE_SHARED_1
    t_vertexSparseShared1 (*m_vertexSparseShared1)[PP_N_VERTEX_SPARSE_SHARED_1];
#endif
#ifdef PP_N_VERTEX_SPARSE_SHARED_2
    t_vertexSparseShared2 (*m_vertexSparseShared2)[PP_N_VERTEX_SPARSE_SHARED_2];
#endif
#ifdef PP_N_VERTEX_SPARSE_SHARED_3
    t_vertexSparseShared3 (*m_vertexSparseShared3)[PP_N_VERTEX_SPARSE_SHARED_3];
#endif
#ifdef PP_N_VERTEX_SPARSE_SHARED_4
    t_vertexSparseShared4 (*m_vertexSparseShared4)[PP_N_VERTEX_SPARSE_SHARED_4];
#endif
#ifdef PP_N_VERTEX_SPARSE_SHARED_5
    t_vertexSparseShared5 (*m_vertexSparseShared5)[PP_N_VERTEX_SPARSE_SHARED_5];
#endif

    // sparse face data for every concurrent forward run
#ifdef PP_N_FACE_SPARSE_PRIVATE_1
    t_faceSparsePrivate1 (*m_faceSparsePrivate1)[PP_N_FACE_SPARSE_PRIVATE_1][N_CRUNS];
#endif
#ifdef PP_N_FACE_SPARSE_PRIVATE_2
    t_faceSparsePrivate2 (*m_faceSparsePrivate2)[PP_N_FACE_SPARSE_PRIVATE_2][N_CRUNS];
#endif
#ifdef PP_N_FACE_SPARSE_PRIVATE_3
    t_faceSparsePrivate3 (*m_faceSparsePrivate3)[PP_N_FACE_SPARSE_PRIVATE_3][N_CRUNS];
#endif
#ifdef PP_N_FACE_SPARSE_PRIVATE_4
    t_faceSparsePrivate4 (*m_faceSparsePrivate4)[PP_N_FACE_SPARSE_PRIVATE_4][N_CRUNS];
#endif
#ifdef PP_N_FACE_SPARSE_PRIVATE_5
    t_faceSparsePrivate5 (*m_faceSparsePrivate5)[PP_N_FACE_SPARSE_PRIVATE_5][N_CRUNS];
#endif
    // shared face data shared among the concurrent forward runs
#ifdef PP_N_FACE_SPARSE_SHARED_1
    t_faceSparseShared1 (*m_faceSparseShared1)[PP_N_FACE_SPARSE_SHARED_1];
#endif
#ifdef PP_N_FACE_SPARSE_SHARED_2
    t_faceSparseShared2 (*m_faceSparseShared2)[PP_N_FACE_SPARSE_SHARED_2];
#endif
#ifdef PP_N_FACE_SPARSE_SHARED_3
    t_faceSparseShared3 (*m_faceSparseShared3)[PP_N_FACE_SPARSE_SHARED_3];
#endif
#ifdef PP_N_FACE_SPARSE_SHARED_4
    t_faceSparseShared4 (*m_faceSparseShared4)[PP_N_FACE_SPARSE_SHARED_4];
#endif
#ifdef PP_N_FACE_SPARSE_SHARED_5
    t_faceSparseShared5 (*m_faceSparseShared5)[PP_N_FACE_SPARSE_SHARED_5];
#endif
#ifdef PP_N_FACE_SPARSE_SHARED_6
    t_faceSparseShared6 (*m_faceSparseShared6)[PP_N_FACE_SPARSE_SHARED_6];
#endif

    // sparse element data for every concurrent forward run
#ifdef PP_N_ELEMENT_SPARSE_PRIVATE_1
    t_elementSparsePrivate1 (*m_elementSparsePrivate1)[PP_N_ELEMENT_SPARSE_PRIVATE_1][N_CRUNS];
#endif
#ifdef PP_N_ELEMENT_SPARSE_PRIVATE_2
    t_elementSparsePrivate2 (*m_elementSparsePrivate2)[PP_N_ELEMENT_SPARSE_PRIVATE_2][N_CRUNS];
#endif
#ifdef PP_N_ELEMENT_SPARSE_PRIVATE_3
    t_elementSparsePrivate3 (*m_elementSparsePrivate3)[PP_N_ELEMENT_SPARSE_PRIVATE_3][N_CRUNS];
#endif
#ifdef PP_N_ELEMENT_SPARSE_PRIVATE_4
    t_elementSparsePrivate4 (*m_elementSparsePrivate4)[PP_N_ELEMENT_SPARSE_PRIVATE_4][N_CRUNS];
#endif
#ifdef PP_N_ELEMENT_SPARSE_PRIVATE_5
    t_elementSparsePrivate5 (*m_elementSparsePrivate5)[PP_N_ELEMENT_SPARSE_PRIVATE_5][N_CRUNS];
#endif
    // sparse element data shared among the concurrent forward runs
#ifdef PP_N_ELEMENT_SPARSE_SHARED_1
    t_elementSparseShared1 (*m_elementSparseShared1)[PP_N_ELEMENT_SPARSE_SHARED_1];
#endif
#ifdef PP_N_ELEMENT_SPARSE_SHARED_2
    t_elementSparseShared2 (*m_elementSparseShared2)[PP_N_ELEMENT_SPARSE_SHARED_2];
#endif
#ifdef PP_N_ELEMENT_SPARSE_SHARED_3
    t_elementSparseShared3 (*m_elementSparseShared3)[PP_N_ELEMENT_SPARSE_SHARED_3];
#endif
#ifdef PP_N_ELEMENT_SPARSE_SHARED_4
    t_elementSparseShared4 (*m_elementSparseShared4)[PP_N_ELEMENT_SPARSE_SHARED_4];
#endif
#ifdef PP_N_ELEMENT_SPARSE_SHARED_5
    t_elementSparseShared5 (*m_elementSparseShared5)[PP_N_ELEMENT_SPARSE_SHARED_5];
#endif
#ifdef PP_N_ELEMENT_SPARSE_SHARED_6
    t_elementSparseShared6 (*m_elementSparseShared6)[PP_N_ELEMENT_SPARSE_SHARED_6];
#endif

    // dense data for every mode and concurrent run
#ifdef PP_N_ELEMENT_MODE_PRIVATE_1
    t_elementModePrivate1 (*m_elementModePrivate1)[PP_N_ELEMENT_MODE_PRIVATE_1][N_ELEMENT_MODES][N_CRUNS];
#endif
#ifdef PP_N_ELEMENT_MODE_PRIVATE_2
    t_elementModePrivate2 (*m_elementModePrivate2)[PP_N_ELEMENT_MODE_PRIVATE_2][N_ELEMENT_MODES][N_CRUNS];
#endif
#ifdef PP_N_ELEMENT_MODE_PRIVATE_3
    t_elementModePrivate3 (*m_elementModePrivate3)[PP_N_ELEMENT_MODE_PRIVATE_3][N_ELEMENT_MODES][N_CRUNS];
#endif

    // dense data for every mode, shared among the concurrent runs
#ifdef PP_N_ELEMENT_MODE_SHARED_1
    t_elementModeShared1 (*m_elementModeShared1)[PP_N_ELEMENT_MODE_SHARED_1][N_ELEMENT_MODES];
#endif
#ifdef PP_N_ELEMENT_MODE_SHARED_2
    t_elementModeShared2 (*m_elementModeShared2)[PP_N_ELEMENT_MODE_SHARED_2][N_ELEMENT_MODES];
#endif
#ifdef PP_N_ELEMENT_MODE_SHARED_3
    t_elementModeShared3 (*m_elementModeShared3)[PP_N_ELEMENT_MODE_SHARED_3][N_ELEMENT_MODES];
#endif

    // dense data for every concurrent run
#ifdef PP_N_ELEMENT_PRIVATE_1
    t_elementPrivate1 (*m_elementPrivate1)[PP_N_ELEMENT_PRIVATE_1][N_CRUNS];
#endif
#ifdef PP_N_ELEMENT_PRIVATE_2
    t_elementPrivate2 (*m_elementPrivate2)[PP_N_ELEMENT_PRIVATE_2][N_CRUNS];
#endif
#ifdef PP_N_ELEMENT_PRIVATE_3
    t_elementPrivate3 (*m_elementPrivate3)[PP_N_ELEMENT_PRIVATE_3][N_CRUNS];
#endif

    // dense data for every element, shared among the runs
#ifdef PP_N_ELEMENT_SHARED_1
    t_elementShared1 (*m_elementShared1)[PP_N_ELEMENT_SHARED_1];
#endif
#ifdef PP_N_ELEMENT_SHARED_2
    t_elementShared2 (*m_elementShared2)[PP_N_ELEMENT_SHARED_2];
#endif
#ifdef PP_N_ELEMENT_SHARED_3
    t_elementShared3 (*m_elementShared3)[PP_N_ELEMENT_SHARED_3];
#endif
#ifdef PP_N_ELEMENT_SHARED_4
    t_elementShared4 (*m_elementShared4)[PP_N_ELEMENT_SHARED_4];
#endif

    // dense data for every face and concurrent run
#ifdef PP_N_FACE_MODE_PRIVATE_1
    t_faceModePrivate1 (*m_faceModePrivate1)[PP_N_FACE_MODE_PRIVATE_1][N_FACE_MODES][N_CRUNS];
#endif
#ifdef PP_N_FACE_MODE_PRIVATE_2
    t_faceModePrivate2 (*m_faceModePrivate2)[PP_N_FACE_MODE_PRIVATE_2][N_FACE_MODES][N_CRUNS];
#endif
#ifdef PP_N_FACE_MODE_PRIVATE_3
    t_faceModePrivate3 (*m_faceModePrivate3)[PP_N_FACE_MODE_PRIVATE_3][N_FACE_MODES][N_CRUNS];
#endif

    //! vertex characteristics
    t_vertexChars *m_vertexChars;

    //! face characteristics
    t_faceChars *m_faceChars;

    //! element characteristics
    t_elementChars *m_elementChars;

    // connectivity information
    t_connect m_connect;

    /**
     * Constructor.
     **/
    Internal(): m_initScratch(false), m_initDense(false), m_initSparse(false) {}

    /**
     * Initializes the scratch memory.
     **/
    void initScratch() {
      // init scratch memory if requested
#ifdef PP_SCRATCH_MEMORY
      m_initScratch = true;

#ifdef PP_USE_OMP
#pragma omp parallel for
#endif
      for( int l_td = 0; l_td < parallel::g_nThreads; l_td++ )  {
#ifdef PP_USE_OMP
#pragma omp critical
#endif
        if( parallel::g_scratchMem == nullptr ) parallel::g_scratchMem = (t_scratchMem*) common::allocate( sizeof(t_scratchMem),
                                                                                                           ALIGNMENT.BASE.HEAP,
                                                                                                           m_memTypes.hbw.scratchMem,
                                                                                                           m_memTypes.huge.scratchMem );
      }

#endif
    }

    /**
     * Initializes the dense data structures.
     *
     * @param i_nVertices number of vertices.
     * @param i_nFaces number of faces.
     * @param i_nElements number of elements.
     **/
    void initDense( int_el i_nVertices,
                    int_el i_nFaces,
                    int_el i_nElements ) {
      m_initDense = true;

      m_nVertices  = i_nVertices;
      m_nFaces     = i_nFaces;
      m_nElements  = i_nElements;

      // check for non-empty mesh
      EDGE_CHECK( m_nVertices > 0 );
      EDGE_CHECK( m_nFaces    > 0 );
      EDGE_CHECK( m_nElements > 0 );

      // faces
#ifdef PP_N_FACE_MODE_PRIVATE_1
      size_t l_faceModePrivateSize1      = std::size_t(m_nFaces) * PP_N_FACE_MODE_PRIVATE_1 * N_FACE_MODES    * N_CRUNS       * sizeof(t_faceModePrivate1);
#endif
#ifdef PP_N_FACE_MODE_PRIVATE_2
      size_t l_faceModePrivateSize2      = std::size_t(m_nFaces) * PP_N_FACE_MODE_PRIVATE_2 * N_FACE_MODES    * N_CRUNS       * sizeof(t_faceModePrivate2);
#endif
#ifdef PP_N_FACE_MODE_PRIVATE_3
      size_t l_faceModePrivateSize3      = std::size_t(m_nFaces) * PP_N_FACE_MODE_PRIVATE_3 * N_FACE_MODES    * N_CRUNS       * sizeof(t_faceModePrivate3);
#endif

      // elements
#ifdef PP_N_ELEMENT_MODE_PRIVATE_1
      size_t l_elementModePrivateSize1   = std::size_t(m_nElements) * PP_N_ELEMENT_MODE_PRIVATE_1 * N_ELEMENT_MODES * N_CRUNS * sizeof(t_elementModePrivate1);
#endif
#ifdef PP_N_ELEMENT_MODE_PRIVATE_2
      size_t l_elementModePrivateSize2   = std::size_t(m_nElements) * PP_N_ELEMENT_MODE_PRIVATE_2 * N_ELEMENT_MODES * N_CRUNS * sizeof(t_elementModePrivate2);
#endif
#ifdef PP_N_ELEMENT_MODE_PRIVATE_3
      size_t l_elementModePrivateSize3   = std::size_t(m_nElements) * PP_N_ELEMENT_MODE_PRIVATE_3 * N_ELEMENT_MODES * N_CRUNS * sizeof(t_elementModePrivate3);
#endif


#ifdef PP_N_ELEMENT_MODE_SHARED_1
      size_t l_elementModeSharedSize1    = std::size_t(m_nElements) * PP_N_ELEMENT_MODE_SHARED_1  * N_ELEMENT_MODES           * sizeof(t_elementModeShared1);
#endif
#ifdef PP_N_ELEMENT_MODE_SHARED_2
      size_t l_elementModeSharedSize2    = std::size_t(m_nElements) * PP_N_ELEMENT_MODE_SHARED_2  * N_ELEMENT_MODES           * sizeof(t_elementModeShared2);
#endif
#ifdef PP_N_ELEMENT_MODE_SHARED_3
      size_t l_elementModeSharedSize3    = std::size_t(m_nElements) * PP_N_ELEMENT_MODE_SHARED_3  * N_ELEMENT_MODES           * sizeof(t_elementModeShared3);
#endif

#ifdef PP_N_ELEMENT_PRIVATE_1
      size_t l_elementPrivateSize1       = std::size_t(m_nElements) * PP_N_ELEMENT_PRIVATE_1 * N_CRUNS                        * sizeof(t_elementPrivate1);
#endif
#ifdef PP_N_ELEMENT_PRIVATE_2
      size_t l_elementPrivateSize2       = std::size_t(m_nElements) * PP_N_ELEMENT_PRIVATE_2 * N_CRUNS                        * sizeof(t_elementPrivate2);
#endif
#ifdef PP_N_ELEMENT_PRIVATE_3
      size_t l_elementPrivateSize3       = std::size_t(m_nElements) * PP_N_ELEMENT_PRIVATE_3 * N_CRUNS                        * sizeof(t_elementPrivate3);
#endif

#ifdef PP_N_ELEMENT_SHARED_1
      size_t l_elementSharedSize1        = std::size_t(m_nElements) * PP_N_ELEMENT_SHARED_1                                   * sizeof(t_elementShared1);
#endif
#ifdef PP_N_ELEMENT_SHARED_2
      size_t l_elementSharedSize2        = std::size_t(m_nElements) * PP_N_ELEMENT_SHARED_2                                   * sizeof(t_elementShared2);
#endif
#ifdef PP_N_ELEMENT_SHARED_3
      size_t l_elementSharedSize3        = std::size_t(m_nElements) * PP_N_ELEMENT_SHARED_3                                   * sizeof(t_elementShared3);
#endif
#ifdef PP_N_ELEMENT_SHARED_4
      size_t l_elementSharedSize4        = std::size_t(m_nElements) * PP_N_ELEMENT_SHARED_4                                   * sizeof(t_elementShared4);
#endif

      size_t l_connElVeSize              =  std::size_t(m_nElements)      * C_ENT[T_SDISC.ELEMENT].N_VERTICES                 * sizeof(int_el);
      size_t l_connElVeElSize            =  std::size_t((m_nElements + 1))                                                    * sizeof(int_el*);
      size_t l_connElFaElSize            =  std::size_t(m_nElements)      * C_ENT[T_SDISC.ELEMENT].N_FACES                    * sizeof(int_el);
      size_t l_connElFaSize              =  std::size_t(m_nElements)      * C_ENT[T_SDISC.ELEMENT].N_FACES                    * sizeof(int_el);
      size_t l_connFaElSize              =  std::size_t(m_nFaces)         * 2                                                 * sizeof(int_el);
      size_t l_connFaVeSize              =  std::size_t(m_nFaces)         * C_ENT[T_SDISC.FACE].N_VERTICES                    * sizeof(int_el);
      size_t l_connfIdElFaElSize         =  std::size_t(m_nElements)      * C_ENT[T_SDISC.ELEMENT].N_FACES                    * sizeof(unsigned short);
      size_t l_connvIdElFaElSize         =  std::size_t(m_nElements)      * C_ENT[T_SDISC.ELEMENT].N_FACES                    * sizeof(unsigned short);

      size_t l_vertexCharsSize           =  std::size_t(m_nVertices)                                                          * sizeof(t_vertexChars);
      size_t l_faceCharsSize             =  std::size_t(m_nFaces)                                                             * sizeof(t_faceChars);
      size_t l_elementCharsSize          =  std::size_t(m_nElements)                                                          * sizeof(t_elementChars);

      // faces
#ifdef PP_N_FACE_MODE_PRIVATE_1
      m_faceModePrivate1      = (t_faceModePrivate1   (*)[PP_N_FACE_MODE_PRIVATE_1][N_FACE_MODES][N_CRUNS]       ) common::allocate( l_faceModePrivateSize1,
                                                                                                                                     ALIGNMENT.BASE.HEAP );
#endif
#ifdef PP_N_FACE_MODE_PRIVATE_2
      m_faceModePrivate2      = (t_faceModePrivate2   (*)[PP_N_FACE_MODE_PRIVATE_2][N_FACE_MODES][N_CRUNS]       ) common::allocate( l_faceModePrivateSize2,
                                                                                                                                     ALIGNMENT.BASE.HEAP );
#endif
#ifdef PP_N_FACE_MODE_PRIVATE_3
      m_faceModePrivate3      = (t_faceModePrivate3   (*)[PP_N_FACE_MODE_PRIVATE_3][N_FACE_MODES][N_CRUNS]       ) common::allocate( l_faceModePrivateSize3,
                                                                                                                                     ALIGNMENT.BASE.HEAP );
#endif


      // elements
#ifdef PP_N_ELEMENT_MODE_PRIVATE_1
      m_elementModePrivate1   = (t_elementModePrivate1 (*)[PP_N_ELEMENT_MODE_PRIVATE_1][N_ELEMENT_MODES][N_CRUNS]) common::allocate( l_elementModePrivateSize1,
                                                                                                                                     ALIGNMENT.BASE.HEAP,
                                                                                                                                     m_memTypes.hbw.elementModePrivate1,
                                                                                                                                     m_memTypes.huge.elementModePrivate1 );
#endif
#ifdef PP_N_ELEMENT_MODE_PRIVATE_2
      m_elementModePrivate2   = (t_elementModePrivate2 (*)[PP_N_ELEMENT_MODE_PRIVATE_2][N_ELEMENT_MODES][N_CRUNS]) common::allocate( l_elementModePrivateSize2,
                                                                                                                                     ALIGNMENT.BASE.HEAP,
                                                                                                                                     m_memTypes.hbw.elementModePrivate2,
                                                                                                                                     m_memTypes.huge.elementModePrivate2 );
#endif
#ifdef PP_N_ELEMENT_MODE_PRIVATE_3
      m_elementModePrivate3   = (t_elementModePrivate3 (*)[PP_N_ELEMENT_MODE_PRIVATE_3][N_ELEMENT_MODES][N_CRUNS]) common::allocate( l_elementModePrivateSize3,
                                                                                                                                     ALIGNMENT.BASE.HEAP,
                                                                                                                                     m_memTypes.hbw.elementModePrivate3,
                                                                                                                                     m_memTypes.huge.elementModePrivate3 );
#endif


#ifdef PP_N_ELEMENT_MODE_SHARED_1
      m_elementModeShared1    = (t_elementModeShared1 (*)[PP_N_ELEMENT_MODE_SHARED_1][N_ELEMENT_MODES]           ) common::allocate( l_elementModeSharedSize1,
                                                                                                                                     ALIGNMENT.BASE.HEAP );
#endif
#ifdef PP_N_ELEMENT_MODE_SHARED_2
      m_elementModeShared2    = (t_elementModeShared2 (*)[PP_N_ELEMENT_MODE_SHARED_2][N_ELEMENT_MODES]           ) common::allocate( l_elementModeSharedSize2,
                                                                                                                                     ALIGNMENT.BASE.HEAP );
#endif
#ifdef PP_N_ELEMENT_MODE_SHARED_3
      m_elementModeShared3    = (t_elementModeShared3 (*)[PP_N_ELEMENT_MODE_SHARED_3][N_ELEMENT_MODES]           ) common::allocate( l_elementModeSharedSize3,
                                                                                                                                     ALIGNMENT.BASE.HEAP );
#endif


#ifdef PP_N_ELEMENT_PRIVATE_1
      m_elementPrivate1   = (t_elementPrivate1 (*)[PP_N_ELEMENT_PRIVATE_1][N_CRUNS]) common::allocate( l_elementPrivateSize1,
                                                                                                       ALIGNMENT.BASE.HEAP,
                                                                                                       m_memTypes.hbw.elementPrivate1 );
#endif
#ifdef PP_N_ELEMENT_PRIVATE_2
      m_elementPrivate2   = (t_elementPrivate2 (*)[PP_N_ELEMENT_PRIVATE_2][N_CRUNS]) common::allocate( l_elementPrivateSize2,
                                                                                                       ALIGNMENT.BASE.HEAP,
                                                                                                       m_memTypes.hbw.elementPrivate2 );
#endif
#ifdef PP_N_ELEMENT_PRIVATE_3
      m_elementPrivate3   = (t_elementPrivate3 (*)[PP_N_ELEMENT_PRIVATE_3][N_CRUNS]) common::allocate( l_elementPrivateSize3,
                                                                                                       ALIGNMENT.BASE.HEAP,
                                                                                                       m_memTypes.hbw.elementPrivate3 );
#endif


#ifdef PP_N_ELEMENT_SHARED_1
      m_elementShared1        = (t_elementShared1 (*)[PP_N_ELEMENT_SHARED_1]                                     ) common::allocate( l_elementSharedSize1,
                                                                                                                                     ALIGNMENT.BASE.HEAP,
                                                                                                                                     m_memTypes.hbw.elementShared1,
                                                                                                                                     m_memTypes.huge.elementShared1 );
#endif
#ifdef PP_N_ELEMENT_SHARED_2
      m_elementShared2        = (t_elementShared2 (*)[PP_N_ELEMENT_SHARED_2]                                     ) common::allocate( l_elementSharedSize2,
                                                                                                                                     ALIGNMENT.BASE.HEAP,
                                                                                                                                     m_memTypes.hbw.elementShared2,
                                                                                                                                     m_memTypes.huge.elementShared2 );
#endif
#ifdef PP_N_ELEMENT_SHARED_3
      m_elementShared3        = (t_elementShared3 (*)[PP_N_ELEMENT_SHARED_3]                                     ) common::allocate( l_elementSharedSize3,
                                                                                                                                     ALIGNMENT.BASE.HEAP,
                                                                                                                                     m_memTypes.hbw.elementShared3,
                                                                                                                                     m_memTypes.huge.elementShared3 );
#endif
#ifdef PP_N_ELEMENT_SHARED_4
      m_elementShared4        = (t_elementShared4 (*)[PP_N_ELEMENT_SHARED_4]                                     ) common::allocate( l_elementSharedSize4,
                                                                                                                                     ALIGNMENT.BASE.HEAP,
                                                                                                                                     m_memTypes.hbw.elementShared4,
                                                                                                                                     m_memTypes.huge.elementShared4 );
#endif

      m_connect.elVe            = (int_el         (*)[C_ENT[T_SDISC.ELEMENT].N_VERTICES]                         ) common::allocate( l_connElVeSize,
                                                                                                                                     ALIGNMENT.BASE.HEAP );
      m_connect.elVeEl          = (int_el**                                                                      ) common::allocate( l_connElVeElSize,
                                                                                                                                     ALIGNMENT.BASE.HEAP );
      m_connect.elFaEl          = (int_el         (*)[C_ENT[T_SDISC.ELEMENT].N_FACES]                            ) common::allocate( l_connElFaElSize,
                                                                                                                                     ALIGNMENT.BASE.HEAP );
      m_connect.elFa            = (int_el         (*)[C_ENT[T_SDISC.ELEMENT].N_FACES]                            ) common::allocate( l_connElFaSize,
                                                                                                                                     ALIGNMENT.BASE.HEAP );
      m_connect.faEl            = (int_el         (*)[2]                                                         ) common::allocate( l_connFaElSize,
                                                                                                                                     ALIGNMENT.BASE.HEAP );
      m_connect.faVe            = (int_el         (*)[C_ENT[T_SDISC.FACE].N_VERTICES]                            ) common::allocate( l_connFaVeSize,
                                                                                                                                     ALIGNMENT.BASE.HEAP );
      m_connect.fIdElFaEl       = (unsigned short (*)[C_ENT[T_SDISC.ELEMENT].N_FACES]                            ) common::allocate( l_connfIdElFaElSize,
                                                                                                                                     ALIGNMENT.BASE.HEAP );
      m_connect.vIdElFaEl       = (unsigned short (*)[C_ENT[T_SDISC.ELEMENT].N_FACES]                            ) common::allocate( l_connvIdElFaElSize,
                                                                                                                                     ALIGNMENT.BASE.HEAP );
      m_vertexChars             = (t_vertexChars  *                                                              ) common::allocate( l_vertexCharsSize,
                                                                                                                                     ALIGNMENT.BASE.HEAP );
      m_faceChars               = (t_faceChars    *                                                              ) common::allocate( l_faceCharsSize,
                                                                                                                                     ALIGNMENT.BASE.HEAP );
      m_elementChars            = (t_elementChars *                                                              ) common::allocate( l_elementCharsSize,
                                                                                                                                     ALIGNMENT.BASE.HEAP );
    }

    /**
     * Initializes the sparse data structures.
     *
     * @param i_nVeSp1 number of sparse vertices #1.
     * @param i_nFaSp1 number of sparse faces #1.
     * @param i_nElSp1 number of sparse elements #1.
     * @param i_nVeSp2 number of sparse vertices #2.
     * @param i_nFaSp2 number of sparse faces #2.
     * @param i_nElSp2 number of sparse elements #2.
     * @param i_nVeSp3 number of sparse vertices #3.
     * @param i_nElSp3 number of sparse elements #3.
     * @param i_nfaSp3 number of sparse faces #3.
     * @param i_nVeSp4 number of sparse vertices #4.
     * @param i_nElSp4 number of sparse elements #4.
     * @param i_nfaSp4 number of sparse faces #4.
     * @param i_nVeSp5 number of sparse vertices #5.
     * @param i_nElSp5 number of sparse elements #5.
     * @param i_nfaSp5 number of sparse faces #5.
     * @param i_nVeSp6 number of sparse vertices #6.
     * @param i_nElSp6 number of sparse elements #6.
     * @param i_nfaSp6 number of sparse faces #6.
     **/
    void initSparse( int_el i_nVeSp1 = 0,
                     int_el i_nFaSp1 = 0,
                     int_el i_nElSp1 = 0,
                     int_el i_nVeSp2 = 0,
                     int_el i_nFaSp2 = 0,
                     int_el i_nElSp2 = 0,
                     int_el i_nVeSp3 = 0,
                     int_el i_nFaSp3 = 0,
                     int_el i_nElSp3 = 0,
                     int_el i_nVeSp4 = 0,
                     int_el i_nFaSp4 = 0,
                     int_el i_nElSp4 = 0,
                     int_el i_nVeSp5 = 0,
                     int_el i_nFaSp5 = 0,
                     int_el i_nElSp5 = 0,
                     int_el i_nVeSp6 = 0,
                     int_el i_nFaSp6 = 0,
                     int_el i_nElSp6 = 0 ) {
      m_initSparse = true;

      m_nVeSp1 = i_nVeSp1;
      m_nVeSp2 = i_nVeSp2;
      m_nVeSp3 = i_nVeSp3;
      m_nVeSp4 = i_nVeSp4;
      m_nVeSp5 = i_nVeSp5;
      m_nVeSp6 = i_nVeSp6;

      m_nFaSp1 = i_nFaSp1;
      m_nFaSp2 = i_nFaSp2;
      m_nFaSp3 = i_nFaSp3;
      m_nFaSp4 = i_nFaSp4;
      m_nFaSp5 = i_nFaSp5;
      m_nFaSp6 = i_nFaSp6;

      m_nElSp1 = i_nElSp1;
      m_nElSp2 = i_nElSp2;
      m_nElSp3 = i_nElSp3;
      m_nElSp4 = i_nElSp4;
      m_nElSp5 = i_nElSp5;
      m_nElSp6 = i_nElSp6;

      // vertices
#ifdef PP_N_VERTEX_SPARSE_PRIVATE_1
      size_t l_veSpPrivate1              = std::size_t(m_nVeSp1) * PP_N_VERTEX_SPARSE_1 * N_CRUNS                            * sizeof(t_vertexSparsePrivate1);
#endif
#ifdef PP_N_VERTEX_SPARSE_PRIVATE_2
      size_t l_veSpPrivate2              = std::size_t(m_nVeSp2) * PP_N_VERTEX_SPARSE_2 * N_CRUNS                            * sizeof(t_vertexSparsePrivate2);
#endif
#ifdef PP_N_VERTEX_SPARSE_PRIVATE_3
      size_t l_veSpPrivate3              = std::size_t(m_nVeSp3) * PP_N_VERTEX_SPARSE_3 * N_CRUNS                            * sizeof(t_vertexSparsePrivate3);
#endif
#ifdef PP_N_VERTEX_SPARSE_PRIVATE_4
      size_t l_veSpPrivate4              = std::size_t(m_nVeSp4) * PP_N_VERTEX_SPARSE_4 * N_CRUNS                            * sizeof(t_vertexSparsePrivate4);
#endif
#ifdef PP_N_VERTEX_SPARSE_PRIVATE_5
      size_t l_veSpPrivate5              = std::size_t(m_nVeSp5) * PP_N_VERTEX_SPARSE_5 * N_CRUNS                            * sizeof(t_vertexSparsePrivate5);
#endif

#ifdef PP_N_VERTEX_SPARSE_SHARED_1
      size_t l_veSpShared1               = std::size_t(m_nVeSp1) * PP_N_VERTEX_SPARSE_1                                      * sizeof(t_vertexSparseShared1);
#endif
#ifdef PP_N_VERTEX_SPARSE_SHARED_2
      size_t l_veSpShared2               = std::size_t(m_nVeSp2) * PP_N_VERTEX_SPARSE_2                                      * sizeof(t_vertexSparseShared2);
#endif
#ifdef PP_N_VERTEX_SPARSE_SHARED_3
      size_t l_veSpShared3               = std::size_t(m_nVeSp3) * PP_N_VERTEX_SPARSE_3                                      * sizeof(t_vertexSparseShared3);
#endif
#ifdef PP_N_VERTEX_SPARSE_SHARED_4
      size_t l_veSpShared4               = std::size_t(m_nVeSp4) * PP_N_VERTEX_SPARSE_4                                      * sizeof(t_vertexSparseShared4);
#endif
#ifdef PP_N_VERTEX_SPARSE_SHARED_5
      size_t l_veSpShared5               = std::size_t(m_nVeSp5) * PP_N_VERTEX_SPARSE_5                                      * sizeof(t_vertexSparseShared5);
#endif

      // faces
#ifdef PP_N_FACE_SPARSE_PRIVATE_1
      size_t l_faSpPrivate1              = std::size_t(m_nFaSp1) * PP_N_FACE_SPARSE_PRIVATE_1 * N_CRUNS                       * sizeof(t_faceSparsePrivate1);
#endif
#ifdef PP_N_FACE_SPARSE_PRIVATE_2
      size_t l_faSpPrivate2              = std::size_t(m_nFaSp2) * PP_N_FACE_SPARSE_PRIVATE_2 * N_CRUNS                       * sizeof(t_faceSparsePrivate2);
#endif
#ifdef PP_N_FACE_SPARSE_PRIVATE_3
      size_t l_faSpPrivate3              = std::size_t(m_nFaSp3) * PP_N_FACE_SPARSE_PRIVATE_3 * N_CRUNS                       * sizeof(t_faceSparsePrivate3);
#endif
#ifdef PP_N_FACE_SPARSE_PRIVATE_4
      size_t l_faSpPrivate4              = std::size_t(m_nFaSp4) * PP_N_FACE_SPARSE_PRIVATE_4 * N_CRUNS                       * sizeof(t_faceSparsePrivate4);
#endif
#ifdef PP_N_FACE_SPARSE_PRIVATE_5
      size_t l_faSpPrivate5              = std::size_t(m_nFaSp5) * PP_N_FACE_SPARSE_PRIVATE_5 * N_CRUNS                       * sizeof(t_faceSparsePrivate5);
#endif

#ifdef PP_N_FACE_SPARSE_SHARED_1
      size_t l_faSpShared1               = std::size_t(m_nFaSp1) * PP_N_FACE_SPARSE_SHARED_1                                  * sizeof(t_faceSparseShared1);
#endif
#ifdef PP_N_FACE_SPARSE_SHARED_2
      size_t l_faSpShared2               = std::size_t(m_nFaSp2) * PP_N_FACE_SPARSE_SHARED_2                                  * sizeof(t_faceSparseShared2);
#endif
#ifdef PP_N_FACE_SPARSE_SHARED_3
      size_t l_faSpShared3               = std::size_t(m_nFaSp3) * PP_N_FACE_SPARSE_SHARED_3                                  * sizeof(t_faceSparseShared3);
#endif
#ifdef PP_N_FACE_SPARSE_SHARED_4
      size_t l_faSpShared4               = std::size_t(m_nFaSp4) * PP_N_FACE_SPARSE_SHARED_4                                  * sizeof(t_faceSparseShared4);
#endif
#ifdef PP_N_FACE_SPARSE_SHARED_5
      size_t l_faSpShared5               = std::size_t(m_nFaSp5) * PP_N_FACE_SPARSE_SHARED_5                                  * sizeof(t_faceSparseShared5);
#endif
#ifdef PP_N_FACE_SPARSE_SHARED_6
      size_t l_faSpShared6               = std::size_t(m_nFaSp6) * PP_N_FACE_SPARSE_SHARED_6                                  * sizeof(t_faceSparseShared6);
#endif

      // elements
#ifdef PP_N_ELEMENT_SPARSE_PRIVATE_1
      size_t l_elSpPrivate1              = std::size_t(m_nElSp1) * PP_N_ELEMENT_SPARSE_PRIVATE_1 * N_CRUNS                    * sizeof(t_elementSparsePrivate1);
#endif
#ifdef PP_N_ELEMENT_SPARSE_PRIVATE_2
      size_t l_elSpPrivate2              = std::size_t(m_nElSp2) * PP_N_ELEMENT_SPARSE_PRIVATE_2 * N_CRUNS                    * sizeof(t_elementSparsePrivate2);
#endif
#ifdef PP_N_ELEMENT_SPARSE_PRIVATE_3
      size_t l_elSpPrivate3              = std::size_t(m_nElSp3) * PP_N_ELEMENT_SPARSE_PRIVATE_3 * N_CRUNS                    * sizeof(t_elementSparsePrivate3);
#endif
#ifdef PP_N_ELEMENT_SPARSE_PRIVATE_4
      size_t l_elSpPrivate4              = std::size_t(m_nElSp4) * PP_N_ELEMENT_SPARSE_PRIVATE_4 * N_CRUNS                    * sizeof(t_elementSparsePrivate4);
#endif
#ifdef PP_N_ELEMENT_SPARSE_PRIVATE_5
      size_t l_elSpPrivate5              = std::size_t(m_nElSp5) * PP_N_ELEMENT_SPARSE_PRIVATE_5 * N_CRUNS                    * sizeof(t_elementSparsePrivate5);
#endif

#ifdef PP_N_ELEMENT_SPARSE_SHARED_1
      size_t l_elSpShared1               = std::size_t(m_nElSp1) * PP_N_ELEMENT_SPARSE_SHARED_1                               * sizeof(t_elementSparseShared1);
#endif
#ifdef PP_N_ELEMENT_SPARSE_SHARED_2
      size_t l_elSpShared2               = std::size_t(m_nElSp2) * PP_N_ELEMENT_SPARSE_SHARED_2                               * sizeof(t_elementSparseShared2);
#endif
#ifdef PP_N_ELEMENT_SPARSE_SHARED_3
      size_t l_elSpShared3               = std::size_t(m_nElSp3) * PP_N_ELEMENT_SPARSE_SHARED_3                               * sizeof(t_elementSparseShared3);
#endif
#ifdef PP_N_ELEMENT_SPARSE_SHARED_4
      size_t l_elSpShared4               = std::size_t(m_nElSp4) * PP_N_ELEMENT_SPARSE_SHARED_4                               * sizeof(t_elementSparseShared4);
#endif
#ifdef PP_N_ELEMENT_SPARSE_SHARED_5
      size_t l_elSpShared5               = std::size_t(m_nElSp5) * PP_N_ELEMENT_SPARSE_SHARED_5                               * sizeof(t_elementSparseShared5);
#endif
#ifdef PP_N_ELEMENT_SPARSE_SHARED_6
      size_t l_elSpShared6               = std::size_t(m_nElSp6) * PP_N_ELEMENT_SPARSE_SHARED_6                               * sizeof(t_elementSparseShared6);
#endif


      // vertices
#ifdef PP_N_VERTEX_SPARSE_PRIVATE_1
      m_vertexSparsePrivate1  = (t_vertexSparsePrivate1 (*)[PP_N_VERTEX_SPARSE_PRIVATE_1][N_CRUNS]               ) common::allocate( l_veSpPrivate1,
                                                                                                                                     ALIGNMENT.BASE.HEAP );
#endif
#ifdef PP_N_VERTEX_SPARSE_PRIVATE_2
      m_vertexSparsePrivate2  = (t_vertexSparsePrivate2 (*)[PP_N_VERTEX_SPARSE_PRIVATE_2][N_CRUNS]               ) common::allocate( l_veSpPrivate2,
                                                                                                                                     ALIGNMENT.BASE.HEAP );
#endif
#ifdef PP_N_VERTEX_SPARSE_PRIVATE_3
      m_vertexSparsePrivate3  = (t_vertexSparsePrivate3 (*)[PP_N_VERTEX_SPARSE_PRIVATE_3][N_CRUNS]               ) common::allocate( l_veSpPrivate3,
                                                                                                                                     ALIGNMENT.BASE.HEAP );
#endif
#ifdef PP_N_VERTEX_SPARSE_PRIVATE_4
      m_vertexSparsePrivate4  = (t_vertexSparsePrivate4 (*)[PP_N_VERTEX_SPARSE_PRIVATE_4][N_CRUNS]               ) common::allocate( l_veSpPrivate4,
                                                                                                                                     ALIGNMENT.BASE.HEAP );
#endif
#ifdef PP_N_VERTEX_SPARSE_PRIVATE_5
      m_vertexSparsePrivate5  = (t_vertexSparsePrivate5 (*)[PP_N_VERTEX_SPARSE_PRIVATE_5][N_CRUNS]               ) common::allocate( l_veSpPrivate5,
                                                                                                                                     ALIGNMENT.BASE.HEAP );
#endif

#ifdef PP_N_VERTEX_SPARSE_SHARED_1
      m_vertexSparseShared1   = (t_vertexSparseShared1 (*)[PP_N_VERTEX_SPARSE_SHARED_1]                          ) common::allocate( l_veSpShared1,
                                                                                                                                     ALIGNMENT.BASE.HEAP );
#endif
#ifdef PP_N_VERTEX_SPARSE_SHARED_2
      m_vertexSparseShared2   = (t_vertexSparseShared2 (*)[PP_N_VERTEX_SPARSE_SHARED_1]                          ) common::allocate( l_veSpShared2,
                                                                                                                                     ALIGNMENT.BASE.HEAP );
#endif
#ifdef PP_N_VERTEX_SPARSE_SHARED_3
      m_vertexSparseShared3   = (t_vertexSparseShared3 (*)[PP_N_VERTEX_SPARSE_SHARED_1]                          ) common::allocate( l_veSpShared3,
                                                                                                                                     ALIGNMENT.BASE.HEAP );
#endif
#ifdef PP_N_VERTEX_SPARSE_SHARED_4
      m_vertexSparseShared4   = (t_vertexSparseShared4 (*)[PP_N_VERTEX_SPARSE_SHARED_1]                          ) common::allocate( l_veSpShared4,
                                                                                                                                     ALIGNMENT.BASE.HEAP );
#endif
#ifdef PP_N_VERTEX_SPARSE_SHARED_5
      m_vertexSparseShared5   = (t_vertexSparseShared5 (*)[PP_N_VERTEX_SPARSE_SHARED_1]                          ) common::allocate( l_veSpShared5,
                                                                                                                                     ALIGNMENT.BASE.HEAP );
#endif


      // faces
#ifdef PP_N_FACE_SPARSE_PRIVATE_1
      m_faceSparsePrivate1     = (t_faceSparsePrivate1 (*)[PP_N_FACE_SPARSE_PRIVATE_1] [N_CRUNS]                 ) common::allocate( l_faSpPrivate1,
                                                                                                                                     ALIGNMENT.BASE.HEAP );
#endif
#ifdef PP_N_FACE_SPARSE_PRIVATE_2
      m_faceSparsePrivate2     = (t_faceSparsePrivate2 (*)[PP_N_FACE_SPARSE_PRIVATE_2] [N_CRUNS]                 ) common::allocate( l_faSpPrivate2,
                                                                                                                                     ALIGNMENT.BASE.HEAP );
#endif
#ifdef PP_N_FACE_SPARSE_PRIVATE_3
      m_faceSparsePrivate3     = (t_faceSparsePrivate3 (*)[PP_N_FACE_SPARSE_PRIVATE_3] [N_CRUNS]                 ) common::allocate( l_faSpPrivate3,
                                                                                                                                     ALIGNMENT.BASE.HEAP );
#endif
#ifdef PP_N_FACE_SPARSE_PRIVATE_4
      m_faceSparsePrivate4     = (t_faceSparsePrivate4 (*)[PP_N_FACE_SPARSE_PRIVATE_4] [N_CRUNS]                 ) common::allocate( l_faSpPrivate4,
                                                                                                                                     ALIGNMENT.BASE.HEAP );
#endif
#ifdef PP_N_FACE_SPARSE_PRIVATE_5
      m_faceSparsePrivate5     = (t_faceSparsePrivate5 (*)[PP_N_FACE_SPARSE_PRIVATE_5] [N_CRUNS]                 ) common::allocate( l_faSpPrivate5,
                                                                                                                                     ALIGNMENT.BASE.HEAP );
#endif

#ifdef PP_N_FACE_SPARSE_SHARED_1
      m_faceSparseShared1      = (t_faceSparseShared1 (*)[PP_N_FACE_SPARSE_SHARED_1]                             ) common::allocate( l_faSpShared1,
                                                                                                                                     ALIGNMENT.BASE.HEAP );
#endif
#ifdef PP_N_FACE_SPARSE_SHARED_2
      m_faceSparseShared2      = (t_faceSparseShared2 (*)[PP_N_FACE_SPARSE_SHARED_2]                             ) common::allocate( l_faSpShared2,
                                                                                                                                     ALIGNMENT.BASE.HEAP );
#endif
#ifdef PP_N_FACE_SPARSE_SHARED_3
      m_faceSparseShared3      = (t_faceSparseShared3 (*)[PP_N_FACE_SPARSE_SHARED_3]                             ) common::allocate( l_faSpShared3,
                                                                                                                                     ALIGNMENT.BASE.HEAP );
#endif
#ifdef PP_N_FACE_SPARSE_SHARED_4
      m_faceSparseShared4      = (t_faceSparseShared4 (*)[PP_N_FACE_SPARSE_SHARED_4]                             ) common::allocate( l_faSpShared4,
                                                                                                                                     ALIGNMENT.BASE.HEAP );
#endif
#ifdef PP_N_FACE_SPARSE_SHARED_5
      m_faceSparseShared5      = (t_faceSparseShared5 (*)[PP_N_FACE_SPARSE_SHARED_5]                             ) common::allocate( l_faSpShared5,
                                                                                                                                     ALIGNMENT.BASE.HEAP );
#endif
#ifdef PP_N_FACE_SPARSE_SHARED_6
      m_faceSparseShared6      = (t_faceSparseShared6 (*)[PP_N_FACE_SPARSE_SHARED_6]                             ) common::allocate( l_faSpShared6,
                                                                                                                                     ALIGNMENT.BASE.HEAP );
#endif

     // elements
#ifdef PP_N_ELEMENT_SPARSE_PRIVATE_1
      m_elementSparsePrivate1     = (t_elementSparsePrivate1 (*)[PP_N_ELEMENT_SPARSE_PRIVATE_1][N_CRUNS]         ) common::allocate( l_elSpPrivate1,
                                                                                                                                     ALIGNMENT.BASE.HEAP );
#endif
#ifdef PP_N_ELEMENT_SPARSE_PRIVATE_2
      m_elementSparsePrivate2     = (t_elementSparsePrivate2 (*)[PP_N_ELEMENT_SPARSE_PRIVATE_2][N_CRUNS]         ) common::allocate( l_elSpPrivate2,
                                                                                                                                     ALIGNMENT.BASE.HEAP );
#endif
#ifdef PP_N_ELEMENT_SPARSE_PRIVATE_3
      m_elementSparsePrivate3     = (t_elementSparsePrivate3 (*)[PP_N_ELEMENT_SPARSE_PRIVATE_3][N_CRUNS]         ) common::allocate( l_elSpPrivate3,
                                                                                                                                     ALIGNMENT.BASE.HEAP );
#endif
#ifdef PP_N_ELEMENT_SPARSE_PRIVATE_4
      m_elementSparsePrivate4     = (t_elementSparsePrivate4 (*)[PP_N_ELEMENT_SPARSE_PRIVATE_4][N_CRUNS]         ) common::allocate( l_elSpPrivate4,
                                                                                                                                     ALIGNMENT.BASE.HEAP );
#endif
#ifdef PP_N_ELEMENT_SPARSE_PRIVATE_5
      m_elementSparsePrivate5     = (t_elementSparsePrivate5 (*)[PP_N_ELEMENT_SPARSE_PRIVATE_5][N_CRUNS]         ) common::allocate( l_elSpPrivate5,
                                                                                                                                     ALIGNMENT.BASE.HEAP );
#endif

#ifdef PP_N_ELEMENT_SPARSE_SHARED_1
      m_elementSparseShared1      = (t_elementSparseShared1 (*)[PP_N_ELEMENT_SPARSE_SHARED_1]                    ) common::allocate( l_elSpShared1,
                                                                                                                                     ALIGNMENT.BASE.HEAP );
#endif
#ifdef PP_N_ELEMENT_SPARSE_SHARED_2
      m_elementSparseShared2      = (t_elementSparseShared2 (*)[PP_N_ELEMENT_SPARSE_SHARED_2]                    ) common::allocate( l_elSpShared2,
                                                                                                                                     ALIGNMENT.BASE.HEAP );
#endif
#ifdef PP_N_ELEMENT_SPARSE_SHARED_3
      m_elementSparseShared3      = (t_elementSparseShared3 (*)[PP_N_ELEMENT_SPARSE_SHARED_3]                    ) common::allocate( l_elSpShared3,
                                                                                                                                     ALIGNMENT.BASE.HEAP );
#endif
#ifdef PP_N_ELEMENT_SPARSE_SHARED_4
      m_elementSparseShared4      = (t_elementSparseShared4 (*)[PP_N_ELEMENT_SPARSE_SHARED_4]                    ) common::allocate( l_elSpShared4,
                                                                                                                                     ALIGNMENT.BASE.HEAP );
#endif
#ifdef PP_N_ELEMENT_SPARSE_SHARED_5
      m_elementSparseShared5      = (t_elementSparseShared5 (*)[PP_N_ELEMENT_SPARSE_SHARED_5]                    ) common::allocate( l_elSpShared5,
                                                                                                                                     ALIGNMENT.BASE.HEAP );
#endif
#ifdef PP_N_ELEMENT_SPARSE_SHARED_6
      m_elementSparseShared6      = (t_elementSparseShared6 (*)[PP_N_ELEMENT_SPARSE_SHARED_6]                    ) common::allocate( l_elSpShared6,
                                                                                                                                     ALIGNMENT.BASE.HEAP );
#endif
    }

    /**
     * Finalizes the data structures.
     **/
    void finalize() {
      if( m_initScratch == true ) {
        // free scratch memory
#ifdef PP_SCRATCH_MEMORY

#ifdef PP_USE_OMP
#pragma omp parallel for
#endif
        for( int l_td = 0; l_td < parallel::g_nThreads; l_td++ ) {
#ifdef PP_USE_OMP
#pragma omp critical
#endif
          if( parallel::g_scratchMem != nullptr ) common::release( parallel::g_scratchMem,
                                                                   m_memTypes.hbw.scratchMem,
                                                                   m_memTypes.huge.scratchMem );
          parallel::g_scratchMem = nullptr;
        }

#endif
      }

      if( m_initDense == true ) {
        // faces
#ifdef PP_N_FACE_MODE_PRIVATE_1
        common::release(m_faceModePrivate1);
#endif
#ifdef PP_N_FACE_MODE_PRIVATE_2
        common::release(m_faceModePrivate2);
#endif
#ifdef PP_N_FACE_MODE_PRIVATE_3
        common::release(m_faceModePrivate3);
#endif

       // elements
#ifdef PP_N_ELEMENT_MODE_PRIVATE_1
       common::release(m_elementModePrivate1, m_memTypes.hbw.elementModePrivate1, m_memTypes.huge.elementModePrivate1);
#endif
#ifdef PP_N_ELEMENT_MODE_PRIVATE_2
       common::release(m_elementModePrivate2, m_memTypes.hbw.elementModePrivate2, m_memTypes.huge.elementModePrivate2);
#endif
#ifdef PP_N_ELEMENT_MODE_PRIVATE_3
       common::release(m_elementModePrivate3, m_memTypes.hbw.elementModePrivate3, m_memTypes.huge.elementModePrivate3);
#endif

#ifdef PP_N_ELEMENT_PRIVATE_1
       common::release(m_elementPrivate1, m_memTypes.hbw.elementPrivate1);
#endif
#ifdef PP_N_ELEMENT_PRIVATE_2
       common::release(m_elementPrivate2, m_memTypes.hbw.elementPrivate2);
#endif
#ifdef PP_N_ELEMENT_PRIVATE_3
       common::release(m_elementPrivate3, m_memTypes.hbw.elementPrivate3);
#endif

#ifdef PP_N_ELEMENT_MODE_SHARED_1
       common::release(m_elementModeShared1);
#endif
#ifdef PP_N_ELEMENT_MODE_SHARED_2
       common::release(m_elementModeShared2);
#endif
#ifdef PP_N_ELEMENT_MODE_SHARED_3
       common::release(m_elementModeShared3);
#endif

#ifdef PP_N_ELEMENT_SHARED_1
       common::release(m_elementShared1, m_memTypes.hbw.elementShared1, m_memTypes.huge.elementShared1);
#endif
#ifdef PP_N_ELEMENT_SHARED_2
       common::release(m_elementShared2, m_memTypes.hbw.elementShared2, m_memTypes.huge.elementShared2);
#endif
#ifdef PP_N_ELEMENT_SHARED_3
       common::release(m_elementShared3, m_memTypes.hbw.elementShared3, m_memTypes.huge.elementShared3);
#endif
#ifdef PP_N_ELEMENT_SHARED_4
       common::release(m_elementShared4, m_memTypes.hbw.elementShared4, m_memTypes.huge.elementShared4);
#endif

      common::release(m_connect.elVe);
      common::release(m_connect.elFa);
      common::release(m_connect.faEl);
      common::release(m_connect.faVe);
      common::release(m_connect.elVeEl);
      common::release(m_connect.elFaEl);
      common::release(m_connect.fIdElFaEl);
      common::release(m_connect.vIdElFaEl);

      common::release(m_vertexChars);
      common::release(m_faceChars);
      common::release(m_elementChars);
      }

      if( m_initSparse == true ) {
        // vertices
#ifdef PP_N_VERTEX_SPARSE_PRIVATE_1
        common::release(m_vertexSparsePrivate1);
#endif
#ifdef PP_N_VERTEX_SPARSE_PRIVATE_2
        common::release(m_vertexSparsePrivate2);
#endif
#ifdef PP_N_VERTEX_SPARSE_PRIVATE_3
        common::release(m_vertexSparsePrivate3);
#endif
#ifdef PP_N_VERTEX_SPARSE_PRIVATE_4
        common::release(m_vertexSparsePrivate4);
#endif
#ifdef PP_N_VERTEX_SPARSE_PRIVATE_5
        common::release(m_vertexSparsePrivate5);
#endif
#ifdef PP_N_VERTEX_SPARSE_SHARED_1
        common::release(m_vertexSparseShared1);
#endif
#ifdef PP_N_VERTEX_SPARSE_SHARED_2
        common::release(m_vertexSparseShared2);
#endif
#ifdef PP_N_VERTEX_SPARSE_SHARED_3
        common::release(m_vertexSparseShared3);
#endif
#ifdef PP_N_VERTEX_SPARSE_SHARED_4
        common::release(m_vertexSparseShared4);
#endif
#ifdef PP_N_VERTEX_SPARSE_SHARED_5
        common::release(m_vertexSparseShared5);
#endif

       // faces
#ifdef PP_N_FACE_SPARSE_PRIVATE_1
        common::release(m_faceSparsePrivate1);
#endif
#ifdef PP_N_FACE_SPARSE_PRIVATE_2
        common::release(m_faceSparsePrivate2);
#endif
#ifdef PP_N_FACE_SPARSE_PRIVATE_3
        common::release(m_faceSparsePrivate3);
#endif
#ifdef PP_N_FACE_SPARSE_PRIVATE_4
        common::release(m_faceSparsePrivate4);
#endif
#ifdef PP_N_FACE_SPARSE_PRIVATE_5
        common::release(m_faceSparsePrivate5);
#endif

#ifdef PP_N_FACE_SPARSE_SHARED_1
        common::release(m_faceSparseShared1);
#endif
#ifdef PP_N_FACE_SPARSE_SHARED_2
        common::release(m_faceSparseShared2);
#endif
#ifdef PP_N_FACE_SPARSE_SHARED_3
        common::release(m_faceSparseShared3);
#endif
#ifdef PP_N_FACE_SPARSE_SHARED_4
        common::release(m_faceSparseShared4);
#endif
#ifdef PP_N_FACE_SPARSE_SHARED_5
        common::release(m_faceSparseShared5);
#endif
#ifdef PP_N_FACE_SPARSE_SHARED_6
        common::release(m_faceSparseShared6);
#endif
        // elements
#ifdef PP_N_ELEMENT_SPARSE_PRIVATE_1
        common::release(m_elementSparsePrivate1);
#endif
#ifdef PP_N_ELEMENT_SPARSE_PRIVATE_2
        common::release(m_elementSparsePrivate2);
#endif
#ifdef PP_N_ELEMENT_SPARSE_PRIVATE_3
        common::release(m_elementSparsePrivate3);
#endif
#ifdef PP_N_ELEMENT_SPARSE_PRIVATE_4
        common::release(m_elementSparsePrivate4);
#endif
#ifdef PP_N_ELEMENT_SPARSE_PRIVATE_5
        common::release(m_elementSparsePrivate5);
#endif

#ifdef PP_N_ELEMENT_SPARSE_SHARED_1
        common::release(m_elementSparseShared1);
#endif
#ifdef PP_N_ELEMENT_SPARSE_SHARED_2
        common::release(m_elementSparseShared2);
#endif
#ifdef PP_N_ELEMENT_SPARSE_SHARED_3
        common::release(m_elementSparseShared3);
#endif
#ifdef PP_N_ELEMENT_SPARSE_SHARED_4
        common::release(m_elementSparseShared4);
#endif
#ifdef PP_N_ELEMENT_SPARSE_SHARED_5
        common::release(m_elementSparseShared5);
#endif
#ifdef PP_N_ELEMENT_SPARSE_SHARED_6
        common::release(m_elementSparseShared6);
#endif
      }
    }

};

#endif
