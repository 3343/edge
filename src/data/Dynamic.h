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
 * Dynamic memory allocations, bookeeping, and memory release on destruction.
 **/

#ifndef DYNAMIC_HPP
#define DYNAMIC_HPP

#include "common.hpp"

namespace edge {
  namespace data {
    class Dynamic;
  }
}

class edge::data::Dynamic {
  private:
    //! allocated memory
    std::vector< void* > m_mem;
    //! memory type of allocated memory
    std::vector< bool  > m_hbw;

  public:
   /**
    * Destructor, which frees all allocated memory.
    **/
    ~Dynamic();

    /**
     * Allocates aligned memory of the given size.
     * The memory will be released automatically on destruction of the object.
     * -> Do not free manually.
     *
     * @param i_size size in bytes.
     * @param i_alignment alignment of the base pointer.
     * @param i_hbw if true, function allocates high bandwidth memory if available.
     * @return pointer to memory.
     **/
    void* allocate( std::size_t i_size,
                    std::size_t i_alignment=64,
                    bool        i_hbw=false );

    /**
     * Allocates EDGE's flex data type.
     *   The size of the memory allocation is based on sparse types.
     *   If the given bit-mask is set for an entity, the corresponding memory is allocated.
     *   If more than one bit-bask applies, the sum of required memory is allocated for that entity.
     *   The linear memory layout first stores the base data and then the data for the sparse types.
     *   An array of pointers for the data of base and every sparse type allows access based on entities.
     *
     * Example:
     *   i_nEn:                        5
     *   i_nSp:                        3
     *   i_bSize                       2
     *   i_spTypes:                   (0, 1, 0), (1, 1, 0), (0, 0, 1)
     *   i_spSizes:                    5,         3,         8
     *   .spType of i_enChars:        (0, 0, 1), (0, 1, 0), (1, 1, 0), (1, 0, 0), (1, 1, 1)
     *
     *                [-------base---------|------(0,1,0)------|------(1,1,0)-------|-------(0,0,1)------]
     *     addresses: [ b0, b1, b2, b3, b4 |  sp00, sp01, sp02 |   sp10,     sp11   |    sp20,    sp21   ]
     *     sizes:     [      5*2=10        |       3*5=15      |       2*3=6        |        2*8=16      ]
     *     pointers:         p00: b0             p10: nullptr        p20: nullptr          p30: sp20 
     *                       p01: b1             p11: sp00           p21: nullptr          p31: nullptr
     *                       p02: b2             p12: sp01           p22: sp10             p32: nullptr
     *                       p03: b3             p13: nullptr        p32: nullptr          p33: nullptr
     *                       p04: b4             p14: sp02           p33: sp11             p34: sp21
     *     total size: 10+15+6+16=47
     *
     * @param i_nEn number of entities.
     * @param i_nSp number of sparse types.
     * @param i_bSize base size per entry, which is always accounted for. 1 corresponds to sizeof(TL_T_FL).
     * @param i_spTypes array of sparse types.
     * @param i_spSizes associated sizes of the sparse types. 1 corresponds to sizeof(TL_T_FL).
     * @param i_enChars entity characteristics.
     * @param i_alignment alignment of the base pointer.
     * @param i_hbw if true, function allocates high bandwidth memory (if available).
     * @return array of pointers, which point to the respective start addresses of the entities. Last entry ( (i_nSp+1)*i_nEn+1 in total) is ghost for consistent size computations.
     *
     * @paramt TL_T_FL type of the items in the flex data structure.
     * @paramt TL_T_INT_LID integer type of local ids.
     * @paramt TL_T_INT_SP integer type of the sparse type.
     * @paramt TL_T_CHARS type of the entity characteristics. Offers a public member .spType for comparison with i_spTypes.
     **/
    template< typename TL_T_FL,
              typename TL_T_INT_LID,
              typename TL_T_INT_SP,
              typename TL_T_CHARS >
    TL_T_FL** flex( TL_T_INT_LID    i_nEn,
                    unsigned short  i_nSp,
                    std::size_t     i_bSize,
                    TL_T_INT_SP    *i_spTypes,
                    std::size_t    *i_spSizes,
                    TL_T_CHARS     *i_enChars,
                    std::size_t     i_alignment = 64,
                    bool            i_hbw       = false ) {
      // determine required amount of memory
      std::size_t l_memSize = 0;

      for( TL_T_INT_LID l_en = 0; l_en < i_nEn; l_en++ ) {
        l_memSize += i_bSize;

        for( unsigned short l_st = 0; l_st < i_nSp; l_st++ ) {
          if( (i_enChars[l_en].spType & i_spTypes[l_st]) == i_spTypes[l_st] ) {
            l_memSize += i_spSizes[l_st];
          }
        }
      }
      l_memSize *= sizeof( TL_T_FL );

      // allocate raw memory
      TL_T_FL *l_raw = (TL_T_FL*) allocate( l_memSize,
                                            i_alignment,
                                            i_hbw );

      // allocate memory for the pointers
      TL_T_FL** l_ptrs = (TL_T_FL **) allocate( ( (i_nSp+1)*i_nEn + 1 )*sizeof(TL_T_FL),
                                                i_alignment,
                                                i_hbw );

      // set the base pointers
      for( TL_T_INT_LID l_en = 0; l_en < i_nEn; l_en++ ) {
        l_ptrs[l_en] = l_raw;
        l_raw += i_bSize;
      }

      // set the sparse pointers
      for( unsigned short l_st = 0; l_st < i_nSp; l_st++ ) {
        for( TL_T_INT_LID l_en = 0; l_en < i_nEn; l_en++ ) {
          std::size_t l_ptrId = i_nEn*(l_st+1)+l_en;

          if( (i_enChars[l_en].spType & i_spTypes[l_st]) == i_spTypes[l_st] ) {
            l_ptrs[l_ptrId] = l_raw;
            l_raw += i_spSizes[l_st];
          }
          else {
            l_ptrs[l_ptrId] = nullptr;
          }
        }
      }

      // set ghost
      l_ptrs[ (i_nSp+1)*i_nEn ] = l_raw;

      return l_ptrs;
    }
};
#endif
