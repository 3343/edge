/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2020-2021, Friedrich Schiller University Jena
 * Copyright (c) 2019-2020, Alexander Breuer
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
 * Mesh functions and storage.
 **/
#include "Mesh.h"

#include "../geom/Generic.h"
#include "../geom/Geom.h"
#include "../io/logging.h"

void edge_v::mesh::Mesh::getElFaEl( t_entityType         i_elTy,
                                    t_idx                i_nEls,
                                    t_idx        const * i_faEl,
                                    t_idx        const * i_elFa,
                                    t_idx              * o_elFaEl ) {
  unsigned short l_nElFas = CE_N_FAS( i_elTy );

#ifdef PP_USE_OMP
#pragma omp parallel for
#endif
  for( t_idx l_el = 0; l_el < i_nEls; l_el++ ) {
    for( unsigned short l_fa = 0; l_fa < l_nElFas; l_fa++ ) {
      t_idx l_faId = i_elFa[l_el*l_nElFas + l_fa];

      if( i_faEl[l_faId*2 + 0] == l_el ) {
        o_elFaEl[ l_el*l_nElFas + l_fa ] =  i_faEl[l_faId*2 + 1];
      }
      else {
        EDGE_V_CHECK_EQ( i_faEl[l_faId*2 + 1], l_el );
        o_elFaEl[ l_el*l_nElFas + l_fa ] =  i_faEl[l_faId*2 + 0];
      }
    }
  }
}

edge_v::t_idx edge_v::mesh::Mesh::getAddEntry( t_idx         i_sizeFirst,
                                               t_idx         i_sizeSecond,
                                               t_idx const * i_first,
                                               t_idx const * i_second ) {
  for( t_idx l_se = 0; l_se < i_sizeSecond; l_se++ ) {
    bool l_present = false;

    for( t_idx l_fi = 0; l_fi < i_sizeFirst; l_fi++ ) {
      if( i_second[l_se] == i_first[l_fi] ) {
        l_present = true;
        break;
      }
    }

    if( l_present == false ) return i_second[l_se];
  }

  return std::numeric_limits< t_idx >::max();
}

void edge_v::mesh::Mesh::getEnVeCrds( t_entityType          i_enTy,
                                      t_idx        const  * i_enVe,
                                      double       const (* i_veCrds)[3],
                                      double             (* o_enVeCrds)[3] ) {
  // number of vertices for the entity type
  unsigned short l_nEnVes = CE_N_VES( i_enTy );

  // gather coordinates
  for( unsigned short l_ve = 0; l_ve < l_nEnVes; l_ve++ ) {
    t_idx l_veId = i_enVe[l_ve];

    for( unsigned short l_di = 0; l_di < 3; l_di++ ) {
      o_enVeCrds[l_ve][l_di] = i_veCrds[l_veId][l_di];
    }
  }
}

void edge_v::mesh::Mesh::addSparseTypeEn( t_entityType         i_enTy,
                                          t_idx                i_nEnsDe,
                                          t_idx                i_nEnsSp,
                                          t_idx        const * i_enVeDe,
                                          t_idx        const * i_enVeSp,
                                          t_sparseType         i_spTypeAdd,
                                          t_sparseType       * io_spType ) {
  // get the number of vertices
  unsigned short l_nEnVes = CE_N_VES( i_enTy );

  // current sparse id
  t_idx l_sp = 0;

  // iterate over the dense entities
  for( t_idx l_de = 0; l_de < i_nEnsDe; l_de++ ) {
    // abort if all sparse entities have been processed
    if( l_sp >= i_nEnsSp ) break;

    // check for number of matching vertices
    unsigned short l_nMatch = 0;
    for( unsigned short l_ve = 0; l_ve < l_nEnVes; l_ve++ ) {
      if( i_enVeDe[ l_de * l_nEnVes + l_ve] == i_enVeSp[ l_sp * l_nEnVes + l_ve ] ) {
        l_nMatch++;
      }
    }

    // set type and increase counter for matches
    if( l_nMatch == l_nEnVes ) {
      io_spType[l_de] |= i_spTypeAdd;
      l_sp++;
    }
  }

  // check that we found all sparse entities
  EDGE_V_CHECK_EQ( l_sp, i_nEnsSp );
}

void edge_v::mesh::Mesh::setInDiameter( t_entityType          i_enTy,
                                        t_idx                 i_nEns,
                                        t_idx        const  * i_enVe,
                                        double       const (* i_veCrds)[3],
                                        double              * o_inDia ) {
  // get the number of vertices
  unsigned short l_nEnVes = CE_N_VES( i_enTy );

  // check the size of the stack-memory below
  EDGE_V_CHECK_LE( l_nEnVes, 8 );

#ifdef PP_USE_OMP
#pragma omp parallel for
#endif
  for( t_idx l_en = 0; l_en < i_nEns; l_en++ ) {
    // get vertex coordinates
    double l_veCrds[8][3] = {};
    getEnVeCrds( i_enTy,
                 i_enVe + (l_nEnVes * l_en),
                 i_veCrds,
                 l_veCrds );

    // compute diameter
    o_inDia[l_en] = geom::Geom::inDiameter( i_enTy,
                                            l_veCrds );
  }
}

void edge_v::mesh::Mesh::setPeriodicBnds( t_entityType                  i_elTy,
                                          t_idx                         i_nFas,
                                          int                           i_peBndTy,
                                          int                  const  * i_faBndTys,
                                          t_idx                const  * i_faVe,
                                          double               const (* i_veCrds)[3],
                                          t_idx                       * io_faEl,
                                          t_idx                const  * i_elFa,
                                          t_idx                       * io_elFaEl,
                                          std::vector< t_idx >        & o_pFasGt ) {
  o_pFasGt.resize(0);

  // get boundary faces
  std::vector< t_idx > l_bndFas;
  for( t_idx l_fa = 0; l_fa < i_nFas; l_fa++ ) {
    if( i_faBndTys[l_fa] == i_peBndTy ) {
      EDGE_V_CHECK_EQ( io_faEl[l_fa*2+1], std::numeric_limits< t_idx >::max() );
      l_bndFas.push_back( l_fa );
    }
  }

  unsigned short l_nDis = CE_N_DIS( i_elTy );
  unsigned short l_nFaVes = CE_N_VES( CE_T_FA( i_elTy ) );

  // lambda which derives the dimension,
  // where the number of matching vertex coordinates matches the requested limit
  auto l_eqDi = [ l_nDis, l_nFaVes, i_faVe, i_veCrds ]( unsigned short i_limit,
                                                        t_idx          i_f0,
                                                        t_idx          i_f1 ) {
    for( unsigned short l_di = 0; l_di < l_nDis; l_di++) {
      unsigned short l_nEq = 0;

      for( unsigned short l_ve0 = 0; l_ve0 < l_nFaVes; l_ve0++ ) {
        for( unsigned short l_ve1 = 0; l_ve1 < l_nFaVes; l_ve1++ ) {
          t_idx l_ve0Id = i_faVe[i_f0*l_nFaVes + l_ve0];
          t_idx l_ve1Id = i_faVe[i_f1*l_nFaVes + l_ve1];

          double l_diff = i_veCrds[l_ve0Id][l_di] - i_veCrds[l_ve1Id][l_di];
          if( std::abs(l_diff) < m_tol ) l_nEq++;
        }
      }

      // if the number of vertices sharing coordinates in a dimension matches the request, we are done
      if( l_nEq == i_limit ) return l_di;
    }
    return std::numeric_limits< unsigned short >::max();
  };

  // lambda which derives the number of matching face-vertices assuming a fixed dimensions
  auto l_maVes = [ l_nDis, l_nFaVes, i_faVe, i_veCrds ]( unsigned short i_fixedDi,
                                                         t_idx          i_f0,
                                                         t_idx          i_f1 ) {
    unsigned short l_nMaVes = 0;

    for( unsigned short l_ve0 = 0; l_ve0 < l_nFaVes; l_ve0++ ) {
      for( unsigned short l_ve1 = 0; l_ve1 < l_nFaVes; l_ve1++ ) {
        t_idx l_ve0Id = i_faVe[i_f0*l_nFaVes + l_ve0];
        t_idx l_ve1Id = i_faVe[i_f1*l_nFaVes + l_ve1];

        unsigned short l_nMaDis = 0;
        for( unsigned short l_di = 0; l_di < l_nDis; l_di++) {
          double l_diff = i_veCrds[l_ve0Id][l_di] - i_veCrds[l_ve1Id][l_di];
          if(    l_di != i_fixedDi
              && std::abs(l_diff) < m_tol ) l_nMaDis++;
        }

        if( l_nMaDis == l_nDis-1 ) l_nMaVes++;
      }
    }

    return l_nMaVes;
  };

  // 1) iterate over all face-pairs,
  // 2) for each face, determine the constant dimension (assumption for our periodic boundaries)
  // 3) if the two faces have the same constant dim, check if the faces' vertices share coordinates in another dim
  std::vector< t_idx > l_faPairs;
  for( std::size_t l_f0 = 0; l_f0 < l_bndFas.size(); l_f0++ ) {
    unsigned short l_constDim0 = l_eqDi( l_nFaVes*l_nFaVes, l_bndFas[l_f0], l_bndFas[l_f0] );
    for( std::size_t l_f1 = 0; l_f1 < l_bndFas.size(); l_f1++ ) {
      unsigned short l_constDim1 = l_eqDi( l_nFaVes*l_nFaVes, l_bndFas[l_f1], l_bndFas[l_f1] );

      if( l_constDim0 == l_constDim1 ) {
        unsigned short l_nMaVes = l_maVes( l_constDim0, l_bndFas[l_f0], l_bndFas[l_f1] );

        if( l_f0 != l_f1 && l_nMaVes == l_nFaVes ) {
          l_faPairs.push_back( l_f1 );
        }
      }
    }

    // add dummy if nothing was found (happens for periodic partition boundaries)
    if( l_faPairs.size() != l_f0 + 1 ) {
      l_faPairs.push_back( std::numeric_limits< t_idx >::max() );
    }
  }

  // check that we got a partner for every found periodic face
  EDGE_V_CHECK_EQ( l_faPairs.size(), l_bndFas.size() );
  for( std::size_t l_fa = 0; l_fa < l_faPairs.size(); l_fa++ ) {
    if( l_faPairs[l_fa] != std::numeric_limits< t_idx >::max() ) {
      EDGE_V_CHECK_EQ( l_faPairs[ l_faPairs[l_fa] ], l_fa );
    }
  }

  // insert the missing elements and store the face with the larger adjacent element
  unsigned short l_nElFas = CE_N_FAS( i_elTy );
  for( std::size_t l_f0 = 0; l_f0 < l_faPairs.size(); l_f0++ ) {
    t_idx l_f1 = l_faPairs[l_f0];
    if( l_f1 == std::numeric_limits< t_idx >::max() ) continue;

    t_idx l_f0Id = l_bndFas[ l_f0 ];
    t_idx l_f1Id = l_bndFas[ l_f1 ];

    t_idx l_el0 = io_faEl[l_f0Id*2 + 0];
    t_idx l_el1 = io_faEl[l_f1Id*2 + 0];
    io_faEl[l_f0Id*2 + 1] = l_el1;

    // store faces where the first element has greater id than the second one
    if( l_el0 > l_el1 ) {
      o_pFasGt.push_back( l_f0Id );
    }

    for( unsigned short l_ad = 0; l_ad < l_nElFas; l_ad++ ) {
      if( i_elFa[l_el0*l_nElFas + l_ad] == l_f0Id ) {
        io_elFaEl[l_el0*l_nElFas + l_ad] = l_el1;
      }
    }
  }
}

void edge_v::mesh::Mesh::normOrder( t_entityType          i_elTy,
                                    t_idx                 i_nFas,
                                    t_idx                 i_nEls,
                                    double       const (* i_veCrds)[3],
                                    t_idx               * io_faVe,
                                    t_idx               * io_faEl,
                                    t_idx               * io_elVe,
                                    t_idx               * io_elFa,
                                    t_idx               * io_elFaEl ) {
  unsigned short l_nFaVes = CE_N_VES( CE_T_FA( i_elTy ) );
  unsigned short l_nElVes = CE_N_VES( i_elTy );
  unsigned short l_nElFas = CE_N_FAS( i_elTy );

#ifdef PP_USE_OMP
#pragma omp parallel for
#endif
  for( t_idx l_fa = 0; l_fa < i_nFas; l_fa++ ) {
    std::sort( io_faVe+l_fa*l_nFaVes, io_faVe+(l_fa+1)*l_nFaVes );
    std::sort( io_faEl+l_fa*2,        io_faEl+(l_fa+1)*2 );
  }

#ifdef PP_USE_OMP
#pragma omp parallel for
#endif
  for( t_idx l_el = 0; l_el < i_nEls; l_el++ ) {
    std::sort( io_elVe+l_el*l_nElVes, io_elVe+(l_el+1)*l_nElVes );

    // lambda, which compares two faces lexicographically
    auto l_faLess = [ io_faVe, l_nFaVes ]( t_idx i_faId0,
                                           t_idx i_faId1 ) {
      // get the vertex ids of the faces
      t_idx const * l_faVe0 = io_faVe + i_faId0*l_nFaVes;
      t_idx const * l_faVe1 = io_faVe + i_faId1*l_nFaVes;

      // iterate over the vertices and compare their ids
      for( unsigned short l_ve = 0; l_ve < l_nFaVes; l_ve++ ) {
        // abort only for non-equal
        if( l_faVe0[l_ve] < l_faVe1[l_ve] )
          return true;
        else if( l_faVe0[l_ve] > l_faVe1[l_ve] )
          return false;
      }

      // equal would be valid but not for the use-case below
      EDGE_V_LOG_FATAL;
      return false;
    };

    // sort the faces
    for( unsigned short l_f0 = 0; l_f0 < l_nElFas; l_f0++ ) {
      for( unsigned short l_f1 = l_f0+1; l_f1 < l_nElFas; l_f1++ ) {
        // get the ids of the two faces
        t_idx l_faId0 = io_elFa[l_el*l_nElFas + l_f0];
        t_idx l_faId1 = io_elFa[l_el*l_nElFas + l_f1];

        // check if the id of face at position f1 is lexicographically less than the one at f0
        bool l_less = l_faLess( l_faId1, l_faId0 );

        // move faces around if not in order
        if( l_less ) {
          io_elFa[l_el*l_nElFas + l_f0] = l_faId1;
          io_elFa[l_el*l_nElFas + l_f1] = l_faId0;

          t_idx l_elId0 = io_elFaEl[l_el*l_nElFas + l_f0];
          t_idx l_elId1 = io_elFaEl[l_el*l_nElFas + l_f1];

          io_elFaEl[l_el*l_nElFas + l_f0] = l_elId1;
          io_elFaEl[l_el*l_nElFas + l_f1] = l_elId0;
        }
      }
    }
  }

  for( t_idx l_el = 0; l_el < i_nEls; l_el++ ) {
    // get vertex coordinates
    double l_veCrds[8][3] = {};
    getEnVeCrds( i_elTy,
                 io_elVe + (l_nElVes * l_el),
                 i_veCrds,
                 l_veCrds );

    geom::Geom::normVesFas( i_elTy,
                            l_veCrds,
                            io_elVe+l_el*l_nElVes,
                            io_elFa+l_el*l_nElFas,
                            io_elFaEl+l_el*l_nElFas );
  }
}

edge_v::mesh::Mesh::Mesh( edge_v::io::Gmsh const & i_gmsh,
                          int                      i_periodic ) {
  // get the element type of the mesh
  m_elTy = i_gmsh.getElType();
  t_entityType l_faTy = CE_T_FA( m_elTy );

  // derive element properties
  unsigned short l_nElVes = CE_N_VES( m_elTy );
  unsigned short l_nElFas = CE_N_FAS( m_elTy );
  unsigned short l_nFaVes = CE_N_VES( l_faTy );

  // query mesh
  m_nVes = i_gmsh.nVes();
  EDGE_V_CHECK_GT( m_nVes, 0 );
  m_nEls = i_gmsh.nEls();
  EDGE_V_CHECK_GT( m_nEls, 0 );

  std::vector< t_idx > l_elVe;
  l_elVe.resize( nEls()*l_nElVes );
  i_gmsh.getElVe( l_elVe.data() );

  std::vector< t_idx > l_elFaVe;
  l_elFaVe.resize( m_nEls*l_nElFas*l_nFaVes );
  i_gmsh.getElFaVe( l_elFaVe.data() );

  // convert faces to custom struct
  struct Face {
    t_idx el;
    unsigned short fa;
    t_idx ves[4] = { std::numeric_limits< t_idx >::max(),
                     std::numeric_limits< t_idx >::max(),
                     std::numeric_limits< t_idx >::max(),
                     std::numeric_limits< t_idx >::max() };
  };
  std::vector< Face > l_mapFas;
  l_mapFas.resize( m_nEls*l_nElFas );
  for( std::size_t l_el = 0; l_el < m_nEls; l_el++ ) {
    for( unsigned short l_fa = 0; l_fa < l_nElFas; l_fa++ ) {
      std::size_t l_mapId = l_el*l_nElFas + l_fa;
      l_mapFas[l_mapId].el = l_el;
      l_mapFas[l_mapId].fa = l_fa;
      for( unsigned short l_ve = 0; l_ve < l_nFaVes; l_ve++ ) {
        l_mapFas[l_mapId].ves[l_ve] = l_elFaVe[l_mapId*l_nFaVes + l_ve];
      }
      std::sort( l_mapFas[l_mapId].ves,
                 l_mapFas[l_mapId].ves+4 );
    }
  }

  // sort by vertices
  std::sort( l_mapFas.begin(),
             l_mapFas.end(),
             []( Face const & l_f0, Face const & l_f1 ) -> bool {
               return   std::tie( l_f0.ves[0], l_f0.ves[1], l_f0.ves[2], l_f0.ves[3] )
                      < std::tie( l_f1.ves[0], l_f1.ves[1], l_f1.ves[2], l_f1.ves[3] ); } );

  // derive number of internal faces
  t_idx l_nFasInt = 0;
  for( std::size_t l_fa = 0; l_fa < l_mapFas.size()-1; l_fa++ ) {
    unsigned short l_nVesShared = 0;
    for( unsigned short l_ve = 0; l_ve < l_nFaVes; l_ve++ ) {
      if( l_mapFas[l_fa].ves[l_ve] == l_mapFas[l_fa+1].ves[l_ve] ) {
        l_nVesShared++;
      }
    }
    
    if( l_nVesShared == l_nFaVes ) {
      l_nFasInt++;
    }
  }
  t_idx l_nFasBnd = l_mapFas.size() - l_nFasInt*2;
  m_nFas = l_nFasInt + l_nFasBnd;

  // allocate memory
  t_idx l_size  = m_nVes * 3;
  m_veCrds = (double (*)[3]) new double[ l_size ];

  l_size = m_nFas * l_nFaVes;
  m_faVe = new t_idx[ l_size ];

  l_size = m_nFas * 2;
  m_faEl = new t_idx[ l_size ];

  l_size  = m_nEls * l_nElVes;
  m_elVe = new t_idx[ l_size ];

  l_size = m_nEls * l_nElFas;
  m_elFa = new t_idx[ l_size ];

  l_size = m_nVes;
  m_spTypeVe = new t_sparseType[ l_size ];

  l_size = m_nFas;
  m_spTypeFa = new t_sparseType[ l_size ];

  l_size = m_nEls;
  m_spTypeEl = new t_sparseType[ l_size ];

  l_size  = m_nEls;
  m_inDiasEl = new double[ l_size ];

  l_size = m_nEls * l_nElFas;
  m_elFaEl = new t_idx[ l_size ];

  // init data
  for( t_idx l_ve = 0; l_ve < m_nVes; l_ve++ ) {
    for( unsigned short l_di = 0; l_di < 3; l_di++ ) {
      m_veCrds[l_ve][l_di] = std::numeric_limits< double >::max();
    }
    m_spTypeVe[l_ve] = 0;
  }
  for( t_idx l_fa = 0; l_fa < m_nFas; l_fa++ ) {
    for( unsigned short l_ve = 0; l_ve < l_nFaVes; l_ve++ ) {
      m_faVe[l_fa*l_nFaVes + l_ve] = std::numeric_limits< t_idx >::max();
    }
    for( unsigned short l_sd = 0; l_sd < 2; l_sd++ ) {
      m_faEl[l_fa*2 + l_sd] = std::numeric_limits< t_idx >::max();
    }
    m_spTypeFa[l_fa] = 0;
  }
  for( t_idx l_el = 0; l_el < m_nEls; l_el++ ) {
    for( unsigned short l_ve = 0; l_ve < l_nElVes; l_ve++ ) {
      m_elVe[l_el*l_nElVes + l_ve] = std::numeric_limits< t_idx >::max();
    }
    for( unsigned short l_fa = 0; l_fa < l_nElFas; l_fa++ ) {
      m_elFa[l_el*l_nElFas + l_fa] = std::numeric_limits< t_idx >::max();
      m_elFaEl[l_el*l_nElFas + l_fa] = std::numeric_limits< t_idx >::max();
    }
    m_spTypeEl[l_el] = 0;
  }

  // query gmsh for 1-to-1 data structures
  i_gmsh.getElVe( m_elVe );
  i_gmsh.getVeCrds( m_veCrds );

  // assign face related data
  std::size_t l_mapId = 0;
  for( t_idx l_fa = 0; l_fa < m_nFas; l_fa++ ) {
    // check if this is an internal face
    bool l_isInt = false;
    if( l_mapId < l_mapFas.size()-1 ) {
      unsigned short l_nVesShared = 0;
      for( unsigned short l_ve = 0; l_ve < l_nFaVes; l_ve++ ) {
        if( l_mapFas[l_mapId].ves[l_ve] == l_mapFas[l_mapId+1].ves[l_ve] ) {
          l_nVesShared++;
        }
      }
      if( l_nVesShared == l_nFaVes ) l_isInt = true;
    }

    // assign data accordingly
    for( unsigned short l_ve = 0; l_ve < l_nFaVes; l_ve++ ) {
      m_faVe[l_fa*l_nFaVes + l_ve] = l_mapFas[l_mapId].ves[l_ve];
    }

    m_faEl[l_fa*2 + 0] = l_mapFas[l_mapId].el;
    if( l_isInt ) {
      m_faEl[l_fa*2 + 1] = l_mapFas[l_mapId+1].el;
    }

    m_elFa[l_mapFas[l_mapId].el*l_nElFas + l_mapFas[l_mapId].fa] = l_fa;
    if( l_isInt ) {
      m_elFa[l_mapFas[l_mapId+1].el*l_nElFas + l_mapFas[l_mapId+1].fa] = l_fa;
    }

    if( l_isInt ) {
      m_elFaEl[l_mapFas[l_mapId  ].el*l_nElFas + l_mapFas[l_mapId  ].fa] = l_mapFas[l_mapId+1].el;
      m_elFaEl[l_mapFas[l_mapId+1].el*l_nElFas + l_mapFas[l_mapId+1].fa] = l_mapFas[l_mapId  ].el;
    }

    l_mapId++;
    if( l_isInt ) l_mapId++;
  }

  // derive sparse types
  t_idx l_nPhGrs = i_gmsh.nPhysicalGroupsFa();
  int const * l_phGrs = i_gmsh.getPhysicalGroupsFa();

  for( t_idx l_gr = 0; l_gr < l_nPhGrs; l_gr++ ) {
    t_idx l_nFasPhGr = i_gmsh.nFas( l_phGrs[l_gr] );

    // get faces' vertices of the group
    l_size = l_nFasPhGr * l_nFaVes;
    t_idx * l_faVePhGr = new t_idx[ l_size ];
    i_gmsh.getFaVe( l_phGrs[l_gr],
                    l_faVePhGr );

    // sort vertices
    for( t_idx l_fa = 0; l_fa < l_nFasPhGr; l_fa++ ) {
      std::sort( l_faVePhGr+l_fa*l_nFaVes, l_faVePhGr+(l_fa+1)*l_nFaVes );
    }

    // sort lexicographically
    geom::Generic::sortLex( l_nFasPhGr,
                            l_nFaVes,
                            l_faVePhGr );

    // add physical group to sparse type
    addSparseTypeEn( l_faTy,
                     m_nFas,
                     l_nFasPhGr,
                     m_faVe,
                     l_faVePhGr,
                     l_phGrs[l_gr],
                     m_spTypeFa );


    // free memory
    delete[] l_faVePhGr;
  }

  // perform sanity checks
  for( t_idx l_ve = 0; l_ve < m_nVes; l_ve++ ) {
    for( unsigned short l_di = 0; l_di < 3; l_di++ ) {
      EDGE_V_CHECK_NE( m_veCrds[l_ve][l_di], std::numeric_limits< double >::max() );
    }
  }
  for( t_idx l_fa = 0; l_fa < m_nFas; l_fa++ ) {
    for( unsigned short l_ve = 0; l_ve < l_nFaVes; l_ve++ ) {
      EDGE_V_CHECK_NE( m_faVe[l_fa*l_nFaVes + l_ve], std::numeric_limits< t_idx >::max() );
    }
    EDGE_V_CHECK_NE( m_faEl[l_fa*2 + 0], std::numeric_limits< t_idx >::max() );
  }
  for( t_idx l_el = 0; l_el < m_nEls; l_el++ ) {
    for( unsigned short l_ve = 0; l_ve < l_nElVes; l_ve++ ) {
      EDGE_V_CHECK_NE( m_elVe[l_el*l_nElVes + l_ve], std::numeric_limits< t_idx >::max() );
    }
    for( unsigned short l_fa = 0; l_fa < l_nElFas; l_fa++ ) {
      EDGE_V_CHECK_NE( m_elFa[l_el*l_nElFas + l_fa], std::numeric_limits< t_idx >::max() );
    }
  }

  // compute mesh properties
  setInDiameter( m_elTy,
                 m_nEls,
                 m_elVe,
                 m_veCrds,
                 m_inDiasEl );

  // adjust for periodic boundaries
  std::vector< t_idx > l_pFasGt;
  if( i_periodic != std::numeric_limits< int >::max() ) {
    setPeriodicBnds( m_elTy,
                     m_nFas,
                     i_periodic,
                     m_spTypeFa,
                     m_faVe,
                     m_veCrds,
                     m_faEl,
                     m_elFa,
                     m_elFaEl,
                     l_pFasGt );
  }

  normOrder( m_elTy,
             m_nFas,
             m_nEls,
             m_veCrds,
             m_faVe,
             m_faEl,
             m_elVe,
             m_elFa,
             m_elFaEl );

  // reverse the order of larger-element periodic faces for consistent normals
  for( std::size_t l_pf = 0; l_pf < l_pFasGt.size(); l_pf++ ) {
    t_idx l_fa = l_pFasGt[l_pf];

    t_idx l_tmpEl = m_faEl[l_fa*2+0];
    m_faEl[l_fa*2+0] = m_faEl[l_fa*2+1];
    m_faEl[l_fa*2+1] = l_tmpEl;
  }
}

edge_v::mesh::Mesh::~Mesh() {
  delete[] m_veCrds;
  delete[] m_faVe;
  delete[] m_faEl;
  delete[] m_elVe;
  delete[] m_elFa;
  delete[] m_elFaEl;
  delete[] m_spTypeVe;
  delete[] m_spTypeFa;
  delete[] m_spTypeEl;
  delete[] m_inDiasEl;
  if( m_volFa != nullptr    ) delete[] m_volFa;
  if( m_volEl != nullptr    ) delete[] m_volEl;
  if( m_cenEl != nullptr    ) delete[] m_cenEl;
  if( m_normals != nullptr  ) delete[] m_normals;
  if( m_tangents != nullptr ) delete[] m_tangents;
}

void edge_v::mesh::Mesh::printStats() const {
  EDGE_V_LOG_INFO << "mesh stats:";
  EDGE_V_LOG_INFO << "  #vertices: " << m_nVes;
  EDGE_V_LOG_INFO << "  #faces:    " << m_nFas;
  EDGE_V_LOG_INFO << "  #elements: " << m_nEls;
}

double const (* edge_v::mesh::Mesh::getCentroidsEl())[3] {
  // only work on this once
  if( m_cenEl == nullptr ) {
    // allocate memory
    t_idx l_size  = m_nEls * 3;
    m_cenEl = (double (*)[3]) new double[ l_size ];

    unsigned short l_nElVes = CE_N_VES( m_elTy );

#ifdef PP_USE_OMP
#pragma omp parallel for
#endif
    for( t_idx l_el = 0; l_el < m_nEls; l_el++ ) {
      // get vertex coordinates
      double l_veCrds[8][3] = {};
      getEnVeCrds( m_elTy,
                   m_elVe + (l_nElVes * l_el),
                   m_veCrds,
                   l_veCrds );

      // compute centroid
      geom::Geom::centroid( m_elTy,
                            l_veCrds,
                            m_cenEl[l_el] );
    }
  }

  return m_cenEl;
}

double const * edge_v::mesh::Mesh::getAreasFa() {
  // only work on this once
  if( m_volFa == nullptr ) {
    // allocate memory
    m_volFa = new double[ m_nFas ];

    t_entityType l_faTy = CE_T_FA( m_elTy );
    unsigned short l_nFaVes = CE_N_VES( l_faTy );

#ifdef PP_USE_OMP
#pragma omp parallel for
#endif
    for( t_idx l_fa = 0; l_fa < m_nFas; l_fa++ ) {
      // get vertex coordinates
      double l_veCrds[4][3] = {};
      getEnVeCrds( l_faTy,
                   m_faVe + (l_nFaVes * l_fa),
                   m_veCrds,
                   l_veCrds );

      // compute "volume"
      m_volFa[l_fa] = geom::Geom::volume( l_faTy,
                                          l_veCrds );
    }
  }

  return m_volFa;
}

double const * edge_v::mesh::Mesh::getVolumesEl() {
  // only work on this once
  if( m_volEl == nullptr ) {
    // allocate memory
    m_volEl = new double[ m_nEls ];

    unsigned short l_nElVes = CE_N_VES( m_elTy );

#ifdef PP_USE_OMP
#pragma omp parallel for
#endif
    for( t_idx l_el = 0; l_el < m_nEls; l_el++ ) {
      // get vertex coordinates
      double l_veCrds[8][3] = {};
      getEnVeCrds( m_elTy,
                   m_elVe + (l_nElVes * l_el),
                   m_veCrds,
                   l_veCrds );

      // compute volume
      m_volEl[l_el] = geom::Geom::volume( m_elTy,
                                          l_veCrds );
    }
  }

  return m_volEl;
}

double const (* edge_v::mesh::Mesh::getNormalsFa() )[3] {
  // only work on this once
  if( m_normals == nullptr ) {
    // allocate memory
    m_normals = (double (*)[3]) new double[ m_nFas*3 ];

    t_entityType l_faTy = CE_T_FA( m_elTy );
    unsigned short l_nElVes = CE_N_VES( m_elTy );
    unsigned short l_nFaVes = CE_N_VES( l_faTy );

#ifdef PP_USE_OMP
#pragma omp parallel for
#endif
    for( t_idx l_fa = 0; l_fa < m_nFas; l_fa++ ) {
      // get vertex coordinates
      double l_veCrds[4][3] = {};
      getEnVeCrds( l_faTy,
                   m_faVe + (l_nFaVes * l_fa),
                   m_veCrds,
                   l_veCrds );

      // normal point from the first element
      t_idx l_el = m_faEl[l_fa*2];
      t_idx l_np = getAddEntry( l_nFaVes,
                                l_nElVes,
                                m_faVe+(l_nFaVes*l_fa),
                                m_elVe+(l_nElVes*l_el) );
      EDGE_V_CHECK_NE( l_np, std::numeric_limits< t_idx >::max() );

      // compute normals
      geom::Geom::normal( l_faTy,
                          l_veCrds,
                          m_veCrds[l_np],
                          m_normals[l_fa] );
    }
  }

  return m_normals;
}

double const (* edge_v::mesh::Mesh::getTangentsFa() )[2][3] {
  // only work on this once
  if( m_tangents == nullptr ) {
    // allocate memory
    m_tangents = (double (*)[2][3]) new double[ m_nFas*2*3 ];

    t_entityType l_faTy = CE_T_FA( m_elTy );
    unsigned short l_nElVes = CE_N_VES( m_elTy );
    unsigned short l_nFaVes = CE_N_VES( l_faTy );

#ifdef PP_USE_OMP
#pragma omp parallel for
#endif
    for( t_idx l_fa = 0; l_fa < m_nFas; l_fa++ ) {
      // get vertex coordinates
      double l_veCrds[4][3] = {};
      getEnVeCrds( l_faTy,
                   m_faVe + (l_nFaVes * l_fa),
                   m_veCrds,
                   l_veCrds );

      // normal point from the first element
      t_idx l_el = m_faEl[l_fa*2];
      t_idx l_np = getAddEntry( l_nFaVes,
                                l_nElVes,
                                m_faVe+(l_nFaVes*l_fa),
                                m_elVe+(l_nElVes*l_el) );
      EDGE_V_CHECK_NE( l_np, std::numeric_limits< t_idx >::max() );

      // compute tangents
      geom::Geom::tangents( l_faTy,
                            l_veCrds,
                            m_veCrds[l_np],
                            m_tangents[l_fa] );
    }
  }

  return m_tangents;
}

void edge_v::mesh::Mesh::getFaIdsAd( t_idx                  i_nFas,
                                     t_idx                  i_elOff,
                                     t_idx          const * i_el,
                                     unsigned short const * i_fa,
                                     unsigned short       * o_faIdsAd ) const {
  edge_v::geom::Geom::getFaIdsAd( m_elTy,
                                  i_nFas,
                                  i_elOff,
                                  i_el,
                                  i_fa,
                                  m_elFaEl,
                                  o_faIdsAd );
}


void edge_v::mesh::Mesh::getVeIdsAd( t_idx                  i_nFas,
                                     t_idx                  i_elOff,
                                     t_idx          const * i_el,
                                     unsigned short const * i_fa,
                                     unsigned short       * o_veIdsAd ) const {
  edge_v::geom::Geom::getVeIdsAd( m_elTy,
                                  i_nFas,
                                  i_elOff,
                                  i_el,
                                  i_fa,
                                  m_elVe,
                                  m_elFaEl,
                                  m_veCrds,
                                  o_veIdsAd );
}