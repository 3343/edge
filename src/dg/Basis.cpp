/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2016-2017, Regents of the University of California
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
 * Basis of an element.
 **/

#include "Basis.h"
#include "linalg/Mappings.hpp"
#include "QuadraturePoints.h"
#include <cassert>
#include <cmath>
#include <io/logging.h>
#include "data/common.hpp"
#include "monitor/instrument.hpp"

namespace edge {
  namespace pre {
    namespace dg {
      // mass matrix
      extern double const * g_massRaw;
      extern std::size_t const g_massSize;

      // stiffness matrices, premultiplied by the inverse mass matrix
      extern double const * g_stiffVRaw;
      extern std::size_t const g_stiffVSize;

      // transposed stiffness matrices, premultiplied by the inverse mass matrix (after transpose)
      extern double const * g_stiffTRaw;
      extern std::size_t const g_stiffTSize;
    }
  }
}

void edge::dg::Basis::initMassMatrix() {
  // check that the size matches
  EDGE_CHECK_EQ( pre::dg::g_massSize, m_nBaseFuncs*m_nBaseFuncs );

  real_base *l_mass;
  l_mass = new real_base[ pre::dg::g_massSize ];

  // set values
  for( std::size_t l_va = 0; l_va < pre::dg::g_massSize; l_va++ )
    l_mass[l_va] = pre::dg::g_massRaw[l_va];

  // convert to coordinate format
  linalg::Matrix::denseToCrd( m_nBaseFuncs,
                              m_nBaseFuncs,
                              l_mass,
                              m_mass );

  delete[] l_mass;
}

void edge::dg::Basis::computeFluxMatricesLine() {
  // determine number of flux matrices
  unsigned int l_nFluxMats = C_ENT[m_entType].N_FACES + C_ENT[m_entType].N_FACES * C_ENT[m_entType].N_FACES;

  m_flux.resize(l_nFluxMats);

  // iterate over faces
  for( int_md l_fa1 = 0; l_fa1 < C_ENT[m_entType].N_FACES; l_fa1++ ) {
    // local quadrature points
    std::vector< real_mesh > l_localQpsXi, l_localQpsEta, l_localQpsZeta, l_weights;

    l_localQpsXi.push_back(   C_REF_ELEMENT.VE.LINE[0][l_fa1] );
    l_localQpsEta.push_back(  C_REF_ELEMENT.VE.LINE[1][l_fa1] );
    l_localQpsZeta.push_back( C_REF_ELEMENT.VE.LINE[2][l_fa1] );
    l_weights.push_back(1);

    // compute local flux matrix
    for( unsigned int l_ro = 0; l_ro < m_nBaseFuncs; l_ro++ ) {
      for( unsigned int l_co = 0; l_co < m_nBaseFuncs; l_co++ ) {
        // value of the resp flux matrix
        real_base l_val = 0;

        for( unsigned int l_qp = 0; l_qp < l_localQpsXi.size(); l_qp++ ) {
          // evaluate basis at quad point
          real_base l_bValRo = 0;
          real_base l_bValCo = 0;

          evalBasisLine( l_ro, l_localQpsXi[l_qp], l_bValRo );
          evalBasisLine( l_co, l_localQpsXi[l_qp], l_bValCo );

          l_val += ( l_bValRo * l_bValCo ) * l_weights[l_qp];
        }

        if( std::abs( l_val ) > TOL.BASIS ) {
          m_flux[l_fa1].nz.push_back(l_val);
          m_flux[l_fa1].co.push_back(l_co);
          m_flux[l_fa1].ro.push_back(l_ro);
        }
      }
    }

    // compute neighboring flux matrices
    for( unsigned int l_fa2 = 0; l_fa2 < C_ENT[m_entType].N_FACES; l_fa2++ ) {
      std::vector< real_mesh > l_neighQpsXi, l_neighQpsEta, l_neighQpsZeta;

      l_neighQpsXi.push_back(   C_REF_ELEMENT.VE.LINE[0][l_fa2] );
      l_neighQpsEta.push_back(  C_REF_ELEMENT.VE.LINE[1][l_fa2] );
      l_neighQpsZeta.push_back( C_REF_ELEMENT.VE.LINE[2][l_fa2] );

      // map quad points to ref. element's face of adjacent element
      for( unsigned int l_qp = 0; l_qp < l_localQpsXi.size(); l_qp++ ) {
        linalg::Mappings::mapFaceCoords( LINE,
                                         l_fa1,
                                         l_fa2,
                                         0,
                                         l_localQpsXi[l_qp],
                                         l_localQpsEta[l_qp],
                                         l_localQpsZeta[l_qp],
                                         l_neighQpsXi[l_qp],
                                         l_neighQpsEta[l_qp],
                                         l_neighQpsZeta[l_qp] );
      }

      for( unsigned int l_ro = 0; l_ro < m_nBaseFuncs; l_ro++ ) {
        for( unsigned int l_co = 0; l_co < m_nBaseFuncs; l_co++ ) {
          // value of the flux matrix respective
          real_base l_val = 0;

          for( unsigned int l_qp = 0; l_qp < l_localQpsXi.size(); l_qp++ ) {
            // evaluate basis at quad point
            real_base l_bValRo = 0;
            real_base l_bValCo = 0;

            evalBasisLine( l_ro, l_neighQpsXi[l_qp], l_bValRo );
            evalBasisLine( l_co, l_localQpsXi[l_qp], l_bValCo );

            l_val += ( l_bValRo * l_bValCo ) * l_weights[l_qp];
          }
          if( std::abs( l_val ) > TOL.BASIS ) {
            m_flux[C_ENT[m_entType].N_FACES + l_fa1*C_ENT[m_entType].N_FACES+l_fa2].nz.push_back(l_val);
            m_flux[C_ENT[m_entType].N_FACES + l_fa1*C_ENT[m_entType].N_FACES+l_fa2].co.push_back(l_co);
            m_flux[C_ENT[m_entType].N_FACES + l_fa1*C_ENT[m_entType].N_FACES+l_fa2].ro.push_back(l_ro);
          }
        }
      }
    }

  }
}

void edge::dg::Basis::computeFluxMatrices2d() {
  // determine number of flux matrices
  // 1) local and neighboring faces
  unsigned int l_nFluxMats = C_ENT[m_entType].N_FACES * C_ENT[m_entType].N_FACES;

  // 2) local flux matrices stored upfront
               l_nFluxMats += C_ENT[m_entType].N_FACES;

  m_flux.resize(l_nFluxMats);

  // iterate over faces
  for( unsigned int l_fa1 = 0; l_fa1 < C_ENT[m_entType].N_FACES; l_fa1++ ) {
    // local quadrature points
    std::vector< real_mesh > l_localQpsXi, l_localQpsEta, l_localQpsZeta, l_weights;

    // get nodes of this face
    real_mesh l_faVes[3][2];
    for( unsigned int l_dim = 0; l_dim < 3; l_dim++ ) {
      l_faVes[l_dim][0] = C_REF_ELEMENT.FA.ENT[m_entType][
                             l_dim * C_ENT[m_entType].N_FACES * C_ENT[m_entType].N_FACE_VERTICES
                           + l_fa1                            * C_ENT[m_entType].N_FACE_VERTICES
                           + 0 ];
      l_faVes[l_dim][1] = C_REF_ELEMENT.FA.ENT[m_entType][
                             l_dim * C_ENT[m_entType].N_FACES * C_ENT[m_entType].N_FACE_VERTICES
                           + l_fa1                            * C_ENT[m_entType].N_FACE_VERTICES
                           + 1 ];
    }

    // get quad points of this face
    dg::QuadraturePoints::getQpts( LINE,
                                   m_order,
                                   l_faVes,
                                   l_localQpsXi,
                                   l_localQpsEta,
                                   l_localQpsZeta,
                                   l_weights );

    // sum the weights
    real_mesh l_sum = 0;
    for( unsigned int l_qp = 0; l_qp < l_weights.size(); l_qp++ ) l_sum += l_weights[l_qp];

    // normalize to face length to 1, this is what we are assuming in the mappings
    for( unsigned int l_qp = 0; l_qp < l_weights.size(); l_qp++ ) l_weights[l_qp] /= l_sum;

    // compute local flux matrices
    for( unsigned int l_ro = 0; l_ro < m_nBaseFuncs; l_ro++ ) {
      for( unsigned int l_co = 0; l_co < m_nBaseFuncs; l_co++ ) {
        // value of the flux matrix respective
        real_base l_val = 0;

        for( unsigned int l_qp = 0; l_qp < l_localQpsXi.size(); l_qp++ ) {
          // evaluate basis at quad point
          real_base l_bValRo = 0;
          real_base l_bValCo = 0;

          evalBasis( l_ro,
                     m_entType,
                     l_bValRo,
                     l_localQpsXi[l_qp],
                     l_localQpsEta[l_qp],
                     0, // zeta
                     -1, m_order );

          evalBasis( l_co,
                     m_entType,
                     l_bValCo,
                     l_localQpsXi[l_qp],
                     l_localQpsEta[l_qp],
                     0, // zeta
                     -1, m_order );

          l_val += ( l_bValRo * l_bValCo ) * l_weights[l_qp];
        }
        if( std::abs( l_val ) > TOL.BASIS ) {
          m_flux[l_fa1].nz.push_back(l_val);
          m_flux[l_fa1].co.push_back(l_co);
          m_flux[l_fa1].ro.push_back(l_ro);
        }
      }
    }

    // compute neighboring flux matrices
    for( unsigned int l_fa2 = 0; l_fa2 < C_ENT[m_entType].N_FACES; l_fa2++ ) {
      std::vector< real_mesh > l_neighQpsXi, l_neighQpsEta, l_neighQpsZeta;

      l_neighQpsXi.resize(   m_nBaseFuncs );
      l_neighQpsEta.resize(  m_nBaseFuncs );
      l_neighQpsZeta.resize( m_nBaseFuncs );

      // map quad points to ref. element's face of adjacent element
      for( unsigned int l_qp = 0; l_qp < l_localQpsXi.size(); l_qp++ ) {
        linalg::Mappings::mapFaceCoords( m_entType,
                                         l_fa1,
                                         l_fa2,
                                         0,
                                         l_localQpsXi[l_qp],
                                         l_localQpsEta[l_qp],
                                         l_localQpsZeta[l_qp],
                                         l_neighQpsXi[l_qp],
                                         l_neighQpsEta[l_qp],
                                         l_neighQpsZeta[l_qp] );
      }

      for( unsigned int l_ro = 0; l_ro < m_nBaseFuncs; l_ro++ ) {
        for( unsigned int l_co = 0; l_co < m_nBaseFuncs; l_co++ ) {
          // value in the respective flux matrix
          real_base l_val = 0;

          for( unsigned int l_qp = 0; l_qp < l_localQpsXi.size(); l_qp++ ) {
            // evaluate basis at quad point
            real_base l_bValRo = 0;
            real_base l_bValCo = 0;

            evalBasis( l_ro,
                       m_entType,
                       l_bValRo,
                       l_neighQpsXi[l_qp],
                       l_neighQpsEta[l_qp],
                       0, // zeta
                       -1, m_order );

            evalBasis( l_co,
                       m_entType,
                       l_bValCo,
                       l_localQpsXi[l_qp],
                       l_localQpsEta[l_qp],
                       0, // zeta
                       -1, m_order );

            l_val += ( l_bValRo * l_bValCo ) * l_weights[l_qp];
          }
          if( std::abs( l_val ) > TOL.BASIS ) {
            m_flux[C_ENT[m_entType].N_FACES + l_fa1*C_ENT[m_entType].N_FACES+l_fa2].nz.push_back(l_val);
            m_flux[C_ENT[m_entType].N_FACES + l_fa1*C_ENT[m_entType].N_FACES+l_fa2].co.push_back(l_co);
            m_flux[C_ENT[m_entType].N_FACES + l_fa1*C_ENT[m_entType].N_FACES+l_fa2].ro.push_back(l_ro);
          }
        }
      }
    }

  }

}

void edge::dg::Basis::computeFluxMatricesHex() {
  assert( C_ENT[m_entType].N_FACES == 6 );

  // fixed size of 12 flux matrices
  m_flux.resize( 12 );

  // iterate ofver faces
  for( unsigned int l_fa1 = 0; l_fa1 < 6; l_fa1++ ) {
    // local quadrature points
    std::vector< real_mesh > l_localQpsXi, l_localQpsEta, l_localQpsZeta, l_weights;

    // get nodes of this face
    real_mesh l_faVes[3][4];
    for( unsigned int l_dim = 0; l_dim < 3; l_dim++ ) {
      for( unsigned int l_ve = 0; l_ve < 4; l_ve++ ) {
        l_faVes[l_dim][l_ve] = C_REF_ELEMENT.FA.HEX[l_dim][l_fa1][l_ve];
      }
    }

    // get quad points of this face
    dg::QuadraturePoints::getQpts( QUAD4R,
                                   m_order,
                                   l_faVes,
                                   l_localQpsXi,
                                   l_localQpsEta,
                                   l_localQpsZeta,
                                   l_weights );

    // sum the weights
    real_mesh l_sum = 0;
    for( unsigned int l_qp = 0; l_qp < l_weights.size(); l_qp++ ) l_sum += l_weights[l_qp];

    // normalize to face length to 1, this is what we are assuming in the mappings
    for( unsigned int l_qp = 0; l_qp < l_weights.size(); l_qp++ ) l_weights[l_qp] /= l_sum;

    // compute local flux matrices
    for( unsigned int l_ro = 0; l_ro < m_nBaseFuncs; l_ro++ ) {
      for( unsigned int l_co = 0; l_co < m_nBaseFuncs; l_co++ ) {
        // value of the flux matrix respective
        real_base l_val = 0;

        for( unsigned int l_qp = 0; l_qp < l_localQpsXi.size(); l_qp++ ) {
          // evaluate basis at quad point
          real_base l_bValRo = 0;
          real_base l_bValCo = 0;

          evalBasis( l_ro,
                     m_entType,
                     l_bValRo,
                     l_localQpsXi[l_qp],
                     l_localQpsEta[l_qp],
                     l_localQpsZeta[l_qp],
                     -1, m_order );

          evalBasis( l_co,
                     m_entType,
                     l_bValCo,
                     l_localQpsXi[l_qp],
                     l_localQpsEta[l_qp],
                     l_localQpsZeta[l_qp],
                     -1, m_order );

          l_val += ( l_bValRo * l_bValCo ) * l_weights[l_qp];
        }
        if( std::abs( l_val ) > TOL.BASIS ) {
          m_flux[l_fa1].nz.push_back(l_val);
          m_flux[l_fa1].co.push_back(l_co);
          m_flux[l_fa1].ro.push_back(l_ro);
        }
      }
    }

    // derive neighboring face of regular meshes
    unsigned int l_fa2;
    if(      l_fa1 == 0 ) l_fa2 = 5;
    else if( l_fa1 == 1 ) l_fa2 = 3;
    else if( l_fa1 == 2 ) l_fa2 = 4;
    else if( l_fa1 == 3 ) l_fa2 = 1;
    else if( l_fa1 == 4 ) l_fa2 = 2;
    else if( l_fa1 == 5 ) l_fa2 = 0;

    std::vector< real_mesh > l_neighQpsXi, l_neighQpsEta, l_neighQpsZeta;

    l_neighQpsXi.resize(   m_nBaseFuncs );
    l_neighQpsEta.resize(  m_nBaseFuncs );
    l_neighQpsZeta.resize( m_nBaseFuncs );

    // map quad points to ref. element's face of adjacent element
    for( unsigned int l_qp = 0; l_qp < l_localQpsXi.size(); l_qp++ ) {
      linalg::Mappings::mapFaceCoords( m_entType,
                                       l_fa1,
                                       l_fa2,
                                       0,
                                       l_localQpsXi[l_qp],
                                       l_localQpsEta[l_qp],
                                       l_localQpsZeta[l_qp],
                                       l_neighQpsXi[l_qp],
                                       l_neighQpsEta[l_qp],
                                       l_neighQpsZeta[l_qp] );
    }

    for( unsigned int l_ro = 0; l_ro < m_nBaseFuncs; l_ro++ ) {
      for( unsigned int l_co = 0; l_co < m_nBaseFuncs; l_co++ ) {
        // value in the respective flux matrix
        real_base l_val = 0;

        for( unsigned int l_qp = 0; l_qp < l_localQpsXi.size(); l_qp++ ) {
          // evaluate basis at quad point
          real_base l_bValRo = 0;
          real_base l_bValCo = 0;

          evalBasis( l_ro,
                     m_entType,
                     l_bValRo,
                     l_neighQpsXi[l_qp],
                     l_neighQpsEta[l_qp],
                     l_neighQpsZeta[l_qp],
                     -1, m_order );

          evalBasis( l_co,
                     m_entType,
                     l_bValCo,
                     l_localQpsXi[l_qp],
                     l_localQpsEta[l_qp],
                     l_localQpsZeta[l_qp],
                     -1, m_order );

          l_val += ( l_bValRo * l_bValCo ) * l_weights[l_qp];
        }
        if( std::abs( l_val ) > TOL.BASIS ) {
          m_flux[C_ENT[m_entType].N_FACES + l_fa1].nz.push_back(l_val);
          m_flux[C_ENT[m_entType].N_FACES + l_fa1].co.push_back(l_co);
          m_flux[C_ENT[m_entType].N_FACES + l_fa1].ro.push_back(l_ro);
        }
      }
    }

  }
}

void edge::dg::Basis::computeFluxMatricesTet() {
  assert( m_entType == TET4 );
  assert( C_ENT[m_entType].N_FACES == 4 );
  assert( C_ENT[m_entType].N_FACE_VERTICES == 3 );

  // fixed size of 52 flux matrices
  m_flux.resize( 52 );

  // get triangular quadrature points
  std::vector< real_mesh > l_faQpsChi1, l_faQpsChi2, l_dummy, l_faWeights;

  dg::QuadraturePoints::getQpts( TRIA3,
                                 m_order,
                                 C_REF_ELEMENT.VE.ENT[TRIA3],
                                 l_faQpsChi1, l_faQpsChi2, l_dummy, l_faWeights );

  // iterate over local faces
  for( unsigned int l_fa1 = 0; l_fa1 < 4; l_fa1++ ) {

    // compute local flux matrices
    for( unsigned int l_ro = 0; l_ro < m_nBaseFuncs; l_ro++ ) {
      for( unsigned int l_co = 0; l_co < m_nBaseFuncs; l_co++ ) {
        // value of the flux matrix respectively
        real_base l_val = 0;

        // iterate over quadrature points
        for( unsigned int l_qp = 0; l_qp < l_faWeights.size(); l_qp++  ) {
          // volume coords
          real_mesh l_volPt[3];

          real_mesh l_faPt[2];
          l_faPt[0] = l_faQpsChi1[l_qp];
          l_faPt[1] = l_faQpsChi2[l_qp];

          // get volume coords
          linalg::Mappings::faToVolRefTet4( l_fa1,
                                            l_faPt,
                                            l_volPt );

          // evaluate basis at quad point
          real_base l_bValRo = 0;
          real_base l_bValCo = 0;

          evalBasis( l_ro,
                     TET4,
                     l_bValRo,
                     l_volPt[0],
                     l_volPt[1],
                     l_volPt[2] );

          evalBasis( l_co,
                     TET4,
                     l_bValCo,
                     l_volPt[0],
                     l_volPt[1],
                     l_volPt[2] );

          l_val += ( l_bValRo * l_bValCo ) * l_faWeights[l_qp];
        }
        if( std::abs( l_val ) > TOL.BASIS ) {
          m_flux[l_fa1].nz.push_back(l_val);
          m_flux[l_fa1].co.push_back(l_co);
          m_flux[l_fa1].ro.push_back(l_ro);
        }
      }
    }

    // iterate over neighboring faces
    for( unsigned int l_fa2 = 0; l_fa2 < 4; l_fa2++ ) {
      // iterate over the neighboring face's vertices
      for( unsigned int l_ve = 0; l_ve < 3; l_ve++ ) {
        // compute neighboring flux matrices
        for( unsigned int l_ro = 0; l_ro < m_nBaseFuncs; l_ro++ ) {
          for( unsigned int l_co = 0; l_co < m_nBaseFuncs; l_co++ ) {
            // value of the flux matrix respectively
            real_base l_val = 0;

            // iterate over quadrature points
            for( unsigned int l_qp = 0; l_qp < l_faWeights.size(); l_qp++  ) {
              // local and neighboring volume coords
              real_mesh l_loVol[3], l_neVol[3];

              // local and neighboring face coordinates; ne coords are expressed in terms of the local coords
              real_mesh l_loFa[2], l_neFa[2];

              // get local face coords
              l_loFa[0] = l_faQpsChi1[l_qp];
              l_loFa[1] = l_faQpsChi2[l_qp];

              // get 'neighboring' coords
              linalg::Mappings::faLocToFaNeiTet4( l_ve, l_loFa, l_neFa );

              // get volume coords
              linalg::Mappings::faToVolRefTet4( l_fa1, l_loFa, l_loVol );
              linalg::Mappings::faToVolRefTet4( l_fa2, l_neFa, l_neVol );


              // evaluate basis at quad point
              real_base l_bValRo = 0;
              real_base l_bValCo = 0;

              evalBasis( l_ro,
                         TET4,
                         l_bValRo,
                         l_neVol[0],
                         l_neVol[1],
                         l_neVol[2] );

              evalBasis( l_co,
                         TET4,
                         l_bValCo,
                         l_loVol[0],
                         l_loVol[1],
                         l_loVol[2] );

              l_val += ( l_bValRo * l_bValCo ) * l_faWeights[l_qp];
            }

            // derive index of flux matrix
            unsigned short l_fId = 4 + l_fa1 * 4 * 3 + l_fa2 * 3 + l_ve;

            if( std::abs( l_val ) > TOL.BASIS ) {
              m_flux[l_fId].nz.push_back(l_val);
              m_flux[l_fId].co.push_back(l_co);
              m_flux[l_fId].ro.push_back(l_ro);
            }
          }
        }
      }
    }
  }
}

void edge::dg::Basis::computeFluxMatrices() {
  if(      m_entType == LINE   ) computeFluxMatricesLine();
  else if( m_entType == QUAD4R ) computeFluxMatrices2d();
  else if( m_entType == TRIA3  ) computeFluxMatrices2d();
  else if( m_entType == HEX8R  ) computeFluxMatricesHex();
  else if( m_entType == TET4   ) computeFluxMatricesTet();
  else assert( false );
}

edge::dg::Basis::Basis( t_entityType   i_entityType,
                        unsigned short i_order ):
 m_entType(i_entityType),
 m_order(i_order) {
  m_nBaseFuncs = CE_N_ELEMENT_MODES( i_entityType, i_order );

  // init the mass matrix
  initMassMatrix();

  // compute the flux matrices
  computeFluxMatrices();

  // pre-compute basis at quad points
  initEvalQpRefEl( i_order+1 );
}

void edge::dg::Basis::print() const {
  EDGE_VLOG(2) << "shared matrices? alright:";
  EDGE_VLOG(2) << "  mass matrix:";
  if( EDGE_VLOG_IS_ON(2) ) {
    linalg::Matrix::printMatrixCrd( m_nBaseFuncs,
                                    m_nBaseFuncs,
                                    m_mass );
  }
  for( unsigned int l_flux = 0; l_flux < m_flux.size(); l_flux++ ) {
    EDGE_VLOG(2) << "  flux matrix #" << l_flux << ": ";
    if( EDGE_VLOG_IS_ON(2) ) {
      linalg::Matrix::printMatrixCrd( m_nBaseFuncs,
                                      m_nBaseFuncs,
                                      m_flux[l_flux]);
    }
  }
}

void edge::dg::Basis::initEvalQpRefEl( unsigned short i_order ) {
  m_qpEval.resize( i_order );

  // iterate over degree of quadrature
  for( unsigned short l_pq = 0; l_pq < i_order; l_pq++ ) {
    // get quad points
    dg::QuadraturePoints::getQpts( m_entType,
                                   l_pq+1,
                                   C_REF_ELEMENT.VE.ENT[m_entType],
                                   m_qpEval[l_pq].xi1,
                                   m_qpEval[l_pq].xi2,
                                   m_qpEval[l_pq].xi3,
                                   m_qpEval[l_pq].weights );

    // get #quad points
    unsigned int l_nQps = m_qpEval[l_pq].weights.size();

    m_qpEval[l_pq].basis.resize( l_nQps );

    // iterate over quad points
    for( unsigned int l_lo = 0; l_lo < l_nQps; l_lo++ ) {
      m_qpEval[l_pq].basis[l_lo].resize( i_order );

      // iterate over degree of basis
      for( unsigned short l_pb = 0; l_pb < i_order; l_pb++ ) {
        // get #modes
        unsigned int l_nModes = CE_N_ELEMENT_MODES( m_entType, l_pb+1 );

        m_qpEval[l_pq].basis[l_lo][l_pb].val.resize( l_nModes );
        for( unsigned short l_di = 0; l_di < C_ENT[m_entType].N_DIM; l_di++ ) {
          m_qpEval[l_pq].basis[l_lo][l_pb].valD[l_di].resize( l_nModes );
        }

        for( int_md l_md = 0; l_md < l_nModes; l_md++ ) {
          evalBasis( l_md, m_entType,
                     m_qpEval[l_pq].basis[l_lo][l_pb].val[l_md],
                     m_qpEval[l_pq].xi1[l_lo], m_qpEval[l_pq].xi2[l_lo], m_qpEval[l_pq].xi3[l_lo],
                     -1, l_pb+1 );

          for( unsigned short l_di = 0; l_di < C_ENT[m_entType].N_DIM; l_di++ ) {
            evalBasis( l_md, m_entType,
                       m_qpEval[l_pq].basis[l_lo][l_pb].valD[l_di][l_md],
                       m_qpEval[l_pq].xi1[l_lo], m_qpEval[l_pq].xi2[l_lo], m_qpEval[l_pq].xi3[l_lo],
                       l_di, l_pb+1 );
          }
        }
      }
    }
  }
}

void edge::dg::Basis::qpts2modal( const real_base    *i_evalF,
                                        unsigned int  i_orderQr,
                                        real_base    *o_modes ) const {
  PP_INSTR_FUN("qp2modal")

  EDGE_CHECK( i_orderQr > 0 );
  EDGE_CHECK( i_orderQr <= m_qpEval.size() );

  // iterate over modes
  for( int_md l_md = 0; l_md < m_nBaseFuncs; l_md++ ) {
    // reset mode
    o_modes[l_md] = 0;

    // iterate over quad points
    for( unsigned int l_qp = 0; l_qp < m_qpEval[i_orderQr-1].weights.size(); l_qp++ ) {
      // add contribution of the quad point
      o_modes[l_md] +=   i_evalF[l_qp]
                       * m_qpEval[i_orderQr-1].basis[l_qp][m_order-1].val[l_md]
                       * m_qpEval[i_orderQr-1].weights[l_qp];
    }

    // divide by matching mass entry
    assert( m_mass.nz.size() > l_md && (m_mass.ro[l_md] == m_mass.co[l_md]) );
    o_modes[l_md] /= m_mass.nz[l_md];
  }
}

real_base edge::dg::Basis::modal2ptval( const real_mesh  i_pt[3],
                                        const real_base *i_modes ) const {

  real_base l_val = 0;

  for( int_md l_md = 0; l_md < m_nBaseFuncs; l_md++ ) {
    real_base l_bVal;

    // evaluate basis function
    evalBasis( l_md,
               m_entType,
               l_bVal,
               i_pt[0],
               i_pt[1],
               i_pt[2] );

    // add contribution
    l_val += l_bVal * i_modes[l_md];
  }

  return l_val;
}

real_base edge::dg::Basis::modal2refPtVal(       unsigned short i_orderQr,
                                                 unsigned int   i_quadPoi,
                                           const real_base *i_modes ) const {
  real_base l_val = 0;
  for( int_md l_md = 0; l_md < m_nBaseFuncs; l_md++ ) {
    // read evaluated basis
    real_base l_bVal = m_qpEval[i_orderQr-1].basis[i_quadPoi][m_order-1].val[l_md];

    // add contribution
    l_val += l_bVal * i_modes[l_md];
  }

  return l_val;
}

void edge::dg::Basis::getMassInvDense( int_md     i_nModes,
                                       real_base *o_matrix,
                                       bool       i_rowMajor ) const {
  // check that we have enough basis functions
  CHECK( i_nModes <= m_nBaseFuncs );

  // non-zero entry of stiffness matrix
  unsigned int l_nzId = 0;

  // iterate mass matrix entries
  for( unsigned int l_ro = 0; l_ro < i_nModes; l_ro++ ) {
    for( unsigned int l_co = 0; l_co < i_nModes; l_co++ ) {
      // check for a match
      if( m_mass.nz.size() > l_nzId &&
          m_mass.ro[l_nzId] == l_ro &&
          m_mass.co[l_nzId] == l_co ) {
        EDGE_CHECK( l_ro == l_co );

        o_matrix[l_ro*i_nModes + l_co] = 1 / m_mass.nz[l_nzId];
        l_nzId++;
      }
      else
        o_matrix[l_ro*i_nModes + l_co] = 0;
    }
  }

  // change from row major to column major if requested
  if( i_rowMajor == false ) {
    linalg::Matrix::transposeDense( i_nModes, o_matrix );
  }
}

void edge::dg::Basis::getStiffMm1Dense( int_md      i_nModes,
                                        real_base  *o_matrices,
                                        bool        i_tStiff ) const {
  if( i_tStiff == false ) {
    // check that the size matches
    EDGE_CHECK_EQ( pre::dg::g_stiffVSize, i_nModes*i_nModes*C_ENT[m_entType].N_DIM );

    // set values
    for( std::size_t l_va = 0; l_va < pre::dg::g_stiffVSize; l_va++ )
      o_matrices[l_va] = pre::dg::g_stiffVRaw[l_va];
  }
  else {
    // check that size matches
    EDGE_CHECK_EQ( pre::dg::g_stiffTSize, i_nModes*i_nModes*C_ENT[m_entType].N_DIM );

    // set values
    for( std::size_t l_va = 0; l_va < pre::dg::g_stiffTSize; l_va++ )
      o_matrices[l_va] = pre::dg::g_stiffTRaw[l_va];
  }
}

void edge::dg::Basis::getFluxMm1Dense( int_md     i_nModes,
                                       real_base *o_matrices,
                                       bool       i_rowMajor ) const {
  // check that we have enough basis functions
  CHECK( i_nModes <= m_nBaseFuncs );

  // iterate over flux matrices
  for( unsigned int l_fm = 0; l_fm < m_flux.size(); l_fm++ ) {
    // non-zero entry of the flux matrix
    unsigned int l_nzId = 0;

    // iterate dense entries
    for( unsigned int l_ro = 0; l_ro < i_nModes; l_ro++ ) {
      for( unsigned int l_co = 0; l_co < i_nModes; l_co++ ) {
        unsigned int l_mId = l_fm * (i_nModes*i_nModes) + l_ro * i_nModes + l_co;

        // check for a match
        if( m_flux[l_fm].nz.size() > l_nzId &&
            m_flux[l_fm].ro[l_nzId] == l_ro &&
            m_flux[l_fm].co[l_nzId] == l_co ) {
          o_matrices[l_mId] = m_flux[l_fm].nz[l_nzId];
          l_nzId++;
        }
        else
          o_matrices[l_mId] = 0;
      }
    }

    // multiply with inverse mass matrix from the right.
    for( unsigned int l_ro = 0; l_ro < i_nModes; l_ro++ ) {
      // check for a diagonal mass matrix
      CHECK( m_mass.nz.size() >= i_nModes && m_mass.ro[l_ro] == m_mass.co[l_ro] && m_mass.ro[l_ro] == l_ro );

      for( unsigned int l_co = 0; l_co < i_nModes; l_co++ ) {
        // divide by entry of diagonal mass matrix
        assert( m_mass.nz[l_ro] > TOL.BASIS );

        unsigned int l_mId = l_fm * (i_nModes*i_nModes) + l_ro * i_nModes + l_co;
        o_matrices[l_mId] /= m_mass.nz[l_co];
      }
    }

    // change from row major to column major if requested
    if( i_rowMajor == false ) {
      linalg::Matrix::transposeDense( i_nModes, o_matrices+(l_fm * (i_nModes*i_nModes)) );
    }
  }
}

void edge::dg::Basis::evalBasisLine( unsigned int  b,
                                     real_mesh     xi,
                                     real_base    &val ) {
#include "generated/basis_line.inc"
}

void edge::dg::Basis::evalBasisDerLine( unsigned int  b,
                                        real_mesh     xi,
                                        real_base    &valDxi ) {
#include "generated/basis_line_der.inc"
}

void edge::dg::Basis::evalBasisQuad( unsigned int  i_b,
                                     real_mesh     i_xi,
                                     real_mesh     i_eta,
                                     real_base    &o_val,
                                     int           i_der,
                                     unsigned int  i_order ) {
  // TODO: move to a hierachical formulation

  // derive tensor indices of rows
  unsigned int l_tIdXi  = i_b % i_order;
  unsigned int l_tIdEta = i_b / i_order;

  // evaluate the 1D basis accordingly
  real_base l_baseEval1D[2];
  if( i_der == 0 ) {
    evalBasisDerLine( l_tIdXi,  i_xi,  l_baseEval1D[0] );
  }
  else {
    evalBasisLine( l_tIdXi,  i_xi,  l_baseEval1D[0] );
  }

  if( i_der == 1 ) {
    evalBasisDerLine( l_tIdEta, i_eta, l_baseEval1D[1] );
  }
  else {
    evalBasisLine( l_tIdEta, i_eta, l_baseEval1D[1] );
  }

  o_val = l_baseEval1D[0] * l_baseEval1D[1];
}

void edge::dg::Basis::evalBasisTria( unsigned int  b,
                                     real_mesh     xi_1,
                                     real_mesh     xi_2,
                                     real_base    &o_val,
                                     int           i_der ) {
  // temporary vars
  double val, valDxi1, valDxi2;
  val = valDxi1 = valDxi2 = std::numeric_limits< double >::max();

  if( i_der == 0 ) {
#include "generated/basis_tria_der_xi1.inc"
    o_val = valDxi1;
  }
  else if( i_der == 1 ) {
#include "generated/basis_tria_der_xi2.inc"
    o_val = valDxi2;
  }
  else {
#include "generated/basis_tria.inc"
    o_val = val;
  }
}

void edge::dg::Basis::evalBasisTet( unsigned int  b,
                                    real_mesh     xi_1,
                                    real_mesh     xi_2,
                                    real_mesh     xi_3,
                                    real_base    &o_val,
                                    int           i_der ) {
  // temporary vars
  double val, valDxi1, valDxi2, valDxi3;
  val = valDxi1 = valDxi2 = std::numeric_limits< double >::max();

  if( i_der == 0 ) {
#include "generated/basis_tet_der_xi1.inc"
    o_val = valDxi1;
  }
  else if( i_der == 1 ) {
#include "generated/basis_tet_der_xi2.inc"
    o_val = valDxi2;
  }
  else if( i_der == 2 ) {
#include "generated/basis_tet_der_xi3.inc"
    o_val = valDxi3;
  }
  else {
#include "generated/basis_tet.inc"
    o_val = val;
  }
}

void edge::dg::Basis::evalBasisHex( unsigned int  i_b,
                                    real_mesh     i_xi,
                                    real_mesh     i_eta,
                                    real_mesh     i_zeta,
                                    real_base    &o_val,
                                    int           i_der,
                                    unsigned int  i_order ) {
  // TODO: move to a hierachical formulation

  // derive tensor indices of rows
  unsigned int l_tIdXi   =   i_b %  i_order;
  unsigned int l_tIdEta  = ( i_b % (i_order*i_order) ) / i_order;
  unsigned int l_tIdZeta =   i_b / (i_order*i_order);

  // evaluate the 1D basis accordingly
  real_base l_baseEval1D[3];
  if( i_der == 0 ) {
    evalBasisDerLine( l_tIdXi,  i_xi,  l_baseEval1D[0] );
  }
  else {
    evalBasisLine( l_tIdXi,  i_xi,  l_baseEval1D[0] );
  }

  if( i_der == 1 ) {
    evalBasisDerLine( l_tIdEta, i_eta, l_baseEval1D[1] );
  }
  else {
    evalBasisLine( l_tIdEta, i_eta, l_baseEval1D[1] );
  }

  if( i_der == 2 ) {
    evalBasisDerLine( l_tIdZeta, i_zeta, l_baseEval1D[2] );
  }
  else {
    evalBasisLine( l_tIdZeta, i_zeta, l_baseEval1D[2] );
  }

  o_val = l_baseEval1D[0] * l_baseEval1D[1] *  l_baseEval1D[2];
}

void edge::dg::Basis::evalBasis( unsigned int  i_b,
                                 t_entityType  i_entType,
                                 real_base    &o_val,
                                 real_mesh     i_xi,
                                 real_mesh     i_eta,
                                 real_mesh     i_zeta,
                                 int           i_der,
                                 unsigned int  i_order ) {
  PP_INSTR_FUN("eval_basis")

  if( i_entType == LINE )
    if( i_der != 1 ) evalBasisLine( i_b, i_xi, o_val );
    else             evalBasisDerLine( i_b, i_xi, o_val );
  else if( i_entType == QUAD4R )
    evalBasisQuad( i_b,
                   i_xi,
                   i_eta,
                   o_val,
                   i_der,
                   i_order );
  else if( i_entType == TRIA3 )
    evalBasisTria( i_b,
                   i_xi,
                   i_eta,
                   o_val,
                   i_der );
  else if( i_entType == HEX8R ) {
    evalBasisHex( i_b,
                  i_xi,
                  i_eta,
                  i_zeta,
                  o_val,
                  i_der,
                  i_order );
  }
  else if( i_entType == TET4 ) {
    evalBasisTet( i_b,
                  i_xi,
                  i_eta,
                  i_zeta,
                  o_val,
                  i_der );
  }
  else
    EDGE_LOG_FATAL << "element type not supported";
}

void edge::dg::Basis::evalBasis( unsigned int        i_b,
                                 t_entityType        i_enType,
                                 real_base          &o_val,
                                 real_mesh*   const  i_pt,
                                 int                 i_der,
                                 unsigned int        i_order ) {
  PP_INSTR_FUN("eval_basis")

  if( i_enType == LINE )
    if( i_der != 1 ) evalBasisLine( i_b, i_pt[0], o_val );
    else             evalBasisDerLine( i_b, i_pt[0], o_val );
  else if( i_enType == QUAD4R )
    evalBasisQuad( i_b,
                   i_pt[0],
                   i_pt[1],
                   o_val,
                   i_der,
                   i_order );
  else if( i_enType == TRIA3 )
    evalBasisTria( i_b,
                   i_pt[0],
                   i_pt[1],
                   o_val,
                   i_der );
  else if( i_enType == HEX8R ) {
    evalBasisHex( i_b,
                  i_pt[0],
                  i_pt[1],
                  i_pt[2],
                  o_val,
                  i_der,
                  i_order );
  }
  else if( i_enType == TET4 ) {
    evalBasisTet( i_b,
                  i_pt[0],
                  i_pt[1],
                  i_pt[2],
                  o_val,
                  i_der );
  }
  else
    EDGE_LOG_FATAL << "element type not supported";
}
