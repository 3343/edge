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
 * Extrudes a given 2D rectangular surface mesh.
 **/
#include "Extrude.h"
#include <limits>
#include <algorithm>
#include "io/logging.hpp"

void edge_cut::mesh::Extrude::addSide( std::size_t                  i_nSteps,
                                       double                       i_zTarget,
                                       std::vector< point > const & i_veCrdsBnd,
                                       std::vector< point >       & o_veCrds ) {
  for( std::size_t l_st = 0; l_st < i_nSteps+1; l_st++ ) {
    for( std::size_t l_ve = 0; l_ve < i_veCrdsBnd.size(); l_ve++ ) {
      double l_z = std::get<2>( i_veCrdsBnd[l_ve] );
      double l_dx = (i_zTarget - l_z) / i_nSteps;
             l_z += l_st * l_dx;

      o_veCrds.push_back( std::make_tuple( std::get<0>( i_veCrdsBnd[l_ve] ),
                                           std::get<1>( i_veCrdsBnd[l_ve] ),
                                           l_z ) );
    }
  }
}

unsigned short edge_cut::mesh::Extrude::writeTriasBlTr( std::size_t    i_off,
                                                        std::size_t    i_strideX,
                                                        std::size_t    i_strideY,
                                                        std::ostream & o_tria ) {
  o_tria << "3 " << i_off
         <<  " " << i_off + i_strideX
         <<  " " << i_off + i_strideX + i_strideY
         << "\n";

  o_tria << "3 " << i_off
         <<  " " << i_off + i_strideX + i_strideY
         <<  " " << i_off +             i_strideY
         << "\n";

  return 2;
}

unsigned short edge_cut::mesh::Extrude::writeTriasBrTl( std::size_t    i_off,
                                                        std::size_t    i_strideX,
                                                        std::size_t    i_strideY,
                                                        std::ostream & o_tria ) {
  o_tria << "3 " << i_off
         <<  " " << i_off + i_strideX
         <<  " " << i_off +             i_strideY
         << "\n";
  o_tria << "3 " << i_off + i_strideX
         <<  " " << i_off + i_strideX + i_strideY
         <<  " " << i_off +             i_strideY
         << "\n";

  return 2;
}

unsigned short edge_cut::mesh::Extrude::writeTriasCr( std::size_t    i_off,
                                                      std::size_t    i_strideX,
                                                      std::size_t    i_strideY,
                                                      std::ostream & o_tria ) {
  o_tria << "3 " << i_off
         <<  " " << i_off + 2*i_strideX
         <<  " " << i_off + 2*i_strideX + 1*i_strideY
         << "\n";
  o_tria << "3 " << i_off
         <<  " " << i_off + 2*i_strideX + 1*i_strideY
         <<  " " << i_off +               2*i_strideY
         << "\n";
  o_tria << "3 " << i_off + 2*i_strideX + 1*i_strideY
         <<  " " << i_off + 2*i_strideX + 2*i_strideY
         <<  " " << i_off +               2*i_strideY
         << "\n";

  return 3;
}

unsigned short edge_cut::mesh::Extrude::writeTriasCl( std::size_t    i_off,
                                                      std::size_t    i_strideX,
                                                      std::size_t    i_strideY,
                                                      std::ostream & o_tria ) {
  o_tria << "3 " << i_off
         <<  " " << i_off + 2*i_strideX
         <<  " " << i_off +               1*i_strideY
         << "\n";
  o_tria << "3 " << i_off + 2*i_strideX
         <<  " " << i_off + 2*i_strideX + 2*i_strideY
         <<  " " << i_off +               1*i_strideY
         << "\n";
  o_tria << "3 " << i_off +               1*i_strideY
         <<  " " << i_off + 2*i_strideX + 2*i_strideY
         <<  " " << i_off +               2*i_strideY
         << "\n";

  return 3;
}

unsigned short edge_cut::mesh::Extrude::writeTriasCt( std::size_t    i_off,
                                                      std::size_t    i_strideX,
                                                      std::size_t    i_strideY,
                                                      std::ostream & o_tria ) {
  o_tria << "3 " << i_off
         <<  " " << i_off + 1*i_strideX + 2*i_strideY
         <<  " " << i_off +               2*i_strideY
         << "\n";
  o_tria << "3 " << i_off
         <<  " " << i_off + 2*i_strideX
         <<  " " << i_off + 1*i_strideX + 2*i_strideY
         << "\n";
  o_tria << "3 " << i_off + 2*i_strideX
         <<  " " << i_off + 2*i_strideX + 2*i_strideY
         <<  " " << i_off + 1*i_strideX + 2*i_strideY
         << "\n";

  return 3;
}

unsigned short edge_cut::mesh::Extrude::writeTriasCb( std::size_t    i_off,
                                                      std::size_t    i_strideX,
                                                      std::size_t    i_strideY,
                                                      std::ostream & o_tria ) {
  o_tria << "3 " << i_off
         <<  " " << i_off + 1*i_strideX
         <<  " " << i_off +               2*i_strideY
         << "\n";
  o_tria << "3 " << i_off + 1*i_strideX
         <<  " " << i_off + 2*i_strideX + 2*i_strideY
         <<  " " << i_off               + 2*i_strideY
         << "\n";
  o_tria << "3 " << i_off + 1*i_strideX
         <<  " " << i_off + 2*i_strideX
         <<  " " << i_off + 2*i_strideX + 2*i_strideY
         << "\n";

  return 3;
}

std::size_t edge_cut::mesh::Extrude::writeTriaOff( std::size_t      i_nx,
                                                   std::size_t      i_ny,
                                                   std::ostream   & o_tria,
                                                   unsigned short   i_nLevels ) {
  std::size_t l_nTria = 0;

  // derive outermost stride and total buffer
  std::size_t l_stride = 1;
  for( unsigned short l_lv = 0; l_lv < i_nLevels; l_lv++ ) {
    l_stride *= 2;
  }

  std::size_t l_buffer = 0;
  for( unsigned short l_lv = 0; l_lv < i_nLevels; l_lv++ ) {
    // bottom left corner
    l_nTria += writeTriasBlTr( l_buffer*i_ny + l_buffer,
                               l_stride*i_ny,
                               l_stride,
                               o_tria );

    // bottom right corner
    l_nTria += writeTriasBrTl( (i_nx-1-l_buffer-l_stride)*i_ny + l_buffer,
                               l_stride*i_ny,
                               l_stride,
                               o_tria );

    // top left corner
    l_nTria += writeTriasBrTl( l_buffer*i_ny + (i_ny-1-l_buffer-l_stride),
                               l_stride*i_ny,
                               l_stride,
                               o_tria );

    // top right corner
    l_nTria += writeTriasBlTr( (i_nx-1-l_buffer-l_stride)*i_ny + (i_ny-1-l_buffer-l_stride),
                               l_stride*i_ny,
                               l_stride,
                               o_tria );

    // add left edge
    for( std::size_t l_vy = l_buffer+l_stride; l_vy < i_ny-l_buffer-l_stride-1; l_vy += l_stride ) {
      l_nTria += writeTriasCr( l_buffer*i_ny + l_vy,
                               (l_stride/2)*i_ny,
                               (l_stride/2),
                               o_tria );
    }

    // add right edge
    for( std::size_t l_vy = l_buffer+l_stride; l_vy < i_ny-l_buffer-l_stride-1; l_vy += l_stride ) {
      l_nTria += writeTriasCl( (i_nx-1-l_stride-l_buffer)*i_ny + l_vy,
                               (l_stride/2)*i_ny,
                               (l_stride/2),
                               o_tria );
    }

    // add bottom edge
    for( std::size_t l_vx = (l_buffer+l_stride)*i_ny; l_vx < (i_nx-l_buffer-l_stride-1)*i_ny; l_vx += l_stride*i_ny ) {
      l_nTria += writeTriasCt( l_buffer + l_vx,
                               (l_stride/2)*i_ny,
                               (l_stride/2),
                               o_tria );
    }

    // add top edge
    for( std::size_t l_vx = (l_buffer+l_stride)*i_ny; l_vx < (i_nx-l_buffer-l_stride-1)*i_ny; l_vx += l_stride*i_ny ) {
      l_nTria += writeTriasCb( i_ny-1-l_buffer-l_stride + l_vx,
                               (l_stride/2)*i_ny,
                               (l_stride/2),
                               o_tria );
    }

    l_buffer += l_stride;
    l_stride /= 2;
  }

  for( std::size_t l_vx = l_buffer; l_vx < i_nx-1-l_buffer; l_vx++ ) {
    for( std::size_t l_vy = l_buffer; l_vy < i_ny-1-l_buffer; l_vy++ ) {
      l_nTria += writeTriasBlTr( l_vx*i_ny + l_vy,
                                 i_ny,
                                 1,
                                 o_tria );
    }
  }

  return l_nTria;
}

void edge_cut::mesh::Extrude::writePtsOff( std::vector< point > const & i_pts,
                                           std::ostream               & o_pts ) {
  for( std::size_t l_pt = 0; l_pt < i_pts.size(); l_pt++ ) {
    o_pts <<        std::get<0>( i_pts[l_pt] )
          << " " << std::get<1>( i_pts[l_pt] )
          << " " << std::get<2>( i_pts[l_pt] )
          << "\n";
  }
}

edge_cut::mesh::Extrude::Extrude( std::size_t          i_nVes,
                                  double      const (* i_veCrds)[3],
                                  double               i_zTarget,
                                  unsigned short       i_nLevels,
                                  double               i_epsilon ) {
  m_eps = i_epsilon;
  m_zTarget = i_zTarget;
  m_nLevels = i_nLevels;

  // init min and max
  for( unsigned short l_di = 0; l_di < 3; l_di++ ) {
    m_veMinSurf[l_di] = std::numeric_limits< double >::max();
    m_veMaxSurf[l_di] = std::numeric_limits< double >::lowest();
  }

  // derive minimum and maximum vertex coordinates, store locally
  for( std::size_t l_ve = 0; l_ve < i_nVes; l_ve++ ) {
    for( unsigned short l_di = 0; l_di < 3; l_di++ ) {
      m_veMinSurf[l_di] = std::min( m_veMinSurf[l_di], i_veCrds[l_ve][l_di] );
      m_veMaxSurf[l_di] = std::max( m_veMaxSurf[l_di], i_veCrds[l_ve][l_di] );
    }

    m_veCrdsSurf.push_back( std::make_tuple( i_veCrds[l_ve][0],
                                             i_veCrds[l_ve][1],
                                             i_veCrds[l_ve][2] ) );
  }

  // sort
  std::sort( m_veCrdsSurf.begin(), m_veCrdsSurf.end() );

  // derive number of points in x-, y-direction
  for( std::size_t l_ve = 0; l_ve < i_nVes-1; l_ve++ ) {
    double l_diff = std::get<0>( m_veCrdsSurf[l_ve] ) - std::get<0>( m_veCrdsSurf[l_ve+1] );
    l_diff = std::abs( l_diff );
    if( l_diff > m_eps) {
      m_ny[0] = l_ve+1;
      break;
    }
  }
  EDGE_CUT_CHECK_NE( m_ny[0], std::numeric_limits< std::size_t >::max() );
  EDGE_CUT_CHECK_EQ( i_nVes%m_ny[0], 0 );
  m_nx[0] = i_nVes / m_ny[0];

  // derive refinement ratio nad number of coarse cells
  for( unsigned short l_le = 0; l_le < i_nLevels; l_le++)
    m_refRatio *= 2;

  EDGE_CUT_CHECK_EQ( (m_nx[0]-1)%m_refRatio, 0 );
  EDGE_CUT_CHECK_EQ( (m_ny[0]-1)%m_refRatio, 0 );

  if( m_refRatio > 1 ) {
    m_nx[1] = m_nx[0] / m_refRatio + 1;
    m_ny[1] = m_ny[0] / m_refRatio + 1;
  }
  else {
    m_nx[1] = m_nx[0];
    m_ny[1] = m_ny[0];
  }

  // derive minimum dx and dy in the coarse mesh
  for( std::size_t l_ve = 0; l_ve < m_ny[1]-1; l_ve++ ) {
    std::size_t l_veNext = (l_ve+1) * m_refRatio;
    double l_diff = std::get<1>( m_veCrdsSurf[l_veNext] ) - std::get<1>( m_veCrdsSurf[l_ve] );
    m_dMin[1] = std::min( m_dMin[1], l_diff );
  }
  for( std::size_t l_ve = 0; l_ve < (m_nx[0]-m_refRatio)*m_ny[0]; l_ve+=m_ny[0]*m_refRatio ) {
    std::size_t l_veNext = l_ve+m_ny[0]*m_refRatio;
    double l_diff = std::get<0>( m_veCrdsSurf[l_veNext] ) - std::get<0>( m_veCrdsSurf[l_ve] );
    m_dMin[0] = std::min( m_dMin[0], l_diff );
  }

  // derive number of steps
  double l_zMaxDist = m_veMaxSurf[2] - m_zTarget;
  EDGE_CUT_CHECK_GT( l_zMaxDist, 0 );

  m_nSteps = l_zMaxDist / std::min( m_dMin[0], m_dMin[1] );
             m_nSteps += 1;

  // add front
  std::vector< point > l_tmpVes;
  for( std::size_t l_ve = 0; l_ve < m_nx[0]*m_ny[0]; l_ve += m_ny[0]*m_refRatio ) {
    l_tmpVes.push_back( m_veCrdsSurf[l_ve] );
  }
  addSide( m_nSteps,
           m_zTarget,
           l_tmpVes,
           m_veCrdsFront );

  // add back
  l_tmpVes.resize(0);
  for( std::size_t l_ve = 0; l_ve < m_nx[0]*m_ny[0]; l_ve += m_ny[0]*m_refRatio ) {
    l_tmpVes.push_back( m_veCrdsSurf[m_ny[0]-1 + l_ve] );
  }
  addSide( m_nSteps,
           m_zTarget,
           l_tmpVes,
           m_veCrdsBack );

  // add left
  l_tmpVes.resize(0);
  for( std::size_t l_ve = 0; l_ve < m_ny[0]; l_ve+=m_refRatio ) {
    l_tmpVes.push_back( m_veCrdsSurf[l_ve] );
  }
  addSide( m_nSteps,
           m_zTarget,
           l_tmpVes,
           m_veCrdsLeft );

  // add right
  l_tmpVes.resize(0);
  for( std::size_t l_ve = 0; l_ve < m_ny[0]; l_ve+=m_refRatio ) {
    l_tmpVes.push_back( m_veCrdsSurf[(m_nx[0]-1)*m_ny[0] + l_ve] );
  }
  addSide( m_nSteps,
           m_zTarget,
           l_tmpVes,
           m_veCrdsRight );
}

void edge_cut::mesh::Extrude::writeOff( std::ostream & io_left,
                                        std::ostream & io_right,
                                        std::ostream & io_front,
                                        std::ostream & io_back,
                                        std::ostream & io_bottom,
                                        std::ostream & io_top ) const {
  std::string l_header = "OFF\n";
              l_header += "#\n";
              l_header += "# EDGEcut generated surface mesh (version: ";
              l_header += PP_EDGE_VERSION;
              l_header += ")\n";
              l_header += "# EDGEcut is available from https://dial3343.org\n";
              l_header += "#\n";

  // left
  io_left << l_header;
  io_left << m_veCrdsLeft.size() << " " << (m_ny[1]-1)*m_nSteps*2 << " 0\n";
  writePtsOff( m_veCrdsLeft,
               io_left );
  writeTriaOff( m_nSteps+1,
                m_ny[1],
                io_left );

  // right
  io_right << l_header;
  io_right << m_veCrdsRight.size() << " " << (m_ny[1]-1)*m_nSteps*2 << " 0\n";
  writePtsOff( m_veCrdsRight,
               io_right );
  writeTriaOff( m_nSteps+1,
                m_ny[1],
                io_right );

  // front
  io_front << l_header;
  io_front << m_veCrdsFront.size() << " " << (m_nx[1]-1)*m_nSteps*2 << " 0\n";
  writePtsOff( m_veCrdsFront,
               io_front );
  writeTriaOff( m_nSteps+1,
                m_nx[1],
                io_front );

  // back
  io_back << l_header;
  io_back << m_veCrdsBack.size() << " " << (m_nx[1]-1)*m_nSteps*2 << " 0\n";
  writePtsOff( m_veCrdsBack,
               io_back );
  writeTriaOff( m_nSteps+1,
                m_nx[1],
                io_back );

  // bottom
  std::vector< point > l_bottom;
  for( std::size_t l_vx = 0; l_vx < m_nx[0]; l_vx+=m_refRatio ) {
    for( std::size_t l_vy = 0; l_vy < m_ny[0]; l_vy+=m_refRatio ) {
      std::size_t l_ve = l_vx * m_ny[0] + l_vy;
      l_bottom.push_back( std::make_tuple( std::get<0>( m_veCrdsSurf[l_ve] ),
                                          std::get<1>( m_veCrdsSurf[l_ve] ),
                                          m_zTarget ) );
    }
  }
  io_bottom << l_header;
  io_bottom << l_bottom.size() << " " << (m_nx[1]-1)*(m_ny[1]-1)*2 << " 0\n";
  writePtsOff( l_bottom,
               io_bottom );
  writeTriaOff( m_nx[1],
                m_ny[1],
                io_bottom );

  // top
  std::stringstream l_top;
  writePtsOff( m_veCrdsSurf,
               l_top );
  std::size_t l_nTria = writeTriaOff( m_nx[0],
                                      m_ny[0],
                                      l_top,
                                      m_nLevels );
  io_top << l_header;
  io_top << m_veCrdsSurf.size() << " " << l_nTria << " 0\n";
  io_top << l_top.str();
}