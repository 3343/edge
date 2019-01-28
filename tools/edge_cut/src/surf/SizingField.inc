/**
 * @file This file is part of EDGE.
 *
 * @author David Lenz (dlenz AT ucsd.edu)
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
 * Class describing how mesh refinement varies throughout the domain
 **/
 template< class Mesh_domain >
 edge_cut::surf::SizingField<Mesh_domain>::SizingField( FT i_innerVal,
                                                        FT i_scale,
                                                        Point_3 i_center,
                                                        FT i_innerRad,
                                                        FT i_outerRad,
                                                        std::map< FT, FT, std::greater<FT> > i_layers ) :
  m_innerVal( i_innerVal ),
  m_scale( i_scale ),
  m_center( i_center ),
  m_innerRad( i_innerRad ),
  m_outerRad( i_outerRad ),
  m_layers( i_layers )
{
  EDGE_CHECK_NE( m_scale, 0 );
  
  if ( m_innerRad != m_outerRad ) {
    m_widthInv = 1 / ( m_outerRad - m_innerRad );
  }
  else {
    m_widthInv = 0;
  }
}


template< class Mesh_domain >
typename edge_cut::surf::SizingField<Mesh_domain>::FT
edge_cut::surf::SizingField<Mesh_domain>::operator()( Point_3 const & i_p,
                                                      int     const,
                                                      Index   const &     ) const
{
  FT l_distance = std::sqrt( std::pow( i_p.x()-m_center.x(), 2 ) + std::pow( i_p.y()-m_center.y(), 2 ) );
  FT l_depthScale = getDepthScale( i_p.z() );

  if ( l_distance <= m_innerRad ) {
    return l_depthScale * m_innerVal;
  }
  else if ( l_distance >= m_outerRad ) {
    return m_scale * l_depthScale * m_innerVal;
  }
  else {
    // Control should never reach here if inner and outer radii are equal
    EDGE_CHECK_NE( m_outerRad, m_innerRad );
    return ( 1 + ( l_distance - m_innerRad ) * ( m_scale - 1) * m_widthInv ) * l_depthScale * m_innerVal;
  }
}


template< class Mesh_domain >
typename edge_cut::surf::SizingField<Mesh_domain>::FT
edge_cut::surf::SizingField<Mesh_domain>::getDepthScale( FT i_depth ) const {
  FT l_depthScale;
  if ( m_layers.size() == 0 ) {
    l_depthScale = 1;
  }
  else {
    typename std::map< FT, FT, std::greater<FT> >::const_iterator l_layerIt = m_layers.begin();
    while ( i_depth <= l_layerIt->first && l_layerIt != m_layers.end() ) {
      l_layerIt++;
    }

    if ( l_layerIt == m_layers.end() ) {
      l_depthScale = (--m_layers.end())->second;
    }
    else {
      l_depthScale = l_layerIt->second;
    }
  }

  return l_depthScale;
}
