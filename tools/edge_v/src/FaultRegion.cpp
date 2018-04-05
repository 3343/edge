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
 * This file contains implementations for the FaultRegion class, and for derived
 * classes that model specific dynamic rupture benchmarks. It is intended that
 * more derived classes will be added to this file in order to model a variety
 * of benchmarks.
 *
 * The base class "FaultRegion" represents a rectangular region of a planar
 * fault. To represent a region for a specific benchmark, create a derived class
 * like those below.  YOU MUST define getQuantity() to whatever makes sense
 * for your example, since it is a pure virtual method in FaultRegion.
 **/

#include "FaultRegion.h"

FaultRegion::FaultRegion( double i_xMin,
                          double i_xMax,
                          double i_yMin,
                          double i_yMax ) :
  m_xMin( i_xMin ),
  m_xMax( i_xMax ),
  m_yMin( i_yMin ),
  m_yMax( i_yMax )
{ }

// Assumes that i_pt lies within the fault plane
bool FaultRegion::contains( xyz_point_t i_pt ){
  bool l_contained = false;

  if( i_pt.x >= m_xMin && i_pt.x <= m_xMax &&
      i_pt.y >= m_yMin && i_pt.y <= m_yMax    )
    l_contained = true;

  return l_contained;
}

std::string FaultRegion::getAllQuantities( xyz_point_t i_pt ){
  std::stringstream l_qtyString;
  for( const auto& l_qName : m_qtyNames ){
    l_qtyString << getQuantity( l_qName, i_pt ) << " ";
    // l_qtyString << m_quantities.find( l_qName )->second << " ";
  }

  return l_qtyString.str();
}


// *~~~~~~~~~~~ FR_static ~~~~~~~~~~~*
void FR_static::insert( std::string i_s, double i_d ){
  m_qtyNames.push_back( i_s );
  m_quantities.insert( std::pair< std::string, double >( i_s, i_d ) );
}

double FR_static::getQuantity( std::string i_qName, xyz_point_t i_pt ){
  return m_quantities.find( i_qName )->second;
}


// *~~~~~~~~~~~ FR_tpv5 ~~~~~~~~~~~*
FR_tpv5Nucleate::FR_tpv5Nucleate( double i_xMin,
                                  double i_xMax,
                                  double i_yMin,
                                  double i_yMax ) :
  FR_static( i_xMin, i_xMax, i_yMin, i_yMax )
{
  insert( "S_STRESS_STK",  81.6E6 );
  insert( "S_STRESS_DIP",   0.0   );
  insert( "N_STRESS",    -120.0E6 );
}

FR_tpv5LowStress::FR_tpv5LowStress( double i_xMin,
                                    double i_xMax,
                                    double i_yMin,
                                    double i_yMax ) :
  FR_static( i_xMin, i_xMax, i_yMin, i_yMax )
{
  insert( "S_STRESS_STK",   62.0E6 );
  insert( "S_STRESS_DIP",    0.0   );
  insert( "N_STRESS",     -120.0E6 );
}

FR_tpv5HighStress::FR_tpv5HighStress( double i_xMin,
                                      double i_xMax,
                                      double i_yMin,
                                      double i_yMax ) :
  FR_static( i_xMin, i_xMax, i_yMin, i_yMax )
{
  insert( "S_STRESS_STK",   78.0E6 );
  insert( "S_STRESS_DIP",    0.0   );
  insert( "N_STRESS",     -120.0E6 );
}
