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
 * This file contains headers for the FaultRegion class, and derived classes
 * that model specific dynamic rupture benchmarks. It is intended that more
 * derived classes will be added to this file in order to model a variety of
 * benchmarks.
 *
 * The base class "FaultRegion" represents a rectangular region of a planar
 * fault. To represent a region for a specific benchmark, create a derived class
 * like those below.  YOU MUST define getQuantity() to whatever makes sense
 * for your example, since it is a pure virtual method in FaultRegion.
 **/

#ifndef EDGE_V_FAULT_REGION_H
#define EDGE_V_FAULT_REGION_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <map>
#include <cmath>
#include <vector>


// x direction is along strike
// y direction is along dip, into Earth
// z is normal to plane
typedef struct xyz_point_t{
  double x, y, z;
} xyz_point_t;


// Abstract base class
class FaultRegion{
public:
  double m_xMin, m_xMax, m_yMin, m_yMax, m_faultAngle;
  std::vector< std::string > m_qtyNames;

  FaultRegion( double, double, double, double );
  bool contains( xyz_point_t );
  std::string getAllQuantities( xyz_point_t );

private:
  virtual double getQuantity( std::string, xyz_point_t ) = 0;
};


// *************** Derived Classes Below *******************
// *********************************************************
// (these represent fault regions for particular benchmarks)

// Basic region where all quantities are constant, and are set at runtime
class FR_static : public FaultRegion{
public:
  using FaultRegion::FaultRegion;
  void insert( std::string, double );

private:
  std::map< std::string, double > m_quantities;

  double getQuantity( std::string, xyz_point_t );
};

// Nucleation Patch for TPV5
// http://scecdata.usc.edu/cvws/tpv5docs.html
class FR_tpv5Nucleate : public FR_static{
public:
  FR_tpv5Nucleate( double, double, double, double );
};

// Lower Stress Patch for TPV5
// http://scecdata.usc.edu/cvws/tpv5docs.html
class FR_tpv5LowStress : public FR_static{
public:
  FR_tpv5LowStress( double, double, double, double );
};

// Higher Stress Patch for TPV5
// http://scecdata.usc.edu/cvws/tpv5docs.html
class FR_tpv5HighStress : public FR_static{
public:
  FR_tpv5HighStress( double, double, double, double );
};

#endif // FAULT_REGION_H
