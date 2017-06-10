##
# @file This file is part of EDGE.
#
# @author Alexander Breuer (anbreuer AT ucsd.edu)
#
# @section LICENSE
# Copyright (c) 2016, Regents of the University of California
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# @section DESCRIPTION
# Projects from WGS84 to geocentric coordinates.
##

import logging
import pyproj
import argparse
import netCDF4
import math
import numpy

# set up logger
logging.basicConfig( level=logging.DEBUG,
                     format='%(asctime)s - %(name)s - %(levelname)s - %(message)s' )
# command line arugments
l_parser    = argparse.ArgumentParser( description='Projects the given file to geocentric coordinates.' )
l_parser.add_argument( '--origin',
                       dest     = "origin",
                       nargs    = 3,
                       type     = float,
                       required = True,
                       help     = "Origin in the target coordinate sytem.",
                       metavar  = "LON LAT HEIGHT" )
l_parser.add_argument( '--in',
                       dest     = "in",
                       required = True,
                       help     = "Input file in netCDF format and lon/lat coords (WGS84).",
                       metavar  = "IN_FILE" )
l_parser.add_argument( '--out',
                       dest     = "out",
                       required = True,
                       help     = "Output file (xyz). Result will be written in geocentric coords (meters).",
                       metavar  = "OUT_FILE" )

l_args = vars(l_parser.parse_args())

# open file
logging.info( 'reading ' + l_args['in'] )
l_ncIn = netCDF4.Dataset( l_args['in'], 'r' )
for l_name in l_ncIn.ncattrs():
  l_attr = str( getattr( l_ncIn, l_name ) )
  logging.info( '  Global attribute ' + l_name + ' = ' +  l_attr )

# read data
l_dataIn = {}
l_dataIn['lat'] = l_ncIn.variables['lat'][:]
l_dataIn['lon'] = l_ncIn.variables['lon'][:]
l_dataIn['z']   = l_ncIn.variables['z'][:][:]

# setup the projection
l_projWgs84      = pyproj.Proj('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs') ##4326
l_projGeocentric = pyproj.Proj('+proj=geocent +datum=WGS84 +units=m     +no_defs') ##4978

# derive coords in meters for the origin
l_orig = numpy.zeros(3)
l_orig[0], l_orig[1], l_orig[2] = pyproj.transform( l_projWgs84,
                                                    l_projGeocentric,
                                                    l_args['origin'][0],
                                                    l_args['origin'][1],
                                                    l_args['origin'][2] )

# define basis for a local east-north-up coordinate system
l_locBas = numpy.zeros((3,3))
l_locBas[0][2], l_locBas[1][2], l_locBas[2][2] = pyproj.transform( l_projWgs84,
                                                                   l_projGeocentric,
                                                                   l_args['origin'][0],
                                                                   l_args['origin'][1],
                                                                   l_args['origin'][2] + 6378137 )
l_locBas[0][0], l_locBas[1][0], l_locBas[2][0] = pyproj.transform( l_projWgs84,
                                                                   l_projGeocentric,
                                                                   l_args['origin'][0]+0.001,
                                                                   l_args['origin'][1],
                                                                   l_args['origin'][2] )
l_locBas[0][1], l_locBas[1][1], l_locBas[2][1] = pyproj.transform( l_projWgs84,
                                                                   l_projGeocentric,
                                                                   l_args['origin'][0],
                                                                   l_args['origin'][1]+0.001,
                                                                   l_args['origin'][2] )

# normalize basis to length 1
for l_b in xrange(3):
  l_locBas[0][l_b] = l_locBas[0][l_b] - l_orig[0]
  l_locBas[1][l_b] = l_locBas[1][l_b] - l_orig[1]
  l_locBas[2][l_b] = l_locBas[2][l_b] - l_orig[2]
  l_norm = math.sqrt(l_locBas[0][l_b]**2 + l_locBas[1][l_b]**2 + l_locBas[2][l_b]**2 )

  for l_dim in xrange(3):
    l_locBas[l_dim][l_b] = l_locBas[l_dim][l_b] / l_norm

# transformation from geocentric basis in local coordinate system
l_trafo = numpy.linalg.inv( l_locBas )

# create output file
l_xyzOut = open( l_args['out'], 'w' )

logging.info( 'doing the ' +
              str(len(l_dataIn['lat'])*len(l_dataIn['lon'])) +
              ' projections..' )

# iterate over data
for l_lat in xrange( len(l_dataIn['lat']) ):
  for l_lon in xrange( len(l_dataIn['lon']) ):
    # do the projection
    l_coords = numpy.zeros(3)
    l_coords[0], l_coords[1], l_coords[2] = pyproj.transform( l_projWgs84,
                                                              l_projGeocentric,
                                                              l_dataIn['lon'][l_lon],
                                                              l_dataIn['lat'][l_lat],
                                                              l_dataIn['z'][l_lat][l_lon] )

    # express relative to origin
    l_coords = l_coords - l_orig

    # change to local coordinate system
    l_coords = numpy.dot( l_trafo, l_coords )

    # write the data
    l_xyzOut.write( str(l_coords[0]) + ' ' + str(l_coords[1]) + ' ' + str(l_coords[2]) + '\n' )


l_ncIn.close()
l_xyzOut.close()

logging.info( 'done' )
