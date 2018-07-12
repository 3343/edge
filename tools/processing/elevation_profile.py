#!/usr/bin/env python
##
# @file This file is part of EDGE.
#
# @author Alexander Breuer (anbreuer AT ucsd.edu)
#
# @section LICENSE
# Copyright (c) 2018, Regents of the University of California
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
# Extracts elevation profiles from a DEM and writes them as CSV or gmsh-mesh.
##
import math
import netCDF4
import argparse
import logging

##
# Searches for the next id in the grid.
#
# @param i_x x-coordinate.
# @param i_y y-coordinate.
# @param i_gridX grid coordinates in x-direction.
# @param i_gridY grid coordinates in y-direction.
#
# @return x- and y-id.
##
def searchIds( i_x, i_y, i_gridX, i_gridY ):
  for l_x in range( len(i_gridX) ):
    if i_x < i_gridX[l_x]:
      break
  for l_y in range( len(i_gridY) ):
    if i_y < i_gridY[l_y]:
      break
  
  return l_x, l_y

# command line arguments
l_parser = argparse.ArgumentParser( description='Extracts elevation profiles from a DEM and writes them as CSV or gmsh-mesh.' )

l_parser.add_argument( '--dem',
                       dest     = 'dem',
                       required = True,
                       type     = str,
                       help     = 'Digital elevation map in netCDF format.' )

l_parser.add_argument( '--left',
                       dest     = 'left',
                       required = True,
                       nargs    = 2,
                       type     = float,
                       help     = 'Location of the left point (w.r.t. to the generated free surface) in the DEM.' )

l_parser.add_argument( '--right',
                       dest     = 'right',
                       required = True,
                       nargs    = 2,
                       type     = float,
                       help     = 'Location of the right point (w.r.t. to the generated free surface) in the DEM.' )

l_parser.add_argument( '--depth',
                       dest     = 'depth',
                       required = True,
                       type     = float,
                       help     = 'Depth of the mesh.' )

l_parser.add_argument( '--sponge',
                       dest     = 'sponge',
                       required = True,
                       type     = float,
                       help     = 'Size of the sponge zone on both sides (added to the left and right boundaries).' )

l_parser.add_argument( '--sampling',
                       dest     = 'sampling',
                       required = True,
                       nargs    = 2,
                       type     = float,
                       help     = 'Sampling of the topography at the surface in domain of interest (bounded by left and right) and sponge zone.' )

l_parser.add_argument( '--offset',
                       dest     = 'offset',
                       default  = [0, 0],
                       nargs    = 2,
                       type     = float,
                       help     = 'Offset, applied to the coordinates of the final mesh. If the offset is (0,0) (or not given), the left boundary is zero in x-direction. Zero in y-direction is equivalent to DEM (typically sea level).' )

l_parser.add_argument( '--output',
                       dest     = 'output',
                       required = True,
                       nargs    = '+',
                       type     = str,
                       help     = 'Path to output files. \'.geo\' file ending will generate a gmsh-script. \'.csv\' file ending will generate a CSV with the elevation profile' )

l_args = vars(l_parser.parse_args())

# set up logger
logging.basicConfig( level=logging.INFO,
                     format='%(asctime)s - %(name)s - %(levelname)s - %(message)s' )

# length of the given segment
l_lengthSeq = math.sqrt( (l_args['right'][0] - l_args['left'][0])**2 + (l_args['right'][1] - l_args['left'][1])**2 )

# length including the sponge zone
l_lengthSpo = l_lengthSeq + 2 * l_args['sponge']

# derive unit vector from left to right
l_unit = []
for l_di in range(2):
  l_unit = l_unit + [ l_args['right'][l_di] - l_args['left'][l_di] ]
l_unit = [ l_en / l_lengthSeq for l_en in l_unit ]

# derive the corners
l_cornersSurf = [ l_args['left'], l_args['right'], l_args['right'], l_args['left']  ]
for l_di in range(2):
  l_cornersSurf[0][l_di] = l_cornersSurf[0][l_di] - l_args['sponge']*l_unit[l_di]
  l_cornersSurf[1][l_di] = l_cornersSurf[1][l_di] + l_args['sponge']*l_unit[l_di]


# print startup info
logging.info( 'generating mesh..')
logging.info( 'length of segment: ' + str(l_lengthSeq) )
logging.info( 'length including sponge zone: ' + str(l_lengthSpo) )
logging.info( 'left corner at the surface (including sponge):  ' + str(l_cornersSurf[0]) )
logging.info( 'right corner at the surface (including sponge): ' + str(l_cornersSurf[1]) )

# open netCDF file
l_topo = netCDF4.Dataset( l_args['dem'] )

# gmsh geometry
l_gmshGeo = ""
l_csv = "x,y\n"

l_pt = 0
l_pt1d = 0

while( l_pt1d <= l_lengthSpo ):
  l_ptSrc = [0,0]
  l_ptSrc[0] = l_cornersSurf[0][0] + l_unit[0] * l_pt1d
  l_ptSrc[1] = l_cornersSurf[0][1] + l_unit[1] * l_pt1d

  l_x, l_y = searchIds( l_ptSrc[0], l_ptSrc[1], l_topo['x'], l_topo['y'] )

  # add to gmsh's geo file
  l_gmshGeo = l_gmshGeo + "Point(" + str(l_pt+1) + ") = { " + str(l_pt1d + l_args['offset'][0]) + ", " + str(l_topo['z'][l_x][l_y] + l_args['offset'][1]) + ", 0 };\n"
  l_csv     = l_csv + str(l_pt1d + l_args['offset'][0]) + ", " + str(l_topo['z'][l_x][l_y] + l_args['offset'][1]) + '\n'

  # update point coordinates
  if( l_pt1d < l_args['sponge'] or l_pt1d >= l_lengthSeq+l_args['sponge'] ):
    l_ptCheck = l_pt1d + l_args['sampling'][1]
  else:
    l_ptCheck = l_pt1d + l_args['sampling'][0]
  l_pt = l_pt + 1

  # add bottom corners in last iteration
  if(  l_ptCheck > l_lengthSpo ):
    l_gmshGeo = l_gmshGeo + "Point(" + str(l_pt+1) + ") = { " + str(l_pt1d + l_args['offset'][0]) + ", " + str(l_args['depth'] + l_args['offset'][1]) + ", 0 };\n"
    l_gmshGeo = l_gmshGeo + "Point(" + str(l_pt+2) + ") = { " + str(l_args['offset'][0])          + ", " + str(l_args['depth'] + l_args['offset'][1]) + ", 0 };\n"

  # update 1d point
  l_pt1d = l_ptCheck

# add the lines
# 1) surface
l_gmshGeo = l_gmshGeo + "BSpline(0) = { 1:" + str(l_pt) + " };\n"
# 2) remainder
l_gmshGeo = l_gmshGeo + "Line(1) = { "+str(l_pt  ) + ',' + str(l_pt+1)  + " };\n"
l_gmshGeo = l_gmshGeo + "Line(2) = { "+str(l_pt+1) + ',' + str(l_pt+2)  + " };\n"
l_gmshGeo = l_gmshGeo + "Line(3) = { "+str(l_pt+2) +                    ",1 };\n"

# write output
for l_path in l_args['output']:
  if( 'geo' in l_path ):
    with open(l_path, 'w' ) as l_file:
      l_file.write( l_gmshGeo )
  elif( 'csv' in l_path ):
    with open(l_path, 'w' ) as l_file:
      l_file.write( l_csv )
  else:
    logging.error( 'don\'t know what to do with' + l_path + ', ignoring..' )
