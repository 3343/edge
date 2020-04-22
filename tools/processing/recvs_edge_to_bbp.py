#!/usr/bin/env python3
##
# @file This file is part of EDGE.
#
# @author Alexander Breuer (breuer AT mytum.de)
#
# @section LICENSE
# Copyright (c) 2020, Alexander Breuer
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
# Converts one of EDGE's seismograms to the format of the broadband platform (*.bbp).
##
import logging
import argparse
import pandas

# setup logger
logging.basicConfig( level=logging.INFO,
                     format='%(asctime)s - %(name)s - %(levelname)s - %(message)s' )

# command line arguments
l_parser = argparse.ArgumentParser( formatter_class=argparse.RawTextHelpFormatter,
                                    description='Converts one a seismogram in EDGEs native format to BBP.' )

l_parser.add_argument( '-i', '--input',
                       dest     = 'input',
                       required = True,
                       type     = str,
                       help     = 'Input seismogram in EDGEs format.' )

l_parser.add_argument( '-l', '--locations',
                       dest = 'locations',
                       required = False,
                       type = str,
                       help = 'White space separated file containing station locations in format \'name lat lon\'' )

l_parser.add_argument( '-c', '--columns',
                       dest = 'columns',
                       required = False,
                       type = str,
                       nargs = 3,
                       default = ['Q6_C0', 'Q7_C0', 'Q8_C0'],
                       help = 'Columns which are converted' )

l_parser.add_argument( '-t', '--time_shift',
                       dest = 'time_shift',
                       required = False,
                       type = float,
                       default = 0.0,
                       help = 'Time shift applied to the seismogram. E.g., if set to 0.3 all entries smaller than 0.3 are removed with time 0.0 being the old 0.3.' )

l_parser.add_argument( '-s', '--scale',
                       dest = 'scale',
                       required = False,
                       type = float,
                       default = 100.0,
                       help = 'Scales the seismograms by the given value, e.g., for m/s to cm/s by using a factor of 100.' )

l_parser.add_argument( '-o', '--output',
                       dest     = 'output',
                       required = True,
                       type     = str,
                       help     = 'Output seismogram in BBP format.' )

l_args = vars( l_parser.parse_args() )

logging.info('parsing input: ' + l_args['input'] )
# read and store header
l_header= []
l_name = 'unknown'
with open( l_args['input'], 'r' ) as l_fileIn:
  l_lines = l_fileIn.readlines()
  for l_li in l_lines:
    if( l_li[0] == '#' ):
      l_header.append( l_li )

      if( 'receiver name:' in l_li ):
        l_name = l_li.split(': ')[1].strip()
        logging.info( 'found receiver name: ' + l_name  )

# derive lon/lat location of the receiver
l_lon = 'unknown'
l_lat = 'unknown'
if( l_args['locations'] ):
  l_locations = pandas.read_csv( l_args['locations'],
                                 sep = ' ',
                                 names = ['name', 'lat', 'lon'],
                                 dtype = str )
  l_locations = l_locations[ l_locations['name'] == l_name ]
  assert( len(l_locations) <= 1 )
  if( len(l_locations) == 1 ):
    l_lon = l_locations['lon'].iloc[0]
    l_lat = l_locations['lat'].iloc[0]
    logging.info( 'found lon / lat: ' + l_lon + ' / ' + l_lat )

# read data
l_data = pandas.read_csv( l_args['input'],
                          comment = '#',
                          float_precision = "high" )

# shift seismogram in time
l_data = l_data[ l_data['time'] >= l_args['time_shift'] ]
l_data['time'] -= l_args['time_shift']

# scale the particle velocities
for l_co in l_args['columns']:
 l_data[l_co] *= l_args['scale']

# open output file
logging.info( 'writing output: ' + l_args['output'] )
l_file = open( l_args['output'], 'w' )
l_file.write( '#     Station= ' + l_name.split('_')[1] + '\n' )
l_file.write( '#        time= 00/00/00,0:0:0.0 UTC\n' )
l_file.write( '#         lon= ' + l_lon + '\n')
l_file.write( '#         lat= ' + l_lat + '\n' )
l_file.write( '#       units= cm/s\n' )
l_file.write( '# orientation= 0.0,90.0,up\n' )
l_file.write( '#\n' )
l_file.write( '# Data fields are TAB-separated\n' )
l_file.write( '# Column 1: Time (s)\n' )
l_file.write( '# Column 2: H1 component ground velocity (+ is 0.0)\n' )
l_file.write( '# Column 3: H2 component ground velocity (+ is 90.0)\n' )
l_file.write( '# Column 4: V component ground velocity (+ is up)\n' )
l_file.write( '#\n' )
l_file.write( '###\n' )
for l_he in l_header:
  l_file.write( l_he )
l_file.write( '###\n' )

# write data
l_data.to_csv( l_file,
               sep = ' ',
               columns = ['time'] + l_args['columns'],
               header = False,
               index = False,
               float_format = '%5.6e' )

logging.info( 'finished, exiting now' )
