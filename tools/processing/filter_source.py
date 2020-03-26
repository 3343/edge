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
# Filters a given kinematic source.
##
import logging
import argparse
import h5py
import numpy
import obspy.signal.filter

# setup logger
logging.basicConfig( level=logging.INFO,
                     format='%(asctime)s - %(name)s - %(levelname)s - %(message)s' )

# command line arguments
l_parser = argparse.ArgumentParser( formatter_class=argparse.RawTextHelpFormatter,
                                    description='Filters the source data.' )

l_parser.add_argument( '-i', '--input',
                       dest     = 'input',
                       required = True,
                       type     = str,
                       help     = 'Input HDF5-file.' )

l_parser.add_argument( '-d', '--delay',
                       dest     = 'delay',
                       required = True,
                       type     = float,
                       help     = 'Delay of the point sources in seconds.' )

l_parser.add_argument( '-f', '--low_pass_frequency',
                       dest     = 'low_pass',
                       required = True,
                       type     = float,
                       help     = 'Low-pass frequency.' )

l_parser.add_argument( '-s', '--n_samples',
                       dest     = 'n_samples',
                       required = True,
                       type     = int,
                       help     = 'Number of output samples per point source.' )

l_parser.add_argument( '-o', '--output',
                       dest     = 'output',
                       required = True,
                       type     = str,
                       help     = 'Output HDF5-file with filtered time-series.' )

l_args = vars( l_parser.parse_args() )

logging.info('opening input: ' + l_args['input'] )
l_input = h5py.File( l_args['input'], 'r' )
logging.info('opening output: ' + l_args['output'] )
l_output = h5py.File( l_args['output'], 'w' )

# iterate over all points and adjust accordingly
l_nPts = len( l_input['time_parameters'] )
assert( len(l_input['time_pointers']) == l_nPts+1 )

logging.info( 'adjusting time series of ' + str(l_nPts) + ' point sources using a low-pass filter of ' + str(l_args['low_pass']) + ' Hz' )

# init time series data
l_output['time_series'] = numpy.zeros( l_nPts*l_args['n_samples'], dtype='float32' )

# write filtered time series data
for l_pt in range(l_nPts):
  # derive frequency
  l_freq = 1.0 / l_input['time_parameters'][l_pt, 1]
  l_freq = int(l_freq+0.5)

  # derive delay in number of samples
  l_nDeSas = int(l_freq * l_args['delay'] + 0.5)

  # get data which is filtered
  l_idIn0 = int(l_input['time_pointers'][l_pt])
  l_idIn1 = int(l_input['time_pointers'][l_pt+1])

  l_filtered = numpy.zeros( l_args['n_samples'], dtype='float32' )
  l_filtered[l_nDeSas:l_nDeSas+l_idIn1-l_idIn0][:] = l_input['time_series'][l_idIn0:l_idIn1]

  # do a sanity check on the filtered source
  l_absMax = numpy.max( numpy.abs( l_filtered ) )
  assert( abs(l_filtered[ 0]) < l_absMax * 0.01 )
  assert( abs(l_filtered[-1]) < l_absMax * 0.01 )

  # filter data
  l_filtered = obspy.signal.filter.lowpass( data      = l_filtered, 
                                            freq      = l_args['low_pass'],
                                            corners   = 4,
                                            df        = l_freq,
                                            zerophase = True )

  # write data back
  l_idOut0 = l_pt*l_args['n_samples']
  l_idOut1 = (l_pt+1)*l_args['n_samples']
  l_output['time_series'][l_idOut0:l_idOut1] = l_filtered

# set remaining data
l_output['time_pointers'] = numpy.array( range( 0,
                                                (l_nPts+1)*l_args['n_samples'],
                                                l_args['n_samples'] ), dtype='uint64' )

l_output['points'] = l_input['points'][:]
l_output['scalings'] = l_input['scalings'][:]
l_output['time_parameters'] = l_input['time_parameters'][:]

logging.info('finished writing ' + l_args['output'] )
