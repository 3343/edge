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
# Converts a kinematic source description in NRF-format to EDGE's HDF5 format.
##
import logging
import argparse
import netCDF4
import h5py
import numpy

# setup logger
logging.basicConfig( level=logging.INFO,
                     format='%(asctime)s - %(name)s - %(levelname)s - %(message)s' )

# command line arguments
l_parser = argparse.ArgumentParser( formatter_class=argparse.RawTextHelpFormatter,
                                    description='Visualizes the training data.' )

l_parser.add_argument( '-i', '--input_nrf',
                       dest     = 'input_nrf',
                       required = True,
                       type     = str,
                       help     = 'Input NRF-file.' )

l_parser.add_argument( '-o', '--output_hdf5',
                       dest     = 'output_hdf5',
                       required = True,
                       type     = str,
                       help     = 'Output HDF5-file.' )

l_parser.add_argument( '-m', '--mu',
                       dest     = 'mu',
                       required = False,
                       type     = float,
                       help     = 'Parameter mu, which is used if non-positive values are found in NRF.' )

l_parser.add_argument( '-l', '--lambda',
                       dest     = 'lambda',
                       required = False,
                       type     = float,
                       default  = 0,
                       help     = 'Parameter lambda, which is used (defaults to zero).' )

l_args = vars( l_parser.parse_args() )

logging.info( 'reading: ' + l_args['input_nrf'] )

l_ncDs = netCDF4.Dataset( l_args['input_nrf'] )

# determine the number of dimensions
if( 'sliprates3' in l_ncDs.variables ):
  l_nDis = 3
else:
  l_nDis = 2
logging.info( 'assuming ' + str(l_nDis) + ' dimensions in the NetCDF-file.' )

# read the centres from NetCDF
l_nPts = len( l_ncDs.variables['centres'] )

# init points
l_points = numpy.zeros( (l_nPts, l_nDis), dtype='float32' )

logging.info( 'reading ' + str(l_nPts) + ' NetCDF centres' )
l_ncCens = l_ncDs.variables['centres'][:]
logging.info( 'storing points' )
for l_pt in range(l_nPts):
  for l_di in range(l_nDis):
    l_points[l_pt][l_di] = l_ncCens[l_pt][l_di]

# read subfaults from NetCDF
logging.info( 'reading ' + str(l_nPts) + ' subfaults' )
l_ncSf = l_ncDs.variables['subfaults']
assert( len(l_ncSf) == l_nPts)

# init time parameters
logging.info( 'storing time parameters' )
l_timePars = numpy.zeros( (l_nPts, 2), dtype='float32' )
for l_pt in range(l_nPts):
  l_timePars[l_pt][0] = l_ncSf[l_pt][0]
  l_timePars[l_pt][1] = l_ncSf[l_pt][1]

# read the NetCDF time pointers
logging.info( 'reading sroffsets' )
l_ncPtrs = l_ncDs.variables['sroffsets']

assert( len(l_ncPtrs) == l_nPts + 1 )

# story the pointers
logging.info( 'storing time pointers' )
l_timePointers = numpy.zeros( l_nPts+1, dtype='uint64' )
for l_pt in range(l_nPts+1):
  l_timePointers[l_pt] = l_ncPtrs[l_pt][0]

# read sliprate smaples
logging.info( 'reading NetCDF sliprates' )
l_ncSr = l_ncDs.variables['sliprates1']
assert( len(l_ncSr) == l_timePointers[-1] )

# store time series
logging.info( 'storing time_series' )
l_timeSeries = numpy.zeros( l_timePointers[-1], dtype='float32' )
l_timeSeries[:] = l_ncSr[:]

# derive scalings
logging.info( 'deriving scalings' )
if( l_nDis == 2 ): l_nQts = 5
else:              l_nQts = 9
l_scalings = numpy.zeros( (l_nPts, l_nQts), dtype='float32' )

##
# Slip coefficients:
#   delta u * [ lambda * l_k * n_k * delta_{ij} + mu * ( l_i * n_j + l_j * n_i )
# 
#   delta u: slip cofficient is the term multiplied with it.
#   lambda: Lame parameter lambda
#   delta_{ij}: Kronecker delta
#   mu: Lame parameter mu
#   l: unit vector in slip-direction  l = ( l_1, l_2 ) or l = ( l_1, l_2, l_3 )
#   n: unit normal n = (n_1, n_2) or n = (n_1, n_2, n_3)
#  
#   Reference:
#     Source Mechanisms of Earthquakes, 2014
#     Theory and Practice
#     Udias, Madariaga, Buforn
#     Eq. (5.26)
##
for l_pt in range(l_nPts):
  # read area and mu
  l_A =  l_ncSf[l_pt][3]
  l_lam = float(l_args['lambda'])
  l_mu = l_ncSf[l_pt][2]
  if( l_mu <= 0 ):
    if( l_args['mu'] != None ):
      l_mu = l_args['mu']
    else:
      logging.error( 'found negative mu in NetCDF file for point ' + str(l_pt) + ', aborting' )
      exit(1)

  # reads slip directions
  l_sds = numpy.zeros( (l_nDis, l_nDis) )
  for l_d0 in range(l_nDis):
    for l_d1 in range(l_nDis):
      # normal goes first, but is last in NetCDF-file
      if( l_d0 == 0 ):
        l_sds[l_d0][l_d1] = l_ncSf[l_pt][-1][l_d1]
      elif( l_d0 == 1 ):
        l_sds[l_d0][l_d1] = l_ncSf[l_pt][-3 + 3 - l_nDis][l_d1]
      elif( l_d0 == 2 ):
        l_sds[l_d0][l_d1] = l_ncSf[l_pt][-2][l_d1]
        
      else:
        l_sds[l_d0][l_d1] = l_ncSf[l_pt][-l_nDis + l_d0][l_d1]

  # iterate over slip directions
  for l_sd in range(l_nDis):

    # scalar product with normal
    l_spN = 0
    for l_di in range(l_nDis):
      l_spN += l_sds[l_sd][l_di]*l_sds[0][l_di]

    # contribution if slip is in normal direction
    for l_di in range(l_nDis):
      l_scalings[l_pt][l_di] += l_A * l_lam * l_spN

    # normal stress
    l_scalings[l_pt][l_di] += l_A * 2.0 * l_mu * l_sds[l_sd][l_di] * l_sds[0][l_di]

    # moment tensor indices
    l_mId = [ [0,1], [1,2], [0,2] ]

    # shear stresses
    if( l_nDis == 2 ):
      l_nShear = 1
    else:
      l_nShear = 3

    for l_sh in range(l_nShear):
      l_m1 = l_mId[l_sh][0]
      l_m2 = l_mId[l_sh][1]

      l_scalings[l_pt][l_di] += l_A * l_mu * l_sds[l_sd][l_m1] * l_sds[0][l_m2] + \
                                l_A * l_mu * l_sds[l_sd][l_m2] * l_sds[0][l_m1]

# write the HDF5
logging.info( 'writing HDF5 output to ' + l_args['output_hdf5'] )

l_hdf5 = h5py.File( l_args['output_hdf5'], 'w' )

l_hdf5['points'] = l_points
l_hdf5['scalings'] = l_scalings
l_hdf5['time_parameters'] = l_timePars
l_hdf5['time_pointers'] = l_timePointers
l_hdf5['time_series'] = l_timeSeries

logging.info( 'finished' )