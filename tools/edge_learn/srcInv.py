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
# Learns about kinematic sources models.
##
import logging
l_logo =["##########################################################################",
         "##############   ##############            ###############  ##############",
         "##############   ###############         ################   ##############",
         "#####            #####       #####      ######                       #####",
         "#####            #####        #####    #####                         #####",
         "#############    #####         #####  #####                  #############",
         "#############    #####         #####  #####      #########   #############",
         "#####            #####         #####  #####      #########           #####",
         "#####            #####        #####    #####        ######           #####",
         "#####            #####       #####      #####       #####            #####",
         "###############  ###############         ###############   ###############",
         "###############  ##############           #############    ###############",
         "#####################################################################learn"];

# setup logger
logging.basicConfig( level=logging.INFO,
                     format='%(asctime)s - %(name)s - %(levelname)s - %(message)s' )

# print our logo
for l_li in l_logo:
  logging.info( l_li )

import argparse
import h5py
import numpy
import pandas
import keras
import time
import tensorflow
import edge_learn.io.Synthetics
import edge_learn.io.Ruptures
import edge_learn.models.Rupture

# command line arguments
l_parser = argparse.ArgumentParser( formatter_class=argparse.RawTextHelpFormatter,
                                    description='Visualizes the training data.' )

l_parser.add_argument( '-i', '--seismograms_hdf5',
                       dest     = 'seismograms_hdf5',
                       required = False,
                       type     = str,
                       help     = 'HDF5 file, containing the synthetics' )

l_parser.add_argument( '-d', '--seismograms_ds',
                       dest     = 'seismograms_ds',
                       required = False,
                       default  = 'synthetics',
                       type     = str,
                       help     = 'Name of the synthetic data set.' )

l_parser.add_argument( '-s', '--slice',
                       dest     = 'slice',
                       required = False,
                       type     = int,
                       nargs    = 3,
                       default  = [0,0,1],
                       help     = 'Slices the seismograms before usage. First: start, second: relative end (offset from back), third: step. Default: [0,1,0]' )

l_parser.add_argument( '-r', '--rups_hdf5',
                       dest     = 'rups_hdf5',
                       required = False,
                       type     = str,
                       help     = 'HDF5 file, containing the rupture descriptions.' )

l_parser.add_argument( '-q', '--rups_ds',
                       dest     = 'rups_ds',
                       required = False,
                       nargs    = '+',
                       default  = ['rupture_time', 'rise_time', 'total_slip'],
                       type     = str,
                       help     = 'Name of the rupture model data sets.' )

l_parser.add_argument( '-b', '--batch_size',
                       dest     = 'batch_size',
                       required = False,
                       type     = int,
                       default  = 16,
                       help     = 'Used batch size in the training.')

l_parser.add_argument( '-e', '--epochs',
                       dest     = 'epochs',
                       required = False,
                       type     = int,
                       default  = 100,
                       help     = 'Number of epochs used in the training.')

l_parser.add_argument( '-t', '--val_split',
                       dest     = 'val_split',
                       required = False,
                       type     = float,
                       default  = 0.2,
                       help     = 'Relative amount of data in [0,1) assigned to validation data.')

l_parser.add_argument( '-n', '--noise',
                       dest     = 'noise',
                       required = False,
                       type     = float,
                       default  = 0.0,
                       help     = 'Relative noise at to each data point (scaled by 1.0+p*n, where p is normally distributed in [-1/2, 1/2] and n the given noise value.' )

l_parser.add_argument( '-p', '--plot_model',
                       dest     = 'plot_model',
                       required = False,
                       type     = str,
                       default  = None,
                       help     = 'Plots the Keras model to the given location.' )

l_parser.add_argument( '-j', '--in_model',
                       dest     = 'in_model',
                       required = False,
                       type     = str,
                       default  = None,
                       help     = 'Path from which the trained Keras model will be read (no training will be performed).' )

l_parser.add_argument( '-k', '--out_model',
                       dest     = 'out_model',
                       required = False,
                       type     = str,
                       default  = None,
                       help     = 'Output path to which the trained Keras model will be written.' )


l_args = vars( l_parser.parse_args() )

if( l_args['in_model'] == None ):
  # read the labels
  logging.info( 'reading rupture models' )
  l_rups = edge_learn.io.Ruptures.Ruptures( l_args['rups_hdf5'],
                                            l_args['rups_ds'] )
  logging.info( 'shape of the ruptures: ' + str( l_rups.getData().shape)  )
else:
  l_rups = None

# read the synthetics
logging.info( 'reading synthetic seismograms' )
l_seismograms = edge_learn.io.Synthetics.Synthetics( l_args['seismograms_hdf5'],
                                                     l_args['seismograms_ds'] )

# slice the synthetics
logging.info( 'slicing synthetic seismograms' )
l_start = l_args['slice'][0]
l_end   = l_seismograms.m_data.shape[-1] - l_args['slice'][1]
l_step  = l_args['slice'][2]
l_seismograms.slice( l_start,
                     l_end,
                     l_step )
logging.info( 'shape of the sliced synthetics: ' + str( l_seismograms.getData().shape)  )


if( l_args['in_model'] != None ):
  logging.info( 'loading trained model: ' + l_args['in_model'] )

l_rupInv = edge_learn.models.Rupture.Rupture( l_args['val_split'],
                                              l_args['batch_size'],
                                              l_seismograms,
                                              l_rups,
                                              l_args['in_model'] )

if( l_args['in_model'] == None ):
  logging.info( 'compiling the model' )
  l_rupInv.compile()

  if( l_args['plot_model'] != None ):
    logging.info( 'plotting the model' )
    l_rupInv.plot( l_args['plot_model'] )

  logging.info( 'training the model' )
  l_rupInv.train( l_args['epochs'] )


if( l_args['out_model'] != None ):
  logging.info( 'saving trained model to: ' + l_args['out_model'] )
  l_rupInv.save( l_args['out_model'] )