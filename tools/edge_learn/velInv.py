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
# Learns about velocity models.
##
import argparse
import logging
import h5py
import numpy
import pandas
import keras
import time
import tensorflow
import edge_learn.io.Config
import edge_learn.io.Synthetics
import edge_learn.models.VelocityClosedForm

# command line arguments
l_parser = argparse.ArgumentParser( formatter_class=argparse.RawTextHelpFormatter,
                                    description='Visualizes the training data.' )

l_parser.add_argument( '-c', '--config',
                       dest     = 'config',
                       required = True,
                       type     = str,
                       help     = 'Configurations of the forward runs.' )

l_parser.add_argument( '-d', '--data_set',
                       dest     = 'data_set',
                       required = False,
                       default  = 'synthetics',
                       type     = str,
                       help     = 'Name of the data set.' )

l_parser.add_argument( '-5', '--dir_hdf5',
                       dest     = 'dir_hdf5',
                       required = True,
                       type     = str,
                       help     = 'Path of the directories containing the HDF5 file (split by batches).' )

l_parser.add_argument( '-s', '--slice',
                       dest     = 'slice',
                       required = False,
                       type     = int,
                       nargs    = 3,
                       default  = [0,0,1],
                       help     = 'Slices the synthetics before traning. First: start, second: relative end (offset from back), third: step. Default: [0,1,0]' )

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

l_parser.add_argument( '-t', '--split',
                       dest     = 'split',
                       required = False,
                       type     = float,
                       default  = 0.8,
                       help     = 'Split between training data and validation data.')

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

l_args = vars( l_parser.parse_args() )

# setup logger
logging.basicConfig( level=logging.DEBUG,
                     format='%(asctime)s - %(name)s - %(levelname)s - %(message)s' )

# read the CSV config
l_config = edge_learn.io.Config.Config( l_args['config'],
                                        [ 'fault_angle_01', 'fault_offset_01' ] )

# read the synthetics
l_synthetics = edge_learn.io.Synthetics.Synthetics( l_args['dir_hdf5'], l_args['data_set'] )

# slice the synthetics
l_start = l_args['slice'][0]
l_end   = l_synthetics.m_data.shape[-1] - l_args['slice'][1]
l_step  = l_args['slice'][2]
l_synthetics.slice( l_start,
                    l_end,
                    l_step )


logging.info( 'shape of the sliced training data: ' + str( l_synthetics.getData().shape)  )


l_veInv = edge_learn.models.VelocityClosedForm.VelocityClosedForm( l_args['split'],
                                                                   l_args['batch_size'],
                                                                   l_config,
                                                                   l_synthetics,
                                                                   l_args['noise'] )

logging.info( 'compiling the model' )
l_veInv.compile()
if( l_args['plot_model'] != None ):
  logging.info( 'plotting the model' )
  l_veInv.plot( l_args['plot_model'] )

logging.info( 'training the model' )
l_veInv.train( l_args['epochs'] )

logging.info( 'printing model, applied to training data' )
print( l_veInv.predictVal() )