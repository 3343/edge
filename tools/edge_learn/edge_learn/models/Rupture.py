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
# Performs an inversion for the rupture process.
##
import numpy
import keras
import copy

class Rupture:
  def __init__( self,
                i_valSplit=None,
                i_batchSize=None,
                i_synths=None,
                i_ruptures=None,
                i_path=None ):
    """Initializes the model for the rupture inversion.
       Input arguments are either (i_valSplit, i_batchSize, i_synths and i_ruptures) or (i_path).
    
    Arguments:
      i_valSplit {float} -- Relative amount of data in [0,1) assigned to validation data.
      i_batchSize {int} -- Used batch size.
      i_synths {io.Synthetics} -- Synthetics of the forward simulations.
      i_ruptures {io.Ruptures} -- Ruptures, associated to the synthetics.
      i_path {string} -- path to a computed Keras model, stored as HDF5.
    """
    # ensure valid arguments
    assert( (i_valSplit != None and i_batchSize != None and i_synths != None and i_ruptures != None) or i_path != None )

    # ensure a valid split
    assert( i_valSplit >= 0 and i_valSplit < 1 )

    if i_path == None:
      # transpose the synthetics (velocities/channels as second slowest dim)
      l_x = copy.deepcopy( i_synths.getData() )
      l_x = l_x.transpose( 0, 2, 1, 3 )
      # reference labels
      l_y = i_ruptures.getData()

      # store the split
      self.m_split = i_valSplit

      self.m_samples = {}
      self.m_samples['x'] = l_x
      self.m_samples['y'] = l_y

      # store the batch size
      self.m_batchSize = i_batchSize
    else: self.load( i_path )

  def load( self,
            i_path ):
    """Load the Keras model from the given HDF5 file.
    
    Arguments:
      i_path {string} -- path to HDF5 file, containing the Keras model.
    """
    self.m_model = keras.models.load_model( i_path )

  def save( self,
            i_path ):
    """Saves the underlying Keras model to the given HDF5 file.
    
    Arguments:
      i_path {string} -- Path to HDF5 file
    """
    self.m_model.save( i_path )

  def compile( self ):
    """Compiles the CNN.
    """
    self.m_model = keras.Sequential()

    # perform convolutions in time
    for l_co in range(5):
      self.m_model.add( keras.layers.Conv2D( input_shape=( self.m_samples['x'].shape[1:] if l_co == 0 else []),
                                             filters=2**(l_co+1),
                                             kernel_size=(1, 7),
                                             data_format='channels_first',
                                             activation='relu' ) )
      self.m_model.add( keras.layers.MaxPooling2D( pool_size=(1, 2),
                                                   data_format='channels_first' ) )

    # perform a single convolution in space
    self.m_model.add( keras.layers.Conv2D( filters=2**6,
                                          kernel_size=(self.m_samples['x'].shape[2], 1),
                                          data_format='channels_first',
                                          activation='relu' ) )

    self.m_model.add( keras.layers.Flatten() )

    # add flattened output layers
    for _ in range(3):
      self.m_model.add( keras.layers.Dense( numpy.prod( self.m_samples['y'].shape[1:] ) ) )

    # compile the model
    self.m_model.compile( loss='mean_absolute_error',
                          metrics=['mean_squared_error'],
                          optimizer='adam' )

  def plot( self,
            i_path ):
    """Plots the compiled model.
    
    Arguments:
      i_path {string} -- Path to the output file.
    """
    keras.utils.plot_model( self.m_model,
                            i_path,
                            show_shapes=True )

  def train( self,
             i_epochs,
             i_logDir=None ):
    """Trains the CNN.
    
    Arguments:
      i_epochs {int} -- Number of epochs used in the training.
      i_logDir {string} -- Path to which the Tensorboard logs are written (optional)
    """

    # init adaptive reduction in learning rate
    l_cbs = [ keras.callbacks.ReduceLROnPlateau( monitor='val_loss',
                                                 factor=0.8,
                                                 patience=10,
                                                 verbose=1,
                                                 min_lr=1.0E-9 ) ]

    if( i_logDir != None ):
      l_cbs.append( TensorBoard( log_dir=i_logDir+'/{}'.format( 'vel_closed' ) ) )

    self.m_model.fit( x                   = self.m_samples['x'],
                      y                   = self.m_samples['y'].reshape(-1, numpy.prod( self.m_samples['y'].shape[1:] ) ),
                      validation_split    = self.m_split,
                      epochs              = i_epochs,
                      verbose             = 1,
                      callbacks           = l_cbs )


  def predict( self,
               i_synths ):
    """Applies the trained model to the given data.
    Arguments:
      i_synths {array of float} -- Synthetic seismograms.
    Returns:
      {numpy array} -- predicted labels.
    """
    l_x = copy.deepcopy( i_synths.getData() )
    l_x = l_x.transpose( 0, 2, 1, 3 )

    return self.m_model.predict( l_x )