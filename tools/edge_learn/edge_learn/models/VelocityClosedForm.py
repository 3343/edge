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
# Performs an inversion for the velocity model, which has a closed form description.
##
import logging
import edge_learn.io.Config
import edge_learn.io.Synthetics
import edge_learn.augment.Single
import keras

class VelocityClosedForm:
  def __init__( self,
                i_split,
                i_batchSize,
                i_config,
                io_synths,
                i_noise=0 ):
    """Initializes the closed form velocity inversion.
      The initialization will modify the underlying data of the synthetics for memory efficiency.
      If this is not intended, pass a copy.
    
    Arguments:
      i_split {float} -- split in training and validation data.
      i_batchSize {int} -- size of the batches.
      i_config {io.Config} -- Config of the forward simulations.
      io_synths {io.Synthetics} -- Synthetics of the forward simulations.
      i_noise relative noise, which is added by scaling the true numerical value with (1+s*p), where s is the given scaling and p normally distributed in [1/2,-1/2]
    """
    logging.info( 'initializing CNN data for velocity inversion' )

    # reorder data, "image" is 3D: sim, recv, time
    io_synths.transpose( (0, 3, 1, 2, 4) )
    self.m_shape = io_synths.getData().shape

    logging.info( 'shape of the CNN data: ' + str(self.m_shape) )

    # generate training and validation data generators
    assert( i_split > 0 and i_split <= 1 )
    l_split = int( i_split * self.m_shape[0] )

    logging.info( 'using ' + str(l_split) + ' samples as training data and ' + str(self.m_shape[0]-l_split) + ' samples for validation'  )

    self.m_genTra = edge_learn.augment.Single.Single( i_batchSize,
                                                      io_synths.getData()[0:l_split],
                                                      i_config.getData()[0:l_split],
                                                      True,
                                                      i_noise )

    self.m_genVal = edge_learn.augment.Single.Single( i_batchSize,
                                                      io_synths.getData()[l_split:-1],
                                                      i_config.getData()[l_split:-1],
                                                      False )

  def compile( self ):
    """Compiles the CNN.
    """
    self.m_model = keras.Sequential()
    # perform convolutions in time
    for l_co in range(5):
      self.m_model.add( keras.layers.Conv3D( input_shape=( self.m_shape[1:] if l_co == 0 else []),
                                             filters=16,
                                             kernel_size=(1, 1, 7),
                                             data_format='channels_first',
                                             activation='relu' ) )
      self.m_model.add( keras.layers.MaxPooling3D( pool_size=(1, 1, 2),
                                                   data_format='channels_first' ) )

    # perform convolution over receivers and source
    for _ in range(2):
      self.m_model.add( keras.layers.Conv3D( filters=8,
                                             kernel_size=(1, 7, 1),
                                             data_format='channels_first',
                                             activation='relu' ) )
      self.m_model.add( keras.layers.MaxPooling3D( pool_size=(1, 2, 1),
                                                   data_format='channels_first' ) )

    for _ in range(2):
      self.m_model.add( keras.layers.Conv3D( filters=8,
                                             kernel_size=(3, 1, 1),
                                             data_format='channels_first',
                                             activation='relu' ) )
      self.m_model.add( keras.layers.MaxPooling3D( pool_size=(2, 1, 1),
                                                   data_format='channels_first' ) )

    self.m_model.add( keras.layers.Flatten() )

    self.m_model.add( keras.layers.Dense(2) )

    self.m_model.compile( loss='mean_squared_error',
                          metrics=['mean_absolute_error'],
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

    self.m_model.fit_generator( generator=self.m_genTra,
                                validation_data=self.m_genVal,
                                epochs=i_epochs,
                                verbose=1,
                                use_multiprocessing=True,
                                callbacks=l_cbs )

  def predictVal( self ):
    """Simply predicts the validation data.
    Returns:
      {numpy array} -- predicted labels.
    """
    return self.m_model.predict_generator( self.m_genVal )