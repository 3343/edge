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
# Data generator, which augments the single seismograms by noise and time shifts.
##
import keras.utils.data_utils
import numpy.random
import copy

class Single( keras.utils.data_utils.Sequence ):
  def __init__( self,
                i_batchSize,
                i_data,
                i_labels,
                i_shuffle,
                i_noise=0.0 ):
    """Initializes the data generated for direct mappings of data and receivers.
       Data augmentation is only done on a per-datum basis.
    
    Arguments:
      i_batchSize {int} -- Batch size.
      i_data {numpy array} -- Training and validation data, where the slowest dimension is used for the batches.
      i_labels {numpy array} -- Labels, where the slows dimension is used for the batches.
      i_shuffle {bool} -- Shuffles the dataset at construction and after every epoch.
      i_noise {float} -- adds uniform random noise of th given magnitude to all data points, when returning batch samples.
    """

    # check for sound input
    assert( len(i_data) == len(i_labels) )
    assert( i_batchSize < len(i_data) )

    self.m_data = i_data
    self.m_labels = i_labels
    self.m_batchSize = i_batchSize

    self.m_ids = list( range( len(i_data) ) )
    self.m_shuffle = i_shuffle
    self.m_noise = i_noise

  def on_epoch_end( self ):
    """Shuffles the data (if set) at the very beginning and after every epoch.
    """
    if( self.m_shuffle ):
      numpy.random.shuffle( self.m_ids )

  def __len__( self ):
    """Returns the number of samples in the data set.
    
    Returns:
      int -- number of samples in the data set.
    """

    l_len = numpy.ceil( len( self.m_data ) / float( self.m_batchSize ) )
    return int( l_len )

  def addNoise( self,
                i_data ):
    l_data = copy.deepcopy( i_data )

    # only continue for positive noise (avoids numpy errors)
    if( self.m_noise > 0.0 ):
      # derive scaling for the noise
      l_noise = numpy.random.normal( loc=1,
                                     scale=self.m_noise,
                                     size=i_data.shape )

      l_data = numpy.multiply( l_data, l_noise )

    return l_data

  def __getitem__( self,
                   i_bId ):
    """Gets the samples and labels for the batch with the given id.
    
    Arguments:
      i_bId {int} -- Id of the batch.
    
    Returns:
      two numpy arrays -- Sample data and labels.
    """

    l_ids = self.m_ids[ i_bId * self.m_batchSize :  (i_bId + 1) * self.m_batchSize ]

    l_data = copy.deepcopy( self.m_data[ l_ids ] )
    self.addNoise( l_data )

    return l_data, self.m_labels[ l_ids ]