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
# Deconvolutions and convolutions.
import h5py
import numpy

def convSrc( i_nSynth,
             i_srcNew,
             i_dec ):
  """Convolves the deconvolved signal with a new source
  
  Arguments:
    i_nSynth number of entries in the generated receivers, the convolved signal will be shortened/zero-passed accordingly.
    i_srcNew {array of float} -- source signal, used in the convolution.
    i_dec {arary of float} -- deconvolved signal.


  Returns:
    [array of float] -- Convolved signal.
  """

  # assemble source sampling, including zero padding
  l_src = i_srcNew
  if( len(l_src)%2 == 0 ):
    l_src = numpy.append( l_src, [0] )

  # convolve with new source
  l_conv = numpy.convolve( i_dec, i_srcNew, mode='same' )

  # pad zeros to account for the time shift
  l_conv = numpy.append( numpy.zeros( int( len(i_srcNew)/2 ) ), l_conv )

  # get rid of overhead data or add zero
  l_up = i_nSynth - len(l_conv)

  # get rid of overhead data or add zeros
  if( l_up < 0 ):
    l_conv = l_conv[0:l_up]
  elif( l_up > 0 ):
    l_conv = numpy.append( l_conv, numpy.zeros(l_up) )

  return l_conv

class DeConvolve:
  def __init__( self,
                i_path='',
                i_dataSet='deconv' ):
    """Initializes the deconvolved data by reading the respective HDF5 files.
    
    Arguments:
      i_path {string} -- Path to the data (optional).
      i_dataSet {string} -- Dataset in the HDF5 file(s), which will be read.
    """
    if( i_path != '' ):
      l_hdf5 = h5py.File( i_path, 'r' )
      self.m_data = numpy.array( l_hdf5[ i_dataSet ] )
      l_hdf5.close()

  def write( self,
             i_path,
             i_dataSet='deconv' ):
    """Writes the deconvolved data.
    
    Keyword Arguments:
      i_path {string} -- Output path.
      i_dataSet {string} -- Data set within the HDF5 file (default: {'deconv'}).
    """
    with h5py.File( i_path, 'w') as l_h5:
      l_h5.create_dataset( i_dataSet,
                            self.m_data.shape,
                            dtype='float32' )
      l_h5[i_dataSet][:] = self.m_data[:]
      l_h5.close()

  def getData( self ):
    """Gets the raw data of the deconvolved signals.
    
    Returns:
      [array of float] -- [raw data]
    """
    return self.m_data

  def deconvSrc( self,
                 i_src,
                 i_data ):
    """Deconvolves the synthetics with the given source sampling.
    
    Arguments:
      i_src {array of float} -- Sampling of the source in the same frequency as the synthetics.
      i_data {array of float} -- Raw data of the synthetics.
    """

    import copy
    import scipy.signal

    # store src sampling locally
    l_src = copy.deepcopy(i_src)

    # produce a symmetric stencil
    if( len(i_src)%2 == 0 ):
      l_src = numpy.append( l_src, [0] )

    # alloc memory for the deconvolutions
    self.m_data = numpy.zeros( list(i_data.shape[:-1]) + [ i_data.shape[-1] - len(l_src) + 1  ] )

    # iterate over synthetics
    for l_id in numpy.ndindex( i_data.shape[0:-1] ):
      # deconvolve the signal
      self.m_data[l_id], _ = scipy.signal.deconvolve( i_data[l_id], l_src )