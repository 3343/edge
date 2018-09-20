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
# Representation of raw synthetics.
##
import os
import h5py
import numpy
import functools

class Synthetics:
  def __init__( self,
                i_path,
                i_dataSet='synthetics' ):
    """Initializes the synthetic data type by reading the respective HDF5 files.
    
    Arguments:
      i_path {string} -- Path to the directory or file, containing the batches.
                         Directory: All files named synthetics.hdf5 in numeric subdirectories will be read.
      i_dataSet {string} -- Dataset in the HDF5 file(s), which will be read.
    """
    if os.path.isdir( i_path ):
      # assemble sub-directories
      l_subDirs = []
      for l_di in os.listdir( i_path ):
        try:
          # check if this is an integer
          l_subDirs.append( int( l_di ) )
        except ValueError:
          # nothing here
          continue

      # sort the sub-directories
      l_subDirs.sort()

      # ensure that everything is there
      for l_sd in range( len(l_subDirs[1:] ) ):
        assert( l_sd == l_subDirs[l_sd] )

      # assemble paths to HDF5 files
      l_files = [ i_path + '/' + str(l_sd) + '/synthetics.hdf5' for l_sd in l_subDirs ]

      # read the data
      self.m_data = []
      for l_fi in l_files:
        l_hdf5 = h5py.File( l_fi, 'r' )
        self.m_data.append( l_hdf5[ i_dataSet ] )

      # convert to numpy array
      self.m_data = numpy.concatenate( self.m_data )
    else:
      l_hdf5 = h5py.File( i_path, 'r' )
      self.m_data = numpy.array( l_hdf5[ i_dataSet ] )
      l_hdf5.close()

  def getData( self ):
    """Gets the raw data of the synthetics.
    
    Returns:
      [array of float] -- [raw data]
    """
    return self.m_data

  def clearData( self ):
    """Deletes the raw data.
    """
    self.m_data = None

  def write( self,
            i_path,
            i_dataSet='synthetics' ):
    """Writes the synthetics.
    
    Keyword Arguments:
      i_path {string} -- Output path.
      i_dataSet {string} -- Data set within the HDF5 file.
    """
    with h5py.File( i_path, 'w') as l_h5:
      l_h5.create_dataset( i_dataSet,
                            self.m_data.shape,
                            dtype='float32' )
      l_h5[i_dataSet][:] = self.m_data[:]
      l_h5.close()

  def slice( self,
             i_start,
             i_end,
             i_step ):
    """Slices the synthetics by extracting the values (fastest dimension) from start to end by the given sampling.
    
    Arguments:
      i_start {integer} -- First value in the fast dimension, which is extracted.
      i_size {integer} -- Number of values, which are extracted.
                          If this does not match a multiple of step size, values will be missing.
      i_step {integer} -- Step size, used in the extraction of values.
    """
    # reshape the data to a simple 2D format
    l_shape = self.m_data.shape
    self.m_data = self.m_data.reshape( (functools.reduce(lambda x, y: x*y, l_shape[0:-1]), l_shape[-1]) )
    # slice the data
    self.m_data = self.m_data[ :, i_start:i_end:i_step ]

    # transfer back
    self.m_data = self.m_data.reshape( list(l_shape[0:-1]) + [ self.m_data.shape[-1] ] )

  def transpose( self, i_perm ):
    """Transposes the synthetics.
    
    Arguments:
      i_perm {tuple of integers} Permutation of the dimensions.
    """
    self.m_data = numpy.transpose( self.m_data, i_perm )

  def lowpass( self,
               i_sampling,
               i_filter ):
    """Low pass filters the synthetics
    
    Arguments:
      i_sampling {float} -- Sampling frequency of the synthetics, used in the filter.
      i_filter {float} -- Lowpass filter frequency.
    """
    import obspy.signal.filter

    # iterate over synthetics and apply filter
    for l_id in numpy.ndindex( self.m_data.shape[0:-1] ):
      self.m_data[l_id] = obspy.signal.filter.lowpass( self.m_data[l_id],  i_filter, i_sampling )
