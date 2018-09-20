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
# Representation of rupture models.
##
import h5py
import copy
import numpy

class Ruptures:
  def __init__( self,
                i_path,
                i_dataSets=['rupture_time', 'rise_time', 'total_slip'] ):
    """Initializes the rupture model by reading the respective HDF5 file.
    
    Arguments:
      i_path {string} -- Path to HDF5 file.
      i_dataSets {string} -- Datasets in the HDF5 file(s), which will be read.
    """
    assert( len(i_dataSets) > 0 )

    l_hdf5 = h5py.File( i_path, 'r' )

    # store quantities
    self.m_qts = copy.deepcopy( i_dataSets )

    # assemble shape of the data
    l_shape = [ len(self.m_qts) ] + list( l_hdf5[ self.m_qts[0] ].shape )

    # allocate memory
    self.m_data = numpy.zeros( l_shape, dtype='float32' )

    # read data
    for l_qt in range( len(self.m_qts) ):
      # ensure that all quantities have the same shape
      assert( list(l_hdf5[ self.m_qts[l_qt] ].shape) == l_shape[1:] )

      # read the data
      self.m_data[l_qt] = l_hdf5[ self.m_qts[l_qt] ][:]
    
    # store samples as slowest dim, not quantities
    self.m_data = numpy.swapaxes( self.m_data, 0, 1 )


  def getData( self ):
    """Gets the raw data of the ruptures.
    
    Returns:
      [array of float] -- [raw data]
    """
    return self.m_data