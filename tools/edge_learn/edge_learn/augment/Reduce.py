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
# Reduces the given raw synthetics.
##
def sum( i_freq,
         i_nSrc,
         i_nSynth,
         i_scale,
         i_shift,
         i_stretch,
         i_dec ):
  """Sums the given time series data.
     Assume a given 2D reduction of shape (220, 240).
     Now, assume that a total of 30 time series in a 2D array is given: (10, 3, 1500).
     Each isolated time series has therefore 1500 samples in [0, 1500*i_freq].
     In that case, the entire input data would have shape: (220, 240, 10, 30, 1500).

     This function support three parameters: "scale", "shift", "stretch".
     Each of the parameters is given for each of the reduction entries of shape (220, 240).
     Additionally, we support the application of multiple reductions in one call.
     This reduces memory pressure, as the input data is only streamed once.
     Assume, a batch size of 16, performing 16 reductions.
     In this case, the shape of the three parameters would be (16, 220, 240).
     Analogue, the reduced output has shape (16, 10, 30)

     The parameter "scale" scales the amplitude of each time series.
     The parameter "shift" shifts the respective time series in time.
     For example, with shift=0.5, the contributions would be sampled after 0.5s.
     The parameter "stretch" stretches (or squeezes) the time series.
  
  Arguments:
    i_freq {float} -- frequency of the input data.
    i_nSynth number of entries in the generated receivers, the convolved signal will be shortened/zero-passed accordingly.
    i_nSrc {integer} -- number of entires in the original boxcar source.
    i_scale {array of float} -- scaling of the summands.
    i_shift {array of float} -- time shift of the summands.
    i_stretch {array of float} -- stretch of the summands.
    i_dec {array of float} -- deconvolved signals.
  
  Returns:
    [array of float] -- reduced data.
  """
  import numpy
  import edge_learn.augment.DeConvolve
  import scipy.signal

  l_data = i_dec.getData()

  # sanity checks on the input dimensions
  assert( i_scale.shape == i_shift.shape == i_stretch.shape )
  assert( l_data.shape[:(len(i_scale.shape)-1)] == i_scale.shape[1:] )

  # output data
  l_out = numpy.zeros( [i_scale.shape[0]]+list(l_data.shape[len(i_scale.shape)-1:-1])+[i_nSynth] )

  # the point in time of the original data 
  l_time = numpy.linspace( 0, i_nSynth*(1.0/i_freq), i_nSynth )

  # iterate over the individual contributions (excluding batches)
  for l_i0 in numpy.ndindex( i_scale.shape[1:] ):
    # points which are queried from each time series
    l_query = numpy.zeros( [ i_scale.shape[0], i_nSynth ] )

    # assemble the query data for each batch
    l_srcsNew = []
    for l_ba in range( i_scale.shape[0] ):
      l_start = -i_shift[l_ba][l_i0]
      l_stop = l_start + i_nSynth / i_freq
      l_query[l_ba] = numpy.linspace( l_start, l_stop, i_nSynth )

      # derive the new source signal
      l_nSrcNew = int(i_nSrc*i_stretch[l_ba][l_i0]+0.5)
      l_sca = float(i_nSrc) / l_nSrcNew
      l_srcsNew.append( numpy.repeat( [l_sca], l_nSrcNew ) )

    # iterate over time series data
    for l_i1 in numpy.ndindex( l_data[l_i0].shape[:-1] ):
      # iterate over batch and do the interpolation
      for l_ba in range( i_scale.shape[0] ):
        # ignore zero strengths
        if( i_scale[l_ba][l_i0] == 0 ): continue

        # convolve with new source
        l_signal = edge_learn.augment.DeConvolve.convSrc( i_nSynth,
                                                          l_srcsNew[l_ba],
                                                          l_data[l_i0][l_i1] )

        # add contribution
        l_out[l_ba][l_i1] = l_out[l_ba][l_i1] + i_scale[l_ba][l_i0] * numpy.interp( l_query[l_ba], l_time, l_signal )
  
  return l_out