##
# @file This file is part of EDGE.
#
# @author Alexander Breuer (anbreuer AT ucsd.edu)
#
# @section LICENSE
# Copyright (c) 2017, Regents of the University of California
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
# Plots matrices.
##
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import MaxNLocator
import numpy

##
# Plots the sparsity pattern of a given matrix
#
# @param i_mat matrix, which is plotted.
# @param i_out output path.
##
def sparsity( i_mat, i_out ):
  # open pdf
  with PdfPages(i_out) as l_pdf:
    # create local copy of matrix
    l_mat = numpy.matrix( i_mat )

    # create image by including color channels
    l_pat = numpy.ones( l_mat.shape + (3,) )

    # assign black to non-zeros
    l_pat[ abs(l_mat) > 0 ] = [0,0,0]

    # create figure
    l_ax =  matplotlib.pyplot.figure().gca()

    # plot the image
    l_ax.imshow( l_pat,  interpolation='nearest' )

    # add number of non-zeros to title
    l_nnz = numpy.count_nonzero( l_mat )
    matplotlib.pyplot.title( '#non-zero: '+str(l_nnz) )

    # enforce int labels
    l_ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    l_ax.yaxis.set_major_locator(MaxNLocator(integer=True))

    # save figure
    l_pdf.savefig()

    # shutdown pyplot
    matplotlib.pyplot.close()