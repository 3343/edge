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
# Compares two receivers in csv-format.
##
import logging
import argparse
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
from matplotlib.backends.backend_pdf import PdfPages

##
# Reads the given column from the file. Header columns are ignored.
#
# @param i_file file which is read.
# @param i_col column which is read (0 is first column).
##
def readCsv( i_file, i_col ):
  # values
  l_vals = []

  # read file
  with open( i_file, 'r' ) as l_fi:
    l_csv = csv.reader( l_fi )

    for l_ro in l_csv:
      if len( l_ro ) > 0:
        # ignore non-numbers
        try:
          l_vals = l_vals + [ float(l_ro[i_col]) ]
        except ValueError: continue

  return l_vals

##
# Plots the two time series.
#
# @parma i_recvNames names of the recivers in the plot.
# @param i_first first time series. Two-element list of lists. First list is time, second values.
# @param i_second second time series. Two-element list of lists. First list is time, second values.
# @param i_outPdf path to outpfule pdf-file which will be written.
##
def plot( i_recvNames, i_first, i_second, i_outPdf ):
  # create pdf
  l_pdf = PdfPages( i_outPdf )

  # create figure
  l_fig = matplotlib.pyplot.figure( figsize=(10, 7) )
  # add data
  matplotlib.pyplot.plot( i_first[0],  i_first[1],  label=i_recvNames[0] )
  matplotlib.pyplot.plot( i_second[0], i_second[1], label=i_recvNames[1] )
  matplotlib.pyplot.legend()

  # create legends
  matplotlib.pyplot.xlabel('time')
  matplotlib.pyplot.ylabel('values')

  # save figure
  l_pdf.savefig( l_fig )

  # close pdf
  l_pdf.close()

# set up logger
logging.basicConfig( level=logging.DEBUG,
                     format='%(asctime)s - %(name)s - %(levelname)s - %(message)s' )

# command line arguments
l_parser = argparse.ArgumentParser( description='Converts plain receivers (1 value per line) to csv-receiver data.' )

l_parser.add_argument( '--in_csv',
                       dest     = 'in_csv',
                       required = True,
                       nargs    = 2,
                       help     = 'Paths of the two receivers',
                       metavar  = ('RECV_0', 'RECV_1') )

l_parser.add_argument( '--cols',
                       dest     = 'cols',
                       required = True,
                       nargs    = 2,
                       help     = 'Columns in the two files which are compared. Counting starts at 0, thus the first column (typically time) as id 0.',
                       metavar  = ('COL_0', 'COL_1') )

l_parser.add_argument( '--out_pdf',
                       dest     = 'out_pdf',
                       required = True,
                       help     = 'Paths to output pdf' )

l_parser.add_argument( '--recv_names',
                       dest     = 'recv_names',
                       required = False,
                       nargs    = 2,
                       type     = str,
                       default  = ['first', 'second'],
                       metavar  = ('FIRST_NAME', 'SECOND_NAME'),
                       help     = 'Paths to output pdf' )

l_args = vars(l_parser.parse_args())

logging.info( 'reading receivers '+l_args['in_csv'][0]+' '+l_args['in_csv'][1] )

# read receiver info
l_first  = [ readCsv( l_args['in_csv'][0], 0 ),
             readCsv( l_args['in_csv'][0], int(l_args['cols'][0]) ) ]
l_second = [ readCsv( l_args['in_csv'][1], 0 ),
             readCsv( l_args['in_csv'][1], int(l_args['cols'][1]) ) ]

plot( l_args['recv_names'], l_first, l_second, l_args['out_pdf'] )

logging.info( 'got it, bye bye' )
