#!/usr/bin/env python3
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
import numpy
import copy

##
# Parses the control data of the prorgam TF_MISFIT_GOF_CRITERIA.
#
# @param i_file path to the file containing the control data.
# @return 1) boolean if the control data is valid, 2) parsed control data.
##
def tfMisContDat( i_file ):
  # open file
  with open( i_file, 'r' ) as l_fh:
    # create csv reader
    l_reader = csv.reader( l_fh, delimiter=' ', skipinitialspace=True )

    # config
    l_conf = []

    # parse file
    for l_ro in l_reader:
      l_conf = l_conf + [ l_ro ]

    l_valid = True

    # check if we have valid data
    if any( ('NaN' in li) or ('Infinity' in li) for li in l_conf ):
      l_valid = False

    # parse control data
    assert( len(l_conf) >= 10)

    l_cdat = { 'fmin':         l_conf[0][0],
               'fmax':         l_conf[0][1],
               'nf_tf':        l_conf[1][0],
               'mt':           l_conf[1][1],
               'dt':           l_conf[2][0],
               'nc':           l_conf[2][1],
               'em':           [],
               'pm':           [],
               'eg':           [],
               'pg':           []
             }

    # convert string to int
    for l_ke in ['nf_tf', 'mt', 'nc']:
      l_cdat[l_ke] = int( l_cdat[l_ke] )

    # extract single valued misfits
    for l_co in range(l_cdat['nc']):
      l_cdat['em'] = l_cdat['em'] + [ l_conf[4 + l_co][0] ]
      l_cdat['pm'] = l_cdat['pm'] + [ l_conf[4 + l_co][1] ]

      l_cdat['eg'] = l_cdat['eg'] + [ l_conf[4 + l_cdat['nc'] + l_co][0] ]
      l_cdat['pg'] = l_cdat['pg'] + [ l_conf[4 + l_cdat['nc'] + l_co][1] ]


    # convert to float and int if valid
    if l_valid:
      for l_ke in [ 'fmin', 'fmax', 'dt', 'nc' ]:
        l_cdat[l_ke] = float( l_cdat[l_ke] )

      for l_ke in [ 'em', 'pm', 'eg', 'pg' ]:
        for l_co in range(len(l_cdat[l_ke])):
          l_cdat[l_ke][l_co] = float( l_cdat[l_ke][l_co] )

    # return data
    return l_valid, l_cdat

##
# Generate a 2D, time-frequency plot.
#
# @param i_gs grid specs for the image and the colorbar.
# @param i_data numpy array containing the raw data.
# @param i_style style of the plot ( 'tfrs1', 'tfrs2', 'tfem', 'tfpm', 'tfeg', or 'tfpg' ).
# @param i_rangeFreq frequency range.
# @param i_rangeTime time range.
##
def timeFreq( i_gs,
              i_data,
              i_style,
              i_rangeFreq,
              i_rangeTime ):
  # dictionary containing titles
  l_titles = { 'tfrs1': 'Time-Frequency Representation of the input Signal (TFRS1)',
               'tfrs2': 'Time-Frequency Representation of the reference Signal (TFRS2)',
               'tfem': 'Time-Frequency Envelope Misfit (TFEM)',
               'tfpm': 'Time-Frequency Phase Misfit (TFPM)',
               'tfeg': 'Time-Frequency Envelope Goodness-of-fit (TFEG)',
               'tfpg': 'Time-Frequency Phase Goodness-of-fit (TFPG)' }

  # number of frequency samples
  l_nfs, l_nts = i_data.shape

  # copy over data
  l_data = copy.deepcopy( i_data )

  # default value of the color basr
  l_cbarAno = ''

  # set style-specific properties
  if i_style in ['tfeg', 'tfpg']:
    # derive color scale
    l_vmin = 0
    l_vabs = 10.0
    l_map = matplotlib.pyplot.get_cmap('hot')
  elif i_style in ['tfrs1', 'tfrs2']:
    l_cbarAno = '%'
    l_data = l_data * 100
    l_vabs = abs( max( l_data.min(), l_data.max(), key=abs) )
    l_vmin = 0
    l_map = matplotlib.pyplot.get_cmap('binary')
  else:
    l_cbarAno = '%'
    l_data = l_data * 100
    l_vabs = abs( max( l_data.min(), l_data.max(), key=abs) )
    l_vmin = -l_vabs
    l_map = matplotlib.pyplot.get_cmap('seismic')

  # extend of the plot
  l_extend = [  i_rangeTime[0],  i_rangeTime[1],  i_rangeFreq[0],  i_rangeFreq[1] ]

  # produce plot
  matplotlib.pyplot.subplot( i_gs[0] )
  l_ax = matplotlib.pyplot.imshow( l_data, cmap=l_map, origin='lower', vmin=l_vmin, vmax=l_vabs, aspect='auto', extent=l_extend  )

  # annotate plot
  matplotlib.pyplot.ylabel( 'frequency (Hz)' )

  if i_style in l_titles.keys():
    matplotlib.pyplot.title( l_titles[i_style] )

  l_axes = matplotlib.pyplot.gca()
  l_axes.grid(which='major', color='black', linestyle=':', linewidth=1)

  # add color bar
  l_ax = matplotlib.pyplot.subplot( i_gs[1] )
  l_cb = matplotlib.pyplot.colorbar( label=l_cbarAno, cax=l_ax )

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
    # detect dialect
    l_sniff = csv.Sniffer().sniff( l_fi.read(), [',', ' '] )
    l_fi.seek(0)

    l_csv = csv.reader( l_fi, l_sniff, skipinitialspace=True )

    for l_ro in l_csv:
      if len( l_ro ) > 0:
        # ignore comments
        if l_ro[0].lstrip()[0] == '#':
          continue

        # ignore indices out of range and non-numbers
        try:
          l_vals = l_vals + [ float(l_ro[i_col]) ]
        except IndexError: continue
        except ValueError: continue

  return l_vals

##
# Plots the two time series.
#
# @param i_title title of the plot.
# @param i_recvNames names of the recivers in the plot.
# @param i_first first time series. Two-element list of lists. First list is time, second values.
# @param i_second second time series. Two-element list of lists. First list is time, second values.
##
def plot( i_title, i_recvNames, i_first, i_second ):
  # determine time range
  l_tMin = max( min(i_first[0]), min(i_second[0]) )
  l_tMax = min( max(i_first[0]), max(i_second[0]) )

  # add data
  matplotlib.pyplot.plot( i_first[0],  i_first[1],  label=i_recvNames[0] )
  matplotlib.pyplot.plot( i_second[0], i_second[1], label=i_recvNames[1] )
  matplotlib.pyplot.legend()

  # set limit
  matplotlib.pyplot.xlim( l_tMin, l_tMax )

  # add title
  matplotlib.pyplot.title( i_title )

  # create legends
  matplotlib.pyplot.xlabel('time (s)')

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

l_parser.add_argument( '--in_tmgc',
                       dest     = 'in_tmgc',
                       required = False,
                       help     = 'Paths to direcory conaining output of program TF-MISFIT_GOF_CRITERIA.' )

l_parser.add_argument( '--cols',
                       dest     = 'cols',
                       required = True,
                       nargs    = 2,
                       help     = 'Columns in the two files which are compared. Counting starts at 0, thus the first column (typically time) as id 0.',
                       metavar  = ('COL_0', 'COL_1') )

l_parser.add_argument( '--tmgc_comp',
                       dest     = 'tmgc_comp',
                       required = False,
                       type     = int,
                       default  = 1,
                       help     = 'Component (as parsed TF-MISFIT_GOF-CRITERIA) to use for the tmgc-plots.' )

l_parser.add_argument( '--out_pdf',
                       dest     = 'out_pdf',
                       required = True,
                       help     = 'Paths to output pdf' )

l_parser.add_argument( '--shift',
                       dest     = 'shift',
                       required = False,
                       nargs    = 2,
                       type     = float,
                       default  = [0, 0],
                       help     = 'Shifts the seismogram in time.' )

l_parser.add_argument( '--title',
                       dest     = 'title',
                       required = False,
                       type     = str,
                       default  = 'Receiver Comparison',
                       help     = 'Title of the plot' )

l_parser.add_argument( '--recv_names',
                       dest     = 'recv_names',
                       required = False,
                       nargs    = 2,
                       type     = str,
                       default  = ['first', 'second'],
                       metavar  = ('FIRST_NAME', 'SECOND_NAME'),
                       help     = 'Paths to output pdf' )

l_args = vars(l_parser.parse_args())

logging.info( 'parsing receivers '+l_args['in_csv'][0]+' '+l_args['in_csv'][1] )

# create pdf
l_pdf = PdfPages( l_args['out_pdf'] )

# read TF_MISFIT_GOF_CRITERIA INFO
if l_args['in_tmgc']:
  l_valid, l_cDat = tfMisContDat( l_args['in_tmgc'] + '/MISFIT-GOF.DAT' )
else: l_valid = False

# create figure
if l_valid:
  l_fig = matplotlib.pyplot.figure( figsize=(10, 21) )
else:
  l_fig = matplotlib.pyplot.figure( figsize=(10, 7) )


# read receiver info
l_first  = [ readCsv( l_args['in_csv'][0], 0 ),
             readCsv( l_args['in_csv'][0], int(l_args['cols'][0]) ) ]
l_second = [ readCsv( l_args['in_csv'][1], 0 ),
             readCsv( l_args['in_csv'][1], int(l_args['cols'][1]) ) ]

# shift receivers in time
l_first[0]  = [ l_t + l_args['shift'][0] for l_t in l_first[0] ]
l_second[0] = [ l_t + l_args['shift'][1] for l_t in l_second[0] ]

# use subplots if time-frequency are part of the figure
if( l_valid ):
  l_r = 3.5
  l_gs = matplotlib.gridspec.GridSpec(7, 2, height_ratios=[6,l_r,l_r,l_r,l_r,l_r,l_r], width_ratios=[25,1], wspace=0.1, hspace=0.5)
  matplotlib.pyplot.subplot( l_gs[0] )

# plot receiver
plot( l_args['title'], l_args['recv_names'], l_first, l_second )

# read and plot TF-MISFIT_GOF_CRITERIA data, if available
if l_args['in_tmgc']:
  if l_valid:
    # assemble frequency and time range
    l_fr = ( l_cDat['fmin'], l_cDat['fmax'] )
    l_tr = ( 0, (l_cDat['mt']-1)*l_cDat['dt'])

    # id of the subplot
    l_spId = 2

    for l_en in [ ['tfrs1', 'TFRS1_' +str( l_args['tmgc_comp'] )+ '.DAT' ],
                  ['tfrs2', 'TFRS2_' +str( l_args['tmgc_comp'] )+ '.DAT' ],
                  ['tfem',  'TFEM'   +str( l_args['tmgc_comp'] )+ '.DAT' ],
                  ['tfpm',  'TFPM'   +str( l_args['tmgc_comp'] )+ '.DAT' ],
                  ['tfeg',  'TFEG'   +str( l_args['tmgc_comp'] )+ '.DAT' ],
                  ['tfpg',  'TFPG'   +str( l_args['tmgc_comp'] )+ '.DAT' ] ]:
      # add dir
      l_en[1] = l_args['in_tmgc'] + '/' + l_en[1]

      logging.info( 'parsing time-frequency data '+l_en[1] )

      # read data
      l_data = numpy.genfromtxt( l_en[1] )

      # reshape data
      l_data.reshape( l_cDat['nf_tf'], l_cDat['mt'] )

      # generate plot
      timeFreq( [ l_gs[l_spId], l_gs[l_spId+1] ], l_data, l_en[0], l_fr, l_tr )

      # update sub-plot id
      l_spId = l_spId + 2

# save figure
l_pdf.savefig( l_fig, bbox_inches='tight' )

# close pdf
l_pdf.close()

logging.info( 'got it, bye bye' )
