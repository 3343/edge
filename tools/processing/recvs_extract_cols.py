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
# Extracts the given columns from a receiver.
##
import logging
import argparse
import csv

# set up logger
logging.basicConfig( level=logging.INFO,
                     format='%(asctime)s - %(name)s - %(levelname)s - %(message)s' )

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
# Writes a CSV file.
#
# @param i_data data which is written (list of lists).
# @param i_file file which is written.
##
def writeCsv( i_data, i_file ):
  with open(i_file, 'w') as l_fi:
    # zip list
    l_rows = zip(*i_data)

    l_csv = csv.writer( l_fi, delimiter=',' )
    l_csv.writerows( l_rows )

# command line arguments
l_parser = argparse.ArgumentParser( description='Extracts the given columns from a receiver. All header lines are ignored' )

l_parser.add_argument( '--in_csv',
                       dest     = 'in_csv',
                       required = True,
                       type     = str,
                       help     = 'Paths of the input receiver.' )
l_parser.add_argument( '--cols',
                       dest     = 'cols',
                       required = True,
                       nargs    = '+',
                       type     = int,
                       help     = 'Columns which will be extracted (first column has id 0).',
                       metavar  = ('COL_0', 'COL_1') )

l_parser.add_argument( '--out_csv',
                       dest     = 'out_csv',
                       required = True,
                       help     = 'Path to output receiver.' )
l_args = vars(l_parser.parse_args())

# read data
l_data = []

for l_co in l_args['cols']:
  logging.info( 'reading column ' + str(l_co) +' of receiver '+l_args['in_csv'] )
  l_data = l_data + [ readCsv( l_args['in_csv'], l_co ) ]

writeCsv( l_data, l_args['out_csv'] )

logging.info( 'done with extracting' )
