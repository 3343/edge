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
# Converts plain receivers (1 value per line) to csv-receiver data.
##
import logging
import argparse
import csv

##
# Reads the plain receiver input from the given file.
#
# @param i_file input file which is parsed.
# @return list of lists, where the first list contains the time values and each of the following lists represents the values of a single receiver.
##
def readPlain( i_file ):
  with open( i_file ) as l_in:
    # read header
    l_header = l_in.readline().split()

    # check for three entries
    assert( len(l_header) == 3 )

    # parse header
    l_nSamples = int(   l_header[0].translate(None, '#') )
    l_freq     = float( l_header[1]                      )
    l_nRecvs   = int(   l_header[2]                      )

    # read receiver data as float
    l_dataRaw = []
    for l_ln in l_in.readlines():
      l_ln = l_ln.strip()
      # read only for non-empty lines
      if l_ln != '': l_dataRaw = l_dataRaw + [ float( l_ln ) ]

    # check the the total number of samples matches the data
    assert( len(l_dataRaw) % l_nSamples == 0 )

    # receiver data, initialized with time dimension
    l_dataRecv = [ [ x * l_freq for x in range( 0, l_nSamples ) ] ]

    # iterate over receiver
    for l_re in xrange( l_nRecvs ):
      # add new list
      l_dataRecv = l_dataRecv + [[]]

      # iterate over this receiver's entries
      for l_en in xrange( l_nSamples ):
        # get id in receiver data
        l_id = l_en + l_re*l_nSamples

        # add datum
        l_dataRecv[-1] = l_dataRecv[-1] + [ l_dataRaw[l_id] ]

    return l_dataRecv

##
# Reads the plain receiver input from the given file.
#
# @param i_data list of lists, where the first list contains the time values and each of the following lists represents the values of a single receiver.
# @param i_basePath base path to which each of the receivers is written.
##
def writeCsv( i_data, i_basePath ):
  # iterate over receivers
  for l_re in range( 1, len(i_data) ):

    l_out = i_basePath+'_'+str(l_re)+'.csv'
    logging.info( "writing csv-file "+l_out )

    # write data
    with open( l_out, 'w' ) as l_file:
      l_csv = csv.writer( l_file )

      l_rows = zip( i_data[0], i_data[l_re] )

      for l_ro in l_rows:
        l_csv.writerow( l_ro )

# set up logger
logging.basicConfig( level=logging.DEBUG,
                     format='%(asctime)s - %(name)s - %(levelname)s - %(message)s' )

# command line arguments
l_parser = argparse.ArgumentParser( description='Converts plain receivers (1 value per line) to csv-receiver data.' )

l_parser.add_argument( '--in_file',
                       dest     = 'in_file',
                       required = True,
                       help     = 'Path to the plain-file which will be converted to csv.' )

l_parser.add_argument( '--out_base',
                       dest     = 'out_base',
                       required = True,
                       help     = 'Base-path to the csv-files which will be generated from the plain-file. The base path will be extended with \'_RECVID.csv\', where RECVID is, starting from 1, the id of the receiver.' )

l_arguments = vars(l_parser.parse_args())

logging.info( "parsing "+l_arguments['in_file'] )
l_data = readPlain( l_arguments['in_file'] )

writeCsv( l_data, l_arguments['out_base'] )

logging.info( "all done, see you later" )
